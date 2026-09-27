# 共识引擎 v6.6f 配套清理：8 数据集逐级相容性对照与 Phylum-Class 定性

> **版本状态（2026-09-16 补注）**：本文 §1–§7 记录的是 **v6.6f 候选（`/tmp/rc66_multi`，未上线）与当时线上 v6.6 世代成品** 的对照，以及由此提出的 A/B/C 三方案。**§6.5 的核心论断「半数票永不触发淘汰」已被 v6.7 取代**：2026-09-16 大王拍板把淘汰门槛语义改为「半数票也淘汰」（`share >= CASCADE_MIN_SHARE`，默认 0.5，`TAX_GATE_VERSION="6.7"`），8 数据集已全量重跑，见新增 **§8**。读本文时请注意「线上」列 = 刷新前线上成品（v6.6 世代），不是 v6.7。

日期：2026-08-30
脚本：`/tmp/chimera_multi.py`（参照库口径，判据与 `/tmp/chk_chimera.py` 完全一致）、
`/tmp/diag_phylum_class.py`、`/tmp/quant_phylum_class_fix.py`、`/tmp/cmp_bad_contig_size.py`
输入：8 数据集线上成品 `05_Taxonomy/Votus.integrated/final_integrated_classification.tsv` vs v6.6f 产物 `/tmp/rc66_multi/<数据集>/final_integrated_classification.tsv`
参照：`~/database/taxonomy/rankedlineage.dmp`

## 1. 8 数据集参照口径对照（相邻阶元矛盾行）

| 数据集 | 行数 | 线上矛盾行 | 线上矛盾率 | v6.6f 矛盾行 | v6.6f 矛盾率 |
|---|---|---|---|---|---|
| RNA-Alternaria_alternata_out | 1728 | 356 | 20.60% | 62 | 3.59% |
| RNA-Aphis_gossypii_out | 2448 | 418 | 17.08% | 59 | 2.41% |
| RNA-Fusarium_nematophilum_out | 1637 | 290 | 17.72% | 49 | 2.99% |
| RNA-Lycium_amarum_out | 957 | 167 | 17.45% | 9 | 0.94% |
| RNA-Lycium_barbarum_out | 20892 | 2738 | 13.11% | 465 | 2.23% |
| RNA-Lycium_chinense_out | 11007 | 1258 | 11.43% | 232 | 2.11% |
| RNA-Lycium_ruthenicum_out | 15181 | 1838 | 12.11% | 322 | 2.12% |
| onekp-virus | 535103 | 71295 | 13.32% | 21541 | 4.03% |

barbarum 用旧脚本复算一致（2738 / 465、Phylum-Class 97 / 273），判据等价性已验证。

## 2. 唯一恶化的边界：Phylum <-> Class

| 数据集 | 线上 | v6.6f |
|---|---|---|
| Alternaria | 4 | 41 |
| Aphis | 5 | 33 |
| Fusarium | 5 | 31 |
| amarum | 10 | 6 |
| barbarum | 97 | 273 |
| chinense | 28 | 117 |
| ruthenicum | 32 | 200 |
| onekp | 14276 | 12676 |

6/8 上升，2/8 下降（amarum、onekp）。其余边界一律大幅下降（barbarum 科-属 1273 -> 9、属-种 1011 -> 62）。

## 3. 成因链（实测，非推断）

样例 `CRR1126134_clean_NODE_21_length_1900_cov_0.549535`（barbarum，工具原始报告）：

| 工具 | Realm | Kingdom | Phylum | Class | Order | Family |
|---|---|---|---|---|---|---|
| mmseqs | Varidnaviria | Bamfordvirae | Nucleocytoviricota | Megaviricetes | Pimascovirales | Pithoviridae |
| VITAP | Duplodnaviria | Heunggongvirae | Uroviricota | Caudoviricetes | NA | NA |

v6.6f 成品该行：Realm=Duplodnaviria(VITAP) Kingdom=Heunggongvirae(VITAP) Phylum=Uroviricota(VITAP)
**Class=Megaviricetes(mmseqs)** Order=Pimascovirales(mmseqs) Family=Pithoviridae(mmseqs)。

机制（对应 `virus_classifier_analysis.R` md5 `c77c43e2468ad3de48466da3e89fff9e`）：

1. Phylum 两工具 1-1 平票，`share = 0.5`，`share > CASCADE_MIN_SHARE` 不成立（L110 注释明确「并列时永不触发」）-> 不淘汰任何人（规则 #4「无明确多数不淘汰」）。
2. Class 两候选仍是 1-1 平票，按下一条已定规则「先用下一级一致性」（L727-739）：`next_score = modal/n_tools`。
   mmseqs 有 Order=Pimascovirales -> 1/1 = 1；VITAP 的 Order 为 NA -> 0/1 = 0 -> mmseqs 的 Megaviricetes 胜出。
3. 结果：粗阶元留在 VITAP 阵营，Class 起翻回 mmseqs 阵营，形成跨阵营世系。

即 Phylum-Class 拼接率上升是「平票不淘汰」+「平票用下一级一致性裁决」两条已定规则的联合后果，不是新回归。

## 4. 273 行全覆盖量化（barbarum）

`/tmp/quant_phylum_class_fix.py` 结果：

- Phylum<->Class 参照矛盾行：273
- 其中「Phylum 与 Class 的支持工具完全不重叠」的跨阵营翻转：**273（100%）**
- 其中「Phylum 阵营工具自己报的 Class 与 Phylum 参照相容」（存在可替代值）：**273（100%）**
- 其中「Class 阵营工具自己报的 Phylum 与参照一致」（错在粗阶元一侧）：**273（100%）**

## 5. 矛盾行的 contig 规模（barbarum）

| 集合 | n | length 中位 | 四分位 | cov 中位 | cov<1 占比 |
|---|---|---|---|---|---|
| 全体行 | 20229 | 1589 | 809/2369 | 3.21 | 7.3% |
| v6.6f P<->C 矛盾行 | 263 | 2138 | 1560/2748 | 3.43 | 3.8% |
| 线上 P<->C 矛盾行 | 89 | 1800 | 933/2469 | 3.00 | 5.6% |

矛盾行并不集中在低覆盖短 contig，反而略长于全体中位。属两工具对同一 contig 给出不同分类的实质分歧，不是噪声行。

## 6. 三个可选补丁（待决策）

| 方案 | 做法 | 预期覆盖 barbarum 273 行 | 代价 |
|---|---|---|---|
| A 平票也淘汰 | 平票时按 `TIE_BREAK_ORDER` 定音并淘汰异议工具（改规则 #4） | 273 | 被淘汰工具在该行的 Order/Family 一并丢失（本例丢 Pimascovirales/Pithoviridae），与「避免误杀」取向冲突 |
| B Phylum-Class 参照闸门 | 投票时或后处理用 rankedlineage 校验，Class 收敛到 Phylum 阵营 | 273 | 同样向下波及 Order/Family，本质保粗弃细 |
| C 粗阶元服从细阶元 | Class 与 Phylum 参照矛盾时，Phylum := ref_phylum_of_class(Class) | 273 | 新增一条「矛盾时粗阶元由细阶元+参照决定」的规则；保留更具体的科属信息，成品与参照完全相容 |

三者都能把 barbarum 矛盾行从 465 降到约 192（0.92%）。C 是唯一既消矛盾又不丢信息的方案，但改动语义最大；A/B 属于保粗弃细。

## 6.5 为什么 v6.6f 矛盾率还不是 0（三层拆解，全部实测）

### L1 跨阵营拼接：占残余行 95.9%~100%（引擎侧）

`/tmp/decomp_residual.py`（8 数据集）：按 `*_agree` 支持工具集合是否相交把矛盾行分桶。

| 数据集 | 矛盾行 | 跨阵营行 | 同阵营行 | 跨阵营占比 | 主导跨阵营边界 |
|---|---|---|---|---|---|
| Alternaria | 62 | 61 | 1 | 98.4% | Phylum-Class 41/69 |
| Aphis | 59 | 59 | 0 | 100.0% | Phylum-Class 33/60 |
| Fusarium | 49 | 47 | 2 | 95.9% | Phylum-Class 31/53 |
| amarum | 9 | 9 | 0 | 100.0% | Phylum-Class 6/9 |
| barbarum | 465 | 456 | 9 | 98.1% | Phylum-Class 273/494 |
| chinense | 232 | 229 | 3 | 98.7% | Phylum-Class 117/247 |
| ruthenicum | 322 | 318 | 4 | 98.8% | Phylum-Class 200/340 |
| onekp | 21541 | 21240 | 301 | 98.6% | Phylum-Class 12676/22099 |

跨阵营矛盾对里，上层阶元的票首占比分布（`/tmp/decomp_layers.py`）：

| 上层票首占比 | barbarum（494 对） | onekp（22099 对） |
|---|---|---|
| 1/2 = 0.50 | 194（39.3%） | 18285（82.7%） |
| 2/4 = 0.50 | 155（31.4%） | 1839（8.3%） |
| 3/6 = 0.50 | 51（10.3%） | 22（0.1%） |
| 1/3 | 33（6.7%） | 1232（5.6%） |
| 1/1 | 24（4.9%） | 380（1.7%） |
| 2/5 | 20（4.0%） | 249（1.1%） |
| 其他 | 17 | 92 |
| 本应触发淘汰（> 0.5） | 25（5.1%） | 381（1.7%） |

即：跨阵营矛盾对里，上层阶元票首**恰好等于一半**（1/2、2/4、3/6）的占 barbarum 81.0%、onekp 91.1%。规则要求 `share > CASCADE_MIN_SHARE`（0.5），半数票不满足，永不触发淘汰（L110 注释），下层阶元便可能落到另一阵营。

> ⚠️ **已被 v6.7 取代（2026-09-16）**：现行实现为 `share >= CASCADE_MIN_SHARE`，半数票（`share == 0.5`）**会**触发淘汰。本段描述的「半数票永不淘汰」只适用于 v6.6f 及更早版本。见 §8。

### L1 的实证：门槛降到 0.499

`MMPV_CASCADE_MIN_SHARE=0.499`（`/tmp/thr_exp_0499.sh`，无代码改动，只用已上线脚本的环境变量），barbarum：

| 口径 | 矛盾行 | 矛盾率 | Phylum-Class |
|---|---|---|---|
| 线上 | 2738 | 13.11% | 97 |
| v6.6f（门槛 0.5） | 465 | 2.23% | 273 |
| 门槛 0.499 | **87** | **0.42%** | **13** |

`/tmp/cmp_thr_vs_v66f.py`：Phylum<->Class 参照相容性 both_ok 20619 / **fixed 260** / **broken 0** / both_bad 13。
代价（同一批行）：Order 有值->空 263 行、Family 有值->空 216 行；Class 取值变化 278 行、Order 334、Family 323、Genus 282、Species 252；反向（空->有值）为 0。

### L2 缺值不算异议（引擎侧）

上层阶元只有 1 个工具报值（1/1，barbarum 4.9%、onekp 1.7%）时没有异议方可淘汰；未在该阶元报值的工具不属异议、不被淘汰，仍可在更细阶元投票，于是仍能产出与已定上层不相容的下层值。5.1% / 1.7% 那些「票首已过 0.5」却仍矛盾的矛盾对，部分也源于此。

### L3 工具自身与参照不相容（工具/参照侧，占残余行 0.9%~4.1%）

同阵营矛盾行，barbarum 10 行、onekp 314 行。反复出现的形态：

- Family=Adenoviridae 且 Genus=Organic Lake phycodnavirus（参照库该属父级为 Phycodnaviridae），支持工具 CAT+diamond_lca、metabuli，如 `CRR732652_clean_NODE_613`、`SRR23107139_clean_NODE_4887`。
- Realm=Monodnaviria 且 Kingdom=Sangervirae（参照库父级为 Floreoviria），仅 genomad 单工具，如 `SRR23107136_clean_NODE_28958`。

这类是工具所据分类体系与 NCBI 参照库不一致（或参照滞后），任何投票规则都消不掉，除非让参照库成为裁决者。

### 分母口径注

参照库未收录的名字直接跳过、不计入矛盾（barbarum 14072 个 name 命中 13991，81 个未命中，0.58%；onekp 比例更高）。因此工具侧真实不一致率被轻微低估。

### 结论

矛盾率不为 0 是「裁决者是工具票数」这一前提的必然结果：只要票数是唯一权威，就总会有与参照不相容的票首组合留下拼接。要压到 0 只能让参照库对每一级都做裁决，那时矛盾的定义与消除同源，等于用参照库而非共识定义分类，与「共识优先」取向冲突。

> **v6.7 修正（2026-09-16）**：本段原写「只在严格多数时淘汰异议方」是 v6.6f 的门槛语义。改 `>=` 后残余矛盾率从 12.11%~20.60% 降到 0.18%~0.52%，但**不为零**，L2（缺值不算异议）与 L3（工具自身与参照不相容）两层仍在，结论方向不变，量级已变。见 §8。

## 7. 本轮已修 / 待定

已修并验证：
- rescue 报告层两处真错误（`prepass`/`failed`/`C blastn_completo` 行、md 分母自相矛盾），8 数据集报告回填，branch 合计 == 候选数全部核对通过。
- Stage 10 总报告元数据区块从不出现：`rglob("**/metadata_association/**")` 只返回目录，改 `**/*`；同时 `_stage_files` 标题改相对路径。刷新后图嵌入 135/135、`data:image` 432 次、metadata 关键词 406 次。
- 孤儿副本 `utils/auto_known_virus.py` 归档为 `archive/auto_known_virus.orphan_utils_20260915.py`。

待定：
- `05_Taxonomy/Votus.integrated/calibration_20260914/` 是否纳入成品（口径 A 校准副产物，输入为线上成品 `12e8f171…`，blank 1246 / review 278 / conflict_one_side 26）。
- rescue md 数字变化（barbarum 合计 308->286、未拯救 748->770、过短 864->749+21）是否同步论文正文。
- 上面 A/B/C 三方案是否上线。

---

## 8. v6.7 定稿与 8 数据集重跑（2026-09-16）

### 8.1 定稿口径

- 淘汰门槛：`share >= CASCADE_MIN_SHARE`（默认 `0.5`），即**半数票也淘汰**；不写魔法数 0.499。等价性已用 md5 证明：新默认产物与 `MMPV_CASCADE_MIN_SHARE=0.499` 实验产物逐字节相同（`a3ead5b5bc4fe290a0e9922497d2fc35`，barbarum）。
- `TAX_GATE_VERSION` 由 `6.6f` bump 到 **`6.7`**（纯数字，`_ver()` 按 `.` split 后 `int()`）。
- 官方脚本：`/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R`，md5 **`f0683ac6af1a9ee1395a100d128365ca`**，备份 `.bak_thr0499_20260916`。
- 执行器：`scripts/audit/refresh_v67.sh`（服务器 `/tmp/refresh_v67.sh`，md5 `d832c4868a2b8298d7bfca44be0dcb18`）。
- 未跑 `--stage analysis`，因此成品未引入 `Nucleic_acid` 列（超范围 schema 变化）。

### 8.2 8 数据集矛盾率（参照口径 R1，本报表口径）

| 数据集 | 数据行 | 旧（刷新前线上）矛盾行 | 旧率 | 新（v6.7）矛盾行 | 新率 | 降幅 |
|---|---|---|---|---|---|---|
| RNA-Alternaria_alternata_out | 1728 | 356 | 20.60% | 6 | 0.35% | 59× |
| RNA-Lycium_amarum_out | 957 | 167 | 17.45% | 2 | 0.21% | 84× |
| RNA-Aphis_gossypii_out | 2448 | 418 | 17.08% | 7 | 0.29% | 60× |
| RNA-Lycium_barbarum_out | 20892 | 2738 | 13.11% | 87 | 0.42% | 31× |
| RNA-Lycium_chinense_out | 11007 | 1258 | 11.43% | 57 | 0.52% | 22× |
| RNA-Fusarium_nematophilum_out | 1637 | 290 | 17.72% | 3 | 0.18% | 97× |
| onekp-virus | 535103 | 71295 | 13.32% | 1460 | 0.27% | 49× |
| RNA-Lycium_ruthenicum_out | 15181 | 1838 | 12.11% | 51 | 0.34% | 36× |

8 个数据集全部：`gate_version 6.7`、`script_md5 f0683ac6…`、`dual_ref_conflict 0 PASS`、`species_dual_ref_conflict 0 PASS`、行数与 contig_id 集合不变、旧产物备份目录存在。

### 8.3 为什么还不是 0（v6.7 下重述三层）

改 `>=` 后，L1 的「半数票不淘汰」主因被消掉（barbarum 2738→87、97 行 P↔C → 13 行），残余由两层构成：

- **L2 缺值不算异议**：上层阶元仅 1 个工具报值（1/1）时无异议方可淘汰，未在该阶报值的工具仍可在更细阶元投票，于是仍能产出与已定上层不相容的下层值。
- **L3 工具自身与参照不相容**：如 Family=Adenoviridae 且 Genus=Organic Lake phycodnavirus（参照属 Phycodnaviridae），Realm=Monodnaviria 且 Kingdom=Sangervirae（参照父级 Floreoviria），单工具形态。工具分类体系与 NCBI 参照库不一致，任何投票规则都消不掉，除非让参照库成为裁决者。
- 新产物里 barbarum 的 Phylum↔Class 相容性为 `both_ok 20783 / fixed 96 / broken 12 / both_bad 1`：13 行 P↔C 矛盾中 12 行是**新产生的**跨阵营拼接（细阶元服从多数阵营后，粗阶元残留另一阵营），1 行旧产物也坏。这是改 `>=` 的结构性代价，量级远小于收益。

### 8.4 参照口径 R1 与 R2

`/tmp/indep_chimera.py`（不复用管线代码的独立实现）同时给出两口径，barbarum 实测：

| 口径 | 定义 | 旧 | 新 |
|---|---|---|---|
| R1 | 按 `tax_name` 在 `rankedlineage.dmp` **首次出现**行取整条谱系（与本文全部既有数字同口径） | 2738（13.11%） | 87（0.42%） |
| R2 | 按阶元列取**多父级集合**，任一父级相容即不算矛盾 | 2585（12.37%） | 53（0.25%） |

R1 与既有报告逐格一致（7 对分布、合计）。R2 数值更低，差异集中在 Family–Genus 与 Genus–Species，成因是同名多义（同一 name 在 dmp 中以不同阶元出现）与多父级集合包含行内高阶元。**方法章口径待定，倾向 R1**（与既有全部报告一致、可复现）。

### 8.5 逐数据集审计报告

`scripts/audit/refresh_v67/{barbarum,Alternaria,amarum,Aphis,chinense,Fusarium,ruthenicum,onekp}.md`，每份含 gate 原文、独立重算的矛盾率与 7 对分布、逐文件归因、备份路径与未动清单；`indep_verify.md` 为 8 数据集 R1/R2 的独立交叉验证。

### 8.6 下游同步（本轮已做）

- rescue 报告层 `rescue_pipeline.py` 修复（md5 `2960e4fe3496ad9d72b51699d2cd95c3`）+ 8 数据集回填。
- Stage 10 元数据区块修复 `generate_pipeline_report.py`（md5 `12c239956c2e22b22d6c2c6c3eb1ead7`）。
- 孤儿副本 `utils/auto_known_virus.py` 归档为 `archive/auto_known_virus.orphan_utils_20260915.py`。

### 8.7 §6 的 A/B/C 三方案状态

方案 A（平票也淘汰）已由 v6.7 实施，但**只做淘汰，不做「粗阶元服从细阶元」的后处理**；方案 B/C 未上线。v6.7 后 barbarum 残余 87 行，若还要再降需另起方案，目前停在「L2/L3 结构性残余」的定性上。
