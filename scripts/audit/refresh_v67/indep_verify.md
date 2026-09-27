# refresh_v67 独立交叉验证报告（分类矛盾率 R1/R2）

- 生成时间：2026-09-16（本地会话）
- 任务性质：只读 + 本地报告；服务器临时文件仅置于 `/tmp/indep_verify/`
- 纪律：未删除任何文件；未改动管线脚本 / R 脚本 / `refresh_v67.sh` / 既有 `/tmp/refresh_v67/**`；未运行 `virome_pipeline.py` 任何 stage；产物目录只读
- 本报告全部数字来自本人独立脚本的实跑；**未采信任何既有子代理结论**

---

## 1. 脚本与指纹

| 对象 | 路径 | md5 | 备注 |
|---|---|---|---|
| 独立复算脚本（任务提供，未改动） | `/tmp/indep_chimera.py` | `73315578d6160ea6359ac9e7b7cd23cd` | 121 行。自写 csv 解析、自建 R1/R2 两套索引，`grep` 全文件无 `import` 任何管线模块，确认不复用管线代码 |
| 参照库 | `/home/zhangwenda/database/taxonomy/rankedlineage.dmp` | `f6b86ae7000837d638dd18035c7c41d0` | 2,831,927 行，375,904,413 bytes |
| 管线参考实现（仅阅读，未使用） | `/tmp/chimera_multi.py` | `e5f5961e77c6eb016d7f32ff05ed1c7d` | 用于确认 R1 口径 |
| 本人辅助：逐格比对 | `/tmp/indep_verify/compare_chimera.py` | `8be2b4ff59aad2dde0eec91ad78420bb` | R1 vs chimera.txt |
| 本人辅助：R1/R2 分歧明细+抽样 | `/tmp/indep_verify/detail.py` | `3c86e45c2105c987f9e8cf8bcb274d51` | |
| 本人辅助：R1-only 标签集中度 | `/tmp/indep_verify/freq_r1only.py` | `40d8912c3a1cdaa1ad187dd737db69b2` | |
| 本人辅助：整行抓取 | `/tmp/indep_verify/fetch_rows.py` | `103c68ea2c85bb648c7a76a65bfc31d1` | |
| 本人辅助：按名定位实例 | `/tmp/indep_verify/find_one.py` | — | |
| 本人辅助：schema/contig 集合 | `/tmp/indep_verify/schema_check.sh` | `cd9b4cb3c0873ce6b1d816a3f9c40af4` | |
| 运行日志 | `/tmp/indep_verify/logs/<label>.{new,old}.log` | — | 16 份 |

R1 口径（与 `chimera_multi.py` 一致）：按 `tax_name` 在 dmp 中**首次出现**行取整条谱系，查低阶元名在该高阶层元的祖先值，与行内高阶元比较，不等即"该行该对上矛盾"。
R2 口径：按 `(低阶元, 低阶元值)` 聚合所有"该值出现在低阶元列"的行的**高阶取值集合**，行内高阶值不在集合中即矛盾。

---

## 2. R1 与既有 `chimera.txt` 逐格对齐

`compare_chimera.py` 自动比对 8 数据集 × (online=新产物, v66f=旧产物) = 16 行，每行的行数、7 对阶元计数、矛盾行数全部比对：

```
ALL_CELLS_MATCH     （16/16 行、每格一致）
```

**没有一个数据集对不上**，因此无需列出差异行。

---

## 3. 汇总表（旧→新；R1、R2 两口径并列）

矛盾率 = 矛盾行 / 数据行数。

| 数据集 | 行数 | R1 旧产物 | R1 新产物 | R2 旧产物 | R2 新产物 |
|---|---|---|---|---|---|
| Alternaria | 1728 | 356 (20.60%) | **6 (0.35%)** | 335 (19.39%) | 5 (0.29%) |
| amarum | 957 | 167 (17.45%) | **2 (0.21%)** | 166 (17.35%) | 2 (0.21%) |
| Aphis | 2448 | 418 (17.08%) | **7 (0.29%)** | 411 (16.79%) | 5 (0.20%) |
| barbarum | 20892 | 2738 (13.11%) | **87 (0.42%)** | 2585 (12.37%) | 53 (0.25%) |
| chinense | 11007 | 1258 (11.43%) | **57 (0.52%)** | 1237 (11.24%) | 34 (0.31%) |
| Fusarium | 1637 | 290 (17.72%) | **3 (0.18%)** | 289 (17.65%) | 1 (0.06%) |
| onekp | 535103 | 71295 (13.32%) | **1460 (0.27%)** | 70488 (13.17%) | 915 (0.17%) |
| ruthenicum | 15181 | 1838 (12.11%) | **51 (0.34%)** | 1805 (11.89%) | 46 (0.30%) |

旧产物矛盾率区间 11.43%–20.60%，新产物降至 0.18%–0.52%（R1）；两个口径给出的"刷新有效"结论一致。

---

## 4. R1 与 R2 的分歧

### 4.1 规模

| 数据集 | 新产物 R1-only 行 | 旧产物 R1-only 行 | R2-only 行（新/旧） |
|---|---|---|---|
| Alternaria | 1 | 21 | 0 / 0 |
| amarum | 0 | 1 | 0 / 0 |
| Aphis | 2 | 7 | 0 / 0 |
| barbarum | 34 | 153 | 0 / 0 |
| chinense | 23 | 21 | 0 / 0 |
| Fusarium | 2 | 1 | 0 / 0 |
| onekp | 545 | 807 | 0 / 0 |
| ruthenicum | 5 | 33 | 0 / 0 |

**16 个文件里 R2-only 全部为 0，观测上 R2 的矛盾集合是 R1 的严格子集**（R2 ⊆ R1）。分歧 = R1-only 行，且**只出现在最低两对阶元**（`Family-Genus`、`Genus-Species`），其余 5 对（Realm-Kingdom … Order-Family）两口径**逐个数字完全相同**。新产物逐对示例：

| 数据集 | Family-Genus R1/R2 | Genus-Species R1/R2 |
|---|---|---|
| Alternaria | 0 / 0 | 2 / 1 |
| Aphis | 0 / 0 | 3 / 1 |
| barbarum | 8 / 0 | 37 / 11 |
| chinense | 3 / 0 | 27 / 7 |
| Fusarium | 1 / 0 | 1 / 0 |
| onekp | 300 / 0 | 291 / 43 |
| ruthenicum | 2 / 0 | 9 / 6 |

### 4.2 原因（量化）

`detail.py` 对每一例 R1-only 做归因。**全部 16 个文件的 R1-only 计数都落在同一个原因 `A`，没有任何一例来自"多父级集合包含行内高阶元"**（即任务中假设的 C 情形一次都没出现），也没有 B/D/E/F/G/H。

原因 A 的本质：
- `R1` 用**名字索引**（只看 `tax_name`），只要名字在库里出现就能给出它在对应高阶层元的祖先值；
- `R2` 用**阶元列索引** `(低阶元, 值)`，只有当该名字**确实以该低阶元身份出现在库里**（即它被当作祖先、其子代行里带有它）时才有键。对"错阶名"和"终端名"没有键，于是**直接弃权、不判矛盾**。

叠加一个重要的参照库事实：本 `rankedlineage.dmp` 是 **self-exclusive（自排除）** 形态——任何分类单元自身所在阶元列为空，谱系列只写严格高于自身的祖先。例：`Homo sapiens` 行 species 列空、genus=Homo；`root`/`Riboviria` 全空；species 列非空仅 290,168 行（都是种下单元）。因此在 R2 里，一个"叶子种名"连自身那行都不会进 species 列，几乎注定无键。

### 4.3 抽样实例（新产物，每条含 contig / 阶元值 / 参照命中）

| 数据集 | contig_id | 对 | 行内值 | R1 参照祖先 | R2 | 归因 |
|---|---|---|---|---|---|---|
| Alternaria | `SRR23083390_clean_NODE_4577...` | Genus-Species | Genus=Moumouvirus, Species=Pandoravirus quercus | Genus=Pandoravirus | 无键 | A：Pandoravirus quercus 无种下子代，species 列无该名 |
| Aphis | `SRR36562187_clean_NODE_4196...` | Genus-Species | Genus=Moumouvirus, Species=Pandoravirus quercus | Genus=Pandoravirus | 无键 | 同上 |
| Fusarium | `contig_869` | Family-Genus | Family=Adenoviridae, Genus=Organic Lake phycodnavirus | Family=Phycodnaviridae | 无键 | A：该名非属级（Phycodnaviridae 下种级），genus 列无该名 |
| barbarum | `CRR732652_clean_NODE_613...` | Family-Genus | Family=Adenoviridae, Genus=Organic Lake phycodnavirus | Family=Phycodnaviridae | 无键 | 同上 |
| barbarum | `CRR1440122_clean_NODE_16...` | Genus-Species | Genus=Medusavirus, Species=Pandoravirus inopinatum | Genus=Pandoravirus | 无键 | A |
| chinense | `SRR23107134_clean_NODE_2245...` | Family-Genus | Family=Adenoviridae, Genus=Organic Lake phycodnavirus | Family=Phycodnaviridae | 无键 | 同上 |
| chinense | `SRR23107134_clean_NODE_15588...` | Genus-Species | Genus=Theiavirus, Species=Pandoravirus inopinatum | Genus=Pandoravirus | 无键 | A |
| ruthenicum | `CRR527056_clean_NODE_4696...` | Family-Genus | Family=Adenoviridae, Genus=Organic Lake phycodnavirus | Family=Phycodnaviridae | 无键 | 同上 |
| onekp | `ERR2040124_clean_NODE_2483...` | Family-Genus | Family=Adenoviridae, Genus=Organic Lake phycodnavirus | Family=Phycodnaviridae | 无键 | 同上 |
| onekp | `ERR2040533_clean_NODE_997...` | Genus-Species | Genus=Varicellovirus, Species=Pandoravirus salinus | Genus=Pandoravirus | 无键 | A：Pandoravirus salinus 无种下子代 |
| onekp | `ERR2040185_clean_NODE_4893...` | Genus-Species | Genus=Medusavirus, Species=Pandoravirus inopinatum | Genus=Pandoravirus | 无键 | A |

表中"参照命中"列即 R1 对低阶元名首次出现行取到的高阶层元祖先值。

参照库端的具体落位（已核）：
- `Organic Lake phycodnavirus`(tax_id 938083)：谱系为 `Phycodnaviridae / Algavirales / Megaviricetes / Nucleocytoviricota / Bamfordvirae / Varidnaviria`，**是种级名，不在 Adenoviridae**。行里却配 `Family=Adenoviridae`，属真正的跨谱系拼接。
- `Pandoravirus inopinatum`(tax_id 1605721)：`Genus=Pandoravirus`（family/order 等在库中为空）。行里却配 `Genus=Medusavirus/Mastadenovirus/Theiavirus…`。

### 4.4 残余矛盾的集中度（新产物）

`freq_r1only.py` 显示，新产物的 R1-only 并非弥散噪声，而是被**极少数具体标签**主导：

| 数据集 | R1-only（对-次） | 最主要标签 |
|---|---|---|
| onekp | 548 | `Family=Adenoviridae + Genus=Organic Lake phycodnavirus` **300 次**；其后多为 `Species=Pandoravirus salinus/inopinatum/...` 配错属 |
| barbarum | 34 | `Species=Pandoravirus salinus + Genus=Alphahydrivirus` 11；`Adenoviridae+Organic Lake` 8；余为 Pandoravirus 种配错属 |
| chinense | 23 | `Pandoravirus salinus+Alphahydrivirus` 7；`Adenoviridae+Organic Lake` 3；余同型 |
| ruthenicum | 5 | `Adenoviridae+Organic Lake` 2；`Pandoravirus japonicus` 配错属 2 |
| Fusarium | 2 | `Pandoravirus inopinatum+Mastadenovirus` 1；`Adenoviridae+Organic Lake` 1 |
| Aphis | 2 | 均为 `Pandoravirus quercus/macleodensis` 配错属 |
| Alternaria | 1 | `Pandoravirus quercus+Moumouvirus` |
| amarum | 0 | — |

一句话：新产物的全部残余矛盾基本来自两类反复出现的错配——**①分类器把 `Organic Lake phycodnavirus`（本质是藻病毒）当属、却配了腺病毒科；②把 `Pandoravirus` 属的种配到了其它病毒属**。

### 4.5 判断：哪个口径更适合写进方法（**以下为我的判断，非事实**）

我倾向以 **R1 作为方法正文的主口径**，理由：
1. R1 与管线自身既有的一致性对照实现（`chimera_multi.py`）同源，且 `taxonomy_gate_check.tsv` 的 gate 也是按这一路参照判据建设，口径自洽；
2. R1 对"名字可解析"的行从不弃权，度量的是"该行低阶元名在参照里的归位是否与行内高阶元自洽"，更贴近"分类矛盾"想表达的东西；
3. R2 更保守，但它的"少"在此数据上主要来自**弃权**（R2-only=0、且 100% 落在原因 A），把"名字不在该阶元/无子代"跳过，而非真正接受了一个合法的另一父级。把这样一个系统性偏保守、且叠加了 self-exclusive 参照效应而额外弃权的数当成"更干净的结果"来宣传，风险较大。

同时应指出：**判据换成 R2 不改变结论本身**——旧 11–21% → 新 0.2–0.5% 的量级和方向在两个口径下一致，所以"刷新大幅降低了分类矛盾"这一主结论对口径选择是稳健的。若审稿人质疑"低阶阶元名存在同名录/多父级"，可以补报 R2 作为敏感性分析，但需说明 R2 的弃权机制。

---

## 5. gate 与 schema 核对结论

### 5.1 gate 文件（8 个数据集，全部一致）

`<ROOT>/05_Taxonomy/Votus.integrated/taxonomy_gate_stamp.tsv`：

| 字段 | 值 |
|---|---|
| gate_version | `6.7` |
| script_path | `/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R` |
| script_md5 | `f0683ac6af1a9ee1395a100d128365ca` |
| ref_file | `/home/zhangwenda/database/taxonomy/genus_family_ref.tsv` |
| ref_md5 | `865283b4b1a3a2f3c4c3ac702403c554` |
| rows | 与各产物行数一致 |

`taxonomy_gate_check.tsv`：8 个数据集的 `dual_ref_conflict` 均为 `n=0 / expected=0 / PASS`，`species_dual_ref_conflict` 亦 `0 / 0 / PASS`。附带注意到各集有少量 `single_ref_conflict`（如 onekp 300、barbarum 8）标为 `unfixable_by_design / INFO`，非 PASS 项。

### 5.2 schema / 行数 / contig 集合（新 vs 旧）

| 数据集 | 行数(新=旧) | 列数(新=旧) | 表头相同 | 含 `Nucleic_acid` | contig 集合相同 | id 唯一 |
|---|---|---|---|---|---|---|
| Alternaria | 1728=1728 | 20=20 | 是 | 否 | 是 | 是 |
| amarum | 957=957 | 20=20 | 是 | 否 | 是 | 是 |
| Aphis | 2448=2448 | 20=20 | 是 | 否 | 是 | 是 |
| barbarum | 20892=20892 | 20=20 | 是 | 否 | 是 | 是 |
| chinense | 11007=11007 | 20=20 | 是 | 否 | 是 | 是 |
| Fusarium | 1637=1637 | 20=20 | 是 | 否 | 是 | 是 |
| onekp | 535103=535103 | 20=20 | 是 | 否 | 是 | 是 |
| ruthenicum | 15181=15181 | 20=20 | 是 | 否 | 是 | 是 |

- 表头 20 列：`contig_id, primary_tool, completeness, confidence, Realm..Species, Realm_agree..Species_agree`，**无 `Nucleic_acid` 列**（符合本轮预期）。
- 新/旧 contig_id 集合**逐个相同**（排序后 diff 为空），且均无重复 id。
- 附注：产物 TSV 的**数据行字段是双引号包裹**的（表头不裹）。这对 `csv` 读取无碍；若后续用裸 `cut/split` 处理 contig_id 需注意去引号。

---

## 6. 我认为不确定 / 需要你决策的点

1. **R2 在本参照库上的语义偏移（最需要留意）**：`rankedlineage.dmp` 是 self-exclusive 的，导致 R2 对"终端种名"几乎必然弃权。如果 R2 原本是想验证"多父级/同名录下的另外一种合法父级"，那它在本数据上**并没有测到那种情形**（C 情形 0 例）。是否要在报告中保留 R2、以及是否需要用一份 self-inclusive 的参照重跑 R2，请你定。
2. **gate `dual_ref_conflict=0 PASS` 与独立参照下仍有残余矛盾的关系**：gate 判的是管线自己的 `genus_family_ref`，而本报告用 `rankedlineage.dmp`；新 onekp 在后者下仍有 1460 行（0.27%）被判矛盾。两者参照不同，不构成冲突，但若对外表述"零矛盾"需谨慎，建议以"新产物对该参照的矛盾率 0.2–0.5%"表述。
3. **新产物残余错配的处置**：`Adenoviridae + Organic Lake phycodnavirus`（onekp 300 例）与 `Pandoravirus` 种配错属（多数据集）是**成体系的标签错配**，不是随机噪声。是否要回查分类器（diamond_lca 的 genus 赋值）并做一次针对性修正，超出本次只读范围，请你决策。
4. **一个正面旁证**：onekp 旧产物中 `ERR2040122_clean_NODE_23705...`（Genus=Hepacivirus / Species=Hepacivirus Q）在新产物已变为 Genus=Orthohepacivirus、与参照自洽，说明刷新确实修掉了一部分错配；但新产物仍保留了 Organic Lake / Pandoravirus 两类系统错配（见 4.4），建议一并处置。

---

## 7. 复现方式

服务器上（全部只读，临时文件在 `/tmp/indep_verify/`）：

```bash
# 1) 16 份表的 R1/R2（3 并发）
nohup bash /tmp/indep_verify/run_all.sh > /tmp/indep_verify/run_all.out 2>&1 &
# 2) R1 与既有 chimera.txt 逐格比对
python3 /tmp/indep_verify/compare_chimera.py
# 3) R1/R2 分歧明细 + 逐例归因/抽样
nohup bash /tmp/indep_verify/detail_run.sh > /tmp/indep_verify/detail2.out 2>&1 &
# 4) 新产物 R1-only 标签集中度
python3 /tmp/indep_verify/freq_r1only.py > /tmp/indep_verify/freq.out 2>&1
# 5) schema / contig 集合（新 vs 旧）
bash /tmp/indep_verify/schema_check.sh
```
