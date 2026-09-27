# 8 数据集门槛 v6.7 全量重跑汇总总表

- 生成时间：2026-09-16 02:5x（北京时间）
- 服务器：`zhangwenda@202.119.189.246`
- 官方脚本：`/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R`
  - md5 **`f0683ac6af1a9ee1395a100d128365ca`**（备份 `.bak_thr0499_20260916`）
  - `TAX_GATE_VERSION = "6.7"`，`CASCADE_MIN_SHARE` 默认 `0.5`，淘汰条件 `share >= CASCADE_MIN_SHARE`
- 参照表：`/home/zhangwenda/database/taxonomy/genus_family_ref.tsv`
  - md5 **`865283b4b1a3a2f3c4c3ac702403c554`**（8 数据集全部一致）
- 执行器：`scripts/audit/refresh_v67.sh`（服务器 `/tmp/refresh_v67.sh`，md5 `d832c4868a2b8298d7bfca44be0dcb18`）
- 参照谱系库：`rankedlineage.dmp`（列序 `tax_id|tax_name|species|genus|family|order|class|phylum|kingdom|superkingdom`）

---

## 0. 一句话结论

把级联淘汰门槛从「票首严格过半」改成「票首达到半数即淘汰」后，8 个数据集全部刷新成功：共识成品行数与 contig 集合一行未变，`dual_ref_conflict` / `species_dual_ref_conflict` 全为 0 PASS，级间（参照口径 R1）矛盾行从合计 **78360 / 588953（13.31%）** 降到 **1673 / 588953（0.28%）**，降幅 47 倍，且没有任何一个数据集出现 gate 自检失败或行数漂移。

> 本轮范围原则（最高优先）：只维护**分类链路**。非自身引起的底层破损（断链软链、被外部写坏的候选 fasta、僵尸进程、0 字节空候选）只记录不动刀；自身重跑造成的覆盖一律复原。复原与结案明细见 §6.1。

---

## 1. 8 数据集主表

| label | 数据行 | 旧矛盾行 | 旧率 | 新矛盾行 | 新率 | 降幅 | 旧成品 md5 | 新成品 md5 | gate 戳时间 |
|---|---|---|---|---|---|---|---|---|---|
| Alternaria | 1728 | 356 | 20.60% | 6 | 0.35% | 59× | `78f6677f850b7e7d009d52fb2cb187cd` | `79ff593bc4a06c16f003b89ea0fc9edf` | 01:07:19 |
| amarum | 957 | 167 | 17.45% | 2 | 0.21% | 84× | `8d26262525fe390b07f5e7bf1515606c` | `868e4248bac3f878a20d9b31cbaa29b5` | 01:08:20 |
| Aphis | 2448 | 418 | 17.08% | 7 | 0.29% | 60× | `64491217d3292b794c313cc18525ab2d` | `b72315fe25cddcf76a2ec89778d653e3` | 01:08:12 |
| barbarum | 20892 | 2738 | 13.11% | 87 | 0.42% | 31× | `12e8f1715922ee967db9b10cb70e523a` | `a3ead5b5bc4fe290a0e9922497d2fc35` | 01:02:21 |
| chinense | 11007 | 1258 | 11.43% | 57 | 0.52% | 22× | `da1488531f5e00f031edcebec7da5103` | `977bb21b8f0d745c2c74dae33aa1f60c` | 01:06:37 |
| Fusarium | 1637 | 290 | 17.72% | 3 | 0.18% | 97× | `ff119f84fd1ffdd11c03bbfbf308c621` | `974ebd044386c029f6e09952e9d8e3b8` | 01:06:45 |
| onekp | 535103 | 71295 | 13.32% | 1460 | 0.27% | 49× | `028e85b2aeea9a3977fea8566b6d72c7` | `02a2b8446ee748721e4c5b6c60f79bf0` | 01:18:35 |
| ruthenicum | 15181 | 1838 | 12.11% | 51 | 0.34% | 36× | `8d143dc3a684bd52845dc1d6d21954a4` | `2416cc1374a88c3643288b351e9109a5` | 01:08:28 |
| **合计** | **588953** | **78360** | **13.31%** | **1673** | **0.28%** | **47×** | | | |

- 「旧」= 刷新前线上成品（2026-08-12 生成，v6.5/v6.6 世代）；「新」= 2026-09-16 v6.7 重跑产物。
- 每个数据集均满足：行数不变、`contig_id` 集合与顺序不变（集合比对，非仅计数）、无重复、`dual_ref_conflict = 0 PASS`、`species_dual_ref_conflict = 0 PASS`、`script_md5 f0683ac6…`、`ref_md5 865283b4…`、备份目录存在。

### 命名警告（重要）

备份目录名 `Votus.integrated.bak_v66f_20260916` / `10_Reports.bak_v66f_20260916` 里的 `v66f` 是执行器硬编码字符串，**不代表内容版本**。其内容 = 刷新前线上成品（8 月 12 日生成，barbarum 矛盾 2738 行），即 `REPORT_chimera_multi_20260830.md` 的「线上」列，**不是**该报告里的 `v6.6f` 候选列（465 行 / 2.23%，只在 `/tmp/rc66_multi`，从未上线）。

顺带一提，现行 R 脚本头注写的是「barbarum 2.23% -> 0.42%」，比较对象是 v6.6f 候选版本。这个数字对候选基线是对的，但与服务器上真实部署过的成品（13.11% -> 0.42%）不是同一个起点。若论文或报告要引用「改前」数字，需要显式说明基线是哪一版。

---

## 2. 7 对相邻阶元分布（旧 → 新）

单位为「该阶元对上的矛盾出现次数」，同一行可能同时计入多对，故各对之和大于矛盾行数。

| 数据集 | Realm-Kingdom | Kingdom-Phylum | Phylum-Class | Class-Order | Order-Family | Family-Genus | Genus-Species |
|---|---|---|---|---|---|---|---|
| Alternaria | 6 → 1 | 5 → 0 | 4 → 2 | 43 → 2 | 81 → 0 | 167 → 0 | 192 → 2 |
| amarum | 1 → 1 | 1 → 0 | 10 → 1 | 17 → 0 | 64 → 0 | 99 → 0 | 48 → 0 |
| Aphis | 68 → 3 | 4 → 1 | 5 → 0 | 90 → 0 | 185 → 0 | 193 → 0 | 111 → 3 |
| barbarum | 57 → 27 | 79 → 2 | 97 → 13 | 536 → 1 | 369 → 0 | 1273 → 8 | 1011 → 37 |
| chinense | 16 → 14 | 22 → 3 | 28 → 8 | 220 → 3 | 180 → 1 | 691 → 3 | 292 → 27 |
| Fusarium | 2 → 1 | 7 → 0 | 5 → 0 | 67 → 0 | 82 → 0 | 164 → 1 | 99 → 1 |
| onekp | 1177 → 375 | 14586 → 0 | 14276 → 304 | 29935 → 48 | 10195 → 192 | 23612 → 300 | 572 → 291 |
| ruthenicum | 69 → 29 | 42 → 1 | 32 → 8 | 430 → 3 | 177 → 2 | 878 → 2 | 644 → 9 |

两点读数：

1. 降幅最大的是 Kingdom-Phylum 与 Class-Order（onekp 14586 → 0、29935 → 48），这两对是改 `>=` 直接收割的「半数票并列」区。
2. 残余最集中的是 onekp 与 barbarum 的 Realm-Kingdom、Phylum-Class、Genus-Species，形态是「单工具自成一派」或「参照库口径与工具体系不一致」，不是票数能消掉的部分。

---

## 3. 阶元取值变化行数（旧 → 新）

| 数据集 | Realm | Kingdom | Phylum | Class | Order | Family | Genus | Species |
|---|---|---|---|---|---|---|---|---|
| Alternaria | 59 | 60 | 67 | 71 | 92 | 149 | 216 | 349 |
| amarum | 3 | 3 | 5 | 13 | 37 | 103 | 71 | 145 |
| Aphis | 20 | 51 | 55 | 80 | 158 | 309 | 393 | 460 |
| barbarum | 69 | 106 | 112 | 232 | 624 | 905 | 2248 | 3472 |
| chinense | 100 | 110 | 120 | 120 | 308 | 394 | 1064 | 1811 |
| Fusarium | 48 | 49 | 47 | 56 | 106 | 164 | 216 | 281 |
| onekp | 18399 | 18688 | 5037 | 18380 | 40262 | 39530 | 60148 | 68985 |
| ruthenicum | 100 | 69 | 68 | 124 | 468 | 618 | 1467 | 2446 |

### Order / Family 信息量净损失

| 数据集 | Order 有值→空 | Family 有值→空 | Order 空→有值 | Family 空→有值 |
|---|---|---|---|---|
| Alternaria | 64 | 58 | 0 | 0 |
| amarum | 16 | 11 | 0 | 0 |
| Aphis | 60 | 44 | 0 | 0 |
| barbarum | 397 | 425 | 0 | 0 |
| chinense | 176 | 179 | 0 | 0 |
| Fusarium | 61 | 50 | 0 | 0 |
| onekp | 23184 | 25995 | 0 | 0 |
| ruthenicum | 294 | 275 | 0 | 0 |

这是改 `>=` 的**已知代价**：票首恰好半数时，异议工具在更细阶元的取值一并让位，于是出现「细阶元原本有值、现在留空」的行，且没有一行反向回填。方向是拿信息量换自洽性，属口径选择。

### 各数据集 `taxonomy_fill_stats.tsv` 抽查（barbarum）

| Rank | 旧 Filled | 旧率 | 新 Filled | 新率 |
|---|---|---|---|---|
| Realm | 19015 | 0.9102 | 19015 | 0.9102 |
| Kingdom | 19142 | 0.9162 | 19142 | 0.9162 |
| Phylum | 19039 | 0.9113 | 19038 | 0.9113 |
| Class | 19459 | 0.9314 | 19450 | 0.9310 |
| Order | 18162 | 0.8693 | 17765 | 0.8503 |
| Family | 17919 | 0.8577 | 17494 | 0.8374 |
| Genus | 16022 | 0.7669 | 14630 | 0.7003 |
| Species | 18056 | 0.8643 | 17201 | 0.8233 |

---

## 4. Phylum↔Class 参照相容性（仅 barbarum 有 both_bad / onekp 有较大 broken）

| 数据集 | both_ok | fixed | broken | both_bad |
|---|---|---|---|---|
| Alternaria | 1722 | 4 | 2 | 0 |
| amarum | 946 | 10 | 1 | 0 |
| Aphis | 2443 | 5 | 0 | 0 |
| barbarum | 20783 | 96 | 12 | 1 |
| chinense | 10971 | 28 | 8 | 0 |
| Fusarium | 1632 | 5 | 0 | 0 |
| onekp | 520758 | 14041 | 69 | 235 |
| ruthenicum | 15141 | 32 | 8 | 0 |

`fixed` = 旧产物与参照冲突、新产物相容；`broken` = 相反。barbarum 上 `96 fixed / 12 broken` 与 §2 的 97 → 13 完全自洽：13 行残余 P↔C 矛盾里 12 行是改门槛后新产生的（细阶元服从多数阵营、粗阶元残留另一阵营），1 行旧产物也坏。收益 96 行、代价 12 行。

---

## 5. rescue 层（08_Rescue/Plant，全 8 数据集）

| 数据集 | 候选 | A（含 prepass） | B | C | D | B+C | final（去冗余） | fail |
|---|---|---|---|---|---|---|---|---|
| Alternaria | 33 | 0 | 2 | 26 | 0 | 28 | 27 | 5 |
| amarum | 76 | 0 | 0 | 8 | 0 | 8 | 8 | 68 |
| Aphis | 50 | 2 | 1 | 27 | 0 | 28 | 28 | 22 |
| barbarum | 1056 | 22 | 13 | 273 | 0 | 286 | 276 | 770 |
| chinense | 542 | 9 | 7 | 117 | 0 | 124 | 122 | 418 |
| Fusarium | 50 | 0 | 2 | 34 | 0 | 36 | 36 | 14 |
| onekp | 20517 | 130 | 50 | 9212 | 0 | 9262 | 8717 | 11255 |
| ruthenicum | 1271 | 28 | 21 | 274 | 0 | 295 | 272 | 976 |
| **合计** | **23595** | **191** | **96** | **9971** | **0** | **10067** | **9486** | **12828** |

- `fail` 的构成主要是 `<2000bp`（过短）与少数属组失败，各数据集 `rescue_summary.md` 已按「占未拯救」重算分母。
- **rescue 数字与门槛不构成因果关系**：多条 `fail` → `C / blastn_completo / BLASTN qcov≥98% 直达` 的翻转来自 rescue 报告生成器（`regen_rescue_reports.py` / `report_pipeline.py`）自身的策略修复（C 分支 Completo 直达未计数、A 分支把 prepass 混入分母），本轮只是顺带重跑。Aphis 隔离出步骤 6 的净效应只有 6 行。
- 论文正文 §3.5 现写 *L. barbarum* 284 contigs、*L. ruthenicum* 294 contigs、*L. chinense* 130 contigs，与本次运行产物（286/295/124，或去冗余 276/272/122）三者互不相等，需要按定稿口径重取，见 §8 第 2 条。

### 5.1 救援计数四套口径的出处（2026-09-16 实测）

| 口径 | barbarum | ruthenicum | chinense | 出处 | 含 HMM |
|---|---|---|---|---|---|
| A+B+C（8/19 版口径） | 308 | 323 | 133 | `08_Rescue/Plant/rescue_report.tsv.bak_v661` + `rescue_summary.md.bak_v661` | 否 |
| B+C（候选表内，今晚重生成） | 286 | 295 | 124 | 当前 `08_Rescue/Plant/rescue_report.tsv` / `rescue_summary.md` | 否 |
| KEEP+REVIEW（09b 证据整合） | 227（pre-HMM 208） | 217 | 76 | `09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv` 第 26 列 `verdict` | **仅 barbarum** |
| 去冗余 HQ vOTU | 276 | 272 | 122 | 同上 `rescue_summary.md` | 否 |

- 308→286 / 323→295 / 133→124 的差正好是 A 支（CheckV ≥90%，22/28/9 条）。今晚重生成报告时把 A 移出候选表、单列「免拯救 prepass」，与 HMM 无关。
- `rescue_pipeline.py`（md5 `2960e4fe3496ad9d72b51699d2cd95c3`）全文零处 `hmm` 字样；`rescue_report.tsv` 第 5 列 `method` 只有 `vsi_extended`（B）与 `blastn_completo`（C），故 **286/295/124 不含 HMM 判据**。
- HMM 判据落点在 09b：`integrate_rescue_evidence.py`（mtime 9/15 00:34）五层证据链里，CDD 域是唯一的「门」，CT3 病毒 HMM 全长 profile 是「补票 2」（`hmm_ct3_evidence.py`，mtime 9/14 21:53；默认参与裁决：无 CDD 病毒域但 CT3 有非噬菌体病毒命中时把 DROP 抬到 REVIEW，`--no-hmm-rescue` 可关闭）。产物 `rescue_evidence_scored.tsv` 末列 `hmm_rescue` + 13 列 `ct3_*`。
- **只有 barbarum 跑过 HMM**：`hmm_ct3_hits.tsv` 310 命中（三库 RDRP_HMMs / Useful_Annotation_HMMs / Virion_HMMs），42 条被抬（19 `DROP->KEEP` + 23 `REVIEW->KEEP`），KEEP 144→186、DROP 83→64、REVIEW 64→41。
- ruthenicum / chinense 的 `rescue_evidence_scored.tsv` 与 8/25 备份逐字节一致（92322 / 39814 B），`hmm_rescue` 列全空，且两数据集没有 `hmm_ct3_hits.tsv`；HMM 判据尚未在这两个数据集落地。
- 另一路 `hmm_rvdb_evidence.py`（mtime 9/14 17:56）是 RVDB-prot-HMM 的第三路证据，按 LCA 谱系门控分六桶，只认病毒桶，供 9b 与 CDD / blast 并列，不参与 08_Rescue 计数。

---

## 6. 下游产物与已修项

| 项 | 状态 |
|---|---|
| `05_Taxonomy/Votus.integrated/final_integrated_classification.tsv` | 8 数据集全部刷新（新 md5 见 §1） |
| `10_Reports/**` | 8 数据集全部重生成，`diff -rq` 行数（内容差异 + 仅新有）：Alternaria 17、amarum 18、Aphis 17、barbarum 21（19+2）、chinense 17、Fusarium 16、onekp 13、ruthenicum 18 |
| `08_Rescue/Plant/rescue_report.tsv`、`rescue_summary.md` | 8 数据集全部重生成 |
| `rescue_pipeline.py` | 修复，md5 `2960e4fe3496ad9d72b51699d2cd95c3`，8 数据集回填完成 |
| `generate_pipeline_report.py` | Stage 10 元数据区块修复，md5 `12c239956c2e22b22d6c2c6c3eb1ead7` |
| 孤儿副本 `utils/auto_known_virus.py` | 归档 `archive/auto_known_virus.orphan_utils_20260915.py` |
| `REPORT_chimera_multi_20260830.md` | 已加版本状态批注、§6.5 两处加取代标记、新增 §8（v6.7 定稿口径与 8 数据集结果） |

10_Reports 的变化分三类，不能全算到门槛头上：

1. 门槛直接派生：`final_integrated_classification.tsv`、`taxonomy_composition.tsv`、`taxonomy_fill_stats.tsv`、`taxonomy_sunburst.html`、两张 sankey、`plant_final_taxonomy.tsv`、`plant_virus_summary.tsv`、`All_plant.viruses_info.tsv`、`HQ_plant_viruses_info.tsv`、`all_plant_viruses_genus_summary.tsv`、`upset_data.json`
2. rescue 层重生成：`rescue_report.tsv`、`rescue_evidence_scored.tsv`、`keep_summary.tsv`
3. 环境 / 生成器：`directory_tree.txt`（纳入新备份目录）、`pipeline_report.html`、`stage_summary.tsv`（路径从已不存在的 `/home/zhangwenda/virus/data-2026/data-test/...` 更正为当前数据集树）、`filter_summary.tsv`（新值全 0，属本轮重跑打坏，**已复原**，见 §6.1）

另有一类不在 diff 里、但确实没刷新的文件，只在 onekp 暴露：`classification_sankey.html`、`classification_sankey_plant.html`、`plant_final_taxonomy.tsv` 三者 mtime 仍是 `2026-09-12 15:41`。根因写在 onekp 的 `report.log` 里：`[WARN] Sankey 生成失败: ... taxonomic_sankey.py ... timed out after 120 seconds`。536k 行的输入在 120 秒硬超时内跑不完，脚本降级跳过；goji 的 7 个数据集体量小，sankey 正常重生成。详见 `onekp.md` §6.5。该三文件已于本轮补跑回填，见 §6.1。

### 6.1 本轮「自己造成的覆盖」复原记录

范围原则（本轮最高优先）：只维护分类链路；非自身引起的底层破损只记录不动刀；自身重跑造成的覆盖一律复原。

**（1）`10_Reports/filter_summary.tsv` 全 0 退化已复原**

- 成因：`02b_Filter` 候选 fasta 是指向 `/home/zhangwenda/goji-virome/01_data/data-test/...`（`goji-virome` 整个目录已不存在）的断链符号链接，`_count_fasta`（`report_pipeline.py` L65-68）用 `os.path.isfile` 跟随链接读到 0；写行条件 `if fm == 'comb' or nf > 0`（L288-307）于是每样本只剩一行 `comb=0`（barbarum 1179 → 397 行、全 0、9731B）。
- 处理：先 `cp -a` 把退化版备份为 `filter_summary.tsv.bak_v67zero_20260916`，再从 `10_Reports.bak_v66f_20260916/` 复原；复原后 md5 与备份逐一对齐：

| 数据集 | 复原后 md5（= 备份 md5） | 行数 |
|---|---|---|
| Alternaria | `594d1f200b29a5729d8e905f447d99bb` | 40 |
| amarum | `87f9468e322970f4c8fda4c586b719aa` | 10 |
| Aphis | `280a99f7f87e5458ea16a3420ca44172` | 22 |
| barbarum | `ab48df363db2327af6c1d88df02fabad` | 1179 |
| chinense | `9ecc2915ecbbb74168bd0aa562f006fc` | 133 |
| Fusarium | `e0cebcd6fe2731cfe11619cff1813fe8` | 19 |
| onekp | `0a004bc472c018501d911c728b4c08e2` | 4027 |
| ruthenicum | `a133b12b871cc361175a5c09d070e87d` | 1145 |

- 断链本身（`02b_Filter` 1560 条 + `00b_HostDepletion` 849 条，`repair_symlinks.py` dry-run manifest 实测）按范围原则**不修**，故该文件在下次全量重跑时会再次退化，属已知行为而非本次遗漏。

**（2）onekp 三个分类产物已回填 `10_Reports/`**

- 源 `/tmp/sankey_fix/out/`，日志 `/tmp/sankey_fix/logs/onekp.log` 实测 `EXIT_ALL=0 ELAPSED_ALL=122`（旧硬超时 120s 只差 2 秒）、`EXIT_PLANT=0 ELAPSED_PLANT=5`、`ELAPSED_BUILD=4`、`plant_ids 20054`；旧版（9/12 15:41）备份为 `.bak_v67presankey_20260916`，回填后源目 md5 逐一对齐：

| 文件 | 回填后 md5 | 字节 | 旧字节 |
|---|---|---|---|
| `classification_sankey.html` | `6ac308fe1bd05a653044a473771d8861` | 4947250 | 5079507 |
| `classification_sankey_plant.html` | `45853abbcdda7512f430af0b813a9e63` | 4888776 | 4901797 |
| `plant_final_taxonomy.tsv` | `804f90650e2e7db18d49dc2751791d0b` | 6998984 | 7046137 |

- 配套补丁（`report_pipeline.py`，线上 md5 `5e51804253d0ca5423ddaecff2181122`，原版备份 `report_pipeline.py.bak_sankey_20260916`）：`SANKEY_TIMEOUT_DEFAULT = 1800` + `_sankey_timeout()` + `_run_sankey()`（L1011-1023）、filter_summary 全 0 守卫 WARN、rescue 拷贝改确定性优先 `Plant`（L702-708）。

**（3）两项结案（非 v6.7 不一致）**

- onekp `upset_data.json`：生成者是 R 脚本 VIS 段（`logs/r_consensus.log` 有 `导出 upset_data.json (549934 元素)`），`report_pipeline.py` 仅 `shutil.copy2` 拷贝。该文件内容为工具投票交互集，与共识门槛无关；本次 mtime `2026-09-16 01:08`，`05_Taxonomy/Votus.integrated/`、`10_Reports/` 与备份三处 md5 同为 `ca59b07714cfc5e875d066cd3e4100c2`，逐字节相同属正常（barbarum 新旧不同 `d31f068c...` / `be63de0d...`，为数据集间差异）。
- `keep_summary.tsv`：源 `09b_Analysis_Verify/keep_summary.tsv` mtime `2026-09-15 00:59`、49569B，早于本轮刷新（备份那份 8/25 19:14、50192B），属 9/15 那批 `_kill_20260915` 工作，与 v6.7 门槛无关。

### 环境层已知问题（与门槛无关，未修）

1. `02b_Filter` / `00b_HostDepletion` 断链符号链接（共 2409 条，全指向已删除的 `goji-virome/01_data/data-test/...`）：按范围原则不修，只影响 `filter_summary.tsv` 一类报告层计数，本轮已用复原方式回避。
2. 刷新前 `10_Reports` 有部分文件来自旧根目录 `/home/zhangwenda/virus/data-2026/data-test/...`（已不存在），所以 diff 里混着「对齐到当前数据集树」的生成器漂移。
3. `10_Reports/rescue_report.tsv` 与 `08_Rescue/Plant/rescue_report.tsv` 内容不一致：`report_pipeline.py` 步骤 5 的 `copy2` 早于步骤 6 的重生成，拷贝滞后一步（脚本既定顺序，本轮未改）。

---

## 7. 参照口径 R1 / R2

独立实现 `/tmp/indep_chimera.py`（md5 `73315578d6160ea6359ac9e7b7cd23cd`，121 行，自写 csv 解析与两套索引，全文件无管线模块 import）同时输出两口径；参照库 `rankedlineage.dmp`（md5 `f6b86ae7000837d638dd18035c7c41d0`，2,831,927 行）。

- **R1**：按 `tax_name` 在 dmp **首次出现**行取整条谱系，比对该低阶名在该高阶层元的祖先值与行内高阶元。与既有全部报告同口径。
- **R2**：按 `(低阶元, 值)` 聚合所有该值出现在该低阶元列的行的高阶取值集合，行内高阶值不在集合即矛盾。

### 7.1 两口径 8 数据集全表

| 数据集 | 行数 | R1 旧 | R1 新 | R2 旧 | R2 新 |
|---|---|---|---|---|---|
| Alternaria | 1728 | 356 (20.60%) | 6 (0.35%) | 335 (19.39%) | 5 (0.29%) |
| amarum | 957 | 167 (17.45%) | 2 (0.21%) | 166 (17.35%) | 2 (0.21%) |
| Aphis | 2448 | 418 (17.08%) | 7 (0.29%) | 411 (16.79%) | 5 (0.20%) |
| barbarum | 20892 | 2738 (13.11%) | 87 (0.42%) | 2585 (12.37%) | 53 (0.25%) |
| chinense | 11007 | 1258 (11.43%) | 57 (0.52%) | 1237 (11.24%) | 34 (0.31%) |
| Fusarium | 1637 | 290 (17.72%) | 3 (0.18%) | 289 (17.65%) | 1 (0.06%) |
| onekp | 535103 | 71295 (13.32%) | 1460 (0.27%) | 70488 (13.17%) | 915 (0.17%) |
| ruthenicum | 15181 | 1838 (12.11%) | 51 (0.34%) | 1805 (11.89%) | 46 (0.30%) |

旧 11.43%–20.60% → 新 0.18%–0.52%（R1）。两个口径给出的量级与方向完全一致，主结论对口径选择稳健。

### 7.2 R1 与既有 `chimera.txt` 逐格一致

自动比对 8 数据集 × (新产物, 旧产物) = 16 行，每行的行数、7 对阶元计数、矛盾行数全部比对，结果 `ALL_CELLS_MATCH`。这是本次门槛变更最硬的一条一致性证据。

### 7.3 R2 是 R1 的严格子集，分歧只有一种成因

16 份表里 **R2-only 恒为 0**，即 R2 ⊆ R1；分歧全部是 R1-only，且**只出现在最低两对阶元**（Family-Genus、Genus-Species），其余 5 对两口径逐格相同。

R1-only 的成因归因（`detail.py`）在 16 个文件里**全部落在同一类（原因 A）**，任务设想的「多父级集合包含行内高阶元」一例都没出现。原因 A 的机制：R1 用名字索引，名字在库里出现即可定祖先；R2 用阶元列索引，只有当该名字确实以该低阶元身份出现在库里才有键，对「错阶名」与「终端名」**直接弃权**。所以 R2 的「低」来自弃权，而非接受了另一个合法父级。

叠加一个参照库侧的事实：本 `rankedlineage.dmp` 是 **self-exclusive（自排除）** 形态，任何分类单元自身所在阶元列为空，谱系列只写严格高于自身的祖先（`Homo sapiens` 行 species 空、genus=Homo；species 列非空仅 290,168 行，都是种下单元）。于是叶子种名在 R2 里几乎注定无键。

### 7.4 残余矛盾不是弥散噪声，被少数标签成体系贡献

| 数据集 | R1-only（对-次） | 最主要标签 |
|---|---|---|
| onekp | 548 | `Family=Adenoviridae + Genus=Organic Lake phycodnavirus` **300 次**；其后是 `Pandoravirus` 属的种配错属 |
| barbarum | 34 | `Species=Pandoravirus salinus + Genus=Alphahydrivirus` 11；`Adenoviridae+Organic Lake` 8 |
| chinense | 23 | `Pandoravirus salinus+Alphahydrivirus` 7；`Adenoviridae+Organic Lake` 3 |
| ruthenicum 5 / Fusarium 2 / Aphis 2 / Alternaria 1 / amarum 0 | | 同型 |

参照库落位已核：`Organic Lake phycodnavirus`（tax_id 938083）谱系为 `Phycodnaviridae / Algavirales / Megaviricetes / Nucleocytoviricota / Bamfordvirae / Varidnaviria`，**是种级名，不在 Adenoviridae**，行里却配 `Family=Adenoviridae`，是真的跨谱系拼接；`Pandoravirus inopinatum`（tax_id 1605721）属为 `Pandoravirus`，行里却配 `Genus=Medusavirus / Mastadenovirus / Theiavirus / Moumouvirus` 等。

一句话：残余矛盾基本来自两类反复出现的错配，一是分类器把藻病毒当属却配了腺病毒科，二是把 `Pandoravirus` 属的种配到其它病毒属。

### 7.5 口径建议

**建议方法正文用 R1 为主口径**，理由：

1. R1 与管线自身的对照实现（`chimera_multi.py`）同源，`taxonomy_gate_check.tsv` 的 gate 也按这一路判据建设，口径自洽；
2. R1 对可解析的名字从不弃权，度量的是「该行低阶名在参照中的归位是否与行内高阶元自洽」，更贴近「分类矛盾」想表达的东西；
3. R2 更保守，但其「少」在本文数据上主要来自弃权（R2-only = 0，100% 原因 A），把「名字不在该阶元或无子代」跳过，并非真正接受另一合法父级。把这样一个系统性偏保守、还叠加了 self-exclusive 参照效应的数字当「更干净的结果」对外报，审稿人按名字口径复算会得到更高数字，风险较大。

若审稿人质疑同名多义或多父级，可补报 R2 作敏感性分析，但需同时说明弃权机制。

### 7.6 残余错误的可验证清单（2026-09-16 复核）

**（1）8 数据集现状重新点验**：gate 戳全部 `6.7`、`script_md5 f0683ac6af1a9ee1395a100d128365ca`、`ref_md5 865283b4b1a3a2f3c4c3ac702403c554`，`dual_ref_conflict` 与 `species_dual_ref_conflict` 全 0 PASS，成品 md5 与 §1 表逐一对齐（Alternaria `79ff593b…`、Aphis `b72315fe…`、Fusarium `974ebd04…`、amarum `868e4248…`、barbarum `a3ead5b5…`、chinense `977bb21b…`、ruthenicum `2416cc13…`、onekp `02a2b844…`）。门槛引发的工具间矛盾属闭环状态。

**（2）`single_ref_conflict` 是真错配，不是「不可判」**。判据见 `virus_classifier_analysis.R` L508-512：

```r
n_single_conf <- cur[dual == FALSE & fam_l != fam_l_ref, .N]
add("single_ref_conflict", n_single_conf, "unfixable_by_design", "INFO")
```

`dual == FALSE` = 该属在参照表里只有 NCBI 或 VMR 单侧来源，此时产物 `Family` 与那唯一来源不一致即计入，但状态是 INFO 而非 FAIL。8 数据集实测：

| 数据集 | single_ref_conflict | species_single_ref_conflict |
|---|---|---|
| onekp | **300** | **223** |
| barbarum | 8 | 10 |
| chinense | 3 | 4 |
| ruthenicum | 2 | 4 |
| Fusarium | 1 | 1 |
| Aphis | 0 | 1 |
| amarum | 0 | 1 |
| Alternaria | 0 | 0 |
| **合计** | **314** | **244** |

onekp 那 300 行经逐行实测：`Family="Adenoviridae"` + `Genus="Organic Lake phycodnavirus"`，300/300 全部如此。而参照表 `genus_family_ref.tsv` 该属写 `Phycodnaviridae`（`VMR_Family` 空，故归单源桶），NCBI `rankedlineage.dmp` tax_id 938083 谱系为 `Phycodnaviridae / Algavirales / Megaviricetes / Nucleocytoviricota / Bamfordvirae / Varidnaviria`。参照侧两处一致、产物错，性质上可修，只是 gate 把它放在 INFO 桶里放行了。同型行散布：barbarum 8、chinense 3、ruthenicum 2、Fusarium 1。

**（3）跨阶元启发式检查（新增脚本，只读）**：`scripts/audit/rank_consistency_check.py`（服务器 `/tmp/rank_consistency_check.py`，md5 `970e5b75823eb453e13b70d53e68de9e`）比对「Species 首词」与 `Genus` 列：

| 数据集 | total | Genus+Species 皆非空 | cross_rank_mismatch | 占比 | Genus 列含空格（物种级名） |
|---|---|---|---|---|---|
| Alternaria | 1728 | 1111 | 55 | 4.95% | 0 |
| Aphis | 2448 | 1541 | 111 | 7.20% | 23 |
| Fusarium | 1637 | 1040 | 63 | 6.06% | 1 |
| amarum | 957 | 621 | 56 | 9.02% | 2 |
| barbarum | 20892 | 13448 | 1025 | 7.62% | 66 |
| chinense | 11007 | 7185 | 694 | 9.66% | 20 |
| ruthenicum | 15181 | 10286 | 1181 | 11.48% | 70 |
| onekp | 535103 | 202279 | 17454 | 8.63% | 1060 |

**这张表的占比不能当错误率用**：物种名本就常以宿主 / 地点 / 形态开头，`Avian leukosis virus` 属 `Alpharetrovirus`、`Paramecium bursaria Chlorella virus 1` 属 `Chlorovirus`、`Bathycoccus sp. virus` 属 `Prasinovirus` 都是合法归位，全部被启发式误报。它真正指出的是一类阶元归属问题：`Genus` 列被填进物种级或株级名字，onekp 1060 行（`Cotesia sesamiae bracovirus` 475、`Organic Lake phycodnavirus` 300、`Cotesia vestalis bracovirus` 123、`Human papillomavirus` 122 …）、barbarum 66、ruthenicum 70、chinense 20、Aphis 23、amarum 2、Fusarium 1。这类与 `Organic Lake` 同源（名字落在错的阶元），但成因是参照表无该属记录时工具回填自由文本，不是参照冲突。

**（4）归纳**：门槛解决的是「工具之间对不上」，残余里能继续压低的是「名字落错阶元 / 单源参照被忽略」这一族，量级 314 行（family 级）+ 244 行（species 级）≈ 0.09%，集中在 onekp。修与不修要大王定，见 §9 第 12 项。

---

## 8. 论文正文同步建议

目标文件：`MMPV-RNA-paper/MMPV_Lycium_virome_manuscript_EN_v2_citations.md`（v1 与 v2 需一并处理）。

### 8.1 必须改：共识引擎描述（三份稿件各一处）

现状（v2_citations §2.3，v1 第 54 行、v2 第 56 行同款）：

> ...integrated through a **weighted-voting consensus engine in R** across the eight-rank ICTV hierarchy.

这句与现行实现冲突。引擎自 v6.5 起已改为**计数制**（每工具每阶元一票，`MMPV_VOTE_WEIGHT=weighted` 仅作回退），v6.6 起平票裁决显式两级化，v6.7 起门槛为 `>= 半数`。建议替换为：

> Taxonomy was assigned by an eight-tool ensemble (geNomad, Metabuli, CAT/BAT, DIAMOND LCA, VITAP, MMseqs2, ACVirus and vConTACT3) integrated by a **count-based cascade consensus implemented in R** across the eight-rank ICTV hierarchy. At each rank, proceeding from Realm to Species and counting one vote per tool per rank, the modal value among tools not yet eliminated was retained; when the modal value reached **one half or more of the votes cast at that rank**, tools that voted otherwise were eliminated from all finer ranks, so that finer assignments stay nested within the surviving higher-rank lineage. Unresolved values were never treated as dissent, and an eliminated tool could not re-enter. Ties were resolved first by consistency with the next finer rank and then by a fixed tool order (`ACVirus > CAT > VITAP > diamond_lca > genomad > metabuli > mmseqs`). The consensus engine version (`6.7`) and the script and reference-table MD5 sums are written to `taxonomy_gate_stamp.tsv` alongside each output.

要点：删掉 `weighted-voting`；写明「每工具每阶元一票」；写明淘汰门槛是「达到半数」；写明 NA 不是异议票；写明版本戳与 md5 可追溯。

### 8.2 必须核对：rescue 数字（§3.5）

现文：*L. barbarum* **284** rescued contigs、*L. ruthenicum* **294**、*L. chinense* **130**。
本次运行产物：barbarum 286（B+C）/ 276（去冗余 final）、ruthenicum 295 / 272、chinense 124 / 122。

三者互不相等，说明正文数字来自更早的 vintage，且「rescued contigs」用的是哪一种计数（拯救条目数 B+C，还是去冗余后的 HQ vOTU 数）未注明。建议：先定口径（一个数只对应一种计数定义），再从刷新后的 `08_Rescue/Plant/rescue_report.tsv` + 去冗余结果重取，然后在正文写清「rescued contigs (deduplicated HQ vOTUs)」之类的限定语。

### 8.3 建议补：可复现性数字

- 方法段可以点名 gate 版本 `6.7` 与 `share >= 1/2` 的淘汰语义，正文加一句即可，不必展开票数例子。
- 8 数据集矛盾率表适合放进补充材料，作为「共识引擎与参照库一致性」的自检证据（13.31% → 0.28%，并说明残余两类成因）。
- 若正文或补充材料出现「改前」矛盾率，必须写清基线是哪一版（真实部署的 13.11%，还是未上线的 v6.6f 候选 2.23%）。

### 8.4 其余文档（不属于论文，但会互相引用）

按 `docs_sweep` 只读横切的结论，风险最高三处：`REPORT_cascade_v64_20260830.md §16.4` 英文方法段（`strict majority` + gate v6.6）、`REPORT_cascade_v64_20260830.md §16.1` 定稿表、`doc/10-taxonomy-host.md` / `doc/METHODS_TEMPLATE.md` / `mmpv-*` skill（仍写「加权投票 / 深度权重」）。`REPORT_chimera_multi_20260830.md` 已处理。`TAX_GATE_REQUIRED="6.2"` 与文件名里含 `499` 的脚本属无需改动项。

---

## 9. 待决项清单

| # | 事项 | 影响 |
|---|---|---|
| 1 | `prevalence_full_table.tsv` 是否重跑 | 可复现调用路径未定位，输入指向 `~/virus/data-2026/data-test` |
| 2 | `09_Virome_Analysis/**`、`09b_Analysis_Verify/**` 是否重跑 | 跑 `--stage analysis` 会就地补 `Nucleic_acid` 列，属超范围 schema 变化 |
| 3 | `05_Taxonomy/Votus.integrated/calibration_20260914/` 是否纳入成品 | 当前仅作参考目录保留 |
| 4 | rescue 数字（284/294/130 vs 286/295/124 vs KEEP+REVIEW 227/217/76 vs 去冗余 276/272/122）取哪一套口径同步论文 | 决定 §3.5 是否重写，四套口径出处见 §5.1 |
| 5 | 方法章用 R1 还是 R2 口径（倾向 R1，见 §7.5） | 决定正文残余矛盾率写 0.42% 还是 0.25%（barbarum） |
| 5b | `Organic Lake phycodnavirus` 被当属并配 `Adenoviridae`（onekp 300 例）、`Pandoravirus` 属的种被配到其它属，是否回查分类器 genus 赋值 | 这是残余矛盾的两大来源，也是唯一能继续压低的实际入口；根因已核到 `single_ref_conflict` 桶（§7.6），参照侧写 `Phycodnaviridae`、产物写 `Adenoviridae` |
| 6 | `REPORT_chimera_multi §6` 的 A/B/C 三方案：A 已由 v6.7 实施，B/C 是否上线 | B/C 会把残余再降一档，但会进一步牺牲信息量 |
| 7 | `10_Reports` 与 `08_Rescue` 的 rescue 报告不同步是否修 | 下游若以 `10_Reports` 为准会读到偏旧一版 |
| 8 | `.taxonomy.ok` / `.analysis.ok` / `.report.ok` / `.rescue.ok` 时间戳未更新 | flag 不代表本次刷新状态，本轮未 touch |
| ~~9~~ | ~~`filter_summary.tsv` 全 0~~ | **已复原**（§6.1），8 数据集 md5 与备份对齐；断链本身按范围原则不修 |
| ~~10~~ | ~~onekp 两个 sankey 与 `plant_final_taxonomy.tsv` 补跑~~ | **已补跑回填**（§6.1），`report_pipeline.py` 超时已从 120s 提到 1800s |
| ~~11~~ | ~~onekp `upset_data.json` 生成器是否被调用~~ | **结案**：生成器（R 脚本 VIS 段）已调用，输出与门槛无关故逐字节相同 |
| 12 | `single_ref_conflict` 314 行（onekp 300 = `Organic Lake phycodnavirus` 配 `Adenoviridae`）是否修：让单源属的 `Family` 采用唯一来源，或至少把 gate 状态从 INFO 升为可修桶并列台账 | 唯一可继续压低的实际入口，修则 5 个数据集需重跑（onekp / barbarum / chinense / ruthenicum / Fusarium） |
| 13 | `Genus` 列含物种级名字（onekp 1060 行、barbarum 66、ruthenicum 70、chinense 20、Aphis 23）是否收紧 | 属阶元归属问题，与 12 同源但成因不同（参照表无属记录时工具回填自由文本） |
| 14 | HMM 判据（`hmm_ct3_evidence.py` + `integrate_rescue_evidence.py --hmm-report`）是否补齐到其余 7 数据集 | 不补齐则「含 HMM 判据」只对 barbarum 成立，方法段与数据不一致 |

---

## 10. 证据与文件索引

### 脚本与报告（本地）

- 执行器：`scripts/audit/refresh_v67.sh`
- 独立口径脚本：`scripts/audit/indep_chimera.py`（用法 `python3 indep_chimera.py <tsv> [dmp]`，输出 R1/R2 + `SUMMARY` 行）
- 逐数据集审计报告：`scripts/audit/refresh_v67/{barbarum,Alternaria,amarum,Aphis,chinense,Fusarium,ruthenicum}.md`
- 独立交叉验证：`scripts/audit/refresh_v67/indep_verify.md`（8 数据集 R1/R2、R1 与 `chimera.txt` 逐格比对、分歧归因与标签集中度、gate 与 schema 核对）
- onekp 专项复核：`scripts/audit/refresh_v67/onekp.md`（成品与 gate 原文、R1/R2、阶元净损失、P↔C、rescue、10_Reports 逐文件四类归因、应刷未刷三文件清单）
- 成因源报告：`scripts/audit/REPORT_chimera_multi_20260830.md`（含 §8）

### 服务器产物与日志

- 每数据集：`<ds>/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv`、同目录 `taxonomy_gate_stamp.tsv`、`taxonomy_gate_check.tsv`
- 备份：`05_Taxonomy/Votus.integrated.bak_v66f_20260916/`、`10_Reports.bak_v66f_20260916/`、`08_Rescue/Plant/*.bak_v661`（`.bak_v661` 为 8-19 既有，非本轮建立）
- 日志与中间产物：`/tmp/refresh_v67/<label>/{manifest.tsv,old.tsv,r.log,report.log,rescue.log,chimera.txt,cmp.txt}`，nohup `/tmp/refresh_v67_<label>.nohup`

### 复核命令（可原样重跑）

```
md5sum /home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R
# 期望 f0683ac6af1a9ee1395a100d128365ca

D=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out
cat $D/05_Taxonomy/Votus.integrated/taxonomy_gate_stamp.tsv
cat $D/05_Taxonomy/Votus.integrated/taxonomy_gate_check.tsv

python3 /tmp/indep_chimera.py $D/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv
# 期望 SUMMARY 20892 87 0.4164 53 0.2537
python3 /tmp/indep_chimera.py $D/05_Taxonomy/Votus.integrated.bak_v66f_20260916/final_integrated_classification.tsv
# 期望 SUMMARY 20892 2738 13.1055 2585 12.3732
```

### 本轮全程未动

未删任何文件、未改 R/py 管线脚本逻辑（除 `rescue_pipeline.py`、`generate_pipeline_report.py` 两处已申报修复）、未跑 `--stage analysis`、未触碰 `prevalence_full_table.tsv`、`09/**`、`09b/**`、`.ok` flag。
