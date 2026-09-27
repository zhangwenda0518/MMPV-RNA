# barbarum 门槛 v6.7 重跑审计报告

- **数据集 label**：`barbarum`
- **目录**：`/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out`（下称 `$D`）
- **执行器**：`/tmp/refresh_v67.sh`（md5 `d832c4868a2b8298d7bfca44be0dcb18`，未修改）
- **窗口**：START `01:01:13` → DONE `01:04:31`（R 01:02:22 完成，report/rescue 随后）
- **审计性质**：矛盾率、7 对阶元分布、阶元取值变化均由本人独立实现 `/tmp/indep_chimera.py`（不复用 `chimera_multi.py`）重算，并与脚本产物 `/tmp/refresh_v67/barbarum/chimera.txt` 逐格对照。
- **命名提示**：备份目录名 `Votus.integrated.bak_v66f_20260916` 由执行器统一命名，其内容是**刷新前线上成品**（8 月 12 日生成，矛盾 2738 行 / 13.11%），即 `REPORT_chimera_multi_20260830.md` 表里的「线上」列，**不是**该报告里的 `v6.6f` 候选列（465 行 / 2.23%，只存在于 `/tmp/rc66_multi`，从未上线）。

---

## 1. 共识成品

| 项 | 旧（刷新前线上） | 新（v6.7） |
|---|---|---|
| 路径 | `05_Taxonomy/Votus.integrated.bak_v66f_20260916/final_integrated_classification.tsv`（副本 `/tmp/refresh_v67/barbarum/old.tsv`） | `05_Taxonomy/Votus.integrated/final_integrated_classification.tsv` |
| md5 | `12e8f1715922ee967db9b10cb70e523a` | `a3ead5b5bc4fe290a0e9922497d2fc35` |
| 行数 | 20893（数据行 20892 + 表头） | 20893（数据行 20892 + 表头） |
| mtime | 2026-08-12 01:08 | 2026-09-16 01:02 |

- 共同 contig 数 **20892 / 20892**，无增无减；列结构未变（未引入 `Nucleic_acid`，本轮未跑 `--stage analysis`）。

## 2. Gate 自检（原文）

`taxonomy_gate_stamp.tsv` 第 2 行：

```
6.7	/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R	f0683ac6af1a9ee1395a100d128365ca	/home/zhangwenda/database/taxonomy/genus_family_ref.tsv	865283b4b1a3a2f3c4c3ac702403c554	20892	2026-09-16 01:02:21
```

`taxonomy_gate_check.tsv`：

```
rows_total	20892	-	INFO
family_filled	17494	-	INFO
genus_filled	14630	-	INFO
genus_ref_known	13559	-	INFO
dual_ref_consistent	12583	-	INFO
dual_ref_conflict	0	0	PASS
single_ref_conflict	8	unfixable_by_design	INFO
species_dual_ref_conflict	0	0	PASS
species_ref_known	12048	-	INFO
species_single_ref_conflict	10	unfixable_by_design	INFO
```

`r.log` 对应行：`科-属自检通过: 双参照一致冲突 0 条; 单侧冲突 8 条; 双参照一致相容 12583 条`、`科-种自检通过: 双参照互斥 0 条; 单侧冲突 10 条`。

## 3. 矛盾率

判据同参考实现：对每对相邻阶元 (H 高, L 低)，用 `rankedlineage.dmp` 查 L 的定型父级，L 有参照父级且行内 H 有值时才比。barbarum 需查参照 name 3723、命中 3700。

| 口径 | 矛盾行 | 矛盾率 |
|---|---|---|
| 旧 | **2738** | **13.11%** |
| 新 | **87** | **0.42%** |

7 对阶元（旧 → 新）：

| 阶元对 | 旧 | 新 |
|---|---|---|
| Realm–Kingdom | 57 | 27 |
| Kingdom–Phylum | 79 | 2 |
| Phylum–Class | 97 | 13 |
| Class–Order | 536 | 1 |
| Order–Family | 369 | 0 |
| Family–Genus | 1273 | 8 |
| Genus–Species | 1011 | 37 |
| **合计** | **2738** | **87** |

本人独立实现结果：`SUMMARY 20892 2738 13.1055 2585 12.3732`（旧）、`SUMMARY 20892 87 0.4164 53 0.2537`（新），R1 与 `chimera.txt` 两行**逐格一致**（含 7 对分布）。第二列 R2 为另一种参照口径，见末节。

## 4. 阶元取值变化与丢失（本人重算）

| 阶元 | 取值变化行数 |
|---|---|
| Realm | 69 |
| Kingdom | 106 |
| Phylum | 112 |
| Class | 232 |
| Order | 624 |
| Family | 905 |
| Genus | 2248 |
| Species | 3472 |

- **有值 → 空**：Order **397**、Family **425**
- **空 → 有值**：Order 0、Family 0
- `taxonomy_fill_stats.tsv` 印证：Order 18162 → 17765、Family 17919 → 17494、Genus 16022 → 14630、Species 18056 → 17201；Realm/Kingdom/Phylum 基本不动（Phylum 19039 → 19038）。
- Phylum↔Class 参照相容性（旧 → 新）：`both_ok 20783 / fixed 96 / broken 12 / both_bad 1`。即新产物里 13 行 P↔C 矛盾中 12 行是新增的跨阵营拼接（旧产物里它们参照相容），只剩 1 行是旧产物也坏的。这是「半数票也淘汰」把细阶元改成服从多数阵营后，粗阶元阵营残留所致。

## 5. 下游产物

### 10_Reports

`diff -rq $D/10_Reports $D/10_Reports.bak_v66f_20260916` → **19 个文件内容差异 + 2 个仅新有**（合计 21 行 diff），体积 85M → 93M。

- 门槛派生：`final_integrated_classification.tsv`、`taxonomy_composition.tsv`、`taxonomy_fill_stats.tsv`、`taxonomy_sunburst.html`、`classification_sankey.html`、`classification_sankey_plant.html`、`plant_final_taxonomy.tsv`、`plant_virus_summary.tsv`、`All_plant.viruses_info.tsv`、`HQ_plant_viruses_info.tsv`、`all_plant_viruses_genus_summary.tsv`、`upset_data.json`
- rescue 层重生成：`rescue_report.tsv`、`rescue_evidence_scored.tsv`、`keep_summary.tsv`
- 环境/生成器因素：`directory_tree.txt`（纳入新备份目录）、`pipeline_report.html`、`stage_summary.tsv`（路径从已不存在的 `/home/zhangwenda/virus/data-2026/data-test/...` 更正为当前数据集树）、`filter_summary.tsv`（新值全 0，根因见 §7）
- 仅新有：`hmm_ct3_hits.tsv`、`hmm_ct3_rescue_candidates.tsv`

### 08_Rescue / 10_Reports rescue

`rescue.log` 摘要行：

```
候选=1056 A=22(含prepass 22) B=13 C=273 D=0 final=276
新报告 branch 分布: {'fail': 770, 'C': 273, 'B': 13}
```

- 拯救 = B + C = **286 / 1056 = 27.1%**；未拯救 **770**（其中 <2000bp 749、属组失败 21）
- 与旧版对比：合计 308 → 286、未拯救 748 → 770、过短 864 → 749 + 21
- `08_Rescue/Plant/rescue_report.tsv` 与 `rescue_report.tsv.bak_v661` 行数均 1057（1056 数据行），差异集中在 `contig_2246/2638/2790/3456 …` 一类的 `fail/未通过任何分支` → `C/blastn_completo` 翻转，属 rescue 生成器策略变化，不是门槛直接后果
- `10_Reports/rescue_report.tsv` 与 `08_Rescue/Plant/rescue_report.tsv` 内容不一致：`report_pipeline.py` 先 `copy2`（步骤 5）再重生成源（步骤 6），拷贝滞后一步，脚本既定顺序，本轮未改

## 6. 备份

- `05_Taxonomy/Votus.integrated.bak_v66f_20260916/`（旧成品，md5 `12e8f171…`）
- `10_Reports.bak_v66f_20260916/`（旧报告树）
- `08_Rescue/Plant/rescue_report.tsv.bak_v661`（既有基线，非本轮建立）
- 日志：`/tmp/refresh_v67/barbarum/{r.log,report.log,rescue.log,chimera.txt,cmp.txt,old.tsv}`

全程未删任何文件；未改 R/py 管线脚本。

## 7. 需要知会的问题

1. `filter_summary.tsv` 新值全 0：`02b_Filter` 下候选 fasta 是指向 `/home/zhangwenda/goji-virome/01_data/data-test/...` 的断链符号链接，`_count_fasta` 读 0。与门槛无关。
2. 刷新前 `10_Reports` 部分文件是旧根目录 `/home/zhangwenda/virus/data-2026/data-test/...`（已不存在）下的陈旧产物，故 diff 里混有「对齐到当前数据集树」的生成器漂移，不能全算到门槛头上。
3. 未动（待决）：`prevalence_full_table.tsv`、`09_Virome_Analysis/**`、`09b_Analysis_Verify/**`。

## 8. 两种参照口径（R1 / R2）

`/tmp/indep_chimera.py` 同时输出两口径：

- **R1**（按 `tax_name` 在 dmp **首次出现**行取整条谱系）：旧 2738（13.11%）、新 87（0.42%）
- **R2**（按阶元列取多父级集合，任一父级相容即不算矛盾）：旧 **2585**（12.37%）、新 **53**（0.25%）

R1 与既有报告口径一致；R2 数值更低，差异集中在 Family–Genus（旧 1273→1265、新 8→0）与 Genus–Species（旧 1011→838、新 37→11），成因是同一名字在 dmp 中按不同阶元出现（同名多义）或多父级集合里含行内高阶元。方法章用哪一口径待定，倾向 R1（与既有全部报告一致、可复现）。
