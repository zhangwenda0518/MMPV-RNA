# refresh_v67 重跑与独立复核报告 — Fusarium

- 生成时间：2026-09-16 01:2x（服务器时间，+0800）
- 执行者：Hana（在大王授权下，经 ssh 到 202.119.189.246）
- 任务性质：线上 R 共识脚本改为「半数票也淘汰」（命中判定 `share >= 0.5`，默认门槛仍 0.5，`TAX_GATE_VERSION = "6.7"`）后，重跑单数据集共识产物并刷新派生报告，随后**独立重算复核**。

## 1. 数据集与目录

| 项 | 值 |
|---|---|
| label | `Fusarium` |
| 目录 D | `/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Fusarium_nematophilum_out` |
| 输入 | `$D/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv`（750367 B，2026-08-12 02:50） |
| 成品 | `$D/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv` |
| R 脚本 | `/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R`，md5 `f0683ac6af1a9ee1395a100d128365ca`（与任务给定一致 ✓） |
| 执行器 | `/tmp/refresh_v67.sh`，md5 `d832c4868a2b8298d7bfca44be0dcb18`（未改动） |

## 2. 指纹对比（旧 / 新）

| | 旧（v6.6f） | 新（v6.7） |
|---|---|---|
| md5 | `ff119f84fd1ffdd11c03bbfbf308c621` | `974ebd044386c029f6e09952e9d8e3b8` |
| 行数（含表头 wc -l） | 1638 | 1638 |
| 数据行 | 1637 | 1637 |
| contig_id 唯一数 | 1637 | 1637 |

旧 md5 的独立来源有两处，互相印证：`$D/05_Taxonomy/Votus.integrated.bak_v66f_20260916/final_integrated_classification.tsv` 与 `/tmp/refresh_v67/Fusarium/old.tsv`，均等于 `ff119f84…`。

## 3. gate 自检（`05_Taxonomy/Votus.integrated/`）

`taxonomy_gate_stamp.tsv` 原文：

```
gate_version	script_path	script_md5	ref_file	ref_md5	rows	timestamp
6.7	/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R	f0683ac6af1a9ee1395a100d128365ca	/home/zhangwenda/database/taxonomy/genus_family_ref.tsv	865283b4b1a3a2f3c4c3ac702403c554	1637	2026-09-16 01:06:45
```

`taxonomy_gate_check.tsv` 原文（节选，含关键两项）：

```
check	n	expected	status
...
dual_ref_conflict	0	0	PASS
single_ref_conflict	1	unfixable_by_design	INFO
species_dual_ref_conflict	0	0	PASS
species_single_ref_conflict	1	unfixable_by_design	INFO
```

结论：版本 6.7 ✓，script_md5 与给定一致 ✓，`dual_ref_conflict = PASS` ✓（`species_dual_ref_conflict` 亦 PASS）。

## 4. 行数与集合（自写 python 用集合比对，非只比计数）

```
行数相等: True (old 1637 / new 1637)
contig_id 集合完全一致: True
仅旧有 0 个: []
仅新有 0 个: []
```

即：**行数不变，contig_id 集合逐元素完全一致**，无增无删。

## 5. 矛盾率（独立重算，判据对齐 `/tmp/chimera_multi.py`）

判据：对每对相邻阶元 (H 高, L 低)，用 `rankedlineage.dmp` 查 L 的定型父级（列序 `tax_id|tax_name|species|genus|family|order|class|phylum|kingdom|superkingdom`，即 `Realm→superkingdom`、`Kingdom→kingdom`），低阶元有参照父级且行内高阶元有值时才比，不等即矛盾。参照库命中 name 919 / 929。

| 口径 | 行数 | Realm-Kingdom | Kingdom-Phylum | Phylum-Class | Class-Order | Order-Family | Family-Genus | Genus-Species | 矛盾行 | 矛盾率 |
|---|---|---|---|---|---|---|---|---|---|---|
| OLD (v6.6f) | 1637 | 2 | 7 | 5 | 67 | 82 | 164 | 99 | **290** | **17.72%** |
| NEW (v6.7) | 1637 | 1 | 0 | 0 | 0 | 0 | 1 | 1 | **3** | **0.18%** |

与脚本侧 `/tmp/refresh_v67/Fusarium/chimera.txt` 对照：**逐项完全一致，无出入**。改前 17.72% 与任务给定锚点吻合；改后 0.18%（低于「1% 上下」的预期，方向一致、更彻底）。

## 6. 阶元取值变化行数 / Order-Family 有值变空

共同 contig 1637（与 `/tmp/refresh_v67/Fusarium/cmp.txt` 逐项一致）：

| 阶元 | 变化行数 |
|---|---|
| Realm | 48 |
| Kingdom | 49 |
| Phylum | 47 |
| Class | 56 |
| Order | 106 |
| Family | 164 |
| Genus | 216 |
| Species | 281 |

- Order 有值→空：**61**；Family 有值→空：**50**
- Order 空→有值：0；Family 空→有值：0
- Phylum↔Class 参照相容性：both_ok 1632，fixed 5，broken 0，both_bad 0

变化集中在 Genus/Species（末端更严命中），高阶元变化较少；信息只减不增（无「空→有值」）。

## 7. 10_Reports 变化清单（对照 `$D/10_Reports.bak_v66f_20260916`）

### 7a. 内容变化的文件（16 个）

**本轮 R/report_pipeline 直接重写，mtime = 2026-09-16 01:06–01:07（12 个）**

`classification_sankey.html`、`classification_sankey_plant.html`、`directory_tree.txt`、`filter_summary.tsv`、`final_integrated_classification.tsv`、`pipeline_report.html`、`plant_virus_summary.tsv`、`stage_summary.tsv`、`taxonomy_composition.tsv`、`taxonomy_fill_stats.tsv`、`taxonomy_sunburst.html`、`upset_data.json`

**内容变化但保留旧 mtime（report_pipeline 用 `shutil.copy2` 从上游阶段拷入，非本轮共识驱动；4 个）**

| 文件 | 备份侧 mtime | 新侧 mtime | 拷入来源 |
|---|---|---|---|
| all_plant_viruses_genus_summary.tsv | 2026-09-15 00:59:08 | 2026-09-15 14:26:11 | `09_Virome_Analysis/all_plant_analysis/` |
| All_plant.viruses_info.tsv | 2026-09-15 00:59:03 | 2026-09-15 14:26:11 | `09_Virome_Analysis/all_plant_analysis/` |
| HQ_plant_viruses_info.tsv | 2026-09-15 00:59:03 | 2026-09-15 14:26:11 | `09_Virome_Analysis/HQ_analysis/` |
| rescue_report.tsv | 2026-09-15 00:59:03 | 2026-09-15 23:49:43 | `08_Rescue/Plant/`（本步骤早于 rescue 重跑，故拷入的是上一代） |

> 说明：这 4 个是 report_pipeline 第 744–768 行 `copy_map` / 第 692 行 `shutil.copy2` 的拷贝行为；`copy2` 保留源文件 mtime，所以 mtime 看着"旧"，内容是上游 09 阶段（09-15 14:26）或 08_Rescue 上一代（09-15 23:49）的。它们与「备份前状态（09-15 00:59）」不同，属既存漂移在本轮被同步，而不是 6.7 共识造成的。

### 7b. mtime 刷新但内容未变（13 个）

`assembly_summary.tsv`、`checkv_confidence.tsv`、`checkv_summary.tsv`、`cluster_pipeline_reduction.tsv`、`cluster_size_distribution.tsv`、`cobra_summary.tsv`、`host_decision_method.tsv`、`host_distribution.tsv`、`ident_summary.tsv`、`plant_final_taxonomy.tsv`、`taxonomy_agreement_stats.tsv`、`taxonomy_consistency_summary.tsv`、`tool_filter_summary.tsv`

（其中 `plant_final_taxonomy.tsv` 在共识变化后内容仍字节不变，值得留意但非异常。）

### 7c. 未触碰（mtime 未变，15 个）

`cdd_evidence_report.tsv`、`chart.min.js`、`frequency_table.tsv`、`host_ictv_classification_summary.tsv`、`host_ictv_confidence.tsv`、`keep_summary.tsv`、`plant_virus_cluster_summary.tsv`、`rescue_detection_summary.tsv`、`rescue_evidence_scored.tsv`、`taxonomy_images.json`、`Viroid.all_info.tsv`、`Viroid.per_sample.tsv`、`Viroid.per_virus.tsv`、`Viroid.species_info.tsv`

无新增文件、无缺失文件（`diff -rq` 无 "Only in"）。

## 8. rescue 报告（`$D/08_Rescue/Plant/`，对照 `.bak_v661`）

| 文件 | 旧（.bak_v661，2026-08-20 05:31） | 新（2026-09-16 01:07:10） |
|---|---|---|
| rescue_report.tsv | 5399 B | 6000 B |
| rescue_summary.md | 753 B | 697 B |

内容实质变化：`rescue_report.tsv` 34 行改动，大量行从 `fail` 变为 `C / blastn_completo / "BLASTN qcov≥98% 直达"`（如 Velarivirus、Potyvirus、Tymovirus、Marafivirus、Deltapartitivirus、Carlavirus 等）。

branch 分布：旧 `{B: 2, fail: 48}` → 新 `{B: 2, C: 34, fail: 14}`。候选 50，A=0（含 prepass 0），B=2，C=34，D=0，final=36。

`rescue_summary.md` 同时被结构化清理：旧的「未通过任何分支 20 / 无属分类 10 / 长度达标但未被拯救 3」等超额占比（>100%）条目被移除，新表只剩「过短 (<2000bp) 13」并合计 14（100%）；表头改为「合计 (候选表内) 36 / 未拯救 (候选表内) 14」。这符合 `regen_rescue_reports.py` 声明的两处修复（C 分支 Completo 直达 pass 未计数；A 计数含 prepass 混入分母）。

注意：`$D/10_Reports/rescue_report.tsv` 是 report_pipeline 在第 5 步从 `08_Rescue/Plant/` 拷入的**上一代**（09-15 23:49）版本，与第 6 步刚刷新出的 `08_Rescue/Plant/rescue_report.tsv`（01:07:10）不同步——两者差 64 B。属步骤先后顺序产物，非本轮共识问题。

## 9. 备份路径与日志

- `$D/05_Taxonomy/Votus.integrated.bak_v66f_20260916/`（旧成品完整目录）
- `$D/10_Reports.bak_v66f_20260916/`（旧报告目录）
- `$D/08_Rescue/Plant/rescue_report.tsv.bak_v661`、`rescue_summary.md.bak_v661`（08-20 既有，未新建）
- 运行日志：`/tmp/refresh_v67/Fusarium/`（`old.tsv`、`r.log`、`chimera.txt`、`cmp.txt`、`report.log`、`rescue.log`、`manifest.tsv`）；stdout：`/tmp/refresh_v67_Fusarium.nohup`
- 独立复核脚本：本地 `scripts/audit/refresh_v67/audit_fusarium.py`，服务器 `/tmp/audit_fusarium.py`（md5 `bf062740558493276aa57a1a7c218d59`）

## 10. 耗时

服务器 01:06:06 启动 → 01:07:10 完成，共 **64 秒**（R 重跑约 20 s，报告+rescue 约 25 s，其余为检查/备份/对照）。R rc=0，report rc=0，rescue rc=0，无 FATAL。

## 11. 问题与不确定处

1. **改后矛盾率 0.18% 低于预期 1%**：方向正确且更彻底（barbarum 是 0.42%），非异常，但说明不同数据集对门槛语义的敏感度差异较大。
2. **10_Reports 中 4 个文件 mtime 保留旧值**：已定位为 report_pipeline 的 `shutil.copy2` 拷贝上游 09 阶段/08_Rescue 旧代产物，**不是本轮共识变更驱动**，也不影响最终结论，但看 mtime 会误判，需留意。
3. **`$D/10_Reports/rescue_report.tsv` 与 `$D/08_Rescue/Plant/rescue_report.tsv` 不同步**：前者是步骤 5 拷入的上一代（09-15 23:49），后者是步骤 6 新刷（01:07:10）。若下游以 10_Reports 为准会读到偏旧版本，建议后续统一口径。
4. **待决项（本轮未动，按纪律保留）**：
   - `prevalence_full_table.tsv` 与 `09/09b` 派生表本轮未重跑（其内容仍停在 09-15 14:26 一代）；
   - 未运行 `virome_pipeline.py --stage analysis`（`run_analysis` 会给分类表就地补 `Nucleic_acid` 列，属 schema 变化，超出本次范围）。
5. 参照库对 10 个 name 未命中（929 查 / 919 中），与参考实现同源同结果，未影响判据。
6. 所有 ssh 内联命令均以单引号包裹发送，未触发 PowerShell 吞引号问题。
