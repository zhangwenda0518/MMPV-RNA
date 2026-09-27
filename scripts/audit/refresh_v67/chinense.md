# chinense 门槛 v6.7 重跑审计报告

- **数据集 label**：`chinense`
- **目录**：`/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_chinense_out`（下称 `$D`）
- **服务器**：`zhangwenda@202.119.189.246`
- **执行器**：`/tmp/refresh_v67.sh`（md5 `d832c4868a2b8298d7bfca44be0dcb18`，未修改）
- **审计时间**：2026-09-16 01:05–01:12（北京时间）
- **执行耗时**：START `01:05:43` → DONE `01:07:29`，墙钟约 **106 秒（1 分 46 秒）**；R rc=0、report rc=0、rescue rc=0，无 FATAL。
- **审计性质**：全部关键数字均由本人独立重算（自写 `/tmp/verify_chinense.py`，未复用 `chimera_multi.py`），仅将脚本产物用于交叉对照。

---

## 1. 共识成品（05_Taxonomy/Votus.integrated/final_integrated_classification.tsv）

| 项 | 旧（v6.6 线上成品） | 新（v6.7） |
|---|---|---|
| md5 | `da1488531f5e00f031edcebec7da5103` | `977bb21b8f0d745c2c74dae33aa1f60c` |
| 行数（`wc -l`） | 11008（数据行 11007 + 表头） | 11008（数据行 11007 + 表头） |
| 数据行 | 11007 | 11007 |
| mtime | 2026-08-12 02:14:59 | 2026-09-16 01:06:37 |

- **contig_id 集合**：两版完全相同（11007），`only_old=0 / only_new=0`（集合比对，非仅计数）。
- `10_Reports/final_integrated_classification.tsv` 已同步为新 md5（旧值为 `da1488…` 存于备份目录）。

---

## 2. Gate 自检（原文）

**`taxonomy_gate_stamp.tsv` 第 2 行**：

```
6.7	/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R	f0683ac6af1a9ee1395a100d128365ca	/home/zhangwenda/database/taxonomy/genus_family_ref.tsv	865283b4b1a3a2f3c4c3ac702403c554	11007	2026-09-16 01:06:37
```

- 版本 = `6.7` ✅；`script_md5 = f0683ac6af1a9ee1395a100d128365ca`（与要求一致）✅；rows = 11007。

**`taxonomy_gate_check.tsv` 原文**：

```
check	n	expected	status
rows_total	11007	-	INFO
family_filled	9531	-	INFO
genus_filled	7854	-	INFO
genus_ref_known	7394	-	INFO
dual_ref_consistent	6780	-	INFO
dual_ref_conflict	0	0	PASS
single_ref_conflict	3	unfixable_by_design	INFO
species_dual_ref_conflict	0	0	PASS
species_ref_known	6482	-	INFO
species_single_ref_conflict	4	unfixable_by_design	INFO
```

- `dual_ref_conflict` = 0 / expected 0 → **PASS** ✅；`species_dual_ref_conflict` → PASS ✅。

---

## 3. 矛盾率（本人独立重算，判据与参考实现一致）

判据：对每对相邻阶元 (H 高, L 低)，用 `rankedlineage.dmp`（列序 `tax_id|tax_name|species|genus|family|order|class|phylum|kingdom|superkingdom`）查 L 的定型父级，低阶元有参照父级且行内高阶元有值时才比，不等即矛盾。参照名命中 `2330 / 需查 2344`。

| 口径 | 行数 | 矛盾行 | 矛盾率 |
|---|---|---|---|
| 旧（v6.6f 备份） | 11007 | **1258** | **11.43%** |
| 新（v6.7 产物） | 11007 | **57** | **0.52%** |

**7 对相邻阶元分布**（旧 → 新）：

| 阶元对 | 旧 | 新 |
|---|---|---|
| Realm–Kingdom | 16 | 14 |
| Kingdom–Phylum | 22 | 3 |
| Phylum–Class | 28 | 8 |
| Class–Order | 220 | 3 |
| Order–Family | 180 | 1 |
| Family–Genus | 691 | 3 |
| Genus–Species | 292 | 27 |

- 基线锚点校验：旧矛盾率 **11.43%** 与预期（线上 v6.6 值 ≈11.43%）吻合 ✅；新值 **0.52%**，落在「1% 上下」区间 ✅。
- 与 `/tmp/refresh_v67/chinense/chimera.txt` **逐格完全一致**（含 7 对分布、矛盾行、矛盾率），无出入。

**Phylum↔Class 参照相容性**（`cmp.txt`，本人复核）：`both_ok 10971 / fixed 28 / broken 8 / both_bad 0` —— 修好 28、弄坏 8，净改善。

---

## 4. 各阶元取值变化行数 / 有值变空

**各阶元取值变化行数**（旧 vs 新，共同 11007 行）：

| Realm | Kingdom | Phylum | Class | Order | Family | Genus | Species |
|---|---|---|---|---|---|---|---|
| 100 | 110 | 120 | 120 | 308 | 394 | 1064 | 1811 |

- **Order 有值 → 空：176 行；Family 有值 → 空：179 行**
- Order 空 → 有值：0；Family 空 → 有值：0（单向信息收敛，无新增填充）
- 与 `/tmp/refresh_v67/chinense/cmp.txt` 一致。

---

## 5. 10_Reports 报告层

- 备份：`$D/10_Reports.bak_v66f_20260916`（`cp -a`，创建于 01:06:50，目录 mtime 继承源码 09-15 00:59，未嵌套）。
- 逐文件 `diff -q` 结果：**17 个文件内容变化**；另有 12 个文件被今日重写但字节完全相同（见下）。
- **mtime 并非可靠的"今日刷新"信号**：4 个变化文件保留了旧 mtime（见注）。

**变化文件清单（17）**，按派生性质分三类：

1. **今日从 v6.7 taxonomy 重算/重绘（mtime 为今日）**：
   - `final_integrated_classification.tsv`（成品同步）
   - `plant_final_taxonomy.tsv`、`plant_virus_summary.tsv`
   - `taxonomy_composition.tsv`、`taxonomy_fill_stats.tsv`
   - `taxonomy_sunburst.html`、`classification_sankey.html`、`classification_sankey_plant.html`、`pipeline_report.html`
   - `upset_data.json`、`stage_summary.tsv`、`filter_summary.tsv`、`directory_tree.txt`

2. **由 report_pipeline 从上游"保留时间戳 copy"同步、尚未反映今日 v6.7 taxonomy（mtime 仍为旧值）**：
   - `all_plant_viruses_genus_summary.tsv`、`All_plant.viruses_info.tsv`、`HQ_plant_viruses_info.tsv`
     → 三者与 `$D/09_Virome_Analysis/`（`all_plant_analysis/`、`HQ_analysis/`）下同名文件**字节与 mtime 完全一致**（09-15 14:26），即来源是 09 阶段旧产物。
   - `rescue_report.tsv` → 与运行前 `08_Rescue/Plant/rescue_report.tsv` 逐字节一致（09-15 23:49），是从 regen 之前的旧版 copy。

3. **被今日重写但内容字节相同（12，未进入 diff 清单）**：
   `assembly_summary.tsv`、`checkv_confidence.tsv`、`checkv_summary.tsv`、`cluster_pipeline_reduction.tsv`、`cluster_size_distribution.tsv`、`cobra_summary.tsv`、`host_decision_method.tsv`、`host_distribution.tsv`、`ident_summary.tsv`、`taxonomy_agreement_stats.tsv`、`taxonomy_consistency_summary.tsv`、`tool_filter_summary.tsv`

---

## 6. rescue 报告

| 文件 | 新 | `.bak_v661` 备份 |
|---|---|---|
| `08_Rescue/Plant/rescue_report.tsv` | 55950 B，md5 `74b78c4afd75977987d2fc324f66985b`，mtime 09-16 01:07:29 | 52923 B，mtime 08-19 17:27 |
| `08_Rescue/Plant/rescue_summary.md` | 821 B，mtime 09-16 01:07:29 | 756 B，mtime 08-19 17:27 |

- `rescue_report.tsv`：新旧均 543 行；`diff` 输出 **608 行**（变化显著）。
- `rescue_summary.md`：分支 A 由 `9` 变为 `0`（prepass 免拯救移出候选表分母），候选表内合计 `133 → 124`，未拯救 `409 → 418`；旧版原因分布存在 `>100%` 的计数错误（过短 440 / 占比 107.6%），新版已修正（过短 398 / 占比 95.2%）。
- **重要口径差异**：`.bak_v661` 是 **08-19** 版本，早于上一代（09-15 23:49）rescue 报告，故 608 行 diff 混合了两轮变更，无法单独隔离今日变化。用 `10_Reports/rescue_report.tsv`（保留的运行前 09-15 23:49 副本，55597 B，md5 `37c53857191b25ee1316462a396fec2c`）近似隔离今日变化，得 **296 行 diff**。

---

## 7. 备份与日志路径

- 共识目录备份：`$D/05_Taxonomy/Votus.integrated.bak_v66f_20260916`（`cp -a` 运行前整目录）
- 报告目录备份：`$D/10_Reports.bak_v66f_20260916`（`cp -a` 运行前整目录）
- rescue 备份：`$D/08_Rescue/Plant/rescue_report.tsv.bak_v661`、`rescue_summary.md.bak_v661`（08-19，**本轮未更新**——regen 仅在备份缺失时才写）
- 日志：`/tmp/refresh_v67/chinense/`（`old.tsv`、`r.log`、`chimera.txt`、`cmp.txt`、`manifest.tsv`、`report.log`、`rescue.log`）
- nohup：`/tmp/refresh_v67_chinense.nohup`
- 审计脚本：`/tmp/verify_chinense.py`（md5 `740aa99e69e12d0c80fc5509c1dadc06`）

---

## 8. 本轮 R 在 Votus.integrated 下的文件增删

- 未见删除任何文件。
- 相对运行前（备份清单 + 侦察清单双重确认），R 新增 4 个文件：
  `taxonomy_gate_stamp.tsv`、`taxonomy_gate_check.tsv`、`tool_weights.tsv`、`species_containment_blanked.tsv`。

---

## 9. 问题与不确定处（均未自行处置，按纪律只记录）

1. **`10_Reports/rescue_report.tsv` 落后一代**：执行器 step5（report_pipeline）先于 step6（regen rescue），导致 `10_Reports` 内这份是从"regen 之前"的 08_Rescue 版本 copy 而来，与现在的 `08_Rescue/Plant/rescue_report.tsv` 相差 296 行。是否需要在 rescue 刷新后重跑一次报告层，待定。
2. **3 个 plant 表未反映 v6.7 taxonomy**：`all_plant_viruses_genus_summary.tsv`/`All_plant.viruses_info.tsv`/`HQ_plant_viruses_info.tsv` 来自未重跑的 09 阶段（09-15 14:26）。与 `prevalence_full_table.tsv`（08-25）、`09/09b` 派生表一起列为**待决项**（本轮按纪律未跑 `virome_pipeline.py --stage analysis`，未动 `prevalence_full_table.tsv` 与 09 目录；`find $D -newermt 2026-09-16` 对二者零命中，确认未被触碰）。
3. **mtime 语义**：report_pipeline 对部分文件采用"保留源时间戳"的 copy，mtime 不能作为"是否已按 v6.7 刷新"的充分判据，判断时须以内容 md5/diff 为准。
4. **`filter_summary.tsv` 变化非 taxonomy 明显派生**（4193 B → 1122 B），疑为 report_pipeline 自身行为差异，未深究，标记为待确认。
5. **`directory_tree.txt` 变化**含新建备份目录（`.bak_v66f_20260916`）条目，属结构性变化而非纯 taxonomy 变化。
6. **`.bak_v661` 基线偏老**：见第 6 节，608 行 diff 含两轮变更，已在报告中用 296 行近似隔离今日变化。
7. **与 barbarum 的差异**：barbarum 行数为 20892，chinense 为 11007（不同数据集）；两数据集 md5 各自独立，仅"行数不变、集合不变、矛盾率大幅下降"的模式一致。
8. **环境小插曲**：首次轮询因 PowerShell 吞掉内联双引号导致远程 bash 语法错误，改用无引号写法后正常；对任务本身无影响。

---

## 10. 结论

- gate 版本 `6.7`、script_md5 `f0683ac6…` 校验通过；`dual_ref_conflict` PASS。
- 新旧行数一致、contig_id 集合完全一致；矛盾率 **11.43% → 0.52%**（旧值命中锚点），7 对阶元分布与参考实现逐格吻合。
- 报告层已刷新（17 文件内容变化）；rescue 报告已刷新（今日变化 296 行，对照 `.bak_v661` 608 行）。
- 遗留 3 个 plant 表与 `10_Reports/rescue_report.tsv` 的时间线/来源问题、`prevalence_full_table.tsv` 与 09/09b 派生表为待决项，见第 9 节。
