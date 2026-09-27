# refresh_v67 独立复核报告 — Alternaria

- 数据集 label：`Alternaria`（RNA-Alternaria_alternata_out）
- 数据集目录：`/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Alternaria_alternata_out`
- 执行器：`/tmp/refresh_v67.sh`（md5 `d832c4868a2b8298d7bfca44be0dcb18`，未修改）
- 关键脚本：`virus_classifier_analysis.R`（md5 `f0683ac6af1a9ee1395a100d128365ca`，TAX_GATE_VERSION 6.7）
- 复核人：独立重算（自写 `/tmp/audit_alternaria.py`，md5 `e7f7f490ad4fbb0299465f5ec9ada847`，未复制参考脚本代码，仅对齐判据语义）
- 服务器：`zhangwenda@202.119.189.246`

## 1. 关键数字

| 项目 | 旧（v6.6f / 线上成品） | 新（v6.7 重跑） |
|---|---|---|
| 成品文件 | `05_Taxonomy/Votus.integrated/final_integrated_classification.tsv` | 同 |
| md5 | `78f6677f850b7e7d009d52fb2cb187cd` | `79ff593bc4a06c16f003b89ea0fc9edf` |
| 行数 | 1728 数据行（1729 含表头） | 1728 数据行（1729 含表头） |
| contig_id 集合 | — | 与旧**完全相同**（集合相等，且顺序一致；仅旧有 0 / 仅新有 0） |
| 相邻阶元矛盾行 | 356 行 | 6 行 |
| 矛盾率 | **20.60%**（基线锚点，命中） | **0.35%** |
| 可比格子数 | 9988 | 9713 |

矛盾率锚点：任务给定旧基线 ≈20.60%，实测 20.60%，对上。

### 7 对相邻阶元矛盾数（独立重算，与 `/tmp/refresh_v67/Alternaria/chimera.txt` 逐格一致）

| 阶元对 | 旧 (v66f) | 新 (v6.7) |
|---|---|---|
| Realm-Kingdom | 6 | 1 |
| Kingdom-Phylum | 5 | 0 |
| Phylum-Class | 4 | 2 |
| Class-Order | 43 | 2 |
| Order-Family | 81 | 0 |
| Family-Genus | 167 | 0 |
| Genus-Species | 192 | 2 |
| **矛盾行** | **356** | **6** |
| **矛盾率** | **20.60%** | **0.35%** |

（矛盾计数总和 旧 498 / 新 7，因单行可同时命中多对，故与矛盾行数不等。）

### 各阶元取值变化行数（旧→新，与 `cmp.txt` 一致）

| 阶元 | 变化行数 |
|---|---|
| Realm | 59 |
| Kingdom | 60 |
| Phylum | 67 |
| Class | 71 |
| Order | 92 |
| Family | 149 |
| Genus | 216 |
| Species | 349 |

- Order 有值→空：**64**；空→有值：0
- Family 有值→空：**58**；空→有值：0

### Phylum↔Class 参照相容性（旧→新）

| 类别 | 行数 |
|---|---|
| both_ok | 1722 |
| fixed（旧矛盾→新相容） | 4 |
| broken（旧相容→新矛盾） | 2 |
| both_bad | 0 |

### 判据敏感性

dmp 同名「首次出现」与「最后一次出现」两种取名策略，旧产物矛盾行均为 356，结果不敏感。

## 2. gate 自检原文

`05_Taxonomy/Votus.integrated/taxonomy_gate_stamp.tsv`：

```
gate_version	script_path	script_md5	ref_file	ref_md5	rows	timestamp
6.7	/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R	f0683ac6af1a9ee1395a100d128365ca	/home/zhangwenda/database/taxonomy/genus_family_ref.tsv	865283b4b1a3a2f3c4c3ac702403c554	1728	2026-09-16 01:07:19
```

`05_Taxonomy/Votus.integrated/taxonomy_gate_check.tsv`：

```
check	n	expected	status
rows_total	1728	-	INFO
family_filled	1409	-	INFO
genus_filled	1188	-	INFO
genus_ref_known	1071	-	INFO
dual_ref_consistent	971	-	INFO
dual_ref_conflict	0	0	PASS
single_ref_conflict	0	unfixable_by_design	INFO
species_dual_ref_conflict	0	0	PASS
species_ref_known	1000	-	INFO
species_single_ref_conflict	0	unfixable_by_design	INFO
```

结论：版本 6.7 ✓，script_md5 = `f0683ac6af1a9ee1395a100d128365ca` ✓，`dual_ref_conflict` = **PASS** ✓。

R 日志要点（`/tmp/refresh_v67/Alternaria/r.log`）：
`逐级淘汰: 阶元计票 11995 格次, 淘汰工具投票权 1321 格次 (票首门槛 0.50)`；`最终分类输出: 1728 个序列`；`全流程完美执行完毕`。

## 3. 10_Reports 变化清单

基线：`10_Reports.bak_v66f_20260916`（由脚本在刷新前 `cp -a` 生成，即 09-15 00:59 的状态）。
`diff -rq` 结果：**无仅存于一方的文件**，内容变化 17 个：

| 文件 | 旧/新 行数 | 旧/新 字节 | 性质 |
|---|---|---|---|
| final_integrated_classification.tsv | 1729 / 1729 | 612583 / 603427 | taxonomy 派生（＝新成品，md5 同步为 `79ff593b…`） |
| taxonomy_composition.tsv | 893 / 836 | 25252 / 23470 | taxonomy 派生 |
| taxonomy_fill_stats.tsv | 9 / 9 | 165 / 166 | taxonomy 派生 |
| taxonomy_sunburst.html | 3887 / 3887 | 4872992 / 4873193 | taxonomy 派生（可视化） |
| plant_final_taxonomy.tsv | 33 / 33 | 11529 / 11380 | taxonomy 派生 |
| plant_virus_summary.tsv | 25 / 29 | 5176 / 6303 | taxonomy 派生 |
| all_plant_viruses_genus_summary.tsv | 19 / 20 | 396 / 804 | taxonomy 派生 |
| All_plant.viruses_info.tsv | 31 / 31 | 7488 / 8127 | taxonomy 派生 |
| HQ_plant_viruses_info.tsv | 25 / 25 | 3890 / 4403 | taxonomy 派生 |
| classification_sankey.html | 3887 / 3887 | 4886920 / 4881493 | taxonomy 派生（可视化） |
| classification_sankey_plant.html | 3887 / 3887 | 4866789 / 4866726 | taxonomy 派生（可视化） |
| upset_data.json | 0 / 1 | 172189 / 160608 | taxonomy 派生（工具一致性 upset 数据） |
| rescue_report.tsv | 30 / 34 | 3051 / 3905 | rescue 层（见 §4） |
| directory_tree.txt | 475 / 507 | 12808 / 13836 | **非 taxonomy**：目录树纳入新建的备份目录 |
| pipeline_report.html | 933 / 933 | 23184223 / 23156208 | **非 taxonomy**：整体重渲染 |
| filter_summary.tsv | 40 / 14 | 1333 / 364 | **非 taxonomy**：报告生成器结构差异（旧为 filter/strict/comb 三模式逐样本，新为 comb 单模式） |
| stage_summary.tsv | 48 / 48 | 3859 / 4003 | **非 taxonomy**：内嵌绝对路径由旧 `/data-test/…` 更正为本目录 `/home/…/goji-virome/…` |

未变化的 taxonomy 相关文件：`taxonomy_agreement_stats.tsv`、`taxonomy_consistency_summary.tsv`（内容逐字节相同）。
被重写但内容相同的文件（新 mtime 01:07，`diff` 无差异）：`checkv_summary.tsv`、`checkv_confidence.tsv`、`ident_summary.tsv`、`cobra_summary.tsv`、`host_*.tsv`、`assembly_summary.tsv`、`cluster_*` 等，说明流水线对其重算结果稳定。

## 4. rescue 报告

`08_Rescue/Plant/`（基线为 `*.bak_v661`，2026-08-20 10:58）：

| 文件 | 旧 mtime | 新 mtime | 旧/新 行数 | 旧/新 字节 |
|---|---|---|---|---|
| rescue_report.tsv | 2026-08-20 10:58:49 | **2026-09-16 01:07:48** | 34 / 34 | 3458 / 3791 |
| rescue_summary.md | 2026-08-20 10:58:49 | **2026-09-16 01:07:48** | 25 / 23 | 692 / 637 |

两者相对 `.bak_v661` 内容均**有变化**（`diff -q` 返回差异）。

注意：该目录在本次刷新前曾于 2026-09-15 23:49 被刷新过一次（旧字节 3905 / 637），`*.bak_v661` 是 8/20 的更老基线，故「相对 bak_v661」的差异**混入了 8/20→9/15 的旧改动**，不能全部归因于本轮。用 10_Reports 快照隔离本轮增量：

`10_Reports/rescue_report.tsv`（09-15 快照 30 行 → 本轮 34 行，48 行发生变化）显示，一批条目由 `fail / 未通过任何分支` 变为 `C / blastn_completo / BLASTN qcov≥98% 直达`，`rescue_summary.md` 的「未拯救原因分布」表也由旧的多行（含 380% 等越界占比）改为干净的 100% 口径。**该 rescue 变化来自 rescue 报告生成器（`regen_rescue_reports.py` / `report_pipeline.py`）自身的策略与版本，而非 v6.7 taxonomy gate 直接导致**，属需知会的范围外增量。

`10_Reports/rescue_report.tsv`（本轮 md5 `0ba08be8…`）与 `08_Rescue/Plant/rescue_report.tsv`（本轮 md5 `c12c1e04…`）内容不相同，两者由不同步骤产出。

## 5. 备份路径（本轮全部只读校验确认，未删除任何文件）

- `05_Taxonomy/Votus.integrated.bak_v66f_20260916/`（含旧成品 md5 `78f6677f…`）
- `10_Reports.bak_v66f_20260916/`
- `08_Rescue/Plant/rescue_report.tsv.bak_v661`、`rescue_summary.md.bak_v661`（既有）
- 脚本内旧成品副本：`/tmp/refresh_v67/Alternaria/old.tsv`（md5 `78f6677f…`）

日志：`/tmp/refresh_v67/Alternaria/`（r.log、report.log、rescue.log、chimera.txt、cmp.txt、manifest.tsv、old.tsv）；nohup：`/tmp/refresh_v67_Alternaria.nohup`。

## 6. 耗时

- 脚本 START `2026-09-16 01:06:40` → DONE `2026-09-16 01:07:48`，约 **68 秒**。
- 其中 R `build_consensus` 合计 10.1 秒；10_Reports 重生成至 01:07:47；rescue 至 01:07:48.2。

## 7. 未触碰项（按纪律保留为待决）

- `prevalence_full_table.tsv`：mtime `2026-08-25 13:13`，未变。
- `09_Virome_Analysis/`（08-25）、`09b_Analysis_Verify/`（09-15 00:59）：未变。
-`.analysis.ok`（08-26 13:58）、`.report.ok`（08-25 16:21）：未变。
- 未运行 `virome_pipeline.py --stage analysis`，故分类表未新增 `Nucleic_acid` 列（属本轮范围外的 schema 变化）。

## 8. 问题与不确定处

1. **重复触发（已查明，无副作用）**：`/tmp/refresh_v67_Alternaria.nohup` 同时含 `FATAL 备份已存在` 与 `DONE`。经查是两个进程写同一文件：我于 01:06:40 启动的 run A 正常跑完（DONE 01:07:48）；另有一个 run B 于 01:07:41 启动，在**前置检查阶段**即因备份已存在而 `FATAL` 退出，**未做任何改动**（备份就在 step 0 之前被拦截）。nohup 中的 NUL 空洞来自两进程共用 fd 的写入偏移错位，非数据损坏。run B 来源未确认，怀疑与服务器上的并行编排有关（同期 `onekp`、`ruthenicum` 亦在跑其他 label）。若为编排器统一派发，属预期内重复。
2. **rescue 报告变化非 gate 驱动**（见 §4），建议确认是否符合预期。
3. **备份目录被纳入报告产物**：`Votus.integrated.bak_v66f_20260916/` 建在 `05_Taxonomy` 内，导致 `directory_tree.txt`、`pipeline_report.html` 出现备份目录条目（轻微美观污染，无功能影响）。
4. **10_Reports 生成器版本漂移**：`filter_summary.tsv`（40→14 行，comb 值全 0）与 `stage_summary.tsv`（绝对路径 `/data-test/…`→本目录）显示旧 10_Reports 由更早的生成器/在别处（data-test）产出，本轮刷新统一到当前生成器口径。属报告层结构变化，非 taxonomy 数据变化，建议知会。
5. `upset_data.json` 被重写为单行（160608 B），内容为新一致性数据，格式压缩无影响。
