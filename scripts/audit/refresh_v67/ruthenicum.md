# ruthenicum —— 门槛新语义（v6.7, share >= 0.5 半数票也淘汰）重跑独立复核

复核人：Hanako（独立重算，不转述脚本输出）
复核时间：2026-09-16 01:07–01:15 (GMT+8)
执行器：`/tmp/refresh_v67.sh`（md5 `d832c4868a2b8298d7bfca44be0dcb18`，未改动）
线上 R 脚本：`/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R`
　　复核时 md5 = `f0683ac6af1a9ee1395a100d128365ca`，与 TAX_GATE_VERSION=6.7 要求一致。

## 0. 结论摘要

| 项目 | 结果 |
|---|---|
| 执行 | DONE，无 FATAL（START 01:07:24 → DONE 01:12:36，耗时 5m12s） |
| gate 版本 | `6.7`，script_md5 = `f0683ac6af1a9ee1395a100d128365ca` ✅ |
| gate 自检 | `dual_ref_conflict = 0 → PASS`，`species_dual_ref_conflict = 0 → PASS` ✅ |
| 行数 | 旧 15181 = 新 15181，`contig_id` 集合与顺序完全一致 ✅ |
| 矛盾率（自算） | 旧 **12.11%**（1838 行）→ 新 **0.34%**（51 行） |
| 与执行器对照 | 完全一致（chimera.txt 逐格元命中） ✅ |
| 基线锚点 | 旧值 12.11% 与任务给定 v6.6 锚点吻合 ✅ |
| 新值 vs 预期 | 低于「1% 上下」的粗预期，与 barbarum 的 0.42% 同量级（见 §11） |

## 1. 数据集与目录

- label：`ruthenicum`
- 目录 `D` = `/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_ruthenicum_out`
- 前置侦察确认存在：
  - `D/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv`（44101 行，md5 `557effbb7074e37701f87b0ffa8d2730`，未改动）
  - `D/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv`
- 改前时间戳：`10_Reports/` mtime `Sep 15 00:59`；`.analysis.ok` = `2026-08-26 13:59:46`；`.report.ok` = `2026-08-25 16:58:40`；`.taxonomy.ok` = `2026-08-12 01:47:46`
- 注：本次运行**未改动**任何 `.ok` 戳（重跑只刷新产物与报告，不触碰 stage 完成戳）。

## 2. 终态指纹

| | 旧（v6.6 线上成品，改前） | 新（v6.7） |
|---|---|---|
| 文件 | `05_Taxonomy/Votus.integrated/final_integrated_classification.tsv` | 同左（就地覆盖） |
| md5 | `8d143dc3a684bd52845dc1d6d21954a4` | `2416cc1374a88c3643288b351e9109a5` |
| 数据行 | 15181 | 15181 |
| 字节 | 5484654 | 5434317 |
| mtime | Aug 12 01:47 | Sep 16 01:08 |

- 与背景中 barbarum 的 md5（旧 `12e8f171…` → 新 `a3ead5b5…`）不同属正常，各数据集内容不同。
- `10_Reports/final_integrated_classification.tsv` md5 与新产物**完全相同**（`2416cc13…`，5434317 B），报告层副本已同步。

## 3. gate 自检原文

`05_Taxonomy/Votus.integrated/taxonomy_gate_stamp.tsv`：

```
gate_version	script_path	script_md5	ref_file	ref_md5	rows	timestamp
6.7	/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R	f0683ac6af1a9ee1395a100d128365ca	/home/zhangwenda/database/taxonomy/genus_family_ref.tsv	865283b4b1a3a2f3c4c3ac702403c554	15181	2026-09-16 01:08:28
```

`05_Taxonomy/Votus.integrated/taxonomy_gate_check.tsv`：

```
check	n	expected	status
rows_total	15181	-	INFO
family_filled	13008	-	INFO
genus_filled	11114	-	INFO
genus_ref_known	10381	-	INFO
dual_ref_consistent	9768	-	INFO
dual_ref_conflict	0	0	PASS
single_ref_conflict	2	unfixable_by_design	INFO
species_dual_ref_conflict	0	0	PASS
species_ref_known	8744	-	INFO
species_single_ref_conflict	4	unfixable_by_design	INFO
```

R 日志门槛行（确认新语义生效）：
`[01:08:14] INFO: 逐级淘汰: 阶元计票 106365 格次, 淘汰工具投票权 9048 格次 (票首门槛 0.50); 平票 3303 格次: 下一级一致性定 430, 工具顺序兜底 2873`

## 4. 行数与集合核对（自写脚本，非计数）

用集合而非计数比对：

- old data rows = 15181（唯一 id 15181）；new data rows = 15181（唯一 id 15181）
- 仅旧有 id：0；仅新有 id：0
- **id 集合相同：True**；**行顺序亦完全相同：True**

即本次变更为纯「同集合内的取值改写」，无增删、无重排。

## 5. 矛盾率（独立重算）

判据与参考实现一致：对每对相邻阶元 (H 高, L 低)，用
`/home/zhangwenda/database/taxonomy/rankedlineage.dmp`（列序 `tax_id|tax_name|species|genus|family|order|class|phylum|kingdom|superkingdom`）查 L 的定型父级（父级位于 H 对应的 dmp 列），仅当 L 有可查父级且行内 H 有值时才比较，不等即该行矛盾。参照 name 命中 2907/2928。

| 阶元对 (H-L) | 旧(v6.6) 矛盾数 | 新(v6.7) 矛盾数 |
|---|---|---|
| Realm–Kingdom | 69 | 29 |
| Kingdom–Phylum | 42 | 1 |
| Phylum–Class | 32 | 8 |
| Class–Order | 430 | 3 |
| Order–Family | 177 | 2 |
| Family–Genus | 878 | 2 |
| Genus–Species | 644 | 9 |
| **矛盾行合计** | **1838** | **51** |
| **矛盾率** | **12.1072% (≈12.11%)** | **0.3359% (≈0.34%)** |

- 与执行器 `/tmp/refresh_v67/ruthenicum/chimera.txt` **逐格元完全一致**，无出入。
- 旧值 12.11% = 任务给的 v6.6 基线锚点 ✅。

## 6. 各阶元取值变化 / 有值变空

自写脚本按 `contig_id` 对齐 15181 条共同行：

| 阶元 | 取值变化行数 | 有值→空 | 空→有值 |
|---|---|---|---|
| Realm | 100 | 0 | 0 |
| Kingdom | 69 | 0 | 0 |
| Phylum | 68 | 0 | 0 |
| Class | 124 | 2 | 0 |
| Order | 468 | 294 | 0 |
| Family | 618 | 275 | 0 |
| Genus | 1467 | 850 | 0 |
| Species | 2446 | 722 | 0 |

- **Order/Family 有值变空：Order 294 / Family 275**（空→有值均为 0，无回填）。
- 与执行器 `cmp.txt` 逐项一致。
- 填充率自洽核验（报告层 `taxonomy_fill_stats.tsv`，与 `gate_check` 的 family_filled/genus_filled 互证）：

| 阶元 | 旧 Filled (Rate) | 新 Filled (Rate) | Δ |
|---|---|---|---|
| Realm | 14044 (0.9251) | 14044 (0.9251) | 0 |
| Kingdom | 14107 (0.9293) | 14107 (0.9293) | 0 |
| Phylum | 14029 (0.9241) | 14029 (0.9241) | 0 |
| Class | 14331 (0.9440) | 14329 (0.9439) | −2 |
| Order | 13150 (0.8662) | 12856 (0.8468) | −294 |
| Family | 13283 (0.8750) | 13008 (0.8569) | −275 |
| Genus | 11964 (0.7881) | 11114 (0.7321) | −850 |
| Species | 13375 (0.8810) | 12653 (0.8335) | −722 |

## 7. 10_Reports 变化清单

`10_Reports/` 目录 mtime 已刷新（多数文件落 Sep 16 01:07–01:12）；两侧文件数均 56，文件集合一致（无新增/删除）。
与备份 `D/10_Reports.bak_v66f_20260916` 逐文件 `diff -rq`，**共 18 个文件内容变化**：

**A. taxonomy 派生、本次真刷新（mtime 01:07–01:12）——13 个**
1. `final_integrated_classification.tsv`（= 05 产物副本，md5 `2416cc13…`）
2. `taxonomy_fill_stats.tsv`（01:08）
3. `taxonomy_composition.tsv`（01:12）
4. `taxonomy_sunburst.html`（01:12）
5. `upset_data.json`（01:07，= `05_Taxonomy/Votus.integrated/upset_data.json`）
6. `pipeline_report.html`（01:12）
7. `classification_sankey.html`（01:12）
8. `classification_sankey_plant.html`（01:12）
9. `stage_summary.tsv`（01:12）
10. `directory_tree.txt`（01:12）
11. `plant_final_taxonomy.tsv`（01:12）
12. `plant_virus_summary.tsv`（01:12）
13. `filter_summary.tsv`（01:11）

**B. copy2 从 09/09b 舞台目录同步而来（mtime 是**被保留的旧时间**，内容早于本次重跑）——4 个**
14. `All_plant.viruses_info.tsv`　mtime `Sep 15 14:26`；源 `09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv`
15. `HQ_plant_viruses_info.tsv`　mtime `Sep 15 14:26`；源 `09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv`
16. `all_plant_viruses_genus_summary.tsv`　mtime `Sep 15 14:26`；源 `09_Virome_Analysis/all_plant_analysis/…`
17. `keep_summary.tsv`　mtime `Sep 15 00:59`；源 `09b_Analysis_Verify/keep_summary.tsv`

（依据：`report_pipeline.py:745–762` 用 `shutil.copy2` 复制这些舞台文件，故 mtime 被保留；它们的内容**不是**本次 v6.7 结果，属 §12 待决项。）

**C. rescue_report.tsv（见 §8 的时序问题）**
18. `rescue_report.tsv`　mtime `Sep 15 23:49`（被保留），内容 = 重跑前的 08 报告

**mtime 已刷新但内容未变**（未列入 diff）：`assembly_summary.tsv`、`ident_summary.tsv`、`cobra_summary.tsv`、`checkv_summary.tsv`、`checkv_confidence.tsv`、`host_distribution.tsv`、`host_decision_method.tsv`、`tool_filter_summary.tsv`、`cluster_pipeline_reduction.tsv`、`cluster_size_distribution.tsv`、
以及 `taxonomy_agreement_stats.tsv`、`taxonomy_consistency_summary.tsv`（二者 `diff` 判定 IDENTICAL）。
`Viroid.*`、`plant_virus_cluster_summary.tsv` 等未变（未列入 diff）。

## 8. rescue 报告

| 文件 | 改前 (Sep 15 23:49 / bak) | 改后 | mtime |
|---|---|---|---|
| `08_Rescue/Plant/rescue_report.tsv` | 128930 B（Sep 15 23:49） | **129139 B**，md5 `ba62ad556287bdf93d8454ae2fcf906e` | **Sep 16 01:12** |
| `08_Rescue/Plant/rescue_summary.md` | 824 B | 824 B（内容更新） | **Sep 16 01:12** |

对照 `.bak_v661`（Aug 19 17:29，即更早的原始版）：
- `rescue_report.tsv`：分支 C 的 `blastn_completo` 直达行由 `fail` 改为 `C + BLASTN qcov≥98% 直达`；同时多行 Genus 被新分类改写（如 `Pseudovirus → Caulimovirus`、`Potexvirus → (空)`、`Ampelovirus/…` 走直达）。
- `rescue_summary.md`：分支 A 由「28 / 2.2%」改为「0（另含 28 条免拯救 prepass，候选表外）」；合计口径改为「候选表内 295 / 23.2%」；未拯救 976；未拯救原因重算为「过短 941 / 属组失败 35」；最终无冗余 vOTU = 272。
- 本次新报告 branch 分布：`{'fail': 976, 'C': 274, 'B': 21}`；候选 1271，A=28(含 prepass 28)、B=21、C=274、D=0，final=272。

⚠️ **时序问题（重要）**：`report_pipeline.py`（步骤 5）在 `regen_rescue_reports.py`（步骤 6）**之前**执行，并复制当时的 08 报告到 `10_Reports/`。因此
- `10_Reports/rescue_report.tsv` = 重跑前（Sep 15 23:49）那版，md5 `cf438eaa275d491ab0f1855dedb5a6ae`；
- `08_Rescue/Plant/rescue_report.tsv` = 刚重生成版，md5 `ba62ad556287bdf93d8454ae2fcf906e`；
- 二者**当前不一致（10_Reports 版落后一步）**。

## 9. 备份路径

- 产物备份：`D/05_Taxonomy/Votus.integrated.bak_v66f_20260916/`（cp -a，mtime 保留 Aug 24 02:18）
- 报告备份：`D/10_Reports.bak_v66f_20260916/`（cp -a，56 项）
- 旧成品快照：`/tmp/refresh_v67/ruthenicum/old.tsv`（md5 `8d143dc3a684bd52845dc1d6d21954a4`，mtime 保留 Aug 12 01:47）
- rescue 旧版备份：`D/08_Rescue/Plant/rescue_report.tsv.bak_v661` 与 `rescue_summary.md.bak_v661`（Aug 19 17:29，**预先存在，未被覆盖**）
- 日志：`/tmp/refresh_v67/ruthenicum/`（`r.log`、`chimera.txt`、`cmp.txt`、`report.log`、`rescue.log`、`manifest.tsv`、`old.tsv`）；nohup `/tmp/refresh_v67_ruthenicum.nohup`

## 10. 耗时

- START `2026-09-16 01:07:24` → DONE `01:12:36` = **5 分 12 秒**（R 重跑主体约 01:07:4x–01:08:30，report_pipeline “收集 stage 数据”约占 207 s）。
- 期间服务器上另有并发的同脚本任务（`chinense`/`onekp`/`Alternaria`/`Aphis`/`Fusarium`，由其他会话发起），本数据集进程独立、无报错。

## 11. 遇到的问题与不确定处

1. **新矛盾率 0.34% 低于「1% 上下」的粗预期。** 方向正确、与 barbarum 的 0.42% 同量级；残余 51 行中 Realm–Kingdom 占 29、Genus–Species 占 9，为主要来源。是否满足验收需按项目方的实际阈值判断。
2. **`10_Reports/rescue_report.tsv` 落后 `08_Rescue/Plant/rescue_report.tsv` 一步**（§8 时序问题）。如需两处一致，报告层需在 rescue 重生成之后再同步一次。
3. **4 个 copy2 同步文件（§7-B）内容早于本次重跑**，其 mtime 被保留（Sep 15 14:26 / Sep 15 00:59），不能当作「本次已刷新」的证据。
4. `taxonomy_agreement_stats.tsv`、`taxonomy_consistency_summary.tsv` mtime 更新但内容未变，属正常（仅依赖工具一致性，不依赖门槛结果）。
5. `.taxonomy.ok / .analysis.ok / .report.ok` 时间戳未变，脚本设计如此（不触碰 stage 戳），不代表报告未刷新。
6. gate_check 中 `single_ref_conflict = 2`、`species_single_ref_conflict = 4` 为 `unfixable_by_design`（单侧参照、设计内不动），非失败。

## 12. 待决项（本轮按纪律未动，仅记录）

- `prevalence_full_table.tsv`（Aug 25 13:12）与 `09_Virome_Analysis/`、`09b_Analysis_Verify/` 派生表：本轮**未重跑**。因此经 §7-B 同步进 `10_Reports` 的 `All_plant.viruses_info.tsv`、`HQ_plant_viruses_info.tsv`、`all_plant_viruses_genus_summary.tsv`、`keep_summary.tsv` 仍是旧门槛口径内容（其挂靠的 genus/family 判定可能滞后于 v6.7）。
- 本轮**未**运行 `virome_pipeline.py --stage analysis`（避免 `run_analysis` 给分类表补 `Nucleic_acid` 列的 schema 变化）。

---

### 复核方法说明（可复现）
- 自写复核脚本：`/tmp/audit_ruthenicum.py`（本地源 `scripts/audit/refresh_v67/audit_ruthenicum.py`），独立实现同一判据，未复用执行器的 `chimera_multi.py`。
- 调用：`python3 /tmp/audit_ruthenicum.py /tmp/refresh_v67/ruthenicum/old.tsv D/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv`
- 全程只读比对 + 仅由执行器写入本数据集目录；未删除任何文件；未改动 R/py 管线脚本；未触碰其他数据集。
