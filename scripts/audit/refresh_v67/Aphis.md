# 门槛变更（TAX_GATE 6.7，「半数票也淘汰」）重跑复核 —— Aphis

- 数据集 / label：`Aphis`
- 目录：`/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Aphis_gossypii_out`
- 执行器：`/tmp/refresh_v67.sh`（md5 `d832c4868a2b8298d7bfca44be0dcb18`，未修改）
- 共识脚本：`/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R`（md5 `f0683ac6af1a9ee1395a100d128365ca`，运行前核对一致）
- 运行窗口：2026-09-16 01:07:30 → 01:08:37（耗时 **67 秒**），脚本自报 `=== Aphis DONE ===`
- 本轮未跑 `virome_pipeline.py --stage analysis`；未改任何 R/py 管线脚本；未删文件；只动本数据集目录

## 1. 成品 md5 / 行数

| 项 | 旧（v6.6 线上产物） | 新（gate 6.7） |
|---|---|---|
| md5 | `64491217d3292b794c313cc18525ab2d` | `b72315fe25cddcf76a2ec89778d653e3` |
| 数据行数 | 2448 | 2448 |
| 文件行数（含表头） | 2449 | 2449 |

旧成品 md5 与题面给出的 barbarum 基线值（`12e8f171...`）不同，属正常：那是 barbarum 数据集的值，本数据集有自己的旧 md5。

## 2. gate 自检（原文）

`05_Taxonomy/Votus.integrated/taxonomy_gate_stamp.tsv` 第 2 行：

```
6.7	virus_classifier_analysis.R	f0683ac6af1a9ee1395a100d128365ca	genus_family_ref.tsv	865283b4b1a3a2f3c4c3ac702403c554	2448	2026-09-16 01:08:12
```

`taxonomy_gate_check.tsv`（原文）：

```
check	n	expected	status
rows_total	2448	-	INFO
family_filled	2105	-	INFO
genus_filled	1714	-	INFO
genus_ref_known	1603	-	INFO
dual_ref_consistent	1512	-	INFO
dual_ref_conflict	0	0	PASS
single_ref_conflict	0	unfixable_by_design	INFO
species_dual_ref_conflict	0	0	PASS
species_ref_known	1401	-	INFO
species_single_ref_conflict	1	unfixable_by_design	INFO
```

- gate_version = **6.7** ✓，script_md5 = **f0683ac6af1a9ee1395a100d128365ca** ✓
- `dual_ref_conflict` = **PASS**（n=0）✓；`species_dual_ref_conflict` = **PASS**（n=0）✓

## 3. 行数与集合（自写 setchk.py 复核，非只比计数）

- A(旧=old.tsv) rows=2448 unique=2448 dup=0；B(新) rows=2448 unique=2448 dup=0
- A_only=0，B_only=0 → **contig_id 集合完全相同**（`SET_EQUAL=True`），且顺序也完全相同（`ORDER_IDENTICAL=True`）

## 4. 矛盾率（自写 chi.py 重算，判据对齐 /tmp/chimera_multi.py）

参照库 `rankedlineage.dmp`；对每对相邻阶元 (H 高, L 低) 用 L 的定型父级比对，低阶元有参照父级且行内高阶元有值时才比。参照命中 1180 / 1190 个 name（与参考实现 REF_HITS 一致）。

| 口径 | Realm-Kingdom | Kingdom-Phylum | Phylum-Class | Class-Order | Order-Family | Family-Genus | Genus-Species | 矛盾行 | 矛盾率 |
|---|---|---|---|---|---|---|---|---|---|
| 旧（v6.6） | 68 | 4 | 5 | 90 | 185 | 193 | 111 | **418** | **17.0752%** |
| 新（6.7） | 3 | 1 | 0 | 0 | 0 | 0 | 3 | **7** | **0.2859%** |

- 基线锚点核对：旧值 17.0752% ≈ 题面 17.08% ✓；改动前实测锚点吻合
- 改后 0.29%，低于「1% 上下」的预期区间上沿（barbarum 为 0.42%，本数据集更低）
- 与 `/tmp/refresh_v67/Aphis/chimera.txt` 逐格一致（realm 3 / kingdom 1 / … / species 3 / 7 / 0.29%），无出入
- 残留 7 行中：Realm-Kingdom 3 + Kingdom-Phylum 1 + Genus-Species 3；gate_check 的 `species_single_ref_conflict=1` 属设计内不可修（INFO）

## 5. 各阶元取值变化行数 / 有值变空（自写 chg.py 复核）

共同 contig 2448（old 2448 / new 2448）。各阶元取值发生变化的行数：

| Realm | Kingdom | Phylum | Class | Order | Family | Genus | Species |
|---|---|---|---|---|---|---|---|
| 20 | 51 | 55 | 80 | 158 | 309 | 393 | 460 |

- Order 有值→空：**60**；Family 有值→空：**44**；Order/Family 空→有值：**0 / 0**（无信息回填）
- 与脚本 `/tmp/refresh_v67/Aphis/cmp.txt` 完全一致
- 交叉印证 `taxonomy_fill_stats`：Order 有值 2143→2083（−60），Family 2149→2105（−44），与上表丢失数吻合
- Phylum↔Class 参照相容性：both_ok 2443 / fixed 5 / broken 0 / both_bad 0（修复 5，无新引入矛盾）

## 6. 10_Reports 变化清单

刷新前快照：`10_Reports.bak_v66f_20260916`（本轮 report_pipeline 执行前一刻 `cp -a`，因此其 diff 为纯刷新效应）。全部报告文件 mtime 已刷新到 2026-09-16 01:07–01:08。

`diff -rq` 共 17 个文件内容改变：

**taxonomy / 门槛直接驱动**
1. `final_integrated_classification.tsv`（新成品副本）
2. `taxonomy_fill_stats.tsv`（Order 2143→2083、Family 2149→2105、Genus 1920→1714、Species 2115→1973）
3. `taxonomy_composition.tsv`（各阶元计数重分布，如 Realm Varidnaviria 1365→1360、Duplodnaviria 318→329）
4. `upset_data.json`（R 重跑时重写 `Votus.integrated/upset_data.json`，260822→244188 字节）
5. `taxonomy_sunburst.html`、`classification_sankey.html`、`classification_sankey_plant.html`
6. `plant_final_taxonomy.tsv`（植物病毒行重分类，如 NODE_622 metabuli→ACVirus、物种改判）
7. `pipeline_report.html`（聚合页，随上述内容重生成）

**上游产物陈旧回填 / 环境驱动（非门槛驱动）**
8. `all_plant_viruses_genus_summary.tsv`（外加 Genome_Type/Structure/Source 三列，属 09 分析新版 schema）
9. `All_plant.viruses_info.tsv`、`HQ_plant_viruses_info.tsv`（来自 09_Virome_Analysis，9-15 的新版）
10. `plant_virus_summary.tsv`（26→31 行）
11. `rescue_report.tsv`（report_pipeline 从 08_Rescue 拷贝，见 §7 顺序说明）
12. `filter_summary.tsv`（**变为全 0**，见「问题」）
13. `directory_tree.txt`、`stage_summary.tsv`（根路径由旧 `/home/zhangwenda/virus/data-2026/data-test/...` 切换为当前目录）

**被重写但内容不变**（未出现在 diff 中）：`taxonomy_agreement_stats.tsv`、`taxonomy_consistency_summary.tsv`、`assembly_summary.tsv`、`ident_summary.tsv`、`tool_filter_summary.tsv`、`cobra_summary.tsv`、`cluster_pipeline_reduction.tsv`、`cluster_size_distribution.tsv`、`host_distribution.tsv`、`checkv_confidence.tsv`、`checkv_summary.tsv`、`host_decision_method.tsv`。

**未刷新**：`taxonomy_images.json`（保持 8-20）、`chart.min.js`（6-15，静态资源）。

## 7. rescue 报告变化

文件：`08_Rescue/Plant/rescue_report.tsv`、`rescue_summary.md`（mtime 均刷新到 2026-09-16 01:08:37）。

对照 `.bak_v661`（**注意：该备份是 2026-08-19 的版本，不是本次运行前的直接前驱**，见「问题」）：
- `rescue_report.tsv`：多条由 `fail / 未通过任何分支` 改为 `C / blastn_completo / BLASTN qcov≥98% 直达`（branch C 直达修复），len 列由 `-` 补为实际长度。
- `rescue_summary.md`：A 分支 2→0（另计 prepass）、合计 30→28、未拯救 20→22，「原因」分母改为「占未拯救」并行（过短 21/95.5%、属组失败 1/4.5%）。

**隔离本次步骤 6 的净效应**（用 `10_Reports/rescue_report.tsv`＝步骤 5 拷贝的「步骤 6 之前」快照，对比当前 08_Rescue 版本）：
- 仅 **6 行**差异，全部是 branch C 行（SRR36562186_NODE_1698 / SRR36562183_NODE_2490 / SRR36562186_NODE_5525 / contig_1802 / SRR36562186_NODE_1598 / SRR36562184_NODE_899）的 `genus` / `genus_avg_len` / `pct_of_genus` 由有值改为空。
- 即：fail→C 的大批改动发生在 **9-15 之前**，不是本轮引入；本轮步骤 6 只动了这 6 行 + summary。

## 8. 备份路径（保留，未删除）

- `05_Taxonomy/Votus.integrated.bak_v66f_20260916/`（含旧成品 876495B / 8-12 02:42）
- `10_Reports.bak_v66f_20260916/`
- `08_Rescue/Plant/rescue_report.tsv.bak_v661`、`rescue_summary.md.bak_v661`（8-19 既有备份，脚本未覆盖）

## 9. 遇到的问题与不确定处

1. **filter_summary.tsv 变全 0**：`02b_Filter/<sample>/*.fasta` 是指向 `/home/zhangwenda/goji-virome/01_data/data-test/RNA-Aphis_gossypii_out/...` 的**符号链接，目标已不存在（断链）**，`_count_fasta` 读到 0。属报告刷新的环境副作用，与门槛变更无关；刷新前那份非零值是 8-26 生成的旧文件。不影响 taxonomy 产物。
2. **旧 10_Reports 是旧根目录产物**：刷新前 `stage_summary.tsv` 里的 Details 指向 `/home/zhangwenda/virus/data-2026/data-test/RNA-Aphis_gossypii_out`（该路径已不存在）。所以本次 17 个变化文件里，有一部分是「报告对齐到当前数据集树」的陈旧回填（all_plant/HQ/plant_virus_summary/stage_summary/directory_tree），而非门槛驱动。逐文件归因已见 §6。
3. **10_Reports/rescue_report.tsv 与 08_Rescue 版本不一致**：`report_pipeline.py`（步骤 5，第 691 行 `shutil.copy2`）先把 08_Rescue 的**旧版**拷进 10_Reports，步骤 6 才覆盖 08_Rescue，故 10_Reports 的 rescue 拷贝滞后一步（5512B / mtime 9-15 23:49 vs 5414B / 9-16 01:08）。脚本既定顺序所致，未改脚本；如需一致需在重生成后再拷一次（超出本轮范围，未做）。
4. **rescue 备份非直接前驱**：`regen_rescue_reports.py` 仅在 `.bak_v661` 不存在时才备份，而该文件自 8-19 起已存在，故 9-15 版本的 rescue 报告**未被备份**。因此 §7 对照 `.bak_v661` 的 diff 混入了 9-15 的既有改动；已用步骤内快照隔离出本轮净效应（6 行）作为补充。
5. **taxonomy_images.json 未刷新**：其写入分支被「`upset_data.json` 已存在」短路，保持 8-20 内容；`pipeline_report.html` 可能内嵌该旧图。属既有行为，非本轮引入。
6. **.ok 标志未被触碰**：`.taxonomy.ok`(8-12)、`.analysis.ok`(8-26)、`.report.ok`(8-25)、`.rescue.ok`(8-19) 均未更新，符合「不动 schema / 不跑 analysis」的约束，但意味着这些 flag 不代表本次刷新状态；本轮未 touch。
7. **待决项（本轮明确不动，仅登记）**：`prevalence_full_table.tsv`（8-25）与 `09_Virome_Analysis` / `09b_Analysis_Verify` 派生表均未变更；`run_analysis` 会给分类表补 `Nucleic_acid` 列，属超出本次范围的 schema 变化，未执行。

## 10. 结论

门槛由 `share > 0.5` 改为 `share >= 0.5`（gate 6.7）后，Aphis 共识产物行数与 contig 集合完全不变，矛盾率 **17.0752% → 0.2859%**（418 → 7 行），gate 自检 PASS，脚本自检与独立复算逐格一致。10_Reports 已刷新（17 文件变化，其中 taxonomy 直接派生项由门槛驱动，其余为上游陈旧回填/环境因素），rescue 报告已刷新（本轮净改动 6 行 + summary 重写）。以上问题均为报告层/环境层既有现象，未影响本次门槛语义的正确重跑。
