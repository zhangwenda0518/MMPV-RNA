# onekp（onekp-virus）门槛 v6.7 刷新复核

- 数据集根：`/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus`
- 时间：成品与 gate 戳 `2026-09-16 01:18`，10_Reports 派生与 rescue 落在 `01:37`–`01:40`
- 脚本：`virus_classifier_analysis.R` md5 `f0683ac6af1a9ee1395a100d128365ca`（gate 6.7）
- 参照表：`genus_family_ref.tsv` md5 `865283b4b1a3a2f3c4c3ac702403c554`
- 日志：`/tmp/refresh_v67/onekp/{r.log,report.log,rescue.log,chimera.txt,cmp.txt,old.tsv,manifest.tsv}`

本数据集 535,103 行，是 8 个数据集中最大的一个（占合计 588,953 行的 91%），因此单列一份复核。

---

## 1. 成品

| 项 | 值 |
|---|---|
| 行数 | 535103（含表头 535104），旧新一致 |
| 旧 md5 | `028e85b2aeea9a3977fea8566b6d72c7` |
| 新 md5 | `02a2b8446ee748721e4c5b6c60f79bf0` |
| contig_id 集合 | 完全相同、无重复、顺序一致 |
| 文件大小 | 175,998,835 bytes |

### gate 自检原文（`taxonomy_gate_stamp.tsv` / `taxonomy_gate_check.tsv`）

```
6.7  f0683ac6af1a9ee1395a100d128365ca  865283b4b1a3a2f3c4c3ac702403c554  rows 535103  2026-09-16 01:18:35
rows_total               535103  -                      INFO
family_filled            409038  -                      INFO
genus_filled             257651  -                      INFO
genus_ref_known          235296  -                      INFO
dual_ref_consistent      216905  -                      INFO
dual_ref_conflict             0  0                      PASS
single_ref_conflict         300  unfixable_by_design    INFO
species_dual_ref_conflict     0  0                      PASS
species_ref_known        185098  -                      INFO
species_single_ref_conflict 223  unfixable_by_design    INFO
```

`dual_ref_conflict = 0 PASS`、`species_dual_ref_conflict = 0 PASS`，两项按设计不可修的单一参照冲突 300 / 223 属 INFO，不阻断。

---

## 2. 分类矛盾率（R1 口径，参照谱系库 `rankedlineage.dmp`）

| 口径 | Realm-Kingdom | Kingdom-Phylum | Phylum-Class | Class-Order | Order-Family | Family-Genus | Genus-Species | 矛盾行 | 矛盾率 |
|---|---|---|---|---|---|---|---|---|---|
| 旧 | 1177 | 14586 | 14276 | 29935 | 10195 | 23612 | 572 | 71295 | 13.32% |
| 新 | 375 | 0 | 304 | 48 | 192 | 300 | 291 | 1460 | 0.27% |

- 独立实现（`/tmp/indep_chimera.py`）与本数据集 `chimera.txt` 逐格一致，已由第三方交叉验证脚本自动比对确认（16 行全格 `ALL_CELLS_MATCH`）。
- R2 口径对照：旧 70488（13.17%）→ 新 915（0.17%）。
- 各对之和（1510）大于矛盾行数（1460），因同一行可同时计入多对。

## 3. 阶元取值变化与信息量净损失

| 阶元 | 变化行数 | | 有值→空 | 空→有值 |
|---|---|---|---|---|
| Realm | 18399 | Order | 23184 | 0 |
| Kingdom | 18688 | Family | 25995 | 0 |
| Phylum | 5037 | | | |
| Class | 18380 | | | |
| Order | 40262 | | | |
| Family | 39530 | | | |
| Genus | 60148 | | | |
| Species | 68985 | | | |

Order 与 Family 的填充率各降约 4.3 与 4.9 个百分点（23184 / 25995 / 535103）。这是 `>=` 门槛的已知代价，本数据集承担了其中约 87%（23916 中的 23184、27037 中的 25995），原因是它的行数基数最大。

## 4. Phylum↔Class 参照相容性

| both_ok | fixed | broken | both_bad |
|---|---|---|---|
| 520758 | 14041 | 69 | 235 |

`fixed : broken` = 14041 : 69，收益远大于代价。与新产物的 Phylum-Class 残余 304 行、旧产物 14276 行自洽。`both_bad` 235 行属参照库本身与行内双方都不相容的情形，是 R 头注里写明的「参照库口径与工具体系不相容」类，门槛动不了。

## 5. rescue 层

| 项 | 数值 | 占比 |
|---|---|---|
| 候选 | 20517 | 100% |
| A（CheckV ≥90% 蛋白完整） | 0（另含 130 条免拯救 prepass，不在候选表内） | 0.0% |
| B（VSI reads 延伸 + genus 回退） | 50 | 0.2% |
| C（BLASTN + RGA/ragtag） | 9212 | 44.9% |
| D（genus_len 属水平兜底） | 0 | 0.0% |
| 合计（B+C） | 9262 | 45.1% |
| 未拯救 | 11255 | 54.9% |
| 最终无冗余 vOTU | 8717 | |

未拯救原因：过短（<2000bp）10329（91.8%）、属组处理失败 924（8.2%）、未通过任何分支判定 1、无属分类信息 1。

`08_Rescue/Plant/rescue_report.tsv` mtime `2026-09-16 01:40`，已随本轮刷新。与门槛的因果关系同其他数据集：rescue 报告的翻转主要来自生成器策略，本轮只是顺带重跑。

## 6. 10_Reports 逐文件归因

`diff -rq` 共 13 个文件内容变化。按成因分四类：

### 6.1 门槛直接派生（6 个）

`final_integrated_classification.tsv`（mtime 01:18）、`taxonomy_composition.tsv`（01:37）、`taxonomy_fill_stats.tsv`、`taxonomy_sunburst.html`（01:37）、`plant_virus_summary.tsv`、`All_plant.viruses_info.tsv`

### 6.2 9/15 上游拷贝（2 个，与门槛无关）

`all_plant_viruses_genus_summary.tsv`、`HQ_plant_viruses_info.tsv`，两者 mtime 均为 `2026-09-15 14:26`。这是 `report_pipeline.py` 用 `copy2` 从 `09_Virome_Analysis` 拷来的产物（copy2 保留源 mtime），而 `09_Virome_Analysis` 本轮未重跑。它们出现在 diff 里，是因为备份是刷新前一刻的快照，而这两份是 9/15 才刚被拷进来的。

### 6.3 rescue 重生成（1 个）

`rescue_report.tsv`（01:40）

### 6.4 生成器 / 环境（4 个）

`directory_tree.txt`（纳入新备份目录）、`pipeline_report.html`、`stage_summary.tsv`（路径更正）、`filter_summary.tsv`（全 0，断链软链，见总表 §6）

### 6.5 应刷未刷（3 个，本数据集独有，需补跑）

| 文件 | mtime | 根因 |
|---|---|---|
| `classification_sankey.html` | 2026-09-12 15:41 | `report.log` 原文：`[WARN] Sankey 生成失败: ... taxonomic_sankey.py ... timed out after 120 seconds`，随后 `Done — 120s` |
| `classification_sankey_plant.html` | 2026-09-12 15:41 | 同上 |
| `plant_final_taxonomy.tsv` | 2026-09-12 15:41 | 同一批次未重写 |

即 `taxonomic_sankey.py` 对本数据集（535k 行、输入 176MB）在 120 秒硬超时内跑不完，脚本按 WARN 降级跳过，两个 sankey 停在 9/12 的 v6.6 前版本。goji 的 7 个数据集体量小，sankey 正常重生成，所以此问题只在 onekp 暴露。

`plant_final_taxonomy.tsv`（7,046,137 bytes，9/12）同样是 v6.7 前版本，需确认它的生成入口是否与 sankey 同批。

`upset_data.json` 内容与备份相同（mtime 01:08，早于 gate 戳 01:18），未随本次 taxonomy 变化而改变，也值得单独确认生成器是否被调用。

---

## 7. 备份与未动

- 备份：`05_Taxonomy/Votus.integrated.bak_v66f_20260916/`、`10_Reports.bak_v66f_20260916/`（备份名里的 `v66f` 是执行器硬编码，内容为刷新前线上成品，不是 v6.6f 候选）
- 未动：`prevalence_full_table.tsv`、`09_Virome_Analysis/**`、`09b_Analysis_Verify/**`、`.ok` flag；未跑 `--stage analysis`；未引入 `Nucleic_acid` 列
- 未删任何文件

## 8. 本数据集新增待决项

| # | 事项 | 说明 |
|---|---|---|
| 10 | onekp 的两个 sankey 与 `plant_final_taxonomy.tsv` 补跑 | 需提高 `taxonomic_sankey.py` 超时或改分段/抽样生成，否则该数据集报告层永久停在 v6.6 前版本 |
| 11 | `upset_data.json` 生成器是否被调用 | 内容与刷新前相同，需确认是数值巧合还是未执行 |
