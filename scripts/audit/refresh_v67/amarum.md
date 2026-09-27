# amarum 门槛 v6.7 重跑与独立复核报告

- 数据集 label: `amarum`
- 目录 `D` = `/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_amarum_out`
- 执行器: `/tmp/refresh_v67.sh`（md5 `d832c4868a2b8298d7bfca44be0dcb18`，未修改）
- R 共识脚本: `/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R`
  - md5 `f0683ac6af1a9ee1395a100d128365ca`（与本轮要求一致）
  - 第 96 行 `TAX_GATE_VERSION <- "6.7"`，第 112 行门槛语义 `share >= CASCADE_MIN_SHARE`
- 输入 combined: `$D/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv`（422730 B, 2026-08-12 02:24）
- 复核时间: 2026-09-16 01:0x–01:1x (+0800)
- 个人独立复核脚本: `/tmp/audit_v67.py`（md5 `e8648ab6ad811d2e05650a32090b28d0`）、`/tmp/residual_amarum.py`；判据与参考实现 `/tmp/chimera_multi.py` 一致

---

## 1. 成品：旧 / 新

| 口径 | md5 | 数据行数 | mtime |
|---|---|---|---|
| 旧（线上 v6.6 产物） | `8d26262525fe390b07f5e7bf1515606c` | 957 | 2026-08-12 02:25:22 |
| 新（v6.7 重跑） | `868e4248bac3f878a20d9b31cbaa29b5` | 957 | 2026-09-16 01:08:20 |

- 文件为 `$D/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv`，文件总行 958（含表头），数据行 957。
- 旧 md5 三处一致：改动前 live 文件、`/tmp/refresh_v67/amarum/old.tsv`、`05_Taxonomy/Votus.integrated.bak_v66f_20260916/final_integrated_classification.tsv` 均为 `8d2626…`。
- 注意：背景里给的 `12e8f171…` → `a3ead5b5…`、20892 行是 **barbarum** 的数据，与 amarum 无关，不作为本轮锚点。
- 10_Reports 内的成品副本 md5 = `868e4248bac3f878a20d9b31cbaa29b5`，与 `05_Taxonomy` 下一致。

**集合比对（自写脚本，非计数比对）**：`old_only=0  new_only=0  common=957  identical_set=True`。
即 contig_id 集合完全相同，行数不变。

---

## 2. gate 自检原文

`taxonomy_gate_stamp.tsv`（2026-09-16 01:08:20，cat -A 确认 tab 分隔、行尾正常）：

```
gate_version	script_path	script_md5	ref_file	ref_md5	rows	timestamp
6.7	/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R	f0683ac6af1a9ee1395a100d128365ca	/home/zhangwenda/database/taxonomy/genus_family_ref.tsv	865283b4b1a3a2f3c4c3ac702403c554	957	2026-09-16 01:08:20
```

- `gate_version=6.7` ✓
- `script_md5=f0683ac6af1a9ee1395a100d128365ca` ✓

`taxonomy_gate_check.tsv`：

```
check	n	expected	status
rows_total	957	-	INFO
family_filled	814	-	INFO
genus_filled	676	-	INFO
genus_ref_known	637	-	INFO
dual_ref_consistent	592	-	INFO
dual_ref_conflict	0	0	PASS
single_ref_conflict	0	unfixable_by_design	INFO
species_dual_ref_conflict	0	0	PASS
species_ref_known	568	-	INFO
species_single_ref_conflict	1	unfixable_by_design	INFO
```

- `dual_ref_conflict = 0  PASS` ✓
- 注：`species_single_ref_conflict = 1`，状态为 `unfixable_by_design / INFO`（设计内不可修，非失败）。
- 补充：脚本内联的那句 `grep -e gate_version -e script_md5 …` 只匹配到了表头（数据行不含这两个字面量），故 nohup 日志里那段是空内容；stamp 的真实数据行由我单独读取，如上。

---

## 3. 矛盾率（独立重算）

判据：对每对相邻阶元 (H 高, L 低)，用 `rankedlineage.dmp` 查行内 L 的定型条目、取其 H 阶父级名，与行内 H 比较；L 有参照父级、行内 H 有值时才比，不等即矛盾；一行任一对矛盾即计矛盾行。

| 口径 | 行数 | Realm-Kingdom | Kingdom-Phylum | Phylum-Class | Class-Order | Order-Family | Family-Genus | Genus-Species | 矛盾行 | 矛盾率 |
|---|---|---|---|---|---|---|---|---|---|---|
| 旧 v6.6 | 957 | 1 | 1 | 10 | 17 | 64 | 99 | 48 | **167** | **17.4504%** |
| 新 v6.7 | 957 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | **2** | **0.2090%** |

- 与 `/tmp/refresh_v67/amarum/chimera.txt` 逐项比对：`online(new)` = 1/0/1/0/0/0/0，矛盾 2，0.21%；`v66f(old)` = 1/1/10/17/64/99/48，矛盾 167，17.45%。**完全一致，无出入。**
- 基线锚点校验：在改动前的 live 成品上先独立算过一次 = 167 / **17.4504%**，与约定锚点 17.45% 吻合。
- 改动前后 live 文件与 `old.tsv` 的旧数值两次独立重算一致。
- 改后落到 0.21%，低于"1% 上下"的预期上限。

**新产物残余 2 条矛盾的定位**（`/tmp/residual_amarum.py`）：

| contig_id | 矛盾对 | 详情 |
|---|---|---|
| `CRR2703952_clean_NODE_839_length_1425_cov_2.764793` | Realm-Kingdom | 行 `Realm=Riboviria`，参照 Kingdom 定型 Realm=`Floreoviria` |
| `CRR2703952_clean_NODE_9055_length_563_cov_1.630612` | Phylum-Class | 行 `Phylum=Nucleocytoviricota`，参照 Class 定型 Phylum=`Lenarviricota` |

两者均为旧产物中**同样存在**的矛盾（旧 Realm-Kingdom=1、Phylum-Class=10，新 Phylum-Class 只剩这 1 条），属于单参照、门槛语义无法消除的残余，非本轮引入。

---

## 4. 各阶元取值变化行数（old → new）

| 阶元 | 变化行数 |
|---|---|
| Realm | 3 |
| Kingdom | 3 |
| Phylum | 5 |
| Class | 13 |
| Order | 37 |
| Family | 103 |
| Genus | 71 |
| Species | 145 |

- **Order 有值→空 16 行；Family 有值→空 11 行；反向（空→有值）Order 0 / Family 0。**
- 属/种有值数下降可由 fill 统计交叉印证：Order 826→810、Family 825→814、Genus 735→676、Species 836→797（Realm/Kingdom/Phylum/Class 不变）。
- Phylum↔Class 参照相容性：`both_ok 946 / fixed 10 / broken 1 / both_bad 0`。
- 上述三组数字与 `/tmp/refresh_v67/amarum/cmp.txt` 完全一致。

---

## 5. 10_Reports 刷新情况

- `$D/10_Reports` 目录内文件 mtime 已刷新到 2026-09-16 01:08:0x–01:08:40；report_pipeline rc=0，输出 39M。
- 与 `$D/10_Reports.bak_v66f_20260916` 逐个 `diff -rq`，**内容变化文件 18 个**：

**A. taxonomy 派生（预期随共识变化）— 14 个**
`final_integrated_classification.tsv`、`plant_final_taxonomy.tsv`、`taxonomy_composition.tsv`、`taxonomy_fill_stats.tsv`、`taxonomy_sunburst.html`、`pipeline_report.html`、`classification_sankey.html`、`classification_sankey_plant.html`、`plant_virus_summary.tsv`、`All_plant.viruses_info.tsv`、`HQ_plant_viruses_info.tsv`、`all_plant_viruses_genus_summary.tsv`、`upset_data.json`、`rescue_report.tsv`

**B. 非 taxonomy 但属刷新副作用 — 2 个**
- `stage_summary.tsv`：路径由旧的 `/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_amarum_out/…` 更正为当前根；taxonomy 计数同步更新（Known 836→797、NewSp 46→55、Rank Fill Order 826→810 / Family 825→814 / Genus 735→676）。
- `directory_tree.txt`：目录清单刷新，纳入本轮新建的 `.bak_v66f_20260916` 目录等。

**C. 需注意 / 与门槛无关 — 2 个**
- `filter_summary.tsv`：由旧的真实计数（CRR2703952 1854→197 等）变为 `comb 0 0 0.0`。**根因**：`02b_Filter/*/` 下的候选 fasta 全是**悬空软链**，指向已不存在的 `/home/zhangwenda/goji-virome/01_data/data-test/…`；report_pipeline 忠实读到 0。这是数据迁移遗留，与门槛语义无关（`readlink -f` 与 `test -r` 均显示目标缺失）。
- `host_distribution.tsv`：仅行序变化（`Algae 69` 行位置移动），内容等价。

**D. mtime 被刷新但内容未变（cmp 相同）**
`taxonomy_agreement_stats.tsv`、`taxonomy_consistency_summary.tsv`、`assembly_summary.tsv`、`checkv_summary.tsv`、`cobra_summary.tsv`（这二者来自 per-tool 一致性/阶段摘要，不依赖共识结果）。
另外 `10_Reports` 内的 `*.tsv.bak_kill_20260915` 等历史备份文件未变动。

---

## 6. rescue 报告

对照 `$D/08_Rescue/Plant/` 与旧备份 `.bak_v661`：

| 文件 | 旧 | 新 |
|---|---|---|
| `rescue_report.tsv` | 7393 B, 2026-08-13 20:08:26 | 7591 B, 2026-09-16 01:08:41 |
| `rescue_summary.md` | 603 B, 2026-08-13 20:08:26 | 637 B, 2026-09-16 01:08:41 |

内容变化：
- `rescue_report.tsv`：8 行由 `fail`（note=`<2000bp` / 未通过任何分支）变为 `C / blastn_completo`（note=`BLASTN qcov≥98% 直达`）；branch 分布 fail 68 + C 8（候选 76，拯救 8，10.5%）。
- `rescue_summary.md`：修正了旧版的百分比口径错误（旧"过短 (<2000bp) 75 / 110.3%"→新"68 / 100.0%"，并标注分母为"候选表内 / 占未拯救"）。

**需注意**：`$D/10_Reports/rescue_report.tsv`（mtime 2026-09-15 23:49，本次未刷新）与 `08_Rescue/Plant/rescue_report.tsv` md5 不同，差异恰为 2 行 genus 列（`10_Reports` 副本仍有 `Biavirus`，08_Rescue 新报告为空白）。原因：执行器先跑 report 层、后跑 rescue 重生成，而 report_pipeline 用 `shutil.copy2` 复制（保留源 mtime），所以 `10_Reports` 里拷到的是**重生成之前**的上一版。该 2 行的属名变化已由成品印证：`NODE_1623`、`NODE_5886` 旧 `Genus=Biavirus`，新 `Genus=NA`（Species 亦改为 Guapo/Berere partitivirus）。

---

## 7. 备份路径

| 内容 | 路径 |
|---|---|
| 旧 Votus.integrated 整目录 | `$D/05_Taxonomy/Votus.integrated.bak_v66f_20260916` |
| 旧 10_Reports 整目录 | `$D/10_Reports.bak_v66f_20260916` |
| 旧成品单文件 | `/tmp/refresh_v67/amarum/old.tsv` |
| rescue 旧版（早于本轮，已存在） | `$D/08_Rescue/Plant/rescue_report.tsv.bak_v661` / `rescue_summary.md.bak_v661` |
| 本轮日志 | `/tmp/refresh_v67/amarum/`（r.log, chimera.txt, cmp.txt, report.log, rescue.log, manifest.tsv） |

未删除任何文件；未改动 R/py 管线脚本；未触碰其他数据集。

---

## 8. 耗时

- 执行器自身：START 2026-09-16 01:07:43 → DONE 2026-09-16 01:08:41，约 **58 s**（含 R 重跑、gate、对照、report、rescue 全流程）。
- 启动（`nohup` 下发）到 DONE 约 1–2 分钟。

---

## 9. 问题与不确定处（待决）

1. **filter_summary.tsv 归零**（悬空软链）——与门槛无关，但会使 10_Reports 里该表失真；建议后续统一 02b_Filter 软链目标（本轮未动，只记录）。
2. **10_Reports/rescue_report.tsv 落后一代**——执行器固定"先 report 后 rescue"顺序 + `copy2` 保留 mtime 所致；若要求报告层一致，需在 rescue 重生成后再补一次 report，或重排步骤（本轮未改脚本）。
3. **旧 10_Reports 基于旧根生成**（`data-2026/data-test`）——本轮 report 刷新把路径与计数一并纠正，故 diff 中的路径类变化属"修正"而非"回归"。
4. **待决项（本轮明确未动）**：`prevalence_full_table.tsv` 与 `09_Virome_Analysis` / `09b_Analysis_Verify` 派生表未参与本轮刷新；另未运行 `virome_pipeline.py --stage analysis`（避免就地补 `Nucleic_acid` 列的 schema 变化）。
5. **新产物仍有 2 条参照矛盾**（NODE_839 Realm、NODE_9055 Phylum），单参照、设计内不可修，非本轮引入。
6. `species_single_ref_conflict = 1`（INFO/unfixable_by_design），属 gate 自检的已知类别，非失败。
7. **并发环境**：巡检时同机另有 `onekp`、`ruthenicum` 的 refresh_v67 进程在跑（其他会话）；共享只读 `rankedlineage.dmp`，未见对 amarum 产生干扰，仅记录。

---

## 10. 结论

amarum 本轮重跑在版本与语义上均达标：gate 6.7 + 正确脚本 md5，`dual_ref_conflict PASS`；行数与 contig 集合完全不变；矛盾率 17.4504% → 0.2090%（167 → 2 行），7 对阶元分布逐项与参考实现一致；代价是 Order 16 行、Family 11 行的分类信息被淘汰为空（反向为 0）。报告层与 rescue 层已刷新，遗留 2 处需人工留意的报告层问题（filter_summary 软链、10_Reports rescue 副本落后一代），均与门槛语义无关，已在上面单列。
