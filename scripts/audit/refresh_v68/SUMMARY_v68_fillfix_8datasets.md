# 8 数据集 fill 层修复全量重跑汇总（gate 6.7 + fill-fix1，2026-09-16）

配套根因报告：`REPORT_fill_rootcause_20260916.md`（同目录）。

## 0. 结论

上游标准化 fill 层的两个语义错误（同名歧义劫持 + 覆盖式写入 + 重复回填）已修复并全量重跑 8 数据集。family 级参照一致性冲突 314 → 0，「Phycodnaviridae 被写成 Adenoviridae」类错误定点消除；成品行数合计 588,953 → 589,441（+488），非 NA 单元格净增 6,826。878 个「值→NA」丢失格里 25 格是撤销旧 fill 伪造值，853 格由下游既定规则/门槛/参照库缺口产生（其中 22 格属可修的占位值不回落）。另发现一个独立缺陷：`KNOWN_REALMS` 白名单漏收 ICTV 2025 新 realm，静默清空 12,797 行 Realm。

## 1. 主表

| 数据集 | 成品 md5（新） | 行数（旧→新） | gate | gate 戳时间 |
|---|---|---|---|---|
| Alternaria alternata | `17cb6c36d51d350ce9cd8ebb7c62fa30` | 1728 → 1729 (+1) | 6.7 | 2026-09-16 05:44:52 |
| Aphis gossypii | `e7757d0cf1cb714368ac20d43941ef5c` | 2448 → 2456 (+8) | 6.7 | 05:46:27 |
| Fusarium nematophilum | `9212a871d5c101b284adb132834de985` | 1637 → 1640 (+3) | 6.7 | 05:48:07 |
| Lycium amarum | `8bfe99a29e9c7763c2fb8febb5af8495` | 957 → 957 (0) | 6.7 | 05:49:37 |
| Lycium barbarum | `b5cc2aee0891ca31f2443db9e9971390` | 20892 → 20935 (+43) | 6.7 | 05:51:50 |
| Lycium chinense | `0a2f6c7f05f8ee86bd3335361e48166a` | 11007 → 11020 (+13) | 6.7 | 05:55:43 |
| Lycium ruthenicum | `95a79b9a6cff1cf1a53ee0b79f76c422` | 15181 → 15219 (+38) | 6.7 | 05:58:19 |
| onekp（1KP） | `d46cbf1bfec0d3d67b1f0af3c6d5f400` | 535103 → 535485 (+382) | 6.7 | 06:17:40（R 步） |
| 合计 | | 588,953 → 589,441 (+488) | | |

版本三件套（8 数据集一致）：gate `6.7`；R 脚本 `virus_classifier_analysis.R` md5 `f0683ac6af1a9ee1395a100d128365ca`（1439 行，`CASCADE_MIN_SHARE` 默认 0.5）；Python 层 `virus_classifier.py` 48433 B / mtime 2026-09-16 05:39（fill-fix1）。

成品列数 20 未变，键（contig_id）唯一，无新增/丢失键类型（`verify_pairs.sh` 8 数据集全跑，日志 `/tmp/vp_*.log`）。

**收尾核对（2026-09-16 11:47，`/tmp/status_v68.sh`）**：8 数据集 gate 戳行数与成品行数逐一相等（1729 / 2456 / 1640 / 957 / 20935 / 11020 / 15219 / 535485）；onekp 报告步 06:49:37 结束（rc=0，`10_Reports` 419M，`pipeline_report.html`、两张 sankey、`plant_final_taxonomy.tsv` 均为 06:46–06:49 本轮产物）；刷新进程全部退出，无残留。

## 2. 单元格变化（按键对比，非按行 zip）

| 数据集 | 值→NA | NA→值 | 值→不同值 |
|---|---|---|---|
| Alternaria | 4 | 7 | 184 |
| Aphis | 9 | 44 | 375 |
| Fusarium | 7 | 14 | 224 |
| amarum | 2 | 10 | 139 |
| barbarum | 81 | 125 | 2751 |
| chinense | 50 | 48 | 1235 |
| ruthenicum | 62 | 101 | 2169 |
| onekp | 663 | 7355 | 84283 |
| 合计 | **878** | **7704** | **91360** |

## 3. 丢失格（878 格）归因

### 3.1 机制轴（`explain_lost_cells.py`）

| 桶 | 含义 | goji | onekp | 合计 | 性质 |
|---|---|---|---|---|---|
| A | 祖先链不相容：值在参照库中的祖先与行内粗阶元冲突 | 172 | 487 | 659 | 设计使然（逐级相容性约束） |
| A2 | 参照库对该名无祖先记录 | 7 | 152 | 159 | 参照库缺口 |
| B | 占位/伪名串取胜后不回落 | 20 | 2 | 22 | **可修（待决 ②）** |
| C | 逐级淘汰（门槛 `share >= 0.5` 的直接后果） | 15 | 5 | 20 | 设计使然（门槛语义） |
| D | 未解释 → 手工分解 | 1 | 17 | 18 | 全部归因，见下 |

D 桶明细：chinense 1 格 = Realm 白名单（清 Floreoviria）；onekp 17 格 = 6 格 Realm 白名单 + 4 格票面为空（旧值是 fill 构造）+ 7 格真种 vs `environmental samples` 1:1 平票、占位串取胜后清空（B 的同机制变体）。

### 3.2 投票背书轴（`probe_lost_backing.py`：旧值是否真有工具投票）

| 数据集 | 丢失格 | 真丢失（旧值在新 combined 仍被投票） | 撤销 fill 伪造值 |
|---|---|---|---|
| 7 个 goji | 215 | 213 | 2（barbarum Family 1 格票面为空 + chinense Realm 1 格仅旧票面有值） |
| onekp | 663 | 640 | 23（17 格票面为空 + 6 格仅旧票面有值） |
| 合计 | 878 | **853** | **25** |

即 fill 修复本身没有丢真数据；25 格「丢失」是撤销旧 fill 写入的伪造值。

## 4. 收益与代价

### 4.1 收益

- family 级 `single_ref_conflict`：修前合计 314（onekp 300 / barbarum 8 / chinense 3 / ruthenicum 2 / Fusarium 1）→ 修后 **8 数据集全 0**。
- 定点：onekp `Organic Lake phycodnavirus` 行 Family 由 `Adenoviridae`（300 行）改为 `Phycodnaviridae`（该 Genus 下现 674 行 Family=Phycodnaviridae），与参照表一致。
- `dual_ref_conflict`、`species_dual_ref_conflict`：前后均 0 PASS。
- 单元格净增 7704 − 878 = **+6826** 个非 NA 值（NA→值 7704 与值→NA 878 是同一对比口径的两个方向；基线为 gate 6.7 产物 `.bak_fillfix_20260916`）。

### 4.2 代价 / 需观察

- `species_single_ref_conflict`：修前 244 → 修后 **247**（Alternaria 0 / Aphis 2 / Fusarium 1 / amarum 1 / barbarum 10 / chinense 4 / ruthenicum 4 / onekp 225）。属残余类，与 v67 结论一致（种级占位/伪名与参照库种名质量差异），非本轮引入。
- 相邻阶元矛盾行（`chimera_rows.py`，基线同为 `.bak_fillfix_20260916`）：goji 合计 213 → 262（Alternaria 6→8、Aphis 7→12、Fusarium 3→5、amarum 2→6、barbarum 87→102、chinense 57→63、ruthenicum 51→66）；onekp 1460 → **1344**（净 −116）。逐对计数 goji 220→267、onekp 1510→1393（同一行可命中多对，故不等于矛盾行数）。该指标双向变化：修正粗阶元会暴露原本被劫持值「吸收」的细阶元错配，也会消除劫持造成的错配。goji 净 +49 的性质待逐条判（待决 ⑦）。
- 行数 +488：旧版因整块被写成同一谱系而产生行合并，修复后还原独立行。

## 5. 新发现：`KNOWN_REALMS` 白名单静默清空 12,797 行

详见根因报告 §5。要点：R 脚本 L33 白名单 6 个旧 realm，L215 把白名单外的 Realm 一律清成 NA；参照库 `rankedlineage.dmp`（2026-05-23）病毒侧 realm 为 Riboviria / Duplodnaviria / Floreoviria / Varidnaviria / Ribozyviria / Adnaviria / Singelaviria，**Floreoviria 14,881 条与 Singelaviria 62 条不在白名单**，而已被 2025 释放淘汰的 Monodnaviria 仍在白名单里。受影响成品行 12,797（含 onekp 11,924）。数据里实际出现的越界值只有 Floreoviria。

## 6. 待决项（需拍板）

① **`KNOWN_REALMS` 补新 realm**（R 一行常量）：建议 `+ Floreoviria + Singelaviria`，Monodnaviria 保留（工具自带库仍在用）。影响面 12,797 行 Realm 由 NA 恢复为真值（含 6 格丢失格）。改动需与「只重跑 R」统一执行。

② **B 机制补回落**（R 层）：占位/伪名串票数最高或并列取胜后不回落，直接留 NA。goji 20 格 + onekp 约 9 格（2 格票首占位 + 7 格平票取胜）。修法：占位串胜出时按剩余候选依次尝试，都不可用才 NA；或平票时优先非占位候选。

③ **Genus 列混入物种级名**（R 层）：Genus 含空格的行数（新成品实测）：Alternaria 0 / Aphis 23 / Fusarium 1 / amarum 2 / barbarum 69 / chinense 21 / ruthenicum 73 / onekp 1112 = **1301 行**。典型值 `Cotesia sesamiae bracovirus`（onekp 475+165）、`Organic Lake phycodnavirus`（onekp 311）、`Human papillomavirus`（onekp 122）、`Small anellovirus`。ICTV 属名均为单词，建议加规则：Genus 含空格且非参照库属名 → 清 NA（或下沉到 Species 列）。注意这会与 `Organic Lake phycodnavirus` 的 Family 修正产生联动（其 Genus 位置本无 ICTV 属）。

④ 报告层：`report_pipeline.py` 步骤 5 `copy2` 早于步骤 6 的重生成，导致 `10_Reports` 与 `08_Rescue` 的 `rescue_report.tsv` 不同步。本轮实测：8 个数据集里 **6 个一致**，`barbarum`（10_Reports 01:04 vs 08_Rescue 05:54）与 `onekp`（01:40 vs 06:26）两个整整落后一版。修法二选一：改脚本步骤顺序后只补跑这两个数据集的 report 步，或直接以 08_Rescue 为准同步两份文件。

⑤ `prevalence_full_table.tsv`、`09_Virome_Analysis/**`（补 `Nucleic_acid` 列）、`calibration_20260914/` 是否重跑/纳入成品。（HMM 判据已补齐，8/8 全部有 09b 证据，见 §7.2）

⑥ `REPORT_chimera_multi §6` B/C 方案是否上线。

⑦ goji 相邻阶元矛盾行净 +49 的逐条性质判定（是否全为「修正粗阶元后暴露的细阶元错配」）。

⑧ **rescue 口径定调**：论文 §3.5 现文 284 / 294 / 130 与三套实测口径都不等；建议正文改用去冗余 HQ vOTU（276 / 272 / 122，barbarum / ruthenicum / chinense）并写明限定语。四套口径出处见 §7.2。

## 7. 论文同步建议

### 7.1 §2.3（分类共识与门槛）

现稿若描述 weighted voting（权重 × 工具偏好 × 检出率 × 种名质量），需改写为当前的 count 制：

- 逐级投票按工具计数（`RANK_DEPTH_WEIGHTS` 在阶元内恒为常数、且 count 模式下已停用；`TOOL_BIAS`、检出率校准、种名质量降权同停用，仅保留为历史分支）；
- 票首占比 `>= 0.5` 即淘汰持异议工具（2026-08-30 拍板，`CASCADE_MIN_SHARE = 0.5`，环境变量可覆写）；
- `NA` 不计为异议票，也不参与计数；
- 平票两级裁决：先看下一（更细）级一致性，再按 `TIE_BREAK_ORDER`（ACVirus > CAT > VITAP > diamond_lca > genomad > metabuli > mmseqs）兜底；
- 逐级相容性约束：细阶元值在参照谱系上必须与行内粗阶元相容，否则清 NA（这是本轮 659 格「值→NA」的主因，属设计行为，应在方法里写明）；
- 上游标准化层：按「同名多谱系零冲突 + 唯一最优才回填 + 只填 NA 不覆盖」校验，避免 NCBI 同名歧义劫持（fill-fix1）；
- 版本戳：正文/补充材料写 `gate_version 6.7` + R 脚本 md5 `f0683ac6af1a9ee1395a100d128365ca`，成品 20 列、键 contig_id。

### 7.2 §3.5（rescue 数字）

本轮 refresh 同步重生成 `08_Rescue/Plant/rescue_summary.md`（barbarum 09-16 05:54 / ruthenicum 06:00 / chinense 05:56），数字与上一轮一致，即 rescue 计数不受 fill 修复与门槛影响。四套口径（barbarum / ruthenicum / chinense）：

| 口径 | barbarum | ruthenicum | chinense | 出处 |
|---|---|---|---|---|
| 论文现文 | 284 | 294 | 130 | 更早 vintage，计数定义未注明 |
| B+C（候选表内） | **286** | **295** | **124** | `08_Rescue/Plant/rescue_summary.md` |
| 去冗余 HQ vOTU | **276** | **272** | **122** | 同上（末行「最终无冗余 vOTU」） |
| KEEP+REVIEW（09b，含 HMM） | **227** | **228** | **84** | `09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv`（hmm_rescue 生效 42 / 43 / 28） |
| A+B+C（8/19 版） | 308 | 323 | 133 | `rescue_report.tsv.bak_v661` |

建议正文用「去冗余 HQ vOTU」并在句中写明限定语（rescued contigs, deduplicated），补充材料给全口径与生成时间。A 支（CheckV ≥90%）已移出候选表、单列为免拯救 prepass（22/28/9 条），这解释了 308→286 / 323→295 / 133→124。

**09b 侧 HMM 判据现已覆盖 8/8 数据集**（此前只有 barbarum）：`hmm_ct3_hits.tsv` 与 `rescue_evidence_scored.tsv` mtime 分别为 barbarum 09-15 00:50 / 其余 7 个 09-16 03:56–06:26（onekp 06:26）。命中数：Alternaria 29、Aphis 31、Fusarium 37、amarum 9、barbarum 311、chinense 134、ruthenicum 312、onekp 8916。因此 KEEP+REVIEW 口径的 ruthenicum（217→228）与 chinense（76→84）数字已变，与 pre-HMM vintage 不可互比；若采用该口径须从刷新后的 `rescue_evidence_scored.tsv` 重取。

## 8. 证据索引

- 根因 diff：`cd ~/MMPV-RNA/virome_discovery_pipeline && diff -u virus_classifier.py.bak_fillname_20260916 virus_classifier.py`
- 刷新日志：`/tmp/refresh_v68_all.log`（每数据集一段：备份路径 / 矛盾行表 / report 步）
- 键对比日志：`/tmp/vp_*.log`（8 数据集）
- 背书探针输出：`/tmp/backing_onekp.log`
- 探针脚本本地：`D:\桌面\延伸基因组\MMPV-RNA\scripts\audit\fix_fill\`
  `{verify_pairs.sh, explain_lost_cells.py, probe_lost_cells.py, probe_lost_backing.py, probe_d_cells.py, probe_realm_whitelist.py, probe_realm_inventory.py, count_genus_multiword.py, diff_prod_by_key.py, chimera_rows.py}`
- 备份：成品 `<ds>/05_Taxonomy/Votus.integrated.bak_fillfix_20260916/`（v67 产物）、`.bak_v66f_20260916/`；combined `Votus_combined_taxonomy.tsv.bak_fillfix_20260916`；报告 `10_Reports.bak_fillfix_20260916`
- 上一轮汇总：`scripts/audit/refresh_v67/SUMMARY_v67_8datasets.md`
