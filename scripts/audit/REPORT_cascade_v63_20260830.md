# cascade v6.3 上线前实测报告

数据：`RNA-Lycium_barbarum_out` Votus 层，20,892 contig，7 工具（ACVirus/metabuli/mmseqs 等）
参照：VMR MSL41（`~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv`，18,833 个 (rank,value) 谱系键）
补丁：`scripts/audit/patch/virus_classifier_analysis.R.cascade_v63`
md5：`76e84e4100d129eb37dc60c1ba35c9aa`（本地与服务器 `/tmp/` 一致，LF 无 CRLF）
A/B 产物：`/tmp/rc_legacy/`、`/tmp/rc_cascade/`（未触碰线上任何产物）
对账脚本：`scripts/audit/cmp_products.py`、`probe_stale.py`、`chk_ph.py`

---

## 1. 计时：23.5 分钟假设被实测推翻

| 运行 | 起止 | 耗时 | `build_consensus` 内部 |
|---|---|---|---|
| legacy | 15:52:28 → 15:53:18 | **50 s** | 5.7 s |
| cascade | 15:56 → 15:56:22 | **53 s** | 9.1 s |

分阶段（cascade）：阶段1-2 长表+权重 0.3 s｜阶段3 逐级淘汰 **3.2 s**｜阶段4 0.1 s｜阶段5 主工具 1.0 s｜闸门-由种提属 0.1 s｜科属校准 0.2 s｜跨参照相容 2.6 s｜阶段6 agree+清理 1.5 s。

08-12 日志里的「00:44:54 → 01:08:31 = 23.5 分钟」两次独立重跑都复现不出来，差两个数量级。结论：不是代码热点，与是否开 cascade 无关。

## 2. 三向对账（VMR 口径）

| 指标 | 盘上 08-12 产物 | legacy（现行代码+闸门） | cascade |
|---|---|---|---|
| 科属不匹配（Genus 侧） | 1509 / 14444 = **10.4 %** | 0 / 13535 = **0.0 %** | 0 / 13239 = **0.0 %** |
| 科种不匹配（Species 侧） | 1382 / 12047 = **11.5 %** | 1382 / 12047 = **11.5 %** | 30 / 11047 = **0.3 %** |
| Genus 非空 | 16022 | 14869 | 14692 |
| Species 非空 | 18056 | 17905 | 17310 |
| Family 非空 | 17919 | 17919 | 17512 |
| 占位值残留 | Genus 1 + Species 151 | 0 | 0 |

归因清楚：

- **科属 10.4 % → 0 %** 全部来自现有闸门（`harmonize_family_genus` + 跨参照相容性）。legacy 已经做到 0，cascade 再加 0 收益。
- **科种 11.5 % → 0.3 %** 全部来自 cascade。legacy 的科种数字与盘上旧产物**逐位相同**（1382/12047），说明现行闸门对 Species 侧零作用，这一栏从来没被管过。
- 三份产物 contig 集合完全相同（20892 / 20892，无新增、无丢失）。

逐格 diff（legacy → cascade，仅统计共同 contig）：

| 阶元 | 改变 | 变空 | 由空变有 |
|---|---|---|---|
| Realm | 0 | 0 | 0 |
| Kingdom | 28 | 0 | 0 |
| Phylum | 96 | 1 | 0 |
| Class | 97 | 7 | 0 |
| Order | 491 | 368 | 0 |
| Family | 825 | 407 | 0 |
| Genus | 870 | 341 | 164 |
| Species | 2250 | 595 | 0 |
| 合计 | **4657** | 1719 | 164 |

`Genus 由空变有 164` 是闸门「由种提属」把 cascade 新选中的二名法 Species 派生成属，属侧 VMR 一致性仍为 0 不匹配，风险低，但属于「投票之外的信息通道」，建议抽查几条再定。

旧产物 → legacy 只动了 2163 格：Genus 2012（变空 1153 = 闸门置空 1126 + 占位 27，与日志**逐格吻合**）、Species 151（全是占位值）。

## 3. 我自己补丁里的缺陷（已修）

阶段 4「最完整工具行」兜底在 **cascade 下实测填补 2410 格**：`best_data` 取的是全工具里 completeness 最高的那一行，不区分该工具是否已在粗阶元被淘汰，等于让被淘汰工具从后门复活，与「淘汰不复活」正面冲突。legacy 下恒 0 格（`best_data` 有值的阶元 `wide` 必然也有值），所以它一直是死代码，只在 cascade 下复活。

已改为 cascade 下跳过。重跑后日志变为「最完整工具行兜底: cascade 模式已停用」，且科属校准动作从「替换 308 行 / 置空 1046 行」降到「替换 139 行 / 置空 67 行」（此前是被兜底灌进去的值触发的）。

发现它的唯一途径是那行计数埋点。埋点保留。

## 4. 占位值清理修复（已含在补丁）

原实现只做整串精确匹配，`PLACEHOLDER_TAXA` 里只有 `environmental`、`samples` 两个独立整串，多词形式全部漏网。改为「整串命中 ∪ 限定词子串命中」（限定词只用 `environmental / uncultured / unclassified / unidentified / unassigned / not assigned / no hit / unknown / undefined`，泛词仍走整串，避免误伤 Narnavirus 这类正常属名）。

结果：legacy 清空 Genus 27 + Species 151，cascade 清空 Genus 28 + Species 173，两版产物 8 个阶元占位残留**全部为 0**。

## 5. cascade 残留的 30 条科种不匹配（要你定）

全部是「科来自多工具共识、种来自另一支系」的物种级嵌合，例：

| Family | Species | VMR 里该种的真实科 |
|---|---|---|
| Mimiviridae | Simplexvirus macacinealpha1 | Orthoherpesviridae |
| Mimiviridae | Coltivirus dermacentoris | Spinareoviridae |
| Pithoviridae | Chlorovirus conductrix | Phycodnaviridae |
| Nudiviridae | Coccolithovirus huxleyi | Phycodnaviridae |
| Schizomimiviridae | Lentivirus humimdef1 | Retroviridae |
| Phycodnaviridae | Plazymidvirus ZM41 | Peduoviridae |

现行闸门只有「属-科」约束，没有「种-科」约束。把 `enforce_rank_containment` 的判据扩到 Species（复用同一份 `genus_family_ref`），预计把这 30 条压到接近 0，代价是再置空几十格。这一项独立于 cascade，legacy 下同样能装（legacy 现在这 1382 条也是无人管的）。

## 6. 四项口径问题（要你定改不改）

1. **`*_agree` 列在闸门之前计算**。legacy 下 Genus 1153 行、Species 151 行「值已被置空，但 `_agree` 仍写着 3/3: ACVirus,metabuli,mmseqs」。cascade 下 Genus 95 / Species 173。下游若拿 `_agree` 当置信度，读到的是已经不存在的支持。
2. **`completeness` / `confidence` 同样在闸门之前算**（脚本 649 行，早于 627-632 行的闸门）。legacy 1297 / 20892 行与实际非空阶元数不符，其中 1105 行写着 `completeness=8` 实际只有 7 个阶元；cascade 265 行。
3. **`compute_tool_weights` 第 276 行 `rate > 0` 才乘惩罚因子**：`rate=0` 的工具完全免罚，把「该工具在该阶元没有数据（应为 1.0）」与「该工具全错（应重罚）」混成一个条件，逻辑倒挂。现权重只轻罚 0.5-0.9 的中间态，放过了最差的。
4. **占位值清理只覆盖 Family/Genus/Species**，Realm/Kingdom/Phylum/Class/Order 未清；Kingdom 层也没有 `-virae` 后缀校验。本数据集实测这两个缺口**零影响**（8 阶元残留 0，Kingdom 8 个取值全部以 -virae 结尾、1750 个空），属潜在漏洞而非当前损失。

## 7. 上线与回滚（等你拍板）

1. 备份 `~/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R` → `.bak_cascade_20260830`
2. 覆盖为 `.cascade_v63`，`TAX_GATE_VERSION` 6.2 → 6.3 会自动让 8 个项目重跑该 stage（不会 SKIP）
3. 回滚两条路：恢复备份，或直接 `MMPV_CONSENSUS_MODE=legacy`（不必改文件）
4. 建议顺序：先只重跑 barbarum，人工过一遍 30 条种级残留与 4657 格改变，确认后再铺开 8 个项目

`CASCADE_MIN_SHARE` 当前 0.5，只淘汰「票首占比 > 50 % 时的异议工具」。想更保守可以调到 0.6 / 0.7（淘汰变少、更接近 legacy），需要的话我再跑一轮对照。
