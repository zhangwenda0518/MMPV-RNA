# 分类共识 v6.4 验收报告（barbarum 单项目实测）

日期：2026-08-30
补丁：`scripts/audit/patch/virus_classifier_analysis.R.cascade_v64`
当前文件 md5 本地 `8E51CF21B8EF801A1539E09272DD0CB1` = 服务器 `/tmp` 副本 `8e51cf21b8ef801a1539e09272dd0cb1`（Rscript parse 通过）
上一版 md5 `7ADE86FE97E290CA22EF31F71431B0B9`（已由复核修正版取代，差异仅为注释与两条日志标签，见第七节）
对照基线：`REPORT_cascade_v63_20260830.md`（v6.3 只做了科-属级，本版补齐 5 个口子）
输入：`RNA-Lycium_barbarum_out/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv`（20892 contig × 7 工具）
产物：`/tmp/rc64b_legacy/`、`/tmp/rc64b_cascade/`（由当前 md5 的补丁跑出；上一版同输入产物 `/tmp/rc64_*` 数据逐字节相同）
日志：`/tmp/rc64b_legacy.log`、`/tmp/rc64b_cascade.log`；执行脚本 `scripts/audit/rc64b_run.sh`

## 一、这一版修的 5 件事，逐条实测

| # | 问题（v6.3 实测） | v6.4 改法 | 实测结果 |
|---|---|---|---|
| 1 | 种级无人管：科-属校准会把自洽工具报的属换进来，Species 里残留的旧属主张谁也看不到（科种不匹配 1382 行 / 11.5%） | 新增 `enforce_species_containment()`：取 Species 首词作属名，双参照（NCBI + VMR）定型科都与行内 Family 互斥时置空 Species | legacy 置空 1283 行、cascade 58 行；科种不匹配 1382 → 213（2.0%）→ 2（0.0%） |
| 2 | `*_agree` 是闸门前的口径：值被置空仍声称「N 个工具支持」（legacy Genus 1153 行 / Species 151 行） | 整段挪到收尾清理之后，对**最终值**重算；`n_tot` 口径不变 | 探针 [A] 两种模式均 **0 行** |
| 3 | `completeness`/`confidence` 在阶段 4 就算死，之后闸门/清理还在改值（不符 1297/20892 行） | 末端唯一赋值点重算 | 探针 [B] **不符 0 / 20892**（两种模式） |
| 4 | 权重与表现倒挂：`if (rate > 0)` 让 rate == 0 的工具免罚拿全额，rate = 0.01 的反而被削到 1% | 有 rate 就乘 `max(rate, RATE_FLOOR=0.25)`，只有 NULL/NA 才视作无数据 | **7 格触底**（VITAP Family/Genus/Species、genomad Genus/Species、mmseqs Genus/Species）；权重落盘 `tool_weights.tsv` 可追溯 |
| 5 | 占位值只清 Family/Genus/Species 三阶元；Kingdom 无格式校验 | 清理扩到全 8 阶元；Kingdom 强制 `-virae$` 结尾，其余阶元只记后缀诊断 | 8 阶元占位残留 **全 0**；Kingdom 合规 **0 置空**（现有取值本已全部合规）；Genus 后缀异常 16 格 / cascade 15 格仅记录 |

修复 4 的反推证据（`probe_weights.py`，由落盘权重除以 base = RANK_DEPTH_WEIGHTS × TOOL_BIAS 反解 rate）：

```
  VITAP Genus  base=70.4  实际=17.6  -> rate <= 0.25（旧代码 rate=0 时会拿满 70.4）
  mmseqs Genus base=64.0  实际=16.0  -> rate <= 0.25（旧代码同样拿满）
  genomad Species / VITAP Species / mmseqs Species / VITAP Family / genomad Genus 同样触底
```

## 二、四方对账（VMR MSL41，18833 谱系键）

| 指标 | 线上旧产物 | v6.3 legacy | v6.4 legacy | v6.4 cascade |
|---|---|---|---|---|
| 科-属不匹配（Genus 侧） | 1509 / 14444（10.4%） | 0 / 13535 | 0 / 13531 | 0 / 13239 |
| 科-种不匹配（Species 侧） | 1382 / 12047（11.5%） | 1382 / 12047（11.5%） | 213 / 10905（2.0%） | **2 / 11047（0.0%）** |
| 占位残留 | Genus 1 + Species 151 | 0 | 0 | 0 |
| Genus 有值 | 16022 | 14869 | 14868 | 14694 |
| Species 有值 | 18056 | 17905 | 16623 | 17259 |
| 科-属自检 | 无 | 0 | 0 PASS | 0 PASS |
| 科-种自检 | 无 | 无 | 0 PASS | 0 PASS |

逐格 diff：

- 旧产物 → v6.3 legacy：Genus 2012 格（其中变空 1153）、Species 151 格（全部是占位值清理）
- v6.3 legacy → v6.4 legacy：Family 18、Genus 52（变空 7 / 由空变有 6）、Species 1394（变空 1282）。前两项来自修复 4（权重变了，科-属校准从「替换 959 / 置空 1126」变到「替换 974 / 置空 1130」），Species 来自修复 1
- v6.4 legacy → v6.4 cascade：合计 4622 格。属以上各阶元 28 至 491 格，Family 828（变空 407）、Genus 886（变空 340 / 由空变有 166）、Species 2196（变空 394 / 由空变有 1030）。cascade 的净效应是「严进严出」：置空更多，但由空变有也显著增加（低票工具被淘汰后，剩余工具的票首反而能填上空缺）

## 三、残留在哪里，为什么不动

科种不匹配残留 213 行（legacy）全部落在「单侧判据」：

```
Mimiviridae          Orthopestivirus bovis        VMR 科=['Pestiviridae']
Phycodnaviridae      Orthohepacivirus hominis     VMR 科=['Hepaciviridae']
Mimiviridae          Alphahydrivirus permafrostis VMR 科=['Hydriviridae']
```

即 VMR 认识这个属、NCBI 不认识（或反是），只有一张参照表能发声。这与科-属闸门是同一取舍：单侧不一致多半是分类学版本漂移（如 Alphasatellitidae 改属后科名未同步），贸然置空等于拿一张表的漂移去覆盖另一张表。cascade 下这类只剩 2 行。

参数化同类项（`probe_sp_gate.py`）：cascade 残留 11 行是「双参照彼此不一致且只有一个与行内相符」，2 行单侧，3130 行属名不在参照表，2 行种名非双名法。**现判据可治 0 行，加严到 A∪B 也是 0 行**，说明种级已经做到底，不需要再加码。

另外一个可解释的口径：legacy 保留 16623 个 Species 值，cascade 保留 17259 个，cascade 反而更多。原因是 cascade 在 Genus 层淘汰了低票工具，Family/Genus 的共识值更集中，「Species 首词与行内 Family 相容」的通过率随之提高。

## 四、自检台账（跑完即产，随产物落盘）

`taxonomy_gate_check.tsv`（两模式均 PASS）：

```
v6.4 legacy     v6.4 cascade    含义
rows_total                  20892            20892
family_filled               17919            17512
genus_filled                14868            14694
genus_ref_known             13918            13589
dual_ref_consistent         12829            12659
dual_ref_conflict               0 PASS           0 PASS   科-属级：双参照都反驳的行
single_ref_conflict            11                8        科-属级：单侧，设计内不动
species_dual_ref_conflict       0 PASS           0 PASS   科-种级：同口径
species_ref_known           11786            12164
species_single_ref_conflict   235               13        科-种级：单侧
```

`species_containment_blanked.tsv` 同步落盘（legacy 1283 行 / cascade 58 行，含表头），逐行可复核。

## 五、待拍板的部署决策

1. **cascade 是否上线**：科种不匹配从 213 行降到 2 行，Family 侧变空 407 格、Species 由空变有 1030 格。上线 = `MMPV_CONSENSUS_MODE` 保持默认 cascade 并把 `TAX_GATE_VERSION` 提到 6.4（管线会按戳强制重跑 05，8 个项目各自约 1 分钟）。
2. **`CASCADE_MIN_SHARE = 0.5` 是否维持**：这是「票首过半才淘汰异议方」的门槛。降到 0.4 会更激进，升到 0.6 更保守。
3. **是否先只重跑 barbarum 看全链路**：其余 7 个项目（含 DNA 组）暂不动。
4. 若上线，`enforce_rank_containment` 里 `blank_species = TRUE` 这条老路径与新的种级闸门在语义上有重叠（前者随科-属级一起置空 Species），当前实测无冲突，但值得在上线后复查一次；`06` 层 `*_agree` 栏的下游消费方尚未验证。

## 六、复现命令

```bash
# 双模式重跑（约 125 秒）
bash /tmp/rc64b_run.sh
# 全量验收（日志 + 台账 + 四方对账 + 四项探针）
bash /tmp/rc64_verify.sh
# 独立复算 *_agree（不依赖 R，从 standardized_*.tsv + 产物重建）
python3 /tmp/chk_agree.py /tmp/rc64b_legacy   /tmp/chk_agree_legacy.txt
python3 /tmp/chk_agree.py /tmp/rc64b_cascade  /tmp/chk_agree_cascade.txt
# 参照表多义性 + 闸门首行/并集口径差异
python3 /tmp/ref_gate_audit.py /home/zhangwenda/database/taxonomy/genus_family_ref.tsv \
        /tmp/rc64b_legacy/species_containment_blanked.tsv /tmp/rc64b_cascade/species_containment_blanked.tsv
# 上一版与当前版的产物字节比对
bash /tmp/cmp_rc64b.sh
```

`cmp_products.py` / `probe_stale.py` / `chk_ph.py` / `probe_sp_gate.py` / `probe_weights.py` / `chk_agree.py` / `ref_gate_audit.py` 均随本报告归档在 `scripts/audit/`。

## 七、从头复核（v6.3 → v6.4 → v6.4b）

复核方式：不看结论，只沿着 19 个改动块逐块追调用方与被调用方，再用三份独立证据压住「口径没被偷偷改掉」这个最难的断言。

### 7.1 `*_agree` 重算的等价性自证

v6.4 把 agree 从「闸门前共识值」换成「最终值」，必须证明变的只是快照，不是口径：

- v6.3 口径 = `match_dt[match_weight > 0]`，其中 `match_dt` 由 `long × cons_long` 合并后 `match_weight := fcase(tolower(Taxon)==tolower(cons_Taxon), rw, default=0)`；v6.4 新块（`:862-868`）用同构的 `long × cons_final` 加同一条字符串相等判据，只有 `cons_*` 的来源从闸门前 `wide` 换成最终 `wide`。
- 关键前提是 `rw` 不为 NA：`rank_weights_dt <- data.table(Rank = names(RANK_DEPTH_WEIGHTS), ...)`（`:769`）覆盖全部 8 个阶元，所以 `match_weight > 0` 与「字符串相等」严格同义，不存在「匹配上了却因为缺权重被漏掉」的角落。
- `n_tot` 仍取自同一个 `long`（`long[, .N, by=.(contig_id, Rank)]`），未改。

### 7.2 独立复算：167136 格逐格重建（`chk_agree.py`）

不调 R、不读 R 的中间产物，直接从 7 个 `standardized_*.tsv`（各含 `contig_id` + 8 阶元）与最终产物重建 `long`、`is_valid_value_vec`、工具 ASCII 排序、`%d/%d: tools` 格式，逐格比对：

```
legacy : 比对 167136 格，一致 167136 格，不一致 0 格
cascade: 比对 167136 格，一致 167136 格，不一致 0 格
格式自检：非 a/b 结构的格子 0 个
```

两种模式、全部 8 阶元、全部 20892 contig，无一格不吻合。这条同时把「等价的旧口径」与「新增的最终值口径」都钉了一遍：被闸门置空的阶元一律 `0/N`，未动的阶元数值与原口径相同。

### 7.3 参照表首行语义：实测无损（`ref_gate_audit.py`）

`load_genus_family_ref_dual()` 里 `unique(ref_dt[...], by = "g_l")` 每个属只留第一行，该语义从 v6.3 原样继承（`norm_tax` → `tax_norm`、`ref_dt` 构造与 `unique` 全部逐行搬移，无改写）。风险在于「一属多科」时首行会掩盖另一行的相容证据，实测：

```
参照表 234804 行 / 去重后属数 NCBI 侧 234526、VMR 侧 2929
[A] 一属多科：NCBI 侧 0 个属（0.000%），VMR 侧 0 个属（0.000%）
[B] 置空清单回查：legacy 1283 行、cascade 58 行
    并集口径会豁免的行 0（NCBI 侧任一科相容 0 / VMR 侧 0）-> 首行口径潜在误伤上限 0.00%
```

一属多科在这个参照表里根本不存在，1283 / 58 行的置空在两个口径下完全一致，闸门动作没有一行是「本可以相容却被首行截断」造成的。

### 7.4 逐点确认的其余断言

| 复核点 | 结论 | 证据 |
|---|---|---|
| `verify_family_genus_consistency` 里 `sp_ck`/`n_sp_conf` 作用域 | 安全 | 赋值与使用同处 `if (!is.null(ref_dual)) {}` 块（`:466` / `:486`），R 函数级作用域，无跨分支裸用 |
| `completeness`/`confidence` 是否有陈旧读者 | 无 | 全文件仅 `:857-858` 赋值、`:1298` 输出列、`:1313` 绘图读取 |
| `final_result` 在 `build_consensus` 之后是否又被改 | 未改 | `:1287` 赋值后只被自检、版本戳、`fwrite` 读取 |
| `analyze_consistency_optimized` 的 `agreement_stats.tsv` | 不受影响 | 其 `Consensus` 由工具自身众数现算（`:1029`），是独立诊断量，不读产物 `*_agree` |
| 占位清理扩到 8 阶元是否误伤 | 无 | 子串白名单不含 `virus`/`other`/`na` 等泛词，8 阶元残留全 0 |
| `opt$output` 是目录 | 是 | 三处落盘均 `file.path(opt$output, ...)`，实测三份台账同目录产出 |
| 计时标签笔误 `闸门-跳参照相容性` | 已修 | 见 7.5 |

### 7.5 复核改掉的两处（纯注释 / 日志串）

1. `enforce_rank_containment` 文档注释仍写「`*_agree` 列保留校准前的投票记录，便于回溯」，与 v6.4 行为直接矛盾，改为「自 v6.4 起对最终值重算，需回溯单工具原始值直接读 `standardized_*.tsv`」。
2. 计时标签 `闸门-跳参照相容性` → `闸门-科属相容性`，`闸门-种科相容性` → `闸门-科种相容性`（与日志正文 `(科-属)`/`(科-种)` 统一）。

改动后重跑（`rc64b_run.sh`，legacy 61s / cascade 62s，均 exit 0）并与上一版产物逐字节比对（`cmp_rc64b.sh`）：

```
legacy / cascade 两模式：22 个数据产物全部 SAME（含 final_integrated_classification.tsv、
  tool_weights.tsv、species_containment_blanked.tsv、taxonomy_gate_check.tsv、7 份 standardized_*、
  agreement_stats.tsv、comparison_* 等）
DIFF 仅 3 处：analysis_summary.txt（内含时间戳）、taxonomy_gate_stamp.tsv（脚本 md5 + 时间）、
  运行日志（时间戳 + 上述两条标签）
```

顺带拿到了一个附加结论：**同输入两次独立运行数据产物逐字节相同**，R 端共识链路是确定性的（无依赖随机或哈希遍历顺序）。

### 7.6 复核后仍然保留的差异项（知情保留，非遗漏）

| 项 | 现状 | 判断 |
|---|---|---|
| `primary_tool` 仍由闸门前共识打分 | 未动（`:768-773`） | 语义是「谁最贴合投票共识」，属溯源标签；若要它反映最终值需另拍板，会改动下游可见列 |
| `rate` 取自闸门前朴素多数 | 未动 | 避免「用最终共识评估工具」的循环，作为先验量反而更稳 |
| `RATE_FLOOR` 对 legacy 同样生效 | legacy 不再是 v6.2 纯复刻 | 实测侧证：科-属校准从「替换 959 / 置空 1126」变到「替换 974 / 置空 1130」，占位清理 Genus 由 1 格变 24 格。若需要 v6.2 原样复刻做基线，必须用旧补丁 |
| `CASCADE_MIN_SHARE` 无上下界钳制 | `RATE_FLOOR` 有 0..1 钳制，它没有 | 环境变量填 >1 会静默退化成逐阶元独立投票；上线前建议一并钳到 (0,1] |
| `load_genus_family_ref_dual()` 被调 4 次（`:465/:506/:622` + 科-属闸门内） | 每次重读 234804 行 | 单次约 2.5–3.0s，占 `build_consensus` 约 1/5；8 项目铺开影响可忽略，可 memo 但不急 |
| `enforce_rank_containment` 的 `blank_species = TRUE` 与种级闸门重叠 | 实测无冲突 | 前者随科-属级一起置空，后者独立判据，置空集合不互相放大 |

## 八、逻辑层验证（单函数人造数据穷举 + 真实数据量化）

前面七节的证据都是「同一份真实输入，换实现、比产物」。这一节换刀口：不看真实数据，把规则拆到单函数层，用人造数据把边界逐一钉住；再用探针量一下新发现的口子在真实数据上到底有多大。

### 8.1 做法

`scripts/audit/logic_selftest.R`（355 行，md5 `1bc0363fe5a8d32c89fb5bd2ac08a724`）：

- 不跑真实数据。只从补丁里 `parse` 出全部 `function` 定义与 12 个常量（`TAX_LEVELS` / `TAX_V2_COLS` / `KNOWN_REALMS` / `SUBRANK_SUFFIX` / `PLACEHOLDER_TAXA` / `RANK_DEPTH_WEIGHTS` / `TOOL_BIAS` / `GENUS_FAMILY_REF` / `TAX_GATE_VERSION` / `CONSENSUS_MODE` / `CASCADE_MIN_SHARE` / `RATE_FLOOR`），eval 定义但不执行主流程。
- 期望值全部手写推算，不从代码反推。`chk` 记 PASS/FAIL，`inf` 记观测项（只记录不判定），有 FAIL 则 `quit(status=1)`。
- 六组：A 投票引擎与阶元归一 10 项 / B 科-属与科-种闸门 6 项 / C 权重计算 6 项 / D Species 质量打分 3 项 / E 端到端合成用例 7 项 / F 淘汰规则临界 3 项。
- 合成数据必须自带合规后缀（Kingdom `-virae$`、Phylum `viricota$`、Class `viricetes$`、Order `virales$`、Family `viridae$`），否则会被后缀闸门整格置空，测不到目标分支。

| 文件 | 行数 | md5 |
|---|---|---|
| `scripts/audit/logic_selftest.R`（= 服务器 `/tmp/logic_selftest_lf.R`） | 355 | `1bc0363fe5a8d32c89fb5bd2ac08a724` |
| `scripts/audit/run_selftest.sh` | 10 | `9cdbfdd591e064a18933bdd1f07df8f4` |
| `scripts/audit/probe_ref_ws_and_genus.py` | 128 | `f698b4adb86e5cfc23173ff22351aafe` |
| `scripts/audit/probe_genus_conflict_breakdown.py` | 93 | `937ff48680eb23aab4118273d9463634` |
| 被测补丁 `/tmp/virus_classifier_analysis.R.cascade_v64` | 1344 | `8e51cf21b8ef801a1539e09272dd0cb1` |

### 8.2 结果

服务器实跑 `Rscript logic_selftest.R /tmp/virus_classifier_analysis.R.cascade_v64`：

```
合计 PASS=72 FAIL=0 INFO=11
== Rscript exit = 0 ==
```

输出 230 行（`/tmp/logic_selftest.out`，已随报告归档）。分组落点：A 组 10/10、B 组 6/6（含 3 项 INFO）、C 组 6/6、D 组 3/3（含 1 项 INFO）、E 组满项（含 2 项 INFO）、F 组 6/6。11 项 INFO 里值得记的是：B4/B5 闸门处置行清单、E7 由种提属被科-属校准改回的实例、D3 `Simplexvirus sp.` = 0.3×1.5 = 0.45 的叠加打分。

### 8.3 首轮 9 个 FAIL 的归因：期望值错，代码对

首轮 PASS=55 FAIL=9，全部归因为我的期望值写错，代码行为本身正确。

**C 组 6 个（C1 三格 + C2 + C3 + C5）**：工具名 `A`/`B`/`C`/`ZZZ` 都没登记在 `TOOL_BIAS`，实际走默认 bias 0.8，我却按已登记工具（ACVirus 1.2）手算。观测值恰好是期望值的 1/1.5，六个数 `12.8 / 102.4 / 10.24 / 0.8 / 51.2 / 5.12` 全部自洽。改期望为 0.8 后全 PASS，并补 C6 用真名反向钉住「查表没失效」：ACVirus 1.2 / VITAP 1.1 / vcontact3 0.7 / contigtax 0.6 → 76.8 / 70.4 / 44.8 / 38.4。

**E5 3 个**：我把 A 的 Phylum 写成 `Peploviricota`、B/C 写 `Uroviricota` 而权重同为 1，A 在 Phylum 层 1:2 就已经出局，根本走不到我设计的 Order 层权重对比（日志「淘汰工具投票权 1 格次」正是 A@Phylum）。重设计为：粗阶元三家报同一值，淘汰精确发生在 Order（A 3 : B 1 : C 1），B/C 在 Family 的权重故意给 5（A 只有 1）。

这轮 9/9 都归因到期望侧，说明自测测的是「文档与代码的差」，不是「代码与代码的差」。

### 8.4 F 组：把「淘汰优先于权重」钉死

E5 在端到端层给出 `Herpesviridae`，但端到端会被闸门与校准加工。F 组绕过闸门直测 `build_cascade_winners`：

- **F1** 淘汰发生在 Order（3:1:1），B/C 在 Family 的权重故意给 5（A 只有 1）：Family 落 `Herpesviridae`，证明被淘汰者不复活、且压过 5 倍权重；粗阶元 0.6 定音。
- **F2** 票首 share 恰 0.5 不淘汰（门槛是严格大于），细阶元回到权重多数（`Mimiviridae`）；粗阶元按行序先到（`Herpesvirales`）。
- **F3** 1:2 时少数 A 被淘汰（`Betaovirales`），幸存工具 C 在细阶元照常定音（`Gammaviridae`）。

同一合成输入在 legacy 侧的 E5 对照给出权重驱动的 `Mimiviridae / Mimivirus / Mimivirus beta`，两模式的差异被压成一对可复现断言。

### 8.5 真实数据量化一：归一化口径不对称（潜伏风险，当前零触发）

背景：`load_genus_family_ref_dual()` 里 NCBI 侧 `tolower(tax_norm(ref$NCBI_Family))` 保留内部空白，VMR 侧 `vmr_pad` 与行侧 `fam_c = gsub("[[:space:]]+", "", ...)` 都去内部空白。若参照表出现带空白（或分号分隔多科）的科名，两侧容忍度不一致会造成误判置空。

`probe_ref_ws_and_genus.py` 实测（参照表 234804 行，列 `Genus,NCBI_Family,NCBI_n,VMR_Family,VMR_n,Domain`）：

```
NCBI_Family 含内部空白 0 行
VMR_Family  含内部空白 0 行
VMR_Family  含 ';'（一属多科）0 行
结论: 不对称当前无实际触发，属潜伏风险
```

即：这个不对称今天是零触发的，不构成正确性问题，可作为清理项排期，不必为它单独重跑。

### 8.6 真实数据量化二：最终 Genus 与 Species 首词（独立属主张）不一致

双名法 Species 的首词本身构成一条独立属主张。最终 Genus 与它不一致，说明 `harmonize_genus_species`（由种提属）先跑、`harmonize_family_genus`（科-属校准）后跑把 Genus 改回去了，而两道闸门都不复核这层关系。

`probe_ref_ws_and_genus.py` 三表实测：

| 产物 | 总行 | 双名法 Species 行 | 属主张 ≠ 最终 Genus | 占比 | 其中 Species_agree 仍写 1/1 |
|---|---|---|---|---|---|
| live v6.2 产物 | 20892 | 14764 | 1246 | 8.44% | 107 |
| v6.4 legacy | 20892 | 13528 | 301 | 2.23% | 60 |
| v6.4 cascade | 20892 | 14138 | 97 | 0.69% | 25 |

`probe_genus_conflict_breakdown.py` 按「首词属的定型科 vs 行内 Family」分层：

| 产物 | 冲突行 | A 同科不同属 | B 跨科 | C 行 Family 空 | D 首词不在参照表 | 带占位限定词 | S_agree=1/1 |
|---|---|---|---|---|---|---|---|
| live v6.2 | 1246 | 532 (42.7%) | 497 (39.9%) | 27 (2.2%) | 190 (15.2%) | 64 | 107 |
| v6.4 legacy | 301 | 16 (5.3%) | 99 (32.9%) | 0 | 186 (61.8%) | 17 | 60 |
| v6.4 cascade | 97 | 16 (16.5%) | **1 (1.0%)** | 0 | 80 (82.5%) | 24 | 25 |

读法：

1. cascade 把真跨科从 497（live）/ 99（legacy）压到 1 行，真实数据侧再次印证「从粗到细逐级淘汰」确实清掉了科层面的异议工具。残留那 1 行是 `SRR23107136_clean_NODE_29173`：F=`Orpheoviridae` / G=`Alphaorpheovirus` / S=`Alphahydrivirus permafrostis`。
2. 「同科不同属」两模式都是 16 行、完全相同，不是投票差异造成，而是科-属校准把 Genus 换成「自洽工具集票首属」的结果，与 Species 首词本就出自两个独立来源。
3. 残留主体（80/97 = 82.5%）是 D 类：首词属在 234804 行的参照表里查不到定型科，闸门按 fail-safe 不判定。里面既有真新属（`Llyrvirus`、`Thornevirus`、`Tokyovirus`），也有工具报的种名首词根本不是属名（`Virus NIOZ-UU157`）。
4. 24 行的 Species 里含 `sp.`/`cf.`/`strain` 等占位限定词，名称本身不完整，首词的属主张天然不可靠（这类 Species 在打分里已被 0.45 系数降权，但没有被排除）。
5. 25 行 Species_agree 仍写 1/1：这一栏的口径是「最终 Species 值与工具上报值的一致性」（见 7.1），不是「Species 与 Genus 的一致性」。它写 1/1 与上面这条失配并不矛盾，但读表的人容易误读。属标签语义问题，不是数值错误。

### 8.7 第八节新增的待拍板项

1. 是否在 `harmonize_family_genus` 之后再加一道「Species 首词属 == 最终 Genus」的复核。可治 16 行同科不同属 + 1 行跨科；D 类 80 行参照表查不到、无法判。修法二选一：把 Genus 拉回首词属，或把这类 Species 降级为单名/置空。
2. `Species_agree` 是否补一个「Species 与 Genus 一致」的并行列，或至少补口径说明（避免 8.6 第 5 条那种误读）。
3. 参照表覆盖：D 类占 82.5% 说明 `genus_family_ref.tsv` 对新属覆盖不足，这是闸门判定的天花板，是否补新版。

### 8.8 第八节复现命令

```bash
bash /tmp/run_selftest.sh                                   # Rscript 自测（PASS/FAIL + 退出码）
python3 /tmp/probe_ref_ws_and_genus.py                      # 参照表空白 + 三表属主张失配
python3 /tmp/probe_genus_conflict_breakdown.py \
        /tmp/rc64b_cascade/final_integrated_classification.tsv \
        /tmp/rc64b_legacy/final_integrated_classification.tsv \
        ~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv
```

一个自查记录：探针三首版的 `QUAL` 正则写成 `\b(sp\.|…)\b`，句点是非单词字符，后面接 `\b` 永不匹配，占位限定词统计恒为 0；已改为 `(?:^|\s)(…)(?=\s|$)`，上表的 24 / 17 / 64 是修正后的数字。

## 九、cascade 之后为什么还有科属打架（闸门关闭实测）

方法：从 v6.4 补丁生成闸门关闭变体 `/tmp/virus_classifier_analysis.R.pregate`（74341 B，md5 `8a4b7a55b5bc7204a8103e416fc01e93`，只注释 `harmonize_family_genus` / `enforce_rank_containment` / `enforce_species_containment` 三处调用，生成器 `scripts/audit/make_pregate_variant.py`）。cascade 模式重跑 55 s，计票与正式版逐字相同（144229 格次计票 / 13681 次工具淘汰事件），产物 `/tmp/pregate_cascade/`，日志 `/tmp/pregate_cascade.log`。

### 9.1 机制：逐级比较比的是「同一阶元的各家取值」，不是「相邻阶元配不配套」

cascade 逐级独立定音：每一级取票首，并把这一级报了别家取值的工具踢出局。它从头到尾不检查相邻阶元是否配套，因此科可以由 A 家工具定、属由 B 家工具定，行内组合不成立。

20892 行的裸共识结果（闸门介入前）：

| 状态 | 行数 |
|---|---|
| 自洽（存在单一工具同时报出这对科属） | 13475 |
| 有科无属 | 3838 |
| 科属全空 | 2300 |
| 有属无科 | 1080 |
| **不配套** | **199** |

199 行的三种来路：

| 成因 | 行数 | 说明 |
|---|---|---|
| 属无任何工具支持 | 85 | 属是「由种提属」从 Species 首词提出、或兜底填出的，与定科的工具无关 |
| 沉默逃逸 | 81 | 报属的工具在科级没有任何值，不构成科级异议因此没被淘汰，留着资格在属级成为票首 |
| 各阶元赢家不同工具 | 33 | 报属的工具在科级报了别的科，它自己的科级主张早被淘汰 |

样例两条：

- `CRR1126135_clean_NODE_29_length_3464_cov_4.068711`：科 `Kyanoviridae` 由 VITAP/genomad 定，属 `Llyrvirus` 由 ACVirus 定，而 ACVirus 科级为空（沉默逃逸）。
- `CRR1440126_clean_NODE_914_length_1792_cov_3.429319`：科 `Marseilleviridae` 由 CAT/mmseqs 定，属 `Fadolivirus` 由 ACVirus 定，而 ACVirus 报的科是 `Mimiviridae`（各阶元赢家不同工具）。

### 9.2 参照侧裁决力与科分布（199 行）

| 参照侧结论 | 行数 |
|---|---|
| 属不在参照表（参照不可判） | 135 |
| 双参照一致反驳行内科 | 38 |
| 双参照不排斥行内科（该属一类多科） | 15 |
| 参照单侧且有反驳 | 6 |
| 参照单侧且相容 | 5 |

共识科分布 top10：Mimiviridae 43、Phycodnaviridae 30、Marseilleviridae 29、Retroviridae 10、Kyanoviridae 8、Mamonoviridae 6、Herelleviridae 6、Pithoviridae 6、Partitiviridae 5、Autographiviridae 5。主角是巨型病毒，原因两重：工具在这类群的科级命中普遍少于属级（沉默逃逸的土壤），参照表又缺这些属（`Llyrvirus`、`Tokyovirus`、`Pandoravirus` 在 234804 行的参照表里查不到）。

### 9.3 闸门动作与「199」的对账

闸门开关两版逐行相比（`scripts/audit/diag_gate_effect_bridge.py`）：Family 变化 **0 行**（v6.4 不含以属修科，科从不被改写）；Genus 被替换 133、被补上 6（合计 = 日志「替换 139」）、被置空 66（日志 67，差 1 行本已为空）；仅 Species 变化 21 行。

对账：199（科属都非空且工具层不自洽）+ 6（属空缺被候选补上）= 205 行净变化，再加 1 行空到空 = 日志的 206。

### 9.4 三个改法在这批数据上的真实覆盖面

| 改法 | 可作用行数 | 实质 |
|---|---|---|
| 清空越权属（现行） | 205 | 管得动，但 135 行属不在参照表，等于「说不清的一律清掉」 |
| 以属修科 | 38 | 只有双参照一致反驳的 38 行改得动，其余属查不到，无从改 |
| 按证据强弱分派 | 38 + 161 | 双参照一致时以属修科，其余走替换/置空 |

### 9.5 复现命令

```bash
python3 /tmp/make_pregate_variant.py                                        # 生成闸门关闭变体
MMPV_CONSENSUS_MODE=cascade Rscript /tmp/virus_classifier_analysis.R.pregate \
    --combined <Votus_combined_taxonomy.tsv> --output /tmp/pregate_cascade    # 55 s
python3 /tmp/diag_rank_winner_mismatch.py                                    # 199 行成因 + 参照侧裁决力
python3 /tmp/diag_gate_effect_bridge.py                                      # 闸门动作与 199 的对账
```

一个坑记录：参照表 `genus_family_ref.tsv` 的 Genus 列是小写（`mimivirus`），探针首版按 ICTV 大写拼写做字典键，199 行全部落到「属不在参照表」，把参照侧裁决力误判为 0。R 层 `load_genus_family_ref_dual` 用 `tolower` 归一，闸门本身没有这个问题。凡离线复核该表，键必须先 `lower()`。

## 十、逐阶元投票实录：科属为什么会由两家工具分别赢下

把 cascade 的投票引擎用 Python 复刻（`scripts/audit/trace_cascade_engine.py`），逐阶元打印「活跃工具的取值与权重 → 票首占比 → 是否淘汰谁 → 落定值」。

验证口径：复刻结果与 `/tmp/pregate_cascade/final_integrated_classification.tsv` 逐格比对，Realm / Kingdom / Phylum / Class / Order / Family **全部 0 差异**；Genus 123 行差异、Species 863 行差异，来自复刻时未跟的两步（`harmonize_genus_species` 由种提属、Species 的 `species_quality_score_vec` 加权），不影响 Family 层的结论。

### 10.1 样例一 `CRR1440126_clean_NODE_914`（7 工具原始取值）

| 阶元 | ACVirus | CAT | VITAP | diamond_lca | genomad | metabuli | mmseqs |
|---|---|---|---|---|---|---|---|
| Realm | Varidnaviria | Varidnaviria | Duplodnaviria | Varidnaviria | Varidnaviria | | Varidnaviria |
| Kingdom | Bamfordvirae | Bamfordvirae | Heunggongvirae | Bamfordvirae | Bamfordvirae | | Bamfordvirae |
| Phylum | Nucleocytoviricota | Nucleocytoviricota | Uroviricota | Nucleocytoviricota | Nucleocytoviricota | | Nucleocytoviricota |
| Class | Megaviricetes | Megaviricetes | Caudoviricetes | Megaviricetes | Megaviricetes | | Megaviricetes |
| Order | Imitervirales | Pimascovirales | | | Algavirales | | Pimascovirales |
| Family | Mimiviridae | Marseilleviridae | | | Phycodnaviridae | | Marseilleviridae |
| Genus | Fadolivirus | | | | | | |
| Species | Fadolivirus algeromassiliense | Marseillevirus LCMAC201 | | | | | |

逐阶元投票：

| 阶元 | 取值与权重 | 票首占比 | 决定 |
|---|---|---|---|
| Realm | Varidnaviria[1.1+0.8+0.7+0.7+0.7] > Duplodnaviria[1.1] | 83.3% | 采用，淘汰 VITAP |
| Kingdom / Phylum / Class | 一致 | 100% | 采用，无异议 |
| Order | Pimascovirales[10.1+8.2] > Imitervirales[12.5] > Algavirales[9.6] | 45.3% | 采用 Pimascovirales，**不淘汰任何人** |
| Family | Marseilleviridae[15.2+14.8] > Mimiviridae[22.3] > Phycodnaviridae[16.5] | 43.6% | 采用 Marseilleviridae，**不淘汰任何人** |
| Genus | Fadolivirus[50.9]（只有 ACVirus 报了值） | 100% | 采用 Fadolivirus |
| Species | Fadolivirus algeromassiliense[53.3] > Marseillevirus LCMAC201[15.2] | 62.2% | 采用前者，淘汰 CAT |

读法：ACVirus 自己拿出来的是 Imitervirales / Mimiviridae / Fadolivirus 一条自洽谱系。它在目级（12.5 对 18.3）和科级（22.3 对 30.0）都输了，但两级票首占比都 ≤ 50%，按规则不淘汰任何人，ACVirus 依然活跃。到属级只有它一个报了值，而属级权重 50.9 远高于其余工具（15.7 至 25.7），属直接归它。于是行内落成「科 = Marseilleviridae（另一家定的）+ 属 = Fadolivirus（ACVirus 定的）」。

另一个细节：这一行连目级也不配套。Order 判给了 Pimascovirales，Family = Marseilleviridae 属于 Pimascovirales（自洽），而 Genus = Fadolivirus 属于 Mimiviridae / Imitervirales（不自洽）。八级各投各的票，跨级组合无人检查。

### 10.2 样例二 `CRR1126135_clean_NODE_29`（沉默逃逸）

七个工具分成两派：Caudoviricetes 派（ACVirus / VITAP / genomad）与 Megaviricetes 派（CAT / diamond_lca / mmseqs）。Realm 级 Duplodnaviria[1.1+0.7+0.7] 以 54.8% 过半，直接把 CAT、diamond_lca、mmseqs 淘汰，此后 Kingdom 到 Order 全部 100% 通过。

关键在 Family 与 Genus：Family 级 Kyanoviridae 由 VITAP(8.8) + genomad(16.5) 投出（ACVirus 在科级**没有值**），100% 通过；Genus 级只有 ACVirus 报 Llyrvirus（50.9）、VITAP 报 Acionnavirus（17.6），Llyrvirus 以 74.3% 胜出。ACVirus 从头到尾没在科级说过话，不构成「科级异议」，不会被淘汰，留着资格到属级把属赢了。

### 10.3 样例三 `CRR1126135_clean_NODE_499`（属在投票之后被改写）

Family 级 Mimiviridae 五个工具 100% 通过；Genus 级只有 ACVirus 报 Fadolivirus（100%）。此时行内科属是配套的（Fadolivirus 属于 Mimiviridae）。但 Species 级 metabuli 的 `Mimivirus LCMiAC01`（权重 85.6）以 61.6% 压掉 ACVirus 的 `Fadolivirus algeromassiliense`，随后 `harmonize_genus_species`（由种提属）把 Genus 改成 Species 首词 `Mimivirus`，已淘汰工具的属主张就此作废。

结论：这类行的科属不配套**不是投票引擎造成的**，而是投票之后由种提属把属单独改了一次，而这一步不回头检查科。

### 10.4 199 行按「属是谁带来的」重分解

| 属的来源 | 行数 |
|---|---|
| cascade 属级投票选出来的 | 111 |
| cascade 之后由种提属从 Species 首词提上来的 | 88 |

再按科级投票状态切（全 199 行）：

| 科级这一关的情况 | 行数 |
|---|---|
| 科级票首过半，已把在科级报别家答案的工具淘汰 | 155 |
| 科级票首不过半，按规则不淘汰任何人 | 44 |

155 行里，属级赢家仍有 87 行在科级是沉默的（沉默不算反调，因此没被淘汰）。

### 10.5 机制总结

cascade 淘汰的是「在某个阶元唱反调的工具」，不是「谱系不配套的工具」。级与级之间只传递一件事：你在这一级唱了反调，更细的级就别投了。于是：

1. 票首没过半的级不淘汰任何人，输家照样参加下级投票，靠该级权重翻盘（44 行，属级权重 ACVirus 50.9 对 CAT 15.7 / VITAP 17.6 / genomad 14.4 / mmseqs 16.0）；
2. 某级沉默的工具永远不被淘汰，可以在它没表态的级被别人定音、又在更细的级成为唯一发言者（87 行）；
3. 属级落定之后，由种提属还能把属单独换掉而不回头查科（88 行）。

复现：`python3 /tmp/trace_cascade_engine.py`、`python3 /tmp/show_contig_votes.py <contig_id>`。

## 十一、v6.5：逐级投票改为「只看出现次数」（A/B 实测）

### 11.1 补丁改动

补丁：`/tmp/virus_classifier_analysis.R.cascade_v65`（1360 行，md5 `bda4ffd41dbb15af45929c705f8b6c9f`）。开关：`MMPV_VOTE_WEIGHT`，取值 `count`（默认）/ `weighted`。

| 位置 | 改动 |
|---|---|
| `:56` | 新增 `VOTE_WEIGHT_MODE`（env `MMPV_VOTE_WEIGHT`，默认 `count`） |
| `:295` | `compute_tool_weights` 阶元循环顶部：count 模式 → 每工具每阶元权重恒为 1，直接 `next` |
| `:725` | `build_consensus` 种名质量降权只在 `weighted` 模式生效 |
| `:1300` | 启动日志打印 `投票计权模式: count` |
| `:75` | `TAX_GATE_VERSION` 6.4 → 6.5 |

`RANK_DEPTH_WEIGHTS` 在阶元内是常数，对占比无影响，count 模式下彻底不用；`TOOL_BIAS`、`consensus_stats` 的 rate 校准、`species_quality_score_vec` 三者在 count 模式下全部停用。`build_cascade_winners` 本身不读该开关（只对 Weight 列求和），因此 `logic_selftest.R` 的 72 条断言不受影响。

### 11.2 回归验证：weighted 口径零扰动

同一补丁加 `MMPV_VOTE_WEIGHT=weighted` 重跑，与 v6.4 正式产物**逐字节相同**：

```
9ab1f1d7680955acdf82c4200c6efef6  /tmp/rc65_weighted/final_integrated_classification.tsv
9ab1f1d7680955acdf82c4200c6efef6  /tmp/rc64b_cascade/final_integrated_classification.tsv
d0e1bd8d69322799da9391a345a8d4ac  tool_weights.tsv （两边同）
```

差异全部来自计权模式本身。

### 11.3 实测对照（同一输入 20892 contig）

| 指标 | v6.4 加权 | v6.5 只数次数 | 变化 |
|---|---|---|---|
| 阶元计票格次 | 144229 | 144771 | 持平 |
| 工具淘汰事件（contig × 工具 对，补丁日志 `n_prune` 口径） | 13681 | 8523 | -38% |
| 其中「过半且有异议」的计票格次（只量了 count 侧，见 §13） | 未测 | 7519 | |
| 票首平票（多值同分）格次 | 0 | 9789 | 新出现 |
| 科-属校准 替换 / 置空 | 139 / 67 | 296 / 192 | 约 2 倍 |
| 科-种闸门置空 Species | 58 | 87 | 1.5 倍 |
| 裸共识「科属不同源」行 | 199 | 473 | 2.4 倍 |
| 成品 Family 有值行 | 17512 | 17710 | +198 |
| 成品 Genus 有值行 | 14694 | 14607 | -87 |
| 成品 vs 加权版逐行差异 | | Family 366 / Genus 466 / Species 1334 行 | |

> 口径更正：本表原先把 13681 / 8523 记作「淘汰异议工具格次」是错的。补丁里 `n_prune <- n_prune + nrow(kill)`，`kill` 是 (contig_id, Tool) 对，所以这个数是**工具淘汰事件数**，不是格次数。原先与之并列的 132826 / 137252 来自另一个早期探针，粒度和这边不同，已删掉不用。端正的格次级分类见 §13。

成品自检两边均 PASS（双参照一致冲突 0 / 物种双参照冲突 0）。

### 11.4 为什么等权会让淘汰机制失效

7 家工具等权时票首占比只能是 1/7 的整数倍；两派 3:3、2:2 对撞恰好卡在 50%，按「> 50% 才淘汰」的规则不淘汰任何人。图景变成：每一级都是独立计票，cascade 退回 legacy 行为。

样例 `CRR1126135_clean_NODE_29`：

| 模式 | Realm 计票 | 是否淘汰 | 最终 Family |
|---|---|---|---|
| v6.4 加权 | Duplodnaviria 54.8% | 淘汰 CAT、diamond_lca、mmseqs | Kyanoviridae |
| v6.5 计数 | Duplodnaviria 3 : Varidnaviria 3 = 50.0% | 不淘汰 | Marseilleviridae（Family 层 2:2 平票，按行序取 CAT 的值） |

新增的平票还带来一个以前从不生效的裁决点：票首并列时取「`names(data_list)` 里最早的工具所报的值」（ACVirus → CAT → VITAP → diamond_lca → genomad → metabuli → mmseqs）。加权模式下权重是连续数、几乎不会平票，这套行序裁决几乎从未生效；等权后它决定了 9789 个格次（占计票格次 6.8%）。

### 11.5 待定

1. 门槛是否维持 `> 50%`：等权后「过半且有异议」的计票格次只有 7519（5.2%），工具淘汰事件 8523 次，实际上把淘汰关了。格次级实测与四种放宽方案的代价见 §13。
2. 平票裁决规则：继续用表格行序（ACVirus 等前置工具系统性获益），还是改为「平票则本阶元不淘汰、且交由上一级已定值的谱系约束」。
3. `species_quality_score_vec` 是否恢复：本次按「不看权重」一并停用，Species 变化 1334 行里含这部分影响，需要时可用单开关单独切出来量。

## 十二、v6.6：平票改成两级裁决（下一级一致性优先，工具顺序兜底）

### 12.1 规则（2026-08-30 定）

平票即同一阶元上多个取值权重相等。v6.5 的裁决是「取表格行序最早的工具所报的值」，行序一改结论就变。v6.6 把它改成显式两级：

1. **第一级：看下一（更细）级的谐调度**。平票各候选的支持工具里，谁在下一级上取值更齐整（同值工具数 / 支持工具数），谁胜出。
2. **第二级：工具顺序兜底**。仍分不出的（Species 本无下一级、两候选下一级都沉默、两候选谐调度一样），按 `TIE_BREAK_ORDER` 取先到者：`ACVirus > CAT > VITAP > diamond_lca > genomad > metabuli > mmseqs`（未列入的工具按输入原序追在末尾）。

补丁：`/tmp/virus_classifier_analysis.R.cascade_v66`（1414 行，md5 `dfb1cd4fbf5c5de34b186c04cd0f9b3e`）。v65 保留为「纯计数 + 隐式行序」基线，两份并存便于 A/B。

改动点：

| 位置 | 改动 |
|---|---|
| `:42-50` | 新增常量 `TIE_BREAK_ORDER` + 工具分序函数 `order_tools_by_tiebreak()` |
| `build_cascade_winners` | `order_pos`（long 行序）→ `tool_rank = match(Tool, tool_order)`；`agg` 新增 `is_top` / `n_top` / `pick`；平票时先算下一级谐调度再定 winner |
| `:734` | `build_cascade_winners(long, tool_order)`，`tool_order` 由 `TIE_BREAK_ORDER` 显式排定 |
| legacy 分支 | `votes` 带 `t_rank`，`setorder(..., -total_weight, t_rank)`（与显式规则一致，当前两者结果相同） |
| `harmonize_family_genus` | 属级候选同分同样按 `TIE_BREAK_ORDER` 显式排定 |
| `:85` | `TAX_GATE_VERSION` 6.5 → 6.6，版本注补 6.5 / 6.6 两行 |

### 12.2 weighted 口径零扰动（硬回归）

同一补丁加 `MMPV_VOTE_WEIGHT=weighted` 重跑，与 v6.4 正式产物逐字节相同（改前改后一样）：

```
9ab1f1d7680955acdf82c4200c6efef6  /tmp/rc66_weighted/final_integrated_classification.tsv
9ab1f1d7680955acdf82c4200c6efef6  /tmp/rc64b_cascade/final_integrated_classification.tsv
d0e1bd8d69322799da9391a345a8d4ac  /tmp/rc66_weighted/tool_weights.tsv
d0e1bd8d69322799da9391a345a8d4ac  /tmp/rc64b_cascade/tool_weights.tsv
```

原因：加权模式权重连续，平票恒 0；而工具分序在 7 家工具上恰好与输入列序一致，所以 `order_pos` 换成显式序也不动结果。这条同时证明「工具顺序」这个变量已经被隔离干净。

### 12.3 平票去向实测（`/tmp/rc66_count`，62 秒）

```
逐级淘汰: 阶元计票 144771 格次, 淘汰工具投票权 8523 格次 (票首门槛 0.50);
         平票 9789 格次: 下一级一致性定 1879, 工具顺序兜底 7910
科-属校准(自洽工具集): 替换 311 行, 置空 141 行
科-种相容性约束(科-种): 双参照一致反驳 85 行已置空 Species
科-属自检通过: 双参照一致冲突 0 条; 单侧冲突 9 条（设计内不动）
```

平票 19.2% 由下一级一致性分辨，80.8% 落到工具顺序。按阶元看平票分布：Species 2443 > Family 1448 > Order 1148 > Class 1084 > Genus 1027 > Phylum 958 > Kingdom 854 > Realm 827。属以上的平票里，下一级一致性是有效信息；Species 的 2443 格全走工具顺序，因为种级下面没有级了。

### 12.4 v6.5 隐式行序 → v6.6 显式两级裁决：成品逐格差异

| 阶元 | 差异行 |
|---|---|
| Realm | 3 |
| Kingdom | 8 |
| Phylum | 12 |
| Class | 292 |
| Order | 71 |
| Family | 59 |
| Genus | 62 |
| Species | 16 |

粗阶元的差异沿链条向下传递（Class 一级改了 292 行，下级再各自收敛，Family 只剩 59 行不同）。

### 12.5 独立复刻核验（`trace_cascade_engine.py`）

Python 仿真同步换成两级规则后，对 `/tmp/pregate66_count` 逐格比对：

```
Realm 0 / Kingdom 0 / Phylum 0 / Class 0 / Order 0 / Family 0
Genus 285 / Species 202（投要之外的后置闸门，仿真未跟）
```

两份实现（R 与 Python）在 144771 个计票格次上逐格相同，包括 9789 个平票的两级裁决与归属（Python 侧同样得出 1879 / 7910）。

跨版本对比（同一诊断脚本、同一 pregate 口径）：裸共识「科属不同源」v6.4 199 行 → v6.5 473 行 → v6.6 **437 行**。

### 12.6 待定（更新）

1. ~~平票裁决规则~~：已定，就是上面两级。11.5 的第 2 条关掉。
2. 门槛 `> 50%` 是否维持：等权后「过半且有异议」的计票格次只有 7519（5.2%），工具淘汰事件 8523 次，实际上把淘汰关了。已量出「唯一票首但不过半」719 格次与四种放宽方案的代价，见 §13，仍未拍板。
3. cascade 是否上线 8 个项目、是否先只重跑 barbarum。

### 12.7 规则层自测：78 条断言全过，且能区分新旧规则

自测脚本 `scripts/audit/logic_selftest.R`（服务器 `/tmp/logic_selftest_lf.R`，md5 `0e26a0165969ebc195e5928f5ce6d313`，输出 `/tmp/selftest_v66.out`）：

```
[环境] VOTE_WEIGHT_MODE=count | TIE_BREAK_ORDER=ACVirus>CAT>VITAP>diamond_lca>genomad>metabuli>mmseqs | TAX_GATE_VERSION=6.6
合计 PASS=78 FAIL=0 INFO=11
```

**关键：旧基线断言在 v66 下确实会失效，说明这组测试真能盖住规则改动。** 旧基线输出 `/tmp/logic_selftest.out`（那次默认补丁是 v64，`PASS=72 FAIL=0`；v65 的平票同样按行序，这三条行为不变）里的三条正是旧行序规则：

| 断言（旧版） | 旧基线结果 | v66 结果 | 现在写成 |
|---|---|---|---|
| `A1 Kingdom 1/2=0.5 不淘汰取行序先到` | `Orthornavirae` PASS | `Shotokuvirae` | A 无 Phylum、C 有 Uroviricota -> 下一级一致性定 C |
| `A2 0.5 取行序先到` | `Genusalpha` PASS | `Genusbeta` | B 报了 Species、A 没有 -> 下一级一致性定 B |
| `A5 换行序后翻转为 B` | `Genusbeta` PASS | `Genusalpha` | 新版应「不翻转」，断言改成稳定性 |

本次同步做了三件事：

1. 自测环境补上 `VOTE_WEIGHT_MODE`（`Sys.getenv` 派生常量不在框架的提取范围，以前跑到 C 节会 `object not found`）。`scripts/audit/fix_selftest_env.py`。
2. `win()` 助手按生产调用点接上 `order_tools_by_tiebreak()`，平票用例真测规则而不是测行序。新增 A5b（平票兜底按 TIE_BREAK_ORDER，行序反着给也改不了赢家）、A5c（下一级一致性压倒工具顺序：mmseqs 报了种，赢过更靠前的 ACVirus）、A5d（两级拼谐度相同 -> 回到工具顺序）、A5e（分序函数本身：登记名按次序排，未登记名按原序追尾，正序幂等）。
3. C 节显式切 `weighted`（那 10 条断言本来就是加权口径），另加 C7 钉住新默认：count 下权重恒 1，bias / rate / RATE_FLOOR 全部停用。

至此「先下一级一致性、最后按固定工具顺序」两级规则同时被四道绳拴住：R 实跑产物、Python 独立复刻、规则层单测、weighted 字节级回归。

## 十三、门槛 `> 50%` 到底裁掉了什么（格次级实测）

探针 `scripts/audit/probe_share_threshold.py`（服务器 `/tmp/probe_share_threshold_lf.py`，输出 `/tmp/probe_share_out.txt`），输入 `/tmp/rc66_count`（count 模式成品，contig 21630）。

### 13.1 计票格次四分类

| 类别 | 格次 | 占比 | 含义 |
|---|---|---|---|
| E0 全票一致 | 126744 | 87.55% | 活跃工具取值全相同，无人可踢 |
| T 并列票首 | 9789 | 6.76% | n_top ≥ 2；数学上并列票首占比恒 ≤ 50%，永不触发淘汰，走两级平票裁决 |
| E1 过半且有异议 | 7519 | 5.19% | 现状真正会踢人的格次 |
| U 唯一票首但不过半 | 719 | 0.50% | 票势其实明确，被门槛拦下 |
| 合计 | 144771 | | 与补丁日志「阶元计票 144771 格次」一致 |

把「有异议」的格次加起来（E1 7519 + T 9789 + U 719 = **18027**）才是真正需要裁决的部分，占全部计票格次的 12.45%。**U 占其中的 4.0%。**

交叉验证：本探针在同一遍仿真里累加被淘汰的 (contig, 工具) 对，得 **8523**，与补丁日志 `n_prune` 严格相等，说明活跃集动态与 R 实现一致。

### 13.2 U 类的形状：几乎全是「精确卡线」

| 票型（有值工具数） | 格次 | 占比 | 票首占比 |
|---|---|---|---|
| 2:1:1（4 家） | 574 | 79.8% | 恰好 50% |
| 3:2:1（6 家） | 106 | 14.7% | 恰好 50% |
| 2:1:1:1（5 家） | 25 | 3.5% | 40% |
| 3:1:1:1（6 家） | 14 | 1.9% | 恰好 50% |

719 个 U 格次里，**精确卡在 50% 的有 694 个**，严格小于 50% 的只有 25 个（全是 2:1:1:1）。按阶元：Class 112 / Phylum 107 / Species 103 / Family 102 / Realm 92 / Kingdom 90 / Order 88 / Genus 25。

样例 `CRR1440126_clean_NODE_709`：Family 2:1:1 票首 Hepaciviridae 50%、Genus 2:1:1 票首 Orthohepacivirus 50%、Species 2:1:1 票首 `Orthohepacivirus hominis` 50%，三级都卡在门口。

### 13.3 四种放宽方案的代价（同一引擎只换淘汰判据，全量仿真）

| 方案 | 判据 | 变动格次 | 涉及 contig |
|---|---|---|---|
| A 现状 | `share > 0.5` | 基线 | 0 |
| B 降门槛 | `share > 0.40` | Realm 0 / Kingdom 8 / Phylum 20 / Class 279 / Order 335 / Family 324 / Genus 251 / Species 283 | **487** |
| C 唯一票首即定音 | `n_top == 1` | Realm 0 / Kingdom 1 / Phylum 4 / Class 3 / Order 31 / Family 30 / Genus 81 / Species 66 | **103** |
| D 门槛 0.40 且票首唯一 | `share > 0.40 and n_top == 1` | 同 C 的阶元形状（Genus 78 / Species 63） | **100** |
| E 唯一票首且 `share >= 0.5` | `share >= 0.5 and n_top == 1` | 同 D | **100** |

Genus/Species 含复刻未跟的后置闸门，数字仅作量级参考；Realm 到 Family 可信。

读法：

1. **B 的影响面大得奇怪，原因在于它把并列也当成了过半。** 2:2 = 50% > 0.40，门槛一降，平票裁决就被绕过，并列的败方直接被踢。这与「平票走两级裁决」的定案自相矛盾，不建议。
2. C 只骑在 719 个 U 格次上，最终只动 103 个 contig（样本的 0.48%）。
3. D 与 E 结果完全一致（100 contig），说明 25 个 share 40% 的格次只贡献了 3 个 contig 差异。**只想救回卡线格次，最小改动是 E：把 `>` 改成 `>=` 并加一个 `n_top == 1` 前置。**
4. 三种放宽都不动 Realm，因为 Realm 层的 18096 个格次里 14742 个是全票一致。
5. 表里的「变动格次」是**全量重仿真**的最终值差异（轨迹含回馈）。但读「E 放行 694 个格次」这句话时要注意粒度：694 是在 v6.6 轨迹上数出来的卡线格次，不是 E 上线后的净增格次。真实净增见 §14.3。

### 13.4 附：为什么「3:3 平票」在这份数据里不是主角

等权后数学上会出现 3:3 = 50% 的平票（T 类，9789 格次）。但 U 类（唯一票首卡线）只有 719 格次，两者不是同一回事：前者谁也说不服谁，后者是「2 家报 A、1 家报 B、1 家报 C」，**票势已经明确，只是被「必须过半」的规则挡住**。大王问的就是这一类。

## 十四、v6.7：E 方案落地实跑与代价核算（2026-08-30）

### 14.1 补丁与三跑台账

补丁 `scripts/audit/patch/virus_classifier_analysis.R.cascade_v67`（由 v6.6 复制，5 处改动）：版本注加 6.7 行 + `TAX_GATE_VERSION <- "6.7"`；新增 `n_cell <- 0L`；win 表带出 `win_ntop`；淘汰判据改为 `if (any(win$share >= min_share & win$win_ntop == 1L))` 且 `prune_win <- win[share >= min_share & win_ntop == 1L, ...]`；日志改为「逐级淘汰: 阶元计票 %d 格次, 淘汰格次 %d, 淘汰工具投票权 %d 次 (唯一票首且票首 >= %.2f)」。运行脚本 `scripts/audit/rc67_run.sh`。

| 跑次 | 日志关键行 |
|---|---|
| count（正式口径） | 阶元计票 144619 格次, 淘汰格次 7698, 淘汰工具投票权 8995 次 (唯一票首且票首 >= 0.50); 平票 9637 格次: 下一级一致性定 1851, 工具顺序兜底 7786 |
| weighted（回归） | 阶元计票 144229 格次, 淘汰格次 11466, 淘汰工具投票权 13681 次; 平票 0 格次 |
| pregate（闸门关闭） | 供诊断复刻用，产物 `/tmp/pregate67_count` |

后置闸门联动：科-属校准 替换 296 行 / 置空 88 行（v6.6 为 311 / 141）；科-种 双参照一致反驳置空 Species 51 行（v6.6 为 85）。两项都降低，说明递上来的行更自洽了。

### 14.2 weighted 硬回归：字节级零扰动

| 产物 | v6.7(weighted) | v6.4 正式产物 |
|---|---|---|
| `final_integrated_classification.tsv` | `9ab1f1d7680955acdf82c4200c6efef6` | 同 |
| `tool_weights.tsv` | `d0e1bd8d69322799da9391a345a8d4ac` | 同 |

原因不是巧合：只要计权方式让「唯一票首恰好等于 50%」成为零概率事件，`>=` 与 `>` 就是恒等变换。weighted 口径因此继续充当冻结参照。

### 14.3 净效果不是 +694，而是 +179 格次

静态读法（在 v6.6 轨迹上数 U 格）会说「E 放行 694 个卡线格次」。同一引擎换判据全量重跑（轨迹含回馈）后，真实数字是：

| 口径 | 淘汰格次 | 淘汰工具投票权 |
|---|---|---|
| A `share > 0.5` | 7519 | 8523 |
| E 唯一票首且 `share >= 0.5` | 7698 | 8995 |
| 净增 | **+179** | +472 |

原因是淘汰具有**吸收性**：粗阶元多踢一次，下游格次的活跃集就变了，很多卡线格次在新轨迹里根本不再出现。同理计票格次总量 144771 → 144619（−152），平票 9789 → 9637（−152），都是活跃集收缩的副产品。

交叉验证：探针改成全量重仿真后逐位复现 R 日志（A 7519/8523，E 7698/8995），说明量的是同一件事。

### 14.4 成品逐格差异（count 口径，v6.6 → v6.7）

行数 20892 / 20892，列名一致。非空行数（NA 视作空）：

| 阶元 | 非空 v66 | 非空 v67 | 差 | 有→空 | 空→有 | 值变 |
|---|---|---|---|---|---|---|
| Realm | 19015 | 19015 | 0 | 0 | 0 | 0 |
| Kingdom | 19142 | 19142 | 0 | 0 | 0 | 1 |
| Phylum | 19039 | 19039 | 0 | 0 | 0 | 4 |
| Class | 19457 | 19457 | 0 | 0 | 0 | 3 |
| Order | 18028 | 18007 | **−21** | 21 | 0 | 10 |
| Family | 17710 | 17691 | **−19** | 19 | 0 | 11 |
| Genus | 14658 | 14643 | −15 | 19 | 4 | 8 |
| Species | 17239 | 17230 | −9 | 23 | 14 | 6 |

涉及 contig 57（20892 的 0.27%）。总计：失 82 个赋值 / 得 18 个 / 改值 35 个，净少 64 个赋值。**代价全部落在 Order 及以下，粗阶元只改值不置空。**

典型变化（`/tmp/cmp_v66_v67.out` 有 57 行全清单）：

- `CRR1440136_clean_NODE_104`：v66 是 Riboviria + Nucleocytoviricota + Imitervirales + Mimiviridae/Mimivirus，v67 改为 Riboviria + Pararnavirae + Artverviricota + Ortervirales + Retroviridae/Gammaretrovirus。v66 那条链把 Riboviria（RNA 病毒域）和 Nucleocytoviricota（核质巨 DNA 病毒门）拼在一起，本身矛盾；v67 换来一条自洽的逆转录病毒链。
- `CRR1440127_clean_NODE_1826`：Phylum 由 Uroviricota 改为 Nucleocytoviricota，代价是 Family/Genus/Species 置空。
- `CRR527046_clean_NODE_12`：清掉与 Caudoviricetes 矛盾的 Imitervirales/Mimiviridae，保住 Genus/Species（Casadabanvirus）。

### 14.5 科属自洽性：437 → 369（3.15% → 2.68%）

判据：科、属都非空，且找不到任何单一工具同时报了这两个值（工具层不自洽）。诊断跑在闸门关闭的 `pregate66/67` 上，用 `trace_cascade_engine.py`，`MMPV_RULE=gt/ge`。

**先看成品层（闸门已开启），同一个判据：**

| 口径 | 科属都非空（分母） | 工具层不自洽 |
|---|---|---|
| v6.6 `rc66_count` | 13732 | **0（0.00%）** |
| v6.7 `rc67_count` | 13711 | **0（0.00%）** |

成品层的「科属打架」被闸门清到 0，两版一样。下面的 437/369 是**闸门跑之前**的中间态，它衡量的是「闸门要动多少手术」，不是成品质量。

| 指标 | A `gt` | E `ge` |
|---|---|---|
| 科属都非空行数（分母） | 13858 | 13784 |
| 其中工具层不自洽 | **437（3.15%）** | **369（2.68%）** |
| ├ A 属就是投票结果（不匹配由投票本身造成） | 225 | 173 |
| └ B 属在 cascade 之后被改写（属 == Species 首词） | 212 | 196 |
| 科级票首达阈值且唯一并淘汰了异议工具 | 306 | 304 |
| 科级票首未达阈值或并列 → 不淘汰，属级重新投票 | **131** | **65** |
| 属级赢家在科级报了别的科/部分沉默 | **68** | **15** |

后两行是矛盾的直接产地：科级没淘汰，属级就换了另一批工具来投，于是科属由两家说出。E 把这两个桶砍掉一半到四分之三。换算成闸门手术量：科-属校准置空 141 → 88 行，科-种相容性置空 85 → 51 行，合计少动 87 行。复刻对账：两轮 Realm..Family 与 R 产物 0 差异（Genus 285/202、266/201 的残差是复刻未跟的后置闸门），平票去向 1879+7910 与 1851+7786 也和 R 日志逐位相同。

### 14.6 平票是谁说了算（ACVirus 归属实测）

表头是「赢家值由哪个工具提供」（按票首值归属到兜底序最早的支持工具）：

| 工具 | 顺序定 | 下一级定 | 合计（gt） | gt 占比 | 合计（ge） |
|---|---|---|---|---|---|
| ACVirus | 5365 | 1195 | **6560** | 67.0% | 6427 |
| CAT | 869 | 227 | 1096 | 11.2% | 1089 |
| metabuli | 484 | 212 | 696 | 7.1% | 693 |
| VITAP | 561 | 77 | 638 | 6.5% | 631 |
| diamond_lca | 475 | 67 | 542 | 5.5% | 541 |
| genomad | 156 | 7 | 163 | 1.7% | 163 |
| mmseqs | 0 | 94 | 94 | 1.0% | 93 |
| 合计 | 7910 | 1879 | 9789 | | 9637 |

读法：`TIE_BREAK_ORDER` 把 ACVirus 放在首位，平票格次里有 67% 最终取了 ACVirus 报的值；单看「顺序兜底」那一档，ACVirus 占 5365/7910 = 67.8%。这批格次占全部计票格次（144771）的 3.7%，是当前结果里「顺序即先验」的集中体现。mmseqs 排最后，顺序档 0 次，但在「下一级一致性」档反而拿下 94 次。

### 14.7 结论

E 的收益：科属不自洽率 3.15% → 2.68%（分子 −68），「科级不淘汰导致属级另起炉灶」131 → 65，「属级赢家在科级报别的科」68 → 15，科-属校准的置空 141 → 88 行。代价：净少 64 个赋值（失 82 / 得 18）、35 个格次改值，影响 57 个 contig（0.27%），weighted 口径零扰动。

待大王拍板：E（v6.7）是否上线；若上线，cascade 是八个项目一起换，还是先只重跑 barbarum 看产物流水与论文取数是否要跟着改。

**但这个「收益」的说法需要收紧。** 上一栏的 3.15%/2.68% 是闸门关闭时的中间态，成品层同一个判据两版都是 0 行：闸门把科属打架清干净了，E 只是让闸门动手更少（少置空 87 行）。所以 E 的真实面目是：

1. **口径变化，不是 bug 修复。** 旧判据 `share > 0.5` 在「票首唯一」时等价于「过半」，在「并列」时永不触发；E 把「唯一票首恰好 50%」也放行，于是 2:1:1 的格次里两家反对者被踢。这是从 **majority rule 转到 plurality rule**（相对多数决），要在方法章里辩护「为什么 2 票可以踢掉两个各 1 票的异见」。
2. **净账接近零。** 闸门少置空 87 行（带外部参照证据的纠错）换来自己净少 64 个赋值、35 个格次改值（无外部证据的程序性沉默）。一个换一个。
3. **两版出货都自洽。** 成品层自洽率 100% 对 100%，57 个 contig 的差别不改任何结论。

## 十五、建议：不上线 E，改做三件事（2026-08-30）

### 15.1 建议不上线 `share >= 0.5`

依据就是 §14.7 那三条。为一个成品层等价、净账为零、却要把「多数决」改成「相对多数决」的改动动主线数据，方法章的辩护成本大于收益。判据描述不精确这点（`share > 0.5` 实际是「非并列且过半」）在注释里写清楚即可，不必要改行为。

主线 `05_Taxonomy/Votus.integrated/final_integrated_classification.tsv`（md5 `12e8f1715922ee967db9b10cb70e523a`）不动，八个项目不重跑。

### 15.2 该做的三件事

**一、把「投票定值 / 闸门管自洽」的分工固化。** 现在的代码本来就是这个分工，但注释与文档没说清，于是每轮都在讨论「门槛该不该动」。动作：补丁头加一段架构说明，报告固化。零风险、不改行为、不需要重跑。

**二、量化 `TIE_BREAK_ORDER` 的隐性先验（优先级最高）。** 实测：平票格次 9789，占全部计票格次（144771）的 6.76%；其中顺序兜底 7910，ACVirus 独占 5365，即**全部计票格次的 3.7%** 由一个排序约定决定。门槛争议只涉及 0.5%（719 格次），顺序先验是它的 **7 倍**。

做法：跑两次对照（ACVirus 从首位挪到末位；或按工具在域内的历史一致性排序），看成品变动多少 contig。结论分两支：变动 < 1% → 方法章写一句「平票按固定顺序兜底，顺序敏感性 <1%」；变动大 → 平票改为「弃权留空」。

**三、把「平票留空」摆上桌（等第二件事的结果）。** 平票时该阶元不猜、留空，代价是空值变多（7910 格次，约占计票 5.5%），换来结果不被排序约定绑架。与「宁空不猜」的取向一致，但要先看论文能接受多少空值。

### 15.3 如果仍要上线 E，最低条件

只 barbarum；补丁 `TAX_GATE_VERSION` 与管线 `virome_pipeline.py:51 TAX_GATE_REQUIRED` 同步 6.2 → 6.7（版本戳不达标会强制重跑 05，不需要手工清产物）；产物先落备份目录；57 个 contig 逐条过一遍再决定是否铺开八个项目。

## 十六、定稿（2026-08-30 晚，本节为终局结论）

### 16.1 三个决定

| 项 | 决定 | 依据（实测） |
|---|---|---|
| 共识引擎 | **上线 v6.6 行为**，文件为 `patch/virus_classifier_analysis.R.cascade_v66f` | 与 v6.6 逐字节同产物（`4c31be09bf2f2b035dfb1df0cd08cf4f`），只是把注释写实 |
| 方案 E（`share >= 0.5`） | **不上线** | §14.7 / §15.1：成品层等价、净账为零、性质是多数决改相对多数决 |
| 平票留空（v6.9） | **不上线** | 代价见 16.3；ACVirus 恰是实测质量最高的工具，置空等于丢掉最可信的一票 |

定稿脚本 md5 `c77c43e2468ad3de48466da3e89fff9e`（LF 版 `/tmp/v66f_lf.R`，原补丁 `/tmp/virus_classifier_analysis.R.cascade_v66f`）。

### 16.2 上线影响（vs 线上成品 `12e8f171…`）

| 阶元 | 线上非空 | 定稿非空 | 变动 | 有→空 | 值变 |
|---|---|---|---|---|---|
| Realm | 19015 | 19015 | 69 | 0 | 69 |
| Kingdom | 19142 | 19142 | 112 | 0 | 112 |
| Phylum | 19039 | 19039 | 100 | 0 | 100 |
| Class | 19459 | 19457 | 396 | 2 | 394 |
| Order | 18162 | 18028 | 355 | 134 | 221 |
| Family | 17919 | 17710 | 658 | 209 | 449 |
| Genus | 16022 | 14658 | 2311 | 1364 | 947 |
| Species | 18056 | 17239 | 3434 | 817 | 2617 |

涉及 contig **4354**（20.8%）。粗阶元（Realm/Kingdom/Phylum）只改值不改覆盖，各自变动 0.3% ~ 0.6%；细阶元覆盖下降集中在 Genus（−1364，8.5%）与 Species（−817，4.5%）。

这批减少是**设计意图**：旧口径在无共识时回落到「最完整工具行」，该行由单一工具撑起整条世系；新口径删除该兜底、且被淘汰工具不复活，于是这些「只有一个工具说话」的格子改为留空。丢掉的正是本次改造想要停止信任的那部分赋值。

### 16.3 平票：为什么保留固定顺序，代价是多少

**留空代价（v6.9 实测，`tb_count`）**：顺序兜底档 7910 格次置空后，Realm 掉 822、Kingdom 830、Phylum 932、Class 558、Order 973、Family 635、Genus 723（含 6 个回填 / 18 个改值）、Species 2383 个非空格，涉及 3415 contig（16.3%）。Species 覆盖掉 13.8%，而这正是「已知 vs 新病毒」这篇要用的精度层。

**ACVirus 的质量实测（`measure_tool_consistency.py`，属可查时与参照库比对）**：

| 工具 | 报科数 | 属可查 | 科一致率 | 域一致率 | 种首词=属 |
|---|---|---|---|---|---|
| ACVirus | 8564 | 8564 | 100.0% | 100.0% | 100.0% |
| metabuli | 9062 | 6421 | 100.0% | 100.0% | 91.1% |
| VITAP | 3279 | 2263 | 100.0% | 100.0% | 100.0% |
| CAT | 6968 | 3077 | 100.0% | 100.0% | 89.1% |
| mmseqs | 4928 | 1947 | 100.0% | 100.0% | 99.6% |
| diamond_lca | 7292 | 3457 | 99.7% | 100.0% | 89.0% |
| genomad | 842 | 228 | 100.0% | 100.0% | 100.0% |

七工具的科一致率几乎全是 100%，域一致率全 100%，**度量没有区分度**，无法据此排出一个有意义的平票顺序（按科一致率排出来的顺序与现行顺序基本重合，只换 diamond_lca 与 mmseqs）。ACVirus 三项全满且属可查率 100%，把平票判给它有实测支撑；顺序敏感性（S1 2374 / S2 3357 个 contig）以披露处理。

**敏感性形态**：S1/S2 下各阶元非空数几乎不变（Realm 19015→19015，Species 17239→17069 / 17044），变化形式是候选之间换边。最坏情形下域级组成变动不超过 2.5% 的 contig。

### 16.4 方法章可直接粘贴的段落

> Taxonomic assignments from the seven classifiers were consolidated by a rank-wise cascade consensus (taxonomy gate v6.6). At each rank, from realm to species, only classifiers that had not been overruled at a coarser rank were allowed to vote, and each classifier cast at most one vote per rank; vote weighting was disabled, so a rank was accepted for a contig when a single value held the strict majority of the votes cast at that rank. Classifiers reporting any other value at an accepted rank were overruled and excluded from all finer ranks, were never reinstated, and contigs for which no classifier reported a value at a rank were left unassigned at that rank without penalising any classifier. Ties (6.8% of all rank-contig votes) were resolved in two steps: first by the agreement of the values the tied classifiers reported at the next, finer rank, and, only when that was uninformative, by a fixed tool priority (ACVirus > CAT > VITAP > diamond_lca > geNomad > Metabuli > MMseqs). Permuting this priority changed 11-16% of individual contig assignments (2,374 contigs for the most adverse single-tool move, 3,357 for the fully reversed order) while leaving the number of non-empty assignments at every rank essentially unchanged, so the priority redistributes assignments among plausible alternatives rather than altering the amount of information recovered. Lineage consistency was then enforced by reference-backed gates that harmonise family-genus pairs and require species epithets to be compatible with the assigned genus.

最后一句对应的是「投票定值 / 闸门管自洽」的分工，已在补丁注释里固化：引擎只负责同一阶元内多取值之间定一个值，跨阶元的自洽交给后置闸门，调门槛属口径选择而非 bug 修复。

### 16.5 交付物与上线步骤

新增/定稿文件（本地 `scripts/audit/`）：

| 文件 | 用途 |
|---|---|
| `patch/virus_classifier_analysis.R.cascade_v66f` | 定稿补丁（注释固化版 v6.6） |
| `deploy_cascade_v66f.sh` | 校验 → 备份 → 安装 → 复验，一条命令 |
| `rollback_cascade.sh` | 还原最近一次备份（也可指定文件） |
| `make_tieblank_variant.py` | 生成 v6.9「平票留空」变体（留档，不上线） |
| `make_tiebreak_variant.py` | 生成 S1/S2 顺序敏感性变体 |
| `measure_tool_consistency.py` | 工具一致率度量（16.3 表） |
| `sens_tiebreak_run.sh` | 敏感性两跑 |

步骤：

```bash
# 1. 传补丁（本地）
scp scripts/audit/patch/virus_classifier_analysis.R.cascade_v66f zhangwenda@202.119.189.246:/tmp/
# 2. 服务器上装（备份+校验+复验，不重跑任何项目）
bash /tmp/deploy_cascade_v66f_lf.sh
# 3. 回滚（如需）
bash /tmp/rollback_cascade_lf.sh
```

装完后版本戳报 6.6，已 ≥ 管线 `TAX_GATE_REQUIRED = "6.2"`，因此**已有产物不会被强制重跑**，八个项目保持现状；新建项目自动走新引擎。要让某个已有项目按新口径重出，显式重跑它的 05 阶段；要让全部八个项目一起换，把 `virome_pipeline.py:51` 提到 `6.6` 即可（版本戳不达标会强制重跑，不需要手工清产物）。

**上线执行记录（2026-09-15 23:15，已按大王指令完成）**

| 项 | 值 |
|---|---|
| 安装前线上脚本 | `3e6c0cd81dff1a61714f28ea310df949`（v6.2, 55874 B；与本地 `snapshot/virus_classifier_analysis.R` 一致） |
| 备份 | `virome_discovery_pipeline/virus_classifier_analysis.R.bak_cascade_v66f_20260915_231529` |
| 安装后线上脚本 | `c77c43e2468ad3de48466da3e89fff9e`（v6.6f, 81100 B, LF） |
| 复验 | `parse` 通过；`TAX_GATE_VERSION <- "6.6"` |
| 管线守卫 | `_ver()` 为整数元组比较，(6,6) > (6,2) → 达标，不强制重跑 |
| 回滚 | `bash /tmp/rollback_cascade_lf.sh`（或 `cp -p <备份> <线上脚本>`） |

论文取数仍是旧口径产物（`Votus.integrated/final_integrated_classification.tsv`，md5 `12e8f171…`）。要让 barbarum 按新口径重出（单次约 67 s，产物先落备份目录比对），或要八个项目一起换（把 `virome_pipeline.py:51 TAX_GATE_REQUIRED` 提到 `6.6`），需另行指令。

### 16.6 唯一保留的开关

若 Genus 掉 1364 格（8.5%）不可接受，唯一的杠杆是**允许被淘汰工具在「该阶元全体活跃工具都无值」时回填**（即把「不复活」放宽为「无票时才复活」），一行改动。这与「被淘汰工具不能复活」的原约定相冲突，是否要改由大王定；不改则按 16.1 定稿。

