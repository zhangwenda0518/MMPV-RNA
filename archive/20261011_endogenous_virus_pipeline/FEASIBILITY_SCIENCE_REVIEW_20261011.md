# EVE 内源性病毒管线 · 可行性 + 科学性审查（2026-10-11）

**审查对象**：`archive/20261011_endogenous_virus_pipeline/`（原 `endogenous_virus_pipeline/`，平台 ⑥）全链——
正向筛查 `eve_screen.py` / `eve_genome_scan.py` / `eve_scan_core.py`、元件级后处理 `eve_element_db.py`、
判别链 `eve_distinguish/`（s1/s2/s2b/s2c/s3/s4 + run_all.sh）、QC 与报告 `eve_qc_contam.py` / `eve_report.py`，
共 22 个 Python（≈16.3k 行）+ 1 个 485 行 shell。

**审查方式**：静态读码 + **最小实测**（可复跑的探针见文末"验证清单"）。
每条结论标注 **✅已实测** / **推断**。审查为**只读**，未修改任何源文件。

**与既有报告的关系**：本报告独立于 `outputs/MMPV_8PIPELINE_REVIEW_20261010/REVIEWER_AUDIT_20261010_05_endogenous_virus.md`
（3 P1 + 4 R-major）。凡与旧报告重合处均注明；凡旧报告条目经本轮复核**不成立**的，在附录 A 明确撤回。

---

## 0. 总体评价

**先说做对的**（这部分很扎实，是这套代码真正的资产）：

- **口径纪律好**：默认值从 `eve_scan_core` 单一来源 import（`eve_screen.py:224-231`），表头抽成模块级常量由测试直接 import，不两头各抄一份。
- **失败不静默**：`eve_screen.py:357` 批次有 error 即返回非零；`eve_scan_core.py:475-478` 对不在 id2div 表的 sseqid 主动告警；`eve_distinguish/run_all.sh` 用 `run_or_die` 打印 `FAILED(rc)`+完整命令行再中断。
- **判据被"量"过**：`eve_distinguish/verify/` 有整套校准——整条洗牌零假阳性、阈值逐个扫描、组件判定用 hmmscan 独立复核。
- **主动披露局限**：README 如实写"分布式退化仅 1.21× 富集""12 个阈值里 4 个对结论零影响""组件复核未定论"。这种自我披露在科研代码里很少见。

**但结论是：归档决定（"方法不成熟"）成立，且本轮找到了一个旧报告未能量化到位的核心判据缺陷。**

系统性模式一句话：**"判据写对了、代码按判据执行了"这两件事被反复验证过，但"判据依赖的坐标语义"没被验证过**——
`s1` 的读框换算就是一个典型：它自带了完整的"坐标换算"单测（`tests/test_eve_distinguish_logic.py:203-229`），
但那个单测是**照实现重抄一遍公式**，证明的是"抄对了"，不是"公式对"。真值一测就翻。

---

## 1. P0 — 会导致结论错误

### P0-1 ✅已实测 · `s1` 把读框偏移**用了两次**，约 2/3 的候选被读成错误读框

**位置**：`eve_distinguish/s1_decay_scan.py:218`（`offset = abs(fnum) - 1`）、`:221` / `:224`（`scan_s = qs-1` / `L-qe`）、
`:228`（`scan_frame_stops(sseq, scan_s, scan_e, offset)`）、`:137`（`i = start + offset`）

**真值（本轮实测，两个独立工具）**：用 NCBI BLAST 与 DIAMOND 2.2.6 分别对同一组"把同一蛋白放进 6 个读框"的最小查询跑
`--outfmt 6 ... qstart qend qframe`：

| qframe | +1 | +2 | +3 | −1 | −2 | −3 |
|---|---|---|---|---|---|---|
| qstart（1-based） | 1 | 2 | 3 | 660 | 660 | 660 |
| qend | 660 | 661 | 662 | 1 | 1 | 1 |

两者**完全一致**：`qstart - 1 ≡ |qframe| - 1 (mod 3)`。也就是说 **`qs-1` 本来就已经落在该读框的密码子栅格上**。
代码在此基础上又 `start + offset`，等于把栅格整体挪了一位/两位：

- `offset=0`（±1 框）→ 恰好正确；
- `offset=1`（±2 框）→ 实际扫的是 **±3 框**的密码子；
- `offset=2`（±3 框）→ 实际扫的是 **±2 框**。

而且 `:230-233` 的"另两框基线"用的是同一批错位起点，于是**真实病毒读框被算进了 `other` 基线**——
分子（`stops`）虚高、分母（`other_mean`）被真值拉低，两个方向的偏差叠加。

**端到端实证**（把 DIAMOND 真实输出行喂进 `hsp_metrics`）：

```
查询     真读框  真终止数   | 代码 premature_stops_region | 代码 decay_class
L_f1      1       0       | 0                          | coding_intact
L_f2      2       0       | 2                          | coding_intact
L_f3      3       0       | 35                         | distributed_decay   ← 反了
L_neg1   -1       0       | 0                          | coding_intact
L_neg2   -2       0       | 2                          | coding_intact
L_neg3   -3       0       | 35                         | distributed_decay   ← 反了
```

另一条随机 400 密码子、**全程无终止密码子**的干净编码段（真实读框 +2）：
代码报 `premature_stops_region=17`、`enrichment=1.714`、`decay_class=distributed_decay`；
正确口径应为 `0` / `0.5` / `coding_intact`。

**为什么这会改结论**：`coding_intact ↔ distributed_decay` 正是 `s3_verdict.py` 二维表的主轴——

- `coding_intact` + 位点齐全 → **`virus_candidate`（KEEP_virus，要留的真 DNA 病毒）**；
- `distributed_decay` + 位点齐全 → **`EVE_STRONG_provirus` / `EVE_suspect` / `ancient_EVE`（MOVE_EVE）**。

所以这个 bug 的方向是**把真病毒判成内源化石**，与判别链想回答的问题正好相反。
它同时解释了作者自己校准里那个反常数字：`distributed_decay` 只有 **1.21×** 富集——
一部分"退化"本来就是错读框读出来的。

**验收**：`qs-1` / `L-qe` 已经是栅格起点，`scan_frame_stops` 不应再加 `offset`；为正负链 × 3 读框 × 区间裁切建立
"人工可核算"的用例（**期望值必须来自 blastx 实测坐标，不能照实现重抄**）；修后 `coding_intact` 与 `distributed_decay`
的分布会显著移动，`eve_distinguish/README.md` 与 `doc.md` 里全部 v6.x 数字需整体重算（作者已在 README 里写明
"面板换代与 s1 修复都会改归类，不能只换局部数字"，这条纪律要照做）。

### P0-2 ✅已实测 · `s3` 对 evidence 表按 ID 静默连接，**缺行 → 走唯一会删数据的分支**

**位置**：`eve_distinguish/s3_verdict.py:160-163`（`for cid, d1 in s1.items(): ... e = ev.get(cid, {})`）、
`:198`（`upstream_viral = fam not in UPSTREAM_NA`）、`:82-83`（`UPSTREAM_NA` 含空串 `""`）、`:239-250`（REMOVE 分支）

**实测**（同一候选，只换 evidence 表）：

| evidence 表 | 结果 verdict | action |
|---|---|---|
| **缺该候选行** | `host_contamination_likely` | **REMOVE_host_contamination** |
| 有该行且 `tax_family=Caulimoviridae` | `host_conflict_review` | REVIEW |

即：候选集由 `s1` 驱动（`:160`），evidence 只是 `dict.get` 补字段；**缺行不报错、不告警**，
`fam` 落空串 → `upstream_viral=False` → 若该候选又有自身物种寄主同源（≥90%/0.80）且无 MP/CP/AP，
就落进 `host_contamination_likely`——**判别链里唯一允许删数据的分支**。

**触发**：evidence 表来自发现管道的 `10_Reports/rescue_evidence_scored.tsv`（或 `09b_...` 原位文件）。
`run_all.sh:201-207` 只检查 query/evidence 文件非空、且允许 `-q`/`-e` 分别显式指定 → 误传旧批次、不同候选集合、
或 ID 命名不一致的表，都会**静默**把候选变成"可删的寄主污染"。本轮未发现实际发生过，但没有任何机制阻止它发生。

**验收**：启动前比较候选 ID 与 evidence ID 覆盖率；ID 不匹配时不得静默产出可执行 `REMOVE`/`MOVE_EVE`
（至少显式标 `missing_evidence` 并保持 REVIEW）；记录 query/evidence 的来源标识或内容指纹。
（与旧报告 E03 独立一致。）

---

## 2. P1 — 结果完整性 / 健壮性

### P1-1 推断（静态路径，与旧报告 E02 独立一致）· Stage2 零命中位点进不了 Stage3，汇总却显示 `undetermined`

- `eve_scan_core.py:435-471`：`verdict_loci` 的输出域 = **S2 原始 blastx 命中里出现过的 locus**（`for locus in sorted(hits)`）。
- `eve_scan_core.py:538`：`summarize_genome` 对缺失 verdict 默认 `["", "undetermined", ...]`。
- `eve_scan_core.py:922-930`：`stage3_annotate` 的 `keep` **只从 verdict 表已有行**收集 `viral_supported`/`undetermined`。

→ Stage1 有命中、但 S2 双侧库零命中的位点：**汇总里被写成"未判定"，实际却拿不到 RVDB 深度归属**。
"汇总的 undetermined 计数"与"实际送 RVDB 的候选集"不是同一个集合。验收：S2 以 Stage1 全部 locus 为输出域，
未命中显式写 `undetermined`；或由 Stage3 以 Stage1 locus 列表补全集。

### P1-2 ✅已实测 · 归档后 `eve_report.py` 无法运行

`eve_report.py:33-37` 工作区版本导入 `mmpv_common.report_style`，回退路径是 `Path(__file__).resolve().parent.parent`
（= `archive/`），该层没有 `mmpv_common`。实测 `python eve_report.py --help` 直接 `ModuleNotFoundError`。
管线在仓库根时能跑，**归档把这条依赖打断了**。恢复时把目录移回仓库根即可（与 `ARCHIVED_20261011.md` 的恢复步骤一致），
但"归档即不可运行"这一点应写进归档说明。

### P1-3 ✅已实测 · 不变量测试不通过（1/274）

`python -m pytest tests/ -q` → **273 passed / 1 failed / 14 skipped**。
失败项 `tests/test_eve_invariants.py::TestSourceLineEndings::test_no_crlf_in_module_sources`：
`eve_report.py` 工作区副本有 **292 处 CRLF**（其余 4 个模块 0）。
`git show :archive/.../eve_report.py` 是纯 LF，且 `.gitattributes` 已声明 `*.py text eol=lf`——
属工作区副本问题（`git status` 对该文件为 `AM`：作者在暂存之后又改过它），修复成本低，但**当前套件不是全绿**。

### P1-4 ✅已实测 · `credit_components` 的可用区间闸门是"全有全无"

`eve_distinguish/s2_domain_scan.py:103-118`：`avail >= nt_min` 一次性放行，而封顶 `cap = aln // COMP_AA_MIN`
只按**本条 HSP 的比对长度**折算，不看"扣掉已占用区间后还剩多少"。
实测：一条 300 aa 的融合条目比对，在只剩 600 nt 未占用区间时仍记 **3 个组件**（撑起 3 个蛋白本需 900 nt）。
影响面：`best_locus_ncomp` 在 HSP 重叠严重的位点会偏高。**缓解**：作者自己的阈值扫描显示
`PROVIRUS_NCOMP`/`PROVIRUS_SPAN` 在 1,027 条上改 0 条结论（`README.md` 校准第 3 条），故实际 verdict 影响受限，
属口径瑕疵而非结论错误。

### P1-5 推断（子代理静态发现，本轮未逐条复核）· 元件层的两个计数口径

- `eve_element_db.py` `family_quant.tsv` 的 `copies` 直接等于 RepeatMasker **命中行数**、`covered_bp` 为逐行求和；
  注释称"碎片校正交给 regions 阶段"，但 `stage_rm()` 未使用 `element_regions.tsv`。（同旧报告 R01）
- `n_genomes` 按 genome **名称**集合计算，而 `eve_screen.py:264-274` 只拒绝"清洗后重名"，
  不比对路径/accession/内容 → 同一 assembly 换名重复登记会伪增跨基因组支持数。（同旧报告 R02）

---

## 3. P2 — 遗留物 / 一致性

| # | 位置 | 内容 |
|---|---|---|
| P2-1 | `eve_screen.py:366` | `r.get('elapsed')` 在 `skip` 状态下可能打出 `Nones`（`status=='skip'` 时 elapsed 未赋值） |
| P2-2 | `doc.md:176` vs `eve_scan_core.py:360` | 文档说 `s1_hits.bed` 是"BED3 + 第4列内嵌原始行, **共4列**"，实际是 **15 列**（BED3 + 12 列原始行展开）；`assign_best_hits` 的 `len(c) < 15` 守卫与真实格式一致，是**文档写错** |
| P2-3 | `eve_distinguish/run_all.sh` | 已存在的 `panel.fasta` 原样沿用（作者已声明）；面板换代后旧 OUT 目录会**静默**用旧面板 → 新旧数字不可比，靠人记住换 `-o` |
| P2-4 | `eve_distinguish/s1_decay_scan.py:134-142` | `scan_frame_stops` 的 `pos.append(i - start)` 与 P0-1 同源；修 P0-1 时一并理清"相对偏移"的记账 |

---

## 4. 文档 vs 代码一致性核对

| 文档承诺 | 代码实际 | 判定 |
|---|---|---|
| 完成标记按**请求的阶段**逐个查产物（`doc.md` 第二轮修复表） | `eve_genome_scan.py` 有 `STAGE_OUTPUTS` 逐阶段判定 | ✅ 一致 |
| Stage3 = `viral_supported + undetermined` 位点 vs RVDB | 候选集只来自 S2 判定表**已有行**，零命中位点缺失 | ❌ 不一致（P1-1） |
| 汇总 12 列（不含 genome）、汇总表 13 列 | `summarize_genome` 写 12 列、`merge_all` 补 `genome` | ✅ 一致 |
| `s1_hits.bed` 4 列 | 实际 15 列 | ❌ 文档错（P2-2） |
| 归档"原样保留，未做删改" | 目录结构一致，但 `eve_report.py` 工作区版本相对 git 有未提交改动（`AM`，+36/−16） | ⚠ 部分 |
| `pytest tests/ -q` 全绿（文档称 132 用例） | 实测 274 用例、1 失败 | ❌ 不一致（P1-3） |

---

## 5. 修复优先级建议

**第一档（会改变结论，必须先修）**
1. **P0-1 读框偏移**：这是全链唯一"系统性把真病毒推向 EVE"的缺陷。修完必须用
   `eve_distinguish/verify/rerun_diff.py` 在服务器对 1,027 条重算，产出 verdict 迁移矩阵，
   并整体重写 README/doc 里的 v6.x 数字。**判据回归用例的期望值要从 blastx 实测坐标来，不能照实现重抄。**
2. **P0-2 evidence 连接校验**：在产出 `REMOVE` 之前必须有覆盖率闸门。这条直接影响"删不删数据"。

**第二档（结果完整性 / 可运维）**
3. P1-1 Stage2 输出域补齐到 Stage1 全 locus。
4. P1-2 归档说明补一句"`eve_report.py` 依赖 `mmpv_common`，归档位置不可运行"。
5. P1-3 修 `eve_report.py` 行尾（或让该不变量测试豁免归档目录）。

**第三档（口径与可维护性）**
6. P1-4 把 `cap` 改为按 `avail` 折算；P1-5 区分"命中片段数"与"拷贝数/非冗余覆盖度"，并记录 assembly 指纹。

**投稿前口径（不修代码也要写进 Methods/限制）**
7. 沿用作者已披露的边界：`distributed_decay` 富集仅 1.21×，退化通道不能单独作 EVE 定论；
   `host_sample_map` 仅覆盖 23/1000+ 样本，寄主否决天花板很低；`--fast` 丢约 54% 位点；
   `family_of` 是关键词粗分类；loci 是无方向区间（`restore_coordinates` 不请求 qframe）。

---

## 6. 值得保留的设计

- **同链另两读框做本地基线**的思路是对的（P0-1 修好后这个设计才真正发挥作用）——把"序列组成本来就容易读出终止"这一项扣掉。
- **`check_extracted()`**：抓住了 samtools 遇越界区域**退出码仍是 0** 这个真实陷阱（本机 samtools 1.21 实测确认），
  否则脏序列会静默流到 Stage2/3。
- **`locus_key()`** 统一 `:`↔`_`：抽序列工具改写分隔符会让 Stage2 判定被静默丢弃，这条护栏必要且正确。
- **寄主否决三分支**：把"跨物种旁证""两通道矛盾""确证寄主序列"分开，且只有最后一支允许 REMOVE——
  这是"不要用一个结论覆盖两种证据状态"的正确做法。
- **`eve_distinguish/verify/` 的校准套件**：整条洗牌零假阳性、阈值逐个扫描、hmmscan 独立复核，
  是这套代码里最有科研价值的部分，建议随方法一起保留。

---

## 附录 A · 复核记录（撤回的误报）

按"子代理结论必须自己复核"的纪律，以下条目经本轮**实测/读源码复核后不成立，予以撤回**：

1. **"`assign_best_hits` 的双指针会把命中错配到相邻位点（P0）"** —— **不成立（不可达）**。
   双指针在"位点互不重叠且有序"时等价于暴力匹配：按真实路径 `merge_loci → assign_best_hits`
   随机对拍 400 组（766 个位点）**0 不一致**。构造性反例（一条命中跨越两个位点）在真实流程里
   **不可能出现**——任何跨越两个位点的命中都会让 `merge_loci`（`:333`，`s <= cur_e + merge_d`）把它们并成一个位点。
   结论：仅当 `loci.bed` 被外部手工编辑成重叠区间时才可能触发；作为**理论健壮性注记**保留，不作为缺陷。
2. **"`VIRAL_BS` 被传给了植物侧判定"** —— **不成立**。`verdict_loci` 全仓只有 1 个生产调用点
   （`eve_scan_core.py:898-899`），位置参数 `(cfg.host_bs, cfg.viral_bs)` 与形参
   `(host_bs, viral_bs)` 一一对应。
3. **"`s2c` 的 `host_dominant` 分支不可达"** —— **不成立**。`s2c_refgenome_scan.py:182` 明确
   `call = "host_dominant"`，且 `:129` 注释说明它"先于一切"判定。

> 记下这三条的原因：它们都是"看着像 bug 的代码形状"，但一个是**契约保证下的不可达**，
> 两个是**没 grep 调用点就下结论**。下次同类审查先跑 `grep -rn <函数名>` 确认调用链，再判可达性。

---

## 附录 B · 本轮验证做到哪一步

**已实测（可复跑）**
- 全 22 个 `.py` 语法编译：0 错误；`bash -n run_all.sh`：通过。
- 测试套件：`pytest tests/ -q` → 273 passed / 1 failed / 14 skipped（失败项见 P1-3）。
- 入口可用性：`eve_screen.py` / `eve_genome_scan.py` / `eve_element_db.py` / `eve_qc_contam.py` 的 `--help` 正常；`eve_report.py` 失败（P1-2）。
- **读框真值**：NCBI BLAST（本机 2.x）与 DIAMOND 2.2.6 双工具、6 个读框、同款 `--outfmt`，结论一致（`qstart-1 ≡ |qframe|-1 mod 3`）。
- **P0-1 端到端**：DIAMOND 真实输出行 → `s1_decay_scan.hsp_metrics` → 6 读框逐条对拍；另造 400 密码子无终止 ORF 复现 `distributed_decay` 误判。
- **P0-2 端到端**：构造 s1/s2/s2b/evidence 四表，只切换 evidence 表内容，`s3_verdict.py` 产出 `host_contamination_likely` vs `host_conflict_review`。
- **P1-4**：直接调 `credit_components`，构造重叠 HSP 观察封顶行为。
- `assign_best_hits` 可达性：按真实路径随机对拍 400 组 / 766 位点。
- `eve_report.py` CRLF：字节级计数 + `git show` 对比暂存 blob。

**未做（限制）**
- 未跑完整 pipeline，未接真实数据库（`ncbi-virus_ref.pep.dmnd` / `plant_virus_final.dmnd` / RVDB / id2div）与真实基因组。
- 未在服务器对 1,027 条真实候选重算 P0-1 的影响面（需原始 `panel_hits.tsv` 等中间产物）。
- P1-1、P1-5 为静态路径推断（P1-1 与旧报告 E02 独立一致；P1-5 来自子代理，未逐条复核）。
- 未评估 `--fast` / `--overlap` 的真实召回损失量级（作者给的是 −54% 位点，未独立复算）。

**副作用**：仅运行了 `pytest`（更新了本就存在的 `.pytest_cache/`，已被 git 忽略）与内存内探针；
**未修改任何源文件**。探针脚本在 `%TEMP%/eve_probe/`，不在仓库内。
