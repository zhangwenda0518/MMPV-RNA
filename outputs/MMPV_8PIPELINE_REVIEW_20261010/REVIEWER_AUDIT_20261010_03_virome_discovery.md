# virome_discovery_pipeline 审查报告（2026-10-10）

## 0. 结论与审查边界

**结论：当前代码包含数项可导致阶段失败被上层报告为成功、过滤条件绕过或断点续跑丢失已完成结果的 P1 流程缺陷；此外，分支 C 的“completo”判据不足以支持参考基因组完整性结论，Methods 对 30/30/40 分值与最终裁决关系的表述也与实现不一致。** 这些风险不等于已证明既有样本结果全部错误；报告区分确定的代码路径与真实队列中的影响范围。

审查范围仅限 `virome_discovery_pipeline/` 本次已确认的源码路径及团队此前已取得的静态/AST/受控桩证据。本轮未新增探针、测试运行、依赖安装、数据库/真实样本运行或外部搜索。代码行号按当前工作树读取。

证据等级：
- **静态高可信**：由当前源码可达调用链直接确定；不代表已运行外部生信程序。
- **既有受控证据**：团队此前的 AST/代理桩结果，本报告仅引用既有结论，不声称本轮重跑。
- **科研条件风险**：代码中的机制确定，但对真实队列及论文结论的实际影响尚未量化。

严重级别：**P1**=流程/输出可靠性缺陷；**R-major**=可能改变科研结论、完整性解释或方法可复现性。未将下列科研建议写成“必然造成全部结果错误”。

## 1. P1：实现与流程可靠性

### V01 · P1：默认容错继续，但主入口最终仍可返回成功状态

- **证据与调用链**：`--stop-on-error` 是可选开关，help 明确默认“仅警告，继续执行” (`virome_pipeline.py:3140-3141`)。阶段循环捕获 `SystemExit` 与一般异常，默认只记录 warning 并将阶段名加入 `failed_stages` (`virome_pipeline.py:3906-3919`)；循环结束仅打印失败阶段 (`virome_pipeline.py:3923-3924`)，没有按 `failed_stages` 设置非零退出码。
- **触发条件**：未传 `--stop-on-error`，任一阶段以 `sys.exit(1)` 或抛异常失败；主进程继续后续阶段，最终可正常返回。
- **后果**：调度器/批处理系统可能把不完整的目录当作成功产物；下游可能在旧文件、缺失文件或部分结果上继续，失败信息只留在日志中。此问题不是“某个子程序返回非零就一定漏检”：单阶段函数可以显式 `sys.exit`，但主循环会捕获它。
- **可信度**：静态高可信，默认参数及异常处理链清楚；本报告未新增运行探针。
- **修复及验收**：允许继续执行以收集诊断，但只要 `failed_stages` 非空，主入口最后必须以非零码退出，并在摘要中列出失败阶段。以一个抛 `SystemExit(1)` 和一个抛 `RuntimeError` 的阶段桩验证：默认继续后最终非零；全部成功时为零；显式 stop-on-error 仍立即非零退出。

### V02 · P1：cluster 阶段的产物校验恒真，旧 `.cluster.ok` 可掩盖缺失结果

- **证据**：`_stage_validated('cluster')` 的兜底分支写成 `bool(d['root'].glob('*/4_centroids/final_centroids.fasta'))` (`virome_pipeline.py:433-440`)。`Path.glob()` 返回迭代器对象，迭代器本身布尔值恒真，并不表示存在匹配文件。阶段跳过逻辑在 `.ok` 存在且验证函数为真时跳过 (`virome_pipeline.py:381-400`)。
- **触发条件**：输出根目录有旧 `.cluster.ok`，但相应 `4_centroids/final_centroids.fasta` 已不存在或从未生成，且契约路径检查也没有命中。
- **后果**：cluster 被错误判为已完成并跳过；后续 taxonomy/host/rescue 可能失败，或误消费其他历史位置的文件。
- **可信度**：高；另有团队既有证据：cluster 第 439 行的 AST/受控空迭代器检查确认空匹配时 `bool(Path.glob())` 仍为真。本轮未重做该探针。
- **修复及验收**：对 glob 结果做实际存在性判断，例如 `any(...)` 并检查文件大小/FASTA 有效性。测试无匹配、仅目录有匹配、真实有效文件三种状态：只有有效目标产物可使阶段跳过。

### V03 · P1：单样本或 `--no-merge` 时 COBRA 可能绕过 filter，静默回退原始鉴定候选

- **证据与真实调用链**：`run_filter()` 在候选只有一个或显式 `--no-merge` 时调用逐样本过滤 (`virome_pipeline.py:849-854`)；逐样本输出位于 `d_filter/<sample>/` (`virome_pipeline.py:927-955`)。过滤脚本把 COBRA 镜像写到其 `outdir/candidates/<sample>/...` (`filter_virus.py:533-576`)，因此逐样本路径会嵌套在 `d_filter/<sample>/candidates/...`。但 `_cobra_candidate_source()` 只在 `d_filter/candidates/*/*_virus.all.candidate.fasta` 查找 filter 镜像 (`virome_pipeline.py:990-994`)；graph 候选不可用时直接返回 identification 目录且 warning 为 `None` (`virome_pipeline.py:995-1007`)。`run_cobra()` 随后把该来源传给 COBRA (`virome_pipeline.py:1045-1059`)。
- **触发条件**：单个待过滤样本（即使没有显式 `--no-merge`，长度为 1 也走 per-sample），或多样本显式 `--no-merge`；graph 候选未提供替代镜像。
- **后果**：filter 阶段虽执行并产出过滤结果，COBRA 输入却可能仍来自 02a 原始鉴定目录；strict/filter 的筛选不会约束 COBRA 扩展。由于回退分支不附带 `src_warn`，日志不一定提示实际绕过。
- **可信度**：静态高可信；路径布局、镜像位置和候选源查找层级均由源码直接确定。
- **修复及验收**：统一逐样本与合并模式的镜像契约，或让候选源定位器递归/按样本目录解析，并在 fallback 前显式警告。以单样本默认过滤和多样本 `--no-merge` 各跑一个轻量目录桩，断言 COBRA 收到的 FASTA 与 filter 最终候选一致；filter 结果缺失时必须告警或失败，不能静默回退。

### V04 · P1：strict 模式生成的 CDD strict 结果未进入最终 UniProt/候选筛选链

- **证据与调用链**：`filter_fa_path` 初始为 `None`，只有 `args.mode in ('filter','both')` 时才赋为 CDD filter 文件 (`filter_virus.py:653-674`)。strict 模式另行生成 `<base>_cdd_strict.fasta` (`filter_virus.py:675-681`)，但随后的 UniProt 输入是 `filter_fa_path or current_input` (`filter_virus.py:683-695`)；strict 下 `filter_fa_path` 仍为 `None`，所以 UniProt 对原始输入执行。最终 `split_output` 仅把 `filter_fa_path` 作为 CDD 结果传入，严格输出则优先使用 UniProt 结果 (`filter_virus.py:697-710`)。
- **触发条件**：`--filter-mode strict`，且 CDD 与 UniProt 均启用（strict 模式默认如此）。
- **后果**：代码虽计算 CDD strict FASTA，但其通过/剔除结果不约束最终 strict 候选集；实现效果不是 help 所称的“UniProt(仅TaxID)+CDD(仅viral)”联合严格筛选 (`virome_pipeline.py:3178-3180`)。被 CDD strict 拒绝但通过 UniProt 的序列仍可能进入最终输出。
- **可信度**：静态高可信；不依赖 CDD/DIAMOND 软件行为，差异来自输入路径选择。
- **修复及验收**：明确 strict 的组合逻辑，并把 CDD strict 输出接入 UniProt 输入或将两套结果按声明规则取交集。构造一条 CDD strict 排除、UniProt strict 通过的序列及一条两者均通过的序列，验证最终 FASTA 分别被排除/保留；同时核对 merged 与 per-sample 两种模式。

### V05 · P1：branch C 组级断点续跑没有把既有 pass 序列恢复到内存集合

- **证据**：断点阶段会从 `branchC_pass.fasta` 和 `branchC_checkpoint.txt` 把已有 ID 加进 `done_groups` (`rescue_pipeline.py:1052-1064`)，据此从 `pending_groups` 排除已完成组 (`:1073-1078`)；但随后 `pass_records` 仅初始化为本轮重新识别的 `completed_pass` (`:1065-1067`)，没有把已有 pass FASTA 中完成组的序列读回。成功组仅将新序列追加到当前进程的 `pass_records` (`:1253-1258`)，结尾以该内存列表重写 `branchC_pass.fasta` (`:1269-1270`)。
- **触发条件**：已有若干属组完成并写入 pass/checkpoint；恢复运行时还有其他待处理组（或直接 pass contig），后续产生新的 pass 记录并执行最终写回。
- **后果**：输出文件被本次内存中的直接通过与新完成组覆盖，之前已完成但本轮跳过的组序列消失；checkpoint 却仍可能保留这些组的完成标记，后续恢复不会再补回。
- **可信度**：静态高可信；涉及的断点加载、跳过、追加、重写路径均明确。本报告未新增模拟中断实验。
- **修复及验收**：恢复时把已有 pass FASTA 记录载入并与本轮新记录按稳定 ID 合并，或采用可校验的组级独立产物汇总，避免仅重写内存增量。验收：完成组 A 后中断，继续处理 B；最终 FASTA 同时含 A、B 且无重复；连续恢复多次结果集合不变。

## 2. R-major：科研解释与方法口径

### R01 · R-major：branch C 将 query coverage ≥98% 直接解释为参考基因组“completo”并通过

- **证据与调用链**：每条 BLAST hit 的 `qcov` 是 `alen / qlen`，其中 `qlen` 来自待判 contig 的自身长度；读出的 hit query length 变量未用于判定 (`rescue_pipeline.py:981-994`)。最佳 hit 按 bitscore 选取后，只要 `best_qcov >= 0.98` 就标为 `completo` (`:998-1007`)，没有同时要求参考覆盖率或最低 identity。此状态在 Phase 2 直接加入 pass (`:1030-1038`)，branch C pass 后续进入 A/B/C/D 合并及最终去冗余 (`:1937-1970`)。
- **科学风险**：一个短片段若几乎全长对齐到较长参考，query coverage 仍可接近 100%，但这不能证明参考基因组主体已覆盖；低 identity 的全 query 局部命中也可能达到该门槛。因此“completo”可能把局部/不完整序列当作完整序列直达 pass，影响完整性与 HQ catalog 解释。是否已有真实序列命中这一边界，尚未量化。
- **可信度**：判定机制静态高可信；实际命中比例为科研条件风险，未用真实数据验证。
- **改进及验收**：分开报告 query coverage、subject/reference coverage、identity 与比对区间；若目标是完整基因组，应结合参考长度/预期节段长度或独立完整性证据，避免仅凭 query coverage 命名 `completo` 并直达 pass。用“短 query 完整落在长 reference”及“近全长、高一致性 reference 对齐”两个对照验证：前者不得标完整，后者按预先定义门槛通过；报告中将 query 覆盖与基因组完整度分开表述。

### R02 · R-major：Methods 声称 30/30/40 加权分数产生 verdict，代码实际按门控规则裁决、分数只排序

- **证据**：Methods 写“weighted evidence scoring … yielding KEEP/REVIEW/DROP verdicts” (`METHODS_Virome_Discovery.md:48`)。实现计算 BLASTN、BLASTX 与 domain 的 30/30/40 总分 (`integrate_rescue_evidence.py:345-354`)，但 KEEP/REVIEW/DROP 在随后独立的 CDD 证据、BLAST identity/coverage 与长度条件分支中赋值 (`:355-422`)，并仅按 `total_score` 对输出行排序 (`:483-484`)。CT3-HMM 还存在单独的 DROP→REVIEW rescue 及显式启用时提升 KEEP 的规则 (`:424-453`)。
- **影响与边界**：论文方法文字会使读者误以为 verdict 是加权总分阈值的结果，无法从当前 Methods 复现实际门控；但“分数仅用于排序、CDD 为门、BLAST 辅助定性”可以是有意设计。本条**不要求**总分必须直接决定 verdict，也不把门控本身判作缺陷。
- **可信度**：实现与文字差异静态高可信；对已发表/待投稿结论的影响取决于结果是否按该 Methods 描述解释。
- **改进及验收**：若门控是目标设计，修订 Methods 为与实现一致的 CDD/BLAST/HMM 分层规则，并说明 30/30/40 分数仅排序；若总分本来应参与裁决，则需另行给出预设阈值及完整决策表。验收：用固定输入逐条映射 Methods 到代码分支，KEEP/REVIEW/DROP、HMM rescue/promote 与分数排序的职责无歧义。

## 3. 未列为缺陷的已确认策略与边界

- **Graph 默认关闭不报缺陷**：主流程、参数 help 和配置均明确为可选子步骤；是否启用属于分析策略，不是实现意外 (`virome_pipeline.py:762-765, 857-860, 3066-3068`)。
- **CDD 无结果保留不报缺陷**：无 CDD 命中/未分类候选的保留行为是筛选策略；本报告没有把“没有 CDD 结果”本身写成错误，也没有要求把所有未命中序列删除。
- **30/30/40 不直接裁决不单独报 bug**：仅记录 Methods 与实现的口径不一致；验收允许保留当前门控设计，不要求分数必须决定 verdict。
- COBRA 的 `virus_mode=raw` 单独看不是绕过滤器：主链预期由 orchestrator 先定位过滤后的候选镜像，再以 raw 模式读取；V03 指向的是单样本/per-sample 镜像实际路径与定位器预期不一致。
- 不根据缺少真实数据库、未运行第三方软件或未量化真实样本影响提出推断性 bug。

## 4. 优先级摘要

| 级别 | 项目 | 主要验收目标 |
|---|---|---|
| P1 | V01 默认失败仍最终成功退出 | 阶段失败后最终退出码非零 |
| P1 | V02 cluster glob 产物校验恒真 | 无真实 centroid 时不得跳过 |
| P1 | V03 单样本/不合并绕过 filter 输入 | COBRA 实际消费 filter 镜像或明确失败 |
| P1 | V04 strict CDD 结果未进入最终链路 | strict 候选符合声明的 CDD+UniProt 组合规则 |
| P1 | V05 branch C checkpoint 丢已完成组 | 中断恢复前后 pass ID 集合一致 |
| R-major | R01 query coverage 被当参考完整度 | 短片段不得因 query 覆盖高而称完整 |
| R-major | R02 分数裁决的 Methods 与实现不一致 | Methods 可复现实际门控，且不误称分数决定 verdict |

**本报告未声称**运行了完整端到端分析，也未判定真实 HQ 序列集合已有多少比例受影响。