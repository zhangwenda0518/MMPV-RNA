# 内源病毒流水线审查（2026-10-10）

## 0. 审查边界与结论

审查对象：`endogenous_virus_pipeline/` 的位点筛查、DNA 病毒/EVE 判别及元件级后处理。共列 **7 条**：3 条条件性实现/流程缺陷，4 条科研解释或数据口径风险。本文沿用主代理已取得的 AST probe 结果，只核对当前源码及行号；**未重新搜索、运行测试、重算或加载真实数据**。因此，凡涉及实际队列命中规模的内容均标为条件风险或既有校准结果，不推断影响数量。

级别：P1=可能改变候选归类或漏掉下游证据；R-major=影响科学解释、计数或独立性。可信度区分为既有 probe、静态代码路径与条件性影响。

**总体判断：** 当前流程的结构判定有一项已有 probe 揭示的高影响读框坐标问题；此外，Stage2 零命中与 Stage3 候选集合衔接不完整。即便修复实现问题，RepeatMasker 行数不能直接解释为拷贝数，候选/证据文件也缺少一致性校验；EVE 退化类判据的既有空模型校准显示区分度有限。`eve_distinguish/README.md` 已披露部分科学局限，本文将其作为投稿口径边界，而非重复报成代码 bug。

## 1. 实现与流程缺陷

### E01 · P1：BLASTX 的 qframe 偏移被重复应用，可能把完整 ORF 判成分布式退化

- **证据**：`eve_distinguish/s1_decay_scan.py:216-228` 从 qframe 计算 `offset=abs(fnum)-1`，正链扫描起点又取 `qs-1`（负链也按 HSP 坐标换算）；`scan_frame_stops()` 在 `:134-141` 再以 `start + offset` 开始扫描。即在 HSP 坐标已定位到的扫描区间上再叠加读框偏移。
- **触发与影响**：已有 AST probe 对合法 HSP/完整 600 nt 编码段验证：原序列没有终止密码子，但该实现计出 199 个并归为 `distributed_decay`。这可令原本 `KEEP_virus` 或 `REVIEW` 的候选沿 `s3_verdict.py:277-286`、`s4_filter.py:41-48` 转成 EVE 方向的 verdict/action（包括 `MOVE_EVE`）。该结论来自既有 probe，本次未重跑。
- **可信度**：既有 AST probe + 当前源码路径，高；影响的是符合所测坐标条件的非零 qframe 命中，真实候选数未知。
- **验收**：明确 qstart/qend 与 qframe 的坐标契约，保证读框偏移只用于一次坐标换算；为正/负链、qframe ±1/±2/±3 和区间裁切建立人工可核算的用例。合法完整 ORF 应得到 0 个 premature stop，且不再因坐标偏移进入 `distributed_decay`。

### E02 · P1（条件）：Stage2 零命中位点虽在汇总里记作 undetermined，却不会进入 Stage3

- **证据**：`eve_scan_core.py:435-471` 的 `verdict_loci()` 只遍历 Stage2 原始结果中实际出现的 query（`sorted(hits)`），没有为“零命中”query 写 `undetermined` 行。相反，`summarize_genome()` 在 `:524-544` 对缺失 verdict 使用默认值 `undetermined`。Stage3 的 `stage3_annotate()` 在 `:922-941` 只从已有的 S2 verdict 行收集 `viral_supported` / `undetermined`，缺行的候选不会被加入集合。
- **触发与影响**：某 Stage1 位点对植物/病毒拆分库没有任何 S2 命中时，单基因组汇总会显示其为 `undetermined`，但 S3 RVDB 检索候选中缺少它；汇总的“未判定”与实际深度注释对象不一致。该路径是否影响具体运行取决于 S2 零命中位点数。
- **可信度**：静态代码路径，高；真实触发频率未知。
- **验收**：S2 应以 Stage1 全部 locus 为输出域，未命中显式写 `undetermined`；或由 Stage3 以 Stage1 locus 列表补全集。构造一条 S1 存在、S2 零命中的候选，确认它在汇总与 Stage3 输入中都可追踪，计数一致。

### E03 · P1（条件）：S3 对候选与外部 evidence 表按 ID 静默连接，不校验覆盖或来源

- **证据**：`run_all.sh:201-207` 只检查 query/evidence 文件非空，允许分别指定；`s3_verdict.py:143-166` 分别读入 S1/S2/S2b/evidence 后只遍历 S1 的 ID，并以 `ev.get(cid,{})` 取证据；`:196-206` 缺行时 `tax_family`、蛋白/核苷酸证据都退化为空值，没有不匹配告警或失败闸门。
- **触发与影响**：误传旧批次、不同候选集合或 ID 命名不一致的 `rescue_evidence_scored.tsv` 时，缺失的候选被当成“没有上游病毒科/同源证据”，仍继续给 verdict/action；若有不相同但同 ID 的旧表，则可能将错误分类证据套到当前序列。当前没有证据表明本次实际输入发生过错配，故不据此声称现有结果已错。
- **可信度**：静态高可信；触发与实际影响为条件风险。
- **验收**：启动前比较候选 ID 与 evidence ID 覆盖率，要求每个候选恰有符合契约的证据行（如确需缺省，显式标记 `missing_evidence` 并保持 REVIEW）；记录 query/evidence 的来源标识或内容指纹。ID 不匹配时不得静默产出可执行 MOVE/REMOVE。

## 2. 科研设计与计数口径

### R01 · R-major：RepeatMasker 命中行被直接当作拷贝数，且覆盖碱基数为逐行求和

- **证据**：`eve_element_db.py:488-512` 解析 RepeatMasker `.out` 的每一行；`accumulate_quant()` 在 `:515-527` 对每个命中行直接 `copies += 1`、`covered_bp += end-begin+1`。注释称碎片校正交给 regions 阶段（`:518-519`），但 `stage_regions()` 在 `:634-660` 从筛查 summary 独立生成 `element_regions.tsv`，没有被 `stage_rm()`（`:806-821`）用于计算 `family_quant.tsv`。
- **触发与影响**：同一共识元件在一个位点被 RepeatMasker 分段命中，或同家族命中区间重叠时，行数/逐行碱基和不等于独立整合拷贝数/非冗余覆盖度；当前 `family_quant.tsv` 的 `copies` 仍是 alignment-row count。真实偏差规模尚未量化。
- **可信度**：静态实现确定；实际偏差取决于 RepeatMasker 输出的碎片化与重叠情况。
- **验收**：保留 query/contig 与坐标，按可追溯规则合并同一位点的碎片并对覆盖区间取并集；分别输出 `hit_segments`、`loci/copies`、`covered_bp_union`，并在已知单拷贝但多段命中的对照上验证拷贝数仍为 1。若暂不校正，至少将列名/报告改为“命中片段数”。

### R02 · R-major（条件）：同一 assembly 以不同名称重复输入时，可伪增跨基因组支持数

- **证据**：`eve_screen.py:264-274` 只拒绝清洗后重名的 batch name，没有比对路径、assembly accession 或序列内容。元件聚类在 `eve_element_db.py:361-379` 以 metadata 中的 genome 名称集合计算 `n_genomes`；因此同一 fasta 用不同名称登记时会被计成多个 genome。
- **触发与影响**：batch 中同一 assembly 被别名重复登记，或同一生物样本/组装被多个名称重复计入时，跨 genome 簇支持数、科×基因组分布及基于 genome 名的支持度被抬高；不能将命名条目数自动解释为独立物种/独立生物学样本数。没有检查实际 batch，故为条件性科学风险。
- **可信度**：代码路径高；是否触发需核对实际 manifest。
- **验收**：记录 assembly accession 与文件 hash；启动前报告相同 hash/ accession 的重复项，需显式说明是否保留；汇总同时给出条目数、唯一 assembly 数与唯一生物学样本/物种数。

### R03 · R-major（范围受限且已披露）：污染清理是事后且靶向有限，未自动阻止污染候选流入后处理

- **证据**：`eve_qc_contam.py:15-18` 明确只针对 phiX174 / 大肠杆菌噬菌体并将其它噬菌体排除在脚本范围外；检测规则见 `eve_scan_core.py:991-1000`，报告生成 `eve_qc_contam.py:219-269` 产出建议清单。与此同时，元件收集 `eve_element_db.py:324-340,555-590` 按 verdict/bit score 收集候选，没有读取该 QC 清单作为过滤条件。
- **触发与影响**：若筛查库中的噬菌体命中来自建库 spike-in、细菌污染或其他非目标来源，且未由用户在下游应用建议清单/额外 QC，仍可能进入跨基因组聚类和元件库。**不应据此把所有噬菌体命中视为污染**；问题在于该后处理是窄范围提示，不是已接入的普适污染排除闸门。
- **可信度**：流程边界静态高；真实污染占比未重新估计。该范围限制已在模块说明中披露，属于解释/流程风险，不是新发现的隐藏行为。
- **验收**：报告必须区分已筛查靶标、未覆盖污染源及未执行筛查；对可能用于物种组成/元件库的 Phage 类结果提供细菌/载体交叉核验或独立来源支持。过滤需留存规则与名单，不能自动把未知噬菌体全删。

### R04 · R-major（既有校准结果）：`distributed_decay` 对随机组成的增量区分度有限，不宜单独作为强 EVE 定论

- **证据**：现有 `eve_distinguish/README.md:408-430` 记录的 1,027 条输入空模型/同源区置换结果中，`distributed_decay` 真实比例 45.4%、空模型 37.6%，富集 1.21×；文档估算真实触发中约 17% 超出随机组成，并指出依赖退化通道的结论证据弱于架构/编码保留通道。判别与动作路径为 `s3_verdict.py:277-286` 及 `s4_filter.py:41-48`，其中 `ancient_EVE` 映射为 `MOVE_EVE`，而 `EVE_suspect` 为 REVIEW。
- **科学影响**：该校准不证明某条候选不是 EVE，但表明退化分类仅略高于该置换零模型；将退化类和具有整合位点/结构证据的类合并表述，或仅凭 `ancient_EVE` 的自动 `MOVE_EVE` 就称为已确认 EVE，会夸大现有证据。该限制已在 README 中披露，不能冒充未披露的实现 bug。
- **可信度**：转引仓库既有校准结果；本次未复算，也未独立评估其零模型设计。
- **验收/投稿口径**：在独立验证集上报告 specificity、false-positive rate 与置信区间；比较替代空模型/阴性集合，并对退化类单列。未取得外部验证前，将退化驱动的 `MOVE_EVE` 明确标为自动筛出、需复核的候选，避免与宿主侧翼/独立整合证据确认的 EVE 等同。

## 3. 优先级摘要

| 优先级 | 条目 | 核心处理 |
|---|---|---|
| 首先修复 | E01 | 解决 qframe 与 query 坐标重复偏移；按正负链及 6 框回归 |
| 结果完整性 | E02、E03 | 补齐零命中 S2 行；拒绝或显式标注 evidence 错配 |
| 结果口径 | R01 | 区分命中片段、拷贝数与非冗余覆盖度 |
| 投稿前核查 | R02–R04 | 审计重复 assembly；声明污染筛查覆盖；退化类分层验证与谨慎措辞 |

## 4. 已有证据与限制

- E01 沿用主代理提供的既有 AST probe 摘要；其摘要还报告探针环境错误已排除、29 个 Python 文件编译通过。本次没有重跑 probe 或编译。
- 本报告仅对当前源码位置做了只读核对。R02 是否发生需对真实 batch/assembly manifest 去重；E02 的影响量需有真实 Stage2 零命中输出；E03 需拿 query 与 evidence 原始 ID/来源对拍；R01 偏差大小需对 RepeatMasker 结果按位点复核。
- R04 的数值是仓库既有验证文档所载结果，未独立审核置换模型、数据抽样或统计不确定性。
