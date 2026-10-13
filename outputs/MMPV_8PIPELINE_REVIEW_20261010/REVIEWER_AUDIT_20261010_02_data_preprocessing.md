# 数据预处理管线审稿视角审查（2026-10-10）

## 总体结论与范围

范围：`data_preprocessing_pipeline/` 的 clean、host depletion、可选 BBNorm 与 report 入口，并与 ⑧ 宿主表达的输入契约交叉核对。**结论：主流程的核心步骤与风险有清晰文档，但有几处可能丢失样本/lane、销毁可复用输入或把不完整结果写成交接成功的问题；统计留存率和读段配对也需按“实际 read 数”而非样本文件数验证。**

本报告仅使用当前源码行号和本轮已取得的内存桩/AST探针。没有执行 fastp、Kraken2、Bowtie2/HISAT2、RiboDetector 或 BBNorm，也没有运行真实数据。严重级别 P1=可能破坏样本/阶段数据契约或误导下游；R-major=投稿前需补充的科研验证/解释，不等于已证明真实队列受影响。

## 1. 实现与流程问题

### PP01 · P1：去宿主后默认删除 clean reads，与宿主表达管线输入契约冲突

- **证据**：`data_preprocessing.py:394-395` 每当 deplete 正常返回就调用 `cleanup_00a_after_deplete()`；该函数在 `data_preprocessing.py:248-298` 依据去宿主目录中出现的样本名删除 `1.fastp_tmp/`、`2.fasta/`、`3.clumpify/` 中同名文件。已取得的 AST 函数探针以虚拟文件系统执行原清理函数，记录到 R1/R2 四次 `unlink` 请求；未触碰真实数据。
- **跨管线调用链**：`host_expression_pipeline/README.md:56-60,158-163,169-172` 与 `doc/IO_LAYOUT_DESIGN.md:140` 将 ② 的 clean reads 作为宿主转录组 Salmon 定量输入。去宿主 reads 不含宿主表达所需 reads。
- **影响**：执行 ② `--stage all` 后，clean FASTQ/FASTA 可能不可用于后续 ⑧；若使用者以去宿主 reads 替代，会改变宿主表达估计对象，无法回答原研究问题。该问题作用于采用默认自动清理且之后需要 clean inputs 的工作流；不是说去宿主本身算法错误。
- **可信度**：调用路径静态确认；原函数以虚拟文件系统实测请求删除。
- **建议/验收**：清理设置须可关闭，且与下游契约协调（例如保留 clean 原始交接，或要求显式确认已归档/不使用 ⑧）。验收主链完成后 `01_CleanData` 中必要输入仍可供 `host_expression.py --project-root` 读取，或系统明确以可恢复方式拒绝清理。

### PP02 · P1（多 lane 条件）：同一 sample name 只选首个 R1/R2 文件

- **证据**：`clean-data.py:135-144` 先按 `extract_sample_name()` 聚合文件，再以 `next(...)` 选第一个符合 R1/R2 标签的文件；每侧其它匹配文件不进入 `tasks`（`149-153`）。`clean-data.py:441` 主流程消费该任务列表。
- **触发与后果**：同一生物样本有多个 lane 且命名规范化后归为同一 sample 时，只有排序后首个 R1/R2 被处理，其他 lane 不合并，测序深度与病毒检出能力下降。此条是条件风险；真实输入多 lane 比例未核实。
- **可信度**：源码路径直接可见；报告不声称真实项目必含多 lane。
- **建议/验收**：将同样本全部 lane 按对应 read end 合并/逐 lane 处理后再汇总；dry-run 显示每个样本接入的全部文件。用两个 lane 的 PE fixture 验证输入 read 数守恒且 R1/R2 配对一致。

### PP03 · P1：报告交接目录合并扫描去宿主与 clean 读段，未实现文档所称优先回退

- **证据**：`preprocess_report.py:173-195` 声称“去宿主优先、回退 clean”，但实现将 deplete 根与 clean 的 `3.clumpify/2.fasta` 全部加入 `roots` 并累积匹配文件；`build_rows()` 在 `:202-207` 直接把这个结果交给后续汇总。
- **影响**：同一 sample 可能同时出现去宿主与未去宿主输入；双端半成品或不同格式的文件也可能一起进入最终列表。交接表中的“最终 reads”不再唯一对应真实下游输入语义，容易重复纳入或把含宿主的 clean reads 当作已去宿主 reads。
- **可信度**：静态路径确认；不声称任何真实运行已生成重复条目。
- **建议/验收**：严格按每 sample 先选完整去宿主 PE/SE，否则才回退 clean；分别校验 R1/R2 完整性，并在 TSV 标出来源阶段。fixture 同时放 clean 与 hostdep，输出必须只选一套；只有 R1 的 PE 不得标为完整 PE。

### PP04 · P1：标准 Bowtie2 对齐率日志格式可能无法被报告解析

- **证据**：`preprocess_report.py:143-160` 的 regex 是 `overall alignment rate\D+([\d.]+)`，即预期百分比出现在“overall alignment rate”之后；Bowtie2 常见 stderr 字面格式为 `95.50% overall alignment rate`（百分比在前）。既有标准输出内存串探针导致该正则未捕获比对率；没有运行外部 Bowtie2。
- **影响**：report 中 `Aligner_Rate(%)` 为空，映射率/QC报告不完整，不能将“没提取到”解释为工具报告了 0%。
- **可信度**：解析器当前正则静态确认；探针为内存合成标准格式文本，不是外部工具实测。
- **建议/验收**：支持已记录的 Bowtie2/HISAT2 日志格式并在无法解析时显式写 `UNPARSED`；将真实版本日志加入 fixture，正确解析 95.50，而不是空值。

### PP05 · P1（多 lane 条件）：seqkit 阶段 reads 数被后续行覆盖而非汇总

- **证据**：`preprocess_report.py:86-113` 在每个 seqkit 行按 Sample/Stage 得到 key，再执行 `out[sample][key] = n`；同样 sample/stage 多 lane 或 PE 两端会覆盖前一行，没有累加。报告计算留存率使用这些单一值（`:202-245`）。
- **影响**：基数和阶段 read 数可能只取最后一个文件/端的 count，留存率在不同 lane 数或 PE/SE 样本间不可比；同一报告会将两端总 raw reads 与单行 seqkit 值混在一起。
- **可信度**：源码路径静态确认；实际偏差取决于 summary 文件行粒度，未用真实日志量化。
- **建议/验收**：明确每个 summary 行的单位与 R1/R2/lane 聚合规则；在样本级汇总配对端 read 数，并用两个 lane、双端的手算 fixture 确认分子分母同口径。

## 2. 科研设计与投稿解释（R-major）

1. **宿主/rRNA过滤与病毒召回需要已知阳性、阴性及 spike-in 验证。** Kraken2 confidence、病毒 taxid rescue、aligner 的 `SAM.view_filter` 配对策略和 RiboDetector 都会决定保留哪些 reads。`host_depletion.py:798-800` 与工具调用链可确认过滤发生，但本次没有外部工具对拍或真实队列召回率数据。应提供病毒低丰度 spike-in、宿主阴性和已知病毒 reads 对照，报告按病毒类群/深度分层的灵敏度与宿主残留率。README/METHODS 的双端处理描述不能替代召回验证。
2. **clean 转 FASTA 后无法再做质量值依赖的变异证据。** 当前 clean 输出语义和 host expression 使用 FASTA 的设计已公开；这不影响 Salmon 伪比对本身，但如后续把同一 FASTA 当作变异 caller 输入，碱基质量证据已丢失。投稿应区分用于定量/组装的 FASTA 与用于变异调用、保留原始质量值的 FASTQ/BAM。

## 3. 值得保留的设计与限制

- Kraken2 宿主税onomic node 解析未误报：已取得的标准 6 列 Kraken report 小样例正确得到 `{10239,12345}`；不列为缺陷。
- `SAM.view_filter` 的 `-f 12` 是严格要求 pair 两端均 unmapped 的策略，并非语法错误；它可能减少病毒 read 召回，应以 PP/R-major 的真实阳性验证量化，而不是直接改成宽松过滤。
- 对外部工具的实际版本、数据库内容与性能未做端到端验证；本报告没有据此断言外部工具本身失效。

## 4. 优先级

| 优先级 | 条目 | 核心处理 |
|---|---|---|
| 先处理跨管线可用性 | PP01 | clean 产物不得在未协调下游前被默认销毁 |
| 输入完整性 | PP02、PP03 | lane 全量处理；交接列表阶段/PE 完整性唯一明确 |
| 结果可信度 | PP04、PP05 | 对齐率真实可解析；计数单位和留存率口径一致 |
| 投稿补证 | R-major 1–2 | spike-in/阴性对照；区分定量/组装 FASTA 与质量感知变异输入 |

## 5. 覆盖边界

本报告覆盖 7 个管线 Python 文件的语法编译及所列相关函数的静态审查；收到的已有虚拟文件系统探针验证了 clean 自动清理，但没有重新运行任何测试。未执行 fastp、Kraken2、Bowtie2、HISAT2、RiboDetector 或 BBNorm；未统计真实队列中多 lane/PE 文件数量与实际漏失 reads 数，也未检查服务器数据库/工具版本。上述科研风险须以真实样本统计和控制实验补证，不是对当前全部输出失效的断言。
