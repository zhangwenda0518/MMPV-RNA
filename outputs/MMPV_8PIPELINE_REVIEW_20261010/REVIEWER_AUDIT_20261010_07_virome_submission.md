# Virome submission pipeline 审稿视角审查（2026-10-10）

## 0. 总结与边界

**总体判断：当前提交流程存在至少一条已由既有 probe 复现的 CDS 坐标错误，以及多项可造成注释、元数据或导出文件失真的缺陷；相关分支在修正并通过下列验收前，不宜把生成文件直接视为可提交终稿。** 这不等同于全部提交产物都错误，影响取决于所走分支及输入条件。

- 审查对象：`virome_submission_pipeline/` 中统一 discovery/analysis 提交入口、Sequin 构建器、已知病毒参考引导分支、hypothetical 注释、同步回写及 HTML 元数据导出。
- 证据范围：按任务要求停止扩展探针/搜索，仅使用已取得的 probe 结果，并核对当前源代码行号和调用链；没有重新运行外部提交工具或读取真实提交数据。
- 证据分级：**实测**=已有 probe 实际运行结果；**静态高可信**=当前源码中可达路径直接支持；**条件风险**=机制确定但对真实样本的影响范围未测。
- 级别：P0=可直接破坏提交 CDS/序列注释结果；P1=提交文件、元数据回写或导出可靠性问题；R-major=需限制科学解释或审稿口径的问题。
- 本报告不包含实际序列/样本 ID。既有 probe 报告的 450 条 ID 变化仅记数量，不展示值。
- `report_html.py` 当前工作区含用户既有修改，本次只读核对并保留，未修改该文件。

## 1. P 级实现缺陷

### S01 · P0（条件）参考引导 CDS 坐标出现偏移

- **证据与位置**：既有 probe 在一个参考注释样例中观察到应为 `7..15` 的 assembly CDS 坐标写成 `4..12`。参考 CDS 坐标进入对齐列的映射见 `known_virus_submit.py:252-268`；assembly 坐标起点按对齐前非 gap 数累加的实现见 `known_virus_submit.py:313-328`。主流程 `main → gen_features` 可达该逻辑：`known_virus_submit.py:536-542`。
- **后果**：feature table 中 CDS 与实际核酸序列错位，可能改变 CDS 翻译、氨基酸序列与下游注释；这是直接的序列注释完整性问题，不只是格式问题。
- **可信度**：probe 实测，高；当前实现位置与受影响链路静态确认。此处结论限定于已测的映射样例，不据此推定所有参考/对齐均偏移相同数值。
- **验收**：对含参考前端插入、缺失、无 gap、正负链的已知小型对齐样例，逐个独立推导 assembly 起止坐标并核对输出 `.tbl`；所有坐标必须覆盖预期 CDS 碱基，提取序列应与独立期望序列一致，并验证翻译一致。

### S02 · P1 hypothetical 重写丢失未建模 CDS qualifiers

- **证据与位置**：解析器的 CDS 数据结构只保留 `product/inferences/notes`：`analyze_hypothetical.py:65-82`；属性解析仅识别这三类字段：`analyze_hypothetical.py:339-348`；重写时 `to_tbl_lines()` 也只输出位置、product、inference、note：`analyze_hypothetical.py:96-112`，其结果用于每个 CDS 的重写：`analyze_hypothetical.py:1162-1192`。既有 probe 实测 hypothetical 更新输出丢失 `gene`、`locus_tag`、`transl_table`、`codon_start`、partial qualifiers。
- **调用链**：统一入口在发现模式且指定 RVDB 时运行 hypothetical：`submission_pipeline.py:1176-1183`；随后如果更新表存在，优先把该表作为提交 feature table：`submission_pipeline.py:1247-1255`。
- **后果**：即使只修改假定蛋白 product，其他 CDS 的结构/标识注释也会随整表重写而丢失，提交注释不完整，并可能影响后续 feature 对应与翻译解释。
- **可信度**：实测，高；静态调用链确认更新表可进入提交步骤。
- **验收**：输入至少含上述 qualifier 的未命中 CDS 和命中 hypothetical CDS；运行后所有原 qualifier 必须逐 feature 保留，仅预期 product/inference/note 可变；再通过 table2asn/NCBI validator 对照。

### S03 · P1（条件）空 KEEP 产物会回退提交 HQ 全集

- **证据与位置**：`_detect_keep_fasta()` 只接受大于 100 字节的 `class_KEEP.fasta`；未命中后会继续回退 `HQ_plant_viruses.fasta`：`submission_pipeline.py:199-224`。既有 probe 实测：KEEP 文件为空时选择了 HQ 全集。
- **调用链**：discovery 参数缺省时调用该检测器：`submission_pipeline.py:1078-1084`。
- **后果**：当 9b 已运行但 KEEP 结果为空/过小，流程把“KEEP 为空”与“9b 未运行”混为一谈，可能绕过高置信子集筛选，把 HQ 中不属于 KEEP 的序列带入提交。
- **可信度**：实测，高；真实业务中空 KEEP 是否应停止取决于上游语义，属条件风险。
- **验收**：区分“KEEP 文件不存在”和“KEEP 文件存在但为空”：前者才允许按明确策略回退；后者应停止并提示，不得静默提交全集。用不存在、空文件、有效 KEEP 三例验证选择结果。

### S04 · P1 Sequin `.src` 表头与每条数据行列数不一致

- **证据与位置**：Sequin builder 写出的 `.src` 表头为 11 列：`sequin_builder.py:464-468`；数据行由 `:469-476` 写出。既有 probe 实测表头 11 列、数据行 12 列，且中间插入的空字段使后续值错列。
- **调用链**：统一入口在 `--mode sequin/both` 时调用 `sequin_builder.py`：`submission_pipeline.py:1560-1579`；每种病毒调用 `_write_src()`：`sequin_builder.py:336-338`。
- **后果**：source modifier 字段可能错配到相邻列或被工具忽略，BioProject、BioSample、SRA 等值无法按声明的列名可靠解释。
- **可信度**：probe 实测，高；写表代码位置静态确认。
- **验收**：所有输出行的字段数必须等于表头；按列名读取数据，确认每个字段值落在正确列；覆盖有/无位置数据及多序列病毒，再由实际 Sequin/table2asn 输入检查。

### S05 · P1 CSV 回写不能把规范化序列 ID 映射回 `.src`

- **证据与位置**：ID 规范化会替换非字母数字/下划线字符，并对长 ID 截断加摘要：`fix_sqn.py:128-134`；生成 `source_norm.src` 时应用该转换：`fix_sqn.py:275-288`。同步脚本只去除 `lcl|` 前缀，再用余下 ID 查 CSV 原始键：`sync_sqn_from_csv.py:240-256`。既有 probe 观察到所测 CSV 的 450 条 ID 全部发生规范化变化，实际值不在本报告披露。
- **调用链**：同步入口先以 `source_filled.src` ID 集校验 CSV：`sync_sqn_from_csv.py:230-238`，之后分别更新 `source_filled.src` 与供 table2asn 使用的 `source_norm.src`：`:240-256`。因此原始表更新成功并不意味着规范化表也被更新。
- **后果**：用户编辑的元数据可能未写入 table2asn 实际使用的规范化 source 文件；脚本可继续重建 `.sqn`，造成“编辑已同步”的假象。
- **可信度**：既有 probe 实测，高；两个 ID 变换/查找路径静态一致。
- **验收**：在规范化时保存原始 ID↔规范化 ID 双向映射并按映射回写；对全部 sequence_name 证明每条恰好匹配一条源记录，零未命中/重复；用含特殊字符和超长 ID 的 fixture 验证最终 `source_norm.src` 中修改值与 CSV 一致。

### S06 · P1 同步脚本可能在 validator 报错时输出成功状态

- **证据与位置**：既有 probe 的同步命令输出 3 条 table2asn `Error:` 诊断但进程返回码为 0。脚本统计错误数后，失败条件只看非零返回码、`.sqn` 不存在或文件小于 10 KB，并未将 `err_n` 纳入成功判定：`sync_sqn_from_csv.py:274-289`。
- **边界说明**：本条只报告同步脚本的实测行为，不泛化成“统一主入口恒 exit 0”。统一入口通用 `run()` 会检查子进程返回码：`submission_pipeline.py:143-171`。
- **后果**：在 probe 所见的“工具 rc=0 但输出含 validator Error”情形下，`.sqn` 达到大小阈值即可打印 `[OK]`；错误诊断数虽展示，却不阻止成功状态。
- **可信度**：特定命令 probe 实测，高；成功条件静态确认。不能据此断言所有 table2asn 版本都以同样方式返回码。
- **验收**：使用含已知 validator Error 的受控输出且 rc=0 的桩测试，必须失败或进入人工复核状态；只有错误为零并满足所需产物/验证条件时才报告成功。真实工具版本的 `.val`/stdout/stderr 结果也应进入判定。

### S07 · P1 `--reference-metadata` CLI 参数未进入元数据加载链

- **证据与位置**：统一入口确实在 `core14` 可用时传入 `--reference-metadata`：`submission_pipeline.py:1610-1617`。CLI 将它保存为 `args.ref_meta`：`unified_metadata.py:829-830`；但主函数只用 `args.metadata` 调用 `load_metadata()`：`unified_metadata.py:861-864`，而不是 `args.ref_meta`。既有 probe 证实此参数无消费引用。
- **后果**：通过统一提交入口传入的公共参考元数据不会用于 CSV 填充；日期、地点、BioSample、平台等字段可能留空/保留占位符，用户容易误以为元数据已经读取。
- **可信度**：静态高可信，既有零引用检查与实参/形参链一致。
- **验收**：从统一入口传入带唯一哨兵值的 reference metadata，生成 CSV 必须反映对应字段；不传入时不得出现该哨兵值；验证日志显示实际读取文件及匹配数。

### S08 · P1 HTML 编辑/CSV 导出会静默截断超过 80 字符的字段

- **证据与位置**：`safe_row()` 默认最大长度 80 并将长值裁成前 77 字符加省略号：`report_html.py:164-172`；源表及 cmt/bs 表都调用该函数构造 HTML 行：`report_html.py:233-249`。浏览器 CSV 导出基于表格 DOM 当前内容：`report_html.py:68-99`，因此导出的值也是被截断后的呈现文本。既有 probe 观察到从 120 字符原值导出为 80 字符展示值。
- **后果**：HTML 中编辑/下载的 CSV 不是原 CSV 的无损副本；长 note、地点或其他元数据可能被截短后再次写回，造成不可见的数据损失。
- **可信度**：probe 实测，高；当前实现静态确认。
- **验收**：至少 80、81、120 字符和含逗号/引号/非 ASCII 的值，下载 CSV 后须逐字节/逐字段等值；如 UI 要限长，应保留完整原值供导出，并明确标示显示省略。

## 2. R-major 科研数据/解释风险

### SR01 · R-major 已知病毒分支不能识别 ERR/DRR run accession

- **证据与位置**：已知病毒分支的 run 提取仅匹配 `[SC]RR\d+`，不匹配 ERR/DRR：`known_virus_submit.py:175-179`；元数据生成路径再次以 `[SC]RR\d+` 识别：`known_virus_submit.py:480-493`。既有 probe 实测 ERR/DRR 输入落入 `UNKNOWN`。
- **调用链**：`main → gen_metadata → run_suvtk_sqn`：`known_virus_submit.py:525-555`；生成 source 行并填入 SRA/来源元数据见 `:440-521`。
- **科研后果**：该分支处理 ENA/其他归档来源的 ERR/DRR run 时，样本与真实 run 元数据的关联可能丢失，输出 `UNKNOWN` 或使用回退值；提交记录中的原始数据来源和样本字段因此不可靠。真实队列中此类 accession 的比例未量化。
- **可信度**：probe 实测，高；覆盖影响属条件风险。
- **验收**：将 SRR/CRR/ERR/DRR accession 作为独立测试，所有合法前缀都应提取并与 metadata 表准确联接；未知格式必须显式告警且不得冒充完整填充；输出 coverage 按 accession 类别汇总。

### SR02 · R-major 已知病毒独立分支仍把作者/平台等提交信息固定为默认值

- **证据与位置**：`known_virus_submit.py` 生成 submit block 时在代码中直接构造作者/联系/题名字段，调用接口未暴露逐项目作者参数：`known_virus_submit.py:348-375`；MIUVIG/Assembly 表也直接写入测序平台、测序方法、软件及方法：`:377-393`、`:506-520`。这些值不一定适用于每个输入项目。此处不抄录代码中的具体个人或单位信息。
- **调用链**：已知病毒主函数完成 CDS 后直接生成 metadata 并制作 `.sqn`：`known_virus_submit.py:536-555`；该分支独立于统一入口的通用参数。
- **科研后果**：多项目/不同作者或测序平台输入可能被提交为错误归属/方法元数据；这是记录真实性与可追溯性风险，而非单纯格式占位。
- **可信度**：静态高可信；真实多项目的错误覆盖范围未测试。
- **验收**：作者、联系人、平台、建库/组装方法从必需配置/输入显式获得，并在缺失时中止或标为待人工确认；至少用两套不同配置确认生成的 `.sbt`、MIUVIG/Assembly 文件分别反映各自输入。

### SR03 · R-major 内部终止密码子不会被提交前 gate 检出

- **证据与位置**：gate 对 CDS 检查长度为 3 的倍数、边界及起止密码子，但只检查首/末密码子：`submission_pipeline.py:321-333`、`:405-433`。既有 probe 以内部 TGA 的 CDS 测试，gate 返回通过。该段落没有将序列按遗传密码表翻译后检查内部 stop。
- **职责判断**：不将它简单定性为“必须拒绝所有内部 TGA 的 validator bug”：部分病毒可能存在程序性读穿/重编码等特殊情况，需按生物学背景判断。但对按标准遗传密码表解释的 CDS，内部 TGA 是需人工核查的强信号；当前 gate 给出的“通过”不能证明 CDS 翻译完整。
- **科研后果**：内部 stop 可能代表错误 CDS 边界、移码/组装错误或特殊翻译机制；未经核查的提交会传播不完整或错误蛋白注释。
- **可信度**：gate 对内部 TGA 未报警为 probe 实测，高；其是否为真实注释错误取决于序列与病毒生物学，属条件风险。
- **验收**：按明确遗传密码表翻译 CDS 并报告内部 stop；常规情况至少 WARN，按项目策略阻断或要求人工确认；对程序性读穿/重编码提供显式例外与可追溯依据。验收 fixture 覆盖标准 CDS、内部 stop、partial CDS、负链及例外注释。

## 3. 证据限制与不作过度结论的事项

1. 已有 probe 观察到 `sync_sqn_from_csv.py` 的工具输出含 3 个 Error 但 rc=0；本报告只据此评价该同步脚本的成功判定，不宣称整个统一入口总是错误退出。
2. 已知 `.val` 文件没有被该提交链显式解析/汇总，但当前证据不足以断言实际 `.val` 内容存在错误，故不把它单列为已证实缺陷。应在验收中验证 `.val` 是否生成、错误/警告数如何与退出状态及最终成功状态关联。
3. 未重新下载/运行 suvtk、table2asn、Sequin 或 NCBI validator；未用真实提交集估计 ERR/DRR 覆盖比例，也未推断所有 `.sqn` 都不可提交。
4. 本报告不审查 `report_html.py` 之外的其他报告模板渲染行为；`report_html.py` 保留工作区原有改动。

## 4. 投稿前优先验收清单

| 优先级 | 条目 | 最低验收结果 |
|---|---|---|
| 立即阻断受影响分支 | S01 | 参考引导 CDS 坐标及提取序列与独立期望一致 |
| 修复后才可用相关输出 | S02–S07 | qualifiers 保留；空 KEEP 不静默扩集；`.src` 列数和值对齐；规范化 ID 可回写；错误诊断不能报成功；reference metadata 真正加载 |
| 防止数据二次损坏 | S08 | HTML 下载不丢失长字段 |
| 提交前研究记录核实 | SR01–SR03 | accession 可联接，作者/平台来自项目事实，内部 stop 经生物学核查 |
| 证据补齐而非预设错误 | `.val` | 明确 `.val` 解析/人工复核策略，并与实际工具输出对拍 |
