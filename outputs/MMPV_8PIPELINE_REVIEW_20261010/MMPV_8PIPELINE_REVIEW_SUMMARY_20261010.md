# MMPV-RNA 八条管线审稿视角综合审查

**日期：2026-10-10**  
**范围：** public metadata、data preprocessing、virome discovery、virome analysis、EVE、phylogeny、submission、host expression 八条管线。  
**定位：** 这是各管线独立审查报告的综合导航和跨管线判断；逐项证据、行号、触发条件及验收要求以各目录报告为准。

## 0. 先给结论

这套代码有不少值得保留的工程基础：阶段化编排、配置继承、统一 I/O 布局、部分输入/输出契约、可选模块隔离、失败日志、校验和科学方法局限的文档披露。但从审稿角度，**主要风险不是“脚本跑不跑”，而是多个地方可能让结果口径、样本集合或完成状态与最终表述不一致**。

**在修正或解释前，不宜把相关结果直接写成确定性科学结论的优先问题：**

1. **统计量有可复现的数值错误。** 系统发育模块的 Tajima’s D、Hd、Fst 及 FEL/MEME 阳性判据存在问题；分析模块的 Pi 不是通常定义的 pairwise nucleotide diversity；宿主富集 ORA 的检验族可能使 q 值翻转。
2. **某些样本/候选会被漏掉、误筛或被当成成功。** 公共检索默认参数不兼容；GSA 双端下载可能成功后立刻停止、少下一个 mate；预处理会删除仍被宿主表达管线需要的 clean reads；发现管线会把失败 stage 记录下来却仍可返回成功。
3. **方法门控与执行链未完全相连。** BEAST 可在 DRT 未通过/未执行时继续；部分过滤结果没有真正进入后续候选源；空 high-confidence 表可回退到更宽松结果；checkpoint/`--force` 不总能让旧结果失效。
4. **科学证据和代码字段语义有错位。** Salmon EM 估计数被填进 `Uniq_Reads`；用相对 TPM 与 pooled AF 得到的相关/“单倍型”不能自动解释为唯一 reads、绝对剂量或真实病毒相位单倍型；EVE 某读框偏移可把完整 ORF 判为退化候选。
5. **提交链存在可改变序列注释的缺陷。** 已知病毒参考引导 CDS 坐标样例发生偏移；hypothetical 重写会丢失部分 CDS qualifiers；source 表列错位及 ID 回写失配可令最终提交元数据不等于用户编辑内容。

上述结论**不代表所有历史分析、所有样本或所有提交都错误**。每项都有具体触发条件；未对真实历史队列重新量化的内容已标为条件风险/科研改进，不外推影响比例。

## 1. 八条管线报告索引

| 管线 | 独立报告 | 条目数与本报告摘要的级别 | 首要复核点 |
|---|---|---|---|
| ① 公共数据与宿主参考 | [`REVIEWER_AUDIT_20261010_01_public_metadata.md`](REVIEWER_AUDIT_20261010_01_public_metadata.md) | 10 项：高优先级 3、中优先级 6、科研可追溯性 1 | `--detailed` 参数失配、PE mate 漏下、`--force` 未传入子建库 |
| ② 数据预处理 | [`REVIEWER_AUDIT_20261010_02_data_preprocessing.md`](REVIEWER_AUDIT_20261010_02_data_preprocessing.md) | 5 项 P1、2 项 R-major | 去宿主后删除 clean reads、multi-lane 首文件选择、交接目录混源 |
| ③ 病毒发现 | [`REVIEWER_AUDIT_20261010_03_virome_discovery.md`](REVIEWER_AUDIT_20261010_03_virome_discovery.md) | 5 项 P1、2 项 R-major | stage 失败仍可成功结束、cluster 检查点恒真、严格过滤未贯穿 |
| ④ 已知病毒分析 | [`REVIEWER_AUDIT_20261010_04_virome_analysis.md`](REVIEWER_AUDIT_20261010_04_virome_analysis.md) | 1 项 P0、6 项 P1、3 项 R-major | Pi 定义、缺 MD 变 0、空高置信筛选回退、AF callability |
| ⑤ 系统发育与分子进化 | [`REVIEWER_AUDIT_20261010_06_virome_phylo.md`](REVIEWER_AUDIT_20261010_06_virome_phylo.md) | 8 项 P0、1 项 P1 | MAFFT 假成功、多个群体统计量、DRT/BEAST 门控 |
| ⑥ 内源性病毒 EVE（模块已归档 2026-10-11） | [`REVIEWER_AUDIT_20261010_05_endogenous_virus.md`](REVIEWER_AUDIT_20261010_05_endogenous_virus.md) | 3 项 P1、4 项 R-major | qframe 偏移导致完整 ORF 被判退化、Stage2/3 集合断裂 |
| ⑦ GenBank 提交 | [`REVIEWER_AUDIT_20261010_07_virome_submission.md`](REVIEWER_AUDIT_20261010_07_virome_submission.md) | 1 项 P0、7 项 P1、3 项 R-major | CDS 坐标、qualifier 保留、KEEP 空集合、CSV 回写 |
| ⑧ 宿主表达（模块已归档 2026-10-11） | [`REVIEWER_AUDIT_20261010_08_host_expression.md`](REVIEWER_AUDIT_20261010_08_host_expression.md) | 3 项 P0、4 项 P1、5 项 R-major（其中富集占位符同时属实现/统计空间问题） | 网络/索引缓存失效、ORA 检验族、报告虚报显著 |

八份独立报告合计 **73 个审查条目**（按报告条目数相加、未跨管线去重）。不同报告的级别命名并非完全同一套分档，因此不将它们误加成一个精确的全局 P0/P1 统计。P 级表示实现/流程问题，R-major 表示研究设计或解释限制；后者不等同于已发现代码 bug。

## 2. 跨管线共同模式

### A. “有文件/有标记”不等于“当前输入的正确产物”

- 宿主表达网络切换病毒标签、改参考或 `--force` 后，旧 network/index 仍可能被 `--resume` 或目录扫描复用。
- discovery 的 cluster `.ok` 校验将 `bool(Path.glob(...))` 当文件存在性检查，空迭代器对象仍为真。
- analysis 空过滤输出未落盘时，`get_summary()` 退回未过滤 `best_summary`。
- public hostdb 外层强制重跑未同步清除子进程 checkpoint。

**审稿影响：** 日志/参数显示本次运行的配置，文件却可能来自旧样本、旧参考或旧筛选口径。应以输入文件指纹、参数与数据库版本共同定义缓存有效性，且结果报告写出实际消费路径与哈希。

### B. 失败、缺失、未检出和零值在不同节点被混为一谈

- discovery stage 的异常默认记录后继续，最后可正常结束；主进程未基于失败列表返回非零。
- EVE Stage2 的零命中 query 可能不写 verdict 行，Stage3 只看已有行，汇总却将缺失 verdict 显示为 `undetermined`。
- AF fallback 无 AF 文件的样本被跳过；未记录位点按参考态编码是显式策略，但 QC `missing=0` 不能证明可调用性。
- analysis 缺 MD 标签时异常可被吞并为零统计；未检测的 run 可能被标签为病毒阴性。
- public 数据中 CollectionDate 用 release/load date 回填后，Core14 丢失来源标记；系统发育消费者优先读 CollectionDate（`data_collector.py:158-164`），存在把代理日期当采样日期的链路风险。

**审稿影响：** 研究报告应至少区分 `observed zero / not tested / missing / tool failed / imputed`，不得用一个 `0` 或 `success` 代替全部状态。

### C. 输入/输出边界的契约有错位

- 预处理自动清理 clean reads，而宿主表达管线需要它们；report 同时收集 clean 与 host-depleted reads，不保证只交接一套完整 PE。
- public 下载链接循环读取 `url_1/url_2`，任一个成功即停；若是 R1/R2 链接，成功状态不等于双端齐全。
- discovery 逐样本 filter 镜像目录与 COBRA resolver 期待路径错位，可能回退识别阶段原始候选。
- submission CSV 规范化 ID 与同步脚本的反查 ID 不对齐，编辑过的 source 值未必进入 table2asn 使用的规范化文件。
- 系统发育默认删除重组子后，部分 CAPHEINE/BEAST 消费者仍可能读原/旧 alignment。

**审稿影响：** 论文 Methods 所称某一筛选/过滤/去重应用于“后续全部分析”时，需要按消费者逐个列样本 ID 与输入文件；不能依靠阶段顺序推定。

### D. 统计名词和代码字段名未必表示同一个量

- Salmon `NumReads`（EM 分配估计数）被填入 `Uniq_Reads` 并与 `Multi_Reads=0` 配套。
- 已知病毒分析 Pi 代码计算的是 `1 - p_max`，不是 pairwise diversity。
- phylo Tajima’s D 混用每位点 π 与绝对 θW；排除参考后 Hd 仍遍历未排除集合；Fst gap 剔除后群体标签错位。
- ORA 只对出现查询命中的 term 做 BH 校正，改变检验族；HTML 还将非显著 term 标作显著。
- QST AF token clustering 不含位点相位证据，不能直接称真实 haplotype/准种组成。
- EVE RepeatMasker alignment 行数不等于独立拷贝数；`distributed_decay` 空模型富集有限，不应单独当 EVE 定论。

**审稿影响：** 需要让列名、图标题、检验定义和 Methods 同步；任何指标重命名/计算修正后，应评估受影响历史结果是否重算。

## 3. 最先处理的“改变结论”清单

按对序列集合或统计结论的直接影响排序（行号见对应管线报告）：

1. **修复直接算错的统计量/筛选判据：** phylo 的 Tajima’s D、Hd、Fst、FEL/MEME q 判据；analysis Pi；host ORA 的 BH 检验族。建议先用独立手算 fixture 定义期望值，再重跑受影响统计。
2. **修复定年和比对前置门控：** phylo batch 原始 FASTA 必须实际比对；DRT 低随机化次数不可用排名百分位伪装显著；BEAST 要受 DRT 结论/状态约束。
3. **修复真实序列判别链：** submission 的参考引导 CDS 坐标；EVE 的 qframe 重复偏移；discovery strict CDD 与后续最终候选的衔接。对 EVE 受影响候选及提交 feature table，必须回看原始核酸序列和注释。
4. **修复把门槛绕开的空结果语义：** high-conf 为空不得回退未筛选 summary；空 KEEP 与尚未运行 KEEP 阶段必须区分；公共 PE Run 只有两个 mate 都齐全才算完成。
5. **修复阶段与缓存状态：** stage 失败最终非零退出；cluster 产物校验真实枚举 glob；输入/参数变化使各管线 checkpoint 失效；明确 XML-only、报告跳过、待人工状态，不把它们统一显示为完成。

## 4. 投稿前必须补的设计与证据

1. **样本单位与独立性：** 检查 Run/BioSample/BioProject 重复；技术重复、组织、物种、项目和测序批次要进入分析设计表。每组“30 个 Run”不自动等于 30 个独立生物学个体。
2. **选择偏倚与检测阴性：** 对预筛、去宿主、rRNA、覆盖度阈值、read rescue 和 HQ 筛选，报告各步分子/reads 流失；低丰度 spike-in、阴性对照、不同阈值敏感性分析比只展示默认阈值更能回答审稿疑问。
3. **LLM 元数据处理：** 保留原始字段、规范化值、推断/回填来源、模型/提示词版本及逐 Run 修改记录；不应覆盖非空原始值后只输出规范化 Core14。
4. **病毒宿主表达：** unknown 检测与真阴性分开；co-infection 需声明；同一队列先选 DEG/模块再对同一队列做载量相关不能写成独立验证；组织内与项目内混杂应由匹配、协变量、holdout 或独立队列处理。
5. **群体遗传和系统发育：** 对 AF fallback 明示 missing/reference 编码、callability 数据、indel SNP 域；说明 pooled frequency 不携带相位；重组 drop/mask 策略与每个消费者的样本集需可追溯。
6. **EVE 与新病毒声明：** 将算法自动分类、人工复核与独立整合/结构证据分层报告；不把候选筛选结果等同实验验证/已确认整合或完整新病毒基因组。

## 5. 本轮验证到底做到哪一步

- 全仓轻量语法编译覆盖 **438 个 Python 文件**；唯一报告的语法错误位于 `virome_analysis_pipeline/archive/_patch_html.py:243`，是归档路径，不据此推定现役入口失败。该扫描不是功能测试。
- 多条重点结论来自先前 AST 抽取/内存桩/手算对照，例如 cleanup 虚拟 FS、align 不调用 MAFFT、ORA 检验族反例、EVE qframe 反例、submission 坐标/列数反例及 phylo 的统计量反例。相关管线报告注明具体可信度。
- 未跑完整 pipeline、真实服务器数据库、NCBI/GSA/Salmon/DIAMOND/MAFFT/IQ-TREE/BEAST/table2asn 等外部工具，也未联网核实外部数据库版本；未估计每项缺陷命中真实样本/结果的比例。
- 本轮为只读审查：**未修改既有 Python/YAML/文档源文件**。新增产物为八份逐管线报告和本综合报告；仓库原有工作区改动不归因于本审查。

## 6. 建议执行顺序

- **第一批：** 修复 phylo/analysis 的统计定义、DRT→BEAST 门控、提交 CDS 坐标、EVE qframe 与宿主 ORA 检验族；为每项保留旧/新实现对照和最小回归用例。
- **第二批：** 修复 CLI 参数契约、双端下载、clean 保留/交接、discovery failure/strict/cluster 检查点、各模块输入指纹与缓存失效。
- **第三批：** 补真实数据敏感性分析：按项目/宿主/组织/独立 BioSample 分层；spike-in/阴性控；CALLABILITY/缺失率；不同阈值及参考版本；对筛选后样本集输出逐阶段名单与哈希。
- **发布论文前：** 将实际运行 config、软件/数据库版本、sample manifest、排除/缺失/回填决策、可复跑测试结果归档；逐句对照 Methods 与实际运行记录，不要从默认配置直接复制成论文方法。
