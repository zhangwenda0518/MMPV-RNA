# 宿主表达管线审稿视角审查（2026-10-10）

## 0. 总体评价与审查边界

审查对象：`host_expression_pipeline/`，33 个 Python 文件、8,526 行。已先读 `README.md`、`teagcn/PROVENANCE.md` 和根目录 `doc/METHODS_TEMPLATE.md:188-194`，再精读主管线、定量/标签/DE/共表达/模块/载量/富集/报告链路。仅新增本报告；未修改源码、旧报告或记忆，未联网、装包、运行 Salmon/TEA-GCN 重计算、加载真实科研数据或执行会产生缓存的测试。

**总体判断：可作为探索性宿主关联分析模块，但当前结果不宜未经修正直接作为“病毒特异、独立验证、显著重布线或因果剂量响应”的投稿证据。** 做得好的地方是独立输入输出边界、clean 未去宿主 reads、明确 TEA-GCN/DE/载量三分支、计数与 TPM 模型分离、样本量/失衡门槛、自边剔除、跨网络按 Target 交集对齐、NaN-safe BH、网络 QC 和上游出处留痕。系统性问题主要是：参数变更没有可靠失效旧产物；公共 run 的独立性/组织项目混杂没有进入正式设计；探索性选择后的关联被当验证；ORA 检验族及注释缺失处理破坏统计口径。

证据等级：
- **实测-AST**：从当前源码 AST 原样抽取函数、语句或表达式，标准库命名空间/只读桩执行；不等同于整条生信端到端测试。
- **静态高可信**：当前主管线可达调用链逐段追踪。
- **科学条件风险**：代码机制已确认，但真实队列是否命中、影响多少样本/基因，未用真实数据量化。

严重级别：P0=能改变结论的实现错误，P1=流程/输出可靠性错误；R-major=影响研究设计或统计解释，R-minor=补充说明/敏感性分析。以下共 12 条，未发现可据实宣称的全新“无依赖即 bug”。

源码快照：按相对路径排序，对每个 `.py` 的相对路径字节及原文件字节依次 SHA-256：`10e0d3c35f9a8a0b26406193beb8419dea1f270593c4c1199d96f23ddd12669b`。

## 1. 实现缺陷

### H01 · P0：换分组/输入无缓存失效；`--force` 仍可能保留旧 RAW 网络

- **证据**：`host_expression.py:233-236` 只判产物存在与 `.ok`；`:426-429` 子步亦无输入/参数指纹；`:624-625` 网络命令无条件带 `--resume`。`teagcn/main/Run_TEA-GCN_mp.py:115-123` 从已有文件个数生成 offset；全基因完成时 offset 等于基因数，不提交任何新任务。`build_coexp_matrices.py:71-73` 自身也仅靠文件存在跳过。
- **原代码摘要**：`if not args.force and ck.is_file() and all(exists): return False`；网络命令末尾固定 `--resume`。
- **触发链**：首次 coexp 完成 → 同输出根换 `--virus`/labels/QC/TPM/阈值 → 无 force 时顶层或子步直接复用。加 force 时矩阵/分区重跑，但如果新旧 selected k 相同，RAW 网络目录名相同，network resume 仍把旧全套文件认作完成；Rank_transform 再把旧 RAW 变成“新”Completed。若 k 改变则产生另一网络目录，可能由唯一目录检查报错，而非可靠重跑。
- **后果**：日志与 CLI 展示新比较，Δz/模块实际上来自旧队列；这与 README `:48-50`“切换口径会覆盖产物”的行为承诺冲突。
- **可信度**：实测-AST，高。真实 `run_networks` 在 `force=True` 仍输出 `--resume`；真实 `_coexp_substep` 给新 virus 且旧产物齐全，返回 False、执行次数 0。
- **建议/验收**：以矩阵、标签、QC名单、比较水平、关键算法参数/版本生成指纹；变更失效全部依赖。force 路径不 resume 旧计算，应写独立新目录并原子替换。验收：同 k、翻转标签、改 TPM 三个对照，各证明 RAW 被重新计算，而非仅 mtime/日志变化；新病毒与旧病毒不可共享完成标记。

### H02 · P0（条件）：新宿主参考可自动消费旧/旁边任意 Salmon index，跨参考矩阵缺失补零

- **证据**：`salmon_auto_pipeline.py:366-374`“有效 index”仅检查 `versionInfo.json`；`:390-408` 先取输出目录已有 index，不核参考；`:416-424` 扫一级子目录的第一个有效 index。主管线 `host_expression.py:352-357` 传新 reference，但 force 仅 `--force-quant`，不传 `--force-index`。定量汇总 `salmon_auto_pipeline.py:733-736` 对不同文件特征取并集，`:697` 对缺失特征补 `0.0`。
- **触发链**：同 quant 输出根切宿主/CDS，或参考父目录存在另一物种 index → stage_quant → 自动 find_existing_index → 对旧 index 定量 → 上游 run_summary 却仍登记请求的新 reference。不同参考的历史 quant 混合时汇总并集补零。
- **后果**：物种/基因注释错配；“未在该参考中”的特征变生物学零表达；DE/coexp/ORA 都可能错。不是“所有首次运行一定错”，是旧索引存在或多索引情况下的确定可达风险。
- **可信度**：实测-AST + 静态高可信。真实 find_existing_index 在给 `new_species/cds.fa` 时返回 `out/salmon_index`。
- **建议/验收**：索引与参考序列 hash、特征ID集绑定，不匹配即拒绝/明确重建；quant.sf 特征集/Length 必须一致，不应把跨参考缺失补零。验收：A/B不同参考、A残留index，B请求不能无提示接受A；异构 quant 集应拒绝。

### H03 · P1：默认并行入口断点续跑仍按完成文件“数量”跳过，而不是文件名集合

- **证据**：`TEA_mp.py:131` 实际为 `tasks = [(gene_idx, gene) ... if gene_idx >= offset]`；`:136` 使用 `imap_unordered`。`TEA_fast.py:89-94` 同样按 gene_idx/offset 后线程乱序完成。`PROVENANCE.md:36-44` 却明确称已修为按文件名集合。
- **触发链**：默认 mp → 中断后完成文件不构成顺序前缀 → 重发同命令（主管线固定 resume）→ 缺失早期基因永远不提交。即使所有原始输出最终顺序完整，也可能单个文件写入中断而仅靠名字被视作完成。
- **后果**：续跑不能修复断点，可能反复失败。**已纠正误报边界**：不能把此问题宣称为“当前 rank 默认静默漏基因”。`rank_transform_mp.py:126-128` 按现有文件建 gene_dict，`:77` 对未完成目标查字典；标准构造缺 B 而 A文件含B时，实测 `KeyError('B')`。因此常规后续是明确阻断，而非静默完整成功。
- **可信度**：实测-AST，高。当前任务表达式对 `PCC_genes=['A','B','C'], offset=2` 只得 `[(2,'C')]`，漏 B；真实 `_extract_worker` 读桩文件确实 KeyError。
- **建议/验收**：按文件名且完整内容/预期目标集判断完成，临时文件写后原子 rename；不以个数作为游标。验收：预置A/C完、B缺，必须只补B；截断A必须重算；更新 provenance 使声明与真实 AST 对齐。

### H04 · P1：`--strata-col` 对 DE 的“共用”承诺未接线

- **证据**：`host_expression.py:945-950` 称 DE 分层按 `--strata-col`，但 `:700-701` 实际传 `--strata args.de_strata`；`run_deseq2.py:226-229` 找不到指定列仅 warning 后跳过。LoadScan `host_expression.py:733` 才真正传 args.strata_col。
- **触发链**：用户按 help 用 `--de-strata enabled --strata-col Batch`，预期按 Batch 防混杂 → DE 收到 `--strata enabled` → 缺列跳过，仍产生 all 结果且整个 DE 正常完成。README `:187` 通用示例的裸 `--de-strata --strata-col Batch` 还会因 `--de-strata` 必须值而 argparse 拒绝。
- **后果**：用户以为执行了项目/组织分层，实际未执行。`--de-strata Tissue` 当前确实能用，应保留这一事实，不笼统称“DE分层不可用”。
- **可信度**：实测-AST，高。真实 parse_args 和 stage_de 捕获命令 `--strata enabled`，没有 Batch。
- **建议/验收**：要么把 de-strata 定义成 boolean并传 strata_col，要么明确它本身就是列名并改 help/示例；显式请求缺失层列应错误退出。验收 Batch 示例真实命令列名正确，并含每层输出/样本数；缺列不可伪装成功。

### H05 · P0：ORA 对观察到命中的 term 才 BH，检验族由查询结果决定

- **证据**：`enrich_analysis.py:173-174` 先获得背景规模符合门槛的 bg_terms；`:189` 却 `tested=[t for t in bg_terms if any(g in gene_set ...)]`；`:192-197` 仅该命中集合进入 hypergeom/BH。文档 `:27-28` 写“BH校正全部检验term”，未说明结果依赖缩小检验族。
- **触发链**：DE/hub/module → stage_enrich → 已知背景term → 先看查询命中再选检验族 → q 被减小。
- **后果**：零命中 term 的 p=1被删，按声明的背景预定义检验族校正时显著性可翻转。标准库组合数反例：N=100、M=10、n=3、k=2，p=0.0257884972；预定义100个符合门槛term且另99个零命中，完整BH q=1；当前命中单term BH q=0.0257884972，错误显示显著。该例99个term可具有重叠背景成员，不需互斥。
- **可信度**：实测-AST + 精确组合数，高。真实 tested 表达式给一个命中和一个零命中 term，实际只返回命中term；不依赖 scipy 或心算。
- **建议/验收**：在查询结果之外固定全部 eligible term 检验族，零命中也以p=1进入BH；若采用有效独立过滤必须给出设计/控制证据。验收上述对照不显著；N_Terms_Tested与预定义族一致，并与独立完整BH对拍。

### H06 · P1 / R-major：eggNOG `-` 占位符被当 term，破坏“可注释空间”修复

- **证据**：`enrich_analysis.py:80-83` 只 dropna 后按逗号/空格拆分，直接 `rows.append((gene,t,ns))`，无 `-` 排除/GO和KO格式校验；`:165-166` annot所有gene进入universe；`:170` `-` 也获得term背景集合。
- **触发链**：常见 emapper GOs/KEGG_ko='-' → parse_emapper → 伪term '-'；若一列有真实term、另一列为'-'也同时进入，占位term混合跨namespace；后续run_ora。
- **后果**：没有任何真实功能注释的基因也进入“可注释空间”，重现注释覆盖率差异引起的偏倚；'-'若规模5–500还被实际检验，可能输出“缺注释富集”。不是仅UI显示瑕疵。
- **可信度**：实测-AST，高。真实token循环输入g1:'-'、g2:'GO:1'输出 `[('g1','-','GO'),('g2','GO:1','GO')]`。
- **建议/验收**：明确过滤'-'/空值与末尾注释行；GO/KO至少词汇格式校验；按namespace定义可注释宇宙并记录覆盖率。验收双列'-'基因不计N_Annotated/N_BG；任何输出Term不为'-'，混合注释空间与手算一致。

### H07 · P1：HTML 把非显著 ORA 前8条展示成“显著 term”

- **证据**：`host_expression_report.py:239-246` 从每个 enrich_*.tsv 读取rows，只要非空就显示标题“显著term (Top8)”并输出rows[1:9]，不判 `Significant`/`Padj`。`enrich_analysis.py:259-262` 输出的是全部检验term，不是仅显著。
- **触发链**：ORA有合法term但无显著结果 → 生成包含Padj>0.05的文件 → report标题作显著宣称。
- **后果**：最终可分享HTML产生不存在的显著通路结论；README `:273`“报告自动附每集合显著term Top8”未兑现。
- **可信度**：实测-AST，高。执行真实report富集循环，输入 `NO_SIGNAL/Padj=0.9/Significant=False`，仍出现在“DE显著term”标题下。
- **建议/验收**：严格筛选Significant为真或Padj<声明阈值，无显著时明示0；可另给“全部term”。验收q=.9不出现在显著段，q=.01出现；汇总显著数与显示行一致。

## 2. 科研设计与解释问题

### H08 · R-major：run独立性、BioSample/BioProject与组织/物种混杂未成为分析契约

- **证据**：`derive_virus_labels.py:208-217` 用Run/Sample合并，仅保留Tissue/ScientificName/CollectionDate，丢BioSample/BioProject/LibraryStrategy等；`run_deseq2.py:96-105` metadata只含group，设计仅`~group`。`build_coexp_matrices.py:81-106` 门槛数的是TPM列/run；不统计唯一BioSample，不限制单物种，不匹配组织/项目。`viral_load_scan.py:192-209` 仅组织内扫描；同组织可仍混项目/宿主/批次。
- **触发链**：公共队列同一BioSample多run，或阳性多来自一个项目/组织/宿主而阴性多来自另一个 → labels → DE/coexp/scan；数目和比率门槛可能全部PASS。
- **科学后果**：技术重复被当生物重复，30/30不代表30独立个体；项目/组织/品种表达差被叫病毒响应。PCA标准化/k-means不能凭无批次变量识别消除这种混杂；DE组织分层也不能解决组织内项目混杂。README `:223` 把组织平衡理念映射到人数闸门+另分支分层，不等价于主共表达网络组织平衡。
- **可信度**：静态高可信，真实队列影响为科学条件风险；未宣称本项目确有重复/多物种。
- **建议/验收**：保留BioSample/BioProject/宿主/测序策略；技术重复先按明确规则合并；网络按独立样本计门槛，单宿主分层、组织和项目匹配/残差化；DE设计加可识别协变量并检测完全混杂。验收提供run→BioSample映射、去重前后计数与组×Tissue×Project×Species表；完全混杂不得报告可辨识感染效应；leave-one-project-out或独立队列方向一致。

### H09 · R-major：CDS/转录本特征直接叫“基因”，未gene-level汇总/长度offset处理

- **证据**：`host_expression.py:275` gffread `-x`提CDS；`salmon_auto_pipeline.py:667-675` 以quant.sf Name原样为特征；`:755-759` NumReads原样汇总；`run_deseq2.py:96` round().astype(int).T，`:100-106`未传EffectiveLength/offset。主管线生成四矩阵但DE只消费NumReads `host_expression.py:691-694`。没有tx2gene/GFF Parent汇总步骤。
- **触发链**：参考含多转录本/重复CDS/共享序列，或样本间片段长度/异构体比例变 → transcript estimated counts → DE及模块 → 注释/ORA却使用“gene”概念。
- **科学后果**：结果是参考feature层面的检验，不必然gene层；多异构体同基因可被多计，EM共享read不确定性未传播，样本依赖有效长度会影响原始估计计数。转录本counts取整本身是公开设计，不应单独报必然错误；单基因单CDS且长度稳定时风险明显较低。
- **可信度**：静态高可信，影响依参考设计而定。
- **建议/验收**：要么明确“CDS/transcript-level”，导出唯一gene映射并证明一一对应；要么tx2gene汇总、适当使用长度校正的估计counts/offset，再一致映射DE/module/ORA。验收多异构体模拟同基因一行，gene总量稳定仅isoform switch时不产生伪gene-DE；输出 feature→transcript→gene表与歧义序列说明。

### H10 · R-major：未完成/未覆盖检测的run被补阴性；单病毒0组也包含其他病毒感染

- **证据**：`host_expression.py:335-342,396-398` 使用宿主QC名单，而非④检测成功名单；`derive_virus_labels.py:159-171` 无样本/病毒长表行就pos/load=0并Virus_any='-'；`:129` 缺覆盖也fillna(0)；`:109-128` 多行覆盖/载量取max。单病毒`host_expression.py:445-460` 只选 `<virus>_positive` 1/0，未排除 N_virus>0 的其他病毒。
- **触发链**：④仅处理部分run、检出表截断、coverage缺失 → 宿主QC仍通过 → 赋“阴性”；或目标病毒阴性但其他病毒阳性 → 被当目标比较的0组。
- **科学后果**：unknown与真正完成检测未检出混淆；缺失机制/深度可能与项目和表达相关，组差异偏倚。单病毒比较本来允许co-infection是可成立设计，问题是不能称“纯该病毒感染 vs 健康”或把目标特异效应当因果。分节段max最多证明某一参考/节段覆盖高，不必然全病毒感染证据。
- **可信度**：实测-AST标签行为高；队列中unknown/co-infection比例未量化。真实样本循环给RUN_NEVER_TESTED赋 `Virus_any='-', positive=0, load=0`。
- **建议/验收**：对齐④每run完成状态，三态positive/negative/unknown，未完成排除；区别覆盖不达阈值、未检出、未检测。报告单感染/共感染矩阵，加入其他病毒载量/状态或单感染敏感性分析，分节段标签写明证据标准。验收未测试run不进阴性组，检测成功零命中才进阴性；目标0组其他病毒比例显式展示。

### H11 · R-major：同数据DEG选module再做module-load，非独立剂量验证；样本域还不同

- **证据**：`host_expression.py:472-495,713-717` 默认用同队列all DEGs限定网络，并再传同TPM与virus_labels；`run_deseq2.py:224` 始终emit(all)，组织分层输出不是模块默认输入。`coexp_modules.py:63-82` DEG边选网络；`:139-152` 载量关联只取长表曾出现的sample；`:154-180` 对每组module均使用全selected samples，按病毒缺行再排除，不限定网络组、组织、项目，也不使用holdout；`:188-189` 对选后p做BH。
- **触发链**：同样本病毒二分类选DEGs → 单组网络学习模块 → 同样本载量与profile相关。载量与二分类高度关联，因此被选择的响应gene自然更易载量相关。模块Load使用长表交集且各病毒缺行删除，而全基因组 `viral_load_scan.py:151` 缺行补零；两翼统计域不同。
- **科学后果**：选择后p/FDR未涵盖DEG筛选与社区学习，不提供独立验证；阳性网络module的Group只是网络出处，不是关联样本组，容易误解。公共队列即使加DE组织分层，模块仍由汇总all DEG筛选并跨组织相关。缺行删除改变零载量/阴性样本组成，“互为校验”不能自动成立。
- **可信度**：静态高可信；选择偏倚幅度未真实重采样量化。
- **建议/验收**：明确exploratory/post-selection association；独立BioProject/BioSample留出训练DEG+模块，固定成员后holdout关联；或全流程嵌套置换纳入选择。明确 all/组内/组织内估计目标和完整检测样本域，两翼输出n及样本名单。验收heldout数据不参与任何模块选择；调试例中载量N与预注册domain一致；关联p不能冒称独立验证。

### H12 · R-major：TPM秩相关“不受变换影响/与计数等价”宣称不成立；病毒Load也非绝对剂量

- **证据**：`README.md:151-154`、`doc/METHODS_TEMPLATE.md:194` 写“TPM（秩相关下等价）”；载量scan `viral_load_scan.py:61-81,147-161` 用每gene跨样本TPM；模块 `coexp_modules.py:161-167` 先各gene跨样本z-score再均值。标签 `derive_virus_labels.py:46,118,127` 优先Asm_TPM并取max；README `:19-22`写“是否渐进、可归因于载量”。
- **原理/反例**：Spearman只对同一变量所有观测统一的单调变换不变；TPM每个样本的归一化分母不同，并非统一变换。固定特征长度，counts=[10,20,30]，样本length-normalized总量=[100,1000,10000]，TPM=[100000,20000,3000]；与载量秩[1,2,3]的rho由+1变-1（标准库实算）。gene z-score后相加也不与原TPM均值/特征均值秩等价。
- **科学后果**：宿主整体表达组成/高丰度转录本变化即可制造TPM负相关；Asm_TPM是另一分母下的相对丰度而非可直接解释的绝对感染剂量。相关、组织内相关或跨分支一致不排除反向因果/项目混杂，不能据此归因病毒载量。
- **可信度**：数学反例实测，高；真实composition bias规模未知。
- **建议/验收**：删“秩相关下等价”和因果措辞，改“相对丰度的探索性关联”。注明Load分母与检测样本域、区分abs/relative load；比较适当归一化counts/VST或稳健log-ratio/残差与TPM，必要时外源标定/宿主分母载量。验收上述反例纳入文档/单测，真实数据提供跨归一化/项目分层方向与排名敏感性；模块谱正文应写标准化gene TPM均值而非简单TPM均值。

## 3. 文档与代码对拍

| 文档承诺 | 当前代码 | 判断 |
|---|---|---|
| README `:48-50` 切virus覆盖05产物 | 顶层/子步存在缓存；force仍resume，H01 | 不一致，必须修复后才可承诺 |
| PROVENANCE `:36-44` 按名字集合修resume | TEA_mp `:131`、fast `:90`仍offset | 不一致；原始实际AST探针确证 |
| README `:46-47,136` 每组30/ratio3，all不过跳过 | build_matrices `:98-106`，host `:516,585-588` | 一致，是run数量门槛，不是独立样本/组织平衡 |
| README `:48-49` virus只切coexp，DE不变 | host `:445-460` vs `:696` | 一致，有意设计，不报参数传递bug；病毒特异解释须受限 |
| help `:945-950` DE与Load共用strata-col；README `:187`示例 | DE传de_strata本身，Load才strata_col | H04，不一致；`--de-strata Tissue`可用 |
| README `:145-150`分位非校准、Louvain无bootstrap | diff分位；module固定seed/res；HTML `:169-174,197-199`确有警示 | 已披露且报告保留，不重报未知算法bug |
| README `:153`/Methods `:194` TPM秩等价 | 样本特异分母与标准化均值，H12 | 数学错误，与是否采用TPM设计无关 |
| README `:273`只展示显著富集term | HTML `:239-246`显示全部前8 | H07，不一致 |
| README `:111`/报告称“显著对” | 有明确非校准警示，但标题/字段仍significant | R-minor：可改“尾部候选对”，不能当统计显著 |
| README `:224`GO/KEGG/MapMan未纳入 vs §9有ORA | §9/实际enrich实现已有 | R-minor：旧映射表应说明未纳入的是注释生成/专用流程，不是ORA缺失 |
| README `:152`模块谱TPM均值 | module `:161-167`gene逐行z-score均值 | R-minor：明确谱定义与参与样本域 |

## 4. 科研改进（不与已知限制重复计缺陷）

1. **先提供设计表而非先加统计复杂度**：唯一样本/项目/宿主/组织的交叉计数、感染状态与测序深度/策略分布；若完全混杂，明确不能估计感染效应。主网络的组织平衡不能由另分支DE分层替代。
2. **网络QC是诊断不是验证**：`coexp_network_qc.py:90-96` Status仅取Max_z，mean/share虽记录但不参与gate；PASS并不证明网络整体可比、完整或有充分功效。`host_expression.py:981-984`默认warn确有文档声明，不报未知bug；投稿运行建议error/skip并保留失败组，不能因为某一个极值超过1.7就写“Δz可信”。同时报告selected k、每分区n、组内有效相关数/缺失边比例。
3. **分位/模块稳定性**：已知分位非检验、单res无bootstrap应继续保留。尾部候选数可作工程指标，不作生物学显著计数。可先做组大小匹配子采样/多个resolution排序一致性，真正推断须纳入整网重建的零模型。
4. **小样本载量p值和检验族**：`viral_load_scan.py:99,65-81`最少n=8、t近似；module最少n=8且zero/ties常见。数值与scipy同式不等于p值在该样本规模已校准。报告零载量/ ties比例，适当置换；scan当前每病毒×stratum各BH（`:162`），不等于全项目跨病毒/组织整体FDR，应预定义检验族，模块FDR同样不包含模型学习。
5. **真正的复现清单**：索引hash、参考ID和gene映射、每分支精确sample名单、输入指纹、算法及依赖版本、selected k、Louvain/分区seed、组织/项目分层定义；生成HTML不应仅依赖目录里“看见什么就展示什么”。无模块/无关联重跑时 `coexp_modules.py:230-266,279-285` 不写新的组表/assoc，旧文件可残留；应为当前run写空状态而不是让富集扫描历史module表。
6. **参考定位边界**：`host_expression.py:164-170`先取all.genome.uniq.fasta，`:182-185`另取第一个GFF，没有组装accession/contig集合配对验收。当前未证实真实①产物有错配，故列为需要真实参考manifest核验的条件风险，不单列P0。多宿主项目应一宿主一参考一分析根。

## 5. 值得保留的设计与撤回/不报项目

- clean未去宿主reads与宿主表达目标一致；FASTA支持明确，不以质量值剥除本身报定量bug。
- optional ⑧不改主链；边界manifest为良好留痕基础，但路径不等于内容指纹。
- `--virus`只影响coexp、DE默认Virus_any有明确说明，**不是未知接线错误**；需要写清各证据分支的估计目标。
- 阳/阴共同基因集文件核验、Target交集对齐、自边按ID剔除、退化Δz退出码3都是合理防护。缺失边丢弃需统计数量，但不把已声明取交集误报为行序bug。
- QC缺meta为fail；QC失败样本从下游可通过名单排除，门槛可调；默认warn网络QC是公开策略，风险应以投稿闸门/状态解释处理。
- Salmon有样本失败会写failed_samples并 `salmon_auto_pipeline.py:1295-1300`退出2，主管线check=True阻断，**不报“失败后总是退出0”**。
- Louvain种子42、k-means种子42明确；未bootstrap已披露，非未知可复现性漏洞。
- 并行resume漏早期基因已实测，但当前rank通常查缺gene时报错；**撤回“必然静默完成不完整网络”说法**。
- 未实际安装/运行pydeseq2 0.4/0.5，不把环境依赖缺失或API兼容猜测当已证bug。

## 6. 修复与投稿优先级

| 优先档 | 条目 | 验收核心 |
|---|---|---|
| 必须先修、可直接翻结论 | H01/H02/H05/H06/H07 | 输入/参数变化不复用旧网；index-reference绑定；完整检验族BH；占位符排除；非显著term不作显著展示 |
| 流程可靠性与文档兑现 | H03/H04 | 乱序断点只补缺失且验证完整；显式分层真正执行，缺列不静默降级 |
| 投稿设计必须交代/独立验证 | H08/H09/H10/H11/H12 | 独立生物样本与项目混杂处理；feature层级/长度校正；检测unknown三态；选择后推断与样本域；相对丰度非因果解释 |

## 7. 实测记录、覆盖限制与可复核方式

- 工具：环境实际托管Python报告 `Python 3.13.14`，所有命令用 `python -B -c`，仅stdlib（ast/pathlib/types/io/math/re/argparse/typing等）。33/33 AST parse与 `compile(source, path, 'exec')` 内存编译通过；没有导入分析模块，也没有pyc/test/cache写入。
- 实测结果：`ACTUAL_RESUME_TASKS [(2,'C')]`；`force=True`真实网络命令仍 `--resume`；新virus旧缓存执行0次；旧index在新reference请求仍返回；缺B目标触发真实rank worker `KeyError('B')`；`--strata enabled`而非Batch；ORA hit-only族并精确BH翻转反例；emapper'-'变GO term；报告展示q=.9/False为“显著term”；未测试run自动阴性；TPM反例rho+1→-1。
- 一次 `_coexp_substep`抽取探针最初因受控命名空间缺`Sequence`出现NameError；补入typing所需绑定后通过。该失败是**审查探针环境问题，非管线缺陷**，没有拿它作缺陷证据。
- 新报告可由逐条 `文件:行号`与上述实际AST机制复核。未创建额外探针文件，避免违反仅报告写入要求。所有实际验证为内存表达式/语句与只读桩，不是假称运行了依赖齐全的真实生信链。
- 覆盖：主编排、Salmon建index/quant/汇总、标签、样本闸门、Δz、网络QC、Louvain模块、DE、载量scan、ORA、HTML均精读；TEA并行/分区/rank重点追踪。全部33py只做语法广度覆盖，未宣称每个相关系数底层和跨宿主保守性算法完成同等深审。
- 未执行全tests（其可产生临时文件、cache和第三方依赖要求），未联网重新核上游commit及此前39k文件md5声明，未测真实Salmon/DESeq2/网络性能、真实组数/ID碰撞/技术重复/多宿主频率。既有`__pycache__`在开始ls时已存在，本审查未增删。本文不是对真实结果已全面失效的结论，而是对当前可达实现与投稿解释的证据分级审查。
