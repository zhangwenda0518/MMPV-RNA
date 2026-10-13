# MMPV-RNA 全管线全面审查报告（2026-10-13）

> 审查对象：工作树现状（含 226 文件未提交重构：⑥⑧归档、根级 GUI 移除、AI 层未跟踪入库）。
> 方法：7 路并行只读审查（①公共数据+②清洗 / ③发现前半 / ③发现后半 / ④分析 / ⑤进化 / ⑦提交+AI层 / 横切文档配置产物）+ 全量测试 + 两个常态审计器实跑。
> 纪律：本轮未修改任何仓库文件；所有行号基于 2026-10-13 工作树。

---

## 一、总体结论

**测试 254 passed / 0 failed；主链默认路径（standard 布局、graph 关闭、全部依赖就位）未发现必然崩溃的 P0；但发现 4 个 P0 级问题（全部集中在 ⑤ 系统发育管线，其中 1 个是本次未提交改动引入的回归）和 21 个 P1。** 最大的结构性风险是：**工作区里 `virome_phylo_pipeline/phylo_pipeline.py` 相对 HEAD 丢失了已提交的退出码/blocked/Incomplete 记账整套逻辑**（详见 P0-1），意味着当前工作树不能直接提交；其次是 ③ 编排器"失败仍 exit 0 + 失败仍写 .ok"的一族问题，会让批处理把失败运行当成功、resume 永久跳过失败阶段。

数量统计：P0 ×4，P1 ×21，P2 ×40+（正文全列），审计器误报澄清 ×2 组，验证通过的关键契约 ×30 项。

---

## 二、P0（会产错结果 / 相对已提交版本回归）

### P0-1 【⑤·未提交回归】工作区丢失 HEAD 的「受阻 stage 记账 + 退出码 + Incomplete」整套语义
- `virome_phylo_pipeline/phylo_pipeline.py`（工作区 md5 `f562c44e…`）对比 `git show HEAD`：`pipeline_exit_code()`、五态 `_skip` 记账（blocked）、`SystemExit` 三分支退出码、CSV 的 `Incomplete/Stages_blocked` 列**全部被删**。`st = {"done","skipped","failed"}` 三态（:4188），`_skip` 只打日志（:4287-4289），main() 三个模式末尾无任何 `raise SystemExit`。
- 后果：⑤ 所有模式退出码恒 0（`&&` 链/批处理把失败当成功）；geo/geo_paths/gene_dating 缺席不再记账；HTML PARTIAL/Incomplete/Stages_blocked 全部消失；`utils/report_builder.py:374-375,427` 读 `stages_blocked` 成无写入方的死代码。RUN_GUIDE.md:298-301 仍按旧契约承诺验收。
- 建议：从 `git show HEAD:virome_phylo_pipeline/phylo_pipeline.py` 取回该段合入新 stage（treetime/phymapr），补单测；修好前 RUN_GUIDE 该节标注"待恢复"。
- **这是本轮最高优先级：决定工作树能不能提交。**

### P0-2 【⑤·AUDIT-1009 P0-1 未修】`align` 阶段 batch 模式从不执行 MAFFT
- `phylo_pipeline.py:1955-1987`：batch 模式 `prep["combined_fasta"]` 已被 clean 就地改成未比对的 `clean.fasta`（:4443），`_usable_alignment` 只判"非空+含>`"（:367-386）→ 把未比对序列当 provided alignment 早退返回；`run_mafft` 调用点 :1984 不可达。`STAGE_REFERENCE.md:311-312` 仍承诺 `mafft --auto`。
- 后果：rdp5/tree/clock/beast/gene_dating/splitstree 全链吃**未比对**输入。
- 建议：给 prep 增加 `provided_alignment` 标志（仅 `--virus`/work_dir 的 `--alignment` 允许早退），batch 必须落到 :1984。

### P0-3 【⑤·AUDIT-1009 P0-3 未修】publication_figures 收到未比对 fasta + 空 tree，进程整体失败
- `phylo_pipeline.py:5894-5896`：`--rtt_fasta` 传已被改写的 clean.fasta、`--rtt_tree` 传不存在的键（恒 `""`）。`publication_figures.py:111` `AlignIO.read` 对未比对输入抛 `ValueError`，`build_all_publication_outputs`(:423) 与 CLI(:587-620) 均无 try → Fig2/Fig3/Table S1 全部缺失，`run_postprocess` 只留一行 WARNING。
- 建议：改传 `results["outputs"]["alignment"]/["tree"]`；或让 `plot_temprest` 捕获异常按 `_declare` 语义登记 None。

### P0-4 【⑤·AUDIT-1009 P0-2 未修】`clean` 以「metadata 无此行」为由静默删除整批序列
- `utils/seq_clean.py:358-381`：`row is None` → 无条件计入 issues → REMOVE；fail-loud 护栏（:288-302）不覆盖"FASTA id 与 metadata 对不上"；`n_keep==0` 时 `process_virus` 记 `_skip('clean')` 并**降级用原始序列**（:4449），无任何标记。
- 建议：仅当 `drop_no_date/-location` 打开时才把"无此行"计 issues；收敛 `_norm_header` 与 `seq_ids.base_sample_id`。

---

## 三、P1（静默错误 / 契约断裂 / 真实崩溃）

### ① 公共数据 + ② 清洗
| # | 位置 | 问题 | 建议 |
|---|---|---|---|
| P1-1 | `public_metadata_pipeline/preprocess_unified.py:263-278` | 调 host_depletion 未传 `--logs_dir` → 统计产物落 CWD `host_depletion_logs/`，② report 阶段找不到 → `After_Kraken2/After_HostAlign/Retained/Aligner_Rate` 全空 + 污染启动目录 | stage_deplete 补 `--logs_dir {deplete}/logs`（一行修） |
| P1-2 | `data_preprocessing_pipeline/preprocess_report.py:120-128` | 按 file 列去重对 host_depletion 的临时随机路径无效 → 失败重跑样本的 After_Kraken2/After_HostAlign **双倍累加**（append 模式 + temp 路径每次不同） | 按 (sample,stage) 取最后一批而非按 file 求和 |
| P1-3 | `data_preprocessing_pipeline/clean-data.py:140-151` | 组内 >2 文件两个静默数据问题：多 lane 输入整对丢弃（`X_L001_R1/R2` + `X_L002_R1/R2` 只留 L001）；只有 R1 到货时同样本名拆多任务并发写同一输出（jobs>1 交错损坏）。内存复现确认，系旧问题非 f3286ed 引入 | 同 sample 多文件显式告警/任务号区分输出名 |

### ③ 发现管线（前半）
| # | 位置 | 问题 | 建议 |
|---|---|---|---|
| P1-4 | `virome_discovery_pipeline/virome_pipeline.py:3913-3925` | **编排器失败恒 exit 0**：阶段失败只 warning（默认容错模式），与 ⑤ 三档退出码契约、① 失败即非零的口径都不一致 | `if failed_stages: sys.exit(1)`（或分档 1/2/0） |
| P1-5 | `virome_pipeline.py:750-753,1077-1083` + `virus_identification.py:1339-1357` + `cobra_pipeline.py:1005-1033,1122-1260` | ident/cobra 子脚本对逐任务失败恒 exit 0 → 编排器 `.ok` 照写 → `--resume` **永久跳过失败阶段**（feb352f 只堵了"空壳目录"，没堵"带日志的失败目录"，`_stage_validated` 判据是目录非空） | cobra/ident 全部任务失败时退出非零；`_stage_validated` 改要求真实产物（`*cobra*.fa` 等） |
| P1-6 | `graph_propagation.py:390-399` ↔ 未提交 `filter_virus.py` | graph demote 读 `cdd/cdd_calls.tsv`，但未提交改动已删 filter 的 CDD 层 → 文件不存在 → demote 集合恒空、审计表 cdd_call 列恒空；单测自造文件所以仍绿 | graph 回退读 `{sv}/uniprot/uniprot.unclassified.id`（数据已存在）；提交 CDD 删除时同步 graph+单测 |
| P1-7 | `filter_virus.py:381-385` ↔ `virome_pipeline.py:990-994` | 逐样本 filter（单样本项目/`--no-merge`）时 COBRA 镜像落 `<d_filter>/<S>/candidates`，编排器查顶层 `<d_filter>/candidates` 落空 → 回退吃 02a 原始候选，**"filter 重新成为主链闸门"在该模式不成立**（graph 开启时恰好掩盖） | `build_cobra_mirror` 用 filter 根定位；补断链测试 |
| P1-8 | `virome_pipeline.py:929-935,946-958` + `filter_virus.py:59-62` | 逐样本模式编排器断点文件名错配（查从不产出的 `{sample}.virus.candidate_filtered.fasta`）+ filter 无 `--force` → `--force` 无法强制重算 UniProt 层，切 filter-mode 复用旧产物 | 断点判据改 `{sample}_virus.all.candidate_filtered.fasta`；filter 加 --force 并透传 |

### ③ 发现管线（后半）
| # | 位置 | 问题 | 建议 |
|---|---|---|---|
| P1-9 | `virome_pipeline.py:1511-1513,1568-1570` | rescue 两个早退分支（无目标宿主 centroids / 无 vclust novel）在写 `HQ_plant_viruses.fasta`（:1661+）**之前** return → 产物直接消失，但 `.rescue.ok` 照写；9a(:2185)/9b(:2481)/报告(:837-845)/⑦ HQ 回退全链静默空。触发条件是真实输入类别（某宿主全部 known 或全部 CheckV≥90% 直通） | HQ 汇总提到早退之前；`_stage_validated['rescue']` 补 HQ 文件存在性 |

### ④ 分析管线
| # | 位置 | 问题 | 建议 |
|---|---|---|---|
| P1-10 | HEAD `generate_pipeline_report.py:2399` ↔ `auto_known_virus.py:931-933` | **HEAD 提交态**：编排器传 `--ai-api-key`，报告只认 `--ai-key` → 环境里有 API key 时 Stage 9 报告必崩（argparse exit 2）。**工作区未提交改动已修**（:2491 加别名）——需尽快提交固化 | 提交该别名修复；`--report_ai` 真值路径加测试 |
| P1-11 | `virus-full.py:846-847,826-828` + 241/269-271,1029,1083 | **de novo 组装全失败时静默用"参考基因组本体"顶替**（Step4 把参考当融合骨架、Step8 用 REFERENCE_BASE 回填、Step9/10 copy 兜底）→ 终产物=参考序列，且被 `check_assembly_quality` 判 100% 覆盖 → 写 `Fully-assembled.ok`，报告/漏斗计入"全长组装成功"，无任何字段区分 | Evolution_Stats 加 `origin` 列（de_novo/fast_track/ref_fallback）；f2/f3 无效禁止写 ok |
| P1-12 | `generate_pipeline_report.py:1192-1193`（未提交引入） | `plots_root` 单根选择：batch_virus_depth 总是创建 `plots/`，编排器却把 Sample Distribution/Co-abundance 图写到 `01_Detection/` 顶层 → 报告取 `plots/sample_distribution` 必然落空，4 张图恒缺 | 两个候选根都遍历 |

### ⑤ 系统发育（除 P0 外）
| # | 位置 | 问题 | 建议 |
|---|---|---|---|
| P1-13 | `gene_partition_dating.py:599-605` + `phylo_pipeline.py:2958-2961` | 缺 BEAST 时假成功（子进程 exit 0 + 外层按 returncode 判成功写 checkpoint）→ 永久跳过；"BEAST 不可用"的 stdout 被 capture 吞掉 | 早退 `sys.exit(3)` 或落 sentinel；成功判据改为"确有条链日志/XML" |
| P1-14 | `phylo_pipeline.py:4021,4830-4835` | 相对 HEAD 丢 `_beast_ran()`/日志行数闸门 → `--skip_beast_run` 假 done；近空 `merged.log` 解析出**假 TMRCA** 进汇总（AUDIT P1-4 的落地实现被删） | 恢复 `_beast_ran()` 与 <10 行拒绝复用闸门 |
| P1-15 | `phylo_pipeline.py:2694-2695,5759` + `gene_align_generator.py:115,122-148,258` | 硬编码 legacy 模块目录残留（capheine 的 `data/ncbi_ref`、work_dir 分支 `time/beast`、gene_split/参考 CDS）→ standard 默认布局下 GBK 参考找不到、Table S1/Fig2/3 面板 M-P 静默缺失；同文件两种口径并存 | 三处统一 `ph_dir()`；加"根目录+模块名字面量"静态检查 |
| P1-16 | `recombination_analysis.py:25,49,96,113` | `UTILS = Path(__file__).parent` 指向根而非 `utils/` → `rdp5_validate.py`/`rdp5_bdp_plot.py` 两个子步骤以"无法打开文件"退出码 2 **静默失效**（`--rdp5_stage all` 时）。这是全部断引用里少数真缺失之一 | `UTILS = Path(__file__).parent / "utils"`；调用前断言存在 |
| P1-17 | `utils/seq_clean.py:331-340,381` | `--strip-illegal` 无效：先无条件记 issue 再剥离、状态只看 issues → 开开关也被 REMOVE（AUDIT P1-1 未修） | 剥离成功不计 issue |

### ⑦ 提交 + AI 层
| # | 位置 | 问题 | 建议 |
|---|---|---|---|
| P1-18 | `virome_submission_pipeline/submission_pipeline.py:299` vs `:380,:1692` | `_gate_locate(workdir, explicit, names, rel_dirs)` 第二形参叫 `explicit`，但 `validate_submission_gate` 以 `src_path=` 关键字传入 → 必然 `TypeError`，调用点无 try → **⑦ 主流程跑到提交闸门即崩**（一行修） | 形参改名 `src_path` |
| P1-19 | `mmpv_common/copilot.py:161-173` + `advisor.py:198-211` | 命令草稿 spec 与 5 条管线真实 argparse **全面漂移**（submission 缺 `discovery` 子命令；preprocessing 7 个 flag 不存在；analysis `-t/-j/--resume` 不存在；public 8 个不存在；discovery 3 个不存在）→ plan_chain/agent 生成的步骤必然 validate 失败或 rc=2 被当 blocked | spec 加 `subcommand` 字段并从各 `--help` 重新生成 |
| P1-20 | `virome_analysis_pipeline/utils/tbl2gb.py:120-126` | 负链坐标**两端各缩 1bp**（`start-1` 后负链分支未重新归一：`(1179,1101)`→应为 `(1100,1180)`）→ 负链 CDS 全部错位（fdb417a 修复只做了一半）。另：顶层 `virome_analysis_pipeline/tbl2gb.py:122` 是未修复旧副本（无调用方但可被 import） | 负链分支 `FeatureLocation(min-1, max, strand=-1)`；删顶层副本或改转发 |
| P1-21 | `mmpv_common/copilot.py:223` | 解释器白名单 token 过宽：`pythonXxx`/`pythonevil` 等 basename 被当解释器剥掉 → 配合 agent.py 直 subprocess 执行 | 收紧为 `re.fullmatch(r"python\d*(\.\d+)?(\.exe)?", base)` |

### 横切
| # | 位置 | 问题 | 建议 |
|---|---|---|---|
| P1-22 | `report_pipeline.py:962,978,1012` + `scripts/audit/check_io_handoffs.py:100-106` | **未提交改动让常态审计闸门从 0 FAIL 变 3 FAIL**：stage 显示标签从"09b_分析验证"改成 `"09b_Analysis_Verify"`，与 legacy 目录名字面量相同 → 审计器字面判 FAIL（实为标签非路径，`_resolve_acv_dir` 双名解析正常）。doc/IO_LAYOUT_DESIGN.md:312、wiki/09:7、wiki/00:93 声称的 "0 FAIL" 已不成立 | 标签改中性名（如 `09b Verify`）或审计器加语境判定；同步三处文档 |
| P1-23 | `doc/MMPV_dependencies.md:125-144` | A12 表缺 6 类 ⑤ 实际调用的外部工具：`trimal`、`raxml-ng`、`lsd2`、`mrbayes`、`phylobayes`、`splitstree6-tools`+`Xvfb`（splitstree 还是正式 stage）→ 换机按 A12 部署必缺件 | A12 增行并标注可选/降级行为 |

---

## 四、P2（精选 24 条，全量见各管线分报告）

**配置与布局**
1. `pipeline_config.yaml:4-11` `io_layout:` 键是死配置（无代码读它；实际链是 CLI > env > standard），用户写 legacy 会被静默忽略。
2. `pipeline_config.yaml:77-81,149-160` `deplete:` 与 `downstream.cleanup/skip_*` 段无任何消费者（`--profile downstream` 退化为 default）。
3. `virome_pipeline.py:763-764` vs `:858-859`：`graph_enabled` 兜底默认一处 True 一处 False（当前 yaml 恒赋值故无实害）。
4. `mmpv_common/io_layout.py:240` resolve_layout docstring 写"env > legacy"与 `DEFAULT_LAYOUT="standard"`（:37）矛盾（②④⑤ 三路独立发现）。
5. `io_layout.py:222-231` 未定义 ⑦ 提交输出根（实际产 `submission_<run_title>/`，process_report 只能 glob 猜）。
6. `data_preprocessing.py:95-96` ensure_hostdb 硬编码 `host_reference`，与 `preprocess_unified.py:434` 布局感知落点（05_HostRef）不一致。
7. `preprocess_unified.py:451-455` `if not str(fasta_input)` 恒假 → "clean 无输出即终止"是死代码。

**③ 报告与状态**
8. `report_pipeline.py:962/978/1012`（除标签问题外）standard 布局下阶段状态表显示 legacy 名——显示层混杂。
9. `report_pipeline.py:604-613` cluster_skip 读 `<d_cluster>/centroids`，实际写方在 `4_centroids` → 块恒不执行（死代码）。
10. `report_pipeline.py:1141-1152` `_find_final_taxonomy` 只取第一份 `*.integrated`，多份 taxonomy 时报告只统计一份（89c7a77 已把"逐样本"定为契约，build_ref_info 已合并全部，报告侧没跟上）。
11. `virome_pipeline.py:441-451` `_stage_validated['taxonomy']` 的逐样本 glob `*final_...` 不跨目录 → 非 Votus 样本名时 .ok 每轮被删、R 共识每轮重跑。
12. `report_pipeline.py:1023-1052` 报告硬编码 `Plant` 宿主子目录，多宿主 `Plant_Fungi` 时 plant summary 静默缺列。
13. `report_pipeline.py:409-416` filter_summary 表 `filter/comb` 两行同一文件（CDD 移除后的陈旧语义）。
14. `virome_pipeline.py:1646-1652` rescue 硬失败只 warning 不 raise → `.rescue.ok` 照写（与本文件 750/893/1072 三处自订契约相悖）。
15. `discovery2analysis.py:78` completeness=NA 时桥直接 traceback（官方入口无防御）。
16. `virus_classifier.py:411-433,878-884` PhaGCN3 无法被选中也不进合并表，但 --help/文档列为可用（"8 工具"实为 7+1）。
17. `virus_classifier.py:451,464` 文本 IO 无 encoding → Windows cp936 下可能 UnicodeError（④ 侧已全量 utf-8，此处未同步）。
18. `virome_pipeline.py:1982-1988` Rscript 调用 `VAR=值 cmd` bash 专有语法 + 路径未引号。
19. `filter_virus.py:204 等` `contig.split('_',1)[0]` 样本名含下划线时拆分目录名截断 → graph/COBRA 查找全部错配。
20. `virus_identification.py:439-445,1101-1103` 多组装工具并发同样本写同一目录（竞态）；串行模式无 per-sample 容错；鉴定空结果不清理旧候选（stale 产物被下游消费）。
21. `pipeline_config.yaml:118-121` 承诺 `--stage graph` 但 CLI 无此选择项（graph 只能改 YAML 开启）。
22. `virome_pipeline.py:117-150` shell=True 无 executable + 63 处 f-string 未 quote → 路径含空格即断。

**④ 分析与报告**
23. Stage2 过滤表把"过滤前统计"标成 Post-Filter（`filter_summary.py:164-165` 在过滤前写 summary）；过滤 0 通过仍写 run.ok 且 `get_summary()` 静默回落 best 层 → 下游悄悄跑在 best 层；`auto_known_virus.py:527` `args.filter` 默认 True 使任何单阶段运行附带跑 filter。
24. `extract_full_fasta.py:257-261` `--no_dna_fill` 死参数 + fill 模式 min_len/max_n 失效；accession 正则漏单字母前缀（`Potato_virus_Y_U13908.2` 提取失败）→ N-fill 静默跳过 + 报告统计丢条目；节段病毒 `Accession` 逗号拼接 → 报告 π/Shannon 永不命中；相似性面板硬编码 `Mode_Filter`、失败仍 mark ok；`batch_virus_depth.py:786-807` 空批次落盘 + 残留高编号批次混入陈旧样本；`batch_virus_full.py:119` 硬编码 `python3` + 超时分支死代码不清孙进程；`virus-full.py:942-945` 缺 rnaviralspades 静默跳过；`batch_virema_dvg.py:570` 缺 Rscript 阶段崩溃、全失败仍 exit 0；报告对 legacy 存量项目只认 standard 目录名 + `--report_ai` 死开关；方法学章节仍写 capheine/"10 stages"（实为 9 阶段）；`virus_metadata_plot.py:590` `--legacy` 引用不存在的 legacy 脚本；根目录与 utils/ 两套 `extract_full_fasta.py` 并存（文档只写裸文件名，易改错文件）；viral_maftools.R 运行时联网装包（离线集群挂）。

**⑤ 进化（除 P0/P1）**：beast checkpoint 三口径不一（beast1.xml vs merged.log vs mcc.tree）；README/utils 计数 51 vs 实际 65、两个根目录脚本列在 utils 下；treetime/phymapr 文档承诺 blocked 实代码 `_skip('na')`；splitstree bridge 缺失被记硬失败（文档说可选）；`_unweight.py` 死回退分支；原生布局固定端口 8001 + pgrep 复用任意同名进程；stage state 非原子写；popgen_r_backend.R 的 nh/Hd 是占位假值（已自我披露）；geo 坐标粒度/来源信息消费侧丢失（fuzzy→'resolved'、coord_source 有写无读）。

**⑦ 提交 + AI**
25. Core14 死键：`unified_metadata.py:345,402,964` 读 `age_stage`/`taxid`，但 load_metadata 的 lookup 不含这两个键 → `src-Dev_stage` 与 `samp_taxon_id` **永远取不到值**（Core14 有 Age_GrowthStage/TaxID 列，两行修）。
26. `load_metadata:183-187` 用 CenterName（机构）推断采样地 `'China'` —— 机构≠采样地的臆测残留。
27. `report_html.py:234-236` suvtk features CDS 计数恒 0（`split()[-1:] == ['CDS']` 判不出行尾 CDS）。
28. `copilot.py:105-129` COMMAND_SPECS 两个 `"hostref"` 键 → build_host_pipeline 规格被覆盖成死数据。
29. MCP `run_digest` 有一个 `write=True` 写盘口（与"只建议不执行"红线字面冲突，建议默认关闭）；`requires_confirmation` 把"校验失败"也算"需确认"（语义误导）；LLM 输出未转义直入报告 HTML（ai.py:459）；`plan_with_ai` 把用户自然语言原样写入 ai_log.jsonl（可能随 provenance 发布）；缓存键不含 api_key。
30. 管线内 `submission_gui.py` 列定义是第二副本（与 unified_metadata 逐字段一致但会二次漂移），建议改 import。

**横切文档与产物**
31. "必需 25 / 可选 25"实为 **24/26**（6 处文档：README、DB_DEPLOYMENT、wiki×3、Home）；DATABASE_SETUP.md §7 CheckV DB 写 ~3GB 实为 ~6.5GB、MEGAN map（~20GB）全仓 0 引用。
32. SOFTWARE_VERSIONS.txt 停在 2026-08-26（仅 7 行 pixi global 工具 vs 现约 114 依赖）；pixi.toml 与 report_pipeline.py:29 仍写 **v2.3**（报告标题输出 v2.3）。
33. wiki 三处矛盾：Home.md:46 "legacy 默认"（实为 standard 默认）；"467 脚本" vs 08-脚本参数参考 "198 核心"；"0 FAIL/26 WARN" vs 现跑 3 FAIL/23 WARN。
34. `doc/14-pipeline-diagrams.md`（v2.3，legacy 名）与 `MMPV_LOGIC_FRAMEWORK.md`（旧①-⑤编号）无"已过时"横幅仍挂在索引。
35. `check_script_refs.py` 解析器缺陷确认：不支持链式 `SCRIPT_DIR/"utils"/"x.py"`、`os.path.dirname(abspath(__file__))`（58 处）、`.parents[N]`、缩进赋值（`^` 锚定）→ 43 对断引用中 **~36 对误报（≈84%）**；真/待确认 7 对：build_virus_db 的 3 个 db_build 外部脚本（文档已披露待收编）、biosoft/VirPhyKit 3 处文件名不符（`Insert_node_numbers.py` 实为 `NodeNumbers_Insert.py`）、`virus_metadata_plot_legacy.py`（有 .exists() 保护，--legacy 必报错）。
36. `outputs/MMPV_8PIPELINE_REVIEW_20261010/` 总结索引 8 条链接全部断链（相对路径指向不存在位置）。
37. `.gitignore` 未覆盖 `_build/`、`_trash_20261010/`（19MB）、`.tmp_probe/`、`outputs/` → git status 噪音 + 误提交风险。
38. 旧产物 `reports/METAGENOMIC_out_run_20261003/pipeline_report.html`：9 个内部锚点悬空；版本 v2.3 **现在重跑也不会变**（report_pipeline.py:29 写死）；无乱码。
39. 10-10 八管线审查"改变结论清单"①统计量（phylo Tajima's D/Hd/Fst、FEL/MEME q、analysis Pi、host ORA BH）与 ③ submission/EVE 项**未被后续修复声称覆盖**——仍开放。
40. `AUDIT_REPORT_20261009.md`（⑤ 自查）逐条核对：P0-1~4、P1-1~5 **全部未修**（本轮 P0-2/3/4 即其 P0-1/3/2），P1-4 还发生回退（闸门实现被删）。

---

## 五、误报澄清（审计器口径修正，避免误修）

1. **check_io_handoffs 3 条 FAIL**：`report_pipeline.py:962/978/1012` 的 `"09b_Analysis_Verify"` 是 `_add()` 的 stage 显示标签，不参与路径拼接；真实目录走 `_resolve_acv_dir` → `_D['d_verify']`（standard=`12_AnalysisVerify`/legacy=`09b_Analysis_Verify` 双名回退），已实测两布局解析正确。→ 属审计器把标签当路径的误报（但因未提交改动才触发，见 P1-22）。
2. **check_script_refs 54 条断引用**：43 对去重后 ~36 对为误报（解析器不支持链式路径/`.resolve()` 盲区/外部命令名），7 对真/待确认（见 P2-35）。`virome_analysis.py→analyze_hypothetical/viral_topology`、`batch_virema_dvg→DI-tector_06.py`、`generate_pipeline_report→capheine_pipeline.py`（已迁⑤存在）等均为误报。
3. ③ 编排器 5 处 legacy 名 warn（:1005/1090/2207/3058/3097）逐行确认全部是日志/help/docstring 文本，非路径字面量。

---

## 六、验证通过的关键契约（实测核过）

- **测试**：254 passed / 0 failed（mmpv_common 8 件 + ①②③ tests，10.44s）；23 个在审 .py `ast.parse` 全绿（累计 36+ 文件）。
- **AI 层红线主体成立**：`mcp_server._DISPATCH` 无任何执行工具（有测试断言）；API key 只以掩码形式出现，ai_log 只记 prompt 哈希/usage；无 key 降级完整；agent.py 是唯一 subprocess 点且 token 白名单 + 危险命令强制跳过 + `--set` 不可注入（P1-21 是唯一口子）。
- **报告模板层**：Jinja2 autoescape+StrictUndefined fail-loud；mmpv_common 8 个测试与实现 AST 级零漂移；④ 报告未迁移模板（f-string + `_esc`，不存在"StrictUndefined 在旧数据上崩"的问题）。
- **io_layout 注册表**：standard/legacy 键集完全一致；35 个调用方全部经 import，**0 处硬编码布局目录字面量**；③→④ 桥 `ref_info.tsv` 17 列表头与 doc 逐字一致。
- **Core14 列同步**：① `gsa_sra.info.py FINAL_14`（14 列）= ⑦ `unified_metadata` 文档列 = 报告/GUI 消费列，无 14↔13 漂移（P2-25 的死键是 lookup 构造问题，非列定义漂移）。
- **② 配对修复 f3286ed**：`GS-2_S2_L001_R1/R2`、`QH-1`、`GS-2_1/_2`、单端 `241-1.fq.gz` 内存实测全部正确；`test_clean_pairing.py` 与代码同步。
- **③ graph（开启时）闭环**：输入/输出/救回降级逻辑完整；`graph.enabled=false` 时 run_filter 与 8fea067 之前逐行等价（唯一非等价是 COBRA 候选源改镜像，系 10-07 裁定的有意变更）。
- **prescreen/virus_profile 契约**：stat_prescreen DEC_COLS / stat_virus_profile / ① 报告三方列名一致；`stage_funnel_summary.py` 属 ④（只吃 ④ 五阶段产物），与 prescreen 无契约关系（任务前提不成立，无对齐问题）。
- **panvirome 分母**已收敛到 `read_sample_total()` 单点，全仓无第二处 prevalence 计算；**stage_funnel** 五阶段输入对上真实产物（含 `.tmp_demo_funnel` 实跑产物核对）。
- **fdb417a 报告合并**：被删 379 行是旧精简版，保留版是超集，无功能丢失。
- **③→④ 桥**：`locate_centroids/locate_taxonomy` 三代候选名兜底实测命中；台账落 `90_Handoff/handoff_manifest.tsv`。
- **⑦ sequin 修复**：mol_type 受控词表、12 列对齐、作者两步拆分、非法 CDS 行切断，方向全部正确；table2asn 缺失降级链完整（`.sqn >1000B` 成功判据）。
- **GUI 移除**：全仓 grep 确认零活跃死引用（仅历史文档提及）；pixi.toml 归档 task 已注释；check_io_handoffs/check_script_refs 的 SCOPES 已同步移除 ⑥⑧/GUI。
- **⑤ 20-stage 体系自洽**：STAGE_ORDER/DEPS/GROUPS/CHECKPOINT_FILES 程序化核对无缺项；未跟踪新 utils（phymapr/treetime/splitstree_native_viz/unweight/popgen_plots/overview_fig/sci_palette）**全部已接线且签名匹配**；"坐标不臆造"在 geo 全链贯彻（miss→None→跳过/null+标注）。
- **capheine 迁移**：④ 侧调用已清理，⑤ 入口输入契约与 FIX_REPORT 一致。

---

## 七、建议处理顺序（醒来后按此修）

1. **先决断（提交前必做）**：⑤ 工作区 phylo_pipeline.py 与 HEAD 的 diff——恢复 P0-1 丢失的退出码/记账逻辑（连同 P1-14 的 `_beast_ran` 闸门），否则整棵未提交改动不能整体提交。
2. **一行修三连**：⑦ `_gate_locate` 形参（P1-18）；① `--logs_dir` 透传（P1-19 前的 P1-1）；④ 提交 `--ai-key/--ai-api-key` 别名（P1-10）。
3. **③ .ok 语义批修**：ident/cobra 全失败退出非零 + `_stage_validated` 收紧（P1-5）+ rescue HQ 早退（P1-9）+ 编排器 `sys.exit(1)`（P1-4）——这一族决定"跑完算不算成功"的机器可读性。
4. **④ virus-full ref_fallback**：加 origin 列 + 禁止零证据写 Fully-assembled.ok（P1-11，直接关系论文数据可信度）。
5. **跨模块两处**：graph demote 数据源切 uniprot（P1-6，随 CDD 删除提交一起）；逐样本 filter 镜像落点（P1-7）。
6. **⑤ align/MAFFT**（P0-2）与 publication_figures 传参（P0-3）、clean 元数据护栏（P0-4）、gene_dating 假成功（P1-13）、UTILS 路径（P1-16）、硬编码目录三处（P1-15）。
7. **AI 层**：copilot spec 重新生成（P1-19）、解释器白名单收紧（P1-21）、MCP write 口默认关、LLM 输出转义、Core14 死键两行修（P2-25）。
8. **审计器与文档批量**：check_io_handoffs 标签白名单、check_script_refs 解析器补链式/`.resolve()`/dirname（消掉 84% 误报）；24/26 计数、v2.3 版本号、wiki 三处矛盾、A12 六类工具、10-10 索引断链、.gitignore 四目录。

---

## 八、审查覆盖与方法说明

- 7 路并行子代理共读约 90 个核心文件全文 + 300+ 处定点核对；两个常态审计器实跑；`mmpv_db.py list` 实跑（50 项清单核对）；`layout_dirs`/`build_discovery_dirs` 两布局实测；配对逻辑与正则边界用纯内存复现（未落盘）。
- 各路分报告（含 P2 全量与逐条证据）由 7 个审查会话产出；本文件为汇总裁决版。误报判定均给出判定依据（审计器解析器源码 + 目标文件实际存在性 + 字符串语境）。
- 未覆盖：未实跑任何管线阶段（无 241/246 环境与外部工具）；⑤ 无正式 tests/（唯一测试实体在 reports/_phylo_fix_20261011/tests，且针对旧快照版本）； biosoft/ 第三方代码只处理被点名条目；archive/ 未审（按定义）。
