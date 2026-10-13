# MMPV-RNA 全管线修复报告（2026-10-13）

> 对应审查报告：`outputs/MMPV_FULL_REVIEW_20261013/REVIEW_SUMMARY_20261013.md`（同日 7 路只读审查的裁决版）。
> 本轮修复：P0 ×4 全修、P1 ×23 全修、P2 修复 50+ 项；少数 P2 因需产品决策/超范围/第三方目录而跳过（见第四节，逐条含原因）。
> 修复方式：6 路并行修复代理（①②/③编排/③filter+graph/③报告/④/⑦+AI，文件域互斥）+ ⑤ 系统发育管线与横切文档由主会话亲手修复。

---

## 一、验收指标（修复后实测）

| 检查项 | 结果 |
|---|---|
| 全量测试（mmpv_common + ①②③） | **282 passed / 0 failed**（修复前 254 passed；净增 28 个回归用例） |
| `mmpv_common/tests`（AI 层） | 164 passed（新增 15 条 spec 漂移/子命令/别名/白名单/缓存盐回归） |
| `scripts/audit/check_io_handoffs.py` | **PASS — 0 FAIL / 23 WARN**（回到基线；修复中途曾引入 2 FAIL，已改由 io_layout 注册表解析候选名后消除） |
| `scripts/check_script_refs.py` | 214 处引用 / 54 告警（与修复前同量级；**~84% 为审计器自身解析器盲区误报**，逐条判定见审查报告第五节。真缺失 7 对中：rdp5 两脚本路径已修，其余为已披露的外部收编缺口/第三方目录内部命名） |
| `mmpv_db.py list` | 50 项（必需 24 / 可选 26），与全部文档计数已对齐 |
| 改动文件语法 | 所有 touched .py `ast.parse` 全绿；`viral_maftools.R` Rscript parse 通过 |
| ⑤ 退出码机制 | 冒烟实测：failed→1 / blocked→2 / 干净→0；CSV `Incomplete/Stages_blocked` 两列恢复；`count_beast_log_rows` 闸门恢复 |

---

## 二、P0 修复（4/4，全部 ⑤ 系统发育管线）

| # | 问题 | 修法 | 位置 |
|---|---|---|---|
| P0-1 | 工作区相对 HEAD 丢失整套「退出码/blocked/Incomplete」记账（提交阻断项） | 从 HEAD 恢复并适配新代码：`pipeline_exit_code()`、五态 `_skip`（checkpoint/na/pending/degraded/blocked）、blocked 兜底记账、`stages_blocked/stages_skip_reasons/incomplete`、受阻与 pending 双警告、CSV `Incomplete/Stages_blocked` 两列、work_dir 模式独立出口 `SystemExit(_wcode)`、main 末尾三分支退出码 | `virome_phylo_pipeline/phylo_pipeline.py`（10 处调用点同步升级：metadata/clean→degraded，host→na/blocked 分流，beast/phylogeo→pending，gene_dating/geo/geo_paths→blocked） |
| P0-2 | `align` 阶段 batch 模式把未比对 clean.fasta 当"已有比对"早退，MAFFT 从不执行 | 新增 `prep["provided_alignment"]` 标志：仅 work-dir 模式（`--alignment` 语义）允许早退；batch/datasets 模式必落到 MAFFT | 同上 :1955-1985、work_dir prep :5876 |
| P0-3 | publication_figures 收到未比对 fasta + 空 tree → AlignIO 崩、Fig2/3/TableS1 全灭 | batch 模式按病毒名从记账结果取 `outputs["alignment"]/["tree"]` 真实产物（回退旧取值）；datasets 报告模式优先已产出的 mafft.aln.fasta | 同上 :6021-6054、:5839-5842 |
| P0-4 | `clean` 以"metadata 无此行"静默整批 REMOVE，上层降级用原始序列无人知晓 | 仅当 drop_no_date/-location 开启时才判缺行；未开门槛仅告警保留 | `utils/seq_clean.py:356-372` |

## 三、P1 修复（23/23）

**①公共数据 + ②清洗（代理 A，7 项中 3 项 P1）**
- P1-1 `preprocess_unified.py:273` 补 `--logs_dir`，去宿主统计不再落 CWD、report 阶段汇总表不再全空（monkeypatch 冒烟验证命令串）。
- P1-2 `preprocess_report.py` 改按 (sample,stage) 取"最后一批"，重跑样本不再双倍累加；新增 2 个回归用例 + CLI 端到端实测（1900/900/45% vs 修复前 3700/1900/95%）。
- P1-3 `clean-data.py` 多 lane 感知配对（跨 lane 样本键加 `_Lxxx` 后缀），歧义时 ValueError fail-loud；新增 5 个用例，241 命名（GS-2 等）无回归。

**③ 发现管线（代理 B1/B2/C，8 项 P1）**
- P1-4 编排器 `failed_stages` 非空 → `sys.exit(1)`（:3996-4002）。
- P1-5a/5b identification 全样本失败 / cobra 全任务失败退出非零；P1-5c `_stage_validated` filter/cobra/rescue 判据全部收紧到真实产物级（纯日志空壳不再通过，rescue 要求 HQ 文件存在）。
- P1-6 graph demote 双源回退：cdd_calls.tsv 缺失时读 uniprot unclassified 清单；新增 3 用例（含 cdd 优先级判别）。
- P1-7 COBRA 镜像统一契约：per-sample 模式也落 `<d_filter>/candidates/<S>/`；测试断言同步 + 新增 merged 回归用例 + 端到端核验。
- P1-8 filter 新增 `--force`（编排器 per-sample 分支透传）+ 编排器断点文件名改真实产物名。
- P1-9 rescue HQ 汇总提取为 `_write_hq_summary()`，三个早退 return 前均调用，HQ 产物任何路径不丢。
- P1-22 报告器 3 处 `"09b_Analysis_Verify"` 标签 → `"09b Analysis Verify"`（`09b` 前缀保留、`_stage_key_of` 映射实测不变），审计器回到 0 FAIL。

**④ 分析管线（代理 D，3 项 P1）**
- P1-10 `--ai-key/--ai-api-key` 别名（工作区既有修复，保留并纳入提交建议）。
- P1-11 virus-full：`origin` 列（de_novo/fast_track/ref_fallback）写入 Evolution_Stats；de novo 无有效产物**不写 Fully-assembled.ok**；batch 侧 ref_fallback ⚠️ 且不计成功数。
- P1-12 报告 Stage1 图区双根查找（`detect_dir` 与 `detect_dir/plots`），Sample Distribution/Co-abundance 不再恒缺。
- P1-20 `utils/tbl2gb.py` 负链改 `FeatureLocation(min-1, max+1, strand=-1)`（0-based 半开，比审查建议的字面公式更准确——代理以生产方 known_virus_submit.py 口径 + 合成 tbl 实测 1179/1101→79bp 验证）；顶层 tbl2gb.py 改 importlib 转发薄壳。

**⑤ 系统发育（主会话修，5 项 P1）**
- P1-13 gene_dating 缺 BEAST → `sys.exit(3)`（约定码=工具不可用）；wrapper 对 3 显式处理 + returncode 0 但目录无产物不写 checkpoint。
- P1-14 `_beast_ran()` + `<10 行` 日志闸门恢复（`count_beast_log_rows` 一并从 HEAD 取回——该函数也在工作区重写中丢失）；`--skip_beast_run` 不再假 done、近空日志不再产出假 TMRCA。
- P1-15 硬编码 legacy 模块目录 ×3 全换布局感知：capheine `data/ncbi_ref`→`ph_dir('data')`、work_dir beast→`ph_dir('time')`、`gene_align_generator.py` 新增带回退的 `_ph_dir()` 贯穿 7 处。
- P1-16 `recombination_analysis.py` `UTILS` 指向 utils/ —— `--rdp5_stage all` 的 validate/plot 子步骤恢复可用。
- P1-17 `--strip-illegal` 落实（剥离成功撤销 issue；剥空仍 REMOVE）。

**⑦ 提交 + AI 层（代理 E，4 项 P1）**
- P1-18 `_gate_locate` 形参 `explicit`→`src_path`（消除关键字调用隐患；注：原报 TypeError 未复现——内部实为位置传参，属预防性修复）。
- P1-19 copilot/advisor 命令 spec 全面重建：CommandSpec 支持 `subcommands/aliases/any_required`；7 条管线 spec 以 `--help` 实跑 + AST 提取核对，**幻影 flag=0、未登记=0**（phylo 仅最小修，见未修清单）；submission 草稿补 `discovery` 子命令。
- P1-21 解释器白名单 `re.fullmatch(r"python\d*(\.\d+)*(\.exe)?|py(\.exe)?")`，`pythonXxx` 类被拒。
- （P1-20 修复归 ④ 侧，tbl2gb 在其目录内）

**横切（主会话修，2 项 P1）**
- P1-22 见上（与 ③ 报告器同项）。
- P1-23 `doc/MMPV_dependencies.md` A12 补录 6 类 ⑤ 实际调用的外部工具（trimal/raxml-ng/lsd2/mrbayes/phylobayes/splitstree6-tools+Xvfb，含各自缺失时行为）+ VCF2Dis/VCF2PCACluster 提示。

## 四、P2 修复（50+ 项，按域归并）

- **①②**：ensure_hostdb 布局感知、`if not str(fasta_input)` 死代码、down 阶段失败回非零、比对率样本键 rsplit。
- **③ 编排**：graph_enabled 兜底统一 False、`--graph-enabled` CLI 开关新增、taxonomy 校验 glob 跨样本、rescue 全败 raise、`--tax_tools` help 摘除 PhaGCN3、ident 空 results 清 stale 候选 + 串行分支容错 + 同样本多组装工具任务合并（消除竞态）、（顺带）report_pipeline seqkit 同表口径问题已登记待议。
- **③ filter/graph**：死参数与死函数清理、docstring 现状化、`--sample-names` 最长前缀匹配防样本名截断。
- **③ 报告**：cluster_skip 三候选、`_find_final_taxonomy` 返回全部并合并去重（三调用方同步）、宿主子目录解析多宿主化、comb 行删除、桥 completeness NA 防御、virus_classifier 8 处 UTF-8、PhaGCN3 help 标注预留、virus_distribution 候选补 'centroids'、`VERSION v2.3→v3.1`。
- **④**：Stage2 Post-Filter 改读 passed 表（0 通过写表头空表、get_summary 回落告警、`args.filter` 门控修正）、extract `--no_dna_fill` 落实 + min_len/max_n 生效 + accession 正则放宽 4 形态 + ref_not_found 记账、S3 节段病毒多键、Mode_* 优先级 + similarity 全败不写 ok、批次空表不落盘 + 残留批次清理、`sys.executable` + TimeoutExpired 杀进程、rnaviralspades 缺失告警、DVG Rscript 预检 + run.partial、报告 legacy 三代候选（**验收时改为 io_layout 注册表解析**，消除验收抓到的 2 个新 FAIL）、`--report_ai` 删除、方法学 9 阶段口径、根级 `--legacy` 死分支删除、bed 防覆盖、R 脚本去联网安装。
- **⑦+AI**：Core14 死键 age_stage/taxid 补全、CenterName→China 臆测删除、CDS 计数等价重写（原"恒 0"前提经两份真实 .tbl 复核不成立，按更直白写法落地）、`hostref/hostref_build` 双键可匹配、MCP write 默认关 + digest 单次构建、`requires_confirmation` 语义修正、LLM 输出 html.escape、ai_log 只记白名单 flag 键、缓存键加盐、管线内 GUI 列定义改 import 单源。
- **横切**：pipeline_config（io_layout 死键说明、deplete/downstream 未接线标注、graph 开启口径）、.gitignore 收口（_build/_trash_*/.tmp_*，outputs/ 留作审查交付物显式保留）、pixi 版本 v2.3.0→v3.1.0、io_layout docstring 与 DEFAULT_LAYOUT 对齐、wiki Home（脚本数 198、⑤ 20 stage、布局默认 standard）、wiki 09/00/10 + README + DB_DEPLOYMENT 计数（必需 24/可选 26、0FAIL/23WARN、213 处）、10-10 审查索引 8 条断链修复、两篇历史文档加过时横幅、DATABASE_SETUP（CheckV 6.5GB、MEGAN 标注未引用）、SOFTWARE_VERSIONS 追记、⑤ STAGE_REFERENCE treetime/phymapr 失败语义与代码对齐（skip(na) 不阻断）、⑤ AUDIT_REPORT_20261009 顶部加逐条修复状态注记。

## 五、未修清单（逐条含原因）

| 项 | 原因 |
|---|---|
| ③ `virome_pipeline.py` 全文件 63 处命令串未 quote + run_cmd 无 executable | 改动面大且管线目标环境为 Linux 服务器；建议单独立项统一 `shlex.quote` |
| ⑤ AUDIT P1-3（clock `summary_rows` 非空即 mark）与 P1-5 | 涉及 RTT/DRT 门控的产品级口径选择，需作者拍板（已在 AUDIT_REPORT 注记中标"未修"） |
| ⑤ copilot spec 的 phylo 127 个 flag 仅登记 7 个 | 未登记只"误拒"草稿不误放行；⑤ CLI 正在演进，全登记易产生假告警（漂移测试会在变更时报警） |
| ⑦ io_layout 未登记 ⑦ 提交输出根（submission_root）+ stage_tools 提交后置条件借位 | 需注册表键位设计决策（涉及 run_digest 域），建议与下次 io_layout 版本一起做 |
| ④ 根级 `filter_summary.py`/`extract_full_fasta.py` 旧副本与 utils/ 漂移 | 建议按 tbl2gb 模式改转发壳；本轮只修了编排器实际调用的 utils/ 版 |
| ④ DNA `--max_n` 现为"填补前"判定 | 语义选择（填补前 vs 填补后），当前与过滤模式口径一致；需作者拍板再调 |
| `report_pipeline.py` seqkit 汇总 PE 取 max（非两端合计） | 与 ② 合计口径的取舍需作者拍板（② 代理范围外发现，已登记） |
| biosoft/VirPhyKit 内部 3 处文件名不符、db_build 3 个外部脚本待收编 | 第三方目录/已披露的外部收编缺口（README 与 db_build/README 已有指引），不动第三方 |
| 10-10 审查报告内列的统计量类遗留 | 属统计方法学修正，需逐项设计方案；**2026-10-13 复核后逐条状态如下表** |

## 六、统计量类遗留：二次核对 → 已全部修复（2026-10-13 当日）

先逐条对到当前工作树代码复核（当日第二轮），确认 5 条真实存在后，**按用户指示全部改为调用成熟库实现、不再手写公式**（引入 scikit-allel 1.3.13，已验证与 numpy 2.4.4 兼容；登记进 pixi.toml，`pixi install` 时会自动更新 lock）：

| 条目 | 修法（调用包） | 位置 | 数值验证 |
|---|---|---|---|
| ⑤ P0-3 Tajima's D 量纲混用 | 交由 `allel.stats.diversity.tajima_d`（内部为绝对 θπ−θW 标准式），输入 = ACGT 列的等位基因计数矩阵（pairwise-deletion 语义与 π/S 一致） | `popgen_analysis.py` `_diversity_core` + 新 helper `_ac_matrix` | 测试断言 == allel 参考值；且 fixture 上旧式分子与新式分子符号相反（旧实现会翻转解读方向） |
| ⑤ P0-5 位点法 Fst gap/标签错位 | 从「样本×位点×等位基因」计数矩阵聚合（gap 序列天然缺席，标签不可能错位），估计量换 **Hudson (1992)/Bhatia (2013)**（`allel.hudson_fst`，比值和口径）；全局 = 各群体两两聚合。`fst_pairwise` 热图同改 Hudson。选 Hudson 而非 Weir&Cockerham 是因为 allel 的 wc 实现硬性只支持二倍体，病毒是单倍体 | `popgen_analysis.py` `fst()`/`fst_pairwise()` | 测试：测试内**独立**重建组计数矩阵直调 allel 对拍一致；两两矩阵对称、对角 0；两组时全局==两两（自洽） |
| ⑤ P0-6 FEL/MEME q 不参与判定 | 阳性判定改 `q < fdr`（FEL 另加 β>α），`--fdr` 真实生效；BH 校正继续用仓库既有 `bh_fdr`（标准单调实现，10-10 未质疑） | `selection_analysis.py` `parse_fel_meme` | 测试：fdr=0.01 时 4×p=0.001 (q=0.005) 阳性、16×p=0.04 (q=0.04) 排除——旧实现全判阳性 |
| ④ V03 Pi=1−p_max | 交由 `allel.mean_pairwise_difference`（n/(n−1)·(1−Σp²)，样本量校正的逐位点成对差异）对有效覆盖位点取均值；N 读段按缺失；有效深度 <2 记 NaN。抽成模块级 `_site_pi_mean()` 便于无 pysam 环境单测 | `batch_virus_variants.py` | 测试：手算对照 0.5444 精确命中（旧 1−p_max 口径为 0.45，必不相等）；N/低深度位点剔除语义验证 |
| ④ V04 缺 MD 静默全零 | 两层：(a) 去掉 `with_seq=True`——代码只用比对位置不用参考碱基，MD 从此非必需；(b) 裸 `except: pass` → logger.warning + Pi/Shannon/ANI 以 NaN 初始化——"不可计算"与"真的是 0"在结果里可区分（聚合 np.average 如实传播 nan） | `batch_virus_variants.py` `worker_calc_metrics` | 结构性修复 + 语法/导入验证；`test_variant_stats.py` 覆盖 Pi 数值面 |

⑧ H05（ORA 检验族）仍随归档搁置；⑤ P0-4（Hd）此前已由 10-12 工作区改动修掉（本轮补了回归锁测试）。

**顺带排掉一颗历史地雷**：⑤ 测试与 ① 测试同会话时暴露 `gsa_sra.info.py:46` 的 `builtins.print = safe_print` 全局劫持——它让 `print.__module__` 变成 `gsa_sra_info`，同进程随后任何 numba 导入（scikit-allel 的依赖链）直接 `KeyError`。已改为模块级 `print = safe_print` 遮蔽（本脚本内所有函数经 globals 解析，防进度条撕裂的行为完全不变，不再污染进程）。

**结果数字会变的口径说明（论文侧需知）**：⑤ `tajima_D` 与 `fst_site`（含热图）的数值系统改用标准估计量，与旧版不可比；④ `Pi/Pi_avr` 从"非主等位占比"变为真正的核苷酸多样性（旧值系统性含义错误）；FEL/MEME 阳性位点集合改由 `--fdr` 控制（默认 0.1，因 BH 的 q≤p 性质，默认档比旧 p<0.05 略宽，需更严可传 `--fdr 0.05`）。

**新增测试**：`virome_phylo_pipeline/tests/test_popgen_stats.py`（4 用例）+ `virome_analysis_pipeline/tests/test_variant_stats.py`（3 用例，含 pysam 桩注入）；全量测试 **289 passed / 0 failed**，IO 审计器保持 PASS（0 FAIL / 23 WARN）。

## 七、提交建议

工作树现在包含：你此前的 226 文件重构 + 本轮 ~85 处修复。建议按两步提交（或合并为一次）：
1. ⑤ `phylo_pipeline.py` 的记账恢复（P0-1）与其余 ⑤ 修复单独成一个 commit（提交阻断项优先落库）；
2. 其余管线修复 + 文档/配置收口一个 commit。
未修清单第五节各项如需继续，建议按"统计量修正 → shell 引号化 → io_layout ⑦ 根键"顺序立项。
