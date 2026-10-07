# ⑧ host_expression_pipeline — 宿主表达分析（可选补充模块）

> **定位：可选补充模块，不在 MMPV ①→⑤ 主链上。**
> 抽取自 goji-virome 差异表达分析项目（枸杞病毒宿主响应）的核心链路，以零侵入
> 方式并入 MMPV：不改动任何既有管线的行为与产物；输入全部取自主链现有产物；
> 输出落独立根 `08_HostExpression/`（standard 布局）。不跑它，主链完全不受影响。
>
> 编号说明：管线号为 ⑧（wiki 的 ①-⑦ 编号中 ⑦ 是提交管线）；输出根 `08_HostExpression`
> 遵循"根号 = 管线号"惯例（07 号空缺 —— ⑦ 提交管线消费 ④ 产物、无独立输出根）。

主链回答"有哪些病毒"，本模块回答"**宿主如何响应病毒感染**"，分三个层面：

| 层面 | 阶段 | 方法 | 定位 |
|---|---|---|---|
| **基因间协调关系**（主分析） | 05_Coexp | TEA-GCN 阳/阴双网络 → ΔzScore(MR) → z-QC 闸门 → DEG 限定模块 | 核心分析（goji 正稿主图） |
| 基因表达量 | 06_DE | DESeq2 阳性 vs 阴性（可选组织分层） | 交叉参考 |
| 基因剂量响应 | 08_LoadScan | 病毒载量 × 全基因组 Spearman + BH-FDR | 剂量验证（默认含于 `all`，`--no-with-loadscan` 可关） |

三条分支共享上游、互为校验：二分差异（coexp/DE）回答"有没有差"，剂量扫描回答
"是否渐进、可归因于载量"。goji 的方法学结论：**组织平衡与剂量建模防的是不同的
伪差**——公共队列上汇总口径受组织构成混杂影响（ERAD 教训），结论以组织内分层
与剂量响应的交叉验证为准。

---

## 1. 快速开始

```bash
# 全流程 (输入自动定位; 参考由 ① 的 genome+GFF3 现场提取)
python host_expression_pipeline/host_expression.py --stage all \
    --project-root <项目根> -t 60 -j 8

# 预览将执行的命令 (不落盘、不建目录)
python host_expression_pipeline/host_expression.py --stage all --dry-run \
    --project-root <项目根>

# 单阶段重跑 (断点续传; --force 强制)
python host_expression_pipeline/host_expression.py --stage coexp \
    --output_dir <项目>/08_HostExpression --force
```

**一核两翼定位**：核心 = 05_Coexp（TEA-GCN Δz → z-QC → DEG 限定模块）；
06_DE 是核心的输入提供者（模块检测的 DEG 表 + 方向拆分富集）兼交叉参考；
08_LoadScan 为剂量验证翼（默认跑，`--no-with-loadscan` 可关）。
公共队列（样本量充足）建议追加：`--de-strata Tissue`（DE 组织分层）。
coexp 闸门默认每组 ≥30 样本且两组比 ≤3×（`--coexp-min-group/--coexp-max-imbalance`
可调；`--stage all` 下闸门不过会跳过 coexp 让后续阶段继续，显式跑则报错）。
指定病毒单独分析：`--virus GCVA` → coexp 分组切为 `GCVA_positive` 1 vs 0
（默认 Virus_any 阳性 vs 阴性；06_DE 不受影响；切换口径会覆盖 05_Coexp 现有产物，
需保留时先归档）。

输出根结构（standard 布局 `08_HostExpression/`；两布局同名）：

```
08_HostExpression/
├── 01_HostRef/     host_cds.fasta              宿主 CDS (genome+GFF3 → gffread -x)
├── 02_Quant/       quant/<样本>/quant.sf        salmon 定量
│                   matrix/{TPM,NumReads,FPKM,EffectiveLength}_matrix.tsv
├── 03_QC/          qc_summary.tsv  qc_passed_samples.txt
├── 04_Labels/      virus_labels.tsv  sample_labels.tsv  virus_columns.tsv
├── 05_Coexp/       matrices/{positive,negative}.tsv + group_summary.tsv (闸门记录)
│                   teagcn_{positive,negative}/ (分区/RAW_GCN/Completed_GCN)
│                   network_qc.tsv (网络 z 分布 QC, Fig S8 闸门)
│                   diff_coexp_{candidates,significant,gene_aggregate,summary}.tsv
│                   modules/{module_table_<组>.tsv, module_summary.tsv,
│                             module_load_assoc.tsv} (DEG 限定模块, Fig5 流程)
├── 06_DE/          deseq2_all_results.tsv  deseq2_degs.tsv  volcano_plot.png
│                   deseq2_strata/<组织>.tsv (--de-strata 时)
├── 08_LoadScan/    <病毒>_genome_scan.tsv  <病毒>__<组织>_genome_scan.tsv
│                   load_scan_summary.tsv (stratum 列区分 all/组织)
├── 07_Enrich/      enrich_summary.tsv  enrich_<集合>.tsv (DE/dz_hub/module_*)
│                   background_genes.txt  sets/<集合>.txt
├── 09_Report/      host_expression_report.html
├── .checkpoints/   断点标记 (含 coexp 7 个子步)
└── host_expression.log
```

## 2. 阶段与脚本

| 阶段 | 目录 | 脚本 | 说明 |
|---|---|---|---|
| hostref | `01_HostRef/` | 编排器内联（gffread） | genome+GFF3 → CDS；`--reference` 直给则跳过 |
| quant | `02_Quant/` | `salmon_auto_pipeline.py` | 自动 SE/PE、自动建索引、并行 quant、汇总 4 矩阵 |
| qc | `03_QC/` | `quant_qc.py` | 比对率/reads 双阈值，剔除 16S 等垃圾样本 |
| labels | `04_Labels/` | `derive_virus_labels.py` | ④ 检出表 → 阳性/阴性 + 连续载量；合并 Core14 Tissue |
| **coexp** | `05_Coexp/` | `build_coexp_matrices.py` → `teagcn/` 引擎 → `coexp_network_qc.py` → `diff_coexp.py` → `coexp_modules.py` | 见 §3 |
| de | `06_DE/` | `run_deseq2.py` | pydeseq2（兼容 0.4/0.5+ 两代 API）+ 可选 `--de-strata` 组织分层 |
| loadscan | `08_LoadScan/` | `viral_load_scan.py` | 剂量验证层（默认跑，`--no-with-loadscan` 可关）：向量化 Spearman + BH-FDR，Tissue 存在时自动加做组织内扫描 |
| enrich | `07_Enrich/` | `enrich_analysis.py` | 功能富集 ORA（超几何+BH）；注释缺省时**自动跳过**，见 §9 |
| report | `09_Report/` | `host_expression_report.py` | 单页 HTML（各段缺省自动跳过） |

## 3. coexp：TEA-GCN 差异共表达（主分析）

数据流（`teagcn/` 为 vendored 上游引擎，见 `teagcn/PROVENANCE.md`）：

```
02_Quant TPM 矩阵 ─┐
03_QC 通过名单 ────┤→ build_coexp_matrices (按 04_Labels 分组切两个子矩阵 + 闸门)
04_Labels 分组 ────┘      │
                          ▼
        Generate_partitions ×2 (k-means 扫 k, silhouette 选 k)
                          ▼
        Run_TEA-GCN_mp ×2 (PCC/SCC/bicor → nanmax → RAvg; fork 多进程)
                          ▼
        Rank_transform ×2 (MR = √(rank_ij·rank_ji) → zScore, 越大越强)
                          ▼
        coexp_network_qc (每组网络: Max_z / z≥1.5·1.7 份额 → PASS/LOW/FAIL;
                          FAIL ⇒ Δz 闸门 warn/skip/error, 默认 warn)
                          ▼
        diff_coexp (Δz = z(阳)−z(阴); 每基因 Top-N 邻居;
                    0.5%/99.5% 分位阈值 → 显著对 + hub 聚合)
                          ▼
        coexp_modules (z≥1.7 ∩ DEG 的边 → Louvain res=1.0 → 模块表/hub;
                       模块表达谱 × 病毒载量 Spearman + BH-FDR;
                       DEG 表在 06_DE 之后 → coexp 阶段先跳过, de 阶段自动补跑)
```

**KingdomStress 增量（2026-10，源自 bioRxiv 10.64898/2026.04.17.719161 方法移植）**：
① `coexp_network_qc.py` = 其 Fig S8 物种剔除规则的工程化——样本量不足会整体压低
MR z-score 分布（36 物种中 5 个 max_z<1.7 被剔除），Δz 是两个网络各自 z 标准化的
差值，任一侧分布截断都会让 Δz 被量纲差主导，故 Δz 前强制体检；② `coexp_modules.py`
= 其 Fig 5 模块检测流程（z≥1.7 且两端均为 DEG 的边 → Louvain resolution=1.0，
seed=42）；③ `coexp_module_conservation.py` = 其跨物种模块保守性分析（OrthoFinder
直系组映射，OC=共享OG/min≥0.5 且共享≥10 判保守，物种内基因→OG 置换 100 次做经验
p 值）——多宿主项目各跑一份 ⑧ 后独立调用，不在编排器链上。

**相对源项目脚本的改进**（`diff_coexp.py` 重写自 `diff_coexp_strict.py + _post.py`）：
全量行序校验（源版只抽查前 3 个基因）、自边按基因名精确剔除（源版固定 dz[0]=0）、
Target 输出真实基因名（源版 0 基索引曾致 7,547/3,908 口径误用）、Top-N/采样步长参数化。

**TEA-GCN 引擎补丁**（输出逐字节不变，已验证）：① 死代码跳过（full=False 时
从不写入的字符串拼接占单基因 ~80% 耗时 → 5× 提速）；② Rank_transform 的
Windows 路径修复。mp 入口（fork 多进程，Linux）与串行版在真数据上
39,224/39,224 文件 md5 一致（goji 验证）。

**样本量要求（方法内禀）**：每组 ≥30、两组比 ≤3×（闸门强制）；且分区样本数
n/k 决定相关系数功效——goji 为 239/19≈12.6。样本不足时请用 06_DE +
08_LoadScan（闸门信息里也会给出此建议）。

**磁盘预估**：双网络 RAW+Completed ≈ `2×70B×基因数²`（39k 基因 ≈ 200 GB）；
编排器会在 >50GB 时警告。Rank_transform 已由 `rank_transform_mp.py` fork 并行
（输出与串行逐字节一致，39k 基因 5.8h→~25min；共享数组内存 ≈ 12 GB）。

**已知限制（审稿人视角，如实声明）**：
- Δz 的 0.5%/99.5% 分位阈值是**排序/截断装置而非校准的检验**——候选对为每基因
  Top-N 极值，相对全总体分位天然偏松（最大值统计效应），显著对数量不宜作生物学
  计数；结论请基于 hub 排序与下游模块/富集。原则性替代：置换零模型（组内打乱
  样本标签重建网络的经验分布），计算昂贵，暂未实现。
- Louvain 模块为单分辨率（res=1.0, seed=42），未做 bootstrap 稳定性；可用
  `--module-resolution` 扫描核验跨分辨率一致性。
- 模块×载量关联以 TPM 模块均值谱做 Spearman（秩相关，对单调变换不敏感）；
  组织混杂防护依赖 QC 后分组构成，公共队列建议配合 `--de-strata` 交叉验证。
- DE 用 pydeseq2 原始计数模型；载量扫描与模块谱用 TPM（秩相关下等价），
  两类检验不混用计数模型。

## 4. 输入边界（跨管线交接，B8）

| 来源 | standard 路径 | 本模块参数 |
|---|---|---|
| ② clean reads（**未去宿主**） | `02_Preprocessing/01_CleanData/{3.clumpify,2.fasta}` | `--reads_dir` |
| ① 宿主基因组 + GFF3 | `01_PublicData/05_HostRef/genome/` | `--genome` / `--gff` |
| ① Core14 元数据（Tissue 列） | `01_PublicData/02_Unified/Global_Unified_Metadata_Core14.{csv,tsv}` | `--metadata` |
| ④ 病毒检出表 | `04_Analysis/01_Detection/summary/all_viruses.best.summary.tsv` | `--detect` |

以上经 `--project-root` 自动定位（standard 优先、legacy 兜底；目录名统一走
`mmpv_common/io_layout.py`）。每次运行的消费/产出登记到
`90_Handoff/handoff_manifest.tsv`（boundary=`hostexpr`）。

**为什么用 ② 的 clean reads**：分析的是宿主基因表达，需要宿主 reads；② clean
已做接头修整且与 ④ 同为"未去宿主"输入语义。**FASTA 说明**：② 产物为 FASTA
（质量值已剥除），salmon 不使用质量值，定量不受影响；引擎已扩展识别 `.fa/.fna/
.fasta`（含 .gz）。

## 5. 通用性（不止病毒项目）

模块的引擎层（定量/QC/coexp/DE/报告）是**物种与分组无关**的；病毒特异性只集中在
两处——04_Labels 的检出表推导与 08_LoadScan 的载量扫描（缺 virus_labels 时自动
跳过）。三种典型用法：

**① MMPV 病毒组学模式（默认）**：`--project-root` 自动定位 ①②④ 产物，
分组 = `Virus_any`（+/-），`--group-col <病毒>_positive` 可切单病毒对比。

**② 通用二分组模式（非病毒/任意比较）**：自备分组表（`Sample` + 分组列，可含
任意元数据列），跳过 ④ 检出表；分组/分层列全部可换：

```bash
python host_expression_pipeline/host_expression.py --stage coexp de report     --output_dir out/08_HostExpression     --reads_dir clean_reads/ --reference host_cds.fasta     --labels-file my_labels.tsv     --group-col Condition --ref-level control --alt-level treated     --pos-level treated --neg-level control     --de-strata --strata-col Batch     --coexp-min-group 10
```

适用：处理 vs 对照、病/健组织、品种对比、任意物种（宿主参考由 `--reference`
或 genome+GFF3 提供）。TEA-GCN 引擎本身物种无关（吃任意 TPM 矩阵）。

**③ 连续变量场景**：分组无意义时跳过 coexp/DE，仅用 08_LoadScan 需病毒载量；
纯"表达×性状"关联属未来扩展（源项目范式，接口已留：sample_labels 任意列）。

闸门与默认值的来源说明：`--coexp-min-group 30`、`--coexp-max-imbalance 3.0`、
`--min-mapped 50` 等是从 goji 大队列经验固化的**保守默认**，全部 CLI 可调；
小样本/失衡场景的降级路径（06_DE + 08_LoadScan）在闸门信息中自动给出。

## 6. 依赖

| 组件 | 用途 | 获取 |
|---|---|---|
| salmon ≥ 1.10 / gffread ≥ 0.12 | 定量 / CDS 提取 | pixi 已含 |
| pandas / numpy / scipy / matplotlib | 各分析阶段 | pixi 已含 |
| pydeseq2 | 06_DE | pixi 已含（`[dependencies]`） |
| **einsumt** | 05_Coexp（TEA-GCN bicor 预计算） | pixi 已含（conda 无包，走 `[pypi-dependencies]`） |

teagcn 引擎与主环境版本范围兼容（上游 pin 的 numpy==2.0.1 等均在 MMPV 范围内）；
在非上游 pin 环境重跑可能与 goji 原结果存在浮点级差异（rank 并列翻转），首次
真机运行建议做一次串行 vs mp 的 md5 冒烟对比（`--teagcn-entry serial`）。

## 7. 与源项目（goji-virome 差异表达分析）的对应

| 源项目组件 | 本模块处置 |
|---|---|
| `salmon_auto_pipeline.py` | **原样抽取**（仅加 FASTA 后缀 + 出处注明） |
| 病毒标签（Covered%≥90%） | 重写为 `derive_virus_labels.py`（列名嗅探/分节段聚合） |
| TEA-GCN 四步 + 多进程加速版 | **vendor 进 `teagcn/`**（上游逐字节一致 + 2 个输出不变补丁 + 你们的 mp/fast 入口） |
| `diff_coexp_strict.py` + `_post.py` | 重写为 `diff_coexp.py`（见 §3 改进清单） |
| `deseq2_cross.py` / 精简工作流 `04_deseq2_de.py` | 合并重写为 `run_deseq2.py`（+ `--strata` 组织分层） |
| `viral_load_scan.py` | 重写为向量化版（+ 自动组织内扫描、stratum 汇总列） |
| `tissue_balanced.py`（组织平衡共表达检验） | 理念吸收进 coexp 闸门 + DE 分层 + LoadScan 组织内扫描 |
| GO/KEGG/MapMan 富集 | **未纳入**（需宿主功能注释源：eggNOG DB 或 Mercator4） |
| ggNetView 网络可视化 / NLP 分区注释 | **未纳入**（见源项目目录与上游 Colab） |

## 7b. KingdomStress 论文增量对应（2026-10）

| 论文组件（bioRxiv 10.64898/2026.04.17.719161） | 本模块处置 |
|---|---|
| Fig S8 网络 z 分布检查（36 物种剔 5） | `coexp_network_qc.py`（qc 子步 + `--coexp-qc-gate` Δz 闸门） |
| Fig 5 模块检测（z≥1.7 ∩ DEG → Louvain res=1.0） | `coexp_modules.py`（modules 子步；de 阶段自动补跑） |
| Fig 5 跨物种模块保守性（OC + 置换零模型） | `coexp_module_conservation.py`（独立脚本，多宿主项目用） |
| 模块×胁迫关联思路 | 模块表达谱 × 病毒载量 Spearman（`module_load_assoc.tsv`） |
| Fig 4 通路共现网络（小样本替代路线） | **未纳入**（可基于 workflow/05_mapman_enrich.py 扩展） |
| SD2/SD3 公开数据（13.6M DE 调用 / 36 物种阈值网络） | 外部验证资源（含番茄），未内置 |

## 8. 测试

```bash
python host_expression_pipeline/tests/test_host_expression.py
# 或: python -m pytest host_expression_pipeline/tests/ -q
```

41 用例：标签推导、QC、BH-FDR 已知值、Spearman 与 scipy 逐值一致、**coexp 切分
与闸门**（小样本/失衡/skip 模式）、**Δz 逻辑**（自边剔除/分位阈值/行序校验拒绝）、
**网络 z-QC**（PASS/LOW/FAIL 判级、error 退出码）、**模块检测**（双植入团精确回收、
非 DEG 排除、FAIL 组跳过、载量关联信号）、**跨宿主保守性**（OC/置换零模型/未映射
基因）、**编排器 modules 子步接线**（真实子进程 + 断点双条件）、salmon FASTA 识别、
输入定位、HTML 报告降级、**功能富集**（emapper/通用注释解析、植入 term 显著、诱饵不显著、基因集收集）。端到端验证（合成数据 + 真 TEA-GCN 引擎）：40+40 样本植入
"仅阳性组共表达"的模块 → Δz 找回其内部显著对 **12/0 全部增强**（方向性完美回收）；
DE 分层、LoadScan 组织内自动扫描、续跑/断点、台账同步验证通过。

## 9. 功能富集（07_Enrich，可选）

对 ⑧ 自动产出的基因集（06_DE 显著 DEG、Δz hub Top-N、每个共表达模块）做
超几何过表达分析 + BH-FDR；背景 = TPM 矩阵全部基因（自动写入）。**注释由用户
事先准备**（物种无关，跑一次即可），两种格式自动嗅探：

```bash
# ① eggNOG-mapper (推荐; 对 01_HostRef 的 CDS 翻译/或直接蛋白序列跑一次)
emapper.py -i host_pep.fasta -m diamond --cpu 60 -o host  # 产出 *.emapper.annotations
# ② 或 Mercator4 (MapMan bin) / 任意注释: 通用映射 TSV, 必含 gene,term 列
#    可选 namespace,name 列 (报告显示 term 名称)

python host_expression.py --stage all --project-root <项目根> \
    --annotation host.emapper.annotations      # 不给则富集阶段自动跳过
```

注意：注释文件的基因 ID 必须与 TPM 矩阵基因 ID 同名（同一 GFF 提取的 CDS
天然一致）；eggNOG 走 GO + KEGG_KO 两类 term，通用表任意词汇；`--enrich-hub-n`
控制 Δz hub 取前 N（默认 50）。产物 `enrich_<集合>.tsv` 含全部检验 term
（Significant 列过滤），报告自动附每集合显著 term Top 8。
