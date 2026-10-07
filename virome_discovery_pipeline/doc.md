# virome_discovery_pipeline — 病毒组发现流程文档

## 数据流总览

```
Raw FASTQ
  │
[00a_CleanData]  ← data_preprocessing_pipeline/clean-data.py     (fastp → seqkit → clumpify)
  │
[00b_HostDepletion] ← data_preprocessing_pipeline/host_depletion.py (kraken2 → align → rrna)
  │
  │
[01_Assembly] ← assembly_pipeline.py (megahit / rnaviralspades / penguin)
  │
[02_Identification] ← virus_identification.py (10 工具并行)
  │
[02b_Filter/03_Filter] ← filter_virus.py (UniProt + CDD 两层过滤; 直通 COBRA 镜像)
  │
[03_Filter/graph] ← graph_propagation.py (配对读段链接图 + ensemble 种子 + Louvain 标签传播;
  │              对 filter 无证据保留结果二次裁决 (救回/降级); 产出 COBRA 输入镜像)
  │
[03_COBRA]  ← cobra_pipeline.py     (bwa-mem2 → coverm → cobra-meta;
  │            候选源 = d_filter/candidates 镜像, 恢复 filter→COBRA 主链契约)
  │
[04_CLUSTER] ← cluster_pipeline.py  (seqkit → cd-hit → vclust)
  │              ├── centroids/final_centroids.fasta
  │              ├── vclust_clusters.tsv
  │              └── split_fastas/
  │
[05_Taxonomy] ← virus_classifier.py + virus_classifier_analysis.R
  │              └── final_integrated_classification.tsv
  │
[06_HostPrediction] ← run_host_prediction.py
  │                    └── ensemble_host_summary.tsv
  │
[07_Checkv] ← checkv (直接调用)
  │
[08_Rescue] ← rescue_pipeline.py  (CheckV → Virseqimprover → BLASTN)
  │              └── final_centroids.fasta (HQ vOTUs)
  │
[09_Reports] ← report_pipeline.py (TSV + HTML + Sankey)
```

---

## 1. virome_pipeline.py — 主编排器

**用途:** 端到端宏病毒组全自动主控流水线，依次调用所有阶段脚本，管理数据流和检查点。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--input_reads` | — | 原始 FASTQ 目录 |
| `--output_dir` / `-o` | 必需 | 项目输出根目录 |
| `--stage` | `all` | 阶段选择: all,clean,deplete,assembly,identification,cobra,cluster,taxonomy,host,checkv,rescue,report |
| `--config` / `--profile` | `default` | YAML 配置文件和 profile |
| `--coassembly` | False | 合并所有样本 reads 进行共组装 |
| `--assembler` | `penguin` | 组装工具: megahit, rnaviralspades, penguin |
| `--aligner` | `bowtie2` | 宿主比对工具: bowtie2, hisat2, minimap2 |
| `--host-filter` | `Plant` | 宿主过滤类型 (Plant/Animal/Fungi/Bacteria) |
| `--threads` / `-t` | 20 | 线程数 |
| `--jobs` / `-j` | 2 | 并行任务数 |
| `--tax_jobs` | 1 | 分类并行数 |
| `--min-length` | 500 | 最小 contig 长度 |
| `--ani` | 0.95 | 聚类 ANI 阈值 |
| `--qcov` | 0.85 | 聚类覆盖度阈值 |
| `--checkv_threshold` | 90.0 | CheckV 完整性阈值 |
| `--skip_clean` / `--skip_depletion` / `--force` | — | 流程控制 |
| `--host_db` / `--virus_db` / `--checkv_db` 等 | — | 数据库路径 |

### 输入
原始 FASTQ 目录；各类参考数据库路径

### 输出
```
{output_dir}/
  00a_CleanData/         ← 清洗后数据
  00b_HostDepletion/     ← 去宿主后 reads
  01_Assembly/           ← {sample}/{sample}_{tool}.contig.fasta
  02_Identification/     ← {sample}_virus.all.candidate.fasta
  03_COBRA/              ← {sample}.{tool}.cobra.fa
  04_CLUSTER/centroids/  ← final_centroids.fasta
  05_Taxonomy/integrated/← final_integrated_classification.tsv
  06_HostPrediction/     ← ensemble_host_summary.tsv
  07_Checkv/             ← CheckV 质量报告
  08_Rescue/Plant/centroids/ ← HQ vOTU catalog
  09_Reports/            ← pipeline_report.html + TSV 汇总
```

---

## 2. clean-data.py — 数据清洗（`data_preprocessing_pipeline/`）

**用途:** RNA-seq 数据清洗: fastp (QC+去接头) → seqkit (FASTQ→FASTA) → clumpify (聚类重排，提升压缩率)。自动识别 PE/SE。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input` | 必需 | 输入 FASTQ 目录 |
| `-o` / `--output` | `./clean_out` | 输出目录 |
| `-j` / `--jobs` | 2 | 并行样本数 |
| `-t` / `--fastp-threads` | 4 | 单样本线程数 |
| `--dedup` | False | fastp 去重 |
| `--skip-clumpify` | False | 跳过 clumpify |
| `--clumpify-memory` | `10g` | Java 堆内存 |
| `--force` / `--dry-run` | — | 流程控制 |

### 输入
包含 `.fq/.fastq/.fq.gz/.fastq.gz` 的目录

### 输出
```
{output}/
  1.fastp_tmp/        ← 临时清洗 FASTQ
  2.fasta/            ← 转换后的 FASTA
  3.clumpify/         ← 聚类重排后 FASTA
  logs/               ← 日志 + fastp HTML/JSON
  .clean_checkpoints  ← 断点续传
```

**工具依赖:** fastp, seqkit, clumpify.sh

---

## 3. host_depletion.py — 宿主去除（`data_preprocessing_pipeline/`）

**用途:** 三阶段混合去宿主: Kraken2 分类 → 比对工具 (Bowtie2/HISAT2/Minimap2) 精细过滤 → Ribodetector rRNA 去除。支持 FASTA/FASTQ, PE/SE。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--tool` | 必需 | bowtie2 / hisat2 / minimap2 |
| `--seq-type` | 必需 | dna-short / rna-short / nanopore / pacbio |
| `-k` / `--kraken2_index` | 必需 | Kraken2 数据库路径 |
| `-x` / `--step2_index` | 必需 | 比对索引路径 |
| `-I` / `--input-dir` | — | 输入目录 (自动扫描) |
| `-O` / `--outdir` | — | 输出目录 |
| `--jobs` | 1 | 并行样本数 |
| `-t` / `--threads` | 4 | 线程数 |
| `--confidence` | 0.4 | Kraken2 置信度 |
| `--rrna` | False | 启用 rRNA 去除 |
| `--rrna_tool` | ribodetector | ribodetector / silva |
| `--steps` | kraken2,align,rrna | 步骤选择 |
| `-f` / `--filter` | true | true=移除宿主, false=提取宿主 |

### 输入
FASTA/FASTQ 目录；Kraken2 宿主库；Bowtie2/HISAT2/Minimap2 宿主索引

### 输出
- 去宿主后 reads (压缩 FASTA/FASTQ)
- Kraken2 分类报告
- seqkit stats 汇总 + 可视化图表
- 资源使用统计

**工具依赖:** kraken2, bowtie2/hisat2/minimap2, samtools, seqkit, ribodetector_cpu

---

## 4. BBNorm 覆盖度归一化（不在本管线，跨管线引用）

> 本管线（`virome_pipeline.py`）**没有 bbnorm 阶段**（`:3580` stage_order 无该项），
> 脚本 `run_bbnorm.py` 已迁至 `data_preprocessing_pipeline/run_bbnorm.py`。
> 作为数据预处理管线的**可选阶段**，两个入口均默认不跑、输出统一为 `00c_BBnorm/`：
> - `data_preprocessing_pipeline/data_preprocessing.py --stage bbnorm`（`--stage all` 不含）
> - `public_metadata_pipeline/preprocess_unified.py --bbnorm`（从 .sra 起步的公共数据入口）
>
> 用途: 多样本深度差异大时、共组装前的覆盖度归一化（target=70, mindepth=2）。
> 产出目录可直接作为本管线 `--input_reads`；若与预处理共用输出根目录，
> `report_pipeline.py` 会在报告里顺带展示该阶段（仅当目录存在）。

---

## 5. assembly_pipeline.py — 序列组装

**用途:** 统一宏转录组组装。支持 MEGAHIT, rnaviralSPAdes, Penguin。自动 PE/SE 检测，refineC 拆分合并后处理，资源使用追踪。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-t` / `--tool` | 必需 | megahit / rnaviralspades / penguin / all |
| `-i` / `--input` | 必需 | 输入文件或目录 |
| `-l` / `--length` | 0 | 最小 contig 长度 |
| `-n` / `--threads` | 8 | 线程数 |
| `-m` / `--memory` | 64 | 内存 (GB) |
| `-j` / `--jobs` | 1 | 并行任务数 |
| `--refineC_split` / `--refineC_merge` | — | refineC 后处理 |
| `--refineC_min_id` | 0.95 | refineC 最小相似度 |
| `--refineC_min_cov` | 0.50 | refineC 最小覆盖度 |
| `-o` / `--output-dir` | `./results` | 输出目录 |

### 输入
FASTA/FASTQ 文件或目录

### 输出
```
{output}/{sample}/
  {sample}_{tool}.contig.fasta    ← 组装 contig (加样本前缀)
  {sample}_{tool}.scaffolds.fasta ← rnaviralspades scaffolds
  *.{log, time.mem.log}           ← 日志 + 资源使用
  {sample}_{tool}_refineC/        ← refineC 拆分输出
```

**工具依赖:** megahit, rnaviralspades.py, penguin, refineC

---

## 6. virus_identification.py — 10 工具并行病毒鉴定

**用途:** 全自动病毒鉴定管线 v16.0。集成 10 个工具: Genomad, Diamond BLASTX (VP+NR+UniProt), RdrpCatch, ViraLM, VirBot, VirSorter2, ViralVerify, VirHunter, Metabuli, Viroid BLASTN。产出统一候选 FASTA + Venn/UpSet 可视化。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input` | 必需 | FASTA 文件或目录 |
| `-o` / `--output` | `5.virus_identification` | 输出目录 |
| `--identify_tools` | `all` | 工具选择 |
| `--db_dir` | `~/database/virus-db` | 数据库根目录 |
| `--virus_protein_db` / `--nr_db` / `--uniprot_db` | — | Diamond 数据库 |
| `--blast_evalue` | `1e-5` | E-value 阈值 |
| `--blast_mode` | `filter` | strict / filter / both |
| `--viralm_path` | `utils/viralm_cpu.py` | ViraLM 脚本路径 |
| `--virbot_path` | `../biosoft/VirBot/VirBot.py` | VirBot 路径 |
| `--virhunter_path` | `../biosoft/virhunter/predict_cpu.py` | VirHunter 路径 |
| `--virsorter_group` | `dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae` | VirSorter2 组 |
| `--jobs` | 1 | 并行样本数 |
| `--threads` | 20 | 线程数 |

### 输入
组装 contig FASTA 目录

### 输出
```
{output}/{sample}/
  {sample}.virus.candidate.fasta          ← ≥500bp 病毒候选
  {sample}.viroids.candidate.fasta        ← 200-1000bp 类病毒候选
  {sample}_virus.{tool}.result.id         ← 各工具结果 ID
  {tool}_output/                          ← 各工具详细输出
  {sample}_virus.all.candidate.fasta      ← 合并候选
  comparison_plots/                       ← Venn/UpSet/柱状图
  {sample}.processing.done                ← 完成标记
```

**工具依赖:** seqkit, genomad, diamond, blastn, virsorter, viralverify, metabuli, taxonkit, conda (rdrpcatch/viralm envs)

---

## 6b. graph_propagation.py — 组装图标签传播（graph 阶段；**默认不运行**）

**默认状态:** `graph.enabled: false`（2026-10-07 裁定）——`--stage all` 自动跳过本阶段，COBRA 吃 filter 直通镜像；需要时 `--stage graph` 显式运行或配置改 true。

**用途:** 在组装图拓扑层面重估每条 contig 的 viral/non-viral 归属（graph 阶段，插在 identification 与 filter 之间）。思想源自 GMW (Chen et al. 2026, Sci China Life Sci 69:1910-1917)，本实现为其配对读段改编（不使用 GMW 代码）：组装图 GFA → 配对读段链接图（assembler 无关）；Kraken2 单源种子 → 本管线 10 工具 ensemble 阳性清单；Louvain 社区传播 + GMW ≥10 倍多数比规则；断边沿用 GMW 口径（GC 差 >20% / 深度比 >20:1）；
链接=双端都比对上的主比对（跨 contig pair 在 bwa 语义下永远 non-proper，默认不要求 proper flag）。

### 判定语义

| verdict | 含义 |
|---------|------|
| `viral_seed` | 被 ≥`seed_min_votes`(默认2) 个工具判为病毒（种子，不可变） |
| `nonviral_seed` | len ≥`nonviral_min_len`（被工具检查过）且 0 票且非 viroid 阳性（种子，不可变） |
| `graph_viral` / `graph_nonviral` | unknown 落在 ≥10 倍多数比的社区 |
| `conflicting` | 社区两种种子都有但均未达多数比（保持 unknown） |
| `unknown` | 其余（无链接/无种子/短 contig 未被工具检查过） |

### 参数（编排器透传；默认值见 pipeline_config.yaml `graph:` 块）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--asm-dir` / `--ident-dir` / `--reads-dir` | io_layout | 输入三件套（d_asm / d_ident / d_hostdep） |
| `-o` / `--output-dir` | `d_graph` | 输出目录 |
| `-t` / `-j` | 8 / 1 | 单样本映射线程 / 并行样本数 |
| `--min-link-pairs` | 3 | 建边最少跨 contig 配对读段数 |
| `--gc-diff` / `--depth-ratio` | 0.20 / 20.0 | 断边阈值（GMW 口径） |
| `--require-proper` | off（勿开） | 仅收 proper pair — bwa 的 0x2 只在同 contig 时置位，跨 contig pair 永远 non-proper，开启即零链接（246 实测 0/196079） |
| `--resolution` | 1.0 | Louvain 分辨率（GMW 论文用 2） |
| `--seed-min-votes` / `--nonviral-min-len` | 2 / 1000 | 种子阈值（票数<阈值的弱阳性归 unknown 交图裁决） |
| `--majority-ratio` | 10.0 | 社区多数比 |
| `--max-reads` | 0 | 每端取前 N 条（确定性前缀截取，规模兜底） |
| `--map-min-len` | 1000 | 映射参考只保留 ≥该长度的 contigs（0=全量；短 contig 无可执行判定，剔除省 30-60% 映射时间，unknown 裁决不受损） |
| `--keep-bam` | off | 保留映射 BAM 供调试 |
| `--no-candidates` | off | 不写派生候选目录 |

### 输出

```
d_graph/{S}/
  {S}.graph_verdict.tsv     ← 逐 contig 判定+审计列 (contig/source_tool/length/gc/depth/
  │                            votes/seed_label/community_*/verdict/reason)
  {S}.graph_summary.tsv     ← 社区级汇总
  {S}.graph.done            ← 样本级断点标记
d_graph/candidates/{N}/
  {N}_virus.all.candidate.fasta  ← 派生候选 = ident 候选 − graph_nonviral
                                   (子目录命名逐镜像 d_ident, 供 COBRA raw 模式发现)
```

**消费点:** ① filter_virus.py 末端 pass（rescue/demote，见 6c 节）；② filter 镜像：COBRA raw 候选源优先吃 `d_filter/candidates`（filter 最终候选集的 02a 兼容镜像，2026-10-07 恢复 filter→COBRA 主链契约）；graph 派生候选降为 filter 未跑时的兜底。cobra_pipeline.py 内部零改动）。

**假阳性削减量化（蒙特卡洛，`tests/montecarlo_fp_reduction.py`，参数按植物宏转录组校准）：**
候选池假阳性按票数三层分解——≥2 票成为 viral_seed（不可变，结构性不由本模块削减，
由下游 CDD/rescue 兜底）；=1 票弱阳性为可削区；0-evidence 保留者归 filter 层图裁决。
模拟结果（600 contigs × 15 重复 × 6 参数组合，tool FP 率 0.5%–2%）：**全 FP 池削减 68–79%，
可削区（1 票）收割率 76–82%，真病毒候选误伤 0**（含病毒-宿主嵌合边扰动）。
真实数字服务器实跑后直接可数：verdict TSV 的 `graph_nonviral` 行数（COBRA 池削减）+
filter `graph_calls.tsv` 的 `demoted` 行数（实削数）。

**246 真实数据实跑（2026-10-02，N20220441-06-01 项目 7 样本，173,882 contigs；计票 bug 修正后数字）：**
① 候选池真实票数结构 = **≥2 票 565 (8%) + 1 票弱阳性 6,440 (92%)**；修复计票 bug 后
**候选层削减 105 条 (1.50%)、救回 78 条，TP 误伤 0**。初版报告的"候选池全强票、FP 削减 0%"
是 all.result.id 计票 bug 的伪象（并集文件被当独立工具，每条候选虚增 1 票）。
② RVDB 核酸 blastn 对本项目候选层**不具判别力**：viral_seed 抽样 70/70 零命中
（工具经蛋白/HMM 级鉴定的病毒，核酸层 e<=1e-5 比对不上属预期）；被降级 0/105 与被救回
0/78 对称 → 无危害信号但也无法用该代理度量收益。真值裁决以 mock benchmark 为准。
③ **mock benchmark（真值三臂，tests/score_mock_benchmark.py）**：wgsim 模拟 4 植物 RNA
病毒（TMV 100×/LYMov 50×/SPFMV 20×/AMV 10×）spike-in 进 NX-4 真实去宿主背景（447,867 对）
→ megahit → 6 工具 ensemble → graph。候选池 36 条（精确率 13.9%，1 票弱阳性 25 条）：
graph 臂削减 3 条背景 FP（精确率 15.2%）、**0 TP 误伤、四病毒召回完全不变**
（TMV 100.3% / SPFMV 99.8% / AMV 98.4% / LYMov 67.2%）；0 票 graph_viral 若注入会拉低
精确率（12.8%）→ "不注入"边界设计获真值支持。真实削减幅度（候选层 1.5%/9.7%FP）远低于
蒙特卡洛的 68-79%，原因是真实组装中孤立 contig 占比高（无链接上下文的 FP 图不可达，
保守放行）——蒙特卡洛高估了社区连通比例，此为重要校准。
④ **map_min_len 加速实测（246 mock A/B）**：参考裁剪 8,710→4,339 条后映射 92s→75s
（真实组装短 contig 占比更高时收益更大）；候选层 36 条裁决中 34 条完全一致、
4 条 unknown→graph_nonviral 且经真值确认全部为背景 FP（裁剪更优）；真病毒召回不变。
代价：丢失短 contig 桥接证据（连接性差的样本可设 0 回退全量）。

**已知边界:** SE 样本/无 reads → 全 unknown 优雅退化；graph 重跑不自动失效 filter/cobra 的 `.ok`（编排器会打显式 warning，需 `--force` 重消费）；短 contig 的 graph_viral 仅记录在 TSV，不注入候选池。

**工具依赖:** bwa-mem2, samtools, coverm, networkx（均在 pixi 环境）

**独立运行与验证:** 本模块可脱离编排器单独运行（所有目录均为 CLI 参数，io_layout 仅作兜底）；
`tests/validate_graph_propagation_standalone.py` 提供无外部映射工具的两场景验证
（A: 合成混合感染社区 16/16 contig 判定全对齐、五基因组社区零混合；B: GMW 论文示例真实 SPAdes
组装图 377 段 + 本机 blastn 亚型真值：断边阈值对真实图温和（边存活 100%）、拓扑救回 4 段全部
经 blastn 确认为流感、Louvain 社区 54/57 单亚型），用法
`python virome_discovery_pipeline/tests/validate_graph_propagation_standalone.py --gmw-repo <gmw仓库>`。
**与 GMW 真工具头对头（2026-10-02, Windows 本机）:** 本机安装 python-louvain + 重建 blast 库后
GMW 全流程可跑（`--disable_taxon_unfold` 绕开 Kraken2）。同图对比：GMW 参考展开剪掉 67% 连接
(482→158, 用参考信息)、本模块三规则剪边 0%（保连通性交给 Louvain+种子）——两种互补哲学。
注意：GMW 示例图 377 段经 blastn 全部为流感（H3N2 362/H5N6 15），**无污染可去**，它是
混合感染分型基准而非去污染基准；单目标模式（H5N6=目标, 6 个长段做种子）下本模块传播精度
77.8%/召回 63.6%，2 个假阳性恰为论文所述保守段重组嵌合体（GMW 靠参考展开解决，新病毒无参考
场景下该差异是两工具定位分歧点）。出版级去污染对比需按 GMW Methods 用 ART 重建其
Paramecium/Klebsiella 模拟数据集（仓库未附带，服务器任务）。

---

## 6c. filter_virus.py — UniProt + CDD 两层过滤（含图判定末端 pass）

**用途:** 对 identification 候选做两层过滤：CDD tiering（tier1=tier 病毒域 → viral；tier2/3 → ambiguous；仅宿主/非毒 → non_viral；零命中 → unclassified）→ UniProt/DIAMOND blastx（病毒 taxid / 关键词 / 弱比对抢救 / 全局 no-hit 保留）。**末端 graph pass**（仅 `-m filter`）消费 graph 阶段判定 TSV：demote = 零正面证据保留（cdd unclassified 且无 blastx hit）且判定 `graph_nonviral` → 移出；rescue = 被整链丢弃且判定 `graph_viral` → 加回。分层符号链接（uniprot/cdd/02a 原始）不动，审计链保留。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input` | 必需 | 候选 FASTA（或 `--input-dir` 合并模式） |
| `-o` / `--outdir` | 必需 | 输出目录（d_filter） |
| `-m` / `--filter-mode` | `filter` | filter / strict / raw |
| `-t` / `--threads` | env THREADS | DIAMOND 线程 |
| `--uniprot-db` / `--virus-taxid` / `--cdd-*` | DEFAULTS | 数据库路径 |
| `--graph-verdict-dir` | io_layout d_graph | graph 判定目录（缺失则 pass 无操作） |
| `--no-graph-demote` / `--no-graph-rescue` | off | 关闭图判定降级/救回 |

### 输出（split 模式逐样本）

```
d_filter/{S}/
  {S}_virus.all.candidate.fasta        ← 02a 原始链 (不动)
  {S}.virus.candidate_cdd_filtered.fasta / ..._uniprot_filtered.fasta  ← 分层链 (不动)
  {S}.virus.candidate_filtered.fasta   ← 最终链 (graph 调整后重指向 graph/ 产物)
  graph/{base}_graph_adjusted.fasta    ← 图判定调整后保留集 (仅有调整时写)
  graph/{base}_graph_calls.tsv         ← 全量审计表 (contig/cdd_call/blastx_hit/verdict/action)
  candidates/{S}/{S}_virus.all.candidate.fasta ← COBRA 主链镜像 (02a 兼容结构;
                                   filter 最终集原样, 恢复主链契约)
  cdd/ uniprot/                        ← CDD 调用表 + blastx 明细 + unclassified.id 清单
```

**工具依赖:** mmseqs2 (CDD), diamond (UniProt)

---

## 7. cobra_pipeline.py — COBRA 序列延伸

**用途:** COBRA 病毒序列批量延伸。运行 BWA-MEM2 比对 → CoverM 覆盖度 → COBRA-Meta 延伸。每样本/组装工具独立运行。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--mode` | `mix` | virus / other / mix / all |
| `-s` / `--samples` | — | 样本列表 (逗号分隔或文件) |
| `-o` / `--output-dir` | `10.cobra_result/` | 输出目录 |
| `--jobs` | 1 | 并行样本数 |
| `--threads` | 20 | 线程数 |
| `--mink` | 21 | COBRA 最小 k-mer |
| `--maxk` | 141 | COBRA 最大 k-mer |
| `--resume` | True | 断点续传 |

### 输入
- Reads 目录 (去宿主后)
- Contigs 目录 (组装后)
- 病毒鉴定目录 (候选序列)

### 输出
```
{output}/{sample}/cobra_{tool}_result/
  {task_id}.sorted.bam + .bai    ← BWA-MEM2 比对结果
  {task_id}.coverage.txt         ← CoverM 覆盖度
  {sample}.{mode}.{tool}.cobra.fa ← COBRA 延伸结果
```

**工具依赖:** bwa-mem2, samtools, coverm, cobra-meta

---

## 8. cluster_pipeline.py — 病毒基因组聚类去冗余

**用途:** 三支路病毒基因组聚类 v3.0: (1) seqkit 长度过滤, (2) CD-HIT 参考引导预聚类 (合并 ICTV/NCBI 参考), (3) vclust Leiden 聚类新序列。输出 centroids + 分簇文件。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--input-fasta` / `-i` | 必需 | 输入 FASTA |
| `--output-dir` / `-o` | 必需 | 输出目录 |
| `-t` / `--threads` | 64 | 线程数 |
| `--min-length` | 500 | 最小长度 |
| `--ani` | 0.95 | ANI 阈值 |
| `--qcov` | 0.85 | 覆盖度阈值 |
| `--ref-genomes` | — | 参考基因组 FASTA 路径 (可多个) |
| `--skip-vclust` | False | 跳过 vclust |
| `--resume` | False | 断点续传 |

### 输入
COBRA 延伸后病毒 FASTA；可选的 ICTV/NCBI 参考基因组

### 输出
```
{output}/
  1_seqkit/virus.candidate.fasta         ← 长度过滤后
  2_cdhit/                                ← CD-HIT 聚类
    known_centroids.fasta                 ← 参考引导 centroids
    known_association.tsv                 ← contig-参考映射
    novel_contigs.fasta                   ← 新序列 (给 vclust)
  3_vclust/                               ← vclust 聚类
    vclust_clusters.tsv                   ← 聚类结果
    all.cluster.ref.fasta                 ← 所有代表序列
    split_fastas/                         ← 分簇 FASTA
  centroids/
    final_centroids.fasta                 ← 最终 centroids
    known_ids.txt
```

**工具依赖:** seqkit, vclust (prefilter, align, cluster, deduplicate)

---

## 9. virus_classifier.py — 8 工具病毒分类

**用途:** 病毒分类整合 v4.2。支持 8 工具: Genomad, Metabuli, CAT, Diamond LCA, VITAP, MMseqs2, ACVirus, vConTACT3。合并为 8 级 taxonomy，加权投票生成 combined_taxonomy.tsv。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-g` / `--genomes` | — | 输入 FASTA 文件 |
| `-i` / `--input-dir` | — | 批量输入目录 |
| `-t` / `--tools` | `all` | 工具选择 (逗号分隔) |
| `-o` / `--output-dir` | `./classify_output` | 输出目录 |
| `-p` / `--threads` | 20 | 线程数 |
| `-j` / `--jobs` | 1 | 并行样本数 |
| `--db-dir` | `~/database/virus-db` | 数据库根目录 |

### 输入
Centroids FASTA；各工具数据库

### 输出
```
{output}/{sample}.classed/
  {sample}_{tool}_taxonomy.tsv       ← 各工具原始结果
  {sample}_combined_taxonomy.tsv     ← 8 级合并 taxonomy
  {tool}_output/                     ← 各工具工作目录
```

**工具依赖:** genomad, metabuli, diamond, CAT_pack, VITAP, mmseqs, ACVirus, vcontact3, taxonkit

---

## 10. virus_classifier_analysis.R — R 共识整合

**用途:** 读取 virus_classifier.py 输出的 combined_taxonomy.tsv，运行加权投票共识引擎，生成 final_integrated_classification.tsv (含 completeness, confidence, per-rank agreement)。

### 参数

| 参数 | 说明 |
|------|------|
| `--combined` | virus_classifier.py 的 combined_taxonomy.tsv |
| `--mmseqs` / `--vcontact3` / `--vitap` 等 | 各工具独立输出 (可选) |
| `-o` / `--output` | 输出目录 (默认 `.`) |
| `--cores` | 线程数 |

### 输入
combined_taxonomy.tsv；可选各工具独立输出

### 输出
```
{output}/
  final_integrated_classification.tsv   ← 加权投票共识 (8 级 + agree + confidence)
  agreement_stats.tsv                   ← 各工具各层级一致性
  consistency_summary.tsv               ← 一致性汇总
  intersection_venn.pdf / upset.pdf     ← 工具交集可视化
  agreement_rates.pdf                   ← 一致率柱状图
```

**R 包依赖:** optparse, data.table, ggplot2, VennDiagram, ggVennDiagram, UpSetR, cowplot, patchwork

---

## 11. run_host_prediction.py — 宿主预测

**用途:** 集成宿主预测: RNAVirHost → PhaBOX2 (CHERRY) → ICTV C9 分类查找。决策树整合 (ICTV > RNAVirHost > PhaBOX2)。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input` | 必需 | 输入 FASTA |
| `--tax` | 必需 | taxonomy TSV (来自 virus_classifier) |
| `-o` / `--output-dir` | `host_out` | 输出目录 |
| `-t` / `--threads` | 40 | 线程数 |
| `--phabox-db` | `~/database/virus-db/phabox_db_v2_2` | PhaBOX2 数据库 |
| `--mode` | `all` | all / ICTV / RNAVirHost / PhaBOX2 |

### 输入
Centroids FASTA + taxonomy TSV

### 输出
```
{output}/
  RVH_result/           ← RNAVirHost 输出
  phabox2_output/       ← PhaBOX2 输出
  C9_ICTV_result/       ← ICTV C9 分类
  ensemble_host_summary.tsv  ← 集成宿主预测
  host_classified_fasta/     ← 按宿主分类的 FASTA
```

---

## 12. rescue_pipeline.py — 四支路病毒抢救

**用途:** A 支路: CheckV 质控 (completeness > 90% 直通)。B 支路: Virseqimprover reads 延伸。C 支路: BLASTN 参考比对 + ragtag 参考引导延伸。D 支路: genus_len 属水平长度拯救 (±15%)。四路 HQ 序列 vclust 去重。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--centroids` / `-c` | 必需 | Centroids FASTA |
| `--clusters-tsv` | 必需 | vclust clusters.tsv |
| `--split-dir` | 必需 | split_fastas/ 目录 |
| `--output-dir` / `-o` | 必需 | 输出目录 |
| `--fastq-dir` / `-fq` | 必需 | 原始 reads 目录 |
| `--checkv-db` / `-cv` | 必需 | CheckV 数据库 |
| `--blast-db` | — | BLAST 数据库 |
| `--checkv-threshold` | 90.0 | CheckV 完整性阈值 |
| `--threads` / `-t` | 64 | 线程数 |
| `--jobs` / `-j` | 4 | 并行任务 |

### 输入
Centroids FASTA + clusters.tsv + split_fastas/ + reads 目录

### 输出
```
{output}/
  branch_a/ ← CheckV 直通
  branch_b/ ← Virseqimprover 延伸
  branch_c/ ← BLASTN 延伸
  merged/all_HQ.fasta                ← 三路合并
  centroids/final_centroids.fasta    ← 去重后 HQ vOTU catalog
```

**工具依赖:** checkv, blastn, Virseqimprover.py, salmon, ragtag, makeblastdb, blastdbcmd, vclust

---

## 13. report_pipeline.py — 管线报告生成

**用途:** 收集所有阶段数据，生成 TSV 汇总 + Sankey 分类图 + 交互式 HTML 报告。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-o` / `--output-dir` | 必需 | 管线输出根目录 |
| `--skip-sankey` | False | 跳过 Sankey 图 |
| `--skip-html` | False | 跳过 HTML 报告 |
| `--ai-summary` | False | AI 总结 |
| `--ai-provider` | `openai` | AI 提供商 |

### 输入
管线输出根目录 (含 00a-08 各阶段子目录)

### 输出
```
09_Reports/
  pipeline_report.html               ← 交互式 HTML
  stage_summary.tsv                  ← 各阶段状态/指标
  data_summary.tsv                   ← QC 统计
  assembly_summary.tsv               ← 组装统计
  ident_summary.tsv                  ← 鉴定统计
  checkv_summary.tsv                 ← CheckV 质量
  sankey_*.png / sankey_*.html       ← Sankey 图
```

---

## 14. 辅助脚本

| 脚本 | 用途 | 类型 |
|------|------|------|
| `data_preprocessing_pipeline/preprocess.py` | 清洗+去宿主轻量包装 | 独立运行 |
| `data_preprocessing_pipeline/data_preprocessing.py` | 清洗+去宿主可配置版 | 独立运行 |
| `viroid_circular_detect.py` | 类病毒环状检测 (self-BLASTN) | 独立运行 |
| `Virseqimprover.py` | reads 级别迭代延伸病毒基因组 | 被 rescue 调用 |
| `../virome_submission_pipeline/virome_submission_pipeline.py` | GenBank 提交准备 | 独立运行 |
| `utils/classify_contigs.py` | ICTV 参考分类查找 | 被 run_host_prediction 调用 |
| `utils/discovery2analysis.py` | centroids → analysis 格式转换 | 桥接脚本 |
| `utils/validate_novel_viruses.py` | 基于分类层级判断新病毒 | 独立运行 |
| `terminal_repeats.py` | 末端重复检测 (DTR/ITR) + 环形判断 + 起点归一 | 独立运行 |

### 09b 证据链脚本（analysis_verify 阶段）

「一道门 + 两次补票 + 整合」结构：

| 环节 | 脚本 | 说明 |
|---|---|---|
| 门 | `validate_rescue_cdd.py` | CDD curated 结构域验证（覆盖窄、特异性高） |
| 补票 1 | `gen_final_judgement.py` | blast 身份注释（known / new） |
| 补票 2 | `hmm_ct3_evidence.py` | CT3 五库全长 profile（覆盖宽），模型名六桶分级 + `VIRAL_CLUSTER_TOKENS` 支系识别 |
| 整合 | `integrate_rescue_evidence.py` | 五层证据加权 → `rescue_evidence_scored.tsv`；打分 `0.3*blastn+0.3*blastx+0.4*domain`，域证据层 `domain = max(cdd, hmm)`（HMM 进分默认开，`--no-hmm-in-score` 回退纯 CDD；同为 profile 证据取 max 防双重计分，CDD 盲区由 HMM 补上） |

> ⚠️ **CT3 不是纯病毒库，不能「命中即病毒」**：`phrogs_for_ct` 是噬菌体专属；`Useful_Annotation` 含大量
> `cdd_cluster`（COG/PHA 噬菌体结构蛋白）与 `baits_missed`（噬菌体尾部/门户）；`DNA_rep` 含
> `circular_genetic_element` / `plasmid` 这类非病毒环状元件。因此按**模型名**分六桶：
> VIRAL_FAMILY / VIRAL_CLUSTER / VIRAL_FUNCTION / AMBIGUOUS / PHAGE / OTHER，前三类合起来才算非噬菌体病毒证据。
>
> 早期曾用 RVDB-prot-HMM（`hmm_rvdb_evidence.py` + `_tools/build_famtax.py`）作补票 2，已由 CT3 路径取代
> （2026-09-14）：CT3 五库覆盖面更广，且两者一为二进制 `.h3m`、一为已 press 的文本，物理合并需转文本重压
> 且混两套命名体系，得不偿失。RVDB 两个脚本归档保留于 `archive/server_local_diff_20260922/server_only_scripts/`，**不接入流程**。


### terminal_repeats.py 说明

逐条 contig 检测末端重复: 5' 端拷贝出现在后半段且一直延伸到 3' 端 → 末端重复
(DTR, 环状装配的标记), 修剪掉冗余拷贝; 两端互为反向互补 → 反向重复 (ITR)。
判定逻辑与 Cenote-Taker3 的 `terminal_repeats.py` 逐行一致 (已实测 120 条合成
contig 600 格零差异), 另加容错 DTR 模式与基于完整 ORF 的起点归一。

```bash
python terminal_repeats.py -i contigs.fasta -o out_dir -t 8 \
    --method auto --rotation repeat
```

关键参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--method` | auto | `exact` 只有完全一致; `auto` 先 exact 再容错 (min-identity 0.90); `fuzzy` 只容错 |
| `--rotation` | repeat | `repeat` 修剪后起点即末端重复处 (位移 0); `orf` 起点归一到最长完整 ORF 起始密码子 (CT3 实际做法); `none` 只修剪 |
| `--min-dtr` / `--max-itr` | 20 / 1000 | DTR 最短长度 / ITR 最长长度 |
| `--max-length` | 1e6 | 超过此长度的 contig 不做 DTR 判定 (置 NA), 与 CT3 一致 |
| `--circ-file` | 无 | 用户直接声明环状的 contig 名单 (id / description / 完整 header 均可) |
| `--header` | full | 输出 FASTA 写完整 header; `id` 只写 id |

输出: `trimmed_contigs.fasta` (修剪/旋转后序列)、`terminal_repeats_summary.tsv`
(前 5 列与 CT3 同名同序, 追加 dtr_length/dtr_identity/dtr_method/itr_length/
topology/rotated/rotation_offset/rotation_mode)、`circular_contigs.txt`。

注意: ITR 单独不作为环状证据 (`topology` 只取 `circular_DTR` / `circular_user` /
`undetermined`)。CT3 的 `--circ-file` 只比对 description 去掉 id (见其源码
`fmt_desc`), 与其帮助文本「names of contigs (header line sans '>')」相悖, 本实现
同时接受 id / description / 完整 header。

自测: `python tests/test_terminal_repeats.py` (44 项断言)。
