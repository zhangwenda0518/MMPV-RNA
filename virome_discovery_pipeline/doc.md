# virome_discovery_pipeline — 病毒组发现流程文档

## 数据流总览

```
Raw FASTQ
  │
[00a_CleanData]  ← clean-data.py     (fastp → seqkit → clumpify)
  │
[00b_HostDepletion] ← host_depletion.py (kraken2 → align → rrna)
  │
[00c_BBnorm]  ← run_bbnorm.py      (可选, co-assembly 前归一化)
  │
[01_Assembly] ← assembly_pipeline.py (megahit / rnaviralspades / penguin)
  │
[02_Identification] ← virus_identification.py (10 工具并行)
  │
[03_COBRA]  ← cobra_pipeline.py     (bwa-mem2 → coverm → cobra-meta)
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
| `--stage` | `all` | 阶段选择: all,clean,deplete,bbnorm,assembly,identification,cobra,cluster,taxonomy,host,checkv,rescue,report |
| `--config` / `--profile` | `default` | YAML 配置文件和 profile |
| `--coassembly` | False | 合并所有样本 reads 进行共组装 |
| `--bbnorm` | False | 启用 BBNorm 归一化 |
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
  00c_BBnorm/            ← 归一化 reads (可选)
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

## 2. clean-data.py — 数据清洗

**用途:** RNA-seq 数据清洗: fastp (QC+去接头) → seqkit (FASTQ→FASTA) → clumpify (光学去重)。自动识别 PE/SE。

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
  3.clumpify/         ← 去重后 FASTA
  logs/               ← 日志 + fastp HTML/JSON
  .clean_checkpoints  ← 断点续传
```

**工具依赖:** fastp, seqkit, clumpify.sh

---

## 3. host_depletion.py — 宿主去除

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

## 4. run_bbnorm.py — 覆盖度归一化

**用途:** BBNorm k-mer 覆盖度归一化，降低高覆盖度序列偏差。自动识别 PE/SE。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input` | 必需 | 输入 reads 目录 |
| `-o` / `--output` | 必需 | 输出目录 |
| `-t` / `--threads` | 16 | 线程数/样本 |
| `-j` / `--jobs` | 1 | 并行样本数 |
| `--target` | 70 | 目标覆盖度 |
| `--mindepth` | 2 | 最低 k-mer 深度 |

### 输入
去宿主后 reads 目录

### 输出
`{output}/*_norm_R1.fq.gz`, `*_norm_R2.fq.gz`, `*_norm_SE.fq.gz`

**工具依赖:** bbnorm.sh

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

## 12. rescue_pipeline.py — 三支路病毒抢救

**用途:** A 支路: CheckV 质控 (completeness > threshold)。B 支路: Virseqimprover reads 延伸。C 支路: BLASTN 参考引导延伸。三路 HQ 序列 vclust 去重。

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
| `preprocess.py` | 清洗+去宿主轻量包装 | 独立运行 |
| `data_preprocessing.py` | 清洗+去宿主可配置版 | 独立运行 |
| `viroid_circular_detect.py` | 类病毒环状检测 (self-BLASTN) | 独立运行 |
| `Virseqimprover.py` | reads 级别迭代延伸病毒基因组 | 被 rescue 调用 |
| `../virome_submission_pipeline/virome_submission_pipeline.py` | GenBank 提交准备 | 独立运行 |
| `utils/classify_contigs.py` | ICTV 参考分类查找 | 被 run_host_prediction 调用 |
| `utils/discovery2analysis.py` | centroids → analysis 格式转换 | 桥接脚本 |
| `utils/validate_novel_viruses.py` | 基于分类层级判断新病毒 | 独立运行 |
