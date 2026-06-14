# MMPV-RNA v2.3 全脚本解读

> **Meta-transcriptomic Mining of Plant Virome from RNA-seq data**
> 宏转录组植物病毒端到端发现管道 — 从原始 RNA-seq reads 到 HQ vOTU 目录

---

## 一、项目总览

MMPV-RNA 是一个完整的宏转录组病毒发现与分析框架，包含 **20 个 Python/R 脚本**和 **13 个文档**，覆盖从原始 FASTQ 到最终病毒目录的全流程。

### 核心数据流

```
Raw FASTQ → clean → deplete → assembly → identification → COBRA → cluster
                                                                      ↓
         ┌──────────────────────────────────────────────────────────────┘
         ↓            ↓           ↓
    taxonomy       host        checkv
         ↓            ↓           ↓
         └────────────┼───────────┘
                      ↓
                   rescue (A→C→D 三支路)
                      ↓
                  validate
```

### 并行支路 (已知病毒分析)

```
HostDepletion → auto_known_virus
                  ├─ detect (batch_virus_depth40)
                  ├─ variants (batch_virus_variants)
                  └─ full (batch_virus_full)
```

---

## 二、各脚本详细解读

### 1. `virome_pipeline.py` — 主控编排器 (1786 行)

**角色**: 整个框架的中央调度器，管理 10 个阶段的串行/独立执行。

**核心类**: `ViromePipeline`

**架构特点**:
- **双通道日志**: 控制台 INFO + 文件 DEBUG
- **断点续传**: 每个子脚本支持 `--resume/--force`
- **模糊样本匹配**: `fuzzy_match()` 处理 `_clean` 等命名字缀差异
- **自动数据库推导**: `--host_db` 自动查找 kraken2/bowtie2/hisat2 子目录
- **多工具组装合并**: ≥2 工具自动启用 `refineC split/merge`
- **宿主过滤**: rescue 阶段前按 Final_Host 分离 centoids
- **CD-HIT known 免拯救**: known + CheckV(≥90%) 直接输出

**10 个阶段映射**:

| 阶段 | 方法 | 子脚本 |
|------|------|--------|
| clean | `run_clean()` | clean-data.py |
| deplete | `run_depletion()` | host_depletion.py |
| assembly | `run_assembly()` | assembly_pipeline.py |
| identification | `run_identification()` | virus_identification16.py |
| cobra | `run_cobra()` | cobra_pipeline.py |
| cluster | `run_cluster()` | cluster_pipeline.py |
| taxonomy | `run_taxonomy()` | virus_classifier2.py + R |
| host | `run_host()` | run_host_prediction.py |
| checkv | `run_checkv_stage()` | checkv completeness |
| rescue | `run_rescue()` | rescue_pipeline.py |

**关键设计决策**:
- `--stage` 支持独立/串行执行，非 reads 依赖阶段自动 skip_clean/skip_depletion
- COBRA 阶段 jobs 减半 (`max(1, self.args.jobs // 2)`) 因为内存密集
- assembly 阶段自动检测已有工具输出，避免重复传递

---

### 2. `clean-data.py` — 数据清洗 (458 行)

**角色**: 原始 FASTQ → 清洁 FASTA，三步流水线。

**流程**: Fastp 质控 → Seqkit FASTQ→FASTA 转换 → Clumpify 光学去重

**关键类**:
- `UI`: 终端彩色进度条，线程安全的 `print_lock`
- `CheckpointManager`: 断点文件 `.clean_checkpoints` 持久化
- `FastqScanner`: 智能配对 R1/R2，支持 `_R1.` / `_R1_` / `_1.` 多种命名
- `CleaningPipeline`: 每个样本独立 run_single

**技术细节**:
- Fastp: `--qualified_quality_phred 20 --length_required 50 -g --poly_g_min_len 10`
- Seqkit: `fq2fa -w 0` 转换 FASTA (不换行)
- Clumpify: BBMap 去重，`reorder dedupe subs=0`
- 内存探针: Linux 下 `/proc/{pid}/statm` 实时监控
- 中间 FASTQ 自动清理 (节省磁盘)
- `--remove-raw` 危险选项：成功后删除原始数据

---

### 3. `host_depletion.py` — 去宿主 + 去 rRNA (918 行)

**角色**: Kraken2 分类 → 精准比对去宿主 → rRNA 剔除，三阶段混合管道。

**核心类**: `HybridPipeline`, `Aligners`, `SAM`, `CheckpointManager`

**工具支持矩阵**:

| 阶段 | 工具 | 说明 |
|------|------|------|
| Kraken2 | `kraken2 --confidence` | 物种分类标记 |
| 精准比对 | Bowtie2 / HISAT2 / Minimap2 | 可配置 `--tool` |
| rRNA 去除 | Ribodetector / SILVA Bowtie2 | 可配置 `--rrna_tool` |

**关键设计**:
- **消融实验支持**: `--steps kraken2,align,rrna` 可任意组合
- **SILVA 模式**: Bowtie2 `--very-sensitive-local` 比对 SILVA 库，`--un-conc-gz` 直接输出
- **资源监控**: 每样本记录耗时/内存/CPU 到 TSV
- **Seqkit 阶段统计**: 每步跟踪 reads 数量变化 (Raw → Kraken2 → Host Filtered → rRNA Filtered)
- **自动配对**: `auto_pair_files()` 支持 PE/SE 自动识别
- **可视化**: matplotlib 绘制 reads 数量变化折线图/柱状图

**SAM 操作链**:
```
align → samtools view -f 12 -F 256 → samtools sort -n → samtools fastq -n
```

---

### 4. `assembly_pipeline.py` — 宏转录组组装 (约 2000+ 行)

**角色**: 三工具并行组装 + refineC 后处理。

**支持工具**: Penguin (宏转录组专用) / MEGAHIT / rnaviralSPAdes

**关键类**: `AssemblyPipeline`

**核心功能**:
- **自动数据类型检测**: PE/SE, FASTQ/FASTA 自适应
- **refineC 集成**: split (拆分嵌合体) + merge (合并重叠 contig)
- **多工具并行**: `--assembler all` 同时运行三工具
- **资源监控**: 每工具记录 wall_sec/cpu_sec/mem_mb
- **优雅退出**: SIGINT/SIGTERM 信号处理

**refineC merge 参数**:
- `--refineC_min_id 0.97` (最小序列一致性)
- `--refineC_min_cov 0.50` (最小覆盖度)
- `--refineC_frag_min_len 1000` (最小片段长度)

---

### 5. `virus_identification16.py` — 病毒序列鉴定 (约 2000+ 行)

**角色**: 6 工具并行病毒鉴定 + Venn/Upset 图可视化。

**6 工具矩阵**:

| 工具 | 方法 | 适用 |
|------|------|------|
| Genomad | 深度学习 (病毒/质粒/前病毒) | 全基因组 |
| Diamond BLASTX | 蛋白比对 (RefSeq + NR) | 蛋白水平 |
| VirSorter2 | HMM 隐马尔可夫模型 | 全基因组 |
| ViralVerify | 病毒蛋白 HMM 验证 | 蛋白水平 |
| VirHunter | 机器学习 | 全基因组 |
| Metabuli | k-mer 分类 | 序列水平 |

**关键特性**:
- **多库对抗**: Diamond 多数据库 (RefSeq/UniProt/NR) 分层比对
- **高维 Venn 图**: `venn` 库支持 2-6 集合，`UpSet` 库支持 >6 集合
- **资源记录**: 每样本/工具独立 resource.tsv，最终合并汇总
- **Blast 救援机制**: 主库无结果时自动切换备用库

---

### 6. `cobra_pipeline.py` — COBRA 批量延伸 (1240 行)

**角色**: BWA-MEM2 → CoverM 覆盖度 → COBRA 重叠延伸 → CheckV 评估

**核心类**: `CobraPipeline`

**三种模式**: `virus` / `other` / `mix`，通过后缀 (`.unmapped.virus` 等) 自动识别

**病毒序列来源** (`--virus-mode`):
- `raw`: 原始鉴定结果
- `filter`: UniProt 过滤
- `strict`: 严格过滤 (默认)

**处理流水线**:
1. 查找 contig + virus + reads 三元组
2. 标准化 contig/virus 序列名称
3. BWA-MEM2 索引 + 比对 + Samtools sort
4. CoverM 计算覆盖度 (covered_fraction/mean/rpkm)
5. COBRA 重叠延伸 (`cobra-meta`)
6. 合并 COBRA_category_*.fasta → `.cobra.fa`

**COBRA 参数**: mink=21, maxk=141, linkage_mismatch=2

**断点管理**: JSON checkpoint 文件 (`checkpoint_status.json`)

---

### 7. `cluster_pipeline.py` — 聚类去冗余 (554 行)

**角色**: CD-HIT 参考引导预聚类 + vclust Leiden 聚类，产出 centroids。

**三步流程**:

1. **seqkit 长度过滤**: `≥min-length bp` (默认 500)
2. **CD-HIT 参考引导预聚类** (可选 `--ref-genomes`):
   - vclust deduplicate 去重参考
   - contig + 参考合并 (加 `ref|`/`our|` 前缀)
   - vclust cd-hit 聚类
   - 拆分 known (代表以 `ref|` 开头) / novel
3. **vclust Leiden 聚类**: 仅 novel 部分
   - prefilter → align → cluster 三件套

**输出结构**:
```
04_CLUSTER/
├── 1_seqkit/virus.candidate.fasta
├── 2_cdhit/
│   ├── cdhit_combined.fasta
│   ├── known_centroids.fasta
│   ├── known_association.tsv
│   ├── novel_contigs.fasta
│   └── known_clusters/cluster_*.all.fasta
├── 3_vclust/
│   ├── vclust_clusters.tsv
│   ├── split_fastas/cluster_*.all.fasta
│   └── cluster_summary.tsv
└── centroids/
    ├── final_centroids.fasta
    └── known_association.tsv
```

---

### 8. `virus_classifier2.py` — 病毒分类 (约 1000+ 行)

**角色**: 9 工具并行分类 + 8 级 taxonomy 整合。

**9 工具**:
- Genomad, Metabuli, CAT, Diamond LCA, MMseqs2, VITAP, ACVirus, vConTACT3, PhaGCN3

**lineage_to_ranks 算法**:
- 解析 `;` 分隔的 lineage 字符串
- 过滤 subrank (`viricotina` 等)、unclassified、skip 词
- 8 级标准化: Realm/Kingdom/Phylum/Class/Order/Family/Genus/Species

---

### 9. `virus_classifier_analysis14.R` — R 共识整合

**角色**: 多工具投票共识，优先级 `vcontact3 > vitap > acvirus > mmseqs > genomad`。

---

### 10. `run_host_prediction.py` — 宿主预测 (约 400+ 行)

**角色**: 三种工具宿主预测 + 决策树集成。

**三工具**:
1. **RNAVirHost**: 两步法 `classify_order` + `predict`
2. **PhaBOX2**: CHERRY 宿主预测
3. **ICTV (C9)**: 分类库查找宿主

**决策树规则**: `ICTV > RNAVirHost > PhaBOX2`
- ICTV==RVH → 直接采用
- 分歧时 PB2 决胜
- 全分歧 → ICTV > RVH > PB2

---

### 11. `C9_classify_contigs.py` — ICTV 宿主分类

**角色**: 基于概率查找表的级联宿主分类。

**分类级联**: Species → Genus → Family → Order (最特异优先)

**输出**: 按宿主类别分组的 FASTA + TSV + 置信度报告

---

### 12. `rescue_pipeline.py` — 三支路级联拯救 (594 行)

**角色**: 独立拯救脚本，接收已聚类的 centroids，执行三支路拯救。

**三支路**:
- **分支 A**: CheckV 并行评估 → completeness ≥ 90% → pass
- **分支 C**: Virseqimprover reads 延伸 (cluster 多样本聚合) → CheckV
- **分支 D**: BLASTN + CheckV + VSI 最后拯救

**免拯救**: CD-HIT known + CheckV pass(≥90%)

**关键优化**:
- 分支 C 聚合 cluster 内所有样本 reads 提升覆盖
- 分支 D BLASTN megablast 快速搜索
- 最终 vclust 去重合并

**VSI 调用**:
```bash
python Virseqimprover.py -1 reads_R1 -2 reads_R2 -scaffold ref.fa -o out -salmon salmon -t N
```

---

### 13. `Virseqimprover.py` — 迭代延伸引擎 (约 800+ 行)

**角色**: Salmon 定量 → BBMap 提取 → SPAdes 组装 → 循环迭代。

**迭代停止条件**:
- CheckV completeness > 90% (自动停止)
- 序列长度不再增长
- 达到最大迭代次数

**参数**: k-mer 自动/手动, 环状检测 (minOverlap=5000, minIdentity=95%)

---

### 14. `validate_novel_viruses.py` — 新颖性判断 (489 行)

**角色**: 基于分类层级判断病毒新颖性，无需 BLASTN。

**分类规则**:
```
Species ≠ NA              → ★ known (已知)
Genus ≠ NA, Species = NA  → ★★ novel_species (新种)
Family ≠ NA, Genus = NA   → ★★ novel_genus (新属)
Order/Class ≠ NA          → ★★★ novel_family/order (新科/目)
全是 NA                    → ★★★ truly_novel (全新)
```

**输出**:
- `novel_viruses.annotated.tsv`: 全 vOTU 注释表
- `final_virus_catalog.fasta`: 病毒目录 (ID 含分类标记)
- `validation_report.html`: Plotly.js 交互 HTML 报告

---

### 15. `auto_known_virus.py` — 已知病毒分析主控 (332 行)

**角色**: 已知病毒的检测 → 变异 → 全长组装三阶段。

**三阶段**:
1. **detect** (batch_virus_depth40): 快速检测 (Salmon/Kallisto/Bowtie2)
2. **variants** (batch_virus_variants): 变异分析 (FreeBayes/iVar/LoFreq + SnpEff + SnpGenie)
3. **full** (batch_virus_full): 全长组装 (SPAdes/IVA 等)

---

### 16. `batch_virus_depth40.py` — 已知病毒快速检测 (849 行)

**角色**: 伪比对/传统比对 + Poisson Ratio 去假阳性 + 双轨过滤。

**核心类**: `UnifiedVirusPipeline`

**双引擎**:
- **伪比对**: Salmon/Kallisto (极速)
- **传统比对**: Bowtie2/BWA/Minimap2/StrobeAlign

**双轨过滤系统**:
- **A 轨** (全基因组): coverage ≥ 阈值 AND Poisson_Ratio ≥ 阈值
- **B 轨** (基因区转录覆盖): gene_total_cov + gene_avr_cov
- RNA 病毒: A+B 双轨; DNA 病毒: 仅 A 轨 (转录组不适用全长)

**输出指标**: EM_Reads, CPM, FPKM, TPM, Avg_Read_ANI, Poisson_Ratio, Pi

**Poisson Ratio 原理**: `(Coverage%/100) / Predicted_Support` 判断覆盖均匀性

---

### 17. `batch_virus_variants.py` — 变异分析 (1230 行)

**角色**: 三段独立解耦的变异分析管道。

**类**: `PostProcessPipeline`

**三个独立模块**:
1. **变异检测** (`worker_call_variants`):
   - FreeBayes/LoFreq/iVar 三选一
   - 动态 VCF 过滤 (QUAL>20, DP/SAF/AF 阈值随深度自适应)
   - bcftools filter PASS/FAIL 标记

2. **SnpEff 注释** (`worker_run_snpeff`):
   - 自动 NCBI 下载 GenBank 构建本地数据库
   - 智能断点：检测 `snpEffectPredictor.bin` 存在则跳过

3. **SNPGenie 进化分析** (`worker_run_snpgenie`):
   - dN/dS 选择压力计算
   - 自动展平 SNPGenie_Results 子目录

**关键优化**:
- 虚拟 GTF 生成 (无 GenBank 注释时兜底)
- iVar TSV → VCF 反转录
- 等位基因频率提取 (AF 字段解析)

---

### 18. `batch_virus_full.py` — 全长组装 (253 行)

**角色**: 桥接 VAP 产出与 OmniVirusAssembler 的批量引擎。

**目录结构**: `{Taxonomy}_{Accession}/{Sample}_{Accession}/`

**Reads 来源优先级**:
1. 靶向提取 reads (extract_reads)
2. 原始 Clean Data (raw reads)

**断点续传**: 检测 `final.fasta/scaffolds.fasta/contigs.fasta` 存在则跳过

---

### 19. `preprocess.py` — 数据预处理 (159 行)

**角色**: clean-data + host_depletion 二合一快捷脚本。

**两阶段**: 清洗 (Fastp+Seqkit+Clumpify) → 去宿主 (Kraken2+Align+rRNA)

输出 reads 指针自动从前一阶段流转到下一阶段。

---

### 20. `viroid_circular_detect.py` — 类病毒环状检测

**角色**: 检测类病毒 (viroid) 的环状 RNA 特征。

---

## 三、技术架构总结

### 依赖层次

```
第一层 (基础设施):
  Python 3.10+, Polars, BioPython, pysam, pandas, tqdm, colorlog

第二层 (核心算法):
  fastp, seqkit, clumpify.sh, kraken2, bowtie2, hisat2, minimap2,
  samtools, bwa-mem2, coverm, cobra-meta, checkv, blastn, blastx,
  salmon, kallisto, megahit, spades, penguin, vclust, cd-hit,
  genomad, virsorter2, freebayes, lofreq, ivar, snpEff, snpgenie,
  rnavirhost, phabox2, R

第三层 (编排调度):
  subprocess/ProcessPoolExecutor/ThreadPoolExecutor 多进程并行
  Spawn 强制启动方法 (避免 Polars/C 扩展 Fork 死锁)
```

### 数据流设计模式

- **断点续传**: 所有脚本通过 checkpoint/--resume 支持
- **目录约定**: `00a_CleanData/`, `00b_HostDepletion/`, `01_Assembly/`...
- **样本追踪**: `extract_base_sample()` + `fuzzy_match()` 处理命名差异
- **资源监控**: 每样本/工具记录 wall_sec/cpu_sec/mem_mb
- **日志分离**: 主控 INFO → 控制台, DEBUG → 文件, 子进程独立日志
- **多进程安全**: spawn 强制启动, 文件锁 (Lock), 独立日志文件避免竞态

### 核心创新点

1. **CD-HIT 参考引导预聚类**: 将碎片化 contig 关联到 ICTV/NCBI 完整基因组
2. **三支路级联拯救**: CheckV → VSI(多样本聚合) → BLASTN+VSI 渐进提升 HQ 产出
3. **Poisson Ratio 建模**: 假阳性精准剔除 (覆盖均匀性检验)
4. **双轨过滤 (A/B)**: RNA 病毒全基因组+基因区双验证, DNA 病毒单轨
5. **宿主过滤优化**: 拯救前按宿主分类过滤, 节省 70%+ 计算量
6. **基于分类层级的新颖性判断**: 不依赖 BLASTN, 直接从 taxonomy 完整性判断
7. **SILVA rRNA 去除**: Bowtie2 比对替代 Ribodetector, 不依赖特定工具
```

---

## 四、命令行速查

```bash
# 全流程
python virome_pipeline.py --stage all --input_reads /data/raw/ --output_dir /data/out/ ...

# 已知病毒
python auto_known_virus.py --stage detect     --reads_dir ... --ref_info ... --reference ...
python auto_known_virus.py --stage variants   --reads_dir ... --ref_info ... --reference ... --snpeff --snpgenie
python auto_known_virus.py --stage full       --reads_dir ...

# 独立阶段
python virome_pipeline.py --stage clean       --input_reads ... --output_dir ...
python virome_pipeline.py --stage deplete     --output_dir ... --kraken2_db ... --host_align_db ...
python virome_pipeline.py --stage assembly    --output_dir ...
python virome_pipeline.py --stage identification --output_dir ... --virus_db ...
python virome_pipeline.py --stage cobra       --output_dir ...
python virome_pipeline.py --stage cluster     --output_dir ...
python virome_pipeline.py --stage taxonomy    --output_dir ... --virus_db ...
python virome_pipeline.py --stage host        --output_dir ... --phabox-db ...
python virome_pipeline.py --stage checkv      --output_dir ... --checkv_db ...
python virome_pipeline.py --stage rescue      --output_dir ... --checkv_db ...

# 交叉验证
python validate_novel_viruses.py -i centroids.fasta --taxonomy taxonomy.tsv --cdhit-known known.tsv ...
```

---

## 五、最终产出目录

```
out/
├── 00a_CleanData/           # 清洗后 reads
├── 00b_HostDepletion/       # 去宿主清洁 reads
├── 01_Assembly/{sample}/    # 三工具 contig
├── 02_Identification/{sample}/ # 候选病毒
├── 03_COBRA/{sample}/       # COBRA 延伸
├── 04_CLUSTER/              # 聚类 + centroids
│   ├── centroids/final_centroids.fasta
│   ├── centroids/known_association.tsv
│   └── 3_vclust/vclust_clusters.tsv
├── 05_Taxonomy/integrated/  # 分类注释
├── 06_HostPrediction/       # 宿主预测
├── 07_Checkv/               # CheckV 预评估
├── 08_Rescue/               # 三支路拯救 ★ HQ vOTU
├── known_viruses/           # 已知病毒分析
└── 09_Validation/           # 新颖性验证
```

---

*文档生成时间: 2026-06-14 | MMPV-RNA v2.3*
