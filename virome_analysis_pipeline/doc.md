# virome_analysis_pipeline — 已知病毒分析流程文档

## 数据流总览

```
reads_dir/ + reference.fasta + ref_info.tsv
  │
[10] report ← generate_pipeline_report.py
  ▲
  │(reads all directories)
  │
[1] detect → batch_virus_depth.py       →  1_FastViromeExplorer/
[2] filter → utils/filter_summary.py    →  high_conf.summary.tsv
[3] variants → batch_virus_variants.py  →  2_Virus_variants_Results/
[4] full → batch_virus_full.py          →  3_Virus_assemblies_final/
[5] extract → utils/extract_full_fasta.py → 4_assemblies_clean/
[6] post → 6-script suite               →  5_post_analysis/
[7] capheine → capheine_pipeline.py     →  6_capheine/
[8] similarity → virus_auto_pipeline.py →  7_similarity/
[9] dvg → batch_virema_dvg.py           →  8_virema_dvg/
```

---

## 1. auto_known_virus.py — 分析编排器

**用途:** 10 阶段已知病毒分析编排器。顺序调用所有阶段脚本，管理检查点、日志、并发和数据流。对应 `auto_known_virus.py --stage all` 运行全流程。

### 10 阶段说明

| 阶段 | 脚本 | 输出目录 | 产出 |
|------|------|----------|------|
| detect | batch_virus_depth.py | `1_FastViromeExplorer/` | best.summary.tsv |
| filter | utils/filter_summary.py | `1_FastViromeExplorer/summary/` | high_conf.summary.tsv |
| variants | batch_virus_variants.py | `2_Virus_variants_Results/` | VCF + SnpEff + SNPGenie |
| full | batch_virus_full.py | `3_Virus_assemblies_final/` | 全长组装 |
| extract | utils/extract_full_fasta.py | `4_assemblies_clean/` | 最长 contig |
| post | 6 脚本并行 | `5_post_analysis/` | VCF viz + PCA + MAF + SnpGenie |
| capheine | capheine_pipeline.py | `6_capheine/` | HyPhy 正选择 |
| similarity | virus_auto_pipeline.py | `7_similarity/` | 相似性热图 |
| dvg | batch_virema_dvg.py | `8_virema_dvg/` | DVG/重组 |
| report | generate_pipeline_report.py | `Pipeline_Summary_Report.html` | 交互报告 |

### 关键参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--reads_dir` | — | 清洗后 reads 目录 |
| `--output_dir` / `-o` | — | 输出根目录 |
| `--ref_info` | — | 参考信息 TSV |
| `--reference` | — | 参考基因组 FASTA |
| `--stage` | `all` | 阶段选择 |
| `--threads` | 60 | 总线程数 |
| `--jobs` | 4 | 并行任务数 |
| `--tool` | `salmon` | 定量引擎 (salmon/kallisto/bowtie2/bwa/...) |
| `--variant_caller` | `ivar` | 变异 caller (freebayes/ivar/lofreq) |
| `--coverage` | 10.0 | 覆盖度阈值 |
| `--ratio` | 0.3 | Poisson ratio 阈值 |
| `--snpeff` / `--snpgenie` | False | 启用 SnpEff/SNPGenie |

### 输入
- 清洗后 reads 目录 (host-depleted FASTQ)
- 参考基因组 FASTA
- 参考信息 TSV (含 Accession, Taxid, Species, Segment 列)

### 输出
10 个子目录 + Pipeline_Summary_Report.html

---

## 2. batch_virus_depth.py — 快速定量 (Stage 1)

**用途:** 双引擎病毒定量管线。伪比对 (Salmon/Kallisto) 或传统比对 (Bowtie2/BWA/Minimap2/HISAT2)。Poisson Ratio 建模过滤假阳性，双轨过滤 (全基因组 + 基因区)，ANI/Pi 进化测算。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input_dir` | 必需 | FASTQ/FASTA 目录 |
| `-r` / `--reference` | 必需 | 参考 FASTA |
| `--ref_info` | 必需 | 参考信息 TSV |
| `-o` / `--output_dir` | `./virus_out` | 输出目录 |
| `--tool` | `bowtie2` | 比对/定量工具 |
| `-t` / `--threads` | 8 | 线程数 |
| `--coverage` | 10.0 | 覆盖度阈值 |
| `--ratio` | 0.3 | Poisson ratio 阈值 |
| `--meandepth` | 0.5 | 平均深度阈值 |
| `--min_tpm` | 1.0 | 最小 TPM |
| `--genes_cov` | — | 基因覆盖度文件 (可选双轨过滤) |
| `--min_aln_len` | 80 | 最小比对长度 |
| `--min_pid` | 0.90 | 最小相似度 |

### 输入
FASTQ/FASTA 目录；参考 FASTA；参考信息 TSV

### 输出
```
{output}/
  bam/           ← 比对文件
  stat/          ← pandepth 深度文件
  summary/       ← all_viruses.best.summary.tsv (最终检出表)
  index/         ← 参考索引
  plots/         ← 深度频率图 (Rscript)
  batches/       ← Parquet 检查点
  logs/          ← 资源使用报告
```

**工具依赖:** salmon/kallisto/bowtie2/bwa/bwa-mem2/hisat2/minimap2/strobealign, samtools, pandepth, coverm(可选), Rscript

### 关键输出解读

`all_viruses.best.summary.tsv`:
- **Rep_Accession:** 参考代表株 accession
- **Rep_Coverage(%):** 基因组覆盖度
- **Rep_MeanDepth:** 平均测序深度
- **Poisson_Ratio:** 泊松比 (<0.3 可能为局部堆叠假阳性)
- **Asm_TPM/CPM/RPM:** 标准化丰度
- **Uniq_Reads/Multi_Reads:** 唯一/多重比对 reads
- **Avg_Read_ANI:** 平均 reads 一致性

---

## 3. batch_virus_variants.py — 变异检测+注释 (Stage 3)

**用途:** 三阶段并行: (1) 变异检出 (FreeBayes/iVar/LoFreq) + 动态 VCF 过滤, (2) SnpEff 功能注释 (自动构建病毒数据库), (3) SNPGenie 进化分析 (dN/dS, piN/piS)。含病毒 reads 提取 + 共识序列构建。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--summary` | 必需 | 上游检出表 TSV |
| `--info` | 必需 | 参考信息 TSV |
| `--reference` | 必需 | 参考 FASTA |
| `--fastq` / `--bam` | 二选一 | reads 来源 |
| `--variant_caller` | `freebayes` | 变异检出工具 |
| `--snpeff` | False | 启用 SnpEff |
| `--snpeff_jar` / `--snpeff_config` | biosoft/snpEff/ | SnpEff 路径 |
| `--snpgenie` | False | 启用 SNPGenie |
| `-q` / `--vc_qual` | 20 | 最低质量 |
| `-d` / `--vc_depth` | 5 | 最低深度 |
| `-f` / `--vc_freq` | 0.5 | 最低频率 (%) |
| `--threads` | 8 | 线程数 |
| `--jobs` | 4 | 并行病毒数 |

### 输入
- 检出 summary TSV (Stage 1 产出)
- 参考信息/基因组
- FASTQ 或 BAM 目录

### 输出
```
{output}/
  virus-fasta/          ← 提取的参考序列
  virus-bam/            ← 病毒特异比对
  virus-consensus/      ← 共识序列
  virus-variants/       ← raw + filtered VCF
  virus-SnpEff/         ← 注释 VCF + TSV
  virus-SNPGenie/       ← population_summary.txt 等
  summary/              ← all_summary.tsv + 共感染矩阵
  logs/
```

**工具依赖:** samtools, bowtie2/bowtie2-build, pigz, viral_consensus, freebayes/ivar/lofreq, bcftools, java (snpEff), snpgenie.pl

### 结果解读

- **VCF:** `*.filtered.vcf` — 经过质量/深度/频率过滤的变异位点
- **SnpEff:** `*.ann.vcf` — 含 ANN 功能注释字段 (HIGH/MODERATE/LOW/MODIFIER)
- **SNPGenie:** `population_summary.txt` — dN/dS, πN/πS, Tajima's D

---

## 4. batch_virus_full.py — 全长组装调度 (Stage 4)

**用途:** 桥接变异分析到 OmniVirusAssembler。按病毒/样本配对，提取 reads 并启动 virus-full.py 12 步组装。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--downstream_dir` | 必需 | Stage 3 输出目录 |
| `--summary` | 必需 | 检出/变异 summary TSV |
| `--clean_data` | 必需 | 清洗后 reads 目录 |
| `--min_covered` | 0.0 | 最低覆盖度 |
| `--assembly_tools` | `all` | 组装工具 |
| `--jobs` | 4 | 并行数 |

### 输入
Stage 3 输出 + 清洗后 reads

### 输出
```
{outdir}/{Taxonomy}_{Accession}/{Sample}_{Accession}/
  (virus-full.py 12 步完整输出)
```

---

## 5. virus-full.py — OmniVirusAssembler (Stage 4 引擎)

**用途:** 12 步全长病毒基因组组装 v9.0。MEGAHIT/SPAdes → refineC → Divine Fusion → PVGA 延伸 → 迭代抛光 → 双引擎 Gap Filling → 环化检测 → 覆盖度可视化。

### 12 步骤

| 步骤 | 内容 |
|------|------|
| 1-2 | MEGAHIT/SPAdes/Penguin 组装 + refineC 合并 |
| 3 | Shiver-like 清理 |
| 4 | Divine Fusion 参考合并 |
| 5 | PVGA 延伸 + gap 分割 |
| 6-7 | 预融合合并 + rmDup 去重 |
| 8 | 最终 Divine Fusion |
| 9 | 迭代 reads 级别抛光 (viral_consensus) |
| 10 | 双引擎 Gap Filling (gmcloser + abyss-sealer) |
| 11 | 环化闭合检测 |
| 12 | 覆盖度可视化 |

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-r` / `--reference` | 必需 | 参考基因组 |
| `--assembly_reads` | 必需 | 组装用 reads |
| `-o` / `--output` | `Ultimate_Result` | 输出目录 |
| `-t` / `--threads` | 8 | 线程数 |
| `-n` / `--iter` | 3 | 抛光迭代次数 |
| `--assembly_tools` | `all` | 组装工具 |
| `--chk-len` | 0.98 | 完美判定长度比 |
| `--chk-n50` | 0.95 | 完美判定 N50 比 |

### 输出
```
{output}/{sample}/
  1.DeNovo_Assembly/ 到 12.Coverage_Visualization/
  11.Ultimate_Circular_Result.fasta  ← 最终环化结果
  Fully-assembled.ok                 ← 成功标记
```

**工具依赖:** megahit, rnaviralspades, penguin, refineC, mafft, minimap2, viral_consensus, blastn, samtools, bbmerge.sh, pvga, gmcloser, abyss-sealer, utils/genome_rmDuplicates.pl

---

## 6. 后处理脚本组 (Stage 6 — 6 脚本并行)

### virus_variants_analyzer.py — 9 图变异可视化

| 参数 | 说明 |
|------|------|
| `-i` / `--input` | VCF/TSV 目录 |
| `-o` / `--output` | 输出目录 |
| `-d` / `--min-depth` | 最低深度 (50) |
| `-f` / `--min-af` | 最低等位基因频率 (0.05) |

**产出 9 图:** mutation landscape, Ts/Tv 饼图, 蛋白功能饼图, 聚类热图, AF 频谱, 变异密度+基因轨道, AF 提琴图, 群体遗传滑窗 (π + Tajima's D), PCA 谱系聚类

### virus_vcf_pipeline.py — VCF 合并+PCA+距离

| 参数 | 说明 |
|------|------|
| `-d` / `--dir` | VCF 目录 |
| `-p` / `--pattern` | `**/*.filtered.vcf` |
| `-o` / `--out_dir` | 输出目录 |

**流程:** bcftools merge → biallelic SNP 过滤 → VCF2PCACluster (PCA) → VCF2Dis (距离矩阵)

### snpeff_analysis.py — SnpEff 宏观分析

| 参数 | 说明 |
|------|------|
| `--miner` | SnpEff VCF 目录 |
| `--outdir` | 输出目录 |

**产出:** 6 CSV 矩阵 + 4 Python 图 + 1 R OncoPrint

### snpeff2maf.py + viral_maftools.R — MAF 转换+瀑布图

- `snpeff2maf.py`: VCF → MAF 格式转换
- `viral_maftools.R`: oncoplot 瀑布图 + lollipop 蛋白结构图 + TiTv 汇总

### snpgenie_master.py — 11 图进化分析

| 参数 | 说明 |
|------|------|
| `-i` / `--input` | SNPGenie 输出目录 |
| `-o` / `--output` | 输出目录 |

**产出:** VAF 频谱, iSNV 密度, dN/dS jointplot, πN/πS jointplot, 单基因 dN/dS (Kruskal-Wallis), bootstrap CI, 突变频谱, 无义突变热点, 双轨滑窗, Trinity 全景, 2D/3D PCA+KMeans

---

## 7. capheine_pipeline.py — 正选择分析 (Stage 7)

**用途:** HyPhy 正选择管线。(1) 去终止密码子, (2) 按基因拆分, (3) cawlign 密码子感知比对, (4) 过滤歧义序列, (5) HyPhy CLN 去重, (6) IQ-TREE 建树, (7) 标记前景/参考分支, (8) HyPhy FEL/MEME/PRIME/BUSTED/CONTRASTFEL/RELAX, (9) DRHIP 聚合, (10) MultiQC 报告。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--reference` / `-r` | 必需 | 参考 CDS FASTA |
| `--unaligned` / `-u` | 必需 | 未比对序列 FASTA |
| `--outdir` / `-o` | 必需 | 输出目录 |
| `--foreground_list` | — | 前景 taxa 列表 |
| `--code` | `1` | 遗传密码 |
| `--workers` | 2 | 并行 worker 数 |
| `--cpus_iqtree` | 6 | IQ-TREE 线程 |
| `--cpus_hyphy` | 16 | HyPhy 线程 |

### 输出
```
{outdir}/
  cawlign/           ← 密码子比对
  hyphy/FEL/MEME/PRIME/BUSTED/  ← 各选择检验
  iqtree/            ← 系统发育树
  drhip/combined_sites.csv  ← 正选择位点汇总
  multiqc/           ← MultiQC 报告
```

**工具依赖:** cawlign, hyphy, iqtree, drhip, multiqc

---

## 8. virus_auto_pipeline.py — 相似性全景 (Stage 8)

**用途:** 两阶段相似性分析: Phase 1 从参考 GenBank 提取基因区, Phase 2 构建两两相似性矩阵 (NT/AA), 层次聚类, SDT 热图, CD-HIT 去重。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input` | 必需 | 共识 FASTA 文件/目录 |
| `-g` / `--genbank` | — | GenBank 参考 |
| `-o` / `--output_dir` | `./pipeline_results` | 输出 |
| `--mode` | `filter` | strict/filter/fill/all |
| `--align_method` | `pairwise` | pairwise/mafft |
| `--cdhit` | False | 启用 CD-HIT 去重 |
| `--max_n_genome` | 5.0 | 最大基因组 N% |
| `--threads` | 4 | 线程 |

### 输出
每 mode 子目录含:
- `Full_Dataset/02_similarity_matrices/` — 相似性矩阵 CSV
- `Full_Dataset/overall_*` / `gene_*` — 热图 + 分布图
- `Dedup_Dataset/` — 去重后分析

**工具依赖:** cd-hit-est, mafft (可选), NCBI Entrez

---

## 9. batch_virema_dvg.py — DVG/重组检测 (Stage 9)

**用途:** ViReMa DVG 检测。PE reads 合并 → FASTQ 清理 → ViReMa 运行 → 收集 SAM/BED/BEDPE → R 报告生成。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-s` / `--summary` | 必需 | 检出 summary TSV |
| `-r` / `--reference` | 必需 | 参考 FASTA |
| `-d` / `--reads_dir` | 必需 | Reads 目录 |
| `-v` / `--virema_script` | 必需 | ViReMa.py 路径 |
| `--seed` | 25 | ViReMa seed 长度 |
| `--mindel` | 15 | 最小 deletion |
| `--min_cov` | 0.0 | 覆盖度过滤 |
| `-j` / `--jobs` | 4 | 并行 |
| `--shm` | False | 使用 /dev/shm |

### 输出
```
{output}/
  virema_results/          ← 每病毒 ViReMa 输出
  virema_sams/             ← SAM 文件
  all_gathered_beds/       ← BED/BEDPE 汇总
  Summary_Analysis_Report/ ← R 生成图 (scatter + Circos)
```

**工具依赖:** bbmerge.sh, bowtie2-build, python (ViReMa.py), Rscript

---

## 10. generate_pipeline_report.py — 交互式报告 (Stage 10)

**用途:** 扫描管线全输出，生成自包含 HTML 交互式报告。含侧边栏病毒导航、嵌入图表、数据表、AI 写作建议。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-d` / `--dir` | 必需 | 管线输出根目录 |
| `-o` / `--output` | — | 输出 HTML 路径 |
| `--no_images` | False | 跳过图片嵌入 |

### 输入
管线输出根目录 (自动定位各阶段产物)

### 输出
`Pipeline_Summary_Report.html` (自包含，可离线查看)

---

## 11. 辅助脚本

| 脚本 | 位置 | 用途 |
|------|------|------|
| `utils/filter_summary.py` | Stage 2 | 高置信度过滤 (覆盖度/深度/reads/TPM/Poisson) |
| `utils/extract_full_fasta.py` | Stage 5 | 提最长 contig |
| `utils/gbk_extractor.py` | Stage 7 | GenBank CDS 提取 → capheine 输入 |
| `utils/visual_codon_miner.py` | Stage 7 | DRHIP 正选择位点密码子频率图 |
| `utils/consensus_extract.py` | 辅助 | NCBI 共识序列提取 |
| `virema_summary_report.R` | Stage 9 | DVG 可视化 (scatter + Circos 4-track) |
| `snpeff_build.py` | Stage 3 内部 | SnpEff 数据库构建 |
