# MMPV-RNA v2.3

**Meta-transcriptomic Mining of Plant Virome from RNA-seq data**

宏转录组植物病毒端到端发现管道 — 从原始 RNA-seq reads 到 HQ vOTU 目录，含已知病毒定量、新病毒发现、宿主预测、交叉验证。

---

## 架构

```
Raw FASTQ
    │
    ▼ clean          Fastp + Seqkit + Clumpify
00a_CleanData/
    │
    ▼ deplete        Kraken2 + Bowtie2 + rRNA (Ribodetector / SILVA Bowtie2)
00b_HostDepletion/
    │
    ├──────────────────────────┐
    ▼                          ▼
assembly                   auto_known_virus
Penguin + SPAdes + MEGAHIT   ├─ detect  Salmon/Bowtie2 快速检测
01_Assembly/                 ├─ variants FreeBayes+SnpEff+SnpGenie
    │                        └─ full    全长组装
    ▼
identification               known_viruses/
Genomad+Diamond+VirSorter2+
ViralVerify+VirHunter+Metabuli
02_Identification/
    │
    ▼
cobra                         BWA-MEM2 + COBRA + CheckV
03_COBRA/
    │
    ▼
cluster                       CD-HIT 参考引导 + vclust Leiden
04_CLUSTER/                   centroids + clusters + split_fastas
    │
    ├──→ taxonomy             5工具分类 + R共识
    │    05_Taxonomy/
    │
    ├──→ host                 ICTV > RNAVirHost > PhaBOX2 决策树
    │    06_HostPrediction/
    │
    ├──→ checkv               按宿主 CheckV 预评估
    │    07_Checkv/
    │
    └──→ rescue               宿主过滤 + 三支路级联拯救
         08_Rescue/
              │
              ▼
         validate              分类层级判断 ★known/★★novel/★★★truly
         09_Validation/
```

## 快速开始

```bash
# 克隆仓库
git clone https://github.com/zhangwenda0518/MMPV-RNA.git
cd MMPV-RNA

# 安装依赖
conda create -n mmpr-rna python=3.10
conda activate mmpr-rna
conda install -c bioconda fastp seqkit bowtie2 samtools kraken2 blast checkv salmon spades megahit vclust
pip install polars biopython pandas tqdm

# 准备参考数据库
mkdir -p db/
# ├── kraken2_db/          Kraken2 宿主库
# ├── host_align_db/       Bowtie2 宿主比对索引
# ├── virus_db/            病毒鉴定/分类数据库
# ├── checkv_db/           CheckV 数据库
# ├── ref.fasta            已知病毒参考序列
# ├── ref_info.tsv         参考元数据 (Accession, Taxid, Species)
# ├── ICTV_plant_viruses.fasta  ICTV/NCBI 植物病毒完整基因组
# ├── phabox_db_v2_2/      PhaBOX2 数据库
# ├── cross_analysis/      ICTV 宿主概率表
# └── silva_index/         SILVA Bowtie2 索引 (可选)

# 一键运行
python virome_pipeline.py --stage all \
    --input_reads /data/raw/ \
    --output_dir /data/out/ \
    --kraken2_db    /db/kraken2_db/ \
    --host_align_db /db/host_align_db/ \
    --virus_db      /db/virus_db/ \
    --checkv_db     /db/checkv_db/ \
    --phabox-db     /db/phabox_db_v2_2/ \
    --prob-dir      /db/cross_analysis/ \
    --ref-genomes   /db/ICTV_plant_viruses.fasta \
    --assembler all --aligner bowtie2 --seq_type rna-short \
    --rrna --host-filter Plant \
    --min-length 500 --ani 0.95 --qcov 0.85 \
    -t 120 -j 20 -m 256
```

## 11 个阶段

| # | 阶段 | 脚本 | 核心功能 |
|---|------|------|----------|
| 1 | `clean` | clean-data.py | Fastp 质控 + Seqkit 统计 + Clumpify 去重 |
| 2 | `deplete` | host_depletion.py | Kraken2 分类 + Bowtie2/HISAT2 去宿主 + rRNA 去除 |
| 3 | `assembly` | assembly_pipeline.py | Penguin / MEGAHIT / rnaviralSPAdes 三工具组装 |
| 4 | `identification` | virus_identification16.py | 6工具并行病毒鉴定 |
| 5 | `cobra` | cobra_pipeline.py | COBRA 批量延伸 |
| 6 | `cluster` | cluster_pipeline.py | CD-HIT 参考引导 + vclust Leiden 聚类 |
| 7 | `taxonomy` | virus_classifier2.py + R | 5工具分类 + R 共识整合 |
| 8 | `host` | run_host_prediction.py | ICTV>RNAVirHost>PhaBOX2 宿主预测 |
| 9 | `checkv` | (内置) | 按宿主 CheckV 完整性预评估 |
| 10 | `rescue` | rescue_pipeline.py + Virseqimprover.py | 三支路级联拯救 (A:CheckV → C:VSI → D:BLASTN+VSI) |
| 11 | — | validate_novel_viruses.py | 基于分类层级的新颖性判断 + 病毒频率统计 |

### 分步运行

```bash
python virome_pipeline.py --stage clean       --help
python virome_pipeline.py --stage deplete     --help
python virome_pipeline.py --stage assembly    --help
python virome_pipeline.py --stage identification --help
python virome_pipeline.py --stage cobra       --help
python virome_pipeline.py --stage cluster     --help
python virome_pipeline.py --stage taxonomy    --help
python virome_pipeline.py --stage host        --help
python virome_pipeline.py --stage checkv      --help
python virome_pipeline.py --stage rescue      --help
```

每个阶段支持 `--stage X --help` 查看详细参数和子脚本 CLI。

## 已知病毒分析 (独立并行)

```bash
# 与主管道并行运行, 共用 00b_HostDepletion/

# 快速检测
python auto_known_virus.py --stage detect \
    --reads_dir /data/out/00b_HostDepletion/ \
    --output_dir /data/out/known_viruses/ \
    --ref_info /db/ref_info.tsv --reference /db/ref.fasta \
    --tool salmon --threads 40 --jobs 4

# 变异分析
python auto_known_virus.py --stage variants \
    --reads_dir /data/out/00b_HostDepletion/ \
    --output_dir /data/out/known_viruses/ \
    --ref_info /db/ref_info.tsv --reference /db/ref.fasta \
    --variant_caller freebayes --snpeff --snpgenie \
    --threads 40 --jobs 4

# 全长组装
python auto_known_virus.py --stage full \
    --reads_dir /data/out/00b_HostDepletion/ \
    --output_dir /data/out/known_viruses/ \
    -j 4 -t 40
```

## 交叉验证

```bash
python validate_novel_viruses.py \
    -i /data/out/08_Rescue/Plant/centroids/final_centroids.fasta \
    --taxonomy /data/out/05_Taxonomy/integrated/final_integrated_classification.tsv \
    --cdhit-known /data/out/04_CLUSTER/centroids/known_association.tsv \
    --clusters-tsv /data/out/04_CLUSTER/3_vclust/vclust_clusters.tsv \
    --host /data/out/06_HostPrediction/ensemble_host_summary.tsv \
    -o /data/out/09_Validation/
```

**分类规则:**

```
Species ≠ NA              → ★ known          已知病毒
Genus ≠ NA, Species = NA  → ★★ novel_species 新种
Family ≠ NA, Genus = NA   → ★★ novel_genus   新属
Order/Class ≠ NA          → ★★★ novel_family 新科
全是 NA                    → ★★★ truly_novel  全新
```

## 三支路级联拯救

```
centroids
    │
    ▼
分支 A: CheckV 并行评估 → completeness ≥ 90% → pass
    │ fail
    ▼
分支 C: Virseqimprover reads 迭代延伸 (cluster 多样本 reads 聚合)
    Salmon 定量 → BBMap 提取 → SPAdes 组装 → CheckV
    │ fail
    ▼
分支 D: BLASTN 参考搜索 + CheckV + VSI 最后拯救
    │
    ▼
合并 + vclust 最终去重 → HQ vOTU

免拯救 = CD-HIT 已知 + CheckV pass (≥90%)
```

## 核心特性

- **CD-HIT 参考引导预聚类**: 引入 ICTV/NCBI 完整基因组作为参考，vclust cd-hit 贪婪增量聚类将碎片化 contig 关联到已知物种
- **三支路级联拯救**: CheckV → Virseqimprover (多样本 reads 聚合) → BLASTN+VSI，三支路渐进式提升 HQ vOTU 产出
- **宿主过滤**: 在拯救前按宿主分类过滤，节省 70%+ 计算量
- **SILVA rRNA 去除**: 支持 Bowtie2 比对 SILVA 库去除 rRNA (替代 Ribodetector)
- **断点续传**: 全部脚本支持 --resume/--force，大规模生产环境可安全中断恢复
- **基于分类层级的新颖性判断**: 不依赖 BLASTN，直接从 taxonomy 完整性判断新病毒

## 最终产出

```
out/
├── 00a_CleanData/                        clean: 清洗后 reads
├── 00b_HostDepletion/                    deplete: 去宿主清洁 reads
├── 01_Assembly/{sample}/                 assembly: 三工具 contig
├── 02_Identification/{sample}/           identification: 候选病毒
├── 03_COBRA/{sample}/                    cobra: COBRA 延伸
├── 04_CLUSTER/                           cluster: 聚类 + centroids
│   ├── 2_cdhit/known_association.tsv     contig→参考映射
│   ├── 3_vclust/vclust_clusters.tsv      聚类结果
│   ├── 3_vclust/split_fastas/            per-cluster 拆分
│   └── centroids/final_centroids.fasta   全部代表序列
├── 05_Taxonomy/                          taxonomy: 分类注释
│   └── integrated/final_integrated_classification.tsv
├── 06_HostPrediction/                    host: 宿主预测
│   ├── ensemble_host_summary.tsv
│   └── host_classified_fasta/
├── 07_Checkv/                            checkv: 按宿主预评估
│   └── checkv_pass_ids.txt
├── 08_Rescue/                            rescue: 三支路拯救
│   ├── Plant/centroids/final_centroids.fasta  ★ HQ vOTU
│   └── checkv/                           CheckV 质量报告
├── known_viruses/                        auto_known_virus
│   ├── 1_FastViromeExplorer/
│   ├── 2_Virus_variants_Results/
│   └── 3_Virus_assemblies_final/
└── 09_Validation/                        validate_novel_viruses
    ├── novel_viruses.annotated.tsv
    ├── final_virus_catalog.fasta
    └── validation_report.html            Plotly.js 交互报告
```

## 参数速查

| 参数 | 默认 | 说明 |
|------|------|------|
| `--input_reads` | — | 原始 FASTQ 目录 |
| `--output_dir` | **必需** | 项目输出根目录 |
| `--stage` | `all` | 运行阶段 |
| `--kraken2_db` | — | Kraken2 宿主库 |
| `--host_align_db` | — | 宿主比对索引 |
| `--virus_db` | — | 病毒数据库根目录 |
| `--checkv_db` | — | CheckV 数据库 |
| `--ref-genomes` | — | ICTV/NCBI 参考基因组 |
| `--phabox-db` | — | PhaBOX2 数据库 |
| `--prob-dir` | — | ICTV 宿主概率表 |
| `--host-filter` | `Plant` | 目标宿主 |
| `--rrna` | 关闭 | 开启 rRNA 剔除 |
| `--rrna_tool` | `ribodetector` | rRNA 工具 (ribodetector/silva) |
| `--silva_index` | — | SILVA Bowtie2 索引 |
| `--assembler` | `penguin` | 组装工具 |
| `--aligner` | `bowtie2` | 去宿主比对工具 |
| `--seq_type` | `rna-short` | 测序类型 |
| `-t, --threads` | `20` | 单任务线程 |
| `-m, --memory` | `64` | 内存 GB |
| `-j, --jobs` | `2` | 并行任务数 |
| `--min-length` | `500` | 病毒最小长度 bp |
| `--ani` | `0.95` | 聚类 ANI 阈值 |
| `--qcov` | `0.85` | 聚类 QCOV 阈值 |
| `--force` | 关闭 | 强制重跑 |

## 引用

Zhang W. et al. MMPV-RNA: Meta-transcriptomic Mining of Plant Virome from RNA-seq data. 2026.
