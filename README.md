# MMPV-RNA v2.3

**Meta-transcriptomic Mining of Plant Virome from RNA-seq data**

> 宏转录组植物病毒端到端发现管道 — End-to-end plant virome discovery pipeline from raw RNA-seq reads to HQ vOTU catalog, including known virus quantification, novel virus discovery, host prediction, and cross-validation.

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![BioConda](https://img.shields.io/badge/bioconda-supported-brightgreen.svg)](https://bioconda.github.io/)

---

## Architecture / 架构

```
Raw FASTQ
    │
    ▼ clean          Fastp + Seqkit + Clumpify
00a_CleanData/       (QC → FASTA conversion → dedup)
    │
    ▼ deplete        Kraken2 + Bowtie2 + rRNA (Ribodetector / SILVA Bowtie2)
00b_HostDepletion/   (host classification → alignment → rRNA removal)
    │
    ├──────────────────────────┐
    ▼                          ▼
assembly                   auto_known_virus
Penguin + SPAdes + MEGAHIT   ├─ detect  Salmon/Bowtie2 rapid detection
01_Assembly/                 ├─ variants FreeBayes+SnpEff+SnpGenie
    │                        └─ full    de novo full-length assembly
    ▼
identification               known_viruses/
Genomad+Diamond+VirSorter2+
ViralVerify+VirHunter+Metabuli
02_Identification/
    │
    ▼
cobra                         BWA-MEM2 + COBRA + CheckV
03_COBRA/                     (contig overlap-based re-assembly)
    │
    ▼
cluster                       CD-HIT reference-guided + vclust Leiden
04_CLUSTER/                   centroids + clusters + split_fastas
    │
    ├──→ taxonomy             9-tool classification + R consensus
    │    05_Taxonomy/
    │
    ├──→ host                 ICTV > RNAVirHost > PhaBOX2 decision tree
    │    06_HostPrediction/
    │
    ├──→ checkv               per-host CheckV pre-evaluation
    │    07_Checkv/
    │
    └──→ rescue               host-filtered 3-branch cascade rescue
         08_Rescue/
              │
              ▼
         validate              taxonomy-level novelty judgment
         09_Validation/        ★known / ★★novel / ★★★truly
```

---

## Quick Start / 快速开始

```bash
# Clone repository / 克隆仓库
git clone https://github.com/zhangwenda0518/MMPV-RNA.git
cd MMPV-RNA

# Install dependencies / 安装依赖
conda create -n mmpr-rna python=3.10
conda activate mmpr-rna
conda install -c bioconda fastp seqkit bowtie2 samtools kraken2 blast checkv salmon spades megahit vclust
pip install polars biopython pandas tqdm

# Prepare reference databases / 准备参考数据库
mkdir -p db/
# ├── kraken2_db/          Kraken2 host database / Kraken2 宿主库
# ├── host_align_db/       Bowtie2 host alignment index / Bowtie2 宿主比对索引
# ├── virus_db/            Virus identification/classification database / 病毒鉴定/分类数据库
# ├── checkv_db/           CheckV database / CheckV 数据库
# ├── ref.fasta            Known virus reference sequences / 已知病毒参考序列
# ├── ref_info.tsv         Reference metadata (Accession, Taxid, Species) / 参考元数据
# ├── ICTV_plant_viruses.fasta  ICTV/NCBI complete plant virus genomes / 植物病毒完整基因组
# ├── phabox_db_v2_2/      PhaBOX2 database / PhaBOX2 数据库
# ├── cross_analysis/      ICTV host probability tables / ICTV 宿主概率表
# └── silva_index/         SILVA Bowtie2 index (optional) / SILVA Bowtie2 索引 (可选)

# One-shot full pipeline / 一键运行
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

---

## Pipeline Stages / 流水线阶段 (11 stages)

| # | Stage / 阶段 | Script / 脚本 | Core Function / 核心功能 |
|---|------|------|----------|
| 1 | `clean` | clean-data.py | Fastp QC + Seqkit stats + Clumpify dedup |
| 2 | `deplete` | host_depletion.py | Kraken2 + Bowtie2/HISAT2 host depletion + rRNA removal |
| 3 | `assembly` | assembly_pipeline.py | 3-tool assembly (Penguin/MEGAHIT/rnaviralSPAdes) |
| 4 | `identification` | virus_identification.py | 6-tool parallel virus identification |
| 5 | `cobra` | cobra_pipeline.py | COBRA batch extension |
| 6 | `cluster` | cluster_pipeline.py | CD-HIT ref-guided + vclust Leiden clustering |
| 7 | `taxonomy` | virus_classifier.py + R | 9-tool classification + R consensus |
| 8 | `host` | run_host_prediction.py | ICTV > RNAVirHost > PhaBOX2 host prediction |
| 9 | `checkv` | (built-in) | Per-host CheckV completeness pre-evaluation |
| 10 | `rescue` | rescue_pipeline.py + Virseqimprover.py | 3-branch cascade rescue (A:CheckV → C:VSI → D:BLASTN+VSI) |
| 11 | — | validate_novel_viruses.py | Taxonomy-level novelty judgment + virus frequency stats |

### Step-by-step / 分步运行

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

Each stage supports `--stage X --help` for detailed parameters and sub-script CLI.
/ 每个阶段支持 `--stage X --help` 查看详细参数和子脚本 CLI。

---

## Known Virus Analysis / 已知病毒分析 (independent parallel track)

```bash
# Runs in parallel with main pipeline, sharing 00b_HostDepletion/
# 与主管道并行运行, 共用 00b_HostDepletion/

# Rapid detection / 快速检测
python auto_known_virus.py --stage detect \
    --reads_dir /data/out/00b_HostDepletion/ \
    --output_dir /data/out/known_viruses/ \
    --ref_info /db/ref_info.tsv --reference /db/ref.fasta \
    --tool salmon --threads 40 --jobs 4

# Variant analysis / 变异分析
python auto_known_virus.py --stage variants \
    --reads_dir /data/out/00b_HostDepletion/ \
    --output_dir /data/out/known_viruses/ \
    --ref_info /db/ref_info.tsv --reference /db/ref.fasta \
    --variant_caller freebayes --snpeff --snpgenie \
    --threads 40 --jobs 4

# Full-length assembly / 全长组装
python auto_known_virus.py --stage full \
    --reads_dir /data/out/00b_HostDepletion/ \
    --output_dir /data/out/known_viruses/ \
    -j 4 -t 40
```

---

## Cross-validation / 交叉验证

```bash
python validate_novel_viruses.py \
    -i /data/out/08_Rescue/Plant/centroids/final_centroids.fasta \
    --taxonomy /data/out/05_Taxonomy/integrated/final_integrated_classification.tsv \
    --cdhit-known /data/out/04_CLUSTER/centroids/known_association.tsv \
    --clusters-tsv /data/out/04_CLUSTER/3_vclust/vclust_clusters.tsv \
    --host /data/out/06_HostPrediction/ensemble_host_summary.tsv \
    -o /data/out/09_Validation/
```

**Classification rules / 分类规则:**

```
Species ≠ NA              → ★ known          已知病毒
Genus ≠ NA, Species = NA  → ★★ novel_species 新种
Family ≠ NA, Genus = NA   → ★★ novel_genus   新属
Order/Class ≠ NA          → ★★★ novel_family 新科
All NA                    → ★★★ truly_novel  全新
```

---

## Three-Branch Cascade Rescue / 三支路级联拯救

```
centroids
    │
    ▼
Branch A: CheckV parallel evaluation → completeness ≥ 90% → pass
    │ fail
    ▼
Branch C: Virseqimprover iterative extension (multi-sample read aggregation)
    Salmon quantification → BBMap extraction → SPAdes assembly → CheckV
    │ fail
    ▼
Branch D: BLASTN reference search + CheckV + VSI final rescue
    │
    ▼
Merge + vclust final dedup → HQ vOTU

Skip rescue = CD-HIT known + CheckV pass (≥90%)
免拯救 = CD-HIT 已知 + CheckV pass (≥90%)
```

---

## Key Features / 核心特性

- **CD-HIT reference-guided pre-clustering / CD-HIT 参考引导预聚类**: Incorporates ICTV/NCBI complete genomes as references; vclust cd-hit greedy incremental clustering associates fragmented contigs with known species.
- **Three-branch cascade rescue / 三支路级联拯救**: CheckV → Virseqimprover (multi-sample reads aggregation) → BLASTN+VSI; progressively boosts HQ vOTU yield.
- **Host filtering / 宿主过滤**: Pre-filter by host category before rescue, saving 70%+ compute.
- **SILVA rRNA removal / SILVA rRNA 去除**: Supports Bowtie2 alignment to SILVA database for rRNA depletion (alternative to Ribodetector).
- **Checkpoint-resume / 断点续传**: All scripts support `--resume/--force`; safe interruption and recovery for production-scale runs.
- **Taxonomy-level novelty judgment / 基于分类层级的新颖性判断**: No BLASTN dependency; determines novel viruses directly from taxonomy completeness.
- **Orchestrator self-control / 编排器自身控制**: `--dry-run`, `--log-level`, `--stop-on-error` for both `virome_pipeline.py` and `auto_known_virus.py`.

---

## Final Output / 最终产出

```
out/
├── 00a_CleanData/                        clean: cleaned reads / 清洗后 reads
├── 00b_HostDepletion/                    deplete: host-depleted clean reads / 去宿主清洁 reads
├── 01_Assembly/{sample}/                 assembly: 3-tool contigs / 三工具 contig
├── 02_Identification/{sample}/           identification: candidate viruses / 候选病毒
├── 03_COBRA/{sample}/                    cobra: COBRA extensions / COBRA 延伸
├── 04_CLUSTER/                           cluster: clustering + centroids / 聚类 + centroids
│   ├── 2_cdhit/known_association.tsv     contig→ref mapping / contig→参考映射
│   ├── 3_vclust/vclust_clusters.tsv      cluster results / 聚类结果
│   ├── 3_vclust/split_fastas/            per-cluster splits / per-cluster 拆分
│   └── centroids/final_centroids.fasta   all representative sequences / 全部代表序列
├── 05_Taxonomy/                          taxonomy: classification / 分类注释
│   └── integrated/final_integrated_classification.tsv
├── 06_HostPrediction/                    host: host prediction / 宿主预测
│   ├── ensemble_host_summary.tsv
│   └── host_classified_fasta/
├── 07_Checkv/                            checkv: per-host pre-evaluation / 按宿主预评估
│   └── checkv_pass_ids.txt
├── 08_Rescue/                            rescue: 3-branch rescue / 三支路拯救
│   ├── Plant/centroids/final_centroids.fasta  ★ HQ vOTU
│   └── checkv/                           CheckV quality reports / CheckV 质量报告
├── known_viruses/                        auto_known_virus / 已知病毒
│   ├── 1_FastViromeExplorer/
│   ├── 2_Virus_variants_Results/
│   └── 3_Virus_assemblies_final/
└── 09_Validation/                        validate_novel_viruses / 新颖性验证
    ├── novel_viruses.annotated.tsv
    ├── final_virus_catalog.fasta
    └── validation_report.html            Plotly.js interactive report / 交互报告
```

---

## Parameter Reference / 参数速查

| Param / 参数 | Default / 默认值 | Description / 说明 |
|------|------|------|
| `--input_reads` | — | Raw FASTQ directory / 原始 FASTQ 目录 |
| `--output_dir` | **required / 必需** | Project output root / 项目输出根目录 |
| `--stage` | `all` | Run stage(s) / 运行阶段 |
| `--dry-run` | off | Scan samples & show config only / 仅扫描并显示配置 |
| `--log-level` | `INFO` | Log verbosity: DEBUG/INFO/WARNING/ERROR |
| `--stop-on-error` | off | Abort on sub-script failure / 子脚本失败时立即终止 |
| `--kraken2_db` | — | Kraken2 host database / Kraken2 宿主库 |
| `--host_align_db` | — | Host alignment index / 宿主比对索引 |
| `--virus_db` | — | Virus database root / 病毒数据库根目录 |
| `--checkv_db` | — | CheckV database / CheckV 数据库 |
| `--ref-genomes` | — | ICTV/NCBI reference genomes / 参考基因组 |
| `--phabox-db` | — | PhaBOX2 database / PhaBOX2 数据库 |
| `--prob-dir` | — | ICTV host probability tables / ICTV 宿主概率表 |
| `--host-filter` | `Plant` | Target host(s) / 目标宿主 |
| `--rrna` | off | Enable rRNA removal / 开启 rRNA 剔除 |
| `--rrna_tool` | `ribodetector` | rRNA tool: ribodetector/silva |
| `--silva_index` | — | SILVA Bowtie2 index / SILVA 索引 |
| `--assembler` | `penguin` | Assembly tool / 组装工具 |
| `--aligner` | `bowtie2` | Host depletion aligner / 去宿主比对工具 |
| `--seq_type` | `rna-short` | Sequencing type / 测序类型 |
| `-t, --threads` | `20` | Threads per task / 单任务线程 |
| `-m, --memory` | `64` | Memory GB / 内存 GB |
| `-j, --jobs` | `2` | Parallel tasks / 并行任务数 |
| `--min-length` | `500` | Min virus length bp / 病毒最小长度 |
| `--ani` | `0.95` | Clustering ANI threshold / 聚类 ANI 阈值 |
| `--qcov` | `0.85` | Clustering QCOV threshold / 聚类 QCOV 阈值 |
| `--force` | off | Force re-run / 强制重跑 |

---

## Documentation / 文档索引

| Document / 文档 | Description / 说明 |
|------|------|
| [doc/00-scripts-interpretation.md](doc/00-scripts-interpretation.md) | Full 20-script interpretation / 全脚本解读 |
| [doc/14-pipeline-diagrams.md](doc/14-pipeline-diagrams.md) | 8 Mermaid diagrams (flowcharts/architecture/sequence) / 流程框架图 |
| [doc/01-virome_pipeline.md](doc/01-virome_pipeline.md) | Orchestrator details / 编排器详解 |
| [doc/02-clean-data.md](doc/02-clean-data.md) | Clean stage / 数据清洗 |
| [doc/03-host_depletion.md](doc/03-host_depletion.md) | Deplete stage / 去宿主 |
| [doc/04-assembly_pipeline.md](doc/04-assembly_pipeline.md) | Assembly stage / 组装 |
| [doc/05-virus_identification.md](doc/05-virus_identification.md) | Identification stage / 病毒鉴定 |
| [doc/06-cobra_pipeline.md](doc/06-cobra_pipeline.md) | COBRA stage / 延伸 |
| [doc/07-cluster_pipeline.md](doc/07-cluster_pipeline.md) | Cluster stage / 聚类 |
| [doc/08-rescue_pipeline.md](doc/08-rescue_pipeline.md) | Rescue stage / 拯救 |
| [doc/09-Virseqimprover.md](doc/09-Virseqimprover.md) | Virseqimprover / 迭代延伸 |
| [doc/10-taxonomy-host.md](doc/10-taxonomy-host.md) | Taxonomy + Host / 分类+宿主 |
| [doc/11-auto_known_virus.md](doc/11-auto_known_virus.md) | Known virus pipeline / 已知病毒 |
| [doc/12-validate.md](doc/12-validate.md) | Validation / 验证 |
| [doc/13-utilities.md](doc/13-utilities.md) | Utility scripts / 工具脚本 |

---

## Citation / 引用

Zhang W. et al. **MMPV-RNA: Meta-transcriptomic Mining of Plant Virome from RNA-seq data**. 2026.

> If you use MMPV-RNA in your research, please cite the above reference.
> 如果您在研究中使用了 MMPV-RNA，请引用上述文献。
