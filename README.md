# MMPV-RNA v3.0

**Macro-metavirome Plant Virus Discovery & Analysis Platform**

> 端到端植物病毒组发现与分析平台 — 从原始 RNA-seq reads 到 HQ vOTU catalog 的完整闭环，含公共数据挖掘、新病毒发现、已知病毒定量/变异/进化深度分析。

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.2-blue.svg)](https://www.r-project.org/)
[![pixi](https://img.shields.io/badge/pixi-enabled-green.svg)](https://pixi.sh/)
[![BioConda](https://img.shields.io/badge/bioconda-supported-brightgreen.svg)](https://bioconda.github.io/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## 三大管线 / Three Pipelines

```
public_metadata_pipeline/     → 公共数据获取
       │
       ├── SRA/GSA 搜索 → 元数据清洗 → 批量下载 → 可视化
       └── 宿主参考基因组下载 → Kraken2/Bowtie2/HISAT2/Minimap2 索引
              │
              ▼
virome_discovery_pipeline/    → 新病毒发现
       │
       ├── 清洁 → 去宿主 → 组装 → 10工具并行鉴定 → COBRA延伸
       ├── CD-HIT参考聚类 → vclust去冗余 → 8工具分类 → 宿主预测
       └── CheckV评估 → 三支路级联拯救 → HQ vOTU catalog
              │
              ▼
virome_analysis_pipeline/     → 已知病毒深度分析
       │
       ├── 快速定量 (Salmon/Bowtie2) → 变异检测 (FreeBayes/iVar)
       ├── SnpEff注释 → SNPGenie进化 → 12步全长组装
       ├── HyPhy正选择 → 相似性全景 → DVG/重组检测
       └── 交互式HTML综合报告
```

---

## 快速开始 / Quick Start

### 1. 一键部署

```bash
git clone https://github.com/zhangwenda0518/MMPV-RNA.git
cd MMPV-RNA

# 安装全部 100 个依赖 (pixi)
pixi install

# 安装 ViraLM 独立环境 (可选)
conda env create -f envs/viralm.yaml -n viralm

# 下载参考数据库 (见下方数据库清单)
```

### 2. 发现管线

```bash
pixi run discovery-downstream \
    --output_dir /data/out/ \
    --input_reads /data/00b_HostDepletion/ \
    --host_db /db/hostdb/ \
    --virus_db /db/virus-db/ \
    --checkv_db /db/checkv-db-v1.7/ \
    --host-filter Plant \
    --coassembly \
    -t 120 -j 20
```

### 3. 分析管线

```bash
# 快速定量
python virome_analysis_pipeline/auto_known_virus.py --stage detect \
    --reads_dir /data/out/00b_HostDepletion/ \
    --reference /db/ref.fasta --ref_info /db/ref_info.tsv \
    --tool salmon -t 40 -j 4

# 全流程
python virome_analysis_pipeline/auto_known_virus.py --stage all \
    --reads_dir /data/out/00b_HostDepletion/ \
    --reference /db/ref.fasta --ref_info /db/ref_info.tsv \
    --snpeff --snpgenie -t 40 -j 4
```

### 4. 公共数据管线

```bash
# 搜索某物种的公共 SRA/GSA 数据
python public_metadata_pipeline/public_data_pipeline.py \
    --species "Solanum lycopersicum" --taxid 4081 \
    --stage search info down plot

# 构建宿主参考基因组索引
python public_metadata_pipeline/build_host_pipeline.py \
    --species "Solanum lycopersicum" --taxid 4081 \
    --stage all --threads 30
```

---

## 发现管线阶段 / Discovery Pipeline Stages

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 0a | `clean` | clean-data.py | Fastp QC + Seqkit FASTA转换 + Clumpify去重 |
| 0b | `deplete` | host_depletion.py | Kraken2 → Bowtie2/HISAT2/Minimap2 → rRNA去除 |
| 0c | `bbnorm` | run_bbnorm.py | BBNorm覆盖度归一化 (共组装前可选) |
| 1 | `assembly` | assembly_pipeline.py | MEGAHIT / rnaviralSPAdes / Penguin 组装 |
| 2 | `identification` | virus_identification.py | 10工具并行鉴定 (Genomad+Diamond+VirSorter2+VirHunter+Metabuli+...) |
| 3 | `cobra` | cobra_pipeline.py | BWA-MEM2 → CoverM → COBRA-Meta 延伸 |
| 4 | `cluster` | cluster_pipeline.py | CD-HIT参考引导预聚类 + vclust Leiden聚类 |
| 5 | `taxonomy` | virus_classifier.py + R | 8工具分类 + 加权投票共识 (8级 taxonomy) |
| 6 | `host` | run_host_prediction.py | ICTV > RNAVirHost > PhaBOX2 决策树宿主预测 |
| 7 | `checkv` | (内置) | CheckV 完整性评估 |
| 8 | `rescue` | rescue_pipeline.py | 三支路级联拯救 (CheckV → Virseqimprover → BLASTN) |
| 9 | `report` | report_pipeline.py | TSV汇总 + Sankey图 + 交互式HTML报告 |

---

## 分析管线阶段 / Analysis Pipeline Stages

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `detect` | batch_virus_depth.py | Salmon/Kallisto/Bowtie2 快速定量 + Poisson过滤 |
| 2 | `filter` | utils/filter_summary.py | 高置信度过滤 |
| 3 | `variants` | batch_virus_variants.py | FreeBayes/iVar/LoFreq + SnpEff + SNPGenie |
| 4 | `full` | batch_virus_full.py → virus-full.py | 12步全长组装 |
| 5 | `extract` | utils/extract_full_fasta.py | 最长contig提取 |
| 6 | `post` | 6脚本并行 | VCF可视化 + PCA + MAF + SnpGenie分析 |
| 7 | `capheine` | capheine_pipeline.py | HyPhy FEL/MEME/BUSTED/PRIME正选择 |
| 8 | `similarity` | virus_auto_pipeline.py | 全长相似性热图 + 层次聚类 |
| 9 | `dvg` | batch_virema_dvg.py | ViReMa DVG/重组检测 + Circos图 |
| 10 | `report` | generate_pipeline_report.py | 交互式HTML综合报告 |

---

## 目录结构 / Directory Structure

```
MMPV-RNA/
├── virome_discovery_pipeline/      # 新病毒发现 (18+4脚本)
│   ├── virome_pipeline.py          # 主编排器
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 辅助工具
│
├── virome_analysis_pipeline/       # 已知病毒深度分析 (17+11脚本)
│   ├── auto_known_virus.py         # 分析编排器
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 辅助工具
│
├── public_metadata_pipeline/       # 公共数据获取 (8+2脚本)
│   ├── public_data_pipeline.py     # 公共数据编排器
│   ├── build_host_pipeline.py      # 宿主库构建编排器
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 共享工具模块
│
├── biosoft/                        # 第三方工具 (脚本/JAR, 无需编译)
│   ├── VirBot/VirBot.py
│   ├── virhunter/predict_cpu.py + weights/
│   ├── ViReMa/ViReMa.py
│   └── snpEff/snpEff.jar + config + scripts/
│
├── stats/                          # 辅助统计/可视化
├── doc/                            # 历史文档 (17篇)
├── envs/viralm.yaml                # ViraLM 独立conda环境
├── pixi.toml                       # 一键部署 (100个依赖)
├── pipeline_config.yaml            # 管线配置文件
├── SOFTWARE_VERSIONS.txt           # 全部第三方软件版本
└── README.md
```

---

## 核心特性 / Key Features

- **10工具并行病毒鉴定**: Genomad + Diamond BLASTX + VirSorter2 + ViralVerify + VirHunter + Metabuli + RdrpCatch + ViraLM + VirBot + Viroid BLASTN
- **CD-HIT参考引导预聚类**: 整合ICTV/NCBI完整基因组为参考，关联碎片contig到已知物种
- **三支路级联拯救**: CheckV → Virseqimprover (多样本reads聚合) → BLASTN+VSI，逐步提升HQ vOTU产出
- **宿主过滤**: 拯救前按宿主类别预过滤，节省70%+计算量
- **pixi一键部署**: 100个conda包精确版本管理，`pixi install`即可
- **断点续传**: 全部脚本支持 `--resume/--force`，大规模运行安全中止和恢复
- **基于分类层级的新颖性判断**: 无需BLASTN依赖，直接从taxonomy completeness判断新病毒
- **完整闭环**: 从公共数据挖掘到投稿图表，一站式完成

---

## 数据库清单 / Required Databases

| 数据库 | 用途 | 来源 |
|--------|------|------|
| CheckV DB | 病毒完整性评估 | https://bitbucket.org/berkeleylab/checkv/ |
| geNomad DB | 病毒鉴定 | https://zenodo.org/records/14828026 |
| RVDB | 病毒参考序列 | https://fzer.github.io/rvdbtools/ |
| NR (diamond) | 蛋白过滤 | NCBI nr |
| NCBI virus ref | BLAST参考 | NCBI virus |
| ViralVerify HMM | HMM验证 | ViralVerify |
| VirSorter2 DB | 病毒分类 | VirSorter2 |
| ViraLM DB | DNABERT-2鉴定 | Google Drive (gdown) |
| Kraken2 + align DB | 宿主去除 | 本管线构建 (build_host_pipeline.py) |
| PhaBOX2 DB | 宿主预测 | PhaBOX2 |

---

## 文档 / Documentation

| 文档 | 说明 |
|------|------|
| `virome_discovery_pipeline/doc.md` | 发现管线全部脚本 — 参数/输入输出/结果解读 |
| `virome_analysis_pipeline/doc.md` | 分析管线全部脚本 — 参数/输入输出/结果解读 |
| `public_metadata_pipeline/doc.md` | 公共数据管线 + 宿主库构建 |
| `doc/14-pipeline-diagrams.md` | Mermaid流程图/架构图 |
| `SOFTWARE_VERSIONS.txt` | 全部第三方软件版本记录 |
| `pixi.toml` | 依赖配置 (可直接查看) |

---

## 引用 / Citation

Zhang W. et al. **MMPV-RNA: Macro-metavirome Plant Virus Discovery and Analysis Platform**. 2026.

> If you use MMPV-RNA in your research, please cite the above reference.
