# MMPV-RNA v3.0

**Macro-metavirome Plant Virus Discovery & Analysis Platform**

> 端到端植物病毒组发现与分析平台 — 从原始 RNA-seq reads 到 HQ vOTU catalog 的完整闭环，含公共数据挖掘、新病毒发现、已知病毒定量/变异/进化深度分析。

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.2-blue.svg)](https://www.r-project.org/)
[![pixi](https://img.shields.io/badge/pixi-enabled-green.svg)](https://pixi.sh/)
[![BioConda](https://img.shields.io/badge/bioconda-supported-brightgreen.svg)](https://bioconda.github.io/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## 四大管线 / Four Pipelines

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
              │
              ▼
virome_submission_pipeline/   → 数据提交
       │
       ├── 拓扑判断 → 元数据模板 → 假定蛋白注释
       └── suvtk tbl2asn / Sequin tbl2asn 双模式 → .sqn 提交文件
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

### 5. 提交管线

```bash
# 从 08_Rescue 一键生成 NCBI 提交文件
python virome_submission_pipeline/submission_pipeline.py \
    --work-dir $OUT/08_Rescue/ \
    --run-title my_plant_virome \
    --mode both \
    --suvtk-db ~/database/virus-db/suvtk_db/ \
    -t 40

# 启动提交 GUI 桌面应用
python submission_gui/submission_gui.py
```

---

## 公共数据管线阶段 / Public Metadata Pipeline Stages

### SRA/GSA 公共数据获取

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `search` | gsa_sra.search.py | NCBI SRA + CNCB GSA 双引擎物种检索 → SRA_GSA_Merged_Final.csv |
| 2 | `info` | gsa_sra.info.py | SRA XML解析 + GSA爬虫 + AI元数据清洗 → Global_Unified_Metadata_Core13.csv (13列) |
| 3 | `down` | gsa_sra.down.py | aria2c/wget/prefetch 双协议下载 → 原始 FASTQ/SRA |
| 4 | `plot` | gsa_sra.plot.py | 时间/组织/地区/机构 SCI级6面板可视化 |

### 宿主参考数据库构建

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `genome-down` | download_host_genome.py | NCBI datasets 下载参考基因组 + GFF3 + 细胞器基因组 |
| 2 | `hostdb` | build_hostbase.py | Kraken2 + Bowtie2 + HISAT2 + Minimap2 四种索引构建 |

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

## 数据提交管线阶段 / Submission Pipeline Stages

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `topology` | viral_topology.py | BWA-MEM2 末端比对 → 病毒基因组拓扑判断 (circular/linear) |
| 2 | `metadata` | unified_metadata.py | 统一元数据模板生成 (source.src + features + organism) |
| 3 | `hypothetical` | analyze_hypothetical.py | 假定蛋白功能注释 (HHsuite/Diamond/DeepLoc/PSORTb/TMHMM 5工具) |
| 4 | `sequin` | sequin_builder.py | Sequin .tbl 格式构建 (Cenote-Taker3 风格) |
| 5 | `submit` | submission_pipeline.py | 端到端编排: suvtk tbl2asn / Sequin tbl2asn 双模式 → .sqn |
| 6 | `gui` | submission_gui/submission_gui.py | PySide6 桌面 GUI: 交互式编辑/验证/导出提交文件 |
| 7 | `report` | report_html.py | 交互式HTML全表编辑报告 |

---

## 目录结构 / Directory Structure

```
MMPV-RNA/
├── virome_discovery_pipeline/      # 新病毒发现 (18+4脚本)
│   ├── virome_pipeline.py          # 主编排器
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 辅助工具 (组装/鉴定/COBRA统计, Sankey)
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
├── metadata_gui/                   # 元数据管理桌面应用 (PyQt6)
│   ├── main.py                     # GUI 主入口
│   ├── controllers/                # 搜索桥接 / AI补全 / 元数据控制
│   ├── models/                     # 数据存储模型
│   ├── views/                      # 主窗口 / 搜索视图 / 表格 / 可视化 / 详情面板
│   └── utils/                      # 辅助工具
│
├── virome_submission_pipeline/     # 提交管线 (GenBank/CNCB)
│   ├── submission_pipeline.py      # 主编排器
│   ├── sequin_builder.py           # Sequin 构建器
│   └── ...                         # 元数据/报告/拓扑分析
│
├── submission_gui/                 # 提交桌面 GUI
│   └── submission_gui.py           # PySide6 交互式编辑/验证/导出
│
├── metadata_gui/                   # 元数据管理 GUI
│
├── biosoft/                        # 第三方工具 (脚本/JAR, 无需编译)
│   ├── VirBot/VirBot.py
│   ├── virhunter/predict_cpu.py + weights/
│   ├── ViReMa/ViReMa.py
│   └── snpEff/snpEff.jar + config + scripts/
│
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
- **元数据 GUI**: PyQt6 桌面应用, 支持 SRA/GSA 元数据可视化搜索、过滤、表格浏览和图表分析
- **提交管线**: GenBank/CNCB 序列提交自动化, 支持 Sequin 构建和假设蛋白分析

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
| `virome_submission_pipeline/` | 提交管线脚本头部均有详细 docstring |
| `doc/14-pipeline-diagrams.md` | Mermaid流程图/架构图 |
| `SOFTWARE_VERSIONS.txt` | 全部第三方软件版本记录 |
| `pixi.toml` | 依赖配置 (可直接查看) |

---

## 引用 / Citation

Zhang W. et al. **MMPV-RNA: Macro-metavirome Plant Virus Discovery and Analysis Platform**. 2026.

> If you use MMPV-RNA in your research, please cite the above reference.
