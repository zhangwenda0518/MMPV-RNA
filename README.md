# MMPV v3.1 — Massive Meta-mining of Plant Viruses

**大规模宏植物病毒挖掘与分析平台**

**中文** | [English](README_EN.md)

> 从公共测序数据中高通量鉴定植物病毒及其他病毒的完整闭环 — 含公共数据挖掘、病毒从头发现（de novo）、已知病毒定量/变异/进化深度分析。

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.2-blue.svg)](https://www.r-project.org/)
[![pixi](https://img.shields.io/badge/pixi-enabled-green.svg)](https://pixi.sh/)
[![BioConda](https://img.shields.io/badge/bioconda-supported-brightgreen.svg)](https://bioconda.github.io/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## 📖 平台 Wiki / Documentation

**20 页完整文档** (架构总览 · 七管线逐段解读 · 198 个核心脚本 `--help` 参数库 · 审计体系 · 部署运维):
👉 [平台 Wiki](https://github.com/zhangwenda0518/MMPV-RNA/wiki) — 从 [Home](https://github.com/zhangwenda0518/MMPV-RNA/wiki/Home) 进入

---

## 七大管线 / Seven Pipelines

```
public_metadata_pipeline/     → 数据获取 (测序数据检索/下载/转换 ＋ 宿主参考基因组获取/建库)
       │
       ├── 测序数据: SRA/GSA 搜索 → 元数据清洗 → 批量下载 → 可视化 → SRA→FASTQ.GZ
       └── 宿主参考: 多通道下载·或用户 FASTA → Kraken2/Bowtie2/HISAT2/Minimap2 四索引
              │
              ├──▶ endogenous_virus_pipeline/ → 内源性病毒（EVE）正向筛查
              │      滑窗 blastx vs 病毒参考蛋白 → 位点合并 →
              │      植物/病毒双侧 bitscore 判定 → RVDB 深度归属
              │      (输入是宿主基因组 FASTA, 与 reads 管线相互独立)
              │
              ▼
data_preprocessing_pipeline/  → 数据清洗 (质控 Fastp→Seqkit→Clumpify → 去宿主 Kraken2→比对→Ribodetector ＋可选 BBNorm 归一化)
              │
              ▼
virome_discovery_pipeline/    → 病毒发现（de novo）
       │
       ├── 组装 → 10工具并行鉴定 → 分层过滤 → COBRA延伸
       ├── 合并共组装 → CD-HIT参考聚类 → vclust去冗余 → 8工具分类 → 宿主预测
       └── CheckV评估 → 四支路级联拯救 (A/B/C/D) → HQ vOTU catalog
              │
              ▼
virome_analysis_pipeline/     → 已知病毒深度分析
       │
       ├── 快速定量 (Salmon/Bowtie2) → 变异检测 (FreeBayes/iVar)
       ├── SnpEff注释 → SNPGenie进化 → 12步全长组装
       ├── 全长相似性全景 (SDT) → DVG/重组检测 (ViReMa)
       └── 交互式HTML综合报告
              │
              ▼
virome_submission_pipeline/   → 数据提交
       │
       ├── 拓扑判断 → 元数据模板 → 假定蛋白注释
       └── suvtk tbl2asn / Sequin tbl2asn 双模式 → .sqn 提交文件
```

---

## 统一 I/O 目录布局 / Unified I/O Layout (v3.1)

五条管线默认沿用 v3.0 目录名（`legacy` 布局，checkpoint 兼容）。新项目可切
`standard` 布局——**每条管线一个独立输出根，根内各自独立编号**：

```
<项目>/ 01_PublicData/ 02_Preprocessing/ 03_Discovery/ 04_Analysis/ 05_Phylo/ 90_Handoff/
```

管线边界交接可自动定位（③→④ 桥一条命令；①→⑤ Core14 原生直读）：

```bash
export MMPV_IO_LAYOUT=standard          # 或各编排器 --io-layout standard
python virome_discovery_pipeline/utils/discovery2analysis.py \
    --from-discovery <项目>/03_Discovery --output_prefix <项目>/90_Handoff/analysis
```

完整项目树、五边界契约、legacy↔standard 映射表与迁移规则见
**`doc/IO_LAYOUT_DESIGN.md`**。

---

## 快速开始 / Quick Start

### 1. 一键部署

```bash
git clone https://github.com/zhangwenda0518/MMPV-RNA.git
cd MMPV-RNA

# 安装全部依赖 (pixi, 已逐包核查渠道存在性)
pixi install

# 注意: penguin / refineC / flye / PhaBOX2 / CHERRY / VirHunter / VirBot / ViraLM
#       为手动安装 (bioconda 无包), 见 pixi.toml 手动安装区与 doc/MMPV_dependencies.md
# 注意: virome_phylo_pipeline 的 BEAST/RDP5/PAML 等需按 doc/MMPV_dependencies.md A12 单独部署

# 安装 ViraLM 独立环境 (可选)
conda env create -f envs/viralm.yaml -n viralm

# 下载参考数据库 (见下方数据库清单)
```

> 宿主基因组获取支持多通道：`--download-source auto|datasets|ngd|ftp|gget`（auto = datasets 失败自动回退 ncbi-genome-download / FTP 直连；gget 走 Ensembl/Ensembl Plants，适合 NCBI 缺物种的非模式植物），自有基因组直接 `--genome-fasta` 跳过下载。

### 1b. 数据预处理（宿主参考建库 ＋ 质控 → 去宿主，两条管线）

```bash
# 数据质控: Fastp → Seqkit → Clumpify
python data_preprocessing_pipeline/clean-data.py -i raw_fastqs/ -o out/00a_CleanData/ -j 20 -t 16

# 宿主参考: NCBI 多通道下载 → 四索引 (或 --genome-fasta 用自有基因组)
python public_metadata_pipeline/build_host_pipeline.py --species "Solanum lycopersicum" --taxid 4081 --stage all

# 去宿主 + 去rRNA: Kraken2 → 比对 → Ribodetector
python data_preprocessing_pipeline/host_depletion.py     -k kraken2_db/ -x bowtie2_idx/ -I out/00a_CleanData/ -O out/00b_HostDepletion/     --tool bowtie2 --seq-type rna-short --rrna -t 40 -j 10

# 报告 + 交接清单: preprocessing_summary.tsv + assembly_ready.list + HTML
# (一条龙 --stage all = clean → deplete → report 自动收尾)
python data_preprocessing_pipeline/data_preprocessing.py --stage report --output_dir out/

# 可选: 覆盖度归一化 (仅共组装前、深度差异大时才需要; 默认不跑 → out/00c_BBnorm/)
python data_preprocessing_pipeline/data_preprocessing.py --stage bbnorm --output_dir out/ -t 16 -j 4

# 更省的做法: 公共数据管道 --stage all 已包含 hostref(宿主基因组下载)+hostdb(四索引建库) 与 report,
# 检索物种即宿主物种, 一次跑通 检索→下载→转换→建库→报告:
# python public_metadata_pipeline/public_data_pipeline.py --species "..." --taxid ... --stage all

# 公共数据从 .sra 一条命令到 clean reads (可选 --bbnorm 加做归一化):
# python public_metadata_pipeline/preprocess_unified.py --sra-dir ./sra/ --outdir out/ --bbnorm \
#     --kraken2-db kraken2_db/ --step2-index bowtie2_idx/host -t 40 -j 4
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

### 3. 内源性病毒（EVE）筛查

```bash
# 单基因组: 宿主参考基因组 → EVE 位点 (数据库路径从 pipeline_config.yaml 读)
python endogenous_virus_pipeline/eve_screen.py \
    -g host_reference/genome/all.genome.uniq.fasta \
    -n Solanum_lycopersicum -o out/eve/ -t 40

# 批量: batch.tsv 每行 NAME<TAB>基因组路径 (.tar.gz 包也可直接填, Stage1 自动解包;
#       NAME 必须唯一, 重名开跑前报错 — 重名会共用 01_Loci/<NAME>/ 产物)
python endogenous_virus_pipeline/eve_screen.py \
    -B batch.tsv -o out/eve/ -t 40 -J 5 --fast
# -J 并行基因组数, -t 每基因组线程数, 总核数 = J x t (256 核推荐 -J5 -t40)

# 断点续传: 按请求的阶段判定 (阶段产物齐全才算完成), 重发同命令即可续跑
# 汇总: kingdom_summary.tsv / family_by_genome.tsv (--merge 单独重跑汇总)
# 抽序列: samtools faidx (需 samtools >= 1.11, 版本不够开跑前体检即报错)
#          索引写到输出目录, 不依赖参考基因组目录可写
# 防挂死: --cmd-timeout 7200 (单条外部命令超时秒数, 0=不限)
```

判明"是宿主基因还是内源病毒化石"，以及判别病毒组组装出的候选 contig 是真 DNA 病毒
还是退化 EVE，另有两个模块：`eve_screen.py` 做宿主基因组侧的位点发现，
`eve_distinguish/` 做候选侧的判别（跑完发现管道后单独调用）：

```bash
# DNA 病毒候选 vs EVE 判别 (详见 endogenous_virus_pipeline/eve_distinguish/README.md)
bash endogenous_virus_pipeline/eve_distinguish/run_all.sh \
    -D /path/to/discovery_out -A /path/to/1kp/assemblies -t 32
#   产物: dna_vs_eve_filter.tsv (每条候选一个 action, 下游按列过滤即可)
```

### 4. 分析管线

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

### 5. 公共数据管线

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

### 6. 提交管线

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
| 2 | `info` | gsa_sra.info.py | SRA XML解析 + GSA爬虫 + AI元数据清洗 → Global_Unified_Metadata_Core14.csv (14列) |
| 3 | `down` | gsa_sra.down.py | aria2c/wget/prefetch 双协议下载 → 原始 FASTQ/SRA |
| 4 | `convert` | sra2fastx.py | fasterq-dump 批量 SRA → FASTQ.GZ（失败回退 fastq-dump） |
| 5 | `plot` | gsa_sra.plot.py | 时间/组织/地区/机构 SCI级6面板可视化 |
| 6 | `hostref` | download_host_genome.py | 宿主基因组下载 (datasets/ngd/FTP/gget 四通道; `--host-fasta` 用户基因组跳过下载) → `host_reference/genome/` |
| 7 | `hostdb` | build_hostbase.py | Kraken2/Bowtie2/HISAT2/Minimap2 四索引构建 → `host_reference/hostdb/`（默认复用 --species/--taxid） |
| 8 | `report` | generate_report.py | 交互式 HTML 报告 + `sample_handoff.csv` 交接清单 (Run→FASTQ 路径+元数据, 交棒数据清洗管道) |

### 宿主参考数据库构建

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `genome-down` | public_metadata_pipeline/download_host_genome.py | 多通道下载参考基因组 (datasets/ngd/FTP/gget · 或 --genome-fasta 用户提供) + GFF3 + 细胞器基因组 |
| 2 | `hostdb` | public_metadata_pipeline/build_hostbase.py | Kraken2 + Bowtie2 + HISAT2 + Minimap2 四种索引构建 |

### 统一预处理入口 (从 .sra 一条命令到 clean reads)

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `convert` | sra2fastx.py | .sra → FASTQ.GZ（可选，已有 FASTQ 时 `--skip-convert`） |
| 2 | `clean` | data_preprocessing_pipeline/clean-data.py | Fastp → Seqkit → Clumpify → `00a_CleanData/` |
| 3 | `hostref` | download_host_genome.py + build_hostbase.py | 可选：宿主基因组四通道下载 + 四索引构建 |
| 4 | `deplete` | data_preprocessing_pipeline/host_depletion.py | Kraken2 → 比对 → Ribodetector → `00b_HostDepletion/` |
| 5 | `bbnorm` | data_preprocessing_pipeline/run_bbnorm.py | **可选**（`--bbnorm`，默认不跑）：覆盖度归一化 → `00c_BBnorm/` |

入口脚本: `python public_metadata_pipeline/preprocess_unified.py --sra-dir ... --outdir ...`
（产出目录名与 `data_preprocessing_pipeline` / 发现管线一致，可直接作为 `--input_reads`）

---

## 发现管线阶段 / Discovery Pipeline Stages

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `clean` | data_preprocessing_pipeline/clean-data.py | Fastp QC + Seqkit FASTA转换 + Clumpify 聚类重排（提升压缩率） |
| 2 | `deplete` | data_preprocessing_pipeline/host_depletion.py | Kraken2 → Bowtie2/HISAT2/Minimap2 → rRNA去除 |
| 3 | `assembly` | assembly_pipeline.py | MEGAHIT / rnaviralSPAdes / Penguin 组装 |
| 4 | `identification` | virus_identification.py | 10工具并行鉴定 (Genomad+Diamond BLASTX+RdrpCatch+ViraLM+VirBot+VirSorter2+ViralVerify+VirHunter+Metabuli+Viroid BLASTN) |
| 5 | `filter` | filter_virus.py | UniProt-strict 高置信过滤 |
| 6 | `cobra` | cobra_pipeline.py | BWA-MEM2 → CoverM → COBRA-Meta 单样本延伸 |
| 7 | `merge` | (内置, Flye) | 多样本 Flye 共组装 (co-assembly 模式) |
| 8 | `cluster` | cluster_pipeline.py | CD-HIT参考引导预聚类 + vclust Leiden聚类 |
| 9 | `taxonomy` | virus_classifier.py + R | 8工具分类 (Genomad+Metabuli+CAT+Diamond LCA+VITAP+MMseqs2+ACVirus+vConTACT3) 加权投票 → 8级 taxonomy |
| 10 | `host` | run_host_prediction.py | ICTV > RNAVirHost > PhaBOX2 三级宿主预测 + 植物专属后过滤 |
| 11 | `checkv` | (内置) | CheckV 完整性评估 |
| 12 | `rescue` | rescue_pipeline.py | 四支路级联拯救 (CheckV → Virseqimprover → BLASTN/PlantVirusDB 参考重建 → genus-length 兜底) |
| 13 | `analysis` | (内置) | 整合拯救证据 + 病毒组下游分析 |
| 14 | `analysis_verify` | (内置) | 09b 分析复核 (HMM/CT3 证据扫描) |
| 15 | `report` | report_pipeline.py | TSV汇总 + Sankey图 + 交互式HTML报告 |

> 注: `clean` / `deplete` 两阶段的脚本独立为顶层管线 `data_preprocessing_pipeline/`
> （含 clean→deplete 一条龙入口）；宿主参考基因组获取/建库在 `public_metadata_pipeline/`
> （`virome_pipeline.py` 仍可整链编排调用）。
>
> 注: **BBNorm 覆盖度归一化已不属于发现管线阶段**（供共组装前可选用途），
> 归 `data_preprocessing_pipeline/run_bbnorm.py`，两个入口均默认不跑：
> `data_preprocessing_pipeline/data_preprocessing.py --stage bbnorm`（`--stage all` 不含）
> 与公共数据入口 `public_metadata_pipeline/preprocess_unified.py --bbnorm`，
> 输出统一为 `00c_BBnorm/`，可直接作为 `virome_pipeline.py --input_reads`。

---

## 分析管线阶段 / Analysis Pipeline Stages

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `detect` | batch_virus_depth.py | Salmon/Kallisto/Bowtie2 快速定量 + Poisson过滤 |
| 2 | `filter` | utils/filter_summary.py | 高置信度过滤 |
| 3 | `variants` | batch_virus_variants.py | FreeBayes/iVar/LoFreq 变异检出 + 共识序列构建 |
| 4 | `post` | virus_vcf_pipeline.py | VCF可视化 + SnpEff宏 + MAF 瀑布 + SnpGenie dN/dS |
| 5 | `full` | batch_virus_full.py → virus-full.py | **全长基因组构建 (OmniVirusAssembler V9.0)**: 12步精炼 — SHIVER-like Divine Fusion 实心骨架 → Reads级迭代抛光 (单碱基精度) → gmcloser + abyss-sealer 双引擎补洞 → 环化检测 (Circular=True) |
| 6 | `extract` | utils/extract_full_fasta.py | 最长contig提取 + 参考N-fill (RNA: N<5% 才填) |
| 7 | `similarity` | virus_auto_pipeline.py | 全长相似性热图 + 层次聚类 (MAFFT --auto) |
| 8 | `dvg` | batch_virema_dvg.py | ViReMa DVG/重组检测 + Circos图 |
| 9 | `report` | generate_pipeline_report.py | 交互式HTML综合报告 + AI 解读提示词 |

> 注: 正选择分析 (`capheine` / HyPhy FEL·MEME·BUSTED·PRIME) 已于 v3.0 迁移至
> `virome_phylo_pipeline/` (STAGE_ORDER 中的 `capheine` 阶段)，本管线 9 阶段定义见
> `auto_known_virus.py` 文件头 docstring。

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
├── data_preprocessing_pipeline/    # 数据清洗 (质控 clean-data.py → 去宿主 host_depletion.py; 含一条龙入口)
├── virome_discovery_pipeline/      # 病毒发现 (de novo, 15 阶段)
│   ├── virome_pipeline.py          # 主编排器
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 辅助工具 (组装/鉴定/COBRA统计, Sankey)
│
├── virome_analysis_pipeline/       # 已知病毒深度分析 (9 阶段)
│   ├── auto_known_virus.py         # 分析编排器
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 辅助工具
│
├── public_metadata_pipeline/       # 公共数据获取 (8 阶段)
│   ├── public_data_pipeline.py     # 公共数据编排器 (检索/元数据/下载/可视化)
│   ├── build_host_pipeline.py      # 宿主参考编排 (多通道下载/用户FASTA → 四索引)
│   ├── preprocess_unified.py       # 全链路预处理入口 (convert→clean→hostref→deplete→[bbnorm])
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 共享工具模块
│
├── metadata_gui/                   # 元数据管理桌面应用 (PySide6)
│   ├── main.py                     # GUI 主入口
│   ├── controllers/                # 搜索桥接 / AI补全 / 元数据控制
│   ├── models/                     # 数据存储模型
│   ├── views/                      # 主窗口 / 搜索视图 / 表格 / 可视化 / 详情面板
│   └── utils/                      # 辅助工具
│
├── virome_phylo_pipeline/          # 系统发育/进化/群体遗传 (18 stage, 8 大模块)
│   ├── phylo_pipeline.py           # 总管线 (BEAST 定年 + pypopart + 正选择)
│   └── doc/ 权威体系 (RUN_GUIDE/STAGE_REFERENCE/METHODS)
├── virome_submission_pipeline/     # 提交管线 (GenBank/CNCB)
│   ├── submission_pipeline.py      # 主编排器
│   ├── sequin_builder.py           # Sequin 构建器
│   └── ...                         # 元数据/报告/拓扑分析
│
├── endogenous_virus_pipeline/      # 内源性病毒 (EVE): 宿主基因组筛查 + 候选判别
│   ├── eve_screen.py               # 筛查主编排器 (三阶段: discover→verdict→annotate→merge)
│   ├── eve_genome_scan.py          # 单基因组 worker (按阶段产物断点续传)
│   ├── eve_scan_core.py            # 纯逻辑核心 (坐标还原/位点合并/判定/汇总)
│   ├── tests/                      # 142 个单测 (核心逻辑/判别逻辑/产物不变量)
│   ├── doc.md                      # 完整流程文档
│   └── eve_distinguish/            # ★ 后运行脚本: DNA 病毒候选 vs EVE 判别 (v6.2)
│       ├── run_all.sh              # 一键入口 (0 面板 → 1 建库 → 2 blastx → 3 结构/域 → 4 寄主 → 5/6 判定)
│       ├── s1_decay_scan.py        # 结构退化 (终止子富集 + 分布剖面)
│       ├── s2_domain_scan.py       # 域架构 (MP/CP/AP/RT/RH 组件记功)
│       ├── s2b_locus_scan.py       # 候选 blastn 回寄主 → 基因座架构
│       ├── s3_verdict.py           # 退化 × 架构 二维判定 + 寄主否决
│       ├── s4_filter.py            # verdict → action (REMOVE/MOVE_EVE/KEEP_virus/REVIEW)
│       ├── build_panel.py          # 参考面板构建 (随模块提交 panel.fasta/baits.fa)
│       └── README.md               # 方法与踩坑记录 (含逐条实测数字)
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
├── doc/                            # 文档 (17篇 + 方法章总模板 + 依赖清单)
│   ├── METHODS_TEMPLATE.md         # 五大管道论文级方法章合并版 (SCI投稿直接取用)
│   └── MMPV_dependencies.md        # 外部依赖清单 (工具/数据库/路径)
├── scripts/                        # 开发期一次性脚本归档 (非核心管线, 见 scripts/README.md)
│   ├── blacklist/                  # 宿主黑名单构建/同步/验证 (bl_layer1~4)
│   ├── db_query/                   # 结果库临时查询 (q_*)
│   ├── patch/                      # 已合入管线的一次性补丁留档
│   ├── audit/                      # 结果审计/交叉核对
│   ├── pilot_rvdb/                 # RVDB 试点分析
│   ├── mapper_bench/               # 比对器基准测试
│   ├── utils/                      # md2docx / 服务器端启动器
│   └── check_script_refs.py        # 管线脚本引用完整性检查
├── paper/                          # SCI 论文稿 (英文全流程初稿)
├── reports/                        # 生成的 HTML 报告输出 (gitignored)
├── envs/viralm.yaml                # ViraLM 独立conda环境
├── pixi.toml                       # 一键部署 (conda 可装依赖; 手动工具见文件尾注释)
├── pipeline_config.yaml            # 管线配置文件
├── SOFTWARE_VERSIONS.txt           # 全部第三方软件版本
└── README.md
```

---

## 核心特性 / Key Features

- **10工具并行病毒鉴定**: Genomad + Diamond BLASTX + VirSorter2 + ViralVerify + VirHunter + Metabuli + RdrpCatch + ViraLM + VirBot + Viroid BLASTN
- **CD-HIT参考引导预聚类**: 整合ICTV/NCBI完整基因组为参考，关联碎片contig到已知物种
- **四支路级联拯救**: A: CheckV≥90% 直通 → B: Virseqimprover reads 延伸 → C: BLASTN 参考引导重建 → D: genus-length 兜底，逐步提升 HQ vOTU 产出
- **宿主过滤**: 拯救前按宿主类别预过滤，节省70%+计算量
- **pixi一键部署**: 100个conda包精确版本管理，`pixi install`即可
- **断点续传**: 全部脚本支持 `--resume/--force`，大规模运行安全中止和恢复
- **基于分类层级的新颖性判断**: 无需BLASTN依赖，直接从taxonomy completeness判断新病毒
- **完整闭环**: 从公共数据挖掘到投稿图表，一站式完成
- **元数据 GUI**: PySide6 桌面应用, 支持 SRA/GSA 元数据可视化搜索、过滤、表格浏览和图表分析
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
| 植物/病毒拆分库 (hybrid) | EVE 双侧判定 | SwissProt 拆分 + id2div 映射 (endogenous_virus_pipeline) |

**完整清单与部署方案**：全管线共依赖 50 项数据库（必需 25 / 可选 25，全量约 710GB）。
逐库的下载与建索引命令见 `DATABASE_SETUP.md`；**可复现的集中部署方案**（共享库根 + 软链接
兼容层 + 环境变量 + 校验）见 `doc/DB_DEPLOYMENT.md`，配套工具：

```bash
python scripts/db_setup/mmpv_db.py list                                 # 打印 50 项依赖清单
python scripts/db_setup/mmpv_db.py check                                # 部署前体检
python scripts/db_setup/mmpv_db.py verify                               # 校验现有库是否齐全
python scripts/db_setup/mmpv_db.py down --apply                          # 下载/构建可自动获取的库
python scripts/db_setup/mmpv_db.py adopt --db-from ~/database --apply   # 迁移到共享库根
python scripts/db_setup/mmpv_db.py link --apply                          # 建旧路径兼容软链接
```

换服务器时用 `down` 拉取；其中 RVDB 派生的一整套库（diamond / mmseqs / metabuli / CAT）
由仓库自带的 `virome_discovery_pipeline/build_virus_db.py` 从上游 RVDB fasta 现场构建，
详见部署方案 §3.2。

---

## 文档 / Documentation

| 文档 | 说明 |
|------|------|
| `virome_discovery_pipeline/doc.md` | 发现管线全部脚本 — 参数/输入输出/结果解读 |
| `virome_analysis_pipeline/doc.md` | 分析管线全部脚本 — 参数/输入输出/结果解读 |
| `public_metadata_pipeline/doc.md` | 公共数据管线 + 宿主参考构建（含四通道下载） |
| `virome_submission_pipeline/` | 提交管线脚本头部均有详细 docstring |
| `doc/14-pipeline-diagrams.md` | Mermaid流程图/架构图 |
| `doc/DB_DEPLOYMENT.md` | 数据库部署方案（共享库根/软链接/环境变量/校验，可复现） |
| `SOFTWARE_VERSIONS.txt` | 全部第三方软件版本记录 |
| `pixi.toml` | 依赖配置 (可直接查看) |

---

## 引用 / Citation

Zhang W. et al. **MMPV: Massive Meta-mining of Plant Viruses**. 2026.

> If you use MMPV in your research, please cite the above reference.
