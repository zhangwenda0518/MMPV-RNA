# public_metadata_pipeline — 公共数据分析流程文档

## 数据流总览

```
public_data_pipeline.py (主编排器)
  │
[Stage 1] search  → gsa_sra.search.py    → SRA_GSA_Merged_Final.csv
[Stage 2] info    → gsa_sra.info.py      → Global_Unified_Metadata_Core14.csv
[Stage 3] down    → gsa_sra.down.py      → downloaded SRA
[Stage 4] convert → sra2fastx.py         → FASTQ.GZ
[Stage 5] plot    → gsa_sra.plot.py      → SCI 级可视化图表
[Stage 6] hostref → download_host_genome.py → host_reference/genome/all.genome.uniq.fasta (四通道/用户FASTA)
[Stage 7] hostdb   → build_hostbase.py     → host_reference/hostdb/ 四索引
[Stage 8] report   → generate_report.py    → Pipeline_Summary_Report.html + sample_handoff.csv

宿主参考构建 (与测序数据获取同管道, build_host_pipeline.py 编排)
  │
[Stage 1] genome-down  → download_host_genome.py  → all.genome.uniq.fasta (四通道: datasets/ngd/FTP/gget)
[Stage 2] hostdb       → build_hostbase.py        → kraken2/bowtie2/hisat2/minimap2 索引

统一预处理入口 (preprocess_unified.py, 按需单独调用; 从 .sra 起步的公共数据 on-ramp)
  │
convert → clean → [hostref] → deplete → [bbnorm]
  → 00_Converted/ 00a_CleanData/ 00b_HostDepletion/ [00c_BBnorm/]
```

---

## 1. public_data_pipeline.py — 公共数据主编排器

**用途:** 八阶段公共数据获取调度器 (`search → info → down → convert → plot → hostref → hostdb → report`)。按物种拉丁名和 NCBI TaxID，搜索 NCBI SRA + CNCB GSA 数据库，清洗元数据，下载原始测序数据，转换 FASTQ，生成出版级可视化，下载宿主参考基因组并构建四索引，收尾产出 HTML 报告与下游交接清单。检查点机制避免重复执行。各 stage 直接调用工作脚本，不嵌套子管道。

### 8 阶段说明

| 阶段 | 脚本 | 输入 | 输出 |
|------|------|------|------|
| search | gsa_sra.search.py | — | `search/SRA_GSA_Merged_Final.csv` |
| info | gsa_sra.info.py | `sra.list` | `info/Global_Unified_Metadata_Core14.csv` |
| down | gsa_sra.down.py | `sra.list` | `down/{sample}/*.sra` |
| convert | sra2fastx.py | `down/*.sra` | `down/{sample}/*.fastq.gz` |
| plot | gsa_sra.plot.py | 上述 CSV | `plot/Combined_Landscape_Full.pdf` |
| hostref | download_host_genome.py | --host-species/--host-taxid (默认复用 --species/--taxid)；`--host-fasta` 跳过下载 | `host_reference/genome/all.genome.uniq.fasta` |
| hostdb | build_hostbase.py | 上一步基因组 FASTA | `host_reference/hostdb/{kraken2,bowtie2,hisat2,minimap2}/` |
| report | generate_report.py | 上述全部产物 | `Pipeline_Summary_Report.html` + `sample_handoff.csv` |

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--species` | 必需 | 物种拉丁学名 |
| `--taxid` | 必需 | NCBI Taxonomy ID |
| `--stage` | 必需 | search / info / down / convert / plot / all |
| `--deepseek-api` | 环境变量 | DeepSeek API Key (AI 元数据清洗) |
| `--ncbi-api` | 环境变量 | NCBI API Key |
| `--work-dir` | `./public_data_pipeline_output` | 输出根目录 |
| `--source-type` | `TRANSCRIPTOMIC` | 测序类型过滤 |
| `--detailed` | True | 详细模式 (提取 Tissue/Age/Location) |
| `--ngdc-method` | `aria2c` | 下载工具 |
| `--threads` | 4 | 线程数 |

### 输入
无 (所有数据从 NCBI E-utilities + CNCB GSA 在线获取)

### 输出
```
{work_dir}/
  search/SRA_GSA_Merged_Final.csv              ← 搜索结果
  info/sra.list                                ← Run accession 列表
  info/Global_Unified_Metadata_Core14.csv      ← 14 列统一元数据
  down/{sample_name}/{accession}.sra           ← 下载的测序数据 (原始 SRA)
  down/{sample_name}/{accession}.fastq.gz      ← 转换后的测序数据 (convert 阶段)
  plot/from_search/Combined_Landscape_Full.*   ← 搜索数据可视化
  plot/from_info/Combined_Landscape_Full.*     ← 元数据可视化
```

### 上下游衔接
- **→ virome_discovery_pipeline:** 下载的 FASTQ 作为 `--input_reads` 输入
- **→ virome_analysis_pipeline:** 元数据 CSV 供样本注释

---

## 2. gsa_sra.search.py — 双引擎物种检索

**用途:** 同时搜索 NCBI SRA (E-utilities) 和 CNCB GSA (网页爬虫)，合并去重输出统一 CSV。detailed 模式下解析 XML/Excel 提取 Tissue, Age, Location 等深层特征，可选 DeepSeek AI 清洗。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-q` / `--query` | 必需 | 物种拉丁名 |
| `-s` / `--source` | — | 测序类型 (如 TRANSCRIPTOMIC) |
| `-o` / `--outdir` | `./Global_Species_Results` | 输出目录 |
| `--detailed` | False | 详细模式 |
| `--ncbi-api` | — | NCBI API Key |
| `--deepseek-api` | — | DeepSeek API Key |

### 输入
无 (在线 API 调用)

### 输出
```
{outdir}/
  SRA_GSA_Merged_Final.csv          ← 合并结果 (Database, Run, BioProject, BioSample, Tissue, Age, Location, ...)
  GSA_Results/0_web_cache/          ← GSA 网页缓存
  GSA_Results/1_xls_cache/          ← GSA Excel 缓存 (detailed 模式)
```

### 结果解读
- `Database` 列: `SRA` (NCBI) 或 `GSA` (CNCB)
- `Run` 列: SRR/ERR/DRR/CRR accession
- `ReleaseDate` 列: 数据发布日期
- `Tissue`, `Age_GrowthStage`, `Location` 列: detailed 模式下提取的深层特征

---

## 3. gsa_sra.info.py — 全局元数据统一引擎

**用途:** 从 SRA/GSA accession 列表出发，下载并解析 SRA XML 元数据 (esearch|efetch)，爬取 GSA 网页和 Excel，可选 AI (DeepSeek/Kimi) 推理/仲裁，BioProject→PubMed 文献追溯，分类学解析，产生 14 列标准化元数据 (CSV+TSV) 和交互式 datavzrd HTML 报告。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input` | 必需 | accession 列表文件 |
| `-o` / `--outdir` | `./Global_Metadata_Results` | 输出目录 |
| `-m` / `--mode` | `both` | local / api / both |
| `-t` / `--threads` | 4 | 线程数 |
| `--deepseek-api` | — | DeepSeek API Key |
| `--ncbi-api` | — | NCBI API Key |
| `--use-scholar` | False | 允许 Google Scholar |
| `--fill-date` | False | 时间兜底填充 |

### 输入
文本文件，每行一个 SRA/GSA accession

### 输出
```
{outdir}/
  Global_Unified_Metadata_Core14.csv    ← 14 列核心元数据
  Global_Unified_Metadata_Core14.tsv
  Global_Unified_Metadata_Full.csv      ← 全部字段
  SRA_Results/                          ← SRA 处理中间文件
    1_raw_xml/{srr}.xml
    2_full_json/{srr}_full.json
    3_3_ai_arbitrated/{srr}_arb.json   ← AI 仲裁结果
    SRA_Ultimate_Merged.csv
  GSA_Results/                          ← GSA 处理中间文件
    0_web_cache/{acc}.json
    2_xlsx/{CRA}.xlsx
    GSA_Ultimate_Merged.csv
  BioProject_Results/                   ← 文献追溯
    BioProject_Trace_Details.csv
  Report_Global_Unified_Metadata_Full/  ← datavzrd 交互报告
```

### 核心 14 列
`Run`, `ReleaseDate`, `CollectionDate`, `Location`, `Source`, `Tissue`, `Age_GrowthStage`, `ScientificName`, `TaxID`, `LibrarySource`, `CenterName`, `BioProject`, `BioSample`, `PMID`

**工具依赖:** esearch, efetch (NCBI E-utilities CLI), datavzrd

---

## 4. gsa_sra.down.py — SRA 智能下载器

**用途:** 从 NGDC (aria2c/wget/requests 双协议 FTP/HTTP 回退) 和 NCBI (prefetch) 下载原始测序数据。支持跳过列表、TSV 矩阵组织、进度追踪、失败列表。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--list` | 三选一 | accession 列表文件 (TXT) |
| `--srr` | 三选一 | 单个 accession |
| `--tsv` | 三选一 | TSV 矩阵 (行=样本, 列=accession) |
| `--skip-list` | — | 跳过列表 |
| `--ngdc-method` | `aria2c` | requests / aria2c / wget |
| `--ngdc-concurrency` | 1 | NGDC 并发数 |
| `--prefetch-concurrency` | 4 | Prefetch 并发数 |
| `-o` / `--output` | `./sra_data` | 输出目录 |

### 输入
accession 列表文件 (每行一个 Run ID)

### 输出
```
{output}/
  {sample_name}/{accession}.sra      ← 下载的原始数据 (SRA)
  download_report_{timestamp}.csv    ← 下载状态报告
  failed_sra_{timestamp}.txt         ← 失败列表 (供重试)
```

**工具依赖:** aria2c, wget, prefetch (SRA Toolkit)

### 下载策略
- **CRR** (CNCB): NGDC FTP → HTTP 回退 → aria2c/wget/requests
- **SRR/ERR/DRR** (NCBI): prefetch（SRA 下载；fasterq-dump 转换见 convert 阶段）

---

## 5. gsa_sra.plot.py — SCI 级宏观可视化

**用途:** 读取合并元数据 CSV，生成 3×2 面板出版级图表: (A) 时间分布折线图, (B) 数据库占比环形图, (C) Top 10 机构, (D) Top 研究组织, (E) Top 采样地区, (F) Top 发育阶段/年龄。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-i` / `--input` | 必需 | CSV 文件路径 |
| `-o` / `--outdir` | `SCI_Figures_Output` | 输出目录 |

### 输入
含 `Run`, `ReleaseDate`, `Organization_CenterName`, `Tissue`, `Location`, `Age_GrowthStage` 列的 CSV

### 输出
```
{outdir}/
  Combined_Landscape_Full.pdf   ← 矢量图
  Combined_Landscape_Full.png   ← 光栅图 (150 dpi)
```

---

## 6. build_host_pipeline.py — 宿主参考编排器

**用途:** 两阶段编排: `genome-down` (多通道下载参考基因组: datasets → ngd → FTP → gget 自动回退；或 `--genome-fasta` 直接用自有 FASTA) → `hostdb` (build_hostbase.py 构建 Kraken2/Bowtie2/HISAT2/Minimap2 四索引)。全检查点化。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--species` | 必需 | 物种拉丁学名 |
| `--taxid` | 必需 | NCBI Taxonomy ID |
| `--stage` | 必需 | genome-down / hostdb / all |
| `--work-dir` | `./build_host_pipeline_output` | 输出目录 |
| `--genome-fasta` | — | 已有基因组 FASTA (跳过下载) |
| `--download-source` | `auto` | 下载通道: auto / datasets / ngd / ftp / gget |
| `--hostdb-tools` | `kraken2,bowtie2,hisat2,minimap2` | 建库工具 |
| `--seq-types` | `dna-short,rna-short,nanopore,pacbio` | Minimap2 序列类型 |
| `--k2-libs` | `archaea,bacteria,plasmid,fungi,protozoa,UniVec` | Kraken2 标准库 (编排调用会覆盖 build_hostbase 自身默认 7 库) |
| `--taxonomy-dir` | — | Kraken2 本地 Taxonomy 目录 |
| `--ncbi-api` | — | NCBI API Key |
| `--threads` | 30 | 线程数 |
| `--force` / `--dry-run` | False | 强制重跑 / 仅预览 |

### 输入
无 (或 `--genome-fasta` 跳过下载)

### 输出
```
{work_dir}/
  genome/all.genome.uniq.fasta       ← 合并去重基因组
  hostdb/
    kraken2/                         ← Kraken2 数据库 (含 taxonomy + standard libraries)
    bowtie2/host.*.bt2               ← Bowtie2 索引
    hisat2/host.*.ht2                ← HISAT2 索引
    minimap2/host_{type}.mmi         ← Minimap2 索引
```

### 上下游衔接
**→ data_preprocessing_pipeline:** `hostdb/` 直接作为 `host_depletion.py -k/-x` 输入与发现管线 `--host_db`；
**→ preprocess_unified.py:** 其 hostref 阶段自动调用本脚本。

**工具依赖:** datasets (NCBI CLI), ncbi-genome-download (pip), gget (pip), requests, unzip, build_hostbase.py

---

## 7. build_hostbase.py — 竞争性宿主数据库构建器

**用途:** 单脚本构建四种宿主索引: Kraken2 (含 taxonomy 下载+标准库+TaxID 注入), Bowtie2, HISAT2, Minimap2。支持检查点续传、软链接挂载 taxonomy、补充自定义 FASTA。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--tool` | 必需 | kraken2 / bowtie2 / hisat2 / minimap2 (逗号分隔) |
| `-i` / `--input` | 必需 | 宿主 FASTA |
| `-o` / `--output` | 必需 | 输出目录 |
| `-t` / `--threads` | `os.cpu_count()` | 线程数 |
| `--taxonomy` | — | Kraken2 本地 Taxonomy 目录 |
| `--taxid` | — | 非模式物种 TaxID |
| `--k2-libs` | `archaea,bacteria,plasmid,fungi,protozoa,UniVec,UniVec_Core` | Kraken2 标准库 (独立运行默认 7 库; 经编排器传入时被覆盖为 6 库) |
| `--add-library` | — | 补充自定义 FASTA |
| `--force` | False | 强制重建 |

### 输入
宿主参考基因组 FASTA

### 输出
```
{output}/
  kraken2/*.k2d                       ← Kraken2 编译数据库
  bowtie2/host.{1..4}.bt2, rev.*.bt2  ← Bowtie2 索引
  hisat2/host.{1..8}.ht2              ← HISAT2 索引
  minimap2/host_{type}.mmi            ← Minimap2 索引
```

**工具依赖:** kraken2-build, k2, bowtie2-build, hisat2-build, minimap2

---

## 8. download_host_genome.py — 参考基因组多通道下载器

**用途:** 宿主参考基因组获取，四通道自动回退: `datasets` (NCBI datasets CLI, 官方首选) → `ngd` (ncbi-genome-download, 按 TaxID 并行) → `ftp` (E-utilities 定位 + HTTPS 直连) → `gget` (Ensembl/Ensembl Plants, NCBI 缺物种的非模式植物补充)。下载核基因组 + GFF3 + seq-report，自动解压、合并多 FASTA、去重序列名，可选下载叶绿体/线粒体基因组。

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--species` | 必需 | 物种拉丁学名 |
| `--outdir` | 必需 | 输出目录 |
| `--source` | `auto` | 下载通道: auto / datasets / ngd / ftp / gget |
| `--taxid` | — | NCBI Taxonomy ID (ngd 通道必需) |
| `--ensembl-release` | — | Ensembl 版本 (gget 通道) |
| `--ncbi-api` | — | NCBI API Key (提升速率) |
| `--include-organelles` | False | 下载叶绿体和线粒体 |
| `--min-length` | 0 | 最小序列长度 |
| `--verify-only` | False | 仅验证现有下载 |
| `--skip-datasets` | False | 跳过 NCBI datasets 步骤 |

### 输入
无

### 输出
```
{outdir}/
  genome_down.zip                ← 原始下载
  extracted/ncbi_dataset/        ← 解压内容
  all.genome.uniq.fasta          ← 合并去重基因组
  genome_report.txt              ← 摘要报告
  chloroplast/chloroplast.fasta  ← (可选) 叶绿体
  mitochondrion/mitochondrion.fasta ← (可选) 线粒体
```

**工具依赖:** datasets (NCBI CLI), ncbi-genome-download (pip), gget (pip), requests, unzip

---

## 9. 共享工具模块

### utils/pipeline_utils.py

三个共享类:

- **UI** — 终端彩色输出 (banner, stage, ok, warn, err, info)
- **Checkpoint** — JSON 检查点管理 (`{work_dir}/.checkpoints/state.json`)
- **run_cmd(cmd, stage_name, log_dir)** — 安全 shell 命令执行，含日志记录、密钥脱敏、超时控制

### utils/extract_sra_list.py

从 CSV 提取 SRA accession 列表，支持按列过滤:
```
python extract_sra_list.py --input merged.csv --output sra.list --filter-db SRA --filter-col Tissue --filter-val leaf
```

---

## 10. preprocess_unified.py — 统一预处理入口（convert → clean → [hostref] → deplete → [bbnorm]）

**用途:** 从 `.sra` 一条命令跑到"宿主剔除后 clean reads"的公共数据 on-ramp。
按阶段调用上游脚本（不嵌套子管道），产出目录名与 `data_preprocessing_pipeline/`
及发现管线一致，因此输出可直接作为下游 `--input_reads` / `--reads_dir`。
之所以放在本目录，是因为其中的 `convert`（`sra2fastx.py`）与 `hostref`
（`download_host_genome.py` + `build_hostbase.py`）三个脚本只在本目录。

### 5 阶段说明

| # | 阶段 | 输入 → 输出 | 说明 |
|---|------|-------------|------|
| 1 | `convert` | `.sra` → `00_Converted/` | `sra2fastx.py`（已有 FASTQ 时 `--skip-convert`，并可改用 `--fastq-dir`） |
| 2 | `clean` | FASTQ → `00a_CleanData/` | `data_preprocessing_pipeline/clean-data.py`（fastp → seqkit → clumpify） |
| 3 | `hostref` | → `host_reference/hostdb/` | **可选**：给了 `--host-species`/`--host-fasta` 才跑；四通道下载 + 四索引构建 |
| 4 | `deplete` | `00a_CleanData/` → `00b_HostDepletion/` | `host_depletion.py`（Kraken2 → 比对 → Ribodetector） |
| 5 | `bbnorm` | `00b_HostDepletion/` → `00c_BBnorm/` | **可选**：`--bbnorm` 显式启用，默认不跑；`run_bbnorm.py`（target=70, mindepth=2） |

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--sra-dir` / `--fastq-dir` | — | 输入二选一（.sra 目录 或 FASTQ.GZ 目录） |
| `--outdir` / `-o` | `./preprocessed` | 输出根目录 |
| `--skip-convert` / `--skip-clean` / `--skip-deplete` | False | 跳过对应阶段（跳过 deplete 时要求宿主库已存在） |
| `--bbnorm` | False | 可选覆盖度归一化 → `00c_BBnorm/`；需 deplete 产物 |
| `--kraken2-db` + `--step2-index` | — | 已有宿主库（与 hostref 二选一） |
| `--host-species` / `--host-taxid` / `--host-fasta` | — | 现场构建宿主参考（`--host-taxid` 必填） |
| `--host-work-dir` | `<outdir>/host_reference` | 宿主参考构建输出目录 |
| `--host-download-source` | `auto` | 基因组下载通道 datasets→ngd→FTP 回退；gget 走 Ensembl |
| `--tool` | `bowtie2` | 比对工具 bowtie2 / hisat2 / minimap2 |
| `--seq-type` | `rna-short` | dna-short / rna-short / nanopore / pacbio |
| `--dedup` / `--skip-clumpify` / `--clumpify-memory` | — / False / `10g` | clean 阶段（fastp 可选去重·默认关、跳过 clumpify、堆内存） |
| `--threads` / `-t`、`--jobs` / `-j` | 40、4 | 总线程、并行任务数 |
| `--tmp` | — | 临时目录 |

### 输出

```
{outdir}/
  00_Converted/       ← .sra → .fastq.gz
  00a_CleanData/      ← fastp + clumpify + FASTA
  00b_HostDepletion/  ← 去宿主 + 去 rRNA 后 reads
  00c_BBnorm/         ← 覆盖度归一化 reads (仅 --bbnorm)
  preprocess.log      ← 全流程日志
  host_reference/     ← 仅 hostref 阶段 (genome/ + hostdb/ 四索引)
```

### 上下游衔接

- **← 上游:** `data_preprocessing_pipeline` / 本目录 `hostdb/` 提供宿主库（`--kraken2-db`/`--step2-index`），
  或由 `--host-species` 现场构建
- **→ 下游:** `00b_HostDepletion/`（或加了 `--bbnorm` 时的 `00c_BBnorm/`）作为
  `virome_discovery_pipeline/virome_pipeline.py --input_reads` 或
  `virome_analysis_pipeline/auto_known_virus.py --reads_dir`
- **说明:** 覆盖度归一化仅用于多样本深度差异大时的共组装前处理（默认关闭）；
  与发现管线共用输出根目录时，`report_pipeline.py` 会在报告里展示 `00c_BBnorm` 阶段
