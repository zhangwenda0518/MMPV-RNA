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

## 六大管线 / Six Pipelines

```
public_metadata_pipeline/     → 数据获取 (测序数据检索/下载/转换 ＋ 宿主参考基因组获取/建库)
       │
       ├── 测序数据: SRA/GSA 搜索 → 元数据清洗 → 批量下载 → 可视化 → SRA→FASTQ.GZ
       └── 宿主参考: 多通道下载·或用户 FASTA → Kraken2/Bowtie2/HISAT2/Minimap2 四索引
              │
              ├──▶ public_metadata_pipeline/eve_genome/ → 基因组 EVE 分析
              │      hostdb 前: 可选 EVE 掩码 (-–eve-mask)
              │      hostdb 后: 扫描 → 双侧判定 → RVDB → 元件库 → HTML (eve 阶段)
              │
              ▼
data_preprocessing_pipeline/  → 数据清洗 (质控 Fastp→Seqkit→Clumpify → 去宿主 Kraken2→比对→Ribodetector ＋可选 BBNorm 归一化)
              │
              ▼
virome_discovery_pipeline/    → 病毒发现（de novo）
       │
       ├── 组装 → 10工具并行鉴定 → 分层过滤 → COBRA延伸
       ├── 合并共组装 → CD-HIT参考聚类 → vclust去冗余 → 8工具分类 → 宿主预测
       ├── CheckV评估 → 两支路级联拯救 (A: 直通 → B: VSI延伸; C/D 已归档) → HQ vOTU catalog
       └── 09b KEEP 前: 候选+科标签 → eve_discriminate (16科 ML + 可选 Δbit)
                        EVE_like → KEEP 降 REVIEW
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

**可选补充模块（已归档 2026-10-11、2026-10-27 移出主仓，方法不成熟——
独立归档于 `差异表达分析` 项目，主链与各管线不受影响）**：

```
host_expression_pipeline/     → ⑧ 宿主表达分析（宿主侧：病毒阳性/阴性分组）【已归档】
       ├── 宿主 CDS 提取 (genome+GFF3 → gffread) → Salmon 定量 → 定量 QC
       ├── 病毒标签推导 (④ 检出表 → 阳性/阴性 + 连续载量; 合并 Core14 组织)
       ├── TEA-GCN 差异共表达 (主分析: 阳/阴双网络 → ΔzScore → 网络 z-QC 闸门 → DEG 限定 Louvain 模块; 内置样本量/失衡闸门)
       ├── DESeq2 差异表达 (交叉参考, 可选组织分层) → 载量×全基因组扫描 (剂量验证)
       └── HTML 汇总报告     ← 抽取自 goji-virome 差异表达分析核心链路
```

> 归档原因与恢复方法见该独立归档内的 `README.md` / `ARCHIVED_20261011.md`。输入全部取自主链现有产物
> （② clean reads、① 宿主参考/Core14、④ 检出表），输出落 `08_HostExpression/`。

---

## 统一 I/O 目录布局 / Unified I/O Layout (v3.1)

五条管线默认采用 `standard` 布局——**每条管线一个独立输出根，根内各自独立编号**。
存量项目（已有 checkpoint / 正在跑的任务）显式切回 `legacy` 布局（v3.0 目录名）：

```
<项目>/ 01_PublicData/ 02_Preprocessing/ 03_Discovery/ 04_Analysis/ 05_Phylo/ 90_Handoff/
        (EVE 产物: ① <host_reference>/eve/ · ③ <发现根>/09b_Analysis_Verify/virus_validation/eve/;
         旧 06_EVE/ 与 08_HostExpression/ 已随对应模块归档 2026-10-11)
```

管线边界交接可自动定位（③→④ 桥一条命令；①→⑤ Core14 原生直读）：

```bash
export MMPV_IO_LAYOUT=legacy            # 仅存量项目需要; 新项目默认 standard
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

# 同一流程换 Kun-peng 初筛 (自动把 Kraken2 宿主库 hashshard 成双后端库)
# python data_preprocessing_pipeline/host_depletion.py --classifier kunpeng     -k kraken2_db/ -x bowtie2_idx/ -I out/00a_CleanData/ -O out/00b_HostDepletion/     --tool bowtie2 --seq-type rna-short --rrna -t 40 -j 10

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

### 3. 内源病毒（EVE）分析（并入 ① 与 ③，不单独成管线）

EVE 拆成两段，分别挂在宿主库构建和病毒发现的 DNA 候选分析上；结果为各自管线的 HTML 增加章节。

**① 基因组侧** — `public_metadata_pipeline/eve_genome/`

| 环节 | 位置 | 说明 |
|---|---|---|
| 掩码 | `hostdb` **之前**（`--eve-mask`） | 先把 EVE 区间屏蔽，宿主索引里就不含 EVE 序列 |
| 分析汇总 | `hostdb` **之后**（`eve` 阶段，`--eve-scan`） | 扫描 → 双侧判定 → RVDB；可选元件库/consensus + 六框蛋白库 |
| 产物 | `<host_reference>/eve/` | `04_Summary/eve_report.html`、`05_Elements/eve_consensus_library.fa`、`.sixframe.faa` |
| 报告 | `Pipeline_Summary_Report.html` | **S7: EVE (genome)** 直接内嵌明细（卡片 + 判定分布/逐基因组/病毒科/位点/元件表）；`eve_report.html` 仅作可选独立导出 |

```bash
# 九阶段公共数据管线：hostref → hostdb → eve → report
python public_metadata_pipeline/public_data_pipeline.py --species <拉丁名> --taxid <ID> \
    --stage hostref,hostdb,eve,report --eve-scan --eve-elements \
    --work-dir <项目>/01_PublicData

# 或独立的宿主构建管线
python public_metadata_pipeline/build_host_pipeline.py --species <拉丁名> --taxid <ID> \
    --stage all --eve-scan --eve-elements --work-dir out/host
```

**③ 转录组侧** — `virome_discovery_pipeline/eve_discriminate/`

| 环节 | 位置 | 说明 |
|---|---|---|
| 判别 | `analysis_verify`（09b）内，CDD 验证之后、五层证据整合**之前** | 候选 fasta + taxonomy 科标签 → 16 科 ML（默认开） |
| 否决 | `integrate_rescue_evidence.py --eve-calls` | `EVE_like` 且 verdict=KEEP → 降为 REVIEW（只降不删） |
| 可选第二路 | `--eve-faa` + `--eve-exo-dmnd` | Δbit-score 双库竞争；EVE 库可直接用 ① 产出的 `.sixframe.faa` |
| 产物 | `09b_Analysis_Verify/virus_validation/eve/` | `eve_calls.tsv`、`eve_report.html` |
| 报告 | `10_Reports/pipeline_report.html` | 09b 卡片新增 **EVE discrimination (before KEEP)** 区块（模型数/EVE_like/否决数 + 链接） |

```bash
python virome_discovery_pipeline/virome_pipeline.py --stage analysis_verify \
    --output_dir <发现项目根> \
    --eve-faa <host_reference>/eve/05_Elements/eve_consensus_library.sixframe.faa \
    --eve-exo-dmnd <ICTV 蛋白库>.dmnd        # 可选：启用 Δbit-score
```

**为什么在 KEEP 前而不是后**：KEEP 由 `integrate_rescue_evidence.py` 判定，`class_KEEP.fasta`
直接从它抽取。放在 KEEP 之后只能给已进目录的 EVE 贴标签；放在之前，EVE 证据才能参与裁决
（`EVE_like` → REVIEW）。只降不删，保留人工复核空间。

**模型/数据**：16 科模型与阈值在 `database/eve_ml/`；训练脚本在 `database/eve_ml/build/`。
**基准可视化**：`database/eve_ml/benchmark/EVE_BENCHMARK.html`（4 主图 + 7 补充图 + Tier A/B/Δbit 表）。

> 旧 EVE 模块（`endogenous_virus_pipeline/` + `eve_distinguish/`）已于 2026-10-11
> 归档至 `archive/20261011_endogenous_virus_pipeline/`，保留用于历史结果复现。

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

# 提交 GUI 已移交外部本地平台 (D:/桌面/植物病毒分析平台); 管线内仍有独立 GUI:
# python virome_submission_pipeline/submission_gui.py
```

---

## 公共数据管线阶段 / Public Metadata Pipeline Stages

### SRA/GSA 公共数据获取

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `search` | gsa_sra.search.py | NCBI SRA + CNCB GSA 双引擎物种检索 → SRA_GSA_Merged_Final.csv |
| 2 | `info` | gsa_sra.info.py | SRA XML解析 + GSA爬虫 + AI元数据清洗 → Global_Unified_Metadata_Core14.csv (14列) |
| 3 | `down` | gsa_sra.down.py | aria2c/wget/prefetch 双协议下载 → 原始 FASTQ/SRA |
| 4 | `convert` | sra2fastx.py | 批量 SRA → FASTQ.GZ：默认 sracha（快 5–13x，严格校验失败自动降级 `--no-strict`），并可回退 fasterq-dump/fastq-dump |
| 5 | `plot` | gsa_sra.plot.py | 时间/组织/地区/机构 SCI级6面板可视化 |
| 6 | `hostref` | download_host_genome.py | 宿主基因组下载 (datasets/ngd/FTP/gget 四通道; `--host-fasta` 用户基因组跳过下载) → `host_reference/genome/` |
| 6 | `hostdb` | build_hostbase.py | Kraken2/Bowtie2/HISAT2/Minimap2 四索引构建 → `host_reference/hostdb/`（`--tool kraken2,kunpeng` 可追加 Kun-peng 分片库；默认由 Kraken2 库 hashshard 原地转换） |
| 7 | `report` | generate_report.py | 交互式 HTML 报告 + `sample_handoff.csv` 交接清单 (Run→FASTQ 路径+元数据, 交棒数据清洗管道) |

### 宿主参考数据库构建

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `genome-down` | public_metadata_pipeline/download_host_genome.py | 多通道下载参考基因组 (datasets/ngd/FTP/gget · 或 --genome-fasta 用户提供) + GFF3 + 细胞器基因组 |
| 2 | `hostdb` | public_metadata_pipeline/build_hostbase.py | Kraken2/Bowtie2/HISAT2/Minimap2 四种索引构建 (+ 可选 Kun-peng 分片) |

### 统一预处理入口 (从 .sra 一条命令到 clean reads)

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `convert` | sra2fastx.py | .sra → FASTQ.GZ（可选，已有 FASTQ 时 `--skip-convert`） |
| 2 | `clean` | data_preprocessing_pipeline/clean-data.py | Fastp → Seqkit → Clumpify → `00a_CleanData/` |
| 3 | `hostref` | download_host_genome.py + build_hostbase.py | 可选：宿主基因组四通道下载 + 四索引构建 |
| 4 | `deplete` | data_preprocessing_pipeline/host_depletion.py | 分类 (Kraken2/Kun-peng) → 比对 → Ribodetector → `00b_HostDepletion/` |
| 5 | `bbnorm` | data_preprocessing_pipeline/run_bbnorm.py | **可选**（`--bbnorm`，默认不跑）：覆盖度归一化 → `00c_BBnorm/` |

入口脚本: `python public_metadata_pipeline/preprocess_unified.py --sra-dir ... --outdir ...`
（产出目录名与 `data_preprocessing_pipeline` / 发现管线一致，可直接作为 `--input_reads`）

---

## 发现管线阶段 / Discovery Pipeline Stages

| # | 阶段 | 脚本 | 功能 |
|---|------|------|------|
| 1 | `clean` | data_preprocessing_pipeline/clean-data.py | Fastp QC + Seqkit FASTA转换 + Clumpify 聚类重排（提升压缩率） |
| 2 | `deplete` | data_preprocessing_pipeline/host_depletion.py | 分类 (Kraken2/Kun-peng) → Bowtie2/HISAT2/Minimap2 → rRNA去除 |
| 3 | `assembly` | assembly_pipeline.py | MEGAHIT / rnaviralSPAdes / Penguin 组装 |
| 4 | `identification` | virus_identification.py | 10工具并行鉴定 (Genomad+Diamond BLASTX+RdrpCatch+ViraLM+VirBot+VirSorter2+ViralVerify+VirHunter+Metabuli+Viroid BLASTN) |
| 5 | `filter` | filter_virus.py + graph_propagation.py | UniProt-strict 高置信过滤 (CDD 层已移除, CDD 仅用于 09b rescue 验证) + 可选 graph 子步骤: 对无证据保留结果二次裁决 (rescue/demote) → COBRA 输入; 思想源自 GMW, 2026 (默认关闭, graph.enabled 开启) |
| 6 | `cobra` | cobra_pipeline.py | BWA-MEM2 → CoverM → COBRA-Meta 单样本延伸 (候选源 = filter 主链镜像, 恢复 filter→COBRA 契约) |
| 7 | `merge` | (内置, Flye) | 多样本 Flye 共组装 (co-assembly 模式) |
| 8 | `cluster` | cluster_pipeline.py | CD-HIT参考引导预聚类 + vclust Leiden聚类 |
| 9 | `taxonomy` | virus_classifier.py + R | 8工具分类 (Genomad+Metabuli+CAT+Diamond LCA+VITAP+MMseqs2+ACVirus+vConTACT3) 加权投票 → 8级 taxonomy |
| 10 | `host` | run_host_prediction.py | ICTV > RNAVirHost > PhaBOX2 三级宿主预测 + 植物专属后过滤 |
| 11 | `checkv` | (内置) | CheckV 完整性评估 |
| 12 | `rescue` | rescue_pipeline.py | 两支路级联拯救 (CheckV ≥90% 直通 → Virseqimprover reads 延伸); 原分支 C/D 已归档 (archive/20261011_rescue_branch_cd/) |
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
| 3 | `extract` | batch_virus_variants.py --align_only + batch_extract_reads.py | **比对与制备**: virus-fasta/bam 索引比对 + all_summary 重算 + consensus + 靶向 reads 提取 (`virus_reads/`) |
| 4 | `variants` | batch_virus_variants.py --variants_only | **纯变异三件套**: FreeBayes/iVar/LoFreq VCF + SnpEff 注释 + SnpGenie dN/dS (复用 extract 的 BAM) |
| 5 | `post` | virus_vcf_pipeline.py | VCF可视化 + SnpEff宏 + MAF 瀑布 + SnpGenie dN/dS |
| 6 | `full` | batch_virus_full.py → virus-full.py (+ utils/extract_full_fasta.py) | **全长基因组构建 (OmniVirusAssembler V9.0)**: 12步精炼 — SHIVER-like Divine Fusion 实心骨架 → Reads级迭代抛光 (单碱基精度) → gmcloser + abyss-sealer 双引擎补洞 → 环化检测 (Circular=True); 尾段默认集成原 extract 的最长contig提取+参考N-fill (`--no_full_extract` 可关) |
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
| 6 | `gui` | submission_gui.py (管线内) | PySide6 桌面 GUI: 交互式编辑/验证/导出提交文件 (外部平台承接桌面交互) |
| 6 | `report` | report_html.py | 交互式HTML全表编辑报告 |

---

## 目录结构 / Directory Structure

```
MMPV-RNA/
├── data_preprocessing_pipeline/    # 数据清洗 (质控 clean-data.py → 去宿主 host_depletion.py; 含一条龙入口)
├── virome_discovery_pipeline/      # 病毒发现 (de novo, 15 阶段)
│   ├── virome_pipeline.py          # 主编排器 (09b KEEP 前接 EVE 判别)
│   ├── eve_discriminate/           # EVE 判别模块: no_ref_genome/ (ML/Δbit/无监督) + with_ref_genome/ (整合证据) + 报告
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 辅助工具 (组装/鉴定/COBRA统计, Sankey)
│
├── virome_analysis_pipeline/       # 已知病毒深度分析 (9 阶段)
│   ├── auto_known_virus.py         # 分析编排器
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 辅助工具
│
├── public_metadata_pipeline/       # 公共数据获取 (9 阶段)
│   ├── public_data_pipeline.py     # 公共数据编排器 (检索/元数据/下载/可视化 + hostdb/eve)
│   ├── build_host_pipeline.py      # 宿主参考编排 (多通道下载/用户FASTA → 四索引 → EVE 报告)
│   ├── eve_genome/                 # 基因组 EVE: eve_screen/element_db/report + 六框蛋白库
│   ├── preprocess_unified.py       # 全链路预处理入口 (convert→clean→hostref→deplete→[bbnorm])
│   ├── doc.md                      # 完整流程文档
│   └── utils/                      # 共享工具模块
│
├── virome_phylo_pipeline/          # 系统发育/进化/群体遗传 (18 stage, 8 大模块)
│   ├── phylo_pipeline.py           # 总管线 (BEAST 定年 + pypopart + 正选择)
│   └── doc/ 权威体系 (RUN_GUIDE/STAGE_REFERENCE/METHODS)
├── virome_submission_pipeline/     # 提交管线 (GenBank/CNCB)
│   ├── submission_pipeline.py      # 主编排器
│   ├── sequin_builder.py           # Sequin 构建器
│   └── ...                         # 元数据/报告/拓扑分析
│
├── database/eve_ml/                # EVE ML 资产: 16 科模型 + index.json + FEATURE_SPEC + build/ 构建脚本
│
├── archive/20261011_endogenous_virus_pipeline/  # 【已归档 2026-10-11】旧 EVE 模块 (被 ①③ 内的 EVE 模块取代, 保留历史复现)
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
├── biosoft/                        # 第三方工具 (脚本/JAR, 无需编译)
│   ├── VirBot/VirBot.py
│   ├── virhunter/predict_cpu.py + weights/
│   ├── ViReMa/ViReMa.py
│   └── snpEff/snpEff.jar + config + scripts/
│
├── doc/                            # 文档 (17篇 + 方法章总模板 + 依赖清单)
│   ├── METHODS_TEMPLATE.md         # 五大管道论文级方法章合并版 (SCI投稿直接取用)
│   └── MMPV_dependencies.md        # 外部依赖清单 (工具/数据库/路径)
├── scripts/                        # 脚本引用检查 (一次性脚本已归档到 archive/20261025_cleanup/scripts_oneoff/)
│   ├── README.md
│   └── check_script_refs.py        # 管线脚本引用完整性检查
├── MMPV-RNA-paper/                 # SCI 论文稿与引用库 (citations/ + figures/)
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
- **两支路级联拯救**: A: CheckV≥90% 直通 → B: Virseqimprover reads 延伸，只保留有质量证据的产出 (原 C 直通筛选/D 零产出已归档: archive/20261011_rescue_branch_cd/)
- **宿主过滤**: 拯救前按宿主类别预过滤，节省70%+计算量
- **pixi一键部署**: 100个conda包精确版本管理，`pixi install`即可
- **断点续传**: 全部脚本支持 `--resume/--force`，大规模运行安全中止和恢复
- **基于分类层级的新颖性判断**: 无需BLASTN依赖，直接从taxonomy completeness判断新病毒
- **完整闭环**: 从公共数据挖掘到投稿图表，一站式完成
- **提交管线**: GenBank/CNCB 序列提交自动化, 支持 Sequin 构建和假设蛋白分析
- **AI 驱动层 (2026-10)**: 四层递进架构 — run digest 状态数字化 → AI 解读 (五报告统一接线 + **过程面板: 报告不止读数据, 兼载全过程与 AI 参与轨迹**) → 建议/候选 triage 实战标注 (`--preset ident` 鉴定证据矩阵, 只写 ai_* 列) → 自然语言副驾驶 + MCP server + 受控 agent 全链路编排 (`mmpv_ai.py agent`) + 项目级主报告 (`mmpv_ai.py report`, 含报告互链/digest 演化图/token 用量) + 溯源导出 (`mmpv_ai.py provenance` → 中英文方法学段落) + 事后补 AI (`mmpv_ai.py enrich`: 无 key 先跑, 配 key 后一键补解读, 幂等)。AI 只建议不执行, 离线完整可用; 入口 `python mmpv_ai.py --help`, 论证见 `doc/AI_DRIVEN_DESIGN.md`
- **AI 集成二期 (2026-10-28)**: stage 级自动摘要 (`ai_hooks` → `90_Handoff/ai/`)、`summarize/watch/backfill/doctor/eval` 子命令、失败诊断 (12 条错误签名 + 白名单校验的修复草稿)、参数建议、文档 BM25 检索 (`file:line`)、MCP 14 个只读工具、6 个 `skills/mmpv-*` 工作流包。摘要默认模板兜底, 无 key/断网可用; 数字 grounding 与 prompt 版本全程留痕。详见 `doc/AI_INTEGRATION.md`
- **报告模板层 (2026-10)**: 报告 HTML 由 Jinja2 统一渲染 (`mmpv_common/report_template.py` + `templates/`: base 骨架 + 组件宏), 过程面板/项目主报告/⑦ 提交报告已迁移; autoescape 修复元数据特殊字符破坏标记, StrictUndefined 缺变量 fail-loud, 离线零外部资源

---

## 数据库清单 / Required Databases

| 数据库 | 用途 | 来源 |
|--------|------|------|
| CheckV DB | 病毒完整性评估 | https://bitbucket.org/berkeleylab/checkv/ |
| geNomad DB | 病毒鉴定 | https://zenodo.org/records/14828026 |
| RVDB | 病毒参考序列 | https://fzer.github.io/rvdbtools/ |
| NR (diamond) | 蛋白过滤 | NCBI nr |
| NCBI virus ref | BLAST参考 | NCBI virus |
| ViralVerify HMM + NBC 模型表 | HMM 验证 + Naive Bayes 判定 | ViralVerify |
| VirSorter2 DB | 病毒分类 | VirSorter2 |
| ViraLM DB | DNABERT-2鉴定 | Google Drive (gdown) |
| Kraken2 + align DB | 宿主去除 | 本管线构建 (build_host_pipeline.py) |
| PhaBOX2 DB | 宿主预测 | PhaBOX2 |
| 植物/病毒拆分库 (hybrid) | `eve_genome` 双侧判定 | SwissProt 拆分 + id2div 映射 (`database/setup/mmpv_db_manifest.tsv`) |

> ViralVerify 的 NBC 模型表 `classifier_table_uniq.txt` 默认与 HMM 同目录；缺失时
> pyhmmer 路径自动退回 HMM-hit 判定（下载见 `DATABASE_SETUP.md` §2.4）。

**完整清单与部署方案**：全管线共依赖 50 项数据库（必需 24 / 可选 26，全量约 710GB）。
逐库的下载与建索引命令见 `DATABASE_SETUP.md`；**可复现的集中部署方案**（共享库根 + 软链接
兼容层 + 环境变量 + 校验）见 `doc/DB_DEPLOYMENT.md`，配套工具：

```bash
python database/setup/mmpv_db.py list                                 # 打印 50 项依赖清单
python database/setup/mmpv_db.py check                                # 部署前体检
python database/setup/mmpv_db.py verify                               # 校验现有库是否齐全
python database/setup/mmpv_db.py down --apply                          # 下载/构建可自动获取的库
python database/setup/mmpv_db.py adopt --db-from ~/database --apply   # 迁移到共享库根
python database/setup/mmpv_db.py link --apply                          # 建旧路径兼容软链接
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
| `doc/AI_DRIVEN_DESIGN.md` | AI 驱动四层架构: 设计论证 / 红线 / 上手指南 |
| `doc/DB_DEPLOYMENT.md` | 数据库部署方案（共享库根/软链接/环境变量/校验，可复现） |
| `SOFTWARE_VERSIONS.txt` | 全部第三方软件版本记录 |
| `pixi.toml` | 依赖配置 (可直接查看) |

---

## 引用 / Citation

Zhang W. et al. **MMPV: Massive Meta-mining of Plant Viruses**. 2026.

> If you use MMPV in your research, please cite the above reference.
