# MMPV 统一 I/O 目录布局设计（standard 布局规范 · 管线独立根版）

> 2026-09-27 定稿（同日第二版：由"项目单树全链编号"改为**管线独立根**方案）。
> 解决三个问题：①五管线目录命名体系各自为政（③ 的 `00a/02a/09a` 字母后缀、
> ① 的散落小写目录）；②②的产物寄居在③的输出根里，管线不独立；③改名会砸
> checkpoint。
>
> **结论：双布局并存，legacy（现行名）为默认，standard 经开关启用。
> standard 下每条管线一个独立输出根（`01_PublicData … 06_EVE`），
> 根内部各自独立编号（⑤ per-virus 模块按输出顺序 01-08 编号），互不占号段。
> 不开开关 = 一切照旧。**

---

## 1. 设计原则

1. **管线独立**：standard 下五条管线各有自己的输出根目录，`--output_dir`/
   `--work-dir` 直接指向本管线自己的根；管线之间只通过边界产物路径衔接
   （§3 契约），任何管线的输出目录都可以单独搬走、单独删除、单独重跑。
2. **独立顺序**：每个管线根内部从 `01` 起独立编号（② `01_CleanData…`，
   ③ `01_Assembly…13_Report`，④ `01_Detection…09_Report`，① 按 stage 顺序
   `01_Search…05_HostRef`），两管线之间号段互不相干。
3. **零破坏**：`legacy` 布局与 v3.0 目录名逐字节一致（`mmpv_common/tests/`
   冻结测试保护）。所有正在跑的任务、已有 checkpoint、GUI、文档默认不受影响。
4. **单一解析点**：目录名只活在 `mmpv_common/io_layout.py` 一张注册表里。
5. **布局传播靠环境**：编排器入口 `normalize_layout_env()` 把解析结果写回
   `MMPV_IO_LAYOUT`，子进程自动继承。
6. **边界留痕**：跨界产物登记 `90_Handoff/handoff_manifest.tsv`（项目根）。

## 2. standard 布局：六个独立输出根

```
<project_root>/
├── 01_PublicData/                 ← ① 公共数据管线的输出根 (--work-dir 指向这里)
│   ├── 01_Search/                   SRA_GSA_Merged_Final.csv + sra_raw/ gsa_raw/
│   ├── 02_Unified/                  Global_Unified_Metadata_Core14.csv (★全平台元数据权威)
│   ├── 03_RawData/                  ERR…_1.fastq.gz … (down/convert 产物)
│   ├── 04_Plots/                    检索可视化 6 面板
│   └── 05_HostRef/
│       ├── genome/                  all.genome.uniq.fasta (+GFF3/细胞器)
│       └── hostdb/                  Kraken2/Bowtie2/HISAT2/Minimap2 四索引
│
├── 02_Preprocessing/              ← ② 数据清洗管线的输出根 (--output_dir 指向这里)
│   ├── 01_CleanData/                1.fastp/ 2.fasta/ 3.clumpify/ logs/
│   ├── 02_HostDepleted/             ★ 去宿主去rRNA reads (全下游 reads 唯一来源)
│   ├── 03_BBnorm/                   (可选)
│   ├── preprocessing_summary.tsv  assembly_ready.list
│   └── preprocessing_report.html  preprocessing.log
│
├── 03_Discovery/                  ← ③ 病毒发现管线的输出根 (--output_dir 指向这里)
│   ├── 01_Assembly/                 {sample}/{sample}_{tool}.contig.fasta
│   ├── 02_Identification/           10 工具并行鉴定
│   ├── 03_Filter/                   UniProt-strict 过滤
│   ├── 04_COBRA/                    单样本延伸
│   ├── 05_CoAssembly/               Flye 共组装
│   ├── 06_CLUSTER/                  centroids/final_centroids.fasta (★③→④交接物)
│   ├── 07_Taxonomy/                 integrated/final_integrated_classification.tsv (★)
│   ├── 08_HostPrediction/  09_CheckV/
│   ├── 10_Rescue/                   HQ_plant_viruses.fasta + final_judgement_table.tsv (★)
│   ├── 11_ViromeAnalysis/  12_AnalysisVerify/  13_Report/
│   └── (③ 若在本管线内跑 clean/deplete 阶段, 产物落 01_CleanData/02_HostDepleted,
│        与 legacy 同语义; 推荐做法仍是 ② 独立跑完后 --input_reads 接入)
│
├── 04_Analysis/                   ← ④ 已知病毒分析管线的输出根 (--output_dir 指向这里)
│   ├── 01_Detection/                summary/all_viruses.best.summary.tsv
│   ├── 02_Filtering/                high_conf.summary.tsv
│   ├── 03_Variants/  04_PostAnalysis/
│   ├── 05_Assembly/  06_Extraction/  ← 06_Extraction/<病毒>/<样本>.full.fasta (★④→⑤交接物)
│   ├── 07_Similarity/  08_DVG/  09_Report/
│   └── tmp/
│
├── 05_Phylo/                      ← ⑤ 系统发育管线的输出根 (--output_dir 指向这里)
│   └── <virus>/                     per-virus 工作目录, 内部 8 模块按输出顺序编号:
│     ├── 01_data/                   prep+metadata+online+clean (dates.csv/三通道表/clean.fasta)
│     ├── 02_phylogeny/              align+splitstree+tree (mafft.aln.fasta, iqtree.treefile)
│     ├── 03_popgen/                 popgen+host (pypopart, host_analysis/)
│     ├── 04_recomb/                 rdp5 (事件CSV + masked 比对)
│     ├── 05_select/                 capheine (cawlign→iqtree→HyPhy→DRHIP)
│     ├── 06_time/                   clock+beast+gene_dating (时间信号门禁/定年)
│     ├── 07_geography/              phylogeo+geo+geo_paths (CTMC 系统地理/传播路径)
│     ├── 08_report/                 phylo_report.html + phylo_summary.csv
│     └── .checkpoints/              断点标记 (不编号)
│
├── 06_EVE/                        ← ⑥ 内源性病毒管线的输出根 (eve_screen -o 指向这里;
│   │                                内部编号两布局一致, 本就符合编号惯例)
│   ├── 01_Loci/<NAME>/              Stage1 滑窗 blastx 发现 (fna/chunks/s1_hits/loci.bed/loci.fa)
│   ├── 02_Verdict/<NAME>/           Stage2 植物/病毒双侧 bitscore 判定 (s2_verdict)
│   ├── 03_RVDB/<NAME>/              Stage3 RVDB 深度归属 (cand3.bed/s3_rvdb)
│   ├── 04_Summary/                  <NAME>_eve_summary.tsv + kingdom_summary.tsv +
│   │                                family_by_genome.tsv + eve_report.html (HTML 报告)
│   └── (eve_distinguish 的判别产物建议放 03_Discovery/14_EVE_Distinguish/, 见 §3 B7)
│
└── 90_Handoff/                    ← 跨管线交接层 (项目根, 两布局统一位置)
    ├── analysis.reference.fasta     桥脚本 ③→④ 产物
    ├── analysis.ref_info.tsv
    └── handoff_manifest.tsv         边界台账 (boundary/key/path/produced_by/produced_at)
```

命名规则：管线根 = `0N_PascalCase`（N=管线号 1-5）；根内 = `NN_PascalCase`
独立编号；⑤ 的内部结构沿用其 per-virus 惯例不编号。

## 3. 五边界契约（管线根之间的衔接）

| # | 边界 | 生产者 → 产物（standard 路径） | 消费者入口参数 |
|---|------|------------------------------|----------------|
| B1 | ①→② | `01_PublicData/03_RawData/` + `01_PublicData/02_Unified/Core14` + `01_PublicData/05_HostRef/hostdb/` | `data_preprocessing.py --input_reads --host_db` |
| B2 | ②→③ | `02_Preprocessing/02_HostDepleted/`（+ `assembly_ready.list`） | `virome_pipeline.py --input_reads` |
| B3 | ③→④ | `03_Discovery/06_CLUSTER/centroids/final_centroids.fasta` + `07_Taxonomy/integrated/final_integrated_classification.tsv` → 桥 → `90_Handoff/analysis.reference.{fasta,tsv}` | `auto_known_virus.py --reference --ref_info --reads_dir <②的02_HostDepleted>` |
| B4 | ④→⑤ | `04_Analysis/06_Extraction/` | `phylo_pipeline.py --extract_dir` |
| B5 | ①→⑤ | `01_PublicData/02_Unified/Global_Unified_Metadata_Core14.tsv` | `phylo_pipeline.py --metadata`（`data_collector.load_sample_metadata` 原生直读，无需转换器） |
| B6 | ①→⑥ | `01_PublicData/05_HostRef/genome/all.genome.uniq.fasta` | `eve_screen.py -g <genome.fa>` 或批量 `-B batch.tsv`（EVE 与 reads 主链相互独立） |
| B7 | ③↔⑥ | ⑥ 的判别器吃 ③ 的发现输出，产出 `dna_vs_eve_filter.tsv`（建议落 `03_Discovery/14_EVE_Distinguish/`） | `eve_distinguish/run_all.sh -D <项目>/03_Discovery -A <1kp_assemblies>` |

③→④ 桥（两布局通吃，无需开关）：

```bash
python virome_discovery_pipeline/utils/discovery2analysis.py \
    --from-discovery <项目>/03_Discovery \
    --output_prefix  <项目>/90_Handoff/analysis
```

## 4. legacy ↔ standard 映射总表

standard 相对路径均相对**本管线自己的根**。

| 逻辑键 | legacy（相对共享根） | standard（相对管线根） | 归属根 |
|---|---|---|---|
| m_search | `search/` | `01_Search/` | 01_PublicData |
| m_info | `info/` | `02_Unified/` | 01_PublicData |
| m_down | `down/` | `03_RawData/` | 01_PublicData |
| m_plot | `plot/` | `04_Plots/` | 01_PublicData |
| h_genome | `host_reference/genome/` | `05_HostRef/genome/` | 01_PublicData |
| h_hostdb | `host_reference/hostdb/` | `05_HostRef/hostdb/` | 01_PublicData |
| d_clean | `00a_CleanData/` | `01_CleanData/` | 02_Preprocessing（③内跑则为 03_Discovery） |
| d_hostdep | `00b_HostDepletion/` | `02_HostDepleted/` | 同上 |
| d_bbnorm | `00c_BBnorm/` | `03_BBnorm/` | 02_Preprocessing |
| d_asm | `01_Assembly/` | `01_Assembly/` | 03_Discovery |
| d_ident | `02a_Identification/` | `02_Identification/` | 03_Discovery |
| d_filter | `02b_Filter/` | `03_Filter/` | 03_Discovery |
| d_cobra | `03a_COBRA/` | `04_COBRA/` | 03_Discovery |
| d_merge | `03b_MergeSamples/` | `05_CoAssembly/` | 03_Discovery |
| d_cluster | `04_CLUSTER/` | `06_CLUSTER/` | 03_Discovery |
| d_centroids | `04_CLUSTER/4_centroids/` | `06_CLUSTER/centroids/` | 03_Discovery |
| d_taxonomy | `05_Taxonomy/` | `07_Taxonomy/` | 03_Discovery |
| d_host_pred | `06_HostPrediction/` | `08_HostPrediction/` | 03_Discovery |
| d_checkv | `07_Checkv/` | `09_CheckV/` | 03_Discovery |
| d_rescue | `08_Rescue/` | `10_Rescue/` | 03_Discovery |
| d_analysis | `09a_Virome_Analysis/` | `11_ViromeAnalysis/` | 03_Discovery |
| d_verify | `09b_Analysis_Verify/` | `12_AnalysisVerify/` | 03_Discovery |
| d_reports | `10_Reports/` | `13_Report/` | 03_Discovery |
| a_detect | `01_detection/` | `01_Detection/` | 04_Analysis |
| a_filter | `02_filtering/` | `02_Filtering/` | 04_Analysis |
| a_variants | `03_variants/` | `03_Variants/` | 04_Analysis |
| a_post | `04_post_analysis/` | `04_PostAnalysis/` | 04_Analysis |
| a_assembly | `05_assembly/` | `05_Assembly/` | 04_Analysis |
| a_extract | `06_extraction/` | `06_Extraction/` | 04_Analysis |
| a_similarity | `07_similarity/` | `07_Similarity/` | 04_Analysis |
| a_dvg | `08_dvg/` | `08_DVG/` | 04_Analysis |
| a_report | `09_report/` | `09_Report/` | 04_Analysis |
| p_root | `phylo_results/` | `05_Phylo/` | 项目根 |
| ph_data | `data/` | `01_data/` | 05_Phylo/<virus> |
| ph_phylogeny | `phylogeny/` | `02_phylogeny/` | 05_Phylo/<virus> |
| ph_popgen | `popgen/` | `03_popgen/` | 05_Phylo/<virus> |
| ph_recomb | `recomb/` | `04_recomb/` | 05_Phylo/<virus> |
| ph_select | `select/` | `05_select/` | 05_Phylo/<virus> |
| ph_time | `time/` | `06_time/` | 05_Phylo/<virus> |
| ph_geography | `geography/` | `07_geography/` | 05_Phylo/<virus> |
| ph_report | `report/` | `08_report/` | 05_Phylo/<virus> |
| e_loci | `01_Loci/` | `01_Loci/`（同） | 06_EVE |
| e_verdict | `02_Verdict/` | `02_Verdict/`（同） | 06_EVE |
| e_rvdb | `03_RVDB/` | `03_RVDB/`（同） | 06_EVE |
| e_summary | `04_Summary/` | `04_Summary/`（同） | 06_EVE |
| x_handoff | `90_Handoff/` | `90_Handoff/`（同） | 项目根 |

> ⑤ 模块编号顺序 = 主线分析顺序（与 `STAGE_REFERENCE.md` 8 大模块表一致）：
> data(prep→clean) → phylogeny(align,splitstree,tree) → popgen(popgen,host) →
> recomb(rdp5) → select(capheine) → time(clock,beast,gene_dating) →
> geography(phylogeo,geo,geo_paths) → report。⑤ 的两个信息字典
> （`STAGE_MODULE_DIR`/`STAGE_GROUPS`）保存的是逻辑模块名（即 `ph_dir()` 的键），
> 两布局通用，不随布局改名。

## 5. 启用方式与传播机制

```bash
export MMPV_IO_LAYOUT=standard        # 或各编排器 --io-layout standard (CLI 优先)

# ① --work-dir 直接指 ① 自己的根
python public_metadata_pipeline/public_data_pipeline.py \
    --work-dir <项目>/01_PublicData --stage all ...

# ② --output_dir 直接指 ② 自己的根
python data_preprocessing_pipeline/data_preprocessing.py \
    --input_reads <项目>/01_PublicData/03_RawData \
    --host_db <项目>/01_PublicData/05_HostRef/hostdb \
    --output_dir <项目>/02_Preprocessing --stage all

# ③ --output_dir 直接指 ③ 自己的根, reads 从 ② 根接入
python virome_discovery_pipeline/virome_pipeline.py \
    --input_reads <项目>/02_Preprocessing/02_HostDepleted \
    --output_dir <项目>/03_Discovery --stage assembly…report

# ③→④ 桥 + ④
python virome_discovery_pipeline/utils/discovery2analysis.py \
    --from-discovery <项目>/03_Discovery --output_prefix <项目>/90_Handoff/analysis
python virome_analysis_pipeline/auto_known_virus.py \
    --reference <项目>/90_Handoff/analysis.reference.fasta \
    --ref_info <项目>/90_Handoff/analysis.ref_info.tsv \
    --reads_dir <项目>/02_Preprocessing/02_HostDepleted \
    --output_dir <项目>/04_Analysis --stage all

# ⑤
python virome_phylo_pipeline/phylo_pipeline.py \
    --extract_dir <项目>/04_Analysis/06_Extraction \
    --metadata <项目>/01_PublicData/02_Unified/Global_Unified_Metadata_Core14.csv \
    --output_dir <项目>/05_Phylo --stage all

# ⑥ EVE (与 reads 主链相互独立, 输入=① 的宿主参考基因组)
python endogenous_virus_pipeline/eve_screen.py \
    -g <项目>/01_PublicData/05_HostRef/genome/all.genome.uniq.fasta \
    -n <物种名> -o <项目>/06_EVE -t 40

# ⑥ 的 HTML 报告 (在 --merge 之后; 补齐与 ①-⑤ 一致的末阶段 HTML 闭环)
python endogenous_virus_pipeline/eve_report.py -o <项目>/06_EVE

# ⑥ 的后运行判别 (吃 ③ 发现输出, 建议产物落 03_Discovery/14_EVE_Distinguish/)
bash endogenous_virus_pipeline/eve_distinguish/run_all.sh \
    -D <项目>/03_Discovery -A <1kp_assemblies> -t 32
```

优先级：CLI `--io-layout` > 环境变量 `MMPV_IO_LAYOUT` > `legacy`。

## 6. 迁移与回滚

- **不迁移**：legacy 是默认，老项目继续 legacy 跑完；新旧项目共存。
- **切换时机**：只允许在**新项目开跑前**切。跑到一半切布局 = checkpoint
  判定失效（产物目录名变了），必须 `--force` 重跑或手工搬目录。
- **回滚**：不开开关即回滚；代码层回滚 = revert `mmpv_common/io_layout.py`。
- **手工搬旧项目**：按 §4 映射表先按管线归组（`00a/00b…→02_Preprocessing/`，
  `01_Assembly…10_Reports→03_Discovery/`，`01_detection…→04_Analysis/`），
  再改各自根内目录名；文件内容零依赖目录名。

## 7. 实现清单（2026-09-27 落地，含第二版修订）

| 文件 | 改动 |
|---|---|
| `mmpv_common/io_layout.py` | **新增**：目录注册表（standard=管线独立根+根内独立编号）+ 解析/规范化 + 三个构造器 + 边界定位/台账 |
| `mmpv_common/tests/test_io_layout.py` | **新增**：15 用例（legacy 冻结 / 键集一致 / 优先级 / 构造器两布局期望 / 定位 / 台账幂等） |
| `data_preprocessing_pipeline/data_preprocessing.py` | 目录名 → `layout_dir_name()`；`--io-layout`；report 阶段显式传 clean/deplete 目录（preprocess_report 默认值是 legacy 名） |
| `virome_discovery_pipeline/virome_pipeline.py` | `self.d` → `build_discovery_dirs()`；`--io-layout`；rescue 阶段 00b 硬编码检查改布局感知 |
| `virome_discovery_pipeline/utils/discovery2analysis.py` | `--from-discovery` 自动定位（两布局 + 三代旧名）+ 台账登记 |
| `virome_analysis_pipeline/auto_known_virus.py` | 9 目录 → `build_analysis_dirs()`；`--io-layout` |
| `public_metadata_pipeline/public_data_pipeline.py` | dirs/宿主参考路径 → 布局；`--io-layout` |
| `public_metadata_pipeline/preprocess_unified.py` | 00a/00b/00c → `layout_dir_name()`；`--io-layout` |
| `virome_phylo_pipeline/phylo_pipeline.py` | 69 处 `os.path.join(out_dir,"模块",…)` → `ph_dir()`（布局感知模块名）；`--io-layout`；`STAGE_MODULE_DIR`/`STAGE_GROUPS` 保留逻辑名不动 |
| `endogenous_virus_pipeline/eve_scan_core.py` | `STAGE_DIRS` → `io_layout` e_* 键（带脱离仓库兜底回退）；standard 输出根约定 06_EVE |
| `endogenous_virus_pipeline/eve_report.py` | **新增**：⑥ 的 HTML 报告生成器（merge 产物 → 自包含 HTML），补齐六管线"末阶段出 HTML"闭环 |
| 衔接审计（2026-09-27 第二轮） | **13 个 worker 脚本补接线**：③ report_pipeline/build_ref_info/filter_virus/rescue_pipeline/virome_analysis/virome_pipeline(cobra 回退+日志目录)、④ generate_pipeline_report(含 33 处相对子路径+陈旧 09_dvg 修复)/generate_pipeline_report_full/panvirome_×2/virus_detection_summary、② preprocess_report、① generate_report/preprocess_unified(hostref 根)、⑤ phylo(--output_dir 默认)。审计器：`scripts/audit/check_io_handoffs.py`（0 fail / 23 warn，warn 均为帮助文本措辞） |
| 衔接审计（2026-09-28 第三轮，含 report/HTML 链路实测） | **⑦ 提交管线纳入审计范围并接线**：submission_pipeline(suvtk 输出探测/KEEP FASTA 探测/ref_info 定位)、fix_sqn(08_Rescue×3)。**⑤ utils 接线**：clock_analysis(6 处 time/phylogeny 路径)、data_collector(prep 的 data 目录)、report_builder(MODULES 注册表 7 个 dir 字段 + 注册表重排为编号序)、beast1_bridge(IQ-TREE 报告探测)。**实测**：③ report_pipeline 空 standard 树实跑产出 `13_Report/pipeline_report.html`（目录引用全为 standard 名）；④ generate_pipeline_report 合成 summary 实跑，standard 路径读取入 HTML；② preprocess_report 合成 fastp/seqkit 数据实跑出 HTML。审计器 SCOPES 补 ⑦ |

验证：`python mmpv_common/tests/test_io_layout.py`（15/15 OK）+ 全部改动文件
`py_compile` 通过 + `discovery2analysis --from-discovery` 端到端冒烟通过。
