# endogenous_virus_pipeline — 内源性病毒元件 (EVE) 筛查流程文档

## 数据流总览

```
eve_screen.py (主编排器)
  │
  ├─ discover  (Stage1) → 滑窗 50kb → diamond blastx vs 病毒参考蛋白
  │                     → 坐标还原 → 300nt 间隔合并位点 → 抽 locus 序列
  ├─ verdict   (Stage2) → locus vs 植物/病毒拆分库 (SwissProt 混合)
  │                     → 两侧 bitscore 比较 → viral_supported / host_like / undetermined
  ├─ annotate  (Stage3) → viral_supported + undetermined vs RVDB (深度科/属归属)
  ├─ viroid    (可选)   → locus vs 类病毒库 blastn
  └─ merge              → 跨基因组汇总 (kingdom_summary / family_by_genome)
```

输入是**宿主参考基因组 FASTA**（不是测序 reads）——与其余 6 条管线的数据流截然不同，
因此独立成顶级管线。宿主基因组来自 `public_metadata_pipeline --stage hostref`
（`host_reference/genome/all.genome.uniq.fasta`）或用户直接提供。

## 与 virome_discovery_pipeline 的分工（含子模块 `eve_distinguish/`）

| | 主管线 `eve_screen.py` | 子模块 `eve_distinguish/` |
|---|---|---|
| 输入 | 完整宿主参考基因组 | 病毒组 de novo 候选 contig |
| 手段 | 正向搜索：病毒蛋白库 blastx 扫基因组 | 结构退化 × 单位点基因座架构双通道判别 |
| 产出 | 基因组坐标级 EVE 位点 + 判定 | contig 级 verdict + `dna_vs_eve_filter.tsv`（action 列） |
| 运行方式 | `eve_screen.py` 编排的 Stage1→3 | **独立后运行脚本**，不是对方的 stage |

两条线**互补不重复**：主管线回答"宿主基因组里哪里有病毒插入、是什么"；
`eve_distinguish/` 回答"病毒组组装的候选 contig 是真 DNA 病毒还是退化 EVE"。
下游可交叉验证（`eve_distinguish` 的 review 类的升级路径之一即宿主基因组比对）。

### 为什么 `eve_distinguish` 不做成管线的 stage

v3 之前它挂在 `virome_pipeline.py --stage eve` 上，输入/输出全靠主管道用环境变量注入
（`QUERY`/`EVIDENCE`/`OUT_DIR`/`ASM_ROOT`）。实际用起来两个问题：

1. 主流程多一条"可选但影响结论"的岔路——忘记开就整类假阳性静默漏掉
   （`host_contamination_likely` 检测不到时不报错）；
2. 判别逻辑散在两处（verdict→action 映射在 python 里，参数默认值在 shell 里），
   改一头忘另一头。

v5 起它搬进本管线，作为**跑完发现管道再单独调用**的后运行脚本，自带 CLI 与完整参数校验
（缺 `-A` 直接报错而不是静默降级）；`virome_pipeline.py` 里的 eve 阶段接入已整体移除。

```bash
# 指向一次发现管道的输出目录: 候选/证据表按阶段目录约定自动定位, 产物落 <DISCOVERY>/08b_EVE_Distinguish
bash endogenous_virus_pipeline/eve_distinguish/run_all.sh \
    -D /path/to/onekp-virus -A /path/to/1kp/assemblies -t 32
```

输入约定（`-D` 时自动找，也可 `-q`/`-e` 显式给）：`10_Reports/eve_candidates.fasta`
（或 `eve_screen/evescan_query.fasta`）与 `10_Reports/rescue_evidence_scored.tsv`。
方法与踩坑记录见 [`eve_distinguish/README.md`](eve_distinguish/README.md)。

---

## 1. eve_screen.py — 主编排器

**用途:** 三阶段 EVE 筛查调度（`discover → verdict → annotate` + `merge`）。
支持单基因组 (`-g`)、批量 TSV (`-B`，路径可直接填 `.tar.gz`/`.gz` 基因组包)；
基因组级 ProcessPoolExecutor 并行（总核数 = `-J × -t`）；断点续传以
`04_Summary/<NAME>_eve_summary.tsv` 存在为完成标记，已完成基因组自动跳过；
单基因组失败不影响整批（结果以 ok/error 汇总返回）。

### 阶段说明

| 阶段 | 输入 | 输出 | 说明 |
|------|------|------|------|
| discover (1) | 基因组 FASTA | `01_Loci/<NAME>/` | 见下 |
| verdict (2) | loci + hits BED | `02_Verdict/<NAME>/s2_verdict.tsv` | 双侧 bitscore 判定 |
| annotate (3) | 候选序列 | `03_RVDB/<NAME>/s3_rvdb.tsv` | 无候选也写空文件（断点标记） |
| merge | 全部 summary | `kingdom_summary.tsv` 等 | 全阶段跑完自动执行，`--no-merge` 关 |

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-g/--genome` | — | 单基因组 FASTA（`.fa/.fasta/.fna/.gz/.tar.gz`） |
| `-n/--name` | 文件名 | 基因组名（自动清洗非法字符） |
| `-B/--batch` | — | 批量 TSV：每行 `NAME<TAB>/path/to/genome.fa`（路径可为 `.tar.gz`/`.gz` 包；**NAME 必须唯一**，重名开跑前报错中止） |
| `--list` | 关 | 干跑：只列任务不执行 |
| `-o/--outdir` | `eve_results` | 输出根目录 |
| `-t/--threads` | 40 | 每基因组线程数 |
| `-J/--jobs` | 1 | 并行基因组数（总核数 = jobs × threads；256 核机推荐 `-J5 -t40`） |
| `--stage` | `all` | `discover/verdict/annotate` 逗号组合或 `all`（别名 `1,2,3`） |
| `--merge` | 关 | 只跑跨基因组汇总 |
| `--no-merge` | 关 | 批量完成后不自动汇总 |
| `--force` | 关 | 忽略断点全量重跑 |
| `--fast` | 关 | Stage1 用 diamond `--fast`（约快 6 倍，位点数约 −54%，千种全扫推荐） |
| `--overlap` | 关 | 切块步长减半（25kb 重叠，修跨界截断；重点基因组复扫用） |
| `--skip-s3` | 关 | 跳过 RVDB 层 |
| `--blastn-viroid` | 关 | 附加类病毒 blastn 层（需 blastn 可执行 + `--viroids-db`） |
| `--cleanup` | 关 | 删除已被下游产物取代的中间文件（chunks/raw/hits/cand3） |
| `--samtools` | 取 yaml `tools.samtools` | 抽序列用 `samtools faidx` 的可执行文件（空 = PATH 里的 samtools） |
| `--cmd-timeout` | 0（不限） | 单条外部命令超时秒数（防 diamond/seqkit 挂死白占并行核） |
| `--config/--profile` | 自动查找 | YAML 配置文件与预设 |

数据库路径三级优先级：**CLI > pipeline_config.yaml（profile）> 环境变量默认**
（`MMPV_VIRUS_DB` / `MMPV_DB_ROOT` / `MMPV_EVE_PV_DB` / `MMPV_EVE_RVDB_DB`）。

### 输出

```
{outdir}/
├── 01_Loci/<NAME>/                    Stage1 产物
│   ├── <NAME>.fna                     归一化后的基因组 FASTA (tar.gz 自动解出)
│   ├── <NAME>.chunks.fna              50kb 滑窗切块 (seqkit sliding -s 50kb)
│   ├── <NAME>.s1_raw.tsv              diamond blastx 原始输出 (块内坐标)
│   ├── <NAME>.s1_hits.bed             坐标还原后的命中 (bed 6列 + 原始行内嵌第4列)
│   ├── <NAME>.loci.bed                300nt 间隔合并后的 EVE 位点
│   ├── <NAME>.loci.fa                 位点序列 (samtools faidx 按区域文件抽取)
│   ├── <NAME>.loci.fa.regions.txt     抽序列用的区域文件 (chr:FROM-TO, 1基闭区间)
│   ├── <NAME>.genome.fai              基因组索引 (固定写输出目录, 不依赖参考目录可写)
│   └── <NAME>.s1_best.tsv             每位点最佳命中 (科级粗分类 ref_family)
├── 02_Verdict/<NAME>/
│   ├── <NAME>.s2_raw.tsv              locus vs 植物/病毒拆分库 blastx 原始输出
│   └── <NAME>.s2_verdict.tsv          判定: locus / verdict / 两侧最佳 id+bitscore
├── 03_RVDB/<NAME>/
│   ├── <NAME>.cand3.bed               送 RVDB 的候选 (viral_supported + undetermined)
│   ├── <NAME>.cand3.fa                候选序列
│   └── <NAME>.s3_rvdb.tsv             RVDB 深度注释 (每位点最佳 hit)
├── 04_Summary/
│   ├── <NAME>_eve_summary.tsv         单基因组汇总 (12列, 断点/完成标记)
│   └── <NAME>_viroid.tsv              [可选] 类病毒层结果
├── logs/                              eve_screen.log (批次) + <NAME>.log (单基因组)
├── kingdom_summary.tsv                merge: 每基因组判定计数
├── kingdom_loci_all.tsv.gz            merge: 全部位点汇总 (13列, 含 genome 列)
└── family_by_genome.tsv               merge: ref_family × genome 矩阵
```

`<NAME>_eve_summary.tsv` 列（与 `kingdom_loci_all.tsv.gz` 只差首列 `genome`）：

| 列 | 含义 |
|---|---|
| `genome` | 基因组名（仅汇总表有） |
| `locus` | 位点标识，如 `chr1:99199-100700` |
| `verdict` | `viral_supported` / `host_like` / `undetermined` |
| `ref_family` | Stage1 最佳命中科级粗分类（轻量关键词映射，非精确分类） |
| `ref_sseqid` / `ref_bitscore` / `ref_qcov` | Stage1 最佳命中详情 |
| `best_viral_sp` / `best_viral_bs` / `best_plant_sp` / `best_plant_bs` | Stage2 两侧最佳 |
| `rvdb_stitle` / `rvdb_bitscore` | Stage3 深度注释 |

---

## 2. eve_genome_scan.py — 单基因组 worker

**用途:** 对单个基因组跑 Stage1→2→3 + summary。由编排器并行调用，也可独立运行
（单基因组重跑/续跑，不必启动整个批次调用池）：

```bash
python eve_genome_scan.py -g genome.fa -n NAME -o OUTDIR -t 40 \
    --ref-db virus_ref.pep.dmnd --pv-db plant_virus.dmnd \
    --id2div id2div.tsv --rvdb-db RVDB.dmnd [--fast] [--blastn-viroid] \
    [--stage 2,3] [--cmd-timeout 7200]
```

断点语义与编排器一致：阶段产物文件存在即完成（空文件 = 完成且零结果），
`--force` 删除旧产物重跑。`--stage`（别名 `--stages`）与编排器同名同义；
数据库走单跑 CLI 指定（批量跑由编排器按 `pipeline_config.yaml` 注入）。
输入 `.tar.gz`/`.gz` 而 `--stage` 不含 `discover` 时，自动改用 Stage1 落盘的
`01_Loci/<NAME>/<NAME>.fna`（须先跑过 `discover` 解包）。

## 3. eve_scan_core.py — 纯逻辑核心

**用途:** 全部可单测的纯函数（不调外部工具）：

| 函数 | 职责 |
|------|------|
| `restore_coordinates` | `_sliding:START-END` 块内坐标 → 基因组坐标，写 hits BED（原始行内嵌第 4 列，Stage2 复用免重解析） |
| `merge_loci` | 300nt 间隔内命中合并为位点 |
| `assign_best_hits` | 每位点取 bitscore 最高命中（双指针线性扫描） |
| `verdict_loci` | 两侧 bitscore ≥50 比较 → viral_supported / host_like / undetermined |
| `family_of` | stitle 关键词 → 科级粗分类（已知局限见源码注释：只匹配属级词根） |
| `summarize_genome` / `merge_all` | 单基因组 12 列汇总 / 跨基因组三表汇总 |
| `sliding_cmd` | `seqkit sliding` 命令构造（显式 `--greedy` 布尔开关 + 基因组位置参数） |
| `bed_to_regions` / `extract_cmd` | BED(0基半开) → samtools 区域文件(1基闭区间) / `samtools faidx --fai-idx --region-file --output` 命令构造 |
| `check_extracted` | 抽序列产物核对（区域 ↔ 记录一一对应、长度相等）；samtools 遇越界区域退出码仍是 0，只往 stderr 写一行，不查就是静默空序列/截断序列 |
| `samtools_too_old` / `check_tools` | 开跑前体检：samtools < 1.11 缺 `--fai-idx`，拦在批处理之前，而不是每个基因组跑到 Stage1 才炸 |
| `stage1/2/3_discover…` | 阶段编排（`seqkit sliding` + `diamond blastx` + `samtools faidx` 抽取） |

## 4. tests/ — 单元测试

`cd endogenous_virus_pipeline && python -m unittest discover -s tests`

73 个用例，两个文件：

- `test_eve_core.py`（54 个）覆盖坐标还原（含负链空位）、位点合并边界（≤300 合并 /
  >300 分开）、两侧判定三分支、科级映射关键词与已知缺口、summary/merge 产物格式、
  阶段别名解析、跨表 join 键 `:`/`_` 兼容（`TestLocusKeyJoin`）、seqkit 命令构造
  （`TestSlidingCmd`：`--greedy` 不吃值、基因组必须位置参数）、工具体检
  （`TestCheckTools`：samtools 必查、blastn 仅类病毒层需要、缺 diamond/seqkit
  要报出来）、samtools 版本探针（`TestSamtoolsVersionProbe`：1.21 放行、1.9/1.10
  拦下、htslib 版本行不误判、探测失败不拦）、BED→samtools 区域文件换算
  （`TestBedToRegions`：0 基半开 → 1 基闭区间、起止倒挂/列数不足/注释行跳过、
  空 bed 得 0 行）、抽取命令（`TestExtractCmd`：索引走 `--fai-idx` 且不是 `--fai`、
  区域走文件不占命令行）、**抽序列产物核对**（`TestCheckExtracted`：真 samtools
  对越界区域 exit 仍是 0，这里兜住空序列/截断/丢记录，头写成 `_` 也认——
  和下游 join 同一把 `locus_key` 尺子）、
  外部命令超时杀进程（`TestRunTimeout`）、batch.tsv 解析与重名中止
  （`TestBatchJobs`）、跳过 Stage1 时 fasta 推断（`TestResolveGenomeFa`）。
- `test_eve_distinguish.py`（19 个）覆盖子模块的 verdict→action 映射（12 种 verdict
  全覆盖、表外新值按 REVIEW 兜底且返回非零、10 个搬运列逐列对齐）与 `run_all.sh` 的
  装配：参数校验（缺 `-A` 必须报错、`--no-host` 才允许不带、空 `-A` 目录同样报错）、
  模块面板复用、五个 python 阶段的文件交接、`-D` 按发现管道目录约定定位输入与默认
  OUT（`-o` 显式覆盖）、`--no-host` 与缺 `HOST_MAP` 时 `locus_architecture.tsv` 仍是
  完整 16 列表头、缺 `baits.dmnd` 不中断、假 diamond 以 rc=137 崩掉时必须留下
  `FAILED` 行而不是静默退出且不落半成品表。用假 `diamond`（只落空文件）跑端到端，
  不需要真实外部工具。

本地无 diamond/seqkit/samtools，单测只覆盖纯逻辑；另用**离线 mock E2E**
（进程内打桩 `run()`，mock 对选项名做白名单校验、未知选项直接报错，复刻真
samtools 的 unrecognized option）验证过全流程：滑窗 → blastx → 合并 →
抽序列 → 双侧判定 → RVDB → summary → merge / 跨基因组汇总 / 断点重跑 skip /
单 genome worker。两种抽序列行为（原样写头、`:→_` 改写）结论一致。

**已在服务器核实**（246 核 Linux，samtools 1.21 / seqkit v2.11.0，取一份 32M
真实植物组副本，只读 `/tmp` 不碰正式数据）：BED→区域→`samtools faidx` 抽出的
10 条边界区域（染色体首/末、单碱基、奇偶起止、7kb 长区域、多染色体）与纯
Python 0 基切片、`seqkit subseq --bed` 两条独立参照逐条零差异；`--fai-idx` 把
索引写到指定路径，参考基因组所在目录设成只读也不报错、不被污染；输出头就是
区域串，每条记录一行（实测 1.21 的 `--output` 不折行）。另外实测确认 samtools
对越界区域退出码仍是 0（只往 stderr 写 `[faidx] Zero length/Truncated
sequence`），故 `check_extracted()` 做产物核对——这条不查就是静默空序列。
仍未覆盖：真实 blastx 命中数与原脚本逐位点相同（需真库+真数据全批跑），以及
万级位点规模下的耗时。

---

## 上下游衔接

- **← public_metadata_pipeline:** `--stage hostref` 产出宿主参考基因组
  （`host_reference/genome/all.genome.uniq.fasta`），作为本管线 `-g`/`-B` 输入
- **← 用户提供:** 自备组装基因组 FASTA（`.fa/.fasta/.fna`，可 `.gz`/`.tar.gz`）
- **→ virome_discovery_pipeline:** `04_Summary` + `kingdom_*` 全基因组 EVE 清单，
  可与本管线子模块 `eve_distinguish/` 的 contig 级判定交叉验证；viral_supported 位点可入
  rescue/CheckV 复评通道。发现管道跑完再单独调 `eve_distinguish/run_all.sh`
  （见上文"为什么不做成 stage"），拿到 `dna_vs_eve_filter.tsv` 后按 `action` 列过滤
- **→ virome_analysis_pipeline:** `family_by_genome.tsv` 供科级组成统计与可视化

## 移植说明 (来源 hi-fever eve_kingdom.py)

- 硬编码路径全部参数化（`EveConfig` + CLI + `pipeline_config.yaml` 三级覆盖）
- 抽序列改为 `samtools faidx` + 区域文件（原脚本强依赖 `~/biosoft/faextract`）。
  四个为此做的决定：① BED（0 基半开）→ 区域文件（1 基闭区间 `chr:FROM-TO`）
  先由 `bed_to_regions()` 换算，坐标错一位是全批位点序列静默错位，必须有单测盯；
  ② 区域走 `--region-file` 而不是逐个命令行参数 —— 大基因组位点数以万计，
  命令行长度会撞 `ARG_MAX`；③ `--fai-idx` 指向输出目录里的 `<NAME>.genome.fai` ——
  samtools 默认把 `.fai` 写在参考旁边，共享只读库目录（如 `/db`）下会直接失败。
  需 samtools ≥ 1.11：`--region-file`/`--output` 从 1.10 起，但 `--fai-idx` 到 1.11
  才有；`check_tools()` 会跑一次 `samtools --version` 把版本不够拦在开跑前。
  选项名必须写全 `--fai-idx` 而不是 `--fai` —— 后者在某些版本上能被无歧义前缀
  匹配到、看起来能跑，但那是未文档化的解析行为，不能依赖；
  ④ `_extract()` 抽完必须过 `check_extracted()` —— 实测 samtools 1.21 遇到越界
  区域退出码仍是 0，只在 stderr 写一行 `[faidx] Zero length/Truncated sequence`，
  输出里留下空序列或少一截的序列；位点坐标是还原出来的，还原一错正好是这个
  症状，不核对就是脏序列静默流到 Stage2/3
- Stage2 原始 blastx 输出并判解析复用（原脚本两轮 blastx 重打）
- 修了原关键词表 1 个实际 bug：`caulimo` 匹配不到 Cauliflower mosaic virus
  （`cauli-f-lower` 无 `m`），CaMV 标题会落 `Other_viral` → 补 `cauliflower` 词条
- 数值列统一整型化输出（`1.2e+02` → `120`，便于下游数值解析）
- `seqkit sliding` 改为显式 `--greedy` + 基因组位置参数：原写法把基因组路径塞给
  布尔开关 `-g` 当值，靠 pflag"漏值成位置参数"碰巧跑通 —— greedy 名义上由参数
  控制实际恒开、意图全丢，且路径以 `-` 开头即失败（`sliding_cmd()` 锁定）
- 跨表 join 键归一化 `locus_key()`：`loci.bed` 第 4 列写作 `chr1:50-699`，而抽序列
  工具写 FASTA 头、diamond 回显 qseqid 时是否把 `:` 改写成 `_` 取决于具体工具
  （samtools faidx 原样写区域名，历史版本不保证）。Stage2 `s2_verdict`、
  Stage3 候选集、`summarize_genome` 三处原来都靠"两边完全一致"隐式成立，一旦某个
  工具改写分隔符，Stage2 判定会被**静默丢弃**（summary 全变 `undetermined`、
  Stage3 候选为 0、RVDB 层静默清空，结论错误但不报错）。现在两侧都过
  `locus_key()`（仅 join 用，输出仍是 `chr1:50-699`）
- 新增断点续传、并行、`.tar.gz` 基因组包解包、清理与汇总表
- 通用性加固（246 核无人值守跑批场景）：batch.tsv `NAME` 重名开跑前中止（重名会
  共用 `01_Loci/<NAME>/` 互相覆盖）；`--blastn-viroid` 前置体检查 `blastn` 本身；
  跳过 Stage1 时压缩输入自动指向 Stage1 解压产物（`resolve_genome_fa()`，samtools
  读不了 `.tar.gz`）；数据文件读写统一 utf-8（服务器 C locale / 开发机 GBK 下不炸）；
  外部命令可设超时（`--cmd-timeout`）防挂死占核
