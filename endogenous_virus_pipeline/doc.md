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

可选的下游一层（同一输出根目录，产物落 `05_Elements/`）：

```
eve_element_db.py (元件级后处理层)
  ├─ collect → cluster → regions → consensus   viral_supported 位点 → 跨基因组元件
  │                                            → 元件区 → consensus 文库
  └─ rm / cdd / phylo                          → 拷贝数定量 / 结构域注释 / 系统发育输入
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

输入约定（`-D` 时自动找，也可 `-q`/`-e` 显式给）：

| 输入 | 自动定位路径 | 说明 |
|---|---|---|
| 候选 fasta | `<DISCOVERY>/eve_screen/evescan_query.fasta` | EVE 筛查产出的查询集（实测 1,027 条） |
| 候选 fasta（备选） | `<DISCOVERY>/10_Reports/eve_candidates.fasta` | 历史命名，只有早期目录里才有 —— **不是发现管道的产物**，别当成主路径 |
| 证据表 | `<DISCOVERY>/10_Reports/rescue_evidence_scored.tsv` | `report` 阶段拷进 `10_Reports` 的副本 |
| 证据表（备选） | `<DISCOVERY>/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv` | 只跑到 `analysis_verify` 的树里用这份原位文件 |

方法与踩坑记录见 [`eve_distinguish/README.md`](eve_distinguish/README.md)。

**v5.1 改了四个会改下游数字的 bug**（服务器 1,027 条真实候选复跑确认）：无寄主命中的
候选不再被并成伪位点、s2 组件并集不再按字符撕裂、寄主否决区分自身物种与跨物种、
`locus_architecture.tsv` 补回表头行。

**v6 / v6.1 又改了四处**（同一条 1,027 条输入，同一套 `panel_hits.tsv` 复跑，逐条对比）：

| 症状 | 根因 | v6.1 实测 |
|---|---|---|
| 一条 104 aa 的比对被打成 4 个组件（MP+CP+RT+RH），于是 `best_locus_ncomp=4` 直接进 `EVE_STRONG_provirus` | 组件标签是**面板条目名**的事实，比对本身只覆盖蛋白的一小段。旧版把条目名直接当成本条证据 | `credit_components()` 要求每条组件有独立 query 区间 + 至少 `COMP_AA_MIN=100` aa 的长度支撑；新增 `credit_hsps` 列记录真正贡献了证据的 HSP 条数。`best_locus_ncomp>=4` 48→18，`provirus_scale` 12→11 |
| `multi_locus_arch` 20 条 TRUE，但逐条查都是短比对堆出来的 | 同上：70–90 aa 的比对也能算一个组件，凑够 3 个就成"多基因座" | 20→**0**（见 README "multi_locus_arch 现状"） |
| `locus_full` 与 `provirus_scale` 取或：6 条候选没有任何寄主位点，却凭 s2 单基因座进了 `EVE_STRONG_provirus` | s2b 是"这个基因座在寄主基因组里成立吗"的唯一判据，s2 的 contig 级规模只是旁证 | 纵向集成判定只认 `locus_full`；`locus_ncomp` 改取 s2b 的值 |
| 53 条 `host_contamination_likely`（会删数据）里，逐条查上游全带病毒科注释 | 上游给了病毒科 = 病毒侧有证据，与"这就是寄主基因"直接矛盾，旧版不理会照删 | 寄主否决拆成三个分支，其中只有 `host_contamination_likely` 会 `REMOVE`；矛盾的一律转 `host_conflict_review` |

verdict / action 的整体迁移（同一份输入，HEAD=仓库已提交版本）：

| | review | EVE_suspect | host_conflict_review | host_homology_cross_species | host_contamination_likely | EVE_STRONG | virus_candidate | → REMOVE | → MOVE_EVE | → KEEP_virus | → REVIEW |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HEAD | 530 | 108 | 0 | 0 | 53 | 9 | 4 | 53 | 179 | 4 | 791 |
| v6 | 530 | 113 | 31 | 18 | 4 | 4 | 2 | 4 | 174 | 2 | 847 |
| v6.1 | 530 | 113 | 33 | 18 | **2** | 4 | 2 | **2** | 174 | 2 | 849 |

v6.1 与 v6 的差别只有一条口径：上游"病毒侧有证据"原本只数蛋白通道（`aa_pident>=95%`
且有 `aa_species`）和病毒科注释，漏了 blastn 通道。补上"`nt_pident>0` 且 `nt_species`
非空"这一条后，2 条只带核苷酸命中的候选从"判寄主序列"退回矛盾复核（71.1% 一致到
`PP728250.1`、72.1% 一致到 `gb|HQ633072.1|`）。**注意 `host_contamination_likely`
现在只剩 2 条**，且两条的 reason 都自带"面板为 Cauli 专用、零命中是面板不适用而非
确认无病毒结构基因"的注脚——下游不能把这一列当成"已确证的寄主基因"直接清。

三个新 verdict 名与上下游对旧名的引用有关，同步改动：

- `host_homology_cross_species`：样本没有自身物种映射，命中的是别的物种。同源强度再高
  也只是跨物种旁证，**不能**据此判寄主序列。
- `host_conflict_review`：自身物种命中 + 上游任一通道有病毒证据 → 两通道矛盾，contig 级
  不仲裁，转人工。
- `host_contamination_likely`：自身物种命中 + 三条通道（病毒科 / 蛋白 / 核苷酸）全无记录，
  是唯一允许删数据的分支。

消费 `eve_distinguish` 产物时需要知道的两个新口径：

- `locus_architecture.tsv` 多出 `host_scope` / `sample_flag` / `n_host_scaffolds` /
  `locus_key` 四列，`locus_arch` 新增取值 `no_host_locus`（该候选没有任何寄主命中，
  架构未判定）。**没有寄主位点 ≠ 不是 EVE**，只是当前取不到基因座级证据；
  s3 对这批最多给到 `EVE_suspect`。
- 自身物种覆盖很有限（实测 1,027 条里只有 73 条落在 `host_sample_map.tsv` 的 23 个
  物种内），扩样本映射表比调阈值更能提升判别力。

---

## 1. eve_screen.py — 主编排器

**用途:** 三阶段 EVE 筛查调度（`discover → verdict → annotate` + `merge`）。
支持单基因组 (`-g`)、批量 TSV (`-B`，路径可直接填 `.tar.gz`/`.gz` 基因组包)；
基因组级 ProcessPoolExecutor 并行（总核数 = `-J × -t`）；断点续传按**请求的阶段**判定：
每个请求阶段的产物文件都在才算完成，已完成基因组自动跳过；
单基因组失败不影响整批（结果以 ok/error 汇总返回，worker 被 OOM killer 硬杀也只记一行
error 而不中断整批）。

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
| `--cleanup` | 关 | 删除已被下游产物取代的中间文件（chunks/raw/hits/cand3 六个文件，各自需护栏产物在才删；终表与 loci.bed/loci.fa/s1_best.tsv 永不删） |
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
│   ├── <NAME>.s1_hits.bed             坐标还原后的命中 (BED3 + 第4列内嵌原始行, 共4列)
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
│   ├── <NAME>_eve_summary.tsv         单基因组汇总 (12列, 不含 genome) + 完成标记之一
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

`cd endogenous_virus_pipeline && python -m pytest tests/ -q`

132 个用例，三个文件：

- `test_eve_core.py`（68 个）覆盖坐标还原（含负链空位）、位点合并边界（≤300 合并 /
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
  完整 20 列表头、缺 `baits.dmnd` 不中断、假 diamond 以 rc=137 崩掉时必须留下
  `FAILED` 行而不是静默退出且不落半成品表。用假 `diamond`（只落空文件）跑端到端，
  不需要真实外部工具。
- `test_eve_distinguish_logic.py`（45 个）覆盖子模块四个判定阶段本身：s1 读框换算与
  同批基线（另两框是哪两个读框、空表/全停止边界、`stop_enrichment` 越大越脏）、
  s2 位点合并（`h["comps"]` 必须按整词并集——历史上这里按字符迭代，把 `AP+RT`
  撕成 `A`+`P`，产线上 334 条 `best_locus_comps` 变成 `RT+A+P` 这种垃圾值）与
  组件记功（短比对撑不起多组件、重叠 HSP 不得堆组件、单条长比对按自身长度封顶、
  `loci_detail` 按 CANON 序输出且不含 set repr）、
  s2b 位点归组（无寄主命中的候选不得跨 contig 并成伪位点、`locus_key` 与
  `no_host_locus` 不变量、确定性行序、两次运行字节相同）、s3 二维表与寄主否决
  （`provirus` 只由 `locus_full` 立、跨物种命中与自身物种命中分开记、三支否决按
  证据通道分流）。表头一律从被测模块 import，不两头各抄一份。

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

## 5. eve_element_db.py — 元件级后处理层（`05_Elements/`）

**用途:** 前四个阶段收敛的是**位点**（viral_supported / host_like / undetermined），
本层把位点进一步收敛成**元件**：跨基因组同源聚类 → 元件区归并 → consensus 文库 →
拷贝数 / 结构域 / 系统发育输入。移植自 hi-fever 的 `eve_kingdom_post.py`
（千种植物 EVE 下游分析定版脚本），算法参数逐项保留，只改输入发现与工具/断点约定。

```bash
# 只跑不需要额外参考库的四个阶段（默认）
python eve_element_db.py -o eve_results --stages collect,cluster,regions,consensus

# 全量（rm 要 -B 与 RepeatMasker；cdd 要 mmseqs profile 库；phylo 出命令脚本）
python eve_element_db.py -o eve_results --stages all -t 40 -J 8 -B batch.tsv \
    --cdd-db ~/database/cdd/cdd-db/cdd_db --rt-ref RT_refs.fa --run-phylo
```

### 阶段

| 阶段 | 输入 | 输出 | 额外依赖 |
|------|------|------|----------|
| collect | `04_Summary/<NAME>_eve_summary.tsv` + `01_Loci/<NAME>/<NAME>.loci.fa` | `all_candidates.fa`、`candidates_meta.tsv` | 无 |
| cluster | `all_candidates.fa`、`candidates_meta.tsv` | `element_table.tsv`、`cluster_members.tsv`、`_members_C%06d.txt` | mmseqs |
| regions | `04_Summary/*` | `element_regions.tsv` | 无 |
| consensus | `element_table.tsv` + `_members_*` | `eve_consensus_library.fa` | mafft（仅 `--mafft-consensus` 时） |
| rm | `eve_consensus_library.fa` + `-B` 的基因组清单 | `rm/<NAME>/<NAME>.out`、`family_quant.tsv` | RepeatMasker（`--rm-bin` 或 `RM_BIN`）、`-B` |
| cdd | `eve_consensus_library.fa` | `consensus_domains.tsv` | mmseqs + `--cdd-db`（profile 库前缀） |
| phylo | `element_table.tsv` + `all_candidates.fa` | `phylo/<家族>_reps.fa`、`phylo_commands.sh` | `--run-phylo` 时还要 `--rt-ref` 与 mafft/trimal/iqtree |

`--stages` 收 `collect/cluster/regions/consensus/rm/cdd/phylo` 的逗号组合（别名
`1..7`）或 `all`；`--rm`、`--cdd-db`、`--run-phylo` 三个开关会把对应阶段加进本次运行
（源脚本就是这么控制的）。**`all` 会连带 rm/cdd/phylo**，需要相应参数，缺了在开跑前
一次性报出来 —— 源脚本那三处是静默跳过（退出码 0、不产出任何表），拿到空表才发现。

### 输入为什么取 `04_Summary` 而不是 `02_Verdict` / `03_RVDB`

summary 是唯一把三张表 join 好的单表：Stage1 的 `ref_family`/`ref_sseqid`/
`ref_bitscore`、Stage2 的 `verdict`、Stage3 的 `rvdb_stitle`/`rvdb_bitscore` 都在。
源脚本的 `evidence_bs` 取 `max(ref_bitscore, rvdb_bitscore)` 两条证据通道的强者，
`ref_family` 又要喂给聚类的家族投票与 phylo 分组；只读 `02_Verdict` 会丢
`ref_family`/`rvdb_bitscore`，只读 `03_RVDB` 会丢 `verdict`（且 03 只收
viral_supported + undetermined 两类候选，host_like 位点根本不在里面）。

locus 与 `loci.fa` 头两侧都过 `locus_key()` 再 join（与 Stage2/3 同一个坑：
抽序列工具把区域串的 `:` 改写成 `_` 时，join 静默落空、位点序列全丢而退出码仍是 0）。
默认口径：`verdict = viral_supported` 且 `evidence_bs >= 50`（`--verdicts` /
`--min-bs` 可调）；**regions 阶段不设证据阈值**（只看 verdict），因为它是计数校正用的
元数据表，阈值一变元件区就跟着变、拷贝数校正也就不可比。

### 输出

```
{outdir}/05_Elements/                  （文件名沿用源脚本 post/ 下的名字）
├── all_candidates.fa                  位点序列, 头 = <基因组名>__<locus>
├── candidates_meta.tsv                genome/locus/family/verdict/ref_sseqid/evidence_bs
├── cluster_members.tsv                mmseqs linclust 原始簇表 (rep<TAB>member)
├── _clu*_seq.fasta, _clu_tmp/         mmseqs 中间产物 (可删; 大基因组下 _clu_tmp 不小)
├── element_table.tsv                  元件主表 (cluster_id/rep/n_members/n_genomes/
│                                      family/chimera_flag/best_member/best_bs)
├── _members_C000000.txt …             每簇成员清单 (consensus/phylo 复用)
├── element_regions.tsv                元件区 (genome/region_id/contig/start/end/
│                                      n_families/families), ±5kb 同家族归并
├── eve_consensus_library.fa           consensus 文库, 头 = <cluster_id>_<家族>
├── rm/<NAME>/<NAME>.out               RepeatMasker 逐基因组产物 (断点粒度到基因组)
├── family_quant.tsv                   每家族 x 每基因组 拷贝数 / 覆盖碱基
├── _cdd_mm/                           mmseqs 中间库 (可删, 每次重跑前清掉)
├── consensus_domains.tsv              cdd 命中 (query/target/evalue/bits/q…/tlen,
│                                      convertalis 直出, 无表头)
├── phylo/<家族>_reps.fa               家族代表序列 (每簇最长成员)
└── phylo_commands.sh                  mafft → trimal → iqtree 命令脚本 (bash)
```

表头在模块里是导出常量（`CAND_META_HDR` / `ELEMENT_TABLE_HDR` /
`ELEMENT_REGIONS_HDR` / `FAMILY_QUANT_HDR`），写盘与测试共用同一份字面量；下游
`build_eve_db.py` 认的就是这几个文件名与列序，改名等于改契约。

**断点语义:** 与 `eve_screen.py` / `eve_genome_scan.py` 一致 —— **该阶段声明的产物
文件全部存在即视为完成**（清单见 `STAGE_OUTPUTS`），重发同命令即从断点继续，
`--force` 重跑请求的阶段（rm 阶段连带忽略每基因组的 `.out` 标记）。不用 `.done` 标记：
源脚本那套 `.done` 与产物文件两套状态并存，删了标记就重跑、留着标记就不跑，容易对不上。
每个阶段完成时打一行 `[done] <阶段>: <计数> -> <产物>`，`grep '\[done\]'` 即得阶段清单。

### 与源脚本的差异（都是有理由的偏差，不是重新调参）

- **RepeatMasker `.out` 列下标修正**（源脚本 `eve_kingdom_post.py:387-391` 取错一列）：
  见下面"移植时发现的源脚本 bug"第 1 条。算法参数（`-no_is -norna -pa`）原样保留。
- **`--run-phylo` 改为交给 bash 执行 `phylo_commands.sh`**：源脚本把生成的行按空白
  拆成 argv 直接执行，`set -euo pipefail` 会被当成可执行文件、mafft 那行的 `>` 重定向与
  `cat` 的引号都会变成普通参数 —— 那条路本来就跑不通（第 2 条）。
- **缺工具/参数开跑前报错**而不是静默跳过（第 3 条）；`phylo` 不带 `--rt-ref` 时仍只出
  代表序列与注释行（与源脚本一致，这是正常用法）。
- **确定性**：成员表、簇顺序、家族列表一律排序（并列时按名字定序）。源脚本的簇序与元件
  ID 依赖 mmseqs 输出行序，同一份输入换个版本就可能整体搬家，无法做逐字节回归。
- **IO**：全部显式 `utf-8` + `newline="\n"`（源脚本用默认编码与默认换行：本机 GBK 下
  读非 ASCII 名字会炸，Windows 上写出的表带 CRLF，bash 侧按行读会多出 `\r`）。
- **rm 单个基因组失败不连带整批**：逐个 future 兜住，失败的基因组名一起报出来；
  任何一个失败都不写 `family_quant.tsv`（不拿半份数据当完整结果），下次重发只补失败的。

测试：`tests/test_eve_element_db.py`（70 个用例）三层 —— 纯逻辑（locus 解析 / 元件区
归并 / 候选筛选 / consensus 输入装配 / 比对列投票 / 批量表解析 / `.out` 解析）、
端到端（假 mmseqs / mafft / RepeatMasker 跑通 collect→cluster→regions→consensus→rm，
逐条核对表头常量、产物文件、`[done]` 行与**工具收到的参数**）、失败可见性（工具崩掉留
`[fail]` 且不留半成品表；rm 缺 `-B`、cdd 缺 `--cdd-db`、`--run-phylo` 缺 `--rt-ref`
必须在开跑前报出来）。输入表不是手抄的：由 `eve_scan_core.summarize_genome()` 生成，
上游改列名这里立刻失败。

#### 移植时发现的源脚本 bug（已在本层修正，源脚本未动）

1. `eve_kingdom_post.py:387-391`（RepeatMasker 定量）**取错一列**：`.out` 的
   0 基第 4 列是 query 序列名，第 5/6 列才是 query begin/end；源脚本拿 `c[4]`/`c[5]`
   当坐标并作 `isdigit()` 守卫 —— `c[4]` 对真实数据行永远不是数字，于是**每一行命中
   都被当成表头跳过**，`family_quant.tsv` 只剩表头而退出码是 0；反过来若 contig 名恰好
   是纯数字（`'1'`、`'2'`… 很多组装就是这么命名的），守卫放行，又按
   `int(c[5]) - int(c[4]) + 1` 算出"起点 − contig 号 + 1"的垃圾 `covered_bp`。
   本层按官方列序取 5/6，并加了"第 8 列必须是链"的守卫。
2. `eve_kingdom_post.py:443-444,473-477`（`--run-phylo`）**执行方式跑不通**：
   `set -euo pipefail` 不是可执行文件（bash 内建），`mafft … > out.fa` 的 `>` 与
   `cat` 行的引号都会被 `line.split()` 当成普通参数。本层直接 `bash phylo_commands.sh`。
3. `eve_kingdom_post.py:526,529,473` **三处静默跳过**：`rm` 在 stages 里但没有 `-g`、
   `cdd` 没有 `--cdd-db`、`--run-phylo` 没有 `--rt-ref` 时那一阶段什么都不做而退出码 0。
   本层改成开跑前体检并指名报缺。
4. `eve_kingdom_post.py:173-175`：mmseqs 退出码 0 但没产出 `_clu_cluster.tsv` 时，
   `cluster_members.tsv` 不会被创建，下一行 `opf()` 抛的是指向中间文件的
   `FileNotFoundError`（看不出是 mmseqs 没产出）。本层显式检查并报出期望路径。
5. `eve_kingdom_post.py:322-325,451-454`：成员序列已不在 `all_candidates.fa` 里的簇被
   **静默丢弃**（文库比元件表少几行没人知道）；`--mafft-consensus` 时更会拿空集合去跑
   mafft 直接失败。本层改为按"可用序列数"判断并打 WARN 计数。
6. `eve_kingdom_post.py:49-54`：`opf()` 的 `".." in rp.parts` 是死代码（`resolve()`
   已经把 `..` 消解），与本仓库对 `safe_dir()` 的判断一致 —— 无害，但别指望它拦路径。

另有一个**值对不上**的问题在下游消费者 `build_eve_db.py`（不在本层，未改动）：
`GT` 表（`build_eve_db.py:27-52`）的关键词只匹配属/科级词根，而上游 `family_of()`
会产出 `NCLDV`、`Picornavirales`、`Nodaviridae`、`Furoviridae`、`Mycovirus` 这五个
**一个关键词都匹配不到**的值（`"nclcv"` 看着就是 `"ncldv"` 的笔误，而上游写的是大写
缩写 `NCLDV`）→ 它们的 `genome_type` 落 `Unclassified`、`tier` 落 `unclassified`，
于是**被排除在 `core_viral_elements.tsv` 之外**。实测 31 个可能取值里有 5 个如此。
与已修的那两处（列索引/字段含义混用）同属"写死的字面量与真实值对不上"。

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
- **→ 元件级后处理层（本管线内，可选）:** `05_Elements/` 的 `element_table.tsv` /
  `eve_consensus_library.fa` / `family_quant.tsv` 是跨基因组**元件**目录，供科级组成、
  拷贝数与系统发育分析（见上文"5. eve_element_db.py"）

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

## 判据校准：四个工具实测（v6.2，`eve_distinguish/verify/`）

此前所有轮次证明的都是"代码按定义执行"。判据本身站不站得住，由 verify/ 下四个工具量。
工具与完整数据见 `eve_distinguish/verify/README.md`。四条结论：

1. **对随机序列零假阳性**：1,027 条候选整条单核苷酸洗牌（保组成）后重跑全流程，
   diamond 0 命中、s2 组件 0 条、s3 实质判定 0 条，全部落 `review`。
2. **火力来源不均**：只洗同源区时 `coding_intact` 富集 11.1×（真信号），而
   `distributed_decay` 只有 **1.21×** —— 化石退化里约 **83%** 与随机组成无法区分。
   所以依赖退化通道的 `EVE_suspect`(113) / `ancient_EVE`(15) 证据强度明显低于
   `EVE_STRONG_provirus`(4) / `virus_candidate`(2)，对外表述不可混为一谈。
3. **12 个阈值里 4 个对结论零影响**：`PROVIRUS_SPAN` / `PROVIRUS_NCOMP` /
   `MIN_HEAD_STOPS` / `AA_VIRAL_PID` 扫到底变更 0 条；真正决定结论的只有
   `MIN_PID`(11.5%) 与 `COMP_AA_MIN`(8.2%)。**"前病毒规模"已不再是判据**（与
   `multi_locus_arch` 同性质，都是信息列）。
4. **组件判定独立复核未定论**：hmmscan+Pfam 对 RT 的对照灵敏度 100%（可用）、
   AP 20% / RH 0%（不可用）。可解读的只有 RT：一致率 26.6%（下界），另有 98 条候选
   oracle 看到 RT 而 s2 未记 —— 这是提高召回率的具体线索。

### 复用派生表的坑

`EVE_VERIFY_DIR` 与任何复核都必须用**当前代码重算**的派生表。实测复用旧
`s2_domains.tsv`（早于 `credit_hsps` 那次修改）与旧 `locus_architecture.tsv`，会得到
`EVE_STRONG_provirus` 11 条 / `virus_candidate` 5 条；而从原始命中重算是 **4 / 2**，
与文档记载的 v6.1 数字逐条一致。差异方向正好是"结论变了"，不看就会以为代码又坏了。
`s1`/`s2`/`s2b` 都是纯计算（`s2b` 复用缓存 `locus_blastn.tsv`），只有 blastx 与 host
blastn 贵 —— 那两份原始命中不用重跑，派生表则要重算。

## 第二轮修复：三处静默产出错文件的问题（实测复现）

上面那些改完之后又对着代码审查了一轮，修掉三处**会静默产出错文件**的问题，外加一处
文档级错误。三处都在 `eve_scan_core.py` / `eve_genome_scan.py`（筛查侧），不在
`eve_distinguish/`（判别侧）。

| # | 症状 | 根因 | 修后 |
|---|---|---|---|
| 1 | `kingdom_loci_all.tsv.gz` 表头 13 个名字、数据行 14 个字段 —— 用 pandas / `csv.DictReader` 按列名读时每列右移一位（`locus` 列显示基因组名、`verdict` 列显示位点） | `summarize_genome` 写 summary 时首列已经写了 `genome`（13 列），`merge_all` 汇总时又拼一遍 | summary 回到契约里的 12 列（不含 genome），汇总表 13 列严格对齐表头；汇总时逐行校验列数，不齐就补空/截尾并告警，且兼容盘上已有的 13 列旧 summary |
| 2 | `--stage 1` 跑过之后，补跑 `--stage 2,3` **一句"跳过"就返回 0、什么都不跑**；那份只有 `undetermined` 的 summary 被 `merge` 当成真结果统计，`viral_supported` 恒为 0 | 完成标记只看 `04_Summary/<NAME>_eve_summary.tsv` 是否存在，而这份 summary 在只跑了部分阶段时**也会写** | 完成判定改为按**请求的阶段**逐个查产物（discover→`loci.bed`、verdict→`s2_verdict.tsv`、annotate→`s3_rvdb.tsv`、类病毒层→`_viroid.tsv`），`--skip-s3` 时不要 RVDB 产物；补跑缺失阶段即生效 |
| 3 | `--cleanup` 六条规则里**只有两条真的会删**；`s1_raw.tsv` 与 `s1_hits.bed`（剩余体积最大的两个中转文件）从未被删掉 | 键是"逻辑名"而落盘名不同（`loci.bed` 的主体是 `loci` 不是 `loci_bed`，`cand3.fa` 的是 `cand3` 不是 `cand3_fa`），且护栏扩展名写死成 `.tsv` 而实际是 `.bed` | victim/guard 两侧都按磁盘真实文件名解析；新增用例逐条对照 stage 函数真正写的文件名，防再次写错 |

另有三处健壮性加固与一处文档修正：

- `verdict_loci` 对不在 `id2div` 表里的 `sseqid` 现在会告警并给出条数 —— 原先这些命中
  被静默忽略，若 id2div 与 `pv_dmnd` 版本不配套或分隔符不对，全部位点会落成
  `undetermined`、Stage3 候选为 0，而批次退出码仍是 0。
- worker 被 OOM killer 硬杀（`BrokenProcessPool`）时逐 future 兜住并记成一行 error，
  不再让剩余的基因组整批消失；文档承诺的"单基因组失败不影响整批"这才真正成立。
- 两个 CLI 的 `WINDOW/EVALUE/MERGE_D/HOST_BS/VIRAL_BS` 默认值改为从 `eve_scan_core`
  import，不再各抄一份字面量（此前三处数值一致但随时可能分叉）。
- `restore_coordinates` 的 docstring 把 seqkit 的块起点写成 0-based，**代码是对的**
  （1-based 闭区间），但照着注释"修"代码会引入 1 个碱基的静默偏移，而唯一的护栏
  `check_extracted` 只比名字和长度、偏移保长度查不出来 —— 已改成把换算逐行写清。

测试：三套共 **132** 个用例（原 118 + 本轮 14 个回归用例），全部可离线跑。