# phylo_pipeline.py 运行指南（v3 · 三入口 + 17 stage）

> **📌 文档定位：本文件是「怎么跑」的权威** —— CLI 参数、config 分区、输出目录布局、
> checkpoint 机制一律以本文件为准。stage 的逐项内部细节见 `STAGE_REFERENCE.md`；
> 数据流/模块归属见 `PIPELINE_FLOW.md`；论文方法学表述见 `METHODS.md`；
> 从哪份文档读起见 `README.md`。
>
> **同步状态**：2026-09-16 按代码逐项重写（旧版 2026-08-27 的整体过期内容已全部替换）。
> 改动前备份：`_docconv_20260916/backup/RUN_GUIDE.md`。
> 机器校验：`_consistency_check_20260916/verify_docs_vs_code.py`（文档 ↔ 代码对拍）。

---

## 一、目录结构

```
virome_phylo_pipeline/
├── phylo_pipeline.py           # 主编排器（stage 流程化，96 个 CLI 参数）
├── config.yaml                 # 默认运行参数（优先级：CLI 显式 > config.yaml > 代码默认）
├── datasets.yaml               # 数据集注册表（--virus 模式用：每病毒路径/配色/基因坐标）
├── capheine_pipeline.py        # 正选择分析（stage `capheine` 的实现）
├── gene_partition_dating.py    # 分基因 BEAST 定年（stage `gene_dating` 的实现）
├── run_phylogeo.py             # BEAST 统一入口（XML 生成 → 提交 → 状态登记）
├── popgen_analysis.py          # 群体遗传（stage `popgen` 的实现）
├── recombination_analysis.py   # RDP5 重组检测（stage `rdp5` 的实现）
├── selection_suite.py          # 正选择下游（独立 CLI，未进总管线）
├── gbk_extractor.py            # GBK → CDS（分基因流的上游）
├── utils/                      # 依赖模块（51 个；beast1_bridge / virphy_bridge / metadata_channels …）
└── *.md                        # 文档（权威关系见 README.md）
```

---

## 二、三种入口模式

| 模式 | 命令 | 适用 |
|---|---|---|
| **批量** | `--extract_dir 06_extraction/ --metadata <Core14.tsv>` | 多病毒全流程（自动扫描提取目录） |
| **单病毒注册名** | `--virus GCVA` | 已在 `datasets.yaml` 注册的病毒，自动定位比对/元数据 |
| **work-dir** | `--work_dir <dir> --alignment ... --tree ... --meta ...` | 单病毒，已有比对+树，只补下游 |

三者走**同一套 stage 调度**。差异只在「输入从哪来」：

- 批量模式：`data/dates.csv` / `data/sample_metadata.csv` / `data/host.csv` 三通道由 `prep` 写出；
- `--virus` 模式：三通道指向 `datasets.yaml` 登记的元数据；
- work_dir 模式：`--meta` 给时间通道，`--host_meta` 给宿主通道（**省略则由 `--meta` 拆出**），地理通道由 `--meta` 的 location 列派生。

---

## 三、Stage 体系

### 3.1 18 个 stage（`--stage` 子步骤名，与 `STAGE_ORDER` 逐字对应）

| # | stage | 依赖（自动补齐） | 输出（`<virus>/` 内） | checkpoint 标志产物 |
|---|---|---|---|---|
| 1 | `prep` | — | 合并序列 + `data/{dates,sample_metadata,host}.csv` | 不参与 checkpoint |
| 2 | `metadata` | prep | `data/metadata/metadata_std.csv`、`metadata_report.tsv`；三通道回填 | `.checkpoints/metadata.done`（仅标记） |
| 3 | `online` | prep, metadata | `data/ncbi_ref/`、`data/ncbi_harvest/` | `data/ncbi_ref` |
| 4 | `clean` | prep | `data/clean/clean.fasta`（+`removed.fasta`/`clean_report.tsv`） | `.checkpoints/clean.done`（仅标记） |
| 5 | `align` | prep, clean | `phylogeny/mafft.aln.fasta` + `phylogeny/saturation_out/` | `phylogeny/mafft.aln.fasta` + **输入指纹** |
| 6 | `splitstree` | align | `phylogeny/splitstree/neighbornet.nexus` | `phylogeny/splitstree/neighbornet.nexus` |
| 7 | `rdp5` | align | `recomb/rdp5/`（事件 CSV + 干净比对 + `validation/`） | `recomb/rdp5/rdp5_dropped.realn.fasta`（`--rdp5_mask` 时 `masked/masked_N.fasta`） |
| 8 | `tree` | align | `phylogeny/iqtree.treefile` + `iqtree.support.json` | `phylogeny/iqtree.treefile` + **输入指纹（含 IQ-TREE 参数）** |
| 9 | `popgen` | align | `popgen/popgen_report.txt`（+ 图/网络/HTML） | `popgen/popgen_report.txt` |
| 10 | `host` | tree | `popgen/host_analysis/tree_host.pdf` | `popgen/host_analysis/tree_host.pdf` |
| 11 | `capheine` | prep | `select/capheine/`（cawlign/drhip/multiqc） | `select/capheine/drhip` |
| 12 | `clock` | tree | `time/treetime_rtt/`、`time/treedater_ltt/`、`time/temporal_signal/` | `time/treedater_ltt/Phylogeny_dated.pdf` |
| 13 | `beast` | tree | `time/beast/beast1.xml` + `.log`/`.trees`/`postprocess/` | `time/beast/beast1.xml` |
| 14 | `gene_dating` | clock, capheine | `time/gene_dating/`（+`time/gene_split/`） | `time/gene_dating/gcva_watch_parse.sh` |
| 15 | `phylogeo` | tree | `geography/phylogeography/`（XML + `merged/mcc.tree` + `visualization/`） | `geography/phylogeography/merged/mcc.tree` |
| 16 | `geo` | tree | `geography/geo_analysis/mantel_test.pdf`（+`extended/`） | `geography/geo_analysis/mantel_test.pdf` |
| 17 | `geo_paths` | phylogeo | `geography/geo_analysis/pathways/`（+`dispersion/`，可选 `visualization/transmission_map.*`） | `geography/geo_analysis/pathways/pathway_qc.json` |
| 18 | `report` | — | `report/phylo_summary.csv` + `report/phylo_report.html` | `report/phylo_summary.csv` |

> ⚠️ **`beast` 恒为 BEAST 1.x**（产物名固定 `beast1.xml`）。`--beast_version` 只影响
> `phylogeo` 走 BEAST1 还是 BEAST2 路径，**不影响 `beast` stage**。

### 3.2 8 个大模块别名（`--stage data,time` 这样用）

| 别名 | 展开 |
|---|---|
| `data` | prep, online, clean |
| `popgen` | popgen, host |
| `recomb` | rdp5 |
| `select` | capheine |
| `phylogeny` | align, splitstree, tree |
| `time` | clock, beast, gene_dating |
| `geography` | phylogeo, geo, geo_paths |
| `report` | report |

> `metadata` **不在任何别名里**（由 `--metadata` 入口参数驱动、在 `prep` 之后自动执行）；
> 需单独重跑时写 `--stage metadata`。
> `rdp5` 只在 `recomb` 组，不在 `phylogeny` 组。

**依赖自动补齐**：`--stage beast` → 实际执行 `prep, metadata, clean, align, tree, beast`。
用 `--dry` 可先预览真实执行清单。

### 3.3 真正关闭一个 stage：`--disable`

`config.yaml` 里各 `enabled: true` 只保证「不被 `--disable` 关掉」，**不会**把 stage 追加进清单。
要关掉（即使 `--stage all` 也会展开它）用：

```bash
python phylo_pipeline.py --stage all --disable beast,rdp5
```

### 3.4 已废弃 stage（保留兼容，不参与自动编排）

`saturation` / `rtt` / `temporal` / `genes` / `concat` —— 能力已分别并入
`align`（饱和）/ `clock`（RTT+LTT+DRT）/ —— 。**单独调用仍可运行**，但任何大模块都不会展开它们。
历史说明见 `STAGE_REFERENCE.md` 附录 A。

---

## 四、CLI 参数（96 个，按 help 分组）

### Config / Input / Output

| 参数 | 说明 |
|---|---|
| `--config` | `config.yaml` 路径（默认：pipeline 根目录） |
| `--extract_dir` | `06_extraction` 目录（批量模式） |
| `--metadata` | `Global_Unified_Metadata_Core14.tsv` 路径（批量的规范源） |
| `--info_dir` / `--variants_dir` | 上游 `public_metadata_pipeline` 的 `info/` / `03_variants` 目录 |
| `--virus` | 单个病毒名（`datasets.yaml` 注册名，如 `GCVA`/`PSTVD`） |
| `--work_dir` | 单个病毒工作目录（已有 alignment+tree 模式） |
| `--alignment` / `--tree` | work_dir 模式下已有的 MAFFT 比对 / Newick 树 |
| `--meta` | work_dir 模式下的 `dates.csv`（时间通道） |
| `--host_meta` | 宿主元数据 CSV `(name,host[,tissue])` —— **独立成宿主通道 `data/host.csv`，只收有 host 的样品**；省略则由 `--meta` 拆出 |
| `--output_dir` / `-o` | 输出根目录（默认 `phylo_results/`） |

### Runtime

| 参数 | 说明 |
|---|---|
| `--stage` | 大模块名或子步骤名，逗号分隔 |
| `--disable` | 真正关闭 stage（逗号分隔，从最终清单剔除） |
| `--threads` / `-t` | 线程数 |
| `--min_samples` | 最小样本数门槛（低于此的病毒跳过） |
| `--dedup` / `--no-dedup` | `(date, location, seq)` 三元组克隆去重（默认开） |
| `--force` | 忽略所有 checkpoint 重跑（产物存在时复用，不删除） |
| `--skip_msa` / `--skip_tree` | 跳过比对 / 跳过建树 |
| `--max_viruses` | 最多处理几个病毒（0 = 不限） |
| `--iqtree_bin` | IQ-TREE 可执行名（默认 `iqtree2`） |
| `--iqtree_boot` | UFBoot 重复数（`-B`）；默认 **1000**，传 **0** 关闭 |
| `--iqtree_alrt` | SH-aLRT 重复数（`-alrt`）；默认 **1000**，传 **0** 关闭 |
| `--iqtree_bnni` | 给 UFBoot 加 `--bnni`（近缘序列多时压低高估，稍慢） |
| `--iqtree_timeout` | IQ-TREE **单次尝试**超时秒数；默认 **3600**，大比对需上调 |
| `--dry` | 只打印执行计划，不运行 |

> ⚠️ `-B` / `-alrt` 在 IQ-TREE 里**硬性要求 `>= 1000`**，传小于 1000 的值会被静默过滤；
> 想关支持值只能写 0。三段降级阶梯与取证 sidecar 见 `STAGE_REFERENCE.md` §8。

### clean（比对前序列清洗）

`--clean_min_length_ratio`(0.9) / `--clean_max_length_ratio`(1.5) / `--clean_max_n`(0.05) /
`--clean_reference` / `--clean_allow_rna` / `--clean_strip_illegal` / `--clean_allow_gap_chars` /
`--clean_metadata` / `--clean_drop_no_date` / `--clean_drop_no_location`

> `--clean_drop_no_date` / `--clean_drop_no_location` 需要 metadata 表：未显式给 `--clean_metadata` 时
> 会回退到 `prep['meta_csv']`；两者都拿不到时 **`clean` 整段降级为「用原序列」并告警**（不会清空整批）。

### metadata（入口治理：时间 + 地理 检查与矫正）

`--metadata_mode {mid,start}` / `--metadata_date_col` / `--metadata_loc_col` / `--metadata_name_col` /
`--metadata_keep_unknown_geo` / `--metadata_no_derive_coords` / `--metadata_no_blank_placeholders` /
`--metadata_drop_placeholder_rows` / `--metadata_drop_no_date` / `--metadata_drop_no_location`

> 默认口径（2026-09-16 老师拍板）：**坐标先用地名推算，推不出才清空字段；占位符只清字段不删行**。
> 删 `Unknown` 层级后只剩纯国家名的记录（如 `China, Unknown, Unknown_AI`）判为 `country_only`
> → **不造坐标**、location 字段清空；要整行去掉用 `--metadata_drop_no_location`。

### 各分析模块参数

| 模块 | 参数 |
|---|---|
| rdp5 | `--rdp5`、`--rdp5_script`、`--rdp5_genes`、`--mask_min_methods`、`--rdp5_mask`、`--rdp5_stage {run,all}` |
| beast | `--beast`、`--beast_bin`、`--beast_chain`、`--beast_burnin`、`--beast_timeout`、`--beast_chains`、`--skip_beast_run` |
| gene_dating | `--gene_dating`、`--genes_dir`、`--gene_chain`、`--gene_chain_big`、`--big_genes` |
| saturation | `--sat_replicates` |
| popgen | `--popgen_group`、`--popgen_min_n`、`--popgen_perm`、`--popgen_exclude_reference`、`--popgen_exclude_ids` |
| capheine | `--capheine_ref`、`--capheine_unaligned`、`--capheine_code`、`--capheine_workers`、`--capheine_cpus_iqtree`、`--capheine_cpus_hyphy`、`--capheine_mpi/--no-capheine_mpi` |
| phylogeo | `--phylogeo`、`--phylogeo_clock {strict,ucln}`、`--phylogeo_prior {constant,skyline,bdsky,auto}`、`--phylogeo_chain`、`--phylogeo_bf`、`--phylogeo_chains`、`--skip_phylogeo_run`、`--geo_subsample`、`--treeannotator_bin`、`--beast_version {1,2}`、`--rrt_randomized_dir`、geo_paths(2026-09-16): `--no_pathways`、`--direct_snp`、`--indirect_snp`、`--censor_years`、`--n_perm`、`--geo_map`、`--geo_gif`、`--geo_gif_frames`、`--geo_tree_map`、`--calibrate {auto,manual}` |
| clock | `--check_temporal`、`--drt_randomizations`、`--bets`、`--bets_steps`、`--bets_chain` |
| report | `--report_enabled`、`--no-report` |

---

## 五、config.yaml 分区 → CLI dest 映射

优先级：**CLI 显式传参 > config.yaml > argparse 代码默认**。
`runtime` 分区是「key 与 dest 同名」的直通区，其余分区必须逐键登记（`CONFIG_SECTION_ARG`）。

| config 分区 | 键 |
|---|---|
| `runtime` | `stage`, `threads`, `min_samples`, `max_viruses`, `force`, `skip_msa`, `skip_tree`, `iqtree_bin`, **`iqtree_boot`, `iqtree_alrt`, `iqtree_bnni`, `iqtree_timeout`**, `output_dir` |
| `rdp5` | `enabled`(=--rdp5), `script`, `genes`, `mask_min_methods`, `stage` |
| `beast` | `enabled`(=--beast), `bin`, `chain`, `burnin`, `timeout_hours`, `skip_run`, `version`, **`chains`** |
| `gene_dating` | `enabled`(=--gene_dating), `genes_dir`, `chain`, `chain_big`, `big_genes` |
| `saturation` | `replicates` |
| `popgen` | `group`, `min_n`, `perm` |
| `capheine` | `ref`, `unaligned`, `code`, `workers`, `cpus_iqtree`, `cpus_hyphy`, `use_mpi` |
| `phylogeo` | `enabled`(=--phylogeo), `clock`, `prior`, `chain`, `bf`, `treeannotator_bin`, **`chains`**, `rrt_randomized_dir` |
| `geo_paths` | `direct_snp`, `indirect_snp`, `censor_years`, `n_perm`, `map`(=--geo_map), `gif`(=--geo_gif), `gif_frames`(=--geo_gif_frames), `tree_map`(=--geo_tree_map), `calibrate`(=--calibrate) |
| `temporal` | `check`(=--check_temporal), `drt_randomizations`, `bets`, `bets_steps`, `bets_chain` |
| `report` | `enabled` |

> ⚠️ **「键没登记就悄悄丢」是这个文件的历史坑**：`phylogeo.chains` 与 `beast.chains` 都曾
> 写在 `config.yaml` 里却**不在映射中** → 被静默忽略，实际一直用 argparse 默认值（3）。
> 两处均已于 2026-09-16 登记修复。
> **防复发**：`config_orphan_keys()` 会在启动时对未登记的键**显式告警**；
> `verify_docs_vs_code.py` 把它做成硬断言。加 config 键时必须同步登记。

---

## 六、输出目录布局（每个病毒一个 work_dir）

```
<output_dir>/<virus>/
├── data/                          # 模块 data
│   ├── metadata/
│   │   ├── metadata_std.csv       # 标准化后 metadata（下游唯一输入）
│   │   └── metadata_report.tsv    # 逐行体检报告（日期精度/地理层级/删除/小数年）
│   ├── dates.csv                  # 时间通道 (name,date)
│   ├── sample_metadata.csv        # 地理通道 (name,date,location)
│   ├── host.csv                   # 宿主通道 (name,host,tissue) —— 只含**有 host** 的样品
│   ├── ncbi_ref/                  # online: efetch 单参考 (.gbk + .ref.fasta)
│   ├── ncbi_harvest/              # online: SeqHarvester 公共序列 + metadata.csv + 分组图
│   └── clean/                     # clean: clean.fasta + removed.fasta + clean_report.tsv
├── phylogeny/                     # 模块 phylogeny
│   ├── mafft.aln.fasta            # align
│   ├── saturation_out/            # align 内联的替换饱和分析 (Iss)
│   ├── splitstree/neighbornet.nexus
│   ├── qc/                        # 比对后 QC 图 (align_qc)
│   ├── iqtree.treefile            # tree
│   ├── iqtree.iqtree / .log / .ckp.gz
│   ├── iqtree.support.json        # tree: 支持值元数据 sidecar（取证用）
│   └── gene_analysis/             # (已废弃 stage `genes` 的历史产物)
├── recomb/rdp5/                   # rdp5
│   ├── <virus>.csv                #   7 方法事件检测结果
│   ├── rdp5_dropped.realn.fasta   #   默认策略: 删重组子 + 重新 MAFFT
│   ├── masked/masked_N.fasta      #   --rdp5_mask 策略: 重组区置 N
│   └── validation/                #   事件解析 / 拓扑切换 / SimPlot
├── popgen/                        # 模块 popgen
│   ├── popgen_report.txt          # popgen
│   └── host_analysis/tree_host.pdf# host
├── select/capheine/               # capheine (cawlign/ drhip/ multiqc)
├── time/                          # 模块 time
│   ├── treetime_rtt/              # clock: RootToTip_regression.pdf + timetree
│   ├── treedater_ltt/             # clock: Phylogeny_dated.pdf + LTT.pdf
│   ├── temporal_signal/           # clock: drt_results.pdf
│   ├── beast/                     # beast: beast1.xml + beast1.log + beast1.trees + postprocess/ + beauti/
│   ├── gene_split/                # gene_dating: 切出的分基因比对
│   └── gene_dating/               # gene_dating: 每基因 XML + watchdog 脚本
├── geography/                     # 模块 geography
│   ├── phylogeography/            # phylogeo: phylogeo_beast1.xml + merged/mcc.tree
│   ├── geo_analysis/              # geo: mantel_test.pdf (+ extended/)
│   ├── prior_comparison/          # phylogeo: 模型/先验比较
│   ├── tempmig/                   # TempMig 时序迁移追踪
│   └── visualization/             # SpreaD3 / skyline HTML
├── report/                        # 模块 report
│   ├── phylo_summary.csv
│   └── phylo_report.html
└── .checkpoints/<stage>.done      # 断点续传标记（JSON，含完成时间 + 可选输入指纹）
```

---

## 七、断点续传与 stage 完成检测（checkpoint）

**机制**：每个病毒 work_dir 下维护 `.checkpoints/<stage>.done`（JSON）。

- **stage 成功完成** → 写 `.done`
- **重跑同命令** → 已完成 stage 直接跳过，日志显示 `[<stage>] checkpoint OK, skip`
- **产物兼容**：无 `.done` 但标志产物存在（见 §3.1 表）→ 同样视为已完成（兼容历史数据）
- **`--force`** → 忽略所有 checkpoint 重跑（产物存在时复用，不删除；无产物时重新计算）
- **失败 stage** → 不写 `.done`，重跑自动重试；下游因缺前置自动跳过

**输入指纹（P1-1，2026-09-15）**：`align` 与 `tree` 会把**输入文件指纹**（路径+大小+mtime 的 sha1）
记进 `.done`。换了 MAFFT 参数 / 加了新序列 / 改了 IQ-TREE 参数 → 指纹不符 → **自动重跑**，
不会拿旧产物静默复用。`tree` 的指纹里额外并入 `iqtree_param_token()`（建树参数 token）。

**状态追踪**：每次运行日志输出 `stages: done=[...] skipped=[...] failed=[...]`；
`report/phylo_summary.csv` 含 `Stages_done` / `Stages_skipped` / `Stages_failed` 三列。

```bash
# 中断后继续（已完成的 stage 自动跳过）
python phylo_pipeline.py --extract_dir 06_extraction/ --metadata metadata.tsv --stage all
# 强制重跑某阶段
python phylo_pipeline.py --work_dir gcva_results/ --alignment mafft.aln.fasta \
    --tree iqtree.treefile --meta dates.csv --stage tree --force
```

---

## 八、典型案例

```bash
# 1. 批量全流程
python phylo_pipeline.py --extract_dir 06_extraction/ --metadata metadata.tsv --stage all -t 20

# 2. 批量核心（比对 → 建树 → 时钟），自动补依赖
python phylo_pipeline.py --extract_dir 06_extraction/ --metadata metadata.tsv --stage tree,clock -t 20

# 3. 只预览执行清单（强烈建议先跑）
python phylo_pipeline.py --stage all --dry

# 4. work-dir 补下游（时间 + 定年 + 地理，含宿主通道）
python phylo_pipeline.py --work_dir gcva_results/ \
    --alignment mafft.aln.fasta --tree iqtree.treefile \
    --meta dates.csv --host_meta host.csv \
    --stage clock,beast,phylogeo,geo,host

# 5. 单病毒注册名模式
python phylo_pipeline.py --virus GCVA --stage phylogeny,time

# 6. 分基因定年（吃 capheine 产物，--genes_dir 缺省自动找）
python phylo_pipeline.py --work_dir gcva_results/ --alignment mafft.aln.fasta \
    --tree iqtree.treefile --meta dates.csv --stage beast,gene_dating

# 7. 历史 capheine 旧路径
python phylo_pipeline.py --work_dir gcva_results/ --alignment mafft.aln.fasta \
    --tree iqtree.treefile --meta dates.csv --stage gene_dating \
    --genes_dir 07_capheine/Cytorhabdovirus_sp._lycii_OR489165.1_v3/cawlign

# 8. BEAST 只生成 XML 不跑（--beast 恒为 BEAST 1.x）
python phylo_pipeline.py --work_dir gcva_results/ --alignment mafft.aln.fasta \
    --tree iqtree.treefile --meta dates.csv --stage beast --skip_beast_run

# 9. 关掉某个 stage（覆盖 --stage all 的展开）
python phylo_pipeline.py --extract_dir 06_extraction/ --metadata metadata.tsv \
    --stage all --disable rdp5,beast
```

---

## 九、独立工具调用（不进总管线）

```bash
cd ~/MMPV-RNA/virome_phylo_pipeline

# 正选择下游（消费 capheine 产物）
python -m utils.selection_suite fubar --capheine-dir <virus>/select/capheine/

# 重组检测 / mask（stage rdp5 的手工等价调用）
python -m utils.recombination_analysis --fasta aln.fasta --prefix my --outdir rec_out
python -m utils.mask_recombination --fasta aln.fasta --rdp5 rec_out/my.csv --mask --outdir mask_out

# 群体遗传
python -m utils.popgen_analysis --fasta X.fasta --metadata Y.csv --group location

# 链合并（phylogeo/beast 多链收尾的手工入口）
python -m utils.merge_results --virus GCVA
```

> **两个孤儿 CLI（未被 `phylo_pipeline.py` 引用，属预期）**：`utils/glm_predictors.py`、
> `utils/import_export.py` —— 它们是需要手工调用的独立工具，不是导入 bug。

---

## 十、与代码的同步声明

本文件的每一张表都由 `_consistency_check_20260916/verify_docs_vs_code.py` 对拍：

| 断言 | 对拍对象 |
|---|---|
| stage 清单 / 数量（17） | `STAGE_ORDER` |
| 分组别名 | `STAGE_GROUPS` |
| 依赖关系 | `STAGE_DEPS` |
| checkpoint 标志产物 | `STAGE_CHECKPOINT_FILES` |
| CLI 参数名 | `argparse` 实际注册的选项 |
| config 分区键 | `CONFIG_SECTION_ARG` + `config.yaml`（含孤儿键） |
| IQ-TREE 默认支持值 | `IQTREE_BOOT_DEFAULT` / `IQTREE_ALRT_DEFAULT` |
| 输出目录布局 | 代码里的 `os.path.join(out_dir, '<mod>/...')` 字面量 |

任何一项与代码不符 → 校验脚本非零退出。**改代码时请同步改本文件与 `STAGE_REFERENCE.md`。**
