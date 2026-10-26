# 17 · 服务器端脚本与本仓版本比对（2026-10-22）

> 目的：回答「本地新加的 eve-mask 阶段，会不会和服务器上实际在跑的生产脚本打架」。
> 方法：把服务器 `/home/zhangwenda/MMPV-RNA` 的关键脚本 scp 下来，与本仓 `HEAD`（`0a52f29`）和工作区逐文件比对（`git diff --no-index --ignore-cr-at-eol`、`git hash-object`、`git cat-file`）。
>
> **修订提示（2026-10-22 晚）**：本文只比了 3 个文件，范围窄；全仓 1,491 vs 2,058 个脚本的逐文件比对见
> [`doc/18-full-script-diff-20261022.md`](18-full-script-diff-20261022.md)。本文 §2.1/§2.2/§2.3 里原写作「服务器独有旧版 / 服务器独有内容」的
> 若干条目经复核实为 **服务器 = 本仓 `24b32f7` 旧版**（内容在本仓对象库中存在，只是 CRLF 归一后哈希不同），已在下方就地更正。

## 0. 结论（先看这个）

1. **文件层面不冲突。** 本次改动只碰 7 个文件（4 个新增 + 3 个修改），其中与服务器**同路径**的只有 `public_metadata_pipeline/build_host_pipeline.py` 一个；服务器上这份文件是 **本仓 `24b32f7` 的旧版**（blob `58eb7c7d`，16,029 B，CRLF），不是服务器侧独立改动，所以不存在「两边同时改同一个文件」这种冲突。
2. **版本层面落后 3 周。** 服务器停在 `b8a73f9`（2026-09-28，「sync from local 24b32f7」，685 文件），本地 `main` 已到 `0a52f29`（2026-10-13），中间 **37 个提交**；服务器上那份 `build_host_pipeline.py` 还是更旧的变体（见 §2.1）。
3. **真正的风险是「整树同步」，不是 eve-mask 本身。** 服务器工作区是一份**未提交的大改**：
   - **3421 个已跟踪文件被删除**（`virome_analysis_pipeline` 1864、`virome_phylo_pipeline` 1267、`search_results` 244、`metadata_gui` 21…）
   - **153 个文件被修改**（phylo 66、discovery 29、analysis 27、submission 10、public_metadata 6…）
   - **196 个未跟踪文件**（phylo 66、analysis 36、discovery 31、`data_preprocessing_pipeline` 4…）
   - 其中 `generate_report.py`、`gsa_sra.plot.py`、`gsa_sra.search.py`、`preprocess_report.py`、`mmpv_common/io_layout.py` 有**服务器独有内容**（本仓 git 对象库里没有对应 blob，⚠ 经 doc/18 复核：io_layout 的差异实为本地注释/docstring 改动，见 doc/18 §4.8）；而 `pipeline_config.yaml`、`host_depletion.py`、`preprocess.py`、`data_preprocessing.py`、`build_host_pipeline.py` 只是**本仓 `24b32f7` 的旧版**（CRLF 归一后哈希一致，非服务器独有）。
   → 直接 `git pull` / `git checkout .` 会把这些现场冲掉。**先存档，再谈同步。**
4. **一处行为差异值得你确认**：服务器上那份 `host_depletion.py`（7 月 12 日 mtime、服务器侧 untracked，内容实为本仓 `24b32f7` 旧版，LF blob `57cb46b5`）在 `--keep-taxids` 打开时**无条件跳过 align**（默认 `10239` → 默认不跑 bowtie2/hisat2）；本仓版本是「只有没给 `-x/--step2_index` 时才跳过」。而三个入口脚本（`data_preprocessing.py` / `preprocess.py` / `preprocess_unified.py`）都只传比对索引、**不传** `--keep-taxids`，所以在服务器那版上，这一路实际只跑 Kraken2 一层。这会影响 `doc/16` 里「Stage 4 比对删掉 EVE 区 reads」的结论是否适用于已跑的 1,348 个样本（见 §3）。

## 1. 版本锚点

| | 本地 | 服务器 |
|---|---|---|
| HEAD | `0a52f292` · 2026-10-13 · `docs: README/wiki/doc 全量同步至 v3.1 …` | `b8a73f907` · 2026-09-28 · `sync from local 24b32f7 … (685 文件)` |
| 上一提交 | `96bd30a` 论文图脚本 | `e21d9cb` pre-sync backup（59 个服务器版本打包 `_bak_pre_sync_20260928/pre_sync_67.tar.gz`） |
| 同源点 | `24b32f7` · 2026-09-28 · v3.1 统一 I/O 布局 | 服务器从本地 `24b32f7` 同步而来 |

关系：`服务器 ≈ 本地 24b32f7 的快照 + 服务器侧未入库改动`；本地在其后又走了 37 个提交。所以**服务器不是另一个分支，是一份过期快照**。

## 2. 两个目标管线逐文件比对

### 2.1 `public_metadata_pipeline/`

| 文件 | 服务器 | 本地 HEAD | 本地工作区 | 差异性质 | 处置 |
|---|---|---|---|---|---|
| `build_host_pipeline.py` | 16,029 B · `147eb2b`（CRLF；等价 LF blob `58eb7c7d`） | 15,889 B · `682ce453` | 37,740 B（+443/−17，含 eve-mask 阶段） | **服务器 = 本仓 `24b32f7` 旧版**（非服务器独有）：①`datasets download …` 少了 `--reference`；②把 `# Merged host genome` / `# Created` 注释行**直接写进合并 FASTA**（本仓已改成旁车 `.meta.txt`，代码注释写明 `#` 行会让 `bowtie2-build` 拒绝） | **用本地新版整体替换**；不要把 eve-mask 补丁往服务器旧版上打 |
| `generate_report.py` | 69,008 B · `212efa0a` | 67,721 B · `7accd16a` | 72,140 B | 两边都大改，服务器内容本仓没有 | **三方合并**，禁止直接覆盖 |
| `gsa_sra.plot.py` | 17,160 B · `aa47647f` | 15,489 B · `7e4bdb91` | 17,936 B | 两边都改 | 三方合并 |
| `gsa_sra.search.py` | 33,854 B · `71cb7cb4` | 33,081 B · `8f3693b9` | 47,767 B | 两边都改 | 三方合并 |
| `gsa_sra.info.py` | 9c018702 | `c6c94316` | — | 服务器 = **24b32f7 旧版**（内容旧，不是新改） | 跟随后续同步即可 |
| `sra2fastx.py` / `utils/pipeline_utils.py` | 与 HEAD 同内容（`e58dd192` / `0ced7bdc`） | 同左 | — | 只有 `644 → 755` 的 mode 变化 | 无实质差异 |
| `doc.md` | 18,246 B · `b89de561` = 24b32f7 版 | 18,448 B · `9CCFB09` | 21,610 B（含 eve-mask 文档节） | 服务器旧 | 随同步更新 |

### 2.2 `data_preprocessing_pipeline/`

| 文件 | 服务器 | 本地 HEAD | 本地工作区 | 差异性质 | 处置 |
|---|---|---|---|---|---|
| `host_depletion.py` | 63,007 B · `e2b952f9`（CRLF，untracked；等价 LF blob `57cb46b5` = `24b32f7`，mtime 2026-07-12） | 64,107 B · `b347004` | 65,236 B（+20 行 `figure_export` 钩子，2026-10-20 未提交） | **服务器 = 本仓 `24b32f7` 旧版**（非独有），但 `--keep-taxids` 语义不同（§3） | 需先确认生产跑的是哪版，再决定是否替换 |
| `preprocess.py` | 7,339 B · `a81793f7`（untracked） | 7,161 B · `3391f66b` | 9,078 B | 服务器 = `24b32f7` 旧版（内容在本仓对象库中存在） | 随同步更新 |
| `preprocess_report.py` | 27,421 B · `f6f7cccb` | 26,675 B · `cf3b43e3` | 53,317 B | 两边都大改（本地 2 倍扩写） | 三方合并 |
| `README.md` | 5,483 B · `369da947`（untracked） | 5,401 B · `94d2cc03` | 5,483 B | 尺寸相同，疑似仅换行/局部差异 | 低风险，人工看一眼 |
| `data_preprocessing.py` | 19,339 B · `23444de1` | 19,383 B · `201a1965` | 19,816 B | 服务器 = `24b32f7` 旧版（非独有，+14/−8 为本地修复） | 以本地为准覆盖 |
| `run_bbnorm.py` | `d4c31e4a` = 24b32f7 版 | 同左 | 同左 | 无差异 | — |
| `build_hostbase.py`（public 侧） | `c2360448` = 24b32f7 版 | 同左 | 同左 | 无差异 | — |
| `doc/03-host_depletion.md` | 2,794 B（旧） | 3,331 B | 5,192 B（含「宿主库版本与 EVE 掩码库」节） | 服务器旧 | 随同步更新 |

### 2.3 服务器上其它未入库改动（不属本次范围，但会阻塞同步）

- 删除 3421：`virome_analysis_pipeline` 1864、`virome_phylo_pipeline` 1267、`search_results` 244、`metadata_gui` 21、`stats` 7、`submission_gui` 6…
- 修改 153：`virome_phylo_pipeline` 66、`virome_discovery_pipeline` 29、`virome_analysis_pipeline` 27、`virome_submission_pipeline` 10、`public_metadata_pipeline` 6…
- 未跟踪 196：`virome_phylo_pipeline` 66、`virome_analysis_pipeline` 36、`virome_discovery_pipeline` 31、`biosoft` 7、`scripts` 4、`data_preprocessing_pipeline` 4…
- **有服务器独有内容（本仓对象库查不到 blob）**：`generate_report.py`、`gsa_sra.plot.py`、`gsa_sra.search.py`、`preprocess_report.py`（⚠ `mmpv_common/io_layout.py` 经 doc/18 §4.8 复核为本地注释改动，不算独有）。
  全仓这类「两边都改过」的文件共 **43 个**（清单见 [`doc/18-full-script-diff-20261022.md`](18-full-script-diff-20261022.md) §4）。
- **看起来独有、实为本仓 `24b32f7` 旧版**（CRLF 归一后哈希一致）：`build_host_pipeline.py`、`host_depletion.py`、`preprocess.py`、`data_preprocessing.py`、`pipeline_config.yaml` —— 直接用本地版覆盖即可，**不需要三方合并**。

## 3. `--keep-taxids` 行为差异（重点）

服务器版（2026-07-12）`host_depletion.py` 的 `main()`：

```python
# --- 智能步骤调整: --keep-taxids 启用时自动跳过 align ---
if args.keep_taxids:
    steps_list = [s for s in args.steps.split(',') if s != 'align']
    sys.stderr.write("💡 --keep-taxids 已启用，自动跳过 align 阶段 …")
```

本仓 `HEAD`（`0a52f29`）：

```python
# 自建宿主库(提供 -x/--step2_index) → Kraken2 + align 双层去除;
# 通用库如 k2_pluspfp(无比对索引) → 仅 Kraken2 去除。
if args.keep_taxids and not args.step2_index:
    ... 去掉 align ...
```

两版 `--keep-taxids` 默认值都是 `'10239'`（Viruses）、`--steps` 默认 `kraken2,align,rrna`。

**调用方证据**（2026-10-22 直接在服务器 `~/MMPV-RNA` 上 grep 三个入口）：

| 入口脚本 | 传给 `host_depletion.py` 的参数 | 传 `--keep-taxids`？ |
|---|---|---|
| `data_preprocessing_pipeline/data_preprocessing.py` | `--step2_index {host_align_db}` | 否（`grep -c keep.taxids` = **0**） |
| `data_preprocessing_pipeline/preprocess.py` | `--step2_index {args.host_align_db}` | 否 |
| `public_metadata_pipeline/preprocess_unified.py` | `-x {step2_index}` | 否 |

→ 三个入口都**只给比对索引、不给 `--keep-taxids ''`**，于是服务器 7 月版里 `if args.keep_taxids:` 恒为真、align 被无条件剔除。**若 1,348 样本确实由这三个入口在服务器上跑的，则生产去宿主只有 Kraken2 一层**，稿件 Methods 的「Bowtie2 precise alignment」需要按实际运行的脚本版本改口径（或说明该步在某些批次被跳过）。

**后果**：
- 服务器版：默认参数下 **align 永远不执行**（除非显式 `--keep-taxids ''`，那样也同时失去病毒救援）→ 去宿主只剩 Kraken2 一层（`--confidence` 默认 0.4）。
- 本仓版：给了 `-x/--step2_index` 就会跑 Kraken2 + bowtie2/hisat2 双层，EVE 同源 reads 会在 Stage 4 被当宿主删掉（`doc/16` §2 的三个假阴入口正是基于这个版本）。
- 如果 1,348 个样本是用服务器 7 月版跑的，那么稿件 Methods 里「Bowtie2 precise alignment」需要核对是哪一次/哪一版跑的；EVE 假阴的主因也会从「Stage 4 删除」变成「Kraken2 分类层误判」。

**怎么确认（在样本输出目录执行）**：

```bash
# 7 月版会打这句
grep -rl "已启用，自动跳过 align" <样本输出根>/logs 2>/dev/null | head
# 本仓版会打这句
grep -rl "未提供比对索引" <样本输出根>/logs 2>/dev/null | head
```

> 注（2026-10-22 实测）：服务器 `~/MMPV-RNA` 下没有 `02_Preprocessing`；`/data`、`/alldata`、`/mnt` 均为空目录或不存在；磁盘只有 `/dev/sda2` 500G（已用 117G）挂在 `/`。样本输出**不在这台登录机上**，可能在 `241` 或其它生产节点，需要你给出实际路径。

## 4. 本次改动清单（本地未提交，供核对）

| 文件 | 状态 | 说明 |
|---|---|---|
| `public_metadata_pipeline/build_host_pipeline.py` | 修改 +443/−17 | 新增可选 `eve-mask` 阶段，默认关闭；`build_parser()/plan_stages()/main()` 拆分；`db_provenance.json` |
| `public_metadata_pipeline/mask_host_genome.py` | 新增 408 行 | BED / verdict TSV → 区间并集 → 掩码（N 或小写），长度守恒 |
| `public_metadata_pipeline/tests/test_mask_host_genome.py` | 新增 | 14 tests |
| `public_metadata_pipeline/tests/test_eve_mask_stage.py` | 新增 | 21 tests（合计 35 passed） |
| `doc/16-host-eve-masking.md` | 新增 | 决策文档（结论 / 假阴入口 / P0-P2 / 验收 / 稿件措辞） |
| `doc/03-host_depletion.md` | 修改 | 追加「宿主库版本与 EVE 掩码库」节 |
| `public_metadata_pipeline/doc.md` | 修改 | 三阶段编排 + eve 参数 + 上下游衔接 |

服务器上：`mask_host_genome.py`、两个测试、`doc/16` 均不存在；`build_host_pipeline.py` 是旧版（§2.1）。

## 5. 部署方案

### 方案 B′（推荐，最小侵入）：只搬 eve-mask

```bash
# 1) 服务器备份旧版（只为对照；它是旧版，不建议保留）
cp public_metadata_pipeline/build_host_pipeline.py \
   public_metadata_pipeline/build_host_pipeline.py.bak_20261022

# 2) 本地 → 服务器（只这两个文件）
scp public_metadata_pipeline/build_host_pipeline.py \
    public_metadata_pipeline/mask_host_genome.py \
    success:~/MMPV-RNA/public_metadata_pipeline/

# 3) 冒烟
python public_metadata_pipeline/build_host_pipeline.py --help
python public_metadata_pipeline/build_host_pipeline.py --stage all --dry-run      # 期望 eve mask: OFF，plan 不含 eve-mask
python public_metadata_pipeline/build_host_pipeline.py --stage eve-mask --dry-run  # 期望报"需要掩码源"
```

注意事项：
- **不要盲目同步 `pipeline_config.yaml`**：服务器那份是本仓 `24b32f7` 旧版（不是服务器独有），本地版本更新（新增 `--graph-enabled` 说明、EVE 段未接线注释等）；同步前先确认生产 DB 路径没在服务器侧被改过（`git diff b8a73f9 -- pipeline_config.yaml` 对着看）。
- 新版 `run_genome_down()` 带 `--reference`；若当初在服务器删掉它是因为某物种没有 reference assembly 而下载失败，需要按物种试一次（可改成「先带 `--reference`，失败回退」）。
- 掩码库只写 `<work>/hostdb_evemask/`，不覆盖 `hostdb/`；下游用 `--host_db <work>/hostdb_evemask` 显式指定。
- **EVE 管线在服务器上是旧布局**：服务器为 `~/MMPV-RNA/endogenous_virus_pipeline/`（顶层，v3.1 归档前的位置），本仓是 `archive/20261011_endogenous_virus_pipeline/`。`build_host_pipeline.py --eve-scan` 已做双路径自动定位（`archive/...` → `endogenous_virus_pipeline/`），两处都没有时用 `--eve-scan-script` 显式给。
- 服务器上**没有**任何既有的 EVE 掩码实现（`grep -l mask` 只命中 `build_host_pipeline.py` 的 API-key masking 与 `virome_phylo_pipeline/mask_recombination.py` 重组检测）→ eve-mask 与服务器现有代码无功能重叠。

### 方案 A（全量同步，暂不推荐）

前提：先把服务器现场固化，再合并。

```bash
cd ~/MMPV-RNA
git checkout -b server-wip-20261022
git add -A && git commit -m "wip: server-side scripts before v3.1 sync (2026-10-22)"
```

然后与本地 `main` 做合并，**必须人工过**的文件：`preprocess_report.py`、`generate_report.py`、`gsa_sra.plot.py`、`gsa_sra.search.py`、`host_depletion.py`（`--keep-taxids` 语义），以及 discovery 侧 2026-10-18/19 的 `benchfix` 文件（`virome_pipeline.py`、`report_pipeline.py`、`cobra_pipeline.py`）；`io_layout.py` 已确认只是本地注释改动、`pipeline_config.yaml` 只是文档性说明，可直接以本地为准。`virome_analysis_pipeline` / `virome_phylo_pipeline` 的 3421 删除 + 130 未跟踪需要先判断是「有意重排」还是「误删」。全量清单见 [`doc/18-full-script-diff-20261022.md`](18-full-script-diff-20261022.md) §3–§4。

## 6. 待确认

1. **1,348 样本用的是哪版 `host_depletion.py`**（§3 的 grep 命令）→ 决定 `doc/16` §2 的假阴入口是否要改写。
2. **服务器版 `build_host_pipeline.py` 为何删掉 `--reference`**：是有意（某物种下载失败）还是旧版没跟上？找当时的报错日志。
3. **`metadata_gui/**` 等 3421 个删除是否有意**：若有意，b8a73f9 之后的同步会把这些文件全部恢复。
4. **`preprocess_report.py` 的服务器版改动是否要回灌本仓**（本地已 2 倍扩写到 53 KB；服务器版可能有本仓缺的字段名修复）；`io_layout.py` / `pipeline_config.yaml` 已核实为「本地注释/文档改动 vs 服务器 24b32f7 旧版」，回灌意义不大。

## 附：复核命令（本机即可）

```bash
# 服务器文件 hash（在服务器执行）
git hash-object public_metadata_pipeline/build_host_pipeline.py

# 本仓某个 blob 是否存在/属于哪个版本
git cat-file -t <hash> && git log --oneline --all --find-object=<hash> | head

# 两个文件逐行差异（忽略 CRLF）
git diff --no-index --ignore-cr-at-eol base.py server.py
```
