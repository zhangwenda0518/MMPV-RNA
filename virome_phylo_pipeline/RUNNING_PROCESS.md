# phylo_pipeline 完整运行过程（叙事版）

> **📌 文档定位：本文件是「跑起来会发生什么」的叙事说明**，不重复表格。
> 参数、输出路径、config 分区、checkpoint 规则一律**以 `RUN_GUIDE.md` 为准**；
> stage 的逐项内部细节见 `STAGE_REFERENCE.md`。
>
> **同步状态**：2026-09-16 与代码逐项复核后重写。旧版（2026-08-27）把 `RUN_GUIDE.md`
> 的表**整份复制**了一遍，两份同源文档随后各自漂移（旧版还在讲 `rtt`/`temporal`/`genes`
> 三个已废弃 stage、`beast2.xml`、`temporal_signal/` 平铺路径）。**本次删掉全部重复表**，
> 只留叙事，从结构上消除「两份文档互相打脸」的可能。
> 改动前备份：`_docconv_20260916/backup/RUNNING_PROCESS.md`。

---

## 一、三种入口，同一套调度

| 模式 | 命令 | 适用 |
|---|---|---|
| **批量** | `--extract_dir 06_extraction/ --metadata Core14.tsv` | 多病毒全流程 |
| **单病毒注册名** | `--virus GCVA` | 已在 `datasets.yaml` 注册 |
| **work-dir** | `--work_dir <dir> --alignment ... --tree ... --meta ...` | 已有比对+树，补下游 |

三者共用 `process_virus()` 的同一份 stage 分发与 checkpoint 逻辑，行为一致。

---

## 二、启动阶段（所有模式共有）

1. 读 `config.yaml`（默认参数）→ CLI 显式传参覆盖。
2. **config 孤儿键自检**：`config.yaml` 里登记不到 dest 的键会打一行 stderr 告警
   （`phylogeo.chains` / `beast.chains` 就是这么被发现的），不再静默失效。
3. 解析 `--stage` → 展开大模块别名 → 按 `STAGE_DEPS` 补齐依赖 → 按 `STAGE_ORDER` 排序。
4. `--disable` 从最终清单里真正剔除。
5. 检测工具依赖（treetime / biopython / mafft / iqtree / Rscript）。
6. `--dry` 只打印执行计划后退出。

---

## 三、批量模式：逐 stage 发生了什么

> 括号内是该 stage 的模块归属（决定产物落在哪个目录）；路径细节见 `RUN_GUIDE.md` §六。

### `prep`（data）数据收集
`load_sample_metadata` 读全局元数据 → `scan_extraction_dir` 扫病毒目录 →
`prepare_virus_inputs` 为每个病毒构建：合并序列 fasta、采样日期、参考序列。
**同时写出三张通道表的初版**：`data/dates.csv`（时间）、`data/sample_metadata.csv`（地理）、
`data/host.csv`（宿主，只收有 host 的样品）。按样本数排序 → `--min_samples` 过滤 → `--max_viruses` 截断。

### `metadata`（data）入口治理
`govern_table()` 逐行体检：日期格式白名单 + 精度分档；地理层级校验 + 删 `Unknown` 层级；
占位符**先补算坐标、补不出才清空字段**。产出 `data/metadata/metadata_std.csv` 与
`metadata_report.tsv`（原件不动）。
随后 `rebuild_channels()` **回填三通道表**（不增不删行，骨架 = 原表 name 列 + 行顺序），
带来源指纹 checkpoint —— 来源未变直接 skip。若该病毒**没有任何宿主来源**，
`host.csv` **不落地**（fail-safe：绝不用空表覆盖既有文件）。
删除 `Unknown` 层级后只剩纯国家名的记录判为 `country_only` → 不造坐标、location 清空。

### `online`（data）NCBI 参考补充
① 按目录名里的 accession `efetch` 单参考 → `data/ncbi_ref/`；
② SeqHarvester 按学名搜公共序列 → `data/ncbi_harvest/`（含 metadata.csv 与 location/host 分组图）。

### `clean`（data）比对前清洗
长度/N 占比/大小写/非法字符/元数据缺失逐项过滤。产物 `data/clean/clean.fasta`，
**并把 `prep['combined_fasta']` 指向它** —— 后续 `align` 用的是清洗后的序列。
开了 `--clean_drop_no_date/-location` 但拿不到 metadata 表时，本 stage 整段降级为「用原序列」并告警。

### `align`（phylogeny）MAFFT 比对 + 饱和分析
`mafft --auto --thread N --reorder` → `phylogeny/mafft.aln.fasta`。
**替换饱和分析（Iss）内联在本 stage**，产物 `phylogeny/saturation_out/`。
checkpoint 带**输入指纹**：换了 MAFFT 参数/加了序列会自动重跑。

### `splitstree`（phylogeny）分裂网络
SplitsTree6 NeighborNet → `phylogeny/splitstree/neighbornet.nexus`，用于判别树状/网状信号。

### `rdp5`（recomb）重组检测
默认策略是 **删高置信重组子 + 重新 MAFFT**（产物 `recomb/rdp5/rdp5_dropped.realn.fasta`）；
`--rdp5_mask` 改为重组区置 N（产物 `masked/masked_N.fasta`）。
本 stage 的产物**替换 `align` 的输出喂给 `tree`** —— 这是「重组必须在建树前」的实现位置。

### `tree`（phylogeny）IQ-TREE 建树
默认 `-m MFP -B 1000 -alrt 1000 --quiet`（双支持值）。失败时按 **full → alrt_only → no_support**
三段降级，每次尝试与实际生效的支持值都写进 `phylogeny/iqtree.support.json`。
复用旧树前会**核对支持值是否齐全**；checkpoint 指纹含 IQ-TREE 参数 token，参数一变即重建。

### `popgen`（popgen）群体遗传
pypopart：单倍型鉴定 → π/θ/S/Tajima's D/Fu's Fs → MJN 网络 + VirNA 有向网络 → Fst（位点法+permutation）。
产物 `popgen/popgen_report.txt` 及其图表。

### `host`（popgen）宿主分化
**只读 `data/host.csv`** —— 没有该文件（即该病毒样本全都无 host）时**直接 skip 并记日志**，
不做任何宿主分析。产物 `popgen/host_analysis/tree_host.pdf`。

### `capheine`（select）正选择
从 GBK 提 CDS → CAWLign 参考锚定密码子比对 → IQ-TREE 基因树 → HyPhy（FEL/MEME/PRIME/BUSTED/
CONTRASTFEL/RELAX）→ DRHIP → MultiQC。产物 `select/capheine/`（cawlign/drhip）。

### `clock`（time）分子钟联合分析
把旧 `rtt` + `temporal` + `genes` 三个 stage 合并：TreeTime 根到梢回归 → `time/treetime_rtt/`；
TreeDater LTT → `time/treedater_ltt/`；
DRT 日期随机化检验 → `time/temporal_signal/`。DRT 结论（`passed`/`drt_status`）落进 metrics。

### `beast`（time）全长 BEAST 定年
**恒走 BEAST 1.x 路径**，产物名固定 `time/beast/beast1.xml`。
XML 用 UCLN + skyline、`discretize_locations=False`（**纯分子钟，无 location CTMC/BSSVS**）。
`--skip_beast_run` 只出 XML。多链（`--beast_chains`，config `beast.chains`）后台并行，
全部完成后 `utils/merge_results` 做 logcombiner + treeannotator。ESS 检查与 BSP/RSPP 图在 `postprocess/`。

### `gene_dating`（time）分基因定年
输入是 **capheine 产物**（`--genes_dir`，缺省自动找 `select/capheine/cawlign`）。
基因命名适配 → 逐基因 BEAST XML → `time/gene_dating/`（中间产物 `time/gene_split/`）。

### `phylogeo`（geography）系统地理
以采样地点为离散性状的 CTMC + BSSVS。`--beast_version` **只在这里生效**（1 或 2）。
多链跑完合并出 MCC 树 `geography/phylogeography/merged/mcc.tree`，SpreaD3 交互地图落 `visualization/`。
`--rrt_randomized_dir` 给了才跑 VirPhyKit 口径的 RRT（树注释法，需自备随机化 MCC 树）。

### `geo`（geography）地理分析
Mantel 检验（遗传 × 地理距离，默认 9999 次置换）+ 树地理着色 + VirSpaceTime →
`geography/geo_analysis/mantel_test.pdf`（扩展结果在 `extended/`）。**不需要 BEAST 产物**。

### `report`（report）汇总
`report/phylo_summary.csv`（含 `Stages_done/Skipped/Failed` 三列 + headline 指标，来源逐列标注）
+ `report/phylo_report.html`（按模块扫产物）。`--no-report` 可整体关闭。

### 全程日志
`<output_dir>/phylo_pipeline.log`，每病毒每 stage 带时间戳；
每个病毒结束打一行 `stages: done=[...] skipped=[...] failed=[...]`。

---

## 四、work-dir 模式：与批量的差异

1. 直接构造 `prep`：`combined_fasta` 取 `--alignment`、`tree_file` 取 `--tree`、
   `dates_csv` 取 `--meta`、宿主通道取 `--host_meta`（省略则由 `--meta` 拆出）。
2. 随后走**完全相同的 stage 调度**：`align`/`tree` 因产物已在而 checkpoint skip，不重算。
3. 每病毒结束后同样跑 `report`。

```bash
python phylo_pipeline.py --work_dir gcva_results/ \
    --alignment mafft.aln.fasta --tree iqtree.treefile \
    --meta dates.csv --host_meta host.csv \
    --stage clock,beast,phylogeo,geo,host -t 8
```

---

## 五、断点续传

规则、标志产物对照表、输入指纹机制见 **`RUN_GUIDE.md` §七**（本文件不再重复，避免二次漂移）。

要点：`.checkpoints/<stage>.done` 存在即跳过；`--force` 忽略标记重跑；失败 stage 不写标记；
`align`/`tree` 带输入指纹，输入或建树参数一变即自动重跑。

---

## 六、历史运行记录（2026-08-27，**已过期，仅存档**）

> ⚠️ 以下日志用的是**当时的 stage 命名**（`rtt` / `temporal` / `genes` 后已并入 `clock`），
> 产物路径也是当时的平铺布局（后已改为 8 大模块目录）。**不要照抄**。

```
Stages: prep, align, tree, rtt, temporal, beast, geo, host, genes
[align] checkpoint OK, skip
[tree]  running IQ-TREE on mafft.aln.fasta...   (20 seq, 3s)
[rtt]   β=0.00064, R²=0.1811
[geo]   Mantel r=-0.0141 (n=20); 7 regions
[host]  3 hosts
[temporal] DRT: Real R²=0.1811 > 100% randomized → PASSED
stages: done=[tree,rtt,geo,host,genes,temporal,beast] skipped=[align] failed=[]
```

当时的命令（**已废弃写法**，仅作对照）：

```bash
# ❌ 旧写法：rtt/temporal/genes 已并入 clock；--host_meta dates.csv 是当时的临时做法
python phylo_pipeline.py --work_dir gcva_results/ \
    --alignment mafft.aln.fasta --tree iqtree.treefile \
    --meta dates.csv --host_meta dates.csv \
    --stage tree,rtt,temporal,geo,host,genes --beast --skip_beast_run \
    --check_temporal --drt_randomizations 5
```

现行等价命令见 `RUN_GUIDE.md` §八。
