# virome_phylo_pipeline 全面审查报告

> **【修复状态 2026-10-13】** 本报告所列问题已在 10-13 修复轮处理：
> P0-1 (MAFFT 假成功/align 早退) ✅ 已修 — `provided_alignment` 门控, batch 必落到 MAFFT；
> P0-2 (clean 静默删序列) ✅ 已修 — `utils/seq_clean.py` 仅 drop 门槛开启时判缺行；
> P0-3 (publication_figures 空 tree/未比对输入) ✅ 已修 — batch 按病毒名取 align/tree 真实产物；
> P0-4 (--virus 退出码) ✅ 已修 — 且 10-12 重写中丢失的整套退出码/blocked/Incomplete 记账已恢复 (含 CSV 两列与 `_beast_ran` 闸门)；
> P1-1 (--strip-illegal 无效) ✅ 已修；P1-2 部分 ✅ (capheine ncbi_ref / work_dir beast / gene_align_generator 三处硬编码目录已换 ph_dir)；
> P1-3 未修 (clock 的 summary_rows 非空即 mark)；P1-4 ✅ 已修 (行数闸门恢复)；P1-5 未动。
> 详见 `outputs/MMPV_FULL_REVIEW_20261013/REVIEW_SUMMARY_20261013.md`。

- 审查日期：2026-10-09
- 审查对象：`D:\桌面\延伸基因组\MMPV-RNA\virome_phylo_pipeline`
- 审查方式：**只读**。本次未修改任何 `.py` / `.R` / `.sh` / `.md` 源文件。
- 规模：活代码 83 个 Python 文件、约 49.6k 行；主编排器 `phylo_pipeline.py` 4456 行
- 环境：`C:\Python312\python.exe`（依赖齐全）；本机缺 `mafft` / `iqtree` / `beast` / `hyphy`

**证据标注约定**：`✅已实测` = 本次实际执行并观察到结果（附命令或日志行号）；`深读` = 源码控制流推断，未构造数据跑。

---

## 0. 总体评价

**做对了什么（先说好的）**

这套流水线的"失败可检测性"是我见过的同类科研代码里做得相当好的一档：

- **checkpoint 带输入指纹**（`phylo_pipeline.py:251-306`）：不只是判"文件在不在"，还记录输入文件的 `路径+大小+mtime` sha1；`tree` 额外并入 IQ-TREE 参数 token（`:3074-3075`）。换了 MAFFT/IQ-TREE 参数会自动重跑，不会拿旧产物静默复用。
- **`_skip` 五态 + 三档退出码**（`:3091-3108`、`:325-340`）：`checkpoint / na / pending / degraded / blocked` 分开记账，`degraded`/`blocked` 汇总进 `Stages_blocked` 并置 `Incomplete=True`，退出码 failed→1、blocked→2、干净→0。日志实证这套机制真的在工作：`:699` `[geo_paths] skip(blocked)` → `:702` `EXIT 2: 结果不完整`。
- **`validate_stage_graph()`**（`:422-447`）：因为真实执行顺序 = `if '<stage>' in stages:` 代码块的**物理顺序**，它用 `STAGE_ORDER` 的下标钉住"依赖必须排在使用者之前"。这是个很有洞察力的守卫。
- **`config_orphan_keys()`**（`:509-535`）：把"config 键没登记映射就被静默丢弃"变成启动期显式告警。实测当前 `config.yaml` 孤儿键 = `[]`，干净。
- **`_beast_ran()` 二次判定**（`:3468-3472`）：防止 `--skip_beast_run`（只出 XML）被当成"已跑过 BEAST"。
- **`seq_clean` 的 fail-loud 护栏**（`utils/seq_clean.py:295-302`）：开了 `--drop-no-date` 但没识别出日期列时**拒绝执行**而不是整批误删。
- **`--explicit_cli_dests`**（`:538-565`）：让"CLI 显式传参 > config"这条文档承诺真正生效（旧版只要传 `--config` 就会被 yaml 静默盖掉）。

**系统性问题（本报告的主线）**

> **同一份数据存在多条路径，其中一条路径被静默绕过。**

三处 P0 都是同一个形状：代码里有一段"针对 A 模式写的守卫/分支"，但它**没有把 A 模式限定住**，于是在 B 模式下误触发，而下游全部按"成功"继续。

| 编号 | 一句话 | 后果 |
|---|---|---|
| P0-1 | `align` 把输入文件当"比对产物"返回 → **MAFFT 从不执行** | 下游全在未比对序列上跑 |
| P0-2 | `clean` 把"metadata 无此行"当删除理由 → **静默删序列** | 序列集被无声裁剪 |
| P0-3 | batch 收尾把未比对序列当比对喂 `publication_figures` | 发表级图件静默缺失 |
| P0-4 | `--virus` 模式**退出码恒为 0** | 自动化误判成功 |

配套原因：文档声明的机器校验脚本 `_consistency_check_20260916/verify_docs_vs_code.py` **已不存在**，所以文档漂移无人拦截（见 §4）。

---

## 1. P0 — 会导致结论错误或静默失效

### P0-1 ✅已实测：`align` 不跑 MAFFT，`run_mafft` 调用点不可达

**证据链**

`phylo_pipeline.py:1229-1232`：
```python
provided = prep.get("combined_fasta", "")
if provided and os.path.exists(provided) and provided != aln_file:
    log.info(f"  [align] {vname}: using provided alignment: {Path(provided).name}")
    return {"success": True, "alignment": provided}      # ← 任何已存在的 combined_fasta 都被当"比对"
```
`phylo_pipeline.py:1243-1245`：
```python
if os.path.abspath(provided) == os.path.abspath(aln_file):
    return {"success": False, "error": "alignment file empty/corrupt (work-dir mode, no raw input to recompute)"}
```
`phylo_pipeline.py:1249`：`ok = run_mafft(provided, aln_file, args.threads, log)`

两条守卫合起来封死了 `:1249`：
- `provided` 存在且 ≠ `aln_file` → `:1230` 早退；
- `provided` 不存在 → `:1238` 早退；
- `provided` == `aln_file` → `:1243` 早退。

**→ `run_mafft(...)`（`:1249`）没有任何可达路径。MAFFT 在本管线中只被 `rdp5` 的"删重组子后重比对"（`:1326`）使用。**

**实测**（造 5 条长度不等的未比对序列，走 `clean → align`，即 `process_virus` 的真实调用序）：
```
[clean] success= True clean_fasta= clean.fasta
[align] probe: using provided alignment: clean.fasta
[align] 返回: {'success': True, 'alignment': 'clean.fasta'}
[align] 返回的 alignment 是否就是 clean.fasta(未比对)? True
   mafft.aln.fasta 存在? False
```

**为什么这改变结论**

`--extract_dir`（batch）是三个入口里**唯一真正需要 MAFFT 的**（`prepare_virus_inputs` 写出的是各样本全长原始序列，`utils/data_collector.py:330-353`），而它恰好在这一模式下被跳过：

- `phylogeny/mafft.aln.fasta` **永不产生**；
- `tree`（IQ-TREE）、`clock`（TreeTime/Treedater）、`popgen`、`splitstree`、`distance_matrix`、`align_qc` 全部吃**未比对**输入；
- `STAGE_CHECKPOINT_FILES['align'] = 'phylogeny/mafft.aln.fasta'`（`:227`）也永不成立，align 的 checkpoint 完全依赖 `.done` 标记。

真实运行日志已复现该后果（`phylo_results/phylo_pipeline.log`）：
```
:553  INFO   [align] smoke_20261003: using provided alignment: clean.fasta
:569  ERROR TreeTime-RTT failed: ValueError: Sequences must all be the same length
:607  INFO   [smoke_20261003] stages: done=['align', 'saturation', 'tree', 'clock'] ...
```
即：`align` 报成功 → `clock` 报成功，而 TreeTime 因序列不等长直接失败。

**与文档冲突**：`STAGE_REFERENCE.md:310-312` 明确承诺「输入 = clean 的 `clean.fasta`；处理 = `mafft --auto --thread N --reorder`；输出 = `phylogeny/mafft.aln.fasta`」。**代码与文档在此处正相反。**

**旁证**：`--virus` 模式下 `combined_fasta` 来自 `datasets.yaml` 的 `mafft.aln.fasta`（**已是比对文件**，见 `datasets.yaml:67`），此时早退是**正确**的——这正是该守卫的本意。问题是它没被限定在 `--virus` / `work_dir` 模式内。

**附带**：`clean` 会用清洗产物覆盖 `combined_fasta`（`:3188`），于是 `work_dir` 模式下用户传的 `--alignment` 也会被换成 `clean.fasta`，再被 `:1230` 当作"比对"返回。

---

### P0-2 ✅已实测：`clean` 以「metadata 无此行」为由静默删除序列

**证据链**

`utils/seq_clean.py:356-361`：
```python
row = meta.get(_norm_header(sid))
has_date = has_loc = None
if metadata_path:
    if row is None:
        has_date = has_loc = False
        issues_this.append("metadata 无此行")     # ← 与 drop_no_date/location 无关
```
`utils/seq_clean.py:381`：`status = "REMOVE" if issues_this else "KEEP"`

而 `run_stage_clean` **默认就会传** metadata（`phylo_pipeline.py:1173-1175`）：
```python
meta_path = (getattr(args, "clean_metadata", None) or prep.get("meta_csv") or None)
```

**实测**（FASTA id 带 accession 后缀 + metadata 用裸 id，即 `--virus` / `work_dir` 的常见形态）：
```
[用例1] 保留=0 剔除=2 理由=['metadata 无此行']
```

**为什么这改变结论**

- 触发条件不是"数据缺日期/缺地理"，而是**FASTA 的 id 与 metadata 的 `name` 对不上**（`_norm_header` 规范化后仍对不上）。这与用户是否开 `--clean_drop_no_date/-location` **完全无关**。
- `--virus` 模式（比对文件与 metadata 各自独立）与 `work_dir` 模式（用户自备 `--alignment` + `--meta`）最容易被命中；命中即整批删除。
- 且 `n_keep == 0` 时 `run_stage_clean` 返回失败（`:1204-1207`），`process_virus` 记 `degraded` 并**降级使用原始序列**（`:3192`）——于是"既没清洗、也没人知道清洗被跳过了"。
- batch 模式因为 `prepare_virus_inputs` 用同一个 `sample_id` 同时写 fasta 与三张通道表（`utils/data_collector.py:336/362/384`），天然对齐，所以这个坑长期没暴露。

**根因**：`prep` 侧的 ID 规范化走 `utils/seq_ids.py::base_sample_id`（`data_collector.py:311`），`clean` 侧另有一套 `_norm_header`。**同一个"ID 是否同一行"的判据有两处实现**——正是本项目自己反复记录过的那类问题。

---

### P0-3 ✅已实测：batch 模式 `publication_figures` 必然失败

**证据链**

`phylo_pipeline.py:4422-4424`（batch 收尾）：
```python
_cmd = ["--rtt_fasta", vq.get("combined_fasta", ""),
        "--rtt_dates", v_meta,
        "--rtt_tree", vq.get("tree_file", "")]
```
- `vq` 就是 `prep` 同一对象；`process_virus` 已把它改写成 `clean.fasta`（`:3188`，清洗产物 = **未比对**），**从未**指向比对产物；
- `prepare_virus_inputs` 返回的 dict **没有 `tree_file` 键**（`utils/data_collector.py:409-425`）→ `--rtt_tree` 恒为空串。

`publication_figures.py:72`：`aln = AlignIO.read(fasta_file, "fasta")` —— 要求所有序列等长；`:451` 的 help 也写着 `"FASTA alignment for TempEst"`。

**实测**：
```
未比对输入  → ValueError: Sequences must all be the same length  (publication_figures.py:72)
已比对输入  → 正常产出 Fig1_temprest.pdf
```

**对照**：`work_dir` 分支传的是 `_resolve(args.alignment)`（已比对，`:4301`）、`--virus` 分支传 `fasta`（datasets.yaml 的比对文件，`:4239`）——**只有 batch 模式坏**。

**为什么这改变结论**：batch 模式下 TempEst RTT 图、MCC 树图、参数附表（Table S1）**静默缺失**，只留一行 WARNING（`run_postprocess` 失败仅告警，`:2901-2903`）。`report/phylo_summary.csv` 仍会写、`Success` 仍可能为 True。

---

### P0-4 ✅已实测（静态）：`--virus` 模式退出码恒为 0

`pipeline_exit_code()`（`:325`）只在 **work_dir** 分支被调用（`:4316`）；batch 分支内联了等价逻辑（`:4440-4452`）；而 `--virus` 分支以裸 `return` 结束（`:4253`），既不调用也不 `raise SystemExit`。

**→ `--virus` 模式下即使有 stage `failed`，进程退出码也是 0**，与 `RUN_GUIDE.md:288-290`「`main` 退出码：有 failed→1，有 blocked/pending→2，干净→0」的承诺不符。`&&` 链、批处理驱动、CI 都会把失败当成功。

---

## 2. P1 — 健壮性 / 可移植性 / 可复现性

### P1-1 ✅已实测：`--strip-illegal` 选项无效
`utils/seq_clean.py:331-332` 在**无条件**把"非法字符"记入 `issues_this`，`:336-340` 才执行剥离，而 `:381` 只看 `issues_this` 是否非空。
实测：`strip_illegal=True` + 序列含 `*` → `保留 0 / 剔除 1 理由=['非法字符 * x1']`。该开关等于没有。

### P1-2 ✅已实测：默认 I/O 布局下多处绕过 `ph_dir()` 的硬编码目录
`mmpv_common/io_layout.py` 的 `DEFAULT_LAYOUT = "standard"`，映射为 `ph_data→01_data`、`ph_time→06_time`、`ph_report→08_report`。以下位置写死了 legacy 名：

| 位置 | 硬编码 | 应为 | 后果 |
|---|---|---|---|
| `phylo_pipeline.py:1682-1683` | `Path(out_dir,"data","ncbi_ref")` | `ph_dir('data')` | capheine 找不到 `online` 产出的 GBK 参考（`online` 写的是 `:1069` 的 `ph_dir('data')/ncbi_ref`） |
| `phylo_pipeline.py:4281` | `str(work_dir),"report"` | `ph_dir('report')` | work_dir 模式报告落错目录（`--virus` 分支 `:4213` 用的是 `ph_dir`） |
| `phylo_pipeline.py:4291` | `str(work_dir),"time","beast"` | `ph_dir('time')` | work_dir 模式探测不到 beast 日志/MCC → Table S1、MCC 图静默缺失 |
| `gene_align_generator.py:115,122-124,133,139,147-148,258` | `"time"` / `"select"` / `"data"` | 同左 | `gene_dating` 的 `gene_split` 与参考 CDS 查找在 standard 布局下失配 |

**同一函数内两种口径**（最能说明问题）：batch 收尾里 `:4380` 用硬编码 `"report"` 写汇总，`:4397` 用 `ph_dir('report')` 写 publication/sample_plot → 同一批产物落到两个不同目录。

### P1-3 ✅已实测（日志）：`clock` 在 TreeTime-RTT 失败后仍标记完成
`phylo_results/phylo_pipeline.log:569` `ERROR TreeTime-RTT failed: ValueError: ...` → `:607` 仍 `done=[..., 'clock']`。
RTT 在方法学上是**参考口径**（主判据是 DRT，`utils/temporal_signal.py:290` `passed = percentile >= 95.0`，与 `METHODS.md:25-26` 一致），所以不算结论错误；但汇总表 `Beta / R_squared / P_value` 会缺值，而 stage 声称"已完成"。

### P1-4 深读：退化 `merged.log` 不被向上传播
`utils/merge_results.py:219-221` 对 `<10` 行数据只 `log("❌ ...")`，而 `merge()` **无返回值**（函数在 `:309-310` 写 `run_status.json` 后结束）；调用方只判 `os.path.exists`（`phylo_pipeline.py:2795`、`:2582`）；`STAGE_CHECKPOINT_FILES['beast']` 与 `_beast_ran()`（`:3468-3472`）也只判"存在且非空"。
→ burnin 超界/链未跑完造成的退化 `merged.log` 仍会被登记为可用产物。仅"TMRCA 复用"分支有 `<10 行` 闸门（`:3493-3497`）。
（注：`merge_results` 内部有 ❌ 闸门，所以**不是**"完全静默"；问题在闸门结论没有传导到调用方与 checkpoint。）

### P1-5 深读：`_refresh_metadata_channels` 失败不记账
`phylo_pipeline.py:3147-3155` 只在成功时更新 `results["outputs"]`；失败路径仅由被调函数内部 `logger.warning`（`:2996`）。`process_virus` 不 `_skip` / `_fail` → 三通道拆分失败不会出现在 `Stages_*` 列里。

### P1-6 深读：必需外部工具缺失不阻断启动
`check_tools`（`:611-623`）只打印 `✗ NOT FOUND`。这与流水线"fail-soft"的整体策略一致，但意味着 batch 模式在缺 `mafft`/`iqtree` 时会一路跑到 stage 失败才发现（本机日志每次启动都打印 `mafft [必需]: ✗ NOT FOUND` 却继续）。

---

## 3. P2 — 死代码 / 重复实现 / 遗留物

- **不可达代码**：`phylo_pipeline.py:1249` 的 `run_mafft(...)`（见 P0-1）。
- **纯信息常量**：`STAGE_MODULE_DIR`（`:172-181`）注释自称"无代码引用"。
- **重复实现（漂移风险）**：
  - 小数年：`utils/decimal_year.py` 与 `utils/metadata_governance.decimal_year` 两份。✅已实测 `mode='mid'` 与统一实现在 5 个用例上**逐位一致**（`2010`/`2010-06`/`2010-06-15`/`2010-12-31`/`2012-02-29`），当前无差异，但属结构性隐患。
  - ID 规范化：`utils/seq_ids.py` 与 `utils/seq_clean.py::_norm_header`（P0-2 的根因）。
  - 同名文件：根目录与 `utils/` 各一份 `batch_draw_pymol.py`；`run_postprocess` 只会选 `utils.` 版（`:2896-2897`）。
- **未接入主流程（作者自己的审计工具，属预期）**：`orphan_audit.py`、`dependency_audit.py`、`patch_pypopart_haplotype.py`。✅已实测运行 `orphan_audit.py`：`活代码 83 | 管线可达 70 | 预期孤儿 10 | 真孤儿 0`；`dependency_audit.py`：项目内 import **无断链**。
- **嵌套 git 仓库**：`utils/RDP5_RBDP_Rgrapher/.git` —— 第三方仓库连 `.git` 一起入库，会干扰上层 `git status`/子模块语义，建议移除该 `.git`。
- **遗留目录**：`.ruff_cache/`（可直接删）；`cache/geo_cache.json` 是 `utils/geo_resolver.py:14` 声明的 L4 持久缓存，**属设计内，不要删**；`_scratch_tests/` 已被 `.gitignore:69` 覆盖。
- **无测试**：全目录无 `test_*.py`。

---

## 4. 文档与代码不一致（本流程的独特价值）

| # | 文档承诺 | 代码实际 | 判定 |
|---|---|---|---|
| 1 | `RUN_GUIDE.md:1` 「三入口 + **17 stage**」；`:373` 「stage 清单/数量（17）」 | `STAGE_ORDER`（`phylo_pipeline.py:132-138`）= **18** 项 | 文档自相矛盾（`:52` 写的是 18） |
| 2 | `RUN_GUIDE.md:68` beast checkpoint = `time/beast/beast1.xml` | `:234` = `time/beast/merged/merged.log`（2026-10-03 改） | 文档未同步 |
| 3 | `RUN_GUIDE.md:57,59` metadata/clean checkpoint = 「仅标记 `.done`」 | `:224-225` 另有产物回退 `data/metadata/metadata_std.csv` / `data/clean/clean.fasta` | 文档未同步 |
| 4 | `RUN_GUIDE.md:115` 「CLI 参数（**96** 个）」 | `main()` 内 `add_argument` 实测 **108** 次 | 数字过期 |
| 5 | `README.md:47`、`RUN_GUIDE.md:28` 「utils/ 依赖模块（**51** 个）」 | 实测 `utils/*.py` = **58** | 数字过期 |
| 6 | `README.md:51-52` 把 `popgen_analysis.py` / `build_haplo_outputs.py` 列在 `utils/` 下 | 两者都在**根目录**，`utils/` 下不存在 | 路径错误 |
| 7 | `RUN_GUIDE.md:10,369`、`README.md:23`、`phylo_pipeline.py:171` 均引用 `_consistency_check_20260916/verify_docs_vs_code.py`（「不符即非零退出」「硬断言守门」） | **该文件与目录均不存在**（全库 glob 无结果） | **校验机制已失效** —— 这正是 1-6 能长期存在的原因 |
| 8 | `STAGE_REFERENCE.md:310-312` align 处理 = `mafft --auto`、输出 = `phylogeny/mafft.aln.fasta` | `align` 不跑 MAFFT（P0-1） | **文档与代码正相反** |
| 9 | `RUN_GUIDE.md` §六 输出布局用 legacy 名（`data/`、`phylogeny/`、`time/`…） | 默认布局是 standard（`01_data/`、`02_phylogeny/`、`06_time/`…） | 表格未标注布局前提 |
| 10 | `--metadata_mode` help（`:3829-3833`）「mid=月中/**年中**」 | `utils/decimal_year.py:51-52` 只给年 → 取 **6 月 15 日**（实测 `2010 → 2010.4556`） | 命名与实现差约 16 天（0.456 vs 0.5） |

**已核对为「一致」的项（不是问题，列出以免误判）**：

- DRT 主判据：`utils/temporal_signal.py:290` `percentile >= 95.0` == `METHODS.md:25-26`。
- Mantel 置换数：`utils/virphy_bridge.py:1274/1934` 默认 `9999` == `METHODS.md:38`。
- BEAST 链数：`config.yaml` `beast.chains: 5` / `phylogeo.chains: 5` 已登记进 `CONFIG_SECTION_ARG`（`:376`、`:404`），实测 `config_orphan_keys(cfg) == []`。
- IQ-TREE 支持值默认：`IQTREE_BOOT_DEFAULT = IQTREE_ALRT_DEFAULT = 1000`（`:652-653`）== `RUN_GUIDE.md:145-146`。
- `STAGE_GROUPS` 八个别名（`:184-193`）== `RUN_GUIDE.md:78-89`。
- 日期口径分叉**是有意的**：`utils/virphy_bridge.py:714-751` 明确注释「年份/年月两级显式保留上游 VirPhyKit 原口径（年初/月初）」，仅服务 Treedater LTT 路径；其余调用方（`beast1_bridge` / `beast_bridge` / `phylogeo_bridge` / `import_export`）统一到 `utils/decimal_year`。**不是缺陷，无需改。**

---

## 5. 修复优先级建议

**第一档 — 改变结论（先修这三个）**

1. **P0-1 MAFFT 不可达**：把 `:1230` 的守卫限定为"仅当 `combined_fasta` 确实来自 `--virus`/`work_dir` 的 `--alignment`"（例如新增 `prep['provided_alignment']=True` 标志），batch 模式必须走 `run_mafft`。
   ⚠️ 这是**判定语义改动**，且修完会导致所有历史 batch 产物失效（checkpoint 指纹变化）。方向请老师定。
2. **P0-2 静默删序列**：`utils/seq_clean.py:359-361` 的"metadata 无此行"应只在 `drop_no_date or drop_no_location` 打开时才计入 `issues_this`；否则最多记 warning。同时把 `_norm_header` 与 `utils/seq_ids.base_sample_id` 收敛为一处。
3. **P0-3 batch 收尾输入**：`phylo_pipeline.py:4422` 改为传 `results["outputs"]["alignment"]`（或 `aln_file`），`:4424` 改为传 `results["outputs"]["tree"]`。

**第二档 — 部署/运维**

4. P0-4：给 `--virus` 分支补 `pipeline_exit_code` + `raise SystemExit`（照抄 `:4316-4319`）。
5. P1-2：四处硬编码改 `ph_dir(...)`；补一条"全库不得出现 `os.path.join(<root>, 'data'|'time'|...)` 字面量"的静态检查。
6. P1-1：`--strip-illegal` 要么修好（剥离成功则不计 issue），要么从 help 里删掉。
7. P1-4/P1-5：让 `merge()` 返回状态；`_refresh_metadata_channels` 失败时 `_skip(..., 'degraded')`。

**第三档 — 可维护性**

8. **恢复 `verify_docs_vs_code.py`**（或等价脚本）并把 §4 表格做成断言。没有它，文档漂移会继续累积。
9. 同步 §4 的 1-6、9-10 六处文档数字/路径。
10. 删 `utils/RDP5_RB_Rgrapher/.git` 与 `.ruff_cache/`；`cache/geo_cache.json` 保留。

---

## 6. 值得保留的设计（建议加注释保护，不要被后续重构抹掉）

1. **`validate_stage_graph()`**（`:422-447`）—— 把"执行顺序 = 代码块物理顺序"这个隐形契约变成启动期硬断言。这是本文件里最有价值的守卫。
2. **`_inputs_fingerprint` + `MISSING` 稳定值**（`:251-274`）—— 输入指纹让"换参数不重跑"变成不可能；`MISSING` 作为稳定值避免文件不存在导致每次重跑。
3. **`_skip` 五态 + `Incomplete` 列**（`:3091-3108`、`:3602-3603`）—— 把"整块静默缺席"变成可检测、可汇总、可反映到退出码。日志里 `geo_paths(blocked)` 的实际表现证明了它的价值。
4. **`config_orphan_keys()`**（`:509-535`）—— 针对"键没登记就悄悄丢"这一类静默失效的通用解法，值得推广到其它管线。
5. **`_beast_ran()` 二次判定**（`:3468-3472`）—— "XML 存在 ≠ 跑过 MCMC"这个区分做得很准。
6. **`seq_clean` 的 fail-loud 护栏**（`utils/seq_clean.py:295-302`）—— 宁可拒绝执行也不整批误删，是删数据场景的正确姿态。
7. **`--explicit_cli_dests`**（`:538-565`）—— 让"CLI > config > 代码默认"真正可执行，且实现不依赖 deepcopy parser。
8. **`merge_results` 的近空产物闸门**（`:219-221`、`:285-288`）—— 对"文件存在但内容退化"这类最隐蔽的假成功有针对性。

---

## 附：本次审查的副作用与自证

- **未修改任何源文件**。`git status --porcelain -- virome_phylo_pipeline` 为空（无 tracked 改动）。
- 运行 `--dry` 时 `setup_logger` 的固有行为向 `phylo_results/phylo_pipeline.log` 追加了 2 条 dry 记录（12:32:33、12:32:36）—— 日志本身是运行产物，未回滚。
- 审查过程中子代理在 `_scratch_tests/rtt_test/` 造的测试夹具**已由其自行清理**；探针脚本全部落在系统临时目录并已删除。
- 未触碰父仓库其它被修改文件（`README.md`、`doc/*` 等为作者本人的工作区改动）。

## 附：本轮纠正的误报

- 子代理报告「`utils/decimal_year` 把纯年份当 **7 月 1 日**」——**我上一版据此写错，实测为 6 月 15 日**（`to_decimal_year('2010') = 2010.4556`，`utils/decimal_year.py:51-52` 明写"只给年 → 取 6 月 15 日"）。该条已改写为 §4 第 10 行的"命名与实现差约 16 天"。
- 子代理报告「burnin 超界导致 `merged.log` 近空且**完全静默**」——实测 `utils/merge_results.py:219-221` 有显式 ❌ 闸门，**不是静默**。已改写为 P1-4「闸门结论未传导到调用方与 checkpoint」。
