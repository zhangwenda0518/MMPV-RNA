# MMPV-RNA 系统发育与进化分析 · 完整流程

> **📌 文档定位：本文件是「模块归属 / 数据流」的权威**（概念层）—— 主线阶段与散置进化分析模块
> （RDP5 重组、capheine 正选择、pypopart 群体遗传、分基因定年）各自坐在数据流的哪个位置，
> 以及当前是否已挂入总管线 `phylo_pipeline.py`。
> 「怎么跑」见 `RUN_GUIDE.md`；stage 逐项细节见 `STAGE_REFERENCE.md`；论文表述见 `METHODS.md`；
> 文档总索引见 `README.md`。机器校验：`_consistency_check_20260916/verify_docs_vs_code.py`。
>
> **同步状态（2026-09-16）**：路径与开关名已按代码校正（模块化目录 `data/`、`select/capheine/`、
> `--rdp5_mask`、三通道含 `host.csv`）。改动前备份 `_docconv_20260916/backup/PIPELINE_FLOW.md`。

---

## 一、核心结论（先看这个）

整条 phylo 分析存在**两条数据流分叉**，而不是一条直线：

- **A 线 · 全长序列流**：全长组装序列 → MAFFT 全长比对 → 〔群体遗传 pypopart · 重组 RDP5 分叉〕 → IQ-TREE 建树 → 定年/地理/可视化
- **B 线 · 分基因流**：GBK 注释 → gbk_extractor 提 CDS → capheine（CAWLign 密码子比对 + 基因树）→ 〔正选择 HyPhy/selection_suite · 分基因定年 gene_partition_dating 分叉〕

两条线的**汇合点**只有一个：capheine 产出的分基因比对，既喂给正选择，又喂给分基因定年。其余各自独立。

---

## 二、完整数据流图

```
════════════════════════════════════════════════════════════════════
【数据准备】                                                    (主线)
────────────────────────────────────────────────────────────────────
 Phase 1  上游全长组装序列 + 元数据收集
 Phase 2  公共序列收集 + 采样日期匹配
 Phase 2b NCBI 参考补充 (SeqHarvester)
════════════════════════════════════════════════════════════════════

            ┌────────────────────────────────────────────┐
            │         MAFFT 全长比对 (Phase 3)           │
            └────────────────────────────────────────────┘
                          │                    │
                          │                    │
        ┌─────────────────┴──────┐             │
        │                        │             │
        ▼                        ▼             │
 ┌──────────────┐      ┌──────────────────┐    │
 │  Phase P     │      │  Phase R         │    │
 │  群体遗传学  │      │  重组检测 RDP5   │    │
 │  pypopart    │      │  (7 方法)        │    │
 │ π/θ/Tajima/Fst│      │  └─ mask 重组区  │    │
 └──────────────┘      └────────┬─────────┘    │
        (输出报告)               │ 干净比对      │
                                ▼              │
                    ┌──────────────────────┐   │
                    │  IQ-TREE 建树        │   │
                    │  (Phase 3 收尾)     │◄──┘
                    └──────────────────────┘
                          │
      ┌─────────┬─────────┼─────────────┬─────────────┐
      ▼         ▼         ▼             ▼             ▼
 Phase 0a   Phase 0b   Phase 5      Phase 6       Phase 7
 DRT+RTT   BETS       BEAST 定年   系统地理      可视化
 (TreeTime)(path samp) (全长)     (CTMC+BSSVS)  (SpreaD3+HTML)

════════════════════════════════════════════════════════════════════
【分基因流 · 独立于 A 线】                                    (旁支)
────────────────────────────────────────────────────────────────────
 GBK 注释 (来自 analysis 或独立来源)
        │
        ▼
 gbk_extractor.py  ── 提取 CDS (每条=一个基因)
        │
        ▼
 ┌──────────────────────────────────────────────┐
 │  Phase S · capheine_pipeline.py              │
 │  CAWLign 参考锚定密码子比对 → IQ-TREE 基因树  │
 │  → HyPhy FEL/MEME/PRIME/BUSTED/CONTRASTFEL/  │
 │        RELAX → DRHIP → MultiQC                │
 └──────────────────────────────────────────────┘
        │
        ├──────────────────────────────┐
        ▼                              ▼
 ┌────────────────────┐      ┌──────────────────────────┐
 │ 正选择下游         │      │  Phase G · 分基因定年    │
 │ selection_suite    │      │  gene_partition_dating.py │
 │  (FUBAR/codeml/    │      │  --genes-dir 吃 capheine  │
 │   dnasp/timetree/  │      │  产物 → BEAST 分基因 MRCA │
 │   network/structure│      └──────────────────────────┘
 └────────────────────┘
```

---

## 三、阶段总表

| 阶段 | 名称 | 输入 | 核心工具 | 关键输出 | 挂载状态 / 对应 stage |
|---|---|---|---|---|---|
| Phase 1 | 上游全长组装序列+元数据收集 | 06_extraction | data_collector | combined.fasta + 三通道表（dates/sample_metadata/host） | ✅ `prep` |
| Phase 2 | 公共序列收集 + 采样日期匹配 | Phase 1 | SeqHarvester + data_collector | ncbi_harvest/ 公共序列 + 日期匹配 | ✅ `prep`(日期) + `online`(公共序列) ⚠️ **公共序列仅入报告，不进分析**（见下方警示） |
| Phase 2b | NCBI 参考补充 | Phase 2 | efetch + SeqHarvester | `data/ncbi_ref/` + 参考序列元数据 | ✅ `online` |
| **Phase P** | **群体遗传学** | **MAFFT 全长比对** | **pypopart** | **π/S/θ/Hd/Tajima/Fst/MJN** | ✅ 主线 (`popgen`) |
| **Phase R** | **重组检测** | **MAFFT 全长比对 + genes.csv** | **RDP5 (7 方法)** | **`recomb/rdp5/<virus>.csv` + 干净比对** | ✅ 开关（`--rdp5`） |
| Phase 3 | 比对 + 建树 | Phase 2 输出 | MAFFT → IQ-TREE | mafft.aln + iqtree.treefile | ✅ 主线 |
| Phase 0a | DRT 时间信号 | 比对 + dates | TreeTime | DRT p 值 | ✅ 开关 |
| Phase 0b | BETS 时间信号 | 比对 + meta | path_sampling | log BF | ⛔ 已归档止损（开关 `--bets` 仍在，但结论不纳入结果） |
| Phase 5 | BEAST 全长定年 | 比对 + dates | BEAST **1.x**（`time/beast/beast1.xml`） | 全长 MRCA | ✅ 开关（`--beast`） |
| **Phase G** | **分基因定年** | **capheine 产物** | **gene_partition_dating --genes-dir** | **分基因 MRCA** | ✅ 开关（`--gene_dating`） |
| **Phase S** | **正选择** | **GBK CDS** | **capheine → HyPhy** | **FEL/MEME/BUSTED/DRHIP** | ✅ 主线 (`capheine`) |
| Phase 6 | 系统地理 | 比对 + meta + tree | discrete CTMC + BSSVS | SpreaD3 地图 | ✅ 开关 |
| Phase 7 | 可视化 + 报告 | 各阶段产物 | SpreaD3 + report_builder | 汇总 HTML | ✅ 主线 |

> **挂载状态说明**：✅ 主线 = 已由 `phylo_pipeline.py --stage` 编排；✅ 开关 = 已挂载但需 flag 触发（`--beast`/`--phylogeo`/`--check_temporal`/`--rdp5`/`--gene_dating`，可组合，`--beast`+`--gene_dating` 同开即全长+分基因都定年）；⚠️ 独立 CLI = 未进总管线，需 `python -m xxx` 手工调用。

---

## 四、数据流依赖（关键，决定模块位置）

### 4.1 谁吃 MAFFT 全长比对

`recombination_analysis.py` 和 `popgen_analysis.py` 的**共同上游是 MAFFT 全长比对**。因此它们在数据流上并列，都位于「比对之后」：

- **Phase P（pypopart）**：`PopGenAnalyzer(mafft_aln, metadata_csv)`，不改变比对，只产出群体遗传统计与 MJN 网络。
- **Phase R（RDP5）**：`--fasta mafft.aln.fasta --genes genes.csv`，检出重组事件后**默认策略是删除高置信重组子 + 重新 MAFFT**（`recomb/rdp5/rdp5_dropped.realn.fasta`）；加 `--rdp5_mask` 则改为重组区置 N（`masked/masked_N.fasta`）。两种产物**都会覆盖 `align` 的输出喂给 IQ-TREE**。

> 这就是「重组必须在建树前」的准确含义：RDP5 插在 **MAFFT 之后、IQ-TREE 之前**，而不是在 MAFFT 之前。

### 4.2 谁吃 GBK CDS（分基因流）

`capheine_pipeline.py` 不经过全长 MAFFT，它从 **GBK 注释提取的 CDS** 出发，内部用 **CAWLign 参考锚定**完成密码子感知比对 + IQ-TREE 基因树。这条线独立于 A 线。

### 4.3 两条线的汇合点

capheine 产出的**分基因比对目录**是唯一汇合点：

```
capheine 产物目录
   ├─> selection_suite.py（FUBAR/codeml/dnasp/timetree/network/structure）
   └─> gene_partition_dating.py --genes-dir（分基因 BEAST 定年）
```

因此 **Phase S（capheine）是 Phase G（分基因定年）的前置依赖**，也是正选择下游（selection_suite）的前置依赖。

---

## 五、迁移后待统一的接口点

1. **capheine 产物目录约定**：实际路径为 `<virus>/select/capheine/`（**模块化目录**，2026-09-16 校正；旧文档写的 `<virus>/capheine/` 已作废）。`selection_suite.py --capheine-dir <virus>/select/capheine/`；`gene_partition_dating.py --genes-dir` 指向其内部基因比对目录（`cawlign/`）并自动适配命名。
2. **gbk_extractor 的输入来源**：迁移后 phylo 需明确 GBK 注释从哪里取（独立提供 vs 从 analysis 的 `virus-annotations/` 拷贝），避免 phylo 回读 analysis 目录。
3. **RDP5 与 mask 的插入**：已完成。实际位置是 **`align` 之后、`tree` 之前**（stage `rdp5`，`STAGE_DEPS` 里 `rdp5` 依赖 `align`、`tree` 消费 rdp5 产物）。默认删重组子+重比对，`--rdp5_mask` 改为置 N，`--mask_min_methods`（默认 3）控制"高置信"门槛。
4. **capheine 版本确认**：迁移前实测本地与服务器四个文件字节完全一致（capheine 21431 / gbk_extractor 5586 / visual_codon_miner 12654 / batch_draw_pymol 21589），无漂移；archive 中 21440 为历史记录值，非主目录现状。

---

## 六、Phase 编号与现有代码对应关系

| 本文档编号 | `phylo_pipeline.py` 现状 |
|---|---|
| Phase 1 + Phase 2（日期匹配） | stage `prep`（`data_collector`：load_sample_metadata / scan_extraction_dir / prepare_virus_inputs） |
| Phase 2（公共序列收集） | stage `online`（`seq_harvester` 按物种名收集 NCBI 公共序列 + metadata） |
| Phase 2b | stage `online`（`Entrez.efetch` 单参考 + SeqHarvester + seq_grouper 分组图） |
| Phase 3 | v2 已拆分为 stage：`align`（MAFFT）+ `rdp5`（插入点）+ `tree`（IQ-TREE）+ `clock`（TreeTime/TreeDater/DRT）。v1 的单一 `run_phase()` 已不存在 |
| Phase P / S | 已挂主线 stage（`popgen` / `capheine`；实现分别在根目录 `popgen_analysis.py`、`capheine_pipeline.py`） |
| Phase R / G | 已挂主线开关（`--rdp5` / `--gene_dating`） |
| Phase 0a / 0b | `--check_temporal` / `--bets` 开关 |
| Phase 5 / 6 / 7 | 原 Phase 5 / 6 / 7（`run_phase_beast` / phylogeo / viz） |

> ✅ 编号冲突已定夺（2026-08-27）：数据准备第三段定名 **Phase 2b**（对齐代码 `run_phase_data_online`），**Phase 3** 专指「比对+建树+初步定年」。数据准备三段 = Phase 1 / Phase 2 / Phase 2b，无编号冲突。

---

### ⚠️ 两条容易踩的接线事实（2026-09-18 逐接口核查后补记）

**① `ncbi_harvest/` 的公共序列目前"只进报告，不进分析"。**
`online` 里 SeqHarvester 采到的 `ncbi_sequences.fasta` + `ncbi_metadata.csv` + `seqgroup_*.csv`
+ 分组饼图，全库**只有 `report_builder` 引用**（MODULES 的 data 组 summaries/plots glob）——
**不进 `combined.fasta`、不进比对/建树/定年**。`prep` 只从 `06_extraction` 收集序列
（参考序列另走 `--variants_dir`）。
⇒ 别指望"跑一遍 `online` 就自动扩了样本"。要真扩样本：把
`<work>/data/ncbi_harvest/ncbi_sequences.fasta` 并入 `06_extraction/<病毒目录>/` 后重跑 `prep`，
并确认这些样本在 metadata 表里（否则会被 `clean` 的日期/地理门槛剔除）。
（对比：`ncbi_ref/` 的单参考 `*.gbk` **是通的** —— `capheine` 从 `data/ncbi_ref/*.gbk` 读取。）

**② 真实执行顺序 = `STAGE_ORDER`（这条曾不成立，2026-09-18 已修）。**
`process_virus` 里每个 stage 是一段 `if '<stage>' in stages:` 硬编码块，**没有
`for stage in stages` 循环** → 执行顺序取决于**块在文件里的物理位置**，而
`resolve_stages` 的依赖闭包只保证"被选中"，**不保证先后**。
历史事故：`geo_paths` 声明依赖 `phylogeo`（吃它的 `merged/mcc.tree`）却排在其**之前**
→ 首次全新跑时 mcc.tree 尚不存在 → 整段**静默跳过**，且 `_skip` 不计失败 → **仍报成功**。
现已把物理块序对齐到 `STAGE_ORDER`，并加 `validate_stage_graph()` 在启动时拒绝
「依赖排在使用者之后」的改动。**改 `STAGE_ORDER` 时必须同步移动 `process_virus` 里对应的
代码块**，否则顺序不会变。
