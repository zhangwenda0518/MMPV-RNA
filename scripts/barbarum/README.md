# scripts/barbarum — 宁夏枸杞 (Lycium barbarum) 专项脚本

单物种已知病毒分析的复刻/重绘工具。**不属于 MMPV 主流水线**，是针对
`03_known_virus/known_virus_barbarum` 的项目级脚本。

## 为什么需要它们

`known_virus_barbarum` 最初是`known_virus_all_v2`（跨物种 + 宏基因组混合队列，
1400 样本）的软链镜像。整合 42 个缺失样本后已改为独立目录，但 01/02 的
若干产物仍是 v2 的 symlink，检出率口径与单物种 barbarum 不符。

## 权威数据来源

```
02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly/   <- 分母的权威来源
```
该目录 446 个条目 = 445 个 run + 1 个 `logs`。**所有检出率的分母 = 445。**

`01_detection/summary/all_viruses.best.summary.tsv` = 536 条记录 / 262 个run /26 种病毒。

| 口径 | 数值 |
|---|---|
| 分母（队列总数） | 445 |
| 检出样本 | 262（58.9%） |
| 未检出 | 183（41.1%） |
| 单一感染 | 96 |
| 合并感染（≥2） | 166 |
| 检出病毒种类 | 26 |

高流行率病毒（≥5%）：GCVA 179例 40.2% / PSTVd 174 例 39.1% / TCDVd 91 例 20.4%。

## regen_barbarum_detection_report.py

重生成 `01_detection/summary/analysis_report.txt`。

复用 `virome_analysis_pipeline/batch_virus_depth.py::generate_report_txt`
（第 755–773 行）的报告结构与统计口径（平均 CPM / FPKM / 覆盖度 /
Poisson Ratio / 群体检出率），但分母改为 barbarum 权威run 数，并额外补充
队列构成、共感染结构、覆盖度质量分级。

```bash
python scripts/barbarum/regen_barbarum_detection_report.py \
  -v.../known_virus_barbarum/01_detection/summary/all_viruses.best.summary.tsv \
  -d .../known_virus_barbarum/01_detection/summary \
  --denominator-sources .../02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly \
  [--dry-run]
```

**安全护栏**：写盘前若目标是 symlink 会先 `unlink`，绝不写穿到 v2。

## redraw_barbarum_detection_figs.py

重绘 `Fig1_study_design.{png,pdf}` 与 `Fig2_prevalence.{png,pdf}`。

原图是 v2 的跨物种版本（标题写死 `Total 1,156` / `820 Lycium RNA`，
含 `DNA metagenome 313`），且**原生成脚本已不存在**——在 MMPV-RNA 全库 +
git 全历史 + bash history + `/tmp` 穷尽搜索 `Fig1_study_design` 零命中。
`panvirome_suite.R` / `panvirome_landscape.py`（git 已删）出的是
`Plot1_..Plot6_` 病毒级图，不是本图系。

新图内容：
- **Fig1** (a) 单物种研究流程  (b) 队列构成与检出状态（445 / 262 / 183）
- **Fig2** (a) 病毒流行率降序（RNA virus / Viroid 配色区分）
         (b) Top4 病毒的组织嗜性  (c) 地理分布

```bash
export LD_LIBRARY_PATH=/home/zhangwenda/mambaforge/lib:$LD_LIBRARY_PATH
python scripts/barbarum/redraw_barbarum_detection_figs.py \
  -v .../01_detection/summary/all_viruses.best.summary.tsv \
  -m .../01_data/sra_rna.data4/global_metadata/Global_Metadata_Final_with_Coordinates.tsv \
  -d .../01_detection/summary \
  --assembly-root .../02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly
```

matplotlib 3.9.4（与原图同版本），300 dpi。**配色不写在本文件里**——按显式
路径加载 `virome_discovery_pipeline/utils/sci_palette.py`（见下节）。

### 画布宽 = 通栏宽（2026-10-12 修订，重要）

**印刷字号 ≈ 源字号 × (面板宽 / 源画布宽)**。旧版 Fig1 画布 13.0 in、Fig2 16.5 in，
排进通栏（`COL_W_1 = 9.32 in`）缩放只有 **0.717 / 0.565** —— Fig2 的 8.6 pt
种名印出来只有 **4.9 pt**，不可读。

现改为：

| 图 | 画布 | 布局 | PDF 实测 | 通栏缩放 |
|----|------|------|----------|---------|
| Fig1 | `9.32 × 5.60 in` | 1×2 | 9.32 × 5.60 in | **1.000** |
| Fig2 | `9.32 × 8.90 in` | 2 行：(a) 通栏 / (b)(c) 并排 | 9.33 × 8.46 in | **0.999** |

字号按**印刷尺寸**直接设定：面板标题 9.5 / 框标题 8.5 / 副标题 6.5 / 刻度 7.5 pt。
验收手段：`pypdf` 读 MediaBox 算 `9.32 / width`，要求落在 **0.9–1.1**。

### 实现要点（踩坑记录）

1. **`tight_layout()` 不可用于 Fig1**——面板 (a) 用 `axis("off")` +
   数据坐标 0~1 画流程框，与 tight_layout 冲突，会抛
   `margins cannot be made made large enough`，并使 `bbox_inches="tight"`
   算出 80397 px 的超宽图（`2^16` 上限）。改用固定 `subplots_adjust`。
2. **文字高度必须按字号现算到数据坐标**，不能用固定偏移：
   ```python
   pt_per_unit  = FIG1_H * (FIG1_TOP - FIG1_BOTTOM) * 72.0 / (ylim[1] - ylim[0])
   title_line_h = TITLE_PT * 1.20 / pt_per_unit * (title.count("\n") + 1)   # ⚠ 要乘行数
   sub_line_h   = SUB_PT   * LS   / pt_per_unit
   ```
   配两条 `assert` 硬守卫：① 副标题底不得低于框底；② `PAD_TOP + content_h < BOX_H`
   （内容整体塞得进框）。溢出即报错、不产图。
3. **多行标题是隐蔽坑**：`Full-length assembly\n& variant analysis` 只按一行算
   → 标题第二行压副标题。现已压成一行 `Assembly & variant calling`。
4. **图例不要复用 barh 的容器句柄**（会继承 `height=0.72` 的奇怪比例）；
   Fig2(b) 图例必须放**左上**（右上被 `flower` 列 PSTVd ≈100 % 柱顶穿），
   且用短代号（长种名会把图例撑出轴框）。
5. **Fig2 画布改 2 行后 `bbox_inches="tight"` 会外扩**——(a) 的长种名超出
   `add_axes` 左边界会把总宽顶到 9.78 in。把 (a) 左边界从 0.215 调到 **0.265**
   即收回 9.33 in。

## make_discovery_sci_figs.py

**discovery 侧第一张静态出版级主图** `Fig_GojiVirome_Discovery`（**v18：7 面板 A/D/I/J/K/L/M**
+ 组合图 + 7 个单面板导出 + `fig_manifest.json`）。此前 discovery 侧
**只有交互式 HTML 报告，没有任何静态主图**。

> **v18 补全 A 面板链条（2026-10-13 第 12 轮，用户指令「补充 filter 后的病毒的数量」）**：
> 面板 A 的漏斗由 **10 条**补成 **13 环节全链**（全部纳入 + 逐项标注来源）：
>
> | # | 环节 | 值 | 来源 | 口径 |
> |---|---|---|---|---|
> | 1 | Assembled contigs | 14,387,168 | `01_Assembly` contigs | all-host |
> | 2 | Viral sequences (10 tools) | 556,773 | `02_Identification` viral_seqs | all-host |
> | 3 | **Filtered candidates** | **94,826** | `02b_Filter/<样本>/uniprot/all_candidates_uniprot.fasta` 计数汇总 | all-host |
> | 4 | **COBRA queries (extended)** | **35,348** | `cobra_summary.tsv` `Total_Queries` 汇总 | all-host |
> | 5 | Merged / re-assembled | 80,637 | `03b_Merge` combined | all-host |
> | 6 | vOTU clusters | 37,803 | `cluster_pipeline_reduction.tsv` | all-host |
> | 7 | Taxonomically classified | 20,935 | `05_Taxonomy` classified | all-host |
> | 8 | Host-assigned | 20,892 | `06_HostPrediction` classified | all-host |
> | 9 | Plant-host (subset) | 1,063 | 同上 `Plant=` 行 | plant |
> | 10 | **Plant HQ vOTUs** | **63** | `checkv_summary.tsv` **Plant 行** `HQ` 列（33+30） | plant |
> | 11 | Rescued (VSI, plant subset) | 276 | `08_Rescue` rescued | plant |
> | 12 | **CheckV high-confidence** | **45** | `keep_summary.tsv` `checkv_confidence=="high"` | plant |
> | 13 | **Final HQ plant viruses** | **291** | `10_Virome_Analysis` `HQ viruses:` | plant |
>
> ⚠ **口径分两块，不可读成单一漏斗**：前 8 条 = 全宿主（all hosts）；后 5 条 = 植物子集
> （是 `Host-assigned` 的**切片**，不是它的下一级）。图上用**右侧 scope 括号**
> （`all hosts` / `plant subset`，轴分数坐标）+ 括注标清，避免被读成"20,892 → 1,063 掉 95 %"。
>
> ⚠ **新增数据源 + 三个坑**：
> ① `Filtered candidates` **不能**读 `filter_summary.tsv`（本项目里它是 **396 行全 0 的空表**：
> `ident_stats.py` 按 `uniprot_filter_output*` 目录名找产物，与实际 `uniprot/` 不匹配 → 静默全 0）
> → 改从各样本 `uniprot/all_candidates_uniprot.fasta` 的 `>` 计数汇总，落盘
> `10_Reports/filter_per_sample.tsv`（395 行，`Uniprot_Filtered` 汇总 = 94,826，权威可复核）。
> ② `Plant HQ` 用 `checkv_summary.tsv` **Plant 行**（分母 **1,090** ≠ 本图其它处的 **1,063**，
> 属已知 **R6** 两套植物口径分歧），用户 2026-10-13 明确选定此口径（不取 `07_CheckV` 全宿主 200）。
> ③ 服务器 `02b_Filter` 下的 fasta 是**断裂软链**（指向已移走的数据盘）→ 必须从其**子目录实体文件**取数。
>
> ⚠ **新增配色键**：`SEMANTIC["stage"]["filter"]` = `#999999`（中性灰，表达"它是 identification
> 的过滤产物"，而非又一条独立口径）。**未在 `SEMANTIC` 里的键会静默退化为 `#999999`**。
>
> ⚠ **注记一律用轴分数坐标**（`ax.transAxes`）：早前按 `max(vals)×N` 估数据坐标，在 log 轴上
> 要反复反解 xlim，实测把注记顶出画布（体检报 `TEXT-OFF-CANVAS`）；改轴分数后与数据范围无关。

> **v13 重构（2026-10-12）**：原 12 面板 **12/12 全是条形图**，图型编码单一，
> 且存在硬重复。重构为 **9 面板 / 4 种图型**：
>
> | 新 | 内容 | 图型 | 来源 |
> |---|---|---|---|
> | A | 发现管线漏斗（10 条） | 条形 (log) | 保留 |
> | B | 鉴定工具贡献 | 条形 (log) | 保留 |
> | C | vOTU 簇大小分布 | 条形 (log) | 保留 |
> | **D** | **分类谱系流向** | **桑基（6 层）** | 新（替换原条形 D） |
> | E | 新病毒分层 | 条形 | 原 F |
> | **F** | **宿主组：总数 vs HQ** | **成对条形** | 原 G+H **合并** |
> | G | 植物病毒属组成 | 条形 | 原 I |
> | H | 终判类别 | 条形 | 原 J |
> | **I** | **基因组型 + 类病毒** | **成对条形/并排** | 原 K+L **合并** |
>
> **下线**：原 E（分类填充率，7 条全 80–100 %，仅 Genus 69.7 % 有信息）→ 补充材料。
> **合并判据**：原 G 与 H 是**同一组 9 个宿主类目、同一单位**的两张横条（HQ ⊂ total）
> → 硬重复；原 D 被原生生物巨型 DNA 病毒占满，对植物病毒论文是干扰项。

> **v14 借参考图补 2 面板（2026-10-12 第 8 轮）**：用户给 `gabrielvpina/dataViz`
> （Theobroma cacao 病毒组论文）的 `data_report.R` 做**版式参考**。原 9 面板**全缺
> "证据质量分布"**（只有计数类条形图，看不出每个判定有多可靠）→ 借参考图 C / F：
>
> | 新 | 内容 | 图型 | 参考图对应 |
> |---|---|---|---|
> | **J** | nt identity × coverage 散点（+ 顶部/右侧边缘箱线，着色 = verdict） | 散点 + 边缘箱线 | 参考图 **C**（`geom_point` + `ggMarginal(type="boxplot")`） |
> | **K** | 按科的 nt identity 箱线（`n ≥ 3`，按中位数排序） | 箱线 | 参考图 **F**（`geom_boxplot` + `coord_flip`） |
>
> 数据源 `rescue_evidence_scored.tsv`（`nt_pident` / `nt_qcovs` / `tax_family` / `verdict`）。
> **不能**改用 `All_plant.viruses_info.tsv`：它的 `confidence` 是**证据分**（1/0.88/0.75…），
> 不是 % identity，画成"identity 散点"是张冠李戴。

> **v15 复现参考论文 Fig2 版式，再补 2 面板（2026-10-12 第 9 轮）**：用户给
> `Viruses 2025, 17, 624`（*Theobroma cacao* Virome，PMC12115555）的 Fig2（元病毒组
> **6 面板 A–F**）做**结果对齐**参考。逐面板比对后发现仍缺两块"证据质量分布"，
> 遂再补 2 面板，主图升为 **13 面板 A–M**：
>
> | 新 | 内容 | 图型 | 论文对应 |
> |---|---|---|---|
> | **L** | 保守域完整性三分（Complete / Fragment / No conserved domain） | 饼图 | 论文 **Fig2-E**（`complete` / `fragmented` / `absent`） |
> | **M** | 比对靶标（HMM 命中的最佳病毒域按功能归类） | 横向条 | 论文 **Fig2-B**（病毒 contig 的比对靶标） |
>
> **论文 Fig2 六面板对齐总表**（写入 `reports/DISCOVERY_SCI_FIGS_20261012.md` §6.2c）：
>
> | 论文 | 我们的对应 | 状态 |
> |---|---|---|
> | A（contig → Meta/Host/Viral 三分） | 面板 A（五级漏斗） | ⚠ 语义近似（归属 vs 阶段），**不改动、记差异** |
> | B（比对靶标） | **面板 M（本轮新增）** | ✅ 已补 |
> | C（identity × coverage） | 面板 J | ✅ |
> | D（family × 遗传物质，known/novel） | 面板 I 左半 + G 着色 | ⚠ 部分对齐（**缺 known/novel 维**，无法合成） |
> | E（保守域完整性） | **面板 L（本轮新增）** | ✅ 已补 |
> | F（科内 nt 相似度） | 面板 K | ✅ |
>
> **新增数据源两条（用错 = 静默错误）**：
> ① **Fig2-E 结构域完整性只能用 `cdd_evidence_report.tsv` 的 `evidence` 列**
>    （`PASS_CORE`/`PASS_CORE_FAM`→Complete、`PASS_VIRAL`→Fragment、其余→No conserved domain）；
>    **不能**用 `checkv_summary.tsv`（那是**基因组完整性**，不是保守域完整性，语义完全不同）。
> ② **Fig2-B 比对靶标只能从 `hmm_ct3_hits.tsv` 的 `ct3_best_viral_desc` 按关键词派生**（无现成功能注释列）。
>
> **v15 版式参数（已废）**：`figsize=(8.55, 14.4)`、`GridSpec(7, 3, hspace=0.80, wspace=0.62,
> height_ratios=[1.00,1.25,1.05,1.30,1.15,1.05,0.80], left=0.085, right=0.975,
> top=0.945, bottom=0.038)`、`suptitle(y=0.983)`。面板 **M 跨 3 列**（末行独占），
> L 放在第 5 行第 3 列（半栏位）。

> **v16 删 4 面板（2026-10-13 第 10 轮，用户指令「EFGH都没啥意义删了呢」）**：
> 从 13 面板降为 **9 面板**，画布 14.4 → **10.6 in**，`GridSpec` 7 行 → **5 行**。
>
> | 下线 | 原内容 | 下线理由 |
> |---|---|---|
> | **E** | Novelty tier（新颖性分层） | 与 A 漏斗的 novelty 分段**语义重复**，且只 3–4 条 |
> | **F** | Host group: predicted vs HQ | HQ 是总数的**切片**（HQ ⊂ total），且与 I / H 的宿主口径交叉 |
> | **G** | Plant-associated virus genera | 与 **D（植物桑基）**+ M 的家族分布部分重叠 |
> | **H** | Final verdict | 与 I 的基因组型 / 类病毒构成**语义接近** |
>
> ⚠ **字母标识保持稳定（跳号 A/B/C/D/I/J/K/L/M），不重排**：文档 §6.2c 对照表、
> `fig_manifest.json` 的 key、脚本内 `WIDE` 集合与 `panel_letter` 分支**全部锚定字母**
> → 删减不重排，避免牵动其它逻辑与既有引用。四个函数体保留（加「已下线」docstring）。
>
> **v16 版式参数（已废）**：`figsize=(8.55, 10.6)`、`GridSpec(5, 3, hspace=0.80, wspace=0.62,
> height_ratios=[1.00,1.25,1.15,1.05,0.80], left=0.085, right=0.975,
> top=0.945, bottom=0.038)`、`suptitle(y=0.983)`。行内配伍：A(跨 2)+B ｜ C+D(跨 2)
> ｜ I(跨 2)+J ｜ K(跨 2)+L ｜ M(跨 3，末行独占)。
>
> ⚠ **`build_panels()` 只写不删** → 面板下线后旧 `panel_X.{png,pdf}` 会原地留存，
> 极易被误读为「仍在用」。v16 起加 **`prune_stale_panels(outdir, keep)`**：每次
> `--panels-only` 后自动删除不在 `PANELS` 里的 `panel_*.{png,pdf}`（仅动该形态文件，
> 不碰 `phylo_*`）。首跑实测删掉 E/F/G/H 共 8 个陈旧文件；v17 再删 B/C 共 4 个。

> **v17 再删 2 面板（2026-10-13 第 11 轮，用户指令「BC 也删了呢，没太大用」）**：
> 从 9 面板降为 **7 面板**，画布 10.6 → **9.4 in**，`GridSpec` 5 行 → **4 行**。
>
> | 下线 | 原内容 | 下线理由 |
> |---|---|---|
> | **B** | Identification tools（多工具贡献条形） | 工具贡献是**方法学过程量**，与 A 漏斗的 `Viral sequences (10 tools)` 同源；读者已由 A 得知规模 |
> | **C** | vOTU cluster size distribution | 簇统计属**中间产物 QC**（singleton 占比高是宏病毒组常规），不承载结论 |
>
> **v17 版式参数**：`figsize=(8.55, 9.4)`、`GridSpec(4, 3, hspace=0.72, wspace=0.62,
> height_ratios=[1.05,1.30,1.15,0.90], left=0.085, right=0.975,
> top=0.945, bottom=0.042)`、`suptitle(y=0.983)`。
> 行内配伍：**A(跨 2)+J** ｜ **D(跨 2)+L** ｜ **I(跨 2)+K** ｜ **M(跨 3，末行独占)**。
> 字母跳号 **A/D/I/J/K/L/M**（继续不重排）。
>
> **v18（2026-10-13）**：面板 A 漏斗补 5 项 → 13 环节（组装→…→final HQ）。
> **v19（2026-10-13）**：D 升**通栏**（跨 3 列，7.61 in），M 跨 2 列 + I 占 1 列同处末行；
> 画布 (8.55, 12.6)，`GridSpec(5,3, hspace=0.62, wspace=0.62)`，
> 行权 `[1.05, 2.10, 1.25, 0.95, 0.02]` → 面板 A/J/I/K/L/M 全部挤到很小
> （实测 A 高仅 ~1.6 in、J 宽仅 1.8 in）。
>
> ⭐⭐⭐ **v20（2026-10-13 晚，用户指令「A图展开高些，B图(J)也是大些，LKIM 放附图吧」）**：
> 主图从 **7 面板** 砍到 **3 面板 = A + J + D**，**1 列 × 3 行全部通栏**；
> **L/K/I/M 移入新建附图 `Fig_GojiVirome_Discovery_Supp`**（见下节）。
>
> | 项 | v19 | v20 |
> |---|---|---|
> | 主图面板 | A/J/D/L/K/I/M | **A / J / D** |
> | 网格 | 3 列 × 5 行 | **1 列 × 3 行（全通栏）** |
> | 行权 | `[1.05,2.10,1.25,0.95,0.02]` | **`[1.45,1.25,1.00]`** |
> | 画布 | (8.55, 12.6) | **(8.55, 11.6)** |
> | hspace / wspace | 0.62 / 0.62 | **0.34 / —** |
> | 实测 A | ~4.7 × 1.6 in（A 跨 2 列） | **7.61 × 3.37 in**（通栏） |
> | 实测 J | ~1.8 × 1.6 in（J 占 1 列） | **7.61 × 2.91 in**（通栏） |
> | 实测 D | 7.61 × 2.95 in | 7.61 × 2.33 in |
>
> ⚠⚠ **v20 第一版踩过的坑（勿重蹈）**：首版把 A/J 并排（A 跨 2 列、J 跨 1 列）+
> 画布压到 10.2 in → 实测 **J 仅 1.80 × 5.26 in**（宽高比 0.34，散点被压成竖线），
> 且 J 的 `get_position().width < 3.2` 触发 `panel_J_scatter` 的**内嵌边缘箱线**分支
> → inset 标题 `Distribution (box = IQR)` 与 J 主标题**叠在一起（实测仅 2 px 间距）**。
> → **正解 = 三面板纵向堆叠、全部通栏**（A/J 既拿满宽又各自拿到足够行高）。
>
> ⚠⚠ **外置边缘箱线时 `ax_top` 绝不能再 `set_title(...)`**：`ax_top` 起于轴顶上方
> 0.02 轴分数（≈0.06 in），其 title 与主轴 `set_title(loc="left")` 竖直只差半行 → 叠字。
> v20 修法 = 外置分支把说明文字改画在**箱线右端外侧**（`ax_top.text(1.01, 0.5, ...)`），
> 内嵌分支才用 `set_title`。
>
> ⚠ 字母偏移：通栏面板（轴占画布 ~91 % 宽）用默认 `dx=-0.075`（≈-0.57 in）仍落在画布内
> （轴左缘在画布 8.5 % ≈ 0.73 in）；**待检项 = 无 `TEXT-OFF-CANVAS`**（体检实测 0 例）。
>
> ⚠ 画布高 11.6 in 是**反解**出来的：row 权 1.45/1.25/1.00 + 间隙 → A≈3.37 / J≈2.91 /
> D≈2.33 in。D 高 2.33 in 下 37 个 L5 节点各得 ~4.5 pt（历史判据 ≥ 4.2 pt，勉强达标）。
> ⚠ 实际产出因 `bbox_inches="tight"` 会略变（实测 10.15 × 11.29 in）。
>
> ⭐ **本主图的三张伴生附图**（由 `main()` 自动调用，失败只 warn 不拖垮主流程）：
> ① `Fig_GojiVirome_NoveltySankey_Supp`（09b 明细，`fig_novelty_sankey_supp.py`）；
> ② **`Fig_GojiVirome_Discovery_Supp`（v20 新增，`fig_discovery_supp.py`：L/K/I/M 四面板）**。
>
> ⭐ **行高按「宽面板源图宽高比」定，目标是面板比 ≈ 源图比**（否则 `imshow(aspect='equal')`
> 留白）。源码实测比：D=2.87 / I=2.85 / A=2.30 / M=2.25（跨 3 列）→ 行 1 取 1.30 使
> **D 面板 = 4.70 × 1.63 in，比 2.89 ≈ 源 2.87**。同时把桑基 plotly `margin` 的
> `r` 从 110 收到 60、渲染尺寸改 **1760×610（≈2.885）**，D 区右侧留白 **20.7 % → 13.7 %**。
>
> ⚠⚠ **两个坑（本轮实测）**：
> ① **`panel_D_sankey` 的标题有两份副本**（直渲分支 + 旁路分支）→ 改标题时只改一处，
>    图上仍是旧文字。**两处必须逐字一致**。
> ② **旁路资产不会随渲染参数自动更新**：改了 `to_image(width=…)` / `margin` 后，
>    `_assets/sankey_flow_plant.png` 仍是旧的 → 必须**先删该资产**再跑，否则静默沿用旧图。

**配套独立图 `Fig_GojiVirome_Phylo`**（对齐论文 Fig3–5；见 `make_goji_virome_phylo_figs.py` 节）。

```bash
export LD_LIBRARY_PATH=/home/zhangwenda/mambaforge/lib:$LD_LIBRARY_PATH
export MMPV_IO_LAYOUT=legacy          # 跑 legacy 结果树必须
python scripts/barbarum/make_discovery_sci_figs.py \
  --root   ~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out \
  --outdir ~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/10_Reports/S10_SCI_Discovery
```

产出（v20）：主图组合图 **A+J+D 三面板**（png 实测 10.15 × 11.29 in，`bbox_inches="tight"` 后）
+ `panels/panel_{A,D,I,J,K,L,M}.{png,pdf}`（**7 个单面板**；I/J/K/L/M 虽不在主图里，
但仍是独立交付件 → 保留、不 prune）+ `_build/fig_manifest.json`
+ `_assets/sankey_flow_plant.png`
+ **三张伴生附图**：`Fig_GojiVirome_NoveltySankey_Supp.{png,pdf}`（09b 明细）
与 `Fig_GojiVirome_Discovery_Supp.{png,pdf}`（v20 新增，L/K/I/M）。

⚠ `main()` 里两张附图都**用 try/except 包住**（失败只 warn）：主图/面板已落盘，
附图失败不应让整条命令报错退出（用户拿不到主图更糟）。

### ⚠ 桑基图与服务器（kaleido 无 Chrome）

matplotlib 无原生桑基 → 用 **plotly + kaleido** 渲染 PNG 再 `ax.imshow` 贴图。
**服务器没有 Chrome**（`Kaleido requires Google Chrome`，且无外网装不上）→
脚本加**旁路**：若能在 `outdir/`、`outdir/_build/`、`outdir/_assets/`、
脚本同级 `_assets/` 里找到 **`sankey_flow_plant.png`**（v14 起为植物子集版；
旧名 `sankey_flow.png` 为全宿主版），就直接贴它（本地渲好上传）。
全都找不到时才退化为 `_panel_D_bars`（Phylum top-12 横条），图面不空。

**桑基数据源**：必须用 `final_integrated_classification.tsv`（逐 contig 全谱系
20,935 行）；**不能**用 `taxonomy_composition.tsv`（只有 `Rank/Taxon/Count`、
各阶元独立汇总、**无 lineage** → 桑基的层级边无从谈起）。
**v14 起面板 D 换植物子集**：以 `ensemble_host_summary.tsv` 的 `Final_Host == "Plant"`
做 contig 过滤（1,063 条，与 `plant_denom()` 同源）→ 渲染尺寸 **`1760 × 610`**
（v17 由 `1780 × 620` 收窄，配 `margin r=60`；改后 D 区右侧留白 20.7 % → **13.7 %**）。
⚠ **过滤键必须用 `ensemble_host_summary.tsv`**，不能用会漂移的 `plant_final_taxonomy.tsv`
（旧快照写 1,092）。

### 实现要点（踩坑记录）

> ⚠ **字母口径提示（v17, 2026-10-13）**：下列踩坑记录按**当时的字母**记述，
> 其中若干面板已下线或改编号（E/F/G/H、B/C 已删；D 自 v14 起 = 植物桑基而非条形）。
> **记录保留原文不改**（是真实历史与判据来源），读时按下表折算：
> 旧 F→E(删) ｜ 旧 G+H→F(删) ｜ 旧 I→G(删) ｜ 旧 J→H(删) ｜ **现行 = A/D/I/J/K/L/M**。

1. **`└ Novel Rank` 是子行**：`stage_summary.tsv` 里 `Known=/NewSp=/NewGe=/NewFa=`
   在**以 `└ ` 开头的子行**上，主行只有 `"20935 条, ★17282 已知, ★★1125 新种"`。
   `parse_stage_summary()` 会跳过 `└ ` 行 → 面板 F 曾整块空图。修法 =
   扫全表（含子行）正则 `Known=(\d+)\s*NewSp=(\d+)\s*NewGe=(\d+)\s*NewFa=(\d+)`。
2. **`All_Candidate` 是并集不是工具**：`ident_summary.tsv` 的 `All_Candidate`
   = 556,773（多工具并集），当工具画会插出一根超高条。修法 = 从工具列表剔除，
   改画成 `axvline` 参考线 + 箭头注记；同时过滤 `n > 0`
   （`virsorter2`/`virhunter`/`viralm` 为 0）。
3. **面板 A 的植物块不是漏斗延续**：`Plant-host (subset) 1,063` 是**植物子集**
   （`06_HostPrediction` 的 `Plant=` 行），与前面的 20,892 全宿主 vOTU **不同轨**。
   直接连成漏斗会显示 20,892 → 1,063 的"假漏斗"（掉 95 %）。修法 = v18 起用
   **右侧 scope 括号**（`all hosts` / `plant subset`，轴分数坐标）+ 括注标清；
   括注数字从 `06_HostPrediction` / `checkv_summary.tsv` / `08_Rescue` /
   `keep_summary.tsv` / `10_Virome_Analysis` **动态取，不硬编码**。
4. **面板 H 不能用 100 % 堆叠**：HQ 只占 200/20,889，堆叠后 HQ 段被压平。
   改为 HQ **绝对数**横向条 + 组内百分比注记。
5. **`Genome_Type` 有 7 行空值** → 图例标签会变成 `"undefined"`。
   收集时先归一 `g = gt if gt else "none"`。
6. 面板 C 标题换行会把 D 面板标题顶穿 → 标题回单行，注记移到 xlabel。
7. 面板 E 的 `Realm`/`Kingdom` 在 45° 旋转下会侵入零轴 → 改 `rotation=90`。
8. **分类面板必须标口径**：面板 **D 是植物宿主子集**（v17 起：1,063 条，
   `final_integrated_classification.tsv` 按 `ensemble_host_summary.tsv` 的
   `Final_Host=="Plant"` 过滤后画桑基），而面板 **I 是植物宿主子集**（1,090 条，
   另一套口径 —— 见 R6）。两者并排不标口径会被读成"枸杞里的病毒科"（用户实际提出过此疑问）。
   现 D 的标题 = `Taxonomy flow · plant-host subset (n = 1,063)`、
   I 的 xlabel = `Contigs · plant-host subset (n = 1,090)`
   （n 由 `all_plant_viruses_genus_summary.tsv` 的 `n_contigs` 动态求和，不硬编码）。
9. **⚠ 面板 D 只有 133.6 pt 宽，而它原有标题就已 162.4 pt**（本就伸进 73.5 pt 列间隙）。
   加长口径**只能挂 xlabel**（可轻微越出轴框且不撑 `bbox_inches="tight"`）；
   写进标题会到 177.1 pt，虽不重叠但更挤。**改标题/xlabel 前先实测宽度**：
   `ax.set_title(...)` 后 `t.get_window_extent(fig.canvas.get_renderer()).width`。
   改完核对 PDF 尺寸未变（`9.26 × 11.46 in`）即证明未撑大画布。
10. **⚠ `stage_summary.tsv` 会被 `report_pipeline.py` 本地化重写**：该脚本内置
    `_en()` 中译英（`report_pipeline.py:139`），重跑后 `20935 条` → `20935 seqs`、
    `09b_分析验证` → `09b_Analysis_Verify`（实测 **31 行差异**）。`parse_stage_summary()`
    的取数正则必须**同时认中英文**（现为 `([\d,]+)\s*(?:条|seqs)`），否则面板 A 会
    **静默**丢 `Taxonomically classified` / `Host-assigned` 两条（`A.n_items` 下掉 2）。
    **每次重跑 `report_pipeline.py` 后都要复核主图**（判据：manifest 的 `A.n_items == 13`）。
11. **面板 A 补全为 13 环节全链（v18，2026-10-13 第 12 轮）**：用户要求
    「组装-病毒-过滤-cobra-merge-cluster-分类-宿主-植物-HQ-rescue-high-HQ」的数目。
    `funnel()` 现为 **13 条**，分**两块口径**（全宿主 8 条 + 植物子集 5 条）：

    | 条目 | 值 | 出处 | 语义色 |
    |---|---|---|---|
    | `Assembled contigs` | 14,387,168 | `01_Assembly` contigs | `assembly=#0072B2` |
    | `Viral sequences (10 tools)` | 556,773 | `02_Identification` viral_seqs | `identification=#56B4E9` |
    | `Filtered candidates` | 94,826 | `filter_per_sample.tsv`（02b_Filter uniprot fasta 计数） | **`filter=#999999`（v18 新键）** |
    | `COBRA queries (extended)` | 35,348 | `cobra_summary.tsv` `Total_Queries` | `cobra=#009E73` |
    | `Merged / re-assembled` | 80,637 | `03b_Merge` combined | `cobra=#009E73` |
    | `vOTU clusters` | 37,803 | `cluster_pipeline_reduction.tsv` | `cluster=#E69F00` |
    | `Taxonomically classified` | 20,935 | `05_Taxonomy` classified | `taxonomy=#CC79A7` |
    | `Host-assigned` | 20,892 | `06_HostPrediction` classified | `host=#D55E00` |
    | `Plant-host (subset)` | **1,063** | 同上 `Plant=`（正则 `Plant=([\d,]+)`） | **`host=#D55E00`（故意与 `Host-assigned` 同色）** |
    | `Plant HQ vOTUs` | **63** | `checkv_summary.tsv` **Plant 行** `HQ`（33+30） | `checkv=#F0E442` |
    | `Rescued (VSI, plant subset)` | 276 | `08_Rescue` rescued | `rescue=#444444` |
    | `CheckV high-confidence` | 45 | `keep_summary.tsv` `checkv_confidence=="high"` | `checkv=#F0E442` |
    | `Final HQ plant viruses` | 291 | `10_Virome_Analysis` `HQ viruses:` | `final=#000000` |

    ⚠ `Plant-host (subset)` 与 `Host-assigned` **同色是刻意的**：让读者一眼看出前者是后者
    的切片，而不是又一条独立口径。若给不同色反而会被读成"另一路统计"。
    ⚠ **`filter` 必须显式进 `SEMANTIC["stage"]`**（取值 `#999999`）：未登记的键会
    **静默**退化为 `#999999`，与"有意用中性灰"无法区分 —— 登记后注释可说明意图。
    ⚠ **`Plant HQ` 口径**：`checkv_summary.tsv` Plant 行分母是 **1,090**（≠ 本图其它处的
    1,063）→ 属已知 **R6** 两套植物口径分歧；用户 2026-10-13 明确选定此口径
    （不取 `07_CheckV` 全宿主 200）。图上按用户要求**逐项标注来源**。
    ⚠ **`Filtered candidates` 不能读 `filter_summary.tsv`**：本项目里它是 **396 行全 0
    的空表**（`ident_stats.py` 按 `uniprot_filter_output*` 目录名找产物，与实际 `uniprot/`
    不匹配 → 静默全 0）。改从各样本 `uniprot/all_candidates_uniprot.fasta` 的 `>` 计数
    汇总，落盘 `10_Reports/filter_per_sample.tsv`（395 行）。
    ⚠ 服务器 `02b_Filter` 下的 fasta 是**断裂软链**（指向已移走的数据盘，
    `grep -c "^>"` 全 0）→ 必须从**子目录实体文件**取数。
12. **⚠ 带箭头的注记别用固定倍数 → 版式体检假阳性；v18 起改「轴分数坐标」根治**：
    （历史）面板 A 的注记箭头尾端取 `vals[-1]*3.0` 时，体检报
    `[OVERLAP 22.0x4.7px] '1,063' <-> 'plant track (subset of Host-assign…'`。
    **实测 bbox 定位**（`ax.annotate(...)` 后 `t.get_window_extent(ren)`）：
    注记 bbox `x[162.0,469.4]`，而**文字实际左缘是 323.2 px**——`x0=162.0` 正是
    **箭头尾端**（`xy`）。即真墨迹与 `1,063` 标签**并不相碰**，纯 bbox 并集造成的假阳性。
    （历史修法）把箭头尾端挪到数值标签右侧 `xy=(ph*4.0, y_anchor)`。
    ⚠ **v18 起不再用箭头**：改画**右侧 scope 括号 + 括注**，且**一律用 `ax.transAxes`
    轴分数坐标**（`x_brace=0.735`、注记 `0.755`）。理由：早前按 `max(vals)×N` 估数据坐标，
    在 log 轴上要**反复反解 xlim**，实测注记被顶出画布 → 体检报 `TEXT-OFF-CANVAS`
    （`bbox=(641,86,743,139)` 而 canvas 只到 700）。轴分数与数据范围无关，且
    **括号线左侧放竖排组名**（`x_brace-0.028`），注记放右侧，两者不抢位（实测 0 对）。
    唯一的数据坐标动作 = 把 `xlim` 右界由 `2.2×` 撑到 `4.5×`，给条形右侧留横向净空。
13. **⚠ `Annotation.get_window_extent()` 会把箭头并进 bbox**（保留教训）：带箭头的注记
    不要用"贴着末条"的固定倍数；若仍用箭头，尾端要显式让开**数值标签**的 x 区间。

### v13 桑基与版式（本轮新增坑）

13. **⚠ plotly `node.color` 不接受 `None`**：给"只有部分节点要灰"写 `[None, None, …]`
    会直接抛 `ValueError: Invalid element(s) received for the 'color' property
    of sankey.node`。修法 = 每个节点都算出一个实色（灰类给 `"#B0B0B0"`）。
14. **⚠ 手写 `node.x / node.y` 会破坏 `arrangement="snap"` 的去重叠**：实测传了显式坐标后
    节点**挤成一列**、标签互相压住、首末列被裁。修法 = **不传 x/y**，让 plotly 自动布局
    （只在数据层按秩给顺序），颜色仍按层保留。
15. **⚠ 桑基 `margin` 太小会裁掉首末列标签**：`margin=dict(l=2,r=2)` 时只看得到色块、
    看不到 `Realm`/`Genus` 名字。修法 = `margin=dict(l=110, r=110, t=10, b=10)`。
16. **⚠ 面板 F 的百分比标签三坑**：① 画在**条形末端**（`q+eps`）→ 小值全挤在起点叠字；
    ② 改成"统一右列"取 `1.02·xmax` → 撞最长条的值标签（`7,595`）；③ 图例放图内右下/中右
    → 碰末排百分比标签。**修法** = `xlim=(0, 1.75·xmax)` + 百分比列固定 `1.46·xmax`
    + 图例 `bbox_to_anchor=(0.0, -0.22)` 放**图外下方**。
17. **⚠ 图例 `lower right` 右缘会伸出轴框、撞相邻面板的 y 刻度**（跨面板 18.4×1.8 px
    `unassigned` ↔ `Fungal virus/contam.`）。修法 = `wspace 0.42→0.62` +
    `bbox_to_anchor=(0.995, 0.02)` 收进轴内。
18. **⚠ 面板 I 用 `axis("off")` + `inset_axes` → `ax.has_data()` 恒为 `False`**：
    (a) 字母标号不会被画（`if ax.has_data() or ax.texts` 判据失效）→ 须显式加
    `or letter == "I"`；(b) 且 I 的轴几乎占满画布（`[0.06,0.16,0.90,0.72]`）→ 用默认
    `dx=-0.20` 会把字母推到**画布外**（实测 bbox `x0=-187`）→ I 用专用 `dx=-0.055`。
19. **标题里的硬编码字母前缀会与 `panel_letter` 重复**：`_hbar_labels(title="A  Discovery…")`
    + `panel_letter(ax, "A")` → 图上出现 `A   A Discovery…`。修法 = 标题只写内容
    （`"Discovery pipeline reduction"`），字母交给 `panel_letter` 统一画。
20. **⚠ 服务器 kaleido 无 Chrome**：`fig.to_image()` 抛 `RuntimeError: Kaleido requires
    Google Chrome`，`plotly.io._kaleido get_chrome` 也装不上（无外网）。修法 = 加
    `_find_asset("sankey_flow.png", D)` 旁路（查 `outdir/` → `outdir/_build/` →
    `outdir/_assets/` → 脚本同级 `_assets/`），本地渲好上传；全找不到才退化横条。
21. **⚠ `vert=` 在 matplotlib 3.11 弃用、3.13 移除，但旧版不认 `orientation=`**：
    本地 mpl 3.11（告警）与服务器旧版（不认新参数）**互斥** → 写死任一个都会在另一端
    报错/告警。修法 = `_box()` 版本探测后择一传参。
22. **⚠ 散点图的边缘箱线位置必须按**轴宽**探测**：组合图里 J 只占 1 列（≈2.5 in）时，
    外侧 `inset_axes([1.02, 0, 0.16, 1])` 会伸进**相邻面板 I 的 y 轴标签**。修法 =
    `ax.get_position().width × figsize` < 3.2 in 时把边缘箱线改放**内侧**。
23. **⚠ 箱线图 x 轴不能写死 `0–102`**：本数据所有科的中位数都在 78–100 %，固定
    0–102 会让左 2/3 全空白、箱线被压成细条（实测）。修法 = 按实际范围留边（但保留
    80 % 参考线的可见性：`x0 < 80 < x1` 才画）。
24. **⚠ 参考线/竖排注记必须在 `set_xlim`/`set_ylim` **之后**画**：先画时 `get_ylim()`
    还是默认 `(0,1)` → 注记落到画布外被裁（实测）。且应取 `lo_y`（数据底）而非 0。
25. **⚠ 跨平台字形高度差会造出**服务器独有**的重叠**：G 面板 15 行、行距 ~7 px，
    `val_fs=5.6` 在本机不碰，服务器上数值标签竖向叠 **0.8 px**（14 对 OVERLAP）。
    **教训**：版式体检必须**双端各跑一遍**，本机 0 对 ≠ 服务器 0 对。修法 = G 的
    `val_fs` 单独降到 4.8。
26. **identity 散点/箱线的数据源只能是 `rescue_evidence_scored.tsv`**：它的
    `nt_pident`/`nt_qcovs` 才是 % identity/coverage；`All_plant.viruses_info.tsv` 的
    `confidence`（1/0.88/0.75…）是**证据分**，拿来画 identity 是张冠李戴。
27. **⚠ 结构域完整性 ≠ 基因组完整性**（v15 面板 L）：`checkv_summary.tsv` 的
    Complete/Fragment 是**基因组**完整度；论文 Fig2-E 要的是**保守域**完整度 →
    唯一正确来源是 `cdd_evidence_report.tsv` 的 `evidence` 列
    （`PASS_CORE`/`PASS_CORE_FAM`→Complete、`PASS_VIRAL`→Fragment、其余→No conserved domain）。
    实测 n=291（Complete 112 / Fragment 34 / No conserved domain 145），
    **与 `checkv_summary` 的 92/108/… 完全不同**——两者是不同的问题。
28. **⚠ 无现成"功能注释"列时按关键词派生**（v15 面板 M）：`hmm_ct3_hits.tsv` 只有
    `ct3_best_viral_desc`（自由文本描述）→ 用 `_TARGET_RULES` 有序规则表归类
    （RdRp / RT·RNase H / Capsid / Nucleoprotein / Helicase / Protease / Movement /
    Replicase / Integrase），**并把未命中数写在 xlabel 上**（`n=167, 5 unmapped`），
    而不是静默丢弃——否则读者会以为总数就是 167。
29. **⚠ `ax.set_axis_off()` 不隐藏 axis label**（phylo 脚本）：`set_axis_off()` 只去轴线/
    刻度/spine，**不动 xlabel/ylabel** → Bio.Phylo 自加的 `branch length`(xlabel)/
    `taxa`(ylabel) 仍会与相邻子图判为重叠（实测 **226 对**）。修法 = `ax.set_xlabel("")`
    + `ax.set_ylabel("")` 显式清空 → 226 → 6 对。
30. **⚠ Bio.Phylo 会给标签加前导空格**（phylo 脚本）：`Phylo.draw` 把标签写成
    `" " + label` → `t.get_text()` 带空格，与 `_short()` 的无空格值做 `in` 判定**恒 False**
    （表现 = 错行代码"看起来没生效"、bbox 纹丝不动）。修法 = 三处判定
    （`label_colors` lambda / 字号设置 / 错行）全部 `.strip()`。
31. **⚠ 自有 contig 缩写必须带样本前缀**（phylo 脚本）：不同样本各有 `NODE_1`，
    只取 `NODE_{n}` 会撞名（判为重叠）。修法 = `_short()` 取样本后 6 位 + `:NODE_{n}`
    （如 `191322:NODE_1`）。另：**Biopython 1.87 起 `label_colors` 回调收到的是标签
    字符串**（旧版是 clade 对象）→ 不能再调 `c.is_terminal()`（`AttributeError`）。

32. **⚠⚠ 自有 contig 判定漏认 `contig_*` 命名（2026-10-12，用户截图暴露）**：
    `_OWN_RE` 原先只写 `r"^(?:CRR|SRR)\d+_clean_NODE_\d+"` → **只认一种命名**，
    漏掉 `contig_*`。实测后果：`Benyviridae` 树里明明有 `contig_4185`，却判为
    `own = 0`（管线 HTML 报告里它是**有高亮的**）→ 该面板被当成"无自有序列"；
    `Closteroviridae` 同样少算 1 个。**修正**：正则改
    `r"^(?:(?:CRR|SRR)\d+_clean_NODE_\d+|contig_\d+)"`，`_short()` 对 `contig_N`
    **原样返回**（`contig_N` 已全局唯一，**不要**硬套 `NODE_` 前缀）。
    复核：`Benyviridae` 0 → **1**、`Closteroviridae` 2 → **3**，其余科不变
    （`Caulimoviridae` / `Virgaviridae` 仍为 0，是**真**数据事实）。
    **教训**：判"树里有没有自有序列"**不能只按一种命名 pattern 匹配** ——
    先 `grep -o` 一遍 treefile 把实际出现的命名样式列全，再定正则；
    且要把管线报告已渲染的图当**对照物**（不一致时先怀疑自己的判据）。

**验收（双端一致，各跑一遍）**：`_check_discovery_overlap.py` → 主图（**v18，7 面板**）
**16 figures / 0 overlapping pairs**、系统发育图 **44 figures / 0 overlapping pairs**
（本地与服务器**各 72 figures / 0 对**）；组合图文本层可检索全部 7 个面板标题、CJK = 0。
脚本 md5：`make_discovery_sci_figs.py` 本地 **`4762851bc8656468a838bf7bf588fd40`**
（**v21.3：行权 `[0.6,1.0]` + 画布 `(8.55,14.0)` + `hspace 0.28` + 单面板 D `7.61×6.87`**）
/ 服务器 **`15da7e3d7547021032cbcd4e2db3b4e1`**（v17，**尚未同步**）；
`fig_novelty_sankey.py` 本地 **`12b18e3f8463ee1eeccc0db240b25f8b`**
（**逐层 `level_fs`/`node_fs_scale` + `_nw()` + `--all-fs`（含独立图）**；**未上服务器**）；
`fig_discovery_supp.py` 本地 **`37b1dc4a6d5772ca1883a61abbac7767`**（v20，**全新、未上服务器**）；
`fig_novelty_sankey_supp.py` 本地 **`636afee9bc159f68c3189bb992095819`**（**未上服务器**）；
`sci_palette.py` 本地 **`f8fbe4dad8c00ef75175bb8acb5159bc`**（v18，新增 `stage["filter"]`）
/ 服务器 **`1aa6769178f001bd99e2add79b238d2e`**（v17，**尚未同步**）。
⚠ **本轮（v21 + L1-L4）改动只在本机验证**；服务器需同步 `make_discovery_sci_figs.py` +
`fig_novelty_sankey.py` + `fig_discovery_supp.py` + `sci_palette.py`
+ 新数据源 `10_Reports/filter_per_sample.tsv` 后**双端各跑一遍**才算闭环
（⚠ 桑基三脚本禁整篇覆盖，见 R8+R13）。
⚠ **同步 `fig_novelty_sankey.py` 时注意**：仓库里存在**并发会话**改过该文件的记录
（见 §「并发会话」）→ **先 `diff` 再 `scp`**，禁止整篇覆盖。
⚠ `Fig_GojiVirome_PhyloSynteny` 是**图片嵌入型**图（`AxesImage`，Text 极少）→
`_check_discovery_overlap.py` 对它**不适用**（它只测 Text bbox）；该图靠
`phylosynteny_manifest.json` 的 `height_capped` / `print_scale` 复核版式。

## make_goji_virome_phylo_figs.py

**枸杞病毒组系统发育独立图** `Fig_GojiVirome_Phylo`（2026-10-12 新建，对齐参考论文
`Viruses 2025, 17, 624` 的 **Fig3–5**：Kitrinoviricota RdRp / Pisuviricota RdRp /
Caulimoviridae RT 的 ML 树）。

**产出**：组合图 `Fig_GojiVirome_Phylo.{png,pdf}`（**8.886 × 26.467 in**、PDF 文本层
4,838 chars、**CJK = 0**）+ `panels/phylo_{Family}.{png,pdf}`（**21 个科**各一张）
+ `_build/phylo_manifest.json`。

**输入**：服务器 `09b_Analysis_Verify/acvirus_trees/{Family}/phylogeny.treefile`
（**21 个科级 ML 树**）。树叶两类：`NCBI accession` = 参考序列；
`CRR*/SRR*_clean_NODE_*` = **我们自己的 contig**（判"新病毒"的依据）。
⚠ **自有 contig 有两种命名**（实测 2026-10-12）：`CRR*/SRR*_clean_NODE_*` **与**
`contig_*`（Benyviridae / Closteroviridae 用后者）。正则须同时认两种 —— 见"坑 32"。

```bash
export LD_LIBRARY_PATH=/home/zhangwenda/mambaforge/lib:$LD_LIBRARY_PATH
python scripts/barbarum/make_goji_virome_phylo_figs.py \
  --trees  <.../09b_Analysis_Verify/acvirus_trees> \
  --outdir <.../S10_SCI_Discovery> \
  --virus-name "Lycium barbarum"
```

**结构要点**：
- `PHYLUM_FAMILIES` 把 21 个科按门分组（Kitrinoviricota ×11 / Pisuviricota ×4 /
  Artverviricota(RT) ×3 / Cressdnaviricota ×1 / Negarnaviricota ×1 / Unclassified ×1），
  `PANEL_META` 给每组一个标题后缀（如 `Kitrinoviricota · RdRp (ML)`，对齐论文 Fig3/4/5）。
- **组合图 = 2 列网格**，叶数 **≥ `WIDE_AT=60`** 的科（Rhabdoviridae 173、Caulimoviridae 94、
  Draupnirviridae 65）**跨 2 列独占一行**；画布高 `1.0 + 1.30 × 行权`。
- **自有 contig 橙色加粗高亮 + 字号加大**（own 5.0 pt bold / ref 3.4 pt），
  对齐论文"新种红标"的语义。
- `_declash_own_labels()` 只对**自有标签**按 y0 升序做**链式**下移（`min_gap_px=6.0`），
  参考 accession 不动。

**⚠ 数据事实（非画法问题）**：`Caulimoviridae`（94 叶）、`Virgaviridae`（33 叶）
树内**自有 contig = 0** → 对应论文 Fig5 的 Caulimoviridae 面板**无自有序列可标橙**；
`Caulimoviridae` 的 `analysis_sequences.fasta` 只有 94 条**参考**基因组。这是管线数据事实。
（⚠ `Benyviridae` **不**属于这一类 —— 它有 `contig_4185`，见坑 32。）

## make_goji_virome_phylosynteny_fig.py — 「ML 树 + Synteny 两联」独立图

**背景 / 目的**：用户指令「给系统发育产物**补 Synteny 右联**」。

**⭐ 关键决策（用户已定，2026-10-12）**：
- **嵌入管线原图，不重绘** —— 管线 `virome_discovery_pipeline/utils/acvirus_tree_pro.py`
  的 `draw_tree_synteny_figure()` 已为**每个科**产出
  `{Family}/Tree_Synteny_Composite.{png,pdf}`（版式 = 左 `Maximum-Likelihood Phylogeny`
  + 右 `Comparative Synteny & Homology`：含 `Genome Coordinates (bp)` 轴、
  `Protein Identity (%)` 30–100 色标、`Sample_Contig` 蓝底高亮、属级背景色块、支持值）。
  这正是参考论文 / 管线报告 HTML 的版式。
  项目硬规则「整页独立图只能自绘」**已被否决** → **原图直接嵌入，不许重绘/降采样**。
- **保留纯树版**：`Fig_GojiVirome_Phylo`（上面那个脚本）**不动**；本脚本**另出**
  `Fig_GojiVirome_PhyloSynteny`。

**产出**：`Fig_GojiVirome_PhyloSynteny.{png,pdf}`（21 行全宽，一科一行）
+ `panels_phylosynteny/phylosynteny_{Family}.{png,pdf}`（21 个科各一张，原图整幅嵌入）
+ `_build/phylosynteny_manifest.json`（每科像素尺寸 / 行高 / 缩放 / 是否被限高）。

```bash
python scripts/barbarum/make_goji_virome_phylosynteny_fig.py \
  --trees  <.../09b_Analysis_Verify/acvirus_trees> \
  --outdir <.../S10_SCI_Discovery> \
  --virus-name "Lycium barbarum"
```

**版式要点**：
- 源图一律 **4477 px 宽**（Potyviridae 4513）→ **通栏**（`COL_W_1 = 9.32 in`）→
  面板净宽 `9.32 - 0.2 = 9.12 in`，`imshow(aspect="auto")` 填满按原比例算好的框。
- **行高 = 原比例高度，∈ [叶数 × 0.145 in, 9.0 in]**：下限保证多叶科不被压扁，
  上限防 `Rhabdoviridae`（原比例 39 in）一科独占全图。超限则**整幅等比缩小**
  （**不裁切、不重排**），manifest 标 `height_capped`。
- 每行上方一条 `Family (n=…, own=…) · Phylum` 题注带（`CAP_IN = 0.30 in`）。
- ⚠ **suptitle 的 `y` 是 figure 分数、留白是英寸** → 必须显式换算
  （`y_title = 1 - TOP_IN*0.22/H_FIG`），否则图变高后两行标题会压住首行图。

## fig_discovery_supp.py — Discovery 补充附图（主图 v20 拆出的 L/K/I/M 四面板）

**背景 / 目的**：用户指令（2026-10-13 晚）「A图展开高些，B图(J)也是大些，**LKIM 放附图吧**」。
主图 v19 的 7 面板把 A（漏斗）与 J（散点）挤到极小（A 高 ~1.6 in、J 宽 1.8 in）→
把 **L/K/I/M 四面板移入独立附图**，主图只留 A/J/D（见 §make_discovery_sci_figs.py v20）。

⭐ **口径零漂移的关键设计 = 复用主图模块的面板函数**：
```python
m = _load_main()                       # 按显式路径 exec make_discovery_sci_figs.py
m.apply_sci_style(base_size=7.5, grid=False)   # ⚠ 必须：否则吃进程级 rcParams
D = m.Data(root, virus_name); D._outdir = outdir
CELLS = [("A",0,0,1,2, m.panel_M_targets), ("B",1,0,1,2, m.panel_K_family),
         ("C",2,0,1,1, m.panel_L_domains), ("D",2,1,1,1, m.panel_I_pair)]
```
→ 同一指标**只写一份实现**（历史教训：两处各写一遍迟早口径分叉）。

**附图版式**（编号 A–D **独立重编**，与主图字母无对应关系）：

| 附图字母 | 内容 | 复用主图 | 布局 |
|---|---|---|---|
| **A** | Alignment targets of viral contigs | `panel_M_targets` | 通栏（跨 2 列） |
| **B** | Identity by virus family（n≥3） | `panel_K_family` | 通栏（跨 2 列） |
| **C** | Conserved-domain completeness（n=291） | `panel_L_domains` | 1 列（饼图） |
| **D** | Genome type + Viroid species | `panel_I_pair` | 1 列（双轴） |

`GridSpec(3, 2, hspace=0.44, wspace=0.30, height_ratios=[1.10,1.30,0.95],
left=0.135, right=0.985, top=0.905, bottom=0.045)`，`figsize=(7.6, 7.5)`。

⚠ **字母两套体系不得互引**：附图内 A/B/C/D 是**排版位置标号**；主图字母
（A/D/I/J/K/L/M）是**稳定资产 ID**（manifest key / `WIDE` / `panel_letter` 锚定）。
同一张图里 "A" 在两处指不同面板 —— 文档与题注**必须写清是主图还是附图**。

⚠ **`panel_I_pair` 用 `ax.get_position()` 算 inset**：在附图（半栏 ~3.5 in）与主图
（v19 的 1.79 in）下面板宽不同 → inset 分数不变但英寸宽不同，块间净空随之变。
附图里类病毒种名仍截 14 字符（`s[:13]+"…"`）→ ✅ 实测无重叠（体检 0 对）。

**产出**：`Fig_GojiVirome_Discovery_Supp.{png,pdf}`（png 2709×2212 = 9.03 × 7.37 in）
+ `Fig_GojiVirome_Discovery_Supp.manifest.json`。

**调用**：由 `make_discovery_sci_figs.py` 的 `main()` **自动调用**（在 09b 附图之后），
失败只 `warn` 不拖垮主流程；也可手工：
```bash
export MMPV_IO_LAYOUT=legacy
python scripts/barbarum/fig_discovery_supp.py --root <RNA-Lycium_barbarum_out>
```

**体检**：`_check_discovery_overlap.py --module fig_discovery_supp -- --root <root>`
→ 判据 **0 对重叠** + 无 `TEXT-OFF-CANVAS`（实测 0 / 0，canvas 760×750 px）。

## make_discovery_novelty_figs.py — 「新/已知 · 科属种 · 新种 · 命名」独立图

**背景 / 目的**：用户指令（2026-10-13）「根据 09b 的表格补绘类似的图，
覆盖 新病毒 / 已知病毒 · DNA / RNA / 类病毒 · 科 / 属 / 种 · 对属的新种 · 属→goji 的命名」，
参考图为 ① 新/已知物种的**科→种 bipartite 通路图**、② **identity×coverage 散点 + 边缘箱线**。

**产出**：`Fig_GojiVirome_Novelty.{png,pdf}`（三联 a/b/c，`8.6 × 12.0 in`）
+ `panels_novelty/panel_{a,b,c}.{png,pdf}`
+ `_build/novelty_manifest.json`（每个面板的族/种计数、**截断比**、命名清单）。

```bash
python scripts/barbarum/make_discovery_novelty_figs.py \
  --root  <RNA-Lycium_barbarum_out> \
  --outdir <...>/S10_SCI_Discovery \
  --virus-name "Lycium barbarum"
```

**面板口径（⭐ 两张表口径不同，图上必须分别标清）**：

| 面板 | 数据源 | 分档列 | 科/属/种列 |
|---|---|---|---|
| **a** | `09b/keep_summary.tsv`（101 条 keep 集） | `novelty_9a` | `Family`/`Genus`/`Species` |
| **b** | `09b/virus_validation/rescue_evidence_scored.tsv`（291 条） | `category` | `tax_family`/`tax_genus`/`tax_species` |
| **c** | 同 b | — | `nt_pident` × `nt_qcovs` |

**三档语义色**（本图专用，已写进 manifest；全取自 OKABE_ITO）：
`Novel` = `#D55E00`（vermillion）/ `Known` = `#0072B2`（blue）/ `Contaminant` = `#999999`（grey）。

**新种判据**（用户 2026-10-13 选定）：
`len_ratio > 1`（长于属平均）**或** `nt identity ∈ [75, 95] %` 且 `qcov ≥ 70 %`。
panel c 用**空心环**标出（82 个点中 24 个命中）。

**「属→goji 的命名」**：`genus_goji()` 按 `lycii`/`goji`、`-like virus`、`associated`、
`sp.` 关键词从 `Species` 列派生，逐条写进 manifest（panel a 命中 2 条、b 命中 11 条；
其中 `Betacytorhabdovirus lycii` 即枸杞学名直命名的实例）。⚠ 这是**派生量**，
不是策展注释，图面/ manifest 均如实标注规则。

### 本轮新增的硬规则（踩坑记录）

1. ⚠⚠ **bipartite 面板的坐标系必须用「英寸」**：本面板要**按文字物理宽度**决定
   色块宽/截断，而 `_text_w_in()` 返回英寸。若坐标系是 0–1 的轴分数，二者
   **不同量纲** → `_fit_label` 会把「Betaflexiviridae」截成「Be…」（实测首版）。
   用英寸坐标系后 `W = 标签宽 + 两侧内边距` 即为真值。
2. ⚠⚠ **谱系图的标题/注记必须走 `transform=ax.transAxes`**：数据坐标是英寸（0..H），
   漏写 transform 的 `ax.text(0, 1.008, …)` 会落到图**正中**（实测把通路图劈两半）。
3. ⚠ **族数 × 每族物种行数是版式硬约束**：盒子矮了必然把末尾几族挤出画布
   （实测 2.4 in 高时 Kitaviridae 掉出）。解法 = `top_families` 上限 +
   **底部强制对齐**的压缩循环，且**截断比如实写进图面**
   （`⚠ Showing top 10 of 17 families`）。
4. ⚠ **档块竖排名 / 轴下图例都会撞**：竖排档名伸长后与底部图例相撞
   （实测 "Contaminant" × 图例 "Novel" 2.8×1.0 px）。→ 档名与图例**一律放进
   独立的说明轴**（`Reading notes` 卡片），与通路图元不同轴 → 结构上不可能重叠。
5. ⚠ 面板字母与标题**分两个 text 横向错开**（字母 x=0、标题 x=0.030），
   注记也从 x=0.045 起 —— 否则字母 bbox 与注记 bbox 竖向擦边（实测 9×3 px）。

## fig_novelty_sankey.py — 09b 五层桑基「已知/新 → DNA/RNA/Viroid → 科 → 属 → 种」

2026-10-13 新建。用户指令：「`2_Fig_GojiVirome_Discovery.png` 补充 09b 的结果，桑基图，
已知病毒 / 新病毒 · DNA / RNA / viroid · 科属种 · 新种（属-goji 命名，这个属多个的话则 -1 -2）」。

**数据层** = `make_novelty_sankey_data.py`（口径常量集中在此，图脚本不再另写一套）→
`novelty_sankey_data/{contigs,nodes,links,new_species}.tsv` + `manifest.json`
+ ⭐ `novel_dedup_full.tsv`（**附图数据**：去冗余前的全部 Novel 明细）。

**产出**：`Fig_GojiVirome_NoveltySankey.{png,pdf}`（**单张**五层桑基，`10.65 × 3.10 in` 资产尺寸，CJK = 0）
+ `Fig_GojiVirome_NoveltySankey_Supp.{png,pdf}`（**附图**，`7.6 × 9.2 in`，4 面板明细）。

```bash
"C:/Python312/python.exe" scripts/barbarum/make_novelty_sankey_data.py
"C:/Python312/python.exe" scripts/barbarum/fig_novelty_sankey.py --outdir <...>/S10_SCI_Discovery
"C:/Python312/python.exe" scripts/barbarum/fig_novelty_sankey_supp.py --outdir <...>/S10_SCI_Discovery  # 附图
"C:/Python312/python.exe" scripts/barbarum/_check_discovery_overlap.py --root <...> --outdir <...>      # 版式体检
```

**口径沿革**：895（原始）→ 296 → 277 → 190（去 Other + 去 Unassigned）→ 121（+ 只物种层去长尾）
→ **61**（+ **同属 Novel 去冗余 cap=1/属**，2026-10-13 最终定案）。

| 项 | 取值 |
|---|---|
| 主体 | **61** = 40 Known（35 HQ 植物病毒 + 5 类病毒种）+ **21 Novel** |
| L1 tier | `已知 40` / `新 21`（其他 87 条已删） |
| L2 | **RNA 51 / Viroid 5 / DNA 5** |
| L3/L4/L5 | family / genus / species（**17 / 25 / 37** 个节点，无聚合块；边 84） |
| 新种 | `category == "新种候选"` 的 **11** 条（cap=1 后由 21 → 11） |

### ⭐⭐ 2026-10-13 晚间两条口径变更（本轮关键）

**① 同属 Novel 去冗余 cap=1/属**（用户：「同属的新病毒只展示一个，这样避免新病毒过多
引起桑基图不好看，**详细的结果放附图**」）

- 两步：① 同属内**同一物种名**（`<属>-goji-N`）只留 1 条（取 `prevalence_pct` 最高，
  并列取 `contig_id` 字典序最小）；② 若去重后仍 > 1 条，按发生率砍到 **1** 条。
- **Known 完全不动**；**命名不变**（只下采样，不做重新编号）。
- 实测：Novel 81 → **21**（实删 64 条）；`Potyvirus` 24 → 1、`Betacytorhabdovirus` 14 → 1…
- ⚠⚠ **去冗余后不再做二次去长尾**：去冗余会把很多属压到 1 个物种，若再按「单位数 < 2 剔除」，
  会把这 **26 个真实新病毒属整个删掉**（实测 61 → 58），与「展示较完整」直接矛盾。
  → **去冗余本身就是「减少视觉冗余」的手段**，长尾剔除不再叠加。
- 关闭方式：`--novel-dedup-per-genus 0`。

**② 去长尾范围收窄：三层 → 只 species 层**（用户：「只对物种层去长尾」）

- 原「family/genus/species 三层逐层迭代剔除 < 2 单位」→ 改为**只在 species 层剔除**。
- 科/属层保留单例，保证每个病毒都能追到完整的 科–属–种 链。

**③ ⭐ 附图机制**（信息零丢失）

- `novel_dedup_full.tsv`（15 列，含 `shown_in_sankey` True/False）记录**去冗余前的全部 Novel
  明细**（85 条 = 展示 21 + 折叠 64），**不参与桑基几何**，纯作「详细结果」附图数据源。
- 附图 `Fig_GojiVirome_NoveltySankey_Supp` 4 面板：A 汇总（shown/folded 按 category）、
  B 属级（每属 shown+folded 堆叠）、C 折叠明细散点（发生率 vs 长度，空心=折叠）、
  D 折叠 contig 清单（按属分组，前 24/64 条）。
- ⚠ 本附图**图中文本一律 ASCII**（rcParams 无 CJK 字体；category 用英文映射
  `CAT_EN`：新种候选→candidate new species 等）。

**源表与关系**（本地镜像 `reports/_discovery_sci_20261012/data/`）：

| 表 | 行 | 关键列 |
|---|---|---|
| `10_Reports/HQ_plant_viruses_info.tsv` | 291 | `Accession`（== rescue 的 `contig_id`，**291/291** 一致）/ `Family` / `Genus` / `Species` / `Genome_Type` |
| `09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv` | 291 | `contig_id` / `category`（8 类中文） |
| `10_Reports/Viroid.all_info.tsv` | 604 | `contig_id`（**第 2 列**）/ `family` / `genus` / `species` |
| `09b_Analysis_Verify/keep_summary.tsv` | 101 | `contig_id`（**第 2 列**）；`keep ⊂ rescue ⊂ HQ`，仅做自检 |

**自检**（不过就 `sys.exit`）：① `set(HQ) != set(rescue)` → 退出；
② `viroid ∩ HQ ≠ ∅` → 退出；③ 开去冗余时**跳过** EXPECTED_LT 精确比对（口径已变），
   关去冗余时 `len(rows) != 121` → 退出；④ 去冗余后 Novel 条数**必须严格减少** + 保留集
   必须是原集子集 + 无重复 contig_id（三重断言）。

### ⭐ 为什么现在是**单张**（口径收紧后的设计约束变化）

旧版拆双面板是因为类病毒 **604/895 = 67.5%** 把科/属/种三层挤成各 1 个节点，
291 条病毒要在 32.5% 高度里摊给 152 个种 → `n=1` 节点只有 ≈2 px，标签物理不可读。
经「去 Unassigned + 去 Other + 只物种层去长尾 + 同属 Novel 去冗余 cap=1」后，
总量降到 **61**（k_max = 37 节点），单节点像素充足 → **单张五层即可承载**。

- L1 列节点数 = **2**（Known / Novel）
- L2 = **3**（DNA / RNA / Viroid）
- L3 = **17** family / L4 = **25** genus / L5 = **37** species（**无任何 `Other ...` 聚合块**）
- 每列 Σn 都 = **61**；边数 **84**
- 11 条「新种候选」物种标签**加粗**；L5 标签按 tier 着色（绿 = Known / 红 = Novel）

图面注记显式声明「Node height and ribbon width are strictly proportional to unit count
(every column sums to 61)」与「No long-tail aggregation (all three levels are complete)」。

### 本轮新增的硬规则（踩坑记录）

1. ⚠⚠ **`_check_discovery_overlap.py` 对本图会「假通过」**：该脚本剔除 `ax.axison == False`
   的轴，而桑基图**全部轴都是 `axis("off")`** → 会报「0 对重叠」但实际一个 text 都没检查。
   → 另写 `_check_novelty_sankey_overlap.py`，遍历**所有** Text 图元（含 `fig.texts` / `fig.legends`）。
2. ⚠ **列内排序必须按上一层的分组序**，否则彩带大面积交叉：family 按「其成员的主 genome_type」
   排、genus 按「其成员的主 family」排、species 按「其成员的主 genus」排。
3. ⚠ **y 轴从 0 向上 → 首节点落在底部**。要「大者在上」必须做 `Y - y` 翻转（首版没翻，图是倒的）。
4. ⚠⚠ **小节点标签必须去碰撞**：`n ≥ 5` 的节点在属/种层仍可能只有 ~10 px 高，
   而 4.05 pt 标签需 ~22 px。→ 同层标签按 y 降序做**贪心下推**（最小行距
   `fontsize × 1.30 / 72 / 面板高_in × Y`），越界整组上移；位移 > 0.45×行距的**补灰色引导线**。
5. ⚠ **面板标题 / 注记 / 图例三者会互压**：`fig.text` 大标题在 y=0.975 与 `ax.transAxes` 1.075
   的面板标题只差 0.03 图高 → 首版直接叠字。→ 大标题上移 + `top` 降到 0.902。
6. ⚠ 图例把 `Tier`（方块）与 `Genome type`（圆点）**拆成两个 legend 并加粗标题** ——
   合成一条时 `Known` 与 `DNA` 同为 `#0072B2`，读者会误以为是同一维度。
7. ⚠⚠ **主图面板 D 的标题数字曾硬编码「121」** → 口径改 61 后**面板静默说谎**
   （图看着正常，数字是旧的）。→ 改为从 `novelty_sankey_data/manifest.json` 的 `total` **读**，
   兜底值也是 61。⚠ 同一类坑：`fig_novelty_sankey.py` 里的「新种候选数断言」也曾硬编码 21。
8. ⚠⚠ **去长尾范围与去冗余必须成对设计**：只改去冗余（cap=1）而不收窄去长尾
   → 26 个真实新病毒属被整个删掉（61 → 58），且**零报错**。
9. ⚠ 图脚本自检在**开了去冗余**时必须**跳过**「物种层 1 单位」检查 ——
   cap=1 后每属只剩 1 个物种是**预期行为**（代表），不是漏剔。

### ⭐⭐ 2026-10-13 深夜「图D的前面 L1-L4 变大些」（本轮关键，用户指令）

**用户原文**：`图D的前面L1-L4 变大些`；经澄清确认口径 = **标签字号放大 + 节点块加高（两者都要）**。

**背景实测（先量化，再改）**：主图 D 面板贴的是 `sankey_novelty_embed.png`（`--embed-asset-plain`）。
旧资产 `3130×927`（比 3.376）下，L1 的两个块实测只有：

| 层 | 节点数 | 旧块高（px） | 旧占比 |
|---|---|---|---|
| L1 Known（40 单位） | 2 | **83** | 9.0 % 画布高 |
| L1 Novel（21 单位） | 2 | **43** | 4.6 % |
| L2 DNA/Viroid（各 5 单位） | 3 | ≈11 | 1.2 % |

**根因（三条，必须一起看）**：

1. **Y 被「节点最多的列」绑架**：`solve_geometry` 的 `Y = tot/(1 − c·(k_max−1))`，
   而 `k_max = 37`（species 列）→ Y 被撑到 **406.7**，L1 的 61 个单位只占其中 15 %。
2. **gap 已在地板值，压不动**：gap 是**物理量** = `fontsize × 1.30/72` 英寸，
   4.05 pt @300 dpi → **21.94 px**；旧资产实测行距 22.97 px 只余 5 % 余量。
   → 任何「压小 gap 换大块」的尝试都会让 L5 的 37 个标签互相压住。
3. ⚠⚠ **试过且失败的两种「局部放大」**（都会静默破坏版式，已回退）：
   - **给 L1–L4 的 gap 单独乘 boost**：L3(17 节点)/L4(25 节点) 的列高暴涨到
     **157 % / 228 % of Y** → 直接冲出画布（`列高 > Y`，零报错）；
   - **压小全局 gap**：行距掉到 9–19 px < 21.94 px 地板 → L5 标签互压。

**✅ 唯一正解 = 加高嵌入画布**：`1 单位像素 = 画布高 / Y`，而 **Y 只由 gap 地板决定、与画布高无关**
→ 画布越高，单位像素**线性增长**：

| 画布高 | Y | 1 单位 | Known 块 | Novel 块 | 行距 |
|---|---|---|---|---|---|
| 3.10 in（最初） | 406.7 | 2.29 px | 91 px | 48 px | 22.97 px |
| 4.89 in（第 1 轮） | 163.3 | 8.08 px | 323 px | 170 px | 22.97 px（**恒定**） |
| **6.48 in（现行）** | 122.3 | 13.49 px | **540 px** | **283 px** | 22.97 px（**恒定**） |

⚠ **行距恒 22.97 px** 是关键 —— 加高**不会造成标签重叠**（物理量不随画布缩放）。

**实现（`fig_novelty_sankey.py`）**：

```bash
# 主图 D 的贴图资产：加高画布 + 前四层字号/块宽放大
"...fig_novelty_sankey.py" --outdir <out> \
    --embed-asset-plain <out>/sankey_novelty_embed.png \
    --embed-size 10.65x6.48 \      # ⚠ 从 10.65x3.10 两次加高（唯一有效手段）
    --front-fs 5.4 \               # L1–L4 标签字号（FS_NODE=4.05 → 5.4）
    --front-node-scale 1.5         # L1–L4 节点块宽 ×1.5
```

- `solve_geometry` **保持纯几何**（不再引入 boost），加高由 `--embed-size` 的 H 控制；
  注释里写死了上述实测结论，防止后人再走"局部 boost"的弯路。
- `draw()` 新增两个**逐层**参数：`level_fs={L:fs}`（逐层字号，且**逐层重算行距地板**）、
  `node_fs_scale={L:k}`（逐层块宽，经 `_nw(L, scale)`）。
- ⚠⚠ **`node_fs_scale` 只改「块宽 `NW`」，绝不动「块高」** —— 图注写的
  「Node height and ribbon width are strictly proportional to unit count」**必须继续成立**
  （块高恒 = `n` 个单位 × 统一比例）。
- ⚠ 前四层的 label 偏移随块宽走（`off = _w + 0.05`），否则放大块宽后标签会压在块上。

### ⭐ 2026-10-13 21:16「再宽高些」→ 澄清「D 面板 L1–L2 块再高些，不用宽了」

**⚠⚠ 块高只能靠「抬高 D 所在行」实现，不能靠字号/块宽** —— `--front-node-scale` 只动块宽。
故本轮改的是**主图行权与画布高**（`make_discovery_sci_figs.py`）：

| 项 | v21 | **v21.1（现行）** |
|---|---|---|
| `LAYOUT_ROW_RATIOS` | `[1.28, 1.00]` | **`[0.95, 1.00]`** |
| `COMPOSITE_FIGSIZE` | `(8.55, 10.0)` | **`(8.55, 11.5)`** |
| 上排（A/J）实测高 | 4.40 in | **4.32 in**（几乎不变） |
| 下排（D）实测高 | 3.43 in | **4.55 in（+33 %）** |
| D 面板内 Known 块 | 220 px | **291 px（+32 %）** |
| 资产内 Known 块 | 388 px | **678 px** |
| 资产尺寸 / 比 | 3130×1412 / 2.217 | **3130×1868 / 1.676** |

⚠ **行权刻意压到 `0.95`（< 1.00）**：这样 D 行增幅远大于上排。
若沿用 `1.28`，要让 D 到 4.55 in 就得把画布撑到 **17.4 in**（上排跟着涨到 6.9 in，纯浪费）。

⚠ **资产比须重算并钉到新 D 面板比**：D 4.55 in / 宽 7.61 in → 比 **1.673**；
反推 `--embed-size` 的画布高 = `tight_w_target / 目标比 / 0.9626` → **6.48 in**
（一次命中，实测 1.676 vs 1.673）。

⚠ 主图输出从 9.28×9.82 in 变为 **9.28×11.3 in**（高 +15 %）。
⚠ 独立图 `Fig_GojiVirome_NoveltySankey` **不受本轮影响**（它不走 `--embed-size`，md5 未变）。

**验收（实测）**：

| 项 | 旧 | 新 |
|---|---|---|
| 资产尺寸 | 3130×927（比 3.376） | **3130×1868（比 1.676）** |
| Known 块高 | 83 px | **388 px（4.7×）** |
| Novel 块高 | 43 px | **203 px（4.7×）** |
| 主图 D 内 Known 块 | ~83 px | **220 px** |
| 资产比 vs D 面板比 | 3.376 vs 2.216（**差 34 %**） | **2.217 vs 2.216（差 0 %）** |
| 独立图体检 | — | `_check_novelty_sankey_overlap.py`：**0 对重叠 / 无越界** |
| 主图体检 | — | `_check_discovery_overlap.py`：**20 figures, 0 overlapping pairs** |

⚠⚠ **资产宽高比必须钉到 D 面板比（2.216）**：主图贴图用 `imshow(aspect="auto")` 会强行拉伸，
比不匹配 = 桑基被压扁（旧版差 34 % 时是明显的纵向挤压）。反推方法 = 先渲一版量 tight 后的
实际比，再按 `画布高_new = 画布高_old × (目标 tight 高 / 实测 tight 高)` 迭代一次。

⚠ **本轮（第 1 轮，4.89 in）同时改了独立图 `Fig_GojiVirome_NoveltySankey`**（同一 `draw()` 被
两条路径共用）—— 独立图 9.32×8.4 in 竖直空间充足，放大后视觉更清晰，体检 0 对。
**若只想改主图 D 而不动独立图，需给两条 draw 调用传不同的 `level_fs`/`node_fs_scale`**（当前实现是共用的）。
⚠ **第 2 轮（6.48 in）只动主图**（行权 + 画布高 + `--embed-size`），独立图 md5 **未变**。

### ⭐⭐⭐ 2026-10-13 21:32「D图都看不清字体大些」→ v21.2（**先量后判**）

**⚠⚠⚠ 本轮查清了「D 面板字号物理硬上限」——这是本项目最容易反复踩的认知陷阱。**

#### 第一步：量化（不要靠感觉）

| 量 | 实测值 |
|---|---|
| `panel_D.png` 画布 | 2290×1035 px = 7.63×3.45 in |
| **桑基带**（剔标题/注记后） | 6.900 × **2.480 in**（= D 行高 × **0.545**） |
| 面板内 "DNA" 字形高 | **6 px** → 等效字号 **~2.0 pt** |
| 资产 3130×1868（10.43×6.23 in）贴入 | 水平缩放 **0.729** → 5.4 pt 被压成 **3.94 pt** |

#### 第二步：推公式（唯一瓶颈 = L5 行数）

因「**块高 ∝ 单位数**」（图注的严格比例声明）且 gap 是**物理量**，L5 的 **37 行**独占带宽：

```
字号上限[pt] = 桑基带高[in] × 300 / L5行数 / 1.30
```
本轮 `2.48×300/37/1.30 = `**`3.71 pt`**（实测 2.0 pt，还叠加各向异性缩放）。

#### ⚠⚠⚠ 第三步：三条「直觉修法」全部实测否掉（写死在代码注释里防复发）

| 死路 | 为什么无效 | 实测 |
|---|---|---|
| **① 加宽画布** | 比值约束 `tightW = tightH × 面板比` 是**死的**；tight 高饱和在 ~4.8–5.0 in → 加宽 W **不改缩放** | 8.8→16.0 in，`s` 0.905 → **0.441** |
| **② 加高画布** | 反使 `tightH` 涨 → `s` 更小 | 15.0×10.20 → 等效 **1.65 pt（最差）** |
| **③ 抬 `--all-fs`** | 抬字号**同时抬 `Y`**（gap 地板）→ 单位像素更小，**互相拆台** | 资产 6.5 pt 贴入只等效 **2.0 pt** |

⚠⚠ **关键澄清**：D 面板**桑基带**的真实比是 **2.78**（矮而宽），
不是 v21.1 算法里用的 1.673；而资产自然比只到 **~2.58** → 贴图必然有拉伸。
`eff = 资产字号 × 面板宽 / (画布高 × 常数)` —— **必须解方程，不能试参数**。

#### ✅ 第四步：用户拍板「只抬 D 行（有限改善）」

```python
# make_discovery_sci_figs.py
LAYOUT_ROW_RATIOS = [0.75, 1.00]     # 原 [0.95, 1.00]
COMPOSITE_FIGSIZE = (8.55, 13.0)     # 原 (8.55, 11.5)
```
⚠⚠ **`build_panels` 里 D 的单面板画布是硬编码 `7.61×2.98`** → 本轮**同步改 `7.61×3.35`**。
**这是本项目一个隐藏耦合点：单面板尺寸与组合图行高完全分离，改行权必须 `grep` 同步。**

```python
# fig_novelty_sankey.py 新增
--all-fs <pt>    # 覆盖**含 L5** 的全部五层（--front-fs 只管 L1–L4）
```
⚠ `solve_geometry` 调用点改用 `_geo_fs = max(FS_NODE, front_fs, all_fs)` ——
**抬 L5 字号必须同步抬行距地板**，否则 37 个标签互压（**零报错**）。
⚠ 资产参数**维持 v21.1**（`10.65x6.48 --front-fs 5.4 --front-node-scale 1.5`）：本轮不动资产。

#### 验收（实测）

| 项 | v21.1 | **v21.2** |
|---|---|---|
| 行权 / 画布高 | `[0.95,1.00]` / 11.5 | **`[0.75,1.00]` / 13.0** |
| 上排 A/J | 4.32 in | **3.94 in** |
| D 行 | 4.55 in | **5.25 in（+15 %）** |
| 主图 | 9.28×11.3 in | **9.28×12.78 in** |
| D 内 Known 块 | 296 px | **442 px（+49 %）** |
| 主图 D 内字形高 | **6 px** | **14–18 px（2.3–3×）** |
| 单面板桑基带 | 2.480 in | **2.650 in** |
| 主图体检 | 20 fig / 0 重叠 | **20 fig / 0 重叠** |
| 桑基体检 | 0 重叠 / 无脱画布 | **0 重叠 / 无脱画布** |

⚠⚠ **仍未到 6 pt 可读线**（上限 3.97 pt）：要 6 pt 需桑基带 ≥ **4.01 in**，
而 D 行即便 6.18 in（figH=14）也只 3.37 in ⇒ **在当前五层结构下不可达**。
⇒ 真正解法 = **减层数**（去掉 species → 最大 25 行）/ **只标重点 species** /
**主图只放缩略 + 明细看独立图**（已作备选告知用户）。
⚠ 本轮交付 md5：`Fig_GojiVirome_Discovery.png` = **`4406a34b0e277c206193e8e459c7f9ad`**。

### ⭐⭐⭐ 2026-10-13 22:28「字体再大些还是看不清」→ v21.3（**推翻并纠正 v21.2 的分析**）

#### ⚠⚠⚠ 先认错：v21.2 的「硬上限 3.97 pt」是错的

v21.2 里我把 **`panel_D.png` 上量到的「最大连续墨迹段 2.48 in」当成了「桑基带高」**，
于是算出「字号上限 3.71 pt、6 pt 不可达」。**这是错的。**

**正确的量法**：桑基的 `Y` 是**铺满整个 D 行**的，不是只占 2.48 in。判据 =
用 **Known 块高 ÷ 它在 `Y` 中的占比** 反推 `Y`：

```
Y_in_panel = Known块高(in) / (40 / Y_units)      # Y_units = tot/(1-36c)
```
v21.2 实测 Known = 1.473 in、`Y_units = 141` → **`Y_in_panel = 5.19 in ≈ D 行高 5.249 in`** ✅
（我原来量的 2.48 in 只是「**墨迹最长的连续段**」，被 L5 的稀疏标签打断了。）

⇒ **真实上限 = `Y_in_panel × 72 / 37 / 1.30`** = `5.19×72/37/1.30 = `**`7.77 pt`**（不是 3.97 pt）。

#### ✅ 正确的关系链（三条杠杆，可复用）

```
Y_in_asset  = 0.927 × canvas_H          # ylim span 1.03Y 映射到 0.955×canvas_H
Y_in_panel  ≈ D_row                     # 桑基铺满整行
paste_scale = Y_in_panel / Y_in_asset
panel_font  = fs × paste_scale = fs × D_row / (0.927 × canvas_H)
Y_units     = tot / (1 − 36c),  c = fs × 0.01890 / canvas_H
```
⇒ **两个独立杠杆**：① **抬 `D_row`**（线性抬字号）；② **抬 `fs/canvas_H` 比**（抬字号，但
`Y_units` 变大 → 块变薄）。`fs/canvas_H → 1.468` 时 `Y_units → ∞`（几何极限，块消失）。

#### 本轮三处改动（v21.3）

| 项 | v21.2 | **v21.3（现行）** |
|---|---|---|
| `LAYOUT_ROW_RATIOS` | `[0.75, 1.00]` | **`[0.6, 1.0]`** |
| `COMPOSITE_FIGSIZE` | `(8.55, 13.0)` | **`(8.55, 14.0)`** |
| `hspace` | `0.42` | **`0.28`** |
| 资产 `--embed-size` | `10.65x6.48` | **`6.25x5.20`**（**大幅缩短**） |
| 资产 `--all-fs` | 未用（`--front-fs 5.4`） | **`5.52`（全五层）** |
| 单面板 D 画布 | `7.61×3.35` | **`7.61×6.87`**（axes 6.18 in 对齐行高） |

⚠⚠ **`hspace` 0.42 → 0.28 是「白捡的高度」**：0.42 × 平均行高 = **2.08 in 的空带**
（图上肉眼可见大片留白）→ 压到 0.28 后同页高下 D 行 6.18 → **6.56 in**、上排 3.71 → **3.94 in**
（**两排都变大**）。⚠ 判据 = 上排 x 轴标签与 D 标题不得相撞（`_check_discovery_overlap` 0 对）。

⚠⚠⚠ **资产画布要「缩短」而不是「加高」** —— 这与 v21.1/v21.2 的直觉相反：
`panel_font = fs × D_row/(0.927 × canvas_H)`，**`canvas_H` 在分母**。前两轮把 canvas 从
3.10 一路加高到 6.48 in，其实是在**把字号往下压**（6.48 → 5.20 让字号涨 25 %）。

#### ⚠⚠ 顺带抓到的一个静默坑：`--all-fs` 漏了独立图

`--front-fs`/`--front-node-scale` 在代码里**同时作用于独立图**（`_lfs`/`_nsc`，L517–524）。
本轮只传 `--all-fs`（不传 `--front-fs`）→ 独立图的 L1–L4 **回落到 `FS_NODE` 4.05**，
比上一轮的 5.4 还小 —— **即「给主图放大字号」反而让独立图变小**（md5 `7cdd6c9c` → `68c16fb9`）。
→ 已修：`--all-fs` 现在同样作用于独立图的 `_lfs`/`_nsc`/`_base_fs`。
⚠ **一般化：凡是「多产出共用同一 `draw()`」的脚本，新增的逐层参数必须逐一确认覆盖到每一条调用路径**；
判据 = **改完比两条产物的 md5 + 量各自字形像素高**（只看主图会漏）。

#### 验收（实测）

| 项 | v21.2 | **v21.3** |
|---|---|---|
| 上排 A/J | 3.94 in | **3.94 in（不变）** |
| D 行 | 5.25 in | **6.56 in（+25 %）** |
| 主图 | 9.28×12.78 in | **9.28×13.80 in** |
| **主图 D 内 L5 字形高** | **6 px（≈2.0 pt）** | **24 px（cap 5.76 pt → 字号≈8.23 pt）** |
| 主图 D 内 Known 块 | 442 px | **364 px** |
| 资产 tight | 10.43×6.23 in | **6.55×5.00 in** |
| 独立图 L5 字号 | 4.05 pt | **7.54 pt** |
| 主图体检 | 20 fig / 0 重叠 | **20 fig / 0 重叠** |
| 桑基体检 | 0 重叠 / 无脱画布 | **0 重叠 / 无脱画布** |

⚠ 字号 **2.0 → 8.23 pt（4.1×）**，代价 = Known 块 442 → 364 px（**−18 %**）。
交付 md5：主图 **`8245fea2d1d782ee34cbfccabf7c587b`**、独立图 **`a9002e8f424033db0e60e61e2e5c2f32`**。

## UpSet 系列（R / ComplexUpset）— 工具重叠 · 分类证据 · 宿主预测

2026-10-13 新建。用 **R + ComplexUpset** 出三张 UpSet 图（与 Python 主图并存，不改动主图链路）。
共用样式与配色在 `upset_common.R`（`OKABE_ITO` 顺序与 `sci_palette.py` 同源）。

| 脚本 | 产出 | 数据源 | 规模 |
|------|------|--------|------|
| `fig_upset_tools.R` | `Fig_GojiVirome_UpsetTools.{png,pdf}` | `05_Taxonomy/Votus.integrated/intersection_matrix.tsv` | **21,630** vOTU × **7** 工具；111 个非空交集，显示前 18 |
| `make_agree_matrix.py` | `upset_data/agree_matrix_<Rank>.tsv` | `10_Reports/final_integrated_classification.tsv` 的 `<Rank>_agree` 列 | 逐 contig 赞成票矩阵（见下「阶元参数化」） |
| `fig_upset_classification.R` | `Fig_GojiVirome_UpsetClassification[_<Rank>].{png,pdf}` | 同上 | **Species** 17,282 vOTU（49 交集，显示前 **16**）；**Family** 17,512 vOTU（77 交集，显示前 **18**） |
| `make_cls_tool_timing_fig.R` | `upset_data/cls_tool_timing_fig[_<Rank>].png`（**预渲，供上图左上角嵌入**） | `05_Taxonomy/Votus.classed/Votus_resource_usage.tsv` | 7 分类工具**全库**耗时（15.8–786.0 s） |
| `fig_host_prediction.R` | `Fig_GojiVirome_HostPrediction.{png,pdf}`（**4 面板**） | `06_HostPrediction/` | 20,892 vOTU；预测器重叠 **20,443** 条（≥1 命中） |

#### ⭐ 阶元参数化（2026-10-13 追加）

**用户指令**：「图B 分类的一致性」→ 澄清确认要**科级（Family）**的赞成票一致性。

- `make_agree_matrix.py --rank <Rank>`：从 `<Rank>_agree` 列（格式 `"3/4: ACVirus,CAT,diamond_lca"`）
  取**冒号后的工具集合**→ 映射成 7 列 0/1 矩阵。
  ⚠ 与图①**不是同一口径**：图① = 「谁检出了」（21,630 全集）；本矩阵 = 「谁在该 rank 投了赞成票」（子集）。
  ⭐ 对拍验证：`--rank Species` 重现 **17,282** 条 + 分布 `{1:13316, 2:2827, 3:844, 4:189, 5:99, 6:7}`，
  与既有 `agree_matrix_Species.tsv` **逐位一致**（零回归）。
- `fig_upset_classification.R --rank <Rank>`：默认仍 Species（零回归），
  `--rank Family` 出科级版，输出名自动加 `_Family` 后缀。

**族级 vs 种级（实测对比，写论文可引）**：

| 阶元 | 入图 vOTU | 无票行 | 交集数 | ≥3 票 | ≥4 票 | 最高票 |
|---|---|---|---|---|---|---|
| Species | 17,282 | 3,653 | 49 | 1,139 | 295 | 6 |
| **Family** | **17,512** | 3,423 | **77** | **4,414** | **1,931** | 6 |

→ **科级一致性显著高于种级**（≥3 票：4,414 vs 1,139，约 3.9×）：符合分类学规律
（科级可比对的特征更保守、工具更易达成共识）。

#### ⚠⚠ 本轮修掉的两个「零报错」坑

1. **左侧 set-size 面板 y 轴标签变乱码**（实测显示 `6 1 4 2 7 5 3`）：
   工具名含下划线（`diamond_lca`）不是合法 R 名 → ComplexUpset 默认 `encode_sets = TRUE`
   把集合**编码成整数 id**。修法 = **`encode_sets = FALSE`**（与 `fig_ident_filter_stacked_all.R` 同一坑）。
2. **轴标题显示为 `group`**：ComplexUpset 的 `upset_set_size()` 默认 x 轴标题就是 `group`（通用措辞）。
   修法 = 显式 `xlab("Set size")` + `theme(axis.title.y = element_blank())`
   （集合名已由矩阵面板统一显示，不重复）。

#### ⚠⚠⚠ 柱顶数字标到柱子里面（用户 2026-10-13「第一条柱子的数字标在了柱子里面」）

**现象**：第一根柱（全图最高）的数字压在柱内（白字压深色柱），其余柱正常。

**根因（`ggplot_build()` 实测读回，Family 版）**：ComplexUpset 的
`intersection_size(counts = TRUE)` 对「**全图最大值那根柱**」特殊处理 ——
把文字放在 **y = 4,626.55 = 柱高 5,443 的 85 %** 处（即柱内），
而其余柱（2,711 / 1,446 / …）都用**真实柱高**。

⚠ **试过且无效**：把 `scale_y_continuous` 的 `expand` 从 `0.16` 提到 `0.26`
（顶部净空 82 px → 210 px）—— 文字仍在柱内，因为偏移是**柱高比例**不是像素量。

**正解** = **不用内置文字（`counts = FALSE`），自己 `annotate()` 按真实柱高标注**：
```r
intersection_size(counts = FALSE, mapping = aes(fill = as.factor(count))) +
annotate("text", x = seq_len(nrow(sig_top)), y = sig_top$n,
         label = format(sig_top$n, big.mark = ","), vjust = -0.5, size = 2.4)
```
副作用（顺带改进）：所有数字统一在柱顶上方（深色）、加千分位逗号。

**⚠⚠ 自绘带来的新约束：`--top` 不能随便设**

自绘数字依赖「**前 N 个交集大小互不相同**」—— 有并列时 ComplexUpset 的排序先后
与 `sort(sizes)` 的表顺序可能不一致 → **数字整体错位，零报错**。
脚本已加自检（并列即 `stop()`），实测无并列的安全值：

| 阶元 | 安全 top | 首个并列在哪 |
|---|---|---|
| **Species** | **16** | 前 17 起（68 并列） |
| **Family** | **18** | 前 20 内无 |

⚠ 历史：Species **曾用 18** —— 旧版靠 ComplexUpset 内置 `counts=TRUE`，它用**自己的内部排序**
所以不会错位；改自绘后必须降到 **16**。

```bash
cd <project root>
Rscript scripts/barbarum/fig_upset_tools.R                       # 图①
python  scripts/barbarum/make_agree_matrix.py --rank Family       # 先出矩阵（Species 同理）
Rscript scripts/barbarum/make_cls_tool_timing_fig.R               # 图② 左上角耗时图（预渲）
Rscript scripts/barbarum/fig_upset_classification.R               # 图② Species（默认）
# 换阶元时：矩阵与耗时图都要按该阶元的行序重出
python  scripts/barbarum/make_agree_matrix.py --rank Family
Rscript scripts/barbarum/make_cls_tool_timing_fig.R --rank Family \
        --out scripts/barbarum/upset_data/cls_tool_timing_fig_Family.png
Rscript scripts/barbarum/fig_upset_classification.R --rank Family \
        --timing-png scripts/barbarum/upset_data/cls_tool_timing_fig_Family.png
Rscript scripts/barbarum/fig_host_prediction.R                   # 图③（4 面板）
# 依赖：Rscript 4.6 + ComplexUpset 1.3.3 + patchwork + data.table + ggplot2 + png + grid
```

#### ⭐ 图② 左上角「工具耗时」条形图（用户 2026-10-13「时间的图呢，左上角没看到呀」）

**做法** = IdentFilter 同款「预渲 PNG → rasterGrob → gtable+grid 叠」（唯一零 reflow 的画法）。

**框位置（本版干净基线逐像素实测，画布 2796×1800 px @300dpi）**：

| 元素 | 实测 y（px） |
|---|---|
| 标题 | 26..66 |
| 副标题 | 91..118 |
| **可用空区** | **119..527**（x 0..880） |
| set-size 面板刻度左端 | 528 起 |

→ 实际取框 **y 137..513**（上留 18 px、下留 14 px 净空）= **880 × 376 px = 2.9333 × 1.2533 in**。
⚠⚠ **必须用「本次改动后」的干净基线重测**，禁沿用上一版坐标（README 硬规则）。

**⚠⚠ 口径差异（图注必须声明）**：

| 图 | 口径 | 工具集 |
|---|---|---|
| 本图（图② 左上） | **一次性跑全库 21,630 vOTU 的总 wall time** | 7 个**分类**工具 |
| IdentFilter 的 `tool_timing.tsv` | **单样本平均**（395 样本） | 7 个**鉴定**工具 |

两者**只有 metabuli / genomad 两个工具重叠**，**不可比、不可混画**。

**⭐ 行序自动对齐**：`make_cls_tool_timing_fig.R` **不硬编码顺序** —— 它读同一份
`agree_matrix_<Rank>.tsv`、按各工具 set size 降序算出 `ORDER`，与 ComplexUpset 点阵行序同源。
（换阶元后行序会变：Species = metabuli > ACVirus > diamond_lca > CAT > mmseqs > VITAP > genomad；
Family = … > VITAP > genomad 末两位互换 → 脚本自动跟随。）

⚠ 踩过的两个「改行序」静默坑（与 IdentFilter 同源）：
① `as.numeric(factor)` 给的是 levels 序号、**不是轴标签序** → 必须显式 `ypos` 连续坐标；
② `scale_y_continuous` 的 y **向上增大** → `ypos` 要**反相**（`N+1-match`）才能让第 1 名在最上。

**⚠⚠ 2026-10-13 二次微调（用户「时间的图 位置往下些，往右边些」）**：

框从 `x 0..880 / y 137..513`（2.9333 × 1.2533 in）改为
**`x 100..810 / y 162..525`（2.3667 × 1.2100 in，即 `TIM_L=100/2796, TIM_R=810/2796, TIM_T=1-162/1800, TIM_B=1-525/1800`）**。

**为什么是这个框（逐像素实测的「障碍物」清单）**：右下移的约束**只有两处**，
都在 `x 815..1009`：

| 障碍 | 位置 | 性质 |
|---|---|---|
| 刻度数字 `6000` | y 284..305 | 柱面板 y 轴刻度 |
| 刻度数字 `4000` | y 527..541 | 柱面板 y 轴刻度 |

→ 框右界收到 **810**（离最近障碍 5 px，且障碍在**下方**区域，垂直不重叠）；
框底压到 **525**（`4000` 刻度在 527 起，留 2 px）。**最终四邻域实测全 0 墨迹**：

| 邻域 | 墨迹 |
|---|---|
| 上（y110..160，副标题底 y112） | 1115（= 副标题自身，间隙 **49 px**） |
| 下（y526..580） | **0** |
| 左（x40..99） | **0** |
| 右（x811..870） | **0** |

⚠⚠ **改框尺寸后必须重跑 `make_cls_tool_timing_fig.R`**（预渲 PNG 不随参数更新，会静默嵌旧尺寸图）；
且该脚本的 **`--w`/`--h` 默认值已同步改为 `2.3667`/`1.2100`**，避免下次不带参数跑时尺寸错配。

⚠⚠⚠ **`--out` 默认名已改为「按 `--rank` 推导」**：Species → `cls_tool_timing_fig.png`（无后缀，向后兼容）；
其他阶元 → `cls_tool_timing_fig_<Rank>.png`。
**修这个坑的原因**：原先两个 rank 都默认写同一个 `cls_tool_timing_fig.png`，
不带 `--out` 依次跑 Species + Family 会**后者静默覆盖前者**（实测 `--rank Family` 仍写无后缀名，零报错）→
随后 `fig_upset_classification.R --rank Species` 会嵌到 **Family 行序的错误预渲图**。

#### ⚠⚠⚠ ComplexUpset：改「交集柱面板」的 theme **必须写在 `base_annotations` 里**（踩坑 3 次）

本轮为「去除 upset 图 Y 轴的 Agreed vOTUs」连续试了 3 条路，**前两条全部静默失效、零报错**：

| # | 尝试 | 结果 | 原因 |
|---|---|---|---|
| ① | 删掉 `ylab("Agreed vOTUs")` | ✗ 无效 | 该标题**不由 `ylab()` 产出**，是 `upset()` 内部 L176 `annotation_plot + ylab(name)` 补写的 |
| ② | `themes = upset_modify_themes(list("intersections" = ...))` | ✗ 静默失效 | **key 名不匹配任何面板** —— `upset_modify_themes()` 用 `setdiff()` 对未知 key 只是「附加」，**不报错** |
| ③ | `themes = upset_modify_themes(list("Agreed vOTUs" = ...))` | ✗ 仍无效 | key 对了（探针实测该 theme 确实挂到 `patches$plots[[2]]`），但**应用顺序**输给 `base_annotations` |
| ✅ | 写在 **`base_annotations` 的 list 元素 theme 里** | ✓ 生效 | 见下 |

**根因（`upset()` L169-170 的应用顺序）**：

```r
annotation_plot = annotation_plot + selected_theme + do.call(theme, user_theme)
#                                  ↑ themes= 来的        ↑ base_annotations 里跟的 theme_upset(8)
```

`user_theme` 排在**后面** → 它把 `axis.title.y` 从 `element_blank()` 覆盖回 `element_text`。
**所以想在交集柱面板上关掉任何轴元素，必须写在 `base_annotations` 那一侧的 theme 里**（最后应用，必胜）。

**内置 theme key 全集**（`ComplexUpset:::upset_themes`，可 `Rscript` 读出来核对）：
`intersections_matrix` / **`Intersection size`** / `overall_sizes` / `default`
—— ⚠ 注意 `"Intersection size"`（**大写 + 空格**）才是**默认**交集柱面板名；
我们自定义了 list 名（`"Agreed vOTUs"`），故 key 也要跟着自定义名走。

**验证手法（可复用）**：直接枚举运行时的 patchwork 子图，看 theme 到底有没有挂上 ——

```r
p <- upset(...)
for (i in seq_along(p$patches$plots)) {
  pl <- p$patches$plots[[i]]
  if (inherits(pl, "ggplot"))
    cat(i, deparse(pl$labels$y), class(pl$theme$axis.title.y), "\n")
}
```

⚠⚠ **判「轴标题是否真去掉」必须量像素**：本例成品图上文字位于 **x 815..847 / y 527..743**
（宽 33、高 217 的竖排块），改完判据 = **该区域墨迹 = 0**。只看"出图成功"会漏判。

**口径（判据优先，勿混）**：

- **图①「7 集合」** = `diamond_lca / mmseqs / CAT / metabuli / genomad / VITAP / ACVirus`
  （**不是** 10 个 —— R 脚本 `virus_classifier_analysis.R` 支持最多 10 个工具入口，
  但本项目 05_Taxonomy 实际只跑到 **7 个**；判据 = `Votus.integrated/standardized_*.tsv` 有 7 个）。
  degree 分布：1→9,220 / 2→2,266 / 3→2,473 / 4→2,735 / 5→2,704 / 6→1,656 / **7（全票）→576**。
- **图②** 取 `<Rank>_agree` 列的**赞成票集合**（`"3/4: ACVirus,CAT,diamond_lca"` → 括号内工具）。
  本图用 **Species rank**。⚠ 与图①**不是同一口径**：图① = "谁检出了这条 contig"，
  图② = "谁在物种级投了赞成票"（是子集，故 17,282 < 21,630）。
- **图③ 四个预测器** = `ICTV`（`classify_contigs.py` 查库，C9）/ `RVH`（`rnavirhost predict`）/
  `PB2`（`phabox2 --task cherry`）/ `Rule`（**分类级回退**，`Decision_Method` 含 `Rule_Class_Taxonomy`）。
  ⚠ **`Rule` 不是独立序列比对** → 面板 b/c/d 与整图副标题均已标注。
  命中集合大小：ICTV 19,891 / RVH 8,394 / PB2 3,236 / Rule 1,773；
  degree：1→12,077 / 2→7,012 / 3→1,305 / **4（四路全票）→49**。

**⚠ 本系列 5 条硬规则（R / ComplexUpset 特有）**：

1. ⚠⚠ **复杂 UpSet 必须限量**：7 集合最多 **127** 个非空交集，全画 x 轴标签会糊成一片
   （实测首版 111 个交集 + 自动转置 → 标题被挤到图中央）。→ 用 `min_size = sizes[top_n]`
   截断，**并在副标题如实写 "showing top N of M intersections"**。
2. ⚠⚠ **`upset()` 返回 patchwork，不能直接 `+ labs(title=)`**：`labs()` 会落到嵌套子图、
   标题跑进图中央（实测图① 首版）。→ 用 `patchwork::plot_annotation(title=, subtitle=)`；
   若是**嵌套**（UpSet 再进 `patchwork` 布局），再套一层 `wrap_elements(textGrob(...))`
   作独立题头条（图③ 面板 b 用的就是这招）。
3. ⚠ **Rscript 不带 `sys.frame(1)$ofile`**：`source()` 兄弟脚本要自己从 `commandArgs` 解析
   `--file=`。⚠⚠ 且在**含中文的路径**下 `normalizePath()` 会把路径编坏
   （`D:/f!i"/e;6d<8e:e g;/...`）→ **不要 `normalizePath`**，直接用 `dirname(file)`。
4. ⚠ **R 的 `sprintf` 不支持 `%,d`**（那是 Python 语法）→ 千分位用 `format(x, big.mark=",")`。
5. ⚠ **ComplexUpset 的 `min_size` 是"交集大小阈值"不是"交集个数"** → 要显示前 N 个交集，
   得先自己 `sort(table(sig), decreasing=TRUE)` 取第 N 名的 size 当阈值。

## fig_ident_filter.R — 「鉴定 / 过滤」双簇分组条形图

2026-10-13 新建。用户指令「鉴定和过滤的分组条形图」→ 二次确认选定
**「鉴定与过滤双簇」+「过滤只展示保留总数」**。与 UpSet 系列共用 `upset_common.R` 的配色与主题。

| 项 | 内容 |
|----|------|
| 脚本 | `fig_ident_filter.R`（绘图）+ `make_identfilter_data.py`（数据中间件） |
| 产出 | `<S10_SCI_Discovery>/Fig_GojiVirome_IdentFilter.{png,pdf}`（**9.32 × 5.4 in 通栏**） |
| 数据源 | `10_Reports/ident_summary.tsv`（TOTAL 行）+ `10_Reports/filter_per_sample.tsv` |
| 中间件 | `identfilter_data/ident_filter_summary.tsv` + `ident_filter_meta.json` |

```bash
cd <project root>
C:/Python312/python.exe scripts/barbarum/make_identfilter_data.py   # 汇总（含口径自检）
Rscript scripts/barbarum/fig_ident_filter.R                          # 出图
```

**作图数值**（y 轴 = **log10**，线性轴下 604/2,229/3,219 只有 0.1–0.6% 高度、肉眼不可见）：

| 簇 | 条目 | 数值 |
|----|------|------|
| 鉴定 | viralverify / metabuli / blast / genomad / virbot / rdrpcatch | 452,904 / 144,041 / 95,573 / 11,526 / 3,219 / 2,229 |
| 鉴定 | **viroid_blastn**（`*`，**独立类病毒轨道**） | 604 |
| 过滤 | Retained (02b_Filter) | **94,826 = 并集的 17.0%** |
| 参考线 | 鉴定并集 `All_Candidate`（= 过滤输入） | **556,773** |

**口径（判据优先，勿混）**：

1. ⚠⚠ **`viroid_blastn` 不在 556,773 并集里**。2026-10-13 在 `02a_Identification/<样本>/`
   用 `comm -23 <(sort viroids.result.id) <(sort virus.all.result.id)` 逐样本验证：
   viroid id **全部**落在 `virus.all.result.id` 之外（CRR1126134 = 3/3、SRR23107142 = 5/5）。
   全库 604 条、**245 个样本非零**。→ 图上给 `*` 标记 + 副标题显式声明，**不与病毒工具同色**。
2. ⚠⚠ **「10 个工具」是真的 10 个**（`report_pipeline.py:385` 的 `tools_list`）：
   `genomad / blast / metabuli / virsorter2 / viralverify / virhunter / virbot / viralm /
   rdrpcatch / viroid_blastn`。其中 **3 个零命中**（`virsorter2` / `virhunter` / `viralm`）——
   在 log 轴上无法画 0，故只在副标题里声明。
3. ⚠ **各工具之和 ≠ 并集**：Σ = **709,492**，并集 = **556,773**，重叠 **152,719**。
   → 不要拿 Σ 当"鉴定总数"。
4. ⚠ **过滤侧只画 Uniprot 保留总数**（用户选定）。`filter_per_sample.tsv` 里的
   `Cdd_Filtered`（473,523）**不在图上**；且 ⚠ `filter_summary.tsv` 是 **396 行全 0 的空表**，禁用。
5. ⚠ **两个样本集不同**：鉴定 438 样本 vs 过滤 395 样本（副标题只写鉴定侧的 438）。

**⚠ 本图 5 条硬规则（ggplot2 特有，见报告 §9 规则 42–46）**：

1. ⚠⚠ **`geom_col()` 在 `scale_y_log10()` 下会把柱子整批静默移除**：`geom_col` 的
   `ymin = 0` → `log10(0) = -Inf` → 只剩一句 `Removed N rows ... (geom_col())` 的 warning，
   **图上什么柱子都没有**。→ 改用 `geom_rect(ymin = <轴下限>)`。
2. ⚠⚠ **图内文本必须全 ASCII**：副标题里的 `U+00B7 (·)` / `U+2192 (→)` 在 Windows 上触发
   `mbcsToSbcs` 字形转换失败（`dot substituted for <e2>`），**静默渲染成点**。→ 用 `;` 和 `->`。
3. ⚠⚠ **数值标签必须放进柱内**：柱顶 452,904 与并集线 556,773 在对数轴上只差 **0.09 decade**，
   标签放柱外会被并集虚线**横穿**（首版实测碰撞）。→ `vjust = 1.55` 放柱内，
   白字配蓝/绿柱、黑字配橙柱。
4. ⚠ **参考线的文字用 `vjust = 0` 让底边落在该 y 上**，否则线同样会横穿文字。
5. ⚠⚠ **读服务器汇总表前先 `diff` 表头**：本地 `ident_summary.tsv` 是 **10 列**旧副本
   （缺 `viroid_blastn`），服务器是 **11 列** → 直接用会静默漏掉整条类病毒轨道。
   消费端已加表头自检（缺列即 `sys.exit`）。

## fig_ident_filter_upset.R — 「7 鉴定轨道 × 过滤结局」堆叠 UpSet

2026-10-13 新建。用户指令「堆叠的 UpSet 图，展示 7 工具的鉴定和过滤」→ 二次确认选定
**「鉴定阶段 7 轨道」×「柱内分段」**（每根交集柱按过滤结局堆叠）。与 UpSet 系列共用 `upset_common.R`。

| 项 | 内容 |
|----|------|
| 脚本 | `fig_ident_filter_upset.R`（绘图）+ `make_ident_filter_upset_data.py`（数据中间件，**在服务器上跑**） |
| 产出 | `<S10_SCI_Discovery>/Fig_GojiVirome_UpsetIdentFilter.{png,pdf}`（**9.32 × 6.8 in 通栏**） |
| 数据源 | `02a_Identification/<S>/<S>_virus.{tool}.result.id` + `<S>_viroids.result.id` + `02b_Filter/<S>/{cdd,uniprot}/*.fasta` |
| 中间件 | `identfilter_data/ident_filter_upset_signature.tsv` + `_meta.json`（**64** 个签名 × 状态计数，仅 1.6 KB） |
| 集合（7） | viralverify / metabuli / blast / genomad / virbot / rdrpcatch + **viroid_blastn** |

```bash
# 1) 服务器端汇总（读 395 个样本的逐序列产物，约 1 分钟）
ssh <246> 'export LD_LIBRARY_PATH=$HOME/mambaforge/lib; ~/mambaforge/bin/python ~/make_ident_filter_upset_data.py'
scp <246>:~/ident_filter_upset/ident_filter_upset_* scripts/barbarum/identfilter_data/
# 2) 本地出图（服务器同样可跑，已双端验证）
Rscript scripts/barbarum/fig_ident_filter_upset.R
```

**过滤链（实测，链式嵌套成立）**：`union 556,950 ⊇ CDD 保留 473,523 ⊇ Uniprot 保留 94,826`
（= union 的 **17.03%**）。入图 557,439 条（union ∪ viroid-only），状态**互斥**：
Retained 94,818 / Removed by UniProt 378,669 / Removed by CDD 83,436 / Viroid track 516 /
无 02b 产物 27（**剔除不入图**，副标题已声明）。

**图上数值**：显示前 **20 / 64** 个交集（覆盖 **99.58%** 序列）。柱顶 % = 该交集的**最终保留率**：
单工具 viralverify **13%** / metabuli 21% / blast 24% / genomad 50%；4 工具 70%；5 工具 96–99%；6 工具 **98%**。
→ **结论：工具共识越高，保留率越高**（单工具 13% vs 六工具全票 98%）。

**⚠⚠ 本次修正一条旧口径（重要）**：上面 `fig_ident_filter.R` 章节写的
「`viroid_blastn` 全部落在 `virus.all.result.id` 之外」**只在 2 个样本上抽验过**
（CRR1126134 = 3/3、SRR23107142 = 5/5）。本次扫**全部 395 个样本**：
**`viroid ∩ union = 88 ≠ 0`** → 两条轨道**确有**重叠（签名 `0100001`=74、`1100001`=13、`1000001`=1）。
图内已声明 "overlaps the viral tools in 88 contigs"。
⚠ 推论：状态判定必须**先判是否在 union 内**，只有 union **之外**的 viroid 命中才算独立轨道；
若把 viroid 判在前面，这 88 条的过滤结局会被静默吞掉。

**⚠ 本图 6 条硬规则**：

1. ⚠⚠ **`signature` 是 7 位 0/1 串，`fread` 会当整数读、吃掉前导零**（`"0100000"` → `100000`）
   → 必须 `colClasses = list(character = "signature")`，否则位序整体错位（首版即 `nchar` 校验失败）。
2. ⚠⚠ **ComplexUpset 会把 `intersect` 给的集合顺序上下翻转** → 要最大集在最上面，得传 `rev(TRACKS)`。
3. ⚠⚠ **集合标签必须单行**：两行标签（`"viralverify\n(n=452,904)"`）总高超过点阵行距 → 与上一行互相压字。
   改用单行 `sprintf("%s (%s)", ...)` + `set_sizes = FALSE`（ComplexUpset 默认那个左侧面板只显示
   **内部编号裸数字**，本系列其他图也有，但语义不清 → 本图去掉，把集合大小写进标签）。
4. ⚠⚠ **柱顶标注依赖交集顺序**：位置用序号 `1..N`，前提是**前 N 个交集大小互不相同**。
   脚本已加 `anyDuplicated()` 自检，一旦并列立即 `stop()`（否则 % 会整体错位）。
5. ⚠ **`intersection_size(mapping = aes(fill = status))` 才能出堆叠柱**；`counts = TRUE` 的计数标签
   在堆叠下位置会乱 → 用 `counts = FALSE`，自己 `annotate("text", x = 序号, y = 柱高, ...)`。
6. ⚠ **`set_sizes = FALSE` 后左上会留一块空白**（柱面板的 y 轴区）—— 正常，别去"填"它。

**未查清的小差异（0.03%，不影响图）**：id 文件口径 Σunion = **556,950**，而 `ident_summary.tsv`
的 `All_Candidate` TOTAL = **556,773**，差 **177**。逐样本比对指向 `CRR1126132`：
其 `_virus.all.result.id` 有 **177** 行，但 `ident_summary` 记 `All_Candidate = 0`
（同时该样本的 `_virus.all.candidate.fasta` 只有 2 字节 = 空）。→ 两个口径在个别样本上不一致，待查。

**双端复现**：本机与 246 各跑一次均成功（PNG 体积 132 KB / 344 KB —— 字形 hinting 差异，非内容差异）。

## fig_ident_filter_by_tool.R — 「按鉴定工具」过滤结局 100% 堆叠柱

2026-10-13 新建，是上面 `fig_ident_filter_upset.R` 的**伴随图**。用户指令「右侧的堆叠图呢，同样展示」
→ 二次确认选定 **「按工具的堆叠柱」**：把 x 轴由「交集」换成「鉴定工具」，沿用同一份数据、同一套配色。
（「右侧」指 UpSet 图右上那块按结局堆叠的交集柱面板。）

| 项 | 内容 |
|----|------|
| 脚本 | `fig_ident_filter_by_tool.R`（**纯绘图，无新数据中间件**） |
| 产出 | `<S10_SCI_Discovery>/Fig_GojiVirome_IdentFilterByTool.{png,pdf}`（**9.32 × 5.2 in 通栏**） |
| 数据源 | 与 UpSet 图**同一份** `identfilter_data/ident_filter_upset_signature.tsv`（不引入新口径） |
| 轴 | x = 7 条鉴定轨道（6 病毒按命中量降序 + 间隔位 + viroid_blastn）；y = **该工具命中的过滤结局构成（归一化 100%）** |

```bash
Rscript scripts/barbarum/fig_ident_filter_by_tool.R        # 本机
ssh <246> 'export LD_LIBRARY_PATH=$HOME/mambaforge/lib; cd ~/MMPV-RNA/scripts/barbarum && Rscript fig_ident_filter_by_tool.R'
```

**逐工具数值**（`n_hits` 含结局未知者；柱内四段不含）：

| 工具 | n_hits | Retained | Removed by UniProt | Removed by CDD | Viroid | 保留率 |
|------|-------:|---------:|-------------------:|---------------:|-------:|-------:|
| viralverify | 452,904 | 69,817 | 308,789 | 74,279 | 0 | **15.4%** |
| metabuli | 144,041 | 30,425 | 103,252 | 10,359 | 0 | 21.1% |
| blast | 95,573 | 26,811 | 68,543 | 211 | 0 | 28.1% |
| genomad | 11,526 | 6,918 | 2,580 | 2,028 | 0 | 60.0% |
| virbot | 3,219 | 2,682 | 510 | 27 | 0 | **83.3%** |
| rdrpcatch | 2,229 | 1,514 | 618 | 96 | 0 | 68.0% |
| viroid_blastn | 604 | 11 | 60 | 17 | 516 | 1.8% |

Σ(工具命中) = **710,096**（重叠计数，**不是**并集 556,950）。viroid_blastn 的 11+60+17 = **88**
即 `viroid ∩ union`，其余 516 条为 union 之外的类病毒专有命中。
→ **结论：工具越"专一"（命中越少），保留率越高**（viralverify 15.4% → virbot 83.3%）。

**⚠ 本图 4 条硬规则**：

1. ⚠⚠ **必须归一化到 100%，不能画绝对堆叠柱**：7 条轨道量级差 **750 倍**（452,904 vs 604）。
   同一线性轴上后 4 根只有 viralverify 的 0.1%–2.5% 高 → 肉眼不可见。绝对数写在工具名下方 + 柱顶标保留率，信息不丢。
2. ⚠⚠ **`scale_y_continuous(limits = c(0, 1))` 会把柱顶标注（y=1.035）静默变 NA 并丢掉**
   （`limits` 是**数据范围过滤器**，不是 `coord_cartesian` 的视窗裁剪）→ 首版 7 个 % 标签**全部消失**（仅一条 warning）。
   修法：**去掉 `limits`**，只留 `expand = expansion(mult = c(0, 0.085))`。
3. ⚠⚠ **柱内分段顺序**：ggplot2 因子**首档默认在顶** → 要 `Retained` 在**柱底**必须 `position_stack(reverse = TRUE)`。
4. ⚠⚠ **副标题必须折行**：单行 196 字符会被画布**右侧裁掉**（首版实测在 "share retained after" 处断掉）。
   改为 4 行、每行 ≤ 110 字符。**判据 = 量墨迹范围**：本图 ink bbox 左/右边距均 34 px、外侧 12 px 带内零墨。

**双端复现**：本机与 246 各跑一次均成功；同尺寸 2796 × 1560，ink bbox 差 ≤ 2 px、ink 像素数比 **1.0000**、
主导色**逐字节一致**（#009E73 / #BFBFBF / #7A7A7A / #E69F00）；差异仅 4–5.8% 像素、集中在**文字带**（字形 hinting）。

## fig_ident_filter_combined.R — 「堆叠 UpSet（左）+ 按工具堆叠柱（右）」双联同图

2026-10-13 新建。用户指令：「不应该和堆叠的 upset 图一起展示放到右侧嘛」
→ 把上面两张图**合成一张**：左 = UpSet 图体，右 = 按工具 100 % 堆叠柱，共用一套配色与**一个**图例。

| 项 | 内容 |
|----|------|
| 脚本 | `fig_ident_filter_combined.R`（自带数据装载 + 两个面板构建，**有意不改**那两张已验证的单图脚本，以免破坏其零回归；口径变动需三处同步） |
| 产出 | `<S10_SCI_Discovery>/Fig_GojiVirome_IdentFilterCombined.{png,pdf}`（**12.6 × 8.4 in**） |
| 数据源 | 与两张单图**同一份** `identfilter_data/ident_filter_upset_signature.tsv` |
| 版式 | `wrap_elements(title_grob) / (pA \| pB) / guide_area()` + `plot_layout(heights = c(0.15, 1, 0.05), widths = c(2.3, 1), guides = "collect")` |

```bash
Rscript scripts/barbarum/fig_ident_filter_combined.R     # 默认 --top 16（覆盖 99.39%）
```

**⚠ 本图 4 条硬规则**：

1. ⚠⚠ **图例只能出一份**：`patchwork` 的 `guides = "collect"` **进不到 `upset()` 这个嵌套 patchwork 里面**
   → 若两侧都留图例，会得到**两份一模一样**的图例（首版实测）。修法 = 在 UpSet 的 `base_annotations` 主题里
   显式 `legend.position = "none"`，只留右侧面板的图例，再由 `guide_area()` 收到底部居中。
2. ⚠⚠ **收尾不能用 `& theme(legend.position = …)`**：`&` 会把主题应用到**所有**子图（含 pA）→ 会把刚设的
   `"none"` 覆盖掉 → 图例又变两份。图例样式只能写在 pB 自己的 `theme()` 里。
3. ⚠ **「上排堆叠柱有多高」由 `height_ratio` 控制**（ComplexUpset 默认 **0.5**），**不是画布高**：
   只把画布从 7.2 → 8.8 in，柱面板**一点没变高** —— 矩阵被 `coord_equal` 锁死，多出的高度不会给柱面板。
   本图用 `height_ratio = 0.36`（0.30 会把矩阵行距压得太挤）。
4. ⚠ **`--top` 由 20 降到 16**：右侧面板占掉约 1/3 宽度后，20 根柱的柱顶 % 会挤在一起；
   16 根（覆盖 **99.39 %**）在 300 dpi 下标签清晰不叠。**判据 = 裁切 + 2× LANCZOS 放大目视**。

**双端复现**：本机与 246 各跑一次均成功；同尺寸 3780 × 2520，ink bbox 差 ≤ 2 px、主导色**逐字节一致**。
⚠ 4× box 降采样后仍有 **9.5 %** 像素差 > 40/255 —— 这是**字形 hinting 的基线，不是版式差异**：
同一测法下已验证的单面板图（ink 像素数比 1.0000）也有 **7.2 %**；本图文字量约为其 2 倍，故略高。

## fig_ident_filter_stacked_all.R — 「一体式堆叠 UpSet」：set-size 面板 + 交集柱**两侧都堆叠**

2026-10-13 新建。用户指令：「不能拼在一起嘛？类似这种，只是都是堆叠的」（附一张原生 ComplexUpset 参考图：
左侧 set-size 横条 + 右上堆叠交集柱 + 右下点阵 + 右侧图例）。
→ 上一条的 `fig_ident_filter_combined.R` 是 **patchwork 左右拼两块**（两套坐标系、两条 x 轴）；
本图是**原生 UpSet 本体版式**：**只有一套坐标系**，set-size 条、交集柱、点阵三者共享同一组交集位置。
**左侧 set-size 面板也按过滤结局堆叠** —— 这就是把「按工具堆叠柱」并进 UpSet 的做法。
后续七次追加修订：① 「保留的放里面呢」→ 两侧绿段改贴坐标轴（见硬规则 8）；
② 「左侧也标保留率」→ 左侧 7 根条外端各加一个保留率标签（见硬规则 9–11 + 专节）；
③ 「图例放到左上角呢」→ 图例从交集柱面板底部**搬到左上角空白格**（见硬规则 12）；
④ 「绘制 6 个软件的 venn 图」→ 交集柱面板空白区**内嵌 6 工具真 Venn**（见硬规则 13–15 + 专节）；
⑤ 「venn 图放大些 **放中间**」→ Venn **放大 ~31 %（2.18×2.26 → 2.85×2.96 in）且移到面板中部**；
（用户随后更正 ④ 里的「右上角」是口误，真意就是中间；缺口径 = **过滤后**数据 + **只算 6 个病毒工具**。）
⚠ **第 5 次修订同时暴露出 ④ 的一个严重缺陷**：`inset_element()` 作用在 **patchwork 对象**上时
框 npc 被**静默误解析**（详见硬规则 15）→ 上一版 Venn 实际只渲染出 **2.18 in**（应为 2.60+），
"上一版看着放得下"其实是被静默缩了一半。
⑥ 「这一部分空白的，可以补充每个软件的时间消耗的条形图嘛」+「**单样本平均耗时含 viroid_blastn**」
→ 左上格图例右/下方空白块（实测 **x400..1245 y280..870**，墨迹 0 像素）内嵌 **7 轨道耗时条形图**
（见硬规则 16–19 + 专节）。

| 项 | 内容 |
|----|------|
| 脚本 | `fig_ident_filter_stacked_all.R` |
| 产出（默认 `--set-size stack`） | `<S10_SCI_Discovery>/Fig_GojiVirome_IdentFilterStackedAll.{png,pdf}`（**11.4 × 7.4 in** / 3420 × 2220 px） |
| 产出（`--set-size fill`） | 同上加后缀 `_norm` |
| 数据源 | 与图④/⑥/⑦**同一份** `identfilter_data/ident_filter_upset_signature.tsv` |
| 耗时图数据源 | `identfilter_data/tool_timing.tsv`（由 `make_tool_timing_data.py` 从服务器 `02a_Identification/<sample>/<sample>_resource.tsv` + viroid mtime 差分生成） |
| 关键开关 | `set_sizes = upset_set_size(geom = geom_bar(aes(fill = status), position = "stack"\|"fill"))` + `encode_sets = FALSE` + `--no-venn`（关 Venn）+ `--no-timing`（关耗时图）+ `--no-legend-ov`（关图例 overlay） |
| 图内标注 | 左条外端逐工具保留率 + 右交集柱顶保留率（同一套语义：绿段 / 该柱整根） |
| 图例位置 | ⭐ **Venn 右侧一列**（`LEG_BOX = x2890..3350 / y430..1290`，**竖排 4 行一列**，由 **overlay** 画，非 patchwork 格子；见硬规则 19–20）。⚠ 第 8 次修订曾在 `x2560..3320/y435..640`（压到 Venn 顶部）→ 第 9 次修订改正 |
| **Y 轴题位置** | ⭐ **交接柱面板左侧、紧贴刻度文字**（overlay，中心 `x1220 / y975`，90° 旋转；**不用原生 `ylab`** —— 后者脱轴 386 px 且 `labs(y=NULL)` 关不掉，见硬规则 21）。第 10 次修订（2026-10-13 21:03） |
| Venn 位置 | **交集柱面板中部**（6 病毒工具真 Venn，`ggVennDiagram`，过滤后 Retained 数据；**预渲 PNG + `rasterGrob` 叠画**，见硬规则 13–15） |
| 耗时图位置 | ⭐ **左列中部大块**（`x51..1156 / y545..1494` = **3.6833 × 3.1633 in**，比例 1.164；7 轨道单样本平均耗时横条，对数刻度，**行序 = 右点阵同序**；**预渲 PNG + `rasterGrob`**，见硬规则 16–18、22、22b、24）。⚠ 沿革：第 7–10 次 `x400..1245 / y280..870` → 第 11 次下移 `x180..1250 / y1080..1560`（3.57×1.60 in）→ 第 12 次按红框放大到 `x50..1162 / y424..1381`（3.71×3.19 in）→ **第 13 次修正 Y 映射偏移（−110 px）整块下移 121 px** 到现框 |

```bash
Rscript scripts/barbarum/fig_ident_filter_stacked_all.R --set-size stack   # 绝对条长（忠实原生语义）
Rscript scripts/barbarum/fig_ident_filter_stacked_all.R --set-size fill    # 归一化 100%（推荐看图用）
# ⚠ 耗时图是**预渲**的：改了 tool_timing.tsv 需先重出 PNG，再跑主脚本
python scripts/barbarum/make_tool_timing_data.py        # → identfilter_data/tool_timing.tsv
Rscript scripts/barbarum/make_tool_timing_fig.R         # → identfilter_data/tool_timing_fig.png
```

### 两种尺度的取舍（本图的核心设计约束）

7 条鉴定轨道的命中数相差 **750 倍**（452,885 vs 604）。在约 2.1 in 宽的 set-size 面板里：

| 轨道 | 命中 | 绝对轴条长（@300 dpi） |
|---|---:|---:|
| viralverify | 452,885 | 全宽 |
| metabuli | 144,036 | ≈ 32 % |
| blast | 95,565 | ≈ 21 % |
| genomad | 11,526 | ≈ 15 px |
| virbot | 3,219 | ≈ 4 px |
| rdrpcatch | 2,228 | ≈ 3 px |
| viroid_blastn | 604 | **< 1 px（不可见）** |

→ **`stack`** 忠实于「set-size 面板 = 显示集合大小」的原生语义，代价是后 3 根的**段完全不可读**；
**`fill`** 每根等宽、构成完全可读，绝对数由**行标签**（`virbot ( 3,219)`）承载。
**本数据推荐 `fill`**；两种都已产出并双端验收。

### 本轮新增硬规则（24 条踩坑记录）

1. ⚠⚠ **`upset_set_size()` 的 `mapping=` 形参不可用**。传 `mapping = aes(fill = x)` 会在 `compute_aesthetics`
   阶段报 `object 'x' not found`（该 mapping 被当作**图层之外**的表达式求值）。
   **必须把 aes 写进 geom 里**：`geom_bar(aes(fill = x), …)`。
2. ⚠⚠ **必须 `encode_sets = FALSE`**。集合名含空格 / 括号 / 逗号（`"viralverify (452,885)"`）不是合法 R 名
   → 默认 `encode_sets = TRUE` 会把集合**编码成整数 id**，于是 set-size 面板自己的离散轴标签
   变成 `5 3 1 2 6 4 7` 这样的乱码数字（首版实测，按字母序编码）。
   关掉后名称原样保留；再配 `themes$overall_sizes` 里 `axis.text.y = element_blank()`，
   集合名由**矩阵面板统一显示一次**（与参考版式一致，不重复）。
3. ⚠⚠ **set-size 面板不共享 `fill` scale**。它是独立 ggplot → 必须**单独**补
   `scale_fill_manual(values = STATUS_COLORS, drop = FALSE)`，否则用 ggplot 默认色轮
   （红 / 橄榄 / 青 / 紫），与图例完全对不上（首版实测）。
4. ⚠ **`coord_flip()` 后主题元素名指「渲染后」的轴**，不是逻辑轴：
   `axis.text.y` / `axis.title.y` = **垂直**轴（集合名），`axis.text.x` / `axis.title.x` = **水平**轴（数值）。
   首版把 `axis.title.x` 置空 → 把刚写的 `ylab("Hits per identification track")` 一起干掉了。
5. ⚠ **set-size 面板的统计口径 = 展开表行数，不是记录口径 Σ`n_total`**。两者必然不等：
   展开表**不含**无 02b 产物的那部分 → 逐工具差 19/5/8/0/0/1/0，合计 **33**。
   ⚠ **33 ≠ 27**：27 是**去重后**的 nofilter 序列数，33 是**按工具重复计数**（一条序列可被多工具命中）。
   → **行标签必须用展开口径**，否则标签数字与条长对不上。脚本内有自检表打印两种口径。
6. ⚠⚠ **题注必须双端「量宽 + 量高」自检，且判据按最宽的一端**。服务器字形比本机**宽约 12 %**
   → 首版副标题第 2 行 268 字符，本机占画布 **96 %**（看着没超）、服务器**被右边缘裁掉**
   （断在 `the absolute hit cou`）。**只在本机看"没超"不算通过。**
   修法：① 每行 ≤ **170 字符**；② 脚本内 `convertWidth(grobWidth(textGrob(sub_txt, gp)))` 实测，
   超 `W_FIG − 0.30 in` 直接 `stop()`。
   ⚠ 同时**必须量高**：多行副标题若沿用系列里「`y` 当中心点」的写法，块高的一半会落到锚点下方
   → 行数一多就被**标题带下边缘裁掉**（实测 4 行版第 4 行只剩 **11 px** 高）。
   修法 = 标题与副标题改 **`vjust = 1`（顶对齐）**，并实测 `grobHeight`。
7. ⚠⚠ **`plot_layout(heights = c(a, 1))` 的 `a` 是相对权重，不是画布高占比**：
   标题带实际高度 = `H_FIG × a/(a+1)`。首版按 `H_FIG × a` 算自检 → 可用高度**偏大 ~18 %** → **假通过**
   （`a = 0.18` 时按 1.33 in 算、实际只有 1.13 in）。自检公式必须用 `a/(a+1)`。
8. ⚠⚠ **堆叠顺序：`Retained` 必须贴坐标轴（内侧）**（用户 2026-10-13 追加「保留的放里面」）。
   默认 `position_stack()` / `position_fill()` 把**第 1 个因子水平放在最外端**
   → 竖向柱绿段贴顶、横向条绿段贴左（远离轴），**绿色长度各自终点对齐、无法横向比长短**。
   修法 = 两侧都改用 `position_*_reverse(TRUE)`：
   - 竖向交集柱：`intersection_size(..., position = position_stack(reverse = TRUE))`
     （⚠ `intersection_size` **有 `position` 形参**，别自己拼 `geom_col`）
   - 横向 set-size 条：`geom_bar(..., position = position_fill(reverse = TRUE))`
     （该面板 x 轴是 `scale_y_reverse()` 翻转的，0 在**右**端 → 绿段贴右 = 贴轴）
   → 改后绿色段共享同一条基线，**保留率可直接比长短**；`Viroid track` 被推到最外端，
   不与病毒三态混叠。
9. ⚠⚠ **左侧逐工具保留率标签：取列必须用 `LABELS[[t]]`，不能写 `df[[t]]`**（用户 2026-10-13 追加
   「左侧也标保留率」）。脚本上文已经 `setnames(df, TRACKS, unname(LABELS[TRACKS]))` 把列名换成
   带计数的标签（`"viralverify (452,885)"`）→ 再用**裸轨道名** `df[["viralverify"]]` 取列会拿到
   `NULL`，`sum(NULL & …)` **静默得 0**（首版实测 7 个标签全打印 `0`、`0.00%`，无任何警告）。
   ⚠ 一般化：**同一脚本里对数据框做过整体改名后，所有下游取列都必须走同一张映射表。**
10. ⚠⚠⚠ **`scale_y_reverse()` 下 `expansion(mult = c(a, b))` 的方向是反的**（本图最大坑）。
    两个数按**变换后**的区间分配：`trans = "reverse"` 把数据 [0, 1] 映成 [-1, 0]
    → **`a` 落在变换后下界 = 数据上界 `y = max`**（set-size 条的**外端 / 左侧**），
    `b` 落在数据 0（**内侧 / 右侧**）。
    → 标签贴在 `y = 外端`，**余量必须写在第一个位置**。
    ⚠ 写成直觉上的 `c(0.004, 0.13)` 会把余量全留到**右侧空处**，左侧外端只剩 0.4 % →
    7 个标签**全部被面板左缘裁掉**，而且：
    - **不报任何错**；
    - `ggplot_build()` 里 `x` / `y` / `label` 三列**齐全、坐标也对**（实测 `x = 1..7, y = -1`）；
    - 图层在 `$layers` 里**位置正确**（`upset()` 对 set_sizes 是
      `layers = c(stripes, geom, user_layers)` → 用户层在**最后**，不存在被盖的问题）。
    → **只能靠量像素发现**。实测 A/B：`c(0.004, 0.13)` → 标签色 **0 px**；
      `c(0.13, 0.004)` → **97 px**；加到 `c(0.30, 0.004)` 仍 **97 px**（方向对了之后余量已够）。
11. ⚠ **余量是「面板宽的百分比」而标签宽度是「物理尺寸」** → 必须按**服务器端**字形算。
    服务器字形比本机宽 ~12–14 %。本机实测（fill 模式，3420 px 宽画布）：
    条长 737 px ↔ 数据 1 单位，面板宽 ≈ 857 px，`0.16 × 737 = 118 px` 余量区；
    最长标签 `"15.4%"` 实宽 **64 px**（服务器 **73 px**）→ 左端富余 47 px（服务器 34 px）。
    ⚠ 曾把标签宽**误估成 ~33 px**（实为 64 px）——**标签宽度必须实测，不能凭字号推算**。
    ⚠ 另：加了余量后 fill 模式的自动 `breaks` 会从 `0/25/50/75/100` 漂成 `0/30/60/90`
    → 已显式钉 `breaks = seq(0, 1, by = 0.25)`。
12. ⚠⚠ **图例搬到左上角 = 替换 patchwork 的 `spacer` 格**（用户 2026-10-13 追加「图例放到左上角呢」）。
    `upset()` 返回的是 **patchwork**，grid 布局 = **2 列 × 2 行**：
    ```
    widths  = c(0.3, 0.7)   →  [set-size 列 | 主列]
    heights = c(1, 0.5)     →  [交集柱行 | 矩阵行]
    +------------------+------------------+
    | spacer ← 放图例   | 交集柱 (Sequences)|
    +------------------+------------------+
    | set-size (Share) | 矩阵              |
    +------------------+------------------+
    ```
    左上格 (行1, 列1) 是一个 **`spacer`** —— 就是主图左上那片大白区，正好放图例：
    ```r
    idx <- which(vapply(p_body$patches$plots, function(x) inherits(x, "spacer"), logical(1)))
    stopifnot(length(idx) == 1)               # ⚠ 断言，版式一变立刻报错而不是静默错位
    p_body$patches$plots[[idx]] <- legend_plot
    ```
    三个必须记住的点：
    - ⚠⚠ **`widths` 是显式值** → 该列宽固定，放进去的对象**超宽会被裁掉且不报错** → 必须量宽自检。
    - ⚠⚠ **量宽不能写 `grobWidth(ggplotGrob(legend_plot))`** —— `ggplotGrob` 的宽高是 **null 单位**，
      转 in **恒得 `0.00`** → 自检**恒假通过**（首版实测打印 `0.00 × 0.00 in`）。
      正确量法：把图例改成**固定位置**（`legend.position = "right"`）后从 gtable 里抠出
      `guide-box` 再量（尺寸与数值定位版完全一致，只是摆法不同）：
      ```r
      g  <- ggplotGrob(legend_plot + theme(legend.position = "right"))
      gi <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
      leg_w <- convertWidth(grobWidth(g$grobs[[gi]]), "in", valueOnly = TRUE)
      ```
    - ⚠ **图例本体必须留一个 layer**：`ggplot()` 无图层会报 `no layers in plot`。
      用 `geom_col(alpha = 0)` 让它画不出来，再 `guide_legend(override.aes = list(alpha = 1))`
      把**图例色块**改回不透明 —— 不 override 的话图例也跟着透明，**整块图例看不见**。
    - ⚠ 搬走后必须在交集柱面板主题里显式 `legend.position = "none"`，否则同一份图例出两次。
    - 定位：`legend.position = c(0.02, 0.98)` + `legend.justification = c(0, 1)`（npc of 本格，左上角内缩 2 %）。
    实测：图例 **1.05 × 0.41 in**，左上格可用宽 **3.42 in** → 无裁切；
    双端图例包围盒起点**逐像素一致**（key 列 x 67..92，首行 y 474/475）。
13. ⚠⚠ **内嵌子图的框要按「目标中心反推」且必须先扫空白区**（用户 2026-10-13 追加）。
    - ⚠⚠ 若用 `inset_element()`，它在框内「**按纵横比等比缩放 + 居中**」，不是拉伸。
      首版给了个**宽扁**的框（x 1495..3362 / y 400..660 = **6.2 × 0.87 in**）→ 近方形的 Venn
      被缩到框高 0.87 in，**数字小到不可读**。→ 必须**先扫空白区、取近似正方形的块**再反推框。
    - ⚠⚠ **`inset_element(clip = FALSE)` 超框既不报错也不裁**（直接画到框外）→
      任何「量框尺寸」的自检都会**自欺**。**唯一可靠判据 = 截图差分**：
      `np.abs(base − cur).sum(axis=2) > 12` 抠出新增墨迹的包围盒，
      再在该区域检查 **base（关掉 Venn 的那版）content pixels == 0**（证明无碰撞）。
    - ⚠⚠ **量尺寸只能用「渲成 PNG 数像素」**：`grobWidth(ggplotGrob(p))` 恒 `0.00`（null 单位）；
      **即使 `png()` + `pushViewport()` + `grid.draw()` 之后回读，仍是 `0.00 × 0.00 in`**（实测）。
      ```r
      ggsave(tmpf, sub_plot, width = 3, height = 3, dpi = 300, bg = "white")
      vb  <- png::readPNG(tmpf)
      ink <- (vb[,,1] < .98) | (vb[,,2] < .98) | (vb[,,3] < .98)
      frac_w <- (diff(range(which(colSums(ink) > 0))) + 1) / (3 * 300)   # in
      ```
    - ⚠ **本图 `patches$plots` 只有 3 个**：`[1] spacer`、`[2]` **交集柱面板本体**、`[3]` set-size。
      交集柱所占据的格子**不是空格** → **不能像硬规则 12 那样换格子**，只能叠画在上面。
14. ⚠⚠⚠ **`ggVennDiagram` 的 `label_alpha` 是「文字透明度」，不是标签底色**（同上修订）。
    - `label_alpha = 0` → **63 个区域标签全部完全透明消失**，而 `ggplot_build()` 里
      `x` / `y` / `label` 三列**齐全**、图层顺序也对 → **又一次零报错的静默丢失**。必须 `label_alpha = 1`。
    - ⚠ **不显式给 `scale_fill_gradient()` 会把整片椭圆染成近黑**，标签淹在里面 →
      补 `scale_fill_gradient(low = "white", high = "#009E73", guide = "none")`。
    - ⚠ **`ggVennDiagram:::process_data()` 是 S4 泛型**，必须传 `Venn(x)` 对象；
      传裸 `list` → `unable to find an inherited method for function 'process_data' for signature 'venn = "list"'`。
      正解：`ggVennDiagram:::process_data(ggVennDiagram::Venn(vsets), shape_id = "601")$regionLabel`。
    - ⚠ **自动标注会给「全部 63 个区域」都标**（实测其中 **5 个 `count = 0`**）→ 一堆 `0` 白占空间；
      且六椭圆有大量**小到放不下数字**的区域，自动标注会**与椭圆弧线重叠**（实测 `"44687"` 被弧压住）。
      → 改 `label = "none"`，自取 `$regionLabel`，**只标 `count > 0` 的 58 个**，并用
      **「白描边粗体 + 正色」两层 `geom_text` 叠绘**（halo），数字压在椭圆线上也读得出来
      （正色 `#1A5C3A`，与左侧保留率 / 右侧柱顶 `%` 同一套绿）。
    - ⚠ **6 集合只有两种内置形状**（`get_shapes()`）：`601`（多边形/椭圆，经典 Venn 观感）与
      `601t`（三角/放射星）。用户要「椭圆 Venn」→ **`601`**。
15. ⚠⚠⚠ **`patchwork::inset_element()` 作用在 patchwork 对象上时，框的 npc 被静默误解析**
    —— 本轮最严重的一处「零报错缺陷」（用户 2026-10-13 第 6 次修订「放大些放中间」时暴露）。
    - **现象**：`layout + inset_element(child, l, b, r, t, align_to = "full")` 里，
      `left/bottom/right/top` **不以整画布为基准**。标记块实测：目标中心 `x=2426 y=722`
      → 实际落 **`x=2056 y=598`**；"框"报成 **4.31 in**（真值 2.95）。**不报错、不警告。**
      真图上表现为 Venn 只有 **1.40 × 1.37 in**（应 2.60 × 2.68）→
      **上一版「能放下」纯属被静默缩了一半**。
    - ⚠ **`align_to` 三值（`"full"`/`"plot"`/`"panel"`）全都一样错**（逐一实测），
      所以「换 align_to 就好了」是错判。
    - ⚠ **`wrap_elements(full = layout) + inset_element(...)` 也不好**：触发 patchwork
      **二次 reflow**（整图重排，不是刚体平移）→ 墨迹从 `x182..3215` 变 `x203..3194`，
      对齐搜索到 ±6 px 边界仍有 9.88 残差。**wrap `p_body` 同理**（真图左半 mean diff 12.937）。
    - ✓✓ **正解 = 预渲 PNG + `rasterGrob` 叠画**（绕开所有 inset/单位陷阱）：
      ```r
      ggsave(pre, child, width = w_in, height = h_in, dpi = 600, bg = "white")   # 2× 过采样
      rg <- grid::rasterGrob(png::readPNG(pre), interpolate = TRUE,
                             width = unit(w_in, "in"), height = unit(h_in, "in"))
      grid::grid.draw(patchwork:::patchworkGrob(layout))          # 基底：gtable 直画，零 reflow
      grid::pushViewport(grid::viewport(x = cx, y = cy, width = w, height = h, just = "center"))
      grid::grid.draw(rg); grid::popViewport()
      ```
      ⚠ 为什么不用 `grid.draw(ggplotGrob(child))` 直接画：child 会**锚在 vp 左下角**
      且 `grobWidth(ggplotGrob(p))` 恒 **0.00**（null 单位）→ `just=` 给什么都无效
      （实测红块永远落在 2058~2270 之间，到不了 2426）。
      ⚠ **`patchwork:::patchworkGrob(layout)` + `grid.draw`** 本身是**零 reflow** 的
      （基线/加图后墨迹逐字节一致），可放心当底座。
      ⚠ 代价 = PDF 里子图是位图 → 故预渲 **600 dpi**（框 3.19 in → 1912 px）。
    - **验收判据（双端）**：全图 `|diff|>12` 像素 = **0**（本机）、左半 `x<1490` mean diff = **0.0000**；
      落位包围盒两端一致，中心差 ≤ 12 px。

16. ⚠⚠⚠ **在「同一个 R 进程」里 `ggsave()` 一张 `theme_minimal` 图，会让主图产生 ~18 K px 的
    子像素重栅格噪声** —— 第 7 次修订（耗时图）踩到的一个**零报错、零警告、肉眼几乎看不出**的坑。
    - **现象**：主脚本内新增 `ggsave(tim_pre, timing_plot, ...)` 后再画 patchwork，
      整图 `|diff|>12` 从 **0 跳到 158,679**；改动包围盒 `x38..3344 y546..2119`（几乎全图）。
      **但 `Δink` 只有 −273 px（−0.07 %）**，且遮挡区逐行墨迹数几乎不变
      → 是**字形落在不同子像素相位**造成的**重栅格**，**不是版式位移**（截图肉眼完全一样）。
    - **逐项排除（都失败）**：✗ 换 device（`ggsave` → 显式 `png()`+`print()`）✗ 加 `dev.off()`
      ✗ 换渲染尺寸（2.82 in → 3.19 in）✗ 把渲染调用移到 `upset()` 之前
      ✗ 用 `system2` 起子进程渲（gtable 已构好，重画反而失真）。
    - ✓ **定位手法**：把 `timing_plot` 换成 **`theme_void()` 玩具图** → **差异归 0**；
      只 `fread()` 数据不建 ggplot → **差异归 0**；`want_timing=FALSE` → **差异归 0**。
      → 锁定：**`theme_minimal` 图（含 `axis.text` / `panel.grid`，需量字形）的 `ggplot_build()`
      改写了「进程级字形度量缓存」**，而 ComplexUpset 的标签/格子尺寸依赖同一缓存 →
      后续 `grid.draw` 被改写。
    - ✓✓ **正解 = 该图「独立成脚本 + 独立进程」预渲 PNG，主脚本只 `readPNG` 贴**：
      ```bash
      Rscript scripts/barbarum/make_tool_timing_fig.R     # → identfilter_data/tool_timing_fig.png
      ```
      主脚本内**绝不构建**该 ggplot（连 ggplot 对象都不建），只 `rasterGrob(readPNG(...))`。
      → 验证结果：**排除叠加框后差异像素 = 0，Δink = +0**（与纯基线逐字节一致）。
    - ⚠ **推论（通用）**：凡「给 patchwork 主图叠一张自带坐标轴/字形的 ggplot」，
      **一律走「独立进程预渲 PNG + `rasterGrob`」**；不要试图在同进程 `ggsave` 后接着画主图。

17. ⚠⚠ **`scale_x_log10()` + `geom_col()` 会把整批条形静默丢光**（只发 warning）。
    - `geom_col()` 隐含 `xmin = 0`，而 `log10(0) = −Inf` → 全部落在轴范围外 →
      控制台一句 `Removed 7 rows containing missing values or values outside the scale range`，
      **图照样出，但一根条都没有**。
    - ✓ 修法：用 **`geom_rect()` 显式给 `xmin`**（= 轴下限，如 3），或用 `geom_bar(stat="identity")`
      配 `coord_cartesian` 手动换算。本项目采用前者（见 `make_tool_timing_fig.R`）。
    - ⚠ 刻度也要显式给：`breaks = c(3,10,30,100,300)` + `expand = c(0,0)`，
      否则对数轴会把「下限到第一档」留出巨大空白，条形看起来像被截断。

18. ⚠⚠⚠ **改「y 轴行序」时的两个静默错位坑**（第 8 次修订，图⑧ 耗时图踩到）。
    - **坑 A：`as.numeric(factor)` 定位 ≠ 轴标签顺序。** 把行序改成
      `factor(levels = rev(ORDER))` 后若仍用 `y = as.numeric(tool)` 画条，
      而 `scale_y_continuous(breaks = 1:7, labels = ORDER)` 按**原序**贴标签 →
      两套方向相反 → **条长与标签错位整行**（实测 viralverify 显示成 7.5 s，真值 30.0 s）。
      **零报错、零警告**，必须逐行对着 tsv 核才看得出。
      → 正解：显式 `ypos` 连续坐标，**条与标签共用同一个 `ypos`**。
    - **坑 B：`ypos` 取向反相。** `scale_y_continuous` 的 y **向上增大** → `ypos = 7` 在**最上**。
      若直接 `ypos := match(tool, ORDER)`（ORDER[1]→1），ORDER[1] 会落**最下** → 上下颠倒
      （实测：日志打印"上→下 = ORDER[1] > …"，图上却是 ORDER[末] 在最上）。
      → 正解：`ypos := N + 1 - match(tool, ORDER)`（反相），并加自检
      `stopifnot(最上那条 == ORDER[1])`。
    - ⚠ 附：**`identical(sort(int), as.numeric(1:7))` 恒 FALSE**（整数 vs 双精度）→
      守卫用 `all(a == b)`。此守卫**真的拦住了坑 A 的错位图**，别删。
    - **纪律**：改 y 轴类版式后**必须渲染出来逐行核对「标签 ↔ 条长」**，
      不能只看「行都在、顺序对」就收工。

19. ⚠⚠⚠ **多个 overlay 的叠加顺序 = 语义的一部分（白底 raster 会覆盖先画的）**。
    - 第 8 次修订把图例挪到右上角，其框 `x2560..3320` 与 Venn 框 `x1947..2905`
      **横向重叠 x2560..2905（245 px）** —— 两框在右上区**必然交叠**（该区是 Venn 右半
      与面板右缘之间的唯一留白带），**无法靠挪框避开**。
    - 而 **Venn 的预渲 PNG 是白底**（`bg = "white"`）→ 若 `overlays` 里图例排在 Venn **之前**，
      **Venn 的白色矩形会把图例左半整块盖掉**，图面只剩右半两条 —— **零报错**，
      看着像"图例被裁"，实际是**被覆盖**。
    - → 正解 = **Venn 在下、图例在上**（`overlays` 里图例 append 在 Venn 之后）。
      通用纪律：**给白底 raster 排顺序时，必须先想清"谁盖谁"**；
      若哪天改成透明底 PNG，顺序才不再敏感（但现在**不可依赖**）。
      诊断手法：`--no-venn` 单独跑一次 → 若图例完好，即为覆盖而非裁剪。
    - ⭐ **第 9 次修订补充**：图例改到 `x2890..3350`（Venn 右缘 2858 之右）后
      **横向已不重叠**，但**顺序仍保持"Venn 在下、图例在上"** —— 因为两框纵向仍有
      重叠（y430..1290 vs Venn y406..1294），且这是**零成本保险**：Venn 将来若再放大，
      顺序对了才不会静默盖掉图例。⚠ **挪框能解决的只是"当前这一版"，
      顺序正确性才是长效保险，两者都要做。**

20. ⚠⚠⚠ **"空白带"必须用「本次改动后的干净基线」实测，不能沿用上一版的测量**。
    - 第 8 次修订把图例放 `x2560..3320 / y435..640`，依据是「该带墨迹 = 0」——但那次
      测量是在**上一版**（图例还在左上）基线上做的；等图例真搬过去，**它自己就成了
      该带的内容**，而 Venn 顶部墨迹（y 从 406 起、x 到 2856）**一直在那** →
      视觉上就是「图例挡到 Venn」。
    - → 纪律：**每次挪 overlay，都要重出一次"关掉该 overlay"的干净基线**再扫描空档
      （`--no-legend-ov --no-timing` / `--no-venn`），**不能复用历史坐标**。
    - 第 9 次修订的实测（干净基线 3420×2220）：Venn 墨迹右缘 **x=2858**；面板右缘 **x≈3362**；
      **`x2890..3350` 从 y0 到 y1470 全空**（唯一墨迹从 y1478 起，是下方 x 轴条）→
      取 **`x2890..3350 / y430..1290`**（460×860 px = **1.53 × 2.87 in**）作右侧列框。
    - ⚠ 扫描口径：`ink = (R<235)|(G<235)|(B<235)`；逐 40 px 列带 + 逐 100 px 行带各出一遍，
      **两向都为空**才算真空。

21. ⚠⚠⚠ **ComplexUpset 的 `base_annotations` **列表名**会被无条件写成该面板的 y 轴题，
    且 `labs(y = NULL)` 关不掉它**（第 10 次修订，图⑧ 踩到）。
    - **现象**：用户报「upset 的 Y 轴 Sequences 不显示」。核查发现题注**其实渲染了，
      但位置错**：`upset(base_annotations = list("Sequences" = intersection_size(...)))` 的
      intersections 面板在本 patchwork 里带**巨大的左侧 margin**（set-size 列 + spacer 列），
      `axis.title.y` 按「整个面板含左 margin」定位 → 题注飘到 **x≈952..984**，
      而 y 轴刻度文字在 **x≈1257..1499** → **离轴 386 px（≈1.3 in）**，
      孤零零浮在耗时图下方空白里 → 用户视觉判为「不显示」。
    - **根因**（读 `ComplexUpset::upset()` 源码第 168–171 行）：
      ```r
      if (is.null(annotation_plot$labels$y) ||
          (!is.null(annotation_plot$default_y) && annotation_plot$default_y == annotation_plot$labels$y)) {
        annotation_plot = annotation_plot + ylab(name)      # ← name = 列表名 "Sequences"
      }
      ```
      **`labs(y = NULL)` → `labels$y` 为 NULL → 条件成立 → 它又把 `ylab("Sequences")` 加回来。**
      ⚠ **`axis.title.y = element_blank()` 同样关不掉** —— theme 只控制「怎么画」，
      控制不了「标签内容是否被写入」（`labels$y` 是不是 NULL 才是判据）。
      → 两个都试过后实测**出图两个题注并存**（原生在 x≈952 + overlay 在新位）。
    - ✓✓ **正解 = 给一个「非空、且 ≠ 默认 y 映射」的哨兵串**：
      `labs(y = " ")`（单个空格）→ `labels$y` 非 NULL 且 ≠ `default_y`
      → ComplexUpset **跳过** `ylab(name)`；再在 `themes` 里用
      `axis.title.y = element_blank()` 把空格题注**隐藏**。**两道缺一不可**。
      ```r
      # base_annotations 内
      labs(y = " ") + theme_upset(8) + theme(axis.title.y = element_blank(), ...)
      ```
    - **题注位实测**（3420×2220 干净扫）：y 轴刻度文字左缘 **x≈1257**；
      紧贴其左的 **`x1190..1250` 竖带在 y249..1582 全空** → 题注中心取 **x≈1220**；
      纵向中心 = 交集柱面板主区中心（柱顶 ≈y480、基线 ≈y1470）→ **y≈975**；
      框 44 × 300 px（0.15 × 1.00 in），在 overlay 里 90° 旋转绘制 → 零 reflow 同 Venn/耗时/图例。
    - **验收（双端）**：老题注位 `x900..1050 / y860..1090` 墨迹 **2042 → 162 px**
      （残余 162 = 耗时图 x 轴题 `Mean runtime per sample…` 的降部，**非**题注）；
      新题注位 `x1180..1260 / y700..1300` 墨迹 **1918 px（本机）/ 1841 px（服务器）**；
      全图 `|diff|>12` 共 **3768 px，且全部落在两个题注位之一**（1888 新位 + 1880 老位）
      → **除题注外零回归**。

---

**22. ⚠⚠ 耗时条形图「整块下移」——框比例一改，`make_tool_timing_fig.R` 的 `--w/--h` 必须同步重渲，否则静默拉伸**

**触发（用户 2026-10-13 第 11 次修订）**：用户给了一张带**红框 + 右下箭头**的截图，
框内是**当前耗时图所在位置**，箭头指向其**下方的大片空白**（图题 `这个位置不行嘛`）。
→ 诉求不是「改内容」，而是**把耗时图整块下移到更空的地方**。

- ⭐⭐ **先量再动**：干净基线（3420×2220）像素扫描，
  下方空白带 `x150..1250 / y1080..1560` **墨迹 = 0**，
  且**不与任何元素相交**（`Sequences` 题注 `y700..1080` 在其上、Venn `x2260..` 在其右）
  → 新框取 **`x180..1250 / y1080..1560`**，中心下移约 745 px。
- ⚠⚠⚠ **核心坑：改框 = 改宽高比，预渲 PNG 必须按新比例重渲**。
  旧框 `x400..1245 / y280..870` → **2.8167 × 1.9667 in（比 1.43）**；
  新框 `x180..1250 / y1080..1560` → **3.5667 × 1.6000 in（比 2.23）**。
  框比例 1.43 → 2.23（**+56 % 更扁**）而 PNG 不重渲 → `rasterGrob(..., width=unit(w,"in"), height=unit(h,"in"))`
  会**把原图按新框硬拉伸**（文字横向压扁 56 %），**且零报错、零警告**。
  → 必须改 `make_tool_timing_fig.R` 的默认参数并**重跑预渲**：
  ```r
  W_IN <- as.numeric(get_arg("--w", "3.5667"))
  H_IN <- as.numeric(get_arg("--h", "1.6000"))
  ```
  **判据 = 出图前先 `identify` 预渲 PNG 的像素宽高**，应 ≈ `3.5667×600 = 2140` /
  `1.6000×600 = 960`；对不上就是没重渲。
- ⚠ **图例是唯一改框后不需要重渲的 overlay** —— 它是 `gtable` 直绘（矢量），
  只有 raster 资产（耗时图 PNG、题注 PNG）才吃宽高比。**同类改动先问「这个 overlay 是 raster 还是 vector」。**
- **验收（双端）**：
  旧槽 `x400..1245 / y280..870` 墨迹 **0**（**完全腾空**）；
  新槽 `x180..1250 / y1080..1560` 墨迹 **107,005 px**；
  `Sequences` 题注仍在 `y700..1080`、**位于新框之上**（不打架）；
  `_check_discovery_overlap.py` → 全场 **0 对重叠**（20 张图）；
  自检行双端一致：`耗时图框：x 0.053..0.365 y 0.297..0.514 npc｜3.57 × 1.60 in（预渲 PNG）`。
- ⚠ **同步清单**（缺一即不一致）：① 本脚本 `TIM_L/R/T/B` 四个常量；
  ② `make_tool_timing_fig.R` 的 `--w/--h` 默认值；③ **重跑预渲**（用第 ① 步脚本或单独跑）；
  ④ README 汇总表「耗时图位置」行 + 本硬规则；⑤ MEMORY.md `⑧` 条目。

**22b. ⚠⚠ 第 12 次修订：耗时图「整块放大到用户红框」（2026-10-13 21:56）**

**触发**：用户截图（1161×612）在**左上角**画了一个红框，框内只有耗时图（占红框下半部），
红框**上半部大片空白**；题注 `再次调整大小，符合红色区域的大小，使其与预期效果一致，
并明确说明希望的具体尺寸或调整方向，避免重复修改。`
→ 诉求 = **把耗时图放大到红框那么大**（不是移动，是**改大小**）。

- ⭐⭐ **红框 → 画布 px 的映射（截图可能被裁，必须拟合）**：
  ⚠ 该截图**裁掉了顶部 341 px**（主标题 + 副标题均未入镜），故**不能**直接用 `图宽比` 换算。
  正解 = **行/列锚点线性拟合**：用截图的「x 轴刻度带」与「底部题注带」两段的行坐标，
  对齐画布的对应两段（`y2094..2118` 与 `y2130..2165`），再用墨迹列首尾对齐定横向偏移。
  拟合结果（已用**裁剪对拍**验证：截红框内 vs 画布同区，7 个工具名与数值逐一吻合）：
  ```
  act_x = (shot_x - 14.25) / 0.3355        # 等比 0.3355（= 1/2.981）
  act_y = (shot_y + 114.4) / 0.3355
  ```
  红框边框外缘 `shot x31..404 / y28..349` → 画布 **`x50..1162 / y424..1381`**。
- **新框 = `3.7067 × 3.1900 in`（比例 1.162）**；旧框 `3.5667 × 1.6000 in`（比例 2.23）。
  → 框比 **2.23 → 1.162（−48 %，变方）** ⇒ **必须重渲**：
  ```r
  W_IN <- as.numeric(get_arg("--w", "3.7067"))
  H_IN <- as.numeric(get_arg("--h", "3.1900"))
  ```
  判据 = 预渲 PNG 像素应 ≈ `3.7067×600 = 2224` / `3.1900×600 = 1914`（实测 2224×1914 ✓）。
- ⭐⭐ **净空必须用 `--no-timing` 底图判，不能用成图**：
  成图里旧框区域**本身就被耗时 overlay 占着**，拿它当"基线"会把 overlay 内容误判成底图元素。
  实测底图（`--outdir <临时目录> --no-timing`）：
  | 区域 | 墨迹 |
  |---|---|
  | 新框全区 `x50..1162 / y424..1381` | **0** |
  | 旧框缩小后**新暴露**的右条 `x1162..1250 / y1080..1381` | **0** |
  | 旧框缩小后**新暴露**的底条 `x50..1162 / y1381..1560` | **0** |
  → **放大不压元素、缩小不露元素，双向安全**。
  ⚠ 上界 `y424` 与 `Sequences` 题注（`x1198..1242`）横向不重叠（`1162 < 1198`，间隙 41 px）✓
- **验收（双端，两版均同）**：
  框内墨迹 **225,006 px**（已填满）；上邻 `y380..423` = **0**、下邻 `y1382..1450` = **0**、
  左邻 `x0..49` = **0**、右邻 `x1163..1300` = 2,821（= `Sequences` 题注，**在框外**）；
  四边各 4 px 全 **0**（无 `TEXT-OFF-CANVAS`）。
  **零回归判据**：新旧差分 `|Δ|>30` 共 299,612 px，bbox `x55..1078 / y428..1559`
  = **恰好等于耗时图新旧框的并集**（PNG 右侧白边不计）→ **除耗时图外零变化** ✓
  自检行：`耗时图框：x 0.015..0.340 y 0.378..0.809 npc｜3.71 × 3.19 in（预渲 PNG）`。
- ⚠ **本轮同时踩到并已修的静默 bug**：`make_cls_tool_timing_fig.R`（分类图那个）
  的 `--out` 默认名**不随 `--rank` 变化** → 两个阶元互相覆盖。详见硬规则 23。
- ⚠ **同步清单**（同硬规则 22，五项）：① `fig_ident_filter_stacked_all.R` 的 `TIM_*`；
  ② `make_tool_timing_fig.R` 的 `--w/--h` 默认值；③ **重跑预渲**；④ 本 README 汇总表 + 本节；
  ⑤ MEMORY.md `⑧` 条目。

**23. ⚠ 同族预渲脚本：`--out` 默认名必须随「区分参数」变化，否则静默互相覆盖**

`make_cls_tool_timing_fig.R`（分类赞成票图的耗时预渲）与 `make_tool_timing_fig.R` 同族。
它的 `--rank`（`Species` / `Family`）会改变**行序**（末两位互换），若两者都写默认名
`cls_tool_timing_fig.png` → `--rank Family` 那次**静默覆盖** `Species` 的产物（零报错）。
→ 正解：默认名按 `rank` 推导（`get_arg("--out", "")` 为空时用 `sprintf(...rank...)`），
并在脚本头写明该约束。**通用判据**：任何「同一脚本多模式、共享输出目录」的预渲，
先问「默认输出名是否唯一标识该模式」。

**24. ⚠⚠⚠ 第 13 次修订：红框→画布映射的 Y 偏移错了 110 px（耗时图偏高 121 px）——锚点必须用矩阵点阵**

**触发**：用户 2026-10-13 21:55 再次发同一张截图 + 红框（题注
`再次调整红色区域的大小，使其与预期效果一致，并明确说明希望的具体尺寸或调整方向，避免重复修改。`）。
第 12 次修订已按该红框放大过，但**位置仍不对**。

- ⭐⭐⭐ **根因 = 第 12 次修订的映射 Y 偏移错 ≈110 px**（X 也偏 +8..+20 px）。
  它用的锚点是「x 轴刻度带 / 底部题注带」两段行坐标 —— 这两段**太细、太靠下**，
  截图的缩放 + 裁剪误差被放大。实测残差：

  | 锚点 | shot_y | 画布实际 | 第 12 次预测 | 本版预测 |
  |---|---:|---:|---:|---:|
  | 矩阵第 1 行 | 383 | 1594 | 1483（**−111**） | **1594（0）** |
  | 矩阵第 7 行 | 531.5 | 2033 | 1925（−108） | **2033（0）** |
  | 左面板首条 | 388.5 | 1585 | 1499（−86） | 1610（+25） |
  | Identification track | 570 | 2146 | 2040（−106） | **2147（+1）** |

- ⭐⭐ **正解锚点 = UpSet 点阵的 7 行棕点（+ 16 个点列）**。理由：它们是**等距、贯穿全宽、颜色独特**
  的硬基准，且在两张图里都能量准。**识别指纹 = 逐行棕点数递减**（画布 12/11/13/8/4/1/1、
  截图 8/7/8/6/3/1/1）—— 这是 UpSet 点阵「上排集合参与交集多、下排少」的必然形态，
  两侧同现即可确认这 7 行是同一组。
  ```
  act_x = -41.1 + 2.9618 * shot_x        # 16 个点列标定，残差 0 px
  act_y =  462.5 + 2.9545 * shot_y       #  7 个点行标定，残差 0 px
  ```
  红框外缘 `shot x31..404 / y28..349` → 画布 **`x51..1156 / y545..1494`**。
  ⚠ 第 12 次修订「截图裁掉顶部 341 px」的判断**正确**，保留。
- **新框 = `3.6833 × 3.1633 in`（比例 1.164）**；上版 `3.7067 × 3.1900 in`（比例 1.162）。
  → 框比几乎不变（1.162 → 1.164），但仍**按纪律重渲**：
  ```r
  W_IN <- as.numeric(get_arg("--w", "3.6833"))
  H_IN <- as.numeric(get_arg("--h", "3.1633"))
  ```
  判据 = 预渲 PNG ≈ `3.6833×600 = 2210` / `3.1633×600 = 1898`（实测 2209×1897 ✓）。
- ⭐ **净空实测**（画布逐像素扫）：
  | 区域 | 墨迹 |
  |---|---|
  | 新占带 `x51..1156 / y1382..1494`（下移后新增） | **0** → 可安全移入 |
  | 让出带 `x51..1156 / y424..545` | 30,661（现图顶部，移走后腾空） |
  | 与下面板之间 `y1495..1555` | **0**（左面板上缘 `y≈1555`） |
  | 右缘外窄带 `x1157..1190 / y545..1494` | **0**（`Sequences` 题注左缘 1198） |
- **验收（双端）**：耗时条行段 `y573..1396`；框内墨迹 bbox `x56..976 / y550..1490`
  —— **完整落在目标框 `x51..1156 / y545..1494` 内**（四边各余 4–5 px）；
  自检行双端一致：`耗时图框：x 0.015..0.338 y 0.327..0.755 npc｜3.68 × 3.16 in（预渲 PNG）`。
- ⭐⭐ **方法论（本轮最值钱的一条）**：
  **映射基准要选「等距、贯穿、颜色独特」的硬结构（点阵/网格），不要选「细带/文字行」。**
  且**必须报告残差**——第 12 次修订只给了映射式、没给残差，110 px 的系统偏差就漏过去了。
  ⭐ 另一条：用户「再画一次同样的框」= **上一版没对齐**，先怀疑自己的映射，别怀疑用户的意图。
- ⚠ **同步清单**（同 22/22b，五项）：① `TIM_L/R/T/B`；② 预渲 `--w/--h`；③ **重跑预渲**；
  ④ 本 README 汇总表 + 本节；⑤ MEMORY.md `⑧` 条目。

### 逐工具保留率标签（用户 2026-10-13 追加）

左侧 7 根条的**外端**各标一个保留率，语义与右侧交集柱顶的 `%` **完全一致**
（= 绿段 / 该工具整根条）。实测值：

| 轨道 | retained | hits | 保留率 |
|---|---:|---:|---:|
| viralverify | 69,817 | 452,885 | **15.42 %** |
| metabuli | 30,425 | 144,036 | **21.12 %** |
| blast | 26,811 | 95,565 | **28.06 %** |
| genomad | 6,918 | 11,526 | **60.02 %** |
| virbot | 2,682 | 3,219 | **83.32 %** |
| rdrpcatch | 1,514 | 2,228 | **67.95 %** |
| viroid_blastn | 11 | 604 | **1.82 %** |

⚠ **`viroid_blastn` 的 1.82 % 是「绿段 / 整根条（604）」，不是「占进入过滤链的 88 条」**
（后者 = 12.5 %）。副标题已声明 `88 of its 604 hits enter the filter chain`，两处口径不可混。

双端像素验收（归一化版）：

| 项 | 本机 | 服务器 |
|---|---|---|
| 标签数 | 7 / 7 | 7 / 7 |
| 标签宽 | 62–64 px | 72–73 px |
| 标签右端与条外端间隙 | 8 px | 10 px |
| 标签左端距面板左缘 | 47 px | 34 px |
| 是否被裁 | 否 | 否 |

### 中部 6 工具 Venn（用户 2026-10-13 追加；第 6 次修订改「放大 + 居中」）

用户原话「绘制 6 个软件的饼图」→ 随即改为「**绘制 6 个软件的 venn 图**」
（第 5 次修订时口述为「右上角」，**第 6 次已更正为口误，真意 = 面板中部**）。
原数据是 7 轨道，直接画 6 椭圆 Venn 会因**大量区域为空**而不可读 → 经二次确认口径：

| 确认项 | 用户选定 |
|---|---|
| 数据范围 | **过滤后（Retained）** 的签名，**不是**全部命中 |
| 集合 | **只算 6 个病毒工具**，**排除 `viroid_blastn`**（独立类病毒轨道） |

**口径换算**（脚本内 `make_tool_venn()`，并带一致性自检）：

- 取 `n_retained > 0` 且 **signature 第 7 位 == `"0"`**（= `viroid_blastn` 未命中）的签名；
- 用 `mapply` 按 `n_retained` **展开成逐序列 id**，再按位取 6 个集合；
- **自检**：Venn 并集必须 == 「6 工具任一命中」的 Retained 总数（`sum(n) − n["000000"]`），
  不符直接 `stop()`（口径一致性闸门）。

| 项 | 值 |
|---|---:|
| 并集（≥1 工具命中） | **94,807** |
| viralverify | 69,817 |
| metabuli | 30,414 |
| blast | 26,811 |
| genomad | 6,918 |
| virbot | 2,682 |
| rdrpcatch | 1,514 |

按 **k（命中工具数）** 分布：`k=1..6` = **65,760 / 18,343 / 8,578 / 1,179 / 422 / 525**
→ 以单工具命中为主（69 %），但六个工具**确实两两可辨**，Venn 各区域都有内容（63 区域中 58 个非零）。

**空白区测绘（决定 Venn 落位，全部在关掉 Venn 的 base 图上量）**：

| 项 | 范围（3420 × 2220 px） |
|---|---|
| 交集柱面板 | x **1352..3362** |
| 主横网格线 | y = 670 / 948 / 1225 / 1503 |
| 柱基线 / 最高柱顶 | y = 2051 / y = 546 |
| 16 根柱 pitch | ≈ 124 px（首根中心 x ≈ 1425，末根 ≈ 3288） |
| **最左（最高）柱** | x **1370..1481** → **x ≥ 1490 起全空**（左界） |
| **柱顶 `%` 标签带** | 从 **y ≈ 1330** 起（"13 %"/"21 %"… 实测 y1330..1490 有内容）→ 下界 |
| 面板顶 / 右侧刻度标签区 | y ≈ 430（再上是 y 轴刻度「3e+05」的**空白**） |
| 内部唯一障碍 | y = 670 一条浅灰网格线（235,235,235）→ Venn 内部为白，**可忽略** |
| **可用块** | **x 1490..3362、y 430..1320（1872 × 890 px = 6.24 × 2.97 in）→ 高度是瓶颈** |

→ 第 6 次修订最终框：**中心 x = 2426、y = 850，边长 3.19 in（957 px）**（框 x 1948..2905 / y 372..1328）。
先实测 Venn 自然墨迹占比 **90 % × 93 %** → 框内墨迹 **2.86 × 2.96 in**，
落 **x 1996..2856 / y 406..1294**（实测 bbox 2004..2859 / 409..1293）：
**左缘距最左柱 514 px、下缘距标签带 37 px、右侧不越面板缘** → 无碰撞 ✅。
（对比：改前的「右上角」版只有 **2.18 × 2.26 in**，本次 **+31 % / +30 %**。）

**双端验收**（2026-10-13，第 6 次修订，两种模式均通过）：

| 指标 | norm（fill） | stack |
|---|---|---|
| 画布 | 3420 × 2220 两版同 | 同 |
| Venn 框 | **3.19 × 3.19 in**（中心 2426, 850 px）| 同 |
| Venn 自然墨迹占比 | 90 % × 93 %（本机）/ 91 % × 93 %（服务器）| 同 |
| Venn 落位 bbox（本机） | x 2004..2859 / y 409..1293 → **2.85 × 2.95 in** | 同 |
| Venn 落位 bbox（服务器） | x 2006..2844 / y 447..1229 → **2.80 × 2.61 in** | 同 |
| Venn 中心（本机 / 服务器） | 2432,851 / 2425,838（差 ≤ 12 px）| 同 |
| 左缘距最左柱 | **514 px** | 同 |
| 下缘距柱顶 `%` 标签带 | **37 px** | 同 |
| 全图 `\|diff\|>12`（对比 `--no-venn`）| **0 px**（零 reflow）| 同 |
| 左半 `x<1490` mean diff | **0.0000** | 同 |
| 标注区域数 | 58 / 63（`count = 0` 的 5 个不标）| 同 |
| 脚本 md5（双端一致） | `63348b9f2082b67a93434e67d89709b1` | 同 |

⚠ 服务器字形比本机宽 → 服务器版 Venn 墨迹高 **2.61 in**（本机 2.95），
即同一框内**墨迹更紧凑**（占比 93 % 不变 → 是本机字形更小、横向撑得更开所致）；
两端**中心与左缘一致**，版式无差异。本版一律**预渲 600 dpi**，故 PDF 中 Venn 为位图。
⚠ **改前「右上角」版的 2.18 in 是被 `inset_element` 静默缩掉的结果**（见硬规则 15），
并非真实可用尺寸 —— 旧验收表（缩放比 0.797/0.809）随之作废。

### 数值验收（像素实测）

归一化版逐条测量 set-size 条的段占比，与记录口径**逐条吻合**：

| 轨道 | 绿（Retained）像素实测 | 记录保留率 | 橙（Viroid）像素实测 |
|---|---:|---:|---:|
| viralverify | 15.4 % | 15.4 % | 0 |
| metabuli | 21.0 % | 21.1 % | 0 |
| blast | 28.0 % | 28.1 % | 0 |
| genomad | 59.9 % | 60.0 % | 0 |
| virbot | 83.3 % | 83.3 % | 0 |
| rdrpcatch | 67.9 % | 68.0 % | 0 |
| viroid_blastn | 12.1 %（= 11/88） | — | **85.4 %**（= 516/604） |

→ 左右两侧结论**互相印证**：工具越「专一」（命中越少），保留率越高（15.4 % → 83.3 %）。

**双端复现**（2026-10-13，逐工具标签 + 图例移到左上角后重测）：本机与 246 各跑两种模式均成功；
脚本 md5 双端一致（`0b74e59e9f084705b7e4d338d66aa1e8`）；两版产物均 3420 × 2220 @300 dpi。

| 指标 | norm（fill） | stack |
|---|---|---|
| 4× 降采样均差 | 7.19 | 6.78 |
| `> 40/255` 像素占比 | 5.98 % | 5.66 % |
| 墨迹比（服务器 / 本机） | 0.988 | 0.988 |
| 图例包围盒起点 | x 67 / y 474（本机）= x 67 / y 475（服务器） | 同 |
| 图例宽 本机 / 服务器 | 309 px / 340 px（上限 1026 px） | 同 |

→ 均在「字形 hinting 基线 ~7.2 %」以内。⚠ **不要用文件大小判「是否同一张图」**：
本机 `png()` 出**调色板图**（全图仅 **349 色 / 161 KB**），服务器真彩（**1496 色 / 355 KB**）
→ 同一张图差 **2.2×**。判据用「尺寸 + 4× 降采样均差 + 墨迹比」。

**双端复现（第 6 次修订：Venn 放大 + 居中，2026-10-13）** —— 上表为第 3–4 次修订基线；
本轮脚本 md5 双端一致（**`63348b9f2082b67a93434e67d89709b1`**），两版产物均 3420 × 2220 @300 dpi，
逐项数字见上文「中部 6 工具 Venn」的双端验收表（Venn bbox 本机 2004..2859 / 409..1293、
服务器 2006..2844 / 447..1229，中心差 ≤ 12 px，全图 `|diff|>12` = 0）。
⚠ 第 5 次修订的旧数字（缩放比 0.797/0.809、墨迹比 0.9998）**作废** —— 那是 `inset_element`
静默误解析后的产物（见硬规则 15）。

⚠️ 部署纪律：服务器脚本先备份 `.bak_20261013_pre_midvenn`（35,539 B），再 `scp`；
`md5sum` 双端校验一致。**README 同步**同样先备份（`.bak_20261013_pre_venn`），
本轮改动全部为**追加 / 就地改写**，不改动并发会话的 `## fig_novelty_sankey.py` 段落。

⚠ 脚本会打印 4 条 `Ignoring unknown aesthetics: fill`（图④ 同样打印 2 条）——
来自 ComplexUpset 内部对不支持 `fill` 的几何对象套用同一套 aes，**无害**。

### 左上格 7 轨道耗时条形图（用户 2026-10-13 第 7 次修订追加）

用户原话：「这一部分空白的，可以补充每个软件的时间消耗的条形图嘛」
+ 「你在 `~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out` 看看**单样本平均耗时含 viroid_blastn**」。

**目标框（实测）**：`--no-venn` 干净基线上逐像素扫，图例 4 行占 **x67..376 / y474..592**；
其右侧 **x400..1245 / y280..870 墨迹 = 0 像素** → **845 × 590 px = 2.82 × 1.97 in**。
右界 1245 因交集面板 y 轴刻度标签（`3e+05`…）自 x1253 起排；下界 870 因面板 y 轴题
`Sequences` 落在 x952..984 / y896..1057。→ 一块**横条**，正合适 7 根横条。

**数据（单样本平均 wall_sec，对数刻度）**：

| 轨道 | n | 均值 | 中位数 |
|---|---:|---:|---:|
| genomad | 395 | **235.2** | 188.9 |
| rdrpcatch | 395 | **157.1** | 57.2 |
| virbot | 395 | **31.0** | 14.1 |
| viralverify | 394 | **30.0** | 17.7 |
| metabuli | 395 | **14.2** | 11.8 |
| blast | 395 | **9.2** | 6.3 |
| viroid_blastn | 436 | **7.5** | 6.0 |

口径（⭐ 三个必须踩对的点，详见 `make_tool_timing_data.py` 头部注释）：

1. ⚠⚠ **`viroid_blastn` 根本不在 `_resource.tsv` 里**。源码 `virus_identification.py:1282`：
   ```python
   if has_viroid and "viroid" in target_tools:
       t0 = time.time()                       # ← 算了 t0
       run_viroid_identification(...)
       log_step("类病毒")                      # ← 却**从不** log_resource()
   ```
   其余 6 个工具都走 `run_timed()` 包装 → 有行；viroid 是**裸 `time.time()`** → **永无行**。
   → 用 **文件 mtime 差分**重建：`viroid_output/viroid_blastn.log`（= blastn 进入时创建、
   **全树 436/436 均为 0 字节**）到 `<sample>_viroids.result.id`（= 该步全部结束）的时间差。
   实测单样本 3–8 s（CRR1126132 = 7.3 s），全样本均值 7.5 s。
2. ⚠⚠ **同一 (sample, tool) 会有多行** —— `_resource.tsv` 是**逐次追加**的，行序 = 执行序
   （`genomad → blast → virbot → rdrpcatch → viralverify → metabuli`），断点续跑 / 二 pass
   在**末尾**再追加。→ **最终口径 = 每 (sample, tool) 取最后一行**。
   实测行数：virbot **851**（二 pass 普遍）、genomad 514、rdrpcatch 455、blast 399、
   viralverify 396、metabuli 396。⚠ 早期若用 `sort -u` 会**打乱追加序** → 必须按原序读。
3. ⚠ **样本分母**：有 `_resource.tsv` 的样本 = **395**；viroid 覆盖 **436**（⊇395，多出的 41 个
   是病毒支无产出、故未写 resource 的样本）。→ 两者 n 不同，图上按各自 n 计均值。

**画法**：`geom_rect` + `scale_x_log10`（7 条跨度 31 倍，线性会把后 4 条压成不可见 ——
用户确认用对数刻度）+ 条末数值标签 + 耗时降序排列。色用 OKABE_ITO 蓝 `#0072B2`
（与主图绿/灰/橙语义色**不冲突**）。

**⚠⚠ 渲染纪律**：本图**必须独立进程预渲**（`make_tool_timing_fig.R`）→ `tool_timing_fig.png`，
主脚本只 `readPNG` 贴。**同进程内 `ggsave` 这张 `theme_minimal` 图会给整张主图引入
18,173 px 的子像素重栅格噪声**（详见硬规则 16）。

**双端验收（2026-10-13，两种模式均通过）**：

| 项 | 本机 | 服务器 |
|---|---|---|
| 蓝条包围盒（框内） | x 574..1064 / y 285..796 | x 574..1063 / y 284..795 |
| 中心差 | — | **≤ 1 px** |
| 排除叠加框后差异像素 | **0**（vs 纯基线） | — |
| `Δink`（框外） | **+0** | — |
| 落位是否越界 | 否（目标 400..1245 / 280..870） | 否 |

→ 因为是**同一张预渲 PNG**（md5 `3b706e6086f848f2773ebdecadc56d6b` 双端一致），
两端落位差 ≤ 1 px。

⚠️ 部署纪律：服务器脚本先备份 `.bak_20261013_pre_timing`，再 `scp`；
新脚本 md5 双端一致（`ca11ee6431980c762c86fe4fb90a6f68`）；数据文件 md5 双端一致
（`tool_timing.tsv` = `d4d73634ea1a969ee020f5e74797cbd5`、
`tool_timing_fig.png` = `3b706e6086f848f2773ebdecadc56d6b`）。

### 第 8 次修订（用户 2026-10-13 晚间）：图例 → 右上角 + 耗时图行序与点阵对齐

用户原话：「你觉得位置对吗？**图例放右上角吧**」（并圈出左上传例与中部 Venn 右侧空白区）。

#### ① 耗时图行序改为「与右下点阵严格同序」

原先按**耗时降序**（genomad 最上）；现改为 **`TRACKS` 序**（= union 内命中数降序）：
`viralverify → metabuli → blast → genomad → virbot → rdrpcatch → viroid_blastn`，
与右下点阵 7 行**同一水平线**，便于上下对读。
⚠ 代价：条长不再单调递减 —— 但这正是「**命中多 ≠ 耗时高**」的可见证据
（blast 命中 95,573 只花 9.2 s；genomad 命中 11,526 却花 235.2 s = 全图最贵）。

⭐⭐ **本轮踩到的三个"改行序"专属坑**（全部零报错，务必记住）：

1. ⚠⚠⚠ **`as.numeric(factor)` 定位 ≠ 轴标签顺序**。改成 `factor(levels = rev(ORDER))`
   后若仍沿用 `y = as.numeric(tool)` 画条，而 `scale_y_continuous(breaks = 1:7,
   labels = ORDER)` 按**原序**贴标签 → **两套方向相反** → 条长与标签**错位整行**
   （实测图面 viralverify 显示成 7.5 s，真值 30.0 s）。**零报错、零警告**，
   必须逐行对着 tsv 核才看得出。
   → **正解**：显式 `ypos` 连续坐标，条与标签共用同一个 `ypos`。
2. ⚠⚠⚠ **`ypos` 的取向反相**。`scale_y_continuous` 的 y **向上增大** → `ypos = 7` 在**最上**。
   若直接 `ypos := match(tool, ORDER)`（viralverify→1），viralverify 会落**最下**、
   viroid_blastn 落**最上** —— 正好上下颠倒（实测：日志打印"上→下 = viralverify > …"，
   图上却是 viroid_blastn 在最上）。
   → **正解**：`ypos := N + 1 - match(tool, ORDER)`（反相），并加**自检**
   `stopifnot(最上那条 == ORDER[1])`。
3. ⚠ **`identical(sort(int), as.numeric(1:7))` 恒 FALSE**（整数 vs 双精度类型不同）→
   守卫用 `all(a == b)`。（此守卫**真的拦住了坑 1 的错位图**，别删。）

**纪律**：改 y 轴类版式后，**必须渲染出来逐行核对「标签 ↔ 条长」**，
不能只看「7 行都在、顺序对」就收工 —— 错位正是发生在"顺序对、配对错"的时候。

#### ② 图例 → 右上角（用 overlay 通道，不换 patchwork 格子）

⚠ 右上角那块"空白区"**不是 patchwork 的格子** —— 它在交集柱面板（x1490..3362 的
ggplot 面板）**内部**；`p_body$patches$plots[...]` 换不到（可见格子只有 4 个，全已占用）。
→ 沿用既有的、已验证**零 reflow** 的通道：**右上角图例 = 第三个 overlay**
（与耗时图 / Venn 同一套 `save_stacked_all(overlays = )` 机制）。

- **框**：`LEG_BOX = x2560..3320 / y435..640`（画布 px，760 × 205 px ≈ 2.53 × 0.68 in）。
  ⚠ 该框与 Venn 框（x1947..2905）**横向重叠 x2560..2905（245 px）** —— 两框在右上区
  **必然交叠**（该区是 Venn 右半与面板右缘之间的唯一留白带），无法靠挪框避开。
- **排布**：`legend.direction = "horizontal"` + **`guide_legend(nrow = 2)`**
  → 2 行 × 2 列。实测三种排布的自然墨迹（像素法）：
  horizontal 一行 = **3.19 × 0.10 in**（一行放不下 2.31 in 框宽）；
  **horizontal nrow = 2 = 1.82 × 0.23 in** ✓ 最贴合框比；
  vertical（4 行）= 0.93 × 0.33 in（竖长，浪费横向空间）。
- **左侧 spacer 格仍留占位**（`geom_blank` + `theme_void`）：若把该格清空，
  ComplexUpset 的 set-size 列宽可能变化 → **引发 reflow**（破坏"零回归"）。
  真正的图例改由 overlay 画在右上角（`--no-legend-ov` 可关）。

⚠⚠⚠ **叠加顺序坑（本轮最关键的一个）**：图例框与 Venn 框横向交叠 **x2560..2905**，
而 **Venn 的预渲 PNG 是白底**（`bg = "white"`）→ 若先画图例、后画 Venn，
**Venn 的白色矩形会把图例左半（"Retained"/"Removed by UniProt" 两行）整块盖掉**，
实测图面只剩右半两条 —— 且**零报错**（看着像"图例被裁"，实际是**被覆盖**）。
→ **正解 = 调整叠加顺序：Venn 在下、图例在上**（`overlays` 里图例排在 Venn 之后）。
⚠ 副产品：若将来 Venn 改**透明底** PNG，顺序不再敏感 —— 但**现在不可依赖**。

#### ③ 双端验收

| 项 | 本机 | 服务器 |
|---|---|---|
| 耗时图行序（上→下） | viralverify…viroid_blastn ✓ | 同 ✓ |
| 图例实测尺寸 | 1.77 × 0.23 in（框 2.53 × 0.68） | 同 ✓（字形宽 ~12% 仍余量充足） |
| 图例 4 条目全可见 | 是 | 是 |
| 未改动区 diff > 12 px | **0**（Venn 几何 / 右点阵 / 左 set-size / 标题带） | — |
| 改动区 diff > 12 px | 仅耗时图 / 图例旧位 / 图例新位 三处 | — |

→ 新增开关：`--no-legend-ov`（关图例 overlay）。
部署：两脚本先备份 `.bak_20261013_before_roworder` / `.bak_20261013_before_legendTR`，
再 `scp`；`make_tool_timing_fig.R` md5 双端一致（`ef390a04721bb53fd644c4c192ae5a9a`），
`fig_ident_filter_stacked_all.R` `diff` 双端为**空**（742 行逐字节一致，
⚠ 直接比对 md5 会因 CRLF 误报不一致 → 先 `tr -d '\r'`）。

### 第 9 次修订（用户 2026-10-13 20:25）：图例 → **Venn 右侧一列**

**用户原话**：「图例挡到 venn 了，一列右侧」「时间的条形图 在位置左上空白处」。

#### ① 问题定性

第 8 次修订把图例放 `x2560..3320 / y435..640`，当时的依据是「该带墨迹 = 0」——
但那次测量是在**上一版**（图例还在左上角）的干净基线上做的。等图例真搬过去，
**它自己就成了该带的内容**；而 Venn 的顶部墨迹（**x 到 2856、y 从 406 起**）一直在那
→ 实际图面上图例左段正压在 Venn 顶部两条标签上。靠"图例在下层之上"的叠加顺序
只是让它**可见**，视觉上仍是「挡到 Venn 了」。

→ 教训（已升为**硬规则 20**）：**每次挪 overlay，都必须重出"关掉该 overlay"的
干净基线**再扫描空档，**不能复用历史坐标**。

#### ② 做法：右侧列框 + 竖排一列

在干净基线（`--no-legend-ov --no-timing`，3420×2220）上实测：

| 量 | 值 |
|---|---|
| Venn 墨迹右缘 | **x = 2858** |
| 交集柱面板右缘 | **x ≈ 3362** |
| 右侧空档 | **`x2890..3350` 从 y0 到 y1470 全空**（唯一墨迹从 y1478 起 = 下方 x 轴条） |
| 取用框 | **`LEG_BOX = x2890..3350 / y430..1290`** = 460 × 860 px = **1.53 × 2.87 in** |

⚠ 扫描口径：`ink = (R<235)|(G<235)|(B<235)`；**逐 40 px 列带 + 逐 100 px 行带各出一遍**，
两向都为空才算真空（单方向扫会漏）。

图例本体改 **竖排一列**：

```r
guide = guide_legend(override.aes = list(alpha = 1), nrow = 4)   # ← 显式 nrow=4 = 一列 4 行
theme(legend.direction = "vertical",
      legend.key.size  = unit(0.30, "cm"),
      legend.text      = element_text(size = 7.0))
```

⚠ 为什么不缩小框去贴内容：框小 → 留白落在**框外**，仍可能压到东西；
**框取满该列、图例在框内居中**才能保证"这一列是干净的"。

| `key.size` | `text.size` | 实测自然尺寸 |
|---|---|---|
| 0.22 cm | 6.2 | 0.91 × 0.35 in |
| **0.30 cm** | **7.0** | **1.06 × 0.47 in** ✓ 采用 |
| 0.60 cm | 10.0 | 1.56 × 0.94 in（仍远小于框高，说明框有余量） |

#### ③ 叠加顺序：**仍然保持** Venn 在下、图例在上

挪框后两框**横向已不重叠**（2890 > 2858），但**顺序一字未改** —— 因为：
① 纵向仍重叠（y430..1290 vs Venn y406..1294）；
② 这是**零成本保险**：Venn 将来若再放大，顺序对了才不会静默盖掉图例。
⚠ 纪律：**挪框解决的是"当前这一版"，顺序正确性才是长效保险，两者都要做。**

#### ④ 双端验收

| 项 | 本机 | 服务器 |
|---|---|---|
| 图例实测尺寸 | **1.06 × 0.47 in**（框 1.53 × 2.87） | **同**（自检通过） |
| 图例排布 | 竖排 1 列 × 4 行 | 同 |
| 图例 4 条目全可见 | 是 | 是 |
| 图例与 Venn 是否相交 | **否**（2890 > Venn 右缘 2858） | 否 |
| 耗时图位置 | 左列中部大块（`x51..1156 / y545..1494` = **3.68 × 3.16 in**，比例 1.164）—— 第 13 次修订修正 Y 映射偏移后下移 121 px（沿革：第 9 次 `x400..1245 / y280..870` → 第 11 次下移 → 第 12 次放大 `x50..1162 / y424..1381` → 第 13 次下移，见硬规则 22b、24） | 同 |
| **未改动区** `|diff|>12` | **0**（对拍 `--no-legend-ov` 基线；差异仅落在 x2904..3349 / y669..1226 = 图例框内） | — |

→ 开关不变：`--no-legend-ov`（关图例 overlay）。
部署：备份 `.bak_20261013_before_legendR` 后 `scp`；
`fig_ident_filter_stacked_all.R` **753 行双端 `diff` 为空**（⚠ md5 直比因 CRLF 必误报 → 先 `tr -d '\r'`）。

## _check_discovery_overlap.py — 版式碰撞自动检测（通用）

劫持 `plt.Figure.savefig` → `fig.canvas.draw()` → 遍历每个 Text 图元
（轴刻度 / xlabel / ylabel / title / `ax.texts` / legend）的**显示坐标 bbox**
两两求交。剔除 `ax.axison == False` 的轴与超出 view limit 的刻度。
**判据：0 overlapping pairs**（`TOL = 0.5 px`），另报越出画布的字母。

```bash
cd scripts/barbarum
# discovery 主图（默认 module）
python _check_discovery_overlap.py --tol 0.5 --inrange-px 0 -- \
  --root <RNA-Lycium_barbarum_out> --outdir <.../S10_SCI_Discovery>
# 任意作图模块（--module）
python _check_discovery_overlap.py --module redraw_barbarum_detection_figs --tol 0.5 --inrange-px 0 -- \
  -v <best.summary.tsv> -m <metadata.tsv> -d <summary dir> --assembly-root <01_Assembly>
```

`--inrange-px 0` 收紧到「标签中心必须落在轴框内」。退出码 = 有重叠则 1。
判据与 `_check_main_overlap.py` / `_check_var_rich_overlap.py` 完全一致。

## 配色唯一真源：virome_discovery_pipeline/utils/sci_palette.py

2026-10-12 新建，与 phylo 侧 `virome_phylo_pipeline/utils/sci_palette.py`
**同源**（`OKABE_ITO` 顺序、`OKABE_BRANCH`、rcParams 键值、`SEMANTIC` 逐项比对一致）。
内容：`OKABE_ITO`(9 色) / `OKABE_BRANCH`(5 色) / `SEQ_CMAP=viridis` /
`DIV_CMAP=RdBu_r` / 语义色 `C_UP`…`C_EXP` / `apply_sci_style()` /
`sci_rcparams()` / `colors()` / `branch_colors()` / `figsize_for()` /
`despine()` / `sci_legend()` / `COL_W_2=4.58` / `COL_W_1=9.32`
+ **新增 `SEMANTIC`** 业务语义色表（verdict / genome / novelty / quality / stage，
全部取自 `OKABE_ITO`，无游离 hex）。

> ⚠ **必须按显式文件路径加载**，不能用 `from utils.sci_palette import`：
> 服务器 `~/mambaforge/lib/python3.10/site-packages/utils/` 存在**常规包**
> （有 `__init__.py`），优先级高于任意 `sys.path` 上的命名空间包
> `virome_discovery_pipeline/utils/`，**与路径顺序无关**地必然失败。
> ```python
> import importlib.util as _ilu
> spec = _ilu.spec_from_file_location("_sci_palette", cand)
> mod  = _ilu.module_from_spec(spec)
> sys.modules["_sci_palette"] = mod
> spec.loader.exec_module(mod)
> ```

## 02_filtering/metadata_association


**不用本目录脚本**，直接调用主流水线引擎：

```bash
export LD_LIBRARY_PATH=/home/zhangwenda/mambaforge/lib:$LD_LIBRARY_PATH
python virome_analysis_pipeline/virus_metadata_plot.py \
  -v .../02_filtering/high_conf.meta_adapter.tsv \
  -m .../01_data/sra_rna.data4/global_metadata/Global_Metadata_Final_with_Coordinates.tsv \
  -o .../02_filtering/metadata_association
```

调用链：`virome_analysis_pipeline/auto_known_virus.py:729–753`
（v2 引擎 = 仓库根 `virus_metadata_plot.py`；回退v1 = `utils/virus_metadata_plot.py`）。

引擎的分母设计是**整个 meta 队列（935）**，不是检出样本数——这是
docstring 第 1 条明示的意图（用真实流行率而非原始计数，因为公共数据里
各物种测序量严重不均衡）。v2 原版同样是这个口径，不要改成单物种分母。

元数据匹配率：250/262 检出 run 命中（95.4%），`ScientificName` /
`Location` / `Tissue` 100% 覆盖。

## HTML 报告（`virome_discovery_pipeline/report_pipeline.py`）配色

2026-10-12 起 chart 数据色统一走 `sci_palette`：
分类(无序)→`OKABE_ITO`、序数(有序)→`SEQ_CMAP=viridis`（`_ramp(n)` 3–8 档，
色值已固化字面量，本模块不依赖 matplotlib）、语义→`SEMANTIC`。
~60 处游离 hex（`#8B2C1F`/`#537D96`/`#4A6B4A`/`#ffa726`/`#00838f`/`#90caf9` …）已清零。

**有意保留**：CSS 主题变量与阶段徽章底色（`stage_defs` 第 5 列）——不承担数据编码，
且白字对比度约束下 Okabe-Ito 亮色（黄/橙/天蓝）不可用。

> ⚠ **跑 legacy 结果树必须 `export MMPV_IO_LAYOUT=legacy`**。`report_pipeline.py`
> 默认 **standard** 布局（`d_reports = 13_Report`），而枸杞结果树是 legacy
> （`10_Reports`）。裸跑会**静默**把降级报告（591 KB，正常 66 MB）写进新建的
> `13_Report/`，而真正的 `10_Reports/` 纹丝不动 —— 看起来"跑了但没效果"。
> 判据：看日志里 `Reports:` 那一行指向哪个目录。

> ⚠ **`taxonomic_sankey.py` 有两份副本**：`stats/taxonomic_sankey.py`
> （服务器，**未纳入 git**）与 `virome_discovery_pipeline/utils/taxonomic_sankey.py`。
> `generate_sankey()` 优先用 `stats/` 那份 → **只改 utils 不生效**。
> 两份 md5 必须一致。本次给两份都加了 `okabe_ito` 调色板
> （`COLOR_PALETTES` 字面量镜像 `sci_palette.OKABE_ITO`），
> 并把 Sankey 调用从 `--palette set3` 改为 `--palette okabe_ito`。

> ⚠ **核对 Sankey HTML 配色要用宽松正则**：plotly 把颜色序列化成
> `rgba(230, 159, 0, 0.595)`（**逗号后带空格**）。用 `#E69F00` /
> `rgb(230,159,0)` / `230,159,0` 三种 grep 全部 0 命中，会误判为"未生效"。
> 正确形式：`rgba?\(\s*\d+\s*,\s*\d+`。

## 尚未处理的遗留

`01_detection/summary/` 下仍有 5 个 symlink：

- `Fig3_dating_dispersal.png` — GCVA MCC 树 + skyline + BSSVS 迁移矩阵，
  数据源在 `09_phylo_results`，不属于 detection 阶段
- `Fig4_selection_interface.png` — dN/dS，正交选择
- `Fig5_recombination_variation.png` — 重组变异
- `Fig6_host_response.png` — GCVA vs PSTVd 组织嗜性 + GO/KEGG 富集，
  数据源在 `05_expression` 转录组
- `all_viruses.best.summary.tsv.bak_20260828_seg` — v2 未过滤版备份

前四张按决策**不重画**（跨阶段产物，硬套单物种口径会造假数据）。
建议移出 `summary/` 或加README 说明其真实数据源，避免误读为detection 产物。

**其余已知遗留**（详见 `reports/DISCOVERY_SCI_FIGS_20261012.md`）：

| # | 项 | 级别 |
|---|----|------|
| R1 | `report_pipeline.py` 无 `--io-layout` CLI 旗标，legacy 树只能靠环境变量；建议加"布局疑似不符"守卫 | 高 |
| R2 | `taxonomic_sankey.py` 两份副本不同步；建议纳入版本控制或加 md5 一致性告警 | 中 |
| R3 | `Fig_KnownVirome_Composite` 10.96 in 宽 → 通栏缩放 **0.850**（略低于 0.9 下限）；修法 = 画布收 9.32 + 字号 ×1.18 + 重跑 `_check_main_overlap.py` | 低 |
| R5 | `Fig3..Fig6` 软链口径（跨物种混合队列）与 `Fig1/Fig2`（445 单物种）矛盾，并排引用会造成分母口径冲突 | 中 |
| R11 | `10_Reports/filter_per_sample.tsv` **服务器上不存在**（仅本地生成）→ 服务器重跑报告链时该表会缺失 | 中 |
| R12 | id 文件口径 Σunion **556,950** vs `ident_summary.tsv::All_Candidate` TOTAL **556,773**，差 **177**；逐样本比对指向 `CRR1126132`（其 `_virus.all.result.id` 有 177 行，但 summary 记 0、`_virus.all.candidate.fasta` 仅 2 字节 = 空） | 低 |
