# reports/virome_analysis_pipeline — 宁夏枸杞 (Lycium barbarum) 专项绘图

> **⚠ 2026-10-14 迁移**：本目录原为 `MMPV-RNA/scripts/barbarum/`，已整体迁入
> `reports/virome_analysis_pipeline/`（用户指令：「直接在 reports 创建 virome_analysis_pipeline，
> 整理到这里，不是别的目录」）。同批并入的还有：
> `reports/sci_composite/` 的成品图（→ `figures/`）、其历史版本与临时产物（→ `_archive/fig_revisions/`）、
> 脚本 `.bak`/探针（→ `_archive/script_backups/`）、仓库根的 `_g2/`、`_t1/`、`Rplots.pdf`（→ `_archive/root_junk/`）。
>
> **目录约定**
>
> | 子目录 | 内容 |
> |---|---|
> | 本层 `*.py` / `*.R` | 绘图脚本 + 数据准备 + 版式体检（**扁平**，与原 `scripts/barbarum` 同级） |
> | `_assets/` `_append_src/` `identfilter_data/` `novelty_sankey_data/` `upset_data/` | 输入数据与预渲资产（脚本的**兄弟目录**） |
> | `figures/` | **成品图**（当前版本） |
> | `_archive/` | 历史版本 / 临时 / 备份（**只增不删**） |
>
> **路径不变性（为什么迁移几乎零改动）**：`reports/virome_analysis_pipeline/x.py` 与
> `scripts/barbarum/x.py` **距仓库根同为 3 层** → 所有 `Path(__file__).parents[2]`、
> `_HERE.parents[1]`、`REPO_DIR = dirname×3` 的仓库根定位**逐字未改**；数据目录仍是脚本的兄弟目录，
> `HERE/upset_data` 一类引用也未改。**仅 4 处硬编码**动过：
> `fig_novelty_sankey.py` / `fig_novelty_sankey_supp.py` 的 `DATA`（`ROOT/"reports"/"virome_analysis_pipeline"/…`）、
> `fig_host_prediction.R` / `fig_host_prediction_b_panel.R` / `fig_upset_classification.R` 的默认数据目录
> （改用 `script_dir()`，与 CWD 无关）。
>
> **零回归验证（2026-10-14）**：主图 PNG 与迁移前 rev.u **逐字节相同**
> （`985206d6e67e881106ef47c3b316fcbc`，3289 × 4280 px）；版式体检 **0 对重叠**；
> 24 个 `.py` 全部 `py_compile` 通过。
> ⚠⚠ **PDF 的 md5 不可作回归判据** —— matplotlib 内嵌 `/CreationDate`，同一脚本两次运行 md5 必不同
> （实测 `9eebc403…` vs `77ce559c…`，仅时间戳差 35 秒）；**判据只认 PNG**。
>
> **服务器镜像**：`~/MMPV-RNA/reports/virome_analysis_pipeline/`（2026-10-14 同步，`_archive/` 不上传）；
> 旧的 `~/MMPV-RNA/scripts/barbarum/` 已改名为 `scripts/_barbarum_old_20261014/` 留底。

---

# 原 `scripts/barbarum` README（正文，路径已批量改为新位置）

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
python reports/virome_analysis_pipeline/regen_barbarum_detection_report.py \
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
python reports/virome_analysis_pipeline/redraw_barbarum_detection_figs.py \
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
python reports/virome_analysis_pipeline/make_discovery_sci_figs.py \
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
脚本 md5：`make_discovery_sci_figs.py` 本地 **`0687a0a8d34738636c527839f5e553e9`**
（**v21.4：画布 `(10.0,21.48)` + 行权 `[0.5798,1.0,0.6978]` + 新面板 N（⑧ 整合）**）
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
⚠ **追加（2026-10-14 ⑧ 第 15 次修订）**：`fig_ident_filter_stacked_all.R` 本地已到
**`36ac4529295143193ad44b3ea81f9cea`**（行名搬到最左 + 死键修复 + `TIM_SHIFT_PX 246→113`）
→ **服务器仍是旧版**，须与 `make_discovery_sci_figs.py`(v21.4) 一并同步（先 `diff` 再 `scp`），
并**重跑 ⑧ 两版 + 主图 + 双端体检**才算闭环。
⚠ **同步 `fig_novelty_sankey.py` 时注意**：仓库里存在**并发会话**改过该文件的记录
（见 §「并发会话」）→ **先 `diff` 再 `scp`**，禁止整篇覆盖。
⚠ `Fig_GojiVirome_PhyloSynteny` 是**图片嵌入型**图（`AxesImage`，Text 极少）→
`_check_discovery_overlap.py` 对它**不适用**（它只测 Text bbox）；该图靠
`phylosynteny_manifest.json` 的 `height_capped` / `print_scale` 复核版式。

### ⭐⭐ v22（2026-10-14）：三张独立图整合进主图（3 行 → 6 行）

用户指令：「`Fig_GojiVirome_HostPrediction_BPanel.png` 放进主图 **以及宿主分类的那个条形图**」
+ 追加「`Fig_GojiVirome_UpsetClassification_Family.png` 也放进主图」。

| 新面板 | 内容 | 源图（独立图，R） | 裁标题块 | 行高 (in) |
|---|---|---|---|---|
| **O** | 宿主分类分布（9 类，n = 20,892） | `fig_host_prediction_a_panel.R`（**本轮新建**；与 `fig_host_prediction.R` 面板 a 同 ggplot 代码） | 283 px | 4.828 |
| **P** | 宿主预测器证据（耗时 / Venn / UpSet） | `fig_host_prediction_b_panel.R` | 185 px | 5.141 |
| **Q** | 科级分类一致性 UpSet（7 工具，n = 17,512） | `fig_upset_classification.R --rank family` | 212 px | 5.055 |

- `LAYOUT` 顺序 = `A+J / D / O / P / N / Q`（宿主预测两面板插在 D 与 N 之间，Q 垫底）；
  字母延续跳号传统（B/C/E/F/G/H 是历史下线 ID，不回收）→ O/P/Q 全新零冲突。
- 画布 `(10.0, 44.0)`、行权 `[0.5798, 1.0000, 0.7105, 0.7566, 0.6978, 0.7440]`；
  行高仍按「**资产比反解**」→ `ax.get_position()` 实测各贴图行**各向异性 ≤ 0.02 %**
  （D −0.007 % / O +0.006 % / P −0.013 % / N −0.008 % / Q −0.020 %）。
- 新主图（PNG）**3174 × 12610 px**（10.58 × 42.03 in，tight 裁剪后），
  md5（前 16）= **`577cdba2e8d78d5b`**（v22.2 配色协调后）；v22.1 含 Venn 版 = `a04df3e0aed4b344`；
  v22 无 Venn 版 = `1a1d9120c2b2daa6`（存 `.bak_20261014_pre_venn`）；
  v22.2 前 = `a04df3e0aed4b344`（存 `.bak_20261014_pre_color`）；
  旧版（3174 × 6202）存 `.bak_20261014_pre_v22`。
- 体检：`_check_discovery_overlap.py --composite-only` → **2 figures / 0 overlapping pairs**；
  A+J / D / N 三行与旧版**行高逐像素吻合**（1405/1406、1988/1987、1484/1484 px），
  仅因画布高度变化产生**亚像素重采样**差异（±2 px 平移扫描无法对齐 → 非内容改动）。
- ⚠⚠⚠ **新硬规则（本轮踩到，已并入 05 档案 README 第 6 条）**：
  `np.asarray(PIL图)` 前**必须 `convert("RGB")`** —— R 的 `ggsave` / `grDevices::png`
  （Windows GDI+）会写出**调色板 PNG（PIL mode='P'）**，此时 `np.asarray` 返回 **(H,W) 索引**
  而非 (H,W,3) → `imshow` 把它当**标量场**套上当前 `image.cmap`（本项目 = viridis）→
  **白底变亮黄、柱色全错、零报错**。实测：APanel = mode 'P'（中招）；
  BPanel / Family = mode 'RGB'（正常）。**判据** = 行内颜色直方图
  `viridis 黄 (253,231,37) == 0 px` 且 `OKABE_ITO 橙 (230,159,0) > 0`。
- ⚠⚠ **R 脚本的中文 `--outdir` 必须走 ASCII 中转**：Git Bash 把参数里的中文转码损坏
  （实测 `cannot create dir 'D:\f!i"'`）→ 先 `--outdir C:/tmp/<ascii>` 再 `cp` 回中文目录。
  （内部路径经 `script_dir()` 解析的**不受影响** —— 只有命令行参数中招。）
- ⏳ **服务器未同步**：`make_discovery_sci_figs.py`(v22, `78d98648ddef2e7a534891ca1e91df7d`)
  + 新建 `fig_host_prediction_a_panel.R`(`62c32dfe31cc95f80a22a9ad0dca4348`)
  + `fig_host_prediction_b_panel.R`(`139174b23a3e4ea351a0a0af7de1c5da`)
  + `fig_upset_classification.R`(`8100d0c277850d0587f49c4f050cb500`)
  → 须一并 `diff` → `scp`（**禁整篇覆盖**）后再双端各跑一遍体检。

### ⭐⭐ v22.1（2026-10-14）：Q 面板右上角补 7 集 Venn

用户指令：「Q 也绘制 venn 补充到右上角的空白处」，参考
<https://zhuanlan.zhihu.com/p/434441906>（R 语言 7 组韦恩图）。

**改动范围 = 只有 `fig_upset_classification.R`**（主图脚本零改动 → Q 行几何不变）。
脚本 md5 `657b727c7fcee5f1925f7aebd7debcc5`（原 `8100d0c277850d0587f49c4f050cb500`）。

- **框位置**（px @300dpi，画布 2796×1800）：`x 1390..2500, y 220..920` = **1110×700 px = 3.700×2.333 in**。
  - 判据：框内 `g<215` 的墨迹 **= 0**（`g<245` 的 13k 像素全是 4 条极浅网格线，`g<220` 即为 0）。
  - ⚠⚠ **必须落在 Q 的标题裁切线（`UF_CROP_TOP_PX = 212`）之下** —— 否则被主图**静默裁掉**。
  - 该框 **Species 与 Family 两版同时为空**（含 14 px 净空）→ 一套坐标两版共用。
- **实现**：`ggVennDiagram`（`shape_id="701"`，与 ③ BPanel 的 Venn **同族同风格**）→ `ggsave` PNG
  → `rasterGrob` → 叠在 `patchworkGrob` 上（**禁 `inset_element`**：对 patchwork 静默误算 npc）。
- ⚠⚠⚠ **7 集 Venn 标签三坑**（逐个实测）：
  1. 全标 127/77 个区域 → **完全糊死**（`venn` 包与 `ggVennDiagram` 两版皆不可读）；
  2. 只留质心 + 阈值过滤 → 数字仍**压花瓣线**；
  3. ggrepel 强 repel（`box.padding` 0.4）→ 标签被推离质心、**撞上集合名**。
  ✅ **正解 = `geom_label_repel`（`fill="white"` + `linewidth=0`）**：白底遮住穿过的轮廓线，
     `force=1.2` 轻度位移 + 细引线回指质心。
  ⚠ `stroke` / `stroke.colour` **不是 `geom_text` 参数**（静默忽略，只 warning）；
     `label.size` 在 ggplot2 ≥3.5 **已弃用** → 必须用 `linewidth`。
- ⭐ **标注阈值取「前 N 大」（`--venn-top`，默认 15）**，不用固定数值 ——
  Species（49 非空，top12=165）与 Family（77 非空，top12=302）分布差异大，
  固定阈值会让两版标注数差一倍；N=15 → Species 阈值 94 / Family 阈值 211，两版均 15 标签。
- ⚠ 中心区域（`1/2/3/4/5/6/7` 全票）实测 count = **0** → 不标注，无误导。
- **开关**：`--no-venn` / `--venn-top N` / `--venn-color #RRGGBB`（与既有 `--no-timing` 同族）。
- ✅ **零回归三重验证**（本轮实测）：
  1. **框外差异像素 = 0**（Family / Species 两版，改前 vs 改后对拍）；
  2. `--no-venn` 重生成的 Family 资产 md5 = `e013592e2e2d12c9ecdb6040a9fbedb6`
     **与 pre_venn 备份逐字节一致**；
  3. 用 v22 资产（无 Venn）重出主图 → md5 **精确回到 `1a1d9120c2b2daa688c26bccc2500d8f`**。
- **新主图**：PNG md5 = `a04df3e0aed4b34404bad2766c98d33f`（3,295,959 B）；Q 行各向异性仍 **0.0201 %**。
- ⚠⚠ **面板内 inset 的框必须落在该面板「标题裁切线」之下**（已并入 05 档案 README 硬规则 7）。
- ⏳ **服务器未同步**：`fig_upset_classification.R`(`657b727c7fcee5f1925f7aebd7debcc5`)
  + 两张 Q 产物（Species / Family 各 png+pdf）。

### ⭐⭐⭐ v23（2026-10-15）：主图改 **2 列 × 4 行** —— 从 1:4.4 长条到接近 A4

用户：「Fig_GojiVirome_Discovery.png 2*3 还是3*3」→ 经 AskUserQuestion 两轮选定
**2 列 × 4 行** + **面板宽锁 ~8.9 in**（字号零损失）。

**逐行**：行0 `A|J`、行1 `D|O`、行2 `P|N`、行3 `Q`（占左半栏，右下留空）。

| 口径 | 2 列时面板宽 | 字号 | D 桑基标签 |
|---|---|---|---|
| 画布宽固定 8.9 in（期刊双栏） | 3.85 in | **43 %** | 4.56 → **1.97 pt 不可读** ✗ |
| **面板宽锁 ~8.9 in**（采纳） | 8.9 in | **100 %** | **4.80 pt** ✓ |

⚠⚠⚠ **为什么不能「画布宽固定」**：7 个面板全是**横宽型**（资产比 1.3098–2.2350），
加列 = 面板变窄 = **字号等比缩小**。

⚠⚠⚠ **坑 1：`GridSpec` 等宽列必然形变**（首版踩，零报错）——实测 **J +35.2 % / D +33.5 % /
P +15.5 % / Q +13.5 %**。根因：两面板资产比不同（J 1.5525 vs A 2.2350），等宽必然拉伸一个。
→ 正解 = **显式坐标定位**（弃用 GridSpec）：行高 = 该行可用宽 / Σ(资产比)，
同行两面板**列宽不等、高相同**。判据 = 逐面板 `get_position()` 反推 → 形变 **全部 +0.00 %**。

⚠⚠ **坑 2：单面板面板字母 `TEXT-OFF-CANVAS`** —— O/P/Q 原走 `else` 分支（`dx=-0.34` 轴分数），
单面板画布 8.99 → 8.87–11.63 in 变宽后算出 **-386 / -294 / -299 px** 出界。
→ 并入 D/N 同档（`dx=0.0` 贴轴左缘 + `dy=1.12`）。

**几何实测**：画布 **20.00 × 26.28 in**（比 0.761）；PNG **6003 × 7673**（比 0.7824）。

| 面板 | 宽 in | 高 in | 比 | 形变 |
|---|---|---|---|---|
| A | 10.504 | 4.700 | 2.2350 | +0.00 % |
| J | 7.296 | 4.700 | 1.5525 | +0.00 % |
| D | 7.394 | 5.645 | 1.3098 | +0.00 % |
| O | 10.406 | 5.645 | 1.8432 | +0.00 % |
| P | 8.540 | 4.933 | 1.7313 | +0.00 % |
| N | 9.260 | 4.933 | 1.8771 | +0.00 % |
| Q | 8.686 | 4.933 | 1.7607 | +0.00 % |

**验证**：版式体检 **28 figures / 0 overlapping pairs**、`TEXT-OFF-CANVAS` **0**；
D 面板字高实测 **4.80 pt**（PASS）；内容主色与旧版一致。

**产物 md5**：`Fig_GojiVirome_Discovery.png` **`8591443bcd19cde3`**；
`.pdf` **`24207115ab1d2c00`**。旧版（v22.2）PNG `577cdba2e8d78d5b` 存 `.bak_20261015_pre_v23`。
**四处副本已同步**（review root / 桌面无前缀 / 桌面 `2_` / 档案库 `05_`）→ 逐字节一致。
⏳ **服务器未同步**：`make_discovery_sci_figs.py`（**先 `diff` 再 `scp`**）。

### ⭐⭐⭐ v24（2026-10-15）：主图新增面板 **R** = Alphaflexiviridae 物种比对

用户：「Alphaflexiviridae科的物种比对 单独一张，合并到主图，其他作为附图」。

**改动**：v23 的第 4 行原为**单面板**（Q 占左半栏、右下留空）→ 把 R 填进右下，该行变**双面板**。

**来源** = discovery ③ 的 09b 阶段产物
`<...>/_phylo/acvirus_trees/Alphaflexiviridae/Tree_Synteny_Composite.png`
（`virome_discovery_pipeline/utils/acvirus_tree_pro.py` 直出；归档副本在
`MMPV主图档案库_20261014/14_MLtree_科级树与共线性/单科面板/管线原图/Alphaflexiviridae/`）。

⚠⚠ **贴管线原图，不重绘、不裁标题块**（项目硬规则「整页独立图只能自绘」已被否决）：
- 与 O/P/Q 同源纪律 = `imread → imshow`，**不许重绘 / 降采样**。
- **与 O/P/Q 的差异**：O/P/Q 源图顶部 title+subtitle 与主图 `ax.set_title()` 重复 → 必须裁；
  R 顶部只有两个**子面板**标题（「Maximum-Likelihood Phylogeny」/「Comparative Synteny &
  Homology」），是图内元素，裁掉没法读 → **保留**。实测顶部空白仅 ~55 px（内容自 y≈73 起），
  无可裁余量。

⚠⚠ **副作用（用户已知悉并选定）**：行高由 v23 的「单面板沿用上一行高」改为「按 Σ资产比反解」
→ 行 3 高 `4.9331 → 5.8793 in`，**Q 面板由 8.685 → 10.352 in 宽**（等比放大 19 %，零形变）。
⚠ **行 0/1/2 完全不变（零回归）**；画布 `20.0 × 26.28 → 27.23 in`。

⚠ **资产比 = `4477 / 3534 = 1.266844`**（与 `_V23_RATIO["R"]` **必须同源**，
否则 `imshow(aspect="auto")` 静默形变 —— 同 v21.4 的 D 面板 13.2 % 压扁事故）。

**几何（实测，8 面板形变全部 +0.000 %）**：

| 面板 | 宽 (in) | 高 (in) | 面板比 | 形变 |
|---|---|---|---|---|
| A | 10.504 | 4.700 | 2.2350 | +0.00 % |
| J | 7.296 | 4.700 | 1.5525 | +0.00 % |
| D | 7.394 | 5.645 | 1.3098 | +0.00 % |
| O | 10.406 | 5.645 | 1.8432 | +0.00 % |
| P | 8.540 | 4.933 | 1.7313 | +0.00 % |
| N | 9.260 | 4.933 | 1.8771 | +0.00 % |
| Q | 10.352 | 5.879 | 1.7607 | +0.00 % |
| R | 7.448 | 5.879 | 1.2668 | +0.00 % |

**验证**：
- 版式体检 **30 figures / 0 overlapping pairs**；
- 面板字母定位 **10 个字母 / 0 off-canvas**（Q `x0=8.4`、R `x0=1113.3`，画布 2000 px 宽）；
- 8 面板各向异性 **全部 +0.000 %**。

**产物 md5**：`Fig_GojiVirome_Discovery.png` **`3b395c68d8daa09b`**（6027 × 7966 px）；
`.pdf` **`1ec0e8c59960bb3c`**。旧版（v23）PNG `8591443bcd19cde3` 存 `.bak_20261015_pre_v24`。
**四处副本已同步**（review root / 桌面无前缀 / 桌面 `2_` / 档案库 `05_`）→ 逐字节一致。
⏳ **服务器未同步**：`make_discovery_sci_figs.py`（**先 `diff` 再 `scp`**）。

### ⭐⭐⭐ v25（2026-10-15）：新增第 5 行 —— Potexvirus SDT + 覆盖度

用户：「Potexvirus SDT 和 序列的覆盖度图 呢 两个图继续补充到主图」。

**改动**：垫底新增**第 5 行** = **S | T**；行 0–3 **完全不变（零回归）**；
画布 `20.0 × 27.23 → 32.04 in`。

| 面板 | 内容 | 资产（裁后 px） | 资产比 |
|---|---|---|---|
| **S** | Potexvirus 成对同一性矩阵（SDT, NT） | `2327 × 2193` | **1.061104** |
| **T** | 代表序列覆盖度（reads depth vs genome position） | `2379 × 572` | **4.159091** |

⚠⚠ **S / T 均为贴管线原图**（`imread → imshow`），**裁掉源图标题块**：
- S：`POTEX_SDT_CROP_TOP_PX = 232`（源图 `2327 × 2425`，顶部**两行**粗体标题
  墨迹带 y 31–97 / 107–173）。
- T：同样裁顶部标题块（源图 `2379 × 617`）。
- 来源 = `09b_Analysis_Verify/SDT_matrix/Potexvirus/`（SDT 热图）
  与 Potexvirus 覆盖度产物；T **只放 1 条代表** = `SRR8191317_NODE_2`
  （9,471 bp / mean depth **117,291×**，即 6 条自有 Potexvirus contig 中平均深度最高者）。

⚠⚠ **本行的几何张力（重要）**：S 近正方（1.061:1）而 T 极扁（4.159:1）。
在「T 铺满 + 两面板等高 + 零留白」三条件下，`h = 17.80 / Σ比 = 17.80 / 5.220194 = 3.4098 in`
→ **S 只拿到 3.6182 in 宽**（≈ 1/5 栏）、T 拿 14.1818 in。
这是**几何唯一解**，不是取舍失误。⚠ 用户已在二选一（选项描述明写「SDT 被压到 3.55 in 宽」）
后选定「只放 1 条代表」。⚠ 若日后要放大 S，只能改「S 独占一行」或「T 不铺满（留白）」。

**几何（实测，10 面板形变全部 +0.000 %）**：

| 面板 | 宽 (in) | 高 (in) | 面板比 | 形变 |
|---|---|---|---|---|
| A | 10.504 | 4.700 | 2.2350 | +0.00 % |
| J | 7.296 | 4.700 | 1.5525 | +0.00 % |
| D | 7.394 | 5.645 | 1.3098 | +0.00 % |
| O | 10.406 | 5.645 | 1.8432 | +0.00 % |
| P | 8.540 | 4.933 | 1.7313 | +0.00 % |
| N | 9.260 | 4.933 | 1.8771 | +0.00 % |
| Q | 10.352 | 5.879 | 1.7607 | +0.00 % |
| R | 7.448 | 5.879 | 1.2668 | +0.00 % |
| **S** | **3.618** | **3.410** | **1.0611** | **+0.00 %** |
| **T** | **14.182** | **3.410** | **4.1591** | **+0.00 %** |

**验证**：
- 版式体检 **34 figures / 0 overlapping pairs**；
- 面板字母定位 **12 个字母 / 0 off-canvas**（S `x0=58.9`、T `x0=389.5`，画布 2000 px 宽）；
- 10 面板各向异性 **全部 +0.000 %**。

**产物 md5**：`Fig_GojiVirome_Discovery.png` **`f6a8e3725c3f5a00`**（6027 × 9386 px）；
`.pdf` **`6195a0baba62f359`**。旧版（v24）PNG `3b395c68d8daa09b` 存 `.bak_20261015_pre_v25`。
**四处副本已同步**（review root / 桌面无前缀 / 桌面 `2_` / 档案库 `05_`）→ 逐字节一致。
⏳ **服务器未同步**：`make_discovery_sci_figs.py`（**先 `diff` 再 `scp`**）。

### ⭐⭐⭐ v26（2026-10-15）：底行补 **SDT 折线图**（新面板 U）+ 两处列宽对齐

用户：「SDT的折线图图呢？ 这样E和F才能平分一行，不显得F过大 另外D 和 A 长度对齐，J和O对齐」。

> ⚠ **「E / F」是历史面板编号**（v20 之前的主图字母），本图现无 E/F 面板。
> 经三轮二选一确认，用户实际指的是**底行 S / T**（「底行 S 和 T」）、
> 并要求「**零形变 + 行内留白**」「**三图一行**」「**行 2 换序 O\|D**」。

**改动 A — 新增面板 U（SDT 成对同一性「分布」折线图）**

| 面板 | 内容 | 源资产 | 裁后 px | 资产比 |
|---|---|---|---|---|
| **U** | Potexvirus 成对同一性**分布**（SDT, NT，折线） | `SDT_Potexvirus_distribution.png` `3025 × 1711` | `3025 × 1578` | **1.916983** |

- 来源 = 09b `01_data/SDT_matrix/barbarum__Potexvirus/SDT_Potexvirus_distribution.png`
  （与面板 S 同一 27-taxon pool）；落到 `S10_SCI_Discovery/Fig_GojiVirome_Potexvirus_SDTdist.png`。
- **裁标题块**：源图粗体标题墨迹带实测 **y 30–96**（单行
  `Distribution - Genus Potexvirus: 8 user contigs + 19 db refs (09b ACVirus) - NT`），
  下一处墨迹 = 坐标轴**顶脊 y 143–145**（1865 px 宽实线）→ `POTEX_DIST_CROP_TOP_PX = 133`。
  ⚠ 与 S 的 232 px 不同源：S 是**两行**标题、U 是**单行**。
- 结果：**底行由 `S\|T` 改三图 `S\|U\|T`** → T **14.18 → 10.09 in**（−28.8 %，「不显得 F 过大」达成）、
  S 3.62 → **2.58 in**、U = **4.65 in**。⚠ S 变窄是用户明示接受的代价。

**改动 B — 行 1 换序 `D|O` → `O|D`（列宽对齐）**

| 行 | 面板 | 分栏线占比 | 分栏线绝对 x |
|---|---|---|---|
| 0 | `A \| J` | 10.5038/(10.5038+7.2962) = **59.01 %** | 0.86 + 10.5038 = **11.3638 in = 3409.1 px** |
| 1（换序后） | `O \| D` | 1.84324/(1.84324+1.30979) = **58.46 %** | 0.86 + 10.4058 = **11.2658 in = 3379.7 px** |

⇒ 绝对差仅 **29.4 px（0.098 in）** / 6027 px 画布 → **零留白、零形变、无任何面板改尺寸**
（O 仍 10.4058 in、D 仍 7.3942 in）。代价 = D 移到右栏（与 J 对齐）、O 到左栏（与 A 对齐）。
⚠ 换序**只改 `_V23_ROWS` 的第 2 行元素顺序**，`_V23_RATIO` / 行高反解公式零改动。

**几何（v26 实测，画布 20.0000 × 31.0541 in，11 axes，形变全部 +0.000 %）**：

| 行 | 行高 (in) | 面板宽 (in) |
|---|---|---|
| 0 | 4.6997 | A 10.5038 / J 7.2962 |
| 1 | 5.6454 | **O 10.4058 / D 7.3942** |
| 2 | 4.9330 | P 8.5404 / N 9.2596 |
| 3 | 5.8793 | Q 10.3518 / R 7.4482 |
| 4 | **2.4267** | **S 2.5750 / U 4.6520 / T 10.0930** |

⚠⚠ **`_v23_geometry()` 未改**：`h = (usable_w − col_gap×(n−1)) / Σ比` 天然支持三面板行（自动扣 2 个 gap）。

⚠⚠ **三面板行的独立图不缩水**：`S` / `U` / `T` 的**单面板成品图**画布高仍锁 `H = 3.7887 in`
（`S = 3.6548` / `U = 6.6026` / `T = 14.3250`）—— 主图第 4 行虽变矮，**独立交付图尺寸不变**，
**只要求比例与主图面板一致**（尺寸自由）。

**验证**：
- 版式体检 **36 figures / 0 overlapping pairs**（v25 为 34）；
- 面板字母 **11 个字母 + 2 个刻度误报（"63" / "45"）/ 0 off-canvas**；
- 11 面板各向异性 **全部 +0.000 %**（O 有 +0.007 % 舍入）；
- 面板字母档：`O/P/Q/R/S/T/U` 全部 `panel_letter(ax, L, dx=0.0, dy=1.12)`
  （`dx=-0.34` 会算出 O −386 / P −294 / Q −299 / R −206 / S −116 / U −209 px **出界**）。

**产物 md5**：`Fig_GojiVirome_Discovery.png` **`0e33ce6180fc93799b7c8550528a414d`**
（**6027 × 9083 px**，4,193,102 B）；`.pdf` 3,506,107 B（md5 每次不同）。
旧版（v25）PNG `f6a8e3725c3f5a00` / 6027 × 9386 存 `.bak_20261015_pre_v26`。
**四处副本已同步**（review root / 桌面无前缀 / 桌面 `2_` / 档案库 `05_`）→ 逐字节一致。
A/B 对比图存 `05_.../过程留档/主图布局v26_改前改后对比.png`（md5 `8cc3552eb84c8ad6`）。
⏳ **服务器未同步**：`make_discovery_sci_figs.py`（**先 `diff` 再 `scp`**）。

### ⭐⭐ v22.2（2026-10-14）：三张 UpSet 配色协调 —— 类别色 vs 数量色分离

用户：「三个 upset 的颜色可以调整下嘛？尤其是左侧的条形图」+「也一并协调顶部交集柱」。

**采纳方案**（经 AskUserQuestion 两轮确认）：左侧 set-size 条 = **身份色**；
顶部交集柱 = **三图统一的数量色阶**；N 保持不动。

| 图 | 左侧条（改后） | 顶部柱（改后） |
|---|---|---|
| **P** `fig_host_prediction_b_panel.R` | `HOSTPRED_COLORS` 四预测器（与左上耗时 inset 同色） | `QUANT_COLORS` 1–4 档 |
| **N** `fig_ident_filter_stacked_all.R` | **不动**（彩色堆叠 = 过滤结局） | **不动**（类别维度，改单色会丢结论） |
| **Q** `fig_upset_classification.R` | `TOOL_COLORS` 七工具 | `QUANT_COLORS` 1–7 档 |

⚠⚠⚠ **为什么必须分离**（改前实际撞色）：Q 左侧若用 `TOOL_COLORS`、顶部仍用
`DEGREE_COLORS` → **4 个完全同色**（`#56B4E9`/`#F0E442`/`#0072B2`/`#D55E00`）；
P 撞 2 个（`#0072B2`/`#D55E00`）→ 图例串味。⇒ 正解 = **类别用彩色、数量用单色系渐变**。

⚠⚠ **数量色阶跨图同值同色**：三图分别 1–7 / 1–4 档，各自 `viridis(n)` 插值会让同一个
「2」在不同图里颜色不同 → 用固定档位表 `QUANT_COLORS[k]`（`upset_common.R`）**按档位索引**。

**新增共用件**（`upset_common.R`）：`QUANT_COLORS` / `QUANT_BLUES` / `quant_colors()` /
`set_size_colors()` / `set_size_colors_filtered()`。

⚠⚠⚠ **三条硬坑（全部实测，已并入 05 档案 README 硬规则 9–11）**：
1. **左侧条只能用颜色向量**：`aes(fill = ...)` 三种写法在 ComplexUpset 1.3.3 全报错
   （`object 'x'/'set' not found`、`Cannot add <ggproto> objects together`）；
   ✅ `geom_bar(fill = <向量>)`，且**向量按传入 `sets` 顺序**（与 `sort_sets` 渲染顺序无关）。
   ⚠ 条宽保持 **0.6**（写 0.62 → 左侧 3,855 px 差异）。
2. **`min_size` 静默裁集合 → 向量长度必须过滤**：Species 的 `genomad` 仅 13 条
   vs `min_size=68` → 被裁 → 7 色 vs 6 行 → `Aesthetics must be either length 1 or the
   same as the data (6)`。✅ 用 `set_size_colors_filtered(sets, pal, sig_top$sig)`。
   ⚠⚠ **只在 Species 触发（Family 的 7 工具全在）→ 必须两个 rank 都跑**。
3. **色阶方向**：viridis 1→7 = 深紫→黄（**大 = 亮**），与旧 `DEGREE_COLORS`
   （1=浅黄…7=橙红，**大 = 深**）**方向相反** → 属用户可见变更，题注已写明。

**零回归三重验证**（旧口径 `--plain-setsize --quant-scheme degree` 重跑，全部逐字节）：

| 图 | 重跑 md5 | 与改前 |
|---|---|---|
| Q Species | `a583f631ffe13d9a` | ✓ 完全一致 |
| Q Family | `f0ff728a1b29beff` | ✓ 完全一致 |
| P BPanel | `3f2dd9e5dad6d664` | ✓ 完全一致 |
| N IdentFilter_norm | `e1f7174c1e1c176c` | ✓ 未受影响 |

**新版差异分区实测**：P / Q 差异**只**在「左侧条 + 顶部柱」，**Venn / 点阵 / inset = 0 px**。

**产物 md5（v22.2）**：
- 主图 PNG `577cdba2e8d78d5b904eec4155a8da5f`（3,174×12,610，尺寸与 v22.1 一致 → 几何零位移）
- Q Family `17856b5d0d650cf0cb9b9f8431980bbf`｜Q Species `86ee82428a0ea5d89a2487e9b17d7b3c`
- P BPanel `84c8c1c63fb72c890aeae1fa83773073`
- 版式体检 → **2 figures / 0 overlapping pairs**
- 脚本 md5：`upset_common.R` `36f485f3…`→（改后见下）；`fig_upset_classification.R` `657b727c…`→改后
- **四处副本已同步**（review 根 / 桌面 / 档案库 05+09；写作目录无副本），旧版备份 `.bak_20261014_pre_color`
- ⏳ **服务器未同步**：`upset_common.R` + `fig_upset_classification.R` + `fig_host_prediction_b_panel.R` + P/Q 产物

**验收对比图**（`3×2` 布局，2026-10-14）：
- 主网格 **3 列 × 2 行**（列 = Q / P / N，行 = 上改前 v22.1 / 下改后 v22.2）
  + 左侧 set-size 条形图**放大 ×2.5** 细节带（同样 3×2）+ 底部配色图例
- 三张源图宽高比均 **1.54–1.55** → 3 列并排整幅接近 4:3；且改动核心是**三图之间的协调**，
  必须**横向并排**才看得出「左侧条彼此区分 + 顶部柱三图同色阶」
- 归档 `MMPV主图档案库_20261014/09_UpSet_工具分类宿主/过程留档/三张UpSet配色v22.2_验收对比_3x2.png`
  （2,970×2,480）+ 排版脚本 `_三张UpSet配色v22.2_排版脚本.py`
- ⚠ 细节带裁框**必须含左侧行标签**：Q/P = `(30,1130,1360,1670)`、N = `(20,1508,1500,2090)`。
  首版左边界取 `100` → P 列行名被切成 `TV / VH / B2 / ule`（零报错）
- ⚠ **判据坑**：`TOOL_COLORS` 与 `DEGREE_COLORS` 同源 OKABE_ITO → 直接数「旧色残留像素」
  会**假阳性**（`#0072B2`/`#D55E00`/`#F0E442` 两套都有）→ 正解 = **A/B 逐像素差异分区**

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
python reports/virome_analysis_pipeline/make_goji_virome_phylo_figs.py \
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
python reports/virome_analysis_pipeline/make_goji_virome_phylosynteny_fig.py \
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
python reports/virome_analysis_pipeline/fig_discovery_supp.py --root <RNA-Lycium_barbarum_out>
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
python reports/virome_analysis_pipeline/make_discovery_novelty_figs.py \
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
"C:/Python312/python.exe" reports/virome_analysis_pipeline/make_novelty_sankey_data.py
"C:/Python312/python.exe" reports/virome_analysis_pipeline/fig_novelty_sankey.py --outdir <...>/S10_SCI_Discovery
"C:/Python312/python.exe" reports/virome_analysis_pipeline/fig_novelty_sankey_supp.py --outdir <...>/S10_SCI_Discovery  # 附图
"C:/Python312/python.exe" reports/virome_analysis_pipeline/_check_discovery_overlap.py --root <...> --outdir <...>      # 版式体检
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

### ⭐⭐⭐ 2026-10-14 00:0x「`Fig_GojiVirome_IdentFilterStackedAll_norm.png` 整合到主图」→ v21.4

**用户指令**：把 ⑧（`Fig_GojiVirome_IdentFilterStackedAll` 的 **`_norm`** 变体）整合进主图。
**澄清两问**：① 版式 → **「底部通栏 · 主图加宽」**；② ⑧ 自带标题块 → **「裁掉」**。

#### ⚠⚠⚠ 为什么必须「加宽主图」而不是「把 ⑧ 缩窄」

| 事实 | 值 |
|---|---|
| ⑧ 源图 | **11.4 × 7.4 in**（`W_FIG`/`H_FIG`），比 1.5405 |
| 主图旧可用宽 | **7.61 in** → ⑧ 必被缩 **0.667×**（最小字 6.4 pt → **4.3 pt**，糊） |
| ⑧ 能否重渲成窄版 | **不能** —— overlay 坐标是按 `W_FIG = 11.4` 写死的**分数**（如 `leg_w_in <- (LEG_R-LEG_L)*W_FIG`），改 `W_FIG` 会让预渲 inset 的**物理尺寸跟着变** → Venn / 图例 / 耗时条 / 序列标题四处 overlay 全部错位 |

⇒ 唯一出路 = **把主图加宽**。画布宽 **8.55 → 10.0 in**（可用宽 **7.61 → 8.90 in**）→ ⑧ 缩放 **0.781×**（最小字 6.4 → **5.0 pt**）。

#### ⚠⚠ 裁掉标题块：裁切线是「量」出来的，不是拍的

| 量 | 值 |
|---|---|
| 标题块（1 行主标题 + 4 行说明） | y **58..253 px**（第 4 行底 = 253） |
| **最高内容元素 = 右上 Venn** | 顶边 y = **407 px**（1.357 in） |
| 各 x 带「首个内容块」 | 左(耗时条) **573**、中(交集柱) **522**、右(Venn) **407**、最右(图例) **787** |
| ⇒ 裁切线 | **y = 398 px**（Venn 之上留 9 px 净空）→ 裁后 **3420 × 1822 px = 11.4 × 6.073 in（比 1.87706）** |

⚠⚠ **裁切在内存里做**（`panel_N_identfilter` 里 `Image.crop`），**不落中间资产**
→ 源图重出后自动跟随，不会留陈旧副本（本项目吃过「预渲资产不随参数更新」的亏）。
⚠⚠ **裁切线与 `H_FIG = 7.4` 强耦合** → 代码里加了**高度自检**：源图高 ≠ 2220 px 时
**回退为「不裁」并打印告警**，绝不静默裁错位置（会把 Venn 顶切掉或留半行标题）。

#### ⭐⭐ 行权是「反解」的 —— 让两块贴图**各向异性 = 0**

`imshow(aspect="auto")` 会把资产**强行拉伸**到轴框 → **面板比 ≠ 资产比 = 静默形变**。
故两个通栏行的高度必须由**资产比反推**，不能拍：

```
D 行高 = 可用宽 / 桑基资产比 = 8.900 / 1.30979 = 6.795 in
N 行高 = 可用宽 / ⑧ 裁后资产比 = 8.900 / 1.87706 = 4.741 in
上排   = 3.940 in（A/J 不动）
合计 15.476 in ×(1 + 2×hspace/3) = 18.366 in → 画布高 = 18.366/0.855 = 21.48 in
```

| 项 | v21.3 | **v21.4（现行）** |
|---|---|---|
| `COMPOSITE_FIGSIZE` | `(8.55, 14.0)` | **`(10.0, 21.48)`** |
| `LAYOUT_ROW_RATIOS` | `[0.6, 1.0]` | **`[0.5798, 1.0000, 0.6978]`** |
| `LAYOUT` | A(11)+J(9) / D(20) | **A(11)+J(9) / D(20) / N(20)** |
| 面板字母 | A, D, J | **A, D, J, N**（新字母 **N**，不复用下线的 B/C/E/F/G/H） |
| A / J 面板 | 4.120×3.938 / 3.344×3.938 in | **4.819×3.940 / 3.912×3.940 in**（宽了） |
| **D 面板** | 7.610×6.562 in（比 **1.1595**） | **8.900×6.795 in（比 1.30979）** |
| **N 面板** | — | **8.900×4.741 in（比 1.87706）** |
| **D 各向异性** | **13.2 %（横向压扁）** | **0.05 %** ✅ |
| **N 各向异性** | — | **~0 %**（实测 sy = 1373/1759 = 0.7805） |
| 主图 | 9.28×13.80 in | **10.58×20.67 in** |

#### ⚠⚠⚠ 顺手修掉一个存量缺陷：D 一直被横向压扁 13.2 %

v21.3 及以前 **D 面板比 1.1595 ≠ 桑基资产比 1.30979** → `aspect="auto"` 把桑基
**横向压扁 13.2 %**。**实测证据**（不靠推算）：资产内 Known 块 277 px、主图内 364 px
→ `scale_y = 1.3141`；而 `scale_x = 7.610/6.553 = 1.1612` → 各向异性 13.2 %。
本轮按「面板比 = 资产比」反解行高 → **归零**。
⚠ **副作用（须告知用户）**：D 的观感会**略变宽**（不再被压扁），这是修正不是回归。

#### 实现要点

- 新面板函数 `panel_N_identfilter(ax, D)`：`_find_asset` 找 `..._norm.png` → 内存裁顶 → `imshow(aspect="auto")` + `axis("off")` + 标题/题注 + 写 `manifest["N"]`。
- **字母偏移复用 D 的那套**（`elif letter in ("D", "N")` → `dx=0.0, dy=1.12`）：两者同为「轴占满通栏 + `loc="left"` 标题」，另起一套必掉回标题行。
- `GridSpec` 行数改用 **`len(LAYOUT_ROW_RATIOS)`** 而**不写死 2** —— 历史踩坑：行号与 `height_ratios` 下标错位会**零报错**地让面板拿到错的权（末行空白被 `bbox_inches="tight"` 裁掉）。
- 单面板 D 同步改为 **`(8.99, 7.55)`**（axes = 8.900 × 6.795）、新增单面板 N **`(8.99, 5.268)`**（axes = 8.900 × 4.741）—— 组合图与单面板必须**同源**（硬规则 87）。

**验收**：主图体检 **22 figures / 0 重叠**、**0 处 `TEXT-OFF-CANVAS`**；桑基体检 **0 重叠 / 无脱画布**。
交付 md5：主图 **`a1df6acd3f3c5343a1ec2a624216df42`**。
⚠ ⑧ 自身（`Fig_GojiVirome_IdentFilterStackedAll{,_norm}`）**未改动**，桌面副本不动。

### ⭐⭐⭐ 2026-10-14 01:0x「这个行名放最左边可以吗」→ ⑧ 第 15 次修订（含"死键"缺陷修复）

**用户指令**：红框圈住**图中间**那 7 个行名（`viralverify (452,885)` … `viroid_blastn ( 604)`）→ 放到**最左边**。

#### 根因（两处；其中一处是**静默死键**）

1. 行名本来由**矩阵面板**画：`upset()` 的 `intersections_matrix` 面板 `axis.text.y` 出集合名，
   该面板在 patchwork 里位于**中间列** ⇒ 名字必然落在「集合条右端」与「点阵」之间。
2. ⚠⚠⚠ **旧版主题键名写的是 `"intersections"` —— 一个永不生效的「死键」**。
   `names(ComplexUpset:::upset_themes)` 实测**只有 4 个合法键**：
   `intersections_matrix` / `Intersection size` / `overall_sizes` / `default`；
   而 `upset_modify_themes()` 是**按名合并、对不认识的 key 只附加、不报错**
   ⇒ 旧版那句 `axis.text.y = element_blank()` **从未生效**。
   证据链：① `upset_themes[["intersections"]]` 为 NULL；② 合并后 `names()` 里
   `"intersections"` 与 `"intersections_matrix"` **并存**；③ A/B 只改那个键 → 图中名字**纹丝不动**。
   ℹ ⑨（`fig_upset_classification.R`）注释里记着**同一个坑** —— 两图同源同病。

#### 修法

```r
themes = upset_modify_themes(list(
  # 矩阵那份不画名字（必须拼回内置两段！见下）
  "intersections_matrix" = Reduce(`+`, ComplexUpset::upset_themes$intersections_matrix) +
    theme(axis.text.y = element_blank()),
  # 名字只保留「最左一处」= set-size 面板
  "overall_sizes" = theme_upset(8) + theme(axis.text.y = element_text(size = 7.5, colour = "black"), ...)
))
```

- ⚠⚠⚠ **必须 `Reduce(`+`, ComplexUpset::upset_themes$intersections_matrix) + theme(...)`**：
  内置值是**两段 list**（`[[1]]` = 144 元素完整主题、`[[2]]` = 3 元素微调），而 `themes[[key]]`
  是**整体替换**语义。写成 `theme_upset(8) + theme(...)` 会踩两坑：
  ① 丢 `[[2]]` 的 `axis.text.x` blank → 矩阵底部冒出 **16 列 × 7 集合名挤成一团的黑字**
  （还把 "Identification track" 题注挤花，**零报错**）；
  ② `theme_upset(8)` 的 `axis.text` 同时打开 `.x/.y`，把行名**又放回来**。
- ⚠ 翻转后主题元素名指**渲染后**的轴：`axis.text.y` = 左侧垂直轴（本面板 = 集合名）、
  `axis.text.x` = 底部水平轴（数值）。标定法：`coord_flip()` 里把 `.x` 设红 / `.y` 设蓝 → 渲染后蓝在左、红在底。

#### ⭐⭐ 与 ⑨ 的正解**恰好相反**（勿照抄 ⑨ 的 `transparent`！实测 A/B，画布 3420 px）

| | 最左面板原本画名字？ | `element_blank()` | `colour="transparent"` |
|---|---|---|---|
| **⑨** | **是** → 只需把矩阵那份「藏起来」 | 文字槽**凭空消失**、矩阵左移拉宽（柱1 中心 1103→999）→ 用户否决 | ✅ 正确 |
| **⑧** | **否** → 必须**新开**一条名字槽（实测 284 px） | ✅ **搬移**：col2 让出 284 → col1 用掉 299，总宽**近似守恒** | ✗ **复制**：col2 保留 284 且 col1 再占 299，总宽凭空 **+284** |

实测（矩阵点阵首列 / 列距）：基线 `1407 / 124.2`；`element_blank` **`1434 / 122.5`（Δ+27、−1.37 %）**；
`transparent` `1610 / 111.2`（Δ+203、−10.5 %）⇒ **⑧ 取 `element_blank()`**。
残余 +27 px 无法在「不改画布宽」下消掉；已评估加宽画布 ~300 px（⑧ → 12.4 in ⇒ 面板 N 降采样
0.781→0.721、6.4 pt 副标题降 4.6 pt < 5 pt 下限）⇒ **放弃**，残余在主图上 ≈19 px 不可辨。

#### ⭐⭐ 耗时图 overlay 重新对齐（口径**改过一次**，记录清楚）

行名搬到最左后集合条面板整体右移（实测 `166..905` → `465..1194`），而耗时图是**固定画布分数** overlay ⇒ 必须同步右移。

| 口径 | `TIM_SHIFT_PX` | 耗时名字右缘 | 集合条名字右缘 | 名字差 | 绘图区左缘差 | 绘图区右缘差 |
|---|---|---|---|---|---|---|
| 先做（对齐最长条 465 px） | `246` | 454 | 330 | **124 px ✗** | 123 px | **3 px** |
| ✅ **正解（对齐标签列）** | **`113`** | **321** | **330** | **9 px ✓** | **10 px ✓** | 136 px |

⇒ 取 **`TIM_SHIFT_PX <- 113`**：**标签列 + 绘图区左缘同时对齐**（对上两项）。
- ⚠ 曾按 246 交付并全过自检/体检，**是目视才发现 124 px 标签台阶** —— 正是用户这次要修掉的东西。
- 代价（已知并接受）：耗时图绘图区右缘（x=300）1061 vs 集合条 1197 → 短 136 px。
  两者 x 轴标度本就无关（log 秒 vs 命中数）；耗时图是**独立预渲 inset**，比主面板窄属正常。
  ⚠ **若强行同宽须改 `make_tool_timing_fig.R` 的 `W_IN`（绘图区加宽 ~126 px）—— 但用户在第 13 次修订里
  亲自用红框钉死了这块面板的尺寸（「再次调整大小，符合红色区域的大小」）⇒ 不擅自改其绘图区宽度。**
- 副产物（正面）：框内容右缘 1221 → **1088**，与「Sequences」题注左缘 1229 拉开 **141 px**，
  消除旧版「只剩 8 px」的擦边风险。
- ⚠⚠ **历史坑（守卫保留）**：第 13 次修订曾取右移 280 px → 框内数值标签右缘到 1255 px，
  而「Sequences」题注是**白底不透明 raster**，必须夹在它与交集柱 y 轴刻度左缘 1281 之间
  ⇒ 只剩 26 px < 8 pt 旋转字高（≈33 px）→ **题注白底把 "235.2" 横切成 "23"**（实测）。
  ⇒ **不透明 raster overlay 必须做「邻居净空」自检**（已内置 `TIM_CONTENT_R_PX` + 题注自检两条 `cat`）。

#### 交付 / 验收（全部实测）

| 产物 | md5 |
|---|---|
| `fig_ident_filter_stacked_all.R` | **`36ac4529295143193ad44b3ea81f9cea`**（`R 语法 OK`） |
| `Fig_GojiVirome_IdentFilterStackedAll.png` | **`44d989732d364d168316e62a6b65c077`** |
| `Fig_GojiVirome_IdentFilterStackedAll_norm.png` | **`e1f7174c1e1c176c577aacbad0488e56`** |
| `Fig_GojiVirome_Discovery.png`（主图，panels=AJDN） | **`0be7a55376ec7ee2ea2fdaf9936fa208`** |

- ⭐ 规范产物与 `/tmp` 试渲**逐字节一致**（`e1f7174c…` 双端同 md5）。
- 主图与改前备份逐像素对拍：差异 **仅 y 4797..6081（面板 N）**，**A/J/D 行 0 差异**；差异像素 2.33 %。
- 版式体检：`exit=0`、**22 figures / 0 overlapping pairs**、**0 处 `TEXT-OFF-CANVAS`**、主图 `texts=65`。
- **桌面同步**（`/d/桌面/MMPV管线主图_20261013/`）：`Fig_GojiVirome_IdentFilterStackedAll{,_norm}.{png,pdf}`
  + `2_Fig_GojiVirome_Discovery.{png,pdf}` → 覆盖后 md5 与规范一致 ✓。
  ⚠ 核对发现**桌面旧版 = 规范「改前」备份**（逐字节一致）⇒ 桌面此前**从未同步过 v21.4**。
  `6_…_NoveltySankey_Supp` / `9_…_Discovery_Supp` md5 **未变**（无需同步）。
- ⭐ 顺带清掉一个**陈旧同名副本**：`reports/_discovery_sci_20261012/Fig_GojiVirome_Discovery.{png,pdf}`
  （3361×4208 = 旧版式，mtime 10-12 23:39，全仓库**无引用**）→ 备份 `.bak_20261012_oldlayout` 后刷新为新主图。

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
Rscript reports/virome_analysis_pipeline/fig_upset_tools.R                       # 图①
python  reports/virome_analysis_pipeline/make_agree_matrix.py --rank Family       # 先出矩阵（Species 同理）
Rscript reports/virome_analysis_pipeline/make_cls_tool_timing_fig.R               # 图② 左上角耗时图（预渲）
Rscript reports/virome_analysis_pipeline/fig_upset_classification.R               # 图② Species（默认）
# 换阶元时：矩阵与耗时图都要按该阶元的行序重出
python  reports/virome_analysis_pipeline/make_agree_matrix.py --rank Family
Rscript reports/virome_analysis_pipeline/make_cls_tool_timing_fig.R --rank Family \
        --out reports/virome_analysis_pipeline/upset_data/cls_tool_timing_fig_Family.png
Rscript reports/virome_analysis_pipeline/fig_upset_classification.R --rank Family \
        --timing-png reports/virome_analysis_pipeline/upset_data/cls_tool_timing_fig_Family.png
Rscript reports/virome_analysis_pipeline/fig_host_prediction.R                   # 图③（4 面板）
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

**⚠⚠ 2026-10-14 三次微调（用户「左侧的时间图 稍微再往下些」）**：

框**整块下移 60 px** → **`x 100..810 / y 222..585`**（`TIM_T=1-222/1800, TIM_B=1-585/1800`）。
**尺寸不变**（仍 710 × 363 px = 2.3667 × 1.2100 in）→ ⚠⚠ **预渲 PNG 无需重渲**（只改位置、不改尺寸）。

**基准 = 右侧柱面板最高柱顶 `y 222`**（实测 `x1073..1134`，即 6,311 那根）→ 耗时条顶与之齐平。
副标题底 y118 与耗时条顶的间隙由 **44 px → 104 px**（原版几乎贴住副标题，视觉上"顶着"标题行）。

⚠⚠⚠ **上一节「框底压到 525（`4000` 刻度在 527 起，留 2 px）」的约束已作废**：
`4000` 刻度实测只占 **`x 971..1009`**，而本框**右缘 810** → **x 方向本就不重叠**，下移根本不受它约束。
（上节把障碍写成 `x 815..1009` 是**过度宽泛**的包围盒 —— 那是「框右缘到最远刻度」的并集，不是刻度文字自身范围。
⚠ **教训：障碍包围盒必须按「元素自身 bbox」量，不能用「从框边到元素」的跨度**，否则会凭空造出一条假约束。）

**净空判据（必须用 `--no-timing` 底图，不能用成图）**：

| 区域 | 墨迹 |
|---|---|
| 新框全区 `x100..810 / y222..585` | **0** |
| 新占带 `x100..810 / y526..585` | **0** |
| 让出带 `x100..810 / y162..221` | 16,811（移走后腾空） |
| 下方到 set-size 面板上缘 `y1177` | **0**（余 **592 px**） |

**验收（双版均同）**：四边各 4 px 全 **0**（无 `TEXT-OFF-CANVAS`）；框四邻域 上/下/左/右 = **0/0/0/0**；
自检行双版一致：`耗时图框：x 0.036..0.290 y 0.675..0.877 npc｜2.37 × 1.21 in（预渲 PNG）`。

**⚠⚠ 2026-10-14 四次微调（用户「把耗时条形图 再次调整红色区域的大小，使其与预期效果一致」）**：

框**整块放大到用户红框**：`x100..810 / y222..585` → **`x94..946 / y265..1076`**
（**852 × 811 px = 2.8400 × 2.7033 in**，比例 **1.96 → 1.05** —— 由"扁"变"近正方"）。

**红框映射（截图 1136×695 → 画布 2796×1800）** —— 用**柱面板 4 个 y 轴刻度**标定：
```
x_canvas =  46.7 + 2.424 * x_shot
y_canvas = 164.8 + 2.424 * y_shot
```
（`6000` 刻度 shot y53.5 ↔ 画布 294.5；`0` 刻度 shot y389 ↔ 画布 1106.5；
四刻度间距 shot **111.8** / 画布 **271** → **s = 2.424**，逐刻度残差 ≤ 1 px）
红框外缘 shot `x19.5..371.5 / y41.5..376.5` → 画布 **`x94..946 / y265..1076`**。

⚠⚠⚠ **尺寸变了 → 预渲 PNG 必须重渲**（`--w 2.8400 --h 2.7033`），否则硬拉伸零报错（硬规则 22）。
**净空（`--no-timing` 底图实测）**：红框全区 **0 墨迹**；
右缘 946 距柱面板 y 轴刻度（`x971..1009`）**25 px**；底缘 1076 距 set-size 面板（`y1174`）**98 px**。

**验收（双版）**：四边各 4 px 全 **0**；框四邻域 上/下/左/右 = **0/0/0/0**（Species 版；
Family 版下邻的 1,374 墨迹在 `y1174+`，即 **set-size 面板自身**，在框外 98 px）；
自检行双版一致：`耗时图框：x 0.034..0.338 y 0.402..0.853 npc｜2.84 × 2.70 in（预渲 PNG）`。

⚠ **字号未变**（仍 6.0–6.4 pt）—— 框高放大 2.2× 后字**相对**变小，但实测在 2.84 in 宽的插图里
仍清晰可读（印刷约 7.2 × 6.9 cm）。如需再放大，上限受**框宽**约束
（x 轴题 `Total runtime, whole library (s, log)` 在 2.84 in 内最多约 9 pt）。

⚠⚠ **改框尺寸后必须重跑 `make_cls_tool_timing_fig.R`**（预渲 PNG 不随参数更新，会静默嵌旧尺寸图）；
且该脚本的 **`--w`/`--h` 默认值已同步改为 `2.8400`/`2.7033`**，避免下次不带参数跑时尺寸错配。

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

### ⭐⭐⭐ 2026-10-14 00:4x「（点阵矩阵的行标签）放到左侧可以吗」→ 行标签只留最左一处

**需求**：用户截图红框圈的是**矩阵面板那一列重复的集合名**（画布 x826..1030），要求
「行标签移到最左」。实测确认 set-size 面板（x68..237）与矩阵面板**各显示一遍完全相同的集合名**。

**同一条需求踩了 3 次，前两次全部静默失效 / 静默劣化**：

| # | 尝试 | 结果 | 原因 |
|---|---|---|---|
| ① | `matrix = intersection_matrix(...) + theme(axis.text.y = element_blank())` | ✗ **无效**（标签照旧） | 该 theme 先被 `upset()` **L98-102** 的 `+ themes$intersections_matrix` 追加；内置 `[[1]]` 是 **144 元素的完整主题**（显式含 `axis.text.y = element_text`）→ **原样覆盖回去** |
| ② | `themes[["intersections_matrix"]] = theme_upset(8) + theme(axis.text.y = element_blank())` | ✗ **劣化** | `themes[[key]]` 是**整体替换**语义 → 丢掉 `[[2]]` 的内置 blank → 底部冒出**挤成一团的 x 轴集合名 + 竖排 `group`**（零报错） |
| ③ | `themes[["intersections_matrix"]] = Reduce(`+`, upset_themes$intersections_matrix) + theme(axis.text.y = element_blank())` | △ **生效但被否决** | 拼回内置两段后确实关掉了文字，**但**该槽宽度归还给面板 → **点阵/柱整体左移并拉宽**（见下表）→ user 2026-10-14 第二轮：「**点阵整体[平移]保持不对，只把行标签[移到]最左处**」 |
| ✅ | 同上，但写 **`axis.text.y = element_text(colour = "transparent")`** | ✓ 采纳 | **文字仍占位**（grob 宽度照常计算）→ 所有 x 几何与原来**逐像素一致**，只是看不见 |

**内置主题的真实结构**（`ComplexUpset::upset_themes`，务必按此改，别再猜）：

| key | `[[1]]` | `[[2]]` 的 `element_blank` |
|---|---|---|
| `intersections_matrix` | 144 元素完整主题 | `axis.title.y` / **`axis.text.x`** / `axis.ticks.x` ← ⚠ **不含 `axis.text.y`**（所以矩阵行标签默认是显示的） |
| `overall_sizes` | 同上 | `axis.title.y` / **`axis.text.y`** / `axis.ticks.y` ← 本图左侧标签是**我们自己的 `theme_upset(8)` 覆盖回来的** |
| `Intersection size` | 同上 | `axis.title.x` / `axis.text.x` |
| `default` | 同上 | `axis.title.x` / `axis.text.x` |

**⭐ 同源对拍（唯一可靠的判「几何有没有被改」手法）**：同一脚本、只差待验证的一处改动，
各渲一张 `--no-timing` 干净底图，再跑 `_check_upset_geometry.py` 做**逐列差异剖面**：

| 指标 | ②（删文字） | ✅（transparent） |
|---|---|---|
| 交集柱 1 中心 | 999 | **1103**（= 原值） |
| 点阵首列 / 末列 x | 984 / 2434 | **1103 / 2438**（= 原值） |
| 点阵列间距 | 95.7 px | **81..86 px**（= 原值） |
| 与改前的差异列区间 | 全宽（x310..2455） | **只有 x825..1031**（= 标签槽本身），上半区 **0 差异** |

⚠ 柱与点**在两种方案下都逐列对齐**（patchwork 保证）→ 所以**不能靠"柱点是否对齐"判优劣**，
必须看**绝对 x 位置**。
⚠ ③/② 两种坏方案的差别只有「点阵是否平移」，肉眼在成品图上极易误判 →
**判据 = 量像素**。

**新增体检脚本 `_check_upset_geometry.py`**（本图已迭代 5 次，固化成 4 条判据）：

```bash
# 必须喂 --no-timing 干净底图；给参考图则额外打印逐列差异剖面
"C:/Python312/python.exe" reports/virome_analysis_pipeline/_check_upset_geometry.py <fig.png> [<参考图.png>]
```

- **P1** 柱↔点 逐列对齐残差 = 0 px
- **P2** 点阵绝对 x 位置 + 列间距（判「是否被整体平移/拉伸」）
- **P3** 矩阵行标签槽 `x820..1040` 内深色文字段 = **0 处**（判「行标签只剩最左」）
- **P4** 画布四边各 4 px 零墨迹（无 TEXT-OFF-CANVAS）
- ⚠⚠ 脚本会自检「左上是否残留耗时条」：残留时 **P1/P2 不可信**（实测耗时条蓝柱被当成 x541 的"交集柱"）。
- ⚠ 文字的**次像素抗锯齿**会产生暖色边缘 → 用"棕色"掩膜找点阵时会被标签文字误命中，
  故点阵检测限定 `x≥960` 且段宽 20..45 px。
- ⭐ 2026-10-14 01:xx **第 7 次修订**（「红色区域翻转 + 中间空白太多」）把上面的
  `transparent` 改回了 `element_blank()`，并配合 `width_ratio = 0.312` 把点阵**钉回逐像素原位**
  → **本节的"✅ 采纳"结论已被取代，但机制表全部仍然成立**。详见 **规则 28**。

#### ⚠⚠⚠ 顺带查出的真缺陷：`--timing-png` 默认名不随 `--rank` 变化 → 静默嵌错阶元的耗时图

**怎么发现的**：上一条做完后拿 `.bak_20261014_before_rowlabel` 对拍，**上半区**竟出现
`x137..245` / `x659..776` 两处差异（那是耗时条里的工具名与数值标签）→ 顺着查下去。

- **根因**：`fig_upset_classification.R` 原先的默认是
  `file.path("scripts","barbarum","upset_data","cls_tool_timing_fig.png")` —— **恒为 Species 那份**。
  于是 `--rank Family` 只要**没显式传** `--timing-png`，就会把 **Species 行序**的耗时图嵌进 Family 版。
  ⚠ **零报错、零告警**，耗时条照常画出来，肉眼完全看不出。
- **为什么两个 rank 的行序不同**：两图的集合行序都是 **set size 降序**（`sort_sets` 默认 TRUE），
  而两个阶元的 set size 排序在末两位互换：

  | rank | 降序（上→下） | 末两位 set size |
  |---|---|---|
  | Species | … > **VITAP(604)** > **genomad(13)** | 604 > 13 |
  | Family | … > **genomad(654)** > **VITAP(607)** | 654 > 607 |

- **判据（必须量像素）**：耗时条 7 根蓝条（`#0072B2` 精确匹配）的**长度序**：
  - ✅ Family 正确 = `447/1200/169/346/334/`**`874/1027`**
  - ✅ Species 正确 = `447/1200/169/346/334/`**`1027/874`**（末两位互换）
  - ❌ 错嵌 = 与另一个 rank 相同
  ⚠ 只看「耗时条有没有画出来」永远查不出；⚠ 也不能靠"脚本打印的 `ORDER`"自证 ——
  打印的是**计算值**，嵌进去的却是**磁盘上的旧 PNG**，两者必须分别核对。
- ✅ **修法**：默认名按 `--rank` 推导（Species → 无后缀，其余 → `_<Rank>`），
  与 `make_cls_tool_timing_fig.R` 的 `--out` 推导规则**完全一致**（同族硬规则 23）。
  ```r
  tim_png <- get_arg("--timing-png",
    file.path("scripts","barbarum","upset_data",
              if (identical(rank,"Species")) "cls_tool_timing_fig.png"
              else sprintf("cls_tool_timing_fig_%s.png", rank)))
  ```
- **修后验收**：两版与各自 `.bak_20261014_before_rowlabel` 对拍，**上半区 0 差异**，
  下半区只有标签槽 `x825..1031`；桌面 PNG md5 双端一致。

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
C:/Python312/python.exe reports/virome_analysis_pipeline/make_identfilter_data.py   # 汇总（含口径自检）
Rscript reports/virome_analysis_pipeline/fig_ident_filter.R                          # 出图
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
scp <246>:~/ident_filter_upset/ident_filter_upset_* reports/virome_analysis_pipeline/identfilter_data/
# 2) 本地出图（服务器同样可跑，已双端验证）
Rscript reports/virome_analysis_pipeline/fig_ident_filter_upset.R
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
Rscript reports/virome_analysis_pipeline/fig_ident_filter_by_tool.R        # 本机
ssh <246> 'export LD_LIBRARY_PATH=$HOME/mambaforge/lib; cd ~/MMPV-RNA/reports/virome_analysis_pipeline && Rscript fig_ident_filter_by_tool.R'
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
Rscript reports/virome_analysis_pipeline/fig_ident_filter_combined.R     # 默认 --top 16（覆盖 99.39%）
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
Rscript reports/virome_analysis_pipeline/fig_ident_filter_stacked_all.R --set-size stack   # 绝对条长（忠实原生语义）
Rscript reports/virome_analysis_pipeline/fig_ident_filter_stacked_all.R --set-size fill    # 归一化 100%（推荐看图用）
# ⚠ 耗时图是**预渲**的：改了 tool_timing.tsv 需先重出 PNG，再跑主脚本
python reports/virome_analysis_pipeline/make_tool_timing_data.py        # → identfilter_data/tool_timing.tsv
Rscript reports/virome_analysis_pipeline/make_tool_timing_fig.R         # → identfilter_data/tool_timing_fig.png
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

### 本轮新增硬规则（25 条踩坑记录）

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
      Rscript reports/virome_analysis_pipeline/make_tool_timing_fig.R     # → identfilter_data/tool_timing_fig.png
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
- ⚠ **同步清单**（同硬规则 22，五项）：① `fig_ident_filter_stacked_all.R` 的 `TIM_*`
  （⚠ **含 `TIM_SHIFT_PX`** —— 它随「集合条面板左缘」走，改行名槽/名字字号/画布宽都要重核，
  口径见「⑧ 第 15 次修订 → 耗时图 overlay 重新对齐」）；
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

**25. ⚠⚠⚠ Venn 的填充色不能落在「语义色」上 —— `ggVennDiagram` 默认把 fill 映射到区域 count**

**触发**：用户 2026-10-13 22:53「venn图的颜色是怎么回事」。

- ⭐ **现象**：6 工具 Venn 里**只有最大那块**（44,637 = viralverify 独有区）是**纯深绿**，
  其余迅速褪成白 → 看着像"某个集合被高亮/被选中"。
- ⭐⭐ **根因**：`ggVennDiagram()` 把 `fill` 映射到**区域 count**（**连续**标度）。
  原写法 `scale_fill_gradient(low = "white", high = "#009E73", guide = "none")`
  → 渐变高端正好是 **`#009E73`**，而它**就是图例「Retained」色块的同一个 hex**
  （像素实测：Venn 最深区 26,385 px / 图例色块 961 px / 交集柱 14,697 px，**三者均 `#009E73`**）。
- ⚠⚠⚠ **为什么这是错的（语义撞色，最要命）**：Venn 只画 **Retained** 序列
  （代码 `n_retained > 0` 且排除 `viroid_blastn`）→ **每个区域本来就都是 Retained**。
  填充深浅表达的是"计数大小"，**不是**分类归属。而图例就贴在 Venn 右侧 →
  读者必然把这块绿读成"Retained" → **色义被读反**。
- ⚠ **另两个副作用**：
  ① 渐变**暗示排序**（44,637 那块像被"高亮"，其实只是单工具最大值）；
  ② **对比度不足**：数字用 `#1A5C3A` 画在填充上 → 落在 `#009E73` 时只有 **2.35:1**
     （WCAG 正文需 ≥4.5:1），"44,637"/"1,606" 恰好在最难读处。
- ✅ **正解 = 去掉渐变，统一极浅中性填充**：
  ```r
  scale_fill_gradient(low = "#F2F2F2", high = "#F2F2F2", guide = "none")
  ```
  `#F2F2F2` 明显浅于图例的 `#BFBFBF` / `#7A7A7A` 两个灰 → 不会被误读为分类色；
  区域只靠 `#333333` 椭圆描边 + 计数标签区分。
  改后数字 `#1A5C3A` on `#F2F2F2` 对比度 = **7.38:1**（达标）→ **标签色不动**
  （保留它与柱顶 `%` / 左侧保留率的同色语义）。
- ⚠⚠ **实现坑**：`fill` 是**连续**标度（`count`）→ **不能**换 `scale_fill_manual`
  （报 `Continuous value supplied to discrete scale`）；
  必须用**两端同色**的 `scale_fill_gradient` 才能得到"平填充"。
- **验收（双端）**：Venn 区 `#009E73` 像素 **26,385 → 0**（本机）/ **0**（服务器）；
  `#F2F2F2` = 470,471 px（本机）/ 458,240 px（服务器）；
  与上一版对拍差异 bbox `x2006..2774 / y447..901` —— **完全落在 Venn 框内**
  （`x1996..2856 / y406..1294`）→ **除 Venn 外零回归**。
- ⭐⭐ **通用判据（本条最值钱）**：**给任何图层加"数值渐变"前，先确认渐变两端不落在
  本图已用的「语义色」上。** 本图语义色表 = `STATUS_COLORS`
  （Retained `#009E73` / Removed by UniProt `#BFBFBF` / Removed by CDD `#7A7A7A` /
  Viroid track `#E69F00`）→ **Venn 填充不许碰这四色，连灰都不行**。
  ⭐ 推论：**"同一张图里，颜色要么全是语义的，要么全是数值的，不能混。"**
  混用必然产生"读者按语义读、作者按数值画"的错位。

**26. ⚠⚠⚠ 想关掉 ComplexUpset 某面板的轴元素，只能改 `themes[[key]]` —— 但 `themes[[key]]` 是「整体替换」，且「关掉文字」会静默改变整块几何**

**触发**：用户 2026-10-14 00:4x「（点阵矩阵的行标签）放到左侧可以吗」→ 二轮纠正
「**点阵整体[平移]保持不对，只把行标签[移到]最左处**」。

- ⭐⭐ **先弄清内置主题到底 blank 了哪些**（`Rscript -e 'str(ComplexUpset::upset_themes)'` 逐项列 `element_blank`）：

  | key | `[[2]]` blank 掉的 | 含义 |
  |---|---|---|
  | `intersections_matrix` | `axis.title.y` / `axis.text.x` / `axis.ticks.x` | ⚠ **不含 `axis.text.y`** → 矩阵行标签**默认就显示** |
  | `overall_sizes` | `axis.title.y` / `axis.text.y` / `axis.ticks.y` | set-size 面板的集合名**默认被关**；本图是我们自己 `theme_upset(8)` 覆盖回来的 |
  | `Intersection size` / `default` | `axis.title.x` / `axis.text.x` | |

- ⚠⚠⚠ **两条静默失效路（都零报错）**：
  ① 写在**图层链**（`matrix = intersection_matrix(...) + theme(axis.text.y = element_blank())`）→ **无效**：
     该 theme 先被 `upset()` **L98-102** 的 `+ themes$intersections_matrix` 追加，
     而内置 `[[1]]` 是 **144 元素的完整主题**（显式含 `axis.text.y = element_text`）→ **原样覆盖回去**。
  ② `themes[["intersections_matrix"]] = theme_upset(8) + theme(...)` → **整体替换**语义 →
     丢掉 `[[2]]` 的内置 blank → 底部冒出**挤成一团的 x 轴集合名 + 竖排 `group`**。
- ✅ **正解**：`themes[["intersections_matrix"]] = Reduce(`+`, ComplexUpset::upset_themes$intersections_matrix) + theme(...)`
  —— **先拼回内置两段**（`upset_modify_themes` 只做命名 list 合并，值是单个 theme 完全合法），再叠加自己的改动。
- ⚠⚠⚠ **本条最值钱的部分：`element_blank()` 会把该槽宽度归还给面板 → 整块几何被改**。
  实测同源对拍（画布 2796×1800，`--no-timing`）：

  | 指标 | 删文字（`element_blank`） | ✅ `element_text(colour = "transparent")` |
  |---|---|---|
  | 交集柱 1 中心 | 999 | **1103**（= 原值） |
  | 点阵首列 / 末列 x | 984 / 2434 | **1103 / 2438**（= 原值） |
  | 点阵列间距 | 95.7 px | **81..86 px**（= 原值） |
  | 与改前的差异列区间 | 全宽 `x310..2455` | **只有 `x825..1031`**（标签槽本身），上半区 **0 差异** |

  ⚠ 柱与点**两种方案下都逐列对齐**（patchwork 保证）→ **不能靠"柱点是否对齐"判优劣**，
  必须看**绝对 x 位置**；⚠ 用 `"transparent"` 而非 `"white"`（白字在非白底会现形）。
- ⭐ **可复用判据**：任何「关掉某轴文字」的改动，都要跑**同源对拍**（同一脚本、只差这一处改动，
  各渲一张 `--no-timing` 底图）看**逐列差异剖面** —— 差异只落在**被改的那一槽**才算零回归。
  工具 = `reports/virome_analysis_pipeline/_check_upset_geometry.py`（P1 柱↔点对齐 / P2 点阵绝对位置 /
  **P3 标签槽内文字段 = 0** / P4 四边 4 px 零墨迹；带参考图时自动出差异剖面）。
  ⚠ 脚本会自检「左上是否残留耗时条」，残留时 P1/P2 不可信；⚠ 文字**次像素抗锯齿**的暖色边缘
  会污染"棕点"掩膜 → 点阵检测须限定 `x≥960` 且段宽 20..45 px。
  ⚠⚠⚠ **2026-10-14 01:xx 第 7 次修订推翻了本条的"✅ 结论"（机制部分仍然成立）**：
  为满足「吃掉中间空白」的新需求，最终**改回 `element_blank()`**，并**同时把
  `width_ratio` 从 0.26 提到 0.312** 来补偿"面板变宽 + 左移" → 矩阵面板**回到
  `1054..2488` 逐像素原位**，点阵首/末列 `Species 1103/2438`、`Family 1100/2440`
  **与 transparent 版完全一致**，而 set-size 面板净加宽 147 px。→ 详见**规则 28**。

**27. ⚠⚠⚠ 预渲脚本的「消费者侧」默认路径也必须随区分参数变化 —— 否则静默嵌错阶元的资产（规则 23 的另一半）**

**触发**：2026-10-14 做行标签修订时，用 `.bak` 对拍发现**上半区**（耗时条）出现差异 → 顺藤摸瓜。

- 规则 23 管的是**生产者**（`make_cls_tool_timing_fig.R --out` 要随 `--rank` 变）；
  本条管**消费者**：`fig_upset_classification.R` 的 `--timing-png` **默认值原先恒为**
  `cls_tool_timing_fig.png`（Species 那份）→ `--rank Family` 只要不显式传参，
  就把 **Species 行序**的耗时图嵌进 Family 版。**零报错、零告警、图照常出**。
- ⚠ **两个 rank 的行序为什么不同**：行序 = **set size 降序**，而末两位在两个阶元互换 ——
  Species `VITAP(604) > genomad(13)`；Family `genomad(654) > VITAP(607)`。
- ⭐ **判据（只能量像素）**：耗时条 7 根 `#0072B2` 蓝条的**长度序** ——
  Family 正确 `…/`**`874/1027`**、Species 正确 `…/`**`1027/874`**（末两位互换）。
  ⚠ **不能靠脚本打印的 `ORDER` 自证**：打印的是**计算值**，嵌进去的是**磁盘上的 PNG**，
  两者必须分别核对（这次就是"打印对、嵌错图"）。
- ✅ **修法**：消费者侧默认路径按 `--rank` 推导，与生产者侧规则**逐字一致**：
  ```r
  if (identical(rank, "Species")) "cls_tool_timing_fig.png"
  else sprintf("cls_tool_timing_fig_%s.png", rank)
  ```
- ⭐ **通用教训**：凡「预渲资产 + 参数化宿主」的组合，**生产者与消费者的命名规则必须成对**。
  只改一边 = 静默串味。排查手法 = **同源对拍 + 逐列差异剖面**，差异落在**非目标区域**就是串味。

**28. ⚠⚠⚠ 「两列之间的空白」压不掉，除非同时动「标签槽」和 `width_ratio` —— 单调任一个都无效**

**触发**：用户 2026-10-14 01:xx 截图（红框圈 set-size 柱条区）「**红色区域翻转呢，另外中间空白太多了不**」。
AskUserQuestion 答复 = 「**柱条左右镜像**」+「**加宽 set-size 面板填满它**」。

**(a) 柱条镜像 = 绕开 ComplexUpset 默认的 `scale_y_reverse()`**

- ⚠⚠ 机制（`upset()` **L217-230**）：对**左侧** set-size 面板，`default_scale = scale_y_reverse()`
  → 再 `coord_flip()` → 横轴反向 → 柱从**右缘往左长**。想镜像，必须**绕开**它。
- ⭐⭐ 关键在 **L228**：`scale_if_missing(set_sizes, axis = "y", scale = default_scale)`，
  而 `scale_if_missing()` = 「`get_scale(annotation, "y")` 非空就返回 `list()`」
  → **只要自己挂一个 y 标度，默认的 `scale_y_reverse()` 就完全不加**。正解：
  ```r
  set_sizes = upset_set_size() +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    xlab("Set size") + theme_upset(8) + theme(axis.title.y = element_blank()),
  ```
  ⚠ `expansion(mult = c(0, 0.05))`：**非 reverse 标度下按常规方向**（0 端不留白、高端留 5 %）。
  写反 → 柱起点离开刻度 0（= 规则 10「`scale_y_reverse` 下 expansion 方向反向」的反面场景）。
- **验收（量像素）**：改前柱条**右端**对齐（Family 7 根右缘全在 734）；
  改后**左端**对齐 `x≈254`、右端各不相同（`769/720/718/643/516/297/294`）。

**(b) 「中间空白」由两块拼成，`width_ratio` 单独调**一点也减不掉

- ⚠⚠ 先量清楚 `width_ratio` 是什么：`upset()` **L261-264** `plot_layout(widths = c(w, 1-w))`
  → 是**两列「面板」宽之比**（不含坐标轴区）。实测 `w=0.26` → 面板 542 : 1543 px = 0.26 : 0.74 ✓。
- ⭐⭐ **空白的两个来源**（画布 2796×1800，逐行扫非白像素段实测）：

  | 来源 | 宽度 | 能否压掉 |
  |---|---:|---|
  | ① **列 2 的左侧坐标轴区** —— 顶部交集柱面板的 y 轴刻度文字（`5,000`…`0`）住在那儿 | ~150 px | ❌ **结构性**：那是柱图的 y 轴，必须留 |
  | ② 矩阵那串**行标签槽** —— 只有 `element_blank()` 才归零（`transparent` 仍占位） | ~205 px（transparent 时） | ✅ 可压 |

- ⭐⭐⭐ **四组同源对拍**（Family，`--no-timing`；量 set-size / 矩阵两面板的灰条带左右缘）：

  | 组合 | set-size 面板 | 矩阵面板 | 中间空白 |
  |---|---|---|---:|
  | `transparent` + wr 0.26（改前交付版） | 254..757 | 1054..2488 | **297 px** |
  | `transparent` + wr 0.45 | 254..1192 | 1342..2488 | **297 px（不变！）** |
  | `element_blank` + wr 0.26 | 254..796 | 945..2488 | **149 px**（但点阵左移 109 px） |
  | ✅ **`element_blank` + wr 0.312** | **254..904** | **1054..2488（原位）** | **150 px** |

  ⚠⚠ **`width_ratio` 变大只会把"柱条 + 矩阵"整体右推，空白宽度纹丝不动** —— 这是本条最反直觉的点。
- ✅ **正解 = 两件事一起做**：`element_blank()` 消掉 205 px 的透明标签槽；
  再把 `width_ratio` 从 0.26 提到 **0.312**，把矩阵面板**推回原位**（补偿 blank 造成的"变宽 + 左移"）
  → 最终 **点阵逐像素不动、set-size 面板净加宽 147 px（+29 %）、空白 297 → 150 px（−50 %）**。
- ⚠ **试过且无效**：把 set-size 面板右 margin / 矩阵左 margin 归零（`margin(5,0,5,5)`）
  → 只省 **9 px**（150 → 141）。因为列 2 的左轴区由**顶部柱面板**的 y 轴决定，矩阵自己的 margin 不是瓶颈。
- ⭐ **验收（同源对拍，唯一可靠）**：新版 `--no-timing` vs 改前交付版，逐像素差 `|Δ|>25`：

  | 区域 | 差异像素 |
  |---|---:|
  | 矩阵带 `y1176..1600` | **59591，且差异列区间只有 `x254..873`**（= set-size 面板自身） |
  | 矩阵面板本体 `x≥874` | **0** |

  → **点阵完全没动**；`_check_upset_geometry.py` 的 P2 首/末列两 rank 均与改前一致。
- ⚠⚠ **改布局后必须复核的两处**：① `_check_upset_geometry.py` 的 `X_SETSIZE_END` / `X_TAGSLOT`
  常量（set-size 面板右缘 757 → 904，否则 P3 会把灰柱当文字误报）；
  ② **左上耗时图框净空**（框 `x94..946 / y265..1076`）→ 已固化成 **P5 判据**（框内墨迹 = 0）。
  实测框右缘 946 距柱图刻度文字左缘 **971（25 px）**，**与改前相同**（刻度文字锚在柱面板上，不受 wr 影响）。
- ⚠⚠ **P1 的假阳性**：顶部柱图 y 轴刻度文字（`x≈998..1028`）经 **ClearType 子像素抗锯齿**会产出
  `(204,139,95)` / `(95,139,204)` 这类**彩边** → `max-min>60` 的"饱和色"判据会把它当成
  一根**假交集柱**（实测凭空多出 `x≈1010` 那根）。→ 柱检测下界必须 ≥1050，并把
  「列内饱和像素数」阈值从 20 降到 3（否则最矮的几根柱 `cs<20` 被漏掉，P1 永远"数量不等"）。


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
cd reports/virome_analysis_pipeline
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
