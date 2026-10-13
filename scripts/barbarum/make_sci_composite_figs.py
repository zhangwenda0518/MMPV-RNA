#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_sci_composite_figs.py
==========================
宁夏枸杞 (Lycium barbarum) 已知病毒模块 —— SCI 正文组合主图

把多个分析模块的结果合成一张正文主图 (7 个面板 / A-G):

    A  检测 + 组装   Detection & genome-recovery funnel (9 species, 检出=高置信/四段堆叠条)
    B  群体内变异   Intraspecific variation landscape (GCVA, OR489165.1)
    C  共感染结构   Co-infection structure
    D  序列同一性   Pairwise NT identity distribution (只出分布图, 与 A 同行)
    E  变异·选择    Self-drawn: dN-vs-dS scatter + per-gene dN/dS (通栏, 双子图)
    F  变异·准种    Self-drawn: allele-frequency spectrum + substitution spectrum (通栏, 双子图)
    G  组织嗜性     Tissue tropism (virus x organ)

**E/F 是什么内容**（`--ef` 切，**rev.h 起默认 `self`**）：
  `self`（**默认**）= **自绘** GCVA 变异分析，数据来自两张可解析的表（见下）：
      E = `dN vs dS` 双对数散点(＋中性线 `dN=dS`) │ 逐基因 `dN/dS`(抖动散点＋中位＋IQR＋中性线 1.0)；
      F = 等位频率谱 AFS(低频准种 vs 固定) │ 12 类替换谱(Ts/Tv，**剔掉 self 柱**)。
      与 A–D/G 共用同一套 `plt.rcParams` + Okabe-Ito + "只留左/下 spine、无网格" 语言，
      **是唯一可用于投稿的形态**。
  `variation` = 嵌入 `04_post_analysis/<GCVA>/` 里已有的变异分析成品图
      （E 默认 `boot_dnds`，F 默认 `afs`；`--panel-e`/`--panel-f` 可换任意两张）。
      **外观与主图不统一，仅作对照留存，不推荐投稿** —— 理由见下一节。
  `deg` = rev.d 的行为：自绘宿主 DEG 火山图 + GO/KEGG 富集。
      宿主内容不在主图里之后，可加 `--host-figure` 另存一张。
  `none` = **rev.j：主图不要 E/F 两行**（只留 A / D / B / C / G 五个面板，≈1.0 页）。
      变异分析内容全部交给另出的 `--var-rich` 组合图（`_VariationRich`，9 面板），
      两图之间**无面板重复**。落地原因：rev.h 主图 20.90 in ≈ 1.4 页，
      而通栏面板的可用宽度（9.775 in）决定了双栏子图的最小高度，几何上无解，
      只能靠"把变异内容移出主图"来压高度。

**rev.k（2026-10-12）**：`--ef none --b-last` —— 把「群体基因组景观」面板挪到 **C/G 下方**，
并把面板字母按**读序**重排为 **A–E**（A 检出 → B 同一性 → C 共感染 → D 组织嗜性 → E 群体景观）。
同一次改动还修了该面板 π 曲线的**口径 bug**：旧版画的是 `PI` 列（窗口内各位点 π 的**求和**）
却把轴标成 "π (×10⁻³)"，实际 = 真 per-site π × 0.590（**偏低 41 %**）；现 ÷ 窗口长度，
与 `_VariationRich` 面板 B 同口径。⚠ 重排后正文里的「面板 X」引用必须同步改。
不加 `--b-last` 则完全保持 rev.j 的 A→D→B→C→G 旧读序（字母不动）。

**`--append-figs`（2026-10-12，用户指定）**：在主图**下方追加拼图行**，字母按读序顺延。
`;` 分隔行、`,` 分隔同一行内的面板。不带值 = 用户 2026-10-12 定的两行（见 `CHOSEN_ROWS`）：
`mafsummary|circos`（两张按原图嵌入的 PPT 图）+ `E|L`（两个自绘面板）。
（rev.s：**位图行在上、自绘行在下**；字母按读序顺延 → 新 E|F = 位图行、新 G|H = 自绘行。）
**位图按原图嵌入、不重绘** —— 都是整页独立图，按列宽嵌入后印刷字号只有 4.5–5.2 pt
（< 6 pt 可读线），用户明确要求就用这几张，故照办并把字号如实打进日志 / README §10.12。

**rev.h 为什么放弃"嵌入"（2026-10-10 定，本轮的关键转折）**：
  04 的成品图由 snpgenie / vcf_viz **各自的脚本**产出 —— **白底灰网格、红点蓝星号、
  紫色填充、各自的默认字体字号**；而主图 A–D/G 是统一后的"只留左/下 spine、无网格、
  Okabe-Ito、标题 7.3 pt `loc='left'`" 语言。两者拼在同一张图上，**一眼就是"拼凑"**，
  过不了期刊对同一张图子面板风格一致的要求。
  **"嵌入"这条路本身不成立，与挑哪两张无关**（rev.g 挑得再准也救不了）——
  所以 rev.h 改为**自绘**：数据其实都在，不是只有截图 ——
  `snpgenie/**/product_results.txt`（15 列数值）与
  `vcf_viz/Table1_Unified_HighQuality_Variants.csv`（44,357 行）都是可解析的表。

**自绘 E/F 的数据来源**（`--dnds-table` / `--var-table` 可覆盖默认路径）：
  * E ← `<out>/_data/gcva_dnds.tsv`，由 `scripts/barbarum/extract_snpgenie_dnds.py`
    汇总 **161** 个 GCVA 样本的 `virus-SNPGenie/**/product_results.txt`
    （**966 行 = 161 样本 × 6 基因**）。
    product→基因 的映射取自**权威注释** `virus-annotations/OR489165.1.gtf`
    （`N/P/P4/M/G/L` = product `1..6`；CDS 长度 1515/1548/720/861/627/6783 nt
    与 `N_sites + S_sites` 反推的 1512/1545/717/858/624/6780 互相印证）。
    注意 GTF 坐标与本脚本 `GCVA_GENES` 的坐标是**两套注释**，但基因顺序一致，
    故配色 `OKABE[i]` 可直接复用。
  * F ← `04_post_analysis/<GCVA>/vcf_viz/Table1_Unified_HighQuality_Variants.csv`
    （**144 样本 × 44,357 位点**；`REF`/`ALT` 直接给出替换类型，**Ts/Tv = 2.59**）。
  * **口径差异必须写进图注**：E 用 **161** 个样本（要求位点在 CDS 内、SNPGenie 能定 dN/dS），
    F 用 **144** 个样本（不要求位点在 CDS 内）。**两者不是同一批样本、不可混算** ——
    实测交集只有 **129** 个（32 个只在 snpgenie 侧、15 个只在 vcf_viz 侧）。
  * ⚠ **务必取 `known_virus_barbarum/` 下的表**。`known_virus_all_v2/` 里有一张**同名**表，
    是 **7 月旧运行**（175 样本 / 55,768 位点 / Ts/Tv 2.64）；rev.h 期间我曾误用它做统计，
    数字全部对不上。本图其余面板（A/B/C/D/G）都读 `known_virus_barbarum`，F 必须一致。

**dN/dS 为什么画 `dN vs dS` 散点而不是比值的箱线（`self` 的核心取舍）**：
  dN/dS 的分母 dS 可以 ≈ 0（实测 6 基因里 `dS = 0` 的"样本×基因"对共 **49** 个：
  M 11 / G 28 / L 7 / N 1 / P4 2 / P 0）。一旦按比值画箱线，
  **单个 dS≈0 的伪影就能把 y 轴拉到 29** —— 这正是 04 `Fig05` 翻车的原因，其余基因被压成一条线。
  **`dN vs dS` 双对数散点 + `dN=dS` 中性线**天然免疫这个坑，并直接读出关键结论：
  **917 个有效"样本×基因"对里只有 1.5 %（14 个）落在中性线之上 → 净化选择压倒性。**
  右图的逐基因比值**只对 `dS > 0` 的子集计算**，并把因 `dS = 0` 被剔除的对数打进日志。

**自绘 E 的实测结果**（161 样本，**排除 dS = 0 后**；可直接写进正文/图注）：
    gene   n    dS 中位    dN 中位    dN/dS 中位   Q3      max      dS=0
    N     160   0.02672    0.00255     0.096     0.109   0.286      1
    P     161   0.02070    0.00331     0.186     0.239   0.465      0
    P4    159   0.02472    0.00180     0.073     0.097   0.294      2
    M     150   0.02584    0.00000     0.000     0.058   0.291     11
    G     133   0.00772    0.00406     0.442     0.532   1.866     28
    L     154   0.02203    0.00227     0.092     0.106   0.567      7
  **G 基因：dS 只有他者的 1/3（0.0077）× dN 最高（0.0041）× 比值最松（0.442 = P 的 2.4 倍）**
  —— 与"G 受的约束最松"这一叙事一致。
  **排除 `dS = 0` 后全体最大值仅 1.866**，彻底证实 04 `Fig05` 里那个 **29 是伪影**，不是生物学信号。

**`variation` 模式的选图依据（保留作对照，rev.g 实测，勿删）**：
  把 45 张候选**全部下载逐张看过**之后，发现"文件名看着合适"与"实际能不能用"是两回事：
  * `dnds_gene`(Fig05) 的 y 轴被**单个 dS≈0 的伪影（29.1）拉到 31** → **不能进正文**；
  * `boot_dnds`(Fig06) 用 10,000× bootstrap 的 CI + 中性线 1.0 的框架，免疫上述坑；
  * `mut_spectrum`(Fig07) 的 x 轴前 4 根是 **`A>A/T>T/C>C/G>G`（与参考一致的位点，不是替换）**，
    把 y 轴顶到 8000，真正的 Ts/Tv 被压到底部 → **需重绘才可用**（自绘 F 右图就是重绘版）；
  * `afs`(Figure5) 直接把"低频准种 vs AF=1.0 固定"的双峰画出来，**无干扰项、无伪影**。

  判定要**同时**满足：「内容干净」+「通栏高 ≲ 6 in」+「印刷字号 ≥ 7 pt」。
  下表是**裁白边 + 切源图标题之后的实测值**（由出图日志抄回，不是按源图估的）：

    key              裁后 px      宽高比  通栏高   印刷字号  内容判定
    boot_dnds        2617x1337    1.96    4.99 in  8.3 pt   ✓ 干净 + 逐基因 + 有统计
    afs              2910x1632    1.78    5.48 in  7.5 pt   ✓ 双峰清楚
    dnds_gene        2732x1538    1.78    5.50 in  7.9 pt   △ 尺寸好但 y 轴被伪影压扁
    mut_spectrum     2914x1637    1.78    5.49 in  7.4 pt   △ 尺寸好但有 self 柱干扰
    hotspots         2917x1645    1.77    5.51 in  7.4 pt   ○ 可用（面板 B 的细节补充）
    pca_2d           2617x1941    1.35    7.25 in  8.3 pt   ○ 偏高但可读
    intrahost_div    4703x1634    2.88    3.40 in  **4.6 pt** ✗ **源图是 16 in/4800 px 画布**
    interhost_dnds   (未裁)                10.04 in           △ 内容最好但**放不进正文**
    intrahost_pi     (未裁)                10.04 in           △ 同上（与上图是对偶）

  `intrahost_div` 是个**反直觉的坑**：通栏高只有 3.40 in（最省版面），但源图是
  4800 px / 16 in 画布，缩到通栏只剩 **4.6 pt** —— **版面最省的那张恰恰是最读不了的**。
  凡 `EF_LIBRARY` 里 4800 px 那批（`intrahost_div` / `manhattan` / `var_density` /
  `lollipop` / `var_landscape` / `var_landscape_b`）都有这个问题，**选它当面板前先看日志字号**。
  即：**内容最好的两张（Fig03/Fig04 的 dN-vs-dS 散点）尺寸不合格，
  尺寸最好的两张（Fig05/Fig07）内容有缺陷** —— 所以 `variation` 模式下正文只能落在
  Fig06 + Figure5。**要同时拿到"干净内容 + 合格尺寸 + 够大字号"且风格统一，只能自绘。**
  全部 45 张的处置见 `reports/sci_composite/README_SCI_composite.md` §九。

**为什么正文最多只能放 2 张变异图 + `--var-figure` 的定位**（2026-10-10 定）：
  `variation` 模式的源画布 10 in，**只能通栏**（半栏降到 ≈3.5 pt）；通栏高 ≈5.5 in，
  而 A4 正文区约 9.5 in 高 → **正文上限 ≈1.6 张**。所以 E/F 各占一张，
  其余有意义的收进 `--var-figure` 单独出的 `_Variation`（补充材料）。
  **注**：`_Variation` 目前仍是嵌图（同一批成品图），按 rev.h 的标准它同样"风格不统一"，
  待定夺（选项：改自绘 / 转补充材料并接受外观 / 废弃）。

版面（rev.d 起，2026-10-10 定）
----------------------------
    row0 = A 漏斗堆叠条 (左, fig 0.155-0.655) | D 同一性分布 (右窄栏, 0.710-0.945)
    row1 = B 群体内变异景观（通栏）
    row2 = C 共感染结构 | G 组织嗜性
    row3 = E 变异分析 2（通栏）
    row4 = F 变异分析 3（通栏）
面板字母沿用 v3、不因位置调整而重排（避免正文里的「面板 X」引用错位）；
D 挪到 A 右侧后，读图顺序是 A → D → B → C → G → E → F，改字母需与正文同步。

面板 A 的**两条安全阀**（都是实测文字盒驱动，硬编码值踩过坑）
-----------------------------------------------------------
1. 标注列溢出 → 顶 xlim（原就有）；
2. 图例与 x 刻度标签叠字 → 图例逐步下压 + `xlabel.labelpad` 自适应撑到图例下方。
   锚点硬编码 `-0.028` 时，`ncol`/段数一变就与刻度的 "50 / 150" 叠在一起（rev.d 实拍发现）。
   `panel_A` 会把两条安全阀的实测结果打进日志（含 y 刻度标签最宽值与左端 x），
   版式左边距只按这个实测值定，不再用估算。

面板 A 的口径（rev.d：检出 = filter 阶段的高置信结果）
---------------------------------------------------
* 横向堆叠条，**条长 = 该物种的 raw 命中 run 数**；由左至右四段
  「完整组装 (N% <= 5) → 组装未达标 (N% > 5) → 高置信检出未送组装 → 仅 raw 命中(未过高置信)」，
  主线三级 raw → 检出 → 完整组装 的数值与率写在条**右侧的标注列**里
  （率的分母统一 n = 445）。
* **「检出」用 `02_filtering/high_conf.summary.tsv`**（不是 detect 阶段的 `best`）：
  依据 `virome_analysis_pipeline/auto_known_virus.py::get_summary()` 的约定 ——
  「filter 运行后自动切换到 high_conf」，其后的 variants / assembly / dvg / report 全用它。
  于是行集 = 至少通过高置信过滤的物种（barbarum 实测 9 种），
  漏斗 726 → 349 → 349 → 303（**严格单调**：所有高置信检出都送进了组装，46 条卡 N% > 5）。
  注意 detect 阶段自己的图（Fig1/Fig2）仍按 `best` 画（262 run / 58.9%），二者口径不同。
* 「组装失败」与「根本没送组装」分成两段，是合并原 D 面板之后才需要的区分；
  本口径下「高置信检出未送组装」恒为 0（全 0 段自动不进图例）。
* 三层源表：raw `01_detection/summary/all_viruses.raw.tsv` /
  检出 `02_filtering/high_conf.summary.tsv` /
  组装 `06_extraction/extraction_status.tsv`（`status == extracted` 即 N% <= 5）；
  组装层先按 run 与检出层取交集（`_run_of()` 归一 `CRR..._OR...` 这类样本名）。
* 数字放在条**外**，所以 6.3 pt 就够读 —— 这是它比气泡点阵更适合密集多行的原因。

v4 气泡点阵（已实测弃用，勿再尝试）
-----------------------------------
曾把面板 A 换成 `virome_analysis_pipeline/stage_funnel_summary.py:plot_matrix()` 的
「物种 × 五阶段留存气泡点阵」（它在 `.tmp_demo_funnel/preview.html` 的「Stage Funnel
图型对比」里被标注为备选1 / 已选定为正式版）。
问题出在几何：官方 figsize 是 `(8.2, max(3.2, 0.62*nr+1.6))`，nr=26 时高达 17.7 in ——
那是给整页 / 滚动浏览设计的。塞进组合图的一栏后行距只剩 18.3 pt，
泡内字号 6.75 pt；当时临时加的 `--emit-variants` 诊断（本次回退已删）实测
26 行 = 5.8 pt ✗ / 12 行 = 6.5 pt ✓ / 9 行 = 6.9 pt ✓ —— 想达标就只能砍行数。
结论：**分类轴上的气泡不适合本图的行密度**。
（素材保留：`.tmp_demo_funnel/`、`virome_analysis_pipeline/stage_funnel_summary.py`、
`_diag/PanelA_bubble_rowcounts.png`。）

设计约束
--------
* 全部面板由**底层数据表重绘**（不用已渲染的 PNG），保证字体/配色/线宽统一；
* Okabe-Ito 色盲友好配色；Arial -> Liberation Sans -> DejaVu Sans 字体回退；
* 300 dpi PNG + 矢量可编辑 PDF (`pdf.fonttype=42`)；
* 分母口径固定为 **445**（`02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly` run 去重）。
* 保存必须 `bbox_inches="tight"`（画布自适应，避免轴标签被裁；同时它是
  「有元素画到画布外」的探测器 —— 画布异常变大就说明有元素越界）。

用法
----
    export LD_LIBRARY_PATH=/home/zhangwenda/mambaforge/lib:$LD_LIBRARY_PATH
    python make_sci_composite_figs.py \
        --root /home/zhangwenda/MMPV-paper/goji-virome \
        --out  <输出目录>
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.gridspec import GridSpec
from matplotlib.textpath import TextPath
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------
# 全局绘图风格（与 scripts_20260825/figures/*.py 的房屋风格一致）
# --------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.edgecolor": "#333333",
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 2.2,
    "ytick.major.size": 2.2,
})

# Okabe-Ito 色盲友好离散色 —— 与 virome_analysis_pipeline/virus_metadata_plot.py 同成员
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000", "#999999"]
# 语义命名常量: 显式赋值, 不依赖 OKABE 的排列顺序 (重排 OKABE 时不会静默错配)
ORANGE, SKY, GREEN, YELLOW = "#E69F00", "#56B4E9", "#009E73", "#F0E442"
BLUE, VERM, PURPLE = "#0072B2", "#D55E00", "#CC79A7"
GREY = "#999999"      # 中性灰 = `virome_phylo_pipeline/utils/sci_palette.py` 的 `C_NEUTRAL`
# ↑ rev.t：旧值 `#8F867B` 是**偏暖的棕灰**（R>G>B），不在 Okabe-Ito / 项目调色板里，
#   与全图其它中性灰并列时明显"发黄"。统一到调色板的 `#999999`。
# 位置色序沿用本脚本既有约定 (BLUE 打头), 保持既有图件配色不位移
OKABE = [BLUE, ORANGE, VERM, GREEN, PURPLE, SKY, YELLOW, GREY]
DN, UP = BLUE, VERM          # 下调 / 上调
# 连续色标: 顺序型统一 viridis 族 (浅(低) -> 深(高), 与原 YlOrRd 方向一致)
CMAP_SEQ = "viridis"
CMAP_SEQ_R = "viridis_r"

DENOM_DEFAULT = 445          # 权威分母 (barbarum run 去重)
A_PITCH_PT = 14.0            # 面板 A 行距 (pt)：9 行 -> 轴高 (9+2.0)*14/72 = 2.14 in
                             # （旧口径 26 行时是 11.0；行数降到 9 后放宽，两行标签也放得下）


# --------------------------------------------------------------------------
# 小工具
# --------------------------------------------------------------------------
def panel_letter(fig, ax, letter, dx=-0.062, dy=1.045, size=12, anchor_x=None):
    """在轴框左上角外挂面板字母。

    `dx` 是**轴分数**（默认 `-0.062` = 轴框左外 6.2 % 轴宽）。
    `anchor_x`（**figure 分数**）给定时改用**绝对 x** 定位 —— 供「同一列里轴框被
    `_shrink_x` 内缩过的自绘面板」与「轴框 = 整格的位图面板」**字母对齐**用
    （见 `panel_E_land` / `panel_L` 的 `pack_cell`）。

    ⚠ 为什么需要它：`panel_var`（位图行）的轴框就是整格 → 字母落在
    `cell_left − 0.062×cell_w`；而自绘面板的轴框被内缩了 8.5–9.5 % 列宽 →
    沿用 `dx=-0.062` 会把字母**整体右移 0.34–0.48 in**，与上一行的字母错位。
    """
    if anchor_x is not None:
        _p = ax.get_position()
        if _p.width > 1e-9:
            dx = (anchor_x - _p.x0) / _p.width
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="top", ha="left", color="#333333")


SHORT_NAME = {
    "Grapevine wood holobiome associated mycobunyavirales-like virus 2":
        "Grapevine wood holobiome assoc. mycobunya-like virus 2",
    "Camellia japonica associated betaflexivirus 1": "Camellia japonica betaflexivirus 1",
    "Camellia ringspot associated virus 1": "Camellia ringspot assoc. virus 1",
    "Cucumber mosaic virus satellite RNA": "Cucumber mosaic virus satRNA",
    "Grapevine-associated RNA virus 4": "Grapevine-assoc. RNA virus 4",
    "Cytorhabdovirus sp. 'lycii'": "Cytorhabdovirus sp. 'lycii' (GCVA)",
    "Cytorhabdovirus_sp._lycii_OR489165.1": "Cytorhabdovirus sp. 'lycii' (GCVA)",
}


def pretty_virus(name: str, maxlen: int = 56) -> str:
    """把长物种名压成图上可读的短名。"""
    n = str(name)
    if n in SHORT_NAME:
        return SHORT_NAME[n]
    n = n.replace("sp_lycii", "sp. 'lycii'").replace("_", " ")
    n = n.replace(" associated ", " assoc. ").replace("-associated", "-assoc.")
    n = n.replace("mycobunyavirales-like", "mycobunya-like")
    for acc in (".1", ".2", ".3"):
        n = n.replace(acc, "")
    n = n.replace("NC ", "NC_").strip()
    return n if len(n) <= maxlen else n[: maxlen - 1] + "\u2026"


def read_table(path, **kw):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return pd.read_csv(path, **kw)


def shorten(term, n=42):
    t = str(term)
    return t if len(t) <= n else t[: n - 1] + "\u2026"


def wrap_name(s, width=30, maxlines=2, merge_short=8):
    """物种名**按词折行、不截断**：最多 maxlines 行，最后一行吸收余下全部词。

    末行很短（<= merge_short 字符）时并回上一行 —— 否则会出现
    「Cytorhabdovirus sp. 'lycii'」/「(GCVA)」这种孤行。
    """
    n = pretty_virus(s, 400)
    if len(n) <= width or maxlines <= 1:
        return n
    lines, cur = [], ""
    for w in n.split(" "):
        cand = f"{cur} {w}".strip()
        if cur and len(cand) > width and len(lines) < maxlines - 1:
            lines.append(cur)
            cur = w
        else:
            cur = cand
    if cur:
        lines.append(cur)
    if len(lines) == 2 and len(lines[1]) <= merge_short:
        lines = [f"{lines[0]} {lines[1]}"]
    return "\n".join(lines)


# --------------------------------------------------------------------------
# 面板 E / F（变异分析）—— 嵌入 04_post_analysis 的成品图
# --------------------------------------------------------------------------
# 为什么"嵌入"而不是重绘：`snpgenie` 那套（逐基因 dN/dS、iSNV 密度、突变谱）**没有
# 落任何 CSV**，重算要重跑整套 SNPGenie；`vcf_viz` / `snpeff_macro` 的 CSV 只能还原
# 其中一部分。嵌入是唯一一次拿全的做法。
#
# **两条实测硬约束**（决定只能通栏、不能半栏，与组合图设计宽无关）：
#   源图 3000 px 宽 = 300 dpi 下的 10 in 画布，标题 12 pt 渲染成 50 px。
#     通栏展示：印刷字号 = 12 × (D/10) × (7.09/D) = 12 × 0.709 ≈ **8.5 pt** ✓
#     半栏展示：12 × 0.354 ≈ **4.25 pt** ✗
#   （D = 组合图设计宽，7.09 in = 期刊 180 mm 双栏；两项里的 D 约掉了，
#     所以这是"比性质"，把组合图整体放大也不会变好。）
#   4800 px 那批（16 in 画布）连通栏也只有 ≈5.3 pt，故**优先选 3000 px 那批**。
#
# 裁白边是必须的：源图自带整页白边 + 自带 "Figure N: ..." 标题栏，不裁会白占
# 大片纵向空间、而且那个 "Figure N" 会和面板字母打架。

# 仓库根目录（`<repo>/scripts/barbarum/本文件` 往上三级）。
# `EF_LIBRARY` 里以 `{repo}/` 开头的路径按它展开 —— 用于**随仓库一起走的固定素材**
# （如用户亲手挑的 4 张 PPT 成品图），既不受 `04_post_analysis` 结构影响，
# 也不依赖 `--root`，保证任何机器上克隆下来都能复现。
REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))

EF_LIBRARY = {
    # key              : (相对 04_post_analysis/<virus>/ 的路径, 面板标题)
    # =============================================================================
    # 2026-10-10 把 45 张候选**全部下载逐张看过**之后的实测判定。
    # 通栏高 = 9.775 in × h/w（用**源图未裁**的像素算是个上界；实际嵌入前会
    #   _crop_white 裁白边 + 切源图自带标题，裁后值看出图日志 "通栏高 x.xx in"）。
    # 两条硬约束：**内容干净**（无伪影 / 无干扰柱 / 不是长度的代理变量）
    #   + **通栏高 ≲ 6 in**。正文只放 2 张（见 docstring），其余进补充或淘汰。
    # 完整论据见 `reports/sci_composite/README_SCI_composite.md` §九。
    # =============================================================================

    # ---- ✓ 正文可用：内容干净、尺寸合格 -----------------------------------------
    # ★E 默认｜10,000× bootstrap 的逐基因 dN/dS + 95% CI，中性线画在 1.0；
    #   6 个基因全部 *** 显著低于 1 → 净化选择。**天然免疫"比值分母趋零"的伪影**。
    "boot_dnds":     ("snpgenie/Fig06_Bootstrapped_dNdS.png",
                      "Per-gene dN/dS with 95% bootstrap CI (neutral line at 1.0)"),
    # ★F 默认｜intra-host 等位频率谱：**双峰**（低频准种 ~3.9k vs AF=1.0 固定 ~14k）。
    "afs":           ("vcf_viz/Figure5_AFS.png",
                      "Intra-host allele frequency spectrum (low-frequency quasispecies vs fixed)"),

    # ---- △ 有真结论但需重绘 / 或通栏太高 ----------------------------------------
    # 逐基因 dN/dS 的**分布**（小提琴+箱线）。结论对（G 最松），但 y 轴被一个
    # dS≈0 的伪影（29.1）拉到 31，其余 5 个基因压成一条线 → 重绘（log y 或截断）后才可进正文。
    "dnds_gene":     ("snpgenie/Fig05_Gene_dNdS_Stats.png",
                      "Per-gene dN/dS distribution (y-axis dominated by one dS~0 outlier)"),
    # dN vs dS 散点 + dN=dS 中性线：**净化选择最直接的可视化**（所有点在线的下方一侧）。
    # 但源图 1811x1861（含 marginal 直方图）→ 通栏 10.04 in，**放不进正文**，只能单占一页。
    "interhost_dnds": ("snpgenie/Fig03_InterHost_dNdS.png",
                       "Inter-host dN vs dS with dN=dS neutral line"),
    # 上一张的**宿主内对偶版**（piN vs piS），同样全部在中性线下方、同样 10.04 in。
    "intrahost_pi":  ("snpgenie/Fig04_IntraHost_pi.png",
                      "Intra-host piN vs piS with piN=piS neutral line"),
    # 逐样本「准种(5-80%) vs 共识(>80%)」变异构成。尺寸最省（通栏仅 3.67 in），
    # 但 x 轴 ~110 个样本名密到不可读；作为紧凑补充图可以。
    "intrahost_div": ("snpeff_macro/Figure_3_IntraHost_Diversity.png",
                      "Per-sample quasispecies vs consensus variant composition"),
    # PCA 2D（PC1 88.8% / PC2 9.6%，K=3 演化类群）。干净可读，但通栏 7.57 in 偏高。
    "pca_2d":        ("snpgenie/Fig11a_PCA_2D.png",
                      "PCA of viral populations (K=3 evolutionary classes)"),
    # 基因组景观一览。**2026-10-11 修：改用原始版**。
    # 原指向 `Fig10_Trinity_Landscape_SCI.png`，实测该"SCI 重制版"把基因组长度标签
    # (`OR489165.1 / 14812 nt`) 压在 N 基因色块上、π 轴 "Neutral" 注释被截断 —— 反而更差。
    # 原始版基因盒与文字不重叠，含正确的逐基因 intra-host dN/dS（G 0.318 … M 0.038）。
    # 缺点：与主图面板 B 重叠 3/4（基因图 + 滑窗 SNP + pi 都重合），所以只能进补充。
    "trinity_sci":   ("snpgenie/Fig10_Trinity_Landscape.png",
                      "GCVA evolutionary landscape and per-gene intra-host dN/dS (original)"),
    # 共享变异 AF 层次聚类。**2026-10-11 修：改用原始版**。
    # 原指向 `Figure4_Clustermap_SCI.png`，实测 220 样本行标签糊成一片（3531x4165 位图）。
    # 原始版标签清晰可读。
    "clustermap_sci": ("vcf_viz/Figure4_Clustermap.png",
                       "Top-50 shared variants AF clustermap (original)"),
    # 逐位点热点（Top intra-host variants）。尺寸合适，可作面板 B 的细节补充。
    "hotspots":      ("snpgenie/Fig08_Top_Hotspots.png",
                      "Top intra-host variant hotspots"),

    # ---- ✗ 不建议：信息被量纲压扁 / 是长度的代理 / 术语错位 ----------------------
    # 前 4 根柱是 A>A / T>T / C>C / G>G —— **与参考一致的位点，不是替换**，
    # 把 y 轴顶到 8000，真正的 Ts/Tv 被压到底部 1/5 → 重绘剔除 self 柱后才可用。
    "mut_spectrum":  ("snpgenie/Fig07_Mutational_Spectrum.png",
                      "Mutational spectrum and codon buffer bias"),
    # 逐基因累计变异数 ≈ **基因长度**（L 21k vs P4 1.5k），未归一化 → 无独立信息。
    "evol_load":     ("snpeff_macro/Figure_2_Gene_Payload.png",
                      "Per-gene evolutionary load (variants by IMPACT)"),
    # 标题称 depth-adjusted，**但图上没有深度轴** → "与深度无关"这个自证不成立。
    "isnv_density":  ("snpgenie/Fig02_Adjusted_iSNV_Density.png",
                      "Depth-adjusted intra-host iSNV density by viral gene"),

    # ---- 同义重复 / 与主图面板 B 重叠 / 其它（登记只为复现，不推荐）--------------
    "vaf":           ("snpgenie/Fig01_VAF_Spectrum.png",
                      "Variant allele frequency spectrum of iSNVs"),
    "af_violin":     ("vcf_viz/Figure7_AF_Violin.png",
                      "Allele-frequency distribution per sample"),
    "afs_merge":     ("vcf_merge/figs/afs.png", "Allele frequency spectrum (merged)"),
    "titv_pie":      ("vcf_viz/Figure2_TsTv_Pie.png",
                      "Transition/transversion spectrum (pie)"),
    "titv_summary":  ("maftools/maftools_06_TiTv_Summary.png",
                      "Transition/transversion summary (maftools)"),
    "trinity":       ("snpgenie/Fig10_Trinity_Landscape.png",
                      "The Trinity evolutionary pan-landscape (original)"),
    "manhattan":     ("snpeff_macro/Figure_1_Manhattan_Mut_Landscape.png",
                      "Genomic mutational landscape of viral quasispecies"),
    "var_landscape": ("vcf_viz/Figure1A_All_Variants_Landscape.png",
                      "High-quality variants: genome position x allele frequency"),
    "var_landscape_b": ("vcf_viz/Figure1B_Top50_Variants_Landscape.png",
                        "Top-50 variants: genome position x allele frequency"),
    "var_density":   ("vcf_viz/Figure6_Variant_Density.png",
                      "Variant density across the genome"),
    "lollipop":      ("maftools/maftools_03c_WholeGenome_Lollipop_Map.png",
                      "Non-synonymous mutation map"),
    "lollipop_all":  ("maftools/maftools_03b_Lollipop_All_Stitched.png",
                      "Per-gene lollipop plots (all stitched; 16.3 in tall)"),
    "oncoplot":      ("maftools/maftools_02_Oncoplot.png",
                      "Mutational co-occurrence matrix (TCGA oncoplot)"),
    "quasi_oncoplot": ("snpeff_macro/Figure_5_Ultimate_Quasispecies_OncoPrint.png",
                       "Quasispecies co-occurrence matrix (OncoPrint)"),
    "maf_summary":   ("maftools/maftools_01_mafSummary_TCGA.png",
                      "Variant summary (maftools)"),
    "pca":           ("vcf_viz/Figure9_PCA_Lineages.png", "PCA of viral populations"),
    "pca_3d":        ("vcf_viz/Figure9_PCA_Lineages_3D.png",
                      "PCA of viral populations (3D)"),
    "pca_merge":     ("vcf_merge/figs/pca.png", "PCA of viral populations (merged)"),
    "clustermap":    ("vcf_viz/Figure4_Clustermap.png",
                      "Population structure clustermap (original)"),
    "lineage_clust": ("snpeff_macro/Figure_4_Lineage_Clustermap.png",
                      "Lineage clustermap"),
    "epinet":        ("vcf_merge/figs/epistatic_network.png",
                      "Epistatic network of co-occurring variants"),
    "mantel_ibd":    ("vcf_merge/figs/mantel_ibd.png",
                      "Isolation-by-distance (Mantel test)"),

    # ---- ★ 用户 2026-10-12 亲自指定的 4 张（`--append-figs` 的默认值）------------
    # 这 4 张**不是** `04_post_analysis/` 的流水线产出，而是用户从论文章节目录
    # `gcva_ppt_figs/` 里**亲自挑的成品图**（用户原话：「我特别选的，就听我的」）。
    # 故拷进仓库固定路径 `reports/sci_composite/_append_src/`（`{repo}/` 前缀，
    # 见 `resolve_ef`），不依赖 `04_post_analysis` 或 `--root`，保证可复现。
    #
    # ⚠ **不要**把它们挂到既有 key 上：`oncoplot` / `lollipop` / `trinity` 在流水线里
    #   另有同名条目指向 `maftools/` 与 `snpgenie/` 的产出，**二者不是同一张图**
    #   （标题、画布、数值都不同）。挂过去会**误伤 `--ef variation` / `--var-panels`
    #   的既有行为**，故此处一律另起 `app_` 前缀。
    #
    # 来源（本地 `论文章节目录/gcva_ppt_figs/`，md5 已核对）：
    #   app_oncoplot  ← 14_2_Maftools_Oncoplot_01.png              11520×10368
    #                   "Altered in 180 (100%) of 180 samples." 6 基因 P/N/P4/M/L/G
    #   app_circos    ← 15_DVG_Circos_4Track.png                    3000×3000
    #                   "Recombination Mutational Landscape"（中心 Recombination Jump）
    #   app_lollipop  ← 14_3_WholeGenome_Lollipop_01.png           18432×6912
    #                   "Holo-Genome Architectural Polypeptide Mutation Hubs"（轴 0–15000）
    #   app_trinity   ← 14_6_Super_Evolutionary_Landscape.png       6342×3958
    #                   "Genomic & Evolutionary Pan-Landscape: ...OR489165.1"
    #
    # ⚠ 源图都是整页高分辨率位图（3000–18432 px 宽），按 9.775 in 通栏嵌入后，
    #   若按"源为 300 dpi"折算，印刷字号偏小（日志里的 `N.N pt` 即此口径）。
    #   用户明确要求**按原图嵌入、不要重绘**，故此处**不做重绘、也不降采样**；
    #   字号口径如实记在 README §10.12，不在此处"优化"。
    "app_oncoplot":  ("{repo}/reports/sci_composite/_append_src/"
                      "oncoplot_maftools_180samples.png",
                      "Mutational co-occurrence matrix (oncoplot, 180 samples)"),
    "app_circos":    ("{repo}/reports/sci_composite/_append_src/"
                      "circos_recombination_landscape.png",
                      "Recombination mutational landscape"),
    "app_lollipop":  ("{repo}/reports/sci_composite/_append_src/"
                      "lollipop_wholegenome_mutation_hubs.png",
                      "Holo-genome architectural polypeptide mutation hubs"),
    "app_trinity":   ("{repo}/reports/sci_composite/_append_src/"
                      "trinity_genomic_evolutionary_panlandscape.png",
                      "Genomic and evolutionary pan-landscape"),
    # 2026-10-12 追加：maftools 变异汇总 dashboard（6 子图）。
    # ⚠ **必须用 PPT 版（12672×9216，4-27）**，不能用流水线版
    #   `maftools/maftools_01_mafSummary_TCGA.png`（3300×2400，10-11）—— 二者
    #   **不是同一张图**：PPT 版 Variant Classification 只有 4 类、x 轴 0–25000、
    #   Top10 是 L 94%/P 100%/N 99%/G 91%/M 97%/P4 99%；流水线版有 8 类
    #   （多 5'Flank/3'Flank/In_Frame_Del/In_Frame_Ins）、x 轴 0–8000、
    #   L 80%/P 94%/N 95%/G 72%/M 81%/P4 90%。**变异数差约 3 倍**（180 vs 144 样本口径）。
    #   用户 2026-10-12 贴的截图 = PPT 版，故取 PPT 版（与另外 4 张同属 4-27 那套）。
    "app_mafsummary": ("{repo}/reports/sci_composite/_append_src/"
                       "mafsummary_dashboard.png",
                       "Mutational summary dashboard (maftools)"),
}
EF_IMG_KEYS = tuple(EF_LIBRARY)

# **源图自带标题不能切**的 key：`_crop_white(drop_title=True)` 的判据是
# "段高小 + 横向覆盖率低 → 像一行文本就切"。maftools dashboard 的**第一段墨迹就是
# 三个子图标题**（Variant Classification / Variant Type / SNV Class），切了会丢真内容，
# 所以这一张显式关掉 `drop_title`。
KEEP_SRC_TITLE = {"app_mafsummary"}

# 用户 2026-10-12 亲挑的图 + 拼图行结构（`--append-figs` 的默认值）。
# 语法：`;` 分隔**行**，`,` 分隔**同一行内的面板**。
# 两个**自绘 token**（不读位图，走 `panel_*` 现画）：
#   `E` = 「GCVA 群体基因组景观」**三层**（π+D 双轴 / 逐基因 Miss-Syn / 基因条带；
#         见 `panel_E_land`）
#   `L` = 「全基因组突变枢纽棒棒糖图」（重绘自 PPT 原图，y = 单位点**去重**携带样本数）
# 面板字母按**读序**（行内从左到右、行间从上到下）顺延；自绘 token 也占字母。
#   rev.l/m 旧序（每张一行）：oncoplot / circos / lollipop / trinity → 已废弃
#   rev.n 旧序（E|lollipop|trinity + oncoplot|mafsummary + circos）→ 已被 rev.o 取代
#   rev.o 旧序（「EFG 三者重绘，放一起展示」+「E 与 G 合并成一张」）：
#         row1 = E(合并景观) | L(自绘棒棒糖)   row2 = oncoplot | mafsummary   row3 = circos
#   rev.r（用户 2026-10-12「去除 G，H I 在一行」）：
#         row1 = E(三层景观) | L(自绘棒棒糖)   row2 = mafsummary | circos
#         → **oncoplot（原 G）不再进主图**，字母重排为 A–H（原 H→G、原 I→H）。
#   rev.s（用户 2026-10-12「GH 和 EF 位置重新换一下」）：
#         row1 = mafsummary | circos   row2 = E(三层景观) | L(自绘棒棒糖)
#         → 把**位图行**与**自绘行**上下互换。字母按**读序**顺延，故互换后
#           新 E|F = 位图行（mafsummary|circos）、新 G|H = 自绘行（E|L）——
#           字母跟着**位置**走、不跟内容走（这是 `--append-figs` 的既定约定）。
#         ⚠ 两行行高之和不变（max(3.26,4.05) + max(3.05,3.05) = 4.05 + 3.05 = 7.10 in）
#           → `H_extra` / `fig_h` / 输出像素尺寸**与 rev.r 完全相同**，只有行序变了。
CHOSEN_ROWS = (("app_mafsummary", "app_circos"),
               ("E", "L"))
CHOSEN_FIGS = ";".join(",".join(_r) for _r in CHOSEN_ROWS)

# ---- 拼图行里**自绘面板**的显示高度（英寸）--------------------------------
# 位图面板的高度是「列宽 × 高/宽」算出来的，自绘面板的宽高比不固定（层数、图例行数
# 都是版式定的）→ 必须显式给高度。两者不同高，行内取 max 后矮的那个**垂直居中**。
PACK_E_H = 3.05     # `E` = GCVA 群体基因组景观（**三层**）
# ↑ 与 `PACK_L_H` **同高是免费的**（行高取 max），所以直接给满 —— 三层各分到的高度
#   比 rev.q 的四层宽裕（rev.q 最矮的层只有 ~0.63 in → ylabel 必须用短版；rev.r 三层
#   是 1.66 / 0.86 / 0.33 in）。条带那 0.33 in 只是轴框，x 刻度与 xlabel 画在框外。
PACK_L_H = 3.05     # `L` = 全基因组突变枢纽棒棒糖图（自带一行 4 列图例）
# ↑ 与 `PACK_E_H` **同高是免费的**：行高取 max，L 矮了也不会让整幅变矮，白白浪费
#   纵向空间（棒棒糖的 hub 标签要 5 级错落，越高越好）。


def _crop_white(path, tol=245, pad=8, drop_title=True, title_frac=None):
    """把 04 出的整页图裁到墨迹外接框；`drop_title` 再把自带的 "Figure N: ..." 切掉。

    返回 PIL.Image。裁切同时**压高度**（去掉整页白边）与**提字号**（同样的像素宽
    塞进更窄的框 = 放大），是嵌入前的必要一步。

    **2026-10-10 修：源图标题可能是两行**。实测 `vcf_viz/Figure5_AFS` 的标题是
    "Figure 5: Allele Frequency Spectrum" + "(OR489165.1)" 两行；旧版只切第一段墨迹，
    于是第二行 `(OR489165.1)` 留在面板上，**字号比面板自己的标题还大**，看着像重复。
    现改成**按"段"循环切**，判据是"这一段像不像一行文本"：
        段高 < `_kmax` px（文本行）**且** 横向覆盖率 < 0.85
    轴框的上边线是一条长横线（覆盖率 ≈ 1.0）、图形本体的段高也远大于 `_kmax`，
    所以循环会在图形本体处自然停下。旧版的"只切第一段 + `i <= 130`"已被替换。

    **2026-10-12 修：`_kmax` 原本硬编码 100 px（= 12 pt @300 dpi 的两倍）**，
    对 `04_post_analysis` 那批 300 dpi 图够用，但对 `--append-figs` 的**高分辨率
    PPT 导出图**（宽 3000–18432 px）会失效 —— 实测标题带高
    `oncoplot 181 px / lollipop 241 px`，都 > 100 被当成"图形本体"，于是**源图标题
    原样留在面板上**（lollipop 那张的标题与我们的图注一字不差 → 明显重复）。
    新增 `title_frac`：给定时 `_kmax = max(100, title_frac × 图高)`，即**按图高比例**
    定阈值，与分辨率无关。**默认 `None` 时行为与旧版完全一致**（`_kmax = 100`），
    所以只有显式传 `title_frac` 的调用方（`--append-figs`）受影响 ——
    `--ef self/variation/deg` 的嵌入结果**构造上逐像素不变**。
    """
    from PIL import Image
    im = Image.open(path).convert("RGB")
    a = np.array(im.convert("L")) < tol
    ys = np.where(a.any(axis=1))[0]
    xs = np.where(a.any(axis=0))[0]
    if len(ys) == 0:
        raise ValueError(f"整张全白，无法裁切: {path}")
    y0, y1 = int(ys[0]), int(ys[-1])
    if drop_title:
        _kmax = 100 if title_frac is None else max(100, int(title_frac * im.height))
        # 一个"段" = 一段连续有墨迹的行。段落像文本行就切掉，像图形就停。
        for _ in range(4):                      # 最多切 4 段（主标题 + 副标题 + 容错）
            prof = a[y0:y1 + 1].any(axis=1)
            k = 0
            while k < len(prof) and prof[k]:    # 当前段的高度
                k += 1
            if k == 0 or k > _kmax:             # 到图形本体（或已无墨迹）→ 停
                break
            if a[y0:y0 + k].any(axis=0).mean() > 0.85:   # 像一条长横线（轴框顶）→ 停
                break
            j = k
            while j < len(prof) and not prof[j]:         # 段后的空白缝
                j += 1
            if j >= len(prof) or (j - k) < 6:   # 缝太薄 → 不是"标题|图形"的分界
                break
            y0 += j
    box = (max(0, int(xs[0]) - pad), max(0, y0 - pad),
           min(im.width, int(xs[-1]) + 1 + pad), min(im.height, y1 + 1 + pad))
    return im.crop(box)


def resolve_ef(key, ref, post_dir, log=print, root=None):
    """把 `EF_LIBRARY` 的 key（或直接给一个 .png 路径）解析成 (绝对路径, 标题)。

    路径前缀有两种，按序展开：
      `{root}/...` → `--root`（= 结果根目录）；用于不在 `04_post_analysis/<virus>/`
                     下的流水线产物（如 `01_data/.../Figure9_Circos_Landscape.png`）。
      `{repo}/...` → 仓库根（`REPO_DIR`）；用于**随仓库一起走的固定素材**
                     （`reports/sci_composite/_append_src/` 里用户亲挑的 4 张）。
    其余相对路径按 `post_dir`（`04_post_analysis/<virus>/`）拼接。
    """
    if ref in EF_LIBRARY:
        rel, cap = EF_LIBRARY[ref]
    elif str(ref).lower().endswith(".png"):
        rel, cap = str(ref), os.path.basename(str(ref))
    else:
        raise KeyError(f"--panel-{key.lower()} 不认得 {ref!r}；可用 key: {sorted(EF_LIBRARY)}")
    if rel.startswith("{root}/"):
        if not root:
            raise ValueError(f"{key!r} 的路径 {rel!r} 需要 `--root`，但调用方没给")
        rel = os.path.join(root, rel[len("{root}/"):])
    elif rel.startswith("{repo}/"):
        rel = os.path.join(REPO_DIR, rel[len("{repo}/"):])
    p = rel if os.path.isabs(rel) else os.path.join(post_dir, rel)
    if not os.path.exists(p):
        raise FileNotFoundError(p)
    return p, cap


def panel_var(fig, spec, png, caption, letter, log=print, im=None,
              disp_w=9.775, title_wrap=None):
    """把 04_post_analysis 的一张成品图当作组合图的一个面板（通栏，或拼图行的一列）。

    行高由调用方按裁切后的宽高比算好（见 `main()`），所以这里不会出现留白信件框。
    `im` 可直接给一张**已裁好**的 PIL 图（`--var-figure` 需要先量高度再建 figure，
    不能等到建完再裁）；给了 `im` 就不再看 `png`，`png` 只用于日志显示。

    `disp_w` / `title_wrap`（rev.n 新增，**默认值 = 旧行为**）：
      `disp_w` = 该面板**实际**占的显示宽度（英寸）。通栏时是 9.775；拼图行里是列宽。
               只用于日志里的"印刷字号"，不影响绘制。
      `title_wrap` = 图注按多少字符折行。窄栏里长图注（如
                "Holo-genome architectural polypeptide mutation hubs" 50 字符）
                会横向溢出、压到左邻面板上 —— 传了就 `textwrap.fill`。
    """
    ax = fig.add_subplot(spec)
    if im is None:
        im = _crop_white(png)
    ax.imshow(im, interpolation="antialiased")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if title_wrap:
        import textwrap
        caption = textwrap.fill(caption, title_wrap)
    ax.set_title(caption, fontsize=8.5, fontweight="bold", loc="left", pad=5)
    panel_letter(fig, ax, letter, dy=1.045)
    # 印刷字号：源图按 300 dpi 出，正文基准 12 pt；裁后像素宽铺满 `disp_w` in，
    # 再乘期刊 180 mm / 设计宽 11.5 in 的印前缩放。>= 6 pt 才算能进正文。
    pt = 12.0 * disp_w * 300.0 / im.width * (7.09 / 11.5)
    log(f"  [{letter}] 嵌入成品图: {os.path.basename(png) if png else '(预裁图)'} "
        f"({im.width}x{im.height} px, 宽高比 {im.width / im.height:.2f}, "
        f"显示宽 {disp_w:.2f} in → 高 {disp_w * im.height / im.width:.2f} in, "
        f"印刷字号 ≈ {pt:.2f} pt@180mm{'  OK' if pt >= 6 else '  ** 偏小 **'})")
    return ax


# --------------------------------------------------------------------------
# 面板 A —— 检测 & 组装漏斗（四段堆叠条；「检出」= filter 阶段的高置信结果）
# --------------------------------------------------------------------------
# 条长 = 该物种的 raw 命中 run 数；由左至右四段：
#   完整组装 (N% <= 5) → 送组装但 N% > 5 → 高置信检出未送组装 → 仅 raw 命中(未过高置信过滤)
# 主线三级 raw → 检出 → 完整组装 的数值与率写在条**右侧的标注列**里，
# 两个率的分母统一为队列 n = denom (445)。数字放在条外是这张图能用 6.3 pt
# 还读得清的关键 —— 气泡点阵把数字塞进泡里，26 行时只能给到 5.8 pt（见 docstring）。
# 「检出」的定义见 load_funnel_table()：用 02_filtering/high_conf.summary.tsv，
# 依据 auto_known_virus.py::get_summary()（filter 跑完下游一律用 high_conf）。
SEG_FULL = GREEN            # 完整组装 (N% <= 5)
SEG_PART = ORANGE           # 送组装但 N% > 5
SEG_DET = SKY               # 高置信检出但未送组装
SEG_RAW = "#DDDDDD"         # 仅 raw 命中（未过高置信过滤）

_ACC_TAIL = re.compile(r"([A-Z]{2}_?\d+\.\d+)$")


def _acc_of(v):
    """`extraction_status.tsv` 的 virus 列是 `_` 连接形式，取尾部 accession。"""
    m = _ACC_TAIL.search(str(v))
    return m.group(1) if m else None


def _run_of(s):
    """`CRR1126132_OR489165.1` → `CRR1126132`（提取产物把 run 与 accession 用 `_` 拼在一起）。"""
    return str(s).split("_")[0]


def load_funnel_table(bar_dir, log=print):
    """三层源表 → 逐物种四列计数（面板 A / 文字结论共用同一份口径）。

    返回 `(df, meta)`：
      df    index = 物种名，列 = `raw` / `det` / `asm_row` / `asm_ok`，按 raw 降序；
      meta  诊断用计数（raw 表规模、可用检出 run 数、送组装行数等）。

    **口径（2026-10-10 rev.d 定）：`det` = filter 阶段的高置信结果**
      raw      01_detection/summary/all_viruses.raw.tsv    `Species` 去重 `Sample`
      det      02_filtering/high_conf.summary.tsv          `Adjusted_Species` 去重 `Sample`
      asm_row  06_extraction/extraction_status.tsv         `status ∈ {extracted, skipped_rna_n_over}`
      asm_ok   同上                                        `status == extracted`（N% <= 5）

    依据：`virome_analysis_pipeline/auto_known_virus.py::get_summary()` 的约定 ——
    「filter 运行后自动切换到 high_conf」，其后的 variants / assembly / dvg / report
    全用它；所以正文图里「检出」应当是高置信过滤后的结果，而不是 `best`（Poisson 过滤）。
    行集 = 至少通过 high_conf 的物种（barbarum 实测 9 种）。

    **组装层与检出层按 (物种, run) 对取交集**：`extraction_status.sample` 是 `<run>_<accession>`
    形式，先经 `_run_of()` 归一到 run；**必须按对（不是按 run）取交**，否则一个 run 在别的
    病毒上过了高置信，就会把它在本物种的组装记录也带进来（实测 TCDVd 40>39、GARV4 5>3）。
    barbarum 实测：366 行 extraction 记录里 349 行命中 high_conf 的 (物种, run)，
    余下 17 行（12 个 run）全是 `skipped_rna_n_over`（N%/长度未达标），故取交集不影响
    `asm_ok`（303 条一条不少），只让漏斗保持严格单调。
    物种归属：raw 用自身 `Species` 列（实测与 `Adjusted_Species` 逐字同名，可直接 join）；
    `extraction_status` 用 virus 列尾部 accession 经 `Rep_Accession → Adjusted_Species` 映射
    （以 high_conf 为准、best 兜底 —— 有 2 个 acc（OR701334.1 / NC_003347.1）只在 best 里）。
    """
    raw = read_table(os.path.join(bar_dir, "01_detection/summary/all_viruses.raw.tsv"),
                     sep="\t", dtype=str)
    hc = read_table(os.path.join(bar_dir, "02_filtering/high_conf.summary.tsv"),
                    sep="\t", dtype=str)
    best = read_table(os.path.join(bar_dir, "01_detection/summary/all_viruses.best.summary.tsv"),
                      sep="\t", dtype=str)
    ext = read_table(os.path.join(bar_dir, "06_extraction/extraction_status.tsv"),
                     sep="\t", dtype=str)

    hc = hc.assign(run=[_run_of(s) for s in hc["Sample"]])
    n_det = hc.groupby("Adjusted_Species")["run"].nunique()
    n_raw = raw.groupby("Species")["Sample"].nunique()

    acc2sp = (best.drop_duplicates("Rep_Accession")
              .set_index("Rep_Accession")["Adjusted_Species"].to_dict())
    acc2sp.update(hc.drop_duplicates("Rep_Accession")
                  .set_index("Rep_Accession")["Adjusted_Species"].to_dict())
    ext = ext.assign(Species=[acc2sp.get(_acc_of(v)) for v in ext["virus"]])
    n_unmapped = int(ext["Species"].isna().sum())
    if n_unmapped:
        log(f"  [warn] extraction_status: {n_unmapped}/{len(ext)} 行无法映射到物种，已忽略")
        ext = ext[ext["Species"].notna()]

    ext = ext.assign(run=[_run_of(s) for s in ext["sample"]])
    n_ext_all = int(len(ext))
    hc_pairs = set(zip(hc["Adjusted_Species"], hc["run"]))
    keep = [(s, r) in hc_pairs for s, r in zip(ext["Species"], ext["run"])]
    dropped = ext[~np.array(keep)]
    ext = ext[np.array(keep)]                          # 组装层 ∩ 检出层（按对）
    n_dropped = n_ext_all - int(len(ext))
    if n_dropped:
        st = dropped["status"].value_counts().to_dict()
        log(f"  [A] 组装层∩检出层丢弃 {n_dropped} 行（{st}）")

    sent = ext[ext["status"].isin(("extracted", "skipped_rna_n_over"))]
    n_sent = sent.groupby("Species")["run"].nunique()
    n_ok = ext[ext["status"] == "extracted"].groupby("Species")["run"].nunique()

    df = pd.DataFrame({"det": n_det})
    df["raw"] = n_raw.reindex(df.index)
    df["asm_row"] = n_sent.reindex(df.index)
    df["asm_ok"] = n_ok.reindex(df.index)
    df = df.fillna(0).astype(int)

    bad = ((df["raw"] < df["det"]) | (df["det"] < df["asm_row"])
           | (df["asm_row"] < df["asm_ok"]))
    if bool(bad.any()):
        log(f"  [warn] {int(bad.sum())} 个物种违反 raw >= 检出 >= 送组装 >= 完整组装，"
            f"已按上级夹取（否则堆叠出现负宽）")
        df["det"] = np.minimum(df["det"], df["raw"])
        df["asm_row"] = np.minimum(df["asm_row"], df["det"])
        df["asm_ok"] = np.minimum(df["asm_ok"], df["asm_row"])

    df = df.sort_values(["raw", "det"], ascending=False)
    log(f"  [A] 检出=high_conf(filter)  species={len(df)}/{int(hc['Adjusted_Species'].nunique())}"
        f"  raw={int(df['raw'].sum())}/{len(raw)} 条（raw 表共 {int(raw['Species'].nunique())} 种）"
        f"  det={int(df['det'].sum())}  sent={int(df['asm_row'].sum())}  "
        f"assembled={int(df['asm_ok'].sum())}  "
        f"（组装层∩检出层丢弃 {n_dropped}/{n_ext_all} 行；映射失败 {n_unmapped}）")
    meta = dict(raw_rows=int(len(raw)), raw_species=int(raw["Species"].nunique()),
                det_species=int(hc["Adjusted_Species"].nunique()),
                det_runs=int(hc["run"].nunique()),
                sent_rows=int(len(sent)), ext_rows=n_ext_all,
                ext_dropped=n_dropped, unmapped=n_unmapped)
    return df, meta


def panel_A(fig, spec, df, denom, letter="A", log=print):
    """A —— 检测 & 组装漏斗：四段堆叠条 + 右侧 `raw → 高置信检出(率) → 完整组装(率)` 标注列。

    `df` 由 `load_funnel_table()` 给出（main() 里先取一次，是为了按行数定版面高度）。
    """
    ax = fig.add_subplot(spec)
    if df is None or not len(df):
        ax.axis("off")
        ax.set_title("[A] funnel table missing", fontsize=8.5, loc="left")
        return ax

    n = len(df)
    y = np.arange(n)[::-1].astype(float)
    left = np.zeros(n)
    segs = ((df["asm_ok"].to_numpy(float), SEG_FULL, "Fully assembled  N% \u2264 5"),
            ((df["asm_row"] - df["asm_ok"]).to_numpy(float), SEG_PART,
             "Assembled but N% > 5"),
            ((df["det"] - df["asm_row"]).to_numpy(float), SEG_DET,
             "High-conf. detected, not assembled"),
            ((df["raw"] - df["det"]).to_numpy(float), SEG_RAW,
             "Raw hit only (below filter)"))
    handles, labels = [], []
    for w, c, lab in segs:
        ax.barh(y, w, left=left, height=0.70, color=c, linewidth=0, zorder=3)
        # 全 0 的段不进图例（barbarum 实测「检出未送组装」= 0：所有高置信检出都送进了组装）
        if float(w.sum()) > 0:
            handles.append(Rectangle((0, 0), 1, 1, fc=c, ec="none"))
            labels.append(f"{lab}  ({int(w.sum())})")
        left = left + w

    # --- 右侧标注列：主线三级数值 + 两个率（分母统一 denom） ---
    # rev.c/d：面板 A 与 D 同行（轴宽 ~5.7 in），标注列可用宽度 = 0.339 × 轴宽 ≈ 1.95 in，
    # 表头用短版；末尾另有溢出安全阀（见下）。
    xmax = float(df["raw"].max())
    ann = []
    for yy, r, d, a in zip(y, df["raw"].to_numpy(), df["det"].to_numpy(),
                           df["asm_ok"].to_numpy()):
        ann.append(ax.text(xmax * 1.045, yy,
                           f"{int(r)}  \u2192  {int(d)} ({100.0 * int(d) / denom:.1f}%)  \u2192  "
                           f"{int(a)} ({100.0 * int(a) / denom:.1f}%)",
                           ha="left", va="center", fontsize=6.3,
                           color="#222222", zorder=4))
    ann.append(ax.text(xmax * 1.045, n - 0.5 + 0.58,
                       "raw  \u2192  high-conf. (%)  \u2192  assembled (%)",
                       ha="left", va="center", fontsize=5.9, color="#999999", zorder=4))

    ax.set_yticks(y)
    # 长名按词折行（不截断）—— rev.d 起物种名不再用 "…" 截断，改由 wrap_name() 两行呈现。
    # barbarum 实测：这 9 个名字 prettify 后最长 34 字符（"Cytorhabdovirus sp. 'lycii' (GCVA)"），
    # width=30 时其末行 "(GCVA)" 只有 6 字符 → 触发 merge_short 并回上一行 → 实际仍全部单行。
    # 这不是未生效，而是"能单行就不折行"的正确结果；名字再长（>=35 字符）才会真的折两行。
    ax.set_yticklabels([wrap_name(s, 30) for s in df.index],
                       fontsize=7.0, fontstyle="italic")
    ax.set_ylim(-0.75, n - 0.5 + 1.25)
    ax.set_xlim(0, xmax * 1.58)
    ax.set_xticks(np.arange(0, xmax + 1, 50))
    ax.grid(axis="x", color="#DDDDDD", lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # 安全阀：标注列若在更窄的版式下会溢出轴框（文字不裁剪，会压到右邻面板），
    # 就把 xlim 顶到最宽一行的右端。实测 6.3 pt 下需要轴宽 >= 4.31 in。
    try:
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        inv = ax.transData.inverted()
        right = max(inv.transform((t.get_window_extent(renderer=rend).x1, 0))[0]
                    for t in ann)
        if right > ax.get_xlim()[1]:
            ax.set_xlim(0, right * 1.02)
            log(f"  [A] 标注列偏长：xlim 由 {xmax * 1.58:.0f} 扩到 {right * 1.02:.0f}"
                f"（等价 {right * 1.02 / xmax:.2f}×xmax）")
    except Exception as exc:                       # 无 renderer 时跳过，不影响出图
        log(f"  [warn] 标注列溢出检测跳过: {exc}")

    ax.set_title("Detection and genome-recovery funnel (high-confidence set)",
                 fontsize=8.5, fontweight="bold", loc="left", pad=7)
    ax.set_xlabel(f"Runs per species holding a raw detection  (cohort n = {denom})",
                  fontsize=7, labelpad=21)
    # 图例段数自适应（barbarum 实测只有 3 段进图例：四段里的「检出未送组装」恒为 0）。
    # 锚点只是一个起点，随后由安全阀 2 按**实测文字盒**往下推 —— 硬编码的 -0.028
    # 在 rev.d 曾与 x 刻度的 "50 / 150" 叠字（图例与刻度不同长，ncol 一改就翻车）。
    leg_y = -0.115
    leg = ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.0, leg_y),
                    ncol=max(1, len(labels)), frameon=False, fontsize=6.4,
                    handlelength=1.1, columnspacing=1.3, handletextpad=0.5,
                    borderaxespad=0.0)

    # 安全阀 2：图例在轴框下方，与 x 刻度标签同处一条带。用实测文字盒把它压到刻度之下，
    # 再把 xlabel 的 labelpad 撑到图例之下（labelpad 自刻度标签底量起），否则图例与 xlabel 又叠。
    try:
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        tl_y0 = min(t.get_window_extent(renderer=rend).y0 for t in ax.get_xticklabels())
        for _ in range(24):
            if leg.get_window_extent(renderer=rend).y1 < tl_y0 - 2.0:
                break
            leg_y -= 0.02
            leg.set_bbox_to_anchor((0.0, leg_y))
            fig.canvas.draw()
            rend = fig.canvas.get_renderer()
        gap_pt = (tl_y0 - leg.get_window_extent(renderer=rend).y0) / fig.dpi * 72.0
        ax.xaxis.labelpad = max(21.0, gap_pt + 5.0)
        # 顺带把 y 刻度标签的实测宽度/左端位置打出来 —— 版式左边距靠它定，别再拿估算值凑
        ytw = max(t.get_window_extent(renderer=rend).width
                  for t in ax.get_yticklabels()) / fig.dpi
        yx0 = min(t.get_window_extent(renderer=rend).x0
                  for t in ax.get_yticklabels()) / fig.dpi
        log(f"  [A] legend anchor y={leg_y:.3f}  xlabel labelpad={ax.xaxis.labelpad:.1f}pt  "
            f"| ytick max_w={ytw:.2f}in left_in={yx0:+.2f}in "
            f"({'OK' if yx0 >= 0 else '溢出画布'})")
    except Exception as exc:
        log(f"  [warn] 图例/刻度重叠自检跳过: {exc}")

    panel_letter(fig, ax, letter, dy=1.026)
    return ax


# --------------------------------------------------------------------------
# 面板 B —— 群体内变异景观 (GCVA)
# --------------------------------------------------------------------------
GCVA_GENES = [("N", 25, 1939), ("P", 1951, 4043), ("P4", 4049, 5021),
              ("M", 5027, 6337), ("G", 6344, 7363), ("L", 7369, 14792)]
GCVA_LEN = 14812


def panel_B(fig, spec, gcva_dir, letter="B", log=print, compact=False):
    """面板 E —— GCVA 群体基因组景观（基因图 / 滑窗变异+π / 每基因 miss·syn）。

    `compact=False`（默认）＝ 通栏 9.775 in 的原始版式（rev.k / rev.m 用），
    **逐像素与旧版一致**。

    `compact=True`（rev.n 拼图行专用）：E 被塞进 `E|G|H` 拼图行的**窄列**（≈2.9 in），
    通栏版式在这宽度下会散架，故只做**三处收窄**（内容与数据一字不改）：
      ① 标题折两行并去掉 accession / 长度（118 字符的长标题会横跨到右邻列）；
      ② 左两条 ylabel 缩短 + 字号 6.8→6.3 pt（通栏时二者就贴着，窄列里会互撞）；
      ③ 右轴 ylabel 由 "Nucleotide diversity π (per site)" 缩成 "π (per site)"
         （通栏时它长约 1.7 in，窄列里顶端会蹭到右邻面板的字母）。
    图例与数据标注保持原样，窗口长度 590 bp 移进左 ylabel（"SNP count / 590 bp"）。
    """
    # 标题：通栏用全称；拼图窄列折两行、去掉 accession/长度（见 docstring ①）
    _TITLE = ("GCVA (Cytorhabdovirus sp. 'lycii', OR489165.1, 14,812 bp)  \u2014  "
              "population-genomic landscape (144 variant-carrying samples)")
    _YL_LAND = "SNP count per 590-bp window"
    _YL_SEL = "Missense / synonymous ratio"
    _YL_PI = "Nucleotide diversity \u03c0 (per site)"
    _FS = 6.8
    if compact:
        _TITLE = ("GCVA population-genomic landscape\n"
                  "(144 variant-carrying samples)")
        _YL_LAND = "SNP count / 590 bp"
        _YL_SEL = "Miss / Syn ratio"
        _YL_PI = "\u03c0 (per site)"
        _FS = 6.3
    sub = spec.subgridspec(3, 1, height_ratios=[0.20, 1.00, 0.52], hspace=0.10)
    ax_map = fig.add_subplot(sub[0, 0])
    ax_land = fig.add_subplot(sub[1, 0])
    ax_sel = fig.add_subplot(sub[2, 0])

    # --- B0: 基因图 -------------------------------------------------------
    ax_map.set_xlim(0, GCVA_LEN)
    ax_map.set_ylim(0, 1)
    ax_map.axis("off")
    for i, (g, s, e) in enumerate(GCVA_GENES):
        ax_map.add_patch(Rectangle((s, 0.22), e - s, 0.56,
                                   facecolor=OKABE[i % 7], edgecolor="#333333", lw=0.5))
        ax_map.text((s + e) / 2, 0.50, g, ha="center", va="center",
                    fontsize=6.5, color="white", fontweight="bold")
    ax_map.plot([0, GCVA_LEN], [0.13, 0.13], color="#333333", lw=0.8)
    ax_map.text(-0.004, 0.5, "Genome", transform=ax_map.transAxes, rotation=90,
                ha="center", va="center", fontsize=6.2, color="#333333")
    ax_map.set_title(_TITLE, fontsize=8.5, fontweight="bold", loc="left", pad=8)
    # 标题宽度自检（rev.n）：窄列里长标题会横跨到右邻列，而体检脚本**查不到**
    # （标题挂在 `axis("off")` 的 ax_map 上；`_check_main_overlap.py` 已补收 title，
    #  但那只报「文字撞文字」，撞到邻列**空白处**仍不报）。8.5 pt 粗体 ≈ 0.065 in/字符。
    _bw = ax_map.get_position().width * fig.get_figwidth()
    _tw = max(len(_ln) for _ln in _TITLE.split("\n")) * 0.065
    if _tw > _bw:
        log(f"  [warn] 面板 {letter} 标题估宽 {_tw:.2f} in > 列宽 {_bw:.2f} in，"
            f"会横跨到右邻面板 —— 请缩短标题或折行")

    # --- B1: 滑窗变异密度 + π --------------------------------------------
    # ⚠ 口径（**2026-10-12 修正**）：`PI` 列是**窗口内各位点 π 的求和**，不是 per-site π
    #   （`virome_analysis_pipeline/virus_variants_analyzer.py:248`
    #    `pi_window = in_window['Pi'].sum()`；对照 `virus_vcf_pipeline.py:1402` 是求均值，
    #    两套实现口径不同）。旧版**直接画 `PI`** 并把右轴标成 "π (×10⁻³)" →
    #   实际画出来的是 per-site π × 0.590，**系统性偏低 41 %**（实测比值 = 0.590）。
    #   现改为 ÷ 窗口长度，与 `panel_V_pi_tajima`（`_VariationRich` 面板 B）**同口径**：
    #   per-site π = 0.0034–0.0346（均值 0.0134）；末窗按实际覆盖长度截断，避免低估。
    win = read_table(os.path.join(gcva_dir, "vcf_viz/Table2_Sliding_Window_PopGen_Stats.csv"))
    mid = win["BIN_MID"].values
    _bs = pd.to_numeric(win["BIN_START"], errors="coerce").values.astype(float)
    _be = pd.to_numeric(win["BIN_END"], errors="coerce").values.astype(float)
    bw = float(np.nanmedian(_be - _bs))
    _Lw = np.minimum(_be, GCVA_LEN) - _bs + 1.0
    _Lw = np.where(np.isfinite(_Lw) & (_Lw > 0), _Lw, bw + 1.0)
    pi = pd.to_numeric(win["PI"], errors="coerce").values.astype(float) / _Lw  # -> per-site π
    ax_land.bar(mid, win["N_VARIANTS_S"], width=bw * 0.88, color=SKY,
                edgecolor="none", alpha=0.85, zorder=2, label="SNP count (window)")
    ax_land.set_ylabel(_YL_LAND, fontsize=_FS, color="#0072B2")
    ax_land.tick_params(axis="y", labelcolor="#0072B2")
    ax_land.set_xlim(0, GCVA_LEN)
    ax_land.set_ylim(0, win["N_VARIANTS_S"].max() * 1.18)

    ax_pi = ax_land.twinx()
    ax_pi.plot(mid, pi, color="#222222", lw=1.0, zorder=5,
               label="Nucleotide diversity \u03c0 (per site)")
    ax_pi.set_ylabel(_YL_PI, fontsize=_FS)
    ax_pi.set_ylim(0, float(np.nanmax(pi)) * 1.18)
    ax_pi.spines["top"].set_visible(False)
    log(f"  [{letter}] 滑窗 π: {len(mid)} 窗 × {bw + 1:.0f} bp; per-site π "
        f"[{np.nanmin(pi):.4f}, {np.nanmax(pi):.4f}] (均值 {np.nanmean(pi):.4f}); "
        f"原始求和列 PI [{win['PI'].min():.2f}, {win['PI'].max():.2f}] 已 ÷ 窗口长度")

    for _, s, e in GCVA_GENES:
        ax_land.axvline(s, color="#DDDDDD", lw=0.5, zorder=1)
    h1, l1 = ax_land.get_legend_handles_labels()
    h2, l2 = ax_pi.get_legend_handles_labels()
    ax_land.legend(h1 + h2, l1 + l2, loc="upper left", frameon=False, ncol=2)
    ax_land.set_xticklabels([])
    ax_land.spines["top"].set_visible(False)
    ax_pi.spines["right"].set_linewidth(0.7)
    ax_land.spines["right"].set_visible(False)

    # --- B2: 每基因 miss/syn ---------------------------------------------
    sel = read_table(os.path.join(gcva_dir, "snpeff_macro/Matrix_04_Gene_Selection_Pressure.csv"))
    order = {g: i for i, (g, _, _) in enumerate(GCVA_GENES)}
    sel = sel[sel["GENE"].isin(order)].copy()
    sel["_o"] = sel["GENE"].map(order)
    sel = sel.sort_values("_o")
    for _, r in sel.iterrows():
        g, s, e = GCVA_GENES[order[r["GENE"]]]
        ax_sel.bar((s + e) / 2, r["Miss_Syn_Ratio"], width=(e - s) * 0.72,
                   color=GREEN, edgecolor="#333333", lw=0.5, zorder=3)
        ax_sel.text((s + e) / 2, r["Miss_Syn_Ratio"] + 0.03, f"{r['Miss_Syn_Ratio']:.3f}",
                    ha="center", va="bottom", fontsize=5.8, color="#222222")
    ax_sel.axhline(1.0, color=VERM, lw=0.8, ls="--", zorder=4)
    ax_sel.text(GCVA_LEN * 0.90, 1.03, "Miss/Syn = 1", ha="right", va="bottom",
                fontsize=6, color=VERM)
    ax_sel.set_ylabel(_YL_SEL, fontsize=_FS)
    ax_sel.set_xlabel("Genomic coordinate (bp)", fontsize=7)
    ax_sel.set_xlim(0, GCVA_LEN)
    ax_sel.set_ylim(0, max(1.22, sel["Miss_Syn_Ratio"].max() * 1.22))
    ax_sel.set_xticks(np.arange(0, GCVA_LEN + 1, 2000))
    for s in ("top", "right"):
        ax_sel.spines[s].set_visible(False)

    panel_letter(fig, ax_map, letter, dx=-0.062, dy=1.30)
    return ax_land


# --------------------------------------------------------------------------
# 面板 E（rev.r）—— GCVA 群体基因组景观（三层，共用一条基因组 x 轴）
# rev.o 合并版 → rev.p 只留原三层 → rev.q 加 Tajima's D 成四层
# → rev.r 按用户「上部分换成 `_VariationRich` 面板 B、保留下部分、条带移到底部」重排
# --------------------------------------------------------------------------
def _shrink_x(fig, axes, pad_l, pad_r):
    """把一组轴**左右各内缩** `pad_l`/`pad_r`（占该列宽的比），把版面留给 ylabel/刻度。

    为什么必须内缩：`fig.add_subplot(subspec_cell)` 拿到的轴框 = 整格宽，而 ylabel、
    ytick、twinx 的右轴刻度与 ylabel **全部画在轴框之外**。拼图行里相邻两列的
    "左列右侧装饰" + "右列左侧装饰" 会挤进同一段行内缝 → 撞在一起。
    内缩后装饰落在本列自己的留白里，行内缝只需承担纯空白。
    """
    for _ax in axes:
        p = _ax.get_position()
        _ax.set_position([p.x0 + pad_l * p.width, p.y0,
                          p.width * (1.0 - pad_l - pad_r), p.height])


def _nat_h(kind, im, col_w):
    """拼图行里一个面板的**自然高度**（英寸）。

    * 自绘面板（`E` 合并景观 / `L` 棒棒糖）高度由我们**显式指定** —— 它们的宽高比
      不固定（层数、图例行数都是版式决定的），不能像位图那样按宽高比反算；
    * 位图面板高度 = `列宽 × 高/宽`（`imshow(aspect="equal")` 不留信件框的前提）。
    """
    if kind == "E":
        return PACK_E_H
    if kind == "L":
        return PACK_L_H
    return col_w * im.height / im.width


def panel_E_land(fig, spec, gcva_dir, letter="E", log=print, disp_w=9.775,
                 pack_cell=None):
    """rev.r 面板 E —— GCVA 群体基因组景观（**三层**，全自绘矢量）。

    `pack_cell`（rev.t 新增，**默认 `None` = 旧行为**）：拼图行里本面板所占的
    **整格** `(left_fig, width_fig)`（figure 分数）。给了就把**字母与图注锚到格左缘**
    而不是内缩后的轴框左缘 —— 否则会与上一行（位图行，轴框 = 整格）错位 0.48 in。
    详见 `panel_letter` 的 `anchor_x` 与 README §10.20。

    沿革：rev.o 四层（并进 trinity 的逐基因 dN/dS）→ rev.p 三层（基因盒 + 滑窗 SNP/π
    + 逐基因 Miss/Syn）→ rev.q 加第③层 Tajima's D 成四层 →
    **rev.r（2026-10-12）按用户「E 的上部分换成 `_VariationRich` 面板 B、保留下部分、
    基因组条带改到下面和 F 对齐」重构**。

    层（自上而下）：
      ① 滑窗 π(per site) + Tajima's D **双轴**（= `_VariationRich` 面板 B 的画法：
         π 走左轴「浅蓝柱 + 蓝线」，D 走右轴「橙线 + `D=0` 虚线 + neutrality 注释」）
      ② 逐基因 Miss/Syn 比 + 中性线 1.0
      ③ 基因条带（N/P/P4/M/G/L）+ x 轴刻度与标签 —— **从顶部挪到最下面**
    ⚠ rev.r 起**不再画「590 bp 滑窗 SNP 数」**（用户 2026-10-12 选定「不保留，完全按
      B 的双轴」）；位点沿基因组的分布仍可在 `_VariationRich` 面板 A 里读。

    ⚠ 条带为什么放最下面：`panel_L`（F 列）的基因盒也画在**自身轴框底部**（`y<0`），
      两列的条带因此落在同一水平带 —— 这就是用户说的「和 F 对齐」。x 轴刻度与
      `Genomic coordinate (bp)` 标签随条带一起下移（原来在 Miss/Syn 层）。

    口径（与 `panel_B` / `_VariationRich` 面板 B 完全一致，不另起一套）：
      * π(per site) = `vcf_viz/Table2_Sliding_Window_PopGen_Stats.csv::PI`
        **÷ 窗口长度**（`PI` 是**窗口内各位点 π 的求和**，见
        `virus_variants_analyzer.py:248`；末窗按实际覆盖截断）——
        实测 per-site π **0.0040–0.0346**（均值 0.0141）；
      * Tajima's D = 同表 `TAJIMA_D`（实测 −2.79…−0.36，**51/51 全为负**）；
      * Miss/Syn = `snpeff_macro/Matrix_04_Gene_Selection_Pressure.csv`。

    `disp_w` = 该面板实际占的显示宽度（英寸）；< 3.5 时自动切窄列版式（标题折行、字号降档）。
    """
    compact = disp_w < 3.5

    # ylabel 用**短版 / 单行版**：三层里 ② 只有 ~0.86 in，长 ylabel 旋转 90° 后会垂直
    # 溢出到邻层（rev.q 实测 2 对重叠）。
    # ⚠ π 那条**不能折成两行**：两行时靠轴的那一行会贴到 y 刻度标签上，且 `labelpad`
    #   对多行 ylabel 的锚点行为不一致 —— 实测 labelpad=4 → 4.3 px 重叠，labelpad=9 →
    #   **11.3 px**（越调越糟）。改用**单行长版**（与 `panel_V_pi_tajima` 逐字同款），
    #   字号降到 6.4 → 实测宽 ≈ 1.46 in < 顶轴高 1.66 in，装得下。
    _YL_PI = "Nucleotide diversity \u03c0 (per site)"
    _YL_TD = "Tajima's $D$"
    _YL_SEL = "Miss/Syn"
    if compact:
        _TITLE = ("GCVA population-genomic landscape\n"
                  "(144 variant-carrying samples)")
        _FS, _LS = 6.3, 5.6
    else:
        _TITLE = ("GCVA (Cytorhabdovirus sp. 'lycii', OR489165.1, 14,812 bp)\n"
                  "population-genomic landscape (144 variant-carrying samples)")
        _FS, _LS = 6.8, 6.4

    # 三层（自上而下）：① π + D 双轴 / ② 逐基因 Miss/Syn / ③ 基因条带 + x 轴
    # ⚠ 条带**必须是最下面一层**（用户 2026-10-12：「基因组条带改到下面，和 F 对齐」）。
    sub = spec.subgridspec(3, 1, height_ratios=[1.00, 0.52, 0.20], hspace=0.10)
    ax_top = fig.add_subplot(sub[0, 0])
    ax_sel = fig.add_subplot(sub[1, 0])
    ax_map = fig.add_subplot(sub[2, 0])
    # 左右各留出 ylabel / 刻度 / twinx 右轴的位置（见 `_shrink_x` 的说明）
    _shrink_x(fig, (ax_top, ax_sel, ax_map), 0.095, 0.075)

    # ---- ① 滑窗 π(per site) + Tajima's D（双轴 = `_VariationRich` 面板 B 的画法）---
    win = read_table(os.path.join(gcva_dir,
                                  "vcf_viz/Table2_Sliding_Window_PopGen_Stats.csv"))
    mid = win["BIN_MID"].values
    _bs = pd.to_numeric(win["BIN_START"], errors="coerce").values.astype(float)
    _be = pd.to_numeric(win["BIN_END"], errors="coerce").values.astype(float)
    bw = float(np.nanmedian(_be - _bs))
    _Lw = np.minimum(_be, GCVA_LEN) - _bs + 1.0
    _Lw = np.where(np.isfinite(_Lw) & (_Lw > 0), _Lw, bw + 1.0)
    pi = pd.to_numeric(win["PI"], errors="coerce").values.astype(float) / _Lw
    td = pd.to_numeric(win["TAJIMA_D"], errors="coerce").values.astype(float)

    ax_top.bar(mid, pi, width=bw * 0.90, color=SKY, edgecolor="none",
               alpha=0.45, zorder=2)
    ax_top.plot(mid, pi, color=BLUE, lw=1.0, zorder=3,
                label="Nucleotide diversity $\\pi$")
    for _, s, _e in GCVA_GENES:
        ax_top.axvline(s, color="#DDDDDD", lw=0.5, zorder=1)
    ax_top.set_xlim(0, GCVA_LEN)
    ax_top.set_ylim(0, float(np.nanmax(pi)) * 1.22)
    ax_top.set_ylabel(_YL_PI, fontsize=max(5.6, _FS - 0.4), color=BLUE)
    ax_top.tick_params(axis="y", labelcolor=BLUE)
    ax_td = ax_top.twinx()
    # ⚠ 必须显式同步 position：`twinx()` 复制的是 parent 的 **`_originalPosition`（内缩前）**，
    #   而 `_shrink_x` 只改了 `_position` → ax_td 会比 ax_top 宽（实测 w=0.39457 vs 0.32749）。
    #   后果：第①层的**右轴 spine、右轴刻度与 `D=0` 灰虚线**全部画到数据区之外，
    #   与第②③层（ax_sel / ax_map）的左右边界对不齐 —— 视觉上就是「上下宽度没对齐」。
    ax_td.set_position(ax_top.get_position())
    ax_td.plot(mid, td, color=VERM, lw=1.1, zorder=4, label="Tajima's $D$")
    ax_td.axhline(0.0, color="#999999", lw=0.8, ls="--", zorder=1)
    ax_td.set_ylabel(_YL_TD, fontsize=_FS, color=VERM)
    ax_td.tick_params(axis="y", labelcolor=VERM)
    _tlo, _thi = float(np.nanmin(td)), float(np.nanmax(td))
    ax_td.set_ylim(min(_tlo * 1.18, -0.4), max(_thi * 1.6, 0.4))
    # `D=0` 注释**贴在线上**（不能写固定轴分数 0.03 —— 右轴 ylim 是 (-3.29, 0.4)，
    #  `D=0` 在顶部，标 0.03 会飘到底部、与线脱节）。与 `panel_V_pi_tajima` 同法反算。
    _d0f = (0.0 - ax_td.get_ylim()[0]) / (ax_td.get_ylim()[1] - ax_td.get_ylim()[0])
    ax_td.text(0.997, min(_d0f + 0.02, 0.995), "neutrality  $D=0$",
               transform=ax_td.transAxes, ha="right", va="bottom",
               fontsize=5.4, color="#999999")
    ax_td.spines["top"].set_visible(False)
    ax_td.spines["right"].set_linewidth(0.7)
    h1, l1 = ax_top.get_legend_handles_labels()
    h2, l2 = ax_td.get_legend_handles_labels()
    ax_top.legend(h1 + h2, l1 + l2, loc="upper left", frameon=False, ncol=2,
                  fontsize=_LS)
    ax_top.set_xticklabels([])
    for _s in ("top", "right"):
        ax_top.spines[_s].set_visible(False)
    # rev.t：图注**锚到整格左缘**（`pack_cell`），与位图行的 `panel_var` 同一起点。
    # `loc="left"` 默认 x=0（= 内缩后的轴框左缘）→ 会把图注右移 0.478 in。
    _tx = 0.0
    if pack_cell is not None:
        _pp = ax_top.get_position()
        if _pp.width > 1e-9:
            _tx = (pack_cell[0] - _pp.x0) / _pp.width
    ax_top.set_title(_TITLE, fontsize=8.5, fontweight="bold", loc="left", pad=8,
                     x=_tx)
    _bw_in = ax_top.get_position().width * fig.get_figwidth()
    _avail_in = (pack_cell[1] * fig.get_figwidth()) if pack_cell is not None else _bw_in
    _tw = max(len(_ln) for _ln in _TITLE.split("\n")) * 0.065
    if _tw > _avail_in:
        log(f"  [warn] 面板 {letter} 标题估宽 {_tw:.2f} in > 可用宽 {_avail_in:.2f} in，"
            f"可能横跨到右邻面板")
    _n_neg = int(np.sum(td < 0))
    _n_ok = int(np.isfinite(td).sum())
    log(f"  [{letter}] 滑窗 π + Tajima's D: {len(mid)} 窗 × {bw + 1:.0f} bp; per-site π "
        f"[{np.nanmin(pi):.4f}, {np.nanmax(pi):.4f}] (均值 {np.nanmean(pi):.4f}); "
        f"原始求和列 PI [{win['PI'].min():.2f}, {win['PI'].max():.2f}] 已 ÷ 窗口长度; "
        f"D 范围 [{_tlo:.2f}, {_thi:.2f}], D<0 {_n_neg}/{_n_ok}")

    # ---- ② 逐基因 Miss/Syn ------------------------------------------------
    sel = read_table(os.path.join(gcva_dir,
                                  "snpeff_macro/Matrix_04_Gene_Selection_Pressure.csv"))
    order = {g: i for i, (g, _, _) in enumerate(GCVA_GENES)}
    sel = sel[sel["GENE"].isin(order)].copy()
    sel["_o"] = sel["GENE"].map(order)
    sel = sel.sort_values("_o")
    for _, r in sel.iterrows():
        g, s, e = GCVA_GENES[order[r["GENE"]]]
        ax_sel.bar((s + e) / 2, r["Miss_Syn_Ratio"], width=(e - s) * 0.72,
                   color=GREEN, edgecolor="#333333", lw=0.5, zorder=3)
        ax_sel.text((s + e) / 2, r["Miss_Syn_Ratio"] + 0.03, f"{r['Miss_Syn_Ratio']:.3f}",
                    ha="center", va="bottom", fontsize=5.8, color="#222222")
    ax_sel.axhline(1.0, color=VERM, lw=0.8, ls="--", zorder=4)
    ax_sel.text(GCVA_LEN * 0.995, 1.04, "Miss/Syn = 1", ha="right", va="bottom",
                fontsize=5.8, color=VERM)
    ax_sel.set_ylabel(_YL_SEL, fontsize=_FS)
    ax_sel.set_xlim(0, GCVA_LEN)
    ax_sel.set_ylim(0, max(1.24, sel["Miss_Syn_Ratio"].max() * 1.24))
    # x 轴刻度/标签已下移到条带层（③）→ 本层刻度连同刻度线一起去掉
    ax_sel.set_xticks([])
    for s in ("top", "right"):
        ax_sel.spines[s].set_visible(False)

    # ---- ③ 基因条带（挪到最下面，自带 x 轴）--------------------------------
    # ⚠ 不能用 `axis("off")` —— 这一层要承担整列的 x 轴刻度与 `Genomic coordinate (bp)`。
    ax_map.set_xlim(0, GCVA_LEN)
    ax_map.set_ylim(0, 1)
    for i, (g, s, e) in enumerate(GCVA_GENES):
        ax_map.add_patch(Rectangle((s, 0.26), e - s, 0.54,
                                   facecolor=OKABE[i % 7], edgecolor="#333333", lw=0.5))
        ax_map.text((s + e) / 2, 0.53, g, ha="center", va="center",
                    fontsize=6.5, color="white", fontweight="bold")
    ax_map.plot([0, GCVA_LEN], [0.18, 0.18], color="#333333", lw=0.8)
    ax_map.text(-0.004, 0.53, "Genome", transform=ax_map.transAxes, rotation=90,
                ha="center", va="center", fontsize=6.2, color="#333333")
    for _s in ("top", "right", "left", "bottom"):
        ax_map.spines[_s].set_visible(False)
    ax_map.set_yticks([])
    ax_map.set_xticks(np.arange(0, GCVA_LEN + 1, 2000))
    ax_map.tick_params(axis="x", labelsize=6.5, length=2.5, pad=2)
    ax_map.set_xlabel("Genomic coordinate (bp)", fontsize=7)

    # 字母画在轴框左上角（见 `panel_L` 的说明）。纵向偏移按**本层实际高度**换算 ——
    # 顶轴从 rev.q 的 0.26 in（条带层）涨到 ~1.66 in，沿用固定 `dy=1.30` 会把字母
    # 甩到标题上方。这里按"距轴框顶 0.30 in"反算轴分数。
    _axh_in = max(0.2, ax_top.get_position().height * fig.get_figheight())
    # rev.t：字母也锚到整格左缘（与位图行同一起点）—— 见 `panel_letter` 的 `anchor_x`
    _anc = (pack_cell[0] - 0.062 * pack_cell[1]) if pack_cell is not None else None
    panel_letter(fig, ax_top, letter, dx=-0.062, dy=1.0 + 0.30 / _axh_in,
                 anchor_x=_anc)
    return ax_top


# --------------------------------------------------------------------------
# 面板 F（rev.o 自绘）—— 全基因组突变枢纽棒棒糖图
# --------------------------------------------------------------------------
# 变异类别 → 图例名 + 颜色。类别映射与 maftools 口径一致；非编码的
# upstream/downstream 变体**剔除**（它们不在基因盒内，画上去只会污染 x 轴两端）。
_LOLLI_CLS = (
    ("Frame_Shift_Del", VERM),
    ("Frame_Shift_Ins", ORANGE),
    ("In_Frame_Del", PURPLE),
    ("In_Frame_Ins", YELLOW),
    ("Missense_Mutation", GREEN),
    ("Nonsense_Mutation", BLUE),
    ("Silent", "#999999"),
)
_LOLLI_COLOR = dict(_LOLLI_CLS)


def _lolli_class(effect, ref, alt):
    """`Matrix_01` 的 `EFFECT`（+ REF/ALT 长短）→ 图例类别；非编码返回 `None`。"""
    e = str(effect)
    if e == "synonymous_variant":
        return "Silent"
    if e == "missense_variant":
        return "Missense_Mutation"
    if e == "stop_gained":
        return "Nonsense_Mutation"
    if e == "frameshift_variant":
        return "Frame_Shift_Del" if len(str(alt)) < len(str(ref)) else "Frame_Shift_Ins"
    if e in ("disruptive_inframe_deletion", "conservative_inframe_deletion"):
        return "In_Frame_Del"
    if e in ("conservative_inframe_insertion", "disruptive_inframe_insertion"):
        return "In_Frame_Ins"
    return None


def _lolli_label(vid):
    """`N_639_p.Val68Val` → `N: p.Val68Val`（与源图 hub 标注**逐字符**同形：冒号后带空格）。"""
    s = str(vid)
    if "_p." in s:
        return s.split("_p.")[0].split("_")[0] + ": p." + s.split("_p.", 1)[1]
    return s


def panel_L(fig, spec, gcva_dir, letter="L", log=print, disp_w=9.775,
            pack_cell=None):
    """rev.o 面板 F —— 全基因组突变枢纽棒棒糖图（**自绘矢量**，重绘自用户 PPT 原图）。

    `pack_cell`（rev.t 新增，**默认 `None` = 旧行为**）：拼图行里本面板所占的
    **整格** `(left_fig, width_fig)`（figure 分数）。给了就把**字母与图注锚到格左缘**
    而不是内缩后的轴框左缘 —— 否则会与上一行（位图行，轴框 = 整格）错位 0.34 in。
    详见 `panel_letter` 的 `anchor_x` 与 README §10.20。

    用户 2026-10-12：「EFG 三者重绘，放一起展示」；y 轴口径经追问选定
    **「携带样本数（去重）」** —— 同一个位点即使多个样本都变异，只计一次。
    **rev.r**：图例从 cell 内的行改成 `fig.legend` 挂到 x 轴标签下方，让轴框独占
    cell（2.21 → 3.05 in）、基因盒与 E 列的条带落到同一水平带（详见函数体里的注释）。

    数据：`snpeff_macro/Matrix_01_Global_Melted_Variants.csv`（44,210 行 / **144 样本**）。
    y = 该位点上**去重后的携带样本数**（实测单位点最高 **127**，位点 3646）；
    茎按该位点的**优势类别**着色，点按各自的类别着色（同一位点可有多类）。

    ⚠ 与用户 PPT 原图是**内容级差异**（不是分辨率）：PPT 版 y 轴到 ~200，而当前队列只有
    **144** 个样本 → 去重后上限 127。流水线自带的
    `maftools/maftools_03c_WholeGenome_Lollipop_Map.png` 也是同一量级（~120）
    → **PPT 版是旧队列的产物**，本函数按当前队列重绘。
    """
    p = os.path.join(gcva_dir, "snpeff_macro/Matrix_01_Global_Melted_Variants.csv")
    m = read_table(p)
    m["POS"] = pd.to_numeric(m["POS"], errors="coerce")
    m = m[m["POS"].notna()].copy()
    m["POS"] = m["POS"].astype(int)
    m["_cls"] = [_lolli_class(e, r, a)
                 for e, r, a in zip(m["EFFECT"], m["REF"], m["ALT"])]
    _n_raw = len(m)
    m = m[m["_cls"].notna()].copy()
    n_samp = int(m["Sample"].nunique())
    n_pos = int(m["POS"].nunique())

    g = m.groupby(["POS", "_cls"])["Sample"].nunique().rename("n").reset_index()
    tot = m.groupby("POS")["Sample"].nunique().rename("n_tot").reset_index()
    ymax = float(tot["n_tot"].max())
    # 纵向留白要给**顶上的 hub 标签**腾地方（最多 5 级错落，见下）：
    #   最高点 = ymax + (0.05 + 0.10×4)×ymax = 1.45×ymax，再加一个标签高 ≈ 0.09×ymax
    #   → ytop = 1.72×ymax 才不顶到轴框。
    ytop = ymax * 1.72

    # ⚠ rev.r（用户 2026-10-12「改 F，让两列自然靠近」）：**图例不再占 cell 内的行**。
    #   旧版 `subgridspec(2,1,[1.00,0.16])` 把 0.35 in 图例 + 0.42 in 行距压在轴框下面
    #   → 轴框底被抬高 0.77 in → 基因盒（画在 `y<0`，贴着轴框底）比 E 列的条带高 0.85 in，
    #   两列条带对不上（实测 F 2545–2603 px vs E 2799–2850 px）。
    #   现在轴框**独占整个 cell**（高 2.21 → 3.05 in），图例改挂 `fig.legend` 放到
    #   **x 轴标签下方**（图外）—— 条带随之落到与 E 同一水平带，且 F 的绘图区白赚 38 %。
    sub = spec.subgridspec(1, 1)
    ax = fig.add_subplot(sub[0, 0])
    _shrink_x(fig, (ax,), 0.085, 0.030)

    # 基因盒画在 y<0（与源图一致：盒子挂在 0 线下），0 线自己画一条，底部 spine 关掉。
    # ⚠ rev.r：盒子的**厚度与高度**对齐 E 列条带。E 的条带轴高 0.3325 in、盒子占其
    #   0.26–0.80 → 距 cell 底 **0.0865–0.266 in（厚 0.18 in）**。轴框涨到 3.05 in 后，
    #   旧值 `(-0.150, -0.040)×ytop` 会画成 **0.27 in 厚**（比 E 粗 1.5 倍，并排看很显眼）。
    #   按「轴框底 = -0.150×ytop」反算：盒子取 `(-0.118, -0.051)×ytop` → 0.085–0.263 in。
    #   `set_ylim` 因此**不再写 `_b0*1.25`** —— 否则改 `_b0` 会连带把轴底一起挪走。
    _b0, _b1 = -ytop * 0.118, -ytop * 0.051
    for i, (gg, s, e) in enumerate(GCVA_GENES):
        # rev.t（用户 2026-10-12「H 的条带调成和作图条带一样的颜色」）：条带配色与 G 列
        #   （`panel_E_land` 的 `ax_map`，见上方同名循环）**逐字一致** —— 同为 `OKABE[i % 7]`
        #   + `edgecolor="#333333"` + `lw=0.5`。旧值是一律 `#4D4D4D` 深灰，
        #   与 G 的 Okabe-Ito 彩色条带上下并排时明显不是一套配色。
        ax.add_patch(Rectangle((s, _b0), e - s, _b1 - _b0, facecolor=OKABE[i % 7],
                               edgecolor="#333333", lw=0.5, zorder=4))
        ax.text((s + e) / 2, (_b0 + _b1) / 2, gg, ha="center", va="center",
                fontsize=6.2, color="white", fontweight="bold", zorder=5)

    # **每个（位点 × 类别）各一根茎**，茎色 = 该类别色 —— 与源图一致：
    # 高位点的茎是 Missense 的绿色，Silent 只贡献底下一小撮灰茎。
    # （旧写法是"一位点一根、按优势类别着色"，因 Silent 行数最多 → 满图灰茎，
    #   和源图观感差很远。）
    ax.vlines(g["POS"].values, 0, g["n"].values,
              colors=[_LOLLI_COLOR[c] for c in g["_cls"]], lw=0.45, alpha=0.40, zorder=2)
    for cls, _c in _LOLLI_CLS:
        s = g[g["_cls"] == cls]
        if s.empty:
            continue
        ax.scatter(s["POS"], s["n"], s=2.8, color=_c, edgecolor="none",
                   alpha=0.85, zorder=3, label=cls)
    ax.axhline(0, color="#333333", lw=0.9, zorder=1)

    # top-10 hub 标注（优先非 Silent）。
    # ⚠ 不能用「上一标签 > 1500 bp 就回 0 级」的**顺序**错落，也不能只用一维区间着色：
    #   rev.o 体检实测 —— 顺序 3 级 → 4 对互压；一维区间着色 → 还剩 2 对，
    #   因为相邻层的**纵向**间距（0.10×ymax）是按"两点同高"算的，而各 hub 的
    #   `n_tot` 差可达十几 → 上一层的标签会从下面顶进下一层。
    #   改成**真正的二维贪心**：逐个 hub 试第 0..N−1 层，只有当"横向区间"与
    #   "纵向区间"**同时**与已放置的标签相交时才换下一层。
    #   文字尺寸不估：用 `TextPath` 按当前字号精确量宽/高，再按轴框换算成数据单位。
    _vmap = m.drop_duplicates(["POS", "_cls"]).set_index(["POS", "_cls"])["Variant_ID"]
    # top-12（源图 PPT 版就是 12 个标注；`n_tot` 用**该位点全部类别**的去重样本数）
    _hub = (g[g["_cls"] != "Silent"].sort_values("n", ascending=False)
            .drop_duplicates("POS").merge(tot, on="POS")
            .sort_values("n_tot", ascending=False).head(12).sort_values("POS"))
    _FS_HUB = 5.2
    _axw_in = max(0.5, disp_w * (1.0 - 0.085 - 0.030))     # 轴框实际宽度（英寸）
    _axh_in = max(0.3, ax.get_position().height * fig.get_figheight())
    _bp_per_in = GCVA_LEN / _axw_in
    _y_per_in = (ytop * 1.150) / _axh_in      # 轴的数据跨度 = ytop − (−0.150×ytop)
    # 6 级：源图 PPT 版标注 12 个 hub，N/P 两基因里就挤了 8 个，5 级会打满。
    # 最高层底边 = 0.05+0.10×5 = 0.55×ymax，再加标签高 ≈ 0.07×ymax → 0.62×ymax，
    # 仍在 `ytop = 1.72×ymax` 的留白里。
    _LVL_N = 6
    _placed = []                                            # [(x, y_bot, half_bp, h_units)]
    _n_capped = 0
    for _, r in _hub.iterrows():
        _x = int(r["POS"])
        _txt = _lolli_label(_vmap.loc[(_x, r["_cls"])])
        _ext = TextPath((0, 0), _txt, size=_FS_HUB).get_extents()
        _half = 0.5 * _ext.width / 72.0 * 1.08 * _bp_per_in     # ×1.08 给渲染留余量
        _h = _ext.height / 72.0 * _y_per_in
        _ybot = None
        for _j in range(_LVL_N):
            _yb = r["n_tot"] + (0.05 + 0.10 * _j) * ymax
            _yt = _yb + _h
            if all(not (abs(_x - _px) < (_half + _ph)
                        and _yb < _pyb + _phh and _yt > _pyb)
                   for _px, _pyb, _ph, _phh in _placed):
                _ybot = _yb
                break
        if _ybot is None:                                   # N 级都放不下 → 最高层保底
            _n_capped += 1
            _ybot = r["n_tot"] + (0.05 + 0.10 * (_LVL_N - 1)) * ymax
        _placed.append((_x, _ybot, _half, _h))
        ax.text(_x, _ybot, _txt, ha="center", va="bottom",
                fontsize=_FS_HUB, color="#333333", zorder=6)
    if _n_capped:
        log(f"  [warn] 面板 {letter}: {_n_capped} 个 hub 标签在 {_LVL_N} 级内仍放不下，"
            f"已压到最高层（可能仍有轻微互压）")

    ax.set_xlim(0, GCVA_LEN)
    ax.set_ylim(-ytop * 0.150, ytop)
    ax.set_ylabel("Mutation count in cohort\n(distinct samples)", fontsize=6.4)
    ax.set_xlabel("Whole genomic position (bp)", fontsize=7)
    ax.set_xticks(np.arange(0, GCVA_LEN + 1, 2000))
    ax.tick_params(axis="x", length=0, pad=3)
    ax.spines["bottom"].set_visible(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    # 标题口径与 A–D/E 一致（8.5 pt 粗体、左对齐、pad=6），字母走 `panel_letter` 的默认 dy
    # rev.t：图注**锚到整格左缘**（`pack_cell`），与位图行的 `panel_var` 同一起点。
    _tx = 0.0
    if pack_cell is not None:
        _pp = ax.get_position()
        if _pp.width > 1e-9:
            _tx = (pack_cell[0] - _pp.x0) / _pp.width
    ax.set_title("Whole-genome mutational hubs", fontsize=8.5, fontweight="bold",
                 loc="left", pad=6, x=_tx)
    # 图例挂到**图外、本列 x 轴标签下方**：`loc="upper center"` + 轴框底再下移
    # `0.34 in`（要跨过 x 刻度 6.5 pt + tick pad + xlabel 7 pt 那 ~0.26 in）。
    # 用 `fig.legend`（不是 `ax.legend`）—— 图例必须画在轴框**之外**，`ax.legend`
    # 会被算进轴框的 tight bbox 反而把轴挤矮。
    h, l = ax.get_legend_handles_labels()
    _p = ax.get_position()
    fig.legend(h, l, loc="upper center",
               bbox_to_anchor=((_p.x0 + _p.x1) / 2.0,
                               _p.y0 - 0.34 / fig.get_figheight()),
               ncol=4, frameon=False, fontsize=5.7,
               handletextpad=0.35, columnspacing=1.1, markerscale=2.4)
    log(f"  [{letter}] 自绘棒棒糖: {_n_raw:,} 行 → 编码变异 {len(m):,} 行 / "
        f"{n_samp} 样本 / {n_pos:,} 位点；单位点最高 {ymax:.0f}（去重样本数）；"
        f"hub 标签用 {_LVL_N} 级区间着色（轴宽 {_axw_in:.2f} in, "
        f"{_bp_per_in:.0f} bp/in）")
    # rev.t：字母也锚到整格左缘（与位图行同一起点）—— 见 `panel_letter` 的 `anchor_x`
    _anc = (pack_cell[0] - 0.062 * pack_cell[1]) if pack_cell is not None else None
    panel_letter(fig, ax, letter, dx=-0.062, dy=1.045, anchor_x=_anc)
    return ax


# --------------------------------------------------------------------------
# 面板 C —— 共感染结构
# --------------------------------------------------------------------------
def panel_C(fig, spec, bar_dir, denom, letter="C", log=print):
    ax = fig.add_subplot(spec)
    up = read_table(os.path.join(bar_dir, "01_detection/summary/all_viruses.best.summary.tsv"), sep="\t")
    per_sample = up.groupby("Sample")["Adjusted_Species"].nunique()

    maxk = int(per_sample.max())
    counts = [denom - len(per_sample)]                      # 0 = 未检出
    counts += [int((per_sample == k).sum()) for k in range(1, maxk + 1)]
    ks = np.arange(len(counts))
    pct = [100 * c / denom for c in counts]

    cols = ["#DDDDDD"] + [plt.cm.viridis(x) for x in np.linspace(0.15, 0.9, maxk)]
    ax.bar(ks, counts, color=cols, edgecolor="white", lw=0.6, width=0.72, zorder=3)
    for k, c, p in zip(ks, counts, pct):
        ax.text(k, c + denom * 0.012, f"{c}\n{p:.1f}%", ha="center", va="bottom",
                fontsize=5.9, color="#222222")
    ax.set_xticks(ks)
    ax.set_xticklabels([str(k) for k in ks])
    ax.set_xlabel("Number of virus species per sample", fontsize=7)
    ax.set_ylabel("Number of samples", fontsize=7)
    ax.set_ylim(0, max(counts) * 1.22)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#DDDDDD", lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    n_det = int((per_sample >= 1).sum())
    n_multi = int((per_sample >= 2).sum())
    ax.set_title("Co-infection structure", fontsize=8.5, fontweight="bold", loc="left", pad=6)
    ax.text(0.98, 0.94,
            f"detected {n_det} ({100*n_det/denom:.1f}%)\n"
            f"single {int((per_sample==1).sum())} ({100*int((per_sample==1).sum())/denom:.1f}%)\n"
            f"co-infected (\u22652) {n_multi} ({100*n_multi/denom:.1f}%)",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.2, color="#222222",
            bbox=dict(boxstyle="round,pad=0.32", fc="#F5F5F5", ec="#DDDDDD", lw=0.5))
    panel_letter(fig, ax, letter)
    return ax


# --------------------------------------------------------------------------
# 面板 D —— 序列同一性 (SDT)：只出「同一性分布」，与面板 A 同行右侧窄栏
# --------------------------------------------------------------------------
def _read_sdt(path):
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    m = df.values.astype(float)
    iu = np.triu_indices_from(m, k=1)
    return m, m[iu][~np.isnan(m[iu])]


def panel_D(fig, spec, bar_dir, letter="D", log=print):
    """D —— 4 个病毒的成对全基因组 NT 同一性分布（窄栏，与 A 同行）。

    **只画分布，不画 GCVA 的 128×128 热图**（v3→rev.c 按需删）：
    该矩阵 96–100% 区间几乎饱和，结构信息已由分布曲线 + 80/90% 分界线承载；
    热图塞进窄栏也读不出结构，反而挤掉分布图的宽度。

    **数据源（2026-10-10 确认）**：只有 **`known_virus_barbarum` 树**（老编号 7 stage）
    下才有 `07_similarity/`；`known_virus_all_v2` 是 `08_similarity`。
    所以本面板的条数是 **109 / 140 / 3 / 3**（barbarum 树），不是 all_v2 的 128 / 172 / 3 / 5。

    **这 4 个病毒怎么来的**：`07_similarity/<病毒>/Mode_Filter/Full_Dataset/` 下
    只有「过 `trim_ends()` 后 `len >= 150` 且 `N% <= 5.0` 的共识 **>= 3 条**」
    的病毒才写了 `overall_WholeGenome_NT.fasta` 与 `*_sdt_exact.csv`
    （`virus_auto_pipeline.py` 的 `:1057 if len(over_recs) >= 3:`）。
    barbarum 树 9 个病毒目录里只有这 4 个过闸；TCDVd / GARV4 之所以只有 3 条，
    是上游 `03_variants/virus-consensus`（ivar：深度 < 5 的位点写 N）里
    42 / 2 条共识 N% 达 11–56% 被挡。**详见 README §八。**
    """
    ax_d = fig.add_subplot(spec)

    targets = [
        ("Cytorhabdovirus_sp._lycii_OR489165.1", "GCVA", OKABE[0]),
        ("Potato_spindle_tuber_viroid_NC_002030.1", "PSTVd", ORANGE),
        ("Tomato_chlorotic_dwarf_viroid_NC_000885.1", "TCDVd", VERM),
        ("Grapevine-associated_RNA_virus_4_MW648525.1", "GARV4", GREEN),
    ]
    base = "07_similarity/{}/Mode_Filter/Full_Dataset/02_similarity_matrices/overall_NT_only_sdt_exact.csv"

    dists, dims, labels, colors = {}, {}, [], []
    for d, lab, col in targets:
        p = os.path.join(bar_dir, base.format(d))
        if not os.path.exists(p):
            log(f"  [warn] 缺 SDT 矩阵: {p}")
            continue
        _m, v = _read_sdt(p)
        dists[lab] = v
        dims[lab] = int(_m.shape[0])          # 矩阵维度 = 参与 SDT 的**序列条数**
        labels.append(lab)
        colors.append(col)

    # 图例**双记法**（rev.f）：矩阵维度（序列条数）与成对数分开写。
    # 2026-10-10 追查「D 怎么就 3 对」时发现：只写 pairs 会被读成"就 3 条序列"，
    # 而 n 条序列 → C(n,2) 对，两个数量级在 GCVA 上是 109 vs 5,886，差 54 倍。
    for lab, col in zip(labels, colors):
        v = dists[lab]
        d_n = dims[lab]
        # n 很小（TCDVd / GARV4 各 3 对）时**不画直方图**：3 个点套 0.5 宽的 bin
        # 会画出两根针尖（占纵轴 55% / 33%），高度由 bin 宽度而非数据量决定，
        # 看上去像"这里有一个分布"—— 是误导。改成离散刻度线如实表达"只有 3 个值"。
        if len(v) < 10:
            ax_d.vlines(v, 0.0, 6.0, color=col, lw=2.0, alpha=0.9, zorder=4)
            ax_d.plot([], [], color=col, lw=2.0,             # 图例句柄
                      label=f"{lab} ({d_n} seqs \u2192 {len(v)} pairs, "
                            f"median {np.median(v):.1f}%)")
            continue
        bins = np.arange(80, 100.25, 0.5)
        h, edges = np.histogram(v, bins=bins, density=True)
        ax_d.plot(edges[:-1] + 0.25, 100 * h * 0.5, lw=1.2, color=col,
                  label=f"{lab} ({d_n} seqs \u2192 {len(v):,} pairs, "
                        f"median {np.median(v):.1f}%)")
    ax_d.set_xlabel("Pairwise whole-genome NT identity (%)", fontsize=7)
    ax_d.set_ylabel("Sequence pairs (%)", fontsize=7)
    ax_d.set_xlim(80, 100)
    ax_d.set_ylim(0, 58)
    ax_d.axvspan(80, 90, color="#F5F5F5", zorder=0)
    ax_d.text(85.0, 20, "intra-species\nsequence space", ha="center", va="top",
              fontsize=6, color=GREY)
    for thr, labt in ((90, "viroid demarcation 90%"), (80, "virus demarcation 80%")):
        ax_d.axvline(thr, color=GREY, lw=0.7, ls=":")
        ax_d.text(thr + 0.25, 26, labt, rotation=90, fontsize=5.6, color=GREY, va="top")
    for s in ("top", "right"):
        ax_d.spines[s].set_visible(False)
    ax_d.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.015, 0.985),
                handlelength=1.3, labelspacing=0.24, fontsize=6.0)
    ax_d.set_title("Pairwise nucleotide identity (SDT)", fontsize=8.5,
                   fontweight="bold", loc="left", pad=6)
    panel_letter(fig, ax_d, letter, dx=-0.10, dy=1.045)
    log(f"  [D] identity distribution only; series={labels}")
    return ax_d


# --------------------------------------------------------------------------
# 面板 E —— 宿主差异表达 (volcano)
# --------------------------------------------------------------------------
def panel_E(fig, spec, expr_dir, letter="E", log=print):
    ax = fig.add_subplot(spec)
    p = os.path.join(expr_dir, "deseq2_virus_vs_non.tsv")
    d = read_table(p, sep="\t")
    d = d.dropna(subset=["padj", "log2FoldChange"])
    d["padj"] = d["padj"].clip(lower=1e-300)
    d["nlp"] = -np.log10(d["padj"])
    sig = d["padj"] < 0.05
    lfc = d["log2FoldChange"]
    cls = np.where(sig & (lfc > 0), "up", np.where(sig & (lfc < 0), "down", "ns"))
    cmap = {"ns": "#DDDDDD", "up": UP, "down": DN}
    for k in ("ns", "down", "up"):
        m = cls == k
        ax.scatter(lfc[m], d["nlp"][m], s=2.2, c=cmap[k], alpha=0.55,
                   linewidths=0, zorder=3 if k != "ns" else 2, label=f"{k} ({int(m.sum()):,})")
    ax.axhline(-np.log10(0.05), color="#333333", lw=0.7, ls="--", zorder=4)
    ax.axvline(0, color="#333333", lw=0.6, ls=":", zorder=4)
    ax.set_xlabel(r"log$_2$ fold change (virus-positive vs. virus-negative)", fontsize=7)
    ax.set_ylabel(r"$-\log_{10}$(adjusted $P$)", fontsize=7)
    ymax = float(np.nanmax(d["nlp"]))
    ax.set_ylim(0, min(ymax * 1.06, 340))
    ax.set_xlim(np.nanpercentile(lfc, 0.05) * 1.1, np.nanpercentile(lfc, 99.95) * 1.1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="upper left", markerscale=3.2)
    ax.set_title("Host transcriptome response (DESeq2)", fontsize=8.5,
                 fontweight="bold", loc="left", pad=6)
    ax.text(1.0, 1.012, f"n = {len(d):,} genes", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.2, color=GREY)
    panel_letter(fig, ax, letter)
    return ax


# --------------------------------------------------------------------------
# 面板 F —— 功能富集 (GO / KEGG dot)
# --------------------------------------------------------------------------
def _enrich_dot(ax, csv, title, ntop=9):
    d = read_table(csv)
    d = d[d["p.adjust"] < 0.05].copy()
    if d.empty:
        ax.axis("off")
        ax.set_title(title, fontsize=7.5, fontweight="bold", loc="left")
        return
    d = d.sort_values("p.adjust").head(ntop).iloc[::-1]
    y = np.arange(len(d))
    gr = d["GeneRatio"].astype(str).apply(lambda s: float(s.split("/")[0]) / float(s.split("/")[1])
                                          if "/" in s else np.nan)
    sc = ax.scatter(gr, y, s=8 + 40 * (d["Count"] / max(d["Count"].max(), 1)),
                    c=-np.log10(d["p.adjust"]), cmap=CMAP_SEQ_R,
                    edgecolor="#333333", linewidth=0.4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([shorten(t, 34) for t in d["Description"]], fontsize=5.8)
    ax.set_xlabel("Gene ratio", fontsize=6.5)
    ax.set_ylim(-0.8, len(d) - 0.2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#DDDDDD", lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=7.5, fontweight="bold", loc="left", pad=4)
    return sc


def panel_F(fig, spec, expr_dir, letter="F", log=print):
    sub = spec.subgridspec(2, 2, hspace=0.80, wspace=0.52)
    axes = [fig.add_subplot(sub[i, j]) for i in range(2) for j in range(2)]
    sets = [("DE_up_GO.csv", "GO \u2014 up-regulated"),
            ("DE_down_GO.csv", "GO \u2014 down-regulated"),
            ("DE_up_KEGG.csv", "KEGG \u2014 up-regulated"),
            ("DE_down_KEGG.csv", "KEGG \u2014 down-regulated")]
    sc = None
    for ax, (f, t) in zip(axes, sets):
        p = os.path.join(expr_dir, f)
        if os.path.exists(p):
            sc = _enrich_dot(ax, p, t) or sc
        else:
            ax.axis("off")
            ax.set_title(t + " (missing)", fontsize=7)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=4.5, mfc="#DDDDDD",
                          mec="#333333", label="Gene count"),
               plt.Line2D([], [], marker="o", ls="", ms=4.5, mfc=UP,
                          mec="#333333", label="More significant")]
    axes[0].legend(handles=handles, frameon=False, loc="lower right",
                   handletextpad=0.2, borderpad=0.2)
    panel_letter(fig, axes[0], letter, dx=-0.22, dy=1.42)
    return axes


# --------------------------------------------------------------------------
# 面板 G —— 病毒 × 组织部位关联（组织嗜性）
# --------------------------------------------------------------------------
ORGANS = ["Fruit", "Leaf", "Flower", "Root", "Stem"]

_ORGAN_ALIAS = {
    "fruit": "Fruit", "berry": "Fruit", "root": "Root",
    "leaf": "Leaf", "leaf of wolfberry": "Leaf",
    "flower": "Flower", "pistil": "Flower", "stamen": "Flower", "pollen": "Flower",
    "stem": "Stem", "stem node": "Stem", "stem tip": "Stem",
    "rhizosphere": "Rhizosphere", "rhizosphere soil": "Rhizosphere", "soil": "Rhizosphere",
    "all plant": "Unresolved",
    # 花+果混合取样，无法归到单一器官
    "mixed sample of flower and fruit": "Other",
}


def norm_organ(t) -> str:
    """把元数据里 41 种自由文本组织名归一到器官类。

    规则（顺序敏感）：
      1) 先剔除匿名/表型/非器官（tissue1..12、gall、all plant）
      2) 精确别名表
      3) 关键词兜底 —— 注意 **先把 "wolfberry" 换成 "lycium"**，
         否则 "leaf of wolfberry" 里的 "berry" 会把它误判成 Fruit。
    """
    if not isinstance(t, str) or not t.strip():
        return "Unknown"
    s = t.strip().lower()
    if re.fullmatch(r"tissue\d+", s):
        return "Unresolved"          # 匿名化的 tissue1..tissue12
    if "gall" in s:
        return "Gall"                # 虫瘿表型，不是器官
    if s in _ORGAN_ALIAS:
        return _ORGAN_ALIAS[s]
    s = s.replace("wolfberry", "lycium")
    if "berry" in s:
        return "Fruit"               # "berry under drought stress in the stage 2" 等自由文本
    if "leaf" in s:                  # "leaf under drought stress in the stage 2"
        return "Leaf"
    if "flower" in s:                # "flower under drought stress in the stage 2"
        return "Flower"
    if s.startswith("stem"):
        return "Stem"
    return "Other"


def cohort_runs(asm_dir: str, log=print) -> set:
    """从权威组装目录解析 445 个 run（条数 = run + logs）。"""
    runs = set()
    for e in os.listdir(asm_dir):
        if os.path.isdir(os.path.join(asm_dir, e)):
            runs |= set(re.findall(r"(CRR\d+|SRR\d+|ERR\d+|DRR\d+)", e))
    return runs


def panel_G(fig, spec, bar_dir, meta_path, asm_dir, letter="G", log=print):
    sub = spec.subgridspec(1, 2, width_ratios=[1.0, 0.034], wspace=0.030)
    ax = fig.add_subplot(sub[0, 0])
    cax = fig.add_subplot(sub[0, 1])

    runs = cohort_runs(asm_dir, log)
    meta = read_table(meta_path, sep="\t", dtype=str, low_memory=False)
    sub = meta[meta["Sample"].isin(runs)].copy()
    sub["organ"] = sub["Tissue"].map(norm_organ)
    organ_runs = {o: set(sub.loc[sub["organ"] == o, "Sample"]) for o in ORGANS}
    all_org = set().union(*organ_runs.values())
    n_org = {o: len(organ_runs[o]) for o in ORGANS}
    total_org = len(all_org)

    up = read_table(os.path.join(bar_dir, "01_detection/summary/all_viruses.best.summary.tsv"), sep="\t")
    det = {}
    for r, v in zip(up["Sample"], up["Adjusted_Species"]):
        det.setdefault(r, set()).add(v)

    counts = up.groupby("Adjusted_Species")["Sample"].nunique().sort_values(ascending=False)
    # 器官子集内检出 >= 3 例才成行（全 0 / 单例的行没有信息量），最多 10 行
    cand = [(v, sum(1 for r in all_org if v in det.get(r, ()))) for v in counts.index]
    keep = [v for v, k in cand if k >= 3][:10]
    log(f"  [G] organ rows: {[(pretty_virus(v, 22), k) for v, k in cand if k >= 3]}")

    prev = np.zeros((len(keep), len(ORGANS)))
    qv = np.ones_like(prev)
    raw_p = []
    for i, v in enumerate(keep):
        tot = sum(1 for r in all_org if v in det.get(r, ()))
        for j, o in enumerate(ORGANS):
            rs = organ_runs[o]
            n, k = len(rs), sum(1 for r in rs if v in det.get(r, ()))
            prev[i, j] = 100.0 * k / n if n else 0.0
            # 该器官 vs 其余器官
            _, p = fisher_exact([[k, n - k], [tot - k, (total_org - n) - (tot - k)]])
            raw_p.append(p)
    pv = multipletests(raw_p, method="fdr_bh")[1].reshape(prev.shape)
    qv = pv

    im = ax.imshow(prev, cmap=CMAP_SEQ_R, vmin=0, vmax=100,
                   extent=(-0.5, len(ORGANS) - 0.5, len(keep) - 0.5, -0.5),
                   aspect="auto")
    ax.set_xlim(-0.5, len(ORGANS) + 1.05)
    ax.set_ylim(len(keep) - 0.5, -0.5)

    def stars(q):
        return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else ""

    for i in range(len(keep)):
        for j in range(len(ORGANS)):
            st = stars(qv[i, j])
            ax.text(j, i, f"{prev[i, j]:.1f}{st}", ha="center", va="center",
                    fontsize=6.0, fontweight="bold" if st else "normal",
                    color="white" if prev[i, j] > 60 else "#222222")
        overall = 100.0 * sum(1 for r in all_org if keep[i] in det.get(r, ())) / total_org
        ax.text(len(ORGANS) + 0.15, i, f"{overall:.1f}", ha="left", va="center",
                fontsize=6.0, color=GREY)

    ax.set_yticks(np.arange(len(keep)))
    ax.set_yticklabels([pretty_virus(v, 28) for v in keep], fontsize=6.2, fontstyle="italic")
    ax.set_xticks(np.arange(len(ORGANS)))
    ax.set_xticklabels([f"{o}\n(n={n_org[o]})" for o in ORGANS], fontsize=6.4)
    ax.text(len(ORGANS) + 0.15, -0.72, "All\norgans", ha="left", va="center",
            fontsize=6.0, color=GREY, fontweight="bold")
    ax.axvline(len(ORGANS) - 0.5, color="#DDDDDD", lw=0.7)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)

    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Detection rate (%)", fontsize=6.2, labelpad=2)
    cb.ax.tick_params(labelsize=5.8, pad=1)
    cb.outline.set_linewidth(0.5)

    # 本面板与 C 同行，单元格较窄 —— 长标题会溢出画布，故副标题下移为 xlabel
    ax.set_title("Tissue tropism of the virome", fontsize=8.5,
                 fontweight="bold", loc="left", pad=6)
    ax.set_xlabel(f"{total_org} organ-annotated runs;  significance: BH-FDR $q$ (Fisher exact)\n"
                  f"* $q$ < 0.05   ** $q$ < 0.01   *** $q$ < 0.001",
                  fontsize=6.4)
    panel_letter(fig, ax, letter, dx=-0.42, dy=1.10)
    log(f"  [G] organ-annotated n={total_org} {n_org}; viruses={len(keep)}; "
        f"significant cells={(qv < 0.05).sum()}")
    return ax


# --------------------------------------------------------------------------
# 面板 E/F（自绘）—— GCVA 变异分析，与 A–D/G 共用同一套 Okabe-Ito 语言
#
# rev.h 的决定：**不再嵌入** `04_post_analysis/` 的成品图。
# 那些图由 snpgenie / vcf_viz 各自的脚本产出（白底灰网格、红点蓝星号、紫色填充、
# 各自的默认字体），与本文 A–D/G 的设计语言完全不同 —— 拼进同一张主图，
# 一眼就是"拼凑"，过不了 SCI 的图件关。数据其实都在（不是只有截图），所以自绘。
#
# 数据来源（都是 03/04 的中间产物，不是截图）：
#   dN/dS  : `_data/gcva_dnds.tsv`  ← scripts/barbarum/extract_snpgenie_dnds.py
#            161 样本 × 6 基因 = 966 行。
#            **product 1..6 -> N/P/P4/M/G/L 的映射取自权威 GTF**
#            `virus-annotations/OR489165.1.gtf`（CDS 436-1950 / 2499-4046 / 4305-5024 /
#            5480-6340 / 6742-7368 / 7562-14344），不是按长度猜的。
#   变异谱 : `04_post_analysis/<GCVA>/vcf_viz/Table1_Unified_HighQuality_Variants.csv`
#            **144 样本 × 44,357 个位点**（含 REF/ALT/ALT_FREQ/Molecular_Type）。
#
# ⚠ 口径必须写进图注（2026-10-10 实测，别再抄错）：
#   * E 用 SNPGenie 侧的 **161** 个样本（位点须在 CDS 内才能定 dN/dS）；
#     F 用 vcf_viz 侧的 **144** 个样本（不要求位点在 CDS 内）。**不是同一批样本**：
#     两侧交集只有 **129** 个，另有 32 个只在 snpgenie、15 个只在 vcf_viz。
#   * **务必用 `known_virus_barbarum/` 下的表**。`known_virus_all_v2/` 里那张同名
#     表是 **7 月旧运行**（175 样本 / 55,768 位点 / Ts/Tv 2.64），曾在 rev.h 期间
#     被误当成本次数据源引用过 —— 数字全对不上，因为它不是同一批样本。
#     本图其余面板（A/B/C/D/G）都读 `known_virus_barbarum`，F 必须一致。
# --------------------------------------------------------------------------
GCVA_ORDER = ["N", "P", "P4", "M", "G", "L"]     # 与 GCVA_GENES 同序 -> 基因配色自动一致
GENE_COLOR = {g: OKABE[i % 7] for i, g in enumerate(GCVA_ORDER)}
NEUTRAL_C = VERM

# 12 类替换（Ts 在前、Tv 在后），键与 REF/ALT 的拼法一致
SUB_TS = ["A>G", "G>A", "C>T", "T>C"]
SUB_TV = ["A>C", "A>T", "C>A", "C>G", "G>C", "G>T", "T>A", "T>G"]


def read_dnds_table(path):
    """读 dN/dS 汇总表（161 样本 × 6 基因）。"""
    df = read_table(path, sep="\t")
    for c in ("mean_dN_vs_ref", "mean_dS_vs_ref", "piN", "piS", "N_sites", "S_sites"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def panel_E_var(fig, spec, dnds, letter="E", log=print):
    """E：选择压力。左 = dN vs dS 散点（双对数 + 中性线 dN=dS）；
    右 = 逐基因 dN/dS（样本散点 + 中位 + IQR + 中性线 1.0）。

    **为什么用 dN-vs-dS 散点而不是 dN/dS 比值的箱线**：
    比值的分母 dS 可以 ≈ 0，会造出 30+ 的假"正选择"点 ——
    04 的 `Fig05_Gene_dNdS_Stats` 就是这么被毁的（y 轴被一个 29.1 的伪影拉到 31，
    其余 5 个基因全被压成一条线）。散点 + 中性线没有这个问题，
    而且把"每个点落在中性线的哪一侧"直接交给读者判断。
    """
    sub = spec.subgridspec(1, 2, width_ratios=[1.36, 1.00], wspace=0.32)
    ax1 = fig.add_subplot(sub[0, 0])
    ax2 = fig.add_subplot(sub[0, 1])

    dN, dS = "mean_dN_vs_ref", "mean_dS_vs_ref"

    # ---- E1: dN vs dS（双对数）------------------------------------------
    # 三类计数要分开，否则图上的标注与正文对不上：
    #   n_def  = dS > 0（比值有定义）—— 这才是"能否算 dN/dS"的人群
    #   n_above= dS > 0 且 dN >= dS  —— 真的偏中性/正选择的对
    #   n_dn0  = dS > 0 但 dN == 0   —— 比值 = 0，但对数轴画不出（dN=0 → -inf）
    #   n_plot = dN > 0 且 dS > 0    —— 真正落在图上的点
    # rev.h 首版只报了 n_plot（809），与正文引用的 917 差 108 个 dN=0 的对，
    # 数字对不上。现在两者都报，图上写 n_def 口径，另用小字注明 off-scale 的对数。
    xlo, xhi = 3.0e-3, 0.34
    ylo, yhi = 3.0e-4, 0.075
    n_def = n_above = n_dn0 = n_plot = 0
    for g in GCVA_ORDER:
        s = dnds[(dnds["gene"] == g) & dnds[dN].notna() & dnds[dS].notna()]
        n_def += int((s[dS] > 0).sum())
        n_above += int(((s[dS] > 0) & (s[dN] >= s[dS])).sum())
        n_dn0 += int(((s[dS] > 0) & (s[dN] == 0)).sum())
        p = s[(s[dN] > 0) & (s[dS] > 0)]
        if not len(p):
            continue
        n_plot += len(p)
        ax1.scatter(p[dS], p[dN], s=5.2, color=GENE_COLOR[g], alpha=0.60,
                    edgecolor="none", zorder=3, label=g)
    ax1.plot([ylo, xhi], [ylo, xhi], color=NEUTRAL_C, lw=1.0, ls="--", zorder=2)
    ax1.text(0.115, 0.052, "neutral\n$d_N=d_S$", ha="right", va="top",
             fontsize=6.0, color=NEUTRAL_C, linespacing=1.2, zorder=6)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(xlo, xhi)
    ax1.set_ylim(ylo, yhi)
    ax1.set_xlabel("$d_S$  synonymous divergence (vs. OR489165.1)", fontsize=6.8)
    ax1.set_ylabel("$d_N$  non-synonymous divergence", fontsize=6.8)
    ax1.legend(loc="lower right", frameon=False, ncol=2, fontsize=5.8, title="gene",
               title_fontsize=5.8, handletextpad=0.28, columnspacing=0.85,
               borderpad=0.15, labelspacing=0.26)
    ax1.text(0.0055, 0.055,
             f"only {n_above} / {n_def} $d_N/d_S$-definable pairs\n"
             f"({100.0 * n_above / max(n_def, 1):.1f} %) have $d_N\\geq d_S$",
             ha="left", va="top", fontsize=6.0, color="#222222", linespacing=1.25, zorder=6)
    ax1.text(0.0055, 0.018,
             f"{n_dn0} pairs with $d_N=0$ below axis range",
             transform=ax1.transAxes, ha="left", va="bottom", fontsize=5.4,
             color=GREY, zorder=6)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.tick_params(which="minor", length=1.3, width=0.5)
    ax1.set_title("Purifying selection dominates the GCVA genome", fontsize=7.3,
                  loc="left", pad=6, color="#333333")

    # ---- E2: 逐基因 dN/dS ------------------------------------------------
    # **y 轴封顶**（ycap）：6 个基因的中位全在 0.44 以下，只有 G 有 14 个点越过 1.0
    # （最大 1.87）。若按 0–1.95 画，整个点云被压进下半区区，N/P/P4/M/L 之间的差异
    # 就看不出来了 —— 而这个差异恰恰是本图要传达的。所以把轴封到 1.12，
    # 越顶的点用**同色空心标记钉在封顶线**上，并在角落注明条数与真实最大值。
    YCAP = 1.12
    n_drop = n_over = 0
    vmax = 0.0
    for i, g in enumerate(GCVA_ORDER):
        s = dnds[(dnds["gene"] == g) & dnds[dN].notna() & dnds[dS].notna()]
        n_drop += int((s[dS] <= 0).sum())
        v = (s.loc[s[dS] > 0, dN] / s.loc[s[dS] > 0, dS]).values
        if not len(v):
            continue
        vmax = max(vmax, float(v.max()))
        jit = (np.random.RandomState(11 + i).rand(len(v)) - 0.5) * 0.62
        ok = v <= YCAP
        ax2.scatter(i + jit[ok], v[ok], s=4.2, color=GENE_COLOR[g], alpha=0.40,
                    edgecolor="none", zorder=3)
        if not ok.all():
            n_over += int((~ok).sum())
            ax2.scatter(i + jit[~ok], np.full(int((~ok).sum()), YCAP), s=7.0,
                        facecolors="none", edgecolors=GENE_COLOR[g], linewidths=0.5,
                        zorder=4)
        med = float(np.median(v))
        q1, q3 = np.percentile(v, [25, 75])
        ax2.plot([i - 0.32, i + 0.32], [med, med], color="#333333", lw=1.4,
                 zorder=5, solid_capstyle="butt")
        ax2.plot([i, i], [q1, q3], color="#333333", lw=0.9, zorder=4)
        # 中位标签贴在横线右侧；但 M 的中位 = 0，直接写在 y=0 会被 x 轴压掉半行，
        # 所以给它一个最小抬升（0.035 ≈ 半个字高）。
        ax2.text(i + 0.40, max(med, 0.035), f"{med:.2f}", ha="left", va="center",
                 fontsize=5.9, color="#333333", fontweight="bold", zorder=6)
    ax2.axhline(1.0, color=NEUTRAL_C, lw=1.0, ls="--", zorder=2)
    ax2.text(5.55, 1.02, "neutral", ha="right", va="bottom", fontsize=6.0,
             color=NEUTRAL_C, zorder=6)
    ax2.set_xticks(range(len(GCVA_ORDER)))
    ax2.set_xticklabels(GCVA_ORDER)
    ax2.set_xlim(-0.60, len(GCVA_ORDER) - 0.12)
    ax2.set_ylim(0.0, YCAP * 1.13)
    ax2.set_ylabel("$d_N/d_S$", fontsize=6.8)
    ax2.set_xlabel("Viral gene (5$'\\rightarrow$3$'$)", fontsize=6.8)
    if n_over:
        # 放在封顶标记的**正上方**（G 列上方是空的）——放在右下角会压到 L 的数据点
        ax2.text(0.615, 0.995,
                 f"$\\circ$ {n_over} pairs above axis limit\n(max $d_N/d_S$ = {vmax:.2f})",
                 transform=ax2.transAxes, ha="left", va="top", fontsize=5.4,
                 color=GREY, linespacing=1.25, zorder=6)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.set_title("Per-gene $d_N/d_S$", fontsize=7.3, loc="left", pad=6, color="#333333")

    panel_letter(fig, ax1, letter, dx=-0.082, dy=1.14)
    log(f"  [E] self-drawn: dN-vs-dS scatter (log-log) — "
        f"{n_above}/{n_def} dN/dS-definable pairs above neutral "
        f"({100.0 * n_above / max(n_def, 1):.1f} %); plotted {n_plot}, "
        f"off-scale dN=0 {n_dn0} | per-gene dN/dS ({n_drop} pairs with dS=0 excluded; "
        f"{n_over} pairs above ycap={YCAP}, max {vmax:.2f})")
    return ax1


def panel_F_var(fig, spec, var_csv, letter="F", log=print):
    """F：宿主内变异谱。左 = 等位频率谱 AFS（低频准种 vs 固定，双峰）；
    右 = 12 类替换谱（Ts 蓝 / Tv 橙，标 Ts/Tv）—— 相当于把 04 的
    `Fig07_Mutational_Spectrum` **剔掉 `A>A/T>T/C>C/G>G` 那四根 self 柱后重画**：
    那四根根本不是替换，却把 y 轴顶到 8000、把真正的 Ts/Tv 压到底部 1/5。"""
    sub = spec.subgridspec(1, 2, width_ratios=[1.30, 1.00], wspace=0.30)
    ax1 = fig.add_subplot(sub[0, 0])
    ax2 = fig.add_subplot(sub[0, 1])

    df = read_table(var_csv)
    af = pd.to_numeric(df["ALT_FREQ"], errors="coerce").dropna()
    n_samp = int(df["Sample_ID"].nunique()) if "Sample_ID" in df.columns else -1

    # ---- F1: AFS ---------------------------------------------------------
    # AF 下限实测 = 0.05（vcf_viz 自己的过滤阈值），bin 从 0.05 起 → 不丢数据、不留空柱。
    # 换物种/换批次后先看日志里的 `off-bin=`：>0 就说明有数据落在直方图外，
    # **必须**把 edges 的起点降到实测 min，否则会静默少画（最隐蔽的一类图错）。
    #
    # 起点 0.05 + 步长 0.02 的一个**刻意的副作用**：bin 边界恰好落在 0.97 上，
    # 于是「末两个 bin（mid 0.98 / 1.00）」的计数之和 == `AF >= 0.97` 的位点数 ——
    # 图上的橙色柱与标注里的 n 严格同口径，不会出现"柱子加起来≠标注数字"。
    # 上界必须 > 1.00（取 1.0201）才能把 **AF 恰好 = 1.0** 的那一批（4,040 个）收进 bin；
    # rev.h 首版取 1.0401，末 bin 是空的 [1.02,1.04)，于是 `cnt[-1] == 0`，
    # `xy=(0.99, cnt[-1]*0.92)` 把箭头钉在 y=0（对数轴上越界）→ **"fixed" 标注整条被静默丢弃**。
    edges = np.arange(0.05, 1.0201, 0.02)
    n_off = int((af.values < edges[0]).sum())
    cnt, _ = np.histogram(af.values, bins=edges)
    mid = edges[:-1] + 0.01
    cols = [GREEN if m < 0.97 else ORANGE for m in mid]
    ax1.bar(mid, np.maximum(cnt, 0.7), width=0.019, color=cols, edgecolor="none", zorder=3)
    ax1.set_yscale("log")
    ax1.set_ylim(0.7, cnt.max() * 2.6)
    ax1.set_xlim(0.045, 1.035)
    ax1.set_xlabel("Intra-host allele frequency (AF)", fontsize=6.8)
    ax1.set_ylabel("Variant sites (log scale)", fontsize=6.8)
    n_low = int((af < 0.5).sum())
    # 与橙色柱**同口径**：`AF >= 0.97` == mid >= 0.97 的那两个 bin 之和（见上面的 bin 边界说明）
    n_fix = int(cnt[mid >= 0.97].sum())
    n_mid = len(af) - n_low - n_fix
    i_low = int(np.argmin(np.abs(mid - 0.20)))
    ax1.annotate(f"low-frequency quasispecies\n"
                 f"AF < 0.5:  n={n_low:,} ({100.0 * n_low / len(af):.0f} %)",
                 xy=(mid[i_low], cnt[i_low]), xytext=(0.15, 0.85),
                 textcoords="axes fraction",
                 fontsize=6.0, color=GREEN, ha="left", va="center", linespacing=1.25,
                 arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.7, shrinkA=1, shrinkB=3))
    ax1.annotate(f"fixed\n"
                 f"AF $\\geq$ 0.97:  n={n_fix:,} ({100.0 * n_fix / len(af):.0f} %)",
                 xy=(mid[-1], cnt[-1] * 0.95), xytext=(0.70, 0.92),
                 textcoords="axes fraction",
                 fontsize=6.0, color=ORANGE, ha="center", va="center", linespacing=1.25,
                 arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.7, shrinkA=1, shrinkB=3))
    # 「中间态」（0.5 <= AF < 0.97）不单独标注：柱体自底部堆到柱顶，坐标内除柱顶上方
    # 没有可用空白，再加一个标签必然压到柱子上。三档占比写进图注即可（日志已打）。
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.set_title("Allele-frequency spectrum: quasispecies + fixed variants",
                  fontsize=7.3, loc="left", pad=6, color="#333333")

    # ---- F2: 12 类替换谱 -------------------------------------------------
    sig = df["REF"].astype(str).str.upper() + ">" + df["ALT"].astype(str).str.upper()
    order = SUB_TS + SUB_TV
    vc = sig.value_counts()
    vals = [int(vc.get(k, 0)) for k in order]
    cols2 = [BLUE if k in SUB_TS else ORANGE for k in order]
    xs = np.arange(len(order))
    ax2.bar(xs, vals, width=0.74, color=cols2, edgecolor="#333333", lw=0.35, zorder=3)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(order, fontsize=5.6, rotation=90)
    ax2.set_ylim(0, max(vals) * 1.28)
    ax2.set_ylabel("Variant sites", fontsize=6.8)
    ax2.set_xlabel("Substitution type", fontsize=6.8)
    ts, tv = sum(vals[:4]), sum(vals[4:])
    n_other = len(df) - ts - tv          # indel / 多核苷酸位点，12 类替换谱里放不下
    ax2.text(0.035, 0.965, f"Ts/Tv = {ts / max(tv, 1):.2f}", transform=ax2.transAxes,
             ha="left", va="top", fontsize=6.6, color="#333333", fontweight="bold")
    if n_other > 0:
        ax2.text(0.035, 0.895, f"{n_other:,} indel sites not shown",
                 transform=ax2.transAxes, ha="left", va="top", fontsize=5.4, color=GREY)
    h = [plt.Line2D([], [], marker="s", ls="", ms=4.2, mfc=BLUE, mec="#333333", mew=0.35),
         plt.Line2D([], [], marker="s", ls="", ms=4.2, mfc=ORANGE, mec="#333333", mew=0.35)]
    ax2.legend(h, ["Transition", "Transversion"], loc="upper right", frameon=False,
               fontsize=5.9, handletextpad=0.3, borderpad=0.15)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.set_title("Substitution spectrum (reference-identical sites excluded)",
                  fontsize=7.3, loc="left", pad=6, color="#333333")

    panel_letter(fig, ax1, letter, dx=-0.082, dy=1.14)
    log(f"  [F] self-drawn: AFS ({len(af):,} sites / {n_samp} samples; "
        f"low-freq AF<0.5 {n_low:,} [{100.0 * n_low / len(af):.0f} %], "
        f"fixed AF>=0.97 {n_fix:,} [{100.0 * n_fix / len(af):.0f} %], "
        f"intermediate {n_mid:,}; off-bin={n_off}) "
        f"+ substitution spectrum (Ts/Tv={ts / max(tv, 1):.2f}, "
        f"{n_other:,} indel sites excluded)")
    return ax1


# ==========================================================================
# 遗传变异分析组合图（自绘 · 9 面板）—— 2026-10-11
# ==========================================================================
# 为什么单独成图：
#   主图 A–D/G 已 20.90 in ≈ 1.4 页（README §五）。把变异分析的 9 个核心面板塞进主图
#   会到 1.7 页以上；而变异分析本身有 9 个「互不重复、各自有结论」的面板。
#   单独成图既保住每面板 ≥ 2.5 in 的可读高度，又让正文能按需裁剪。
#
# 为什么自绘而不是嵌 04 成品图：
#   与 rev.h 同一个理由（README §七）。补充一条**实测**依据 —— 你文档里那 4 张「主图级」
#   候选，按通栏（9.775 in）嵌入后印刷字号为：
#       Fig10_Trinity  源画布 14.8 in -> 4.9 pt  ✗
#       Figure8_PopGen 源画布 15.8 in -> 4.6 pt  ✗
#       maftools_Oncoplot 源画布 14.0 in -> 5.2 pt ✗
#       Fig11a_PCA_2D  源画布  8.9 in -> 8.2 pt  ✓
#   即「整页独立图」的设计塞进组合图后全部低于 6 pt 可读线 —— 只能自绘。
#
# 面板清单（数据源见 --var-data-dir）：
#   A  全基因组变异景观   Figure1A_All_Variants_Landscape_Data.csv   POS×Mean_AF 气泡
#   B  滑窗 π + Tajima's D Table2_Sliding_Window_PopGen_Stats.csv    多样性 + 选择动态
#   C  群体结构 PCA        Matrix_ML_Clusters_3D.csv                 K=3, PC1 88.8%
#   D  准种 πN vs πS       Matrix_ML_Clusters_3D.csv                 净化选择（多态性口径）
#   E  dN vs dS 散点       gcva_dnds.tsv                            净化选择（分歧口径）
#   F  逐基因 dN/dS        gcva_dnds.tsv                            G 基因约束最松
#   G  Top 12 变异热点     Matrix_03_Population_Hotspots.csv         热点定位
#   H  等位频率谱 AFS      Table1_Unified_HighQuality_Variants.csv   准种 vs 固定
#   I  12 类替换谱         Table1_Unified_HighQuality_Variants.csv   Ts/Tv 偏倚

VAR_RICH_TITLE = ("Within-host genetic variation of Cytorhabdovirus sp. 'lycii' "
                  "(GCVA, OR489165.1): genome-wide landscape, population structure, "
                  "selection and quasispecies")

# Fig11a 的 PC1/PC2 方差解释率：**只能从 04 的成品图抄**（坐标表里没有特征值）。
# 换病毒/换批次时务必回 `snpgenie/Fig11a_PCA_2D.png` 的轴标签核对，否则图上数字是错的。
PCA_PC1_PCT, PCA_PC2_PCT = 88.8, 9.6


def _v_track(fig, spec):
    """基因轨道（与 panel_B 的 B0 同画法）。"""
    ax = fig.add_subplot(spec)
    ax.set_xlim(0, GCVA_LEN)
    ax.set_ylim(0, 1)
    ax.axis("off")
    for i, (g, s, e) in enumerate(GCVA_GENES):
        ax.add_patch(Rectangle((s, 0.24), e - s, 0.52, facecolor=OKABE[i % 7],
                               edgecolor="#333333", lw=0.4))
        if (e - s) > 700:
            ax.text((s + e) / 2, 0.50, g, ha="center", va="center",
                    fontsize=6.2, color="white", fontweight="bold")
    ax.plot([0, GCVA_LEN], [0.10, 0.10], color="#333333", lw=0.8)
    ax.text(-0.004, 0.5, "Genome", transform=ax.transAxes, rotation=90,
            ha="center", va="center", fontsize=6.0, color="#333333")
    return ax


def _v_track_split(fig, spec, row_h, gap_in=0.34):
    """把一行拆成『基因轨道 + 主图』两格，返回 (track_ax, plot_ax)。

    ⚠ 版式坑（2026-10-11 实测）：`subgridspec(hspace=...)` 的 hspace 是**平均格高**
    （= row_h/2）的比例，不是绝对英寸。要拿到 `gap_in` 英寸的间隙必须
    `hspace = 2*gap_in/row_h`。原先写 0.08 → 实际只有 0.09 in，而 7.3 pt 标题
    + pad 6 pt 需要 ≈0.19 in → **标题压在轨道基线上**（A/B 两个面板都中招）。
    """
    h = 2.0 * gap_in / row_h if row_h else 0.08
    sub = spec.subgridspec(2, 1, height_ratios=[0.17, 1.0], hspace=h)
    return _v_track(fig, sub[0, 0]), fig.add_subplot(sub[1, 0])


def _v_finish(ax, title, letter=None, tsize=7.3):
    """统一收尾：只留左/下 spine + 面板字母 + 左对齐标题。

    ⚠ 面板字母**必须画在轴框左上角（transAxes 的 (0,1)）而不是用 `panel_letter`
    往左偏**（2026-10-11 实测）：本图是 6 行版式，上一行是**通栏**面板，它的
    xlabel 居中在整幅中点；右列面板的字母若按 dx=-0.135 往左偏，正好落进
    上一行 xlabel 的横向范围里 —— 实测 `B 的 'Genomic coordinate (bp)'` 与
    `D 的字母` 重叠 12×7 px。画在轴框左上角则**必然落在本行的水平带内**，
    同时把标题右移到 x=0.048 给字母让位。
    """
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if letter:
        ax.text(0.0, 1.0, letter, transform=ax.transAxes, fontsize=12,
                fontweight="bold", va="bottom", ha="left", color="#333333")
    if title:
        ax.set_title(title, fontsize=tsize, loc="left",
                     x=0.048 if letter else 0.0, pad=6, color="#333333")


def panel_V_landscape(fig, spec, var_csv, letter="A", log=print, row_h=None):
    """A：全基因组宿主内变异景观（逐位点 AF 气泡；大小 = 携带样本数，色 = Ts/Tv）。"""
    _track, ax = _v_track_split(fig, spec, row_h)
    df = read_table(var_csv)
    for c in ("POS", "Mean_AF", "Sample_Count"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["POS", "Mean_AF", "Sample_Count"])
    is_ts = df["Molecular_Type"].astype(str).str.contains("Transition")
    smax = float(df["Sample_Count"].max()) or 1.0
    size = 1.6 + 30.0 * (df["Sample_Count"] / smax) ** 1.6
    for m, c, lab in ((is_ts, BLUE, "Transition (Ts)"),
                      (~is_ts, ORANGE, "Transversion (Tv)")):
        ax.scatter(df.loc[m, "POS"], df.loc[m, "Mean_AF"], s=size[m], color=c,
                   alpha=0.50, edgecolor="none", zorder=3, label=lab)
    ax.set_xlim(0, GCVA_LEN)
    ax.set_ylim(-0.03, 1.06)
    ax.set_xlabel("Genomic coordinate (bp)", fontsize=6.8)
    ax.set_ylabel("Mean intra-host allele frequency", fontsize=6.8)
    n_sites = int(df["POS"].nunique())
    n_samp = int(smax)
    _v_finish(ax, f"Genome-wide intra-host variant landscape "
                  f"({n_sites:,} sites, {len(df):,} sample\u00d7site records; "
                  f"bubble size = no. of carrying samples)", letter)
    hs = [plt.Line2D([], [], marker="o", ls="",
                     ms=float(np.sqrt(1.6 + 30.0 * (v / smax) ** 1.6)),
                     mfc="#999999", mec="none")
          for v in (smax * 0.25, smax * 0.55, smax)]
    # 这一格是**满幅密集散点**（4,428 条记录铺满 0~1 的 AF 区间），
    # 按全项目惯例的 frameon=False 会让图例文字压在气泡上读不出来 —— 例外加半透明底。
    _lgkw = dict(frameon=True, framealpha=0.86, facecolor="white",
                 edgecolor="none", fancybox=False)
    lg1 = ax.legend(loc="upper right", fontsize=5.9, ncol=2,
                    handletextpad=0.3, columnspacing=0.9, title="substitution",
                    title_fontsize=5.9, **_lgkw)
    ax.add_artist(lg1)
    ax.legend(hs, [f"{int(smax * 0.25)}", f"{int(smax * 0.55)}", f"{n_samp}"],
              loc="upper left", fontsize=5.6, ncol=3,
              handletextpad=0.4, columnspacing=0.8, title="samples",
              title_fontsize=5.9, **_lgkw)
    log(f"  [A] landscape: {n_sites:,} sites / {len(df):,} records, max carriers={n_samp}")
    return ax


def panel_V_pi_tajima(fig, spec, win_csv, letter="B", log=print, row_h=None):
    """B：滑窗核苷酸多样性 π + Tajima's D（590 bp 窗）。

    ⚠ 口径陷阱（2026-10-11 实测）：`Table2_Sliding_Window_PopGen_Stats.csv` 的 `PI` 列
    **不是 per-site π**，而是窗口内**各位点 π 的求和** —— 见
    `virome_analysis_pipeline/virus_variants_analyzer.py:248`
    `pi_window = in_window['Pi'].sum()`（对照：同仓库 `virus_vcf_pipeline.py:1402`
    的实现是 `/ valid_sites` 求均值，**两套实现口径不同**，本表列名 `N_VARIANTS_S`
    表明它来自 `virus_variants_analyzer.py`）。所以原始值是 2.0–20.4 这种量级。
    本面板**除以窗口长度换成 per-site π**（= 对全部位点求均值，未分离位点 π=0），
    这样才与其他图的 π 可比、也才是文献里报的 π；Tajima's D 仍用求和尺度（其分子
    `π − S/a1` 必须与 S 同尺度），故两者 y 轴量级不同是**正确**的。
    换算只在作图时做，CSV 不动；日志会同时打两种口径。
    """
    _track, ax = _v_track_split(fig, spec, row_h)
    w = read_table(win_csv)
    bs = pd.to_numeric(w["BIN_START"], errors="coerce").values
    be = pd.to_numeric(w["BIN_END"], errors="coerce").values
    mid = pd.to_numeric(w["BIN_MID"], errors="coerce").values
    bw = float(np.nanmedian(be - bs)) + 1.0          # 窗口长度（bp），含两端
    # 末窗可能越过基因组末端 -> 分母按实际覆盖长度截断，避免低估
    Lw = np.minimum(be, GCVA_LEN) - bs + 1.0
    Lw = np.where(np.isfinite(Lw) & (Lw > 0), Lw, bw)
    pi_sum = pd.to_numeric(w["PI"], errors="coerce").values
    pi = pi_sum / Lw                                  # -> per-site π
    td = pd.to_numeric(w["TAJIMA_D"], errors="coerce").values
    ax.bar(mid, pi, width=bw * 0.90, color=SKY, edgecolor="none", alpha=0.45, zorder=2)
    ax.plot(mid, pi, color=BLUE, lw=1.0, zorder=3, label="Nucleotide diversity $\\pi$")
    for _, s, _e in GCVA_GENES:
        ax.axvline(s, color="#DDDDDD", lw=0.5, zorder=1)
    ax.set_xlim(0, GCVA_LEN)
    ax.set_ylim(0, float(np.nanmax(pi)) * 1.22)
    ax.set_ylabel("Nucleotide diversity $\\pi$ (per site)", fontsize=6.8, color=BLUE)
    ax.tick_params(axis="y", labelcolor=BLUE)
    ax.set_xlabel("Genomic coordinate (bp)", fontsize=6.8)
    ax2 = ax.twinx()
    ax2.plot(mid, td, color=VERM, lw=1.1, zorder=4, label="Tajima's $D$")
    ax2.axhline(0.0, color="#999999", lw=0.8, ls="--", zorder=1)
    ax2.set_ylabel("Tajima's $D$", fontsize=6.8, color=VERM)
    ax2.tick_params(axis="y", labelcolor=VERM)
    lo, hi = float(np.nanmin(td)), float(np.nanmax(td))
    ax2.set_ylim(min(lo * 1.18, -0.4), max(hi * 1.6, 0.4))
    # ⚠ 标注要贴在 `D=0` 那条线上，不能用固定轴分数 0.03 —— 右轴 ylim 是
    #   (-3.29, 0.4)，`D=0` 在**顶部**（轴分数 0.89），标在 0.03 会飘到底部、
    #   与线脱节。这里按 ylim 反算轴分数，并把文字放在线上方（那片是空的）。
    _d0f = (0.0 - ax2.get_ylim()[0]) / (ax2.get_ylim()[1] - ax2.get_ylim()[0])
    ax2.text(0.997, min(_d0f + 0.02, 0.995), "neutrality  $D=0$",
             transform=ax2.transAxes, ha="right", va="bottom",
             fontsize=5.7, color="#999999")
    ax2.spines["top"].set_visible(False)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", frameon=False, ncol=2, fontsize=5.9)
    n_neg = int(np.sum(td < 0))
    _v_finish(ax, f"Sliding-window diversity and selection "
                  f"({len(mid)} windows of {bw:.0f} bp; Tajima's $D<0$ in "
                  f"{n_neg}/{int(np.isfinite(td).sum())} windows)", letter)
    log(f"  [B] pi+Tajima: {len(mid)} windows x {bw:.0f} bp; "
        f"per-site pi [{np.nanmin(pi):.4f}, {np.nanmax(pi):.4f}] "
        f"(原始求和列 PI [{np.nanmin(pi_sum):.2f}, {np.nanmax(pi_sum):.2f}] / 窗口长度); "
        f"Tajima D range [{lo:.2f}, {hi:.2f}], D<0 {n_neg}/{int(np.isfinite(td).sum())}")
    return ax


def panel_V_pca(fig, spec, ml_csv, letter="C", log=print):
    """C：群体结构 PCA（K=3 无监督聚类，PC1 88.8%）。"""
    ax = fig.add_subplot(spec)
    d = read_table(ml_csv)
    x = pd.to_numeric(d["PC1"], errors="coerce").values
    y = pd.to_numeric(d["PC2"], errors="coerce").values
    cl = d["ML_Cluster"].astype(str).values
    classes = sorted(set(cl), key=lambda s: (len(s), s))
    n_out = 0
    for i, c in enumerate(classes):
        m = cl == c
        ax.scatter(x[m], y[m], s=13, color=OKABE[i % 7], alpha=0.80,
                   edgecolor="#333333", linewidths=0.25, zorder=3, label=c)
        if m.sum() >= 3:                      # 协方差椭圆（<3 点无意义）
            cov = np.cov(np.vstack([x[m], y[m]]))
            vals, vecs = np.linalg.eigh(cov)
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:, order]
            ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
            # **1 SD（68 %）而不是 95 %**：实测 Evo-Class 3 有个 (PC1≈2, PC2≈3.2) 的
            # 离群样本，95 % 椭圆（×2.447）会被它撑成横跨整个象限的大色块，
            # 反而读不出簇的形状；1 SD 稳健得多且是 PCA 图的通行画法。
            wid, hei = 2.0 * 1.0 * np.sqrt(np.maximum(vals, 0))
            ax.add_patch(Ellipse((x[m].mean(), y[m].mean()), wid, hei, angle=ang,
                                 facecolor=OKABE[i % 7], alpha=0.12,
                                 edgecolor=OKABE[i % 7], lw=0.7, zorder=1))
            # 落在椭圆外的点数（按马氏距离）—— 让读者知道有离群样本
            inv = np.linalg.inv(cov)
            d2 = np.array([(np.array([xi - x[m].mean(), yi - y[m].mean()]) @ inv
                            @ np.array([xi - x[m].mean(), yi - y[m].mean()]))
                           for xi, yi in zip(x[m], y[m])])
            n_out += int((d2 > 1.0 ** 2 * 2).sum())
    ax.set_xlabel(f"PC1 ({PCA_PC1_PCT:.1f} %)", fontsize=6.8)
    ax.set_ylabel(f"PC2 ({PCA_PC2_PCT:.1f} %)", fontsize=6.8)
    ax.legend(loc="upper right", frameon=False, fontsize=5.9, title="cluster",
              title_fontsize=5.9, handletextpad=0.3, borderpad=0.2)
    _v_finish(ax, f"Population structure: PCA of {len(d)} samples "
                  f"(unsupervised K=3; 1 SD ellipse)", letter)
    log(f"  [C] PCA: {len(d)} samples, " +
        ", ".join(f"{c}={int((cl == c).sum())}" for c in classes) +
        f"; {n_out} points outside 1 SD ellipse")
    return ax


def panel_V_pin_pis(fig, spec, ml_csv, letter="D", log=print):
    """D：准种多样性 πN vs πS（多态性口径的净化选择证据）。"""
    ax = fig.add_subplot(spec)
    d = read_table(ml_csv)
    xs = pd.to_numeric(d["piS"], errors="coerce")
    ys = pd.to_numeric(d["piN"], errors="coerce")
    ok = xs.notna() & ys.notna() & (xs > 0) & (ys > 0)
    x, y = xs[ok].values, ys[ok].values
    ax.scatter(x, y, s=9, color=PURPLE, alpha=0.55, edgecolor="none", zorder=3)
    lim = [min(x.min(), y.min()) * 0.55, max(x.max(), y.max()) * 1.8]
    ax.plot(lim, lim, color=NEUTRAL_C, lw=1.0, ls="--", zorder=2)
    ax.text(lim[1] * 0.62, lim[1] * 0.30, "neutral\n$\\pi_N=\\pi_S$", ha="right", va="top",
            fontsize=5.9, color=NEUTRAL_C, linespacing=1.2, zorder=6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("$\\pi_S$  synonymous diversity", fontsize=6.8)
    ax.set_ylabel("$\\pi_N$  non-synonymous diversity", fontsize=6.8)
    n_above = int((y >= x).sum())
    ratio = y / x
    ax.text(0.035, 0.965,
            f"$\\pi_N\\geq\\pi_S$ in {n_above} / {len(x)} samples "
            f"({100.0 * n_above / max(len(x), 1):.1f} %)\n"
            f"$\\pi_N/\\pi_S$ = {np.min(ratio):.3f}\u2013{np.max(ratio):.3f} "
            f"(median {np.median(ratio):.3f})",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.0,
            color="#222222", linespacing=1.25, zorder=6)
    _v_finish(ax, "Purifying selection at the polymorphism level "
                  "($\\pi_N$ vs $\\pi_S$)", letter)
    log(f"  [D] piN-vs-piS: {len(x)} samples, {n_above} above neutral "
        f"({100.0 * n_above / max(len(x), 1):.1f} %); "
        f"piN/piS {np.min(ratio):.3f}-{np.max(ratio):.3f} (median {np.median(ratio):.3f})")
    return ax


def _draw_dnds_scatter(ax, dnds):
    """E 左：dN vs dS 双对数散点 + 中性线。返回 (n_def, n_above, n_dn0, n_plot)。"""
    dN, dS = "mean_dN_vs_ref", "mean_dS_vs_ref"
    xlo, xhi = 3.0e-3, 0.34
    ylo, yhi = 3.0e-4, 0.075
    n_def = n_above = n_dn0 = n_plot = 0
    for g in GCVA_ORDER:
        s = dnds[(dnds["gene"] == g) & dnds[dN].notna() & dnds[dS].notna()]
        n_def += int((s[dS] > 0).sum())
        n_above += int(((s[dS] > 0) & (s[dN] >= s[dS])).sum())
        n_dn0 += int(((s[dS] > 0) & (s[dN] == 0)).sum())
        p = s[(s[dN] > 0) & (s[dS] > 0)]
        if not len(p):
            continue
        n_plot += len(p)
        ax.scatter(p[dS], p[dN], s=5.2, color=GENE_COLOR[g], alpha=0.60,
                   edgecolor="none", zorder=3, label=g)
    ax.plot([ylo, xhi], [ylo, xhi], color=NEUTRAL_C, lw=1.0, ls="--", zorder=2)
    ax.text(0.115, 0.052, "neutral\n$d_N=d_S$", ha="right", va="top",
            fontsize=6.0, color=NEUTRAL_C, linespacing=1.2, zorder=6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel("$d_S$  synonymous divergence (vs. OR489165.1)", fontsize=6.8)
    ax.set_ylabel("$d_N$  non-synonymous divergence", fontsize=6.8)
    ax.legend(loc="lower right", frameon=False, ncol=2, fontsize=5.8, title="gene",
              title_fontsize=5.8, handletextpad=0.28, columnspacing=0.85,
              borderpad=0.15, labelspacing=0.26)
    # ⚠ 原先把 "n pairs with dN=0" 单独画在**轴分数** (0.0055, 0.018) —— 实测压在左下角的
    #   数据点上（那里正好有一簇 dS≈0.008–0.02、dN≈3e-4 的点）；而上面那条说明块用的是
    #   **数据坐标** (0.0055, 0.055)，所以跑到了左上角。两处坐标口径不一致本身就容易踩坑，
    #   统一改成轴分数并把 "dN=0" 并入说明块：左上角（中性线之上）是这张图唯一稳定留白的区域。
    ax.text(0.02, 0.93,
            f"only {n_above} / {n_def} $d_N/d_S$-definable pairs\n"
            f"({100.0 * n_above / max(n_def, 1):.1f} %) have $d_N\\geq d_S$\n"
            f"{n_dn0} pairs with $d_N=0$ not shown",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.0,
            color="#222222", linespacing=1.25, zorder=6)
    ax.tick_params(which="minor", length=1.3, width=0.5)
    return n_def, n_above, n_dn0, n_plot


def _draw_dnds_gene(ax, dnds, ycap=1.12):
    """E 右：逐基因 dN/dS（抖动散点 + 中位 + IQR + 中性线）。"""
    dN, dS = "mean_dN_vs_ref", "mean_dS_vs_ref"
    n_drop = n_over = 0
    vmax = 0.0
    for i, g in enumerate(GCVA_ORDER):
        s = dnds[(dnds["gene"] == g) & dnds[dN].notna() & dnds[dS].notna()]
        n_drop += int((s[dS] <= 0).sum())
        v = (s.loc[s[dS] > 0, dN] / s.loc[s[dS] > 0, dS]).values
        if not len(v):
            continue
        vmax = max(vmax, float(v.max()))
        jit = (np.random.RandomState(11 + i).rand(len(v)) - 0.5) * 0.62
        ok = v <= ycap
        ax.scatter(i + jit[ok], v[ok], s=4.2, color=GENE_COLOR[g], alpha=0.40,
                   edgecolor="none", zorder=3)
        if not ok.all():
            n_over += int((~ok).sum())
            ax.scatter(i + jit[~ok], np.full(int((~ok).sum()), ycap), s=7.0,
                       facecolors="none", edgecolors=GENE_COLOR[g], linewidths=0.5, zorder=4)
        med = float(np.median(v))
        q1, q3 = np.percentile(v, [25, 75])
        ax.plot([i - 0.32, i + 0.32], [med, med], color="#333333", lw=1.4,
                zorder=5, solid_capstyle="butt")
        ax.plot([i, i], [q1, q3], color="#333333", lw=0.9, zorder=4)
        ax.text(i + 0.40, max(med, 0.035), f"{med:.2f}", ha="left", va="center",
                fontsize=5.9, color="#333333", fontweight="bold", zorder=6)
    ax.axhline(1.0, color=NEUTRAL_C, lw=1.0, ls="--", zorder=2)
    ax.text(5.55, 1.02, "neutral", ha="right", va="bottom", fontsize=6.0,
            color=NEUTRAL_C, zorder=6)
    ax.set_xticks(range(len(GCVA_ORDER)))
    ax.set_xticklabels(GCVA_ORDER)
    ax.set_xlim(-0.60, len(GCVA_ORDER) - 0.12)
    ax.set_ylim(0.0, ycap * 1.13)
    ax.set_ylabel("$d_N/d_S$", fontsize=6.8)
    ax.set_xlabel("Viral gene (5$'\\rightarrow$3$'$)", fontsize=6.8)
    if n_over:
        ax.text(0.615, 0.995,
                f"$\\circ$ {n_over} pairs above axis limit\n(max $d_N/d_S$ = {vmax:.2f})",
                transform=ax.transAxes, ha="left", va="top", fontsize=5.4,
                color=GREY, linespacing=1.25, zorder=6)
    return n_drop, n_over, vmax


def panel_V_dnds_scatter(fig, spec, dnds, letter="E", log=print):
    """E：dN vs dS 散点（分歧口径的净化选择）。"""
    ax = fig.add_subplot(spec)
    n_def, n_above, n_dn0, n_plot = _draw_dnds_scatter(ax, dnds)
    _v_finish(ax, "Purifying selection dominates the GCVA genome", letter)
    log(f"  [E] dN-vs-dS: {n_above}/{n_def} pairs above neutral "
        f"({100.0 * n_above / max(n_def, 1):.1f} %); plotted {n_plot}, dN=0 {n_dn0}")
    return ax


def panel_V_dnds_gene(fig, spec, dnds, letter="F", log=print):
    """F：逐基因 dN/dS（y 轴封顶 1.12 + 越顶空心圈）。"""
    ax = fig.add_subplot(spec)
    n_drop, n_over, vmax = _draw_dnds_gene(ax, dnds)
    _v_finish(ax, "Per-gene $d_N/dS$", letter)
    log(f"  [F] per-gene dN/dS: {n_drop} pairs with dS=0 excluded; "
        f"{n_over} above cap, max={vmax:.2f}")
    return ax


_AA3TO1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
           "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
           "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
           "Tyr": "Y", "Val": "V", "Ter": "*", "Stop": "*"}


def _short_hotspot(v):
    """`P_3646_p.Ala383Val` -> `P 3646 A383V`（y 刻度标签必须短，否则挤出画布左缘）。

    实测：原样 19 字符 @5.5pt ≈ 1.2 in，比左边距（0.86 in）还宽 →
    整个 tight bbox 被撑宽、版面右移、字母 G 跑到画布外。
    """
    p = str(v).split("_")
    if len(p) >= 3:
        g, pos, rest = p[0], p[1], "_".join(p[2:])
        rest = rest.replace("p.", "").replace("c.", "")
        m = re.match(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3})$", rest)
        if m and m.group(1) in _AA3TO1 and m.group(3) in _AA3TO1:
            rest = f"{_AA3TO1[m.group(1)]}{m.group(2)}{_AA3TO1[m.group(3)]}"
        return f"{g} {pos} {rest}"
    return str(v)


def panel_V_hotspots(fig, spec, hs_csv, letter="G", log=print, n_top=12):
    """G：Top N 变异热点（按携带样本数）。"""
    ax = fig.add_subplot(spec)
    d = read_table(hs_csv)
    d["Sample_Count"] = pd.to_numeric(d["Sample_Count"], errors="coerce")
    d = d.dropna(subset=["Sample_Count"]).sort_values("Sample_Count", ascending=False)
    d = d.head(n_top).iloc[::-1]
    ypos = np.arange(len(d))
    genes = [str(g) for g in d["GENE"]]
    cols = [GENE_COLOR.get(g, GREY) for g in genes]
    ax.barh(ypos, d["Sample_Count"], height=0.72, color=cols,
            edgecolor="#333333", lw=0.35, zorder=3)
    for yi, (_, r) in zip(ypos, d.iterrows()):
        ax.text(r["Sample_Count"] + d["Sample_Count"].max() * 0.015, yi,
                f"{int(r['Sample_Count'])}", va="center", ha="left",
                fontsize=5.7, color="#222222")
    ax.set_yticks(ypos)
    ax.set_yticklabels([_short_hotspot(v) for v in d["Variant_ID"]], fontsize=5.5)
    # 基因名已写在 y 刻度标签里 -> 直接把标签染成基因色，省掉图例。
    # （原先的图例放在 lower right，与最短那根柱的数值标签重叠 5×6 px。）
    for lab, g in zip(ax.get_yticklabels(), genes):
        lab.set_color(GENE_COLOR.get(g, "#222222"))
    ax.set_xlabel("No. of carrying samples", fontsize=6.8)
    ax.set_xlim(0, d["Sample_Count"].max() * 1.16)
    _v_finish(ax, f"Top {len(d)} intra-host variant hotspots", letter)
    log(f"  [G] hotspots: top {len(d)} by sample count, "
        f"max={int(d['Sample_Count'].max())}; "
        f"genes={dict(zip(*np.unique(genes, return_counts=True)))}")
    return ax


def _draw_afs(ax, var_csv, ytop=1.9):
    """F 左：AFS 直方图（低频准种 vs 固定双峰）。返回 (n_low, n_fix, n_mid, n_off, n_af, n_samp)。

    `ytop` = y 上限相对最高柱的倍数。注意这是**对数轴**：倍数要留够标注空间，
    但也不能太大（原先 2.6 → 图上 60 % 的高度是空的）。
    ⚠ 本函数目前只服务 `_VariationRich` 的 H 面板；`panel_F_var` 里还有一份**同源的
    内联副本**（主图已定稿，暂不动它，避免改出回归）。两者口径必须保持一致。
    """
    df = read_table(var_csv)
    af = pd.to_numeric(df["ALT_FREQ"], errors="coerce").dropna()
    n_samp = int(df["Sample_ID"].nunique()) if "Sample_ID" in df.columns else -1
    edges = np.arange(0.05, 1.0201, 0.02)
    n_off = int((af.values < edges[0]).sum())
    cnt, _ = np.histogram(af.values, bins=edges)
    mid = edges[:-1] + 0.01
    cols = [GREEN if m < 0.97 else ORANGE for m in mid]
    ax.bar(mid, np.maximum(cnt, 0.7), width=0.019, color=cols, edgecolor="none", zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.7, cnt.max() * ytop)
    ax.set_xlim(0.045, 1.035)
    ax.set_xlabel("Intra-host allele frequency (AF)", fontsize=6.8)
    ax.set_ylabel("Variant sites (log scale)", fontsize=6.8)
    n_low = int((af < 0.5).sum())
    n_fix = int(cnt[mid >= 0.97].sum())
    n_mid = len(af) - n_low - n_fix
    i_low = int(np.argmin(np.abs(mid - 0.20)))
    # 两条标注都拉到顶部空白区：AF≈1.0 那根橙柱很高，右对齐放会撞上它。
    ax.annotate(f"low-frequency quasispecies\nAF < 0.5:  n={n_low:,} "
                f"({100.0 * n_low / len(af):.0f} %)",
                xy=(mid[i_low], cnt[i_low]), xytext=(0.05, 0.93),
                textcoords="axes fraction", fontsize=6.0, color=GREEN,
                ha="left", va="top", linespacing=1.25,
                arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.7, shrinkA=1, shrinkB=3))
    ax.annotate(f"fixed\nAF $\\geq$ 0.97:  n={n_fix:,} "
                f"({100.0 * n_fix / len(af):.0f} %)",
                xy=(mid[-1], cnt[-1] * 0.95), xytext=(0.63, 0.93),
                textcoords="axes fraction", fontsize=6.0, color=ORANGE,
                ha="center", va="top", linespacing=1.25,
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.7, shrinkA=1, shrinkB=3))
    return n_low, n_fix, n_mid, n_off, len(af), n_samp


def panel_V_afs(fig, spec, var_csv, letter="H", log=print):
    """H：等位频率谱 AFS（准种 vs 固定）。"""
    ax = fig.add_subplot(spec)
    n_low, n_fix, n_mid, n_off, n_af, n_samp = _draw_afs(ax, var_csv)
    _v_finish(ax, "Allele-frequency spectrum: quasispecies + fixed variants", letter)
    log(f"  [H] AFS: {n_af:,} sites / {n_samp} samples; low {n_low:,} "
        f"[{100.0 * n_low / n_af:.0f} %], fixed {n_fix:,} "
        f"[{100.0 * n_fix / n_af:.0f} %], mid {n_mid:,}; off-bin={n_off}")
    return ax


def _draw_spectrum(ax, var_csv):
    """F 右：12 类替换谱。返回 (ts, tv, n_other)。"""
    df = read_table(var_csv)
    sig = df["REF"].astype(str).str.upper() + ">" + df["ALT"].astype(str).str.upper()
    order = SUB_TS + SUB_TV
    vc = sig.value_counts()
    vals = [int(vc.get(k, 0)) for k in order]
    cols2 = [BLUE if k in SUB_TS else ORANGE for k in order]
    xs = np.arange(len(order))
    ax.bar(xs, vals, width=0.74, color=cols2, edgecolor="#333333", lw=0.35, zorder=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(order, fontsize=5.6, rotation=90)
    ax.set_ylim(0, max(vals) * 1.28)
    ax.set_ylabel("Variant sites", fontsize=6.8)
    ax.set_xlabel("Substitution type", fontsize=6.8)
    ts, tv = sum(vals[:4]), sum(vals[4:])
    n_other = len(df) - ts - tv
    ax.text(0.035, 0.965, f"Ts/Tv = {ts / max(tv, 1):.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=6.6, color="#333333", fontweight="bold")
    if n_other > 0:
        ax.text(0.035, 0.895, f"{n_other:,} indel sites not shown",
                transform=ax.transAxes, ha="left", va="top", fontsize=5.4, color=GREY)
    h = [plt.Line2D([], [], marker="s", ls="", ms=4.2, mfc=BLUE, mec="#333333", mew=0.35),
         plt.Line2D([], [], marker="s", ls="", ms=4.2, mfc=ORANGE, mec="#333333", mew=0.35)]
    ax.legend(h, ["Transition", "Transversion"], loc="upper right", frameon=False,
              fontsize=5.9, handletextpad=0.3, borderpad=0.15)
    return ts, tv, n_other


def panel_V_spectrum(fig, spec, var_csv, letter="I", log=print):
    """I：12 类替换谱（Ts/Tv）。"""
    ax = fig.add_subplot(spec)
    ts, tv, n_other = _draw_spectrum(ax, var_csv)
    _v_finish(ax, "Substitution spectrum (reference-identical sites excluded)", letter)
    log(f"  [I] spectrum: Ts/Tv={ts / max(tv, 1):.2f}, {n_other:,} indel sites excluded")
    return ax


VAR_TABLES = {
    "var_csv":  "Table1_Unified_HighQuality_Variants.csv",
    "dnds_csv": "gcva_dnds.tsv",
    "land_csv": "Figure1A_All_Variants_Landscape_Data.csv",
    "win_csv":  "Table2_Sliding_Window_PopGen_Stats.csv",
    "ml_csv":   "Matrix_ML_Clusters_3D.csv",
    "hs_csv":   "Matrix_03_Population_Hotspots.csv",
}


def build_var_rich(data_dir, out_dir, stem, log=print):
    """生成自绘的「遗传变异分析」9 面板组合图（PNG + 矢量 PDF）。

    `data_dir` 可以是**一个目录**，也可以是**目录列表**（按序查找）——
    6 张表分散在 `04_post_analysis/<GCVA>/` 的三个子目录
    （`vcf_viz/` 3 张、`snpgenie/` 1 张、`snpeff_macro/` 1 张）+ `09_report/.../_data/` 1 张，
    所以服务器上直接用目录列表，不必先拼一个扁平目录。
    """
    dirs = [data_dir] if isinstance(data_dir, str) else list(data_dir)
    paths = {}
    for key, name in VAR_TABLES.items():
        hit = next((os.path.join(d, name) for d in dirs
                    if d and os.path.exists(os.path.join(d, name))), None)
        if not hit:
            raise FileNotFoundError(
                f"variation figure 缺数据表 {name}；已找过:\n  " +
                "\n  ".join(str(d) for d in dirs))
        paths[key] = hit
        log(f"  [var-rich] {name}  <-  {hit}")
    var_csv, dnds_csv = paths["var_csv"], paths["dnds_csv"]
    land_csv, win_csv = paths["land_csv"], paths["win_csv"]
    ml_csv, hs_csv = paths["ml_csv"], paths["hs_csv"]
    dnds = read_dnds_table(dnds_csv)

    # 行高（in）：两个通栏 + 三行双栏 + 一个通栏尾行
    H1, H2, HR, H3 = 2.25, 2.45, 3.05, 1.95
    # 行间距必须**装得下上一行的 x 刻度标签 + xlabel + 本行的面板标题**：
    #   刻度标签 6.8pt(0.094) + tick pad(0.049) + labelpad(0.056) + xlabel 6.8pt(0.094)
    #   = 0.293 in（上一行轴框以下）
    #   标题 pad 6pt(0.083) + 7.3pt 字高(0.101) = 0.184 in（本行轴框以上）
    #   => 0.477 in + 视觉余量 → 取 0.55 in。**再小就会标题压 xlabel**（实测 0.40 时中招）。
    hsp = 0.55
    rows = [H1, H2, HR, HR, HR, H3]
    H = sum(rows) + hsp * (len(rows) - 1)
    fh = 0.62 + H + 0.34
    # ⚠ hspace 是"平均行高"的比例，不是绝对英寸。想让每道间隙正好 = hsp 英寸，
    #   需要 hspace = hsp * n / Σrows（推导：cell_h = (tot − hsp(n−1))/n = Σrows/n）。
    #   原先写成 hsp/(H/(n−1)) = hsp(n−1)/Σrows → 间隙只有 0.31 in 而非 0.40 in，
    #   于是**每一行的面板标题都压到上一行的 xlabel 上**（A/B 的标题还压住基因轨道）。
    hspace = hsp * len(rows) / sum(rows)
    log(f"[layout] var-rich 9 panels: rows={rows} in, gap={hsp} in "
        f"(hspace={hspace:.4f}), H={H:.2f} in, fig_h={fh:.2f} in")
    fig = plt.figure(figsize=(11.5, fh))
    gs = fig.add_gridspec(
        len(rows), 2, left=0.075, right=0.955,
        top=1.0 - 0.62 / fh, bottom=0.34 / fh,
        height_ratios=rows, hspace=hspace,
        width_ratios=[1.0, 1.0], wspace=0.30)

    panel_V_landscape(fig, gs[0, :], land_csv, "A", log, row_h=H1)
    panel_V_pi_tajima(fig, gs[1, :], win_csv, "B", log, row_h=H2)
    panel_V_pca(fig, gs[2, 0], ml_csv, "C", log)
    panel_V_pin_pis(fig, gs[2, 1], ml_csv, "D", log)
    panel_V_dnds_scatter(fig, gs[3, 0], dnds, "E", log)
    panel_V_dnds_gene(fig, gs[3, 1], dnds, "F", log)
    panel_V_hotspots(fig, gs[4, 0], hs_csv, "G", log)
    panel_V_afs(fig, gs[4, 1], var_csv, "H", log)
    panel_V_spectrum(fig, gs[5, :], var_csv, "I", log)

    fig.text(0.5, 1.0 - 0.20 / fh, VAR_RICH_TITLE, ha="center", va="bottom",
             fontsize=10, fontweight="bold", color="#333333")

    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, f"{stem}_VariationRich.png")
    pdf = os.path.join(out_dir, f"{stem}_VariationRich.pdf")
    for p in (png, pdf):
        if os.path.islink(p):
            os.unlink(p)
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    from PIL import Image as _Img
    with _Img.open(png) as _im:
        log(f"[size] variation-rich {_im.size[0]} x {_im.size[1]} px @300dpi")
    log(f"[ok ] {png}")
    log(f"[ok ] {pdf}")
    return png


# --------------------------------------------------------------------------
# 主程序
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="生成 barbarum 已知病毒 SCI 正文组合主图")
    ap.add_argument("--root", default="/home/zhangwenda/MMPV-paper/goji-virome")
    ap.add_argument("--out", default=None)
    ap.add_argument("--denominator", type=int, default=DENOM_DEFAULT)
    ap.add_argument("--stem", default="Fig_KnownVirome_Composite")
    ap.add_argument("--metadata", default=None,
                    help="全局元数据 TSV（含 Tissue 列），默认 01_data/sra_rna.data4/global_metadata/")
    ap.add_argument("--assembly-root", default=None,
                    help="权威组装目录（分母来源），默认 02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly")
    ap.add_argument("--ef", choices=("self", "variation", "deg", "none"), default="self",
                    help="E/F 两个通栏面板放什么。"
                         "self（默认，rev.h）= **自绘** GCVA 变异分析（dN-vs-dS 散点 + 逐基因 dN/dS | "
                         "AFS + 替换谱），与 A–D/G 共用同一套 Okabe-Ito 语言；"
                         "variation = 嵌入 04_post_analysis 的成品图（**外观与主图不统一，不推荐用于投稿**）；"
                         "deg = 自绘宿主 DEG 火山图 + GO/KEGG 富集（rev.d 行为）；"
                         "none（rev.j）= **去掉 E/F 两行**，主图只剩 A/D/B/C/G 五个面板（≈1.0 页），"
                         "变异内容改由 --var-rich 的 _VariationRich 承担（无重复）")
    ap.add_argument("--b-last", action="store_true",
                    help="（rev.k，仅与 `--ef none` 连用）把「群体基因组景观」面板移到 C/G 两栏"
                         "**下方**，并把面板字母按**读序**重排为 A–E："
                         "A 检出漏斗 → B 序列同一性 → C 共感染 → D 组织嗜性 → E 群体变异景观"
                         "（即原 A→A / D→B / C→C / G→D / B→E）。"
                         "⚠ 正文里的「面板 X」引用必须同步改。不加此旗标则保持 rev.j 的 "
                         "A→D→B→C→G 旧读序（字母不动）。")
    ap.add_argument("--append-figs", nargs="?", const=CHOSEN_FIGS, default=None,
                    help="在主图**下方**追加面板（EF_LIBRARY key，或 .png 绝对路径）。"
                         "语法：`;` 分隔**行**、`,` 分隔**同一行内的面板**；不含 `;` 时"
                         "每个逗号项各占一个通栏行（旧行为）。两个**自绘 token**（不读位图）："
                         "**`E`** = GCVA 群体基因组景观（三层；仅 `--ef none` 有效，"
                         "且只能作行首）；**`L`** = 全基因组突变枢纽棒棒糖图（y = 单位点去重"
                         "携带样本数）。面板字母按**读序**顺延（自绘 token 也占字母）。"
                         f"不带值 = 用户 2026-10-12 定的拼图：`{CHOSEN_FIGS}`"
                         "（素材固定在 `reports/sci_composite/_append_src/`，走 `{repo}/` 前缀）。"
                         "传空串可显式关闭。"
                         "⚠ 追加图按**原图**嵌入（不重绘、不降采样），源图分辨率极高的那些"
                         "按 300 dpi 折算的印刷字号会偏小 —— 只如实记日志，不做优化。")
    ap.add_argument("--dnds-table", default=None,
                    help="自绘 E 用的 dN/dS 汇总表，默认 <out>/_data/gcva_dnds.tsv；"
                         "由 scripts/barbarum/extract_snpgenie_dnds.py 生成")
    ap.add_argument("--var-table", default=None,
                    help="自绘 F 用的变异总表，默认 <04_post_analysis>/<GCVA>/vcf_viz/"
                         "Table1_Unified_HighQuality_Variants.csv")
    ap.add_argument("--panel-e", default="boot_dnds",
                    help=f"E 用哪张 04 图（key 或 .png 路径）。可用: {sorted(EF_LIBRARY)}")
    ap.add_argument("--panel-f", default="afs",
                    help="F 用哪张 04 图（key 或 .png 路径），同上")
    ap.add_argument("--post-dir", default=None,
                    help="04_post_analysis/<病毒>/ 目录，默认 barbarum 的 "
                         "Cytorhabdovirus_sp._lycii_OR489165.1")
    ap.add_argument("--host-figure", action="store_true",
                    help="额外把宿主响应（火山图 + GO/KEGG）单独出一张图，"
                         "这样 E/F 让给变异分析时不丢内容")
    ap.add_argument("--var-figure", action="store_true",
                    help="额外出一张 `_Variation` 变异分析组合图（默认 4 个通栏面板，"
                         "供补充材料；正文只能容 1-2 张通栏，见 README §九）")
    ap.add_argument("--var-panels", default="dnds_gene,pca_2d,hotspots,af_violin",
                    help="`--var-figure` 的 4 个面板，逗号分隔的 EF_LIBRARY key"
                         "（也可给 .png 绝对路径）。默认 = 补充材料那四张：逐基因 dN/dS 分布 / "
                         "群体 PCA / 变异热点 / 逐样本 AF 小提琴。"
                         "**选面板前先看日志的『印刷字号』**：源图是 4800 px 那批"
                         "（intrahost_div / manhattan / var_density / lollipop 等）缩到通栏只有 4.6–5.3 pt，不合格。")
    ap.add_argument("--var-rich", action="store_true",
                    help="另出自绘的 **9 面板**遗传变异分析组合图 "
                         "`<stem>_VariationRich.png/.pdf`（A 全基因组变异景观 / B 滑窗 π+Tajima's D / "
                         "C 群体 PCA(K=3) / D πN-vs-πS / E dN-vs-dS 逐位点 / F 逐基因 dN-dS / "
                         "G 变异热点 / H 等位基因频率谱 / I 替换谱）。"
                         "与 `--var-figure` 的区别：那个是**嵌 04 成品图**（外观不统一、只能通栏），"
                         "这个是**自绘**（与主图 A–D/G 同一套 Okabe-Ito + 只留左/下 spine 语言）。"
                         "数据表默认取 `<out>/_diag/data/`。")
    ap.add_argument("--var-data-dir", default=None,
                    help="`--var-rich` 的数据表目录（可给多个，用 `%s` 分隔）。"
                         "默认按序找：<out>/_diag/data、<out>/_data、"
                         "<out>/../_diag/data，以及 04_post_analysis/<GCVA>/ 下的 "
                         "vcf_viz / snpgenie / snpeff_macro 三个子目录。"
                         % os.pathsep)
    args = ap.parse_args()

    root = args.root
    bar_dir = os.path.join(root, "03_known_virus/known_virus_barbarum")
    gcva_dir = os.path.join(bar_dir, "04_post_analysis/Cytorhabdovirus_sp._lycii_OR489165.1")
    expr_dir = os.path.join(root, "05_expression/coexpression/barbarum/back_up")
    meta_path = args.metadata or os.path.join(
        root, "01_data/sra_rna.data4/global_metadata/Global_Metadata_Final_with_Coordinates.tsv")
    asm_dir = args.assembly_root or os.path.join(
        root, "02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly")
    out_dir = args.out or os.path.join(bar_dir, "09_report/S9_SCI_Composite")
    os.makedirs(out_dir, exist_ok=True)

    def log(*a):
        print(*a, flush=True)

    log(f"[root] {root}")
    log(f"[out ] {out_dir}")

    # ---- 版面用「英寸」算，再由总高换算成 figure 分数 -----------------------
    # A 面板高度按行数自适应，目标行距 A_PITCH_PT：
    #   pitch_pt = 72 * H_ax / (n + 2.0)   （+2.0 = ylim 留白：下 0.75 + 上 1.25）
    #   => H_ax = (n + 2.0) * A_PITCH_PT / 72
    # 标题与 xlabel / 图例都画在轴框之外，由上方 0.60 in 与下方 A_B_GAP 的留白承接。
    try:
        funnel, _fmeta = load_funnel_table(bar_dir, log=log)
        n_a = len(funnel)
    except Exception as exc:                       # 源表缺失不致命，A 留空即可
        log(f"  [warn] 漏斗表读取失败（面板 A 将留空）: {exc}")
        funnel, n_a = None, 26
    gsA_h = (n_a + 2.0) * A_PITCH_PT / 72.0
    # REST_H = 图题 0.60 + A/B 间隔 0.90 + B~F 区 + 下边距 0.53
    #   A/B 间隔 0.90 要装下 A 轴框下方整摞：x 刻度标签 0.15 in + 图例（1 行）0.17 in
    #   + 自适应 labelpad 后的 xlabel 0.10 in ≈ 0.52-0.56 in，再给 B 的面板标题留 ~0.20 in。
    #   实测余量 ~0.30 in；若把 A_PITCH_PT 或图例段数改大，必须回看 panel_A 打的
    #   "legend anchor / labelpad" 日志与 fig_h，别只看总高没变就以为没动版。
    A_B_GAP = 0.90
    # ---- B~F 区：四行的高度（in）显式给出，再换算成 height_ratios ----------------
    # 约定：height_ratios 直接填"英寸"，配 bottom/top 反算出等单位 u，于是
    #   H_bf = Σd × (1 + hspace×(n-1)/n)
    # （hspace 是"平均行高"的比例，n=4 → ×3/4）。rev.d 实测 d=[2.414,1.971,2.052,2.856]、
    # hspace=0.46 → H_bf=12.499 ≈ 12.50，与旧写死的常数一致，说明该换算是对的。
    d_b, d_cg = 2.414, 1.971
    post_dir = args.post_dir or gcva_dir

    # ---- 先解析 `--append-figs` 的**行结构**：决定 `E` 是否被挪进拼图行 -----------
    # 必须在算 `ratio_bf` / `H_bf` 之前做 —— E 一旦挪走，主 gridspec 就只剩 C|D 一行。
    _spec_raw = (args.append_figs or "").strip()
    PACKED = ";" in _spec_raw                       # 含 `;` = 拼图模式（一行多张）
    if _spec_raw:
        if PACKED:
            rows_spec = [[t.strip() for t in _r.split(",") if t.strip()]
                         for _r in _spec_raw.split(";")]
        else:
            rows_spec = [[t.strip()] for t in _spec_raw.split(",") if t.strip()]
        rows_spec = [_r for _r in rows_spec if _r]
    else:
        rows_spec = []
    # 两个**自绘 token**（不读位图，走 `panel_E_land` / `panel_L` 现画）：
    #   `E` = GCVA 群体基因组景观（**三层** = π+D 双轴 + 逐基因 Miss/Syn + 基因条带）
    #   `L` = 全基因组突变枢纽棒棒糖图（重绘自 PPT 原图）
    # 两者都**只在 `--ef none` 的拼图档**下用（字母表按 A–I 重排，见下）。
    _SELF_TOKENS = ("E", "L")
    E_MOVED = PACKED and any("E" in _r for _r in rows_spec)
    if E_MOVED and args.ef != "none":
        raise SystemExit("[err] `--append-figs` 里的 `E` 只在 `--ef none` 下有意义"
                         "（其它档的 E 是 04 成品图 / 自绘 dN-dS，不叫这个名）")
    if E_MOVED and any(_r.count("E") > 1 or "E" in _r[1:] for _r in rows_spec):
        raise SystemExit("[err] `E` 只能作为拼图行里**每行第一个**位置出现一次")
    for _r in rows_spec:
        if _r.count("L") > 1:
            raise SystemExit("[err] 每个拼图行里 `L`（自绘棒棒糖）最多出现一次")
    # dN/dS 汇总表：**只给 `--ef self` 的自绘 E/F 用**（rev.o 起 `E` token 不再需要它）。
    # 默认找 `<out>/_data/`；找不到再退到**正式产物目录**的同名文件 —— 这样
    # `--out /tmp/...` 试版也能直接出图，不必每次先拷 90 KB 的表。
    if args.dnds_table:
        dnds_path = args.dnds_table
    else:
        _cand_dnds = [os.path.join(out_dir, "_data/gcva_dnds.tsv"),
                      os.path.join(bar_dir, "09_report/S9_SCI_Composite/_data/gcva_dnds.tsv")]
        dnds_path = next((_p for _p in _cand_dnds if os.path.exists(_p)), _cand_dnds[0])
    e_img = f_img = None
    dnds_df = var_csv = None
    if args.ef == "variation":
        # E/F 的高度必须由嵌入图的**实际宽高比**决定，否则 imshow(aspect=equal)
        # 会在格子里留信件框白边（纵向白白多出几英寸）。
        try:
            e_png, e_cap = resolve_ef("E", args.panel_e, post_dir, log)
            f_png, f_cap = resolve_ef("F", args.panel_f, post_dir, log)
            e_img, f_img = _crop_white(e_png), _crop_white(f_png)
            EF_W = (0.945 - 0.095) * 11.5            # 通栏可用宽度 = 9.775 in
            d_e = EF_W * e_img.height / e_img.width
            d_f = EF_W * f_img.height / f_img.width
        except Exception as exc:
            log(f"  [warn] 变异面板解析失败，回退到自绘 DEG 版: {exc}")
            args.ef, e_img, f_img = "deg", None, None
    if args.ef == "self":
        # 自绘：行高由我们定（双面板，宽高比 ≈ 2.06），不吃嵌入图的宽高比。
        # `dnds_path` 已在上面解析好（含正式目录回退），这里不重复赋值。
        var_csv = args.var_table or os.path.join(
            post_dir, "vcf_viz/Table1_Unified_HighQuality_Variants.csv")
        if not os.path.exists(dnds_path):
            raise SystemExit(
                f"[err] 缺少 dN/dS 汇总表 {dnds_path}\n"
                f"      先用 scripts/barbarum/extract_snpgenie_dnds.py 生成，"
                f"或改回 --ef deg / --ef variation")
        if not os.path.exists(var_csv):
            raise SystemExit(f"[err] 缺少变异总表 {var_csv}")
        dnds_df = read_dnds_table(dnds_path)
        d_e = d_f = 4.75
        log(f"  [self] dN/dS 表 {os.path.basename(dnds_path)} ({len(dnds_df)} 行, "
            f"{dnds_df['sample'].nunique()} 样本 × {dnds_df['gene'].nunique()} 基因)")
        # E 与 F 的样本集**不同**（E 要求位点在 CDS 内），图注必须写明。
        # 这里把两条口径的重叠实测出来打进日志，避免正文/图注写错样本数。
        try:
            _v = read_table(var_csv)
            _s_e = set(dnds_df["sample"].astype(str))
            _s_f = set(_v["Sample_ID"].astype(str))
            log(f"  [self] F 表 {os.path.basename(var_csv)} ({len(_v):,} 位点, "
                f"{len(_s_f)} 样本) | E∩F 样本 {len(_s_e & _s_f)} "
                f"(E 独有 {len(_s_e - _s_f)} / F 独有 {len(_s_f - _s_e)})")
        except Exception as exc:
            log(f"  [warn] F 表预检失败（不影响出图）: {exc}")
    elif args.ef == "deg":
        d_e, d_f = 2.052, 2.856                      # rev.d 的自绘几何
    else:                                            # --ef none：只留 B / C|G 两行
        d_e = d_f = None
    if args.b_last and args.ef != "none":
        raise SystemExit("[err] `--b-last` 目前只支持与 `--ef none` 连用"
                         "（E/F 两行还在时「B 放到底部」的读序未定义）")
    # ---- 面板字母（rev.k：--b-last 时按**读序**重排为 A–E）-------------------
    # 不加 --b-last 则沿用 v3/rev.j 的旧读序 A→D→B→C→G（字母不动），
    # 以免正文里的「面板 X」引用错位。
    if args.b_last:
        L_A, L_D, L_C, L_G, L_B = "A", "B", "C", "D", "E"
    else:
        L_A, L_D, L_C, L_G, L_B = "A", "D", "C", "G", "B"
    if args.ef == "none":
        # --b-last：row0 = C|G，row1 = B（通栏）；否则 row0 = B（通栏），row1 = C|G
        # E_MOVED（rev.n 拼图）：E 挪进下方拼图行 → 主 gridspec 只剩 C|D 一行
        if E_MOVED:
            ratio_bf = [d_cg]
        else:
            ratio_bf = [d_cg, d_b] if args.b_last else [d_b, d_cg]
    else:
        ratio_bf = [d_b, d_cg, d_e, d_f]
    n_bf = len(ratio_bf)
    mean_bf = sum(ratio_bf) / n_bf
    # 行间距固定成"英寸"更稳：hspace 是按平均行高成比例的，行一变高间距就会失控地涨
    hspace = 0.46 if args.ef == "deg" else min(0.46, 0.95 / mean_bf)
    H_bf = sum(ratio_bf) * (1 + hspace * (n_bf - 1.0) / n_bf)
    if args.ef == "none":
        _ef_geo = "CG/B  (rev.k, --b-last)" if args.b_last else "B/CG  (rev.j)"
    else:
        _ef_geo = f"B/CG/E/F = {d_b:.2f}/{d_cg:.2f}/{d_e:.2f}/{d_f:.2f}"
    log(f"[layout] ef={args.ef}{'  --b-last' if args.b_last else ''}  {_ef_geo} in  "
        f"hspace={hspace:.3f}  H_bf={H_bf:.2f} in")
    log(f"[layout] 面板字母 A={L_A}  D={L_D}  C={L_C}  G={L_G}  B={L_B}")
    # ---- 追加 / 拼图行（`--append-figs`）---------------------------------------
    # 用户明确要求：这几张**按原图**进主图，不要重绘。行高由**裁白边后**的宽高比决定
    # （否则 `imshow(aspect="equal")` 会在格子里留信件框白边），所以必须在建 figure
    # **之前**就裁好、量好 —— 与 `--ef variation` 同一套做法。
    # rev.n：`;` 让**一行能放多张**（用户要求「E|G|H 一行、F|附图 一行」）。
    # rev.o：自绘 token 由 `E` 一个扩到 `E`/`L` 两个，字母改为**按读序顺延**（见下）。
    EF_W = (0.945 - 0.095) * 11.5                   # 通栏可用宽度 = 9.775 in
    EXTRA_GAP = 0.55                                # **拼图行之间**的英寸间距（紧凑）
    # ↑ rev.t（用户 2026-10-12「全图的尺寸可以再调整下嘛」→ 选「紧凑去留白」）：0.95 → **0.55**。
    #   0.95 是 rev.n 为「E 是自绘、右侧有 ylabel/ytick」留的富余，但那时两行之间夹着的是
    #   **整块留白**；现在两行都是「面板 + 图注」，0.55 in 已足够（体检判据 0 对，见 §10.19）。
    GAP_MAIN = 0.80                                 # **主区块 ↔ 拼图块**的英寸间距（更大）
    # ↑ rev.t 踩坑后新增：`EXTRA_GAP` 原先**同时**管两种边界，但二者需要的余量不同 ——
    #   ① 拼图行**之间**：上行的位图铺满轴框、下行的自绘面板字母只外挂 0.045×行高，
    #      0.55 in 够；
    #   ② **主区块 ↔ 拼图块**：主区块最下一行（`--b-last` 时是 C|G）的 **xlabel 挂在轴框
    #      下方**（组织嗜性面板是**两行** xlabel：`319 organ-annotated runs; …` +
    #      `* q<0.05 ** q<0.01 *** q<0.001`），拼图块最上一行（位图行）的**图注又挂在
    #      轴框上方**（`panel_var` 的 `set_title(..., pad=5)`）—— **两者都在"声明行高"之外**，
    #      实测两侧墨迹合计 67 px，0.55 in = 55 px 时**重叠 12 px**（体检报 1 对）。
    #      0.80 in = 80 px → 余量 13 px。这一处比行内缝宽 0.25 in ≈ 75 px。
    PACK_GAP = 0.70                                 # **同一行内**面板之间的水平缝（英寸）
    # ↑ 0.30 太窄：E 是**自绘**的、右侧有 ylabel + ytick（`Nucleotide diversity π`、
    #   `0.03`…），会捅进右邻列的图注/字母里。实测 `_check_main_overlap.py` 在 0.30 时
    #   报 2 处新碰撞（ytick `0.03`↔字母 `F`、ylabel↔F 的标题），0.55 后归零。
    # `PACK_E_H` / `PACK_L_H` 是模块级常量（见 `CHOSEN_ROWS` 下方）—— `_nat_h` 要读。
    extra_rows = []                                 # [[(L,key,path,cap,im,kind)], ...]
    # 字母表：主 gridspec 里已经用掉的字母之后的第一个空位。`E_MOVED` 时 `panel_B`
    # 不再画在主网格里 → 字母 `E` 空出来给合并景观面板（这正是它的名字）。
    _used = {L_A, L_D, L_C, L_G}
    if not E_MOVED:
        _used.add(L_B)
    if args.ef != "none":
        _used |= {"E", "F"}
    _next = max(ord(_c) for _c in _used) + 1

    # ---- 列宽：**等高解**（rev.t）--------------------------------------------------
    # 位图面板的高度 = 列宽 × 高/宽，而行高取行内 max → 行内较矮的那张必然在上下留出死白
    # （rev.s 实测 E 列：mafsummary 自然高 3.26 in 挤在 4.05 in 的行里 → 上下各 ~0.40 in 白）。
    # 解法：让**同一行内各列的显示高度相等** —— 列宽按各图的**宽高比 (w/h)** 分配，
    # 于是每列高度都 = `可用宽 / Σ(w/h)`，行内零死白、行高从 4.05 降到 ~3.62 in。
    # ⚠ 这一组列宽**所有拼图行共用** —— 否则上下两行之间的列缝会错位，并排看很显眼。
    _col_fracs = None
    if PACKED:
        _nset = {len(_r) for _r in rows_spec}
        if len(_nset) == 1 and next(iter(_nset)) > 1:
            for _row in rows_spec:
                if any(_k in _SELF_TOKENS for _k in _row):
                    continue                     # 自绘面板高度固定 → 不参与等高解
                _asp = []
                for _k in _row:
                    try:
                        _pp, _ = resolve_ef("probe", _k, post_dir,
                                            lambda *a, **k: None, root=root)
                        _ii = (_crop_white(_pp, drop_title=False)
                               if _k in KEEP_SRC_TITLE
                               else _crop_white(_pp, title_frac=0.05))
                    except Exception:
                        _asp = []
                        break
                    _asp.append(_ii.width / float(_ii.height))
                if len(_asp) == len(_row) and all(a > 0 for a in _asp):
                    _s = float(sum(_asp))
                    _col_fracs = [a / _s for a in _asp]
                break
    if _col_fracs:
        log(f"[layout] 拼图列宽按**等高解**分配（宽高比 "
            + " / ".join(f"{a:.3f}" for a in _col_fracs)
            + f"）→ 各列显示高相等 = {EF_W - PACK_GAP * (len(_col_fracs) - 1):.2f} / "
            f"Σ(w/h) in")

    for _ri, _row in enumerate(rows_spec):
        _n = len(_row)
        _avail = EF_W if not PACKED else EF_W - PACK_GAP * (_n - 1)
        _fracs = (_col_fracs if (_col_fracs is not None and len(_col_fracs) == _n)
                  else [1.0 / _n] * _n)
        _ws = [_avail * _f for _f in _fracs]
        _items = []
        for _k, _w in zip(_row, _ws):
            _L = chr(_next)
            _next += 1
            if _k in _SELF_TOKENS:                   # `E` / `L`：自绘，不读位图
                _items.append((_L, _k, None, None, None, _k))
                continue
            try:
                _p, _cap = resolve_ef(f"append{_ri}_{_L}", _k, post_dir, log, root=root)
                # `title_frac=0.05`：追加图多为**高分辨率 PPT 导出**（宽 3000–18432 px），
                # 源图自带标题的带宽远超 `_crop_white` 默认的绝对阈值 100 px
                # （实测 oncoplot 181 px / lollipop 241 px）→ 不传比例阈值的话标题会
                # 原样留在面板上，与我们的图注重复。
                # 例外：`KEEP_SRC_TITLE` 里的图**首段墨迹就是真正的子图标题**，不能切。
                _im = (_crop_white(_p, drop_title=False)
                       if _k in KEEP_SRC_TITLE else _crop_white(_p, title_frac=0.05))
            except Exception as exc:
                raise SystemExit(f"[err] `--append-figs` 解析 {_k!r} 失败: {exc}")
            _items.append((_L, _k, _p, _cap, _im, None))
        # 行高 = 行内各面板自然高的最大值（自绘面板用固定高度 PACK_E_H / PACK_L_H）
        _row_h = max(_nat_h(_kind, _im, _w)
                     for (_, _, _, _, _im, _kind), _w in zip(_items, _ws))
        extra_rows.append((_items, _ws, _row_h))
        for (_L, _k, _p, _cap, _im, _kind), _w in zip(_items, _ws):
            if _kind == "E":
                log(f"  [{_L}] 自绘「GCVA 群体基因组景观」（三层，"
                    f"列宽 {_w:.2f} in, 高 {PACK_E_H:.2f} in）")
                continue
            if _kind == "L":
                log(f"  [{_L}] 自绘「全基因组突变枢纽棒棒糖图」（列宽 {_w:.2f} in, "
                    f"高 {PACK_L_H:.2f} in）")
                continue
            _nh = _w * _im.height / _im.width
            _pt = 12.0 * _w * 300.0 / _im.width * (7.09 / 11.5)
            log(f"  [{_L}] {_k}: {os.path.basename(_p)} ({_im.width}x{_im.height} px, "
                f"列宽 {_w:.2f} in → 自然高 {_nh:.2f} in, "
                f"印刷字号 ≈ {_pt:.2f} pt@180mm"
                f"{'  OK' if _pt >= 6 else '  ** 偏小 **'})")
    # 缝隙清单（自上而下）：第 0 条是**主区块 ↔ 拼图块**（`GAP_MAIN`，更宽），
    # 其余是**拼图行之间**（`EXTRA_GAP`）。见 `GAP_MAIN` 的注释。
    _gaps = [GAP_MAIN] + [EXTRA_GAP] * max(0, len(extra_rows) - 1)
    H_extra = sum(_rh for _, _, _rh in extra_rows) + sum(_gaps)
    if extra_rows:
        log(f"[layout] 追加 {len(extra_rows)} 行 / "
            f"{sum(len(_it) for _it, _, _ in extra_rows)} 个面板  H_extra={H_extra:.2f} in"
            f"（含 主↔拼图 {GAP_MAIN} in × 1"
            f" + 行间 {EXTRA_GAP} in × {max(0, len(extra_rows) - 1)}"
            + (f" + 行内缝 {PACK_GAP} in × {sum(max(0, len(_it)-1) for _it, _, _ in extra_rows)}"
               if PACKED else "") + "）")

    REST_H = 0.60 + A_B_GAP + H_bf + H_extra + 0.53
    fig_h = gsA_h + REST_H

    fig = plt.figure(figsize=(11.5, fig_h))
    top = 1.0 - 0.60 / fig_h
    a_bot = top - gsA_h / fig_h
    gs_top = a_bot - A_B_GAP / fig_h

    # ---- 顶行 = A | D（两个 gridspec：A 要宽左边距放物种名，D 只要窄栏） ----
    # A 左边距 0.155 fig ≈ 1.78 in；panel_A 会打印 y 刻度标签的实测最宽值
    # （barbarum 这 9 名 prettify 后 7.0 pt 斜体约 1.6 in）与左端 x 坐标，
    # 换物种集/换字号后以那个日志为准复核，不要沿用这里的数字。
    # 右侧标注列 1.46 in、表头 1.5 in → xlim=1.58×xmax 下 A 轴宽需 ≥ 4.4 in（现 5.75 in，余量足）。
    gsA = fig.add_gridspec(1, 1, left=0.155, right=0.655, top=top, bottom=a_bot)
    gsD = fig.add_gridspec(1, 1, left=0.710, right=0.945, top=top, bottom=a_bot)
    # B..F —— n_bf 行 x 2 列（n_bf = 4，除 `--ef none` 时 = 2）：
    #   ef=none 且**不加** --b-last（rev.j）：row0 = B 变异景观(通栏) | row1 = C 共感染 | G 组织嗜性
    #   ef=none 且**加** --b-last  （rev.k）：row0 = C 共感染 | G 组织嗜性 | row1 = B 变异景观(通栏)
    #   其它 ef：row0 = B 变异景观(通栏)      row1 = C 共感染 | G 组织嗜性
    #            row2 = E 变异分析 1(通栏)    row3 = F 变异分析 2(通栏)   ← ef=variation 时是 04 成品图
    #            （ef=deg 时 row2/3 = 火山图 / 功能富集）
    # height_ratios 直接填"英寸"（见上面的换算说明）；hspace 在变异模式下按"目标行距 0.95 in"
    # 反算，否则嵌入图那两行一变高、比例式的 hspace 会把行距拉到 1.8 in 以上。
    # wspace=0.62 是为了让第 2 列的 G 面板（长行标签 + 色标）不撞到第 1 列。
    gs = fig.add_gridspec(n_bf, 2, left=0.095, right=0.945, top=gs_top,
                          bottom=(0.53 + H_extra) / fig_h,
                          height_ratios=ratio_bf,
                          width_ratios=[1.0, 1.0],
                          hspace=hspace, wspace=0.62)
    log(f"[layout] fig_h={fig_h:.2f}in  gsA={gsA_h:.2f}in "
        f"(rows={n_a}, pitch={A_PITCH_PT}pt)  top row = A[0.155-0.655] | D[0.710-0.945]")

    panel_A(fig, gsA[0, 0], funnel, args.denominator, letter=L_A, log=log)
    panel_D(fig, gsD[0, 0], bar_dir, letter=L_D, log=log)
    if args.b_last:
        # rev.k：row0 = C 共感染 | G 组织嗜性，row1 = B 群体基因组景观（通栏）
        # rev.n（E_MOVED）：E 已被挪进下方拼图行 → 这里只画 C|D，不画 E
        panel_C(fig, gs[0, 0], bar_dir, args.denominator, letter=L_C, log=log)
        panel_G(fig, gs[0, 1], bar_dir, meta_path, asm_dir, letter=L_G, log=log)
        if not E_MOVED:
            panel_B(fig, gs[1, :], gcva_dir, letter=L_B, log=log)
    else:
        panel_B(fig, gs[0, :], gcva_dir, letter=L_B, log=log)
        panel_C(fig, gs[1, 0], bar_dir, args.denominator, letter=L_C, log=log)
        panel_G(fig, gs[1, 1], bar_dir, meta_path, asm_dir, letter=L_G, log=log)
    if args.ef == "variation":
        panel_var(fig, gs[2, :], e_png, e_cap, "E", log=log)
        panel_var(fig, gs[3, :], f_png, f_cap, "F", log=log)
    elif args.ef == "self":
        panel_E_var(fig, gs[2, :], dnds_df, log=log)
        panel_F_var(fig, gs[3, :], var_csv, log=log)
    elif args.ef == "deg":
        panel_E(fig, gs[2, :], expr_dir, log=log)
        panel_F(fig, gs[3, :], expr_dir, log=log)
    # --ef none：E/F 两行整体不存在（变异内容由 --var-rich 的 _VariationRich 承担）

    # ---- 追加 / 拼图行：从下往上摆，最下一行贴 0.53 in 下边距 ------------------
    # 每行内部：列宽 `_ws`（rev.t 起按**等高解**分配，见上；无解时退回等宽），
    # 行内水平缝 `PACK_GAP`；各面板**保持自身宽高比**，比行高矮的在行内**垂直居中**
    # （`_b` 偏移），不拉伸、不留信件框。等高解生效时 `_h == _row_h`，居中偏移自然为 0。
    # ---- 拼图块的左边距：**与主网格左缘对齐**（rev.t，两件事）-----------------------
    # ① 坐标系坑（用户 2026-10-12「你都没对齐呀」的真因）：
    #    `gs` 的 `left=0.095` 是 **figure 分数**（= 0.095 × 11.5 = **1.0925 in**），
    #    而 `_PACK_L` 是**英寸**。早期代码两处都写成字面量 `0.095`，于是拼图块比主网格
    #    整整左移 `1.0925 − 0.095 = 0.998 in`；rev.t 又把 `_PACK_L` 提到 0.362 in，
    #    错位仍有 `1.0925 − 0.362 = 0.73 in`（= 219 px）—— 上下两块**同宽（9.775 in）
    #    却左右错开**，右缘也跟着短 0.73 in。现在统一换算成英寸再比。
    # ② 字母留位：`panel_var` / `panel_E_land` / `panel_L` 都用 `dx = -0.062`（轴分数）
    #    把字母画在**轴框左外 6.2 % 列宽**处；列宽 ~5 in 时是 ~0.31 in。若左边距不够，
    #    `bbox_inches="tight"` 会把画布向左撑大、字母孤零零挂在整幅左缘外（rev.t 实测
    #    tightbbox x0 = −0.217 in，整幅从 3289 撑到 3429）。所以取**两者的大者**。
    _GRID_L = 0.095 * 11.5                            # 主网格左缘（英寸）= 1.0925 in
    _max_col_w = max((_w for _, _ws, _ in extra_rows for _w in _ws), default=0.0)
    _PACK_L = _GRID_L
    if PACKED:
        _PACK_L = max(_PACK_L, 0.062 * _max_col_w + 0.05)
    if extra_rows:
        log(f"[layout] 拼图块左边距 = **{_PACK_L:.3f} in**（主网格左缘 {_GRID_L:.3f} in"
            + (f"，最宽列 {_max_col_w:.2f} in 的字母外挂量 {0.062 * _max_col_w:.3f} in"
               if PACKED else "") + "）")
        _ltr = _PACK_L - 0.062 * _max_col_w if PACKED else _PACK_L
        log(f"[layout] 拼图块字母 x = {_ltr:.3f} in，右缘 = {_PACK_L + EF_W:.3f} in"
            f"（主网格 {_GRID_L:.3f} – {_GRID_L + EF_W:.3f} in）")

    _bot = 0.53
    for _gi, (_items, _ws, _row_h) in enumerate(reversed(extra_rows)):
        _x = _PACK_L
        for (_L, _k, _p, _cap, _im, _kind), _w in zip(_items, _ws):
            _lf = _x / 11.5
            _rf = (_x + _w) / 11.5
            _h = _nat_h(_kind, _im, _w)
            _b = _bot + (_row_h - _h) / 2.0
            _g = fig.add_gridspec(1, 1, left=_lf, right=_rf,
                                  bottom=_b / fig_h, top=(_b + _h) / fig_h)
            # rev.t：本格的整格几何（figure 分数）—— 传给自绘面板，让它们的**字母与图注
            # 锚到格左缘**而不是内缩后的轴框左缘，从而与位图行的 `panel_var`（轴框 = 整格）
            # 对齐。见 `panel_letter` 的 `anchor_x` / README §10.20。
            _cell = (_lf, _w / 11.5)
            if _kind == "E":
                # rev.r：三层景观（π+D 双轴 + 逐基因 Miss/Syn + 基因条带在**最下面**）
                panel_E_land(fig, _g[0, 0], gcva_dir, letter=_L, log=log, disp_w=_w,
                             pack_cell=_cell)
            elif _kind == "L":
                # rev.o：全基因组突变枢纽棒棒糖图（y = 单位点**去重**携带样本数）
                panel_L(fig, _g[0, 0], gcva_dir, letter=_L, log=log, disp_w=_w,
                        pack_cell=_cell)
            else:
                # 窄栏里长图注会横向溢出、压到左邻面板 → 按列宽折行
                # （8.5 pt 粗体平均字宽 ≈ 0.062 in，留 10 % 余量）
                panel_var(fig, _g[0, 0], _p, _cap, _L, log=log, im=_im,
                          disp_w=_w, title_wrap=max(18, int(_w / 0.068)))
            _x += _w + PACK_GAP
        # 缝隙与 `_gaps` 对齐：最后摆的一行（`_gi == len-1`，即最上面那行）吃 `_gaps[0]`
        _bot += _row_h + _gaps[len(extra_rows) - 1 - _gi]

    _parts = ["detection and genome recovery", "tissue tropism"]
    if args.ef != "none":
        _parts.append("population variation")
    if extra_rows:
        _parts.append("recombination and mutational landscape")
    _parts.append("sequence identity")
    title = ("The known RNA virome of Lycium barbarum: "
             + ", ".join(_parts[:-1]) + " and " + _parts[-1])
    if args.ef == "deg":
        title += " and host response"
    fig.text(0.5, 1.0 - 0.18 / fig_h, title,
             ha="center", va="bottom", fontsize=10, fontweight="bold", color="#333333")

    png = os.path.join(out_dir, f"{args.stem}.png")
    pdf = os.path.join(out_dir, f"{args.stem}.pdf")
    for p in (png, pdf):
        if os.path.islink(p):
            os.unlink(p)
    # bbox_inches="tight": 让画布自适应，避免右/下溢出的轴标签被裁
    # （手工调 right/bottom 边距估不准 —— 本项目 plot_fig1.py 也是这么做的）
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    from PIL import Image as _Img
    with _Img.open(png) as _im:
        log(f"[size] {_im.size[0]} x {_im.size[1]} px @300dpi")
    log(f"[ok ] {png}")
    log(f"[ok ] {pdf}")

    # ---- 可选：把宿主响应（火山图 + GO/KEGG）单独成图 --------------------------
    # 主图的 E/F 让给变异分析后，宿主侧内容就从主图里消失了。这里用**同一批 panel_E/F**
    # 另存一张，保证解耦不丢内容；两行给的高度比主图宽裕，因为不再跟别的面板抢版面。
    if args.host_figure:
        d_v, d_g = 3.30, 3.60                     # 火山图 / 富集四点图
        h2 = 0.30
        H = (d_v + d_g) * (1 + h2 / 2.0)
        fh = 0.55 + H + 0.35
        f2 = plt.figure(figsize=(11.5, fh))
        g2 = f2.add_gridspec(2, 1, left=0.095, right=0.945,
                             top=1.0 - 0.55 / fh, bottom=0.35 / fh,
                             height_ratios=[d_v, d_g], hspace=h2)
        panel_E(f2, g2[0, 0], expr_dir, letter="A", log=log)
        panel_F(f2, g2[1, 0], expr_dir, letter="B", log=log)
        f2.text(0.5, 1.0 - 0.16 / fh,
                "Host transcriptome response to virus infection in Lycium barbarum: "
                "differentially expressed genes and functional enrichment",
                ha="center", va="bottom", fontsize=10, fontweight="bold", color="#333333")
        png2 = os.path.join(out_dir, f"{args.stem}_HostResponse.png")
        pdf2 = os.path.join(out_dir, f"{args.stem}_HostResponse.pdf")
        for p in (png2, pdf2):
            if os.path.islink(p):
                os.unlink(p)
        f2.savefig(png2, dpi=300, bbox_inches="tight", pad_inches=0.03)
        f2.savefig(pdf2, bbox_inches="tight", pad_inches=0.03)
        plt.close(f2)
        with _Img.open(png2) as _im:
            log(f"[size] host figure {_im.size[0]} x {_im.size[1]} px @300dpi")
        log(f"[ok ] {png2}")
        log(f"[ok ] {pdf2}")

    # ---- 可选：GCVA 变异分析**另成一张组合图**（补充材料候选） -------------------
    # 为什么需要它：04 出的成品图**只能通栏**（§七 的两条硬约束：源画布 10 in，
    # 半栏缩到 ≈3.5 pt 全废）。而一张 A4 正文区约 9.5 in 高，除以通栏高 5.5 in
    # → 正文**最多塞得下 1.6 张**变异图。于是「有意义的」变很多，能进正文的极少：
    # 本图把第 2 梯队（AFS / iSNV 密度 / 演化负荷 / 突变谱）收进一张独立图，
    # 放补充材料，不跟正文抢版面。
    if args.var_figure:
        keys = [k.strip() for k in args.var_panels.split(",") if k.strip()]
        pngs, caps, imgs = [], [], []
        for k in keys:
            p, c = resolve_ef("V", k, post_dir, log)
            pngs.append(p)
            caps.append(c)
            imgs.append(_crop_white(p))
        EF_W = (0.945 - 0.095) * 11.5
        d_all = [EF_W * im.height / im.width for im in imgs]
        mean_d = sum(d_all) / len(d_all)
        h3 = min(0.46, 0.95 / mean_d)
        H = sum(d_all) * (1 + h3 * (len(d_all) - 1) / len(d_all))
        fh = 0.55 + H + 0.30
        log(f"[layout] var-figure panels={keys}")
        log(f"[layout] var-figure each_h={[round(d, 2) for d in d_all]} in  "
            f"hspace={h3:.3f}  H={H:.2f} in  fig_h={fh:.2f} in")
        f3 = plt.figure(figsize=(11.5, fh))
        g3 = f3.add_gridspec(len(d_all), 1, left=0.095, right=0.945,
                             top=1.0 - 0.55 / fh, bottom=0.30 / fh,
                             height_ratios=d_all, hspace=h3)
        for i, (im, cap, p) in enumerate(zip(imgs, caps, pngs)):
            panel_var(f3, g3[i, 0], p, cap, "ABCDEFGH"[i], log=log, im=im)
        f3.text(0.5, 1.0 - 0.16 / fh,
                "Within-host variation of Cytorhabdovirus sp. 'lycii' (GCVA): "
                "selection, mutational spectrum and allele-frequency structure",
                ha="center", va="bottom", fontsize=10, fontweight="bold", color="#333333")
        png3 = os.path.join(out_dir, f"{args.stem}_Variation.png")
        pdf3 = os.path.join(out_dir, f"{args.stem}_Variation.pdf")
        for p in (png3, pdf3):
            if os.path.islink(p):
                os.unlink(p)
        f3.savefig(png3, dpi=300, bbox_inches="tight", pad_inches=0.03)
        f3.savefig(pdf3, bbox_inches="tight", pad_inches=0.03)
        plt.close(f3)
        with _Img.open(png3) as _im:
            log(f"[size] variation figure {_im.size[0]} x {_im.size[1]} px @300dpi")
        log(f"[ok ] {png3}")
        log(f"[ok ] {pdf3}")

    # ---- 可选：**自绘**的 9 面板遗传变异分析组合图（正文主图候选） -----------------
    # 为什么不复用 `--var-figure`：那张是嵌 04 成品图，外观与主图 A–D/G 不统一
    # （默认字体/白底灰网格/位图 vs 只留左/下 spine/无网格/Okabe-Ito），
    # 且 04 那批源画布 14–16 in，缩到通栏只剩 4.6–5.2 pt，印刷不可读 —— 实测见 README。
    # 数据全在，所以自绘；版式与字号依据见 README_SCI_composite.md 的 _VariationRich 一节。
    if args.var_rich:
        dirs = ([d for d in args.var_data_dir.split(os.pathsep) if d]
                if args.var_data_dir else [])
        dirs += [
            os.path.join(out_dir, "_diag", "data"),
            os.path.join(out_dir, "_data"),
            os.path.join(os.path.dirname(out_dir), "_diag", "data"),
            os.path.join(gcva_dir, "vcf_viz"),
            os.path.join(gcva_dir, "snpgenie"),
            os.path.join(gcva_dir, "snpeff_macro"),
        ]
        dirs = [d for d in dirs if os.path.isdir(d)]
        if not dirs:
            raise SystemExit("[err] --var-rich 一个候选数据目录都不存在，"
                             "用 --var-data-dir 显式指定")
        log(f"[var-rich] search dirs ({len(dirs)}):")
        for d in dirs:
            log(f"            {d}")
        build_var_rich(dirs, out_dir, args.stem, log=log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
