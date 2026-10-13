#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sci_palette.py — virome_discovery_pipeline 统一 SCI 出版级配色与 matplotlib 样式。

⚠ 与 `virome_phylo_pipeline/utils/sci_palette.py` **同源**。
   两处共用同一套 OKABE_ITO 顺序、同一组语义色、同一份 `sci_rcparams()`。
   任何一侧改了调色板/字号/dpi，另一侧必须同步（否则主图与补充图会串色）。
   基线（2026-10-12 固化）：phylo 侧 `sci_palette.py` md5 `77420a00d7f76fb59c3f3572da5340af`。

口径（项目约定，2026-10-12 固化）
--------------------------------
分类色   OKABE_ITO     色盲安全 8+1 色，**顺序固定**（见下），不得随意重排。
树/分支  OKABE_BRANCH  排除最浅的 YELLOW / SKY，保证 0.5 pt 细线在纸上可辨。
连续色   viridis       （需要反向用 viridis_r）
发散色   RdBu_r

硬规则
------
1. **未显式指定颜色的图必须接管 `axes.prop_cycle`**（`apply_sci_style()` 会做），
   否则 matplotlib 默认 tab10 会混进来 —— 这是配色失控的根因。
2. seaborn 的 `violinplot/boxplot` 必须传 `saturation=1`，否则颜色被去饱和。
3. 颜色检查脚本必须先 `tokenize` 剔除注释，否则注释里的色值会误报。
4. **画布尺寸 = 最终排版宽度**：2 栏面板 4.58 in、通栏 9.32 in。
   源画布远大于排版宽度时字号会被等比压小（源 10 in → 排进 4.58 in，10 pt 只剩 4.6 pt）。
5. **不要用裸 hex 写"业务色"**：verdict / genome type / host / quality 一律走本模块的
   `SEMANTIC`（见下）。discovery 侧报告 HTML（`report_pipeline.py`）与静态 SCI 图
   必须取同一组色值。

用法
----
    from utils.sci_palette import apply_sci_style, colors, SEMANTIC, sci_legend
    apply_sci_style(base_size=8)
    ax.bar(x, y, color=[SEMANTIC["verdict"][v] for v in verdicts])
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── 调色板 ─────────────────────────────────────────────────────────
# 项目统一顺序（与 MEMORY.md 一致，勿改）
OKABE_ITO = [
    "#E69F00",   # orange
    "#56B4E9",   # sky blue
    "#009E73",   # bluish green
    "#F0E442",   # yellow
    "#0072B2",   # blue
    "#D55E00",   # vermillion
    "#CC79A7",   # reddish purple
    "#000000",   # black
    "#999999",   # grey
]

# 树/分支图专用：去掉最浅的 yellow / sky，保证细线可辨
OKABE_BRANCH = [
    "#0072B2",   # blue
    "#D55E00",   # vermillion
    "#009E73",   # green
    "#E69F00",   # orange
    "#CC79A7",   # purple
]

# 语义色（固定含义，跨图一致）
C_UP = "#D55E00"        # 上升 / 显著 / 观测值
C_DOWN = "#0072B2"      # 下降 / 不显著 / 期望值
C_NEUTRAL = "#999999"   # 参考线 / 灰底
C_REF = "#444444"       # 阈值线 / 参考线
C_OBS = "#D55E00"       # 观测统计量
C_EXP = "#0072B2"       # 随机化/期望分布

SEQ_CMAP = "viridis"
DIV_CMAP = "RdBu_r"

# Fig 3 实测栏宽（in）：2 栏面板 / 通栏面板
COL_W_2 = 4.58
COL_W_1 = 9.32

# ── discovery 侧业务语义色 ────────────────────────────────────────
# ⚠ 一律从 OKABE_ITO 取，保证与主图同色系；不要新增游离 hex。
SEMANTIC = {
    # 终判（final_judgement_table.category）
    "verdict": {
        "KEEP": "#009E73",          # bluish green —— 保留
        "REVIEW": "#E69F00",        # orange —— 待核
        "DROP": "#D55E00",          # vermillion —— 丢弃
        "已知植物病毒": "#0072B2",
        "新种候选": "#009E73",
        "远缘病毒候选": "#E69F00",
        "远缘(未定)": "#F0E442",
        "未定": "#999999",
        "宿主基因污染": "#CC79A7",
        "昆虫病毒/污染": "#D55E00",
        "真菌病毒/污染": "#56B4E9",
    },
    # 基因组型（hq_genome_type_summary / all_plant_viruses_genus_summary）
    "genome": {
        "DNA": "#0072B2",
        "RNA": "#D55E00",
        "Viroid": "#009E73",
        "Unassigned": "#999999",
        "none": "#999999",
    },
    # 新病毒等级（05_Taxonomy Novel Rank）
    "novelty": {
        "Known": "#999999",
        "NewSp": "#E69F00",
        "NewGe": "#D55E00",
        "NewFa": "#CC79A7",
    },
    # CheckV 质量分级
    "quality": {
        "Complete": "#0072B2",
        "High-quality": "#009E73",
        "Medium-quality": "#E69F00",
        "Low-quality": "#D55E00",
        "Not-determined": "#CCCCCC",
        # Fig2-E（参考论文 Viruses 2025, 17, 624）保守域完整性三分：
        #   Complete 复用 CheckV Complete 的蓝；Fragment 用 Medium 橙；
        #   absent 用浅灰（"无信息"而非"差"）。
        "Fragment": "#E69F00",
        "No conserved domain": "#CCCCCC",
    },
    # 管线阶段（漏斗 / 流程）
    "stage": {
        "reads": "#999999",
        "assembly": "#0072B2",
        "identification": "#56B4E9",
        # 02b_Filter 高置信过滤（UniProt-strict）。
        # ⚠ 取 **grey `#999999`**（与 `reads` 同）：它是 identification 的**下游
        #   缩减**，不是新口径 —— 用中性灰可让读者一眼看出"这是前一档的过滤产物"，
        #   而不是又引入一个独立语义色。切勿改回 `identification` 的 sky blue
        #   （会与上一档同色而看不出是过滤步骤）。
        "filter": "#999999",
        "cobra": "#009E73",
        "cluster": "#E69F00",
        "taxonomy": "#CC79A7",
        "host": "#D55E00",
        "checkv": "#F0E442",
        "rescue": "#444444",
        "final": "#000000",
    },
}

__all__ = [
    "OKABE_ITO", "OKABE_BRANCH", "SEQ_CMAP", "DIV_CMAP", "SEMANTIC",
    "C_UP", "C_DOWN", "C_NEUTRAL", "C_REF", "C_OBS", "C_EXP",
    "COL_W_2", "COL_W_1",
    "apply_sci_style", "colors", "branch_colors", "sci_legend", "despine",
    "figsize_for", "sci_rcparams",
]


def colors(n: int, start: int = 0):
    """取 n 个分类色（OKABE_ITO 循环，顺序固定）。"""
    n = max(int(n), 0)
    return [OKABE_ITO[(start + i) % len(OKABE_ITO)] for i in range(n)]


def branch_colors(n: int, start: int = 0):
    """取 n 个树/分支色（OKABE_BRANCH 循环）。"""
    n = max(int(n), 0)
    return [OKABE_BRANCH[(start + i) % len(OKABE_BRANCH)] for i in range(n)]


def figsize_for(ncol: int = 1, ratio: float = 0.62, width: float = None):
    """按排版栏宽给画布尺寸。

    ncol=1 → 通栏 9.32 in；ncol=2 → 单栏 4.58 in。
    ratio  = 高/宽（默认 0.62）。
    """
    w = width if width else (COL_W_1 if ncol >= 2 else COL_W_2)
    return (w, w * ratio)


def sci_rcparams(base_size: float = 8.0, *,
                 grid: bool = True,
                 despine_top_right: bool = True,
                 font_family: str = "DejaVu Sans",
                 title_size: float = None,
                 label_size: float = None,
                 tick_size: float = None,
                 legend_size: float = None) -> dict:
    """返回 SCI 样式对应的**完整** rcParams 字典（不写入全局）。

    ⚠ 这是全项目样式的**唯一真源**。`apply_sci_style()` 与任何"自带一套字号体系、
    因此不调用 `apply_sci_style`"的模块都必须拿**同一个** dict 去
    `rcParams.update(...)`。

    为什么必须共享：matplotlib 的 rcParams 是**进程级全局量**。若某模块只钉死
    其中一部分键，它的输出就取决于"同一进程里是否先跑过别的模块" —— 实测同一
    脚本单独跑与跟在别的任务后跑会得到不同 dpi / 字体类型 / 文件大小。
    **只钉一部分仍会残留差异**：`axes.spines.top|right`、`axes.prop_cycle`、
    `xtick.major.size`、`lines.linewidth`、`patch.linewidth`、`legend.*`、
    `grid.*` 都是差异源。所以这里把整份 dict 交给调用方，谁都不许自己另写一份。
    """
    title_size = title_size if title_size is not None else base_size + 1.5
    label_size = label_size if label_size is not None else base_size
    tick_size = tick_size if tick_size is not None else base_size - 0.5
    legend_size = legend_size if legend_size is not None else base_size - 0.5

    rcp = {
        # 字体：英文出版用 DejaVu Sans（图内不出现中文）
        "font.family": "sans-serif",
        "font.sans-serif": [font_family, "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base_size,
        "axes.titlesize": title_size,
        "axes.labelsize": label_size,
        "xtick.labelsize": tick_size,
        "ytick.labelsize": tick_size,
        "legend.fontsize": legend_size,
        "figure.titlesize": title_size + 1,
        "mathtext.fontset": "dejavusans",
        # 颜色：接管 prop_cycle —— 硬规则 1
        "axes.prop_cycle": plt.cycler(color=OKABE_ITO),
        # 版面
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#333333",
        "axes.labelcolor": "#222222",
        "text.color": "#222222",
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "#BBBBBB",
        "legend.facecolor": "white",
        "legend.borderpad": 0.35,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.3,
        "lines.linewidth": 1.2,
        "lines.markersize": 4.0,
        "patch.linewidth": 0.6,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "pdf.fonttype": 42,      # TrueType 嵌入 —— 出版要求可编辑文本
        "ps.fonttype": 42,
        # 网格：**无论开关都写全**，否则 grid.* 会残留上一个模块的值
        "axes.grid": bool(grid),
        "grid.color": "#DDDDDD",
        "grid.linewidth": 0.5,
        "grid.alpha": 0.7,
        "axes.axisbelow": True,
        # 上/右边框：**无论开关都写全**（默认 matplotlib 是 True，
        # apply_sci_style 是 False —— 不写全就是顺序依赖）
        "axes.spines.top": not despine_top_right,
        "axes.spines.right": not despine_top_right,
    }
    return rcp


def apply_sci_style(base_size: float = 8.0, *,
                    grid: bool = True,
                    despine_top_right: bool = True,
                    font_family: str = "DejaVu Sans",
                    title_size: float = None,
                    label_size: float = None,
                    tick_size: float = None,
                    legend_size: float = None):
    """统一 SCI 样式。

    **必须**在绘图前调用：它接管 `axes.prop_cycle`，避免 matplotlib 默认 tab10 混入。
    """
    plt.rcParams.update(sci_rcparams(
        base_size, grid=grid, despine_top_right=despine_top_right,
        font_family=font_family, title_size=title_size,
        label_size=label_size, tick_size=tick_size, legend_size=legend_size))


def despine(ax, top: bool = True, right: bool = True):
    """去掉指定轴的边框（默认去上/右）。"""
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)
    return ax


def _legend_keep_previous(ax, **kw):
    """`ax.legend()` 会**删掉上一条 legend**（`ax.legend_` 被覆盖后旧对象被 remove）。

    多图例场景（颜色图例 + 尺寸图例挂在同一个 ax 上）必须把前一条重新挂回
    `artists`，否则**后加的图例会把先加的静默吃掉**。
    """
    prev = getattr(ax, "legend_", None)
    leg = ax.legend(**kw)
    if prev is not None and prev is not leg:
        try:
            ax.add_artist(prev)
        except Exception:
            pass
    return leg


def sci_legend(ax, *, title: str = None, loc: str = "best",
               ncol: int = 1, frameon: bool = True, fontsize: float = None,
               outside: bool = False, **kw):
    """统一图例样式（白底细框、紧凑、可读）。"""
    if outside:
        kw.setdefault("bbox_to_anchor", (1.02, 1.0))
        loc = "upper left"
    leg = _legend_keep_previous(ax, loc=loc, ncol=ncol, frameon=frameon,
                                fontsize=fontsize, title=title, **kw)
    if title:
        leg.get_title().set_fontsize(
            (fontsize if fontsize is not None else plt.rcParams["legend.fontsize"]))
        leg.get_title().set_fontweight("bold")
    return leg
