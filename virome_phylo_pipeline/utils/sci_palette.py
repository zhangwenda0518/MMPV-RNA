#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sci_palette.py — 全项目统一 SCI 出版级配色与 matplotlib 样式。

口径（项目约定，2026-10-12 固化）
--------------------------------
分类色   OKABE_ITO     色盲安全 8+1 色，**顺序固定**（见下），不得随意重排。
树/分支  OKABE_BRANCH  排除最浅的 YELLOW / SKY，保证 0.5 pt 细线在纸上可辨。
连续色   viridis       （需要反向用 viridis_r）
发散色   RdBu_r

硬规则
------
1. **未显式指定颜色的图必须接管 `axes.prop_cycle`**（`apply_sci_style()` 会做），
   否则 matplotlib 默认 tab10 会混进来 —— 这是 2026-10-12 审计发现的配色失控根因。
2. seaborn 的 `violinplot/boxplot` 必须传 `saturation=1`，否则颜色被去饱和。
3. 颜色检查脚本必须先 `tokenize` 剔除注释，否则注释里的色值会误报。
4. **画布尺寸 = 最终排版宽度**：2 栏面板 4.58 in、通栏 9.32 in（Fig 3 的实测栏宽）。
   源画布远大于排版宽度时，字号会被等比压小（源 10 in → 排进 4.58 in，10 pt 只剩 4.6 pt）。

用法
----
    from utils.sci_palette import apply_sci_style, colors, branch_colors, sci_legend
    apply_sci_style(base_size=8)
    ax.plot(x, y, color=colors(1)[0])
    sci_legend(ax, title="Location")
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

__all__ = [
    "OKABE_ITO", "OKABE_BRANCH", "SEQ_CMAP", "DIV_CMAP",
    "C_UP", "C_DOWN", "C_NEUTRAL", "C_REF", "C_OBS", "C_EXP",
    "COL_W_2", "COL_W_1",
    "apply_sci_style", "colors", "branch_colors", "sci_legend", "despine",
    "figsize_for",
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
    因此不调用 `apply_sci_style`"的模块（如 `utils/pipeline_overview_fig.py` 的
    Fig 3）都必须拿**同一个** dict 去 `rcParams.update(...)`。

    为什么必须共享：matplotlib 的 rcParams 是**进程级全局量**。若某模块只钉死
    其中一部分键，它的输出就取决于"同一进程里是否先跑过别的模块" —— 实测
    `pipeline_overview_fig` 单独跑 = 150 dpi / Type3 / 2.37 MB，跟在别的任务后
    跑 = 300 dpi / Type42 / 5.97 MB。**只钉一部分仍会残留差异**：实测把字号/
    字体/dpi/颜色钉死后，`axes.spines.top|right`（面板细边框的显隐）、
    `axes.prop_cycle`、`xtick.major.size`、`lines.linewidth`、`patch.linewidth`、
    `legend.*`、`grid.*` 都还是差异源（残留 mean|Δ| = 0.19、0.46 % 像素差 >8）。
    所以这里把整份 dict 交给调用方，谁都不许自己另写一份。
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

    多图例场景（颜色图例 + 气泡尺寸图例挂在同一个 ax 上）必须把前一条重新挂回
    `artists`，否则**后加的图例会把先加的静默吃掉** —— 2026-10-12 实测：
    面板 A/B 的地点颜色图例就是这样消失的（图上再也查不到"颜色=地点"）。
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
