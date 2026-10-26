#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_goji_virome_phylo_figs.py — 枸杞 (Lycium barbarum) 科级 ML 系统发育图（独立产出）

背景 / 目的
-----------
复现参考论文 **Viruses 2025, 17, 624**（*Theobroma cacao* 病毒组，PMC12115555）
Fig3–Fig5 的版式：

    Fig3  Kitrinoviricota 的 RdRp 最大似然树（IQ-TREE, 1000 bootstrap），
          **新病毒序列红色高亮**，外群定根；
    Fig4  Pisuviricota 的同款 RdRp ML 树；
    Fig5  Caulimoviridae 的 RT（反转录酶）ML 树。

我们的对应数据 = `09b_Analysis_Verify/acvirus_trees/<Family>/phylogeny.treefile`
（21 个科级 IQ-TREE ML 树；树叶里 `NCBI accession` = 参考序列、
`CRR*/SRR*_clean_NODE_*` = **我们自己的 contig** → 判别「新病毒」的依据）。

⚠ **数值不与论文相同、也不应相同**：论文 = 109 可可文库 / 34 条手策序列；
我们 = 445 枸杞文库 / 20,889 vOTU。本次复现的是 **版式与面板语义**。

⚠ **本脚本只出"门级合并大图"**（把同门各科的树竖排拼成一张），不重跑 IQ-TREE。

产出
----
    <outdir>/Fig_GojiVirome_Phylo.{png,pdf}        三栏（Kitrinoviricota / Pisuviricota / Caulimoviridae）
    <outdir>/panels/phylo_<Phylum>.{png,pdf}       单门面板
    <outdir>/_build/phylo_manifest.json            每个门用了哪些科、叶数、自有 contig 数

用法
----
    python reports/virome_analysis_pipeline/make_goji_virome_phylo_figs.py \
      --trees <dir containing <Family>/phylogeny.treefile> \
      --outdir <...>/S10_SCI_Discovery \
      [--virus-name "Lycium barbarum"] [--dpi 300]
"""
from __future__ import annotations

# ── 图件四格式导出 (SVG / PNG / PDF / PNG@1200dpi) — 默认开启 ──────────
# 2026-10-20: 所有出图脚本统一口径; 用 MMPV_FIG_EXPORT=0 可整体关闭,
# MMPV_FIG_FORMATS=svg,png,pdf,hires 可裁剪, MMPV_FIG_HIRES_DPI 可改高清 dpi。
try:
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
except ImportError:  # 仓库根不在 sys.path (脚本被单独执行) → 向上找 mmpv_common/
    import os as _os
    import sys as _sys
    _d = _os.path.dirname(_os.path.abspath(__file__))
    while _d and not _os.path.isdir(_os.path.join(_d, "mmpv_common")):
        _parent = _os.path.dirname(_d)
        if _parent == _d:
            break
        _d = _parent
    if _d and _d not in _sys.path:
        _sys.path.insert(0, _d)
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
_install_fig_hook()
# ──────────────────────────────────────────────────────────────────────

import argparse
import csv
import json
import re
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np

from Bio import Phylo

# ── 统一配色（显式文件路径加载，理由同 make_discovery_sci_figs.py）──
import importlib.util as _ilu

_HERE = Path(__file__).resolve()
_REPO = _HERE.parents[2]

_PAL = _REPO / "virome_discovery_pipeline" / "utils" / "sci_palette.py"
_spec = _ilu.spec_from_file_location("sci_palette", _PAL)
_sp = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_sp)

OKABE_ITO = _sp.OKABE_ITO
apply_sci_style = _sp.apply_sci_style
despine = _sp.despine

DEFAULT_DPI = 300

# ── 门 → 该门下的科（来自 final_integrated_classification.tsv 的实测映射）──
# ⚠ 只放**我们实际建了树**的科（`acvirus_trees/` 里有 `phylogeny.treefile` 的）。
#   映射依据 = 10_Reports/final_integrated_classification.tsv 的 Phylum→Family 实测，
#   不是凭记忆写的。`Draupnirviridae` 在本数据集未落到任何门（Unclassified）→ 单列。
PHYLUM_FAMILIES = {
    "Kitrinoviricota": [
        "Alphaflexiviridae", "Betaflexiviridae", "Bromoviridae",
        "Closteroviridae", "Kitaviridae", "Mayoviridae", "Tymoviridae",
        "Tombusviridae", "Benyviridae", "Virgaviridae", "Botourmiaviridae",
    ],
    "Pisuviricota": [
        "Potyviridae", "Secoviridae", "Solemoviridae", "Iflaviridae",
    ],
    "Artverviricota (RT)": [     # 论文 Fig5 = Caulimoviridae 的 RT 树（属 Artverviricota 旁支）
        "Caulimoviridae", "Metaviridae", "Pseudoviridae",
    ],
    "Cressdnaviricota": [
        "Geminiviridae",
    ],
    "Negarnaviricota": [
        "Rhabdoviridae",
    ],
    "Unclassified": [            # 本数据集未落到任何门
        "Draupnirviridae",
    ],
}

# 面板标题里显示的「门 / 家族块」名与论文对应
PANEL_META = {
    "Kitrinoviricota": ("Kitrinoviricota · RdRp (ML)", "Fig3"),
    "Pisuviricota": ("Pisuviricota · RdRp (ML)", "Fig4"),
    "Artverviricota (RT)": ("Artverviricota (Caulimoviridae/Metaviridae/Pseudoviridae) · RT (ML)", "Fig5"),
    "Cressdnaviricota": ("Cressdnaviricota · Rep/CP (ML)", "—"),
    "Negarnaviricota": ("Negarnaviricota · RdRp (ML)", "—"),
    "Unclassified": ("Unclassified families (ML)", "—"),
}

# 自有 contig 的命名模式。⚠ **两种**（实测 2026-10-12）：
#   ① `CRR*/SRR*_clean_NODE_*`  —— 宏病毒组组装出的 node（多数科）
#   ② `contig_*`                —— 部分科（Benyviridae / Closteroviridae）用这套命名
# 只认 ① 会**静默漏报**：实测 Benyviridae 树里明明有 `contig_4185`，
# 旧正则判为 0 → 该面板被当成"无自有序列"，而管线 HTML 报告里它是有高亮的。
_OWN_RE = re.compile(r"^(?:(?:CRR|SRR)\d+_clean_NODE_\d+|contig_\d+)")


def _is_own(name: str) -> bool:
    """判断叶名是不是「我们自己的 contig」——这是「新病毒红色高亮」的判据。"""
    return bool(_OWN_RE.match(name or ""))


def _short(name: str, maxlen: int = 30) -> str:
    """叶标签缩写。

    ⚠ 自有 contig 必须带**样本前缀**：不同样本各有 `NODE_1`、`NODE_5`… → 只写
      `own NODE_1` 会撞名（实测 Rhabdoviridae 里两处 `own NODE_1` 重叠）。
      参考 accession 原样（已唯一）。⚠ `contig_N` 命名**本身已唯一**（无样本概念）
      → 原样返回，不要硬套 `NODE_` 前缀。
    """
    if _is_own(name):
        m = re.match(r"^((?:CRR|SRR)\d+)_clean_NODE_(\d+)", name)
        if m:
            return f"{m.group(1)[-6:]}:NODE_{m.group(2)}"
        m2 = re.search(r"_clean_NODE_(\d+)", name)
        if m2:
            return f"own NODE_{m2.group(1)}"
        return name if len(name) <= maxlen else name[:maxlen - 1] + "…"
    return name if len(name) <= maxlen else name[:maxlen - 1] + "…"


def read_tree(path: Path):
    try:
        return Phylo.read(str(path), "newick")
    except Exception:
        try:
            return Phylo.read(str(path), "nexus")
        except Exception:
            return None


def ladderize(tree):
    try:
        tree.ladderize()
    except Exception:
        pass
    return tree


def draw_tree(ax, tree, *, own_color="#D55E00", ref_color="#333333",
              max_ref_labels=24):
    """在 ax 上画一棵 ladderized 树。

    标签策略（关键）：**只给自有 contig 打标签**，参考 accession 只在叶数少时显示
    若干。理由：Rhabdoviridae 有 173 叶、Caulimoviridae 有 94 叶 —— 全打标签会糊成
    一片黑（实测 100 % 重叠）。论文里也是「树形承载结构、红字只标新种」。
    """
    tree = ladderize(tree)
    terminals = tree.get_terminals()
    n = len(terminals)
    own_labels = {_short(t.name) for t in terminals if _is_own(t.name)}
    # 参考叶：叶数少（≤ max_ref_labels）时显示，否则只显示自有叶。
    show_ref = n <= max_ref_labels

    def _lab(c):
        if not c.is_terminal():
            return ""
        s = _short(c.name)
        if _is_own(c.name):
            return s
        return s if show_ref else ""

    _draw_kw = dict(axes=ax, do_show=False, label_func=_lab,
                    label_colors=(lambda lab: own_color
                                  if lab.strip() in own_labels else ref_color),
                    branch_labels=(lambda c: "" if c.is_terminal() else ""),
                    show_confidence=False)
    try:
        Phylo.draw(tree, **_draw_kw)
    except TypeError:
        # 老版本 `Phylo.draw` 不认 `label_colors` / `show_confidence` → 逐个剔除重试
        for k in ("label_colors", "show_confidence", "branch_labels"):
            _draw_kw.pop(k, None)
        Phylo.draw(tree, **_draw_kw)
    # Bio.Phylo 的 `label_fs` 关键字在部分版本不存在 → 统一在画完后按 Text 对象改字号
    # ⚠ `t.get_text()` 带 Bio.Phylo 加的前导空格 → 必须 strip 后再与 own_labels 比。
    for t in ax.texts:
        _own = t.get_text().strip() in own_labels
        t.set_fontsize(5.0 if _own else 3.4)
        if _own:
            t.set_fontweight("bold")
    ax.set_axis_off()
    ax.set_title("")
    # ⚠ Bio.Phylo 会自加 `branch length` (xlabel) 与 `taxa` (ylabel)；`set_axis_off()`
    #   只隐藏轴线/刻度/spine，**不隐藏 axis label** → 相邻子图的 xlabel 会互相
    #   "重叠"（实测组合图 226 对，全是这两类标签）。必须显式清空。
    ax.set_xlabel("")
    ax.set_ylabel("")
    _declash_own_labels(ax, own_labels)
    return n, len(own_labels)


def _declash_own_labels(ax, own_labels, min_gap_px=6.0):
    """只对**自有 contig 标签**做竖向错行（参考 accession 不改）。

    ⚠ 为什么只动自有标签：参考叶有几十上百个，全动会变"竖直引线贯穿全图"
      （splitstree 那次四次修法全被否的教训 —— 视觉最差）。自有 contig 只有 1–8 个，
      错行范围小、且它们才是读者要看的对象。

    做法 = 取所有自有标签，按 y 升序；自下而上依次保证与前一个的**竖向净空 ≥
    `min_gap_px`**，不够就把当前这个往下挪。整条链一次扫完（不是两两独立）。
    """
    try:
        ax.figure.canvas.draw()
        rend = ax.figure.canvas.get_renderer()
    except Exception:
        return
    tags = []
    for t in ax.texts:
        # ⚠ Bio.Phylo 会在标签前**加一个空格**（`" " + label`）→ `t.get_text()` 带前导
        #   空格，而 `own_labels` 是不带空格的 `_short()` 值 → 直接 `in` 判定恒为 False
        #   （这就是错行"看起来没生效"的原因）。必须 strip 后再比。
        if t.get_text().strip() in own_labels:
            b = t.get_window_extent(renderer=rend)
            tags.append([t, b.y0, b.y1, b])
    if len(tags) < 2:
        return
    tags.sort(key=lambda z: z[1])                     # 按 y0 升序
    inv = ax.transData.inverted()
    prev_top = None
    for tag in tags:
        t, y0, y1, b = tag
        if prev_top is not None and y0 < prev_top + min_gap_px:
            dy_px = (prev_top + min_gap_px) - y0
            xd, yd = t.get_position()
            # px 位移 → 数据坐标位移（y 方向）
            _, y_a = inv.transform((b.x0, b.y0))
            _, y_b = inv.transform((b.x0, b.y0 + dy_px))
            t.set_position((xd, yd + (y_b - y_a)))
            y0, y1 = y0 + dy_px, y1 + dy_px
        prev_top = y1


def _tree_stats(tree):
    terms = tree.get_terminals()
    own = [t.name for t in terms if _is_own(t.name)]
    return len(terms), len(own)


def panel_for_phylum(ax, trees_dir: Path, phylum: str, families: list):
    """把一个门下的各科树竖排画在同一个 ax 里。

    返回 (used_families, total_leaves, total_own)。
    """
    # 收集可用的科树
    avail = []
    for fam in families:
        p = trees_dir / fam / "phylogeny.treefile"
        if p.is_file():
            t = read_tree(p)
            if t is not None and len(t.get_terminals()) >= 4:
                avail.append((fam, t, p))
    if not avail:
        ax.axis("off")
        ax.text(0.5, 0.5, f"(no tree for {phylum})", ha="center", va="center",
                fontsize=7, color="#999999", transform=ax.transAxes)
        return [], 0, 0

    n_fam = len(avail)
    # 每个科的 x 范围独立（树的枝长尺度各异）→ 用「子轴纵排」而非共享 x
    # 这里用 GridSpecFromSubplotSpec 不方便（ax 已给定）→ 手工按比例切 ax 的位置
    pos = ax.get_position()
    fig = ax.figure
    ax.set_axis_off()
    x0, y0, w, h = pos.x0, pos.y0, pos.width, pos.height
    # ⚠ 顶部预留一条「门标题带」：否则第一个科的科名会与门标题相撞（实测）。
    top_band = 0.055 * h
    h_eff = h - top_band
    gap = 0.10 / max(n_fam, 1)
    hh = (h_eff - gap * (n_fam - 1)) / n_fam
    total_leaves = total_own = 0
    used = []
    for i, (fam, t, p) in enumerate(avail):
        sub = fig.add_axes([x0, y0 + (n_fam - 1 - i) * (hh + gap), w, hh])
        n, no = _tree_stats(t)
        total_leaves += n
        total_own += no
        used.append({"family": fam, "leaves": n, "own_contigs": no,
                     "treefile": str(p)})
        draw_tree(sub, t)
        # 科名放**该子轴的右下角内侧**：既不压枝、也不与门标题/上一棵树相撞。
        sub.text(1.0, -0.02,
                 f"{fam}  (n={n}" + (f", own={no}" if no else "") + ")",
                 transform=sub.transAxes, fontsize=5.0, fontweight="bold",
                 va="bottom", ha="right", color="#222222",
                 path_effects=_stroke())
    return used, total_leaves, total_own


def _stroke(lw=2.2):
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=lw, foreground="white")]


# ══════════════════════════════════════════════════════════════════
# 组装
# ══════════════════════════════════════════════════════════════════
def save(fig, outdir: Path, stem: str, dpi: int):
    outdir.mkdir(parents=True, exist_ok=True)
    png, pdf = outdir / f"{stem}.png", outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    fig.savefig(pdf, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    return png, pdf


def _phylum_weight(trees_dir: Path, ph: str) -> float:
    """门内「可用科的叶数总和 / 该门科数」→ 用于按信息量分配画布高度。

    优先按**叶数和中位**估，避免个别超大科（Rhabdoviridae 173 叶）吃掉全部高度。
    """
    counts = []
    for fam in PHYLUM_FAMILIES[ph]:
        p = trees_dir / fam / "phylogeny.treefile"
        if p.is_file():
            t = read_tree(p)
            if t is not None:
                counts.append(len(t.get_terminals()))
    if not counts:
        return 1.0
    counts.sort()
    med = counts[len(counts) // 2]
    # 高度 ∝ 科数 × log(中位叶数+1)，对数压缩避免超级科一家独大
    return max(1.0, len(counts) * np.log10(med + 10))


def _family_row_height(trees_dir: Path, fam: str) -> float:
    """一个科的画布高度权重 ∝ log10(叶数)。"""
    p = trees_dir / fam / "phylogeny.treefile"
    t = read_tree(p) if p.is_file() else None
    n = len(t.get_terminals()) if t is not None else 10
    return max(0.6, np.log10(n + 8))


def build_composite(trees_dir: Path, outdir: Path, dpi: int, virus_name: str):
    """2 列网格，一科一格；**叶数 ≥ 60 的科跨 2 列独占一行**（叶标签才排得开）。

    ⚠ 版式理由：Kitrinoviricota 有 11 个科、Negarnaviricota 有 1 个 173 叶的大科；
       "一个门挤成一竖条"会把 11 个科的叶标签压到一起（实测不可读）。
       且 173 叶的 Rhabdoviridae 若只占半页宽，自有 contig 标签会与相邻叶相叠
       （实测 3 对 OVERLAP）→ 大科跨全宽。
    """
    # 展开成 (phylum, family) 列表，跳过没有树的科
    rows = []
    for ph, fams in PHYLUM_FAMILIES.items():
        for fam in fams:
            p = trees_dir / fam / "phylogeny.treefile"
            if p.is_file():
                t = read_tree(p)
                if t is not None and len(t.get_terminals()) >= 4:
                    rows.append((ph, fam, t, p))
    if not rows:
        sys.exit("[fatal] no usable trees found")

    WIDE_AT = 60                       # 叶数 ≥ 此值 → 跨 2 列
    # 先按「两列填格 + 大科独占一行」排成 (row, col, colspan) 序列
    slots = []
    ncol = 2
    cur_col = 0
    for ph, fam, t, p in rows:
        n = len(t.get_terminals())
        span = ncol if n >= WIDE_AT else 1
        if span == ncol:
            if cur_col != 0:           # 半格未填满 → 先补空位，避免错行
                slots.append(None)
                cur_col = 0
            slots.append((ph, fam, t, p, ncol))
            cur_col = 0
        else:
            slots.append((ph, fam, t, p, 1))
            cur_col = (cur_col + 1) % ncol
            if cur_col == 0:
                pass
    if cur_col != 0:
        slots.append(None)             # 补齐最后一行的空格
    nrow = (len(slots) + ncol - 1) // ncol

    # 行高权重
    row_w = []
    i = 0
    while i < len(slots):
        s = slots[i]
        if s is None:
            row_w.append(0.8)
            i += 1
            continue
        spans = []
        j = i
        while j < len(slots) and j < i + ncol:
            if slots[j] is not None:
                spans.append(_family_row_height(trees_dir, slots[j][1]))
            j += 1
        row_w.append(max(spans) if spans else 0.8)
        i += ncol

    fig = plt.figure(figsize=(9.2, 1.0 + 1.30 * sum(row_w)))
    gs = GridSpec(nrow, ncol, figure=fig, hspace=0.30, wspace=0.15,
                  left=0.030, right=0.985, top=0.966, bottom=0.028,
                  height_ratios=row_w)
    manifest = {}
    for i, s in enumerate(slots):
        if s is None:
            continue
        ph, fam, t, p, span = s
        r, c = divmod(i, ncol)
        ax = fig.add_subplot(gs[r:r + 1, c:c + span])
        n, no = _tree_stats(t)
        draw_tree(ax, t)
        lab = f"{fam}  (n={n}" + (f", own={no}" if no else "") + ")"
        ax.text(0.0, 1.01, lab, transform=ax.transAxes, fontsize=6.2,
                fontweight="bold", va="bottom", ha="left", color="#111111")
        manifest[fam] = {"phylum": ph, "leaves": n, "own_contigs": no,
                         "colspan": span, "treefile": str(p),
                         "paper_fig": PANEL_META.get(ph, ("", "—"))[1]}
    fig.legend(
        handles=[Line2D([], [], color="#D55E00", lw=1.8,
                        label="own contig (putative novel virus)"),
                 Line2D([], [], color="#333333", lw=1.8,
                        label="reference sequence (NCBI)")],
        loc="lower center", ncol=2, fontsize=6.4, frameon=False,
        bbox_to_anchor=(0.5, 0.002))
    fig.suptitle(
        f"Maximum-likelihood phylogenies of plant-associated virus families in the "
        f"{virus_name} RNA virome\n"
        f"per-family IQ-TREE trees; own contigs highlighted in orange",
        fontsize=8.6, fontweight="bold", y=0.989)
    png, pdf = save(fig, outdir, "Fig_GojiVirome_Phylo", dpi)
    plt.close(fig)
    return png, pdf, manifest


def build_panels(trees_dir: Path, outdir: Path, dpi: int):
    """单科面板：每科一张（宽画布），供排版单独取用。"""
    pdir = outdir / "panels"
    pdir.mkdir(parents=True, exist_ok=True)
    made = []
    for ph, fams in PHYLUM_FAMILIES.items():
        for fam in fams:
            p = trees_dir / fam / "phylogeny.treefile"
            if not p.is_file():
                continue
            t = read_tree(p)
            if t is None or len(t.get_terminals()) < 4:
                continue
            n, no = _tree_stats(t)
            h = max(2.6, 0.062 * n)
            fig = plt.figure(figsize=(7.0, h))
            ax = fig.add_axes([0.03, 0.02, 0.94, 0.92])
            draw_tree(ax, t)
            fig.suptitle(f"{fam}  (n={n}" + (f", own={no}" if no else "") + ")",
                         fontsize=8.0, fontweight="bold", y=0.988)
            save(fig, pdir, f"phylo_{fam}", dpi)
            plt.close(fig)
            made.append(fam)
    return made


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True,
                    help="含 <Family>/phylogeny.treefile 的目录")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--virus-name", default="Lycium barbarum")
    ap.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    ap.add_argument("--panels-only", action="store_true")
    ap.add_argument("--composite-only", action="store_true")
    a = ap.parse_args()

    trees_dir = Path(a.trees).expanduser().resolve()
    if not trees_dir.is_dir():
        sys.exit(f"[fatal] --trees not a directory: {trees_dir}")
    outdir = Path(a.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    apply_sci_style(base_size=7.5, grid=False)

    print(f"[in]  trees dir = {trees_dir}")
    man = {"trees_dir": str(trees_dir), "virus_name": a.virus_name, "panels": {}}
    if not a.panels_only:
        png, pdf, m = build_composite(trees_dir, outdir, a.dpi, a.virus_name)
        man["panels"] = m
        print(f"[out] composite = {png.name} / {pdf.name}")
    if not a.composite_only:
        made = build_panels(trees_dir, outdir, a.dpi)
        print(f"[out] panels    = {len(made)} files in {outdir / 'panels'}")

    mpath = outdir / "_build" / "phylo_manifest.json"
    mpath.parent.mkdir(parents=True, exist_ok=True)
    mpath.write_text(json.dumps(man, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[out] manifest  = {mpath}")


if __name__ == "__main__":
    main()
