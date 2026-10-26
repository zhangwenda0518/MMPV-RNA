#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_goji_virome_phylosynteny_fig.py — 枸杞病毒「ML 树 + Synteny 两联」独立图

背景 / 目的
-----------
用户的明确指令：给系统发育产物 **补 Synteny 右联**。

⚠ **口径（用户已定，2026-10-12）**：
  · **嵌入管线原图，不重绘** —— 管线 `utils/acvirus_tree_pro.py`
    的 `draw_tree_synteny_figure()` 已为**每个科**产出
    `Tree_Synteny_Composite.{png,pdf}`（版式 = 左 `Maximum-Likelihood Phylogeny`
    + 右 `Comparative Synteny & Homology`，含 Genome Coordinates (bp) 轴、
    Protein Identity (%) 色标 30–100、Sample_Contig 高亮、属级背景色块）。
    这**正是**参考论文 / 管线报告 HTML 里的版式。
    项目硬规则：「整页独立图只能自绘」**已被否决** → **原图直接嵌入，不许重绘/降采样**。
  · **保留纯树版**：现有 `Fig_GojiVirome_Phylo`（`make_goji_virome_phylo_figs.py`
    画的 21 科纯树）**不动**；本脚本**另出** `Fig_GojiVirome_PhyloSynteny`。

产出
----
    <outdir>/Fig_GojiVirome_PhyloSynteny.{png,pdf}   21 行全宽，一科一行（左树+右 Synteny 原图）
    <outdir>/panels_phylosynteny/<Family>.{png,pdf}  单科面板（原图整幅嵌入，供单独取用）
    <outdir>/_build/phylosynteny_manifest.json       每科的像素尺寸、缩放、来源路径

版式依据（项目硬规则）
----------------------
  ⭐ **面板宽必须匹配源画布宽**：源图一律 **4477 px 宽**（Potyviridae 4513）→
     属**通栏**（≥ 8.8 in 级）→ 本图按通栏 9.32 in 排气。
  ⭐ **不许拉伸**：`imshow(..., aspect="equal")` 保持原图长宽比，行高按
     `h_src / w_src × w_panel` 算，**不做任何非等比缩放**。
  ⚠ 源图高度差异极大（Mayoviridae 1715 px → Rhabdoviridae 19142 px，11×）→
     若按原比例直排，单张 Rhabdoviridae 会占据整页 → 对**超高图**设上限
     `MAX_ROW_IN`，超限则该行给足上限高度、整幅**等比缩小**并记 manifest
     （**不裁切、不重排**）。

用法
----
    python reports/virome_analysis_pipeline/make_goji_virome_phylosynteny_fig.py \
      --trees <dir containing <Family>/Tree_Synteny_Composite.png> \
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
import json
import re
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import importlib.util as _ilu

from Bio import Phylo

_HERE = Path(__file__).resolve()
_REPO = _HERE.parents[2]

# ── 统一配色（显式文件路径加载：服务器 site-packages/utils 是常规包，
#    `from utils.sci_palette import` 必失败）──
_PAL = _REPO / "virome_discovery_pipeline" / "utils" / "sci_palette.py"
_spec = _ilu.spec_from_file_location("sci_palette", _PAL)
_sp = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_sp)

OKABE_ITO = _sp.OKABE_ITO
apply_sci_style = _sp.apply_sci_style
despine = _sp.despine

DEFAULT_DPI = 300

# 通栏宽（in）：与 make_discovery_sci_figs.py / 项目版式规则同源
COL_W_1 = getattr(_sp, "COL_W_1", 9.32)

# ── 门 → 科 映射：与 make_goji_virome_phylo_figs.py **完全一致**（同源实测映射）──
#    本图的行序 = 该映射的展开序，与纯树版对齐，便于两图并排阅读。
PHYLUM_FAMILIES = {
    "Kitrinoviricota": [
        "Alphaflexiviridae", "Betaflexiviridae", "Bromoviridae",
        "Closteroviridae", "Kitaviridae", "Mayoviridae", "Tymoviridae",
        "Tombusviridae", "Benyviridae", "Virgaviridae", "Botourmiaviridae",
    ],
    "Pisuviricota": [
        "Potyviridae", "Secoviridae", "Solemoviridae", "Iflaviridae",
    ],
    "Artverviricota (RT)": [
        "Caulimoviridae", "Metaviridae", "Pseudoviridae",
    ],
    "Cressdnaviricota": [
        "Geminiviridae",
    ],
    "Negarnaviricota": [
        "Rhabdoviridae",
    ],
    "Unclassified": [
        "Draupnirviridae",
    ],
}

# 自有 contig 命名模式：与纯树脚本同源（两种命名，缺一即静默漏报）
_OWN_RE = re.compile(r"^(?:(?:CRR|SRR)\d+_clean_NODE_\d+|contig_\d+)")


def png_size(p: Path):
    """读 PNG 的 IHDR 拿 (w, h)（不依赖 PIL，服务器上也稳）。"""
    with open(p, "rb") as f:
        d = f.read(33)
    if d[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return (int.from_bytes(d[16:20], "big"), int.from_bytes(d[20:24], "big"))


def tree_stats(trees_dir: Path, fam: str):
    """(叶数, 自有 contig 数) —— 从树的 treefile 数；缺 treefile 时返回 (None, None)。"""
    p = trees_dir / fam / "phylogeny.treefile"
    if not p.is_file():
        return None, None
    try:
        t = Phylo.read(str(p), "newick")
    except Exception:
        try:
            t = Phylo.read(str(p), "nexus")
        except Exception:
            return None, None
    terms = t.get_terminals()
    return len(terms), sum(1 for x in terms if _OWN_RE.match(x.name or ""))


def save(fig, outdir: Path, stem: str, dpi: int):
    outdir.mkdir(parents=True, exist_ok=True)
    png, pdf = outdir / f"{stem}.png", outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    fig.savefig(pdf, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    return png, pdf


def _collect(trees_dir: Path, exclude=frozenset()):
    """按 PHYLUM_FAMILIES 的顺序收集 (phylum, family, img, w, h, n, own)。

    ⚠ `exclude`（2026-10-15 v24 新增）：跳过的科名集合。用于出**补充附图** ——
      Alphaflexiviridae 已被提升进主图（`Fig_GojiVirome_Discovery` 面板 R），
      其余 20 科走 `Fig_GojiVirome_PhyloSynteny_Supp`，**同一科不得两处重复交付**。
    """
    rows = []
    for ph, fams in PHYLUM_FAMILIES.items():
        for fam in fams:
            if fam in exclude:
                continue
            img = trees_dir / fam / "Tree_Synteny_Composite.png"
            if not img.is_file():
                continue
            sz = png_size(img)
            if not sz:
                continue
            n, own = tree_stats(trees_dir, fam)
            rows.append((ph, fam, img, sz[0], sz[1], n, own))
    return rows


def build_composite(trees_dir: Path, outdir: Path, dpi: int, virus_name: str,
                    exclude=frozenset(), stem: str = "Fig_GojiVirome_PhyloSynteny",
                    supp_note: str = ""):
    """21 行 × 1 列，一科一行，全幅嵌入管线原图（不重绘、不拉伸）。

    ⚠ v24：`exclude` / `stem` / `supp_note` 三个参数用于出**补充附图** ——
      `exclude={"Alphaflexiviridae"}` → 20 科；`stem` 换成 `..._Supp`；
      `supp_note` 追加到 suptitle 第二行，说明本图是补充图、主图已收录哪一科。
    """
    rows = _collect(trees_dir, exclude)
    if not rows:
        sys.exit("[fatal] no Tree_Synteny_Composite.png found under --trees")

    W_FIG = COL_W_1                       # 通栏宽（in）
    M_L = M_R = 0.10                      # 左右边距（in）→ 面板净宽
    w_panel = W_FIG - M_L - M_R           # ≈ 9.12 in

    # ── 行高：按原图长宽比等比算 → 需的高度 ──
    #   h_in = h_src / w_src * w_panel
    # ⚠ 源图高度差 11×（Mayoviridae 1715 → Rhabdoviridae 19142 px）。
    #   两个约束取折中：
    #     · **下限** = 叶数 × MIN_ROW_IN → 保证叶标签可读（不许把多叶科压扁）；
    #     · **上限** = MAX_ROW_IN → 单行不超过一页，否则 Rhabdoviridae 一张
    #       （原比例 39 in）会独占整张图。
    #   上限触发时**整幅等比缩小**（不裁切/不重排），并在 manifest 标 `height_capped`。
    MIN_ROW_IN = 0.145                    # 单叶最小行距（in）→ 173 叶需 ~25 in
    MAX_ROW_IN = 9.0                      # 单行画布上限（in）
    H_MIN_ROW = 0.9                       # 单行下限
    h_rows = []
    meta_scale = []
    for ph, fam, img, w_src, h_src, n, own in rows:
        h_nat = h_src / float(w_src) * w_panel
        h_floor = MIN_ROW_IN * (n or 10)
        h_use = min(max(h_nat, h_floor, H_MIN_ROW), MAX_ROW_IN)
        h_rows.append(h_use)
        meta_scale.append({"natural_h_in": round(h_nat, 3),
                           "floor_h_in": round(h_floor, 3),
                           "capped": bool(h_nat > MAX_ROW_IN + 1e-9),
                           "expanded": bool(h_nat < h_floor - 1e-9),
                           "used_h_in": round(h_use, 3)})

    CAP_IN = 0.30                         # 每行顶部题注净空（in）
    hs = 0.34                             # 行间距（in）
    TOP_IN = 0.95                         # 图顶 → 首行的净空（in），给两行 suptitle 让位
    H_FIG = sum(h + CAP_IN for h in h_rows) + hs * (len(rows) - 1) + TOP_IN + 0.12
    fig = plt.figure(figsize=(W_FIG, H_FIG))

    # 自顶向下逐行放：用绝对坐标（figure fraction）排，避免 GridSpec 行高换算误差
    y_in = H_FIG - TOP_IN                 # 首行「题注带」的上沿
    x0 = M_L / W_FIG
    wp = w_panel / W_FIG
    manifest = {}
    for k, (ph, fam, img, w_src, h_src, n, own) in enumerate(rows):
        h_use = h_rows[k]
        # ── 题注带（紧贴本行图上方）──
        cap_top = y_in
        y_in -= CAP_IN
        top = y_in / H_FIG
        bottom = (y_in - h_use) / H_FIG
        cap = f"{fam}   (n={n}" + (f", own={own}" if own else "") + f")   ·   {ph}"
        ax_cap = fig.add_axes([x0, top, wp, CAP_IN / H_FIG])
        ax_cap.axis("off")
        ax_cap.text(0.0, 0.0, cap, fontsize=6.8, fontweight="bold",
                    va="bottom", ha="left", color="#111111")
        ax = fig.add_axes([x0, bottom, wp, h_use / H_FIG])
        try:
            im = plt.imread(img)
        except Exception:
            traceback.print_exc()
            ax.axis("off")
            ax.text(0.5, 0.5, f"{fam}: image load failed",
                    ha="center", va="center", color="#999999")
        else:
            ax.imshow(im, aspect="auto")   # aspect="auto" 填满已按原比例算好的框
            ax.set_axis_off()
        # 面板宽/字号换算说明（印刷字号 = 源字号 × 面板宽 / 源画布宽）
        scale = (w_panel * dpi) / float(w_src)
        manifest[fam] = {
            "phylum": ph, "leaves": n, "own_contigs": own,
            "src_png": str(img), "src_px": [w_src, h_src],
            "row_h_in": round(h_use, 3),
            "natural_h_in": meta_scale[k]["natural_h_in"],
            "height_capped": meta_scale[k]["capped"],
            "print_scale": round(scale, 3),
        }
        y_in -= h_use + hs
    _ = cap_top  # 保留变量语义（题注带上沿），便于日后调试

    # ⚠ suptitle 的 y 是 **figure 分数**，而 TOP_IN 是**英寸** → 必须显式换算，
    #   否则 H_FIG 一变（本脚本 H_FIG 随叶数变化，可达数十英寸）两行标题就会压住首行图。
    #   在留白（TOP_IN）的中部落第一行、第二行在其下。
    y_title = 1.0 - (TOP_IN * 0.22) / H_FIG
    fig.suptitle(
        f"Comparative synteny and homology of plant-associated virus families in the "
        f"{virus_name} RNA virome\n"
        f"per-family maximum-likelihood phylogeny (left) with pairwise protein-identity "
        f"co-linearity (right); own contigs highlighted"
        + (f"\n{supp_note}" if supp_note else ""),
        fontsize=8.6, fontweight="bold", y=y_title)
    png, pdf = save(fig, outdir, stem, dpi)
    plt.close(fig)
    return png, pdf, manifest


def build_panels(trees_dir: Path, outdir: Path, dpi: int, exclude=frozenset()):
    """单科面板：原图**整幅嵌入**（不缩放不裁切），供排版单独取用。

    ⚠ 面板画布 = 原图本身的尺寸（等比），不做任何重绘。
    ⚠ `exclude`（2026-10-15 v24）：跳过的科（已进主图者），避免重复交付。
    """
    pdir = outdir / "panels_phylosynteny"
    pdir.mkdir(parents=True, exist_ok=True)
    made = []
    for ph, fams in PHYLUM_FAMILIES.items():
        for fam in fams:
            if fam in exclude:
                continue
            img = trees_dir / fam / "Tree_Synteny_Composite.png"
            if not img.is_file():
                continue
            sz = png_size(img)
            if not sz:
                continue
            w_src, h_src = sz
            w_in = COL_W_1
            h_in = h_src / float(w_src) * w_in
            n, own = tree_stats(trees_dir, fam)
            fig = plt.figure(figsize=(w_in, h_in))
            ax = fig.add_axes([0, 0, 1, 1])
            ax.imshow(plt.imread(img), aspect="auto")
            ax.set_axis_off()
            save(fig, pdir, f"phylosynteny_{fam}", dpi)
            plt.close(fig)
            made.append(fam)
    return made


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True,
                    help="含 <Family>/Tree_Synteny_Composite.png 的目录")
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
        print(f"[out] panels    = {len(made)} files in {outdir / 'panels_phylosynteny'}")

    mpath = outdir / "_build" / "phylosynteny_manifest.json"
    mpath.parent.mkdir(parents=True, exist_ok=True)
    mpath.write_text(json.dumps(man, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[out] manifest  = {mpath}")


if __name__ == "__main__":
    main()
