#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_discovery_supp.py — 枸杞病毒组发现的 **补充附图**（附图 S-Disc）。

为什么存在（用户 2026-10-13 晚明确指令）
  > 「A图展开高些，B图(J)也是大些，**LKIM 放附图吧**」
  主图 `Fig_GojiVirome_Discovery` 原为 7 面板（3 列 × 5 行）。面板数过多导致
  A/J 被挤到极小（实测 A 高仅 ~1.6 in、J 宽仅 1.8 in）。用户要求把
  **A/J 放大**、同时把 L/K/I/M 四面板**移入独立附图**。
  → 主图 v20 = A + J + D（三行通栏）；本附图承载 L/K/I/M。

  ⚠ 本附图**不是新数据**：四面板全部**复用主图模块的面板函数**
    （`make_discovery_sci_figs.panel_L_domains / panel_K_family /
      panel_I_pair / panel_M_targets`），只换版式容器 → 口径零漂移。
    （历史教训：同一指标在两处各写一遍实现，迟早口径分叉。）

附图版式（本脚本设计，编号 A–D **重新编号**、与主图字母无对应关系）
  ┌ A  Alignment targets of viral contigs（M，横向条，通栏） ┐
  ├ B  Identity by virus family（K，14 科横箱线，通栏）      ┤
  ├ C  Conserved-domain completeness（L，饼图） │ D  Genome  │
  └                                             │   type + │
                                                │   Viroid  ┘
  ⚠ 附图内部字母**独立重编**（A/B/C/D）；**主图字母保持跳号稳定**
    （A/D/I/J/K/L/M），两套字母体系不得互引。

口径（与主图同源，勿混读）
  · 数据源 = 与主图**完全相同的 `Data(root)`**（同一 root、同一 445 calls 口径）。
  · L = `cdd_evidence_report.tsv` 的 domain 完整性（n=291），**不是** CheckV 基因组完整性；
  · K = `All_plant.viruses_info.tsv` 的 nt identity 按科分组（n≥3 的科）；
  · I = `hq_genome_type_summary.tsv` + `Viroid.species_info.tsv`（合计 895 基因组型，
    与主图 A 的 "HQ vOTUs = 200" **不是一回事**）；
  · M = `hmm_ct3_hits.tsv` 的 `ct3_best_viral_desc` 关键词归类（派生量）。

用法
  python fig_discovery_supp.py --root <RNA-Lycium_barbarum_out> [--outdir DIR]
  （通常由 `make_discovery_sci_figs.py` 的 main() 自动调用，无需手工跑）

⚠ 踩过的坑（沿用主图纪律）
  · 复用面板函数 → **必须提前 `apply_sci_style`**，否则吃进程级 rcParams
    （150 dpi/Type3 vs 300 dpi/Type42，不可复现缺陷）；
  · 退化守卫：任一面板数据算不出 → 该面板画 `_empty` 占位，不产半张空图；
  · 面板字母用 `ax.transAxes`，避免 `bbox_inches="tight"` 下的像素偏移。
"""

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
import importlib.util
import json
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]                      # reports/virome_analysis_pipeline -> <repo>


def _load_main():
    """按显式路径加载主图模块（**不要用 import**：文件名含下划线但非包）。"""
    spec = importlib.util.spec_from_file_location(
        "make_discovery_sci_figs", _HERE / "make_discovery_sci_figs.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build(root: Path, outdir: Path, dpi: int = 300,
          virus_name: str = "Lycium barbarum",
          stem: str = "Fig_GojiVirome_Discovery_Supp"):
    m = _load_main()
    m.apply_sci_style(base_size=7.5, grid=False)

    D = m.Data(root, virus_name)
    D._outdir = outdir

    # 附图中每个面板的**字母**独立重编（A–D），与主图字母无对应
    #   A = M 靶标（通栏）  B = K 科箱线（通栏）  C = L 饼图  D = I 双轴
    CELLS = [
        ("A", 0, 0, 1, 2, m.panel_M_targets),
        ("B", 1, 0, 1, 2, m.panel_K_family),
        ("C", 2, 0, 1, 1, m.panel_L_domains),
        ("D", 2, 1, 1, 1, m.panel_I_pair),
    ]
    # 行权：M(8 横条) 与 K(14 科) 都给足；底行 L/I 是矮面板
    ROW_RATIOS = [1.10, 1.30, 0.95]
    FIG = (7.6, 7.5)

    fig = plt.figure(figsize=FIG)
    gs = GridSpec(3, 2, figure=fig, hspace=0.44, wspace=0.30,
                  height_ratios=ROW_RATIOS,
                  left=0.135, right=0.985, top=0.905, bottom=0.045)
    made = []
    for letter, r, c, rs, cs, fn in CELLS:
        ax = fig.add_subplot(gs[r:r + rs, c:c + cs])
        try:
            fn(ax, D)
        except Exception as e:                        # noqa: BLE001
            traceback.print_exc()
            m._empty(ax, f"panel {letter} failed: {type(e).__name__}")
        # panel_I_pair 走 axis("off")+inset → has_data() 恒 False，须显式补字母
        if ax.has_data() or ax.texts or letter == "D":
            if letter == "A":
                # 通栏轴占满画布 → 偏移必须很小（≈ 0），dy 抬到标题上方
                m.panel_letter(ax, letter, dx=0.0, dy=1.10)
            elif letter == "B":
                m.panel_letter(ax, letter, dx=-0.045, dy=1.05)
            elif letter == "C":
                # 饼图轴占满左缘 → 偏移接近 0，dy 抬到标题上方
                m.panel_letter(ax, letter, dx=-0.005, dy=1.14)
            else:
                # I 走 axis("off")+inset，轴左缘在画布 ~2% 处 → 偏移必须极小
                m.panel_letter(ax, letter, dx=0.0, dy=1.06)
        made.append(letter)

    fig.suptitle(
        f"Supplementary \u00b7 {virus_name} virome discovery: evidence quality, "
        f"genome composition and alignment targets",
        x=0.012, y=0.985, ha="left", fontsize=9.0, fontweight="bold")
    fig.text(0.012, 0.952,
             "Companion to Fig_GojiVirome_Discovery (main figure: A funnel / J identity-vs-coverage "
             "/ D novelty sankey). Panels here were split off the main figure to give A and J more room.",
             ha="left", fontsize=5.4, color="#555555")

    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    fig.savefig(pdf, dpi=dpi, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    plt.close(fig)
    print(f"OK  {png}")
    print(f"OK  {pdf}")

    man = {
        "fig": stem,
        "root": str(root),
        "virus_name": virus_name,
        "note": "Supplementary figure split off Fig_GojiVirome_Discovery (v20). "
                "Panels reuse the main module's functions (zero口径 drift). "
                "⚠ Appendix letters A-D are renumbered independently of the main "
                "figure letters (A/D/I/J/K/L/M).",
        "panels": {
            "A": "alignment targets of viral contigs (main-fig M)",
            "B": "nt identity by virus family (main-fig K)",
            "C": "conserved-domain completeness (main-fig L)",
            "D": "genome type + viroid species (main-fig I)",
        },
    }
    (outdir / f"{stem}.manifest.json").write_text(
        json.dumps(man, ensure_ascii=False, indent=2), encoding="utf-8")
    return png


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="discovery 输出根目录（含 10_Reports/）")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--virus-name", default="Lycium barbarum")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--stem", default="Fig_GojiVirome_Discovery_Supp")
    a = ap.parse_args()

    root = Path(a.root).expanduser().resolve()
    if not root.is_dir():
        sys.exit(f"[fatal] --root not a directory: {root}")
    outdir = Path(a.outdir).expanduser().resolve() if a.outdir \
        else (root / "10_Reports" / "S10_SCI_Discovery")
    build(root, outdir, a.dpi, a.virus_name, a.stem)


if __name__ == "__main__":
    main()
