#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_novelty_sankey_supp.py — 09b **附图**：新颖性分类的**完整明细**（附图 S-09b）。

为什么存在（用户 2026-10-13 明确指令）
  > 「同属的新病毒只展示一个，这样避免新病毒过多引起桑基图不好看，**详细的结果放附图**」
  主图桑基（`Fig_GojiVirome_NoveltySankey`）为可读性，Novel 每属只展示 **1 条代表**
  （取发生率最高）→ 85 条 Novel 里 **64 条被折叠**（未进桑基几何）。
  信息零丢失：全部 85 条的完整明细已由数据层写到
  `novelty_sankey_data/novel_dedup_full.tsv`（15 列，含 `shown_in_sankey` 标记）。
  **本脚本把这份明细可视化** = 用户要的「附图」。

附图版式（本脚本设计）
  ┌ A 汇总计数条 ─────────────────────────────────────────────┐
  │  Novel 全量 85 = 桑基展示 21 + 折叠 64（按 category/tier 分层）│
  ├ B 属级折叠条（横向堆叠） ──────────────────────────────────┤
  │  每属：展示 vs 折叠（Potyvirus 1+23、Betacytorhabdovirus 1+13…）│
  ├ C 折叠明细点阵（发生率 vs 长度，标注被折叠者） ──────────────┤
  │  85 点；实心=展示 / 空心=折叠；颜色按 category              │
  └ D 折叠 contig 清单（按属分组，前 N 条） ────────────────────┘

口径
  · 数据源唯一 = `novel_dedup_full.tsv`（数据层产物，勿另起炉灶）
  · 与主图桑基**共用同一份口径**：`shown_in_sankey` 由 `make_novelty_sankey_data.py`
    在**最终保留集**上打出，本脚本只读不算（防两处口径漂移）。
  · 与主图/面板 A/I 是不同口径，不得混读（A = 445 calls、I = 1,063 植物宿主池）。
  · 图中文本一律 ASCII（rcParams 无 CJK 字体）。

⚠ 本脚本踩过的坑（沿用主图纪律）
  · 不用 `TextPath` 量宽（本机字体度量与服务器差 ~12–14%）→ 只用 `_text_w_in` 保守估宽；
  · 面板字母/标题用 `ax.transAxes`，避免 `bbox_inches="tight"` 下的像素偏移；
  · 空图守卫：明细为空直接退出且不产图。
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
import csv
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "reports" / "virome_analysis_pipeline" / "novelty_sankey_data"

# ── 调色板：按显式文件路径加载（项目纪律，勿改成 import）──
_pal = ROOT / "virome_discovery_pipeline" / "utils" / "sci_palette.py"
_spec = importlib.util.spec_from_file_location("sci_palette", _pal)
sp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sp)

# 语义色：展示 vs 折叠（与主图 tier 语义一致：Known 蓝 / Novel 橙）
C_SHOWN = sp.SEMANTIC["genome"]["DNA"]      # #0072B2（深蓝）
C_FOLDED = "#E8A33D"                        # 折叠（浅橙，与 Novel 橙同族但更淡）
C_FOLDED2 = "#F4C98A"
# category 色（三分类）——取自 Okabe-Ito 固定序，保证色盲友好
# ⚠ 图中文本一律 ASCII（rcParams 无 CJK 字体，否则静默吞字并刷 warning）
CAT_COLOR = {
    "新种候选": sp.SEMANTIC.get("verdict", {}).get("KEEP", "#009E73"),
    "远缘病毒候选": "#CC79A7",
    "远缘(未定)": "#999999",
}
CAT_EN = {
    "新种候选": "candidate new species",
    "远缘病毒候选": "divergent virus candidate",
    "远缘(未定)": "divergent (unassigned)",
}
CAT_ORDER = ["新种候选", "远缘病毒候选", "远缘(未定)"]

FS_TITLE = 7.4
FS_SUB = 5.6
FS_TICK = 5.2
FS_NODE = 5.0
FS_NOTE = 4.6


def _load_detail(path: Path):
    if not path.is_file():
        sys.exit(f"❌ 明细数据不存在：{path}\n"
                 f"   先跑 make_novelty_sankey_data.py（会写 novel_dedup_full.tsv）")
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"❌ 明细数据为空：{path}")
    for r in rows:
        r["_shown"] = (str(r.get("shown_in_sankey", "")).strip() == "True")
        for k in ("length", "nt_pident", "nt_qcovs", "n_samples"):
            try:
                r["_" + k] = float(r[k]) if r.get(k) not in (None, "", "NA") else None
            except (TypeError, ValueError):
                r["_" + k] = None
        try:
            r["_prev"] = float(r["prevalence_pct"]) if r.get("prevalence_pct") not in (None, "") else None
        except (TypeError, ValueError):
            r["_prev"] = None
    return rows


def panel_A_summary(ax, rows):
    """A 汇总：Novel 全量 / 展示 / 折叠（按 category 分层横向堆叠）。"""
    n_all = len(rows)
    n_shown = sum(1 for r in rows if r["_shown"])
    n_fold = n_all - n_shown

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # 总量条
    ax.barh(0.80, 1.0, height=0.20, left=0.0, color=C_FOLDED2,
            edgecolor="white", linewidth=0.6)
    ax.barh(0.80, n_shown / n_all, height=0.20, left=0.0, color=C_SHOWN,
            edgecolor="white", linewidth=0.6)
    ax.text(0.012, 0.80, f"shown in sankey  n={n_shown}", va="center", ha="left",
            fontsize=FS_NODE, color="white", fontweight="bold")
    ax.text(1.0 - 0.012, 0.80, f"folded (detail-only)  n={n_fold}",
            va="center", ha="right", fontsize=FS_NODE, color="#7A4A06", fontweight="bold")

    # 按 category 分层
    y = 0.52
    cats = [c for c in CAT_ORDER
            if any(r["category"] == c for r in rows)]
    for c in cats:
        sub = [r for r in rows if r["category"] == c]
        s = sum(1 for r in sub if r["_shown"])
        f = len(sub) - s
        col = CAT_COLOR.get(c, "#999999")
        ax.barh(y, len(sub) / n_all, height=0.15, left=0.0, color=col,
                alpha=0.30, edgecolor=col, linewidth=0.5)
        ax.barh(y, s / n_all, height=0.15, left=0.0, color=col,
                edgecolor="white", linewidth=0.5)
        ax.text(-0.008, y, CAT_EN.get(c, c), va="center", ha="right",
                fontsize=FS_TICK, color="#333333")
        # ⚠⚠ 类别计数标签**不能放在 `1.0 + ε`**：panel A 轴占满画布左缘、右缘即画布右缘，
        #   任何 `x > 1.0` 的文本都会越出画布（体检判 `TEXT-OFF-CANVAS`，实测 3 例：
        #   '11 + 10' / '7 + 17' / '3 + 37'，bbox x1=779~784 > canvas 760）。
        #   → 首选**紧跟条右端**（`bar_end + 0.012`）；若会越界（> 0.965）则**回落到条内右端**，
        #     并加白色描边保证在任何条色上都可读。
        _bar_end = len(sub) / n_all
        _lbl_x = _bar_end + 0.012
        if _lbl_x > 0.965:
            ax.text(0.965, y, f"{s} + {f}", va="center", ha="right",
                    fontsize=FS_TICK, color=col, fontweight="bold",
                    path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])
        else:
            ax.text(_lbl_x, y, f"{s} + {f}", va="center", ha="left",
                    fontsize=FS_TICK, color=col, fontweight="bold")
        y -= 0.155

    ax.text(0.0, 0.06,
            "left segment = shown (unique representative per genus)   "
            "right = folded (kept in this table, omitted from the sankey)",
            transform=ax.transAxes, fontsize=FS_NOTE, color="#555555", ha="left")
    ax.set_title("A  Novel contigs \u00b7 shown vs folded", fontsize=FS_TITLE,
                 fontweight="bold", loc="left")


def panel_B_genus(ax, rows):
    """B 属级：每属展示 vs 折叠（横向堆叠条，按总量降序）。"""
    g = defaultdict(lambda: [0, 0])
    for r in rows:
        g[r["genus"]][0 if r["_shown"] else 1] += 1
    items = sorted(g.items(), key=lambda kv: (-(kv[1][0] + kv[1][1]), kv[0]))

    ax.set_xlim(0, max((s + f) for _, (s, f) in items) * 1.18)
    ax.set_ylim(-0.8, len(items) - 0.2)
    ax.axis("off")

    for i, (gen, (s, f)) in enumerate(items):
        ax.barh(i, s, left=0, height=0.62, color=C_SHOWN,
                edgecolor="white", linewidth=0.5)
        ax.barh(i, f, left=s, height=0.62, color=C_FOLDED,
                edgecolor="white", linewidth=0.5)
        lab = gen if gen != "(Unassigned)" else "Unassigned"
        ax.text(-0.4, i, lab, va="center", ha="right", fontsize=FS_TICK,
                color="#333333")
        ax.text(s + f + ax.get_xlim()[1] * 0.012, i, f"{s}+{f}",
                va="center", ha="left", fontsize=FS_TICK, color="#555555")

    ax.set_title("B  Per-genus: shown + folded", fontsize=FS_TITLE,
                 fontweight="bold", loc="left")
    ax.plot([0, 0], [-0.6, len(items) - 0.4], color="#CCCCCC", lw=0.5, zorder=0)


def panel_C_scatter(ax, rows):
    """C 折叠明细点阵：发生率(x) vs 长度(y)，实心=展示 / 空心=折叠。"""
    pts = [r for r in rows if r["_prev"] is not None and r["_length"] is not None]
    if not pts:
        ax.text(0.5, 0.5, "no length/prevalence data", ha="center", va="center",
                transform=ax.transAxes, color="#999999", fontsize=FS_SUB)
        ax.axis("off")
        ax.set_title("C  Folded detail (prevalence vs length)", fontsize=FS_TITLE,
                     fontweight="bold", loc="left")
        return

    for r in pts:
        col = CAT_COLOR.get(r["category"], "#999999")
        if r["_shown"]:
            ax.scatter(r["_prev"], r["_length"], s=16, facecolor=col,
                       edgecolor="white", linewidth=0.35, zorder=3)
        else:
            ax.scatter(r["_prev"], r["_length"], s=13, facecolor="none",
                       edgecolor=col, linewidth=0.55, zorder=2, alpha=0.9)

    ax.set_xlabel("prevalence (% of 312 samples)", fontsize=FS_TICK)
    ax.set_ylabel("contig length (nt)", fontsize=FS_TICK)
    ax.tick_params(labelsize=FS_TICK - 0.4, length=2)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
    ax.grid(True, color="#EEEEEE", lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    # ⚠ 主体长度都在 600–3600 nt，少数几个 15–17 k nt 的离群值会把其余 80 点
    #   压成一条线 → 纵轴截到 6 k，离群值单独标注条数（不丢信息）。
    _ys = sorted(r["_length"] for r in pts)
    _ylim = 6000
    _n_off = sum(1 for v in _ys if v > _ylim)
    ax.set_ylim(0, _ylim)
    if _n_off:
        ax.text(0.995, 0.02, f"+ {_n_off} contigs >6,000 nt (off scale)",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=FS_NOTE - 0.6, color="#999999")
    ax.set_title("C  Folded detail (prevalence vs length)", fontsize=FS_TITLE,
                 fontweight="bold", loc="left")

    # 图例
    h = [Line2D([], [], marker="o", ls="", mfc=C_SHOWN, mec="white", ms=3.4,
                label="shown"),
         Line2D([], [], marker="o", ls="", mfc="none", mec="#777777", ms=3.4,
                label="folded")]
    for c in CAT_ORDER:
        if any(r["category"] == c for r in pts):
            h.append(Line2D([], [], marker="s", ls="", mfc=CAT_COLOR[c],
                            mec="none", ms=3.4, label=CAT_EN.get(c, c)))
    ax.legend(handles=h, fontsize=FS_TICK - 0.6, frameon=False,
              loc="upper right", handletextpad=0.4, labelspacing=0.28,
              borderpad=0.2)


def panel_D_list(ax, rows, topn=24):
    """D 折叠 contig 清单（按属分组，组内按发生率降序；只列前 topn 条）。

    ⚠ 行高按**面板实际英寸高**换算（FS_NOTE 5 pt × 1.5 行距 ≈ 0.104 in），
      否则在 2.4 权重行里 22 行会互相压字（实测首版重叠）。
    """
    folded = [r for r in rows if not r["_shown"]]
    folded.sort(key=lambda r: (r["genus"], -(r["_prev"] or 0.0), r["contig_id"]))
    show = folded[:topn]

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title(f"D  Folded contigs (top {len(show)} of {len(folded)}, "
                 f"by genus / prevalence)", fontsize=FS_TITLE,
                 fontweight="bold", loc="left")

    if not show:
        ax.text(0.5, 0.5, "nothing folded", ha="center", va="center",
                transform=ax.transAxes, color="#999999", fontsize=FS_SUB)
        return

    n = len(show)
    # 面板在该 gridspec 下约 2.0 in 高；表头占 top 8%，行高均分剩余
    y_top, y_bot = 0.94, 0.015
    dy = (y_top - y_bot) / max(n, 1)
    # 表头
    hdr = [("genus", 0.000), ("species", 0.175), ("type", 0.435),
           ("len", 0.500), ("prev%", 0.598), ("contig_id", 0.700)]
    for txt, x in hdr:
        ax.text(x, 1.0, txt, fontsize=FS_NOTE - 0.5, color="#777777",
                va="top", ha="left", fontweight="bold")
    # 属分组分隔线 + 组名只写一次
    last_gen = None
    for i, r in enumerate(show):
        yy = y_top - (i + 0.5) * dy
        col = CAT_COLOR.get(r["category"], "#999999")
        gen = r["genus"] if r["genus"] != "(Unassigned)" else "Unassigned"
        if gen != last_gen:
            ax.text(0.000, yy, gen, fontsize=FS_NOTE - 0.5, color="#333333",
                    va="center", ha="left", fontweight="bold")
            last_gen = gen
        ax.text(0.175, yy, r["species"], fontsize=FS_NOTE - 0.5, color=col, va="center")
        ax.text(0.435, yy, r["genome_type"], fontsize=FS_NOTE - 0.5,
                color="#333333", va="center")
        ln = f"{int(r['_length'])}" if r["_length"] else "-"
        ax.text(0.500, yy, ln, fontsize=FS_NOTE - 0.5, color="#333333", va="center")
        pv = f"{(r['_prev'] or 0):.1f}" if r["_prev"] is not None else "-"
        ax.text(0.598, yy, pv, fontsize=FS_NOTE - 0.5, color="#333333", va="center")
        cid = r["contig_id"]
        cid = cid if len(cid) <= 40 else cid[:37] + "..."
        ax.text(0.700, yy, cid, fontsize=FS_NOTE - 1.0, color="#888888", va="center",
                family="monospace")


def build(detail_tsv: Path, outdir: Path, stem="Fig_GojiVirome_NoveltySankey_Supp"):
    rows = _load_detail(detail_tsv)
    n_all = len(rows)
    n_shown = sum(1 for r in rows if r["_shown"])
    print(f"   明细 {n_all} 条（展示 {n_shown} / 折叠 {n_all - n_shown}）")

    plt.rcParams.update(sp.sci_rcparams(base_size=8))
    fig = plt.figure(figsize=(7.6, 9.2))

    gs = fig.add_gridspec(
        4, 1,
        height_ratios=[0.95, 2.5, 2.3, 2.6],
        hspace=0.42, left=0.16, right=0.985, top=0.915, bottom=0.045)

    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])
    axC = fig.add_subplot(gs[2])
    axD = fig.add_subplot(gs[3])

    panel_A_summary(axA, rows)
    panel_B_genus(axB, rows)
    panel_C_scatter(axC, rows)
    panel_D_list(axD, rows)

    fig.suptitle(
        "Supplementary \u00b7 Goji virome novelty taxonomy: complete Novel-contig detail",
        x=0.012, y=0.985, ha="left", fontsize=9.4, fontweight="bold")
    fig.text(0.012, 0.958,
             f"All {n_all} Novel contigs from the 09b pool. The sankey (main figure) shows only "
             f"{n_shown} genus representatives (highest prevalence per genus); "
             f"the remaining {n_all - n_shown} are listed here.",
             ha="left", fontsize=FS_SUB, color="#555555")

    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"OK  {png}")
    print(f"OK  {pdf}")

    # manifest（供主图/体检脚本消费，防口径漂移）
    man = {
        "fig": stem,
        "data_source": str(detail_tsv.relative_to(ROOT)).replace("\\", "/"),
        "n_novel_total": n_all,
        "n_shown_in_sankey": n_shown,
        "n_folded": n_all - n_shown,
        "note": "Companion detail figure for Fig_GojiVirome_NoveltySankey. "
                "Sankey shows 1 representative per genus; this figure lists all.",
        "panels": {"A": "summary shown/folded by category",
                   "B": "per-genus shown+folded",
                   "C": "prevalence vs length scatter",
                   "D": "folded contig table"},
    }
    (outdir / f"{stem}.manifest.json").write_text(
        json.dumps(man, ensure_ascii=False, indent=2), encoding="utf-8")
    return png


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detail-tsv", default=str(DATA / "novel_dedup_full.tsv"))
    ap.add_argument("--outdir", default=str(
        ROOT / "reports" / "_discovery_sci_20261012" / "_lt" / "10_Reports" / "S10_SCI_Discovery"))
    ap.add_argument("--stem", default="Fig_GojiVirome_NoveltySankey_Supp")
    args = ap.parse_args()
    build(Path(args.detail_tsv), Path(args.outdir), args.stem)


if __name__ == "__main__":
    main()
