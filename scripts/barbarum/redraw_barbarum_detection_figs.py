#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""redraw_barbarum_detection_figs.py — 宁夏枸杞 (Lycium barbarum) 01_detection 检出图重绘

背景:
  known_virus_barbarum/01_detection/summary/Fig1_study_design.{png,pdf} 与
  Fig2_prevalence.{png,pdf} 原为指向 known_virus_all_v2 的 symlink, 内容是
  跨物种混合队列 (Total 1,156 / 820 Lycium RNA / 含 DNA metagenome 313),
  与单物种 barbarum 不符。原生成脚本已不存在 (MMPV-RNA 全库 + git 全历史零命中)。

本脚本按 L. barbarum 单物种口径重绘:
  Fig1 (a) 研究流程示意 (b) 队列构成与检出状态
  Fig2 (a) 病毒流行率 (b) 组织嗜性 (c) 地理分布

口径 (与 analysis_report.txt 一致):
  分母 = 445 (novel_virus/01_Assembly 权威 run 数)
  检出 = best.summary 中的 run (262, 58.9%)
  元数据 = Global_Metadata_Final_with_Coordinates.tsv

配色 (2026-10-12 修订):
  **一律走 `virome_discovery_pipeline/utils/sci_palette.py`**（全项目唯一真源）。
  旧版本在文件内自带一份 `OKABE_ITO`，且**顺序与真源不同**
  （旧: `#0072B2,#E69F00,...`；真源: `#E69F00,#56B4E9,...`）→ 跨图串色。
  已改为按**显式文件路径**加载真源（服务器 site-packages 里有同名第三方 `utils`
  包会盖住命名空间包，`from utils.sci_palette import` 会失败）。

用法:
  python redraw_barbarum_detection_figs.py \
     -v.../01_detection/summary/all_viruses.best.summary.tsv \
     -m .../Global_Metadata_Final_with_Coordinates.tsv \
     -d .../01_detection/summary \
     --assembly-root .../02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly
"""
import argparse
import collections
import csv
import importlib.util as _ilu
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

# ── 统一配色（唯一真源，按显式路径加载） ────────────────────────────
_REPO = Path(__file__).resolve().parents[2]      # MMPV-RNA/


def _load_palette():
    for rel in ("virome_discovery_pipeline/utils/sci_palette.py",
                "virome_phylo_pipeline/utils/sci_palette.py"):
        cand = _REPO / rel
        if cand.is_file():
            spec = _ilu.spec_from_file_location("_sci_palette", cand)
            mod = _ilu.module_from_spec(spec)
            sys.modules["_sci_palette"] = mod
            spec.loader.exec_module(mod)
            return mod
    raise ImportError(f"sci_palette.py not found under {_REPO}")


_P = _load_palette()
OKABE_ITO = _P.OKABE_ITO
apply_sci_style = _P.apply_sci_style
COL_W_1 = getattr(_P, "COL_W_1", 9.32)   # 通栏宽 (in)
COL_W_2 = getattr(_P, "COL_W_2", 4.58)   # 半栏宽 (in)

# 语义色：从真源 OKABE_ITO 取，禁止游离 hex
C_RNA = OKABE_ITO[4]        # #0072B2 blue   —— RNA 病毒
C_VIROID = OKABE_ITO[0]     # #E69F00 orange —— 类病毒
C_TOTAL = "#DDDDDD"         # 队列总量（中性灰）
C_NEG = "#999999"           # 未检出（深一档灰）
C_INK = "#333333"           # 文字 / 描边
C_ACCENT = OKABE_ITO[4]     # 流程框强调色（与 RNA 病毒同蓝，全图一致）


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    margin = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / denom
    return (max(0.0, centre - margin), min(1.0, centre + margin))


def count_assembly_runs(root):
    p = Path(root)
    if not p.is_dir():
        return None
    return len([d for d in p.iterdir()
                if d.is_dir() and d.name.lower() not in ("logs", "log")])


def short_virus(name):
    return {
        "Cytorhabdovirus sp. 'lycii'": "Cytorhabdovirus sp. 'lycii' (GCVA)",
        "Potato spindle tuber viroid": "Potato spindle tuber viroid (PSTVd)",
        "Tomato chlorotic dwarf viroid": "Tomato chlorotic dwarf viroid (TCDVd)",
        "Citrus exocortis Yucatan viroid": "Citrus exocortis Yucatan viroid (CEYVd)",
        "Grapevine-associated RNA virus 4": "Grapevine-associated RNA virus 4 (GARV-4)",
    }.get(name, name)


def short_code(name):
    """面板 (b) 图例用的短代号（长种名会把图例撑到轴上）。"""
    return {
        "Cytorhabdovirus sp. 'lycii'": "GCVA",
        "Potato spindle tuber viroid": "PSTVd",
        "Tomato chlorotic dwarf viroid": "TCDVd",
        "Citrus exocortis Yucatan viroid": "CEYVd",
        "Grapevine-associated RNA virus 4": "GARV-4",
        "Potato virus H": "PVH",
        "Potato virus M": "PVM",
    }.get(name, name)


def mol_class(name):
    return "Viroid" if "viroid" in name.lower() else "RNA virus"


# --------------------------------------------------------------------------- Fig 1
# ⚠ 画布宽 = **通栏宽 COL_W_1 = 9.32 in**（2026-10-12 修订）。
#   旧版 13.0 in 宽 → 排进通栏时缩放 9.32/13.0 = 0.717，10.5 pt 框标题印出来
#   只剩 7.5 pt；若按半栏排更是 45 折。现画布 = 通栏宽 ⇒ 印刷缩放 1.000，
#   源字号即印刷字号（框标题 8.5 / 副标题 6.5 / 刻度 7.5 pt）。
FIG1_W, FIG1_H = COL_W_1, 5.60
FIG1_TOP, FIG1_BOTTOM = 0.90, 0.12


def fig1_study_design(out_stem, denom, detected_runs, meta, run_meta):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(FIG1_W, FIG1_H),
                                   gridspec_kw={"width_ratios": [1.05, 1]})

    # ---- (a) 工作流 ----
    steps = [
        ("Public data", f"SRA + GSA public RNA-seq\n({denom} L. barbarum runs)"),
        ("QC & host depletion", "fastp · Kraken2 · rRNA removal"),
        ("Known-virus detection",
         f"salmon · PlantVirusDB\nlayered filtering\n({len(detected_runs)} positives)"),
        ("Assembly & variant calling", "contig consensus · SNP/indel"),
        ("Population evolution", "GCVA · PSTVd\ndating / selection / recombination"),
        ("Host response", "TEA-GCN dual networks\ntissue tropism"),
    ]
    N = len(steps)
    STEP = 1.0
    BOX_H = 0.94
    PAD_TOP = 0.06          # 框顶 → 标题顶
    PAD_GAP = 0.03          # 标题底 → 副标题顶
    TITLE_PT, SUB_PT, LS = 8.5, 6.5, 1.40

    ylim = (-0.75, N - 1 + 0.75)
    # ⚠ 文字高度必须换算到**数据坐标**才能保证「不出框」：
    #   1 数据单位 = 轴高(in)×72 / ylim 跨度(pt)。旧版用固定 0.16/0.20 的
    #   数据坐标偏移，3 行副标题实际占 0.54 单位 → 溢出框底 0.12 单位，
    #   压到下一个框的标题上（2026-10-12 实测：'(262 positives)' ×
    #   'Full-length assembly'）。这里改成按字号现算。
    ax_h_in = FIG1_H * (FIG1_TOP - FIG1_BOTTOM)
    pt_per_unit = ax_h_in * 72.0 / (ylim[1] - ylim[0])
    title_line_h = TITLE_PT * 1.20 / pt_per_unit
    sub_line_h = SUB_PT * LS / pt_per_unit

    for i, (title, sub) in enumerate(steps):
        y = N - 1 - i
        n_title = title.count("\n") + 1
        n_sub = sub.count("\n") + 1
        # ⚠ 内容高度必须**整段**塞进框内：标题行数 × 行高 + 间隙 + 副标题行数 × 行高
        content_h = title_line_h * n_title + PAD_GAP + sub_line_h * n_sub
        assert PAD_TOP + content_h < BOX_H, (
            f"step {title!r}: content {content_h:.3f} + pad {PAD_TOP:.2f} "
            f"> box height {BOX_H:.2f} (title {n_title} line(s), sub {n_sub} line(s))")
        axA.add_patch(mpatches.FancyBboxPatch(
            (0.04, y - BOX_H / 2), 0.92, BOX_H,
            boxstyle="round,pad=0.010,rounding_size=0.02",
            facecolor="#EAF2F8" if i % 2 == 0 else "#DCE9F5",
            edgecolor=C_ACCENT, linewidth=0.9))
        t_top = y + BOX_H / 2 - PAD_TOP
        axA.text(0.5, t_top, title, ha="center", va="top",
                 fontsize=TITLE_PT, color=C_ACCENT, fontweight="bold",
                 linespacing=1.20)
        # 副标题顶 = 标题顶 − 标题**全部行**的高度 − 间隙（旧版只减了一行 → 多行标题压副标题）
        s_top = t_top - title_line_h * n_title - PAD_GAP
        axA.text(0.5, s_top, sub, ha="center", va="top",
                 fontsize=SUB_PT, color=C_INK, linespacing=LS)
        # 溢出体检：副标题底不得低于框底
        s_bottom = s_top - sub_line_h * n_sub
        assert s_bottom > y - BOX_H / 2, \
            f"step {title!r}: subtitle overflows box by {(y - BOX_H / 2) - s_bottom:.3f}"
        if y > 0:
            axA.annotate("", xy=(0.5, y - BOX_H / 2 - 0.055),
                         xytext=(0.5, y - BOX_H / 2 - 0.005),
                         arrowprops=dict(arrowstyle="-|>", color=C_ACCENT, lw=1.2))
    axA.set_xlim(0, 1)
    axA.set_ylim(*ylim)
    axA.axis("off")
    axA.set_title("(a) Study workflow", fontsize=9.5, fontweight="bold", loc="left")

    # ---- (b) 队列构成与检出状态 ----
    n_pos = len(detected_runs)
    n_neg = denom - n_pos
    bars = [("Total\n(processed)", denom, C_TOTAL),
            ("Virus-positive", n_pos, C_RNA),
            ("Virus-negative", n_neg, C_NEG)]
    xs = range(len(bars))
    vals = [b[1] for b in bars]
    axB.bar(xs, vals, width=0.58, color=[b[2] for b in bars],
            edgecolor=C_INK, linewidth=0.8)
    for x, v in zip(xs, vals):
        axB.text(x, v + denom * 0.018, f"{v}\n({v/denom*100:.1f}%)",
                 ha="center", va="bottom", fontsize=8.0, fontweight="bold",
                 linespacing=1.4, color=C_INK)

    axB.set_xticks(list(xs))
    axB.set_xticklabels([b[0] for b in bars], fontsize=7.5)
    axB.tick_params(axis="x", pad=5)
    axB.set_ylabel("Number of runs", fontsize=8.0)
    axB.set_ylim(0, denom * 1.22)
    axB.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
    axB.set_title(f"(b) Cohort composition (L. barbarum, n = {denom})",
                  fontsize=9.5, fontweight="bold", loc="left")
    axB.text(0.98, 0.95,
             f"prevalence = {n_pos}/{denom} = {n_pos/denom*100:.1f}%",
             transform=axB.transAxes, ha="right", va="top", fontsize=7.5,
             color=C_INK,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#F5F5F5",
                       edgecolor="#999999", linewidth=0.8))

    # Fig1 (a) 用 axis("off") + 数据坐标 0~1 画流程框, 与 tight_layout 冲突
    # (会报 "margins cannot be made large enough" 并导致 bbox 横向溢出),
    # 因此改用固定 subplots_adjust。
    fig.subplots_adjust(left=0.02, right=0.98, top=FIG1_TOP,
                        bottom=FIG1_BOTTOM, wspace=0.18)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300)
    plt.close(fig)
    return f"{out_stem}.png/.pdf"


# --------------------------------------------------------------------------- Fig 2
def _label_a(nm, width=42):
    """面板 (a) 的 y 轴标签：种名 + 缩写；过长截断（超长名会把绘图区挤没）。"""
    s = short_virus(nm)
    return s if len(s) <= width else s[:width - 1].rstrip() + "…"


def fig2_prevalence(out_stem, by_virus, denom, run_meta, min_prev=0.0, top_n=25):
    # ⚠ 画布宽 = **通栏宽 COL_W_1 = 9.32 in**（2026-10-12 修订）。
    #   旧版 16.5 in 宽 → 排进通栏缩放 9.32/16.5 = 0.565，8.6 pt 的种名印出来
    #   只有 4.9 pt（低于可读下限）。改成 **2 行布局**（(a) 通栏 / (b)(c) 并排）
    #   后各面板宽都 ≈ 通栏宽量级 ⇒ 印刷缩放 1.000，源字号即印刷字号。
    fig = plt.figure(figsize=(COL_W_1, 8.90))
    axA = fig.add_axes([0.265, 0.578, 0.710, 0.358])   # (a) 流行率（通栏）
    axB = fig.add_axes([0.078, 0.070, 0.352, 0.378])   # (b) 组织嗜性
    axC = fig.add_axes([0.560, 0.070, 0.415, 0.378])   # (c) 地理分布

    order = sorted(by_virus.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    order = [(nm, g) for nm, g in order if len(g) / denom >= min_prev][:top_n]

    # ---- (a) 流行率 ----
    names = [_label_a(nm) for nm, _ in order]
    ns = [len(g) for _, g in order]
    rates = [n / denom * 100 for n in ns]
    bar_colors = [C_RNA if mol_class(nm) == "RNA virus" else C_VIROID
                  for nm, _ in order]
    ypos = list(range(len(names)))
    axA.barh(ypos, rates, color=bar_colors, edgecolor=C_INK,
             linewidth=0.5, height=0.70)
    axA.set_yticks(ypos)
    axA.set_yticklabels(names, fontsize=7.2)
    axA.invert_yaxis()
    axA.set_xlabel("Prevalence in L. barbarum cohort (%)", fontsize=8.0)
    axA.set_title(f"(a) Virome prevalence (n = {denom})", fontsize=9.5,
                  fontweight="bold", loc="left")
    for y, r, n in zip(ypos, rates, ns):
        axA.text(r + max(rates) * 0.012 + 0.25, y, f"{r:.1f}%  (n={n})",
                 va="center", fontsize=6.8, color=C_INK)
    axA.set_xlim(0, max(rates) * 1.30)
    axA.xaxis.set_major_locator(MaxNLocator(nbins=6))
    axA.tick_params(axis="x", labelsize=7.5)
    # ⚠ 图例句柄不要复用 barh 的容器（会继承 height=0.72 的奇怪比例）
    axA.legend(handles=[mpatches.Patch(facecolor=C_RNA, edgecolor=C_INK,
                                       linewidth=0.5, label="RNA virus"),
                        mpatches.Patch(facecolor=C_VIROID, edgecolor=C_INK,
                                       linewidth=0.5, label="Viroid")],
               loc="lower right", fontsize=7.5)

    # ---- (b) 组织嗜性 (仅 top4 病毒) ----
    top4 = [(nm, g) for nm, g in order[:4]]
    tissues = ["fruit", "leaf", "flower", "root", "stem"]
    width = 0.8 / max(len(top4), 1)
    xs = list(range(len(tissues)))
    for i, (nm, g) in enumerate(top4):
        inf = {s for s in (r["Sample"] for r in g)}
        row = []
        for t in tissues:
            grp = [s for s in inf
                   if (run_meta.get(s, {}).get("Tissue") or "").strip().lower() == t]
            den = sum(1 for s, m in run_meta.items()
                      if (m.get("Tissue") or "").strip().lower() == t)
            row.append(len(grp) / den * 100 if den else 0.0)
        axB.bar([x + i * width for x in xs], row, width=width * 0.92,
                color=OKABE_ITO[i], edgecolor=C_INK, linewidth=0.5,
                label=short_code(nm))
    axB.set_xticks([x + 0.4 - width / 2 for x in xs])
    axB.set_xticklabels(tissues, fontsize=7.5)
    axB.tick_params(axis="y", labelsize=7.5)
    axB.set_ylabel("Prevalence within tissue (%)", fontsize=8.0)
    axB.set_title("(b) Tissue tropism (top viruses)", fontsize=9.5,
                  fontweight="bold", loc="left")
    # ⚠ 图例移到**左上**：右上会被 flower 列（PSTVd ~100 %）顶穿；
    #   并改用短代号，长种名会把图例撑出轴框。
    axB.legend(fontsize=7.0, loc="upper left", framealpha=0.92,
               borderpad=0.3, labelspacing=0.25)
    axB.set_ylim(0, 118)
    axB.set_yticks([0, 25, 50, 75, 100])

    # ---- (c) 地理分布 ----
    prov = collections.Counter()
    for s, m in run_meta.items():
        loc = (m.get("Location") or "").strip()
        if not loc:
            continue
        parts = [p.strip() for p in loc.split(",")]
        p = parts[1] if len(parts) > 1 and parts[0].lower().startswith("china") \
            else parts[0]
        prov[p] += 1
    top = prov.most_common(8)[::-1]
    axC.barh([t[0] for t in top], [t[1] for t in top],
             color=C_RNA, edgecolor=C_INK, linewidth=0.5, height=0.68)
    for y, (_, v) in enumerate(top):
        axC.text(v + max(t[1] for t in top) * 0.015, y, str(v),
                 va="center", fontsize=7.0, color=C_INK)
    axC.set_xlabel("Runs in metadata", fontsize=8.0)
    axC.set_title("(c) Geographic origin", fontsize=9.5, fontweight="bold", loc="left")
    axC.set_xlim(0, max(t[1] for t in top) * 1.18)
    axC.tick_params(axis="both", labelsize=7.5)

    for ext in ("png", "pdf"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return f"{out_stem}.png/.pdf"


def main():
    ap = argparse.ArgumentParser(description="重绘 L. barbarum 01_detection 检出图")
    ap.add_argument("-v", "--virus-summary", required=True)
    ap.add_argument("-m", "--meta", required=True)
    ap.add_argument("-d", "--out-dir", required=True)
    ap.add_argument("--assembly-root", default=None,
                    help="novel_virus/01_Assembly, 用于统计权威分母")
    ap.add_argument("--denominator", type=int, default=None)
    ap.add_argument("--top-n", type=int, default=25)
    ap.add_argument("--min-prev", type=float, default=0.0)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    vsum, meta_p = Path(a.virus_summary), Path(a.meta)
    for p in (vsum, meta_p):
        if not p.is_file():
            raise SystemExit(f"[FATAL] 不存在: {p}")

    with open(vsum, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise SystemExit("[FATAL] best.summary 为空")

    by_virus = collections.defaultdict(list)
    for r in rows:
        by_virus[(r.get("Adjusted_Species") or "").strip()].append(r)
    detected = sorted({(r.get("Sample") or "").strip() for r in rows} - {""})

    run_meta = {}
    with open(meta_p, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            k = (r.get("Sample") or "").strip()
            if k:
                run_meta[k] = r

    n_asm = count_assembly_runs(a.assembly_root) if a.assembly_root else None
    denom = a.denominator or n_asm or len(detected)
    if denom < len(detected):
        raise SystemExit(f"[FATAL] 分母 {denom} < 检出样本 {len(detected)}，口径矛盾")

    print(f"[data] 分母={denom}  检出样本={len(detected)} ({len(detected)/denom*100:.1f}%)")
    print(f"[data] 病毒={len(by_virus)}  元数据匹配={sum(1 for s in detected if s in run_meta)}/{len(detected)}")

    if a.dry_run:
        print("[dry-run] 未写盘")
        return

    outdir = Path(a.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    apply_sci_style(base_size=8, grid=False)

    made = []
    for stem, fn in (
        ("Fig1_study_design",
         lambda s: fig1_study_design(s, denom, detected, meta_p, run_meta)),
        ("Fig2_prevalence",
         lambda s: fig2_prevalence(s, by_virus, denom, run_meta,
                                   min_prev=a.min_prev, top_n=a.top_n)),
    ):
        target = outdir / stem
        for ext in (".png", ".pdf"):
            t = Path(str(target) + ext)
            if t.is_symlink():
                t.unlink()
                print(f"[guard] 解除 symlink: {t.name}")
        made.append(fn(str(target)))

    for m in made:
        print(f"[OK] {m}")


if __name__ == "__main__":
    main()
