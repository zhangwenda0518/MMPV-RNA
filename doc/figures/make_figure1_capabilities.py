# -*- coding: utf-8 -*-
"""MMPV Figure 1 — 四大集成能力 + 十三阶段发现管道 + 三条设计策略 (论文用)

用法:  python doc/figures/make_figure1_capabilities.py
输出:  doc/figures/MMPV_Figure1.{png,pdf,svg}
画幅:  190x138 mm, Okabe-Ito 配色, Arial
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

mpl.use("Agg")
mpl.rcParams.update({
    "font.family": "Arial",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

W, H = 190.0, 138.0
OUT = Path(__file__).resolve().parent

C1, C2, C3, C4 = "#0072B2", "#E69F00", "#009E73", "#CC79A7"   # 四能力
INK, SUB, LINE = "#222222", "#4A4A4A", "#C9D2DA"
PILL_FILL, PILL_EC = "#FDF2E2", "#E3B36F"

fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.set_aspect("equal"); ax.axis("off")


def rbox(x, y, w, h, fc, ec, lw=0.8, r=1.4, z=3):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle=f"round,pad=0,rounding_size={r}",
                 fc=fc, ec=ec, lw=lw, zorder=z))


def arrow(p0, p1, color="#6B7A88", lw=1.2, ms=8, dashed=False, z=5):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms,
                 color=color, lw=lw, linestyle=(0, (2.2, 1.6)) if dashed else "-",
                 shrinkA=0, shrinkB=0, zorder=z))


# ══════════════ 输入条 ══════════════
ax.text(W / 2, 135.2,
        "INPUT — Public archives: NGDC GSA · NCBI SRA  |  self-generated data",
        fontsize=5.8, color="#8A97A3", ha="center", va="center",
        fontweight="bold")
arrow((31, 133.4), (31, 130.6), color="#8A97A3", lw=1.0, ms=6)

# ══════════════ 能力 1 (左上) ══════════════
rbox(4, 88, 54, 42, "white", C1, lw=1.0)
ax.add_patch(plt.Rectangle((4, 128.2), 54, 1.8, fc=C1, ec="none", zorder=4))
ax.text(6.2, 126.6, "1 · Public Data Acquisition", fontsize=6.9,
        color=C1, fontweight="bold", ha="left", va="center", zorder=5)
ax.text(6.2, 123.4, "public_metadata_pipeline", fontsize=4.9, color="#9AA6B2",
        style="italic", family="monospace", ha="left", va="center", zorder=5)
c1_items = ["NGDC GSA + NCBI SRA retrieval", "AI-driven metadata curation",
            "Literature-traced completion", "Parallel download (aria2c · prefetch)",
            "6 SCI-grade landscape plots"]
fy = 119.2
for it in c1_items:
    ax.text(6.4, fy, "·", fontsize=6.5, color=C1, ha="left", va="center",
            zorder=5, fontweight="bold")
    ax.text(8.4, fy, it, fontsize=5.8, color=SUB, ha="left", va="center", zorder=5)
    fy -= 4.5
rbox(6.2, 90.2, 49.6, 6.2, "#EAF2F8", "#9DC3DE", lw=0.7, r=1.0, z=4)
ax.text(31, 93.3, "Desktop GUI — manual editing & visualization",
        fontsize=5.2, color="#2E5E8C", ha="center", va="center",
        zorder=5, fontweight="bold")

# ══════════════ 能力 2 (右上, 主角) ══════════════
rbox(61, 88, 125, 42, "white", C2, lw=1.2)
ax.add_patch(plt.Rectangle((61, 128.2), 125, 1.8, fc=C2, ec="none", zorder=4))
ax.text(63, 126.6, "2 · Benchmark-Driven Discovery Pipeline", fontsize=6.9,
        color="#9C6500", fontweight="bold", ha="left", va="center", zorder=5)
ax.text(63, 123.4, "virome_discovery_pipeline", fontsize=4.9, color="#9AA6B2",
        style="italic", family="monospace", ha="left", va="center", zorder=5)
ax.text(63.5, 120.2, "13 stages · every default parameter benchmark-validated",
        fontsize=5.4, color=INK, ha="left", va="center", style="italic", zorder=5)

stages = ["Quality\ncontrol", "Host\nfiltering", "Depth\nnormalization",
          "De novo\nassembly", "Multi-tool\nidentification", "UniProt-strict\nfiltering",
          "COBRA\nextension", "Flye\nco-assembly",
          "Ref-guided\nclustering", "Multi-tool\ntaxonomy", "3-tier host\nprediction",
          "CheckV\ncompleteness", "4-branch\ncascade rescue"]
OUT_STAGE = "Integrated\nreport"
PX, PW, PG, PH = 65.0, 14.5, 2.6, 6.5
ROW1_Y, ROW2_Y = 110.0, 100.0


def pill(cx, cy, text, output=False):
    rbox(cx, cy, PW, PH, "#E8F5EE" if output else PILL_FILL,
         "#5BAF7D" if output else PILL_EC, lw=0.8, r=1.0, z=5)
    ax.text(cx + PW / 2, cy + PH / 2, text, fontsize=4.6,
            color="#1F5C3D" if output else INK, ha="center", va="center",
            zorder=6, linespacing=1.05)


for j in range(7):                                   # row1: 左→右 (stage 1-7)
    x = PX + j * (PW + PG)
    pill(x, ROW1_Y, stages[j])
    if j < 6:
        arrow((x + PW + 0.35, ROW1_Y + PH / 2), (x + PW + PG - 0.35, ROW1_Y + PH / 2),
              color="#B08640", lw=1.0, ms=5)
arrow((PX + 6 * (PW + PG) + PW / 2, ROW1_Y),                # 换行: 7→8
      (PX + 6 * (PW + PG) + PW / 2, ROW2_Y + PH), color="#B08640", lw=1.0, ms=5)
for j in range(6, -1, -1):                           # row2: 右→左 (stage 8-13 + 产出)
    x = PX + j * (PW + PG)
    idx = 13 - j
    pill(x, ROW2_Y, stages[idx] if idx < 13 else OUT_STAGE, output=(idx >= 13))
    if j > 0:                                        # 从右列指向左列
        arrow((x - 0.35, ROW2_Y + PH / 2), (x - PG + 0.35, ROW2_Y + PH / 2),
              color="#B08640", lw=1.0, ms=5)

bench = ["6-D benchmark · 20+ tools", "1,536 identification combinations",
         "6 assemblers × 54 depth × background conditions"]
bx = 66.5
for t in bench:
    tw = 2.3 + 0.86 * len(t)
    rbox(bx, 91.0, tw, 5.0, "#FDF2E2", "#E3B36F", lw=0.7, r=1.0, z=5)
    ax.text(bx + tw / 2, 93.5, t, fontsize=4.9, color="#7A5200",
            ha="center", va="center", zorder=6, fontweight="bold")
    bx += tw + 3.0

# ══════════════ 能力 3 / 4 (中排) ══════════════
rbox(4, 56, 88, 26, "white", C3, lw=1.0)
ax.add_patch(plt.Rectangle((4, 80.2), 88, 1.8, fc=C3, ec="none", zorder=4))
ax.text(6.2, 78.4, "3 · Rapid Known-Virus Analysis", fontsize=6.9, color=C3,
        fontweight="bold", ha="left", va="center", zorder=5)
ax.text(6.2, 75.2, "virome_analysis_pipeline", fontsize=4.9, color="#9AA6B2",
        style="italic", family="monospace", ha="left", va="center", zorder=5)
c3L = ["Immediate detection", "Quantification (TPM · RPM · depth)",
       "Variant detection"]
c3R = ["Evolutionary characterization", "Co-infection & association",
       "Pan-virome landscape"]
for col, items in ((6.4, c3L), (50.0, c3R)):
    fy = 70.6
    for it in items:
        ax.text(col, fy, "·", fontsize=6.5, color=C3, ha="left", va="center",
                zorder=5, fontweight="bold")
        ax.text(col + 2.0, fy, it, fontsize=5.7, color=SUB, ha="left",
                va="center", zorder=5)
        fy -= 4.5

rbox(95, 56, 91, 26, "white", C4, lw=1.0)
ax.add_patch(plt.Rectangle((95, 80.2), 91, 1.8, fc=C4, ec="none", zorder=4))
ax.text(97.2, 78.4, "4 · Semi-automatic GenBank Submission", fontsize=6.9,
        color=C4, fontweight="bold", ha="left", va="center", zorder=5)
ax.text(97.2, 75.2, "virome_submission_pipeline", fontsize=4.9, color="#9AA6B2",
        style="italic", family="monospace", ha="left", va="center", zorder=5)
c4_items = ["Topology triage · MIUVIG metadata templates",
            "Sequin / suvtk packaging → .sqn"]
fy = 70.6
for it in c4_items:
    ax.text(97.4, fy, "·", fontsize=6.5, color=C4, ha="left", va="center",
            zorder=5, fontweight="bold")
    ax.text(99.4, fy, it, fontsize=5.7, color=SUB, ha="left", va="center", zorder=5)
    fy -= 4.5
rbox(99.4, 58.2, 62, 5.6, "#F7EDF4", "#D39EC2", lw=0.7, r=1.0, z=4)
ax.text(130.4, 61.0, "Dedicated desktop GUI", fontsize=5.2, color="#7C4A6B",
        ha="center", va="center", zorder=5, fontweight="bold")

# ══════════════ 主流程箭头 ══════════════
arrow((58, 109), (61, 109))                    # 1 → 2
arrow((70, 88), (70, 82.4))                    # 2 → 3
arrow((92, 69), (95, 69))                      # 3 → 4

# ══════════════ 三条设计策略 (底层) ══════════════
ax.add_patch(FancyBboxPatch((4, 4), 182, 40,
             boxstyle="round,pad=0,rounding_size=2.0",
             fc="#F6F7F8", ec=LINE, lw=0.7, zorder=1))
ax.text(6.6, 42.2, "CROSS-CUTTING DESIGN STRATEGIES", fontsize=5.8,
        color="#8A97A3", ha="left", va="top", fontweight="bold", zorder=2)

strategies = [
    (C1, "Data Reduction", "99.8%", "host-filtering efficiency",
     ["Reference-guided pre-clustering", "Host-based sequence partitioning"]),
    (C2, "Quality Enhancement", "4-branch", "cascade genome rescue",
     ["COBRA resolves de Bruijn breakpoints", "Flye co-assembly merges cross-sample",
      "BLASTN → PlantVirusDB rebuild", "bwa-mem2 + viral_consensus / ragtag",
      "Genus-length rescue as final fallback"]),
    (C3, "Reliability Enhancement", "F1 = 0.942", "identification (UniProt-strict)",
     ["UniProt-strict ID ≈ 5-tool ensemble", "Family accuracy 94.1% (+72.3 pp)",
      "3-tier host + plant-specific filter"]),
]
SX, SW, SG, SY, SH = 6.5, 57.7, 3.0, 6.0, 32.0
for i, (acc, name, stat, statlab, lines) in enumerate(strategies):
    x = SX + i * (SW + SG)
    rbox(x, SY, SW, SH, "white", "#D8DCE0", lw=0.7)
    ax.add_patch(plt.Rectangle((x, SY), 1.6, SH, fc=acc, ec="none", zorder=4))
    ax.text(x + 3.4, SY + SH - 3.2, name, fontsize=6.2, color=acc,
            fontweight="bold", ha="left", va="center", zorder=5)
    ax.text(x + 3.4, SY + SH - 9.2, stat, fontsize=8.6, color=acc,
            fontweight="bold", ha="left", va="center", zorder=5)
    ax.text(x + 3.4 + 1.7 * len(stat) + 1.8, SY + SH - 9.2, statlab,
            fontsize=4.8, color="#9AA6B2", ha="left", va="center",
            zorder=5, style="italic")
    fy = SY + SH - 15.2
    for ln in lines:
        ax.text(x + 3.6, fy, "·", fontsize=6.0, color=acc, ha="left",
                va="center", zorder=5, fontweight="bold")
        ax.text(x + 5.4, fy, ln, fontsize=5.0, color=SUB, ha="left",
                va="center", zorder=5)
        fy -= 4.0

for xa in (35, 100, 160):                       # 策略支撑上层 (落在卡片底边内)
    arrow((xa, 44.4), (xa, 55.6), color="#9AA6B2", lw=0.9, ms=6, dashed=True, z=2)

for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
    fig.savefig(OUT / f"MMPV_Figure1.{ext}", **kw)
print("saved:", *(str(OUT / f"MMPV_Figure1.{e}") for e in ("png", "pdf", "svg")), sep="\n  ")
