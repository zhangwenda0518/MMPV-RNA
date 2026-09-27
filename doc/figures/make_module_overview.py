# -*- coding: utf-8 -*-
"""MMPV 模块总览图 (论文 Figure 1 / 图形摘要) — 分层架构 + 层内数据流

用法:  python doc/figures/make_module_overview.py
输出:  doc/figures/MMPV_module_overview.{png,pdf,svg}
配色:  Okabe-Ito (与管线内已发表图风格一致); 画幅 190x128 mm
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

mpl.use("Agg")
mpl.rcParams.update({
    "font.family": "Arial",
    "svg.fonttype": "none",          # SVG 保留可编辑文本
    "pdf.fonttype": 42,              # PDF 内嵌 TrueType (可编辑)
})

W, H = 190.0, 128.0                 # mm
OUT = Path(__file__).resolve().parent

# Okabe-Ito
C1, C2, C3, C4, C5 = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"
INK, SUB, LINE = "#222222", "#555555", "#C9D2DA"

fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.set_aspect("equal"); ax.axis("off")


def band(x, y, w, h, fill, tag):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0,rounding_size=2.0",
                 fc=fill, ec=LINE, lw=0.7, zorder=1))
    ax.text(x + 2.6, y + h - 2.6, tag, fontsize=5.8, color="#8A97A3",
            ha="left", va="top", fontweight="bold", zorder=2)


def card(x, y, w, h, title, tcolor, lines, tfs=7.0, bfs=6.1, dir_name=None,
         fc="white", ec="#B7C2CC", center=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0,rounding_size=1.4",
                 fc=fc, ec=ec, lw=0.8, zorder=3))
    cy = y + h - 3.6
    ax.text(x + w / 2 if center else x + 2.4, cy, title,
            fontsize=tfs, color=tcolor, fontweight="bold",
            ha="center" if center else "left", va="center", zorder=4)
    cy -= 4.0
    if dir_name:
        ax.text(x + w / 2 if center else x + 2.4, cy, dir_name,
                fontsize=5.0, color="#9AA6B2", style="italic", family="monospace",
                ha="center" if center else "left", va="center", zorder=4)
        cy -= 3.6
    for ln in lines:
        ax.text(x + w / 2 if center else x + 2.8, cy, ln,
                fontsize=bfs, color=INK if center else SUB,
                ha="center" if center else "left", va="center", zorder=4)
        cy -= 4.35


# ══════════════════ 交互层 ══════════════════
band(4, 103, 182, 21, "#F2F6FA", "INTERACTIVE LAYER")
card(30, 105.5, 62, 13.5, "Metadata GUI", "#3E5C76",
     ["Search · edit · visualize unified metadata"],
     tfs=7.2, bfs=6.0, dir_name="metadata_gui", center=True)
card(98, 105.5, 62, 13.5, "Submission GUI", "#3E5C76",
     ["Guided metadata entry · GenBank packaging"],
     tfs=7.2, bfs=6.0, dir_name="submission_gui", center=True)

for x0, x1, lab in ((61, 56, "metadata in"), (129, 169, "assists module 5")):
    ax.add_patch(FancyArrowPatch((x0, 105.5), (x1, 96.5), arrowstyle="-|>",
                 mutation_scale=7, linestyle=(0, (2.2, 1.6)), color="#9AA6B2",
                 lw=0.8, shrinkA=0, shrinkB=0, zorder=2))
    ax.text((x0 + x1) / 2 + 2.5, 101.2, lab, fontsize=4.8, color="#9AA6B2",
            style="italic", ha="left", va="center", zorder=2)

# ══════════════════ 管线层 ══════════════════
band(4, 22, 182, 75, "#FBFCFE", "PIPELINE LAYER · CORE WORKFLOW")

modules = [
    (C1, "Public Data\nAcquisition", "public_metadata_pipeline",
     ["SRA / GSA dual search", "AI-assisted query terms",
      "Metadata retrieval", "Batch download & convert",
      "Host genome DB build"],
     "→ FASTQs + metadata"),
    (C2, "Novel Virus\nDiscovery", "virome_discovery_pipeline",
     ["QC & host depletion", "De novo assembly",
      "10-tool identification", "Taxonomy voting (8-way)",
      "Host prediction", "Rescue & report"],
     "→ contigs + taxonomy"),
    (C3, "Known-Virus\nAnalysis", "virome_analysis_pipeline",
     ["Quantification (TPM/RPM)", "Variant calling",
      "snpEff annotation", "Co-infection analysis",
      "Pan-virome landscape"],
     "→ variants + quant tables"),
    (C4, "Phylodynamics", "virome_phylo_pipeline",
     ["Recombination (RDP5)", "Selection (HyPhy)",
      "Population genetics", "BEAST dating · MCC",
      "Transmission maps"],
     "→ trees + dynamics figs"),
    (C5, "Data\nSubmission", "virome_submission_pipeline",
     ["Topology triage", "MIUVIG metadata",
      "Sequin / suvtk", "tbl2asn .sqn build",
      "Hypothetical protein QA"],
     "→ submission-ready .sqn"),
]

MX, MW, GAP, MY, MH = 6.5, 33.5, 3.0, 25.5, 64.5
for i, (col, title, dirname, funcs, out_chip) in enumerate(modules):
    x = MX + i * (MW + GAP)
    ax.add_patch(FancyBboxPatch((x, MY), MW, MH,
                 boxstyle="round,pad=0,rounding_size=1.4",
                 fc="white", ec=col, lw=0.9, zorder=3))
    hh = 11.5
    ax.add_patch(FancyBboxPatch((x, MY + MH - hh), MW, hh,
                 boxstyle="round,pad=0,rounding_size=1.4",
                 fc=col, ec=col, lw=0.9, zorder=4))
    ax.add_patch(plt.Rectangle((x, MY + MH - hh), MW, 2.0, fc=col, ec="none", zorder=4))
    ax.text(x + MW / 2, MY + MH - 3.0, f"MODULE {i + 1}", fontsize=5.2,
            color="white", alpha=0.85, ha="center", va="center",
            fontweight="bold", zorder=5)
    ax.text(x + MW / 2, MY + MH - 7.6, title, fontsize=7.3, color="white",
            ha="center", va="center", fontweight="bold", zorder=5, linespacing=1.05)
    ax.text(x + MW / 2, MY + MH - hh - 3.4, dirname, fontsize=4.9, color="#9AA6B2",
            style="italic", family="monospace", ha="center", va="center", zorder=5)
    fy = MY + MH - hh - 9.0
    for fn in funcs:
        ax.text(x + 2.0, fy, "·", fontsize=6.5, color=col, ha="left",
                va="center", zorder=5, fontweight="bold")
        ax.text(x + 4.0, fy, fn, fontsize=6.1, color=SUB, ha="left",
                va="center", zorder=5)
        fy -= 5.55
    # 卡片底部: 产物输出芯片 (填充纵向空间 + 传达模块间数据流)
    ax.plot([x + 2.0, x + MW - 2.0], [MY + 8.2, MY + 8.2], color=LINE,
            lw=0.6, zorder=5)
    ax.text(x + 2.0, MY + 5.4, out_chip, fontsize=5.0, color=col,
            style="italic", ha="left", va="center", zorder=5, fontweight="bold")
    if i < 4:
        ax.add_patch(FancyArrowPatch((x + MW + 0.4, MY + MH - hh / 2),
                     (x + MW + GAP - 0.4, MY + MH - hh / 2), arrowstyle="-|>",
                     mutation_scale=8, color="#6B7A88", lw=1.3,
                     shrinkA=0, shrinkB=0, zorder=5))

# ══════════════════ 设施层 ══════════════════
band(4, 3, 182, 15, "#F5F4EF", "SHARED INFRASTRUCTURE")

infra = [
    ("Tool Library", "30+ wrapped tools: Kraken2 · Bowtie2 · snpEff\nVirHunter · ViReMa · BEAST · RDP5"),
    ("Reference Databases", "ICTV VMR taxonomy · plant-virus\nreference DB · host genomes"),
    ("Metadata Governance", "Time/location parsing, normalization\n& three-channel split"),
    ("Reproducible Environment", "pixi-pinned environments · stage\ncheckpoints · one-command reruns"),
]
IX, IW, IGAP, IY, IH = 6.5, 43.0, 2.3, 4.3, 9.2
for i, (t, body) in enumerate(infra):
    x = IX + i * (IW + IGAP)
    ax.add_patch(FancyBboxPatch((x, IY), IW, IH,
                 boxstyle="round,pad=0,rounding_size=1.2",
                 fc="white", ec="#D5D2C8", lw=0.7, zorder=3))
    ax.text(x + 2.2, IY + IH - 2.6, t, fontsize=6.0, color="#5A5648",
            fontweight="bold", ha="left", va="center", zorder=4)
    ax.text(x + 2.2, IY + 2.9, body, fontsize=5.0, color=SUB, ha="left",
            va="center", zorder=4, linespacing=1.35)

for xa in (48, 95, 142):
    ax.add_patch(FancyArrowPatch((xa, 21.6), (xa, 18.4), arrowstyle="-|>",
                 mutation_scale=6, linestyle=(0, (2.2, 1.6)), color="#B9B4A6",
                 lw=0.8, shrinkA=0, shrinkB=0, zorder=2))

for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
    fig.savefig(OUT / f"MMPV_module_overview.{ext}", **kw)
print("saved:", *(str(OUT / f"MMPV_module_overview.{e}") for e in ("png", "pdf", "svg")), sep="\n  ")
