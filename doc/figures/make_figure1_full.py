# -*- coding: utf-8 -*-
"""MMPV Figure 1 全量版 — 七面板展示全部模块功能 (v2, 按代码 stage 定义核准)

  (a) 公共数据获取   (b) 基准驱动发现管道(13阶段+4支路拯救+基准)
  (c) 已知病毒分析(9阶段, stage 定义取自 auto_known_virus.py;
      全长基因组构建 = OmniVirusAssembler; HyPhy 已迁至 phylo)
  (d) 系统发育动力学(18 阶段, STAGE_ORDER 取自 phylo_pipeline.py)
  (e) GenBank 提交   (f) 三条设计策略   (g) 共享基础设施

用法:  python doc/figures/make_figure1_full.py
输出:  doc/figures/MMPV_Figure1_full.{png,pdf,svg}   190x228 mm
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

mpl.use("Agg")
mpl.rcParams.update({"font.family": "Arial", "svg.fonttype": "none",
                     "pdf.fonttype": 42})

W, H = 190.0, 228.0
OUT = Path(__file__).resolve().parent
INK, SUB, LINE = "#222222", "#4A4A4A", "#C9D2DA"
CA, CB, CC, CD, CE = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00"
TINT = {CA: ("#EAF2F8", "#9DC3DE"), CB: ("#FDF2E2", "#E3B36F"),
        CC: ("#E8F5EE", "#5BAF7D"), CD: ("#F7EDF4", "#D39EC2"),
        CE: ("#FBEAE3", "#E09B7A")}

fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.set_aspect("equal"); ax.axis("off")


def rbox(x, y, w, h, fc, ec, lw=0.8, r=1.2, z=3):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle=f"round,pad=0,rounding_size={r}",
                 fc=fc, ec=ec, lw=lw, zorder=z))


def arrow(p0, p1, color="#6B7A88", lw=1.2, ms=7, dashed=False, z=5):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms,
                 color=color, lw=lw,
                 linestyle=(0, (2.2, 1.6)) if dashed else "-",
                 shrinkA=0, shrinkB=0, zorder=z))


def panel(x, y, w, h, letter, title, color, dirname):
    rbox(x, y, w, h, "white", color, lw=1.0, r=1.6, z=2)
    ax.add_patch(plt.Rectangle((x, y + h - 1.6), w, 1.6, fc=color, ec="none", zorder=3))
    ax.text(x + 2.4, y + h - 4.6, f"({letter})", fontsize=8.0, color=INK,
            fontweight="bold", ha="left", va="center", zorder=4)
    tx = x + 2.4 + 1.32 * len(letter) + 3.4
    ax.text(tx, y + h - 4.6, title, fontsize=6.9, color=color,
            fontweight="bold", ha="left", va="center", zorder=4)
    tx2 = tx + 1.45 * len(title) + 3.0
    ax.text(tx2, y + h - 4.6, dirname, fontsize=4.7, color="#9AA6B2",
            style="italic", family="monospace", ha="left", va="center", zorder=4)


def chip(x, y, text, color, fs=4.9, h=5.0, bold=True, pad=2.2):
    w = pad + 0.86 * len(text)
    rbox(x, y, w, h, TINT[color][0], TINT[color][1], lw=0.7, r=1.0, z=5)
    ax.text(x + w / 2, y + h / 2, text, fontsize=fs, color=color if bold else SUB,
            ha="center", va="center", zorder=6, fontweight="bold" if bold else "normal")
    return w


def stage_pill(x, y, w, h, text, color, out=False, fs=4.5):
    rbox(x, y, w, h, "#E8F5EE" if out else TINT[color][0],
         "#5BAF7D" if out else TINT[color][1], lw=0.8, r=1.0, z=5)
    ax.text(x + w / 2, y + h / 2, text, fontsize=fs,
            color="#1F5C3D" if out else INK, ha="center", va="center",
            zorder=6, linespacing=1.05)


def elbow(x_from, x_to, y_row1_bottom, y_row2_top, color):
    """行间换行连接: row1 末列底部 → 间隙中线 → 横穿 → 进入 row2 首列顶部"""
    y_mid = (y_row1_bottom + y_row2_top) / 2
    ax.plot([x_from, x_from, x_to, x_to],
            [y_row1_bottom, y_mid, y_mid, y_row2_top + 0.8],
            color=color, lw=1.0, solid_capstyle="round", zorder=5)
    arrow((x_to, y_row2_top + 0.8), (x_to, y_row2_top), color=color, lw=1.0, ms=5)


# ═══════════ 输入条 ═══════════
ax.text(W / 2, 225.3,
        "INPUT — Public archives: NGDC GSA · NCBI SRA  |  self-generated data",
        fontsize=5.8, color="#8A97A3", ha="center", va="center", fontweight="bold")

# ═══════════ (a) 公共数据获取 y200-224 ═══════════
panel(4, 200, 182, 24, "a", "Public Data Acquisition & Curation", CA,
      "public_metadata_pipeline")
pills_row = ["GSA + SRA\ndual search", "AI metadata\ncuration (DeepSeek)",
             "Batch download\naria2c · prefetch", "SRA → FASTQ\nconversion",
             "6-panel\nlandscape plots"]
pw5, pg5 = 32.0, 3.0
for i, t in enumerate(pills_row):
    x = 6.5 + i * (pw5 + pg5)
    stage_pill(x, 206.5, pw5, 6.8, t, CA, fs=4.8)
    if i < 4:
        arrow((x + pw5 + 0.35, 206.5 + 3.4), (x + pw5 + pg5 - 0.35, 206.5 + 3.4),
              color=TINT[CA][1], lw=1.0, ms=5)
chip(6.5, 200.8, "Host genome DB: Kraken2 · Bowtie2 · HISAT2 · Minimap2 indexes", CA)
chip(118.0, 200.8, "Desktop GUI — AI completion · visualize · export", CA)

# ═══════════ (b) 发现管道 y150-196 ═══════════
panel(4, 150, 182, 46, "b", "Benchmark-Driven Discovery Pipeline", CB,
      "virome_discovery_pipeline")
ax.text(6.5, 185.2, "13 stages · every default parameter benchmark-validated",
        fontsize=5.3, color=INK, style="italic", ha="left", va="center", zorder=5)
row1 = ["Quality\ncontrol", "Host\nfiltering", "Depth\nnormalization", "De novo\nassembly",
        "Multi-tool\nidentification", "UniProt-strict\nfiltering", "COBRA\nextension"]
row2 = ["Flye\nco-assembly", "Ref-guided\nclustering", "Multi-tool\ntaxonomy",
        "3-tier host\nprediction", "CheckV\ncompleteness", "4-branch\ncascade rescue",
        "Integrated\nreport"]
PW, PG, PH = 14.5, 2.6, 6.8
R1Y, R2Y = 173.6, 163.4
for j, t in enumerate(row1):
    x = 6.5 + j * (PW + PG)
    stage_pill(x, R1Y, PW, PH, t, CB)
    if j < 6:
        arrow((x + PW + 0.35, R1Y + PH / 2), (x + PW + PG - 0.35, R1Y + PH / 2),
              color=TINT[CB][1], lw=1.0, ms=5)
for j, t in enumerate(row2):
    x = 6.5 + j * (PW + PG)
    out = (j == 6)
    stage_pill(x, R2Y, PW, PH, t, CB, out=out)
    if j < 6:
        arrow((x + PW + 0.35, R2Y + PH / 2), (x + PW + PG - 0.35, R2Y + PH / 2),
              color="#5BAF7D" if j == 5 else TINT[CB][1], lw=1.0, ms=5)
elbow(6.5 + 6 * (PW + PG) + PW / 2, 6.5 + PW / 2, R1Y, R2Y, TINT[CB][1])
bx = 6.5
for t in ["6-D benchmark · 20+ tools", "1,536 identification combinations",
          "6 assemblers × 54 depth × background"]:
    bx += chip(bx, 152.6, t, CB) + 2.6
# 4-branch rescue inset
rbox(126.0, 152.6, 57.5, 27.7, "#FFFBF5", "#E3B36F", lw=0.9, r=1.2, z=4)
ax.text(127.6, 177.4, "4-branch cascade rescue", fontsize=5.4, color="#9C6500",
        fontweight="bold", ha="left", va="center", zorder=6)
branches = [("A", "CheckV ≥ 90% complete → pass"),
            ("B", "Virseqimprover extension"),
            ("C", "BLASTN → PlantVirusDB rebuild:\nbwa-mem2 + viral_consensus / ragtag"),
            ("D", "Genus-length rescue (fallback)")]
by = 173.6
for tag, txt in branches:
    rbox(127.6, by - 1.7, 3.4, 3.4, CB, CB, lw=0, r=0.8, z=6)
    ax.text(129.3, by, tag, fontsize=4.6, color="white", ha="center",
            va="center", zorder=7, fontweight="bold")
    nl = txt.count("\n") + 1
    ax.text(132.4, by + (nl - 1) * 1.7, txt, fontsize=4.7, color=SUB,
            ha="left", va="top", zorder=6, linespacing=1.25)
    by -= 4.6 if nl == 1 else 6.4

# ═══════════ (c) 已知病毒分析 y112-146 (9 stages, auto_known_virus.py) ═══════════
panel(4, 112, 182, 34, "c", "Rapid Known-Virus Analysis", CC,
      "virome_analysis_pipeline")
crow1 = ["Known-virus\ndetection & quant.", "High-confidence\nfiltering",
         "Variant calling\nSnpEff · SnpGenie", "Full-length genome\nconstruction",
         "N-fill extraction\n(extract_full_fasta)"]
crow2 = ["Post-analysis\nVCF · MAF", "Similarity panorama\n(SDT)",
         "DVG & recomb.\n(ViReMa)", "HTML report\n+ AI interpretation"]
CW, CG, CH_ = 21.3, 2.6, 6.8
C1Y, C2Y = 128.6, 119.8
for j, t in enumerate(crow1):
    x = 6.5 + j * (CW + CG)
    stage_pill(x, C1Y, CW, CH_, t, CC, fs=4.5)
    if j < 4:
        arrow((x + CW + 0.35, C1Y + CH_ / 2), (x + CW + CG - 0.35, C1Y + CH_ / 2),
              color=TINT[CC][1], lw=1.0, ms=5)
for j, t in enumerate(crow2):
    x = 6.5 + j * (CW + CG)
    out = (j == 3)
    stage_pill(x, C2Y, CW, CH_, t, CC, out=out, fs=4.5)
    if j < 3:
        arrow((x + CW + 0.35, C2Y + CH_ / 2), (x + CW + CG - 0.35, C2Y + CH_ / 2),
              color="#5BAF7D" if j == 2 else TINT[CC][1], lw=1.0, ms=5)
elbow(6.5 + 4 * (CW + CG) + CW / 2, 6.5 + CW / 2, C1Y, C2Y, TINT[CC][1])
# OmniVirusAssembler inset (stage 4 `full`)
rbox(126.0, 119.8, 57.5, 15.6, "#F2FBF6", "#5BAF7D", lw=0.9, r=1.2, z=4)
ax.text(127.6, 133.2, "OmniVirusAssembler (stage `full`)", fontsize=5.2,
        color="#1F5C3D", fontweight="bold", ha="left", va="center", zorder=6)
omni = ["12-step refinement · SHIVER-like fusion skeleton",
        "Reads-level iterative polishing — single-base precision",
        "Dual-engine gap closing: gmcloser + abyss-sealer",
        "Circularization detection → [Circular=True]"]
oy = 129.6
for ln in omni:
    ax.text(127.8, oy, "·", fontsize=6.0, color=CC, ha="left", va="center",
            zorder=6, fontweight="bold")
    ax.text(129.6, oy, ln, fontsize=4.6, color=SUB, ha="left", va="center", zorder=6)
    oy -= 2.95
cx = 6.5
for t in ["Pan-virome: prevalence · landscape · oncoprint",
          "Co-abundance network", "snpeff annotation"]:
    cx += chip(cx, 113.2, t, CC) + 2.4

# ═══════════ (d) 系统发育动力学 y84-108 (18 stages, phylo_pipeline.py) ═══════════
panel(4, 84, 182, 24, "d", "Phylodynamics & Evolution", CD,
      "virome_phylo_pipeline")
pills_d = ["Data prep\ngovernance·online", "Phylogeny\nalign·tree·SplitsTree",
           "Recombination\nRDP5", "Selection\nFEL·MEME·BUSTED·PRIME",
           "PopGen\nhost · DNASP", "Time\nBEAST·clock·dating",
           "Geography\nphylogeo·MCC", "Transmission\npaths·maps·GIF"]
for i, t in enumerate(pills_d):
    x = 6.5 + i * 22.0
    stage_pill(x, 91.2, 20.0, 6.8, t, CD, fs=4.5, out=(i == 7))
    if i < 7:
        arrow((x + 20.0 + 0.3, 91.2 + 3.4), (x + 22.0 - 0.3, 91.2 + 3.4),
              color=TINT[CD][1], lw=1.0, ms=5)
dx = 6.5
for t in ["Metadata governance: time/location · 3-channel split",
          "Genetic subsampling · GLM predictors"]:
    dx += chip(dx, 85.6, t, CD) + 2.4

# ═══════════ (e) GenBank 提交 y62-78 ═══════════
panel(4, 62, 182, 16, "e", "Semi-automatic GenBank Submission", CE,
      "virome_submission_pipeline")
pills_e = ["Topology\ncircular / linear", "Metadata\ntemplates (MIUVIG)",
           "Hypothetical\n5-tool QA", "Sequin\n.tbl build",
           "suvtk / Sequin\ndual → .sqn", "HTML\nreport"]
for i, t in enumerate(pills_e):
    x = 6.5 + i * 21.0
    stage_pill(x, 65.0, 19.0, 6.8, t, CE, fs=4.5, out=(i == 5))
    if i < 5:
        arrow((x + 19.0 + 0.3, 65.0 + 3.4), (x + 21.0 - 0.3, 65.0 + 3.4),
              color=TINT[CE][1], lw=1.0, ms=5)
chip(140.0, 65.0, "Desktop GUI (PySide6)", CE)

# ═══════════ 主流程箭头 ═══════════
for y0, y1 in ((200, 196.4), (150, 146.4), (112, 108.4), (84, 80.4)):
    arrow((31, y0), (31, y1), lw=1.3, ms=7)

# ═══════════ (f) 设计策略 y23-59 ═══════════
ax.add_patch(FancyBboxPatch((4, 23), 182, 36,
             boxstyle="round,pad=0,rounding_size=2.0",
             fc="#F6F7F8", ec=LINE, lw=0.7, zorder=1))
ax.text(6.6, 57.3, "(f)  CROSS-CUTTING DESIGN STRATEGIES", fontsize=5.8,
        color="#8A97A3", ha="left", va="top", fontweight="bold", zorder=2)
strategies = [
    (CA, "Data Reduction", "99.8%", "host-filtering efficiency",
     ["Reference-guided pre-clustering", "Host-based sequence partitioning"]),
    (CB, "Quality Enhancement", "4-branch", "cascade genome rescue",
     ["COBRA resolves de Bruijn breakpoints", "Flye co-assembly merges cross-sample",
      "BLASTN → PlantVirusDB · viral_consensus / ragtag",
      "Genus-length rescue as final fallback"]),
    (CC, "Reliability Enhancement", "F1 = 0.942", "identification (UniProt-strict)",
     ["UniProt-strict ID ≈ 5-tool ensemble", "Family accuracy 94.1% (+72.3 pp)",
      "3-tier host + plant-specific filter"]),
]
SX, SW, SG, SY, SH = 6.5, 57.7, 3.0, 26.0, 28.0
for i, (acc, name, stat, statlab, lines) in enumerate(strategies):
    x = SX + i * (SW + SG)
    rbox(x, SY, SW, SH, "white", "#D8DCE0", lw=0.7)
    ax.add_patch(plt.Rectangle((x, SY), 1.6, SH, fc=acc, ec="none", zorder=4))
    ax.text(x + 3.4, SY + SH - 3.4, name, fontsize=6.1, color=acc,
            fontweight="bold", ha="left", va="center", zorder=5)
    ax.text(x + 3.4, SY + SH - 9.4, stat, fontsize=8.4, color=acc,
            fontweight="bold", ha="left", va="center", zorder=5)
    ax.text(x + 3.4 + 1.7 * len(stat) + 1.8, SY + SH - 9.4, statlab,
            fontsize=4.8, color="#9AA6B2", ha="left", va="center",
            zorder=5, style="italic")
    fy = SY + SH - 15.6
    for ln in lines:
        ax.text(x + 3.6, fy, "·", fontsize=6.0, color=acc, ha="left",
                va="center", zorder=5, fontweight="bold")
        ax.text(x + 5.4, fy, ln, fontsize=4.9, color=SUB, ha="left",
                va="center", zorder=5)
        fy -= 3.9
for xa in (35, 100, 160):
    arrow((xa, 59.4), (xa, 61.6), color="#9AA6B2", lw=0.9, ms=6, dashed=True, z=2)

# ═══════════ (g) 共享基础设施 y4-19 ═══════════
ax.add_patch(FancyBboxPatch((4, 4), 182, 15,
             boxstyle="round,pad=0,rounding_size=2.0",
             fc="#F5F4EF", ec=LINE, lw=0.7, zorder=1))
ax.text(6.6, 17.4, "(g)  SHARED INFRASTRUCTURE", fontsize=5.8, color="#8A97A3",
        ha="left", va="top", fontweight="bold", zorder=2)
infra = [
    ("Tool Library", "30+ wrapped: Kraken2 · Bowtie2 · snpEff · VirHunter\nViReMa · BEAST · RDP5 · HyPhy · CheckV"),
    ("Reference Databases", "ICTV VMR taxonomy · PlantVirusDB\nhost genomes · UniProt / VP"),
    ("Metadata Governance", "Time/location parsing & normalization\nthree-channel split"),
    ("Reproducible Environment", "pixi-pinned environments\nstage checkpoints · one-command reruns"),
]
IX, IW, IGAP, IY, IH = 6.5, 43.0, 2.3, 5.4, 9.4
for i, (t, body) in enumerate(infra):
    x = IX + i * (IW + IGAP)
    rbox(x, IY, IW, IH, "white", "#D5D2C8", lw=0.7)
    ax.text(x + 2.2, IY + IH - 2.6, t, fontsize=5.8, color="#5A5648",
            fontweight="bold", ha="left", va="center", zorder=4)
    ax.text(x + 2.2, IY + 2.7, body, fontsize=4.6, color=SUB, ha="left",
            va="center", zorder=4, linespacing=1.3)

for ext, kw in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
    fig.savefig(OUT / f"MMPV_Figure1_full.{ext}", **kw)
print("saved:", *(str(OUT / f"MMPV_Figure1_full.{e}") for e in ("png", "pdf", "svg")), sep="\n  ")
