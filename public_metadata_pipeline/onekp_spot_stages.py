#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
onekp_spot_stages.py — onekp 四阶段 spot 分布：箱线图 + 折线图
==============================================================

四阶段
------
  raw          01_detection/summary/all_viruses.raw.tsv        （原始比对）
  best         01_detection/summary/all_viruses.best.summary.tsv（白名单）
  high_conf    02_filtering/high_conf.summary.tsv             （质量过滤）
  extracted    06_extraction/extraction_status.tsv 中 status=extracted
               （完整组装，spot 回查 high_conf）

口径
----
  spot = `Asm_EM_Reads`（raw 表里叫 `EM_Reads`），即 EM 估计的 reads 数。
  四个阶段都是同一口径，可直接纵向比较。

输出
----
  fig_spot_stages_box.png    箱线图（log10 y 轴，四阶段并列）
  fig_spot_stages_line.png   折线图（各分位数随阶段变化 + 阈值存活曲线）
  spot_stages_stats.tsv      各阶段分位数明细
  spot_min_per_virus.tsv     每个病毒在四阶段的最小值

用法
----
    python onekp_spot_stages.py -o ./fig_spot
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.expanduser("~/MMPV-paper/onekp-virome/onekp-analysis")

STAGES = [
    ("raw", "Raw hits", os.path.join(ROOT, "01_detection/summary/all_viruses.raw.tsv"),
     ("Asm_EM_Reads", "EM_Reads"), "Species"),
    ("best", "Best (whitelist)", os.path.join(ROOT, "01_detection/summary/all_viruses.best.summary.tsv"),
     ("Asm_EM_Reads", "EM_Reads"), "Species_NCBI"),
    ("high_conf", "High-confidence", os.path.join(ROOT, "02_filtering/high_conf.summary.tsv"),
     ("Asm_EM_Reads", "EM_Reads"), "Species_NCBI"),
]

COLORS = ["#B0BEC5", "#64B5F6", "#42A5F5", "#1E88E5"]
THRESH = [10, 100, 1000, 10000]


def read_tsv(path: str) -> list[dict]:
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def get_val(row: dict, names: tuple[str, ...]) -> float | None:
    for n in names:
        v = row.get(n)
        if v not in (None, ""):
            try:
                return float(v)
            except ValueError:
                continue
    return None


def load() -> tuple[dict[str, list[float]], dict[str, dict[str, list[float]]], dict[str, str]]:
    """返回 (阶段 -> spot 列表, 阶段 -> 病毒 -> spot 列表, 病毒 -> 显示名)"""
    by_stage: dict[str, list[float]] = {}
    by_virus: dict[str, dict[str, list[float]]] = {}
    hc_rows = None

    for key, _label, path, cols, spcol in STAGES:
        rows = read_tsv(path)
        if key == "high_conf":
            hc_rows = rows
        vals: list[float] = []
        vv: dict[str, list[float]] = defaultdict(list)
        for r in rows:
            v = get_val(r, cols)
            if v is None:
                continue
            vals.append(v)
            sp = (r.get(spcol) or "").strip()
            if sp:
                vv[sp].append(v)
        by_stage[key] = vals
        by_virus[key] = dict(vv)
        print(f"[in] {key:10s} {len(vals):6d} 条, {len(vv):4d} 种病毒")

    # extracted：从 extraction_status 回查 high_conf
    idx: dict[tuple[str, str], dict] = {}
    for r in hc_rows or []:
        idx[(r["Sample"], r["Rep_Accession"])] = r
        for s in (r.get("Segment_Accessions") or "").split(","):
            s = s.strip().split(":")[-1]
            if s:
                idx.setdefault((r["Sample"], s), r)

    vals = []
    vv = defaultdict(list)
    n_miss = 0
    for e in read_tsv(os.path.join(ROOT, "06_extraction/extraction_status.tsv")):
        if e.get("status") != "extracted":
            continue
        acc, _, rest = e["sample"].partition("_")
        r = idx.get((acc, rest))
        if r is None:
            n_miss += 1
            continue
        v = get_val(r, ("Asm_EM_Reads", "EM_Reads"))
        if v is None:
            continue
        vals.append(v)
        sp = (r.get("Species_NCBI") or "").strip()
        if sp:
            vv[sp].append(v)
    by_stage["extracted"] = vals
    by_virus["extracted"] = dict(vv)
    print(f"[in] {'extracted':10s} {len(vals):6d} 条, {len(vv):4d} 种病毒"
          f"  ({n_miss} 条在 high_conf 中找不到，未计入)")
    return by_stage, by_virus, {}


def stats(vals: list[float]) -> dict[str, float]:
    v = np.array(sorted(vals), dtype=float)
    return {
        "n": len(v), "min": v[0], "p10": np.percentile(v, 10),
        "p25": np.percentile(v, 25), "median": np.median(v),
        "p75": np.percentile(v, 75), "p90": np.percentile(v, 90),
        "max": v[-1], "mean": v.mean(),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outdir", default="./fig_spot")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    by_stage, by_virus, _ = load()
    order = [k for k, *_ in STAGES] + ["extracted"]
    labels = ["Raw hits", "Best\n(whitelist)", "High-\nconfidence", "Complete\nassembly"]

    # ---------------- 统计表 ----------------
    stats_path = os.path.join(args.outdir, "spot_stages_stats.tsv")
    with open(stats_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["stage", "n", "min", "p10", "p25", "median", "p75", "p90",
                    "max", "mean"] + [f"ge{t}" for t in THRESH])
        for k in order:
            s = stats(by_stage[k])
            row = [k, s["n"], f"{s['min']:.1f}", f"{s['p10']:.1f}", f"{s['p25']:.1f}",
                   f"{s['median']:.1f}", f"{s['p75']:.1f}", f"{s['p90']:.1f}",
                   f"{s['max']:.1f}", f"{s['mean']:.1f}"]
            row += [sum(1 for x in by_stage[k] if x >= t) for t in THRESH]
            w.writerow(row)
    print(f"[out] {stats_path}")

    # ---------------- 每病毒最小值 ----------------
    mv_path = os.path.join(args.outdir, "spot_min_per_virus.tsv")
    viruses = sorted(by_virus["extracted"], key=lambda x: -len(by_virus["extracted"][x]))
    with open(mv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["virus", "n_extracted", "min_raw", "min_best", "min_high_conf",
                    "min_extracted", "median_extracted"])
        for v in viruses:
            def mn(stage):
                a = by_virus.get(stage, {}).get(v)
                return f"{min(a):.0f}" if a else ""
            med = by_virus["extracted"][v]
            w.writerow([v, len(med), mn("raw"), mn("best"), mn("high_conf"),
                        mn("extracted"), f"{np.median(med):.0f}"])
    print(f"[out] {mv_path}")

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "figure.dpi": 150,
    })

    # ============ 图 1：箱线图 ============
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    data = [by_stage[k] for k in order]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=True,
                    flierprops=dict(marker="o", markersize=2.2, alpha=0.35,
                                    markerfacecolor="#78909C", markeredgecolor="none"),
                    medianprops=dict(color="#C62828", linewidth=1.6),
                    whiskerprops=dict(color="#546E7A", linewidth=0.9),
                    capprops=dict(color="#546E7A", linewidth=0.9))
    for patch, c in zip(bp["boxes"], COLORS):
        patch.set_facecolor(c)
        patch.set_edgecolor("#37474F")
        patch.set_linewidth(0.8)

    ax.set_yscale("log")
    ax.set_xticklabels(labels)
    ax.set_ylabel("Spot (EM-estimated reads, log scale)")
    ax.set_title("onekp: spot distribution across four stages", fontsize=11, pad=10)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6, which="major")
    ax.set_axisbelow(True)

    # 阈值参考线
    for t, ls in ((10, "--"), (100, "-"), (1000, ":")):
        ax.axhline(t, color="#EF6C00", linestyle=ls, linewidth=0.9, alpha=0.7, zorder=1)
        ax.text(4.62, t, f" {t}", va="center", ha="left", fontsize=8, color="#EF6C00")

    # 各阶段最小值：单独一行标在图下方，避免与箱体/轴标签重叠
    ax.annotate("min", xy=(-0.62, -0.155), xycoords=("data", "axes fraction"),
                ha="right", va="center", fontsize=8, color="#B71C1C",
                fontweight="bold")
    for i, k in enumerate(order, start=1):
        m = min(by_stage[k])
        ax.annotate(f"{m:,.0f}", xy=(i, -0.155), xycoords=("data", "axes fraction"),
                    ha="center", va="center", fontsize=8, color="#B71C1C",
                    fontweight="bold")
    ax.set_xlim(0.5, 5.15)
    fig.tight_layout()
    p1 = os.path.join(args.outdir, "fig_spot_stages_box.png")
    fig.savefig(p1, dpi=300, bbox_inches="tight")
    fig.savefig(p1.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[out] {p1}")

    # ============ 图 2：折线图（双面板）============
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6))

    # 左：分位数随阶段变化
    x = np.arange(len(order))
    qspecs = [("min", "#B71C1C", "o", 1.9), ("p10", "#EF6C00", "s", 1.4),
              ("p25", "#FDD835", "^", 1.4), ("median", "#2E7D32", "D", 1.9),
              ("p75", "#1565C0", "v", 1.4), ("p90", "#6A1B9A", "P", 1.4),
              ("max", "#455A64", "*", 1.9)]
    for name, color, marker, lw in qspecs:
        ys = [stats(by_stage[k])[name] for k in order]
        ax1.plot(x, ys, marker=marker, color=color, linewidth=lw,
                 markersize=5.5, label=name, markeredgecolor="white",
                 markeredgewidth=0.6, zorder=3)
    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=9)
    ax1.set_ylabel("Spot (log scale)")
    ax1.set_title("Quantiles across stages", fontsize=10.5)
    ax1.grid(linestyle=":", linewidth=0.6, alpha=0.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8, ncol=2, loc="lower right")
    for i, k in enumerate(order):
        ax1.annotate(f"{stats(by_stage[k])['median']:,.0f}",
                     xy=(i, stats(by_stage[k])["median"]), xytext=(0, 9),
                     textcoords="offset points", ha="center", fontsize=7.5,
                     color="#2E7D32", fontweight="bold")

    # 右：阈值存活曲线
    for t, color in zip(THRESH, ["#90CAF9", "#42A5F5", "#1E88E5", "#0D47A1"]):
        ys = []
        for k in order:
            v = by_stage[k]
            ys.append(100.0 * sum(1 for z in v if z >= t) / len(v) if v else 0.0)
        ax2.plot(x, ys, marker="o", color=color, linewidth=1.9, markersize=6,
                 label=f"spot ≥ {t:,}", markeredgecolor="white", markeredgewidth=0.7)
        for xi, yi in zip(x, ys):
            if yi > 0:
                ax2.annotate(f"{yi:.0f}%", xy=(xi, yi), xytext=(0, 8),
                             textcoords="offset points", ha="center",
                             fontsize=7, color=color)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("Fraction of records passing threshold (%)")
    ax2.set_ylim(-3, 108)
    ax2.set_title("Survival rate by spot threshold", fontsize=10.5)
    ax2.grid(linestyle=":", linewidth=0.6, alpha=0.6)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, fontsize=8, loc="center left")

    fig.suptitle("onekp virome pipeline: spot evidence through four filtering stages",
                 fontsize=11.5, y=1.01)
    fig.tight_layout()
    p2 = os.path.join(args.outdir, "fig_spot_stages_line.png")
    fig.savefig(p2, dpi=300, bbox_inches="tight")
    fig.savefig(p2.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[out] {p2}")

    # ---------------- 控制台摘要 ----------------
    print("\n=== 四阶段 spot 分位 ===")
    print("%-12s %6s %10s %10s %10s %10s %10s" %
          ("stage", "n", "min", "median", "p90", "max", ">=100"))
    for k in order:
        s = stats(by_stage[k])
        ge = sum(1 for x in by_stage[k] if x >= 100)
        print("%-12s %6d %10.0f %10.0f %10.0f %10.0f %6d (%.0f%%)" %
              (k, s["n"], s["min"], s["median"], s["p90"], s["max"],
               ge, 100.0 * ge / s["n"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
