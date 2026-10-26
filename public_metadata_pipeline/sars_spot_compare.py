#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sars_spot_compare.py — SARS-CoV-2 spot 分布：人类队列 vs 植物队列
=================================================================

输入
----
  --human   human_sars_hits.tsv        （human_sars_spot.py 产物，433,183 行）
  --plant   declared_plant_virome_rows.parquet 里的 SARS-CoV-2 行

输出
----
  fig_sars_compare.png   双面板：分布对比 + 阈值存活对比
  sars_compare_stats.tsv 两队列分位数对照表
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import duckdb
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V = os.path.expanduser("~/plant_virus_mining/10_core/12_virome/"
                       "declared_plant_virome_rows.parquet")
TARGET = "Severe acute respiratory syndrome coronavirus 2"
THRESH = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000, 100000]


def qs(v: np.ndarray) -> dict[str, float]:
    return {"n": len(v), "min": v[0], "p10": np.percentile(v, 10),
            "p25": np.percentile(v, 25), "median": np.median(v),
            "p75": np.percentile(v, 75), "p90": np.percentile(v, 90),
            "p95": np.percentile(v, 95), "p99": np.percentile(v, 99),
            "max": v[-1], "mean": v.mean()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--human", default=os.path.expanduser(
        "~/MMPV-paper/pv_spot_work/human_sars/human_sars_hits.tsv"))
    ap.add_argument("-o", "--outdir", default=os.path.expanduser(
        "~/MMPV-paper/pv_spot_work/sars_compare"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ---- 人类 ----
    hv = []
    with open(args.human, "r", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                hv.append(float(r["self_count"]))
            except (KeyError, ValueError):
                continue
    human = np.array(sorted(hv), dtype=float)
    print(f"[in] human: {len(human):,} runs")

    # ---- 植物 ----
    con = duckdb.connect()
    con.execute("SET threads TO 16")
    pv = [float(r[0]) for r in con.execute(
        f"SELECT self_count FROM read_parquet('{V}') WHERE name=? ORDER BY 1",
        [TARGET]).fetchall()]
    plant = np.array(sorted(pv), dtype=float)
    print(f"[in] plant: {len(plant):,} runs")

    # ---- 对照表 ----
    sp = os.path.join(args.outdir, "sars_compare_stats.tsv")
    with open(sp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "human", "plant"])
        h, p = qs(human), qs(plant)
        for k in ("n", "min", "p10", "p25", "median", "p75", "p90", "p95",
                  "p99", "max", "mean"):
            w.writerow([k, f"{h[k]:,.1f}", f"{p[k]:,.1f}"])
        for t in THRESH:
            hn = int((human >= t).sum()); pn = int((plant >= t).sum())
            w.writerow([f"ge{t}_n", hn, pn])
            w.writerow([f"ge{t}_pct", f"{100.0*hn/len(human):.2f}",
                        f"{100.0*pn/len(plant):.2f}"])
    print(f"[out] {sp}")

    # ---- 图 ----
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "figure.dpi": 150,
    })
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.5))

    # 左：分档占比对比（并列柱）
    BINS = [(1, 1), (2, 4), (5, 9), (10, 49), (50, 99), (100, 499),
            (500, 999), (1000, 4999), (5000, 99999), (100000, None)]
    labels, hshare, pshare = [], [], []
    for lo, hi in BINS:
        if hi is None:
            hn = int((human >= lo).sum()); pn = int((plant >= lo).sum())
            lab = f"≥{lo:,}"
        elif hi == lo:
            hn = int((human == lo).sum()); pn = int((plant == lo).sum())
            lab = str(lo)
        else:
            hn = int(((human >= lo) & (human <= hi)).sum())
            pn = int(((plant >= lo) & (plant <= hi)).sum())
            lab = f"{lo:,}–{hi:,}"
        labels.append(lab)
        hshare.append(100.0 * hn / len(human))
        pshare.append(100.0 * pn / len(plant))

    x = np.arange(len(labels))
    wd = 0.38
    ax1.bar(x - wd / 2, hshare, wd, label=f"Human (n={len(human):,})",
            color="#1E88E5", edgecolor="#37474F", linewidth=0.5)
    ax1.bar(x + wd / 2, pshare, wd, label=f"Plant (n={len(plant):,})",
            color="#66BB6A", edgecolor="#37474F", linewidth=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax1.set_ylabel("% of runs in cohort")
    ax1.set_xlabel("SARS-CoV-2 spot")
    ax1.set_title("Spot distribution: two independent cohorts", fontsize=10.5)
    ax1.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8.5, loc="upper right")

    # 右：阈值存活曲线
    hp = [100.0 * (human >= t).sum() / len(human) for t in THRESH]
    pp = [100.0 * (plant >= t).sum() / len(plant) for t in THRESH]
    ax2.plot(THRESH, hp, marker="o", color="#1E88E5", linewidth=2.2,
             markersize=6, markeredgecolor="white", markeredgewidth=0.8,
             label=f"Human (n={len(human):,})")
    ax2.plot(THRESH, pp, marker="s", color="#66BB6A", linewidth=2.2,
             markersize=6, markeredgecolor="white", markeredgewidth=0.8,
             label=f"Plant (n={len(plant):,})")
    ax2.set_xscale("log")
    ax2.set_xlabel("Spot threshold")
    ax2.set_ylabel("% of runs passing threshold")
    ax2.set_title("Survival curve: near-identical noise floor", fontsize=10.5)
    ax2.grid(linestyle=":", linewidth=0.6, alpha=0.6)
    ax2.set_axisbelow(True)
    ax2.set_ylim(-4, 112)
    ax2.legend(frameon=False, fontsize=8.5)
    for t in (100, 1000):
        ax2.axvline(t, color="#EF6C00", linestyle="--", linewidth=0.9, alpha=0.65)
        ax2.text(t, 60, f" {t}", rotation=90, fontsize=7.5, color="#EF6C00",
                 va="center")
    for t, a, b in zip(THRESH, hp, pp):
        if t in (10, 100, 1000, 10000):
            ax2.annotate(f"{a:.1f}", xy=(t, a), xytext=(0, -14),
                         textcoords="offset points", ha="center", fontsize=7,
                         color="#1E88E5")
            ax2.annotate(f"{b:.1f}", xy=(t, b), xytext=(0, 7),
                         textcoords="offset points", ha="center", fontsize=7,
                         color="#2E7D32")

    fig.suptitle("SARS-CoV-2 spot in human vs plant SRA cohorts: "
                 "a shared contamination baseline", fontsize=11.5, y=1.02)
    fig.tight_layout()
    png = os.path.join(args.outdir, "fig_sars_compare.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(png.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[out] {png}")

    # ---- 摘要 ----
    print("\n=== 对照 ===")
    print("%-10s %14s %14s" % ("metric", "human", "plant"))
    for k in ("n", "median", "p90", "p99", "max"):
        print("%-10s %14s %14s" % (k, f"{h[k]:,.0f}", f"{p[k]:,.0f}"))
    print("\n=== 阈值通过率 ===")
    print("%-10s %12s %12s" % ("thr", "human", "plant"))
    for t, a, b in zip(THRESH, hp, pp):
        print("%-10s %11.2f%% %11.2f%%" % (format(t, ","), a, b))
    return 0


if __name__ == "__main__":
    sys.exit(main())
