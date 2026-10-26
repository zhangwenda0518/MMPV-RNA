#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sars_spot_dist.py — SARS-CoV-2 在植物队列中的 spot 分布（污染案例）
===================================================================

数据来源
--------
  plant_virus_mining/10_core/12_virome/declared_plant_virome_rows.parquet
  种级行：name = 'Severe acute respiratory syndrome coronavirus 2'
  队列：声明为植物的 SRA Run（1,910,691 个有病毒 STAT 行的 run）

口径
----
  spot = self_count（该种节点直接命中）。
  该 taxid 在所有 run 里 self_count == total_count（无下游株系节点）。

输出
----
  sars_spot_summary.tsv    分位数 / 分档 / 阈值通过率
  sars_spot_per_year.tsv   逐年 run 数 + spot 分位数
  fig_sars_spot_dist.png   三联图：直方图 + 阈值存活 + 逐年
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

B = os.path.expanduser("~/plant_virus_mining")
V = B + "/10_core/12_virome/declared_plant_virome_rows.parquet"
M = B + "/00_data/03_queue_master/plant_runs_master.parquet"
TARGET = "Severe acute respiratory syndrome coronavirus 2"

THRESH = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000]
BINS = [(1, 1), (2, 4), (5, 9), (10, 49), (50, 99), (100, 499),
        (500, 999), (1000, 4999), (5000, None)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outdir", default="./sars_spot")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    con = duckdb.connect()
    con.execute("SET threads TO 16")

    n_rows, n_runs = con.execute(
        f"SELECT COUNT(*), COUNT(DISTINCT acc) FROM read_parquet('{V}') WHERE name=?",
        [TARGET]).fetchone()
    print(f"[in] SARS-CoV-2: {n_rows:,} 行 / {n_runs:,} run")

    vals = np.array([r[0] for r in con.execute(
        f"SELECT self_count FROM read_parquet('{V}') WHERE name=? ORDER BY 1",
        [TARGET]).fetchall()], dtype=float)
    print(f"[in] self_count: min={vals[0]:.0f} median={np.median(vals):.0f} max={vals[-1]:,.0f}")

    # ---- 汇总表 ----
    sp = os.path.join(args.outdir, "sars_spot_summary.tsv")
    with open(sp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in (("n_rows", n_rows), ("n_runs", n_runs),
                     ("min", vals[0]), ("p10", np.percentile(vals, 10)),
                     ("p25", np.percentile(vals, 25)), ("median", np.median(vals)),
                     ("p75", np.percentile(vals, 75)), ("p90", np.percentile(vals, 90)),
                     ("p95", np.percentile(vals, 95)), ("p99", np.percentile(vals, 99)),
                     ("max", vals[-1]), ("mean", vals.mean()),
                     ("sum", vals.sum())):
            w.writerow([k, f"{v:,.1f}" if isinstance(v, float) else v])
        for t in THRESH:
            n = int((vals >= t).sum())
            w.writerow([f"ge{t}", n])
            w.writerow([f"ge{t}_pct", f"{100.0*n/len(vals):.2f}"])
    print(f"[out] {sp}")

    # ---- 逐年 ----
    yr_rows = con.execute(f"""
        SELECT m.theyear AS yr, COUNT(*) AS n,
               CAST(quantile_cont(v.self_count, 0.5) AS BIGINT) AS med,
               CAST(quantile_cont(v.self_count, 0.9) AS BIGINT) AS p90,
               MAX(v.self_count) AS mx
        FROM read_parquet('{V}') v
        JOIN read_parquet('{M}') m ON v.acc = m.acc
        WHERE v.name = ? AND m.theyear IS NOT NULL AND m.theyear BETWEEN 2010 AND 2026
        GROUP BY 1 ORDER BY 1""", [TARGET]).fetchall()
    yp = os.path.join(args.outdir, "sars_spot_per_year.tsv")
    with open(yp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["year", "n_runs", "median_spot", "p90_spot", "max_spot"])
        w.writerows(yr_rows)
    print(f"[out] {yp}")

    # ---- 图 ----
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "figure.dpi": 150,
    })
    fig = plt.figure(figsize=(15, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.05], wspace=0.28)

    # 左：分档柱状图（对数 y）
    ax1 = fig.add_subplot(gs[0, 0])
    labels, counts = [], []
    for lo, hi in BINS:
        if hi is None:
            c = int((vals >= lo).sum()); lab = f"≥{lo:,}"
        elif hi == lo:
            c = int((vals == lo).sum()); lab = str(lo)
        else:
            c = int(((vals >= lo) & (vals <= hi)).sum()); lab = f"{lo:,}–{hi:,}"
        labels.append(lab); counts.append(c)
    colors = ["#B71C1C" if l.startswith(("1", "2")) and "–" not in l and l != "1"
              else "#EF9A9A" if l in ("2–4", "5–9") else "#64B5F6" if l in ("10–49", "50–99")
              else "#1E88E5" for l in labels]
    bars = ax1.bar(range(len(labels)), counts, color=colors, edgecolor="#37474F", linewidth=0.6)
    ax1.set_yscale("log")
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax1.set_ylabel("Runs (log scale)")
    ax1.set_xlabel("SARS-CoV-2 spot")
    ax1.set_title(f"Spot distribution (n = {n_runs:,} runs)", fontsize=10, pad=12)
    ax1.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    ax1.set_axisbelow(True)
    ax1.set_ylim(0.5, max(counts) * 12)
    for i, c in enumerate(counts):
        if c > 0:
            ax1.text(i, c * 1.3, f"{c:,}", ha="center", fontsize=7, rotation=90,
                     va="bottom", color="#263238")

    # 中：阈值存活曲线
    ax2 = fig.add_subplot(gs[0, 1])
    pcts = [100.0 * (vals >= t).sum() / len(vals) for t in THRESH]
    ax2.plot(THRESH, pcts, marker="o", color="#1565C0", linewidth=2,
             markersize=6, markeredgecolor="white", markeredgewidth=0.7)
    ax2.set_xscale("log")
    ax2.set_xlabel("Spot threshold")
    ax2.set_ylabel("% of runs passing")
    ax2.set_title("Survival curve", fontsize=10)
    ax2.grid(linestyle=":", linewidth=0.6, alpha=0.6)
    ax2.set_axisbelow(True)
    for t, p in zip(THRESH, pcts):
        if p > 2:
            ax2.annotate(f"{p:.1f}%", xy=(t, p), xytext=(0, 9),
                         textcoords="offset points", ha="center", fontsize=7.5,
                         color="#1565C0")
    ax2.set_ylim(-3, 112)
    # 标注关键阈值
    for t, note in ((100, "10.1%"), (1000, "3.4%")):
        ax2.axvline(t, color="#EF6C00", linestyle="--", linewidth=0.9, alpha=0.7)
        ax2.text(t, 50, f" {t}", rotation=90, fontsize=7.5, color="#EF6C00",
                 va="center")

    # 右：逐年
    ax3 = fig.add_subplot(gs[0, 2])
    if yr_rows:
        years = [r[0] for r in yr_rows]
        ns = [r[1] for r in yr_rows]
        meds = [r[2] if r[2] else 0 for r in yr_rows]
        colors3 = ["#C62828" if y >= 2020 else "#78909C" for y in years]
        ax3.bar(years, ns, color=colors3, edgecolor="#37474F", linewidth=0.5, alpha=0.85)
        ax3.set_ylabel("Runs with SARS-CoV-2", color="#37474F")
        ax3.set_xlabel("Release year")
        ax3.set_xticks(years[::2])
        ax3.set_xticklabels([str(y) for y in years[::2]], rotation=45, fontsize=8)
        ax3.set_title("Detections by year (red = 2020+)", fontsize=10)
        ax3.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
        ax3.set_axisbelow(True)
        ax3b = ax3.twinx()
        ax3b.plot(years, meds, marker="s", color="#1E88E5", linewidth=1.6,
                  markersize=4.5, label="median spot")
        ax3b.set_ylabel("Median spot", color="#1E88E5")
        ax3b.tick_params(axis="y", colors="#1E88E5")
        ax3b.spines["right"].set_visible(True)
        ax3b.spines["right"].set_color("#1E88E5")
        ax3.set_zorder(ax3b.get_zorder() + 1)
        ax3.patch.set_visible(False)

    fig.suptitle("SARS-CoV-2 detections in plant-declared SRA runs: "
                 "a contamination signature", fontsize=12, y=1.03)
    fig.subplots_adjust(left=0.06, right=0.94, top=0.84, bottom=0.22, wspace=0.32)
    png = os.path.join(args.outdir, "fig_sars_spot_dist.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(png.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[out] {png}")

    # ---- 控制台摘要 ----
    print("\n=== SARS-CoV-2 spot 分位 ===")
    print(f"  n = {n_runs:,} runs")
    for k in ("min", "p25", "median", "p75", "p90", "p99", "max"):
        v = {"min": vals[0], "p25": np.percentile(vals, 25), "median": np.median(vals),
             "p75": np.percentile(vals, 75), "p90": np.percentile(vals, 90),
             "p99": np.percentile(vals, 99), "max": vals[-1]}[k]
        print(f"  {k:8s} {v:>12,.0f}")
    print("\n=== 阈值通过 ===")
    for t in THRESH:
        n = int((vals >= t).sum())
        print(f"  >= {t:<6} {n:>8,} ({100.0*n/len(vals):5.2f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
