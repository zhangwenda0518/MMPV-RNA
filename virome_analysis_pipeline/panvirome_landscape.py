#!/usr/bin/env python3
"""
panvirome_landscape.py — 泛病毒全基因组进化景观构建器
=======================================================
从 panvirome_aggregator.py 的输出中读取每病毒的合并数据，
生成三轨全景图: 基因注释 + SNP/π 密度 + 基因级 dN/dS。

输入:
  --merged-dir  Aggregator 输出目录 (每病毒子目录含 gene_coordinates.csv, merged_site_pi.csv, merged_product_results.csv)
输出:
  每病毒子目录下: Plot6_Super_Evolutionary_Landscape.pdf / .png
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator

# 出版级配置
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["axes.linewidth"] = 1.3

# 配色
COLOR_SNP = "#41B6C4"
COLOR_PI = "#737373"
COLOR_DNDS = "#E31A1C"
COLOR_PURIFY = "#74C476"
NEUTRAL_LINE = "#969696"
GENE_COLORS = ["#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99"]
WINDOW_SIZE = 50


def safe_read(filepath):
    if not os.path.exists(filepath):
        return pd.DataFrame()
    try:
        return pd.read_csv(filepath)
    except Exception as e:
        print(f"  [WARN] {filepath}: {e}")
        return pd.DataFrame()


def parse_args():
    p = argparse.ArgumentParser(description="泛病毒进化景观构建")
    p.add_argument("--merged-dir", required=True,
                   help="Aggregator 输出目录 (每病毒子目录)")
    p.add_argument("--virus", default=None,
                   help="仅处理指定病毒 (默认: 全部)")
    return p.parse_args()


def main():
    args = parse_args()
    base = args.merged_dir
    if not os.path.isdir(base):
        print(f"[ERROR] 目录不存在: {base}")
        sys.exit(1)

    subdirs = [os.path.join(base, d) for d in os.listdir(base)
               if os.path.isdir(os.path.join(base, d))]
    if args.virus:
        subdirs = [d for d in subdirs if args.virus in os.path.basename(d)]

    print(f"  共 {len(subdirs)} 个病毒")

    for current_dir in subdirs:
        virus_name = os.path.basename(current_dir)
        print(f"\n  ▶ {virus_name}")

        gene_file = os.path.join(current_dir, "gene_coordinates.csv")
        site_file = os.path.join(current_dir, "merged_site_pi.csv")
        prod_file = os.path.join(current_dir, "merged_product_results.csv")

        genes_df = safe_read(gene_file)
        site_df = safe_read(site_file)
        prod_df = safe_read(prod_file)

        if genes_df.empty:
            print("    [SKIP] 无基因坐标 (类病毒/非编码)")
            continue
        if site_df.empty or prod_df.empty:
            print("    [SKIP] 缺 π 或 dN/dS 数据")
            continue

        # 基因级 dN/dS
        prod_summary = prod_df.groupby("product")[["piN", "piS"]].mean().reset_index()
        prod_summary["dNdS"] = prod_summary["piN"] / (prod_summary["piS"] + 1e-5)
        max_dnds = prod_summary["dNdS"].max()
        ylim_top = max(1.3, min(max_dnds * 1.3, 10.5))

        # 全基因组阵列
        max_len = int(max(genes_df["end"].max(), site_df["site"].max())) + 100
        pi_array = np.zeros(max_len)
        snp_array = np.zeros(max_len)

        for _, row in site_df.iterrows():
            pos = int(row["site"]) - 1
            if pos < max_len:
                pi_array[pos] = row["pi"]
                snp_array[pos] = 1

        x_coords = np.arange(1, max_len + 1)
        roll_snps = pd.Series(snp_array).rolling(window=WINDOW_SIZE, center=True).sum().fillna(0)
        roll_pi = pd.Series(pi_array).rolling(window=WINDOW_SIZE, center=True).mean().fillna(0)

        # 三轨画布
        fig, (ax_struct, ax_div, ax_sel) = plt.subplots(
            3, 1, figsize=(16, 10), sharex=True,
            gridspec_kw={"height_ratios": [1, 3.5, 2]})
        fig.subplots_adjust(hspace=0.1)

        # 轨1: 基因注释
        ax_struct.set_xlim(0, max_len)
        ax_struct.set_ylim(-0.8, 0.8)
        ax_struct.plot([0, max_len], [0, 0], color="black", linewidth=1.5, zorder=1)
        for idx, row in genes_df.iterrows():
            gs, ge, gn = row["start"], row["end"], row["gene"]
            gw = ge - gs
            c = GENE_COLORS[idx % len(GENE_COLORS)]
            rect = patches.Rectangle((gs, -0.4), gw, 0.8, facecolor=c, edgecolor="black", linewidth=1.2, zorder=2)
            ax_struct.add_patch(rect)
            ax_struct.text(gs + gw / 2, 0, gn, ha="center", va="center", fontweight="bold", fontsize=12)
        ax_struct.axis("off")
        ax_struct.set_title(f"Genomic & Evolutionary Pan-Landscape: {virus_name}", fontsize=18, fontweight="bold", pad=15)

        # 轨2: SNP + π
        ax_div.plot(x_coords, roll_snps, color=COLOR_SNP, alpha=0.9, linewidth=1.8)
        ax_div.fill_between(x_coords, roll_snps, color=COLOR_SNP, alpha=0.25)
        ax_div.set_ylabel(f"SNPs / {WINDOW_SIZE}nt", color=COLOR_SNP, fontsize=13, fontweight="bold")
        ax_div.tick_params(axis="y", labelcolor=COLOR_SNP)
        ax_div.set_ylim(bottom=0)

        ax_twin = ax_div.twinx()
        ax_twin.plot(x_coords, roll_pi, color=COLOR_PI, linestyle="--", linewidth=1.5)
        ax_twin.set_ylabel("Nucleotide Diversity (π)", color=COLOR_PI, fontsize=13, fontweight="bold")
        ax_twin.tick_params(axis="y", labelcolor=COLOR_PI)
        ax_twin.set_ylim(bottom=0)

        # 轨3: dN/dS
        ax_sel.axhline(1.0, color=NEUTRAL_LINE, linestyle="--", linewidth=1.2, zorder=1)
        for idx, row in genes_df.iterrows():
            gs, ge, gn = row["start"], row["end"], row["gene"]
            gw = ge - gs
            dnds_row = prod_summary.loc[prod_summary["product"].str.contains(gn, na=False, case=False)]
            if not dnds_row.empty:
                dnds_val = min(dnds_row["dNdS"].values[0], 10.0)
                color = COLOR_DNDS if dnds_val > 1.0 else COLOR_PURIFY
                ax_sel.bar(gs + gw / 2, dnds_val, width=gw, facecolor=color, alpha=0.85, edgecolor="black", linewidth=1.0, zorder=2)
                ax_sel.text(gs + gw / 2, dnds_val + ylim_top * 0.04, f"{dnds_val:.3f}",
                            ha="center", va="bottom", fontsize=10, fontweight="bold", color="black")

        ax_sel.set_ylabel("Orthologous dN/dS", fontsize=13, fontweight="bold")
        ax_sel.set_xlabel("Genomic Coordinates (bp)", fontsize=15, fontweight="bold")
        ax_sel.set_ylim(0, ylim_top)
        ax_sel.spines["top"].set_visible(False)
        ax_sel.spines["right"].set_visible(False)
        ax_sel.xaxis.set_major_locator(MaxNLocator(integer=True, prune="both", nbins=12))

        out_pdf = os.path.join(current_dir, "Plot6_Super_Evolutionary_Landscape.pdf")
        out_png = os.path.join(current_dir, "Plot6_Super_Evolutionary_Landscape.png")
        plt.tight_layout()
        plt.savefig(out_pdf, bbox_inches="tight", dpi=300)
        plt.savefig(out_png, bbox_inches="tight", dpi=400)
        plt.close()
        print(f"    → {os.path.basename(out_pdf)}")

    print(f"\n  完成")


if __name__ == "__main__":
    main()
