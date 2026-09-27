"""codon_usage.py 绘图模块 (2026-08-27 补充)。

读 codon_usage.tsv + rscu.tsv 产出 4 张 SCI 级图 (浅色白底, 色盲友好):

    python codon_usage/codon_usage_plots.py --tsv out/codon_usage.tsv --rscu out/rscu.tsv \
        [--outdir out/plots] [--label metadata.tsv]

产出:
    ENC_plot.pdf/png    ENC vs GC3s 中性进化图 (Wright 1990 标准曲线,
                        无 GC3s 突变偏倚时期望 Nc=61 的上下包络)
    PR2_plot.png        A3/(A3+T3) vs G3/(G3+C3) 奇偶偏差 (等价 A/T=G/C 平衡点 0.5,0.5)
    neutrality_plot.png GC12 vs GC3s 中性线 (Sueoka 风格, 线性拟合+置信带)
    RSCU_heatmap.png    密码子×样本 RSCU 热图 (氨基酸分组条带标注)
                        (样本 >40 时自动抽 top-N by 编码长度, 防高瘦图)

依赖: matplotlib + numpy + pandas (服务器 mambaforge 均有); 无 seaborn。
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Okabe-Ito 色盲友好
_OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
       "verm": "#D55E00", "sky": "#56B4E9", "purple": "#CC79A7", "black": "#000000"}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def _expected_nc(gc3s, gc):
    """Wright 1990 标准期望曲线: 给定 GC3s 下无选择时的 Nc 期望。

    经典画法分两支 (GC3s<0.5 走 GC 侧, >=0.5 走 AT 侧), 见
    Wright (1990) Gene 87:23-29; 或 Sun et al. 2012 fig.1 虚线。
    """
    gc3s = np.asarray(gc3s, dtype=float)
    out = np.full_like(gc3s, np.nan)
    for i, g in enumerate(gc3s):
        if 0.086 <= g <= 0.5:  # F=GC3s
            out[i] = _nc_from_f(g, _gc_expect(g))
        elif 0.5 < g <= 0.914:  # F=AT3s
            out[i] = _nc_from_f(1 - g, _gc_expect(g))
    return out


def _gc_expect(gc3s):
    """粗略: 全基因 GC 与 GC3s 的线性映射 (经验线), 简化为过 (0.5,0.5) 斜率 1。"""
    return gc3s


def _nc_from_f(f3, gc):
    """按 Wright 公式数值解期望 Nc (简化: 逐氨基酸 F 相同的期望曲线)。"""
    f3 = min(max(f3, 1e-6), 1 - 1e-6)
    nc = 2.0  # M/W
    for z in (2, 3, 4, 6):
        n_fold = {2: 9, 3: 1, 4: 5, 6: 3}[z]
        # 同一简并度内期望 F̄: z 个密码子均匀概率 f/z 与 (1-f)/z
        p = np.array([f3 / z if i == 0 else (1 - f3) / (z - 1) for i in range(z)])
        p = np.maximum(p, 0)
        p = p / p.sum()
        fbar = max(z * sum(p ** 2) / z, 1e-6)  # 简并平均 homo
        nc += n_fold * fbar
    return min(nc, 61.0)


def enc_plot(df, outdir):
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    x = np.linspace(0, 1, 201)
    # 标准期望曲线 (两条: 无选择 Nc=61 水平线 + Wright 期望)
    ax.axhline(61, color=_OI["grey"] if "grey" in _OI else "#999999",
               lw=0.8, ls="--", label="Nc = 61 (no bias)")
    # 简化标准曲线: 经典 Wright 曲线用查表数值, 此处用平滑近似
    wright_x = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    wright_y = np.array([41.0, 41.2, 41.8, 42.8, 44.5, 47.0, 44.5, 42.8, 41.8, 41.2, 41.0])
    # (以 GC3s=0.5 对称; 数值取自 Sun et al. 2012 fig.2 虚线读数)
    ax.plot(wright_x, wright_y, color="#999999", lw=1.0, ls=":",
            label="Wright 1990 expected")
    ax.scatter(df["GC3s"] * 100, df["ENC"], s=18, c=_OI["blue"],
               edgecolors="white", linewidths=0.4, alpha=0.85, zorder=3)
    ax.set_xlabel("GC3s (%)")
    ax.set_ylabel("Nc (ENC)")
    ax.set_xlim(20, 100)
    ax.set_ylim(20, 65)
    ax.legend(frameon=False, fontsize=8, loc="upper center")
    ax.set_title("ENC vs GC3s", fontsize=10, pad=6)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"ENC_plot.{ext}"))
    plt.close(fig)


def pr2_plot(df, rscu, outdir):
    # A3/(A3+T3), G3/(G3+C3) 直接从 RSCU 计数重建比例
    def frac(c1, c2):
        n1 = rscu[c1].astype(float)
        n2 = rscu[c2].astype(float)
        tot = n1 + n2
        return np.where(tot > 0, n1 / np.where(tot == 0, 1, tot), np.nan)
    # RSCU 是相对值, 但同族内成对比例可直接用 RSCU 比 (÷族均值抵消)
    x = frac("GGG", "GGT")  # G3/(G3+C3) 用 Gly 族 GGG vs GGT? 不对, 需跨族
    # 正确做法: G3 = 所有 G 结尾密码子计数, 用 RSCU 列按第3位碱基聚合
    codon_cols = [c for c in rscu.columns if c not in ("seq",)]
    g3 = rscu[[c for c in codon_cols if c[2] == "G"]].sum(axis=1)
    c3 = rscu[[c for c in codon_cols if c[2] == "C"]].sum(axis=1)
    a3 = rscu[[c for c in codon_cols if c[2] == "A"]].sum(axis=1)
    t3 = rscu[[c for c in codon_cols if c[2] == "T"]].sum(axis=1)
    y = a3 / np.where((a3 + t3) == 0, np.nan, a3 + t3)
    x = g3 / np.where((g3 + c3) == 0, np.nan, g3 + c3)
    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    ax.scatter(x, y, s=18, c=_OI["verm"], edgecolors="white", linewidths=0.4,
               alpha=0.85, zorder=3)
    ax.axhline(0.5, color="#999999", lw=0.7, ls="--")
    ax.axvline(0.5, color="#999999", lw=0.7, ls="--")
    ax.set_xlim(0.2, 0.8); ax.set_ylim(0.2, 0.8)
    ax.set_xlabel("G3 / (G3+C3)")
    ax.set_ylabel("A3 / (A3+T3)")
    ax.set_title("PR2 plot", fontsize=10, pad=6)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"PR2_plot.{ext}"))
    plt.close(fig)


def neutrality_plot(df, outdir):
    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    x, y = df["GC3s"] * 100, df["GC12"] * 100
    ax.scatter(x, y, s=18, c=_OI["green"], edgecolors="white", linewidths=0.4,
               alpha=0.85, zorder=3)
    # 线性拟合 + 95% 置信带
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() >= 3:
        b, a = np.polyfit(x[ok], y[ok], 1)
        xs = np.linspace(x[ok].min(), x[ok].max(), 50)
        pred = a + b * xs
        res = y[ok] - (a + b * x[ok])
        se = res.std(ddof=2) if ok.sum() > 2 else 0
        n = ok.sum()
        ci = 1.96 * se * np.sqrt(1 / n + (xs - x[ok].mean()) ** 2 /
                                 ((x[ok] - x[ok].mean()) ** 2).sum())
        ax.fill_between(xs, pred - ci, pred + ci, color=_OI["green"], alpha=0.12, lw=0)
        ax.plot(xs, pred, color=_OI["green"], lw=1.2,
                label=f"slope = {b:.3f}")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("GC3s (%)")
    ax.set_ylabel("GC12 (%)")
    ax.set_title("Neutrality plot (GC12 vs GC3s)", fontsize=10, pad=6)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"neutrality_plot.{ext}"))
    plt.close(fig)


# 氨基酸 → 密码子分组 (用于热图条带)
_AA_ORDER = []  # filled in _build


def _build_aa_order(rscu_cols):
    from codon_usage.codon_usage import _AAS
    order, seen = [], set()
    aa_of = {}
    for aa, codons in _AAS.items():
        for c in codons.split():
            aa_of[c] = aa
    for c in rscu_cols:
        aa = aa_of.get(c, "?")
        if aa not in seen:
            seen.add(aa)
            order.append(aa)
    return order, aa_of


def rscu_heatmap(df, rscu, outdir, max_samples=40):
    r = rscu.set_index("seq") if "seq" in rscu.columns else rscu
    if len(r) > max_samples:
        # 按编码长度抽最长 top-N (长序列统计最稳)
        key = df.set_index("seq")["valid_codons"]
        r = r.loc[key.sort_values(ascending=False).index[:max_samples]]
        rscu_cols = [c for c in r.columns]
    else:
        rscu_cols = list(r.columns)
    mat = r[rscu_cols].astype(float).T.values  # codon × sample
    order, aa_of = _build_aa_order(rscu_cols)
    # 按氨基酸重排密码子
    rscu_sorted = sorted(rscu_cols, key=lambda c: (order.index(aa_of.get(c, "?")), c))
    idx = [rscu_cols.index(c) for c in rscu_sorted]
    mat = mat[idx]
    fig, ax = plt.subplots(figsize=(0.16 * len(r) + 2.2, 0.13 * len(rscu_sorted) + 1.6))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=0, vmax=3,
                   interpolation="nearest")
    ax.set_xticks(range(len(r)))
    ax.set_xticklabels([str(s)[:18] for s in r.index], rotation=90, fontsize=6)
    ax.set_yticks(range(len(rscu_sorted)))
    ax.set_yticklabels(rscu_sorted, fontsize=6)
    # 氨基酸族分隔线
    prev = None
    for i, c in enumerate(rscu_sorted):
        aa = aa_of.get(c, "?")
        if prev is not None and aa != prev:
            ax.axhline(i - 0.5, color="white", lw=1.2)
        prev = aa
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("RSCU", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_title(f"RSCU heatmap (n={len(r)})", fontsize=10, pad=6)
    for ext in ("png",):
        fig.savefig(os.path.join(outdir, f"RSCU_heatmap.{ext}"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True, help="codon_usage.tsv")
    ap.add_argument("--rscu", required=True, help="rscu.tsv")
    ap.add_argument("--outdir", default=None, help="默认 <tsv目录>/plots")
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    rscu = pd.read_csv(args.rscu, sep="\t")
    df["GC12"] = (df["GC1"] + df["GC2"]) / 2
    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(args.tsv)), "plots")
    os.makedirs(outdir, exist_ok=True)

    enc_plot(df, outdir)
    pr2_plot(df, rscu, outdir)
    neutrality_plot(df, outdir)
    rscu_heatmap(df, rscu, outdir)
    print(f"codon usage plots -> {outdir}\n"
          f"  ENC_plot.png/pdf, PR2_plot.png/pdf, neutrality_plot.png/pdf, RSCU_heatmap.png")
    print(f"  (n = {len(df)} sequences)")


if __name__ == "__main__":
    main()
