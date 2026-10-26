"""MCMC 收敛诊断绘图（合并 log 的 ESS 可视化）。

用法:
    python -m utils.mcmc_diag_plot --log merged/merged.log --outdir merged \
        [--burnin 0.1] [--log2 chain2.log] [--params treeModel.rootHeight,ucld.mean]

产出:
    mcmc_diagnostics.pdf  — 关键参数 trace + density + ESS 柱状图
    chain_scatter.pdf     — 双链收敛对比 (y=x 散点 + 相关系数 R, 可选)

ESS 算法: utils.ess (Geyer 1992 initial monotone sequence, 与 Tracer 一致)。
颜色约定同 Tracer 分级、配色用项目 OKABE: ESS>=200 绿, >=100 橙, <100 朱红。
"""
import os
import csv
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.ess import calculate_ess
from utils.ess_metrics import ess_bulk, ess_tail, ess_evolution

# 2026-10-12: 统一 SCI 配色/样式。脚本既能 `import utils.mcmc_diag_plot`，
# 也能被直接 `python utils/mcmc_diag_plot.py` 执行（此时 sys.path[0]=utils/，
# 包路径不存在）→ 双路径 import。
try:
    from utils.sci_palette import (
        apply_sci_style, OKABE_ITO, COL_W_1, COL_W_2, despine, sci_legend)
    # 2026-10-16: 图件双格式 (PNG+PDF 同主名)
    from utils.plot_io import save_both
except ImportError:                                    # pragma: no cover
    from sci_palette import (
        apply_sci_style, OKABE_ITO, COL_W_1, COL_W_2, despine, sci_legend)
    from plot_io import save_both


def _hdi(x, prob=0.95):
    """最高密度区间 (Tracer/arviz 口径, 与 utils.ess.calculate_95hpd 同源; ceil 口径)。"""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    k = max(int(np.ceil(prob * n)), 1)
    if n <= k:
        return float(x[0]), float(x[-1])
    widths = x[k:] - x[:n - k]
    i = int(np.argmin(widths))
    return float(x[i]), float(x[i + k - 1])


def read_log_columns(log_file, burnin_frac=0.1):
    """读 BEAST tab 分隔 log, 返回 (列名, {col: post-burnin np.array})。"""
    rows = []
    header = None
    with open(log_file, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = [p.strip() for p in parts]
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue  # 注释/非数据行
    if not header or not rows:
        raise ValueError(f"No data in {log_file}")
    n_total = len(rows)
    n_skip = int(n_total * float(burnin_frac))
    data = np.array(rows[n_skip:])
    cols = {}
    for i, name in enumerate(header):
        try:
            cols[name] = data[:, i]
        except IndexError:
            pass
    return header, cols


def pick_params(header, params_arg=None, max_params=6):
    """优先用户指定; 否则选常见关键参数, 再补方差最大的列。"""
    if params_arg:
        wanted = [p.strip() for p in params_arg.split(",") if p.strip()]
        return [p for p in wanted if p in header][:max_params]
    key_candidates = [
        "treeModel.rootHeight", "treeHeight", "height",
        "ucld.mean", "clock.rate", "branchRates",
        "constant.popSize", "skyline.popSize", "b_POPSize",
        "kappa", "location.rates",
    ]
    picked = [k for k in key_candidates if k in header]
    # 历史坑: beast1 实际列名是 skyline.popSize1..N (带后缀), 精确匹配永不命中 → skyline 系列不进 ESS 图
    if not any(k.startswith("skyline.popSize") for k in picked):
        picked += [c for c in header if c.startswith("skyline.popSize")][:3]
    if not picked:
        # 方差最大的数值列 (排除第一列 state)
        return []
    return picked[:max_params]


def _ess_color(ess):
    """ESS 分级配色（Tracer 口径 → 项目 OKABE 语义色）。

    2026-10-12: 原为手写 #7fbf7f/#eede91/#ee9191（淡绿/淡黄/淡粉），
    在 300 dpi 印刷下饱和度不足、且不属于项目调色板。现改用
    OKABE 绿(#009E73)/橙(#E69F00)/朱红(#D55E00)。
    """
    if ess >= 200:
        return OKABE_ITO[2]      # bluish green — 充分
    if ess >= 100:
        return OKABE_ITO[0]      # orange       — 偏低
    return OKABE_ITO[5]          # vermillion   — 不足


def plot_diagnostics(log_file, out_pdf, burnin_frac=0.1, params=None):
    # 2026-10-12: 统一 SCI 样式 + 画布按**排版栏宽**（Fig 3 里 N/S 面板占 2 栏 = 4.58 in）。
    # 旧画布固定 16 in 宽 → 排进 4.58 in 只剩 0.29 倍，9 pt 标题变 2.6 pt，完全不可读。
    apply_sci_style(base_size=6.5)
    header, cols = read_log_columns(log_file, burnin_frac)
    picked = pick_params(header, params)
    if not picked:
        raise ValueError("No plottable parameters found")
    n = len(picked)
    fig, axes = plt.subplots(n, 3, figsize=(COL_W_2, COL_W_2 * 0.185 * n))
    if n == 1:
        axes = axes.reshape(1, 3)
    ess_values = []
    for i, p in enumerate(picked):
        v = cols[p]
        # trace
        ax = axes[i, 0]
        ax.plot(v, lw=0.35, color=OKABE_ITO[4])
        burn_x = int(len(v) * burnin_frac)
        ax.axvline(burn_x, color=OKABE_ITO[8], ls="--", lw=0.7)
        ax.set_ylabel(p, fontsize=5.5)
        ax.set_title("Trace", fontsize=6.5, loc="left")
        despine(ax)
        # density: 直方图 + 95% HPD 阴影标注 + 统计盒
        ax = axes[i, 1]
        ax.hist(v, bins=60, density=True, color=OKABE_ITO[4], alpha=0.75)
        lo, hi = _hdi(v, 0.95)
        ax.axvspan(lo, hi, color=OKABE_ITO[5], alpha=0.10)
        for xv in (lo, hi):
            ax.axvline(xv, color=OKABE_ITO[5], ls="--", lw=0.8)
        ax.text(0.02, 0.97,
                f"mean={np.mean(v):.4g}\nmedian={np.median(v):.4g}\n"
                f"95% HPD=[{lo:.4g}, {hi:.4g}]\ntailESS={ess_tail(v):.0f}",
                transform=ax.transAxes, va="top", fontsize=5,
                bbox=dict(fc="white", ec=OKABE_ITO[8], alpha=0.85, pad=1.4))
        ax.set_title("Density (95% HPD)", fontsize=6.5, loc="left")
        despine(ax)
        # autocorrelation 前 100 lag
        ax = axes[i, 2]
        vc = v - v.mean()
        var = vc.var(ddof=1)
        if var > 0:
            lags = min(100, len(v) // 2)
            acf = [1.0] + [float(np.sum(vc[:-k] * vc[k:]) / (var * len(vc)))
                           for k in range(1, lags)]
            ax.bar(range(len(acf)), acf, color=OKABE_ITO[4], width=0.9)
            ax.axhline(0, color="#222222", lw=0.5)
        ax.set_title("ACF", fontsize=6.5, loc="left")
        despine(ax)
        ess_values.append((p, float(calculate_ess(v)), len(v),
                           float(ess_bulk(v)), float(ess_tail(v))))

    # ESS 柱状图: bulk + tail 双柱 (Tracer 1.7+ 双指标)
    fig2, ax = plt.subplots(figsize=(COL_W_2, 0.22 * len(ess_values) + 0.95))
    names = [e[0] for e in ess_values][::-1]
    bulk = [e[3] for e in ess_values][::-1]
    tail = [e[4] for e in ess_values][::-1]
    tot = [e[2] for e in ess_values][::-1]
    y = np.arange(len(names))
    ax.barh(y + 0.2, bulk, height=0.38, color=[_ess_color(v) for v in bulk],
            edgecolor="#333333", linewidth=0.4, label="ESS bulk")
    ax.barh(y - 0.2, tail, height=0.38, color=[_ess_color(v) for v in tail],
            edgecolor="#333333", linewidth=0.4, alpha=0.55, label="ESS tail")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=5.5)
    for j, (b, t, tt) in enumerate(zip(bulk, tail, tot)):
        ax.text(b, j + 0.2, f" {b:.0f}/{tt}", va="center", fontsize=4.8)
        ax.text(t, j - 0.2, f" {t:.0f}", va="center", fontsize=4.8)
    ax.axvline(200, color=OKABE_ITO[2], ls=":", lw=0.8)
    ax.axvline(100, color=OKABE_ITO[0], ls=":", lw=0.8)
    ax.set_xlabel("ESS   (green $\\geq$200, orange $\\geq$100, red <100)", fontsize=6)
    # 图例用**中性灰** Patch：柱子颜色编码的是 ESS 分级，若图例沿用柱色会被误读成
    # “红=bulk / 浅红=tail”。分级含义已写在 x 轴标签里。
    from matplotlib.patches import Patch
    _h = [Patch(facecolor="#888888", edgecolor="#333333", linewidth=0.4,
                label="ESS bulk"),
          Patch(facecolor="#888888", edgecolor="#333333", linewidth=0.4,
                alpha=0.55, label="ESS tail")]
    sci_legend(ax, handles=_h, loc="lower right", fontsize=5.5)
    ax.set_title(f"ESS bulk/tail — {os.path.basename(log_file)} "
                 f"(burnin={burnin_frac:.0%})", fontsize=6.5)
    despine(ax)
    fig2.tight_layout()
    save_both(fig2, out_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig2)

    # trace+density 也存一份
    fig.tight_layout()
    trace_pdf = out_pdf.replace(".pdf", "_trace.pdf")
    save_both(fig, trace_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)

    # ESS evolution 面板: bulk/tail 随链长增长 (回答"链长够不够")
    evo_pdf = out_pdf.replace(".pdf", "_evolution.pdf")
    n_par = len(picked)
    fig3, axes3 = plt.subplots(n_par, 1, figsize=(COL_W_2, 0.46 * n_par))
    if n_par == 1:
        axes3 = [axes3]
    for ax3, p in zip(axes3, picked):
        fr, bl, tl = ess_evolution(cols[p])
        ax3.plot(fr * len(cols[p]), bl, color=OKABE_ITO[4], lw=1.2, label="ESS bulk")
        ax3.plot(fr * len(cols[p]), tl, color=OKABE_ITO[5], lw=1.2, ls="--",
                 label="ESS tail")
        ax3.axhline(200, color=OKABE_ITO[2], ls=":", lw=0.8)
        ax3.axhline(100, color=OKABE_ITO[0], ls=":", lw=0.8)
        ax3.set_ylabel(p, fontsize=5.5)
        sci_legend(ax3, loc="lower right", fontsize=5)
        ax3.grid(alpha=0.3)
        despine(ax3)
    axes3[-1].set_xlabel("Samples (cumulative)", fontsize=6)
    fig3.suptitle("ESS evolution — curves should plateau above thresholds",
                  fontsize=7, fontweight="bold")
    fig3.tight_layout()
    save_both(fig3, evo_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig3)
    return ess_values, trace_pdf, evo_pdf


def plot_chain_scatter(log1, log2, out_pdf, burnin_frac=0.1, params=None):
    """链间一致性 **QQ 图**: 两条链同参数的分位数是否落在 y=x 上。

    2026-10-12 修正
    --------------
    原实现是「按行号逐样本对齐的散点 + Pearson r」。但 BEAST 各链**独立采样**，
    相同 state 的样本之间**没有对应关系** → 点云天然是两团云而非沿对角线，
    r 恒被稀释到 ≈0（GCVA 实测 0.22 / −0.21 / 0.07 …），读者会把它误读成
    "两链相关性低"，而它其实什么也没说明。
    改为分位数-分位数对比：两条链的边际后验一致 → 点落 y=x；系统性偏离 → 未收敛。
    判据量取 **归一化最大偏差** max|Q1−Q2| / 数据跨度（0 = 完全一致）。
    """
    apply_sci_style(base_size=6.5)
    h1, c1 = read_log_columns(log1, burnin_frac)
    h2, c2 = read_log_columns(log2, burnin_frac)
    common = [k for k in pick_params(h1, params) if k in h2]
    if not common:
        common = [k for k in h1 if k in h2 and k.lower() not in ("state", "iterations")][:6]
    if not common:
        raise ValueError("No common numeric columns between the two logs")
    n = len(common)
    # 2026-10-12: 原为 1×n 横排（n=6 → 每格 0.53 in，标签全糊）。
    # 改为 3 列网格；画布 = Fig 3 的**通栏**宽度（9.32 in，P 面板独占一行）。
    ncol = min(3, n)
    nrow = -(-n // ncol)
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(COL_W_1, COL_W_1 * 0.30 * nrow))
    axes = np.atleast_1d(axes).ravel()
    rs = []
    for ax, p in zip(axes, common):
        x, y = np.sort(c1[p]), np.sort(c2[p])
        m = max(min(len(x), len(y)), 2)
        q = np.linspace(0.0, 1.0, m)
        xs, ys = np.quantile(x, q), np.quantile(y, q)
        lo = min(xs.min(), ys.min())
        hi = max(xs.max(), ys.max())
        span = (hi - lo) or 1.0
        dev = float(np.max(np.abs(xs - ys)) / span)
        rs.append((p, dev))
        ax.scatter(xs, ys, s=2.6, alpha=0.5, color=OKABE_ITO[4], edgecolors="none")
        ax.plot([lo, hi], [lo, hi], "--", color="#222222", lw=0.7)
        ax.set_xlabel(f"chain 1  {p}", fontsize=6.5)
        ax.set_ylabel(f"chain 2  {p}", fontsize=6.5)
        # 2026-10-12（审稿人视角）：dev 是**已归一化**的量（max|Δ分位数| ÷ 数据跨度），
        # 旧标签只写 "max dev = 0.645" 读者无从判断量纲，也没法和别的参数比。
        # 现写明"占该参数数据跨度的百分比"。
        ax.set_title(f"max deviation = {dev * 100:.1f}% of range", fontsize=7)
        # 刻度过密会互相压字（ucld.mean 一档 7 位小数）→ 限 3 个刻度
        ax.locator_params(axis="x", nbins=3)
        ax.locator_params(axis="y", nbins=4)
        despine(ax)
    for ax in axes[len(common):]:
        ax.axis("off")
    fig.suptitle(f"Chain agreement (QQ): {os.path.basename(log1)} vs {os.path.basename(log2)}"
                 f"   — points on the dashed line = identical marginal posteriors;"
                 f" deviation is relative to each parameter's own range",
                 fontsize=8, fontweight="bold")
    fig.tight_layout()
    save_both(fig, out_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return rs, out_pdf


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--log", required=True, help="合并后 BEAST log (tab 分隔)")
    ap.add_argument("--log2", default=None, help="第二条链 log (画 y=x 收敛对比)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--burnin", type=float, default=0.1)
    ap.add_argument("--params", default=None, help="逗号分隔参数名; 缺省自动选关键参数")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_pdf = os.path.join(args.outdir, "mcmc_diagnostics.pdf")
    ess_values, trace_pdf, evo_pdf = plot_diagnostics(args.log, out_pdf, args.burnin, args.params)
    print(f"ESS summary ({os.path.basename(args.log)}):")
    for p, ess, tot, eb, et in ess_values:
        flag = "OK" if (eb >= 200 and et >= 200) else ("LOW" if (eb >= 100 and et >= 100) else "FAIL")
        print(f"  {p:32s} bulk={eb:7.0f} tail={et:7.0f} / {tot}  [{flag}]")
    print(f"plots: {out_pdf} + {trace_pdf} + {evo_pdf}")

    if args.log2 and os.path.exists(args.log2):
        rs, sc_pdf = plot_chain_scatter(
            args.log, args.log2,
            os.path.join(args.outdir, "chain_scatter.pdf"), args.burnin, args.params)
        print(f"chain scatter: {sc_pdf} (r: " + ", ".join(f"{p}={r:.3f}" for p, r in rs) + ")")


if __name__ == "__main__":
    main()
