"""MCMC 收敛诊断绘图（合并 log 的 ESS 可视化）。

用法:
    python -m utils.mcmc_diag_plot --log merged/merged.log --outdir merged \
        [--burnin 0.1] [--log2 chain2.log] [--params treeModel.rootHeight,ucld.mean]

产出:
    mcmc_diagnostics.pdf  — 关键参数 trace + density + ESS 柱状图
    chain_scatter.pdf     — 双链收敛对比 (y=x 散点 + 相关系数 R, 可选)

ESS 算法: utils.ess (Geyer 1992 initial monotone sequence, 与 Tracer 一致)。
颜色约定同 Tracer: ESS>=200 绿, >=100 黄, <100 红。
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
    if ess >= 200:
        return "#7fbf7f"
    if ess >= 100:
        return "#eede91"
    return "#ee9191"


def plot_diagnostics(log_file, out_pdf, burnin_frac=0.1, params=None):
    header, cols = read_log_columns(log_file, burnin_frac)
    picked = pick_params(header, params)
    if not picked:
        raise ValueError("No plottable parameters found")
    n = len(picked)
    fig, axes = plt.subplots(n, 3, figsize=(16, 2.8 * n))
    if n == 1:
        axes = axes.reshape(1, 3)
    ess_values = []
    for i, p in enumerate(picked):
        v = cols[p]
        # trace
        ax = axes[i, 0]
        ax.plot(v, lw=0.4, color="#4e699a")
        burn_x = int(len(v) * burnin_frac)
        ax.axvline(burn_x, color="gray", ls="--", lw=0.8)
        ax.set_ylabel(p, fontsize=8)
        ax.set_title("Trace", fontsize=9, loc="left")
        # density: 直方图 + 95% HPD 阴影标注 + 统计盒
        ax = axes[i, 1]
        ax.hist(v, bins=60, density=True, color="#4e699a", alpha=0.75)
        lo, hi = _hdi(v, 0.95)
        ax.axvspan(lo, hi, color="red", alpha=0.1)
        for xv in (lo, hi):
            ax.axvline(xv, color="red", ls="--", lw=0.9)
        q05, q95 = np.quantile(v, [0.05, 0.95])
        ax.text(0.02, 0.97,
                f"mean={np.mean(v):.4g}\nmedian={np.median(v):.4g}\n"
                f"95% HPD=[{lo:.4g}, {hi:.4g}]\ntailESS={ess_tail(v):.0f}",
                transform=ax.transAxes, va="top", fontsize=7,
                bbox=dict(fc="white", ec="#999", alpha=0.85, pad=2))
        ax.set_title("Density (95% HPD)", fontsize=9, loc="left")
        # autocorrelation 前 100 lag
        ax = axes[i, 2]
        vc = v - v.mean()
        var = vc.var(ddof=1)
        if var > 0:
            lags = min(100, len(v) // 2)
            acf = [1.0] + [float(np.sum(vc[:-k] * vc[k:]) / (var * len(vc)))
                           for k in range(1, lags)]
            ax.bar(range(len(acf)), acf, color="#4e699a", width=0.9)
            ax.axhline(0, color="black", lw=0.5)
        ax.set_title("ACF", fontsize=9, loc="left")
        ess_values.append((p, float(calculate_ess(v)), len(v),
                           float(ess_bulk(v)), float(ess_tail(v))))

    # ESS 柱状图: bulk + tail 双柱 (Tracer 1.7+ 双指标)
    fig2, ax = plt.subplots(figsize=(9, 0.6 * len(ess_values) + 1.5))
    names = [e[0] for e in ess_values][::-1]
    bulk = [e[3] for e in ess_values][::-1]
    tail = [e[4] for e in ess_values][::-1]
    tot = [e[2] for e in ess_values][::-1]
    y = np.arange(len(names))
    ax.barh(y + 0.2, bulk, height=0.38, color=[_ess_color(v) for v in bulk],
            edgecolor="#333", label="ESS bulk")
    ax.barh(y - 0.2, tail, height=0.38, color=[_ess_color(v) for v in tail],
            edgecolor="#333", alpha=0.55, label="ESS tail")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    for j, (b, t, tt) in enumerate(zip(bulk, tail, tot)):
        ax.text(b, j + 0.2, f" {b:.0f}/{tt}", va="center", fontsize=7)
        ax.text(t, j - 0.2, f" {t:.0f}", va="center", fontsize=7)
    ax.axvline(200, color="green", ls=":", lw=0.8)
    ax.axvline(100, color="orange", ls=":", lw=0.8)
    ax.set_xlabel("ESS (green>=200, yellow>=100, red<100)")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title(f"ESS bulk/tail — {os.path.basename(log_file)} (burnin={burnin_frac:.0%})")
    fig2.tight_layout()
    fig2.savefig(out_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig2)

    # trace+density 也存一份
    fig.tight_layout()
    trace_pdf = out_pdf.replace(".pdf", "_trace.pdf")
    fig.savefig(trace_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig)

    # ESS evolution 面板: bulk/tail 随链长增长 (回答"链长够不够")
    evo_pdf = out_pdf.replace(".pdf", "_evolution.pdf")
    n_par = len(picked)
    fig3, axes3 = plt.subplots(n_par, 1, figsize=(8, 2.2 * n_par))
    if n_par == 1:
        axes3 = [axes3]
    for ax3, p in zip(axes3, picked):
        fr, bl, tl = ess_evolution(cols[p])
        ax3.plot(fr * len(cols[p]), bl, color="#4e699a", lw=1.5, label="ESS bulk")
        ax3.plot(fr * len(cols[p]), tl, color="#c46a6a", lw=1.5, ls="--", label="ESS tail")
        ax3.axhline(200, color="green", ls=":", lw=0.8)
        ax3.axhline(100, color="orange", ls=":", lw=0.8)
        ax3.set_ylabel(p, fontsize=8)
        ax3.legend(fontsize=7, loc="lower right")
        ax3.grid(alpha=0.3)
    axes3[-1].set_xlabel("Samples (cumulative)")
    fig3.suptitle("ESS evolution — curves should plateau above thresholds")
    fig3.tight_layout()
    fig3.savefig(evo_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig3)
    return ess_values, trace_pdf, evo_pdf


def plot_chain_scatter(log1, log2, out_pdf, burnin_frac=0.1, params=None):
    """YR-MPE tracecomp 式双链收敛对比: y=x 散点 + 相关系数。"""
    h1, c1 = read_log_columns(log1, burnin_frac)
    h2, c2 = read_log_columns(log2, burnin_frac)
    common = [k for k in pick_params(h1, params) if k in h2]
    if not common:
        common = [k for k in h1 if k in h2 and k.lower() not in ("state", "iterations")][:6]
    if not common:
        raise ValueError("No common numeric columns between the two logs")
    n = len(common)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.4))
    if n == 1:
        axes = [axes]
    rs = []
    for ax, p in zip(axes, common):
        # 长度对齐 (min)
        m = min(len(c1[p]), len(c2[p]))
        x, y = c1[p][:m], c2[p][:m]
        r = float(np.corrcoef(x, y)[0, 1])
        rs.append((p, r))
        lims = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.scatter(x, y, s=6, alpha=0.4, color="#4e699a", edgecolors="none")
        ax.plot(lims, lims, "k--", lw=0.8)
        ax.set_xlabel(f"chain1 {p}", fontsize=8)
        ax.set_ylabel(f"chain2 {p}", fontsize=8)
        ax.set_title(f"r={r:.3f}", fontsize=9)
    fig.suptitle(f"Chain convergence: {os.path.basename(log1)} vs {os.path.basename(log2)}")
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight", dpi=150)
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
