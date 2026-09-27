#!/usr/bin/env python3
"""
publication_figures.py — 发表级三件套：TempEst 图 + MCC 树 + 附表
==============================================================
生成 SCI 期刊标准的三张关键图和一份附表：

  Fig 1: TempEst root-to-tip 回归图 (时间信号)
  Fig 2: MCC 时间树 (按地点着色 + 节点后验概率)
  Table S1: 参数估计表 (mean, 95% HPD, ESS)

用法:
  python publication_figures.py \
      --rtt_fasta mafft.aln.fasta --rtt_dates dates.csv --rtt_tree tree.nwk \
      --mcc_tree mcc.tree --mcc_locations Ningxia,Beijing,Neimenggu \
      --beast_log beast.log \
      --output figures/
"""

import csv, json, os, sys, re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle  # unused but kept for reference
from matplotlib.lines import Line2D

import matplotlib.patheffects as pe

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector
from utils.ess import calculate_ess, calculate_95hpd

# ═══════════════════════════════════════════════════════════════════
# Color palette (Nature-style)
# ═══════════════════════════════════════════════════════════════════

NATURE_COLORS = ['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F',
                 '#8491B4','#91D1C2','#DC0000','#7E6148']

# ═══════════════════════════════════════════════════════════════════
# Fig 1: TempEst Root-to-Tip Regression
# ═══════════════════════════════════════════════════════════════════

def plot_temprest(
    fasta_file: str,
    dates_csv: str,
    output_path: str,
    tree_file: Optional[str] = None,
    title: str = "Root-to-Tip Regression",
    figsize: Tuple[float, float] = (10, 7),
    log: Optional[LogCollector] = None,
):
    """
    生成 TempEst 风格的 root-to-tip 回归图。

    直接计算遗传距离 vs 采样日期的线性回归，不依赖 BEAST。
    使用 Biopython 距离计算而非 TreeTime 内部方法（避免递归问题）。

    Returns dict with R², slope, intercept, p_value
    """
    from Bio import AlignIO, Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from scipy import stats as scipy_stats

    log = log or LogCollector()

    # Load data
    aln = AlignIO.read(fasta_file, "fasta")

    # Load dates
    dates = {}
    with open(dates_csv, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('name', '').strip()
            date = row.get('date', '').strip()
            if name and date:
                try:
                    dates[name] = float(date.replace('-', '.').split('.')[0])
                except ValueError:
                    pass

    # Build or load tree
    if tree_file and os.path.exists(tree_file):
        tree = Phylo.read(tree_file, "newick")
    else:
        # Build NJ tree from alignment
        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

    # Compute root-to-tip distances
    tree.root_with_outgroup(tree.get_terminals()[0])
    # 2026-09-16 修: `Bio.Phylo.depths()` 返回的 dict **以 Clade 对象为键**(不是名字) ——
    # 原写法用字符串 `rtt_dists[n]` 查, 命中率恒为 0 → `valid_tips` 恒空 →
    # 这张论文级 Fig 1 (TempEst root-to-tip 回归图) **一直静默不产出**
    # (只在日志里留一句 "Only 0 dated tips")。改为先建 name→depth 映射。
    # 注: 主树加了分支支持值后, 内部节点也有 name, 故只取 terminal。
    rtt_dists = tree.depths(unit_branch_lengths=True)
    rtt_by_name = {c.name: d for c, d in rtt_dists.items() if c.is_terminal() and c.name}

    # Match tips with dates
    tip_names = [c.name for c in tree.get_terminals()]
    valid_tips = [(n, rtt_by_name[n], dates[n]) for n in tip_names
                  if n in rtt_by_name and n in dates and dates[n] is not None]

    if len(valid_tips) < 3:
        log.warning(f"Only {len(valid_tips)} dated tips, cannot plot")
        return

    names, dists, times = zip(*valid_tips)
    dists = np.array(dists)
    times = np.array(times)

    # Linear regression
    slope, intercept, r, p, se = scipy_stats.linregress(times, dists)
    r2 = r ** 2
    tmrca = -intercept / slope if slope != 0 else float('inf')

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # Scatter
    ax.scatter(times, dists, c=NATURE_COLORS[0], s=80, alpha=0.85,
              edgecolors='white', linewidth=1.5, zorder=3)

    # Regression line
    x_range = np.array([times.min() - 0.5, times.max() + 0.5])
    y_pred = slope * x_range + intercept
    ax.plot(x_range, y_pred, '--', color=NATURE_COLORS[3], linewidth=2.5,
           alpha=0.8, zorder=2, label=f'R² = {r2:.3f}')

    # TMRCA vertical line
    if tmrca > times.min() and tmrca < times.max() + 5:
        ax.axvline(tmrca, color='#888888', linestyle=':', linewidth=1.5,
                  alpha=0.7, label=f'TMRCA ≈ {tmrca:.1f}')

    # Labels
    ax.set_xlabel('Sampling date (year)', fontsize=14, fontweight='medium')
    ax.set_ylabel('Root-to-tip distance (subs/site)', fontsize=14, fontweight='medium')
    ax.set_title(title, fontsize=16, fontweight='bold', pad=15)

    # Stats box
    stats_text = (
        f"Rate: {slope:.2e} subs/site/yr\n"
        f"R²: {r2:.3f}\n"
        f"TMRCA: {tmrca:.1f}\n"
        f"n = {len(valid_tips)}"
    )
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top',
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                     edgecolor='#cccccc', alpha=0.9))

    ax.legend(fontsize=11, loc='lower right', framealpha=0.9)
    ax.tick_params(labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, format='pdf', facecolor='white',
               edgecolor='none', bbox_inches='tight')
    plt.close(fig)

    log.emit(f"TempEst: {output_path} (R²={r2:.3f}, rate={slope:.2e}, n={len(valid_tips)})")
    return {"R2": r2, "slope": slope, "tmrca": tmrca, "p_value": p, "n": len(valid_tips)}


# ═══════════════════════════════════════════════════════════════════
# Fig 2: MCC Time-Tree with Location Colors
# ═══════════════════════════════════════════════════════════════════

def plot_mcc_tree_publication(
    mcc_tree_file: str,
    output_path: str,
    locations: Optional[List[str]] = None,
    tip_location_map: Optional[Dict[str, str]] = None,
    title: str = "Maximum Clade Credibility Tree",
    figsize: Tuple[float, float] = (12, 16),
    show_posterior: bool = True,
    log: Optional[LogCollector] = None,
):
    """
    发表级 MCC 树图：按地点着色 + 节点后验概率气泡。

    使用矩形分支布局（cladogram style），比 Phylo.draw 更适合 120+ taxa。
    """
    from Bio import Phylo
    log = log or LogCollector()

    try:
        tree = Phylo.read(mcc_tree_file, "newick")
    except Exception as e:
        log.warning(f"Cannot read MCC tree: {e}")
        return

    # Extract tip locations from tree annotations
    if tip_location_map is None:
        tip_location_map = {}
        for clade in tree.find_clades():
            if clade.is_terminal() and hasattr(clade, 'comment') and clade.comment:
                comment = clade.comment.strip('[]')
                for part in comment.split(','):
                    if 'location=' in part:
                        loc = part.split('location=')[-1].strip().strip('"')
                        tip_location_map[clade.name] = loc
                        break

    if not locations:
        locations = sorted(set(tip_location_map.values())) if tip_location_map else []

    # Color map
    color_map = {}
    for i, loc in enumerate(locations):
        color_map[loc] = NATURE_COLORS[i % len(NATURE_COLORS)]

    # Draw using rectangular layout
    fig, ax = plt.subplots(figsize=figsize)

    # 历史坑: Phylo.draw 参数是 show_confidence, show_conviction 是无效参数被静默吞掉
    # (节点后验概率永不显示)
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: c.name[:30] if c.name and c.is_terminal() else "",
               show_confidence=show_posterior,  # node labels as posterior probability
               branch_labels=lambda c: "")

    # Color tips
    for txt in ax.texts:
        name = txt.get_text().strip()
        if name in tip_location_map:
            loc = tip_location_map[name]
            txt.set_color(color_map.get(loc, '#333333'))
            txt.set_fontsize(5)
            txt.set_fontfamily('monospace')

    # Legend
    legend_elements = [Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=color_map[loc],
                              markersize=12, label=loc)
                       for loc in locations if loc in color_map]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9,
              title='Location', title_fontsize=10, framealpha=0.9)

    ax.set_title(title, fontsize=16, fontweight='bold', loc='left', pad=20)

    fig.savefig(output_path, dpi=300, format='pdf', facecolor='white',
               edgecolor='none', bbox_inches='tight')
    plt.close(fig)
    log.emit(f"MCC tree: {output_path} ({len(locations)} locations)")


# ═══════════════════════════════════════════════════════════════════
# Table S1: Parameter Estimates
# ═══════════════════════════════════════════════════════════════════

def build_parameter_table(
    log_file: str,
    output_path: str,
    burnin_pct: float = 10.0,
    ess_threshold: float = 200.0,
    log: Optional[LogCollector] = None,
) -> pd.DataFrame:
    """
    生成发表级参数估计表 (CSV + 可转 LaTeX)。

    Returns DataFrame with columns: Parameter, Mean, Median, 95% HPD Lower,
    95% HPD Upper, ESS, Converged.
    """
    log = log or LogCollector()

    if not os.path.exists(log_file):
        log.warning(f"Log not found: {log_file}")
        return pd.DataFrame()

    df = pd.read_csv(log_file, sep="\t", comment="#")
    burnin_idx = int(len(df) * burnin_pct / 100.0)
    df_burnin = df.iloc[burnin_idx:]

    skip_cols = {"Sample", "sample", "state", "State"}
    rows = []

    for col in df_burnin.select_dtypes(include=[np.number]).columns:
        if col in skip_cols:
            continue
        vals = df_burnin[col].dropna().values
        if len(vals) < 10:
            continue

        # ESS (Tracer 同算法, Geyer monotone sequence)
        ess = calculate_ess(vals)
        hpd_lo, hpd_hi = calculate_95hpd(vals)

        rows.append({
            "Parameter": col,
            "Mean": f"{np.mean(vals):.4e}",
            "Median": f"{np.median(vals):.4e}",
            "HPD_95_Lower": f"{hpd_lo:.4e}",
            "HPD_95_Upper": f"{hpd_hi:.4e}",
            "ESS": f"{ess:.1f}",
            "Converged": "Yes" if ess >= ess_threshold else "No",
        })

    table = pd.DataFrame(rows).sort_values("Parameter")

    # Save CSV
    table.to_csv(output_path, index=False)
    log.emit(f"Parameter table: {output_path} ({len(table)} params)")

    # Print key parameters
    key_patterns = ["clock", "ucld", "kappa", "popsize", "rootHeight",
                    "coefficientOfVariation", "posterior", "prior", "likelihood"]
    log.emit("\n  Key Parameters:")
    for _, row in table.iterrows():
        for kp in key_patterns:
            if kp.lower() in row["Parameter"].lower():
                flag = "✅" if row["Converged"] == "Yes" else "⚠️"
                log.emit(f"  {flag} {row['Parameter']:<35s} = {row['Mean']}  "
                        f"[{row['HPD_95_Lower']} – {row['HPD_95_Upper']}]  "
                        f"ESS={row['ESS']}")
                break

    return table


def _calc_ess(x: np.ndarray) -> float:
    """Deprecated 别名: 委托 utils.ess.calculate_ess (Tracer 同算法)。
    旧版逐 lag 正截断 AR(1) 近似会低估 ESS, 已废弃。"""
    return calculate_ess(x)


# ═══════════════════════════════════════════════════════════════════
# All-in-one: Generate all publication outputs
# ═══════════════════════════════════════════════════════════════════

def build_all_publication_outputs(
    rtt_fasta: str,
    rtt_dates: str,
    beast_log: str,
    output_dir: str,
    rtt_tree: Optional[str] = None,
    mcc_tree: Optional[str] = None,
    locations: Optional[List[str]] = None,
    prefix: str = "Fig",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    一次性生成所有发表级输出。

    Returns dict with all output paths.
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    outputs = {}

    # Fig 1: TempEst
    log.emit("=" * 50)
    log.emit("Fig 1: TempEst Root-to-Tip Regression")
    tempest_path = os.path.join(output_dir, f"{prefix}1_temprest.pdf")
    rtt_result = plot_temprest(
        fasta_file=rtt_fasta, dates_csv=rtt_dates,
        output_path=tempest_path, tree_file=rtt_tree,
        title="Root-to-Tip Regression of GCVA Sequences",
        log=log)
    outputs["tempest"] = tempest_path
    outputs["tempest_stats"] = rtt_result

    # Fig 2: MCC Tree
    if mcc_tree and os.path.exists(mcc_tree):
        log.emit("\nFig 2: MCC Time-Tree")
        mcc_path = os.path.join(output_dir, f"{prefix}2_mcc_tree.pdf")
        plot_mcc_tree_publication(
            mcc_tree_file=mcc_tree, output_path=mcc_path,
            locations=locations,
            title="Time-Calibrated Phylogeny of GCVA",
            log=log)
        outputs["mcc_tree"] = mcc_path

    # Table S1: Parameter Estimates
    if beast_log and os.path.exists(beast_log):
        log.emit("\nTable S1: Parameter Estimates")
        table_path = os.path.join(output_dir, f"{prefix}S1_parameters.csv")
        table = build_parameter_table(
            log_file=beast_log, output_path=table_path, log=log)
        outputs["param_table"] = table_path

        # Also generate LaTeX version
        if len(table) > 0:
            latex_path = os.path.join(output_dir, f"{prefix}S1_parameters.tex")
            _export_latex_table(table, latex_path)
            outputs["param_table_latex"] = latex_path

    log.emit(f"\nAll outputs in: {output_dir}")
    return outputs


def _export_latex_table(df, path):
    """Export parameter table as LaTeX."""
    # Map parameter names to LaTeX-friendly versions
    name_map = {
        "clock.rate": r"$\mu$ (clock rate)",
        "ucld.mean": r"$\mu_{UCLD}$ (mean rate)",
        "ucld.stdev": r"$\sigma_{UCLD}$ (rate stdev)",
        "kappa": r"$\kappa$ (Ts/Tv ratio)",
        "constant.popSize": r"$\theta$ (pop. size)",
        "treeModel.rootHeight": r"$T_{MRCA}$ (root height)",
        "alpha": r"$\alpha$ (gamma shape)",
        "coefficientOfVariation": r"CoV (rate variation)",
        "posterior": "Posterior",
        "prior": "Prior",
        "likelihood": "Likelihood",
    }

    with open(path, 'w') as f:
        f.write(r"\begin{table}[ht]" + "\n")
        f.write(r"\centering" + "\n")
        f.write(r"\caption{Parameter estimates from BEAST UCLN analysis.}" + "\n")
        f.write(r"\label{tab:parameters}" + "\n")
        f.write(r"\begin{tabular}{lrrrrrc}" + "\n")
        f.write(r"\hline" + "\n")
        f.write(r"Parameter & Mean & 95\% HPD Lower & 95\% HPD Upper & ESS & Converged \\" + "\n")
        f.write(r"\hline" + "\n")

        for _, row in df.iterrows():
            param = row["Parameter"]
            latex_name = name_map.get(param, param.replace('_', r'\_'))
            conv = r"\checkmark" if row["Converged"] == "Yes" else ""
            f.write(f"{latex_name} & {row['Mean']} & {row['HPD_95_Lower']} & "
                   f"{row['HPD_95_Upper']} & {row['ESS']} & {conv} \\\\\n")

        f.write(r"\hline" + "\n")
        f.write(r"\end{tabular}" + "\n")
        f.write(r"\end{table}" + "\n")

    return path


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Publication figures generator")
    parser.add_argument("--rtt_fasta", help="FASTA alignment for TempEst")
    parser.add_argument("--rtt_dates", help="Dates CSV")
    parser.add_argument("--rtt_tree", default=None, help="Tree for TempEst")
    parser.add_argument("--mcc_tree", default=None, help="MCC tree")
    parser.add_argument("--mcc_locations", default=None, help="Comma-separated locations")
    parser.add_argument("--beast_log", default=None, help="BEAST log file")
    parser.add_argument("-o", "--output_dir", default="publication_figures")
    parser.add_argument("--prefix", default="Fig")
    args = parser.parse_args()

    lc = LogCollector()
    locations = args.mcc_locations.split(",") if args.mcc_locations else None

    build_all_publication_outputs(
        rtt_fasta=args.rtt_fasta,
        rtt_dates=args.rtt_dates,
        rtt_tree=args.rtt_tree,
        mcc_tree=args.mcc_tree,
        beast_log=args.beast_log,
        output_dir=args.output_dir,
        locations=locations,
        prefix=args.prefix,
        log=lc)
