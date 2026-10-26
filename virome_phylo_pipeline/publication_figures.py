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
# 2026-10-16: 图件双格式 (PNG+PDF 同主名); 统一由 plot_io.save_both 收口
from utils.plot_io import save_both
from utils.ess import calculate_ess, calculate_95hpd

# ═══════════════════════════════════════════════════════════════════
# Color palette (Nature-style)
# ═══════════════════════════════════════════════════════════════════

NATURE_COLORS = ['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F',
                 '#8491B4','#91D1C2','#DC0000','#7E6148']


def read_tree_any(tree_file: str, log: Optional[LogCollector] = None):
    """读一棵树, **按内容自动识别** NEXUS / Newick; 多树取第一棵。读不出抛异常。

    2026-10-11 修 (P0): 旧实现硬编码 `Phylo.read(f, "newick")`, 而 BEAST/TreeAnnotator
    合并产出的 `time/beast/merged/mcc.tree` 是 **NEXUS**
    (`#NEXUS` … `tree TREE1 = [&R] (((…)));`)。按 newick 读 NEXUS 时 Biopython 会把
    注释里的分号当成树分隔符 → 抛
        "There are multiple trees in this file; use parse() instead."
    → `plot_mcc_tree_publication` 直接 return → **Fig 2 MCC 树静默缺失**
    (实测 PSTVd: report/publication/ 里只有 Fig1_temprest.pdf + FigS1_*)。

    实测该文件按 nexus 解析 = 1 棵树 / 129 叶; 按 newick 解析 = 8 棵 (全是残片)。
    """
    from Bio import Phylo
    log = log or LogCollector()
    try:
        with open(tree_file, "r", errors="replace") as fh:
            head = fh.read(4096).lstrip()
    except OSError as e:
        raise ValueError(f"树文件读不了: {tree_file} ({e})")
    first = "nexus" if head[:6].lower().startswith("#nexus") else "newick"
    order = (first, "newick" if first == "nexus" else "nexus")
    for fmt in order:
        try:
            trees = list(Phylo.parse(tree_file, fmt))
        except Exception:                       # noqa: BLE001 — 换格式再试
            continue
        if not trees:
            continue
        if fmt != first:
            log.emit(f"  [tree] {os.path.basename(tree_file)}: 按 {first} 解析失败, "
                     f"改按 {fmt} 成功")
        if len(trees) > 1:
            log.emit(f"  [tree] {os.path.basename(tree_file)} 含 {len(trees)} 棵树, 取第一棵")
        return trees[0]
    raise ValueError(f"无法解析树文件 (试过 {'/'.join(order)}): {tree_file}")


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
    # 2026-10-11 修 (P1): 原写法 `float(date.replace('-','.').split('.')[0])` 把
    # `2024-04` 一律截成 `2024.0` —— **月份/日期信息全丢**, 所有样本挤在整数年上,
    # root-to-tip 回归的分辨率被人为拉低 (而 clock_summary.tsv 用的是小数年,
    # 于是 Fig1 与 clock 阶段口径不一致)。改用全管线统一的 utils.decimal_year。
    from utils.decimal_year import to_decimal_year_safe
    dates = {}
    with open(dates_csv, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('name', '').strip()
            date = row.get('date', '').strip()
            if name and date:
                dy = to_decimal_year_safe(date, sample=name)
                if dy is not None:
                    dates[name] = dy

    # Build or load tree
    if tree_file and os.path.exists(tree_file):
        # 2026-10-11: 走格式自动识别 (IQ-TREE 是 newick, 但调用方可能传 NEXUS)
        tree = read_tree_any(tree_file, log)
    else:
        # Build NJ tree from alignment
        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

    # Compute root-to-tip distances
    # 2026-10-11 修 (P1): 原先无条件 `tree.root_with_outgroup(tree.get_terminals()[0])`
    # —— 拿**任意一个 tip** 当外群重定根。root-to-tip 距离完全取决于根的位置,
    # 这种定根没有生物学依据; 更严重的是它与管线 clock 阶段口径不一致
    # (clock_analysis.root_to_tip_distances 直接按 Newick 给定的根算, 从不重定根),
    # 于是同一批序列「Fig1 的 R²」和「clock_summary.tsv 的 R²」对不上。
    # 现改为**沿用输入树既有的根**, 与 clock 阶段严格同源。
    # (需要换根请在 tree 阶段产出已定根的树, 不在这里静默重定。)
    # 2026-09-16 修: `Bio.Phylo.depths()` 返回的 dict **以 Clade 对象为键**(不是名字) ——
    # 原写法用字符串 `rtt_dists[n]` 查, 命中率恒为 0 → `valid_tips` 恒空 →
    # 这张论文级 Fig 1 (TempEst root-to-tip 回归图) **一直静默不产出**
    # (只在日志里留一句 "Only 0 dated tips")。改为先建 name→depth 映射。
    # 注: 主树加了分支支持值后, 内部节点也有 name, 故只取 terminal。
    # 2026-10-12 修 (P1-4): 旧写法 `tree.depths(unit_branch_lengths=True)` ——
    # 该参数语义是"把所有分支长度当 1", 返回值是**从根走下来的边数**(计数),
    # 不是累加枝长; 而 Fig1 的 y 轴写的是 subs/site, 回归得到的 R²/TMRCA 全是
    # 错口径 (实测阶梯树: 错口径 R²=0.028 / slope=-0.077 "边/年", 正确口径
    # R²=0.966 / slope=1.9e-4 subs/site/yr)。改为累加枝长, 与
    # utils/clock_analysis.root_to_tip_distances 严格同源 (实测逐 tip 一致)。
    rtt_dists = tree.depths(unit_branch_lengths=False)
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
    ax.set_xlabel('Sampling date (year)', fontsize=14, fontweight='normal')     # P3-15: 'medium' 无对应字面 → normal
    ax.set_ylabel('Root-to-tip distance (subs/site)', fontsize=14, fontweight='normal')
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
    save_both(fig, output_path, dpi=300, format='pdf', facecolor='white',
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
        # 2026-10-11 修 (P0): 原来硬编码 "newick" —— 而 mcc.tree 是 NEXUS, 于是
        # 必然抛 "There are multiple trees in this file" → Fig 2 静默缺失。
        tree = read_tree_any(mcc_tree_file, log)
    except Exception as e:
        log.warning(f"Cannot read MCC tree: {e}")
        return

    # 2026-10-11 修 (P1): Biopython 的 Phylo.draw 会把 `tree.name` 写成**居中标题**
    # (Bio/Phylo/_utils.py:592 `axes.set_title(name)`), 而 NEXUS 解析出的
    # tree.name 就是树块名 "TREE1"。我们自己的标题用的是 loc='left' —— 二者是
    # **互相独立的 Text artist**, 于是同一行并排出现。实测 pdftotext 读出:
    #     "Time-Calibrated Phylogeny of GCVATREE1"
    # 清空 tree.name 即可: Biopython 只在 `if name:` 为真时才 set_title。
    tree.name = ""

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
    # 2026-10-11 修 (P1): 无 location 注释时 (实测 PSTVd `mcc.tree` 里 `location=`
    # 命中 **0** 次 —— BEAST 合并树只保留 taxon 名, 不带地点) 旧代码仍画一个
    # **空图例** (只剩 "Location" 标题、零条目), 且所有 tip 都不着色 —— 与项目
    # 「退化守卫」纪律相悖 (统计量算不出来就不该产出对应图元)。改为有地点才画。
    if locations:
        legend_elements = [Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor=color_map[loc],
                                  markersize=12, label=loc)
                           for loc in locations if loc in color_map]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=9,
                  title='Location', title_fontsize=10, framealpha=0.9)
    else:
        log.emit("  [mcc] 树中无 location 注释 → 不做地点着色 / 不画空图例")

    ax.set_title(title, fontsize=16, fontweight='bold', loc='left', pad=20)

    save_both(fig, output_path, dpi=300, format='pdf', facecolor='white',
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

    # 2026-10-11 修 (P0): rows 为空时 `pd.DataFrame([])` **没有任何列**,
    # 直接 `.sort_values("Parameter")` 抛 `KeyError: 'Parameter'`。异常冒泡到
    # `build_all_publication_outputs` → **整个 publication 阶段中止**,
    # 连后面的 Table S1 都产不出。实测 GCVA: `report/publication/` 是空目录,
    # 日志只有一行 `[postprocess] publication_figures 失败: KeyError: 'Parameter'`。
    # (空表是正常情形: beast_log 里没有可用的数值列 / 未传 --beast_log。)
    table = pd.DataFrame(rows)
    if len(table) > 0:
        table = table.sort_values("Parameter")

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
    virus_name: str = "",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    一次性生成所有发表级输出。

    Returns dict with all output paths.
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    outputs = {}

    # 2026-10-11 修 (P1): Fig1/Fig2 的标题原先把 "GCVA" **硬编码**在函数体里 ——
    # 于是跑 PSTVd 时两张发表级图的标题都写成了 GCVA。实测 (PSTVd 活副本,
    # 用 PSTVd 的 mcc.tree 直跑) pdftotext 读出:
    #     "Time-Calibrated Phylogeny of GCVATREE1"
    # 现按 virus_name 生成; 未提供时退化为**不带物种名**的通用标题
    # (宁可通用, 也不能张冠李戴)。
    _vlabel = (virus_name or "").strip()
    _fig1_title = (f"Root-to-Tip Regression of {_vlabel} Sequences"
                   if _vlabel else "Root-to-Tip Regression")
    _fig2_title = (f"Time-Calibrated Phylogeny of {_vlabel}"
                   if _vlabel else "Time-Calibrated Phylogeny")

    # ── 2026-10-11 修 (P1): "声明即存在" ──
    # 旧版无条件写 `outputs["tempest"] = tempest_path`, 即使 plot_temprest 中途
    # return (日期不足 3 条 / 读不了比对) —— 于是下游看到的是"有产物"但文件不存在,
    # 排查时只能靠猜。现改为: **文件真的落盘了才登记**, 否则记 None 并告警。
    def _declare(key, path, what):
        if path and os.path.exists(path) and os.path.getsize(path) > 0:
            outputs[key] = path
            return True
        outputs[key] = None
        log.warning(f"  [publication] {what} 未产出 ({path}) —— 不登记为交付物")
        return False

    # Fig 1: TempEst
    log.emit("=" * 50)
    log.emit("Fig 1: TempEst Root-to-Tip Regression")
    tempest_path = os.path.join(output_dir, f"{prefix}1_temprest.pdf")
    rtt_result = plot_temprest(
        fasta_file=rtt_fasta, dates_csv=rtt_dates,
        output_path=tempest_path, tree_file=rtt_tree,
        title=_fig1_title,
        log=log)
    _declare("tempest", tempest_path, "Fig 1 TempEst 回归图")
    outputs["tempest_stats"] = rtt_result

    # Fig 2: MCC Tree
    if mcc_tree and os.path.exists(mcc_tree):
        log.emit("\nFig 2: MCC Time-Tree")
        mcc_path = os.path.join(output_dir, f"{prefix}2_mcc_tree.pdf")
        plot_mcc_tree_publication(
            mcc_tree_file=mcc_tree, output_path=mcc_path,
            locations=locations,
            title=_fig2_title,
            log=log)
        _declare("mcc_tree", mcc_path, "Fig 2 MCC 时间树")
    else:
        # 2026-10-11: 原先**完全静默** —— 论文级三件套少一张图却毫无提示。
        # MCC 树来自 BEAST 多链合并产物 (time/beast/merged/mcc.tree), 多链模式
        # 下若未做合并, 该文件就不存在 (GCVA 实测: beast11/12/13 都跑完了,
        # 但 merged/ 从未生成 → 论文主图缺一张而日志无任何线索)。
        log.warning(f"  [publication] 跳过 Fig 2 (MCC 时间树): 未提供或不存在 "
                    f"--mcc_tree ({mcc_tree})。"
                    f"若 BEAST 多链已跑完, 请先合并链 (time/beast/merged/mcc.tree) 再重跑本步骤。")

    # Table S1: Parameter Estimates
    if beast_log and os.path.exists(beast_log):
        log.emit("\nTable S1: Parameter Estimates")
        table_path = os.path.join(output_dir, f"{prefix}S1_parameters.csv")
        table = build_parameter_table(
            log_file=beast_log, output_path=table_path, log=log)
        _declare("param_table", table_path, "Table S1 参数估计表")

        # Also generate LaTeX version
        if len(table) > 0:
            latex_path = os.path.join(output_dir, f"{prefix}S1_parameters.tex")
            _export_latex_table(table, latex_path)
            _declare("param_table_latex", latex_path, "Table S1 LaTeX 版")
    else:
        log.warning(f"  [publication] 跳过 Table S1 (参数估计表): 未提供或不存在 "
                    f"--beast_log ({beast_log})")

    # ── Fig 3: 管线各阶段代表性结果总览（2026-10-12 新增）──
    # 与 Fig1 (TempEst 时间信号) / Fig2 (MCC 时间树) 同套，补上"管线每跑一个阶段
    # 都产出了什么"。data_dir 由 output_dir 反推：output_dir = <data_dir>/report/publication
    # （legacy 布局 report/；standard 布局是 08_report/ —— 反推对两套都成立）。
    _data_root = os.path.dirname(os.path.dirname(os.path.abspath(output_dir)))
    if os.path.isdir(_data_root):
        log.emit("\nFig 3: Pipeline stage overview")
        try:
            from utils.pipeline_overview_fig import build_pipeline_overview
            _ov = build_pipeline_overview(
                _data_root, output_dir, title=_vlabel or "Virus",
                prefix=prefix, log=log)
            _declare("pipeline_overview_pdf", _ov.get("pdf"),
                     "Fig 3 管线各阶段总览 (PDF)")
            _declare("pipeline_overview_png", _ov.get("png"),
                     "Fig 3 管线各阶段总览 (PNG)")
        except Exception as _e:                       # 总览图失败不得阻断 Fig1/Fig2
            log.warning(f"  [publication] Fig 3 生成失败 (不阻断其它产物): {_e}")
    else:
        log.warning(f"  [publication] 跳过 Fig 3: 由 output_dir 反推的 data_dir 不存在 "
                    f"({_data_root})")

    # ── Fig 4: 核心进化结果（2026-10-13 新增）──
    # 与 Fig3「管线各阶段总览」互补：Fig3 回答"每个阶段产出了什么"，
    # Fig4 只收 5 张最有信息量的结果图，按叙事顺序排成一张可独立阅读的主图
    # （系统发育 → 采样结构 → 重组 → 正选择 → 地理传播）。
    # 与 Fig3 同样由 output_dir 反推 data_dir（legacy/standard 两套布局都成立）。
    if os.path.isdir(_data_root):
        log.emit("\nFig 4: Core evolutionary results")
        try:
            from utils.core_results_fig import build_core_results
            _cr = build_core_results(
                _data_root, output_dir, title=_vlabel or "Virus", log=log)
            _declare("core_results_pdf", _cr.get("pdf"),
                     "Fig 4 核心进化结果 (PDF)")
            _declare("core_results_png", _cr.get("png"),
                     "Fig 4 核心进化结果 (PNG)")
        except Exception as _e:                       # 核心图失败不得阻断其它产物
            log.warning(f"  [publication] Fig 4 生成失败 (不阻断其它产物): {_e}")
    else:
        log.warning(f"  [publication] 跳过 Fig 4: 由 output_dir 反推的 data_dir 不存在 "
                    f"({_data_root})")

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
    # 2026-10-11: 图题里的物种名此前硬编码 GCVA; 现由调用方 (phylo_pipeline) 传入。
    parser.add_argument("--virus_name", default="",
                        help="病毒名, 用于 Fig1/Fig2 标题 (缺省则用不带物种名的通用标题)")
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
        virus_name=args.virus_name,
        log=lc)
