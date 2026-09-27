#!/usr/bin/env python3
"""
phylogeo_figures.py — 发表级系统地理学配图
===========================================
生成 SCI 期刊标准的系统地理学多面板图:

  Panel A: MCC 时间树 (按地点着色 + 节点后验概率)
  Panel B: 迁移路线地图 (SpreaD3 风格)
  Panel C: Bayesian Skyline Plot (有效群体大小变化)
  Panel D: 时钟速率分布 (UCLN mean + stdev)
  Panel E: 迁移 BF 热力图

配色方案: Nature-style deep color palette, colorblind-safe
"""

import csv, json, os, sys, re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from matplotlib import cm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector
# 历史坑: calculate_95hpd 曾只在 plot_skyline/plot_clock_rate 局部 import,
# _draw_skyline_on_ax 三处裸用 → NameError → Panel C 必空。
from utils.ess import calculate_95hpd

# ═══════════════════════════════════════════════════════════════════
# Color palette (Nature-inspired, colorblind-safe)
# ═══════════════════════════════════════════════════════════════════

NATURE_PALETTE = [
    '#E64B35',  # Red
    '#4DBBD5',  # Cyan
    '#00A087',  # Green
    '#3C5488',  # Blue
    '#F39B7F',  # Salmon
    '#8491B4',  # Grey-purple
    '#91D1C2',  # Mint
    '#DC0000',  # Dark red
    '#7E6148',  # Brown
    '#B09C85',  # Tan
]

LOCATION_COLORS = {
    "Beijing": "#E64B35",
    "Ningxia": "#4DBBD5",
    "Neimenggu": "#00A087",
    "China": "#3C5488",
    "Japan": "#F39B7F",
    "Korea": "#8491B4",
    "Europe": "#91D1C2",
    "India": "#DC0000",
    "USA": "#7E6148",
    "Other": "#B09C85",
}


# ═══════════════════════════════════════════════════════════════════
# Panel A: MCC Time-Tree with Location Colors
# ═══════════════════════════════════════════════════════════════════

def plot_mcc_tree(
    mcc_tree_file: str,
    output_path: str,
    locations: Optional[List[str]] = None,
    figsize: Tuple[float, float] = (10, 12),
    log: Optional[LogCollector] = None,
):
    """
    Draw MCC tree with branch colors by inferred location.
    """
    from Bio import Phylo

    log = log or LogCollector()

    try:
        tree = Phylo.read(mcc_tree_file, "newick")
    except Exception as e:
        log.warning(f"Cannot read MCC tree: {e}")
        return

    # Extract location annotations
    tip_locations = {}
    internal_locations = {}
    tip_posteriors = {}

    for clade in tree.find_clades():
        if hasattr(clade, 'comment') and clade.comment:
            comment = clade.comment.strip('[]')
            loc = None
            post = None
            for part in comment.split(','):
                part = part.strip()
                if '=' in part:
                    k, v = part.split('=', 1)
                    k, v = k.strip(), v.strip().strip('"')
                    if k == 'location':
                        loc = v
                    elif k == 'posterior':
                        try:
                            post = float(v)
                        except ValueError:
                            pass
            if loc:
                if clade.is_terminal():
                    tip_locations[clade.name] = loc
                else:
                    internal_locations[id(clade)] = loc
            if post:
                tip_posteriors[clade.name if clade.is_terminal() else id(clade)] = post

    if not locations:
        locations = sorted(set(list(tip_locations.values()) +
                               list(internal_locations.values())))

    # Color map
    color_map = {}
    for i, loc in enumerate(locations):
        color_map[loc] = LOCATION_COLORS.get(loc, NATURE_PALETTE[i % len(NATURE_PALETTE)])

    # Draw
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis('off')

    # Use rectangular layout
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: c.name[:25] if c.name and c.is_terminal() else "",
               branch_labels=lambda c: "",
               show_confidence=False)

    # Color tips by location
    for txt in ax.texts:
        name = txt.get_text()
        if name in tip_locations:
            txt.set_color(color_map.get(tip_locations[name], '#333333'))
            txt.set_fontsize(7)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[loc],
               markersize=10, label=loc)
        for loc in locations if loc in color_map
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9,
              title='Location', title_fontsize=10, framealpha=0.9)

    ax.set_title('A. Maximum Clade Credibility Tree', fontsize=14, fontweight='bold',
                loc='left', pad=10)

    fig.savefig(output_path, bbox_inches='tight', dpi=300, format='pdf',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    log.emit(f"MCC tree: {output_path}")


# ═══════════════════════════════════════════════════════════════════
# Panel C: Bayesian Skyline Plot
# ═══════════════════════════════════════════════════════════════════

def plot_skyline(
    log_file: str,
    output_path: str,
    burnin_pct: float = 10.0,
    figsize: Tuple[float, float] = (8, 5),
    log: Optional[LogCollector] = None,
):
    """
    Draw Bayesian Skyline Plot from BEAST log.

    2026-09-16 (T5): 时间轴改由 `_skyline_series` 还原 (组跨度/距今年), 不再用
    `np.linspace(0, 1, n)` 的编造刻度; 没有 groupSize 列时退回等距并在图上标注。
    """
    log = log or LogCollector()

    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
    except Exception as e:
        log.warning(f"Cannot read log: {e}")
        return

    burnin_idx = int(len(df) * burnin_pct / 100.0)
    df = df.iloc[burnin_idx:]

    ser = _skyline_series(df)
    if ser is None:
        log.warning("No skyline population size columns found")
        return

    fig, ax = plt.subplots(figsize=figsize)

    if ser["edges"] is not None:
        xs = list(ser["edges"])
        ax.set_xlabel('Years before present (0 = youngest sample)', fontsize=12)
    else:
        xs = [0.0] + list(ser["mid"])
        ax.set_xlabel('Relative time (groups, equally spaced)', fontsize=12)

    xl, y_med, y_lo, y_hi = [], [], [], []
    for i in range(ser["n"]):
        xl += [xs[i], xs[i + 1]]
        y_med += [ser["med"][i]] * 2
        y_lo += [ser["lo"][i]] * 2
        y_hi += [ser["hi"][i]] * 2

    ax.fill_between(xl, y_lo, y_hi, alpha=0.25, color='#4DBBD5', step='post')
    ax.step(xl, y_med, where='post', color='#3C5488', linewidth=2.5)

    ax.set_ylabel('Effective population size (Nₑτ)', fontsize=12)
    ax.set_title('C. Bayesian Skyline Plot', fontsize=14, fontweight='bold', loc='left')
    ax.tick_params(labelsize=10)
    if ser["warnings"]:
        ax.text(0.02, 0.02, "⚠ " + "; ".join(ser["warnings"])[:160],
                transform=ax.transAxes, fontsize=7, color='#b9770e', va='bottom')

    fig.savefig(output_path, bbox_inches='tight', dpi=300, format='pdf',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    log.emit(f"Skyline: {output_path}")


# ═══════════════════════════════════════════════════════════════════
# Panel C-2: Bayesian Skyline — 自包含交互 HTML (日历年 + 95% HPD 带)
# T5 (2026-09-16)。平台侧已有同型图 (tools.html skyline Ne), 管线侧缺。
# ═══════════════════════════════════════════════════════════════════

def _skyline_series(df) -> Optional[Dict]:
    """BEAST log (已去 burnin) → skyline 序列 + **真实**时间轴。

    历史坑: 旧 `plot_skyline` / `_draw_skyline_on_ax` 用 `np.linspace(0, 1, n)`
    当 x 轴 —— 刻度是编的, 读者无从判断扩张/收缩发生在哪一年, 也无法与采样年对齐。
    这里:
      · 组跨度取 `skyline.groupSize1..k` (BEAST 里 group 1 = 最靠近现在的组);
      · 量纲自检: sum(groupSize) 应 ≈ treeModel.rootHeight (差 >5% 就在告警里说清);
      · 有采样年时换算成**日历年** (= max_tip_year − 距今年), 无则明说「相对时间」。
    返回 None 表示 log 里没有 skyline 参数列 (调用方负责把原因写出来)。
    """
    pop_cols = [c for c in df.columns
                if ('skyline' in c.lower() and 'pop' in c.lower())]
    if not pop_cols:
        return None

    def _ord(col):
        m = re.search(r'(\d+)\s*$', str(col))
        return int(m.group(1)) if m else 0

    pop_cols = sorted(pop_cols, key=_ord)
    grp_cols = sorted([c for c in df.columns
                       if 'skyline' in c.lower() and 'groupsize' in c.lower()], key=_ord)

    hp = [calculate_95hpd(df[c].values) for c in pop_cols]
    med = [float(np.median(df[c].values)) for c in pop_cols]
    lo = [float(h[0]) for h in hp]
    hi = [float(h[1]) for h in hp]

    warn = []
    root_h = None
    for c in ('treeModel.rootHeight', 'tree.height', 'TreeHeight'):
        if c in df.columns:
            root_h = float(np.median(df[c].values))
            break

    bounds = None
    if len(grp_cols) == len(pop_cols):
        # BEAST 的 group 1 是**最年轻**组 → 距今年边界从 0 起累加
        spans = [float(np.median(df[c].values)) for c in grp_cols]
        tot = sum(spans)
        if root_h and root_h > 0:
            rel = abs(tot - root_h) / root_h
            if rel > 0.05:
                warn.append(f"组跨度之和 {tot:.3g} 与树高 {root_h:.3g} 相差 "
                            f"{rel*100:.1f}% (量纲可能不是年) → 时间轴按组跨度, 标为相对时间")
                bounds = None
            else:
                bounds = np.cumsum([0.0] + spans).tolist()
        else:
            warn.append("log 无 treeModel.rootHeight, 无法校验组跨度量纲 → 用组跨度当时间轴")
            bounds = np.cumsum([0.0] + spans).tolist()
    else:
        warn.append(f"log 无 skyline.groupSize 列 (找到 {len(grp_cols)} 个, popSize "
                    f"{len(pop_cols)} 个) → x 轴只能等距排布, 不代表真实时间")

    n = len(pop_cols)
    if bounds is None:
        # 等距回退: 中点在 [0,1] 均匀分布 —— 显式标注为「等距 (非真实时间)」
        xs = [(i + 0.5) / n for i in range(n)]
        edges = None
    else:
        edges = bounds  # 长度 n+1, 递增 (距今年)
        xs = [(bounds[i] + bounds[i + 1]) / 2.0 for i in range(n)]

    return {"pop_cols": pop_cols, "med": med, "lo": lo, "hi": hi,
            "edges": edges, "mid": xs, "root_h": root_h,
            "warnings": warn, "n": n}


def max_tip_year_from_metadata(meta_csv: str, log=None) -> Optional[float]:
    """从 sample_metadata.csv 取最大采样年 (十进制年)。取不到返回 None (不编数)。"""
    try:
        from utils.decimal_year import to_decimal_year
    except ImportError:
        return None
    if not meta_csv or not os.path.exists(meta_csv):
        return None
    yrs = []
    try:
        with open(meta_csv, encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                d = (row.get('date') or '').strip()
                if not d:
                    continue
                try:
                    yrs.append(to_decimal_year(d, sample=row.get('name', '')))
                except ValueError:
                    continue
    except OSError as e:
        if log is not None:
            log.warning(f"  ⚠ 读取采样元数据失败, 无法换算日历年: {e}")
        return None
    return max(yrs) if yrs else None


def build_skyline_html(log_file: str, output_path: str,
                       max_tip_year: Optional[float] = None,
                       burnin_pct: float = 10.0, height: int = 560,
                       log: Optional[LogCollector] = None) -> Optional[str]:
    """Bayesian skyline (中位数 + 95% HPD 带) → 单文件自包含 HTML, 返回路径或 None。

    max_tip_year: 最大采样年 (十进制年)。给了 → x 轴是**日历年**; 没给 → 距今年,
                  并在图上明说。x 轴意义永远写清楚, 不靠读者猜。
    """
    log = log or LogCollector()
    if not log_file or not os.path.exists(log_file):
        log.warning(f"  skyline HTML 跳过: log 不存在 ({log_file})")
        return None
    try:
        import plotly.graph_objects as go
        from utils.self_contained_html import write_figure
    except ImportError as e:
        log.warning(f"  skyline HTML 跳过: 缺 plotly / self_contained_html ({e})；"
                    f"修复: pip install plotly")
        return None

    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
    except Exception as e:
        log.warning(f"  skyline HTML 跳过: log 读取失败 ({e})")
        return None
    if len(df) < 10:
        log.warning(f"  skyline HTML 跳过: log 只有 {len(df)} 行 (无法去 burnin)")
        return None
    df = df.iloc[int(len(df) * burnin_pct / 100.0):]

    ser = _skyline_series(df)
    if ser is None:
        log.warning("  skyline HTML 跳过: log 无 skyline.popSize 列 → 该 run 的树先验不是 "
                    "skyline (检查 --tree_prior / beast XML 的 generalizedSkyLineLikelihood)")
        return None

    # ── x 轴: 日历年优先 ──
    edges, mid = ser["edges"], ser["mid"]
    if edges is not None:
        bp = list(edges)                                  # 距今年 (0=现在)
        if max_tip_year is not None:
            xs = [max_tip_year - b for b in bp]
            xlab = "采样年（日历年）"
        else:
            xs = bp
            xlab = "距今年代（年，0 = 采样最晚样本）"
            ser["warnings"].append("未提供 max_tip_year → x 轴为距今年, 不是日历年")
    else:
        xs = [0.0] + list(mid)                            # 等距回退
        xlab = "相对时间（等距，非真实时间）"

    # 逐步 (piecewise-constant) 折线: 每个组一条水平段
    xl, y_med, y_lo, y_hi = [], [], [], []
    for i in range(ser["n"]):
        xl += [xs[i], xs[i + 1] if i + 1 < len(xs) else xs[i]]
        y_med += [ser["med"][i]] * 2
        y_lo += [ser["lo"][i]] * 2
        y_hi += [ser["hi"][i]] * 2

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=xl, y=y_hi, mode="lines", line=dict(width=0),
                             hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(x=xl, y=y_lo, mode="lines", line=dict(width=0),
                             fill="tonexty", fillcolor="rgba(77,187,213,0.35)",
                             name="95% HPD", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=xl, y=y_med, mode="lines", line=dict(color="#3C5488", width=2.5),
                             name="中位数 (Nₑτ)"))

    nsig = f"{ser['n']} 个时间组"
    unit = "日历年" if max_tip_year is not None and edges is not None else "相对/距今年代"
    title = f"Bayesian Skyline — {Path(log_file).name}"
    fig.update_layout(
        title=dict(text=title, x=0.02, xanchor="left", font=dict(size=15)),
        xaxis=dict(title=dict(text=xlab), tickfont=dict(size=10), showgrid=True,
                   gridcolor="#eee"),
        yaxis=dict(title=dict(text="有效群体大小 Nₑτ"), tickfont=dict(size=10),
                   showgrid=True, gridcolor="#eee"),
        legend=dict(orientation="h", y=1.02, x=0.98, xanchor="right"),
        margin=dict(l=70, r=30, t=70, b=70), height=height,
        plot_bgcolor="white", font=dict(family="Microsoft YaHei, Segoe UI, Arial", size=11))

    _root_txt = f"{ser['root_h']:.3g}" if ser["root_h"] else "log 未提供 (无法校验量纲)"
    note = (f"burnin {burnin_pct:g}% 后 {len(df)} 条样本 · {nsig} · 时间轴: {unit} · "
            f"树高中位数 {_root_txt} · 95% HPD = Tracer 口径 (utils/ess.calculate_95hpd)")
    if ser["warnings"]:
        note += "<br>⚠ " + "；".join(ser["warnings"])
    if max_tip_year is not None and edges is not None:
        note += f"<br>日历年换算: 采样最晚样本 {max_tip_year:.2f} 年, 根 ~{max_tip_year - edges[-1]:.1f} 年"

    try:
        write_figure(fig, output_path, title=title, note=note, height=height, log=log)
    except RuntimeError as e:
        log.warning(f"  skyline HTML 生成失败: {e}")
        return None
    return output_path


# ═══════════════════════════════════════════════════════════════════
# Panel D: Clock Rate Distribution
# ═══════════════════════════════════════════════════════════════════

def plot_clock_rate(
    log_file: str,
    output_path: str,
    burnin_pct: float = 10.0,
    figsize: Tuple[float, float] = (8, 5),
    log: Optional[LogCollector] = None,
):
    """Draw clock rate posterior density + UCLD stats."""
    log = log or LogCollector()

    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
    except Exception:
        return

    burnin_idx = int(len(df) * burnin_pct / 100.0)
    df = df.iloc[burnin_idx:]

    # Find clock rate column
    # 历史坑: 只认 clock.rate/clockRate/meanRate, BEAST1 UCLN 的列是 ucld.mean → 永不生成
    clock_cols = [c for c in df.columns
                  if c in ('clock.rate', 'clockRate', 'meanRate', 'ucld.mean', 'default.meanRate')
                  or 'clockRate' in c]
    if not clock_cols:
        # Fallback: search
        for c in df.columns:
            if ('clock' in c.lower() and 'rate' in c.lower() and 'ucld' not in c.lower()) \
                    or c == 'ucld.mean':
                clock_cols = [c]
                break

    if not clock_cols:
        return

    clock_vals = df[clock_cols[0]].values

    fig, ax = plt.subplots(figsize=figsize)

    ax.hist(clock_vals, bins=60, color='#E64B35', alpha=0.7, density=True,
            edgecolor='white', linewidth=0.5)
    ax.axvline(np.mean(clock_vals), color='#3C5488', linewidth=2.5,
              linestyle='--', label=f'Mean = {np.mean(clock_vals):.2e}')
    from utils.ess import calculate_95hpd
    _cl, _ch = calculate_95hpd(clock_vals)
    ax.axvline(_cl, color='#888888', linewidth=1.5,
              linestyle=':', alpha=0.8)
    ax.axvline(_ch, color='#888888', linewidth=1.5,
              linestyle=':', alpha=0.8, label=f'95% HPD')

    ax.set_xlabel('Clock rate (substitutions/site/year)', fontsize=12)
    ax.set_ylabel('Posterior density', fontsize=12)
    ax.set_title('D. Molecular Clock Rate', fontsize=14, fontweight='bold', loc='left')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.tick_params(labelsize=10)

    fig.savefig(output_path, bbox_inches='tight', dpi=300, format='pdf',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    log.emit(f"Clock rate: {output_path}")


# ═══════════════════════════════════════════════════════════════════
# Panel E: Migration BF Heatmap
# ═══════════════════════════════════════════════════════════════════

def plot_migration_heatmap(
    migration_csv: str,
    output_path: str,
    locations: List[str],
    figsize: Tuple[float, float] = (6, 5),
    log: Optional[LogCollector] = None,
):
    """Draw migration Bayes Factor heatmap."""
    log = log or LogCollector()

    if not os.path.exists(migration_csv):
        log.warning(f"Migration CSV not found: {migration_csv}")
        return

    df = pd.read_csv(migration_csv)

    n_loc = len(locations)
    matrix = np.zeros((n_loc, n_loc))

    for _, row in df.iterrows():
        from_loc = row.get("from", row.get("From", ""))
        to_loc = row.get("to", row.get("To", ""))
        bf = float(row.get("bf", row.get("BF", 0)))

        if from_loc in locations and to_loc in locations:
            fi = locations.index(from_loc)
            ti = locations.index(to_loc)
            matrix[fi, ti] = np.log10(max(bf, 1.0))

    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto', vmin=0,
                   vmax=max(np.log10(10), matrix.max()))

    # Labels
    ax.set_xticks(range(n_loc))
    ax.set_yticks(range(n_loc))
    ax.set_xticklabels(locations, rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(locations, fontsize=10)
    ax.set_title('E. Migration Routes (log₁₀ BF)', fontsize=14,
                fontweight='bold', loc='left')

    # Annotate significant routes
    for fi in range(n_loc):
        for ti in range(n_loc):
            if fi != ti and matrix[fi, ti] > 0:
                bf_val = 10 ** matrix[fi, ti]
                text = f'{bf_val:.0f}' if bf_val < 100 else f'{bf_val:.0f}'
                color = 'white' if matrix[fi, ti] > np.log10(3) else 'black'
                ax.text(ti, fi, f'BF={bf_val:.0f}', ha='center', va='center',
                       fontsize=8, color=color, fontweight='bold')

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('log₁₀(Bayes Factor)', fontsize=10)

    fig.savefig(output_path, bbox_inches='tight', dpi=300, format='pdf',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    log.emit(f"Migration heatmap: {output_path}")


# ═══════════════════════════════════════════════════════════════════
# Multi-Panel Figure (A+B+C+D+E)
# ═══════════════════════════════════════════════════════════════════

def build_multipanel_figure(
    log_file: str,
    trees_file: str,
    migration_csv: str,
    locations: List[str],
    output_path: str,
    mcc_tree_file: Optional[str] = None,
    burnin_pct: float = 10.0,
    title: str = "Bayesian Phylogeography of Virus X",
    log: Optional[LogCollector] = None,
):
    """
    Generate a 5-panel publication figure:
      A: MCC tree with location colors
      B: Migration map (simplified)  
      C: Bayesian Skyline Plot
      D: Clock rate distribution
      E: Migration BF heatmap
    """
    log = log or LogCollector()
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

    fig = plt.figure(figsize=(18, 22))

    # Layout: A full-width top, B+C+D+E in 2×2 grid below
    gs = fig.add_gridspec(3, 2, height_ratios=[1.5, 1, 1],
                          hspace=0.35, wspace=0.3)

    # ── Panel A: MCC Tree (spans top row, both columns) ──
    ax_a = fig.add_subplot(gs[0, :])
    try:
        from Bio import Phylo
        tree = Phylo.read(trees_file, "newick")
        Phylo.draw(tree, axes=ax_a, do_show=False,
                   label_func=lambda c: c.name[:20] if c.name and c.is_terminal() else "",
                   branch_labels=lambda c: "")
        ax_a.set_title('A. Maximum Clade Credibility Tree', fontsize=14,
                      fontweight='bold', loc='left')
    except Exception as e:
        ax_a.text(0.5, 0.5, f'Tree not available\n{e}', ha='center', va='center',
                 transform=ax_a.transAxes, fontsize=12, color='gray')
        ax_a.set_title('A. Maximum Clade Credibility Tree', fontsize=14,
                      fontweight='bold', loc='left')

    # ── Panel B: Migration Network (circular layout) ──
    # G5 (2026-09-16): 原名 'B. Migration Routes (BF≥5)' 会被读成地图, 实际是
    # 「地点排成圆环 + 有向箭头」的网络图 (无经纬度、无底图, 见 _draw_simple_migration_map)。
    # 标题与脚注都写明是圆周布局; 阈值不写死 5 (实际阈值由 --phylogeo_bf 决定)。
    ax_b = fig.add_subplot(gs[1, 0])
    _draw_simple_migration_map(ax_b, migration_csv, locations)
    ax_b.set_title('B. Migration Network (circular layout)', fontsize=14,
                  fontweight='bold', loc='left')
    ax_b.text(0.5, -0.04, 'nodes = locations on a circle (no geographic projection) · '
                          'arrow width ∝ BF · significant routes only',
              transform=ax_b.transAxes, ha='center', va='top',
              fontsize=8, color='#7f8c8d')

    # ── Panel C: Skyline Plot ──
    ax_c = fig.add_subplot(gs[1, 1])
    _draw_skyline_on_ax(ax_c, log_file, burnin_pct)
    ax_c.set_title('C. Bayesian Skyline Plot', fontsize=14,
                  fontweight='bold', loc='left')

    # ── Panel D: Clock Rate ──
    ax_d = fig.add_subplot(gs[2, 0])
    _draw_clock_on_ax(ax_d, log_file, burnin_pct)
    ax_d.set_title('D. Clock Rate Distribution', fontsize=14,
                  fontweight='bold', loc='left')

    # ── Panel E: Migration Heatmap ──
    ax_e = fig.add_subplot(gs[2, 1])
    _draw_heatmap_on_ax(ax_e, migration_csv, locations)
    ax_e.set_title('E. Migration BF Heatmap', fontsize=14,
                   fontweight='bold', loc='left')

    # Suptitle
    fig.suptitle(title, fontsize=18, fontweight='bold', y=0.98)

    fig.savefig(output_path, bbox_inches='tight', dpi=300, format='pdf',
                facecolor='white', edgecolor='none')
    plt.close(fig)

    log.emit(f"Multi-panel figure: {output_path}")
    return output_path


def _draw_simple_migration_map(ax, migration_csv, locations):
    """Draw simplified migration map with arrows."""
    if not os.path.exists(migration_csv):
        ax.text(0.5, 0.5, 'Migration data not available', ha='center',
               va='center', transform=ax.transAxes)
        return

    # Simple layout: arrange locations in a circle
    n = len(locations)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    positions = {loc: (np.cos(a), np.sin(a)) for loc, a in zip(locations, angles)}

    df = pd.read_csv(migration_csv)
    # 历史坑: 只认 'Yes' 过滤, 但 extract_migration_bf 输出的 significant 列是
    # True/False (bool to_csv) → 恒 False → Panel B 迁移箭头全丢。兼容多种表示。
    sig_mask = None
    for col in ('significant', 'Significant'):
        if col in df.columns:
            sig_mask = df[col].astype(str).str.upper().isin(['TRUE', 'YES', '1', 'T'])
            break
    sig_routes = df[sig_mask] if sig_mask is not None and sig_mask.any() else df

    # Draw nodes
    for i, loc in enumerate(locations):
        x, y = positions[loc]
        color = LOCATION_COLORS.get(loc, NATURE_PALETTE[i % len(NATURE_PALETTE)])
        ax.scatter(x, y, s=300, c=color, edgecolors='white', linewidth=2,
                  zorder=3)
        ax.annotate(loc, (x, y), textcoords="offset points", xytext=(0, 15),
                   ha='center', fontsize=10, fontweight='bold')

    # Draw routes
    max_bf = 1
    for _, row in sig_routes.iterrows():
        bf = float(row.get('BF', row.get('bf', 0)))
        max_bf = max(max_bf, bf)

    for _, row in sig_routes.iterrows():
        fr = row.get('From', row.get('from', ''))
        to = row.get('To', row.get('to', ''))
        bf = float(row.get('BF', row.get('bf', 0)))

        if fr in positions and to in positions:
            x1, y1 = positions[fr]
            x2, y2 = positions[to]
            width = 1 + 5 * bf / max(max_bf, 1)
            alpha = 0.3 + 0.5 * bf / max(max_bf, 1)
            ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                       arrowprops=dict(arrowstyle='->', color='#E64B35',
                                      lw=width, alpha=alpha,
                                      connectionstyle='arc3,rad=0.2'))

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.axis('off')


def _draw_skyline_on_ax(ax, log_file, burnin_pct):
    """Draw skyline on existing axis.

    2026-09-16 (T5): x 轴不再用 `np.linspace(0, 1, n)` —— 那是编出来的刻度。
    改用 `_skyline_series` 的组跨度还原真实时间 (有则距今年, 无则等距并明说)。
    """
    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
        burnin_idx = int(len(df) * burnin_pct / 100.0)
        df = df.iloc[burnin_idx:]

        ser = _skyline_series(df)
        if ser is None:
            # Single constant pop size
            for c in ['constant.popSize', 'popSize.t:aln']:
                if c in df.columns:
                    vals = df[c].values
                    ax.axhline(np.median(vals), color='#3C5488', linewidth=2, alpha=0.5)
                    ax.fill_between([0, 1],
                                   [calculate_95hpd(vals)[0]] * 2,
                                   [calculate_95hpd(vals)[1]] * 2,
                                   alpha=0.2, color='#4DBBD5')
                    ax.set_ylabel('Nₑτ (constant)', fontsize=11)
                    ax.set_xlabel('Time (arbitrary)', fontsize=11)
                    return

        if ser["edges"] is not None:
            xs = list(ser["edges"])
            ax.set_xlabel('Years before present (0 = youngest sample)', fontsize=11)
        else:
            xs = [0.0] + list(ser["mid"])
            ax.set_xlabel('Relative time (groups, equally spaced)', fontsize=11)

        xl, y_med, y_lo, y_hi = [], [], [], []
        for i in range(ser["n"]):
            xl += [xs[i], xs[i + 1]]
            y_med += [ser["med"][i]] * 2
            y_lo += [ser["lo"][i]] * 2
            y_hi += [ser["hi"][i]] * 2

        ax.fill_between(xl, y_lo, y_hi, alpha=0.25, color='#4DBBD5', step='post')
        ax.plot(xl, y_med, color='#3C5488', linewidth=2.5, drawstyle='steps-post')
        if ser["warnings"]:
            ax.text(0.02, 0.02, "⚠ " + "; ".join(ser["warnings"])[:160],
                    transform=ax.transAxes, fontsize=7, color='#b9770e', va='bottom')
        ax.set_ylabel('Effective population size', fontsize=11)
    except Exception:
        ax.text(0.5, 0.5, 'Skyline data unavailable', ha='center', va='center',
               transform=ax.transAxes)


def _draw_clock_on_ax(ax, log_file, burnin_pct):
    """Draw clock rate on existing axis."""
    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
        burnin_idx = int(len(df) * burnin_pct / 100.0)
        df = df.iloc[burnin_idx:]

        for c in ['clock.rate', 'clockRate', 'clockRate.c:aln', 'ucld.mean']:
            if c in df.columns:
                vals = df[c].values
                ax.hist(vals, bins=50, color='#E64B35', alpha=0.7, density=True,
                       edgecolor='white', linewidth=0.5)
                ax.axvline(np.mean(vals), color='#3C5488', linewidth=2, linestyle='--')
                ax.set_xlabel('Substitution rate (subs/site/yr)', fontsize=11)
                ax.set_ylabel('Density', fontsize=11)
                return
    except Exception:
        pass
    ax.text(0.5, 0.5, 'Clock rate unavailable', ha='center', va='center',
           transform=ax.transAxes)


def _draw_heatmap_on_ax(ax, migration_csv, locations):
    """Draw migration heatmap on existing axis."""
    try:
        df = pd.read_csv(migration_csv)
        n = len(locations)
        matrix = np.zeros((n, n))

        for _, row in df.iterrows():
            fr = str(row.get('From', row.get('from', '')))
            to = str(row.get('To', row.get('to', '')))
            bf = float(row.get('BF', row.get('bf', 0)))
            if fr in locations and to in locations:
                matrix[locations.index(fr), locations.index(to)] = np.log10(max(bf, 1))

        im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(locations, rotation=45, ha='right', fontsize=9)
        ax.set_yticklabels(locations, fontsize=9)
        plt.colorbar(im, ax=ax, shrink=0.8, label='log₁₀ BF')
    except Exception:
        ax.text(0.5, 0.5, 'Migration data unavailable', ha='center', va='center',
               transform=ax.transAxes)


def main():
    """CLI 入口: 发表级系统地理多面板图 (5 面板)"""
    import argparse
    ap = argparse.ArgumentParser(description="发表级系统地理配图")
    ap.add_argument("--log", required=True, help="BEAST log")
    ap.add_argument("--trees", required=True, help="BEAST trees")
    ap.add_argument("--migration", default="", help="迁移 BF CSV (可选)")
    ap.add_argument("--locations", default="", help="逗号分隔地点列表")
    ap.add_argument("--mcc-tree", default=None, help="MCC 树 (可选)")
    ap.add_argument("-o", "--output", required=True, help="输出 PDF")
    ap.add_argument("--title", default="Bayesian Phylogeography", help="图标题")
    args = ap.parse_args()
    locs = [l.strip() for l in args.locations.split(",") if l.strip()]
    build_multipanel_figure(
        log_file=args.log, trees_file=args.trees,
        migration_csv=args.migration or os.path.join(os.path.dirname(args.output), "migration_bf.csv"),
        locations=locs, output_path=args.output,
        mcc_tree_file=args.mcc_tree, title=args.title)
    print("PHYLOGEO_FIGURES_DONE")


if __name__ == "__main__":
    main()
