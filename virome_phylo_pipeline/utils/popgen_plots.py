#!/usr/bin/env python3
"""
popgen_plots.py — 群体遗传统计多面板可视化 (π/θ/Hd/Tajima's D/Fu's Fs/Fst)
==========================================================================

2026-10-12 新增。此前 popgen stage 的统计量 (π/θ/Hd/Tajima's D/Fu's Fs/Fst)
只有 popgen_report.txt 纯文本, 没有任何图; pypopart 的 visualization 模块也
只画单倍型网络 (static_plot.create_publication_figure), 不画统计图。

数据源: popgen_analysis --stats-json 产出的结构化 JSON (全局 + 按群体 +
pairwise Fst 矩阵)。渲染 2×3 多面板 (matplotlib, 纯本地, 不依赖 pypopart):

  A  π vs θ-W   按群体分组柱状图   (Bizard: Composition/GroupedBarplot)
  B  Hd         柱状图 (0-1), 柱底标 n
  C  Tajima's D 发散条形图 (0 线基准, 负=群体扩张信号)
  D  Fu's Fs    发散条形图 (同上; Ewens 口径, θ̂ 由单倍型数拟合)
  E  Fst 热图   下三角标注数值 (Bizard: Correlation/Heatmap)
  F  全局统计摘要文本框 (含方法口径注)

图型参考 Bizard 图库:
  https://openbiox.github.io/Bizard/Composition/GroupedBarplot.html
  https://openbiox.github.io/Bizard/Hiplot/043-diverging-scale.html
  https://openbiox.github.io/Bizard/Correlation/Heatmap.html

用法:
  python -m utils.popgen_plots --stats popgen_stats.json \
      --pdf popgen_stats.pdf --png popgen_stats.png [--title GCVA]
"""

import argparse
import json
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Patch, Rectangle
import numpy as np

# 2026-10-16: 图件双格式 (PNG+PDF 同主名); 统一由 plot_io.save_both 收口。
# 双路径 import: `-m utils.popgen_plots` 与 `python utils/popgen_plots.py` 均可用。
try:
    from utils.plot_io import save_both
except ImportError:                              # pragma: no cover
    from plot_io import save_both

# Okabe-Ito 色盲安全调色板 (figure-maker 出版规范)
_OKABE_ITO = ['#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7', '#56B4E9', '#F0E442', '#000000']


def _setup_font():
    """地点/宿主分组值可能含中文, 有 CJK 字体就用。"""
    try:
        from matplotlib import font_manager
        names = {f.name for f in font_manager.fontManager.ttflist}
        for cand in ('Microsoft YaHei', 'SimHei', 'Noto Sans CJK SC',
                     'Source Han Sans SC', 'WenQuanYi Zen Hei', 'PingFang SC'):
            if cand in names:
                matplotlib.rcParams['font.sans-serif'] = [cand, 'DejaVu Sans']
                break
        matplotlib.rcParams['axes.unicode_minus'] = False
    except Exception:  # noqa: BLE001
        pass


def _sig_stars(q: float) -> str:
    """BH 校正后 q 值 → 显著性星号 (借 EasyHap plotting._sig_stars 约定)。"""
    if q is None or (isinstance(q, float) and math.isnan(q)):
        return ''
    return '***' if q < 0.001 else '**' if q < 0.01 else '*' if q < 0.05 else ''


def _short(s, width=10):
    """过长的分组名截断 (图轴标签防重叠)。"""
    s = str(s)
    return s if len(s) <= width else s[:width - 1] + '…'


def render_popgen_stats(stats, out_pdf, out_png=None, title=None, figsize=(17, 9)):
    """渲染 2×3 多面板统计图。

    stats: popgen_analysis --stats-json 的 dict (含 global/populations/fst/fst_pairwise)
    figsize: 画布英寸。默认 (17,9) 为独立整图; 作组合图面板时按目标面板宽出图
             (如 9.32 = Fig3 通栏宽), 保证图内字号 1:1 印刷不缩水。
    返回 {'success': bool, 'panels': {...}, 'error': None}
    """
    _setup_font()
    # ── 字号/间距自适应: 整图 (17in) 与组合图面板 (≈9.32in) 两档 ──
    # 面板宽度下 17in 版的标题字号必然互相重叠 (实测 A/B/C 与 D/E/F 标题互压),
    # 缩比例: 标题 11→8.5, suptitle 13→10, 行距/顶边距加大。
    panel_mode = figsize[0] < 12
    _fs_title = 8.5 if panel_mode else 11
    _fs_supt = 10 if panel_mode else 13
    _top, _hs = (0.86, 0.62) if panel_mode else (0.90, 0.42)
    pops = stats.get('populations') or []
    g = stats.get('global') or {}
    fst_pw = stats.get('fst_pairwise') or {}
    fst_g = stats.get('fst') or {}

    fig, axes = plt.subplots(2, 3, figsize=figsize, facecolor='white')
    fig.subplots_adjust(bottom=0.10, top=_top, wspace=0.28, hspace=_hs)
    panels_drawn = 0

    names = [_short(p.get('pop', '?')) for p in pops]
    x = np.arange(len(pops))
    # 真实数据可能 10+ 群体: 组多时降字号 + 加大旋转, 防标签挤压重叠
    _rot, _tfs = (45, 7.5) if len(pops) > 8 else (30, 8.5)
    span_eps = 0.0  # C 面板零值可见化阈值 (span 确定后更新)

    def _style(ax, panel):
        ax.set_title(panel, fontsize=_fs_title, loc='left', fontweight='bold')
        ax.tick_params(labelsize=8.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # ── A. π vs θ-W 分组柱状 (统一换算为 per-site 百分比, 量纲可比) ──
    # 原始 θ-W = S/a1 是**计数值**, 与 π (per site) 差 L 倍, 直接同图会把
    # π 柱压扁到不可见 → θ-W 除以 L 换算成每位点, 两者同单位 (%)
    ax = axes[0][0]
    pi_v = [(p.get('pi') or 0.0) * 100 for p in pops]
    th_v = [((p.get('theta_w') or 0.0) / p['L'] * 100) if p.get('L') else 0.0 for p in pops]
    if pops:
        w = 0.38
        ax.bar(x - w / 2, pi_v, w, label='π', color=_OKABE_ITO[0])
        ax.bar(x + w / 2, th_v, w, label='θ-W (per site)', color=_OKABE_ITO[1])
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=_rot, ha='right', fontsize=_tfs)
        ax.set_ylabel('diversity (%)', fontsize=9)
        ax.legend(fontsize=8, frameon=False)
        panels_drawn += 1
    else:
        ax.text(0.5, 0.5, 'no population data', ha='center', va='center', fontsize=9)
    _style(ax, 'A. Nucleotide diversity (π vs θ-W, %)')

    # ── B. Hd 柱状 (柱底标 n) ──
    ax = axes[0][1]
    hd_v = [p.get('Hd') or 0.0 for p in pops]
    if pops:
        ax.bar(x, hd_v, 0.6, color=_OKABE_ITO[2])
        ax.set_xticks(x)
        # 双行 (名+n=) 标签只在宽画布且组少时用; 面板模式/组多 → 单行短名
        if len(pops) > 8 or panel_mode:
            ax.set_xticklabels(names, rotation=_rot, ha='right', fontsize=_tfs)
        else:
            ax.set_xticklabels([f"{n}\n(n={p.get('n', '?')})" for n, p in zip(names, pops)],
                               rotation=_rot, ha='right', fontsize=_tfs)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel('haplotype diversity Hd', fontsize=9)
        for xi, v in zip(x, hd_v):
            ax.text(xi, min(v + 0.02, 1.03), f'{v:.3f}', ha='center', fontsize=7)
        panels_drawn += 1
    else:
        ax.text(0.5, 0.5, 'no population data', ha='center', va='center', fontsize=9)
    _style(ax, 'B. Haplotype diversity (Hd)')

    # ── C. Tajima's D 发散条形 ──
    ax = axes[0][2]

    def _num(v):
        return None if v is None or (isinstance(v, float) and not math.isfinite(v)) else v

    td_v = [_num(p.get('tajima_D')) for p in pops]
    if pops and any(v is not None for v in td_v):
        vals = [v if v is not None else 0.0 for v in td_v]
        colors = [_OKABE_ITO[3] if v < 0 else _OKABE_ITO[1] for v in vals]
        ax.bar(x, vals, 0.6, color=colors)
        ax.axhline(0, color='#333333', lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=_rot, ha='right', fontsize=_tfs)
        ax.set_ylabel("Tajima's D", fontsize=9)
        span = max([abs(v) for v in vals] + [0.5])
        # 零值可见化 (2026-10-13, 用户报告 "好多空白"): D≈0 的组画灰色细柱 + 加粗
        # 0.00 标签 —— 否则零高度柱看起来像 "没算"
        span_eps = max(span * 0.015, 1e-9)
        for xi, v in zip(x, vals):
            if abs(v) < span_eps:
                ax.bar(xi, span_eps, 0.6, color='#bbbbbb', zorder=1.5)
        for xi, v in zip(x, td_v):
            if v is not None:
                ax.text(xi, v + (span * 0.03 if v >= 0 else -span * 0.03),
                        f'{v:.2f}', ha='center', fontweight='bold' if abs(v) < span_eps else 'normal',
                        va='bottom' if v >= 0 else 'top', fontsize=7)
        panels_drawn += 1
    else:
        ax.text(0.5, 0.5, "Tajima's D unavailable (S=0 or n<2)",
                ha='center', va='center', fontsize=9)
    _style(ax, "C. Tajima's D (negative → expansion)")

    # ── D. Fu's Fs 发散条形 ──
    ax = axes[1][0]
    fs_v = [_num(p.get('fu_fs')) for p in pops]
    if pops and any(v is not None for v in fs_v):
        vals = [v if v is not None else 0.0 for v in fs_v]
        colors = [_OKABE_ITO[3] if v < 0 else _OKABE_ITO[1] for v in vals]
        ax.bar(x, vals, 0.6, color=colors)
        # None 可见化 (2026-10-13): k≥n (全唯一单倍型) 时 Fs 未定义 —— 画灰色
        # "n/a" 柱代替留白, 与 "没算" 区分
        span_pre = max([abs(v) for v in vals if v is not None] + [0.5])
        for xi, v in zip(x, fs_v):
            if v is None:
                ax.bar(xi, span_pre * 0.02, 0.6, color='#cccccc', zorder=1.5)
                ax.text(xi, span_pre * 0.04, 'n/a', ha='center', va='bottom',
                        fontsize=6.5, color='#666666', rotation=90)
        ax.axhline(0, color='#333333', lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=_rot, ha='right', fontsize=_tfs)
        ax.set_ylabel("Fu's Fs", fontsize=9)
        span = span_pre
        ax.set_ylim(min(min(vals), 0) - span * 0.15, max(max(vals), 0) + span * 0.15)
        for xi, v in zip(x, fs_v):
            if v is not None:
                ax.text(xi, v + (span * 0.03 if v >= 0 else -span * 0.03),
                        f'{v:.1f}', ha='center',
                        va='bottom' if v >= 0 else 'top', fontsize=7)
        panels_drawn += 1
    else:
        ax.text(0.5, 0.5, "Fu's Fs unavailable (n<2 or single haplotype)",
                ha='center', va='center', fontsize=9)
    _style(ax, "D. Fu's Fs, Ewens (negative → expansion)")

    # ── E. pairwise Fst 热图 (下三角) ──
    ax = axes[1][1]
    fp_pops = fst_pw.get('pops') or []
    mat = fst_pw.get('matrix') or []
    if len(fp_pops) >= 2 and mat:
        # 只显示下三角 (含对角 0 置灰): 上三角置 NaN 遮掉, 与标题一致
        m = np.array([[np.nan if (j >= i or v is None) else v
                       for j, v in enumerate(row)] for i, row in enumerate(mat)])
        masked = np.ma.masked_invalid(m)
        im = ax.imshow(masked, cmap='Reds', vmin=0, vmax=max(0.01, float(np.nanmax(m))))
        ax.set_xticks(range(len(fp_pops)))
        ax.set_yticks(range(len(fp_pops)))
        ax.set_xticklabels([_short(p, 8) for p in fp_pops], rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels([_short(p, 8) for p in fp_pops], fontsize=8)
        for i in range(len(fp_pops)):
            for j in range(len(fp_pops)):
                if j < i and not np.isnan(m[i][j]):
                    ax.text(j, i, f'{m[i][j]:.3f}', ha='center', va='center', fontsize=7,
                            color='white' if m[i][j] > 0.6 * float(np.nanmax(m)) else 'black')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        panels_drawn += 1
    else:
        ax.text(0.5, 0.5, 'pairwise Fst unavailable\n(group size < min_n?)',
                ha='center', va='center', fontsize=9)
    # 面板模式宽度下长标题会压到右邻面板标题 (实测 "lower triangle"∩"Summary")
    _style(ax, 'E. Pairwise Fst (lower tri.)' if panel_mode
           else 'E. Pairwise Fst (site-wise, lower triangle)')

    # ── F. 全局统计摘要 ──
    ax = axes[1][2]
    ax.axis('off')
    lines = ['Global (all samples, reference-excluded):', '']
    if g:
        fu = g.get('fu_fs')
        td = g.get('tajima_D')
        lines += [
            f"  samples n = {g.get('n', '?')}   sites L = {g.get('L', '?')}",
            f"  haplotypes = {g.get('haplotypes', '?')}   Hd = {g.get('Hd', '?')}",
            f"  π = {g.get('pi', '?')}   θ-W = {g.get('theta_w', '?')}   S = {g.get('S', '?')}",
            f"  Tajima's D = {td}",
            f"  Fu's Fs = {fu} (Ewens)",
        ]
    if fst_g.get('fst_site') is not None:
        lines += [f"  Fst(site, global) = {fst_g['fst_site']:.4f}  p = {fst_g.get('p')}",
                  f"  Fst(Nei haplotype) = {fst_g.get('fst_haplotype', {})}"]
    lines += ['', f"group column: {stats.get('group_col', '?')}   pops plotted: {len(pops)}",
              'methods: manual ACGT pairwise-deletion;', "Fu's Fs: Ewens (θ̂ from haplotype count)"]
    # ── EasyHap 借鉴补充 (2026-10-12): 私有单倍型 / 相似性聚类 / 富集星号 ──
    pv = stats.get('private_haplotypes') or {}
    if pv.get('n_private') is not None:
        _top = ', '.join(f"{k}={v.get('n_haps', 0)}" for k, v in
                         sorted((pv.get('by_group') or {}).items(),
                                key=lambda kv: -kv[1].get('n_haps', 0))[:3])
        lines.append(f"private haplotypes: {pv['n_private']}/{pv.get('n_total_haps', '?')} ({_top})")
    cl = stats.get('clusters') or {}
    if cl.get('n_clusters') is not None:
        try:
            lines.append(f"similarity clusters (≤{float(cl.get('threshold', 0)):.0%} diff): {cl['n_clusters']}")
        except (TypeError, ValueError):
            pass
    en = stats.get('enrichment') or {}
    sig = [r for r in (en.get('rows') or []) if r.get('star')][:3]
    if sig:
        lines.append('haplotype→location enrichment (BH):')
        for r in sig:
            lines.append(f"  {r['hap_id']} → {r['group']} ({r['in_group']}/{r['group_n']}) "
                         f"q={r['q']:g} {r.get('star', '')}")
    ax.text(0.02, 0.98, '\n'.join(lines), ha='left', va='top', fontsize=8.5,
            family='monospace', transform=ax.transAxes)
    panels_drawn += 1
    _style(ax, 'F. Summary')

    cap = 'Population genetics summary'
    if title:
        cap += f' — {title}'
    if stats.get('group_col'):
        cap += f" (grouped by {stats['group_col']})"
    fig.suptitle(cap, fontsize=_fs_supt, y=0.99 if panel_mode else 0.97)
    # 2026-10-16: 图件双格式 (PNG+PDF 同主名) —— 调用方没传 out_png 时按主名
    # 自动补 .png; 原先 PDF 那次 savefig 参数不变。
    save_both(fig, out_pdf, png=out_png, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'saved: {out_pdf} ({panels_drawn}/6 statistical panels)')
    return {'success': True, 'panels': panels_drawn, 'error': None}


# ── 基因结构带配色 (与 biosoft/EasyHap easyhap/plotting.py 的口径一致, 便于
#    两条管线出的基因结构图观感统一; 病毒基因无 UTR 注释时 UTR 段为空, 不报错) ──
GENE_CDS_COLOR = 'black'
GENE_UTR_COLOR = '#BFBFBF'
GENE_INTRON_COLOR = '#8C8C8C'
GENE_BAND_FACE = '#EEF3F8'
GENE_BAND_EDGE = '#333333'

# 连接线 (基因坐标 → 热图列) 与轨道的配色
LINK_COLOR = '#7F7F7F'
VARIANT_MARK_COLOR = '#C82525'   # 变异位点倒三角 (hastat 论文 Fig.2 构图)

# 群体/群体对配色 (Okabe-Ito 色盲安全; 论文 Fig.2 的三 ecotype 同为红/黄/蓝系)
GROUP_PALETTE = ['#D55E00', '#F0E442', '#0072B2', '#009E73',
                 '#CC79A7', '#56B4E9', '#E69F00', '#000000']
PAIR_PALETTE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7',
                '#56B4E9', '#F0E442', '#E69F00', '#000000']


def _resolve_group_colors(groups, group_colors=None):
    """群体 → 颜色 映射 (缺省按 GROUP_PALETTE 顺序取; 显式传入的优先)。"""
    out = dict(group_colors or {})
    palette_i = 0
    for g in groups:
        if g in out:
            continue
        out[g] = GROUP_PALETTE[palette_i % len(GROUP_PALETTE)]
        palette_i += 1
    return out


def _link_gene_to_heatmap(fig, gene_ax, heat_ax, gene_intervals,
                          gene_y=0.5, link_y=1.0):
    """从基因区间的列向热图对应列拉虚线 (借 hastat HapGene.draw_genotype
    的逐列连接线与 EasyHap _connect_gene_to_axis 的思想)。

    gene_ax 与 heat_ax 共享 x 轴 (列坐标), 用 ConnectionPatch 跨坐标系连线。
    gene_intervals: [(c, c), ...] 单列区间列表 (逐位点引线) 或 [(c0, c1), ...]
        区间 (只连两端); 调用方按位点密度决定传哪种。

    返回画出的线条数。
    """
    n = 0
    for c0, c1 in gene_intervals:
        if c1 < c0:
            continue
        cols = range(c0, c1 + 1) if c1 > c0 else (c0,)
        for cx in cols:
            con = ConnectionPatch(
                xyA=(cx, gene_y), coordsA=gene_ax.transData,
                xyB=(cx, link_y), coordsB=heat_ax.transData,
                color=LINK_COLOR, linewidth=0.6, linestyle=(0, (3, 2)),
                alpha=0.85, zorder=0, clip_on=False)
            fig.add_artist(con)
            n += 1
    return n


def _draw_gene_band(ax, gene_intervals, labels, xlim, strand='+'):
    """基因结构带: 外显子 (黑) / UTR (灰) 矩形 + 内含子箭头连线 + 基因名。

    gene_intervals: [(c0, c1, kind)] —— kind ∈ {'CDS', 'UTR'} (列坐标, 闭区间,
        与热图列对齐); 单基因注释只有起止时 kind='CDS' 整段涂黑。
    labels: 与区间一一对应的显示名 (只取第一个非空作基因名)。
    strand: '+'/'-' —— 内含子箭头方向 (借 hastat draw_introns 的链向感知;
        EasyHap _draw_gene_structure 同样按链翻转坐标, 这里热图列已按链排好序,
        只需箭头朝向正确)。
    """
    if not gene_intervals:
        return
    name = next((l for l in labels if l), 'gene')
    # 基因块垂直居中偏上 (0.32-0.68): 上方留 5′/3′ 标签位 (0.78+), 下方
    # 留变异位点倒三角位 (0.10, 由调用方画) —— 论文 Fig.2 B 面板布局。
    for c0, c1, kind in gene_intervals:
        color = GENE_CDS_COLOR if kind == 'CDS' else GENE_UTR_COLOR
        ax.add_patch(Rectangle((c0 - 0.5, 0.32), max(0.05, c1 - c0 + 1), 0.36,
                               facecolor=color, edgecolor='none',
                               linewidth=0, zorder=3, clip_on=False))
    # 内含子: 相邻区间之间的缺口画成带箭头线段, 方向随链向
    ordered = sorted(gene_intervals, key=lambda t: t[0])
    for (a0, a1, _), (b0, b1, _k2) in zip(ordered, ordered[1:]):
        gap0, gap1 = a1 + 0.5, b0 - 0.5
        if gap1 <= gap0:
            continue
        if strand == '-':
            gap0, gap1 = gap1, gap0
        ax.annotate('', xy=(gap1, 0.5), xytext=(gap0, 0.5),
                    arrowprops=dict(arrowstyle='-|>', color=GENE_INTRON_COLOR,
                                    lw=0.9, shrinkA=0, shrinkB=0,
                                    mutation_scale=7), zorder=2)
    ax.text(ordered[0][0] - 0.5, 0.80, f"5′ {name}", ha='left', va='bottom',
            fontsize=7.5, clip_on=False)
    ax.text(ordered[-1][1] + 0.5, 0.80, '3′', ha='right', va='bottom',
            fontsize=7.5, clip_on=False)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1)
    ax.axis('off')


def render_hap_state_heatmap(hap_rows, out_pdf, out_png=None, title=None,
                             genes=None, max_sites=400, max_token_sites=40,
                             tracks=None, dpi=300, hap_groups=None,
                             group_colors=None, max_track_groups=6):
    """单倍型 × 多态位点等位状态热图 (hastat 论文 Fig.2 构图, 2026-10-22)。

    自上而下五层, 全部共享基因组横轴:
      1. Fst 轨道 (按群体两两): 每条线一个群体对, 右端标对名 (论文 B/C 面板)
      2. π 轨道 (按群体): 每条线一个群体 (论文 B/C 面板)
      3. 基因结构带: CDS 黑块 / UTR 灰块 / 内含子箭头 (链向) + 基因名;
         下缘红色倒三角逐列标出**变异位点** (论文 B 面板的倒三角)
      4. 单倍型等位状态矩阵: 行按频次降序; 左侧行色条按该单倍型的**主导
         群体**着色 (论文 B 面板彩色方块); 位点多时为纯栅格, 少时画碱基字符
      5. 横轴: 多态位点编号 (5′→3′)

    hap_rows: [(hap_id, seq, count)] 按频次降序; 只统计 ACGT 列的多态位点。
    genes:    [(name, start, end)] 参考坐标 1-based 闭区间 (可选; 画基因带,
              并从每个基因两端引连接线到热图列区间)。
    tracks:   {'pi': {group: [(pos, val), ...]},
               'fst': {pair_label: [(pos, val), ...]}} 可选 —— 滑窗轨道,
              pos 与热图列同一坐标基准 (基因组坐标或比对列坐标)。
    hap_groups: {hap_id: group_label} 可选 —— 单倍型主导群体 (行色条与图例);
    group_colors: {group_label: color} 可选 —— 自定义群体配色。
    max_track_groups: 群体/群体对超过此数时轨道只画前 N 个 (防图例爆炸)。
    返回 {'success', 'n_haps', 'n_sites', 'subsampled', 'error'}
    """
    _setup_font()
    try:
        if not hap_rows:
            return {'success': False, 'error': 'no haplotypes', 'panels': 0}
        L = len(hap_rows[0][1])
        n_h = len(hap_rows)
        # 多态位点 (仅 ACGT 列且 ≥2 种碱基)
        poly = []
        for j in range(L):
            col = {s[j] for _, s, _ in hap_rows if s[j] in 'ACGT'}
            if len(col) >= 2:
                poly.append(j)
        subsampled = False
        if len(poly) > max_sites:
            idx = np.unique(np.linspace(0, len(poly) - 1, max_sites).astype(int))
            poly = [poly[i] for i in idx]
            subsampled = True
        ncol = len(poly)
        if ncol == 0:
            return {'success': False, 'error': 'no polymorphic ACGT sites', 'panels': 0}

        base_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        cmap = matplotlib.colors.ListedColormap(
            ['#4daf4a', '#377eb8', '#ff7f00', '#e41a1c', '#dddddd'])  # A C G T other
        mat = np.full((n_h, ncol), 4, dtype=np.int8)
        for r, (_, s, _c) in enumerate(hap_rows):
            for c, j in enumerate(poly):
                mat[r, c] = base_idx.get(s[j], 4)

        # ── 轨道: Fst (按群体两两, 上) + π (按群体, 下) ──
        # 坐标基准统一: sliding_tracks 的 pos 是**窗口中心** (基因组/比对列坐标,
        # 1-based), 而矩阵横轴是**多态位点序号** (0..ncol-1)。两者必须换算到
        # 同一基准才能真正 sharex —— 论文里 hastat 的轨道与矩阵共享横轴, 是因为
        # 它的窗口位置就落在基因组物理坐标且矩阵列按基因组排; 我们位点抽稀后
        # 列号与坐标不再线性对应, 故用 poly 反查: 窗口中心落在哪个多态列上。
        # 找不到 (窗口落在无多态位点的区段) 的点丢掉, 避免轨道画到错误列。
        import bisect
        tracks = tracks or {}

        def _to_col(pos):
            """基因组/比对列坐标 (1-based) → 最近的多态位点列号 (抽稀后)。

            窗口中心很少正好落在多态位点上; 取左侧或右侧最近的一个
            (位点抽稀后列距变大, 但"最近列"仍保持轨道的形状与位置关系)。
            """
            j = bisect.bisect_left(poly, pos - 1)
            cands = [c for c in (j - 1, j) if 0 <= c < ncol]
            if not cands:
                return None
            return min(cands, key=lambda c: abs(poly[c] - (pos - 1)))

        def _map_series(series_dict):
            out = []
            for key, pts in sorted(series_dict.items()):
                mapped = []
                for pos, v in sorted(pts):
                    c = _to_col(pos)
                    if c is not None:
                        mapped.append((c, v))
                if mapped:
                    out.append((key, mapped))
            return out

        fst_series = _map_series(tracks.get('fst') or {})
        pi_series = _map_series(tracks.get('pi') or {})
        if len(fst_series) > max_track_groups:
            fst_series = fst_series[:max_track_groups]
        if len(pi_series) > max_track_groups:
            pi_series = pi_series[:max_track_groups]
        # 群体 (行色条/图例) 汇总: 显式 hap_groups > π 轨道的群体键
        all_groups = []
        if hap_groups:
            all_groups = sorted({g for g in hap_groups.values() if g})
        elif pi_series:
            all_groups = [g for g, _ in pi_series]
        gcolors = _resolve_group_colors(all_groups, group_colors)

        # ── 画布布局: 自上而下 fst / pi / 基因带 / 矩阵, 全 sharex ──
        # 各层的"相对高度"按英寸折算成 ratios (matplotlib gridspec 只认比例):
        #   矩阵每行 row_h in; 轨道每条 track_h_each in; 基因带 band_h in。
        row_h = 0.16 if n_h > 60 else 0.24
        n_tracks = (1 if fst_series else 0) + (1 if pi_series else 0)
        track_h_each = 0.62                       # 每条轨道高度 (英寸)
        band_h = 0.62 if genes else 0.0           # 基因结构带高度 (英寸)
        fig_h = n_h * row_h + n_tracks * track_h_each + band_h + 0.55
        fig_w = min(20.0, max(6.0, ncol * 0.08 + 2))
        # 左侧行色条占 0.22 in (与 y 轴标签同一区域, 用 twinx 式的偏移轴实现)
        bar_w_in = 0.22 if all_groups else 0.0
        fig = plt.figure(figsize=(fig_w + bar_w_in, fig_h), facecolor='white')
        # ratios: 轨道/基因带用固定英寸数, 矩阵用 n_h*row_h (行数主导)
        ratios = ([track_h_each] * n_tracks
                  + ([band_h] if genes else [])
                  + [n_h * row_h])
        gs = fig.add_gridspec(len(ratios), 1, height_ratios=ratios,
                              hspace=0.12 if n_tracks else 0)
        axes = [fig.add_subplot(gs[i]) for i in range(len(ratios))]
        track_axes = axes[:n_tracks]
        ax_rest = axes[n_tracks:]
        gene_ax = ax_rest[0] if genes else None
        ax = ax_rest[-1] if ax_rest else axes[-1]
        for a in axes[:-1]:
            a.tick_params(labelbottom=False, bottom=False)

        # ── 1. Fst 轨道 (按群体两两, 论文 B/C 面板最上层) ──
        if fst_series:
            a = track_axes[0]
            for i, (pair, pts) in enumerate(fst_series):
                pts = sorted(pts)
                if not pts:
                    continue
                color = PAIR_PALETTE[i % len(PAIR_PALETTE)]
                a.plot([p[0] for p in pts], [p[1] for p in pts],
                       color=color, lw=0.9, label=pair)
            a.set_ylabel('Fst', fontsize=7.5)
            # 图例放轴内左上 (标题在轴上方居中, 不冲突); 群体对多时两列
            a.legend(frameon=False, fontsize=5.5,
                     ncol=2 if len(fst_series) > 4 else 1,
                     loc='upper left', handlelength=1.2, columnspacing=0.8,
                     borderaxespad=0.1, labelspacing=0.2)
            a.tick_params(labelsize=6)

        # ── 2. π 轨道 (按群体) ──
        if pi_series:
            a = track_axes[-1]
            for g, pts in pi_series:
                pts = sorted(pts)
                if not pts:
                    continue
                a.plot([p[0] for p in pts], [p[1] for p in pts],
                       color=gcolors.get(g, '#888888'), lw=0.9, label=str(g))
            a.set_ylabel('π (1e-2)', fontsize=7.5)
            a.legend(frameon=False, fontsize=5.5,
                     ncol=2 if len(pi_series) > 4 else 1,
                     loc='upper left', handlelength=1.2, columnspacing=0.8,
                     borderaxespad=0.1, labelspacing=0.2)
            a.tick_params(labelsize=6)

        # ── 3. 基因结构带 + 逐基因连接线 + 变异位点倒三角 ──
        if genes:
            import bisect
            gene_band_ivals = []
            for gname, gs, ge in genes:
                c0 = bisect.bisect_left(poly, gs - 1)
                c1 = bisect.bisect_right(poly, ge - 1) - 1
                if c1 < c0:
                    continue
                gene_band_ivals.append((c0, c1, 'CDS'))
                # 逐位点引线 (论文 Fig.2: 位点↔矩阵列一一对应), 但位点列数
                # 过多时线间距小于线宽会糊成灰幕 → 每基因区间内等间隔抽稀到
                # ≤60 条/基因; 位点很少时全画 (信息最完整)。
                span = c1 - c0 + 1
                if span <= 60:
                    cols = list(range(c0, c1 + 1))
                else:
                    stride = span // 60 + 1
                    cols = list(range(c0, c1 + 1, stride))
                _link_gene_to_heatmap(fig, gene_ax, ax,
                                      [(c, c) for c in cols],
                                      gene_y=0.5, link_y=1.0)
            _draw_gene_band(gene_ax, gene_band_ivals,
                            [g[0] for g in genes], (0 - 0.5, ncol - 0.5))
            gene_ax.add_patch(Rectangle((0 - 0.5, 0.02), ncol, 0.96,
                                        facecolor=GENE_BAND_FACE,
                                        edgecolor=GENE_BAND_EDGE, lw=0.8,
                                        zorder=0, clip_on=False))
            # 变异位点倒三角: 每列下缘一个红三角, 直指矩阵顶 (论文 B 面板)。
            # 画在带内底部 (y=0.10) 且 clip_on=False: 位点列与矩阵列 1:1 对齐,
            # 位点越密三角越密, 直观显示变异分布; 抽稀后三角也跟着抽稀。
            if ncol <= 200:                     # 太密时三角会糊成一片红线
                gene_ax.scatter(range(ncol), [0.10] * ncol, marker='v', s=7,
                                color=VARIANT_MARK_COLOR, zorder=6,
                                clip_on=False, linewidths=0)

        # ── 4. 单倍型矩阵 + 行色条 (主导群体) ──
        ax.imshow(mat, aspect='auto', interpolation='nearest', cmap=cmap,
                  vmin=0, vmax=4)
        # x 轴 tick: ≤40 列全画; 41-200 列每 5 列一个; 更多则纯栅格 (无 tick)
        if ncol <= 40:
            ax.set_xticks(range(ncol))
            ax.set_xticklabels([p + 1 for p in poly], fontsize=5, rotation=90)
        elif ncol <= 200:
            step_t = 5
            keep = list(range(0, ncol, step_t))
            ax.set_xticks(keep)
            ax.set_xticklabels([poly[k] + 1 for k in keep], fontsize=5,
                               rotation=90)
        else:
            ax.set_xticks([])
        ax.set_yticks(range(n_h))
        ax.set_yticklabels([f"{hid} ({c})" for hid, _s, c in hap_rows],
                           fontsize=5.5 if n_h > 60 else 7)
        ax.set_ylabel('haplotypes (by frequency)', fontsize=9)
        ax.set_xlabel('polymorphic sites (5\u2032→3\u2032)' +
                      ('; subsampled' if subsampled else ''), fontsize=9)
        # 行色条: 与矩阵等高、贴在 y 轴标签右侧的窄条
        if all_groups:
            pos = ax.get_position()
            bar_w_frac = min(0.012, 0.22 / max(fig_w, 1e-6))
            bar_ax = fig.add_axes([pos.x0 - bar_w_frac - 0.004, pos.y0,
                                   bar_w_frac, pos.height], sharey=ax)
            bar_ax.set_xlim(0, 1)
            bar_ax.set_xticks([])
            bar_ax.set_yticks([])
            for spine in bar_ax.spines.values():
                spine.set_visible(False)
            bar_ax.patch.set_alpha(0)
            for r, (hid, _s, _c) in enumerate(hap_rows):
                g = (hap_groups or {}).get(hid)
                if g:
                    bar_ax.add_patch(Rectangle((0, r - 0.5), 1, 1,
                                               facecolor=gcolors[g],
                                               edgecolor='none'))
            # 图例 (行色条 ↔ 群体): 挂在矩阵**下方** x 轴标签旁, 避免与基因带
            # 的 5′/3′ 标签抢位置 (2026-10-22 论文构图对齐)
            handles = [Patch(facecolor=gcolors[g], edgecolor='none', label=str(g))
                       for g in all_groups]
            ax.legend(handles=handles, frameon=False, fontsize=6.5,
                      loc='upper left', bbox_to_anchor=(0.0, -0.10),
                      ncol=min(6, len(all_groups)), handlelength=1.0,
                      columnspacing=0.8, borderaxespad=0)

        if title:
            # 标题挂最顶部的轴 (Fst 轨道或基因带), 不压矩阵标题区
            top_ax = track_axes[0] if track_axes else (gene_ax or ax)
            top_ax.set_title(title, fontsize=10, pad=4)
        # 2026-10-16: 同 render_popgen_stats —— PNG 缺失时按主名自动补 (.png 必出)
        save_both(fig, out_pdf, png=out_png, png_dpi=dpi,
                  bbox_inches='tight', facecolor='white')
        plt.close(fig)
        _tag = []
        if genes:
            _tag.append(f"{len(genes)} genes")
        if pi_series:
            _tag.append(f"π×{len(pi_series)}")
        if fst_series:
            _tag.append(f"Fst×{len(fst_series)}")
        print(f'saved: {out_pdf} ({n_h} haps x {ncol} sites'
              + (', subsampled' if subsampled else '')
              + (f", {', '.join(_tag)}" if _tag else '') + ')')
        return {'success': True, 'n_haps': n_h, 'n_sites': ncol,
                'subsampled': subsampled, 'error': None}
    except Exception as e:  # noqa: BLE001
        return {'success': False, 'n_haps': 0, 'n_sites': 0, 'error': str(e)}


def main():
    ap = argparse.ArgumentParser(description='群体遗传统计多面板图 (吃 popgen_analysis --stats-json)')
    ap.add_argument('--stats', required=True, help='popgen_analysis --stats-json 产出的 JSON')
    ap.add_argument('--pdf', required=True)
    ap.add_argument('--png', default=None)
    ap.add_argument('--title', default=None)
    ap.add_argument('--width', type=float, default=17.0,
                    help='画布宽 (英寸); 组合图面板按目标面板宽出图, 字号 1:1')
    ap.add_argument('--height', type=float, default=9.0)
    a = ap.parse_args()
    with open(a.stats, encoding='utf-8') as f:
        stats = json.load(f)
    render_popgen_stats(stats, a.pdf, a.png, a.title, figsize=(a.width, a.height))


if __name__ == '__main__':
    main()
