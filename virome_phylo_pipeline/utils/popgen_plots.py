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
import numpy as np

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
        span = max(abs(v) for v in vals) or 1.0
        ax.set_ylim(min(min(vals), 0) - span * 0.15, max(max(vals), 0) + span * 0.15)
        for xi, v in zip(x, td_v):
            if v is not None:
                ax.text(xi, v + (span * 0.03 if v >= 0 else -span * 0.03),
                        f'{v:.2f}', ha='center',
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
        ax.axhline(0, color='#333333', lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=_rot, ha='right', fontsize=_tfs)
        ax.set_ylabel("Fu's Fs", fontsize=9)
        span = max(abs(v) for v in vals) or 1.0
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
    fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    if out_png:
        fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'saved: {out_pdf} ({panels_drawn}/6 statistical panels)')
    return {'success': True, 'panels': panels_drawn, 'error': None}


def render_hap_state_heatmap(hap_rows, out_pdf, out_png=None, title=None,
                             genes=None, max_sites=400, max_token_sites=40,
                             dpi=300):
    """单倍型 × 多态位点等位状态热图 (借 EasyHap plot_haplotype_heatmap 思想,
    缩放到全基因组尺度: 位点多时为纯栅格; 可选基因结构带 + 连接线;
    位点少时格子内画碱基字符, 字号随格数自适应 —— EasyHap 的
    _cell_text_base_fontsize 思路)。

    hap_rows: [(hap_id, seq, count)] 按频次降序; 只统计 ACGT 列的多态位点。
    genes:    [(name, start, end)] 参考坐标 1-based 闭区间 (可选; 画顶部基因带,
              并从基因两端引连接线到热图列区间 —— EasyHap
              _connect_gene_to_axis 的思想)。
    max_sites: 位点数超限时等距抽稀 (纯栅格展示)。
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

        row_h = 0.16 if n_h > 60 else 0.24
        fig_h = n_h * row_h + (1.1 if genes else 0.7)
        fig, ax = plt.subplots(figsize=(min(17, max(6, ncol * 0.08 + 2)), fig_h),
                               facecolor='white')
        ax.imshow(mat, aspect='auto', interpolation='nearest', cmap=cmap,
                  vmin=0, vmax=4)
        ax.set_xticks([]) if ncol > 80 else ax.set_xticks(range(ncol))
        if ncol <= 80:
            ax.set_xticklabels([p + 1 for p in poly], fontsize=5, rotation=90)
        ax.set_yticks(range(n_h))
        ax.set_yticklabels([f"{hid} ({c})" for hid, _s, c in hap_rows],
                           fontsize=5.5 if n_h > 60 else 7)
        ax.set_ylabel('haplotypes (by frequency)', fontsize=9)
        ax.set_xlabel('polymorphic sites (5\u2032→3\u2032)' +
                      ('; subsampled' if subsampled else ''), fontsize=9)

        # ── 基因结构带 + 连接线 (借 EasyHap _connect_gene_to_axis 思想) ──
        if genes:
            from matplotlib.patches import Rectangle
            import bisect
            span_y = n_h - 0.5   # 热图顶
            strip_h = 0.9
            for gi, (gname, gs, ge) in enumerate(genes):
                c0 = bisect.bisect_left(poly, gs - 1)
                c1 = bisect.bisect_right(poly, ge - 1) - 1
                if c1 < c0:
                    continue
                y0 = span_y + 0.4 + strip_h
                ax.add_patch(Rectangle((c0 - 0.5, y0), c1 - c0 + 1, strip_h * 0.8,
                                       facecolor='#bdd7e7', edgecolor='#333333',
                                       lw=0.8, clip_on=False))
                ax.text((c0 + c1) / 2, y0 + strip_h * 0.4, gname, ha='center',
                        va='center', fontsize=7, clip_on=False)
                # 连接线: 基因两端 → 热图列边界
                for cx in (c0 - 0.5, c1 + 0.5):
                    ax.plot([cx, cx], [span_y + 0.3, y0], color='#666666',
                            lw=0.6, clip_on=False)
        if title:
            ax.set_title(title, fontsize=10, pad=12 if genes else 6)
        fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
        if out_png:
            fig.savefig(out_png, dpi=dpi, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f'saved: {out_pdf} ({n_h} haps x {ncol} sites'
              + (', subsampled' if subsampled else '') + ')')
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
