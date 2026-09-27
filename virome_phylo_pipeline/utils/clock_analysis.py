#!/usr/bin/env python3
"""
clock_analysis.py -- 分子钟联合分析脚本 (整体 + 基因级 RTT / LTT / DRT)

一条命令完成过去 rtt + temporal + genes 三个 stage 的全部工作:
  ① 整体 RTT   : TreeTime root-to-tip 回归 (β, R², p)
  ①b 稳健增强 (2026-08-27, 默认开启, 借鉴 shinyTempSignal 思想):
      - 学生化残差迭代剔除日期离群样本 (|t|>3), 报告剔除前后对比
      - 严格钟 GLS/PIC 回归 (系统发育独立对比, 校正非独立性)
  ② 整体 LTT   : TreeDater 分子定年 + LTT 曲线
  ③ 整体 DRT   : date-randomization 时间信号检验
  ④ 基因级 RTT : GenBank CDS 坐标切分 → 每基因 IQ-TREE + TreeTime-RTT
  ⑤ 基因级 LTT : 每基因 TreeDater (可选 --gene-ltt)
  ⑥ 基因级 DRT : 每基因随机化检验 (次数较低, 默认 10)
  ⑦ 汇总表     : clock_summary.tsv (整体与各基因的 β/R²/p/DRT 结论并排)

也可被 phylo_pipeline 的 clock stage 直接 import 调用 (run_clock_analysis)。

用法 (独立运行):
  python3 utils/clock_analysis.py \
      --alignment phylogeny/mafft.aln.fasta \
      --tree phylogeny/iqtree.treefile \
      --dates data/dates.csv --meta data/meta.csv \
      --gbk reference.gb --out-dir . --threads 16 \
      [--drt-randomizations 10] [--gene-drt-randomizations 10] \
      [--gene-ltt] [--rscript Rscript] [--iqtree iqtree2]

输出目录 (与原三个 stage 完全一致, 下游无感):
  time/treetime_rtt/  time/treedater_ltt/  time/temporal_signal/
  phylogeny/gene_analysis/{alignments,trees,drt}/
  time/clock_summary.tsv
"""
import argparse
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# 跨管线统一 I/O 布局: ⑤ per-virus 模块目录名随 MMPV_IO_LAYOUT 解析
try:
    from mmpv_common.io_layout import ph_dir
except ImportError:  # 脱离包环境直接运行时自举
    import os as _os, sys as _sys
    _REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    if _REPO_ROOT not in _sys.path:
        _sys.path.insert(0, _REPO_ROOT)
    from mmpv_common.io_layout import ph_dir

UTILS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(UTILS_DIR))
sys.path.insert(0, str(UTILS_DIR.parent))

from virphy_bridge import (          # noqa: E402
    LogCollector,
    run_treetime_rtt,
    run_treedater_ltt,
    slice_genes_from_alignment,
    run_gene_rtt_batch,
)
from temporal_signal import date_randomization_test   # noqa: E402


# ══════════════════════════════════════════════════════════════════════
# 稳健时间信号层 (借鉴 shinyTempSignal, 纯 Python 实现, 无新依赖)
#   · 学生化残差剔除离群日期标签 (迭代, |t| > 阈值)
#   · 严格钟 GLS 等价估计: 系统发育独立对比 (PIC) 回归
# ══════════════════════════════════════════════════════════════════════

def _betacf(a, b, x):
    """不完全 beta 连分数 (Lentz), 用于 t 分布 p 值"""
    MAXIT, EPS, FPMIN = 200, 3e-12, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        dele = d * c
        h *= dele
        if abs(dele - 1.0) < EPS:
            break
    return h


def _t_sf(t: float, df: int) -> float:
    """t 分布双尾 p 值 (via incomplete beta)"""
    if df <= 0:
        return float('nan')
    x = df / (df + t * t)
    ib = _betacf(df / 2.0, 0.5, x)
    p = ib * math.exp(math.lgamma(df / 2.0) + math.lgamma(0.5)
                      - math.lgamma((df + 1) / 2.0)) / math.sqrt(df)
    return min(1.0, max(0.0, p))


def _parse_newick(text: str):
    """极简 Newick 解析: 返回 (edges: {node: parent}, blen: {(p,c): len}, tip: {name: node}, root)"""
    text = text.strip().rstrip(';')
    edges, blen, tip = {}, {}, {}
    counter = [0]

    def new_internal():
        counter[0] -= 1
        return counter[0]

    def parse(s, parent):
        s = s.strip()
        if s.startswith('('):
            close = s.rindex(')')
            inner = s[1:close]              # 括号内内容
            label_part = s[close + 1:].strip()  # 节点标签[:blen]
            # 括号内顶层逗号切分 (depth==0)
            parts, depth, start = [], 0, 0
            for j, ch in enumerate(inner):
                if ch == '(':
                    depth += 1
                elif ch == ')':
                    depth -= 1
                elif ch == ',' and depth == 0:
                    parts.append(inner[start:j]); start = j + 1
            parts.append(inner[start:])
            node = new_internal()
            if ':' in label_part:
                try:
                    blen[(parent, node)] = float(label_part.split(':', 1)[1].strip())
                except ValueError:
                    pass
            edges[node] = parent
            for p_ in parts:
                parse(p_, node)
            return node
        else:
            # leaf: name[:blen]
            name, bl = s, None
            for m in re.finditer(r"'([^']*)'", s):
                name = m.group(1)
            if ':' in s:
                nm, bl_s = s.rsplit(':', 1)
                name = nm.strip().strip("'")
                try:
                    bl = float(bl_s.strip())
                except ValueError:
                    bl = None
            name = name.strip().strip("'")
            node = name
            tip[name] = node
            edges[node] = parent
            if bl is not None and parent is not None:
                blen[(parent, node)] = bl
            return node

    root = parse(text, None)
    return edges, blen, tip, root


def root_to_tip_distances(tree: str) -> Optional[Dict[str, float]]:
    """从 Newick 树算每个 tip 的 root-to-tip 距离 (含 branch lengths)"""
    txt = Path(tree).read_text() if os.path.exists(tree) else tree
    if not txt.strip():
        return None
    edges, blen, tip, root = _parse_newick(txt)
    if not tip:
        return None
    out = {}
    for name in tip:
        d, cur, seen = 0.0, name, set()
        while edges.get(cur) is not None and cur not in seen:
            seen.add(cur)
            p = edges[cur]
            d += blen.get((p, cur), 0.0)
            cur = p
        out[name] = d
    return out


def _ols(x, y):
    """OLS: 返回 slope/intercept/r2/p/mse"""
    n = len(x)
    if n < 3:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    if sxx <= 0:
        return None
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    slope = sxy / sxx
    intercept = my - slope * mx
    ss_res = sum((b - (intercept + slope * a)) ** 2 for a, b in zip(x, y))
    ss_tot = sum((b - my) ** 2 for b in y)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    df = n - 2
    mse = ss_res / df if df > 0 else float('nan')
    if mse != mse or mse <= 0:
        return {'slope': slope, 'intercept': intercept, 'r2': r2, 'p': None,
                'mse': mse, 'n': n}
    se = math.sqrt(mse / sxx)
    tval = slope / se if se > 0 else float('nan')
    p = _t_sf(abs(tval), df) if tval == tval else None
    return {'slope': slope, 'intercept': intercept, 'r2': r2, 'p': p,
            'mse': mse, 'n': n}


def greedy_r2_trim(x, y, names, min_gain=0.1, max_drop_frac=0.2):
    """贪心 R² 提升: 每轮尝试删每个点, 删谁后剩余 OLS R² 提升最大且 ≥min_gain 就删谁。
    解决高杠杆错标点对残差的 masking 效应 (学生化残差单独看不到它们)。
    返回 (kept_idx, dropped:[{name,r2_gain}], final_ols)"""
    idx = list(range(len(x)))
    dropped = []
    max_drop = max(1, int(len(x) * max_drop_frac))
    while len(idx) >= 8 and len(dropped) < max_drop:
        base_fit = _ols([x[i] for i in idx], [y[i] for i in idx])
        if not base_fit or base_fit['r2'] != base_fit['r2']:
            break
        best_gain, best_i = 0.0, None
        for i in idx:
            rest = [j for j in idx if j != i]
            if len(rest) < 5:
                continue
            f2 = _ols([x[j] for j in rest], [y[j] for j in rest])
            if f2 and f2['r2'] == f2['r2']:
                gain = f2['r2'] - base_fit['r2']
                if gain > best_gain:
                    best_gain, best_i = gain, i
        if best_i is None or best_gain < min_gain:
            break
        dropped.append({'name': names[best_i], 'r2_gain': round(best_gain, 4),
                        'method': 'greedy_r2'})
        idx.remove(best_i)
    final = _ols([x[i] for i in idx], [y[i] for i in idx]) if len(idx) >= 3 else None
    return idx, dropped, final


def robust_rtt_outliers(x, y, names, threshold=3.0):
    """两层串联: 贪心 R² 提升抓高杠杆 masking 点 + 学生化残差抓普通离群点。
    返回 (kept_idx, dropped:[{...}], final_ols)"""
    idx1, dropped1, _ = greedy_r2_trim(x, y, names)
    x2 = [x[i] for i in idx1]; y2 = [y[i] for i in idx1]
    n2 = [names[i] for i in idx1]
    idx2, dropped2, final, _n = studentized_outlier_trim(x2, y2, n2, threshold=threshold)
    kept = [idx1[i] for i in idx2]
    return kept, dropped1 + dropped2, final


def studentized_outlier_trim(x, y, names, threshold=3.0, max_drop_frac=0.15):
    """迭代学生化残差剔除: |t_i| > threshold 的样本逐个剔除 (每次删最离群一个)
    返回 (kept_idx, dropped:[{name,t,resid}], final_ols, iterations)"""
    idx = list(range(len(x)))
    dropped = []
    max_drop = max(1, int(len(x) * max_drop_frac))
    while len(idx) >= 5 and len(dropped) < max_drop:
        cx = [x[i] for i in idx]; cy = [y[i] for i in idx]
        fit = _ols(cx, cy)
        if not fit or fit['mse'] != fit['mse']:
            break
        n = len(idx); df = n - 2
        mx = sum(cx) / n
        sxx = sum((a - mx) ** 2 for a in cx)
        if sxx <= 0:
            break
        worst_t, worst_i = 0.0, None
        for k, i in enumerate(idx):
            resid = cy[k] - (fit['intercept'] + fit['slope'] * cx[k])
            h = 1.0 / n + (cx[k] - mx) ** 2 / sxx   # 杠杆值
            denom = math.sqrt(fit['mse'] * (1.0 - h))
            if denom <= 0:
                continue
            t_i = abs(resid) / denom
            if t_i > worst_t:
                worst_t, worst_i = t_i, i
        if worst_i is None or worst_t <= threshold:
            break
        k = idx.index(worst_i)
        resid = cy[k] - (fit['intercept'] + fit['slope'] * cx[k])
        dropped.append({'name': names[worst_i], 't': round(worst_t, 3),
                        'resid': round(resid, 5)})
        idx.remove(worst_i)
    final = _ols([x[i] for i in idx], [y[i] for i in idx]) if len(idx) >= 3 else None
    return idx, dropped, final, len(dropped)


def pic_regression(tree: str, dates: Dict[str, float]):
    """严格钟 GLS 等价: 系统发育独立对比 (Felsenstein 1985)
    双性状 (date, root-to-tip) contrasts → 过原点回归斜率 = 演化速率
    返回 {rate, r_pic, n_contrasts}; 思想同 ape::pic + nlme::gls(corBrownian)"""
    rtt = root_to_tip_distances(tree)
    if not rtt:
        return None
    common = [t for t in rtt if t in dates]
    if len(common) < 4:
        return None
    edges, blen, tip, root = _parse_newick(
        Path(tree).read_text() if os.path.exists(tree) else tree)
    children = {}
    for c, p in edges.items():
        if p is not None:
            children.setdefault(p, []).append(c)
    val = {t: (dates[t], rtt[t], 0.0, 1) for t in common}
    contrasts = []

    def collect(node):
        if node in val:
            return val[node]
        sub = [collect(c) for c in children.get(node, [])]
        sub = [s for s in sub if s is not None]
        if len(sub) < 2:
            return sub[0] if len(sub) == 1 else None
        (d1, r1, v1, n1), (d2, r2, v2, n2) = sub[0], sub[1]
        kids = [c for c in children.get(node, [])]
        l1 = blen.get((node, kids[0]), 0.0) + v1
        l2 = blen.get((node, kids[1]), 0.0) + v2
        if l1 + l2 <= 0:
            l1 = l2 = 1e-9
        w1, w2 = l2 / (l1 + l2), l1 / (l1 + l2)
        v_new = l1 * l2 / (l1 + l2)
        if v_new > 0:
            contrasts.append(((d1 - d2) / math.sqrt(v_new),
                              (r1 - r2) / math.sqrt(v_new)))
        val[node] = (w1 * d1 + w2 * d2, w1 * r1 + w2 * r2, v_new, n1 + n2)
        return val[node]

    try:
        collect(root)
    except RecursionError:
        return None
    cs = contrasts
    if len(cs) < 3:
        return None
    xs = [c[0] for c in cs]; ys = [c[1] for c in cs]
    # 过原点回归
    sxx = sum(a * a for a in xs)
    if sxx <= 0:
        return None
    rate = sum(a * b for a, b in zip(xs, ys)) / sxx
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    cov = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    sy = math.sqrt(sum((b - my) ** 2 for b in ys))
    r_pic = cov / (sx * sy) if sx > 0 and sy > 0 else float('nan')
    return {'rate': rate, 'r_pic': r_pic, 'n_contrasts': len(cs)}


def _fmt(v, fmt='{:.6g}'):
    if v is None:
        return 'NA'
    try:
        return fmt.format(v)
    except (TypeError, ValueError):
        return str(v)


def _load_dates(dates_csv: str) -> Optional[Dict[str, float]]:
    """读 dates CSV (name,date 两列, 兼容常见表头), 返回 {tip: decimal_date}"""
    import csv as _csv
    if not dates_csv or not os.path.exists(dates_csv):
        return None
    out = {}
    with open(dates_csv) as f:
        rdr = _csv.DictReader(f)
        cols = rdr.fieldnames or []
        name_col = next((c for c in cols if c.lower() in ('name', 'id', 'taxon', 'accession', 'tip', 'seq', 'sequence')), cols[0] if cols else None)
        date_col = next((c for c in cols if c.lower() in ('date', 'sampling_date', 'date_decimal', 'decimal_date', 'year')), None)
        if not name_col or not date_col:
            return None
        for row in rdr:
            nm, dt = (row.get(name_col) or '').strip(), (row.get(date_col) or '').strip()
            if not nm or not dt:
                continue
            # 兼容 YYYY-MM-DD / YYYY.MM.DD / YYYY/MM/DD → decimal
            m = re.match(r'(\d{4})(?:[-/.](\d{1,2}))?(?:[-/.](\d{1,2}))?', dt)
            if m:
                y = int(m.group(1)); mo = int(m.group(2) or 1); d_ = int(m.group(3) or 1)
                try:
                    import datetime as _dt
                    base = _dt.date(y, 1, 1)
                    dec = y + (base.toordinal() +
                               (_dt.date(y, mo, min(d_, 28)).toordinal() - base.toordinal())) / 365.25 - 1
                    # 简化: 年 + 年内天数比例
                    doy = (_dt.date(y, mo, min(d_, 28)) - base).days
                    dec = y + doy / 365.25
                except ValueError:
                    continue
                out[nm] = dec
            else:
                try:
                    out[nm] = float(dt)
                except ValueError:
                    continue
    return out or None


def _plot_robust_rtt(xs, ys, names, dropped_names, fit_before, fit_after,
                    out_pdf, out_png, dpi=300):
    """稳健 RTT 对比图 (借鉴 shinyTempSignal 图设计):
    灰点+灰虚线=全样本, 蓝/橙点+实线=剔除后, 离群样本红色标出
    SCI 风格: 白底, Okabe-Ito 配色, 全英文标签"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    dropped_set = set(dropped_names)
    fig, ax = plt.subplots(figsize=(5.2, 4.0), dpi=dpi)
    # 全样本 (灰)
    ax.scatter(xs, ys, s=22, c='#BBBBBB', edgecolors='#888888',
               linewidths=0.4, zorder=2, label='All samples')
    if fit_before:
        x0, x1 = min(xs), max(xs)
        y0 = fit_before['intercept'] + fit_before['slope'] * x0
        y1 = fit_before['intercept'] + fit_before['slope'] * x1
        ax.plot([x0, x1], [y0, y1], '--', c='#999999', lw=1.1, zorder=3,
                label=f"All: R²={fit_before.get('r2', float('nan')):.3f}")
    # 剔除后 (蓝) + 离群点 (红)
    kx = [(x, y, n) for x, y, n in zip(xs, ys, names) if n not in dropped_set]
    dx = [x for x, y, n in zip(xs, ys, names) if n in dropped_set]
    dy = [y for x, y, n in zip(xs, ys, names) if n in dropped_set]
    if kx:
        ax.scatter([t[0] for t in kx], [t[1] for t in kx], s=24,
                   c='#0072B2', edgecolors='white', linewidths=0.4,
                   zorder=4, label='Retained')
        if fit_after:
            x0, x1 = min(t[0] for t in kx), max(t[0] for t in kx)
            y0 = fit_after['intercept'] + fit_after['slope'] * x0
            y1 = fit_after['intercept'] + fit_after['slope'] * x1
            ax.plot([x0, x1], [y0, y1], '-', c='#D55E00', lw=1.4, zorder=5,
                    label=f"Trimmed: R²={fit_after.get('r2', float('nan')):.3f}")
    if dx:
        ax.scatter(dx, dy, s=30, c='#D55E00', marker='X',
                   edgecolors='white', linewidths=0.4, zorder=6,
                   label=f'Outliers (n={len(dx)})')
    ax.set_xlabel('Sampling date (decimal year)', fontsize=10)
    ax.set_ylabel('Root-to-tip divergence', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.legend(fontsize=8, frameon=False, loc='best')
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=dpi)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    return True


def run_clock_analysis(alignment: str, tree: str, dates_csv: str,
                       out_dir: str, meta_csv: str = None,
                       gbk: str = None, seq_len: int = 10000,
                       threads: int = 40, iqtree_bin: str = 'iqtree2',
                       rscript: str = 'Rscript',
                       drt_randomizations: int = 20,
                       gene_drt_randomizations: int = 10,
                       gene_ltt: bool = False,
                       skip_genome: bool = False, skip_genes: bool = False,
                       outlier_threshold: float = 3.0,
                       log=None) -> Dict:
    """整体+基因级 RTT/LTT/DRT。返回 {summary_rows, metrics, outputs}。"""
    log = log or LogCollector()
    os.makedirs(out_dir, exist_ok=True)
    rows: List[Dict] = []
    metrics: Dict = {}
    outputs: Dict = {}

    have_dates = dates_csv and os.path.exists(dates_csv)
    if not have_dates:
        log.emit('[clock] 缺 dates CSV, DRT/LTT 将跳过')

    # ── ① 整体 RTT ──
    if not skip_genome:
        rtt_dir = os.path.join(out_dir, ph_dir('time'), 'treetime_rtt')
        if have_dates and tree and os.path.exists(tree):
            log.emit('[clock] ① 整体 RTT (TreeTime)...')
            r = run_treetime_rtt(fasta_file=alignment, tree_file=tree,
                                 dates_file=dates_csv, output_dir=rtt_dir, log=log)
            metrics.update({k: r.get(k) for k in ('beta', 'r_squared', 'p_value')})
            outputs['rtt_plot'] = r.get('plots', {}).get('rtt')
            outputs['timetree'] = r.get('output_tree')

            # ── ①b 稳健层: 学生化残差剔除 + PIC 严格钟 (默认开启) ──
            try:
                rtt_d = root_to_tip_distances(tree)
                dates_map = _load_dates(dates_csv)
                if rtt_d and dates_map:
                    common = [t for t in rtt_d if t in dates_map]
                    if len(common) >= 8:
                        xs = [dates_map[t] for t in common]
                        ys = [rtt_d[t] for t in common]
                        # 两层剔除: 贪心 R² + 学生化残差 (默认开启)
                        idx, dropped, final_fit = robust_rtt_outliers(
                            xs, ys, common, threshold=outlier_threshold)
                        robust_dir = os.path.join(out_dir, ph_dir('time'), 'treetime_rtt', 'robust')
                        os.makedirs(robust_dir, exist_ok=True)
                        rep = {'n_total': len(common), 'n_dropped': len(dropped),
                               'dropped': dropped,
                               'ols_before': _ols(xs, ys),
                               'ols_after': final_fit}
                        if dropped:
                            log.emit(f"[clock] ①b 离群日期样本剔除 {len(dropped)} 条: "
                                     f"{', '.join(d['name'] for d in dropped)}")
                        # PIC 严格钟 (全样本)
                        pic = pic_regression(tree, dates_map)
                        if pic:
                            rep['pic_strict_clock'] = pic
                            metrics['pic_rate'] = pic.get('rate')
                            metrics['pic_r'] = pic.get('r_pic')
                            log.emit(f"[clock] ①b PIC 严格钟: rate={_fmt(pic.get('rate'))} "
                                     f"r={_fmt(pic.get('r_pic'))} "
                                     f"(n={pic.get('n_contrasts')} contrasts)")
                        import json as _json
                        with open(os.path.join(robust_dir, 'robust_report.json'), 'w') as f:
                            _json.dump(rep, f, indent=2, default=str)
                        outputs['robust_report'] = os.path.join(robust_dir, 'robust_report.json')
                        # 稳健 RTT 对比图 (借鉴 shinyTempSignal 设计)
                        ok = _plot_robust_rtt(
                            xs, ys, common,
                            [d['name'] for d in dropped],
                            rep['ols_before'], final_fit,
                            os.path.join(robust_dir, 'robust_rtt_compare.pdf'),
                            os.path.join(robust_dir, 'robust_rtt_compare.png'))
                        if ok:
                            outputs['robust_rtt_plot'] = os.path.join(robust_dir, 'robust_rtt_compare.png')
                            log.emit('[clock] ①b 稳健 RTT 对比图 → robust_rtt_compare.{pdf,png}')
            except Exception as e:
                log.emit(f'[clock] ①b 稳健层失败(不阻断): {e}')

            rows.append({'level': 'genome', 'gene': '-', 'length': '-',
                         'beta': r.get('beta'), 'r_squared': r.get('r_squared'),
                         'p_value': r.get('p_value'), 'drt_passed': '',
                         'drt_pct': ''})

            # ── ② 整体 LTT ──
            ltt_dir = os.path.join(out_dir, ph_dir('time'), 'treedater_ltt')
            log.emit('[clock] ② 整体 LTT (TreeDater)...')
            d = run_treedater_ltt(tree_file=r.get('output_tree') or tree,
                                  metadata_file=dates_csv, seq_len=seq_len,
                                  output_dir=ltt_dir, plot_ltt=True,
                                  threads=threads,
                                  rscript_path=rscript)
            outputs['dated_phylo'] = d.get('plots', {}).get('phylogeny')
            outputs['ltt_plot'] = d.get('plots', {}).get('ltt')
            # ── 2026-08-31 treedater 四功能 metrics: CoV/钟检验/离群谱系/CI ──
            _tdm = d.get('metrics') or {}
            for _k, _lbl in [('cov_rate', 'td_cov_rate'),
                             ('relaxed_clock_p', 'td_relaxed_clock_p'),
                             ('tmrca', 'td_tmrca'),
                             ('mean_rate', 'td_mean_rate'),
                             ('n_outlier_tips', 'td_n_outlier_tips')]:
                if _k in _tdm:
                    metrics[_lbl] = _tdm[_k]
            if 'tmrca_ci' in _tdm:
                metrics['td_tmrca_ci95'] = f"{_tdm['tmrca_ci'][0]:.3f}-{_tdm['tmrca_ci'][1]:.3f}"
            if 'rate_ci' in _tdm:
                metrics['td_rate_ci95'] = f"{_tdm['rate_ci'][0]:.3e}-{_tdm['rate_ci'][1]:.3e}"
            if d.get('outlier_tips'):
                outputs['treedater_outlier_tips'] = d['outlier_tips']

            # ── ③ 整体 DRT ──
            drt_dir = os.path.join(out_dir, ph_dir('time'), 'temporal_signal')
            log.emit('[clock] ③ 整体 DRT (x%d)...' % drt_randomizations)
            g = date_randomization_test(fasta_file=alignment, dates_csv=dates_csv,
                                        output_dir=drt_dir,
                                        n_randomizations=drt_randomizations,
                                        threads=threads,
                                        tree_file=tree, log=log)
            metrics['temporal_signal'] = g.get('conclusion', '')
            metrics['drt_passed'] = g.get('passed', False)
            if rows:
                rows[-1].update({'drt_passed': g.get('passed'),
                                 'drt_pct': g.get('r2_percentile')})
        else:
            log.emit('[clock] 整体级跳过 (缺树或日期)')

    # ── ④⑤⑥ 基因级 ──
    if not skip_genes and gbk and os.path.exists(gbk) and alignment:
        gene_dir = os.path.join(out_dir, ph_dir('phylogeny'), 'gene_analysis')
        log.emit('[clock] ④ 基因级 RTT (切分 + IQ-TREE + TreeTime)...')
        gene_alns = slice_genes_from_alignment(gbk, alignment,
                                               os.path.join(gene_dir, 'alignments'),
                                               log=log)
        if gene_alns:
            gene_results = run_gene_rtt_batch(
                gene_alns, meta_csv or dates_csv, os.path.join(gene_dir, 'trees'),
                iqtree_bin=iqtree_bin, threads=threads, log=log)
            drt_dates = dates_csv if have_dates else None
            # 文件名映射: 切分产物是 {gene}.mafft.fasta, 必须从 gene_alns 拿真实路径
            aln_path_by_gene = {ga['gene']: ga['fasta'] for ga in gene_alns}
            # 逐基因 DRT: 各基因独立, 基因间并行 (2026-08-31 加速)
            def _gene_drt(gr):
                gname = gr.get('gene', '')
                g_aln = aln_path_by_gene.get(gname) or os.path.join(
                    gene_dir, 'alignments', f'{gname}.fasta')
                if not (g_aln and os.path.exists(g_aln)):
                    return gname, None, None, f'skip: 找不到比对文件 ({g_aln})'
                n_seq = sum(1 for l in open(g_aln) if l.startswith('>'))
                if not (10 <= n_seq <= 500):
                    return gname, None, None, f'skip: n={n_seq} 不在 [10,500]'
                g_tree = os.path.join(gene_dir, 'trees', gname, f'{gname}.treefile')
                if not os.path.exists(g_tree):
                    g_tree = os.path.join(gene_dir, 'trees', f'{gname}.treefile')
                try:
                    gd = date_randomization_test(
                        fasta_file=g_aln, dates_csv=drt_dates,
                        output_dir=os.path.join(gene_dir, 'drt', gname),
                        n_randomizations=gene_drt_randomizations,
                        threads=max(1, threads // max(1, len(gene_results))),
                        tree_file=g_tree if os.path.exists(g_tree) else None,
                        log=None)
                    return gname, gd.get('passed'), gd.get('r2_percentile'), None
                except Exception as e:
                    return gname, None, None, f'DRT 失败(不阻断): {e}'

            if drt_dates:
                from concurrent.futures import ThreadPoolExecutor
                _gd_thr = max(1, threads // max(1, len(gene_results)) or 1)
                with ThreadPoolExecutor(max_workers=min(_gd_thr, len(gene_results))) as _ex:
                    _gd_res = list(_ex.map(_gene_drt, gene_results))
                _gd_map = {gn: (pv, pct, msg) for gn, pv, pct, msg in _gd_res}
                for gn, pv, pct, msg in _gd_res:
                    if msg:
                        log.emit(f'[clock] {gn} DRT {msg}')
            else:
                _gd_map = {}
            for gr in gene_results:
                gname = gr.get('gene', '')
                row = {'level': 'gene', 'gene': gname,
                       'length': gr.get('length'),
                       'beta': gr.get('beta'), 'r_squared': gr.get('r_squared'),
                       'p_value': gr.get('p_value'), 'drt_passed': '', 'drt_pct': ''}
                if drt_dates and gr.get('success'):
                    pv, pct, _msg = _gd_map.get(gname, (None, None, None))
                    row['drt_passed'] = pv
                    row['drt_pct'] = pct
                    metrics[f'drt_{gname}_passed'] = pv if pv is not None else False
                    metrics[f'drt_{gname}_pct'] = pct
                # 基因级 LTT (可选)
                if gene_ltt and drt_dates and gr.get('success'):
                    g_tree2 = os.path.join(gene_dir, 'trees', gname, f'{gname}.treefile')
                    if not os.path.exists(g_tree2):
                        g_tree2 = os.path.join(gene_dir, 'trees', f'{gname}.treefile')
                    if os.path.exists(g_tree2):
                        try:
                            run_treedater_ltt(tree_file=g_tree2,
                                              metadata_file=drt_dates, threads=threads,
                                              seq_len=int(gr.get('length') or 1000),
                                              output_dir=os.path.join(gene_dir, 'ltt', gname),
                                              plot_ltt=True, rscript_path=rscript)
                        except Exception as e:
                            log.emit(f'[clock] {gname} LTT 失败(不阻断): {e}')
                rows.append(row)
                if gr.get('r_squared') is not None:
                    metrics[f'R2_{gname}'] = gr['r_squared']
            outputs['gene_analysis_dir'] = gene_dir
        else:
            log.emit('[clock] 无基因比对产出 — 参考 GBK 无 CDS (类病毒不适用) 或 gbk CDS 注释缺失')
    elif not skip_genes:
        log.emit('[clock] 基因级跳过 (未提供 --gbk)')

    # ── ⑦ 汇总表 ──
    summary = os.path.join(out_dir, ph_dir('time'), 'clock_summary.tsv')
    os.makedirs(os.path.dirname(summary), exist_ok=True)
    with open(summary, 'w') as f:
        f.write('level\tgene\tlength\tbeta\tr_squared\tp_value\tdrt_passed\tdrt_r2_percentile\n')
        for r in rows:
            f.write('\t'.join(str(r.get(k, '')) if r.get(k) is not None else 'NA'
                              for k in ('level', 'gene', 'length', 'beta',
                                        'r_squared', 'p_value', 'drt_passed',
                                        'drt_pct')) + '\n')
    outputs['clock_summary'] = summary
    log.emit(f'[clock] 汇总表 → {summary} ({len(rows)} 行)')
    return {'summary_rows': rows, 'metrics': metrics, 'outputs': outputs,
            'success': bool(rows)}


def main():
    ap = argparse.ArgumentParser(description='整体+基因级 RTT/LTT/DRT 联合分析')
    ap.add_argument('--alignment', required=True)
    ap.add_argument('--tree', required=True)
    ap.add_argument('--dates', required=True)
    ap.add_argument('--meta', default=None)
    ap.add_argument('--gbk', default=None, help='GenBank 注释 (CDS 坐标切分基因)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--seq-len', type=int, default=10000)
    ap.add_argument('--threads', type=int, default=40)
    ap.add_argument('--iqtree', default='iqtree2')
    ap.add_argument('--rscript', default='Rscript')
    ap.add_argument('--drt-randomizations', type=int, default=20)
    ap.add_argument('--gene-drt-randomizations', type=int, default=10)
    ap.add_argument('--gene-ltt', action='store_true', help='基因级 LTT (默认关, 耗时)')
    ap.add_argument('--outlier-threshold', type=float, default=3.0,
                    help='学生化残差剔除阈值 |t| (默认 3.0, 0=关闭)')
    ap.add_argument('--skip-genome', action='store_true')
    ap.add_argument('--skip-genes', action='store_true')
    a = ap.parse_args()
    log = LogCollector()
    r = run_clock_analysis(
        alignment=a.alignment, tree=a.tree, dates_csv=a.dates, out_dir=a.out_dir,
        meta_csv=a.meta, gbk=a.gbk, seq_len=a.seq_len, threads=a.threads,
        iqtree_bin=a.iqtree, rscript=a.rscript,
        drt_randomizations=a.drt_randomizations,
        gene_drt_randomizations=a.gene_drt_randomizations,
        gene_ltt=a.gene_ltt, skip_genome=a.skip_genome, skip_genes=a.skip_genes,
        outlier_threshold=a.outlier_threshold, log=log)
    if not r['success']:
        sys.exit(1)


if __name__ == '__main__':
    main()
