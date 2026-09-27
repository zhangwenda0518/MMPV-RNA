#!/usr/bin/env python3
"""
stage_plots.py — 补齐各 stage 的 summary 图 (SCI 风格: 白底/Okabe-Ito/300dpi/pdf+png)

覆盖四个空白:
  ① align   QC 图   : 序列长度分布 + gap 比例 (比对质量一目了然)
  ② tree    可视化  : matplotlib 环形树 (无 ETE 依赖, 全 Python)
  ③ popgen  图      : π/θ/Tajima's D/Fst 柱状图 (从 popgen_report.txt 解析)
  ④ gene_dating 汇总: 各基因 tMRCA (95% HPD) 森林图 (从各基因目录扫 phylo_dates_*.csv)

也可独立运行: python utils/stage_plots.py --kind align_qc --fasta aln.fa --outdir d/
"""
import os
import re
import sys
from pathlib import Path

UTILS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(UTILS_DIR))
sys.path.insert(0, str(UTILS_DIR.parent))

from virphy_bridge import LogCollector  # noqa: E402

# Okabe-Ito
OI = {'blue': '#0072B2', 'orange': '#D55E00', 'green': '#009E73',
      'sky': '#56B4E9', 'yellow': '#E69F00', 'purple': '#CC79A7',
      'grey': '#999999'}


def _finish(fig, out_pdf, out_png, dpi=300):
    fig.savefig(out_pdf, dpi=dpi)
    fig.savefig(out_png, dpi=dpi)
    import matplotlib.pyplot as plt
    plt.close(fig)


def _style_ax(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=9)


# ═══ ① align QC ═══

def plot_align_qc(fasta_file, output_dir, log=None):
    """序列长度分布 + 每序列 gap 比例; 返回 png 路径或 None"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    names, lengths, gaps = [], [], []
    cur_n, cur_s = None, []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                if cur_n is not None:
                    names.append(cur_n); lengths.append(len(cur_s))
                    gaps.append(cur_s.count('-') / max(1, len(cur_s)))
                cur_n = line[1:].split()[0]; cur_s = []
            else:
                cur_s.append(line.strip())
    if cur_n is not None:
        names.append(cur_n); lengths.append(len(cur_s))
        gaps.append(cur_s.count('-') / max(1, len(cur_s)))
    if not names:
        return None
    os.makedirs(output_dir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4), dpi=300)
    axes[0].hist(lengths, bins=30, color=OI['blue'], edgecolor='white', linewidth=0.4)
    axes[0].set_xlabel('Aligned sequence length (bp)', fontsize=10)
    axes[0].set_ylabel('Count', fontsize=10)
    axes[0].set_title(f'n = {len(names)}', fontsize=9)
    _style_ax(axes[0])
    axes[1].hist([g * 100 for g in gaps], bins=30, color=OI['orange'],
                 edgecolor='white', linewidth=0.4)
    axes[1].set_xlabel('Gap fraction (%)', fontsize=10)
    axes[1].set_ylabel('Count', fontsize=10)
    _style_ax(axes[1])
    fig.tight_layout()
    _finish(fig, os.path.join(output_dir, 'align_qc.pdf'),
            os.path.join(output_dir, 'align_qc.png'))
    return os.path.join(output_dir, 'align_qc.png')


# ═══ ② tree 可视化 (纯 Python newick → 环形树) ═══

def _parse_newick(s):
    s = s.strip().rstrip(';')
    pos = [0]

    def node():
        children = []
        if s[pos[0]] == '(':
            pos[0] += 1
            while True:
                children.append(node())
                if s[pos[0]] == ',':
                    pos[0] += 1
                else:
                    break
            pos[0] += 1  # ')'
        # label
        m = re.match(r'[^,:()\[\]]*', s[pos[0]:])
        label = m.group(0)
        pos[0] += len(label)
        length = 0.0
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            m = re.match(r'-?[0-9.eE+-]+', s[pos[0]:])
            if m:
                length = float(m.group(0))
                pos[0] += len(m.group(0))
        return {'label': label, 'length': length, 'children': children}

    return node()


def plot_tree_circular(tree_file, output_dir, log=None):
    """环形系统发育树 (matplotlib, 无 ETE 依赖)"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None
    try:
        txt = open(tree_file).read()
        # 去掉 bootstrap 支持值标签保留分支长: (A:0.1,(B:0.2,C:0.3)95:0.4);
        txt = re.sub(r'\)\s*\d+(\.\d+)?\s*([,):])', r')\2', txt)
        root = _parse_newick(txt)
    except Exception as e:
        (log or LogCollector()).emit(f"[tree_plot] 解析失败: {e}")
        return None

    # 收集 tips
    tips = []

    def collect(n, depth):
        if not n['children']:
            tips.append({'label': n['label'], 'd': depth})
        for c in n['children']:
            collect(c, depth + c['length'])
    collect(root, 0.0)
    if not tips:
        return None
    tips.sort(key=lambda t: t['d'])
    n_tips = len(tips)

    os.makedirs(output_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=300,
                           subplot_kw={'projection': 'polar'})
    max_d = max(t['d'] for t in tips) or 1.0
    for i, t in enumerate(tips):
        theta = 2 * np.pi * i / n_tips
        ax.plot([0, theta], [0, t['d']], color=OI['grey'], lw=0.6)
        if t['label']:
            ax.text(theta, t['d'] + max_d * 0.03, t['label'],
                    fontsize=3.2, ha='center', va='center', rotation=0)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.spines['polar'].set_visible(False)
    ax.set_ylim(0, max_d * 1.18)
    fig.tight_layout(pad=0.4)
    _finish(fig, os.path.join(output_dir, 'tree_circular.pdf'),
            os.path.join(output_dir, 'tree_circular.png'))
    return os.path.join(output_dir, 'tree_circular.png')


# ═══ ③ popgen 图 ═══

def plot_popgen(report_txt, output_dir, log=None):
    """从 popgen_report.txt 解析 π/θ/Tajima's D/Fst 画柱状图"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    if not os.path.exists(report_txt):
        return None
    txt = open(report_txt, encoding='utf-8', errors='ignore').read()
    metrics = {}
    for key, pat in [('π', r"π\s*=\s*([0-9.eE+-]+)"),
                     ('θw', r"θ[wW]?\s*=\s*([0-9.eE+-]+)"),
                     ("Tajima's D", r"Tajima'?s?\s*D\s*=\s*([0-9.eE+-]+)"),
                     ('Fst', r"F[Ss][Tt]\s*=\s*([0-9.eE+-]+)")]:
        m = re.search(pat, txt)
        if m:
            metrics[key] = float(m.group(1))
    if not metrics:
        return None
    os.makedirs(output_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(4.6, 3.4), dpi=300)
    keys = list(metrics.keys())
    vals = [metrics[k] for k in keys]
    colors = [OI['blue'], OI['sky'], OI['orange'], OI['green']][:len(keys)]
    bars = ax.bar(range(len(keys)), vals, color=colors, width=0.55,
                  edgecolor='white', linewidth=0.4)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2,
                v + (0.02 * max(abs(x) for x in vals) or 0.02) * (1 if v >= 0 else -1),
                f'{v:.4g}', ha='center', va='bottom' if v >= 0 else 'top', fontsize=8)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(keys, fontsize=9)
    ax.axhline(0, color='#555555', lw=0.6)
    ax.set_ylabel('Value', fontsize=10)
    _style_ax(ax)
    fig.tight_layout()
    _finish(fig, os.path.join(output_dir, 'popgen_metrics.pdf'),
            os.path.join(output_dir, 'popgen_metrics.png'))
    return os.path.join(output_dir, 'popgen_metrics.png')


# ═══ ④ gene_dating 汇总 (tMRCA 森林图) ═══

def plot_gene_dating_summary(gene_dating_dir, output_dir, log=None):
    """扫描各基因子目录 phylo_dates_*.csv, 提取 tMRCA/HPD 画森林图 + 汇总 TSV"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        return None, None
    if not os.path.isdir(gene_dating_dir):
        return None, None
    rows = []
    for sub in sorted(Path(gene_dating_dir).iterdir()):
        if not sub.is_dir():
            continue
        for csv_f in sorted(sub.glob('phylo_dates_*.csv')):
            try:
                df = pd.read_csv(csv_f)
                for col_cand in ['tmra', 'tmrca', 'TMRCA', 't_mrca']:
                    if col_cand in df.columns:
                        vals = pd.to_numeric(df[col_cand], errors='coerce').dropna()
                        if len(vals):
                            rows.append({'gene': sub.name, 'file': csv_f.name,
                                         'tmrca_mean': vals.mean(),
                                         'tmrca_lower': vals.quantile(0.025),
                                         'tmrca_upper': vals.quantile(0.975)})
                        break
            except Exception:
                continue
    if not rows:
        return None, None
    os.makedirs(output_dir, exist_ok=True)
    # 汇总 TSV
    tsv = os.path.join(output_dir, 'gene_tmrca_summary.tsv')
    with open(tsv, 'w') as f:
        f.write('gene\ttmrca_mean\ttmrca_lower95\ttmrca_upper95\n')
        for r in rows:
            f.write(f"{r['gene']}\t{r['tmrca_mean']:.4f}\t"
                    f"{r['tmrca_lower']:.4f}\t{r['tmrca_upper']:.4f}\n")
    # 森林图
    rows.sort(key=lambda r: r['tmrca_mean'])
    fig, ax = plt.subplots(figsize=(5.6, 1.2 + 0.5 * len(rows)), dpi=300)
    for i, r in enumerate(rows):
        ax.plot([r['tmrca_lower'], r['tmrca_upper']], [i, i],
                color=OI['blue'], lw=2.2, solid_capstyle='round')
        ax.plot(r['tmrca_mean'], i, 'o', color=OI['orange'],
                ms=6, mec='white', mew=0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r['gene'] for r in rows], fontsize=9)
    ax.set_xlabel('tMRCA (years before present)', fontsize=10)
    _style_ax(ax)
    ax.invert_yaxis()
    fig.tight_layout()
    _finish(fig, os.path.join(output_dir, 'gene_tmrca_forest.pdf'),
            os.path.join(output_dir, 'gene_tmrca_forest.png'))
    return tsv, os.path.join(output_dir, 'gene_tmrca_forest.png')


# ═══ CLI ═══

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--kind', required=True,
                    choices=['align_qc', 'tree', 'popgen', 'gene_dating'])
    ap.add_argument('--input', required=True)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()
    if a.kind == 'align_qc':
        print(plot_align_qc(a.input, a.outdir))
    elif a.kind == 'tree':
        print(plot_tree_circular(a.input, a.outdir))
    elif a.kind == 'popgen':
        print(plot_popgen(a.input, a.outdir))
    elif a.kind == 'gene_dating':
        t, p = plot_gene_dating_summary(a.input, a.outdir)
        print(t, p)
