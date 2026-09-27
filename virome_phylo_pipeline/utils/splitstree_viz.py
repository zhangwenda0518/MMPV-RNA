#!/usr/bin/env python3
"""
splitstree_viz.py -- SplitsTree6 NeighborNet nexus 的等角布局 (equal-angle) 出版级绘制
解析 nexus 的 TAXA/CYCLE/SPLITS 块, 用 zonotope (平行四边形 zone) 渲染分裂网络,
视觉等价 SplitsTree 的 equal-angle 视图。

用法: python3 splitstree_viz.py -i neighbornet.nexus -o neighbornet.pdf [--png] [--label-every N]
输出: PDF (矢量, 期刊首选) + 可选 PNG 300dpi; 白底, 全英文标签
"""

import argparse
import math
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_nexus(path):
    text = open(path, errors='ignore').read()
    # TAXLABELS
    labels = re.findall(r"\[\d+\]\s+'([^']+)'", text)
    if not labels:
        m = re.search(r'TAXLABELS(.*?);', text, re.S)
        labels = m.group(1).split() if m else []
    # CYCLE
    m = re.search(r'CYCLE\s+([\d\s]+?);', text, re.S)
    if not m:
        raise ValueError('no CYCLE block (not a circular split system?)')
    cycle = [int(x) for x in m.group(1).split()]
    # SPLITS MATRIX: lines like "[k, size=s] \t w \t t1 t2 ...,"
    m = re.search(r'^BEGIN SPLITS;(.*?)^END;', text, re.S | re.M)
    if not m:
        raise ValueError('no SPLITS block')
    splits = []  # (weight, frozenset taxa-ids)
    for line in m.group(1).splitlines():
        line = line.strip()
        mm = re.match(r'\[\d+,\s*size=\d+\]\s+([\d.eE+-]+)\s+(.*?),?\s*$', line)
        if not mm:
            continue
        w = float(mm.group(1))
        taxa = frozenset(int(x) for x in mm.group(2).split())
        splits.append((w, taxa))
    return labels, cycle, splits


def equal_angle_zones(labels, cycle, splits):
    """每个 split 渲染为一个 zone (平行四边形, 四边均为网络边)。
    zone 半径 d = 包含它的 split 权重之和; zone 高 = 自身权重 w。"""
    n = len(cycle)
    pos = {t: k for k, t in enumerate(cycle)}  # taxon-id -> cycle position
    sector = 2 * math.pi / n

    zones = []
    for w, taxa in splits:
        if w <= 0 or len(taxa) in (0, n):
            continue
        ps = sorted(pos[t] for t in taxa)
        # 圆形 split 系统: 某一侧在 cycle 上连续。取连续的一侧。
        # 判断: 若补集连续或本集连续
        def contiguous(positions):
            if not positions:
                return None
            positions = sorted(positions)
            # 尝试每个断点作为绕回点
            gaps = [(positions[(i + 1) % len(positions)] - positions[i]) % n
                    for i in range(len(positions))]
            if all(g == 1 for g in gaps) and len(positions) < n:
                return positions
            big = max(range(len(positions)),
                      key=lambda i: gaps[i])
            if gaps[big] > 1 and all(g == 1 for i, g in enumerate(gaps) if i != big):
                # 断点在 big 之后, 起点为 positions[big+1]
                start = big + 1
                rotated = positions[start:] + positions[:start]
                return rotated
            return None

        side = contiguous(ps)
        if side is None:
            comp = contiguous(sorted(set(range(1, n + 1)) - set(ps)))
            side = comp
        if side is None:
            continue
        # 连续区间 [side[0], side[-1]] (可能绕回)
        a, b = side[0], side[-1]
        span = (b - a) % n + 1
        zones.append({'taxa': frozenset(taxa), 'set_side': frozenset(taxa),
                      'a': a, 'b': b, 'span': span, 'w': w})
    # containment depth
    for z in zones:
        z['d'] = sum(o['w'] for o in zones
                     if o is not z and z['taxa'] < o['taxa'])
    return zones, n, sector


def render(nexus, out_pdf, out_png=None, label_every=1, title=None):
    labels, cycle, splits = parse_nexus(nexus)
    zones, n, sector = equal_angle_zones(labels, cycle, splits)
    total_w = max((z['d'] + z['w'] for z in zones), default=1.0)
    R = 1.0
    fig, ax = plt.subplots(figsize=(9, 9), facecolor='white')
    ax.set_facecolor('white')

    def pt(angle, r):
        return (r * math.cos(angle), r * math.sin(angle))

    # zone 平行四边形: 内弦(d) / 外弦(d+w) / 两条径向边
    for z in zones:
        th1 = z['a'] * sector            # 区间起始射线
        th2 = (z['b'] + 1) * sector      # 区间结束射线
        d_in = z['d'] / total_w * R
        d_out = (z['d'] + z['w']) / total_w * R
        p1, p2 = pt(th1, d_in), pt(th2, d_in)
        p3, p4 = pt(th2, d_out), pt(th1, d_out)
        quad = [p1, p2, p3, p4, p1]
        xs, ys = zip(*quad)
        ax.plot(xs, ys, color='#333333', lw=0.6, zorder=2, solid_capstyle='round')

    # 分类单元: 放在外圈, 连线到最外层包含它的 zone 外弦中点方向
    for k, t in enumerate(cycle):
        th = (k + 0.5) * sector
        outer = max((z for z in zones if t in z['taxa']),
                    key=lambda z: z['d'] + z['w'], default=None)
        r0 = (outer['d'] + outer['w']) / total_w * R if outer else 0.0
        p1, p2 = pt(th, r0), pt(th, R * 1.03)
        ax.plot(*zip(p1, p2), color='#888888', lw=0.5, zorder=1)
        lx, ly = pt(th, R * 1.10)
        lab = labels[t - 1] if t - 1 < len(labels) else str(t)
        ha = 'left' if math.cos(th) > 0.1 else ('right' if math.cos(th) < -0.1 else 'center')
        if k % label_every == 0:
            ax.text(lx, ly, lab, ha=ha, va='center', fontsize=5.5, color='black')

    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect('equal')
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=10)
    nw = len([z for z in zones])
    fit = None
    import re as _re
    m = _re.search(r'fit=([\d.]+)', open(nexus, errors='ignore').read())
    if m:
        fit = float(m.group(1))
    cap = 'Neighbor-Net (SplitsTree6), n=%d taxa, %d splits' % (n, nw)
    if fit is not None:
        cap += ', fit=%.1f%%' % fit
    fig.text(0.5, 0.03, cap, ha='center', fontsize=9)
    fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    if out_png:
        fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('saved:', out_pdf, '(%d zones)' % len(zones))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', required=True)
    ap.add_argument('-o', '--pdf', required=True)
    ap.add_argument('--png', default=None)
    ap.add_argument('--label-every', type=int, default=1,
                    help='标注每 N 个 taxa (大图降密度)')
    ap.add_argument('--title', default=None)
    a = ap.parse_args()
    render(a.input, a.pdf, a.png, a.label_every, a.title)
