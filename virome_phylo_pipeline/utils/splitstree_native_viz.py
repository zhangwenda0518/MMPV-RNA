#!/usr/bin/env python3
"""
splitstree_native_viz.py -- 用 **SplitsTree6 自己的布局** 出图 (不再自绘 EqualAngle)

背景 (2026-10-12 第 8 轮, 用户: "你这根本不是 splitstree 的结果图呀, 画的太丑了,
我们有 splitstree, 你打开看看"):
  旧 `splitstree_viz.py` 用 Python 重新实现了 Huson & Bryant 的 equal-angle,
  形状与 SplitsTree 本体**不一致** (自绘出来是"星芒辐射扇", SplitsTree 出来是
  宽扁网状)。本轮改为直接消费 SplitsTree 本体的布局坐标:

    ~/biosoft/splitstree6-tools/tools/webserver
      GET /splitstree/draw_distances
          ?matrix=<nexus DISTANCES 块>&algorithm=nnet&layout=radial
          &output=coordinates&width=<W>&height=<H>
      → 文本坐标 (第 1 行 <nNodes>\\t<nEdges>;
                  节点 <id>\\t<x>\\t<y>[\\t<taxonLabel>];
                  边   <id>\\t<src>\\t<tgt>\\t<weight>\\t<splitId>)

  ⚠ 端点是 `/draw_distances` (下划线)。`/drawDistances` 会被 com.sun.net.httpserver
    的**前缀匹配**吃掉, 落到 `/draw` → 报 `Unknown parameter: 'matrix'`。

用法:
  python3 splitstree_native_viz.py --coords nnet_radial.coords -o neighbornet.pdf \
      [--png neighbornet.png] [--nexus neighbornet.nexus] \
      [--meta sample_metadata.csv --group-col location] [--label-rings 3] \
      [--title "Cytorhabdovirus sp. lycii (OR489165.1)"]

等长 + 等粗 (SplitsTree `SplitsTopology` 观感):
  默认 `--unweight auto` + `--lw-preset auto` 即为**等长 + 等粗** —— 输入若还是
  加权坐标, 脚本会先内联调用 `splitstree_unweight.py` 把权重当 1 重放 (等长),
  再把线宽切成 `flat` (等粗)。要复现**旧加权图**给 `--unweight off --lw-preset legacy`。

输出: PDF (矢量) + 可选 PNG 300dpi; 白底, 全英文。
"""

import argparse
import math
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt

try:                                     # 样式真源 (与本目录其它模块共享同一份)
    from utils.sci_palette import sci_rcparams, OKABE_ITO, COL_W_1, COL_W_2
    # 2026-10-16: 图件双格式 (PNG+PDF 同主名)
    from utils.plot_io import save_both
except ImportError:                      # 允许以脚本形式直接跑
    from sci_palette import sci_rcparams, OKABE_ITO, COL_W_1, COL_W_2
    from plot_io import save_both

try:
    from splitstree_viz import (_setup_font, load_group_map, _match_groups,
                               _norm_label, _short_group)
except ImportError:
    def _setup_font():
        return None

    def load_group_map(*a, **k):
        return {}, None

    def _match_groups(labels, groups):
        return {}

    def _norm_label(s):
        return str(s).strip().lower()

    def _short_group(grp):
        return str(grp).split(',')[-1].strip() or str(grp)

# ── 组→色固定映射 (2026-10-14 修复, 全模块唯一定义) ─────────────────────────
# 病根: 本模块原先有**四处独立取色** —— 弧线/弧名按「组大小降序」取 OKABE,
# tip 点/标签/图例按「组名升序」取 OKABE —— 两组排序在宿主数据上不一致
# (barbarum 103 / ruthenicum 14 / chinense 3 → 弧序=b,r,c / 名序=b,c,r),
# 于是面板 D 里 ruthenicum 拿了 chinense 的浅蓝、chinense 拿了 ruthenicum 的绿,
# 与面板 C (sci_palette 固定映射) 直接对撞 —— 同图同物不同色。
# 现钉死: 已知宿主组 (大小写不敏感) → 全线约定色; 其余组 (地点等) → viridis
# 冷色系按名升序采样 (与 classic 树地点点同族, 不与宿主暖色撞)。
_HOST_FIXED = {
    'lycium barbarum': '#E69F00',    # OKABE_ITO[0]
    'lycium chinense': '#56B4E9',    # OKABE_ITO[1]
    'lycium ruthenicum': '#009E73',  # OKABE_ITO[2]
}


def _canonical_group_color(gnames):
    """组名可迭代对象 → {组名: 颜色}。

    已知宿主组用 _HOST_FIXED 固定色; 其余组按名升序从 viridis 采样,
    且避开已占用的固定色 (宿主/地点混出时也不撞)。
    """
    known, unknown, used = {}, [], set()
    for g in gnames:
        c = _HOST_FIXED.get(str(g).strip().lower())
        if c:
            known[g] = c
            used.add(c)
        else:
            unknown.append(g)
    for i, g in enumerate(sorted(unknown)):
        frac = i / max(1, len(unknown) - 1) if len(unknown) > 1 else 0.5
        col = matplotlib.colormaps['viridis'](frac)
        if matplotlib.colors.to_hex(col).upper() in {u.upper() for u in used}:
            col = matplotlib.colormaps['viridis'](
                min(1.0, frac + 0.618 / max(1, len(unknown))))
        known[g] = col
    return known

C_EDGE = '#1a1a1a'          # 网络线 (与文献 NeighborNet 一致: 统一深色, 权重只控线宽)
C_TIP = '#111111'
C_LEADER = '#b8b8b8'
_LBL_STROKE = 2.2

# ── 线宽映射 (split 权重 → 线段粗细) ─────────────────────────────────────────
# ⚠ 2026-10-12 实测: 本数据 829 条边的权重**极度右偏** (p50=0.00026, max=0.01461,
#   跨 56 倍)。旧的线性映射 `0.22+0.85*(w/wmax)` 在这种分布下 → 中位线宽只有
#   0.23 pt, 87 % 的边 < 0.3 pt, 仅 4/829 条 ≥ 0.5 pt → **肉眼完全看不出粗细差异**
#   (虽然代码"有"权重映射, 实际等于没展示枝长)。
#
# 解法 = 先做**幂次压缩**再线性铺到 [lw_min, lw_max]:
#     lw = lw_min + (lw_max - lw_min) * (w / w_ref) ** lw_gamma
#   `lw_gamma < 1` 把右偏分布压平, 小权重也被抬离地板。
#   参考尺度 `w_ref` 取 **p95 分位**而非 max: max 是单条离群极值, 用它做分母会把
#   其余 99 % 的边全压到接近 0。p95 让"最粗那一档"由"显著粗的边"占据。
_LW_MIN = 0.35              # 最细边的线宽 (pt)
_LW_MAX = 1.55              # 最粗边的线宽 (pt)
_LW_GAMMA = 0.45            # 幂次压缩指数 (<1 抬升小权重)
_LW_REF_PCT = 95            # 归一化分母用的分位 (%)

# 预设 (供 CLI `--lw-preset` 选用; 空字符串=用模块默认三参)
LW_PRESETS = {
    'legacy':  (0.22, 1.07, 1.00, 100),   # 旧行为 (线性/max 归一) —— 只为复现旧图
    'std':     (_LW_MIN, _LW_MAX, _LW_GAMMA, _LW_REF_PCT),
    'bold':    (0.50, 1.90, 0.50, 95),    # 更粗, 适合缩到半栏
    'subtle':  (0.30, 1.15, 0.55, 95),    # 更含蓄
    'flat':    (0.55, 0.55, 1.00, 100),   # **等粗** (所有权重 → 同一线宽)
}

# ── 等长布局 (unweighted / topology layout) ──────────────────────────────────
# 用户 2026-10-12: "展示线宽就别有枝长了 ... 取消枝长, 等长, 我们给的图不就是嘛"。
#
# 背景: 参考图 (SplitsTree GUI 截图) 是**等长**的 —— 所有边长度一致 → 视觉上是
#   规则网格; 而 SplitsTree 加权布局把 split 权重烙进几何 → 边长跨 **344 倍**
#   (p50=19.5 / max=1111.8) → 放射状长枝。两者观感差得远。
#
# 解法: 同一套角度, 把每条边的权重**当 1 重放** DFS (源码 EqualAngle:
#   `p_w = p_v + weight_e * (cos, sin)`, `weight_e = useWeights ? getWeight : 1`)。
#   实现在 `_unweight.py` (见其文件头, 含"原权重重放 == 原坐标"自检)。
#
# 本模块**既产出也接收** unweighted 坐标 (2026-10-13 起 `unweight_coords()` 内联),
# 渲染判据:
#   ① 边长 max/min ≈ 1.0x (等长);  ② 线宽走 `flat` 预设 (等粗)。
_UNW_TOL = 1.35        # 判"是否已是等长布局"的 max/min 阈值

# ── 外圈分组弧线 (2026-10-12 第 9 轮) ────────────────────────────────────────
# 目标样式 = 文献里最常见的 NeighborNet 配图: 每组 tips 外面套一条**同色平滑弧线**,
# 组名用同色粗体挂在弧外。⚠ SplitsTree6 本体**没有**这个功能 (源码实测见
# `reports/_phylo_fix_20261011/FIX_REPORT_20261012.md` §⑬): 它只有
# `SplitsDiagramType.Outline`(整张网络**一个**轮廓, 单一填色)、`TraitsFormat`(节点饼图)、
# `view/format/taxmark`(每个 taxon 一个小标记)。所以弧线必须我们自己画。
_ARC_GAP = 0.30        # 同组相邻 tip 的角间隔 > 此值 (rad) → 断开成新的一段弧
_ARC_MIN_RUN = 3       # 一段弧至少含几个 tip, 不足则丢弃 (组间高度交错时抑制碎片)
_ARC_EXT = 0.10        # 弧两端各向外延伸 (rad)
_ARC_WIN = 0.30        # 外缘包络的角窗宽 (rad): 窗内 tips 的最大半径 = 该方向的外缘
_ARC_PAD = 0.045       # 弧相对网络外缘的径向外推量 (占最大半径比例)
_ARC_LW = 1.7          # 弧线宽
_ARC_LABEL_PAD = 0.085  # 组名相对弧的径向外推量
_ARC_LABEL_MIN_D = 0.16  # 组名之间的最小间距 (占最大半径比例), 冲突则继续外推


# ──────────────────────────────────────────────────────────────────────────────
# 读 SplitsTree 原生坐标
# ──────────────────────────────────────────────────────────────────────────────
def load_coords(path):
    """SplitsTree6 webserver /draw_distances?output=coordinates 文本 → (nodes, labels, edges).

    nodes : {nodeId: (x, y)}
    labels: {nodeId: taxonLabel}   (仅带标签的节点)
    edges : [(srcId, tgtId, weight, splitId)]
    """
    with open(path, encoding='utf-8', errors='replace') as f:
        lines = f.read().splitlines()
    if not lines:
        raise ValueError('empty coords file: %s' % path)
    n_n, n_e = (int(v) for v in lines[0].split()[:2])
    nodes, labels = {}, {}
    for i in range(1, 1 + n_n):
        p = lines[i].split('\t')
        nodes[int(p[0])] = (float(p[1]), float(p[2]))
        if len(p) > 3 and p[3].strip():
            labels[int(p[0])] = p[3].strip()
    edges = []
    for i in range(1 + n_n, 1 + n_n + n_e):
        p = lines[i].split('\t')
        if len(p) < 5:
            continue
        edges.append((int(p[1]), int(p[2]), float(p[3]), int(p[4])))
    return nodes, labels, edges


def coords_bbox(nodes):
    xs = [v[0] for v in nodes.values()]
    ys = [v[1] for v in nodes.values()]
    return min(xs), max(xs), min(ys), max(ys)


# ──────────────────────────────────────────────────────────────────────────────
# 外圈分组弧线 (文献 NeighborNet 配图的标准画法)
# ──────────────────────────────────────────────────────────────────────────────
def _angular_runs(theta_sorted, gap):
    """把按角度排好序的 tip 切成若干**连续段** (环形)。

    返回 [[start, end), ...] —— 下标是"展开数组"上的区间 (start 可能 > end 时
    表示跨 0 环绕, 调用方用 `i % n` 取值)。

    为什么必须分段: 本数据 Yinchuan 的 69 个 tip 在角向被其他组切成 19 段
    (第 6 轮实测)。整组只画一条"最短包弧"必然横扫大半个圆, 把不相干的 tip 圈进来。
    """
    n = len(theta_sorted)
    if n == 0:
        return []
    if n == 1:
        return [(0, 1)]
    two_pi = 2 * math.pi
    gaps = [i for i in range(n)
            if (theta_sorted[(i + 1) % n] - theta_sorted[i]) % two_pi > gap]
    if not gaps:
        return [(0, n)]
    runs = []
    for k, g in enumerate(gaps):
        s = (g + 1) % n
        e = gaps[(k + 1) % len(gaps)] + 1
        if e <= s:
            e += n
        runs.append((s, e))
    return runs


def draw_group_arcs(ax, tips, labels, groups, fontsize=8.0, min_sep=None,
                    obstacles=None):
    """在每个分组的 tips 外侧画一条**同色平滑弧线**, 组名用同色粗体挂在弧外。

    tips  : {nodeId: (x, y)}  —— 只含带标签的节点
    labels: {nodeId: taxonLabel}
    groups: {元数据名: 组名}  —— 内部走 `_match_groups` 归一到 nexus 标签口径
    obstacles: [(cx, cy, hx, hy)]  —— **已有的文字**占位框 (数据坐标, 半宽/半高)。
       整页版有 145 个 tip 标签, 组名必须避开它们, 否则核心标注被压
       (实测整页版 24 对重叠全是组名 vs tip 标签)。面板版无 tip 标签 → 传 None。

    几何 (极坐标, 原点取 tips 的包围框中心):
      ① 外缘包络 R(θ) = 角窗 `_ARC_WIN` 内所有 tip 的**最大**半径 —— 这样弧线是
         贴着网络外轮廓走的"衣领", 而不是一条等半径圆 (本数据网络是斜向菱形,
         等半径圆会切进网面);
      ② 同组 tips 按角度分段 (见 `_angular_runs`), 每段两端外延 `_ARC_EXT`;
      ③ 弧半径 = max(包络, 本段最大半径) + `_ARC_PAD`, 再沿角向做滑动平均平滑;
      ④ 组名挂在最大那一段的角中点、弧外 `_ARC_LABEL_PAD` 处, 相互冲突则径向外推。

    返回 (legend_items, pts, pad_needed):
      legend_items = [(组名, 颜色, tip 数)] (按 tip 数降序, 颜色取 OKABE_ITO)
      pts          = 所有弧线/组名用到的点 (供调用方扩画布)
      pad_needed   = 需要的最小边距 (相对 w / h 的比例)
    """
    if not tips or not groups:
        return [], [], 0.0
    lab2node = {}
    for i, lab in labels.items():
        lab2node.setdefault(lab, i)
    gmap = _match_groups(labels.values(), groups) or {}
    if not gmap:
        return [], [], 0.0

    by_group = {}
    for lab, grp in gmap.items():
        i = lab2node.get(lab)
        if i is not None and i in tips:
            by_group.setdefault(_short_group(grp), []).append(i)
    if not by_group:
        return [], [], 0.0

    cx = 0.5 * (min(v[0] for v in tips.values()) + max(v[0] for v in tips.values()))
    cy = 0.5 * (min(v[1] for v in tips.values()) + max(v[1] for v in tips.values()))
    pol = {i: (math.atan2(tips[i][1] - cy, tips[i][0] - cx),
               math.hypot(tips[i][0] - cx, tips[i][1] - cy)) for i in tips}
    r_all = [p[1] for p in pol.values()]
    r_max = max(r_all) if r_all else 1.0
    _all = list(pol.values())

    def _env_all(q):
        """**全部** tips 在角窗内的最大半径 = 该方向上网络的外缘。

        组名用它定位 —— 只贴本组弧线会把名字压在密集簇内部 (本数据 6 个组的
        "最大段"都在中心密集区), 放到全局外缘之外才与参考图一致。
        """
        best = 0.0
        for th, rr in _all:
            d = abs((th - q + math.pi) % (2 * math.pi) - math.pi)
            if d <= _ARC_WIN and rr > best:
                best = rr
        return best

    items = sorted(by_group.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    gcolor = _canonical_group_color(by_group.keys())
    pts, labels_out = [], []
    for gi, (grp, ids) in enumerate(items):
        color = gcolor[grp]
        ids_sorted = sorted(ids, key=lambda i: pol[i][0])
        ths = [pol[i][0] for i in ids_sorted]
        runs = [r for r in _angular_runs(ths, _ARC_GAP) if (r[1] - r[0]) >= _ARC_MIN_RUN]

        def _seg_rmax(s, e):
            return max(pol[ids_sorted[i % len(ths)]][1] for i in range(s, e))

        def _seg_qm(s, e):
            """段角中点。⚠ 跨 0 环绕的段必须先把 e-1 的角度 +2π 展开再取中点 ——
            否则中点会落到**对侧** (实测 "China"/"Xining" 两个组名被放到半径 142
            的圆心附近, 完全不可见)。弧线绘制那段一直是对的, 标签这段漏过。"""
            _q0 = ths[s % len(ths)]
            _q1 = ths[(e - 1) % len(ths)]
            if e - s > len(ths):
                _q1 += 2 * math.pi
            return 0.5 * (_q0 + _q1)

        # 组名锚点段 = **相对暴露度最高**的那一段:
        #     score = 段最大半径 / 该角度的全局外缘
        #   ⚠ 两个更朴素的选法都不行 (实测):
        #     ① 选"tip 数最多的段" → 本数据各组的最大段几乎都在中心密集区、角中点
        #        几乎同向 → 组名只能沿同一条射线外推, 整串被推到右缘外并**被画布切掉**;
        #     ② 选"绝对最远的段" → 会把名字推到**别的组**伸得更远的方向上 (实测
        #        Neimenggu 的段在 143°/r=690, 而 Yinchuan 的 tip 在同方向伸到 1052,
        #        名字被顶到 r≈1189、离自己的弧 500 单位, 完全脱节)。
        #     用相对暴露度则名字落在"本组弧最外露"的方向上, 与参考图一致。
        def _score(r):
            _q = _seg_qm(r[0], r[1])
            _env = _env_all(_q)
            _ro = _seg_rmax(r[0], r[1])
            return ((_ro / _env) if _env > 0 else 1.0, _ro)

        anchor = max(runs, key=_score) if runs else None
        for (s, e) in runs:
            q0 = ths[s % len(ths)]
            q1 = ths[(e - 1) % len(ths)]
            if e - s > len(ths):
                q1 += 2 * math.pi
            a0, a1 = q0 - _ARC_EXT, q1 + _ARC_EXT
            # 包络只用**本段自己的 tips** (不掺入其它组/远处的 tip):
            #   本数据网络是"中心密集 + 长悬枝", 用全局包络会被远端悬枝把弧线拉飞
            #   (实测弧线飞到图外, 且横穿网面)。
            seg = [(pol[ids_sorted[i % len(ths)]][0], pol[ids_sorted[i % len(ths)]][1])
                   for i in range(s, e)]

            def _seg_env(q, _seg=seg):
                best = 0.0
                for th, rr in _seg:
                    d = abs((th - q + math.pi) % (2 * math.pi) - math.pi)
                    if d <= _ARC_WIN and rr > best:
                        best = rr
                if best <= 0.0:                      # 角窗内没有本段 tip → 取最近者
                    best = min(_seg, key=lambda t: abs(
                        (t[0] - q + math.pi) % (2 * math.pi) - math.pi))[1]
                return best

            n_samp = max(24, int((a1 - a0) / 0.02))
            qs = [a0 + (a1 - a0) * k / (n_samp - 1) for k in range(n_samp)]
            rs = [_seg_env(q) * (1.0 + _ARC_PAD) for q in qs]
            # 角向滑动平均 → 平滑弧 (窗口 5)
            sm = [sum(rs[max(0, k - 2):k + 3]) / len(rs[max(0, k - 2):k + 3])
                  for k in range(len(rs))]
            xy = [(r * math.cos(q), r * math.sin(q)) for q, r in zip(qs, sm)]
            ax.plot([p[0] + cx for p in xy], [p[1] + cy for p in xy], '-',
                    color=color, lw=_ARC_LW, solid_capstyle='round', zorder=2.4)
            pts.extend((p[0] + cx, p[1] + cy) for p in xy)
        if anchor is not None:
            s, e = anchor
            qm = _seg_qm(s, e)
            rm = _seg_rmax(s, e) * (1.0 + _ARC_PAD)
            # 至少推到该角度的**全局外缘**之外, 否则名字会落在密集簇内部
            rm = max(rm, _env_all(qm) * (1.0 + _ARC_PAD))
            labels_out.append([grp, color, qm, rm + _ARC_LABEL_PAD * r_max])

    # 组名去重叠: 半径大的先落位, 冲突就沿径向外推
    #   ⚠ 间距必须**随字号/画布缩放** —— 固定 `_ARC_LABEL_MIN_D * r_max` 在窄画布
    #     (面板版 4.58 in / 6.5 pt) 上不够, 实测 6 个组名里有 3 对互相压。
    #   ⚠ 判据必须是**逐轴**的: 两个轴对齐文本框不压的充要条件是
    #     `|dx| >= 字宽` 或 `|dy| >= 字高`。用"中心距 >= 对角线"**不等价**
    #     (实测残留 1 对 12x11 px)。
    if isinstance(min_sep, (tuple, list)):
        sx, sy = float(min_sep[0]), float(min_sep[1])
    else:
        sx = sy = float(min_sep or 0.0)
    sx = max(sx, _ARC_LABEL_MIN_D * r_max)
    sy = max(sy, 0.6 * _ARC_LABEL_MIN_D * r_max)
    # 候选位移按"**切向优先, 径向兜底**"排列 ——
    #   ⚠ 只做径向外推是错的: 本数据 3 个组 (Guangzhou/Wuwei/Beijing) 的最外侧段都在
    #     右侧密集区且几乎同向 (θ≈4~12°), 径向外推只能沿同一条射线排队, 实测整串组名
    #     被推到网络右侧的空白区、与自己的弧线脱节 (面板版还一度被右缘切掉)。
    #     先切向 (±θ) 能把它们沿弧线方向摊开, 仍贴着各自的弧。
    _cands = [(0.0, 0.0)]
    for _m in range(1, 17):
        for _s in (1.0, -1.0):
            _cands.append((_s * _m * 0.045, 0.0))
    for _m in range(1, 13):
        for _s in (1.0, -1.0):
            _cands.append((_s * _m * 0.045, _m * 0.045 * r_max))
    # ⚠⚠ 兜底 (2026-10-12 第 18 轮): 上面两圈候选是"切向 ±0.7 rad + 径向 12 步"。
    #   等长布局把网络压扁 → tip 标签环更密 → 组名可落空间被挤压; 实测 rings=12 时
    #   `Neimenggu` 的候选**全部**被 tip 标签挡死 → `for ... if _clear: break` 无
    #   `else` 分支 → **直接放置** → 残留 3 对 (rings=9 生产默认无此问题)。
    #   修法 = 候选用尽后再做一轮**更密 / 更远**的径向扫描, 宁可把组名外推得远一点,
    #   也不允许压字 (组名压字是硬判据, 外推只是观感代价)。
    for _m in range(1, 41):
        for _s in (1.0, -1.0):
            _cands.append((_s * (0.22 + _m * 0.075), _m * 0.09 * r_max))
            _cands.append((_s * _m * 0.075, 0.0))
    placed = []
    for grp, color, q, rr in sorted(labels_out, key=lambda t: -t[3]):
        lx = ly = 0.0

        def _clear(lx, ly):
            # 与已落位的组名互不压 (逐轴判据; 组名是 ha/va=center, 故中心即盒心)
            for px, py in placed:
                if abs(lx - px) < sx and abs(ly - py) < sy:
                    return False
            # 与**已有文字** (tip 标签) 互不压: 标准轴对齐盒相交判据
            for ox0, ox1, oy0, oy1 in (obstacles or ()):
                if (lx - 0.5 * sx < ox1 and lx + 0.5 * sx > ox0
                        and ly - 0.5 * sy < oy1 and ly + 0.5 * sy > oy0):
                    return False
            return True

        for _dq, _dr in _cands:
            lx = cx + (rr + _dr) * math.cos(q + _dq)
            ly = cy + (rr + _dr) * math.sin(q + _dq)
            if _clear(lx, ly):
                break
        placed.append((lx, ly))
        ax.text(lx, ly, grp, ha='center', va='center', fontsize=fontsize,
                color=color, fontweight='bold', zorder=6,
                path_effects=[pe.withStroke(linewidth=_LBL_STROKE * 0.8,
                                            foreground='white')])
        pts.append((lx, ly))

    legend = [(g, gcolor[g], len(v))
              for k, (g, v) in enumerate(items)]
    pad_needed = 0.0
    if pts:
        px = [p[0] for p in pts]
        py = [p[1] for p in pts]
        w = max(1e-9, max(tips[i][0] for i in tips) - min(tips[i][0] for i in tips))
        h = max(1e-9, max(tips[i][1] for i in tips) - min(tips[i][1] for i in tips))
        pad_needed = max((max(px) - max(tips[i][0] for i in tips)) / w,
                         (min(tips[i][0] for i in tips) - min(px)) / w,
                         (max(py) - max(tips[i][1] for i in tips)) / h,
                         (min(tips[i][1] for i in tips) - min(py)) / h)
    return legend, pts, max(0.0, pad_needed)


# ──────────────────────────────────────────────────────────────────────────────
# 外围引线标签布局 (等价 SplitsTree 的 RadialLabelLayout: 标签推到网络外圈 + 引线)
# ──────────────────────────────────────────────────────────────────────────────
def _border_point(cx, cy, x, y, X0, X1, Y0, Y1):
    """从 (cx,cy) 经 (x,y) 射出, 与矩形 (X0..X1, Y0..Y1) 的交点。"""
    dx, dy = x - cx, y - cy
    if dx == 0 and dy == 0:
        return cx, Y0
    t = []
    if dx > 0:
        t.append((X1 - cx) / dx)
    elif dx < 0:
        t.append((X0 - cx) / dx)
    if dy > 0:
        t.append((Y1 - cy) / dy)
    elif dy < 0:
        t.append((Y0 - cy) / dy)
    t = [v for v in t if v > 0]
    tt = min(t) if t else 1.0
    return cx + dx * tt, cy + dy * tt


def _rect_param(bx, by, X0, X1, Y0, Y1, per):
    """矩形边界上一点的周长参数 t ∈ [0, per)。"""
    w, h = X1 - X0, Y1 - Y0
    eps = 1e-9
    if abs(by - Y0) < eps:                       # 下边
        return bx - X0
    if abs(bx - X1) < eps:                       # 右边
        return w + (by - Y0)
    if abs(by - Y1) < eps:                       # 上边
        return w + h + (X1 - bx)
    return 2 * w + h + (Y1 - by)                 # 左边


def _rect_at(t, X0, X1, Y0, Y1, per):
    """周长参数 → 矩形边界点 + 该处文字的水平/垂直对齐。"""
    w, h = X1 - X0, Y1 - Y0
    t = t % per
    if t < w:
        return X0 + t, Y0, 'center', 'top'
    t -= w
    if t < h:
        return X1, Y0 + t, 'left', 'center'
    t -= h
    if t < w:
        return X1 - t, Y1, 'center', 'bottom'
    t -= w
    return X0, Y1 - t, 'right', 'center'


def _spread(vals, gap, per):
    """一维单调排布: 保证相邻间距 >= gap (周长是环形, 做首尾环绕修正)。"""
    v = list(vals)
    n = len(v)
    if n < 2:
        return v
    for _ in range(400):
        moved = False
        for j in range(1, n):
            d = v[j] - v[j - 1]
            if d < gap:
                shift = (gap - d) / 2.0
                v[j - 1] -= shift
                v[j] += shift
                moved = True
        # 环绕: 末元素与首元素 (跨 per)
        d = (v[0] + per) - v[-1]
        if d < gap:
            shift = (gap - d) / 2.0
            v[-1] -= shift
            v[0] += shift
            moved = True
        if not moved:
            break
    return v


def _spread_rect(pairs, per, X0, X1, Y0, Y1, gap_min):
    """矩形边界上的**逐对精确**排布 (比 `_spread` 的均匀 gap 更省空间)。

    ⚠ 实测 (2026-10-12 第 11 轮) —— `_spread` 用**统一** `gap = max(text_w, text_h)`
      在**周长参数**上推开。这只在「同一 1:1 边 + 同向文字」时正确; 矩形上不同位置
      的「1 单位周长对应的物理宽度」不同 (横边=Δx, 竖边=Δy, 而 x/y 尺度不同), 且角点
      附近文字会跨边 → 实测底部边 `SRR15037xxx` 一批仍残留 **12.5x4.4 px** 重叠 65 对。

    ⚠⚠ 二次修正 (同轮): 初版 `_need` 用 `pairs[i][3]` 里**预存**的边标签, 但那是
      **移动前**的位置 —— 标签被推开后会跨到别的边, 判据随即失效 (实测残留 15 对
      大交叠 `55x4.4px`)。这里改为**每轮迭代按当前 t 动态判边**。

    三次修正: 每字符宽系数由 `0.56` 改为实测 `0.8656` (旧值低估 35%)。

    pairs: [(t, w_text, h_text, _unused), ...]
    返回新的参数列表。
    """
    v = [p[0] for p in pairs]
    n = len(v)
    if n < 2:
        return v
    ww, hh = X1 - X0, Y1 - Y0

    def _need(i, j):
        """第 i、j 个标签之间的最小周长间距 (**按预存边标签**)。

        ⚠ 试过"每轮按当前 t 动态判边", 实测**更差** (ring9: 15 → 51 对):
          判据随迭代变化 → 相邻对的 `need` 反复跳变 → 震荡不收敛。
          故保留"按移动前位置判边"——它稳定, 且语义正确 (决定标签大致在哪条边)。
        """
        wi, hi = pairs[i][1], pairs[i][2]
        wj, hj = pairs[j][1], pairs[j][2]
        ei, ej = pairs[i][3], pairs[j][3]
        horiz = ('b', 't')
        if (ei in horiz) and (ej in horiz):
            # 都沿横边排布: 水平中心距 ≥ 半宽之和
            return 0.5 * (wi + wj)
        if (ei not in horiz) and (ej not in horiz):
            # 都沿竖边排布: 垂直中心距 ≥ 半高之和
            return 0.5 * (hi + hj)
        # 跨角点 (一横一竖): 用对角线保守
        return 0.5 * (math.hypot(wi, hi) + math.hypot(wj, hj))

    for _ in range(800):
        moved = False
        for j in range(n):
            i = j - 1
            ii = i if j > 0 else n - 1
            d = (v[j] - v[i]) if j > 0 else (v[0] + per - v[-1])
            need = max(_need(ii, j), gap_min)
            if d < need:
                shift = (need - d) / 2.0
                if j > 0:
                    v[i] -= shift
                    v[j] += shift
                else:
                    v[-1] -= shift
                    v[0] += shift
                moved = True
        if not moved:
            break
    return v


def draw_peripheral_labels(ax, tips, x0, x1, y0, y1, w, h,
                           fontsize=3.4, n_rings=3, color_of=None,
                           pad=0.030, ring_gap=0.055, data_per_inch=None,
                           label_chars=20):
    """把 taxon 标签推到网络外圈的多层环上, 用细引线连回 tip。

    tips: {nodeId: (x, y)}; color_of: {nodeId: color} 可选。
    data_per_inch: 1 英寸对应多少数据单位 (决定标签的最小物理间距)。
    返回实际画出的标签数。
    """
    if not tips:
        return 0
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    # ⚠⚠ 每字符宽度的**实测系数** (2026-10-12 第 11 轮):
    #   DejaVu Sans 渲染 `SRR17208227_OR489165.1` (22 char, 3.2 pt) 真实宽
    #   = 0.6093 in → 每字符 0.02770 in = **0.8656 * fontsize/72**。
    #   旧代码用 `0.56 * fontsize/72` (记作 text_w / label_chars) → **低估 35 %**,
    #   导致 `_spread_rect` 算出的半宽之和偏小 → 同横边相邻标签仍水平重合
    #   (实测残留 40 对 `dy==0`、水平交叠 40–60 px)。
    #   改用量出来的系数, 并**按每个标签自己的字符数**算 (不统一用 label_chars)。
    _CHW = 0.8656 * fontsize / 72.0
    dpi_data = float(data_per_inch) if data_per_inch else max(w, h) / 9.32
    text_h = 1.55 * fontsize / 72.0 * dpi_data
    text_w = label_chars * _CHW * dpi_data
    # ⚠ 试过把 gap 抬到 1.14*text_h 清 rings=6 的 0.7–3.5 px 细缝 —— **反而更差**
    #   (8→10 对): 那些缝不是 `_spread_rect` 的 margin 造成, 而是**逃逸搜索**落点
    #   的先后顺序造成。调 gap 会改变整条排布的连锁, 不能治本。故保留 1.06。
    gap_min = text_h * 1.06
    _CHAR_W = _CHW * dpi_data

    items = sorted(((math.atan2(y - cy, x - cx), i, x, y)
                    for i, (x, y) in tips.items()), key=lambda t: t[0])
    n = len(items)
    n_rings = max(1, min(n_rings, n))
    drawn = 0
    # ⚠⚠ 环**间**协同 (2026-10-12 第 11 轮): 各环独立排布后, 同一角度方向上被分到
    #   不同环的标签会**横向对撞** —— 实测残留 15 对里 `dx` 呈 133/266/399 px 阶梯
    #   (= `need` 的整数倍, 同边排布特征) 而 `dy` 非零 (= 不同环) 正是此症。
    #   解法: 内层环已放置的标签框作为**外层障碍**, 落位后做一次矩形求交, 撞了就沿
    #   该环周长向外挪一个身位。障碍存 (x0,x1,y0,y1) 数据坐标 —— 同一坐标空间, 无需换算。
    _placed = []          # [(x0, x1, y0, y1), ...] 数据坐标

    def _rect_of(bx, by, ha, va, wtxt, htxt):
        if ha == 'left':
            x0, x1 = bx, bx + wtxt
        elif ha == 'right':
            x0, x1 = bx - wtxt, bx
        else:
            x0, x1 = bx - 0.5 * wtxt, bx + 0.5 * wtxt
        if va == 'bottom':
            y0, y1 = by, by + htxt
        elif va == 'top':
            y0, y1 = by - htxt, by
        else:
            y0, y1 = by - 0.5 * htxt, by + 0.5 * htxt
        return x0, x1, y0, y1

    def _hits(box):
        # ⚠ 试过给判据加净空 `_CLR` 清 rings=6 的细缝 —— **显著更差** (10→23 对):
        #   环本身容量不足时, 加净空会让后续标签被"连锁挤开", 反而制造新重叠。
        #   结论: rings=6 对 145 个标签**容量不足**, 不是判据问题 (生产默认 9 已够)。
        #   保留严格相交判据。
        for o in _placed:
            if (box[0] < o[1] and box[1] > o[0]
                    and box[2] < o[3] and box[3] > o[2]):
                return True
        return False

    for k in range(n_rings):
        ring = items[k::n_rings]
        if not ring:
            continue
        mx = (pad + k * ring_gap) * w
        my = (pad + k * ring_gap) * h
        X0, X1, Y0, Y1 = x0 - mx, x1 + mx, y0 - my, y1 + my
        per = 2 * ((X1 - X0) + (Y1 - Y0))
        recs = []
        for ang, i, x, y in ring:
            bx, by = _border_point(cx, cy, x, y, X0, X1, Y0, Y1)
            recs.append([_rect_param(bx, by, X0, X1, Y0, Y1, per), i, x, y])
        recs.sort(key=lambda r: r[0])

        def _edge_of(t):
            tt = t % per
            ww, hh = X1 - X0, Y1 - Y0
            if tt < ww:
                return 'b'
            if tt < ww + hh:
                return 'r'
            if tt < 2 * ww + hh:
                return 't'
            return 'l'

        pairs = []
        for r in recs:
            lab = tips_label_of(r[1]) or ''
            wi = max(1, len(lab)) * _CHAR_W
            pairs.append((r[0], wi, text_h, _edge_of(r[0])))
        pos = _spread_rect(pairs, per, X0, X1, Y0, Y1, gap_min)
        for (t_old, i, x, y), t_new in zip(recs, pos):
            lab = tips_label_of(i) or ''
            wi = max(1, len(lab)) * _CHAR_W
            # 环间避让 (2026-10-12 第 18 轮修):
            #   ⚠ 旧实现只**正向**挪 `+1.15*wi`, 且 40 次用尽就**放弃**(直接放置)。
            #     等长布局把网络压扁 (canvas 高 575→419 px) → 相邻环的横边 y 更近
            #     → 同 y 对撞更易发生; 实测残留 **1 对** (组名仍 0)。
            #   修法: 双向搜索 (先 +, 再 -), 步长从细到粗, 命中即停。
            #     判据仍是 `_hits` 对**已放置框**求交 (数据坐标, 同一空间)。
            def _ok(tt):
                _bx, _by, _ha, _va = _rect_at(tt, X0, X1, Y0, Y1, per)
                return not _hits(_rect_of(_bx, _by, _ha, _va, wi, text_h))
            if not _ok(t_new):
                _esc = False
                for _step_mult in (1.15, 1.6, 2.4):
                    _stp = _step_mult * wi
                    for _k in range(1, 61):
                        for _t2 in ((t_new + _k * _stp) % per,
                                    (t_new - _k * _stp) % per):
                            if _ok(_t2):
                                t_new = _t2
                                _esc = True
                                break
                        if _esc:
                            break
                    if _esc:
                        break
            bx, by, ha, va = _rect_at(t_new, X0, X1, Y0, Y1, per)
            _placed.append(_rect_of(bx, by, ha, va, wi, text_h))
            col = (color_of or {}).get(i, C_TIP)
            ax.plot([x, bx], [y, by], '-', color=C_LEADER, lw=0.28,
                    solid_capstyle='round', zorder=2)
            ax.text(bx, by, tips_label_of(i), fontsize=fontsize, ha=ha, va=va,
                    color=col, zorder=6,
                    path_effects=[pe.withStroke(linewidth=_LBL_STROKE * 0.45,
                                                foreground='white')])
            drawn += 1
    return drawn


# 占位 (在 render 内被替换成闭包, 保持 draw_peripheral_labels 签名简单)
def tips_label_of(_i):          # pragma: no cover - 运行时被覆盖
    return ''


def draw_tip_labels(ax, tips, x0, x1, y0, y1, w, h, fontsize=3.2,
                    color_of=None, off=0.012):
    """把标签直接挂在 tip 外侧 (SplitsTree 默认观感), 按方位自动对齐。

    tips: {nodeId: (x, y)}; color_of: {nodeId: color} 可选。
    """
    if not tips:
        return 0
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    for i, (x, y) in tips.items():
        dx, dy = x - cx, y - cy
        d = math.hypot(dx, dy) or 1.0
        ux, uy = dx / d, dy / d
        tx, ty = x + ux * off * w, y + uy * off * h
        if ux > 0.30:
            ha = 'left'
        elif ux < -0.30:
            ha = 'right'
        else:
            ha = 'center'
        va = 'bottom' if uy > 0.30 else ('top' if uy < -0.30 else 'center')
        ax.text(tx, ty, tips_label_of(i), fontsize=fontsize, ha=ha, va=va,
                color=(color_of or {}).get(i, C_TIP), zorder=6,
                path_effects=[pe.withStroke(linewidth=_LBL_STROKE * 0.45,
                                            foreground='white')])
    return len(tips)


# ──────────────────────────────────────────────────────────────────────────────
# 主渲染
# ──────────────────────────────────────────────────────────────────────────────
def _text_w_in(s, fontsize):
    """文本的**真实字形宽度** (英寸)。用 TextPath 量, 不用"字符数×系数"估。

    项目硬规则: 题注/标注必须量宽, 不许目测 (见 MEMORY.md 第 7 轮事故)。
    """
    try:
        from matplotlib.textpath import TextPath
        from matplotlib.font_manager import FontProperties
        tp = TextPath((0, 0), str(s), size=fontsize,
                      prop=FontProperties(weight='bold'))
        return tp.get_extents().width / 72.0
    except Exception:  # noqa: BLE001
        return 0.62 * fontsize * len(str(s)) / 72.0


def _weight_ref(ws, pct=_LW_REF_PCT):
    """归一化分母: 权重的 `pct` 分位 (默认 p95)。

    ⚠ 不用 max —— 本数据权重右偏 (p50=0.00026 vs max=0.01461), 用 max 会把
    99 % 的边压到接近地板, 粗细差异消失。p95 让最粗档由"显著粗的边"占据。
    """
    xs = sorted(v for v in ws if v > 0)
    if not xs:
        return 1.0
    k = int(len(xs) * pct / 100.0) - 1
    k = max(0, min(k, len(xs) - 1))
    return xs[k] or 1.0


def _lw_fn(ws, lw_min=_LW_MIN, lw_max=_LW_MAX, gamma=_LW_GAMMA,
           ref_pct=_LW_REF_PCT):
    """返回 `weight -> linewidth(pt)` 闭包 (见文件头"线宽映射"说明)。"""
    ref = _weight_ref(ws, ref_pct)
    def _f(w):
        t = (w / ref) if ref > 0 else 0.0
        if t > 1.0:
            t = 1.0
        if gamma != 1.0:
            t = t ** gamma
        return lw_min + (lw_max - lw_min) * t
    return _f


def edge_length_stats(nodes, edges):
    """几何边长统计 → (n, min, p50, p95, max, max/min)。用于判"是否等长布局"。

    判据: `max/min <= _UNW_TOL` (1.35) → 该坐标已是 unweighted/等长布局。
    """
    ls = sorted(math.hypot(nodes[a][0] - nodes[b][0], nodes[a][1] - nodes[b][1])
                for a, b, *_ in edges if a in nodes and b in nodes)
    ls = [v for v in ls if v > 1e-9]
    if not ls:
        return 0, 0.0, 0.0, 0.0, 0.0, 0.0
    q = lambda p: ls[min(int(len(ls) * p / 100.0), len(ls) - 1)]  # noqa: E731
    return len(ls), ls[0], q(50), q(95), ls[-1], (ls[-1] / ls[0])


def render(coords, out_pdf, out_png=None, nexus=None, title=None,
           groups=None, group_note=None, color_splits=False,
           label_mode='tip', label_rings=3, label_fontsize=3.2,
           label_min_taxa=0, width_in=None, fit_text=True, show_legend=True,
           group_arcs=True, arc_label_fontsize=8.0, lw_preset='auto',
           lw_scale=1.0, canvas_aspect=None):
    # ⚠ 2026-10-13: `lw_preset` 默认由 `'std'` 改为 `'auto'` —— 与 CLI 默认对齐。
    #   旧默认 `'std'` 会让**等长坐标**也走权重线宽 (双重编码): 管线调 `render()`
    #   时未传该参 → 拿到 `'std'` → auto-flat 永不触发 → 图是"等长几何 + 变宽线"
    #   的半成品。`'auto'` 对**加权**输入的行为与 `'std'` **逐点一致**
    #   (ratio > _UNW_TOL → 回落 'std'), 故对既有加权图零回归。
    #
    # ⚠ 2026-10-13: `canvas_aspect` (目标画布比 h/w)。**不传** = 按数据 bbox 出
    #   画布 (旧行为, 图会被裁到紧贴网络)。**传值** = 画布按该比出图, 网络用
    #   `set_aspect('equal')` 在画布内**居中** —— 用于把面板画布比**钉到同排
    #   邻居的比**, 否则 Fig3 里两个并排面板高度不等 → 面板字母错位
    #   (实测 useWeights=false 后面板 D 比从 0.9969 变 0.9105, 字母 D 比 C 低
    #   0.41 in; 对齐到 C 的 1.0 即消除)。
    _setup_font()
    nodes, labels, edges = load_coords(coords)
    x0, x1, y0, y1 = coords_bbox(nodes)
    w, h = x1 - x0, y1 - y0
    asp = (w / h) if h else 1.0

    fig_w = float(width_in) if width_in else COL_W_1
    if canvas_aspect:
        fig_h = fig_w * float(canvas_aspect)
    else:
        fig_h = fig_w / asp if asp else fig_w
    rc = sci_rcparams(base_size=8) if callable(sci_rcparams) else {}
    matplotlib.rcParams.update(rc)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # ── 网络本体 ─────────────────────────────────────────────────────────────
    wmax = max((e[2] for e in edges), default=1.0) or 1.0
    # 等长布局 (unweighted) 自动降级为 **等粗**: 边已等长 → 再让线宽随权重变
    # 就是"用一种编码表达另一种编码", 且权重与几何已解耦 (重放时被当 1), 画出来
    # 只会是噪声。判据 = 几何边长 max/min <= _UNW_TOL。
    _n_e, _emin, _e50, _e95, _emax, _eratio = edge_length_stats(nodes, edges)
    _auto_flat = (_eratio <= _UNW_TOL) if lw_preset == 'auto' else False
    _eff_preset = 'flat' if (_auto_flat or lw_preset == 'flat') else lw_preset
    if lw_preset == 'auto' and not _auto_flat:
        _eff_preset = 'std'
    _preset = LW_PRESETS.get(_eff_preset, LW_PRESETS['std'])
    _lwf = _lw_fn([e[2] for e in edges], _preset[0] * lw_scale,
                  _preset[1] * lw_scale, _preset[2], _preset[3])
    if _auto_flat:
        sys.stderr.write('  [viz] 等长布局检测到 (边长 max/min=%.2f <= %.2f) '
                         '→ 线宽自动切 flat (等粗)\n' % (_eratio, _UNW_TOL))
    cmap = None
    if color_splits:
        cmap = matplotlib.colormaps['viridis']
    for a, b, wt, _sp in edges:
        if a not in nodes or b not in nodes:
            continue
        xa, ya = nodes[a]
        xb, yb = nodes[b]
        if cmap is not None:
            col = cmap(0.12 + 0.88 * (wt / wmax))
            lw = 0.35
        else:
            col = C_EDGE
            lw = _lwf(wt)          # 权重 → 线宽 (幂次压缩, 见文件头说明)
        ax.plot([xa, xb], [ya, yb], '-', color=col, lw=lw,
                solid_capstyle='round', zorder=1)

    # ── tip 点 (可选按分组着色) ──────────────────────────────────────────────
    gmap = {}
    if groups:
        gmap = _match_groups(labels.values(), groups) or {}
    gnames = sorted({v for v in gmap.values()})
    gcolor = _canonical_group_color(gnames)
    tip_color = {}
    for i, lab in labels.items():
        g = gmap.get(lab)
        tip_color[i] = gcolor.get(g, C_TIP) if g else C_TIP

    lx = [nodes[i][0] for i in labels]
    ly = [nodes[i][1] for i in labels]
    if gnames:
        for g in gnames:
            sel = [i for i in labels if gmap.get(labels[i]) == g]
            ax.plot([nodes[i][0] for i in sel], [nodes[i][1] for i in sel],
                    'o', ms=1.5, color=gcolor[g], mec='none', zorder=4)
        others = [i for i in labels if not gmap.get(labels[i])]
        if others:
            ax.plot([nodes[i][0] for i in others], [nodes[i][1] for i in others],
                    'o', ms=1.2, color='#9a9a9a', mec='none', zorder=4)
    else:
        ax.plot(lx, ly, 'o', ms=1.1, color=C_TIP, mec='none', zorder=4)

    # ── 标签 ────────────────────────────────────────────────────────────────
    global tips_label_of
    tips_label_of = lambda i: labels.get(i, '')
    n_drawn = 0
    _rings = int(label_rings) if str(label_mode).lower() == 'ring' else 0
    pad = 0.055 + 0.075 * _rings
    ax_w = (x1 - x0) * (1.0 + 2.0 * pad)
    dpi_data = (ax_w / fig_w) if fig_w else 1.0
    if str(label_mode).lower() != 'none' and len(labels) >= label_min_taxa:
        tips = {i: nodes[i] for i in labels}
        if _rings:
            max_chars = max((len(v) for v in labels.values()), default=16)
            n_drawn = draw_peripheral_labels(
                ax, tips, x0, x1, y0, y1, w, h,
                fontsize=label_fontsize, n_rings=_rings, color_of=tip_color,
                data_per_inch=dpi_data, label_chars=max_chars)
        else:
            n_drawn = draw_tip_labels(ax, tips, x0, x1, y0, y1, w, h,
                                      fontsize=label_fontsize,
                                      color_of=tip_color)

    # ── 图例/题注占位 (必须在分组弧线之前算: 决定坐标区可用高度) ────────────
    ncol = min(4, max(1, len(gnames))) if (gnames and show_legend) else 0
    nrow = (len(gnames) + ncol - 1) // ncol if ncol else 0
    bottom = (0.055 + 0.055 * nrow) if (nrow or fit_text) else 0.02
    top = 0.925 if title else 0.995
    # ⚠ 坐标区**实际**宽度不等于 0.976*fig_w: `set_aspect('equal')` 会按高度把
    #   axes 盒收缩 (默认 adjustable='box')。本图 asp≈1.003 而可用区 9.10x8.08 in
    #   是**高度受限** → 实际宽只有 8.08*1.003 ≈ 8.11 in。按 9.10 算会让"数据/英寸比"
    #   偏小 12% → 障碍框/组名间距全被低估 (实测残留 1 对 3.0x2.8 px 的组名压标签)。
    _axw_eff = min((0.988 - 0.012) * fig_w, (top - bottom) * fig_h * asp)

    # ── 外圈分组弧线 ────────────────────────────────────────────────────────
    # 必须在 set_xlim 之前画: 弧线/组名会伸到网络外, 返回值给出需要的最小边距。
    arc_pad = 0.0
    if group_arcs and groups and labels:
        _tips = {i: nodes[i] for i in labels}
        # 组名最小中心间距 (数据单位) = 组名字框对角线 × **真实**数据/英寸比。
        # ⚠ 真实比值 = w*(1+2*pad)/(0.976*fig_w), 含 pad —— 而 pad 又取决于组名落点,
        #   是个循环。故做**两遍**: 第一遍拿 arc_pad 估出 pad, 第二遍用真实比值重放。
        #   (只做一遍会低估 2 倍 → 实测 6 个组名里仍有 1 对压字。)
        _tw_in = max(_text_w_in(_short_group(g), arc_label_fontsize)
                     for g in set(groups.values()))
        _th_in = 1.25 * arc_label_fontsize / 72.0
        _diag_in = math.hypot(_tw_in, _th_in)
        _n_ln = len(ax.lines)
        # ⚠ 只删**本次新建**的 Text —— `ax.texts` 里还有 145 个 tip 标签,
        #   整列表 remove() 会把标签全删光 (实测踩过)。
        _tx_before = list(ax.texts)

        # pad 的**闭式解**: 组名是文字, 中心点进框不代表文字进画布 —— θ≈0/π 处会被
        #   右/左缘切掉。文字从中心向两侧各伸 _tw_in/2 英寸, 换成"占数据框宽的比例"要
        #   除以**坐标区实际宽** 0.976*fig_w (subplots_adjust left/right = 0.012/0.988),
        #   再乘 (1+2*pad) 把 pad 自身的放大算进去:
        #       pad = (arc_pad + 0.012) + k*(1 + 2*pad),  k = (tw/2)/(0.976*fig_w)
        #   ⇒  pad = (A + k) / (1 - 2k),  A = arc_pad + 0.012
        #   ⚠ 旧写法 `max(pad, arc_pad + 0.012 + (tw/2)/fig_w)` 漏了 (1+2pad)/0.976
        #     因子 → 低估约 2 倍, 实测面板版 "Guangzhou" 被右缘切掉。
        _k = (_tw_in / 2.0) / max(_axw_eff, 0.1)

        def _pad_from(ap):
            _A = ap + 0.012
            return (_A + _k) / max(1e-6, 1.0 - 2.0 * _k)

        # tip 标签的**占位框** (数据坐标): 整页版 145 个标签会与组名抢同一片区域,
        #   实测 24 对重叠全是"组名 vs tip 标签"。把它们的真实 bbox 量出来当障碍物,
        #   比"整体加大组名外推量"精准 —— 后者是**全局**代价 (实测 pad 0.085→0.30
        #   才清零组名, 但 tip 标签互压从 403 涨到 494)。
        #   ⚠⚠ 障碍框**只能存"锚点 + 对齐 + 英寸尺寸"**, 不能存"bbox 中心":
        #     文字尺寸是物理量, 在数据坐标下的大小 = 英寸 × 数据/英寸比。而 ha='left'
        #     的文字, 其 bbox **中心** = 锚点 + 宽/2 —— 中心会随比值漂移, 不是不变量
        #     (实测按中心存 → 障碍中心偏 15 数据单位 → 判据放行了一个实际压字的落点)。
        #     锚点 (ax.text 的 x,y) 才是数据坐标不变量。
        _obs_in = []
        if str(label_mode).lower() != 'none' and n_drawn:
            try:
                fig.canvas.draw()
                for _t in ax.texts:
                    _b = _t.get_window_extent(renderer=fig.canvas.get_renderer())
                    _ax_, _ay_ = _t.get_position()
                    _obs_in.append((float(_ax_), float(_ay_),
                                    _b.width / fig.dpi, _b.height / fig.dpi,
                                    str(_t.get_ha()), str(_t.get_va())))
            except Exception:                   # noqa: BLE001
                _obs_in = []

        def _mk_obs(_dpi):
            """把障碍框按**当轮**的数据/英寸比展开成 (x0, x1, y0, y1) 盒子。"""
            out = []
            for _ax_, _ay_, _w, _h, _ha, _va in _obs_in:
                _w *= _dpi
                _h *= _dpi
                if _ha == 'left':
                    _x0, _x1 = _ax_, _ax_ + _w
                elif _ha == 'right':
                    _x0, _x1 = _ax_ - _w, _ax_
                else:
                    _x0, _x1 = _ax_ - 0.5 * _w, _ax_ + 0.5 * _w
                if _va == 'bottom':
                    _y0, _y1 = _ay_, _ay_ + _h
                elif _va == 'top':
                    _y0, _y1 = _ay_ - _h, _ay_
                else:
                    _y0, _y1 = _ay_ - 0.5 * _h, _ay_ + 0.5 * _h
                out.append((_x0, _x1, _y0, _y1))
            return out

        # ── 迭代到收敛 (不用固定"两遍") ────────────────────────────────────
        # 组名落位 ↔ 需要的最小边距 pad ↔ **数据/英寸比** 三者互为因果:
        #   落位用比值 → 落位决定 pad → pad 决定真实比值。
        # 只做 1~2 遍都留下残差: 实测"两遍"后仍剩 1 对 3.0x2.8 px (组名 vs tip 标签),
        # 因为 pad 在第二遍之后还会被 `_pad_from` 抬高一点 → 真实比值比落位时用的
        # 估计值大 → 障碍框(按比值换算)被低估。
        # 判据: 相邻两轮的比值相对变化 < 0.1% 即认为收敛 (本数据 3 轮内收敛)。
        _dpi_cur = dpi_data
        _pad_cur = pad
        for _it in range(6):
            for _t in list(ax.texts):
                if _t not in _tx_before:
                    _t.remove()
            for _ln in list(ax.lines[_n_ln:]):
                _ln.remove()
            _, _, arc_pad = draw_group_arcs(
                ax, _tips, labels, groups, fontsize=arc_label_fontsize,
                min_sep=(_tw_in * _dpi_cur, _th_in * _dpi_cur),
                obstacles=_mk_obs(_dpi_cur))
            _pad_cur = max(pad, _pad_from(arc_pad))
            _dpi_new = ((x1 - x0) * (1.0 + 2.0 * _pad_cur)
                        / max(_axw_eff, 0.1))
            if abs(_dpi_new - _dpi_cur) <= 1e-3 * max(abs(_dpi_cur), 1e-9):
                _dpi_cur = _dpi_new
                break
            _dpi_cur = _dpi_new
        pad = _pad_cur

    # ── 画布与坐标轴 ────────────────────────────────────────────────────────
    ax.set_xlim(x0 - pad * w, x1 + pad * w)
    ax.set_ylim(y0 - pad * h, y1 + pad * h)
    ax.set_aspect('equal')
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=9.0, fontweight='bold', color='#222222')

    # ── 图例 ────────────────────────────────────────────────────────────────
    # ⚠ 顺序很重要: 先 subplots_adjust 定好边距, 再放 legend / caption。
    #   反过来 (tight_layout 之后再 fig.legend) 会把图例压到坐标区上并叠字。
    #   (ncol / nrow / bottom / top 已在分组弧线之前算好 —— 那里要用它推坐标区实际宽)
    fig.subplots_adjust(left=0.012, right=0.988, top=top, bottom=bottom)
    if ncol:
        from matplotlib.lines import Line2D
        counts = {g: sum(1 for i in labels if gmap.get(labels[i]) == g)
                  for g in gnames}
        handles = [Line2D([0], [0], marker='o', ls='none', ms=3.2,
                          color=gcolor[g], label='%s (n = %d)' % (g, counts[g]))
                   for g in gnames]
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.5, bottom - 0.045 * nrow - 0.012),
                   ncol=ncol, frameon=False, fontsize=6.2, handlelength=0.9,
                   columnspacing=1.4, handletextpad=0.4)

    # ── 题注 ────────────────────────────────────────────────────────────────
    if fit_text:
        fit = None
        n_splits = len({e[3] for e in edges})
        if nexus and os.path.exists(nexus):
            m = re.search(r'fit=([\d.]+)', open(nexus, errors='ignore').read())
            if m:
                fit = float(m.group(1))
        cap = ('Neighbor-Net (SplitsTree6 native layout): n = %d taxa, %d splits'
               % (len(labels), n_splits))
        if fit is not None:
            cap += ', fit = %.1f%%' % fit
        if gnames:
            cap += '; tip colour = %s' % (group_note or 'group')
        fig.text(0.5, 0.008, cap, ha='center', va='bottom',
                 fontsize=6.0, color='#333333')

    # 2026-10-16: 双格式 —— 调用方没给 out_png 时按主名自动补 .png (PNG 必出,
    # 例: Fig3 面板 D 的 neighbornet_panel.pdf 现在也有同名 PNG)
    save_both(fig, out_pdf, png=out_png, dpi=300, facecolor='white')
    plt.close(fig)
    print('saved: %s (%d nodes, %d edges, %d labels, %d group colours)'
          % (out_pdf, len(nodes), len(edges), n_drawn, len(gnames)))


# ── 等长 (unweight) 内联: 让本脚本自己就能出"等长 + 等粗" ────────────────────
# 用户 2026-10-13: "修改我们的脚本绘制等长等粗"。
# 此前 unweight 只做在 `phylo_pipeline.py` 里 (单独 import 转换器), 直接跑本脚本
# 只能画**加权**布局 —— 想画等长得先手动转一次坐标。现在内联进 CLI:
#   `--unweight auto` (默认): 输入非等长 (边长 max/min > _UNW_TOL) 就自动转;
#   `--unweight on`         : 总是转;  `off`: 原样用 (复现旧加权图时配 --lw-preset legacy)。
_UNW_MODULE_CANDIDATES = ('splitstree_unweight.py', '_unweight.py')


def find_unweight_module():
    """定位同目录下的 unweight 转换器 (`splitstree_unweight.py` 优先)。"""
    here = os.path.dirname(os.path.abspath(__file__))
    for name in _UNW_MODULE_CANDIDATES:
        p = os.path.join(here, name)
        if os.path.exists(p):
            return p
    return None


def unweight_coords(src, dst, module_path=None):
    """加权坐标 → **等长**坐标 (同角度把权重当 1 重放 DFS)。

    返回 `(ok, 路径)`。⚠ `splitstree_unweight.main()` **无参**, 从 `sys.argv[1]/[2]`
    取输入/输出 → 这里必须临时替换 `sys.argv` 再恢复 (踩过: 直接 `main([a,b])` 报 TypeError)。
    """
    mp = module_path or find_unweight_module()
    if not mp:
        return False, src
    import importlib.util as _ilu
    spec = _ilu.spec_from_file_location('_unweight_rt', mp)
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    saved = sys.argv
    try:
        sys.argv = [mp, src, dst]
        rc = mod.main()
    finally:
        sys.argv = saved
    ok = (rc == 0 and os.path.exists(dst))
    return ok, (dst if ok else src)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--coords', required=True, help='SplitsTree6 native coords txt')
    ap.add_argument('-o', '--pdf', required=True)
    ap.add_argument('--png', default=None)
    ap.add_argument('--nexus', default=None, help='取 fit 值 (可选)')
    ap.add_argument('--title', default=None)
    ap.add_argument('--meta', default=None, help='元数据 CSV/TSV (含 name + 分组列)')
    ap.add_argument('--group-col', default='auto')
    ap.add_argument('--color-splits', action='store_true',
                    help='按 split 权重 viridis 着色 (默认统一深色线)')
    ap.add_argument('--label-mode', default='tip', choices=['tip', 'ring', 'none'],
                    help='tip=标签挂 tip 外侧 (SplitsTree 观感, 默认); '
                         'ring=推到外圈多层环 + 引线; none=不画标签')
    ap.add_argument('--label-rings', type=int, default=3,
                    help='--label-mode ring 时的环数')
    ap.add_argument('--label-fontsize', type=float, default=3.2)
    ap.add_argument('--width-in', type=float, default=None,
                    help='图宽 (in); 默认通栏 %s' % COL_W_1)
    ap.add_argument('--no-legend', action='store_true',
                    help='不画分组图例 (面板版用, 图例交给外层题注)')
    ap.add_argument('--legend', action='store_true',
                    help='强制画分组图例 (默认: 画了分组弧线就不画图例 —— 组名已在图内)')
    ap.add_argument('--no-group-arcs', action='store_true',
                    help='不画外圈分组弧线 (只按组给 tip 上色 + 图例)')
    ap.add_argument('--arc-label-fontsize', type=float, default=8.0,
                    help='分组弧线旁的组名字号')
    ap.add_argument('--no-caption', action='store_true', help='不画底部题注')
    ap.add_argument('--lw-preset', default='auto', choices=list(LW_PRESETS) + ['auto'],
                    help='线宽映射预设 (auto=按布局自动选, 等长→flat/加权→std; '
                         'flat=等粗; legacy=旧线性行为; std/bold/subtle=幂次压缩)')
    ap.add_argument('--lw-scale', type=float, default=1.0,
                    help='线宽整体倍率 (在预设之上再乘, 半栏缩印可给 1.2)')
    ap.add_argument('--unweight', default='auto', choices=['on', 'off', 'auto'],
                    help='等长转换 (权重当 1 重放 DFS): on=总是转; off=用原坐标; '
                         'auto=输入非等长时自动转 (默认) —— 等长 + --lw-preset auto '
                         '即"等长等粗"')
    ap.add_argument('--unweight-out', default=None,
                    help='等长坐标保存路径 (默认 <coords>_unweighted.coords)')
    a = ap.parse_args()

    # ── 等长转换 (必须在渲染之前) ────────────────────────────────────────────
    coords_in = a.coords
    _n0, _l0, _e0 = load_coords(coords_in)
    _er0 = edge_length_stats(_n0, _e0)[5]
    _need_unw = (a.unweight == 'on') or (a.unweight == 'auto' and _er0 > _UNW_TOL)
    if _need_unw:
        _dst = a.unweight_out or (os.path.splitext(coords_in)[0] + '_unweighted.coords')
        _ok, _new = unweight_coords(coords_in, _dst)
        if _ok:
            coords_in = _new
            _n1, _l1, _e1 = load_coords(coords_in)
            _er1 = edge_length_stats(_n1, _e1)[5]
            print('[viz] 等长转换 OK: 边长 max/min %.3f -> %.3f  (%s)'
                  % (_er0, _er1, coords_in), file=sys.stderr)
        else:
            print('[viz] 等长转换失败 (缺 %s?), 用原坐标 (边长 max/min=%.3f)'
                  % ('/'.join(_UNW_MODULE_CANDIDATES), _er0), file=sys.stderr)
    elif a.unweight == 'auto':
        print('[viz] 输入已是等长布局 (边长 max/min=%.3f), 跳过转换' % _er0,
              file=sys.stderr)

    groups, note = (None, None)
    if a.meta:
        raw, note = load_group_map(a.meta, a.group_col)
        if not raw:
            print('[warn] 元数据无可用分组列 (%s), 只画网络本体' % a.group_col,
                  file=sys.stderr)
        else:
            groups = raw
    # 图例默认策略: 画了分组弧线 → 组名已在图内, 不再叠一个图例; 无弧线 → 画图例
    show_legend = False if a.no_legend else (True if a.legend else a.no_group_arcs)
    render(coords_in, a.pdf, a.png, nexus=a.nexus, title=a.title,
           groups=groups, group_note=note, color_splits=a.color_splits,
           label_mode=a.label_mode, label_rings=a.label_rings,
           label_fontsize=a.label_fontsize, width_in=a.width_in,
           show_legend=show_legend, fit_text=not a.no_caption,
           group_arcs=not a.no_group_arcs,
           arc_label_fontsize=a.arc_label_fontsize,
           lw_preset=a.lw_preset, lw_scale=a.lw_scale)


if __name__ == '__main__':
    main()
