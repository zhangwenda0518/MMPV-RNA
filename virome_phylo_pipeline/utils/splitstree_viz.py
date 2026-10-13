#!/usr/bin/env python3
"""
splitstree_viz.py -- SplitsTree6 NeighborNet nexus 的等角布局 (equal-angle) 出版级绘制
解析 nexus 的 TAXA/CYCLE/SPLITS 块, 用 zonotope (平行四边形 zone) 渲染分裂网络,
视觉等价 SplitsTree 的 equal-angle 视图。

2026-10-12 新增: 外圈彩色分组弧线 (--meta + --group-col, 默认 location→host 自动回退)。
同组 tip 在 CYCLE 上取角度区间, 画细闭合弧线 (内边贴 tip 标签外沿, 外边按弧线嵌套深度
逐层外推), 组名同色标注 —— 即文献里 NeighborNet 按地理/宿主分组的常见画法。

2026-10-12（审稿人视角）三处重做:
  1. **分组弧线全 0 条的真 bug**: SplitsTree 会把 tip 名里的 '.' 规范成 '_'
     (`CRR1126135_OR489165.1` → `CRR1126135_OR489165_1`), 而元数据里还是带点的
     原值 —— 直接 `dict` 查表必然 0 命中。现在两侧都做 `_norm_label()` 归一
     (小写 + 非字母数字→'_' + 折叠), 再回退"去掉后缀只留登录号前缀"的模糊匹配。
  2. **网状信号量化**: 分裂系统里"互不相容(交错)"的 split 对就是 NeighborNet
     上的**箱形 (box) 结构** —— 重组 / 不完全谱系分选的直接证据。只画网络不报
     这个数, 读者无法判断"这张图到底有没有网状信号"。现在题注给出
     `交错对数 / 总对数 (占比)`, 网络本体按 **split 权重**着色 (viridis + 色标)。
     ⚠ 不要"把参与箱形的 split 标红": 实测 242 个 zone 里 241 个都参与过至少
     一对交错 —— 全红等于没标。
  3. **样式真源**: 原先本模块不调 `apply_sci_style`, 于是吃进程级全局 rcParams
     (顺序敏感 → 不可复现, 与 Fig3 同一类缺陷)。现在显式
     `rcParams.update(sci_rcparams(...))`; 分组弧线颜色也由 `tab20` 改为项目
     统一的 OKABE_ITO。

2026-10-12（第 6 轮, 用户: "这有什么意义吗, 我说的 splitstree 的结果呢 / 你 splitstree
的主题没了, 绘制类似的图"）—— 旧图是一颗"灰色太阳"(145 根等长辐条 + 中心空洞),
根本不像分裂网络。逐个定位到**四个真 bug**:

  1. **pendant split 全丢 / 画反** (最致命): SplitsTree 写 nexus 时每个 split 只写
     一侧, 且实测本数据 271 个 split 里**只有 1 个写成单 taxon**, 其余写的都是**长侧**。
     旧版对单 taxon 侧调用 `contiguous()` 得 None → 回退补集 → 把 pendant split 画成
     绕圆一周的巨型 zone, tip 永远没有自己的 zone。→ `equal_angle_zones` 重写为
     "两侧都试, 取较短的一侧作规范表示"。
  2. **补集坐标空间错**: `set(range(1, n + 1)) - set(ps)` 而 ps 是 0-based → 补集恒含
     非法下标 n、恒丢 0。
  3. **tip 被强行拉到统一外圈**: 旧版画一条 `r0 → R*1.03` 的灰线, 于是所有 tip 落在一个
     圆上, 外轮廓永远是圆。→ 改用**自然半径**(tip = 自身最外层 zone 的外弦), 外轮廓成为
     不规则凸多边形 (与文献 NeighborNet 一致); 画布按实际外接框自适应铺满。
  4. **分组弧"整组一条最短包弧"**: 本数据地理组在 CYCLE 上**高度交错**(Yinchuan 69 个
     tip 分成 19 段), 整组包弧必然横扫大半个圆, 把不相干的 tip 圈进来。→ 改为**逐连续段
     画贴边色带**, 组名只挂最大一段并做径向去重叠; 另加**组别图例**与**按组着色的 tip 点**。
     同时修掉组名撞名 (`China, Beijing` 与 `China, Beijing, Beijing` 简称都是 Beijing
     → 旧版画出两条都叫 Beijing 的弧、Xining 消失), 现按简称合并。
  5. **网状信号量化改用教科书写法**: 两个 split 不相容 ⟺ 四个交集全非空。旧版用"两弧
     交错"且弧取的是 nexus 写的长侧 → 报出 626 对 (2.1 %) 的虚高数字; 现为 184/29 403
     (0.6 %), 与侧的选择无关。
  6. 网络线默认改为**统一深色**(权重只控线宽), 与文献 NeighborNet 一致; 想要权重着色 +
     色标用 `--color-splits`。

验证: 1930 个嵌套 split 对 **0 个违反平面性** (d(σ) ≥ d(τ)+w(τ) 对所有 σ⊊τ 成立)。

用法: python3 splitstree_viz.py -i neighbornet.nexus -o neighbornet.pdf [--png out.png]
      [--meta sample_metadata.csv --group-col location] [--color-splits]
输出: PDF (矢量, 期刊首选) + 可选 PNG 300dpi; 白底, 全英文标签 (分组值含中文时自动切 CJK 字体)
"""

import argparse
import math
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np

try:                                    # 样式真源 (与本目录其它模块共享同一份)
    from utils.sci_palette import sci_rcparams, OKABE_ITO, COL_W_2
except ImportError:                     # 允许以脚本形式直接跑
    from sci_palette import sci_rcparams, OKABE_ITO, COL_W_2

# 分组弧线的几何常量 (单位: 网络外接半径 R=1.0 的比例)
_TIP_LINE_R = 1.03     # tip 悬枝终点
_TIP_LABEL_R = 1.10    # tip 标签半径
_ARC_BASE_OUT = 1.06   # 分组色带的最小外半径 (lim 兜底用)
_ARC_MIN_TIPS = 2      # 单 tip 的连续段不画色带
# 贴边分组色带 (2026-10-12 第 6 轮替代同心圆弧)
_RIB_PAD_IN = 0.018    # 内边离 tip 的距离
_RIB_BAND = 0.028      # 色带宽度
_RIB_LW = 1.5
# 组名最小间距 (数据坐标, 网络半径 ~1)
# ⚠ 第 7 轮：0.085 → 0.11。0.085 只是"中心点距"，而组名是横排文字
#   (8 pt 下 "Zhongning" 宽 ≈ 0.55 in ≈ 0.28 数据单位) → 两个中心相距 0.085
#   的文字在水平方向上仍会互相压字。0.11 + 白描边后实测不再相撞。
_LBL_MIN_D = 0.11
_LBL_STROKE = 2.6      # 组名白色描边宽度 (pt) —— 组名常落在网线之上, 无描边读不清

# tip 标签密度阈值: 超过就默认不画标签 (145 条挤一圈必然糊成一片)
_LABEL_AUTO_MAX = 45

C_ARC = OKABE_ITO          # 分组弧线配色


def _setup_font():
    """分组值可能含中文 (location/host 列), 有 CJK 字体就用, 没有保持默认。"""
    try:
        from matplotlib import font_manager
        names = {f.name for f in font_manager.fontManager.ttflist}
        for cand in ('Microsoft YaHei', 'SimHei', 'Noto Sans CJK SC',
                     'Source Han Sans SC', 'WenQuanYi Zen Hei', 'PingFang SC'):
            if cand in names:
                matplotlib.rcParams['font.sans-serif'] = [cand, 'DejaVu Sans']
                break
        matplotlib.rcParams['axes.unicode_minus'] = False
    except Exception:  # noqa: BLE001 — 字体探测失败不致命
        pass


def _norm_label(s):
    """标签归一化: 小写 + 非字母数字→'_' + 折叠重复下划线 + 去首尾下划线。

    SplitsTree 写 nexus 时会把 '.' 等字符换成 '_' (`CRR1126135_OR489165.1`
    → `CRR1126135_OR489165_1`), 而元数据里还是原值 —— 不归一则分组永远 0 命中。
    """
    s = str(s).strip().lower()
    s = re.sub(r'[^a-z0-9]+', '_', s)
    return re.sub(r'_+', '_', s).strip('_')


def _label_keys(name):
    """一个 tip 名 → 一组候选 key (精确归一 + 前缀模糊)。

    前缀模糊: 取第一段 (通常是登录号 SRR/CRR/OR 号), 用于 nexus 标签带了
    额外后缀而元数据没有 (或反之) 的情形。
    """
    n = _norm_label(name)
    keys = {n}
    if '_' in n:
        head = n.split('_')[0]
        if len(head) >= 4:
            keys.add(head)
    return keys


def load_group_map(meta_path, group_col="auto"):
    """元数据 → ({tip名: 组名}, 实际使用的列名)。

    group_col='auto' 时依次试 location → host (管线三通道表列名约定, 大小写不敏感);
    指定列不存在或全空 → ({}, None)。空值/缺失的 tip 归为未分组 (不画弧)。
    """
    try:
        import pandas as pd
        df = pd.read_csv(meta_path, encoding='utf-8-sig', sep=None, engine='python')
    except Exception:  # noqa: BLE001
        return {}, None
    cols = {str(c).lower().strip(): c for c in df.columns}
    name_col = cols.get('name') or df.columns[0]
    candidates = ['location', 'host'] if str(group_col).lower() in ('', 'auto') else [group_col]
    for cand in candidates:
        c = cols.get(str(cand).lower())
        if c is None:
            continue
        s = df[c].fillna('').astype(str).str.strip()
        if (s != '').sum() == 0:
            continue
        gmap = {str(n): v for n, v in zip(df[name_col].astype(str), s) if v}
        if gmap:
            return gmap, str(c)
    return {}, None


def _match_groups(labels, groups):
    """把 {元数据名: 组名} 映射到 nexus 的 labels 上 (归一化 + 前缀回退)。

    返回 {nexus标签: 组名}; 一个 nexus 标签命中多个元数据行时取首个。
    """
    if not groups:
        return {}
    exact, head = {}, {}
    for name, grp in groups.items():
        n = _norm_label(name)
        exact.setdefault(n, grp)
        if '_' in n:
            h = n.split('_')[0]
            if len(h) >= 4:
                head.setdefault(h, grp)
    out = {}
    for lab in labels:
        n = _norm_label(lab)
        if n in exact:
            out[lab] = exact[n]
            continue
        h = n.split('_')[0] if '_' in n else n
        if h in head:
            out[lab] = head[h]
    return out


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


def _contiguous_run(positions, n):
    """0-based cycle 位置集合 → 连续段 (圆形, 绕回起点重排); 不连续返回 None。"""
    if not positions:
        return None
    ps = sorted(set(positions))
    if len(ps) == n or len(ps) == 1:
        return ps
    gaps = [(ps[(i + 1) % len(ps)] - ps[i]) % n for i in range(len(ps))]
    big = max(range(len(ps)), key=lambda i: gaps[i])
    if all(g == 1 for i, g in enumerate(gaps) if i != big):
        start = big + 1
        return ps[start:] + ps[:start]
    return None


def equal_angle_zones(labels, cycle, splits):
    """每个 split 渲染为一个 zone (平行四边形, 四边均为网络边)。
    zone 半径 d = 包含它的 split 权重之和; zone 高 = 自身权重 w。

    ⚠ 2026-10-12 第 6 轮重写 —— 旧版两个致命缺陷, 直接导致整张图退化成"太阳":
      ① **侧的选择错**: SplitsTree 的 nexus 只写 split 的一侧, 且实测本数据
         271 个 split 里**只有 1 个写成单 taxon**(其余写的都是**长侧**)。旧版对
         单 taxon 侧调用 `contiguous()` 返回 None → 回退到补集 → 把 pendant
         split 画成绕圆一周的巨型 zone, tip 也就永远没有自己的 zone。
      ② **补集用错坐标空间**: `set(range(1, n + 1)) - set(ps)`, 而 ps 是
         0-based → 补集恒含非法的 n、恒丢 0。
    现在: 两侧都试, 取**较短的一侧**作规范表示 (圆形 split 系统两侧都连续)。
    """
    n = len(cycle)
    pos = {t: k for k, t in enumerate(cycle)}  # taxon-id -> cycle position
    sector = 2 * math.pi / n
    all_pos = set(range(n))

    zones = []
    for w, taxa in splits:
        if w <= 0 or len(taxa) in (0, n):
            continue
        ps = sorted(pos[t] for t in taxa if t in pos)
        if not ps:
            continue
        cands = []
        for poss in (ps, sorted(all_pos - set(ps))):
            run = _contiguous_run(poss, n)
            if run is not None:
                cands.append(run)
        if not cands:
            continue
        side = min(cands, key=len)          # 规范化为较短的一侧
        a, b = side[0], side[-1]
        zones.append({'taxa': frozenset(cycle[p] for p in side),
                      'a': a, 'b': b, 'span': len(side), 'w': w})
    # containment depth
    for z in zones:
        z['d'] = sum(o['w'] for o in zones
                     if o is not z and z['taxa'] < o['taxa'])
    return zones, n, sector


def _in_arc(x, a, s, n):
    """圆上坐标 x 是否落在"起点 a、长度 s"的弧内 (弧沿逆时针展开)。"""
    return (x - a) % n < s


def flag_incompatible(zones, n, cycle=None):
    """标出互不相容 (incompatible) 的 split 对 —— 即 NeighborNet 的箱形结构。

    两个 split σ=A|A', τ=B|B' 互不相容 ⟺ **四个交集 A∩B, A∩B', A'∩B, A'∩B'
    全部非空** (教科书写法)。不相容的 split 不可能同时存在于一棵树上 —— 这就是
    网状信号 (重组 / 不完全谱系分选) 的直接证据, 也是"数据到底是不是树状"的答案。

    ⚠ 2026-10-12 第 6 轮: 旧版用"两条弧是否交错"判定, 而弧取的是 nexus 里写的
    那一侧 (实测多半是**长侧**) → 两个几乎绕满圆的弧恒判为交错, 报出 626 对
    (2.1%) 的虚高数字。改为四象限判据后与侧的选择无关。
    """
    nz = len(zones)
    if nz < 2:
        return set(), 0, 0
    if cycle is None:
        cycle = list(range(1, n + 1))
    allset = set(cycle)
    S = [(set(z['taxa']), allset - set(z['taxa'])) for z in zones]

    flags = set()
    pairs = 0
    n_pairs = 0
    for i in range(nz):
        ai, ci = S[i]
        for j in range(i + 1, nz):
            aj, cj = S[j]
            n_pairs += 1
            if (ai & aj) and (ai & cj) and (ci & aj) and (ci & cj):
                pairs += 1
                flags.add(i)
                flags.add(j)
    return flags, pairs, n_pairs


# ═══════════════════════════════════════════════════════════════════
# 分组弧线: 文献 NeighborNet 分组画法
# ═══════════════════════════════════════════════════════════════════

def _runs(indices, n):
    """cycle 位置集合 → 连续段列表 (圆形, 首尾相接的段合并)。

    2026-10-12 第 6 轮: 本数据的地理组在 CYCLE 上**高度交错**(Yinchuan 69 个 tip
    分成 19 段、Zhongning 21 个分成 6 段) —— 用"整组一条最短包弧"必然横扫大半个
    圆, 把不相干的 tip 圈进来, 既难看又误导。改为**逐段**画贴边色带。
    """
    if not indices:
        return []
    idx = sorted(set(indices))
    if len(idx) == n:
        return [idx]
    runs = [[idx[0]]]
    for p in idx[1:]:
        if p - runs[-1][-1] == 1:
            runs[-1].append(p)
        else:
            runs.append([p])
    # 环绕: 末段接首段 (0 与 n-1 同属一段)
    if len(runs) > 1 and idx[0] == 0 and idx[-1] == n - 1:
        runs[0] = runs[-1] + runs[0]
        runs.pop()
    return runs


def _catmull_rom_closed(pts):
    """闭合点列 → Catmull-Rom 平滑的 matplotlib Path (三次贝塞尔)。"""
    from matplotlib.path import Path as MPath
    if len(pts) < 3:
        return None
    P = [np.asarray(p, float) for p in pts]
    P = P + P[:3]
    verts, codes = [tuple(P[1])], [MPath.MOVETO]
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        c1 = tuple(p1 + (p2 - p0) / 6.0)
        c2 = tuple(p2 - (p3 - p1) / 6.0)
        verts += [c1, c2, tuple(p2)]
        codes += [MPath.CURVE4] * 3
    return MPath(verts, codes)


def _short_group(grp):
    """组名取最细一级行政区 (`China, Ningxia, Yinchuan` → `Yinchuan`)。

    ⚠ 必须**合并同名**再画: 本数据的 location 列里 `China, Beijing` 与
    `China, Beijing, Beijing` 是两个值、简称都是 `Beijing` —— 旧版直接
    `split(',')[-1]` 就画出**两条都叫 Beijing** 的弧 (PDF 文本层里 Beijing
    出现两次、Xining 消失), 读者无法分辨。
    """
    return str(grp).split(',')[-1].strip() or str(grp)


def draw_group_ribbons(ax, cycle, sector, labels, groups, tip_r):
    """按"连续段"画贴边分组色带 (参考图 NeighborNet 分组画法)。

    返回 (legend_items, r_max_used)。legend_items = [(组名, 颜色, tip 数)]。
    """
    two_pi = 2 * math.pi
    n = len(cycle)
    lab_pos = {}
    for k, t in enumerate(cycle):
        lab = labels[t - 1] if t - 1 < len(labels) else str(t)
        lab_pos.setdefault(lab, []).append(k)

    # 按**简称**聚合 (同名合并, 见 _short_group)
    by_group = {}
    for lab, grp in (groups or {}).items():
        for k in lab_pos.get(lab, []):
            by_group.setdefault(_short_group(grp), []).append(k)
    if not by_group:
        return [], _TIP_LABEL_R

    def _pt(th, r):
        return (r * math.cos(th), r * math.sin(th))

    from matplotlib.patches import PathPatch
    _bbox = []                      # 供调用方自适应画布
    _labels = []
    legend, r_max_used = [], _TIP_LABEL_R
    items = sorted(by_group.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    for gi, (grp, ks) in enumerate(items):
        color = C_ARC[gi % len(C_ARC)]
        runs = [r for r in _runs(ks, n) if len(r) >= _ARC_MIN_TIPS]
        biggest = max(runs, key=len) if runs else None
        for run in runs:
            ths = [(k + 0.5) * sector for k in run]
            rs = [tip_r[k] for k in run]
            # 3 点取大 → 色带跨过低半径凹陷, 视觉上成为"包络"而非锯齿
            rs = [max(rs[max(0, i - 1):i + 2]) for i in range(len(rs))]
            ext = 0.35 * sector
            ths = [ths[0] - ext] + ths + [ths[-1] + ext]
            rs = [rs[0]] + rs + [rs[-1]]
            inner = [_pt(th, r + _RIB_PAD_IN) for th, r in zip(ths, rs)]
            outer = [_pt(th, r + _RIB_PAD_IN + _RIB_BAND)
                     for th, r in zip(ths, rs)]
            path = _catmull_rom_closed(inner + outer[::-1])
            if path is not None:
                ax.add_patch(PathPatch(path, fill=False, edgecolor=color,
                                       lw=_RIB_LW, zorder=1.6, alpha=0.95))
            r_max_used = max(r_max_used, max(rs) + _RIB_PAD_IN + _RIB_BAND)
            _bbox.extend(inner)
            _bbox.extend(outer)
        if biggest is not None:
            # 组名挂在**该段最外**的 tip 方向 (挂中点会压进网面)
            m = max(range(len(biggest)), key=lambda i: tip_r[biggest[i]])
            _labels.append((grp, color, (biggest[m] + 0.5) * sector,
                            tip_r[biggest[m]] + _RIB_PAD_IN + _RIB_BAND + 0.045))
        legend.append((grp, color, len(ks)))

    # ── 组名去重叠: 按半径从大到小落位, 冲突则沿径向逐步外推 ──
    #    (本数据的 7 个地理组在 CYCLE 上交错, 组名自然位置会挤在中心 —— 实测
    #     Wuwei / Beijing 直接压在一起)
    placed = []
    for grp, color, th, rr in sorted(_labels, key=lambda t: -t[3]):
        for _ in range(16):
            lx, ly = _pt(th, rr)
            if all((lx - px) ** 2 + (ly - py) ** 2 >= _LBL_MIN_D ** 2
                   for px, py in placed):
                break
            rr += 0.030
        placed.append((lx, ly))
        # 白描边 (path_effects.withStroke): 组名经常落在黑色网线之上, 没有描边
        # 就是"黑线穿字"。描边画在字形外侧, 不改变字色/字号/位置。
        ax.text(lx, ly, grp, ha='center', va='center', fontsize=8.0,
                color=color, fontweight='bold', zorder=3,
                path_effects=[pe.withStroke(linewidth=_LBL_STROKE,
                                            foreground='white')])
        _bbox.append((lx, ly))
        r_max_used = max(r_max_used, rr)
    return legend, r_max_used, _bbox


def render(nexus, out_pdf, out_png=None, label_every=None, title=None,
           groups=None, group_note=None, color_splits=False):
    # ── 样式真源: 显式钉死全部 rcParams, 避免顺序依赖 (与 Fig3 同一类缺陷) ──
    matplotlib.rcParams.update(sci_rcparams(base_size=8, grid=False,
                                            despine_top_right=False))
    _setup_font()

    labels, cycle, splits = parse_nexus(nexus)
    zones, n, sector = equal_angle_zones(labels, cycle, splits)
    total_w = max((z['d'] + z['w'] for z in zones), default=1.0)
    R = 1.0

    # ── 网状信号量化: 交错(互不相容)的 split 对数 ──
    _box_flags, n_box_pairs, n_pairs = flag_incompatible(zones, n, cycle)

    # ── tip 半径 = 该 taxon 最外层 zone 的外弦半径 (SplitsTree equal-angle 原意) ──
    #    ⚠ 旧版把所有 tip 强行拉到固定外圈 R*1.03 → 145 根等长灰辐条 + 中心空洞,
    #    看着像"太阳"而不是分裂网络。改用自然半径后外轮廓成为不规则凸多边形,
    #    与文献 NeighborNet 一致。
    tip_r = []
    for k, t in enumerate(cycle):
        outer = max((z for z in zones if t in z['taxa']),
                    key=lambda z: z['d'] + z['w'], default=None)
        tip_r.append((outer['d'] + outer['w']) / total_w * R if outer else 0.0)

    fig, ax = plt.subplots(figsize=(COL_W_2 + 0.6, COL_W_2 + 0.6),
                           facecolor='white')
    ax.set_facecolor('white')

    def pt(angle, r):
        return (r * math.cos(angle), r * math.sin(angle))

    # ── 分组弧线先画 (zorder 低于网络线; 弧在网络外圈, 不遮网面) ──
    legend_items, r_arc_max, _rib_bbox = ([], _TIP_LABEL_R, [])
    if groups:
        # ⚠ 第 7 轮发现 (P1, 静默丢标注)：`groups` 有**两种口径**
        #   · `load_group_map()` 返回 {**元数据名**: 组名}  ← 管线自己就是这么传的
        #     (phylo_pipeline.run_stage_splitstree: `_ss_render(..., groups=_groups)`)
        #   · `draw_group_ribbons()` 要 {**nexus 标签**: 组名}
        #   旧实现把前者直接透传 → `lab_pos.get(lab)` 全部落空 →
        #   **管线自产的 neighbornet.pdf 一条色带都没有**，且不报错。
        #   这里统一过一次 `_match_groups`（对已是标签口径的输入是幂等的：
        #   精确命中按归一化标签，回退分支不会触发）。
        _mapped = _match_groups(labels, groups)
        if _mapped:
            groups = _mapped
        try:
            legend_items, r_arc_max, _rib_bbox = draw_group_ribbons(
                ax, cycle, sector, labels, groups, tip_r)
        except Exception as e:  # noqa: BLE001 — 色带失败不拖垮网络本体
            print(f'[warn] 分组色带绘制失败 (仅出网络): {e}', file=sys.stderr)

    # zone 平行四边形: 内弦(d) / 外弦(d+w) / 两条径向边。
    # 颜色按 **split 权重** 着色 (viridis): 深色 = 强支持的分裂, 浅色 = 弱分裂。
    # (曾试过"把参与箱形结构的 split 标红" —— 实测 242 个 zone 里 241 个都参与过
    #  至少一对交错, 全红等于没标; 权重着色才是 NeighborNet 的标准读法。)
    from matplotlib import cm as _cm
    wmax = max((z['w'] for z in zones), default=1.0) or 1.0
    _zone_pts = []
    for z in zones:
        th1 = z['a'] * sector            # 区间起始射线
        th2 = (z['b'] + 1) * sector      # 区间结束射线
        d_in = z['d'] / total_w * R
        d_out = (z['d'] + z['w']) / total_w * R
        p1, p2 = pt(th1, d_in), pt(th2, d_in)
        p3, p4 = pt(th2, d_out), pt(th1, d_out)
        _zone_pts.extend([p1, p2, p3, p4])
        xs, ys = zip(p1, p2, p3, p4, p1)
        t = min(1.0, (z['w'] / wmax) ** 0.4)     # 幂律拉伸: 弱分裂仍可见
        col = _cm.viridis(0.12 + 0.88 * t) if color_splits else '#4A4A4A'
        ax.plot(xs, ys, color=col, lw=0.30 + 0.85 * t,
                alpha=0.92 if color_splits else 0.80, zorder=2,
                solid_capstyle='round')
    sm = None
    if color_splits:
        sm = _cm.ScalarMappable(
            cmap=_cm.viridis,
            norm=matplotlib.colors.PowerNorm(0.4, 0.0, wmax))
        sm.set_array([])

    # 分类单元: 落在**自然半径**处 (不再拉到统一外圈); 按组着色
    show_labels = (label_every is not None and label_every > 0) or \
                  (label_every is None and n <= _LABEL_AUTO_MAX)
    step = label_every if (label_every and label_every > 0) else 1
    grp_color = {g: c for (g, c, _k) in legend_items}
    lab_of_cycle = [labels[t - 1] if t - 1 < len(labels) else str(t)
                    for t in cycle]
    bbox = list(_rib_bbox) + list(_zone_pts)
    for k, t in enumerate(cycle):
        th = (k + 0.5) * sector
        r0 = tip_r[k]
        g = _short_group(groups[lab_of_cycle[k]]) \
            if (groups and lab_of_cycle[k] in groups) else None
        col = grp_color.get(g, '#B9B9B9')
        x, y = r0 * math.cos(th), r0 * math.sin(th)
        ax.plot([x], [y], marker='o', ms=1.9, mfc=col, mec='none', zorder=2.2)
        bbox.append((x, y))
        if show_labels and k % step == 0:
            lx, ly = pt(th, r0 + 0.028)
            ha = 'left' if math.cos(th) > 0.1 else ('right' if math.cos(th) < -0.1 else 'center')
            ax.text(lx, ly, lab_of_cycle[k], ha=ha, va='center',
                    fontsize=4.2, color='#333333')
            bbox.append((lx, ly))

    # ── 自适应画布: 按**实际几何**外接框铺满, 而不是套一个圆 ──
    #    旧版固定 lim=1.35, 而等角布局的外轮廓是不规则凸多边形 → 图缩在角落、
    #    四周大片留白 (实测仅占 ~40 % 画布)。
    bx = [p[0] for p in bbox] or [-1.0, 1.0]
    by = [p[1] for p in bbox] or [-1.0, 1.0]
    x0, x1, y0, y1 = min(bx), max(bx), min(by), max(by)
    pad = 0.075 * max(x1 - x0, y1 - y0, 1e-6)
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_aspect('equal')
    ax.axis('off')
    # 画布纵横比跟随数据外接框, 避免方画布 + 扁形状 → 两侧大留白
    _bw = (x1 - x0) + 2 * pad
    _bh = (y1 - y0) + 2 * pad
    _figw = COL_W_2 + 0.6
    _figh = max(_figw * (_bh / _bw) * (0.96 / 0.72), 1.15 / 0.24)
    fig.set_size_inches(_figw, _figh)
    ax.set_position([0.02, 0.24, 0.96, 0.72])
    if title:
        ax.set_title(title, fontsize=9.0, fontweight='bold', color='#222222')

    # ── 组别图例 (分组在 CYCLE 上交错, 弧线标签只能点最大的一段 → 必须有图例) ──
    if legend_items:
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], color=c, lw=2.2,
                          label='%s (n = %d)' % (g, k))
                   for (g, c, k) in legend_items]
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.5, 0.055), ncol=4, frameon=False,
                   fontsize=7.0, handlelength=1.3, columnspacing=1.4,
                   handletextpad=0.5)

    # ── 色标: 仅 --color-splits 时出 (默认统一深色线, 与文献 NeighborNet 一致) ──
    if sm is not None:
        cax = fig.add_axes([0.36, 0.030, 0.28, 0.010])
        cb = fig.colorbar(sm, cax=cax, orientation='horizontal')
        cb.ax.xaxis.set_label_position('top')
        cb.set_label('split weight (substitutions / site)',
                     fontsize=5.2, labelpad=1.5)
        from matplotlib import ticker as _ticker
        cb.locator = _ticker.MaxNLocator(nbins=3)
        cb.ax.tick_params(labelsize=4.6, length=1.5, pad=1.0)
        cb.update_ticks()
        cb.outline.set_linewidth(0.4)

    fit = None
    m = re.search(r'fit=([\d.]+)', open(nexus, errors='ignore').read())
    if m:
        fit = float(m.group(1))
    frac = (100.0 * n_box_pairs / n_pairs) if n_pairs else 0.0
    cap = 'Neighbor-Net (SplitsTree6): n = %d taxa, %d splits' % (n, len(zones))
    if fit is not None:
        cap += ', fit = %.1f%%' % fit
    cap2 = ('%s/%s split pairs are incompatible (%.1f%%)'
            % (format(n_box_pairs, ','), format(n_pairs, ','), frac))
    if legend_items:
        cap2 += '; %d %s groups (arcs and tip dots)' % (
            len(legend_items), group_note or 'meta')
    fig.text(0.5, 0.010, cap + '\n' + cap2, ha='center', va='bottom',
             fontsize=6.0, color='#333333', linespacing=1.55)
    fig.savefig(out_pdf, dpi=300, facecolor='white')
    if out_png:
        fig.savefig(out_png, dpi=300, facecolor='white')
    plt.close(fig)
    print('saved:', out_pdf, '(%d zones, %d interleaving pairs, %d group arcs)'
          % (len(zones), n_box_pairs, len(legend_items)))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', required=True)
    ap.add_argument('-o', '--pdf', required=True)
    ap.add_argument('--png', default=None)
    ap.add_argument('--label-every', type=int, default=None,
                    help='标注每 N 个 taxa; 不给则 n<=45 全标, 否则不标')
    ap.add_argument('--title', default=None)
    ap.add_argument('--meta', default=None,
                    help='元数据 CSV/TSV (含 name + 分组列); 给了才画分组弧线')
    ap.add_argument('--group-col', default='auto',
                    help="分组列: auto(默认, location→host 自动回退) / location / host / 列名")
    ap.add_argument('--color-splits', action='store_true',
                    help="按 split 权重用 viridis 着色 + 出色标 (默认统一深色线)")
    a = ap.parse_args()
    groups = None
    note = None
    if a.meta:
        raw, note = load_group_map(a.meta, a.group_col)
        if raw:
            labels, _c, _s = parse_nexus(a.input)
            groups = _match_groups(labels, raw)
            if not groups:
                print('[warn] 元数据与 nexus 标签无法匹配, 只画网络本体',
                      file=sys.stderr)
        else:
            print(f'[warn] 元数据无可用分组列 ({a.group_col}), 只画网络本体',
                  file=sys.stderr)
    render(a.input, a.pdf, a.png, a.label_every, a.title,
           groups=groups, group_note=note, color_splits=a.color_splits)
