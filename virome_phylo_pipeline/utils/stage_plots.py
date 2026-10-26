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

# ── 图件四格式导出 (SVG / PNG / PDF / PNG@1200dpi) — 默认开启 ──────────
# 2026-10-20: 所有出图脚本统一口径; 用 MMPV_FIG_EXPORT=0 可整体关闭,
# MMPV_FIG_FORMATS=svg,png,pdf,hires 可裁剪, MMPV_FIG_HIRES_DPI 可改高清 dpi。
try:
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
except ImportError:  # 仓库根不在 sys.path (脚本被单独执行) → 向上找 mmpv_common/
    import os as _os
    import sys as _sys
    _d = _os.path.dirname(_os.path.abspath(__file__))
    while _d and not _os.path.isdir(_os.path.join(_d, "mmpv_common")):
        _parent = _os.path.dirname(_d)
        if _parent == _d:
            break
        _d = _parent
    if _d and _d not in _sys.path:
        _sys.path.insert(0, _d)
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
_install_fig_hook()
# ──────────────────────────────────────────────────────────────────────
import os
import re
import sys
import math
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


def _norm_tip(s):
    """tip 名归一化 (与 splitstree_viz._norm_label 同口径: 小写 + 非字母数字→'_')。"""
    s = str(s).strip().lower()
    s = re.sub(r'[^a-z0-9]+', '_', s)
    return re.sub(r'_+', '_', s).strip('_')


def _find_meta_for_tree(output_dir):
    """从 phylogeny/ 目录反推可用的元数据表 (按优先级)。

    2026-10-13: 原只找 `data/sample_metadata.csv`。现同时找 `data/host.csv` ——
    两者都可能各自缺列（本数据的 `sample_metadata.csv` 只有 location 没有 host，
    `host.csv` 只有 host 没有 location），所以**分开查、各取所需**。
    2026-10-15: 候选目录布局感知（io_layout ph_dir('data'): standard=01_data /
    legacy=data）—— 原只查 legacy `data/`，standard 默认布局下树图 tip 的
    宿主/地点着色静默缺失。
    """
    root = os.path.dirname(os.path.abspath(output_dir))
    try:
        from mmpv_common.io_layout import ph_dir as _ph_dir
        _prefixes = []
        for _p in (_ph_dir('data'), "data"):
            if _p not in _prefixes:
                _prefixes.append(_p)
    except Exception:
        _prefixes = ["data"]
    out = {}
    for key, tail in (("meta", "sample_metadata.csv"), ("host", "host.csv")):
        for pref in _prefixes:
            for rel in (f"{pref}/{tail}",
                        f"{pref}/metadata/{tail}",
                        f"{pref}/clean/{tail}"):
                p = os.path.join(root, rel)
                if os.path.exists(p):
                    out[key] = p
                    break
            if key in out:
                break
    return out


def _load_tip_attr(meta_path, col_candidates, reduce=None):
    """元数据表 → {归一化 tip 名: 属性值}。

    col_candidates: 列名候选（大小写不敏感，取第一个命中的）。
    reduce: 值有多段逗号时如何取（如地点 "China, Ningxia, Yinchuan" → 末段）。
    """
    try:
        import csv as _csv
        with open(meta_path, encoding='utf-8-sig', newline='') as fh:
            rd = _csv.DictReader(fh)
            fn = rd.fieldnames or []
            cols = {str(c).strip().lower(): c for c in fn}
            name_c = cols.get('name') or (fn[0] if fn else None)
            val_c = next((cols[c] for c in col_candidates if c in cols), None)
            if not name_c or not val_c:
                return {}
            out = {}
            for row in rd:
                n = _norm_tip(row.get(name_c, ''))
                v = (row.get(val_c) or '').strip()
                if n and v:
                    out[n] = (reduce(v) if reduce else v)
            return out
    except Exception:  # noqa: BLE001
        return {}


def _load_location_map(meta_path):
    """元数据 → {归一化 tip 名: 地点标签}。地点取最后一段 (China, Ningxia, Yinchuan → Yinchuan)。"""
    return _load_tip_attr(meta_path, ('location',),
                          reduce=lambda v: v.split(',')[-1].strip() or v)


def _load_host_map(meta_path):
    """元数据 → {归一化 tip 名: 宿主物种名}。保留全名（如 Lycium barbarum）。

    ⚠ 宿主名是**双名法**，不能像地点那样取逗号末段 —— 复数宿主会写成
    "Lycium barbarum, Nicotiana benthamiana" 形式，此处取**首段**（主宿主）。
    """
    return _load_tip_attr(meta_path, ('host', 'host_species'),
                          reduce=lambda v: v.split(',')[0].strip() or v)


def _plot_tree_classic(tree_file, output_dir, log=None, meta=None,
                       show_tip_labels=True, dpi=300, branch_color="none"):
    """经典样式环形 ML 树（2026-10-14, 方案 A / ggtree·iTOL 惯例读法）。

    与彩色版 (`plot_tree_circular`) 的分工：
      · 彩色版 = Fig3 半栏嵌入用（无标签 cladogram + 双环色带）；
      · 本函数 = **独立大图**用：单色枝 + 等深 cladogram + 全部 tip 标签 +
        tip 端双色点环（内点 = 宿主物种 / 外点 = 采样地点）。

    branch_color 参数（2026-10-14 补, 用户要求「宿主/地理映射到树枝上」）：
      · "none"      (默认) 枝不带映射 —— 宿主/地理只在 tip 点上；
      · "host"      枝色 = 后代 tip 的**众数宿主**（≥90% 纯度才上色, 混合/未
                    注释枝 = 灰）; 地点点保留, 宿主点撤（避免双重编码）;
      · "location"  枝色 = 众数地点; 宿主点保留, 地点点撤。
      ⚠ 这是把旧彩版「内部枝无含义的统一蓝」修正为有统计含义的映射：
        枝色语义 = 「该枝后代以某组为主」, 图例注明; 不能按 tip 各自着色
        后直接当内部枝的颜色（那就是 2026-10-14 被判「混乱」的旧病）。

    ⚠ 2026-10-14 彩色版被用户判「颜色混乱」的三个根因，本样式逐一规避：
      ① 彩色版内部骨架是统一蓝（`OKABE_ITO[4]`），与图例 *L. chinense* 的
         浅蓝同族 → 读者误读内部枝的宿主含义；本样式**通体单色**，颜色只
         出现在 tip 点上（语义 = 图例，无一例外）。
      ② 彩色版内外双环 + 分支共用 OKABE → 宿主与地点撞色；本样式宿主用
         OKABE、地点用 viridis 冷色系（色相族不同）。
      ③ 无 tip 标签无法锚定；本样式全部 127 个 accession 可读。

    布局 = **等深 cladogram**（每层等宽、tip 一律拉到最外圈）——与用户提供的
    参考图同款：127 个标签在统一半径成环, 相邻标签弧间距 ≈ 12 pt > 字高,
    不叠字。⚠ cladogram 不表达枝长（本树枝长跨度 52.8×, 真实比例会把 120 条
    tip 挤成中团）—— 尺度信息在 treefile/.iqtree 里, 注记中写明。
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None
    from matplotlib.lines import Line2D

    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    try:
        txt = open(tree_file).read()
        # 去掉 bootstrap 支持值标签保留分支长（与彩色版同口径）
        txt = re.sub(r'\)\s*\d+(\.\d+)?\s*([,):])', r')\2', txt)
        root = _parse_newick(txt)
    except Exception as e:
        log.emit(f"[tree_classic] 解析失败: {e}")
        return None

    # ── 零长枝塌缩（与彩色版同算法同口径 1e-6, 只塌内部枝）──
    _n_collapsed = [0]

    def _collapse_zero(n):
        for _c in list(n['children']):
            _collapse_zero(_c)
        _changed = True
        while _changed:
            _changed = False
            _new = []
            for _c in n['children']:
                if _c['children'] and (_c.get('length') or 0.0) <= 1e-6:
                    _new.extend(_c['children'])
                    _n_collapsed[0] += 1
                    _changed = True
                else:
                    _new.append(_c)
            n['children'] = _new

    try:
        _collapse_zero(root)
    except Exception:                                     # noqa: BLE001
        _n_collapsed[0] = 0

    # ── 等深 cladogram 布局: x = 层级深度, tip 一律拉到最外圈 ──
    leaves = []

    def assign(n, depth):
        n['_nd'] = depth
        if not n['children']:
            n['y'] = len(leaves)
            leaves.append(n)
        else:
            for c in n['children']:
                assign(c, depth + 1)
            n['y'] = sum(c['y'] for c in n['children']) / len(n['children'])

    assign(root, 0)
    _max_nd = max(n['_nd'] for n in leaves) or 1

    def _setx(n):
        n['x'] = float(n['_nd'])          # 等深: 每层等宽
        for c in n['children']:
            _setx(c)
    _setx(root)
    for _lf in leaves:
        _lf['x'] = float(_max_nd)         # tip 一律拉到最外圈
    max_x = float(_max_nd)
    n_tips = len(leaves)
    if n_tips < 2:
        log.emit("[tree_classic] tips < 2, 跳过")
        return None

    # ── 元数据: 宿主 (OKABE) + 地点 (Dark2, 与宿主无撞色) ──
    paths = _find_meta_for_tree(output_dir)
    if meta:
        paths = dict(paths)
        paths["host"] = meta
    host_map = _load_host_map(paths.get("host")) or {}
    loc_map = _load_location_map(paths.get("meta")) or {}

    def _attr_of(amap):
        out = {}
        for lf in leaves:
            k = _norm_tip(lf.get('label', ''))
            v = amap.get(k)
            if v is None and '_' in k:
                h = k.split('_')[0]
                v = next((vv for kk, vv in amap.items()
                          if kk.split('_')[0] == h), None)
            if v:
                out[id(lf)] = v
        return out

    try:
        from utils.sci_palette import apply_sci_style, OKABE_ITO
    except ImportError:                                   # 脚本方式直跑
        from sci_palette import apply_sci_style, OKABE_ITO
    host_of = _attr_of(host_map)
    loc_of = _attr_of(loc_map)
    host_groups = sorted(set(host_of.values()))
    loc_groups = sorted(set(loc_of.values()))
    host_color = {g: OKABE_ITO[i] for i, g in enumerate(host_groups)}
    # 地点用 viridis 冷色系采样 —— 与宿主 OKABE（暖橙/浅蓝/绿）无色相碰撞:
    # 2026-10-14 首版用 Dark2, 其橙色 (#D95F02) 与宿主 *L. barbarum* 的
    # #E69F00 同族, 双环同图仍有撞色风险。宿主暖色 / 地点冷色, 语义也顺。
    try:
        _loc_cmap = plt.get_cmap('viridis')
    except AttributeError:                                # matplotlib >= 3.9
        import matplotlib as _mpl
        _loc_cmap = _mpl.colormaps['viridis']
    loc_color = {_g: _loc_cmap(i / max(1, len(loc_groups) - 1))
                 for i, _g in enumerate(loc_groups)}

    # ── 枝色映射 (branch_color): 自底向上众数继承, ≥90% 纯度才上色 ──
    _GREY = '#BBBBBB'
    _AMAP = {"host": host_of, "location": loc_of}.get(branch_color)
    _PALETTE = {"host": host_color, "location": loc_color}.get(branch_color)

    def _majority(amap, threshold=0.9):
        """{id(node): 组名} —— 节点子树内有标注 tip 的众数组, 纯度≥threshold 才记录。"""
        out = {}
        import collections as _co

        def visit(n):
            counts = _co.Counter()
            if not n['children']:
                g = amap.get(id(n))
                if g:
                    counts[g] += 1
            else:
                for c in n['children']:
                    counts.update(visit(c))
            if counts:
                g, k = counts.most_common(1)[0]
                if k / sum(counts.values()) >= threshold:
                    out[id(n)] = g
            return counts

        visit(root)
        return out

    _branch_map = _majority(_AMAP) if (_AMAP and _PALETTE) else None

    # ── 画布: 独立大图 ──
    apply_sci_style(base_size=7, grid=False)
    fig = plt.figure(figsize=(9.6, 9.6), dpi=dpi)
    ax = fig.add_subplot(111, projection='polar')
    two_pi = 2 * np.pi
    phase = np.pi / 2                      # 12 点钟起步, 顺时针

    def theta(y):
        return phase - two_pi * y / n_tips

    # 半径预算（数据单位 = 层级）: tip 环 → 宿主点 → 地点点 → 标签
    _U = max_x * 0.01
    _r_dot_in = max_x + _U * 1.2           # 宿主点环
    _r_dot_out = max_x + _U * 3.0          # 地点点环
    _r_lab = max_x + _U * 5.2              # 标签统一半径

    _C_TREE = '#4D4D4D'
    stack = [root]
    while stack:
        n = stack.pop()
        for c in n['children']:
            if _branch_map is not None:
                _g = _branch_map.get(id(c))
                _cc = _PALETTE.get(_g, _GREY) if _g else _GREY
            else:
                _cc = _C_TREE
            ax.plot([theta(c['y']), theta(c['y'])], [n['x'], c['x']],
                    color=_cc, lw=0.8, solid_capstyle='round', zorder=2)
            stack.append(c)
        if n['children']:
            ys = [c['y'] for c in n['children']]
            t0, t1 = theta(min(ys)), theta(max(ys))
            k = max(2, int(abs(t1 - t0) * 160))
            # 弧段 = 节点自身的枝色 (众数继承), 混合节点灰 —— clade 色块连贯
            if _branch_map is not None:
                _g = _branch_map.get(id(n))
                _cc = _PALETTE.get(_g, _GREY) if _g else _GREY
            else:
                _cc = _C_TREE
            ax.plot(np.linspace(t0, t1, k), [n['x']] * k,
                    color=_cc, lw=0.8, solid_capstyle='round', zorder=2)

    # ── tip 点环 + 径向标签 ──
    # 枝色已表达某维度时, 该维度的点撤掉 (双重编码 = 又回到「混乱」老路)
    _show_host_dots = branch_color != "host"
    _show_loc_dots = branch_color != "location"
    for lf in leaves:
        th = theta(lf['y'])
        if _show_host_dots and host_of.get(id(lf)):
            ax.plot([th], [_r_dot_in], 'o', ms=2.6,
                    color=host_color[host_of[id(lf)]], zorder=4)
        if _show_loc_dots and loc_of.get(id(lf)):
            ax.plot([th], [_r_dot_out], 'o', ms=2.6,
                    color=loc_color[loc_of[id(lf)]], zorder=4)
        if show_tip_labels and lf.get('label'):
            # 径向文字: 旋转 = 方位角本身, 左半圈翻 180° + 锚点换边 → 字头朝外
            _deg = np.degrees(th) % 360.0
            if 90.0 < _deg < 270.0:
                _rot, _ha = _deg - 180.0, 'right'
            else:
                _rot, _ha = _deg, 'left'
            ax.text(th, _r_lab, lf['label'], rotation=_rot,
                    rotation_mode='anchor', ha=_ha, va='center',
                    fontsize=5.8, color='#222222', zorder=5)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines['polar'].set_visible(False)
    # 标签区 ≈ 最长标签(22 字符 × 0.6 × 5.8pt) ≈ 1.06 in; 轴半径 ~3.9 in
    _lab_chars = max((len(lf.get('label') or '') for lf in leaves), default=0)
    _lab_in = _lab_chars * 0.60 * 5.8 / 72.0
    _r_lim = _r_lab + _lab_in / 3.9 * _r_lab
    ax.set_ylim(0, _r_lim * 1.06)
    ax.set_title(f'Maximum-likelihood tree  (n = {n_tips} tips)',
                 fontsize=9, fontweight='bold', pad=14)

    # ── 图例: 枝色图例 (branch_color 时) + 点图例, 左边距分栏 ──
    handles = []
    _n_branch = 0
    if _branch_map is not None:
        # 枝色图例只列**实际出现在枝上的组** (被灰枝吞掉的少数组不列)
        _used = sorted({_g for _g in _branch_map.values()})
        handles += [Line2D([], [], color=_PALETTE[g], lw=1.6, label=g)
                    for g in _used]
        _n_branch = len(_used)
    if _show_host_dots:
        handles += [Line2D([], [], linestyle='none', marker='o', ms=4,
                           color=host_color[g], label=g) for g in host_groups]
    if _show_loc_dots:
        handles += [Line2D([], [], linestyle='none', marker='o', ms=4,
                           color=loc_color[g], label=g) for g in loc_groups]
    if not handles:
        return None
    _n_h = _n_branch if _branch_map is not None else max(1, len(host_groups))
    # 图例小标题: 枝色 / 点色各是什么
    _t1 = ({'host': 'Host plant (branch)',
            'location': 'Sampling location (branch)'}.get(branch_color,
                                                           'Host plant (dot)')
           if handles else '')
    leg1 = ax.legend(handles=handles[:_n_h],
                     loc='upper left', bbox_to_anchor=(-0.145, 1.02),
                     fontsize=6.5, frameon=True, framealpha=0.9,
                     title=_t1, title_fontsize=6.8,
                     labelspacing=0.32, borderpad=0.5)
    ax.add_artist(leg1)
    _rest = handles[_n_h:]
    if _rest:
        _t2 = ('Sampling location (dot)' if branch_color == "host"
               else 'Host plant (dot)' if branch_color == "location"
               else 'Sampling location (dot)')
        ax.legend(handles=_rest, loc='upper left',
                  bbox_to_anchor=(-0.145, 0.76),
                  fontsize=6.5, frameon=True, framealpha=0.9,
                  title=_t2, title_fontsize=6.8,
                  labelspacing=0.32, borderpad=0.5, ncol=1)

    _note_l1 = ('cladogram (branch lengths not to scale; topology from '
                'IQ-TREE ML, see treefile/.iqtree for substitutions per site)')
    if _branch_map is not None:
        _dim_lab = 'host plant' if branch_color == "host" else 'sampling location'
        _note_l2 = (f'branch colour = majority {_dim_lab} among descendant tips '
                    '(≥90% purity; grey = mixed or unassigned)'
                    + (f'; {len(host_of)}/{n_tips} tips with host'
                       if branch_color == "location" else '')
                    + (f'; {_n_collapsed[0]} zero-length internal branches collapsed'
                       if _n_collapsed[0] else ''))
    else:
        _note_l2 = ('dots: inner = host plant, outer = sampling location'
                    f'; {len(host_of)}/{n_tips} tips with host, '
                    f'{len(loc_of)}/{n_tips} with location'
                    + (f'; {_n_collapsed[0]} zero-length internal branches collapsed'
                       if _n_collapsed[0] else ''))
    fig.text(0.012, 0.006, _note_l1 + '\n' + _note_l2, fontsize=6.0,
             color='#666666', ha='left', va='bottom', linespacing=1.4)
    # 左侧 15% 边距给图例（不用 tight_layout, 它会把极轴撑回满幅）
    fig.subplots_adjust(left=0.155, right=0.995, top=0.945, bottom=0.055)
    _suffix = {'host': '_brhost', 'location': '_brloc'}.get(branch_color, '')
    _finish(fig, os.path.join(output_dir, f'tree_circular{_suffix}.pdf'),
            os.path.join(output_dir, f'tree_circular{_suffix}.png'))
    log.emit(f"[tree_classic] {n_tips} tips, host {len(host_of)}, "
             f"loc {len(loc_of)}, branch_color={branch_color} "
             f"-> {os.path.join(output_dir, f'tree_circular{_suffix}.png')}")
    return os.path.join(output_dir, f'tree_circular{_suffix}.png')


def plot_tree_circular(tree_file, output_dir, log=None, meta=None, color_by="host",
                       show_branch_length=False, collapse_zero_len=True,
                       zero_len_tol=1e-6, style="colored", branch_color="none"):
    """环形系统发育图 —— 真正的树结构 (径向分支 + 等半径弧)。

    2026-10-12 重写: 旧实现只画 `[0, theta] → [0, depth]` 的**放射直线**,
    等于把 127 条 tip 画成"星形图", 既没有内部节点也没有分支弧, 不是树;
    且 tip 标签只有 3.2 pt (排进 4.58 in 栏宽后 ~2 pt), 属于视觉噪声。
    现在按标准矩形系统发育图 → 极坐标包裹: 每个节点给 (x=根到节点距离, y=叶序),
    分支画径向段、内部节点画等半径弧。

    2026-10-12 再补: 末端分支着色。本数据 127 条 tip 里 120 条枝长极短,
    单色径向图看上去就是一团墨 —— 必须把某个元数据维度编码进颜色。

    2026-10-13 改口径: 着色维度由**采样地点**改为**宿主植物**（用户指定
    「分支展示不同物种的，不按地区」）。`color_by` 可选：
      · "host"     (默认) → 宿主物种，如 Lycium barbarum / L. ruthenicum / L. chinense
      · "location"        → 采样地点（旧口径，保留以便回退 / 对拍）
      · "none"            → 单色
    数据源按 `color_by` 各自反推：host → `data/host.csv`；location →
    `data/sample_metadata.csv`。对应表缺失时**退回单色**（不报错）。

    2026-10-13 定稿: **默认不展示枝长（cladogram）** —— 本树枝长跨度 52.8×，
    真实枝长会把 120 条 tip 挤在圆心成一团墨；cladogram（层级等距）拓扑清晰得多。
    外圈保留**双环色带**（内环 = 主维度 / 外环 = 另一维度），参考图那种
    「沿圆周分段上色 + 组名」的排法。
    `show_branch_length=True` 可回到真实枝长 + `BL_FLOOR_FRAC=0.06` 下限截断。

    2026-10-13 再改（用户要求「枝长为 0 的去除」）: **零长枝塌缩**
    （`collapse_zero_len=True`, `zero_len_tol=1e-6`）。
    ⚠ 本树 IQ-TREE 输出里**没有字面 `:0`** —— 最短枝长恒为 **9.98e-07**，
      而该值出现 **43 次**（其余所有值都只出现 1 次）⇒ 这是 IQ-TREE 的
      **数值下界**，即这 43 条枝的真实枝长在数值上等同于 0
      （18 条内部枝 + 25 条 tip 枝）。
    ⚠ cladogram 模式下 **tip 枝长不进入图形**（所有 tip 一律拉到最外圈），
      所以真正影响形态的只有 **内部零长枝** —— 它们各自白占一个层级。
      实测塌缩 18 条内部零长枝 → 最大深度 25→23、平均深度 14.2→12.4、
      层级数 26→24，树更紧凑且分支不再被"零长度分化"稀释。
    ⚠ **只塌缩内部节点**：塌缩 tip 等于删掉序列，绝对不行。
    ⚠ 枝长语义：被塌缩节点的 children 提升到祖父，其枝长**不变**
      （因为被合并的那段距离 ~0）。
    ⚠⚠ 实现必须**自底向上递归 + 本层循环到稳定**：早期版本用「walk 快照 +
      逐节点提升」，被提升走的节点仍留在快照里、其孙子被重复检查 ⇒ 移除数
      虚高（实测把 18 条算成 26）。此处每个节点至多被移除一次。

    style 参数（2026-10-14）: "colored"（默认, Fig3 半栏嵌版, 双环+分支着色）
      | "classic"（独立大图: 单色枝+真实枝长+全标签+tip 双色点环, 见
      `_plot_tree_classic` —— 用户 2026-10-14 判彩色版「颜色混乱」后按参考
      图样式新增）。
    """
    if style == "classic":
        return _plot_tree_classic(tree_file, output_dir, log=log,
                                  meta=meta, branch_color=branch_color)
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

    # ── 零长枝塌缩（2026-10-13 用户要求：枝长为 0 的去除后重绘）──
    # ⚠⚠ 本树 IQ-TREE 输出里**没有字面 `:0`** —— 最短枝长恒为 **9.98e-07**，
    #   且该值出现 **43 次**（其余所有值都只出现 1 次）⇒ 这是 IQ-TREE 的
    #   **数值下界**：这 43 条枝的真实枝长在数值上等同于 0
    #   （实测 18 条内部枝 + 25 条 tip 枝）。
    # ⚠ cladogram 模式下 **tip 枝长不进入图形**（所有 tip 一律拉到最外圈），
    #   所以真正影响形态的只有 **内部零长枝** —— 它们各自白占一个层级。
    # ⚠ **只塌缩内部节点**：塌缩 tip 等于删掉序列，绝对不行。
    # ⚠ 枝长语义：被塌缩节点的 children 提升到祖父，其枝长**不变**
    #   （被合并掉的那段距离 ~0）。
    _n_collapsed = [0]

    def _collapse_zero(n):
        """自底向上塌缩枝长 <= zero_len_tol 的内部枝（每个节点至多移除一次）。"""
        for _c in list(n['children']):        # 快照：本层随后会改写 children
            _collapse_zero(_c)                # 先把子树塌缩干净
        _changed = True
        while _changed:                       # 本层循环到稳定
            _changed = False
            _new = []
            for _c in n['children']:
                if _c['children'] and (_c.get('length') or 0.0) <= zero_len_tol:
                    _new.extend(_c['children'])    # 提升孙子，_c 被移除
                    _n_collapsed[0] += 1
                    _changed = True
                else:
                    _new.append(_c)
            n['children'] = _new

    if collapse_zero_len:
        try:
            _collapse_zero(root)
            if _n_collapsed[0]:
                (log or LogCollector()).emit(
                    f"[tree_plot] 塌缩零长枝 {_n_collapsed[0]} 条 "
                    f"(内部枝长 <= {zero_len_tol:g}；IQ-TREE 数值下界 9.98e-07)")
        except Exception as _ce:              # noqa: BLE001
            (log or LogCollector()).emit(f"[tree_plot] 零长枝塌缩失败(已忽略): {_ce}")
            _n_collapsed[0] = 0

    # ── 布局 ──
    leaves = []
    # ⚠ 2026-10-13 定稿：**默认 cladogram（不展示枝长）**。枝长跨度实测 52.8×
    #   （中位 2.7e-4 / max 1.4e-2）→ 真实枝长会把 120 条 tip 全挤在圆心成一团墨。
    #   `show_branch_length=True` 时回退到真实枝长 + 下限截断（每条枝长抬到
    #   ≥`BL_FLOOR_FRAC` × 最长枝），既保相对比例、又给极短枝留出可见长度。
    BL_FLOOR_FRAC = 0.06
    if not show_branch_length:
        # cladogram（不展示枝长）
        # ⚠⚠ 2026-10-13 二次修正：**原实现 `x = 层级深度(整数)` 在环形布局下会
        #   得到「同心迷宫」**，用户反馈「树枝去哪了」。根因（实测本数据）：
        #   127 tips 的层级深度 1–25、均值 14.2，且 251 条枝里 89 条 ≤1e-4、
        #   43 条 ≤1e-6 —— 即**大量近零枝却贡献了真实层级**。于是：
        #     ① tip 各自停在自己的深度半径上（深度 1 的贴圆心、深度 25 的在外圈）
        #        → 内圈孤立点 + 外圈一圈叠一圈的弧 → 看起来像靶纸/迷宫；
        #     ② 每层占等宽半径 → 25 层把径向空间切得太碎。
        #   ✅ 正解（=「等长 tip cladogram」业界标准做法）：
        #     ① **所有 tip 一律拉到最外圈**（等长 tip 段）—— 这是环状树之所以
        #        「像树」的关键：分支密度由**弧长**表达，径向只表达层级；
        #     ② **层级做 sqrt 压缩**（`ND_COMPRESS="sqrt"`）→ 内圈不再被前几层
        #        独占，整棵树由内向外均匀展开。
        #   ⚠ 若想回到旧行为（同心环）：`ND_COMPRESS="raw"`。
        ND_COMPRESS = "sqrt"          # {"sqrt" | "raw" | "log1p"}
        def assign(n, depth):
            n['_nd'] = depth
            if not n['children']:
                n['y'] = len(leaves)
                leaves.append(n)
            else:
                for c in n['children']:
                    assign(c, depth + 1)
                n['y'] = sum(c['y'] for c in n['children']) / len(n['children'])
        assign(root, 0)
        _max_nd = max(n['_nd'] for n in leaves) or 1
        if ND_COMPRESS == "raw":
            _f = lambda d: float(d)
        elif ND_COMPRESS == "log1p":
            _k = 1.6
            _f = lambda d: (math.log1p(d * _k) / math.log1p(_max_nd * _k)
                            * _max_nd) if _max_nd else 0.0
        else:                                              # sqrt（默认）
            _f = lambda d: (math.sqrt(d) / math.sqrt(_max_nd)
                            * _max_nd) if _max_nd else 0.0

        def _setx(n):
            n['x'] = _f(n['_nd'])
            for c in n['children']:
                _setx(c)
        _setx(root)
        for _lf in leaves:                                 # ⭐ tip 一律拉到最外圈
            _lf['x'] = float(_max_nd)
        max_x = float(_max_nd) or 1.0
        _radial_note = ('branch lengths not shown (cladogram; '
                        + ('levels sqrt-compressed' if ND_COMPRESS == 'sqrt'
                           else 'levels log1p-compressed' if ND_COMPRESS == 'log1p'
                           else 'equal depth per level')
                        + (f'; {_n_collapsed[0]} zero-length branches collapsed'
                           if _n_collapsed[0] else '')
                        + ')')
    else:
        # 真实枝长 + 下限截断
        _all_bl = []
        _stack = [root]
        while _stack:
            _n = _stack.pop()
            for _c in _n['children']:
                _all_bl.append(_c.get('length', 0.0) or 0.0)
                _stack.append(_c)
        _max_bl = max(_all_bl) if _all_bl else 1.0
        _floor = BL_FLOOR_FRAC * _max_bl

        def assign(n, x):
            x = x + max(n['length'], _floor)
            n['x'] = x
            if not n['children']:
                n['y'] = len(leaves)
                leaves.append(n)
            else:
                for c in n['children']:
                    assign(c, x)
                n['y'] = sum(c['y'] for c in n['children']) / len(n['children'])
        assign(root, 0.0)
        max_x = max(n['x'] for n in leaves) or 1.0
        _radial_note = (f'branch length 0 \u2013 {_max_bl:.4f} subs/site '
                        f'(min clipped to {BL_FLOOR_FRAC:.0%} of max for legibility)'
                        + (f'; {_n_collapsed[0]} zero-length branches collapsed'
                           if _n_collapsed[0] else ''))

    n_tips = len(leaves)
    if n_tips < 2:
        return None

    # ── tip 按组排序（2026-10-13 关键修复）──
    # ⚠ 首版双环色带的组名在一圈里重复 4~5 次：因为 Newick 的叶序是拓扑序
    #   （跟分组无关），同一组被切成十几段碎弧。参考图的做法是**先把 tip 按组
    #   排好**，环带自然收敛成"一组一段"的大弧。这里在 y 赋值完成后做一次
    #   叶序置换：按 (主维度组, 副维度组) 排序 → 重排 y → 再回写所有内部节点的 y。
    _sort_key = None
    _paths_pre = _find_meta_for_tree(output_dir)
    if meta:                                   # 与渲染段一致：显式路径优先
        _paths_pre = dict(_paths_pre)
        _paths_pre["host" if color_by == "host" else "meta"] = meta
    try:
        _hmap = _load_host_map(_paths_pre.get("host")) or {}
        _lmap = _load_location_map(_paths_pre.get("meta")) or {}
        if _hmap or _lmap:
            def _getv(amap, lf):
                k = _norm_tip(lf.get('label', ''))
                v = amap.get(k)
                if v is None and '_' in k:
                    h = k.split('_')[0]
                    v = next((vv for kk, vv in amap.items()
                              if kk.split('_')[0] == h), None)
                return v or '~'          # '~' 排最后（未分配的 tip）
            if color_by == "host":
                _sort_key = lambda lf: (_getv(_hmap, lf), _getv(_lmap, lf))
            else:
                _sort_key = lambda lf: (_getv(_lmap, lf), _getv(_hmap, lf))
    except Exception:
        _sort_key = None
    if _sort_key is not None:
        try:
            leaves.sort(key=_sort_key)
            for _i, _lf in enumerate(leaves):            # 重排叶序
                _lf['y'] = _i
            _st2 = [root]
            while _st2:                                   # 自底向上回写内部节点 y
                _n = _st2.pop()
                if _n['children'] and _n not in leaves:
                    _n['y'] = sum(c['y'] for c in _n['children']) / len(_n['children'])
                _st2.extend(_n['children'])
        except Exception:
            pass

    os.makedirs(output_dir, exist_ok=True)
    from utils.sci_palette import apply_sci_style, COL_W_2, OKABE_ITO
    apply_sci_style(base_size=7, grid=False)
    fig, ax = plt.subplots(figsize=(COL_W_2, COL_W_2), dpi=300,
                           subplot_kw={'projection': 'polar'})

    two_pi = 2 * np.pi
    phase = np.pi / 2          # 从 12 点钟方向开始

    def theta(y):
        return phase - two_pi * y / n_tips

    # ── 着色: tip → 元数据标签（2026-10-13 起默认 = 宿主，不再是地区）──
    paths = _find_meta_for_tree(output_dir)
    if meta:                       # 显式传入的路径按 color_by 归位
        paths = dict(paths)
        paths["host" if color_by == "host" else "meta"] = meta

    _DIM = {
        "host":     (paths.get("host"), _load_host_map,
                     "host plant", "host plant"),
        "location": (paths.get("meta"), _load_location_map,
                     "sampling location", "sampling location"),
    }

    def _build_attr(dim):
        """dim -> ({id(leaf): 组名}, [组名...], {组名: 颜色}, 说明标签, 短名)"""
        src, loader, label, short = _DIM.get(dim, (None, None, None, None))
        amap = loader(src) if (src and loader) else {}
        if not amap:
            return {}, [], {}, None, None
        aof = {}
        for lf in leaves:
            k = _norm_tip(lf.get('label', ''))
            g = amap.get(k)
            if g is None and '_' in k:
                h = k.split('_')[0]
                g = next((v for kk, v in amap.items()
                          if kk.split('_')[0] == h), None)
            if g:
                aof[id(lf)] = g
        gs = sorted(set(aof.values()))
        # ⚠ 组数多时色板循环会导致两组同色 → 超过 OKABE_ITO 长度就退回 viridis
        if len(gs) <= len(OKABE_ITO):
            gc = {g: OKABE_ITO[i] for i, g in enumerate(gs)}
        else:
            from matplotlib import cm as _cm
            gc = {g: _cm.get_cmap('viridis')(
                i / max(1, len(gs) - 1)) for i, g in enumerate(gs)}
        return aof, gs, gc, label, short

    # 主维度（决定分支颜色 / 图例 / 标题文案）
    attr_of, groups_sorted, grp_color, _label, _short = _build_attr(color_by)
    # 环带维度：外圈双环 = 内环 location + 外环 host（2026-10-13 用户指定）
    _ring_outer = "location" if color_by == "host" else "host"
    ring_in, ring_in_g, ring_in_c, ring_in_lab, _ = _build_attr(color_by)
    ring_out, ring_out_g, ring_out_c, ring_out_lab, _ = _build_attr(_ring_outer)

    col = OKABE_ITO[4]
    stack = [root]
    while stack:
        n = stack.pop()
        for c in n['children']:
            # 径向分支段: 末端分支按元数据维度着色, 内部一律用统一的蓝
            g = attr_of.get(id(c))
            cc = grp_color.get(g, col) if g else col
            lw = 0.75 if g else 0.45
            ax.plot([theta(c['y']), theta(c['y'])], [n['x'], c['x']],
                    color=cc, lw=lw, solid_capstyle='round', zorder=3)
            stack.append(c)
        if n['children']:
            # 等半径弧 (内部节点的"横线")
            ys = [c['y'] for c in n['children']]
            t0, t1 = theta(min(ys)), theta(max(ys))
            k = max(2, int(abs(t1 - t0) * 140))
            ax.plot(np.linspace(t0, t1, k), [n['x']] * k,
                    color=col, lw=0.45, solid_capstyle='round', zorder=3)

    # 淡参考环
    for frac in (0.25, 0.50, 0.75, 1.00):
        ax.plot(np.linspace(0, two_pi, 240), [max_x * frac] * 240,
                color='#DDDDDD', lw=0.4, zorder=1)

    # ── 外圈双环色带（2026-10-13 用户指定：内环宿主 / 外环地区）──
    # 参考图那种「沿圆周分段上色 + 组名」的排法。环带画在**内容实际外缘**之外。
    # ⚠ 2026-10-13 修的 4 个坑（首版双环实测）：
    #   ① 环带半径不能用 `max_x*1.06` —— `max_x` 由**单条最长枝**决定（本树
    #      最长枝 1.4e-2 而中位仅 2.7e-4），绝大多数 tip 半径远小于它 → 环带与
    #      树之间出现巨大空白。改用**叶半径 95 分位**做几何基准（长枝自然穿带，
    #      视觉上就是"超长枝"标记，参考图同样如此）。
    #   ② 组名**同一个组在同一环内只标一次** —— 首版按"相邻同组"分段，而 tip
    #      顺序未按组分块 → 同一组碎成 4~5 段 → 组名重复 4~5 次。现按组聚合。
    #   ③ 标签翻转判据写错（`not (-90 <= _rot <= 90)` 漏掉 ±180 边界）→ 左下/右下
    #      象限出现倒字。改用标准写法 + `va` 随内外环反向，保证字头永远朝外。
    #   ④ 环带外还要再留出组名的位置，否则文字被画布裁掉（`bbox_inches='tight'`
    #      只撑画布、不救被轴裁掉的 text）。
    _r_leaf = [lf['x'] for lf in leaves if lf.get('x')]
    # cladogram 下所有 tip 半径几乎相同（=max_x）；真实枝长下用 95 分位，
    # 否则单条超长枝会把环带推到很远（首版 4.7 倍空白就是这原因）。
    _r_geo = (max_x if not show_branch_length
              else (float(np.percentile(_r_leaf, 95)) if _r_leaf else max_x))
    _RING_W = _r_geo * 0.050         # 单条环带厚度
    _RING_GAP = _r_geo * 0.030       # 两环间距（⚠ 太小会让内外环视觉上糊成一条）
    _r_base = _r_geo * 1.02          # 内环内缘半径（贴住树的最外沿）
    _ring_labels = []                # 收集 (组名, 颜色, 弧中点角, 环半径, 环序)
    for _ri, (_rmap, _rg, _rc, _rlab) in enumerate(
            [(ring_in, ring_in_g, ring_in_c, ring_in_lab),
             (ring_out, ring_out_g, ring_out_c, ring_out_lab)]):
        if not _rmap:
            continue
        r_in = _r_base + _ri * (_RING_W + _RING_GAP)
        # 按叶序扫描，把连续的同一组归成弧段（画的是连续弧段，不是逐 tip 小块）
        order = [_rmap.get(id(lf)) for lf in leaves]
        # ⚠ 环带画成**整圈底色**再叠各组的彩色段 → 未分配（None）的 tip 显示为
        #   浅灰底，不会在环上开天窗（参考图里无数据区同样是浅灰）。
        _RING_LW = 2.6
        ax.plot(np.linspace(0, two_pi, 361), [r_in] * 361,
                color='#E8E8E8', lw=_RING_LW,
                solid_capstyle='butt', zorder=2.4)
        segs = []
        _s = 0
        for _i in range(1, n_tips + 1):
            if _i == n_tips or order[_i] != order[_s]:
                if order[_s]:
                    segs.append((order[_s], _s, _i))
                _s = _i
        # 首尾同组时合并（环是闭合的）
        if len(segs) > 1 and segs[0][0] == segs[-1][0]:
            g0, _s0, _e0 = segs.pop(0)
            _s1, e1 = segs[-1][1], segs[-1][2]
            segs[-1] = (g0, _s1, e1)
        _seg_len = {}
        for _grp, _a, _b in segs:
            # 极角是**顺时针**走的（theta 随 y 递减）→ 取 min/max 保证弧段连续
            t0, t1 = theta(_a), theta(_b)
            _lo, _hi = (t0, t1) if t0 <= t1 else (t1, t0)
            _kk = max(3, int((_hi - _lo) * 160))
            ax.plot(np.linspace(_lo, _hi, _kk), [r_in] * _kk,
                    color=_rc[_grp], lw=2.6, solid_capstyle='butt',
                    zorder=2.5, alpha=0.9)
            _seg_len[_grp] = _seg_len.get(_grp, 0.0) + (_hi - _lo)
        # ② 每个组只取**最长的那一段**标字 → 同名组在一圈内最多出现一次
        for _grp, _a, _b in segs:
            t0, t1 = theta(_a), theta(_b)
            _lo, _hi = (t0, t1) if t0 <= t1 else (t1, t0)
            if (_hi - _lo) < _seg_len.get(_grp, 0.0) - 1e-9:
                continue                    # 不是该组最长的段 → 不标字
            if (_hi - _lo) < 0.06:          # 太窄的字放不下 → 不标，避免叠字
                continue
            _ring_labels.append((_grp, _rc[_grp], (_lo + _hi) / 2.0,
                                 r_in + _RING_W * 0.5, _ri, _hi - _lo))

    # 环带组名：沿切向写；标准翻转判据 + `va` 随内外环反向 → 字头永远朝外
    for _grp, _c, _mid, _rr, _ri, _span in _ring_labels:
        _deg = np.degrees(_mid) % 360.0
        _rot = _deg - 90.0
        if _rot > 90.0:
            _rot -= 180.0
        elif _rot < -90.0:
            _rot += 180.0
        # 弧度按半径缩放：内环半径小 → 同字号占比大，故内环用更小字号
        _fs = 4.6 if _ri == 1 else 4.2
        ax.text(_mid, _rr, _grp, rotation=_rot, rotation_mode='anchor',
                ha='center', va='center', fontsize=_fs, color='#222222',
                zorder=7,
                bbox=dict(boxstyle='round,pad=0.12', fc='white',
                          ec='none', alpha=0.65))

    _outer_edge = _r_base + 2 * _RING_W + _RING_GAP
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines['polar'].set_visible(False)
    ax.set_ylim(0, max(_outer_edge * 1.02, max_x * 1.08))
    _ring_note = ''
    if ring_in and ring_out:
        _ring_note = f'; rings: inner = {ring_in_lab}, outer = {ring_out_lab}'
    elif ring_in:
        _ring_note = f'; ring = {ring_in_lab}'
    ax.set_title(f'Maximum-likelihood tree  (n = {n_tips} tips'
                 + (f'; branch colour = {_label})' if attr_of else ')'),
                 fontsize=7, fontweight='bold', pad=2)
    if attr_of:
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], color=grp_color[g], lw=1.5, label=g)
                   for g in groups_sorted]
        ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, 1.0),
                  fontsize=5.2, frameon=True, framealpha=0.9,
                  handlelength=1.3, handletextpad=0.5, labelspacing=0.28,
                  borderpad=0.4)
    # ⚠⚠ 2026-10-13 修：底部注记原先用**单行 `fig.text(0.01, 0.005, …, ha='left')`**，
    #   而注记串会随「宿主分组 + 双环」增长 → 实测 5.836 in > 可用 4.534 in（超出 1.29×）
    #   → **右端被画布裁掉且零报错**（硬规则：「被面板裁掉的文字零报错，必须量像素」）。
    #   修法 = **量宽折行**（复用 `pipeline_overview_fig._text_w_in/_wrap_caption`），
    #   行数决定 y 起点（多行时整体上移，保持底线对齐）。
    try:
        from utils.pipeline_overview_fig import _text_w_in, _wrap_caption
    except Exception:                                     # 独立脚本模式（PYTHONPATH=.）
        _text_w_in = _wrap_caption = None
    _note_full = (_radial_note
                  + (f'; {len(attr_of)}/{n_tips} tips assigned to '
                     f'{len(groups_sorted)} {_short}' if attr_of else '')
                  + _ring_note)
    _fs_note = 5.5
    _fig_w_in = fig.get_figwidth()
    _avail_in = _fig_w_in - 0.02                          # 左右各留 1% 边距
    if _wrap_caption is not None:
        # ⚠ 先按**分号边界**切成语义块，再贪心地把块装行 —— 纯空格折行会在
        #   `rings: inner =` 后断开、把 "host plant" 拆散（难看且易误读）。
        _segs = [s.strip() for s in _note_full.split('; ') if s.strip()]
        _note_lines = []
        for _i, _seg in enumerate(_segs):
            _cand = (_note_lines[-1] + '; ' + _seg) if _note_lines else _seg
            # ⚠⚠ 首块时 `_note_lines` 为空 → `_cand == _seg`，此时**必须 append**
            #   而非 `[-1]=`（越界 IndexError；实测服务器上直接崩）。
            if not _note_lines:
                _note_lines.append(_cand)
            elif _text_w_in(_cand, _fs_note) <= _avail_in:
                _note_lines[-1] = _cand                    # 能接就接在上一行尾
            elif _text_w_in(_seg, _fs_note) <= _avail_in:
                _note_lines.append(_seg)                   # 整块另起一行
            else:
                _note_lines.extend(_wrap_caption(_seg, _avail_in,
                                                 fs=_fs_note))   # 单块仍超宽 → 折
        if not _note_lines:
            _note_lines = [_note_full]
    else:                                                 # 兜底：按字符比例粗分
        _est = _text_w_in(_note_full, _fs_note) if _text_w_in else \
            (0.62 * _fs_note / 72.0 * len(_note_full))
        _k = max(1, int(_est / _avail_in) + 1)
        _note_lines = [_note_full[i::_k] for i in range(_k)]
    # ⚠ 字体高度按「点 → 英寸」换算（1 pt = 1/72 in），行高取 1.35× 字号
    _lh_in = _fs_note * 1.35 / 72.0
    _n_lines = len(_note_lines)
    fig.text(0.01, 0.005 if _n_lines == 1 else 0.005,
             '\n'.join(_note_lines),
             fontsize=_fs_note, color='#666666', ha='left',
             va='bottom', linespacing=1.35)
    fig.tight_layout(pad=0.3)
    # ⚠⚠ `tight_layout()` **不搬动 `fig.text`** → 必须在它之后**实测反推**是否溢出，
    #   仍溢出则按需要的宽度降字号（零报错缺陷的兜底闸门）。
    try:
        fig.canvas.draw()
        _ren = fig.canvas.get_renderer()
        _bb = [t.get_window_extent(renderer=_ren) for t in fig.texts]
        if _bb:
            _x1 = max(b.x1 for b in _bb) / fig.dpi          # 英寸
            if _x1 > _fig_w_in - 0.01:
                _scale = (_fig_w_in - 0.01) / _x1
                for t in fig.texts:
                    t.set_fontsize(t.get_fontsize() * _scale)
    except Exception:
        pass
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
    ap.add_argument('--style', default='colored', choices=['colored', 'classic'],
                    help='tree 专用: colored=Fig3 半栏彩版 / classic=单色大图+全标签')
    ap.add_argument('--meta', default=None,
                    help='tree 专用: 显式宿主/地点表 CSV (缺省自动发现 data/)')
    ap.add_argument('--branch-color', default='none',
                    choices=['none', 'host', 'location'],
                    help='tree 专用 (classic): 枝色映射 (众数继承, 混合枝灰)')
    a = ap.parse_args()
    if a.kind == 'align_qc':
        print(plot_align_qc(a.input, a.outdir))
    elif a.kind == 'tree':
        print(plot_tree_circular(a.input, a.outdir, style=a.style,
                                 meta=a.meta if a.meta else None,
                                 branch_color=a.branch_color))
    elif a.kind == 'popgen':
        print(plot_popgen(a.input, a.outdir))
    elif a.kind == 'gene_dating':
        t, p = plot_gene_dating_summary(a.input, a.outdir)
        print(t, p)
