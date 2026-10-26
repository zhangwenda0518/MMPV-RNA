#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_novelty_sankey.py — 09b「已知/新 -> 基因组型 -> 科 -> 属 -> 种」**单张**五层桑基图。

口径（用户 2026-10-13 拍板 + 七条追加；数据层 = make_novelty_sankey_data.py）
  · 主体 **61**（= 40 Known + 21 Novel；cap=1/属 去冗余后）
    ⚠ **类病毒按物种去重**（追加 1）：「类病毒太多了，其实种类没那么多，去冗余下」
      604 条 BLAST 命中只落在 5 个 species 上 → 每个 species 只算 1 个单位。
    ⚠ **剔除 Unassigned**（追加 2）：「去除unass结果」→ 删 `family == "(Unassigned)"` 的
      19 条病毒记录（类病毒强制豁免：`Rubber viroid India/2009` 的科本身就是 `(Unassigned)`）。
    ⚠ **剔除 L1 的「其他」**（追加 3）：「去除 L1 的 other」→ 删 tier == 其他 的 87 条，
      L1 只剩 Known / Novel。
    ⚠ **Novel 全部改名**（追加 4）：「为什么 L1 Novel 的结果，物种层面还是注释给了已知病毒的名字
      → 同样的属 goji 1 2 3 这样」→ Novel 的物种层一律 `<属>-goji[-N]`，
      **不再出现任何已知物种名**；Known 分支保持源表真名。
    ⚠ **长尾剔除**（追加 5）：「所有的 other 去掉」→ 图上**不再出现任何 `Other ...` 节点**。
      ⚠⚠ 2026-10-13 晚**范围收窄**：原「family/genus/species 三层逐层迭代」→ 改为
      **只在 species 层剔除**（科/属层保留单例，保证科–属–种链完整）。
      类病毒 + 新种候选豁免。
      ⚠ 消费端阈值必须随之置 0，否则单例又会被聚成 `Other species`（见 `THR_EFF`）。
    ⚠⚠ **同属 Novel 去冗余 cap=1/属**（追加 6，**本轮最关键**）：
      「同属的新病毒只展示一个，这样避免新病毒过多引起桑基图不好看，
       **详细的结果放附图**」→ 数据层两步：① 同属内同一物种名(goji-N)只留 1 条
       （取发生率最高）；② 仍 > 1 则按发生率砍到 1 条。**Known 完全不动**。
       → Novel 81 → 21（64 条折叠）；完整明细另存 `novel_dedup_full.tsv`，
       由附图 `Fig_GojiVirome_NoveltySankey_Supp` 承载（信息零丢失）。
      ⚠ 去冗余后**不再做二次去长尾**（否则会删掉 26 个「单物种新病毒属」，
      与「展示较完整」矛盾 —— 实测 61 vs 58）。
  · L1 = Known(40) / Novel(21)
  · L2 = DNA(5) / RNA(51) / Viroid(5)
  · L3/L4/L5 = family / genus / species（17 / 25 / 37 个节点，无聚合块）
  · 11 条「新种候选」的物种标签**加粗**（cap=1 后由 21 → 11）
  · L5 物种标签按 tier 着色：绿 = Known / 红 = Novel（追加 7，用户：「新病毒标红色字，已知标绿色字」）
  · L5 标签附发生率 `(xx.x%)`（独立图；嵌入主图面板 D 时用 plain 版无后缀）

⚠ 为什么是**单张**（用户追加 3：「ab 不能合一起嘛」）
  旧版拆双面板是因为类病毒 604/895 = 67.5% 把科/属/种三层挤成各 1 个节点，
  291 条病毒要在 32.5% 高度里摊给 152 个种 → n=1 节点只有 ≈2 px，标签物理不可读。
  去冗余 + 去 Unassigned + 去 Other + 去长尾 + cap=1/属后总量降到 61（k_max = 37 节点），
  1 单位节点仍有充足像素（@300 dpi），标签靠去碰撞 + 引导线排布 → 单张图即可承载。

⚠ 版式硬规则（本脚本踩过）：
  · 列内排序必须按**上一层的分组序**（否则彩带大面积交叉）
  · y 轴从 0 向上 → 首节点落在底部；要「大者在上」须做 `Y - y` 翻转
  · **列内 gap 全局统一 + 各列在 [0,Y] 内居中**（`solve_geometry`）：L1 只有 2 个类别，
    若按「每列铺满同一 Y」会把 gap 撑到 800+，图上出现大片空白
  · 小节点标签去碰撞 + 灰色引导线（最小行距按**面板实际英寸高**换算）
  · 注记必须**拆成多行短句**，否则会撞右侧图例
  · 末列标签宽度 = `xlim 右端 - (末列 x + NW)`，必须按字号实测（见 `_text_w_in`）
  · 文本一律 ASCII（项目无 CJK 字体配置；R 侧 mbcsToSbcs 亦有同类静默坑）
  · 图内数字一律取自 `manifest.json`，禁止硬编码（口径切换会静默说谎）

复现：
  "C:/Python312/python.exe" reports/virome_analysis_pipeline/make_novelty_sankey_data.py
  "C:/Python312/python.exe" reports/virome_analysis_pipeline/fig_novelty_sankey.py
  "C:/Python312/python.exe" reports/virome_analysis_pipeline/_check_novelty_sankey_overlap.py

并入主图（用户 2026-10-13：「和主图放一起可以吗」→ 选「主图改通栏」）：
  # ① 无发生率后缀的嵌入资产（主图 D 用；通栏面板里后缀会把标签顶出画布）
  "...fig_novelty_sankey.py" --embed-asset-plain <out>/sankey_novelty_embed.png \
                             --embed-size 14.00x3.28
  # ② 主图：D 面板升为通栏 + 底行重排（见 make_discovery_sci_figs.py::LAYOUT）
  #    资产需落在 <outdir>/（`_find_asset` 首位候选）
"""
from __future__ import annotations

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

import argparse
import csv
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "reports" / "virome_analysis_pipeline" / "novelty_sankey_data"
DEFAULT_OUT = (ROOT / "reports" / "_discovery_sci_20261012" / "_lt" / "10_Reports" / "S10_SCI_Discovery")

# ── 调色板：按显式文件路径加载（项目纪律，勿改成 import）──
_pal = ROOT / "virome_discovery_pipeline" / "utils" / "sci_palette.py"
_spec = importlib.util.spec_from_file_location("sci_palette", _pal)
sp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sp)

C_KNOWN = sp.SEMANTIC["genome"]["DNA"]      # #0072B2
C_NOVEL = sp.SEMANTIC["genome"]["RNA"]      # #D55E00
C_OTHER = "#999999"
TIER_COLOR = {"已知": C_KNOWN, "新": C_NOVEL, "其他": C_OTHER}
TIER_ORDER = ["已知", "新", "其他"]
TIER_LABEL = {"已知": "Known", "新": "Novel", "其他": "Other"}
TYPE_COLOR = dict(sp.SEMANTIC["genome"])
TYPE_ORDER = ["DNA", "RNA", "Viroid", "Unassigned"]

NW = 0.075
FS_NODE = 4.05
FS_TITLE = 6.8
FS_NOTE = 4.4
THR = {"family": 2, "genus": 2, "species": 2}
# 长尾聚合阈值（**仅 species 层**）——仅在**数据层未剔除长尾**（`manifest["drop_longtail"] == False`）时生效。
# 当前默认口径（去 Other + 去 Unassigned + Novel 全 goji 化 + **只 species 层去长尾**
#   + **同属 Novel 去冗余 cap=1**）总 **61** 单位，
# 数据层已把物种层单位数 < 2 的类群删掉 → 阈值被 `THR_EFF` 置 0（不再聚合），
# 列节点数 = 2 / 3 / 17 / 25 / 37，图上无任何 `Other ...` 块。
# 关掉剔除（`--no-drop-longtail`）时回到旧口径：总 190，阈值 2 生效，
#   family 29 唯一 / genus 55 / species 117（单例 91）→ 保留 22 / 31 / 26 个 + 3 个长尾块。


def _text_w_in(s: str, fontsize: float) -> float:
    try:
        from matplotlib.textpath import TextPath
        fp = FontProperties(size=fontsize)
        return TextPath((0, 0), s, prop=fp).get_extents().width / 72.0
    except Exception:
        return len(s) * fontsize * 0.5 / 72.0


def read_tsv(p: Path):
    with p.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ══════════════════════════════════════════════════════════════════════
def build_columns(contigs, thr):
    """返回 cols（每列 [(显示名, n)]，首元素在最上方）、各层 name->display 映射。"""
    cnt = {lv: Counter(r[lv] for r in contigs) for lv in ("tier", "genome_type", "family", "genus", "species")}

    # 每个 family/genus/species 的「主」上层归属（按该上层取值计数最大者）
    def primary(lv, up):
        out = {}
        for r in contigs:
            out.setdefault(r[lv], Counter())[r[up]] += 1
        return out
    fam_gt = primary("family", "genome_type")
    gen_fam = primary("genus", "family")
    sp_gen = primary("species", "genus")

    def split(c, key, th, protect):
        keep = sorted([k for k, v in c.items() if v >= th or k in protect], key=key)
        rest_n = sum(v for k, v in c.items() if v < th and k not in protect)
        rest_k = sum(1 for k, v in c.items() if v < th and k not in protect)
        return keep, rest_n, rest_k

    # ⚠ 两类节点**强制保留**（否则会被长尾聚合整条吞掉）：
    #   ① 类病毒谱系：按物种去重后每条只有 1 个单位，低于任何阈值；
    #   ② 21 条新种候选：合并面板后它们就是「主角」，必须在种层可读。
    prot = {lv: {r[lv] for r in contigs
                 if r["genome_type"] == "Viroid" or r["is_new_species"] == "True"}
            for lv in ("family", "genus", "species")}

    f_keep, f_rn, f_rk = split(cnt["family"],
                               lambda k: (min(TYPE_ORDER.index(g) for g in fam_gt[k]), -cnt["family"][k], k),
                               thr["family"], prot["family"])
    f_disp = {k: k for k in f_keep}
    if f_rn:
        A = f"Other families (n={f_rn:,}; {f_rk} taxa)"
        for k, v in cnt["family"].items():
            if v < thr["family"] and k not in prot["family"]:
                f_disp[k] = A
    fam_idx = {k: i for i, k in enumerate(f_keep)}

    ge_keep, ge_rn, ge_rk = split(
        cnt["genus"],
        lambda k: (fam_idx.get(max(gen_fam[k], key=lambda f: gen_fam[k][f]), 99), -cnt["genus"][k], k),
        thr["genus"], prot["genus"])
    ge_disp = {k: k for k in ge_keep}
    if ge_rn:
        A = f"Other genera (n={ge_rn:,}; {ge_rk} taxa)"
        for k, v in cnt["genus"].items():
            if v < thr["genus"] and k not in prot["genus"]:
                ge_disp[k] = A
    gen_idx = {k: i for i, k in enumerate(ge_keep)}

    sp_keep, sp_rn, sp_rk = split(
        cnt["species"],
        lambda k: (gen_idx.get(max(sp_gen[k], key=lambda g: sp_gen[k][g]), 99), -cnt["species"][k], k),
        thr["species"], prot["species"])
    sp_disp = {k: k for k in sp_keep}
    if sp_rn:
        A = f"Other species (n={sp_rn:,}; {sp_rk} taxa)"
        for k, v in cnt["species"].items():
            if v < thr["species"] and k not in prot["species"]:
                sp_disp[k] = A

    def col_of(keys, agg_label):
        return list(keys) + ([agg_label] if agg_label else [])

    def agg_name(disp):
        return next((v for v in disp.values() if v.startswith("Other")), None)

    f_all = col_of(f_keep, agg_name(f_disp))
    ge_all = col_of(ge_keep, agg_name(ge_disp))
    sp_all = col_of(sp_keep, agg_name(sp_disp))

    cols = [
        [(t, cnt["tier"][t]) for t in TIER_ORDER if cnt["tier"][t]],
        [(g, cnt["genome_type"][g]) for g in TYPE_ORDER if cnt["genome_type"][g]],
        [(n, sum(v for k, v in cnt["family"].items() if f_disp[k] == n)) for n in f_all],
        [(n, sum(v for k, v in cnt["genus"].items() if ge_disp[k] == n)) for n in ge_all],
        [(n, sum(v for k, v in cnt["species"].items() if sp_disp[k] == n)) for n in sp_all],
    ]
    srcs = [lambda r: r["tier"], lambda r: r["genome_type"],
            lambda r: f_disp[r["family"]], lambda r: ge_disp[r["genus"]],
            lambda r: sp_disp[r["species"]]]
    E = []
    for L in range(4):
        c = Counter()
        for r in contigs:
            c[(srcs[L](r), srcs[L + 1](r), r["tier"])] += 1
        E.append([(s, t, n, tier) for (s, t, tier), n in c.items()])
    return cols, E, {0: TIER_LABEL}, {2: f_disp, 3: ge_disp, 4: sp_disp}


# ══════════════════════════════════════════════════════════════════════
def solve_geometry(cols, ax_h_in, label_fs):
    """解出坐标高 Y 与**列内 gap**。

    ⚠ 为什么不用「每列铺满同一 Y」：L1 只有 2 个类别（Known / Novel），
      铺满会把 (Y-tot)/(k-1) 撑到 800+，图上出现大片空白；而 species 列 50+ 个节点
      又需要足够行距。→ 所有列共用一个 gap（= 标签可读行距），各列在 [0, Y] 内居中。
      节点高度仍严格 ∝ 流量，只是列高随该列节点数变化。
    方程：gap = c·Y（c = 行距换算），Y = tot + gap·(k_max-1)  ⇒  Y = tot / (1 - c·(k_max-1))

    ⭐⭐ 2026-10-13（用户「图D的前面 L1-L4 变大些」）**实测结论 —— 只能靠加高画布**：
      · 「L1-L4 显小」的根因 = **Y 被节点最多的列（species，37 个）绑架**。
        gap 是**物理量**（`fs×1.30/72` 英寸），在固定画布高下 L5 的 37 个标签把 gap
        顶到地板值（4.05 pt @300 dpi → **21.94 px**，实测现行 22.97 px 只余 5 %），
        ⇒ **gap 不能再压** ⇒ Y 不能小 ⇒ L1 的 Known 只占 91 px / 927 px 画布。
      · ⚠⚠ **试过且失败的两种「局部放大」**（会静默破坏版式）：
        ① 给 L1-L4 的 gap 单独乘 boost → L3(17 节点)/L4(25 节点) 的列高暴涨到
           **157 % / 228 % of Y**，直接冲出画布；② 压小 gap → L5 的 37 个标签
           行距掉到 9–19 px < 21.94 px 地板 → 互相压住。
      · ✅ **唯一正解 = 加高嵌入画布**：`1 单位像素 = 画布高 / Y`，而 Y 只由 gap 地板
        决定（与画布高**无关**）→ 画布越高，单位像素线性增长：
          3.10 in → 2.3 px/单位（Known 91 px） ｜ 4.40 in → 8.1 px/单位（Known 323 px）
        行距恒 ≈22.97 px（物理量）→ **加高不会造成标签重叠**。
      ⇒ 由调用方（`--embed-size` 的 H）决定放大量级；本函数保持纯几何、不再引入 boost。
    """
    tot = sum(n for _, n in cols[0])
    k_max = max(len(c) for c in cols)
    if k_max <= 1:
        return tot * 1.05, 0.0
    c = label_fs * 1.30 / 72.0 / ax_h_in
    Y = tot / max(0.15, 1.0 - c * (k_max - 1))      # 下限 0.15 防止节点过多时除零/爆炸
    return Y, c * Y


def layout(cols, Y, gap):
    """每列在 [0, Y] 内居中；**首元素在最上方**（y 向上）。"""
    pos = {}
    for L, names in enumerate(cols):
        k = len(names)
        tot = sum(n for _, n in names)
        g = gap if k > 1 else 0.0
        h = tot + g * (k - 1)
        y = (Y - h) / 2.0
        for nm, n in names:
            pos[(L, nm)] = (Y - y - n, Y - y)
            y += n + g
    return pos


def ribbon(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha, zorder):
    xm = (x0 + x1) / 2.0
    verts = [(x0, y0a), (xm, y0a), (xm, y1a), (x1, y1a),
             (x1, y1b), (xm, y1b), (xm, y0b), (x0, y0b), (x0, y0a)]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
             MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=color, edgecolor="none",
                           alpha=alpha, zorder=zorder))


def _nw(L: int, node_fs_scale) -> float:
    """该层节点块宽（默认 `NW`；`node_fs_scale[L]` 放大）。

    ⚠ 只影响**块宽**（视觉醒目度），**不影响块高** —— 块高恒 ∝ 单位数。
    """
    return NW * float((node_fs_scale or {}).get(L, 1.0))


def draw(ax, cols, edges, node_color, ribbon_color, label_side, Y, xmax,
         label_fs=FS_NODE, label_map=None, bold_labels=None, gap=None, color_of=None,
         xmin=-0.80, level_fs=None, node_fs_scale=None):
    """画桑基本体。

    ⭐ 2026-10-13（用户「图D的前面 L1-L4 变大些」）新增两个**逐层**参数：
      · `level_fs`      = `{L: fontsize}`，逐层覆盖标签字号（未给则用 `label_fs`）；
      · `node_fs_scale` = `{L: k}`，逐层放大**节点块**宽度与最小行距（k>1 让前四层更醒目）。
        ⚠⚠ **绝不改变块高 ∝ 单位数** —— 只放大 `NW`（块宽）与 label 行距。
    """
    label_map = label_map or {}
    bold_labels = bold_labels or set()
    color_of = color_of or (lambda L, base: "#111111")
    level_fs = level_fs or {}
    node_fs_scale = node_fs_scale or {}
    ncol = len(cols)
    xs = [float(i) for i in range(ncol)]
    if gap is None:                                   # 兜底：退回「每列铺满同一 Y」
        tot = sum(n for _, n in cols[0])
        gap = (Y - tot) / max(1, len(cols[0]) - 1)
    pos = layout(cols, Y, gap)
    idx = [{nm: i for i, (nm, _) in enumerate(col)} for col in cols]

    for L, E in enumerate(edges):
        cur = {nm: pos[(L, nm)][1] for nm, _ in cols[L]}          # 自上而下
        src = {}
        for e in sorted(E, key=lambda e: (idx[L][e[0]], idx[L + 1][e[1]], TIER_ORDER.index(e[3]))):
            y1 = cur[e[0]]
            src[e] = (y1 - e[2], y1)
            cur[e[0]] = y1 - e[2]
        cur = {nm: pos[(L + 1, nm)][1] for nm, _ in cols[L + 1]}
        for e in sorted(E, key=lambda e: (idx[L + 1][e[1]], idx[L][e[0]], TIER_ORDER.index(e[3]))):
            y1 = cur[e[1]]
            ribbon(ax, xs[L] + _nw(L, node_fs_scale), xs[L + 1], src[e][0], src[e][1],
                   y1 - e[2], y1, ribbon_color(e), 0.46 if L == 0 else 0.32, 1 + L)
            cur[e[1]] = y1 - e[2]

    # ── 节点 + 标签（含去碰撞 + 引导线）──
    labels = []
    for L, names in enumerate(cols):
        for nm, n in names:
            y0, y1 = pos[(L, nm)]
            ax.add_patch(Rectangle((xs[L], y0), _nw(L, node_fs_scale), y1 - y0,
                                   facecolor=node_color(L, nm),
                                   edgecolor="white", linewidth=0.35, zorder=5))
            if label_side.get(L) is not None:
                labels.append([L, (y0 + y1) / 2.0, label_map.get(L, {}).get(nm, nm),
                               (y0 + y1) / 2.0, nm])       # 末位 = **基础名**（加粗判定用，不含发生率后缀）

    h_in = ax.get_position().height * ax.figure.get_figheight()

    def _fs(L):
        return float(level_fs.get(L, label_fs))

    for L in sorted({l[0] for l in labels}):
        # ⚠ 逐层最小行距（`level_fs` 放大过的层要用放大后的字号算，否则大字号标签互相压住）
        mg = (_fs(L) * 1.30 / 72.0) / h_in * Y
        grp = sorted([l for l in labels if l[0] == L], key=lambda l: -l[1])
        for i in range(1, len(grp)):
            if grp[i - 1][1] - grp[i][1] < mg:
                grp[i][1] = grp[i - 1][1] - mg
        if grp and grp[-1][1] < 0.0:                     # 越界则整组上移
            d = -grp[-1][1]
            for l in grp:
                l[1] += d
    for L, ly, txt, ny, base in labels:
        side = label_side[L]
        _w = _nw(L, node_fs_scale)
        off = _w + 0.05 if side == "right" else -0.05 - (_w - NW)
        xl = xs[L] + off
        _mg = (_fs(L) * 1.30 / 72.0) / h_in * Y
        if abs(ly - ny) > _mg * 0.45:                    # 位移过大 → 补引导线
            ax.plot([xs[L] + (_w if side == "right" else 0.0), xl - 0.012],
                    [ny, ly], color="#8A8A8A", lw=0.35, zorder=6, solid_capstyle="butt")
        ax.text(xl, ly, txt, ha="left" if side == "right" else "right", va="center",
                fontsize=_fs(L), color=color_of(L, base), zorder=7,
                fontweight="bold" if base in bold_labels else "normal",
                path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-Y * 0.015, Y * 1.015)
    ax.axis("off")



# ══════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(DEFAULT_OUT))
    ap.add_argument("--stem", default="Fig_GojiVirome_NoveltySankey")
    # ── 嵌入模式（用户 2026-10-13：「和主图放一起可以吗」→ 选「替换 D 面板」）──
    #   为 `Fig_GojiVirome_Discovery` 的 D 面板产出一张**不带注记**的贴图资产。
    #   ⚠ 为什么不把桑基逻辑搬进主图脚本：本项目已有「同一面板两份标题副本」
    #     （直渲 + 旁路）的静默坑 —— 桑基两处实现必然漂移。这里让**唯一实现**
    #     多输出一个资产，主图只负责贴图（与 D 原有的 `sankey_flow_plant.png` 旁路同构）。
    ap.add_argument("--embed-asset", default=None,
                    help="额外输出一张供主图 D 面板贴图的 PNG（无标题/无注记，紧密裁边）")
    ap.add_argument("--embed-asset-plain", default=None,
                    help="同 `--embed-asset`，但**不带发生率后缀** —— 供主图这类「标签必须极短」"
                         "的场景（挤在通栏面板里时后缀会把标签顶出画布）。")
    ap.add_argument("--embed-no-pct", action="store_true",
                    help="单独使用时：`--embed-asset` 也不加发生率后缀。")
    ap.add_argument("--embed-size", default="12.60x4.40",
                    help="嵌入资产的画布英寸尺寸 `WxH`。默认 12.60x4.40 是为把 **tight 裁边后**的"
                         "宽高比钉到 **2.88**（= 主图 D 面板 4.70x1.63 的比 2.883）——"
                         "⚠ 裁边会把两侧留白吃掉，故画布要先给足宽度；直接用面板比会偏窄到 2.2 而留白。")
    # ── ⭐ 2026-10-13（用户「图D的前面 L1-L4 变大些」）──
    #   ⚠⚠ **放大量级由 `--embed-size` 的 H 决定**（见 `solve_geometry` 实测结论）：
    #     L1-L4 显小的根因是 Y 被 species 列（37 节点）绑架，而 gap 已在地板值
    #     （21.94 px @ fs 4.05 pt，300 dpi）→ **只能加高画布**把 1 单位像素放大。
    #   本组开关只做「同一画布高内」的**额外**强化（字号 + 块宽），不碰几何。
    ap.add_argument("--front-fs", type=float, default=0.0,
                    help="L1–L4 标签字号（0 = 用 `FS_NODE`）。⚠ 放大后每列行距按新字号"
                         "自动重算（`draw` 内逐层），但**末列 species 字号不变**。")
    ap.add_argument("--front-node-scale", type=float, default=1.0,
                    help="L1–L4 **节点块宽度**放大系数（视觉醒目度）。⚠ 只改块宽 `NW`，"
                         "**块高恒 ∝ 单位数**（图注的严格比例声明不受影响）。")
    # ── ⭐ 2026-10-13 21:32（用户「D图都看不清字体大些」）──
    #   ⚠⚠ **根因不是「字号设小了」，而是「贴图被缩小」**：资产画布宽 10.43 in 贴进
    #     D 面板 7.61 in → **水平缩放 0.729**，把资产内 5.4 pt 压成面板内等效 **3.94 pt**、
    #     4.05 pt（L5）压成 **2.95 pt** —— 远低于印刷可读下限（≈6 pt）。
    #   ⇒ 两手一起做：① `--all-fs` 抬高**全部五层**字号（含 L5，`--front-fs` 只管 L1–L4）；
    #     ② `--embed-size` 的 **W 加宽**，让 tight 后的资产宽 ≈ 面板宽（缩放比→1.0）。
    ap.add_argument("--all-fs", type=float, default=0.0,
                    help="**全部五层**（含末列 species）标签字号（0 = 用 `FS_NODE`）。"
                         "⚠ 与 `--front-fs` 的区别：`--front-fs` 只覆盖 L1–L4，本参数覆盖 L5。"
                         "⚠ 抬高 L5 字号会**同步抬高 gap 地板**（`solve_geometry` 的 c 项），"
                         "故画布高可能要一起加 —— 见 `--embed-size`。")
    args = ap.parse_args()

    contigs = read_tsv(DATA / "contigs.tsv")
    man = json.loads((DATA / "manifest.json").read_text(encoding="utf-8"))
    TOTAL = int(man["total"])                                   # 296（类病毒按物种去重）
    # 自检三判据：① 总数 ② tier 分层 ③ 类病毒单位数
    assert len(contigs) == TOTAL, f"contigs 应为 {TOTAL} 行，实得 {len(contigs)}"
    tc = Counter(r["tier"] for r in contigs)
    assert dict(tc) == man["tiers"], f"tier 分布不符：{dict(tc)} vs {man['tiers']}"
    gtc = Counter(r["genome_type"] for r in contigs)
    assert dict(gtc) == man["genome_types"], f"genome_type 分布不符：{dict(gtc)} vs {man['genome_types']}"
    n_vir = gtc.get("Viroid", 0)
    assert n_vir == man["viroid_units"], \
        f"类病毒单位数不符：{n_vir} vs manifest {man['viroid_units']}（口径 {man['viroid_unit']}）"
    nov = [r for r in contigs if r["is_new_species"] == "True"]
    # ⚠ 新种候选数**不得硬编码**：2026-10-13 用户追加「同属 Novel 去冗余（每属 cap 3）」后，
    #   属内候选数 >3 的会被砍（Betacytorhabdovirus 4→2、Carlavirus 3→2）→ 21 → 18。
    #   → 一律从 manifest 的 `new_species`（数据层实际保留数）读，防静默说谎。
    n_ns_expect = int(man.get("new_species", len(nov)))
    assert len(nov) == n_ns_expect, \
        f"新种候选应为 {n_ns_expect} 条（manifest），实得 {len(nov)}"
    print(f"口径 = {man['viroid_unit']}  总 {TOTAL}  分层 {dict(tc)}  "
          f"类型 {dict(gtc)}  类病毒单位 {n_vir}")

    plt.rcParams.update(sp.sci_rcparams(base_size=8))
    fig = plt.figure(figsize=(9.32, 8.40))
    axA = fig.add_axes([0.030, 0.118, 0.958, 0.775])     # 单面板（用户 2026-10-13：「ab 不能合一起嘛」）

    # ── 单张五层桑基 ──
    # ⚠⚠ 用户 2026-10-13 追加「所有的 other 去掉」→ 数据层已**逐层剔除**长尾类群
    #    （`manifest["drop_longtail"] == True`）。此时消费端**绝不能再聚合**：
    #    阈值置 0（`v >= 0` 恒真）→ 全部类群各自成节点，`Other ...` 块不再产生。
    #    否则 26 个被豁免的单例（21 新种候选 + 5 类病毒）又会被聚成 `Other species`，等于白做。
    THR_EFF = {k: 0 for k in THR} if man.get("drop_longtail") else dict(THR)
    cols, E, lmap, disp = build_columns(contigs, THR_EFF)
    # ── 自检：图上不得出现任何 `Other ...` 长尾块 ──
    _agg = [n for col in cols for n, _ in col if n.startswith("Other ")]
    if man.get("drop_longtail") and _agg:
        sys.exit(f"❌ 数据层已剔除长尾，但图上仍出现聚合块：{_agg}")
    if not man.get("drop_longtail") and not _agg:
        print("   ⚠ 未开长尾剔除，但也没有聚合块（阈值可能已被改动）")
    # ── 自检：开了剔除后，**仅 species 层**不得残留未被豁免的 1 单位类群 ──
    #   ⚠⚠ 2026-10-13 晚（两次口径变更后）：
    #     ① 同属 Novel 去冗余（cap 1/属）会把很多属压到 1 个物种
    #        （`Soymovirus`/`Closterovirus`/`Potexvirus`…）→ 用户定「**去冗余后不再剔物种层**」
    #        （否则 26 个真实新病毒属整个被删，与「展示较完整」矛盾）。
    #     ② 科/属层的单例是合法保留（保证科–属–种链完整）。
    #   → 本自检只在**去冗余关闭**时才对物种层严格检查；
    #     去冗余开启时，物种层单例是**预期行为**（代表按发生率选出的每属 1 条），跳过。
    _dedup_on = bool(man.get("novel_dedup_per_genus"))
    if man.get("drop_longtail") and not _dedup_on:
        _ltscope = man.get("longtail_scope", ["species"])
        _prot = {lv: {r[lv] for r in contigs
                      if r["genome_type"] == "Viroid" or r["is_new_species"] == "True"}
                 for lv in _ltscope}
        for lv in _ltscope:
            _c = Counter(r[lv] for r in contigs)
            _bad = {k: v for k, v in _c.items() if v < 2 and k not in _prot[lv]}
            if _bad:
                sys.exit(f"❌ {lv} 层仍有未豁免的 1 单位类群（应已在数据层剔除）：{list(_bad)[:5]}")
    elif _dedup_on:
        print("   ℹ 已开同属 Novel 去冗余 → 物种层单例为预期（每属只展示代表），跳过长尾自检")
    node_color = lambda L, nm: (TIER_COLOR.get(nm, "#333333") if L == 0
                                else TYPE_COLOR.get(nm, "#5A5A5A") if L == 1
                                else "#C4C4C4")
    # 21 条新种候选的物种节点标签**加粗**（合并面板后取代 Panel b 的「逐条可读」功能）
    bold = {r["species"] for r in contigs if r["is_new_species"] == "True"}

    # ── L5 物种层标签按 tier 着色（用户 2026-10-13：「新病毒标红色字，已知标绿色字」）──
    #   ⚠⚠ 只有 **L5** 能唯一着色：L2 / L3 / L4 大量 tier 混合（如 `Betaflexiviridae`
    #      已知 20 / 新 16、`Carlavirus` 已知 18 / 新 7）→ 按主导 tier 染色会**误导**，
    #      故只对 L == 4 生效，其余层保持黑色。
    #   ⚠ 旧口径下的长尾块 `Other species (…)` 混合（已知 3 / 新 62）→ 无唯一 tier → 黑字。
    #     2026-10-13 长尾块已整体去除 → **51 个 L5 节点全部 tier 唯一**，全部着色（无黑字残留）。
    #   色值取项目 Okabe-Ito：bluish green `#009E73`（已知）/ vermillion `#D55E00`（新）。
    TIER_TEXT = {"已知": "#009E73", "新": "#D55E00"}
    tier_of = {}
    for r in contigs:
        tier_of.setdefault(r["species"], set()).add(r["tier"])

    def label_color(L, base):
        if L != 4:
            return "#111111"
        ts = tier_of.get(base)
        return TIER_TEXT[next(iter(ts))] if ts and len(ts) == 1 else "#111111"

    # ── L5 发生率附注（用户 2026-10-13：「发生率可以补上嘛」→ 物种层 + 多条取 max）──
    #   ⚠ 发生率是「样本维度」（分母 manifest["prevalence_denom"]），节点高是「contig 维度」
    #     → 只作标签附注，**不参与任何几何计算**。
    DENOM = int(man["prevalence_denom"])
    prev_max = {}
    for r in contigs:
        if r["prevalence_pct"] == "":
            continue
        d = disp[4][r["species"]]                        # 聚合块的成员也归到块名上
        prev_max[d] = max(prev_max.get(d, -1.0), float(r["prevalence_pct"]))
    l5_label = {nm: (f"{nm} ({prev_max[nm]:.1f}%)" if nm in prev_max else nm)
                for nm, _ in cols[4]}
    lmap = dict(lmap); lmap[4] = l5_label

    Y, gap = solve_geometry(cols, axA.get_position().height * fig.get_figheight(),
                            max([FS_NODE] + ([args.front_fs] if args.front_fs > 0 else [])
                                + ([args.all_fs] if args.all_fs > 0 else [])))
    _lfs = {}
    _nsc = {}
    if args.front_fs > 0:
        for _L in range(min(4, len(cols) - 1)):
            _lfs[_L] = args.front_fs
    if args.front_node_scale != 1.0:
        for _L in range(min(4, len(cols) - 1)):
            _nsc[_L] = args.front_node_scale
    # ⚠⚠ 2026-10-13 v21.3：`--all-fs` **必须同样作用于独立图**。
    #   否则会踩一个静默坑：只传 `--all-fs`（不传 `--front-fs`）时，
    #   独立图的 L1–L4 会**回落到 `FS_NODE`（4.05）**，比上一轮的 5.4 还小 ——
    #   即「改主图字号」反而**让独立图变小**（实测 md5 从 `7cdd6c9c` 变 `68c16fb9`）。
    _base_fs = FS_NODE
    if args.all_fs > 0:
        _base_fs = args.all_fs
        for _L in range(len(cols)):
            _lfs[_L] = args.all_fs
        _nsc.update({_L: args.front_node_scale
                     for _L in range(min(4, len(cols) - 1))})
    draw(axA, cols, E, node_color, lambda e: TIER_COLOR[e[3]],
         {0: "left", 1: "right", 2: "right", 3: "right", 4: "right"},
         Y=Y, xmax=4.92, label_fs=_base_fs, label_map=lmap, bold_labels=bold, gap=gap,
         color_of=label_color, level_fs=_lfs, node_fs_scale=_nsc)

    # ── 标题 / 注记（全部数字取自 manifest.json）──
    tier_txt = "  |  ".join(f"{TIER_LABEL[t]} = {man['tiers'][t]}" for t in TIER_ORDER
                            if man["tiers"].get(t))
    n_vir_virus = man["tiers"]["已知"] - man["viroid_units"]
    tail_phrase = ("no long-tail aggregation (all three levels are complete)"
                   if man.get("drop_longtail") else "long tail aggregated per level")
    fig.text(0.030, 0.968, "Goji virome discovery - 09b novelty taxonomy",
             fontsize=9.0, fontweight="bold", ha="left", va="center", color="#111111")
    fig.text(0.030, 0.934, f"Taxonomic flow of the {TOTAL} high-confidence units",
             fontsize=FS_TITLE, fontweight="bold", ha="left", va="center", color="#111111")
    fig.text(0.030, 0.910,
             f"Tiers  {tier_txt}   (Known = {n_vir_virus} viruses + {man['viroid_units']} viroid species)"
             f"        {tail_phrase}; candidate new species are shown in bold",
             fontsize=FS_NOTE, ha="left", va="center", color="#555555")

    n_kept_unass = sum(1 for r in contigs if r["family"] == "(Unassigned)")
    n_kept_phrase = (f"{n_kept_unass} viroid record kept" if n_kept_unass == 1
                     else f"{n_kept_unass} viroid records kept")
    # ⚠ 注记**拆成多行短句**（每行 ≤ ~0.46 图宽），否则会撞右侧图例（x=0.585 起）
    for i, line in enumerate([
        f"Node height and ribbon width are strictly proportional to unit count (every column sums to {TOTAL}).",
        f"Viroids de-duplicated to species level (604 BLAST hits -> {man['viroid_units']} species);"
        f" {man['n_dropped']} Unassigned records dropped ({n_kept_phrase});"
        f" {man['n_dropped_other']} 'Other' tier records removed.",
        f"No long-tail aggregation: {man['n_dropped_longtail']} units in taxa with < {man['min_units']}"
        f" units were removed at the family/genus/species level; every remaining taxon is drawn separately.",
        "Novel branch is named <Genus>-goji[-N]; bold = the 21 candidate new species."
        "  No known species name appears on the Novel branch.",
        "Species labels are coloured by tier: green = Known, red = Novel.",
        f"Prevalence in brackets = n_samples / {DENOM} samples (max within each species node);"
        " label note only, not used in the layout.",
    ]):
        fig.text(0.030, 0.101 - i * 0.0170, line,
                 fontsize=FS_NOTE, color="#777777", ha="left", va="center")

    h_tier = [Line2D([], [], marker="s", ls="", ms=4.4, mfc=TIER_COLOR[t], mec="white",
                     label=TIER_LABEL[t]) for t in TIER_ORDER if man["tiers"].get(t)]
    # ⚠ 图例只列**数据里实际存在**的类别：删掉某类后若仍硬编码，图上会多出一个空类别（误导）。
    types_present = [t for t in TYPE_ORDER if any(r["genome_type"] == t for r in contigs)]
    h_type = [Line2D([], [], marker="o", ls="", ms=4.4, mfc=TYPE_COLOR[t], mec="white",
                     label=t) for t in types_present]
    lg1 = fig.legend(handles=h_tier, loc="lower left", bbox_to_anchor=(0.585, 0.020), ncol=3,
                     frameon=False, fontsize=5.0, handletextpad=0.35, columnspacing=0.9,
                     title="Tier", title_fontsize=5.0)
    lg1.get_title().set_fontweight("bold")
    lg1._legend_box.align = "left"
    lg2 = fig.legend(handles=h_type, loc="lower left", bbox_to_anchor=(0.725, 0.020),
                     ncol=len(types_present),
                     frameon=False, fontsize=5.0, handletextpad=0.35, columnspacing=0.9,
                     title="Genome type", title_fontsize=5.0)
    lg2.get_title().set_fontweight("bold")
    lg2._legend_box.align = "left"

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / f"{args.stem}.png"
    pdf = outdir / f"{args.stem}.pdf"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(pdf, facecolor="white")

    # ── 嵌入资产（供主图 D 面板贴图）：无标题 / 无注记 / 无图例，紧密裁边 ──
    #   ⚠ 只画桑基本体（`axA`），标题与 6 行注记属于**主图 D 自己的**职责。
    #   ⚠ 用 `bbox_inches="tight"` 让裁边后的宽高比由内容决定；主图侧用
    #     `imshow(aspect="equal")` 贴图，比不匹配会留白（D 现有旁路同款处理）。
    embed_png = None
    embed_png_plain = None
    _ESZ = None
    # ⚠ 两个嵌入资产**必须用同一次布局**（同 Y/gap）：若各自 `solve_geometry`，
    #   两张图的 L5 行位置会因画布高不同而不同 → 主图与独立图对不上（静默不一致）。
    for _path, _no_pct in ((args.embed_asset, args.embed_no_pct),
                           (args.embed_asset_plain, True)):
        if not _path:
            continue
        try:
            ew, eh = (float(v) for v in args.embed_size.lower().split("x"))
        except Exception:
            sys.exit(f"❌ --embed-size 格式应为 `WxH`（如 9.60x4.30），实得 {args.embed_size!r}")
        _ESZ = (ew, eh)
        fig2 = plt.figure(figsize=(ew, eh))
        ax2 = fig2.add_axes([0.012, 0.020, 0.976, 0.955])
        # ⚠⚠ 2026-10-13 21:32：`solve_geometry` 的 gap 地板必须按**实际用到的最大字号**算。
        #   `--all-fs` 同时抬高 L5 → 行距地板随之上抬 → 若仍传 `FS_NODE`，L5 的 37 个
        #   标签会互相压住（**零报错**）。这里取 max(FS_NODE, front_fs, all_fs)。
        _geo_fs = max([FS_NODE] + ([args.front_fs] if args.front_fs > 0 else [])
                      + ([args.all_fs] if args.all_fs > 0 else []))
        Y2, gap2 = solve_geometry(cols, ax2.get_position().height * fig2.get_figheight(),
                                  _geo_fs)
        _lm2 = dict(lmap)
        if _no_pct:
            # 去掉 `(xx.x%)` 后缀，只留物种名（发生率信息仍在**独立图**的注记里）
            _lm2[4] = {nm: nm for nm, _ in cols[4]}
        # ⭐ 2026-10-13：前四层放大（仅嵌入资产）—— 用户「图D的前面 L1-L4 变大些」
        _lfs2 = {}
        _nsc2 = {}
        for _L in range(min(4, len(cols) - 1)):
            if args.front_fs > 0:
                _lfs2[_L] = args.front_fs
            if args.front_node_scale != 1.0:
                _nsc2[_L] = args.front_node_scale
        # ⭐⭐ 2026-10-13 21:32（用户「D图都看不清字体大些」）：**全部五层**统一字号。
        #   根因 = 贴图被缩小 0.729×（资产 10.43 in → 面板 7.61 in），把 5.4 pt 压成 3.94 pt。
        #   本参数覆盖 L0–L4（含末列 species）；`_geo_fs` 已同步抬高行距地板。
        _base_fs = FS_NODE
        if args.all_fs > 0:
            _base_fs = args.all_fs
            for _L in range(len(cols)):          # ⚠ 含 L5 —— 这是与 `--front-fs` 的关键区别
                _lfs2[_L] = args.all_fs
            _nsc2.update({_L: args.front_node_scale for _L in range(min(4, len(cols) - 1))})
        draw(ax2, cols, E, node_color, lambda e: TIER_COLOR[e[3]],
             {0: "left", 1: "right", 2: "right", 3: "right", 4: "right"},
             Y=Y2, xmax=4.78, xmin=-0.40, label_fs=_base_fs, label_map=_lm2,
             bold_labels=bold, gap=gap2, color_of=label_color,
             level_fs=_lfs2, node_fs_scale=_nsc2)
        _pp = Path(_path)
        _pp.parent.mkdir(parents=True, exist_ok=True)
        fig2.savefig(_pp, dpi=300, facecolor="white", bbox_inches="tight",
                     pad_inches=0.02)
        plt.close(fig2)
        if _no_pct:
            embed_png_plain = _pp
        else:
            embed_png = _pp
        print(f"OK  {_pp}  ({ew} x {eh} in 画布，tight 裁边，"
              f"发生率后缀 {'无' if _no_pct else '有'})")
    plt.close(fig)

    # ── 自检 ──
    flows = [sum(n for _, n in c) for c in cols]
    assert all(f == TOTAL for f in flows), f"列流量不守恒（应 {TOTAL}）：{flows}"
    # 末列标签宽度实测（**含发生率后缀**，因为后缀是后加的，最容易把标签顶出画布）
    avail_in = (4.92 - (4 + NW + 0.05)) * (0.958 * 9.32) / (4.92 + 0.80)
    widest_txt = max(l5_label.values(), key=lambda s: _text_w_in(s, FS_NODE))
    w_last = _text_w_in(widest_txt, FS_NODE)
    assert w_last < avail_in, f"末列标签「{widest_txt}」实测宽 {w_last:.3f} in ≥ 可用 {avail_in:.3f} in"
    # 发生率自检：每个 L5 节点都必须拿到值（静默漏 join 会让图上出现空百分比）
    miss_prev = [nm for nm, _ in cols[4] if nm not in prev_max]
    assert not miss_prev, f"{len(miss_prev)} 个 L5 节点缺发生率（例：{miss_prev[:3]}）"
    # ── 嵌入资产自检：末列标签（**无后缀版**）必须放得下 ──
    #   ⚠ 只对**实际产出过**的嵌入资产校验，且用各自的 xlim / 画布宽
    #     （独立图用 (-0.80, 4.92)，嵌入资产用 (-0.40, 4.78)，两者不同）。
    if embed_png_plain is not None and _ESZ:
        _ew = _ESZ[0]
        for _nm, _xmax, _xmin in (("plain", 4.78, -0.40),):
            _av = (4.78 - (4 + NW + 0.05)) * (_ew * 0.976) / (4.78 - (-0.40))
            _wd = max((nm for nm, _ in cols[4]),
                      key=lambda s: _text_w_in(s, FS_NODE))
            _wv = _text_w_in(_wd, FS_NODE)
            assert _wv < _av, \
                f"嵌入资产({_nm})末列标签「{_wd}」{_wv:.3f} in ≥ 可用 {_av:.3f} in"
    print(f"OK  {png}")
    print(f"OK  {pdf}")
    print(f"    列节点数 {[len(c) for c in cols]}  列流量 {flows}")
    print(f"    新种候选（加粗）{len(bold)} 个标签；末列最宽「{widest_txt}」{w_last:.3f} in（可用 {avail_in:.3f}）")
    _pv = sorted(prev_max.values(), reverse=True)
    print(f"    发生率（分母 {DENOM}，L5 取 max）：最高 {_pv[0]:.2f}%  中位 {_pv[len(_pv)//2]:.2f}%  "
          f"最低 {_pv[-1]:.2f}%；<1% 的 {sum(1 for v in _pv if v < 1)}/{len(_pv)} 个")


if __name__ == "__main__":
    main()
