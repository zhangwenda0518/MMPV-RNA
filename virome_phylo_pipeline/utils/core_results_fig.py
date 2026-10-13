#!/usr/bin/env python3
"""core_results_fig.py — 「核心进化结果」组合主图（Fig4_core_evolution）

把 5 张最有信息量的结果图合成一张，叙事顺序 = 系统发育 → 采样结构 →
重组 → 正选择 → 地理传播：

    A  ML 系统发育树          phylogeny/tree_circular.pdf                   半栏 4.580 x 4.580
    B  采样时间分布           report/sample_plot/sample_temporal.pdf        半栏 4.691 x 4.223
    C  重组断点分布 (RDP5)    recomb/rdp5/*_publication.pdf                 通栏 16.000 x 6.000
    D  正选择位点-结构映射    select/capheine/figures/Fig_Selection_Structure_Overview.png
                                                                            通栏 13.333 x 15.360
    E  系统地理传播           geography/absorb_test_20260916_full/transmission_map.pdf
                                                                            通栏 8.801 x 7.486

设计纪律（与 `pipeline_overview_fig.py` 完全一致，且**直接复用它的实现**）
────────────────────────────────────────────────────────────────────────
1. **嵌入各阶段原图，不重绘、不降采样**（用户明确要求）。
2. **面板宽必须匹配源画布宽**：源 ≤ 5.2 in → 半栏（w ≈ 4.58）；源 ≥ 8.8 in →
   通栏（w = 9.32）。理由：印刷字号 ≈ 源字号 × (面板宽 / 源画布宽)，把 16 in 宽
   的 RDP5 图或 13.3 in 宽的结构图压进半栏，字号会掉到 4 pt 以下。
3. **题注必须量宽折行**（`_text_w_in` + `_wrap_caption`），行高按**该行最长题注
   的行数**计入 H —— 否则折行会压到下一行面板字母。
4. **rcParams 用唯一真源** `sci_palette.sci_rcparams()`，不许手抄（手抄必漏键，
   漏键必顺序依赖 → 同一份代码两次跑出不同的图）。
5. **PDF 必须显式传 dpi**：不传就吃全局 `savefig.dpi`，输出会取决于同进程里
   是否先跑过别的模块（`apply_sci_style` 会把它改成 300）。

⚠ 面板几何改动后**必须重跑版式体检**（判据 0 对重叠 + 无 TEXT-OFF-CANVAS）。

用法
────
    python utils/core_results_fig.py --root <data_dir> --out <dir> --title "病毒名"
"""
from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

# 复用 Fig3 的排版机制 —— 题注量宽/折行、光栅化、路径解析、RDP5 现生成。
# 不重复实现：这些函数的边界条件（TextPath 度量、布局换算、RDP5 参数）都在
# Fig3 那边踩过坑并固化过，抄一份必然漂移。
try:                                       # 包内导入
    from utils.pipeline_overview_fig import (
        _recomb_rdp5_dirs, _text_w_in, _wrap_caption, maybe_build_rdp5_pdf,
        rasterize, resolve_panel_path,
    )
    from utils.sci_palette import sci_rcparams
except ImportError:                        # 脚本目录内直接运行
    from pipeline_overview_fig import (
        _recomb_rdp5_dirs, _text_w_in, _wrap_caption, maybe_build_rdp5_pdf,
        rasterize, resolve_panel_path,
    )
    from sci_palette import sci_rcparams

Image.MAX_IMAGE_PIXELS = None

# ── 字体：图内**全英文**，直接用 DejaVu Sans（与 Fig3 一致）────────────
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.unicode_minus"] = False

# ── 面板登记表 ────────────────────────────────────────────────────────
# (字母, 英文题注, 相对路径, 类型)
# ⚠ rel 是 **legacy 字面量**（第一段 = 模块名），取值时经 resolve_panel_path()
#   按当前 I/O 布局换算（standard: report→08_report 等），legacy 名兜底。
# ⚠ 题注是自足句：写明每个视觉编码的含义 + 关键统计量，读者不必回正文。
PANELS: List[tuple] = [
    # ⚠ A 的题注**不与图内自带标题重复**（图内已有 "Maximum-likelihood tree
    #   (n = 127 tips; colour = host plant)" 与 "105/127 tips assigned to
    #   3 host plant"）—— 题注只补图内没有的信息：全基因组比对 + 外圈双环含义。
    # ⚠ 2026-10-13 着色口径由「采样地点」改为「宿主植物」（用户指定
    #   「分支展示不同物种的，不按地区」）→ 题注同步，勿漏改。
    # ⚠ 2026-10-13 再改：树枝改为 **cladogram（不展示枝长）** → 原题注里的
    #   "scale bar = substitutions/site" 已不成立（该图**没有标尺**）；
    #   同时外圈新增**双环色带**（内=宿主 / 外=采样地点），一并写进题注。
    ("A", "Whole-genome maximum-likelihood cladogram of GCVA (branch lengths "
          "not shown); tip colour and inner ring = host plant, outer ring = sampling location",
     "phylogeny/tree_circular.pdf", "pdf"),
    # ⚠ 2026-10-13 题注口径更新（用户「图 B 把物种也加上」）：气泡内新增宿主饼
    #   → 原题注只提 "bubble area scales with count"，未说明饼的含义。
    #   **必须与 Fig3 面板 B 逐字同步**（同一张源图，两处题注口径要一致）。
    ("B", "Sampling through time; line = yearly total; bubble area scales with count, "
          "pie slices = host plant composition",
     "report/sample_plot/sample_temporal.pdf", "pdf"),
    # ⚠ RDP5 断点图：优先用树内已有的**排版版**（`*_publication.pdf`，16 x 6 in）；
    #   没有才按排版尺寸现生成（见 maybe_build_rdp5_pdf 的说明）。
    ("C", "Recombination breakpoint distribution: orange = RBN (200-nt window), "
          "blue = 99% CI, black line = RBN, boxes = gene map",
     "@RDP5@", "abs"),
    ("D", "Positive selection (HYPHY FEL/MEME): 10 sites in 4 of 6 genes mapped onto "
          "AlphaFold3 models — orange = pervasive (FEL, P < 0.05), "
          "blue = episodic (MEME only, P < 0.05)",
     "select/capheine/figures/Fig_Selection_Structure_Overview.png", "png"),
    ("E", "Phylogeographic spread; arrow colour = age of transmission branch midpoint "
          "(years before present), labels = posterior support",
     "geography/absorb_test_20260916_full/transmission_map.pdf", "pdf"),
]

# 行分组：窄图并排，宽图独占通栏。
#   A/B 源画布 4.58 / 4.69 in → 半栏并排（缩放 1.00 / 0.98，字号不缩水）；
#   C/D/E 源画布 16.0 / 13.33 / 8.80 in → 各占通栏。
ROWS: List[List[str]] = [["A", "B"], ["C"], ["D"], ["E"]]

# ── 版面几何（英寸） ──
W_FIG = 10.0
M_L, M_R = 0.40, 0.28
M_T, M_B = 1.00, 0.30
WU = W_FIG - M_L - M_R
GAP_X = 0.16
ROW_GAP = 0.30        # 行间距（无组头，行与行之间留白）
H_CAP = 0.42          # 单行题注行高
CAP_FS = 6.8          # 题注字号（必须与 fig.text(fontsize=...) 一致）
CAP_LINE_H = 0.115    # 每多一行题注增加的高度（6.8 pt x 1.2 / 72 ≈ 0.113）
CAP_TOP_OFF = 0.10    # 题注顶端相对面板下沿的偏移
CAP_WRAP_TOL = 0.98
DPI = 300

DEFAULT_STEM = "Fig4_core_evolution"


def build_core_results(root: str, out_dir: str, *, title: str = "",
                       rdp5_pdf: Optional[str] = None,
                       stem: str = DEFAULT_STEM,
                       dpi: int = DPI, log=None) -> Dict[str, Optional[str]]:
    """生成「核心进化结果」组合主图。返回 {'png':..,'pdf':..}（失败项为 None）。"""
    def _log(m, lvl="emit"):
        if log is not None:
            getattr(log, lvl, log.emit)(m)
        else:
            print(m)

    # rcParams 用唯一真源钉死（理由见模块 docstring 第 4 条）。
    #   despine_top_right=False：每个面板带一圈细灰边框（下面会给四条边上色）。
    matplotlib.rcParams.update(sci_rcparams(base_size=8, grid=False,
                                            despine_top_right=False))

    root = os.path.abspath(os.path.expanduser(root))
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    if not stem:
        stem = DEFAULT_STEM

    # ── RDP5 面板：树内已有的排版版优先，否则现生成 ──
    if rdp5_pdf is None:
        import glob
        pre = []
        for rd in _recomb_rdp5_dirs(root):
            pre += glob.glob(os.path.join(rd, "breakpoint_distribution_publication.pdf"))
            pre += glob.glob(os.path.join(rd, "*_publication.pdf"))
        rdp5_pdf = sorted(pre)[0] if pre else None
    if rdp5_pdf is None:
        rdp5_pdf = maybe_build_rdp5_pdf(root, os.path.join(out_dir, "_core_src"),
                                        title or "Virus", log=log)

    # ── 收图 ──
    info = {}
    for letter, en, rel, kind in PANELS:
        p = rdp5_pdf if rel == "@RDP5@" else resolve_panel_path(root, rel)
        if not p or not os.path.exists(p):
            _log(f"  [core] 跳过面板 {letter}: {rel}")
            continue
        try:
            im = rasterize(p, dpi)
        except Exception as e:                       # noqa: BLE001
            _log(f"  [core] 跳过面板 {letter}: 读图失败 {e}", "warning")
            continue
        info[letter] = dict(en=en, path=p, img=im, aspect=im.height / im.width)

    if len(info) < 3:
        _log("  [core] 可用面板不足 3 个 → 不出图", "warning")
        return {"png": None, "pdf": None}

    # ── 分行 + 题注折行 ──
    layout = []
    for letters in ROWS:
        letters = [x for x in letters if x in info]
        if not letters:
            continue
        n = len(letters)
        w_panel = (WU - (n - 1) * GAP_X) / n
        h_row = max(w_panel * info[x]["aspect"] for x in letters)
        caps = [_wrap_caption(info[x]["en"], w_panel) for x in letters]
        n_cap = max(len(c) for c in caps)
        if n_cap > 1:
            _log("  [core] 题注折行(栏宽 %.2f in): %s" % (
                w_panel, ", ".join("%s->%d行" % (ltr, len(c))
                                   for ltr, c in zip(letters, caps) if len(c) > 1)))
        layout.append((letters, w_panel, h_row, caps, n_cap))

    H = M_T + M_B + ROW_GAP * max(0, len(layout) - 1)
    for _letters, _wp, h_row, _caps, n_cap in layout:
        H += h_row + H_CAP + (n_cap - 1) * CAP_LINE_H

    _log(f"  [core] 设计尺寸 {W_FIG:.2f} x {H:.2f} in "
         f"({len(layout)} 行 / {len(info)} 面板)")

    fig = plt.figure(figsize=(W_FIG, H), dpi=150)
    fig.patch.set_facecolor("white")
    fx = lambda v: v / W_FIG   # noqa: E731
    fy = lambda v: v / H       # noqa: E731

    _t = title or "Virus"
    fig.text(fx(M_L), fy(H - 0.42),
             f"{_t} — core evolutionary results",
             fontsize=13, fontweight="bold", ha="left", va="center")
    fig.text(fx(M_L), fy(H - 0.72),
             "Phylogeny, sampling structure, recombination, positive selection "
             "and phylogeographic spread",
             fontsize=9.5, color="#444444", ha="left", va="center")

    y_top = H - M_T
    for _i, (letters, w_panel, h_row, caps, n_cap) in enumerate(layout):
        if _i:
            y_top -= ROW_GAP
        y_pt = y_top - h_row
        for j, letter in enumerate(letters):
            d = info[letter]
            x0 = M_L + j * (w_panel + GAP_X)
            h = w_panel * d["aspect"]
            ax = fig.add_axes([fx(x0), fy(y_pt), fx(w_panel), fy(h)])
            ax.imshow(d["img"], interpolation="antialiased")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_edgecolor("#bbbbbb"); s.set_linewidth(0.6)
            ax.set_anchor("N")
            ax.text(0.0, 1.005, letter, transform=ax.transAxes,
                    fontsize=11, fontweight="bold", ha="left", va="bottom")
            fig.text(fx(x0), fy(y_pt - CAP_TOP_OFF), "\n".join(caps[j]),
                     fontsize=CAP_FS, ha="left", va="top", color="#222222",
                     linespacing=1.2)
        y_top = y_pt - H_CAP - (n_cap - 1) * CAP_LINE_H

    out_png = os.path.join(out_dir, stem + ".png")
    out_pdf = os.path.join(out_dir, stem + ".pdf")
    # ⚠ dpi 必须显式传（理由见模块 docstring 第 5 条）。
    fig.savefig(out_png, dpi=dpi, facecolor="white",
                bbox_inches="tight", pad_inches=0.12)
    fig.savefig(out_pdf, dpi=dpi, facecolor="white",
                bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    _log(f"  [core] -> {out_pdf} ({os.path.getsize(out_pdf)} B)")
    _log(f"  [core] -> {out_png} ({os.path.getsize(out_png)} B)")
    return {"png": out_png, "pdf": out_pdf}


def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="Core evolutionary results figure (Fig4_core_evolution)")
    ap.add_argument("--root", required=True, help="data_dir (含各阶段目录)")
    ap.add_argument("--out", required=True, help="输出目录")
    ap.add_argument("--title", default="", help="图标题里的病毒名")
    ap.add_argument("--rdp5-pdf", default=None, help="现成的 rdp5 断点图 PDF")
    ap.add_argument("--stem", default=DEFAULT_STEM)
    ap.add_argument("--dpi", type=int, default=DPI)
    a = ap.parse_args()
    r = build_core_results(a.root, a.out, title=a.title, rdp5_pdf=a.rdp5_pdf,
                           stem=a.stem, dpi=a.dpi)
    if not r.get("pdf"):
        sys.exit(1)


if __name__ == "__main__":
    main()
