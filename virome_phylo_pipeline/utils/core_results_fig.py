#!/usr/bin/env python3
"""core_results_fig.py — 「核心进化结果」组合主图（Fig4_core_evolution）

把 8 张最有信息量的结果图合成一张，叙事顺序 = 系统发育 → 采样结构 →
重组 → 正选择 → 系统动力学 → 地理传播 → 群体遗传：

    行1  A  ML 系统发育树 (v2)   phylogeny/tree_circular_v2.pdf              源 10.600 x 12.400
         B  采样时间分布         report/sample_plot/sample_temporal.pdf       源  4.740 x  4.267
    行2  C  重组断点分布 (RDP5)  recomb/rdp5/*_publication.pdf               源 16.000 x  6.000
    行3  D  正选择位点-结构映射  select/capheine/figures/Fig_Selection_Structure_Overview.png
                                                                              源 13.333 x 15.360
    行4  E  贝叶斯天际线图       time/beast/postprocess/BSP.pdf               源 10.000 x  6.000
         F  系统地理传播         geography/absorb_test_20260916_full/transmission_map.pdf
                                                                              源  8.801 x  7.486
    行5  G  群体遗传汇总         popgen/popgen_stats_panel.pdf                 源 10.423 x  5.651
         H  单倍型网络           haplo_outputs/gcva_network_location_pub.pdf   源  9.900 x  7.898

⚠ 2026-10-15 第二轮（用户「BSP 和 E 放一行、BSP 在 E 左侧；FG 放一行；AB 放一行；
  CD 单行不变了」）→ 版式由「7 行每行一个面板」改为 **5 行（3 个两栏行 + 2 个单行）**：
  · 面板字母**按阅读顺序整体重排**：旧 E(传播图)→**F**、旧 F(群体遗传)→**G**、旧 G(网络)→**H**；
    新增 **E = BSP**。理由：BSP 插在传播图左侧，沿用旧字母会让图内字母序变成
    `A,B / C / D / H,E / F,G`（字母不按阅读顺序 = 投稿图硬伤）。
    ✅ 已核实无副作用：稿件里的 `Fig. 4A,B` 指的是**另一张图**
    （`metadata_association` 的 Virus–host associations），不是本图，重排不破坏正文引用。
  · **代价（用户明确要求，非缺陷）**：两栏行面板宽 = 4.58 in → 图内字号等比缩小，
    A 树最小字 **2.51 pt**、G 群体遗传 **2.86 pt**、H 网络 **3.24 pt**、
    F 传播图 **4.16 pt**、E BSP **4.58 pt**、B 采样 **4.73 pt**（全部明细见 HALF_COL_MAX_IN 后的注记）。
    单图细节请用桌面独立成品图（如 `7_MLtree_circular_host.png`）。
  · 整图高度 56.85 → **约 32 in**（省 43%）。

⚠ 2026-10-15 第一轮（用户「tree_circular_v2 替换图A；图A 和图B 位置互换；
  补充 popgen_stats 与 gcva_network_location_pub 到主图」）：
  · A 换成 **tree_v2**（源画布 10.6 in，cladogram + 样本号 tip label + 单环外圈=采样地点）；
  · 新增 G/H 两个面板。G 用 `popgen_stats_panel.pdf`（10.42 in）而不是
    `popgen_stats.pdf`（14.81 in）—— 后者即使通栏也只有 4.09 pt（不可读）。
  · ⚠⚠ G 与 H 来自**两套不同口径**（G: n=127 EasyHap；H: n=145 pypopart）→
    题注**必须各自写明 n 与来源**，不得让读者当成同一套分析。

设计纪律（与 `pipeline_overview_fig.py` 完全一致，且**直接复用它的实现**）
────────────────────────────────────────────────────────────────────────
1. **嵌入各阶段原图，不重绘、不降采样**（用户明确要求）。
2. **面板宽必须匹配源画布宽**：源 ≤ 5.2 in → 半栏（w ≈ 4.58）；源 ≥ 8.8 in →
   通栏（w = 9.32）。理由：印刷字号 ≈ 源字号 × (面板宽 / 源画布宽)，把 16 in 宽
   的 RDP5 图或 13.3 in 宽的结构图压进半栏，字号会掉到 4 pt 以下。
   ⚠ **本条在第 2 轮版式（3 个两栏行）里被用户显式让位**：C/D 仍是通栏，
   其余 6 个面板一律按 4.58 in 排 → 字号按 `4.58 / 源画布宽` 缩小（明细见下）。
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
    # ⚠ A 的题注**不与图内自带标题重复**：图内已有
    #   "Maximum-likelihood tree (n = 127 tips, cladogram; 18 zero-length internal
    #    branches collapsed)" 与两条图例（"branch colour = majority host plant
    #   (>=90% purity; grey = mixed)" / "ring cells = sampling location"）——
    #   题注只补图内没有的信息（全基因组比对 + n 与折叠枝数）。
    # ⚠ 历史（2026-10-13，**仅对旧树 tree_circular.pdf 成立，v2 已作废**）：
    #   旧树曾是 "tip 色 + 内环 = 宿主 / 外环 = 采样地点" 的**双环**版，
    #   故旧题注写 "tip colour and inner ring = host plant, outer ring = sampling location"。
    # ⚠ 2026-10-15：A 的源图换成 **tree_v2**（`tree_circular_v2.pdf`，10.600 x 12.400 in）。
    #   v2 的视觉编码与旧树**不同**，题注必须整体重写：
    #     · 旧树：tip 色 + **内环** = 宿主，**外环** = 采样地点（双环）
    #     · v2  ：**branch colour** = majority host plant（>=90% purity，grey = mixed），
    #             外圈**只有一环** = sampling location
    #   实测图内文字（PyMuPDF 提取，非目测）：
    #     "branch colour = majority host plant (>=90% purity; grey = mixed)"
    #     "ring cells = sampling location"
    #     host: Lycium barbarum n=93 / chinense n=3 / ruthenicum n=9（合 105/127）
    #     loc : Yinchuan 66 / Zhongning 15 / Neimenggu 9 / Beijing 8 /
    #           Guangzhou 3 / Wuwei 3 / Xining 1（合 105/127）
    #   ⚠ 地点口径是**城市级**（Yinchuan…），与面板 B 的**省级**口径（Ningxia…）
    #     不同层级 —— 这是两张源图各自的设计，不是数据冲突。
    ("A", "Whole-genome maximum-likelihood cladogram (n = 127 tips; branch lengths "
          "not to scale, 18 zero-length internal branches collapsed); branch colour = "
          "majority host plant (\u2265 90% purity, grey = mixed), outer ring = sampling location",
     "phylogeny/tree_circular_v2.pdf", "pdf"),
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
    # ── 2026-10-15 第二轮（用户「BSP 和 E 放一行，BSP 在 E 左侧；FG 放一行；AB 放一行」）──
    # ⚠⚠ 面板字母**按阅读顺序整体重排**（E/F/G/H 相对上一版各后移一位）：
    #     旧 E(传播图) → 新 **F**；旧 F(群体遗传) → 新 **G**；旧 G(网络) → 新 **H**。
    #   理由：BSP 插在传播图**左侧**，若沿用旧字母，图内字母序会变成
    #   `A,B / C / D / H,E / F,G` —— 字母不按阅读顺序在投稿图里是硬伤。
    #   ✅ 已核实**无副作用**：稿件里的 `Fig. 4A,B` 指的是**另一张图**
    #   （`metadata_association` 的「Virus–host associations」），
    #   不是本图 `Fig4_core_evolution`，故重排不破坏任何正文引用。
    #
    # ⚠ BSP（贝叶斯天际线图）源 `time/beast/postprocess/BSP.pdf`，10.000 x 6.000 in。
    #   **y 轴是对数轴**（不是线性）—— 实测刻度标签像素间距
    #   10→30 = 175.6 pt、30→50 = 81.7 pt，比值 2.149，与 log 理论比
    #   log10(3)/log10(5/3) = 2.151 吻合（**目测会误判成线性，必须量**）。
    #   图内文字：`Bayesian skyline plot` / `95% HPD` / `Posterior median` /
    #   `Time (years before present)` / `Effective population size (Neτ)`。
    ("E", "Bayesian skyline plot (BEAST); log-scale y-axis; line = posterior median, "
          "band = 95% HPD; x = time in years before present",
     "time/beast/postprocess/BSP.pdf", "pdf"),
    ("F", "Phylogeographic spread; arrow colour = age of transmission branch midpoint "
          "(years before present), labels = posterior support",
     "geography/absorb_test_20260916_full/transmission_map.pdf", "pdf"),
    # ── 2026-10-15 第一轮新增（用户「补充 popgen_stats 与 gcva_network_location_pub 到主图」）──
    # ⚠⚠ G 与 H 来自**两套不同口径**（用户裁定「都放 + 题注标明各自 n」）：
    #     G = n = 127（EasyHap, gcva127 集, 按 location 分组, 6 个 n>=3 的群体）
    #     H = n = 145（pypopart, 含 142 单倍型 / 197 边）
    #   两者**不是同一套分析**，题注必须各自写明 n 与工具，否则读者会当成一次运行。
    #
    # ⚠ G 的源图选 **`popgen_stats_panel.pdf`** 而不是用户点名的 `popgen_stats.pdf`：
    #     后者源画布 14.812 in → 通栏缩放 0.629 → 最小字 6.5 pt 印成 **4.09 pt**（不可读）；
    #     `_panel` 版源画布 10.423 in → 缩放 0.894 → 最小字 6.5 pt 印成 **5.81 pt**（可读）。
    #   同源、同内容、同口径（同一 2x3 面板），只是为通栏排版重出过画布。
    #   ⚠ 这一版是**服务器 2026-10-15 03:44** 的重出件（面板 D 补上了 `n/a` 标注，
    #     说明为何只有 Ningxia 有 Fu's Fs 值），本地 `_scratch_tests/` 的 10-13 版没有。
    ("G", "Population genetics by sampling location (n = 127, EasyHap; 6 populations with "
          "n \u2265 3, Qinghai n = 1 excluded): nucleotide diversity \u03c0 vs \u03b8-W, "
          "haplotype diversity Hd, Tajima's D, Fu's Fs (Ewens), pairwise Fst (lower triangle); "
          "global Fst = 0.057 (P = 0.003)",
     "popgen/popgen_stats_panel.pdf", "pdf"),
    ("H", "Haplotype network (n = 145 samples, pypopart; 142 haplotypes, 197 edges, "
          "0 median vectors); node colour = sampling location, node label = haplotype ID, "
          "edge label = number of mutations; scale bar = 1 mutation",
     "haplo_outputs/gcva_network_location_pub.pdf", "pdf"),
]

# 行分组（2026-10-15 第二轮，用户指定）：
#   [A|B] 并排 —— 用户「AB 放一行」；
#   C、D  各自**单行通栏**（用户「CD 单行不变了」）；
#   [E|F] 并排 —— BSP 在传播图**左侧**（用户「BSP 和 E 放一行，放在 E 左侧」）；
#   [G|H] 并排 —— 群体遗传 + 单倍型网络（用户「FG 放一行」）。
# ⚠ 与上一版（7 行每行一个面板）相比：整图高度 56.85 → 约 32 in。
# ⚠ 代价：两栏行里每个面板只有 4.58 in 宽 → 图内字号按 `源画布宽/4.58` 缩小，
#   A(树) 最小字 5.8 pt 印成 **2.51 pt**（详见下方 HALF_COL_MAX_IN 说明后的注记）。
ROWS: List[List[str]] = [["A", "B"], ["C"], ["D"], ["E", "F"], ["G", "H"]]

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

# ── 半栏图「独占一行」时不得被拉伸到通栏宽 ──────────────────────────────
# ⚠⚠ 2026-10-15 第一轮实测（本约束仍保留，用于防「单面板行把窄图拉宽 2 倍」）：
#   布局引擎里 **单面板行** 的面板宽 = 整个可用宽 WU = 9.32 in
#   （`w_panel = (WU - (n-1)*GAP_X)/n`，n=1 时就是 WU）。
#   面板 B 的源画布只有 4.740 in（半栏图）→ 独占一行会被**放大 1.99 倍**，
#   图内字号跟着放大 ~2 倍，与通栏面板（缩放 ~0.9）差一倍以上，肉眼一眼可见。
#   正解 = 源画布 ≤ HALF_COL_MAX_IN 的「半栏图」独占一行时，
#   把面板宽**夹到源画布宽**（缩放 1.00，字号 1:1，符合模块 docstring 第 2 条）。
#   源画布 ≥ HALF_COL_MAX_IN 的通栏图仍填满 WU（轻微上采样 ≤1.06x，无碍）。
# ⚠ 2026-10-15 第二轮（改 2 栏布局）后**没有任何单面板行**，
#   故本约束当前是**空转的保险**（所有行 n=2，w_panel = 4.58 < 5.2 的源宽全部不触发）。
#   保留原因：布局随时可能改回单面板行，届时没有它就会静默把窄图拉宽 2 倍。
HALF_COL_MAX_IN = 5.2

# ⚠⚠ 2026-10-15 第二轮**已知代价**（用户裁定「AB 放一行 / BSP+E 放一行 / FG 放一行」）：
#   两栏行里每个面板宽 = (9.32 - 0.16)/2 = **4.58 in**，
#   图内字号按 `源画布宽 / 4.58` 等比缩小。实测最小字号：
#     A  树 v2    源 10.600 in → 5.80 pt × 0.432 = **2.51 pt**  ← 已低于 5 pt 可读下限
#     B  采样     源  4.740 in → 4.90 pt × 0.966 = 4.73 pt
#     E  BSP      源 10.000 in → 10.0 pt × 0.458 = 4.58 pt
#     F  传播图   源  8.801 in → 8.00 pt × 0.520 = 4.16 pt
#     G  群体遗传 源 10.423 in → 6.50 pt × 0.439 = 2.86 pt
#     H  网络     源  9.900 in → 7.00 pt × 0.463 = 3.24 pt
#   **这是用户明确要求的紧凑版式**（单图细节另有桌面 `7_MLtree_circular_host.png`
#   等独立成品图可看）。若要恢复可读性，正解 = 把该行改回单面板行（`ROWS` 一行一个字母）。

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
        # ⚠ `src_w_in` = 源画布物理宽（英寸）。`rasterize()` 对 PDF 是按 `dpi`
        #   渲染的 → `im.width / dpi` 即源页宽；PNG 则等于「该位图按 `dpi` 输出时
        #   的物理宽」（本图 PNG 源均为 300 dpi 导出，实测对得上）。
        #   用途：判断该面板是不是「半栏图」（见 HALF_COL_MAX_IN）。
        info[letter] = dict(en=en, path=p, img=im, aspect=im.height / im.width,
                            src_w_in=im.width / float(dpi))

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
        # ⚠ 半栏图独占一行 → 夹到源画布宽，禁止放大（理由见 HALF_COL_MAX_IN）。
        #   取本行**最小**源宽做判据：本图所有行都是单面板行，多面板行也不会把
        #   窄图拉宽（保守取小，宁可不放大）。
        _src_w = min(info[x]["src_w_in"] for x in letters)
        if _src_w < HALF_COL_MAX_IN and w_panel > _src_w:
            _log("  [core] 半栏图独占一行: %s 面板宽 %.3f -> %.3f in (源画布宽)"
                 % ("/".join(letters), w_panel, _src_w))
            w_panel = _src_w
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
