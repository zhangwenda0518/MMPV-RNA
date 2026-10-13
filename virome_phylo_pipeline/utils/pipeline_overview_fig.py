#!/usr/bin/env python3
"""pipeline_overview_fig.py — 「系统发育分析管线各阶段代表性结果」组合主图（Fig3）

给 `publication_figures.build_all_publication_outputs()` 调用，也可独立跑：

    python utils/pipeline_overview_fig.py --root <data_dir> --out <dir>

设计纪律（沿用本项目既有约定）
────────────────────────────────
1. **嵌入各阶段原图，不重绘、不降采样**（用户明确要求）。
2. **按源图物理宽度分行**：宽图给通栏整行、窄图并排。理由：印刷字号 ≈ 源字号 ×
   (面板宽 / 源画布宽)，把 9–19 in 的源图压进半栏会掉到 <5 pt。
3. **只收有信息量的图**。经逐张目视确认后剔除退化图（单柱 / 单圆）：
   `seqgroup_location` / `seqgroup_host`（全部序列同地点同宿主 → 一个大圆）、
   `align_qc`（所有序列等长同 gap → 两根针）、`popgen_metrics`（两根柱）。
4. PDF 用 **PyMuPDF 光栅化到 300 dpi** 再嵌入（matplotlib 不能直接放 PDF 页）。

⚠ 面板几何改动后**必须重跑版式体检**（`_check_fig3_overlap.py`，判据 0 overlapping pairs）：
面板字母画在面板上沿外侧，`H_CAP` 必须 ≥ 题注高 + 字母高 + 余量，否则字母会压上一行题注。
"""
from __future__ import annotations

import os
import subprocess
import sys
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from PIL import Image

try:                                       # 包内导入
    from utils.sci_palette import sci_rcparams
except ImportError:                        # 脚本目录内直接运行
    from sci_palette import sci_rcparams

# ── 布局感知 (2026-10-12 修 P0-1) ─────────────────────────────────────
# PANELS 的 rel 全是 **legacy** 模块目录名 (report/、phylogeny/…), 而 standard
# (默认布局) 下模块目录叫 08_report/、02_phylogeny/… —— 旧实现直接
# os.path.join(root, rel), 于是 23 个面板**全部**"跳过", Fig3 永不产出。
# 现统一经 ph_dir() 换算第一段模块名: standard 名优先、legacy 字面量兜底
# (存量 legacy 结果树兼容)。
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import ph_dir  # noqa: E402


def _module_candidates(root: str, mod: str) -> List[str]:
    """模块名 → 候选目录绝对路径 (standard 名优先, legacy 同名字面量兜底)。"""
    cands: List[str] = []
    try:
        cands.append(os.path.join(root, ph_dir(mod)))
    except (KeyError, ValueError):
        pass
    legacy = os.path.join(root, mod)
    if legacy not in cands:
        cands.append(legacy)
    return cands


def resolve_panel_path(root: str, rel: str) -> Optional[str]:
    """面板相对路径 (PANELS 表 legacy 字面量) → 实际存在的绝对路径; 均不存在返回 None。"""
    parts = rel.split("/")
    cands: List[str] = []
    try:
        cands.append(os.path.join(root, ph_dir(parts[0]), *parts[1:]))
    except (KeyError, ValueError):
        pass
    legacy = os.path.join(root, rel)
    if legacy not in cands:
        cands.append(legacy)
    for p in cands:
        if os.path.exists(p):
            return p
    return None


def _recomb_rdp5_dirs(root: str) -> List[str]:
    """recomb/rdp5 的候选目录 (standard 优先, legacy 兜底)。"""
    return [os.path.join(d, "rdp5") for d in _module_candidates(root, "recomb")]


Image.MAX_IMAGE_PIXELS = None

# ── 字体：**图内全英文**（用户 2026-10-12 要求），故直接 DejaVu Sans。
# 原实现为"有 CJK 就用 Noto Sans CJK"以支持中英混排；中英混排标题/组头/题注
# 已全部改英文，用 CJK 字体的拉丁字形反而不如 DejaVu Sans 规整。
_PRIMARY = "DejaVu Sans"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.unicode_minus"] = False

# ── 面板登记表 ─────────────────────────────────────────────
# (字母, 组id, 英文题注, 相对路径, 类型)
# 2026-10-12: 图内文字全英文 → 原双语题注的中文键已删除（确认从未绘制）。
#
# ⚠ rel 是 **legacy 字面量** (第一段 = 模块名), 取值时经 resolve_panel_path()
#   按当前 I/O 布局换算 (standard: report→08_report 等), legacy 名兜底 —— 改面板
#   路径请只动模块名之后的部分, 不要写死带编号的目录名。
#
# 2026-10-12（审稿人视角）：英文题注原为"Sampling locations"这类纯名词短语，
# 读者必须回到正文才能知道图上每个视觉编码代表什么。现改为**自足题注** ——
# 写明编码含义（颜色/大小/线型）与关键统计量（n、检验、置信区间）。
# ⚠ 题注是 `fig.text` 单行不换行：2 栏行的题注必须 ≲ 80 字符，
#   否则会越过列缝压到右列题注上（版式体检会报 overlap）。
#
# ── 2026-10-12（第 5 轮）：第 5–8 组重排 ───────────────────────────────
# 用户反馈「5 6 7 8 核心结果 U 算一张，其他的好多你再好好查查」。
# 逐张实测源画布宽度 + PDF 内字号（`page.get_text("dict")` 的 span size）后定位到
# **真正的根因是版式而非内容**：第 5–8 组有 9 个面板的源画布是 8.8–11.5 in 宽，
# 却被塞进 4.58 in 的半栏 → 印刷字号 = 源字号 × (4.58 / 源画布宽)：
#
#   面板 源画布  源最小字号  旧缩放  旧印刷字号 | 通栏缩放 通栏印刷字号
#   J RTT    8.99 in   8.4 pt   0.51   4.3 pt   |  1.04   8.7 pt
#   M LTT   10.00      10.0     0.46   4.6 pt   |  0.93   9.3 pt
#   N dated 10.00       7.0     0.46   3.2 pt ← 不可读
#   P BSP   10.00      10.0     0.46   4.6 pt   |  0.93   9.3 pt
#   R geneG 10.27      10.5     0.45   4.7 pt   |  0.91   9.5 pt
#   U Mantel 9.09       9.0     0.50   4.5 pt   |  1.03   9.2 pt
#   V trans 8.80        8.0     0.52   4.2 pt   |  1.06   8.5 pt
#   W virsp 9.32        7.0     0.49   3.4 pt ← 不可读
#
# 处置（只改**放置方式**，不改任何源图；被移出的图仍留在结果树里）：
#   • J / K / L(dated) / P(BSP) / Q(QQ) / U(Mantel) / V(transmap) → 各占**通栏**；
#   • ESS-BEAST 与 ESS-phylogeo 源画布都是 4.58 in（= COL_W_2），半栏并排**零缩放**，
#     且两者直接可比（同一种诊断、两次 run）；
#   • **删** RTT-robust（与 J 是同一张散点、同一结论 R²≈0.02）、
#     LTT（与 BSP 同答"谱系/Ne 随时间怎么变"）、
#     VirSpaceTime（把 A 的空间分布与 B 的时间分布再画一遍）。
# 结果：第 5–8 组由 14 张 → 11 张，且没有一张印在 6 pt 以下。
PANELS: List[tuple] = [
    ("A", "g1", "Sampling locations; bubble area scales with isolate count",
     "report/sample_plot/sample_map.pdf", "pdf"),
    ("B", "g1",
     # ⚠ 2026-10-13 题注口径更新（用户「图 B 把物种也加上」）：气泡内新增宿主饼
     #   → 原题注只提 "bubble area scales with count"，未说明饼的含义，必须补。
     #   同时图例已从坐标区下方移入**两面板之间的空隙**，题注不必再提位置。
     "Sampling through time; line = yearly total; bubble area scales with count, "
     "pie slices = host plant composition",
     "report/sample_plot/sample_temporal.pdf", "pdf"),
    # ⚠ 2026-10-13 着色口径由「采样地点」改为「宿主植物」（用户指定）
    #   同时树枝改为 cladogram（不展示枝长，见 stage_plots::plot_tree_circular）。
    ("C", "g2", "Maximum-likelihood cladogram (n = 127 tips; branch lengths not "
                "shown); tips coloured by host plant",
     "phylogeny/tree_circular.pdf", "pdf"),
    # 2026-10-12 新增：管线一直有 `splitstree` 阶段（SplitsTree6 NeighborNet），
    # 但产物 nexus 从未被画出来 —— 主图缺"数据到底是不是树状"这一环。
    # ⚠ 第 7 轮两处修正：
    #   ① 题注必须**折行**：原先单行 148 字符 → 实测右缘 11.643 in >
    #      设计 9.720 in，`bbox_inches="tight"` 把画布从 10.00 撑到 11.42 in
    #      （详见 CAP_* 常量块）。
    #   ② 题注**不与图内自带的题注重复**：`neighbornet.pdf` 自己底部已有
    #      "n = 145 taxa, 243 splits, fit = 99.2% / 184/29,403 split pairs …"，
    #      面板题注改成只写**结论 + 图例含义**，避免同一句话在图里出现两遍。
    # ⚠ 第 8 轮：D 面板改用 **SplitsTree6 自己算的布局**（`neighbornet_panel.pdf`，
    #   由 splitstree_native_viz 按**半栏 4.58 in** 出图、无标题/无图例/无题注、
    #   不画 taxon 标签 —— 符合「面板宽必须匹配源画布宽」，字号不再打 5 折）；
    #   旧的整页 `neighbornet.pdf`（9.32 in 通栏 + tip 标签 + 图例）仍作为独立图保留。
    # ⚠ 第 9 轮：面板 D 的源图换成 **449 节点纯 EqualAngle 辐射布局**（旧版是
    #   514 节点 EqualAngle+ConvexHull 的宽扁"彗星"），并**自绘分组弧线 + 组名**
    #   （SplitsTree6 本体没有分组轮廓功能，已穷举源码确认 → 必须自绘）。
    #   源图尺寸 4.580 x 4.566 in（aspect 0.997，近方形；旧版 4.580 x 2.829 in）。
    #   题注恢复提 "arcs"（本轮图内确实画了分组弧线），并写清颜色/弧线的共同含义。
    ("D", "g2",
     # ⚠ 2026-10-13 着色口径由「采样地点」改为「宿主植物」（与面板 C 统一）：
     #   管线 `--splitstree_group_col` 默认值已同步改为 host。
     "Neighbor-Net (SplitsTree6 native layout): near tree-like genealogy "
     "(fit 99.2%); tip colour and outer arcs = host plant",
     "phylogeny/splitstree/neighbornet_panel.pdf", "pdf"),
    ("E", "g2", "Substitution saturation: Ti/Tv ratio vs p-distance; no saturation (R² = 0.005)",
     "phylogeny/saturation_out/saturation_plot.png", "png"),
    ("F", "g3", "Linkage disequilibrium decay; left = all 1.60 M SNP pairs, right = binned means",
     "popgen/dnasp/figures/ld_decay.png", "png"),
    # 2026-10-13（用户: 「面板 G 新绘制」）: 源图由 dnasp summary（仅全局 π/Hd）
    #   换为 `popgen/popgen_stats_panel.png` —— utils/popgen_plots 的 2×3 多面板
    #   （π vs θ-W(%) / Hd / Tajima's D / Fu's Fs / pairwise Fst 下三角热图 / 摘要）。
    #   数据口径: gcva127.fasta（= G 基因树 127 tips ∩ 全长比对, 与历史报告/论文
    #   n=127 同集）; Fu's Fs 为 Ewens 口径（utils/popgen_plots docstring）;
    #   全局 Fst 取自历史 popgen_report.txt（127 集, p=0.003）—— MJN 在该数据上
    #   中介向量爆炸（实测 348GB 不收敛）, 故 fst 网络路径跳过, 数值沿用历史。
    #   源画布 9.32 × 5.40 in = 通栏宽 → 字号 1:1, 符合「面板宽匹配源画布宽」规则。
    ("G", "g3", "Population genetics by location (n = 127): π vs θ-W, Hd, Tajima's D "
     "(Nei 1987), Fu's Fs, pairwise Fst (Hudson); global Fst = 0.057 (P = 0.003)",
     "popgen/popgen_stats_panel.png", "png"),
    ("H", "g3", "Host differentiation; boxplots and permutation test (9,999 permutations)",
     "popgen/host_analysis/host_differentiation.pdf", "pdf"),
    ("I", "g4",
     "Recombination breakpoint distribution: orange = RBN (200-nt window), "
     "blue = 99% CI, black line = RBN, boxes = gene map",
     "@RDP5@", "abs"),
    # ── g5 时间信号与分子钟（3 张，全部通栏） ──
    ("J", "g5",
     "Root-to-tip regression (TreeTime): rate 1.59×10⁻⁴ subs/site/yr, "
     "R² = 0.019, P = 0.156 — no detectable temporal signal",
     "time/treetime_rtt/RootToTip_regression.pdf", "pdf"),
    ("K", "g5",
     "Date-randomization test (n = 20): observed R² and rate fall inside the "
     "randomized-date null — confirms the negative RTT result",
     "time/temporal_signal/drt_results.pdf", "pdf"),
    ("L", "g5",
     "Time-calibrated phylogeny (TreeDater); node bars = 95% CI from parametric bootstrap",
     "time/treedater_ltt/Phylogeny_dated.pdf", "pdf"),
    # ── g6 贝叶斯系统动力学（4 张 / 3 行） ──
    ("M", "g6",
     "MCMC ESS — BEAST: kappa = 130, the other five parameters < 20",
     "time/beast/merged/mcmc_diagnostics.pdf", "pdf"),
    ("N", "g6",
     "MCMC ESS — phylogeography: 4 of 6 parameters > 200",
     "geography/phylogeography/merged/mcmc_diagnostics.pdf", "pdf"),
    ("O", "g6", "Bayesian skyline plot; log-scale y-axis, band = 95% HPD",
     "time/beast/postprocess/BSP.pdf", "pdf"),
    ("P", "g6",
     "Between-chain agreement (QQ): points on the dashed line = identical marginal "
     "posteriors, so the low-ESS chains still agree",
     "time/beast/merged/chain_scatter.pdf", "pdf"),
    # ── g7 基因定年（2 张并排） ──
    ("Q", "g7", "Gene dating — G gene; time-calibrated tree (substitutions/site)",
     "phylogeny/gene_analysis/trees/G/TimeTree.pdf", "pdf"),
    ("R", "g7", "Gene dating — L gene; time-calibrated tree (substitutions/site)",
     "phylogeny/gene_analysis/trees/L/TimeTree.pdf", "pdf"),
    # ── g8 系统地理学与传播（2 张，全部通栏） ──
    ("S", "g8",
     "Mantel test: r = −0.048 (r² = 0.002, P = 0.0002) — no isolation-by-distance",
     "geography/geo_analysis/mantel_test.pdf", "pdf"),
    ("T", "g8",
     "Phylogeographic spread; arrow colour = age of transmission branch midpoint, "
     "labels = posterior support",
     "geography/absorb_test_20260916_full/transmission_map.pdf", "pdf"),
    # ── g9 正选择（1 张通栏总览）────────────────────────────────────────
    # 2026-10-12 第 5 轮补漏：`STAGE_GROUPS` 里的 `select`(capheine) 这一组
    # **此前没有任何面板** —— 用户指出「缺少正选择位点分析的结果」。
    # 管线实际跑了 HYPHY 8 种方法（BUSTED/CLN/CONTRASTFEL/FEL/LABELTREE/MEME/
    # PRIME/RELAX）× 6 个基因，但**只落 JSON，从不画图**；位点级结果
    # （`selection_sites.tsv`：10 个位点 / 4 个基因）由
    # `reports/_sel_struct_20261012/` 的 PyMOL 渲染链映射到 AlphaFold3 结构上。
    #
    # 第 6 轮（用户「缺少 Fig_Selection_Structure_Overview.png」）：原先 4 张
    # 单基因图（各 15.0 × 6.07 in，共吃 15.1 in 高）换成 **1 张四基因 2×2 总览**
    # —— 总览在每个蛋白上直接点名全部 10 个位点，信息不减、高度省 4.3 in。
    # 源图 13.33 × 15.36 in（aspect 1.152）→ 通栏缩放 9.32/13.33 = 0.70 →
    # 图内 11.5 pt 题注印成 **8.1 pt**（可读）；半栏只有 4.0 pt，故必须通栏。
    # ⚠ 4 张单基因图 + `Fig_Selection_Structure_Zooms.png` 仍留在
    #   `select/capheine/figures/` 作补充图（管线本身不产它们）。
    ("U", "g9",
     "Positive selection (HYPHY FEL/MEME): 10 sites in 4 of 6 genes — "
     "orange = pervasive (FEL), blue = episodic (MEME only)",
     "select/capheine/figures/Fig_Selection_Structure_Overview.png", "png"),
]

GROUPS = [
    ("g1", "Sampling & metadata"),
    ("g2", "Alignment & phylogeny"),
    ("g3", "Population genomics & host"),
    ("g4", "Recombination"),
    ("g5", "Temporal signal & molecular clock"),
    ("g6", "Bayesian phylodynamics"),
    ("g7", "Gene dating"),
    ("g8", "Phylogeography & spread"),
    ("g9", "Positive selection"),
]

# ⚠ 分行规则（2026-10-12 第 5 轮固化）：**面板宽必须匹配源画布宽**。
#   源画布 ≤ 5.2 in → 半栏（w≈4.58 in，缩放 ≈0.9–1.0）；
#   源画布 ≥ 8.8 in → 通栏（w≈9.32 in，缩放 ≈0.9–1.1）。
#   把 9–11.5 in 的源图放进半栏 = 字号打 45 折（实测 3.2–4.7 pt），是本轮返工的根因。
ROWS = [
    ("g1", ["A", "B"]), ("g2", ["C", "D"]), ("g2", ["E"]), ("g3", ["F"]), ("g3", ["G"]),
    ("g3", ["H"]), ("g4", ["I"]),
    ("g5", ["J"]), ("g5", ["K"]), ("g5", ["L"]),
    ("g6", ["M", "N"]), ("g6", ["O"]), ("g6", ["P"]),
    ("g7", ["Q", "R"]),
    ("g8", ["S"]), ("g8", ["T"]),
    # g9 正选择：源画布 13.33 × 15.36 in（aspect 1.152）→ 占**通栏**。
    # 通栏 w = 9.32 in → 缩放 0.70 → 图内 11.5 pt 题注印成 8.1 pt。
    ("g9", ["U"]),
]

# ── 版面几何（英寸） ──
W_FIG = 10.0
M_L, M_R = 0.40, 0.28
M_T, M_B = 1.00, 0.30
WU = W_FIG - M_L - M_R
GAP_X = 0.16
GAP_GROUP = 0.28
H_HEADER = 0.50      # 组头 -0.16 / 细线 -0.28 / 面板顶 -0.50
H_CAP = 0.42         # **单行**题注行的固定高度（题注 -0.10~-0.213 + 下一行字母 0.155 + 余量）
DPI = 300

# ── 面板题注（2026-10-12 第 7 轮：从"单行硬写"改为"按栏宽自动折行"）────────────
# ⚠ 事故：D 面板题注写成单行 148 字符 → 实测右缘 **11.643 in**，而设计右缘
#   `W_FIG - M_R = 9.720`。`bbox_inches="tight"` 于是把整张画布从 10.00 in
#   撑到 **11.42 in**（PDF page rect 实测），右列题注甩出图外、版面失去对齐。
#   修法 = 量宽后折行（不许再靠"目测差不多"）：
#     · `_text_w_in()` 用 `TextPath` 取**真实字形宽度**（同字号同字体，非估算）；
#     · `_wrap_caption()` 折到 `w_panel * CAP_WRAP_TOL` 以内；
#     · 行高按**该行最长题注的行数**计入 `H`（`H_CAP + (n-1) * CAP_LINE_H`），
#       而不是全局写死 —— 否则折行会压到下一行面板字母。
CAP_FS = 6.8         # 题注字号（必须与下面 fig.text(fontsize=...) 一致）
CAP_LINE_H = 0.115   # 单行题注占高（6.8 pt × 1.2 / 72 ≈ 0.113 in）
CAP_TOP_OFF = 0.10   # 题注顶端相对面板下沿的偏移（与下面 fig.text 的 -0.10 一致）
CAP_WRAP_TOL = 0.98  # 折行留 2 % 余量（TextPath 与 hinting 的细微差）

# rdp5 排版参数（见 rd p5_bdp_plot.py 的 --width/--height/--font-scale）
RDP5_W, RDP5_H = 16.0, 6.0
RDP5_FONT, RDP5_GENE_FONT = 1.0, 12.0


# ═══════════════════════════════════════════════════════════════════
# rdp5 断点图（缺图时按排版尺寸现生成）
# ═══════════════════════════════════════════════════════════════════
def _find_rdp5_inputs(root: str):
    """找 (breakpoint csv, positions csv, gene-coords csv)。找不到返回 None。

    2026-10-12 修 (P0-1): rdp5 目录经布局换算 (standard=04_recomb/rdp5,
    legacy=recomb/rdp5); 基因坐标候选 data/ 同样两布局都试。
    """
    import glob
    rd = next((d for d in _recomb_rdp5_dirs(root) if os.path.isdir(d)), None)
    if rd is None:
        return None
    bp = [p for p in glob.glob(os.path.join(rd, "*breakpoint distribution.csv"))
          if not p.endswith("BreakpointPositions.csv")]
    pos = glob.glob(os.path.join(rd, "*BreakpointPositions.csv"))
    if not bp or not pos:
        return None
    # 基因坐标: 优先 root 下的 *genes.csv (gene,start,end,type)，其次 data/ 下同名
    cands = list(glob.glob(os.path.join(root, "*genes.csv")))
    for d in _module_candidates(root, "data"):
        cands += glob.glob(os.path.join(d, "*genes.csv"))
    if not cands:
        return None
    return sorted(bp)[0], sorted(pos)[0], sorted(cands)[0]


def _orfcoords_from_genes(genes_csv: str, out_csv: str) -> bool:
    """`gene,start,end,type` → RDP5_RBDP_Rgrapher 要的 `Start,Stop,Gene symbol`。

    同一基因可能同时有 gene/CDS 两行 → 去重（按 gene 名取第一次出现）。
    """
    import csv as _csv
    seen, rows = set(), []
    try:
        with open(genes_csv, newline="") as f:
            for r in _csv.DictReader(f):
                g = (r.get("gene") or r.get("Gene") or r.get("Gene symbol") or "").strip()
                s = (r.get("start") or r.get("Start") or "").strip()
                e = (r.get("end") or r.get("Stop") or r.get("stop") or "").strip()
                if not g or not s or not e or g in seen:
                    continue
                seen.add(g)
                rows.append((int(float(s)), int(float(e)), g))
    except Exception:
        return False
    if not rows:
        return False
    with open(out_csv, "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["Start", "Stop", "Gene symbol"])
        w.writerows(rows)
    return True


def maybe_build_rdp5_pdf(root: str, work_dir: str, name: str,
                         log=None) -> Optional[str]:
    """用管线自带的 RDP5_RBDP_Rgrapher 生成**按排版尺寸**的断点分布图。

    2026-10-12: `rdp5_bdp_plot.py` 的画布/字号原先是写死的（ggsave width=30/height=8、
    R 侧 fontSizeMultiplier=2），30 in 宽的图塞进通栏只剩 0.31 倍、字号掉到 ~4 pt。
    现该脚本已暴露 `--width/--height/--font-scale/--gene-font-size`，这里按 16×6 生成。
    """
    def _log(m, lvl="emit"):
        if log is not None:
            getattr(log, lvl, log.emit)(m)
        else:
            print(m)

    inp = _find_rdp5_inputs(root)
    if not inp:
        _log("  [overview] 跳过面板 H (rdp5): 缺 breakpoint/positions/genes CSV")
        return None
    bp, pos, genes = inp
    here = os.path.dirname(os.path.abspath(__file__))
    tool = os.path.join(here, "rdp5_bdp_plot.py")
    if not os.path.exists(tool):
        _log("  [overview] 跳过面板 H (rdp5): 未找到 utils/rdp5_bdp_plot.py")
        return None
    os.makedirs(work_dir, exist_ok=True)
    orf = os.path.join(work_dir, "orfcoords.csv")
    if not _orfcoords_from_genes(genes, orf):
        _log(f"  [overview] 跳过面板 H (rdp5): 无法从 {os.path.basename(genes)} 解析基因坐标")
        return None
    out = os.path.join(work_dir, "rdp5_breakpoint.pdf")
    cmd = [sys.executable, tool, "-b", bp, "-g", orf, "-p", pos,
           "-n", name, "--width", str(RDP5_W), "--height", str(RDP5_H),
           "--font-scale", str(RDP5_FONT), "--gene-font-size", str(RDP5_GENE_FONT),
           "-o", out]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except Exception as e:
        _log(f"  [overview] 跳过面板 H (rdp5): 调用失败 {e}", "warning")
        return None
    if os.path.exists(out) and os.path.getsize(out) > 0:
        _log(f"  [overview] 面板 H (rdp5) 已生成: {out} ({os.path.getsize(out)} B)")
        return out
    _log(f"  [overview] 跳过面板 H (rdp5): 生成失败 rc={r.returncode} "
         f"{(r.stderr or '')[-300:]}", "warning")
    return None


# ═══════════════════════════════════════════════════════════════════
# 主入口
# ═══════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════
# 题注量宽 / 折行（第 7 轮新增）
# ═══════════════════════════════════════════════════════════════════
def _text_w_in(s: str, fs: float = CAP_FS) -> float:
    """返回字符串在字号 `fs` 下的**真实渲染宽度**（英寸）。

    用 `TextPath` 走 FreeType 的字形轮廓，与 `fig.text` 用同一套字体度量；
    比"字符数 × 平均宽"可靠（本图里 `W`/`i`/`×`/上标 差异极大）。
    空串返回 0；异常时退回保守估算，绝不返回 0 而让折行失效。
    """
    if not s:
        return 0.0
    try:
        from matplotlib.textpath import TextPath
        from matplotlib.font_manager import FontProperties
        tp = TextPath((0, 0), s, size=fs, prop=FontProperties())
        return float(tp.get_extents().width) / 72.0
    except Exception:
        return 0.62 * fs / 72.0 * len(s)   # 保守：0.62 em/字符


def _wrap_caption(s: str, max_w_in: float, fs: float = CAP_FS,
                  tol: float = CAP_WRAP_TOL):
    """把题注折成若干行，每行宽 ≤ `max_w_in * tol`。返回行列表。

    · `\\n` 是**硬换行**（作者可控），硬换行后仍超宽的段会继续自动折行；
    · ⚠ 单个词超过栏宽时**不硬切词**（会切出 `incompatib/le` 这种），
      而是让该行单独超宽 —— 那种情况说明题注本身要改写，不是折行能救的。
    """
    lim = max_w_in * tol
    lines = []
    for para in (s or "").split("\n"):
        words = para.split()
        if not words:
            continue
        cur = ""
        for w in words:
            cand = w if not cur else cur + " " + w
            if _text_w_in(cand, fs) <= lim or not cur:
                cur = cand
            else:
                lines.append(cur)
                cur = w
        if cur:
            lines.append(cur)
    return lines or [""]


def rasterize(path: str, dpi: int = DPI) -> Image.Image:
    if path.lower().endswith((".png", ".jpg", ".jpeg")):
        return Image.open(path).convert("RGB")
    import fitz
    doc = fitz.open(path)
    pm = doc[0].get_pixmap(matrix=fitz.Matrix(dpi / 72.0, dpi / 72.0), alpha=False)
    im = Image.frombytes("RGB", (pm.width, pm.height), pm.samples)
    doc.close()
    return im


def build_pipeline_overview(root: str, out_dir: str, *, title: str = "",
                            rdp5_pdf: Optional[str] = None,
                            prefix: str = "Fig", stem: str = "",
                            dpi: int = DPI, log=None) -> Dict[str, Optional[str]]:
    # ⚠ 本模块**不**调用 sci_palette.apply_sci_style（Fig3 有自己的一套字号/网格），
    #   因此必须自己把全局 rcParams 钉死，否则输出取决于同一进程里是否先跑过
    #   别的模块（它们会调用 apply_sci_style 改全局 rcParams）。实测同一份代码：
    #     单独跑        = 150 dpi / Type3 字体 / font.size=10 / 2.37 MB
    #     跟在别的任务后 = 300 dpi / Type42 字体 / font.size=8  / 5.97 MB
    #   而且**像素级也不一致** —— 对审稿人来说这就是"同一份数据两次跑出来的图不一样"。
    #   2026-10-12 二次修正：原先只手抄了字号/字体/颜色，仍残留差异
    #   （实测 mean|Δ|=0.19、0.46 % 像素差 >8），漏掉的是
    #   `axes.spines.top|right`（= 每个面板细边框的显隐：apply_sci_style 关掉、
    #   matplotlib 默认打开）、`axes.prop_cycle`、`xtick.major.size`、
    #   `lines.linewidth`、`patch.linewidth`、`legend.*`、`grid.*`。
    #   现改为直接复用**唯一真源** `sci_palette.sci_rcparams()` —— 谁都不许再手抄
    #   一份 rcParams（手抄必然漏键，漏键必然顺序依赖）。
    #   ⚠ 这里显式 `despine_top_right=False`：Fig3 每个面板都带一圈细灰边框
    #   （下面 `for s in ax.spines.values()` 会给四条边上色），去掉上/右会破坏观感。
    matplotlib.rcParams.update(sci_rcparams(base_size=8, grid=False,
                                            despine_top_right=False))
    """生成 Fig3 管道各阶段总览。返回 {'png':..,'pdf':..}（失败项为 None）。"""
    def _log(m, lvl="emit"):
        if log is not None:
            getattr(log, lvl, log.emit)(m)
        else:
            print(m)

    root = os.path.abspath(os.path.expanduser(root))
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    if not stem:
        stem = f"{prefix}3_pipeline_overview"

    # rdp5 面板：先看显式传入 / 树内已有的排版版，否则现生成
    if rdp5_pdf is None:
        import glob
        pre = []
        for rd in _recomb_rdp5_dirs(root):
            pre += glob.glob(os.path.join(
                rd, "breakpoint_distribution_publication.pdf"))
            pre += glob.glob(os.path.join(rd, "*_publication.pdf"))
        rdp5_pdf = sorted(pre)[0] if pre else None
    if rdp5_pdf is None:
        rdp5_pdf = maybe_build_rdp5_pdf(root, os.path.join(out_dir, "_fig3_src"),
                                        title or "Virus", log=log)

    info = {}
    for letter, gid, en, rel, kind in PANELS:
        if rel == "@RDP5@":
            p = rdp5_pdf
        else:
            p = resolve_panel_path(root, rel)
        if not p or not os.path.exists(p):
            _log(f"  [overview] 跳过面板 {letter}: {rel}")
            continue
        try:
            im = rasterize(p, dpi)
        except Exception as e:
            _log(f"  [overview] 跳过面板 {letter}: 读图失败 {e}", "warning")
            continue
        info[letter] = dict(en=en, gid=gid, path=p, img=im,
                            aspect=im.height / im.width)

    if len(info) < 3:
        _log("  [overview] 可用面板不足 3 个 → 不出图", "warning")
        return {"png": None, "pdf": None}

    layout = []
    for gid, letters in ROWS:
        letters = [x for x in letters if x in info]
        if not letters:
            continue
        n = len(letters)
        w_panel = (WU - (n - 1) * GAP_X) / n
        h_row = max(w_panel * info[x]["aspect"] for x in letters)
        # 第 7 轮：题注按**栏宽**折行（见 CAP_* 常量块的说明）。
        #   avail = w_panel（左列起点 M_L + WU - w_panel = W_FIG - M_R - w_panel 恒成立）
        caps = [_wrap_caption(info[x]["en"], w_panel) for x in letters]
        n_cap = max(len(c) for c in caps)
        if n_cap > 1:
            _log("  [overview] 题注折行(栏宽 %.2f in): %s" % (
                w_panel, ", ".join("%s->%d行" % (ltr, len(c))
                                   for ltr, c in zip(letters, caps) if len(c) > 1)))
        layout.append((gid, letters, w_panel, h_row, caps, n_cap))

    H = M_T + M_B
    prev = None
    for gid, letters, w_panel, h_row, caps, n_cap in layout:
        if gid != prev:
            H += H_HEADER + (GAP_GROUP if prev is not None else 0.0)
            prev = gid
        H += h_row + H_CAP + (n_cap - 1) * CAP_LINE_H

    _log(f"  [overview] 设计尺寸 {W_FIG:.2f} x {H:.2f} in "
         f"({len(layout)} 行 / {len(info)} 面板)")

    fig = plt.figure(figsize=(W_FIG, H), dpi=150)
    fig.patch.set_facecolor("white")
    fx = lambda v: v / W_FIG   # noqa: E731
    fy = lambda v: v / H        # noqa: E731

    _t = title or "Virus"
    # 2026-10-12: 全英文（原为主标题 + 组头 + 题注均为中英混排）
    fig.text(fx(M_L), fy(H - 0.42),
             f"{_t} — phylogenetic analysis pipeline",
             fontsize=13, fontweight="bold", ha="left", va="center")
    fig.text(fx(M_L), fy(H - 0.72),
             "Representative output of every analysis stage, "
             "from sampling metadata to phylodynamics and phylogeography",
             fontsize=9.5, color="#444444", ha="left", va="center")

    y_top = H - M_T
    prev = None
    _gidx = 0
    for gid, letters, w_panel, h_row, caps, n_cap in layout:
        if gid != prev:
            if prev is not None:
                y_top -= GAP_GROUP
            e = next((e for g, e in GROUPS if g == gid), "")
            _gidx += 1
            fig.text(fx(M_L), fy(y_top - 0.16), f"{_gidx}.  {e}",
                     fontsize=10.5, fontweight="bold", ha="left", va="center")
            fig.add_artist(plt.Line2D([fx(M_L), fx(W_FIG - M_R)],
                                      [fy(y_top - 0.30), fy(y_top - 0.30)],
                                      color="#999999", linewidth=0.7,
                                      transform=fig.transFigure))
            y_top -= H_HEADER
            prev = gid

        y_pt = y_top - h_row
        for i, letter in enumerate(letters):
            d = info[letter]
            x0 = M_L + i * (w_panel + GAP_X)
            h = w_panel * d["aspect"]
            ax = fig.add_axes([fx(x0), fy(y_pt), fx(w_panel), fy(h)])
            ax.imshow(d["img"], interpolation="antialiased")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_edgecolor("#bbbbbb"); s.set_linewidth(0.6)
            ax.set_anchor("N")
            ax.text(0.0, 1.005, letter, transform=ax.transAxes,
                    fontsize=11, fontweight="bold", ha="left", va="bottom")
            # 第 7 轮：折行后的多行题注；行距用 `linespacing` 与 CAP_LINE_H 对齐
            fig.text(fx(x0), fy(y_pt - CAP_TOP_OFF), "\n".join(caps[i]),
                     fontsize=CAP_FS, ha="left", va="top", color="#222222",
                     linespacing=1.2)
        y_top = y_pt - H_CAP - (n_cap - 1) * CAP_LINE_H

    out_png = os.path.join(out_dir, stem + ".png")
    out_pdf = os.path.join(out_dir, stem + ".pdf")
    fig.savefig(out_png, dpi=dpi, facecolor="white",
                bbox_inches="tight", pad_inches=0.12)
    # ⚠ PDF 必须显式给 dpi：不传就吃全局 rcParams['savefig.dpi']，而该值会被
    #   同一进程里先跑的其他模块（如 mcmc_diag_plot → sci_palette.apply_sci_style）
    #   改成 300。于是"图里内嵌位图的分辨率"取决于**任务执行顺序** ——
    #   同代码单独跑 = 150 dpi / 2.37 MB，跟在别的任务后面跑 = 300 dpi / 5.85 MB。
    #   显式钉死后输出与进程状态无关（审稿人视角的"可复现"要求）。
    fig.savefig(out_pdf, dpi=dpi, facecolor="white",
                bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    _log(f"  [overview] -> {out_pdf} ({os.path.getsize(out_pdf)} B)")
    _log(f"  [overview] -> {out_png} ({os.path.getsize(out_png)} B)")
    return {"png": out_png, "pdf": out_pdf}


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Pipeline stage overview figure (Fig3)")
    ap.add_argument("--root", required=True, help="data_dir (含各阶段目录)")
    ap.add_argument("--out", required=True, help="输出目录")
    ap.add_argument("--title", default="", help="图标题里的病毒名")
    ap.add_argument("--rdp5-pdf", default=None, help="现成的 rdp5 断点图 PDF")
    ap.add_argument("--prefix", default="Fig")
    ap.add_argument("--stem", default="")
    ap.add_argument("--dpi", type=int, default=DPI)
    a = ap.parse_args()
    r = build_pipeline_overview(a.root, a.out, title=a.title, rdp5_pdf=a.rdp5_pdf,
                                prefix=a.prefix, stem=a.stem, dpi=a.dpi)
    if not r.get("pdf"):
        sys.exit(1)


if __name__ == "__main__":
    main()
