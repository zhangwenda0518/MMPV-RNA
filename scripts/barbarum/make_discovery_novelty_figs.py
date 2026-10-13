#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_discovery_novelty_figs.py — 枸杞病毒组「新/已知 · 基因型 · 科属种」通路图

背景
----
用户 2026-10-13 提供两张参考图（`gabrielvpina/dataViz` 风格），要求**基于 09b 的表格**
补绘同类图，覆盖维度：
    新病毒 / 已知病毒 · DNA / RNA / 类病毒 · 科 / 属 / 种 ·
    对属的新种 · 属→goji 的命名

参考图形态
----------
  ① **bipartite 通路图**：左侧两个粗色块（`Novel species` 红 / `Known species` 蓝）
     → 中间科名色块 → 右侧物种名，之间用**灰带**连接（带宽 ∝ 条目数）。
  ② **identity × coverage 散点 + 上下/左右边缘箱线**。

与既有的关系
------------
`make_discovery_sci_figs.py` 的面板 **J** 已有一版 identity×coverage 散点（按 `verdict`
着色）。本脚本的 panel **c** 是它的**改造版**：改为按**新/已知 + 基因型**着色，并叠加
「种界 95 % / 属界 75 %」参考分区 —— 把参考图①的分类信息压进参考图②。

数据口径（**全部来自 09b 的表，脚本不重算统计量**）
--------------------------------------------------
  panel a  `09b_Analysis_Verify/keep_summary.tsv`（101 行）
           `novelty_9a` 三值分档 → Novel / Known
           `Family` / `Genus` / `Species` + `len_ratio`（对属平均长的比）
  panel b  `09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv`（291 行）
           `category` 中文分档 → Novel / Known / Contaminant
           `tax_family` / `tax_genus` / `tax_species` + `nt_pident` / `nt_qcovs`
  panel c  同 panel b 的 `nt_pident` × `nt_qcovs`

⚠ **两张表口径不同，必须在图上分别标清**（用户选「两张表都画」）：
    · keep 表 = 09b 顶层 keep 集（101 条），novelty 三档来自 `novelty_9a`；
    · rescue 表 = `virus_validation` 的 291 条候选，category 来自 `rescue` 复检。
  两者的「科/属」列名不同（`Family` vs `tax_family`）→ 用各自的列，勿交叉。

⚠ **`len_ratio`（对属平均长的比）** 只有 keep 表有（`rescue` 表是 `length`/`genus_avg_len`
  两列，需自行相除）。本项目里 keep 表 83/101 行有值。

⚠ 图内**全英文**（与 discovery 侧其它静态图一致）；PDF 文本层 CJK 必须 = 0。

用法
----
    python scripts/barbarum/make_discovery_novelty_figs.py \
      --root  <RNA-Lycium_barbarum_out> \
      --outdir <...>/10_Reports/S10_SCI_Discovery \
      --virus-name "Lycium barbarum" [--panels-only] [--composite-only] [--dpi 300]

产出
----
    <outdir>/Fig_GojiVirome_Novelty.{png,pdf}         三联图 a/b/c
    <outdir>/panels_novelty/panel_a|b|c.{png,pdf}      单面板
    <outdir>/_build/novelty_manifest.json              面板数值来源
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import traceback
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np

# ── 统一配色（唯一真源）──
# ⚠ **不要**写成 `from utils.sci_palette import ...`：服务器 site-packages 里有同名
# 常规包会静默盖掉 → 必须按**显式文件路径**加载（2026-07-xx 起沿用至今）。
import importlib.util as _ilu

_HERE = Path(__file__).resolve()
_REPO = _HERE.parents[2]


def _load_palette():
    for rel in ("virome_discovery_pipeline/utils/sci_palette.py",
                "virome_phylo_pipeline/utils/sci_palette.py"):
        cand = _REPO / rel
        if cand.is_file():
            spec = _ilu.spec_from_file_location("_sci_palette", cand)
            mod = _ilu.module_from_spec(spec)
            sys.modules["_sci_palette"] = mod
            spec.loader.exec_module(mod)
            return mod
    raise ImportError(f"sci_palette.py not found under {_REPO}")


_P = _load_palette()
apply_sci_style = _P.apply_sci_style
SEMANTIC = _P.SEMANTIC
despine = _P.despine
sci_rcparams = _P.sci_rcparams

DEFAULT_DPI = 300

# ══════════════════════════════════════════════════════════════════
# 配色：novel / known / contaminant 三分（用户要的「新病毒 / 已知病毒」）
# ══════════════════════════════════════════════════════════════════
# ⚠ 这套色是**本图专用**的语义色（不属于 sci_palette 的既有键）→ 在此显式定义并
#   写进 manifest，保证「同义同色」可复核。全部取自 OKABE_ITO，不引入游离 hex。
C_NOVEL = "#D55E00"        # vermillion —— 新病毒 / 新种候选（= SEMANTIC 里"新"系）
C_KNOWN = "#0072B2"        # blue —— 已知病毒
C_CONTAM = "#999999"       # grey —— 污染 / 未定（不承载结论）

NOVEL_COLORS = {
    "Novel": C_NOVEL,
    "Known": C_KNOWN,
    "Contaminant": C_CONTAM,
}

C_BAND = "#4D4D4D"         # 灰带（参考图①的连接带）
C_BAND_ALPHA = 0.30
C_NODE_EDGE = "#333333"

# 档的固定展示顺序（Novel → Known → Contaminant）；所有排序/取档都锚定它
TIER_ORDER = ["Novel", "Known", "Contaminant"]


# ══════════════════════════════════════════════════════════════════
# I/O 辅助
# ══════════════════════════════════════════════════════════════════
def read_tsv(path: Path):
    if not path.is_file():
        return None
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _num(s, default=None):
    if s is None:
        return default
    t = str(s).strip().strip('"')
    if t in ("", "NA", "N/A", "-", "nan", "None"):
        return default
    t = t.replace(",", "").rstrip("%")
    try:
        return float(t)
    except ValueError:
        return default


def _txt(s):
    """列值 → 干净的字符串；`NA`/`-`/空 → ""（便于统一当"无归属"处理）。"""
    t = (s or "").strip()
    return "" if t in ("NA", "N/A", "-", "nan", "None") else t


def _find_09b(root: Path) -> Path | None:
    for name in ("09b_Analysis_Verify", "09b_analysis_verify"):
        p = root / name
        if p.is_dir():
            return p
    return None


def _find_reports(root: Path) -> Path:
    for name in ("10_Reports", "09_Reports", "09_report", "reports"):
        p = root / name
        if p.is_dir():
            return p
    return root


# ══════════════════════════════════════════════════════════════════
# 数据装配
# ══════════════════════════════════════════════════════════════════
class Data:
    def __init__(self, root: Path, virus_name: str):
        self.root = root
        self.virus_name = virus_name
        self.rep = _find_reports(root)
        self.manifest = {}
        self._load()

    def _load(self):
        b = _find_09b(self.root)
        self.b09 = b
        self.keep = self._pick("keep_summary.tsv", b)
        self.rescue = self._pick("virus_validation/rescue_evidence_scored.tsv", b)

    def _pick(self, rel: str, b: Path | None) -> list | None:
        """按 09b → root → reports 顺序找表（本地镜像与服务器路径不同）。"""
        cands = []
        if b is not None:
            cands.append(b / rel)
        cands.append(self.root / rel)
        cands.append(self.rep / rel)
        cands.append(self.rep / Path(rel).name)
        for c in cands:
            rows = read_tsv(c)
            if rows:
                self.manifest.setdefault("_sources", {})[Path(rel).name] = str(c)
                return rows
        return None

    # ── 派生：三档标签 ─────────────────────────────────────────
    @staticmethod
    def _tri_keep(raw: str) -> str:
        """`novelty_9a` → Novel / Known（keep 表没有污染档，全为植物 keep 集）。"""
        t = (raw or "").strip()
        if t.startswith("已知"):
            return "Known"
        if t.startswith("新"):
            return "Novel"
        return "Novel"          # 兜底：keep 集本体即"候选"

    @staticmethod
    def _tri_rescue(raw: str) -> str:
        """`category` → Novel / Known / Contaminant。"""
        t = (raw or "").strip()
        if t in ("已知植物病毒",):
            return "Known"
        if t in ("新种候选", "远缘病毒候选", "远缘(未定)"):
            return "Novel"
        if "污染" in t or t == "未定":
            return "Contaminant"
        return "Novel"

    # ── 派生：新种判据（用户 2026-10-13 选定）───────────────────
    #   「len_ratio > 1 或 (nt_pident 75–95 % + 高覆盖)」
    #   · 长于属平均 = 基因组规模异常 → 提示新种（ICTV 物种界定含基因组长度）
    #   · nt identity 落在 75–95 % = 高于属界（~75 %）但低于种界（~95 %）
    #     + qcov ≥ 70 % = 比对覆盖足够（低覆盖下的 identity 不可信）
    @staticmethod
    def is_new_species(len_ratio=None, pident=None, qcov=None) -> bool:
        if len_ratio is not None and len_ratio > 1.0:
            return True
        if pident is not None and qcov is not None:
            return (75.0 <= pident <= 95.0) and (qcov >= 70.0)
        return False

    # ── 派生：属→goji 命名 ─────────────────────────────────────
    def genus_goji(self, genus: str, species: str) -> str | None:
        """判定该物种名是否为「宿主关联 / 枸杞相关」命名。

        用户要的「属→goji 的命名」= 属名 → 该属内以**宿主/地点**命名的物种
        （xx-like virus / host-associated / 枸杞学名 `lycii`）。返回命中的
        关键词，未命中返回 None。
        """
        s = (species or "").lower()
        g = (genus or "").lower()
        if not s:
            return None
        keys = []
        if "lycii" in s or "goji" in s:
            keys.append("lycii")            # 枸杞学名 Lycium barbarum → lycii
        if "-like virus" in s or "like virus" in s:
            keys.append("virus-like")
        if "associated" in s:
            keys.append("associated")
        if "sp." in s:
            keys.append("sp.")
        return "/".join(keys) if keys else None

    # ── panel a/b 的通用结构 ───────────────────────────────────
    def flow(self, which: str, min_family=2, top_species=14, top_families=14):
        """返回 bipartite 通路图的三列数据。

        `which` = "keep" | "rescue"。
        ⚠ `top_families` 是**版式硬约束**而非取数偏好：族数 × 各族物种行数决定
          所需高度，超出画布高度就会把末尾族挤出（实测）。取不到就如实截断并在
          标题注记里写明（`… showing top N families`）。
        返回 `(tiers, families, species, edges, src)`。
        """
        rows = self.keep if which == "keep" else self.rescue
        if not rows:
            return None
        if which == "keep":
            fam_k, gen_k, sp_k, tri = "Family", "Genus", "Species", self._tri_keep
            n_k = "novelty_9a"
            src = "09b/keep_summary.tsv"
        else:
            fam_k, gen_k, sp_k, tri = "tax_family", "tax_genus", "tax_species", self._tri_rescue
            n_k = "category"
            src = "09b/virus_validation/rescue_evidence_scored.tsv"

        recs = []
        for r in rows:
            t = tri(r.get(n_k) or "")
            f = _txt(r.get(fam_k))
            g = _txt(r.get(gen_k))
            s = _txt(r.get(sp_k))
            if not f and not s:
                continue
            recs.append((t, f or "(unassigned)", g, s or "(unassigned)", r))
        if len(recs) < 5:
            return None

        fam_cnt = defaultdict(int)
        for t, f, g, s, r in recs:
            fam_cnt[f] += 1
        fams = [(f, n) for f, n in fam_cnt.items() if n >= min_family]
        if len(fams) < 3:
            return None
        fams.sort(key=lambda x: -x[1])
        n_fam_total = len(fams)
        fams = fams[:top_families]          # 版式硬约束（见 docstring）
        keep_f = {f for f, _ in fams}

        sp_cnt = defaultdict(int)
        sp_genus = {}
        for t, f, g, s, r in recs:
            if f in keep_f:
                sp_cnt[(f, s)] += 1
                if g and (f, s) not in sp_genus:
                    sp_genus[(f, s)] = g
        sps = sorted(sp_cnt.items(), key=lambda x: (-x[1], x[0][0], x[0][1]))[:top_species]

        # ⚠ 一个科内可能**同时**有 Novel / Known / Contaminant 记录（如 Betaflexiviridae
        #   既有已知 Carlavirus 又有新种候选）→ 不能「后写覆盖」（会得到随机色，
        #   且可能出现某档在 families 里不存在却在 recs 里的 KeyError）。
        #   取**出现次数最多**的档作为该科的展示档，并把计数记进 manifest（可复核）。
        fam_tri = defaultdict(lambda: defaultdict(int))
        for t, f, g, s, r in recs:
            if f in keep_f:
                fam_tri[f][t] += 1
        tier_of = {}
        for f in keep_f:
            d = fam_tri[f]
            tier_of[f] = max(d.items(), key=lambda kv: (kv[1], -TIER_ORDER.index(kv[0])))[0]
        # ⚠ 档的展示顺序固定 Novel → Known → Contaminant（跳号也稳定）
        tiers = [x for x in TIER_ORDER if x in set(tier_of.values())]

        families = [(f, n, tier_of.get(f, "Novel")) for f, n in fams]
        species = [(f"{s}  ({sp_genus[(f, s)]})" if sp_genus.get((f, s)) else s,
                    n, f, tier_of.get(f, "Novel"))
                   for (f, s), n in sps]
        edges = []
        tcnt = defaultdict(int)
        for t, f, g, s, r in recs:
            if f in keep_f:
                tcnt[(t, f)] += 1
        for (t, f), n in tcnt.items():
            edges.append((t, f, n))
        for (f, s), n in sps:
            lab = (f"{s}  ({sp_genus[(f, s)]})" if sp_genus.get((f, s)) else s)
            edges.append((f, lab, n))
        # ── 属 → goji / 宿主关联命名的统计（用户要的「属→goji 的命名」维度）──
        #   ⚠ 必须在 `flow()` 内算（`recs` 是局部变量）；随返回值交给 `_bipartite`
        #     写进 manifest，保证"命名判定"可复核、不靠图面文字自证。
        seen_n = set()
        naming = []
        for t, f, g, s, r in recs:
            if not g or (g, s) in seen_n:
                continue
            key = self.genus_goji(g, s)
            if key:
                seen_n.add((g, s))
                naming.append({"genus": g, "species": s, "family": f,
                               "class": t, "naming_key": key})
        naming.sort(key=lambda d: (d["genus"].lower(), d["species"].lower()))
        return tiers, families, species, edges, src, n_fam_total, naming

    def scatter(self, which="rescue", min_pts=8):
        """identity × coverage 散点数据（panel c）。

        返回 `(pts, xlab, ylab)`，`pts = [(ident, cov, tier, family, is_new), ...]`。
        ⚠ 只取**两坐标都有效**（> 0）的行：`nt_pident = 0` 表示无成对命中，
          画进去会在 x=0 堆一条假竖线（与面板 J 同一条纪律）。
        """
        rows = self.keep if which == "keep" else self.rescue
        if not rows:
            return None
        if which == "keep":
            x_k, y_k, fam_k, lenr = "nt_pident", "nt_qcovs", "Family", "len_ratio"
            tri = self._tri_keep
            n_k = "novelty_9a"
        else:
            x_k, y_k, fam_k, lenr = "nt_pident", "nt_qcovs", "tax_family", None
            tri = self._tri_rescue
            n_k = "category"

        pts = []
        for r in rows:
            x = _num(r.get(x_k))
            y = _num(r.get(y_k))
            if x is None or y is None or x <= 0 or y <= 0:
                continue
            t = tri(r.get(n_k) or "")
            f = _txt(r.get(fam_k)) or "(unassigned)"
            lr = _num(r.get(lenr)) if lenr else None
            new = self.is_new_species(lr, x, y)
            pts.append((x, y, t, f, new))
        if len(pts) < min_pts:
            return None
        return pts, "Nucleotide identity (%)", "Coverage (%)", which


# ══════════════════════════════════════════════════════════════════
# 通用绘制
# ══════════════════════════════════════════════════════════════════
def _empty(ax, msg="(no data)"):
    ax.axis("off")
    ax.text(0.5, 0.5, msg, ha="center", va="center", color="#999999",
            fontsize=7, transform=ax.transAxes)


def _band(ax, x0, y0, x1, y1, lw, color=C_BAND, alpha=C_BAND_ALPHA, z=1):
    """灰带 = 「参考图①」的连接带（贝塞尔曲线 + 圆头 + 半透明）。

    ⚠ 用 `FancyArrowPatch` 的 `connectionstyle="arc3,rad=..."` 而不是直线：
      参考图的带子是**S 形**（左端水平出发、右端水平进入），直线会显得生硬、
      且在密集处把节点"穿"成一团。
    """
    p = FancyArrowPatch((x0, y0), (x1, y1),
                        connectionstyle="arc3,rad=0.22",
                        arrowstyle="-", lw=lw, color=color, alpha=alpha,
                        zorder=z, capstyle="round")
    ax.add_patch(p)


def _curve_band(ax, x0, y0, x1, y1, lw, color=C_BAND, alpha=C_BAND_ALPHA, z=1):
    """S 形带：控制点取两端中点，做成水平出发/进入的样条。"""
    from matplotlib.path import Path as MPath
    from matplotlib.patches import PathPatch
    xm = (x0 + x1) / 2.0
    verts = [(x0, y0), (xm, y0), (xm, y1), (x1, y1)]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
    pp = PathPatch(MPath(verts, codes), facecolor="none", edgecolor=color,
                   lw=lw, alpha=alpha, zorder=z, capstyle="round")
    ax.add_patch(pp)


def _text_w_in(s: str, fontsize: float, weight: str = "normal") -> float:
    """文本真实字宽（in）—— 用于**折行/判断右缘是否出界**（项目既有方法）。"""
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    fp = FontProperties(size=fontsize, weight=weight)
    try:
        tp = TextPath((0, 0), s, prop=fp)
        return tp.get_extents().width / 72.0
    except Exception:
        return len(s) * fontsize * 0.5 / 72.0


def _fit_label(s: str, max_in: float, fontsize: float) -> str:
    """把标签截断到给定的物理宽度内（保留尾部 "(Genus)" 若存在）。"""
    if _text_w_in(s, fontsize) <= max_in:
        return s
    tail = ""
    if s.endswith(")") and "(" in s:
        i = s.rfind("(")
        tail, s = s[i:], s[:i].rstrip()
    ell = "…"
    while s and _text_w_in(s + ell + (" " + tail if tail else ""), fontsize) > max_in:
        s = s[:-1]
    return s + ell + (" " + tail if tail else "")


def _pill(ax, x, y, w, h, text, color, fontsize, text_color="white",
          edge=C_NODE_EDGE, lw=0.5, z=4):
    """圆角色块 + 居中文字（参考图①的科名色块 / 物种名色块）。"""
    from matplotlib.patches import FancyBboxPatch
    pad = h * 0.30
    box = FancyBboxPatch((x, y - h / 2 + pad / 2), w, h - pad,
                         boxstyle=f"round,pad=0,rounding_size={h * 0.30}",
                         facecolor=color, edgecolor=edge, lw=lw, zorder=z,
                         mutation_aspect=1.0)
    ax.add_patch(box)
    ax.text(x + w / 2.0, y, text, ha="center", va="center", fontsize=fontsize,
            color=text_color, zorder=z + 1)


# ══════════════════════════════════════════════════════════════════
# panel a / b：Novel-Known → 科 → 种（bipartite 通路图）
# ══════════════════════════════════════════════════════════════════
def _bipartite(ax, D: Data, which: str, title: str, note: str, letter_tag: str,
               box_in=None, min_family=2, top_species=14, top_families=14):
    """bipartite 通路图（参考图①）。

    ⚠ **坐标系用「英寸」**（`xlim=(0,W)` / `ylim=(0,H)`，W/H 由 `box_in` 给定）：
      唯一的理由是本面板要**按文字物理宽度**决定色块宽度/截断 —— `_text_w_in()`
      返回英寸，若坐标系是 0–1 的轴分数，二者**不同量纲**，`_fit_label` 会把
      「Betaflexiviridae」截成「Be…」（实测），这是本项目版式事故的高频来源。
      用英寸坐标系后，`W = 标签宽 + 两侧内边距` 即为真值，不再需要猜系数。
    ⚠ `box_in=None` 时按轴的实际物理尺寸自适应（组合图里每个槽位大小不同）。
    """
    got = D.flow(which, min_family=min_family, top_species=top_species,
                 top_families=top_families)
    if not got:
        return _empty(ax)
    tiers, families, species, edges, src, n_fam_total, naming_uniq = got
    if box_in is None:
        try:
            pos = ax.get_position()
            fw, fh = ax.figure.get_size_inches()
            box_in = (pos.width * fw, pos.height * fh)
        except Exception:
            box_in = (7.3, 2.55)
    W, H = box_in
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    FS_TITLE = 7.4
    FS_FAM = 5.6
    FS_SP = 5.2
    FS_TIER = 5.6

    # ── x 布局（英寸）──
    X_TIER = 0.03
    W_TIER = 0.17
    X_FAM = 0.90
    X_SP = 3.55
    PAD = 0.055                       # 色块内左右留白
    fam_w = max(_text_w_in(f, FS_FAM, "bold") for f, _, _ in families) + 2 * PAD
    fam_w = min(fam_w, 1.35)
    # 物种标签最长者决定色块宽（但不超过剩余空间）
    sp_w = 0.0
    for lab, n, f, t in species:
        sp_w = max(sp_w, _text_w_in(lab, FS_SP))
    sp_w = min(sp_w + 2 * PAD, W - X_SP - 0.30)
    if sp_w < 1.2:
        sp_w = 1.2

    # ── y 布局（英寸）：按档分组铺科块 ──
    # 先按档分组重排（Novel 在上、Known 次之、Contaminant 最后）
    fams_sorted = sorted(families, key=lambda x: (TIER_ORDER.index(x[2]), -x[1]))
    # ⚠ 每个科的「行高」必须 ≥ 其下挂物种数 × 行距，**否则相邻科的物种行会叠**。
    #   且总高可能超画布 → 先按理想行距算，超了就**同时压缩科间距与行距**，
    #   但行距有下限（SP_MIN），避免压到看不清。
    sp_by_fam = defaultdict(list)
    for lab, n, f, t in species:
        sp_by_fam[f].append((lab, n, t))
    SP_STEP = 0.150                   # 种节点理想行距
    SP_MIN = 0.115                    # 行距下限
    FAM_GAP = 0.075
    HEAD = 0.10                       # 数据区顶部额外留白（in）
    TOP_PAD = 0.10                    # 数据区顶端小留白（标题/注记都在轴**外**，
                                      #   不占数据区高度 → 无需大留白）
    LEG_H = 0.02                      # 数据区底部小留白（图例已移到轴外，见下）
    nfam = len(fams_sorted)
    for _ in range(30):
        total = sum(max(0.19, len(sp_by_fam.get(f, [])) * SP_STEP)
                    for f, _, _ in fams_sorted) + FAM_GAP * (nfam - 1)
        avail = H - HEAD - TOP_PAD - LEG_H
        if total <= avail or SP_STEP <= SP_MIN:
            break
        # 超了就按缺口等比缩，但行距不破下限
        SP_STEP = max(SP_MIN, SP_STEP * avail / total)
        FAM_GAP = max(0.045, FAM_GAP * avail / total)
    # ⚠ 压缩后可能**仍**超（SP_STEP 已到下限）→ 把整体从顶部往下压到 LOW 之上，
    #   强制底部对齐，宁可族间距变小也不越界。
    total = sum(max(0.19, len(sp_by_fam.get(f, [])) * SP_STEP)
                for f, _, _ in fams_sorted) + FAM_GAP * (nfam - 1)
    avail = H - HEAD - TOP_PAD - LEG_H
    if total > avail:
        scale = max(0.35, avail / total)
        SP_STEP = max(SP_MIN * 0.7, SP_STEP * scale)
        FAM_GAP = max(0.030, FAM_GAP * scale)

    fam_slot = {}
    fam_y = {}
    y = H - TOP_PAD
    for f, n, t in fams_sorted:
        h = max(0.19, len(sp_by_fam.get(f, [])) * SP_STEP)
        fam_y[f] = y - h / 2.0
        fam_slot[f] = h
        y -= (h + FAM_GAP)
    # 档块：覆盖其成员科的整段 y 范围
    tier_span = {}
    for t in tiers:
        ys = [(fam_y[f] - fam_slot[f] / 2.0, fam_y[f] + fam_slot[f] / 2.0)
              for f, n, tt in fams_sorted if tt == t]
        if ys:
            tier_span[t] = (min(a for a, _ in ys), max(b for _, b in ys))

    # ── 种节点 y ──
    sp_y = {}
    for f, lst in sp_by_fam.items():
        if f not in fam_y:
            continue
        yc = fam_y[f]
        k = len(lst)
        top = yc + (k - 1) * SP_STEP / 2.0
        for i, (lab, n, t) in enumerate(lst):
            sp_y[(f, lab)] = top - i * SP_STEP

    # ── 画带（先带后块）──
    max_fam = max(n for _, n, _ in families)
    for t, f, n in edges:
        if f not in fam_y:
            continue
        lw = 0.55 + 3.2 * (n / max_fam)
        _curve_band(ax, X_TIER + W_TIER, sum(tier_span[t]) / 2.0, X_FAM, fam_y[f], lw)
    max_sp = max(n for _, n, _, _ in species) if species else 1
    for f, lst in sp_by_fam.items():
        for lab, n, t in lst:
            if f not in fam_y or (f, lab) not in sp_y:
                continue
            lw = 0.5 + 2.4 * (n / max_sp)
            _curve_band(ax, X_FAM + fam_w, fam_y[f], X_SP, sp_y[(f, lab)], lw)

    # ── 画节点 ──
    for t in tiers:
        if t not in tier_span:
            continue
        ylo, yhi = tier_span[t]
        yc = (ylo + yhi) / 2.0
        from matplotlib.patches import FancyBboxPatch
        ax.add_patch(FancyBboxPatch((X_TIER, ylo), W_TIER, yhi - ylo,
                                    boxstyle="round,pad=0.02,rounding_size=0.05",
                                    facecolor=NOVEL_COLORS[t], edgecolor=C_NODE_EDGE,
                                    lw=0.5, zorder=5))
        # ⚠ **不画竖排档名**：档块高度随族数变化，竖排文字要么被裁成半截、
        #   要么伸长到轴外与底部图例相撞（实测 "Contaminant" × 图例 "Novel"
        #   撞 2.8×1.0 px）。档的身份**只靠颜色 + 轴外图例**表达，更稳。
        _ = (ylo, yhi)                       # 保留 span 计算，不再落文字

    for f, n, t in fams_sorted:
        yc = fam_y[f]
        h = 0.175
        ax.add_patch(FancyBboxPatch((X_FAM, yc - h / 2), fam_w, h,
                                    boxstyle="round,pad=0.0,rounding_size=0.05",
                                    facecolor=NOVEL_COLORS[t], edgecolor=C_NODE_EDGE,
                                    lw=0.5, zorder=5))
        ax.text(X_FAM + PAD, yc, _fit_label(f, fam_w - 2 * PAD, FS_FAM),
                ha="left", va="center", fontsize=FS_FAM, color="white", zorder=6,
                fontweight="bold")
        ax.text(X_FAM + fam_w + 0.035, yc, f"n={n}", ha="left", va="center",
                fontsize=4.6, color="#777777", zorder=7)

    for lab, n, f, t in species:
        yc = sp_y.get((f, lab))
        if yc is None:
            continue
        h = 0.125
        ax.add_patch(FancyBboxPatch((X_SP, yc - h / 2), sp_w, h,
                                    boxstyle="round,pad=0.0,rounding_size=0.045",
                                    facecolor=NOVEL_COLORS[t], edgecolor=C_NODE_EDGE,
                                    lw=0.45, zorder=5))
        ax.text(X_SP + PAD, yc, _fit_label(lab, sp_w - 2 * PAD, FS_SP),
                ha="left", va="center", fontsize=FS_SP, color="white", zorder=6)

    # ⚠ **不在轴内/轴缘画图例**：档块的竖排名与轴下图例都会随族数/画布变化而与
    #   族色块相撞（实测 Bromoviridae × 图例 "Known" 撞 10×5 px）。档的颜色语义
    #   统一放进**右侧独立说明卡**（`Reading notes`，一个 `axis("off")` 的轴，
    #   与任何通路图元都不同轴 → 结构上不可能重叠）。

    ax.text(0.0, 1.085, letter_tag, transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="center", ha="left", color="#111111")
    ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", loc="left", x=0.030,
                 y=1.085)
    # ⚠ 注记走 axes 分数坐标（数据坐标是英寸 0..H，漏 transform 会落到图正中）。
    # ⚠ 且**横向从 0.045 起**：面板字母占 x∈[0, 0.030]，注记左移会让二者
    #   bbox 在竖直方向擦边（实测 9×3 px）→ 直接错开 x，与 y 无关地根治。
    # ⚠ 截断必须**如实写进图面**（`TOP N of M families`）：读者不能以为画的就是全部。
    if n_fam_total > len(families):
        note = note + f"  ⚠ Showing top {len(families)} of {n_fam_total} families."
    ax.text(0.045, 1.000, note, transform=ax.transAxes, fontsize=4.8,
            color="#444444", va="bottom", ha="left", linespacing=1.35)
    # ── 属 → goji / 宿主关联命名（已由 `flow()` 算好，写进 manifest）──
    D.manifest[letter_tag] = {
        "source": src,
        "which": which,
        "n_families_shown": len(families),
        "n_families_total": n_fam_total,
        "n_species": len(species),
        "tiers": tiers,
        "colors": dict(NOVEL_COLORS),
        "note": note,
        "host_linked_naming": {
            "rule": "species name contains `lycii`/`goji`, `-like virus`, "
                    "`associated`, or `sp.`",
            "n": len(naming_uniq),
            "items": naming_uniq,
        },
    }
    return None


def panel_a_keep(ax, D: Data):
    """a  Novel / Known 物种 → 科 → 物种（09b keep_summary.tsv）。"""
    note = ("Source: 09b/keep_summary.tsv (n=101 kept contigs); class from `novelty_9a` "
            "(known virus / novel candidate by BLAST protein or domain-only).\n"
            "Grey bands ∝ record count; family label = `Family`, species label = "
            "`Species` (Genus).")
    return _bipartite(ax, D, "keep",
                      "Novel vs known species by family (09b keep set)",
                      note, "a", min_family=2, top_species=3, top_families=10)


def panel_b_rescue(ax, D: Data):
    """b  Novel / Known / Contaminant → 科 → 物种（09b rescue 复检表）。"""
    note = ("Source: 09b/virus_validation/rescue_evidence_scored.tsv (n=291 candidates); "
            "class from `category` (known plant virus / new-species / distant candidate / "
            "contaminant).\nDifferent table & column names from panel a — same layout, "
            "independent record set.")
    return _bipartite(ax, D, "rescue",
                      "Novel vs known species by family (09b rescue set)",
                      note, "b", min_family=3, top_species=3, top_families=12)


def panel_c_scatter(ax, D: Data):
    """c  identity × coverage 散点：按新/已知着色 + 种界/属界分区（参考图② 改造版）。"""
    got = D.scatter("rescue")
    if not got:
        return _empty(ax)
    pts, xlab, ylab, which = got
    x = np.array([p[0] for p in pts], dtype=float)
    y = np.array([p[1] for p in pts], dtype=float)

    ax.axvspan(95, 101, color="#0072B2", alpha=0.045, zorder=0)
    ax.axvspan(0, 75, color="#999999", alpha=0.055, zorder=0)
    ax.axvspan(75, 95, color="#D55E00", alpha=0.055, zorder=0)

    groups = defaultdict(lambda: ([], []))
    new_x, new_y = [], []
    for xx, yy, t, f, isn in pts:
        groups[t][0].append(xx)
        groups[t][1].append(yy)
        if isn:
            new_x.append(xx)
            new_y.append(yy)

    for t in ("Novel", "Known", "Contaminant"):
        gx, gy = groups.get(t, ([], []))
        if not gx:
            continue
        ax.scatter(gx, gy, s=14, alpha=0.70, color=NOVEL_COLORS[t],
                   edgecolors="#333333", linewidths=0.25,
                   label=f"{t} (n={len(gx)})", zorder=3)
    if new_x:
        ax.scatter(new_x, new_y, s=46, facecolors="none", edgecolors="#111111",
                   linewidths=0.7, marker="o", zorder=4,
                   label=f"new-species call (n={len(new_x)})")

    ax.set_xlabel(xlab, fontsize=6.2)
    ax.set_ylabel(ylab, fontsize=6.2)
    ax.set_xlim(0, 104)
    ax.set_ylim(0, 104)
    ax.tick_params(labelsize=5.2)
    for xv, lab in ((75, "75% genus"), (95, "95% species")):
        ax.axvline(xv, color="#666666", lw=0.55, ls=":", zorder=1)
        ax.text(xv, 101.5, lab, fontsize=4.4, color="#444444",
                rotation=90, va="top", ha="right")

    ax.set_title("Identity vs coverage by class (n = {:,})".format(len(pts)),
                 fontsize=7.2, fontweight="bold", loc="left", x=0.030)
    ax.text(0.0, 1.055, "c", transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="center", ha="left", color="#111111")
    ax.legend(loc="lower left", fontsize=4.6, frameon=True, borderpad=0.25,
              labelspacing=0.24, handletextpad=0.35, framealpha=0.92)
    despine(ax)
    D.manifest["c"] = {
        "source": "09b/virus_validation/rescue_evidence_scored.tsv "
                  "(`nt_pident` × `nt_qcovs`)",
        "n_points": len(pts),
        "n_new_species_call": len(new_x),
        "note": "Shaded zones: 0-75% below genus threshold; 75-95% between genus and "
                "species thresholds (candidate new species); 95-100% above species "
                "threshold (known species).\n"
                "Hollow rings = new-species call by rule (len_ratio>1 OR "
                "nt_identity 75-95% with qcov>=70%).",
    }


# ══════════════════════════════════════════════════════════════════
# 组装
# ══════════════════════════════════════════════════════════════════
PANELS = [
    ("a", panel_a_keep),
    ("b", panel_b_rescue),
    ("c", panel_c_scatter),
]

WIDE = {"a", "b"}       # 通路图要宽；散点近方形


def save(fig, outdir: Path, stem: str, dpi: int):
    outdir.mkdir(parents=True, exist_ok=True)
    png, pdf = outdir / f"{stem}.png", outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", pad_inches=0.05,
                facecolor="white")
    fig.savefig(pdf, dpi=dpi, bbox_inches="tight", pad_inches=0.05,
                facecolor="white")
    return png, pdf


def _letter(ax, letter):
    """兜底：仅在面板函数**抛异常**（`_empty` 占位）时补画字母。

    ⚠ 正常路径下字母已写进各面板的标题前缀（见 `_bipartite` / `panel_c_scatter`），
      这里只在失败占位图上补一个，且贴在轴**内**左上角，避免 TEXT-OFF-CANVAS。
    """
    ax.text(0.0, 1.02, letter, transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="bottom", ha="left", color="#111111")


def build_composite(D: Data, outdir: Path, dpi: int):
    # 版式：**3 行 × 2 列**，a/b 各跨 2 列（通路图宽高比 ≈ 2.8），c 占左下、
    #   右下放口径说明卡。
    #   ⚠ 通路图的高度是**硬约束**：族数 × (每族物种行数) 是真实需求，盒子矮了
    #     必然把末尾几个族连同其物种挤出画布（实测：2.4 in 高时 Kitaviridae 掉出）。
    #     故行高比按「a/b 各 3.3 in、c 3.0 in」定，画布总高 ≈ 11.0 in。
    fig = plt.figure(figsize=(8.6, 12.0))
    gs = GridSpec(3, 2, figure=fig, hspace=0.50, wspace=0.30,
                  left=0.055, right=0.985, top=0.935, bottom=0.035,
                  height_ratios=[1.12, 1.20, 0.95])
    made = []
    for i, (letter, fn) in enumerate(PANELS):
        if letter in WIDE:
            ax = fig.add_subplot(gs[i, 0:2])
        else:
            ax = fig.add_subplot(gs[i, 0:1])
        try:
            fn(ax, D)
        except Exception as e:
            traceback.print_exc()
            _empty(ax, f"panel {letter} failed: {type(e).__name__}")
            _letter(ax, letter)
        made.append(letter)
    # 右下角：口径说明卡（纯文本轴）
    axn = fig.add_subplot(gs[2, 1])
    axn.axis("off")
    axn.text(0.0, 0.97, "Reading notes", fontsize=6.6, fontweight="bold",
             va="top", ha="left", color="#111111")
    notes = [
        "• Two independent record sets, same layout:",
        "    a = 09b keep set (101 contigs, `novelty_9a`)",
        "    b = 09b rescue set (291 candidates, `category`)",
        "• Band width ∝ number of records.",
        "• New-species call (panel c rings):",
        "    length > genus average, OR nt identity",
        "    75–95 % with query coverage ≥ 70 %.",
        "• Genus → host-linked naming (`… lycii`, `…-like virus`,",
        "    `… associated`, `sp.`) is flagged in the manifest.",
        "• Panel c zones: <75 % below genus threshold;",
        "    75–95 % between genus and species (new sp.);",
        "    >95 % above species threshold (known sp.).",
    ]
    axn.text(0.0, 0.855, "\n".join(notes), fontsize=4.7, va="top", ha="left",
             color="#333333", linespacing=1.55)
    # 颜色键（放在独立说明轴上，与通路图元不同轴 → 结构上不可能重叠）
    axn.text(0.0, 0.29, "Colour key", fontsize=5.4, fontweight="bold",
             va="top", ha="left", color="#111111")
    from matplotlib.patches import Rectangle as _R
    ky = 0.225
    for t in TIER_ORDER:
        axn.add_patch(_R((0.0, ky - 0.020), 0.045, 0.040,
                         facecolor=NOVEL_COLORS[t], edgecolor=C_NODE_EDGE, lw=0.4,
                         transform=axn.transAxes, clip_on=False))
        axn.text(0.058, ky, t, fontsize=5.0, va="center", ha="left",
                 color="#333333", transform=axn.transAxes)
        ky -= 0.062
    fig.suptitle(
        f"Novelty and taxonomic placement of {D.virus_name} virome contigs",
        fontsize=8.6, fontweight="bold", y=0.982)
    png, pdf = save(fig, outdir, "Fig_GojiVirome_Novelty", dpi)
    plt.close(fig)
    return png, pdf, made


def build_panels(D: Data, outdir: Path, dpi: int):
    pdir = outdir / "panels_novelty"
    pdir.mkdir(parents=True, exist_ok=True)
    made = []
    for letter, fn in PANELS:
        if letter == "c":
            fig = plt.figure(figsize=(5.2, 3.6))
            ax = fig.add_axes([0.12, 0.14, 0.84, 0.76])
        else:
            fig = plt.figure(figsize=(7.6, 4.0))
            ax = fig.add_axes([0.02, 0.02, 0.96, 0.80])
        try:
            fn(ax, D)
        except Exception as e:
            traceback.print_exc()
            _empty(ax, f"panel {letter} failed: {type(e).__name__}")
            _letter(ax, letter)
        save(fig, pdir, f"panel_{letter}", dpi)
        plt.close(fig)
        made.append(letter)
    return made


def prune_stale(outdir: Path, keep: list[str]) -> list[Path]:
    """删掉 panels_novelty/ 下不属于 keep 的旧图（与主图同一条纪律）。"""
    pdir = outdir / "panels_novelty"
    if not pdir.is_dir():
        return []
    valid = {f"panel_{x}.{ext}" for x in keep for ext in ("png", "pdf")}
    removed = []
    for f in sorted(pdir.glob("panel_*.*")):
        if f.name not in valid:
            f.unlink()
            removed.append(f)
    return removed


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--virus-name", default="Lycium barbarum")
    ap.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    ap.add_argument("--panels-only", action="store_true")
    ap.add_argument("--composite-only", action="store_true")
    a = ap.parse_args()

    root = Path(a.root).expanduser().resolve()
    if not root.is_dir():
        sys.exit(f"[fatal] --root not a directory: {root}")
    outdir = Path(a.outdir).expanduser().resolve() if a.outdir \
        else (root / "10_Reports" / "S10_SCI_Discovery")

    # ⚠ 必须显式接管 rcParams（否则沿用进程级旧值 → 同代码两次跑出不同 dpi/字体嵌入）
    matplotlib.rcParams.update(sci_rcparams(base_size=7.5, grid=False))
    apply_sci_style(base_size=7.5, grid=False)

    D = Data(root, a.virus_name)
    print(f"[in]  root      = {root}")
    print(f"[in]  09b dir   = {D.b09}")
    for k, v in (D.manifest.get("_sources") or {}).items():
        print(f"[in]  {k:34s} = {v}")
    if not D.keep:
        print("[warn] keep_summary.tsv not found → panel a will be empty")
    if not D.rescue:
        print("[warn] rescue_evidence_scored.tsv not found → panels b/c will be empty")

    made = []
    if not a.panels_only:
        png, pdf, m = build_composite(D, outdir, a.dpi)
        print(f"[out] composite = {png.name} / {pdf.name}  panels={''.join(m)}")
        made += m
    if not a.composite_only:
        m = build_panels(D, outdir, a.dpi)
        print(f"[out] panels    = {len(m)} panels in {outdir / 'panels_novelty'}")
        stale = prune_stale(outdir, [k for k, _ in PANELS])
        if stale:
            print(f"[prune] removed {len(stale)} stale: "
                  f"{', '.join(f.name for f in stale)}")

    man = outdir / "_build" / "novelty_manifest.json"
    man.parent.mkdir(parents=True, exist_ok=True)
    man.write_text(json.dumps({"root": str(root), "virus_name": a.virus_name,
                               "panels": D.manifest},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[out] manifest  = {man}")


if __name__ == "__main__":
    main()
