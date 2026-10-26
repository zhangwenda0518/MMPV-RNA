#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
stat_virus_profile.py — 病毒层面画像：名录 / 丰度 / 跨样本检出率 + 绘图
======================================================================

为什么需要这个脚本
------------------
`stat_prescreen.py` 回答的是**样本级**问题（这条 run 要不要下载），
但病毒本身也是结果 —— 一个项目里"到底有哪些病毒、各占多少、在多少样本里出现"
是独立且更有科研价值的产出。本脚本把样本级的病毒明细长表**上卷**成病毒层面的画像。

输入（stat_prescreen.py 的产物）
-------------------------------
  stat_viral_taxa.tsv    每 run × 每病毒分类单元（长表，10864 行的量级）
  sample_decision.tsv    可选，用于取分母（有分类数据 / 有病毒信号的样本数）

输出
----
  virus_prevalence.tsv     病毒层面汇总：检出率 / 丰度 / 阶元 / 命中份额
  virus_sample_matrix.tsv  样本 × 病毒 矩阵（长表，供热图与下游统计）
  virus_top_by_sample.tsv  每样本的主要病毒（取自 sample_decision）
  VIRUS_PROFILE_REPORT.md  人读摘要
  figures/*.png            4 张单图 + 1 张 2x2 组合图

⚠ 层级重复是本脚本最核心的处理点
--------------------------------
`stat_viral_taxa.tsv` 里同一个病毒分支会同时出现父节点与子节点，例如某样本同时有
`Pospiviroidae`（科）/ `Pospiviroid`（属）/ `Potato spindle tuber viroid`（种）。
直接按名字计数会让"检出率"和"丰度"都被父节点虚高。

处理方式：**样本内判定末端节点（leaf）**，即该分支上实际检出的最细分类单元：
  * NCBI：`self_count == total_count 且 self_count > 0`
    （自身命中数等于含后代累计数 → 说明没有更细的后代被命中）
  * GSA ：`path` 不是任何其他节点 `path` 的前缀

汇总表同时给出 `n_samples`（该名字出现过）与 `n_leaf_samples`（它作为末端被检出），
下游可按需选口径。

用法
----
    # 基本：只出表，不绘图（无需 matplotlib）
    python stat_virus_profile.py -i ./prescreen -o ./prescreen/virus_profile

    # 出表 + 绘图（需 numpy/matplotlib/seaborn）
    python stat_virus_profile.py -i ./prescreen -o ./prescreen/virus_profile --plot

    # 调整绘图筛选
    python stat_virus_profile.py -i ./prescreen -o out --plot --min-samples 5 --top-n 30

    # 绘图口径：默认 leaf（末端节点，父子不并列）；all 则任意层级都画
    python stat_virus_profile.py -i ./prescreen -o out --plot --prevalence-by all
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
import json
import os
import sys
from collections import defaultdict
from datetime import datetime
from statistics import mean, median
from typing import Any, Dict, List, Optional, Tuple

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# 复用 stat_prescreen 的阶元推断，避免两处逻辑分叉
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from stat_prescreen import rank_guess  # type: ignore
except Exception:  # pragma: no cover
    def rank_guess(name: str) -> str:  # type: ignore
        return "unknown"


# ============================================================================
# 病毒大类：viroid 是植物病毒研究的重点，单独摘出来
# ============================================================================
CLASS_ORDER = ("viroid", "virus", "phage", "virophage", "satellite", "other")
CLASS_COLOR = {
    "viroid":     "#D1495B",   # 类病毒 —— 植物病原，重点
    "virus":      "#2E6F9E",
    "phage":      "#7A9E9F",
    "virophage":  "#B08BBB",
    "satellite":  "#E5A54B",
    "other":      "#9AA0A6",
}


def virus_class(name: str) -> str:
    """按名字粗分病毒大类（只用于分组着色，不做分类学断言）。"""
    n = (name or "").lower()
    if not n:
        return "other"
    if "viroid" in n:
        return "viroid"
    if "virophage" in n:
        return "virophage"
    if "phage" in n:
        return "phage"
    if "satellite" in n:
        return "satellite"
    if "incertae sedis" in n or n in ("viruses", "unclassified"):
        return "other"
    return "virus"


# ============================================================================
# 读取
# ============================================================================
def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _to_int(v: Any, d: int = 0) -> int:
    try:
        return int(float(str(v).strip() or 0))
    except (TypeError, ValueError):
        return d


def _to_float(v: Any, d: float = 0.0) -> float:
    try:
        return float(str(v).strip() or 0)
    except (TypeError, ValueError):
        return d


def mark_leaf(rows: List[Dict[str, str]]) -> None:
    """在**样本内**标记每个病毒节点是否为末端（该分支实际检出的最细单元）。

    直接改 rows（原地加 `_leaf` 键）。
    """
    by_acc: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for r in rows:
        by_acc[r.get("acc", "")].append(r)

    for acc, rs in by_acc.items():
        gsa_paths = [r.get("path", "") for r in rs if r.get("source") == "GSA"]
        for r in rs:
            if r.get("source") == "GSA":
                p = r.get("path", "") or ""
                # path 以 " > " 连接；若没有任何其他节点以它为前缀，它就是末端
                r["_leaf"] = not any(
                    o and o != p and o.startswith(p + " > ") for o in gsa_paths)
            else:
                s = _to_int(r.get("self_count"))
                t = _to_int(r.get("total_count"))
                r["_leaf"] = (s > 0 and s == t)


# ============================================================================
# 聚合
# ============================================================================
def aggregate(rows: List[Dict[str, str]], n_ok: int, n_viral: int) -> List[Dict[str, Any]]:
    """把 run×taxon 长表上卷成「每个病毒分类单元一行」。"""
    agg: Dict[str, Dict[str, Any]] = {}

    for r in rows:
        name = (r.get("org") or "").strip()
        if not name:
            continue
        a = agg.setdefault(name, {
            "virus": name,
            "virus_class": virus_class(name),
            "rank": rank_guess(name),
            "tax_id": (r.get("tax_id") or "").strip(),
            "is_leaf_node": False,
            "samples": set(), "leaf_samples": set(),
            "pcts": [], "hits": 0, "leaf_hits": 0,
            "sources": set(),
        })
        acc = r.get("acc", "")
        a["samples"].add(acc)
        a["sources"].add(r.get("source", ""))
        if r.get("tax_id") and not a["tax_id"]:
            a["tax_id"] = r["tax_id"].strip()
        if r.get("_leaf"):
            a["is_leaf_node"] = True
            a["leaf_samples"].add(acc)
            a["leaf_hits"] += _to_int(r.get("total_count"))
        a["pcts"].append(_to_float(r.get("pct")))
        a["hits"] += _to_int(r.get("total_count"))

    # 命中份额：以「末端节点」的累计命中为分母 —— 避免层级重复把父节点算重
    total_leaf_hits = sum(a["leaf_hits"] for a in agg.values()) or 0

    out: List[Dict[str, Any]] = []
    for name, a in agg.items():
        ns = len(a["samples"])
        nl = len(a["leaf_samples"])
        pcts = [p for p in a["pcts"] if p > 0]
        out.append({
            "virus": name,
            "virus_class": a["virus_class"],
            "rank": a["rank"],
            "tax_id": a["tax_id"],
            "is_leaf_node": "yes" if a["is_leaf_node"] else "no",
            "n_samples": ns,
            "prevalence_pct": round(100.0 * ns / n_ok, 4) if n_ok else 0.0,
            "n_leaf_samples": nl,
            "leaf_prevalence_pct": round(100.0 * nl / n_ok, 4) if n_ok else 0.0,
            "prevalence_in_viral_pct": round(100.0 * ns / n_viral, 4) if n_viral else 0.0,
            "mean_pct": round(mean(pcts), 6) if pcts else 0.0,
            "median_pct": round(median(pcts), 6) if pcts else 0.0,
            "max_pct": round(max(pcts), 6) if pcts else 0.0,
            "total_hits": a["hits"],
            "leaf_hits": a["leaf_hits"],
            "hits_share_pct": (round(100.0 * a["leaf_hits"] / total_leaf_hits, 4)
                               if total_leaf_hits else 0.0),
            "sources": "+".join(sorted(s for s in a["sources"] if s)),
        })

    out.sort(key=lambda x: (-x["n_samples"], -x["total_hits"]))
    return out


def build_matrix(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """样本 × 病毒 长表（供热图与下游统计）。"""
    m = []
    for r in rows:
        name = (r.get("org") or "").strip()
        if not name:
            continue
        m.append({
            "acc": r.get("acc", ""),
            "source": r.get("source", ""),
            "virus": name,
            "virus_class": virus_class(name),
            "rank": rank_guess(name),
            "is_leaf": "yes" if r.get("_leaf") else "no",
            "pct": round(_to_float(r.get("pct")), 6),
            "hits": _to_int(r.get("total_count")),
        })
    m.sort(key=lambda x: (x["virus"], x["acc"]))
    return m


# ============================================================================
# 输出表
# ============================================================================
PREV_COLS = [
    "virus", "virus_class", "rank", "tax_id", "is_leaf_node",
    "n_samples", "prevalence_pct", "prevalence_in_viral_pct",
    "n_leaf_samples", "leaf_prevalence_pct",
    "mean_pct", "median_pct", "max_pct",
    "total_hits", "leaf_hits", "hits_share_pct", "sources",
]


def write_tsv(path: str, cols: List[str], rows: List[Dict[str, Any]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})


# ============================================================================
# 绘图（无 matplotlib 时整体降级，不影响出表）
# ============================================================================
def _setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    # CJK 字体置前 (病毒名/类目可能含中文), 原西文字体保留作回退
    mpl.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial", "Helvetica", "DejaVu Sans", "sans-serif"]
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["axes.spines.top"] = False
    mpl.rcParams["axes.edgecolor"] = "#5A6A7E"
    mpl.rcParams["axes.labelcolor"] = "#1A2332"
    mpl.rcParams["text.color"] = "#1A2332"
    mpl.rcParams["xtick.color"] = "#5A6A7E"
    mpl.rcParams["ytick.color"] = "#5A6A7E"
    mpl.rcParams["savefig.dpi"] = 300
    mpl.rcParams["figure.facecolor"] = "white"
    return plt


def _wrap(s: str, width: int = 34) -> str:
    import textwrap
    return "\n".join(textwrap.wrap(s, width)) or s


def _rel(path: str, start: str) -> str:
    """相对路径；跨盘符（Windows D: → C:）时退回绝对路径，不让 relpath 抛异常。"""
    try:
        return os.path.relpath(path, start)
    except ValueError:
        return path


def _fmt_pct(v: Any, nd: int = 1) -> str:
    """百分数的人类可读格式 —— 避免表里冒出 4.6e-05 这种科学计数法。

    只用于 md 报告; TSV/JSON 保留原始精度给下游分析。
    """
    x = _to_float(v)
    if x == 0:
        return "0"
    if abs(x) < 0.001:
        return "<0.001"
    if abs(x) < 1:
        return f"{x:.4f}".rstrip("0").rstrip(".")
    return f"{x:.{nd}f}"


def plot_prevalence_bar(plt, df, out: str, top_n: int = 25, scope: str = "leaf-node") -> None:
    """图 1：Top N 病毒跨样本检出率（横向条形）。"""
    d = df.head(top_n).iloc[::-1]
    h = max(4.0, 0.34 * len(d) + 1.6)
    fig, ax = plt.subplots(figsize=(9.2, h))
    colors = [CLASS_COLOR.get(c, "#9AA0A6") for c in d["virus_class"]]
    ax.barh(range(len(d)), d["prevalence_pct"], color=colors, height=0.72)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([_wrap(v) for v in d["virus"]], fontsize=8)
    xmax = max(d["prevalence_pct"].max(), 1.0)
    for y, (p, n) in enumerate(zip(d["prevalence_pct"], d["n_samples"])):
        ax.text(p + xmax * 0.012, y, f"{p:.1f}%  (n={int(n)})",
                va="center", fontsize=7.5, color="#33475B")
    ax.set_xlim(0, xmax * 1.22)
    ax.set_xlabel(f"Prevalence  (% of samples with taxonomy data, {scope})", fontsize=9.5)
    ax.set_title(f"Top {len(d)} viruses by cross-sample prevalence", fontsize=11, pad=10)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="x", color="#E8ECF1", lw=0.7)
    ax.set_axisbelow(True)
    # 图例
    from matplotlib.patches import Patch
    present = [c for c in CLASS_ORDER if c in set(d["virus_class"])]
    ax.legend(handles=[Patch(facecolor=CLASS_COLOR[c], label=c) for c in present],
              fontsize=7.5, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_prevalence_vs_abundance(plt, df, out: str, top_n: int = 60,
                                 scope: str = "leaf-node") -> None:
    """图 2：检出率 vs 平均丰度（气泡图，点大小 = 样本数）。"""
    d = df.head(top_n)
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    colors = [CLASS_COLOR.get(c, "#9AA0A6") for c in d["virus_class"]]
    sizes = [18 + 3.2 * n for n in d["n_samples"]]
    ax.scatter(d["prevalence_pct"], d["mean_pct"].clip(lower=1e-5),
               s=sizes, c=colors, alpha=0.72, edgecolors="white", linewidths=0.6)
    ax.set_xlabel(f"Prevalence (% of samples, {scope})", fontsize=9.5)
    ax.set_ylabel("Mean abundance within sample (%, log scale)", fontsize=9.5)
    ax.set_yscale("log")
    ax.set_title("Prevalence vs. abundance", fontsize=11, pad=10)
    ax.grid(color="#EEF1F5", lw=0.7)
    ax.set_axisbelow(True)
    # 标注检出率最高的几个 + 丰度最高的几个
    lab = list(dict.fromkeys(
        list(d.sort_values("prevalence_pct", ascending=False).head(4)["virus"]) +
        list(d.sort_values("mean_pct", ascending=False).head(3)["virus"])))
    for _, r in d.iterrows():
        if r["virus"] in lab:
            ax.annotate(r["virus"], (r["prevalence_pct"], max(r["mean_pct"], 1e-5)),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=7, color="#33475B")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(plt, matrix: List[Dict[str, Any]], out: str,
                 prev_df, top_v: int = 20, top_s: int = 40) -> None:
    """图 3：Top 病毒 × Top 样本 热图（log1p 变换的丰度）。"""
    import numpy as np
    import pandas as pd
    import seaborn as sns

    vs = list(prev_df.head(top_v)["virus"])
    dfm = pd.DataFrame(matrix)
    dfm = dfm[dfm["virus"].isin(vs)]
    # 样本按「检出病毒数」排序，让信息量大的样本排前面
    cnt = dfm.groupby("acc")["virus"].nunique().sort_values(ascending=False)
    ss = list(cnt.head(top_s).index)
    dfm = dfm[dfm["acc"].isin(ss)]

    if dfm.empty:
        return
    mat = dfm.pivot_table(index="virus", columns="acc", values="pct",
                          aggfunc="max").reindex(index=vs, columns=ss)
    mat = np.log10(mat.fillna(0).astype(float) + 1e-5)

    fig, ax = plt.subplots(figsize=(max(8.0, 0.26 * len(ss) + 3.2),
                                    max(5.0, 0.32 * len(vs) + 2.2)))
    sns.heatmap(mat, cmap="rocket_r", ax=ax, cbar_kws={"label": "log10(abundance %)", "shrink": 0.62},
                linewidths=0.25, linecolor="#F2F4F7",
                xticklabels=True, yticklabels=True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
    ax.set_yticklabels([_wrap(t.get_text(), 30) for t in ax.get_yticklabels()], fontsize=7.5)
    ax.set_xlabel(f"Samples (top {len(ss)} by distinct virus count)", fontsize=9.5)
    ax.set_ylabel("")
    ax.set_title("Virus abundance matrix", fontsize=11, pad=10)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_class_composition(plt, df, out: str, scope: str = "leaf-node",
                           min_samples: int = 2) -> None:
    """图 4：病毒大类构成（左：分类单元数；右：检出事件数）。

    注意本图只统计**进入绘图筛选的**分类单元（口径 + 样本数下限），
    与 md 报告里"全部分类单元"的构成表不是同一个分母 —— 标题里写明，
    避免同一份材料里两个数对不上。
    """
    from collections import Counter
    cnt_taxa = Counter(df["virus_class"])
    cnt_ev = Counter()
    for c, n in zip(df["virus_class"], df["n_samples"]):
        cnt_ev[c] += n
    order = [c for c in CLASS_ORDER if c in cnt_taxa]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4))
    for ax, cnt, title, xlabel in (
            (axes[0], cnt_taxa, "Distinct taxa", "Number of distinct viral taxa"),
            (axes[1], cnt_ev, "Detection events", "Sum of samples per taxon")):
        vals = [cnt[c] for c in order]
        cols = [CLASS_COLOR[c] for c in order]
        ax.bar(range(len(order)), vals, color=cols, width=0.66)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, fontsize=9)
        for i, v in enumerate(vals):
            ax.text(i, v + max(vals) * 0.02, str(v), ha="center", fontsize=8.5, color="#33475B")
        ax.set_title(title, fontsize=10.5)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.grid(axis="y", color="#EEF1F5", lw=0.7)
        ax.set_axisbelow(True)
        ax.set_ylim(0, max(vals) * 1.18 if vals else 1)
    fig.suptitle(f"Viral composition by class  "
                 f"({scope}, taxa in ≥{min_samples} samples: n={len(df)})",
                 fontsize=11.5, y=0.99)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_combined(plt, df, matrix, out: str, top_n: int = 18,
                  scope: str = "leaf-node") -> None:
    """组合图（2x2），便于直接放进报告/投稿。"""
    import numpy as np
    import pandas as pd
    import seaborn as sns

    fig = plt.figure(figsize=(14.5, 10.5))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.26,
                          width_ratios=[1.15, 1.0])

    # (A) 检出率
    ax = fig.add_subplot(gs[0, 0])
    d = df.head(top_n).iloc[::-1]
    cols = [CLASS_COLOR.get(c, "#9AA0A6") for c in d["virus_class"]]
    ax.barh(range(len(d)), d["prevalence_pct"], color=cols, height=0.7)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([_wrap(v, 30) for v in d["virus"]], fontsize=7.5)
    xm = max(d["prevalence_pct"].max(), 1.0)
    for y, (p, n) in enumerate(zip(d["prevalence_pct"], d["n_samples"])):
        ax.text(p + xm * 0.015, y, f"{p:.0f}% ({int(n)})", va="center", fontsize=6.8, color="#33475B")
    ax.set_xlim(0, xm * 1.25)
    ax.set_xlabel(f"Prevalence (%, {scope})", fontsize=9)
    ax.set_title("A  Top viruses by prevalence", fontsize=10.5, loc="left", fontweight="bold")
    ax.grid(axis="x", color="#EEF1F5", lw=0.7); ax.set_axisbelow(True)

    # (B) 检出率 vs 丰度
    ax = fig.add_subplot(gs[0, 1])
    d2 = df.head(60)
    ax.scatter(d2["prevalence_pct"], d2["mean_pct"].clip(lower=1e-5),
               s=[16 + 3.0 * n for n in d2["n_samples"]],
               c=[CLASS_COLOR.get(c, "#9AA0A6") for c in d2["virus_class"]],
               alpha=0.72, edgecolors="white", linewidths=0.5)
    ax.set_yscale("log")
    ax.set_xlabel(f"Prevalence (%, {scope})", fontsize=9)
    ax.set_ylabel("Mean abundance (%, log)", fontsize=9)
    ax.set_title("B  Prevalence vs. abundance", fontsize=10.5, loc="left", fontweight="bold")
    ax.grid(color="#EEF1F5", lw=0.7); ax.set_axisbelow(True)

    # (C) 热图
    ax = fig.add_subplot(gs[1, :])
    vs = list(df.head(16)["virus"])
    dfm = pd.DataFrame(matrix)
    dfm = dfm[dfm["virus"].isin(vs)]
    cnt = dfm.groupby("acc")["virus"].nunique().sort_values(ascending=False)
    ss = list(cnt.head(40).index)
    dfm = dfm[dfm["acc"].isin(ss)]
    if not dfm.empty:
        mat = dfm.pivot_table(index="virus", columns="acc", values="pct",
                              aggfunc="max").reindex(index=vs, columns=ss)
        mat = np.log10(mat.fillna(0).astype(float) + 1e-5)
        sns.heatmap(mat, cmap="rocket_r", ax=ax, linewidths=0.25, linecolor="#F2F4F7",
                    cbar_kws={"label": "log10(abundance %)", "shrink": 0.5})
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=5.5)
        ax.set_yticklabels([_wrap(t.get_text(), 28) for t in ax.get_yticklabels()], fontsize=7)
        ax.set_xlabel("Samples", fontsize=9)
    ax.set_title("C  Virus x sample abundance matrix", fontsize=10.5, loc="left", fontweight="bold")

    from matplotlib.patches import Patch
    present = [c for c in CLASS_ORDER if c in set(df["virus_class"])]
    fig.legend(handles=[Patch(facecolor=CLASS_COLOR[c], label=c) for c in present],
               loc="upper right", fontsize=8, frameon=False, ncol=len(present),
               bbox_to_anchor=(0.99, 0.985))
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# 主流程
# ============================================================================
def main() -> int:
    p = argparse.ArgumentParser(
        description="病毒层面画像：名录 / 丰度 / 跨样本检出率 + 绘图",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python stat_virus_profile.py -i ./prescreen -o ./prescreen/virus_profile --plot
  python stat_virus_profile.py -i ./prescreen -o out --plot --min-samples 5 --top-n 30
""")
    p.add_argument("-i", "--indir", required=True,
                   help="stat_prescreen 输出目录（含 stat_viral_taxa.tsv）")
    p.add_argument("-o", "--outdir", default=None, help="输出目录（默认 <indir>/virus_profile）")
    p.add_argument("--decision", default=None,
                   help="sample_decision.tsv 路径（用于取分母；默认在 indir 下自动找）")
    p.add_argument("--plot", action="store_true", help="生成图表（需 numpy/matplotlib/seaborn）")
    p.add_argument("--min-samples", type=int, default=2,
                   help="绘图与名录的样本数下限（默认 2，滤掉单例噪声；表仍保留全部）")
    p.add_argument("--top-n", type=int, default=25, help="条形图展示的病毒数（默认 25）")
    p.add_argument("--prevalence-by", choices=["leaf", "all"], default="leaf",
                   help="绘图口径: leaf=末端节点（实际检出到的最细单元，默认，不会层级重复）；"
                        "all=任意层级（含祖先节点，父子会并列出现）")
    args = p.parse_args()

    indir = os.path.abspath(args.indir)
    outdir = os.path.abspath(args.outdir or os.path.join(indir, "virus_profile"))
    os.makedirs(outdir, exist_ok=True)

    taxa_path = os.path.join(indir, "stat_viral_taxa.tsv")
    if not os.path.isfile(taxa_path):
        print(f"[ERR] 找不到 {taxa_path}")
        return 2

    rows = read_tsv(taxa_path)
    if not rows:
        print(f"[ERR] {taxa_path} 为空")
        return 2
    mark_leaf(rows)

    # ---- 分母：优先从 sample_decision.tsv 取真实统计 ----
    dec_path = args.decision or os.path.join(indir, "sample_decision.tsv")
    n_ok = n_viral = 0
    if os.path.isfile(dec_path):
        for r in read_tsv(dec_path):
            if (r.get("status") or "").strip() == "OK":
                n_ok += 1
                if _to_float(r.get("virus_pct_ident")) > 0:
                    n_viral += 1
    if not n_ok:
        n_ok = len(set(r.get("acc", "") for r in rows))
        n_viral = n_ok

    prev = aggregate(rows, n_ok, n_viral)
    matrix = build_matrix(rows)

    write_tsv(os.path.join(outdir, "virus_prevalence.tsv"), PREV_COLS, prev)
    write_tsv(os.path.join(outdir, "virus_sample_matrix.tsv"),
              ["acc", "source", "virus", "virus_class", "rank", "is_leaf", "pct", "hits"],
              matrix)

    # ---- 每样本主要病毒（若上游有）----
    top_by_sample = []
    if os.path.isfile(dec_path):
        for r in read_tsv(dec_path):
            if _to_float(r.get("virus_pct_ident")) > 0:
                top_by_sample.append({
                    "acc": r.get("acc", ""), "source": r.get("source", ""),
                    "virus_pct": r.get("virus_pct", ""),
                    "virus_pct_ident": r.get("virus_pct_ident", ""),
                    "top_virus": r.get("top_virus", ""),
                    "top_virus_pct": r.get("top_virus_pct", ""),
                    "virus_signal": r.get("virus_signal", ""),
                    "organism_match": r.get("organism_match", ""),
                    "decision": r.get("decision", ""),
                })
        top_by_sample.sort(key=lambda x: -_to_float(x["virus_pct_ident"]))
        write_tsv(os.path.join(outdir, "virus_top_by_sample.tsv"),
                  ["acc", "source", "virus_pct", "virus_pct_ident", "top_virus",
                   "top_virus_pct", "virus_signal", "organism_match", "decision"],
                  top_by_sample)

    # ---- 摘要统计 ----
    filtered = [r for r in prev if r["n_samples"] >= args.min_samples]
    leaf_only = [r for r in prev if r["is_leaf_node"] == "yes"
                 and r["n_samples"] >= args.min_samples]
    by_class: Dict[str, int] = {}
    for r in prev:
        by_class[r["virus_class"]] = by_class.get(r["virus_class"], 0) + 1

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "input": _rel(taxa_path, outdir),
        "n_samples_with_taxonomy": n_ok,
        "n_samples_with_virus": n_viral,
        "n_viral_taxa_total": len(prev),
        "n_viral_taxa_min_samples": len(filtered),
        "n_leaf_taxa_min_samples": len(leaf_only),
        "min_samples": args.min_samples,
        "by_class": by_class,
        "top20_by_prevalence": [
            {"virus": r["virus"], "class": r["virus_class"], "rank": r["rank"],
             "n_samples": r["n_samples"], "prevalence_pct": r["prevalence_pct"],
             "mean_pct": r["mean_pct"], "hits_share_pct": r["hits_share_pct"]}
            for r in prev[:20]],
    }

    # ---- 绘图口径 ----
    # leaf：同一分支只保留**实际检出到的最细单元**，父子节点不再并列（默认）
    # all ：任意层级都画（含祖先节点，父子会同时出现，仅用于看层级结构）
    scope = "leaf-node" if args.prevalence_by == "leaf" else "any-rank"
    summary["prevalence_by"] = args.prevalence_by
    summary["plot_scope"] = scope

    # ---- 绘图 ----
    figures: List[str] = []
    if args.plot:
        try:
            plt = _setup_mpl()
            import pandas as pd
            df = pd.DataFrame(filtered)
            if args.prevalence_by == "leaf":
                # 把「检出率」整体切到末端节点口径；丰度列不变（样本内丰度与层级无关）
                df = df.copy()
                df["n_samples"] = df["n_leaf_samples"]
                df["prevalence_pct"] = df["leaf_prevalence_pct"]
                df = df[df["n_samples"] >= args.min_samples]
                df = df.sort_values(["n_samples", "leaf_hits"],
                                    ascending=False).reset_index(drop=True)
            summary["n_plot_taxa"] = int(len(df))
            if not df.empty:
                figdir = os.path.join(outdir, "figures")
                os.makedirs(figdir, exist_ok=True)
                jobs = [
                    ("01_virus_prevalence_bar.png",
                     lambda pth: plot_prevalence_bar(plt, df, pth, args.top_n, scope)),
                    ("02_prevalence_vs_abundance.png",
                     lambda pth: plot_prevalence_vs_abundance(plt, df, pth, 60, scope)),
                    ("03_virus_sample_heatmap.png",
                     lambda pth: plot_heatmap(plt, matrix, pth, df)),
                    ("04_virus_class_composition.png",
                     lambda pth: plot_class_composition(plt, df, pth, scope,
                                                        args.min_samples)),
                    ("Combined_Virus_Profile.png",
                     lambda pth: plot_combined(plt, df, matrix, pth, 18, scope)),
                ]
                for name, fn in jobs:
                    pth = os.path.join(figdir, name)
                    try:
                        fn(pth)
                        if os.path.isfile(pth):
                            figures.append(_rel(pth, outdir))
                    except Exception as e:
                        print(f"  [WARN] 绘图失败 {name}: {type(e).__name__}: {e}")
            summary["figures"] = figures
        except ImportError as e:
            print(f"  [WARN] 缺少绘图依赖（{e}），已跳过绘图；表已正常产出")
            summary["figures"] = []
    else:
        summary["figures"] = []

    with open(os.path.join(outdir, "virus_profile_meta.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    # ---- 人读报告 ----
    rep = os.path.join(outdir, "VIRUS_PROFILE_REPORT.md")
    with open(rep, "w", encoding="utf-8") as f:
        f.write("# 病毒层面画像（名录 / 丰度 / 跨样本检出率）\n\n")
        f.write(f"- 生成时间: {summary['generated_at']}\n")
        f.write(f"- 输入: `{summary['input']}`\n")
        f.write(f"- 有分类数据的样本: **{n_ok}**\n")
        f.write(f"- 检出病毒的样本: **{n_viral}**\n")
        f.write(f"- 病毒分类单元总数: **{len(prev)}**"
                f"（任意层级口径，样本数 ≥ {args.min_samples} 的 {len(filtered)} 个）\n")
        f.write(f"- 其中**末端节点**（实际检出到的最细单元）: "
                f"**{sum(1 for r in prev if r['is_leaf_node'] == 'yes')}**"
                f"（样本数 ≥ {args.min_samples} 的 {len(leaf_only)} 个）\n\n")

        f.write("## 口径说明（先读这段）\n\n")
        f.write("`stat_viral_taxa.tsv` 里同一分支会同时出现**父节点与子节点**"
                "（如 `Pospiviroidae` 科 / `Pospiviroid` 属 / `Potato spindle tuber viroid` 种）。\n"
                "直接按名字计数会让父节点的检出率与丰度虚高，因此本表同时给出两套口径：\n\n")
        f.write("| 列 | 含义 |\n|---|---|\n")
        f.write("| `n_samples` | 该名字出现在多少个样本里（含作为祖先出现的） |\n")
        f.write("| `n_leaf_samples` | 该名字作为**末端节点**被检出的样本数（更严格） |\n")
        f.write("| `prevalence_pct` | `n_samples` / 有分类数据的样本数 |\n")
        f.write("| `leaf_prevalence_pct` | `n_leaf_samples` / 有分类数据的样本数 |\n")
        f.write("| `prevalence_in_viral_pct` | `n_samples` / **检出病毒**的样本数 |\n")
        f.write("| `mean_pct` / `median_pct` / `max_pct` | 样本内丰度（各来源自身口径） |\n")
        f.write("| `hits_share_pct` | 该病毒在**全部末端命中**中的份额（避免层级重复） |\n")
        f.write("| `rank` | 从学名字形推断的阶元（species/genus/family/order/...） |\n\n")
        f.write("> ⚠️ 两个来源分母不同：NCBI STAT 是 analysed（已鉴定 reads），"
                "GSA 是总 reads。跨来源比较请用 `stat_run_summary.tsv` 的 `virus_pct_ident`。\n\n")

        f.write(f"## 病毒大类构成（共 {len(prev)} 个分类单元）\n\n")
        f.write(f"> 这里是**全部**分类单元的构成；`04_virus_class_composition.png` "
                f"只画进入绘图筛选的那部分（口径 {scope}、样本数 ≥ {args.min_samples}），"
                f"分母更小，两个数不一样是正常的。\n\n")
        f.write("| 大类 | 分类单元数 |\n|---|---|\n")
        for c in CLASS_ORDER:
            if by_class.get(c):
                f.write(f"| {c} | {by_class[c]} |\n")
        f.write("\n")

        f.write(f"## 检出率 Top 30（**任意层级**口径，样本数 ≥ {args.min_samples}）\n\n")
        f.write("> 这一节按 `n_samples`（该名字出现过）排序，"
                "**含只作为祖先出现的节点**（如 `Pospiviroidae` 科），"
                "因此数值会偏高，仅用于看有哪些分支被命中。\n"
                "> 要做病毒名录/流行率判断，请用下一节的**末端节点**口径。\n\n")
        f.write("| 病毒 | 大类 | 阶元 | 检出样本数 | 检出率% | 末端检出数 | 平均丰度% | 命中份额% |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for r in filtered[:30]:
            f.write(f"| {r['virus']} | {r['virus_class']} | {r['rank']} | {r['n_samples']} | "
                    f"{_fmt_pct(r['prevalence_pct'])} | {r['n_leaf_samples']} | "
                    f"{_fmt_pct(r['mean_pct'], 3)} | {_fmt_pct(r['hits_share_pct'], 2)} |\n")
        f.write("\n")

        leaf_rows = sorted((r for r in leaf_only if r["n_leaf_samples"] >= args.min_samples),
                           key=lambda x: (-x["n_leaf_samples"], -x["leaf_hits"]))[:20]
        if leaf_rows:
            f.write("## 末端节点 Top 20（实际检出到的最细分类单元）\n\n")
            f.write("| 病毒 | 阶元 | 末端检出样本数 | 末端检出率% | 平均丰度% |\n|---|---|---|---|---|\n")
            for r in leaf_rows:
                f.write(f"| {r['virus']} | {r['rank']} | {r['n_leaf_samples']} | "
                        f"{_fmt_pct(r['leaf_prevalence_pct'])} | "
                        f"{_fmt_pct(r['mean_pct'], 3)} |\n")
            f.write("\n")

        if figures:
            f.write("## 图\n\n")
            f.write(f"绘图口径：**{scope}**"
                    + ("（只统计实际检出到的最细分类单元，父子节点不会并列，"
                       "因此不会虚高）\n\n" if args.prevalence_by == "leaf"
                       else "（含祖先节点，父子会并列出现，仅用于观察层级结构）\n\n"))
            for g in figures:
                f.write(f"![{os.path.basename(g)}]({g})\n\n")

        f.write("## 注意\n\n")
        f.write("- 这是 **k-mer 快速分类**的结果，不是严谨比对：低丰度病毒会漏、"
                "宿主内源病毒序列会误报。**定位是预筛**，结论请用下游比对/组装复核。\n")
        f.write("- 检出率高不等于该病毒真的流行 —— 也可能是数据库背景"
                "（如 `Human gammaherpesvirus 4`、噬菌体类在 RNA-seq 里常见）。\n")
        f.write("- 判读时建议结合 `virus_class`：`viroid` 是植物病原，"
                "参考价值通常高于噬菌体/人类病毒背景。\n")

    # ---- 控制台 ----
    print()
    print("=" * 60)
    print("  病毒层面画像")
    print("=" * 60)
    print(f"  有分类数据样本 : {n_ok}")
    print(f"  检出病毒样本   : {n_viral}")
    print(f"  病毒分类单元   : {len(prev)} 个（≥{args.min_samples} 样本的 {len(filtered)} 个）")
    print(f"  大类构成       : " + ", ".join(f"{c}={by_class[c]}" for c in CLASS_ORDER
                                            if by_class.get(c)))
    # 控制台头条与绘图保持同一口径，避免"图上一个数、屏上另一个数"
    if args.prevalence_by == "leaf":
        headline = sorted((r for r in prev if r["n_leaf_samples"] >= args.min_samples),
                          key=lambda x: (-x["n_leaf_samples"], -x["leaf_hits"]))[:10]
        hk, hl = "n_leaf_samples", "leaf_prevalence_pct"
    else:
        headline = filtered[:10]
        hk, hl = "n_samples", "prevalence_pct"
    print()
    print(f"  检出率 Top 10（口径 = {scope}）:")
    for r in headline:
        print(f"    {r[hk]:>4} 样本  {r[hl]:>6.1f}%  "
              f"[{r['virus_class']:<8}] {r['virus']}")
    print()
    print(f"  汇总表   -> {os.path.join(outdir, 'virus_prevalence.tsv')}")
    print(f"  矩阵表   -> {os.path.join(outdir, 'virus_sample_matrix.tsv')}")
    if top_by_sample:
        print(f"  样本表   -> {os.path.join(outdir, 'virus_top_by_sample.tsv')}")
    print(f"  元信息   -> {os.path.join(outdir, 'virus_profile_meta.json')}")
    print(f"  报告     -> {rep}")
    if figures:
        print(f"  图       -> {os.path.join(outdir, 'figures')}  "
              f"({len(figures)} 张, 口径={scope})")
    elif args.plot:
        print("  图       -> 未生成（依赖缺失或数据不足）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
