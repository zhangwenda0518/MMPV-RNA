"""出版级静态树图 / UpSet 图封装 (PhyloSuite 内嵌库复用)。

依赖: biosoft/PhyloSuite/PhyloSuite 内嵌的 phytreeviz 与 marsilea 源码
      (本地直接 sys.path 引用, 不复制不安装; 服务器侧需将两库目录随包部署)。

两个入口:
    plot_mcc_tree(mcc_tree, out_pdf, ...)  — matplotlib 原生 MCC 树 (SCI 白底)
    plot_upset(sets_df, out_pdf, ...)      — marsilea 出版级静态 UpSet
                                             (virus×host 交叉分析产物直接可喂)

衔接:
    - plot_mcc_tree 的输入与 publication_figures.plot_mcc_tree 相同
      (Newick MCC 树 + tip_location_map, "[location=xxx]" comment 解析复用),
      但用 TreeViz 渲染, 避开 Bio.Phylo.draw 的排版局限。
    - plot_upset 的输入是 binary DataFrame (index=virus, columns=host/地点),
      与 phylo_cross / 病毒-宿主关联分析的宽表格式一致; 之前的交互式
      UpSet (8/23 集成) 保留用于探索, 本模块出静态 PDF 用于发表。
"""
import os
import sys
import logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# PhyloSuite 内嵌库路径 (相对管线根目录)
_PS_LIB = os.path.normpath(os.path.join(
    os.path.dirname(__file__), "..", "..", "biosoft", "PhyloSuite", "PhyloSuite"))


def _ensure_ps_path():
    if _PS_LIB not in sys.path:
        if not os.path.isdir(_PS_LIB):
            raise ImportError(
                f"PhyloSuite lib not found at {_PS_LIB}; "
                "deploy biosoft/PhyloSuite alongside the pipeline or pip install marsilea/phytreeviz")
        sys.path.insert(0, _PS_LIB)


# ═══════════════════════════════════════════════════════════════════
# MCC tree (phytreeviz.TreeViz)
# ═══════════════════════════════════════════════════════════════════

def parse_tip_locations(tree):
    """从 BEAST MCC newick comment 提取 {tip_name: location}。"""
    loc_map = {}
    for clade in tree.find_clades():
        if clade.is_terminal() and getattr(clade, "comment", None):
            comment = clade.comment.strip("[]")
            for part in comment.split(","):
                if "location=" in part:
                    loc_map[clade.name] = part.split("location=")[-1].strip().strip('"')
                    break
    return loc_map


# SCI 友好色板 (Okabe-Ito, 色盲安全)
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
             "#56B4E9", "#F0E442", "#000000"]


def plot_mcc_tree(mcc_tree, out_pdf, tip_location_map=None, title="",
                  figsize=(10, 6), dpi=150, show_confidence=True,
                  show_scale_bar=True):
    """用 phytreeviz 画 MCC 树, tip 按 location 着色。

    mcc_tree: Newick 文件路径或 Bio.Phylo 树对象
    返回 (out_pdf, loc_map)。
    """
    from Bio import Phylo
    _ensure_ps_path()
    from phytreeviz import TreeViz

    if isinstance(mcc_tree, (str, os.PathLike)):
        tree = Phylo.read(str(mcc_tree), "newick")
    else:
        tree = mcc_tree
    if tip_location_map is None:
        tip_location_map = parse_tip_locations(tree)
    locations = sorted(set(tip_location_map.values()))
    color_map = {loc: OKABE_ITO[i % len(OKABE_ITO)]
                 for i, loc in enumerate(locations)}

    tv = TreeViz(tree)
    if show_confidence:
        tv.show_confidence()
    if show_scale_bar:
        tv.show_scale_bar()
    if title:
        tv.set_title(title, fontproperties={"weight": "bold"})
    # tip label 按 location 着色
    for name, loc in tip_location_map.items():
        try:
            tv.set_node_label_props(
                name, color=color_map.get(loc, "#333333"), size=6)
        except Exception:
            pass  # 名称不完全匹配时跳过着色

    fig = tv.plotfig(dpi=dpi)
    fig.set_size_inches(*figsize)
    fig.savefig(out_pdf, dpi=dpi, format="pdf", facecolor="white",
                edgecolor="none", bbox_inches="tight")
    plt.close(fig)

    # legend 单独补一张 (TreeViz 不直接暴露 legend axes)
    if locations:
        fig_lg, ax_lg = plt.subplots(figsize=(2, 0.35 * len(locations) + 0.5))
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=color_map[loc], markersize=9, label=loc)
                   for loc in locations]
        ax_lg.legend(handles=handles, loc="center left", frameon=False,
                     title="Location")
        ax_lg.axis("off")
        lg_pdf = out_pdf.replace(".pdf", "_legend.pdf")
        fig_lg.savefig(lg_pdf, bbox_inches="tight", facecolor="white")
        plt.close(fig_lg)
    return out_pdf, tip_location_map


# ═══════════════════════════════════════════════════════════════════
# Static UpSet (marsilea)
# ═══════════════════════════════════════════════════════════════════

def plot_upset(sets_df, out_pdf, title="", min_cardinality=1,
               sort_subsets="cardinality", dpi=150, width=None, height=None):
    """marsilea 出版级静态 UpSet。

    sets_df: binary DataFrame — index=virus/条目名, columns=集合名
             (host/地点), 值 0/1。与病毒-宿主关联宽表直接兼容。
    返回 out_pdf。
    """
    _ensure_ps_path()
    from marsilea.upset import UpsetData, Upset

    data = UpsetData(sets_df)
    up = Upset(data, min_cardinality=min_cardinality,
               sort_subsets=sort_subsets,
               color=".15", width=width, height=height)
    up.render()
    if title:
        up.get_main_ax().set_title(title, fontweight="bold", loc="left")
    plt.savefig(out_pdf, dpi=dpi, format="pdf", facecolor="white",
                edgecolor="none", bbox_inches="tight")
    plt.close("all")
    return out_pdf


def plot_upset_from_long(rows, out_pdf, **kw):
    """从 (virus, host) 长表构造 UpSet。rows: iterable of (item, set_name)。"""
    import pandas as pd
    df = pd.DataFrame([(i, s) for i, s in rows], columns=["item", "set"])
    wide = (df.assign(v=1)
              .pivot_table(index="item", columns="set", values="v",
                           fill_value=0, aggfunc="first"))
    return plot_upset(wide, out_pdf, **kw)


def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="Static publication-grade UpSet from a long (item,set) TSV")
    ap.add_argument("--long-tsv", required=True,
                    help="两列 TSV: item<TAB>set")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    args = ap.parse_args()
    rows = []
    with open(args.long_tsv, encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0]:
                rows.append((parts[0], parts[1]))
    plot_upset_from_long(rows, args.out, title=args.title)
    print(f"UpSet: {args.out} ({len(rows)} item-set pairs)")


if __name__ == "__main__":
    main()
