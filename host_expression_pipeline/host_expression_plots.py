#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""host_expression_plots.py — ⑧ 宿主表达结果出版级统计图 (09_Report/figures/)
=========================================================================
从 ⑧ 各阶段 TSV 产物出图, 遵循期刊规范 (figure-maker 约定):
  - Okabe-Ito 色盲安全配色; Arial; 300 dpi PNG + 矢量 PDF
  - 单栏 85 mm / 双栏 170 mm; 去顶右框线; 面板字母加粗

图清单:
  fig01_hub_ranking      Δz hub Top-N 横向棒棒糖图 (Ranking)
  fig02_module_load      模块×病毒载量关联气泡矩阵 (Correlation;
                         size=-log10 Padj, color=Spearman Rho 发散色)
  fig03_enrichment       DE_up/DE_down 富集 dot plot 分面 (Omics)
  fig04_module_network   模块网络布局 (每面: 全网络概览 + Top-3 模块细分;
                         python-igraph drl/FR 布局, 边稀疏化)

用法:
  python host_expression_plots.py --root <08_HostExpression> [--top-n 15]
输出:
  <root>/09_Report/figures/fig0*_*.png|pdf   (缺数据的图自动跳过)
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── 期刊规范常量 ──────────────────────────────────────────────────────
OKABE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
         "#56B4E9", "#F0E442", "#000000"]
MM = 1 / 25.4                      # mm → inch
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def _save(fig, out: Path, name: str) -> Path:
    fig.savefig(out / f"{name}.png")
    fig.savefig(out / f"{name}.pdf")
    plt.close(fig)
    print(f"  ✓ {out / (name + '.png')}")
    return out / f"{name}.png"


def _panel_label(ax, letter: str) -> None:
    ax.text(-0.10, 1.16, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left")


def _label_halo(ax, x, y, text, fontsize=5.2):
    """hub 标签: 白色描边 + 轻微错位避让 (matplotlib path_effects)。"""
    import matplotlib.patheffects as pe
    ax.annotate(text, (x, y), fontsize=fontsize, color="#111111",
                xytext=(3, 3), textcoords="offset points", zorder=5,
                path_effects=[pe.withStroke(linewidth=1.4, foreground="white")])


# ── fig01: Δz hub 排名 ────────────────────────────────────────────────

def fig_hub_ranking(root: Path, out: Path, top_n: int = 15):
    """Δz 全分布直方图 (阈值线 + 饱和天花板注记 + 增强/减弱计数)。

    注: Top-N 对/聚合表在 Top-20 截断方案下均饱和于 Δz 天花板
    (pos 秩1边与 neg 极值的差), 无区分度 —— 分布图如实呈现方法行为。"""
    cand = root / "05_Coexp" / "diff_coexp_candidates_sampled.tsv"
    summ = root / "05_Coexp" / "diff_coexp_summary.tsv"
    if not (cand.is_file() and summ.is_file()):
        print("  - fig01 跳过 (无 diff_coexp_candidates_sampled/summary)")
        return None
    dz = pd.read_csv(cand, sep="	", usecols=["Delta_z"])["Delta_z"]
    dz = dz.dropna()
    sm = pd.read_csv(summ, sep="	").iloc[0]
    fig, ax = plt.subplots(figsize=(95 * MM, 62 * MM))
    ax.hist(dz, bins=80, color="#c9d7e8", edgecolor="#7d97b5",
            linewidth=0.3, zorder=1)
    lo, hi = float(sm["threshold_low"]), float(sm["threshold_high"])
    for x, lab in ((lo, f"p{sm['q_low']:g} = {lo:.2f}"),
                   (hi, f"p{sm['q_high']:g} = {hi:.2f}")):
        ax.axvline(x, color=OKABE[3], lw=0.9, ls="--", zorder=3)
        ax.text(x, ax.get_ylim()[1] * 0.97, lab, rotation=90,
                va="top", ha="right", fontsize=5.4, color=OKABE[3])
    ceil = dz.max()
    if ceil > hi + 0.3:
        ax.axvline(ceil, color="#999999", lw=0.8, zorder=3)
        ax.text(ceil - 0.03, ax.get_ylim()[1] * 0.80,
                "saturation ceiling\n|Δz| = " + f"{ceil:.2f}", rotation=90,
                va="top", ha="right", fontsize=5.2, color="#666666")
    ax.axvspan(lo, hi, color="none", zorder=0)
    ax.set_xlabel("Δz = z(pos network) − z(neg network)")
    ax.set_ylabel("Candidate pairs (sampled)")
    ax.set_title("Δz distribution "
                 f"({int(sm['n_enhanced'])} enhanced / "
                 f"{int(sm['n_weakened'])} weakened significant pairs)",
                 loc="left", fontsize=7.5)
    ax.text(0.50, 0.86,
            "thresholds: 0.5% / 99.5% percentiles\n"
            "counts are threshold artifacts,\nnot biology (see report note)",
            transform=ax.transAxes, ha="center", fontsize=5.4,
            color="#666666")
    fig.tight_layout()
    return _save(fig, out, "fig01_delta_z_distribution")


# ── fig02: 模块×病毒载量关联气泡矩阵 ──────────────────────────────────

def fig_module_load(root: Path, out: Path):
    f = root / "05_Coexp" / "modules" / "module_load_assoc.tsv"
    if not f.is_file():
        print("  - fig02 跳过 (无 module_load_assoc.tsv)")
        return None
    df = pd.read_csv(f, sep="\t")
    df = df[df["Padj"] < 0.05]
    if df.empty:
        print("  - fig02 跳过 (无 Padj<0.05 关联)")
        return None
    df = df.copy()
    df["label"] = df["Group"].str.replace("positive", "pos",
                                          regex=False).str.replace(
        "negative", "neg", regex=False) + " M" + df["Module"].astype(str)
    df["size"] = (-np.log10(df["Padj"].clip(lower=1e-30))).clip(upper=15)
    row_order = (df.groupby("label")["Rho"]
                 .apply(lambda x: x.abs().mean()).sort_values().index.tolist())
    df["label"] = pd.Categorical(df["label"], categories=row_order,
                                 ordered=True)
    df = df.sort_values("label")
    virus_order = sorted(df["Virus"].unique())
    df["vx"] = df["Virus"].map({v: i for i, v in enumerate(virus_order)})
    n_row = df["label"].nunique()
    fig, ax = plt.subplots(figsize=(95 * MM, max(60, 7.5 * n_row) * MM))
    for _, r in df.iterrows():
        ax.scatter(r["vx"], r["label"], s=12 + 16 * r["size"],
                   c=[OKABE[0] if r["Rho"] > 0 else OKABE[1]],
                   edgecolor="white", lw=0.5, alpha=0.95, zorder=2)
        ax.text(r["vx"] + 0.13, r["label"], f"ρ={r['Rho']:+.2f}",
                va="center", fontsize=5.2, color="#444444", zorder=3)
        star = "**" if r["Padj"] < 1e-5 else ("*" if r["Padj"] < 1e-3 else "")
        if star:
            ax.text(r["vx"], r["label"], star, ha="center", va="center",
                    fontsize=5.5, color="white", fontweight="bold", zorder=4)
    ax.set_xticks(range(len(virus_order)))
    ax.set_xticklabels([v[:28] for v in virus_order], rotation=18,
                       ha="right")
    ax.set_xlim(-0.6, len(virus_order) - 0.2)
    ax.set_xlabel("Spearman ρ, module activity vs viral load\n"
                  r"* Padj < $10^{-3}$   ** Padj < $10^{-5}$   "
                  r"bubble size = $-\log_{10}$(Padj)")
    ax.set_ylabel("")
    ax.set_title("Module × viral load association", loc="left")
    h_col = [plt.Line2D([], [], marker="o", ls="", ms=5, color=OKABE[0],
                        label="ρ > 0"),
             plt.Line2D([], [], marker="o", ls="", ms=5, color=OKABE[1],
                        label="ρ < 0")]
    # 图例键直径与数据同标尺: 直径 pt = sqrt(s) = sqrt(12+16v)
    leg1 = ax.legend(handles=h_col, frameon=False, fontsize=5.8,
                     loc="center left", bbox_to_anchor=(1.01, 0.62),
                     handletextpad=0.4, borderaxespad=0)
    ax.add_artist(leg1)                                # 防第二次 legend 覆盖
    h_size = [plt.Line2D([], [], marker="o", ls="",
                         ms=(12 + 16 * v) ** 0.5, mfc="#cccccc",
                         mec="#999999",
                         label=r"$-\log_{10}$Padj = " + str(v))
              for v in (5, 10, 15)]
    ax.legend(handles=h_size, frameon=False, fontsize=5.8,
              loc="center left", bbox_to_anchor=(1.01, 0.28),
              handletextpad=0.6, borderaxespad=0, labelspacing=2.2)
    fig.tight_layout()
    return _save(fig, out, "fig02_module_load")


# ── fig03: 富集 dot plot (DE_up / DE_down 分面) ───────────────────────

def fig_enrichment(root: Path, out: Path, top_n: int = 10):
    frames = []
    for set_name, panel in (("DE_up", "up-regulated"), ("DE_down", "down-regulated")):
        f = root / "07_Enrich" / f"enrich_{set_name}.tsv"
        if not f.is_file():
            continue
        d = pd.read_csv(f, sep="\t").head(top_n)
        if d.empty:
            continue
        d["panel"] = panel
        d["label"] = d["Name"].fillna("").where(
            d["Name"].fillna("") != "", d["Term"])
        d["label"] = d["label"].str.slice(0, 44) + np.where(
            d["label"].str.len() > 44, "…", "")
        d["nlp"] = -np.log10(d["Padj"].clip(lower=1e-30))
        frames.append(d)
    if not frames:
        print("  - fig03 跳过 (无 enrich_DE_up/down.tsv)")
        return None
    df = pd.concat(frames, ignore_index=True)
    panels = list(dict.fromkeys(df["panel"]))
    vmin, vmax = df["nlp"].min(), df["nlp"].max()      # 全图统一色标
    fig, axes = plt.subplots(1, len(panels),
                             figsize=(170 * MM, 78 * MM), sharex=False)
    if len(panels) == 1:
        axes = [axes]
    for ax, panel in zip(axes, panels):
        sub = df[df["panel"] == panel].iloc[::-1]
        y = np.arange(len(sub))
        sc = ax.scatter(sub["Fold"], y, s=14 + 5 * sub["N_Set_in_Term"],
                        c=sub["nlp"], cmap="viridis", vmin=vmin, vmax=vmax,
                        edgecolor="white", lw=0.4, zorder=2)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["label"])
        ax.set_xlabel("Fold enrichment")
        ax.set_title(f"{panel} (top {len(sub)})", loc="left")
        ax.grid(axis="x", lw=0.3, color="#dddddd", zorder=0)
        ax.margins(x=0.10, y=0.10)                     # 大泡不出界
        _panel_label(ax, "A" if panel == panels[0] else "B")
    # 独立 colorbar 轴 (不压数据); size 图例放 colorbar 下方 (图外空白)
    fig.subplots_adjust(right=0.87, wspace=0.72)
    cax = fig.add_axes([0.895, 0.30, 0.012, 0.52])
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label(r"$-\log_{10}$(Padj)", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6)
    h_size = [plt.Line2D([], [], marker="o", ls="", ms=(14 + 5 * v) ** 0.5,
                         mfc="#cccccc", mec="#999999", label=f"{v} genes")
              for v in (5, 15, 25)]
    fig.legend(handles=h_size, frameon=False, fontsize=5.8,
               loc="lower right", bbox_to_anchor=(0.995, 0.02),
               title="Bubble size", title_fontsize=6, labelspacing=1.2)
    return _save(fig, out, "fig03_enrichment")


# ── fig04: 模块网络布局 (python-igraph 布局 + matplotlib 出版渲染) ─────

def _thin_edges(edges: pd.DataFrame, k_per_node: int) -> pd.DataFrame:
    """每节点保留 Z 最大的 k 条边 (并去重), 网络可视化的标准稀疏化。"""
    a = pd.concat([edges.rename(columns={"Gene1": "n", "Gene2": "m"}),
                   edges.rename(columns={"Gene2": "n", "Gene1": "m"})])
    a = a.sort_values("Z", ascending=False).groupby("n", sort=False).head(k_per_node)
    keep = set(zip(np.minimum(a["n"], a["m"]), np.maximum(a["n"], a["m"])))
    out = edges[[tuple(x) in keep for x in zip(edges["Gene1"], edges["Gene2"])]]
    return out.drop_duplicates(subset=["Gene1", "Gene2"])


def _draw_network(ax, nodes: list, edges: pd.DataFrame, node_color,
                  seed: int = 42, layout: str = "drl",
                  label_hubs: int = 0, edge_alpha: float = 0.10):
    import igraph as ig
    import random
    ig.set_random_number_generator(random.Random(seed))
    idx = {g: i for i, g in enumerate(nodes)}
    es = [(idx[u], idx[v]) for u, v in
          zip(edges["Gene1"], edges["Gene2"]) if u in idx and v in idx]
    g = ig.Graph(n=len(nodes), edges=es)
    lo = g.layout_drl() if layout == "drl" else \
        g.layout_fruchterman_reingold(niter=800)
    xy = np.array(lo.coords)
    deg = np.array(g.degree())
    segs = np.array([xy[a] for a, b in es] + [xy[b] for a, b in es])
    if len(segs):
        lc = matplotlib.collections.LineCollection(
            segs.reshape(-1, 2, 2), colors="#8a8a8a", linewidths=0.3,
            alpha=edge_alpha, zorder=1)
        ax.add_collection(lc)
    ax.scatter(xy[:, 0], xy[:, 1], s=3 + 1.8 * np.sqrt(deg),
               c=node_color, edgecolor="white", linewidths=0.15,
               alpha=0.92, zorder=2)
    if label_hubs:
        # 取度最高的 hub; 两 hub 过近时只标度更高的一个 (防文字贴死)
        placed = []
        for hi in list(np.argsort(-deg)[:label_hubs + 2]):
            if len(placed) >= label_hubs:
                break
            pt = xy[hi]
            if any(float(np.hypot(pt[0] - q[0], pt[1] - q[1])) < 0.06
                   for q in placed):
                continue
            placed.append(pt)
            _label_halo(ax, pt[0], pt[1], nodes[hi], fontsize=5.2)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.autoscale()
    ax.margins(0.06)
    return xy, deg


def fig_module_networks(root: Path, out: Path, top_modules: int = 3,
                        seed: int = 42):
    """每面 (positive/negative) 一张 2×2: A 全网络 (按模块着色) + B-D Top-3 模块。"""
    made = []
    for side in ("positive", "negative"):
        tf = root / "05_Coexp" / "modules" / f"module_table_{side}.tsv"
        ef = root / "05_Coexp" / "modules" / f"module_edges_{side}.tsv"
        if not (tf.is_file() and ef.is_file()):
            continue
        table = pd.read_csv(tf, sep="\t")
        edges = pd.read_csv(ef, sep="\t")
        sizes = table.groupby("Module")["Gene"].size().sort_values(
            ascending=False)
        top_ids = sizes.head(top_modules).index.tolist()
        cmap = {m: OKABE[i % len(OKABE)] for i, m in enumerate(
            sizes.index[:6])}
        fig, axes = plt.subplots(2, 2, figsize=(170 * MM, 168 * MM))
        # A: 全网络概览 (稀疏化, 按模块着色)
        gmap = dict(zip(table["Gene"], table["Module"]))
        thin = _thin_edges(edges, 3)
        nodes = sorted(set(thin["Gene1"]) | set(thin["Gene2"]))
        colors = [cmap.get(gmap.get(n), "#aaaaaa") for n in nodes]
        _draw_network(axes[0][0], nodes, thin, colors, seed=seed,
                      layout="drl", edge_alpha=0.10)
        axes[0][0].set_title(f"{side}: DEG network by module "
                             f"({len(nodes)} genes)", loc="left",
                             fontsize=7.5, pad=8)
        _panel_label(axes[0][0], "A")
        from collections import Counter
        drawn = Counter(gmap.get(n) for n in nodes)    # 按所绘子网计数
        h_mod = [plt.Line2D([], [], marker="o", ls="", ms=4.5,
                            color=cmap.get(m), label=f"M{m} ({drawn[m]})")
                 for m in sizes.index[:6] if drawn.get(m)]
        if drawn.get(None):
            h_mod.append(plt.Line2D([], [], marker="o", ls="", ms=4.5,
                                    color="#aaaaaa",
                                    label=f"other ({drawn[None]})"))
        axes[0][0].legend(handles=h_mod, frameon=False, fontsize=5.2,
                          loc="lower right", ncol=2, handletextpad=0.2,
                          columnspacing=0.8)
        # B-D: Top-3 模块细分
        for j, m in enumerate(top_ids):
            members = set(table.loc[table["Module"] == m, "Gene"])
            sub = edges[[u in members and v in members
                         for u, v in zip(edges["Gene1"], edges["Gene2"])]]
            sub = _thin_edges(sub, 6)
            ns = sorted(set(sub["Gene1"]) | set(sub["Gene2"]))
            ax = axes[(1 + j) // 2][(1 + j) % 2]
            _draw_network(ax, ns, sub, cmap.get(m, OKABE[0]), seed=seed + j,
                          layout="fr", label_hubs=2, edge_alpha=0.22)
            ax.set_title(f"module {m} ({len(members)} genes)", loc="left",
                         fontsize=7.5, pad=8)
            _panel_label(ax, "BCD"[j])
        fig.tight_layout()
        made.append(_save(fig, out, f"fig04_module_network_{side}"))
    if not made:
        print("  - fig04 跳过 (无 module_edges_*.tsv)")
    return made


def fig_circular_pathway(root: Path, out: Path, edges_per_node: int = 6,
                         seed: int = 42):
    """TEA 原文同款环形网络: 节点按模块排序上环、实线=正相关/虚线=负相关、
    每节点保留 Top 边防毛球。数据: DEG ∩ ko00196 天线蛋白通路 (阳性网络)。"""
    f = root / "05_Coexp" / "circular_pathway_edges.tsv"
    if not f.is_file():
        print("  - fig05 跳过 (无 circular_pathway_edges.tsv)")
        return None
    e = pd.read_csv(f, sep="	")
    genes = sorted(set(e["Gene1"]) | set(e["Gene2"]))
    # 功能标签 (eggNOG emapper: Preferred_name 优先, 否则 Description 归为 LHC)
    gname = {}
    ff = root / "05_Coexp" / "circular_pathway_funcs.tsv"
    if ff.is_file():
        fu = pd.read_csv(ff, sep="	", dtype=str)
        fu = fu[fu["Gene"].notna() & fu["Gene"].str.startswith("ningxia")]
        for _, r in fu.iterrows():
            nm = (r["Preferred_name"] or "").strip()
            if not nm and "light-harvesting" in (r["Description"] or ""):
                nm = "LHC"
            gname[r["Gene"]] = nm
    gmod = {}
    for _, r in e.iterrows():
        gmod[r["Gene1"]] = r["Module1"]
        gmod[r["Gene2"]] = r["Module2"]
    mod_of = sorted(set(gmod.values()))
    cmap = {m: OKABE[i % len(OKABE)] for i, m in enumerate(mod_of)}
    # 同功能相邻 (复刻原文的家族排序), 再按模块分组
    order = sorted(genes, key=lambda g: (gname.get(g, "zz"), g))
    n = len(order)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2
    pos = {g: (np.cos(a), np.sin(a)) for g, a in zip(order, ang)}
    # 每节点保留 Top 边
    a = pd.concat([e.rename(columns={"Gene1": "n", "Gene2": "m"}),
                   e.rename(columns={"Gene2": "n", "Gene1": "m"})])
    a = a.sort_values("Z", ascending=False).groupby("n", sort=False).head(
        edges_per_node)
    keep = {(min(x, y), max(x, y)) for x, y in zip(a["n"], a["m"])}
    ee = e[[pair in keep for pair in
            zip(e["Gene1"], e["Gene2"])]].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(120 * MM, 120 * MM))
    for _, r in ee.iterrows():
        x1, y1 = pos[r["Gene1"]]
        x2, y2 = pos[r["Gene2"]]
        if r["Cor"] > 0:
            ax.plot([x1, x2], [y1, y2], color="#8a8a8a", lw=0.9,
                    alpha=0.55, zorder=1, solid_capstyle="round")
        else:
            ax.plot([x1, x2], [y1, y2], color="#8a8a8a", lw=0.8,
                    alpha=0.55, ls=(0, (2, 2)), zorder=1)
    deg = pd.concat([ee["Gene1"], ee["Gene2"]]).value_counts()
    for g in order:
        x, y = pos[g]
        ax.scatter(x, y, s=90 + 22 * deg.get(g, 0), c=[cmap[gmod[g]]],
                   edgecolor="white", lw=0.8, zorder=3)
        lx, ly = x * 1.22, y * 1.22
        ha = "left" if x > 0.05 else ("right" if x < -0.05 else "center")
        nm = gname.get(g, "")
        if nm in ("-", "", "nan"):
            nm = ""
        short = g.replace("ningxia", "n")
        lab = short if not nm else short + "\n" + nm
        ax.text(lx, ly, lab, ha=ha, va="center", fontsize=5.8,
                color="#222222")
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.42, 1.42)
    ax.set_aspect("equal")
    ax.axis("off")
    h_mod = [plt.Line2D([], [], marker="o", ls="", ms=5.5, color=cmap[m],
                        label=f"module {m}") for m in mod_of]
    h_edge = [plt.Line2D([], [], color="#8a8a8a", lw=1.0,
                         label="positive co-exp."),
              plt.Line2D([], [], color="#8a8a8a", lw=1.0, ls=(0, (2, 2)),
                         label="negative co-exp.")]
    ax.legend(handles=h_mod + h_edge, frameon=False, fontsize=6,
              loc="center left", bbox_to_anchor=(1.02, 0.5),
              handletextpad=0.4)
    ax.set_title("ko00196 Photosynthesis - antenna proteins "
                 f"({n} DEGs, positive TEA-GCN)", loc="center",
                 fontsize=8)
    fig.tight_layout()
    return _save(fig, out, "fig05_circular_pathway")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="⑧ 结果出版级统计图")
    ap.add_argument("--root", required=True, help="08_HostExpression 输出根")
    ap.add_argument("--out-dir", help="图输出目录 (默认 <root>/09_Report/figures)")
    ap.add_argument("--top-n", type=int, default=15, help="hub/hub 与 term 各取前 N")
    args = ap.parse_args(argv)

    root = Path(args.root)
    out = Path(args.out_dir) if args.out_dir else \
        root / "09_Report" / "figures"
    out.mkdir(parents=True, exist_ok=True)
    print(f"输出: {out}")
    fig_hub_ranking(root, out, args.top_n)
    fig_module_load(root, out)
    fig_enrichment(root, out, args.top_n)
    fig_module_networks(root, out)
    fig_circular_pathway(root, out)
    print("完成")


if __name__ == "__main__":
    main()


def main_testable(argv):
    main(argv)
