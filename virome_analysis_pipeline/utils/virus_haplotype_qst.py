#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
virus_haplotype_qst.py — 准种型 (Quasispecies Type, QST) 单倍型形式化层 + token 化目录

定位
====
EasyHap 思想的 AF 兼容重写 (零 GPL 代码): 把 04_post_analysis 已有的 样本×位点 矩阵
升格为可检验的单倍型对象。现有 clustermap/PCA 只产出"图", 本脚本补上群体遗传学统计层:
模式命名 → 簇频率表 → 多样性 Hd → 物种私有 QST。

方法
====
1) Token 化: 每个样本在每个位点编码为离散 token。
   - 默认二值: 值 >= --af-min 记 '1', 否则 '0', 缺失记 'N'
   - --af-bins "0.05,0.25,0.75": 分级剂量 '0'(低于首阈)/'L'/'M'/'H', 保留准种梯度信息
2) QST 定义: 跨全部受检位点的 token 元组即一条单倍型模式; 相同模式归并。
3) 聚类: union-find 连通分量, 两模式 Hamming 距离 <= --hamming-threshold 即合并;
   距离只在双方均非 'N' 的位点计算 (缺失跳过), 共享位点数为 0 时永不合并。
   复杂度 O(u^2 * L), u 为唯一模式数。
4) 命名: 按携带样本数降序编号 QST001, QST002, ...
5) 统计: 每物种(需 -M 元数据)输出 n、检出 QST 数、单倍型多样性 Hd = n/(n-1)*(1-Σx_i²)、
   私有 QST 数 (成员全部来自单一物种且 >= --private-min)。

输入矩阵布局
============
与 variant_species_association.py 完全一致 (自动识别 snp_matrix.tsv / Matrix_05.csv)。

输出 (-o 目录)
==============
- QST_Definitions.tsv          每个 QST 的 token 模式与成员样本
- Sample_QST_Assignment.tsv    样本 → QST 映射
- QST_by_Host_Species.tsv      QST × 物种 计数表 (无元数据时跳过)
- QST_Diversity_Stats.tsv      各物种 n / QST 数 / Hd / 私有数 / 主导QST占比
- Private_QST_List.tsv         物种私有 QST 清单 (无元数据时跳过)
- Variant_Token_Catalog.tsv    位点级 token 编码目录 (#4 人读表)
- qst_token_heatmap.png/.pdf   样本×位点 token 热图, 按 QST 簇分块排列

用法示例
========
python utils/virus_haplotype_qst.py \
    -m 04_post_analysis/<virus>/vcf_merge/matrix/snp_matrix.tsv \
    -M sra_rna.data4/global_metadata/Global_Unified_Metadata_Core14.tsv \
    -o 04_post_analysis/<virus>/qst --label PSTVd --hamming-threshold 0.15
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

# Okabe-Ito 色板 (SCI 浅色白底规范)
TOKEN_COLORS = {
    "0": "#F5F5F5",
    "L": "#56B4E9",
    "M": "#E69F00",
    "H": "#D55E00",
    "1": "#0072B2",
    "N": "#BDBDBD",
}


# ==========================================
# 输入校验与矩阵装载 (与关联脚本同风格)
# ==========================================
class InputError(ValueError):
    pass


def _err(path, lineno, msg):
    raise InputError(f"{os.path.basename(path)}:{lineno}: {msg}")


def load_matrix(path):
    """读入 样本×位点 数值矩阵, 自动识别两种布局, 返回 (df, note)。"""
    if not os.path.isfile(path):
        raise InputError(f"矩阵文件不存在: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, index_col=0, dtype=str)
    idx_name = str(df.index.name or "").strip().lower()

    if idx_name == "sample":
        note = "layout=samples_rows (snp_matrix 风格)"
    else:
        df = df.T
        df.index.name = "Sample"
        note = "layout=sites_rows (Matrix_05 风格, 已转置)"

    dup_s = df.index[df.index.duplicated()].unique().tolist()
    if dup_s:
        _err(path, 0, f"重复样本名: {dup_s[:5]}")
    dup_v = df.columns[df.columns.duplicated()].unique().tolist()
    if dup_v:
        _err(path, 0, f"重复位点ID: {dup_v[:5]}")
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.shape[1] == 0:
        _err(path, 1, "未解析到任何位点列")
    return df, note


def parse_site_position(label):
    parts = str(label).strip().split("_")
    if len(parts) >= 3 and parts[-3].isdigit():
        return int(parts[-3])
    for p in reversed(parts):
        if p.isdigit():
            return int(p)
    import re
    m = re.findall(r"\d+", str(label))
    return int(m[-1]) if m else 0


def _resolve_host_col(meta, host_col):
    """宿主列自动检测: 显式指定优先, 否则在常见列名中找。"""
    if host_col:
        return host_col
    for cand in ("Host", "ScientificName", "Species", "host", "Adjusted_Species"):
        if cand in meta.columns:
            return cand
    return ""


def load_species_map(path, run_col, host_col, species_list):
    if not os.path.isfile(path):
        raise InputError(f"元数据文件不存在: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    meta = pd.read_csv(path, sep=sep, dtype=str)
    if run_col not in meta.columns:
        _err(path, 1, f"缺少必需列 '{run_col}'; 实际列: {list(meta.columns)[:15]}")
    host_col = _resolve_host_col(meta, host_col)
    if not host_col:
        _err(path, 1, f"未找到宿主列 (尝试 Host/ScientificName/Species); 实际列: {list(meta.columns)[:15]}")
    run_map, seen = {}, {}
    for i, row in meta.iterrows():
        run = str(row[run_col]).strip()
        if not run or run.lower() == "nan":
            continue
        if run in seen:
            _err(path, i + 2, f"Run 重复出现: '{run}' (首次在第 {seen[run]} 行)")
        seen[run] = i + 2
        h_low = str(row[host_col]).strip().lower()
        for name in species_list:
            if name.lower() in h_low:
                run_map[run] = name
                break
    return run_map


def load_group_map(path, run_col, group_col):
    """按任意元数据列分组 (--group-col): 返回 {Run: 列值}。"""
    if not os.path.isfile(path):
        raise InputError(f"元数据文件不存在: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    meta = pd.read_csv(path, sep=sep, dtype=str)
    for col in (run_col, group_col):
        if col not in meta.columns:
            _err(path, 1, f"缺少必需列 '{col}'; 实际列: {list(meta.columns)[:15]}")
    gmap, seen = {}, {}
    for i, row in meta.iterrows():
        run = str(row[run_col]).strip()
        if not run or run.lower() == "nan":
            continue
        if run in seen:
            _err(path, i + 2, f"Run 重复出现: '{run}' (首次在第 {seen[run]} 行)")
        seen[run] = i + 2
        gval = str(row[group_col]).strip()
        if gval and gval.lower() != "nan":
            gmap[run] = gval
    return gmap


# ==========================================
# Token 化
# ==========================================
def tokenize_matrix(M, af_min, af_bins):
    """返回 tokens (list[list[str]], 行=样本), 以及编码说明 legend 字符串。"""
    vals = M.to_numpy(dtype=float)
    if af_bins:
        thr = sorted(float(x) for x in af_bins.split(",") if x.strip())
        if len(thr) > 3:
            raise InputError(f"--af-bins 最多支持 3 个阈值, 收到 {len(thr)} 个")
        grades = ["L", "M", "H", "X"][: len(thr)]
        legend = _bins_legend(thr, grades)
        codes = np.full(vals.shape, grades[0], dtype=object)
        codes[np.isnan(vals)] = "N"
        for k in range(len(thr)):
            codes[vals >= thr[k]] = grades[k]
        return [[c for c in row] for row in codes], legend
    legend = f"二值编码: >={af_min:g}='1', <{af_min:g}='0', 缺失='N'"
    codes = np.where(np.isnan(vals), "N", np.where(vals >= af_min, "1", "0"))
    return [[c for c in row] for row in codes], legend


def _bins_legend(thr, grades):
    segs = [f"<{thr[0]:g}='0'"]
    for lo, hi, g in zip(thr[:-1], thr[1:], grades[1:]):
        segs.append(f"{lo:g}-{hi:g}='{g}'")
    segs.append(f">={thr[-1]:g}='{grades[0]}'")
    return "分级剂量: " + ", ".join(segs)


# ==========================================
# Union-Find 与 Hamming 聚类
# ==========================================
class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def hamming_ok(t1, t2, threshold):
    """共享非 N 位点上的归一化 Hamming; 无共享位点 → 不合并。"""
    shared = mism = 0
    for a, b in zip(t1, t2):
        if a != "N" and b != "N":
            shared += 1
            if a != b:
                mism += 1
    if shared == 0:
        return False
    return (mism / shared) <= threshold


def cluster_patterns(pattern_ids, threshold):
    """pattern_ids: 唯一模式列表 (tuple[str])。返回每个模式的簇号列表。"""
    u = len(pattern_ids)
    dsu = DSU(u)
    for i in range(u):
        pi = pattern_ids[i]
        for j in range(i + 1, u):
            if hamming_ok(pi, pattern_ids[j], threshold):
                dsu.union(i, j)
    roots = {}
    clusters = []
    for i in range(u):
        r = dsu.find(i)
        if r not in roots:
            roots[r] = len(clusters)
        clusters.append(roots[r])
    return clusters


# ==========================================
# 主流程
# ==========================================
def main():
    ap = argparse.ArgumentParser(
        description="QST 准种型单倍型形式化层 (token 化 + Hamming 聚类 + Hd)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-m", "--matrix", required=True,
                    help="样本×位点矩阵 (Matrix_05 或 snp_matrix.tsv)")
    ap.add_argument("-o", "--outdir", required=True, help="输出目录")
    ap.add_argument("-M", "--metadata", default="", help="可选元数据 TSV/CSV (启用物种层统计)")
    ap.add_argument("--run-col", default="Run")
    ap.add_argument("--host-col", default="",
                    help="元数据中宿主列名 (默认自动检测: Host/ScientificName/Species)")
    ap.add_argument("--species", nargs="+",
                    default=["Lycium barbarum", "Lycium ruthenicum",
                             "Nicotiana chinense", "Nicotiana amarum"])
    ap.add_argument("--af-min", type=float, default=0.05, help="二值模式的 ALT 阈值")
    ap.add_argument("--af-bins", default="", help="可选分级阈值, 如 '0.05,0.25,0.75'")
    ap.add_argument("--hamming-threshold", type=float, default=0.15,
                    help="QST 合并的 Hamming 距离上限 (共享位点口径)")
    ap.add_argument("--private-min", type=int, default=2,
                    help="私有 QST 判定的最小成员样本数")
    ap.add_argument("--no-plot", action="store_true", help="跳过热图输出")
    ap.add_argument("--label", default="", help="病毒名标签 (标题/日志)")
    ap.add_argument("--group-col", default="",
                    help="可选: 指定元数据列名作为分组 (默认按物种列表匹配 Host)")
    ap.add_argument("--trait-cols", default="",
                    help="可选: 逗号分隔的数值型性状列 (按 QST 画箱线图 + KW/两两MWU)")
    ap.add_argument("--network-min-count", type=int, default=1,
                    help="网络图节点最小携带样本数")
    ns = ap.parse_args()

    os.makedirs(ns.outdir, exist_ok=True)
    tag = ns.label or os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(ns.matrix))))
    print(f"\n{'=' * 64}\n QST 单倍型形式化 | {tag}\n{'=' * 64}")

    try:
        M, layout = load_matrix(ns.matrix)
        samples = [str(s) for s in M.index]
        site_ids = list(M.columns)
        print(f"[矩阵] {len(samples)} 样本 × {len(site_ids)} 位点 ({layout})")

        species_of = {}
        groups = {}
        unmapped = []
        if ns.metadata:
            if ns.group_col:
                gmap = load_group_map(ns.metadata, ns.run_col, ns.group_col)
                for s in samples:
                    g = gmap.get(s.strip())
                    if g is None:
                        unmapped.append(s)
                    else:
                        species_of[s] = g
                        groups.setdefault(g, []).append(s)
                print(f"[分组] {len(species_of)}/{len(samples)} 样本按 '{ns.group_col}' 分组"
                      + (f"; 未映射 {len(unmapped)} 个" if unmapped else ""))
            else:
                run_map = load_species_map(ns.metadata, ns.run_col, ns.host_col, ns.species)
                for s in samples:
                    sp = run_map.get(s.strip())
                    if sp is None:
                        unmapped.append(s)
                    else:
                        species_of[s] = sp
                        groups.setdefault(sp, []).append(s)
                print(f"[元数据] {len(species_of)}/{len(samples)} 样本映射到物种"
                      + (f"; 未映射 {len(unmapped)} 个" if unmapped else ""))
        else:
            print("[元数据] 未提供, 跳过物种层统计")

        tokens, legend = tokenize_matrix(M, ns.af_min, ns.af_bins)
        print(f"[token] {legend}")

        # ---- 唯一模式 → 聚类 ----
        pat_index = {}
        sample_pat = []
        for row in tokens:
            key = tuple(row)
            if key not in pat_index:
                pat_index[key] = len(pat_index)
            sample_pat.append(pat_index[key])
        unique_pats = list(pat_index.keys())
        print(f"[模式] 唯一 token 模式 {len(unique_pats)} 条, 开始 Hamming<= "
              f"{ns.hamming_threshold:g} 连通聚类 ...")

        # 全 N 模式 (无任何变异信息) 不参与合并也不单独成簇统计
        pat_cluster = cluster_patterns(unique_pats, ns.hamming_threshold)
        cl_patterns, cl_samples = {}, {}
        for pi_idx, cl in enumerate(pat_cluster):
            cl_patterns.setdefault(cl, []).append(pi_idx)
        for i, s in enumerate(samples):
            cl_samples.setdefault(pat_cluster[sample_pat[i]], []).append(s)
        # 按携带样本数降序命名 (EasyHap 语义): 样本最多的簇为 QST001
        ranked = sorted(cl_samples.items(), key=lambda kv: (-len(kv[1]), kv[0]))
        cl_name = {cl: f"QST{k + 1:03d}" for k, (cl, _) in enumerate(ranked)}

        # ---- 样本分配 ----
        assign = pd.DataFrame({
            "Sample": samples,
            "QST_ID": [cl_name[pat_cluster[sample_pat[i]]] for i in range(len(samples))],
            "Species": [species_of.get(s, "") for s in samples],
        })

        # ---- QST 定义表 ----
        defs = []
        for cl, members_samples in ranked:
            example = "|".join(unique_pats[cl_patterns[cl][0]])
            defs.append({
                "QST_ID": cl_name[cl],
                "N_Subpatterns": len(cl_patterns[cl]),
                "Pattern_Example": example,
                "N_Samples": len(members_samples),
                "Member_Samples": ";".join(sorted(members_samples)),
            })
        pd.DataFrame(defs).to_csv(
            os.path.join(ns.outdir, "QST_Definitions.tsv"), sep="\t", index=False)

        assign.to_csv(os.path.join(ns.outdir, "Sample_QST_Assignment.tsv"),
                      sep="\t", index=False)

        # ---- 位点目录 (#4 token 人读表) ----
        pos_list = [parse_site_position(s) for s in site_ids]
        carrier_counts = [
            sum(1 for i in range(len(samples)) if tokens[i][j] not in ("0", "N"))
            for j in range(len(site_ids))]
        catalog = pd.DataFrame({
            "Site_ID": site_ids, "POS": pos_list,
            "N_Samples_Carrying_Allele": carrier_counts,
        }).sort_values("POS")
        catalog.to_csv(os.path.join(ns.outdir, "Variant_Token_Catalog.tsv"),
                       sep="\t", index=False)

        # ---- 物种层统计 ----
        div_rows, priv_rows = [], []
        if groups:
            ct = pd.crosstab(assign["QST_ID"], assign["Species"])
            ordered_cols = [g for g in sorted(groups) if g in ct.columns]
            other_cols = [c for c in ct.columns if c not in ordered_cols]
            ct = ct[ordered_cols + other_cols]
            ct["Total"] = ct.sum(axis=1)
            ct.to_csv(os.path.join(ns.outdir, "QST_by_Host_Species.tsv"), sep="\t")

            for sp in sorted(groups):
                sub = assign[assign["Species"] == sp]
                n = len(sub)
                freqs = sub["QST_ID"].value_counts(normalize=True)
                hd = (n / (n - 1) * (1 - float((freqs ** 2).sum()))
                      if n >= 2 else float("nan"))
                dom_q = freqs.index[0]
                priv = [q for q in freqs.index
                        if set(assign.loc[assign["QST_ID"] == q, "Species"].dropna()) == {sp}
                        and int((assign["QST_ID"] == q).sum()) >= ns.private_min]
                div_rows.append({
                    "Species": sp, "N_Samples": n, "N_QST_Present": len(freqs),
                    "Haplotype_Diversity_Hd": round(hd, 4) if n >= 2 else "NA",
                    "Private_QST_Count": len(priv),
                    "Dominant_QST": dom_q,
                    "Dominant_Share_pct": round(freqs.iloc[0] * 100, 1),
                })
                for q in priv:
                    priv_rows.append({"Species": sp, "QST_ID": q,
                                      "N_Members": int((assign["QST_ID"] == q).sum())})
            pd.DataFrame(div_rows).to_csv(
                os.path.join(ns.outdir, "QST_Diversity_Stats.tsv"), sep="\t", index=False)
            pd.DataFrame(priv_rows, columns=["Species", "QST_ID", "N_Members"]).to_csv(
                os.path.join(ns.outdir, "Private_QST_List.tsv"), sep="\t", index=False)

        top_name, top_members = ranked[0]
        print(f"[结果] {len(ranked)} 个 QST 簇; 最大簇 {cl_name[top_name]}"
              f"(n={len(top_members)}); 唯一模式 {len(unique_pats)} 条")
        if groups and div_rows:
            for r in div_rows:
                print(f"   - {r['Species']}: n={r['N_Samples']}, "
                      f"QST={r['N_QST_Present']}, Hd={r['Haplotype_Diversity_Hd']}, "
                      f"私有={r['Private_QST_Count']}")

        # ---- 热图 ----
        if not ns.no_plot:
            draw_heatmap(tokens, samples, assign, pos_list,
                         cl_name, pat_cluster, sample_pat, tag,
                         ns.outdir, legend)
            print(f"[图]   {os.path.join(ns.outdir, 'qst_token_heatmap.png')}")
            # 分组图 (EasyHap 思路): 物种/分组 × QST 构成
            if groups:
                draw_group_stacked(assign, ns.outdir, tag)
                draw_group_pies(assign, ns.outdir, tag)
            # MST 单倍型网络图
            draw_haplotype_network(samples, tokens, sample_pat, unique_pats,
                                   species_of, ns.network_min_count,
                                   ns.outdir, tag)
            # 性状箱线图 (QST × trait)
            if ns.trait_cols and ns.metadata:
                _draw_traits(ns, assign, tag)
    except InputError as e:
        sys.exit(f"[输入错误] {e}")


def draw_heatmap(tokens, samples, assign, positions, cl_name, pat_cluster,
                 sample_pat, tag, outdir, legend):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.colors import ListedColormap

    size_of_qst = assign["QST_ID"].value_counts().to_dict()
    order = sorted(range(len(samples)),
                   key=lambda i: (-size_of_qst.get(assign.iloc[i]["QST_ID"], 0),
                                  assign.iloc[i]["QST_ID"], samples[i]))
    site_order = sorted(range(len(positions)), key=lambda j: positions[j])

    code_map = {"0": 0, "L": 1, "M": 2, "H": 3, "1": 4, "N": 5}
    mat = np.array([[code_map[tokens[i][j]] for j in site_order] for i in order])
    color_list = [TOKEN_COLORS[c] for c in ["0", "L", "M", "H", "1", "N"]]
    cmap = ListedColormap(color_list)

    fig, ax = plt.subplots(figsize=(min(16, max(8, len(site_order) * 0.18)),
                                    min(14, max(4, len(samples) * 0.22))))
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=5,
              interpolation="nearest")

    # 动态字号 (EasyHap cell_text 思路): 格子物理尺寸可读才画 token 文字
    cell_pt = min(fig.get_figwidth() / max(1, len(site_order)),
                  fig.get_figheight() / max(1, len(samples))) * 72.0
    fs = min(11.5, max(2.0, 0.55 * cell_pt))
    if fs >= 4.5 and len(samples) * len(site_order) <= 6000:
        for i, si in enumerate(order):
            for j, sj in enumerate(site_order):
                ax.text(j, i, tokens[si][sj], ha="center", va="center",
                        fontsize=fs, color="#333333", zorder=5)

    # QST 簇边界线
    qs_seq = [assign.iloc[i]["QST_ID"] for i in order]
    for k in range(1, len(qs_seq)):
        if qs_seq[k] != qs_seq[k - 1]:
            ax.axhline(k - 0.5, color="black", lw=0.9)

    ax.set_xticks(range(len(site_order)))
    ax.set_xticklabels([str(positions[j]) for j in site_order],
                       fontsize=6, rotation=90)
    ax.set_ylabel("Samples (grouped by QST cluster)", fontsize=10)
    ax.set_xlabel("Genomic position (bp)", fontsize=10)
    ax.set_title(f"QST token map — {tag}\n({legend})",
                 fontsize=11, fontweight="bold")
    handles = [Patch(facecolor=TOKEN_COLORS[t], edgecolor="black", linewidth=0.4, label=t)
               for t in ["0", "L", "M", "H", "1", "N"]]
    ax.legend(handles=handles, title="Token", loc="center left",
              bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "qst_token_heatmap.png"),
                dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(os.path.join(outdir, "qst_token_heatmap.pdf"),
                facecolor="white", bbox_inches="tight")
    plt.close(fig)


# ==========================================
# 绘图增强 (EasyHap plotting 思路重写, 零 GPL)
# ==========================================
def _hap_palette(n):
    """Okabe-Ito 循环配色。"""
    base = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9",
            "#CC79A7", "#F0E442", "#999999", "#000000"]
    return [base[i % len(base)] for i in range(n)]


def _sig_stars(q):
    if q < 1e-4:
        return "****"
    if q < 1e-3:
        return "***"
    if q < 1e-2:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def draw_group_stacked(assign, outdir, tag):
    """QST × 分组 构成比例堆叠 bar (EasyHap plot_group_distribution 思路)。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ct = pd.crosstab(assign["QST_ID"], assign["Species"])
    ct = ct.loc[:, [c for c in ct.columns if str(c).strip() != ""]]
    if ct.empty or ct.shape[1] == 0:
        print("[分组bar] 无分组数据, 跳过")
        return
    prop = ct.T.div(ct.T.sum(axis=1), axis=0)  # 分组行 × QST 列
    colors = _hap_palette(len(prop.columns))
    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(prop)), 5.2))
    bottom = np.zeros(len(prop))
    x = np.arange(len(prop))
    for k, q in enumerate(prop.columns):
        ax.bar(x, prop[q].values, bottom=bottom, label=q, color=colors[k],
               edgecolor="black", linewidth=0.3)
        bottom += prop[q].values
    n_per = ct.T.sum(axis=1).astype(int)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\nn={int(n_per.loc[g])}" for g in prop.index],
                       fontsize=9, rotation=15 if len(prop) > 3 else 0)
    ax.set_ylim(0, 1)
    ax.set_ylabel("QST frequency", fontsize=10)
    ax.set_xlabel("Group", fontsize=10)
    ax.set_title(f"QST composition by group — {tag}", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "qst_group_stacked.png"), dpi=300,
                facecolor="white", bbox_inches="tight")
    fig.savefig(os.path.join(outdir, "qst_group_stacked.pdf"),
                facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"[分组bar] {os.path.join(outdir, 'qst_group_stacked.png')}")


def draw_group_pies(assign, outdir, tag):
    """QST × 分组 pie 网格 (每分组一子图, EasyHap plot_group_pie_chart 思路)。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    ct = pd.crosstab(assign["QST_ID"], assign["Species"])
    ct = ct.loc[:, [c for c in ct.columns if str(c).strip() != ""]]
    if ct.empty or ct.shape[1] == 0:
        print("[分组pie] 无分组数据, 跳过")
        return
    colors = _hap_palette(len(ct.index))
    n = ct.shape[1]
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 3.2 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for idx, sp in enumerate(ct.columns):
        ax = axes[idx]
        s = ct[sp]
        s = s[s > 0]
        ax.pie(s.values, colors=[colors[ct.index.get_loc(q)] for q in s.index],
               autopct=lambda p: f"{p:.0f}%" if p >= 4 else "",
               startangle=90, counterclock=False,
               wedgeprops={"edgecolor": "white", "linewidth": 0.8},
               textprops={"fontsize": 7})
        ax.set_title(f"{sp} (n={int(s.sum())})", fontsize=10)
    for idx in range(n, len(axes)):
        axes[idx].axis("off")
    handles = [Patch(facecolor=colors[k], edgecolor="none", label=q)
               for k, q in enumerate(ct.index)]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(0.99, 0.5),
               frameon=False, title="QST", fontsize=8)
    fig.suptitle(f"QST frequency by group — {tag}", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 0.86, 0.94])
    fig.savefig(os.path.join(outdir, "qst_group_pies.png"), dpi=300,
                facecolor="white", bbox_inches="tight")
    fig.savefig(os.path.join(outdir, "qst_group_pies.pdf"),
                facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"[分组pie] {os.path.join(outdir, 'qst_group_pies.png')}")


def draw_trait_boxplots(assign, trait_map, trait, outdir, tag):
    """QST × 性状 箱线图 + Kruskal-Wallis + 两两 Mann-Whitney U (BH 校正)。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import kruskal, mannwhitneyu
    rows = []
    for s in assign["Sample"]:
        v = trait_map.get(str(s))
        if v is not None:
            rows.append({"QST_ID": assign.loc[assign["Sample"] == s, "QST_ID"].iloc[0],
                         "value": float(v)})
    tmp = pd.DataFrame(rows).dropna()
    order = [q for q, n in tmp["QST_ID"].value_counts().items() if n >= 2]
    if len(order) < 2:
        print(f"[性状箱线] {trait}: 有效组 < 2, 跳过")
        return
    data = [tmp.loc[tmp["QST_ID"] == q, "value"].values for q in order]
    fig, ax = plt.subplots(figsize=(max(5.4, 0.6 * len(order) + 4), 5.2))
    bp = ax.boxplot(data, tick_labels=[str(q) for q in order], showfliers=False,
                    patch_artist=True, widths=0.58,
                    medianprops={"color": "black", "lw": 1.5})
    colors = _hap_palette(len(order))
    for b, c in zip(bp["boxes"], colors):
        b.set_facecolor(c)
        b.set_alpha(0.72)
        b.set_edgecolor("0.25")
    for i, vals in enumerate(data, 1):
        ax.text(i, 0.02, f"n={len(vals)}", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=8)
    try:
        overall = kruskal(*data).pvalue
    except Exception:
        overall = float("nan")
    pairs = []
    for a in range(len(order)):
        for b in range(a + 1, len(order)):
            try:
                p = mannwhitneyu(data[a], data[b], alternative="two-sided").pvalue
            except Exception:
                p = float("nan")
            if np.isfinite(p):
                pairs.append([a, b, p])
    m = len(pairs)
    if m:
        ranks = sorted(range(m), key=lambda k: pairs[k][2])
        running = 1.0
        qv = [0.0] * m
        for rk in range(m, 0, -1):
            k = ranks[rk - 1]
            running = min(running, pairs[k][2] * m / rk)
            qv[k] = running
        for k in range(m):
            pairs[k].append(qv[k])
        sig = sorted([p for p in pairs if p[3] < 0.05], key=lambda r: r[3])[:8]
        if sig:
            vals_all = np.concatenate(data)
            ymin, ymax = float(np.nanmin(vals_all)), float(np.nanmax(vals_all))
            span = ymax - ymin or max(abs(ymax), 1.0) * 0.2
            base = ymax + 0.08 * span
            step = 0.10 * span
            for kk, (a, b, p, q) in enumerate(sig):
                y = base + kk * step
                ax.plot([a + 1, a + 1, b + 1, b + 1],
                        [y, y + 0.025 * span, y + 0.025 * span, y],
                        lw=0.8, color="black", clip_on=False)
                ax.text((a + b + 2) / 2, y + 0.025 * span, _sig_stars(q),
                        ha="center", va="bottom", fontsize=9)
            ax.set_ylim(ymin - 0.08 * span, base + (len(sig) + 1.2) * step)
    ax.set_title(f"{trait} by QST — {tag}\n"
                 + (f"Kruskal-Wallis P = {overall:.2e}" if np.isfinite(overall) else ""),
                 fontsize=10, fontweight="bold")
    ax.set_xlabel("QST_ID", fontsize=10)
    ax.set_ylabel(trait, fontsize=10)
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", lw=0.4, alpha=0.25)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, f"qst_{trait}_trait_boxplot.png"), dpi=300,
                facecolor="white", bbox_inches="tight")
    fig.savefig(os.path.join(outdir, f"qst_{trait}_trait_boxplot.pdf"),
                facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"[性状箱线] {os.path.join(outdir, f'qst_{trait}_trait_boxplot.png')}")


def _draw_traits(ns, assign, tag):
    """读取元数据中 --trait-cols 的数值列, 逐个画 QST × 性状箱线图。"""
    import pandas as _pd
    sep = "," if ns.metadata.lower().endswith(".csv") else "\t"
    meta_all = _pd.read_csv(ns.metadata, sep=sep, dtype=str)
    for tc in [x.strip() for x in ns.trait_cols.split(",") if x.strip()]:
        if tc not in meta_all.columns:
            print(f"[性状箱线] 元数据无列 '{tc}', 跳过")
            continue
        trait_map = {}
        for _, r in meta_all.iterrows():
            run = str(r[ns.run_col]).strip()
            v = r[tc]
            if not run or not v or str(v).lower() == "nan":
                continue
            try:
                trait_map[run] = float(v)
            except ValueError:
                pass  # 非数值列跳过
        if trait_map:
            draw_trait_boxplots(assign, trait_map, tc, ns.outdir, tag)


def draw_haplotype_network(samples, tokens, sample_pat, unique_pats, species_of,
                           min_count, outdir, tag):
    """MST 单倍型网络图 (EasyHap plot_haplotype_network 思路重写)。
    节点=唯一 token 模式, 边=共享位点 Hamming 距离, Kruskal 最小生成树;
    节点 pie=宿主分组构成, 大小=携带样本数。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Wedge, Circle
    counts, by_species = {}, {}
    for i, s in enumerate(samples):
        p = sample_pat[i]
        counts[p] = counts.get(p, 0) + 1
        sp = species_of.get(s, "NA")
        by_species.setdefault(p, {}).setdefault(sp, 0)
        by_species[p][sp] += 1
    nodes = [p for p in range(len(unique_pats))
             if counts.get(p, 0) >= max(1, int(min_count))]
    if len(nodes) < 1:
        print("[网络图] 无可绘制模式, 跳过")
        return

    def hdist(i, j):
        t1, t2 = unique_pats[nodes[i]], unique_pats[nodes[j]]
        shared = mism = 0
        for a, b in zip(t1, t2):
            if a != "N" and b != "N":
                shared += 1
                if a != b:
                    mism += 1
        return (mism if shared else 0, shared)

    edges = []
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            d, sh = hdist(i, j)
            if sh > 0:
                edges.append((i, j, d))
    mst = []
    if edges:
        dsu = DSU(len(nodes))
        for a, b, d in sorted(edges, key=lambda e: e[2]):
            if dsu.find(a) != dsu.find(b):
                dsu.union(a, b)
                mst.append((a, b, d))
    try:
        import networkx as nx
        G = nx.Graph()
        for k in range(len(nodes)):
            G.add_node(k)
        for a, b, d in mst:
            G.add_edge(a, b, distance=d)
        pos = (nx.kamada_kawai_layout(G, weight="distance") if len(nodes) > 1
               else {0: np.array([0.0, 0.0])})
        layout = "kamada_kawai"
    except Exception:
        pos = {k: np.array([np.cos(2 * np.pi * k / len(nodes)),
                            np.sin(2 * np.pi * k / len(nodes))])
               for k in range(len(nodes))}
        layout = "circular(fallback)"
    species = sorted({sp for p in nodes for sp in by_species.get(p, {})})
    colors = _hap_palette(len(species))
    sp_color = {sp: colors[i] for i, sp in enumerate(species)}
    fig, ax = plt.subplots(figsize=(9.0, 8.0))
    xs = [float(pos[k][0]) for k in range(len(nodes))]
    ys = [float(pos[k][1]) for k in range(len(nodes))]
    span = max(max(xs) - min(xs) if len(xs) > 1 else 1.0,
               max(ys) - min(ys) if len(ys) > 1 else 1.0, 1.0)
    for a, b, d in mst:
        x1, y1 = pos[a]
        x2, y2 = pos[b]
        ax.plot([x1, x2], [y1, y2], color="0.35", lw=1.2, zorder=1)
        if d > 0:
            ax.text((x1 + x2) / 2, (y1 + y2) / 2, str(d), fontsize=7,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6},
                    zorder=4)
    max_count = max(counts[p] for p in nodes)
    for k, p in enumerate(nodes):
        x, y = pos[k]
        total = max(1, counts[p])
        radius = 0.035 * span * (0.8 + 0.7 * (total / max_count) ** 0.5)
        start = 90.0
        for sp in species:
            c = by_species.get(p, {}).get(sp, 0)
            if c <= 0:
                continue
            t2 = start + 360.0 * c / total
            ax.add_patch(Wedge((x, y), radius, start, t2, facecolor=sp_color[sp],
                               edgecolor="white", lw=0.6, zorder=5))
            start = t2
        ax.add_patch(Circle((x, y), radius, fill=False, edgecolor="black", lw=0.8, zorder=6))
        ax.text(x, y + radius * 1.35, f"P{nodes[k] + 1:03d} (n={total})",
                ha="center", va="bottom", fontsize=7.5, zorder=7)
    pad = 0.22 * span
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    handles = [Patch(facecolor=sp_color[sp], edgecolor="none", label=sp)
               for sp in species]
    if handles:
        ax.legend(handles=handles, title="Host species", frameon=False,
                  bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    ax.set_title(f"Minimum-spanning haplotype network — {tag}\n"
                 f"edge labels = mutational distance (layout: {layout})",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "qst_haplotype_network.png"), dpi=300,
                facecolor="white", bbox_inches="tight")
    fig.savefig(os.path.join(outdir, "qst_haplotype_network.pdf"),
                facecolor="white", bbox_inches="tight")
    plt.close(fig)
    with open(os.path.join(outdir, "QST_Network_Edges.tsv"), "w") as f:
        f.write("Pattern1\tPattern2\tMutationalDistance\n")
        for a, b, d in sorted(mst):
            f.write(f"P{nodes[a] + 1:03d}\tP{nodes[b] + 1:03d}\t{d}\n")
    print(f"[网络图] {os.path.join(outdir, 'qst_haplotype_network.png')} ({layout})")


if __name__ == "__main__":
    main()
