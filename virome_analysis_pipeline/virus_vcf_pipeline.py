#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
病毒VCF合并、质控、下游群体遗传学分析管线  v3.0
=================================================
merge → QC/过滤 → SNP矩阵导出 → 自定义距离 → NJ树 → 可视化
所有新增模块均可通过 CLI flag 独立开启, 默认行为与 v1 兼容。
"""

import os
import sys
import glob
import subprocess
import argparse
import shutil
from collections import defaultdict

# ---- 科学计算 / 可视化 (优雅降级) ----
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_VIZ = True
except ImportError:
    HAS_VIZ = False

try:
    from sklearn.decomposition import PCA
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


# ============================================================
# 工具函数
# ============================================================

def check_dependencies(extra=None):
    """检查命令行工具; extra 为可选 Python 包名列表"""
    tools = ["bcftools", "bgzip", "tabix"]
    if extra is None:
        extra = []
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        print(f"[错误] 找不到工具: {', '.join(missing)}")
        sys.exit(1)
    pkg_map = {"numpy": HAS_NUMPY, "scipy": HAS_SCIPY, "viz": HAS_VIZ, "sklearn": HAS_SKLEARN}
    for pkg in extra:
        if pkg in pkg_map and not pkg_map[pkg]:
            print(f"[错误] 缺少 Python 包 '{pkg}' (pip install {pkg})")
            sys.exit(1)


def run_command(cmd, return_output=False, check=True):
    """执行 shell 命令。return_output=True 时返回 stdout 字符串。"""
    try:
        if return_output:
            return subprocess.check_output(cmd, shell=True, text=True)
        else:
            subprocess.run(cmd, shell=True, check=check)
    except subprocess.CalledProcessError as e:
        print(f"[警告] 命令返回非零: {cmd}")
        if return_output:
            return ""


# ============================================================
# 1. QC: 样本级统计 & 异常过滤
# ============================================================

def qc_sample_stats(merged_vcf):
    """
    从 merged VCF 提取每个样本的:
      - n_variants: 非缺失基因型数
      - missing_rate: 缺失率
      - ts_count / tv_count / ts_tv_ratio
    返回 dict: {sample_name: {...}}
    """
    # 样本列表
    samples_str = run_command(f"bcftools query -l '{merged_vcf}'", return_output=True)
    samples = [s.strip() for s in samples_str.strip().split("\n") if s.strip()]
    if not samples:
        return {}

    # 每个位点的非缺失样本信息: CHROM POS REF ALT GT_ARRAY
    raw = run_command(
        f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' '{merged_vcf}'",
        return_output=True
    )

    n_sites = 0
    sample_variants = defaultdict(int)
    sample_missing = defaultdict(int)
    sample_ts = defaultdict(int)
    sample_tv = defaultdict(int)
    transitions = [{"A", "G"}, {"G", "A"}, {"C", "T"}, {"T", "C"}]

    for line in raw.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 5 + len(samples):
            continue
        ref, alt = parts[2].upper(), parts[3].upper()
        is_ts = {ref, alt} in transitions
        gts = parts[4:]  # one per sample

        n_sites += 1
        for i, gt_raw in enumerate(gts):
            # iVar VCF 的 GT 可能是 "./." 或 "0" 或 "1" 或 "0/1" 等
            if gt_raw in ("./.", ".", "./.", ""):
                sample_missing[samples[i]] += 1
            else:
                # 只要不是缺失就算 variant
                sample_variants[samples[i]] += 1
                if is_ts:
                    sample_ts[samples[i]] += 1
                else:
                    sample_tv[samples[i]] += 1

    stats = {}
    for s in samples:
        n_var = sample_variants.get(s, 0)
        n_mis = sample_missing.get(s, 0)
        total = n_var + n_mis
        ts_c = sample_ts.get(s, 0)
        tv_c = sample_tv.get(s, 0)
        stats[s] = {
            "n_variants": n_var,
            "n_missing": n_mis,
            "missing_rate": round(n_mis / total, 4) if total > 0 else 1.0,
            "ts_count": ts_c,
            "tv_count": tv_c,
            "ts_tv_ratio": round(ts_c / tv_c, 4) if tv_c > 0 else float("inf"),
        }
    return stats


def qc_filter_samples(merged_vcf, out_vcf, stats, max_snps=0, max_missing=0.2):
    """
    根据阈值过滤样本, 输出新的 VCF。
    max_snps > 0: 过滤超过该值的样本
    max_missing > 0: 过滤缺失率超过该值的样本
    返回 (保留的样本列表, 被排除的样本列表)
    """
    excluded = set()
    for s, st in stats.items():
        if max_snps > 0 and st["n_variants"] > max_snps:
            excluded.add(s)
        if max_missing > 0 and st["missing_rate"] > max_missing:
            excluded.add(s)

    kept = [s for s in stats if s not in excluded]

    if excluded:
        print(f"[QC] 过滤掉 {len(excluded)} 个异常样本: {', '.join(sorted(excluded))}")
        # 用 bcftools view 保留子集
        kept_str = ",".join(kept)
        run_command(f"bcftools view -s '{kept_str}' '{merged_vcf}' -O z -o '{out_vcf}'")
        run_command(f"tabix -p vcf '{out_vcf}'")
    else:
        print(f"[QC] 所有 {len(kept)} 个样本通过过滤。")
        if out_vcf != merged_vcf:
            shutil.copy(merged_vcf, out_vcf)
            if not os.path.exists(out_vcf + ".tbi"):
                run_command(f"tabix -p vcf '{out_vcf}'")

    return kept, excluded


def qc_summary_table(stats, excluded, out_path):
    """输出 QC 统计表 TSV"""
    print("\n[QC 提示]: bcftools merge 产生的 ./. 多见于病毒群体 (异质性极高),")
    print("  通常代表该位点为野生型或覆盖不足, 而非测序失败。")
    print("  建议: 优先用 --max-snps 过滤超突变样本, missing_rate 仅作参考。\n")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("Sample\tn_variants\tn_missing\tmissing_rate\tts_count\ttv_count\tts_tv_ratio\tstatus\n")
        for s, st in sorted(stats.items()):
            status = "EXCLUDED" if s in excluded else "PASS"
            f.write(f"{s}\t{st['n_variants']}\t{st['n_missing']}\t{st['missing_rate']}\t"
                    f"{st['ts_count']}\t{st['tv_count']}\t{st['ts_tv_ratio']}\t{status}\n")
    print(f"[QC] 汇总表已保存: {out_path}")


# ============================================================
# 2. SNP 矩阵导出 (样本 × 位点, 0/1/NA)
# ============================================================

def export_snp_matrix(vcf_path, out_tsv):
    """
    从 VCF 提取基因型矩阵: 样本(行) × 位点(列), 值 0/1/NaN。
    使用 bcftools query + numpy 后处理, 无需 cyvcf2。
    """
    assert HAS_NUMPY, "需要 numpy"

    samples_str = run_command(f"bcftools query -l '{vcf_path}'", return_output=True)
    samples = [s.strip() for s in samples_str.strip().split("\n") if s.strip()]
    if not samples:
        print("[SNP矩阵] 无样本, 跳过。")
        return None, None

    # 获取位点信息: CHROM POS REF ALT
    positions = run_command(
        f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' '{vcf_path}'", return_output=True
    )
    pos_info = [p.strip().split("\t") for p in positions.strip().split("\n") if p.strip()]
    pos_labels = [f"{p[0]}_{p[1]}_{p[2]}_{p[3]}" for p in pos_info]

    # 获取 GT 矩阵
    raw_gt = run_command(
        f"bcftools query -f '[\\t%GT]\\n' '{vcf_path}'", return_output=True
    )

    n_sites = len(pos_labels)
    n_samples = len(samples)
    M = np.full((n_samples, n_sites), np.nan, dtype=np.float64)

    for i, line in enumerate(raw_gt.strip().split("\n")):
        if i >= n_sites:
            break
        parts = line.strip().split("\t")
        for j, gt_str in enumerate(parts):
            if j >= n_samples:
                break
            gt_str = gt_str.strip()
            if gt_str in (".", "./.", ".|."):
                continue  # keep NaN
            # 支持 0, 1, 0/1, 1/1, 0|1 等格式
            alleles = set(gt_str.replace("/", "|").split("|"))
            alleles.discard(".")
            if not alleles:
                continue
            # haploid iVar: GT 就是 "1" 代表 alt; 任意非0等位基因视为 alt
            M[j, i] = 0 if alleles == {"0"} else 1

    # 保存 TSV
    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("Sample\t" + "\t".join(pos_labels) + "\n")
        for j, s_name in enumerate(samples):
            row_str = "\t".join(
                str(int(M[j, k])) if not np.isnan(M[j, k]) else "NA"
                for k in range(n_sites)
            )
            f.write(f"{s_name}\t{row_str}\n")

    print(f"[SNP矩阵] {n_samples} x {n_sites} 已导出: {out_tsv}")
    return M, samples, pos_labels


# ============================================================
# 3. 距离矩阵计算 (Jaccard / Hamming — 保留用于 0/1 矩阵)
# ============================================================

def pairwise_distance_matrix(M, metric="jaccard"):
    """
    M: samples × sites (0/1/NaN numpy array)
    metric: "jaccard" | "hamming" | "both"
    返回 dict: {"jaccard": ndarray, "hamming": ndarray}
    """
    assert HAS_NUMPY, "需要 numpy"
    n = M.shape[0]

    D_jacc = np.zeros((n, n))
    D_hamm = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            a, b = M[i], M[j]
            mask = ~np.isnan(a) & ~np.isnan(b)
            shared = mask.sum()
            if shared == 0:
                D_jacc[i, j] = D_jacc[j, i] = np.nan
                D_hamm[i, j] = D_hamm[j, i] = np.nan
                continue

            ai, bi = a[mask].astype(int), b[mask].astype(int)

            # Jaccard distance: 1 - |intersection| / |union|
            inter = np.sum((ai == 1) & (bi == 1))
            union = np.sum((ai == 1) | (bi == 1))
            d_j = 1.0 - (inter / union) if union > 0 else 0.0  # 两样本全 ref → 距离 0
            D_jacc[i, j] = D_jacc[j, i] = d_j

            # Hamming: fraction of mismatching sites
            d_h = np.sum(ai != bi) / shared
            D_hamm[i, j] = D_hamm[j, i] = d_h

    results = {}
    if metric in ("jaccard", "both"):
        results["jaccard"] = D_jacc
    if metric in ("hamming", "both"):
        results["hamming"] = D_hamm
    return results


def save_distance_matrix(D, samples, out_path):
    """保存距离矩阵为 TSV (含行列名)"""
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("Sample\t" + "\t".join(samples) + "\n")
        for i, s in enumerate(samples):
            row = "\t".join(f"{D[i, j]:.6f}" if not np.isnan(D[i, j]) else "NA"
                            for j in range(len(samples)))
            f.write(f"{s}\t{row}\n")
    print(f"[距离矩阵] 已保存: {out_path}")


# ============================================================
# 4. NJ 树 (scipy linkage → Newick)
# ============================================================

def linkage_to_newick(Z, labels):
    """
    将 scipy linkage matrix 转为 Newick 字符串。
    Z: scipy.cluster.hierarchy.linkage 输出
    labels: 叶节点名称列表
    """
    n = len(labels)
    nodes = {i: labels[i] for i in range(n)}

    for k, (c1, c2, dist, _) in enumerate(Z):
        c1, c2 = int(c1), int(c2)
        left = nodes[c1]
        right = nodes[c2]
        # 若子节点为内部节点, 加括号
        if c1 >= n:
            left = f"({left})"
        if c2 >= n:
            right = f"({right})"
        bl = f":{dist:.6f}"
        # length already stored in dist for this join
        nodes[n + k] = f"({left}{bl},{right}{bl})"

    return f"({nodes[2 * n - 2]});"


def build_nj_tree(D, samples, out_path, method="average"):
    """
    从距离矩阵构建 NJ 树 (通过 UPGMA/WPGMA 近似)。
    method: "average" | "ward" | "single" | "complete"
    """
    assert HAS_SCIPY, "需要 scipy"

    # NaN 处理: 无共有位点的样本对 → 赋予群体最大距离 (而非 0)
    max_dist = np.nanmax(D) if not np.isnan(D).all() else 1.0
    D_clean = np.nan_to_num(D, nan=max_dist * 1.1)

    # 对称化
    np.fill_diagonal(D_clean, 0.0)

    condensed = squareform(D_clean, checks=False)
    Z = linkage(condensed, method=method)
    newick = linkage_to_newick(Z, list(samples))

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(newick + "\n")
    print(f"[树] Newick 已保存: {out_path}")
    return Z, newick


# ============================================================
# 5. 可视化
# ============================================================

def compute_ibs_kinship(M):
    """
    IBS Kinship (VCF2PCACluster Method 1: Normalized_IBS / Yang-BaldingNicols).
    M: samples × sites (0/1/NaN).
    返回: (Kinship matrix, 样本数, 样本名)
    """
    n = M.shape[0]
    K = np.zeros((n, n))
    C = np.zeros((n, n))  # 有效位点计数

    for i in range(n):
        for j in range(i, n):
            a, b = M[i], M[j]
            mask = ~np.isnan(a) & ~np.isnan(b)
            shared = mask.sum()
            if shared == 0:
                continue
            ai, bi = a[mask].astype(int), b[mask].astype(int)
            ibs = (ai == bi).sum()  # 相同基因型的位点数
            K[i, j] = ibs
            K[j, i] = ibs
            C[i, j] = shared
            C[j, i] = shared

    # 归一化: Kinship = IBS / count
    valid = C > 0
    K[valid] /= C[valid]
    return K, n


def pca_from_kinship(K, samples, out_path):
    """
    从 IBS Kinship 矩阵做特征分解 → PCA (等价于 VCF2PCACluster 默认方法)。
    """
    assert HAS_SCIPY, "需要 scipy (linalg.eigh)"

    eigvals, eigvecs = np.linalg.eigh(K)
    # 按特征值降序排列
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # 只保留正值特征值对应的 PC
    n_pos = (eigvals > 1e-10).sum()
    n_pcs = max(2, min(n_pos, len(samples) - 1))
    total = eigvals[:n_pcs].sum()

    # 前两个 PC
    pc1 = eigvecs[:, 0] * np.sqrt(eigvals[0]) if eigvals[0] > 0 else eigvecs[:, 0]
    pc2 = eigvecs[:, 1] * np.sqrt(eigvals[1]) if n_pcs > 1 and eigvals[1] > 0 else eigvecs[:, 0]
    r1 = eigvals[0] / total * 100 if total > 0 else 0
    r2 = eigvals[1] / total * 100 if total > 0 and n_pcs > 1 else 0

    # 画图
    plt.figure(figsize=(8, 7))
    plt.scatter(pc1, pc2, s=80, c="#2b83ba", edgecolors="black", alpha=0.85)
    for i, s_name in enumerate(samples):
        plt.annotate(s_name, (pc1[i], pc2[i]), fontsize=7, alpha=0.7,
                     textcoords="offset points", xytext=(4, 4))
    plt.xlabel(f"PC1 ({r1:.1f}%)", fontsize=12, fontweight="bold")
    plt.ylabel(f"PC2 ({r2:.1f}%)", fontsize=12, fontweight="bold")
    plt.title("PCA via IBS Kinship (VCF2PCACluster method)", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[PCA图-IBS] 已保存: {out_path}")


def pca_from_genotypes(M, samples, out_path):
    """
    标准基因型 PCA (仅均值中心化, 不做 StandardScaler)。
    对病毒 haploid 数据, 不做方差缩放, 避免稀有 SNP 权重被放大。
    """
    assert HAS_SKLEARN, "需要 sklearn"

    M_filled = np.where(np.isnan(M), 0.0, M)
    col_std = np.std(M_filled, axis=0)
    M_var = M_filled[:, col_std > 0]
    if M_var.shape[1] < 2:
        print("[PCA] 变异位点太少, 跳过。")
        return

    # 只中心化, 不缩放 (不用 StandardScaler)
    X_centered = M_var - M_var.mean(axis=0)
    pca = PCA(n_components=2)
    X_pc = pca.fit_transform(X_centered)
    var_r = pca.explained_variance_ratio_

    plt.figure(figsize=(8, 7))
    plt.scatter(X_pc[:, 0], X_pc[:, 1], s=80, c="#e74c3c", edgecolors="black", alpha=0.85)
    for i, s_name in enumerate(samples):
        plt.annotate(s_name, (X_pc[i, 0], X_pc[i, 1]), fontsize=7, alpha=0.7,
                     textcoords="offset points", xytext=(4, 4))
    plt.xlabel(f"PC1 ({var_r[0]*100:.1f}%)", fontsize=12, fontweight="bold")
    plt.ylabel(f"PC2 ({var_r[1]*100:.1f}%)", fontsize=12, fontweight="bold")
    plt.title("PCA of SNP Genotypes (center-only)", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[PCA图-Genotype] 已保存: {out_path}")


def plot_pca(M, samples, out_path, method="kinship"):
    """
    PCA 散点图。
    method: "kinship" (VCF2PCACluster 等价) | "genotype" (标准中心化 PCA)
    """
    assert HAS_VIZ, "需要 matplotlib"

    if method == "kinship":
        K, _ = compute_ibs_kinship(M)
        pca_from_kinship(K, samples, out_path)
    else:
        pca_from_genotypes(M, samples, out_path)


def plot_distance_heatmap(D, samples, out_path, title="Sample Distance Heatmap"):
    """距离矩阵热图 + 层次聚类树"""
    assert HAS_VIZ, "需要 matplotlib + seaborn"

    D_df = {samples[i]: {samples[j]: D[i, j] for j in range(len(samples))}
            for i in range(len(samples))}
    import pandas as pd
    df = pd.DataFrame(D_df)

    # 处理 NaN
    mask = df.isna()

    g = sns.clustermap(df, cmap="YlOrRd_r", metric="precomputed",
                       row_cluster=True, col_cluster=True,
                       linewidths=0.5, figsize=(max(10, len(samples)*0.4),
                                                 max(8, len(samples)*0.4)),
                       mask=mask, cbar_kws={"label": "Distance"})
    g.ax_heatmap.set_title(title, fontsize=14, fontweight="bold", pad=20)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[热图] 已保存: {out_path}")


def plot_afs(M, out_path):
    """等位基因频率谱直方图"""
    assert HAS_VIZ, "需要 matplotlib"

    freqs = np.nanmean(M, axis=0)
    freqs = freqs[~np.isnan(freqs)]

    if len(freqs) == 0:
        print("[AFS] 无数据, 跳过。")
        return

    plt.figure(figsize=(8, 5))
    plt.hist(freqs, bins=30, color="#8da0cb", edgecolor="black", alpha=0.8)
    plt.xlabel("Alternate Allele Frequency", fontsize=12, fontweight="bold")
    plt.ylabel("Number of Sites", fontsize=12, fontweight="bold")
    plt.title("Allele Frequency Spectrum (AFS)", fontsize=14, fontweight="bold")

    # 标记稀有变异区
    plt.axvline(x=0.05, color="red", linestyle="--", alpha=0.6, label="Rare (5%)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[AFS图] 已保存: {out_path}")


def plot_dendrogram(Z, samples, out_path, title="Hierarchical Clustering"):
    """层次聚类树状图"""
    assert HAS_VIZ, "需要 matplotlib"

    plt.figure(figsize=(max(8, len(samples)*0.3), 6))
    dendrogram(Z, labels=list(samples), leaf_font_size=9,
               color_threshold=0.7 * max(Z[:, 2]))
    plt.title(title, fontsize=14, fontweight="bold")
    plt.xlabel("Sample", fontsize=12)
    plt.ylabel("Distance", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[树图] 已保存: {out_path}")


# ============================================================
# 5.6  LD (r²) 共突变分析 — 用于寻找上位效应, 不画假 LD decay
# ============================================================

def compute_ld_r2(M, site_labels, min_r2=0.5):
    """
    计算双等位变异的连锁不平衡 (r²), 仅返回 r² >= min_r2 的强相关位点对。
    用于检测 RNA 病毒中的共突变 (上位效应), 而非真核生物的 LD decay。
    M: samples × sites (0/1/NaN)
    """
    assert HAS_NUMPY, "需要 numpy"
    M_sites = M.T  # sites × samples
    n_sites = M_sites.shape[0]

    positions = []
    for label in site_labels:
        try:
            pos = int(label.split("_")[1])
        except (IndexError, ValueError):
            pos = 0
        positions.append(pos)

    strong_pairs = []
    for i in range(n_sites):
        ai = M_sites[i]
        mask_i = ~np.isnan(ai)
        if mask_i.sum() < 3:
            continue
        pi = np.nanmean(ai)
        if pi == 0 or pi == 1:
            continue
        for j in range(i + 1, n_sites):
            bj = M_sites[j]
            mask_j = ~np.isnan(bj)
            mask = mask_i & mask_j
            if mask.sum() < 3:
                continue
            a, b = ai[mask], bj[mask]
            pa, pb = a.mean(), b.mean()
            if pa == 0 or pa == 1 or pb == 0 or pb == 1:
                continue
            cov = ((a - pa) * (b - pb)).mean()
            denom = pa * (1 - pa) * pb * (1 - pb)
            if denom <= 0:
                continue
            r2 = (cov / np.sqrt(denom)) ** 2
            if r2 >= min_r2:
                strong_pairs.append({
                    "Site_1": site_labels[i],
                    "Site_2": site_labels[j],
                    "Distance_bp": abs(positions[j] - positions[i]),
                    "R_squared": round(r2, 4)
                })

    if not strong_pairs:
        print("[共突变] 未找到 r² >= 0.5 的强连锁位点对。")
        return None

    return pd.DataFrame(strong_pairs).sort_values(by="R_squared", ascending=False)


def plot_epistatic_network(df_pairs, site_labels, positions, out_path):
    """上位效应可视化: 网络散点(x=位点1,y=位点2,色=r²) + r²分布直方图"""
    if df_pairs is None or df_pairs.empty:
        return
    assert HAS_VIZ, "需要 matplotlib"
    pos_map = dict(zip(site_labels, positions))
    pts = []
    for _, r in df_pairs.iterrows():
        p1, p2 = pos_map.get(r["Site_1"], 0), pos_map.get(r["Site_2"], 0)
        if p1 and p2:
            pts.append((p1, p2, r["R_squared"]))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    if pts:
        xs, ys, cs, ss = zip(*[(x, y, c, max(3, c * 20)) for x, y, c in pts])
        sc = ax1.scatter(xs, ys, c=cs, s=ss, cmap="YlOrRd",
                         edgecolors="grey", linewidth=0.5, alpha=0.8)
        ax1.plot([min(xs + ys), max(xs + ys)], [min(xs + ys), max(xs + ys)],
                 '--', color='grey', alpha=0.3)
        plt.colorbar(sc, ax=ax1, shrink=0.8).set_label("R²", fontsize=11)
    ax1.set_xlabel("Site 1 Position"); ax1.set_ylabel("Site 2 Position")
    ax1.set_title(f"Epistatic Co-mutation Network ({len(pts)} pairs)", fontweight="bold")

    ax2.hist(df_pairs["R_squared"], bins=20, color="#E74C3C", edgecolor="white")
    ax2.axvline(0.5, color="black", linestyle="--", label="threshold")
    ax2.set_xlabel("R²"); ax2.set_ylabel("Count")
    ax2.set_title("R² Distribution", fontweight="bold"); ax2.legend()

    plt.tight_layout(); plt.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close()
    print(f"[共突变图] 已保存: {out_path}")


# ============================================================
# 5.7  滑动窗 π 与 Tajima's D
# ============================================================

def compute_sliding_window_popgen(M, site_labels, win_size=0, step=0):
    """
    M: samples × sites
    site_labels: site IDs for position parsing
    win_size/step: 窗口大小和步长 (0=自动)
    返回 DataFrame: BIN_START, BIN_END, BIN_MID, N_SITES, PI, TAJIMA_D
    """
    assert HAS_NUMPY, "需要 numpy"
    n_samples = M.shape[0]

    # 解析位置
    positions = []
    for label in site_labels:
        try:
            pos = int(label.split("_")[1])
        except (IndexError, ValueError):
            pos = 0
        positions.append(pos)

    if not positions:
        return None

    pos_arr = np.array(positions)
    max_pos = pos_arr.max()

    # 自动窗口
    if win_size <= 0 or step <= 0:
        if max_pos <= 500:
            win_size, step = 50, 25
        else:
            step = max(50, int(max_pos // 50))
            win_size = step * 2

    # 每个位点的 π_site (haploid pairwise differences)
    M_T = M.T  # sites × samples
    results = []

    for start in range(1, int(max_pos), step):
        end = start + win_size - 1
        mask = (pos_arr >= start) & (pos_arr <= end)
        idx = np.where(mask)[0]
        if len(idx) < 2:
            results.append({
                "BIN_START": start, "BIN_END": end, "BIN_MID": (start + end) / 2,
                "N_SITES": len(idx), "PI": np.nan, "TAJIMA_D": np.nan
            })
            continue

        window_M = M_T[idx]  # sites_in_window × samples

        # π: mean pairwise difference per site
        pi_sum = 0.0
        valid_sites = 0
        for row in window_M:
            valid = ~np.isnan(row)
            n_v = valid.sum()
            if n_v < 2:
                continue
            vals = row[valid].astype(int)
            p = vals.mean()
            # π per biallelic site: 1 - p² - (1-p)² ≡ 2p(1-p)
            # (前提: bcftools view -m2 -M2 已过滤为纯双等位; 多等位时需 1-Σpi²)
            h = 2 * p * (1 - p)
            h_corrected = h * n_v / (n_v - 1) if n_v > 1 else h
            pi_sum += h_corrected
            valid_sites += 1
        pi_window = pi_sum / valid_sites if valid_sites > 0 else 0.0
        S = valid_sites

        # Tajima's D
        if S > 0 and n_samples >= 2:
            a1 = sum(1.0 / i for i in range(1, n_samples))
            a2 = sum(1.0 / (i ** 2) for i in range(1, n_samples))
            b1 = (n_samples + 1.0) / (3.0 * (n_samples - 1.0))
            b2 = 2.0 * (n_samples ** 2 + n_samples + 3.0) / (9.0 * n_samples * (n_samples - 1.0))
            c1 = b1 - 1.0 / a1
            c2 = b2 - (n_samples + 2.0) / (a1 * n_samples) + a2 / (a1 ** 2)
            e1 = c1 / a1
            e2 = c2 / (a1 ** 2 + a2)
            V_d = e1 * S + e2 * S * (S - 1.0)
            taj_d = (pi_window - S / a1) / math.sqrt(V_d) if V_d > 0 else np.nan
        else:
            taj_d = np.nan

        results.append({
            "BIN_START": start, "BIN_END": end, "BIN_MID": (start + end) / 2,
            "N_SITES": len(idx), "PI": pi_window, "TAJIMA_D": taj_d
        })

    return pd.DataFrame(results)


def plot_popgen_windows(df_popgen, out_path):
    """滑动窗 π + Tajima's D 双面板图"""
    assert HAS_VIZ, "需要 matplotlib"

    df = df_popgen.dropna(subset=["PI", "TAJIMA_D"])
    if df.empty:
        print("[PopGen] 无有效窗口数据, 跳过绘图。")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    ax1.plot(df["BIN_MID"], df["PI"], color="#e74c3c", linewidth=1.5)
    ax1.fill_between(df["BIN_MID"], 0, df["PI"], color="#e74c3c", alpha=0.15)
    ax1.set_ylabel("Nucleotide Diversity (π)", fontsize=12, fontweight="bold")
    ax1.set_title("Sliding Window Population Genetics", fontsize=14, fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.4)

    ax2.plot(df["BIN_MID"], df["TAJIMA_D"], color="#2980b9", linewidth=1.5)
    ax2.fill_between(df["BIN_MID"], 0, df["TAJIMA_D"], color="#2980b9", alpha=0.15)
    ax2.axhline(y=0, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax2.set_ylabel("Tajima's D", fontsize=12, fontweight="bold")
    ax2.set_xlabel("Genomic Position (bp)", fontsize=12, fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[PopGen图] 已保存: {out_path}")


# ============================================================
# 5.8  Per-site 注释 (从 per-sample SnpEff 汇总)
# ============================================================

def parse_gtf_gene_spans(gtf_path):
    """从 GTF 文件提取基因坐标区间: {gene_name: (start, end)}"""
    spans = {}
    if not gtf_path or not os.path.exists(gtf_path):
        return spans
    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            start, end = int(parts[3]), int(parts[4])
            # 提取 gene_name
            import re
            m = re.search(r'gene_name "([^"]+)"', parts[8])
            if not m:
                m = re.search(r'gene_id "([^"]+)"', parts[8])
            gname = m.group(1) if m else f"gene_{start}"
            if gname not in spans:
                spans[gname] = [start, end]
            else:
                spans[gname][0] = min(spans[gname][0], start)
                spans[gname][1] = max(spans[gname][1], end)
    return spans


def collect_per_site_annotations(snpeff_dir):
    """
    扫描 per-sample SnpEff annotation_summary.tsv, 按 (POS,REF,ALT) 去重汇总。
    返回 dict: {(pos, ref, alt): {"gene": ..., "effect": ..., "impact": ...}}
    """
    import glob as _glob
    ann_files = _glob.glob(os.path.join(snpeff_dir, "**", "annotation_summary.tsv"), recursive=True)
    if not ann_files:
        ann_files = _glob.glob(os.path.join(snpeff_dir, "**", "*.annotation_summary.tsv"), recursive=True)
    if not ann_files:
        print("[注释] 未找到 annotation_summary.tsv 文件, 跳过。")
        return {}

    site_ann = {}
    for fp in ann_files:
        try:
            with open(fp, "r") as f:
                header = f.readline()
                for line in f:
                    if not line.strip():
                        continue
                    cols = line.strip().split("\t")
                    if len(cols) < 9:
                        continue
                    # 格式: CHROM POS REF ALT GENE EFFECT IMPACT DNA_CHANGE AA_CHANGE
                    chrom, pos, ref, alt = cols[0], cols[1], cols[2], cols[3]
                    gene = cols[4] if cols[4] else "intergenic"
                    effect = cols[5] if cols[5] else "UNKNOWN"
                    impact = cols[6] if len(cols) > 6 else "MODIFIER"
                    key = (pos, ref, alt)
                    if key not in site_ann:
                        site_ann[key] = {"gene": gene, "effect": effect, "impact": impact}
        except Exception:
            pass

    print(f"[注释] 从 {len(ann_files)} 个文件中收集了 {len(site_ann)} 个唯一位点的注释。")
    return site_ann


def build_per_gene_summary(M, site_labels, site_annotations, df_popgen, gene_spans):
    """
    构建 per-gene 汇总表。
    M: samples × sites
    site_labels: ["CHROM_POS_REF_ALT", ...]
    site_annotations: {(pos, ref, alt): {...}}
    df_popgen: sliding window popgen DataFrame (可为 None)
    gene_spans: {gene: (start, end)}
    返回 DataFrame。
    """
    assert HAS_NUMPY, "需要 numpy"

    # 为每个位点分配基因
    site_positions = []
    site_genes = []
    for label in site_labels:
        parts = label.split("_")
        try:
            pos = int(parts[1])
            ref, alt = parts[2], parts[3]
        except (IndexError, ValueError):
            site_genes.append("unknown")
            site_positions.append(0)
            continue
        site_positions.append(pos)
        key = (str(pos), ref, alt)
        ann = site_annotations.get(key, {})
        gene = ann.get("gene", "intergenic")
        # 如果注释没有基因信息，用 GTF spans 反查
        if gene == "intergenic" and gene_spans:
            for gname, (gs, ge) in gene_spans.items():
                if gs <= pos <= ge:
                    gene = gname
                    break
        site_genes.append(gene)

    # 按基因统计
    from collections import defaultdict
    gene_data = defaultdict(lambda: {
        "n_sites": 0, "positions": [], "alt_freqs": [], "top_site": "", "top_freq": 0.0
    })

    for i, gene in enumerate(site_genes):
        gd = gene_data[gene]
        gd["n_sites"] += 1
        gd["positions"].append(site_positions[i])
        # 等位频率
        col = M[:, i]
        valid = col[~np.isnan(col)]
        af = valid.mean() if len(valid) > 0 else 0.0
        gd["alt_freqs"].append(af)
        if af > gd["top_freq"]:
            gd["top_freq"] = af
            gd["top_site"] = site_labels[i]

    rows = []
    for gene, gd in sorted(gene_data.items()):
        row = {
            "gene": gene,
            "n_sites": gd["n_sites"],
            "top_site": gd["top_site"],
            "top_freq": round(gd["top_freq"], 4),
            "mean_af": round(np.mean(gd["alt_freqs"]), 4) if gd["alt_freqs"] else 0.0,
            "cds_start": gene_spans.get(gene, [0, 0])[0],
            "cds_end": gene_spans.get(gene, [0, 0])[1],
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    # 合并 popgen 数据 (如有)
    if df_popgen is not None and len(df_popgen) > 0:
        for _, prow in df.iterrows():
            if prow["cds_start"] == 0 and prow["cds_end"] == 0:
                continue
            gs, ge = prow["cds_start"], prow["cds_end"]
            mask = (df_popgen["BIN_MID"] >= gs) & (df_popgen["BIN_MID"] <= ge)
            pg = df_popgen[mask]
            if len(pg) > 0:
                df.loc[df["gene"] == prow["gene"], "mean_pi"] = round(pg["PI"].mean(), 6)
                df.loc[df["gene"] == prow["gene"], "mean_tajd"] = round(pg["TAJIMA_D"].mean(), 4)
            else:
                df.loc[df["gene"] == prow["gene"], "mean_pi"] = np.nan
                df.loc[df["gene"] == prow["gene"], "mean_tajd"] = np.nan

    return df


# ============================================================
# 5.9  带注释的可视化 (per-gene, 基因轨道)
# ============================================================

def plot_per_gene_variants(gene_df, out_path):
    """Per-gene 变异位点数柱状图"""
    assert HAS_VIZ, "需要 matplotlib"

    df = gene_df[gene_df["n_sites"] > 0].sort_values("n_sites", ascending=True)
    if df.empty:
        return

    plt.figure(figsize=(8, max(4, len(df) * 0.3)))
    colors = ["#e74c3c" if row["mean_af"] > 0.5 else "#3498db" for _, row in df.iterrows()]
    plt.barh(df["gene"], df["n_sites"], color=colors, edgecolor="black")

    for i, (_, row) in enumerate(df.iterrows()):
        label = f"  n={row['n_sites']}, top AF={row['top_freq']:.2f}"
        plt.text(row["n_sites"] + 0.2, i, label, va="center", fontsize=8)

    plt.xlabel("Number of Variant Sites", fontsize=12, fontweight="bold")
    plt.title("Variants per Gene", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[Per-gene图] 已保存: {out_path}")


def plot_popgen_with_genes(df_popgen, gene_spans, out_path):
    """π + Tajima's D 滑动窗 + 基因轨道 (三面板)"""
    assert HAS_VIZ, "需要 matplotlib"

    if df_popgen is None or df_popgen.empty:
        return
    if not gene_spans:
        plot_popgen_windows(df_popgen, out_path)
        return

    df = df_popgen.dropna(subset=["PI", "TAJIMA_D"])
    if df.empty:
        return

    import matplotlib.patches as mpatches

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 9),
                                         gridspec_kw={'height_ratios': [4, 4, 1]},
                                         sharex=True)

    # Panel A: π
    ax1.plot(df["BIN_MID"], df["PI"], color="#e74c3c", linewidth=1.2)
    ax1.fill_between(df["BIN_MID"], 0, df["PI"], color="#e74c3c", alpha=0.12)
    ax1.set_ylabel("Nucleotide Diversity (π)", fontsize=11, fontweight="bold")
    ax1.set_title("Sliding Window Population Genetics", fontsize=14, fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.3)

    # Panel B: Tajima's D
    ax2.plot(df["BIN_MID"], df["TAJIMA_D"], color="#2980b9", linewidth=1.2)
    ax2.fill_between(df["BIN_MID"], 0, df["TAJIMA_D"], color="#2980b9", alpha=0.12)
    ax2.axhline(y=0, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax2.set_ylabel("Tajima's D", fontsize=11, fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.3)

    # Panel C: gene tracks
    max_pos = int(df["BIN_END"].max())
    ax3.set_xlim(0, max_pos)
    ax3.set_ylim(0, 1)
    ax3.set_xlabel("Genomic Position (bp)", fontsize=12, fontweight="bold")
    ax3.get_yaxis().set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    gene_colors = ['#FF9F1C', '#2ECC71', '#3498DB', '#9B59B6', '#E74C3C', '#F1C40F',
                   '#1ABC9C', '#E67E22', '#2980B9', '#8E44AD']
    sorted_genes = sorted(gene_spans.items(), key=lambda x: x[1][0])

    for gi, (gname, (gs, ge)) in enumerate(sorted_genes):
        color = gene_colors[gi % len(gene_colors)]
        rect = mpatches.Rectangle((gs, 0.15), ge - gs, 0.7, linewidth=1,
                                  edgecolor='black', facecolor=color, alpha=0.85)
        ax3.add_patch(rect)
        mid = (gs + ge) / 2
        ax3.text(mid, 0.5, gname, ha='center', va='center', fontsize=9,
                 fontweight='bold', color='white')

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[PopGen+基因图] 已保存: {out_path}")


def plot_annotated_snp_heatmap(M, site_labels, site_annotations, samples, out_path):
    """SNP 矩阵热图，列按基因分组着色"""
    assert HAS_VIZ, "需要 matplotlib + seaborn"

    if M.shape[0] < 2 or M.shape[1] < 2:
        return

    # 为每个位点分配基因
    gene_list = []
    for label in site_labels:
        parts = label.split("_")
        try:
            pos, ref, alt = parts[1], parts[2], parts[3]
        except IndexError:
            gene_list.append("unknown")
            continue
        ann = site_annotations.get((pos, ref, alt), {})
        gene_list.append(ann.get("gene", "intergenic"))

    # 按基因排序位点
    gene_order = sorted(set(gene_list), key=lambda g: (
        0 if g == "intergenic" else 1 if g == "unknown" else 2, g))
    col_order = sorted(range(len(gene_list)),
                       key=lambda i: (gene_order.index(gene_list[i]), i))

    M_ordered = M[:, col_order]
    genes_ordered = [gene_list[i] for i in col_order]

    # 准备 DataFrame
    import pandas as _pd
    heatmap_labels = [site_labels[i] for i in col_order]
    df = _pd.DataFrame(M_ordered, index=samples, columns=heatmap_labels)
    df = df.fillna(0.5)  # missing → grey

    # 行聚类
    g = sns.clustermap(df, cmap="RdBu_r", center=0.5,
                       row_cluster=True, col_cluster=False,
                       figsize=(max(12, len(heatmap_labels) * 0.2),
                                max(8, len(samples) * 0.3)),
                       cbar_kws={"label": "Genotype (0=ref, 1=alt, 0.5=NA)"})
    g.ax_heatmap.set_title("SNP Genotype Matrix (sites ordered by gene)",
                           fontsize=14, fontweight="bold", pad=20)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[注释热图] 已保存: {out_path}")


# ============================================================
# 6. 主流程
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="病毒VCF合并与下游群体遗传学分析管线 v3.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 默认用法 (kinship PCA + Jaccard/Hamming + 树 + 6 图)
  python virus_vcf_pipeline.py -d ./vcfs/ -o ./out/ --prefix my_virus \\
      --dist-metrics both --tree --visualize

  # 完整下游 (QC + LD + PopGen + per-gene 汇总)
  python virus_vcf_pipeline.py -d ./vcfs/ -o ./out/ --prefix my_virus \\
      --qc --max-snps 500 --dist-metrics both --tree --visualize \\
      --ld --popgen-windows \\
      --snpeff-dir ../virus-SnpEff/ --snpgenie --snpgenie-ref ref.fa --snpgenie-gtf ref.gtf \\
      --gene-summary

  # 合并 SnpEff + SNPGenie 一键
  python virus_vcf_pipeline.py -d ./vcfs/ -o ./out/ \\
      --visualize --gene-summary \\
      --snpeff-dir ../virus-SnpEff/ \\
      --snpgenie --snpgenie-ref ref.fa --snpgenie-gtf ref.gtf
"""
    )
    # ---- v1 兼容参数 ----
    parser.add_argument("-d", "--dir", required=True, help="包含样本 VCF 的目录")
    parser.add_argument("-p", "--pattern", default="**/*.filtered.vcf",
                        help="VCF 文件匹配 glob 模式 (默认: **/*.filtered.vcf)")
    parser.add_argument("-o", "--out_dir", default=".",
                        help="输出目录 (默认: 当前目录)")
    parser.add_argument("--prefix", default="merged_virus",
                        help="输出前缀 (默认: merged_virus)")

    # ---- 新参数 ----
    parser.add_argument("--qc", action="store_true",
                        help="开启样本级质控, 输出 qc_summary.tsv")
    parser.add_argument("--max-snps", type=int, default=0,
                        help="过滤超过此突变数的样本 (0=不过滤)")
    parser.add_argument("--max-missing", type=float, default=0.0,
                        help="过滤缺失率超过此值的样本 (0=不过滤; 建议 0.2)")

    parser.add_argument("--snp-matrix", action="store_true",
                        help="导出 SNP 矩阵 (样本×位点, 0/1/NA)")
    parser.add_argument("--dist-metrics", choices=["jaccard", "hamming", "both"],
                        default=None, help="在 Python 中计算距离矩阵 (替代/补充 VCF2Dis)")
    parser.add_argument("--tree", action="store_true",
                        help="从距离矩阵构建 NJ/UPGMA 树")

    parser.add_argument("--visualize", action="store_true",
                        help="生成全套图表 (PCA/热图/AFS/树图/LD/PopGen/per-gene柱状图/注释热图/基因轨道)")
    parser.add_argument("--pca-method", choices=["kinship", "genotype"], default="kinship",
                        help="PCA 方法: kinship (VCF2PCACluster 等价) | genotype (中心化)")

    parser.add_argument("--ld", action="store_true",
                        help="寻找高 r² 共突变位点对 (上位效应), 输出 epistatic_co_mutations.tsv")
    parser.add_argument("--popgen-windows", action="store_true",
                        help="滑动窗计算 π 与 Tajima's D")
    parser.add_argument("--win-size", type=int, default=0,
                        help="滑动窗口大小 (0=自动)")
    parser.add_argument("--win-step", type=int, default=0,
                        help="滑动窗口步长 (0=自动)")

    parser.add_argument("--snpeff-dir", type=str, default=None,
                        help="Per-sample SnpEff 注释目录 (virus-SnpEff/) 用于汇总 per-site 注释")
    parser.add_argument("--snpgenie", action="store_true",
                        help="在合并 VCF 上运行 SNPGenie (需 --snpgenie-ref 和 --snpgenie-gtf)")
    parser.add_argument("--snpgenie-ref", type=str, default=None,
                        help="SNPGenie 参考 FASTA")
    parser.add_argument("--snpgenie-gtf", type=str, default=None,
                        help="SNPGenie GTF 注释文件")
    parser.add_argument("--gene-summary", action="store_true",
                        help="生成 per-gene 汇总表与可视化 (需 --snpeff-dir 或 --snpgenie-gtf)")

    parser.add_argument("--use-vcf2tools", action="store_true",
                        help="使用 VCF2PCACluster 和 VCF2Dis (默认已改用 Python sklearn/numpy)")

    parser.add_argument("--ivar", action="store_true",
                        help="输入为 iVar 产物, 启用 awk 修复缺失的 FORMAT/GT 列 (Freebayes/bcftools 不要加此参数)")

    args = parser.parse_args()

    input_dir = os.path.abspath(args.dir)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    prefix = args.prefix

    # ---- 依赖检查 ----
    extra_pkgs = []
    if args.snp_matrix or args.dist_metrics or args.ld or args.popgen_windows:
        extra_pkgs.append("numpy")
    if args.tree:
        extra_pkgs.append("scipy")
    if args.visualize or args.gene_summary or (args.ld and HAS_VIZ) or (args.popgen_windows and HAS_VIZ):
        extra_pkgs.extend(["viz", "numpy"])
    if args.visualize:
        extra_pkgs.append("sklearn")
    check_dependencies(extra_pkgs)

    if args.use_vcf2tools:
        missing = [t for t in ["VCF2PCACluster", "VCF2Dis"] if shutil.which(t) is None]
        if missing:
            print(f"[错误] --use-vcf2tools 需要: {', '.join(missing)}")
            sys.exit(1)

    # ---- 目录准备 ----
    stats_dir = os.path.join(out_dir, "stats")
    matrix_dir = os.path.join(out_dir, "matrices")
    tree_dir = os.path.join(out_dir, "trees")
    fig_dir = os.path.join(out_dir, "figs")
    for d in [stats_dir, matrix_dir, tree_dir, fig_dir]:
        os.makedirs(d, exist_ok=True)

    # ---- 查找 VCF 文件 ----
    search_pattern = os.path.join(input_dir, args.pattern)
    vcf_files = glob.glob(search_pattern, recursive=True)
    if not vcf_files:
        vcf_files = glob.glob(os.path.join(input_dir, "*.filtered.vcf"))
    if not vcf_files:
        print(f"[错误] 在 {input_dir} 中未找到 VCF 文件。")
        sys.exit(1)

    print(f"\n找到 {len(vcf_files)} 个 VCF 文件。")

    # ---- 工作目录 ----
    work_dir = os.path.join(out_dir, f".{prefix}_workdir")
    os.makedirs(work_dir, exist_ok=True)

    # ============================================================
    # Step 1: VCF 重命名 + bgzip/tabix (标准 VCF 或 iVar VCF)
    # ============================================================
    caller_tag = "iVar" if args.ivar else "Freebayes/bcftools"
    print(f"\n================ [步骤 1] 准备 VCF ({caller_tag}) ================")
    ready_vcfs = []

    for vcf in vcf_files:
        if args.ivar:
            # iVar INFO-only VCF: 从文件名提取样本名, awk 补 FORMAT/GT 列
            filename = os.path.basename(vcf)
            sample_name = filename.split('.')[0].split('_')[0]
            out_vcf_gz = os.path.join(work_dir, f"{sample_name}.vcf.gz")

            awk_cmd = (
                f"awk -v sample='{sample_name}' '"
                r'BEGIN { FS="\t"; OFS="\t" } '
                r'/^##/ { print $0; next } '
                r'/^#CHROM/ { '
                r'  if (NF < 9) { '
                r'      print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; '
                r'      print $0, "FORMAT", sample; '
                r'      add_cols=1 '
                r'  } else { '
                r'      $10=sample; '
                r'      print $0; '
                r'      add_cols=0 '
                r'  }; '
                r'  next '
                r'} '
                r'{ if (add_cols==1) { print $0, "GT", "1/1" } else { print $0 } }'
                f"' '{vcf}' | bgzip -c > '{out_vcf_gz}'"
            )
            run_command(awk_cmd)
        else:
            # 标准 VCF (Freebayes / bcftools): 先读 VCF 已有的样本名
            sample_name = run_command(
                f"bcftools query -l '{vcf}'", return_output=True
            ).strip().split('\n')[0].strip()
            # Freebayes 单样本 VCF 常写 "unknown"，回退到文件名
            if not sample_name or sample_name == "unknown":
                sample_name = os.path.basename(vcf).split('.')[0]
            out_vcf_gz = os.path.join(work_dir, f"{sample_name}.vcf.gz")

            # 重设样本名（如果 VCF 里是 unknown）
            reheader_map = os.path.join(work_dir, f"_rename_{sample_name}.txt")
            orig_name = run_command(
                f"bcftools query -l '{vcf}'", return_output=True
            ).strip().split('\n')[0].strip()
            if orig_name and orig_name != sample_name:
                with open(reheader_map, "w") as rf:
                    rf.write(f"{orig_name} {sample_name}\n")
                run_command(f"bcftools reheader -s '{reheader_map}' '{vcf}' | bgzip -c > '{out_vcf_gz}'")
            elif vcf.endswith(".gz"):
                run_command(f"zcat '{vcf}' | bgzip -c > '{out_vcf_gz}'")
            else:
                run_command(f"bgzip -c '{vcf}' > '{out_vcf_gz}'")

        run_command(f"tabix -p vcf '{out_vcf_gz}'")
        ready_vcfs.append(out_vcf_gz)

    # ============================================================
    # Step 2: Merge (保留缺失为 ./. 不强制填参考)
    # ============================================================
    print("\n================ [步骤 2] 合并 VCF ================")
    raw_merged_vcf = os.path.join(work_dir, "raw_merged.vcf.gz")

    # 使用文件列表而非命令行拼接, 防御超大队列 ARG_MAX 溢出
    if len(ready_vcfs) == 1:
        shutil.copy(ready_vcfs[0], raw_merged_vcf)
        print(f"[信息] 仅 1 个样本, 跳过合并。")
    else:
        vcf_list_file = os.path.join(work_dir, "vcf_merge_list.txt")
        with open(vcf_list_file, "w") as f:
            f.write("\n".join(ready_vcfs) + "\n")
        run_command(f"bcftools merge -l '{vcf_list_file}' -O z -o '{raw_merged_vcf}'")

    # ---- SNP filter ----
    merged_vcf = os.path.join(out_dir, f"{prefix}.vcf.gz")
    run_command(f"bcftools view -m2 -M2 -v snps '{raw_merged_vcf}' -O z -o '{merged_vcf}'")
    run_command(f"tabix -p vcf '{merged_vcf}'")

    snp_count = int(run_command(f"bcftools view -H '{merged_vcf}' | wc -l", return_output=True).strip() or 0)
    if snp_count == 0:
        print(f"\n[警告] 过滤后剩余 0 个 SNP, 终止。")
        shutil.rmtree(work_dir)
        return
    print(f"[信息] 合并后保留 {snp_count} 个 bi-allelic SNP 位点。")

    # ============================================================
    # Step 2.5: QC (如果开启)
    # ============================================================
    qc_merged = merged_vcf  # 默认直接用原始合并结果
    excluded_samples = set()

    if args.qc:
        print("\n================ [步骤 2.5] 样本 QC ================")
        stats = qc_sample_stats(merged_vcf)

        if args.max_snps > 0 or args.max_missing > 0:
            qc_merged = os.path.join(work_dir, "qc_filtered.vcf.gz")
            kept, excluded_samples = qc_filter_samples(
                merged_vcf, qc_merged, stats,
                max_snps=args.max_snps, max_missing=args.max_missing
            )
        else:
            # 只出报告不过滤
            kept = list(stats.keys())
            excluded_samples = set()

        qc_path = os.path.join(stats_dir, "qc_summary.tsv")
        qc_summary_table(stats, excluded_samples, qc_path)

    # ============================================================
    # Step 3: SNP 矩阵导出
    # ============================================================
    M = None
    samples = None
    site_labels = None

    if args.snp_matrix or args.dist_metrics or args.visualize or args.ld or args.popgen_windows or args.gene_summary or args.snpgenie:
        print("\n================ [步骤 3] 导出 SNP 矩阵 ================")
        M, samples, site_labels = export_snp_matrix(
            qc_merged, os.path.join(matrix_dir, "snp_matrix.tsv")
        )
        if M is None:
            print("[错误] SNP 矩阵导出失败。")
            sys.exit(1)


    # ============================================================
    # Step 4: 距离矩阵 (Python)
    # ============================================================
    D_dict = {}
    Z_linkage = None

    if args.dist_metrics and M is not None and len(samples) > 0:
        print("\n================ [步骤 4] 计算距离矩阵 ================")
        D_dict = pairwise_distance_matrix(M, metric=args.dist_metrics)

        for metric_name, D_mat in D_dict.items():
            save_path = os.path.join(matrix_dir, f"distance_{metric_name}.tsv")
            save_distance_matrix(D_mat, samples, save_path)

    # ============================================================
    # Step 5: NJ 树
    # ============================================================
    if args.tree and M is not None and len(samples) > 0:
        print("\n================ [步骤 5] 构建 NJ/UPGMA 树 ================")
        if "jaccard" in D_dict:
            D_tree = D_dict["jaccard"]
        elif "hamming" in D_dict:
            D_tree = D_dict["hamming"]
        else:
            # 用 Jaccard 兜底
            D_dict2 = pairwise_distance_matrix(M, metric="jaccard")
            D_tree = D_dict2["jaccard"]

        tree_path = os.path.join(tree_dir, "tree.newick")
        Z_linkage, _ = build_nj_tree(D_tree, samples, tree_path)

    # ============================================================
    # Step 6: 可视化
    # ============================================================
    if args.visualize and M is not None and len(samples) > 0:
        print("\n================ [步骤 6] 生成可视化图表 ================")

        # PCA
        if M.shape[1] >= 2:
            plot_pca(M, samples, os.path.join(fig_dir, "pca.png"),
                     method=args.pca_method)

        # 距离热图
        if "jaccard" in D_dict:
            plot_distance_heatmap(D_dict["jaccard"], samples,
                                  os.path.join(fig_dir, "distance_clustermap.png"),
                                  title="Jaccard Distance Heatmap")
        elif "hamming" in D_dict:
            plot_distance_heatmap(D_dict["hamming"], samples,
                                  os.path.join(fig_dir, "distance_clustermap.png"),
                                  title="Hamming Distance Heatmap")

        # AFS
        plot_afs(M, os.path.join(fig_dir, "afs.png"))

        # 树图
        if Z_linkage is not None:
            plot_dendrogram(Z_linkage, samples,
                            os.path.join(fig_dir, "dendrogram.png"),
                            title="UPGMA Hierarchical Clustering")

    # ============================================================
    # Step 6a: 共突变分析 (r², 上位效应, 不画 LD decay)
    # ============================================================
    ld_df = None
    if args.ld and M is not None and len(samples) > 0:
        print("\n================ [步骤 6a] 寻找共突变位点对 (上位效应) ================")
        ld_df = compute_ld_r2(M, site_labels, min_r2=0.5)
        if ld_df is not None:
            ld_out = os.path.join(stats_dir, "epistatic_co_mutations.tsv")
            ld_df.to_csv(ld_out, sep="\t", index=False)
            print(f"[共突变] 已导出 {len(ld_df)} 对强连锁变异: {ld_out}")
            # 网络散点 + r² 分布图
            if args.visualize and HAS_VIZ and site_labels is not None:
                pos_arr = [int(l.split("_")[1]) if "_" in l else 0 for l in site_labels]
                plot_epistatic_network(ld_df, site_labels, pos_arr,
                                       os.path.join(fig_dir, "epistatic_network.png"))

    # ============================================================
    # Step 6b: 滑动窗 π + Tajima's D
    # ============================================================
    df_popgen = None
    if args.popgen_windows and M is not None and len(samples) > 0:
        print("\n================ [步骤 6b] 滑动窗 π + Tajima's D ================")
        df_popgen = compute_sliding_window_popgen(
            M, site_labels, win_size=args.win_size, step=args.win_step
        )
        if df_popgen is not None:
            df_popgen.to_csv(os.path.join(stats_dir, "popgen_windows.tsv"),
                             sep="\t", index=False)
            print(f"[PopGen] 滑动窗统计已保存: {stats_dir}/popgen_windows.tsv")
            if args.visualize:
                plot_popgen_windows(df_popgen, os.path.join(fig_dir, "popgen_windows.png"))

    # ============================================================
    # Step 6c: SnpEff per-site 注释汇总
    # ============================================================
    site_annotations = {}
    gene_spans = {}

    if args.snpeff_dir or args.gene_summary:
        if args.snpeff_dir:
            print("\n================ [步骤 6c] 汇总 per-site SnpEff 注释 ================")
            site_annotations = collect_per_site_annotations(args.snpeff_dir)

        # 解析基因坐标
        gtf_for_spans = args.snpgenie_gtf or (
            os.path.join(args.snpeff_dir, "..", "virus-annotations") if args.snpeff_dir else None
        )
        if gtf_for_spans:
            # 尝试从 GTF 目录找对应文件
            if os.path.isdir(gtf_for_spans):
                import glob as _glob2
                gtf_files = _glob2.glob(os.path.join(gtf_for_spans, "*.gtf"))
                if gtf_files:
                    gtf_for_spans = gtf_files[0]
            if os.path.isfile(gtf_for_spans):
                gene_spans = parse_gtf_gene_spans(gtf_for_spans)
                print(f"[基因坐标] 从 GTF 解析了 {len(gene_spans)} 个基因区间。")

    # ============================================================
    # Step 6d: 合并 SNPGenie
    # ============================================================
    snpgenie_products = None
    if args.snpgenie:
        if not args.snpgenie_ref or not args.snpgenie_gtf:
            print("[错误] --snpgenie 需要 --snpgenie-ref 和 --snpgenie-gtf")
        elif not shutil.which("snpgenie.pl"):
            print("[警告] 找不到 snpgenie.pl, 跳过合并 SNPGenie。")
        else:
            print("\n================ [步骤 6d] 合并 SNPGenie 选择压力分析 ================")
            sg_dir = os.path.join(out_dir, "snpgenie_merged")
            os.makedirs(sg_dir, exist_ok=True)

            vcf_work = os.path.join(sg_dir, "merged.vcf")
            fa_work = os.path.join(sg_dir, "ref.fasta")
            gtf_work = os.path.join(sg_dir, "ref.gtf")

            # 解压到工作目录
            run_command(f"bcftools view '{qc_merged}' -Ov -o '{vcf_work}'")
            shutil.copy(args.snpgenie_ref, fa_work)
            shutil.copy(args.snpgenie_gtf, gtf_work)

            run_command(
                f"cd '{sg_dir}' && snpgenie.pl --vcfformat=2 "
                f"--snpreport='{os.path.basename(vcf_work)}' "
                f"--fastafile='{os.path.basename(fa_work)}' "
                f"--gtffile='{os.path.basename(gtf_work)}'",
                check=False
            )

            # 展平 SNPGenie_Results
            sg_results = os.path.join(sg_dir, "SNPGenie_Results")
            if os.path.isdir(sg_results):
                for item in os.listdir(sg_results):
                    src = os.path.join(sg_results, item)
                    dst = os.path.join(sg_dir, item)
                    if os.path.exists(dst):
                        if os.path.isdir(dst):
                            shutil.rmtree(dst)
                        else:
                            os.remove(dst)
                    shutil.move(src, dst)
                shutil.rmtree(sg_results)

            # 解析 product_results.txt
            prod_file = os.path.join(sg_dir, "product_results.txt")
            if os.path.exists(prod_file):
                try:
                    snpgenie_products = pd.read_csv(prod_file, sep="\t")
                    print(f"[SNPGenie] 解析了 {len(snpgenie_products)} 个基因产物的选择压力数据。")
                except Exception as e:
                    print(f"[警告] 解析 product_results.txt 失败: {e}")

            # 汇总 per-gene SNPGenie 数据
            if snpgenie_products is not None and len(snpgenie_products) > 0:
                sg_cols = [c for c in snpgenie_products.columns]
                sg_out = os.path.join(stats_dir, "snpgenie_per_gene.tsv")
                # 提取关键列
                key_cols = ["product"]
                for c in ["N_diffs_vs_ref", "S_diffs_vs_ref", "N_sites", "S_sites",
                           "mean_dN_vs_ref", "mean_dS_vs_ref", "piN", "piS"]:
                    if c in sg_cols:
                        key_cols.append(c)
                sg_sub = snpgenie_products[key_cols].copy()
                # 计算 dN/dS 和 πN/πS
                if "mean_dN_vs_ref" in sg_sub.columns and "mean_dS_vs_ref" in sg_sub.columns:
                    dS = pd.to_numeric(sg_sub["mean_dS_vs_ref"], errors="coerce")
                    dN = pd.to_numeric(sg_sub["mean_dN_vs_ref"], errors="coerce")
                    sg_sub["dNdS"] = np.where(dS > 0, dN / dS, np.nan)
                if "piN" in sg_sub.columns and "piS" in sg_sub.columns:
                    piS = pd.to_numeric(sg_sub["piS"], errors="coerce")
                    piN = pd.to_numeric(sg_sub["piN"], errors="coerce")
                    sg_sub["piN_piS"] = np.where(piS > 0, piN / piS, np.nan)
                sg_sub.to_csv(sg_out, sep="\t", index=False)
                print(f"[SNPGenie] Per-gene 汇总: {sg_out}")

    # ============================================================
    # Step 6e: Per-gene 汇总表
    # ============================================================
    gene_summary_df = None
    if args.gene_summary and M is not None and len(samples) > 0:
        print("\n================ [步骤 6e] 生成 per-gene 汇总表 ================")
        gene_summary_df = build_per_gene_summary(
            M, site_labels, site_annotations, df_popgen, gene_spans
        )

        # 合并 SNPGenie per-gene 数据
        if snpgenie_products is not None and len(snpgenie_products) > 0:
            for _, srow in snpgenie_products.iterrows():
                prod_name = str(srow.get("product", ""))
                # 尝试匹配基因名
                for gcol in [prod_name, prod_name.replace("_CDS", ""), prod_name.split("_")[0]]:
                    if gcol in gene_summary_df["gene"].values:
                        mask = gene_summary_df["gene"] == gcol
                        if "mean_dN_vs_ref" in snpgenie_products.columns:
                            gene_summary_df.loc[mask, "dNdS"] = float(srow.get("mean_dN_vs_ref", np.nan)) / float(srow.get("mean_dS_vs_ref", 1)) if float(srow.get("mean_dS_vs_ref", 0)) > 0 else np.nan
                        break

        gene_path = os.path.join(stats_dir, "per_gene_summary.tsv")
        gene_summary_df.to_csv(gene_path, sep="\t", index=False)
        print(f"[Per-gene] 汇总表: {gene_path}")

    # ============================================================
    # Step 6f: 带注释的可视化
    # ============================================================
    if args.visualize and gene_summary_df is not None and len(gene_summary_df) > 0:
        print("\n================ [步骤 6f] 生成基因级可视化 ================")

        # Per-gene 变异数柱状图
        plot_per_gene_variants(gene_summary_df,
                               os.path.join(fig_dir, "per_gene_variants.png"))

        # SNP 矩阵带注释热图
        if site_annotations:
            plot_annotated_snp_heatmap(
                M, site_labels, site_annotations, samples,
                os.path.join(fig_dir, "snp_heatmap_annotated.png"))

        # PopGen + 基因轨道
        if df_popgen is not None and gene_spans:
            plot_popgen_with_genes(
                df_popgen, gene_spans,
                os.path.join(fig_dir, "popgen_with_genes.png"))

    # ============================================================
    # Step 7: VCF2工具 (仅当 --use-vcf2tools)
    # ============================================================
    if args.use_vcf2tools:
        print("\n================ [步骤 7] VCF2PCACluster / VCF2Dis (旧版兼容) ================")

        # 短 ID 映射
        if len(excluded_samples) > 0:
            vcf_for_tools = qc_merged
        else:
            vcf_for_tools = merged_vcf

        samples_str = run_command(f"bcftools query -l '{vcf_for_tools}'", return_output=True)
        samples_list = samples_str.strip().split('\n')

        mapping_dict = {}
        rename_lines = []
        for i, s_name in enumerate(samples_list):
            short_id = f"S{i:05d}"
            mapping_dict[short_id] = s_name
            rename_lines.append(f"{s_name} {short_id}")

        mapping_file = os.path.join(work_dir, "rename_map.txt")
        with open(mapping_file, "w") as f:
            f.write("\n".join(rename_lines) + "\n")

        short_vcf = os.path.join(work_dir, "short_merged.vcf.gz")
        run_command(f"bcftools reheader -s '{mapping_file}' '{vcf_for_tools}' > '{short_vcf}'")
        run_command(f"tabix -p vcf '{short_vcf}'")

        # PCA
        pca_prefix = os.path.join(out_dir, f"{prefix}_PCA")
        run_command(f"VCF2PCACluster -InVCF '{short_vcf}' -OutPut '{pca_prefix}'")

        # Distance
        dis_mat = os.path.join(out_dir, f"{prefix}_Distance.mat")
        run_command(f"VCF2Dis -InPut '{short_vcf}' -OutPut '{dis_mat}'")

        # 还原样本名
        if os.path.exists(dis_mat):
            with open(dis_mat, 'r') as f:
                lines = f.readlines()
            with open(dis_mat, 'w') as f:
                for line in lines:
                    parts = line.rstrip('\n').split('\t')
                    new_parts = [mapping_dict.get(p, p) for p in parts]
                    f.write('\t'.join(new_parts) + '\n')

        pca_file = f"{pca_prefix}.pca"
        if os.path.exists(pca_file):
            with open(pca_file, 'r') as f:
                lines = f.readlines()
            with open(pca_file, 'w') as f:
                for line in lines:
                    parts = line.strip().split()
                    if not parts:
                        continue
                    parts[0] = mapping_dict.get(parts[0], parts[0])
                    f.write('\t'.join(parts) + '\n')

    # ---- 清理 ----
    print("\n================ 清理中间文件 ================")
    shutil.rmtree(work_dir)

    # ---- 最终摘要 ----
    print(f"\n{'='*55}")
    print(f"  流程完成")
    print(f"{'='*55}")
    print(f"  合并 VCF:     {merged_vcf}")
    if args.qc:
        print(f"  QC 报告:      {stats_dir}/qc_summary.tsv")
    if args.snp_matrix:
        print(f"  SNP 矩阵:     {matrix_dir}/snp_matrix.tsv")
    if args.dist_metrics:
        for m in D_dict:
            print(f"  距离矩阵({m}): {matrix_dir}/distance_{m}.tsv")
    if args.tree:
        print(f"  NJ 树:        {tree_dir}/tree.newick")
    if args.visualize:
        print(f"  PCA 图:       {fig_dir}/pca.png")
        print(f"  热图:         {fig_dir}/distance_clustermap.png")
        print(f"  AFS 图:       {fig_dir}/afs.png")
        if Z_linkage is not None:
            print(f"  树图:         {fig_dir}/dendrogram.png")
        if ld_df is not None:
            print(f"  共突变表:     {stats_dir}/epistatic_co_mutations.tsv")
        if df_popgen is not None:
            print(f"  PopGen 图:    {fig_dir}/popgen_windows.png")
    if args.ld and ld_df is not None:
        print(f"  共突变位点对:  {len(ld_df)} 对 (r² >= 0.5)")
    if args.popgen_windows and df_popgen is not None:
        print(f"  PopGen 表:    {stats_dir}/popgen_windows.tsv")
    if args.snpeff_dir:
        print(f"  Per-site注释:  {len(site_annotations)} 位点")
    if args.snpgenie and snpgenie_products is not None:
        print(f"  SNPGenie:     {stats_dir}/snpgenie_per_gene.tsv")
    if args.gene_summary and gene_summary_df is not None:
        print(f"  Per-gene汇总:  {stats_dir}/per_gene_summary.tsv")
        if args.visualize:
            print(f"  Per-gene图:    {fig_dir}/per_gene_variants.png")
            if site_annotations:
                print(f"  注释热图:      {fig_dir}/snp_heatmap_annotated.png")
            if gene_spans and df_popgen is not None:
                print(f"  PopGen+基因:   {fig_dir}/popgen_with_genes.png")
    if args.use_vcf2tools:
        print(f"  VCF2  PCA:    {out_dir}/{prefix}_PCA.pca")
        print(f"  VCF2 距离:    {out_dir}/{prefix}_Distance.mat")
    print(f"{'='*55}\n")


if __name__ == "__main__":
    main()
