#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
=======================================================================================
 SNPGenie Population Miner (V13.2 OMEGA 大满贯究极整合版)
 融合了统计检验、机器学习与泛病毒三维景观渲染引擎，防断点，全免疫异常。
=======================================================================================
"""

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
import glob
import argparse
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from scipy import stats

# ================================
# 外部依赖检查模块
# ================================
try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    HAS_ML = True
except ImportError:
    HAS_ML = False
    print("  -> [警告] 未安装 scikit-learn，机器学习聚类图将被安全跳过。")

try:
    from Bio import Entrez, SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("  -> [警告] 未安装 biopython，Fig10 基因注释放弃联网，将使用本地数据推测边界。")

warnings.filterwarnings('ignore')

# ================================
# 统一配色 (SCI 发表级 / Okabe-Ito 色盲友好)
#   与 virome_analysis_pipeline/virus_metadata_plot.py 的 OKABE_ITO 同序同成员。
#   2026-10-11: 原文件混用 ColorBrewer Reds/Blues、Set2/Set3、flare 与杂牌 hex。
# ================================
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000", "#999999"]
ORANGE, SKY, GREEN, YELLOW, BLUE, VERM, PURPLE = OKABE_ITO[:7]
BLACK, GREY = "#000000", "#999999"
CMAP_SEQ = "viridis"        # 深(低) -> 浅(高); 原 flare 走这支
CMAP_SEQ_R = "viridis_r"    # 浅(低) -> 深(高); 原 YlOrRd 走这支
CMAP_DIV = "RdBu_r"         # 发散型: 蓝(低) - 白(中) - 红(高)


def okabe_cycle(n):
    """返回 n 个 Okabe-Ito 循环色 (n > 7 时自动循环)。"""
    return [OKABE_ITO[i % 7] for i in range(max(0, int(n)))]


# 分支/树图专用色序: 只用中深色, 排除最浅的 YELLOW(#F0E442) 与 SKY(#56B4E9)。
# 原因: 152 样本 + color_threshold=0.7*max 会产生十几个簇, 颜色循环后大量落在
# 橙/黄/浅蓝上, 白底细线分支对比度严重不足(视觉"发虚")。
# 保留的 5 色 (深蓝/朱红/深绿/橙/玫红) 仍全部是 Okabe-Ito 成员, 色盲友好不变。
OKABE_BRANCH = [BLUE, VERM, GREEN, ORANGE, PURPLE]


def okabe_branch_cycle(n):
    """分支/树图 (scipy.dendrogram、sns.clustermap) 专用循环色。"""
    return [OKABE_BRANCH[i % len(OKABE_BRANCH)] for i in range(max(0, int(n)))]


# ================================
# 全局审美基建 (Publication Ready)
# ================================
def setup_plot_style():
    sns.set_theme(style="ticks", context="paper", font_scale=1.2)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    plt.rcParams['axes.linewidth'] = 1.3
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['axes.labelweight'] = 'bold'
    # 默认离散色环也统一到 Okabe-Ito: 未显式指定颜色的图 (含 scipy/seaborn
    # 树状图分支的 'C0'..'C9') 不再落到 matplotlib 默认 tab10
    # 2026-10-11: 改用 okabe_branch_cycle (排除浅黄 #F0E442 / 浅蓝 #56B4E9)。
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=okabe_branch_cycle(10))

def P_to_stars(p):
    if p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else: return "ns"

# ==========================================
# 1. 核心数据解析与低频质控引擎
# ==========================================
def parse_and_filter_data(input_dir):
    print(f"[*] 开启深空探测阵列，解析目录: {os.path.abspath(input_dir)}")
    pop_list, prod_list, site_list = [], [],[]
    summary_files = glob.glob(os.path.join(input_dir, "**", "population_summary.txt"), recursive=True)
    
    if not summary_files:
        raise FileNotFoundError("【致命错】未找到 population_summary.txt 数据集！")
        
    for s_file in summary_files:
        d_name = os.path.dirname(s_file)
        sample_id = os.path.basename(d_name).replace("SNPGenie_Results", "")
        sample_id = sample_id.strip("_") or os.path.basename(d_name)
        
        try:
            df_pop = pd.read_csv(s_file, sep='\t')
            if not df_pop.empty:
                r = df_pop.iloc[0].copy()
                r['Sample_ID'] = sample_id
                dn, ds = pd.to_numeric(r.get('mean_dN_vs_ref', np.nan), errors='coerce'), pd.to_numeric(r.get('mean_dS_vs_ref', np.nan), errors='coerce')
                pi_n, pi_s = pd.to_numeric(r.get('piN', np.nan), errors='coerce'), pd.to_numeric(r.get('piS', np.nan), errors='coerce')
                r['dN_dS_ratio'] = dn / ds if pd.notna(ds) and ds > 0 else np.nan
                r['piN_piS_ratio'] = pi_n / pi_s if pd.notna(pi_s) and pi_s > 0 else np.nan
                pop_list.append(r)
        except Exception: pass
            
        prod_f = os.path.join(d_name, "product_results.txt")
        if os.path.exists(prod_f):
            try:
                df = pd.read_csv(prod_f, sep='\t', on_bad_lines='skip')
                if not df.empty: 
                    df['Sample_ID'] = sample_id
                    prod_list.append(df)
            except Exception: pass

        site_f = os.path.join(d_name, "site_results.txt")
        if os.path.exists(site_f):
            try:
                df = pd.read_csv(site_f, sep='\t', on_bad_lines='skip')
                if not df.empty: 
                    df['Sample_ID'] = sample_id
                    site_list.append(df)
            except Exception: pass

    df_pop = pd.DataFrame(pop_list)
    df_prod = pd.concat(prod_list, ignore_index=True) if prod_list else pd.DataFrame()
    df_site = pd.concat(site_list, ignore_index=True) if site_list else pd.DataFrame()
    
    print(f"  -> 数据装载完工：抓取到 {len(df_pop)} 个样本, {len(df_prod)} 个产物记录与 {len(df_site)} 处变异印记。")
    return df_pop, df_prod, df_site

def apply_strict_isnv_filters(df_site):
    print("[*] 施加极严苛亚群多态性質控 (Cov >= 100, VAF 2.5% ~ 97.5%)...")
    if df_site.empty: return df_site
    df = df_site.copy()
    b_cols =[c for c in df.columns if c.strip().upper() in list('ACGT')]
    if len(b_cols) < 4: return df
    
    df[b_cols] = df[b_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    df['coverage'] = pd.to_numeric(df['coverage'], errors='coerce').fillna(0)
    df = df[df['coverage'] >= 100]
    
    max_c = df[b_cols].max(axis=1)
    df['Minor_Count'] = df['coverage'] - max_c
    df['VAF'] = df['Minor_Count'] / df['coverage']
    df_filtered = df[(df['VAF'] >= 0.025) & (df['VAF'] <= 0.975)]
    return df_filtered

# ==========================================
# [Fig 1] VAF 频次分布谱图
# ==========================================
def plot_vaf_spectrum(df_filtered, output_dir):
    print("  -> Fig 1: 绘制 VAF 等位基因频数谱...")
    if df_filtered.empty: return
    plt.figure(figsize=(10, 6))
    sns.histplot(data=df_filtered, x='VAF', bins=50, kde=True, color=SKY, edgecolor='black')
    plt.axvline(x=0.05, color=VERM, linestyle='--', label='Rare Variants (5%)')
    plt.title('Fig 1: Variant Allele Frequency (VAF) Spectrum of iSNVs', fontsize=14, fontweight='bold')
    plt.xlabel('Minor Allele Frequency (VAF)', fontweight='bold')
    plt.ylabel('Density', fontweight='bold')
    plt.legend(); sns.despine(); plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Fig01_VAF_Spectrum.png"))
    plt.savefig(os.path.join(output_dir, "Fig01_VAF_Spectrum.pdf"))
    plt.close()

# ==========================================
# [Fig 2] 深度残差回归校正 iSNV 密度箱线图
# ==========================================
def calculate_adjusted_n_per_kb(df_site_filtered, df_prod, output_dir):
    print("  -> Fig 2: 执行覆盖度对数残差混杂矫正...")
    if df_site_filtered.empty or df_prod.empty: return
    df_p = df_prod.copy()
    df_p['gene_length'] = pd.to_numeric(df_p['N_sites'], errors='coerce') + pd.to_numeric(df_p['S_sites'], errors='coerce')
    gene_info = df_p.groupby(['Sample_ID', 'product'])['gene_length'].mean().reset_index()
    
    isnv_cnts = df_site_filtered.groupby(['Sample_ID', 'product']).agg(n_isnvs=('site', 'count'), mean_cov=('coverage', 'mean')).reset_index()
    df_m = pd.merge(gene_info, isnv_cnts, on=['Sample_ID', 'product'], how='left').fillna({'n_isnvs': 0})
    df_m['mean_cov'] = df_m['mean_cov'].fillna(df_m['mean_cov'].median() if not pd.isna(df_m['mean_cov'].median()) else 100)
    df_m['n_per_kb'] = df_m['n_isnvs'] / (df_m['gene_length'] / 1000.0)
    df_m['cov_log'] = np.log2(df_m['mean_cov'] + 1)
    
    v_mask = (df_m['mean_cov'] > 0) & (df_m['n_per_kb'] > 0)
    if v_mask.sum() > 5:
        sl, itc, _, _, _ = stats.linregress(df_m.loc[v_mask, 'cov_log'], df_m.loc[v_mask, 'n_per_kb'])
        df_m['expected'] = sl * df_m['cov_log'] + itc
        df_m['n_per_kb_adjusted'] = np.clip(df_m['n_per_kb'] - df_m['expected'] + df_m['n_per_kb'].mean(), 0, None)
    else:
        df_m['n_per_kb_adjusted'] = df_m['n_per_kb']
        
    plt.figure(figsize=(10, 6))
    # saturation=1: seaborn 默认 0.75 会把 Okabe 色去饱和 (实测 #E69F00 -> #C9941D),
    # 渲染出的就不是精确 Okabe-Ito 了。
    sns.boxplot(data=df_m, x='product', y='n_per_kb_adjusted', saturation=1,
                palette=okabe_cycle(df_m['product'].nunique()), showfliers=False)
    sns.stripplot(data=df_m, x='product', y='n_per_kb_adjusted', color='black', alpha=0.5, jitter=True)
    plt.title('Fig 2: Depth-Adjusted iSNV Density Across Viral Genes', fontsize=14, fontweight='bold')
    plt.xlabel('Viral Product', fontweight='bold'); plt.ylabel('Adjusted # iSNVs per Kb', fontweight='bold')
    sns.despine(); plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Fig02_Adjusted_iSNV_Density.png"))
    plt.savefig(os.path.join(output_dir, "Fig02_Adjusted_iSNV_Density.pdf"))
    plt.close()

# ==========================================
#[Fig 3 & 4] Inter & Intra-host JointPlots
# ==========================================
def plot_joint_dynamics(df_pop, output_dir):
    print("  -> Fig 3 & 4: 渲染宏观与微观演化双动量...")
    # Fig 3 (Inter-host)
    _cols_dnds = ['mean_dN_vs_ref', 'mean_dS_vs_ref']
    if not all(c in df_pop.columns for c in _cols_dnds):
        print("  -> [跳过 Fig 3&4] 未检测到 dN/dS 列 (可能为非蛋白编码序列)")
        return
    df_inter = df_pop.dropna(subset=_cols_dnds).copy()
    if not df_inter.empty:
        g = sns.JointGrid(data=df_inter, x='mean_dS_vs_ref', y='mean_dN_vs_ref', height=6)
        g.plot_joint(sns.scatterplot, alpha=0.7, color=BLUE, s=60, edgecolor='white')
        g.plot_marginals(sns.histplot, kde=True, color=BLUE, alpha=0.6)
        m_x = df_inter['mean_dS_vs_ref'].max() * 1.1
        m_y = df_inter['mean_dN_vs_ref'].max() * 1.2
        m_max = min(m_x, m_y)
        g.ax_joint.plot([0, m_max], [0, m_max], 'r--', label='dN=dS')
        g.ax_joint.legend(loc='upper left')
        g.ax_joint.set_xlabel('dS (Synonymous Divergence)', fontweight='bold')
        g.ax_joint.set_ylabel('dN (Nonsynonymous Divergence)', fontweight='bold')
        plt.suptitle("Fig 3: Inter-host Evolutionary Dynamics", y=1.03, fontsize=14, fontweight='bold')
        plt.savefig(os.path.join(output_dir, "Fig03_InterHost_dNdS.png"), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, "Fig03_InterHost_dNdS.pdf"), bbox_inches='tight'); plt.close()

    # Fig 4 (Intra-host)
    df_intra = df_pop[(df_pop.get('piN', 0) > 0) | (df_pop.get('piS', 0) > 0)].dropna(subset=['piN', 'piS']).copy()
    if not df_intra.empty:
        g2 = sns.JointGrid(data=df_intra, x='piS', y='piN', height=6)
        g2.plot_joint(sns.scatterplot, alpha=0.7, color=VERM, s=60, edgecolor='white')
        g2.plot_marginals(sns.histplot, kde=True, color=VERM, alpha=0.6)
        m_x2 = df_intra['piS'].max() * 1.1
        m_y2 = df_intra['piN'].max() * 1.2
        m_max2 = min(m_x2, m_y2)
        g2.ax_joint.plot([0, m_max2], [0, m_max2], 'k--', label='\u03C0N=\u03C0S')
        g2.ax_joint.legend(loc='upper left')
        g2.ax_joint.set_xlabel('\u03C0S (Synonymous Diversity)', fontweight='bold')
        g2.ax_joint.set_ylabel('\u03C0N (Nonsynonymous Diversity)', fontweight='bold')
        plt.suptitle("Fig 4: Intra-host Quasispecies Dynamics", y=1.03, fontsize=14, fontweight='bold')
        plt.savefig(os.path.join(output_dir, "Fig04_IntraHost_pi.png"), bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, "Fig04_IntraHost_pi.pdf"), bbox_inches='tight'); plt.close()

        # Kruskal-Wallis on \u03C0N/\u03C0S ratio across genes
        if 'product' in df_intra.columns and df_intra['product'].nunique() >= 2:
            df_intra_kw = df_intra.copy()
            df_intra_kw['pi_ratio'] = df_intra_kw['piN'] / (df_intra_kw['piS'] + 1e-10)
            groups = [g.dropna().values for _, g in df_intra_kw.groupby('product')['pi_ratio'] if len(g) >= 3]
            if len(groups) >= 2:
                h, p = stats.kruskal(*groups)
                with open(os.path.join(output_dir, "Stats_Kruskal_Wallis.csv"), "a") as f:
                    genes_str = " vs ".join(df_intra['product'].unique())
                    f.write(f"Kruskal-Wallis_\u03C0N/\u03C0S,{h:.4f},{p:.6g},{genes_str}\n")

# ==========================================
# [Fig 5] dN/dS Wilcoxon 统计小提琴图
# ==========================================
def plot_gene_dnds_with_stats(df_prod, output_dir):
    print("  -> Fig 5: 提取小提琴统计推衍图层...")
    if df_prod.empty or 'mean_dS_vs_ref' not in df_prod.columns:
        print("    -> [跳过 Fig 5] 缺少 dN/dS 列 (非蛋白编码序列)")
        return
    df_valid = df_prod[pd.to_numeric(df_prod['mean_dS_vs_ref'], errors='coerce') > 0].copy()
    if df_valid.empty: return
    df_valid['dNdS'] = pd.to_numeric(df_valid['mean_dN_vs_ref']) / pd.to_numeric(df_valid['mean_dS_vs_ref'])
    
    plt.figure(figsize=(11, 6))
    ax = plt.gca()
    sns.violinplot(data=df_valid, x='product', y='dNdS', hue='product', saturation=1,
                   inner="box", palette=okabe_cycle(df_valid['product'].nunique()),
                   legend=False, ax=ax)
    sns.stripplot(data=df_valid, x='product', y='dNdS', color='black', alpha=0.3, jitter=True, ax=ax)
    plt.axhline(y=1, color=VERM, linestyle='--', label='Neutral Selection')

    # 稳健量纲 + 离群隔离: 实测单个基因 (G) 的 dN/dS 可达 ~30, 既把其余 5 个基因
    # 压成平线, 又会把小提琴 KDE 的带宽撑爆 -> G 整根变成矩形。故:
    #   ① 以 99 分位定显示上限 (只收缩量纲, 不删数据)
    #   ② 小提琴/散点只画未越界数据, 越界点以文字说明个数与最大值
    y_max = float(df_valid['dNdS'].max())
    q99 = float(np.nanpercentile(df_valid['dNdS'].dropna(), 99)) if len(df_valid) else y_max
    y_top = min(y_max * 1.20, max(q99 * 1.6, 1.5))
    df_plot = df_valid[df_valid['dNdS'] <= y_top]
    n_clip = int(len(df_valid) - len(df_plot))

    plt.figure(figsize=(11, 6))
    ax = plt.gca()
    sns.violinplot(data=df_plot, x='product', y='dNdS', hue='product', saturation=1,
                   inner="box", palette=okabe_cycle(df_plot['product'].nunique()),
                   legend=False, ax=ax)
    sns.stripplot(data=df_plot, x='product', y='dNdS', color='black', alpha=0.3, jitter=True, ax=ax)
    plt.axhline(y=1, color=VERM, linestyle='--', label='Neutral Selection')

    # 顶部留出注记带; 星号/KW 文本都放进坐标区内, 避免像旧版那样被甩出画布并压住标题
    ax.set_ylim(0, y_top * 1.18)
    for i, gene in enumerate(df_valid['product'].unique()):
        g_d = df_valid[df_valid['product'] == gene]['dNdS'].dropna()
        if len(g_d) >= 3:
            try:
                _, p_val = stats.wilcoxon([x - 1.0 for x in g_d], alternative='two-sided')
                ax.text(i, 0.880, P_to_stars(p_val), ha='center', va='bottom', fontsize=14,
                        color=VERM, fontweight='bold', transform=ax.get_xaxis_transform())
            except Exception: pass

    # Kruskal-Wallis: compare dN/dS distributions across genes
    gene_groups = [g.dropna().values for _, g in df_valid.groupby('product')['dNdS']]
    gene_groups = [g for g in gene_groups if len(g) >= 3]
    kw_h, kw_p = np.nan, np.nan
    if len(gene_groups) >= 2:
        try:
            kw_h, kw_p = stats.kruskal(*gene_groups)
            kw_text = f"Kruskal-Wallis H={kw_h:.2f}, p={kw_p:.2e}" if kw_p < 0.001 else f"Kruskal-Wallis H={kw_h:.2f}, p={kw_p:.4f}"
            ax.text(0.5, 0.945, kw_text, ha='center', va='bottom', fontsize=10,
                    color=BLUE, fontweight='bold', transform=ax.get_xaxis_transform())
        except Exception: pass

    if n_clip:
        ax.text(0.005, 0.995, f"▲ {n_clip} point(s) above axis (max dN/dS = {y_max:.1f})",
                transform=ax.transAxes, ha='left', va='top', fontsize=8, color=GREY)
    plt.title('Fig 5: Selection Pressures (Wilcoxon per-gene + Kruskal-Wallis across genes)', fontsize=13, fontweight='bold')
    plt.xlabel('Viral Gene', fontweight='bold')
    plt.ylabel('dN/dS Ratio', fontweight='bold')
    plt.legend(); sns.despine(); plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Fig05_Gene_dNdS_Stats.png"))
    plt.savefig(os.path.join(output_dir, "Fig05_Gene_dNdS_Stats.pdf")); plt.close()

    # Save Kruskal-Wallis result
    if not np.isnan(kw_h):
        with open(os.path.join(output_dir, "Stats_Kruskal_Wallis.csv"), "w") as f:
            f.write("Test,Statistic,P_value,Genes_Compared\n")
            genes_str = " vs ".join(df_valid['product'].unique())
            f.write(f"Kruskal-Wallis,{kw_h:.4f},{kw_p:.6g},{genes_str}\n")

# ==========================================
# [Fig 6] 10,000x Bootstrapping
# ==========================================
def perform_10000x_bootstrap_dnds(df_prod, output_dir, bs=10000):
    print(f"  -> Fig 6: 引爆 {bs} 轮蒙特卡洛纯化选择边界计算...")
    needed_cols = ['product', 'N_diffs_vs_ref', 'S_diffs_vs_ref', 'N_sites', 'S_sites']
    if df_prod.empty or not all(c in df_prod.columns for c in needed_cols):
        print("    -> [跳过 Fig 6] 缺少 bootstrap 所需列 (非蛋白编码序列)")
        return
    df = df_prod[needed_cols].copy().apply(pd.to_numeric, errors='ignore')
    results = []
    
    for gene in df['product'].unique():
        df_g = df[df['product'] == gene]
        if len(df_g) < 3: continue
        
        real_n_sites = df_g['N_sites'].sum()
        real_s_sites = df_g['S_sites'].sum()
        r_N = df_g['N_diffs_vs_ref'].sum() / real_n_sites if real_n_sites > 0 else 0
        r_S = df_g['S_diffs_vs_ref'].sum() / real_s_sites if real_s_sites > 0 else 0
        
        matrix = df_g[['N_diffs_vs_ref', 'S_diffs_vs_ref', 'N_sites', 'S_sites']].values
        n_samp = len(matrix)
        boot_dnds, boot_dn_m_ds = np.zeros(bs), np.zeros(bs)
        
        for b in range(bs):
            rx = matrix[np.random.choice(n_samp, n_samp, replace=True)]
            bs_nsites = np.nansum(rx[:, 2])
            bs_ssites = np.nansum(rx[:, 3])
            nd = np.nansum(rx[:, 0])/bs_nsites if bs_nsites > 0 else 0
            sd = np.nansum(rx[:, 1])/bs_ssites if bs_ssites > 0 else 0
            boot_dn_m_ds[b] = nd - sd
            boot_dnds[b]    = nd / sd if sd > 0 else np.nan
            
        valid_ratio = np.sum(~np.isnan(boot_dnds)) / bs
        if np.isnan(boot_dnds).all() or valid_ratio < 0.5:
            continue  # skip genes with too few valid bootstrap replicates

        prop_above = np.mean(boot_dn_m_ds >= 0)
        p_boot = 2 * min(prop_above, 1 - prop_above)
        results.append({
            'Gene': gene,
            'Real': r_N/r_S if r_S > 0 else np.nan,
            'CI_L': np.nanpercentile(boot_dnds, 2.5),
            'CI_U': np.nanpercentile(boot_dnds, 97.5),
            'Sig': P_to_stars(p_boot)
        })
    
    if not results: return
    df_res = pd.DataFrame(results).dropna(subset=['Real'])
    if df_res.empty: return
    
    plt.figure(figsize=(9, 5))
    x = np.arange(len(df_res))
    yerr_lower = np.maximum(0, df_res['Real'] - df_res['CI_L'])
    yerr_upper = np.maximum(0, df_res['CI_U'] - df_res['Real'])
    
    plt.errorbar(x, df_res['Real'], yerr=[yerr_lower, yerr_upper],
                 fmt='o', color=VERM, ecolor='black', elinewidth=2, capsize=5, markersize=8)
    plt.axhline(1.0, color=GREY, linestyle='--')
    plt.xticks(x, df_res['Gene'])
    for i, r in df_res.iterrows(): plt.text(i, r['CI_U']+0.05, r['Sig'], ha='center', fontweight='bold', color=BLUE)
    plt.title('Fig 6: 10,000x Bootstrapping Purifying Selection Validation', fontsize=14, fontweight='bold')
    plt.ylabel('dNdS with 95% Bootstrap CI', fontweight='bold')
    sns.despine(); plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Fig06_Bootstrapped_dNdS.png"))
    plt.savefig(os.path.join(output_dir, "Fig06_Bootstrapped_dNdS.pdf"))
    plt.close()

# ==========================================
# [Fig 7 & 8] 突变频谱与标数高频靶点
# ==========================================
def _alt_nt_from_counts(row, ref_nt):
    """从 SNPGenie site_results 的 A/C/G/T 计数列取「非参考碱基中计数最高者」。

    为什么必须用它：SNPGenie 的 maj_nt 是位点上的**多数碱基**，不是变异碱基。
    位点不变时 maj_nt == ref_nt，若直接用 ref_nt->maj_nt 当替换标签，
    就会产出 A->A / T->T 这类不存在的替换（历史 bug）。
    返回大写碱基；若无任何非参考碱基（即不变位点）则返回 None。
    """
    ref = str(ref_nt).strip().upper()
    best, best_n = None, 0.0
    for b in ('A', 'C', 'G', 'T'):
        if b == ref:
            continue
        try:
            n = float(row.get(b, 0) or 0)
        except (TypeError, ValueError):
            n = 0.0
        if n > best_n:
            best, best_n = b, n
    return best


def plot_signatures_and_hotspots(df_site_filtered, df_site_all, output_dir):
    print("  -> Fig 7 & 8: 输出真实替换频谱与变异热点...")
    # -- Fig 7 (突变频谱) --
    # 只统计「真实存在的替换」：优先 QC 过滤后的多态位点，为空时退回全量但仍排除不变位点；
    # 替换标签一律由 A/C/G/T 计数列推出真实变异等位基因（不再用 maj_nt）。
    needed = ['ref_nt', 'A', 'C', 'G', 'T']
    base = df_site_filtered if not df_site_filtered.empty else (
        df_site_all if df_site_all is not None and not df_site_all.empty else df_site_filtered)
    if base.empty or not all(c in base.columns for c in needed):
        print("    -> [跳过 Fig 7] site_results 缺少 ref_nt 或 A/C/G/T 计数列")
    else:
        df_mut = base[base['ref_nt'].astype(str).str.upper().isin(list('ACGT'))].copy()
        df_mut['alt_nt'] = df_mut.apply(lambda r: _alt_nt_from_counts(r, r['ref_nt']), axis=1)
        df_mut = df_mut.dropna(subset=['alt_nt'])
        df_mut['ref_nt'] = df_mut['ref_nt'].astype(str).str.upper()
        df_mut = df_mut[df_mut['ref_nt'] != df_mut['alt_nt']]   # 双保险：排除不变位点
        if not df_mut.empty:
            # [防崩溃补丁] 安全降级转换类型
            df_mut['pos'] = pd.to_numeric(df_mut.get('position_in_codon', np.nan), errors='coerce')
            df_mut = df_mut.dropna(subset=['pos'])
            df_mut['pos'] = df_mut['pos'].astype(int)

            df_mut['Type'] = df_mut['ref_nt'] + "->" + df_mut['alt_nt']
            pivot = df_mut.groupby(['Type', 'pos']).size().unstack(fill_value=0)

            plot_cols = [c for c in [1, 2, 3] if c in pivot.columns]
            if plot_cols:
                pivot = pivot[plot_cols]
                pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
                plt.figure(figsize=(10, 6))
                pivot.plot(kind='bar', stacked=True, color=okabe_cycle(len(pivot.columns)), ax=plt.gca(), edgecolor='black')
                plt.title('Fig 7: Mutational Spectrum & Codon Buffer Bias', fontsize=14, fontweight='bold')
                plt.legend(title='Codon Pos'); plt.xlabel('Substitution Motif', fontweight='bold'); plt.ylabel('Occurrences', fontweight='bold')
                sns.despine(); plt.tight_layout()
                plt.savefig(os.path.join(output_dir, "Fig07_Mutational_Spectrum.png"))
                plt.savefig(os.path.join(output_dir, "Fig07_Mutational_Spectrum.pdf")); plt.close()

    # -- Fig 8 (趋同热点靶点) --
    if 'class_vs_ref' in df_site_filtered.columns and all(c in df_site_filtered.columns for c in needed):
        df_hot = df_site_filtered[df_site_filtered['class_vs_ref'] == 'Nonsynonymous'].copy()
        df_hot['alt_nt'] = df_hot.apply(lambda r: _alt_nt_from_counts(r, r['ref_nt']), axis=1)
        df_hot = df_hot.dropna(subset=['alt_nt'])
        df_hot['ref_nt'] = df_hot['ref_nt'].astype(str).str.upper()
        df_hot = df_hot[df_hot['ref_nt'] != df_hot['alt_nt']]   # 排除不变位点
        if not df_hot.empty:
            tc = df_hot.groupby(['product', 'site', 'ref_nt', 'alt_nt']).size().reset_index(name='n').sort_values('n', ascending=False).head(12)
            tc['L'] = tc.apply(lambda x: f"{x['product']}:{x['site']} ({x['ref_nt']}->{x['alt_nt']})", axis=1)
            plt.figure(figsize=(10, 6))
            ax = sns.barplot(data=tc, x='n', y='L', palette=CMAP_SEQ, edgecolor='black')
            for c in ax.containers: ax.bar_label(c, padding=4, fontsize=10, fontweight='bold')
            plt.title('Fig 8: Top Nonsynonymous Hotspots (QC Passed)', fontsize=14, fontweight='bold')
            plt.xlim(0, tc['n'].max() * 1.15)
            plt.xlabel('Sample Count', fontweight='bold'); plt.ylabel('Viral Target', fontweight='bold')
            sns.despine(); plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "Fig08_Top_Hotspots.png"))
            plt.savefig(os.path.join(output_dir, "Fig08_Top_Hotspots.pdf")); plt.close()

# ==========================================
#[Fig 9] 双轨滑动多态窗口 (pi_N vs pi_S)
# ==========================================
def plot_dual_track_window(df_site_all, output_dir, window=50):
    print("  -> Fig 9: 分轨滑窗探测(pi_N 激荡红线 vs pi_S 免疫深蓝)...")
    if df_site_all.empty or 'pi' not in df_site_all.columns: return
    df = df_site_all[['site', 'pi', 'product', 'class_vs_ref']].copy().dropna(subset=['site'])
    df['pi'], df['site'] = pd.to_numeric(df['pi'], errors='coerce').fillna(0), pd.to_numeric(df['site']).astype(int)
    
    m_ns = df[df['class_vs_ref'] == 'Nonsynonymous'].groupby('site')['pi'].mean().reset_index()
    m_sy = df[df['class_vs_ref'] == 'Synonymous'].groupby('site')['pi'].mean().reset_index()
    
    min_st, max_st = int(df['site'].min()), int(df['site'].max())
    all_s = pd.DataFrame({'site': range(min_st, max_st + 1)})
    
    mr_n = pd.merge(all_s, m_ns, how='left').fillna(0)['pi'].rolling(window=window, center=True).mean().fillna(0)
    mr_s = pd.merge(all_s, m_sy, how='left').fillna(0)['pi'].rolling(window=window, center=True).mean().fillna(0)
    
    plt.figure(figsize=(14, 5))
    plt.plot(all_s['site'], mr_s, color=BLUE, lw=1.5, label='Synonymous (\u03C0S)')
    plt.plot(all_s['site'], mr_n, color=VERM, lw=1.5, label='Nonsynonymous (\u03C0N)')
    plt.fill_between(all_s['site'], mr_n, color=VERM, alpha=0.2)
    plt.title(f'Fig 9: Dual-Track Diversity Landscape (w={window}nt)', fontsize=15, fontweight='bold')
    plt.xlabel('Genomic Position', fontweight='bold'); plt.ylabel('Smoothed \u03C0', fontweight='bold')
    plt.legend()
    sns.despine(); plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Fig09_DualTrack_Window.png"))
    plt.savefig(os.path.join(output_dir, "Fig09_DualTrack_Window.pdf"))
    plt.close()

# ==========================================
# [Fig 10] 自动连通原生 GenBank 解析物理坐标 (Trinity基建)
# ==========================================
def fetch_gb_coordinates(accession, output_dir):
    """自动联网下载并规避 noncoding 掩盖 Bug"""
    if not HAS_BIOPYTHON: return pd.DataFrame(), 0
    Entrez.email = "viralevolution@ncbi.nlm.nih.gov"
    gb_file = os.path.join(output_dir, f"{accession}.gb")
    
    if not os.path.exists(gb_file):
        try:
            with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as h:
                with open(gb_file, "w") as f_out: f_out.write(h.read())
        except Exception: return pd.DataFrame(), 0
            
    spans =[]
    try:
        record = SeqIO.read(gb_file, "genbank")
        for feat in record.features:
            if feat.type == 'CDS':
                st = int(feat.location.start) + 1  
                ed = int(feat.location.end)
                gn = feat.qualifiers.get('gene', feat.qualifiers.get('product', ['ORF']))[0]
                spans.append({'product': gn, 'st': st, 'ed': ed})
        return pd.DataFrame(spans).sort_values('st').drop_duplicates('product').reset_index(drop=True), len(record.seq)
    except Exception: return pd.DataFrame(), 0

# ==========================================
#[Fig 10] Trinity 终极进化全境图
# ==========================================
def plot_trinity_landscape(df_site, df_prod, output_dir, window=50, ref_acc=""):
    print("  -> Fig 10: 架构封顶！渲染 Trinity 级泛病毒宏伟生态全景大作...")
    if df_site.empty or df_prod.empty: return
    df_s = df_site.copy(); df_s['site'] = pd.to_numeric(df_s['site'], errors='coerce')
    df_s = df_s.dropna(subset=['site'])
    
    # 启用神级自动坐标解析
    spans, max_genome_len = fetch_gb_coordinates(ref_acc, output_dir)
    
    if spans.empty: 
        df_v = df_s[~df_s['product'].str.contains('noncoding', case=False, na=False)]
        spans = df_v.dropna(subset=['product']).groupby('product')['site'].agg(st='min', ed='max').reset_index()
        max_l = int(df_s['site'].max()) + 50
    else: max_l = max_genome_len
    
    if spans.empty: return

    st_pi = df_s.groupby('site')['pi'].mean()
    st_c  = df_s['site'].unique()
    
    pi_a, snp_a = np.zeros(max_l+1), np.zeros(max_l+1)
    for p in st_pi.index: 
        if pd.notna(p) and int(p) <= max_l: pi_a[int(p)] = st_pi[p]
    for p in st_c: 
        if pd.notna(p) and int(p) <= max_l: snp_a[int(p)] = 1
        
    r_snp = pd.Series(snp_a).rolling(window, center=True).sum().fillna(0)
    r_pi  = pd.Series(pi_a).rolling(window, center=True).mean().fillna(0)
    
    if 'mean_dN_vs_ref' not in df_prod.columns: return
    
    # [BUGFIX]: 安全数据类型转换机制
    df_p = df_prod.copy()
    df_p['mean_dN_vs_ref'] = pd.to_numeric(df_p['mean_dN_vs_ref'], errors='coerce')
    df_p['mean_dS_vs_ref'] = pd.to_numeric(df_p['mean_dS_vs_ref'], errors='coerce')
    d_sum = df_p.groupby('product')[['mean_dN_vs_ref', 'mean_dS_vs_ref']].mean().reset_index()
    d_sum['dNdS'] = d_sum['mean_dN_vs_ref'] / (d_sum['mean_dS_vs_ref'] + 1e-6)
    ytop = min(d_sum['dNdS'].max() * 1.3, 10.0) if d_sum['dNdS'].max() > 1 else 1.3

    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(15, 9), sharex=True, gridspec_kw={'height_ratios':[1, 3.5, 2]})
    fig.subplots_adjust(hspace=0.1)

    a1.plot([1, max_l],[0, 0], color='black', lw=1.5, zorder=1)
    a1.set_ylim(-0.8, 0.8); a1.axis('off')
    cs = okabe_cycle(6)
    for i, r in spans.iterrows():
        a1.add_patch(patches.Rectangle((r['st'], -0.4), r['ed']-r['st'], 0.8, facecolor=cs[i%len(cs)], ec='black', lw=1.2))
        a1.text(r['st']+(r['ed']-r['st'])/2, 0, r['product'], ha='center', va='center', fontweight='bold', fontsize=12)
    a1.set_title(f"Fig 10: The Trinity Evolutionary Pan-Landscape ({ref_acc})", fontsize=16, fontweight='bold')

    a2.plot(np.arange(max_l+1), r_snp, color=SKY, lw=1.5, label='SNPs')
    a2.fill_between(np.arange(max_l+1), r_snp, color=SKY, alpha=0.25)
    a2.set_ylabel(f'SNP Counts (w={window})', color=SKY, fontweight='bold')
    a2_t = a2.twinx()
    a2_t.plot(np.arange(max_l+1), r_pi, color=GREY, ls='--', lw=1.5, label='Pi')
    a2_t.set_ylabel(r'Diversity ($\pi$)', color=GREY, fontweight='bold')

    a3.axhline(1.0, color=GREY, ls='--')
    for i, r in spans.iterrows():
        dr = d_sum[d_sum['product'] == r['product']]
        if not dr.empty:
            dv = min(dr['dNdS'].values[0], 10.0)
            c = VERM if dv > 1.0 else GREEN
            a3.bar(r['st']+(r['ed']-r['st'])/2, dv, width=r['ed']-r['st'], facecolor=c, ec='black', alpha=0.8)
            a3.text(r['st']+(r['ed']-r['st'])/2, dv + ytop*0.05, f"{dv:.3f}", ha='center', fontsize=9, fontweight='bold')
    
    a3.set_ylabel('Orthologous dN/dS', fontweight='bold')
    a3.set_xlabel('Genomic Coordinates (bp)', fontweight='bold')
    a3.set_ylim(0, ytop); a3.set_xlim(0, max_l)
    a3.spines['top'].set_visible(False); a3.spines['right'].set_visible(False)
    plt.tight_layout(); plt.savefig(os.path.join(output_dir, "Fig10_Trinity_Landscape.png"), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "Fig10_Trinity_Landscape.pdf"), bbox_inches='tight'); plt.close()

# ==========================================
# ==========================================
# [Fig 11] 无监督机器学习聚类 (PCA + KMeans)
# Feature matrix: [dN, dS, piN, piS] — evolutionary dynamics, NOT genotype
# This is the key difference from traditional pop-gen PCA (plink on 0/1/2)
# ==========================================
def plot_ml_clusters(df_pop, output_dir, metadata_file=None):
    if not HAS_ML: return
    print("  -> Fig 11: 启动机器学习引擎，构建群体遗传学 2D & 3D 降维宇宙...")
    feats = ['mean_dN_vs_ref', 'mean_dS_vs_ref', 'piN', 'piS']
    df = df_pop.dropna(subset=feats).copy()

    if len(df) < 15:
        print("    -> [跳过] 样本数量不足 15 个，PCA 聚类意义不大。")
        return

    X = StandardScaler().fit_transform(df[feats])

    # === AI Auto-K: Silhouette Score maximization ===
    max_k = max(2, min(7, len(X) // 5))
    best_k, best_score = 2, -1.0
    for k in range(2, max_k + 1):
        labels_temp = KMeans(n_clusters=k, random_state=42).fit_predict(X)
        score = silhouette_score(X, labels_temp)
        if score > best_score:
            best_score, best_k = score, k
    print(f"    -> [AI 决策] 自动评估完成：数据在 K={best_k} 时轮廓系数最高 (Score: {best_score:.3f})，已锁定最优亚群数。")

    # === 3-component PCA (not just 2 — PC3 reveals hidden substructure) ===
    pca = PCA(n_components=3)
    p_res = pca.fit_transform(X)
    df['PC1'], df['PC2'], df['PC3'] = p_res[:, 0], p_res[:, 1], p_res[:, 2]

    df['ML_Cluster'] = [f"Evo-Class {i+1}" for i in KMeans(n_clusters=best_k, random_state=42).fit_predict(X)]

    # === Optional: external metadata coloring (geographic origin, host species, etc.) ===
    color_col = 'ML_Cluster'
    title_suffix = f"(Unsupervised K={best_k})"

    if metadata_file and os.path.exists(metadata_file):
        try:
            meta = pd.read_csv(metadata_file)
            df = pd.merge(df, meta, how='left')
            meta_label = meta.columns[1]
            df[meta_label] = df[meta_label].fillna('Unknown')
            color_col = meta_label
            title_suffix = f"(Colored by True {meta_label})"
            print(f"    -> [成功] 已挂载外部元数据 {meta_label} 用于映射真实流行病学分支。")
        except Exception as e:
            print(f"    -> [警告] 元数据融合失败: {e}，回退至 AI 自动聚类着色。")

    df.to_csv(os.path.join(output_dir, "Matrix_ML_Clusters_3D.csv"), index=False)
    palette = okabe_cycle(len(df[color_col].unique()))

    # ================= 2D PCA =================
    fig, ax2d = plt.subplots(figsize=(9, 7))
    sns.scatterplot(data=df, x='PC1', y='PC2', hue=color_col, palette=palette,
                    s=130, edgecolor='white', linewidth=1.2, alpha=0.85, ax=ax2d)

    # Dynamic headroom for frosted-glass info panel
    y_min, y_max = ax2d.get_ylim()
    ax2d.set_ylim(y_min, y_max + (y_max - y_min) * 0.28)

    ax2d.set_title(f'Fig 11a: PCA of Evolutionary Dynamics {title_suffix}', fontsize=14, fontweight='bold')
    ax2d.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontweight='bold')
    ax2d.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", fontweight='bold')
    # 图例标题不暴露内部列名/方法学标签 (旧版直接显示 "ML_Cluster", 投稿时属术语错位)
    legend_title = 'Evolutionary class' if color_col == 'ML_Cluster' else color_col
    ax2d.legend(title=legend_title, loc='upper left')

    # Cluster center profiling (piN, piS, dN)
    cc = df.groupby('ML_Cluster')[feats].mean()
    txt = "Cluster profiling (mean):\n"
    for cl in sorted(cc.index):
        txt += f"  {cl}: piN={cc.loc[cl,'piN']:.4f}, piS={cc.loc[cl,'piS']:.4f} | dN={cc.loc[cl,'mean_dN_vs_ref']:.4f}\n"

    ax2d.text(0.97, 0.98, txt.strip(), transform=ax2d.transAxes, fontsize=10, va='top', ha='right',
              bbox=dict(boxstyle='round,pad=0.6', facecolor='#f9f9f9', edgecolor='black', alpha=0.9))

    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Fig11a_PCA_2D.png"), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "Fig11a_PCA_2D.pdf"), bbox_inches='tight')
    plt.close()

    # ================= 3D PCA =================
    fig = plt.figure(figsize=(10, 8))
    ax3d = fig.add_subplot(111, projection='3d')

    unique_labels = df[color_col].unique()
    color_dict = {label: col for label, col in zip(unique_labels, palette)}

    for label in unique_labels:
        subset = df[df[color_col] == label]
        ax3d.scatter(subset['PC1'], subset['PC2'], subset['PC3'],
                     label=label, c=[color_dict[label]], s=80, edgecolors='w', alpha=0.9)

    ax3d.set_title(f'Fig 11b: 3D PCA Landscape {title_suffix}', fontsize=15, fontweight='bold')
    ax3d.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontweight='bold')
    ax3d.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", fontweight='bold')
    ax3d.set_zlabel(f"PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)", fontweight='bold')

    ax3d.view_init(elev=20, azim=45)
    ax3d.grid(True, linestyle=':', alpha=0.6)
    ax3d.legend(title=legend_title, bbox_to_anchor=(1.15, 0.8), loc='upper left')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "Fig11b_PCA_3D.png"), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "Fig11b_PCA_3D.pdf"), bbox_inches='tight')
    plt.close()

    print(f"    -> [完成] 2D + 3D PCA 双模式图表已生成。")
    print(f"    -> 论文模板: 'To prevent human bias, optimal K was determined by maximizing "
          f"the Silhouette Score. The algorithm converged on K={best_k}, revealing {best_k} "
          f"distinct evolutionary trajectories among viral isolates.'")

# ==========================================
# 主管道调度员
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="SNPGenie Mega-Miner (V13.2 终极整合与防错版)")
    parser.add_argument('-i', '--input', type=str, required=True, help='挂载核心 SNPGenie 输出源')
    parser.add_argument('-o', '--output', type=str, default='./SNPGenie_Omega_Results', help='结果导出域')
    parser.add_argument('-r', '--ref', type=str, default='', help='NCBI accession for coordinate download (required for Fig 10)')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    setup_plot_style()
    
    print("\n" + "="*80)
    print(" SNPGENIE SUPER OMEGA MINER[绝密战列舰级] ".center(80, "="))
    print("="*80 + "\n")

    try:
        # P1: 剥离假象加载数据
        df_pop, df_prod, df_site = parse_and_filter_data(args.input)
        df_site_filtered = apply_strict_isnv_filters(df_site)
        
        print("\n[*] ------------------------[神像级绘图管线全开] ------------------------")
        plot_vaf_spectrum(df_site_filtered, args.output)                   # 1
        calculate_adjusted_n_per_kb(df_site_filtered, df_prod, args.output)# 2
        plot_joint_dynamics(df_pop, args.output)                           # 3 & 4
        plot_gene_dnds_with_stats(df_prod, args.output)                    # 5
        perform_10000x_bootstrap_dnds(df_prod, args.output)                # 6
        # 7 & 8: Fig7 主用 QC 过滤后的多态位点，df_site 仅作兜底（仍排除不变位点）
        plot_signatures_and_hotspots(df_site_filtered, df_site, args.output)
        plot_dual_track_window(df_site, args.output)                       # 9
        if args.ref:
            plot_trinity_landscape(df_site, df_prod, args.output, ref_acc=args.ref) # 10
        else:
            print("  -> [跳过 Fig 10] 未提供 -r/--ref 参数，Trinity 全景观图需要 NCBI accession")
        plot_ml_clusters(df_pop, args.output)                              # 11

        print("\n"+"="*80)
        print(" 🎉 【大满贯闭环完成】 11张具备 Nature/Cell 系列碾压实力的学术底片集结完毕！")
        print(f" 打开检阅: {os.path.abspath(args.output)}")
        print("="*80+"\n")
        
    except Exception as e:
        print(f"\n[阻滞] 超算引擎触发未知死锁: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
