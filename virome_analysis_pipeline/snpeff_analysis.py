#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
【计算病毒群体遗传学终极一体化枢纽】(The Grand Unified Pipeline - Ultimate Edition)
作者：顶级生物信息开发专家
描述：
 1. 兼容所有病毒/类病毒，自动剥离 NCBI 编号（智能正则）。
 2. 具备生物学防护，遇到无编码区(Viroid)智能跳过相关蛋白分析图。
 3. 自动生成并【在后台一键静默拉起执行】R 语言 OncoPrint 瀑布图脚本。
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

import argparse
import os
import glob
import re
import subprocess
import pandas as pd
import numpy as np

# 强制服务器无头模式绘图
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", font_scale=1.1)
plt.rcParams['pdf.fonttype'] = 42

# ============================================================
# 统一配色 (SCI 发表级 / Okabe-Ito 色盲友好)
#   与 virome_analysis_pipeline/virus_metadata_plot.py 的 OKABE_ITO 同序同成员;
#   R 模板 (ComplexHeatmap OncoPrint) 的 col_map 同步改用 Okabe 色相。
#   2026-10-11: 原文件用 ColorBrewer Reds/Blues 与 RdYlBu 派生色, 与管线其余图不一致。
# ============================================================
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


def okabe_shade(hex_color, white_frac):
    """把 Okabe 基色按 white_frac 向白色混合, 生成同色相的浅色 (用于固定/主要/次要分级)。
    R 模板的 OncoPrint 需要 12 个状态色, Okabe 只有 7 个基色, 故按
    「功能类别定色相 + 频率定深浅」展开, 与 alter_fun 的矩形高度编码互相印证。"""
    h = str(hex_color).lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    f = max(0.0, min(1.0, float(white_frac)))
    mix = lambda c: int(round(c + (255 - c) * f))
    return "#%02X%02X%02X" % (mix(r), mix(g), mix(b))


# clustermap 的行列树状图分支颜色取自 matplotlib 默认色环 'C0'..'C9',
# 必须在 sns.set_theme 之后再接管, 否则会退回 tab10。
# 2026-10-11: 改用 okabe_branch_cycle (排除浅黄 #F0E442 / 浅蓝 #56B4E9), 提升细线分支可读性。
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=okabe_branch_cycle(10))

# ==========================================
# 核心质控参数与辅助函数
# ==========================================
MIN_DP = 50     # 最小测序深度阈值 (与 Stage6 --post_min_dp 对齐)
MIN_AF = 0.05   # 最小等位基因频率阈值

def get_titv(ref, alt):
    if len(ref) != 1 or len(alt) != 1: return "Indel/Complex"
    if set([ref, alt]) in [{"A", "G"}, {"C", "T"}]: return "Ti"
    return "Tv"

def determine_oncoprint_state(af, impact):
    if af < MIN_AF: return ""
    freq = "Minor" if af < 0.50 else ("Major" if af < 0.95 else "Fixed")
    if impact == 'HIGH': imp = "Truncating"
    elif impact == 'MODERATE': imp = "Missense"
    elif impact == 'LOW': imp = "Synonymous"
    else: imp = "Regulatory"
    return f"{freq}_{imp}"

# ==========================================
# Phase 1: 极度稳健的多维 VCF 数据张量提取
# ==========================================
def orchestrate_data_mining(input_dir):
    print("\n[Phase 1/4] 🚀 开始跨文件数据挖掘与质控清洗...")
    
    # 动态匹配任何病毒文件目录，不局限于特定的 OR489165.1
    vcf_pattern = os.path.join(input_dir, "*", "*.ann.vcf")
    vcf_files = glob.glob(vcf_pattern)
    
    if not vcf_files:
        raise FileNotFoundError(f"在此目录未发现符合特征的 VCF 文件: {vcf_pattern}")
    
    all_records =[]
    
    for vcf in vcf_files:
        # 智能剥离后缀，兼容类似 _OR489165.1 或 _NC_002030.1 的命名
        folder_name = os.path.basename(os.path.dirname(vcf))
        sample = re.sub(r'_[A-Za-z]{2,}_?\d+\.\d+$', '', folder_name)
        
        with open(vcf, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'): continue
                cols = line.strip().split('\t')
                
                # 质量过滤
                if len(cols) < 8 or cols[6] not in ('PASS', '.'): continue
                
                chrom, pos, ref, alt, info_str = cols[0], cols[1], cols[3], cols[4], cols[7]
                info = dict(item.split('=', 1) if '=' in item else (item, True) for item in info_str.split(';'))
                
                try:
                    dp = int(info.get('DP', 0))
                    af = float(info.get('AF', '0.0').split(',')[0])
                except ValueError: continue
                if dp < MIN_DP or af < MIN_AF: continue
                
                # SnpEff 注释拆解
                ann_str = info.get('ANN', '')
                if ann_str:
                    anns = ann_str.split(',')[0].split('|')
                    effect = anns[1] if len(anns)>1 else 'UNKNOWN'
                    impact = anns[2] if len(anns)>2 else 'UNKNOWN'
                    gene = anns[3] if len(anns)>3 else 'Intergenic'
                    aa_change = anns[10] if len(anns)>10 else ''
                else:
                    effect, impact, gene, aa_change = ('NONE', 'NONE', 'Intergenic', '')
                
                # 变异ID
                mut_marker = aa_change if aa_change else f"{ref}>{alt}"
                if 'upstream' in effect or 'intergenic' in effect:
                    vid = f"Reg_{pos}_{mut_marker}"
                else:
                    vid = f"{gene}_{pos}_{mut_marker}"
                    
                all_records.append({
                    'Sample': sample, 'CHROM': chrom, 'POS': int(pos),
                    'Variant_ID': vid, 'REF': ref, 'ALT': alt,
                    'TiTv': get_titv(ref, alt), 'GENE': gene,
                    'EFFECT': effect, 'IMPACT': impact,
                    'DP': dp, 'AF': af,
                    'Quasispecies': 'Consensus(>80%)' if af >= 0.8 else 'iSNV(5-80%)',
                    'Onco_State': determine_oncoprint_state(af, impact)
                })
                
    return pd.DataFrame(all_records)

# ==========================================
# Phase 2: 六大矩阵统计生成
# ==========================================
def export_statistical_matrices(df, outdir):
    print("[Phase 2/4] 🧮 计算与降维，正在融合输出六大核心发文数据矩阵...")
    
    df.to_csv(os.path.join(outdir, "Matrix_01_Global_Melted_Variants.csv"), index=False)
    
    s_stats = df.pivot_table(index='Sample', columns='IMPACT', values='POS', aggfunc='count', fill_value=0)
    s_stats['Total_Mutations'] = s_stats.sum(axis=1)
    s_stats.to_csv(os.path.join(outdir, "Matrix_02_PerSample_Burden.csv"))
    
    hotspots = df.groupby(['Variant_ID', 'GENE', 'IMPACT', 'POS']).agg(
        Sample_Count=('Sample', 'nunique'), Mean_AF=('AF', 'mean')
    ).reset_index().sort_values('Sample_Count', ascending=False)
    hotspots['Sample_Frequency(%)'] = (hotspots['Sample_Count'] / df['Sample'].nunique()) * 100
    hotspots.to_csv(os.path.join(outdir, "Matrix_03_Population_Hotspots.csv"), index=False)
    
    gene_stats = df.groupby('GENE').agg(
        Total_Observed=('Variant_ID', 'count'),
        Missense_Events=('IMPACT', lambda x: (x=='MODERATE').sum()),
        Synonymous_Events=('IMPACT', lambda x: (x=='LOW').sum())
    )
    gene_stats['Miss_Syn_Ratio'] = gene_stats['Missense_Events'] / (gene_stats['Synonymous_Events'] + 1)
    gene_stats.to_csv(os.path.join(outdir, "Matrix_04_Gene_Selection_Pressure.csv"))
    
    af_matrix = df.pivot_table(index='Variant_ID', columns='Sample', values='AF', fill_value=0.0)
    af_matrix.to_csv(os.path.join(outdir, "Matrix_05_Continuous_AF_Matrix.csv"))
    
    onco_matrix = df.pivot_table(index='Variant_ID', columns='Sample', values='Onco_State',
        aggfunc=lambda x: ';'.join(set(str(v) for v in x if pd.notna(v) and str(v).strip() != "")),
        fill_value="")
    # P2-17: 单样本时 nunique=1 -> max(2, 0)=2 会要求同一变异出现 >=2 个样本,
    # Matrix_06 必空 -> R oncoPrint 静默失败; n<2 时按 min_occ=1 保留该样本全部变异
    n_samples = df['Sample'].nunique()
    min_occ = 1 if n_samples < 2 else max(2, int(n_samples * 0.03))
    valid_vars = hotspots[hotspots['Sample_Count'] >= min_occ]['Variant_ID']
    onco_out = onco_matrix.loc[onco_matrix.index.intersection(valid_vars)]
    if onco_out.empty:
        print(f"   -> ⚠️ Matrix_06 (OncoPrint) 为空: 无变异达到 min_occ={min_occ} (样本数={n_samples})。")
    else:
        print(f"   -> 📄 Matrix_06 (OncoPrint): {len(onco_out)} 变异 x {onco_out.shape[1]} 样本 (min_occ={min_occ})。")
    onco_out.to_csv(os.path.join(outdir, "Matrix_06_Ultimate_OncoPrint_State.csv"))
    
    return hotspots

# ==========================================
# Phase 3: Python 本土的四大顶级图表渲染 (带类病毒智能防护)
# ==========================================
def render_core_figures(df, hotspots, outdir):
    print("[Phase 3/4] 🎨 正在调用内核渲染生信图像...")
    
    # --- 图 1：突变全景曼哈顿图 ---
    plt.figure(figsize=(14, 6))
    pos_stats = df.groupby(['POS', 'IMPACT']).agg(Samples=('Sample','nunique'), MeanAF=('AF','mean')).reset_index()
    pal_imp = {'HIGH': VERM, 'MODERATE': ORANGE, 'LOW': SKY, 'MODIFIER': BLUE, 'UNKNOWN': GREY}
    
    sns.scatterplot(data=pos_stats, x='POS', y='Samples', size='MeanAF', sizes=(20,250), hue='IMPACT', palette=pal_imp, alpha=0.8)
    plt.title("Figure 1: Genomic Mutational Landscape of Viral Quasispecies", fontweight='bold')
    plt.xlabel("Genomic Position (Nucleotide POS)")
    plt.ylabel("Number of Affected Hosts")
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Figure_1_Manhattan_Mut_Landscape.pdf"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(outdir, "Figure_1_Manhattan_Mut_Landscape.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- 图 2：突变效应载荷对比 (添加类病毒判定防空图) ---
    coding_df = df[df['IMPACT'].isin(['LOW', 'MODERATE', 'HIGH'])].copy()
    if not coding_df.empty:
        gene_summary = coding_df.groupby(['GENE', 'IMPACT']).size().unstack(fill_value=0)
        if not gene_summary.empty:
            fig, ax = plt.subplots(figsize=(10, 6))
            colors =[pal_imp.get(c, '#000') for c in gene_summary.columns]
            gene_summary.plot(kind='bar', stacked=True, ax=ax, color=colors, edgecolor='black')
            plt.title("Figure 2: Evolutionary Payload per Viral Gene", fontweight='bold')
            plt.ylabel("Cumulative Discovered Variations")
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "Figure_2_Gene_Payload.pdf"), dpi=300)
            plt.savefig(os.path.join(outdir, "Figure_2_Gene_Payload.png"), dpi=300)
            plt.close()
    else:
        print("   -> ⚠️ 无蛋白编码区 (类病毒特性)，自动跳过 Figure 2 (基因突变图)。")

    # --- 图 3：准种动力学多样性 ---
    if not df.empty:
        div_stats = df.groupby(['Sample', 'Quasispecies']).size().unstack(fill_value=0)
        div_stats['Tot'] = div_stats.sum(axis=1)
        div_stats = div_stats[div_stats['Tot'] > 2].sort_values('Tot', ascending=False).head(50).drop(columns=['Tot'])
        if not div_stats.empty:
            fig, ax = plt.subplots(figsize=(16, 6))
            c_map = {'Consensus(>80%)': PURPLE, 'iSNV(5-80%)': GREEN}
            colors_used =[c_map.get(c, 'grey') for c in div_stats.columns]
            div_stats.plot(kind='bar', stacked=True, color=colors_used, ax=ax, edgecolor='white')
            plt.title("Figure 3: Intra-host Dynamics (Quasispecies iSNVs vs Consensus)", fontweight='bold')
            plt.ylabel("Variant Count")
            plt.xticks(rotation=90, fontsize=8)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "Figure_3_IntraHost_Diversity.pdf"), dpi=300)
            plt.savefig(os.path.join(outdir, "Figure_3_IntraHost_Diversity.png"), dpi=300)
            plt.close()

    # --- 图 4：重点错义聚类网络 ---
    missense_df = df[df['IMPACT'] == 'MODERATE']
    valid_mrk = hotspots[(hotspots['IMPACT'] == 'MODERATE') & (hotspots['Sample_Frequency(%)'] >= 5)]['Variant_ID']
    if not missense_df.empty and len(valid_mrk) >= 2:
        heatmap_df = missense_df[missense_df['Variant_ID'].isin(valid_mrk)].pivot_table(index='Sample', columns='Variant_ID', values='AF', fill_value=0)
        heatmap_df = heatmap_df.loc[(heatmap_df != 0).any(axis=1)] 
        if heatmap_df.shape[0] >= 2 and heatmap_df.shape[1] >= 2:
            n_col = heatmap_df.shape[1]
            fig_w = max(12.0, min(40.0, 4.0 + n_col * 0.22))
            cm = sns.clustermap(heatmap_df, cmap=CMAP_SEQ, figsize=(fig_w, 10), metric='euclidean', method='ward',
                                xticklabels=True, yticklabels=False, dendrogram_ratio=0.10,
                                cbar_kws={'label':'Allele Frequency'})
            cm.ax_heatmap.set_xlabel("Shared missense variant", fontweight='bold')
            cm.ax_heatmap.set_ylabel("Sample", fontweight='bold')
            # 列数多时逐列标签会糊成黑带 -> 字号随列数收缩并旋转 90 度
            plt.setp(cm.ax_heatmap.get_xticklabels(), rotation=90,
                     fontsize=max(4, min(9, int(300.0 / max(n_col, 1)))))
            cm.fig.suptitle("Figure 4: Unsupervised Lineage Clustering via Shared Missense Mutations", fontweight='bold', y=1.02)
            # 必须用 clustermap 自带的 fig 保存并带 bbox_inches:
            # 旧版 plt.savefig(无 bbox_inches) -> suptitle 被裁掉
            cm.fig.savefig(os.path.join(outdir, "Figure_4_Lineage_Clustermap.pdf"), dpi=300, bbox_inches='tight')
            cm.fig.savefig(os.path.join(outdir, "Figure_4_Lineage_Clustermap.png"), dpi=300, bbox_inches='tight')
            plt.close()
    else:
        print("   -> ⚠️ 无足够多错义突变参与分组，自动跳过 Figure 4 (错义分支聚类图)。")

# ==========================================
# Phase 4: 在输出目录自动生成 R 脚本并通过系统桥接运行
# ==========================================
def inject_and_run_r_script(outdir):
    print("[Phase 4/4] 🧬 正在自动生成并后台拉起 R 语言，绘制顶区级瀑布图...")
    outdir_safe = outdir.replace('\\', '/')
    onco_csv = f"{outdir_safe}/Matrix_06_Ultimate_OncoPrint_State.csv"
    out_pdf = f"{outdir_safe}/Figure_5_Ultimate_Quasispecies_OncoPrint.pdf"
    out_png = f"{outdir_safe}/Figure_5_Ultimate_Quasispecies_OncoPrint.png"

    # P2-17: Matrix_06 为空 (或缺失) 时 R 的 read.csv/oncoPrint 必然失败且原先静默:
    # 这里显式跳过并说明原因, 不再拉起注定失败的 R 进程
    try:
        if not os.path.isfile(onco_csv) or pd.read_csv(onco_csv, index_col=0).empty:
            print("   -> ⚠️ 跳过 Figure 5 (OncoPrint): Matrix_06 为空 (无变异达到 min_occ 阈值)。")
            return
    except Exception as e:
        print(f"   -> ⚠️ 跳过 Figure 5 (OncoPrint): Matrix_06 读取失败 ({e})。")
        return

    # OncoPrint 12 状态色 (Okabe 色相): 功能类别定色相, 频率定深浅 (Fixed 最深 -> Minor 最浅),
    # 与 alter_fun 的矩形高度编码 (h / 0.75h / 0.4h) 互相印证。
    _onco_base = {"Missense": VERM, "Synonymous": BLUE,
                  "Regulatory": GREEN, "Truncating": BLACK}
    _onco_state_colors = {}
    for _imp, _base in _onco_base.items():
        for _freq, _frac in (("Fixed", 0.0), ("Major", 0.30), ("Minor", 0.55)):
            _onco_state_colors[f"{_freq}_{_imp}"] = okabe_shade(_base, _frac)
    _col_map_r = ",\n".join(f'    "{k}" = "{v}"' for k, v in _onco_state_colors.items())

    r_code = f"""
options(warn=-1)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(grid))

onco_file <- "{onco_csv}"
onco_mat <- as.matrix(read.csv(onco_file, row.names=1, check.names=FALSE))
# 过绘守卫: 行(变异) 与 列(样本) 双向封顶; 旧版只封行, 200+ 样本 x 200+ 变异时格子全糊
if (nrow(onco_mat) > 120) onco_mat <- onco_mat[1:120, , drop=FALSE]
if (ncol(onco_mat) > 200) onco_mat <- onco_mat[, 1:200, drop=FALSE]
row_fs <- max(5, min(10, round(900 / max(nrow(onco_mat), 1))))
fig_h <- max(6, min(18, 3 + nrow(onco_mat) * 0.06))
fig_w <- max(10, min(30, 4 + ncol(onco_mat) * 0.10))

col_map <- c(
{_col_map_r}
)

alter_fun = list(
    background = function(x,y,w,h) {{ grid.rect(x,y,w-unit(1,"pt"),h-unit(1,"pt"), gp=gpar(fill="#F5F5F5",col=NA)) }},
    Fixed_Missense = function(x,y,w,h) {{ grid.rect(x,y,w,h, gp=gpar(fill=col_map["Fixed_Missense"], col=NA)) }},
    Major_Missense = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.75, gp=gpar(fill=col_map["Major_Missense"], col=NA)) }},
    Minor_Missense = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.4, gp=gpar(fill=col_map["Minor_Missense"], col=NA)) }},
    Fixed_Synonymous = function(x,y,w,h) {{ grid.rect(x,y,w,h, gp=gpar(fill=col_map["Fixed_Synonymous"], col=NA)) }},
    Major_Synonymous = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.75, gp=gpar(fill=col_map["Major_Synonymous"], col=NA)) }},
    Minor_Synonymous = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.4, gp=gpar(fill=col_map["Minor_Synonymous"], col=NA)) }},
    Fixed_Regulatory = function(x,y,w,h) {{ grid.rect(x,y,w,h, gp=gpar(fill=col_map["Fixed_Regulatory"], col=NA)) }},
    Major_Regulatory = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.75, gp=gpar(fill=col_map["Major_Regulatory"], col=NA)) }},
    Minor_Regulatory = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.4, gp=gpar(fill=col_map["Minor_Regulatory"], col=NA)) }},
    Fixed_Truncating = function(x,y,w,h) {{ grid.rect(x,y,w,h, gp=gpar(fill=col_map["Fixed_Truncating"], col=NA)) }},
    Major_Truncating = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.75, gp=gpar(fill=col_map["Major_Truncating"], col=NA)) }},
    Minor_Truncating = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.4, gp=gpar(fill=col_map["Minor_Truncating"], col=NA)) }}
)

out_pdf <- "{out_pdf}"
out_png <- "{out_png}"

ht = oncoPrint(onco_mat,
               alter_fun = alter_fun, col = col_map,
               remove_empty_columns = TRUE, remove_empty_rows = TRUE,
               show_column_names = FALSE,
               row_names_gp = gpar(fontsize = row_fs),
               pct_side = "right", row_names_side = "left",
               column_title = "Figure 5: Mutational Waterfall Landscape of Viral Isolates",
               heatmap_legend_param = list(title = "Fixation State", at=names(col_map)))

pdf(out_pdf, width=fig_w, height=fig_h)
draw(ht)
invisible(dev.off())

png(out_png, width=fig_w, height=fig_h, units="in", res=300)
draw(ht)
invisible(dev.off())
"""

    r_file = os.path.join(outdir, "Execute_Figure_5_R_ComplexHeatmap.R")
    with open(r_file, "w", encoding="utf-8") as f:
        f.write(r_code.strip())

    log_file = os.path.join(outdir, "Rscript_OncoPrint.log")
    try:
        with open(log_file, "w") as rf_log:
            subprocess.run(["Rscript", r_file], check=True, stdout=rf_log, stderr=subprocess.STDOUT)
        print("   -> ✅ 跨语言联用成功！R引擎所绘原生阶梯瀑布 OncoPrint (图5) 已自动落盘！")
    except FileNotFoundError:
        print("   -> ⚠️ 未能在环境中检测到 Rscript，跳过 R语言自动化。")
    except subprocess.CalledProcessError:
        print(f"   -> ⚠️ R 执行失败, 图 5 (OncoPrint) 绘制失败。")
        print(f"   -> 💡 日志: {log_file} (检查是否缺少 ComplexHeatmap 包)")


# ==========================================
# 顶层业务调度执行
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="终极群体病毒发生引擎 - 纯内生自动化整合版")
    parser.add_argument("--miner", required=True, help="输入路径 (包含待挖数据的父级目录)")
    parser.add_argument("--outdir", required=True, help="统一结果存放路径")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    df_master = orchestrate_data_mining(args.miner)
    if df_master.empty:
        print("❌ 核心管线切断：未发现可用突变。")
        exit()
        
    hotspots_df = export_statistical_matrices(df_master, args.outdir)
    render_core_figures(df_master, hotspots_df, args.outdir)
    
    # Python 最后一公里直通 R
    inject_and_run_r_script(args.outdir)
    
    print("\n" + "♛"*40)
    print("      核心分析制图大满贯 全部完成！      ")
    print("♛"*40)
    print(f"➜ 您的所有数据明细矩阵以及5张终极科研图纸已归集至： {os.path.abspath(args.outdir)}\n")
