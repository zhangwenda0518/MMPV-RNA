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
    min_occ = max(2, int(df['Sample'].nunique() * 0.03))
    valid_vars = hotspots[hotspots['Sample_Count'] >= min_occ]['Variant_ID']
    onco_out = onco_matrix.loc[onco_matrix.index.intersection(valid_vars)]
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
    pal_imp = {'HIGH': '#d73027', 'MODERATE': '#fc8d59', 'LOW': '#91bfdb', 'MODIFIER': '#4575b4', 'UNKNOWN': 'gray'}
    
    sns.scatterplot(data=pos_stats, x='POS', y='Samples', size='MeanAF', sizes=(20,250), hue='IMPACT', palette=pal_imp, alpha=0.8)
    plt.title("Figure 1: Genomic Mutational Landscape of Viral Quasispecies", fontweight='bold')
    plt.xlabel("Genomic Position (Nucleotide POS)")
    plt.ylabel("Number of Affected Hosts")
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Figure_1_Manhattan_Mut_Landscape.pdf"), dpi=300)
    plt.savefig(os.path.join(outdir, "Figure_1_Manhattan_Mut_Landscape.png"), dpi=300)
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
            c_map = {'Consensus(>80%)': '#af8dc3', 'iSNV(5-80%)': '#7fbf7b'}
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
            cm = sns.clustermap(heatmap_df, cmap='flare', figsize=(12, 10), metric='euclidean', method='ward',
                                xticklabels=True, yticklabels=False, cbar_kws={'label':'Allele Frequency'})
            cm.fig.suptitle("Figure 4: Unsupervised Lineage Clustering via Shared Missense Mutations", fontweight='bold', y=1.02)
            plt.savefig(os.path.join(outdir, "Figure_4_Lineage_Clustermap.pdf"), dpi=300)
            plt.savefig(os.path.join(outdir, "Figure_4_Lineage_Clustermap.png"), dpi=300)
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

    r_code = f"""
options(warn=-1)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(grid))

onco_file <- "{onco_csv}"
onco_mat <- as.matrix(read.csv(onco_file, row.names=1, check.names=FALSE))

col_map <- c(
    "Fixed_Missense" = "#E31A1C",   "Major_Missense" = "#FC4E2A",   "Minor_Missense" = "#FD8D3C",
    "Fixed_Synonymous" = "#1F78B4", "Major_Synonymous"="#41B6C4",   "Minor_Synonymous" = "#A1DAB4",
    "Fixed_Regulatory" = "#33A02C", "Major_Regulatory"="#74C476",   "Minor_Regulatory" = "#C7E9C0"
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
    Minor_Regulatory = function(x,y,w,h) {{ grid.rect(x,y,w,h*0.4, gp=gpar(fill=col_map["Minor_Regulatory"], col=NA)) }}
)

out_pdf <- "{out_pdf}"
out_png <- "{out_png}"

ht = oncoPrint(onco_mat,
               alter_fun = alter_fun, col = col_map,
               remove_empty_columns = TRUE, remove_empty_rows = TRUE,
               show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 9),
               pct_side = "right", row_names_side = "left",
               column_title = "Figure 5: Mutational Waterfall Landscape of Viral Isolates",
               heatmap_legend_param = list(title = "Fixation State", at=names(col_map)))

pdf(out_pdf, width=16, height=10)
draw(ht)
invisible(dev.off())

png(out_png, width=16, height=10, units="in", res=300)
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
