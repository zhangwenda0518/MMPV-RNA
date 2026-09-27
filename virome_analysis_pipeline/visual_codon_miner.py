#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
@Tool:          Auto Codon Miner Pipeline (DRHIP Integration)
@Description:   Automatically parses drhip results, identifies positive 
                selection sites, queries source alignments, and plots results.
@Author:        World-Class Bioinformatics AI
@Requirements:  biopython, pandas, matplotlib, seaborn
=============================================================================
"""

import os
import sys
import argparse
import csv
import re
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

def get_args():
    parser = argparse.ArgumentParser(
        description="🧬 Auto Codon Miner: 自动读取 drhip 结果，全自动溯源比对文件并出图",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
🔥 推荐使用 (Pipeline 模式):
  python auto_codon_miner.py --drhip drhip/combined_sites.csv --clndir hyphy/CLN/ -o Final_Results

🛠️ 手动使用 (单体模式):
  python auto_codon_miner.py -i hyphy/CLN/ORxxx.part_L-nodups.fasta -r ORxxx.1 -s 1047 -o Results
        """
    )
    
    # 自动管线参数
    parser.add_argument("--drhip", type=str, default=None,
                        help="drhip 结果文件路径 (如: drhip/combined_sites.csv)")
    parser.add_argument("--clndir", type=str, default="hyphy/CLN",
                        help="密码子比对文件所在目录 (默认: hyphy/CLN)")
    
    # 手动模式参数
    parser.add_argument("-i", "--input", type=str, help="输入比对文件 (FASTA)")
    parser.add_argument("-r", "--ref", type=str, help="参考序列 ID")
    parser.add_argument("-s", "--sites", nargs='+', type=int, help="氨基酸位点")
    
    # 公共输出参数
    parser.add_argument("-o", "--outdir", default="./Publication_Plots", type=str,
                        help="输出结果目录 (默认: ./Publication_Plots)")
    
    return parser.parse_args()

def parse_drhip_csv(drhip_file, cln_dir):
    """
    智能解析 DRHIP 文件
    策略：无需提前知道表头结构。逐行读取，只要行内包含 'positive' 字样，
          就提取第 0 列 (Gene Prefix) 和第 1 列 (Site 号)。自动推断关联路径。
    """
    print(f"[*] 正在解析 DRHIP 综合统计报告: {drhip_file}")
    if not os.path.exists(drhip_file):
        raise FileNotFoundError(f"[!] 找不到 DRHIP 文件: {drhip_file}")

    tasks = {}
    positive_count = 0

    with open(drhip_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row: continue # 跳过空行
            
            # 判断该行数据是否命中 'positive' 选择标记
            if "positive" in [str(cell).strip().lower() for cell in row]:
                gene_prefix = row[0] # 例如: OR489165.1.cds-noStopCodons.part_L
                try:
                    site = int(row[1])   # 例如: 1047
                except ValueError:
                    continue # 避免表头干扰
                
                # 智能识别参考基因组 ID (兼容 '.cds' 前缀分隔)
                ref_id = gene_prefix.split('.cds')[0] 
                
                # 构建关联的比对文件路径 (适配了你在 tree 中展示的 -nodups.fasta 格式)
                fasta_path = os.path.join(cln_dir, f"{gene_prefix}-nodups.fasta")
                
                if fasta_path not in tasks:
                    tasks[fasta_path] = {'ref': ref_id, 'sites': set()}
                
                tasks[fasta_path]['sites'].add(site)
                positive_count += 1
                
    print(f"    [+] 成功解析！共发现 {positive_count} 个处于正选择的目标位点，归属 {len(tasks)} 个基因文件。")
    return tasks

def translate_codon(codon_str):
    codon_str = codon_str.upper()
    if "-" in codon_str: return "GAP (-)"
    if not re.match(r'^[ATCGU]{3}$', codon_str): return "AMBIGUOUS"
    try: return str(Seq(codon_str).translate(table=1))
    except Exception: return "UNKNOWN"

def map_reference_codon_to_msa(alignment, ref_id, target_site):
    """
    针对 SRA 组学挖掘结果量身定制的坐标系锚定。
    """
    # 通过 ref_id 精确匹配参考序列，不假设 alignment[0] 一定是参考
    ref_record = None
    for record in alignment:
        if ref_id in record.id:
            ref_record = record
            break
    if ref_record is None:
        raise ValueError(f"Reference sequence '{ref_id}' not found in alignment")

    ref_seq = str(ref_record.seq).upper()
    biol_codon_pos = 0

    for i in range(0, len(ref_seq), 3):
        current_codon = ref_seq[i:i+3]

        # 跳过纯 Gap，累加实际可见的实体密码子
        if current_codon != "---":
            biol_codon_pos += 1
            if biol_codon_pos == target_site:
                return i

    raise IndexError(f"目标位点 {target_site} 越界 (参考尺度的实际有效密码子数: {biol_codon_pos})")


def extract_site_data(alignment, nuc_start_idx):
    records = [{'Amino_Acid': translate_codon(str(r.seq)[nuc_start_idx:nuc_start_idx+3]),
                'Codon': str(r.seq)[nuc_start_idx:nuc_start_idx+3].upper()} 
               for r in alignment]
    df = pd.DataFrame(records)
    total_samples = len(df)
    if total_samples == 0: return None
    
    summary = df.groupby(['Amino_Acid', 'Codon']).size().reset_index(name='Count')
    summary['Frequency (%)'] = (summary['Count'] / total_samples) * 100
    return summary.sort_values(by='Count', ascending=False).reset_index(drop=True)

def plot_publication_figure(df, file_basename, ref_id, site, outdir):
    df['Label'] = df['Amino_Acid'] + "\n(" + df['Codon'] + ")"
    sns.set_theme(style="whitegrid", font_scale=1.1)
    
    # 动态适应柱子的数量，保证图表不拥挤
    plt.figure(figsize=(max(6, len(df)*1.5), 6)) 
    
    ax = sns.barplot(x='Label', y='Frequency (%)', hue='Amino_Acid', data=df, dodge=False, palette="Set2")
    param_gene = file_basename.split('.part_')[-1].split('-')[0] # 智能提取短基因名 (如 L, G)
    
    plt.title(f"Adaptive Evolution at Site {site} of {param_gene} Protein\n(Reference: {ref_id})", 
              fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Amino Acid Substitution (Codon)', fontsize=12, fontweight='bold')
    plt.ylabel('Population Frequency (%)', fontsize=12, fontweight='bold')
    upper = df['Frequency (%)'].max() * 1.15; ax.set_ylim(0, min(100, upper) if upper <= 100 else 100) 
    
    for p, count, freq in zip(ax.patches, df['Count'], df['Frequency (%)']):
        ax.annotate(f"n={int(count)}\n{freq:.1f}%", 
                    (p.get_x() + p.get_width() / 2., p.get_height()), 
                    ha='center', va='bottom', xytext=(0, 5), textcoords='offset points',
                    fontsize=10, color='black', fontweight='bold')
    if ax.get_legend(): ax.get_legend().remove()
    plt.tight_layout()
    
    plot_file = os.path.join(outdir, f"{param_gene}_Protein_Site_{site}_Evolution.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

def process_tasks(tasks_dict, outdir):
    """流水线执行引擎"""
    for fasta_file, meta in tasks_dict.items():
        ref_id = meta['ref']
        sites = sorted(list(meta['sites'])) # 按位点顺序处理
        
        file_basename = os.path.basename(fasta_file)
        print(f"\n[🚀] 加载任务: {file_basename}")
        
        if not os.path.exists(fasta_file):
            print(f"    [!] 找不到对应的比对文件: {fasta_file}，已跳过。")
            continue
            
        try:
            alignment = AlignIO.read(fasta_file, "fasta")
            print(f"    [+] 加载成功，包含序列: {len(alignment)} 条")
        except Exception as e:
            print(f"    [!] 解析 FASTA 失败: {e}"); continue
            
        for site in sites:
            print(f"      -> 挖掘深入位点 : {site}...")
            try:
                # 核心逻辑
                idx = map_reference_codon_to_msa(alignment, ref_id, site)
                df = extract_site_data(alignment, idx)
                if df is None: continue
                
                # 储存与绘图
                short_gene = file_basename.split('.part_')[-1].split('-')[0]
                csv_file = os.path.join(outdir, f"{short_gene}_Protein_Site_{site}.csv")
                df.to_csv(csv_file, index=False)
                plot_publication_figure(df, file_basename, ref_id, site, outdir)
            except Exception as e:
                print(f"      [!] 挖掘失败: {e}")

def main():
    args = get_args()
    print("\n" + "═"*70)
    print("🧬 Auto Codon Miner - DRHIP End-To-End Pipeline 🧬")
    print("═"*70)
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # 模式选择器
    if args.drhip:
        tasks = parse_drhip_csv(args.drhip, args.clndir)
        process_tasks(tasks, args.outdir)
    elif args.input and args.ref and args.sites:
        tasks = {args.input: {'ref': args.ref, 'sites': set(args.sites)}}
        process_tasks(tasks, args.outdir)
    else:
        print("[!] 参数不足。请提供 --drhip 或完整的 (-i, -r, -s) 参数。使用 -h 查看帮助。")
        sys.exit(1)
        
    print("\n" + "═"*70)
    print(f"🎉 全部 Pipeline 运行完毕！\n请打开 '{args.outdir}' 目录核收高水平多态性柱状图及数据表！")
    print("═"*70 + "\n")

if __name__ == "__main__":
    main()
