#!/usr/bin/env python3
import numpy as np
import pandas as pd

# ====== 修复服务器环境下的 xcb 报错 ======
import matplotlib
matplotlib.use('Agg') # 强制使用无头渲染引擎，必须放在导入 pyplot 之前
# ==========================================

import matplotlib.pyplot as plt
import argparse
import os
import gzip
import textwrap
import re

def smooth(y, box_pts):
    """对深度序列进行平滑处理"""
    if len(y) < box_pts:
        return y
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def format_number(value):
    """格式化大数字，添加千位分隔符"""
    return f"{value:,.0f}" if isinstance(value, (int, float)) else value

def create_info_text(v_stat):
    """提取 summary 中的高阶定量与注释信息"""
    info =[
        f"Taxonomy: {textwrap.shorten(str(v_stat.get('taxonomy', 'Unknown')), width=45, placeholder='...')}",
        f"Accession: {v_stat['Virus']}",
        f"Length: {format_number(v_stat['Length'])} bp",
        f"Coverage: {v_stat['Coverage(%)']:.2f}%",
        f"Mean Depth: {v_stat['MeanDepth']:.2f}x",
        "-"*25,
        f"Mapped Reads: {format_number(v_stat['MappedReads'])}",
        f"RPM: {v_stat.get('RPM', 0):.2f}",
        f"FPKM: {v_stat.get('FPKM', 0):.2f}",
        f"TPM: {v_stat.get('TPM', 0):.2f}"
    ]
    return "\n".join(info)

def main():
    parser = argparse.ArgumentParser(description='批量绘制带分类学注释和多维定量指标的全基因组深度分布图')
    parser.add_argument('-d', '--depth_dir', required=True, help='包含 pandepth 输出单碱基深度文件 (.SiteDepth.gz) 的统计目录 (通常为 stat 文件夹)')
    parser.add_argument('-m', '--summary', required=True, help='经过最优筛选的 summary 文件 (all_viruses.best.summary.tsv)')
    parser.add_argument('-n', '--sample', required=False, help='[可选] 指定只绘制单个样本。若不指定，将自动绘制表中的所有样本')
    parser.add_argument('-o', '--output', required=True, help='图片输出的基础目录 (例如 plots/)')
    parser.add_argument('-w', '--window', type=int, default=100, help='平滑窗口大小 (默认: 100)')
    parser.add_argument('-f', '--fontsize', type=float, default=9, help='注释框字体大小 (默认: 9)')
    args = parser.parse_args()

    print(f"📂 正在读取核心汇总文件: {args.summary}")
    try:
        # 智能兼容分隔符
        sep = ',' if args.summary.endswith('.csv') else '\t'
        summary_df = pd.read_csv(args.summary, sep=sep)
    except Exception as e:
        print(f"❌ 读取 summary 文件失败: {e}")
        return

    # 判断是否为批量模式
    if args.sample:
        samples_to_process = [args.sample]
        print(f"🔍 单样本模式：仅处理样本 {args.sample}")
    else:
        samples_to_process = summary_df['Sample'].unique().tolist()
        print(f"🔍 全自动模式：自动探测并处理汇总表中的全部 {len(samples_to_process)} 个样本")

    # 创建基础输出目录
    os.makedirs(args.output, exist_ok=True)
    
    total_plots = 0

    # 开始遍历每一个样本
    for sample in samples_to_process:
        sample_summary = summary_df[summary_df['Sample'] == sample]
        
        if sample_summary.empty:
            print(f"⚠️ 警告: 样本 '{sample}' 在汇总文件中没有有效的检出记录，跳过。")
            continue
            
        valid_viruses = sample_summary['Virus'].tolist()
        
        # 自动寻找对应的 SiteDepth 文件
        depth_file = os.path.join(args.depth_dir, f"{sample}.SiteDepth.gz")
        if not os.path.exists(depth_file):
            # 兼容 pandepth 可能输出的 chr.SiteDepth 命名形式
            alt_depth_file = os.path.join(args.depth_dir, f"{sample}.chr.SiteDepth.gz")
            if os.path.exists(alt_depth_file):
                depth_file = alt_depth_file
            else:
                print(f"⚠️ 警告: 未找到样本 {sample} 的深度文件 ({depth_file})，跳过该样本的绘图。")
                continue

        print(f"\n▶ [{sample}] 找到 {len(valid_viruses)} 个核心代表株，正在读取深度文件...")
        
        try:
            with gzip.open(depth_file, 'rt') as f:
                depth_df = pd.read_csv(f, header=None, sep='\t', names=['chr', 'position', 'depth'])
        except Exception as e:
            print(f"❌ 读取深度文件出错: {e}")
            continue
        
        # 仅保留需要绘制的病毒数据，防止极大数据集 OOM
        depth_df = depth_df[depth_df['chr'].isin(valid_viruses)]
        
        # 创建该样本专属的独立输出目录！
        sample_out_dir = os.path.join(args.output, str(sample))
        os.makedirs(sample_out_dir, exist_ok=True)
        
        # 绘制该样本下每一个病毒的图表
        for _, v_stat in sample_summary.iterrows():
            virus = v_stat['Virus']
            taxonomy = str(v_stat.get('taxonomy', 'Unannotated'))
            
            chrom_df = depth_df[depth_df['chr'] == virus]
            if chrom_df.empty:
                print(f"  - ⚠️ 病毒 {virus} 缺失深度坐标数据，已跳过。")
                continue
                
            x = chrom_df['position']
            y = chrom_df['depth']
            
            # 平滑处理
            y_smooth = smooth(y, args.window)
            mean_depth = v_stat['MeanDepth']
            
            # 初始化画布
            fig, ax = plt.subplots(figsize=(12, 6))
            
            # 绘制主曲线和底色填充 (科技蓝)
            ax.fill_between(x, 0, y_smooth, alpha=0.3, color='#1f77b4')
            ax.plot(x, y_smooth, color='#1f77b4', linewidth=1.5, label=f'Smoothed Depth (Window={args.window})')
            
            # 绘制平均深度红虚线
            ax.axhline(y=mean_depth, color='#d62728', linestyle='--', linewidth=2, 
                       label=f'Mean Depth: {mean_depth:.2f}x')
            
            # 添加左上角数据透视面板
            info_text = create_info_text(v_stat)
            ax.annotate(info_text, 
                         xy=(0.02, 0.96), 
                         xycoords='axes fraction',
                         bbox=dict(boxstyle="round,pad=0.6", fc="#f8f9fa", ec="#ced4da", alpha=0.9),
                         fontsize=args.fontsize, 
                         family='monospace',
                         ha='left', va='top')
            
            # 设置双行直观标题
            short_tax = textwrap.shorten(taxonomy, width=65, placeholder='...')
            ax.set_title(f'Sample: {sample}\n{short_tax} ({virus})', fontsize=14, fontweight='bold', pad=15)
            ax.set_xlabel('Genome Position (bp)', fontsize=12)
            ax.set_ylabel('Sequencing Depth (x)', fontsize=12)
            
            # 图例与网格样式
            ax.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98), framealpha=0.9)
            ax.grid(True, linestyle=':', alpha=0.6)
            
            # 固定坐标轴底端为0
            ax.set_ylim(bottom=0)
            
            plt.tight_layout()
            
            # 获取 taxonomy 并清理非法字符
            tax = str(v_stat.get('taxonomy', 'Unknown'))
            # 允许字母、数字、下划线、连字符、点，其余替换为下划线
            safe_tax = re.sub(r'[^\w\-.]', '_', tax)
            safe_vname = re.sub(r'[^\w\-.]', '_', virus)

            # 构建新文件名：Sample_taxonomy_Virus.pdf
            output_filename = f"{sample}_{safe_tax}_{safe_vname}_depth..pdf"
            output_path = os.path.join(sample_out_dir, output_filename)
            plt.savefig(output_path, bbox_inches='tight', dpi=300)
            plt.close(fig)
            total_plots += 1
            
            print(f"  - ✅ 绘制完成: {virus} -> {os.path.basename(output_path)}")

    print(f"\n🎉 批量绘图任务全部完成，共分类保存了 {total_plots} 张深度分布图。")

if __name__ == "__main__":
    main()
