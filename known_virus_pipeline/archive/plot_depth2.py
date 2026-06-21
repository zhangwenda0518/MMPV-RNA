import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import argparse
import os
import gzip
import textwrap

def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def format_number(value):
    """格式化大数字，添加千位分隔符"""
    return f"{value:,.0f}" if isinstance(value, (int, float)) else value

def create_info_text(chrom_stat):
    """创建格式化的统计信息文本"""
    info = [
        f"Chromosome: {chrom_stat['Chr']}",
        f"Length: {format_number(chrom_stat['Length'])} bp",
        f"Covered Sites: {format_number(chrom_stat['CoveredSite'])}",
        f"Total Depth: {format_number(chrom_stat['TotalDepth'])}",
        f"Coverage: {chrom_stat['Coverage(%)']:.2f}%",
        f"Mean Depth: {chrom_stat['MeanDepth']:.2f}x"
    ]
    return "\n".join(info)

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description='Plot smoothed depth profile with mean depth line')
    parser.add_argument('-s', '--sitedepth', required=True, 
                        help='Single site depth file (gzipped)')
    parser.add_argument('-d', '--depthstat', required=True, 
                        help='Depth statistics file (gzipped)')
    parser.add_argument('-o', '--output', required=True, 
                        help='Output directory for plots')
    parser.add_argument('-w', '--window', type=int, default=100,
                        help='Smoothing window size (default: 100)')
    parser.add_argument('-c', '--coverage_threshold', type=float, default=90.0,
                        help='Coverage percentage threshold (default: 90.0)')
    parser.add_argument('-f', '--fontsize', type=float, default=10,
                        help='Font size for annotations (default: 10)')
    args = parser.parse_args()

    # 读取单碱基深度文件
    with gzip.open(args.sitedepth, 'rt') as f:
        sitedepth_df = pd.read_csv(f, header=None, sep='\t')
    sitedepth_df.columns = ['chr', 'position', 'depth']
    
    # 读取序列覆盖度统计文件
    with gzip.open(args.depthstat, 'rt') as f:
        # 跳过注释行，使用第一行作为列名
        header_line = f.readline().strip()
        # 移除可能存在的注释符号
        if header_line.startswith('#'):
            header = header_line.replace('#', '').split()
        else:
            header = header_line.split()
            f.seek(0)  # 如果第一行不是注释，则重置文件指针
            
        stat_df = pd.read_csv(f, sep='\t', header=None, comment='#', names=header)
    
    # 筛选满足覆盖度阈值的染色体
    filtered_chroms = stat_df[stat_df['Coverage(%)'] >= args.coverage_threshold]['Chr']
    chroms_to_plot = set(sitedepth_df['chr']).intersection(set(filtered_chroms))
    
    if not chroms_to_plot:
        print(f"Warning: No chromosomes found with coverage >= {args.coverage_threshold}%")
        return
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    # 为每个满足条件的染色体创建图表
    for chrom in chroms_to_plot:
        # 提取当前染色体的深度数据
        chrom_df = sitedepth_df[sitedepth_df['chr'] == chrom]
        x = chrom_df['position']
        y = chrom_df['depth']
        
        # 应用平滑处理
        if len(y) >= args.window:
            y_smooth = smooth(y, args.window)
        else:
            print(f"Warning: Chromosome {chrom} has fewer positions ({len(y)}) than window size ({args.window}). Using raw data.")
            y_smooth = y
        
        # 获取当前染色体的统计信息
        chrom_stat = stat_df[stat_df['Chr'] == chrom].iloc[0]
        mean_depth = chrom_stat['MeanDepth']
        coverage = chrom_stat['Coverage(%)']
        
        # 创建图形
        plt.figure(figsize=(12, 6))
        plt.plot(x, y_smooth, label=f'Smoothed Depth (window={args.window})')
        plt.axhline(y=mean_depth, color='r', linestyle='--', 
                    label=f'Mean Depth: {mean_depth:.2f}x')
        
        # 创建并添加统计信息文本（左上角）
        info_text = create_info_text(chrom_stat)
        plt.annotate(info_text, 
                     xy=(0.02, 0.98),  # 左上角位置
                     xycoords='axes fraction',
                     bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", alpha=0.7),
                     fontsize=args.fontsize, 
                     family='monospace',
                     ha='left',  # 水平左对齐
                     va='top')   # 垂直顶对齐
        
        plt.title(f'Depth Profile - {chrom}')
        plt.xlabel('Genome Position')
        plt.ylabel('Depth (x)')
        
        # 添加图例（右上角），自动调整位置避免重叠
        plt.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98), framealpha=0.7)
        
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # 自动调整布局确保所有元素都在坐标系内
        plt.tight_layout(pad=2.0)
        
        # 保存图片
        output_path = os.path.join(args.output, f'{chrom}_depth.pdf')
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
        print(f"Plot saved to: {output_path}")

if __name__ == "__main__":
    main()
