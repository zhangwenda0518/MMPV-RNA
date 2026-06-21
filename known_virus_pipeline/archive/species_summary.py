#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Summarize and visualize unique virus species counts per query')
    parser.add_argument('-i', '--input', required=True,
                        help='Path to annotated_results.tsv file')
    parser.add_argument('-s', '--summary', required=True,
                        help='Output summary file path')
    parser.add_argument('-p', '--plot', 
                        help='Output plot file path (optional)')
    parser.add_argument('-n', '--top_n', type=int, default=10,
                        help='Number of top species to display (default: 10)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Plot resolution (default: 300)')
    parser.add_argument('--width', type=float, default=10.0,
                        help='Plot width in inches (default: 10)')
    parser.add_argument('--height', type=float, default=6.0,
                        help='Plot height in inches (default: 6)')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # 读取注释结果文件
    print(f"Reading input file: {args.input}")
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        sys.exit(f"Error reading input file: {str(e)}")
    
    # 检查必要的列是否存在
    required_columns = ['Query', 'Species']
    for col in required_columns:
        if col not in df.columns:
            sys.exit(f"Error: Required column '{col}' not found in input file")
    
    # 按Query和Species去重 - 每个Query只计一次每个物种
    print("Removing duplicate species per query...")
    unique_df = df.drop_duplicates(subset=['Query', 'Species'])
    
    # 统计每个物种出现的次数（每个物种被多少个不同的Query注释到）
    print("Counting unique species occurrences...")
    species_counts = unique_df['Species'].value_counts().reset_index()
    species_counts.columns = ['Species', 'Unique_Query_Count']
    
    # 保存统计结果
    print(f"Writing summary to: {args.summary}")
    species_counts.to_csv(args.summary, sep='\t', index=False)
    
    # 打印总体统计信息
    total_queries = df['Query'].nunique()
    total_species = len(species_counts)
    print("\nSummary statistics:")
    print(f"  Total unique queries: {total_queries:,}")
    print(f"  Total unique species: {total_species:,}")
    
    # 如果需要绘图
    if args.plot:
        print(f"Generating plot: {args.plot}")
        
        # 准备数据
        plot_df = species_counts.head(args.top_n).copy()
        
        # 创建图表
        plt.figure(figsize=(args.width, args.height))
        
        # 使用seaborn创建水平条形图
        ax = sns.barplot(
            x='Unique_Query_Count', 
            y='Species', 
            data=plot_df,
            palette='viridis',
            edgecolor='black'
        )
        
        # 添加计数标签
        for i, (count, species) in enumerate(zip(plot_df['Unique_Query_Count'], plot_df['Species'])):
            ax.text(count + max(plot_df['Unique_Query_Count'])*0.01, i, 
                    f'{count:,}', va='center', fontsize=9)
        
        # 设置标题和标签
        plt.title(f'Top {args.top_n} Virus Species (Unique Query Count)', fontsize=14)
        plt.xlabel('Number of Unique Queries', fontsize=12)
        plt.ylabel('Virus Species', fontsize=12)
        
        # 优化布局
        plt.tight_layout()
        
        # 保存图表
        plt.savefig(args.plot, dpi=args.dpi)
        print(f"Plot saved to: {args.plot}")
        
        # 显示前10个物种的信息
        print("\nTop species:")
        for i, row in plot_df.iterrows():
            print(f"  {i+1}. {row['Species']}: {row['Unique_Query_Count']:,} unique queries")

if __name__ == "__main__":
    main()
