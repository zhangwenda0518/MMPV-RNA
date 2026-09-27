import pandas as pd
import glob
import os
import argparse
import sys

def parse_args():
    """解析命令行参数并自动生成帮助信息"""
    
    # 描述信息会在输入 -h 时显示在开头
    parser = argparse.ArgumentParser(
        description="""
======================================================================
病毒 VCF 注释结果提取工具 (SnpEff TSV -> OncoPrint Matrix)
----------------------------------------------------------------------
该脚本用于批量读取 SnpEff 生成的 .summary.tsv 文件，
过滤掉非编码区变异，并将数据转换为适用于 R 语言 ComplexHeatmap 
(OncoPrint) 绘图的宽格式 (Wide-format) 矩阵。
======================================================================
        """,
        formatter_class=argparse.RawTextHelpFormatter # 允许在描述中使用换行符
    )

    # 添加输入目录参数 (可选，默认当前目录)
    parser.add_argument(
        "-i", "--input_dir", 
        type=str, 
        default=".", 
        help="包含 SnpEff summary.tsv 文件的目录路径。\n(默认: '.' 当前目录)"
    )

    # 添加文件后缀匹配参数 (可选，默认 .variants.summary.tsv)
    parser.add_argument(
        "-s", "--suffix", 
        type=str, 
        default=".variants.summary.tsv", 
        help="需要匹配的文件后缀名。\n(默认: '.variants.summary.tsv')"
    )

    # 添加输出文件参数 (可选，默认 oncoprint_matrix.csv)
    parser.add_argument(
        "-o", "--output", 
        type=str, 
        default="oncoprint_matrix.csv", 
        help="生成的 OncoPrint 矩阵输出路径和文件名。\n(默认: 'oncoprint_matrix.csv')"
    )

    return parser.parse_args()

def process_virus_variants(args):
    """核心处理逻辑"""
    
    # 1. 构建搜索路径并获取文件
    search_pattern = os.path.join(args.input_dir, f"*{args.suffix}")
    tsv_files = glob.glob(search_pattern)
    
    if not tsv_files:
        print(f"[错误] 在 '{args.input_dir}' 目录下未找到匹配 '{args.suffix}' 的文件！", file=sys.stderr)
        print("提示: 可以使用 -i 指定目录，或使用 -s 修改匹配的后缀名。", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] 找到 {len(tsv_files)} 个匹配的文件，开始处理...")
    
    all_data = []

    # 2. 遍历文件，提取关键信息
    for file in tsv_files:
        # 提取纯文件名（去除路径）
        basename = os.path.basename(file)
        # 从文件名提取样本名，基于第一个点号分割
        sample_name = basename.split('.')[0]
        
        try:
            df = pd.read_csv(file, sep='\t')
        except Exception as e:
            print(f"[警告] 读取文件 {basename} 失败: {e}", file=sys.stderr)
            continue
            
        # 检查是否包含必要的列
        required_columns = {'GENE', 'EFFECT', 'IMPACT'}
        if not required_columns.issubset(df.columns):
            print(f"[警告] 文件 {basename} 缺少必要的列 (GENE, EFFECT, IMPACT)，跳过该文件。", file=sys.stderr)
            continue
        
        # 过滤掉非编码区变异
        coding_df = df[df['IMPACT'].isin(['HIGH', 'MODERATE', 'LOW'])].copy()
        
        if coding_df.empty:
            continue
            
        # 简化突变类型的名称
        effect_mapping = {
            'missense_variant': 'Missense',
            'synonymous_variant': 'Synonymous',
            'nonsense_variant': 'Nonsense',
            'stop_gained': 'Nonsense',
            'frameshift_variant': 'Frameshift'
        }
        
        coding_df['Simplified_Effect'] = coding_df['EFFECT'].map(
            lambda x: effect_mapping.get(x, x.replace('_', ' ').title())
        )
        
        coding_df['Sample'] = sample_name
        subset = coding_df[['Sample', 'GENE', 'Simplified_Effect']]
        all_data.append(subset)

    if not all_data:
        print("[警告] 所有文件中均未提取到有效的编码区变异数据。", file=sys.stderr)
        sys.exit(0)

    # 3. 将所有样本的数据合并
    combined_df = pd.concat(all_data, ignore_index=True)

    # 4. 聚合数据 (同一基因存在多种突变时用分号连接)
    aggregated = combined_df.groupby(['GENE', 'Sample'])['Simplified_Effect'].apply(
        lambda x: ';'.join(x.unique())
    ).reset_index()

    # 5. 长数据转宽数据
    oncoprint_matrix = aggregated.pivot(index='GENE', columns='Sample', values='Simplified_Effect')
    oncoprint_matrix = oncoprint_matrix.fillna("")
    
    # 6. 输出结果
    try:
        # 【新增逻辑】：自动提取输出路径中的文件夹，如果不存在则自动创建
        output_dir = os.path.dirname(args.output)
        if output_dir:  # 确保路径不为空（即不仅是一个纯文件名）
            os.makedirs(output_dir, exist_ok=True)
            
        oncoprint_matrix.to_csv(args.output)
        print(f"[SUCCESS] 处理完成！矩阵已成功保存为: {args.output}")
        print("\n--- 矩阵预览 (前5行) ---")
        print(oncoprint_matrix.head())
    except Exception as e:
        print(f"[错误] 保存文件失败: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # 解析命令行参数
    args = parse_args()
    # 运行主程序
    process_virus_variants(args)
