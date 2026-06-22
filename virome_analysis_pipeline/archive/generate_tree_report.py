import pandas as pd
import argparse
import sys

def generate_report(rna_file, virus_file, output_file):
    print(f"正在读取文件...\n 来源文件: {rna_file}\n 病毒文件: {virus_file}")
    
    # 1. 读取并处理 ncbi.RNA.tab (来源文件)
    try:
        rna_df = pd.read_csv(rna_file, sep='\t')
    except Exception as e:
        print(f"❌ 读取来源文件失败，请检查路径。错误信息: {e}")
        sys.exit(1)

    mapping_df = rna_df.melt(var_name='Source', value_name='Sample')
    mapping_df = mapping_df.dropna(subset=['Sample'])
    mapping_df['Sample'] = mapping_df['Sample'].astype(str).str.strip()
    mapping_df = mapping_df[mapping_df['Sample'] != '']

    # 统计每个 Source 的总样本数
    source_total_samples = mapping_df.groupby('Source')['Sample'].nunique().to_dict()

    # 2. 读取病毒汇总文件并自动判断模式
    try:
        virus_df = pd.read_csv(virus_file, sep='\t')
    except Exception as e:
        print(f"❌ 读取病毒汇总文件失败，请检查路径。错误信息: {e}")
        sys.exit(1)

    virus_df['Sample'] = virus_df['Sample'].astype(str).str.strip()

    # 自动识别数据类型
    if 'Adjusted_Species' in virus_df.columns:
        mode = 'map_full'
        print("✅ 自动识别数据格式: 完整版 Map 定量模式 (包含ANI、Pi等多态性指标)")
        species_col = 'Adjusted_Species'
    elif 'Coverage(%)' in virus_df.columns and 'MeanDepth' in virus_df.columns and 'Virus' in virus_df.columns:
        mode = 'map_summary'
        print("✅ 自动识别数据格式: 简版 Map 定量模式 (包含Coverage, MeanDepth, TPM等)")
        species_col = 'Virus'
    elif 'NumReads' in virus_df.columns and 'Virus' in virus_df.columns:
        mode = 'salmon'
        print("✅ 自动识别数据格式: Salmon 定量模式 (包含 NumReads, TPM)")
        species_col = 'Virus'
    else:
        print("❌ 无法识别病毒文件的表头格式，请确保输入了正确的 Map 或 Salmon 汇总文件。")
        sys.exit(1)

    # 3. 将病毒数据与样品来源进行合并 (Inner Join)
    merged_df = pd.merge(virus_df, mapping_df, on='Sample', how='inner')
    
    if merged_df.empty:
        print("⚠️ 警告：两个文件之间没有匹配到任何相同的 Sample 编号，请检查数据。")
        sys.exit(0)

    # 4. 根据不同的模式定义聚合方法
    if mode == 'map_full':
        agg_funcs = {
            'Sample': 'nunique',                 
            'Asm_CPM': 'mean',                   
            'Asm_RPM': 'mean',                   
            'Asm_FPKM': 'mean',                  
            'Asm_TPM': 'mean',                   
            'Avg_Read_ANI': 'mean',              
            'Avg_Pi': 'mean',                    
            'Rep_Length': 'first',               
            'Rep_Coverage(%)': 'mean',           
            'Rep_MeanDepth': 'mean',             
            'Rep_Accession': lambda x: x.mode()[0] if not x.empty else 'Unknown'
        }
    elif mode == 'map_summary':
        # 针对新提供的 all_viruses.best.summary.tsv 格式
        agg_funcs = {
            'Sample': 'nunique',
            'RPM': 'mean',
            'FPKM': 'mean',
            'TPM': 'mean',
            'Length': 'first',
            'Coverage(%)': 'mean',
            'MeanDepth': 'mean'
        }
    else: # mode == 'salmon'
        agg_funcs = {
            'Sample': 'nunique',
            'NumReads': 'mean',
            'TPM': 'mean'
        }

    # 按照样品来源和病毒种类进行分组计算
    summary = merged_df.groupby(['Source', species_col]).agg(agg_funcs).reset_index()

    # 5. 生成树状格式报告
    output_lines = []
    
    for source, group in summary.groupby('Source'):
        total_samples = source_total_samples.get(source, 1)
        
        header = f"========== 样品来源: {source} (总检测样本数: {total_samples}) ==========\n"
        print(header, end="")
        output_lines.append(header)
        
        # 按检出例数从高到低排序
        group = group.sort_values(by='Sample', ascending=False)
        
        for _, row in group.iterrows():
            species = row[species_col]
            count = row['Sample']
            prevalence = (count / total_samples) * 100
            
            # 构建树状字符串
            if mode == 'map_full':
                block = (
                    f"🎯 {species}: 检出 {count} 例 (群体检出率 {prevalence:.1f}%)\n"
                    f"   ├─ 平均 CPM: {row['Asm_CPM']:.2f} | 平均 RPM: {row['Asm_RPM']:.2f}\n"
                    f"   ├─ 平均 FPKM: {row['Asm_FPKM']:.2f} | 平均 TPM: {row['Asm_TPM']:.2f}\n"
                    f"   ├─ 测定 ANI: {row['Avg_Read_ANI']:.2f}% | 多态性 Pi: {row['Avg_Pi']:.5f}\n"
                    f"   ├─ 代表序列长度: {row['Rep_Length']} bp\n"
                    f"   ├─ 平均覆盖度: {row['Rep_Coverage(%)']:.2f}% | 平均深度: {row['Rep_MeanDepth']:.2f}x\n"
                    f"   └─ 代表株(Rep): {row['Rep_Accession']}\n\n"
                )
            elif mode == 'map_summary':
                block = (
                    f"🎯 {species}: 检出 {count} 例 (群体检出率 {prevalence:.1f}%)\n"
                    f"   ├─ 平均 RPM: {row['RPM']:.2f} | 平均 FPKM: {row['FPKM']:.2f} | 平均 TPM: {row['TPM']:.2f}\n"
                    f"   └─ 序列长度: {row['Length']} bp | 平均覆盖度: {row['Coverage(%)']:.2f}% | 平均深度: {row['MeanDepth']:.2f}x\n\n"
                )
            else: # mode == 'salmon'
                block = (
                    f"🎯 {species}: 检出 {count} 例 (群体检出率 {prevalence:.1f}%)\n"
                    f"   ├─ 平均 NumReads: {row['NumReads']:.2f}\n"
                    f"   └─ 平均 TPM: {row['TPM']:.2f}\n\n"
                )
            
            print(block, end="")
            output_lines.append(block)

    # 6. 保存到文本文件
    try:
        with open(output_file, "w", encoding="utf-8") as f:
            f.writelines(output_lines)
        print(f"\n✅ 报告已成功生成并保存至: {output_file}")
    except Exception as e:
        print(f"❌ 写入输出文件失败。错误信息: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据样品来源统计病毒感染情况，支持多种 Map 和 Salmon 输出结果自适应。")
    parser.add_argument("-m", "--meta", required=True, help="输入的样品来源文件，例如: ncbi.RNA.tab")
    parser.add_argument("-v", "--virus", required=True, help="输入的病毒鉴定汇总文件，例如: all_viruses.best.summary.tsv")
    parser.add_argument("-o", "--output", default="viral_tree_report.txt", help="输出的报告文件名 (默认: viral_tree_report.txt)")
    
    args = parser.parse_args()
    generate_report(args.meta, args.virus, args.output)
