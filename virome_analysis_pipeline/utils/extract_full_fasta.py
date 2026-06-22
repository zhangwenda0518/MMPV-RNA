#!/usr/bin/env python3
"""
extract_full_fasta.py
功能：
1. 从组装结果目录中读取 11.Ultimate_Circular_Result.fasta
2. 自动解析并【仅提取最长的一条序列】（过滤掉未去冗余导致的短序列）
3. 彻底重命名 Header 为 >Sample_Accession (舍弃原始 ID)
4. 输出到全新的分类结果目录中，且 FASTA 序列【不换行】(单行输出)。
"""

import argparse
from pathlib import Path

def get_longest_sequence(fasta_path):
    """
    轻量级 FASTA 解析函数，返回文件中最长的一条序列（纯字符串）
    """
    longest_seq = ""
    current_seq = []
    
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_seq:
                        seq_str = "".join(current_seq)
                        if len(seq_str) > len(longest_seq):
                            longest_seq = seq_str
                    current_seq = []
                else:
                    current_seq.append(line)
        # 处理文件末尾的最后一条序列
        if current_seq:
            seq_str = "".join(current_seq)
            if len(seq_str) > len(longest_seq):
                longest_seq = seq_str
    except Exception as e:
        print(f"读取文件错误 {fasta_path}: {e}")
        
    return longest_seq

def main():
    parser = argparse.ArgumentParser(description="提取最长 FASTA 序列并重命名(单行序列)至独立目录")
    parser.add_argument("-d", "--dir", required=True, help="输入的组装总目录 (例如 7.Virus_assemblies_final)")
    parser.add_argument("-o", "--outdir", required=True, help="输出的全新总目录 (例如 7.Virus_full.result)")
    parser.add_argument("--target_file", default="11.Ultimate_Circular_Result.fasta", help="目标文件名")
    args = parser.parse_args()

    in_base = Path(args.dir)
    out_base = Path(args.outdir)

    if not in_base.exists():
        print(f"❌ 找不到输入目录: {in_base}")
        return

    # 创建独立的输出总目录
    out_base.mkdir(parents=True, exist_ok=True)

    print(f"🔍 正在扫描 {in_base}，提取最长序列(单行格式)至 {out_base}...")
    
    success_count = 0
    missing_count = 0
    empty_count = 0

    # 遍历 Level 1: Taxonomy_Accession 目录
    for tax_dir in in_base.iterdir():
        if not tax_dir.is_dir(): continue

        # 在输出目录中创建对应的 Taxonomy 文件夹
        out_tax_dir = out_base / tax_dir.name
        out_tax_dir.mkdir(parents=True, exist_ok=True)

        # 遍历 Level 2: Sample_Accession 目录
        for sample_dir in tax_dir.iterdir():
            if not sample_dir.is_dir(): continue

            # 向下寻找目标 fasta 文件
            target_files = list(sample_dir.rglob(args.target_file))
            
            if not target_files:
                missing_count += 1
                continue

            for target_file in target_files:
                # 获取样本文件夹名字作为标准 Header，例如: CRR1126132_OR489165.1
                sample_accession = sample_dir.name
                
                # 输出文件名保持 . 分隔的习惯：CRR1126132.OR489165.1.full.fasta
                file_prefix = sample_accession.replace('_', '.')
                out_name = f"{file_prefix}.full.fasta"
                out_path = out_tax_dir / out_name
                
                # 提取最长序列
                longest_seq = get_longest_sequence(target_file)
                
                if not longest_seq:
                    empty_count += 1
                    print(f"⚠️ 警告: {target_file} 中没有找到有效序列。")
                    continue
                
                # 将最长序列写入新文件，彻底重写 Header，且序列单行输出不换行
                try:
                    with open(out_path, 'w') as fout:
                        fout.write(f">{sample_accession}\n")
                        fout.write(f"{longest_seq}\n")  # 一整条序列直接写入
                            
                    success_count += 1
                    print(f"✅ 提取最长 ({len(longest_seq)} bp): {tax_dir.name}/{out_name}")
                except Exception as e:
                    print(f"⚠️ 写入 {out_path} 时出错: {e}")

    print("-" * 60)
    print(f"🎉 提取整理完成！成功处理并生成 {success_count} 个单行(Single-line)格式的最长 FASTA 文件。")
    print(f"📁 完美的纯净版结果目录位于: {out_base.absolute()}")
    
    if missing_count > 0 or empty_count > 0:
        print(f"ℹ️  缺失文件数: {missing_count}，空/无效文件数: {empty_count}")

if __name__ == "__main__":
    main()
