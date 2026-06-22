#!/usr/bin/env python3
import argparse
import sys
import os
from Bio import SeqIO

def extract_sequences(input_gbk, full_out, nucl_out, prot_out, split):
    """
    解析 GenBank 文件并提取请求的序列。
    如果启用了 split，则拆分输出。
    """
    # 仅在非拆分模式下打开全局文件句柄
    f_handle = open(full_out, "w") if (full_out and not split) else None
    n_handle = open(nucl_out, "w") if (nucl_out and not split) else None
    p_handle = open(prot_out, "w") if (prot_out and not split) else None

    # 确保输出目录存在
    for out_path in [full_out, nucl_out, prot_out]:
        if out_path:
            os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    def get_split_filename(base_name, identifier):
        """根据基础文件名和标识符生成拆分后的文件名"""
        if not base_name: return None
        name, ext = os.path.splitext(base_name)
        # 移除非法字符，防止路径错误
        safe_id = str(identifier).replace("/", "_").replace("\\", "_")
        return f"{name}_{safe_id}{ext}"

    try:
        for seq_record in SeqIO.parse(input_gbk, "genbank"):
            print(f"[*] 正在处理 GenBank 记录: {seq_record.id}")

            # 1. 提取全长序列
            if full_out:
                full_fasta = f">{seq_record.id} {seq_record.description}\n{seq_record.seq}\n"
                if split:
                    # 按 record.id 拆分全长序列
                    split_f_name = get_split_filename(full_out, seq_record.id)
                    with open(split_f_name, "w") as sf:
                        sf.write(full_fasta)
                else:
                    f_handle.write(full_fasta)

            # 2. 提取 CDS 核酸和蛋白序列
            if nucl_out or prot_out:
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS":
                        # 获取基因标签，如果没有 locus_tag 则退而求其次
                        locus_tag = seq_feature.qualifiers.get('locus_tag', 
                                    seq_feature.qualifiers.get('gene', ['unknown_locus']))[0]
                        protein_id = seq_feature.qualifiers.get('protein_id', [locus_tag])[0]

                        # -- 提取 CDS 核酸序列 --
                        if nucl_out:
                            cds_seq = seq_feature.extract(seq_record.seq)
                            nucl_fasta = f">{locus_tag} [locus_tag] | {seq_record.id} | CDS_Nucleotide\n{cds_seq}\n"
                            if split:
                                # 按 locus_tag 拆分输出核酸序列
                                split_n_name = get_split_filename(nucl_out, locus_tag)
                                with open(split_n_name, "w") as sn:
                                    sn.write(nucl_fasta)
                            else:
                                n_handle.write(nucl_fasta)

                        # -- 提取 CDS 蛋白翻译序列 --
                        if prot_out:
                            translation = seq_feature.qualifiers.get('translation', [None])[0]
                            if translation:
                                prot_fasta = f">{protein_id} [protein_id] | locus:{locus_tag} | {seq_record.id}\n{translation}\n"
                                if split:
                                    # 按 locus_tag 拆分输出蛋白序列
                                    split_p_name = get_split_filename(prot_out, locus_tag)
                                    with open(split_p_name, "w") as sp:
                                        sp.write(prot_fasta)
                                else:
                                    p_handle.write(prot_fasta)
                            else:
                                print(f"    - 警告: {locus_tag} 没有找到 translation 标签，已跳过。")

        print("\n[+] 提取完成！")

    except Exception as e:
        print(f"[-] 处理过程中发生错误: {e}")
        sys.exit(1)

    finally:
        # 关闭全局文件句柄（如果存在）
        if f_handle: f_handle.close()
        if n_handle: n_handle.close()
        if p_handle: p_handle.close()

def main():
    parser = argparse.ArgumentParser(
        description="从 GenBank (.gb/.gbk) 文件中提取全长序列、CDS 核酸序列及蛋白序列。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("-i", "--input", required=True, 
                        help="输入的 GenBank 文件路径 (必需)")
    parser.add_argument("-f", "--full", 
                        help="输出文件名：提取全长序列 (Fasta 格式)")
    parser.add_argument("-n", "--nucl", 
                        help="输出文件名：提取 CDS 的核酸序列 (Fasta 格式)")
    parser.add_argument("-p", "--prot", 
                        help="输出文件名：提取 CDS 的蛋白翻译序列 (Fasta 格式)")
    
    # 新增 --split 参数
    parser.add_argument("--split", action="store_true", 
                        help="启用拆分模式：将每个基因的 CDS 和蛋白独立输出为单独的文件 (以 locus_tag 命名)")

    args = parser.parse_args()

    if not (args.full or args.nucl or args.prot):
        print("[-] 错误：请至少指定一个输出参数 (-f, -n 或 -p)。\n")
        parser.print_help()
        sys.exit(1)

    extract_sequences(args.input, args.full, args.nucl, args.prot, args.split)

if __name__ == "__main__":
    main()
