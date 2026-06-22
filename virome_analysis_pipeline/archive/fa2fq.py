#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import sys
import os

def parse_fasta(file_handle):
    """生成器：读取FASTA文件，自动处理多行序列拼接"""
    header = None
    seq_lines = []
    for line in file_handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(seq_lines)
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        yield header, ''.join(seq_lines)

def generate_out_filename(in_fasta):
    """根据输入文件名自动推导输出文件名"""
    if in_fasta.endswith('.fasta.gz'):
        return in_fasta.replace('.fasta.gz', '.fastq.gz')
    elif in_fasta.endswith('.fa.gz'):
        return in_fasta.replace('.fa.gz', '.fastq.gz')
    elif in_fasta.endswith('.fasta'):
        return in_fasta.replace('.fasta', '.fastq')
    elif in_fasta.endswith('.fa'):
        return in_fasta.replace('.fa', '.fastq')
    else:
        return f"{in_fasta}.fastq"

# ！！！核心修改在这里：将默认字符改为 'h' (ASCII 104) ！！！
# 这样在 Phred+64 下是 Q40，在 Phred+33 下是 Q71，全平台高质量
def convert_fasta_to_fastq(in_fasta, dummy_qual='h'):
    """核心转换逻辑"""
    out_fastq = generate_out_filename(in_fasta)
    print(f"[*] 正在转换: {in_fasta} -> {out_fastq}")
    
    # 自动判断是否为 gzip 文件并选择合适的打开方式
    open_in = gzip.open if in_fasta.endswith('.gz') else open
    open_out = gzip.open if out_fastq.endswith('.gz') else open
    mode_in = 'rt' if in_fasta.endswith('.gz') else 'r'
    mode_out = 'wt' if out_fastq.endswith('.gz') else 'w'
    
    count = 0
    try:
        with open_in(in_fasta, mode_in) as f_in, open_out(out_fastq, mode_out) as f_out:
            for header, seq in parse_fasta(f_in):
                seq_len = len(seq)
                qual_str = dummy_qual * seq_len
                # 写入 FASTQ 标准的 4 行格式
                f_out.write(f"@{header}\n{seq}\n+\n{qual_str}\n")
                count += 1
    except Exception as e:
        print(f"[!] 处理 {in_fasta} 时出错: {e}", file=sys.stderr)
        sys.exit(1)
        
    print(f"[+] 完成！共转换 {count} 条序列。\n")

def main():
    parser = argparse.ArgumentParser(
        description="高效 FASTA 转 FASTQ 脚本 (伪造质量值 'h': 兼容 Phred33的Q71 和 Phred64的Q40)",
        usage="%(prog)s --paired <R1.fasta.gz> <R2.fasta.gz>"
    )
    
    # 设置 --paired 参数接收两个文件
    parser.add_argument('--paired', nargs=2, metavar=('R1_FASTA', 'R2_FASTA'),
                        help='输入双端配稳的 FASTA 文件 (支持 .gz 压缩格式)')
    
    # 也可以作为一个普通参数处理单端文件，增加脚本泛用性
    parser.add_argument('single_inputs', nargs='*', 
                        help='输入单端 FASTA 文件')

    args = parser.parse_args()

    if not args.paired and not args.single_inputs:
        parser.print_help()
        sys.exit(1)

    # 处理双端数据
    if args.paired:
        for fasta_file in args.paired:
            if not os.path.exists(fasta_file):
                print(f"[!] 文件不存在: {fasta_file}", file=sys.stderr)
                sys.exit(1)
            convert_fasta_to_fastq(fasta_file)

    # 处理额外传入的单端数据
    if args.single_inputs:
        for fasta_file in args.single_inputs:
            if os.path.exists(fasta_file):
                 convert_fasta_to_fastq(fasta_file)

if __name__ == '__main__':
    main()
