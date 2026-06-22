#!/usr/bin/env python3
import re
import argparse

def convert_megahit_to_rnaspades(input_file, output_file):
    """
    将MEGAHIT头信息转换为rnaSPAdes格式
    输入格式：>k141_45489 flag=1 multi=1.0000 len=327
    输出格式：>NODE_36652_length_943_cov_17349.429885_g18276_i0
    """
    pattern = re.compile(r'>(\S+)\s+flag=(\d+)\s+multi=([\d.]+)\s+len=(\d+)')
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                match = pattern.match(line.strip())
                if match:
                    contig_id, flag, multi, length = match.groups()
                    # 提取数字ID部分 (k141_45489 -> 45489)
                    node_id = contig_id.split('_')[-1] if '_' in contig_id else contig_id
                    # 创建rnaSPAdes格式头
                    new_header = f">NODE_{node_id}_length_{length}_cov_{multi}_g{node_id}_i0"
                    fout.write(new_header + '\n')
                else:
                    fout.write(line)
            else:
                fout.write(line)

def convert_rnaspades_to_megahit(input_file, output_file):
    """
    将rnaSPAdes头信息转换为MEGAHIT格式
    输入格式：>NODE_36652_length_943_cov_17349.429885_g18276_i0
    输出格式：>k141_45489 flag=1 multi=1.0000 len=327
    """
    pattern = re.compile(r'>NODE_(\d+)_length_(\d+)_cov_([\d.]+)_g(\d+)_i(\d+)')
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                match = pattern.match(line.strip())
                if match:
                    node_id, length, cov, gene_id, isoform = match.groups()
                    # 创建MEGAHIT格式头 (使用基因ID作为contig ID)
                    new_header = f">g{gene_id}_i{isoform} flag=1 multi={cov} len={length}"
                    fout.write(new_header + '\n')
                else:
                    fout.write(line)
            else:
                fout.write(line)

def main():
    parser = argparse.ArgumentParser(description='转换MEGAHIT和rnaSPAdes组装结果的头信息格式')
    parser.add_argument('-i', '--input', required=True, help='输入FASTA文件')
    parser.add_argument('-o', '--output', required=True, help='输出FASTA文件')
    parser.add_argument('-d', '--direction', required=True, choices=['m2r', 'r2m'],
                        help='转换方向: m2r = MEGAHIT转rnaSPAdes, r2m = rnaSPAdes转MEGAHIT')
    
    args = parser.parse_args()
    
    if args.direction == 'm2r':
        convert_megahit_to_rnaspades(args.input, args.output)
        print(f"转换完成: MEGAHIT -> rnaSPAdes | 输出文件: {args.output}")
    elif args.direction == 'r2m':
        convert_rnaspades_to_megahit(args.input, args.output)
        print(f"转换完成: rnaSPAdes -> MEGAHIT | 输出文件: {args.output}")

if __name__ == "__main__":
    main()#!/usr/bin/env python3
import re
import argparse

def convert_megahit_to_rnaspades(input_file, output_file):
    """
    将MEGAHIT头信息转换为rnaSPAdes格式
    输入格式：>k141_45489 flag=1 multi=1.0000 len=327
    输出格式：>NODE_36652_length_943_cov_17349.429885_g18276_i0
    """
    pattern = re.compile(r'>(\S+)\s+flag=(\d+)\s+multi=([\d.]+)\s+len=(\d+)')
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                match = pattern.match(line.strip())
                if match:
                    contig_id, flag, multi, length = match.groups()
                    # 提取数字ID部分 (k141_45489 -> 45489)
                    node_id = contig_id.split('_')[-1] if '_' in contig_id else contig_id
                    # 创建rnaSPAdes格式头
                    new_header = f">NODE_{node_id}_length_{length}_cov_{multi}_g{node_id}_i0"
                    fout.write(new_header + '\n')
                else:
                    fout.write(line)
            else:
                fout.write(line)

def convert_rnaspades_to_megahit(input_file, output_file):
    """
    将rnaSPAdes头信息转换为MEGAHIT格式
    输入格式：>NODE_36652_length_943_cov_17349.429885_g18276_i0
    输出格式：>k141_45489 flag=1 multi=1.0000 len=327
    """
    pattern = re.compile(r'>NODE_(\d+)_length_(\d+)_cov_([\d.]+)_g(\d+)_i(\d+)')
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                match = pattern.match(line.strip())
                if match:
                    node_id, length, cov, gene_id, isoform = match.groups()
                    # 创建MEGAHIT格式头 (使用基因ID作为contig ID)
                    new_header = f">g{gene_id}_i{isoform} flag=1 multi={cov} len={length}"
                    fout.write(new_header + '\n')
                else:
                    fout.write(line)
            else:
                fout.write(line)

def main():
    parser = argparse.ArgumentParser(description='转换MEGAHIT和rnaSPAdes组装结果的头信息格式')
    parser.add_argument('-i', '--input', required=True, help='输入FASTA文件')
    parser.add_argument('-o', '--output', required=True, help='输出FASTA文件')
    parser.add_argument('-d', '--direction', required=True, choices=['m2r', 'r2m'],
                        help='转换方向: m2r = MEGAHIT转rnaSPAdes, r2m = rnaSPAdes转MEGAHIT')
    
    args = parser.parse_args()
    
    if args.direction == 'm2r':
        convert_megahit_to_rnaspades(args.input, args.output)
        print(f"转换完成: MEGAHIT -> rnaSPAdes | 输出文件: {args.output}")
    elif args.direction == 'r2m':
        convert_rnaspades_to_megahit(args.input, args.output)
        print(f"转换完成: rnaSPAdes -> MEGAHIT | 输出文件: {args.output}")

if __name__ == "__main__":
    main()
