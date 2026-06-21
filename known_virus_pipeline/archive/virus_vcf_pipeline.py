#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import subprocess
import argparse
import shutil
import sys

def check_dependencies():
    """检查所需的命令行工具是否已安装"""
    tools =["bcftools", "bgzip", "tabix", "VCF2PCACluster", "VCF2Dis", "awk"]
    missing_tools =[]
    for tool in tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)
    if missing_tools:
        print(f"[错误] 找不到以下工具: {', '.join(missing_tools)}")
        sys.exit(1)

def run_command(cmd, return_output=False):
    try:
        if return_output:
            return subprocess.check_output(cmd, shell=True, text=True, executable='/bin/bash')
        else:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print(f"[错误] 命令执行失败: {cmd}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="完美解决病毒VCF合并、空SNP崩溃及样本名过长的问题")
    parser.add_argument("-d", "--dir", required=True, help="包含各个样本VCF文件的原始目录路径")
    parser.add_argument("-p", "--pattern", default="**/*.filtered.vcf", help="VCF文件的匹配模式")
    parser.add_argument("-o", "--out_dir", default=".", help="输出目录（默认为当前目录）")
    parser.add_argument("--prefix", default="merged_virus", help="输出文件的前缀（默认为 merged_virus）")
    
    args = parser.parse_args()
    
    # 绝对化输入目录和输出目录
    input_dir = os.path.abspath(args.dir)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)                 # 确保输出目录存在
    
    prefix = args.prefix
    
    check_dependencies()

    search_pattern = os.path.join(input_dir, args.pattern)
    vcf_files = glob.glob(search_pattern, recursive=True)
    if not vcf_files:
        vcf_files = glob.glob(os.path.join(input_dir, "*.filtered.vcf"))
    if not vcf_files:
        print(f"[错误] 在 {input_dir} 中未找到任何 VCF 文件。")
        sys.exit(1)
        
    print(f"\n找到 {len(vcf_files)} 个 VCF 文件准备处理。")
    
    # 临时工作目录放在输出目录下，处理完后删除
    work_dir = os.path.join(out_dir, f".{prefix}_workdir")
    os.makedirs(work_dir, exist_ok=True)
    
    ready_vcfs =[]
    
    # ---------------------------------------------------------
    # 步骤 1: 智能修复 iVar VCF，并重命名样本 (去除 _ 及其后面的部分)
    # ---------------------------------------------------------
    print("\n================[步骤 1] 修复缺失列、去除 _ 后缀并重写 Header ================")
    for vcf in vcf_files:
        filename = os.path.basename(vcf)
        
        # 改名逻辑: 
        # 1. split('.')[0] 去除扩展名, 例: SRR2344926_OR489165.1.filtered.vcf -> SRR2344926_OR489165
        # 2. split('_')[0] 去除下划线后半截, 例: SRR2344926_OR489165 -> SRR2344926
        sample_name = filename.split('.')[0].split('_')[0]
        
        out_vcf_gz = os.path.join(work_dir, f"{sample_name}.vcf.gz")
        
        awk_cmd = (
            f"awk -v sample='{sample_name}' '"
            r'BEGIN { FS="\t"; OFS="\t" } '
            r'/^##/ { print $0; next } '
            r'/^#CHROM/ { '
            r'  if (NF < 9) { '
            r'      print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; '
            r'      print $0, "FORMAT", sample; '
            r'      add_cols=1 '
            r'  } else { '
            r'      $10=sample; '
            r'      print $0; '
            r'      add_cols=0 '
            r'  }; '
            r'  next '
            r'} '
            r'{ if (add_cols==1) { print $0, "GT", "1/1" } else { print $0 } }'
            f"' '{vcf}' | bgzip -c > '{out_vcf_gz}'"
        )
        run_command(awk_cmd)
        run_command(f"tabix -p vcf '{out_vcf_gz}'")
        ready_vcfs.append(out_vcf_gz)

    # ---------------------------------------------------------
    # 步骤 2: 合并 VCF 并严格过滤 (避免 0 SNPs 问题)
    # ---------------------------------------------------------
    print("\n================ [步骤 2] 合并并严格过滤 VCF (使用 -0 修复缺失率) ================")
    raw_merged_vcf = os.path.join(work_dir, "raw_merged.vcf.gz")
    merged_vcf = os.path.join(out_dir, f"{prefix}.vcf.gz")
    
    vcf_list_str = " ".join([f"'{v}'" for v in ready_vcfs])
    # -0 参数解决缺失率过高，视为参考基因型
    run_command(f"bcftools merge -0 {vcf_list_str} -O z -o '{raw_merged_vcf}'")
    
    # 过滤多等位基因和 INDEL
    run_command(f"bcftools view -m2 -M2 -v snps '{raw_merged_vcf}' -O z -o '{merged_vcf}'")
    run_command(f"tabix -p vcf '{merged_vcf}'")

    # ---------------------------------------------------------
    # 步骤 3: 规避 10 字符限制 (临时改名短 ID)
    # ---------------------------------------------------------
    print("\n================ [步骤 3] 映射样本名以规避 VCF2Dis 10字符限制 ================")
    samples_str = run_command(f"bcftools query -l '{merged_vcf}'", return_output=True)
    samples = samples_str.strip().split('\n')
    
    mapping_dict = {}
    rename_lines =[]
    for i, s_name in enumerate(samples):
        short_id = f"S{i:05d}"  # 生成如 S00000 的极短ID，绝对安全
        mapping_dict[short_id] = s_name
        rename_lines.append(f"{s_name} {short_id}")
    
    mapping_file = os.path.join(work_dir, "rename_map.txt")
    with open(mapping_file, "w") as f:
        f.write("\n".join(rename_lines) + "\n")
        
    short_vcf = os.path.join(work_dir, "short_merged.vcf.gz")
    run_command(f"bcftools reheader -s '{mapping_file}' '{merged_vcf}' > '{short_vcf}'")
    run_command(f"tabix -p vcf '{short_vcf}'")

    # ---------------------------------------------------------
    # 步骤 4 & 5: 运行 VCF2PCACluster 和 VCF2Dis
    # ---------------------------------------------------------
    print("\n================ [步骤 4] 运行 PCA 与 距离矩阵计算 ================")
    pca_out_prefix = os.path.join(out_dir, f"{prefix}_PCA")
    cmd_pca = f"VCF2PCACluster -InVCF '{short_vcf}' -OutPut '{pca_out_prefix}'"
    run_command(cmd_pca)

    dis_out_mat = os.path.join(out_dir, f"{prefix}_Distance.mat")
    cmd_dis = f"VCF2Dis -InPut '{short_vcf}' -OutPut '{dis_out_mat}'"
    run_command(cmd_dis)

    # ---------------------------------------------------------
    # 步骤 5: 将结果文件中的短 ID 还原为真实样本名
    # ---------------------------------------------------------
    print("\n================[步骤 5] 还原矩阵文件中的真实样本名 ================")
    
    if os.path.exists(dis_out_mat):
        with open(dis_out_mat, 'r') as f:
            lines = f.readlines()
        with open(dis_out_mat, 'w') as f:
            for line in lines:
                parts = line.rstrip('\n').split('\t')
                new_parts = [mapping_dict.get(p, p) for p in parts]
                f.write('\t'.join(new_parts) + '\n')
                
    pca_file = f"{pca_out_prefix}.pca"
    if os.path.exists(pca_file):
        with open(pca_file, 'r') as f:
            lines = f.readlines()
        with open(pca_file, 'w') as f:
            for line in lines:
                parts = line.strip().split()
                if not parts: continue
                parts[0] = mapping_dict.get(parts[0], parts[0])
                f.write('\t'.join(parts) + '\n')

    # 清理中间文件
    print("\n================ 清理中间文件 ================")
    shutil.rmtree(work_dir)
    print(f"\n【恭喜！全部流程顺利完成】")
    print(f"主 VCF 文件: {merged_vcf}")
    print(f"PCA 结果文件: {pca_file}")
    print(f"距离矩阵文件: {dis_out_mat}")

if __name__ == "__main__":
    main()
