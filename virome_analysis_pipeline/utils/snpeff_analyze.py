#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import gzip
import logging
from multiprocessing import Pool
from functools import partial

# ================= 配置日志 =================
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

def parse_ann_to_tsv(ann_vcf, out_tsv):
    """解析 SnpEff 的 ANN 字段为简易 TSV"""
    opener = gzip.open if ann_vcf.endswith('.gz') else open
    try:
        with opener(ann_vcf, 'rt') as f_in, open(out_tsv, 'w') as f_out:
            f_out.write("CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\tDNA_CHANGE\tAA_CHANGE\n")
            for line in f_in:
                if line.startswith("#"): continue
                cols = line.strip().split("\t")
                if len(cols) < 8: continue
                
                info = cols[7]
                ann_field = [x for x in info.split(";") if x.startswith("ANN=")]
                if not ann_field: continue
                
                # 提取第一个最高权重的注释
                annotations = ann_field[0][4:].split(",")
                first_ann = annotations[0].split("|")
                
                if len(first_ann) > 10:
                    res = [
                        cols[0], cols[1], cols[3], cols[4], 
                        first_ann[3],   # Gene Name
                        first_ann[1],   # Annotation/Effect
                        first_ann[2],   # Impact
                        first_ann[9],   # HGVS.c
                        first_ann[10]   # HGVS.p
                    ]
                    f_out.write("\t".join(res) + "\n")
    except Exception as e:
        logger.error(f"解析 TSV 失败 {ann_vcf}: {e}")

def process_single_vcf(vcf_path, jar, config, db, out_dir, mem):
    """单个 VCF 的完整处理流程：注释 -> 解析"""
    base_name = os.path.basename(vcf_path).replace(".vcf.gz", "").replace(".vcf", "")
    ann_vcf = os.path.join(out_dir, f"{base_name}.ann.vcf")
    out_tsv = os.path.join(out_dir, f"{base_name}.summary.tsv")
    
    logger.info(f"开始处理: {vcf_path}")
    
    # 步骤 1: SnpEff 注释
    cmd = [
        "java", f"-Xmx{mem}", "-jar", os.path.expanduser(jar), 
        "ann", "-c", os.path.expanduser(config), "-noStats", 
        db, vcf_path
    ]
    
    try:
        with open(ann_vcf, "w") as f:
            subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)
        
        # 步骤 2: 转换为 TSV
        parse_ann_to_tsv(ann_vcf, out_tsv)
        logger.info(f"完成处理: {base_name}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"SnpEff 运行出错 {vcf_path}: {e.stderr.decode()}")
        return False
    except Exception as e:
        logger.error(f"处理任务时发生未知错误 {vcf_path}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="高性能并行 VCF 注释分析工具")
    
    # 输入输出
    parser.add_argument("-i", "--input", required=True, help="输入路径 (单个 VCF 文件或包含多个 VCF 的文件夹)")
    parser.add_argument("-d", "--db", required=True, help="SnpEff 数据库名称")
    parser.add_argument("-o", "--out-dir", default="snpeff_results", help="结果输出目录")
    
    # 性能配置
    parser.add_argument("-j", "--jobs", type=int, default=1, help="并行任务数 (默认: 1)")
    parser.add_argument("--mem", default="2g", help="每个任务占用的 Java 内存 (默认: 2g)")
    
    # 环境配置
    parser.add_argument("--jar", default="~/snpEff/snpEff.jar", help="snpEff.jar 路径")
    parser.add_argument("--config", default="~/snpEff/snpEff.config", help="snpEff.config 路径")
    
    args = parser.parse_args()

    # 1. 确定输入文件列表
    input_path = os.path.abspath(args.input)
    vcf_files = []
    
    if os.path.isdir(input_path):
        for f in os.listdir(input_path):
            if f.endswith(".vcf") or f.endswith(".vcf.gz"):
                vcf_files.append(os.path.join(input_path, f))
        logger.info(f"文件夹模式：在 {input_path} 中找到 {len(vcf_files)} 个 VCF 文件")
    elif os.path.isfile(input_path):
        vcf_files.append(input_path)
    else:
        logger.error(f"无效的输入路径: {input_path}")
        sys.exit(1)

    if not vcf_files:
        logger.warning("未找到待处理的 VCF 文件。")
        return

    # 2. 创建输出目录
    os.makedirs(args.out_dir, exist_ok=True)

    # 3. 并行执行
    # 使用 partial 固定不变的参数
    worker_func = partial(
        process_single_vcf, 
        jar=args.jar, 
        config=args.config, 
        db=args.db, 
        out_dir=args.out_dir, 
        mem=args.mem
    )

    logger.info(f"启动并行处理任务，最大进程数: {args.jobs}")
    with Pool(processes=args.jobs) as pool:
        results = pool.map(worker_func, vcf_files)

    # 4. 统计结果
    success_count = sum(1 for r in results if r)
    logger.info(f"--- 任务报告 ---")
    logger.info(f"总计: {len(vcf_files)}, 成功: {success_count}, 失败: {len(vcf_files) - success_count}")
    logger.info(f"所有结果保存在目录: {args.out_dir}")

if __name__ == "__main__":
    main()
