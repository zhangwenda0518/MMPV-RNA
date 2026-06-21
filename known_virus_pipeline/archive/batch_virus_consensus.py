#!/usr/bin/env python3
"""
batch_viral_consensus.py

Generate consensus sequences for virus-positive samples based on
batch_virus_depth.py output.

Usage:
    python batch_viral_consensus.py --result virus_depth_analysis \\
                                    --reference plantvirus.fasta \\
                                    --jobs 10
"""

import os
import sys
import argparse
import math
import subprocess
import logging
import glob
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# 尝试导入 tqdm，若失败则给出友好提示并退出
try:
    from tqdm import tqdm
except ImportError:
    print("Error: tqdm is required. Install with 'pip install tqdm'", file=sys.stderr)
    sys.exit(1)

# 配置日志：只输出到文件，控制台不显示日志
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# 检查外部命令是否可用
REQUIRED_COMMANDS = ['seqkit', 'samtools', 'viral_consensus']
for cmd in REQUIRED_COMMANDS:
    if not subprocess.run(f'which {cmd}', shell=True, capture_output=True).returncode == 0:
        logger.error(f"Required command '{cmd}' not found in PATH")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate viral consensus sequences from batch_virus_depth.py results')
    parser.add_argument('--result', required=True, help='Result directory (e.g., virus_depth_analysis)')
    parser.add_argument('--reference', required=True, help='Reference fasta file containing all viral sequences')
    parser.add_argument('--jobs', type=int, default=4, help='Number of parallel sample processes (default: 4)')
    parser.add_argument('--min-depth', type=int, default=1, help='Minimum depth threshold if MeanDepth < 1 (default: 1)')
    return parser.parse_args()


def setup_file_logging(result_dir):
    """将日志写入 result_dir/consensus.log"""
    log_file = os.path.join(result_dir, 'consensus.log')
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # 不再添加控制台 handler，控制台只显示 tqdm 和最终 print


def read_summary_tasks(summary_file):
    """
    读取 summary 文件（all_viruses.summary.tsv 或单个病毒文件），
    返回任务列表 [(sample, virus, mean_depth_floor), ...]
    """
    tasks = []
    if not os.path.isfile(summary_file):
        logger.warning(f"Summary file not found: {summary_file}")
        return tasks

    with open(summary_file, 'r') as f:
        header = f.readline().strip().split('\t')
        try:
            sample_idx = header.index('Sample')
            virus_idx = header.index('Virus')
            meandepth_idx = header.index('MeanDepth')
        except ValueError as e:
            logger.error(f"Required column not found in {summary_file}: {e}")
            return tasks

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(sample_idx, virus_idx, meandepth_idx):
                continue
            sample = parts[sample_idx]
            virus = parts[virus_idx]
            try:
                meandepth = float(parts[meandepth_idx])
            except ValueError:
                logger.warning(f"Invalid MeanDepth value for {sample} {virus}, skipping")
                continue
            depth_floor = int(math.floor(meandepth))
            if depth_floor < 1:
                depth_floor = args.min_depth
                logger.debug(f"{sample} {virus} MeanDepth={meandepth} <1, using min_depth={depth_floor}")
            tasks.append((sample, virus, depth_floor))
    logger.info(f"Loaded {len(tasks)} tasks from {summary_file}")
    return tasks


def extract_virus_reference(virus_name, ref_fasta, out_dir):
    """
    提取病毒参考序列到 out_dir/{virus_name}.ref.fasta，若已存在则跳过
    返回 (success, ref_path)
    """
    out_file = os.path.join(out_dir, f"{virus_name}.ref.fasta")
    if os.path.isfile(out_file):
        logger.debug(f"Reference for {virus_name} already exists: {out_file}")
        return True, out_file
    cmd = f"seqkit grep -p '{virus_name}' {ref_fasta} -o {out_file}"
    logger.info(f"Extracting reference for {virus_name}")
    ret = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if ret.returncode != 0:
        logger.error(f"Failed to extract reference for {virus_name}: {ret.stderr}")
        return False, None
    return True, out_file


def process_sample(sample, virus, min_depth, bam_dir, virus_ref_file, out_dir):
    """
    处理单个样本-病毒对，生成 consensus 序列
    返回 (success, sample, virus, elapsed_time)
    """
    start_time = time.time()
    fixed_bam = os.path.join(out_dir, f"{virus}.{sample}.fixed.bam")
    consensus_fa = os.path.join(out_dir, f"{virus}.{sample}.consensus.fasta")

    # 如果 consensus 已存在，直接返回成功（耗时记为0）
    if os.path.isfile(consensus_fa):
        logger.info(f"Consensus already exists: {consensus_fa}, skipping")
        return True, sample, virus, 0.0

    bam_file = os.path.join(bam_dir, f"{sample}.sorted.bam")
    if not os.path.isfile(bam_file):
        logger.error(f"BAM file not found: {bam_file}")
        return False, sample, virus, time.time() - start_time

    # 步骤1: 提取该病毒的比对，并过滤 header
    awk_cmd = f"awk '{{if(/^@SQ/ && $2!=\"SN:{virus}\") next; print}}'"
    cmd1 = f"samtools view -h {bam_file} {virus} | {awk_cmd} | samtools view -b -o {fixed_bam}"
    logger.info(f"Generating fixed BAM for {sample} {virus}")
    ret1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    if ret1.returncode != 0:
        logger.error(f"samtools view failed for {sample} {virus}: {ret1.stderr}")
        return False, sample, virus, time.time() - start_time

    # 步骤2: 运行 viral_consensus
    cmd2 = f"viral_consensus -i {fixed_bam} -r {virus_ref_file} -o {consensus_fa} --min_depth {min_depth}"
    logger.info(f"Running viral_consensus for {sample} {virus}")
    ret2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
    if ret2.returncode != 0:
        logger.error(f"viral_consensus failed for {sample} {virus}: {ret2.stderr}")
        if os.path.isfile(consensus_fa) and os.path.getsize(consensus_fa) == 0:
            os.remove(consensus_fa)
        return False, sample, virus, time.time() - start_time

    elapsed = time.time() - start_time
    logger.info(f"Successfully generated consensus: {consensus_fa} (took {elapsed:.2f}s)")
    return True, sample, virus, elapsed


def main():
    global args
    args = parse_arguments()

    result_dir = os.path.abspath(args.result)
    bam_dir = os.path.join(result_dir, 'bam')
    summary_dir = os.path.join(result_dir, 'summary')
    consensus_root = os.path.join(result_dir, 'consensus')

    # 设置文件日志
    setup_file_logging(result_dir)
    logger.info("=== Starting batch_viral_consensus.py ===")
    logger.info(f"Result directory: {result_dir}")
    logger.info(f"Reference: {args.reference}")
    logger.info(f"Jobs: {args.jobs}")

    # 检查必要目录
    for d in [bam_dir, summary_dir]:
        if not os.path.isdir(d):
            logger.error(f"Required directory not found: {d}")
            sys.exit(1)
    os.makedirs(consensus_root, exist_ok=True)

    # 记录总开始时间
    total_start = time.time()

    # 确定 summary 文件：优先使用 all_viruses.summary.tsv
    all_summary = os.path.join(summary_dir, 'all_viruses.summary.tsv')
    if os.path.isfile(all_summary):
        tasks = read_summary_tasks(all_summary)
    else:
        tasks = []
        for sum_file in glob.glob(os.path.join(summary_dir, '*.summary.tsv')):
            if os.path.basename(sum_file) == 'all_viruses.summary.tsv':
                continue
            tasks.extend(read_summary_tasks(sum_file))

    if not tasks:
        logger.error("No tasks found in summary files")
        sys.exit(1)

    viruses = set(task[1] for task in tasks)
    logger.info(f"Found {len(viruses)} viruses: {', '.join(viruses)}")

    # --- 第一步：为每个病毒提取参考序列 ---
    virus_refs = {}
    virus_dirs = {}
    logger.info("Extracting virus reference sequences...")
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_virus = {}
        for virus in viruses:
            virus_dir = os.path.join(consensus_root, f"{virus}_consensus")
            os.makedirs(virus_dir, exist_ok=True)
            virus_dirs[virus] = virus_dir
            future = executor.submit(extract_virus_reference, virus, args.reference, virus_dir)
            future_to_virus[future] = virus

        # 使用 tqdm 显示进度，并输出错误到控制台（通过 pbar.write）
        with tqdm(total=len(future_to_virus), desc="Extracting references", unit="virus", position=0) as pbar:
            for future in as_completed(future_to_virus):
                virus = future_to_virus[future]
                try:
                    success, ref_path = future.result()
                    if success:
                        virus_refs[virus] = ref_path
                    else:
                        pbar.write(f"ERROR: Failed to extract reference for {virus}")
                except Exception as e:
                    pbar.write(f"ERROR: Exception while extracting {virus}: {e}")
                pbar.update(1)

    # 过滤掉参考序列提取失败的任务
    valid_tasks = [t for t in tasks if t[1] in virus_refs]
    logger.info(f"Processing {len(valid_tasks)} valid sample-virus pairs")

    # --- 第二步：并行生成 consensus ---
    success_count = 0
    fail_count = 0
    failed_tasks = []          # 存放 (sample, virus)
    success_times = []          # 存放成功任务的耗时

    logger.info("Generating consensus sequences...")
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_task = {}
        for sample, virus, min_depth in valid_tasks:
            virus_dir = virus_dirs[virus]
            ref_file = virus_refs[virus]
            future = executor.submit(process_sample, sample, virus, min_depth,
                                     bam_dir, ref_file, virus_dir)
            future_to_task[future] = (sample, virus)

        with tqdm(total=len(future_to_task), desc="Consensus generation", unit="sample", position=0) as pbar:
            for future in as_completed(future_to_task):
                sample, virus = future_to_task[future]
                try:
                    success, s, v, elapsed = future.result()
                    if success:
                        success_count += 1
                        success_times.append(elapsed)
                    else:
                        fail_count += 1
                        failed_tasks.append((s, v))
                except Exception as e:
                    pbar.write(f"ERROR: Exception processing {sample} {virus}: {e}")
                    fail_count += 1
                    failed_tasks.append((sample, virus))
                pbar.update(1)

    # --- 最终统计 ---
    total_time = time.time() - total_start
    avg_time = sum(success_times) / len(success_times) if success_times else 0

    # 输出到控制台（此时进度条已结束）
    print("\n" + "="*60)
    print(f"Total runtime: {total_time:.2f} seconds")
    print(f"Successful tasks: {success_count}, Failed tasks: {fail_count}")
    if success_times:
        print(f"Average time per successful task: {avg_time:.2f} seconds")
    if failed_tasks:
        print("\nFailed tasks (Sample, Virus):")
        for s, v in failed_tasks:
            print(f"  {s}\t{v}")
    print("="*60)

    logger.info(f"All done. Success: {success_count}, Failed: {fail_count}")


if __name__ == '__main__':
    main()
