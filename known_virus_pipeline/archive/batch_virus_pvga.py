#!/usr/bin/env python3
"""
批量 pvga 分析流程
自动处理单/双端 FASTA/FASTQ 数据，合并 reads，运行 pvga 病毒基因组组装
"""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from tqdm import tqdm

# -------------------- 工具函数 --------------------
def check_command(cmd):
    """检查命令是否存在于 PATH 中"""
    if shutil.which(cmd) is None:
        raise RuntimeError(f"命令 '{cmd}' 未找到，请确保已安装并加入 PATH")

def is_gzipped(filepath):
    """通过文件头判断是否为 gzip 压缩文件"""
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def open_file(filepath, mode='rt'):
    """智能打开普通文件或 gzip 文件"""
    if is_gzipped(filepath):
        return gzip.open(filepath, mode)
    else:
        return open(filepath, mode)

def run_cmd(cmd, desc=None):
    """运行外部命令，若失败则抛出异常"""
    if desc:
        print(f"[运行] {desc}", file=sys.stderr)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[错误] 命令执行失败: {cmd}", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"命令返回非零状态码: {result.returncode}")
    return result

def convert_fastq_to_fasta(fastq_file, fasta_file):
    """使用 bbmap 的 reformat.sh 将 FASTQ 转为 FASTA"""
    check_command("reformat.sh")
    cmd = f"reformat.sh in={fastq_file} out={fasta_file} overwrite=t"
    run_cmd(cmd, f"转换 FASTQ -> FASTA: {fastq_file}")

def cat_files(input_files, output_file):
    """合并多个文件（支持 gzip 压缩文件）"""
    with open(output_file, 'wb') as out_f:
        for f in input_files:
            with open_file(f, 'rb') as in_f:
                shutil.copyfileobj(in_f, out_f)

# -------------------- 样本识别 --------------------
def find_paired_samples(input_dir, pattern_r1=r'_(R?1|1)(?:_clean)?\.(fast[aq]|f[aq])(?:\.gz)?$',
                        pattern_r2=r'_(R?2|2)(?:_clean)?\.(fast[aq]|f[aq])(?:\.gz)?$'):
    """
    扫描输入目录，识别配对的样本
    返回 dict: {样本前缀: {"r1": path, "r2": path, "format": "fastq"/"fasta", "compressed": bool}}
    """
    input_path = Path(input_dir)
    files = list(input_path.glob("*"))
    samples = defaultdict(dict)

    for f in files:
        if not f.is_file():
            continue
        name = f.name
        # 匹配 R1 或 _1
        match_r1 = re.search(pattern_r1, name, re.IGNORECASE)
        if match_r1:
            prefix = name[:match_r1.start()]
            # 提取格式和压缩信息
            ext = match_r1.group(2).lower()
            fmt = "fastq" if "q" in ext else "fasta"
            compressed = name.endswith(".gz")
            samples[prefix]["r1"] = str(f)
            samples[prefix]["format"] = fmt
            samples[prefix]["compressed"] = compressed
            continue

        # 匹配 R2 或 _2
        match_r2 = re.search(pattern_r2, name, re.IGNORECASE)
        if match_r2:
            prefix = name[:match_r2.start()]
            samples[prefix]["r2"] = str(f)
            continue

        # 单端文件（无 R1/R2 标识）
        # 这里假设单端文件直接作为样本名
        # 但为避免误判，仅当文件匹配常见序列扩展名时才视为单端
        if re.search(r'\.(fast[aq]|f[aq])(?:\.gz)?$', name, re.IGNORECASE):
            single_name = re.sub(r'\.(fast[aq]|f[aq])(?:\.gz)?$', '', name)
            if "r1" not in samples[single_name] and "r2" not in samples[single_name]:
                samples[single_name]["r1"] = str(f)   # 将单端文件视为 R1
                fmt_match = re.search(r'\.(fast[aq]|f[aq])', name, re.IGNORECASE)
                fmt = "fastq" if fmt_match and "q" in fmt_match.group(1).lower() else "fasta"
                samples[single_name]["format"] = fmt
                samples[single_name]["compressed"] = name.endswith(".gz")
                samples[single_name]["single"] = True

    # 清理：仅保留有文件的条目
    valid_samples = {}
    for prefix, data in samples.items():
        if data.get("r1"):
            # 若只有 R1 没有 R2，则为单端
            if "r2" not in data:
                data["single"] = True
            valid_samples[prefix] = data

    return valid_samples

# -------------------- 单个样本处理流程 --------------------
def process_sample(sample_name, files_info, output_dir, pvga_ref, pvga_db, pvga_threads, keep_tmp):
    """
    处理单个样本：
    1. 若双端 FASTQ -> bbmerge -> cat merged+unmerged -> 转 FASTA
    2. 若双端 FASTA -> 直接 cat 两个文件
    3. 若单端 FASTQ -> 转 FASTA
    4. 若单端 FASTA -> 直接使用
    5. 运行 pvga
    """
    sample_out_dir = Path(output_dir) / sample_name
    sample_out_dir.mkdir(parents=True, exist_ok=True)

    r1 = files_info["r1"]
    r2 = files_info.get("r2")
    is_single = files_info.get("single", False)
    fmt = files_info["format"]
    compressed = files_info["compressed"]

    final_fasta = sample_out_dir / f"{sample_name}.final.fasta"
    log_file = sample_out_dir / "pipeline.log"

    try:
        with open(log_file, 'w') as log:
            log.write(f"样本: {sample_name}\n")
            log.write(f"R1: {r1}\n")
            if r2:
                log.write(f"R2: {r2}\n")
            log.write(f"格式: {fmt}, 压缩: {compressed}\n")

            # ---------- 合并步骤 ----------
            if not is_single and r2:
                if fmt == "fastq":
                    # FASTQ 双端：bbmerge
                    check_command("bbmerge.sh")
                    merged_fq = sample_out_dir / "merged.fastq"
                    unmerged1_fq = sample_out_dir / "unmerged1.fastq"
                    unmerged2_fq = sample_out_dir / "unmerged2.fastq"
                    all_fq = sample_out_dir / "all.fastq"

                    cmd_bbmerge = (f"bbmerge.sh in1={r1} in2={r2} out={merged_fq} "
                                   f"outu1={unmerged1_fq} outu2={unmerged2_fq} overwrite=t")
                    run_cmd(cmd_bbmerge, f"bbmerge 合并 {sample_name}")
                    # 合并所有 reads
                    cat_files([merged_fq, unmerged1_fq, unmerged2_fq], all_fq)
                    # 转为 FASTA
                    convert_fastq_to_fasta(all_fq, final_fasta)
                    log.write("合并方式: bbmerge (FASTQ)\n")

                else:  # FASTA
                    # 双端 FASTA 直接拼接
                    cat_files([r1, r2], final_fasta)
                    log.write("合并方式: 直接拼接两个 FASTA 文件\n")

            else:
                # 单端
                if fmt == "fastq":
                    # 单端 FASTQ 转为 FASTA
                    convert_fastq_to_fasta(r1, final_fasta)
                    log.write("单端 FASTQ，已转换为 FASTA\n")
                else:
                    # 单端 FASTA，直接复制或解压
                    if compressed:
                        with open_file(r1, 'rb') as src, open(final_fasta, 'wb') as dst:
                            shutil.copyfileobj(src, dst)
                    else:
                        shutil.copy(r1, final_fasta)
                    log.write("单端 FASTA，直接使用\n")

            # ---------- pvga 运行 ----------
            check_command("pvga")
            pvga_out = sample_out_dir / "pvga"
            cmd_pvga = f"pvga -r {final_fasta} -b {pvga_db} -n {pvga_threads} -o {pvga_out}"
            run_cmd(cmd_pvga, f"运行 pvga: {sample_name}")
            log.write(f"pvga 完成，输出目录: {pvga_out}\n")

            # 清理临时文件（可选）
            if not keep_tmp:
                for tmp in sample_out_dir.glob("*.fastq"):
                    tmp.unlink()
                for tmp in sample_out_dir.glob("unmerged*"):
                    tmp.unlink()

            return sample_name, True, None

    except Exception as e:
        return sample_name, False, str(e)

# -------------------- 主函数 --------------------
def main():
    parser = argparse.ArgumentParser(description="批量运行 pvga 病毒基因组组装流程")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="输入目录，包含测序 reads 文件（支持 .fastq/.fasta/.fa 及 .gz）")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="输出根目录，每个样本将在此目录下创建子文件夹")
    parser.add_argument("-b", "--pvga_db", required=True,
                        help="pvga 的参考数据库文件（-b 参数，例如 HXB2.fa）")
    parser.add_argument("-n", "--threads", type=int, default=10,
                        help="每个 pvga 任务使用的线程数（默认 10）")
    parser.add_argument("--jobs", type=int, default=1,
                        help="并行处理的样本数量（默认 1）")
    parser.add_argument("--keep_tmp", action="store_true",
                        help="保留中间文件（如合并后的 fastq）")
    parser.add_argument("--pattern_r1", default=r'_(R?1|1)(?:_clean)?\.(fast[aq]|f[aq])(?:\.gz)?$',
                        help="R1 文件匹配正则表达式（高级选项）")
    parser.add_argument("--pattern_r2", default=r'_(R?2|2)(?:_clean)?\.(fast[aq]|f[aq])(?:\.gz)?$',
                        help="R2 文件匹配正则表达式（高级选项）")
    args = parser.parse_args()

    # 检查必要命令
    for cmd in ["bbmerge.sh", "reformat.sh", "pvga"]:
        check_command(cmd)

    # 创建输出根目录
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # 扫描样本
    print("正在扫描输入目录，识别样本...")
    samples = find_paired_samples(args.input_dir, args.pattern_r1, args.pattern_r2)
    if not samples:
        print("错误：未找到任何 reads 文件。请检查输入目录和文件名模式。", file=sys.stderr)
        sys.exit(1)
    print(f"共发现 {len(samples)} 个样本")

    # 并行处理
    results = []
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {}
        for sample_name, files_info in samples.items():
            future = executor.submit(
                process_sample,
                sample_name,
                files_info,
                args.output_dir,
                None,               # pvga -r 参数直接使用 final_fasta
                args.pvga_db,
                args.threads,
                args.keep_tmp
            )
            futures[future] = sample_name

        # 使用 tqdm 显示进度
        with tqdm(total=len(futures), desc="处理进度") as pbar:
            for future in as_completed(futures):
                sample_name = futures[future]
                try:
                    name, success, error = future.result()
                    if success:
                        tqdm.write(f"✓ {name} 处理成功")
                    else:
                        tqdm.write(f"✗ {name} 处理失败: {error}")
                except Exception as e:
                    tqdm.write(f"✗ {sample_name} 处理异常: {e}")
                pbar.update(1)
                results.append((name, success))

    # 汇总
    succeeded = sum(1 for _, ok in results if ok)
    failed = len(results) - succeeded
    print(f"\n处理完成：成功 {succeeded} 个，失败 {failed} 个")
    if failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()
