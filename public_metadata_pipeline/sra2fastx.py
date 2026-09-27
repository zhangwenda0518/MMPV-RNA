#!/usr/bin/env python3
"""
SRA → FASTQ.GZ 批量转换器
===========================
递归扫描 .sra 文件，并行调用 fasterq-dump（失败回退 fastq-dump），
压缩为 .fastq.gz，成功后删除原 .sra 文件。

提供 CLI 入口和可编程 API（process_sra / convert_all）。
"""

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
from typing import List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


def convert_all(
    input_dir: str,
    jobs: int = 2,
    threads: int = 6,
    progress: bool = True,
) -> Tuple[int, List[str]]:
    """
    Convert all .sra files under *input_dir* to .fastq.gz.

    Returns:
        (success_count, list_of_error_messages)
    """
    input_path = Path(input_dir)
    sra_files = sorted(input_path.rglob("*.sra"))
    if not sra_files:
        return 0, []

    success_count = 0
    failed_messages: List[str] = []

    if progress:
        print(f"  找到 {len(sra_files)} 个 .sra 文件")
        print(f"  并行任务: {jobs}, 单任务线程: {threads}")
        if HAS_TQDM:
            print(f"  预估 CPU 占用: {jobs * threads} 核\n")

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        future_to_sra = {
            executor.submit(process_sra, sra, threads): sra
            for sra in sra_files
        }
        iterator = (
            tqdm(as_completed(future_to_sra), total=len(sra_files),
                 desc="转换进度", unit="文件")
            if HAS_TQDM and progress
            else as_completed(future_to_sra)
        )
        for future in iterator:
            success, msg = future.result()
            if success:
                success_count += 1
            else:
                failed_messages.append(msg)

    return success_count, failed_messages

def process_sra(sra_path: Path, threads: int) -> Tuple[bool, str]:
    """Convert a single .sra file to .fastq.gz. Returns (success, message)."""
    sra_path = Path(sra_path).resolve()
    sra_dir = sra_path.parent
    sra_prefix = sra_path.stem  # 获取文件名，例如 SRR12345

    # 1. 运行 fasterq-dump，失败回退 fastq-dump
    for tool in ["fasterq-dump", "fastq-dump"]:
        cmd_dump = [
            tool,
            "--split-3",
            "--threads", str(threads) if tool == "fasterq-dump" else "1",
            "-O", str(sra_dir),
            str(sra_path)
        ]
        try:
            subprocess.run(cmd_dump, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            break
        except subprocess.CalledProcessError:
            continue
    else:
        return False, f"转换失败 ({sra_path.name}): fasterq-dump & fastq-dump 均失败"

    # 2. 查找生成的 fastq 文件
    # --split-3 可能会生成 PREFIX.fastq, PREFIX_1.fastq, PREFIX_2.fastq 等
    fastq_files = list(sra_dir.glob(f"{sra_prefix}*.fastq"))
    
    if not fastq_files:
        return False, f"转换失败 ({sra_path.name}): 未找到生成的 fastq 文件。"

    # 3. 压缩 fastq 文件
    # 自动检测系统中是否有更快的 pigz (多线程gzip)，如果没有则使用普通 gzip
    zip_tool = "pigz" if shutil.which("pigz") else "gzip"
    
    for fq in fastq_files:
        cmd_zip = [zip_tool, "-f", str(fq)]
        try:
            subprocess.run(cmd_zip, check=True)
        except subprocess.CalledProcessError:
            return False, f"压缩失败 ({fq.name})"

    # 4. 所有步骤均成功，安全删除原始 SRA 文件
    try:
        sra_path.unlink()
    except OSError as e:
        return False, f"删除SRA文件失败 ({sra_path.name}): {e}"

    return True, f"成功: {sra_path.name}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="批量将SRA转换为fastq.gz。成功后自动删除原SRA文件。"
    )
    parser.add_argument("input_dir", help="包含 SRA 文件的根目录 (支持多级子目录)")
    parser.add_argument("-j", "--jobs", type=int, default=2, help="并行处理的SRA文件数量 (默认: 2)")
    parser.add_argument("-t", "--threads", type=int, default=6, help="每个 fasterq-dump 使用的线程数 (默认: 6)")

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        print(f"错误: 目录 '{input_dir}' 不存在。")
        sys.exit(1)

    success, failures = convert_all(args.input_dir, args.jobs, args.threads, progress=True)
    total = success + len(failures)

    if total == 0:
        print(f"在目录 {input_dir} 中未找到任何 .sra 文件。")
        sys.exit(0)

    print(f"\n任务完成! 成功处理 {success}/{total} 个文件。")
    if failures:
        print("以下文件在处理时遇到错误 (原 SRA 文件已被保留):")
        for err in failures:
            print(f"  - {err}")
        sys.exit(1)

if __name__ == "__main__":
    main()
