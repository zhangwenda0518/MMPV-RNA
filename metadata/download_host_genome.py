#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
📥 宿主参考基因组下载器 v1.0
============================
基于 NCBI datasets CLI，自动下载指定物种的：
  - 核基因组 (genome)
  - 注释文件 (GFF3)
  - 序列报告 (seq-report)

功能特性:
  - 自动解压与合并多文件
  - 序列名去重（处理多基因组合并时的重名问题）
  - 支持叶绿体/线粒体基因组的额外下载
  - 可选基因组大小过滤
  - 输出统一 FASTA 用于下游分析

用法:
  python download_host_genome.py --species "Lycium barbarum" --outdir ./host_genome
  python download_host_genome.py --species "Lycium barbarum" --include-organelles --ncbi-api xxx
"""

import os
import sys
import re
import gzip
import argparse
import subprocess
import hashlib
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Set, Tuple

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


# ==========================================
# 终端UI
# ==========================================
class UI:
    CYAN = '\033[96m'; GREEN = '\033[92m'; YELLOW = '\033[93m'
    RED = '\033[91m'; PURPLE = '\033[95m'; GRAY = '\033[90m'
    BOLD = '\033[1m'; RESET = '\033[0m'

    @staticmethod
    def ok(msg):    print(f"  {UI.GREEN}✓{UI.RESET} {msg}")
    @staticmethod
    def warn(msg):  print(f"  {UI.YELLOW}⚠{UI.RESET} {msg}")
    @staticmethod
    def err(msg):   print(f"  {UI.RED}✗{UI.RESET} {msg}")
    @staticmethod
    def info(msg):  print(f"  {UI.CYAN}→{UI.RESET} {msg}")
    @staticmethod
    def header(msg):
        print(f"\n{UI.PURPLE}{UI.BOLD}{'='*55}{UI.RESET}")
        print(f"{UI.PURPLE}{UI.BOLD} {msg}{UI.RESET}")
        print(f"{UI.PURPLE}{UI.BOLD}{'='*55}{UI.RESET}")


def check_datasets_cli() -> bool:
    """检查 NCBI datasets CLI 是否可用"""
    try:
        result = subprocess.run(
            ['datasets', '--version'],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            UI.ok(f"NCBI datasets CLI 可用: {version}")
            return True
    except FileNotFoundError:
        pass
    except Exception:
        pass

    UI.err("未找到 NCBI datasets CLI")
    UI.info("安装方法: conda install -c conda-forge ncbi-datasets-cli")
    UI.info("或访问: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
    return False


def open_fasta(path: str):
    """智能打开 FASTA 文件（支持 .gz 压缩）"""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, 'r', encoding='utf-8')


def collect_fasta_files(directory: str) -> List[str]:
    """递归收集目录下所有 FASTA 文件"""
    fasta_exts = {'.fna', '.fasta', '.fa', '.fna.gz', '.fasta.gz', '.fa.gz'}
    files = []
    for root, _, filenames in os.walk(directory):
        for f in filenames:
            ext = os.path.splitext(f)[1]
            if f.endswith('.gz'):
                # 检查 .fna.gz 等双重后缀
                base = f[:-3]
                if any(base.endswith(e) for e in ['.fna', '.fasta', '.fa']):
                    files.append(os.path.join(root, f))
            elif ext in {'.fna', '.fasta', '.fa'}:
                files.append(os.path.join(root, f))
    return sorted(files)


def count_sequences(fasta_path: str) -> int:
    """快速统计 FASTA 文件中的序列数"""
    count = 0
    try:
        with open_fasta(fasta_path) as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception:
        pass
    return count


def merge_and_deduplicate(
    fasta_files: List[str],
    output_path: str,
    min_length: int = 0,
    prefix: str = ''
) -> Tuple[int, int]:
    """
    合并多个 FASTA 文件并去重序列名

    返回值: (写入序列数, 重复序列名数)
    """
    UI.header("合并与去重")

    seen_ids: Set[str] = set()
    total_written = 0
    duplicate_count = 0
    total_bp = 0

    with open(output_path, 'w', encoding='utf-8') as out_f:
        out_f.write(f"# Merged host genome\n")
        out_f.write(f"# Created: {datetime.now().isoformat()}\n")
        out_f.write(f"# Source files: {len(fasta_files)}\n")
        out_f.write(f"#\n")

        for fa_path in fasta_files:
            basename = os.path.basename(fa_path)
            UI.info(f"处理: {basename}")
            local_count = 0

            try:
                with open_fasta(fa_path) as in_f:
                    current_seq = []
                    current_id = ''
                    current_full_header = ''

                    for line in in_f:
                        if line.startswith('>'):
                            # 写入上一条序列
                            if current_id and current_seq:
                                seq = ''.join(current_seq)
                                if len(seq) >= min_length:
                                    out_f.write(f">{current_full_header}\n")
                                    for i in range(0, len(seq), 60):
                                        out_f.write(seq[i:i+60] + '\n')
                                    total_written += 1
                                    total_bp += len(seq)

                            # 解析新序列头
                            full_header = line[1:].strip()
                            seq_id = full_header.split()[0]

                            # 序列名去重
                            if seq_id in seen_ids:
                                suffix = 1
                                while f"{seq_id}_dup{suffix}" in seen_ids:
                                    suffix += 1
                                new_id = f"{seq_id}_dup{suffix}"
                                seen_ids.add(new_id)
                                current_full_header = full_header.replace(seq_id, new_id, 1)
                                duplicate_count += 1
                            else:
                                seen_ids.add(seq_id)
                                current_full_header = full_header

                            current_id = seq_id
                            current_seq = []
                            local_count += 1
                        else:
                            current_seq.append(line.strip())

                    # 处理最后一条序列
                    if current_id and current_seq:
                        seq = ''.join(current_seq)
                        if len(seq) >= min_length:
                            out_f.write(f">{current_full_header}\n")
                            for i in range(0, len(seq), 60):
                                out_f.write(seq[i:i+60] + '\n')
                            total_written += 1
                            total_bp += len(seq)

                UI.ok(f"  {local_count} 条序列")

            except Exception as e:
                UI.err(f"  处理 {basename} 时出错: {e}")

    file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
    UI.ok(f"合并完成: {total_written} 条序列, {total_bp:,} bp, {file_size_mb:.1f} MB")
    if duplicate_count:
        UI.warn(f"重命名了 {duplicate_count} 个重复序列名")

    return total_written, duplicate_count


def generate_genome_report(output_dir: str, fasta_files: List[str], merged_fasta: str):
    """生成基因组下载汇总报告"""
    report_path = os.path.join(output_dir, 'genome_report.txt')

    total_seqs = count_sequences(merged_fasta)
    total_size = os.path.getsize(merged_fasta) if os.path.isfile(merged_fasta) else 0

    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("  宿主参考基因组下载报告\n")
        f.write("=" * 60 + "\n")
        f.write(f"  生成时间: {datetime.now().isoformat()}\n")
        f.write(f"  合并文件: {merged_fasta}\n")
        f.write(f"  序列总数: {total_seqs}\n")
        f.write(f"  文件大小: {total_size / (1024*1024):.1f} MB\n")
        f.write(f"  源文件数: {len(fasta_files)}\n")
        f.write("\n  源文件列表:\n")
        for fa in fasta_files:
            f.write(f"    - {fa}\n")
        f.write("=" * 60 + "\n")

    UI.ok(f"报告已保存: {report_path}")


def download_organelle_genome(species: str, organelle: str, out_dir: str,
                               ncbi_api: str = '') -> Optional[str]:
    """
    下载细胞器基因组（叶绿体/线粒体）

    参数:
        species: 物种名
        organelle: 'chloroplast' 或 'mitochondrion'
        out_dir: 输出目录
    """
    UI.info(f"正在下载{organelle}基因组...")

    org_dir = os.path.join(out_dir, organelle)
    os.makedirs(org_dir, exist_ok=True)

    # 使用 NCBI E-utilities 搜索
    if organelle == 'chloroplast':
        search_term = f'"{species}"[Organism] AND chloroplast[filter]'
    else:
        search_term = f'"{species}"[Organism] AND mitochondrion[filter]'

    try:
        import requests
        api_key_param = f'&api_key={ncbi_api}' if ncbi_api else ''
        esearch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={search_term}&retmode=json{api_key_param}'
        resp = requests.get(esearch_url, timeout=30)
        data = resp.json()
        id_list = data.get('esearchresult', {}).get('idlist', [])

        if not id_list:
            UI.warn(f"  未找到 {species} 的 {organelle} 基因组")
            return None

        UI.ok(f"  找到 {len(id_list)} 条 {organelle} 记录")

        # 下载 FASTA
        ids_str = ','.join(id_list[:100])  # 限制100条
        efetch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={ids_str}&rettype=fasta&retmode=text{api_key_param}'
        resp = requests.get(efetch_url, timeout=120)
        fasta_content = resp.text

        out_path = os.path.join(org_dir, f'{organelle}.fasta')
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write(fasta_content)

        seq_count = fasta_content.count('>')
        UI.ok(f"  下载完成: {out_path} ({seq_count} 条序列)")
        return out_path

    except Exception as e:
        UI.err(f"  {organelle} 基因组下载失败: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="📥 宿主参考基因组下载器 — NCBI datasets 封装",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 基础下载
  python download_host_genome.py --species "Lycium barbarum" --outdir ./host_genome

  # 包含细胞器基因组 + API key
  python download_host_genome.py --species "Lycium barbarum" \\
      --outdir ./host_genome --include-organelles --ncbi-api xxx

  # 仅验证/检查已有下载
  python download_host_genome.py --species "Lycium barbarum" \\
      --outdir ./host_genome --verify-only
        """
    )

    parser.add_argument('--species', required=True, help='物种拉丁学名')
    parser.add_argument('--outdir', required=True, help='输出目录')
    parser.add_argument('--ncbi-api', help='NCBI API Key (提升速率)')
    parser.add_argument('--include-organelles', action='store_true',
                        help='同时下载叶绿体和线粒体基因组')
    parser.add_argument('--min-length', type=int, default=0,
                        help='最小序列长度过滤 (bp)')
    parser.add_argument('--verify-only', action='store_true',
                        help='仅验证现有下载')
    parser.add_argument('--skip-datasets', action='store_true',
                        help='跳过 NCBI datasets 步骤 (使用已有文件)')

    args = parser.parse_args()

    UI.header(f"宿主参考基因组下载: {args.species}")

    os.makedirs(args.outdir, exist_ok=True)

    # 验证模式
    if args.verify_only:
        UI.info("验证模式: 检查已有基因组文件...")
        merged_fasta = os.path.join(args.outdir, 'all.genome.uniq.fasta')
        if os.path.isfile(merged_fasta):
            n_seqs = count_sequences(merged_fasta)
            size_mb = os.path.getsize(merged_fasta) / (1024 * 1024)
            UI.ok(f"已存在合并基因组: {n_seqs} 条序列, {size_mb:.1f} MB")
        else:
            extracted_dir = os.path.join(args.outdir, 'extracted')
            fasta_files = collect_fasta_files(extracted_dir)
            if fasta_files:
                UI.ok(f"找到 {len(fasta_files)} 个源 FASTA 文件，但尚未合并")
                for f in fasta_files:
                    UI.info(f"  {f}")
            else:
                UI.warn("未找到任何基因组文件")
        return

    # 步骤 1: 使用 NCBI datasets 下载
    if not args.skip_datasets:
        if not check_datasets_cli():
            UI.warn("datasets CLI 不可用，尝试回退方案...")
        else:
            genome_zip = os.path.join(args.outdir, 'genome_down.zip')

            if os.path.isfile(genome_zip) and os.path.getsize(genome_zip) > 1000:
                UI.ok(f"基因组压缩包已存在: {genome_zip}")
            else:
                UI.header("步骤 1/3: NCBI datasets 下载")
                cmd = [
                    'datasets', 'download', 'genome', 'taxon', args.species,
                    '--filename', genome_zip,
                    '--include', 'genome,gff3,seq-report'
                ]
                if args.ncbi_api:
                    cmd.extend(['--api-key', args.ncbi_api])

                UI.info(f"执行: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

                if result.returncode != 0:
                    UI.err(f"datasets 下载失败: {result.stderr[:500]}")
                    UI.info("尝试使用 E-utilities 直接搜索...")
                    # 回退方案：直接用 E-utilities
                    try:
                        import requests
                        search_term = f'"{args.species}"[Organism] AND refseq[filter]'
                        api_key_param = f'&api_key={args.ncbi_api}' if args.ncbi_api else ''
                        esearch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={search_term}&retmax=5&retmode=json{api_key_param}'
                        resp = requests.get(esearch_url, timeout=30)
                        genome_ids = resp.json().get('esearchresult', {}).get('idlist', [])
                        if genome_ids:
                            UI.ok(f"找到 {len(genome_ids)} 个组装: {genome_ids}")
                        else:
                            UI.err("E-utilities 也未找到组装")
                    except Exception as e:
                        UI.err(f"回退方案失败: {e}")
                else:
                    UI.ok("基因组下载成功")

    # 步骤 2: 解压
    UI.header("步骤 2/3: 解压基因组")
    genome_zip = os.path.join(args.outdir, 'genome_down.zip')
    extracted_dir = os.path.join(args.outdir, 'extracted')

    if os.path.isfile(genome_zip):
        if not os.path.isdir(extracted_dir) or not os.listdir(extracted_dir):
            os.makedirs(extracted_dir, exist_ok=True)
            UI.info("正在解压...")
            result = subprocess.run(
                ['unzip', '-o', genome_zip, '-d', extracted_dir],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                UI.ok("解压完成")
            else:
                UI.err(f"解压失败: {result.stderr[:200]}")
                # 尝试用 Python zipfile
                import zipfile
                try:
                    with zipfile.ZipFile(genome_zip, 'r') as zf:
                        zf.extractall(extracted_dir)
                    UI.ok("解压完成 (Python zipfile)")
                except Exception as e:
                    UI.err(f"Python 解压也失败: {e}")
        else:
            UI.ok("已解压，跳过")
    else:
        UI.warn(f"基因组压缩包不存在: {genome_zip}")

    # 步骤 3: 合并去重
    UI.header("步骤 3/3: 合并与去重")

    fasta_files = collect_fasta_files(extracted_dir)
    if not fasta_files:
        UI.warn("未在解压目录中找到 FASTA 文件，尝试全目录搜索...")
        fasta_files = collect_fasta_files(args.outdir)

    if fasta_files:
        UI.ok(f"找到 {len(fasta_files)} 个 FASTA 文件:")
        for f in fasta_files:
            n_seqs = count_sequences(f)
            UI.info(f"  {os.path.basename(f)} — {n_seqs} 条序列")

        merged_fasta = os.path.join(args.outdir, 'all.genome.uniq.fasta')
        n_written, n_dups = merge_and_deduplicate(
            fasta_files, merged_fasta,
            min_length=args.min_length
        )
        generate_genome_report(args.outdir, fasta_files, merged_fasta)
    else:
        UI.warn("未找到任何 FASTA 文件")

    # 可选: 细胞器基因组
    if args.include_organelles:
        UI.header("额外: 细胞器基因组下载")
        for org in ['chloroplast', 'mitochondrion']:
            org_fasta = download_organelle_genome(
                args.species, org, args.outdir, args.ncbi_api
            )
            if org_fasta:
                # 合并到主文件
                merged_fasta = os.path.join(args.outdir, 'all.genome.uniq.fasta')
                if os.path.isfile(merged_fasta) and os.path.isfile(org_fasta):
                    UI.info(f"合并 {org} 基因组到主文件...")
                    all_files = collect_fasta_files(args.outdir)
                    all_files = [f for f in all_files if f != merged_fasta]
                    merge_and_deduplicate(all_files, merged_fasta + '.tmp',
                                          min_length=args.min_length)
                    os.replace(merged_fasta + '.tmp', merged_fasta)
                    UI.ok(f"{org} 基因组已合并")

    UI.header("下载完成")
    merged_fasta = os.path.join(args.outdir, 'all.genome.uniq.fasta')
    if os.path.isfile(merged_fasta):
        n_seqs = count_sequences(merged_fasta)
        size_mb = os.path.getsize(merged_fasta) / (1024 * 1024)
        UI.ok(f"最终输出: {merged_fasta}")
        UI.ok(f"  {n_seqs} 条序列, {size_mb:.1f} MB")
    else:
        UI.warn("最终合并文件未生成，请检查错误日志")


if __name__ == '__main__':
    main()
