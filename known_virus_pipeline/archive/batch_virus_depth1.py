#!/usr/bin/env python3
"""
完整的病毒分析流程脚本
执行：比对 -> 覆盖度计算 -> 统计 -> FPKM/TPM计算
支持自动识别输入目录中的单双端测序数据
使用tqdm显示进度条，优化输出目录结构
"""

import argparse
import os
import sys
import subprocess
import gzip
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import logging
from datetime import datetime
import multiprocessing
from functools import partial
import tempfile
import re
from tqdm import tqdm
import shutil

# 设置日志
def setup_logging(log_level=logging.INFO, log_file=None, console_only=False):
    """设置日志系统"""
    logger = logging.getLogger()
    logger.setLevel(log_level)
    
    # 清除所有现有的处理器
    logger.handlers.clear()
    
    # 控制台处理器
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', 
                                       datefmt='%H:%M:%S')
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # 文件处理器（如果指定了日志文件且不是仅控制台模式）
    if log_file and not console_only:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_format = logging.Formatter('%(asctime)s - %(levelname)s - %(module)s - %(message)s')
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)
    
    return logger

class VirusAnalyzer:
    """病毒分析器类"""
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.temp_dir = tempfile.mkdtemp(prefix="virus_analysis_")
        self.logger.debug(f"临时目录: {self.temp_dir}")
        
        # 存储中间结果
        self.intermediate_files = []
    
    def check_dependencies(self):
        """检查所需工具是否可用"""
        required_tools = ['strobealign', 'samtools', 'pandepth']
        missing_tools = []
        
        self.logger.info("检查依赖工具...")
        
        for tool in required_tools:
            try:
                result = subprocess.run([tool, '--version'], 
                                      capture_output=True, 
                                      text=True,
                                      check=False)
                if result.returncode == 0:
                    version_line = result.stdout.split('\n')[0] if result.stdout else "unknown"
                    self.logger.info(f"  ✓ {tool}: {version_line}")
                else:
                    missing_tools.append(tool)
            except FileNotFoundError:
                missing_tools.append(tool)
        
        if missing_tools:
            self.logger.error(f"缺少必要工具: {', '.join(missing_tools)}")
            self.logger.error("请确保以下工具已安装并在PATH中:")
            self.logger.error("1. strobealign: https://github.com/ksahlin/strobealign")
            self.logger.error("2. samtools: https://github.com/samtools/samtools")
            self.logger.error("3. pandepth: https://github.com/schatzlab/pandepth")
            return False
        
        return True
    
    def discover_sequencing_files(self, input_dirs):
        """
        自动发现输入目录中的测序文件并配对
        
        Args:
            input_dirs: 输入目录列表
        
        Returns:
            list: 样本列表，每个元素为(样本名, R1路径, R2路径或None)
        """
        samples = []
        self.logger.info("扫描输入目录中的测序文件...")
        
        # 常见的测序文件扩展名
        seq_extensions = ['.fq.gz', '.fastq.gz', '.fq', '.fastq']
        
        all_seq_files = []
        
        # 首先收集所有文件
        with tqdm(total=len(input_dirs), desc="扫描目录", unit="目录") as pbar:
            for input_dir in input_dirs:
                input_path = Path(input_dir)
                if not input_path.exists():
                    self.logger.warning(f"输入目录不存在: {input_dir}")
                    pbar.update(1)
                    continue
                
                # 查找所有测序文件
                for ext in seq_extensions:
                    for file_path in input_path.rglob(f'*{ext}'):
                        all_seq_files.append(file_path)
                
                pbar.update(1)
        
        if not all_seq_files:
            self.logger.error("未找到测序文件")
            return samples
        
        self.logger.info(f"找到 {len(all_seq_files)} 个测序文件")
        
        # 构建文件名到路径的映射
        file_map = {}
        for file_path in all_seq_files:
            file_name = file_path.name
            file_map[file_name] = file_path
        
        # 使用正则表达式识别配对文件
        paired_files = defaultdict(list)
        
        # 定义配对模式
        patterns = [
            # 格式: sample_R1.fastq.gz / sample_R2.fastq.gz
            (r'^(.*?)[._-]R?1[._-]', r'^(.*?)[._-]R?2[._-]'),
            # 格式: sample_1.fastq.gz / sample_2.fastq.gz
            (r'^(.*?)[._-]1[._-]', r'^(.*?)[._-]2[._-]'),
            # 格式: sample.R1.fastq.gz / sample.R2.fastq.gz
            (r'^(.*?)\.R?1\.', r'^(.*?)\.R?2\.'),
        ]
        
        # 标记已处理的文件
        processed_files = set()
        paired_samples = {}
        
        # 先尝试配对双端测序
        for file_name, file_path in file_map.items():
            if file_name in processed_files:
                continue
                
            for r1_pattern, r2_pattern in patterns:
                r1_match = re.match(r1_pattern, file_name, re.IGNORECASE)
                if r1_match:
                    sample_base = r1_match.group(1)
                    
                    # 寻找对应的R2文件
                    for r2_name, r2_path in file_map.items():
                        if r2_name in processed_files:
                            continue
                            
                        if re.match(r2_pattern + r'.*', r2_name, re.IGNORECASE) and r2_name.startswith(sample_base):
                            # 找到配对
                            paired_samples[sample_base] = (str(file_path), str(r2_path))
                            processed_files.add(file_name)
                            processed_files.add(r2_name)
                            self.logger.debug(f"配对成功: {sample_base} - {file_name} / {r2_name}")
                            break
                    break
        
        # 处理未配对的单端文件
        for file_name, file_path in file_map.items():
            if file_name not in processed_files:
                # 提取样本名（去除扩展名）
                sample_base = file_name
                for ext in seq_extensions:
                    if sample_base.endswith(ext):
                        sample_base = sample_base[:-len(ext)]
                        break
                
                # 清理样本名
                sample_base = re.sub(r'[._-](R?[12]|read[12]|fq|fastq)$', '', sample_base, flags=re.IGNORECASE)
                sample_base = sample_base.rstrip('._-')
                
                paired_samples[sample_base] = (str(file_path), None)
                processed_files.add(file_name)
                self.logger.debug(f"单端文件: {sample_base} - {file_name}")
        
        # 构建最终样本列表
        for sample_base, (r1_path, r2_path) in paired_samples.items():
            # 确保样本名唯一
            unique_name = sample_base
            counter = 1
            while unique_name in [s[0] for s in samples]:
                unique_name = f"{sample_base}_{counter}"
                counter += 1
            
            samples.append((unique_name, r1_path, r2_path))
        
        # 统计结果
        paired_count = sum(1 for _, _, r2 in samples if r2 is not None)
        single_count = len(samples) - paired_count
        
        self.logger.info(f"样本发现完成: 共 {len(samples)} 个样本")
        self.logger.info(f"  - 双端测序: {paired_count} 个")
        self.logger.info(f"  - 单端测序: {single_count} 个")
        
        return samples
    
    def align_sample(self, sample_name, r1_path, r2_path, reference, threads=8, output_dir=None):
        """
        对单个样本进行比对
        
        Args:
            sample_name: 样本名称
            r1_path: R1 reads文件路径
            r2_path: R2 reads文件路径（单端测序时为None）
            reference: 参考基因组文件路径
            threads: 线程数
            output_dir: 输出目录
        
        Returns:
            str: 输出BAM文件路径
        """
        # 确定输出目录
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        # 输出文件路径
        bam_file = output_dir / f"{sample_name}.bam"
        
        # 检查文件是否已存在
        if bam_file.exists():
            self.logger.debug(f"BAM文件已存在，跳过比对: {bam_file}")
            return str(bam_file)
        
        # 构建比对命令
        if r2_path:
            # 双端测序
            align_cmd = [
                'strobealign',
                '-t', str(threads),
                str(reference),
                str(r1_path),
                str(r2_path)
            ]
            seq_type = "双端"
        else:
            # 单端测序
            align_cmd = [
                'strobealign',
                '-t', str(threads),
                str(reference),
                str(r1_path)
            ]
            seq_type = "单端"
        
        self.logger.debug(f"比对样本 {sample_name} ({seq_type}测序)")
        
        try:
            # 执行比对和排序管道
            with open(os.devnull, 'w') as devnull:
                align_proc = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=devnull)
                
                # samtools排序命令
                sort_cmd = [
                    'samtools', 'sort',
                    '-@', str(min(threads, 20)),
                    '-o', str(bam_file)
                ]
                sort_proc = subprocess.Popen(sort_cmd, stdin=align_proc.stdout, stderr=devnull)
                
                # 等待进程完成
                align_proc.stdout.close()
                sort_proc.communicate()
            
            if sort_proc.returncode != 0:
                raise subprocess.CalledProcessError(sort_proc.returncode, sort_cmd)
            
            # 创建BAM索引
            index_cmd = ['samtools', 'index', str(bam_file)]
            subprocess.run(index_cmd, check=True, capture_output=True)
            
            # 记录中间文件
            self.intermediate_files.append(str(bam_file))
            self.intermediate_files.append(str(bam_file) + ".bai")
            
            return str(bam_file)
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"比对过程失败: {e}")
            # 清理失败的输出文件
            if bam_file.exists():
                bam_file.unlink()
            raise
    
    def calculate_coverage(self, bam_file, threads=10, output_prefix=None):
        """
        使用pandepth计算覆盖度
        
        Args:
            bam_file: BAM文件路径
            threads: 线程数
            output_prefix: 输出文件前缀
        
        Returns:
            str: 统计文件路径
        """
        bam_path = Path(bam_file)
        sample_name = bam_path.stem
        
        if output_prefix is None:
            output_prefix = bam_path.parent / sample_name
        
        # 检查文件是否已存在
        stat_file = Path(f"{output_prefix}.chr.stat.gz")
        if stat_file.exists():
            self.logger.debug(f"统计文件已存在，跳过计算: {stat_file}")
            return str(stat_file)
        
        # pandepth命令
        pandepth_cmd = [
            'pandepth',
            '-a',
            '-i', str(bam_file),
            '-o', str(output_prefix),
            '-t', str(threads)
        ]
        
        try:
            self.logger.debug(f"计算覆盖度: {sample_name}")
            
            # 执行pandepth
            with open(os.devnull, 'w') as devnull:
                subprocess.run(pandepth_cmd, check=True, stderr=devnull)
            
            # 检查输出文件
            if not stat_file.exists():
                # 尝试其他可能的文件名
                possible_files = [
                    Path(f"{output_prefix}.stat.gz"),
                    Path(f"{output_prefix}.chr.stat"),
                    Path(f"{output_prefix}.stat")
                ]
                
                for file in possible_files:
                    if file.exists():
                        stat_file = file
                        break
                else:
                    raise FileNotFoundError(f"未找到pandepth输出文件: {output_prefix}.*")
            
            # 记录中间文件
            self.intermediate_files.append(str(stat_file))
            
            return str(stat_file)
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"pandepth执行失败: {e}")
            raise
    
    def parse_stat_file(self, stat_file, sample_name, coverage_threshold=90.0, depth_threshold=10.0):
        """
        解析stat文件并计算统计信息
        
        Args:
            stat_file: stat文件路径
            sample_name: 样本名称
            coverage_threshold: Coverage阈值
            depth_threshold: MeanDepth阈值
        
        Returns:
            list: 满足条件的病毒信息列表
        """
        results = []
        
        try:
            # 支持gzip压缩和普通文本文件
            if str(stat_file).endswith('.gz'):
                open_func = gzip.open
                mode = 'rt'
            else:
                open_func = open
                mode = 'r'
            
            with open_func(stat_file, mode) as f:
                for line in f:
                    line = line.strip()
                    
                    # 跳过注释行和空行
                    if not line or line.startswith('#') or line.startswith('##'):
                        continue
                    
                    parts = line.split()
                    if len(parts) < 6:
                        continue
                    
                    virus_name, length, covered_site, total_depth, coverage, mean_depth = parts
                    
                    try:
                        # 转换为适当的类型
                        length_val = int(length)
                        covered_site_val = int(covered_site)
                        total_depth_val = int(total_depth)
                        coverage_val = float(coverage)
                        mean_depth_val = float(mean_depth)
                        
                        # 检查是否满足条件
                        if coverage_val > coverage_threshold and mean_depth_val > depth_threshold:
                            results.append({
                                'Sample': sample_name,
                                'Virus': virus_name,
                                'Length': length_val,
                                'CoveredSite': covered_site_val,
                                'TotalDepth': total_depth_val,
                                'Coverage(%)': coverage_val,
                                'MeanDepth': mean_depth_val
                            })
                    except ValueError:
                        continue
                        
        except Exception as e:
            self.logger.error(f"解析stat文件失败 {stat_file}: {e}")
        
        return results
    
    def calculate_fpkm_tpm(self, stat_file, idxstats_file=None, bam_file=None):
        """
        计算FPKM和TPM值
        
        Args:
            stat_file: stat文件路径
            idxstats_file: samtools idxstats输出文件路径（可选）
            bam_file: BAM文件路径（如果未提供idxstats_file）
        
        Returns:
            tuple: (total_mapped_reads, virus_read_counts)
        """
        virus_reads = {}
        
        try:
            # 获取idxstats文件
            if idxstats_file is None and bam_file is not None:
                idxstats_file = Path(bam_file).with_suffix('.idxstats')
                if not idxstats_file.exists():
                    idxstats_cmd = ['samtools', 'idxstats', bam_file]
                    result = subprocess.run(idxstats_cmd, capture_output=True, text=True, check=True)
                    idxstats_content = result.stdout
                else:
                    with open(idxstats_file, 'r') as f:
                        idxstats_content = f.read()
            elif idxstats_file is not None:
                with open(idxstats_file, 'r') as f:
                    idxstats_content = f.read()
            else:
                return 0, {}
            
            # 解析idxstats文件
            total_mapped = 0
            for line in idxstats_content.strip().split('\n'):
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 4:
                    virus_name = parts[0]
                    if virus_name == '*':  # 跳过未比对的reads
                        continue
                    
                    mapped_reads = int(parts[2])
                    virus_reads[virus_name] = mapped_reads
                    total_mapped += mapped_reads
            
        except Exception as e:
            self.logger.debug(f"计算FPKM/TPM时出错: {e}")
            return 0, {}
        
        return total_mapped, virus_reads
    
    def process_single_sample(self, sample_info, reference, output_dir, 
                              coverage_threshold, depth_threshold, threads,
                              skip_fpkm_tpm=False, keep_intermediate=False):
        """处理单个样本的完整流程"""
        sample_name, r1_path, r2_path = sample_info
        
        try:
            # 1. 比对
            bam_file = self.align_sample(
                sample_name=sample_name,
                r1_path=r1_path,
                r2_path=r2_path,
                reference=reference,
                threads=threads,
                output_dir=output_dir / "alignment"
            )
            
            # 2. 计算覆盖度
            stat_file = self.calculate_coverage(
                bam_file=bam_file,
                threads=threads,
                output_prefix=output_dir / "coverage" / sample_name
            )
            
            # 3. 解析统计文件并计算FPKM/TPM
            results = self.parse_stat_file(
                stat_file=stat_file,
                sample_name=sample_name,
                coverage_threshold=coverage_threshold,
                depth_threshold=depth_threshold
            )
            
            # 4. 计算FPKM和TPM
            total_mapped = 0
            if not skip_fpkm_tpm and results:
                total_mapped, virus_reads = self.calculate_fpkm_tpm(stat_file, bam_file=bam_file)
                
                # 为每个结果添加FPKM和TPM
                if total_mapped > 0 and virus_reads:
                    # 先计算所有病毒的RPK（Reads Per Kilobase）
                    virus_rpk = {}
                    total_rpk = 0.0
                    
                    for virus_name, reads in virus_reads.items():
                        # 从结果中获取长度
                        virus_info = next((r for r in results if r['Virus'] == virus_name), None)
                        if virus_info and virus_info['Length'] > 0:
                            length_kb = virus_info['Length'] / 1000.0
                            rpk = reads / length_kb
                            virus_rpk[virus_name] = rpk
                            total_rpk += rpk
                    
                    # 为每个结果计算FPKM和TPM
                    for result in results:
                        virus_name = result['Virus']
                        reads = virus_reads.get(virus_name, 0)
                        length_kb = result['Length'] / 1000.0
                        
                        # 计算FPKM
                        if total_mapped > 0 and length_kb > 0:
                            result['FPKM'] = (reads * 1e9) / (total_mapped * length_kb)
                        else:
                            result['FPKM'] = 0.0
                        
                        # 计算TPM
                        if total_rpk > 0 and virus_name in virus_rpk:
                            result['TPM'] = (virus_rpk[virus_name] * 1e6) / total_rpk
                        else:
                            result['TPM'] = 0.0
            
            # 5. 清理中间文件（如果不需要保留）
            if not keep_intermediate:
                self.cleanup_intermediate_files()
            
            return sample_name, results, total_mapped
            
        except Exception as e:
            self.logger.error(f"处理样本 {sample_name} 时出错: {e}")
            return sample_name, [], 0
    
    def analyze_samples(self, sample_list, reference, output_dir='./results',
                       coverage_threshold=90.0, depth_threshold=10.0,
                       threads=8, parallel=1, skip_alignment=False,
                       skip_fpkm_tpm=False, keep_intermediate=False):
        """
        分析多个样本
        
        Args:
            sample_list: 样本列表
            reference: 参考基因组文件路径
            output_dir: 输出目录
            coverage_threshold: Coverage阈值
            depth_threshold: MeanDepth阈值
            threads: 每个样本使用的线程数
            parallel: 并行处理的样本数
            skip_alignment: 是否跳过比对步骤
            skip_fpkm_tpm: 是否跳过FPKM/TPM计算
            keep_intermediate: 是否保留中间文件
        """
        output_dir = Path(output_dir)
        
        # 创建目录结构
        self.create_directory_structure(output_dir, keep_intermediate)
        
        self.logger.info(f"开始分析 {len(sample_list)} 个样本")
        self.logger.info(f"输出目录: {output_dir}")
        self.logger.info(f"阈值: Coverage > {coverage_threshold}%, MeanDepth > {depth_threshold}")
        
        all_results = []
        virus_results = defaultdict(list)
        sample_stats = {}
        
        if skip_alignment:
            # 跳过比对步骤，直接分析已有BAM文件
            self.logger.info("跳过比对步骤，直接分析已有BAM文件")
            
            for sample_name, r1_path, r2_path in tqdm(sample_list, desc="分析样本", unit="样本"):
                # 假设r1_path是BAM文件路径
                bam_file = r1_path
                
                try:
                    # 计算覆盖度
                    stat_file = self.calculate_coverage(
                        bam_file=bam_file,
                        threads=threads,
                        output_prefix=output_dir / "coverage" / sample_name
                    )
                    
                    # 解析统计文件
                    results = self.parse_stat_file(
                        stat_file=stat_file,
                        sample_name=sample_name,
                        coverage_threshold=coverage_threshold,
                        depth_threshold=depth_threshold
                    )
                    
                    # 计算FPKM和TPM
                    total_mapped = 0
                    if not skip_fpkm_tpm:
                        total_mapped, virus_reads = self.calculate_fpkm_tpm(stat_file, bam_file=bam_file)
                        
                        # 为每个结果添加FPKM和TPM
                        if total_mapped > 0 and virus_reads and results:
                            # 先计算所有病毒的RPK（Reads Per Kilobase）
                            virus_rpk = {}
                            total_rpk = 0.0
                            
                            for virus_name, reads in virus_reads.items():
                                # 从结果中获取长度
                                virus_info = next((r for r in results if r['Virus'] == virus_name), None)
                                if virus_info and virus_info['Length'] > 0:
                                    length_kb = virus_info['Length'] / 1000.0
                                    rpk = reads / length_kb
                                    virus_rpk[virus_name] = rpk
                                    total_rpk += rpk
                            
                            # 为每个结果计算FPKM和TPM
                            for result in results:
                                virus_name = result['Virus']
                                reads = virus_reads.get(virus_name, 0)
                                length_kb = result['Length'] / 1000.0
                                
                                # 计算FPKM
                                if total_mapped > 0 and length_kb > 0:
                                    result['FPKM'] = (reads * 1e9) / (total_mapped * length_kb)
                                else:
                                    result['FPKM'] = 0.0
                                
                                # 计算TPM
                                if total_rpk > 0 and virus_name in virus_rpk:
                                    result['TPM'] = (virus_rpk[virus_name] * 1e6) / total_rpk
                                else:
                                    result['TPM'] = 0.0
                    
                    sample_stats[sample_name] = total_mapped
                    all_results.extend(results)
                    
                    for result in results:
                        virus_name = result['Virus']
                        virus_results[virus_name].append(result)
                        
                except Exception as e:
                    self.logger.error(f"分析样本 {sample_name} 时出错: {e}")
        
        else:
            # 正常流程：比对 -> 覆盖度计算 -> 统计
            if parallel > 1:
                # 并行处理
                self.logger.info(f"使用 {parallel} 个进程并行处理")
                
                # 创建部分函数
                process_func = partial(
                    self.process_single_sample,
                    reference=reference,
                    output_dir=output_dir,
                    coverage_threshold=coverage_threshold,
                    depth_threshold=depth_threshold,
                    threads=max(1, threads // parallel),
                    skip_fpkm_tpm=skip_fpkm_tpm,
                    keep_intermediate=keep_intermediate
                )
                
                # 使用进程池并行处理
                with multiprocessing.Pool(processes=parallel) as pool:
                    # 使用tqdm显示进度
                    results = list(tqdm(
                        pool.imap(process_func, sample_list),
                        total=len(sample_list),
                        desc="处理样本",
                        unit="样本"
                    ))
                    
                    for sample_name, sample_results, total_mapped in results:
                        sample_stats[sample_name] = total_mapped
                        all_results.extend(sample_results)
                        
                        for result in sample_results:
                            virus_name = result['Virus']
                            virus_results[virus_name].append(result)
            else:
                # 串行处理
                for sample_info in tqdm(sample_list, desc="处理样本", unit="样本"):
                    sample_name, sample_results, total_mapped = self.process_single_sample(
                        sample_info=sample_info,
                        reference=reference,
                        output_dir=output_dir,
                        coverage_threshold=coverage_threshold,
                        depth_threshold=depth_threshold,
                        threads=threads,
                        skip_fpkm_tpm=skip_fpkm_tpm,
                        keep_intermediate=keep_intermediate
                    )
                    sample_stats[sample_name] = total_mapped
                    all_results.extend(sample_results)
                    
                    for result in sample_results:
                        virus_name = result['Virus']
                        virus_results[virus_name].append(result)
        
        # 生成汇总报告
        self.generate_reports(all_results, virus_results, sample_stats, output_dir, skip_fpkm_tpm)
        
        # 清理临时目录（如果不保留中间文件）
        if not keep_intermediate:
            self.cleanup()
        
        return all_results, virus_results, sample_stats
    
    def create_directory_structure(self, output_dir, keep_intermediate=False):
        """创建输出目录结构"""
        directories = [
            output_dir / "summary",
            output_dir / "logs"
        ]
        
        if keep_intermediate:
            directories.extend([
                output_dir / "alignment",
                output_dir / "coverage",
                output_dir / "intermediate"
            ])
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
            self.logger.debug(f"创建目录: {directory}")
    
    def cleanup_intermediate_files(self):
        """清理中间文件"""
        for file_path in self.intermediate_files:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    self.logger.debug(f"删除中间文件: {file_path}")
            except Exception as e:
                self.logger.warning(f"删除文件失败 {file_path}: {e}")
        
        self.intermediate_files.clear()
    
    def cleanup(self):
        """清理临时文件"""
        import shutil
        try:
            if Path(self.temp_dir).exists():
                shutil.rmtree(self.temp_dir)
                self.logger.debug(f"清理临时目录: {self.temp_dir}")
        except Exception as e:
            self.logger.warning(f"清理临时目录失败: {e}")
    
    def generate_reports(self, all_results, virus_results, sample_stats, output_dir, skip_fpkm_tpm):
        """生成报告文件"""
        self.logger.info("生成报告文件...")
        
        if not all_results:
            self.logger.warning("没有满足条件的结果，不生成报告")
            return
        
        summary_dir = output_dir / "summary"
        
        # 1. 所有病毒的汇总表格
        df_all = pd.DataFrame(all_results)
        all_summary_file = summary_dir / "all_viruses.summary.tsv"
        
        # 按Sample和Virus排序
        df_all = df_all.sort_values(['Sample', 'Virus'])
        df_all.to_csv(all_summary_file, sep='\t', index=False, float_format='%.4f')
        
        # 2. 按病毒分开的表格
        for virus_name, records in tqdm(virus_results.items(), desc="生成病毒报告", unit="病毒"):
            df_virus = pd.DataFrame(records)
            
            # 移除Virus列，因为文件名已包含
            if 'Virus' in df_virus.columns:
                df_virus = df_virus.drop('Virus', axis=1)
            
            # 按Sample排序
            df_virus = df_virus.sort_values('Sample')
            
            # 清理病毒名中的特殊字符，避免文件名问题
            safe_virus_name = re.sub(r'[^\w\-_\. ]', '_', virus_name)
            virus_file = summary_dir / f"{safe_virus_name}.summary.tsv"
            df_virus.to_csv(virus_file, sep='\t', index=False, float_format='%.4f')
        
        # 3. 样本统计摘要
        sample_summary = []
        for sample_name, total_mapped in sample_stats.items():
            sample_viruses = [r for r in all_results if r['Sample'] == sample_name]
            if skip_fpkm_tpm:
                sample_summary.append({
                    'Sample': sample_name,
                    'TotalMappedReads': total_mapped,
                    'DetectedViruses': len(sample_viruses),
                    'HighConfidenceViruses': len([v for v in sample_viruses 
                                                if v['Coverage(%)'] > 90 and v['MeanDepth'] > 10])
                })
            else:
                sample_summary.append({
                    'Sample': sample_name,
                    'TotalMappedReads': total_mapped,
                    'DetectedViruses': len(sample_viruses),
                    'HighConfidenceViruses': len([v for v in sample_viruses 
                                                if v['Coverage(%)'] > 90 and v['MeanDepth'] > 10]),
                    'TotalFPKM': sum(v.get('FPKM', 0) for v in sample_viruses),
                    'TotalTPM': sum(v.get('TPM', 0) for v in sample_viruses)
                })
        
        df_sample_summary = pd.DataFrame(sample_summary)
        sample_summary_file = summary_dir / "sample_summary.tsv"
        df_sample_summary.to_csv(sample_summary_file, sep='\t', index=False, float_format='%.4f')
        
        # 4. 病毒检测频率
        virus_freq = []
        for virus_name, records in virus_results.items():
            if skip_fpkm_tpm:
                virus_freq.append({
                    'Virus': virus_name,
                    'DetectionCount': len(records),
                    'DetectionRate(%)': (len(records) / len(sample_stats)) * 100 if sample_stats else 0,
                    'MeanCoverage(%)': np.mean([r['Coverage(%)'] for r in records]) if records else 0,
                    'MeanDepth': np.mean([r['MeanDepth'] for r in records]) if records else 0
                })
            else:
                virus_freq.append({
                    'Virus': virus_name,
                    'DetectionCount': len(records),
                    'DetectionRate(%)': (len(records) / len(sample_stats)) * 100 if sample_stats else 0,
                    'MeanCoverage(%)': np.mean([r['Coverage(%)'] for r in records]) if records else 0,
                    'MeanDepth': np.mean([r['MeanDepth'] for r in records]) if records else 0,
                    'MeanFPKM': np.mean([r.get('FPKM', 0) for r in records]) if records else 0,
                    'MeanTPM': np.mean([r.get('TPM', 0) for r in records]) if records else 0
                })
        
        df_virus_freq = pd.DataFrame(virus_freq)
        df_virus_freq = df_virus_freq.sort_values('DetectionCount', ascending=False)
        virus_freq_file = summary_dir / "virus_frequency.tsv"
        df_virus_freq.to_csv(virus_freq_file, sep='\t', index=False, float_format='%.4f')
        
        # 5. 生成README文件
        self.generate_readme(output_dir, len(all_results), len(virus_results), 
                           len(sample_stats), skip_fpkm_tpm)
        
        self.logger.info(f"报告生成完成，保存在: {summary_dir}")
    
    def generate_readme(self, output_dir, total_records, total_viruses, total_samples, skip_fpkm_tpm):
        """生成README文件"""
        readme_content = f"""# 病毒分析结果报告

## 分析概况
- 分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- 样本数量: {total_samples}
- 检测到的病毒种类: {total_viruses}
- 总记录数: {total_records}

## 目录结构
{output_dir.name}/
├── summary/ # 汇总文件目录
│ ├── all_viruses.summary.tsv # 所有病毒汇总
│ ├── sample_summary.tsv # 样本统计摘要
│ ├── virus_frequency.tsv # 病毒检测频率
│ └── [病毒名称].summary.tsv # 单个病毒详细统计
├── logs/ # 日志文件目录
│ └── analysis.log # 分析日志
├── alignment/ # BAM文件目录（如保留中间文件）
├── coverage/ # 覆盖度统计文件目录（如保留中间文件）
└── README.txt # 本文件

## 文件说明

### 汇总文件 (summary/)
1. **all_viruses.summary.tsv**: 所有检测到的病毒汇总表格
   - 包含每个样本中每个病毒的详细信息
   - 列: Sample, Virus, Length, CoveredSite, TotalDepth, Coverage(%), MeanDepth{', FPKM, TPM' if not skip_fpkm_tpm else ''}

2. **[病毒名称].summary.tsv**: 单个病毒的详细统计
   - 包含该病毒在所有样本中的检测情况

3. **sample_summary.tsv**: 样本统计摘要
   - 每个样本的总比对reads数、检测到的病毒数等

4. **virus_frequency.tsv**: 病毒检测频率统计
   - 每个病毒的检测频率、平均覆盖度、平均深度等

## 分析流程
1. 序列比对: 使用strobealign
2. 文件排序: 使用samtools sort
3. 覆盖度计算: 使用pandepth
4. 统计分析: 基于阈值过滤和FPKM/TPM计算

## 阈值设置
- Coverage阈值: > 90%
- MeanDepth阈值: > 10

---
生成工具: Virus Analysis Pipeline
版本: 1.0.0
"""

        readme_file = output_dir / "README.txt"
        with open(readme_file, 'w', encoding='utf-8') as f:
            f.write(readme_content)

    def __del__(self):
        """析构函数，确保清理"""
        self.cleanup()

def print_banner():
    """打印程序横幅"""
    banner = """
╔══════════════════════════════════════════════════════════╗
║             病毒分析流程工具 v1.0.0                      ║
║          Virus Analysis Pipeline                         ║
╚══════════════════════════════════════════════════════════╝
"""
    print(banner)

def main():
    print_banner()

    parser = argparse.ArgumentParser(
        description='完整的病毒分析流程：比对 -> 覆盖度计算 -> 统计 -> FPKM/TPM计算',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 基本用法：自动发现测序数据
  python virus_analysis.py --reference virus.fasta --input_dirs /path/to/data

  # 处理多个目录
  python virus_analysis.py --reference virus.fasta --input_dirs /path/to/data1 /path/to/data2

  # 指定输出目录和线程数
  python virus_analysis.py --reference virus.fasta --input_dirs /path/to/data --output-dir results --threads 16 --parallel 4

  # 跳过FPKM/TPM计算
  python virus_analysis.py --reference virus.fasta --input_dirs /path/to/data --skip-fpkm-tpm

  # 保留中间文件
  python virus_analysis.py --reference virus.fasta --input_dirs /path/to/data --keep-intermediate
        """
    )

    # 输入输出参数
    input_group = parser.add_argument_group('输入参数')
    input_group.add_argument('--reference', required=True,
                           help='参考基因组文件 (FASTA格式)')
    input_group.add_argument('--input_dirs', nargs='+', required=True,
                           help='输入目录路径（可指定多个），自动识别测序数据')
    input_group.add_argument('--output-dir', default='./virus_analysis_results',
                           help='输出目录 (默认: ./virus_analysis_results)')

    # 分析参数
    analysis_group = parser.add_argument_group('分析参数')
    analysis_group.add_argument('--coverage', type=float, default=90.0,
                              help='Coverage(%%)阈值 (默认: 90)')
    analysis_group.add_argument('--meandepth', type=float, default=10.0,
                              help='MeanDepth阈值 (默认: 10)')
    analysis_group.add_argument('--threads', type=int, default=8,
                              help='每个样本使用的线程数 (默认: 8)')
    analysis_group.add_argument('--parallel', type=int, default=1,
                              help='并行处理的样本数 (默认: 1)')

    # 其他参数
    other_group = parser.add_argument_group('其他参数')
    other_group.add_argument('--skip-alignment', action='store_true',
                           help='跳过比对步骤（直接分析已有BAM文件）')
    other_group.add_argument('--skip-fpkm-tpm', action='store_true',
                           help='跳过FPKM/TPM计算')
    other_group.add_argument('--keep-intermediate', action='store_true',
                           help='保留中间文件（BAM、统计文件等）')
    other_group.add_argument('--log-level', default='INFO',
                           choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                           help='日志级别 (默认: INFO)')
    other_group.add_argument('--log-file',
                           help='日志文件路径')
    other_group.add_argument('--quiet', action='store_true',
                           help='安静模式，减少屏幕输出')

    args = parser.parse_args()

    # 设置日志
    log_level = getattr(logging, args.log_level)

    # 如果指定了日志文件，同时输出到控制台和文件
    if args.log_file:
        logger = setup_logging(log_level, args.log_file, console_only=args.quiet)
    else:
        # 如果没有指定日志文件，只输出到控制台
        logger = setup_logging(log_level, console_only=args.quiet)

    # 检查依赖
    analyzer = VirusAnalyzer(logger)
    if not analyzer.check_dependencies():
        sys.exit(1)

    # 验证参考基因组文件
    ref_path = Path(args.reference)
    if not ref_path.exists():
        logger.error(f"参考基因组文件不存在: {args.reference}")
        sys.exit(1)

    logger.info(f"参考基因组: {ref_path}")

    # 自动发现测序文件
    logger.info(f"扫描 {len(args.input_dirs)} 个输入目录...")

    samples = analyzer.discover_sequencing_files(args.input_dirs)

    if not samples:
        logger.error("未在输入目录中发现测序文件")
        logger.info("支持的测序文件格式:")
        logger.info("  双端测序: sample_R1.fastq.gz / sample_R2.fastq.gz")
        logger.info("  单端测序: sample.fastq.gz")
        logger.info("其他常见格式: *.fq.gz, *.fastq, *.fq")
        sys.exit(1)

    # 显示样本统计
    if not args.quiet:
        print("\n" + "="*60)
        print("样本统计")
        print("="*60)

        paired_count = sum(1 for _, _, r2 in samples if r2 is not None)
        single_count = len(samples) - paired_count

        print(f"总样本数: {len(samples)}")
        print(f"双端测序: {paired_count}")
        print(f"单端测序: {single_count}")

        # 显示前5个样本
        print("\n前5个样本:")
        for i, (sample_name, r1_path, r2_path) in enumerate(samples[:5]):
            seq_type = "双端" if r2_path else "单端"
            print(f"  {i+1}. {sample_name} ({seq_type}): {Path(r1_path).name}")

        if len(samples) > 5:
            print(f"  ... 还有 {len(samples) - 5} 个样本")
        print("="*60 + "\n")

    try:
        # 执行分析
        start_time = datetime.now()

        all_results, virus_results, sample_stats = analyzer.analyze_samples(
            sample_list=samples,
            reference=args.reference,
            output_dir=args.output_dir,
            coverage_threshold=args.coverage,
            depth_threshold=args.meandepth,
            threads=args.threads,
            parallel=args.parallel,
            skip_alignment=args.skip_alignment,
            skip_fpkm_tpm=args.skip_fpkm_tpm,
            keep_intermediate=args.keep_intermediate
        )

        end_time = datetime.now()
        duration = end_time - start_time

        # 输出总结
        if not args.quiet:
            print("\n" + "="*60)
            print("分析完成!")
            print("="*60)

            print(f"分析时间: {duration}")
            print(f"检测到的病毒种类: {len(virus_results)}")
            print(f"总记录数: {len(all_results)}")

            if virus_results:
                print("\n检测频率最高的病毒:")
                for virus_name, records in sorted(virus_results.items(),
                                                key=lambda x: len(x[1]),
                                                reverse=True)[:5]:
                    detection_rate = (len(records) / len(sample_stats)) * 100
                    print(f"  {virus_name}: {len(records)} 个样本 ({detection_rate:.1f}%)")

            print(f"\n详细结果保存在: {args.output_dir}")
            print("="*60)

        logger.info(f"分析完成，耗时: {duration}")
        logger.info(f"结果目录: {args.output_dir}")

    except KeyboardInterrupt:
        logger.info("用户中断分析")
        sys.exit(130)
    except Exception as e:
        logger.error(f"分析过程中出现错误: {e}", exc_info=True)
        sys.exit(1)

if __name__ == '__main__':
    main()
