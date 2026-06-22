#!/usr/bin/env python3
"""
病毒序列比对、深度统计和定量分析完整流程脚本
使用strobealign进行比对，计算深度和定量指标：MeanDepth, FPKM, TPM
"""

import argparse
import os
import sys
import subprocess
import glob
import gzip
import shutil
import time
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import colorlog
from tqdm import tqdm
import pandas as pd
import numpy as np
import re

# 设置彩色日志
def setup_logging(verbose=False):
    """设置日志格式"""
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        '%(log_color)s%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S',
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'bold_red',
        }
    ))
    
    logger = colorlog.getLogger(__name__)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    return logger

logger = setup_logging()

class Timer:
    """计时器类"""
    def __init__(self, name=""):
        self.name = name
        self.start_time = None
        self.end_time = None
    
    def __enter__(self):
        self.start_time = time.time()
        logger.info(f"⏱️  开始: {self.name}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        duration = self.end_time - self.start_time
        if exc_type is None:
            logger.info(f"✅ 完成: {self.name} [耗时: {self.format_duration(duration)}]")
        else:
            logger.error(f"❌ 失败: {self.name} [耗时: {self.format_duration(duration)}]")
    
    @staticmethod
    def format_duration(seconds):
        """格式化时间"""
        if seconds < 60:
            return f"{seconds:.1f}秒"
        elif seconds < 3600:
            minutes = seconds / 60
            return f"{minutes:.1f}分钟"
        else:
            hours = seconds / 3600
            return f"{hours:.1f}小时"

class VirusAnalysisPipeline:
    """病毒分析流程类"""
    
    def __init__(self, args):
        self.args = args
        self.samples = []
        self.results = []
        self.virus_results = defaultdict(list)
        self.output_dir = Path(args.output_dir)
        
        # 检查必要工具
        self.check_tools()
    
    def check_tools(self):
        """检查必要的工具是否可用"""
        required_tools = {
            'samtools': '测序数据处理',
            'pandepth': '深度计算',
            'strobealign': '序列比对'
        }
        
        missing_tools = []
        for tool, description in required_tools.items():
            if shutil.which(tool) is None:
                missing_tools.append(f"{tool} ({description})")
        
        if missing_tools:
            logger.error("❌ 缺少必要的工具:")
            for tool in missing_tools:
                logger.error(f"  - {tool}")
            
            # 提供安装建议
            logger.info("\n💡 安装建议:")
            if 'strobealign' in missing_tools:
                logger.info("  安装strobealign:")
                logger.info("    conda install -c bioconda strobealign")
                logger.info("    或从 https://github.com/ksahlin/strobealign 编译安装")
            if 'pandepth' in missing_tools:
                logger.info("  安装pandepth:")
                logger.info("    conda install -c bioconda pandepth")
            if 'samtools' in missing_tools:
                logger.info("  安装samtools:")
                logger.info("    conda install -c bioconda samtools")
            
            sys.exit(1)
        
        logger.info("✅ 所有必要工具可用")
    
    def setup_output_directory(self):
        """设置输出目录结构"""
        subdirs = ['bam', 'stat', 'summary', 'logs']
        
        with Timer("创建输出目录结构"):
            for subdir in subdirs:
                (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)
                logger.debug(f"创建目录: {self.output_dir / subdir}")
        
        logger.info(f"📁 输出目录: {self.output_dir}")
        return self.output_dir
    
    def find_samples(self):
        """查找样本文件"""
        input_dir = Path(self.args.input_dir)
        # 支持的fastq格式
        fastq_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
        
        samples = []
        
        if self.args.single_end:
            logger.info("🔍 搜索单端测序样本...")
            for ext in fastq_extensions:
                fastq_files = list(input_dir.glob(f'*{ext}'))
                for fastq in tqdm(fastq_files, desc="扫描文件", unit="文件"):
                    sample_name = fastq.name.replace(ext, '').replace('.gz', '')
                    # 清理样本名，去除常见后缀
                    for suffix in ['_1', '_R1', '.R1', "R1_10239",'.1', '_unmapped', '.unmapped', '_trimmed', '.trimmed']:
                        if suffix in sample_name:
                            sample_name = sample_name.split(suffix)[0]
                    samples.append({
                        'name': sample_name,
                        'r1': str(fastq),
                        'r2': None
                    })
        else:
            logger.info("🔍 搜索双端测序样本...")
            # 支持多种命名模式
            patterns = [
                ('*_1.f*q*', '*_2.f*q*'),  # sample_1.fq.gz, sample_2.fq.gz
                ('*_R1*.f*q*', '*_R2*.f*q*'),  # sample_R1.fq.gz, sample_R2.fq.gz
                ('*.R1.*', '*.R2.*'),  # sample.R1.fq.gz, sample.R2.fq.gz
                ('*.1.f*q*', '*.2.f*q*'),  # sample.1.fq.gz, sample.2.fq.gz
                ('*_1_*.f*q*', '*_2_*.f*q*'),  # sample_1_xxx.fq.gz
                ('*_unmapped.R1.fq.gz', '*_unmapped.R2.fq.gz'),  # CRR527054.unmapped.R1.fq.gz
                ('*unmapped.R1.fq.gz', '*unmapped.R2.fq.gz'),  # CRR527054.unmapped.R1.fq.gz
                ('*_1_unmapped.*', '*_2_unmapped.*'),  # sample_1_unmapped.fq.gz
                ('*.unmapped.R1_10239.*', '*.unmapped.R2_10239.*'),
            ]
            
            found_files = set()
            
            for r1_pattern, r2_pattern in tqdm(patterns, desc="匹配文件模式", unit="模式"):
                r1_files = list(input_dir.glob(r1_pattern))
                
                for r1_file in r1_files:
                    if r1_file in found_files:
                        continue
                    
                    # 根据R1文件名生成可能的R2文件名
                    r1_name = r1_file.name
                    r2_name = None
                    
                    # 根据模式替换
                    if '_1.' in r1_name:
                        r2_name = r1_name.replace('_1.', '_2.')
                    elif '_R1' in r1_name:
                        r2_name = r1_name.replace('_R1', '_R2')
                    elif '.R1.' in r1_name:
                        r2_name = r1_name.replace('.R1.', '.R2.')
                    elif '.1.' in r1_name:
                        r2_name = r1_name.replace('.1.', '.2.')
                    elif '_1_' in r1_name:
                        r2_name = r1_name.replace('_1_', '_2_')
                    elif '.unmapped.R1.' in r1_name:
                        r2_name = r1_name.replace('.unmapped.R1.', '.unmapped.R2.')
                    elif 'unmapped.R1.' in r1_name:
                        r2_name = r1_name.replace('unmapped.R1.', 'unmapped.R2.')
                    elif '_1_unmapped.' in r1_name:
                        r2_name = r1_name.replace('_1_unmapped.', '_2_unmapped.')
                    elif '.unmapped.R1_10239.' in r1_name:
                        r2_name = r1_name.replace('.unmapped.R1_10239.', '.unmapped.R2_10239.')
                    
                    if r2_name:
                        r2_file = input_dir / r2_name
                        if r2_file.exists():
                            # 提取样本名
                            sample_name = r1_name
                            # 去除常见后缀
                            for suffix in ['_1', '_R1', '.R1', '.1', '_unmapped', '.unmapped', '_1_unmapped']:
                                if suffix in sample_name:
                                    sample_name = sample_name.split(suffix)[0]
                            
                            # 如果还有扩展名，去掉扩展名
                            sample_name = sample_name.split('.')[0]
                            
                            samples.append({
                                'name': sample_name,
                                'r1': str(r1_file),
                                'r2': str(r2_file)
                            })
                            found_files.add(r1_file)
                            found_files.add(r2_file)
            
            # 如果没有找到，尝试直接搜索所有fastq文件并配对
            if not samples:
                logger.info("尝试直接搜索并配对fastq文件...")
                all_fastq_files = []
                for ext in fastq_extensions:
                    all_fastq_files.extend(input_dir.glob(f'*{ext}'))
                
                # 按前缀分组
                file_groups = defaultdict(list)
                for f in all_fastq_files:
                    name = f.name
                    # 尝试提取前缀
                    prefix = re.sub(r'[._](R?[12]|unmapped\.R?[12])[._].*', '', name)
                    if '.' in prefix:
                        prefix = prefix.split('.')[0]
                    file_groups[prefix].append(str(f))
                
                # 为每个前缀创建样本
                for prefix, files in file_groups.items():
                    if len(files) == 2:
                        # 确定R1和R2
                        r1, r2 = None, None
                        for f in files:
                            if any(pattern in f for pattern in ['_1', '_R1', '.R1.', '.1.', 'unmapped.R1']):
                                r1 = f
                            elif any(pattern in f for pattern in ['_2', '_R2', '.R2.', '.2.', 'unmapped.R2']):
                                r2 = f
                        
                        if r1 and r2:
                            samples.append({
                                'name': prefix,
                                'r1': r1,
                                'r2': r2
                            })
        
        # 如果提供了样本列表文件
        if self.args.sample_list:
            logger.info(f"📋 从文件加载样本列表: {self.args.sample_list}")
            samples = []
            with open(self.args.sample_list, 'r') as f:
                for i, line in enumerate(f, 1):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            sample_name = parts[0]
                            r1_file = parts[1]
                            r2_file = parts[2] if len(parts) > 2 else None
                            
                            # 检查文件是否存在
                            if not Path(r1_file).exists():
                                logger.error(f"文件不存在: {r1_file}")
                                continue
                            if r2_file and not Path(r2_file).exists():
                                logger.error(f"文件不存在: {r2_file}")
                                continue
                            
                            samples.append({
                                'name': sample_name,
                                'r1': r1_file,
                                'r2': r2_file
                            })
            
            logger.info(f"从列表文件加载了 {len(samples)} 个样本")
        
        self.samples = samples
        
        if not samples:
            logger.error("❌ 未找到任何样本文件")
            logger.info("请检查以下可能的问题:")
            logger.info("1. 确保输入目录包含fastq文件")
            logger.info("2. 文件命名应该符合常见模式:")
            logger.info("   - sample_1.fq.gz 和 sample_2.fq.gz")
            logger.info("   - sample_R1.fastq.gz 和 sample_R2.fastq.gz")
            logger.info("   - sample.R1.fq.gz 和 sample.R2.fq.gz")
            logger.info("   - sample.unmapped.R1.fq.gz 和 sample.unmapped.R2.fq.gz")
            logger.info("3. 如果是单端数据，请使用 --single_end 参数")
            logger.info("4. 可以使用 --sample_list 参数指定样本列表文件")
            sys.exit(1)
        
        # 去重样本
        unique_samples = []
        seen_names = set()
        for sample in samples:
            if sample['name'] not in seen_names:
                unique_samples.append(sample)
                seen_names.add(sample['name'])
            else:
                logger.warning(f"重复样本名: {sample['name']}，跳过")
        
        self.samples = unique_samples
        
        # 记录样本信息
        sample_info_file = self.output_dir / 'logs' / 'samples.txt'
        with open(sample_info_file, 'w') as f:
            f.write("样本名\tR1文件\tR2文件\t文件大小(MB)\n")
            for sample in self.samples:
                r1_size = os.path.getsize(sample['r1']) / (1024 * 1024) if os.path.exists(sample['r1']) else 0
                r2_size = os.path.getsize(sample['r2']) / (1024 * 1024) if sample['r2'] and os.path.exists(sample['r2']) else 0
                f.write(f"{sample['name']}\t{sample['r1']}\t{sample['r2'] or 'N/A'}\t{r1_size:.1f}/{r2_size:.1f}\n")
        
        logger.info(f"✅ 找到 {len(self.samples)} 个样本")
        
        # 显示样本统计信息
        total_size = 0
        for i, sample in enumerate(self.samples[:5], 1):
            r1_size = os.path.getsize(sample['r1']) / (1024 * 1024) if os.path.exists(sample['r1']) else 0
            total_size += r1_size
            if sample['r2'] and os.path.exists(sample['r2']):
                r2_size = os.path.getsize(sample['r2']) / (1024 * 1024)
                total_size += r2_size
                logger.debug(f"  {i}. {sample['name']}: R1={r1_size:.1f}MB, R2={r2_size:.1f}MB")
            else:
                logger.debug(f"  {i}. {sample['name']}: R1={r1_size:.1f}MB")
        
        if len(self.samples) > 5:
            logger.debug(f"  ... 还有 {len(self.samples)-5} 个样本")
        
        logger.info(f"📊 总数据量: {total_size:.1f} MB")
    
    def align_sample(self, sample):
        """使用strobealign比对单个样本"""
        try:
            sample_name = sample['name']
            
            # 输出文件路径
            bam_dir = self.output_dir / 'bam'
            sorted_bam = bam_dir / f'{sample_name}.sorted.bam'
            
            # 如果已存在且设置了跳过，则跳过
            if sorted_bam.exists() and self.args.skip_alignment:
                logger.debug(f"⏭️  跳过已存在的比对文件: {sorted_bam}")
                return sorted_bam
            
            # 比对命令
            logger.info(f"🔄 开始比对样本: {sample_name}")
            
            with Timer(f"比对 {sample_name}"):
                # 构建strobealign命令
                align_cmd = [
                    'strobealign',
                    '-t', str(self.args.align_threads),
                    self.args.reference,
                    sample['r1'],
                ]
                
                if sample['r2']:
                    align_cmd.append(sample['r2'])
                
                # 排序命令
                sort_cmd = [
                    'samtools', 'sort',
                    '-@', str(min(4, self.args.threads)),
                    '-o', str(sorted_bam)
                ]
                
                logger.debug(f"比对命令: {' '.join(align_cmd)}")
                logger.debug(f"排序命令: {' '.join(sort_cmd)}")
                
                # 执行管道命令
                align_process = subprocess.Popen(
                    align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                )
                
                sort_process = subprocess.Popen(
                    sort_cmd, stdin=align_process.stdout, stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE, text=True
                )
                
                align_process.stdout.close()
                
                # 等待两个进程完成
                _, align_stderr = align_process.communicate()
                sort_stdout, sort_stderr = sort_process.communicate()
                
                # 检查错误
                if align_process.returncode != 0:
                    logger.error(f"strobealign比对失败: {align_stderr[:500]}")
                    return None
                
                if sort_process.returncode != 0:
                    logger.error(f"排序失败: {sort_stderr[:500]}")
                    return None
                
                # 记录比对统计信息
                if align_stderr:
                    for line in align_stderr.split('\n'):
                        if 'aligned' in line.lower() or 'reads' in line.lower():
                            logger.debug(f"比对统计: {line.strip()}")
            
            # 创建索引
            with Timer(f"索引 {sample_name}"):
                index_cmd = ['samtools', 'index', str(sorted_bam)]
                result = subprocess.run(index_cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.error(f"索引创建失败: {result.stderr}")
                    return None
            
            # 获取比对统计
            flagstat_cmd = ['samtools', 'flagstat', str(sorted_bam)]
            flagstat_result = subprocess.run(flagstat_cmd, capture_output=True, text=True)
            if flagstat_result.returncode == 0:
                total_reads = 0
                mapped_reads = 0
                for line in flagstat_result.stdout.split('\n'):
                    if 'in total' in line:
                        total_reads = int(line.split()[0])
                    elif 'mapped (' in line:
                        mapped_reads = int(line.split()[0])
                
                if total_reads > 0:
                    mapping_rate = (mapped_reads / total_reads) * 100
                    logger.info(f"📊 {sample_name}: 总reads={total_reads:,}, 比对reads={mapped_reads:,}, 比对率={mapping_rate:.1f}%")
            
            logger.info(f"✅ 样本 {sample_name} 比对完成")
            return sorted_bam
        
        except Exception as e:
            logger.error(f"❌ 样本 {sample['name']} 比对失败: {str(e)}")
            return None
    
    def run_pandepth(self, bam_file):
        """运行pandepth计算覆盖度"""
        try:
            bam_path = Path(bam_file)
            sample_name = bam_path.stem.replace('.sorted', '')
            output_prefix = self.output_dir / 'stat' / sample_name
            
            # 如果已存在且设置了跳过，则跳过
            stat_file = output_prefix.with_suffix('.chr.stat.gz')
            if stat_file.exists() and self.args.skip_depth:
                logger.debug(f"⏭️  跳过已存在的覆盖度文件: {stat_file}")
                return stat_file
            
            logger.info(f"📊 计算样本 {sample_name} 的覆盖度...")
            
            with Timer(f"深度计算 {sample_name}"):
                pandepth_cmd = [
                    'pandepth',
                    '-a',  # 输出所有contig
                    '-i', str(bam_file),
                    '-o', str(output_prefix),
                    '-t', str(min(10, self.args.threads))
                ]
                
                logger.debug(f"深度计算命令: {' '.join(pandepth_cmd)}")
                result = subprocess.run(pandepth_cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    logger.error(f"pandepth运行失败: {result.stderr[:500]}")
                    return None
                
                logger.debug(f"深度计算输出:\n{result.stdout[:500]}...")
            
            logger.info(f"✅ 样本 {sample_name} 覆盖度计算完成")
            return stat_file
        
        except Exception as e:
            logger.error(f"❌ pandepth运行出错: {str(e)}")
            return None
    
    def parse_stat_file(self, stat_file, sample_name):
        """解析统计文件并计算FPKM/TPM"""
        try:
            results = []
            virus_stats = []
            
            with gzip.open(stat_file, 'rt') as f:
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
                        length_val = int(length)
                        covered_site_val = int(covered_site)
                        total_depth_val = int(total_depth)
                        coverage_val = float(coverage)
                        mean_depth_val = float(mean_depth)
                        
                        # 存储所有病毒数据
                        virus_stats.append({
                            'Virus': virus_name,
                            'Length': length_val,
                            'TotalDepth': total_depth_val,
                            'Coverage(%)': coverage_val,
                            'MeanDepth': mean_depth_val
                        })
                        
                        # 检查是否满足条件
                        if coverage_val > self.args.coverage and mean_depth_val > self.args.meandepth:
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
            
            # 计算FPKM和TPM
            if virus_stats:
                total_reads = sum([v['TotalDepth'] for v in virus_stats])
                
                # 计算每千碱基读数
                rpkb_values = []
                for v in virus_stats:
                    rpkb = (v['TotalDepth'] * 1000) / v['Length'] if v['Length'] > 0 else 0
                    rpkb_values.append(rpkb)
                
                total_rpkb = sum(rpkb_values)
                
                # 为结果添加FPKM和TPM
                for result in results:
                    virus_name = result['Virus']
                    v_stats = next((v for v in virus_stats if v['Virus'] == virus_name), None)
                    
                    if v_stats and total_reads > 0:
                        fpkm = (v_stats['TotalDepth'] * 1e9) / (total_reads * v_stats['Length']) if v_stats['Length'] > 0 else 0
                        rpkb = (v_stats['TotalDepth'] * 1000) / v_stats['Length'] if v_stats['Length'] > 0 else 0
                        tpm = (rpkb * 1e6) / total_rpkb if total_rpkb > 0 else 0
                        
                        result['FPKM'] = round(fpkm, 2)
                        result['TPM'] = round(tpm, 2)
                    else:
                        result['FPKM'] = 0
                        result['TPM'] = 0
            
            # 应用FPKM和TPM过滤
            filtered_results = [
                r for r in results 
                if r['FPKM'] >= self.args.min_fpkm and r['TPM'] >= self.args.min_tpm
            ]
            
            logger.debug(f"样本 {sample_name}: 找到 {len(filtered_results)} 条满足条件的记录")
            return filtered_results
        
        except Exception as e:
            logger.error(f"❌ 解析统计文件 {stat_file} 失败: {str(e)}")
            return []
    
    def process_sample(self, sample):
        """处理单个样本的完整流程"""
        try:
            logger.info(f"\n📋 处理样本: {sample['name']}")
            
            # 1. 比对
            bam_file = self.align_sample(sample)
            if not bam_file:
                return None
            
            # 2. 计算覆盖度
            stat_file = self.run_pandepth(bam_file)
            if not stat_file:
                return None
            
            # 3. 解析统计文件
            sample_results = self.parse_stat_file(stat_file, sample['name'])
            
            logger.info(f"✅ 样本 {sample['name']} 处理完成")
            return sample_results
        
        except Exception as e:
            logger.error(f"❌ 处理样本 {sample['name']} 失败: {str(e)}")
            return None
    
    def process_samples_parallel(self, max_workers):
        """并行处理样本"""
        all_results = []
        
        logger.info(f"🔀 并行处理样本，最大并行任务数: {max_workers}")
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
            future_to_sample = {
                executor.submit(self.process_sample, sample): sample 
                for sample in self.samples
            }
            
            # 使用tqdm显示进度
            with tqdm(total=len(self.samples), desc="处理样本", unit="样本") as pbar:
                for future in as_completed(future_to_sample):
                    sample = future_to_sample[future]
                    try:
                        result = future.result()
                        if result:
                            all_results.extend(result)
                            pbar.set_postfix_str(f"✅ {sample['name']}")
                        else:
                            pbar.set_postfix_str(f"❌ {sample['name']}失败")
                    except Exception as e:
                        logger.error(f"❌ 样本 {sample['name']} 处理异常: {str(e)}")
                        pbar.set_postfix_str(f"❌ {sample['name']}异常")
                    
                    pbar.update(1)
        
        return all_results
    
    def process_samples_serial(self):
        """串行处理样本"""
        all_results = []
        
        with tqdm(self.samples, desc="处理样本", unit="样本") as pbar:
            for sample in pbar:
                pbar.set_postfix_str(f"正在处理 {sample['name']}")
                result = self.process_sample(sample)
                if result:
                    all_results.extend(result)
                    pbar.set_postfix_str(f"✅ {sample['name']}")
                else:
                    pbar.set_postfix_str(f"❌ {sample['name']}失败")
        
        return all_results
    
    def run_pipeline(self):
        """运行完整分析流程"""
        start_time = datetime.now()
        logger.info("=" * 60)
        logger.info("🚀 开始病毒分析流程 (使用strobealign)")
        logger.info("=" * 60)
        
        # 显示配置信息
        logger.info("📋 分析配置:")
        logger.info(f"  - 比对工具: strobealign")
        logger.info(f"  - 参考基因组: {self.args.reference}")
        logger.info(f"  - 输入目录: {self.args.input_dir}")
        logger.info(f"  - 输出目录: {self.output_dir}")
        logger.info(f"  - 总线程数: {self.args.threads}")
        logger.info(f"  - 比对线程数: {self.args.align_threads}")
        logger.info(f"  - 过滤条件: Coverage > {self.args.coverage}%, MeanDepth > {self.args.meandepth}")
        logger.info(f"  - 定量过滤: FPKM >= {self.args.min_fpkm}, TPM >= {self.args.min_tpm}")
        
        # 设置输出目录
        self.setup_output_directory()
        
        # 查找样本
        self.find_samples()
        
        # 处理样本
        logger.info(f"\n🔄 开始处理 {len(self.samples)} 个样本")
        
        if self.args.parallel and self.args.threads > 1:
            # 计算最大并行任务数
            # 每个样本的比对会使用 align_threads 个线程
            # 我们可以同时运行的样本数 = 总线程数 / 比对线程数
            max_parallel_jobs = max(1, self.args.threads // self.args.align_threads)
            
            # 但不能超过样本数
            max_parallel_jobs = min(max_parallel_jobs, len(self.samples))
            
            # 如果用户指定了并行任务数，使用用户指定的值
            if self.args.parallel_jobs:
                max_parallel_jobs = min(self.args.parallel_jobs, len(self.samples))
                logger.info(f"用户指定并行任务数: {self.args.parallel_jobs}")
            
            logger.info(f"计算得到的并行任务数: {max_parallel_jobs}")
            all_results = self.process_samples_parallel(max_parallel_jobs)
        else:
            logger.info("🔂 使用串行处理")
            all_results = self.process_samples_serial()
        
        # 汇总结果
        if all_results:
            self.save_results(all_results)
        else:
            logger.warning("⚠️  未找到满足条件的记录")
        
        # 计算总耗时
        end_time = datetime.now()
        total_duration = end_time - start_time
        
        logger.info("=" * 60)
        logger.info("✨ 分析流程完成")
        logger.info(f"📊 总计: {len(all_results)} 条记录")
        logger.info(f"⏱️  总耗时: {self.format_duration(total_duration)}")
        
        # 显示结果文件位置
        if all_results:
            summary_file = self.output_dir / 'summary' / f'all_viruses.summary.{self.args.format}'
            logger.info(f"💾 结果文件: {summary_file}")
        
        logger.info("=" * 60)
    
    def save_results(self, all_results):
        """保存结果到文件"""
        logger.info(f"\n💾 保存结果...")
        
        with Timer("汇总和保存结果"):
            # 转换为DataFrame
            df_all = pd.DataFrame(all_results)
            
            # 保存汇总文件
            summary_dir = self.output_dir / 'summary'
            summary_file = summary_dir / f'all_viruses.summary.{self.args.format}'
            
            if self.args.format == 'csv':
                df_all.to_csv(summary_file, index=False)
            else:  # tsv
                df_all.to_csv(summary_file, sep='\t', index=False)
            
            logger.info(f"✅ 汇总表格已保存: {summary_file}")
            
            # 按病毒分组保存
            if 'Virus' in df_all.columns:
                viruses = df_all['Virus'].unique()
                logger.info(f"📁 按病毒分组保存 ({len(viruses)} 个病毒)...")
                
                for virus in tqdm(viruses, desc="保存病毒文件", unit="病毒"):
                    df_virus = df_all[df_all['Virus'] == virus].copy()
                    if 'Virus' in df_virus.columns:
                        df_virus = df_virus.drop('Virus', axis=1)
                    
                    virus_file = summary_dir / f'{virus}.summary.{self.args.format}'
                    
                    if self.args.format == 'csv':
                        df_virus.to_csv(virus_file, index=False)
                    else:
                        df_virus.to_csv(virus_file, sep='\t', index=False)
                    
                    self.virus_results[virus].append(df_virus)
            
            # 生成统计报告
            self.generate_report(df_all)
    
    def generate_report(self, df_all):
        """生成分析报告"""
        report_file = self.output_dir / 'summary' / 'analysis_report.txt'
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("病毒分析流程报告 (使用strobealign比对)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"参考基因组: {self.args.reference}\n")
            f.write(f"总线程数: {self.args.threads}\n")
            f.write(f"比对线程数: {self.args.align_threads}\n")
            f.write(f"过滤条件: Coverage > {self.args.coverage}%, MeanDepth > {self.args.meandepth}\n")
            f.write(f"定量过滤: FPKM >= {self.args.min_fpkm}, TPM >= {self.args.min_tpm}\n\n")
            
            f.write(f"总样本数: {len(self.samples)}\n")
            f.write(f"总记录数: {len(df_all)}\n\n")
            
            # 按病毒统计
            if 'Virus' in df_all.columns:
                virus_counts = df_all['Virus'].value_counts()
                f.write("按病毒统计:\n")
                f.write("-" * 40 + "\n")
                for virus, count in virus_counts.items():
                    virus_data = df_all[df_all['Virus'] == virus]
                    mean_depth = virus_data['MeanDepth'].mean()
                    mean_coverage = virus_data['Coverage(%)'].mean()
                    mean_fpkm = virus_data['FPKM'].mean() if 'FPKM' in virus_data.columns else 0
                    mean_tpm = virus_data['TPM'].mean() if 'TPM' in virus_data.columns else 0
                    
                    f.write(f"{virus}:\n")
                    f.write(f"  样本数: {count}\n")
                    f.write(f"  平均深度: {mean_depth:.2f}\n")
                    f.write(f"  平均覆盖度: {mean_coverage:.2f}%\n")
                    if mean_fpkm > 0:
                        f.write(f"  平均FPKM: {mean_fpkm:.2f}\n")
                    if mean_tpm > 0:
                        f.write(f"  平均TPM: {mean_tpm:.2f}\n")
                    f.write("\n")
            
            # 按样本统计
            if 'Sample' in df_all.columns:
                sample_counts = df_all['Sample'].value_counts()
                f.write("按样本统计:\n")
                f.write("-" * 40 + "\n")
                for sample, count in sample_counts.items():
                    f.write(f"{sample}: {count} 个病毒\n")
        
        logger.info(f"📄 分析报告已保存: {report_file}")
    
    @staticmethod
    def format_duration(td):
        """格式化时间差"""
        total_seconds = int(td.total_seconds())
        hours, remainder = divmod(total_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        parts = []
        if hours > 0:
            parts.append(f"{hours}小时")
        if minutes > 0:
            parts.append(f"{minutes}分钟")
        if seconds > 0 or not parts:
            parts.append(f"{seconds}秒")
        
        return "".join(parts)

def main():
    parser = argparse.ArgumentParser(
        description='病毒序列比对、深度统计和定量分析完整流程 (使用strobealign)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python virus_analysis_pipeline.py --input_dir fastq --reference virus.fasta
  
  # 使用并行处理，指定并行任务数
  python virus_analysis_pipeline.py --input_dir fastq --reference virus.fasta --parallel --parallel_jobs 4
  
  # 自定义过滤条件
  python virus_analysis_pipeline.py --input_dir fastq --reference virus.fasta --coverage 95 --meandepth 20
  
  # 处理单端数据
  python virus_analysis_pipeline.py --input_dir fastq --reference virus.fasta --single_end
  
  # 使用样本列表文件
  python virus_analysis_pipeline.py --reference virus.fasta --sample_list samples.txt
  
  # 使用已有BAM文件，跳过比对步骤
  python virus_analysis_pipeline.py --reference virus.fasta --input_dir bam_dir --skip_alignment
  
  # 查看详细日志
  python virus_analysis_pipeline.py --input_dir fastq --reference virus.fasta --verbose
        """
    )
    
    # 输入输出参数
    input_group = parser.add_argument_group('输入输出参数')
    input_group.add_argument(
        '--input_dir',
        type=str,
        help='输入文件夹路径（包含fastq文件）'
    )
    input_group.add_argument(
        '--reference',
        type=str,
        required=True,
        help='参考基因组fasta文件路径'
    )
    input_group.add_argument(
        '--output_dir',
        type=str,
        default='./virus_analysis',
        help='输出目录路径 [默认: ./virus_analysis]'
    )
    input_group.add_argument(
        '--sample_list',
        type=str,
        help='样本列表文件（每行: 样本名 fastq1 fastq2）'
    )
    input_group.add_argument(
        '--single_end',
        action='store_true',
        help='使用单端测序数据'
    )
    
    # 比对参数
    align_group = parser.add_argument_group('比对参数')
    align_group.add_argument(
        '--align_threads',
        type=int,
        default=8,
        help='每个样本比对的线程数 [默认: 8]'
    )
    
    # 过滤参数
    filter_group = parser.add_argument_group('过滤参数')
    filter_group.add_argument(
        '--coverage',
        type=float,
        default=90.0,
        help='Coverage(%%)阈值 [默认: 90]'
    )
    filter_group.add_argument(
        '--meandepth',
        type=float,
        default=10.0,
        help='MeanDepth阈值 [默认: 10]'
    )
    filter_group.add_argument(
        '--min_fpkm',
        type=float,
        default=0.0,
        help='FPKM最小值过滤 [默认: 0]'
    )
    filter_group.add_argument(
        '--min_tpm',
        type=float,
        default=0.0,
        help='TPM最小值过滤 [默认: 0]'
    )
    
    # 并行处理参数
    parallel_group = parser.add_argument_group('并行处理参数')
    parallel_group.add_argument(
        '--parallel',
        action='store_true',
        help='并行处理样本'
    )
    parallel_group.add_argument(
        '--parallel_jobs',
        type=int,
        help='并行任务数（如果不指定，自动计算）'
    )
    parallel_group.add_argument(
        '--threads',
        type=int,
        default=4,
        help='总线程数 [默认: 4]'
    )
    
    # 流程控制参数
    control_group = parser.add_argument_group('流程控制参数')
    control_group.add_argument(
        '--skip_alignment',
        action='store_true',
        help='跳过比对步骤（使用已有BAM文件）'
    )
    control_group.add_argument(
        '--skip_depth',
        action='store_true',
        help='跳过覆盖度计算步骤（使用已有统计文件）'
    )
    control_group.add_argument(
        '--format',
        type=str,
        choices=['csv', 'tsv'],
        default='tsv',
        help='输出文件格式 [默认: tsv]'
    )
    
    # 其他参数
    other_group = parser.add_argument_group('其他参数')
    other_group.add_argument(
        '--verbose',
        action='store_true',
        help='输出详细日志信息'
    )
    other_group.add_argument(
        '--version',
        action='version',
        version='病毒分析流程 v1.0 (strobealign版)'
    )
    
    args = parser.parse_args()
    
    # 验证参数
    if not args.input_dir and not args.sample_list:
        parser.error("必须提供 --input_dir 或 --sample_list 参数")
    
    # 检查参考基因组文件
    if not os.path.exists(args.reference):
        parser.error(f"参考基因组文件不存在: {args.reference}")
    
    # 检查输入目录
    if args.input_dir and not os.path.exists(args.input_dir):
        parser.error(f"输入目录不存在: {args.input_dir}")
    
    # 设置日志级别
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # 运行分析流程
    pipeline = VirusAnalysisPipeline(args)
    pipeline.run_pipeline()

if __name__ == '__main__':
    main()
