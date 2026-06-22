#!/usr/bin/env python3
"""
病毒序列比引、深度统计和绝对定量分析完整流程脚本
支持多种比对工具，引入三大定量指标：FPKM, RPM (评估绝对载量) 和 TPM (评估相对群落组成)
支持结合 RVDB 分类文件 (RVDB_Taxon_Current.tab) 补充分类注释信息
"""

import argparse
import os
import sys
import subprocess
import glob
import gzip
import shutil
import time
import io
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
        self.index_path = None
        self.taxon_db = None
        
        # 检查必要工具
        self.check_tools()
        # 尝试加载分类信息
        self.load_rvdb_taxonomy()
    
    def check_tools(self):
        """检查必要的工具是否可用"""
        required_tools = {
            'samtools': '测序数据处理',
            'pandepth': '深度计算',
            self.args.tool: '序列比对'
        }
        
        if self.args.tool == 'bowtie2':
            required_tools['bowtie2-build'] = 'Bowtie2索引构建'
        elif self.args.tool == 'hisat2':
            required_tools['hisat2-build'] = 'HISAT2索引构建'
            
        missing_tools = []
        for tool, description in required_tools.items():
            if shutil.which(tool) is None:
                missing_tools.append(f"{tool} ({description})")
        
        if missing_tools:
            logger.error("❌ 缺少必要的工具:")
            for tool in missing_tools:
                logger.error(f"  - {tool}")
            
            logger.info("\n💡 安装建议 (通过 conda):")
            for tool_name in required_tools.keys():
                if any(tool_name in missing for missing in missing_tools):
                    logger.info(f"    conda install -c bioconda {tool_name}")
            sys.exit(1)
        
        logger.info("✅ 所有必要工具可用")

    def load_rvdb_taxonomy(self):
        """加载 RVDB 分类注释文件"""
        if not self.args.rvdb_taxon:
            return
            
        logger.info(f"📂 正在解析 RVDB 分类注释文件: {self.args.rvdb_taxon}")
        if not os.path.exists(self.args.rvdb_taxon):
            logger.error(f"❌ 找不到分类文件: {self.args.rvdb_taxon}")
            sys.exit(1)
            
        try:
            # 安全读取文件，跳过以 #### 开头的注释行
            valid_lines = []
            with open(self.args.rvdb_taxon, 'r', encoding='utf-8') as f:
                for line in f:
                    if not line.startswith('####'):
                        valid_lines.append(line)
            
            df_tax = pd.read_csv(io.StringIO(''.join(valid_lines)), sep='\t', dtype=str)
            
            req_cols = ['accession', 'description', 'taxonomy', 'taxonomy_id']
            if not all(c in df_tax.columns for c in req_cols):
                logger.warning("⚠️ RVDB 分类文件缺少所需列，注释功能将被跳过！")
                return

            taxon_dict = {}
            for _, row in df_tax.iterrows():
                acc = str(row['accession']).strip()
                info = {
                    'description': str(row['description']) if pd.notna(row['description']) else '',
                    'taxonomy': str(row['taxonomy']) if pd.notna(row['taxonomy']) else '',
                    'taxonomy_id': str(row['taxonomy_id']) if pd.notna(row['taxonomy_id']) else ''
                }
                taxon_dict[acc] = info
                
                # 创建一份去掉 .1 版本的冗余映射，增强容错能力
                acc_base = acc.split('.')[0]
                if acc_base not in taxon_dict:
                    taxon_dict[acc_base] = info
                    
            self.taxon_db = taxon_dict
            logger.info(f"✅ 成功加载 {len(self.taxon_db)} 条病毒分类及描述信息")
            
        except Exception as e:
            logger.error(f"❌ 加载 RVDB 分类文件失败: {e}")
            self.taxon_db = None
    
    def setup_output_directory(self):
        """设置输出目录结构"""
        subdirs = ['bam', 'stat', 'summary', 'logs', 'index']
        with Timer("创建输出目录结构"):
            for subdir in subdirs:
                (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)
                logger.debug(f"创建目录: {self.output_dir / subdir}")
        logger.info(f"📁 输出目录: {self.output_dir}")
        return self.output_dir
        
    def build_index(self):
        """构建参考基因组索引"""
        tool = self.args.tool
        ref_path = Path(self.args.reference).resolve()
        
        if tool in ['strobealign', 'minimap2']:
            self.index_path = str(ref_path)
            return
            
        index_dir = self.output_dir / 'index'
        prefix = index_dir / f"{ref_path.stem}_{tool}"
        self.index_path = str(prefix)
        
        if tool == 'bwa':
            if Path(f"{prefix}.bwt").exists(): return
            cmd = ['bwa', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bwa-mem2':
            if Path(f"{prefix}.bwt.2bit.64").exists() or Path(f"{prefix}.0123").exists(): return
            cmd = ['bwa-mem2', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bowtie2':
            if Path(f"{prefix}.1.bt2").exists() or Path(f"{prefix}.1.bt2l").exists(): return
            cmd = ['bowtie2-build', '--threads', str(self.args.threads), str(ref_path), str(prefix)]
        elif tool == 'hisat2':
            if Path(f"{prefix}.1.ht2").exists() or Path(f"{prefix}.1.ht2l").exists(): return
            cmd = ['hisat2-build', '-p', str(self.args.threads), str(ref_path), str(prefix)]
            
        logger.info(f"🏗️ 构建 {tool} 索引...")
        with Timer(f"构建 {tool} 索引"):
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"❌ 索引构建失败: {result.stderr}")
                sys.exit(1)
            logger.info(f"✅ {tool} 索引构建完成")

    def find_samples(self):
        """查找样本文件"""
        input_dir = Path(self.args.input_dir) if self.args.input_dir else None
        fastq_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
        samples = []
        
        if input_dir and self.args.single_end:
            logger.info("🔍 搜索单端测序样本...")
            for ext in fastq_extensions:
                fastq_files = list(input_dir.glob(f'*{ext}'))
                for fastq in tqdm(fastq_files, desc="扫描文件", unit="文件"):
                    sample_name = fastq.name.replace(ext, '').replace('.gz', '')
                    for suffix in ['_1', '_R1', '.R1', "R1_10239",'.1', '_unmapped', '.unmapped', '_trimmed', '.trimmed']:
                        if suffix in sample_name: sample_name = sample_name.split(suffix)[0]
                    samples.append({'name': sample_name, 'r1': str(fastq), 'r2': None})
        elif input_dir:
            logger.info("🔍 搜索双端测序样本...")
            patterns = [
                ('*_1.f*q*', '*_2.f*q*'), ('*_R1*.f*q*', '*_R2*.f*q*'), ('*.R1.*', '*.R2.*'), 
                ('*.1.f*q*', '*.2.f*q*'), ('*_1_*.f*q*', '*_2_*.f*q*'), 
                ('*_unmapped.R1.fq.gz', '*_unmapped.R2.fq.gz'), ('*unmapped.R1.fq.gz', '*unmapped.R2.fq.gz'), 
                ('*_1_unmapped.*', '*_2_unmapped.*'), ('*.unmapped.R1_10239.*', '*.unmapped.R2_10239.*')
            ]
            found_files = set()
            for r1_pattern, r2_pattern in tqdm(patterns, desc="匹配文件模式", unit="模式"):
                r1_files = list(input_dir.glob(r1_pattern))
                for r1_file in r1_files:
                    if r1_file in found_files: continue
                    r1_name = r1_file.name
                    r2_name = None
                    if '_1.' in r1_name: r2_name = r1_name.replace('_1.', '_2.')
                    elif '_R1' in r1_name: r2_name = r1_name.replace('_R1', '_R2')
                    elif '.R1.' in r1_name: r2_name = r1_name.replace('.R1.', '.R2.')
                    elif '.1.' in r1_name: r2_name = r1_name.replace('.1.', '.2.')
                    elif '_1_' in r1_name: r2_name = r1_name.replace('_1_', '_2_')
                    elif '.unmapped.R1.' in r1_name: r2_name = r1_name.replace('.unmapped.R1.', '.unmapped.R2.')
                    elif 'unmapped.R1.' in r1_name: r2_name = r1_name.replace('unmapped.R1.', 'unmapped.R2.')
                    elif '_1_unmapped.' in r1_name: r2_name = r1_name.replace('_1_unmapped.', '_2_unmapped.')
                    elif '.unmapped.R1_10239.' in r1_name: r2_name = r1_name.replace('.unmapped.R1_10239.', '.unmapped.R2_10239.')
                    
                    if r2_name:
                        r2_file = input_dir / r2_name
                        if r2_file.exists():
                            sample_name = r1_name
                            for suffix in ['_1', '_R1', '.R1', '.1', '_unmapped', '.unmapped', '_1_unmapped']:
                                if suffix in sample_name: sample_name = sample_name.split(suffix)[0]
                            sample_name = sample_name.split('.')[0]
                            samples.append({'name': sample_name, 'r1': str(r1_file), 'r2': str(r2_file)})
                            found_files.add(r1_file)
                            found_files.add(r2_file)
            
            if not samples:
                all_fastq_files = []
                for ext in fastq_extensions: all_fastq_files.extend(input_dir.glob(f'*{ext}'))
                file_groups = defaultdict(list)
                for f in all_fastq_files:
                    prefix = re.sub(r'[._](R?[12]|unmapped\.R?[12])[._].*', '', f.name)
                    if '.' in prefix: prefix = prefix.split('.')[0]
                    file_groups[prefix].append(str(f))
                for prefix, files in file_groups.items():
                    if len(files) == 2:
                        r1, r2 = None, None
                        for f in files:
                            if any(pattern in f for pattern in ['_1', '_R1', '.R1.', '.1.', 'unmapped.R1']): r1 = f
                            elif any(pattern in f for pattern in ['_2', '_R2', '.R2.', '.2.', 'unmapped.R2']): r2 = f
                        if r1 and r2: samples.append({'name': prefix, 'r1': r1, 'r2': r2})
        
        if self.args.sample_list:
            logger.info(f"📋 从文件加载样本列表: {self.args.sample_list}")
            with open(self.args.sample_list, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            if not Path(parts[1]).exists(): continue
                            samples.append({
                                'name': parts[0],
                                'r1': parts[1],
                                'r2': parts[2] if len(parts) > 2 and Path(parts[2]).exists() else None
                            })
        
        unique_samples, seen_names = [], set()
        for sample in samples:
            if sample['name'] not in seen_names:
                unique_samples.append(sample)
                seen_names.add(sample['name'])
        
        self.samples = unique_samples
        if not self.samples:
            logger.error("❌ 未找到任何样本文件。请检查输入路径或使用正确的命名模式。")
            sys.exit(1)
            
        logger.info(f"✅ 找到 {len(self.samples)} 个样本")
    
    def align_sample(self, sample):
        """比对单个样本"""
        try:
            sample_name = sample['name']
            bam_dir = self.output_dir / 'bam'
            sorted_bam = bam_dir / f'{sample_name}.sorted.bam'
            
            if sorted_bam.exists() and self.args.skip_alignment:
                logger.debug(f"⏭️  跳过已存在的比对文件: {sorted_bam}")
                return sorted_bam
            
            logger.info(f"🔄 开始比对样本: {sample_name} (使用 {self.args.tool})")
            with Timer(f"比对 {sample_name}"):
                tool = self.args.tool
                align_cmd = []
                
                if tool == 'strobealign':
                    align_cmd = ['strobealign', '-t', str(self.args.align_threads), self.index_path, sample['r1']]
                    if sample['r2']: align_cmd.append(sample['r2'])
                elif tool == 'minimap2':
                    preset = 'sr' if sample['r2'] else 'map-ont'
                    align_cmd = ['minimap2', '-ax', preset, '-t', str(self.args.align_threads), self.index_path, sample['r1']]
                    if sample['r2']: align_cmd.append(sample['r2'])
                elif tool == 'bwa':
                    align_cmd = ['bwa', 'mem', '-v', '1', '-t', str(self.args.align_threads), self.index_path, sample['r1']]
                    if sample['r2']: align_cmd.append(sample['r2'])
                elif tool == 'bwa-mem2':
                    align_cmd = ['bwa-mem2', 'mem', '-v', '1', '-t', str(self.args.align_threads), self.index_path, sample['r1']]
                    if sample['r2']: align_cmd.append(sample['r2'])
                elif tool == 'bowtie2':
                    align_cmd = ['bowtie2', '-p', str(self.args.align_threads), '-x', self.index_path]
                    if sample['r2']: align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']])
                    else: align_cmd.extend(['-U', sample['r1']])
                elif tool == 'hisat2':
                    align_cmd = ['hisat2', '-p', str(self.args.align_threads), '-x', self.index_path]
                    if sample['r2']: align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']])
                    else: align_cmd.extend(['-U', sample['r1']])
                
                sort_cmd = ['samtools', 'sort', '-@', str(min(4, self.args.threads)), '-o', str(sorted_bam)]
                
                align_process = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                sort_process = subprocess.Popen(sort_cmd, stdin=align_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                align_process.stdout.close()
                _, align_stderr = align_process.communicate()
                sort_stdout, sort_stderr = sort_process.communicate()
                
                if align_process.returncode != 0:
                    logger.error(f"{tool} 比对失败: {align_stderr[:500]}")
                    return None
                if sort_process.returncode != 0:
                    logger.error(f"排序失败: {sort_stderr[:500]}")
                    return None
            
            with Timer(f"索引 {sample_name}"):
                result = subprocess.run(['samtools', 'index', str(sorted_bam)], capture_output=True, text=True)
                if result.returncode != 0:
                    logger.error(f"索引创建失败: {result.stderr}")
                    return None
            
            logger.info(f"✅ 样本 {sample_name} 比对完成")
            return sorted_bam
        
        except Exception as e:
            logger.error(f"❌ 样本 {sample['name']} 比对失败: {str(e)}")
            return None
            
    def get_idxstats(self, bam_file):
        """提取 BAM 中的绝对序列总数并计算 Total RPK 以用于后续 TPM/FPKM 计算"""
        try:
            cmd = ['samtools', 'idxstats', str(bam_file)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"samtools idxstats 失败: {result.stderr[:200]}")
                return None
                
            stats = {}
            total_reads_in_library = 0
            total_rpk = 0.0
            
            for line in result.stdout.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 4:
                    ref, length, mapped, unmapped = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
                    stats[ref] = mapped
                    
                    # Total Library Size = 样本中所有比对上 + 未比对上的 reads 总和
                    total_reads_in_library += (mapped + unmapped)
                    
                    # 计算当前病毒的 RPK，并累加到 Total RPK（用于后续算 TPM）
                    if length > 0 and mapped > 0:
                        total_rpk += (mapped * 1000.0) / length
                        
            return {
                'stats': stats, 
                'global_total_reads': total_reads_in_library,
                'total_rpk': total_rpk
            }
        except Exception as e:
            logger.error(f"❌ 提取 idxstats 失败: {str(e)}")
            return None
    
    def run_pandepth(self, bam_file):
        """运行 pandepth 计算覆盖度"""
        try:
            bam_path = Path(bam_file)
            sample_name = bam_path.stem.replace('.sorted', '')
            output_prefix = self.output_dir / 'stat' / sample_name
            stat_file = output_prefix.with_suffix('.chr.stat.gz')
            
            if stat_file.exists() and self.args.skip_depth:
                return stat_file
                
            with Timer(f"深度计算 {sample_name}"):
                pandepth_cmd = ['pandepth', '-a', '-i', str(bam_file), '-o', str(output_prefix), '-t', str(min(10, self.args.threads))]
                result = subprocess.run(pandepth_cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.error(f"pandepth运行失败: {result.stderr[:500]}")
                    return None
            return stat_file
        except Exception as e:
            logger.error(f"❌ pandepth运行出错: {str(e)}")
            return None
    
    def parse_stat_file(self, stat_file, sample_name, idx_data):
        """解析统计文件并计算 TPM, FPKM, RPM 三大指标"""
        if not idx_data:
            return []
            
        global_total_reads = idx_data['global_total_reads']
        total_rpk = idx_data['total_rpk']
        idx_stats = idx_data['stats']
        results = []
        
        try:
            with gzip.open(stat_file, 'rt') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'): continue
                    
                    parts = line.split()
                    if len(parts) < 6: continue
                    virus_name, length, covered_site, total_depth, coverage, mean_depth = parts
                    
                    length_val = int(length)
                    covered_site_val = int(covered_site)
                    total_depth_val = int(total_depth)
                    coverage_val = float(coverage)
                    mean_depth_val = float(mean_depth)
                    
                    # 取出实际比对上的 Reads 数量
                    mapped_reads = idx_stats.get(virus_name, 0)
                    
                    # --------------------------------------------------
                    # 1. RPM 和 FPKM (以整个测序文库的 absolute size 作为分母)
                    # --------------------------------------------------
                    if global_total_reads > 0 and length_val > 0:
                        rpm = (mapped_reads * 1e6) / global_total_reads
                        fpkm = (mapped_reads * 1e9) / (global_total_reads * length_val)
                    else:
                        rpm = 0.0
                        fpkm = 0.0
                        
                    # --------------------------------------------------
                    # 2. TPM (以该样本内所有病毒的 Total RPK 作为分母，用于群落结构比较)
                    # --------------------------------------------------
                    if total_rpk > 0 and length_val > 0:
                        rpk = (mapped_reads * 1000.0) / length_val
                        tpm = (rpk * 1e6) / total_rpk
                    else:
                        tpm = 0.0
                        
                    # 应用阈值过滤
                    if coverage_val > self.args.coverage and mean_depth_val > self.args.meandepth:
                        if fpkm >= self.args.min_fpkm and tpm >= self.args.min_tpm:
                            results.append({
                                'Sample': sample_name,
                                'Virus': virus_name,
                                'Length': length_val,
                                'CoveredSite': covered_site_val,
                                'TotalDepth': total_depth_val,
                                'Coverage(%)': coverage_val,
                                'MeanDepth': mean_depth_val,
                                'MappedReads': mapped_reads,
                                'RPM': round(rpm, 2),
                                'FPKM': round(fpkm, 2),
                                'TPM': round(tpm, 2)
                            })
            
            logger.debug(f"样本 {sample_name}: 找到 {len(results)} 条满足条件的记录")
            return results
        except Exception as e:
            logger.error(f"❌ 解析统计文件 {stat_file} 失败: {str(e)}")
            return []
    
    def process_sample(self, sample):
        """处理单个样本的完整流程"""
        try:
            logger.info(f"\n📋 处理样本: {sample['name']}")
            
            # 1. 比对
            bam_file = self.align_sample(sample)
            if not bam_file: return None
            
            # 2. 提取全局绝对比对统计 (计算 Global Reads 和 Total RPK)
            idx_data = self.get_idxstats(bam_file)
            if not idx_data: return None
            
            # 3. 计算覆盖度
            stat_file = self.run_pandepth(bam_file)
            if not stat_file: return None
            
            # 4. 解析并定量计算
            sample_results = self.parse_stat_file(stat_file, sample['name'], idx_data)
            
            logger.info(f"✅ 样本 {sample['name']} 处理完成")
            return sample_results
        
        except Exception as e:
            logger.error(f"❌ 处理样本 {sample['name']} 失败: {str(e)}")
            return None
    
    def process_samples_parallel(self, max_workers):
        all_results = []
        logger.info(f"🔀 并行处理样本，最大并行任务数: {max_workers}")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_sample = {executor.submit(self.process_sample, s): s for s in self.samples}
            with tqdm(total=len(self.samples), desc="处理样本", unit="样本") as pbar:
                for future in as_completed(future_to_sample):
                    sample = future_to_sample[future]
                    try:
                        result = future.result()
                        if result is not None:
                            all_results.extend(result)
                            pbar.set_postfix_str(f"✅ {sample['name']}")
                        else:
                            pbar.set_postfix_str(f"❌ {sample['name']}失败")
                    except Exception as e:
                        logger.error(f"❌ 样本 {sample['name']} 异常: {str(e)}")
                        pbar.set_postfix_str(f"❌ {sample['name']}异常")
                    pbar.update(1)
        return all_results
    
    def process_samples_serial(self):
        all_results = []
        with tqdm(self.samples, desc="处理样本", unit="样本") as pbar:
            for sample in pbar:
                pbar.set_postfix_str(f"正在处理 {sample['name']}")
                result = self.process_sample(sample)
                if result is not None:
                    all_results.extend(result)
                    pbar.set_postfix_str(f"✅ {sample['name']}")
                else:
                    pbar.set_postfix_str(f"❌ {sample['name']}失败")
        return all_results
    
    def run_pipeline(self):
        start_time = datetime.now()
        logger.info("=" * 60)
        logger.info(f"🚀 开始病毒分析流程 (工具: {self.args.tool} | 支持多维度定量 RPM/FPKM/TPM)")
        logger.info("=" * 60)
        
        self.setup_output_directory()
        self.build_index()
        self.find_samples()
        
        if self.args.parallel and self.args.threads > 1:
            max_parallel_jobs = max(1, self.args.threads // self.args.align_threads)
            max_parallel_jobs = min(max_parallel_jobs, len(self.samples))
            if self.args.parallel_jobs: max_parallel_jobs = min(self.args.parallel_jobs, len(self.samples))
            all_results = self.process_samples_parallel(max_parallel_jobs)
        else:
            all_results = self.process_samples_serial()
        
        if all_results:
            self.save_results(all_results)
        else:
            logger.warning("⚠️  未找到满足条件的记录")
        
        total_duration = datetime.now() - start_time
        logger.info("=" * 60)
        logger.info("✨ 分析流程完成")
        # 🔥🔥🔥 此处修复了方法调用：改用 Timer.format_duration 🔥🔥🔥
        logger.info(f"⏱️  总耗时: {Timer.format_duration(total_duration.total_seconds())}")
        if all_results:
            logger.info(f"💾 汇总结果文件已保存至: {self.output_dir}/summary/")
        logger.info("=" * 60)
    
    def save_results(self, all_results):
        """整合分类信息并保存结果"""
        logger.info(f"\n💾 正在合并分类信息并保存结果...")
        with Timer("汇总和保存结果"):
            df_all = pd.DataFrame(all_results)
            
            # 【合并 RVDB 分类注释信息】
            if self.taxon_db is not None and 'Virus' in df_all.columns:
                def get_tax_info(virus_id, key):
                    if virus_id in self.taxon_db:
                        return self.taxon_db[virus_id][key]
                    base_id = virus_id.split('.')[0]
                    if base_id in self.taxon_db:
                        return self.taxon_db[base_id][key]
                    return "Unannotated"
                    
                df_all['description'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'description'))
                df_all['taxonomy'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy'))
                df_all['taxonomy_id'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy_id'))
                
                # 整理列顺序：把注释列放到 Virus 的后面
                cols = list(df_all.columns)
                for col in ['taxonomy_id', 'taxonomy', 'description']:
                    if col in cols:
                        cols.remove(col)
                        v_idx = cols.index('Virus')
                        cols.insert(v_idx + 1, col)
                df_all = df_all[cols]

            summary_dir = self.output_dir / 'summary'
            summary_file = summary_dir / f'all_viruses.summary.{self.args.format}'
            
            if self.args.format == 'csv': df_all.to_csv(summary_file, index=False)
            else: df_all.to_csv(summary_file, sep='\t', index=False)
            
            if 'Virus' in df_all.columns:
                viruses = df_all['Virus'].unique()
                for virus in tqdm(viruses, desc="生成各病毒独立报表", unit="病毒"):
                    df_virus = df_all[df_all['Virus'] == virus].copy()
                    virus_safe_name = re.sub(r'[\\/*?:"<>|]', "_", virus)
                    virus_file = summary_dir / f'{virus_safe_name}.summary.{self.args.format}'
                    if self.args.format == 'csv': df_virus.to_csv(virus_file, index=False)
                    else: df_virus.to_csv(virus_file, sep='\t', index=False)
            
            self.generate_report(df_all)
    
    def generate_report(self, df_all):
        report_file = self.output_dir / 'summary' / 'analysis_report.txt'
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"病毒丰度分析报告 (多维度定量版: TPM / FPKM / RPM)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"比对工具: {self.args.tool}\n")
            f.write(f"参考基因组: {self.args.reference}\n")
            if self.taxon_db:
                f.write(f"使用RVDB注释文件: {self.args.rvdb_taxon}\n")
            f.write(f"过滤条件: Coverage > {self.args.coverage}%, MeanDepth > {self.args.meandepth}\n")
            f.write(f"定量过滤: FPKM >= {self.args.min_fpkm}, TPM >= {self.args.min_tpm}\n\n")
            
            f.write(f"总样本数: {len(self.samples)}\n")
            f.write(f"总发现记录数: {len(df_all)}\n\n")
            
            if 'Virus' in df_all.columns:
                f.write("常见病毒概览 (平均定量指标):\n")
                f.write("-" * 60 + "\n")
                virus_counts = df_all['Virus'].value_counts()
                for virus, count in virus_counts.items():
                    vd = df_all[df_all['Virus'] == virus]
                    desc = vd['description'].iloc[0] if 'description' in vd.columns else ''
                    mean_tpm = vd['TPM'].mean()
                    mean_fpkm = vd['FPKM'].mean()
                    mean_rpm = vd['RPM'].mean()
                    f.write(f"{virus} ({desc[:40]}...):\n")
                    f.write(f"  样本检出数: {count}\n")
                    f.write(f"  平均 TPM:  {mean_tpm:.2f} (相对群落丰度)\n")
                    f.write(f"  平均 FPKM: {mean_fpkm:.2f} (绝对载量，长度归一化)\n")
                    f.write(f"  平均 RPM:  {mean_rpm:.2f} (绝对载量)\n\n")
        logger.info(f"📄 分析报告已保存: {report_file}")
    
def main():
    parser = argparse.ArgumentParser(
        description='病毒序列比对、深度统计和多维度定量分析流程 (支持多种比对工具与 RVDB 分类注释)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
定量计算说明:
  1. RPM / FPKM: 分母为该样本的 **总测序 Reads 数 (Total Library Size)**。用于横向跨样本比较真实的**病毒绝对载量**。
  2. TPM: 分母为该样本内比对上的**所有病毒 RPK 总和**。用于消除不同样本数据量差异，比较病毒在宏基因组中的**相对结构占比**。

使用示例:
  # 基础分析 (使用 bwa-mem2)
  python batch_virus_depth.py --tool bwa-mem2 --input_dir fastq --reference virus.fasta
  
  # 附加 RVDB 注释并在输出结果中合并
  python batch_virus_depth.py --tool bwa-mem2 --input_dir fastq --reference virus.fasta --rvdb_taxon RVDB_Taxon_Current.tab
        """
    )
    
    # 输入输出参数
    input_group = parser.add_argument_group('输入输出参数')
    input_group.add_argument('--input_dir', type=str, help='输入文件夹路径（包含fastq文件）')
    input_group.add_argument('--reference', type=str, required=True, help='参考基因组fasta文件路径')
    input_group.add_argument('--rvdb_taxon', type=str, help='RVDB 分类注释文件路径 (如: RVDB_Taxon_Current.tab)')
    input_group.add_argument('--output_dir', type=str, default='./virus_analysis', help='输出目录路径 [默认: ./virus_analysis]')
    input_group.add_argument('--sample_list', type=str, help='样本列表文件（每行: 样本名 fastq1 fastq2）')
    input_group.add_argument('--single_end', action='store_true', help='使用单端测序数据')
    
    # 比对参数
    align_group = parser.add_argument_group('比对参数')
    align_group.add_argument('--tool', type=str, choices=['strobealign', 'minimap2', 'bwa', 'bwa-mem2', 'bowtie2', 'hisat2'], default='strobealign', help='比对工具选择 [默认: strobealign]')
    align_group.add_argument('--align_threads', type=int, default=8, help='每个样本比对的线程数 [默认: 8]')
    
    # 过滤参数
    filter_group = parser.add_argument_group('过滤参数')
    filter_group.add_argument('--coverage', type=float, default=90.0, help='Coverage(%%)阈值 [默认: 90]')
    filter_group.add_argument('--meandepth', type=float, default=10.0, help='MeanDepth阈值 [默认: 10]')
    filter_group.add_argument('--min_fpkm', type=float, default=0.0, help='FPKM最小值过滤 [默认: 0]')
    filter_group.add_argument('--min_tpm', type=float, default=0.0, help='TPM最小值过滤 [默认: 0]')
    
    # 并行处理参数
    parallel_group = parser.add_argument_group('并行处理参数')
    parallel_group.add_argument('--parallel', action='store_true', help='并行处理样本')
    parallel_group.add_argument('--parallel_jobs', type=int, help='并行任务数（如果不指定，自动计算）')
    parallel_group.add_argument('--threads', type=int, default=4, help='总线程数 [默认: 4]')
    
    # 流程控制参数
    control_group = parser.add_argument_group('流程控制参数')
    control_group.add_argument('--skip_alignment', action='store_true', help='跳过比对步骤（使用已有BAM文件）')
    control_group.add_argument('--skip_depth', action='store_true', help='跳过覆盖度计算步骤（使用已有统计文件）')
    control_group.add_argument('--format', type=str, choices=['csv', 'tsv'], default='tsv', help='输出文件格式 [默认: tsv]')
    
    # 其他参数
    other_group = parser.add_argument_group('其他参数')
    other_group.add_argument('--verbose', action='store_true', help='输出详细日志信息')
    
    args = parser.parse_args()
    if not args.input_dir and not args.sample_list: parser.error("必须提供 --input_dir 或 --sample_list 参数")
    if not os.path.exists(args.reference): parser.error(f"参考基因组文件不存在: {args.reference}")
    if args.input_dir and not os.path.exists(args.input_dir): parser.error(f"输入目录不存在: {args.input_dir}")
    if args.verbose: logger.setLevel(logging.DEBUG)
    
    pipeline = VirusAnalysisPipeline(args)
    pipeline.run_pipeline()

if __name__ == '__main__':
    main()
