#!/usr/bin/env python3
"""
病毒序列比对、深度统计、多维度定量分析及共识序列构建完整流程
支持多种比对工具，支持极速解析超大型 RVDB 分类注释文件
支持调用 viral_consensus 自动构建共识序列并按病毒种类归档保存
"""

import argparse
import os
import sys
import subprocess
import glob
import gzip
import shutil
import time
import math
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
        if seconds < 60: return f"{seconds:.1f}秒"
        elif seconds < 3600: return f"{seconds/60:.1f}分钟"
        else: return f"{seconds/3600:.1f}小时"

class VirusAnalysisPipeline:
    def __init__(self, args):
        self.args = args
        self.samples = []
        self.results = []
        self.output_dir = Path(args.output_dir)
        self.index_path = None
        
        self.check_tools()
    
    def check_tools(self):
        required_tools = {
            'samtools': '测序数据处理',
            'pandepth': '深度计算',
            self.args.tool: '序列比对'
        }
        
        if self.args.tool == 'bowtie2': required_tools['bowtie2-build'] = 'Bowtie2索引构建'
        elif self.args.tool == 'hisat2': required_tools['hisat2-build'] = 'HISAT2索引构建'
        if self.args.consensus:
            required_tools['viral_consensus'] = '共识序列构建'
            required_tools['awk'] = 'BAM头处理工具'
            
        missing_tools = []
        for tool, description in required_tools.items():
            if shutil.which(tool) is None:
                missing_tools.append(f"{tool} ({description})")
        
        if missing_tools:
            logger.error("❌ 缺少必要的工具:")
            for tool in missing_tools:
                logger.error(f"  - {tool}")
            sys.exit(1)
        logger.info("✅ 所有必要工具可用")

    def fast_load_rvdb_taxonomy(self, found_viruses):
        """流式极速提取 RVDB 数据，解决 1000万 行文件的慢速问题"""
        if not self.args.rvdb_taxon or not found_viruses:
            return None
            
        logger.info(f"📂 正在快速流式扫描 RVDB 注释文件 (匹配 {len(found_viruses)} 个检出病毒)...")
        
        target_exact = set(found_viruses)
        target_base = {v.split('.')[0] for v in found_viruses}
        
        taxon_dict = {}
        found_count = 0
        
        try:
            with open(self.args.rvdb_taxon, 'r', encoding='utf-8') as f:
                header = None
                idx_acc, idx_desc, idx_tax, idx_taxid = -1, -1, -1, -1
                
                for line in f:
                    if line.startswith('####'):
                        continue
                    
                    parts = line.rstrip('\n').split('\t')
                    
                    if header is None:
                        header = parts
                        try:
                            idx_acc = header.index('accession')
                            idx_desc = header.index('description')
                            idx_tax = header.index('taxonomy')
                            idx_taxid = header.index('taxonomy_id')
                        except ValueError:
                            logger.warning("⚠️ RVDB 缺少必要的列头(accession, description等)，跳过注释")
                            return None
                        continue
                    
                    if len(parts) <= max(idx_acc, idx_desc, idx_tax, idx_taxid):
                        continue
                        
                    acc = parts[idx_acc].strip()
                    acc_base = acc.split('.')[0]
                    
                    if acc in target_exact or acc_base in target_base or acc in target_base:
                        info = {
                            'description': parts[idx_desc] if parts[idx_desc] else '',
                            'taxonomy': parts[idx_tax] if parts[idx_tax] else '',
                            'taxonomy_id': parts[idx_taxid] if parts[idx_taxid] else ''
                        }
                        taxon_dict[acc] = info
                        taxon_dict[acc_base] = info
                        found_count += 1
                        
            logger.info(f"✅ 扫描完成，成功提取到 {found_count} 条分类信息")
            return taxon_dict
            
        except Exception as e:
            logger.error(f"❌ 读取 RVDB 失败: {e}")
            return None
    
    def setup_output_directory(self):
        subdirs = ['bam', 'stat', 'summary', 'logs', 'index']
        if self.args.consensus:
            subdirs.append('consensus')
            
        with Timer("创建输出目录结构"):
            for subdir in subdirs:
                (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)
                
            # --- 增加功能: 动态添加 FileHandler 同步保存全量日志 ---
            has_file_handler = any(isinstance(h, logging.FileHandler) for h in logger.handlers)
            if not has_file_handler:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                log_file = self.output_dir / 'logs' / f'analysis_{timestamp}.log'
                file_handler = logging.FileHandler(log_file, encoding='utf-8')
                file_handler.setFormatter(logging.Formatter(
                    '%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S'
                ))
                logger.addHandler(file_handler)
                logger.info(f"📝 运行日志已同步保存至: {log_file}")
                
        logger.info(f"📁 输出目录: {self.output_dir}")
        return self.output_dir
        
    def build_index(self):
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
        input_dir = Path(self.args.input_dir) if self.args.input_dir else None
        fastq_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
        samples = []
        
        if input_dir and self.args.single_end:
            logger.info("🔍 搜索单端测序样本...")
            for ext in fastq_extensions:
                for fastq in input_dir.glob(f'*{ext}'):
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
            for r1_pattern, r2_pattern in patterns:
                for r1_file in input_dir.glob(r1_pattern):
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
                        if len(parts) >= 2 and Path(parts[1]).exists():
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
        try:
            sample_name = sample['name']
            bam_dir = self.output_dir / 'bam'
            sorted_bam = bam_dir / f'{sample_name}.sorted.bam'
            
            if sorted_bam.exists() and self.args.skip_alignment:
                return sorted_bam
            
            align_cmd = []
            tool = self.args.tool
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
            
            align_proc = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            sort_proc = subprocess.Popen(sort_cmd, stdin=align_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            align_proc.stdout.close()
            _, align_stderr = align_proc.communicate()
            sort_stdout, sort_stderr = sort_proc.communicate()
            
            if align_proc.returncode != 0:
                logger.error(f"{tool} 比对失败: {align_stderr[:500]}")
                return None
            if sort_proc.returncode != 0:
                logger.error(f"排序失败: {sort_stderr[:500]}")
                return None
            
            subprocess.run(['samtools', 'index', str(sorted_bam)], capture_output=True)
            return sorted_bam
        except Exception as e:
            logger.error(f"❌ 比对异常: {str(e)}")
            return None
            
    def get_idxstats(self, bam_file):
        try:
            cmd = ['samtools', 'idxstats', str(bam_file)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0: return None
                
            stats = {}
            total_reads_in_library = 0
            total_rpk = 0.0
            for line in result.stdout.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 4:
                    ref, length, mapped, unmapped = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
                    stats[ref] = mapped
                    total_reads_in_library += (mapped + unmapped)
                    if length > 0 and mapped > 0:
                        total_rpk += (mapped * 1000.0) / length
            return {'stats': stats, 'global_total_reads': total_reads_in_library, 'total_rpk': total_rpk}
        except:
            return None
    
    def run_pandepth(self, bam_file):
        try:
            bam_path = Path(bam_file)
            sample_name = bam_path.stem.replace('.sorted', '')
            output_prefix = self.output_dir / 'stat' / sample_name
            stat_file = output_prefix.with_suffix('.chr.stat.gz')
            if stat_file.exists() and self.args.skip_depth: return stat_file
                
            pandepth_cmd = ['pandepth', '-a', '-i', str(bam_file), '-o', str(output_prefix), '-t', str(min(10, self.args.threads))]
            subprocess.run(pandepth_cmd, capture_output=True)
            return stat_file
        except:
            return None
    
    def parse_stat_file(self, stat_file, sample_name, idx_data):
        if not idx_data: return []
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
                    
                    length_val, covered_site_val, total_depth_val = int(length), int(covered_site), int(total_depth)
                    coverage_val, mean_depth_val = float(coverage), float(mean_depth)
                    mapped_reads = idx_stats.get(virus_name, 0)
                    
                    if global_total_reads > 0 and length_val > 0:
                        rpm = (mapped_reads * 1e6) / global_total_reads
                        fpkm = (mapped_reads * 1e9) / (global_total_reads * length_val)
                    else: rpm, fpkm = 0.0, 0.0
                        
                    if total_rpk > 0 and length_val > 0:
                        rpk = (mapped_reads * 1000.0) / length_val
                        tpm = (rpk * 1e6) / total_rpk
                    else: tpm = 0.0
                        
                    if coverage_val > self.args.coverage and mean_depth_val > self.args.meandepth:
                        if fpkm >= self.args.min_fpkm and tpm >= self.args.min_tpm:
                            results.append({
                                'Sample': sample_name, 'Virus': virus_name, 'Length': length_val,
                                'CoveredSite': covered_site_val, 'TotalDepth': total_depth_val,
                                'Coverage(%)': coverage_val, 'MeanDepth': mean_depth_val,
                                'MappedReads': mapped_reads, 'RPM': round(rpm, 2),
                                'FPKM': round(fpkm, 2), 'TPM': round(tpm, 2)
                            })
            return results
        except: return []
    
    def process_sample(self, sample):
        try:
            logger.info(f"📋 处理样本: {sample['name']}")
            bam_file = self.align_sample(sample)
            if not bam_file: return None
            idx_data = self.get_idxstats(bam_file)
            if not idx_data: return None
            stat_file = self.run_pandepth(bam_file)
            if not stat_file: return None
            return self.parse_stat_file(stat_file, sample['name'], idx_data)
        except Exception as e:
            logger.error(f"❌ 样本 {sample['name']} 失败: {str(e)}")
            return None
    
    def run_pipeline(self):
        start_time = datetime.now()
        logger.info("=" * 60)
        logger.info(f"🚀 开始病毒分析流程 (工具: {self.args.tool} | 定量: FPKM/RPM/TPM | 共识构建: {self.args.consensus})")
        logger.info("=" * 60)
        
        self.setup_output_directory()
        self.build_index()
        self.find_samples()
        
        all_results = []
        if self.args.parallel and self.args.threads > 1:
            max_parallel_jobs = min(max(1, self.args.threads // self.args.align_threads), len(self.samples))
            if self.args.parallel_jobs: max_parallel_jobs = min(self.args.parallel_jobs, len(self.samples))
            logger.info(f"🔀 并行处理样本，最大并行任务数: {max_parallel_jobs}")
            
            with ProcessPoolExecutor(max_workers=max_parallel_jobs) as executor:
                future_to_sample = {executor.submit(self.process_sample, s): s for s in self.samples}
                with tqdm(total=len(self.samples), desc="处理样本", unit="样本") as pbar:
                    for future in as_completed(future_to_sample):
                        res = future.result()
                        if res is not None: all_results.extend(res)
                        pbar.update(1)
        else:
            for sample in tqdm(self.samples, desc="处理样本", unit="样本"):
                res = self.process_sample(sample)
                if res is not None: all_results.extend(res)
        
        if all_results:
            self.save_results(all_results)
            # 开启了生成共识序列功能
            if self.args.consensus:
                self.build_consensus_sequences(all_results)
        else:
            logger.warning("⚠️ 未找到满足过滤条件的病毒记录")
        
        total_duration = datetime.now() - start_time
        logger.info("=" * 60)
        logger.info("✨ 分析流程完全结束")
        logger.info(f"⏱️  总耗时: {Timer.format_duration(total_duration.total_seconds())}")
        logger.info("=" * 60)
    
    def save_results(self, all_results):
        with Timer("汇总和保存结果"):
            df_all = pd.DataFrame(all_results)
            
            if self.args.rvdb_taxon and 'Virus' in df_all.columns:
                found_viruses = df_all['Virus'].unique().tolist()
                taxon_db = self.fast_load_rvdb_taxonomy(found_viruses)
                
                if taxon_db:
                    def get_tax_info(vid, key):
                        if vid in taxon_db: return taxon_db[vid][key]
                        base_id = vid.split('.')[0]
                        if base_id in taxon_db: return taxon_db[base_id][key]
                        return "Unannotated"
                        
                    df_all['description'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'description'))
                    df_all['taxonomy'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy'))
                    df_all['taxonomy_id'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy_id'))
                    
                    cols = list(df_all.columns)
                    for col in ['taxonomy_id', 'taxonomy', 'description']:
                        if col in cols: cols.remove(col)
                        cols.insert(cols.index('Virus') + 1, col)
                    df_all = df_all[cols]

            summary_dir = self.output_dir / 'summary'
            summary_file = summary_dir / f'all_viruses.summary.{self.args.format}'
            
            if self.args.format == 'csv': df_all.to_csv(summary_file, index=False)
            else: df_all.to_csv(summary_file, sep='\t', index=False)
            
            if 'Virus' in df_all.columns:
                for virus in df_all['Virus'].unique():
                    df_virus = df_all[df_all['Virus'] == virus].copy()
                    vname = re.sub(r'[\\/*?:"<>|]', "_", virus)
                    vfile = summary_dir / f'{vname}.summary.{self.args.format}'
                    if self.args.format == 'csv': df_virus.to_csv(vfile, index=False)
                    else: df_virus.to_csv(vfile, sep='\t', index=False)
            
            self.generate_report(df_all)

    # ================= 共识序列处理核心模块 =================
    
    def extract_fastas_with_python(self, target_ids, out_dir):
        """Python 内置极速 FASTA 提取器，并为每个病毒建立独立文件夹"""
        logger.info("🧬 正在提取参考 FASTA 并建立独立目录...")
        target_set = set(target_ids)
        found = set()
        
        current_id = None
        current_seq = []
        
        def save_record(vid, seq_list):
            if vid in target_set:
                # --- 增加功能: 为每个病毒建立专属子目录 ---
                vname = re.sub(r'[\\/*?:"<>|]', "_", vid)
                v_dir = out_dir / vname
                v_dir.mkdir(parents=True, exist_ok=True)
                
                with open(v_dir / f"{vname}.ref.fasta", 'w') as outf:
                    outf.write(f">{vid}\n" + "".join(seq_list) + "\n")
                found.add(vid)
                
        with open(self.args.reference, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id: save_record(current_id, current_seq)
                    # 处理带描述的 fasta 头 (只取第一段ID)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id: save_record(current_id, current_seq)
            
        if len(found) < len(target_set):
            logger.warning(f"⚠️ 部分病毒未在 FASTA 中找到序列: {target_set - found}")
            
    def _run_single_consensus(self, task):
        sample, virus, mean_depth = task
        consensus_dir = self.output_dir / 'consensus'
        
        # 使用安全的病毒名称定位专属子目录
        vname = re.sub(r'[\\/*?:"<>|]', "_", virus)
        virus_dir = consensus_dir / vname
        virus_dir.mkdir(parents=True, exist_ok=True)
        
        bam_file = self.output_dir / 'bam' / f"{sample}.sorted.bam"
        ref_fasta = virus_dir / f"{vname}.ref.fasta"
        fixed_bam = virus_dir / f"{sample}.{vname}.fixed.bam"
        out_fasta = virus_dir / f"{sample}.{vname}.consensus.fasta"
        
        if not ref_fasta.exists() or not bam_file.exists():
            return False
            
        # --- 增加功能: 智能截断 --min_depth (最大10，最小1) ---
        min_depth = int(math.floor(mean_depth))
        if min_depth > 10:
            min_depth = 10
        elif min_depth < 1:
            min_depth = 1
        
        try:
            # 1. 过滤冗余 @SQ 头部文件并提取病毒对应的比对信息
            awk_cmd = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
            pipe_cmd = f"samtools view -h '{bam_file}' '{virus}' | {awk_cmd} | samtools view -b -o '{fixed_bam}'"
            subprocess.run(pipe_cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
            
            # 2. 调用 viral_consensus 进行序列推断 (使用动态生成的 min_depth)
            vc_cmd = f"viral_consensus -i '{fixed_bam}' -r '{ref_fasta}' -o '{out_fasta}' --min_depth {min_depth}"
            subprocess.run(vc_cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
            
            # 3. 清理庞大的临时固定 BAM 文件，节省存储空间
            if fixed_bam.exists():
                fixed_bam.unlink()
                
            return True
        except subprocess.CalledProcessError as e:
            return False

    def build_consensus_sequences(self, all_results):
        """生成共识序列的大管家函数"""
        logger.info("\n🧬 开始生成病毒共识序列...")
        consensus_dir = self.output_dir / 'consensus'
        
        with Timer("共识序列构建"):
            # 1. 快速提取所有的参考基因组序列至专属目录
            unique_viruses = list(set(r['Virus'] for r in all_results))
            self.extract_fastas_with_python(unique_viruses, consensus_dir)
            
            # 2. 准备并行任务
            tasks = []
            for r in all_results:
                tasks.append((r['Sample'], r['Virus'], r['MeanDepth']))
                
            # 3. 使用多进程提速
            success_count = 0
            with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_consensus, t): t for t in tasks}
                with tqdm(total=len(tasks), desc="构建进度") as pbar:
                    for future in as_completed(futures):
                        if future.result(): success_count += 1
                        pbar.update(1)
            
            logger.info(f"✅ 共成功生成 {success_count} 个病毒共识序列文件")

    def generate_report(self, df_all):
        report_file = self.output_dir / 'summary' / 'analysis_report.txt'
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"病毒丰度分析报告 (多维度定量版: TPM / FPKM / RPM)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"比对工具: {self.args.tool}\n")
            f.write(f"参考基因组: {self.args.reference}\n")
            if self.args.rvdb_taxon:
                f.write(f"使用RVDB注释文件: {self.args.rvdb_taxon}\n")
            f.write(f"过滤条件: Coverage > {self.args.coverage}%, MeanDepth > {self.args.meandepth}\n")
            f.write(f"定量过滤: FPKM >= {self.args.min_fpkm}, TPM >= {self.args.min_tpm}\n")
            f.write(f"共识序列: {'已生成 (按病毒分子目录保存在 consensus 文件夹下)' if self.args.consensus else '否'}\n\n")
            
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
        description='病毒序列比对、深度统计、多维度定量分析及共识序列生成流程',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
增强功能说明:
  --consensus: 开启共识序列构建。它会解析过滤后的合格病毒记录，按病毒名建立独立目录。
               动态调节 --min_depth (MeanDepth>10设为10，否则设为向下取整)，生成规范的共识序列。

使用示例:
  # 基础分析 + RVDB 注释 + 自动生成带有子目录和智能截断深度的共有序列
  python batch_virus_depth.py --tool bwa-mem2 --input_dir fastq --reference virus.fasta \\
         --rvdb_taxon RVDB_Taxon_Current.tab --consensus
        """
    )
    
    input_group = parser.add_argument_group('输入输出参数')
    input_group.add_argument('--input_dir', type=str, help='输入文件夹路径（包含fastq文件）')
    input_group.add_argument('--reference', type=str, required=True, help='参考基因组fasta文件路径')
    input_group.add_argument('--rvdb_taxon', type=str, help='RVDB 分类注释文件路径 (如: RVDB_Taxon_Current.tab)')
    input_group.add_argument('--output_dir', type=str, default='./virus_analysis', help='输出目录路径 [默认: ./virus_analysis]')
    input_group.add_argument('--sample_list', type=str, help='样本列表文件（每行: 样本名 fastq1 fastq2）')
    input_group.add_argument('--single_end', action='store_true', help='使用单端测序数据')
    
    align_group = parser.add_argument_group('比对及共识序列参数')
    align_group.add_argument('--tool', type=str, choices=['strobealign', 'minimap2', 'bwa', 'bwa-mem2', 'bowtie2', 'hisat2'], default='strobealign', help='比对工具选择 [默认: strobealign]')
    align_group.add_argument('--align_threads', type=int, default=8, help='每个样本比对的线程数 [默认: 8]')
    align_group.add_argument('--consensus', action='store_true', help='对发现的病毒结果按其对应层级自动生成共识序列')
    
    filter_group = parser.add_argument_group('过滤参数')
    filter_group.add_argument('--coverage', type=float, default=90.0, help='Coverage(%%)阈值 [默认: 90]')
    filter_group.add_argument('--meandepth', type=float, default=10.0, help='MeanDepth阈值 [默认: 10]')
    filter_group.add_argument('--min_fpkm', type=float, default=0.0, help='FPKM最小值过滤 [默认: 0]')
    filter_group.add_argument('--min_tpm', type=float, default=0.0, help='TPM最小值过滤 [默认: 0]')
    
    parallel_group = parser.add_argument_group('并行处理参数')
    parallel_group.add_argument('--parallel', action='store_true', help='并行处理样本')
    parallel_group.add_argument('--parallel_jobs', type=int, help='并行任务数（如果不指定，自动计算）')
    parallel_group.add_argument('--threads', type=int, default=4, help='总线程数 [默认: 4]')
    
    control_group = parser.add_argument_group('流程控制参数')
    control_group.add_argument('--skip_alignment', action='store_true', help='跳过比对步骤（使用已有BAM文件）')
    control_group.add_argument('--skip_depth', action='store_true', help='跳过覆盖度计算步骤（使用已有统计文件）')
    control_group.add_argument('--format', type=str, choices=['csv', 'tsv'], default='tsv', help='输出文件格式 [默认: tsv]')
    
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
