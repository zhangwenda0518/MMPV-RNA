#!/usr/bin/env python3
"""
病毒序列比对、深度统计、多维度定量分析、共识序列构建与可视化完整流程
1. 包含 FPKM, RPM, TPM 三大定量体系，及 Relative_Abundance(%) 相对丰度补充
2. 包含 RVDB 极速流式解析 + NCBI Entrez 智能联网补全 (含实时缓存断点续跑)
3. 包含最优同源代表株筛选机制 (Best Representative by Taxonomy)
4. 包含阳性reads提取、变异检测及断点续跑(--resume)
5. 支持 --remap 独立重比对解剖模式，消除数据库竞争，达到极限组装质量
6. 共识模块支持双引擎切换: --consensus_tool {viral_consensus, ivar}
7. 完美兼容 5 大软件的 Unique / Multi 极速分离，源头剔除次优/嵌合假阳性
8. [新增] 报告中增加群体样本检出率统计，TSV 表格中增加 Unique(%) 占比列
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
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import colorlog
from tqdm import tqdm
import pandas as pd
import re

# 设置彩色全局日志
def setup_logging(verbose=False):
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        '%(log_color)s%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S',
        log_colors={'DEBUG': 'cyan', 'INFO': 'green', 'WARNING': 'yellow', 'ERROR': 'red', 'CRITICAL': 'bold_red'}
    ))
    logger = colorlog.getLogger(__name__)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    return logger

logger = setup_logging()

class Timer:
    def __init__(self, name=""): self.name = name
    def __enter__(self):
        self.start_time = time.time()
        logger.info(f"⏱️  开始: {self.name}")
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        duration = time.time() - self.start_time
        if exc_type is None: logger.info(f"✅ 完成: {self.name}[耗时: {self.format_duration(duration)}]")
        else: logger.error(f"❌ 失败: {self.name}[耗时: {self.format_duration(duration)}]")
    @staticmethod
    def format_duration(seconds):
        if seconds < 60: return f"{seconds:.1f}秒"
        elif seconds < 3600: return f"{seconds/60:.1f}分钟"
        else: return f"{seconds/3600:.1f}小时"

class VirusAnalysisPipeline:
    def __init__(self, args):
        self.args = args
        self.samples = []
        self.output_dir = Path(args.output_dir)
        self.index_path = None
        
        if self.args.resume:
            self.args.skip_alignment = True
            self.args.skip_depth = True

        self.check_tools()
        
    def check_tools(self):
        required_tools = {'samtools': '测序数据处理', 'pandepth': '深度计算', self.args.tool: '序列比对'}
        
        if self.args.tool == 'bowtie2': required_tools['bowtie2-build'] = 'Bowtie2索引构建'
        elif self.args.tool == 'hisat2': required_tools['hisat2-build'] = 'HISAT2索引构建'
        elif self.args.tool == 'bwa': required_tools['bwa'] = 'BWA索引构建'
        elif self.args.tool == 'bwa-mem2': required_tools['bwa-mem2'] = 'BWA-MEM2索引构建'
        
        if getattr(self.args, 'extract_reads', False):
            required_tools['pigz'] = '多线程压缩工具'
            
        if self.args.consensus:
            if self.args.consensus_tool == 'viral_consensus':
                required_tools['viral_consensus'] = '共识序列构建'
            elif self.args.consensus_tool == 'ivar':
                required_tools['ivar'] = '共识序列构建 (iVar)'
            required_tools['awk'] = 'BAM头处理工具'
            
        if getattr(self.args, 'call_variants', False):
            required_tools[self.args.variant_caller] = '变异检测工具'
            required_tools['awk'] = 'BAM头处理工具'
            
        missing_tools = [f"{t} ({d})" for t, d in required_tools.items() if shutil.which(t) is None]
        if missing_tools:
            logger.error("❌ 缺少必要的工具:")
            for m in missing_tools: logger.error(f"  - {m}")
            sys.exit(1)
        logger.info("✅ 所有必要工具可用")

    def write_sample_log(self, sample_name, content, mode='a'):
        log_file = self.output_dir / 'logs' / f"{sample_name}.log"
        with open(log_file, mode, encoding='utf-8') as f:
            f.write(content + "\n")

    def fast_load_rvdb_taxonomy(self, found_viruses):
        if not self.args.rvdb_taxon or not found_viruses: return None
        logger.info(f"📂 正在快速流式扫描本地 RVDB 注释文件...")
        target_exact = set(found_viruses)
        target_base = {v.split('.')[0] for v in found_viruses}
        taxon_dict = {}
        found_count = 0
        try:
            with open(self.args.rvdb_taxon, 'r', encoding='utf-8') as f:
                header = None
                idx_acc, idx_desc, idx_tax, idx_taxid = -1, -1, -1, -1
                for line in f:
                    if line.startswith('####'): continue
                    parts = line.rstrip('\n').split('\t')
                    if header is None:
                        header = parts
                        try:
                            idx_acc, idx_desc = header.index('accession'), header.index('description')
                            idx_tax, idx_taxid = header.index('taxonomy'), header.index('taxonomy_id')
                        except ValueError: return None
                        continue
                    if len(parts) <= max(idx_acc, idx_desc, idx_tax, idx_taxid): continue
                    acc = parts[idx_acc].strip()
                    acc_base = acc.split('.')[0]
                    if acc in target_exact or acc_base in target_base or acc in target_base:
                        info = {'description': parts[idx_desc] or '', 'taxonomy': parts[idx_tax] or '', 'taxonomy_id': parts[idx_taxid] or ''}
                        taxon_dict[acc], taxon_dict[acc_base] = info, info
                        found_count += 1
            logger.info(f"✅ 本地扫描完成，成功提取到 {found_count} 条分类信息")
            return taxon_dict
        except Exception as e:
            logger.error(f"❌ 读取 RVDB 失败: {e}")
            return None

    def fetch_ncbi_taxonomy(self, accessions):
        try:
            from Bio import Entrez, SeqIO
        except ImportError:
            logger.warning("⚠️ 找不到 biopython 模块，跳过 NCBI 在线解析。")
            return {}

        ncbi_dict = {}
        cache_file = self.output_dir / 'summary' / 'ncbi_taxonomy_cache.tsv'
        
        if self.args.resume and cache_file.exists():
            try:
                df_cache = pd.read_csv(cache_file, sep='\t', dtype=str)
                for _, row in df_cache.iterrows():
                    ncbi_dict[row['accession']] = {
                        'description': str(row['description']),
                        'taxonomy': str(row['taxonomy']),
                        'taxonomy_id': str(row['taxonomy_id'])
                    }
                logger.info(f"🔄 [断点续跑] 从本地缓存加载了 {len(ncbi_dict)} 条历史 NCBI 获取记录")
            except Exception as e:
                logger.warning(f"⚠️ 读取 NCBI 缓存失败，将重新获取: {e}")

        to_fetch = [acc for acc in accessions if acc not in ncbi_dict]
        
        if not to_fetch:
            logger.info("✅ 所有序列的 NCBI 分类信息已在缓存中找到，跳过联网获取。")
            return ncbi_dict

        Entrez.email = self.args.email
        if self.args.api_key:
            Entrez.api_key = self.args.api_key
            sleep_time = 0.15
        else:
            sleep_time = 0.35

        logger.info(f"🌐 正在通过 NCBI 批量获取 {len(to_fetch)} 个未知病毒分类 (批处理加速模式)...")

        write_header = not cache_file.exists()
        with open(cache_file, 'a', encoding='utf-8') as cf:
            if write_header:
                cf.write("accession\tdescription\ttaxonomy\ttaxonomy_id\n")
            
            # 【核心修复】：批量处理，每 50 个 ID 合并成一次网络请求
            batch_size = 50
            for i in tqdm(range(0, len(to_fetch), batch_size), desc="NCBI 批量解析进度", unit="批次"):
                batch = to_fetch[i:i+batch_size]
                id_list = ",".join(batch)
                
                max_retries = 3
                success = False
                for attempt in range(max_retries):
                    try:
                        # 一次性请求 50 条记录
                        handle = Entrez.efetch(db="nuccore", id=id_list, rettype="gb", retmode="text")
                        # 使用 parse 迭代解析多个结果
                        records = list(SeqIO.parse(handle, "genbank"))
                        handle.close()
                        
                        fetched_accs = set()
                        for record in records:
                            tax_id = "-"
                            for feature in record.features:
                                if feature.type == "source":
                                    for db_xref in feature.qualifiers.get('db_xref', []):
                                        if db_xref.startswith('taxon:'): tax_id = db_xref.split(':')[1]; break
                                    break
                                    
                            desc = record.description
                            tax = record.annotations.get('organism', '-')
                            
                            # 将 Biopython 解析出的 ID 映射回我们请求的原始 ID 列表
                            matched_acc = record.id
                            for b_acc in batch:
                                if b_acc == record.id or b_acc.split('.')[0] == record.id.split('.')[0]:
                                    matched_acc = b_acc
                                    break
                                    
                            ncbi_dict[matched_acc] = {'description': desc, 'taxonomy': tax, 'taxonomy_id': tax_id}
                            cf.write(f"{matched_acc}\t{desc}\t{tax}\t{tax_id}\n")
                            fetched_accs.add(matched_acc)
                            
                        # 标记那些 NCBI 没返回的“死链” ID，防止它们无限重试
                        for missing in batch:
                            if missing not in fetched_accs:
                                ncbi_dict[missing] = {'description': "Unannotated (NCBI Not Found)", 'taxonomy': "Unannotated", 'taxonomy_id': "-"}
                                cf.write(f"{missing}\tUnannotated (NCBI Not Found)\tUnannotated\t-\n")
                                
                        cf.flush()
                        success = True
                        break
                    except Exception as e:
                        time.sleep(2) # 网络失败则等待 2 秒后重试
                        
                if not success:
                    for b_acc in batch:
                        ncbi_dict[b_acc] = {'description': "Unannotated (NCBI Error)", 'taxonomy': "Unannotated", 'taxonomy_id': "-"}
                        cf.write(f"{b_acc}\tUnannotated (NCBI Error)\tUnannotated\t-\n")
                    cf.flush()
                    
                time.sleep(sleep_time)

        return ncbi_dict
                        
    def setup_output_directory(self):
        subdirs = ['bam', 'stat', 'summary', 'logs', 'index']
        if getattr(self.args, 'extract_reads', False): subdirs.append('reads')
        if self.args.consensus: subdirs.append('consensus')
        if getattr(self.args, 'call_variants', False): subdirs.append('variants')
            
        with Timer("创建输出目录结构"):
            for subdir in subdirs: (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)
            has_file_handler = any(isinstance(h, logging.FileHandler) for h in logger.handlers)
            if not has_file_handler:
                log_file = self.output_dir / 'logs' / f'analysis_global_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
                file_handler = logging.FileHandler(log_file, encoding='utf-8')
                file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
                logger.addHandler(file_handler)
        return self.output_dir
        
    def build_index(self):
        tool = self.args.tool
        ref_path = Path(self.args.reference).resolve()
        if tool in ['strobealign', 'minimap2']:
            self.index_path = str(ref_path)
            return
        prefix = self.output_dir / 'index' / f"{ref_path.stem}_{tool}"
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
            
        logger.info(f"🏗️ 构建全局 {tool} 索引...")
        with Timer(f"构建 {tool} 索引"):
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"❌ 索引构建失败: {result.stderr}")
                sys.exit(1)
    


    def find_samples(self):
        input_dir = Path(self.args.input_dir) if self.args.input_dir else None
        fq_exts = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
        samples = []
        if input_dir and self.args.single_end:
            for ext in fq_exts:
                for fastq in input_dir.glob(f'*{ext}'):
                    sname = fastq.name.replace(ext, '').replace('.gz', '')
                    for suf in ['_1', '_R1', '.R1', "R1_10239", '.1', '_unmapped', '.unmapped', '_trimmed', '.trimmed']:
                        if suf in sname: sname = sname.split(suf)[0]
                    samples.append({'name': sname, 'r1': str(fastq), 'r2': None})
        elif input_dir:
            patterns = [('*_1.f*q*', '*_2.f*q*'), ('*_R1*.f*q*', '*_R2*.f*q*'), ('*.R1.*', '*.R2.*'), 
                        ('*.1.f*q*', '*.2.f*q*'), ('*_1_*.f*q*', '*_2_*.f*q*'), 
                        ('*_unmapped.R1.fq.gz', '*_unmapped.R2.fq.gz'), ('*unmapped.R1.fq.gz', '*unmapped.R2.fq.gz'), 
                        ('*_1_unmapped.*', '*_2_unmapped.*'), ('*.unmapped.R1_10239.*', '*.unmapped.R2_10239.*')]
            found_files = set()
            for p1, p2 in patterns:
                for r1_file in input_dir.glob(p1):
                    if r1_file in found_files: continue
                    r1_name, r2_name = r1_file.name, None
                    if '_1.' in r1_name: r2_name = r1_name.replace('_1.', '_2.')
                    elif '_R1' in r1_name: r2_name = r1_name.replace('_R1', '_R2')
                    elif '.R1.' in r1_name: r2_name = r1_name.replace('.R1.', '.R2.')
                    elif '.1.' in r1_name: r2_name = r1_name.replace('.1.', '.2.')
                    elif '_1_' in r1_name: r2_name = r1_name.replace('_1_', '_2_')
                    elif '.unmapped.R1.' in r1_name: r2_name = r1_name.replace('.unmapped.R1.', '.unmapped.R2.')
                    elif 'unmapped.R1.' in r1_name: r2_name = r1_name.replace('unmapped.R1.', 'unmapped.R2.')
                    elif '_1_unmapped.' in r1_name: r2_name = r1_name.replace('_1_unmapped.', '_2_unmapped.')
                    elif '.unmapped.R1_10239.' in r1_name: r2_name = r1_name.replace('.unmapped.R1_10239.', '.unmapped.R2_10239.')
                    
                    if r2_name and (input_dir / r2_name).exists():
                        sname = r1_name
                        for suf in ['_1', '_R1', '.R1', '.1', '_unmapped', '.unmapped', '_1_unmapped']:
                            if suf in sname: sname = sname.split(suf)[0]
                        
                        # 【修复点 1】：去除了 sname = sname.split('.')[0]，保留文件名中的小数点信息
                        samples.append({'name': sname, 'r1': str(r1_file), 'r2': str(input_dir / r2_name)})
                        found_files.update([r1_file, input_dir / r2_name])
            if not samples:
                all_fqs = []
                for ext in fq_exts: all_fqs.extend(input_dir.glob(f'*{ext}'))
                fg = defaultdict(list)
                for f in all_fqs:
                    prefix = re.sub(r'[._](R?[12]|unmapped\.R?[12])[._].*', '', f.name)
                    
                    # 【修复点 2】：同样去除了 fallback 机制中的 if '.' in prefix: prefix = prefix.split('.')[0]
                    fg[prefix].append(str(f))
                for prefix, files in fg.items():
                    if len(files) == 2:
                        r1, r2 = None, None
                        for f in files:
                            if any(p in f for p in ['_1', '_R1', '.R1.', '.1.', 'unmapped.R1']): r1 = f
                            elif any(p in f for p in ['_2', '_R2', '.R2.', '.2.', 'unmapped.R2']): r2 = f
                        if r1 and r2: samples.append({'name': prefix, 'r1': r1, 'r2': r2})
        if self.args.sample_list:
            with open(self.args.sample_list, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2 and not line.startswith('#') and Path(parts[1]).exists():
                        samples.append({'name': parts[0], 'r1': parts[1], 'r2': parts[2] if len(parts)>2 and Path(parts[2]).exists() else None})
        
        unique_samples, seen = [], set()
        for s in samples:
            if s['name'] not in seen: unique_samples.append(s); seen.add(s['name'])
        self.samples = unique_samples
        if not self.samples:
            logger.error("❌ 未找到任何样本文件。")
            sys.exit(1)
        logger.info(f"✅ 找到 {len(self.samples)} 个测序样本")
    
    def align_sample(self, sample):
        sample_name = sample['name']
        bam_dir = self.output_dir / 'bam'
        sorted_bam = bam_dir / f'{sample_name}.sorted.bam'
        err_log = self.output_dir / 'logs' / f'{sample_name}.align.err'
        
        # [安全增强]: 仅当文件存在且大于 1KB 时，断点续跑才生效，无视 0 字节的破损文件
        if sorted_bam.exists() and sorted_bam.stat().st_size > 1024 and self.args.skip_alignment:
            self.write_sample_log(sample_name, f"⏭️[断点续跑] 跳过比对，直接使用已有BAM: {sorted_bam}")
            return sorted_bam
            
        self.write_sample_log(sample_name, f"\n--- 1. 序列比对 (使用 {self.args.tool}) ---")
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
        
        view_cmd = ['samtools', 'view', '-u', '-F', '2304', '-']
        sort_cmd = ['samtools', 'sort', '-@', str(min(4, self.args.threads)), '-o', str(sorted_bam)]
        
        # [逻辑修复]: 使用 shlex 和 bash pipefail 消除 Python 并发死锁与文件描述符崩溃
        import shlex
        align_cmd_str = ' '.join(shlex.quote(str(x)) for x in align_cmd)
        view_cmd_str = ' '.join(shlex.quote(str(x)) for x in view_cmd)
        sort_cmd_str = ' '.join(shlex.quote(str(x)) for x in sort_cmd)
        
        # 使用 set -o pipefail 确保任何一环报错都能截停，错误流安全落盘防死锁
        full_cmd = f"set -o pipefail; ({align_cmd_str} | {view_cmd_str} | {sort_cmd_str}) 2> '{err_log}'"
        
        self.write_sample_log(sample_name, f"Pipeline CMD: {full_cmd}")
        
        proc = subprocess.run(full_cmd, shell=True, executable='/bin/bash')
        
        if proc.returncode != 0:
            err_msg = ""
            if err_log.exists():
                with open(err_log, 'r') as f:
                    err_msg = f.read()[-2000:]
            self.write_sample_log(sample_name, f"❌ 比对与排序失败:\n{err_msg}")
            if sorted_bam.exists(): sorted_bam.unlink()
            return None
            
        self.write_sample_log(sample_name, f"✅ 比对与深度净化(去除次优/嵌合)成功。详细日志见: {err_log.name}")
        
        index_cmd = ['samtools', 'index', str(sorted_bam)]
        subprocess.run(index_cmd, capture_output=True, text=True)
        return sorted_bam

    def get_global_bam_qc(self, bam_file, sample_name):
        qc_file = self.output_dir / 'stat' / f"{sample_name}.global_qc.txt"
        
        if self.args.resume and qc_file.exists():
            self.write_sample_log(sample_name, f"⏭️  [断点续跑] 跳过全局 BAM 统计，直接使用已有文件: {qc_file.name}")
            return qc_file
            
        self.write_sample_log(sample_name, "\n--- 2.1 全局比对质量诊断 (samtools stats) ---")
        cmd = f"samtools stats '{bam_file}' | grep ^SN | cut -f 2- > '{qc_file}'"
        
        try:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            with open(qc_file, 'r') as f:
                lines = f.readlines()
                total, mapped, mq0, unmapped = 0, 0, 0, 0
                for line in lines:
                    if 'raw total sequences' in line: total = line.split('\t')[1].strip()
                    elif line.startswith('reads mapped:'): mapped = line.split('\t')[1].strip()
                    elif line.startswith('reads unmapped:'): unmapped = line.split('\t')[1].strip()
                    elif line.startswith('reads MQ0:'): mq0 = line.split('\t')[1].strip()
                
                log_msg = (
                    f"✅ 全局质控完成:\n"
                    f"  - 总测序 Reads  : {total}\n"
                    f"  - 未比对 (Unmapped) : {unmapped}\n"
                    f"  - 总比对 (Mapped)   : {mapped}\n"
                    f"  - 其中多重比对 (MQ0): {mq0} (占总比对的 {float(mq0)/float(mapped)*100:.1f}% 如果 Mapped>0)\n"
                ) if int(mapped) > 0 else f"✅ 全局质控完成: 0 条 reads 成功比对。"
                self.write_sample_log(sample_name, log_msg)
            return qc_file
        except Exception as e:
            self.write_sample_log(sample_name, f"❌ 全局 BAM 统计失败:\n{e}")
            return None

    def get_idxstats(self, bam_file, sample_name):
        self.write_sample_log(sample_name, "\n--- 2.2 提取测序绝对文库大小 (samtools idxstats) ---")
        cmd = ['samtools', 'idxstats', str(bam_file)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self.write_sample_log(sample_name, f"❌ idxstats 失败:\n{result.stderr}")
            return None
            
        stats, total_reads, total_mapped_reads, total_rpk = {}, 0, 0, 0.0
        for line in result.stdout.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 4:
                ref, length, mapped, unmapped = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
                stats[ref] = mapped
                total_reads += (mapped + unmapped)
                total_mapped_reads += mapped
                if length > 0 and mapped > 0: total_rpk += (mapped * 1000.0) / length
                
        self.write_sample_log(sample_name, f"Total Library Reads: {total_reads} | Total Mapped: {total_mapped_reads} | Total RPK: {total_rpk:.2f}")
        return {'stats': stats, 'global_total_reads': total_reads, 'total_mapped_reads': total_mapped_reads, 'total_rpk': total_rpk}

    def get_unique_multi_stats(self, bam_file, sample_name):
        stat_out = self.output_dir / 'stat' / f"{sample_name}.um_stat.tsv"
        
        if self.args.resume and stat_out.exists():
            self.write_sample_log(sample_name, f"⏭️  [断点续跑] 跳过 Unique/Multi 解析，直接使用已有文件: {stat_out.name}")
            stats = {}
            with open(stat_out, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 3:
                        stats[parts[0]] = {'unique': int(parts[1]), 'multi': int(parts[2])}
            return stats
            
        self.write_sample_log(sample_name, "\n--- 2.5 提取唯一/多重比对精细统计 (Unique/Multi-mapped) ---")
        
        # AWK 逻辑终极优化 (兼容所有比对软件，MAPQ < 10)
        awk_script = """
        {
            ref = $3;
            if (ref == "*") next; 
            
            if ($5 < 10) {
                multi[ref]++
            } else {
                uniq[ref]++
            }
        }
        END {
            for (r in uniq) {
                m = multi[r]+0;
                print r "\\t" uniq[r] "\\t" m;
                delete multi[r];
            }
            for (r in multi) {
                print r "\\t0\\t" multi[r];
            }
        }
        """
        
        cmd = f"samtools view -F 2308 '{bam_file}' | awk '{awk_script}' > '{stat_out}'"
        
        try:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            self.write_sample_log(sample_name, f"✅ Unique/Multi 统计完成 (基于 MAPQ 标准)。")
        except subprocess.CalledProcessError as e:
            self.write_sample_log(sample_name, f"❌ Unique/Multi 统计失败:\n{e}")
            return {}

        stats = {}
        with open(stat_out, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 3:
                    stats[parts[0]] = {'unique': int(parts[1]), 'multi': int(parts[2])}
        return stats

    def run_pandepth(self, bam_file, sample_name):
        output_prefix = self.output_dir / 'stat' / sample_name
        stat_file = output_prefix.with_suffix('.chr.stat.gz')
        if stat_file.exists() and self.args.skip_depth:
            self.write_sample_log(sample_name, f"⏭️  [断点续跑] 跳过深度计算，使用已有文件: {stat_file}")
            return stat_file
            
        self.write_sample_log(sample_name, "\n--- 3. 深度计算 (pandepth) ---")
        pandepth_cmd = ['pandepth', '-a', '-i', str(bam_file), '-o', str(output_prefix), '-t', str(min(10, self.args.threads))]
        result = subprocess.run(pandepth_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self.write_sample_log(sample_name, f"❌ pandepth 失败:\n{result.stderr}")
            return None
        self.write_sample_log(sample_name, "✅ pandepth 运行成功。")
        return stat_file
    
    def parse_stat_file(self, stat_file, sample_name, idx_data):
        if not idx_data: return []
        global_total_reads = idx_data['global_total_reads']
        total_mapped_reads = idx_data.get('total_mapped_reads', 0)
        total_rpk = idx_data['total_rpk']
        idx_stats = idx_data['stats']
        um_stats = idx_data.get('um_stats', {}) 
        
        results = []
        try:
            with gzip.open(stat_file, 'rt') as f:
                for line in f:
                    if not line or line.startswith('#'): continue
                    parts = line.split()
                    if len(parts) < 6: continue
                    virus_name, length, covered_site, total_depth, coverage, mean_depth = parts
                    
                    length_val, covered_site_val, total_depth_val = int(length), int(covered_site), int(total_depth)
                    coverage_val, mean_depth_val = float(coverage), float(mean_depth)
                    mapped_reads = idx_stats.get(virus_name, 0)
                    
                    if mapped_reads == 0: 
                        continue 
                    
                    # 相对丰度计算 (占全部 Mapped Reads 百分比)
                    rel_abundance = (mapped_reads / total_mapped_reads * 100.0) if total_mapped_reads > 0 else 0.0
                    
                    # Unique 和 Multi 读取
                    virus_um = um_stats.get(virus_name, {'unique': 0, 'multi': 0})
                    unique_reads = virus_um['unique']
                    multi_reads = virus_um['multi']
                    
                    # [新增]: 计算 Unique 占比
                    unique_pct = (unique_reads / mapped_reads * 100.0) if mapped_reads > 0 else 0.0
                    
                    rpm = (mapped_reads * 1e6) / global_total_reads if global_total_reads > 0 else 0.0
                    fpkm = (mapped_reads * 1e9) / (global_total_reads * length_val) if global_total_reads > 0 and length_val > 0 else 0.0
                    tpm = ((mapped_reads * 1000.0) / length_val * 1e6) / total_rpk if total_rpk > 0 and length_val > 0 else 0.0
                        
                    results.append({
                        'Sample': sample_name, 'Virus': virus_name, 'Length': length_val,
                        'CoveredSite': covered_site_val, 'TotalDepth': total_depth_val,
                        'Coverage(%)': coverage_val, 'MeanDepth': mean_depth_val,
                        'MappedReads': mapped_reads,
                        'Unique_Reads': unique_reads,  
                        'Multi_Reads': multi_reads,
                        'Unique(%)': round(unique_pct, 2), # 新增的 Unique_Reads/MappedReads 占比
                        'Relative_Abundance(%)': round(rel_abundance, 4), 
                        'RPM': round(rpm, 2),
                        'FPKM': round(fpkm, 2), 'TPM': round(tpm, 2)
                    })
            self.write_sample_log(sample_name, f"✅ 数据解析完毕, 提取出 {len(results)} 条含比对记录的 Raw 数据。")
            return results
        except Exception as e:
            self.write_sample_log(sample_name, f"❌ 统计文件解析异常: {str(e)}")
            return []
    
    def process_sample(self, sample):
        sample_name = sample['name']
        start_time = datetime.now()
        
        log_mode = 'a' if self.args.resume else 'w'
        status_str = "[断点续跑]" if self.args.resume else "[全新运行]"
        self.write_sample_log(sample_name, f"\n{'='*50}\n[{start_time.strftime('%Y-%m-%d %H:%M:%S')}] {status_str} 开始处理测序样本: {sample_name}", mode=log_mode)
        
        try:
            bam_file = self.align_sample(sample)
            if not bam_file: return None
            
            self.get_global_bam_qc(bam_file, sample_name)
            
            idx_data = self.get_idxstats(bam_file, sample_name)
            if not idx_data: return None
            
            um_stats = self.get_unique_multi_stats(bam_file, sample_name)
            if um_stats: idx_data['um_stats'] = um_stats
            
            stat_file = self.run_pandepth(bam_file, sample_name)
            if not stat_file: return None
            
            res = self.parse_stat_file(stat_file, sample_name, idx_data)
            
            end_time = datetime.now()
            duration = end_time - start_time
            self.write_sample_log(sample_name, f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] 🎉 样本整体处理流程结束[耗时: {Timer.format_duration(duration.total_seconds())}]")
            return res
        except Exception as e:
            self.write_sample_log(sample_name, f"❌ 样本执行期间发生致命错误: {str(e)}")
            return None
    
    def run_pipeline(self):
        start_time = datetime.now()
        logger.info("=" * 60)
        logger.info(f"🚀 开始病毒分析管线 (工具: {self.args.tool} | 断点续跑: {self.args.resume})")
        logger.info("=" * 60)
        
        self.setup_output_directory()
        self.build_index()
        self.find_samples()
        
        all_results = []
        if self.args.parallel and self.args.threads > 1:
            max_parallel_jobs = min(max(1, self.args.threads // self.args.align_threads), len(self.samples))
            if self.args.parallel_jobs: max_parallel_jobs = min(self.args.parallel_jobs, len(self.samples))
            logger.info(f"🔀 开始并行处理任务 (最大并行数: {max_parallel_jobs})")
            
            with ProcessPoolExecutor(max_workers=max_parallel_jobs) as executor:
                future_to_sample = {executor.submit(self.process_sample, s): s for s in self.samples}
                with tqdm(total=len(self.samples), desc="处理进度", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(future_to_sample):
                        sample = future_to_sample[future]
                        pbar.set_postfix_str(f"处理: {sample['name']} 中")
                        res = future.result()
                        if res is not None: all_results.extend(res)
                        pbar.update(1)
        else:
            with tqdm(self.samples, desc="处理进度", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                for sample in pbar:
                    pbar.set_postfix_str(f"处理: {sample['name']} 中")
                    res = self.process_sample(sample)
                    if res is not None: all_results.extend(res)
        
        best_summary_file = None
        if all_results:
            best_records, best_summary_file = self.save_results(all_results)
            
            if getattr(self.args, 'extract_reads', False) and best_records:
                self.extract_mapped_reads(best_records)
                
            if getattr(self.args, 'call_variants', False) and best_records:
                self.run_variants_calling(best_records)
                
            if getattr(self.args, 'consensus', False) and best_records:
                self.build_consensus_sequences(best_records)
                
        else:
            logger.warning("⚠️ 未找到任何有 MappedReads 的病毒记录")
        
        total_duration = datetime.now() - start_time
        logger.info("=" * 60)
        logger.info("✨ 分析管线全流程运行结束")
        logger.info(f"⏱️  管线总耗时: {Timer.format_duration(total_duration.total_seconds())}")
        
        if best_summary_file:
            self.generate_report_txt(best_summary_file)
        
        logger.info("=" * 60)

    def save_results(self, all_results):
        with Timer("汇总去重与保存结果"):
            df_all = pd.DataFrame(all_results)
            
            if 'Virus' in df_all.columns:
                found_viruses = df_all['Virus'].unique().tolist()
                df_all['description'] = "Unannotated"
                df_all['taxonomy'] = "Unannotated"
                df_all['taxonomy_id'] = "-"
                
                if self.args.rvdb_taxon:
                    taxon_db = self.fast_load_rvdb_taxonomy(found_viruses)
                    if taxon_db:
                        def get_tax_info(vid, key):
                            if vid in taxon_db: return taxon_db[vid][key]
                            base_id = vid.split('.')[0]
                            if base_id in taxon_db: return taxon_db[base_id][key]
                            return "Unannotated" if key != "taxonomy_id" else "-"
                            
                        df_all['description'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'description'))
                        df_all['taxonomy'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy'))
                        df_all['taxonomy_id'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy_id'))

                unannotated_viruses = df_all[df_all['taxonomy'] == 'Unannotated']['Virus'].unique().tolist()
                if unannotated_viruses:
                    if self.args.email:
                        ncbi_db = self.fetch_ncbi_taxonomy(unannotated_viruses)
                        if ncbi_db:
                            for vid in unannotated_viruses:
                                if vid in ncbi_db:
                                    mask = df_all['Virus'] == vid
                                    df_all.loc[mask, 'description'] = ncbi_db[vid]['description']
                                    df_all.loc[mask, 'taxonomy'] = ncbi_db[vid]['taxonomy']
                                    df_all.loc[mask, 'taxonomy_id'] = ncbi_db[vid]['taxonomy_id']

                cols = list(df_all.columns)
                for col in ['taxonomy_id', 'taxonomy', 'description']:
                    if col in cols: cols.remove(col)
                    cols.insert(cols.index('Virus') + 1, col)
                df_all = df_all[cols]

            summary_dir = self.output_dir / 'summary'
            
            raw_summary_file = summary_dir / f'all_viruses.raw.{self.args.format}'
            if self.args.format == 'csv': df_all.to_csv(raw_summary_file, index=False)
            else: df_all.to_csv(raw_summary_file, sep='\t', index=False)
            self.total_raw_records = len(df_all)
            
            df_filtered = df_all[
                (df_all['Coverage(%)'] > self.args.coverage) &
                (df_all['MeanDepth'] > self.args.meandepth) &
                (df_all['FPKM'] >= self.args.min_fpkm) &
                (df_all['TPM'] >= self.args.min_tpm)
            ].copy()
            
            filtered_summary_file = summary_dir / f'all_viruses.summary.{self.args.format}'
            if self.args.format == 'csv': df_filtered.to_csv(filtered_summary_file, index=False)
            else: df_filtered.to_csv(filtered_summary_file, sep='\t', index=False)
            self.total_filtered_records = len(df_filtered)

            if 'taxonomy' in df_filtered.columns:
                df_filtered['group_tax'] = df_filtered.apply(lambda x: x['taxonomy'] if x['taxonomy'] not in ['Unannotated', '-'] else x['Virus'], axis=1)
            else:
                df_filtered['group_tax'] = df_filtered['Virus']
                
            df_sorted = df_filtered.sort_values(by=['Sample', 'group_tax', 'Coverage(%)', 'MeanDepth'], ascending=[True, True, False, False])
            df_best = df_sorted.drop_duplicates(subset=['Sample', 'group_tax'], keep='first').copy()
            
            df_best.drop(columns=['group_tax'], inplace=True)

            best_summary_file = summary_dir / f'all_viruses.best.summary.{self.args.format}'
            if self.args.format == 'csv': df_best.to_csv(best_summary_file, index=False)
            else: df_best.to_csv(best_summary_file, sep='\t', index=False)
            self.df_best = df_best
            
            logger.info(f"✅ 全量 Mapped 记录已保存至 all_viruses.raw ({self.total_raw_records} 条)")
            logger.info(f"✅ 阈值过滤记录已保存至 all_viruses.summary ({self.total_filtered_records} 条)")
            logger.info(f"✅ 最优代表株结果已保存至 all_viruses.best.summary ({len(df_best)} 条核心记录)")
            
            return df_best.to_dict('records'), best_summary_file

    def extract_fastas_with_python(self, target_ids, virus_to_tax, out_dir):
        target_set = set(target_ids)
        found = set()
        current_id, current_seq = None, []
        
        def save_record(vid, seq_list):
            if vid in target_set:
                vname = re.sub(r'[\\/*?:"<>| ]', "_", vid)
                tax = virus_to_tax.get(vid, 'Unannotated')
                safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", tax)
                if len(safe_tax) > 50: safe_tax = safe_tax[:50]
                
                folder_name = f"{safe_tax}_{vname}"
                v_dir = out_dir / folder_name
                v_dir.mkdir(parents=True, exist_ok=True)
                
                with open(v_dir / f"{folder_name}.ref.fasta", 'w') as outf:
                    outf.write(f">{vid}\n" + "".join(seq_list) + "\n")
                found.add(vid)
                
        with open(self.args.reference, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id: save_record(current_id, current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id: save_record(current_id, current_seq)

    def _run_single_read_extraction(self, sample):
        bam_file = self.output_dir / 'bam' / f"{sample}.sorted.bam"
        reads_dir = self.output_dir / 'reads'
        
        r1 = reads_dir / f"{sample}_virus_1.fastq.gz"
        r2 = reads_dir / f"{sample}_virus_2.fastq.gz"
        rs = reads_dir / f"{sample}_virus_single.fastq.gz"
        
        if self.args.resume and (r1.exists() or r2.exists() or rs.exists()):
            self.write_sample_log(sample, f"\n⏭️ [断点续跑] 跳过 reads 提取，目标文件已存在: {reads_dir.name}/")
            return True
            
        if not bam_file.exists(): return False

        log_block = [f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] === 提取比对上病毒的序列 Reads ==="]
        try:
            cmd = (
                f"samtools fastq -F 4 "
                f"-1 >(pigz -c > '{r1}') "
                f"-2 >(pigz -c > '{r2}') "
                f"'{bam_file}' | pigz -c > '{rs}'"
            )
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash', stderr=subprocess.PIPE)
            
            for f in [r1, r2, rs]:
                if f.exists() and f.stat().st_size <= 50:
                    f.unlink()
                    
            log_block.append(f"✅ {sample} 病毒 reads 提取成功! 存入目录: reads/")
            self.write_sample_log(sample, "\n".join(log_block))
            return True
        except subprocess.CalledProcessError as e:
            err_msg = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or str(e))
            log_block.append(f"❌ reads 提取失败:\n{err_msg}")
            self.write_sample_log(sample, "\n".join(log_block))
            return False

    def extract_mapped_reads(self, best_results):
        logger.info("\n📦 开始并行提取比对上病毒的序列 Reads (仅限阳性样本)...")
        reads_dir = self.output_dir / 'reads'
        reads_dir.mkdir(parents=True, exist_ok=True)
        
        unique_samples = list(set(r['Sample'] for r in best_results))
        success_count = 0
        
        with Timer("阳性病毒 Reads 提取"):
            with ProcessPoolExecutor(max_workers=min(len(unique_samples), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_read_extraction, s): s for s in unique_samples}
                with tqdm(total=len(unique_samples), desc="Reads提取", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(futures):
                        sample = futures[future]
                        pbar.set_postfix_str(f"{sample} 中")
                        if future.result(): success_count += 1
                        pbar.update(1)
            
            logger.info(f"✅ 成功提取了 {success_count} 个阳性样本的病毒序列 Reads")

    def _run_single_variant_calling(self, task):
        sample, virus, taxonomy = task
        variants_dir = self.output_dir / 'variants'
        
        vname = re.sub(r'[\\/*?:"<>| ]', "_", virus)
        safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", taxonomy)
        if len(safe_tax) > 50: safe_tax = safe_tax[:50]
        
        folder_name = f"{safe_tax}_{vname}"
        virus_dir = variants_dir / folder_name
        
        bam_file = self.output_dir / 'bam' / f"{sample}.sorted.bam"
        ref_fasta = virus_dir / f"{folder_name}.ref.fasta"
        fixed_bam = virus_dir / f"{sample}.{folder_name}.fixed.bam"
        
        caller = self.args.variant_caller
        out_ext = 'tsv' if caller == 'ivar' else 'vcf'
        out_file = virus_dir / f"{sample}.{folder_name}.variants.{out_ext}"
        
        if self.args.resume and out_file.exists() and out_file.stat().st_size > 0:
            self.write_sample_log(sample, f"\n⏭️ [断点续跑] 跳过变异检测，目标文件已存在: {out_file.name}")
            return True
            
        if not ref_fasta.exists() or not bam_file.exists(): return False
            
        log_block = [f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] === 变异检测 ({caller}): {taxonomy} ({virus}) ==="]
        
        try:
            if not Path(str(ref_fasta) + '.fai').exists():
                subprocess.run(f"samtools faidx '{ref_fasta}'", shell=True, check=True, stderr=subprocess.PIPE)
                
            awk_cmd = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
            pipe_cmd = f"samtools view -h '{bam_file}' '{virus}' | {awk_cmd} | samtools view -b -o '{fixed_bam}'"
            subprocess.run(pipe_cmd, shell=True, check=True, stderr=subprocess.PIPE)
            subprocess.run(f"samtools index '{fixed_bam}'", shell=True, check=True, stderr=subprocess.PIPE)
            
            if caller == 'freebayes':
                cmd = f"freebayes -p 1 -f '{ref_fasta}' '{fixed_bam}' > '{out_file}'"
                subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE)
                
            elif caller == 'ivar':
                pileup = virus_dir / f"{sample}.pileup.txt"
                prefix = str(out_file).replace('.tsv', '')
                cmd1 = f"samtools mpileup -A -aa -d 0 -Q 0 --reference '{ref_fasta}' '{fixed_bam}' > '{pileup}'"
                cmd2 = f"cat '{pileup}' | ivar variants -r '{ref_fasta}' -p '{prefix}'"
                subprocess.run(cmd1, shell=True, check=True, stderr=subprocess.PIPE)
                subprocess.run(cmd2, shell=True, check=True, stderr=subprocess.PIPE)
                if pileup.exists(): pileup.unlink()
                
            elif caller == 'lofreq':
                tmp_bam = virus_dir / f"{sample}.tmp.lofreq.bam"
                cmd1 = f"lofreq indelqual -f '{ref_fasta}' --dindel -o '{tmp_bam}' '{fixed_bam}'"
                cmd2 = f"lofreq call -f '{ref_fasta}' --call-indels -o '{out_file}' '{tmp_bam}'"
                subprocess.run(cmd1, shell=True, check=True, stderr=subprocess.PIPE)
                subprocess.run(f"samtools index '{tmp_bam}'", shell=True, check=True, stderr=subprocess.PIPE)
                subprocess.run(cmd2, shell=True, check=True, stderr=subprocess.PIPE)
                if tmp_bam.exists(): tmp_bam.unlink()
                if Path(str(tmp_bam) + '.bai').exists(): Path(str(tmp_bam) + '.bai').unlink()
                
            log_block.append(f"✅ {sample} 变异检测 ({caller}) 成功! 存入目录: {folder_name}/")
            
            if fixed_bam.exists(): fixed_bam.unlink()
            if Path(str(fixed_bam) + '.bai').exists(): Path(str(fixed_bam) + '.bai').unlink()
            
            self.write_sample_log(sample, "\n".join(log_block))
            return True
        except subprocess.CalledProcessError as e:
            err_msg = e.stderr if hasattr(e, 'stderr') and e.stderr else str(e)
            log_block.append(f"❌ 变异检测失败:\n{err_msg}")
            self.write_sample_log(sample, "\n".join(log_block))
            return False

    def run_variants_calling(self, best_results):
        logger.info(f"\n🔬 开始并行执行变异检测 (工具: {self.args.variant_caller}, 按 Taxonomy 归档)...")
        variants_dir = self.output_dir / 'variants'
        
        virus_to_tax = {r['Virus']: r.get('taxonomy', 'Unannotated') for r in best_results}
        
        with Timer("核心代表株变异检测"):
            unique_viruses = list(set(r['Virus'] for r in best_results))
            self.extract_fastas_with_python(unique_viruses, virus_to_tax, variants_dir)
            
            tasks = [(r['Sample'], r['Virus'], r.get('taxonomy', 'Unannotated')) for r in best_results]
            success_count = 0
            
            with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_variant_calling, t): t for t in tasks}
                with tqdm(total=len(tasks), desc=f"变异检测", unit="任务", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(futures):
                        task = futures[future]
                        short_tax = str(task[2])[:10]
                        pbar.set_postfix_str(f"{task[0]} ({short_tax}..)中")
                        if future.result(): success_count += 1
                        pbar.update(1)
            
            logger.info(f"✅ 成功提取了 {success_count} 个按分类学命名的核心病毒变异检测结果 (VCF/TSV)")

    def _build_single_virus_index(self, ref_fasta, prefix):
        tool = self.args.tool
        if tool in ['strobealign', 'minimap2']: return True
        if tool == 'bwa' and not Path(f"{prefix}.bwt").exists():
            subprocess.run(['bwa', 'index', '-p', str(prefix), str(ref_fasta)], capture_output=True)
        elif tool == 'bwa-mem2' and not (Path(f"{prefix}.bwt.2bit.64").exists() or Path(f"{prefix}.0123").exists()):
            subprocess.run(['bwa-mem2', 'index', '-p', str(prefix), str(ref_fasta)], capture_output=True)
        elif tool == 'bowtie2' and not (Path(f"{prefix}.1.bt2").exists() or Path(f"{prefix}.1.bt2l").exists()):
            subprocess.run(['bowtie2-build', str(ref_fasta), str(prefix)], capture_output=True)
        elif tool == 'hisat2' and not (Path(f"{prefix}.1.ht2").exists() or Path(f"{prefix}.1.ht2l").exists()):
            subprocess.run(['hisat2-build', str(ref_fasta), str(prefix)], capture_output=True)
        return True

    def _run_single_consensus(self, task):
        sample_dict, virus, mean_depth, taxonomy = task
        sample = sample_dict['name']
        consensus_dir = self.output_dir / 'consensus'
        
        vname = re.sub(r'[\\/*?:"<>| ]', "_", virus)
        safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", taxonomy)
        if len(safe_tax) > 50: safe_tax = safe_tax[:50]
        
        folder_name = f"{safe_tax}_{vname}"
        virus_dir = consensus_dir / folder_name
        
        ref_fasta = virus_dir / f"{folder_name}.ref.fasta"
        out_fasta = virus_dir / f"{sample}.{folder_name}.consensus.fasta"
        
        if self.args.resume and out_fasta.exists() and out_fasta.stat().st_size > 0:
            self.write_sample_log(sample, f"\n⏭️ [断点续跑] 跳过共识提取，目标文件已存在: {out_fasta.name}")
            return True
        if not ref_fasta.exists(): return False
            
        min_depth = int(math.floor(mean_depth))
        if min_depth > 10: min_depth = 10
        elif min_depth < 1: min_depth = 1
        
        mode_str = "[重比对提取解剖模式]" if self.args.remap else "[全库提取模式]"
        tool_str = f"使用 {self.args.consensus_tool}"
        log_block = [f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] === 共识生成 {mode_str} ({tool_str}): {taxonomy} ({virus}) ==="]
        
        fixed_bam = virus_dir / f"{sample}.{folder_name}.fixed.bam"
        try:
            if self.args.remap:
                tool = self.args.tool
                index_path = str(virus_dir / f"{folder_name}.index") if tool not in ['strobealign', 'minimap2'] else str(ref_fasta)
                r1, r2 = sample_dict['r1'], sample_dict['r2']
                
                align_cmd = []
                if tool == 'strobealign':
                    align_cmd = ['strobealign', '-t', '2', index_path, r1]
                    if r2: align_cmd.append(r2)
                elif tool == 'minimap2':
                    preset = 'sr' if r2 else 'map-ont'
                    align_cmd = ['minimap2', '-ax', preset, '-t', '2', index_path, r1]
                    if r2: align_cmd.append(r2)
                elif tool == 'bwa':
                    align_cmd = ['bwa', 'mem', '-v', '1', '-t', '2', index_path, r1]
                    if r2: align_cmd.append(r2)
                elif tool == 'bwa-mem2':
                    align_cmd = ['bwa-mem2', 'mem', '-v', '1', '-t', '2', index_path, r1]
                    if r2: align_cmd.append(r2)
                elif tool == 'bowtie2':
                    align_cmd = ['bowtie2', '-p', '2', '-x', index_path]
                    if r2: align_cmd.extend(['-1', r1, '-2', r2])
                    else: align_cmd.extend(['-U', r1])
                elif tool == 'hisat2':
                    align_cmd = ['hisat2', '-p', '2', '-x', index_path]
                    if r2: align_cmd.extend(['-1', r1, '-2', r2])
                    else: align_cmd.extend(['-U', r1])
                
                view_cmd = ['samtools', 'view', '-u', '-F', '4', '-']
                sort_cmd = ['samtools', 'sort', '-@', '2', '-o', str(fixed_bam)]
                
                align_proc = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                view_proc = subprocess.Popen(view_cmd, stdin=align_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                sort_proc = subprocess.Popen(sort_cmd, stdin=view_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                align_proc.stdout.close()
                view_proc.stdout.close()
                sort_proc.communicate()
            else:
                global_bam_file = self.output_dir / 'bam' / f"{sample}.sorted.bam"
                if not global_bam_file.exists(): return False
                awk_cmd = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
                pipe_cmd = f"samtools view -h '{global_bam_file}' '{virus}' | {awk_cmd} | samtools view -b -o '{fixed_bam}'"
                subprocess.run(pipe_cmd, shell=True, check=True, stderr=subprocess.PIPE, text=True)

            subprocess.run(['samtools', 'index', str(fixed_bam)], check=True, capture_output=True)

            if self.args.consensus_tool == 'viral_consensus':
                vc_cmd = f"viral_consensus -i '{fixed_bam}' -r '{ref_fasta}' -o '{out_fasta}' --min_depth {min_depth}"
                subprocess.run(vc_cmd, shell=True, check=True, capture_output=True, text=True)
                
            elif self.args.consensus_tool == 'ivar':
                prefix = str(out_fasta).replace('.fasta', '')
                ivar_freq = self.args.ivar_freq
                mpileup_cmd = f"samtools mpileup -A -aa -d 600000 -B -Q 0 --reference '{ref_fasta}' '{fixed_bam}'"
                ivar_cmd = f"ivar consensus -p '{prefix}' -t {ivar_freq} -m {min_depth}"
                
                pipe_cmd = f"{mpileup_cmd} | {ivar_cmd}"
                subprocess.run(pipe_cmd, shell=True, check=True, capture_output=True, text=True)
                
                if Path(f"{prefix}.fa").exists():
                    shutil.move(f"{prefix}.fa", out_fasta)

            log_block.append(f"✅ {sample} 共识序列构建成功! 存入目录: {folder_name}/")
            
            if fixed_bam.exists(): fixed_bam.unlink()
            if Path(str(fixed_bam) + '.bai').exists(): Path(str(fixed_bam) + '.bai').unlink()
            
            self.write_sample_log(sample, "\n".join(log_block))
            return True
        except subprocess.CalledProcessError as e:
            err_output = e.stderr if hasattr(e, 'stderr') and e.stderr else str(e)
            log_block.append(f"❌ 共识序列生成失败:\n{err_output}")
            self.write_sample_log(sample, "\n".join(log_block))
            return False

    def build_consensus_sequences(self, best_results):
        tool_str = "iVar" if self.args.consensus_tool == 'ivar' else "viral_consensus"
        mode_str = "(重比对模式)" if self.args.remap else "(全库提取模式)"
        logger.info(f"\n🧬 开始并行生成核心病毒共识序列 [{tool_str}] {mode_str}...")
        
        consensus_dir = self.output_dir / 'consensus'
        virus_to_tax = {r['Virus']: r.get('taxonomy', 'Unannotated') for r in best_results}
        
        with Timer("核心共识序列构建"):
            unique_viruses = list(set(r['Virus'] for r in best_results))
            self.extract_fastas_with_python(unique_viruses, virus_to_tax, consensus_dir)
            
            if self.args.remap:
                logger.info("🛠️ 正在为目标毒株预构建单株独立索引 (Remap)...")
                for vid in unique_viruses:
                    vname = re.sub(r'[\\/*?:"<>| ]', "_", vid)
                    tax = virus_to_tax.get(vid, 'Unannotated')
                    safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", tax)[:50]
                    folder_name = f"{safe_tax}_{vname}"
                    ref_fasta = consensus_dir / folder_name / f"{folder_name}.ref.fasta"
                    index_prefix = consensus_dir / folder_name / f"{folder_name}.index"
                    self._build_single_virus_index(ref_fasta, index_prefix)
            
            sample_map = {s['name']: s for s in self.samples}
            tasks = []
            for r in best_results:
                s_dict = sample_map.get(r['Sample'])
                if s_dict:
                    tasks.append((s_dict, r['Virus'], r['MeanDepth'], r.get('taxonomy', 'Unannotated')))
            
            success_count = 0
            with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_consensus, t): t for t in tasks}
                with tqdm(total=len(tasks), desc="构建进度", unit="文件", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(futures):
                        task = futures[future]
                        short_tax = str(task[3])[:10]
                        pbar.set_postfix_str(f"{task[0]['name']} ({short_tax}..)中")
                        if future.result(): success_count += 1
                        pbar.update(1)
            
            logger.info(f"✅ 成功归档了 {success_count} 个按分类学命名的核心共识序列文件")

    def generate_report_txt(self, best_summary_file):
        report_file = self.output_dir / 'summary' / 'analysis_report.txt'
        df_best = self.df_best
        total_samples = len(self.samples)  # 获取总样本数计算检出率
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"病毒丰度核心分析报告 (基于全测序文库的绝对定量 + 分类去冗余)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"比对工具: {self.args.tool}\n")
            f.write(f"参考基因组: {self.args.reference}\n")
            if self.args.rvdb_taxon: f.write(f"本地注释文件: {self.args.rvdb_taxon}\n")
            f.write(f"基础过滤条件: Coverage > {self.args.coverage}%, MeanDepth > {self.args.meandepth}\n")
            f.write(f"定量过滤要求: FPKM >= {self.args.min_fpkm}, TPM >= {self.args.min_tpm}\n")
            
            consensus_status = "否"
            if self.args.consensus:
                mode_info = "重比对解剖模式(消除同源竞争)" if self.args.remap else "全库提取模式"
                consensus_status = f"已针对核心代表株生成 (引擎: {self.args.consensus_tool}, 策略: {mode_info})"
            f.write(f"共识提取状态: {consensus_status}\n")
            f.write(f"变异检测状态: {'已使用 ' + self.args.variant_caller + ' 进行检测' if getattr(self.args, 'call_variants', False) else '否'}\n\n")
            
            f.write(f"总处理样本数: {total_samples}\n")
            f.write(f"[第一级] Raw 包含全量未过滤记录总数: {self.total_raw_records}\n")
            f.write(f"[第二级] Summary 阈值过滤后记录总数: {self.total_filtered_records}\n")
            f.write(f"[第三级] Best 同源去重后核心记录数: {len(df_best)}\n\n")
            
            if 'Virus' in df_best.columns:
                f.write("同源去重后的最终确诊核心病毒群落概览:\n")
                f.write("-" * 60 + "\n")
                virus_counts = df_best['Virus'].value_counts()
                for virus, count in virus_counts.items():
                    vd = df_best[df_best['Virus'] == virus]
                    desc = vd['description'].iloc[0] if 'description' in vd.columns else ''
                    
                    # [新增]: 检出率与各类均值计算
                    detection_rate = (count / total_samples) * 100.0
                    mean_tpm = vd['TPM'].mean()
                    mean_rel_abund = vd['Relative_Abundance(%)'].mean() if 'Relative_Abundance(%)' in vd.columns else 0
                    mean_unique = vd['Unique_Reads'].mean() if 'Unique_Reads' in vd.columns else 0
                    mean_multi = vd['Multi_Reads'].mean() if 'Multi_Reads' in vd.columns else 0
                    mean_unique_pct = vd['Unique(%)'].mean() if 'Unique(%)' in vd.columns else 0
                    
                    f.write(f"{virus} ({desc[:40]}...):\n")
                    f.write(f"  样本检出数: {count} / {total_samples} (群体检出率: {detection_rate:.1f}%)\n")
                    f.write(f"  核心平均 Unique_Reads: {mean_unique:.1f} (平均占比 {mean_unique_pct:.1f}%, 特异性高置信度)\n")
                    f.write(f"  核心平均 Multi_Reads:  {mean_multi:.1f} (跨物种多重同源比对)\n")
                    f.write(f"  核心平均相对丰度(%):   {mean_rel_abund:.4f}% (占整个建库测序成功比对序列的百分比)\n")
                    f.write(f"  核心平均 TPM:  {mean_tpm:.2f} (长度归一化后的绝对相对占比)\n\n")
        
def main():
    parser = argparse.ArgumentParser(
        description='病毒分析终极版 (带断点续跑、三级过滤、相对丰度、重比对解剖拼接、双共识引擎)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
【重磅增强功能】:
  --remap           : 开启重比对解剖模式。提取单个病毒进行专属建库与原始双端重比对，消除宏基因组多重序列竞争，极大提升 Consensus 质量。
  --consensus_tool  : 支持切换共识生成引擎，可选 viral_consensus 或临床金标准 ivar。
  --ivar_freq       : 仅在选择 ivar 时生效，控制杂合位点阈值 (默认 0.75)。
        """
    )
    
    input_group = parser.add_argument_group('输入输出参数')
    input_group.add_argument('--input_dir', type=str, help='输入文件夹路径（包含fastq文件）')
    input_group.add_argument('--reference', type=str, required=True, help='参考基因组fasta文件路径')
    input_group.add_argument('--rvdb_taxon', type=str, help='RVDB 分类注释文件路径 (如: RVDB_Taxon_Current.tab)')
    input_group.add_argument('--email', type=str, help='[可选] 提供有效邮箱，用于启用 NCBI 联网获取未知病毒注释')
    input_group.add_argument('--api_key', type=str, help='[可选] NCBI API Key，可以提高联网获取的访问频率')
    input_group.add_argument('--output_dir', type=str, default='./virus_analysis', help='输出目录路径 [默认: ./virus_analysis]')
    input_group.add_argument('--sample_list', type=str, help='样本列表文件（每行: 样本名 fastq1 fastq2）')
    input_group.add_argument('--single_end', action='store_true', help='使用单端测序数据')
    
    align_group = parser.add_argument_group('高级流程与变异控制参数')
    align_group.add_argument('--tool', type=str, choices=['strobealign', 'minimap2', 'bwa', 'bwa-mem2', 'bowtie2', 'hisat2'], default='strobealign', help='比对工具选择[默认: strobealign]')
    align_group.add_argument('--align_threads', type=int, default=8, help='每个样本比对的线程数[默认: 8]')
    align_group.add_argument('--extract_reads', action='store_true', help='提取比对上病毒的 reads，自动分离并输出为单双端 fastq.gz')
    
    # 共识模块与 Remap 参数
    align_group.add_argument('--consensus', action='store_true', help='基于去重后的最佳代表株自动生成带分类学归档的共识序列')
    align_group.add_argument('--consensus_tool', type=str, choices=['viral_consensus', 'ivar'], default='viral_consensus', help='共识序列构建工具选择 [默认: viral_consensus]')
    align_group.add_argument('--ivar_freq', type=float, default=0.75, help='iVar 判定主流等位基因的最小频率阈值 (默认 0.75，仅 consensus_tool=ivar 有效)')
    align_group.add_argument('--remap', action='store_true', help='开启重比对解剖模式。在提取 Consensus 时为单株建立独立索引进行无干扰重比对，组装效果翻倍')
    
    align_group.add_argument('--call_variants', action='store_true', help='基于去重后的最佳代表株开启变异检测')
    align_group.add_argument('--variant_caller', type=str, choices=['freebayes', 'ivar', 'lofreq'], default='freebayes', help='变异检测工具选择 [默认: freebayes]')
    align_group.add_argument('--resume', action='store_true', help='开启断点续跑模式，自动跳过已存在的耗时结果及请求')
    
    filter_group = parser.add_argument_group('过滤参数')
    filter_group.add_argument('--coverage', type=float, default=90.0, help='Coverage(%%)阈值[默认: 90]')
    filter_group.add_argument('--meandepth', type=float, default=10.0, help='MeanDepth阈值[默认: 10]')
    filter_group.add_argument('--min_fpkm', type=float, default=0.0, help='FPKM最小值过滤[默认: 0]')
    filter_group.add_argument('--min_tpm', type=float, default=0.0, help='TPM最小值过滤[默认: 0]')
    
    parallel_group = parser.add_argument_group('并行处理参数')
    parallel_group.add_argument('--parallel', action='store_true', help='并行处理样本')
    parallel_group.add_argument('--parallel_jobs', type=int, help='并行任务数（如果不指定，自动计算）')
    parallel_group.add_argument('--threads', type=int, default=4, help='总线程数[默认: 4]')
    
    control_group = parser.add_argument_group('兼容控制参数')
    control_group.add_argument('--skip_alignment', action='store_true', help='跳过比对步骤 (若使用 --resume 则自动包含此功能)')
    control_group.add_argument('--skip_depth', action='store_true', help='跳过深度计算步骤 (若使用 --resume 则自动包含此功能)')
    control_group.add_argument('--format', type=str, choices=['csv', 'tsv'], default='tsv', help='输出文件格式[默认: tsv]')
    
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
