#!/usr/bin/env python3
"""
统一宏转录组组装脚本
支持多种组装工具和单双端数据，支持并行处理
输出目录按样本组织
支持refineC split和merge后处理
增强版：支持各个工具/处理环节的运行时间和峰值内存（MB）统计输出
修复版：
  1. 将 refineC merge 的最终组装结果保留在上一级(样本)目录，辅助日志文件存入子目录。
  2. 简化自动检测阶段的屏幕输出，仅打印文件数、样本总数及单/双端数量。
  3. 自动解压 refineC merge 产生的最终 .gz 文件。
  4. 运行结束后自动删除产生的 *_refineC 中间子目录（除非使用 --keep-temp）。
"""

import os
import sys
import argparse
import tempfile
import shutil
import subprocess
import random
import logging
import glob
import time
import multiprocessing
import re
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import signal
import psutil

class AssemblyPipeline:
    """统一的组装流水线类"""
    
    def __init__(self):
        self.temp_dir = None
        self.output_base = None
        self.sample_output = None
        self.logger = self.setup_logging()
        self._stop_flag = False
        
    def setup_logging(self) -> logging.Logger:
        """设置日志"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.StreamHandler(sys.stdout)]
        )
        return logging.getLogger(__name__)
    
    def setup_signal_handlers(self):
        """设置信号处理器，用于优雅退出"""
        def signal_handler(signum, frame):
            self.logger.info(f"接收到信号 {signum}，正在优雅退出...")
            self._stop_flag = True
            
        signal.signal(signal.SIGINT, signal_handler)
        signal.signal(signal.SIGTERM, signal_handler)
    
    def detect_data_type_and_sample(self, input_path: str) -> Dict[str, Any]:
        """
        自动检测数据类型和样本名
        """
        input_path = Path(input_path).resolve()
        
        if self._stop_flag:
            raise KeyboardInterrupt("用户中断执行")
            
        # 支持处理的序列扩展名集合
        supported_exts = ['.fastq.gz', '.fq.gz', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fasta', '.fa', '10239.fq', '_10239.fq.gz']
            
        if input_path.is_dir():
            seq_patterns = ["*.fastq*", "*.fq*", "*.fasta*", "*.fa*"]
            seq_files = []
            for pattern in seq_patterns:
                seq_files.extend(list(input_path.glob(pattern)))
            
            if not seq_files:
                raise ValueError(f"在目录 {input_path} 中未找到fastq/fq/fasta/fa文件")
            
            # 简化输出，不再打印所有文件名
            self.logger.info(f"在目录 {input_path.name} 中找到 {len(seq_files)} 个序列文件")
            
            sample_groups = {}
            
            for f in seq_files:
                if self._stop_flag:
                    raise KeyboardInterrupt("用户中断执行")
                    
                filename = f.name
                
                # 检测配对标识并提取基础样本名
                if '_R1.' in filename or '_R1_' in filename:
                    if '_R1.' in filename:
                        sample_base = filename.split('_R1.')[0]
                    else:
                        sample_base = filename.split('_R1_')[0]
                    sample_groups.setdefault(sample_base, {})['read1'] = f.resolve()
                    
                elif '_R2.' in filename or '_R2_' in filename:
                    if '_R2.' in filename:
                        sample_base = filename.split('_R2.')[0]
                    else:
                        sample_base = filename.split('_R2_')[0]
                    sample_groups.setdefault(sample_base, {})['read2'] = f.resolve()
                    
                elif '.R1.' in filename or '.R1_' in filename:
                    if '.R1.' in filename:
                        sample_base = filename.split('.R1.')[0]
                    else:
                        sample_base = filename.split('.R1_')[0]
                    sample_groups.setdefault(sample_base, {})['read1'] = f.resolve()
                    
                elif '.R2.' in filename or '.R2_' in filename:
                    if '.R2.' in filename:
                        sample_base = filename.split('.R2.')[0]
                    else:
                        sample_base = filename.split('.R2_')[0]
                    sample_groups.setdefault(sample_base, {})['read2'] = f.resolve()
                elif '1_10239.' in filename:
                    sample_base = filename.split('_1_10239.')[0]
                    sample_groups.setdefault(sample_base, {})['read1'] = f.resolve()
                elif '2_10239.' in filename:
                    sample_base = filename.split('_2_10239.')[0]
                    sample_groups.setdefault(sample_base, {})['read2'] = f.resolve()
                elif '_1.' in filename or '_1_' in filename:
                    if '_1.' in filename:
                        sample_base = filename.split('_1.')[0]
                    else:
                        sample_base = filename.split('_1_')[0]
                    sample_groups.setdefault(sample_base, {})['read1'] = f.resolve()
                elif '_2.' in filename or '_2_' in filename:
                    if '_2.' in filename:
                        sample_base = filename.split('_2.')[0]
                    else:
                        sample_base = filename.split('_2_')[0]
                    sample_groups.setdefault(sample_base, {})['read2'] = f.resolve()
                else:
                    sample_base = filename
                    for ext in supported_exts:
                        if sample_base.endswith(ext):
                            sample_base = sample_base[:-len(ext)]
                            break
                    sample_groups.setdefault(sample_base, {})['reads'] = f.resolve()
            
            valid_samples = {}
            
            for sample_base, files in sample_groups.items():
                if self._stop_flag:
                    raise KeyboardInterrupt("用户中断执行")
                    
                if 'read1' in files and 'read2' in files:
                    valid_samples[sample_base] = {
                        'type': 'paired',
                        'read1': files['read1'],
                        'read2': files['read2']
                    }
                elif 'reads' in files:
                    valid_samples[sample_base] = {
                        'type': 'single',
                        'reads': files['reads']
                    }
                elif 'read1' in files:
                    valid_samples[sample_base] = {
                        'type': 'single', 
                        'reads': files['read1']
                    }
                elif 'read2' in files:
                    valid_samples[sample_base] = {
                        'type': 'single',
                        'reads': files['read2']
                    }
            
            # 统计单双端数量
            paired_count = sum(1 for v in valid_samples.values() if v['type'] == 'paired')
            single_count = sum(1 for v in valid_samples.values() if v['type'] == 'single')

            if len(valid_samples) > 1:
                self.logger.info(f"检测到 {len(valid_samples)} 个样本 (双端: {paired_count}, 单端: {single_count})")
                return {'type': 'batch', 'samples': valid_samples}
            elif len(valid_samples) == 1:
                self.logger.info(f"检测到 1 个样本 (双端: {paired_count}, 单端: {single_count})")
                sample_name = list(valid_samples.keys())[0]
                files_info = valid_samples[sample_name]
                if files_info['type'] == 'paired':
                    return {
                        'type': 'paired', 
                        'sample': sample_name,
                        'read1': files_info['read1'],
                        'read2': files_info['read2']
                    }
                else:
                    return {
                        'type': 'single',
                        'sample': sample_name,
                        'reads': files_info['reads']
                    }
            else:
                raise ValueError("在目录中未解析出任何有效的样本数据。")
        
        else:
            if not input_path.exists():
                raise ValueError(f"输入文件不存在: {input_path}")
            
            filename = input_path.name
            
            # 检测配对标识
            if '_R1.' in filename or '_R1_' in filename:
                if '_R1.' in filename:
                    sample_base = filename.split('_R1.')[0]
                else:
                    sample_base = filename.split('_R1_')[0]
                
                read2_patterns = [
                    input_path.parent / f"{sample_base}_R2.fa.gz",
                    input_path.parent / f"{sample_base}_R2.fasta.gz",
                    input_path.parent / f"{sample_base}_R2.fq.gz",
                    input_path.parent / f"{sample_base}_R2.fastq.gz",
                    input_path.parent / f"{sample_base}_R2.fa",
                    input_path.parent / f"{sample_base}_R2.fasta",
                    input_path.parent / f"{sample_base}_R2.fq",
                    input_path.parent / f"{sample_base}_R2.fastq",
                ]
                
                read2 = None
                for pattern in read2_patterns:
                    if pattern.exists():
                        read2 = pattern.resolve()
                        break
                
                if read2:
                    return {
                        'type': 'paired',
                        'sample': sample_base,
                        'read1': input_path.resolve(),
                        'read2': read2
                    }
                else:
                    return {
                        'type': 'single',
                        'sample': sample_base,
                        'reads': input_path.resolve()
                    }
            
            elif '_R2.' in filename or '_R2_' in filename:
                if '_R2.' in filename:
                    sample_base = filename.split('_R2.')[0]
                else:
                    sample_base = filename.split('_R2_')[0]
                
                read1_patterns = [
                    input_path.parent / f"{sample_base}_R1.fa.gz",
                    input_path.parent / f"{sample_base}_R1.fasta.gz",
                    input_path.parent / f"{sample_base}_R1.fq.gz",
                    input_path.parent / f"{sample_base}_R1.fastq.gz",
                    input_path.parent / f"{sample_base}_R1.fa",
                    input_path.parent / f"{sample_base}_R1.fasta",
                    input_path.parent / f"{sample_base}_R1.fq",
                    input_path.parent / f"{sample_base}_R1.fastq",
                ]
                
                read1 = None
                for pattern in read1_patterns:
                    if pattern.exists():
                        read1 = pattern.resolve()
                        break
                
                if read1:
                    return {
                        'type': 'paired',
                        'sample': sample_base,
                        'read1': read1,
                        'read2': input_path.resolve()
                    }
                else:
                    return {
                        'type': 'single',
                        'sample': sample_base,
                        'reads': input_path.resolve()
                    }
            
            elif '.R1.' in filename or '.R1_' in filename:
                if '.R1.' in filename:
                    sample_base = filename.split('.R1.')[0]
                else:
                    sample_base = filename.split('.R1_')[0]
                
                read2_patterns = [
                    input_path.parent / f"{sample_base}.R2.fa.gz",
                    input_path.parent / f"{sample_base}.R2.fasta.gz",
                    input_path.parent / f"{sample_base}.R2.fq.gz",
                    input_path.parent / f"{sample_base}.R2.fastq.gz",
                    input_path.parent / f"{sample_base}_R2.fa.gz",
                    input_path.parent / f"{sample_base}_R2.fasta.gz",
                    input_path.parent / f"{sample_base}_R2.fq.gz",
                    input_path.parent / f"{sample_base}_R2.fastq.gz",
                ]
                
                read2 = None
                for pattern in read2_patterns:
                    if pattern.exists():
                        read2 = pattern.resolve()
                        break
                
                if read2:
                    return {
                        'type': 'paired',
                        'sample': sample_base,
                        'read1': input_path.resolve(),
                        'read2': read2
                    }
                else:
                    return {
                        'type': 'single',
                        'sample': sample_base,
                        'reads': input_path.resolve()
                    }
            
            elif '.R2.' in filename or '.R2_' in filename:
                if '.R2.' in filename:
                    sample_base = filename.split('.R2.')[0]
                else:
                    sample_base = filename.split('.R2_')[0]
                
                read1_patterns = [
                    input_path.parent / f"{sample_base}.R1.fa.gz",
                    input_path.parent / f"{sample_base}.R1.fasta.gz",
                    input_path.parent / f"{sample_base}.R1.fq.gz",
                    input_path.parent / f"{sample_base}.R1.fastq.gz",
                    input_path.parent / f"{sample_base}_R1.fa.gz",
                    input_path.parent / f"{sample_base}_R1.fasta.gz",
                    input_path.parent / f"{sample_base}_R1.fq.gz",
                    input_path.parent / f"{sample_base}_R1.fastq.gz",
                ]
                
                read1 = None
                for pattern in read1_patterns:
                    if pattern.exists():
                        read1 = pattern.resolve()
                        break
                
                if read1:
                    return {
                        'type': 'paired',
                        'sample': sample_base,
                        'read1': read1,
                        'read2': input_path.resolve()
                    }
                else:
                    return {
                        'type': 'single',
                        'sample': sample_base,
                        'reads': input_path.resolve()
                    }
                    
            elif '_1.' in filename or '_1_' in filename:
                if '_1.' in filename:
                    sample_base = filename.split('_1.')[0]
                else:
                    sample_base = filename.split('_1_')[0]
                
                read2_patterns = [
                    input_path.parent / f"{sample_base}_2.fa.gz",
                    input_path.parent / f"{sample_base}_2.fasta.gz",
                    input_path.parent / f"{sample_base}_2.fq.gz",
                    input_path.parent / f"{sample_base}_2.fastq.gz",
                    input_path.parent / f"{sample_base}_2.fa",
                    input_path.parent / f"{sample_base}_2.fasta",
                    input_path.parent / f"{sample_base}_2.fq",
                    input_path.parent / f"{sample_base}_2.fastq",
                ]
                
                read2 = None
                for pattern in read2_patterns:
                    if pattern.exists():
                        read2 = pattern.resolve()
                        break
                
                if read2:
                    return {
                        'type': 'paired',
                        'sample': sample_base,
                        'read1': input_path.resolve(),
                        'read2': read2
                    }
                else:
                    return {
                        'type': 'single',
                        'sample': sample_base,
                        'reads': input_path.resolve()
                    }
                    
            elif '_2.' in filename or '_2_' in filename:
                if '_2.' in filename:
                    sample_base = filename.split('_2.')[0]
                else:
                    sample_base = filename.split('_2_')[0]
                
                read1_patterns = [
                    input_path.parent / f"{sample_base}_1.fa.gz",
                    input_path.parent / f"{sample_base}_1.fasta.gz",
                    input_path.parent / f"{sample_base}_1.fq.gz",
                    input_path.parent / f"{sample_base}_1.fastq.gz",
                    input_path.parent / f"{sample_base}_1.fa",
                    input_path.parent / f"{sample_base}_1.fasta",
                    input_path.parent / f"{sample_base}_1.fq",
                    input_path.parent / f"{sample_base}_1.fastq",
                ]
                
                read1 = None
                for pattern in read1_patterns:
                    if pattern.exists():
                        read1 = pattern.resolve()
                        break
                
                if read1:
                    return {
                        'type': 'paired',
                        'sample': sample_base,
                        'read1': read1,
                        'read2': input_path.resolve()
                    }
                else:
                    return {
                        'type': 'single',
                        'sample': sample_base,
                        'reads': input_path.resolve()
                    }
            
            else:
                sample_base = filename
                for ext in supported_exts:
                    if sample_base.endswith(ext):
                        sample_base = sample_base[:-len(ext)]
                        break
                return {
                    'type': 'single',
                    'sample': sample_base,
                    'reads': input_path.resolve()
                }
    
    def parse_arguments(self) -> argparse.Namespace:
        """解析命令行参数"""
        parser = argparse.ArgumentParser(
            description='统一宏转录组组装流水线',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
使用示例:
  python assembly_pipeline.py -t megahit -i sample1_1.fa.gz
  python assembly_pipeline.py -t rnaviralspades -i /path/to/fasta_files/ -j 4
  python assembly_pipeline.py -t penguin -i sample.fq.gz -l 0 -n 8
  python assembly_pipeline.py -t all -i sample_R1.fq.gz -j 2
  python assembly_pipeline.py -t all -i sample_1.fa.gz --refineC_split --refineC_merge --log_dirs ./logs
            """
        )
        
        # 必需参数
        parser.add_argument('-t', '--tool', 
                          choices=['megahit', 'rnaviralspades', 'penguin', 'all'],
                          required=True,
                          help='组装工具选择')
        parser.add_argument('-i', '--input',
                          required=True,
                          help='输入文件或目录。自动检测单双端数据和样本名(支持 fastq, fq, fasta, fa)')
        
        # 组装参数
        parser.add_argument('-l', '--length', type=int, default=0,
                          help='contig最小长度阈值 (默认: 0)')
        parser.add_argument('-n', '--threads', type=int, default=8,
                          help='每个任务的线程数 (默认: 8)')
        parser.add_argument('-m', '--memory', type=int, default=64,
                          help='每个任务的内存大小(GB) (默认: 64)')
        parser.add_argument('-j', '--jobs', type=int, default=1,
                          help='并行任务数 (默认: 1，串行运行)')
        
        # refineC后处理参数
        parser.add_argument('--refineC_split', action='store_true',
                          help='运行refineC split后处理，分割重叠的contigs')
        parser.add_argument('--refineC_merge', action='store_true',
                          help='运行refineC merge后处理，合并来自不同工具的contigs（仅当使用多个工具时有效）')
        parser.add_argument('--refineC_threads', type=int, default=None,
                          help='refineC使用的线程数 (默认: 使用--threads参数的值)')
        parser.add_argument('--refineC_frag_min_len', type=int, default=1000,
                          help='refineC split的最小片段长度 (默认: 1000)')
        parser.add_argument('--refineC_min_id', type=float, default=0.95,
                          help='refineC merge的最小序列一致性 (默认: 0.95)')
        parser.add_argument('--refineC_min_cov', type=float, default=0.50,
                          help='refineC merge的最小覆盖度 (默认: 0.50)')
        
        # 输出参数
        parser.add_argument('-o', '--output-dir', default='./results',
                          help='输出目录 (默认: ./results)')
        parser.add_argument('--tmp-dir',
                          help='临时目录 (默认: 系统临时目录)')
        parser.add_argument('--keep-temp', action='store_true',
                          help='保留临时文件和refineC中间目录')
        parser.add_argument('--force', action='store_true',
                          help='强制重新运行，覆盖已有结果')
                          
        # 资源统计输出
        parser.add_argument('--log_dirs', default='./logs',
                          help='运行时间和内存统计日志输出目录 (默认: ./logs)')
        
        return parser.parse_args()

    # --------------- 资源与日志统计相关辅助方法 ---------------
    def _init_stats(self) -> Dict[str, Dict[str, float]]:
        """初始化一个包含各个工具资源的统计字典"""
        return {
            'megahit': {'time': 0.0, 'mem': 0.0},
            'rnaviralspades': {'time': 0.0, 'mem': 0.0},
            'penguin': {'time': 0.0, 'mem': 0.0},
            'refineC_split': {'time': 0.0, 'mem': 0.0},
            'refineC_merge': {'time': 0.0, 'mem': 0.0},
            'total': {'time': 0.0, 'mem': 0.0}
        }

    def _get_max_mem(self, stats: Dict[str, Dict[str, float]]) -> float:
        """从统计结果中得出最高内存峰值"""
        return max([s['mem'] for k, s in stats.items() if k != 'total'], default=0.0)

    def parse_time_mem(self, time_log: Path, py_time: float) -> Tuple[float, float]:
        """解析 .time.mem.log 提取时间和内存"""
        mem_mb = 0.0
        time_s = py_time
        if time_log and time_log.exists():
            try:
                with open(time_log, 'r') as f:
                    content = f.read()
                    time_match = re.search(r'Time:\s*([\d.]+)\s*seconds', content, re.IGNORECASE)
                    if time_match:
                        time_s = float(time_match.group(1))
                    mem_match = re.search(r'Memory:\s*(\d+)\s*KB', content, re.IGNORECASE)
                    if mem_match:
                        mem_mb = float(mem_match.group(1)) / 1024.0
                    cpu_match = re.search(r'CPU:\s*([\d]+)%', content)
                    if cpu_match:
                        cpu_pct = float(cpu_match.group(1))
            except Exception as e:
                self.logger.warning(f"解析日志文件失败 {time_log}: {e}")
        return time_s, mem_mb

    def write_stats_log(self, args: argparse.Namespace, sample: str, stats: Dict[str, Dict[str, float]]) -> None:
        """将单个样本的汇总统计追加到日志文件"""
        if not args.log_dirs:
            return
        log_dir = Path(args.log_dirs)
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / "resource_usage_summary.tsv"
        
        write_header = not log_file.exists()
        
        try:
            with open(log_file, 'a') as f:
                if write_header:
                    headers = [
                        "Sample",
                        "megahit_time(s)", "megahit_mem(MB)",
                        "rnaviralspades_time(s)", "rnaviralspades_mem(MB)",
                        "penguin_time(s)", "penguin_mem(MB)",
                        "M_split_time(s)", "M_split_mem(MB)", "H_split_time(s)", "H_split_mem(MB)", "P_split_time(s)", "P_split_mem(MB)",
                        "refineC_merge_time(s)", "refineC_merge_mem(MB)",
                        "total_time(s)", "total_mem(MB)"
                    ]
                    f.write("\t".join(headers) + "\n")
                
                row = [
                    sample,
                    f"{stats['megahit']['time']:.2f}", f"{stats['megahit']['mem']:.2f}",
                    f"{stats['rnaviralspades']['time']:.2f}", f"{stats['rnaviralspades']['mem']:.2f}",
                    f"{stats['penguin']['time']:.2f}", f"{stats['penguin']['mem']:.2f}",
                    f"{stats.get('refineC_split_megahit', {}).get('time', 0):.2f}", f"{stats.get('refineC_split_megahit', {}).get('mem', 0):.2f}", f"{stats.get('refineC_split_rnaviralspades', {}).get('time', 0):.2f}", f"{stats.get('refineC_split_rnaviralspades', {}).get('mem', 0):.2f}", f"{stats.get('refineC_split_penguin', {}).get('time', 0):.2f}", f"{stats.get('refineC_split_penguin', {}).get('mem', 0):.2f}",
                    f"{stats['refineC_merge']['time']:.2f}", f"{stats['refineC_merge']['mem']:.2f}",
                    f"{stats['total']['time']:.2f}", f"{stats['total']['mem']:.2f}"
                ]
                f.write("\t".join(row) + "\n")
        except Exception as e:
            self.logger.warning(f"写入资源统计日志失败: {e}")
    # -----------------------------------------------------------

    def check_system_resources(self) -> Tuple[bool, str]:
        """检查系统资源"""
        try:
            # 检查内存
            memory_gb = psutil.virtual_memory().total / (1024**3)
            if memory_gb < 4:
                return False, f"系统内存不足: {memory_gb:.1f}GB < 4GB"
            
            # 检查磁盘空间
            disk_usage = psutil.disk_usage('/')
            free_disk_gb = disk_usage.free / (1024**3)
            if free_disk_gb < 10:
                return False, f"磁盘空间不足: {free_disk_gb:.1f}GB < 10GB"
                
            return True, "系统资源充足"
            
        except Exception as e:
            self.logger.warning(f"系统资源检查失败: {e}")
            return True, "资源检查跳过"
    
    def check_tool_dependencies(self, tools: List[str], args: argparse.Namespace) -> Tuple[bool, List[str]]:
        """检查工具依赖"""
        tools_executables = {
            'megahit': 'megahit',
            'rnaviralspades': 'rnaviralspades.py',
            'penguin': 'penguin',
            'refineC': 'refineC'
        }
        
        missing_tools = []
        for tool in tools:
            if tool == 'refineC' and (args.refineC_split or args.refineC_merge):
                tool_executable = tools_executables[tool]
                if not shutil.which(tool_executable):
                    missing_tools.append(tool_executable)
            elif tool != 'refineC':
                tool_executable = tools_executables[tool]
                if not shutil.which(tool_executable):
                    missing_tools.append(tool_executable)
        
        # 检查 bioawk
        if 'rnaviralspades' in tools and not shutil.which('bioawk'):
            missing_tools.append('bioawk')
        
        return len(missing_tools) == 0, missing_tools
    
    def validate_arguments(self, args: argparse.Namespace) -> bool:
        """验证参数有效性"""
        errors = []
        
        # 验证输入路径
        if not os.path.exists(args.input):
            errors.append(f"输入路径不存在: {args.input}")
        
        # 验证数值参数
        if args.threads < 1:
            errors.append("线程数必须大于0")
        
        if args.memory < 1:
            errors.append("内存大小必须大于0")
        
        if args.jobs < 1:
            errors.append("并行任务数必须大于0")
        
        # all 模式自动启用 refineC split + merge
        if args.tool == 'all':
            args.refineC_split = True
            args.refineC_merge = True

        # 检查工具依赖
        tools_to_check = []
        if args.tool == 'all':
            tools_to_check = ['megahit', 'rnaviralspades', 'penguin']
        else:
            tools_to_check = [args.tool]

        # refineC 后处理添加到依赖检查
        if args.refineC_split or args.refineC_merge:
            tools_to_check.append('refineC')

        deps_ok, missing_tools = self.check_tool_dependencies(tools_to_check, args)
        if not deps_ok:
            errors.append(f"缺少必要工具: {', '.join(missing_tools)}")
        
        # 检查系统资源
        resource_ok, resource_msg = self.check_system_resources()
        if not resource_ok:
            errors.append(resource_msg)
        else:
            self.logger.info(resource_msg)
        
        if errors:
            for error in errors:
                self.logger.error(error)
            return False
        
        return True
    
    def setup_directories(self, args: argparse.Namespace, sample: str, data_type: str) -> None:
        """设置工作目录 - 修改为按样本组织"""
        # 输出目录直接按样本组织，不按工具分类
        self.output_base = Path(args.output_dir)
        self.sample_output = self.output_base / sample
        self.sample_output.mkdir(parents=True, exist_ok=True)
        
        # 创建临时目录
        if args.tmp_dir:
            temp_dir_name = f"tmp_{sample}_{random.randint(1000, 9999)}"
            self.temp_dir = Path(args.tmp_dir) / temp_dir_name
        else:
            temp_dir_name = f"assembly_{sample}_"
            self.temp_dir = Path(tempfile.mkdtemp(prefix=temp_dir_name))
        
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"临时目录: {self.temp_dir}")
        self.logger.info(f"样本输出目录: {self.sample_output}")
    
    def prepare_input_files(self, sample_info: Dict[str, Any], sample: str) -> Dict[str, Path]:
        """准备输入文件"""
        input_files = {}
        
        if sample_info['type'] == 'paired':
            input_files['read1'] = sample_info['read1']
            input_files['read2'] = sample_info['read2']
            self.logger.info(f"使用双端数据: {sample_info['read1'].name}, {sample_info['read2'].name}")
        else:
            input_files['reads'] = sample_info['reads']
            self.logger.info(f"使用单端数据: {sample_info['reads'].name}")
        
        return input_files
    
    def build_megahit_command(self, args: argparse.Namespace, input_files: Dict[str, Path]) -> List[str]:
        """构建MEGAHIT命令"""
        cmd = [
            'megahit',
            '--k-list',  '21,31,41,51,61,71,81,91,99',
            '--out-dir', str(self.temp_dir / "assembly"),
            '--out-prefix', 'megahit',
            '--min-contig-len', str(args.length),
            '--tmp-dir', str(self.temp_dir),
            '-t', str(args.threads)
        ]
        
        if 'read1' in input_files and 'read2' in input_files:
            cmd.extend(['-1', str(input_files['read1']), '-2', str(input_files['read2'])])
        else:
            cmd.extend(['-r', str(input_files['reads'])])
        
        return cmd
    
    def build_rnaviralspades_command(self, args: argparse.Namespace, input_files: Dict[str, Path]) -> List[str]:
        """构建rnaviralspades命令"""
        cmd = [
            'rnaviralspades.py',
            '-o', str(self.temp_dir / "assembly"),
            '--tmp-dir', str(self.temp_dir),
            '-t', str(args.threads),
            '-m', str(args.memory)
        ]
        
        if 'read1' in input_files and 'read2' in input_files:
            cmd.extend(['-1', str(input_files['read1']), '-2', str(input_files['read2'])])
        else:
            cmd.extend(['-s', str(input_files['reads'])])
        
        return cmd
    
    def build_penguin_command(self, args: argparse.Namespace, input_files: Dict[str, Path], sample: str) -> List[str]:
        """构建penguin命令"""
        output_file = self.temp_dir / "penguin_out" / f"{sample}.penguin.assembly.fa"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            'penguin',
            'guided_nuclassemble'
        ]
        
        if 'read1' in input_files and 'read2' in input_files:
            cmd.extend([str(input_files['read1']), str(input_files['read2'])])
        else:
            cmd.append(str(input_files['reads']))
        
        cmd.extend([
            '--threads', str(args.threads),
            str(output_file),
            str(self.temp_dir / "penguin_tmp"),
            '--min-contig-len',str(args.length)
        ])
        
        return cmd
    
    def run_assembly(self, args: argparse.Namespace, input_files: Dict[str, Path], sample: str, tool: str, stats: Dict) -> bool:
        """运行组装工具"""
        if self._stop_flag:
            self.logger.info(f"用户中断，停止 {tool} 组装")
            return False
            
        self.logger.info(f"开始 {tool} 组装")
        start_time = time.time()
        
        # 构建命令
        command_builders = {
            'megahit': self.build_megahit_command,
            'rnaviralspades': self.build_rnaviralspades_command,
            'penguin': lambda a, f: self.build_penguin_command(a, f, sample)
        }
        
        cmd = command_builders[tool](args, input_files)
        
        # 准备日志文件
        time_log = self.temp_dir / f"{sample}.{tool}.time.mem.log"
        stdout_log = self.temp_dir / f"{sample}.{tool}.log"
        
        self.logger.info(f"执行命令: {' '.join(cmd)}")
        
        try:
            # 使用time命令记录资源使用
            time_exec = shutil.which('gtime') or shutil.which('/usr/bin/time') or shutil.which('time')
            if not time_exec:
                self.logger.warning("未找到 time 命令，将不记录资源使用情况")
                cmd_to_run = cmd
            else:
                cmd_to_run = [time_exec, '-f', 'Time:%e seconds\nMemory:%M KB\nCPU:%P', '--output', str(time_log)] + cmd
            
            with open(stdout_log, 'w') as log_file:
                process = subprocess.run(
                    cmd_to_run,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=True,
                    timeout=86400  # 24小时超时
                )
            
            elapsed_time = time.time() - start_time
            time_s, mem_mb = self.parse_time_mem(time_log, elapsed_time)
            stats[tool]['time'] = time_s
            stats[tool]['mem'] = mem_mb
            
            self.logger.info(f"{tool} 组装成功完成，耗时: {elapsed_time:.2f} 秒, 内存: {mem_mb:.2f} MB")
            return True
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"组装过程超时 (24小时)")
            return False
        except subprocess.CalledProcessError as e:
            self.logger.error(f"组装过程失败，退出码: {e.returncode}")
            # 记录详细的错误信息
            stdout_log_path = self.temp_dir / f"{sample}.{tool}.log"
            if stdout_log_path.exists():
                with open(stdout_log_path, 'r') as f:
                    error_details = f.read()[-2000:]  # 只显示最后2000字符
                self.logger.error(f"错误日志片段:\n{error_details}")
            return False
        except FileNotFoundError as e:
            self.logger.error(f"找不到命令: {e}")
            return False
        except Exception as e:
            self.logger.error(f"组装过程中发生未知错误: {e}")
            return False
    
    def is_fasta_file_valid(self, fasta_file: Path) -> bool:
        """检查FASTA文件是否有效且非空"""
        try:
            if not fasta_file.exists():
                return False
                
            # 检查文件大小
            if fasta_file.stat().st_size == 0:
                return False
            
            # 如果是压缩文件，检查是否能够读取
            if fasta_file.suffix == '.gz':
                import gzip
                with gzip.open(fasta_file, 'rt') as f:
                    first_line = f.readline()
                    return first_line.startswith('>')
            else:
                with open(fasta_file, 'r') as f:
                    first_line = f.readline()
                    return first_line.startswith('>')
        except Exception as e:
            self.logger.warning(f"检查FASTA文件时出错: {e}")
            return False
    
    def run_refineC_split(self, args: argparse.Namespace, input_fasta: Path, sample: str, tool: str, stats: Dict, prefix_suffix: str = "") -> Optional[Path]:
        """运行refineC split后处理"""
        if self._stop_flag:
            return None
            
        self.logger.info(f"开始 refineC split 处理: {tool}")
        start_time = time.time()
        
        # 设置输出文件名和目录
        split_suffix = prefix_suffix if prefix_suffix else ""
        if prefix_suffix:
            output_prefix = f"{sample}_{tool}.refineC{split_suffix}"
            split_output_dir = self.sample_output / f"{sample}_{tool}_refineC{split_suffix}"
        else:
            output_prefix = f"{sample}_{tool}.refineC"
            split_output_dir = self.sample_output / f"{sample}_{tool}_refineC"
        
        split_output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置refineC线程数
        refineC_threads = args.refineC_threads if args.refineC_threads else args.threads
        
        # 构建refineC split命令
        cmd = [
            'refineC', 'split',
            '--threads', str(refineC_threads),
            '--contigs', str(input_fasta),
            '--prefix', output_prefix,
            '--output', str(split_output_dir),
            '--frag-min-len', str(args.refineC_frag_min_len),
            '--min-id',str(args.refineC_min_id),
            '--min-cov', str(args.refineC_min_cov)

        ]
        
        # 准备日志文件
        time_log = split_output_dir / f"{output_prefix}.time.mem.log"
        stdout_log = split_output_dir / f"{output_prefix}.log"
        
        self.logger.info(f"执行命令: {' '.join(cmd)}")
        
        try:
            # 使用time命令记录资源使用
            time_exec = shutil.which('gtime') or shutil.which('/usr/bin/time') or shutil.which('time')
            if not time_exec:
                self.logger.warning("未找到 time 命令，将不记录资源使用情况")
                cmd_to_run = cmd
            else:
                cmd_to_run = [time_exec, '-f', 'Time:%e seconds\nMemory:%M KB\nCPU:%P', '--output', str(time_log)] + cmd
            
            with open(stdout_log, 'w') as log_file:
                process = subprocess.run(
                    cmd_to_run,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=True,
                    timeout=86400  # 24小时超时
                )
            
            elapsed_time = time.time() - start_time
            time_s, mem_mb = self.parse_time_mem(time_log, elapsed_time)
            
            # 因为 split 可能会对多个 tool 执行，我们取时间累加，内存取最大值
            stats[f'refineC_split_{tool}'] = {'time': time_s, 'mem': mem_mb}
            
            self.logger.info(f"refineC split 成功完成，耗时: {elapsed_time:.2f} 秒, 内存: {mem_mb:.2f} MB")
            
            # 查找输出文件
            base_prefix = f"{sample}_{tool}_refineC"
            possible_locations = [
                self.sample_output / f"{base_prefix}.split.fasta.gz",
                self.sample_output / f"{base_prefix}.split.fasta",
                self.sample_output / f"{sample}_{tool}.refineC.split.fasta.gz",
                self.sample_output / f"{sample}_{tool}.refineC.split.fasta",
                split_output_dir / f"{base_prefix}.split.fasta.gz",
                split_output_dir / f"{base_prefix}.split.fasta",
                split_output_dir / f"{output_prefix}.split.fasta.gz",
                split_output_dir / f"{output_prefix}.split.fasta",
            ]
            
            for location in possible_locations:
                if location.exists():
                    self.logger.info(f"找到split输出文件: {location}")
                    
                    if location.parent == self.sample_output and location.parent != split_output_dir:
                        target_location = split_output_dir / location.name
                        self.logger.info(f"将文件从 {location} 移动到 {target_location}")
                        shutil.move(location, target_location)
                        
                        related_pattern = location.name.replace('.fasta.gz', '.*').replace('.fasta', '.*')
                        for related_file in self.sample_output.glob(related_pattern):
                            if related_file != location:
                                target_related = split_output_dir / related_file.name
                                shutil.move(related_file, target_related)
                                self.logger.info(f"移动相关文件: {related_file.name}")
                        
                        return target_location
                    return location
            
            split_patterns = ["*.split.fasta.gz", "*.split.fasta"]
            for pattern in split_patterns:
                found_files = list(self.sample_output.glob(pattern))
                if found_files:
                    for found_file in found_files:
                        if tool in str(found_file):
                            self.logger.info(f"通过模式匹配找到split文件: {found_file}")
                            return found_file
            
            self.logger.warning(f"未找到任何split输出文件")
            self.logger.info(f"refineC split 没有生成输出fasta文件，将使用原始文件作为split结果")
            
            if input_fasta.suffix == '.gz':
                split_file = split_output_dir / f"{base_prefix}.split.fasta.gz"
                shutil.copy2(input_fasta, split_file)
            else:
                split_file = split_output_dir / f"{base_prefix}.split.fasta"
                if input_fasta.suffix in ['.gz', '.bz2', '.xz']:
                    import gzip
                    with gzip.open(input_fasta, 'rt') as f_in:
                        with open(split_file, 'w') as f_out:
                            f_out.write(f_in.read())
                else:
                    shutil.copy2(input_fasta, split_file)
            
            self.logger.info(f"使用原始文件作为split结果: {split_file}")
            return split_file
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"refineC split 过程超时 (24小时)")
            return None
        except subprocess.CalledProcessError as e:
            self.logger.error(f"refineC split 过程失败，退出码: {e.returncode}")
            if stdout_log.exists():
                with open(stdout_log, 'r') as f:
                    error_details = f.read()[-2000:]
                self.logger.error(f"错误日志片段:\n{error_details}")
            return None
        except FileNotFoundError as e:
            self.logger.error(f"找不到refineC命令: {e}")
            return None
        except Exception as e:
            self.logger.error(f"refineC split 过程中发生未知错误: {e}")
            return None
    
    def run_refineC_merge(self, args: argparse.Namespace, input_fasta: Path, sample: str, stats: Dict, prefix: str = None) -> Optional[Path]:
        """运行refineC merge后处理"""
        if self._stop_flag:
            return None
            
        self.logger.info(f"开始 refineC merge 处理")
        start_time = time.time()
        
        # 设置输出文件名和目录
        if not prefix:
            prefix = f"{sample}.all_tools.refineC"
        
        merge_output_dir = self.sample_output / f"{sample}_all_tools_refineC_merge"
        merge_output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置refineC线程数
        refineC_threads = args.refineC_threads if args.refineC_threads else args.threads
        
        # 构建refineC merge命令
        cmd = [
            'refineC', 'merge',
            '--threads', str(refineC_threads),
            '--contigs', str(input_fasta),
            '--prefix', prefix,
            '--output', str(merge_output_dir),
            '--min-id', str(args.refineC_min_id),
            '--min-cov', str(args.refineC_min_cov),
            '--mnm2-maxtrim','100',
            '--mnm2-overlap','150',
            '--glob-cls-id','0.95'
        ]
        
        # 准备日志文件
        time_log = merge_output_dir / f"{prefix}.time.mem.log"
        stdout_log = merge_output_dir / f"{prefix}.log"
        
        self.logger.info(f"执行命令: {' '.join(cmd)}")
        
        try:
            # 使用time命令记录资源使用
            time_exec = shutil.which('gtime') or shutil.which('/usr/bin/time') or shutil.which('time')
            if not time_exec:
                self.logger.warning("未找到 time 命令，将不记录资源使用情况")
                cmd_to_run = cmd
            else:
                cmd_to_run = [time_exec, '-f', 'Time:%e seconds\nMemory:%M KB\nCPU:%P', '--output', str(time_log)] + cmd
            
            with open(stdout_log, 'w') as log_file:
                process = subprocess.run(
                    cmd_to_run,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=True,
                    timeout=86400  # 24小时超时
                )
            
            elapsed_time = time.time() - start_time
            time_s, mem_mb = self.parse_time_mem(time_log, elapsed_time)
            stats['refineC_merge']['time'] = time_s
            stats['refineC_merge']['mem'] = mem_mb
            
            self.logger.info(f"refineC merge 成功完成，耗时: {elapsed_time:.2f} 秒, 内存: {mem_mb:.2f} MB")
            
            dir_name = merge_output_dir.name
            
            # 处理gz解压的辅助函数
            def finalize_merge_file(file_path: Path) -> Path:
                if file_path.suffix == '.gz':
                    unzipped = file_path.with_suffix('')
                    import gzip
                    self.logger.info(f"正在解压合并结果: {file_path.name} -> {unzipped.name}")
                    with gzip.open(file_path, 'rt') as f_in:
                        with open(unzipped, 'w') as f_out:
                            f_out.write(f_in.read())
                    file_path.unlink()  # 删除原压缩包
                    return unzipped
                return file_path
            
            # 搜索产生的文件位置
            possible_files = [
                self.sample_output / f"{dir_name}.merged.fasta.gz",
                self.sample_output / f"{dir_name}.merged.fasta",
                merge_output_dir / f"{dir_name}.merged.fasta.gz",
                merge_output_dir / f"{dir_name}.merged.fasta",
                merge_output_dir / f"{prefix}.merged.fasta.gz",
                merge_output_dir / f"{prefix}.merged.fasta",
                merge_output_dir / "merged.fasta.gz",
                merge_output_dir / "merged.fasta",
            ]
            
            for possible_file in possible_files:
                if possible_file.exists():
                    self.logger.info(f"找到merge输出文件: {possible_file}")
                    
                    # 1. 如果文件原本就在样本根目录
                    if possible_file.parent == self.sample_output:
                        self.logger.info(f"文件保留在样本根目录作为最终结果: {possible_file}")
                        
                        related_pattern = possible_file.name.replace('.fasta.gz', '.*').replace('.fasta', '.*')
                        for related_file in self.sample_output.glob(related_pattern):
                            if related_file != possible_file and related_file.name != possible_file.name:
                                target_related = merge_output_dir / related_file.name
                                shutil.move(related_file, target_related)
                                self.logger.info(f"移动辅助文件到子目录: {related_file.name}")
                        return finalize_merge_file(possible_file)
                    
                    # 2. 如果文件不幸被写到了子目录中，把它往上提一级
                    elif possible_file.parent == merge_output_dir:
                        target_file = self.sample_output / possible_file.name
                        self.logger.info(f"将最终结果文件从子目录移动到上一级目录: {target_file}")
                        shutil.move(possible_file, target_file)
                        return finalize_merge_file(target_file)
                        
                    return finalize_merge_file(possible_file)
            
            # 如果上面没找到，再执行模糊搜索
            fasta_files_in_sample = list(self.sample_output.glob("*.merged.fasta*"))
            for fasta_file in fasta_files_in_sample:
                if 'merge' in fasta_file.name.lower():
                    self.logger.info(f"在样本目录中找到merge文件，作为最终结果保留: {fasta_file}")
                    return finalize_merge_file(fasta_file)
            
            fasta_files_in_output = list(merge_output_dir.glob("*.fasta*"))
            for fasta_file in fasta_files_in_output:
                if 'merge' in fasta_file.name.lower():
                    self.logger.info(f"在输出目录中找到merge文件: {fasta_file}")
                    target_file = self.sample_output / fasta_file.name
                    shutil.move(fasta_file, target_file)
                    self.logger.info(f"移动到样本根目录作为最终结果: {target_file}")
                    return finalize_merge_file(target_file)
            
            self.logger.warning(f"未找到任何merge输出文件")
            return None
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"refineC merge 过程超时 (24小时)")
            return None
        except subprocess.CalledProcessError as e:
            self.logger.error(f"refineC merge 过程失败，退出码: {e.returncode}")
            if stdout_log.exists():
                with open(stdout_log, 'r') as f:
                    error_details = f.read()[-2000:]
                self.logger.error(f"错误日志片段:\n{error_details}")
            return None
        except FileNotFoundError as e:
            self.logger.error(f"找不到refineC命令: {e}")
            return None
        except Exception as e:
            self.logger.error(f"refineC merge 过程中发生未知错误: {e}")
            return None
    
    def get_output_filename(self, tool: str, sample: str) -> str:
        """获取输出文件名"""
        filenames = {
            'megahit': f'{sample}_{tool}.contig.fasta',
            'rnaviralspades': f'{sample}_{tool}.contig.fasta',
            'penguin': f'{sample}_{tool}.contig.fasta'
        }
        return filenames[tool]
    
    def get_log_filename(self, tool: str, sample: str) -> str:
        """获取日志文件名"""
        return f"{sample}.{tool}.log"
    
    def get_time_log_filename(self, tool: str, sample: str) -> str:
        """获取时间内存日志文件名"""
        return f"{sample}.{tool}.time.mem.log"
    
    def copy_results(self, args: argparse.Namespace, sample: str, tool: str, stats: Dict) -> bool:
        """复制结果文件到样本目录"""
        try:
            # 复制日志文件
            log_files = [
                (self.temp_dir / f"{sample}.{tool}.log", self.get_log_filename(tool, sample)),
                (self.temp_dir / f"{sample}.{tool}.time.mem.log", self.get_time_log_filename(tool, sample))
            ]
        
            for src_file, dst_name in log_files:
                if src_file.exists():
                    dst_file = self.sample_output / dst_name
                    shutil.copy2(src_file, dst_file)
                    self.logger.info(f"日志文件已复制: {dst_file}")
                else:
                    self.logger.warning(f"日志文件不存在: {src_file}")
        
            # 复制组装结果
            output_filename = self.get_output_filename(tool, sample)
        
            if tool == 'penguin':
                source_file = self.temp_dir / "penguin_out" / f"{sample}.penguin.assembly.fa"
                if source_file.exists():
                    dest_file = self.sample_output / output_filename
                    shutil.copy2(source_file, dest_file)
                    self.logger.info(f"组装结果已复制: {dest_file}")
                else:
                    self.logger.error(f"找不到penguin输出文件: {source_file}")
                    return False
            else:
                assembly_dir = self.temp_dir / "assembly"
                if assembly_dir.exists():
                    if tool == 'megahit':
                        original_file = assembly_dir / "megahit.contigs.fa"
                        if original_file.exists():
                            dest_file = self.sample_output / output_filename
                            shutil.copy2(original_file, dest_file)
                            self.logger.info(f"组装结果已复制: {dest_file}")
                        else:
                            self.logger.error(f"找不到组装输出文件: {original_file}")
                            return False
                    else:  # rnaviralspades
                        # 复制 contigs.fasta
                        original_contigs = assembly_dir / "contigs.fasta"
                        if original_contigs.exists():
                            dest_contigs = self.sample_output / output_filename
                            shutil.copy2(original_contigs, dest_contigs)
                            self.logger.info(f"contigs结果已复制: {dest_contigs}")
                        else:
                            self.logger.error(f"找不到contigs输出文件: {original_contigs}")
                            return False
                    
                        # 复制 scaffolds.fasta
                        original_scaffolds = assembly_dir / "scaffolds.fasta"
                        scaffolds_filename = f'{sample}_{tool}.scaffolds.fasta'
                        dest_scaffolds = self.sample_output / scaffolds_filename
                    
                        if original_scaffolds.exists():
                            shutil.copy2(original_scaffolds, dest_scaffolds)
                            self.logger.info(f"scaffolds结果已复制: {dest_scaffolds}")
                        
                            # 在scaffolds头部添加样本名前缀
                            self.add_sample_prefix(dest_scaffolds, sample)
                        
                            # 验证scaffolds文件
                            if not self.validate_output_file(dest_scaffolds):
                                self.logger.warning(f"scaffolds文件验证失败: {dest_scaffolds}")
                        else:
                            self.logger.warning(f"找不到scaffolds输出文件: {original_scaffolds}")
                else:
                    self.logger.error(f"找不到组装目录: {assembly_dir}")
                    return False
        
            # 在contig头部添加样本名前缀
            final_output = self.sample_output / output_filename
            if final_output.exists():
                self.add_sample_prefix(final_output, sample)
                # 验证输出文件
                if self.validate_output_file(final_output):
                    self.logger.info(f"结果文件已保存并验证: {final_output}")
                    
                    # 如果启用了refineC_split，运行refineC split
                    if args.refineC_split:
                        self.logger.info(f"对 {tool} 结果运行 refineC split")
                        split_result = self.run_refineC_split(args, final_output, sample, tool, stats)
                        if split_result:
                            self.logger.info(f"refineC split 完成: {split_result}")
                        else:
                            # 即使refineC split失败，我们仍然认为这个工具是成功的
                            self.logger.warning(f"refineC split 失败或无输出，但组装成功: {tool}")
                else:
                    self.logger.warning(f"输出文件验证失败: {final_output}")
            else:
                self.logger.warning(f"未找到预期的输出文件: {final_output}")
                return False
        
            return True
        
        except Exception as e:
            self.logger.error(f"复制结果时出错: {e}")
            return False

    def validate_output_file(self, fasta_file: Path) -> bool:
        """验证输出FASTA文件"""
        try:
            with open(fasta_file, 'r') as f:
                content = f.read().strip()
            
            if not content:
                self.logger.warning(f"输出文件为空: {fasta_file}")
                return False
            
            lines = content.split('\n')
            if not lines[0].startswith('>'):
                self.logger.warning(f"FASTA文件格式错误: {fasta_file}")
                return False
            
            # 检查是否有序列内容
            has_sequences = any(not line.startswith('>') and line.strip() for line in lines)
            if not has_sequences:
                self.logger.warning(f"FASTA文件没有序列内容: {fasta_file}")
                return False
                
            return True
            
        except Exception as e:
            self.logger.warning(f"验证输出文件时出错: {e}")
            return False
    
    def add_sample_prefix(self, fasta_file: Path, sample: str) -> None:
        """在FASTA序列头部添加样本名前缀"""
        try:
            with open(fasta_file, 'r') as f:
                content = f.read()
            
            # 替换序列头
            content = content.replace('>', f'>{sample}_')
            
            with open(fasta_file, 'w') as f:
                f.write(content)
                
        except Exception as e:
            self.logger.warning(f"添加样本名前缀失败: {e}")
    
    def merge_all_split_results(self, args: argparse.Namespace, sample: str) -> Optional[Path]:
        """合并所有工具的split结果"""
        self.logger.info(f"合并所有工具的split结果: {sample}")
        
        # 查找所有split结果文件
        split_files = []
        tools = ['megahit', 'rnaviralspades', 'penguin']
        
        for tool in tools:
            base_prefix = f"{sample}_{tool}_refineC"
            
            # 可能的文件位置
            possible_files = [
                self.sample_output / f"{base_prefix}.split.fasta.gz",
                self.sample_output / f"{base_prefix}.split.fasta",
                self.sample_output / f"{sample}_{tool}.refineC.split.fasta.gz",
                self.sample_output / f"{sample}_{tool}.refineC.split.fasta",
            ]
            
            split_dir = self.sample_output / f"{sample}_{tool}_refineC"
            if split_dir.exists():
                possible_files.extend([
                    split_dir / f"{base_prefix}.split.fasta.gz",
                    split_dir / f"{base_prefix}.split.fasta",
                    split_dir / f"{sample}_{tool}.refineC.split.fasta.gz",
                    split_dir / f"{sample}_{tool}.refineC.split.fasta",
                ])
            
            found_file = None
            for possible_file in possible_files:
                if possible_file.exists() and self.is_fasta_file_valid(possible_file):
                    found_file = possible_file
                    self.logger.info(f"找到 {tool} 的split文件: {found_file}")
                    break
            
            if not found_file:
                split_patterns = [f"*{tool}*split*.fasta*", f"*{tool}*refineC*.fasta*"]
                for pattern in split_patterns:
                    found = list(self.sample_output.glob(pattern))
                    if found:
                        for file in found:
                            if self.is_fasta_file_valid(file):
                                found_file = file
                                self.logger.info(f"通过模式匹配找到 {tool} 的split文件: {found_file}")
                                break
                    if found_file:
                        break
            
            if found_file:
                split_files.append(found_file)
            else:
                self.logger.warning(f"未找到工具 {tool} 的split结果文件")
        
        if not split_files:
            self.logger.warning("未找到任何split结果文件")
            return None
        
        self.logger.info(f"找到 {len(split_files)} 个split结果文件")
        
        # 创建合并文件
        merged_file = self.sample_output / f"{sample}.all_tools.split_merged.fasta"
        
        try:
            with open(merged_file, 'w') as outfile:
                for split_file in split_files:
                    self.logger.info(f"合并文件: {split_file.name}")
                    
                    if split_file.suffix == '.gz':
                        import gzip
                        with gzip.open(split_file, 'rt') as infile:
                            content = infile.read()
                            outfile.write(content)
                    else:
                        with open(split_file, 'r') as infile:
                            content = infile.read()
                            outfile.write(content)
        
            self.logger.info(f"合并完成: {merged_file}")
            
            if not self.is_fasta_file_valid(merged_file):
                self.logger.warning(f"合并后的文件无效或为空: {merged_file}")
                return None
            
            return merged_file
            
        except Exception as e:
            self.logger.error(f"合并split结果失败: {e}")
            return None
    
    def process_single_sample_with_tool(self, args: argparse.Namespace, sample_info: Dict[str, Any], tool: str) -> Tuple[bool, Dict]:
        """使用指定工具处理单个样本"""
        stats = self._init_stats()
        if self._stop_flag:
            return False, stats
            
        sample = sample_info['sample']
        data_type = sample_info['type']
        start_time = time.time()

        # 先设置目录
        self.setup_directories(args, sample, data_type)

        # 检查最终结果是否存在
        final_output_name = self.get_output_filename(tool, sample)
        final_output_path = self.sample_output / final_output_name
        
        if final_output_path.exists() and not args.force:
            self.logger.info(f"最终结果已存在，跳过处理: {final_output_path}")

            # 1. 捞取组装工具旧日志
            old_time_log = self.sample_output / self.get_time_log_filename(tool, sample)
            if old_time_log.exists():
                t_s, m_mb = self.parse_time_mem(old_time_log, 0.0)
                stats[tool]['time'] = t_s
                stats[tool]['mem'] = m_mb

            # 2. 检查或修复 split 资源
            if args.refineC_split:
                split_dir = self.sample_output / f"{sample}_{tool}_refineC"
                if split_dir.exists():
                    split_files = list(split_dir.glob("*.split.fasta*"))
                    if not split_files:
                        self.logger.info(f"运行refineC split处理已有结果: {tool}")
                        split_result = self.run_refineC_split(args, final_output_path, sample, tool, stats)
                        if split_result:
                            self.logger.info(f"refineC split 完成: {split_result}")
                    else:
                        # ─── 核心修复：如果 split 已经跑完，把 split 的历史资源捞回来 ───
                        old_split_log = split_dir / f"{sample}_{tool}.refineC.time.mem.log"
                        if old_split_log.exists():
                            t_s, m_mb = self.parse_time_mem(old_split_log, 0.0)
                            stats[f'refineC_split_{tool}'] = {'time': t_s, 'mem': m_mb}
                else:
                    self.logger.info(f"运行refineC split处理已有结果: {tool}")
                    split_result = self.run_refineC_split(args, final_output_path, sample, tool, stats)
                    if split_result:
                        self.logger.info(f"refineC split 完成: {split_result}")

            self.cleanup(args)

            # 3. 修正 total 逻辑
            stats['total']['time'] = stats[tool]['time'] + (time.time() - start_time)
            stats['total']['mem'] = self._get_max_mem(stats)
            return True, stats

    def process_single_sample_all_tools(self, args: argparse.Namespace, sample_info: Dict[str, Any]) -> Tuple[bool, Dict]:
        """使用所有工具处理单个样本"""
        stats = self._init_stats()
        if self._stop_flag:
            return False, stats
            
        sample = sample_info['sample']
        tools = ['megahit', 'rnaviralspades', 'penguin']
        success_count = 0
        start_time = time.time()
        
        self.setup_directories(args, sample, sample_info['type'])
        
        for tool in tools:
            if self._stop_flag:
                break
                
            self.logger.info(f"使用 {tool} 处理样本 {sample}")
            
            final_output_name = self.get_output_filename(tool, sample)
            final_output_path = self.sample_output / final_output_name
            
            input_files = self.prepare_input_files(sample_info, sample)
            
            if final_output_path.exists() and not args.force:
                self.logger.info(f"{tool} 结果已存在，跳过组装: {final_output_path}")
                success_count += 1

                # ─── 修复 1: 捞取多工具主日志 ───
                old_time_log = self.sample_output / self.get_time_log_filename(tool, sample)
                if old_time_log.exists():
                    t_s, m_mb = self.parse_time_mem(old_time_log, 0.0)
                    stats[tool]['time'] = t_s
                    stats[tool]['mem'] = m_mb

                if args.refineC_split:
                    split_dir = self.sample_output / f"{sample}_{tool}_refineC"
                    if split_dir.exists():
                        split_files = list(split_dir.glob("*.split.fasta*"))
                        if not split_files:
                            self.logger.info(f"运行refineC split处理已有结果: {tool}")
                            split_result = self.run_refineC_split(args, final_output_path, sample, tool, stats)
                            if split_result:
                                self.logger.info(f"refineC split 完成: {split_result}")
                        else:
                            # ─── 修复 2: 捞取多工具 split 历史资源 ───
                            old_split_log = split_dir / f"{sample}_{tool}.refineC.time.mem.log"
                            if old_split_log.exists():
                                t_s, m_mb = self.parse_time_mem(old_split_log, 0.0)
                                stats[f'refineC_split_{tool}'] = {'time': t_s, 'mem': m_mb}
                    else:
                        self.logger.info(f"运行refineC split处理已有结果: {tool}")
                        split_result = self.run_refineC_split(args, final_output_path, sample, tool, stats)
                        if split_result:
                            self.logger.info(f"refineC split 完成: {split_result}")
                continue
                
            if self.run_assembly(args, input_files, sample, tool, stats):
                if self.copy_results(args, sample, tool, stats):
                    success_count += 1
                else:
                    self.logger.error(f"{tool} 复制结果失败")
            else:
                self.logger.error(f"{tool} 组装失败")
        
        if args.refineC_merge and success_count > 0:
            self.logger.info("开始refineC merge处理")
            
            merged_split_file = self.merge_all_split_results(args, sample)
            if merged_split_file:
                merge_result = self.run_refineC_merge(args, merged_split_file, sample, stats)
                if merge_result:
                    self.logger.info(f"refineC merge 完成: {merge_result}")
                else:
                    self.logger.warning("refineC merge 失败或无输出，但样本组装成功")
            else:
                self.logger.warning("无法合并split结果，跳过refineC merge")
        
        self.cleanup(args)
        
        self.logger.info(f"样本 {sample} 的所有工具处理完成: {success_count}/{len(tools)} 个工具成功")
        
        stats['total']['time'] = time.time() - start_time
        stats['total']['mem'] = self._get_max_mem(stats)
        return success_count > 0, stats

    def process_sample_wrapper(self, task_data: Tuple) -> Tuple[str, bool, Dict]:
        """包装样本处理函数，用于并行执行"""
        args, sample_info, tool = task_data
        sample = sample_info['sample']
        
        pipeline = AssemblyPipeline()
        
        try:
            if tool == 'all':
                success, stats = pipeline.process_single_sample_all_tools(args, sample_info)
            else:
                success, stats = pipeline.process_single_sample_with_tool(args, sample_info, tool)
            return sample, success, stats
        except Exception as e:
            pipeline.logger.error(f"处理样本 {sample} 时发生错误: {e}")
            return sample, False, pipeline._init_stats()
    
    def run_parallel_processing(self, args: argparse.Namespace, sample_info: Dict[str, Any]) -> bool:
        """并行处理多个样本"""
        samples = sample_info['samples']
        total_samples = len(samples)
        
        self.logger.info(f"开始并行处理 {total_samples} 个样本，并行任务数: {args.jobs}")
        
        tasks = []
        for sample, files in samples.items():
            sample_data = {
                'type': files['type'],
                'sample': sample,
            }
            
            if files['type'] == 'paired':
                sample_data['read1'] = files['read1']
                sample_data['read2'] = files['read2']
            else:
                sample_data['reads'] = files['reads']
            
            tasks.append((args, sample_data, args.tool))
        
        success_count = 0
        failed_samples = []
        
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_sample = {
                executor.submit(self.process_sample_wrapper, task): task[1]['sample'] 
                for task in tasks
            }
            
            with tqdm(total=total_samples, desc="处理样本", unit="sample") as pbar:
                for future in as_completed(future_to_sample):
                    sample = future_to_sample[future]
                    try:
                        result_sample, success, stats = future.result()
                        if success:
                            success_count += 1
                            pbar.set_description(f"完成: {sample} (成功)")
                        else:
                            failed_samples.append(sample)
                            pbar.set_description(f"完成: {sample} (失败)")
                        
                        self.write_stats_log(args, sample, stats)
                        pbar.update(1)
                        
                        if self._stop_flag:
                            self.logger.info("接收到停止信号，取消剩余任务...")
                            executor.shutdown(wait=False)
                            for f in future_to_sample:
                                f.cancel()
                            break
                            
                    except Exception as e:
                        self.logger.error(f"处理样本 {sample} 时发生异常: {e}")
                        failed_samples.append(sample)
                        pbar.update(1)
        
        self.logger.info(f"并行处理完成: {success_count}/{total_samples} 个样本成功")
        if failed_samples:
            self.logger.warning(f"失败的样本: {', '.join(failed_samples)}")
        
        return success_count == total_samples
    
    def cleanup(self, args: argparse.Namespace) -> None:
        """清理临时文件和中间目录"""
        if not args.keep_temp and self.temp_dir and self.temp_dir.exists():
            try:
                shutil.rmtree(self.temp_dir)
                self.logger.info("临时文件已清理")
            except Exception as e:
                self.logger.warning(f"清理临时文件失败: {e}")
                
        # 自动清理 refineC 产生的中间子目录
        if not args.keep_temp and self.sample_output and self.sample_output.exists():
            try:
                # 匹配所有带 _refineC 的目录 (例如 *_megahit_refineC, *_all_tools_refineC_merge 等)
                for subdir in self.sample_output.glob('*_refineC*'):
                    if subdir.is_dir():
                        shutil.rmtree(subdir)
                        self.logger.info(f"已自动清理 refineC 中间目录: {subdir.name}")
            except Exception as e:
                self.logger.warning(f"清理 refineC 中间目录失败: {e}")
    
    def run_pipeline(self) -> int:
        """运行完整的组装流水线"""
        args = self.parse_arguments()
        
        self.setup_signal_handlers()
        
        if not self.validate_arguments(args):
            return 1
        
        try:
            sample_info = self.detect_data_type_and_sample(args.input)
            
            if sample_info['type'] == 'batch':
                if args.jobs > 1:
                    success = self.run_parallel_processing(args, sample_info)
                    return 0 if success else 1
                else:
                    success_count = 0
                    total_count = len(sample_info['samples'])
                    
                    for sample, files in tqdm(sample_info['samples'].items(), desc="处理样本", unit="sample"):
                        if self._stop_flag:
                            break
                            
                        sample_data = {
                            'type': files['type'],
                            'sample': sample,
                        }
                        
                        if files['type'] == 'paired':
                            sample_data['read1'] = files['read1']
                            sample_data['read2'] = files['read2']
                        else:
                            sample_data['reads'] = files['reads']
                        
                        if args.tool == 'all':
                            success, stats = self.process_single_sample_all_tools(args, sample_data)
                            if success:
                                success_count += 1
                            else:
                                self.logger.error(f"样本 {sample} 的部分工具处理失败")
                        else:
                            success, stats = self.process_single_sample_with_tool(args, sample_data, args.tool)
                            if success:
                                success_count += 1
                            else:
                                self.logger.error(f"样本 {sample} 处理失败")
                            
                        self.write_stats_log(args, sample, stats)
                    
                    self.logger.info(f"批量处理完成: {success_count}/{total_count} 个样本成功")
                    return 0 if success_count == total_count else 1
                    
            else:
                if args.tool == 'all':
                    success, stats = self.process_single_sample_all_tools(args, sample_info)
                else:
                    success, stats = self.process_single_sample_with_tool(args, sample_info, args.tool)
                
                self.write_stats_log(args, sample_info['sample'], stats)
                return 0 if success else 1
                
        except KeyboardInterrupt:
            self.logger.info("用户中断执行")
            return 1
        except Exception as e:
            self.logger.error(f"流水线执行失败: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return 1

def main():
    """主函数"""
    pipeline = AssemblyPipeline()
    sys.exit(pipeline.run_pipeline())

if __name__ == '__main__':
    main()
