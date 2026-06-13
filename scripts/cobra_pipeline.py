#!/usr/bin/env python3
"""
COBRA病毒序列扩展分析脚本 (Python版) - 批量并行版本

功能特点：
1. 支持多种组装工具和病毒鉴定工具的批量处理
2. 支持自动样本发现模式
3. 支持断点运行和CheckV模式选择
4. 支持COBRA结果检测和重命名
5. 支持virus.fa空文件检测和智能跳过
6. 支持多模式并行处理
"""

import os
import sys
import re
import argparse
import subprocess
import shutil
from pathlib import Path
import logging
from typing import List, Dict, Tuple, Optional, Set
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from tqdm import tqdm
import time
import json
from datetime import datetime
from collections import defaultdict

# ==================== 配置部分 ====================

# 设置主日志记录器
def setup_main_logger():
    """设置主日志记录器"""
    logger = logging.getLogger("cobra_main")
    logger.setLevel(logging.INFO)
    
    # 清除现有处理器
    logger.handlers.clear()
    
    # 控制台处理器 - 只显示重要信息
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
    console_handler.setFormatter(console_formatter)
    
    # 文件处理器 - 记录所有信息
    file_handler = logging.FileHandler('cobra_pipeline.log')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    
    return logger

# 初始化主日志
main_logger = setup_main_logger()

# 组装工具映射
ASSEMBLER_MAP = {
    "megahit": "megahit",
    "rnaviralspades": "rnaviralspades", 
    "penguin": "penguin"
}

# 病毒鉴定工具映射
VIRUS_IDENTIFICATION_MAP = {
    "megahit": "megahit_result",
    "rnaviralspades": "rnaviralspades_result",
    "penguin": "penguin_result"
}

# 模式配置
MODE_CONFIG = {
    "virus": {
        "reads_subdir": "1.virus-reads",
        "contigs_subdir": "1.virus-assembly",
        "virsorter_subdir": "virus",
        "suffix": ".unmapped.virus"
    },
    "other": {
        "reads_subdir": "2.other-reads",
        "contigs_subdir": "2.other-assembly",
        "virsorter_subdir": "other",
        "suffix": ".unmapped.other"
    },
    "mix": {
        "reads_subdir": "3.mix-reads",
        "contigs_subdir": "3.mix-assembly",
        "virsorter_subdir": "mix",
        "suffix": ".unmapped.mix"
    }
}

# ==================== COBRA管道类 ====================

class CobraPipeline:
    def __init__(self, args):
        """初始化COBRA管道"""
        self.args = args
        self.start_time = time.time()
        
        # 初始化断点状态
        self.checkpoint_file = Path(self.args.output_dir) / "checkpoint_status.json"
        self.completed_tasks = self.load_checkpoint()
        
        # 根据模式自动设置目录路径
        self.setup_directories_by_mode()
        
        # 解析样本列表
        self.samples = self.parse_samples(args.samples)
        
        # 解析组装工具列表
        self.assembly_tools = self.parse_assembly_tools(args.assembly_tools)
        
        # 显示配置信息
        self.display_configuration()
        
        # 检查必要工具是否可用
        self.check_tools()
        
        # 验证输入目录
        self.validate_directories()
    
    def display_configuration(self):
        """显示运行配置"""
        config_lines = [
            "=" * 60,
            "COBRA批量分析流程配置",
            "=" * 60,
            f"运行模式: {self.args.mode}",
            f"样本数量: {len(self.samples)}",
            f"组装工具: {', '.join(self.assembly_tools)}",
            f"并行任务数: {self.args.jobs}",
            f"单任务线程数: {self.args.threads}",
            f"断点运行: {'是' if self.args.resume else '否'}",
            f"Reads目录: {self.args.reads_dir}",
            f"Contigs目录: {self.args.contigs_dir}",
            f"病毒鉴定目录: {self.args.virsorter_dir}",
            f"输出目录: {self.args.output_dir}",
            f"COBRA参数: mink={self.args.mink}, maxk={self.args.maxk}, linkage_mismatch={self.args.linkage_mismatch}",
            "=" * 60
        ]
        
        for line in config_lines:
            if line.startswith("="):
                main_logger.info(line)
            else:
                main_logger.info(f"  {line}")
    
    def load_checkpoint(self):
        """加载断点状态"""
        if self.args.resume and self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                main_logger.warning(f"无法读取断点文件: {e}")
                return {}
        return {}
    
    def save_checkpoint(self, task_id):
        """保存断点状态"""
        if task_id not in self.completed_tasks:
            self.completed_tasks[task_id] = {
                "status": "completed",
                "timestamp": datetime.now().isoformat()
            }
            try:
                with open(self.checkpoint_file, 'w') as f:
                    json.dump(self.completed_tasks, f, indent=2)
            except Exception as e:
                main_logger.warning(f"无法保存断点文件: {e}")
    
    def is_task_completed(self, task_id):
        """检查任务是否已完成"""
        if not self.args.resume:
            return False
        return task_id in self.completed_tasks
    
    def setup_directories_by_mode(self):
        """根据模式设置目录路径"""
        mode = self.args.mode
        
        if mode == "all":
            main_logger.info("使用用户指定的目录路径")
            self.args.output_base = Path(self.args.output_dir)
            self.args.output_dir = str(self.args.output_base)
            return
        
        # 根据模式更新默认路径
        config = MODE_CONFIG[mode]
        
        # 如果用户使用的是默认路径，则更新为模式特定的路径
        if self.args.reads_dir == "7.Virus.reads/3.mix-reads/":
            self.args.reads_dir = f"7.Virus.reads/{config['reads_subdir']}/"
            main_logger.info(f"更新reads目录为: {self.args.reads_dir}")
        
        if self.args.contigs_dir == "8.virus_assembly/3.mix-assembly/":
            self.args.contigs_dir = f"8.virus_assembly/{config['contigs_subdir']}/"
            main_logger.info(f"更新contigs目录为: {self.args.contigs_dir}")
        
        if self.args.virsorter_dir == "9.virus_identification/mix/":
            self.args.virsorter_dir = f"9.virus_identification/{config['virsorter_subdir']}/"
            main_logger.info(f"更新virsorter目录为: {self.args.virsorter_dir}")
        
        # 设置输出目录为模式特定子目录的父目录
        self.args.output_base = Path(self.args.output_dir)
        self.args.output_dir = str(self.args.output_base)
    
    def validate_directories(self):
        """验证输入目录是否存在"""
        main_logger.info("验证输入目录...")
        
        dir_checks = [
            ("Reads目录", self.args.reads_dir),
            ("Contigs目录", self.args.contigs_dir),
            ("病毒鉴定目录", self.args.virsorter_dir),
        ]
        
        all_valid = True
        for dir_name, dir_path in dir_checks:
            path = Path(dir_path)
            if not path.exists():
                main_logger.error(f"✗ {dir_name}不存在: {dir_path}")
                all_valid = False
            else:
                main_logger.info(f"✓ {dir_name}存在: {dir_path}")
        
        if not all_valid:
            sys.exit(1)
    
    def parse_samples(self, samples_input: Optional[str]) -> List[str]:
        """解析样本列表"""
        if samples_input:
            # 如果提供了样本输入，按原逻辑处理
            if os.path.exists(samples_input):
                with open(samples_input, 'r') as f:
                    samples = [line.strip() for line in f if line.strip()]
                main_logger.info(f"从文件读取 {len(samples)} 个样本: {samples_input}")
            else:
                samples = [s.strip() for s in samples_input.split(',') if s.strip()]
                main_logger.info(f"从参数读取 {len(samples)} 个样本")
        else:
            samples = self.discover_samples_by_mode()
        
        return samples
    
    def discover_samples_by_mode(self) -> List[str]:
        """根据模式自动发现样本 (支持通用命名回退)"""
        mode = self.args.mode
        reads_dir = Path(self.args.reads_dir)

        if not reads_dir.exists():
            main_logger.error(f"Reads目录不存在: {reads_dir}")
            sys.exit(1)

        samples = set()

        def _collect(scan_dir, suffix):
            """从目录收集样本名"""
            for read_file in scan_dir.glob(f"*{suffix}.R1.*"):
                samples.add(read_file.name.split(".R1.")[0])
            for read_file in scan_dir.glob(f"*{suffix}_1.*"):
                samples.add(read_file.name.split("_1.")[0])
            for read_file in scan_dir.glob(f"*{suffix}_R1.*"):
                samples.add(read_file.name.split("_R1.")[0])

        if mode == "all":
            for sub_mode in ["virus", "other", "mix"]:
                config = MODE_CONFIG[sub_mode]
                suffix = config["suffix"]
                mode_reads_dir = reads_dir.parent / config["reads_subdir"]
                if mode_reads_dir.exists():
                    _collect(mode_reads_dir, suffix)
        else:
            config = MODE_CONFIG[mode]
            suffix = config["suffix"]
            _collect(reads_dir, suffix)

        # 回退: 如果没找到任何样本，尝试不带 suffix 的通用模式
        if not samples:
            main_logger.info("模式匹配无结果，尝试通用命名模式...")
            _collect(reads_dir, "")

        samples_list = sorted(list(samples))

        if not samples_list:
            main_logger.error(f"在目录 {reads_dir} 中未找到任何样本")
            sys.exit(1)

        main_logger.info(f"发现 {len(samples_list)} 个样本")
        return samples_list
    
    def parse_assembly_tools(self, tools_input: str) -> List[str]:
        """解析组装工具列表"""
        if tools_input.lower() == 'all':
            return list(ASSEMBLER_MAP.keys())
        else:
            tools = [t.strip() for t in tools_input.split(',') if t.strip()]
            valid_tools = []
            for tool in tools:
                if tool in ASSEMBLER_MAP:
                    valid_tools.append(tool)
                else:
                    main_logger.warning(f"忽略无效的组装工具: {tool}")
            return valid_tools
    
    def check_tools(self):
        """检查必要的工具是否在PATH中"""
        tools = ['bwa-mem2', 'samtools', 'coverm']
        missing_tools = []
        
        for tool in tools:
            if shutil.which(tool) is None:
                missing_tools.append(tool)
        
        if missing_tools:
            main_logger.error(f"以下工具未找到: {', '.join(missing_tools)}")
            sys.exit(1)
        else:
            main_logger.info("所有必要工具已安装")
    
    def setup_sample_logger(self, sample_name: str, assembler: str, log_file: Path) -> logging.Logger:
        """为样本设置独立的日志记录器"""
        task_id = f"{sample_name}_{assembler}"
        logger_name = f"cobra_{task_id}"
        
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
        file_handler.setFormatter(file_formatter)
        
        logger.addHandler(file_handler)
        
        return logger
    
    def run_command(self, cmd: List[str], desc: str, sample_logger: logging.Logger, 
                   shell: bool = False):
        """运行shell命令并处理错误"""
        sample_logger.info(f"运行: {desc}")
        sample_logger.debug(f"命令: {' '.join(cmd)}")
        
        try:
            if shell:
                cmd_str = cmd[0] if len(cmd) == 1 else ' '.join(cmd)
                result = subprocess.run(
                    cmd_str,
                    shell=True,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True
                )
            else:
                result = subprocess.run(
                    cmd,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True
                )
            
            if result.stdout:
                sample_logger.debug(f"命令输出:\n{result.stdout}")
            
            return result
            
        except subprocess.CalledProcessError as e:
            sample_logger.error(f"命令执行失败: {desc}")
            sample_logger.error(f"错误信息: {e.stderr}")
            raise
    
    def check_fasta_file(self, fasta_file: Path) -> bool:
        """检查FASTA文件是否包含有效序列"""
        if not fasta_file.exists():
            return False
        
        file_size = fasta_file.stat().st_size
        if file_size < 100:
            return False
        
        try:
            with open(fasta_file, 'r') as f:
                content = f.read()
                
                if '>' not in content:
                    return False
                
                sequences = []
                current_seq = ""
                in_sequence = False
                
                for line in content.split('\n'):
                    line = line.strip()
                    if line.startswith('>'):
                        if in_sequence and current_seq:
                            sequences.append(current_seq)
                            current_seq = ""
                        in_sequence = True
                    elif in_sequence and line:
                        current_seq += line
                
                if in_sequence and current_seq:
                    sequences.append(current_seq)
                
                if not sequences:
                    return False
                
                total_length = sum(len(seq) for seq in sequences)
                if total_length < 100:
                    return False
                
                return True
                
        except Exception:
            return False
    
    def normalize_contig_names(self, input_file: Path, output_file: Path, sample_logger: logging.Logger):
        """标准化contig名称"""
        sample_logger.info(f"标准化contig名称: {input_file} -> {output_file}")
        
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    parts = line.split()
                    contig_id = parts[0][1:]
                    if '_' in contig_id:
                        fout.write(f'>{contig_id}\n')
                    else:
                        fout.write(line)
                else:
                    fout.write(line)
    
    def normalize_virus_names(self, input_file: Path, output_file: Path, sample_logger: logging.Logger):
        """标准化病毒序列名称"""
        sample_logger.info(f"标准化病毒序列名称: {input_file} -> {output_file}")
        
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    contig_id = line.split(' ')[0]
                    if '_' in contig_id[1:]:
                        fout.write(f'{contig_id}\n')
                    else:
                        fout.write(line)
                else:
                    fout.write(line)
    
    def get_sample_mode(self, sample: str) -> str:
        """根据样本名获取模式"""
        for mode, config in MODE_CONFIG.items():
            if sample.endswith(config["suffix"]):
                return mode
        return self.args.mode if self.args.mode != "all" else "mix"
    
    def find_contig_file(self, sample: str, assembler: str) -> Path:
        """查找contig文件"""
        sample_mode = self.get_sample_mode(sample)
        
        if self.args.mode == "all":
            config = MODE_CONFIG[sample_mode]
            contigs_dir = Path(self.args.contigs_dir).parent / config["contigs_subdir"]
        else:
            contigs_dir = Path(self.args.contigs_dir)
        
        sample_assembly_dir = contigs_dir / f"{sample}"
        
        if not sample_assembly_dir.exists():
            raise FileNotFoundError(f"样本目录不存在: {sample_assembly_dir}")
        
        contig_patterns = [
            f"{sample}_{ASSEMBLER_MAP[assembler]}.contig.fasta",
            f"{sample}.{ASSEMBLER_MAP[assembler]}.contig.fasta",
            f"final.contigs.fa",
        ]
        
        for pattern in contig_patterns:
            contig_file = sample_assembly_dir / pattern
            if contig_file.exists():
                return contig_file
        
        for file_path in sample_assembly_dir.rglob("*"):
            if ASSEMBLER_MAP[assembler] in file_path.name and file_path.suffix in ['.fa', '.fasta', '.fna']:
                return file_path
        
        raise FileNotFoundError(f"未找到样本 {sample} 的 {assembler} contig文件")
    
    def find_virus_file(self, sample: str, assembler: str) -> Path:
        """查找病毒序列文件 (支持嵌套 {tool}_result/ 和扁平目录)"""
        sample_mode = self.get_sample_mode(sample)

        if self.args.mode == "all":
            config = MODE_CONFIG[sample_mode]
            virsorter_dir = Path(self.args.virsorter_dir).parent / config["virsorter_subdir"]
        else:
            virsorter_dir = Path(self.args.virsorter_dir)

        base_sample = sample
        for mode_config in MODE_CONFIG.values():
            if sample.endswith(mode_config["suffix"]):
                base_sample = sample.replace(mode_config["suffix"], '')
                break

        # 样本名变体 (处理 _clean 差异)
        name_variants = [base_sample]
        if base_sample.endswith('_clean'):
            name_variants.append(base_sample[:-6])
        elif '_clean' not in base_sample:
            name_variants.append(base_sample + '_clean')

        virus_patterns = [
            f"{n}_virus.all.candidate.fasta" for n in name_variants
        ] + [
            f"{n}_virus.fasta" for n in name_variants
        ] + [
            "final-virus-combined.fa",
        ]

        # 1. 先尝试嵌套结构: {virsorter_dir}/{name}/{tool}_result/
        for name in name_variants:
            virus_ident_dir = virsorter_dir / name / VIRUS_IDENTIFICATION_MAP[assembler]
            if virus_ident_dir.exists():
                for pattern in virus_patterns:
                    for vf in virus_ident_dir.glob(pattern):
                        if vf.exists():
                            return vf
                # 模糊搜索
                for fp in virus_ident_dir.rglob("*"):
                    if fp.suffix in ['.fa', '.fasta', '.fna'] and 'virus' in fp.name.lower():
                        return fp

        # 2. 回退到扁平结构: {virsorter_dir}/{name}/
        for name in name_variants:
            flat_dir = virsorter_dir / name
            if flat_dir.exists():
                for pattern in virus_patterns:
                    for vf in flat_dir.glob(pattern):
                        if vf.exists():
                            return vf
                for fp in flat_dir.rglob("*"):
                    if fp.suffix in ['.fa', '.fasta', '.fna'] and 'virus' in fp.name.lower():
                        return fp

        raise FileNotFoundError(
            f"未找到样本 {base_sample} 的 {assembler} 病毒序列文件 "
            f"(已尝试: {name_variants})")
    
    def find_read_files(self, sample: str) -> Tuple[Path, Path]:
        """查找reads文件 (支持 FASTQ/FASTA 及 _clean_ 命名)"""
        sample_mode = self.get_sample_mode(sample)

        if self.args.mode == "all":
            config = MODE_CONFIG[sample_mode]
            reads_dir = Path(self.args.reads_dir).parent / config["reads_subdir"]
        else:
            reads_dir = Path(self.args.reads_dir)

        # 同时尝试原始样本名和去 _clean 后缀的变体
        sample_variants = [sample]
        if sample.endswith('_clean'):
            sample_variants.append(sample[:-6])  # 去掉 _clean
        elif '_clean_' in sample:
            sample_variants.append(sample.replace('_clean', ''))

        read_patterns = [
            # FASTQ
            (f"{s}.R1.fq.gz", f"{s}.R2.fq.gz")
            for s in sample_variants
        ] + [
            (f"{s}_1.fastq.gz", f"{s}_2.fastq.gz")
            for s in sample_variants
        ] + [
            (f"{s}_1.fq.gz", f"{s}_2.fq.gz")
            for s in sample_variants
        ] + [
            (f"{s}_1.fastq", f"{s}_2.fastq")
            for s in sample_variants
        ] + [
            (f"{s}.R1.fastq.gz", f"{s}.R2.fastq.gz")
            for s in sample_variants
        ] + [
            # FASTA
            (f"{s}_1.fa.gz", f"{s}_2.fa.gz")
            for s in sample_variants
        ] + [
            (f"{s}_1.fasta.gz", f"{s}_2.fasta.gz")
            for s in sample_variants
        ] + [
            (f"{s}_1.fa", f"{s}_2.fa")
            for s in sample_variants
        ] + [
            (f"{s}.R1.fa.gz", f"{s}.R2.fa.gz")
            for s in sample_variants
        ] + [
            (f"{s}.R1.fasta.gz", f"{s}.R2.fasta.gz")
            for s in sample_variants
        ]

        for r1_pattern, r2_pattern in read_patterns:
            r1_file = reads_dir / r1_pattern
            r2_file = reads_dir / r2_pattern
            if r1_file.exists() and r2_file.exists():
                return r1_file, r2_file

        raise FileNotFoundError(f"未找到样本 {sample} 的reads文件 (已尝试变体: {sample_variants})")
    
    def process_coverage_file(self, coverage_file: Path, output_file: Path, sample_logger: logging.Logger):
        """处理CoverM生成的文件，提取第一列和第三列"""
        sample_logger.info(f"处理CoverM文件: {coverage_file} -> {output_file}")

        with open(coverage_file, 'r') as fin, open(output_file, 'w') as fout:
            lines = fin.readlines()
            for line in lines[1:]:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    fout.write(f"{parts[0]}\t{parts[2]}\n")
                else:
                    sample_logger.warning(f"CoverM文件行格式异常: {line}")
    def process_cobra_results(self, cobra_output: Path, output_dir: Path, base_sample: str, 
                            sample_mode: str, assembler: str, sample_logger: logging.Logger) -> Optional[Path]:
        """处理COBRA结果文件"""
        sample_logger.info("检查COBRA输出结果...")
        
        # 查找COBRA输出的所有category文件
        category_files = list(cobra_output.glob("COBRA_category_*.fasta"))
        
        if category_files:
            sample_logger.info(f"找到 {len(category_files)} 个COBRA分类文件")
            for file in category_files:
                sample_logger.debug(f"  - {file.name}")
            
            # 合并所有category文件
            new_filename = f"{base_sample}.{sample_mode}.{assembler}.cobra.fa"
            new_filepath = output_dir / new_filename
            
            try:
                with open(new_filepath, 'w') as outfile:
                    for category_file in sorted(category_files):
                        with open(category_file, 'r') as infile:
                            outfile.write(infile.read())
                
                # 检查合并后的文件是否有效
                if not self.check_fasta_file(new_filepath):
                    sample_logger.warning(f"合并后的COBRA结果文件无效或为空: {new_filepath}")
                    return None
                
                sample_logger.info(f"成功合并COBRA结果文件: {new_filepath}")
                return new_filepath
                
            except Exception as e:
                sample_logger.error(f"合并COBRA结果文件失败: {e}")
                return None
        else:
            # 如果找不到category文件，回退到查找contigs.new.fa
            sample_logger.info("未找到COBRA分类文件，尝试查找contigs.new.fa")
            
            cobra_result_patterns = [
                cobra_output / f"{base_sample}.{sample_mode}.{assembler}.contigs.new.fa",
                cobra_output / "*.contigs.new.fa",
                cobra_output / "contigs.new.fa",
                cobra_output / "*.new.fa",
            ]
            
            cobra_result_file = None
            for pattern in cobra_result_patterns:
                if '*' in str(pattern):
                    for file in cobra_output.parent.glob(pattern.name):
                        if file.exists():
                            cobra_result_file = file
                            break
                elif pattern.exists():
                    cobra_result_file = pattern
                    break
                
                if cobra_result_file:
                    break
            
            if cobra_result_file and cobra_result_file.exists():
                if not self.check_fasta_file(cobra_result_file):
                    sample_logger.warning(f"COBRA结果文件无效或为空: {cobra_result_file}")
                    return None
                
                new_filename = f"{base_sample}.{sample_mode}.{assembler}.cobra.fa"
                new_filepath = output_dir / new_filename
                
                sample_logger.info(f"发现COBRA延伸结果: {cobra_result_file}")
                sample_logger.info(f"复制为: {new_filename}")
                
                try:
                    shutil.copy2(cobra_result_file, new_filepath)
                    sample_logger.info(f"COBRA结果文件已保存到: {new_filepath}")
                    return new_filepath
                except Exception as e:
                    sample_logger.error(f"复制COBRA结果文件失败: {e}")
                    return None
            else:
                sample_logger.info("COBRA未生成有效结果文件，可能是由于序列太少延伸失败")
                return None
    def process_sample_tool(self, sample: str, assembler: str) -> Dict:
        """处理单个样本的单个组装工具"""
        sample_mode = self.get_sample_mode(sample)
        base_sample = sample
        for mode_config in MODE_CONFIG.values():
            if sample.endswith(mode_config["suffix"]):
                base_sample = sample.replace(mode_config["suffix"], '')
                break
        
        task_id = f"{base_sample}.{sample_mode}.{assembler}"
        
        if self.is_task_completed(task_id):
            main_logger.debug(f"跳过已处理的任务: {task_id}")
            return {
                "task_id": task_id,
                "sample": sample,
                "assembler": assembler,
                "status": "skipped",
                "output_dir": "已跳过",
                "skip_reason": "断点运行已处理"
            }
        
        main_logger.info(f"开始处理: {task_id}")
        
        try:
            if self.args.mode == "all":
                output_base = Path(self.args.output_dir) / sample_mode / base_sample
            else:
                output_base = Path(self.args.output_dir) / base_sample
            
            if assembler == "megahit":
                output_dir = output_base / "cobra_megahit_result"
            elif assembler == "rnaviralspades":
                output_dir = output_base / "cobra_rnaviralspades_result"
            elif assembler == "penguin":
                output_dir = output_base / "cobra_penguin_result"
            else:
                output_dir = output_base / f"cobra_{assembler}_result"
            
            output_dir.mkdir(parents=True, exist_ok=True)
            
            task_log_file = output_dir / f"{task_id}.log"
            sample_logger = self.setup_sample_logger(base_sample, assembler, task_log_file)
            sample_logger.info(f"开始处理任务: {task_id}")
            
            # 1. 查找输入文件
            sample_logger.info("查找输入文件...")
            contigs_file = self.find_contig_file(sample, assembler)
            virus_file = self.find_virus_file(sample, assembler)
            read1, read2 = self.find_read_files(sample)
            sample_logger.info(f"找到contig文件: {contigs_file}")
            sample_logger.info(f"找到病毒文件: {virus_file}")
            sample_logger.info(f"找到reads文件: {read1}, {read2}")
            
            # 2. 标准化contigs和病毒序列名称
            normalized_contigs = output_dir / f"{task_id}.contigs.fa"
            normalized_virus = output_dir / f"{task_id}.virus.fa"
            
            self.normalize_contig_names(contigs_file, normalized_contigs, sample_logger)
            self.normalize_virus_names(virus_file, normalized_virus, sample_logger)
            
            # 3. 检查virus.fa是否为空或无效
            sample_logger.info("检查virus.fa文件有效性...")
            if not self.check_fasta_file(normalized_virus):
                sample_logger.warning(f"virus.fa文件为空或无效，跳过COBRA运行")

                sample_logger.info(f"任务完成（跳过COBRA）: {task_id}")
                main_logger.info(f"完成处理（跳过COBRA）: {task_id} [virus.fa为空或无效]")

                return {
                    "task_id": task_id,
                    "sample": sample,
                    "assembler": assembler,
                    "status": "success",
                    "cobra_status": "skipped",
                    "skip_reason": "virus.fa为空或无效",
                    "output_dir": str(output_dir),
                    "cobra_file": None
                }
            
            # 4. 构建bwa-mem2索引
            sample_logger.info("构建bwa-mem2索引")
            self.run_command(
                ["bwa-mem2", "index", str(normalized_contigs)],
                "构建bwa-mem2索引",
                sample_logger
            )
            
            # 5. 使用管道进行比对并生成排序BAM
            sample_logger.info("使用管道进行比对并生成排序BAM")
            sorted_bam = output_dir / f"{task_id}.sorted.bam"

            pipe_cmd = (
                f"bwa-mem2 mem -t {self.args.threads} "
                f"{normalized_contigs} {read1} {read2} | "
                f"samtools sort -@ {self.args.threads} -o {sorted_bam} -"
            )
            
            self.run_command(
                [pipe_cmd],
                "bwa比对和samtools排序",
                sample_logger,
                shell=True
            )
            
            # 6. 创建BAM索引
            sample_logger.info("创建BAM索引")
            self.run_command(
                ["samtools", "index", str(sorted_bam)],
                "索引BAM",
                sample_logger
            )
            
            # 7. 计算覆盖率
            sample_logger.info("计算覆盖率")
            coverage_file_raw = output_dir / f"{task_id}.CoverM.txt"
            coverage_file = output_dir / f"{task_id}.coverage.txt"
            
            self.run_command([
                "coverm", "contig",
                "-b", str(sorted_bam),
                "-t", str(self.args.threads),
                "-m", "covered_fraction", "mean", "rpkm",
                "--output-file", str(coverage_file_raw),
                "--min-covered-fraction", "0"
            ], "计算覆盖率", sample_logger)
            
            # 8. 处理CoverM文件
            self.process_coverage_file(coverage_file_raw, coverage_file, sample_logger)
            
            # 9. 运行COBRA
            sample_logger.info("运行COBRA")
            
            assembler_map_for_cobra = {
                "megahit": "megahit",
                "rnaviralspades": "metaspades",
                "penguin": "megahit"
            }
            
            cobra_assembler = assembler_map_for_cobra.get(assembler, "megahit")
            
            cobra_output_name = f"{base_sample}.{sample_mode}.{assembler}.COBRA"
            
            cobra_output = output_dir / cobra_output_name
            cobra_output.parent.mkdir(parents=True, exist_ok=True)
            
            self.run_command([
                "cobra-meta",
                "-f", str(normalized_contigs),
                "-q", str(normalized_virus),
                "-o", str(cobra_output),
                "-c", str(coverage_file),
                "-m", str(sorted_bam),
                "-a", cobra_assembler,
                "-mink", str(self.args.mink),
                "-maxk", str(self.args.maxk),
                "-lm", str(self.args.linkage_mismatch)
            ], "运行COBRA", sample_logger)
            
            # 10. 处理COBRA结果
            cobra_final_file = self.process_cobra_results(
                cobra_output, output_dir, base_sample, sample_mode, assembler, sample_logger
            )

            sample_logger.info(f"任务完成: {task_id}")

            cobra_status = "成功生成延伸序列" if cobra_final_file else "未生成有效延伸序列"
            main_logger.info(f"完成处理: {task_id} [{cobra_status}]")

            return {
                "task_id": task_id,
                "sample": sample,
                "assembler": assembler,
                "status": "success",
                "cobra_status": "success" if cobra_final_file else "no_extension",
                "output_dir": str(output_dir),
                "cobra_file": str(cobra_final_file) if cobra_final_file else None
            }
            
        except Exception as e:
            error_msg = f"处理失败 {task_id}: {str(e)}"
            main_logger.error(error_msg)
            return {
                "task_id": task_id,
                "sample": sample,
                "assembler": assembler,
                "status": "failed",
                "error": str(e)
            }
    
    def run(self):
        """运行整个流程"""
        main_logger.info("=" * 60)
        main_logger.info("开始COBRA批量分析流程")
        main_logger.info("=" * 60)

        tasks = []
        for sample in self.samples:
            for assembler in self.assembly_tools:
                tasks.append((sample, assembler))

        main_logger.info(f"总任务数: {len(tasks)}")

        if tasks and not self.args.resume:
            test_sample, test_assembler = tasks[0]
            main_logger.info(f"测试第一个任务: {test_sample}_{test_assembler}")
            try:
                contigs_file = self.find_contig_file(test_sample, test_assembler)
                virus_file = self.find_virus_file(test_sample, test_assembler)
                read1, read2 = self.find_read_files(test_sample)

                if self.check_fasta_file(virus_file):
                    main_logger.info(f"测试成功，所有文件都能找到，virus.fa有效")
                else:
                    main_logger.warning(f"测试成功，但virus.fa为空或无效，将跳过COBRA运行")

            except Exception as e:
                main_logger.error(f"测试失败: {e}")
                main_logger.error("请检查以上错误并修复后再运行完整流程")
                sys.exit(1)

        results = []

        # 设置tqdm样式 - 添加颜色支持
        tqdm_kwargs = {
            "total": len(tasks),
            "desc": "🚀 处理进度",
            "unit": "任务",
            "bar_format": "{l_bar}{bar:40}{r_bar}",
            "colour": "green",  # 使用青色进度条
            "dynamic_ncols": True,  # 动态调整宽度
            "ascii": False,  # 使用Unicode字符
        }

        with ProcessPoolExecutor(max_workers=self.args.jobs) as executor:
            future_to_task = {
                executor.submit(self.process_sample_tool, sample, assembler): (sample, assembler)
                for sample, assembler in tasks
            }

            # 使用带颜色的进度条
            with tqdm(**tqdm_kwargs) as pbar:
                for future in as_completed(future_to_task):
                    sample, assembler = future_to_task[future]

                    try:
                        result = future.result(timeout=3600)
                        results.append(result)

                        # 主进程保存断点 (避免子进程竞态覆盖)
                        tid = result.get("task_id")
                        if tid and result.get("status") in ("success", "skipped"):
                            self.save_checkpoint(tid)

                        # 根据任务状态设置不同的颜色和图标
                        status_info = []
                        status_color = ""
                        if result.get("status") == "success":
                            cobra_status = result.get("cobra_status", "")
                            if cobra_status == "skipped":
                                status_info.append("✅ 成功(跳过COBRA)")
                                status_color = "yellow"
                            elif cobra_status == "no_extension":
                                status_info.append("⚠️ 成功(无延伸)")
                                status_color = "magenta"
                            else:
                                status_info.append("✅ 成功")
                                status_color = "green"
                        elif result.get("status") == "skipped":
                            status_info.append("⏭️ 跳过")
                            status_color = "blue"
                        else:
                            status_info.append("❌ 失败")
                            status_color = "red"

                        # 更新进度条后置信息
                        pbar.set_postfix_str(
                            f"{sample}_{assembler} ({' '.join(status_info)})",
                            refresh=False
                        )

                    except Exception as e:
                        results.append({
                            "sample": sample,
                            "assembler": assembler,
                            "status": "exception",
                            "error": str(e)
                        })
                        pbar.set_postfix_str(
                            f"{sample}_{assembler} (🔥 异常)",
                            refresh=False
                        )

                    pbar.update(1)

        self.display_summary(results)
    
    def display_summary(self, results: List[Dict]):
        """显示运行摘要"""
        elapsed_time = time.time() - self.start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        summary_stats = defaultdict(int)
        cobra_extension_stats = defaultdict(int)
        skip_reasons = defaultdict(int)
        
        for result in results:
            status = result.get("status", "unknown")
            summary_stats[status] += 1
            
            if status == "success":
                cobra_status = result.get("cobra_status", "unknown")
                cobra_extension_stats[cobra_status] += 1
                
                skip_reason = result.get("skip_reason")
                if skip_reason:
                    skip_reasons[skip_reason] += 1
        
        summary_lines = [
            "=" * 60,
            "COBRA批量分析流程完成",
            "=" * 60,
            f"运行时间: {int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}",
            f"总任务数: {len(results)}",
            f"成功任务: {summary_stats['success']}",
            f"  - 有延伸结果: {cobra_extension_stats.get('success', 0)}",
            f"  - 无延伸结果: {cobra_extension_stats.get('no_extension', 0)}",
            f"  - 跳过COBRA: {cobra_extension_stats.get('skipped', 0)}",
            f"跳过任务: {summary_stats['skipped']}",
            f"失败任务: {summary_stats['failed']}",
            f"异常任务: {summary_stats['exception']}",
            "=" * 60
        ]
        
        main_logger.info("\n")
        for line in summary_lines:
            if line.startswith("="):
                main_logger.info(line)
            else:
                main_logger.info(f"  {line}")
        
        if skip_reasons:
            main_logger.info("\n跳过COBRA原因统计:")
            for reason, count in skip_reasons.items():
                main_logger.info(f"  {reason}: {count} 个任务")
        
        if summary_stats['failed'] > 0 or summary_stats['exception'] > 0:
            main_logger.info("\n失败/异常任务详情:")
            for result in results:
                if result.get("status") in ["failed", "exception"]:
                    main_logger.info(f"  {result['sample']}_{result['assembler']}: {result.get('error', '未知错误')}")
        
        no_extension_tasks = [r for r in results if r.get("cobra_status") == "no_extension"]
        if no_extension_tasks:
            main_logger.info("\n无COBRA延伸结果的任务:")
            for result in no_extension_tasks[:10]:
                main_logger.info(f"  {result['sample']}_{result['assembler']}")
            if len(no_extension_tasks) > 10:
                main_logger.info(f"  还有 {len(no_extension_tasks) - 10} 个任务...")
        
        skipped_cobra_tasks = [r for r in results if r.get("cobra_status") == "skipped"]
        if skipped_cobra_tasks:
            main_logger.info("\n跳过COBRA运行的任务:")
            for result in skipped_cobra_tasks[:10]:
                skip_reason = result.get("skip_reason", "未知原因")
                main_logger.info(f"  {result['sample']}_{result['assembler']} ({skip_reason})")
            if len(skipped_cobra_tasks) > 10:
                main_logger.info(f"  还有 {len(skipped_cobra_tasks) - 10} 个任务...")
        
        main_logger.info(f"\n分析完成! 详细日志请查看各样本目录下的日志文件。")

# ==================== 主函数 ====================

def main():
    parser = argparse.ArgumentParser(
        description="COBRA病毒序列扩展分析脚本 - 批量并行版本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 自动发现mix模式下的所有样本
  python cobra_pipeline.py --mode mix -a all --jobs 4 --threads 8
  
  # 处理所有模式（virus, other, mix）的样本
  python cobra_pipeline.py --mode all -a megahit,penguin
  
  # 手动指定样本
  python cobra_pipeline.py -s "CRR1126132.unmapped.mix,CRR1126133.unmapped.mix" -a megahit,penguin
  
  # 使用断点运行
  python cobra_pipeline.py --mode all -a all --resume

  # 干跑模式，只显示任务列表
  python cobra_pipeline.py --mode all -a all --dry-run
  
输入文件结构:
  7.Virus.reads/1.virus-reads/     # virus模式reads
  7.Virus.reads/2.other-reads/      # other模式reads  
  7.Virus.reads/3.mix-reads/        # mix模式reads
  
  8.virus_assembly/1.virus-assembly/    # virus模式组装结果
  8.virus_assembly/2.other-assembly/    # other模式组装结果
  8.virus_assembly/3.mix-assembly/      # mix模式组装结果
  
  9.virus_identification/virus/         # virus模式鉴定结果
  9.virus_identification/other/         # other模式鉴定结果
  9.virus_identification/mix/           # mix模式鉴定结果
        """
    )
    
    # 模式参数
    parser.add_argument("--mode", choices=["virus", "other", "mix", "all"], default="mix",
                        help="运行模式：virus(病毒), other(其他), mix(混合), all(所有) (默认: mix)")
    
    # 样本参数
    parser.add_argument("-s", "--samples",
                        help="样本列表，可以是逗号分隔的字符串或包含样本列表的文件。如不指定，将根据--mode自动发现样本")
    
    # 目录参数
    parser.add_argument("-c", "--contigs-dir", default="8.virus_assembly/3.mix-assembly/",
                        help="contigs文件目录路径(默认: 8.virus_assembly/3.mix-assembly/)")
    parser.add_argument("-v", "--virsorter-dir", default="9.virus_identification/mix/",
                        help="病毒鉴定结果目录路径(默认: 9.virus_identification/mix/)")
    parser.add_argument("-r", "--reads-dir", default="7.Virus.reads/3.mix-reads/",
                        help="原始reads目录路径(默认: 7.Virus.reads/3.mix-reads/)")
    parser.add_argument("-o", "--output-dir", default="10.cobra_result/",
                        help="输出目录路径(默认: 10.cobra_result/)")
    
    # 工具选择参数
    parser.add_argument("-a", "--assembly-tools", default="all",
                        help="组装工具列表，逗号分隔或'all'(默认: all)")

    # 并行和性能参数
    parser.add_argument("--jobs", type=int, default=1,
                        help="并行任务数(默认: 1)")
    parser.add_argument("--threads", type=int, default=20,
                        help="单个任务的线程数(默认: 20)")
    
    # 断点运行参数
    parser.add_argument("--resume", action="store_true", default=True,
                        help="断点运行: 跳过已成功完成的任务 (默认开启)")
    parser.add_argument("--no-resume", action="store_false", dest="resume",
                        help="禁用断点续传, 强制重跑所有任务")
    
    # COBRA参数
    parser.add_argument("--mink", type=int, default=21,
                        help="最小kmer值(默认: 21)")
    parser.add_argument("--maxk", type=int, default=141,
                        help="最大kmer值(默认: 141)")
    parser.add_argument("--linkage-mismatch", type=int, default=2,
                        help="链接识别不匹配数(默认: 2)")
    
    # 其他参数
    parser.add_argument("--dry-run", action="store_true",
                        help="只显示将要运行的任务，不实际执行")
    parser.add_argument("--verbose", action="store_true",
                        help="显示详细日志")
    
    # 检查是否提供了参数
    if len(sys.argv) == 1:
        parser.print_help()
        print("\n注意: 脚本需要参数才能运行。请至少指定一个参数，或使用示例中的命令。")
        print("例如: python cobra_pipeline.py --mode mix -a all --jobs 4 --threads 8")
        sys.exit(0)
    
    args = parser.parse_args()
    
    # 验证参数
    required_dirs = [
        ("contigs_dir", "Contigs目录"),
        ("virsorter_dir", "病毒鉴定目录"),
        ("reads_dir", "Reads目录")
    ]
    
    for dir_arg, dir_name in required_dirs:
        dir_path = getattr(args, dir_arg)
        if not Path(dir_path).exists():
            main_logger.error(f"{dir_name}不存在: {dir_path}")
            sys.exit(1)
    
    # 创建输出目录
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # 运行流程
    try:
        pipeline = CobraPipeline(args)
        
        if args.dry_run:
            main_logger.info("干跑模式 - 只显示任务列表")
            main_logger.info(f"模式: {args.mode}")
            main_logger.info(f"样本: {len(pipeline.samples)} 个")
            main_logger.info(f"组装工具: {', '.join(pipeline.assembly_tools)}")
            main_logger.info(f"总任务数: {len(pipeline.samples) * len(pipeline.assembly_tools)}")
            
            tasks_to_show = min(5, len(pipeline.samples) * len(pipeline.assembly_tools))
            main_logger.info(f"前 {tasks_to_show} 个任务:")
            count = 0
            for sample in pipeline.samples:
                for assembler in pipeline.assembly_tools:
                    if count < tasks_to_show:
                        main_logger.info(f"  - {sample} [{assembler}]")
                        count += 1
        else:
            pipeline.run()
            
    except KeyboardInterrupt:
        main_logger.info("\n用户中断执行")
        sys.exit(1)
    except Exception as e:
        main_logger.error(f"流程执行失败: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
