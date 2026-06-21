#!/usr/bin/env python3
# virus_identification_pipeline_v6.py
# 病毒鉴定流程自动化脚本 (版本6，简洁清晰的输出)

import os
import sys
import argparse
import subprocess
import shutil
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import glob

# 默认参数
DEFAULT_INPUT_DIRS = "8.virus_assembly"
DEFAULT_READS_TYPE = "all"
DEFAULT_ASSEMBLY_TOOLS = "all"
DEFAULT_IDENTIFY_TOOLS = "all"
DEFAULT_JOBS = 1
DEFAULT_THREADS = 20
DEFAULT_OUTPUT_DIR = "9.virus_identification"
DEFAULT_DB_DIR = os.path.expanduser("~/database/virus-db")

# 组装工具到文件后缀的映射
ASSEMBLY_TOOL_MAPPING = {
    "megahit": "_megahit.contig.fasta",
    "penguinls": "_penguin.contig.fasta",
    "rnaviralspades": "_rnaviralspades.contig.fasta"
}

# reads类型到子目录的映射
READS_TYPE_MAPPING = {
    "virus": "1.virus-assembly",
    "other": "2.other-assembly",
    "mix": "3.mix-assembly"
}

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="病毒鉴定流程自动化脚本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用模式:
  模式1: 指定输入目录，批量处理
    python %(prog)s --input_dirs 8.virus_assembly --reads mix --assembly_tools penguinls --identify_tools all
    
  模式2: 指定单个组装文件，直接处理
    python %(prog)s --assembly_contigs path/to/assembly.fasta --sample_name SAMPLE --identify_tools all

断点续跑功能:
  - 默认自动跳过已完成的样品
  - 使用--force强制重新运行
  - 使用--clean_failed清理失败的任务

注意:
  1. 需要提前安装以下工具和环境:
     - genomad, diamond, blastn, rdrpcatch, viralm, seqkit
  2. 需要设置正确的数据库路径
        """
    )
    
    # 模式选择参数组
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument("--input_dirs", default=DEFAULT_INPUT_DIRS,
                       help="输入目录，包含1.virus-assembly, 2.other-assembly, 3.mix-assembly子目录")
    mode_group.add_argument("--assembly_contigs",
                       help="直接指定组装好的contigs文件路径")
    
    # 通用参数
    parser.add_argument("--reads", choices=["virus", "other", "mix", "all"],
                       default=DEFAULT_READS_TYPE,
                       help="指定reads类型 (当使用--input_dirs时有效)")
    
    parser.add_argument("--assembly_tools", 
                       choices=["megahit", "penguinls", "rnaviralspades", "all"],
                       default=DEFAULT_ASSEMBLY_TOOLS,
                       help="指定组装工具 (当使用--input_dirs时有效)")
    
    parser.add_argument("--identify_tools",
                       choices=["genomad", "blast", "rdrpcatch", "viralm", "all"],
                       default=DEFAULT_IDENTIFY_TOOLS,
                       help="指定鉴定工具")
    
    parser.add_argument("--jobs", type=int, default=DEFAULT_JOBS,
                       help="并行运行任务数")
    
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS,
                       help="每个任务的线程数")
    
    parser.add_argument("--output", default=DEFAULT_OUTPUT_DIR,
                       help="输出目录 (直接模式时，结果直接输出到此目录)")
    
    parser.add_argument("--db_dir", default=DEFAULT_DB_DIR,
                       help="数据库目录")
    
    parser.add_argument("--sample_name", 
                       help="样品名 (当使用--assembly_contigs时必须指定)")
    
    parser.add_argument("--force", action="store_true",
                       help="强制重新运行所有任务，覆盖已有结果")
    
    parser.add_argument("--clean_failed", action="store_true",
                       help="清理失败的任务并重新运行")
    
    parser.add_argument("--dry_run", action="store_true",
                       help="只显示将要运行的任务，不实际执行")
    
    return parser.parse_args()

def validate_args(args):
    """验证参数"""
    # 验证jobs和threads
    if args.jobs < 1:
        raise ValueError("--jobs 必须是正整数")
    
    if args.threads < 1:
        raise ValueError("--threads 必须是正整数")
    
    # 模式1: 使用input_dirs
    if args.input_dirs:
        if not os.path.isdir(args.input_dirs):
            raise ValueError(f"输入目录不存在: {args.input_dirs}")
    
    # 模式2: 使用assembly_contigs
    if args.assembly_contigs:
        if not os.path.isfile(args.assembly_contigs):
            raise ValueError(f"组装文件不存在: {args.assembly_contigs}")
        
        if not args.sample_name:
            raise ValueError("当使用--assembly_contigs时，必须指定--sample_name")
    
    # 验证数据库目录
    if not os.path.isdir(args.db_dir):
        print(f"警告: 数据库目录 {args.db_dir} 不存在", file=sys.stderr)
        response = input("是否继续? (y/n): ")
        if response.lower() != 'y':
            sys.exit(1)
    
    return True

def get_subdirectories_for_reads(args):
    """根据reads类型获取要处理的子目录列表"""
    subdirs = []
    
    if args.reads == "all":
        for reads_type, subdir_name in READS_TYPE_MAPPING.items():
            subdir_path = os.path.join(args.input_dirs, subdir_name)
            if os.path.isdir(subdir_path):
                subdirs.append((reads_type, subdir_path))
    else:
        if args.reads in READS_TYPE_MAPPING:
            subdir_name = READS_TYPE_MAPPING[args.reads]
            subdir_path = os.path.join(args.input_dirs, subdir_name)
            if os.path.isdir(subdir_path):
                subdirs.append((args.reads, subdir_path))
    
    return subdirs

def get_assembly_files_from_subdir(subdir_path, reads_type, assembly_tools):
    """从子目录中获取组装文件列表"""
    assembly_files = []
    
    # 获取该子目录下的所有样品目录
    sample_dirs = []
    for item in os.listdir(subdir_path):
        item_path = os.path.join(subdir_path, item)
        if os.path.isdir(item_path) and not item.startswith("."):
            # 检查是否是样品目录（包含.unmapped.）
            if ".unmapped." in item:
                sample_dirs.append(item)
    
    # 确定要搜索的组装工具后缀
    tool_suffixes = []
    if assembly_tools == "all":
        tool_suffixes = list(ASSEMBLY_TOOL_MAPPING.values())
    else:
        if assembly_tools in ASSEMBLY_TOOL_MAPPING:
            tool_suffixes.append(ASSEMBLY_TOOL_MAPPING[assembly_tools])
    
    # 搜索每个样品目录中的组装文件
    for sample_dir in sample_dirs:
        sample_dir_path = os.path.join(subdir_path, sample_dir)
        
        for suffix in tool_suffixes:
            # 构建可能的文件名模式
            pattern = os.path.join(sample_dir_path, f"{sample_dir}{suffix}")
            
            # 使用glob查找文件
            for file_path in glob.glob(pattern):
                if os.path.isfile(file_path):
                    assembly_files.append(file_path)
    
    return assembly_files

def get_assembly_files_from_input_dirs(args):
    """从输入目录获取组装文件列表"""
    all_files = []
    
    # 获取要处理的子目录
    subdirs = get_subdirectories_for_reads(args)
    
    for reads_type, subdir_path in subdirs:
        # 获取该子目录下的组装文件
        files = get_assembly_files_from_subdir(subdir_path, reads_type, args.assembly_tools)
        all_files.extend([(f, reads_type) for f in files])
    
    return all_files

def get_assembly_files_from_single_contig(args):
    """从单个组装文件获取任务列表"""
    # 直接模式下，不推断reads类型和组装工具
    reads_type = "direct"
    assembly_tool = "direct"
    
    return [(args.assembly_contigs, reads_type, assembly_tool, args.sample_name)]

def get_all_assembly_files(args):
    """获取所有组装文件列表"""
    if args.assembly_contigs:
        # 模式2: 使用单个组装文件
        return get_assembly_files_from_single_contig(args)
    else:
        # 模式1: 使用输入目录
        files_with_reads = get_assembly_files_from_input_dirs(args)
        # 转换为与模式2相同的格式
        result = []
        for filepath, reads_type in files_with_reads:
            sample, _, assembly_tool = extract_info_from_filename(filepath, reads_type)
            result.append((filepath, reads_type, assembly_tool, sample))
        return result

def extract_info_from_filename(filepath, reads_type):
    """从文件名提取信息"""
    filename = os.path.basename(filepath)
    dirname = os.path.dirname(filepath)
    
    # 从目录名获取样品目录名
    sample_dir_name = os.path.basename(dirname)
    
    # 从样品目录名提取样品名
    if sample_dir_name.endswith(f".unmapped.{reads_type}"):
        sample = sample_dir_name[:-len(f".unmapped.{reads_type}")]
    else:
        # 尝试其他可能的格式
        parts = sample_dir_name.split(".")
        if len(parts) >= 3 and parts[1] == "unmapped":
            sample = parts[0]
        else:
            sample = sample_dir_name.split(".")[0]
    
    # 从文件名提取组装工具
    assembly_tool = None
    for tool, suffix in ASSEMBLY_TOOL_MAPPING.items():
        if suffix in filename:
            assembly_tool = tool
            break
    
    if assembly_tool is None:
        # 尝试从文件名推断
        if "_megahit" in filename:
            assembly_tool = "megahit"
        elif "_penguin" in filename:
            assembly_tool = "penguinls"
        elif "_rnaviralspades" in filename:
            assembly_tool = "rnaviralspades"
        else:
            assembly_tool = "unknown"
    
    return sample, reads_type, assembly_tool

def check_tools_available():
    """检查所需工具是否可用"""
    required_tools = ["seqkit"]
    optional_tools = ["genomad", "diamond", "blastn", "rdrpcatch", "viralm"]
    
    print("检查工具可用性:")
    
    available = True
    for tool in required_tools:
        if shutil.which(tool):
            print(f"  ✓ {tool}")
        else:
            print(f"  ✗ {tool} (必需)")
            print(f"错误: 必需工具 {tool} 未找到", file=sys.stderr)
            available = False
    
    for tool in optional_tools:
        if shutil.which(tool):
            print(f"  ✓ {tool}")
        else:
            print(f"  ✗ {tool} (可选)")
    
    return available

def run_command(cmd, log_file=None, timeout=None):
    """运行命令并记录日志"""
    try:
        if log_file:
            with open(log_file, 'w') as f:
                result = subprocess.run(cmd, shell=True, check=True, 
                                       stdout=f, stderr=subprocess.STDOUT,
                                       timeout=timeout)
        else:
            result = subprocess.run(cmd, shell=True, check=True, 
                                   capture_output=True, text=True,
                                   timeout=timeout)
        return True, ""
    except subprocess.CalledProcessError as e:
        error_msg = f"命令执行失败: {cmd}\n返回码: {e.returncode}"
        if log_file and os.path.exists(log_file):
            with open(log_file, 'r') as f:
                error_msg += f"\n错误输出:\n{f.read()[-1000:]}"  # 只显示最后1000字符
        return False, error_msg
    except subprocess.TimeoutExpired:
        error_msg = f"命令执行超时: {cmd}"
        return False, error_msg
    except Exception as e:
        error_msg = f"命令执行异常: {cmd}\n异常: {str(e)}"
        return False, error_msg

def run_genomad(assembly_contigs, sample, output_dir, threads, db_dir, dry_run=False):
    """运行genomad鉴定"""
    print(f"  [{sample}] 运行 genomad...")
    
    if dry_run:
        print(f"    DRY RUN: 将运行 genomad")
        return os.path.join(output_dir, f"{sample}_virus.genomad.result.id")
    
    # 创建输出目录
    genomad_output = os.path.join(output_dir, "genomad_output")
    os.makedirs(genomad_output, exist_ok=True)
    
    # 检查genomad是否可用
    if not shutil.which("genomad"):
        print(f"  [{sample}] 警告: genomad 未找到，跳过genomad鉴定", file=sys.stderr)
        result_file = os.path.join(output_dir, f"{sample}_virus.genomad.result.id")
        open(result_file, 'w').close()
        return result_file
    
    # 构建命令
    cmd = (f"genomad end-to-end --cleanup --threads {threads} "
           f"'{assembly_contigs}' '{genomad_output}' "
           f"'{db_dir}/genomad_db'")
    
    # 运行命令
    log_file = os.path.join(output_dir, "genomad.log")
    success, error_msg = run_command(cmd, log_file, timeout=3600)
    
    if not success:
        print(f"  [{sample}] genomad 执行失败: {error_msg}", file=sys.stderr)
        result_file = os.path.join(output_dir, f"{sample}_virus.genomad.result.id")
        open(result_file, 'w').close()
        return result_file
    
    # 提取结果
    result_file = os.path.join(output_dir, f"{sample}_virus.genomad.result.id")
    virus_summary = None
    
    # 查找virus_summary.tsv文件
    for root, dirs, files in os.walk(genomad_output):
        for file in files:
            if file.endswith("virus_summary.tsv"):
                virus_summary = os.path.join(root, file)
                break
        if virus_summary:
            break
    
    if virus_summary and os.path.exists(virus_summary):
        with open(virus_summary, 'r') as f_in, open(result_file, 'w') as f_out:
            lines = f_in.readlines()
            for line in lines[1:]:
                columns = line.strip().split('\t')
                if columns:
                    f_out.write(columns[0] + '\n')
    else:
        print(f"  [{sample}] 警告: 未找到genomad病毒摘要文件", file=sys.stderr)
        open(result_file, 'w').close()
    
    print(f"  [{sample}] genomad 完成")
    return result_file

def run_blast(assembly_contigs, sample, output_dir, threads, db_dir, dry_run=False):
    """运行blast鉴定"""
    print(f"  [{sample}] 运行 blast...")
    
    if dry_run:
        print(f"    DRY RUN: 将运行 blast (diamond和blastn)")
        return os.path.join(output_dir, f"{sample}_virus.blast.result.id")
    
    # 创建输出目录
    blast_output = os.path.join(output_dir, "blast_output")
    os.makedirs(blast_output, exist_ok=True)
    
    # 检查所需工具
    diamond_cmd = shutil.which("diamond")
    blastn_cmd = shutil.which("blastn")
    
    if not diamond_cmd or not blastn_cmd:
        print(f"  [{sample}] 警告: diamond 或 blastn 未找到，跳过blast鉴定", file=sys.stderr)
        result_file = os.path.join(output_dir, f"{sample}_virus.blast.result.id")
        open(result_file, 'w').close()
        return result_file
    
    # 运行diamond blastx
    diamond_out = os.path.join(blast_output, f"{sample}_virus.diamond_hits.out")
    diamond_cmd_str = (f"diamond blastx "
                      f"-d '{db_dir}/ncbi-virus_ref/ncbi-virus_ref.pep.dmnd' "
                      f"-q '{assembly_contigs}' "
                      f"-o '{diamond_out}' "
                      f"--outfmt 6 --evalue 1e-5 --threads {threads} "
                      f"--block-size 15 --index-chunks 1")
    
    # 运行blastn
    blastn_out = os.path.join(blast_output, f"{sample}_virus.blastn_hits.out")
    blastn_cmd_str = (f"blastn -query '{assembly_contigs}' "
                     f"-db '{db_dir}/ncbi-virus_ref/ncbi-virus_ref.blast.db' "
                     f"-evalue 1e-5 -num_threads {threads} "
                     f"-out '{blastn_out}' -outfmt 6")
    
    # 运行命令
    log_file = os.path.join(output_dir, "blast.log")
    
    # 运行diamond
    success_diamond, error_msg = run_command(diamond_cmd_str, log_file, timeout=7200)
    if not success_diamond:
        print(f"  [{sample}] diamond blastx 执行失败: {error_msg}", file=sys.stderr)
    
    # 运行blastn
    success_blastn, error_msg = run_command(blastn_cmd_str, log_file, timeout=7200)
    if not success_blastn:
        print(f"  [{sample}] blastn 执行失败: {error_msg}", file=sys.stderr)
    
    # 合并结果
    result_file = os.path.join(output_dir, f"{sample}_virus.blast.result.id")
    
    if os.path.exists(diamond_out) and os.path.exists(blastn_out):
        ids = set()
        for infile in [diamond_out, blastn_out]:
            try:
                with open(infile, 'r') as f:
                    for line in f:
                        if line.strip():
                            parts = line.split('\t')
                            if parts:
                                ids.add(parts[0])
            except Exception as e:
                print(f"  [{sample}] 读取文件失败: {infile}: {e}", file=sys.stderr)
        
        with open(result_file, 'w') as f:
            for seq_id in sorted(ids):
                f.write(f"{seq_id}\n")
    else:
        print(f"  [{sample}] 警告: blast输出文件不存在或生成失败", file=sys.stderr)
        open(result_file, 'w').close()
    
    print(f"  [{sample}] blast 完成")
    return result_file

def run_rdrpcatch(assembly_contigs, sample, output_dir, threads, db_dir, dry_run=False):
    """运行rdrpcatch鉴定"""
    print(f"  [{sample}] 运行 rdrpcatch...")
    
    if dry_run:
        print(f"    DRY RUN: 将运行 rdrpcatch")
        return os.path.join(output_dir, f"{sample}_virus.rdrpcatch.result.id")
    
    # 创建输出目录
    rdrpcatch_output = os.path.join(output_dir, "rdrpcatch_output")
    os.makedirs(rdrpcatch_output, exist_ok=True)
    
    # 计算zvalue（contig数量）
    try:
        zvalue_cmd = f"grep -c '^>' '{assembly_contigs}'"
        result = subprocess.run(zvalue_cmd, shell=True, capture_output=True, text=True)
        zvalue = result.stdout.strip()
        if not zvalue.isdigit():
            zvalue = "1000"
    except:
        zvalue = "1000"
    
    # 构建命令
    cmd = (f"conda run -n rdrpcatch rdrpcatch scan -i '{assembly_contigs}' "
           f"-o '{rdrpcatch_output}' "
           f"-db-dir '{db_dir}/rdrp-db/rdrpcatch_dbs/' "
           f"--db-options all --cpus {threads} --zvalue {zvalue} --overwrite")
    
    # 运行命令
    log_file = os.path.join(output_dir, "rdrpcatch.log")
    success, error_msg = run_command(cmd, log_file, timeout=3600)
    
    if not success:
        print(f"  [{sample}] rdrpcatch 执行失败: {error_msg}", file=sys.stderr)
        result_file = os.path.join(output_dir, f"{sample}_virus.rdrpcatch.result.id")
        open(result_file, 'w').close()
        return result_file
    
    # 提取结果
    result_file = os.path.join(output_dir, f"{sample}_virus.rdrpcatch.result.id")
    annotated_file = None
    
    # 查找annotated.tsv文件
    for root, dirs, files in os.walk(rdrpcatch_output):
        for file in files:
            if file.endswith("output_annotated.tsv"):
                annotated_file = os.path.join(root, file)
                break
        if annotated_file:
            break
    
    if annotated_file and os.path.exists(annotated_file):
        ids = set()
        try:
            with open(annotated_file, 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:
                    if line.strip() and not line.startswith("Contig_name"):
                        columns = line.strip().split('\t')
                        if columns:
                            ids.add(columns[0])
            
            with open(result_file, 'w') as f:
                for seq_id in sorted(ids):
                    f.write(f"{seq_id}\n")
        except Exception as e:
            print(f"  [{sample}] 读取rdrpcatch结果文件失败: {e}", file=sys.stderr)
            open(result_file, 'w').close()
    else:
        print(f"  [{sample}] 警告: 未找到rdrpcatch注释文件", file=sys.stderr)
        open(result_file, 'w').close()
    
    print(f"  [{sample}] rdrpcatch 完成")
    return result_file

def run_viralm(assembly_contigs, sample, output_dir, threads, db_dir, dry_run=False):
    """运行viralm鉴定"""
    print(f"  [{sample}] 运行 viralm...")
    
    if dry_run:
        print(f"    DRY RUN: 将运行 viralm")
        return os.path.join(output_dir, f"{sample}_virus.viralm.result.id")
    
    # 创建输出目录
    viralm_output = os.path.join(output_dir, "viralm_result")
    os.makedirs(viralm_output, exist_ok=True)
    
    
    # 构建命令
    cmd = (f"conda run -n viralm  taskset -c 0-60 python ~/biosoft/virus/ViraLM/viralm.py --database '{db_dir}/viralm_db/' "
           f"-i '{assembly_contigs}' "
           f"-o '{viralm_output}' "
           f"--threads {threads} -f ")
    
    # 运行命令
    log_file = os.path.join(output_dir, "viralm.log")
    success, error_msg = run_command(cmd, log_file, timeout=3600)
    
    if not success:
        print(f"  [{sample}] viralm 执行失败: {error_msg}", file=sys.stderr)
        result_file = os.path.join(output_dir, f"{sample}_virus.viralm.result.id")
        open(result_file, 'w').close()
        return result_file
    
    # 提取结果
    result_file = os.path.join(output_dir, f"{sample}_virus.viralm.result.id")
    virus_fasta = None
    
    # 查找virus fasta文件
    for file in os.listdir(viralm_output):
        if file.startswith("virus_") and file.endswith(".fasta"):
            virus_fasta = os.path.join(viralm_output, file)
            break
    
    if virus_fasta and os.path.exists(virus_fasta):
        ids = []
        try:
            with open(virus_fasta, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        seq_id = line[1:].strip().split()[0]
                        ids.append(seq_id)
            
            with open(result_file, 'w') as f:
                for seq_id in ids:
                    f.write(f"{seq_id}\n")
        except Exception as e:
            print(f"  [{sample}] 读取viralm结果文件失败: {e}", file=sys.stderr)
            open(result_file, 'w').close()
    else:
        print(f"  [{sample}] 警告: 未找到viralm病毒fasta文件", file=sys.stderr)
        open(result_file, 'w').close()
    
    print(f"  [{sample}] viralm 完成")
    return result_file

def extract_viral_sequences(assembly_contigs, result_files, sample, output_dir, dry_run=False):
    """提取病毒序列"""
    print(f"  [{sample}] 提取病毒序列...")
    
    if dry_run:
        print(f"    DRY RUN: 将合并结果并提取序列")
        return
    
    # 合并所有结果ID
    all_ids = set()
    for result_file in result_files:
        if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
            try:
                with open(result_file, 'r') as f:
                    for line in f:
                        if line.strip():
                            seq_id = line.strip().split()[0]
                            all_ids.add(seq_id)
            except Exception as e:
                print(f"  [{sample}] 读取结果文件失败 {result_file}: {e}", file=sys.stderr)
    
    # 保存合并的ID
    all_result_file = os.path.join(output_dir, f"{sample}_virus.all.result.id")
    with open(all_result_file, 'w') as f:
        for seq_id in sorted(all_ids):
            f.write(f"{seq_id}\n")
    
    # 使用seqkit提取序列
    seqkit_cmd = shutil.which("seqkit")
    if seqkit_cmd and all_ids:
        candidate_fasta = os.path.join(output_dir, f"{sample}_virus.all.candidate.fasta")
        
        # 创建临时ID文件
        temp_id_file = os.path.join(output_dir, f"{sample}_virus.temp.ids")
        with open(temp_id_file, 'w') as f:
            for seq_id in all_ids:
                f.write(f"{seq_id}\n")
        
        # 使用seqkit提取
        cmd = f"seqkit grep -f '{temp_id_file}' '{assembly_contigs}' > '{candidate_fasta}' 2>&1"
        
        success, error_msg = run_command(cmd, timeout=600)
        
        if success:
            if os.path.exists(candidate_fasta):
                count_cmd = f"grep -c '^>' '{candidate_fasta}'"
                result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
                count = result.stdout.strip()
                if count.isdigit():
                    print(f"  [{sample}] 提取到 {count} 条候选病毒序列")
                else:
                    print(f"  [{sample}] 提取完成")
            else:
                print(f"  [{sample}] 警告: 候选序列文件未生成", file=sys.stderr)
        else:
            print(f"  [{sample}] 警告: seqkit提取序列失败: {error_msg}", file=sys.stderr)
        
        # 清理临时文件
        if os.path.exists(temp_id_file):
            os.remove(temp_id_file)
    else:
        if not seqkit_cmd:
            print(f"  [{sample}] 警告: seqkit 未找到，跳过序列提取", file=sys.stderr)
        elif not all_ids:
            print(f"  [{sample}] 警告: 没有找到候选病毒ID", file=sys.stderr)

def get_output_dir_for_task(filepath, reads_type, assembly_tool, sample, args):
    """获取任务的输出目录"""
    if args.assembly_contigs:
        # 直接模式：输出目录就是args.output
        output_dir = args.output
    else:
        # 目录模式：创建多层子目录
        output_dir = os.path.join(args.output, reads_type, sample, f"{assembly_tool}_result")
    
    return output_dir

def check_task_completion(output_dir, sample):
    """检查任务是否完成"""
    final_result = os.path.join(output_dir, f"{sample}_virus.all.candidate.fasta")
    return os.path.exists(final_result) and os.path.getsize(final_result) > 0

def cleanup_failed_task(output_dir):
    """清理失败的任务目录"""
    if os.path.exists(output_dir):
        try:
            shutil.rmtree(output_dir)
            return True
        except Exception as e:
            print(f"  清理失败目录时出错: {e}", file=sys.stderr)
            return False
    return False

def process_single_sample(filepath, reads_type, assembly_tool, sample, args):
    """处理单个样品"""
    # 获取输出目录
    output_dir = get_output_dir_for_task(filepath, reads_type, assembly_tool, sample, args)
    
    # 检查任务是否已完成
    task_completed = check_task_completion(output_dir, sample)
    
    # 如果是清理失败模式且任务失败，清理目录
    if args.clean_failed and task_completed == False:
        if cleanup_failed_task(output_dir):
            print(f"清理失败任务: {sample}")
        else:
            print(f"无法清理失败任务: {sample}")
    
    # 如果任务已完成且不是强制模式，跳过
    if task_completed and not args.force:
        return sample, True, "已完成，跳过"
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 显示处理信息
    print(f"\n处理样品: {sample}")
    print(f"  组装文件: {filepath}")
    print(f"  输出目录: {output_dir}")
    
    # 如果是dry run模式，只显示信息
    if args.dry_run:
        print("  DRY RUN: 将在此目录运行鉴定流程")
        return sample, True, "dry run"
    
    # 运行指定的鉴定工具
    result_files = []
    
    # 检查需要运行哪些工具
    tools_to_run = []
    if args.identify_tools == "all":
        tools_to_run = ["genomad", "blast", "rdrpcatch", "viralm"]
    else:
        tools_to_run = [args.identify_tools]
    
    for tool in tools_to_run:
        try:
            if tool == "genomad":
                result_file = run_genomad(filepath, sample, output_dir, args.threads, args.db_dir, args.dry_run)
                result_files.append(result_file)
            elif tool == "blast":
                result_file = run_blast(filepath, sample, output_dir, args.threads, args.db_dir, args.dry_run)
                result_files.append(result_file)
            elif tool == "rdrpcatch":
                result_file = run_rdrpcatch(filepath, sample, output_dir, args.threads, args.db_dir, args.dry_run)
                result_files.append(result_file)
            elif tool == "viralm":
                result_file = run_viralm(filepath, sample, output_dir, args.threads, args.db_dir, args.dry_run)
                result_files.append(result_file)
        except Exception as e:
            print(f"  [{sample}] {tool} 运行异常: {e}", file=sys.stderr)
    
    # 合并结果并提取序列
    if result_files and not args.dry_run:
        extract_viral_sequences(filepath, result_files, sample, output_dir, args.dry_run)
    
    print(f"样品 {sample} 处理完成")
    return sample, True, "成功完成"

def main():
    """主函数"""
    # 解析参数
    args = parse_arguments()
    
    # 验证参数
    try:
        validate_args(args)
    except ValueError as e:
        print(f"参数错误: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 检查工具可用性（如果不是dry run）
    if not args.dry_run and not check_tools_available():
        print("错误: 必需工具缺失", file=sys.stderr)
        sys.exit(1)
    
    # 打印参数
    print("=" * 70)
    print("病毒鉴定流程自动化脚本")
    print("=" * 70)
    
    if args.assembly_contigs:
        print(f"运行模式: 直接模式 (使用单个组装文件)")
        print(f"  组装文件: {args.assembly_contigs}")
        print(f"  样品名: {args.sample_name}")
        print(f"  输出目录: {args.output}")
    else:
        print(f"运行模式: 目录模式 (批量处理)")
        print(f"  输入目录: {args.input_dirs}")
        print(f"  reads类型: {args.reads}")
        print(f"  组装工具: {args.assembly_tools}")
        print(f"  输出目录: {args.output}")
    
    print(f"  鉴定工具: {args.identify_tools}")
    print(f"  并行任务: {args.jobs}")
    print(f"  线程数: {args.threads}")
    print(f"  数据库目录: {args.db_dir}")
    print(f"  强制重新运行: {args.force}")
    print(f"  清理失败任务: {args.clean_failed}")
    print(f"  模拟运行: {args.dry_run}")
    print("=" * 70)
    
    # 获取组装文件列表
    print("\n搜索组装文件...")
    assembly_tasks = get_all_assembly_files(args)
    
    if not assembly_tasks:
        print("错误: 未找到任何组装文件", file=sys.stderr)
        print("请检查:")
        print(f"  1. 输入目录是否正确: {args.input_dirs}")
        print(f"  2. reads类型是否正确: {args.reads}")
        print(f"  3. 组装工具是否正确: {args.assembly_tools}")
        sys.exit(1)
    
    print(f"找到 {len(assembly_tasks)} 个组装文件")
    
    # 如果是dry run，显示将要处理的任务
    if args.dry_run:
        print("\n模拟运行 - 将要处理的任务:")
        for i, (filepath, reads_type, assembly_tool, sample) in enumerate(assembly_tasks, 1):
            output_dir = get_output_dir_for_task(filepath, reads_type, assembly_tool, sample, args)
            task_completed = check_task_completion(output_dir, sample)
            status = "已完成" if task_completed else "未完成"
            print(f"{i:3d}. {sample} [{status}]")
            print(f"     文件: {filepath}")
            print(f"     输出目录: {output_dir}")
        print(f"\n总共 {len(assembly_tasks)} 个任务")
        return
    
    # 处理文件
    print(f"\n开始处理样品 (并行数: {args.jobs})...")
    
    results = []
    skipped_count = 0
    completed_count = 0
    
    if args.jobs > 1 and len(assembly_tasks) > 1:
        # 并行处理
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            # 提交任务
            future_to_task = {}
            for task in assembly_tasks:
                filepath, reads_type, assembly_tool, sample = task
                output_dir = get_output_dir_for_task(filepath, reads_type, assembly_tool, sample, args)
                task_completed = check_task_completion(output_dir, sample)
                
                # 如果任务已完成且不是强制模式，跳过
                if task_completed and not args.force:
                    print(f"跳过已完成样品: {sample}")
                    skipped_count += 1
                    results.append((sample, True, "已跳过（已完成）"))
                    continue
                
                future = executor.submit(process_single_sample, filepath, reads_type, assembly_tool, sample, args)
                future_to_task[future] = sample
            
            # 使用tqdm显示进度
            with tqdm(total=len(future_to_task), desc="处理进度") as pbar:
                for future in as_completed(future_to_task):
                    try:
                        sample, success, message = future.result()
                        results.append((sample, success, message))
                        completed_count += 1
                        pbar.update(1)
                        pbar.set_postfix_str(f"已完成: {sample}")
                    except Exception as e:
                        sample = future_to_task[future]
                        print(f"处理样品 {sample} 时出错: {e}", file=sys.stderr)
                        results.append((sample, False, f"处理异常: {str(e)}"))
                        completed_count += 1
                        pbar.update(1)
    else:
        # 串行处理（使用tqdm显示进度）
        for filepath, reads_type, assembly_tool, sample in tqdm(assembly_tasks, desc="处理进度"):
            output_dir = get_output_dir_for_task(filepath, reads_type, assembly_tool, sample, args)
            task_completed = check_task_completion(output_dir, sample)
            
            # 如果任务已完成且不是强制模式，跳过
            if task_completed and not args.force:
                print(f"\n跳过已完成样品: {sample}")
                skipped_count += 1
                results.append((sample, True, "已跳过（已完成）"))
                continue
            
            try:
                sample_result = process_single_sample(filepath, reads_type, assembly_tool, sample, args)
                results.append(sample_result)
                completed_count += 1
            except Exception as e:
                print(f"\n处理样品 {sample} 时出错: {e}", file=sys.stderr)
                results.append((sample, False, f"处理异常: {str(e)}"))
                completed_count += 1
    
    # 打印摘要
    print("\n" + "=" * 70)
    print("运行摘要:")
    print("=" * 70)
    
    total_tasks = len(assembly_tasks)
    failed_tasks = sum(1 for _, success, _ in results if not success)
    successful_tasks = total_tasks - skipped_count - failed_tasks
    
    print(f"总任务数: {total_tasks}")
    print(f"跳过任务: {skipped_count} (已完成的样品)")
    print(f"成功完成: {successful_tasks}")
    print(f"失败任务: {failed_tasks}")
    
    if failed_tasks > 0:
        print("\n失败任务详情:")
        for sample, success, message in results:
            if not success:
                print(f"  {sample}: {message}")
    
    print("\n" + "=" * 70)
    print("病毒鉴定流程完成!")
    print(f"结果保存在: {args.output}")
    print("=" * 70)

if __name__ == "__main__":
    main()
