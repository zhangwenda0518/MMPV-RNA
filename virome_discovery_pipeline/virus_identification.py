#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# virus_identification_pipeline_v16_complete.py
# 病毒鉴定流程自动化脚本 (版本16.0 - 全工具整合版)
# ==============================================================================
# 更新内容:
# 1. 【工具整合】新增 VirSorter2 病毒鉴定模块 (整合自 viralprediction.py)。
# 2. 【工具整合】新增 ViralVerify HMM 病毒鉴定模块 (整合自 viralprediction.py)。
# 3. 【工具整合】新增 VirHunter 深度学习病毒鉴定模块。
# 4. 保留原版 Genomad 简洁调用，不做额外参数控制。
# 5. 【架构复原】全面恢复原版的子目录分类架构。
# 6. 【绘图升级】引入 `venn` 库，支持 2~6 集合的高级 Venn 图绘制！
# 7. 保留进度条置底、样本节点动态追踪 UI、Blastx 多库对抗抢救机制。
# ==============================================================================

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
import resource
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import glob
import time
import re
import warnings
import threading
from Bio import SeqIO
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    psutil = None
    PSUTIL_AVAILABLE = False

# 屏蔽底层库烦人的版本不兼容警告
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# 尝试导入高级 Venn 图库 (支持 2~6 个集合)
try:
    from venn import venn as draw_venn
    VENN_AVAILABLE = True
except ImportError:
    VENN_AVAILABLE = False

# 尝试导入 UpSet 图库 (超过 6 个集合时使用)
try:
    from upsetplot import UpSet, from_contents
    UPSET_AVAILABLE = True
except ImportError:
    UPSET_AVAILABLE = False

# --- 默认参数设置 ---
DEFAULT_JOBS = 1
DEFAULT_THREADS = 20
DEFAULT_OUTPUT_DIR = "5.virus_identification"
DEFAULT_DB_DIR = os.path.expanduser("~/database/virus-db")

LEN_VIRUS_MIN = 500
LEN_VIROID_MIN = 200
LEN_VIROID_MAX = 1000

# 全局资源消耗记录 (由各 process_sample 写入独立文件，最终合并)
# 不在多进程间共享，每样本独立记录

def log_resource(output_dir, sample, tool, wall_sec, cpu_sec=0, mem_mb=0, status="OK"):
    """记录单次工具运行的资源消耗 (写入样本独立文件, 兼容多进程)"""
    tsv_path = os.path.join(output_dir, f"{sample}_resource.tsv")
    header = not os.path.exists(tsv_path)
    with open(tsv_path, 'a') as f:
        if header:
            f.write("sample\ttool\twall_sec\tcpu_sec\tmem_mb\tstatus\n")
        f.write(f"{sample}\t{tool}\t{round(wall_sec, 1)}\t{round(cpu_sec, 1)}\t{round(mem_mb, 1)}\t{status}\n")


def save_resource_summary(output_dir):
    """扫描所有 *_resource.tsv, 生成跨样本资源消耗汇总 TSV + 图表"""
    all_tsvs = glob.glob(os.path.join(output_dir, "**/*_resource.tsv"), recursive=True) + \
               glob.glob(os.path.join(output_dir, "*_resource.tsv"))
    if not all_tsvs:
        return
    dfs = []
    for tsv in all_tsvs:
        try:
            dfs.append(pd.read_csv(tsv, sep='\t'))
        except Exception:
            pass
    if not dfs:
        return
    df = pd.concat(dfs, ignore_index=True)

    # 保存合并记录
    tsv_path = os.path.join(output_dir, "resource_usage.tsv")
    df.to_csv(tsv_path, sep='\t', index=False)

    # 按工具汇总
    summary = df.groupby('tool').agg(
        total_wall_sec=('wall_sec', 'sum'),
        avg_wall_sec=('wall_sec', 'mean'),
        total_cpu_sec=('cpu_sec', 'sum'),
        avg_cpu_sec=('cpu_sec', 'mean'),
        max_mem_mb=('mem_mb', 'max'),
        avg_mem_mb=('mem_mb', 'mean'),
        avg_cpus=('wall_sec', lambda x: (df.loc[x.index, 'cpu_sec'].sum() / x.sum()) if x.sum() > 0 else 0),
        runs=('wall_sec', 'count'),
        failed=('status', lambda x: (x == 'FAIL').sum())
    ).sort_values('total_wall_sec', ascending=False)

    sum_tsv = os.path.join(output_dir, "resource_summary.tsv")
    summary.to_csv(sum_tsv, sep='\t', float_format='%.1f')
    safe_print(f"\n📊 资源消耗记录: {tsv_path}")
    safe_print(f"📊 资源汇总表:   {sum_tsv}")

    # 绘制资源消耗汇总图 (三面板: Wall Time / CPU Time / Max Memory)
    sum_dir = os.path.join(output_dir, "summary_plots")
    os.makedirs(sum_dir, exist_ok=True)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

    tools = summary.index.tolist()
    colors = sns.color_palette("rocket", len(tools))

    ax1.bar(tools, summary['total_wall_sec'] / 60, color=colors)
    ax1.set_title("Total Wall Time (min)", fontsize=12, fontweight='bold')
    ax1.set_ylabel("Minutes")
    ax1.tick_params(axis='x', rotation=45)

    ax2.bar(tools, summary['total_cpu_sec'] / 60, color=colors)
    ax2.set_title("Total CPU Time (min)", fontsize=12, fontweight='bold')
    ax2.set_ylabel("Minutes")
    ax2.tick_params(axis='x', rotation=45)

    ax3.bar(tools, summary['max_mem_mb'] / 1024, color=colors)
    ax3.set_title("Max Memory (GB)", fontsize=12, fontweight='bold')
    ax3.set_ylabel("GB")
    ax3.tick_params(axis='x', rotation=45)

    plt.suptitle(f"Resource Usage Summary ({len(df['sample'].unique())} samples)", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(sum_dir, "Resource_Usage.png"), dpi=300, bbox_inches='tight')
    plt.close()


def set_core_unlimited():
    """设置core dump文件大小无限制"""
    try:
        resource.setrlimit(resource.RLIMIT_CORE, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except Exception:
        pass

def is_file_valid(filepath, min_size=1):
    """检查文件是否存在且大小正常"""
    return os.path.exists(filepath) and os.path.getsize(filepath) > min_size

def safe_print(msg):
    """多进程安全的打印函数，保护底部 tqdm 进度条"""
    try:
        tqdm.write(msg)
    except Exception:
        print(msg)

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="🦠 病毒与类病毒全自动鉴定流水线 (v16.0 全工具整合版)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    io_group = parser.add_argument_group('>>> 输入与输出控制 (I/O Control)')
    io_group.add_argument("-i", "--input", required=True, help="输入路径: 单个 FASTA 文件或目录")
    io_group.add_argument("-ext", "--extension", default=".fasta", help="当输入为目录时搜索的后缀 (例: refineC_merge.merged.fasta)")
    io_group.add_argument("-s", "--sample", help="单文件模式强制指定样品名")
    io_group.add_argument("-o", "--output", default=DEFAULT_OUTPUT_DIR, help="结果输出的根目录")
    io_group.add_argument("--force", action="store_true", help="强制重新运行，覆盖已有结果")
    io_group.add_argument("--skip_run", action="store_true", help="跳过工具运行，仅执行合并+后置过滤 (用于仅修改过滤参数后重新过滤)")
    io_group.add_argument("--clean_failed", action="store_true", help="自动清理运行失败的任务目录")
    io_group.add_argument("--skip_plots", action="store_true", help="跳过图表生成")
    io_group.add_argument("--skip_uniprot_filter", action="store_true", help="跳过 UniProt 后置过滤")
    io_group.add_argument("--skip_nr_filter", action="store_true", help="跳过 NR 后置过滤")

    db_group = parser.add_argument_group('>>> 数据库路径 (Database Paths)')
    db_group.add_argument("--db_dir", default=DEFAULT_DB_DIR, help="数据库根目录")
    db_group.add_argument("--virus_protein_db", help="Diamond 病毒蛋白数据库 (.dmnd)")
    db_group.add_argument("--nr_db", help="Diamond NR 数据库 (.dmnd)")
    db_group.add_argument("--uniprot_db", help="Diamond UniProt 数据库 (.dmnd)")
    db_group.add_argument("--viroids_db", help="类病毒 Blastn 数据库路径")
    db_group.add_argument("--virus_taxid", help="病毒taxid文件路径")

    method_group = parser.add_argument_group('>>> 鉴定方法配置 (Identification Methods)')
    method_group.add_argument("--identify_tools", choices=["genomad", "blast", "rdrpcatch", "viralm", "virbot", "viroid", "virsorter2", "viralverify", "virhunter", "metabuli", "all"], default="all", help="选择鉴定工具")
    method_group.add_argument("--blast_evalue", default="1e-5", help="E-value 阈值")
    method_group.add_argument("--blast_mode", choices=['strict', 'filter', 'both', 'no-filter'], default='filter', help="Blast过滤模式: filter(抢救,默认), strict(严格), both(同时运行两种)")
    method_group.add_argument("--blast_top_n", type=int, default=5, help="Blast对抗验证审查 Top N")
    method_group.add_argument("--virbot_path", default=str(SCRIPT_DIR.parent / "biosoft/VirBot/VirBot.py"), help="VirBot 脚本路径")
    method_group.add_argument("--viralm_path", default=str(SCRIPT_DIR / "utils/viralm_cpu.py"), help="ViraLM CPU 版本脚本路径")
    method_group.add_argument("--virsorter_group", default="dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae",
                       help="VirSorter2 --include-groups (默认: dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae)")
    method_group.add_argument("--virsorter_db",
                       help="VirSorter2 数据库路径 (默认: {db_dir}/db)")
    method_group.add_argument("--viralverify_hmm",
                       help="ViralVerify HMM 模型路径 (默认: {db_dir}/ViralVerify/nbc_hmms.h3m)")
    method_group.add_argument("--virhunter_path",
                       default=str(SCRIPT_DIR.parent / "biosoft/virhunter/predict_cpu.py"),
                       help="VirHunter predict_cpu.py 路径")
    method_group.add_argument("--virhunter_weights",
                       default=str(SCRIPT_DIR.parent / "biosoft/virhunter/weights/generalistic"),
                       help="VirHunter 权重目录 (默认: generalistic)")
    method_group.add_argument("--metabuli_db",
                       help="Metabuli 数据库路径 (默认: {db_dir}/RVDB-v31/RVDB_viroids.metabuli_db)")

    run_group = parser.add_argument_group('>>> 运行资源与其他 (Resources)')
    run_group.add_argument("-j", "--jobs", type=int, default=DEFAULT_JOBS, help="并行样本数")
    run_group.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="单样本线程数")
    run_group.add_argument("--dry_run", action="store_true", help="仅打印流程，不执行命令")
    run_group.add_argument("--no_core_dump", action="store_true", help="禁用 core dump 设置")

    return parser.parse_args()

def validate_args(args):
    """资源校验"""
    total_cores = args.jobs * args.threads
    sys_cores = multiprocessing.cpu_count()
    if total_cores > sys_cores:
        print(f"⚠ 警告: 请求总核心数 ({total_cores}) 超过可用核心数 ({sys_cores})")
        if not args.dry_run:
            time.sleep(1)

def check_tools_available(args):
    """[1/3] 工具依赖检查"""
    target_tools = ["genomad", "blast", "rdrpcatch", "viralm", "virbot", "viroid",
                    "virsorter2", "viralverify", "virhunter", "metabuli"] \
                   if args.identify_tools == "all" else [args.identify_tools]

    print("\n[1/3] 🛠️  相关环境工具检查")
    all_ok = True

    # 必需工具
    if shutil.which("seqkit"):
        print(f"  ✓ seqkit             : 已安装")
    else:
        print(f"  ✗ seqkit             : 未找到 (流程必需)"); all_ok = False

    # 通过 PATH 检查的 CLI 工具
    path_tools = [
        ("genomad", "genomad"),
        ("diamond", "blast"),
        ("blastn", "viroid"),
        ("taxonkit", "blast/Filter"),
        ("virsorter", "virsorter2"),
        ("viralverify", "viralverify"),
        ("metabuli", "metabuli"),
    ]
    for exe, tool_key in path_tools:
        if tool_key in target_tools or args.identify_tools == "all":
            if shutil.which(exe):
                print(f"  ✓ {exe:<18} : 已安装 ({tool_key})")
            else:
                print(f"  - {exe:<18} : 未找到 ({tool_key}, 可选)")

    # conda env 工具
    conda_envs = [
        ("rdrpcatch", "rdrpcatch"),
        ("virhunter", "virhunter"),
        ("viralm", "viralm"),
    ]
    for env, tool_key in conda_envs:
        if tool_key in target_tools or args.identify_tools == "all":
            import subprocess as sp
            r = sp.run(f"conda env list | grep -w {env}", shell=True, capture_output=True, text=True)
            if r.returncode == 0 and env in r.stdout:
                print(f"  ✓ conda env:{env:<12} : 已安装 ({tool_key})")
            else:
                print(f"  - conda env:{env:<12} : 未找到 ({tool_key}, 可选)")

    # Python 脚本路径检查
    if "virbot" in target_tools or args.identify_tools == "all":
        if os.path.exists(args.virbot_path):
            print(f"  ✓ VirBot 脚本        : 找到 ({args.virbot_path})")
        else:
            print(f"  - VirBot 脚本        : 未找到 ({args.virbot_path})")

    if "virhunter" in target_tools or args.identify_tools == "all":
        if os.path.exists(args.virhunter_path):
            print(f"  ✓ VirHunter 脚本     : 找到 ({args.virhunter_path})")
        else:
            print(f"  - VirHunter 脚本     : 未找到 ({args.virhunter_path})")

    if "viralm" in target_tools or args.identify_tools == "all":
        if os.path.exists(args.viralm_path):
            print(f"  ✓ ViralM 脚本        : 找到 ({args.viralm_path})")
        else:
            print(f"  - ViralM 脚本        : 未找到 ({args.viralm_path})")

    # 绘图库
    if not args.skip_plots:
        if VENN_AVAILABLE:
            print(f"  ✓ venn (绘图库)      : 已安装 (支持2~6集合Venn图)")
        else:
            print(f"  - venn (绘图库)      : 未找到 (pip install venn)")
        if UPSET_AVAILABLE:
            print(f"  ✓ upsetplot (绘图库)  : 已安装 (支持>6集合UpSet图)")
        else:
            print(f"  - upsetplot (绘图库)  : 未找到 (pip install upsetplot)")

    if not all_ok:
        sys.exit("\n❌ 致命错误: 缺失必需依赖，流程终止。")
    if not all_ok: sys.exit("\n❌ 致命错误: 缺失必需的依赖工具，流程终止。")

def check_databases_available(args):
    """[2/3] 数据库依赖检查"""
    print("\n[2/3] 🗄️  相关依赖数据库检查")
    all_fatal = False
    target_tools = ["genomad", "blast", "rdrpcatch", "viralm", "virbot", "viroid", "virsorter2", "viralverify", "virhunter", "metabuli"] if args.identify_tools == "all" else [args.identify_tools]

    vp_db = args.virus_protein_db or os.path.join(args.db_dir, "Diamond_VirusProtein_db", "viral_protein.dmnd")
    if not os.path.exists(vp_db): vp_db = os.path.join(args.db_dir, "ncbi-virus_ref", "ncbi-virus_ref.pep.dmnd")
    if os.path.exists(vp_db): print(f"  ✓ 病毒初筛库 (VP)  : 找到 ({vp_db})")
    else:
        if "blast" in target_tools: print(f"  ✗ 病毒初筛库 (VP)  : 未找到 (Blast必需)"); all_fatal = True
        else: print(f"  - 病毒初筛库 (VP)  : 忽略")

    nr_db = args.nr_db or os.path.join(args.db_dir, "Diamond_nr_db", "nr.dmnd")
    if os.path.exists(nr_db): print(f"  ✓ NR 总库          : 找到 ({nr_db})")
    elif "blast" in target_tools and args.blast_mode != 'no-filter': print(f"  ⚠️ NR 总库          : 未找到 (将跳过 NR 库抢救验证)")

    uniprot_db = args.uniprot_db or os.path.join(args.db_dir, "Diamond_uniprot_db", "uniprot.dmnd")
    if os.path.exists(uniprot_db): print(f"  ✓ UniProt 总库     : 找到 ({uniprot_db})")
    elif "blast" in target_tools and args.blast_mode != 'no-filter': print(f"  ⚠️ UniProt 总库     : 未找到 (将跳过 UniProt 库抢救验证)")

    taxid = args.virus_taxid or os.path.join(args.db_dir, "virus.taxid.txt")
    if os.path.exists(taxid): print(f"  ✓ 病毒 TaxID 列表  : 找到 ({taxid})")
    elif "blast" in target_tools and args.blast_mode != 'no-filter': print(f"  ⚠️ 病毒 TaxID 列表  : 未找到 (将由 taxonkit 生成)")

    viroids_db = args.viroids_db or os.path.join(args.db_dir, "viroids-db/viroids.fasta.blast.db")
    has_viroid = any(os.path.exists(f"{viroids_db}{ext}") for ext in ['.nhr', '.nin', '.nsq'])
    if not has_viroid and os.path.exists(viroids_db): has_viroid = any(os.path.exists(f"{viroids_db.replace('.fasta.blast.db', '')}{ext}") for ext in ['.nhr', '.nin', '.nsq'])
    if has_viroid: print(f"  ✓ 类病毒库 (Viroid): 找到 ({viroids_db})")
    else:
        if "viroid" in target_tools: print(f"  ✗ 类病毒库 (Viroid): 未找到"); all_fatal = True

    genomad_db = os.path.join(args.db_dir, "genomad_db")
    if os.path.exists(genomad_db): print(f"  ✓ Genomad 数据库   : 找到 ({genomad_db})")
    else:
        if "genomad" in target_tools: print(f"  ✗ Genomad 数据库   : 未找到"); all_fatal = True

    # VirSorter2 数据库检查
    if "virsorter2" in target_tools:
        vs2_db = args.virsorter_db or os.path.join(args.db_dir, "db")
        if os.path.isdir(vs2_db) or os.path.exists(vs2_db):
            print(f"  ✓ VirSorter2 数据库: 找到 ({vs2_db})")
        else:
            print(f"  ✗ VirSorter2 数据库: 未找到 ({vs2_db})"); all_fatal = True

    # ViralVerify HMM 数据库检查
    if "viralverify" in target_tools:
        vv_hmm = args.viralverify_hmm or os.path.join(args.db_dir, "ViralVerify", "nbc_hmms.h3m")
        if os.path.exists(vv_hmm):
            print(f"  ✓ ViralVerify HMM   : 找到 ({vv_hmm})")
        else:
            print(f"  ✗ ViralVerify HMM   : 未找到 ({vv_hmm})"); all_fatal = True

    # Metabuli 数据库检查
    if "metabuli" in target_tools:
        mb_db = args.metabuli_db or os.path.join(args.db_dir, "RVDB-v31", "RVDB_viroids.metabuli_db")
        if os.path.isdir(mb_db):
            print(f"  ✓ Metabuli 数据库   : 找到 ({mb_db})")
        else:
            print(f"  ✗ Metabuli 数据库   : 未找到 ({mb_db})"); all_fatal = True

    if all_fatal: sys.exit("\n❌ 致命错误: 缺失必需的依赖数据库。")

def run_command(cmd, log_file=None):
    """执行命令并记录日志"""
    try:
        if log_file:
            with open(log_file, 'w') as f:
                subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)
        else:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
        return True, ""
    except subprocess.CalledProcessError as e:
        msg = f"Error code: {e.returncode}"
        if log_file and os.path.exists(log_file):
            try:
                with open(log_file, 'r') as f: msg += f"\nOutput: {f.read()[-500:]}"
            except Exception: pass
        return False, msg
    except Exception as e:
        return False, str(e)

def extract_sample_name(filepath, args):
    """智能提取干净样本名"""
    if args.sample: return args.sample
    parent_dir = os.path.basename(os.path.dirname(filepath))
    if parent_dir and parent_dir not in ['.', '..']:
        clean_name = parent_dir.replace("_clean", "").replace(".unmapped", "")
        if clean_name.lower() not in ["virus", "merged", "results", "mix"]: return clean_name
    base_name = os.path.basename(filepath)
    if base_name.endswith(args.extension): base_name = base_name[:-len(args.extension)]
    for suf in ["_all_tools_refineC_merge", "_megahit.contig", "_penguin.contig", "_rnaviralspades.contig", "_clean", "."]:
        base_name = base_name.replace(suf, "")
    return base_name.strip('_')

def get_file_list(args):
    """递归扫描获取文件列表"""
    tasks = []
    if os.path.isfile(args.input):
        tasks.append((args.input, extract_sample_name(args.input, args)))
        return tasks
    if not os.path.isdir(args.input): sys.exit(f"❌ 错误: 路径 '{args.input}' 不存在！")
    for root, dirs, files in os.walk(args.input):
        for file in files:
            if file.endswith(args.extension):
                filepath = os.path.join(root, file)
                if is_file_valid(filepath, 100):
                    tasks.append((filepath, extract_sample_name(filepath, args)))
    return tasks

def run_seqkit_filtering(input_fasta, base_dir, sample, threads, dry_run=False):
    """SeqKit 长度拆分"""
    virus_out = os.path.join(base_dir, f"{sample}.virus.candidate.fasta")
    viroids_out = os.path.join(base_dir, f"{sample}.viroids.candidate.fasta")
    if dry_run: return virus_out, viroids_out

    cmd_v = f"seqkit seq -g -m {LEN_VIRUS_MIN} -j {threads} '{input_fasta}' -o '{virus_out}'"
    cmd_vd = f"seqkit seq -g -m {LEN_VIROID_MIN} -M {LEN_VIROID_MAX} -j {threads} '{input_fasta}' -o '{viroids_out}'"
    
    success, msg = run_command(f"{cmd_v} && {cmd_vd}", os.path.join(base_dir, "seqkit_filter.log"))
    if not success:
        safe_print(f"  [{sample}] ⚠️ SeqKit 过滤失败: {msg}")
        if os.path.exists(virus_out): os.remove(virus_out)
        if os.path.exists(viroids_out): os.remove(viroids_out)
        return None, None
    return virus_out, viroids_out

def consolidate_virus_results(input_original, result_files, sample, base_dir):
    """提取最终唯一 FASTA"""
    all_ids_file = os.path.join(base_dir, f"{sample}_virus.all.result.id")
    final_fasta = os.path.join(base_dir, f"{sample}_virus.all.candidate.fasta")
    
    ids = set()
    for rf in result_files:
        if is_file_valid(rf, 1):
            with open(rf) as f:
                for line in f:
                    if line.strip(): ids.add(line.strip().split()[0])
    
    with open(all_ids_file, 'w') as f:
        f.write("\n".join(sorted(ids)) + "\n")
            
    if ids:
        os.system(f"seqkit grep -f '{all_ids_file}' '{input_original}' -w 0 -o '{final_fasta}' > /dev/null 2>&1")

def read_virus_taxid_file(path):
    t = set()
    try:
        with open(path, 'r') as f:
            for line in f:
                if tid := line.strip(): t.add(tid)
    except: pass
    return t

# ==========================================================
# 核心鉴定模块：Blast 对抗抢救与各大工具调用 (恢复子目录架构)
# ==========================================================
def run_blast_pipeline(input_fasta, sample, base_dir, args):
    """Blast 主流程 — 仅 VP 初筛，NR/UniProt 过滤后移至合并后统一执行"""
    result_file = os.path.join(base_dir, f"{sample}_virus.blast.result.id")
    if is_file_valid(result_file, 1) and not args.force: return result_file
    blast_dir = os.path.join(base_dir, "blast_output")
    os.makedirs(blast_dir, exist_ok=True)
    prefix = os.path.join(blast_dir, sample)

    vp_db = args.virus_protein_db or os.path.join(args.db_dir, "Diamond_VirusProtein_db", "viral_protein.dmnd")
    if not os.path.exists(vp_db): vp_db = os.path.join(args.db_dir, "ncbi-virus_ref", "ncbi-virus_ref.pep.dmnd")
    if not os.path.exists(vp_db):
        open(result_file, 'w').close()
        return result_file

    vp_out = f"{prefix}.vp.txt"
    blast_cols = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','slen','stitle','salltitles','qcovhsp','nident','staxids']
    cmd_vp = f"diamond blastx -q '{input_fasta}' --db '{vp_db}' --long-reads -o '{vp_out}' -e {args.blast_evalue} --threads {args.threads} --block-size 15 --index-chunks 1 --max-target-seqs 5 --outfmt 6 " + " ".join(blast_cols)
    run_command(cmd_vp, os.path.join(blast_dir, "blast_vp.log"))

    if not is_file_valid(vp_out, 10):
        open(result_file, 'w').close()
        return result_file

    df_vp = pd.read_csv(vp_out, sep='\t', header=None, names=blast_cols)
    all_vp_ids = set(df_vp['qseqid'].unique())

    with open(result_file, 'w') as f:
        f.write("\n".join(sorted(all_vp_ids)) + "\n")
    return result_file


def run_post_filter(db_name, db_path, input_fasta, sample, output_dir, args, is_uniprot=True):
    """后置过滤：diamond blastx + taxid 多数投票 + 关键词白名单 + 指标抢救"""
    sub_dir = os.path.join(output_dir, f"{db_name.lower()}_filter_output")
    os.makedirs(sub_dir, exist_ok=True)

    # ===== 获取病毒 TaxID 列表 =====
    virus_taxid_path = args.virus_taxid or os.path.join(args.db_dir, "virus.taxid.txt")
    if not os.path.exists(virus_taxid_path):
        run_command("taxonkit list --ids 10239 --indent '' > " + virus_taxid_path,
                    os.path.join(sub_dir, "taxonkit.log"))
    virus_taxids = read_virus_taxid_file(virus_taxid_path)
    if not virus_taxids:
        safe_print(f"  [{sample}] 无病毒 taxid，跳过 {db_name} 后置过滤")
        return

    # ===== diamond blastx 比对 =====
    out_file = os.path.join(sub_dir, "blastx_result.txt")
    blast_cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
                  'qstart','qend','sstart','send','evalue','bitscore',
                  'slen','stitle','salltitles','qcovhsp','nident']
    if not is_uniprot:
        blast_cols.append('staxids')

    cmd = (f"diamond blastx -q '{input_fasta}' --db '{db_path}' --long-reads "
           f"-o '{out_file}' -e {args.blast_evalue} --threads {args.threads} "
           f"--block-size 15 --index-chunks 1 --max-target-seqs {args.blast_top_n} "
           f"--outfmt 6 " + " ".join(blast_cols))
    success, _ = run_command(cmd, os.path.join(sub_dir, "blastx.log"))

    if not (success and os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        safe_print(f"  [{sample}] {db_name} 比对无结果")
        return

    # ===== 解析: e-value 预过滤 + 排序 + Top-N (参考 contigtax --top) =====
    df = pd.read_csv(out_file, sep='\t', header=None, names=blast_cols)
    for col in ['pident', 'length', 'slen', 'bitscore', 'evalue']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    # 预过滤: 丢弃不可信比对 (BASTA _check_hit + contigtax 双层)
    df = df[df['evalue'] <= 0.001] if 'evalue' in df.columns else df
    # 排序: bitscore 降 + evalue 升
    df = df.sort_values(by=['qseqid', 'bitscore', 'evalue'], ascending=[True, False, True])
    # Top-N: 每个 query 保留前 blast_top_n 条 best hits
    df_topn = df.groupby('qseqid').head(args.blast_top_n)

    # ===== 提取 TaxID =====
    if is_uniprot:
        # 修复: UniProt 标准 tag 是 OX=, TaxID= 也兼容
        df_topn['taxid'] = df_topn['stitle'].astype(str).str.extract(r'(?:OX|TaxID)=(\d+)', expand=False)
    else:
        # NR: 检查 staxids 中所有 taxid
        def nr_first_viral(taxids_str):
            if pd.isna(taxids_str): return None
            for t in str(taxids_str).split(';'):
                t = t.strip()
                if t in virus_taxids: return t
            return str(taxids_str).split(';')[0].strip()
        df_topn['taxid'] = df_topn['staxids'].apply(nr_first_viral)

    # ===== 病毒关键词白名单 (抢救前噬菌体被误判为宿主) =====
    viral_keywords = ['phage', 'virus', 'virion', 'capsid', 'tail', 'head',
                      'portal', 'terminase', 'integrase', 'prophage']

    def has_viral_keyword(stitle):
        sl = str(stitle if pd.notna(stitle) else '').lower()
        return any(kw in sl for kw in viral_keywords)

    all_fasta_ids = set(SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta")).keys())
    all_df_qids = set(df['qseqid'])

    # ===== 判定逻辑 (Top-N 多数投票) =====
    def _compute_passed(mode):
        """根据指定模式计算通过集合"""
        passed = set()
        for qid, group in df_topn.groupby('qseqid'):
            group_valid = group.dropna(subset=['taxid'])
            any_viral = any(t in virus_taxids for t in group_valid['taxid']) if len(group_valid) > 0 else False

            if mode == 'strict':
                if any_viral:
                    passed.add(qid)
            else:  # filter
                if any_viral:
                    passed.add(qid)          # ① 已知病毒
                elif len(group_valid) == 0:
                    pass  # 无 taxid → 留给全局 no-hit
                elif group_valid['stitle'].apply(has_viral_keyword).any():
                    passed.add(qid)          # ② 关键词抢救
                else:
                    for _, row in group_valid.iterrows():
                        l = row.get('length', 0) or 0
                        p = row.get('pident', 0) or 0
                        e = row.get('evalue', 1.0) or 1.0
                        if l < 50 or p < 30 or e > 1e-5:  # ③ 比对不可信抢救
                            passed.add(qid); break
        # 全局 no-hit 保留
        if mode == 'filter':
            passed |= (all_fasta_ids - all_df_qids)
        return passed

    modes_to_run = ['strict', 'filter'] if args.blast_mode.lower() == 'both' else [args.blast_mode.lower()]

    for mode in modes_to_run:
        passed_ids = _compute_passed(mode)
        suffix = f"_{mode}" if args.blast_mode.lower() == 'both' else ""
        sub_dir_mode = os.path.join(output_dir, f"{db_name.lower()}_filter_output{suffix}")
        os.makedirs(sub_dir_mode, exist_ok=True)

        filtered_id = os.path.join(sub_dir_mode, f"{sample}_virus.{db_name.lower()}_filtered.id")
        filtered_fa = os.path.join(sub_dir_mode, f"{sample}_virus.{db_name.lower()}_filtered.fasta")
        with open(filtered_id, 'w') as f:
            f.write("\n".join(sorted(passed_ids)) + "\n")
        if passed_ids:
            os.system(f"seqkit grep -f '{filtered_id}' '{input_fasta}' -w 0 -o '{filtered_fa}' > /dev/null 2>&1")
        safe_print(f"  [{sample}] {db_name}({mode}) 后置过滤: {len(passed_ids)}/{len(all_fasta_ids)} 通过 → {os.path.basename(filtered_fa)}")

        # 在 both 模式下, 也把 blastx_result.txt 拷贝到 mode 子目录 (前面统一跑过了)
        if args.blast_mode.lower() == 'both' and os.path.exists(os.path.join(sub_dir, "blastx_result.txt")):
            if sub_dir_mode != sub_dir:
                os.system(f"cp '{os.path.join(sub_dir, 'blastx_result.txt')}' '{sub_dir_mode}/' 2>/dev/null")

def run_genomad(input_fasta, sample, base_dir, args):
    """Genomad 鉴定"""
    result_file = os.path.join(base_dir, f"{sample}_virus.genomad.result.id")
    if is_file_valid(result_file, 1) and not args.force: return result_file
    gm_dir = os.path.join(base_dir, "genomad_output")
    cmd = f"genomad end-to-end --cleanup --threads {args.threads} '{input_fasta}' '{gm_dir}' '{args.db_dir}/genomad_db'"

    success, _ = run_command(cmd, os.path.join(base_dir, "genomad.log"))

    virus_summary = None
    for root, _, files in os.walk(gm_dir):
        for f in files:
            if f.endswith("virus_summary.tsv"):
                virus_summary = os.path.join(root, f)
                break
        if virus_summary: break

    with open(result_file, 'w') as fout:
        if success and virus_summary and os.path.exists(virus_summary):
            try:
                with open(virus_summary) as fin:
                    next(fin)
                    for line in fin:
                        if cols := line.split('\t'): fout.write(cols[0] + '\n')
            except: pass
    return result_file

def run_virsorter2(input_fasta, sample, base_dir, args):
    """VirSorter2 病毒鉴定 (整合自 viralprediction.py)"""
    result_file = os.path.join(base_dir, f"{sample}_virus.virsorter2.result.id")
    if is_file_valid(result_file, 1) and not args.force: return result_file
    vs2_dir = os.path.join(base_dir, "virsorter2_output")
    os.makedirs(vs2_dir, exist_ok=True)

    vs2_db = args.virsorter_db or os.path.join(args.db_dir, "db")
    cmd = (f"virsorter run -w '{vs2_dir}' -i '{input_fasta}' "
           f"--include-groups '{args.virsorter_group}' -j {args.threads} "
           f"all --min-score 0.5 --min-length 300 --keep-original-seq "
           f"-d '{vs2_db}'")

    success, _ = run_command(cmd, os.path.join(base_dir, "virsorter2.log"))

    # 从 final-viral-combined.fa 提取 ID (比 TSV 更可靠)
    viral_fa = os.path.join(vs2_dir, "final-viral-combined.fa")
    with open(result_file, 'w') as fout:
        if success and os.path.exists(viral_fa):
            try:
                with open(viral_fa) as fin:
                    for line in fin:
                        if line.startswith('>'):
                            seq_id = line[1:].strip().split()[0].split('||')[0]
                            fout.write(seq_id + '\n')
            except Exception:
                pass
    return result_file


def run_viralverify(input_fasta, sample, base_dir, args):
    """ViralVerify HMM 病毒鉴定 (整合自 viralprediction.py)"""
    result_file = os.path.join(base_dir, f"{sample}_virus.viralverify.result.id")
    if is_file_valid(result_file, 1) and not args.force: return result_file
    vv_dir = os.path.join(base_dir, "viralverify_output")
    os.makedirs(vv_dir, exist_ok=True)

    hmm_db = args.viralverify_hmm or os.path.join(args.db_dir, "ViralVerify", "nbc_hmms.h3m")
    cmd = (f"viralverify -f '{input_fasta}' -o '{vv_dir}' "
           f"--hmm '{hmm_db}' -t {args.threads}")

    success, _ = run_command(cmd, os.path.join(base_dir, "viralverify.log"))

    # 从 Prediction_results_fasta/*_virus.fasta 提取 ID
    viral_fas = glob.glob(os.path.join(vv_dir, "Prediction_results_fasta", "*_virus.fasta"))
    with open(result_file, 'w') as fout:
        if success and viral_fas:
            for vf in viral_fas:
                try:
                    with open(vf) as fin:
                        for line in fin:
                            if line.startswith('>'):
                                fout.write(line[1:].strip().split()[0] + '\n')
                except Exception:
                    pass
    return result_file


def run_virhunter(input_fasta, sample, base_dir, args):
    """VirHunter 深度学习病毒鉴定"""
    result_file = os.path.join(base_dir, f"{sample}_virus.virhunter.result.id")
    if is_file_valid(result_file, 1) and not args.force: return result_file
    vh_dir = os.path.join(base_dir, "virhunter_output")
    os.makedirs(vh_dir, exist_ok=True)

    cmd = (f"conda run -n virhunter python '{args.virhunter_path}' "
           f"--input '{input_fasta}' --weights '{args.virhunter_weights}' "
           f"--cpu {args.threads} --return_viral --length 500 "
           f"--out_dir '{vh_dir}'")

    success, _ = run_command(cmd, os.path.join(base_dir, "virhunter.log"))

    # 从 *_viral.fasta 提取序列 ID (grep ">" | sed 's/>//')
    with open(result_file, 'w') as fout:
        if success:
            viral_fastas = glob.glob(os.path.join(vh_dir, "*_viral.fasta"))
            if viral_fastas:
                for vf in viral_fastas:
                    try:
                        with open(vf) as fin:
                            for line in fin:
                                if line.startswith('>'):
                                    seq_id = line[1:].strip().split()[0]
                                    fout.write(seq_id.split('||')[0] + '\n')
                    except Exception:
                        pass
    return result_file


def run_metabuli(input_fasta, sample, base_dir, args):
    """Metabuli 病毒鉴定 (classify + extract tax-id 10239)"""
    result_file = os.path.join(base_dir, f"{sample}_virus.metabuli.result.id")
    if is_file_valid(result_file, 1) and not args.force: return result_file
    mb_dir = os.path.join(base_dir, "metabuli_output")
    os.makedirs(mb_dir, exist_ok=True)

    mb_db = args.metabuli_db or os.path.join(args.db_dir, "RVDB-v31", "RVDB_viroids.metabuli_db")
    prefix = sample

    cmd_classify = (f"metabuli classify --seq-mode 3 --threads {args.threads} "
                    f"'{input_fasta}' '{mb_db}' '{mb_dir}' '{prefix}'")
    success, msg = run_command(cmd_classify, os.path.join(base_dir, "metabuli.log"))

    if not success:
        return result_file

    class_tsv = os.path.join(mb_dir, f"{prefix}_classifications.tsv")
    if os.path.exists(class_tsv):
        cmd_extract = (f"metabuli extract --seq-mode 3 "
                       f"'{input_fasta}' '{class_tsv}' '{mb_db}' "
                       f"--tax-id 10239 --outdir '{mb_dir}'")
        run_command(cmd_extract, os.path.join(base_dir, "metabuli_extract.log"))

    viral_fnas = glob.glob(os.path.join(mb_dir, "*_10239.fna"))
    with open(result_file, 'w') as fout:
        for vf in viral_fnas:
            try:
                with open(vf) as fin:
                    for line in fin:
                        if line.startswith('>'):
                            seq_id = line[1:].strip().split()[0]
                            fout.write(seq_id + '\n')
            except Exception:
                pass
    return result_file


def run_viroid_identification(input_fasta, sample, base_dir, args):
    """类病毒 Blastn 鉴定 (blastn + jcvi best-hit 过滤)"""
    id_out = os.path.join(base_dir, f"{sample}_viroids.result.id")
    if is_file_valid(id_out, 1) and not args.force: return id_out
    viroid_dir = os.path.join(base_dir, "viroid_output")
    os.makedirs(viroid_dir, exist_ok=True)

    txt_out = os.path.join(viroid_dir, f"{sample}_viroids.blastn.result.txt")

    db_path = args.viroids_db or os.path.join(args.db_dir, "viroids-db/viroids.fasta.blast.db")
    cmd = f"blastn -query '{input_fasta}' -db '{db_path}' -evalue {args.blast_evalue} -num_threads {args.threads} -out '{txt_out}' -outfmt 6"

    success, msg = run_command(cmd, os.path.join(viroid_dir, "viroid_blastn.log"))
    if success and os.path.exists(txt_out) and os.path.getsize(txt_out) > 0:
        # jcvi best-hit 过滤: 每个 query 只保留最优比对 (按 bitscore 降序)
        run_command(f"python -m jcvi.formats.blast best -n 1 '{txt_out}'",
                    os.path.join(viroid_dir, "viroid_besthit.log"))
        # 从 .best 文件提取 ID
        best_file = txt_out + ".best"
        src = best_file if os.path.exists(best_file) else txt_out
        os.system(f"awk '{{print $1}}' '{src}' | sort -u > '{id_out}' 2>/dev/null")
    else:
        safe_print(f"  [{sample}] ⚠️ 类病毒 Blastn 异常: {msg}")

    return id_out

def run_other_tools(tool_name, input_fasta, sample, base_dir, args):
    """统一调用其他的 AI/特殊鉴定工具 (含子目录创建)"""
    res_file = os.path.join(base_dir, f"{sample}_virus.{tool_name}.result.id")
    if is_file_valid(res_file, 1) and not args.force: return res_file

    if tool_name == "rdrpcatch":
        out_subdir = os.path.join(base_dir, "rdrpcatch_output")
        try: zval = subprocess.getoutput(f"grep -c '^>' '{input_fasta}'")
        except: zval = "1000"
        cmd = f"conda run -n rdrpcatch rdrpcatch scan -i '{input_fasta}' -o '{out_subdir}' -db-dir '{args.db_dir}/rdrp-db/rdrpcatch_dbs/' --db-options all --cpus {args.threads} --zvalue {zval} --overwrite"
        run_command(cmd, os.path.join(base_dir, "rdrpcatch.log"))
        
        annotated = os.path.join(out_subdir, f"{sample}.virus.candidate_rdrpcatch_output_annotated.tsv")
        with open(res_file, 'w') as f:
            if os.path.exists(annotated):
                try:
                    with open(annotated) as fin:
                        next(fin)
                        for line in fin:
                            if line.strip(): f.write(line.split('\t')[0] + '\n')
                except: pass

    elif tool_name == "viralm":
        core_dump_cmd = "ulimit -c unlimited && "
        out_subdir = os.path.join(base_dir, "viralm_output")
        cmd = f"{core_dump_cmd} conda run -n viralm taskset -c 0-60 python '{args.viralm_path}' --database '{args.db_dir}/viralm_db/' -i '{input_fasta}' -o '{out_subdir}' --processes {args.threads} -f --batch_size 128 --chunk_size 500 --len {LEN_VIRUS_MIN}"
        run_command(cmd, os.path.join(base_dir, "viralm.log"))
        
        virus_fa_out = None
        if os.path.exists(out_subdir):
            for f in os.listdir(out_subdir):
                if f.startswith("virus_") and f.endswith(".fasta"):
                    virus_fa_out = os.path.join(out_subdir, f)
                    break
        if virus_fa_out: os.system(f"grep '^>' '{virus_fa_out}' | sed 's/>//g' | awk '{{print $1}}' > '{res_file}'")
        else: open(res_file, 'w').close()

    elif tool_name == "virbot":
        out_subdir = os.path.join(base_dir, "virbot_output")
        cmd = f"python '{args.virbot_path}' --input '{input_fasta}' --output '{out_subdir}' --sen --threads {args.threads}"
        run_command(cmd, os.path.join(base_dir, "virbot.log"))
        
        found = False
        for f in glob.glob(os.path.join(out_subdir, "*.vb.fasta")):
            os.system(f"grep '^>' '{f}' | sed 's/>//g' | awk '{{print $1}}' > '{res_file}'")
            found = True
            break
        if not found: open(res_file, 'w').close()

    return res_file

# ==========================================================
# 绘图模块 (≤6集合用Venn, >6集合用UpSet)
# ==========================================================
def generate_comparison_plots(sample, base_dir, tools_to_run, args):
    """绘制单样本维恩图和条形图 (输出到 comparison_plots 子目录)"""
    tool_results = {}
    tool_files = {
        "blast": f"{sample}_virus.blast.result.id",
        "genomad": f"{sample}_virus.genomad.result.id",
        "rdrpcatch": f"{sample}_virus.rdrpcatch.result.id",
        "viralm": f"{sample}_virus.viralm.result.id",
        "virbot": f"{sample}_virus.virbot.result.id",
        "virsorter2": f"{sample}_virus.virsorter2.result.id",
        "viralverify": f"{sample}_virus.viralverify.result.id",
        "virhunter": f"{sample}_virus.virhunter.result.id",
        "metabuli": f"{sample}_virus.metabuli.result.id",
    }
    for tool_name in tools_to_run:
        if tool_name in tool_files:
            file_path = os.path.join(base_dir, tool_files[tool_name])
        else:
            file_path = os.path.join(base_dir, f"{sample}_virus.{tool_name}.result.id")
        if is_file_valid(file_path, 1):
            with open(file_path, 'r') as f:
                tool_results[tool_name] = set(line.strip() for line in f if line.strip())
        else:
            tool_results[tool_name] = set()
            
    valid_tools = [tool for tool in tool_results if len(tool_results[tool]) > 0]
    if len(valid_tools) < 2: return
    
    plot_dir = os.path.join(base_dir, "comparison_plots")
    os.makedirs(plot_dir, exist_ok=True)

    # 1. 集合比较图: <=6 工具用 Venn, >6 工具用 UpSet
    if not args.skip_plots:
        labels_dict = {tool: tool_results[tool] for tool in valid_tools}
        total_unique = len(set().union(*labels_dict.values()))

        if len(valid_tools) <= 6 and VENN_AVAILABLE:
            plt.figure(figsize=(10, 10))
            draw_venn(labels_dict)
            plt.title(f"Virus Identification Comparison - {sample}", fontsize=16, fontweight='bold', pad=20)
            plt.suptitle(f"Total unique sequences: {total_unique}", fontsize=12, style='italic')
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f"{sample}_venn_diagram.png"), dpi=300, bbox_inches='tight')
            plt.close()
        elif len(valid_tools) > 6 and UPSET_AVAILABLE:
            upset_data = from_contents(labels_dict)
            UpSet(upset_data, subset_size='count', show_counts=True).plot()
            plt.title(f"Virus Identification Comparison - {sample}", fontsize=16, fontweight='bold', pad=20)
            plt.suptitle(f"Total unique sequences: {total_unique}", fontsize=12, style='italic')
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f"{sample}_upset_plot.png"), dpi=300, bbox_inches='tight')
            plt.close()
        elif len(valid_tools) > 6 and not UPSET_AVAILABLE:
            safe_print(f"  [{sample}] 提示: 工具数({len(valid_tools)})>6 且 upsetplot 未安装，跳过集合比较图 (pip install upsetplot)")
        else:
            safe_print(f"  [{sample}] 提示: venn 库未安装，跳过 Venn 图绘制。")

    # 2. 条形图绘制
    if not args.skip_plots:
        plt.figure(figsize=(14, 8))
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        tool_counts = {tool: len(tool_results[tool]) for tool in tool_results}
        sorted_tools = sorted(tool_counts.items(), key=lambda x: x[1], reverse=True)
        tools = [t[0] for t in sorted_tools]
        counts = [t[1] for t in sorted_tools]
        
        ax1.bar(tools, counts, color=sns.color_palette("husl", len(tools)))
        ax1.set_title(f"Sequences Identified by Each Tool - {sample}", fontsize=14)
        ax1.tick_params(axis='x', rotation=45)
        
        all_sequences = set()
        for ids in tool_results.values(): all_sequences.update(ids)
        
        tool_coverage = Counter()
        for seq_id in all_sequences: tool_coverage[sum(1 for tool in tool_results if seq_id in tool_results[tool])] += 1
        
        coverages = list(range(1, len(tool_results) + 1))
        coverage_counts = [tool_coverage.get(cov, 0) for cov in coverages]
        ax2.bar([f"{cov} tool{'s' if cov>1 else ''}" for cov in coverages], coverage_counts)
        ax2.set_title(f"Sequence Overlap Across Tools - {sample}", fontsize=14)
        
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f"{sample}_bar_charts.png"), dpi=300, bbox_inches='tight')
        plt.close()

def generate_final_summary_plots(output_dir):
    """跨样本最终统计图 (存入全局的 summary_plots 目录)"""
    all_results = {}
    for root, _, files in os.walk(output_dir):
        for file in files:
            if file.endswith(".all.result.id"):
                sample = file.replace("_virus.all.result.id", "")
                ids_file = os.path.join(root, file)
                if os.path.getsize(ids_file) > 0:
                    with open(ids_file) as f:
                        all_results[sample] = sum(1 for line in f if line.strip())
    
    if len(all_results) < 2: return
    
    sum_dir = os.path.join(output_dir, "summary_plots")
    os.makedirs(sum_dir, exist_ok=True)
    
    plt.figure(figsize=(14, 8))
    samples, counts = zip(*sorted(all_results.items(), key=lambda x: x[1], reverse=True))
    plt.bar(samples, counts, color=sns.color_palette("viridis", len(samples)))
    plt.title("Total Virus Sequences Identified per Sample", fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(sum_dir, "All_Samples_Comparison.png"), dpi=300, bbox_inches='tight')
    plt.close()

# ==========================================================
# 主控制逻辑 (带 UI 追踪器)
# ==========================================================
def process_sample(filepath, sample, args):
    """单样本全流程调度"""
    output_dir = os.path.join(args.output, sample)
    os.makedirs(output_dir, exist_ok=True)
    
    final_flag = os.path.join(output_dir, f"{sample}.processing.done")
    if os.path.exists(final_flag) and not args.force: return sample, True, "已完成"
    
    if args.clean_failed and os.path.exists(output_dir):
        try: shutil.rmtree(output_dir); os.makedirs(output_dir, exist_ok=True)
        except Exception: pass

    tracker = []
    def log_step(step_name):
        tracker.append(f"✓{step_name}")
        safe_print(f"  🧬 {sample:<15} -> [{' '.join(tracker)}]")
        
    virus_fasta, viroids_fasta = run_seqkit_filtering(filepath, output_dir, sample, args.threads, args.dry_run)
    if args.dry_run: return sample, True, "Dry Run"
    if not virus_fasta: return sample, False, "SeqKit 过滤失败"
        
    log_step("预过滤")
    has_virus, has_viroid = is_file_valid(virus_fasta, 10), is_file_valid(viroids_fasta, 10)
    target_tools = ["genomad", "blast", "rdrpcatch", "viralm", "virbot", "viroid", "virsorter2", "viralverify", "virhunter", "metabuli"] if args.identify_tools == "all" else [args.identify_tools]
    tools_ran, virus_results = [], []
    
    def run_timed(tool_name, func, *func_args):
        """执行工具函数并记录资源消耗 (wall time + CPU time + max RSS)"""
        # CPU 时间基线 (子进程累计)
        try:
            r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
            cpu0 = r0.ru_utime + r0.ru_stime
        except Exception:
            cpu0 = None

        # 内存监控线程 (监测当前进程及其子孙的 RSS)
        max_rss = [0]
        stop_mon = threading.Event()
        def monitor_mem():
            try:
                parent = psutil.Process()
                while not stop_mon.is_set():
                    rss = parent.memory_info().rss
                    for child in parent.children(recursive=True):
                        try:
                            rss += child.memory_info().rss
                        except Exception:
                            pass
                    if rss > max_rss[0]:
                        max_rss[0] = rss
                    time.sleep(0.5)
            except Exception:
                pass

        if PSUTIL_AVAILABLE:
            mon_thread = threading.Thread(target=monitor_mem, daemon=True)
            mon_thread.start()

        t0 = time.time()
        result = func(*func_args)
        wall = time.time() - t0

        if PSUTIL_AVAILABLE:
            stop_mon.set()
            mon_thread.join(timeout=2)

        cpu_sec = 0
        if cpu0 is not None:
            try:
                r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
                cpu_sec = (r1.ru_utime + r1.ru_stime) - cpu0
            except Exception:
                pass

        mem_mb = max_rss[0] / 1024 / 1024 if PSUTIL_AVAILABLE else 0
        if wall >= 1.0:   # 跳过断点续跑的假运行 (耗时<1秒 = 工具被跳过)
            log_resource(output_dir, sample, tool_name, wall, cpu_sec, mem_mb)
        return result

    # --skip_run: 跳过工具运行, 直接读取已有 result.id 文件
    if args.skip_run:
        safe_print(f"  [{sample}] --skip_run 模式: 跳过工具运行, 读取已有结果文件")
        if has_virus:
            for tool in [t for t in target_tools if t != "viroid"]:
                rf = os.path.join(output_dir, f"{sample}_virus.{tool}.result.id")
                vt = "blast" if tool == "blast" else tool
                if tool == "blast":
                    rf = os.path.join(output_dir, f"{sample}_virus.blast.result.id")
                elif tool in ("genomad", "virbot", "rdrpcatch", "viralm", "virsorter2", "viralverify", "virhunter", "metabuli"):
                    rf = os.path.join(output_dir, f"{sample}_virus.{tool}.result.id")
                if is_file_valid(rf, 1):
                    virus_results.append(rf)
                    tools_ran.append(tool)
            safe_print(f"  [{sample}] 读取到 {len(virus_results)} 个工具结果: {', '.join(tools_ran)}")
            log_step("读取已有结果")

    elif has_virus:
        if "genomad" in target_tools:
            virus_results.append(run_timed("genomad", run_genomad, virus_fasta, sample, output_dir, args))
            tools_ran.append("genomad"); log_step("Genomad")
        if "blast" in target_tools:
            virus_results.append(run_timed("blast", run_blast_pipeline, virus_fasta, sample, output_dir, args))
            tools_ran.append("blast"); log_step("Blast")
        if "virbot" in target_tools:
            virus_results.append(run_timed("virbot", run_other_tools, "virbot", virus_fasta, sample, output_dir, args))
            tools_ran.append("virbot"); log_step("Virbot")
        if "rdrpcatch" in target_tools:
            virus_results.append(run_timed("rdrpcatch", run_other_tools, "rdrpcatch", virus_fasta, sample, output_dir, args))
            tools_ran.append("rdrpcatch"); log_step("RdRpCatch")
        if "viralm" in target_tools:
            virus_results.append(run_timed("viralm", run_other_tools, "viralm", virus_fasta, sample, output_dir, args))
            tools_ran.append("viralm"); log_step("ViralM")
        if "virsorter2" in target_tools:
            virus_results.append(run_timed("virsorter2", run_virsorter2, virus_fasta, sample, output_dir, args))
            tools_ran.append("virsorter2"); log_step("VirSorter2")
        if "viralverify" in target_tools:
            virus_results.append(run_timed("viralverify", run_viralverify, virus_fasta, sample, output_dir, args))
            tools_ran.append("viralverify"); log_step("ViralVerify")
        if "virhunter" in target_tools:
            virus_results.append(run_timed("virhunter", run_virhunter, virus_fasta, sample, output_dir, args))
            tools_ran.append("virhunter"); log_step("VirHunter")
        if "metabuli" in target_tools:
            virus_results.append(run_timed("metabuli", run_metabuli, virus_fasta, sample, output_dir, args))
            tools_ran.append("metabuli"); log_step("Metabuli")

        if virus_results:
            t0 = time.time()
            consolidate_virus_results(filepath, virus_results, sample, output_dir)
            log_step("合并提取")

            # 后置过滤 (合并后统一对 all.candidate.fasta 执行, 每个工具独立输出过滤ID)
            merged_fasta = os.path.join(output_dir, f"{sample}_virus.all.candidate.fasta")
            if is_file_valid(merged_fasta, 100):
                filter_jobs = []
                uni_db = args.uniprot_db or os.path.join(args.db_dir, "Diamond_uniprot_db", "uniprot.dmnd")
                nr_db = args.nr_db or os.path.join(args.db_dir, "Diamond_nr_db", "nr.dmnd")
                if os.path.exists(uni_db) and not args.skip_uniprot_filter:
                    filter_jobs.append(("UniProt", uni_db, True))
                if os.path.exists(nr_db) and not args.skip_nr_filter:
                    filter_jobs.append(("NR", nr_db, False))
                for db_name, db_path, is_uni in filter_jobs:
                    t0 = time.time()
                    run_post_filter(db_name, db_path, merged_fasta, sample, output_dir, args, is_uniprot=is_uni)
                    log_step(f"{db_name}过滤")
                    # 每个工具的原始结果对过滤集做交集，输出过滤后的 ID
                    suffixes = ["_strict", "_filter"] if args.blast_mode.lower() == 'both' else [""]
                    for sfx in suffixes:
                        sub_dir = os.path.join(output_dir, f"{db_name.lower()}_filter_output{sfx}")
                        filtered_ids_file = os.path.join(sub_dir, f"{sample}_virus.{db_name.lower()}_filtered.id")
                        if os.path.exists(filtered_ids_file):
                            passed = set(open(filtered_ids_file).read().strip().split('\n'))
                            for tool in tools_ran:
                                raw_id_file = os.path.join(output_dir, f"{sample}_virus.{tool}.result.id")
                                if os.path.exists(raw_id_file):
                                    raw_ids = set(open(raw_id_file).read().strip().split('\n'))
                                    tool_filtered = raw_ids & passed
                                    with open(os.path.join(sub_dir, f"{sample}_virus.{tool}.{db_name.lower()}_filtered.id"), 'w') as f:
                                        f.write("\n".join(sorted(tool_filtered)) + "\n")

    if has_viroid and "viroid" in target_tools:
        t0 = time.time()
        run_viroid_identification(viroids_fasta, sample, output_dir, args)
        log_step("类病毒")

    if not args.skip_plots and has_virus:
        t0 = time.time()
        generate_comparison_plots(sample, output_dir, tools_ran, args)
        log_step("绘图")

    open(final_flag, 'w').close()
    return sample, True, "鉴定完毕"

def main():
    set_core_unlimited()
    args = parse_arguments()
    validate_args(args)

    print("\n" + "="*70)
    print("🦠 病毒与类病毒全自动鉴定流水线 (v16.0 全工具整合版)")
    print("="*70)
    print("[配置参数]")
    print(f"  * 输入路径 : {args.input} (扩展名过滤: {args.extension})")
    print(f"  * 输出路径 : {args.output}")
    print(f"  * Blast模式: {args.blast_mode} (Top {args.blast_top_n} 审查)")
    print(f"  * 目标工具 : {args.identify_tools}")
    print(f"  * 并行资源 : {args.jobs} 并发样本 x {args.threads} 线程/样本")

    check_tools_available(args)
    check_databases_available(args)

    print("\n[3/3] 📂 待处理样本数检查")
    tasks = get_file_list(args)
    if not tasks: sys.exit(f"  ❌ 致命错误: 未在 '{args.input}' 找到任何文件。")
    print(f"  ✓ 成功匹配到 {len(tasks)} 个待处理样本。\n")

    print("-" * 70)
    print("🚀 所有前置检查通过，开始执行鉴定流水线...")
    print("-" * 70)

    if args.jobs > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            future_to_sample = {executor.submit(process_sample, fp, smp, args): smp for fp, smp in tasks}
            with tqdm(total=len(tasks), desc="总进度", colour="cyan") as pbar:
                for future in as_completed(future_to_sample):
                    try:
                        _, _, msg = future.result()
                        pbar.set_postfix_str(msg)
                    except Exception as e: safe_print(f"\n❌ 任务异常: {e}")
                    pbar.update(1)
    else:
        for fp, smp in tqdm(tasks, desc="总进度", colour="cyan"):
            process_sample(fp, smp, args)
    
    if not args.skip_plots:
        generate_final_summary_plots(args.output)

    save_resource_summary(args.output)

    print(f"\n✅ 流水线运行结束，结果已保存在: {args.output}")

if __name__ == "__main__":
    main()
