#!/usr/bin/env python3
"""
RNA-seq/Genome 数据清洗与格式转换极致管道
流程: fastp (质控去接头) -> seqkit (转 FASTA) -> 清理中间 FASTQ -> [可选] clumpify
整合: 自动单双端嗅探、统一日志归档、动态进度条 UI、内存峰值探针、断点续传状态机
"""

import os
import re
import sys
import time
import glob
import shutil
import logging
import argparse
import subprocess
import concurrent.futures
from threading import Lock
from collections import defaultdict

# ==========================================
# 0. 终端颜色 UI 与动态进度条
# ==========================================
class UI:
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    GRAY = '\033[90m'
    BOLD = '\033[1m'
    RESET = '\033[0m'

    print_lock = Lock()
    total_tasks = 0
    completed_tasks = 0
    success_tasks = 0
    fail_tasks = 0

    @staticmethod
    def update_progress(success=True, init=False):
        """带有 init 参数，支持初始化 0% 的进度条显示"""
        with UI.print_lock:
            if not init:
                UI.completed_tasks += 1
                if success: UI.success_tasks += 1
                else: UI.fail_tasks += 1
            
            percent = int((UI.completed_tasks / UI.total_tasks) * 100) if UI.total_tasks > 0 else 0
            bar_len = 30
            filled_len = int(bar_len * UI.completed_tasks // UI.total_tasks) if UI.total_tasks > 0 else 0
            bar = '█' * filled_len + '░' * (bar_len - filled_len)
            
            color = UI.GREEN if UI.fail_tasks == 0 else UI.YELLOW
            # 清除当前行并重新绘制进度条
            sys.stdout.write(f"\r\033[K{color}{UI.BOLD}进度: [{bar}] {percent}% ({UI.completed_tasks}/{UI.total_tasks}) | 成功: {UI.success_tasks} | 失败: {UI.fail_tasks}{UI.RESET}")
            sys.stdout.flush()

# ==========================================
# 1. 状态机与断点管理器
# ==========================================
class CheckpointManager:
    def __init__(self, outdir, force=False):
        os.makedirs(outdir, exist_ok=True)
        self.filepath = os.path.join(outdir, '.clean_checkpoints')
        self.completed = set()
        if force and os.path.exists(self.filepath):
            os.remove(self.filepath)
        if os.path.exists(self.filepath):
            with open(self.filepath, 'r') as f:
                self.completed = set(line.strip() for line in f if line.strip())

    def is_done(self, sample_id):
        return sample_id in self.completed

    def mark_done(self, sample_id):
        with UI.print_lock:
            self.completed.add(sample_id)
            with open(self.filepath, 'a') as f:
                f.write(sample_id + '\n')

# ==========================================
# 2. 日志与基础配置
# ==========================================
def configure_logger(name, debug=False):
    logger = logging.getLogger(name)
    if logger.hasHandlers(): logger.handlers.clear()
    format_style = '%(asctime)s [%(levelname)s] %(message)s'
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(format=format_style, level=log_level)
    return logger

logger = logging.getLogger(__name__)

# ==========================================
# 3. 原始文件智能配对扫描
# ==========================================
class FastqScanner:
    @staticmethod
    def extract_sample_name(filename):
        basename = os.path.basename(filename)
        # 支持去除 _fastq.gz 或 .fastq.gz 等各种复杂变体
        patterns = [
            r'[._]f(ast)?q(\.(gz|bz2))?$',
            r'[._-][Rr]?[12](?=[._-]|$)',
            r'_R[12]$', r'[._-][Ll]\d+', r'[._-][Ss]\d+',
            r'[._-](clean|filtered|trimmed|sorted|dedup)',
            r'^[._-]+|[._-]+$'
        ]
        for p in patterns:
            basename = re.sub(p, '', basename, flags=re.IGNORECASE)
        
        if not basename:
            match = re.match(r'^([^_]+)', os.path.basename(filename))
            return match.group(1) if match else os.path.basename(filename).split('.')[0]
        return basename

    @staticmethod
    def scan_and_pair(input_dir):
        patterns = ["**/*fq*", "**/*fastq*"]
        files = set()
        for p in patterns:
            files.update(glob.glob(os.path.join(input_dir, p), recursive=True))
        files = sorted(list(files))
        
        samples = defaultdict(list)
        for f in files:
            samples[FastqScanner.extract_sample_name(f)].append(f)

        tasks = []
        for sample, f_list in samples.items():
            f_list = sorted(f_list)
            # 自动寻找 R1 和 R2 标识
            r1 = next((f for f in f_list if re.search(r'([._-][Rr]?1)(?=[._-]|$)', f, re.I)), None)
            r2 = next((f for f in f_list if re.search(r'([._-][Rr]?2)(?=[._-]|$)', f, re.I)), None)
            
            if len(f_list) == 2 and not r1 and not r2:
                r1, r2 = f_list[0], f_list[1]
                
            if r1 and r2: tasks.append((sample, r1, r2))
            else:
                for f in f_list: tasks.append((sample, f, None))
        
        return tasks

# ==========================================
# 4. 底层执行引擎 (含内存探针)
# ==========================================
def execute(cmds, log_prefix="", logs_dir="", sample_id="", step_name=""):
    os.makedirs(logs_dir, exist_ok=True)
    
    returncodes = []
    procs = []
    log_file_handles = []

    for i, cmd in enumerate(cmds):
        tool_name = os.path.basename(cmd[0]).split('.')[0]
        log_file_path = os.path.join(logs_dir, f"{sample_id}.{step_name}.{i}_{tool_name}.log")
        err_f = open(log_file_path, 'w', encoding='utf-8')
        log_file_handles.append((tool_name, err_f, log_file_path))

        err_f.write(f"=== 时间: {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        err_f.write(f"=== 命令: {' '.join(cmd)} ===\n\n")
        err_f.flush()

        proc = subprocess.Popen(cmd, stdout=err_f, stderr=err_f)
        procs.append(proc)

    peak_mem_mb = 0
    while True:
        all_done = True
        current_mem = 0
        for proc in procs:
            if proc.poll() is None:
                all_done = False
                try:
                    if sys.platform == 'linux':
                        with open(f"/proc/{proc.pid}/statm", 'r') as f:
                            current_mem += int(f.read().split()[1]) * 4 / 1024 
                except: pass
        if current_mem > peak_mem_mb: peak_mem_mb = current_mem
        if all_done: break
        time.sleep(0.5)

    for proc in procs:
        proc.wait()
        returncodes.append(proc.returncode)

    for i, c in enumerate(returncodes):
        tool_name, err_f, log_file_path = log_file_handles[i]
        err_f.close()
        if c != 0:
            with UI.print_lock: 
                sys.stdout.write("\033[2K\r") 
                logger.error(f"{log_prefix} ❌ {tool_name} 失败 (Exit: {c}) 日志: {log_file_path}")
                try:
                    with open(log_file_path, 'r') as f:
                        lines = [line.strip() for line in f.readlines() if line.strip()]
                        if lines: logger.error(f"{log_prefix} 🔍 [探针] -> {' | '.join(lines[-3:])}")
                except: pass

    return returncodes, peak_mem_mb

# ==========================================
# 5. 工具命令构造器
# ==========================================
class Tools:
    @staticmethod
    def fastp(sample, r1, r2, outdir, logsdir, threads, dedup, no_compress):
        # FASTQ 暂存目录 (转换 FASTA 成功后会被清理)
        step_dir = os.path.join(outdir, "1.fastp_tmp")
        os.makedirs(step_dir, exist_ok=True)
        ext = ".fq" if no_compress else ".fq.gz"
        out_r1 = os.path.join(step_dir, f"{sample}_1{ext}")
        
        # 显式添加 -g 和 --poly_g_min_len 10，去除技术假象 Poly-G
        cmd = [
            "fastp", "--thread", str(threads),
            "--qualified_quality_phred", "20",
            "--length_required", "50",
            "-g", "--poly_g_min_len", "10",
            "--html", os.path.join(logsdir, f"{sample}_fastp_report.html"),
            "--json", os.path.join(logsdir, f"{sample}_fastp_report.json")
        ]
        if dedup: cmd.append("--dedup")

        if r2:
            out_r2 = os.path.join(step_dir, f"{sample}_2{ext}")
            cmd += ["-i", r1, "-I", r2, "-o", out_r1, "-O", out_r2, "--detect_adapter_for_pe", "--correction"]
            return cmd, out_r1, out_r2
        else:
            cmd += ["-i", r1, "-o", out_r1]
            return cmd, out_r1, None

    @staticmethod
    def seqkit(sample, in_r1, in_r2, outdir, threads, no_compress):
        step_dir = os.path.join(outdir, "2.fasta")
        os.makedirs(step_dir, exist_ok=True)
        ext = ".fa" if no_compress else ".fa.gz"
        out_r1 = os.path.join(step_dir, f"{sample}_1{ext}")
        
        # -w 0 保证 FASTA 序列不换行输出
        cmd_r1 = ["seqkit", "fq2fa", "-j", str(threads), "-w", "0", in_r1, "-o", out_r1]
        
        if in_r2:
            out_r2 = os.path.join(step_dir, f"{sample}_2{ext}")
            cmd_r2 = ["seqkit", "fq2fa", "-j", str(threads), "-w", "0", in_r2, "-o", out_r2]
            return [cmd_r1, cmd_r2], out_r1, out_r2
        else:
            return [cmd_r1], out_r1, None

    @staticmethod
    def clumpify(sample, in_r1, in_r2, outdir, memory, no_compress):
        step_dir = os.path.join(outdir, "3.clumpify")
        os.makedirs(step_dir, exist_ok=True)
        ext = ".fa" if no_compress else ".fa.gz"
        out_r1 = os.path.join(step_dir, f"{sample}_1{ext}")
        
        cmd = [
            "clumpify.sh", f"in1={in_r1}", f"out1={out_r1}",
            "reorder", "dedupe", "subs=0", f"-Xmx{memory}"
        ]
        
        if in_r2:
            out_r2 = os.path.join(step_dir, f"{sample}_2{ext}")
            cmd += [f"in2={in_r2}", f"out2={out_r2}"]
            return cmd, out_r1, out_r2
        else:
            return cmd, out_r1, None

# ==========================================
# 6. 核心清洗管道流
# ==========================================
class CleaningPipeline:
    def run_single(self, task_id, sample, r1, r2, args):
        start_time = time.time()
        max_mem_used = 0.0
        prefix = f"[{UI.CYAN}任务 {task_id:02d}{UI.RESET} | {UI.BOLD}{sample}{UI.RESET}]"
        
        # 创建统一日志目录
        logs_dir = os.path.join(args.output, "logs")
        os.makedirs(logs_dir, exist_ok=True)

        # [Step 1] fastp 质控
        with UI.print_lock:
            sys.stdout.write("\033[2K\r") 
            logger.info(f"{prefix} -> [1/3] fastp 质控过滤...")
            
        fastp_cmd, fastp_r1, fastp_r2 = Tools.fastp(sample, r1, r2, args.output, logs_dir, args.fastp_threads, args.dedup, args.no_compress)
        ret, mem = execute([fastp_cmd], log_prefix=prefix, logs_dir=logs_dir, sample_id=sample, step_name="01_fastp")
        max_mem_used = max(max_mem_used, mem)
        if any(c != 0 for c in ret): return 1

        # [Step 2] seqkit 格式转换 (FASTQ -> FASTA)
        with UI.print_lock:
            sys.stdout.write("\033[2K\r") 
            logger.info(f"{prefix} -> [2/3] seqkit 格式转换 (FASTQ -> FASTA)...")
            
        seqkit_cmds, fasta_r1, fasta_r2 = Tools.seqkit(sample, fastp_r1, fastp_r2, args.output, args.fastp_threads, args.no_compress)
        ret, mem = execute(seqkit_cmds, log_prefix=prefix, logs_dir=logs_dir, sample_id=sample, step_name="02_seqkit")
        max_mem_used = max(max_mem_used, mem)
        if any(c != 0 for c in ret): return 1

        # [清理机制 A] 转换成功后，立即安全删除 fastp 生成的临时 FASTQ 文件
        try:
            if os.path.exists(fastp_r1): os.remove(fastp_r1)
            if fastp_r2 and os.path.exists(fastp_r2): os.remove(fastp_r2)
            with UI.print_lock:
                sys.stdout.write("\033[2K\r")
                logger.info(f"{prefix} 🗑️ 已清理中间 FASTQ 文件")
        except Exception as e:
            with UI.print_lock:
                sys.stdout.write("\033[2K\r")
                logger.warning(f"{prefix} ⚠️ 清理临时文件失败: {e}")

        # [Step 3] clumpify 光学去重 (可选)
        if not args.skip_clumpify:
            with UI.print_lock:
                sys.stdout.write("\033[2K\r") 
                logger.info(f"{prefix} -> [3/3] clumpify 光学去重...")
                
            clump_cmd, _, _ = Tools.clumpify(sample, fasta_r1, fasta_r2, args.output, args.clumpify_memory, args.no_compress)
            ret, mem = execute([clump_cmd], log_prefix=prefix, logs_dir=logs_dir, sample_id=sample, step_name="03_clumpify")
            max_mem_used = max(max_mem_used, mem)
            if any(c != 0 for c in ret): return 1

        # [清理机制 B] 彻底成功后，删除最初的原始输入文件 (仅在启用 --remove-raw 时)
        if getattr(args, 'remove_raw', False):
            try:
                if r1 and os.path.exists(r1): 
                    os.remove(r1)
                    with UI.print_lock:
                        sys.stdout.write("\033[2K\r")
                        logger.info(f"{prefix} 🗑️ 已删除原始数据: {os.path.basename(r1)}")
                if r2 and os.path.exists(r2): 
                    os.remove(r2)
                    with UI.print_lock:
                        sys.stdout.write("\033[2K\r")
                        logger.info(f"{prefix} 🗑️ 已删除原始数据: {os.path.basename(r2)}")
            except Exception as e:
                with UI.print_lock:
                    sys.stdout.write("\033[2K\r")
                    logger.warning(f"{prefix} ⚠️ 删除原始数据失败: {e}")

        elapsed = time.time() - start_time
        with UI.print_lock:
            sys.stdout.write("\033[2K\r")
            logger.info(f"{prefix} ✅ 流程完成 (耗时: {elapsed:.2f} s | 峰值内存: {max_mem_used:.2f} MB)")

        return 0

# ==========================================
# 7. 主程序入口
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Genome 数据清洗与格式转换管道 (Fastp + Seqkit + Clumpify)")
    
    parser.add_argument("-i", "--input", required=True, help="输入目录 (自动扫描 fq/fastq 及其 gz 包)")
    parser.add_argument("-o", "--output", default="./clean_out", help="清洗结果输出目录")
    
    parser.add_argument("-j", "--jobs", type=int, default=2, help="并发处理的样本数 (默认 2)")
    parser.add_argument("-t", "--fastp-threads", type=int, default=4, help="每个样本使用的工作线程数 (默认 4)")
    
    parser.add_argument("--dedup", action="store_true", help="启用 fastp 自带去重")
    parser.add_argument("--skip-clumpify", action="store_true", help="跳过 clumpify 去重步骤")
    parser.add_argument("--clumpify-memory", default="10g", help="clumpify 分配的 Java 堆内存 (例: 10g)")
    
    parser.add_argument("--no-compress", action="store_true", help="最终结果不使用 gzip 压缩")
    parser.add_argument("--force", action="store_true", help="无视断点记录，强制重新跑所有样本")
    parser.add_argument("--dry-run", action="store_true", help="仅显示识别到的样本，不执行计算")
    parser.add_argument("-d", "--debug", action="store_true", help="显示底层执行命令 Debug 信息")
    
    parser.add_argument("--remove-raw", action="store_true", help="【高危】任务成功后自动删除最初始的测序输入文件")

    args = parser.parse_args()
    global logger
    logger = configure_logger(__name__, debug=args.debug)

    if not os.path.exists(args.input):
        logger.error(f"输入目录不存在: {args.input}")
        sys.exit(1)

    tasks = FastqScanner.scan_and_pair(args.input)
    if not tasks:
        logger.error("在输入目录中未找到任何 fastq/fq 文件！")
        sys.exit(1)

    # 自动统计单双端数目
    pe_count = sum(1 for t in tasks if t[2])
    se_count = sum(1 for t in tasks if not t[2])

    print("\n")
    logger.info("============== 数据清洗与转换管道 ==============")
    logger.info(f"解析任务: 发现 {pe_count} 个双端 (PE) 样本，{se_count} 个单端 (SE) 样本。")
    logger.info(f"执行逻辑: 1.fastp 质控 -> 2.seqkit 转 FASTA -> {'[跳过]' if args.skip_clumpify else '3.clumpify 去重'}")
    logger.info(f"并发配置: {args.jobs} 并发样本 × {args.fastp_threads} 线程/样本 = {args.jobs * args.fastp_threads} 逻辑核心流")
    if args.remove_raw:
        logger.warning("【注意】已开启 --remove-raw，处理成功的原始数据将被自动永久删除！")
    logger.info("================================================\n")

    if args.dry_run:
        logger.info("[Dry-Run] 扫描到的任务清单:")
        for sample, r1, r2 in tasks:
            print(f"  - [{sample}]  {'PE 双端' if r2 else 'SE 单端'} \n    R1: {r1}" + (f"\n    R2: {r2}" if r2 else ""))
        sys.exit(0)

    # 检查环境变量依赖
    missing = [tool for tool in ['fastp', 'seqkit', 'clumpify.sh'] if shutil.which(tool) is None]
    if 'clumpify.sh' in missing and args.skip_clumpify: missing.remove('clumpify.sh')
    if missing:
        logger.error(f"环境缺失！找不到以下命令: {', '.join(missing)}")
        sys.exit(1)

    cm = CheckpointManager(args.output, force=args.force)
    pipeline = CleaningPipeline()
    UI.total_tasks = len(tasks)
    start_time_global = time.time()

    # 初始化进度条
    UI.update_progress(success=True, init=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {}
        for idx, (sample, r1, r2) in enumerate(tasks, 1):
            if cm.is_done(sample):
                with UI.print_lock:
                    sys.stdout.write("\033[2K\r")
                    logger.info(f"{UI.GRAY}⏭ [跳过] 任务 {idx:02d} | {sample} (历史已完成){UI.RESET}")
                UI.update_progress(success=True)
                continue

            future = executor.submit(pipeline.run_single, idx, sample, r1, r2, args)
            futures[future] = sample

        for future in concurrent.futures.as_completed(futures):
            sample = futures[future]
            try:
                if future.result() == 0:
                    cm.mark_done(sample)
                    UI.update_progress(success=True)
                else:
                    UI.update_progress(success=False)
            except Exception as e:
                with UI.print_lock:
                    sys.stdout.write("\033[2K\r")
                    logger.error(f"[{sample}] 发生异常: {e}")
                UI.update_progress(success=False)

    end_time_global = time.time()
    print("\n")
    logger.info("============== 管道执行完毕 ==============")
    logger.info(f"总耗时: {(end_time_global - start_time_global)/60:.2f} 分钟")
    logger.info(f"最终输出: 包含清洗后结果与统一日志的 {os.path.abspath(args.output)} 目录")
    
    sys.exit(1 if UI.fail_tasks > 0 else 0)

if __name__ == '__main__':
    main()
