#!/usr/bin/env python3
"""
独立的数据库与索引批量构建工具 (Kraken2, Bowtie2, HISAT2, Minimap2)
【极致优雅版】利用软链接实现 0 秒挂载 Taxonomy，并在 --clean 前自动安全脱穿隔离源文件。
"""

import os
import re
import sys
import time
import shutil
import logging
import argparse
import subprocess

# ==========================================
# 0. 终端颜色 UI 配置
# ==========================================
class UI:
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    PURPLE = '\033[95m'
    GRAY = '\033[90m'
    BOLD = '\033[1m'
    RESET = '\033[0m'

    @staticmethod
    def print_header(text):
        print(f"\n{UI.PURPLE}{UI.BOLD}{'='*65}{UI.RESET}")
        print(f"{UI.PURPLE}{UI.BOLD} {text} {UI.RESET}")
        print(f"{UI.PURPLE}{UI.BOLD}{'='*65}{UI.RESET}")
        
    @staticmethod
    def print_step(current, total, text):
        print(f"\n{UI.CYAN}{UI.BOLD}▶ [{current}/{total}] {text}{UI.RESET}")

    @staticmethod
    def print_skipped(current, total, text):
        print(f"\n{UI.GRAY}{UI.BOLD}⏭ [{current}/{total}] 跳过: {text} (检测到断点){UI.RESET}")

    @staticmethod
    def print_cmd(cmd_list):
        cmd_str = " ".join(cmd_list)
        print(f"{UI.YELLOW}{UI.BOLD}[执行命令] {UI.RESET}{UI.YELLOW}{cmd_str}{UI.RESET}\n")

# ==========================================
# 1. 状态机与断点管理器
# ==========================================
class CheckpointManager:
    def __init__(self, db_dir, force=False):
        os.makedirs(db_dir, exist_ok=True)
        self.filepath = os.path.join(db_dir, '.checkpoints')
        self.completed = set()
        
        if force and os.path.exists(self.filepath):
            os.remove(self.filepath)
            
        if os.path.exists(self.filepath):
            with open(self.filepath, 'r') as f:
                self.completed = set(line.strip() for line in f if line.strip())

    def is_done(self, step_id):
        return step_id in self.completed

    def mark_done(self, step_id):
        self.completed.add(step_id)
        with open(self.filepath, 'a') as f:
            f.write(step_id + '\n')

# ==========================================
# 2. 日志与基础配置
# ==========================================
def configure_logger(name, debug=False, quiet=False, filename=None):
    logger = logging.getLogger(name)
    if logger.hasHandlers(): logger.handlers.clear()
    log_level = logging.DEBUG if debug else (logging.WARNING if quiet else logging.INFO)
    
    if filename:
        file_handler = logging.FileHandler(filename, mode='a')
        file_handler.setLevel(log_level)
        file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        logger.addHandler(file_handler)
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.WARNING)
    console_handler.setFormatter(logging.Formatter(f'{UI.RED}[%(levelname)s] %(message)s{UI.RESET}'))
    logger.addHandler(console_handler)
    return logger

logger = logging.getLogger(__name__)

# ==========================================
# 3. 安全与环境检查
# ==========================================
def check_dependencies(tools):
    missing = []
    for tool in tools:
        executable = tool
        if tool == 'bowtie2': executable = 'bowtie2-build'
        elif tool == 'hisat2': executable = 'hisat2-build'
        elif tool == 'kraken2': executable = 'kraken2-build'
        if shutil.which(executable) is None: missing.append(executable)
    if missing:
        logger.error(f"环境缺失！找不到以下命令: {', '.join(missing)}")
        sys.exit(1)

def validate_args(args):
    for arg in args:
        if type(arg) == str:
            chars = re.findall(r'[^a-zA-Z0-9\. !#=/_-]', str(arg))
            if not re.match(r'^[a-zA-Z0-9\. !#=/_-]+$|^(?![\s\S])', str(arg)): return False, arg, chars
    return True, '', []

# ==========================================
# 4. 核心执行器
# ==========================================
def execute_realtime(cmd, tool_name="SYS"):
    UI.print_cmd(cmd)
    logger.info(f"执行命令: {' '.join(cmd)}")
    start_time = time.time()
    
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in iter(proc.stdout.readline, ''):
            line = line.strip()
            if line:
                sys.stdout.write(f"\r\033[K[{tool_name}] {line[:120]}") 
                sys.stdout.flush()
                logger.info(f"[{tool_name}] {line}")
        proc.wait()
        sys.stdout.write("\n")
        logger.info(f"命令完成，耗时: {time.time() - start_time:.2f}s, 退出码: {proc.returncode}")
        return proc.returncode
    except KeyboardInterrupt:
        logger.error("\n[SYS] 任务被用户强行中断 (Ctrl+C)！")
        proc.terminate()
        proc.wait()
        return 1

# ==========================================
# 5. 索引构建器核心逻辑 (动态任务队列)
# ==========================================
class IndexBuilders:
    @staticmethod
    def kraken2(fasta, db_dir, threads=1, taxonomy_path=None, k2_libs=None, taxid=None, add_library=None, force=False):
        cm = CheckpointManager(db_dir, force=force)
        actions = [] 

        # ------------------------------------
        # 步骤 1: 软链接挂载 Taxonomy
        # ------------------------------------
        if taxonomy_path:
            abs_tax_path = os.path.abspath(taxonomy_path)
            def link_tax():
                tax_dest = os.path.join(db_dir, 'taxonomy')
                try:
                    if os.path.exists(tax_dest) or os.path.islink(tax_dest):
                        if os.path.islink(tax_dest) or os.path.isfile(tax_dest): os.remove(tax_dest)
                        else: shutil.rmtree(tax_dest)
                    os.symlink(abs_tax_path, tax_dest)
                    print(f"{UI.GREEN}  └── 软链接创建成功！0秒挂载。{UI.RESET}")
                    return 0
                except Exception as e:
                    logger.error(f"创建 Taxonomy 软链接失败: {e}")
                    return 1
            actions.append((f"建立软链接挂载 Taxonomy -> {abs_tax_path}", link_tax, "k2_taxo"))
        else:
            cmd = ['kraken2-build', '--threads', str(threads), '--download-taxonomy', '--db', db_dir]
            actions.append(("联网下载 Taxonomy 分类树", cmd, "k2_taxo"))

        # ------------------------------------
        # 步骤 2: 非模式物种 FASTA 标签注入
        # ------------------------------------
        final_fasta_for_k2 = fasta
        if taxid:
            modified_fasta = os.path.join(db_dir, "host_with_taxid.fasta")
            final_fasta_for_k2 = modified_fasta
            def inject_taxid():
                try:
                    seq_count = 0
                    taxid_tag = f"|kraken:taxid|{taxid}"
                    with open(fasta, 'r') as in_f, open(modified_fasta, 'w') as out_f:
                        for line in in_f:
                            if line.startswith('>'):
                                seq_count += 1
                                parts = line.strip().split(maxsplit=1)
                                new_id = f"{parts[0]}{taxid_tag}"
                                out_f.write(f"{new_id} {parts[1]}\n" if len(parts) == 2 else f"{new_id}\n")
                            else:
                                out_f.write(line)
                    print(f"{UI.GREEN}  └── 成功注入 TaxID: {taxid}，修改 {seq_count} 条序列！{UI.RESET}")
                    return 0
                except Exception as e:
                    logger.error(f"注入 TaxID 标签失败: {e}")
                    return 1
            actions.append((f"为主宿主注入 TaxID: {taxid}", inject_taxid, "k2_taxid_inject"))

        # ------------------------------------
        # 步骤 3: 下载背景标准库
        # ------------------------------------
        if k2_libs and k2_libs[0].lower() != 'none':
            for lib in k2_libs:
                cmd = ['k2', 'download-library', '--library', lib, '--db', db_dir]
                actions.append((f"下载官方标准库: {lib}", cmd, f"k2_lib_{lib}"))

        # ------------------------------------
        # 步骤 4: 补充自定义 FASTA 库
        # ------------------------------------
        if add_library:
            for extra_fa in add_library:
                abs_extra_fa = os.path.abspath(extra_fa)
                basename = os.path.basename(abs_extra_fa)
                cmd_extra = ['kraken2-build', '--threads', str(threads), '--add-to-library', abs_extra_fa, '--db', db_dir]
                step_id = f"k2_add_extra_{basename}"
                actions.append((f"补充自定义 FASTA 库: {basename}", cmd_extra, step_id))

        # ------------------------------------
        # 步骤 5: 混入主宿主
        # ------------------------------------
        cmd_add_host = ['kraken2-build', '--threads', str(threads), '--add-to-library', final_fasta_for_k2, '--db', db_dir]
        actions.append((f"混入主宿主参考基因组", cmd_add_host, "k2_add_host"))

        # ------------------------------------
        # 步骤 6: 编译数据库
        # ------------------------------------
        cmd_build = ['kraken2-build', '--threads', str(threads), '--build', '--db', db_dir]
        actions.append(("编译最终的 Hash 数据库 (极其耗时)", cmd_build, "k2_build"))

        # ------------------------------------
        # 步骤 7: 【神级操作】脱穿隔离（解绑软链接）
        # ------------------------------------
        if taxonomy_path:
            def unlink_tax():
                tax_dest = os.path.join(db_dir, 'taxonomy')
                try:
                    if os.path.islink(tax_dest):
                        os.unlink(tax_dest)
                        print(f"{UI.GREEN}  └── 成功解除软链接，彻底保护源文件！{UI.RESET}")
                    elif os.path.exists(tax_dest):
                        shutil.rmtree(tax_dest)
                    return 0
                except Exception as e:
                    logger.error(f"解除 Taxonomy 软链接失败: {e}")
                    return 1
            actions.append(("安全隔离：解除 Taxonomy 软链接", unlink_tax, "k2_unlink_taxo"))

        # ------------------------------------
        # 步骤 8: 清理缓存
        # ------------------------------------
        cmd_clean = ['kraken2-build', '--threads', str(threads), '--clean', '--db', db_dir]
        actions.append(("清理缓存与冗余序列", cmd_clean, "k2_clean"))

        # === 统一断点执行调度器 ===
        total_steps = len(actions)
        for i, (desc, action, step_id) in enumerate(actions, 1):
            if cm.is_done(step_id):
                UI.print_skipped(i, total_steps, desc)
                logger.info(f"断点跳过: {desc}")
                continue

            UI.print_step(i, total_steps, desc)
            if callable(action):
                if action() != 0: return 1
            else:
                if execute_realtime(action, tool_name="Kraken2") != 0:
                    logger.error(f"Kraken2 建库在 [{desc}] 阶段失败！")
                    return 1
            cm.mark_done(step_id)
        return 0

    @staticmethod
    def bowtie2(fasta, out_prefix, threads=1, force=False):
        dir_path = os.path.dirname(out_prefix)
        os.makedirs(dir_path, exist_ok=True)
        cm = CheckpointManager(dir_path, force=force)
        
        if cm.is_done("bt2_build"):
            UI.print_skipped(1, 1, f"编译 Bowtie2 索引 -> {out_prefix}")
            return 0
            
        UI.print_step(1, 1, f"编译 Bowtie2 索引 -> {out_prefix}")
        if execute_realtime(['bowtie2-build', '--threads', str(threads), fasta, out_prefix], tool_name="Bowtie2") == 0:
            cm.mark_done("bt2_build")
            return 0
        return 1

    @staticmethod
    def hisat2(fasta, out_prefix, threads=1, force=False):
        dir_path = os.path.dirname(out_prefix)
        os.makedirs(dir_path, exist_ok=True)
        cm = CheckpointManager(dir_path, force=force)
        
        if cm.is_done("ht2_build"):
            UI.print_skipped(1, 1, f"编译 HISAT2 索引 -> {out_prefix}")
            return 0
            
        UI.print_step(1, 1, f"编译 HISAT2 索引 -> {out_prefix}")
        if execute_realtime(['hisat2-build', '-p', str(threads), fasta, out_prefix], tool_name="HISAT2") == 0:
            cm.mark_done("ht2_build")
            return 0
        return 1

    @staticmethod
    def minimap2(fasta, out_prefix, seq_type, threads=1, force=False):
        dir_path = os.path.dirname(out_prefix)
        os.makedirs(dir_path, exist_ok=True)
        cm = CheckpointManager(dir_path, force=force)
        
        step_id = f"mm2_{seq_type}"
        if cm.is_done(step_id):
            UI.print_skipped(1, 1, f"编译 Minimap2 索引 [{seq_type}] -> {out_prefix}")
            return 0
            
        cmd = ['minimap2', '-t', str(threads), '-d', out_prefix]
        if seq_type == 'dna-short': cmd += ['-x', 'sr']
        elif seq_type == 'rna-short': cmd += ['-x', 'splice:sr']
        elif seq_type == 'nanopore': cmd += ['-x', 'map-ont']
        elif seq_type == 'pacbio': cmd += ['-x', 'map-pb']
        cmd += [fasta]
        
        UI.print_step(1, 1, f"编译 Minimap2 索引 [{seq_type}] -> {out_prefix}")
        if execute_realtime(cmd, tool_name="Minimap2") == 0:
            cm.mark_done(step_id)
            return 0
        return 1

# ==========================================
# 6. 命令行解析与主干逻辑
# ==========================================
def parse_comma_list(string):
    return [s.strip() for s in string.split(',') if s.strip()]

def main():
    parser = argparse.ArgumentParser(description="一键化竞争性宿主全库构建工具 (支持断点续传与安全隔离)")
    
    # 核心与环境参数
    parser.add_argument('--tool', type=parse_comma_list, required=True, help="建库工具列表，逗号分隔 (如 kraken2,bowtie2)")
    parser.add_argument('-i', '--input', required=True, help="主宿主参考基因组 FASTA 路径")
    parser.add_argument('-o', '--output', required=True, help="输出的主目录 (如 host_db)")
    parser.add_argument('-t', '--threads', type=int, default=os.cpu_count(), help="CPU 线程数")
    
    # 专属建库参数
    parser.add_argument('--seq-type', type=parse_comma_list, default=['dna-short'], help="Minimap2: 预设类型，逗号分隔")
    parser.add_argument('--taxonomy', type=str, help="Kraken2: 本地 Taxonomy 目录，秒级软链接并在 clean 前安全拆卸")
    parser.add_argument('--taxid', type=int, help="Kraken2: (非模式物种专用) 强制为主宿主指定 NCBI TaxID")
    parser.add_argument('--k2-libs', type=parse_comma_list, 
                        default=["archaea", "bacteria", "plasmid", "fungi", "protozoa", "UniVec", "UniVec_Core"],
                        help="Kraken2: 需额外下载的官方库。设为 'none' 则不下载。")
    parser.add_argument('--add-library', type=parse_comma_list, default=[], help="Kraken2: 补充自定义 FASTA (逗号分隔)")
    
    # 运行与调试参数
    parser.add_argument('--force', action='store_true', help="无视历史断点，强制从头重新构建所有库！")
    parser.add_argument('-d', '--debug', action='store_true', help="开启详细日志")
    parser.add_argument('-l', '--log-file', type=str, help="将日志保存到文件")

    args = parser.parse_args()
    global logger
    logger = configure_logger(__name__, debug=args.debug, filename=args.log_file)

    valid_tools, valid_seq = ['kraken2', 'bowtie2', 'hisat2', 'minimap2'], ['dna-short', 'rna-short', 'nanopore', 'pacbio']
    for t in args.tool:
        if t not in valid_tools:
            logger.error(f"不支持的工具: {t}")
            sys.exit(1)
    for st in args.seq_type:
        if st not in valid_seq:
            logger.error(f"不支持的数据类型: {st}")
            sys.exit(1)

    valid, arg, chars = validate_args([args.input, args.output])
    if not valid:
        logger.error(f"输入包含非法路径字符: '{arg}'")
        sys.exit(1)
        
    if not os.path.isfile(args.input):
        logger.error(f"找不到主宿主输入的 FASTA: {args.input}")
        sys.exit(1)
        
    if 'kraken2' in args.tool:
        if args.taxonomy and not os.path.exists(os.path.join(args.taxonomy, 'names.dmp')):
            logger.error("Taxonomy 目录内缺少关键文件 (names.dmp)！")
            sys.exit(1)
        for extra_fa in args.add_library:
            if not os.path.isfile(extra_fa):
                logger.error(f"找不到补充的自定义 FASTA: {extra_fa}")
                sys.exit(1)

    check_dependencies(args.tool)

    tasks = []
    for tool in args.tool:
        if tool == 'minimap2':
            for st in args.seq_type: tasks.append((tool, st))
        else: tasks.append((tool, None))

    start_time_global = time.time()
    
    print("\n")
    print(f"  {UI.BOLD}🎯 All-in-One 竞争性宿主数据库建库向导{UI.RESET}")
    print(f"  {UI.BOLD}📌 主宿主参考 : {args.input}{UI.RESET}")
    print(f"  {UI.BOLD}📁 主存储路径 : {args.output}{UI.RESET}")
    print(f"  {UI.BOLD}🧬 TaxID 注入: {args.taxid if args.taxid else '未启用'}{UI.RESET}")
    if args.add_library:
        print(f"  {UI.BOLD}🧩 自定义补充 : {len(args.add_library)} 个 FASTA 文件{UI.RESET}")
    print(f"  {UI.BOLD}⚡ 断点续传   : {'禁用 (强制重跑)' if args.force else '已启用'}{UI.RESET}")
    print(f"  {UI.BOLD}📦 计划任务   : {len(tasks)} 项架构编译{UI.RESET}")

    results = {}
    base_dir = args.output

    for idx, (tool, st) in enumerate(tasks, 1):
        task_name = f"{tool.upper()}" + (f" ({st})" if st else "")
        UI.print_header(f"任务进度 [{idx}/{len(tasks)}] 🚀 构建 {task_name}")
        
        exit_code = 1
        start_time_task = time.time()
        
        if tool == 'kraken2':
            out_dir = os.path.join(base_dir, 'kraken2')
            exit_code = IndexBuilders.kraken2(
                args.input, out_dir, args.threads, 
                taxonomy_path=args.taxonomy, 
                k2_libs=args.k2_libs, taxid=args.taxid, 
                add_library=args.add_library, force=args.force
            )
            results['kraken2'] = exit_code
            
        elif tool == 'bowtie2':
            out_prefix = os.path.join(base_dir, 'bowtie2', 'host')
            exit_code = IndexBuilders.bowtie2(args.input, out_prefix, args.threads, args.force)
            results['bowtie2'] = exit_code
            
        elif tool == 'hisat2':
            out_prefix = os.path.join(base_dir, 'hisat2', 'host')
            exit_code = IndexBuilders.hisat2(args.input, out_prefix, args.threads, args.force)
            results['hisat2'] = exit_code
            
        elif tool == 'minimap2':
            out_mmi = os.path.join(base_dir, 'minimap2', f'host_{st}.mmi')
            exit_code = IndexBuilders.minimap2(args.input, out_mmi, st, args.threads, args.force)
            results[f'minimap2 ({st})'] = exit_code
            
        elapsed_task = (time.time() - start_time_task) / 60
        if exit_code == 0:
            print(f"\n{UI.GREEN}{UI.BOLD}✔ {task_name} 构建成功！(本次耗时: {elapsed_task:.2f} 分钟){UI.RESET}")
        else:
            print(f"\n{UI.RED}{UI.BOLD}✖ {task_name} 构建失败！(中止前耗时: {elapsed_task:.2f} 分钟){UI.RESET}")
            sys.exit(1)

    end_time_global = time.time()
    
    UI.print_header("建库流水线最终报告")
    for task, code in results.items():
        if code == 0: print(f"  {UI.GREEN}✅ {task.ljust(20)} [SUCCESS]{UI.RESET}")
        else: print(f"  {UI.RED}❌ {task.ljust(20)} [FAILED]{UI.RESET}")
            
    print(f"\n{UI.CYAN}{UI.BOLD}🎉 所有建库任务圆满完成！本次运行总耗时: {(end_time_global - start_time_global)/60:.2f} 分钟{UI.RESET}\n")

if __name__ == '__main__':
    main()
