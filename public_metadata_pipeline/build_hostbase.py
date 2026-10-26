#!/usr/bin/env python3
"""
独立的数据库与索引批量构建工具 (Kraken2, Bowtie2, HISAT2, Minimap2)
【极致优雅版】利用软链接实现 0 秒挂载 Taxonomy，并在 --clean 前自动安全脱穿隔离源文件。
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

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

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
def check_dependencies(tools, kunpeng_bin='kun_peng'):
    missing = []
    for tool in tools:
        executable = tool
        if tool == 'bowtie2': executable = 'bowtie2-build'
        elif tool == 'hisat2': executable = 'hisat2-build'
        elif tool == 'kraken2': executable = 'kraken2-build'
        elif tool == 'kunpeng': executable = kunpeng_bin
        found = shutil.which(executable) is not None or (
            os.path.isfile(executable) and os.access(executable, os.X_OK))
        if not found: missing.append(executable)
    if missing:
        logger.error(f"环境缺失！找不到以下命令: {', '.join(missing)}")
        sys.exit(1)

def validate_args(args):
    for arg in args:
        if type(arg) == str:
            probe = str(arg).replace('\\', '/')   # Windows 路径反斜杠归一
            chars = re.findall(r'[^a-zA-Z0-9\. !#=/:_-]', probe)
            if not re.match(r'^[a-zA-Z0-9\. !#=/:_-]+$|^(?![\s\S])', probe): return False, arg, chars
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
def _inject_kraken_taxid(src_fasta, dst_fasta, taxid):
    """在 FASTA 每条序列 ID 后追加 `|kraken:taxid|N`（Kraken2/Kun-peng 通用）。"""
    seq_count = 0
    taxid_tag = f"|kraken:taxid|{taxid}"
    with open(src_fasta, 'r') as in_f, open(dst_fasta, 'w') as out_f:
        for line in in_f:
            if line.startswith('>'):
                seq_count += 1
                parts = line.strip().split(maxsplit=1)
                new_id = f"{parts[0]}{taxid_tag}"
                out_f.write(f"{new_id} {parts[1]}\n" if len(parts) == 2 else f"{new_id}\n")
            else:
                out_f.write(line)
    return seq_count


def _kunpeng_shards(db_dir):
    return [p for p in glob.glob(os.path.join(db_dir, 'hash_*.k2d'))
            if re.match(r'^hash_\d+\.k2d$', os.path.basename(p))]


def _kunpeng_db_ready(db_dir):
    return all(os.path.isfile(os.path.join(db_dir, f))
               for f in ('hash_config.k2d', 'opts.k2d', 'taxo.k2d')) \
        and bool(_kunpeng_shards(db_dir))


def _reset_kunpeng_db(db_dir):
    """--force: 清掉 library/ 与所有派生文件，避免 add-library 累积重复序列。

    Kun-peng 的 add-library 按 added.md5 去重并 append seqid2taxid.map；
    只删 hash_*.k2d 而不删 library/map 会让重建把旧序列再算一遍。
    """
    for sub in ('library', 'prep'):
        p = os.path.join(db_dir, sub)
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)
    for pat in ('hash_*.k2d', 'chunk_*.k2', 'seqid2taxid.map',
                'taxo.k2d', 'opts.k2d', 'add_seqid2taxid*.map'):
        for p in glob.glob(os.path.join(db_dir, pat)):
            try:
                os.remove(p)
            except OSError:
                pass


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

    # ------------------------------------------------------------------
    # Kun-peng 宿主库
    # ------------------------------------------------------------------
    @staticmethod
    def kunpeng(fasta, db_dir, threads=1, taxonomy_path=None, taxid=None,
                add_library=None, hash_capacity='1G', force=False,
                kunpeng_bin='kun_peng'):
        """原生 Kun-peng 建库: add-library → build-db (estimate+chunk+build)。

        与 Kraken2 的差异:
          · 必须有 <db>/taxonomy/{nodes,names}.dmp (chunk 阶段生成 taxo.k2d 用);
          · FASTA 头同样使用 kraken:taxid|N 标签;
          · 产物是 hash_*.k2d + hash_config.k2d + opts.k2d + taxo.k2d,
            没有 Kraken2 的单文件 hash.k2d。
        """
        os.makedirs(db_dir, exist_ok=True)
        if force:
            _reset_kunpeng_db(db_dir)
        cm = CheckpointManager(db_dir, force=force)

        if _kunpeng_db_ready(db_dir) and not force:
            print(f"{UI.GREEN}  └── Kun-peng 库已就绪，跳过建库: {db_dir}{UI.RESET}")
            return 0

        actions = []

        # 步骤 1: 挂载 taxonomy (copy；Kun-peng chunk 阶段直接读 <db>/taxonomy)
        if not taxonomy_path:
            logger.error("Kun-peng 原生建库需要 --taxonomy (含 nodes.dmp / names.dmp)")
            return 1
        tax_src = os.path.abspath(taxonomy_path)

        def copy_tax():
            tax_dst = os.path.join(db_dir, 'taxonomy')
            os.makedirs(tax_dst, exist_ok=True)
            copied = []
            for name in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
                src = os.path.join(tax_src, name)
                if os.path.isfile(src):
                    shutil.copyfile(src, os.path.join(tax_dst, name))
                    copied.append(name)
            if not {'nodes.dmp', 'names.dmp'} <= set(copied):
                logger.error("taxonomy 目录缺少 nodes.dmp 或 names.dmp")
                return 1
            print(f"{UI.GREEN}  └── taxonomy 已复制: {', '.join(copied)}{UI.RESET}")
            return 0
        actions.append((f"复制 NCBI taxonomy → {db_dir}/taxonomy", copy_tax, "kp_taxo"))

        # 步骤 2: 主宿主 taxid 注入 (可选；已带 taxid 标签的 FASTA 可省)
        final_fasta = fasta
        if taxid:
            prep_dir = os.path.join(db_dir, 'prep')
            os.makedirs(prep_dir, exist_ok=True)
            final_fasta = os.path.join(prep_dir, "host_with_taxid.fasta")

            def inject_taxid():
                try:
                    n = _inject_kraken_taxid(fasta, final_fasta, taxid)
                    print(f"{UI.GREEN}  └── 成功注入 TaxID {taxid}，修改 {n} 条序列！{UI.RESET}")
                    return 0
                except Exception as e:
                    logger.error(f"注入 TaxID 标签失败: {e}")
                    return 1
            actions.append((f"为主宿主注入 TaxID: {taxid}", inject_taxid, "kp_taxid_inject"))

        # 步骤 3: 补充自定义 FASTA
        if add_library:
            for extra_fa in add_library:
                abs_extra_fa = os.path.abspath(extra_fa)
                base = os.path.basename(abs_extra_fa)
                cmd = [kunpeng_bin, 'add-library', '--db', db_dir, '-i', abs_extra_fa]
                actions.append((f"补充自定义 FASTA: {base}", cmd, f"kp_add_extra_{base}"))

        # 步骤 4: 主宿主入库
        cmd_add = [kunpeng_bin, 'add-library', '--db', db_dir, '-i', final_fasta]
        actions.append(("主宿主 FASTA 加入 library", cmd_add, "kp_add_host"))

        # 步骤 5: estimate + chunk + build (Kun-peng 内部自动串联)
        # seqkmer 并行读入要求 n_threads > 2 (0.7.11 实测 -p 2 panic)
        kp_threads = max(3, int(threads))
        if kp_threads != threads:
            logger.warning(f"Kun-peng 要求线程数 > 2, 已从 {threads} 抬到 {kp_threads}")
        cmd_build = [kunpeng_bin, 'build-db', '--db', db_dir,
                     '--hash-capacity', str(hash_capacity), '-p', str(kp_threads)]
        actions.append((f"编译 hash 分片 (hash-capacity={hash_capacity}, 极其耗时)",
                        cmd_build, "kp_build_db"))

        total_steps = len(actions)
        for i, (desc, action, step_id) in enumerate(actions, 1):
            if cm.is_done(step_id):
                UI.print_skipped(i, total_steps, desc)
                logger.info(f"断点跳过: {desc}")
                continue
            UI.print_step(i, total_steps, desc)
            if callable(action):
                if action() != 0:
                    return 1
            else:
                if execute_realtime(action, tool_name="Kun-peng") != 0:
                    logger.error(f"Kun-peng 建库在 [{desc}] 阶段失败！")
                    return 1
            cm.mark_done(step_id)

        if not _kunpeng_db_ready(db_dir):
            logger.error(f"Kun-peng 库校验失败 (缺 hash_config/taxo/opts 或 hash 分片): {db_dir}")
            return 1
        n = len(_kunpeng_shards(db_dir))
        print(f"{UI.GREEN}{UI.BOLD}  ✔ Kun-peng 库就绪: {n} 个 hash 分片 → {db_dir}{UI.RESET}")
        return 0

    @staticmethod
    def kunpeng_from_kraken2(k2_dir, hash_capacity='1G', force=False,
                             kunpeng_bin='kun_peng'):
        """把已建成的 Kraken2 库原地转成 Kun-peng 分片库 (保留 hash.k2d)。

        `kun_peng hashshard` 只追加 hash_*.k2d + hash_config.k2d (并补 taxo/opts)，
        不删除 Kraken2 文件，因此转换后的目录同时支持 kraken2 与 kun_peng。
        """
        if _kunpeng_db_ready(k2_dir) and not force:
            print(f"{UI.GREEN}  └── Kun-peng 分片已存在，跳过转换: {k2_dir}{UI.RESET}")
            return 0
        if not force and os.path.isfile(os.path.join(k2_dir, 'hash_config.k2d')):
            logger.error("检测到半转换状态 (hash_config.k2d 存在但分片不完整); "
                         "请加 --force 重新转换")
            return 1
        if force:
            for p in glob.glob(os.path.join(k2_dir, 'hash_config.k2d')) + \
                     glob.glob(os.path.join(k2_dir, 'hash_*.k2d')):
                try:
                    os.remove(p)
                except OSError:
                    pass
        if not os.path.isfile(os.path.join(k2_dir, 'hash.k2d')):
            logger.error(f"找不到 Kraken2 库 hash.k2d: {k2_dir}")
            return 1
        cmd = [kunpeng_bin, 'hashshard', '--db', k2_dir,
               '--hash-capacity', str(hash_capacity)]
        if execute_realtime(cmd, tool_name="Kun-peng") != 0:
            logger.error("Kun-peng hashshard 转换失败")
            return 1
        if not _kunpeng_db_ready(k2_dir):
            logger.error(f"hashshard 结束但库校验失败: {k2_dir}")
            return 1
        print(f"{UI.GREEN}{UI.BOLD}  ✔ Kraken2 库已转为双后端库 "
              f"({len(_kunpeng_shards(k2_dir))} 个分片): {k2_dir}{UI.RESET}")
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
    parser.add_argument('--tool', type=parse_comma_list, required=True, help="建库工具列表，逗号分隔 (如 kraken2,bowtie2,kunpeng)")
    parser.add_argument('-i', '--input', required=True, help="主宿主参考基因组 FASTA 路径")
    parser.add_argument('-o', '--output', required=True, help="输出的主目录 (如 host_db)")
    parser.add_argument('-t', '--threads', type=int, default=os.cpu_count(), help="CPU 线程数")
    
    # 专属建库参数
    parser.add_argument('--seq-type', type=parse_comma_list, default=['dna-short'], help="Minimap2: 预设类型，逗号分隔")
    parser.add_argument('--taxonomy', type=str, help="Kraken2/Kun-peng: 本地 Taxonomy 目录 (含 nodes.dmp / names.dmp)")
    parser.add_argument('--taxid', type=int, help="Kraken2/Kun-peng: (非模式物种专用) 强制为主宿主指定 NCBI TaxID")
    parser.add_argument('--k2-libs', type=parse_comma_list, 
                        default=["archaea", "bacteria", "plasmid", "fungi", "protozoa", "UniVec", "UniVec_Core"],
                        help="Kraken2: 需额外下载的官方库。设为 'none' 则不下载。")
    parser.add_argument('--add-library', type=parse_comma_list, default=[], help="Kraken2/Kun-peng: 补充自定义 FASTA (逗号分隔)")
    parser.add_argument('--hash-capacity', default='1G',
                        help="Kun-peng: hash 分片容量 (默认 1G ≈ 4GiB/片; 同时用于 Kraken2→Kun-peng 转换)")
    parser.add_argument('--kunpeng-bin', default=os.environ.get('MMPV_KUNPENG_BIN', 'kun_peng'),
                        help="Kun-peng 可执行文件 (默认 $MMPV_KUNPENG_BIN 或 PATH 中的 kun_peng)")
    parser.add_argument('--kunpeng-native', action='store_true',
                        help="同时要求 kraken2,kunpeng 时默认用 hashshard 把 Kraken2 库原地转成分片库 "
                             "(省一半磁盘); 加此开关则从 FASTA 原生重建独立的 kunpeng/ 库")
    
    # 运行与调试参数
    parser.add_argument('--force', action='store_true', help="无视历史断点，强制从头重新构建所有库！")
    parser.add_argument('-d', '--debug', action='store_true', help="开启详细日志")
    parser.add_argument('-l', '--log-file', type=str, help="将日志保存到文件")

    args = parser.parse_args()
    global logger
    logger = configure_logger(__name__, debug=args.debug, filename=args.log_file)

    valid_tools, valid_seq = ['kraken2', 'bowtie2', 'hisat2', 'minimap2', 'kunpeng'], ['dna-short', 'rna-short', 'nanopore', 'pacbio']
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

    kunpeng_native = args.kunpeng_native or 'kraken2' not in args.tool
    need_taxonomy = 'kraken2' in args.tool or ('kunpeng' in args.tool and kunpeng_native)
    if need_taxonomy and args.taxonomy:
        if not os.path.exists(os.path.join(args.taxonomy, 'names.dmp')):
            logger.error("Taxonomy 目录内缺少关键文件 (names.dmp)！")
            sys.exit(1)
        if 'kunpeng' in args.tool and not os.path.exists(os.path.join(args.taxonomy, 'nodes.dmp')):
            logger.error("Kun-peng 原生建库需要 Taxonomy 目录内的 nodes.dmp！")
            sys.exit(1)
    if 'kunpeng' in args.tool and kunpeng_native and not args.taxonomy:
        logger.error("Kun-peng 原生建库需要 --taxonomy (nodes.dmp / names.dmp)")
        sys.exit(1)
    for extra_fa in args.add_library:
        if not os.path.isfile(extra_fa):
            logger.error(f"找不到补充的自定义 FASTA: {extra_fa}")
            sys.exit(1)
    if 'kunpeng' in args.tool and kunpeng_native and args.k2_libs \
            and args.k2_libs[0].lower() != 'none':
        logger.warning("--k2-libs 仅用于 Kraken2 官方库; Kun-peng 原生建库会忽略它 "
                       "(需要背景库请用 --add-library 提供 FASTA)")

    check_dependencies(args.tool, kunpeng_bin=args.kunpeng_bin)

    # kunpeng 的 hashshard 依赖 kraken2 库：无论用户书写顺序如何，先建 kraken2
    if 'kunpeng' in args.tool and 'kraken2' in args.tool and not args.kunpeng_native:
        args.tool = ['kraken2'] + [t for t in args.tool if t != 'kraken2']

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
    if 'kunpeng' in args.tool:
        print(f"  {UI.BOLD}🐧 Kun-peng   : {args.kunpeng_bin} | hash-capacity={args.hash_capacity}"
              f" | {'原生重建' if kunpeng_native else '由 kraken2 库 hashshard 转换'}{UI.RESET}")
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
            
        elif tool == 'kunpeng':
            if 'kraken2' in args.tool and not args.kunpeng_native:
                # 已有 Kraken2 库 → 原地 hashshard, 一个目录同时支持两种后端
                k2_dir = os.path.join(base_dir, 'kraken2')
                exit_code = IndexBuilders.kunpeng_from_kraken2(
                    k2_dir, hash_capacity=args.hash_capacity,
                    force=args.force, kunpeng_bin=args.kunpeng_bin)
                results['kunpeng (from kraken2)'] = exit_code
            else:
                out_dir = os.path.join(base_dir, 'kunpeng')
                exit_code = IndexBuilders.kunpeng(
                    args.input, out_dir, args.threads,
                    taxonomy_path=args.taxonomy, taxid=args.taxid,
                    add_library=args.add_library,
                    hash_capacity=args.hash_capacity,
                    force=args.force, kunpeng_bin=args.kunpeng_bin)
                results['kunpeng'] = exit_code

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
