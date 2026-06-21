#!/usr/bin/env python3
"""
batch_virema_dvg.py — 高通量 DVG 与病毒重组全自动流水线
【大满配终局版：动态内存盘 / 纯净洗涤 / 动态进度轨 / 自动出图 / Coverage脱水 / 结果保全】
"""

import argparse
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import traceback
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import colorlog
import pandas as pd
from tqdm import tqdm

# ==========================================
# 1. 日志与系统工具
# ==========================================
def setup_logging(verbose: bool = False) -> logging.Logger:
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        "%(log_color)s%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
        log_colors={"DEBUG": "cyan", "INFO": "green", "WARNING": "yellow", "ERROR": "red", "CRITICAL": "bold_red"},
    ))
    log = colorlog.getLogger("virema_batch")
    log.handlers.clear()
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    return log

logger = setup_logging()

def check_psutil():
    try:
        import psutil
    except ImportError:
        logger.critical("❌ 为了探测内存安全防爆水位，请先安装: pip install psutil")
        sys.exit(1)

def safe_name(s: str, max_len: int = 100) -> str:
    s = str(s)
    s = re.sub(r'[^A-Za-z0-9\-.]', '_', s)
    s = re.sub(r'_+', '_', s)
    return s.strip('_.')[:max_len]

def run_cmd(cmd: str, log_path: str = None, check: bool = True):
    full_cmd = f"set -o pipefail; {cmd}"
    result = subprocess.run(full_cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    
    if log_path:
        lp = Path(log_path)
        lp.parent.mkdir(parents=True, exist_ok=True)
        with open(lp, "a", encoding="utf-8") as f:
            f.write(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] CMD: {cmd}\n")
            f.write(f"EXIT_CODE: {result.returncode}\n")
            if result.stdout: f.write(f"--- STDOUT ---\n{result.stdout.strip()}\n")
            if result.stderr: f.write(f"--- STDERR ---\n{result.stderr.strip()}\n")
            f.write("-" * 60 + "\n")
            
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    return result


# ==========================================
# 2. 🌟 绝对提纯 FASTQ 隔离海关 (带动态进度仪)
# ==========================================
def sanitize_fastq_stream(files_dict, out_path, is_fasta, sample, max_workers, log_file):
    """
    底层核心海关：
    1. 计算读取物理进度并在终端投影动态加载条。
    2. 无情剿杀一切混有 N 或不明符号的 Reads，阻止 Bowtie C++ 内核爆破。
    3. 赋予唯一标识后缀 _M / _R1 / _R2 避免 ViReMa 字典击穿。
    """
    current_identity = multiprocessing.current_process()._identity
    worker_id = current_identity[0] if current_identity else 1
    pos = (worker_id % max_workers) + 1
    
    # 极其苛刻的纯净序列验证：只要包含非 a,t,g,c,A,T,G,C 的字符（包含 N!），直接视作脏数据火化
    seq_dirty_checker = re.compile(r'[^ATGCatgc]')
    qual_cleaner = re.compile(r'[^\x21-\x7E]') 

    total_bytes = sum(Path(f).stat().st_size for f in files_dict.values() if Path(f).exists())
    
    with open(out_path, 'w') as out_f:
        desc_str = f"[{sample[:10]:<10}] 净化重组"
        with tqdm(total=total_bytes, desc=desc_str, position=pos, leave=False, 
                  unit='B', unit_scale=True, colour='green') as pbar:
            
            for suffix, fpath in files_dict.items():
                if not Path(fpath).exists() or Path(fpath).stat().st_size == 0:
                    continue
                
                opener = gzip.open if str(fpath).endswith('.gz') else open
                try:
                    with opener(fpath, 'rt', encoding='utf-8', errors='replace') as in_f:
                        chunk_bytes = 0
                        if is_fasta:
                            header = ""
                            seq_chunks = []
                            for line in in_f:
                                chunk_bytes += len(line)
                                if chunk_bytes > 500000:
                                    pbar.update(chunk_bytes)
                                    chunk_bytes = 0

                                line = line.strip()
                                if not line: continue
                                if line.startswith('>'):
                                    if header and seq_chunks:
                                        final_seq = "".join(seq_chunks)
                                        # 过滤所有混有 N 且过短的序列
                                        if not seq_dirty_checker.search(final_seq) and len(final_seq) >= 15:
                                            out_f.write(f"{header}\n{final_seq}\n")
                                            
                                    parts = line.split(maxsplit=1)
                                    header = parts[0] + suffix + (" " + parts[1] if len(parts)>1 else "")
                                    seq_chunks = []
                                else:
                                    seq_chunks.append(line)
                                    
                            if header and seq_chunks:
                                final_seq = "".join(seq_chunks)
                                if not seq_dirty_checker.search(final_seq) and len(final_seq) >= 15:
                                    out_f.write(f"{header}\n{final_seq}\n")
                        else:
                            for line in in_f:
                                chunk_bytes += len(line)
                                if chunk_bytes > 500000:
                                    pbar.update(chunk_bytes)
                                    chunk_bytes = 0

                                header = line.strip()
                                if not header.startswith('@'): continue
                                
                                seq = in_f.readline().strip()
                                plus = in_f.readline().strip()
                                qual = in_f.readline().strip()
                                chunk_bytes += len(seq) + len(plus) + len(qual) + 3
                                
                                # 双保险检查等长，同时彻底剔除含 N 的读取
                                if seq and len(seq) == len(qual) and len(seq) >= 15:
                                    if not seq_dirty_checker.search(seq):
                                        qual = qual_cleaner.sub('#', qual)
                                        parts = header.split(maxsplit=1)
                                        new_header = parts[0] + suffix + (" " + parts[1] if len(parts)>1 else "")
                                        out_f.write(f"{new_header}\n{seq}\n+\n{qual}\n")
                        
                        if chunk_bytes > 0: pbar.update(chunk_bytes)
                except Exception as e:
                    with open(log_file, "a") as lf:
                        lf.write(f"⚠️ 解构 {fpath} 中遇大尺度物理损坏，安检强行隔离。 ({str(e)})\n")


# ==========================================
# 3. 核心大闭环：BBMerge -> 洗涤 -> 沙盒 ViReMa
# ==========================================
def worker_virema(args):
    (sample, virus, r1_str, r2_str, is_single, ref_fa, out_dir_str, 
     sam_dir_str, seed, mindel, defuzz, virema_path, threads, resume, keep_temp, log_file, use_shm, max_jobs) = args
    
    out_dir = Path(out_dir_str)
    sam_dir = Path(sam_dir_str)
    out_dir.mkdir(parents=True, exist_ok=True)
    sam_dir.mkdir(parents=True, exist_ok=True)
    
    sam_basename = f"{sample}_{virus}_ViReMa.sam"
    output_tag = f"{sample}_{virus}" # 为该样本的txt和bed打下不灭的姓名烙印
    final_sam_full_path = sam_dir / sam_basename
    
    # 续传检验：判定SAM必须大于 100 Bytes
    if resume and final_sam_full_path.exists() and final_sam_full_path.stat().st_size > 100:
        return True

    r1 = Path(r1_str) if r1_str else None
    r2 = Path(r2_str) if r2_str else None
    is_fasta = True if r1 and any(ext in r1.name.lower() for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']) else False
    fmt_ext = ".fa" if is_fasta else ".fq"
    t_io = min(4, threads)
    
    shm_mount = Path("/dev/shm")
    shm_active = False
    
    import psutil
    if use_shm and shm_mount.exists() and shm_mount.is_mount():
        shm_usage = psutil.disk_usage('/dev/shm')
        if (shm_usage.free / (1024**3)) > 10.0:
            shm_active = True
            tmp_dir = shm_mount / f"virema_ds_{sample}_{virus}"
            with open(log_file, "a") as lf: lf.write(f"\n[⚡ 架构预警] 高速 RAM ({shm_usage.free / 1024**3:.1f}GB 充足)，沙盒迁移至内存盘中...\n")
        else:
            with open(log_file, "a") as lf: lf.write(f"\n[⚠️ 容量预警] RAM盘空间预警 (仅剩 {shm_usage.free / 1024**3:.1f}GB)，安全防爆退回物理系统碟...\n")
            tmp_dir = out_dir / "virema_sandbox"
    else:
        tmp_dir = out_dir / "virema_sandbox"
    
    tmp_dir.mkdir(parents=True, exist_ok=True)
    
    tmp_merged = tmp_dir / f"merged{fmt_ext}"
    tmp_unmerged_r1 = tmp_dir / f"unmerged_R1{fmt_ext}"
    tmp_unmerged_r2 = tmp_dir / f"unmerged_R2{fmt_ext}"
    tmp_input_virema = tmp_dir / f"input_virema{fmt_ext}"
    
    if tmp_input_virema.exists(): tmp_input_virema.unlink()
        
    run_success = False
    try:
        files_to_merge_dict = {}
        # A) 双端物理重叠
        if not is_single and r1 and r2 and r1.exists() and r2.exists():
            cmd_merge = (
                f"bbmerge.sh in1='{r1}' in2='{r2}' "
                f"out='{tmp_merged}' outu1='{tmp_unmerged_r1}' outu2='{tmp_unmerged_r2}' "
                f"t={t_io} 2>>'{log_file}'"
            )
            run_cmd(cmd_merge, log_path=log_file, check=False)
            files_to_merge_dict = {'_M': tmp_merged, '_R1': tmp_unmerged_r1, '_R2': tmp_unmerged_r2}
        else:
            files_to_merge_dict = {'_S': r1}

        # B) 🚀 引导进入最高等级带长进度条的钛合金安检筛查过滤网
        sanitize_fastq_stream(files_to_merge_dict, tmp_input_virema, is_fasta, sample, max_jobs, log_file)

        if not tmp_input_virema.exists() or tmp_input_virema.stat().st_size == 0:
            with open(log_file, "a") as lf: lf.write("⚠️ 数据不合标准被全数销毁，退出后续处理。\n")
            run_success = True
            return True
            
        fasta_flag = "-Fasta " if is_fasta else ""
        
        # C) 运行 ViReMa：写入 -Overwrite 并锁定 --Output_Tag，让产物具备防混识别前缀！
        virema_cmd = (
            f"python '{virema_path}' '{ref_fa}' '{tmp_input_virema}' '{sam_basename}' "
            f"{fasta_flag}--Seed {seed} --MicroInDel_Length {mindel} --Defuzz {defuzz} "
            f"--p {threads} -BED -Overwrite --Output_Tag '{output_tag}' --Output_Dir '{tmp_dir}'"
        )
        run_cmd(virema_cmd, log_path=log_file, check=True)
        
        # D) 暴力提取并收押沙盒中所有的战利品
        # (1) 剥离骨架大文件 SAM
        produced_sam = tmp_dir / sam_basename
        if produced_sam.exists(): shutil.move(produced_sam, final_sam_full_path)
            
        # (2) 深喉提取 BED_Files 目录和外部一切带样本名的 txt 指标
        dest_bed_dir = out_dir / "bed_results"
        dest_bed_dir.mkdir(exist_ok=True)
        
        bed_res_dir = tmp_dir / "BED_Files"
        if bed_res_dir.exists():
            for f in bed_res_dir.glob("*.bed"): shutil.move(f, dest_bed_dir)
            for f in bed_res_dir.glob("*.bedgraph"): shutil.move(f, dest_bed_dir)
            for f in bed_res_dir.glob("*.BEDPE"): shutil.move(f, dest_bed_dir)
            
        for f in tmp_dir.glob("*_Results.txt"): shutil.move(f, dest_bed_dir)
        for f in tmp_dir.glob("*_Insertions.txt"): shutil.move(f, dest_bed_dir)
        for f in tmp_dir.glob("*_Substitutions.txt"): shutil.move(f, dest_bed_dir)
        for f in tmp_dir.glob("*_Micro*.txt"): shutil.move(f, dest_bed_dir)

        run_success = True
        return True
        
    except Exception as e:
        err = f"\n[Python 异常拦截] 样本 {sample} 绝缘崩溃:\n{traceback.format_exc()}\n"
        with open(log_file, "a") as lf: lf.write(err)
        return False
    finally:
        # E) ✅ 无情引爆沙盒，如果启动了内存盘则强行连根拔起释放资源
        if tmp_dir.exists():
            if shm_active or (run_success and not keep_temp):
                shutil.rmtree(tmp_dir)


# ==========================================
# 4. 全局指挥雷达枢纽
# ==========================================
class BatchViReMaPipeline:
    def __init__(self, args):
        self.args = args
        self.out = Path(args.output_dir)
        self.d_refs = self.out / "individual_refs"
        self.d_results = self.out / "virema_results"
        self.d_sams = self.out / "virema_sams"
        self.d_all_beds = self.out / "all_gathered_beds"
        
        for d in [self.d_refs, self.d_results, self.d_sams, self.d_all_beds]:
            d.mkdir(parents=True, exist_ok=True)
        
        self.df = None
        self.ref_map = {}
        self.sample_files = {}

    def run(self):
        check_psutil()
        logger.info("=" * 60)
        mode = "【 内存极速狂飙模式 (自动防爆) 】" if getattr(self.args, 'shm', False) else "【 标准硬盘避险模式 】"
        logger.info(f"🚀 终极霸主·批量 DVG 发掘管线引擎启动！ {mode}")
        logger.info("=" * 60)
        
        self.check_dependencies()
        self.load_summary()
        self.extract_virus_fastas()
        
        logger.info("⚡ 为了避开并发死锁陷阱，向机器内核提前下发建库指令...")
        self.prebuild_all_indexes()

        self.find_reads()
        self.print_dataset_statistics()
        
        self.run_virema()
        self.run_auto_r_report()
        
        logger.info("\n" + "=" * 60)
        logger.info(f"✨ 星舰平稳靠岸，全部科研绘档直落: {self.out}")

    def prebuild_all_indexes(self):
        for vd, fpath in self.ref_map.items():
            if not (Path(fpath).parent / f"{Path(fpath).name}.1.ebwt").exists() and \
               not (Path(fpath).parent / f"{Path(fpath).name}.1.bt2").exists():
                import glob
                base_dir = str(Path(fpath).parent)
                ebwt_files = glob.glob(f"{fpath}*.ebwt") + glob.glob(f"{fpath}*.bt2")
                if not ebwt_files:
                    logger.debug(f"📐 正在为靶向结构 {vd} 前期搭建索引架构...")
                    bcmd = f"bowtie-build '{fpath}' '{fpath}'"
                    run_cmd(bcmd, check=False)

    def print_dataset_statistics(self):
        total_samples = len(self.df['Sample'].unique())
        matched_samples = len(self.sample_files)
        total_viruses = len(self.ref_map)

        total_bytes = sum(
            Path(p).stat().st_size for sf in self.sample_files.values() 
            for p in (sf['r1'], sf['r2']) if p and Path(p).exists()
        )
        single_ends = sum(1 for sf in self.sample_files.values() if sf['is_single'])
        paired_ends = len(self.sample_files) - single_ends

        logger.info("\n📊 " + "=" * 54)
        logger.info(" 🔍 [ 预处理数据扫描快照 ]")
        logger.info(f"   ▶ 图谱导入总控 : {total_samples} 个样本池 / 侦测靶位 {total_viruses} 种病毒")
        logger.info(f"   ▶ 获取实际火控 : 挂载弹药舱 {matched_samples} 个 (双端: {paired_ends} | 单端: {single_ends})")
        if total_samples > matched_samples: logger.warning(f"   ▶ 发现幽灵遗落数据 : {total_samples - matched_samples} 个脱离阵列！")
        logger.info(f"   ▶ 总载重核酸吨位 : {total_bytes / (1024**3):.2f} GB")
        logger.info("=" * 58 + "\n")

    def check_dependencies(self):
        missing = []
        if not shutil.which("bbmerge.sh"): missing.append("bbmerge.sh")
        if not shutil.which("bowtie-build") and not shutil.which("bowtie2-build") and not shutil.which("bwa"): 
            missing.append("bowtie-build (依赖引擎)")
            
        vp = Path(self.args.virema_script)
        if vp.exists() and vp.is_dir():
            py_files = list(vp.glob("ViReMa*.py"))
            if py_files:
                self.args.virema_script = str(py_files[0])
            else:
                missing.append(f"❌ 目录 '{vp}' 内挖不着 ViReMa.py 大脑！")
        elif not vp.exists():
            missing.append(f"❌ ViReMa 脚本本体已遗失: {vp}")
            
        if missing:
            logger.error(f"❌ 核心前置武装组件瘫痪:\n" + "\n".join(f" - {m}" for m in missing))
            sys.exit(1)

    def load_summary(self):
        sep = "," if Path(self.args.summary).suffix.lower() == ".csv" else "\t"
        df = pd.read_csv(self.args.summary, sep=sep, dtype=str)
        df['Sample'] = df['Sample'].apply(lambda x: re.sub(r'(?i)(_clean|_trimmed|_filtered|_val|_fastp)+$', '', str(x).strip()))
        
        for col in ['Rep_Accession', 'Accession', 'Virus']:
            if col in df.columns: 
                df['Virus_Acc'] = df[col]
                break
        for col in ['Taxonomy', 'Adjusted_Species', 'Species_NCBI']:
            if col in df.columns:
                df['Taxonomy'] = df[col]
                break
                
        if 'Virus_Acc' not in df.columns:
            logger.error("❌ Summary 表必须蕴含病毒靶序列(Accession)！")
            sys.exit(1)
            
        if 'Taxonomy' not in df.columns: df['Taxonomy'] = "Unannotated"
        
        # 👑 执行 Coverage 深水过滤
        if self.args.min_cov > 0:
            cov_cols = ['Rep_Coverage(%)', 'Genome coverage', 'Coverage(%)', 'Coverage', 'Cov']
            cov_col_found = next((col for col in cov_cols if col in df.columns), None)
            if cov_col_found:
                df[cov_col_found] = df[cov_col_found].astype(str).str.replace('%', '').apply(pd.to_numeric, errors='coerce')
                filtered_len = len(df) - len(df[df[cov_col_found] >= self.args.min_cov])
                df = df[df[cov_col_found] >= self.args.min_cov]
                if filtered_len > 0: logger.info(f"✂️  Coverage清道夫：移除 {filtered_len} 个丰度低于红线 {self.args.min_cov}% 的标的。")
            else:
                logger.warning(f"⚠️ 表格未搜寻到 Coverage 属性，切除过滤已被跳过。")
                
        self.df = df.drop_duplicates(subset=['Sample', 'Virus_Acc'])

    def extract_virus_fastas(self):
        target_viruses = set(self.df['Virus_Acc'].unique().tolist())
        seq_buf = []; vid_cur = None
        def _save_fasta():
            if vid_cur in target_viruses and seq_buf:
                ref_fa = self.d_refs / f"{safe_name(vid_cur)}.fasta"
                if not ref_fa.exists():
                    with open(ref_fa, "w") as f: f.write(f">{vid_cur}\n" + "".join(seq_buf) + "\n")
                self.ref_map[vid_cur] = str(ref_fa)

        with open(self.args.reference, "r") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    _save_fasta(); vid_cur = line[1:].split()[0]; seq_buf = []
                else: 
                    seq_buf.append(line)
        _save_fasta()

    def find_reads(self):
        rd = Path(self.args.reads_dir)
        all_reads = [f for f in rd.rglob("*") if f.is_file() and any(ext in f.name.lower() for ext in ['.fq', '.fastq', '.fa', '.fasta', '.gz'])]
        
        for sname in self.df['Sample'].unique():
            s_clean = sname.lower()
            matched = []
            for f in all_reads:
                fn_clean = f.name.lower().replace('_', ' ').replace('-', ' ').replace('.', ' ')
                sc_spaced = s_clean.replace('_', ' ').replace('-', ' ').replace('.', ' ')
                if re.search(r'\b' + re.escape(sc_spaced) + r'\b', fn_clean) or f.name.lower().startswith(s_clean + "_") or f.name.lower().startswith(s_clean + "."):
                    matched.append(f)
            
            matched = list(set(matched))
            if not matched: continue
                
            r1, r2, is_single = None, None, True
            for f in matched:
                nl = f.name.lower()
                if any(x in nl for x in ['_r2', '_2.', '.r2', '_2_']): r2 = f
                elif any(x in nl for x in ['_r1', '_1.', '.r1', '_1_']): r1 = f
                else: 
                    if not r1: r1 = f
            
            if r1 and r2: is_single = False
            self.sample_files[sname] = {'r1': str(r1), 'r2': str(r2) if r2 else None, 'is_single': is_single}

    def run_virema(self):
        tasks = []
        for _, r in self.df.iterrows():
            s, v, tax = r["Sample"], r["Virus_Acc"], r["Taxonomy"]
            if s not in self.sample_files or v not in self.ref_map: continue
            
            sf = self.sample_files[s]
            out_dir = str(self.d_results / safe_name(tax) / f"{safe_name(s)}_{safe_name(v)}")
            log_file = str(self.d_results / safe_name(tax) / f"{safe_name(s)}_{safe_name(v)}" / "virema.log")
            
            tasks.append((
                s, v, sf['r1'], sf['r2'], sf['is_single'], self.ref_map[v], out_dir, str(self.d_sams),
                self.args.seed, self.args.mindel, self.args.defuzz, self.args.virema_script,
                self.args.threads, self.args.resume, getattr(self.args, 'keep', False), log_file, getattr(self.args, 'shm', False), self.args.jobs
            ))
            
        logger.info(f"🚀 发力总牵引：锁定全部火力核心并发发射 (Jobs: {self.args.jobs}) >>>")
        
        with tqdm(total=len(tasks), desc="⭕ 巨构阵列流转总体进度", colour='cyan', position=0, maxinterval=1) as pbar_total:
            with ProcessPoolExecutor(max_workers=self.args.jobs) as ex:
                futures = {ex.submit(worker_virema, t): t for t in tasks}
                ok, fail = 0, 0
                for fut in as_completed(futures):
                    sample, virus = futures[fut][:2]
                    try:
                        if fut.result(): ok += 1
                        else: fail += 1; logger.warning(f"DVG task failed: {sample}/{virus}")
                    except Exception as e:
                        fail += 1; logger.error(f"DVG task crashed: {sample}/{virus} - {e}")
                    pbar_total.update(1)
                logger.info(f"DVG tasks: {ok} OK, {fail} failed")

    def run_auto_r_report(self):
        # Look for R script: own directory first, then ViReMa install path
        r_script = Path(__file__).parent / "virema_summary_report.R"
        if not r_script.exists():
            r_script = Path(self.args.virema_script).parent / "virema_summary_report.R"
        if not r_script.exists():
            logger.warning("  virema_summary_report.R not found, skipping R visualization")
            return
            
        logger.info("\n🎨 触发 R 渲染大阵：调用云端分析网并准备发表级绘图...")
        
        # 🚀 修复核心：放开权限！把 bed, bedgraph, BEDPE, txt 统统打包带走！
        gather_files = []
        gather_files.extend(list(self.out.rglob("bed_results/*.bed")))
        gather_files.extend(list(self.out.rglob("bed_results/*.bedgraph")))
        gather_files.extend(list(self.out.rglob("bed_results/*.txt")))
        gather_files.extend(list(self.out.rglob("bed_results/*.BEDPE")))
        
        if not gather_files:
            return
            
        for b in gather_files: 
            shutil.copy(b, self.d_all_beds)
            
        r_out_dir = self.out / "Summary_Analysis_Report"
        r_cmd = f"Rscript '{r_script}' -i '{self.d_all_beds}' -o '{r_out_dir}' --auto_annotate"
        try:
            run_cmd(r_cmd, log_path=str(self.out / "r_summary_run.log"))
            logger.info(f"✅ R 画匠圆满谢幕！成品图集直供于: {r_out_dir}")
            
            if (self.d_sams / "BED_Files").exists(): 
                shutil.rmtree(self.d_sams / "BED_Files")
                
        except Exception as e:
            logger.error(f"❌ R 节点渲染退出，可检视 r_summary_run.log。")

if __name__ == "__main__":
    if sys.platform == "win32":
        import colorama
        colorama.init()

    parser = argparse.ArgumentParser(description="独立批量调度 ViReMa (满配大一统版本)")
    
    parser.add_argument("-s", "--summary", required=True, help="包含 [Sample, Accession] 关系表格")
    parser.add_argument("-r", "--reference", required=True, help="病毒总群 FASTA 指导序列")
    parser.add_argument("-d", "--reads_dir", required=True, help="脱宿主后的原始队列仓")
    parser.add_argument("-v", "--virema_script", required=True, help="ViReMa.py 超级指向器定位")
    parser.add_argument("-o", "--output_dir", default="./virema_out", help="总体产出宇宙")
    
    vr = parser.add_argument_group("核心逻辑阀值与脱水")
    vr.add_argument("--min_cov", type=float, default=0.0, help="【脱水】设定抛弃阈值，如 80 即切除 Coverage 低于 80% 的命中 (默认0,即全通)")
    vr.add_argument("--seed", type=int, default=25, help="ViReMa 切割重组的初始比对种子长")
    vr.add_argument("--mindel", type=int, default=15, help="划定为微缺失(假阳)的最高值界限")
    vr.add_argument("--defuzz", type=int, choices=[0, 3, 5], default=0, help="对齐模糊碱基断点处理(0为中立)")
    
    sys_group = parser.add_argument_group("天火并发阵列调配")
    sys_group.add_argument("-j", "--jobs", type=int, default=4, help="外层全域样本同时发动数 (-j)")
    sys_group.add_argument("-t", "--threads", type=int, default=8, help="内层计算裂变所用核心算力")
    sys_group.add_argument("--resume", action="store_true", help="接防中断的运算防线 (跨越已生成的SAM)")
    sys_group.add_argument("--keep", action="store_true", help="强行扣留原本要被沙盒炸毁的废料现场做人工排查")
    sys_group.add_argument("--shm", action="store_true", help="🔥 [狂飙核心] 开启 /dev/shm 高速内存盘承载 ViReMa 迭代引擎！")
    
    args = parser.parse_args()
    BatchViReMaPipeline(args).run()
