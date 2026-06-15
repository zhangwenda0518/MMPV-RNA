#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Virome Pipeline v2.3 — 宏病毒组端到端全自动主控流水线
===============================================================================

数据流:
  Raw FASTQ ──→ [00a_CleanData] ──→ [00b_HostDepletion] ──→ [01_Assembly]
                                                                   │
   [05_Reports] ←── [04_CLUSTER] ←── [03_COBRA] ←── [02_Identification] ←──┘

依赖脚本 (同目录):
  clean-data.py              → Step 0a: Fastp + Seqkit + Clumpify
  host_depletion.py          → Step 0b: Kraken2 + Align + Ribodetector
  assembly_pipeline.py       → Step 1:  MEGAHIT / rnaviralSPAdes / Penguin
  virus_identification16.py  → Step 2:  Genomad + Blast + VirSorter2 + ...
  cobra_pipeline.py          → Step 3:  BWA-MEM2 + COBRA 批量延伸
  cluster_pipeline.py        → Step 4:  CLUSTER 三支路病毒基因组去冗余
  virus_classifier2.py      → Step 5:  病毒分类注释 (直接调用)
  run_host_prediction.py  → Step 6:  宿主预测

所有 CLI 参数精确匹配底层脚本的真实参数名。
===============================================================================
"""

import os
import sys
import argparse
import subprocess
import logging
import re
import shutil
from pathlib import Path
from collections import defaultdict

import polars as pl
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent


# ═══════════════════════════════════════════════════════════════════
# 1. 日志与工具函数
# ═══════════════════════════════════════════════════════════════════

def setup_logger(output_dir, level='INFO'):
    """配置双通道日志: 控制台 INFO + 文件 DEBUG"""
    logger = logging.getLogger("ViromeOrch")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(output_dir, exist_ok=True)

    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, level.upper(), logging.INFO))
    ch.setFormatter(logging.Formatter(
        '[%(asctime)s] %(levelname)s %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)

    fh = logging.FileHandler(
        os.path.join(output_dir, 'orchestrator.log'), encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '[%(asctime)s] %(name)s %(levelname)s %(message)s'))
    logger.addHandler(fh)

    return logger


def find_cmd(name):
    """查找命令: 优先 PATH, 回退同名"""
    found = shutil.which(name)
    return found if found else name


def run_cmd(cmd, logger, step_name, log_file=None):
    """执行 shell 命令。stdout/stderr 同时输出到终端、主编排日志、可选的阶段日志。"""
    logger.info("[%s] 执行中...", step_name)
    logger.debug("  CMD: %s", cmd)
    lf = open(log_file, 'a', encoding='utf-8') if log_file else None
    try:
        if lf:
            lf.write(f"=== {step_name} ===\nCMD: {cmd}\n\n")
        with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, bufsize=1) as proc:
            for line in proc.stdout:
                line = line.rstrip('\n\r')
                print(line, flush=True)
                logger.debug("  %s", line)
                if lf:
                    lf.write(line + '\n')
                    lf.flush()
            proc.wait()
            if proc.returncode != 0:
                logger.error("[%s] ✗ 失败 (exit=%d)", step_name, proc.returncode)
                return False, f"exit={proc.returncode}"
        logger.info("[%s] ✓ 成功", step_name)
        return True, ""
    except Exception as e:
        logger.error("[%s] ✗ 异常: %s", step_name, e)
        return False, str(e)
    finally:
        if lf:
            lf.write(f"\n=== {step_name} 完成 ===\n\n")
            lf.close()
            # 复制到 09_Reports/logs/ (从 output_dir/orchestrator.log 反推)
            try:
                import shutil
                for h in logger.handlers:
                    if isinstance(h, logging.FileHandler):
                        reports_log = Path(h.baseFilename).resolve().parent / "09_Reports" / "logs"
                        reports_log.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(log_file, str(reports_log / Path(log_file).name))
                        break
            except: pass


# ═══════════════════════════════════════════════════════════════════
# 2. 文件扫描与样本追踪
# ═══════════════════════════════════════════════════════════════════

SEQ_EXTS = ['.fastq', '.fq', '.fasta', '.fa',
            '.fastq.gz', '.fq.gz', '.fasta.gz', '.fa.gz',
            '_fastq.gz', '_fq.gz', '_fasta.gz', '_fa.gz']


def find_seq_files(search_dir):
    """递归查找序列文件，返回去重排序列表"""
    search_dir = Path(search_dir)
    files = []
    for ext in SEQ_EXTS:
        files.extend(search_dir.rglob(f"*{ext}"))
    seen = set()
    return sorted([f for f in files if f.is_file() and not (str(f) in seen or seen.add(str(f)))])


def extract_base_sample(filename):
    """
    提取基础样本名。依次剥离:
      扩展名(.fastq.gz) → R1/R2/_1/_2 标签 → _clean → 组装/鉴定后缀
    """
    name = os.path.basename(str(filename))
    for pat in [
        r'[._]f(ast)?[aq](\.(gz|bz2))?$',
        r'[._-][Rr][12](?=[._-]|\.|$)',
        r'(?<=[._-])[12](?=[._-]|\.|$)',
        r'_S\d+_L\d+',    # Illumina sample/lane tag
        r'_clean',
        r'_megahit\.contig',
        r'_rnaviralspades\.contig',
        r'_penguin\.contig',
        r'_all_tools_refineC_merge',
        r'\.contig',
        r'\.merged',
    ]:
        name = re.sub(pat, '', name, flags=re.IGNORECASE)
    return name.strip('_.-')


def scan_samples_in_dir(d):
    """扫描目录, 返回 {base_sample: {'r1':Path, 'r2':Path|None}}"""
    files = find_seq_files(d)
    if not files:
        return {}
    by_base = defaultdict(list)
    for f in files:
        by_base[extract_base_sample(f)].append(f)

    samples = {}
    for base, flist in by_base.items():
        flist.sort(key=str)
        r1 = r2 = None
        for f in flist:
            n = os.path.basename(str(f))
            if re.search(r'[._-][Rr]1[._-]|[._-]1[._-]', n):
                r1 = f
            elif re.search(r'[._-][Rr]2[._-]|[._-]2[._-]', n):
                r2 = f
        if len(flist) == 2 and (not r1 or not r2):
            r1, r2 = flist[0], flist[1]
        elif len(flist) == 1 and not r1:
            r1 = flist[0]
        if r1:
            samples[base] = {'r1': r1, 'r2': r2}
    return samples


def scan_contig_files(asm_dir, tools):
    """扫描 assembly 输出: {sample_name: {tool: contig_path}}"""
    result = defaultdict(dict)
    for d in Path(asm_dir).iterdir():
        if not d.is_dir():
            continue
        for tool in tools:
            cf = d / f"{d.name}_{tool}.contig.fasta"
            if cf.exists() and cf.stat().st_size > 0:
                result[d.name][tool] = cf
    return dict(result)


def scan_viral_files(ident_dir):
    """扫描 virus_identification 输出: {sample_name: viral_fasta_path}"""
    result = {}
    for d in Path(ident_dir).iterdir():
        if not d.is_dir():
            continue
        for pat in ['*_virus.all.candidate.fasta', '*_virus.candidate.fasta', '*.virus.candidate.fasta']:
            for vf in d.glob(pat):
                if vf.stat().st_size > 50:
                    result[d.name] = vf
                    break
            if d.name in result:
                break
    return result


def fuzzy_match(target, candidates):
    """模糊匹配: 精确 → 包含 → 去 _clean 后匹配"""
    if target in candidates:
        return target
    for c in candidates:
        if target in c or c in target:
            return c
    tc = target.replace('_clean', '')
    for c in candidates:
        if c.replace('_clean', '') == tc:
            return c
    return None


# ═══════════════════════════════════════════════════════════════════
# 3. 核心流水线
# ═══════════════════════════════════════════════════════════════════

class ViromePipeline:
    def __init__(self, args, logger):
        self.args = args
        self.log = logger
        out = Path(args.output_dir).absolute()
        raw = Path(args.input_reads).absolute() if args.input_reads else Path(args.output_dir).absolute()
        script_dir = Path(__file__).parent.resolve()

        # 标准化目录
        self.d = {
            'raw':       raw,
            'root':      out,
            'clean':     out / '00a_CleanData',
            'hostdep':   out / '00b_HostDepletion',
            'asm':       out / '01_Assembly',
            'ident':     out / '02_Identification',
            'cobra':     out / '03_COBRA',
            'cluster':   out / '04_CLUSTER',
            'taxonomy':  out / '05_Taxonomy',
            'host_pred':  out / '06_HostPrediction',
            'checkv_dir': out / '07_Checkv',
            'rescue_dir': out / '08_Rescue',
            'reports':    out / '09_Reports',
        }
        # 仅创建根目录, 各阶段按需创建自己的子目录
        self.d['root'].mkdir(parents=True, exist_ok=True)

        # 脚本路径
        self.sc = {
            'clean':    script_dir / 'clean-data.py',
            'deplete':  script_dir / 'host_depletion.py',
            'assembly': script_dir / 'assembly_pipeline.py',
            'identify': script_dir / 'virus_identification16.py',
            'cobra':    script_dir / 'cobra_pipeline.py',
            'cluster':  script_dir / 'cluster_pipeline.py',
            'rescue':   script_dir / 'rescue_pipeline.py',
            'host_pred': script_dir / 'run_host_prediction.py',
            'classifier': script_dir / 'virus_classifier2.py',
            'classifier_R': script_dir / 'virus_classifier_analysis14.R',
            'c9': script_dir / 'C9_classify_contigs.py',
        }

        # 数据流指针
        self.reads_dir = self.d['raw']

        # 验证
        self._validate()

        # 检测原始样本
        no_reads_stages = {'identification', 'cluster', 'taxonomy', 'host', 'checkv', 'report'}
        if not set(self.args.stage).issubset(no_reads_stages):
            self.orig_samples = scan_samples_in_dir(self.reads_dir)
            if not self.orig_samples:
                logger.error("在 %s 中未找到序列文件!", self.reads_dir)
                sys.exit(1)
            pe = sum(1 for v in self.orig_samples.values() if v['r2'])
            se = len(self.orig_samples) - pe
            logger.info("检测到 %d 个样本 (PE=%d, SE=%d)", len(self.orig_samples), pe, se)
        else:
            self.orig_samples = {}

    def _validate(self):
        for name, path in self.sc.items():
            if not path.exists():
                self.log.error("致命: 找不到 %s", path)
                sys.exit(1)

        stages = set(self.args.stage)
        _all = 'all' in stages
        doing_clean = _all or 'clean' in stages
        need_virus_db = _all or bool(stages & {'identification', 'taxonomy', 'host'})
        need_checkv_db = _all or bool(stages & {'checkv', 'rescue'})

        # deplete 阶段需要 kraken2 + host_align (或 --host_db 自动检测)
        need_deplete_db = _all or 'deplete' in stages or ('clean' in stages and not self.args.skip_depletion)
        if need_deplete_db and not self.args.host_db:
            if not self.args.kraken2_db or not self.args.host_align_db:
                self.log.error("致命: deplete 阶段需要 --kraken2_db/--host_align_db 或 --host_db")
                sys.exit(1)

        # assembly/identification/cobra/cluster/taxonomy/host 阶段需要 virus_db
        if need_virus_db:
            if not self.args.virus_db:
                self.log.error("致命: 需要 --virus_db")
                sys.exit(1)

        # cobra/rescue 阶段需要 checkv_db
        if need_checkv_db:
            if not self.args.checkv_db:
                self.log.error("致命: 需要 --checkv_db")
                sys.exit(1)

        for exe in ['python']:
            if shutil.which(exe) is None:
                self.log.error("致命: 找不到 %s", exe)
                sys.exit(1)
        self.log.info("环境验证通过 (stage=%s)", ','.join(sorted(stages)))

    # ── Step 0a: 数据清洗 ──
    def run_clean(self):
        self.d['clean'].mkdir(parents=True, exist_ok=True)
        if self.args.skip_clean:
            self.log.info("[0a] 跳过 (--skip_clean)")
            return

        self.log.info("=" * 50)
        self.log.info("[0a] Fastp → Seqkit → Clumpify")

        # clean-data.py 完整参数: --input, --output, --fastp-threads, --jobs,
        #   --skip-clumpify, --force, --dedup, --clumpify-memory, --no-compress, --debug
        parts = [
            f"python {self.sc['clean']}",
            f"--input {self.reads_dir}",
            f"--output {self.d['clean']}",
            f"--fastp-threads {self.args.threads}",
            f"--jobs {self.args.jobs}",
        ]
        if self.args.skip_clumpify:
            parts.append("--skip-clumpify")
        if self.args.force:
            parts.append("--force")
        if self.args.dedup:
            parts.append("--dedup")
        if self.args.clumpify_memory:
            parts.append(f"--clumpify-memory {self.args.clumpify_memory}")
        if self.args.no_compress:
            parts.append("--no-compress")
        if self.args.clean_debug:
            parts.append("--debug")

        ok, _ = run_cmd(' '.join(parts), self.log, "Clean", str(self.d['clean'] / "clean.log"))
        if not ok:
            self.log.error("清洗失败, 终止。")
            sys.exit(1)

        # 更新指针: 优先 clumpify 否则 fasta
        cl = self.d['clean'] / '3.clumpify'
        fa = self.d['clean'] / '2.fasta'
        self.reads_dir = cl if (cl.exists() and any(cl.iterdir())) else fa
        self.log.info("  Reads → %s", self.reads_dir)

    # ── Step 0b: 去宿主 ──
    def run_depletion(self):
        self.d['hostdep'].mkdir(parents=True, exist_ok=True)
        if self.args.skip_depletion:
            self.log.info("[0b] 跳过 (--skip_depletion)")
            return

        # 自动从 --host_db 推导子数据库路径
        kraken2_db = self.args.kraken2_db
        host_align_db = self.args.host_align_db
        if self.args.host_db:
            hdb = Path(self.args.host_db)
            if not kraken2_db:
                for sub in ['kraken2', 'kraken2_db', 'kraken']:
                    if (hdb / sub).is_dir():
                        kraken2_db = str(hdb / sub)
                        break
            if not host_align_db:
                aligner = self.args.aligner
                for sub in [aligner, f'{aligner}_index', 'align']:
                    if (hdb / sub).is_dir():
                        # 找索引前缀 (bowtie2: host.1.bt2, hisat2: host.1.ht2, minimap2: host_*.mmi)
                        for prefix in ['host', 'index', 'genome']:
                            test = hdb / sub / prefix
                            if aligner == 'bowtie2' and ((test.parent / (prefix + '.1.bt2')).is_file() or (test.parent / (prefix + '.1.bt2l')).is_file()):
                                host_align_db = str(test); break
                            if aligner == 'hisat2' and (test.parent / (prefix + '.1.ht2')).is_file():
                                host_align_db = str(test); break
                            if aligner == 'minimap2':
                                mmi = test.parent / f'{prefix}_{self.args.seq_type}.mmi'
                                if mmi.is_file():
                                    host_align_db = str(mmi); break
                        if host_align_db:
                            break
            self.log.info("  host_db: %s → kraken2=%s, align=%s", hdb, kraken2_db or '?', host_align_db or '?')

        if not kraken2_db:
            self.log.error("致命: 需要 --kraken2_db 或 --host_db")
            sys.exit(1)
        if not host_align_db:
            self.log.error("致命: 需要 --host_align_db 或 --host_db")
            sys.exit(1)

        self.log.info("=" * 50)
        self.log.info("[0b] Kraken2 → Align → Ribodetector")

        # host_depletion.py 完整参数:
        #   --tool, --seq-type, --kraken2_index, --step2_index,
        #   --input-dir, --outdir, --jobs, --threads, --logs_dir,
        #   --rrna, --rrna_tool, --silva_index, --filter, --confidence,
        #   --tmp, --force, --keep_rrna, --chunk_size, --rrna_report,
        #   --steps, --config, --debug
        parts = [
            f"python {self.sc['deplete']}",
            f"--tool {self.args.aligner}",
            f"--seq-type {self.args.seq_type}",
            f"--kraken2_index {kraken2_db}",
            f"--step2_index {host_align_db}",
            f"--input-dir {self.reads_dir}",
            f"--outdir {self.d['hostdep']}",
            f"--jobs {self.args.jobs}",
            f"--threads {self.args.threads}",
            f"--logs_dir {self.d['hostdep']}/logs",
            "--filter true",
        ]
        if self.args.rrna:
            parts.append("--rrna")
            parts.append(f"--rrna_tool {self.args.rrna_tool}")
            if self.args.silva_index:
                parts.append(f"--silva_index {self.args.silva_index}")
        if self.args.force:
            parts.append("--force")
        if self.args.deplete_tmp:
            parts.append(f"--tmp {self.args.deplete_tmp}")
        parts.append(f"--confidence {self.args.kraken2_confidence}")
        if self.args.keep_rrna:
            parts.append("--keep_rrna")
        parts.append(f"--chunk_size {self.args.rrna_chunk_size}")
        parts.append(f"--rrna_report {self.args.rrna_report}")
        parts.append(f"--steps {self.args.deplete_steps}")
        if self.args.align_config:
            parts.append(f"--config {self.args.align_config}")
        if self.args.deplete_debug:
            parts.append("--debug")

        ok, _ = run_cmd(' '.join(parts), self.log, "Host Depletion", str(self.d['hostdep'] / "hostdep.log"))
        if not ok:
            self.log.error("去宿主失败, 终止。")
            sys.exit(1)

        self.reads_dir = self.d['hostdep']
        self.log.info("  Reads → %s", self.reads_dir)

    # ── Step 1: 组装 ──
    def run_assembly(self):
        self.d['asm'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[1] MEGAHIT / rnaviralSPAdes / Penguin")

        # assembly_pipeline.py 完整参数:
        #   --tool, --input, --length, --threads, --memory, --jobs,
        #   --output-dir, --log_dirs, --refineC_split, --refineC_merge,
        #   --refineC_threads, --refineC_frag_min_len, --refineC_min_id, --refineC_min_cov,
        #   --tmp-dir, --keep-temp, --force
        parts = [
            f"python {self.sc['assembly']}",
            f"--tool {self.args.assembler}",
            f"--input {self.reads_dir}",
            f"--output-dir {self.d['asm']}",
            f"--threads {self.args.threads}",
            f"--length {self.args.contig_length}",
            f"--memory {self.args.memory}",
            f"--jobs {self.args.jobs}",
            f"--log_dirs {self.d['asm']}/logs",
        ]
        # 多工具组装时自动启用 refineC split + merge
        asm_tools = self.args.assembler.split(",") if self.args.assembler != 'all' else ['megahit', 'rnaviralspades', 'penguin']
        if len(asm_tools) >= 2:
            parts.append("--refineC_split --refineC_merge")
        if self.args.force:
            parts.append("--force")
        if self.args.refinec_threads is not None:
            parts.append(f"--refineC_threads {self.args.refinec_threads}")
        parts.append(f"--refineC_frag_min_len {self.args.refinec_frag_min_len}")
        parts.append(f"--refineC_min_id {self.args.refinec_min_id}")
        parts.append(f"--refineC_min_cov {self.args.refinec_min_cov}")
        if self.args.asm_tmp_dir:
            parts.append(f"--tmp-dir {self.args.asm_tmp_dir}")
        if self.args.asm_keep_temp:
            parts.append("--keep-temp")

        ok, _ = run_cmd(' '.join(parts), self.log, "Assembly", str(self.d['asm'] / "assembly.log"))
        if not ok:
            self.log.error("组装失败, 终止。")
            sys.exit(1)

        self.asm_map = scan_contig_files(self.d['asm'], asm_tools)
        self.log.info("  组装完成: %d 样本有 contig 输出", len(self.asm_map))

    # ── Step 2: 病毒鉴定 ──
    def run_identification(self):
        self.d['ident'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[2] Genomad / Blast / VirSorter2 / ViralVerify / VirHunter / Metabuli")

        asm_map = getattr(self, 'asm_map', {})
        if not asm_map:
            # 独立运行时: 优先 --input_assembly, 否则默认 01_Assembly/
            asm_dir = Path(self.args.input_assembly) if self.args.input_assembly else self.d['asm']
            scan_tools = ['megahit', 'rnaviralspades', 'penguin']
            asm_map = scan_contig_files(asm_dir, scan_tools)
        if not asm_map:
            self.log.error("无 contig 文件, 跳过鉴定。")
            return

        # 直接传目录, 子脚本内部 --jobs 并行处理所有样本
        self.log.info("  %d 个样本待鉴定", len(asm_map))

        # virus_identification16.py 完整参数:
        #   --input, --output, --db_dir, --identify_tools, --threads, --jobs,
        #   --blast_mode, --blast_evalue, --blast_top_n, --virsorter_group,
        #   --virus_protein_db, --uniprot_db, --viroids_db, --virsorter_db,
        #   --viralverify_hmm, --metabuli_db, --virus_taxid,
        #   --virhunter_path, --virhunter_weights, --virbot_path, --viralm_path,
        #   --nr_db, --skip_uniprot_filter, --skip_nr_filter,
        #   --skip_plots, --clean_failed, --extension, --force
        parts = [
            f"python {self.sc['identify']}",
            f"--input {asm_dir}",
            f"--output {self.d['ident']}",
            f"--db_dir {self.args.virus_db}",
            f"--identify_tools {self.args.identify_tools}",
            f"--threads {self.args.threads}",
            f"--jobs {self.args.jobs}",
            f"--blast_mode {self.args.blast_mode}",
            f"--blast_evalue {self.args.blast_evalue}",
            f"--blast_top_n {self.args.blast_top_n}",
            f"--virsorter_group {self.args.virsorter_group}",
            f"--extension {self.args.ident_ext}",
        ]
        for arg in ['virus_protein_db', 'uniprot_db', 'viroids_db', 'virsorter_db',
                     'viralverify_hmm', 'metabuli_db', 'virus_taxid',
                     'virhunter_path', 'virhunter_weights', 'virbot_path', 'viralm_path',
                     'nr_db']:
            val = getattr(self.args, arg, None)
            if val:
                parts.append(f"--{arg} {val}")
        if self.args.force:
            parts.append("--force")
        if self.args.skip_uniprot_filter:
            parts.append("--skip_uniprot_filter")
        if self.args.skip_nr_filter:
            parts.append("--skip_nr_filter")
        if self.args.skip_id_plots:
            parts.append("--skip_plots")
        if self.args.clean_failed:
            parts.append("--clean_failed")
        run_cmd(' '.join(parts), self.log, "VirusIdentification", str(self.d['ident'] / "ident.log"))

        self.viral_map = scan_viral_files(self.d['ident'])
        self.log.info("  鉴定完成: %d 样本有病毒候选序列", len(self.viral_map))

    # ── Step 3: COBRA 延伸 ──
    def run_cobra(self):
        self.d['cobra'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[3] COBRA 批量延伸 (BWA-MEM2 + COBRA + CheckV)")

        # cobra_pipeline.py 完整参数:
        #   --mode, --reads-dir, --contigs-dir, --virsorter-dir, --output-dir,
        #   --assembly-tools, --virus-mode, --jobs, --threads,
        #   --mink, --maxk, --linkage-mismatch,
        #   --resume/--no-resume, --verbose
        # 自动检测哪些工具实际有组装输出
        all_tools = ['megahit', 'rnaviralspades', 'penguin']
        asm_tools = []
        for tool in all_tools:
            for d in self.d['asm'].iterdir():
                if d.is_dir() and (d / f"{d.name}_{tool}.contig.fasta").exists():
                    asm_tools.append(tool)
                    break
        if not asm_tools:
            asm_tools = [self.args.assembler] if self.args.assembler != 'all' else all_tools
        self.log.info("  Auto-detect 组装工具: %s", ','.join(asm_tools))

        parts = [
            f"python {self.sc['cobra']}",
            f"--mode mix",
            f"--reads-dir {self.reads_dir}",
            f"--contigs-dir {self.d['asm']}",
            f"--virsorter-dir {self.d['ident']}",
            f"--output-dir {self.d['cobra']}",
            f"--assembly-tools {','.join(asm_tools)}",
            f"--virus-mode {self.args.virus_mode}",
            f"--jobs {max(1, self.args.jobs // 2)}",
            f"--threads {self.args.threads}",
            f"--mink {self.args.cobra_mink}",
            f"--maxk {self.args.cobra_maxk}",
            f"--linkage-mismatch {self.args.cobra_linkage_mismatch}",
        ]
        if self.args.force:
            parts.append("--no-resume")
        if self.args.cobra_verbose:
            parts.append("--verbose")

        ok, _ = run_cmd(' '.join(parts), self.log, "COBRA Pipeline", str(self.d['cobra'] / "cobra.log"))
        if not ok:
            self.log.warning("COBRA 阶段部分任务失败, 检查日志。")

        self.log.info("  COBRA 阶段完成")
        self.log.info("  输出: %s", self.d['cobra'])

    # ── Step 4: CLUSTER 三支路去冗余 ──
    def run_cluster(self):
        self.d['cluster'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[4] CLUSTER 三支路病毒基因组去冗余")

        # 优先使用直接输入, 否则按样本收集 03_COBRA: 有 cobra.fa 用 cobra, 否则 virus.fa
        if self.args.cluster_input and os.path.isfile(self.args.cluster_input):
            cluster_fa = Path(self.args.cluster_input)
            self.log.info("  CLUSTER 直接输入: %s", cluster_fa)
        else:
            cobra_dir = self.d['cobra']
            if not cobra_dir.is_dir():
                self.log.error("COBRA 输出 %s 不存在, 跳过 CLUSTER", cobra_dir)
                return
            cluster_fa = Path(self.d['root']) / "cluster_input.fasta"

            n_cobra, n_virus = 0, 0
            with open(cluster_fa, 'w') as out:
                for sd in cobra_dir.iterdir():
                    if not sd.is_dir():
                        continue
                    # 检查该样本是否有 cobra 延伸结果
                    cobra_files = list(sd.rglob('*.cobra.fa'))
                    if cobra_files:
                        for cf in cobra_files:
                            with open(cf) as inf:
                                out.write(inf.read())
                                n_cobra += 1
                    else:
                        # 无延伸, 回退到原始 virus.fa
                        for vf in sd.rglob('*virus*.fasta'):
                            with open(vf) as inf:
                                out.write(inf.read())
                                n_virus += 1
                        for vf in sd.rglob('*virus*.fa'):
                            with open(vf) as inf:
                                out.write(inf.read())
                                n_virus += 1

            if cluster_fa.stat().st_size == 0:
                self.log.warning("未找到任何输入, 跳过 CLUSTER")
                return
            if n_virus:
                self.log.info("  CLUSTER 输入 (自动收集): %d cobra + %d virus → %s", n_cobra, n_virus, cluster_fa)
            else:
                self.log.info("  CLUSTER 输入 (自动收集): %d cobra → %s", n_cobra, cluster_fa)

        self.log.info("  CLUSTER 输入: %s", cluster_fa)

        # cluster_pipeline.py 完整参数:
        #   -i, -o, -t, --min-length, --ani, --qcov,
        #   --ref-genomes, --cdhit-ani, --cdhit-qcov,
        #   --skip-vclust, --vclust-cluster-file, --resume
        parts = [
            f"python {self.sc['cluster']}",
            f"-i {cluster_fa}",
            f"-o {self.d['cluster']}",
            f"-t {self.args.threads}",
            f"--min-length {self.args.min_length}",
            f"--ani {self.args.ani}",
            f"--qcov {self.args.qcov}",
        ]
        if self.args.ref_genomes:
            parts.append(f"--ref-genomes {' '.join(self.args.ref_genomes)}")
        if self.args.cdhit_ani:
            parts.append(f"--cdhit-ani {self.args.cdhit_ani}")
        if self.args.cdhit_qcov:
            parts.append(f"--cdhit-qcov {self.args.cdhit_qcov}")
        if self.args.skip_vclust:
            parts.append("--skip-vclust")
        if self.args.vclust_cluster_file:
            parts.append(f"--vclust-cluster-file {self.args.vclust_cluster_file}")
        if not self.args.force:
            parts.append("--resume")

        ok, _ = run_cmd(' '.join(parts), self.log, "CLUSTER (vclust only)", str(self.d['cluster'] / "cluster.log"))
        if not ok:
            self.log.error("CLUSTER vclust 阶段失败, 终止。")
            sys.exit(1)

        centroids = Path(self.d['cluster']) / "centroids" / "final_centroids.fasta"
        if centroids.exists():
            n = sum(1 for _ in open(centroids) if _.startswith('>'))
            self.log.info("  CLUSTER 输出: %d 条 centroids → %s", n, centroids)
        else:
            self.log.error("  CLUSTER 未产出 centroids!")
            sys.exit(1)

        self.log.info("  CLUSTER (vclust) 阶段完成 — 三支路拯救交由 rescue 阶段")

    # ── Rescue: 按宿主过滤 → 三支路级联拯救 ──
    def run_rescue(self):
        self.d['rescue_dir'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[Rescue] 宿主过滤 + 三支路级联拯救")

        host_summary = self.d['host_pred'] / "ensemble_host_summary.tsv"
        if not host_summary.exists():
            self.log.error("宿主预测结果 %s 不存在, 请先运行 --stage host", host_summary)
            sys.exit(1)

        centroids_fa = self.d['cluster'] / "centroids" / "final_centroids.fasta"
        if not centroids_fa.exists():
            self.log.error("centroids %s 不存在, 请先运行 --stage cluster", centroids_fa)
            sys.exit(1)

        clusters_tsv = self.d['cluster'] / "3_vclust" / "vclust_clusters.tsv"
        if not clusters_tsv.exists():
            self.log.error("clusters.tsv %s 不存在", clusters_tsv)
            sys.exit(1)

        split_dir = self.d['cluster'] / "3_vclust" / "split_fastas"
        if not split_dir.exists():
            self.log.error("split_fastas %s 不存在, 请重新运行 --stage cluster", split_dir)
            sys.exit(1)

        # 1. 加载宿主预测
        host_df = pl.read_csv(str(host_summary), separator="\t", null_values=["NA", "N/A", ""])
        if "Final_Host" not in host_df.columns or "contig_id" not in host_df.columns:
            self.log.error("宿主 TSV 缺少必需列 (contig_id, Final_Host)")
            sys.exit(1)

        host_filter = self.args.host_filter or ["Plant"]
        if isinstance(host_filter, str):
            host_filter = [h.strip() for h in host_filter.split(",")]

        self.log.info("  目标宿主: %s", ", ".join(host_filter))
        self.log.info("  宿主分布:")
        for h, n in host_df.group_by("Final_Host").agg(pl.len()).sort("len", descending=True).iter_rows():
            self.log.info("    %s: %d", h if h else "Unknown", n)

        # 2. 分离: 目标宿主 / Unknown / 其他
        target_ids = set()
        unknown_ids = set()
        other_ids = {}

        for row in host_df.iter_rows(named=True):
            cid = row["contig_id"]
            fh = row["Final_Host"]
            if fh and fh != "NA" and fh != "Unknown":
                if fh in host_filter:
                    target_ids.add(cid)
                else:
                    other_ids.setdefault(fh, []).append(cid)
            else:
                unknown_ids.add(cid)

        self.log.info("  目标宿主 (%s): %d 条 centroids", ",".join(host_filter), len(target_ids))
        self.log.info("  Unknown:         %d 条 centroids", len(unknown_ids))
        self.log.info("  其他宿主:        %d 条 centroids", sum(len(v) for v in other_ids.values()))

        # 预先加载 centroids 序列 (避免循环内重复解析)
        centroids_map = {}
        for rec in SeqIO.parse(str(centroids_fa), "fasta"):
            centroids_map[rec.id] = rec

        # 3. 输出 Unknown centroids
        if unknown_ids:
            unknown_out = self.d['cluster'] / "centroids" / "unknown_votus.fasta"
            with open(unknown_out, "w") as uf:
                written = 0
                for cid in unknown_ids:
                    if cid in centroids_map:
                        SeqIO.write(centroids_map[cid], uf, "fasta")
                        written += 1
            self.log.info("  Unknown → %s (%d 条)", unknown_out, written)

        # 4. 输出其他宿主 (按宿主分文件)
        for host, ids in other_ids.items():
            if not ids:
                continue
            safe_host = host.replace("/", "_").replace(" ", "_")
            host_out = self.d['cluster'] / "centroids" / f"skipped_{safe_host}.fasta"
            with open(host_out, "w") as hf:
                written = 0
                for cid in ids:
                    if cid in centroids_map:
                        SeqIO.write(centroids_map[cid], hf, "fasta")
                        written += 1
            if written > 0:
                self.log.info("  %s → %s (%d 条)", host, host_out, written)

        if not target_ids:
            self.log.warning("  无目标宿主 centroids, 跳过三支路拯救")
            return

        # 4.5 区分 CD-HIT known vs vclust novel
        known_id_file = self.d['cluster'] / "centroids" / "known_ids.txt"
        cdhit_known_ids = set()
        if known_id_file.is_file():
            with open(known_id_file) as kf:
                for line in kf:
                    cdhit_known_ids.add(line.strip())

        # 读取 CheckV pass IDs (≥90% completeness, 也免拯救)
        checkv_pass_file = self.d['checkv_dir'] / "checkv_pass_ids.txt"
        checkv_pass_ids = set()
        if checkv_pass_file.is_file():
            with open(checkv_pass_file) as pf:
                for line in pf:
                    checkv_pass_ids.add(line.strip())

        # 免拯救 = CD-HIT known + CheckV pass (≥90%)
        skip_rescue_ids = cdhit_known_ids | checkv_pass_ids
        target_known = target_ids & skip_rescue_ids
        target_novel = target_ids - skip_rescue_ids
        n_cdhit = len(target_ids & cdhit_known_ids)
        n_checkv = len(target_ids & checkv_pass_ids - cdhit_known_ids)
        if skip_rescue_ids:
            self.log.info("  免拯救: %d CD-HIT + %d CheckV(≥90%%) = %d 条", n_cdhit, n_checkv, len(target_known))
        self.log.info("  vclust novel: %d 在目标宿主中 (进入三支路拯救)", len(target_novel))

        # 4.6 CD-HIT known + CheckV pass centroids → 直接输出
        if target_known:
            known_out_dir = self.d['rescue_dir'] / "known"; known_out_dir.mkdir(parents=True, exist_ok=True)
            known_final = known_out_dir / "centroids"; known_final.mkdir(parents=True, exist_ok=True)
            known_centroids_fa = known_final / "final_centroids.fasta"

            known_seqs = []
            centroids_map2 = {}
            for rec in SeqIO.parse(str(centroids_fa), "fasta"):
                centroids_map2[rec.id] = rec
            for cid in target_known:
                if cid in centroids_map2:
                    known_seqs.append(centroids_map2[cid])

            if known_seqs:
                SeqIO.write(known_seqs, str(known_centroids_fa), "fasta")
                self.log.info("  CD-HIT known → %s (%d 条, 完整参考基因组)", known_centroids_fa, len(known_seqs))

        if not target_novel:
            self.log.warning("  无 vclust novel centroids, 跳过三支路拯救")
            return

        # 5. 加载 clusters.tsv → 找到目标 novel centroids 所在的 cluster 成员
        target_clusters = set()
        with open(clusters_tsv) as f:
            f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    member_id, cluster_id = parts[0].strip(), parts[1].strip()
                    if member_id in target_novel:
                        target_clusters.add(cluster_id)

        # CD-HIT known 簇 (有 target_known centroids 的簇)
        cdhit_cluster_files = []
        if target_known:
            for fa in split_dir.glob("cdhit_cluster_*.all.fasta"):
                cdhit_cluster_files.append(fa)

        self.log.info("  涉及 cluster: %d vclust + %d cdhit-known", len(target_clusters), len(cdhit_cluster_files))

        # 6. 写入 target_novel centroids FASTA (供 rescue_pipeline.py 使用)
        rescue_centroids = self.d['cluster'] / "rescue_centroids.fasta"
        n_written = 0
        with open(rescue_centroids, "w") as rf:
            for cid in target_novel:
                if cid in centroids_map:
                    SeqIO.write(centroids_map[cid], rf, "fasta")
                    n_written += 1
        self.log.info("  目标 centroids: %d 条 → %s", n_written, rescue_centroids)

        if n_written == 0:
            self.log.error("  无有效 centroids, 跳过三支路拯救")
            return

        # 7. 调用 rescue_pipeline.py (直接使用已有聚类, 不重新 vclust)
        rescue_out = self.d['rescue_dir'] / f"{'_'.join(host_filter)}"
        rescue_out.mkdir(parents=True, exist_ok=True)

        parts = [
            f"python {self.sc['rescue']}",
            f"-c {rescue_centroids}",
            f"--clusters-tsv {clusters_tsv}",
            f"--split-dir {split_dir}",
            f"-o {rescue_out}",
            f"-fq {self.reads_dir}",
            f"-cv {self.args.checkv_db}",
            f"-t {self.args.threads}",
            f"-j {self.args.jobs}",
            f"--ani {self.args.ani}",
            f"--qcov {self.args.qcov}",
        ]
        if self.args.blast_db:
            parts.append(f"-db {self.args.blast_db}")
        vsi_path = self.args.virseqimprover_path or str(self.sc['cluster'].parent / 'Virseqimprover.py')
        parts += [f"--virseqimprover-path {vsi_path}"]
        parts += [f"--salmon-bin {self.args.salmon_bin}"]
        parts += [f"--max-vsi-samples {self.args.max_vsi_samples}"]
        parts += [f"--min-vsi-len {self.args.min_vsi_len}"]
        if hasattr(self.args, 'checkv_threshold'):
            parts += [f"--checkv-threshold {self.args.checkv_threshold}"]
        if not self.args.force:
            parts.append("--resume")

        ok, _ = run_cmd(' '.join(parts), self.log, f"Rescue ({','.join(host_filter)})", str(rescue_out / "rescue.log"))
        if ok:
            final_out = rescue_out / "centroids" / "final_centroids.fasta"
            if final_out.exists():
                n = sum(1 for _ in open(final_out) if _.startswith('>'))
                self.log.info("  Rescue 最终输出: %d 条 HQ vOTU → %s", n, final_out)
            self.log.info("  Rescue 完成 → %s", rescue_out)
        else:
            self.log.warning("  Rescue 部分任务失败")

        # ── 合并 免拯救(known+CheckV-pass) + rescue → 完整病毒集合 ──
        all_plant = self.d['rescue_dir'] / "all_plant_viruses.fasta"
        no_rescue_fa = self.d['rescue_dir'] / "no_rescue" / "centroids" / "final_centroids.fasta"
        if not no_rescue_fa.is_file():
            no_rescue_fa = self.d['rescue_dir'] / "known" / "centroids" / "final_centroids.fasta"
        with open(all_plant, "w") as apf:
            n_no_rescue = 0
            if no_rescue_fa.is_file():
                for line in open(no_rescue_fa): apf.write(line)
                n_no_rescue = sum(1 for l in open(no_rescue_fa) if l.startswith('>'))
            n_rescued_total = 0
            for host in host_filter:
                safe_host = host.replace("/", "_").replace(" ", "_")
                rescue_final = self.d['rescue_dir'] / safe_host / "centroids" / "final_centroids.fasta"
                if rescue_final.is_file():
                    for line in open(rescue_final): apf.write(line)
                    n_rescued_total += sum(1 for l in open(rescue_final) if l.startswith('>'))
        self.log.info("=" * 50)
        self.log.info("  完整植物病毒: %d 条 → %s",
                      n_no_rescue + n_rescued_total, all_plant)
        self.log.info("    CD-HIT known: %d  |  CheckV pass(≥90%%): %d  |  rescued: %d",
                      n_cdhit, n_checkv, n_rescued_total)

        # ── CheckV 质量评估 (按宿主统计) ──
        self.log.info("=" * 50)
        self.log.info("[CheckV] 各宿主质量评估统计")
        self._run_checkv_summary()


    # ── CheckV 辅助 ──

    def _run_checkv_on_fasta(self, fasta_path, out_dir):
        """对单个 FASTA 运行 checkv completeness, 返回 {quality: count} 和总数"""
        fasta_path = Path(fasta_path)
        if not fasta_path.exists() or fasta_path.stat().st_size < 50:
            return {}, 0

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # checkv completeness 产出 completeness.tsv (不是 _quality_summary.tsv)
        completeness_tsv = out_dir / "completeness.tsv"
        if not self.args.force and completeness_tsv.is_file() and completeness_tsv.stat().st_size > 100:
            self.log.debug("  [SKIP] CheckV 已有结果: %s", completeness_tsv)
        else:
            cmd = (
                f"checkv completeness {fasta_path} {out_dir} "
                f"-d {self.args.checkv_db} -t {min(self.args.threads, 16)} "
            )
            ok, _ = run_cmd(cmd, self.log, f"CheckV: {fasta_path.name}", str(self.d['checkv_dir'] / "checkv.log"))
            if not ok:
                self.log.warning("  CheckV 失败: %s", fasta_path.name)
                return {}, 0

        if not completeness_tsv.is_file():
            return {}, 0

        try:
            df = pl.read_csv(str(completeness_tsv), separator="\t", null_values=["NA", "N/A", ""])
        except Exception as e:
            self.log.warning("  解析 CheckV 输出失败: %s", e)
            return {}, 0

        total = df.height
        # 检测 CheckV 质量列 (兼容新旧版本)
        quality_col = None
        for c in ["checkv_quality", "quality"]:
            if c in df.columns:
                quality_col = c
                break

        # 新版 CheckV (v1.x) 无 quality 列, 从 completeness 值推导
        if quality_col is None:
            comp_col = None
            for c in ["completeness", "aai_completeness"]:
                if c in df.columns:
                    comp_col = c
                    break
            if comp_col is None:
                self.log.warning("  CheckV 输出无可用的 quality/completeness 列")
                return {}, total

            dist = {}
            for row in df.iter_rows(named=True):
                raw = row.get(comp_col)
                try:
                    comp = float(raw) if raw and raw != "NA" else None
                except (ValueError, TypeError):
                    comp = None
                if comp is None:
                    key = "Not-determined"
                elif comp >= 90:
                    key = "Complete"
                elif comp >= 50:
                    key = "High-quality"
                elif comp >= 10:
                    key = "Medium-quality"
                else:
                    key = "Low-quality"
                dist[key] = dist.get(key, 0) + 1
        else:
            dist = {}
            for val in df[quality_col].to_list():
                key = val if val and val != "NA" else "Not-determined"
                dist[key] = dist.get(key, 0) + 1

        return dist, total


    def _run_checkv_summary(self):
        """对 rescue 产出 + unknown + skipped 全部运行 CheckV 并打印对比表"""
        host_filter = self.args.host_filter or ["Plant"]
        if isinstance(host_filter, str):
            host_filter = [h.strip() for h in host_filter.split(",")]

        all_stats = {}  # {label: (dist, total)}

        # 1. 目标宿主 rescue 产出
        for host in host_filter:
            safe_host = host.replace("/", "_").replace(" ", "_")
            rescue_final = self.d['rescue_dir'] / safe_host / "centroids" / "final_centroids.fasta"
            checkv_dir = self.d['rescue_dir'] / "checkv" / safe_host
            dist, total = self._run_checkv_on_fasta(rescue_final, checkv_dir)
            all_stats[f"Rescue_{host}"] = (dist, total)

        # 1.5 免拯救: CD-HIT known + CheckV pass (≥90%) — 无需经过 rescue
        no_rescue_fa = self.d['rescue_dir'] / "known" / "centroids" / "final_centroids.fasta"
        if no_rescue_fa.is_file():
            dist, total = self._run_checkv_on_fasta(no_rescue_fa, self.d['rescue_dir'] / "checkv" / "no_rescue")
            all_stats["免拯救(known+≥90%)"] = (dist, total)

        # 2. Unknown centroids
        unknown_fa = self.d['cluster'] / "centroids" / "unknown_votus.fasta"
        if unknown_fa.exists():
            dist, total = self._run_checkv_on_fasta(unknown_fa, self.d['rescue_dir'] / "checkv" / "unknown")
            all_stats["Unknown"] = (dist, total)

        # 3. 跳过的宿主 centroids
        seen_hosts = set()
        for f in (self.d['cluster'] / "centroids").glob("skipped_*.fasta"):
            host_label = f.stem.replace("skipped_", "").replace("_", " ")
            if host_label in seen_hosts:
                continue
            seen_hosts.add(host_label)
            dist, total = self._run_checkv_on_fasta(f, self.d['rescue_dir'] / "checkv" / f"skipped_{host_label}")
            all_stats[f"Skipped_{host_label}"] = (dist, total)

        # 4. 汇总 centroids
        all_centroids = self.d['cluster'] / "centroids" / "final_centroids.fasta"
        if all_centroids.exists():
            dist, total = self._run_checkv_on_fasta(all_centroids, self.d['rescue_dir'] / "checkv" / "all")
            all_stats["All_centroids"] = (dist, total)

        if not all_stats:
            self.log.info("  无 CheckV 结果")
            return

        # ── 统一质量等级排序 ──
        quality_order = ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]

        # 打印表头
        self.log.info("")
        self.log.info("  %-22s %10s %10s %10s %10s %10s %10s" % (
            "", "Complete", "High-qual", "Medium", "Low", "Not-det", "Total"))
        self.log.info("  " + "-" * 82)

        for label, (dist, total) in all_stats.items():
            parts = [f"  {label:<22s}"]
            for q in quality_order:
                parts.append(f"{dist.get(q, 0):>10d}")
            parts.append(f"{total:>10d}")
            self.log.info("".join(parts))

        self.log.info("  " + "-" * 82)

        # 汇总: 按来源分别统计
        self.log.info("  注: Complete ≥90% | High-quality ≥50% | Medium ≥10% | Low <10% | Not-det 无法判断")

        # 免拯救
        no_rescue_dist = all_stats.get("免拯救(known+≥90%)", ({}, 0))
        no_rescue_total = no_rescue_dist[1]
        no_rescue_hq = no_rescue_dist[0].get("Complete", 0) + no_rescue_dist[0].get("High-quality", 0)
        self.log.info("  免拯救 HQ (Complete+High): %d / %d (%.1f%%)",
                     no_rescue_hq, no_rescue_total,
                     no_rescue_hq / no_rescue_total * 100 if no_rescue_total else 0)

        # Rescue 产出
        rescue_keys = [k for k in all_stats if k.startswith("Rescue_")]
        if rescue_keys:
            rescue_total = sum(v[1] for k, v in all_stats.items() if k in rescue_keys)
            rescue_hq = sum(v[0].get("Complete", 0) + v[0].get("High-quality", 0)
                           for k, v in all_stats.items() if k in rescue_keys)
            self.log.info("  Rescue 产出 HQ (Complete+High): %d / %d (%.1f%%)",
                         rescue_hq, rescue_total,
                         rescue_hq / rescue_total * 100 if rescue_total else 0)

        # 全部植物病毒
        plant_total = no_rescue_total + rescue_total
        plant_hq = no_rescue_hq + rescue_hq
        self.log.info("  ★ 全部植物病毒 HQ: %d / %d (%.1f%%)",
                     plant_hq, plant_total,
                     plant_hq / plant_total * 100 if plant_total else 0)

    # ── Step 5: 分类注释 ──
    def run_taxonomy(self):
        self.d['taxonomy'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[5] 病毒分类注释 (virus_classifier2.py + R 共识整合)")

        centroids = self.d['cluster'] / "centroids" / "final_centroids.fasta"
        if not centroids.exists():
            self.log.warning("未找到 centroids, 跳过分类")
            return

        tax_dir = self.d['taxonomy']
        sample = "WVDB_votus"
        out_subdir = tax_dir / f"{sample}.virus_classed"
        combined_tsv = out_subdir / f"{sample}_combined_taxonomy.tsv"
        int_dir = tax_dir / "integrated"
        final_tsv = int_dir / "final_integrated_classification.tsv"

        # virus_classifier2.py 完整参数:
        #   -g, -s, -t, -o, -p, -f, --db-dir,
        #   --genomad-db, --metabuli-db, --cat-db, --cat-tax, --uniprot-db,
        #   --mmseqs-db, --vitap-db, --acvirus-db, --vcontact3-db,
        #   -j, -e, --remove-suffix
        if not self.args.force and combined_tsv.is_file() and combined_tsv.stat().st_size > 100:
            self.log.info("  [SKIP] virus_classifier2 — 已有结果")
        else:
            parts = [
                f"python {self.sc['classifier']}",
                f"-g {centroids}",
                f"-s {sample}",
                f"-t {self.args.tax_tools}",
                f"-o {tax_dir}",
                f"-p {self.args.threads}",
                f"-j {self.args.tax_jobs}",
                f"--db-dir {self.args.virus_db}",
                f"-e {self.args.tax_ext}",
            ]
            if self.args.uniprot_db:
                parts.append(f"--uniprot-db {self.args.uniprot_db}")
            if self.args.metabuli_db:
                parts.append(f"--metabuli-db {self.args.metabuli_db}")
            if self.args.genomad_db:
                parts.append(f"--genomad-db {self.args.genomad_db}")
            if self.args.cat_db:
                parts.append(f"--cat-db {self.args.cat_db}")
            if self.args.cat_tax:
                parts.append(f"--cat-tax {self.args.cat_tax}")
            if self.args.mmseqs_db:
                parts.append(f"--mmseqs-db {self.args.mmseqs_db}")
            if self.args.vitap_db:
                parts.append(f"--vitap-db {self.args.vitap_db}")
            if self.args.acvirus_db:
                parts.append(f"--acvirus-db {self.args.acvirus_db}")
            if self.args.vcontact3_db:
                parts.append(f"--vcontact3-db {self.args.vcontact3_db}")
            if self.args.tax_remove_suffix:
                parts.append(f"--remove-suffix {self.args.tax_remove_suffix}")
            if self.args.force:
                parts.append("-f")
            ok, _ = run_cmd(' '.join(parts), self.log, "virus_classifier2.py", str(self.d['taxonomy'] / "taxonomy.log"))
            if not ok:
                self.log.warning("  virus_classifier2 失败")

        # Step 2: R 共识整合
        if not self.args.force and final_tsv.is_file() and final_tsv.stat().st_size > 100:
            self.log.info("  [SKIP] R consensus — 已有结果")
        elif combined_tsv.is_file():
            int_dir.mkdir(parents=True, exist_ok=True)
            parts = [
                "Rscript", str(self.sc['classifier_R']),
                "--combined", str(combined_tsv),
                "--output", str(int_dir),
            ]
            ok, _ = run_cmd(' '.join(parts), self.log, "R consensus", str(self.d['taxonomy'] / "r_consensus.log"))
            if not ok:
                self.log.warning("  R consensus 失败")

        if final_tsv.is_file():
            n = sum(1 for _ in open(final_tsv)) - 1
            self.log.info("  分类完成: %d 条 → %s", n, final_tsv)
        else:
            self.log.warning("  分类未产出最终结果")

    # ── Step 6: 宿主预测 ──
    def run_host(self):
        self.d['host_pred'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[6] 宿主预测 (RNAVirHost + PhaBOX2 + ICTV, ICTV > RVH > PB2)")

        centroids = self.d['cluster'] / "centroids" / "final_centroids.fasta"
        tax_tsv = self.d['taxonomy'] / "integrated" / "final_integrated_classification.tsv"
        if not centroids.exists():
            self.log.warning("未找到 centroids, 跳过宿主预测")
            return
        if not tax_tsv.exists():
            self.log.error("分类结果 %s 不存在, 请先运行 --stage taxonomy", tax_tsv)
            return

        # run_host_prediction.py 完整参数:
        #   -i, --tax, -o, -t, --phabox-db, --prob-dir,
        #   --mode, -f, --skip-rnavirhost, --skip-phabox, --skip-ictv
        parts = [
            f"python {self.sc['host_pred']}",
            f"-i {centroids}",
            f"--tax {tax_tsv}",
            f"-o {self.d['host_pred']}",
            f"-t {self.args.threads}",
            f"--mode all",
        ]
        if self.args.phabox_db:
            parts.append(f"--phabox-db {self.args.phabox_db}")
        if self.args.prob_dir:
            parts.append(f"--prob-dir {self.args.prob_dir}")
        if self.args.force:
            parts.append("-f")
        if self.args.skip_rnavirhost:
            parts.append("--skip-rnavirhost")
        if self.args.skip_phabox:
            parts.append("--skip-phabox")
        if self.args.skip_ictv:
            parts.append("--skip-ictv")

        ok, _ = run_cmd(' '.join(parts), self.log, "Host Prediction", str(self.d['host_pred'] / "host.log"))
        if ok:
            self.log.info("  宿主预测完成 → %s", self.d['host_pred'])
        else:
            self.log.warning("  宿主预测失败")

    # ── CheckV 预评估: 按宿主分类检查 centroids 完整性 ──
    def run_checkv_stage(self):
        self.d['checkv_dir'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[CheckV] 按宿主分类预评估 centroids 完整性")

        host_summary = self.d['host_pred'] / "ensemble_host_summary.tsv"
        centroids_fa = self.d['cluster'] / "centroids" / "final_centroids.fasta"

        if not host_summary.exists():
            self.log.warning("宿主预测结果不存在, 跳过 CheckV 预评估")
            return
        if not centroids_fa.exists():
            self.log.warning("centroids 不存在, 跳过 CheckV 预评估")
            return

        host_df = pl.read_csv(str(host_summary), separator="\t", null_values=["NA", "N/A", ""])
        centroids_map = {}
        for rec in SeqIO.parse(str(centroids_fa), "fasta"):
            centroids_map[rec.id] = rec

        # 按 Final_Host 分组
        host_groups = {}
        for row in host_df.iter_rows(named=True):
            cid = row["contig_id"]
            fh = row.get("Final_Host", "Unknown")
            if not fh or fh == "NA":
                fh = "Unknown"
            host_groups.setdefault(fh, []).append(cid)

        checkv_dir = self.d['checkv_dir']
        checkv_dir.mkdir(parents=True, exist_ok=True)
        checkv_pass_ids = set()

        self.log.info("  按宿主分组运行 CheckV:")
        for host, ids in sorted(host_groups.items()):
            safe_host = host.replace("/", "_").replace(" ", "_")
            host_fa = checkv_dir / f"{safe_host}.fasta"
            with open(host_fa, "w") as hf:
                for cid in ids:
                    if cid in centroids_map:
                        SeqIO.write(centroids_map[cid], hf, "fasta")

            # 运行 CheckV
            cv_out = checkv_dir / safe_host
            cv_out.mkdir(exist_ok=True)
            cv_tsv = cv_out / "completeness.tsv"
            if not self.args.force and cv_tsv.is_file() and cv_tsv.stat().st_size > 100:
                self.log.debug("    [SKIP] %s 已有结果", host)
            else:
                cmd = (f"checkv completeness {host_fa} {cv_out} "
                       f"-d {self.args.checkv_db} -t {min(self.args.threads, 16)}")
                run_cmd(cmd, self.log, f"CheckV: {host}", str(self.d['checkv_dir'] / "checkv.log"))

            # 解析结果, 标记 pass (>90%)
            if cv_tsv.is_file():
                try:
                    cv_df = pl.read_csv(str(cv_tsv), separator="\t", null_values=["NA", "N/A", ""])
                    comp_col = "aai_completeness" if "aai_completeness" in cv_df.columns else "completeness"
                    n_complete = 0
                    for row in cv_df.iter_rows(named=True):
                        cid = row.get("contig_id", "")
                        val = row.get(comp_col, 0)
                        if val is not None and float(val) >= 90.0:
                            checkv_pass_ids.add(cid)
                            n_complete += 1
                    quality_col = "checkv_quality" if "checkv_quality" in cv_df.columns else None
                    if quality_col:
                        qdist = cv_df.group_by(quality_col).agg(pl.len()).to_dict(as_series=False)
                        qstr = ", ".join(f"{k}={list(v)[0]}" for k, v in zip(qdist[quality_col], qdist["len"]))
                        self.log.info("    %s: %d 条, Complete=%d | %s", host, len(ids), n_complete, qstr)
                    else:
                        self.log.info("    %s: %d 条, Complete(≥90%%)=%d", host, len(ids), n_complete)
                except Exception as e:
                    self.log.warning("    %s: 解析失败 - %s", host, e)

        # 写入 checkv_pass_ids 供 rescue 阶段使用
        pass_file = self.d['checkv_dir'] / "checkv_pass_ids.txt"
        with open(pass_file, "w") as pf:
            for cid in sorted(checkv_pass_ids):
                pf.write(f"{cid}\n")

        self.log.info("  CheckV pass (≥90%%): %d 条 → %s", len(checkv_pass_ids), pass_file)
        self.log.info("  按宿主 CheckV 报告 → %s", checkv_dir)

    # ── 汇总 ──
    @staticmethod
    def _asm_stats(fasta_path):
        """返回 (n, total, max_len, n50, n90, c500, r500, c1000, r1000)"""
        lens = []; seq = ""
        for line in open(fasta_path):
            l = line.strip()
            if l.startswith('>'):
                if seq: lens.append(len(seq))
                seq = ""
            else: seq += l
        if seq: lens.append(len(seq))
        if not lens: return (0,0,0,0,0,0,0,0,0)
        lens.sort(reverse=True); total = sum(lens); cum = 0
        half = total / 2; n90t = total * 0.9; n50 = n90 = 0
        for la in lens:
            cum += la
            if n50 == 0 and cum >= half: n50 = la
            if n90 == 0 and cum >= n90t: n90 = la
        c500 = sum(1 for la in lens if la > 500)
        c1000 = sum(1 for la in lens if la > 1000)
        r500 = round(c500/max(len(lens),1)*100, 1)
        r1000 = round(c1000/max(len(lens),1)*100, 1)
        return (len(lens), total, lens[0], n50, n90, c500, r500, c1000, r1000)

    def run_reports(self):
        """调用独立报告生成脚本 report_pipeline.py"""
        self.log.info("=" * 50)
        self.log.info("[9] Virome Report — 流水线总结报告")
        report_script = SCRIPT_DIR / "report_pipeline.py"
        if not report_script.is_file():
            self.log.error("  report_pipeline.py not found at %s", report_script)
            return
        cmd = [sys.executable, str(report_script), "-o", str(self.d['root'])]
        self.log.info("  → %s", " ".join(cmd))
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if result.stdout:
                for line in result.stdout.strip().split('\n')[-30:]:
                    self.log.info("  %s", line)
            if result.returncode != 0:
                self.log.error("  Report generation failed (rc=%d)", result.returncode)
                if result.stderr:
                    for line in result.stderr.strip().split('\n')[-8:]:
                        self.log.error("  %s", line)
        except Exception as e:
            self.log.error("  Report generation error: %s", e)
def print_help_and_exit():
    print(OVERVIEW)
    sys.exit(0)


def _build_parser(add_help=True):
    """构建 argparse, 可复用"""
    p = argparse.ArgumentParser(
        description='宏病毒组端到端全自动主控流水线', add_help=add_help)

    g = p.add_argument_group('路径配置')
    g.add_argument('--input_reads', help='原始 FASTQ 目录 (--cluster_input 模式可省略)')
    g.add_argument('--output_dir', required=True, help='项目输出根目录')
    g.add_argument('--cluster_input', help='直接输入已合并的病毒 FASTA (跳过 COBRA 收集)')
    g.add_argument('--input_assembly', help='组装结果目录 (identification 阶段, 默认 01_Assembly/)')

    g = p.add_argument_group('流程控制')
    g.add_argument('--stage', default=['all'], nargs='+',
                   choices=['all', 'clean', 'deplete', 'assembly', 'identification', 'cobra', 'cluster', 'taxonomy', 'host', 'checkv', 'rescue', 'report'],
                   help='运行阶段 (可多个, 如: --stage clean deplete)')
    g.add_argument('--host-filter', default='Plant',
                   help='目标宿主 (逗号分隔, rescue 阶段使用, 默认: Plant. Unknown 默认跳过并输出到 unknown_votus.fasta)')
    g.add_argument('--skip_clean', action='store_true', help='跳过数据清洗')
    g.add_argument('--skip_depletion', action='store_true', help='跳过去宿主')
    g.add_argument('--skip_clumpify', action='store_true', help='跳过 Clumpify')
    g.add_argument('--force', action='store_true', help='强制重跑')
    g.add_argument('--dry-run', action='store_true', help='仅扫描样本并显示配置，不实际执行')
    g.add_argument('--log-level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                   help='日志级别 (默认: INFO)')
    g.add_argument('--stop-on-error', action='store_true',
                   help='子脚本失败时立即终止 (默认: 仅警告, 继续执行后续阶段)')

    g = p.add_argument_group('Clean 阶段 (clean-data.py)')
    g.add_argument('--dedup', action='store_true', help='启用 fastp 自带去重')
    g.add_argument('--clumpify_memory', default='10g', help='clumpify Java 堆内存 (默认: 10g)')
    g.add_argument('--no_compress', action='store_true', help='最终结果不使用 gzip 压缩')
    g.add_argument('--clean_debug', action='store_true', help='clean-data.py 详细调试日志')

    g = p.add_argument_group('Deplete 阶段 (host_depletion.py)')
    g.add_argument('--deplete_tmp', help='去宿主临时文件目录')
    g.add_argument('--kraken2_confidence', type=float, default=0.4, help='Kraken2 分类置信度阈值 (默认: 0.4)')
    g.add_argument('--keep_rrna', action='store_true', help='保留分离出的 rRNA 序列到 rrna/ 目录')
    g.add_argument('--rrna_chunk_size', type=int, default=256, help='ribodetector_cpu chunk_size (默认: 256)')
    g.add_argument('--rrna_report', default='ribodetector.report.txt', help='rRNA 统计报告文件名')
    g.add_argument('--deplete_steps', default='kraken2,align,rrna',
                   help='去宿主消融实验步骤 (默认: kraken2,align,rrna)')
    g.add_argument('--align_config', default='', help='透传给比对工具的额外参数')
    g.add_argument('--deplete_debug', action='store_true', help='host_depletion.py 详细调试日志')

    g = p.add_argument_group('Assembly 阶段 (assembly_pipeline.py)')
    g.add_argument('--refinec_threads', type=int, help='refineC 独立线程数 (默认: 使用 --threads)')
    g.add_argument('--refinec_frag_min_len', type=int, default=1000, help='refineC split 最小片段长度 bp (默认: 1000)')
    g.add_argument('--refinec_min_id', type=float, default=0.95, help='refineC merge 最小序列一致性 (默认: 0.95)')
    g.add_argument('--refinec_min_cov', type=float, default=0.50, help='refineC merge 最小覆盖度 (默认: 0.50)')
    g.add_argument('--asm_tmp_dir', help='组装临时文件目录')
    g.add_argument('--asm_keep_temp', action='store_true', help='保留组装临时文件及 refineC 中间目录')

    g = p.add_argument_group('Identification 阶段 (virus_identification16.py)')
    g.add_argument('--nr_db', help='Diamond NR 数据库路径')
    g.add_argument('--skip_uniprot_filter', action='store_true', help='跳过 UniProt 后置过滤')
    g.add_argument('--skip_nr_filter', action='store_true', help='跳过 NR 后置过滤')
    g.add_argument('--skip_id_plots', action='store_true', help='跳过鉴定阶段图表生成')
    g.add_argument('--clean_failed', action='store_true', help='自动清理鉴定失败的任务目录')
    g.add_argument('--ident_ext', default='.fasta', help='输入目录时搜索的后缀 (默认: .fasta)')

    g = p.add_argument_group('COBRA 阶段 (cobra_pipeline.py)')
    g.add_argument('--cobra_mink', type=int, default=21, help='COBRA 最小 kmer (默认: 21)')
    g.add_argument('--cobra_maxk', type=int, default=141, help='COBRA 最大 kmer (默认: 141)')
    g.add_argument('--cobra_linkage_mismatch', type=int, default=2, help='COBRA 链接识别不匹配数 (默认: 2)')
    g.add_argument('--cobra_verbose', action='store_true', help='cobra_pipeline.py 详细日志')

    g = p.add_argument_group('CLUSTER 阶段 (cluster_pipeline.py)')
    g.add_argument('--skip_vclust', action='store_true', help='跳过 vclust 聚类步骤')
    g.add_argument('--vclust_cluster_file', help='复用已有 vclust 聚类 TSV 文件')

    g = p.add_argument_group('Taxonomy 阶段 (virus_classifier2.py)')
    g.add_argument('--tax_tools', default='all', help='分类工具: genomad,metabuli,diamond_lca,VITAP,mmseqs,ACVirus,vcontact3,PhaGCN3,all (默认: all)')
    g.add_argument('--tax_jobs', type=int, default=1, help='分类并行任务数 (默认: 1)')
    g.add_argument('--tax_ext', default='.fasta', help='分类输入文件扩展名 (默认: .fasta)')
    g.add_argument('--tax_remove_suffix', help='分类输入文件去后缀名')

    g = p.add_argument_group('Host 阶段 (run_host_prediction.py)')
    g.add_argument('--skip_rnavirhost', action='store_true', help='跳过 RNAVirHost 宿主预测')
    g.add_argument('--skip_phabox', action='store_true', help='跳过 PhaBOX2 宿主预测')
    g.add_argument('--skip_ictv', action='store_true', help='跳过 ICTV 宿主查找')

    g = p.add_argument_group('数据库路径')
    g.add_argument('--host_db', help='宿主数据库根目录 (自动查找 kraken2/ bowtie2/ hisat2/ minimap2/)')
    g.add_argument('--kraken2_db', help='Kraken2 宿主库 (覆盖 --host_db 自动检测)')
    g.add_argument('--host_align_db', help='宿主比对索引 (覆盖 --host_db 自动检测)')
    g.add_argument('--virus_db', help='病毒鉴定数据库根目录')
    g.add_argument('--checkv_db', help='CheckV 数据库路径')

    g = p.add_argument_group('工具与算法')
    g.add_argument('--aligner', default='bowtie2', choices=['bowtie2', 'hisat2', 'minimap2'])
    g.add_argument('--seq_type', default='rna-short', choices=['dna-short', 'rna-short', 'nanopore', 'pacbio'])
    g.add_argument('--rrna', action='store_true', help='开启 rRNA 剔除')
    g.add_argument('--rrna_tool', default='ribodetector', choices=['ribodetector', 'silva'],
                   help='rRNA 剔除工具: ribodetector (默认) / silva (Bowtie2+SILVA)')
    g.add_argument('--silva_index', help='SILVA Bowtie2 索引前缀 (--rrna_tool silva 时必需)')
    g.add_argument('--assembler', default='penguin', choices=['megahit', 'rnaviralspades', 'penguin', 'all'])
    g.add_argument('--contig-length', '-l', type=int, default=200, help='contig 最小长度 bp (默认 200)')
    g.add_argument('--identify_tools', default='all', help='病毒鉴定工具')
    g.add_argument('--virus_mode', default='strict', choices=['raw', 'filter', 'strict'],
                   help='COBRA 病毒序列来源: raw=原始鉴定, filter=UniProt过滤, strict=严格过滤 (默认: strict)')
    g = p.add_argument_group('鉴定数据库 (virus_identification16.py)')
    g.add_argument('--virus_protein_db', help='病毒蛋白 Diamond DB')
    g.add_argument('--uniprot_db', help='UniProt Diamond DB')
    g.add_argument('--viroids_db', help='类病毒 BLAST DB')
    g.add_argument('--virsorter_db', help='VirSorter2 数据库')
    g.add_argument('--viralverify_hmm', help='ViralVerify HMM 文件')
    g.add_argument('--metabuli_db', help='Metabuli 数据库')
    g.add_argument('--virus_taxid', help='病毒 TaxID 列表')
    g.add_argument('--virhunter_path', help='VirHunter predict_cpu.py 路径')
    g.add_argument('--virhunter_weights', help='VirHunter weights 目录')
    g.add_argument('--virbot_path', help='VirBot.py 路径')
    g.add_argument('--viralm_path', help='viralm_cpu.py 路径')
    g.add_argument('--virsorter_group', default='dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae')
    g.add_argument('--blast_mode', default='both', help='Blast 模式 (默认 both: 同时产出 filter+strict)')
    g.add_argument('--blast_evalue', default='1e-5', help='Blast e-value (默认 1e-5)')
    g.add_argument('--blast_top_n', default='5', help='Blast top N (默认 5)')
    g.add_argument('--phabox-db', help='PhaBOX2 数据库路径 (host 阶段)')
    g.add_argument('--prob-dir', help='ICTV 宿主概率表目录 (host 阶段, 默认: cross_analysis/)')
    g.add_argument('--ref-genomes', nargs='*', help='ICTV/NCBI 参考基因组 FASTA (可多个, CD-HIT 参考引导预聚类)')
    g.add_argument('--blast-db', help='BLAST 参考数据库 (rescue 阶段分支 D, ref.fasta)')

    g = p.add_argument_group('计算资源')
    g.add_argument('--threads', '-t', type=int, default=20, help='线程 (默认 20)')
    g.add_argument('--memory', '-m', type=int, default=64, help='内存 GB (默认 64)')
    g.add_argument('--jobs', '-j', type=int, default=2, help='并行数 (默认 2)')

    g = p.add_argument_group('CLUSTER 参数')
    g.add_argument('--min-length', type=int, default=500, help='病毒最小长度 bp (默认 500)')
    g.add_argument('--ani', type=float, default=0.95, help='vclust ANI 阈值 (默认 0.95)')
    g.add_argument('--qcov', type=float, default=0.85, help='vclust qcov 阈值 (默认 0.85)')
    g.add_argument('--cdhit_ani', type=float, help='CD-HIT ANI 阈值 (默认 0.95, 转录组建议 0.85)')
    g.add_argument('--cdhit_qcov', type=float, help='CD-HIT qcov 阈值 (默认 0.85, 转录组建议 0.50)')
    g.add_argument('--virseqimprover-path', help='Virseqimprover.py 路径')
    g.add_argument('--salmon-bin', default=os.path.expanduser('~/mambaforge/envs/Virseqimprover/bin/salmon'), help='Salmon 二进制路径')
    g.add_argument('--max_vsi_samples', type=int, default=10, help='VSI 最大合并样本数 (0=不限制, 默认: 10)')
    g.add_argument('--min_vsi_len', type=int, default=2000, help='VSI 最小 contig 长度 bp (默认: 2000)')
    g.add_argument('--checkv_threshold', type=float, default=90.0, help='CheckV completeness 通过阈值 (默认90, 植物病毒建议80)')

    g = p.add_argument_group('分类数据库 (virus_classifier2.py)')
    g.add_argument('--genomad_db', help='genomad DB 路径')
    g.add_argument('--cat_db', help='CAT 数据库路径')
    g.add_argument('--cat_tax', help='CAT taxonomy 路径')
    g.add_argument('--mmseqs_db', help='mmseqs 数据库路径')
    g.add_argument('--vitap_db', help='VITAP 数据库路径')
    g.add_argument('--acvirus_db', help='ACVirus 数据库路径')
    g.add_argument('--vcontact3_db', help='vConTACT3 数据库路径')
    return p


def parse_args():
    for i, a in enumerate(sys.argv[1:], 1):
        if a == '--help-all':
            p = _build_parser(add_help=True)
            p.parse_args(['--help'])
            sys.exit(0)
        if a in ('-h', '--help'):
            # 如果同时有 --stage <name> -h, 显示阶段详情
            for j, b in enumerate(sys.argv[1:], 1):
                if b == '--stage' and j < len(sys.argv) - 1:
                    stage = sys.argv[j + 1]
                    if stage in STAGE_HELP:
                        print(STAGE_HELP[stage])
                        sys.exit(0)
            print_help_and_exit()

    if len(sys.argv) == 1:
        print_help_and_exit()

    return _build_parser(add_help=False).parse_args()
def main():
    args = parse_args()
    logger = setup_logger(args.output_dir, level=args.log_level)

    stages = set(args.stage)  # 支持多阶段: --stage clean deplete

    logger.info("=" * 50)
    logger.info("Virome Pipeline v2.3")
    logger.info("  Stage:    %s", ','.join(sorted(stages)))
    logger.info("  Output:   %s", args.output_dir)

    # 根据 stages 自动设置 skip 标志
    needs_reads = bool(stages & {'clean','deplete','assembly','cobra','rescue'})
    if stages == {'clean'}:
        args.skip_clean = False
        args.skip_depletion = True
        logger.info("  Flow:     Clean only (清洗)")
    elif stages == {'deplete'}:
        args.skip_clean = True
        args.skip_depletion = False
        logger.info("  Flow:     Deplete only (去宿主)")
    elif not needs_reads:
        args.skip_clean = True
        args.skip_depletion = True
        logger.info("  Flow:     %s (无需 reads)", ','.join(sorted(stages)))
    elif 'all' not in stages:
        args.skip_clean = 'clean' not in stages
        args.skip_depletion = 'deplete' not in stages
        logger.info("  Flow:     %s", ','.join(sorted(stages)))
    else:  # all
        flow = []
        flow.append('SKIP' if args.skip_clean else 'Clean')
        flow.append('SKIP' if args.skip_depletion else 'Deplete')
        flow.append(f'Assemble({args.assembler})')
        flow.append(f'Identify({args.identify_tools})')
        flow.append('COBRA')
        flow.append('Cluster(vclust)')
        flow.append('Taxonomy')
        flow.append('Host')
        flow.append('CheckV')
        flow.append(f'Rescue({args.host_filter})')
        logger.info("  Flow:     %s", ' → '.join(flow))
    logger.info("  Log Level: %s", args.log_level)
    if args.stop_on_error:
        logger.info("  模式:      遇错即停 (--stop-on-error)")
    else:
        logger.info("  模式:      容错继续 (默认)")
    logger.info("=" * 50)

    # ── --dry-run: 扫描样本后退出 ──
    if args.dry_run:
        logger.info("")
        logger.info("═══ DRY-RUN 模式 — 不执行任何计算 ═══")
        if args.input_reads and Path(args.input_reads).exists():
            samples = scan_samples_in_dir(args.input_reads)
            if samples:
                pe = sum(1 for v in samples.values() if v['r2'])
                logger.info("  输入目录: %s", args.input_reads)
                logger.info("  检测样本: %d (PE=%d, SE=%d)", len(samples), pe, len(samples) - pe)
                for name, info in sorted(samples.items()):
                    tag = "PE" if info['r2'] else "SE"
                    logger.info("    [%s] %s → %s", tag, name, info['r1'])
            else:
                logger.info("  [WARN] 输入目录无序列文件")
        else:
            logger.info("  输入目录: %s (不存在或未指定)", args.input_reads)
        logger.info("  阶段:     %s", stage)
        logger.info("  输出目录: %s", args.output_dir)
        logger.info("═══ DRY-RUN 结束 ═══")
        return

    pipe = ViromePipeline(args, logger)

    needs_reads = bool(stages & {'clean','deplete','assembly','cobra','rescue'})
    if needs_reads:
        # 优先使用 --input_reads 指定的目录
        if args.input_reads and Path(args.input_reads).exists():
            pipe.reads_dir = Path(args.input_reads)
        elif 'deplete' in stages:
            # deplete: 自动使用 clean 输出 (优先 clumpify, 否则 fasta)
            cl = pipe.d['clean'] / '3.clumpify'
            fa = pipe.d['clean'] / '2.fasta'
            pipe.reads_dir = cl if (cl.exists() and any(cl.iterdir())) else fa
        else:
            pipe.reads_dir = pipe.d['hostdep']
        if not pipe.reads_dir.exists() or not any(pipe.reads_dir.iterdir()):
            if stages <= {'cluster', 'taxonomy', 'host', 'checkv', 'report'}:
                logger.info("  无需 reads, 跳过")
            else:
                logger.error("reads 目录 %s 为空", pipe.reads_dir)
                sys.exit(1)
        logger.info("  使用 Reads: %s", pipe.reads_dir)
        pipe.orig_samples = scan_samples_in_dir(pipe.reads_dir)
        if not pipe.orig_samples:
            logger.error("在 %s 中未找到序列文件!", pipe.reads_dir)
            sys.exit(1)
        pe = sum(1 for v in pipe.orig_samples.values() if v['r2'])
        logger.info("  检测到 %d 个样本 (PE=%d, SE=%d)", len(pipe.orig_samples), pe,
                     len(pipe.orig_samples) - pe)

    # ═══ 执行 ═══
    _all = 'all' in stages
    stage_map = {
        'clean': pipe.run_clean, 'deplete': pipe.run_depletion,
        'assembly': pipe.run_assembly, 'identification': pipe.run_identification,
        'cobra': pipe.run_cobra, 'cluster': pipe.run_cluster,
        'taxonomy': pipe.run_taxonomy, 'host': pipe.run_host,
        'checkv': pipe.run_checkv_stage, 'rescue': pipe.run_rescue, 'report': pipe.run_reports,
    }
    # 按流水线顺序排列
    stage_order = ['clean','deplete','assembly','identification','cobra','cluster',
                   'taxonomy','host','checkv','rescue','report']
    stages_to_run = [(s, stage_map[s]) for s in stage_order if _all or s in stages]

    # 准备阶段日志目录
    stage_log_dir = Path(args.output_dir) / "09_Reports" / "logs"
    stage_log_dir.mkdir(parents=True, exist_ok=True)

    failed_stages = []
    for stage_name, stage_func in stages_to_run:
        # 添加阶段独立日志 handler
        stage_handler = logging.FileHandler(str(stage_log_dir / f"{stage_name}.log"), encoding='utf-8')
        stage_handler.setLevel(logging.DEBUG)
        stage_handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s', datefmt='%H:%M:%S'))
        logger.addHandler(stage_handler)
        try:
            stage_func()
        except SystemExit as e:
            if args.stop_on_error:
                logger.error("[%s] 阶段失败 (exit=%d), 终止", stage_name, e.code if e.code else 1)
                sys.exit(e.code if e.code else 1)
            else:
                logger.warning("[%s] 阶段失败 (exit=%d), 继续", stage_name, e.code if e.code else 1)
                failed_stages.append(stage_name)
        except Exception as e:
            if args.stop_on_error:
                logger.error("[%s] 阶段异常: %s, 终止", stage_name, e)
                sys.exit(1)
            else:
                logger.warning("[%s] 阶段异常: %s, 继续", stage_name, e)
                failed_stages.append(stage_name)
        finally:
            logger.removeHandler(stage_handler)

    if failed_stages:
        logger.warning("以下阶段失败: %s", ', '.join(failed_stages))


if __name__ == '__main__':
    main()
