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


def run_cmd(cmd, logger, step_name):
    """执行 shell 命令。stdout/stderr 透传 (进度条可见)。"""
    logger.info("[%s] 执行中...", step_name)
    logger.debug("  CMD: %s", cmd)
    try:
        subprocess.run(cmd, shell=True, check=True)
        logger.info("[%s] ✓ 成功", step_name)
        return True, ""
    except subprocess.CalledProcessError as e:
        logger.error("[%s] ✗ 失败 (exit=%d)", step_name, e.returncode)
        return False, f"exit={e.returncode}"


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

        ok, _ = run_cmd(' '.join(parts), self.log, "Clean")
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

        ok, _ = run_cmd(' '.join(parts), self.log, "Host Depletion")
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

        ok, _ = run_cmd(' '.join(parts), self.log, "Assembly")
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
        run_cmd(' '.join(parts), self.log, "VirusIdentification")

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

        ok, _ = run_cmd(' '.join(parts), self.log, "COBRA Pipeline")
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

        ok, _ = run_cmd(' '.join(parts), self.log, "CLUSTER (vclust only)")
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
        if not self.args.force:
            parts.append("--resume")

        ok, _ = run_cmd(' '.join(parts), self.log, f"Rescue ({','.join(host_filter)})")
        if ok:
            final_out = rescue_out / "centroids" / "final_centroids.fasta"
            if final_out.exists():
                n = sum(1 for _ in open(final_out) if _.startswith('>'))
                self.log.info("  Rescue 最终输出: %d 条 HQ vOTU → %s", n, final_out)
            self.log.info("  Rescue 完成 → %s", rescue_out)
        else:
            self.log.warning("  Rescue 部分任务失败")

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
            ok, _ = run_cmd(cmd, self.log, f"CheckV: {fasta_path.name}")
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

        # 1.5 CD-HIT known centroids (完整参考, 免 rescue)
        known_centroids = self.d['rescue_dir'] / "known" / "centroids" / "final_centroids.fasta"
        if known_centroids.is_file():
            dist, total = self._run_checkv_on_fasta(known_centroids, self.d['rescue_dir'] / "checkv" / "known_ref")
            all_stats["Known_ref"] = (dist, total)

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

        # 汇总: 所有 rescue 产出合并
        if any(k.startswith("Rescue_") for k in all_stats):
            self.log.info("  注: Complete ≥90% | High-quality ≥50% | Medium ≥10% | Low <10% | Not-det 无法判断")
            rescue_total = sum(v[1] for k, v in all_stats.items() if k.startswith("Rescue_"))
            rescue_complete = sum(v[0].get("Complete", 0) + v[0].get("High-quality", 0)
                                 for k, v in all_stats.items() if k.startswith("Rescue_"))
            self.log.info("  Rescue 汇总 HQ (Complete+High): %d / %d (%.1f%%)",
                         rescue_complete, rescue_total,
                         rescue_complete / rescue_total * 100 if rescue_total else 0)

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
            ok, _ = run_cmd(' '.join(parts), self.log, "virus_classifier2.py")
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
            ok, _ = run_cmd(' '.join(parts), self.log, "R consensus")
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

        ok, _ = run_cmd(' '.join(parts), self.log, "Host Prediction")
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
                run_cmd(cmd, self.log, f"CheckV: {host}")

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
        self.log.info("=" * 50)
        self.log.info("[9] Virome Report — 流水线总结报告")
        report_dir = self.d['reports']; report_dir.mkdir(parents=True, exist_ok=True)
        stage_stats = []

        # ── 辅助函数 ──
        def _count_fasta(path):
            if not path or not os.path.isfile(path): return 0
            return sum(1 for _ in open(path) if _.startswith('>'))

        def _count_lines(path):
            if not path or not os.path.isfile(path): return 0
            return sum(1 for _ in open(path))

        def _count_dir(d):
            return sum(1 for _ in Path(d).iterdir() if _.is_dir()) if d and Path(d).is_dir() else 0

        def _file_size(path):
            if not path or not os.path.isfile(path): return "0 B"
            s = os.path.getsize(path)
            for u in ['B', 'KB', 'MB', 'GB']:
                if s < 1024: return f"{s:.1f} {u}"
                s /= 1024
            return f"{s:.1f} TB"

        def _add(stage, status, key_metric="", details=""):
            stage_stats.append({"Stage": stage, "Status": status, "Key_Metric": key_metric, "Details": details})

        # ── 00a_CleanData ──
        clean = self.d['clean']
        if clean.is_dir():
            ns = _count_dir(clean / "3.clumpify") or _count_dir(clean / "2.fasta") or _count_dir(clean)
            _add("00a_CleanData", "✓", key_metric=f"{ns} 样本", details=f"{clean}")
            # fastp summary → data_summary.tsv (所有样本)
            import json
            fp = clean / "logs"
            if fp.is_dir():
                jf_list = list(fp.glob("*_fastp_report.json"))
                if jf_list:
                    n_total_before, n_total_after = 0, 0
                    with open(report_dir / "data_summary.tsv", "w") as ds:
                        hdr = "Sample\tRaw_Reads\tClean_Reads\tRetained(%)\tRaw_Q20(%)\tClean_Q20(%)\tRaw_Q30(%)\tClean_Q30(%)\tLowQ_Reads\tTooShort_Reads\tDup_Rate(%)"
                        ds.write(hdr + "\n")
                        for jf in sorted(jf_list):
                            try:
                                js = json.load(open(jf))
                                sn = jf.name.replace("_fastp_report.json", "")
                                bef = js.get("summary", {}).get("before_filtering", {})
                                aft = js.get("summary", {}).get("after_filtering", {})
                                fil = js.get("filtering_result", {})
                                dup = js.get("duplication", {})
                                nb = bef.get("total_reads", 0)
                                na = aft.get("total_reads", 0)
                                n_total_before += nb; n_total_after += na
                                ds.write(f"{sn}\t{nb}\t{na}\t{round(na/max(nb,1)*100,1)}\t"
                                         f"{round(bef.get('q20_rate',0)*100,1)}\t{round(aft.get('q20_rate',0)*100,1)}\t"
                                         f"{round(bef.get('q30_rate',0)*100,1)}\t{round(aft.get('q30_rate',0)*100,1)}\t"
                                         f"{fil.get('low_quality_reads',0)}\t{fil.get('too_short_reads',0)}\t"
                                         f"{round(dup.get('rate',0)*100,3)}\n")
                            except: pass
                        ds.write(f"TOTAL\t{n_total_before}\t{n_total_after}\t{round(n_total_after/max(n_total_before,1)*100,1)}\n")
                    _add("  └ data_summary", "✓", key_metric=f"reads: {n_total_before:,}→{n_total_after:,}",
                         details=f"{report_dir}/data_summary.tsv")
        else:
            _add("00a_CleanData", "○", details="未运行")

        # ── 00b_HostDepletion ──
        hostdep = self.d['hostdep']
        if hostdep.is_dir():
            ns = _count_dir(hostdep)
            _add("00b_HostDepletion", "✓", key_metric=f"{ns} 样本", details=f"{hostdep}")
        else:
            _add("00b_HostDepletion", "○", details="未运行")

        # ── 01_Assembly ──
        asm = self.d['asm']
        if asm.is_dir():
            ns = _count_dir(asm)
            total_contigs, total_bp = 0, 0
            for d in asm.iterdir():
                if not d.is_dir(): continue
                for f in d.glob("*.contig.fasta"):
                    for line in open(f):
                        if line.startswith('>'): total_contigs += 1
                        else: total_bp += len(line.strip())
            _add("01_Assembly", "✓", key_metric=f"{ns} 样本, {total_contigs:,} contigs, {total_bp/1e6:.1f} Mb",
                 details=f"{asm}")
            # per-sample detail
            for d in sorted(asm.iterdir()):
                if not d.is_dir(): continue
                n_contig = 0
                for f in d.glob("*.contig.fasta"):
                    n_contig = _count_fasta(f)
                _add(f"  └ {d.name}", "✓", key_metric=f"{n_contig} contigs")
            # 生成 assembly_summary.tsv (N50/N90)
            asm_data = []
            with open(report_dir / "assembly_summary.tsv", "w") as af:
                af.write("Sample\tAssembler\tSize(Mb)\tContigs\tMax_Len\tN50\tN90\t>500bp\t>500bp(%)\t>1000bp\t>1000bp(%)\n")
                for d in sorted(asm.iterdir()):
                    if not d.is_dir(): continue
                    for f in d.glob("*.contig.fasta"):
                        n, total, mx, n50, n90, c500, r500, c1000, r1000 = self._asm_stats(str(f))
                        if n == 0: continue
                        at = f.stem.replace(f"{d.name}_", "").replace(".contig", "")
                        af.write(f"{d.name}\t{at}\t{total/1e6:.1f}\t{n}\t{mx}\t{n50}\t{n90}\t{c500}\t{r500}\t{c1000}\t{r1000}\n")
                        asm_data.append({'s': d.name, 'n': n, 'total': total, 'n50': n50})
                if len(asm_data) > 1:
                    t_n = sum(r['n'] for r in asm_data)
                    t_bp = sum(r['total'] for r in asm_data)
                    af.write(f"TOTAL\tall\t{t_bp/1e6:.1f}\t{t_n}\t-\t-\t-\t-\t-\t-\t-\n")
        else:
            _add("01_Assembly", "○", details="未运行")

        # ── 02_Identification ──
        ident = self.d['ident']
        if ident.is_dir():
            ns = _count_dir(ident)
            n_virus = 0
            for d in ident.iterdir():
                if not d.is_dir(): continue
                for f in d.glob("*virus.all.candidate.fasta"):
                    n_virus += _count_fasta(f)
            _add("02_Identification", "✓", key_metric=f"{ns} 样本, {n_virus:,} 病毒序列", details=f"{ident}")
            # 生成 ident_summary.tsv + filter_summary.tsv
            tools_list = ['genomad','blast','metabuli','virsorter2','viralverify',
                          'virhunter','virbot','viralm','rdrpcatch']
            ident_data = []
            with open(report_dir / "ident_summary.tsv", "w") as ids:
                ids.write("Sample\tAll_Candidate\t" + "\t".join(tools_list) + "\n")
                for d in sorted(ident.iterdir()):
                    if not d.is_dir(): continue
                    all_ids = _count_fasta(d / f"{d.name}_virus.all.candidate.fasta")
                    tcounts = {}
                    for tool in tools_list:
                        idf = d / f"{d.name}_virus.{tool}.result.id"
                        tcounts[tool] = _count_lines(idf) if idf.is_file() else 0
                    ids.write(f"{d.name}\t{all_ids}\t" + "\t".join(str(tcounts[t]) for t in tools_list) + "\n")
                    ident_data.append({'Sample': d.name, 'All': all_ids, **tcounts})
                # TOTAL 行
                if len(ident_data) > 1:
                    total_all = sum(r['All'] for r in ident_data)
                    total_tools = {t: sum(r[t] for r in ident_data) for t in tools_list}
                    ids.write(f"TOTAL\t{total_all}\t" + "\t".join(str(total_tools[t]) for t in tools_list) + "\n")
            filter_data = []
            with open(report_dir / "filter_summary.tsv", "w") as fs:
                fs.write("Sample\tMode\tAll_Candidate\tPassed\tRetained(%)\n")
                for d in sorted(ident.iterdir()):
                    if not d.is_dir(): continue
                    all_n = _count_fasta(d / f"{d.name}_virus.all.candidate.fasta")
                    for fm, fd in [('filter','uniprot_filter_output_filter'),
                                   ('strict','uniprot_filter_output_strict'),
                                   ('comb','uniprot_filter_output')]:
                        ff = d / fd / f"{d.name}_virus.uniprot_filtered.fasta"
                        nf = _count_fasta(ff) if ff.is_file() else 0
                        if nf > 0 or (nf == 0 and fd == 'uniprot_filter_output'):
                            fs.write(f"{d.name}\t{fm}\t{all_n}\t{nf}\t{round(nf/max(all_n,1)*100,1)}\n")
                            filter_data.append({'Sample': d.name, 'Mode': fm, 'All': all_n, 'Passed': nf,
                                               'Retained': round(nf/max(all_n,1)*100,1)})
                # TOTAL 行 (按 mode)
                if len(filter_data) > 0:
                    modes_seen = set(r['Mode'] for r in filter_data)
                    for m in sorted(modes_seen):
                        mr = [r for r in filter_data if r['Mode'] == m]
                        t_all = sum(r['All'] for r in mr)
                        t_pass = sum(r['Passed'] for r in mr)
                        fs.write(f"TOTAL\t{m}\t{t_all}\t{t_pass}\t{round(t_pass/max(t_all,1)*100,1)}\n")
            # 多样本屏幕显示
            if len(ident_data) > 1:
                top_tools = sorted(total_tools.items(), key=lambda x: -x[1])[:3]
                best = " | ".join(f"{t}={c}" for t,c in top_tools)
                _add("  └ multi-sample", "✓", key_metric=f"total={total_all}, top={best}")
                for fm in sorted(modes_seen)[:2]:
                    mr = [r for r in filter_data if r['Mode'] == fm]
                    t_p = sum(r['Passed'] for r in mr)
                    _add(f"  └ {fm}", "✓", key_metric=f"{t_p} passed ({len(mr)} samples)")
        else:
            _add("02_Identification", "○", details="未运行")

        # ── 03_COBRA + cobra_summary.tsv ──
        cobra = self.d['cobra']
        if cobra.is_dir():
            ns = _count_dir(cobra)
            n_ext, n_queries, n_orphan = 0, 0, 0
            with open(report_dir / "cobra_summary.tsv", "w") as cs:
                cs.write("Sample\tTotal_Queries\tExtended_Circular\tExtended_Partial\t"
                         "Extended_Failed\tOrphan_End\tExtension_Rate(%)\tOrphan_Rate(%)\t"
                         "Extended_Contigs\tTotal_Gain(bp)\n")
                for sd in sorted(cobra.iterdir()):
                    if not sd.is_dir(): continue
                    sn = sd.name
                    # 统计 cobra.fa
                    sn_ext, sn_gain = 0, 0
                    for cf in sd.rglob("*.cobra.fa"):
                        if cf.stat().st_size > 0:
                            for rec in SeqIO.parse(str(cf), "fasta"):
                                sn_ext += 1
                                sn_gain += len(rec.seq)
                    n_ext += sn_ext
                    # 解析 COBRA log
                    logs = list(sd.rglob("log"))
                    cobra_logs = [f for f in logs if 'COBRA' in str(f.parent.name)]
                    if not cobra_logs: cobra_logs = logs[:1]
                    tq = ec = ep = ef = oe = 0
                    if cobra_logs:
                        try:
                            text = open(str(cobra_logs[0])).read()
                            for line in text.split('\n'):
                                s = line.strip()
                                if s.startswith('# Total queries:'): tq = int(s.split(':')[1].strip())
                                elif 'Self_circular' in s: pass
                                elif 'Extended_circular' in s: ec = int(s.split(':')[1].strip().split()[0])
                                elif 'Extended_partial' in s: ep = int(s.split(':')[1].strip().split()[0])
                                elif 'Extended_failed' in s: ef = int(s.split(':')[1].strip())
                                elif 'Orphan end' in s: oe = int(s.split(':')[1].strip())
                        except: pass
                    n_queries += tq; n_orphan += oe
                    er = round((ec+ep)/max(tq,1)*100, 1)
                    or_ = round(oe/max(tq,1)*100, 1)
                    cs.write(f"{sn}\t{tq}\t{ec}\t{ep}\t{ef}\t{oe}\t{er}\t{or_}\t{sn_ext}\t{sn_gain}\n")
            _add("03_COBRA", "✓", key_metric=f"{ns} 样本, {n_ext} 延伸, {n_queries} query, {n_orphan} orphan",
                 details=f"{report_dir}/cobra_summary.tsv")
        else:
            _add("03_COBRA", "○", details="未运行")

        # ── 04_CLUSTER ──
        cluster = self.d['cluster']
        centroids_fa = cluster / "centroids" / "final_centroids.fasta"
        known_fa = cluster / "2_cdhit" / "known_centroids.fasta"
        if cluster.is_dir():
            n_centroids = _count_fasta(centroids_fa) if centroids_fa.is_file() else 0
            n_known = _count_fasta(known_fa) if known_fa.is_file() else 0
            _add("04_CLUSTER", "✓", key_metric=f"{n_centroids:,} novel + {n_known:,} known centroids",
                 details=f"{cluster}")
            # vclust stats
            ctsv = cluster / "3_vclust" / "vclust_clusters.tsv"
            if ctsv.is_file():
                n_clusters = _count_lines(ctsv) - 1
                # 统计簇大小分布
                sizes = []
                with open(ctsv) as cf:
                    cf.readline()
                    for line in cf:
                        members = line.strip().split('\t')[1] if '\t' in line else ""
                        sizes.append(len(members.split(',')) if members else 0)
                singletons = sum(1 for sz in sizes if sz <= 1)
                max_sz = max(sizes) if sizes else 0
                _add("  └ vclust", "✓", key_metric=f"{n_clusters:,} 簇, {singletons} 单例, 最大簇={max_sz}")
            # CD-HIT known detail
            known_linked_fa = cluster / "2_cdhit" / "known_linked_centroids.fasta"
            n_linked = _count_fasta(known_linked_fa) if known_linked_fa.is_file() else 0
            if n_linked > 0:
                _add("  └ CD-HIT linked", "✓", key_metric=f"{n_linked} 有关联contig的已知簇")
            n_pure = n_known - n_linked
            if n_pure > 0:
                _add("  └ CD-HIT pure", "○", key_metric=f"{n_pure} 纯参考簇 (不进下游)")
        else:
            _add("04_CLUSTER", "○", details="未运行")

        # ── 05_Taxonomy + taxonomy_summary.tsv ──
        tax = self.d['taxonomy']
        int_dir = tax / "integrated"
        final_tax = int_dir / "final_integrated_classification.tsv"
        if final_tax.is_file():
            import csv
            n = _count_lines(final_tax) - 1
            # 1. Novelty 分级
            counts = {"Known": 0, "Novel_Species": 0, "Novel_Genus": 0, "Novel_Family": 0}
            rank_fill = {"Realm":0,"Kingdom":0,"Phylum":0,"Class":0,"Order":0,"Family":0,"Genus":0,"Species":0}
            with open(final_tax) as tf:
                for row in csv.DictReader(tf, delimiter="\t"):
                    sp = row.get("Species", row.get("species", ""))
                    ge = row.get("Genus", row.get("genus", ""))
                    fa = row.get("Family", row.get("family", ""))
                    if sp and sp not in ("NA", "-"): counts["Known"] += 1
                    elif ge and ge not in ("NA", "-"): counts["Novel_Species"] += 1
                    elif fa and fa not in ("NA", "-"): counts["Novel_Genus"] += 1
                    else: counts["Novel_Family"] += 1
                    for rk in ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]:
                        v = row.get(rk, row.get(rk.lower(), ""))
                        if v and v not in ("NA", "-"): rank_fill[rk] += 1
            _add("05_Taxonomy", "✓", key_metric=f"{n} 条, ★{counts['Known']} 已知, ★★{counts['Novel_Species']} 新种", details=f"{tax}")
            _add("  └ Novel Rank", "✓",
                 key_metric=f"Known={counts['Known']} NewSp={counts['Novel_Species']} NewGe={counts['Novel_Genus']} NewFa={counts['Novel_Family']}")
            # 2. 各 rank 填充率
            rk_info = " ".join(f"{r}={rank_fill[r]}" for r in ["Realm","Phylum","Class","Order","Family","Genus","Species"] if rank_fill.get(r,0) > 0)
            if rk_info: _add("  └ Rank Fill", "✓", key_metric=rk_info[:120])

            # 3. R consensus stats
            con_sum = int_dir / "consistency_summary.tsv"
            if con_sum.is_file():
                try:
                    cs = pl.read_csv(str(con_sum), separator="\t", null_values=["NA",""])
                    if "consistency_status" in cs.columns:
                        n_con = cs.filter(pl.col("consistency_status")=="consistent").height
                        n_inc = cs.filter(pl.col("consistency_status")=="inconsistent").height
                        n_tot = cs.height
                        _add("  └ Consensus", "✓", key_metric=f"consistent={n_con} inconsistent={n_inc} total={n_tot}")
                except: pass

            # 4. Tool agreement rate
            ag_tsv = int_dir / "agreement_stats.tsv"
            if ag_tsv.is_file():
                try:
                    ag = pl.read_csv(str(ag_tsv), separator="\t", null_values=["NA",""])
                    if "agreement_rate" in ag.columns:
                        avg_ag = ag["agreement_rate"].mean()
                        _add("  └ Agreement", "✓", key_metric=f"avg pairwise agreement={avg_ag:.1f}%")
                    elif "Agreement_Rate" in ag.columns:
                        avg_ag = ag["Agreement_Rate"].mean()
                        _add("  └ Agreement", "✓", key_metric=f"avg pairwise agreement={avg_ag:.1f}%")
                except: pass

            # 5. Resource usage per tool
            res_tsv = tax / "WVDB_votus.virus_classed" / "WVDB_votus_resource_usage.tsv"
            if res_tsv.is_file():
                try:
                    ru = pl.read_csv(str(res_tsv), separator="\t", null_values=["NA",""])
                    if "tool" in ru.columns and "wall_sec" in ru.columns:
                        total_sec = 0
                        for row in ru.iter_rows(named=True):
                            try: total_sec += float(row.get("wall_sec", 0))
                            except: pass
                        total_min = round(total_sec / 60, 1)
                        n_tools_ran = ru.height
                        _add("  └ Resources", "✓", key_metric=f"{n_tools_ran} tools, {total_min} min total")
                except: pass
        elif tax.is_dir():
            _add("05_Taxonomy", "○", key_metric="已运行但无最终结果", details=f"{tax}")
        else:
            _add("05_Taxonomy", "○", details="未运行")

        # ── 06_HostPrediction ──
        host = self.d['host_pred']
        host_summary = host / "ensemble_host_summary.tsv"
        if host_summary.is_file():
            n = _count_lines(host_summary) - 1
            try:
                hdf = pl.read_csv(str(host_summary), separator="\t", null_values=["NA", "N/A", ""])
                if "Final_Host" in hdf.columns:
                    hcounts = hdf.group_by("Final_Host").agg(pl.len()).sort("len", descending=True)
                    top3 = " | ".join(f"{r[0]}={r[1]}" for r in hcounts.head(3).iter_rows())
                    n_unknown = sum(1 for r in hcounts.iter_rows() if r[0] == "Unknown")
                    _add("06_HostPrediction", "✓", key_metric=f"{n} 条, {top3}", details=f"{host}")
                    # decision tree method stats
                    if "Decision_Method" in hdf.columns:
                        dm_counts = hdf.group_by("Decision_Method").agg(pl.len()).sort("len", descending=True)
                        top_methods = " ".join(f"{r[0]}={r[1]}" for r in dm_counts.head(4).iter_rows())
                        _add("  └ Decision", "✓", key_metric=top_methods)
                else:
                    _add("06_HostPrediction", "✓", key_metric=f"{n} 条", details=f"{host}")
            except Exception as e:
                _add("06_HostPrediction", "✓", key_metric=f"{n} 条", details=f"{host}")

        elif host.is_dir():
            _add("06_HostPrediction", "○", key_metric="已运行但无最终结果", details=f"{host}")
        else:
            _add("06_HostPrediction", "○", details="未运行")

        # ── 07_CheckV + checkv_summary.tsv (按宿主分类) ──
        cv_dir = self.d['checkv_dir']
        QUALITY_ORDER = ["Complete","High-quality","Medium-quality","Low-quality","Not-determined"]
        if cv_dir.is_dir():
            cv_tsvs = list(cv_dir.rglob("completeness.tsv"))
            if cv_tsvs:
                # 按宿主分组统计
                host_qd = {}  # {host: {quality: count}}
                for ct in cv_tsvs:
                    host_name = ct.parent.name  # 目录名即宿主名
                    if host_name not in host_qd:
                        host_qd[host_name] = dict.fromkeys(QUALITY_ORDER, 0)
                    try:
                        cv = pl.read_csv(str(ct), separator="\t", null_values=["NA", "N/A", ""])
                        comp_col = next((c for c in ["aai_completeness", "completeness"] if c in cv.columns), None)
                        if comp_col:
                            for row in cv.iter_rows(named=True):
                                val = row.get(comp_col)
                                try:
                                    v = float(val) if val and val != "NA" else None
                                    if v is None: key = "Not-determined"
                                    elif v >= 90: key = "Complete"
                                    elif v >= 50: key = "High-quality"
                                    elif v >= 10: key = "Medium-quality"
                                    else: key = "Low-quality"
                                except: key = "Not-determined"
                                host_qd[host_name][key] += 1
                    except: pass

                # 全局汇总
                global_qd = dict.fromkeys(QUALITY_ORDER, 0)
                for hqd in host_qd.values():
                    for q in QUALITY_ORDER: global_qd[q] += hqd[q]
                n_total = sum(global_qd.values())
                n_hq = global_qd["Complete"] + global_qd["High-quality"]
                n_eval = n_total - global_qd["Not-determined"]

                # 写入 checkv_summary.tsv
                with open(report_dir / "checkv_summary.tsv", "w") as cvf:
                    cvf.write("Host\t" + "\t".join(QUALITY_ORDER) + "\tTotal\tHQ\n")
                    for h in sorted(host_qd.keys()):
                        hqd = host_qd[h]; t = sum(hqd.values())
                        hq = hqd["Complete"] + hqd["High-quality"]
                        cvf.write(f"{h}\t" + "\t".join(str(hqd[q]) for q in QUALITY_ORDER) + f"\t{t}\t{hq}\n")
                    cvf.write(f"TOTAL\t" + "\t".join(str(global_qd[q]) for q in QUALITY_ORDER) + f"\t{n_total}\t{n_hq}\n")

                _add("07_CheckV", "✓", key_metric=f"{n_hq} HQ / {n_total} total ({n_eval} evaluated)", details=f"{cv_dir}")
                # 按宿主显示 HQ 数 (Top 5)
                host_hq = [(h, host_qd[h]["Complete"]+host_qd[h]["High-quality"], sum(host_qd[h].values()))
                           for h in host_qd]
                host_hq.sort(key=lambda x: -x[2])
                for h, hq, tot in host_hq[:8]:
                    if tot > 0: _add(f"  └ {h}", "✓", key_metric=f"HQ={hq} total={tot}")
                if len(host_hq) > 8: _add(f"  └ ...", "✓", key_metric=f"+{len(host_hq)-8} more hosts")
                # 质量分布摘要
                qparts = [f"{q}={global_qd[q]}" for q in QUALITY_ORDER if global_qd[q] > 0]
                _add("  └ Distribution", "✓", key_metric=" ".join(qparts[:4]))
            else:
                _add("07_CheckV", "✓", key_metric=f"已运行", details=f"{cv_dir}")
        else:
            _add("07_CheckV", "○", details="未运行")

        # ── 08_Rescue ──
        rescue = self.d['rescue_dir']
        if rescue.is_dir():
            rescue_finals = list(rescue.rglob("final_centroids.fasta"))
            n_rescued = 0
            for rf in rescue_finals:
                if "branch" not in str(rf) and "known" not in str(rf):
                    n_rescued += _count_fasta(rf)
            # 各分支统计
            branch_info = []
            for bname, blabel in [("branch_a","CheckV"),("branch_b","VSI"),("branch_c","BLASTN+VSI")]:
                for bd in rescue.rglob(bname):
                    if not bd.is_dir(): continue
                    pass_fa = bd / f"{bname}B" if bname == "branch_b" else bd / f"{bname}_pass.fasta"
                    pass_fa = bd / f"{'branchA' if bname=='branch_a' else 'branchB' if bname=='branch_b' else 'branchC'}_pass.fasta"
                    if not pass_fa.is_file():
                        pass_fa = bd / f"{'branchB' if bname=='branch_b' else 'branchC'}_pass.fasta"
                    if pass_fa.is_file():
                        bp = _count_fasta(pass_fa)
                        if bp > 0: branch_info.append(f"{blabel}={bp}")
            _add("08_Rescue", "✓", key_metric=f"{n_rescued:,} HQ vOTU ({' | '.join(branch_info) if branch_info else '0 pas'})", details=f"{rescue}")
        else:
            _add("08_Rescue", "○", details="未运行")

        # ── 可选: taxonomy Sankey 图 ──
        final_tax = tax / "integrated" / "final_integrated_classification.tsv" if tax.is_dir() else None
        if final_tax and final_tax.is_file():
            sankey_script = Path(__file__).resolve().parent.parent / "analysis" / "taxonomic_sankey.py"
            if sankey_script.is_file():
                try:
                    import importlib.util
                    if not importlib.util.find_spec("plotly"):
                        self.log.info("  安装 plotly...")
                        import subprocess as sp
                        sp.run([sys.executable, "-m", "pip", "install", "plotly", "-q"], check=False)
                    import subprocess as sp
                    sp.run([sys.executable, str(sankey_script),
                            "-i", str(final_tax),
                            "-o", str(report_dir / "classification_sankey.png"),
                            "--format", "png", "--min-flow", "1", "--min-genus-flow", "10",
                            "--palette", "set3", "--height", "1200",
                            "--node-pad", "30", "--label-truncate", "25",
                            "--font-size", "9", "--title-font-size", "16"],
                           capture_output=True, timeout=120)
                    self.log.info("  Sankey 图 → %s", report_dir / "classification_sankey.png")
                except: pass

        # ── 写入报告 ──
        # 1. stage_summary.tsv
        import csv as csv_mod
        summary_tsv = report_dir / "stage_summary.tsv"
        with open(summary_tsv, "w", newline="") as sf:
            w = csv_mod.DictWriter(sf, fieldnames=["Stage", "Status", "Key_Metric", "Details"], delimiter="\t")
            w.writeheader()
            for s in stage_stats: w.writerow(s)

        # 2. 目录树
        tree_file = report_dir / "directory_tree.txt"
        with open(tree_file, "w") as tf:
            for d in sorted(self.d['root'].iterdir()):
                if not d.is_dir(): continue
                tf.write(f"{d.name}/\n")
                for sd in sorted(d.iterdir()):
                    if sd.is_dir():
                        tf.write(f"  {sd.name}/\n")
                        for f in sorted(sd.iterdir())[:5]:
                            tf.write(f"    {f.name}\n")
                        rest = sum(1 for _ in sd.iterdir()) - 5
                        if rest > 0: tf.write(f"    ... +{rest} more\n")
                    else:
                        tf.write(f"  {sd.name}\n")

        # 3. 复制日志
        import shutil as shutil_mod
        log_src = self.d['root'] / "orchestrator.log"
        if log_src.is_file():
            shutil_mod.copy(log_src, report_dir / "orchestrator.log")

        # 4. HTML 报告 (含图表)
        _write_html_report(report_dir, stage_stats)

        # ── 屏幕输出 ──
        self.log.info("=" * 50)
        self.log.info("  流水线完成!")
        self.log.info("")
        self.log.info("  阶段汇总:")
        self.log.info("  %-22s %4s  %s", "Stage", "状态", "关键指标")
        self.log.info("  " + "-" * 70)
        for s in stage_stats:
            if s["Stage"].startswith("  "): continue
            icon = "✓" if s["Status"] == "✓" else "○"
            self.log.info("  %-22s  %s   %s", s["Stage"], icon, s["Key_Metric"])
        self.log.info("")
        self.log.info("  详细报告: %s", report_dir)
        self.log.info("    pipeline_report.html  网页版报告")
        self.log.info("    stage_summary.tsv     阶段汇总表")
        self.log.info("    directory_tree.txt    输出目录树")
        self.log.info("    orchestrator.log      运行日志")
        self.log.info("=" * 50)


# ═══════════════════════════════════════════════════════════════════
# 4. 参数解析
# ═══════════════════════════════════════════════════════════════════

STAGE_HELP = {
    'clean': """
  --stage clean — 数据清洗 (Fastp + Seqkit + Clumpify)

  子脚本: clean-data.py  (支持断点续传: .clean_checkpoints)
    Fastp    : 滑动窗口质控 (Q20, min_len=50, Poly-G 去除)
    Seqkit   : FASTQ → FASTA 格式转换 (--no-compress 控制压缩)
    Clumpify : BBMap 光学去重 (--skip-clumpify 跳过)

  编排器透传参数 (★ 必需 | · 可选):
    ★ --input_reads        输入 FASTQ 目录
    ★ --output_dir         输出根目录 (子目录: 00a_CleanData/)
    · --threads -t         fastp 线程数 (默认 20)
    · --jobs -j            并行样本数 (默认 2)
    · --skip_clumpify      跳过 Clumpify 光学去重
    · --dedup              启用 fastp 自带去重
    · --clumpify_memory    clumpify Java 堆内存 (默认: 10g)
    · --no_compress        最终结果不使用 gzip 压缩
    · --clean_debug        输出 clean-data.py 详细调试日志
    · --force              强制重跑 (清除断点记录)

  输出目录:
    00a_CleanData/
    ├── 1.fastp_tmp/       (fastp 质控后 FASTQ, 转 FASTA 后自动清理)
    ├── 2.fasta/            seqkit 转换的 FASTA
    ├── 3.clumpify/         clumpify 去重 FASTA (最终 reads 指针)
    └── logs/               每样本日志 + fastp HTML/JSON 报告
""",
    'deplete': """
  --stage deplete — 去宿主 + 去rRNA (Kraken2 + 精准比对 + rRNA剔除)

  子脚本: host_depletion.py  (支持断点续传: .checkpoints)
    Kraken2             : 物种分类标记宿主 reads
    Bowtie2/HISAT2/Minimap2 : 精准比对去宿主
    rRNA 去除           : Ribodetector (rna-short) / SILVA Bowtie2 (不限)

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir         输出根目录 (子目录: 00b_HostDepletion/)
    · --host_db            宿主数据库根目录 (自动推导 kraken2/bowtie2 子目录)
    · --kraken2_db         Kraken2 宿主库 (覆盖 --host_db 自动检测)
    · --host_align_db      宿主比对索引前缀 (覆盖 --host_db 自动检测)
    · --aligner            {bowtie2,hisat2,minimap2} (默认 bowtie2)
    · --seq_type           {dna-short,rna-short,nanopore,pacbio} (默认 rna-short)
    · --rrna               开启 rRNA 剔除
    · --rrna_tool          {ribodetector,silva} (默认 ribodetector)
    · --silva_index        SILVA Bowtie2 索引前缀 (--rrna_tool silva 时必需)
    · --kraken2_confidence Kraken2 分类置信度阈值 (默认 0.4)
    · --deplete_steps      消融实验步骤 (默认: kraken2,align,rrna)
    · --keep_rrna          保留分离的 rRNA reads 到 rrna/ 目录
    · --rrna_chunk_size    ribodetector chunk_size (默认 256)
    · --rrna_report        rRNA 统计报告文件名 (默认: ribodetector.report.txt)
    · --align_config       透传给比对工具的额外参数
    · --deplete_tmp        临时文件目录
    · --deplete_debug      输出 host_depletion.py 详细调试日志
    · -t/--threads         每样本线程数 (默认 20)
    · -j/--jobs            并行样本数 (默认 2)
    · --force              强制重跑 (清除断点记录)

  输出目录:
    00b_HostDepletion/
    ├── *_clean_1.fa.gz / *_clean_2.fa.gz  (清洁 reads, 可配对的 fasta)
    ├── rrna/               (--keep_rrna 时保留的 rRNA reads)
    ├── logs/               (每样本子进程日志)
    ├── kraken2_report/     (Kraken2 分类报告 per-sample)
    ├── host_depletion_seqkit_summary.tsv   (各阶段 reads 数量统计)
    ├── host_depletion_resource_usage.tsv   (每样本资源消耗汇总)
    └── host_depletion_plot_*.png           (reads 变化可视化图表)
""",
    'assembly': """
  --stage assembly — 宏转录组组装 (Penguin / MEGAHIT / rnaviralSPAdes)

  子脚本: assembly_pipeline.py  (支持断点续传: 检测已有 .contig.fasta)
    penguin         : 宏转录组专用组装 (guided_nuclassemble, 默认)
    megahit         : 宏基因组组装 (k-list 21..99)
    rnaviralspades  : RNA 病毒专用 SPAdes
    all             : 三工具并行; ≥2 工具自动 refineC split + merge

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir         输出根目录 (子目录: 01_Assembly/)
    · --input_reads        输入 reads 目录 (默认读取 00b_HostDepletion/)
    · --assembler          {penguin,megahit,rnaviralspades,all} (默认 penguin)
    · --contig-length      最小 contig 长度 bp (默认 200)
    · -t/--threads         单任务线程数 (默认 20)
    · -m/--memory          单任务内存 GB (默认 64)
    · -j/--jobs            并行样本数 (默认 2)
    · --refinec_threads    refineC 独立线程数 (默认使用 --threads)
    · --refinec_frag_min_len  refineC split 最小片段长度 bp (默认 1000)
    · --refinec_min_id     refineC merge 最小序列一致性 (默认 0.95)
    · --refinec_min_cov    refineC merge 最小覆盖度 (默认 0.50)
    · --asm_tmp_dir        组装临时文件目录
    · --asm_keep_temp      保留临时文件及 refineC 中间目录 (调试用)
    · --force              强制重跑 (覆盖已有 contig)

  多工具模式 (--assembler all / megahit,rnaviralspades,penguin):
    ≥2 工具时自动启用 refineC split → merge_all → refineC merge 管道
    最终产出: {sample}_all_tools_refineC_merge.merged.fasta

  输出目录:
    01_Assembly/{sample}/
    ├── {sample}_{tool}.contig.fasta          (组装 contig)
    ├── {sample}.{tool}.log                   (stdout/stderr)
    ├── {sample}.{tool}.time.mem.log          (资源使用)
    └── {sample}_all_tools_refineC_merge.merged.fasta  (多工具 merge 产出)
""",
    'identification': """
  --stage identification — 6 工具并行病毒序列鉴定

  子脚本: virus_identification16.py  (支持断点续传 + 资源监控)
    Genomad          : 深度学习病毒/质粒/前病毒分类
    Diamond BLASTX   : RefSeq 病毒蛋白 + UniProt 蛋白比对
    VirSorter2       : 隐马尔可夫模型病毒检测 (--virsorter_group)
    ViralVerify      : 病毒蛋白 HMM 验证
    VirHunter        : 深度学习病毒鉴定
    Metabuli         : k-mer 序列分类
    ─────────────────────────────────────────
    后置过滤         : UniProt 验证 + NR 库对抗过滤
    Venn/Upset 图    : 多工具交集可视化

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir         输出根目录 (子目录: 02_Identification/)
    ★ --virus_db           病毒鉴定数据库根目录
    · --input_assembly     组装结果目录 (默认自动读取 01_Assembly/)
    · --identify_tools     鉴定工具: all / genomad,diamond,metabuli,... (默认 all)
    · --blast_mode         {strict,filter,both,no-filter} (默认 filter)
    · --blast_evalue       BLAST e-value 阈值 (默认 1e-5)
    · --blast_top_n        BLAST 对抗验证 Top N (默认 5)
    · --virsorter_group    VirSorter2 分组 (默认 dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae)
    · --ident_ext          输入文件扩展名 (默认 .fasta)
    · -t/--threads         单样本线程数 (默认 20)
    · -j/--jobs            并行样本数 (默认 2)
    · --force              强制重跑

  鉴定数据库 (全部可选, 自动从 --virus_db 推导):
    · --virus_protein_db  病毒蛋白 Diamond DB
    · --uniprot_db        UniProt Diamond DB
    · --nr_db             NR Diamond DB (对抗过滤)
    · --viroids_db        类病毒 BLAST DB
    · --virsorter_db      VirSorter2 数据库
    · --viralverify_hmm   ViralVerify HMM 文件
    · --metabuli_db       Metabuli 数据库
    · --virus_taxid       病毒 TaxID 列表
    · --virhunter_path    VirHunter predict_cpu.py 路径
    · --virhunter_weights  VirHunter weights 目录
    · --virbot_path       VirBot.py 路径
    · --viralm_path       viralm_cpu.py 路径

  后置过滤控制:
    · --skip_uniprot_filter  跳过 UniProt 后置过滤
    · --skip_nr_filter       跳过 NR 后置过滤
    · --skip_id_plots        跳过 Venn/Upset 图表生成
    · --clean_failed         自动清理失败任务目录

  输出目录:
    02_Identification/{sample}/
    ├── *_virus.all.candidate.fasta         (最终候选病毒)
    ├── uniprot_filter_output_strict/       (UniProt 严格过滤结果)
    ├── uniprot_filter_output_filter/       (UniProt 宽松过滤结果)
    ├── *_resource.tsv                      (资源消耗)
    └── Venn/                               (Venn/Upset 图表)
""",
    'cobra': """
  --stage cobra — COBRA 批量延伸 (BWA-MEM2 + CoverM + COBRA)

  子脚本: cobra_pipeline.py  (支持 JSON 断点续传: checkpoint_status.json)
    流程: BWA-MEM2 比对 → Samtools sort → CoverM 覆盖度 → COBRA 重叠延伸
    自动匹配: reads + contig + virus 三元组
    病毒来源: 支持 raw/filter/strict 三种模式 (--virus_mode)

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir         输出根目录 (子目录: 03_COBRA/)
    · --input_reads        输入 reads 目录 (默认读取 00b_HostDepletion/)
    · --virus_mode         {raw,filter,strict} 病毒序列来源 (默认 strict)
    · --cobra_mink         COBRA 最小 kmer (默认 21)
    · --cobra_maxk         COBRA 最大 kmer (默认 141)
    · --cobra_linkage_mismatch  链接识别不匹配数 (默认 2)
    · --cobra_verbose      输出 cobra_pipeline.py 详细日志
    · -t/--threads         单任务线程数 (默认 20)
    · -j/--jobs            并行任务数 (COBRA 内存密集, jobs 自动减半)
    · --force              禁用断点续传, 强制重跑所有任务

  输出目录:
    03_COBRA/{sample}/
    └── cobra_{tool}_result/
        ├── {sample}.{mode}.{tool}.cobra.fa    (COBRA 延伸结果)
        ├── {sample}.{mode}.{tool}.COBRA/       (COBRA 原始输出)
        └── {sample}.{mode}.{tool}.log          (任务日志)
""",
    'cluster': """
  --stage cluster — 聚类去冗余 (CD-HIT 参考引导 + vclust Leiden)

  子脚本: cluster_pipeline.py  (支持 --resume 断点续传)
    Step 1  seqkit        : 最小长度过滤 (--min-length, 默认 500bp)
    Step 2a CD-HIT 参考引导 : 合并 contig + ICTV/NCBI 参考 → vclust cd-hit
                             → 拆分 known/novel → 产出 association 映射
    Step 2b vclust Leiden  : 仅 novel contig 聚类 (prefilter→align→cluster)
    Step 3  输出            : centroids + per-cluster 拆分 + 统计

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir         输出根目录 (子目录: 04_CLUSTER/)
    · --cluster_input      直接输入已合并 FASTA (跳过 COBRA 收集)
    · --ref-genomes        ICTV/NCBI 参考基因组 FASTA (可多个, CD-HIT 预聚类)
    · --min-length         病毒最小长度 bp (默认 500)
    · --ani                vclust Leiden ANI 阈值 (默认 0.95)
    · --qcov               vclust Leiden QCOV 阈值 (默认 0.85)
    · --cdhit_ani          CD-HIT ANI 阈值 (默认 0.95, 转录组建议 0.85)
    · --cdhit_qcov         CD-HIT QCOV 阈值 (默认 0.85, 转录组建议 0.50)
    · --skip_vclust        跳过 vclust Leiden (仅做 CD-HIT 预聚类)
    · --vclust_cluster_file  复用已有 vclust 聚类 TSV (跳过聚类计算)
    · -t/--threads         线程数 (默认 20)
    · --force              禁用断点续传, 强制重跑

  输出目录:
    04_CLUSTER/
    ├── 1_seqkit/virus.candidate.fasta          (长度过滤后)
    ├── 2_cdhit/
    │   ├── known_centroids.fasta               (已知簇 centroids)
    │   ├── known_association.tsv               (contig→参考映射)
    │   ├── novel_contigs.fasta                 (新颖 contig)
    │   └── known_clusters/                     (per-cluster 已知簇拆分)
    ├── 3_vclust/
    │   ├── vclust_clusters.tsv                 (聚类结果)
    │   ├── cluster_summary.tsv                  (Polars 统计)
    │   └── split_fastas/                       (per-cluster novel 拆分)
    └── centroids/
        ├── final_centroids.fasta               (全部代表序列, 供后续 taxonomy/host)
        ├── known_association.tsv
        └── known_ids.txt
""",
    'taxonomy': """
  --stage taxonomy — 9 工具病毒分类 + R 共识整合

  子脚本: virus_classifier2.py + virus_classifier_analysis14.R
    genomad         : 深度学习全基因组病毒分类
    metabuli        : k-mer 序列分类 + taxonkit lineage
    CAT             : BAT/CAT 蛋白比对分类
    diamond_lca     : Diamond BLASTX LCA 分类
    mmseqs          : 蛋白序列比对分类
    VITAP           : 病毒蛋白分类
    ACVirus         : 古菌病毒分类
    vcontact3       : 蛋白簇网络分类 (基因共享网络)
    PhaGCN3         : 噬菌体 GCN 分类
    R consensus     : 多工具投票共识 (vcontact3>vitap>acvirus>mmseqs>genomad)
                     → 8 级标准 taxonomy (Realm..Species)

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir         输出根目录 (子目录: 05_Taxonomy/)
    ★ --virus_db           病毒分类数据库根目录
    · -t/--threads         线程数 (默认 20)
    · --tax_tools          分类工具, 逗号分隔 (默认: all)
    · --tax_jobs           分类并行任务数 (默认 1)
    · --tax_ext            分类输入扩展名 (默认 .fasta)
    · --tax_remove_suffix  分类输入去后缀名
    · --force              强制重跑

  分类数据库 (全部可选):
    · --genomad_db        genomad 数据库
    · --metabuli_db       Metabuli 数据库
    · --cat_db            CAT 数据库
    · --cat_tax           CAT taxonomy 路径
    · --uniprot_db        UniProt Diamond DB (diamond_lca)
    · --mmseqs_db         MMseqs2 数据库
    · --vitap_db          VITAP 数据库
    · --acvirus_db        ACVirus 数据库
    · --vcontact3_db      vConTACT3 数据库

  输出目录:
    05_Taxonomy/
    ├── WVDB_votus.virus_classed/
    │   └── WVDB_votus_combined_taxonomy.tsv    (9 工具合并分类)
    └── integrated/
        └── final_integrated_classification.tsv  (R 共识最终分类)
""",
    'host': """
  --stage host — 宿主预测 (决策树: ICTV > RNAVirHost > PhaBOX2)

  子脚本: run_host_prediction.py --mode all
    ICTV (C9)       : 官方权威分类库查找宿主 (--prob-dir, 级联: Species→Genus→Family)
    RNAVirHost      : 全生态位模型预测 (plant/animal/fungi/bacteria, 两步法)
    PhaBOX2 CHERRY  : CRISPR+AAI 网络噬菌体宿主预测 (--phabox-db)
    决策树          : Class 硬规则 (Caudoviricetes→Bacteria)
                      ICTV==RVH→直接采用, 分歧时 PB2 决胜, 全分歧→ICT>RVH>PB2

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir         输出根目录 (子目录: 06_HostPrediction/)
    · --phabox-db          PhaBOX2 数据库路径
    · --prob-dir           ICTV 宿主概率表目录 (默认: cross_analysis/)
    · -t/--threads         线程数 (默认 20)
    · --force              强制重跑
    · --skip_rnavirhost    跳过 RNAVirHost 宿主预测
    · --skip_phabox        跳过 PhaBOX2 宿主预测
    · --skip_ictv          跳过 ICTV 宿主查找

  输出目录:
    06_HostPrediction/
    ├── ensemble_host_summary.tsv                    (决策树最终宿主)
    ├── host_classified_fasta/{host}.classified.fasta
    ├── RVH_result/                                 (RNAVirHost 原始结果)
    ├── phabox2_output/                             (PhaBOX2 原始结果)
    └── C9_ICTV_result/                             (ICTV 查找结果)
""",
    'checkv': """
  --stage checkv — 按宿主分类 CheckV 完整性预评估

  功能: 对每个宿主类别的 centroids 分别运行 checkv completeness
        标记 completeness ≥ 90% 的 centroids 为【免拯救】
        (与 CD-HIT known 合并, 在三支路拯救前直接输出)

  调用命令: checkv completeness <fasta> <outdir> -d <db> -t N

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir         输出根目录 (子目录: 07_Checkv/)
    ★ --checkv_db          CheckV 数据库路径
    · -t/--threads         线程数 (默认 20, 自动上限 16)
    · --force              强制重跑

  质量等级推导 (新版 CheckV 无 quality 列时从 completeness 推导):
    Complete       : completeness ≥ 90%
    High-quality   : 50% ≤ completeness < 90%
    Medium-quality : 10% ≤ completeness < 50%
    Low-quality    : completeness < 10%
    Not-determined : 无法判断

  输出目录:
    07_Checkv/
    ├── {host}/completeness.tsv         (各宿主 CheckV 结果)
    ├── {host}.fasta                    (各宿主 centroids)
    └── checkv_pass_ids.txt            (≥90% pass centroids, 供 rescue 阶段)
""",
    'rescue': """
  --stage rescue — 宿主过滤 + 三支路级联拯救 + CheckV 质量报告

  子脚本: rescue_pipeline.py  (支持 --resume 断点续传)
    前置步骤 : 按 Final_Host 过滤 centroids (--host-filter)
               CD-HIT known + CheckV pass(≥90%) → 免拯救, 直接输出
    分支 A    : CheckV 并行评估 centroids (分块, completeness ≥ 90% pass)
    分支 C    : Virseqimprover reads 迭代延伸 (cluster 内多样本 reads 聚合)
                Salmon 定量 → BBMap 提取 → SPAdes 组装 → CheckV 验证
    分支 D    : BLASTN megablast + CheckV + VSI 最后拯救
    合并      : A+C+D pass → vclust 最终去重 → HQ vOTU

  免拯救 = CD-HIT known (参考关联) + CheckV pass (completeness ≥ 90%)

  编排器透传参数 (★ 必需 | · 可选):
    ★ --output_dir             输出根目录 (子目录: 08_Rescue/)
    ★ --checkv_db              CheckV 数据库路径
    · --input_reads            输入 reads 目录 (默认读取 00b_HostDepletion/)
    · --blast-db               BLASTN 参考数据库 (分支 D, ref.fasta)
    · --host-filter            目标宿主, 逗号分隔 (默认 Plant)
    · --virseqimprover-path    Virseqimprover.py 路径 (默认同目录)
    · --salmon-bin             Salmon 二进制路径 (默认: salmon)
    · --min-length             病毒最小长度 bp (默认 500)
    · --ani                    最终 vclust ANI (默认 0.95)
    · --qcov                   最终 vclust QCOV (默认 0.85)
    · -t/--threads             线程数 (默认 20)
    · -j/--jobs                并行任务数 (默认 2)
    · --force                  禁用断点续传

  输出目录:
    08_Rescue/
    ├── known/centroids/final_centroids.fasta     (免拯救 known centroids)
    ├── {host}/
    │   ├── input_centroids.fasta                 (目标宿主 centroids)
    │   ├── branch_a/branchA_pass.fasta           (CheckV pass)
    │   │           branchA_fail.fasta             (CheckV fail)
    │   ├── branch_b/branchB_pass.fasta           (VSI pass)
    │   │           branchB_fail.fasta             (VSI fail)
    │   ├── branch_c/branchC_pass.fasta           (BLASTN+VSI pass)
    │   ├── merged/all_HQ.fasta                   (A+C+D 合并)
    │   └── centroids/final_centroids.fasta       (最终 HQ vOTU ★)
    └── checkv/
        ├── {host}/completeness.tsv               (各宿主 CheckV)
        └── all/completeness.tsv                  (全部 centroids)
""",
}

OVERVIEW = """
═══════════════════════════════════════════════════════════════
  Virome Pipeline v2.3 — 宏病毒组端到端全自动主控
═══════════════════════════════════════════════════════════════

  10 个独立阶段 (--stage):

    clean           清洗 (Fastp + Seqkit + Clumpify)
    deplete         去宿主 (Kraken2 + Bowtie2 + Ribodetector)
    assembly        宏转录组组装 (penguin)
    identification  病毒序列鉴定 (6 工具)
    cobra           COBRA 批量延伸
    cluster         vclust 聚类 (仅聚类, 不拯救)
    taxonomy        分类注释 (5 工具 + R 共识)
    host            宿主预测 (RNAVirHost + PhaBOX2 + ICTV)
    rescue          宿主过滤 + 三支路级联拯救 [NEW]
    all             全流程串行

  新流程 (v2.3):
    cluster → taxonomy → host → rescue(按宿主过滤)
    Unknown 宿主默认跳过, 输出到 unknown_votus.fasta

  示例:
    python virome_pipeline.py --stage all \\
      --input_reads /data/host/ --output_dir /data/out/ \\
      --virus_db /db/virus/ --checkv_db /db/checkv/ \\
      --host-filter Plant,Animal -t 120 -j 20

  查看阶段详情:
    python virome_pipeline.py --stage rescue --help

  查看所有参数:
    python virome_pipeline.py --help
═══════════════════════════════════════════════════════════════
"""


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
    g.add_argument('--salmon-bin', default='salmon', help='Salmon 二进制路径 (默认: salmon)')

    g = p.add_argument_group('分类数据库 (virus_classifier2.py)')
    g.add_argument('--genomad_db', help='genomad DB 路径')
    g.add_argument('--cat_db', help='CAT 数据库路径')
    g.add_argument('--cat_tax', help='CAT taxonomy 路径')
    g.add_argument('--mmseqs_db', help='mmseqs 数据库路径')
    g.add_argument('--vitap_db', help='VITAP 数据库路径')
    g.add_argument('--acvirus_db', help='ACVirus 数据库路径')
    g.add_argument('--vcontact3_db', help='vConTACT3 数据库路径')
    return p


def _write_html_report(report_dir, stage_stats):
    """生成自包含 HTML 流水线报告 (含 Chart.js 图表)"""
    from datetime import datetime
    import re, json as _json

    main_stages = [s for s in stage_stats if not s["Stage"].startswith("  ")]
    sub_stages  = [s for s in stage_stats if s["Stage"].startswith("  ")]

    status_icon = {"✓": "pass", "○": "skip", "✗": "fail"}
    def _esc(v): return str(v).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

    # ── 从 stage_stats 提取图表数据 ──
    chart_scripts = ""
    n_pass = sum(1 for s in main_stages if s["Status"] == "✓")
    n_total = sum(1 for s in main_stages if s["Status"] != "○")
    pct = round(n_pass/max(n_total,1)*100)

    # Helper: 从 Key_Metric 或 Details 提取 "key=value" 对
    def _extract_kv(text, pattern):
        result = {}
        for m in re.finditer(pattern, text):
            result[m.group(1)] = int(m.group(2))
        return result

    # 1. 鉴定工具图 — 读 ident_summary.tsv
    id_tsv = report_dir / "ident_summary.tsv"
    if id_tsv.is_file():
        with open(id_tsv) as f:
            hdr = f.readline().strip().split('\t')
            line = f.readline()  # first sample
            if line:
                parts = line.strip().split('\t')
                tools = hdr[2:]  # skip Sample, All_Candidate
                vals = {}
                for i, t in enumerate(tools):
                    try: vals[t] = int(parts[i+2])
                    except: pass
                if vals:
                    chart_scripts += f"""
    new Chart(document.getElementById('chart_ident'), {{
      type:'bar', data:{{labels:{_json.dumps(list(vals.keys()))},
      datasets:[{{label:'Virus contigs',data:{_json.dumps(list(vals.values()))},
      backgroundColor:'#42a5f5'}}]}},
      options:{{plugins:{{title:{{display:true,text:'Per-Tool Identification'}}}}}}}});
"""

    # 2. 组装 N50 图 — 读 assembly_summary.tsv
    asm_tsv = report_dir / "assembly_summary.tsv"
    if asm_tsv.is_file():
        samples, n50s, contigs = [], [], []
        with open(asm_tsv) as f:
            f.readline()
            for line in f:
                p = line.strip().split('\t')
                if p[0] == 'TOTAL': continue
                samples.append(p[0][:15]); n50s.append(int(p[5])); contigs.append(int(p[3]))
        if samples:
            chart_scripts += f"""
    new Chart(document.getElementById('chart_asm'), {{
      type:'bar', data:{{labels:{_json.dumps(samples)},
      datasets:[{{label:'N50 (bp)',data:{_json.dumps(n50s)},backgroundColor:'#66bb6a',yAxisID:'y'}},
                {{label:'Contigs',data:{_json.dumps(contigs)},backgroundColor:'#ffa726',yAxisID:'y1'}}]}},
      options:{{plugins:{{title:{{display:true,text:'Assembly N50 & Contigs'}}}},
        scales:{{y:{{beginAtZero:true,position:'left'}},y1:{{beginAtZero:true,position:'right',grid:{{drawOnChartArea:false}}}}}}}}}});
"""

    # 3. 过滤图 — 读 filter_summary.tsv
    fil_tsv = report_dir / "filter_summary.tsv"
    if fil_tsv.is_file():
        flabels, fvals = [], []
        with open(fil_tsv) as f:
            f.readline()
            for line in f:
                p = line.strip().split('\t')
                if p[0] == 'TOTAL': continue
                flabels.append(f"{p[0][:10]}-{p[1]}"); fvals.append(int(p[3]))
        if flabels:
            chart_scripts += f"""
    new Chart(document.getElementById('chart_filter'), {{
      type:'bar', data:{{labels:{_json.dumps(flabels)},
      datasets:[{{label:'Passed UniProt filter',data:{_json.dumps(fvals)},
      backgroundColor:{_json.dumps(['#ef5350','#42a5f5','#66bb6a'][:len(flabels)])}}}]}},
      options:{{plugins:{{title:{{display:true,text:'UniProt Filter Passed'}}}}}}}});
"""

    # 4. Host 分布 — 从 stage_stats
    for s in stage_stats:
        if '06_HostPrediction' in s['Stage'] and '=' in s.get('Key_Metric',''):
            host_kv = _extract_kv(s['Key_Metric'], r'(\w+)=(\d+)')
            if host_kv:
                chart_scripts += f"""
    new Chart(document.getElementById('chart_host'), {{
      type:'doughnut', data:{{labels:{_json.dumps(list(host_kv.keys()))},
      datasets:[{{data:{_json.dumps(list(host_kv.values()))},
      backgroundColor:['#42a5f5','#66bb6a','#ffa726','#ef5350','#ab47bc','#26c6da','#7e57c2','#78909c']}}]}},
      options:{{plugins:{{title:{{display:true,text:'Host Prediction Distribution'}}}}}}}});
"""

    # 5. Taxonomy novelty — 从 stage_stats
    for s in stage_stats:
        if 'Novel Rank' in s['Stage'] and 'Known=' in s.get('Key_Metric',''):
            tax_kv = _extract_kv(s['Key_Metric'], r'(Known|NewSp|NewGe|NewFa)=(\d+)')
            if tax_kv:
                chart_scripts += f"""
    new Chart(document.getElementById('chart_tax'), {{
      type:'doughnut', data:{{labels:{_json.dumps(list(tax_kv.keys()))},
      datasets:[{{data:{_json.dumps(list(tax_kv.values()))},
      backgroundColor:['#66bb6a','#42a5f5','#ffa726','#ef5350']}}]}},
      options:{{plugins:{{title:{{display:true,text:'Taxonomy Classification'}}}}}}}});
"""

    # 6. CheckV quality
    checkv_kv = {}
    for s in sub_stages:
        if s['Stage'].strip().startswith('└') and s.get('Key_Metric','').endswith('条'):
            try:
                q = s['Stage'].strip().replace('└ ','')
                v = int(s['Key_Metric'].replace(' 条',''))
                checkv_kv[q] = v
            except: pass
    if checkv_kv and len(checkv_kv) > 1:
        chart_scripts += f"""
    new Chart(document.getElementById('chart_checkv'), {{
      type:'bar', data:{{labels:{_json.dumps(list(checkv_kv.keys()))},
      datasets:[{{label:'Sequences',data:{_json.dumps(list(checkv_kv.values()))},
      backgroundColor:['#66bb6a','#42a5f5','#ffa726','#ef5350','#bdbdbd']}}]}},
      options:{{plugins:{{title:{{display:true,text:'CheckV Quality Distribution'}}}}}}}});
"""

    # 图表 HTML 区域 (仅在有数据时)
    chart_divs = ""
    chart_ids = []
    if id_tsv.is_file(): chart_ids.append('ident')
    if asm_tsv.is_file(): chart_ids.append('asm')
    if fil_tsv.is_file(): chart_ids.append('filter')
    for s in stage_stats:
        if '06_HostPrediction' in s['Stage'] and '=' in s.get('Key_Metric',''): chart_ids.append('host')
        if 'Novel Rank' in s['Stage'] and 'Known=' in s.get('Key_Metric',''): chart_ids.append('tax')
    if checkv_kv and len(checkv_kv) > 1: chart_ids.append('checkv')

    for cid in chart_ids:
        chart_divs += f'<div style="background:#fff;border-radius:10px;padding:20px;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:16px"><canvas id="chart_{cid}" style="max-height:350px"></canvas></div>\n'

    rows_html = ""
    for i, s in enumerate(main_stages):
        cls = status_icon.get(s["Status"], "skip")
        rows_html += f"""<tr class="{cls}"><td style="text-align:center;color:#666;font-size:12px">{i}</td>
          <td><b>{_esc(s['Stage'])}</b></td><td style="text-align:center"><span class="badge badge-{cls}">{s['Status']}</span></td>
          <td>{_esc(s['Key_Metric'])}</td><td style="font-size:12px;color:#888">{_esc(s['Details'])}</td></tr>"""
    sub_rows = ""
    for s in sub_stages:
        sub_rows += f"""<tr><td></td><td style="padding-left:32px;color:#555">{_esc(s['Stage'].strip())}</td>
        <td></td><td>{_esc(s['Key_Metric'])}</td><td style="font-size:12px;color:#888">{_esc(s['Details'])}</td></tr>"""

    html = f"""<!DOCTYPE html>
<html lang="zh"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>MMPV-RNA Pipeline Report</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4"></script>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#f0f2f5;color:#333;line-height:1.6}}
.container{{max-width:1100px;margin:0 auto;padding:20px}}
.header{{background:linear-gradient(135deg,#1a237e,#283593);color:#fff;padding:32px;border-radius:12px;margin-bottom:24px}}
.header h1{{font-size:24px;margin-bottom:8px}}.header p{{opacity:.85;font-size:14px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:16px;margin-bottom:24px}}
.card{{background:#fff;border-radius:10px;padding:20px;box-shadow:0 1px 3px rgba(0,0,0,.08)}}
.card .value{{font-size:28px;font-weight:700;color:#1a237e}}.card .label{{font-size:13px;color:#888;margin-top:4px}}
.charts{{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-bottom:24px}}
@media(max-width:768px){{.charts{{grid-template-columns:1fr}}}}.chart-full{{grid-column:1/-1}}
.table-wrap{{background:#fff;border-radius:10px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:24px}}
table{{width:100%;border-collapse:collapse}}th{{background:#f5f5f5;text-align:left;padding:12px 16px;font-size:13px;color:#666;border-bottom:2px solid #e0e0e0}}
td{{padding:10px 16px;font-size:14px;border-bottom:1px solid #f0f0f0}}
tr.pass td:first-child{{border-left:3px solid #4caf50}}tr.fail td:first-child{{border-left:3px solid #f44336}}tr.skip td:first-child{{border-left:3px solid #ccc}}
.badge{{display:inline-block;padding:2px 10px;border-radius:12px;font-size:12px;font-weight:600}}
.badge-pass{{background:#e8f5e9;color:#2e7d32}}.badge-fail{{background:#fce4ec;color:#c62828}}.badge-skip{{background:#eee;color:#888}}
.footer{{text-align:center;color:#aaa;font-size:12px;padding:20px}}
</style></head><body><div class="container">
<div class="header"><h1>MMPV-RNA v2.3 &mdash; Pipeline Report</h1>
<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} &nbsp;|&nbsp; {n_pass}/{n_total} stages passed ({pct}%)</p></div>
<div class="cards">
<div class="card"><div class="value">{n_pass}/{n_total}</div><div class="label">Stages Passed</div></div>
<div class="card"><div class="value">{len(stage_stats)}</div><div class="label">Total Metrics</div></div>
<div class="card"><div class="value">{pct}%</div><div class="label">Success Rate</div></div>
</div>
<div class="charts">{chart_divs}</div>
<div class="table-wrap"><table><thead><tr><th style="width:40px">#</th><th>Stage</th><th style="width:60px">Status</th><th>Key Metrics</th><th style="width:200px">Details</th></tr></thead><tbody>{rows_html}{sub_rows}</tbody></table></div>
<div class="footer">MMPV-RNA v2.3 &mdash; Generated by virome_pipeline.py</div>
</div>
<script>{chart_scripts}</script>
</body></html>"""

    with open(report_dir / "pipeline_report.html", "w", encoding="utf-8") as hf:
        hf.write(html)

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
    stage_log_dir = Path(args.output_dir) / "09_Virome_Report" / "logs"
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

    pipe.run_reports()


if __name__ == '__main__':
    main()
