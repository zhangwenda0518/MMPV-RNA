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

def setup_logger(output_dir):
    """配置双通道日志: 控制台 INFO + 文件 DEBUG"""
    logger = logging.getLogger("ViromeOrch")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(output_dir, exist_ok=True)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
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
        self.orig_samples = scan_samples_in_dir(self.reads_dir)
        if not self.orig_samples:
            logger.error("在 %s 中未找到序列文件!", self.reads_dir)
            sys.exit(1)

        pe = sum(1 for v in self.orig_samples.values() if v['r2'])
        se = len(self.orig_samples) - pe
        logger.info("检测到 %d 个样本 (PE=%d, SE=%d)", len(self.orig_samples), pe, se)

    def _validate(self):
        for name, path in self.sc.items():
            if not path.exists():
                self.log.error("致命: 找不到 %s", path)
                sys.exit(1)

        stage = self.args.stage
        doing_clean = stage in ('all', 'clean')
        need_virus_db = stage in ('all', 'assembly', 'identification', 'cobra', 'taxonomy', 'host')
        need_checkv_db = stage in ('all', 'checkv', 'rescue')

        # deplete 阶段需要 kraken2 + host_align (或 --host_db 自动检测)
        need_deplete_db = stage in ('all', 'deplete') or (stage == 'clean' and not self.args.skip_depletion)
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
        self.log.info("环境验证通过 (stage=%s)", stage)

    # ── Step 0a: 数据清洗 ──
    def run_clean(self):
        self.d['clean'].mkdir(parents=True, exist_ok=True)
        if self.args.skip_clean:
            self.log.info("[0a] 跳过 (--skip_clean)")
            return

        self.log.info("=" * 50)
        self.log.info("[0a] Fastp → Seqkit → Clumpify")

        # clean-data.py: --input, --output, --fastp-threads, --jobs,
        #                 --skip-clumpify, --force, --dedup, --clumpify-memory
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

        # host_depletion.py: --tool, --seq-type, --kraken2_index, --step2_index,
        #   --input-dir, --outdir, --jobs, --threads, --logs_dir, --rrna,
        #   --filter, --confidence
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

        # 收集所有 contig 文件
        all_contigs = []
        for sample_name, tools in asm_map.items():
            for tool, contig_path in tools.items():
                all_contigs.append(contig_path)

        self.log.info("  %d 个 contig 文件待鉴定 (内部 --jobs %d 并行)",
                     len(all_contigs), self.args.jobs)

        # virus_identification16.py: --input (文件或目录), --output, --db_dir,
        #   --identify_tools, --threads, --jobs, --force
        for contig in all_contigs:
            parts = [
                f"python {self.sc['identify']}",
                f"--input {contig}",
                f"--output {self.d['ident']}",
                f"--db_dir {self.args.virus_db}",
                f"--identify_tools {self.args.identify_tools}",
                f"--threads {self.args.threads}",
                f"--jobs {self.args.jobs}",
                f"--blast_mode {self.args.blast_mode}",
                f"--blast_evalue {self.args.blast_evalue}",
                f"--blast_top_n {self.args.blast_top_n}",
                f"--virsorter_group {self.args.virsorter_group}",
            ]
            for arg in ['virus_protein_db', 'uniprot_db', 'viroids_db', 'virsorter_db',
                         'viralverify_hmm', 'metabuli_db', 'virus_taxid',
                         'virhunter_path', 'virhunter_weights', 'virbot_path', 'viralm_path']:
                val = getattr(self.args, arg, None)
                if val:
                    parts.append(f"--{arg} {val}")
            if self.args.force:
                parts.append("--force")
            run_cmd(' '.join(parts), self.log, f"VirusID: {contig.name}")

        self.viral_map = scan_viral_files(self.d['ident'])
        self.log.info("  鉴定完成: %d 样本有病毒候选序列", len(self.viral_map))

    # ── Step 3: COBRA 延伸 ──
    def run_cobra(self):
        self.d['cobra'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[3] COBRA 批量延伸 (BWA-MEM2 + COBRA + CheckV)")

        # cobra_pipeline.py 自动从 reads 目录发现样本，
        # 从 contigs 目录查找组装结果，从 virsorter 目录查找病毒序列。
        # CLI:
        #   --mode, --reads-dir, --contigs-dir, --virsorter-dir, --output-dir,
        #   --assembly-tools, --checkv-db, --checkv-mode, --jobs, --threads,
        #   --mink, --maxk, --linkage-mismatch
        asm_tools = (['megahit', 'rnaviralspades', 'penguin']
                     if self.args.assembler == 'all'
                     else [self.args.assembler])

        parts = [
            f"python {self.sc['cobra']}",
            f"--mode mix",
            f"--reads-dir {self.reads_dir}",
            f"--contigs-dir {self.d['asm']}",
            f"--virsorter-dir {self.d['ident']}",
            f"--output-dir {self.d['cobra']}",
            f"--assembly-tools {','.join(asm_tools)}",
            f"--jobs {max(1, self.args.jobs // 2)}",
            f"--threads {self.args.threads}",
        ]

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

        # 优先使用直接输入, 否则自动收集 COBRA 结果
        if self.args.cluster_input and os.path.isfile(self.args.cluster_input):
            cluster_fa = Path(self.args.cluster_input)
            self.log.info("  CLUSTER 直接输入: %s", cluster_fa)
        else:
            cluster_input = self.d['cobra']
            if not Path(cluster_input).exists():
                self.log.error("COBRA 输出 %s 不存在, 跳过 CLUSTER", cluster_input)
                return
            cluster_fa = Path(self.d['root']) / "cluster_input.fasta"
            with open(cluster_fa, 'w') as out:
                for f in Path(cluster_input).rglob('*.cobra.fa'):
                    with open(f) as inf:
                        out.write(inf.read())
            if cluster_fa.stat().st_size == 0:
                self.log.warning("未找到任何输入, 跳过 CLUSTER")
                return
            self.log.info("  CLUSTER 输入 (自动收集): %s", cluster_fa)

        self.log.info("  CLUSTER 输入: %s", cluster_fa)

        parts = [
            f"python {self.sc['cluster']}",
            f"-i {cluster_fa}",
            f"-o {self.d['cluster']}",
            f"-t {self.args.threads}",
            f"--min-length {self.args.min_length}",
            f"--ani {self.args.ani}",
            f"--qcov {self.args.qcov}",
            f"--stop-after-vclust",
        ]
        if self.args.ref_genomes:
            parts.append(f"--ref-genomes {' '.join(self.args.ref_genomes)}")

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
                f"--remove_tmp"
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
        checkv_col = None
        for c in ["checkv_quality", "quality"]:
            if c in df.columns:
                checkv_col = c
                break

        if checkv_col is None:
            self.log.warning("  CheckV 输出缺少 quality 列")
            return {}, total

        dist = {}
        for val in df[checkv_col].to_list():
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

        # Step 1: virus_classifier2.py
        if not self.args.force and combined_tsv.is_file() and combined_tsv.stat().st_size > 100:
            self.log.info("  [SKIP] virus_classifier2 — 已有结果")
        else:
            parts = [
                f"python {self.sc['classifier']}",
                f"-g {centroids}",
                f"-s {sample}",
                f"-t genomad,mmseqs,VITAP,ACVirus,vcontact3",
                f"-o {tax_dir}",
                f"-p {self.args.threads}",
                f"--db-dir {self.args.virus_db}",
            ]
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
                       f"-d {self.args.checkv_db} -t {min(self.args.threads, 16)} --remove_tmp")
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
    def run_reports(self):
        self.log.info("=" * 50)
        self.log.info("[Reports] 流水线汇总")

        cobra_dir = self.d['cobra']

        # 统计 COBRA 输出 (仅当目录存在)
        sample_count = 0
        cobra_file_count = 0
        if cobra_dir.is_dir():
            for sample_dir in cobra_dir.iterdir():
                if not sample_dir.is_dir():
                    continue
                sample_count += 1
                for result_dir in sample_dir.iterdir():
                    if result_dir.is_dir() and result_dir.name.startswith('cobra_'):
                        for cf in result_dir.glob('*.cobra.fa'):
                            if cf.stat().st_size > 0:
                                cobra_file_count += 1

        self.log.info("=" * 50)
        self.log.info("  流水线完成!")
        self.log.info("")
        self.log.info("  输出目录结构:")
        self.log.info("    数据清洗:     %s", self.d['clean'] if self.d['clean'].is_dir() else '(未运行)')
        self.log.info("    去宿主:       %s", self.d['hostdep'] if self.d['hostdep'].is_dir() else '(未运行)')
        self.log.info("    组装结果:     %s", self.d['asm'] if self.d['asm'].is_dir() else '(未运行)')
        self.log.info("    病毒鉴定:     %s", self.d['ident'] if self.d['ident'].is_dir() else '(未运行)')
        if cobra_file_count > 0:
            self.log.info("    COBRA 延伸:   %s (%d 样本, %d 个延伸结果)",
                         cobra_dir, sample_count, cobra_file_count)
        self.log.info("    CLUSTER 去冗余:  %s", self.d['cluster'] if self.d['cluster'].is_dir() else '(未运行)')
        self.log.info("    运行日志:     %s", self.d['root'] / 'orchestrator.log')
        self.log.info("=" * 50)


# ═══════════════════════════════════════════════════════════════════
# 4. 参数解析
# ═══════════════════════════════════════════════════════════════════

STAGE_HELP = {
    'clean': """
  --stage clean — 数据清洗 (Fastp + Seqkit + Clumpify)

  子脚本: clean-data.py
    Fastp    滑动窗口质控 (--fastp-threads)
    Seqkit   统计 reads 数和碱基数
    Clumpify BBMap 光学去重 (--skip-clumpify 跳过)

  调用命令: python clean-data.py --input <dir> --output <dir> --fastp-threads N --jobs N

  参数:
    --input_reads    输入 FASTQ 目录 [必需]
    --output_dir     输出根目录 [必需]
    --skip_clumpify  跳过 Clumpify 光学去重
    -t, --threads    线程数 (默认 20)
    -j, --jobs       并行数 (默认 2)
    --force          强制重跑

  输出: 00a_CleanData/
""",
    'deplete': """
  --stage deplete — 去宿主 + 去 rRNA (Kraken2 + Bowtie2/HISAT2 + rRNA)

  子脚本: host_depletion.py
    Kraken2         物种分类标记宿主 reads (--confidence 置信度)
    Bowtie2/HISAT2  精确比对去宿主 (--tool, --seq-type)
    rRNA 去除       Ribodetector (默认) 或 SILVA Bowtie2

  调用命令: python host_depletion.py --tool bowtie2 --seq-type rna-short
             --kraken2_index <dir> --step2_index <dir>
             --input-dir <dir> --outdir <dir> --jobs N --threads N
             [--rrna --rrna_tool silva --silva_index <dir>]

  参数:
    --input_reads    输入目录 (默认读取 00a_CleanData/) [必需]
    --output_dir     输出根目录 [必需]
    --host_db        宿主数据库根目录 (自动查找 kraken2/bowtie2/hisat2/minimap2 子目录)
    --kraken2_db     Kraken2 宿主库 (覆盖 --host_db 自动检测)
    --host_align_db  宿主比对索引 (覆盖 --host_db 自动检测)
    --aligner        {bowtie2,hisat2,minimap2} 默认 bowtie2
    --seq_type       {dna-short,rna-short,nanopore,pacbio} 默认 rna-short
    --rrna           开启 rRNA 剔除
    --rrna_tool      {ribodetector,silva} 默认 ribodetector (仅rna-short)
    --silva_index    SILVA Bowtie2 索引 (--rrna_tool silva 时必需)
    -t, --threads    线程数 (默认 20)
    -j, --jobs       并行数 (默认 2)
    --force          强制重跑

  输出: 00b_HostDepletion/
""",
    'assembly': """
  --stage assembly — 宏转录组组装 (Penguin / MEGAHIT / rnaviralSPAdes)

  子脚本: assembly_pipeline.py
    penguin         宏转录组专用组装 (默认)
    megahit         宏基因组组装
    rnaviralspades  RNA病毒组装
    all             三种工具并行, ≥2工具自动 refineC split/merge

  子脚本完整参数:
    --tool / -t           {megahit,rnaviralspades,penguin,all} 组装工具
    --input / -i          输入文件或目录
    --length / -l         contig最小长度 (默认 200)
    --threads / -n        线程数 (默认 8)
    --memory / -m         内存 GB (默认 64)
    --jobs / -j           并行任务数 (默认 1)
    --output-dir / -o     输出目录
    --log_dirs            日志目录
    --refineC_split       运行 refineC split
    --refineC_merge       运行 refineC merge
    --refineC_threads     refineC 线程
    --refineC_frag_min_len  refineC split 最小片段长度 (默认 1000)
    --refineC_min_id      refineC merge 最小序列一致性 (默认 0.97)
    --refineC_min_cov     refineC merge 最小覆盖度 (默认 0.50)
    --tmp-dir             临时目录
    --keep-temp           保留临时文件
    --force               强制重跑

  编排器透传参数:
    --input_reads    输入目录 (默认读取 00b_HostDepletion/) [必需]
    --output_dir     输出根目录 [必需]
    --assembler      {penguin,megahit,rnaviralspades,all} 默认 penguin
    -t, --threads    线程数 (默认 20)
    -m, --memory     内存 GB (默认 64)
    -j, --jobs       并行数 (默认 2)
    --force          强制重跑

  其余参数请直接调用子脚本: python scripts/assembly_pipeline.py -h

  输出: 01_Assembly/{sample}/{sample}_{tool}.contig.fasta
""",
    'identification': """
  --stage identification — 6 工具病毒序列鉴定

  子脚本: virus_identification16.py
    Genomad         深度学习病毒/质粒/前病毒分类
    Diamond BLASTX  RefSeq 病毒蛋白 + NR 库比对
    VirSorter2      隐马尔可夫模型病毒检测
    ViralVerify     病毒蛋白验证
    VirHunter       机器学习病毒鉴定
    Metabuli        基于 k-mer 的分类

  调用命令: python virus_identification16.py --input <fasta> --output <dir>
             --db_dir <dir> --identify_tools all --threads N --jobs N

  参数:
    --output_dir     输出根目录 [必需]
    --input_assembly 组装结果目录 (默认自动读取 01_Assembly/)
    --virus_db       病毒鉴定数据库根目录 [必需]
    --identify_tools 鉴定工具 (默认 all, 或逗号分隔: genomad,diamond,...)
    -t, --threads    线程数 (默认 20)
    -j, --jobs       并行数 (默认 2)
    --force          强制重跑

  输出: 02_Identification/{sample}/*_virus.all.candidate.fasta
""",
    'cobra': """
  --stage cobra — COBRA 批量延伸 (BWA-MEM2 + Contig Overlap Re-Assembly)

  子脚本: cobra_pipeline.py
    自动匹配 reads + contig + virus 三元组
    BWA-MEM2 比对 → CoverM 覆盖度 → COBRA 重叠延伸 → CheckV 评估

  调用命令: python cobra_pipeline.py --mode mix --reads-dir <dir>
             --contigs-dir <dir> --virsorter-dir <dir> --output-dir <dir>
             --assembly-tools <tools> --jobs N --threads N

  参数:
    --input_reads    项目根目录 (自动读取 00b_HostDepletion/) [必需]
    --output_dir     输出根目录 [必需]
    --assembler      {penguin,megahit,rnaviralspades,all} 默认 penguin
    -t, --threads    线程数 (默认 20)
    -j, --jobs       并行数 (默认 2)
    --force          强制重跑

  输出: 03_COBRA/{sample}/cobra_{tool}_result/*.cobra.fa
""",
    'cluster': """
  --stage cluster — 聚类 (CD-HIT 参考引导 + vclust Leiden, 仅聚类不拯救)

  子脚本: cluster_pipeline.py --stop-after-vclust
    自动收集 03_COBRA/**/*.cobra.fa, 或通过 --cluster_input 直接指定
    seqkit          最小长度过滤 (--min-length, 默认 500bp)
    CD-HIT          参考引导预聚类 (可选 --ref-genomes, ANI 95%, QCOV 85%)
                    vclust deduplicate 去重参考 → cd-hit 聚类 → 拆分 known/novel
    vclust Leiden   novel contig 的 Leiden 聚类 (ANI 95%, QCOV 85%)

  子脚本完整参数:
    -i / --input-fasta        输入 FASTA [必需]
    -o / --output-dir         输出目录 [必需]
    -t / --threads            线程数 (默认 64)
    --min-length              最小长度 bp (默认 500)
    --ani                     vclust ANI (默认 0.95)
    --qcov                    vclust QCOV (默认 0.85)
    --skip-vclust             跳过 vclust (复用已有)
    --vclust-cluster-file     复用已有聚类 TSV
    --stop-after-vclust       聚类后停止 (编排器使用)
    --ref-genomes             ICTV/NCBI 参考基因组 FASTA
    --cdhit-ani               CD-HIT ANI (默认 0.95)
    --cdhit-qcov              CD-HIT QCOV (默认 0.85)
    --resume                  断点续传

  编排器透传参数:
    --output_dir     输出根目录 [必需]
    --cluster_input  直接输入已合并 FASTA (跳过 COBRA 收集)
    --ref-genomes    ICTV/NCBI 参考基因组 FASTA
    --min-length     病毒最小长度 bp (默认 500)
    --ani            vclust ANI (默认 0.95)
    --qcov           vclust QCOV (默认 0.85)
    -t, --threads    线程数 (默认 20)

  输出:
    04_CLUSTER/centroids/final_centroids.fasta
    04_CLUSTER/centroids/known_association.tsv
    04_CLUSTER/3_vclust/vclust_clusters.tsv
    04_CLUSTER/3_vclust/split_fastas/
""",
    'taxonomy': """
  --stage taxonomy — 5 工具分类 + R 共识整合

  子脚本: virus_classifier2.py + virus_classifier_analysis14.R
    genomad        深度学习病毒分类
    mmseqs         蛋白序列比对分类
    VITAP          病毒蛋白分类
    ACVirus        古菌病毒分类
    vcontact3      蛋白簇网络分类
    R consensus    多工具投票共识 (vcontact3>vitap>acvirus>mmseqs>genomad)

  调用命令: python virus_classifier2.py -g <fasta> -s <sample>
             -t genomad,mmseqs,VITAP,ACVirus,vcontact3 -o <dir>
             -p N --db-dir <dir> [-f]
             Rscript virus_classifier_analysis14.R --combined <tsv> --output <dir>

  参数:
    --input_reads    项目根目录 [必需]
    --output_dir     输出根目录 [必需]
    --virus_db       病毒分类数据库根目录 [必需]
    -t, --threads    线程数 (默认 20)
    --force          强制重跑

  输出: 05_Taxonomy/integrated/final_integrated_classification.tsv
""",
    'host': """
  --stage host — 宿主预测 (决策树: ICTV > RNAVirHost > PhaBOX2)

  子脚本: run_host_prediction.py --mode all
    ICTV (C9)      官方权威分类库查找宿主 (--prob-dir)
    RNAVirHost     全生态位模型预测 (plant/animal/fungi/bacteria)
    PhaBOX2        CRISPR+AAI网络噬菌体宿主预测 (--phabox-db)
    决策树          Class硬规则 (Caudoviricetes→Bacteria)
                   ICTV==RVH→直接采用, 分歧时PB2决胜, 全分歧→ICT>RVH>PB2

  调用命令: python run_host_prediction.py -i <fasta> --tax <tsv> -o <dir>
             -t N --mode all [--phabox-db <dir>] [--prob-dir <dir>] [-f]

  参数:
    --input_reads    项目根目录 [必需]
    --output_dir     输出根目录 [必需]
    --phabox-db      PhaBOX2 数据库路径
    --prob-dir       ICTV 宿主概率表目录
    --virus_db       病毒数据库根目录 [必需]
    -t, --threads    线程数 (默认 20)
    --force          强制重跑

  输出:
    06_HostPrediction/ensemble_host_summary.tsv
    06_HostPrediction/host_classified_fasta/{host}.classified.fasta
""",
    'checkv': """
  --stage checkv — 按宿主分类 CheckV 预评估 (新增)

  功能: 对每个宿主类别的 centroids 运行 checkv completeness
        标记 completeness ≥ 90% 的 centroids 为免拯救 (与 CD-HIT known 合并)

  调用命令: checkv completeness <fasta> <outdir> -d <db> -t N --remove_tmp

  参数:
    --input_reads    项目根目录 [必需]
    --output_dir     输出根目录 [必需]
    --checkv_db      CheckV 数据库 [必需]
    -t, --threads    线程数 (默认 20)
    --force          强制重跑

  输出:
    07_Checkv/{host}/completeness.tsv           (各宿主 CheckV 结果)
    07_Checkv/checkv_pass_ids.txt              (≥90% pass centroids)
""",
    'rescue': """
  --stage rescue — 宿主过滤 + 三支路级联拯救 + CheckV 质量报告

  子脚本: rescue_pipeline.py
    分支 A   CheckV 并行评估 centroids (completeness >90% pass)
    分支 C   Virseqimprover reads 迭代延伸 (cluster 多样本 reads 聚合)
    分支 D   BLASTN 参考搜索 + CheckV + VSI 最后拯救
    合并     A+C+D pass → vclust 最终去重 → HQ vOTU

  免拯救 = CD-HIT known + CheckV pass (≥90%)

  子脚本完整参数:
    -c / --centroids            输入 centroids FASTA [必需]
    --clusters-tsv              vclust 聚类结果 TSV [必需]
    --split-dir                 per-cluster 拆分目录 [必需]
    -o / --output-dir           输出目录 [必需]
    -fq / --fastq-dir           reads 目录 [必需]
    -cv / --checkv-db           CheckV 数据库 [必需]
    -db / --blast-db            BLAST 数据库 (分支 D, 可选)
    --virseqimprover-path       Virseqimprover.py 路径
    --salmon-bin                Salmon 路径 (默认 salmon)
    -t / --threads              线程数 (默认 64)
    -j / --jobs                 VSI 并行数 (默认 4)
    --ani                       最终 vclust ANI (默认 0.95)
    --qcov                      最终 vclust QCOV (默认 0.85)
    --resume                    断点续传

  编排器透传参数:
    --input_reads       项目根目录 [必需]
    --output_dir        输出根目录 [必需]
    --checkv_db         CheckV 数据库 [必需]
    --blast-db          BLAST 数据库 (分支 D, 可选)
    --host-filter       目标宿主, 逗号分隔 (默认 Plant)
    --min-length        病毒最小长度 bp (默认 500)
    --ani               vclust ANI (默认 0.95)
    --qcov              vclust QCOV (默认 0.85)
    -t, --threads       线程数 (默认 20)
    -j, --jobs          并行数 (默认 2)

  输出:
    08_Rescue/{host}/centroids/final_centroids.fasta
    08_Rescue/checkv/
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
    g.add_argument('--stage', default='all',
                   choices=['all', 'clean', 'deplete', 'assembly', 'identification', 'cobra', 'cluster', 'taxonomy', 'host', 'checkv', 'rescue'],
                   help='运行阶段')
    g.add_argument('--host-filter', default='Plant',
                   help='目标宿主 (逗号分隔, rescue 阶段使用, 默认: Plant. Unknown 默认跳过并输出到 unknown_votus.fasta)')
    g.add_argument('--skip_clean', action='store_true', help='跳过数据清洗')
    g.add_argument('--skip_depletion', action='store_true', help='跳过去宿主')
    g.add_argument('--skip_clumpify', action='store_true', help='跳过 Clumpify')
    g.add_argument('--force', action='store_true', help='强制重跑')

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
    g.add_argument('--blast_mode', default='filter', help='Blast 模式 (默认 filter)')
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
    g.add_argument('--virseqimprover-path', help='Virseqimprover.py 路径')
    g.add_argument('--salmon-bin', default='salmon', help='Salmon 二进制路径 (默认: salmon)')
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
    logger = setup_logger(args.output_dir)

    stage = args.stage

    readless = stage in ('cluster', 'taxonomy', 'host', 'checkv')
    logger.info("=" * 50)
    logger.info("Virome Pipeline v2.3")
    logger.info("  Stage:    %s", stage)
    if not readless:
        logger.info("  Input:    %s", args.input_reads or args.output_dir)
    logger.info("  Output:   %s", args.output_dir)

    # 根据 stage 自动设置 skip 标志
    downstream = stage in ('assembly', 'cobra', 'rescue')
    if stage == 'clean':
        args.skip_clean = False
        args.skip_depletion = True   # clean 只做清洗, 不做去宿主
        logger.info("  Flow:     Clean only (清洗)")
    elif stage == 'deplete':
        args.skip_clean = True
        args.skip_depletion = False  # deplete 只做去宿主
        logger.info("  Flow:     Deplete only (去宿主)")
    elif downstream:
        args.skip_clean = True
        args.skip_depletion = True
        logger.info("  Flow:     %s", stage)
    elif stage in ('identification', 'cluster', 'taxonomy', 'host', 'checkv'):
        args.skip_clean = True
        args.skip_depletion = True
        logger.info("  Flow:     %s (无需 reads)", stage)
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
    logger.info("=" * 50)

    pipe = ViromePipeline(args, logger)

    if downstream or stage == 'deplete':
        # 优先使用 --input_reads 指定的目录
        if args.input_reads and Path(args.input_reads).exists():
            pipe.reads_dir = Path(args.input_reads)
        elif stage == 'deplete':
            # deplete: 自动使用 clean 输出 (优先 clumpify, 否则 fasta)
            cl = pipe.d['clean'] / '3.clumpify'
            fa = pipe.d['clean'] / '2.fasta'
            pipe.reads_dir = cl if (cl.exists() and any(cl.iterdir())) else fa
        else:
            pipe.reads_dir = pipe.d['hostdep']
        if not pipe.reads_dir.exists() or not any(pipe.reads_dir.iterdir()):
            if stage in ('cluster', 'taxonomy', 'host', 'checkv'):
                logger.info("  %s 阶段无需 reads", stage)
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
    if stage in ('all', 'clean'):
        pipe.run_clean()

    if stage in ('all', 'deplete'):
        pipe.run_depletion()

    if stage in ('all', 'assembly'):
        pipe.run_assembly()

    if stage in ('all', 'identification'):
        pipe.run_identification()

    if stage in ('all', 'cobra'):
        pipe.run_cobra()

    if stage in ('all', 'cluster'):
        pipe.run_cluster()

    if stage in ('all', 'taxonomy'):
        pipe.run_taxonomy()

    if stage in ('all', 'host'):
        pipe.run_host()

    if stage in ('all', 'checkv'):
        pipe.run_checkv_stage()

    if stage in ('all', 'rescue'):
        pipe.run_rescue()

    pipe.run_reports()


if __name__ == '__main__':
    main()
