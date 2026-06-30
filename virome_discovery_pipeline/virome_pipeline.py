#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Virome Pipeline v2.3 — 宏病毒组端到端全自动主控流水线
===============================================================================

数据流:
  Raw FASTQ ──→ [00a_CleanData] ──→ [00b_HostDepletion] ──→ [01_Assembly]
                                                                   │
   [05_Reports] ←── [04_CLUSTER] ←── [03b_MergeSamples] ←── [03a_COBRA] ←── [02_Identification] ←──┘

依赖脚本 (同目录):
  clean-data.py              → Step 0a: Fastp + Seqkit + Clumpify
  host_depletion.py          → Step 0b: Kraken2 + Align + Ribodetector
  assembly_pipeline.py       → Step 1:  MEGAHIT / rnaviralSPAdes / Penguin
  virus_identification.py  → Step 2:  Genomad + Blast + VirSorter2 + ...
  cobra_pipeline.py          → Step 3:  BWA-MEM2 + COBRA 批量延伸
  cluster_pipeline.py        → Step 4:  CLUSTER 三支路病毒基因组去冗余
  virus_classifier.py      → Step 5:  病毒分类注释 (直接调用)
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
                              text=True, bufsize=1, env={**os.environ, 'PYTHONUNBUFFERED':'1'}) as proc:
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
            # 复制到 10_Reports/logs/ (从 output_dir/orchestrator.log 反推)
            try:
                import shutil
                for h in logger.handlers:
                    if isinstance(h, logging.FileHandler):
                        reports_log = Path(h.baseFilename).resolve().parent / "10_Reports" / "logs"
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
            'bbnorm':    out / '00c_BBnorm',
            'asm':       out / '01_Assembly',
            'ident':     out / '02_Identification',
            'cobra':         out / '03_COBRA',
            '_cobra_dirs':   [out / '03a_COBRA', out / '03_COBRA'],
            'merge_samples': out / '03b_MergeSamples',
            'cluster':       out / '04_CLUSTER',
            'centroids':     out / '04_CLUSTER' / '4_centroids',  # 旧集群用 4.centroids
            '_centroids_v1': out / '04_CLUSTER' / '4.centroids',
            '_centroids_v2': out / '04_CLUSTER' / 'centroids',    # 新集群用 centroids
            'taxonomy':  out / '05_Taxonomy',
            'host_pred':  out / '06_HostPrediction',
            'checkv_dir': out / '07_Checkv',
            'rescue_dir': out / '08_Rescue',
            'analysis':  out / '09_Virome_Analysis',
            'reports':    out / '10_Reports',
        }
        # 仅创建根目录, 各阶段按需创建自己的子目录
        self.d['root'].mkdir(parents=True, exist_ok=True)

        # 脚本路径
        self.sc = {
            'clean':    script_dir / 'clean-data.py',
            'deplete':  script_dir / 'host_depletion.py',
            'assembly': script_dir / 'assembly_pipeline.py',
            'identify': script_dir / 'virus_identification.py',
            'cobra':    script_dir / 'cobra_pipeline.py',
            'cluster':  script_dir / 'cluster_pipeline.py',
            'rescue':   script_dir / 'rescue_pipeline.py',
            'host_pred': script_dir / 'run_host_prediction.py',
            'classifier': script_dir / 'virus_classifier.py',
            'classifier_R': script_dir / 'virus_classifier_analysis.R',
            'c9': script_dir / 'utils/classify_contigs.py',
            'bbnorm_script': script_dir / 'run_bbnorm.py',
        }

        # 数据流指针
        self.reads_dir = self.d['raw']

        # 验证
        self._validate()

        # 检测原始样本
        no_reads_stages = {'identification', 'merge', 'cluster', 'taxonomy', 'host', 'checkv', 'report', 'rescue'}
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

    # ── Step 0c: BBNorm 覆盖度归一化 ──
    def run_bbnorm(self):
        self.d['bbnorm'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[0c] BBNorm — 覆盖度归一化 (target=70 mindepth=2)")
        bbnorm_script = self.sc.get('bbnorm_script', SCRIPT_DIR / "run_bbnorm.py")
        parts = [
            f"python {bbnorm_script}",
            f"-i {self.reads_dir}",
            f"-o {self.d['bbnorm']}",
            f"-t {self.args.threads}",
            f"-j {getattr(self.args, 'jobs', 4)}",
        ]
        ok, _ = run_cmd(' '.join(parts), self.log, "BBNorm", str(self.d['bbnorm'] / "bbnorm.log"))
        if ok:
            self.reads_dir = self.d['bbnorm']
            self.log.info("  Reads → %s", self.reads_dir)
        else:
            self.log.warning("  BBNorm 部分失败, 继续使用未归一化 reads")

    # ── Step 1: 组装 ──
    def run_assembly(self):
        self.d['asm'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[1] MEGAHIT / rnaviralSPAdes / Penguin")

        # Co-assembly 模式: 合并所有样本 reads → 单次组装
        # (如有 BBNorm, 应在 --stage bbnorm 阶段先运行, reads 已归一化至此)
        asm_input = str(self.reads_dir)
        if getattr(self.args, 'coassembly', False):
            self.log.info("  [co-assembly] 合并所有样本 reads → 单次组装")
            merged_dir = self.d['asm'] / "coassembly_merged"
            merged_dir.mkdir(exist_ok=True)
            r1_files = sorted(Path(asm_input).glob("*_R1*.fastq.gz")) + \
                       sorted(Path(asm_input).glob("*_R1*.fq.gz")) + \
                       sorted(Path(asm_input).glob("*_1.fastq.gz")) + \
                       sorted(Path(asm_input).glob("*_1.fq.gz")) + \
                       sorted(Path(asm_input).glob("*_1.fa.gz"))
            r2_files = sorted(Path(asm_input).glob("*_R2*.fastq.gz")) + \
                       sorted(Path(asm_input).glob("*_R2*.fq.gz")) + \
                       sorted(Path(asm_input).glob("*_2.fastq.gz")) + \
                       sorted(Path(asm_input).glob("*_2.fq.gz")) + \
                       sorted(Path(asm_input).glob("*_2.fa.gz"))

            if r1_files:
                self.log.info("  合并 %d 个 R1 文件", len(r1_files))
                with open(merged_dir / "ALL_merged_R1.fq.gz", "wb") as out:
                    for f in r1_files:
                        with open(f, "rb") as inf: out.write(inf.read())
            if r2_files:
                self.log.info("  合并 %d 个 R2 文件", len(r2_files))
                with open(merged_dir / "ALL_merged_R2.fq.gz", "wb") as out:
                    for f in r2_files:
                        with open(f, "rb") as inf: out.write(inf.read())
            asm_input = str(merged_dir)
            # 更新 reads_dir → rescue 阶段使用合并后的 reads
            self.coassembly_merged_dir = str(merged_dir)

        # assembly_pipeline.py 完整参数:
        #   --tool, --input, --length, --threads, --memory, --jobs,
        #   --output-dir, --log_dirs, --refineC_split, --refineC_merge,
        #   --refineC_threads, --refineC_frag_min_len, --refineC_min_id, --refineC_min_cov,
        #   --tmp-dir, --keep-temp, --force
        parts = [
            f"python {self.sc['assembly']}",
            f"--tool {self.args.assembler}",
            f"--input {asm_input}",
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

        asm_dir = Path(self.args.input_assembly) if self.args.input_assembly else self.d['asm']
        asm_map = getattr(self, 'asm_map', {})
        if not asm_map:
            scan_tools = ['megahit', 'rnaviralspades', 'penguin']
            asm_map = scan_contig_files(asm_dir, scan_tools)
        if not asm_map:
            self.log.error("无 contig 文件, 跳过鉴定。")
            return

        # 直接传目录, 子脚本内部 --jobs 并行处理所有样本
        self.log.info("  %d 个样本待鉴定", len(asm_map))

        # virus_identification.py 完整参数:
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

    # ── Step 3a: COBRA 延伸 ──
    # ── 自动探测 COBRA 目录 (兼容 03a_COBRA / 03_COBRA) ──
    def _resolve_cobra_dir(self):
        for d in self.d.get('_cobra_dirs', [self.d['root'] / '03a_COBRA']):
            if d.is_dir():
                self.d['cobra'] = d
                return d
        self.d['cobra'] = self.d['root'] / '03a_COBRA'
        return self.d['cobra']

    # ── Step 3a: COBRA 延伸 ──
    def run_cobra(self):
        self._resolve_cobra_dir()
        self.d['cobra'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[3a] COBRA 批量延伸 (BWA-MEM2 + COBRA + CheckV)")

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
            f"--jobs {self.args.cobra_jobs or max(1, self.args.jobs // 2)}",
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

    # ── Step 3b: 合并多样本 COBRA 结果 + 可选 Flye 共组装 ──
    def run_merge_samples(self):
        self._resolve_cobra_dir()
        """
        收集 03a_COBRA 下所有样本的结果, 合并为单个 FASTA。
        默认对合并结果跑 Flye --subassemblies; --skip-flye 可跳过。
        最终输出 Flye 延伸版 + 原始合并版的拼接文件。
        产出: 03b_MergeSamples/all_sample_virus.fasta (或 _combined.fasta)
        """
        merge_dir = self.d['merge_samples']
        merge_dir.mkdir(parents=True, exist_ok=True)

        self._resolve_cobra_dir()
        self.log.info("=" * 50)
        self.log.info("[3b] 合并多样本 COBRA 结果")

        cobra_dir = self.d['cobra']
        if not cobra_dir.is_dir():
            self.log.error("COBRA 输出 %s 不存在, 跳过 merge", cobra_dir)
            sys.exit(1)

        # ── 收集合并 ──
        merged_fa = merge_dir / "all_sample_virus.fasta"
        n_cobra, n_virus = 0, 0
        with open(merged_fa, 'w') as out:
            for sd in cobra_dir.iterdir():
                if not sd.is_dir():
                    continue
                cobra_files = list(sd.rglob('*.cobra.fa'))
                if cobra_files:
                    for cf in cobra_files:
                        with open(cf) as inf:
                            out.write(inf.read())
                            n_cobra += 1
                else:
                    for vf in sd.rglob('*virus*.fasta'):
                        with open(vf) as inf:
                            out.write(inf.read())
                            n_virus += 1
                    for vf in sd.rglob('*virus*.fa'):
                        with open(vf) as inf:
                            out.write(inf.read())
                            n_virus += 1

        if merged_fa.stat().st_size == 0:
            self.log.error("合并后 FASTA 为空, 未找到任何 COBRA 或 virus 结果")
            sys.exit(1)

        n_input = n_cobra + n_virus
        total_bp = sum(len(rec.seq) for rec in SeqIO.parse(str(merged_fa), "fasta"))
        if n_virus:
            self.log.info("  合并: %d cobra + %d virus = %d 条, %.1f Mb → %s", n_cobra, n_virus, n_input, total_bp / 1e6, merged_fa)
        else:
            self.log.info("  合并: %d cobra, %.1f Mb → %s", n_cobra, total_bp / 1e6, merged_fa)

        # ── Flye 共组装 (可选) ──
        if getattr(self.args, 'skip_flye', False):
            self.log.info("  merge 阶段完成 (跳过 Flye)")
            return

        flye_dir = merge_dir / "flye_coassembly"
        flye_dir.mkdir(parents=True, exist_ok=True)
        flye_out = flye_dir / "flye_output"

        min_ovlp = getattr(self.args, 'flye_min_overlap', 500)
        read_err = getattr(self.args, 'flye_read_error', 0.005)

        self.log.info("-" * 40)
        self.log.info("  Flye --subassemblies: min_overlap=%d, read_error=%.3f", min_ovlp, read_err)

        flye_cmd = (
            f"flye --subassemblies {merged_fa}"
            f" -t {self.args.threads}"
            f" --meta"
            f" --read-error {read_err}"
            f" -m {min_ovlp}"
            f" -o {flye_out}"
        )

        ok, _ = run_cmd(flye_cmd, self.log, "Flye co-assembly", str(flye_dir / "flye.log"))
        if not ok:
            self.log.warning("Flye 运行失败, merge 输出为原始合并结果")
            return

        flye_assembly = flye_out / "assembly.fasta"
        if not flye_assembly.exists() or flye_assembly.stat().st_size == 0:
            self.log.warning("Flye 未产出有效 contig, merge 输出为原始合并结果")
            return

        n_flye = sum(1 for _ in open(flye_assembly) if _.startswith('>'))
        flye_bp = sum(len(rec.seq) for rec in SeqIO.parse(str(flye_assembly), "fasta"))
        self.log.info("  Flye 产出: %d 条 contig, %.1f Mb", n_flye, flye_bp / 1e6)

        # 溯源报告
        trace_script = SCRIPT_DIR / "utils" / "flye_trace_native.py"
        if trace_script.exists():
            run_cmd(
                f"python {trace_script} -i {flye_out} -o {flye_dir / 'flye_mapping.tsv'} -s {flye_dir / 'flye_summary.txt'} -p {flye_dir / 'flye_report.png'} --full",
                self.log, "Flye trace", str(flye_dir / "trace.log")
            )

        # 合并: Flye + 原始 → combined
        combined_fa = merge_dir / "all_sample_virus_combined.fasta"
        with open(combined_fa, 'w') as out:
            with open(flye_assembly) as inf:
                out.write(inf.read())
            with open(merged_fa) as inf:
                out.write(inf.read())

        n_combined = sum(1 for _ in open(combined_fa) if _.startswith('>'))
        combined_bp = sum(len(rec.seq) for rec in SeqIO.parse(str(combined_fa), "fasta"))
        self.log.info("  合并: Flye %d + 原始 %d = %d 条, %.1f Mb → %s",
                      n_flye, n_input, n_combined, combined_bp / 1e6, combined_fa)
        self.log.info("  merge 阶段完成 (含 Flye 共组装)")

    # ── Step 4: CLUSTER 三支路去冗余 ──
    def run_cluster(self):
        self._resolve_cobra_dir()
        self.d['cluster'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[4] CLUSTER 三支路病毒基因组去冗余")

        # 确定 cluster 输入
        if self.args.cluster_input and os.path.isfile(self.args.cluster_input):
            cluster_fa = Path(self.args.cluster_input)
            self.log.info("  CLUSTER 直接输入: %s", cluster_fa)
        else:
            # 优先取 merge stage 的输出
            merge_out = self.d['merge_samples']
            combined_fa = merge_out / "all_sample_virus_combined.fasta"
            raw_fa = merge_out / "all_sample_virus.fasta"

            if combined_fa.exists() and combined_fa.stat().st_size > 0:
                cluster_fa = combined_fa
                self.log.info("  CLUSTER 输入 (merge+Flye): %s", cluster_fa)
            elif raw_fa.exists() and raw_fa.stat().st_size > 0:
                cluster_fa = raw_fa
                self.log.info("  CLUSTER 输入 (merge): %s", cluster_fa)
            else:
                # 向后兼容: 现场收集
                self.log.info("  未找到 merge 输出, 现场收集 COBRA 结果...")
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
                        cobra_files = list(sd.rglob('*.cobra.fa'))
                        if cobra_files:
                            for cf in cobra_files:
                                with open(cf) as inf:
                                    out.write(inf.read())
                                    n_cobra += 1
                        else:
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
                self.log.info("  CLUSTER 输入 (现场收集): %d cobra + %d virus", n_cobra, n_virus)

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
        if getattr(self.args, 'skip_rmdup', False):
            parts.append("--skip-rmdup")
        if getattr(self.args, 'rmdup_length', None):
            parts.append(f"--rmdup-length {self.args.rmdup_length}")
        if not self.args.force:
            parts.append("--resume")

        ok, _ = run_cmd(' '.join(parts), self.log, "CLUSTER (vclust only)", str(self.d['cluster'] / "cluster.log"))
        if not ok:
            self.log.error("CLUSTER vclust 阶段失败, 终止。")
            sys.exit(1)

        centroids = self.d['centroids'] / "final_centroids.fasta"
        if not centroids.is_file():
            for alt_key in ['_centroids_v1', '_centroids_v2']:
                centroids = self.d[alt_key] / "final_centroids.fasta"
                if centroids.is_file(): break
        if centroids.is_file():
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

        centroids_fa = self.d['centroids'] / "final_centroids.fasta"
        if not centroids_fa.is_file():
            centroids_fa = self.d['centroids'] / "final_centroids.fasta"
        for alt_key in ['_centroids_v1', '_centroids_v2']:
            if centroids_fa.is_file(): break
            centroids_fa = self.d[alt_key] / "final_centroids.fasta"
        if not centroids_fa.is_file():
            self.log.error("centroids 不存在 (试了 4_centroids, 4.centroids, centroids), 请先运行 --stage cluster")
            sys.exit(1)

        clusters_tsv = self.d['cluster'] / "3_vclust" / "vclust_clusters.tsv"
        if not clusters_tsv.is_file():
            self.log.error("clusters.tsv %s 不存在", clusters_tsv)
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
            unknown_out = self.d['centroids'] / "unknown_votus.fasta"
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
            host_out = self.d['centroids'] / f"skipped_{safe_host}.fasta"
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
        known_id_file = self.d['centroids'] / "known_ids.txt"
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
            f"-fq {getattr(self, 'coassembly_merged_dir', None) or self.reads_dir}",
            f"-cv {self.args.checkv_db}",
            f"-t {self.args.threads}",
            f"-j {self.args.jobs}",
            f"--ani {self.args.ani}",
            f"--qcov {self.args.qcov}",
        ]
        # 分支 C 默认用 centroids + cdhit_combined 自动建库
        # 如需外部 BLAST DB, 在 run_config.json 中加 "extra_blast_db" 即可
        if getattr(self.args, 'extra_blast_db', None):
            parts.append(f"-db {getattr(self.args, 'extra_blast_db')}")
        vsi_path = self.args.virseqimprover_path or str(self.sc['cluster'].parent / 'Virseqimprover.py')
        parts += [f"--virseqimprover-path {vsi_path}"]
        parts += [f"--salmon-bin {self.args.salmon_bin}"]
        parts += [f"--max-vsi-samples {self.args.max_vsi_samples}"]
        parts += [f"--min-vsi-len {self.args.min_vsi_len}"]
        parts += [f"--checkv-threshold {getattr(self.args, 'checkv_threshold', 90.0)}"]
        # Taxonomy + genus_len (分支 D + VSI genus_avg_len 备选截止)
        sample = getattr(self.args, 'tax_sample_name', None) or "Votus"
        tax_tsv = self.d['taxonomy'] / f"{sample}.integrated" / "final_integrated_classification.tsv"
        if tax_tsv.is_file():
            parts.append(f"--taxonomy-tsv {tax_tsv}")
        genus_len_path = Path(os.path.expanduser("~/database/virus-db/db/genus_lens"))
        if not genus_len_path.is_file():
            genus_len_path = self.script_dir.parent / "database" / "genus_lens"
        if genus_len_path.is_file():
            parts.append(f"--genus-len {genus_len_path}")
        if not self.args.force:
            parts.append("--resume")
        # cdhit_combined.fasta (供分支C自动建库)
        cdhit_fa = self.d['cluster'] / "2_cdhit" / "cdhit_combined.fasta"
        if cdhit_fa.is_file():
            parts.append(f"--cdhit-fa {cdhit_fa}")

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
        unknown_fa = self.d['centroids'] / "unknown_votus.fasta"
        if unknown_fa.exists():
            dist, total = self._run_checkv_on_fasta(unknown_fa, self.d['rescue_dir'] / "checkv" / "unknown")
            all_stats["Unknown"] = (dist, total)

        # 3. 跳过的宿主 centroids
        seen_hosts = set()
        for f in (self.d['centroids']).glob("skipped_*.fasta"):
            host_label = f.stem.replace("skipped_", "").replace("_", " ")
            if host_label in seen_hosts:
                continue
            seen_hosts.add(host_label)
            dist, total = self._run_checkv_on_fasta(f, self.d['rescue_dir'] / "checkv" / f"skipped_{host_label}")
            all_stats[f"Skipped_{host_label}"] = (dist, total)

        # 4. 汇总 centroids
        all_centroids = self.d['centroids'] / "final_centroids.fasta"
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
        self.log.info("[5] 病毒分类注释 (virus_classifier.py + R 共识整合)")

        centroids = self.d['centroids'] / "final_centroids.fasta"
        if not centroids.exists():
            self.log.warning("未找到 centroids, 跳过分类")
            return

        tax_dir = self.d['taxonomy']
        sample = getattr(self.args, 'tax_sample_name', None) or "Votus"
        out_subdir = tax_dir / f"{sample}.classed"
        combined_tsv = out_subdir / f"{sample}_combined_taxonomy.tsv"
        int_dir = tax_dir / f"{sample}.integrated"
        final_tsv = int_dir / "final_integrated_classification.tsv"

        # virus_classifier.py 完整参数:
        #   -g, -s, -t, -o, -p, -f, --db-dir,
        #   --genomad-db, --metabuli-db, --cat-db, --cat-tax, --uniprot-db,
        #   --mmseqs-db, --vitap-db, --acvirus-db, --vcontact3-db,
        #   -j, -e, --remove-suffix
        if not self.args.force and combined_tsv.is_file() and combined_tsv.stat().st_size > 100:
            self.log.info("  [SKIP] virus_classifier — 已有结果")
        else:
            out_subdir.mkdir(parents=True, exist_ok=True)
            parts = [
                f"python {self.sc['classifier']}",
                f"-g {centroids}",
                f"-s {sample}",
                f"-t {self.args.tax_tools}",
                f"-o {out_subdir}",
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
            ok, _ = run_cmd(' '.join(parts), self.log, "virus_classifier.py", str(self.d['taxonomy'] / "taxonomy.log"))
            if not ok:
                self.log.warning("  virus_classifier 失败")

        # Step 2: R 共识整合
        if not self.args.force and final_tsv.is_file() and final_tsv.stat().st_size > 100:
            self.log.info("  [SKIP] R consensus — 已有结果")
        elif combined_tsv.is_file():
            int_dir.mkdir(parents=True, exist_ok=True)
            parts = [
                f"cd {int_dir} &&",
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

        centroids = self.d['centroids'] / "final_centroids.fasta"
        sample = getattr(self.args, 'tax_sample_name', None) or "Votus"
        tax_tsv = self.d['taxonomy'] / f"{sample}.integrated" / "final_integrated_classification.tsv"
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
            f"--mode {self.args.host_mode}",
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
        centroids_fa = self.d['centroids'] / "final_centroids.fasta"

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

    # ── Step 09: 病毒基因组下游分析 + GenBank 提交 ──
    def run_analysis(self):
        """[9] 病毒基因组下游分析 + suvtk 注释 (合并阶段)

        1. suvtk taxonomy → features → hypothetical → tbl2gb → topology
        2. integrate_taxonomy_suvtk (05_Taxonomy + suvtk 交叉补充)
        3. virome_analysis.py (丰度/多样性等)
        4. integrated_summary.py (三源交叉验证)
        """
        self.d['analysis'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[9] Virome Analysis — 下游分析 + suvtk 注释")

        all_plant = self.d['rescue_dir'] / "all_plant_viruses.fasta"
        if not all_plant.is_file():
            self.log.warning("  all_plant_viruses.fasta 不存在, 跳过")
            return

        # 复制到 analysis 目录
        import shutil
        dest_fasta = self.d['analysis'] / "all_plant_viruses.fasta"
        if not dest_fasta.exists():
            shutil.copy(all_plant, dest_fasta)

        n = sum(1 for _ in open(all_plant) if _.startswith('>'))
        self.log.info("  输入: %d 条植物病毒", n)

        t = getattr(self.args, 'threads', 40)
        suvtk_db = getattr(self.args, 'suvtk_db', None) or os.path.expanduser('~/database/virus-db/suvtk_db')
        rvdb = getattr(self.args, 'rvdb_db', None) or os.path.expanduser('~/database/virus-db/RVDB-v31')

        # ── Part 1: suvtk taxonomy + features ──
        tax_out = self.d['analysis'] / "suvtk.taxonomy_output"
        tax_tsv = tax_out / "taxonomy.tsv"
        if not tax_tsv.is_file() or tax_tsv.stat().st_size < 100:
            cmd = f"suvtk taxonomy -i {dest_fasta} -o {tax_out} -d {suvtk_db} -s 0.7 -t {t}"
            run_cmd(cmd, self.log, "suvtk taxonomy", str(self.d['analysis'] / "suvtk_taxonomy.log"))
        else:
            self.log.info("  suvtk taxonomy — 已存在, 跳过")

        feat_out = self.d['analysis'] / "suvtk.features_output"
        tbl = feat_out / "featuretable.tbl"
        if not tbl.is_file() or tbl.stat().st_size < 100:
            cmd = f"suvtk features -i {dest_fasta} -o {feat_out} -d {suvtk_db} --coding-complete"
            if tax_tsv.is_file():
                cmd += f" --taxonomy {tax_tsv}"
            cmd += f" -t {t}"
            run_cmd(cmd, self.log, "suvtk features", str(self.d['analysis'] / "suvtk_features.log"))
        else:
            self.log.info("  suvtk features — 已存在, 跳过")

        # ── Part 2: analyze_hypothetical ──
        hypo_out = self.d['analysis'] / "analyze_hypothetical"
        updated_tbl = hypo_out / "featuretable_updated.tbl"
        hypo_script = SCRIPT_DIR.parent / "virome_submission_pipeline" / "analyze_hypothetical.py"
        if not hypo_script.is_file():
            hypo_script = SCRIPT_DIR.parent / "virome_analysis_pipeline" / "analyze_hypothetical.py"
        if tbl.is_file() and (feat_out / "proteins.faa").is_file():
            if not updated_tbl.is_file() or updated_tbl.stat().st_size < 100:
                rvdb_fasta = Path(rvdb) / "U-RVDBv31.0-prot.EX.acc.fasta"
                rvdb_hmm = Path(rvdb) / "U-RVDBv31.0-prot.hmm"
                rvdb_annot = Path(rvdb) / "U-RVDBv31.0-prot.info.tab"
                cmd = (
                    f"python {hypo_script} "
                    f"-t {tbl} -f {feat_out / 'proteins.faa'} "
                    f"-o {hypo_out} "
                    f"--blast {hypo_out / 'merged_blast.txt'} "
                    f"--diamond -d {rvdb_fasta} "
                    f"--hmmer --hmmer-db {rvdb_hmm} "
                    f"--annot {rvdb_annot} "
                    f"--threads {t}"
                )
                run_cmd(cmd, self.log, "analyze_hypothetical",
                        str(self.d['analysis'] / "hypothetical.log"))
            else:
                self.log.info("  analyze_hypothetical — 已存在, 跳过")

        # ── Part 3: tbl2gb (GenBank) ──
        tbl2gb_script = SCRIPT_DIR.parent / "virome_analysis_pipeline" / "tbl2gb.py"
        if not tbl2gb_script.is_file():
            tbl2gb_script = SCRIPT_DIR.parent / "virome_submission_pipeline" / "tbl2gb.py"
        nuc_fna = feat_out / "reoriented_nucleotide_sequences.fna"
        if not nuc_fna.is_file():
            nuc_fna = feat_out / "nucleotide_sequences.fasta"
        gb_out = self.d['analysis'] / "virus-annotations"
        if tbl2gb_script.is_file() and updated_tbl.is_file() and nuc_fna.is_file():
            if not gb_out.is_dir() or not any(gb_out.glob("*.gb")):
                gb_out.mkdir(exist_ok=True)
                cmd = f"python {tbl2gb_script} --tbl {updated_tbl} --fasta {nuc_fna} -o {gb_out}"
                run_cmd(cmd, self.log, "tbl2gb (GenBank)", str(self.d['analysis'] / "tbl2gb.log"))
                n_gb = len(list(gb_out.glob("*.gb"))) if gb_out.is_dir() else 0
                self.log.info("  tbl2gb — 生成 %d 个 .gb 文件", n_gb)
            else:
                self.log.info("  tbl2gb — 已存在, 跳过")
        else:
            self.log.warning("  tbl2gb — 跳过 (缺 tbl/fna/script)")

        # ── Part 4: viral_topology ──
        topo_out = self.d['analysis'] / "topology.tsv"
        if not topo_out.is_file() or topo_out.stat().st_size < 50:
            topo_script = SCRIPT_DIR.parent / "virome_submission_pipeline" / "viral_topology.py"
            cmd = f"python {topo_script} --taxonomy {tax_tsv} --fasta {dest_fasta} -o {topo_out}"
            run_cmd(cmd, self.log, "viral_topology", str(self.d['analysis'] / "topology.log"))

        # ── Part 5: 整合 05_Taxonomy + suvtk ──
        self._integrate_taxonomy_suvtk(dest_fasta, tax_tsv, hypo_out)

        # ── Part 6: 病毒组下游分析 ──
        analysis_script = SCRIPT_DIR / "virome_analysis.py"
        cmd = f"python {analysis_script} -i {dest_fasta} -o {self.d['analysis']} -t {t}"
        ok, _ = run_cmd(cmd, self.log, "Virome Analysis",
                        str(self.d['analysis'] / "analysis.log"))
        if ok:
            integ_script = SCRIPT_DIR / "integrated_summary.py"
            if integ_script.is_file():
                run_cmd(f"python {integ_script} -o {self.d['root']}",
                        self.log, "Integrated Summary",
                        str(self.d['analysis'] / "integrated_summary.log"))
        else:
            self.log.warning("  分析部分失败, 检查日志")

        self.log.info("  分析完成 → %s", self.d['analysis'])

    def run_suvtk_annotate(self):
        """[已废弃] 合并到 run_analysis()"""
        self.run_analysis()

    def _integrate_taxonomy_suvtk(self, fasta_path, tax_tsv, hypo_dir):
        """整合 05_Taxonomy 共识分类 + suvtk → 丰富化输出

        产出:
          miuvig_taxonomy_enriched.tsv  — pred_genome_type/structure 用 consensus 修正
          integrated_summary.tsv        — 两源对比表
          topology.tsv (enriched)       — 增加 taxonomy 列
        """
        import pandas as pd
        from Bio import SeqIO

        # ── 加载数据 ──
        # 05_Taxonomy 共识
        consensus_path = self.d['taxonomy'] / "Votus.integrated" / "final_integrated_classification.tsv"
        consensus = {}
        if consensus_path.exists():
            try:
                cdf = pd.read_csv(consensus_path, sep='\t', quotechar='"')
                for _, row in cdf.iterrows():
                    cid = str(row.get('contig_id', row.iloc[0]))
                    consensus[cid] = {
                        'Realm': str(row.get('Realm', '')), 'Kingdom': str(row.get('Kingdom', '')),
                        'Phylum': str(row.get('Phylum', '')), 'Class': str(row.get('Class', '')),
                        'Order': str(row.get('Order', '')), 'Family': str(row.get('Family', '')),
                        'Genus': str(row.get('Genus', '')), 'Species': str(row.get('Species', '')),
                        'confidence': str(row.get('confidence', '')),
                        'primary_tool': str(row.get('primary_tool', '')),
                    }
            except Exception:
                pass

        # suvtk taxonomy
        suvtk = {}
        if tax_tsv and tax_tsv.exists():
            sdf = pd.read_csv(tax_tsv, sep='\t')
            for _, row in sdf.iterrows():
                suvtk[str(row.iloc[0])] = str(row.iloc[1])

        # suvtk miuvig_taxonomy
        miuvig_tax = {}
        miuvig_path = self.d['analysis'] / "suvtk.taxonomy_output" / "miuvig_taxonomy.tsv"
        if miuvig_path.exists():
            mdf = pd.read_csv(miuvig_path, sep='\t')
            for _, row in mdf.iterrows():
                miuvig_tax[str(row.iloc[0])] = (
                    str(row.get('pred_genome_type', '')), str(row.get('pred_genome_struc', ''))
                )

        # suvtk topology
        topo_path = self.d['analysis'] / "topology.tsv"
        topology = {}
        if topo_path.exists():
            tdf = pd.read_csv(topo_path, sep='\t')
            for _, row in tdf.iterrows():
                topology[str(row.iloc[0])] = {
                    'final_topology': str(row.get('final_topology', '')),
                    'evidence': str(row.get('evidence', '')),
                }

        # 序列长度 & CDS 数
        seq_lens = {}
        for rec in SeqIO.parse(str(fasta_path), 'fasta'):
            seq_lens[rec.id] = len(rec.seq)

        cds_counts = {}
        tbl_path = hypo_dir / "featuretable_updated.tbl"
        if not tbl_path.exists():
            tbl_path = self.d['analysis'] / "suvtk.features_output" / "featuretable.tbl"
        if tbl_path.exists():
            cur = None
            with open(tbl_path) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>Feature '):
                        cur = line.split('>Feature ')[1].strip()
                        cds_counts.setdefault(cur, 0)
                    elif line and line.split()[-1:] == ['CDS'] and cur:
                        cds_counts[cur] += 1

        # ── hypothetical status (from analyze_hypothetical) ──
        hypo_status = {}
        updated_tbl = hypo_dir / "featuretable_updated.tbl"
        if updated_tbl.exists():
            cur_contig = None; hypo_n = 0; anno_n = 0
            for line in open(updated_tbl):
                line = line.strip()
                if line.startswith('>Feature '):
                    if cur_contig and (hypo_n + anno_n) > 0:
                        hypo_status[cur_contig] = f'{anno_n} annotated, {hypo_n} hypothetical'
                    cur_contig = line.split('>Feature ')[1].strip()
                    hypo_n = 0; anno_n = 0
                elif line.startswith('\t\t\tproduct\t') and cur_contig:
                    if 'hypothetical' in line.lower():
                        hypo_n += 1
                    else:
                        anno_n += 1
            if cur_contig and (hypo_n + anno_n) > 0:
                hypo_status[cur_contig] = f'{anno_n} annotated, {hypo_n} hypothetical'

        # ── 生成 enriched miuvig_taxonomy ──
        enriched_miuvig = self.d['analysis'] / "suvtk.taxonomy_output" / "miuvig_taxonomy_enriched.tsv"
        with open(enriched_miuvig, 'w') as f:
            f.write("contig\tpred_genome_type\tpred_genome_struc\t"
                    "suvtk_taxonomy\tconsensus_Family\tconsensus_Genus\tconsensus_Species\t"
                    "confidence\tprimary_tool\n")
            for cid in seq_lens:
                mt = miuvig_tax.get(cid, ('uncharacterized', 'undetermined'))
                sv = suvtk.get(cid, '')
                cs = consensus.get(cid, {})
                f.write(f"{cid}\t{mt[0]}\t{mt[1]}\t"
                        f"{sv}\t{cs.get('Family','')}\t{cs.get('Genus','')}\t{cs.get('Species','')}\t"
                        f"{cs.get('confidence','')}\t{cs.get('primary_tool','')}\n")
        self.log.info("  miuvig_taxonomy_enriched.tsv (%d 条)", len(seq_lens))

        # ── 生成 integrated_summary (两源对比) ──
        int_summary = self.d['analysis'] / "integrated_summary.tsv"
        with open(int_summary, 'w') as f:
            f.write("contig_id\tlength\tcds_count\t"
                    "suvtk_taxonomy\tconsensus_Family\tconsensus_Genus\tconsensus_Species\t"
                    "topology\tevidence\tconfidence\tprimary_tool\n")
            for cid, slen in seq_lens.items():
                sv = suvtk.get(cid, '')
                cs = consensus.get(cid, {})
                tp = topology.get(cid, {})
                f.write(f"{cid}\t{slen}\t{cds_counts.get(cid, 0)}\t"
                        f"{sv}\t{cs.get('Family','')}\t{cs.get('Genus','')}\t{cs.get('Species','')}\t"
                        f"{tp.get('final_topology','')}\t{tp.get('evidence','')}\t"
                        f"{cs.get('confidence','')}\t{cs.get('primary_tool','')}\n")
        self.log.info("  integrated_summary.tsv (%d 条)", len(seq_lens))

        # ── 生成 ref_info.tsv ──
        ref_out = self.d['analysis'] / "ref_info.tsv"
        with open(ref_out, 'w') as f:
            f.write("Accession\tLength\tSpecies\tGenus\tFamily\tRealm\t"
                    "Kingdom\tClass\tOrder\t"
                    "suvtk_Species\tsuvtk_Genus\tsuvtk_Family\t"
                    "CDS_Count\tPrimary_Tool\tConfidence\t"
                    "Molecule_type\tMolecule_Type2\tSegment\t"
                    "topology\ttopology_evidence\thypothetical_status\n")
            for cid, slen in seq_lens.items():
                cs = consensus.get(cid, {})
                cds = cds_counts.get(cid, 0)
                mt = miuvig_tax.get(cid, ('uncharacterized', 'undetermined'))
                mol_type = mt[0]  # ssRNA(+)/ssRNA(-)/dsRNA
                mol_type2 = 'RNA' if 'RNA' in mol_type else ('DNA' if 'DNA' in mol_type else '')
                segment = 'Unsegmented' if mt[1] != 'segmented' else 'Segmented'
                tp = topology.get(cid, {})
                hypo = hypo_status.get(cid, '')
                f.write(f"{cid}\t{slen}\t"
                        f"{cs.get('Species','')}\t{cs.get('Genus','')}\t{cs.get('Family','')}\t{cs.get('Realm','')}\t"
                        f"{cs.get('Kingdom','')}\t{cs.get('Class','')}\t{cs.get('Order','')}\t"
                        f"{cs.get('Species','')}\t{cs.get('Genus','')}\t{cs.get('Family','')}\t"
                        f"{cds}\t{cs.get('primary_tool','MMPV-RNA+suvtk')}\t{cs.get('confidence','')}\t"
                        f"{mol_type}\t{mol_type2}\t{segment}\t"
                        f"{tp.get('final_topology','')}\t{tp.get('evidence','')}\t{hypo}\n")
        self.log.info("  ref_info.tsv (%d 条)", len(seq_lens))

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

STAGE_HELP = {
    'clean':         '00a: Fastp质控 + Seqkit转FASTA + Clumpify去冗余',
    'deplete':       '00b: Kraken2 + Bowtie2/HISAT2/Minimap2 + rRNA去除',
    'bbnorm':        '00c: BBNorm覆盖度均一化 (可选, co-assembly前)',
    'assembly':      '01:  MEGAHIT / rnaviralSPAdes / Penguin 三工具组装',
    'identification':'02:  10工具并行病毒鉴定 (Genomad/Diamond/VirSorter2/...)',
    'cobra':         '03a: BWA-MEM2 + COBRA-Meta 重叠群延伸',
    'merge':         '03b: 合并多样本COBRA结果 + 可选Flye共组装',
    'cluster':       '04:  CD-HIT参考预聚类 + vclust Leiden聚类',
    'taxonomy':      '05:  8工具分类 + R加权投票共识',
    'host':          '06:  ICTV > RNAVirHost > PhaBOX2 宿主预测',
    'checkv':        '07:  CheckV 完整性评估',
    'rescue':        '08:  三支路级联拯救 (CheckV -> Virseqimprover -> BLASTN)',
    'analysis':      '09:  suvtk + integrated summary -> GenBank准备',
    'report':        '10:  TSV汇总 + Sankey图 + 交互式HTML报告',
}

STAGE_ARGS = {
    "clean":         ["clean-data"],
    "deplete":       ["host_depletion"],
    "bbnorm":        [],
    "assembly":      ["assembly_pipeline"],
    "identification":["virus_identification"],
    "cobra":         ["cobra_pipeline"],
    "merge":         [],
    "cluster":       ["cluster_pipeline"],
    "taxonomy":      ["virus_classifier"],
    "host":          ["run_host_prediction"],
    "checkv":        [],
    "rescue":        ["cluster_pipeline"],
    "analysis":      [],
    "report":        [],
}

OVERVIEW = """
╔══════════════════════════════════════════════════════════════╗
║   MMPV-RNA — 宏病毒组端到端全自动分析流水线          ║
╚══════════════════════════════════════════════════════════════╝

Stage 流程 (顺序执行):
  clean         → 原始数据质控 (fastp)
  deplete       → 宿主去除 (bowtie2)
  assembly      → 组装 (rnaviralSPAdes)
  identification→ 病毒鉴定 (geNomad)
  cobra         → 跨样本聚类 (COBRA)
  cluster       → vOTU 聚类 (CD-HIT + vclust)
  taxonomy      → 分类学注释 (mmseqs + genomad)
  host          → 宿主预测
  checkv        → 完整性评估 (CheckV)
  rescue        → 三支路级联拯救 (CheckV → VSI → BLASTN → genus_len)
  analysis      → 病毒组下游分析
  suvtk_annotate→ suvtk 注释 + tbl2gb (GenBank)
  report        → HTML 报告生成

用法:
  python virome_pipeline.py --stage rescue --output_dir <DIR>
  python virome_pipeline.py --stage rescue,analysis,suvtk_annotate,report --output_dir <DIR>
  python virome_pipeline.py --stage all --output_dir <DIR> --input_reads <FASTQ_DIR>
"""

def print_help_and_exit():
    print(OVERVIEW)
    sys.exit(0)


def _build_parser(add_help=True):
    """构建 argparse, 可复用"""
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='宏病毒组端到端全自动主控流水线', add_help=add_help)

    g = p.add_argument_group('路径配置')
    g.add_argument('--input_reads', help='原始 FASTQ 目录 (--cluster_input 模式可省略)')
    g.add_argument('--output_dir', required=True, help='项目输出根目录')
    g.add_argument('--cluster_input', help='直接输入已合并的病毒 FASTA (跳过 COBRA 收集)')
    g.add_argument('--input_assembly', help='组装结果目录 (identification 阶段, 默认 01_Assembly/)')

    g = p.add_argument_group('流程控制')
    g.add_argument('--config', default=None, help='YAML 配置文件路径 (默认: 自动查找 pipeline_config.yaml)')
    g.add_argument('--profile', default='default', help='配置预设 (默认: default, 可选: downstream/plant)')
    g.add_argument('--dump-config', action='store_true', help='仅打印配置摘要并退出 (不运行)')
    g.add_argument('--stage', default=['all'], nargs='+',
                   choices=['all', 'clean', 'deplete', 'bbnorm', 'assembly', 'identification', 'cobra', 'merge', 'cluster', 'taxonomy', 'host', 'checkv', 'rescue', 'analysis', 'suvtk_annotate', 'report'],
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

    g = p.add_argument_group('Identification 阶段 (virus_identification.py)')
    g.add_argument('--nr_db', default='~/database/nr_diamond/nr.dmnd', help='Diamond NR 数据库路径')
    g.add_argument('--skip_uniprot_filter', action='store_true', help='跳过 UniProt 后置过滤')
    g.add_argument('--skip_nr_filter', action='store_true', help='跳过 NR 后置过滤')
    g.add_argument('--skip_id_plots', action='store_true', help='跳过鉴定阶段图表生成')
    g.add_argument('--clean_failed', action='store_true', help='自动清理鉴定失败的任务目录')
    g.add_argument('--ident_ext', default='.fasta', help='输入目录时搜索的后缀 (默认: .fasta)')

    g = p.add_argument_group('COBRA 阶段 (cobra_pipeline.py)')
    g.add_argument('--cobra_mink', type=int, default=21, help='COBRA 最小 kmer (默认: 21)')
    g.add_argument('--cobra_maxk', type=int, default=141, help='COBRA 最大 kmer (默认: 141)')
    g.add_argument('--cobra_jobs', type=int, default=0, help='COBRA 并行数 (默认: jobs//2, 即主 jobs 的一半)')
    g.add_argument('--cobra_linkage_mismatch', type=int, default=2, help='COBRA 链接识别不匹配数 (默认: 2)')
    g.add_argument('--cobra_verbose', action='store_true', help='cobra_pipeline.py 详细日志')

    g = p.add_argument_group('CLUSTER 阶段 (cluster_pipeline.py)')
    g.add_argument('--skip_vclust', action='store_true', help='跳过 vclust 聚类步骤')
    g.add_argument('--vclust_cluster_file', help='复用已有 vclust 聚类 TSV 文件')
    g.add_argument('--skip-rmdup', action='store_true',
                   help='跳过 genome_rmDuplicates 去冗余 (vclust 聚类后)')
    g.add_argument('--rmdup-length', type=int, default=1000,
                   help='genome_rmDuplicates 短序列阈值 bp (默认: 1000)')
    g.add_argument('--skip-flye', action='store_true',
                   help='跳过 Flye 共组装延伸, 仅合并多样本 COBRA 结果')
    g.add_argument('--flye-min-overlap', type=int, default=500,
                   help='Flye 最小重叠长度 bp (默认: 500, 因为输入 contig 通常较短)')
    g.add_argument('--flye-read-error', type=float, default=0.005,
                   help='Flye 读长错误率 (默认: 0.005)')

    g = p.add_argument_group('Taxonomy 阶段 (virus_classifier.py)')
    g.add_argument('--tax_tools', default='all', help='分类工具: genomad,metabuli,diamond_lca,VITAP,mmseqs,ACVirus,vcontact3,PhaGCN3,all (默认: all)')
    g.add_argument('--tax_jobs', type=int, default=1, help='分类并行任务数 (默认: 1)')
    g.add_argument('--tax_ext', default='.fasta', help='分类输入文件扩展名 (默认: .fasta)')
    g.add_argument('--tax_sample_name', default='Votus', help='分类样本名 (输出目录前缀, 默认: Votus)')
    g.add_argument('--tax_remove_suffix', help='分类输入文件去后缀名')

    g = p.add_argument_group('Host 阶段 (run_host_prediction.py)')
    g.add_argument('--skip_rnavirhost', action='store_true', help='跳过 RNAVirHost 宿主预测')
    g.add_argument('--skip_phabox', action='store_true', help='跳过 PhaBOX2 宿主预测')
    g.add_argument('--skip_ictv', action='store_true', help='跳过 ICTV 宿主查找')

    g = p.add_argument_group('数据库路径')
    g.add_argument('--host_db', default='~/database/host_db/', help='宿主数据库根目录')
    g.add_argument('--kraken2_db', help='Kraken2 宿主库 (覆盖 --host_db 自动检测)')
    g.add_argument('--host_align_db', help='宿主比对索引 (覆盖 --host_db 自动检测)')
    g.add_argument('--virus_db', default='~/database/virus-db/', help='病毒鉴定数据库根目录')
    g.add_argument('--checkv_db', default='~/database/virus-db/checkv-db-v1.7', help='CheckV 数据库路径')
    g.add_argument('--blast-db', default='~/database/virus-db/ncbi-virus_ref/ncbi-virus_ref.blast.db', help='BLAST 参考数据库 (rescue 阶段)')

    g = p.add_argument_group('工具与算法')
    g.add_argument('--aligner', default='bowtie2', choices=['bowtie2', 'hisat2', 'minimap2'])
    g.add_argument('--seq_type', default='rna-short', choices=['dna-short', 'rna-short', 'nanopore', 'pacbio'])
    g.add_argument('--rrna', action='store_true', help='开启 rRNA 剔除')
    g.add_argument('--rrna_tool', default='ribodetector', choices=['ribodetector', 'silva'],
                   help='rRNA 剔除工具: ribodetector (默认) / silva (Bowtie2+SILVA)')
    g.add_argument('--silva_index', help='SILVA Bowtie2 索引前缀 (--rrna_tool silva 时必需)')
    g.add_argument('--coassembly', action='store_true', help='Co-assembly 模式: 合并所有样本 reads 进行单次组装')
    g.add_argument('--bbnorm', action='store_true', help='co-assembly 前 BBNorm 归一化 (target=70 mindepth=2)')
    g.add_argument('--assembler', default='penguin', choices=['megahit', 'rnaviralspades', 'penguin', 'all'])
    g.add_argument('--contig-length', '-l', type=int, default=200, help='contig 最小长度 bp (默认 200)')
    g.add_argument('--identify_tools', default='all', help='病毒鉴定工具')
    g.add_argument('--virus_mode', default='strict', choices=['raw', 'filter', 'strict'],
                   help='COBRA 病毒序列来源: raw=原始鉴定, filter=UniProt过滤, strict=严格过滤 (默认: strict)')
    g = p.add_argument_group('鉴定数据库 (virus_identification.py)')
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
    g.add_argument('--host-mode', default='all', choices=['all','ICTV','RNAVirHost','PhaBOX2'], help='宿主预测模式 (默认: all)')
    g.add_argument('--prob-dir', help='ICTV 宿主概率表目录 (host 阶段, 默认: cross_analysis/)')
    g.add_argument('--ref-genomes', nargs='*', help='ICTV/NCBI 参考基因组 FASTA (可多个, CD-HIT 参考引导预聚类)')

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

    g = p.add_argument_group('分类数据库 (virus_classifier.py)')
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
                    stages_help = []
                    k = j + 1
                    while k < len(sys.argv) and sys.argv[k] not in ('-h', '--help'):
                        if sys.argv[k] in STAGE_HELP:
                            stages_help.append(sys.argv[k])
                        k += 1
                    if stages_help:
                        import argparse as _ap
                        always_show = ['path config', 'pipeline control', 'Resources']
                        all_keys = set()
                        for s in stages_help:
                            for kk in STAGE_ARGS.get(s, []):
                                all_keys.add(kk)
                        print()
                        print(f"  {'='*50}")
                        print(f"  Stages: {', '.join(stages_help)}")
                        print(f"  {'='*50}")
                        for s in stages_help:
                            print(f"    {s}: {STAGE_HELP[s]}")
                        print()
                        if all_keys:
                            print("  [Relevant arguments]")
                            print()
                            pp = _build_parser(add_help=False)
                            for grp in pp._action_groups:
                                title = grp.title
                                if any(kk in title for kk in all_keys) or title in always_show:
                                    print(f"  --- {title} ---")
                                    for action in grp._group_actions:
                                        opts = ', '.join(action.option_strings) if action.option_strings else action.dest
                                        if action.help:
                                            h = action.help
                                            if action.default is not None and action.default is not _ap.SUPPRESS and action.default is not False:
                                                if 'default' not in h and '默认' not in h:
                                                    h += f' (default: {action.default})'
                                            print(f"    {opts:42s} {h}")
                                    print()
                        print(f"  [Usage]")
                        print(f"    single:  --output_dir DIR --stage {stages_help[0]}")
                        print(f"    multi:   --output_dir DIR --stage {' '.join(stages_help)}")
                        print(f"    all:     --output_dir DIR --stage all")
                        print()
                        sys.exit(0)
            print_help_and_exit()

    if len(sys.argv) == 1:
        print_help_and_exit()

    return _build_parser(add_help=False).parse_args()


def _load_config(args):
    """加载 YAML 配置, CLI 参数覆盖配置文件值。返回更新后的 args。"""
    config_path = args.config
    if not config_path:
        for p in [SCRIPT_DIR.parent / "pipeline_config.yaml",
                  SCRIPT_DIR / "pipeline_config.yaml",
                  Path.cwd() / "pipeline_config.yaml"]:
            if p.is_file(): config_path = str(p); break
    if not config_path:
        return args  # 无配置文件, 纯 CLI 模式

    try:
        import yaml
    except ImportError:
        print("[WARN] PyYAML 未安装, 跳过配置文件加载 (pip install pyyaml)")
        return args

    with open(config_path, encoding='utf-8') as cf:
        config = yaml.safe_load(cf)

    profiles = config.get('profiles', {})
    profile = profiles.get(args.profile, profiles.get('default', {}))
    if not profile:
        print(f"[WARN] 配置 profile '{args.profile}' 未找到, 使用 CLI 参数")
        return args

    # 数据库路径
    db = profile.get('databases', {})
    for key in ['checkv_db','genomad_db','mmseqs_db','virus_db','uniprot_db',
                'host_db','blast_db','nr_db','db_dir',
                'viralverify_hmm','virsorter_db','metabuli_db','virus_taxid']:
        if getattr(args, key, None) is None and key in db:
            setattr(args, key, db[key])

    # 工具路径
    tools = profile.get('tools', {})
    for key in ['salmon','diamond','ragtag']:
        attr_map = {'salmon': 'salmon_bin', 'diamond': None, 'ragtag': None}
        ak = attr_map.get(key, key)
        if ak and getattr(args, ak, None) is None and key in tools:
            setattr(args, ak, os.path.expanduser(tools[key]))

    # 运行参数 (CLI 默认值与 YAML 值比较, 取较大的)
    rt = profile.get('runtime', {})
    _defaults = {'threads': 20, 'jobs': 2, 'tax_jobs': 1}
    for key in ['threads','jobs','tax_jobs']:
        if key in rt:
            current = getattr(args, key, _defaults[key])
            if current == _defaults[key]:  # CLI 未修改, 用 profile 值
                setattr(args, key, int(rt[key]))

    # assembly
    asm_cfg = profile.get('assembly', {})
    if hasattr(args, 'assembler') and 'assembler' in asm_cfg:
        if getattr(args, 'assembler', 'megahit') == 'megahit':
            setattr(args, 'assembler', asm_cfg['assembler'])

    # identification
    id_cfg = profile.get('identification', {})
    for key in ['virus_mode','blast_mode']:
        if getattr(args, key, None) is None and key in id_cfg:
            setattr(args, key, id_cfg[key])

    # cluster
    cl_cfg = profile.get('cluster', {})
    for key in ['min_length','ani','qcov']:
        if key in cl_cfg:
            current = getattr(args, key, None)
            if current is None:
                setattr(args, key, cl_cfg[key])
    if 'ref_genomes' in cl_cfg and not getattr(args, 'ref_genomes', None):
        setattr(args, 'ref_genomes', cl_cfg['ref_genomes'])

    # cobra
    cobra_cfg = profile.get('cobra', {})
    if 'cobra_jobs' in cobra_cfg and getattr(args, 'cobra_jobs', 0) == 0:
        setattr(args, 'cobra_jobs', int(cobra_cfg['cobra_jobs']))

    # taxonomy
    tax_cfg = profile.get('taxonomy', {})
    if 'tax_sample_name' in tax_cfg and getattr(args, 'tax_sample_name', 'Votus') == 'Votus':
        setattr(args, 'tax_sample_name', tax_cfg['tax_sample_name'])

    # host
    host_cfg = profile.get('host', {})
    if 'host_mode' in host_cfg and getattr(args, 'host_mode', 'all') == 'all':
        setattr(args, 'host_mode', host_cfg['host_mode'])
    if 'host_filter' in host_cfg and getattr(args, 'host_filter', 'Plant') == 'Plant':
        setattr(args, 'host_filter', host_cfg['host_filter'])

    # rescue
    res_cfg = profile.get('rescue', {})
    for key in ['checkv_threshold','max_vsi_samples','min_vsi_len']:
        if key in res_cfg:
            current = getattr(args, key, None)
            default_map = {'checkv_threshold': 90.0, 'max_vsi_samples': 10, 'min_vsi_len': 2000}
            if current is None or current == default_map.get(key):
                setattr(args, key, res_cfg[key])

    return args


def _validate_config(args, logger):
    """验证数据库和工具是否存在, 打印配置摘要, 保存 run_config.json"""
    import json
    checks = []

    # 数据库检查 (BLAST DB 查 .nin 文件, 其余查目录存在)
    db_keys = ['checkv_db','genomad_db','mmseqs_db','virus_db','host_db','blast_db','nr_db']
    for key in db_keys:
        val = getattr(args, key, None)
        if val:
            p = Path(os.path.expanduser(str(val)))
            ok = False
            if key == 'blast_db':
                ok = (p.name + '.nin') in set(os.listdir(str(p.parent))) if p.parent.is_dir() else False
            else:
                ok = p.exists()
            status = '✓' if ok else '✗ MISSING'
            checks.append(('DB', key, str(p), status))

    # 工具检查
    tool_checks = [('salmon', getattr(args, 'salmon_bin', None)),
                   ('diamond', getattr(args, 'virseqimprover_path', None))]
    for name, path in tool_checks:
        if path:
            p = Path(os.path.expanduser(str(path)))
            status = '✓' if p.exists() else '✗ MISSING'
            checks.append(('TOOL', name, str(p), status))

    # 打印摘要
    logger.info("=" * 60)
    logger.info("Configuration Summary")
    logger.info("  Profile: %s", getattr(args, 'profile', 'default'))
    logger.info("  " + "-" * 40)
    for cat, name, val, status in checks:
        logger.info("  [%s] %-20s %s  %s", cat, name, status, val)
    logger.info("  " + "-" * 40)
    missing = [c for c in checks if 'MISSING' in c[3]]
    if missing:
        logger.warning("  %d 个资源未找到 (阶段运行时会报错)", len(missing))
    else:
        logger.info("  所有资源验证通过 ✓")
    logger.info("=" * 60)

    # 保存 run_config.json
    run_cfg = {
        "profile": getattr(args, 'profile', 'default'),
        "stage": getattr(args, 'stage', ['all']),
        "output_dir": str(getattr(args, 'output_dir', '')),
        "threads": getattr(args, 'threads', 20),
        "jobs": getattr(args, 'jobs', 2),
    }
    for key in db_keys:
        run_cfg[key] = str(getattr(args, key, None))
    cfg_path = Path(args.output_dir) / "run_config.json"
    try:
        cfg_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cfg_path, 'w') as cf:
            json.dump(run_cfg, cf, indent=2, ensure_ascii=False)
        logger.info("  run_config.json → %s", cfg_path)
    except: pass


def main():
    args = parse_args()
    args = _load_config(args)
    logger = setup_logger(args.output_dir, level=args.log_level)
    _validate_config(args, logger)

    if getattr(args, 'dump_config', False):
        logger.info("  --dump-config 模式, 仅打印配置摘要")
        return

    stages = set(args.stage)  # 支持多阶段: --stage clean deplete

    logger.info("=" * 50)
    logger.info("Virome Pipeline v2.3")
    logger.info("  Stage:    %s", ','.join(sorted(stages)))
    logger.info("  Output:   %s", args.output_dir)

    # 根据 stages 自动设置 skip 标志
    needs_reads = 'all' in stages or bool(stages & {'clean','deplete','assembly','cobra','rescue'})
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
        args.skip_clean = False
        args.skip_depletion = False
        flow = []
        flow.append('SKIP' if args.skip_clean else 'Clean')
        flow.append('SKIP' if args.skip_depletion else 'Deplete')
        flow.append(f'Assemble({args.assembler})')
        flow.append(f'Identify({args.identify_tools})')
        flow.append('COBRA')
        if getattr(args, 'skip_flye', False):
            flow.append('MergeSamples')
        else:
            flow.append('MergeSamples+Flye')
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

    needs_reads = 'all' in stages or bool(stages & {'clean','deplete','assembly','cobra','rescue'})
    if needs_reads:
        # 优先使用 --input_reads 指定的目录, 但 rescue 阶段推荐 00b_HostDepletion
        if args.input_reads and Path(args.input_reads).exists():
            pipe.reads_dir = Path(args.input_reads).absolute()
            hostdep_dir = Path(args.output_dir) / "00b_HostDepletion"
            if 'rescue' in stages and hostdep_dir.is_dir() and '00b' not in str(pipe.reads_dir):
                logger.warning("  ⚠ rescue 阶段建议 --input_reads out/00b_HostDepletion/ (当前: %s)", pipe.reads_dir)
        elif 'deplete' in stages:
            # deplete: 自动使用 clean 输出 (优先 clumpify, 否则 fasta)
            cl = pipe.d['clean'] / '3.clumpify'
            fa = pipe.d['clean'] / '2.fasta'
            pipe.reads_dir = cl if (cl.exists() and any(cl.iterdir())) else fa
        else:
            pipe.reads_dir = pipe.d['hostdep']
        if not pipe.reads_dir.exists() or not any(pipe.reads_dir.iterdir()):
            if stages <= {'merge', 'cluster', 'taxonomy', 'host', 'checkv', 'report', 'rescue', 'analysis'}:
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
        'bbnorm': pipe.run_bbnorm,
        'assembly': pipe.run_assembly, 'identification': pipe.run_identification,
        'cobra': pipe.run_cobra, 'merge': pipe.run_merge_samples, 'cluster': pipe.run_cluster,
        'taxonomy': pipe.run_taxonomy, 'host': pipe.run_host,
        'checkv': pipe.run_checkv_stage, 'rescue': pipe.run_rescue,
        'analysis': pipe.run_analysis, 'suvtk_annotate': pipe.run_suvtk_annotate, 'report': pipe.run_reports,
    }
    # 按流水线顺序排列
    stage_order = ['clean','deplete','bbnorm','assembly','identification','cobra','merge','cluster',
                   'taxonomy','host','checkv','rescue','analysis','report']
    stages_to_run = [(s, stage_map[s]) for s in stage_order if _all or s in stages]

    # 准备阶段日志目录
    stage_log_dir = Path(args.output_dir) / "10_Reports" / "logs"
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
