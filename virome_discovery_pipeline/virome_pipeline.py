#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Virome Pipeline v2.3 — 宏病毒组端到端全自动主控流水线
===============================================================================

数据流:
  Raw FASTQ ──→ [00a_CleanData] ──→ [00b_HostDepletion] ──→ [01_Assembly]
                                                                   │
   [05_Reports] ←── [04_CLUSTER] ←── [03b_MergeSamples] ←── [03a_COBRA] ←── [02b_Filter] ←── [02a_Identification] ←──┘

依赖脚本 (同目录):
  ../data_preprocessing_pipeline/clean-data.py     → Step 0a: Fastp + Seqkit + Clumpify
  ../data_preprocessing_pipeline/host_depletion.py → Step 0b: Kraken2 + Align + Ribodetector
  assembly_pipeline.py       → Step 1:  MEGAHIT / rnaviralSPAdes / Penguin
  virus_identification.py  → Step 2a: Genomad + Blast + VirSorter2 + ...
  filter_virus.py           → Step 2b: UniProt + CDD 分层过滤
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
import csv
import subprocess
import logging
import re
import json
import shutil
from pathlib import Path
from collections import defaultdict

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import (build_discovery_dirs, normalize_layout_env,
                                   dir_name as layout_dir_name)

import polars as pl
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent

# ── 05 分类阶段: 科-属闸门守卫 ─────────────────────────────────────────
# R 共识脚本每次跑完会写两份台账到集成目录:
#   taxonomy_gate_stamp.tsv  (gate_version / script_md5 / ref_md5 / rows / timestamp)
#   taxonomy_gate_check.tsv  (check / n / expected / status, 其中 dual_ref_conflict 必须 PASS)
# 只有「版本戳达标 且 自检 PASS」的产物才能被 SKIP 复用。否则强制重跑 R 共识。
# 背景: 闸门逻辑上线前的旧产物同样叫 final_integrated_classification.tsv,
#       单看文件存在就 SKIP 会让属服从科 / 逐级相容性两道闸门永久失效 —— 历史错误的规模来源。
TAX_GATE_REQUIRED = "6.2"
TAX_GATE_REF_DEFAULT = "~/database/taxonomy/genus_family_ref.tsv"
TAX_GATE_STAMP = "taxonomy_gate_stamp.tsv"
TAX_GATE_CHECK = "taxonomy_gate_check.tsv"


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


def _count_fasta(path):
    """fasta 序列条数 (读不动或不存在时返回 '?')。"""
    try:
        return sum(1 for line in open(path, encoding='utf-8', errors='replace')
                   if line.startswith('>'))
    except OSError:
        return '?'


def tsv_has_col(path, col):
    """TSV 表头是否含指定列。用于判断旧产物是否需按新表头重建。"""
    try:
        with open(path, encoding='utf-8', errors='replace') as f:
            return col in f.readline().split('\t')
    except OSError:
        return False


def run_cmd(cmd, logger, step_name, log_file=None, env=None, cwd=None):
    """执行 shell 命令。stdout/stderr 同时输出到终端、主编排日志、可选的阶段日志。

    返回 (ok, detail)。env/cwd 用于需要在指定目录并注入环境变量的阶段
    (如独立后运行脚本用变量驱动, 见各阶段 --help)。
    """
    logger.info("[%s] 执行中...", step_name)
    logger.debug("  CMD: %s", cmd)
    lf = open(log_file, 'a', encoding='utf-8') if log_file else None
    try:
        if lf:
            lf.write(f"=== {step_name} ===\nCMD: {cmd}\n\n")
        with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, bufsize=1, cwd=cwd,
                              env={**os.environ, 'PYTHONUNBUFFERED': '1', **(env or {})}) as proc:
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
                        reports_log = Path(h.baseFilename).resolve().parent / layout_dir_name('d_reports') / "logs"
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
        # read 编号 _1/_2: 只匹配下划线前缀 (样本名里的 -1/-2 连字符不匹配)
        r'(?<![.-])[._][12](?=[._]|\.|$)',
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

        # 布局解析 (CLI > MMPV_IO_LAYOUT > legacy) 并写回环境, 子进程继承同一布局;
        # self.d 由 mmpv_common.io_layout 构造, legacy 分支与 v3.0 字面量逐键一致
        self.layout = normalize_layout_env(getattr(args, 'io_layout', None))
        self.d = build_discovery_dirs(out, raw, self.layout)
        # 仅创建根目录, 各阶段按需创建自己的子目录
        self.d['root'].mkdir(parents=True, exist_ok=True)

        # 脚本路径（质控/去宿主已独立为顶层管线目录）
        self.sc = {
            'clean':    script_dir.parent / 'data_preprocessing_pipeline' / 'clean-data.py',
            'deplete':  script_dir.parent / 'data_preprocessing_pipeline' / 'host_depletion.py',
            'assembly': script_dir / 'assembly_pipeline.py',
            'identify': script_dir / 'virus_identification.py',
            'filter':   script_dir / 'filter_virus.py',
            'cobra':    script_dir / 'cobra_pipeline.py',
            'cluster':  script_dir / 'cluster_pipeline.py',
            'rescue':   script_dir / 'rescue_pipeline.py',
            'host_pred': script_dir / 'run_host_prediction.py',
            'classifier': script_dir / 'virus_classifier.py',
            'classifier_R': script_dir / 'virus_classifier_analysis.R',
            'c9': script_dir / 'utils/classify_contigs.py',
        }

        # 数据流指针
        self.reads_dir = self.d['raw']
        # co-assembly 合并 reads 目录 (run_assembly 内赋值; 预置 None 防止
        # 单独续跑 --stage rescue 时 getattr 回退到逐样本 reads_dir)
        self.coassembly_merged_dir = None

        # 验证
        self._validate()
        self._check_db_requirements()

        # 检测原始样本
        no_reads_stages = {'identification', 'filter', 'cobra', 'merge', 'cluster', 'taxonomy', 'host', 'checkv', 'report', 'rescue', 'analysis', 'analysis_verify'}
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

    def _check_db_requirements(self):
        """按 stage 检查必需数据库参数, 缺失即 fail-fast。(2026-09-04 修复: 原逻辑为死代码)"""
        stages = set(self.args.stage)
        _all = 'all' in stages
        doing_clean = _all or 'clean' in stages
        need_virus_db = _all or bool(stages & {'identification', 'taxonomy', 'host'})
        need_checkv_db = _all or bool(stages & {'checkv', 'rescue'})

        # deplete 阶段需要 kraken2 (host_align 仅在 --deplete_steps 含 align 时必需)
        need_deplete_db = _all or 'deplete' in stages or ('clean' in stages and not self.args.skip_depletion)
        need_align = need_deplete_db and 'align' in self.args.deplete_steps.split(',')
        if need_deplete_db and not self.args.host_db:
            if not self.args.kraken2_db or (need_align and not self.args.host_align_db):
                self.log.error("致命: deplete 阶段需要 --kraken2_db (--deplete_steps 含 align 时还需 --host_align_db) 或 --host_db")
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

    def _stage_marker(self, stage):
        """返回阶段完成标记文件路径"""
        return self.d['root'] / f".{stage}.ok"

    def _stage_done(self, stage):
        """标记阶段已完成"""
        self._stage_marker(stage).write_text("ok")
        self.log.debug("  标记 %s 完成", stage)

    def _stage_skip(self, stage):
        """检查是否跳过 (已完成且非 force)"""
        if self.args.force:
            self._stage_marker(stage).unlink(missing_ok=True)
            return False
        if not self._stage_marker(stage).is_file():
            return False
        # 标记存在, 但验证关键产物也真实存在 (防止空跑标记)
        if not self._stage_validated(stage):
            self.log.warning("[%s] .oks 标记存在但产物为空, 删除标记重新执行", stage)
            self._stage_marker(stage).unlink(missing_ok=True)
            return False
        self.log.info("[%s] 跳过 (已完成, .%s.ok 存在)",
                      STAGE_HELP.get(stage, stage).split(':')[0], stage)
        return True

    def _stage_validated(self, stage):
        """验证阶段关键产物存在 (防止空跑标记)。返回 True=产物存在, False=空跑。"""
        import glob as _vg
        d = self.d
        checks = {
            'clean':    lambda: d['clean'].is_dir() and any(d['clean'].iterdir()),
            'deplete':  lambda: d['hostdep'].is_dir() and any(d['hostdep'].iterdir()),
            'assembly': lambda: len(_vg.glob(os.path.join(str(d['asm']), '*', '*.contig.fasta'))) > 0,
            'identification': lambda: len(_vg.glob(os.path.join(str(d['ident']), '*', '*virus.all.candidate.fasta'))) > 0,
            'filter':   lambda: any(
                p.is_dir() and p.name not in ('cdd','uniprot')
                for p in d['filter'].iterdir()
            ) if d['filter'].is_dir() else False,
            'cobra':    lambda: (
                # co-assembly 模式下 COBRA 被有意跳过, 信任标记
                getattr(self.args, 'coassembly', False)
                or any(
                    (cobra_d and cobra_d.is_dir() and any(
                        p.is_dir() for p in cobra_d.iterdir()
                        if p.name not in ('cobra.log','checkpoint_status.json')
                    ))
                    for cobra_d in d.get('_cobra_dirs', [d['cobra']])
                )
            ),
            'merge':    lambda: (
                (d['merge_samples'] / 'all_sample_virus_combined.fasta').exists()
                or ((d['merge_samples'] / 'all_sample_virus.fasta').is_file()
                    and (d['merge_samples'] / 'all_sample_virus.fasta').stat().st_size > 0)
            ),
            'cluster':  lambda: (
                (d['centroids'] / 'final_centroids.fasta').is_file()
                or (d.get('_centroids_v1') and (d['_centroids_v1'] / 'final_centroids.fasta').is_file())
                or (d.get('_centroids_v2') and (d['_centroids_v2'] / 'final_centroids.fasta').is_file())
            ),
            'taxonomy': lambda: (
                (d['taxonomy'] / 'Votus.integrated' / 'final_integrated_classification.tsv').is_file()
                and (d['taxonomy'] / 'Votus.integrated' / 'final_integrated_classification.tsv').stat().st_size > 100
            ),
            'host':     lambda: (
                (d['host_pred'] / 'ensemble_host_summary.tsv').is_file()
                and (d['host_pred'] / 'ensemble_host_summary.tsv').stat().st_size > 100
            ),
            'checkv':   lambda: len(_vg.glob(str(d['checkv_dir'] / '*' / 'completeness.tsv'))) > 0,
            'rescue':   lambda: d['rescue_dir'].is_dir() and any(d['rescue_dir'].iterdir()),
        }
        if stage in checks:
            return checks[stage]()
        # 未列入验证的阶段 (如 analysis/report) 信任标记
        return True


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
        if not host_align_db and 'align' in self.args.deplete_steps.split(','):
            self.log.error("致命: --deplete_steps 含 align 时需要 --host_align_db 或 --host_db")
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
        ]
        if host_align_db:
            parts.append(f"--step2_index {host_align_db}")
        parts += [            f"--input-dir {self.reads_dir}",
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
        if self.args.vmtouch:
            parts.append("--vmtouch")
        if self.args.keep_taxids:
            parts.append(f"--keep-taxids {self.args.keep_taxids}")
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

        # Co-assembly 模式: 合并所有样本 reads → 单次组装
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
        asm_tmp = self.args.asm_tmp_dir or str(self.d['asm'] / "tmp")
        parts.append(f"--tmp-dir {asm_tmp}")
        if self.args.asm_keep_temp:
            parts.append("--keep-temp")

        ok, _ = run_cmd(' '.join(parts), self.log, "Assembly", str(self.d['asm'] / "assembly.log"))
        if not ok:
            self.log.error("组装失败, 终止。")
            sys.exit(1)

        self.asm_map = scan_contig_files(self.d['asm'], asm_tools)
        self.log.info("  组装完成: %d 样本有 contig 输出", len(self.asm_map))
        # 清理临时目录
        asm_tmp_path = Path(asm_tmp)
        if not self.args.asm_keep_temp and asm_tmp_path.exists():
            shutil.rmtree(asm_tmp_path, ignore_errors=True)
            self.log.info("  临时目录已清理: %s", asm_tmp)

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
        parts.append("--skip_uniprot_filter")  # UniProt 过滤移至 02b_Filter
        if self.args.skip_nr_filter:
            parts.append("--skip_nr_filter")
        if self.args.skip_id_plots:
            parts.append("--skip_plots")
        if self.args.clean_failed:
            parts.append("--clean_failed")
        run_cmd(' '.join(parts), self.log, "VirusIdentification", str(self.d['ident'] / "ident.log"))

        self.viral_map = scan_viral_files(self.d['ident'])
        self.log.info("  鉴定完成: %d 样本有病毒候选序列", len(self.viral_map))

    # ── Step 2b: UniProt + CDD 分层过滤 ──
    def run_filter(self):
        self.d['filter'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[2b] UniProt + CDD tiered filter")

        # 收集鉴定产出的 candidate fasta
        # 注意: 优先用 *_virus.all.candidate.fasta (工具鉴定的真病毒候选),
        #       避免用 *.virus.candidate.fasta (纯长度过滤, 含大量非病毒序列)
        candidates = []
        for sample_dir in sorted(self.d['ident'].iterdir()):
            if not sample_dir.is_dir(): continue
            sample = sample_dir.name
            cand = None
            for pat in [f"{sample}_virus.all.candidate.fasta",
                        f"*_virus.all.candidate.fasta",
                        f"{sample}.virus.candidate.fasta"]:
                hits = list(sample_dir.glob(pat))
                if hits and hits[0].stat().st_size > 50:
                    cand = hits[0]
                    break
            if cand:
                candidates.append((sample, str(cand)))

        if not candidates:
            self.log.warning("  未找到 virus.candidate.fasta，跳过过滤")
            return

        self.log.info("  %d 个样本待过滤", len(candidates))

        # 默认合并运行: 合并所有 candidate -> 一次性过滤 -> 拆分回样本
        if not getattr(self.args, 'no_merge', False) and len(candidates) > 1:
            self._run_filter_merged(candidates)
        else:
            self._run_filter_per_sample(candidates)

        self.log.info("  过滤完成, 输出: %s", self.d['filter'])

    def _run_filter_merged(self, candidates):
        """合并模式: 所有样本 candidate 合并 → 一次过滤 → 拆分"""
        merged_fa = self.d['filter'] / 'all_candidates.fasta'
        # 合并时总线程 = threads * jobs
        total_threads = min(self.args.threads * max(1, getattr(self.args, 'jobs', 1)), 64)
        self.log.info("  合并模式: %d 样本 → %s (%d 线程)", len(candidates), merged_fa, total_threads)

        # 合并
        with open(merged_fa, 'w') as fout:
            for sample, cand_fa in candidates:
                with open(cand_fa) as fin:
                    fout.write(fin.read())

        # 构建参数
        parts = [
            f"python {self.sc['filter']}",
            f"-m {getattr(self.args, 'filter_mode', 'filter')}",
            f"-i {merged_fa}",
            f"-o {self.d['filter']}",
            f"-t {total_threads}",
            f"--split-by-prefix",
        ]
        uniprot_db = getattr(self.args, 'uniprot_db', None)
        if uniprot_db:
            parts.append(f"--uniprot-db {uniprot_db}")
        virus_taxid = getattr(self.args, 'virus_taxid', None)
        if virus_taxid:
            parts.append(f"--virus-taxid {virus_taxid}")

        run_cmd(' '.join(parts), self.log, "Filter-merged",
                str(self.d['filter'] / "filter_merged.log"))

        # 清理合并文件
        merged_fa.unlink(missing_ok=True)

    def _run_filter_per_sample(self, candidates):
        """逐样本模式: 每样本独立运行 filter_virus.py"""

        # 构建基础参数
        parts_base = [
            f"python {self.sc['filter']}",
            f"-m {getattr(self.args, 'filter_mode', 'filter')}",
        ]
        uniprot_db = getattr(self.args, 'uniprot_db', None)
        if uniprot_db:
            parts_base.append(f"--uniprot-db {uniprot_db}")
        virus_taxid = getattr(self.args, 'virus_taxid', None)
        if virus_taxid:
            parts_base.append(f"--virus-taxid {virus_taxid}")

        # 并行处理: 每样本 threads 线程, jobs 个样本并行, 支持样本级断点续传
        jobs = max(1, getattr(self.args, 'jobs', 1))
        self.log.info("  并行: %d 样本, 每样本 %d 线程", jobs, self.args.threads)

        from concurrent.futures import ThreadPoolExecutor, as_completed
        try: from tqdm import tqdm
        except ImportError: tqdm = lambda x, **kw: x

        # 断点续传: 跳过已完成的样本
        pending = []
        skipped = 0
        for sample, cand_fa in candidates:
            sample_out = self.d['filter'] / sample
            final_fa = sample_out / f"{sample}.virus.candidate_filtered.fasta"
            if not self.args.force and final_fa.is_file():
                self.log.info("  Filter-%s: 跳过 (已完成, %s 存在)", sample, final_fa.name)
                skipped += 1
                continue
            pending.append((sample, cand_fa))

        if skipped:
            self.log.info("  断点续传: 跳过 %d 已完成, 待处理 %d", skipped, len(pending))
        if not pending:
            self.log.info("  全部已完成, 跳过过滤")
            return

        def filter_one(sample, cand_fa):
            sample_out = self.d['filter'] / sample
            os.makedirs(sample_out, exist_ok=True)
            parts = parts_base + [
                f"-i {cand_fa}",
                f"-o {sample_out}",
                f"-t {self.args.threads}",
            ]
            log_file = sample_out / "filter.log"
            cmd = ' '.join(parts)
            result = subprocess.run(cmd, shell=True, executable='/bin/bash',
                                    capture_output=True, text=True)
            log_file.write_text(result.stdout + '\n' + result.stderr)
            return sample, result.returncode

        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futures = {ex.submit(filter_one, s, f): s for s, f in pending}
            pbar = tqdm(as_completed(futures), total=len(pending), desc="  过滤", unit="sample")
            for fut in pbar:
                sample, rc = fut.result()
                status = "OK" if rc == 0 else f"FAIL({rc})"
                pbar.set_postfix_str(f"{sample[:30]} {status}")
                self.log.info("  Filter-%s: %s", sample, status)

        self.log.info("  过滤完成, 输出: %s", self.d['filter'])

    # ── Step 3a: COBRA 延伸 ──
    # ── 自动探测 COBRA 目录 (兼容 03a_COBRA / 03_COBRA) ──
    def _resolve_cobra_dir(self):
        for d in self.d.get('_cobra_dirs', [self.d['root'] / layout_dir_name('d_cobra', self.layout)]):
            if d.is_dir():
                self.d['cobra'] = d
                return d
        self.d['cobra'] = self.d['root'] / layout_dir_name('d_cobra', self.layout)
        return self.d['cobra']

    # ── Step 3a: COBRA 延伸 ──
    def run_cobra(self):
        self._resolve_cobra_dir()
        # co-assembly 模式: reads 已合并为单一 ALL_merged, 无逐样本延伸意义,
        # 直接跳过 (merge 阶段会改从 02a/02b 收集病毒候选作为合并输入)
        if getattr(self.args, 'coassembly', False):
            self.log.info("=" * 50)
            self.log.info("[3a] COBRA 批量延伸 — [co-assembly] 跳过逐样本延伸")
            self.d['cobra'].mkdir(parents=True, exist_ok=True)
            return
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

        # 病毒候选来源: 逐级选择 (cdd > uniprot > raw)
        virus_mode = getattr(self.args, 'virus_mode', 'cdd')
        if virus_mode == 'raw':
            virsorter_src = self.d['ident']
        elif virus_mode == 'uniprot':
            virsorter_src = self.d['filter']
        else:  # cdd (default)
            virsorter_src = self.d['filter'] if self.d['filter'].is_dir() else self.d['ident']
        self.log.info("  病毒候选来源 (%s): %s", virus_mode, virsorter_src)
        # COBRA 内部只用 raw 模式 (我们的目录结构已解决文件查找问题)
        cobra_virus_mode = 'raw'

        parts = [
            f"python {self.sc['cobra']}",
            f"--mode mix",
            f"--reads-dir {self.reads_dir}",
            f"--contigs-dir {self.d['asm']}",
            f"--virsorter-dir {virsorter_src}",
            f"--output-dir {self.d['cobra']}",
            f"--assembly-tools {','.join(asm_tools)}",
            f"--virus-mode {cobra_virus_mode}",
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

        # co-assembly 模式: 无 COBRA 产物, 改从鉴定/过滤结果收集伪样本
        # ALL_merged 的病毒候选 FASTA 作为合并输入, 后续 Flye 共组装照常
        if getattr(self.args, 'coassembly', False):
            merged_fa = merge_dir / "all_sample_virus.fasta"
            candidates = []
            for src in [self.d['filter'], self.d['ident']]:
                if not src.is_dir():
                    continue
                d_all = src / 'ALL_merged'
                if d_all.is_dir():
                    candidates = sorted(d_all.rglob('*virus*.fasta')) + \
                                 sorted(d_all.rglob('*virus*.fa'))
                    if candidates:
                        break
            if not candidates:
                self.log.error("[co-assembly] 未在 %s / %s 下找到 ALL_merged 病毒候选", self.d['filter'], self.d['ident'])
                sys.exit(1)
            n_in = 0
            with open(merged_fa, 'w') as out:
                for cf in candidates:
                    out.write(cf.read_text())
                    n_in += 1
            total_bp = sum(len(rec.seq) for rec in SeqIO.parse(str(merged_fa), "fasta"))
            self.log.info("=" * 50)
            self.log.info("[3b] [co-assembly] 合并 ALL_merged 病毒候选: %d 文件, %.1f Mb → %s", n_in, total_bp / 1e6, merged_fa)
            return self._merge_flye_stage(merged_fa)

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

        return self._merge_flye_stage(merged_fa)


    # ── Step 3b 内部: Flye 共组装段 (per-sample 与 co-assembly 两路共用) ──
    def _merge_flye_stage(self, merged_fa):
        merge_dir = self.d['merge_samples']
        # ── Flye 共组装 (可选) ──
        if getattr(self.args, 'skip_flye', False):
            self.log.info("  merge 阶段完成 (跳过 Flye)")
            return

        flye_dir = merge_dir / "flye_coassembly"
        flye_dir.mkdir(parents=True, exist_ok=True)
        flye_out = flye_dir / "flye_output"

        min_ovlp = getattr(self.args, 'flye_min_overlap', 1000)
        read_err = getattr(self.args, 'flye_read_error', 0.005)

        self.log.info("-" * 40)
        self.log.info("  Flye --subassemblies: min_overlap=%d, read_error=%.3f", min_ovlp, read_err)

        flye_cmd = (
            f"flye --subassemblies {merged_fa}"
            f" -t {min(max(self.args.threads * self.args.jobs, self.args.threads), 128)}"
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
        mapping_tsv = flye_dir / 'flye_mapping.tsv'
        sample_json = flye_dir / 'flye_sample_map.json'
        if trace_script.exists():
            run_cmd(
                f"python {trace_script} -i {flye_out} -o {mapping_tsv} -s {flye_dir / 'flye_summary.txt'} --full",
                self.log, "Flye trace", str(flye_dir / "trace.log")
            )
        # TSV → JSON: Flye contig → 来源样本 (供 rescue 阶段精确溯源)
        if mapping_tsv.is_file():
            flye_map = {}
            with open(mapping_tsv) as f:
                header = f.readline().strip().split('\t')
                try:
                    input_col = header.index("input_reads")
                except ValueError:
                    input_col = None
                for line in f:
                    cols = line.strip().split('\t')
                    if input_col is not None and len(cols) > input_col:
                        contig_id = cols[0]  # output_contig
                        raw_reads = cols[input_col]
                        if raw_reads:
                            # 从 contig 名提取样本: ERR2040118_clean_NODE_2 → ERR2040118_clean
                            samples = set()
                            for rid in re.split(r'[;,]', raw_reads):
                                rid = rid.strip()
                                m = re.match(r'^(\w+?_\w+?)_', rid)  # 取前两段
                                if m:
                                    samples.add(m.group(1))
                            if samples:
                                flye_map[contig_id] = sorted(samples)
            if flye_map:
                json.dump(flye_map, open(sample_json, 'w'), ensure_ascii=False, indent=2)
                self.log.info("  Flye 样本映射: %d 条 → %s", len(flye_map), sample_json)

        # 合并: Flye + 原始 → combined
        combined_fa = merge_dir / "all_sample_virus_combined.fasta"
        with open(combined_fa, 'w') as out:
            with open(flye_assembly) as inf:
                out.write(inf.read())
            with open(merged_fa) as inf:
                out.write(inf.read())

        n_combined = sum(1 for _ in open(combined_fa) if _.startswith('>'))
        combined_bp = sum(len(rec.seq) for rec in SeqIO.parse(str(combined_fa), "fasta"))
        self.log.info("  合并: Flye %d + 原始合并文件 → %d 条, %.1f Mb → %s",
                      n_flye, n_combined, combined_bp / 1e6, combined_fa)
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
                # 向后兼容: 现场收集 (co-assembly 模式则从鉴定/过滤结果收集 ALL_merged)
                if getattr(self.args, 'coassembly', False):
                    self.log.info("  [co-assembly] 未找到 merge 输出, 从鉴定/过滤结果现场收集...")
                    cluster_fa = Path(self.d['root']) / "cluster_input.fasta"
                    n_in = 0
                    with open(cluster_fa, 'w') as out:
                        for src in [self.d['filter'], self.d['ident']]:
                            d_all = src / 'ALL_merged' if hasattr(src, '__truediv__') else None
                            if d_all is None or not d_all.is_dir():
                                continue
                            for cf in sorted(d_all.rglob('*virus*.fasta')) + sorted(d_all.rglob('*virus*.fa')):
                                with open(cf) as inf:
                                    out.write(inf.read()); n_in += 1
                            if n_in:
                                break
                    if cluster_fa.exists() and cluster_fa.stat().st_size > 0:
                        self.log.info("  CLUSTER 输入 ([co-assembly] 现场收集): %d 文件 → %s", n_in, cluster_fa)
                    else:
                        self.log.error("[co-assembly] 未找到 ALL_merged 病毒候选, 无法进行聚类")
                        return
                else:
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
        # 线程预算: 单进程阶段用满 threads*jobs (总预算), taxonomy 按 7 工具均分
        _budget = max(self.args.threads * self.args.jobs, self.args.threads)
        parts = [
            f"python {self.sc['cluster']}",
            f"-i {cluster_fa}",
            f"-o {self.d['cluster']}",
            f"-t {_budget}",
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
            for alt_key in ['_centroids_v1', '_centroids_v2']:
                centroids_fa = self.d[alt_key] / "final_centroids.fasta"
                if centroids_fa.is_file():
                    break
        if not centroids_fa.is_file():
            self.log.error("centroids 不存在 (试了 4_centroids, 04_centroids, 4.centroids), 请先运行 --stage cluster")
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

        # 将预通过的 CheckV pass IDs 写入 rescue 目录 (供分支A报告使用)
        checkv_pass_plant = target_ids & checkv_pass_ids - cdhit_known_ids
        if checkv_pass_plant:
            rescue_dir = self.d['rescue_dir'] / "_".join(host_filter) if host_filter else self.d['rescue_dir'] / "Plant"
            rescue_dir.mkdir(parents=True, exist_ok=True)
            pre_pass_file = rescue_dir / "branch_a_prepass.txt"
            with open(pre_pass_file, "w") as ppf:
                for cid in sorted(checkv_pass_plant):
                    ppf.write(f"{cid}\n")

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
        # Flye 样本映射 (让 rescue 精确找到 Flye contig 的来源样本)
        flye_map = self.d['merge_samples'] / "flye_coassembly" / "flye_sample_map.json"
        if flye_map.is_file():
            parts.append(f"--flye-sample-map {flye_map}")

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
        all_plant = self.d['rescue_dir'] / "HQ_plant_viruses.fasta"
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
                f"-d {self.args.checkv_db} -t {min(self.args.threads * self.args.jobs, 16)} "
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

    # ── 05 闸门守卫: 判断已有产物是否带达标闸门版本戳 + 自检 PASS ──
    def _tax_gate_state(self, int_dir):
        """返回 (ok, why)。ok=False 时调用方必须重跑 R 共识，不得 SKIP。

        判据: 版本戳存在且 gate_version >= TAX_GATE_REQUIRED，自检台账里
        dual_ref_conflict == PASS。任何缺失/降级都算不达标（含参照表缺失的
        SKIP_REF_MISSING，那种情况下闸门实际没跑）。
        """
        def _ver(v):
            try:
                return tuple(int(x) for x in str(v).strip().split("."))
            except Exception:
                return (0,)

        stamp_path = int_dir / TAX_GATE_STAMP
        check_path = int_dir / TAX_GATE_CHECK
        if not stamp_path.is_file():
            return False, f"缺闸门版本戳 {TAX_GATE_STAMP} (旧产物, 闸门未跑过)"

        stamp_ver = ""
        try:
            lines = stamp_path.read_text(encoding="utf-8", errors="replace").splitlines()
        except Exception as exc:
            return False, f"闸门版本戳不可读: {exc}"
        for line in lines:
            fields = line.split("\t")
            if len(fields) >= 2 and fields[0].strip() == "gate_version":
                stamp_ver = fields[1].strip()
                break
        if not stamp_ver:
            return False, "闸门版本戳无 gate_version 字段"
        if _ver(stamp_ver) < _ver(TAX_GATE_REQUIRED):
            return False, f"闸门版本过低 {stamp_ver} < {TAX_GATE_REQUIRED}"

        if not check_path.is_file():
            return False, f"缺自检台账 {TAX_GATE_CHECK}"
        status = ""
        try:
            for line in check_path.read_text(encoding="utf-8", errors="replace").splitlines():
                fields = line.split("\t")
                if len(fields) >= 4 and fields[0].strip() == "dual_ref_conflict":
                    status = fields[3].strip()
        except Exception as exc:
            return False, f"自检台账不可读: {exc}"
        if status != "PASS":
            return False, f"科-属自检未通过 dual_ref_conflict={status or 'missing'}"
        return True, f"闸门 {stamp_ver} 达标且自检 PASS"

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
            # 线程预算: taxonomy 7 工具并行, 每工具均分 threads*jobs 总预算
            _tax_budget = max((self.args.threads * self.args.jobs) // 7, 1)
            parts = [
                f"python {self.sc['classifier']}",
                f"-g {centroids}",
                f"-s {sample}",
                f"-t {self.args.tax_tools}",
                f"-o {out_subdir}",
                f"-p {_tax_budget}",
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
        # 闸门守卫: 产物不仅要存在, 还要带 >= TAX_GATE_REQUIRED 的闸门版本戳且自检 PASS,
        # 否则旧产物会让「属服从科 + 逐级相容性」两道闸门永久失效。
        gate_ok, gate_why = self._tax_gate_state(int_dir)
        if not self.args.force and gate_ok and final_tsv.is_file() and final_tsv.stat().st_size > 100:
            self.log.info("  [SKIP] R consensus — 已有结果 (%s)", gate_why)
        elif combined_tsv.is_file():
            if final_tsv.is_file() and not gate_ok:
                self.log.warning("  重跑 R consensus: %s", gate_why)
            int_dir.mkdir(parents=True, exist_ok=True)
            ref_path = os.environ.get("MMPV_GENUS_FAMILY_REF") or os.path.expanduser(TAX_GATE_REF_DEFAULT)
            if not os.path.isfile(ref_path):
                self.log.error("  科-属参照表不存在: %s —— 闸门会降级为 SKIP, 该产物不可直接出货", ref_path)
            parts = [
                f"cd {int_dir} &&",
                f"MMPV_GENUS_FAMILY_REF={ref_path}",
                "Rscript", str(self.sc['classifier_R']),
                "--combined", str(combined_tsv),
                "--output", str(int_dir),
            ]
            ok, _ = run_cmd(' '.join(parts), self.log, "R consensus", str(self.d['taxonomy'] / "r_consensus.log"))
            if not ok:
                self.log.warning("  R consensus 失败")
            # 跑完读自检台账: 不通过必须 fail loud, 不允许静默出货
            post_ok, post_why = self._tax_gate_state(int_dir)
            if post_ok:
                self.log.info("  R 共识闸门自检: %s", post_why)
            else:
                self.log.error("  R 共识闸门未通过: %s —— 该分类表带科属冲突, 不得作为合格产物出货", post_why)

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
        # 线程预算: host 单进程, 用满 threads*jobs 总预算
        parts = [
            f"python {self.sc['host_pred']}",
            f"-i {centroids}",
            f"--tax {tax_tsv}",
            f"-o {self.d['host_pred']}",
            f"-t {max(self.args.threads * self.args.jobs, self.args.threads)}",
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
                       f"-d {self.args.checkv_db} -t {min(self.args.threads * self.args.jobs, 16)}")
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
        """[9] 病毒基因组下游分析 + 结果整合

        子目录:
          HQ_analysis/        - HQ 植物病毒信息 (05_Taxonomy)
          all_plant_analysis/ - 全部植物病毒物种级聚合
          viroid_analysis/    - 类病毒批量分析
        """
        self.d['analysis'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[9] Virome Analysis — 下游分析 + 结果整合")

        all_plant = self.d['rescue_dir'] / "HQ_plant_viruses.fasta"
        if not all_plant.is_file():
            self.log.warning("  HQ_plant_viruses.fasta 不存在, 跳过")
            return

        import shutil
        dest_fasta = self.d['analysis'] / "HQ_plant_viruses.fasta"
        if not dest_fasta.exists():
            shutil.copy(all_plant, dest_fasta)

        # ── 软链接: 统一入口文件 ──
        # All_plant.viruses.fasta → 06 宿主分类的全部植物病毒
        plant_classified = self.d['host_pred'] / "host_classified_fasta" / "Plant.classified.fasta"
        all_plant_link = self.d['analysis'] / "All_plant.viruses.fasta"
        if plant_classified.is_file() and not all_plant_link.exists():
            all_plant_link.symlink_to(os.path.relpath(plant_classified, self.d['analysis']))
            self.log.info("  链接: All_plant.viruses.fasta → %s", plant_classified.name)
        # Viroid.all.fasta → 02 鉴定产出的全部类病毒候选
        ident_dir = self.d['ident']
        viroid_all = self.d['analysis'] / "Viroid.all.fasta"
        if ident_dir.is_dir() and not viroid_all.exists():
            viroid_files = list(ident_dir.rglob('*viroids.candidate.fasta'))
            if viroid_files:
                viroid_all.symlink_to(os.path.relpath(viroid_files[0], self.d['analysis']))
                self.log.info("  链接: Viroid.all.fasta → 02a_Identification")

        n = sum(1 for _ in open(all_plant) if _.startswith('>'))
        self.log.info("  输入: %d 条植物病毒", n)

        logs_dir = self.d['analysis'] / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)

        # 子目录
        hq_dir = self.d['analysis'] / "HQ_analysis"; hq_dir.mkdir(parents=True, exist_ok=True)
        plant_dir = self.d['analysis'] / "all_plant_analysis"; plant_dir.mkdir(parents=True, exist_ok=True)
        viroid_dir = self.d['analysis'] / "viroid_analysis"; viroid_dir.mkdir(parents=True, exist_ok=True)

        t = getattr(self.args, 'threads', 40)
        sample_name = getattr(self.args, 'tax_sample_name', 'Votus')
        consensus_tax = self.d['taxonomy'] / f"{sample_name}.integrated" / "final_integrated_classification.tsv"

        # ── 核酸类型 (DNA/RNA) 标注 → 就地补一列 Nucleic_acid 进分类表 ──
        # 依据 VMR 的 Genome 列按谱系判定, 详见 utils/annotate_nucleic_acid.py
        # 幂等: 重跑先剔旧列再追加; 首次自动备份 .bak_nucleic_<日期>
        nucleic_script = str(SCRIPT_DIR / "utils" / "annotate_nucleic_acid.py")
        if not os.path.isfile(nucleic_script):
            self.log.info("  核酸类型标注 — 跳过 (annotate_nucleic_acid.py 不存在)")
        elif not consensus_tax.is_file():
            self.log.info("  核酸类型标注 — 跳过 (缺 final_integrated_classification.tsv)")
        else:
            ncmd = f"{sys.executable} {nucleic_script} --tax {consensus_tax} --inplace"
            vdb = getattr(self.args, 'virus_db', None)
            if vdb:
                ncmd += f" --virus-db {vdb}"
            ok_n, msg_n = run_cmd(ncmd, self.log, "Nucleic Acid Annotation",
                                  str(logs_dir / "nucleic_acid.log"))
            if not ok_n:
                self.log.warning("  核酸类型标注失败 (%s), 继续后续分析", msg_n)

        # ── Part 1: HQ_plant_viruses_info.tsv → HQ_analysis/ ──
        ref_out = hq_dir / "HQ_plant_viruses_info.tsv"
        if (not ref_out.is_file() or ref_out.stat().st_size < 50
                or not tsv_has_col(ref_out, "Nucleic_acid")):
            self._write_ref_info(dest_fasta, ref_out, consensus_tax)
        else:
            self.log.info("  HQ_plant_viruses_info.tsv — 已存在, 跳过")

        # ── Part 2: All_plant.viruses_info.tsv (物种级聚合) → all_plant_analysis/ ──
        plant_info_out = plant_dir / "All_plant.viruses_info.tsv"
        if (not plant_info_out.is_file() or plant_info_out.stat().st_size < 50
                or not tsv_has_col(plant_info_out, "Nucleic_acid")):
            plant_info_script = str(SCRIPT_DIR / "utils" / "build_plant_virus_info.py")
            plant_fa = self.d['host_pred'] / "host_classified_fasta" / "Plant.classified.fasta"
            if os.path.isfile(plant_info_script) and plant_fa.is_file() and consensus_tax.is_file():
                cobra_dir = self.d['cobra']
                rescue_fa = self.d['rescue_dir'] / "HQ_plant_viruses.fasta"
                clusters_tsv = self.d['cluster'] / "3_vclust" / "vclust_clusters.tsv"
                cmd = (f"python {plant_info_script} --output-dir {plant_dir} "
                       f"--taxonomy {consensus_tax} --plant-fasta {plant_fa} --cobra-dir {cobra_dir}")
                if clusters_tsv.is_file():
                    cmd += f" --clusters {clusters_tsv}"
                if rescue_fa.is_file():
                    cmd += f" --rescue-fa {rescue_fa}"
                run_cmd(cmd, self.log, "Plant Virus Info", str(logs_dir / "plant_virus_info.log"))
            else:
                self.log.info("  All_plant.viruses_info.tsv — 跳过 (缺脚本/植物fasta/taxonomy)")
        else:
            self.log.info("  All_plant.viruses_info.tsv — 已存在, 跳过")

        # ── Part 3: Viroid batch analysis ──
        viroid_batch = str(SCRIPT_DIR / "utils" / "analyze_viroid_batch.py")
        if os.path.isfile(viroid_batch):
            ident_dir = self.d["ident"]
            vdb = getattr(self.args, "viroids_db", None)
            vtax = Path(str(vdb)).parent / "viroids.taxonomy_info.tsv" if vdb else None
            vfasta = Path(str(vdb)).parent / "viroids.fasta" if vdb else None
            if vdb and vtax and vtax.is_file() and vfasta and vfasta.is_file():
                viroid_out = viroid_dir
                cmd = f"python {viroid_batch} --ident-dir {ident_dir} --taxonomy {vtax} --ref-fasta {vfasta} -o {viroid_out}"
                run_cmd(cmd, self.log, "Viroid Batch Analysis",
                        str(logs_dir / "viroid_analysis.log"))
            else:
                self.log.info("  类病毒分析 — 跳过 (缺 viroids_db)")
        else:
            self.log.info("  类病毒分析 — 跳过 (analyze_viroid_batch.py 不存在)")

        # ── Part 4: Rescue 候选病毒检测 (salmon 定量, 得到检出频次) ──
        # 对全部 HQ 候选 (dest_fasta) 做 salmon 检测, 得到检出频次 (n_samples)
        # 作为证据整合的输入 (跨样本支持度)
        batch_depth_script = SCRIPT_DIR.parent / "virome_analysis_pipeline" / "batch_virus_depth.py"
        frequency_tsv = None
        if batch_depth_script.is_file() and dest_fasta.is_file():
            detect_out = self.d['analysis'] / "rescue_detection"
            detect_out.mkdir(parents=True, exist_ok=True)
            detect_summary = detect_out / "summary" / "all_viruses.raw.tsv"
            if not detect_summary.is_file():
                ref_info = detect_out / "rescue_ref_info.tsv"
                if not ref_info.is_file():
                    self._write_detect_ref_info(dest_fasta, ref_info)
                def _has_seq(d):
                    # 目录存在且含序列文件才可用 (00a_CleanData 可能只剩日志)
                    if not d.is_dir():
                        return False
                    for f in d.iterdir():
                        if f.is_file() and any(
                                f.name.lower().endswith(e) for e in SEQ_EXTS):
                            return True
                    return False
                reads_dir = (self.d['clean'] if _has_seq(self.d['clean'])
                             else (self.d['hostdep'] if _has_seq(self.d['hostdep'])
                                   else self.d['raw']))
                dcmd = (f"{sys.executable} {batch_depth_script} "
                        f"--input_dir {reads_dir} "
                        f"--output_dir {detect_out} "
                        f"--ref_info {ref_info} "
                        f"--reference {dest_fasta} "
                        f"--tool salmon --threads {t} --align_threads 8 --batch_size 20 "
                        f"--coverage 10.0 --ratio 0.3 --meandepth 0.5 --min_tpm 1.0 --min_uniq_reads 10 "
                        f"--sp_thresh 95.0 --resume")
                run_cmd(dcmd, self.log, "Rescue Detection (salmon)",
                        str(logs_dir / "rescue_detection.log"))
            else:
                self.log.info("  rescue 检测 — 已存在, 跳过")

            # 生成频次表 (contig_id \t n_samples), 供 Part 5 证据整合用
            if detect_summary.is_file():
                frequency_tsv = detect_out / "frequency_table.tsv"
                self._write_frequency_table(detect_summary, frequency_tsv, dest_fasta,
                                            self.d['root'] / "prevalence_full_table.tsv")
        else:
            self.log.info("  rescue 检测 — 跳过 (batch_virus_depth.py 或 HQ fasta 不存在)")

        # ── Part 5 已移至 09b (virus_validation) ──
        # CDD 验证 + 五层证据整合 + verdict 判定 → 09b_ACVirus_Analysis/virus_validation/
        self.log.info("  [9a] Part 5 (CDD/blast 验证/verdict) 已转移 → 09b virus_validation")

        self.log.info("  分析完成 → %s", self.d['analysis'])

    def _write_detect_ref_info(self, fasta_path, out_path):
        """生成 rescue 检测 ref_info (Accession=fasta 头, Species 空, 由 batch_virus_depth 用 fasta 头注释)"""
        from Bio import SeqIO
        n = 0
        with open(out_path, 'w') as f:
            f.write('Accession\tTaxid\tSpecies_NCBI\n')
            for rec in SeqIO.parse(str(fasta_path), 'fasta'):
                f.write(f'{rec.id}\t\t\n')
                n += 1
        self.log.info('  rescue ref_info (%d 条)', n)

    def _write_frequency_table(self, summary_path, out_path, fasta_path=None, prevalence_path=None):
        """从检测 summary 统计检出样本数, 生成频次表 (contig_id + n_samples)
        合并 prevalence_full_table 的 checkv 质量字段, 覆盖全部候选 (含 0 检出)
        """
        from collections import defaultdict
        from Bio import SeqIO
        freq = defaultdict(set)
        with open(summary_path, errors='ignore') as f:
            header = f.readline().rstrip('\n').split('\t')
            col = {h: i for i, h in enumerate(header)}
            acc_col = col.get('Rep_Accession', col.get('Accession', None))
            sample_col = col.get('Sample', 0)
            if acc_col is None:
                self.log.warning("  频次表 — 找不到 Rep_Accession/Accession 列")
                return
            for line in f:
                p = line.split('\t')
                if len(p) <= max(acc_col, sample_col):
                    continue
                acc = p[acc_col]
                sample = p[sample_col]
                if not acc or not sample.startswith(('SRR', 'CRR', 'ERR')):
                    continue
                freq[acc].add(sample)
        # 总样本数 (检测分母, 用于计算 prevalence_%)
        total_samples = 0
        for _acc, _samps in freq.items():
            total_samples = max(total_samples, len(_samps))
        # 更准确: 从 summary 的 Sample 列去重 (含无检出的样本)
        try:
            with open(summary_path, errors='ignore') as _sf:
                _hdr = _sf.readline().rstrip('\n').split('\t')
                _sc = {h: i for i, h in enumerate(_hdr)}.get('Sample', 0)
                _sf.seek(0); _sf.readline()
                _all_samps = set()
                for _line in _sf:
                    _pp = _line.split('\t')
                    if len(_pp) > _sc and _pp[_sc].startswith(('SRR', 'CRR', 'ERR')):
                        _all_samps.add(_pp[_sc])
                if _all_samps:
                    total_samples = len(_all_samps)
        except Exception:
            pass
        # 读 prevalence_full_table 的 checkv 质量字段 (合并进频次表, 避免质量层丢失)
        cv_info = {}
        if prevalence_path and os.path.exists(prevalence_path):
            import csv as _csv
            for r in _csv.DictReader(open(prevalence_path), delimiter='\t'):
                cid = r.get('contig_id', '')
                if cid:
                    cv_info[cid] = r
        # 全部候选 (含 0 检出)
        all_ids = set()
        if fasta_path and os.path.exists(fasta_path):
            for rec in SeqIO.parse(str(fasta_path), 'fasta'):
                all_ids.add(rec.id)
        else:
            all_ids = set(freq.keys())
        with open(out_path, 'w') as f:
            f.write('contig_id\tn_samples\tcheckv_completeness\tcheckv_confidence\tgenus_avg_len\tavg_ANI\tprevalence_%\n')
            for cid in sorted(all_ids):
                cv = cv_info.get(cid, {})
                n_s = len(freq.get(cid, set()))
                # prevalence_% 用检测频次重新计算 (与 n_samples 配套)
                pct = f'{n_s / total_samples * 100:.2f}' if total_samples > 0 else ''
                f.write(f'{cid}\t{n_s}\t'
                        f'{cv.get("checkv_completeness","")}\t{cv.get("checkv_confidence","")}\t'
                        f'{cv.get("genus_avg_len","")}\t{cv.get("avg_ANI","")}\t{pct}\n')
        detected = sum(1 for cid in all_ids if len(freq.get(cid, set())) > 0)
        self.log.info('  频次表: %d/%d 候选检出 (→ %s)', detected, len(all_ids), out_path.name)

    def _write_ref_info(self, fasta_path, out_path, consensus_path):
        """生成 HQ_plant_viruses_info.tsv: 整合 05_Taxonomy 分类"""
        import pandas as pd
        from Bio import SeqIO

        # ── 05_Taxonomy ──
        consensus = {}
        if consensus_path.is_file():
            try:
                cdf = pd.read_csv(consensus_path, sep='\t', quotechar='"')
                for _, row in cdf.iterrows():
                    cid = str(row.get('contig_id', row.iloc[0]))
                    consensus[cid] = {k: str(row.get(k, '')) for k in
                        ['Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species','confidence','primary_tool','Nucleic_acid']}
            except Exception: pass

        # ── 序列长度 ──
        seq_lens = {}
        for rec in SeqIO.parse(str(fasta_path), 'fasta'):
            seq_lens[rec.id] = len(rec.seq)

        with open(out_path, 'w') as f:
            f.write('Accession\tLength\tSpecies\tGenus\tFamily\tRealm\t'
                    'Kingdom\tClass\tOrder\t'
                    'CDS_Count\tPrimary_Tool\tConfidence\tNucleic_acid\n')
            for cid, slen in seq_lens.items():
                cs = consensus.get(cid, {})
                f.write(f'{cid}\t{slen}\t'
                        f'{cs.get("Species","")}\t{cs.get("Genus","")}\t{cs.get("Family","")}\t{cs.get("Realm","")}\t'
                        f'{cs.get("Kingdom","")}\t{cs.get("Class","")}\t{cs.get("Order","")}\t'
                        f'\t{cs.get("primary_tool","")}\t{cs.get("confidence","")}\t'
                        f'{cs.get("Nucleic_acid","")}\n')
        self.log.info('  HQ_plant_viruses_info.tsv (%d 条)', len(seq_lens))

    def run_analysis_verify(self):
        """[9b] 分析验证 — KEEP 序列的 ACVirus 重新分类 + 新病毒鉴定 + 按科进化树

        输出到 09b_Analysis_Verify/:
          class_KEEP.fasta            - 从 09a 提取的 KEEP 序列
          acvirus_classify/           - ACVirus classify 结果 (科/属级分类 + 置信度)
          acvirus_identify/           - ACVirus identify 结果 (新病毒 Novelty 打分)
          acvirus_trees/<Family>/     - 各科进化树 (Tree-Pro, 属级分色 + 共线性)
        """
        self.d['analysis_verify'].mkdir(parents=True, exist_ok=True)
        self.log.info("=" * 50)
        self.log.info("[9b] 分析验证 — 病毒验证/新病毒鉴定/进化树")
        import subprocess as _sp
        import shutil as _sh

        acv_script = SCRIPT_DIR / "utils" / "acvirus_tree_pro.py"
        if not acv_script.is_file():
            acv_script = Path.home() / "bin" / "acvirus_tree_pro.py"  # fallback
        acv_dir = self.d['analysis_verify']
        src_analysis = self.d['analysis']

        # ── Part V: virus_validation (CDD 验证 + blast 身份 + 五层证据整合, 由 9a 迁入) ──
        dest_fasta = src_analysis / "HQ_plant_viruses.fasta"
        if not dest_fasta.is_file():
            self.log.warning("[9b] 缺 HQ_plant_viruses.fasta, 跳过")
            return
        validate_out = acv_dir / "virus_validation"
        validate_out.mkdir(parents=True, exist_ok=True)
        (self.d['analysis'] / 'logs').mkdir(parents=True, exist_ok=True)  # mix 等旧结构可能缺 logs
        t = getattr(self.args, 'threads', 0) or 60
        sample_name = getattr(self.args, 'tax_sample_name', None) or "Votus"
        consensus_tax = self.d['taxonomy'] / f"{sample_name}.integrated" / "final_integrated_classification.tsv"
        frequency_tsv = self.d['analysis'] / "rescue_detection" / "frequency_table.tsv"
        validate_script = SCRIPT_DIR / "validate_rescue_cdd.py"
        if validate_script.is_file():
            cdd_report = validate_out / "cdd_evidence_report.tsv"
            if not cdd_report.is_file() or cdd_report.stat().st_size < 100:
                vcmd = (f"{sys.executable} {validate_script} "
                        f"--input {dest_fasta} "
                        f"--taxonomy {consensus_tax} "
                        f"--outdir {validate_out} "
                        f"--threads {t} --mode filter")
                run_cmd(vcmd, self.log, "Virus CDD Validation",
                        str(self.d['analysis'] / 'logs' / "virus_cdd_validation.log"))
            else:
                self.log.info("  [9b] virus CDD 验证 — 已存在, 跳过")

            gen_script = SCRIPT_DIR / "gen_final_judgement.py"
            cdd_hits_tsv = validate_out / "cdd_hits.tsv"
            freq_gen = frequency_tsv if (frequency_tsv and frequency_tsv.is_file()) else (self.d['root'] / "prevalence_full_table.tsv")
            if gen_script.is_file() and not (self.d['root'] / "final_judgement_table.tsv").is_file():
                import re as _re
                root_name = Path(self.d['root']).name
                m = _re.search(r'(Lycium_[a-z]+|Alternaria_[a-z]+|Aphis_[a-z]+|Fusarium_[a-z]+|Neoceratitis_[a-z]+)', root_name)
                sample_name2 = m.group(1).lower() if m else root_name
                gcmd = (f"{sys.executable} {gen_script} "
                        f"--fasta {dest_fasta} "
                        f"--freq {freq_gen} "
                        f"--cdd {cdd_hits_tsv} "
                        f"--cdd-evidence {cdd_report} "
                        f"--sample {sample_name2} "
                        f"--threads {t} "
                        f"--out {self.d['root'] / 'final_judgement_table.tsv'}")
                run_cmd(gcmd, self.log, "Generate final_judgement_table",
                        str(self.d['analysis'] / 'logs' / "virus_gen_judgement.log"))
            elif gen_script.is_file():
                self.log.info("  final_judgement_table.tsv — 已存在, 跳过")

            # ── 第三路探针: CT3 病毒 HMM 库 (pyhmmer 六框翻译 + 六桶分级) ──
            # 一道门 (CDD) 加两次补票 (blast 身份 / HMM 全长 profile)。
            # CT3 库是二进制 .h3m, 五库分别扫描后按 contig 取并集; 缺脚本或缺库时优雅跳过。
            hmm_script = SCRIPT_DIR / "hmm_ct3_evidence.py"
            hmm_dir = Path(getattr(self.args, 'ct3_hmm_dir', None)
                           or (Path.home() / 'database' / 'virus-db' / 'ct3_DBs' / 'hmmscan_DBs' / 'v3.1.1'))
            hmm_summary = validate_out / "hmm_ct3_hits.tsv"
            hmm_report_arg = None
            if getattr(self.args, 'skip_ct3_hmm', False):
                self.log.info("  [9b] CT3-HMM 第三路探针 — 已按 --skip-ct3-hmm 跳过")
            elif not hmm_script.is_file():
                self.log.info("  [9b] CT3-HMM 第三路探针 — 跳过 (缺 hmm_ct3_evidence.py)")
            elif not hmm_dir.is_dir():
                self.log.info("  [9b] CT3-HMM 第三路探针 — 跳过 (缺 HMM 库目录 %s)", hmm_dir)
            else:
                if getattr(self.args, 'force', False) or not hmm_summary.is_file() or hmm_summary.stat().st_size < 100:
                    hcmd = (f"{sys.executable} {hmm_script} "
                            f"--fasta {dest_fasta} "
                            f"--ct3-dir {hmm_dir} "
                            f"--outdir {validate_out} "
                            f"--threads {t}")
                    run_cmd(hcmd, self.log, "CT3-HMM Third-Probe Validation",
                            str(self.d['analysis'] / 'logs' / "virus_ct3_hmm.log"))
                else:
                    self.log.info("  [9b] CT3-HMM 第三路探针 — 已存在, 跳过")
                if hmm_summary.is_file():
                    hmm_report_arg = hmm_summary

            integrate_script = SCRIPT_DIR / "integrate_rescue_evidence.py"
            judgement_tsv = self.d['root'] / "final_judgement_table.tsv"
            prevalence_tsv = self.d['root'] / "prevalence_full_table.tsv"
            freq_input = frequency_tsv if (frequency_tsv and frequency_tsv.is_file()) else prevalence_tsv
            if integrate_script.is_file() and judgement_tsv.is_file():
                scored = validate_out / "rescue_evidence_scored.tsv"
                if getattr(self.args, 'force', False) or not scored.is_file() or scored.stat().st_size < 100:
                    icmd = (f"{sys.executable} {integrate_script} "
                            f"--judgement {judgement_tsv} "
                            f"--cdd-report {cdd_report} "
                            f"--prevalence {freq_input} "
                            f"--taxonomy {consensus_tax} "
                            f"--fasta {dest_fasta} "
                            f"--outdir {validate_out}")
                    if hmm_report_arg:
                        icmd += f" --hmm-report {hmm_report_arg}"
                    # rescue 默认开 (integrate_rescue_evidence.py 侧 default True)
                    if getattr(self.args, 'hmm_min_specific', 1) != 1:
                        icmd += f" --hmm-min-specific {self.args.hmm_min_specific}"
                    if getattr(self.args, 'hmm_require_family', False):
                        icmd += " --hmm-require-family"
                    if getattr(self.args, 'hmm_count_ambiguous', False):
                        icmd += " --hmm-count-ambiguous"
                    if getattr(self.args, 'no_hmm_rescue', False):
                        icmd += " --no-hmm-rescue"
                    if getattr(self.args, 'no_hmm_in_score', False):
                        icmd += " --no-hmm-in-score"
                    # 升 KEEP 默认开 (integrate_rescue_evidence.py 侧 default True)
                    if getattr(self.args, 'no_hmm_promote_keep', False):
                        icmd += " --no-hmm-promote-keep"
                    if getattr(self.args, 'hmm_promote_min', 1) != 1:
                        icmd += f" --hmm-promote-min {self.args.hmm_promote_min}"
                    run_cmd(icmd, self.log, "Virus Evidence Integration",
                            str(self.d['analysis'] / 'logs' / "virus_evidence_integration.log"))
                else:
                    self.log.info("  五层证据整合 — 已存在, 跳过")
            else:
                self.log.info("  五层证据整合 — 跳过 (缺 integrate_rescue_evidence.py 或 final_judgement_table.tsv)")
        else:
            self.log.info("  virus CDD 验证 — 跳过 (validate_rescue_cdd.py 不存在)")

        # 读 09b virus_validation 的 KEEP (兼容旧 09a/rescue_validation)
        scored = validate_out / "rescue_evidence_scored.tsv"
        if not scored.is_file():
            scored = src_analysis / "rescue_validation" / "rescue_evidence_scored.tsv"
        hq = dest_fasta
        if not scored.is_file() or not hq.is_file():
            self.log.warning("[9b] 缺 scored.tsv 或 HQ fasta, 跳过")
            return
        try:
            import csv as _csv
            from Bio import SeqIO as _sio
            MIN_KEEP_LEN = getattr(self.args, 'keep_min_len', 0)  # KEEP 最小长度 bp (0=不过滤)
            keep_ids = set()
            hmm_keep_ids = set()  # CT3-HMM 比对命中提升的 KEEP: 不受长度门限制
            with open(scored) as fl:
                for row in _csv.DictReader(fl, delimiter='\t'):
                    if row.get('verdict') == 'KEEP':
                        keep_ids.add(row.get('contig_id', ''))
                        if (row.get('hmm_rescue', '') or '').strip().endswith('->KEEP'):
                            hmm_keep_ids.add(row.get('contig_id', ''))
            keep_fa = acv_dir / "class_KEEP.fasta"
            n_dropped_short = 0
            n_short_kept = 0
            with open(keep_fa, 'w') as ko, open(hq) as hqi:
                for rec in _sio.parse(hqi, 'fasta'):
                    if rec.id in keep_ids:
                        if len(rec.seq) >= MIN_KEEP_LEN or rec.id in hmm_keep_ids:
                            if len(rec.seq) < MIN_KEEP_LEN:
                                n_short_kept += 1
                            _sio.write(rec, ko, 'fasta')
                        else:
                            n_dropped_short += 1
            n = sum(1 for _ in _sio.parse(str(keep_fa), 'fasta'))
            self.log.info('  [9b] 提取 KEEP %d 条 → class_KEEP.fasta (筛掉 <%dbp %d 条; HMM 比对豁免长度门 %d 条)',
                          n, MIN_KEEP_LEN, n_dropped_short, n_short_kept)
        except Exception as e:
            self.log.warning('[9b] 提取 KEEP 失败: %s', e)
            return

        if not acv_script.is_file():
            self.log.warning('[9b] 缺 acvirus_tree_pro.py (%s), 跳过构建', acv_script)
            return
        db = Path.home() / 'database' / 'virus-db' / 'acvirus_db'
        db_taxa = db / 'taxa.txt'
        db_fasta = db / 'all_virus.fasta'
        if not db_taxa.is_file() or not db_fasta.is_file():
            self.log.warning('[9b] 缺 ACVirus 参考库, 跳过')
            return

        # classify (收集科/属 - 用 classify 结果确定按哪些科建树)
        cls_out = acv_dir / 'acvirus_classify'
        cls_out.mkdir(exist_ok=True)
        try:
            r = _sp.run(['python3', str(Path.home()/'MMPV-RNA'/'biosoft'/'ACVirus'/'cli.py'),
                         'classify', '--contigs', str(keep_fa), '--data_path', str(db),
                         '--out', str(cls_out), '-t', str((getattr(getattr(self,'args',None),'threads',0) or 60))],
                        capture_output=True, text=True, timeout=1800)
            self.log.info('  [9b] classify rc=%d', r.returncode)
        except Exception as e:
            self.log.warning('[9b] classify 失败: %s', e)

        # identify (新病毒鉴定)
        try:
            idt_out = acv_dir / 'acvirus_identify'
            idt_out.mkdir(exist_ok=True)
            r = _sp.run(['python3', str(Path.home()/'MMPV-RNA'/'biosoft'/'ACVirus'/'cli.py'),
                         'identify', '--classify_out', str(cls_out), '--data_path', str(db),
                         '--out', str(idt_out/'identification_result_score.csv'),
                         '--min_cluster', '3', '-t', str((getattr(getattr(self,'args',None),'threads',0) or 60))],
                        capture_output=True, text=True, timeout=1800)
            self.log.info('  [9b] identify rc=%d', r.returncode)
        except Exception as e:
            self.log.warning('[9b] identify 失败: %s', e)

        # 按科建树 (Tree-Pro macro 模式, 属级分色) - 并发执行
        try:
            from concurrent.futures import ThreadPoolExecutor, as_completed
            tree_root = acv_dir / 'acvirus_trees'
            tree_root.mkdir(exist_ok=True)
            # 从 classify 结果确定主要科, 并按科拆分 KEEP (每棵树只用科内序列, 避免全量陪跑)
            cls_res = cls_out / 'final_result_with_confidence.tsv'
            if cls_res.is_file():
                import csv as _csv2
                fam_of = {}
                fams = set()
                with open(cls_res) as cr:
                    for row in _csv2.DictReader(cr, delimiter='\t'):
                        f = (row.get('Family', '') or '').strip()
                        cid = (row.get('Nucleotide', '') or row.get('Contig', '') or '').strip()
                        if f and f != 'NA':
                            fams.add(f)
                            if cid: fam_of[cid] = f
                # ── 建树过滤: len_ratio>=0.7 (接近属平均长度) 或属长缺失才保留 ──
                # genus_lens 库 (g__Genus -> 平均长度); 优先 v2 表 (含节段属 ref_min)
                _gl = {}
                _gmin = {}   # 节段属: 最短完整单元 (单段 contig 的 ratio 基准)
                for _glfile in [Path.home() / 'database' / 'virus-db' / 'db' / 'genus_lens_v2.tsv',
                                Path.home() / 'database' / 'virus-db' / 'db' / 'genus_lens',
                                SCRIPT_DIR.parent / 'database' / 'genus_lens_v2.tsv',
                                SCRIPT_DIR.parent / 'database' / 'genus_lens']:
                    if _glfile.is_file():
                        for _line in open(_glfile):
                            _p = _line.rstrip().split('\t')
                            if len(_p) >= 2 and _p[0].startswith('g__'):
                                try: _gl[_p[0][3:]] = float(_p[1])
                                except Exception: pass
                                if len(_p) >= 5 and _p[2].strip().lower() in ('yes', 'true', '1'):
                                    try: _gmin[_p[0][3:]] = float(_p[3])
                                    except Exception: pass
                        break
                # classify Genus -> 属平均长度
                _genus_of = {}
                with open(cls_res) as _cr2:
                    for _row in _csv2.DictReader(_cr2, delimiter='\t'):
                        _cid = (_row.get('Nucleotide','') or '').strip()
                        _g = (_row.get('Genus','') or '').strip()
                        if _cid and _g: _genus_of[_cid] = _g
                # 拆分 class_KEEP.fasta → acvirus_trees/<Fam>/keep_<Fam>.fasta (仅 len_ratio>=0.7 或属长缺失)
                MIN_RATIO = getattr(self.args, 'tree_min_ratio', 0.7)
                fam_fa = {}
                n_dropped_ratio = 0
                for rec in _sio.parse(str(keep_fa), 'fasta'):
                    f = fam_of.get(rec.id) or fam_of.get(rec.id.split()[0])
                    if not f: continue
                    # len_ratio 判定 (节段属以最短完整单元为基准, 单段 contig 不被误滤)
                    _g = _genus_of.get(rec.id, '')
                    _glen = _gl.get(_g, 0.0) or 0.0
                    if _g in _gmin:
                        _glen = _gmin[_g]
                    if _glen > 0:
                        _ratio = len(rec.seq) / _glen
                        if _ratio < MIN_RATIO:
                            n_dropped_ratio += 1
                            continue
                    fam_fa.setdefault(f, []).append(rec)
                for f, recs in fam_fa.items():
                    d = tree_root / f; d.mkdir(parents=True, exist_ok=True)
                    with open(d / f'keep_{f}.fasta', 'w') as fo:
                        for rc_ in recs: _sio.write(rc_, fo, 'fasta')
                self.log.info('  [9b] 按科拆分 KEEP: %s (筛掉 len_ratio<%.1f %d 条)', ', '.join(f'{k}×{len(v)}' for k,v in sorted(fam_fa.items(), key=lambda x:-len(x[1]))[:8]), MIN_RATIO, n_dropped_ratio)
                total_th = getattr(getattr(self,'args',None),'threads',0) or 60
                fam_threads = max(4, total_th // 4)
                n_par = max(1, min(3, len(fams)))
                self.log.info('  [9b] 并行建树 %d 科: %d 路 x %d threads/科', len(fams), n_par, fam_threads)
                def _build_one(fam):
                    fam_out = tree_root / fam
                    try:
                        # P3: 已有树则跳过 (重启时复用, 免数小时重跑)
                        if list(fam_out.glob('*.treefile')):
                            return fam, 4, 'tree exists'
                        # P1: 无本科 KEEP 子集则跳过 (无自研 contig; 全量兼底会让 MAFFT 超时)
                        fam_keep = fam_out / f'keep_{fam}.fasta'
                        if not fam_keep.exists():
                            return fam, 3, 'no family-specific KEEP contig'
                        use_fa = fam_keep
                        # P4: 输入过长则跳过 (巨基因组参考会让 MAFFT 超 4h 超时;
                        #     植物病毒基因组 <=30kb, 单条 >200kb 或总量 >5Mb 视为病态输入)
                        try:
                            _tot = 0; _mx = 0
                            for _rc0 in _sio.parse(str(use_fa), 'fasta'):
                                _l = len(_rc0.seq); _tot += _l
                                if _l > _mx: _mx = _l
                            if _mx > 200000 or _tot > 5000000:
                                return fam, 5, 'alignment too large (max %.1fkb, total %.1fMb)' % (_mx/1e3, _tot/1e6)
                        except Exception:
                            pass
                        r = _sp.run(['python3', str(acv_script),
                                     '--mode', 'macro', '--target_name', fam, '--target_rank', 'Family',
                                     '--contigs', str(use_fa), '--db_taxa', str(db_taxa),
                                     '--db_fasta', str(db_fasta), '--outdir', str(fam_out),
                                     '--threads', str(fam_threads)],
                                    capture_output=True, text=True, timeout=14400)
                        return fam, r.returncode, (r.stderr or '')[-200:]
                    except Exception as e:
                        return fam, -1, str(e)[:200]
                with ThreadPoolExecutor(max_workers=n_par) as ex:
                    futs = [ex.submit(_build_one, fam) for fam in sorted(fams)]
                    for fu in as_completed(futs):
                        fam, rc, err = fu.result()
                        if rc == 0:
                            self.log.info('  [9b] tree %s ✓', fam)
                        elif rc in (3, 4, 5):
                            self.log.info('  [9b] tree %s ⊘ %s', fam, err)
                        else:
                            self.log.warning('  [9b] tree %s ✗ rc=%d %s', fam, rc, err)
        except Exception as e:
            self.log.warning('[9b] 按科建树失败: %s', e)

        # ── SDT 属级序列一致性矩阵 (接近属平均长度的 KEEP) ──
        if not getattr(self.args, 'skip_sdt', False):
            try:
                self._run_sdt_matrices(acv_dir)
            except Exception as e:
                self.log.warning('[9b] SDT 矩阵失败: %s', e)

        try:
            self._write_keep_len_filter_table(acv_dir)
        except Exception as e:
            self.log.warning('[9b] KEEP 过滤表生成失败: %s', e)

        self.log.info('  09b 分析验证 完成 → %s', acv_dir)

    def _write_keep_len_filter_table(self, acv_dir):
        """[9b] 生成 KEEP 长度/len_ratio 过滤判定表 → 09b/KEEP_len_filter_table.tsv
        反映两层: 提取时 <keep-min-len 真过滤 + 建树时 len_ratio>=tree_min_ratio 精简版。
        属长缺失不过滤。仅作记录, 不改动 class_KEEP.fasta。"""
        import csv as _c
        from collections import Counter as _ctr
        ksum = acv_dir / 'keep_summary.tsv'
        if not ksum.is_file():
            self.log.info('  [9b] KEEP 过滤表跳过 (无 keep_summary.tsv)')
            return
        _gl = {}
        _gmin = {}   # 节段属: 最短完整单元 (ratio 基准)
        for _pgl in [Path.home() / 'database' / 'virus-db' / 'db' / 'genus_lens_v2.tsv',
                     Path.home() / 'database' / 'virus-db' / 'db' / 'genus_lens',
                     SCRIPT_DIR.parent / 'database' / 'genus_lens_v2.tsv',
                     SCRIPT_DIR.parent / 'database' / 'genus_lens']:
            if _pgl.is_file():
                for _line in open(_pgl):
                    _pp = _line.rstrip().split('\t')
                    if len(_pp) >= 2 and _pp[0].startswith('g__'):
                        try:
                            _gl[_pp[0][3:]] = float(_pp[1])
                        except Exception:
                            pass
                        if len(_pp) >= 5 and _pp[2].strip().lower() in ('yes', 'true', '1'):
                            try:
                                _gmin[_pp[0][3:]] = float(_pp[3])
                            except Exception:
                                pass
                break
        MIN_RATIO = getattr(self.args, 'tree_min_ratio', 0.7)
        rows = list(_c.DictReader(open(ksum), delimiter='\t'))
        for r in rows:
            try:
                L = float(r.get('length', '') or 0)
            except Exception:
                L = 0.0
            genus = (r.get('Genus', '') or '').strip()
            gal = r.get('genus_avg_len', '').strip()
            avg = None
            if gal:
                try:
                    avg = float(gal)
                except Exception:
                    pass
            if avg is None:
                avg = _gl.get(genus)
                if genus in _gmin:      # 节段属: ratio 以最短完整单元为基准
                    avg = _gmin[genus]
            ratio = None
            if avg and avg > 0 and L > 0:
                ratio = L / avg
            if avg is None or avg <= 0 or ratio is None:
                status = 'KEEP_no_genus_len'
            elif ratio < MIN_RATIO:
                status = 'FILTERED'
            else:
                status = 'KEEP'
            r['len_ratio_recalc'] = ('%.3f' % ratio) if ratio is not None else ''
            r['filter_status'] = status
        newcols = list(rows[0].keys()) + ['len_ratio_recalc', 'filter_status']
        _order = {'FILTERED': 0, 'KEEP_no_genus_len': 1, 'KEEP': 2}
        rows.sort(key=lambda r: _order.get(r['filter_status'], 3))
        dst = acv_dir / 'KEEP_len_filter_table.tsv'
        with open(dst, 'w', newline='', encoding='utf-8') as fo:
            w = _c.DictWriter(fo, fieldnames=newcols, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            for r in rows:
                w.writerow(r)
        cc = _ctr(r['filter_status'] for r in rows)
        self.log.info('  [9b] KEEP 过滤表 %s: 总%d | <%.1f精简剔除%d | 无属长保留%d | 保留%d',
                      dst.name, len(rows), MIN_RATIO,
                      cc.get('FILTERED', 0), cc.get('KEEP_no_genus_len', 0), cc.get('KEEP', 0))

    def _run_sdt_matrices(self, acv_dir):
        """[9b] SDT 属级序列一致性矩阵 — 对接近属平均长度的 KEEP (len_ratio>=sdt_min_ratio) 按属绘制
        调 virome_analysis_pipeline/sdt_genus_matrix.py, 输出 09b/SDT_matrix/<属>/"""
        import csv as _c
        import subprocess as _sp  # 独立方法作用域, 必须自行导入 (曾缺此报 NameError)
        sdt_script = SCRIPT_DIR.parent / 'virome_analysis_pipeline' / 'sdt_genus_matrix.py'
        cls_res = acv_dir / 'acvirus_classify' / 'final_result_with_confidence.tsv'
        if not sdt_script.is_file() or not cls_res.is_file():
            self.log.info('  [9b] SDT 跳过 (缺 sdt_genus_matrix.py 或 classify 结果)')
            return
        min_ratio = getattr(self.args, 'sdt_min_ratio', 0.0)
        sdt_out = acv_dir / 'SDT_matrix'
        sdt_out.mkdir(parents=True, exist_ok=True)
        # 属 → 用户 contig 列表 (从 classify)
        from collections import defaultdict as _dd
        gen_of = _dd(list)
        keep_ids = set()
        with open(acv_dir / 'class_KEEP.fasta') as _f:
            for _line in _f:
                if _line.startswith('>'): keep_ids.add(_line[1:].strip())
        # 长度/len_ratio: 从 keep_summary 取 (若存在)
        ratio_map = {}
        ksum = acv_dir / 'keep_summary.tsv'
        if ksum.is_file():
            for r in _c.DictReader(open(ksum), delimiter='\t'):
                try: ratio_map[r['contig_id']] = float(r['len_ratio'] or 0)
                except Exception: pass
        with open(cls_res) as _f:
            for row in _c.DictReader(_f, delimiter='\t'):
                cid = (row.get('Nucleotide','') or '').strip()
                g = (row.get('Genus','') or '').strip()
                if cid in keep_ids and g and g != 'NA':
                    gen_of[g].append(cid)
        sdt_th = getattr(self.args, 'sdt_threads', 8)
        # SDT 是 O(N^2) 全对全精确比对 (每对一次外部 MAFFT), 属内序列过多时
        # 组合数爆炸 (519 条 → 13.4 万对), 必须设上限。超限属跳过并记日志。
        sdt_max = getattr(self.args, 'sdt_max_seqs', 80)
        n_run = 0
        n_skip_big = 0
        for g, cids in sorted(gen_of.items(), key=lambda x: -len(x[1])):
            user_ok = [c for c in cids if ratio_map.get(c, 1.0) >= min_ratio]
            if len(user_ok) < 2:
                continue
            if len(user_ok) > sdt_max:
                n_skip_big += 1
                self.log.info('  [9b] SDT %s 跳过: %d 条 > 上限 %d (O(N^2) 组合 %d 对)',
                              g, len(user_ok), sdt_max, len(user_ok)*(len(user_ok)-1)//2)
                continue
            od = sdt_out / g
            done = any((od / f'SDT_{g}.png').exists() for _ in [0]) or (od / f'SDT_{g}.pdf').exists()
            if done:
                self.log.info('  [9b] SDT %s 已存在, 跳过', g)
                continue
            od.mkdir(parents=True, exist_ok=True)
            self.log.info('  [9b] SDT %s: %d 条用户 contig (ratio>=%.1f)', g, len(user_ok), min_ratio)
            try:
                r = _sp.run(['python3', str(sdt_script),
                             '--analysis_dir', str(acv_dir), '--genus', g, '-o', str(od),
                             '--threads', str(sdt_th), '--short_labels', '--palette', 'cividis',
                             '--resume', '--other_count', '0', '--target_cap', '20'], capture_output=True, text=True, timeout=14400)
                self.log.info('  [9b] SDT %s rc=%d', g, r.returncode)
                if r.returncode != 0 and r.stderr:
                    self.log.warning('  [9b] SDT %s stderr: %s', g, r.stderr[-200:])
                n_run += 1
            except Exception as e:
                self.log.warning('  [9b] SDT %s 异常: %s', g, e)
        self.log.info('  [9b] SDT 矩阵: %d 个属 → %s (超上限跳过 %d 个)',
                      n_run, sdt_out, n_skip_big)


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
            import time as _t
            _t0 = _t.time()
            # 报告需遍历全流程产物 (含数千目录), 600s 不够 → 3600s
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            self.log.info("  report_pipeline 耗时 %.1fs (rc=%d)", _t.time() - _t0, result.returncode)
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
    'clean':         '00a: Fastp质控 + Seqkit转FASTA + Clumpify聚类重排',
    'deplete':       '00b: Kraken2 + Bowtie2/HISAT2/Minimap2 + rRNA去除',
    'assembly':      '01:  MEGAHIT / rnaviralSPAdes / Penguin 三工具组装',
    'identification':'02a: 10工具并行病毒鉴定 (Genomad/Diamond/VirSorter2/...)',
    'filter':        '02b: UniProt + CDD 分层过滤 (保留viral/viral_plant/viral_bacteria)',
    'cobra':         '03a: BWA-MEM2 + COBRA-Meta 重叠群延伸',
    'merge':         '03b: 合并多样本COBRA结果 + 可选Flye共组装',
    'cluster':       '04:  CD-HIT参考预聚类 + vclust Leiden聚类',
    'taxonomy':      '05:  8工具分类 + R加权投票共识',
    'host':          '06:  ICTV > RNAVirHost > PhaBOX2 宿主预测',
    'checkv':        '07:  CheckV 完整性评估',
    'rescue':        '08:  三支路级联拯救 (CheckV -> Virseqimprover -> BLASTN)',
    'analysis':      '09:  病毒组下游分析 + 结果整合',
    'analysis_verify':'09b: 分析复核 (HMM / CDD-CT3 证据)',
    'report':        '10:  TSV汇总 + Sankey图 + 交互式HTML报告',
}

STAGE_ARGS = {
    "clean":         ["clean-data"],
    "deplete":       ["host_depletion"],
    "assembly":      ["assembly_pipeline"],
    "identification":["virus_identification"],
    "filter":        ["filter_virus"],
    "cobra":         ["cobra_pipeline"],
    "merge":         [],
    "cluster":       ["cluster_pipeline"],
    "taxonomy":      ["virus_classifier"],
    "host":          ["run_host_prediction"],
    "checkv":        [],
    "rescue":        ["cluster_pipeline"],
    "analysis":      [],
    "analysis_verify":[],
    "report":        [],
}

OVERVIEW = """
╔══════════════════════════════════════════════════════════════╗
║   MMPV-RNA — 宏病毒组端到端全自动分析流水线          ║
╚══════════════════════════════════════════════════════════════╝

Stage 流程 (顺序执行, 共15个, 与 STAGE_ORDER 一致):
  clean           → 原始数据质控 (fastp)
  deplete         → 宿主去除 (bowtie2)
  assembly        → 组装 (rnaviralSPAdes)
  identification  → 病毒鉴定 (10工具并行)
  filter          → 高置信过滤 (UniProt-strict)
  cobra           → 单样本延伸 (COBRA-Meta)
  merge           → 多样本共组装 (Flye, co-assembly 模式)
  cluster         → vOTU 聚类 (CD-HIT + vclust)
  taxonomy        → 分类学注释 (8工具加权投票)
  host            → 宿主预测 (三级决策树)
  checkv          → 完整性评估 (CheckV)
  rescue          → 四支路级联拯救 (CheckV → VSI → BLASTN/PlantVirusDB → genus_len)
  analysis        → 病毒组下游分析
  analysis_verify → 分析复核 (HMM/CT3 证据)
  report          → HTML 报告生成

用法:
  python virome_pipeline.py --stage rescue --output_dir <DIR>
  python virome_pipeline.py --stage rescue,analysis,report --output_dir <DIR>
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
    g.add_argument('--io-layout', choices=['legacy', 'standard'], default=None,
                   help='I/O 目录布局: legacy=v3.0 现行目录名 (默认); '
                        'standard=统一编号布局 (doc/IO_LAYOUT_DESIGN.md)。'
                        '亦可经环境变量 MMPV_IO_LAYOUT 设置, 子进程自动继承')
    g.add_argument('--cluster_input', help='直接输入已合并的病毒 FASTA (跳过 COBRA 收集)')
    g.add_argument('--input_assembly', help='组装结果目录 (identification 阶段, 默认 01_Assembly/)')

    g = p.add_argument_group('流程控制')
    g.add_argument('--config', default=None, help='YAML 配置文件路径 (默认: 自动查找 pipeline_config.yaml)')
    g.add_argument('--profile', default='default', help='配置预设 (默认: default, 可选: downstream/plant)')
    g.add_argument('--dump-config', action='store_true', help='仅打印配置摘要并退出 (不运行)')
    g.add_argument('--stage', default=['all'], nargs='+',
                   choices=['all', 'clean', 'deplete', 'assembly', 'identification', 'filter', 'cobra', 'merge', 'cluster', 'taxonomy', 'host', 'checkv', 'rescue', 'analysis', 'analysis_verify', 'report'],
                   help='运行阶段 (可多个, 如: --stage clean deplete)')
    g.add_argument('--host-filter', default='Plant',
                   help='目标宿主 (逗号分隔, rescue 阶段使用, 默认: Plant. Unknown 默认跳过并输出到 unknown_votus.fasta)')
    g.add_argument('--skip-sdt', action='store_true', help='跳过 09b 的 SDT 属级序列一致性矩阵绘制')
    g.add_argument('--sdt-min-ratio', type=float, default=0.0,
                   help='SDT 只跑 len_ratio>=该值的 KEEP 所属属 (0=全部, 0.7=仅接近属平均长度)')
    g.add_argument('--sdt-threads', type=int, default=8, help='SDT 逐对比对线程数')
    g.add_argument('--sdt-max-seqs', type=int, default=80,
                   help='SDT 属级矩阵单属序列上限 (超过则跳过, 防 O(N^2) 组合爆炸)')
    g.add_argument('--tree-min-ratio', type=float, default=0.7, help='建树 KEEP 最小 len_ratio (接近属平均长度比例, 0=不过滤; 属长缺失不过滤)')
    g.add_argument('--keep-min-len', type=int, default=0, help='KEEP 最小长度 bp (默认 0=不过滤; 因 RNA 数据里 DNA 病毒常只拼到基因级, 短 contig 也应保留。设 >0 才启用长度门; HMM 比对命中的仍豁免)')
    g.add_argument('--ct3-hmm-dir', default=None,
                   help='CT3 病毒 HMM 库目录 (默认 ~/database/virus-db/ct3_DBs/hmmscan_DBs/v3.1.1, 内含 RDRP/DNA_rep/Virion/Useful_Annotation/phrogs 五个 .h3m; 缺则 09b 第三路探针优雅跳过)')
    g.add_argument('--skip-ct3-hmm', action='store_true', help='跳过 09b 的 CT3-HMM 第三路探针')
    g.add_argument('--no-hmm-rescue', action='store_true',
                   help='关闭第三路探针参与裁决 (默认开: 无 CDD 病毒域但 CT3-HMM 有病毒命中时 DROP→REVIEW)')
    g.add_argument('--no-hmm-in-score', action='store_true',
                   help='关闭 HMM 进分 (默认开: 09b 整合时域证据层取 max(CDD, CT3-HMM), CDD 盲区由 HMM 补; '
                        '关闭则域证据层只看 CDD, 回退旧 30/30/40 纯 CDD 口径)')
    g.add_argument('--hmm-min-specific', type=int, default=1,
                   help='触发 HMM rescue 所需的最少病毒命中数 (默认 1)')
    g.add_argument('--no-hmm-promote-keep', action='store_true',
                   help='关闭 CT3-HMM 升 KEEP (默认开: 有非噬菌体病毒命中时 DROP/REVIEW 直接升 KEEP)')
    g.add_argument('--hmm-promote-min', type=int, default=1,
                   help='触发升 KEEP 所需的最少病毒命中数 (默认 1; 配合 --hmm-require-family 时按科级计数)')
    g.add_argument('--hmm-require-family', action='store_true',
                   help='收紧 HMM rescue 口径: 只认命中 ICTV 分类单元的模型 (ct3_n_family)')
    g.add_argument('--hmm-count-ambiguous', action='store_true',
                   help='放宽 HMM rescue 口径: 把 VOG 等无法判定的 AMBIGUOUS 命中也算进病毒证据')
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
    g.add_argument('--kraken2_confidence', type=float, default=0.2, help='Kraken2 分类置信度阈值 (默认: 0.2)')
    g.add_argument('--keep_rrna', action='store_true', help='保留分离出的 rRNA 序列到 rrna/ 目录')
    g.add_argument('--keep_taxids', default='10239', help='Kraken2 保留的 taxID (逗号分隔), 默认 10239=Viruses, 设空字符串关闭')
    g.add_argument('--vmtouch', action='store_true', help='用 vmtouch 预热 Kraken2 数据库到 page cache (配合 --memory-mapping, 多样本并行时推荐)')
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
    g.add_argument('--nr_db', default=None, help='Diamond NR 数据库路径 (默认: 从 pipeline_config.yaml 读取)')
    g.add_argument('--skip_uniprot_filter', action='store_true', help='跳过 UniProt 后置过滤')
    g.add_argument('--skip_nr_filter', action='store_true', help='跳过 NR 后置过滤')
    g.add_argument('--skip_id_plots', action='store_true', help='跳过鉴定阶段图表生成')
    g.add_argument('--clean_failed', action='store_true', help='自动清理鉴定失败的任务目录')
    g.add_argument('--ident_ext', default='.contig.fasta', help='输入目录时搜索的后缀 (默认: .contig.fasta, 避开 scaffold 避免重复鉴定)')

    g = p.add_argument_group('Filter 阶段 (filter_virus.py)')
    g.add_argument('--filter-mode', choices=['filter','strict','raw'], default='filter',
                   help='过滤模式: raw=不过滤, filter=UniProt(关键词抢救)+CDD(viral+viral_plant+viral_bacteria), strict=UniProt(仅TaxID)+CDD(仅viral)')
    g.add_argument('--no-merge', action='store_true', help='逐样本过滤 (默认合并后拆分)')

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
    g.add_argument('--flye-min-overlap', type=int, default=1000,
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
    g.add_argument('--host_db', default=None, help='宿主数据库根目录 (默认: 从 pipeline_config.yaml 读取)')
    g.add_argument('--kraken2_db', help='Kraken2 宿主库 (覆盖 --host_db 自动检测)')
    g.add_argument('--host_align_db', help='宿主比对索引 (覆盖 --host_db 自动检测)')
    g.add_argument('--virus_db', default=None, help='病毒鉴定数据库根目录 (默认: 从 pipeline_config.yaml 读取)')
    g.add_argument('--checkv_db', default=None, help='CheckV 数据库路径 (默认: 从 pipeline_config.yaml 读取)')
    g.add_argument('--blast-db', default=None, help='BLAST 参考数据库 (rescue 阶段, 默认: 从 pipeline_config.yaml 读取)')

    g = p.add_argument_group('工具与算法')
    g.add_argument('--aligner', default='bowtie2', choices=['bowtie2', 'hisat2', 'minimap2'])
    g.add_argument('--seq_type', default='rna-short', choices=['dna-short', 'rna-short', 'nanopore', 'pacbio'])
    g.add_argument('--rrna', action='store_true', help='开启 rRNA 剔除')
    g.add_argument('--rrna_tool', default='ribodetector', choices=['ribodetector', 'silva'],
                   help='rRNA 剔除工具: ribodetector (默认) / silva (Bowtie2+SILVA)')
    g.add_argument('--silva_index', help='SILVA Bowtie2 索引前缀 (--rrna_tool silva 时必需)')
    g.add_argument('--coassembly', action='store_true', help='Co-assembly 模式: 合并所有样本 reads 进行单次组装')
    g.add_argument('--assembler', default=None, choices=['megahit', 'rnaviralspades', 'penguin', 'all'],
                    help='组装工具 (默认: 从 pipeline_config.yaml 读取, 未配置则用 rnaviralspades)')
    g.add_argument('--contig-length', '-l', type=int, default=200, help='contig 最小长度 bp (默认 200)')
    g.add_argument('--identify_tools', default='all', help='病毒鉴定工具')
    g.add_argument('--virus_mode', default='filter', choices=['raw', 'filter', 'strict'],
                   help='COBRA 病毒序列来源: raw=原始鉴定, filter=UniProt过滤, strict=严格过滤 (默认: strict)')
    g = p.add_argument_group('鉴定数据库 (virus_identification.py)')
    g.add_argument('--virus_protein_db', help='病毒蛋白 Diamond DB')
    g.add_argument('--uniprot_db', help='UniProt Diamond DB')
    g.add_argument('--viroids_db', default=os.path.expanduser('~/database/virus-db/viroids-db/viroids.fasta.blast.db'),
                   help='类病毒 BLAST DB FASTA')
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
    g.add_argument('--prob-dir', help='ICTV 宿主概率表目录 (host 阶段, 默认: database/cross_analysis/)')
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


def _expand_env_config(obj):
    """递归展开配置里的 ${VAR} / ${VAR:-默认}（换环境可用环境变量覆盖, 不必改 yaml）。

    与 virome_phylo_pipeline/utils/dataset_config.py 语义一致：
    默认值允许一层 {占位符} 嵌套；变量未定义且无默认 → 空串。
    """
    pat = re.compile(r"\$\{([A-Za-z_]\w*)(?::-((?:[^{}]|\{[^{}]*\})*))?\}")

    def _s(v):
        def _sub(m):
            name, default = m.group(1), m.group(2)
            val = os.environ.get(name)
            if val:
                return val
            return default if default is not None else ""

        # 迭代展开以支持嵌套（如 ${A:-${B:-x}/y}）：re.sub 单次替换不会
        # 再扫描替换结果，故需多轮直到稳定
        for _ in range(10):
            new = pat.sub(_sub, v)
            if new == v:
                break
            v = new
        return v

    if isinstance(obj, str):
        return _s(obj)
    if isinstance(obj, dict):
        return {k: _expand_env_config(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env_config(v) for v in obj]
    return obj


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
    # 展开 ${VAR} / ${VAR:-默认}（换环境可用环境变量覆盖数据库/工具路径, 不必改 yaml）
    config = _expand_env_config(config)

    profiles = config.get('profiles', {})
    profile = profiles.get(args.profile, profiles.get('default', {}))
    if not profile:
        print(f"[WARN] 配置 profile '{args.profile}' 未找到, 使用 CLI 参数")
        return args

    # 数据库路径
    db = profile.get('databases', {})
    for key in ['checkv_db','genomad_db','mmseqs_db','virus_db','virus_protein_db','uniprot_db',
                'host_db','blast_db','nr_db','db_dir',
                'viralverify_hmm','virsorter_db','metabuli_db','viroids_db','virus_taxid']:
        if getattr(args, key, None) is None and key in db:
            setattr(args, key, db[key])

    # 工具路径
    tools = profile.get('tools', {})
    tool_path_map = {'salmon': 'salmon_bin', 'virbot': 'virbot_path', 'virhunter': 'virhunter_path', 'viralm': 'viralm_path'}
    for key in ['salmon','diamond','ragtag','virbot','virhunter','viralm']:
        ak = tool_path_map.get(key, key)
        if ak and getattr(args, ak, None) is None and key in tools:
            setattr(args, ak, os.path.expanduser(tools[key]))

    # 运行参数: CLI 显式指定 > YAML profile > argparse 默认值
    rt = profile.get('runtime', {})
    _cli_set = set()
    for i, a in enumerate(sys.argv):
        if a in ('-t', '--threads'): _cli_set.add('threads')
        if a in ('-j', '--jobs'): _cli_set.add('jobs')
        if a == '--tax_jobs': _cli_set.add('tax_jobs')
        if a == '--assembler': _cli_set.add('assembler')
        if a == '--rrna': _cli_set.add('rrna')
    for key in ['threads','jobs','tax_jobs']:
        if key in rt and key not in _cli_set:  # CLI 没显式传 → 用 YAML 值
            setattr(args, key, int(rt[key]))

    # assembly
    asm_cfg = profile.get('assembly', {})
    _asm_cli = 'assembler' in _cli_set
    if 'assembler' in asm_cfg and not _asm_cli:
        setattr(args, 'assembler', asm_cfg['assembler'])
    if getattr(args, 'assembler', None) is None:
        setattr(args, 'assembler', 'rnaviralspades')  # 最终兜底

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

    # 按需检查: 只验证当前阶段需要的数据库
    stages = set(getattr(args, 'stage', ['all']))
    _all = 'all' in stages
    need_id = _all or bool(stages & {'identification', 'taxonomy', 'host'})
    need_deplete = _all or bool(stages & {'clean', 'deplete'})
    need_checkv = _all or bool(stages & {'checkv', 'rescue'})

    db_keys = []
    if need_deplete:
        db_keys += ['host_db']
    if need_id:
        db_keys += ['virus_db', 'genomad_db', 'virus_protein_db', 'uniprot_db', 'nr_db', 'mmseqs_db']
    if need_checkv:
        db_keys += ['checkv_db', 'blast_db']
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


def _record_provenance(args, logger):
    """
    记录运行环境快照到 provenance.json。
    纯附加: 任何步骤失败都不影响管道执行。
    """
    import datetime
    import platform
    import importlib.metadata

    prov = {
        "pipeline": "MMPV-RNA v2.3",
        "timestamp": datetime.datetime.now().isoformat(),
        "hostname": platform.node(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "tools": {},
        "databases": {},
        "packages": {},
    }

    # ── 工具版本探测 ──
    TOOL_VERSION_ARGS = {
        "checkv":       ["checkv", "--version"],
        "genomad":      ["genomad", "--version"],
        "blastn":       ["blastn", "-version"],
        "diamond":      ["diamond", "--version"],
        "kraken2":      ["kraken2", "--version"],
        "bowtie2":      ["bowtie2", "--version"],
        "megahit":      ["megahit", "--version"],
        "seqkit":       ["seqkit", "version"],
        "fastp":        ["fastp", "--version"],
        "salmon":       ["salmon", "--version"],
        "mmseqs":       ["mmseqs", "version"],
        "samtools":     ["samtools", "--version"],
        "vclust":       ["vclust", "--version"],
        "bwa-mem2":     ["bwa-mem2", "version"],
        "minimap2":     ["minimap2", "--version"],
        "rnaviralspades": ["rnaspades.py", "--version"],
        "ragtag":       ["ragtag.py", "--version"],
    }

    for tool, cmd in TOOL_VERSION_ARGS.items():
        try:
            exe = shutil.which(cmd[0])
            if not exe:
                continue
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            out = (result.stdout or "").strip()
            err = (result.stderr or "").strip()
            # 很多工具把版本号打到 stderr
            version_line = (out or err).split("\n")[0][:200]
            prov["tools"][tool] = version_line
            logger.debug("  provenance: %s → %s", tool, version_line[:80])
        except Exception:
            pass

    # ── 数据库元数据 ──
    db_keys = ["checkv_db", "genomad_db", "mmseqs_db", "virus_db",
               "virus_protein_db", "uniprot_db", "blast_db", "nr_db", "host_db"]
    for key in db_keys:
        val = getattr(args, key, None)
        if not val:
            continue
        try:
            p = Path(os.path.expanduser(str(val)))
            if p.exists():
                st = p.stat()
                prov["databases"][key] = {
                    "path": str(p),
                    "size_gb": round(st.st_size / (1024**3), 2) if p.is_file() else None,
                    "modified": datetime.datetime.fromtimestamp(st.st_mtime).isoformat(),
                }
            elif p.parent.is_dir() and key == "blast_db":
                # BLAST DB 是多个文件
                prov["databases"][key] = {"path": str(p), "note": "multi-file BLAST DB"}
        except Exception:
            pass

    # ── ICTV MSL 版本检测 ──
    try:
        msl_candidates = [
            Path(os.path.expanduser("~/database/virus-db/VMR_MSL41.v1.20260320.xlsx")),
            Path(__file__).resolve().parent.parent / "database" / "VMR_MSL41.v1.20260320.xlsx",
        ]
        for msl in msl_candidates:
            if msl.is_file():
                m = re.search(r"MSL(\d+)", msl.name)
                prov["ictv_msl"] = f"MSL{m.group(1)}" if m else msl.name
                prov["ictv_msl_file"] = str(msl)
                break
    except Exception:
        pass

    # ── Python 包版本 ──
    PKG_NAMES = ["biopython", "pandas", "numpy", "polars", "matplotlib",
                 "seaborn", "scipy", "Bio", "tqdm", "psutil"]
    for pkg in PKG_NAMES:
        try:
            prov["packages"][pkg] = importlib.metadata.version(pkg)
        except Exception:
            try:
                mod = __import__(pkg)
                prov["packages"][pkg] = getattr(mod, "__version__", "installed")
            except Exception:
                pass

    # ── conda / pixi 环境信息 ──
    try:
        result = subprocess.run(["conda", "info", "--json"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            import json as _json
            info = _json.loads(result.stdout)
            prov["conda_env"] = info.get("active_prefix_name", "unknown")
    except Exception:
        pass

    try:
        pixi_toml = Path(__file__).resolve().parent.parent / "pixi.toml"
        if pixi_toml.is_file():
            prov["pixi_toml"] = str(pixi_toml)
    except Exception:
        pass

    # ── 写入 provenance.json ──
    try:
        prov_path = Path(args.output_dir) / "provenance.json"
        prov_path.parent.mkdir(parents=True, exist_ok=True)
        with open(prov_path, "w", encoding="utf-8") as pf:
            json.dump(prov, pf, indent=2, ensure_ascii=False)
        logger.info("  provenance.json → %s (%d tools, %d databases, %d packages)",
                    prov_path, len(prov["tools"]), len(prov["databases"]), len(prov["packages"]))
    except Exception as e:
        logger.debug("  provenance.json 写入失败: %s", e)


def main():
    args = parse_args()
    args = _load_config(args)
    logger = setup_logger(args.output_dir, level=args.log_level)
    _validate_config(args, logger)
    _record_provenance(args, logger)

    # 进程守护: 停止编排器时连带杀掉整条流水线的所有子孙进程 (各阶段脚本及其 spades/kraken2 等)
    try:
        import process_guard
        process_guard.install(logger)
    except Exception:
        pass

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
        logger.info("  阶段:     %s", ','.join(sorted(stages)))
        logger.info("  输出目录: %s", args.output_dir)
        logger.info("═══ DRY-RUN 结束 ═══")
        return

    pipe = ViromePipeline(args, logger)

    needs_reads = 'all' in stages or bool(stages & {'clean','deplete','assembly','cobra','rescue'})
    if needs_reads:
        # 优先使用 --input_reads 指定的目录, 但 rescue 阶段推荐去宿主产物目录
        if args.input_reads and Path(args.input_reads).exists():
            pipe.reads_dir = Path(args.input_reads).absolute()
            hostdep_dir = Path(args.output_dir) / layout_dir_name('d_hostdep', pipe.layout)
            if 'rescue' in stages and hostdep_dir.is_dir() and hostdep_dir.resolve() != pipe.reads_dir.resolve():
                logger.warning("  ⚠ rescue 阶段建议 --input_reads out/%s/ (当前: %s)",
                               layout_dir_name('d_hostdep', pipe.layout), pipe.reads_dir)
        elif 'deplete' in stages:
            # deplete: 自动使用 clean 输出 (优先 clumpify, 否则 fasta)
            cl = pipe.d['clean'] / '3.clumpify'
            fa = pipe.d['clean'] / '2.fasta'
            pipe.reads_dir = cl if (cl.exists() and any(cl.iterdir())) else fa
        else:
            pipe.reads_dir = pipe.d['hostdep']
        if not pipe.reads_dir.exists() or not any(pipe.reads_dir.iterdir()):
            if stages <= {'merge', 'cluster', 'taxonomy', 'host', 'checkv', 'rescue',
                          'analysis', 'analysis_verify', 'report'}:
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
        'filter': pipe.run_filter,
        'cobra': pipe.run_cobra, 'merge': pipe.run_merge_samples, 'cluster': pipe.run_cluster,
        'taxonomy': pipe.run_taxonomy, 'host': pipe.run_host,
        'checkv': pipe.run_checkv_stage, 'rescue': pipe.run_rescue,
        'analysis': pipe.run_analysis, 'analysis_verify': pipe.run_analysis_verify,
        'report': pipe.run_reports,
    }
    # 按流水线顺序排列
    stage_order = ['clean','deplete','assembly','identification','filter','cobra','merge','cluster',
                   'taxonomy','host','checkv','rescue','analysis','analysis_verify','report']
    # co-assembly 模式: 03a(COBRA) 无逐样本延伸意义, 03b(merge+Flye) 也一并跳过
    # (cluster 的 coassembly 兑底已改从 02a/02b 收集 ALL_merged 候选)
    if getattr(args, 'coassembly', False):
        skipped_03 = [s for s in ('cobra', 'merge') if _all or s in stages]
        if skipped_03:
            logger.info("[co-assembly] 跳过 Stage 03 (COBRA + merge): %s", ','.join(skipped_03))
        stage_order = [s for s in stage_order if s not in ('cobra', 'merge')]
    stages_to_run = [(s, stage_map[s]) for s in stage_order if _all or s in stages]

    # 准备阶段日志目录
    stage_log_dir = Path(args.output_dir) / layout_dir_name('d_reports') / "logs"
    stage_log_dir.mkdir(parents=True, exist_ok=True)

    failed_stages = []
    for stage_name, stage_func in stages_to_run:
        # 断点续传: 已完成则跳过
        if pipe._stage_skip(stage_name):
            # 跳过时仍需更新 reads_dir 指针
            if stage_name == 'clean':
                cl = pipe.d['clean'] / '3.clumpify'
                fa = pipe.d['clean'] / '2.fasta'
                pipe.reads_dir = cl if (cl.exists() and any(cl.iterdir())) else fa
            elif stage_name == 'deplete' and pipe.d['hostdep'].is_dir():
                pipe.reads_dir = pipe.d['hostdep']
            continue
        # 添加阶段独立日志 handler
        stage_handler = logging.FileHandler(str(stage_log_dir / f"{stage_name}.log"), encoding='utf-8')
        stage_handler.setLevel(logging.DEBUG)
        stage_handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s', datefmt='%H:%M:%S'))
        logger.addHandler(stage_handler)
        try:
            stage_func()
            pipe._stage_done(stage_name)
            # 更新 reads_dir 指针
            if stage_name == 'clean':
                cl = pipe.d['clean'] / '3.clumpify'
                fa = pipe.d['clean'] / '2.fasta'
                pipe.reads_dir = cl if (cl.exists() and any(cl.iterdir())) else fa
            elif stage_name == 'deplete' and pipe.d['hostdep'].is_dir() and any(pipe.d['hostdep'].iterdir()):
                pipe.reads_dir = pipe.d['hostdep']
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
