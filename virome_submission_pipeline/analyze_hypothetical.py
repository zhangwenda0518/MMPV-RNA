#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_hypothetical.py — 假定蛋白功能注释工具 v2.1
====================================================

从 suvtk features 产出的 featuretable.tbl 中识别 hypothetical protein，
提取对应蛋白序列，运行 blastp 获取功能注释，并更新特征表和蛋白文件。

流程:
  1. 解析 featuretable.tbl → 找出所有 hypothetical protein CDS
  2. 匹配 proteins.faa → 提取对应蛋白序列
  3. 运行 blastp → 搜索功能同源
  4. 输出 blast 结果
  5. 输出更新后的 featuretable 和 proteins 文件

用法:
  python analyze_hypothetical.py \\
      -t featuretable.tbl \\
      -f proteins.faa \\
      --blast my_blast_results.txt \\
      -d /path/to/blast/db \\
      -o ./output/ \\
      --threads 8
"""

import argparse
import os
import re
import sys
import subprocess
import logging
import tempfile
import time
from collections import OrderedDict
from pathlib import Path

from tqdm import tqdm
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez


# ══════════════════════════════════════════════════════════════
# 日志
# ══════════════════════════════════════════════════════════════

def setup_logger():
    logger = logging.getLogger("AnalyzeHypo")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)
    return logger

LOG = setup_logger()


# ══════════════════════════════════════════════════════════════
# 数据结构
# ══════════════════════════════════════════════════════════════

class CdsEntry:
    """featuretable.tbl 中的一条 CDS 记录"""
    __slots__ = ('seq_id', 'start', 'end', 'strand',
                 'cds_index', 'product', 'inferences', 'notes')

    def __init__(self, seq_id, start, end, strand, cds_index):
        self.seq_id = seq_id
        self.start = start
        self.end = end
        self.strand = strand
        self.cds_index = cds_index
        self.product = None
        self.inferences = []
        self.notes = []

    @property
    def coord_min(self):
        return min(self.start, self.end)

    @property
    def coord_max(self):
        return max(self.start, self.end)

    @property
    def is_hypothetical(self):
        return self.product and 'hypothetical' in self.product.lower()

    def to_tbl_lines(self):
        """生成该 CDS 的 .tbl 文本行"""
        lines = []
        # 位置行: 负链先写大坐标
        if self.strand == -1:
            lines.append(f"{self.end}\t{self.start}\tCDS")
        else:
            lines.append(f"{self.start}\t{self.end}\tCDS")
        # product
        lines.append(f"\t\t\tproduct\t{self.product or 'hypothetical protein'}")
        # inference
        for inf in self.inferences:
            lines.append(f"\t\t\tinference\t{inf}")
        # notes
        for note in self.notes:
            lines.append(f"\t\t\tnote\t{note}")
        return lines


class ProteinEntry:
    """proteins.faa 中的一条蛋白记录"""
    __slots__ = ('header', 'seq_id', 'orig_seq_id', 'start', 'end',
                 'strand', 'record_idx', 'cds_idx', 'sequence')

    def __init__(self, header, seq_id, orig_seq_id, start, end, strand,
                 record_idx, cds_idx, sequence):
        self.header = header
        self.seq_id = seq_id
        self.orig_seq_id = orig_seq_id
        self.start = start
        self.end = end
        self.strand = strand
        self.record_idx = record_idx
        self.cds_idx = cds_idx
        self.sequence = sequence

    @property
    def coord_min(self):
        return min(self.start, self.end)

    @property
    def coord_max(self):
        return max(self.start, self.end)


class BlastHit:
    """一条 blastp 命中结果"""
    __slots__ = ('query', 'sseqid', 'pident', 'evalue', 'bitscore', 'stitle')

    def __init__(self, query, sseqid, pident, evalue, bitscore, stitle):
        self.query = query
        self.sseqid = sseqid
        self.pident = float(pident)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
        self.stitle = stitle


# ══════════════════════════════════════════════════════════════
# 步骤 1: 解析 featuretable.tbl
# ══════════════════════════════════════════════════════════════

def _count_lines(filepath):
    """快速统计文件行数"""
    count = 0
    with open(filepath, 'rb') as f:
        for _ in f:
            count += 1
    return count


def parse_featuretable(tbl_path):
    """解析 featuretable.tbl，返回 {seq_id: [CdsEntry, ...]}"""
    if not os.path.exists(tbl_path):
        LOG.error("featuretable.tbl 不存在: %s", tbl_path)
        sys.exit(1)

    total_lines = _count_lines(tbl_path)
    records = OrderedDict()
    current_seq = None
    current_cds = None
    seq_count = 0
    cds_total = 0

    with open(tbl_path, 'r', encoding='utf-8') as f:
        pbar = tqdm(f, total=total_lines, desc="[1] 解析 featuretable.tbl",
                    unit="lines", ncols=100, bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]')

        for raw_line in pbar:
            line = raw_line.rstrip('\n')

            # >Feature 行 → 新序列
            if line.startswith('>Feature '):
                current_seq = line.split('>Feature ', 1)[1].strip()
                records.setdefault(current_seq, [])
                current_cds = None
                seq_count += 1
                continue

            # CDS 行: 数字 tab 数字 tab CDS
            stripped = line.strip()
            parts = stripped.split()
            if len(parts) >= 3 and parts[2] == 'CDS':
                current_cds = None
                try:
                    s = int(parts[0].replace('>', '').replace('<', ''))
                    e = int(parts[1].replace('>', '').replace('<', ''))
                except ValueError:
                    continue

                strand = 1 if s <= e else -1
                idx = len(records[current_seq]) + 1
                current_cds = CdsEntry(current_seq, s, e, strand, idx)
                records[current_seq].append(current_cds)
                cds_total += 1
                continue

            # 特征属性行 (三个 tab 开头)
            if current_cds and line.startswith('\t\t\t'):
                attr = stripped
                if attr.startswith('product\t'):
                    current_cds.product = attr.split('\t', 1)[1] if '\t' in attr else ''
                elif attr.startswith('inference\t'):
                    current_cds.inferences.append(attr.split('inference\t', 1)[1])
                elif attr.startswith('note\t'):
                    current_cds.notes.append(attr.split('note\t', 1)[1])

    hypo_count = sum(1 for v in records.values() for c in v if c.is_hypothetical)
    LOG.info("[1] 解析 featuretable.tbl: %d 条序列, %d 个 CDS, %d 个 hypothetical",
             seq_count, cds_total, hypo_count)
    return records


# ══════════════════════════════════════════════════════════════
# 步骤 2: 解析 proteins.faa
# ══════════════════════════════════════════════════════════════

def parse_proteins_faa(faa_path):
    """解析 proteins.faa，返回 {(seq_id, coord_min, coord_max): ProteinEntry} 字典

    头部格式:
      >{seq_id} # {start} # {end} # {strand} # ID={record}_{cds};...
    """
    if not os.path.exists(faa_path):
        LOG.error("proteins.faa 不存在: %s", faa_path)
        sys.exit(1)

    total_lines = _count_lines(faa_path)
    proteins = {}
    header, seq_lines = None, []
    parsed_count = 0

    with open(faa_path, 'r', encoding='utf-8') as f:
        pbar = tqdm(f, total=total_lines, desc="[2] 解析 proteins.faa  ",
                    unit="lines", ncols=100, bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]')

        for raw_line in pbar:
            line = raw_line.strip()
            if line.startswith('>'):
                if header:
                    prot = _parse_faa_header(header, ''.join(seq_lines))
                    if prot:
                        key = (prot.seq_id, prot.coord_min, prot.coord_max)
                        proteins[key] = prot
                        parsed_count += 1
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # 最后一条
        if header:
            prot = _parse_faa_header(header, ''.join(seq_lines))
            if prot:
                key = (prot.seq_id, prot.coord_min, prot.coord_max)
                proteins[key] = prot
                parsed_count += 1

    LOG.info("[2] 解析 proteins.faa: %d 条蛋白序列", parsed_count)
    return proteins


def _parse_faa_header(header, sequence):
    """解析单条 proteins.faa 头部

    faa 头部格式: >{record_id}_{gene_index} # {start} # {end} # {strand} # ID=X_Y;...
    其中 pyrodigal 在 record_id 后追加了 _{gene_index} 后缀
    """
    h = header.lstrip('>')

    parts = h.split(' # ')
    if len(parts) < 4:
        return None

    full_seq_id = parts[0].strip()

    # 解析 ID=X_Y
    record_idx, cds_idx = 1, 1
    for part in parts:
        if part.startswith('ID='):
            id_val = part.split(';')[0].replace('ID=', '')
            try:
                r, c = id_val.split('_')
                record_idx = int(r)
                cds_idx = int(c)
            except ValueError:
                pass
            break

    # 提取原始 record_id: 从末尾删除 pyrodigal 追加的 _{cds_idx} 后缀
    suffix = f"_{cds_idx}"
    if full_seq_id.endswith(suffix):
        orig_seq_id = full_seq_id[:-len(suffix)]
    else:
        orig_seq_id = re.sub(r'_\d+$', '', full_seq_id)

    try:
        start = int(parts[1].strip())
        end = int(parts[2].strip())
        strand = int(parts[3].strip())
    except (ValueError, IndexError):
        return None

    return ProteinEntry(header.lstrip('>'), full_seq_id, orig_seq_id,
                        start, end, strand, record_idx, cds_idx, sequence)


# ══════════════════════════════════════════════════════════════
# 步骤 3: 匹配 hypothetical protein → 蛋白序列
# ══════════════════════════════════════════════════════════════

def match_hypotheticals(records, proteins):
    """将 featuretable 中的 hypothetical CDS 匹配到 proteins.faa 中的序列

    匹配: (seq_id, coord_min, coord_max) — seq_id + CDS 坐标区间

    返回: [{seq_id, cds_index, protein_header, protein_sequence, ...}, ...]
    """
    # 构建索引: (orig_seq_id, coord_min, coord_max) → ProteinEntry
    by_ocoord = {}
    for prot in tqdm(proteins.values(), desc="[3] 构建蛋白索引    ",
                     unit="seqs", ncols=100,
                     bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]'):
        by_ocoord[(prot.orig_seq_id, prot.coord_min, prot.coord_max)] = prot

    matched = []
    unmatched_cds = []

    # 收集所有 hypothetical CDS
    all_hypo = []
    for seq_id, cdses in records.items():
        for cds in cdses:
            if cds.is_hypothetical:
                all_hypo.append(cds)

    for cds in tqdm(all_hypo, desc="[3] 匹配 hypothetical  ",
                    unit="cds", ncols=100,
                    bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]'):
        prot = by_ocoord.get((cds.seq_id, cds.coord_min, cds.coord_max))

        if prot:
            matched.append({
                'seq_id': cds.seq_id,
                'cds_index': cds.cds_index,
                'start': cds.start,
                'end': cds.end,
                'strand': cds.strand,
                'product': cds.product,
                'protein_header': prot.header,
                'protein_sequence': prot.sequence,
            })
        else:
            unmatched_cds.append(f"{cds.seq_id}[CDS#{cds.cds_index}]")

    LOG.info("[3] 匹配 hypothetical protein: %d/%d 成功",
             len(matched), len(matched) + len(unmatched_cds))
    if unmatched_cds:
        LOG.warning("  未匹配 %d 条: %s ...", len(unmatched_cds),
                     ", ".join(unmatched_cds[:5]))
    return matched


# ══════════════════════════════════════════════════════════════
# 步骤 4a: 运行在线 BLAST (NCBI)
# ══════════════════════════════════════════════════════════════

def run_online_blast(hypotheticals, output_file, database="nr", delay=3,
                     email=None, api_key=None):
    """对 hypothetical 蛋白序列逐个提交 NCBI 在线 BLAST

    参数:
      hypotheticals : 匹配后的 hypothetical 蛋白列表
      output_file  : 结果输出文件 (与本地 blast --outfmt 6 格式一致)
      database     : NCBI 数据库名 (默认 nr)
      delay        : 每次提交间隔秒数 (默认 3s; 有 API key 可降至 1s)
      email        : NCBI 要求, 用于追踪和限流
      api_key      : NCBI API key, 提升请求频率上限 (3/s → 10/s)
    """
    if not hypotheticals:
        LOG.warning("[4] 没有 hypothetical protein 需要 BLAST, 跳过")
        return

    # 配置 NCBI 认证
    if email:
        Entrez.email = email
        LOG.info("[4] NCBI 邮箱: %s", email)
    if api_key:
        Entrez.api_key = api_key
        LOG.info("[4] NCBI API key: %s...%s", api_key[:4], api_key[-4:])

    total = len(hypotheticals)
    LOG.info("[4] NCBI 在线 BLAST: %d 条序列 → database=%s (间隔=%ds, 预计耗时 ~%d min)",
             total, database, delay, total * (delay + 15) // 60)

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)

    blast_hits = []  # 收集命中结果
    fail_count = 0

    pbar = tqdm(enumerate(hypotheticals, 1), total=total,
                desc="[4] NCBI 在线 BLAST ",
                unit="query", ncols=100,
                bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]')

    for i, hp in pbar:
        seq_id = hp['seq_id']
        seq = hp['protein_sequence']

        pbar.set_postfix_str(f"query={seq_id[:30]}...")

        try:
            # 提交 BLAST
            result_handle = NCBIWWW.qblast(
                program="blastp",
                database=database,
                sequence=seq,
                expect=1e-3,
                hitlist_size=5,
            )

            # 解析结果
            blast_record = NCBIXML.read(result_handle)

            if blast_record.alignments:
                for alignment in blast_record.alignments[:3]:
                    hsp = alignment.hsps[0]
                    identity_pct = (hsp.identities / hsp.align_length * 100) if hsp.align_length > 0 else 0

                    blast_hits.append({
                        'query': hp['protein_header'].lstrip('>'),
                        'sseqid': alignment.accession or alignment.hit_id.split('|')[1] if '|' in alignment.hit_id else alignment.hit_id,
                        'pident': f"{identity_pct:.2f}",
                        'length': str(hsp.align_length),
                        'mismatch': str(hsp.align_length - hsp.identities),
                        'gapopen': str(hsp.gaps or 0),
                        'qstart': str(hsp.query_start),
                        'qend': str(hsp.query_end),
                        'sstart': str(hsp.sbjct_start),
                        'send': str(hsp.sbjct_end),
                        'evalue': f"{hsp.expect:.1e}",
                        'bitscore': f"{hsp.score:.1f}",
                        'stitle': alignment.title,
                    })
            else:
                LOG.debug("  %s: 无命中", seq_id[:50])

        except Exception as e:
            fail_count += 1
            LOG.warning("  %s: BLAST 失败 - %s", seq_id[:50], str(e)[:100])

        # NCBI 限速控制
        if i < total:
            time.sleep(delay)

    pbar.close()

    # 写 tabular 输出 (与 blastp -outfmt 6 格式一致)
    with open(output_file, 'w', encoding='utf-8') as f:
        for hit in blast_hits:
            f.write('\t'.join([
                hit['query'], hit['sseqid'], hit['pident'], hit['length'],
                hit['mismatch'], hit['gapopen'], hit['qstart'], hit['qend'],
                hit['sstart'], hit['send'], hit['evalue'], hit['bitscore'],
                hit['stitle']
            ]) + '\n')

    LOG.info("[4] 在线 BLAST 完成: %d 条命中, %d 条失败 → %s",
             len(blast_hits), fail_count, output_file)

# ══════════════════════════════════════════════════════════════
# 步骤 4b: 运行 diamond blastp (高速本地比对)
# ══════════════════════════════════════════════════════════════

def _ensure_diamond_db(db_path, threads, diamond_bin):
    """检查 diamond 数据库是否存在, 若不存在自动从同名前缀的 FASTA 构建

    db_path 可以是:
      - .dmnd 文件路径
      - .fasta/.fa 文件路径 (自动推导 .dmnd 路径)
    """
    # 标准化: 计算 dmnd 路径和 fasta 路径
    if db_path.endswith('.dmnd'):
        dmnd_path = db_path
        fasta_path = None
        # 尝试从 dmnd 路径反推 fasta
        base = db_path[:-5]  # 去掉 .dmnd
        for ext in ['.fasta', '.fa', '.faa', '.EX.acc.fasta']:
            candidate = base + ext
            if os.path.exists(candidate):
                fasta_path = candidate
                break
    elif db_path.endswith(('.fasta', '.fa', '.faa')):
        fasta_path = db_path
        dmnd_path = re.sub(r'\.(fasta|fa|faa)$', '.dmnd', db_path)
    else:
        # 无后缀, 加上 .dmnd 试试
        dmnd_path = db_path + '.dmnd'
        fasta_path = None
        for ext in ['.fasta', '.fa', '.faa', '.EX.acc.fasta']:
            candidate = db_path + ext
            if os.path.exists(candidate):
                fasta_path = candidate
                break

    if os.path.exists(dmnd_path):
        return dmnd_path

    # dmnd 不存在, 需要从 fasta 构建
    if fasta_path is None:
        # 尝试找同名前缀的 fasta
        base = dmnd_path.replace('.dmnd', '')
        for ext in ['.fasta', '.fa', '.faa', '.EX.acc.fasta']:
            if os.path.exists(base + ext):
                fasta_path = base + ext
                break
            # 也尝试去掉 .EX.acc 等后缀
            for sfx in ['.EX.acc', '.prot']:
                candidate = base.replace(sfx, '') + ext
                if os.path.exists(candidate):
                    fasta_path = candidate
                    break

    if fasta_path is None or not os.path.exists(fasta_path):
        LOG.error("[4] diamond 数据库不存在, 且未找到 FASTA 源文件")
        LOG.error("    dmnd 路径: %s", dmnd_path)
        return None

    LOG.info("[4] diamond 数据库不存在, 自动构建: %s → %s", fasta_path, dmnd_path)
    cmd = (
        f"{diamond_bin} makedb "
        f"--in {fasta_path} "
        f"-d {dmnd_path.replace('.dmnd', '')} "
        f"--threads {threads} "
        f">/tmp/diamond_makedb.log 2>&1"
    )
    ret = os.system(cmd)
    if ret != 0 or not os.path.exists(dmnd_path):
        LOG.error("[4] diamond makedb 失败 (exit=%d)", ret)
        return None

    LOG.info("[4] diamond 数据库构建完成: %s", dmnd_path)
    return dmnd_path


def run_diamond(hypotheticals, db_path, output_file, threads=8, evalue=1e-3,
                diamond_bin="diamond"):
    """对 hypothetical 蛋白序列运行 diamond blastp (自动构建 .dmnd 如不存在)"""
    if not hypotheticals:
        LOG.warning("[4] 没有 hypothetical protein 需要 BLAST, 跳过")
        return

    # 自动检测/构建 diamond 数据库
    db_path = _ensure_diamond_db(db_path, threads, diamond_bin)
    if db_path is None:
        return

    # 写入临时 FASTA
    fd, tmp_faa = tempfile.mkstemp(suffix='.faa', prefix='hypo_')
    os.close(fd)
    with open(tmp_faa, 'w', encoding='utf-8') as f:
        for hp in tqdm(hypotheticals, desc="[4] 写入临时 FASTA   ",
                       unit="seqs", ncols=100,
                       bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]'):
            f.write(f">{hp['protein_header']}\n{hp['protein_sequence']}\n")

    LOG.info("[4] diamond blastp: %d 条序列 → db=%s (threads=%d)",
             len(hypotheticals), db_path, threads)

    # diamond 平衡参数: 不加 --sensitive (太慢), 不加 --fast (可能漏)
    # 用 --block-size/--index-chunks 优化内存
    outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'
    cmd = (
        f"{diamond_bin} blastp "
        f"--query {tmp_faa} "
        f"--db {db_path} "
        f"--out {output_file} "
        f"--outfmt {outfmt} "
        f"--threads {threads} "
        f"--evalue {evalue} "
        f"--max-target-seqs 3 "
        f"--block-size 8 "
        f"--index-chunks 1 "
        f"--tmpdir {tempfile.gettempdir()}"
    )

    start_time = time.time()
    LOG.info("[4] diamond 正在运行, 请等待... (提示: nr库较大, 预计1-3分钟)")

    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)

        # 简单计时器, diamond 批量写入结果, 输出文件不会实时增长
        # 通过检查进程是否还在运行来判断
        pbar = tqdm(desc="[4] diamond 搜索进度  ",
                    unit="s", ncols=100,
                    bar_format='{desc}: {elapsed} |{bar}| running...')

        last_elapsed = 0
        while proc.poll() is None:
            elapsed = int(time.time() - start_time)
            if elapsed > last_elapsed:
                pbar.update(elapsed - last_elapsed)
                last_elapsed = elapsed
            time.sleep(1.0)

        pbar.close()
        stdout, stderr = proc.communicate()
        elapsed = time.time() - start_time

        if proc.returncode != 0:
            LOG.error("[4] diamond 失败 (exit=%d)", proc.returncode)
            if stderr:
                LOG.error("  stderr: %s", stderr[:800])
            if stdout:
                LOG.error("  stdout: %s", stdout[:500])
        else:
            hit_lines = sum(1 for _ in open(output_file)) if os.path.exists(output_file) else 0
            if hit_lines == 0 and stderr:
                LOG.warning("  diamond 无命中结果, stderr: %s", stderr[:500])
            LOG.info("[4] diamond 完成 → %s (耗时 %.1fs, %d hits)",
                     output_file, elapsed, hit_lines)
    except Exception as e:
        LOG.error("[4] diamond 异常: %s", str(e)[:200])
    finally:
        if os.path.exists(tmp_faa):
            os.unlink(tmp_faa)


# ══════════════════════════════════════════════════════════════
# 步骤 4c: 运行 hmmscan (HMM 结构域搜索, 补刀)
# ══════════════════════════════════════════════════════════════

def _ensure_hmm_pressed(hmm_db):
    """检查 HMM 数据库是否已 hmmpress, 若未压缩则自动执行"""
    for suf in ['.h3f', '.h3i', '.h3m', '.h3p']:
        if os.path.exists(hmm_db + suf):
            return True

    LOG.info("[4c] HMM 数据库未压缩, 自动 hmmpress: %s", hmm_db)
    ret = os.system(f"hmmpress -f {hmm_db} >/dev/null 2>&1")
    if ret == 0:
        LOG.info("[4c] hmmpress 完成")
        return True
    LOG.error("[4c] hmmpress 失败 (exit=%d)", ret)
    return False


def run_hmmscan(hypotheticals, hmm_db, output_file, threads=8, evalue=1e-3,
                hmmscan_bin="hmmscan"):
    """对 hypothetical 蛋白序列运行 hmmscan 搜索 HMM 结构域"""
    if not hypotheticals:
        LOG.warning("[4c] 没有 hypothetical protein 需要 hmmscan, 跳过")
        return

    # 自动 hmmpress
    if not _ensure_hmm_pressed(hmm_db):
        return

    # 写入临时 FASTA
    fd, tmp_faa = tempfile.mkstemp(suffix='.faa', prefix='hypo_')
    os.close(fd)
    with open(tmp_faa, 'w', encoding='utf-8') as f:
        for hp in tqdm(hypotheticals, desc="[4c] 写入临时 FASTA  ",
                       unit="seqs", ncols=100,
                       bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]'):
            f.write(f">{hp['protein_header']}\n{hp['protein_sequence']}\n")

    LOG.info("[4c] hmmscan: %d 条序列 → db=%s (threads=%d)",
             len(hypotheticals), hmm_db, threads)

    # hmmscan --tblout 输出格式: target acc query acc evalue score
    fd_tbl, tmp_tbl = tempfile.mkstemp(suffix='.tblout', prefix='hmmer_')
    os.close(fd_tbl)

    cmd = (
        f"{hmmscan_bin} "
        f"--tblout {tmp_tbl} "
        f"--cpu {threads} "
        f"--noali "
        f"-E {evalue} "
        f"{hmm_db} "
        f"{tmp_faa} "
        f">/dev/null"
    )

    start_time = time.time()
    pbar = tqdm(total=len(hypotheticals), desc="[4c] hmmscan 搜索进度  ",
                unit="query", ncols=100,
                bar_format='{desc}: {elapsed} |{bar}| hmmscan running...')

    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)

        while proc.poll() is None:
            time.sleep(1.0)
            pbar.update(0)  # keep alive

        stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            LOG.error("[4c] hmmscan 失败 (exit=%d)", proc.returncode)
            if stderr:
                LOG.error("  stderr: %s", stderr[:500])
        else:
            # 解析 hmmscan --tblout → 转换为 diamond outfmt 6 格式
            hit_count = _convert_hmmer_tblout(tmp_tbl, output_file)
            elapsed = time.time() - start_time
            LOG.info("[4c] hmmscan 完成 → %s (耗时 %.0fs, %d hits)",
                     output_file, elapsed, hit_count)
    except Exception as e:
        LOG.error("[4c] hmmscan 异常: %s", str(e)[:200])
    finally:
        pbar.close()
        if os.path.exists(tmp_faa):
            os.unlink(tmp_faa)
        if os.path.exists(tmp_tbl):
            os.unlink(tmp_tbl)


def _convert_hmmer_tblout(tblout_path, output_path):
    """将 hmmscan --tblout 转换为 diamond outfmt 6 格式

    hmmscan tblout 列:
      target_name  acc  query_name  acc  evalue  score  bias  ...
    """
    hit_count = 0
    with open(tblout_path, 'r') as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            # target=parts[0], query=parts[2], evalue=parts[4], score=parts[5]
            target = parts[0]
            target_acc = parts[1] if parts[1] != '-' else target
            query = parts[2]
            evalue = parts[4]
            score = parts[5]

            # 转为 diamond outfmt 6: qseqid sseqid pident ... evalue bitscore stitle
            # HMM 没有 pident, 用 score 替代; 没有坐标信息填 0
            fout.write('\t'.join([
                query, target_acc, '0.0', '0', '0', '0',
                '0', '0', '0', '0',
                evalue, score, target
            ]) + '\n')
            hit_count += 1
    return hit_count


# ══════════════════════════════════════════════════════════════
# 步骤 4d: 运行本地 blastp
# ══════════════════════════════════════════════════════════════

def run_blastp(hypotheticals, db_path, output_file, threads=8, evalue=1e-3):
    """对 hypothetical 蛋白序列运行 blastp"""
    if not hypotheticals:
        LOG.warning("[4] 没有 hypothetical protein 需要 BLAST, 跳过")
        return

    # 写入临时 FASTA
    fd, tmp_faa = tempfile.mkstemp(suffix='.faa', prefix='hypo_')
    os.close(fd)
    with open(tmp_faa, 'w', encoding='utf-8') as f:
        for hp in tqdm(hypotheticals, desc="[4] 写入临时 FASTA   ",
                       unit="seqs", ncols=100, bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]'):
            f.write(f">{hp['protein_header']}\n{hp['protein_sequence']}\n")

    LOG.info("[4] 运行 blastp: %d 条序列 → db=%s (threads=%d)",
             len(hypotheticals), db_path, threads)

    outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'
    cmd = (
        f"blastp "
        f"-query {tmp_faa} "
        f"-db {db_path} "
        f"-out {output_file} "
        f"-outfmt '{outfmt}' "
        f"-num_threads {threads} "
        f"-evalue {evalue} "
        f"-max_target_seqs 3"
    )

    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)
        pbar = tqdm(total=len(hypotheticals), desc="[4] blastp 搜索进度  ",
                    unit="query", ncols=100, bar_format='{desc}: {elapsed} |{bar}| blastp running...')
        start_time = time.time()

        while proc.poll() is None:
            elapsed = time.time() - start_time
            pbar.set_postfix_str(f"elapsed={elapsed:.0f}s")
            time.sleep(0.5)

        pbar.close()
        stdout, stderr = proc.communicate()

        if proc.returncode != 0:
            LOG.error("[4] blastp 失败 (exit=%d)", proc.returncode)
            if stderr:
                LOG.error("  stderr: %s", stderr[:500])
        else:
            elapsed = time.time() - start_time
            LOG.info("[4] blastp 完成 → %s (耗时 %.0fs)", output_file, elapsed)
    except subprocess.CalledProcessError as e:
        LOG.error("[4] blastp 失败 (exit=%d)", e.returncode)
    finally:
        if os.path.exists(tmp_faa):
            os.unlink(tmp_faa)


# ══════════════════════════════════════════════════════════════
# 工具: 合并 diamond + hmmscan 结果
# ══════════════════════════════════════════════════════════════

def _merge_blast_files(diamond_file, hmmer_file, merged_file):
    """合并 diamond 和 hmmscan 结果: diamond 优先, hmmscan 补刀

    策略:
      - diamond 有命中 → 保留 (序列同源优先)
      - diamond 无命中 → 用 hmmscan 的结构域命中
      - 同一 query 两个都有 → diamond 优先
    """
    diamond_queries = set()
    with open(diamond_file, 'r') as f:
        for line in f:
            if line.strip():
                diamond_queries.add(line.split('\t')[0])

    with open(merged_file, 'w') as fout:
        # 先写 diamond 结果
        with open(diamond_file, 'r') as f:
            for line in f:
                if line.strip():
                    fout.write(line)

        # hmmscan 补刀: 只写 diamond 没命中的 query
        added = 0
        with open(hmmer_file, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                query = line.split('\t')[0]
                if query not in diamond_queries:
                    fout.write(line)
                    added += 1

    LOG.info("[合并] diamond=%d 条, hmmscan 补充=%d 条 → %s",
             len(diamond_queries), added, merged_file)


# ══════════════════════════════════════════════════════════════
# 步骤 5a: 解析 blast 结果
# ══════════════════════════════════════════════════════════════

def parse_blast_results(blast_file, annot_table=None):
    """解析 blastp 输出, 返回 {(seq_id, cds_index): BlastHit}

    优先选择非 hypothetical 命中
    如果提供了 annot_table, 用真实蛋白名替换序列 ID
    """
    # 加载注释表: protein_id → protein_name
    annot_map = {}
    if annot_table and os.path.exists(annot_table):
        LOG.info("[5a] 加载注释表: %s", annot_table)
        with open(annot_table, 'r') as f:
            header = f.readline()  # skip header
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= 3:
                    annot_map[cols[0]] = cols[2]  # cols[0]=Protein_id, cols[2]=Protein_name
        LOG.info("[5a] 注释表加载: %d 条映射", len(annot_map))

    results = {}
    all_hits = {}

    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        LOG.warning("[5a] blast 结果文件为空或不存在")
        return results

    total_lines = _count_lines(blast_file)

    with open(blast_file, 'r', encoding='utf-8') as f:
        pbar = tqdm(f, total=total_lines, desc="[5a] 解析 blast 结果  ",
                    unit="hits", ncols=100, bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]')
        for line in pbar:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 12:
                continue

            query = parts[0]
            if parts[1] == '*' or parts[2] == '-1':
                continue
            key = _query_to_key(query)
            if key is None:
                continue

            # 用注释表替换 stitle
            sseqid = parts[1]
            stitle = parts[12] if len(parts) > 12 else sseqid
            if annot_map:
                real_name = annot_map.get(sseqid)
                if real_name:
                    stitle = real_name

            hit = BlastHit(
                query=query,
                sseqid=sseqid,
                pident=parts[2],
                evalue=parts[10],
                bitscore=parts[11],
                stitle=stitle,
            )

            all_hits.setdefault(key, []).append(hit)

            if key not in results:
                results[key] = hit
            elif _is_hypothetical(results[key].stitle) and not _is_hypothetical(hit.stitle):
                results[key] = hit

    # 对仍命中 "hypothetical" 的, 且 annot 中有对应名称的, 用 annot 替换
    if annot_map:
        for key, hit in results.items():
            real_name = annot_map.get(hit.sseqid)
            if real_name and _is_hypothetical(hit.stitle):
                hit.stitle = real_name

    hypo_improved = sum(
        1 for k in results
        if len(all_hits.get(k, [])) > 1
        and _is_hypothetical(all_hits[k][0].stitle)
        and not _is_hypothetical(results[k].stitle)
    )
    total_with_hits = len(results)
    LOG.info("[5a] 解析 blast 结果: %d 条命中 (%d 条跳过假想蛋白取其下一条)",
             total_with_hits, hypo_improved)
    return results


_HYPOTHETICAL_WORDS = re.compile(
    r'hypothetical|uncharacterized|unnamed|predicted|putative|unknown|DUF\d+|domain of unknown function',
    re.IGNORECASE
)


def _is_hypothetical(title):
    """判断 blast 命中是否仍是假想蛋白"""
    return bool(_HYPOTHETICAL_WORDS.search(title))


def _query_to_key(query_header):
    """从 blast query / faa 头部提取 (orig_seq_id, cds_index) 匹配键

    orig_seq_id 去除了 pyrodigal 追加的 _{gene} 后缀, 与 featuretable.tbl 中的 seq_id 对齐

    支持两种头部格式:
      - 完整格式: >seq_id_1 # start # end # strand # ID=X_Y;...
      - diamond 截断格式: seq_id_1  (diamond 截取 > 到第一个空格)
    """
    h = query_header.lstrip('>').strip()
    parts = h.split(' # ')

    # 完整格式: 从 ID= 字段提取 cds_idx
    if len(parts) >= 4:
        full_seq_id = parts[0].strip()
        cds_idx = 1
        for part in parts:
            if part.startswith('ID='):
                id_val = part.split(';')[0].replace('ID=', '')
                try:
                    _, cds_idx = id_val.split('_')
                    cds_idx = int(cds_idx)
                except ValueError:
                    pass
                break
    else:
        # diamond 截断格式: 只有 seq_id, 从末尾推断 cds_idx
        full_seq_id = h
        m = re.search(r'_(\d+)$', full_seq_id)
        cds_idx = int(m.group(1)) if m else 1

    # 去除 pyrodigal 追加的 _{cds_idx} 后缀, 恢复原始 record_id
    suffix = f"_{cds_idx}"
    if full_seq_id.endswith(suffix):
        orig_seq_id = full_seq_id[:-len(suffix)]
    else:
        orig_seq_id = re.sub(r'_\d+$', '', full_seq_id)

    return (orig_seq_id, cds_idx)


# ══════════════════════════════════════════════════════════════
# 步骤 5b: 输出更新后的 featuretable.tbl
# ══════════════════════════════════════════════════════════════

def write_updated_featuretable(records, blast_results, output_path):
    """写更新后的 featuretable.tbl"""
    updated_count = 0
    total_cds = sum(len(cdses) for cdses in records.values())

    with open(output_path, 'w', encoding='utf-8') as f:
        pbar = tqdm(total=total_cds, desc="[5b] 更新 featuretable ",
                    unit="cds", ncols=100, bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]')

        for seq_id, cdses in records.items():
            f.write(f">Feature {seq_id}\n")
            for cds in cdses:
                key = (cds.seq_id, cds.cds_index)
                blast_hit = blast_results.get(key)

                if cds.is_hypothetical and blast_hit:
                    product_name = _sanitize_product_name(blast_hit.stitle)
                    cds.product = product_name
                    inf = (f"alignment:BLASTP:2.16.0:"
                           f"evalue={blast_hit.evalue:.1e};"
                           f"pident={blast_hit.pident:.1f}%;"
                           f"{blast_hit.sseqid}")
                    cds.inferences.append(inf)
                    updated_count += 1

                for line in cds.to_tbl_lines():
                    f.write(line + '\n')
                pbar.update(1)

    LOG.info("[5b] 更新 featuretable → %s (%d 条 annotated)", output_path, updated_count)
    return updated_count


def _sanitize_product_name(raw_title, max_len=80):
    """从 blast 命中标题提取干净的功能名称"""
    # 移除 OS=, OX=, GN=, PE=, SV= 等 UniProt 后缀
    for sep in [' OS=', ' OX=', ' GN=', ' PE=', ' SV=', ' [', ' (Fragment)']:
        idx = raw_title.find(sep)
        if idx > 0:
            raw_title = raw_title[:idx]
    raw_title = raw_title.strip()
    if len(raw_title) > max_len:
        raw_title = raw_title[:max_len - 3] + '...'
    return raw_title


# ══════════════════════════════════════════════════════════════
# 步骤 5c: 输出更新后的 proteins.faa
# ══════════════════════════════════════════════════════════════

def write_updated_proteins(faa_path, blast_results, output_path):
    """写更新后的 proteins.faa (在头部追加 blast 注释)"""
    if not os.path.exists(faa_path):
        LOG.warning("proteins.faa 不存在, 跳过更新")
        return 0

    updated_count = 0
    total_lines = _count_lines(faa_path)

    with open(faa_path, 'r', encoding='utf-8') as fin, \
         open(output_path, 'w', encoding='utf-8') as fout:
        pbar = tqdm(fin, total=total_lines, desc="[5c] 更新 proteins.faa ",
                    unit="lines", ncols=100, bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}]')

        for line in pbar:
            if line.startswith('>'):
                key = _query_to_key(line.strip())
                blast_hit = blast_results.get(key) if key else None
                if blast_hit:
                    product_name = _sanitize_product_name(blast_hit.stitle)
                    line = line.rstrip('\n') + f" BLAST={product_name}\n"
                    updated_count += 1
            fout.write(line)

    LOG.info("[5c] 更新 proteins.faa → %s (%d 条 annotated)", output_path, updated_count)
    return updated_count


# ══════════════════════════════════════════════════════════════
# 工具: 打印输出摘要
# ══════════════════════════════════════════════════════════════

def print_summary(blast_results, hypotheticals):
    """打印 blast 注释摘要"""
    if not blast_results:
        return

    print(f"\n{'='*70}")
    print(f"BLAST 注释摘要 (共 {len(blast_results)} 条命中)")
    print(f"{'='*70}")
    print(f"{'Contig':<50} {'CDS':>4} {'Best Hit':<40} {'E-value':>10}")
    print(f"{'-'*50} {'-'*4} {'-'*40} {'-'*10}")

    for hp in hypotheticals:
        key = (hp['seq_id'], hp['cds_index'])
        bh = blast_results.get(key)
        if bh:
            contig = hp['seq_id'][:47] + '...' if len(hp['seq_id']) > 50 else hp['seq_id']
            hit = _sanitize_product_name(bh.stitle)[:37] + '...' if len(bh.stitle) > 40 else bh.stitle[:40]
            print(f"{contig:<50} {hp['cds_index']:>4} {hit:<40} {bh.evalue:>10.1e}")

    print(f"{'='*70}\n")


# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="analyze_hypothetical.py — 假定蛋白 blastp 功能注释工具 v2.1",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 本地 blastp (快速, 推荐批量)
  python analyze_hypothetical.py \\
      -t featuretable.tbl -f proteins.faa \\
      --blast my_blast_results.txt \\
      -d /db/nr \\
      -o ./output/ \\
      --threads 16

  # NCBI 在线 BLAST (无需本地数据库)
  python analyze_hypothetical.py \\
      -t featuretable.tbl -f proteins.faa \\
      --blast my_blast_results.txt \\
      --online \\
      -o ./output/
"""
    )
    parser.add_argument('-t', '--tbl', required=True,
                        help='featuretable.tbl 路径')
    parser.add_argument('-f', '--faa', required=True,
                        help='proteins.faa 路径')
    parser.add_argument('--blast', required=True,
                        help='blastp 结果输出文件路径')
    parser.add_argument('-d', '--db',
                        help='blastp 数据库路径 (如 /db/nr 或 /db/uniprot)')
    parser.add_argument('-o', '--output', default='./hypothetical_output',
                        help='输出目录 (默认: ./hypothetical_output)')
    parser.add_argument('-tbl-out', '--tbl-output',
                        help='更新后的 featuretable 输出路径 (默认: {output}/featuretable_updated.tbl)')
    parser.add_argument('-faa-out', '--faa-output',
                        help='更新后的 proteins.faa 输出路径 (默认: {output}/proteins_updated.faa)')
    parser.add_argument('--threads', type=int, default=8,
                        help='blastp 线程数 (默认: 8, 仅本地模式)')
    parser.add_argument('--evalue', type=float, default=1e-3,
                        help='blastp E-value 阈值 (默认: 1e-3)')
    parser.add_argument('--online', action='store_true',
                        help='使用 NCBI 在线 BLAST (无需本地数据库, 序列多时较慢)')
    parser.add_argument('--online-db', default='nr',
                        help='NCBI 在线 BLAST 数据库 (默认: nr, 可选: uniprotkb, swissprot, pdb)')
    parser.add_argument('--delay', type=int, default=3,
                        help='在线 BLAST 每次提交间隔秒数 (默认: 3; 有 API key 可降至 1)')
    parser.add_argument('--email',
                        help='NCBI 邮箱 (强烈推荐, NCBI 可能限制无邮箱的请求)')
    parser.add_argument('--ncbi-api-key',
                        help='NCBI API key (提升请求频率上限, https://ncbi.nlm.nih.gov/account/)')
    parser.add_argument('--diamond', action='store_true',
                        help='使用 diamond blastp 本地高速比对 (配合 -d 指定 .dmnd 数据库)')
    parser.add_argument('--diamond-bin', default='diamond',
                        help='diamond 二进制路径 (默认: diamond, 需在 PATH 中)')
    parser.add_argument('--hmmer', action='store_true',
                        help='使用 hmmscan 搜索 HMM 结构域 (配合 --hmmer-db 指定 .hmm 数据库)')
    parser.add_argument('--hmmer-db',
                        help='hmmscan HMM 数据库路径 (如 U-RVDBv31.0-prot.hmm)')
    parser.add_argument('--hmmer-bin', default='hmmscan',
                        help='hmmscan 二进制路径 (默认: hmmscan)')
    parser.add_argument('--annot',
                        help='功能注释表 (如 U-RVDBv31.0-prot.info.tab), 将序列ID映射为真实蛋白名')
    parser.add_argument('--skip-blast', action='store_true',
                        help='跳过 blastp 运行 (直接使用已有的 --blast 文件)')
    args = parser.parse_args()

    # 输出目录
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    tbl_out = args.tbl_output or str(out_dir / 'featuretable_updated.tbl')
    faa_out = args.faa_output or str(out_dir / 'proteins_updated.faa')
    blast_file = args.blast

    LOG.info("=" * 60)
    LOG.info("analyze_hypothetical.py v2.1")
    LOG.info("  featuretable: %s", args.tbl)
    LOG.info("  proteins.faa: %s", args.faa)
    LOG.info("  output dir:   %s", out_dir)
    LOG.info("=" * 60)

    # ── 步骤 1: 解析 featuretable.tbl ──
    records = parse_featuretable(args.tbl)

    # ── 步骤 2: 解析 proteins.faa ──
    proteins = parse_proteins_faa(args.faa)

    # ── 步骤 3: 匹配 hypothetical protein ──
    hypotheticals = match_hypotheticals(records, proteins)

    if not hypotheticals:
        LOG.warning("未找到任何 hypothetical protein, 无需 BLAST")
        return

    # ── 步骤 4: 运行 blast ──
    if not args.skip_blast:
        if args.online:
            run_online_blast(hypotheticals, blast_file, args.online_db,
                              args.delay, args.email, args.ncbi_api_key)
        elif args.diamond and args.hmmer:
            # 双模式: diamond + hmmscan 互补
            if not args.db or not args.hmmer_db:
                LOG.error("双模式需要 -d (diamond) 和 --hmmer-db 同时指定")
                sys.exit(1)
            # 用临时文件避免与 --blast 路径冲突
            diamond_tmp = str(out_dir / '.diamond_tmp.txt')
            hmmer_tmp = str(out_dir / '.hmmer_tmp.txt')
            run_diamond(hypotheticals, args.db, diamond_tmp, args.threads,
                        args.evalue, args.diamond_bin)
            run_hmmscan(hypotheticals, args.hmmer_db, hmmer_tmp, args.threads,
                        args.evalue, args.hmmer_bin)
            # 合并: diamond 优先, hmmer 补刀 → 写入用户指定的 blast_file
            _merge_blast_files(diamond_tmp, hmmer_tmp, blast_file)
            for tmp in [diamond_tmp, hmmer_tmp]:
                if os.path.exists(tmp): os.unlink(tmp)
        elif args.diamond and args.db:
            run_diamond(hypotheticals, args.db, blast_file, args.threads,
                        args.evalue, args.diamond_bin)
        elif args.hmmer and args.hmmer_db:
            run_hmmscan(hypotheticals, args.hmmer_db, blast_file, args.threads,
                        args.evalue, args.hmmer_bin)
        elif args.db:
            run_blastp(hypotheticals, args.db, blast_file, args.threads, args.evalue)
        else:
            LOG.error("需要指定 BLAST 方式: --online / --diamond -d / --hmmer / -d")
            sys.exit(1)
    else:
        LOG.info("[4] 跳过 blastp (--skip-blast), 使用已有结果: %s", blast_file)

    # ── 步骤 5a: 解析 blast 结果 ──
    blast_results = parse_blast_results(blast_file, args.annot)

    # ── 步骤 5b: 输出更新后的 featuretable ──
    n_tbl = write_updated_featuretable(records, blast_results, tbl_out)

    # ── 步骤 5c: 输出更新后的 proteins.faa ──
    n_faa = write_updated_proteins(args.faa, blast_results, faa_out)

    # ── 打印摘要 ──
    print_summary(blast_results, hypotheticals)

    # ── 最终报告 ──
    LOG.info("=" * 60)
    LOG.info("完成!")
    LOG.info("  hypothetical 总数: %d", len(hypotheticals))
    LOG.info("  blast 命中:       %d", len(blast_results))
    LOG.info("  featuretable:     %s (%d updated)", tbl_out, n_tbl)
    LOG.info("  proteins.faa:     %s (%d updated)", faa_out, n_faa)
    if len(hypotheticals) - len(blast_results) > 0:
        LOG.info("  ⚠️  %d 条 hypothetical 无 blast 命中",
                 len(hypotheticals) - len(blast_results))
    LOG.info("=" * 60)


if __name__ == '__main__':
    main()
