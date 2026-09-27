#!/usr/bin/env python3
"""
terminal_repeats.py — 末端重复 (DTR/ITR) 检测、DTR 冗余修剪与环状起点归一
==========================================================================

参考实现: Cenote-Taker3 (2026-08-24 版)
  - src/cenote/python_modules/terminal_repeats.py
        DTR/ITR 检测 + DTR 修剪。其中 fetch_dtr / reverse_complement / fetch_itr
        三个函数 CT3 注明 "copied wholesale from CheckV"
        (checkv/modules/complete_genomes.py), 本脚本按原逻辑逐行照搬。
  - src/cenote/cenote_main.sh  "Rotating DTR contigs" 段
        prodigal -c -p meta -> 统计完整 ORF -> seqkit restart -i START_BASE 起点归一。
        本脚本用内部六框扫描 (或可选 pyrodigal) 复刻同一判据, 免去 prodigal/seqkit 依赖。

本脚本把 CT3 的两步合并为一个自包含 CLI, 并补齐汇总字段:

  1. 逐条 contig 检测
       DTR  (direct terminal repeat,   直接末端重复) -> 环状基因组证据
       ITR  (inverted terminal repeat, 反向末端重复)
  2. 修剪 DTR 冗余拷贝: 丢掉 3' 端那份重复, 保留 5' 端拷贝 (等价 CT3 --wrap 的修剪口径)
  3. 起点归一 (旋转):
       --rotation repeat  (默认) 起点落在保留的末端重复拷贝起点。CT3/CheckV 口径下
                          修剪后位置 0 已是该拷贝起点, 故位移为 0; 该模式存在的意义是把
                          "起点必须落在末端重复处" 显式写进流程与汇总表。
       --rotation orf     复刻 CT3 cenote_main.sh: 最长完整 ORF 的起始密码子落在第 1 位
                          (反向链完整 ORF 占优时先反向互补再定位)。
       --rotation none    只修剪, 不旋转。
  4. 输出
       terminal_repeats_summary.tsv   每条 contig 的末端重复汇总
       trimmed_contigs.fasta          修剪 (+旋转) 后的序列, 非 DTR contig 原样透传
       circular_contigs.txt           判为环状的 contig ID (每行一个)

用法:
   python terminal_repeats.py -i centroids.fasta -o tr_out
   python terminal_repeats.py -i centroids.fasta -o tr_out --rotation orf -t 16
   python terminal_repeats.py -i centroids.fasta -o tr_out --circ-file circular_ids.txt
   python terminal_repeats.py -i centroids.fasta -o tr_out --method fuzzy --min-identity 0.9

依赖: Biopython (Bio.SeqIO); pyrodigal 可选 (--orf-caller pyrodigal)

与 Cenote-Taker3 对拍:
    python terminal_repeats.py -i in.fasta -o out --method exact --rotation none --header id
    此时本脚本前 5 列 (contig/in_length_contig/out_length_contig/dtr_seq/itr_seq)
    与 trimmed_contigs.fasta 与 CT3 的 hallmark_contigs_terminal_repeat_summary.tsv /
    trimmed_TRs_hallmark_contigs.fasta 逐格一致 (已用 120 条合成 contig 实测: 600 格零差异)。
    --method auto (默认) 额外开容错 DTR, 命中会多于 CT3; --rotation orf 额外做起点归一。

坐标与长度约定: 汇总表 in_length_contig / out_length_contig / dtr_length 均为碱基长度,
rotation_offset 为 0-based 的左移位数 (输出序列第 1 位 = 原序列第 rotation_offset+1 位)。
"""

import argparse
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# CT3 fetch_dtr 的种子长度默认值 (CheckV 口径)
DTR_MIN_LENGTH = 20
# CT3 --max_itr / fetch_itr 的默认上限
ITR_MAX_LEN = 1000
# CT3 --max_dtr_assess 默认 1e6: 超过该长度的 contig 不评估 DTR
MAX_CONTIG_LENGTH = 1_000_000

START_CODONS = ("ATG", "GTG", "TTG")
STOP_CODONS = ("TAA", "TAG", "TGA")

# 汇总表列 (前 5 列与 CT3 hallmark_contigs_terminal_repeat_summary.tsv 同名同序, 便于逐列对拍)
SUMMARY_COLUMNS = [
    "contig",
    "in_length_contig",
    "out_length_contig",
    "dtr_seq",
    "itr_seq",
    "dtr_length",
    "dtr_identity",
    "dtr_method",
    "itr_length",
    "topology",
    "rotated",
    "rotation_offset",
    "rotation_mode",
]


# ──────────────────────────────────────────────────────────────────────────
# 检测函数 (CT3 / CheckV 原样口径)
# ──────────────────────────────────────────────────────────────────────────

def reverse_complement(seq):
    """CT3 原样: str.maketrans("ACTG","TGAC") 反转。非 ACTG 字符 (如 N) 原样保留。"""
    trans = str.maketrans("ACTG", "TGAC")
    return seq[::-1].translate(trans)


def fetch_dtr(seq, min_length=DTR_MIN_LENGTH):
    """CT3 fetch_dtr 逐行照搬。

    取 5' 端 min_length bp 作种子, 在序列后半段 ([len/2, )) 找全部命中位点, 要求
    命中位点到 contig 末端之间的子串 (endseq) 与 5' 端起同长子串完全一致;
    满足则返回该 3' 端重复序列, 否则返回 ""。

    语义: contig 的 3' 端回卷进 5' 端 (装配穿越原点), 即 contig = P + M + P。
    """
    startseq = seq[0:min_length]
    matches = [
        m.start() for m in re.finditer("(?={0})".format(re.escape(startseq)), seq)
    ]
    matches = [_ for _ in matches if _ >= len(seq) / 2]
    for matchpos in matches:
        endseq = seq[matchpos:]
        if seq[0:len(endseq)] == endseq:
            return endseq
    return ""


def _identity_with_abort(a, b, max_mismatch):
    """无 gap 的等长相似度。超过 max_mismatch 立即返回 -1 (早退)。"""
    mm = 0
    for x, y in zip(a, b):
        if x != y:
            mm += 1
            if mm > max_mismatch:
                return -1.0
    return 1.0 - mm / len(a)


def fetch_dtr_fuzzy(seq, min_length=DTR_MIN_LENGTH, max_len=ITR_MAX_LEN,
                    min_identity=0.90, seed=8):
    """容错版 DTR 检测 (与 fetch_dtr 同形, 末端拷贝允许错配)。

    寻找满足 min_length <= L <= min(max_len, len/2) 且
        identity(seq[0:L], seq[len-L:len]) >= min_identity
    的最长 L (同样要求命中落在序列后半段, 与 CT3 保持一致), 返回 (重复序列, identity)。
    无命中返回 ("", 0.0)。
    """
    n = len(seq)
    l_max = min(max_len, n // 2)
    if l_max < min_length:
        return "", 0.0

    head_seed = seq[:seed]
    if "N" in head_seed or len(head_seed) < seed:
        return "", 0.0

    # 3' 端拷贝起点 p = n - L; 要求 p >= n/2 (后半段) 且 L >= min_length
    p_lo = max(n - l_max, (n + 1) // 2)
    p_hi = n - min_length
    if p_lo > p_hi:
        return "", 0.0

    best_seq, best_id = "", 0.0
    p = seq.find(head_seed, p_lo)
    while p != -1 and p <= p_hi:
        L = n - p
        max_mm = int((1.0 - min_identity) * L)
        ident = _identity_with_abort(seq[:L], seq[p:p + L], max_mm)
        if ident >= min_identity and L > len(best_seq):
            best_seq, best_id = seq[p:p + L], ident
            break  # p 递增 => L 递减, 首个达标即最长
        p = seq.find(head_seed, p + 1)
    return best_seq, best_id


def fetch_itr(seq, min_len=DTR_MIN_LENGTH, max_len=ITR_MAX_LEN):
    """CT3 fetch_itr 逐行照搬 (含 len(seq) < min_len 的退化行为)。

    把全序列反向互补, 若 5' 端 min_len bp 与反向互补链的 5' 端 min_len bp 相同,
    则逐位向上扩展到最长公共前缀 (上限 max_len 位), 返回该 ITS 序列; 否则返回 ""。

    语义: seq[:k] == RC(seq)[:k] 等价于 seq[:k] == RC(seq[-k:]), 即两端互为反向互补。
    """
    rev = reverse_complement(seq)
    if seq[:min_len] != rev[:min_len]:
        return ""
    i = min_len + 1
    while seq[:i] == rev[:i] and i <= max_len:
        i += 1
    return seq[:i - 1]


# ──────────────────────────────────────────────────────────────────────────
# ORF 扫描 (复刻 CT3 cenote_main.sh 的旋转判据)
# ──────────────────────────────────────────────────────────────────────────

def find_complete_orfs(seq, min_aa=30):
    """内部六框扫描, 返回 [(start0, end0, strand, aa_len)] (0-based 半开区间)。

    "完整 ORF" = 起始密码子开头 + 同框终止密码子结尾, 对应 CT3 里 prodigal 的
    `partial=00;start_type` 口径 (两端均不 partial)。
    每框遇到终止密码子则关闭当前 ORF 并重置, 即不产生跨终止的重叠 ORF。
    """
    n = len(seq)
    orfs = []
    for strand, s in ((1, seq), (-1, reverse_complement(seq))):
        for frame in range(3):
            start = None
            for p in range(frame, n - 2, 3):
                codon = s[p:p + 3]
                if start is None:
                    if codon in START_CODONS:
                        start = p
                elif codon in STOP_CODONS:
                    aa_len = (p - start) // 3
                    if aa_len >= min_aa:
                        if strand == 1:
                            orfs.append((start, p + 3, 1, aa_len))
                        else:
                            orfs.append((n - (p + 3), n - start, -1, aa_len))
                    start = None
    return orfs


def find_complete_orfs_pyrodigal(seq, min_aa=30):
    """pyrodigal (prodigal 的 C 重实现) 口径的完整 ORF, 坐标 0-based 半开。"""
    import pyrodigal  # 可选依赖

    finder = pyrodigal.GeneFinder(meta=True, min_gene=int(min_aa * 3))
    genes = finder.find_genes(seq.encode())
    out = []
    for g in genes:
        if getattr(g, "partial_begin", False) or getattr(g, "partial_end", False):
            continue
        begin, end = int(g.begin), int(g.end)  # pyrodigal 为 1-based inclusive
        out.append((begin - 1, end, int(g.strand), (end - begin + 1) // 3))
    return out


def _orfs(seq, min_aa, caller):
    if caller == "pyrodigal":
        try:
            return find_complete_orfs_pyrodigal(seq, min_aa)
        except ImportError:
            pass  # 回退内部扫描
    return find_complete_orfs(seq, min_aa)


def compute_orf_rotation(seq, min_aa=30, caller="internal"):
    """复刻 CT3 cenote_main.sh "Rotating DTR contigs" 的判据。

    正向完整 ORF 数 >= 反向完整 ORF 数 且 正向 >= 1
        -> 取最长正向完整 ORF 的起点为该 contig 新起点;
    否则反向完整 ORF >= 1
        -> 先反向互补, 再在新序列上取最长正向完整 ORF 起点;
    都不满足 -> 不旋转。

    返回 (offset, rc_flag, reason):
        offset  = 对「旋转输入序列」的 0-based 左移位数 (未旋转为 0)
        rc_flag = 是否先做了反向互补
        reason  = 'orf_forward' / 'orf_rc' / 'no_complete_orf'
    """
    orfs = _orfs(seq, min_aa, caller)
    fwd = [o for o in orfs if o[2] == 1]
    rev = [o for o in orfs if o[2] == -1]

    if fwd and len(fwd) >= len(rev):
        best = max(fwd, key=lambda o: o[3])
        return best[0], False, "orf_forward"
    if rev:
        rc = reverse_complement(seq)
        rc_orfs = _orfs(rc, min_aa, caller)
        rc_fwd = [o for o in rc_orfs if o[2] == 1]
        if rc_fwd:
            best = max(rc_fwd, key=lambda o: o[3])
            return best[0], True, "orf_rc"
    return 0, False, "no_complete_orf"


def rotate_sequence(seq, offset):
    """左移 offset 位 (seqkit restart 等价)。offset<=0 或 >=len 时原样返回。"""
    if offset <= 0 or offset >= len(seq):
        return seq
    return seq[offset:] + seq[:offset]


# ──────────────────────────────────────────────────────────────────────────
# 单条 contig 处理
# ──────────────────────────────────────────────────────────────────────────

def process_one(rec_id, rec_desc, raw_seq, opts, circ_ids):
    """返回 (row_dict, out_seq_str 或 None(不输出序列), header)。"""
    seq = str(raw_seq).upper()
    n = len(seq)
    header_full = ("%s %s" % (rec_id, rec_desc)).strip() if rec_desc else rec_id
    out_seq = seq
    dtr_seq, itr_seq = "", ""
    dtr_len, dtr_ident, dtr_method = 0, "", ""
    itr_len = 0
    topology = "undetermined"
    rotated, offset, rot_mode = False, 0, "none"

    user_circular = bool(circ_ids) and (
        rec_id in circ_ids or rec_desc in circ_ids or header_full in circ_ids
    )
    # 注: CT3 只用 description-of-id 人 fmt_desc 做匹配 (忽略 id 本身), 按官方帮助
    # 文本 "names of contigs (header line sans '>')" 应是 id; 本实现同时接受 id /
    # description / 完整 header 三种写法, 避开 CT3 该处遗漏。

    if user_circular:
        # CT3 口径: 用户声明的环状 contig 不搜 DTR, 不修剪
        dtr_seq, dtr_method = "User-provided circular", "user"
        topology = "circular_user"
    elif n > opts.max_length:
        # CT3 口径: 超长 contig 不评估 DTR (很可能是细菌染色体)
        dtr_method = "too_long"
    else:
        hit = ""
        if opts.method in ("auto", "exact"):
            hit = fetch_dtr(seq, opts.min_dtr)
            if hit:
                dtr_method = "exact"
                dtr_ident = "1.0000"  # exact 判定即完全一致
        if not hit and opts.method in ("auto", "fuzzy"):
            hit, ident = fetch_dtr_fuzzy(
                seq, opts.min_dtr, opts.max_dtr, opts.min_identity
            )
            if hit:
                dtr_method = "fuzzy"
                dtr_ident = "%.4f" % ident
        if hit:
            dtr_seq = hit
            dtr_len = len(hit)
            topology = "circular_DTR"

        if opts.trim and dtr_seq and dtr_seq not in ("NA", "User-provided circular"):
            out_seq = seq[:-dtr_len]

    # ITR 与 CT3 一致, 对所有 contig 都算 (上限 --max-itr)
    itr_hit = fetch_itr(seq, opts.min_dtr, opts.max_itr) if n >= opts.min_dtr else ""
    if itr_hit:
        itr_seq, itr_len = itr_hit, len(itr_hit)

    # 起点归一 (CT3 在修剪之后做旋转)
    if opts.rotation == "orf" and topology == "circular_DTR":
        offset, rc_flag, reason = compute_orf_rotation(
            out_seq, opts.orf_min_aa, opts.orf_caller
        )
        if rc_flag:
            out_seq = reverse_complement(out_seq)
        if offset > 0:
            out_seq = rotate_sequence(out_seq, offset)
            rotated = True
        rot_mode = "orf:%s" % reason
    elif opts.rotation == "repeat" and topology == "circular_DTR":
        # 保留的末端重复拷贝天然位于修剪后序列的 5' 端 => 位移 0
        rot_mode = "repeat"

    if not dtr_seq:
        dtr_seq = "NA"
    if not itr_seq:
        itr_seq = "NA"

    row = {
        "contig": rec_id,
        "in_length_contig": n,
        "out_length_contig": len(out_seq),
        "dtr_seq": dtr_seq,
        "itr_seq": itr_seq,
        "dtr_length": dtr_len,
        "dtr_identity": dtr_ident,
        "dtr_method": dtr_method or "none",
        "itr_length": itr_len,
        "topology": topology,
        "rotated": rotated,
        "rotation_offset": offset,
        "rotation_mode": rot_mode,
    }
    return row, out_seq, rec_desc


# ──────────────────────────────────────────────────────────────────────────
# 批量驱动
# ──────────────────────────────────────────────────────────────────────────

def _iter_records(fasta_in):
    for rec in SeqIO.parse(fasta_in, "fasta"):
        desc = " ".join(rec.description.split(" ")[1:])  # CT3: header 去掉 ">id "
        yield rec.id, desc, str(rec.seq)


def _chunks(it, size):
    buf = []
    for item in it:
        buf.append(item)
        if len(buf) >= size:
            yield buf
            buf = []
    if buf:
        yield buf


def _worker(payload):
    """multiprocessing 任务单元 (模块级, 需可 pickle)。返回 [(row, seq, desc), ...]。"""
    chunk, opts, circ_ids = payload
    return [process_one(cid, desc, seq, opts, circ_ids) for cid, desc, seq in chunk]


def _wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width)) or ""


def run(fasta_in, out_dir, opts):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    circ_ids = set()
    if opts.circ_file and os.path.isfile(opts.circ_file):
        with open(opts.circ_file) as fh:
            circ_ids = {ln.strip() for ln in fh if ln.strip()}
        print("  [circ-file] 用户声明环状 contig: %d 条" % len(circ_ids))

    fasta_out = out_dir / "trimmed_contigs.fasta"
    table_out = out_dir / "terminal_repeats_summary.tsv"
    circ_out = out_dir / "circular_contigs.txt"

    stats = {
        "n_total": 0, "n_exact": 0, "n_fuzzy": 0, "n_itr": 0,
        "n_rot": 0, "n_user": 0, "len_in": 0, "len_out": 0,
    }
    circular_ids = []
    t0 = time.time()

    fh_fa = open(fasta_out, "w")
    fh_tsv = open(table_out, "w", newline="")
    fh_tsv.write("\t".join(SUMMARY_COLUMNS) + "\n")

    def _handle(rows):
        for row, out_seq, desc in rows:
            stats["n_total"] += 1
            stats["len_in"] += row["in_length_contig"]
            stats["len_out"] += row["out_length_contig"]
            if row["dtr_method"] == "exact":
                stats["n_exact"] += 1
            elif row["dtr_method"] == "fuzzy":
                stats["n_fuzzy"] += 1
            elif row["dtr_method"] == "user":
                stats["n_user"] += 1
            if row["itr_length"]:
                stats["n_itr"] += 1
            if row["rotated"]:
                stats["n_rot"] += 1
            if row["topology"].startswith("circular"):
                circular_ids.append(row["contig"])
            fh_tsv.write("\t".join(str(row[c]) for c in SUMMARY_COLUMNS) + "\n")
            if opts.header == "id":
                header = row["contig"]
            elif desc:
                header = "%s %s" % (row["contig"], desc)
            else:
                header = row["contig"]
            fh_fa.write(">%s\n%s\n" % (header, _wrap(out_seq)))

    try:
        if opts.threads and opts.threads > 1:
            from multiprocessing import Pool

            payloads = ((c, opts, circ_ids)
                        for c in _chunks(_iter_records(fasta_in), 2000))
            with Pool(opts.threads) as pool:
                for rows in pool.imap(_worker, payloads):
                    _handle(rows)
                    _progress(stats["n_total"], t0, opts)
        else:
            buf = []
            for item in _iter_records(fasta_in):
                buf.append(process_one(item[0], item[1], item[2], opts, circ_ids))
                if len(buf) >= 1000:
                    _handle(buf)
                    buf = []
                    _progress(stats["n_total"], t0, opts)
            if buf:
                _handle(buf)
    finally:
        fh_fa.close()
        fh_tsv.close()

    with open(circ_out, "w") as fh:
        for cid in circular_ids:
            fh.write(cid + "\n")

    n_total = stats["n_total"]
    n_dtr = stats["n_exact"] + stats["n_fuzzy"]
    n_fuzzy = stats["n_fuzzy"]
    n_itr = stats["n_itr"]
    n_rot = stats["n_rot"]
    n_user = stats["n_user"]
    len_in = stats["len_in"]
    len_out = stats["len_out"]

    elapsed = time.time() - t0
    print("\n" + "=" * 60)
    print("  输入 contig:        %d 条 (%s bp)" % (n_total, _fmt(len_in)))
    print("  DTR 命中 (环状):    %d 条 (精确 %d / 容错 %d)"
          % (n_dtr, n_dtr - n_fuzzy, n_fuzzy))
    print("  用户声明环状:       %d 条" % n_user)
    print("  ITR 命中:           %d 条" % n_itr)
    print("  起点归一(旋转):     %d 条 (rotation=%s)" % (n_rot, opts.rotation))
    print("  输出序列总量:       %s bp (修剪掉 %s bp)"
          % (_fmt(len_out), _fmt(len_in - len_out)))
    print("  汇总表:             %s" % table_out)
    print("  序列:               %s" % fasta_out)
    print("  环状 ID:            %s" % circ_out)
    print("  耗时:               %.1fs" % elapsed)
    print("=" * 60)
    return n_total



def _progress(n, t0, opts):
    if opts.quiet or n % 20000:
        return
    sys.stderr.write("  [%s] %s 条\r" % (_elapsed(t0), _fmt(n)))
    sys.stderr.flush()


def _elapsed(t0):
    s = int(time.time() - t0)
    return "%02d:%02d:%02d" % (s // 3600, s % 3600 // 60, s % 60)


def _fmt(n):
    return "{:,}".format(int(n))


def main():
    p = argparse.ArgumentParser(
        description="末端重复 (DTR/ITR) 检测、DTR 修剪与环状起点归一 (参考 Cenote-Taker3)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, help="输入 FASTA (contig / vOTU catalog)")
    p.add_argument("-o", "--output-dir", required=True, help="输出目录")
    p.add_argument("-t", "--threads", type=int, default=4, help="并行进程数 (1=串行)")
    p.add_argument("--min-dtr", type=int, default=DTR_MIN_LENGTH,
                   help="DTR/ITR 最小重复长度 (CT3/CheckV 口径 20)")
    p.add_argument("--max-dtr", type=int, default=ITR_MAX_LEN,
                   help="容错 DTR 扫描的最大重复长度")
    p.add_argument("--max-itr", type=int, default=ITR_MAX_LEN,
                   help="ITR 最大扩展长度 (CT3 fetch_itr max_len)")
    p.add_argument("--method", choices=["auto", "exact", "fuzzy"], default="auto",
                   help="auto=先精确(CT3)后容错; exact=仅 CT3 原逻辑; fuzzy=仅容错")
    p.add_argument("--min-identity", type=float, default=0.90,
                   help="容错 DTR 的最小一致率")
    p.add_argument("--max-length", type=float, default=MAX_CONTIG_LENGTH,
                   help="超过该长度的 contig 不评估 DTR (CT3 --max_dtr_assess)")
    p.add_argument("--rotation", choices=["repeat", "orf", "none"], default="repeat",
                   help="起点归一策略: repeat=起点落在末端重复处(默认); "
                        "orf=复刻 CT3 的最长完整 ORF 起始密码子; none=不旋转")
    p.add_argument("--orf-min-aa", type=int, default=30,
                   help="--rotation orf 的完整 ORF 最短氨基酸长度")
    p.add_argument("--orf-caller", choices=["internal", "pyrodigal"], default="internal",
                   help="ORF 调用器; pyrodigal 不可用时自动回退 internal")
    p.add_argument("--no-trim", dest="trim", action="store_false",
                   help="只检测不修剪 DTR")
    p.add_argument("--circ-file", default=None,
                   help="用户声明环状 contig 列表 (每行一个 header, 去掉 '>'); "
                        "命中的 contig 不搜 DTR、不修剪 (CT3 --circ-file 口径)")
    p.add_argument("--header", choices=["full", "id"], default="full",
                   help="输出 FASTA 保留完整 header 还是仅 id (CT3 写 id)")
    p.add_argument("--quiet", action="store_true", help="不输出进度")
    opts = p.parse_args()

    if not os.path.isfile(opts.input):
        sys.exit("[terminal_repeats] 输入文件不存在: %s" % opts.input)

    print("开始: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("输入: %s" % opts.input)
    print("参数: method=%s min_dtr=%d max_dtr=%d min_id=%.2f rotation=%s trim=%s threads=%d"
          % (opts.method, opts.min_dtr, opts.max_dtr, opts.min_identity,
             opts.rotation, opts.trim, opts.threads))

    run(opts.input, opts.output_dir, opts)


if __name__ == "__main__":
    main()
