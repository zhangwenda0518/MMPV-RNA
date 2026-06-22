#!/usr/bin/env python3
"""
viral_topology.py — 病毒基因组拓扑判定 (DTR/ITR → circular/linear)
============================================================
基于 Cenote-Taker3 terminal_repeats.py 的 CheckV DTR 检测算法

输入: FASTA 文件
输出: topology_summary.tsv (contig | length | dtr | itr | topology)

用法:
  python viral_topology.py -i virus.fasta -o topology_summary.tsv
  python viral_topology.py -i virus.fasta -o topology_summary.tsv --min-dtr 20 --max-itr 1000
"""

import argparse, os, re, sys
from Bio import SeqIO
import pandas as pd


def fetch_dtr(fullseq, min_length=20):
    """CheckV DTR 检测: 序列两端的直接重复"""
    startseq = fullseq[:min_length]
    matches = [
        m.start() for m in re.finditer("(?={0})".format(re.escape(startseq)), fullseq)
    ]
    matches = [m for m in matches if m >= len(fullseq) / 2]
    for matchpos in matches:
        endseq = fullseq[matchpos:]
        if fullseq[:len(endseq)] == endseq:
            return endseq
    return ""


def reverse_complement(seq):
    trans = str.maketrans("ACTGactg", "TGACtgac")
    return seq[::-1].translate(trans)


def fetch_itr(seq, min_len=20, max_len=1000):
    """ITR 检测: 末端反向互补重复"""
    rev = reverse_complement(seq)
    if seq[:min_len] == rev[:min_len]:
        i = min_len + 1
        while seq[:i] == rev[:i] and i <= max_len:
            i += 1
        return seq[:i - 1]
    return ""


def main():
    p = argparse.ArgumentParser(description="病毒基因组拓扑判定")
    p.add_argument("-i", "--input", required=True, help="输入 FASTA")
    p.add_argument("-o", "--output", default="topology_summary.tsv", help="输出 TSV")
    p.add_argument("--min-dtr", type=int, default=20, help="DTR 最小长度 bp (默认 20)")
    p.add_argument("--max-itr", type=int, default=1000, help="ITR 最大检测长度 bp (默认 1000)")
    p.add_argument("--max-length", type=int, default=1000000, help="最大检测长度 bp (超限跳过DTR)")
    p.add_argument("--circ-file", help="用户指定的环状 contig 列表 (每行一个)")
    args = p.parse_args()

    circ_ids = set()
    if args.circ_file and os.path.isfile(args.circ_file):
        with open(args.circ_file) as cf:
            for line in cf:
                circ_ids.add(line.strip())

    results = []
    for rec in SeqIO.parse(args.input, "fasta"):
        cid = rec.id
        seq = str(rec.seq)
        seq_len = len(seq)

        if cid in circ_ids:
            dtr = "User-provided circular"
            itr = ""
        elif seq_len > args.max_length:
            dtr = "NA (too long)"
            itr = ""
        else:
            dtr = fetch_dtr(seq, args.min_dtr) or "NA"
            itr = fetch_itr(seq, args.min_dtr, args.max_itr) or "NA"

        topology = "circular" if (dtr and dtr != "NA") else "linear"
        dtr_len = len(dtr) if dtr and dtr != "NA" and "User" not in str(dtr) else 0
        itr_len = len(itr) if itr and itr != "NA" else 0

        results.append({
            "contig_id": cid,
            "length": seq_len,
            "dtr": dtr if dtr else "NA",
            "dtr_length": dtr_len,
            "itr": itr if itr else "NA",
            "itr_length": itr_len,
            "topology": topology,
        })

    df = pd.DataFrame(results)
    df.to_csv(args.output, sep="\t", index=False)

    n_circ = sum(1 for r in results if r["topology"] == "circular")
    n_lin = len(results) - n_circ
    print(f"拓扑判定: {len(results)} 序列 → 环状={n_circ} 线状={n_lin} → {args.output}")


if __name__ == "__main__":
    main()
