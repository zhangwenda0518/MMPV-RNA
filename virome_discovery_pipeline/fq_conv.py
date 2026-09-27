#!/usr/bin/env python3
"""fq_conv.py — FASTA(任意折行) -> FASTQ(合成I质量), 用于 VSI 烟雾测试数据制备。
可选第三个参数: 最多转换的记录数 (R1/R2 分别用相同值截断以保持配对)"""
import sys

def convert(src, dst, max_records=None):
    out = []
    with open(src) as f:
        h, s = None, []
        for line in f:
            if max_records is not None and len(out) // 4 >= max_records:
                break
            line = line.strip()
            if line.startswith(">"):
                if h is not None:
                    emit(h, s, out)
                h, s = line[1:], []
            elif line:
                s.append(line)
        if h is not None:
            emit(h, s, out)
    with open(dst, "w") as w:
        w.write("".join(out))
    print(dst, len(out) // 4, "records")

def emit(h, s, out):
    seq = "".join(s)
    out.append("@%s\n%s\n+\n%s\n" % (h, seq, "I" * len(seq)))

convert(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else None)

