#!/usr/bin/env python3
"""从 query.fasta 抽指定 id 的序列（服务器用）"""
import sys

ids, fa = sys.argv[1], sys.argv[2]
want = set(open(ids).read().split())
name, keep, seq = None, False, []
n = 0
out = open(fa, "w")
for line in open(sys.argv[3] if len(sys.argv) > 3 else "t3_query.fasta"):
    if line.startswith(">"):
        if keep and name:
            out.write(">" + name + "\n" + "\n".join(seq[j:j+60] for j in range(0, len(seq), 60)) + "\n")
            n += 1
        name = line[1:].split()[0]
        keep = name in want
        seq = []
    else:
        seq.append(line.strip())
if keep and name:
    out.write(">" + name + "\n" + "\n".join(seq[j:j+60] for j in range(0, len(seq), 60)) + "\n")
    n += 1
out.close()
print("written", n)
