"""勘明: 06/07 stage 原始产物中哪些文件含 662 条待清理记录。"""
import csv
import os
import subprocess
from collections import Counter

KILL = set(l.strip() for l in open("/tmp/blacklist_kill_ids.txt") if l.strip())
print(f"662 清单: {len(KILL)} 条\n")

TARGETS = [
    ("06_HostPrediction/ensemble_host_summary.tsv", "tsv"),
    ("06_HostPrediction/host_classified_fasta/Plant.classified.fasta", "fasta"),
    ("07_Checkv/Plant.fasta", "fasta"),
]

for rel, kind in TARGETS:
    p = f"/home/zhangwenda/data-test/out/{rel}"
    if not os.path.exists(p):
        print(f"[缺失] {rel}")
        continue
    if kind == "fasta":
        n = tot = 0
        with open(p, errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    tot += 1
                    if line[1:].split()[0] in KILL:
                        n += 1
        print(f"[fasta] {rel}")
        print(f"    总 {tot} 条, 命中 662 清单 {n} 条")
    else:
        # TSV: 找 contig_id 列
        with open(p, errors="replace") as f:
            rd = csv.reader(f, delimiter="\t")
            hdr = next(rd)
            rows = list(rd)
        # 找哪列是 id
        id_cols = [i for i, h in enumerate(hdr)
                   if "contig" in h.lower() or h.lower() in ("id", "seq_id", "qseqid")]
        print(f"[tsv] {rel}")
        print(f"    总 {len(rows)} 行, 表头前 6: {hdr[:6]}")
        for c in id_cols:
            hits = sum(1 for r in rows if len(r) > c and r[c].strip() in KILL)
            print(f"    ID列 '{hdr[c]}'(#{c}): 命中 {hits}")

# 07_Checkv/Plant 目录内容
pd = "/home/zhangwenda/data-test/out/07_Checkv/Plant"
print(f"\n[dir] 07_Checkv/Plant/ 内容:")
if os.path.isdir(pd):
    for fn in sorted(os.listdir(pd))[:15]:
        fp = os.path.join(pd, fn)
        print(f"    {fn}  ({os.path.getsize(fp)} bytes)")
