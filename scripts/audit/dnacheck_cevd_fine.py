#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CEVd 阳性位点精查（按位点独立统计）：
1) {sp}_viroid_word11.tsv 中 pid>=90 且 aln>=80 的命中，±200nt 合并成位点；
2) samtools faidx 抽位点 ±10kb 上下文；
3) 每段上下文跑 word 7 + dust no blastn；
4) 每位点统计：CEVd 371nt 覆盖率 / 最长连续覆盖 / 最好 pid / 片段数。
产出: cevd_fine/{sp}_fine_hits.tsv, fine_summary.tsv, {sp}_loci_context.fa"""
import os, subprocess
from collections import defaultdict

BASE = "/home/zhangwenda/eve_screen_goji3_20260929"
OUT = f"{BASE}/cevd_fine"
VIROIDS = "/home/zhangwenda/database/virus-db/viroids-db/viroids.fasta"
GEN = {
    "zhonghua": "/home/zhangwenda/Lycium_EVE/00_edta/genomes/zhonghua.fa",
    "ningxia": "/home/zhangwenda/Lycium_EVE/00_edta/genomes/ningxia.clean.fa",
    "heiguo": "/home/zhangwenda/Lycium_EVE/00_edta/genomes/heiguo.clean.fa",
}
FLANK, MERGE_GAP, CEVD_L = 10000, 200, 371

os.makedirs(OUT, exist_ok=True)
summary = open(f"{OUT}/fine_summary.tsv", "w")
summary.write("sp\tlocus_id\tlocus\tcevd_cov_pct\tbest_pid\tmax_run_nt\tn_frag\tcontext_len\n")

for sp, g in GEN.items():
    hits = defaultdict(list)
    for line in open(f"{BASE}/{sp}_viroid_word11.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 13:
            continue
        if float(p[4]) >= 90 and int(p[5]) >= 80:
            s, e = sorted((int(p[8]), int(p[9])))
            hits[p[2]].append((s, e))
    merged = []
    for c in sorted(hits):
        iv = sorted(hits[c])
        cs, ce = iv[0]
        for s, e in iv[1:]:
            if s - ce <= MERGE_GAP:
                ce = max(ce, e)
            else:
                merged.append((c, cs, ce))
                cs, ce = s, e
        merged.append((c, cs, ce))
    fa_path, regions = f"{OUT}/{sp}_loci_context.fa", f"{OUT}/{sp}_regions.txt"
    meta = {}
    with open(fa_path, "w") as fo, open(regions, "w") as fr:
        for i, (c, s, e) in enumerate(sorted(merged), 1):
            lid = f"{sp}_L{i:02d}_{c}:{s}-{e}"
            rs, re_ = max(1, s - FLANK), e + FLANK
            meta[lid] = (c, s, e)
            meta[f"{c}:{rs}-{re_}"] = lid  # region 串 -> lid
            fr.write(f"{c}:{rs}-{re_}\n")
            r = subprocess.run(["samtools", "faidx", g, "--fai-idx", f"{g}.fai",
                                f"{c}:{rs}-{re_}"], check=True, capture_output=True, text=True)
            # samtools 头行为 ">chr:s-e"，追加为统一命名
            lines = r.stdout.splitlines()
            fo.write(">" + lid + "\n")
            fo.write("\n".join(lines[1:]) + "\n")
    faidx, name, seq = {}, None, []
    for line in open(fa_path):
        if line.startswith(">"):
            if name:
                faidx[name] = "".join(seq)
            name, seq = line[1:].strip().split()[0], []
        else:
            seq.append(line.strip())
    if name:
        faidx[name] = "".join(seq)
    fout = open(f"{OUT}/{sp}_fine_hits.tsv", "w")
    fout.write("lid\tqseqid\tpident\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n")
    per = defaultdict(list)
    for lid, sub_seq in faidx.items():
        tmp = f"{OUT}/_tmp_subject.fa"
        with open(tmp, "w") as w:
            w.write(f">{lid}\n" + "\n".join(sub_seq[i:i+60] for i in range(0, len(sub_seq), 60)) + "\n")
        r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "7", "-dust", "no",
                            "-query", VIROIDS, "-subject", tmp, "-evalue", "1e-4",
                            "-outfmt", "6 qseqid pident length qstart qend sstart send evalue bitscore"],
                           capture_output=True, text=True)
        for line in r.stdout.splitlines():
            p = line.split("\t")
            fout.write(lid + "\t" + line + "\n")
            per[lid].append((float(p[1]), int(p[2]), int(p[3]), int(p[4])))
    fout.close()
    for lid in sorted(k for k, v in meta.items() if isinstance(v, tuple)):
        c, s, e = meta[lid]
        frags = per.get(lid, [])
        spans = set()
        best_pid = max((f[0] for f in frags), default=0)
        for pid, ln, qs, qe in frags:
            for x in range(min(qs, qe), max(qs, qe) + 1):
                spans.add((x - 1) % CEVD_L)
        run = mx = 0
        prev = None
        for x in sorted(spans):
            run = run + 1 if (prev is not None and x == prev + 1) else 1
            mx = max(mx, run)
            prev = x
        summary.write(f"{sp}\t{lid}\t{c}:{s}-{e}\t{len(spans)/CEVD_L*100:.0f}\t"
                      f"{best_pid:.1f}\t{mx}\t{len(frags)}\t{len(faidx.get(lid, ''))}\n")
    print(f"[{sp}] 位点 {len(meta)} 个, word7 精查完成")
summary.close()
print("done ->", OUT)
