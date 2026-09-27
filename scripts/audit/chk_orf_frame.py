#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_orf_frame.py — 决定性对照：ACVirus 的 Prodigal ORF 是否被翻译错读框
  A) diamond blastp : ACVirus contigs.faa 里的 ORF（就是 ACVirus 实际拿去搜的蛋白）
  B) diamond blastx : 原始 contig 六读框（mmseqs/diamond_lca 实际走的路）
两者打同一个库 U-RVDBv31.0-prot_unique.dmnd
"""
import os, subprocess, sys

DS = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
CENT = os.path.join(DS, "04_CLUSTER/4_centroids/final_centroids.fasta")
FAA = os.path.join(DS, "05_Taxonomy/Votus.classed/ACVirus_results/Votus.acvirus/contigs.faa")
DB = "/home/zhangwenda/database/virus-db/RVDB-v31/RVDB_viroids.diamond_db/U-RVDBv31.0-prot_unique.dmnd"
DIA = "/home/zhangwenda/biosoft/binary/diamond"

T = [
 "CRR1440136_clean_NODE_2764", "CRR1440137_clean_NODE_812",
 "CRR2703946_clean_NODE_2830", "CRR2703948_clean_NODE_611",
 "CRR527041_clean_NODE_304", "CRR527042_clean_NODE_477",
 "CRR527046_clean_NODE_441", "CRR527046_clean_NODE_560",
 "SRR22742702_clean_NODE_373", "SRR23107136_clean_NODE_15609",
 "SRR23215179_clean_NODE_2046", "SRR24305226_clean_NODE_2281",
 "SRR33536650_clean_NODE_3026",
]


def extract(src, dst, want_orf=False):
    if not os.path.exists(src):
        print("  [缺失] %s" % src); return {}
    meta, cur, hdr = {}, None, None
    n = 0
    with open(src, "r", encoding="utf-8", errors="replace") as fi, \
         open(dst, "w", encoding="utf-8") as fo:
        for line in fi:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                hdr = line.rstrip("\n")
                if any(cur.startswith(t) for t in T):
                    fo.write(line); n += 1
                    meta[cur] = hdr
                    keep = True
                else:
                    keep = False
            elif keep:
                fo.write(line)
    print("  写入 %s ：%d 条" % (dst, n))
    return meta


print("=" * 130)
print("A. 抽取序列")
morf = extract(FAA, "/tmp/orf13.faa", True)
mcen = extract(CENT, "/tmp/c13.fasta")
for k in sorted(morf):
    parts = morf[k].split("#")
    ln = None
    print("   ORF %-52s 坐标 %s-%s 链 %s" % (k, parts[1].strip() if len(parts) > 1 else "?",
                                          parts[2].strip() if len(parts) > 2 else "?",
                                          parts[3].strip() if len(parts) > 3 else "?"))

FMT = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "evalue", "bitscore"]


def run(mode, q, out):
    cmd = [DIA, mode, "-q", q, "-d", DB, "--outfmt", "6"] + FMT + \
          ["--max-target-seqs", "4", "--quiet", "-o", out]
    print("\n   $ " + " ".join(cmd))
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print("   !! 退出码 %d：%s" % (r.returncode, r.stdout.decode()[:400]))
        return []
    rows = []
    if os.path.exists(out):
        for line in open(out):
            rows.append(line.rstrip("\n").split("\t"))
    return rows


print("\n" + "=" * 130)
print("B. diamond blastp（ACVirus 实际查询的 ORF 蛋白）")
bp = run("blastp", "/tmp/orf13.faa", "/tmp/orf13_bp.tsv")
for r in bp:
    print("   %-44s -> %-14s pid=%-6s len=%-5s q=%s-%s ev=%-11s bits=%s"
          % (r[0][:44], r[1], r[2], r[3], r[4], r[5], r[6], r[7]))

print("\n" + "=" * 130)
print("C. diamond blastx（原始 contig 六读框）")
bx = run("blastx", "/tmp/c13.fasta", "/tmp/c13_bx.tsv")
for r in bx:
    print("   %-44s -> %-14s pid=%-6s len=%-5s q=%s-%s ev=%-11s bits=%s"
          % (r[0][:44], r[1], r[2], r[3], r[4], r[5], r[6], r[7]))

print("\n" + "=" * 130)
print("D. 汇总 subject（用 eutils 标名称）")
ids = sorted({r[1] for r in bp + bx})
print("   %d 个 subject: %s" % (len(ids), ",".join(ids)))
if ids:
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id="
           + ",".join(ids) + "&retmode=json")
    try:
        import json, urllib.request
        d = json.load(urllib.request.urlopen(url, timeout=60))["result"]
        for k, v in d.items():
            if k == "uids":
                continue
            print("   %-14s | %-38s | %s" % (k, v.get("organism", "?")[:38], v.get("title", "")[:80]))
    except Exception as e:
        print("   [eutils 失败] %s" % e)
