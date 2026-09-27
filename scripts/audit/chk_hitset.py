#!/usr/bin/env python3
"""命中集到底是什么：参考基因组身份 + HG999358 被命中基因的注释"""
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

DB = Path("/home/zhangwenda/database/virus-db/acvirus_db")
W = Path("/tmp/bia_check")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def sh(c):
    p = subprocess.run(c, shell=True, capture_output=True, text=True)
    return (p.stdout or "") + (p.stderr or "")


rows = [l.split("\t") for l in open(W / "dia.tsv") if l.strip()]
by_q = defaultdict(list)
for r in rows:
    by_q[r[0]].append(r)

print("=" * 120)
print("1. 每个 contig 的 top5 命中（accession / pid / alen / bits）")
print("=" * 120)
accs = Counter()
for q, rs in by_q.items():
    rs.sort(key=lambda r: -float(r[5]))
    print("\n%s" % q)
    for r in rs[:5]:
        acc = r[1].rsplit("_", 1)[0]
        accs[acc] += 1
        print("    %-16s pid=%-6s alen=%-5s bits=%-6s e=%s" % (acc, r[2], r[3], r[5], r[4]))
print("\n命中出现的参考 accession 全集 (%d 个):" % len(accs))
for a, n in accs.most_common():
    print("    %-16s 出现在 %d 个 contig 的 top10 里" % (a, n))

print("\n" + "=" * 120)
print("2. 这些参考在 ACVirus 库里的分类归属")
print("=" * 120)
# taxa.txt / vmr.tsv 是 tab 分隔；按 accession 反查
for f in ["vr.tsv", "vmr.tsv", "taxa.txt", "processed_accessions_b.tsv", "processed_accessions_b.fa_names.tsv"]:
    p = DB / f
    if not p.exists():
        continue
    print("\n--- %s ---" % f)
    hits = 0
    for ln in open(p, encoding="utf-8", errors="replace"):
        if any(a in ln for a in accs):
            hits += 1
            print("    %s" % ln.rstrip()[:230])
            if hits >= 18:
                print("    ... (截断)")
                break
    if hits == 0:
        print("    (无命中)")

print("\n" + "=" * 120)
print("3. 未在库里反查到的 accession，去 NCBI 拿标题")
print("=" * 120)
blob = sh('curl -s -g "%s/esummary.fcgi?db=nuccore&id=%s&retmode=json"'
          % (EUTILS, ",".join(accs.keys())))
try:
    res = json.loads(blob)["result"]
    for k, v in res.items():
        if k == "uids":
            continue
        print("    %-16s %-10s %s bp | %s" % (v.get("accessionversion"), v.get("slen"), v.get("slen"),
                                              (v.get("title") or "")[:80]))
except Exception as e:
    print("    esummary 失败: %s | %s" % (str(e)[:60], blob[:200]))

print("\n" + "=" * 120)
print("4. HG999358 被命中的基因是什么蛋白（取 GenBank feature table 注释）")
print("=" * 120)
# 从 header 里抽出坐标: >HG999358.1_933 # 984498 # 986468 # 1 # ...
need = {}
for ln in open(W / "dia.tsv"):
    pass
coords = set()
for r in rows:
    if r[1].startswith("HG999358"):
        m = re.search(r"#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?\d+)", r[1])
        if m:
            coords.add((r[1].split()[0], int(m.group(1)), int(m.group(2)), m.group(3)))
ft = sh('curl -s -g "%s/efetch.fcgi?db=nuccore&id=HG999358&rettype=ft&retmode=text"' % EUTILS)
print("feature table: %d 字节" % len(ft))
blocks = ft.split("\n     CDS")
orfs = []
for i, b in enumerate(ft.split("\n     ")):
    m = re.search(r"CDS\s+(\d+)\.\.(\d+)", b)
    if not m:
        continue
    prod = re.search(r'/product="([^"]+)"', b)
    gene = re.search(r'/gene="([^"]+)"', b)
    orfs.append((int(m.group(1)), int(m.group(2)), prod.group(1) if prod else "-",
                 gene.group(1) if gene else "-", b.replace("\n", " ")[:150]))
print("解析到 %d 个 CDS" % len(orfs))
for tag, s, e, strand in sorted(coords):
    hit = [o for o in orfs if abs(o[0] - s) < 2000]
    print("  %-18s %d-%d (%s) -> %s" % (tag.split("_")[-1] if "_" in tag else tag, s, e, strand,
                                        (hit[0][2] + " | gene=" + hit[0][3]) if hit else "未匹配"))
