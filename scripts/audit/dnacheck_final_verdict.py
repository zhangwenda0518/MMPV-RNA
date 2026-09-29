#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""枸杞 DNA 病毒真伪终裁：合并 minimap2 / dc-megablast / 经典 blastn 三路证据。
基因组同源证据分强/中/弱三档；无同源者为真 DNA 病毒候选，再看检出广度。"""
import csv, json
from collections import Counter, defaultdict

OUT = "/home/zhangwenda/goji_dnacheck_20260929"

# ---------- 读 manifest（500 条） ----------
rows = []
with open(f"{OUT}/dna_manifest.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        r["q"] = f"{r['Cohort']}__{r['contig_id']}"
        rows.append(r)

# ---------- minimap2 round1 ----------
mm = {}
with open(f"{OUT}/genome_match_round1.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        mm[f"{r['Cohort']}__{r['contig_id']}"] = r

# ---------- blast 两轮 best ----------
def best_hits(path, has_eval=False):
    best = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        q, t = p[0], p[1]
        pid, aln = float(p[2]), int(p[3])
        sc = float(p[8])
        cur = best.get(q)
        if cur is None or (pid, aln) > (cur[1], cur[2]):
            best[q] = (t, pid, aln, sc)
    return best

dc = best_hits(f"{OUT}/t3_all_dcmega.tsv")
cl = best_hits(f"{OUT}/t3_all_classic.tsv")

SAME = {"Lycium_barbarum": "ningxia", "Lycium_chinense": "zhonghua",
        "Lycium_ruthenicum": "heiguo"}
NOSPECIES = {"Lycium_amarum", "Aphis_gossypii", "Fusarium_nematophilum", "Alternaria_alternata"}


def merge_evidence(q):
    """返回 (级别, 描述)：EVE_strong / EVE_moderate / EVE_weak / none"""
    cands = []
    m = mm.get(q)
    if m and m.get("genome_hit") not in ("none", ""):
        cands.append((float(m["pid"]), 10**6, f"minimap2 {m['genome_hit']} {m['hit_locus']}"))
    for tag, d in (("dc", dc), ("bl", cl)):
        if q in d:
            t, pid, aln, sc = d[q]
            cands.append((pid, aln, f"{tag} {t} pid={pid:.1f} aln={aln}"))
    if not cands:
        return "none", ""
    # 取证据最强的一条（先按 pid 分档，再按 aln）
    def key(c):
        lvl = 3 if c[0] >= 90 else 2 if c[0] >= 80 else 1 if c[0] >= 55 else 0
        return (lvl, min(c[1], 10**5))
    cands.sort(key=key, reverse=True)
    pid, aln, desc = cands[0]
    if pid >= 90 and aln >= 200:
        return "EVE_strong", desc
    if pid >= 80 and aln >= 200:
        return "EVE_moderate", desc
    if pid >= 70 and aln >= 300:
        return "EVE_moderate", desc
    if pid >= 55 and aln >= 150:
        return "EVE_weak", desc
    if pid >= 55 and aln >= 100:
        return "EVE_weak", desc
    return "none", desc if desc else (";".join(c[2] for c in cands[:1]))


for r in rows:
    lvl, desc = merge_evidence(r["q"])
    r["genome_evidence"] = lvl
    r["evidence_detail"] = desc
    if r["Cohort"] in NOSPECIES:
        r["verdict"] = "无自身参考基因组-无法用基因组排除"
    elif lvl == "EVE_strong":
        r["verdict"] = "基因组强同源(近期整合EVE/同源拷贝)"
    elif lvl == "EVE_moderate":
        r["verdict"] = "基因组中度同源(疑似EVE)"
    elif lvl == "EVE_weak":
        r["verdict"] = "基因组弱同源(远缘EVE信号)"
    else:
        r["verdict"] = "无基因组同源(真病毒候选)"

# 检出广度标签
for r in rows:
    try:
        n = int(float(r["n_samples"]))
    except Exception:
        n = -1
    r["n"] = n

with open(f"{OUT}/final_dna_verdict.tsv", "w") as fo:
    cols = ["Cohort", "contig_id", "Genus", "Species", "Family_true", "Category", "n_samples",
            "length", "genome_evidence", "verdict", "evidence_detail"]
    fo.write("\t".join(cols) + "\n")
    for r in rows:
        fo.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

print("=== 总体判定 ===")
for k, v in Counter(r["verdict"] for r in rows).most_common():
    print(f"  {k}: {v}")
print()
lyc = [r for r in rows if r["Cohort"] not in NOSPECIES]
print(f"=== 三个有参考基因组枸杞种 ({len(lyc)} 条) ===")
for k, v in Counter(r["verdict"] for r in lyc).most_common():
    print(f"  {k}: {v}")
print()
print("=== 真病毒候选（无基因组同源）按科 ===")
cand = [r for r in lyc if r["verdict"] == "无基因组同源(真病毒候选)"]
for k, v in Counter(r["Family_true"] for r in cand).most_common():
    print(f"  {k}: {v}")
print()
print("=== 真病毒候选按队列 ===")
for k, v in Counter(r["Cohort"] for r in cand).most_common():
    print(f"  {k}: {v}")
print()
print("=== 真病毒候选按属 ===")
for k, v in Counter(r["Genus"] for r in cand).most_common():
    print(f"  {k}: {v}")
print()
ns = [r for r in cand if r["n"] >= 2]
print(f"真病毒候选中 n_samples>=2 的: {len(ns)}")
for r in ns[:20]:
    print(f"  {r['Cohort']} {r['contig_id'][:60]} n={r['n_samples']} {r['Genus']} {r['Species'][:40]} len={r['length']}")
