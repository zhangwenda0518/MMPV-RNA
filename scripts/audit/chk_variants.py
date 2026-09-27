#!/usr/bin/env python3
"""四条口径对撞: base / R1-last / R2-last / R1-first / R2-first
指标: 属空缺数、科属核酸矛盾行数、相对 base 的属改动、不变量违反行数
"""
import csv
import importlib.util
from collections import defaultdict
from pathlib import Path

ORIG = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
            "/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
FNAME = "final_integrated_classification.tsv"

TABLES = [
    ("Aug12表", ORIG),
    ("base(现行)", Path("/tmp/rc_base") / FNAME),
    ("R1-last", Path("/tmp/rc_patch") / FNAME),
    ("R2-last", Path("/tmp/rc_v_last_noblank") / FNAME),
    ("R1-first", Path("/tmp/rc_v_first_blank") / FNAME),
    ("R2-first", Path("/tmp/rc_v_first_noblank") / FNAME),
]
TOOLS = ["ACVirus", "CAT", "diamond_lca", "genomad", "metabuli", "mmseqs", "VITAP"]

spec = importlib.util.spec_from_file_location("ann", PIPE / "utils" / "annotate_nucleic_acid.py")
ann = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ann)
IDX = ann.load_vmr(VMR)
B = ann.is_blank


def broad(row, rank):
    v = row.get(rank, "")
    if B(v):
        return None
    return ann.resolve_rank(v, rank, IDX)


def load(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


per_tool = defaultdict(dict)
for t in TOOLS:
    with open(Path("/tmp/rc_base") / ("standardized_%s.tsv" % t), newline="",
              encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            per_tool[r.get("contig_id")][t] = r

tabs = {}
for name, p in TABLES:
    if not p.exists():
        print("缺失: %s" % p)
        continue
    tabs[name] = load(p)
print("装载: %s\n" % ", ".join("%s=%d行" % (k, len(v)) for k, v in tabs.items()))

# 以 base 的共识 Family 建「自洽工具集候选属」(Family 各表应一致)
base = tabs["base(现行)"]
cand = {}
for r in base:
    cid, F = r["contig_id"], r.get("Family", "")
    if B(F):
        continue
    c = set()
    for t, d in per_tool.get(cid, {}).items():
        if (d.get("Family") or "") == F and not B(d.get("Genus") or ""):
            c.add(d["Genus"].lower())
    cand[cid] = c

hdr = "%-12s %8s %8s %8s %9s %9s %9s"
print(hdr % ("口径", "属非空", "属空缺", "科空缺", "科属矛盾", "属改base", "族改base"))
print("-" * 72)
rows_out = {}
for name, rows in tabs.items():
    g_ok = sum(1 for r in rows if not B(r.get("Genus", "")))
    f_ok = sum(1 for r in rows if B(r.get("Family", "")))
    contra = 0
    for r in rows:
        bf, bg = broad(r, "Family"), broad(r, "Genus")
        if bf and bg and bf != bg:
            contra += 1
    gdiff = fdiff = 0
    if name != "base(现行)":
        bm = {r["contig_id"]: r for r in base}
        for r in rows:
            br = bm.get(r["contig_id"])
            if br is None:
                continue
            if (r.get("Genus") or "") != (br.get("Genus") or ""):
                gdiff += 1
            if (r.get("Family") or "") != (br.get("Family") or ""):
                fdiff += 1
    else:
        bm = {r["contig_id"]: r for r in tabs["Aug12表"]}
        for r in rows:
            br = bm.get(r["contig_id"])
            if br and (r.get("Genus") or "") != (br.get("Genus") or ""):
                gdiff += 1
            if br and (r.get("Family") or "") != (br.get("Family") or ""):
                fdiff += 1
    print(hdr % (name, g_ok, len(rows) - g_ok, f_ok, contra, gdiff, fdiff))
    rows_out[name] = rows

print()
print("=== 不变量检查: 属非空但无任何「认同该行科」的工具报出该属 ===")
for name, rows in rows_out.items():
    bad = 0
    for r in rows:
        cid, G, F = r["contig_id"], r.get("Genus", ""), r.get("Family", "")
        if B(G) or B(F):
            continue
        cs = cand.get(cid, set())
        if G.lower() not in cs:
            bad += 1
    print("  %-12s 违反 %5d" % (name, bad))

print()
print("=== 306 打架行 (科属核酸相反) 的终态 ===")
clash = set()
for r in tabs["Aug12表"]:
    bf, bg = broad(r, "Family"), broad(r, "Genus")
    if bf and bg and bf != bg:
        clash.add(r["contig_id"])
print("  Aug12 表判定为打架的 contig: %d" % len(clash))
for name, rows in rows_out.items():
    still = 0
    for r in rows:
        if r["contig_id"] not in clash:
            continue
        bf, bg = broad(r, "Family"), broad(r, "Genus")
        if bf and bg and bf != bg:
            still += 1
    print("  %-12s 仍打架 %5d" % (name, still))
