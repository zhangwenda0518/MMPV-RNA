#!/usr/bin/env python3
"""Biavirus 记录全字段 + 各工具原始命中（这个 'Guapo partitivirus' 到底是什么）"""
import csv
from pathlib import Path

DT = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
CLASSED = DT / "05_Taxonomy/Votus.classed"
TARGETS = [
    "CRR1440136_clean_NODE_2764_length_1347_cov_11.565149",
    "CRR527046_clean_NODE_560_length_1834_cov_106.714935",
]


def clean(v):
    return (v or "").strip().strip('"')


print("=" * 100)
print("【1】VMR 中 Biavirus 那条记录的 29 个字段全值")
print("=" * 100)
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        if clean(r.get("Genus")) == "Biavirus":
            for k, v in r.items():
                print("  %-32s %s" % (k, clean(v)))
print()

print("=" * 100)
print("【2】VMR 中含 raunefjorden / Schizomimivirus 的记录")
print("=" * 100)
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        blob = " ".join(clean(v) for v in r.values()).lower()
        if "raunefjorden" in blob or "schizomimi" in blob:
            print("  Family=%-20s Genus=%-16s Species=%-34s accession=%s" %
                  (clean(r.get("Family")), clean(r.get("Genus")), clean(r.get("Species")),
                   clean(r.get("Virus GENBANK accession"))))
print()

for tf in sorted(CLASSED.glob("Votus_*_taxonomy.tsv")):
    name = tf.name.replace("Votus_", "").replace("_taxonomy.tsv", "")
    with open(tf, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        try:
            hdr = next(rd)
        except StopIteration:
            continue
        rows = [r for r in rd if r and any(t in r[0] for t in TARGETS)]
    if not rows:
        continue
    print("=" * 100)
    print("【3】%s 原始表 (%d 列): %s" % (name, len(hdr), " | ".join(hdr)))
    print("=" * 100)
    for r in rows:
        print("  行: %s" % " | ".join(r))
    print()

print("=" * 100)
print("【4】全库检索 'Guapo'（VMR 之外，看是否出现在别的参考里）")
print("=" * 100)
for p in [DT / "05_Taxonomy/Votus.classed", DT / "05_Taxonomy/Votus.integrated"]:
    for tf in sorted(p.glob("*.tsv")):
        with open(tf, encoding="utf-8", errors="replace") as f:
            hits = [ln.rstrip("\n") for ln in f if "guapo" in ln.lower()]
        if hits:
            print("  %s : %d 行" % (tf.name, len(hits)))
            for h in hits[:3]:
                print("      %s" % (h[:400]))
