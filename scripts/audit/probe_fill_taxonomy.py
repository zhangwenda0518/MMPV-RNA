#!/usr/bin/env python3
# 探针：定位 onekp 中 Family=Adenoviridae + Genus=Organic Lake phycodnavirus 的成因
# 用法: python3 /tmp/probe_fill_taxonomy.py
import sys, os, subprocess

sys.path.insert(0, "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
import virus_classifier as vc

RANK = list(vc.RANK_NAMES)
HDR = ["seq_name", "tool"] + RANK

DIA = ("Viruses; Varidnaviria; Bamfordvirae; Nucleocytoviricota; Megaviricetes; "
       "Algavirales; Phycodnaviridae; environmental samples; Organic Lake phycodnavirus")
CAT = ("Viruses;Varidnaviria;Bamfordvirae;Nucleocytoviricota;Megaviricetes;"
       "Algavirales;Phycodnaviridae;-;Organic Lake phycodnavirus")


def row(cid, tool, lin):
    r = vc.lineage_to_ranks(lin)
    return [cid, tool] + [r.get(k, "NA") for k in RANK]


print("== 1 lineage_to_ranks 单独结果")
print("RANK_NAMES:", RANK)
print("dia:", "\t".join(row("c1", "diamond_lca", DIA)))
print("cat:", "\t".join(row("c1", "CAT", CAT)))

p = "/tmp/probe_combined.tsv"
with open(p, "w") as f:
    f.write("\t".join(HDR) + "\n")
    f.write("\t".join(row("c1", "diamond_lca", DIA)) + "\n")
    f.write("\t".join(row("c1", "CAT", CAT)) + "\n")
print("== 2 fill 之前")
print(open(p).read().strip())
vc.fill_taxonomy_na(p, p)
print("== 3 fill 之后")
print(open(p).read().strip())

print("== 4 dmp 中同名 environmental samples 的分布")
taxdb = os.path.expanduser("~/database/taxonomy/fullnamelineage.dmp")
out = subprocess.run(["grep", "-F", "environmental samples", taxdb],
                     capture_output=True, text=True).stdout
lines = [l for l in out.split("\n") if l.strip()]
exact = []
for l in lines:
    parts = l.split("\t|\t")
    if len(parts) >= 3 and parts[1].strip().rstrip(";| ") == "environmental samples":
        exact.append(parts)
print("grep 命中行数:", len(lines), " 精确同名行数:", len(exact))
agg = {}
for parts in exact:
    lin = parts[2].strip().rstrip(";| ")
    agg[lin] = agg.get(lin, 0) + 1
for i, (lin, n) in enumerate(sorted(agg.items(), key=lambda x: -x[1])[:10]):
    print("  同名谱系 x%d: %s" % (n, lin[:150]))

viral = [p2 for p2 in exact if "Viruses" in p2[2]]
print("  其中谱系含 Viruses 的同名行数:", len(viral))
if viral:
    for p2 in viral[:5]:
        r = vc.lineage_to_ranks(p2[2].strip().rstrip(";| "))
        print("   taxid=%s -> %s" % (p2[0].strip(), {k: v for k, v in r.items() if v != "NA"}))
    last = viral[-1]
    r = vc.lineage_to_ranks(last[2].strip().rstrip(";| "))
    print("  最后一个（会覆盖 n2r 键）taxid=%s -> %s" % (last[0].strip(), r))

print("== 5 单名 organic lake phycodnavirus 的映射")
out2 = subprocess.run(["grep", "-F", "Organic Lake phycodnavirus", taxdb],
                      capture_output=True, text=True).stdout
for l in out2.split("\n")[:5]:
    if not l.strip():
        continue
    parts = l.split("\t|\t")
    if len(parts) >= 3:
        r = vc.lineage_to_ranks(parts[2].strip().rstrip(";| "))
        print("  taxid=%s name=%s -> %s" % (parts[0].strip(), parts[1].strip(), r))
