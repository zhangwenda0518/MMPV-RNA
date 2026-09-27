#!/usr/bin/env python3
# 名字歧义探针：对若干"行内最深名"分别列 dmp 精确同名行、其 Viruses 谱系 tuple、以及 n2r 最后写入的值
# 用法: python3 /tmp/probe_name_ambiguity.py
import sys, os, subprocess

sys.path.insert(0, "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
import virus_classifier as vc

TAXDB = os.path.expanduser("~/database/taxonomy/fullnamelineage.dmp")
NAMES = [
    "environmental samples", "unclassified viruses", "uncultured virus",
    "Alphachlorovirus", "Betachlorovirus", "Gammachlorovirus",
    "Chlorovirus newyorkense", "Chlorovirus americanus",
    "Bracoviriform", "Polydnaviriformidae", "Zimmerviridae",
    "Monodnaviria", "Floreoviria", "Organic Lake phycodnavirus",
    "Pantevenvirales", "Cotesia vestalis bracovirus",
]


def exact_lines(name):
    out = subprocess.run(["grep", "-F", name, TAXDB], capture_output=True, text=True).stdout
    res = []
    for l in out.split("\n"):
        if not l.strip():
            continue
        p = l.split("\t|\t")
        if len(p) >= 3 and p[1].strip().rstrip(";| ") == name:
            res.append((p[0].strip(), p[2].strip().rstrip(";| ")))
    return res


def tup(lin):
    r = vc.lineage_to_ranks(lin)
    return tuple(r[k] for k in vc.RANK_NAMES)


for name in NAMES:
    ex = exact_lines(name)
    viral = [(t, l) for t, l in ex if "Viruses" in l]
    tups = {}
    for t, l in viral:
        tups.setdefault(tup(l), []).append(t)
    print("### %s  精确同名行=%d  含Viruses=%d  不同tuple=%d" % (name, len(ex), len(viral), len(tups)))
    for tp, ids in list(tups.items())[:6]:
        print("    x%-4d %s" % (len(ids), " | ".join(x if x != "NA" else "." for x in tp)))
    if len(tups) > 1:
        print("    !!! 歧义：n2r 最终保留最后一行 = %s" % (viral[-1][1][:110] if viral else "none"))
    print()
