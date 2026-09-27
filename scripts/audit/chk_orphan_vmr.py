#!/usr/bin/env python3
"""让 VMR 裁决被置空的属是不是真信息: 属在 VMR 里归哪个科, 与共识科一致否"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
FNAME = "final_integrated_classification.tsv"
BASE = Path("/tmp/rc_base") / FNAME
R1 = Path("/tmp/rc_patch") / FNAME


def rd(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def blank(v):
    return v is None or v.strip() in ("", "NA", "na", "N/A", "nan", "None")


hdr = [c for c in csv.reader(open(VMR, newline="", encoding="utf-8", errors="replace"),
                            delimiter="\t").__next__()]
print("VMR 列: %s" % hdr[:20])
fam_col = "Family" if "Family" in hdr else None
gen_col = "Genus" if "Genus" in hdr else None
gtype_col = "Genome" if "Genome" in hdr else None
assert fam_col and gen_col, "VMR 缺 Family/Genus 列"

g2fam = defaultdict(Counter)          # genus(lower) -> Counter(family)
g2genome = defaultdict(Counter)       # genus(lower) -> Counter(Genome)
g_disp = {}
vmr_rows = 0
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        vmr_rows += 1
        g = (r.get(gen_col) or "").strip()
        if blank(g):
            continue
        k = g.lower()
        g_disp.setdefault(k, g)
        fm = (r.get(fam_col) or "").strip()
        g2fam[k][fm if not blank(fm) else "(空科)"] += 1
        if gtype_col:
            g2genome[k][(r.get(gtype_col) or "").strip()] += 1
print("VMR 行数 %d, 有属行 %d, 唯一属 %d\n" % (vmr_rows, sum(sum(c.values()) for c in g2fam.values()), len(g2fam)))

base = {r["contig_id"]: r for r in rd(BASE)}
r1 = {r["contig_id"]: r for r in rd(R1)}

blanked, replaced, same, na_before = [], [], 0, 0
for cid, rb in base.items():
    gb = rb.get("Genus", "")
    rr = r1.get(cid)
    if rr is None:
        continue
    gr = rr.get("Genus", "")
    if blank(gb):
        na_before += 1
        continue
    if blank(gr):
        blanked.append((cid, rb, rr))
    elif gb.strip().lower() != gr.strip().lower():
        replaced.append((cid, rb, rr))
    else:
        same += 1

print("base 属非空 %d -> 置空 %d / 替换 %d / 不变 %d\n"
      % (len(base) - na_before, len(blanked), len(replaced), same))

print("=== A. 被置空的属: VMR 怎么判 ===")
cat = Counter()
detail = defaultdict(Counter)
for cid, rb, rr in blanked:
    g = rb["Genus"].strip()
    cfam = (rb.get("Family") or "").strip()
    c = g2fam.get(g.lower())
    if not c:
        cat["① VMR 未收录该属 -> 无据可依"] += 1
        detail["①"][cfam] += 1
    elif any(f.lower() == cfam.lower() for f in c):
        cat["② VMR 认可该属且归在共识科 -> 真信息丢失"] += 1
        detail["②"][cfam] += 1
    else:
        cat["③ VMR 认可该属但归在别的科 -> 矛盾被证实"] += 1
        top = c.most_common(1)[0][0]
        detail["③"]["%s -> VMR: %s" % (cfam, top)] += 1
for k in sorted(cat):
    print("  %-38s %5d" % (k, cat[k]))
print("\n  ③ 组 top10 (共识科 -> VMR 归科):")
for k, v in detail["③"].most_common(10):
    print("      %-52s %4d" % (k, v))
print("\n  ① 组 top8 共识科:")
for k, v in detail["①"].most_common(8):
    print("      %-52s %4d" % (k, v))
print("\n  ② 组 top8 共识科:")
for k, v in detail["②"].most_common(8):
    print("      %-52s %4d" % (k, v))

print("\n  被置空属 VMR genome 类型分布 top8:")
gk = Counter()
for cid, rb, rr in blanked:
    c = g2genome.get(rb["Genus"].strip().lower())
    if c:
        for k, v in c.items():
            gk[k] += 1
for k, v in gk.most_common(8):
    print("      %-24s %4d" % (k or "(空)", v))

print("\n=== B. 被替换的属: 新属在 VMR 里归哪个科 ===")
cat2 = Counter()
for cid, rb, rr in replaced:
    ng = rr["Genus"].strip()
    cfam = (rb.get("Family") or "").strip()
    c = g2fam.get(ng.lower())
    if not c:
        cat2["① VMR 未收录新属"] += 1
    elif any(f.lower() == cfam.lower() for f in c):
        cat2["② VMR: 新属归在共识科 (改对了)"] += 1
    else:
        cat2["③ VMR: 新属归别的科 (仍不符)"] += 1
for k in sorted(cat2):
    print("  %-34s %5d" % (k, cat2[k]))

print("\n=== C. 被置空行的 Species 是否仍写着该属 (双名首词) ===")
rec = 0
for cid, rb, rr in blanked:
    sp = (rb.get("Species") or "").strip()
    g = rb["Genus"].strip()
    if not blank(sp) and " " in sp and sp.split(" ", 1)[0].lower() == g.lower():
        rec += 1
print("  Species 首词 == 被置空属: %d / %d  (%.2f%%)" % (rec, len(blanked), 100.0 * rec / max(1, len(blanked))))
