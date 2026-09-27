"""植物表内科属标签打架的行 + 全库 Plant->Protist 13 行归因。只读。"""
import csv
import os
from collections import Counter

CA = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
RUN = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
PLANT = f"{RUN}/06_HostPrediction/C9_ICTV_result/Plant.classified.tsv"
FULLHOST = f"{RUN}/06_HostPrediction/C9_ICTV_result/classification_result.tsv"
NATIVE = {"Protist", "Algae"}


def rows_of(path):
    with open(path, errors="replace") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            yield {(k.strip() if k else k): (v.strip().strip('"') if isinstance(v, str) else v)
                   for k, v in d.items()}


def load_prob(level):
    out = {}
    for r in rows_of(f"{CA}/{level.lower()}_host_probability.tsv"):
        n = r.get(level, "")
        if n:
            out[n] = r
    return out


fam_p, gen_p = load_prob("Family"), load_prob("Genus")


def lab(tb, name):
    return tb.get(name, {}).get("Predicted_Host", "")


print("=" * 90)
print("[A] 全库 Plant->Protist / Protist->Plant 类不一致的 13 行归因")
print("=" * 90)
buckets = Counter()
rows13 = []
for r in rows_of(FULLHOST):
    fh = lab(fam_p, r.get("Family", ""))
    gh = lab(gen_p, r.get("Genus", ""))
    if fh and gh and fh != gh and (fh in NATIVE or gh in NATIVE):
        buckets[(fh, gh, r.get("Family", ""), r.get("Genus", ""))] += 1
for (fh, gh, fam, gen), v in buckets.most_common(20):
    print(f"  {v:>5} 行  科={fam:<24}({fh})  属={gen:<22}({gh})")

print("\n" + "=" * 90)
print("[B] 植物表 1,065 行内部: 科标签 vs 属标签")
print("=" * 90)
c = Counter()
no_gen = no_fam = 0
for r in rows_of(PLANT):
    fh = lab(fam_p, r.get("Family", ""))
    gh = lab(gen_p, r.get("Genus", ""))
    if not fh:
        no_fam += 1
    if not gh:
        no_gen += 1
    c[(fh or "缺科条目", gh or "缺属条目")] += 1
tot = sum(c.values())
print(f"  总行 {tot}; 缺科级条目 {no_fam}; 缺属级条目 {no_gen}")
print(f"  {'科标签':<14} {'属标签':<14} {'行数':>6}  占比")
for (fh, gh), v in c.most_common():
    print(f"  {fh:<14} {gh:<14} {v:>6}  {v / tot * 100:5.1f}%")

n_conf = sum(v for (fh, gh), v in c.items() if fh != "缺科条目" and gh != "缺属条目" and fh != gh)
n_ok = sum(v for (fh, gh), v in c.items() if fh != "缺科条目" and gh != "缺属条目" and fh == gh)
print(f"\n  可比对(两侧都有条目) {n_conf + n_ok} 行: 一致 {n_ok}, 不一致 {n_conf}")
print("  不一致且涉原生/藻:",
      sum(v for (fh, gh), v in c.items() if fh != gh and (fh in NATIVE or gh in NATIVE)))

print("\n  不一致明细 (前 25 行):")
n = 0
for r in rows_of(PLANT):
    fh = lab(fam_p, r.get("Family", ""))
    gh = lab(gen_p, r.get("Genus", ""))
    if fh and gh and fh != gh:
        n += 1
        if n <= 25:
            print(f"    {r.get('contig_id','')[:40]:<40} 科={r.get('Family',''):<22}({fh:<10}) "
                  f"属={r.get('Genus',''):<20}({gh})")
print(f"    ... 共 {n} 行")
