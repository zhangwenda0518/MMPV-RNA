"""原生(Protist/Algae)与植物病毒能否按科/属机械分开 —— 覆盖率与盲区实测。

数据源:
  1. database/cross_analysis/{order,family,genus,species}_host_probability.tsv (C7 表)
  2. 06_HostPrediction/C9_ICTV_result/Plant.classified.tsv  (植物分支实际输入)
  3. 05_Taxonomy/Votus.integrated/final_integrated_classification.tsv (全库)
只读, 不改任何产物。
"""
import csv
import os
from collections import Counter

CA = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
RUN = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
PLANT = f"{RUN}/06_HostPrediction/C9_ICTV_result/Plant.classified.tsv"
FULL = f"{RUN}/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"

CATS = ["Algae", "Animal_other", "Arachnida", "Archaea", "Aves", "Bacteria",
        "Fungi", "Human", "Insecta", "Mammalia", "Oomycetes", "Plant", "Protist"]
NATIVE = {"Protist", "Algae"}


def rows_of(path):
    with open(path, errors="replace") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            yield {(k.strip() if k else k): (v.strip().strip('"') if isinstance(v, str) else v)
                   for k, v in d.items()}


def load_prob(level):
    out = {}
    p = f"{CA}/{level.lower()}_host_probability.tsv"
    if not os.path.exists(p):
        return out
    for r in rows_of(p):
        n = r.get(level, "")
        if n:
            out[n] = r
    return out


def fl(r, k, d=0.0):
    try:
        return float(r.get(k) or d)
    except Exception:
        return d


def line(ch="-", n=78):
    print(ch * n)


fam_p = load_prob("Family")
gen_p = load_prob("Genus")
sp_p = load_prob("Species")
ord_p = load_prob("Order")

print("=" * 78)
print("[1] C7 概率表规模与宿主类别")
line("=")
print(f"  Order {len(ord_p):>6}   Family {len(fam_p):>6}   Genus {len(gen_p):>6}   Species {len(sp_p):>6}")
for lvl, tb in (("Family", fam_p), ("Genus", gen_p)):
    c = Counter(v["Predicted_Host"] for v in tb.values())
    print(f"  {lvl:<7} Predicted_Host 分布: " +
          "  ".join(f"{k}={c[k]}" for k in CATS if c.get(k)) +
          (f"  Other={sum(v for k, v in c.items() if k not in CATS)}" if any(k not in CATS for k in c) else ""))

# ---- 原生/藻类在概率表里的名单 ----
nat_gen = sorted(k for k, v in gen_p.items() if v["Predicted_Host"] in NATIVE)
nat_fam = sorted(k for k, v in fam_p.items() if v["Predicted_Host"] in NATIVE)
print(f"\n  属级条目判为原生/藻类的属: {len(nat_gen)} 个")
print(f"  科级条目判为原生/藻类的科: {len(nat_fam)} 个 -> {', '.join(nat_fam)}")
print(f"  Biavirus 在属级表中: {'有' if 'Biavirus' in gen_p else '无 (缺失)'}"
      f"   在科级表中: {'有' if 'Biavirus' in fam_p else '无'}")
print(f"  Partitiviridae 科级: host={fam_p.get('Partitiviridae', {}).get('Predicted_Host')} "
      f"P_Max={fam_p.get('Partitiviridae', {}).get('P_Max')} "
      f"records={fam_p.get('Partitiviridae', {}).get('Total_Records')} "
      f"P(Plant)={fam_p.get('Partitiviridae', {}).get('P(Plant)')} "
      f"P(Protist)={fam_p.get('Partitiviridae', {}).get('P(Protist)')} "
      f"P(Fungi)={fam_p.get('Partitiviridae', {}).get('P(Fungi)')}")

# ---- [2] 植物表分层 ----
prows = list(rows_of(PLANT))
print("\n" + "=" * 78)
print("[2] Plant.classified.tsv 分层")
line("=")
print(f"  行数: {len(prows)}")
det = Counter(r.get("Determination_Level", "?") for r in prows)
print("  Determination_Level:")
for k, v in det.most_common():
    print(f"    {k:<22} {v:>6}  ({v / len(prows) * 100:.1f}%)")
print("  Confidence_Level:")
for k, v in Counter(r.get("Confidence_Level", "?") for r in prows).most_common():
    print(f"    {k:<22} {v:>6}")

n_nogen = sum(1 for r in prows if not r.get("Genus"))
n_in_gen = sum(1 for r in prows if r.get("Genus") in gen_p)
n_in_fam = sum(1 for r in prows if r.get("Family") in fam_p)
print(f"\n  有属级条目(按属可判): {n_in_gen}  ({n_in_gen / len(prows) * 100:.1f}%)")
print(f"  属名为空: {n_nogen}")
print(f"  有科级条目: {n_in_fam}  ({n_in_fam / len(prows) * 100:.1f}%)")
print(f"  => 只能靠科级兜底的(属无条目): {len(prows) - n_in_gen}  ({(len(prows) - n_in_gen) / len(prows) * 100:.1f}%)")

# 属级有条目但宿主不是 Plant
conf_rows = []
for r in prows:
    g = r.get("Genus", "")
    if g in gen_p and gen_p[g]["Predicted_Host"] != "Plant":
        conf_rows.append((g, gen_p[g]["Predicted_Host"], r.get("Family", ""), r.get("contig_id", "")))
print(f"\n  属级判为非植物、却被放进 Plant 表: {len(conf_rows)}")
cg = Counter((g, h) for g, h, _, _ in conf_rows)
for (g, h), v in cg.most_common(15):
    print(f"    {g:<26} 属级host={h:<14} {v:>4} 行")

# 原生属/科直接出现在植物表
hit_nat_gen = Counter(r.get("Genus") for r in prows if r.get("Genus") in nat_gen)
hit_nat_fam = Counter(r.get("Family") for r in prows if r.get("Family") in nat_fam)
print(f"\n  植物表里属级原生名单命中: {sum(hit_nat_gen.values())} 行 {dict(hit_nat_gen)}")
print(f"  植物表里科级原生名单命中: {sum(hit_nat_fam.values())} 行 {dict(hit_nat_fam)}")

# ---- [3] 弱证据科级兜底 ----
print("\n" + "=" * 78)
print("[3] 科级兜底的证据强度 (Integrated_Confidence 分档)")
line("=")
bins = [(0.0, 0.4, "极弱 <0.4"), (0.4, 0.6, "弱 0.4-0.6"), (0.6, 0.8, "中 0.6-0.8"), (0.8, 1.01, "强 >=0.8")]
for lo, hi, lab in bins:
    n_all = sum(1 for r in prows if lo <= fl(r, "Integrated_Confidence") < hi)
    n_fam = sum(1 for r in prows if lo <= fl(r, "Integrated_Confidence") < hi
                and str(r.get("Determination_Level", "")).startswith("Family"))
    print(f"  {lab:<12} 全部 {n_all:>5}   其中靠 Family 判定 {n_fam:>5}")

# ---- [4] 跨宿主科 ----
print("\n" + "=" * 78)
print("[4] 跨宿主科: 科内既有植物记录又有原生/藻/真菌记录")
line("=")
fam_rows = Counter(r.get("Family") for r in prows)
mixed = []
for fam, n in fam_rows.items():
    v = fam_p.get(fam)
    if not v:
        continue
    pp, pr, al, fu = fl(v, "P(Plant)"), fl(v, "P(Protist)"), fl(v, "P(Algae)"), fl(v, "P(Fungi)")
    if pp > 0 and (pr + al + fu) > 0:
        mixed.append((n, fam, v["Predicted_Host"], pp, pr + al, fu, v["Total_Records"], fl(v, "Shannon_Entropy")))
mixed.sort(reverse=True)
print(f"  {'行数':>5} {'科':<26} {'科级host':<10} {'P(Plant)':>8} {'P(原生藻)':>9} {'P(Fungi)':>8} {'记录':>7} {'Shannon':>8}")
for n, fam, h, pp, pa, fu, rec, sh in mixed[:25]:
    print(f"  {n:>5} {fam:<26} {h:<10} {pp:>8.3f} {pa:>9.3f} {fu:>8.3f} {rec:>7} {sh:>8.3f}")
print(f"  合计跨宿主科: {len(mixed)}, 涉及植物表 {sum(m[0] for m in mixed)} 行 "
      f"({sum(m[0] for m in mixed) / len(prows) * 100:.1f}%)")

# ---- [5] 全库属级覆盖 ----
print("\n" + "=" * 78)
print("[5] 全库 20,892 行的属级覆盖(参照)")
line("=")
nf = ng = nall = 0
for r in rows_of(FULL):
    nall += 1
    if r.get("Genus"):
        ng += 1
    if r.get("Genus") in gen_p:
        nf += 1
print(f"  行数 {nall}   Genus 非空 {ng} ({ng / nall * 100:.1f}%)   命中属级表 {nf} ({nf / nall * 100:.1f}%)")
