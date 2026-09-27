"""科属不一致现状实测: 当前植物表里"科判原生/藻 + 属判植物"的行还在不在, 有多少."""
import csv
import os
from collections import Counter

C7 = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
K = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus"
G = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"


def load(tbl, col):
    d = {}
    with open(f"{C7}/{tbl}", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[(r.get(col) or "").strip()] = (r.get("Predicted_Host") or "").strip()
    return d


fam_h = load("family_host_probability.tsv", "Family")
gen_h = load("genus_host_probability.tsv", "Genus")
PA_F = {k for k, v in fam_h.items() if v in ("Protist", "Algae")}
PA_G = {k for k, v in gen_h.items() if v in ("Protist", "Algae")}
PLANT_G = {k for k, v in gen_h.items() if v == "Plant"}
PLANT_F = {k for k, v in fam_h.items() if v == "Plant"}
print(f"C7: 原生/藻科 {len(PA_F)} | 原生/藻属 {len(PA_G)} | 植物属 {len(PLANT_G)} | 植物科 {len(PLANT_F)}")

for label, root in (("onekp", K), ("goji", G)):
    p = f"{root}/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"
    print("\n" + "=" * 92)
    print(f"[{label}] {p}")
    print("=" * 92)
    if not os.path.exists(p):
        print("  不存在")
        continue
    import time as _t
    print("  mtime:", _t.strftime("%Y-%m-%d %H:%M", _t.localtime(os.path.getmtime(p))),
          "| 字节:", os.path.getsize(p))
    rows = list(csv.DictReader(open(p, errors="replace"), delimiter="\t"))
    print("  行数:", len(rows))
    print("  列:", list(rows[0].keys())[:26])
    n_fam_pa = 0
    n_conflict = 0
    combo = Counter()
    case = Counter()
    for r in rows:
        fm = (r.get("Family") or "").strip()
        g = (r.get("Genus") or "").strip()
        if fm not in PA_F:
            continue
        n_fam_pa += 1
        if g in PLANT_G:
            n_conflict += 1
            combo[(fm, g)] += 1
        elif g in PA_G:
            case["科=原生/藻, 属也是原生/藻(自洽)"] += 1
        elif not g or g in ("NA", "nan", "Unclassified"):
            case["科=原生/藻, 属缺失"] += 1
        else:
            case["科=原生/藻, 属=其他"] += 1
    print(f"  科字段∈原生/藻科(16): {n_fam_pa} 行")
    print(f"  其中 科属打架(科判原生/藻 + 属∈植物病毒属): {n_conflict} 行")
    print("  其余形态:", dict(case))
    print("  打架 top 组合:")
    for k, v in combo.most_common(12):
        print(f"     {k[0]:<20} x {k[1]:<24} {v}")

    # 反方向: 植物科 + 属判非植物
    rev = 0
    revc = Counter()
    for r in rows:
        fm = (r.get("Family") or "").strip()
        g = (r.get("Genus") or "").strip()
        if fm in PLANT_F and g in gen_h and gen_h[g] not in ("Plant",):
            rev += 1
            revc[(fm, g, gen_h[g])] += 1
    print(f"  反方向(植物科 + 属判非植物): {rev} 行")
    for k, v in revc.most_common(6):
        print(f"     {k[0]:<20} x {k[1]:<22} 属判{k[2]:<10} {v}")
