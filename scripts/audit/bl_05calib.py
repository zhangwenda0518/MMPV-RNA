"""校准表 calib_* 列交叉核对 (review 行是否保留 Genus/Species)."""
import csv
from collections import Counter

P = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
     "/05_Taxonomy/Votus.integrated/calibration_20260914/final_integrated_classification.calibrated.tsv")


def v(x):
    return (x or "").strip().strip('"')


rows = list(csv.DictReader(open(P, errors="replace"), delimiter="\t"))
print("总行数:", len(rows))
ct = Counter()
for r in rows:
    act = v(r.get("calib_action"))
    gen = v(r.get("Genus"))
    sp = v(r.get("Species"))
    fam = v(r.get("Family"))
    key = (act or "(none)",
           "Genus空" if gen in ("", "NA") else "Genus有值",
           "Species空" if sp in ("", "NA") else "Species有值",
           "Family空" if fam in ("", "NA") else "Family有值")
    ct[key] += 1
print("\naction | Genus | Species | Family  分布:")
for k, n in sorted(ct.items(), key=lambda x: (-x[1], x[0])):
    print(f"  {k[0]:<18} {k[1]:<10} {k[2]:<12} {k[3]:<11} {n}")

print("\nreview 行的属名样例(前 8):")
n = 0
for r in rows:
    if v(r.get("calib_action")) == "review":
        print("   ", v(r.get("Family")), "|", v(r.get("Genus")), "|", v(r.get("calib_flag_single_ref")))
        n += 1
        if n >= 8:
            break

print("\nblank 行的属名(前 5, 应已置空, 原值在 calib_prev_genus):")
n = 0
for r in rows:
    if v(r.get("calib_action")) == "blank":
        print("   Fam=", v(r.get("Family")), "| Genus=", repr(v(r.get("Genus"))),
              "| prev=", v(r.get("calib_prev_genus")))
        n += 1
        if n >= 5:
            break
