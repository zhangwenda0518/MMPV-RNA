"""版本漂移定性: 7 个不在 MSL41 的科名是谁; 251 vmr_only 与 1246 blank 的真实成因."""
import csv
import os
from collections import Counter, defaultdict

MSL = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"
BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
        "/05_Taxonomy/Votus.integrated")
CAL = f"{BASE}/calibration_20260914/final_integrated_classification.calibrated.tsv"

msl_fam = set()
msl_gen2fam = defaultdict(set)
with open(MSL, errors="replace") as f:
    for r in csv.DictReader(f):
        fam = (r.get("Family") or "").strip().lower()
        gen = (r.get("Genus") or "").strip().lower()
        if fam:
            msl_fam.add(fam)
        if gen and fam:
            msl_gen2fam[gen].add(fam)

print("=== MSL41 是否保留这些科 ===")
for k in ["mimiviridae", "hydriviridae", "phycodnaviridae", "marseilleviridae", "iridoviridae",
          "pithoviridae", "schizomimiviridae", "hepaciviridae", "flaviviridae", "pestiviridae",
          "ourmiaviridae", "botourmiaviridae", "parahypoviridae", "hypoviridae", "partitiviridae"]:
    print(f"  {k:<22} {'在' if k in msl_fam else '不在'} MSL41")


def val(x):
    x = (x or "").strip().strip('"')
    return None if x in ("", "NA", "N/A", "-") else x


rows = list(csv.DictReader(open(CAL, errors="replace"), delimiter="\t"))

# ── 1. 全部 7 个非 MSL41 科名 (含无候选者) ──────────────────────
cnt = Counter()
for r in rows:
    fam = val(r.get("Family"))
    if fam and fam.lower() not in msl_fam:
        cnt[fam] += 1
print(f"\n=== Family 非 MSL41 科名 全清单 ({len(cnt)} 个) ===")
for k, v in cnt.most_common():
    print(f"  {k:<24} {v} 行")

# ── 2. 251 vmr_only 行: 科名 vs 属的 MSL41 定型科 ────────────────
pair = Counter()
for r in rows:
    if val(r.get("calib_flag_single_ref")) != "vmr_only":
        continue
    fam = val(r.get("Family"))
    g = val(r.get("Genus"))
    newf = msl_gen2fam.get((g or "").lower(), set())
    pair[(fam, ";".join(sorted(newf)) or "(属不在MSL41)")] += 1
print(f"\n=== 251 vmr_only: 行内科名 -> 属的 MSL41 定型科 ===")
for (a, b), n in pair.most_common(12):
    print(f"  {str(a):<22} -> {b:<28} {n} 行")
same = sum(n for (a, b), n in pair.items() if a and a.lower() in msl_fam)
print(f"  其中行内科名本身存在于 MSL41 的: {same} 行 (即科名未过时, 是属放错科)")
print(f"  行内科名不在 MSL41 的: {251 - same} 行 (真·科名过时)")

# ── 3. 1246 blank 行: 两参照对属的定型科是否一致 ─────────────────
ref = {}
with open(os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv"), errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        g = (r.get("Genus") or "").strip().lower()
        if g and g not in ref:
            ref[g] = ((r.get("NCBI_Family") or "").strip().lower(),
                      (r.get("VMR_Family") or "").strip().lower())
st = Counter()
for r in rows:
    if val(r.get("calib_action")) != "blank":
        continue
    g = val(r.get("calib_prev_genus"))
    ncbi, vmr = ref.get((g or "").lower(), ("", ""))
    if ncbi and vmr:
        st["两参照都有: 一致" if ncbi == vmr else "两参照都有: 不一致"] += 1
    elif ncbi or vmr:
        st["仅一侧有"] += 1
    else:
        st["两参照都没有"] += 1
print("\n=== 1246 blank 行: 属的定型科参照情况 ===")
for k, v in st.most_common():
    print(f"  {k:<20} {v} 行")

# ── 4. blank 行若是「属对科错」, 改写科的可行性 ────────────────────
tri = Counter()
for r in rows:
    if val(r.get("calib_action")) != "blank":
        continue
    fam = val(r.get("Family")) or "(行内科为空)"
    g = val(r.get("calib_prev_genus"))
    ncbi, vmr = ref.get((g or "").lower(), ("", ""))
    tri[(fam, ncbi or "-", vmr or "-")] += 1
print("\n=== blank 行 top: 行内科 | NCBI定型科 | VMR定型科 ===")
for (a, b, c), n in tri.most_common(12):
    print(f"  {a:<24} | {b:<20} | {c:<20} {n} 行")
