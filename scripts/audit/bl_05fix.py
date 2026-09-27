"""三条边界的可修复性量化.

边界一: Family 从不改写 -> 有多少行的 Family 是 MSL41 里查不到的老科名, 能否映射到新科名
边界二: 251 行版本漂移 -> 老科名 -> 新科名的映射关系是否 1:1 可枚举
边界三: 已置空的 1,246 行 -> 科名映射后能否把属救回来; 947 行有属无科 -> 能否补科
"""
import csv
import os
from collections import Counter, defaultdict

MSL = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"
REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")
BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
        "/05_Taxonomy/Votus.integrated")
CAL = f"{BASE}/calibration_20260914/final_integrated_classification.calibrated.tsv"
ORIG = f"{BASE}/final_integrated_classification.tsv"

# ── MSL41 ──────────────────────────────────────────────────────────────
msl_fam = set()
msl_gen2fam = defaultdict(set)
msl_gen_rows = 0
with open(MSL, errors="replace") as f:
    for r in csv.DictReader(f):
        fam = (r.get("Family") or "").strip().lower()
        gen = (r.get("Genus") or "").strip()
        if fam:
            msl_fam.add(fam)
        if gen:
            msl_gen_rows += 1
            if fam:
                msl_gen2fam[gen.lower()].add(fam)
print(f"MSL41: 科 {len(msl_fam)} 个 | 属 {len(msl_gen2fam)} 个 (属行 {msl_gen_rows})")

# ── 参照表 genus -> (NCBI科, VMR科) ──────────────────────────────────────
ref = {}
with open(REF, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        g = (r.get("Genus") or "").strip().lower()
        if g and g not in ref:
            ref[g] = ((r.get("NCBI_Family") or "").strip().lower(),
                      (r.get("VMR_Family") or "").strip().lower())


def val(x):
    x = (x or "").strip().strip('"')
    return None if x in ("", "NA", "N/A", "-") else x


rows = list(csv.DictReader(open(CAL, errors="replace"), delimiter="\t"))
print(f"\n枸杞校准表: {len(rows)} 行")

fam_all = Counter()
fam_not_msl = Counter()
old2new = defaultdict(Counter)
rows_fam_not_msl = 0
unmappable = 0
for r in rows:
    fam = val(r.get("Family"))
    if fam:
        fam_all[fam] += 1
        if fam.lower() not in msl_fam:
            fam_not_msl[fam] += 1
            rows_fam_not_msl += 1
            g = val(r.get("Genus")) or val(r.get("calib_prev_genus"))
            cands = msl_gen2fam.get((g or "").lower(), set())
            if cands:
                for c in cands:
                    old2new[fam][c] += 1
            else:
                unmappable += 1

print(f"  非空 Family: {sum(fam_all.values())} 行, 去重 {len(fam_all)} 个科名")
print(f"  MSL41 查不到的科名(老科名候选): {len(fam_not_msl)} 个, 覆盖 {rows_fam_not_msl} 行")
print(f"  其中属也不在 MSL41(无法映射): {unmappable} 行")
print("\n  老科名 -> MSL41 新科名 (按行数):")
for old, c in sorted(old2new.items(), key=lambda x: -sum(x[1].values()))[:15]:
    tot = sum(c.values())
    top = ", ".join(f"{k}({v})" for k, v in c.most_common(3))
    print(f"    {old:<28} {fam_not_msl[old]:>5} 行 -> {top}")

# ── 边界三 A: 被置空的属能否靠科名映射救回 ──────────────────────────────
rescued = 0
blank_rows = 0
for r in rows:
    if val(r.get("calib_action")) != "blank":
        continue
    blank_rows += 1
    fam = val(r.get("Family"))
    pg = val(r.get("calib_prev_genus"))
    if not fam or not pg:
        continue
    cands = msl_gen2fam.get(pg.lower(), set())
    if cands and fam.lower() not in cands:
        # 若把行内老科名换成该属在 MSL41 的定型科, 是否就相容了
        rescued += 1
print(f"\n  blank {blank_rows} 行: 属在 MSL41 有定型科且与行内科名不符 -> 可考虑改写科名后保留属: {rescued} 行")

# ── 边界三 B: 有属无科的 947 行能否补科 ────────────────────────────────
fill_ok = 0
fill_conf = 0
nofam = 0
for r in rows:
    fam = val(r.get("Family"))
    g = val(r.get("Genus"))
    if fam or not g:
        continue
    nofam += 1
    ncbi, vmr = ref.get(g.lower(), ("", ""))
    if ncbi and vmr and ncbi == vmr:
        fill_ok += 1
    elif ncbi or vmr:
        fill_conf += 1
print(f"  有属无科 {nofam} 行: 两参照科名一致可直接补 -> {fill_ok} 行 | 仅一侧或冲突 -> {fill_conf} 行")

# ── 亲代表对照: 老科名规模 ─────────────────────────────────────────────
o = list(csv.DictReader(open(ORIG, errors="replace"), delimiter="\t"))
o_not = 0
for r in o:
    fam = val(r.get("Family"))
    if fam and fam.lower() not in msl_fam:
        o_not += 1
print(f"\n  亲代表 {len(o)} 行中 Family 非 MSL41 科名: {o_not} 行")
