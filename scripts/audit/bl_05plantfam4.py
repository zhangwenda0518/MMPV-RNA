"""修正版: 用已知病毒科全集做植物科判定; 7 个非 MSL41 科逐行拆开看属."""
import csv
import os
import re
from collections import Counter, defaultdict

MSL = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"
REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")
PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
        "/05_Taxonomy/Votus.integrated")
CAL = f"{BASE}/calibration_20260914/final_integrated_classification.calibrated.tsv"

FAMS7 = ["Hepaciviridae", "Pestiviridae", "Zimmerviridae", "Ambiguiviridae",
         "Autographiviridae", "Ourmiaviridae", "Parahypoviridae"]

# 已知病毒科/属全集
known_fam, known_genus = set(), set()
msl_fam = set()
with open(MSL, errors="replace") as f:
    for r in csv.DictReader(f):
        fam = (r.get("Family") or "").strip().lower()
        g = (r.get("Genus") or "").strip().lower()
        if fam:
            known_fam.add(fam)
            msl_fam.add(fam)
        if g:
            known_genus.add(g)
with open(REF, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        g = (r.get("Genus") or "").strip().lower()
        if g:
            known_genus.add(g)
        for c in ("NCBI_Family", "VMR_Family"):
            v = (r.get(c) or "").strip().lower()
            if v:
                known_fam.add(v)
print(f"已知科全集 {len(known_fam)} | MSL41 科 {len(msl_fam)} | 已知属 {len(known_genus)}")

# Plant.tsv: 植物病毒库认的科与属
plant_fam, plant_genus, plant_virus_pairs = set(), set(), 0
with open(PLANT, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        plant_virus_pairs += 1
        toks = [t.strip().lower() for t in re.split(r"[;,]", r.get("Virus_lineage") or "") if t.strip()]
        for t in toks:
            if t in known_fam:
                plant_fam.add(t)
            if t in known_genus:
                plant_genus.add(t)
print(f"植物病毒库认的科 {len(plant_fam)} 个 | 属 {len(plant_genus)} 个 | 宿主对 {plant_virus_pairs}")

print("\n=== 7 个非 MSL41 科的植物归属 (三重口径) ===")
print(f"{'科名':<22}{'在MSL41':>8}{'植物库认':>10}{'宿主黑名单':>12}")
RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
import ast
src = open(RUN, encoding="utf-8", errors="replace").read()
t = ast.parse(src)
keep = [n for n in t.body if isinstance(n, (ast.Import, ast.ImportFrom, ast.Assign, ast.FunctionDef, ast.AnnAssign))]
ns = {"__file__": RUN}
exec(compile(ast.Module(body=keep, type_ignores=[]), "<rhp>", "exec"), ns)
BF = set(ns.get("NON_PLANT_FAMILIES_FALLBACK", []) or [])
for k in FAMS7:
    print(f"  {k:<20}{'是' if k.lower() in msl_fam else '否':>8}"
          f"{'是' if k.lower() in plant_fam else '否':>10}{'在' if k in BF else '不在':>12}")
print(f"  (参考: Botourmiaviridae 在 MSL41 {'是' if 'botourmiaviridae' in msl_fam else '否'}, "
      f"植物库 {'是' if 'botourmiaviridae' in plant_fam else '否'})")


def val(x):
    x = (x or "").strip().strip('"')
    return None if x in ("", "NA", "N/A", "-") else x


rows = list(csv.DictReader(open(CAL, errors="replace"), delimiter="\t"))

# ── 7 个科的行: 属种是什么, 属是不是植物属 ───────────────────────────
print("\n=== 7 个非 MSL41 科的 180 行: 科 x 属 拆解 ===")
agg = Counter()
for r in rows:
    fam = val(r.get("Family"))
    if fam not in FAMS7:
        continue
    g = val(r.get("Genus")) or "(空)"
    pg = val(r.get("calib_prev_genus"))
    agg[(fam, g, "植物属" if g.lower() in plant_genus else "非植物属")] += 1
for (fam, g, tag), v in sorted(agg.items(), key=lambda x: (-x[1], x[0][0])):
    print(f"  {fam:<20} 属={g:<26} {tag:<8} {v} 行")

# ── 修正后的: 科是病毒科但不属植物库, 而属是植物属 ───────────────────
mis = Counter()
for r in rows:
    fam = (val(r.get("Family")) or "").lower()
    g = (val(r.get("Genus")) or "").lower()
    if fam and g and fam in known_fam and fam not in plant_fam and g in plant_genus:
        mis[(val(r.get("Family")), val(r.get("Genus"))) ] += 1
print(f"\n=== 科是已知病毒科但不被植物库收录, 而属是植物病毒属: {sum(mis.values())} 行 ===")
for (a, b), v in mis.most_common(20):
    print(f"  科={a:<22} 属={b:<26} {v} 行")

# ── 251 vmr_only 漂移行: 属是不是植物属 ─────────────────────────────
d = Counter()
for r in rows:
    if val(r.get("calib_flag_single_ref")) != "vmr_only":
        continue
    g = val(r.get("Genus"))
    d["属是植物病毒属" if (g or "").lower() in plant_genus else "属非植物/未知"] += 1
print("\n=== 251 vmr_only 漂移行: 属的植物归属 ===")
for k, v in d.items():
    print(f"  {k:<22} {v} 行")
