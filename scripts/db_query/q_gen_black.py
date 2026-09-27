import csv
import collections
import glob
import os

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

auth_gen = set()
auth_fam = set()
with open(PLANT, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        for p in (row.get("Virus_lineage", "") or "").split(";"):
            p = p.strip()
            if not p:
                continue
            if p.endswith("viridae"):
                auth_fam.add(p)
            elif p.endswith("virus"):
                auth_gen.add(p)

base_dir = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"
gen_stat = collections.Counter()
gen_fam = {}
gen_det = {}
for tsv in glob.glob(os.path.join(base_dir, "*.classified.tsv")):
    with open(tsv, errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if (row.get("Predicted_Host", "") or "").strip() != "Plant":
                continue
            g = (row.get("Genus", "") or "").strip()
            fm = (row.get("Family", "") or "").strip()
            d = (row.get("Determination_Level", "") or "").split("(")[0].strip()
            gen_stat[g] += 1
            gen_fam.setdefault(g, collections.Counter())[fm] += 1
            gen_det.setdefault(g, collections.Counter())[d] += 1

# 属不在权威库 & 其所属科也不在整科黑名单 (即跨宿主科内的属级误判)
print("=== C9 判 Plant 但权威库无此属 (跨宿主科的属级误判候选) ===")
print(f"{'属':<26} {'条数':<7} {'主科':<22} {'科级判定':<8}")
rows = []
for g, n in gen_stat.items():
    if not g or g in ("NA", "nan") or g in auth_gen:
        continue
    top_fam = gen_fam[g].most_common(1)[0][0] if gen_fam[g] else ""
    dets = gen_det.get(g, {})
    fam_lvl = dets.get("Family", 0) + dets.get("Order", 0) + dets.get("None", 0)
    rows.append((g, n, top_fam, fam_lvl))
rows.sort(key=lambda x: -x[1])
for g, n, fm, fl in rows:
    print(f"{g:<26} {n:<7} {fm:<22} {fl:<8}")
print(f"\n共 {len(rows)} 个属, 合计 {sum(r[1] for r in rows)} 条")
