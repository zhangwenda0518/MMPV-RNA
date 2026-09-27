"""黑名单第二问: 开火之后行到底去哪了 (RVH/PB2 复活路径).
数据: ensemble_host_summary.tsv (含 Host_ICTV / pred|L1 / evidence / Host lineage / Final_Host / Decision_Method)
     + C9 classification_result.tsv (取 Determination_Level)
"""
import ast
import csv
import os
from collections import Counter, defaultdict

RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
C7 = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
TREES = {
    "goji-Lycium": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}

src = open(RUN, encoding="utf-8", errors="replace").read()
t = ast.parse(src)
sets = [n for n in t.body if isinstance(n, ast.Assign)
        and any(isinstance(x, ast.Name) and x.id in (
            "NON_PLANT_GENERA", "NON_PLANT_FAMILIES_FALLBACK", "TRUSTED_LEVELS")
            for x in n.targets)]
fns = [n for n in t.body if isinstance(n, ast.FunctionDef)
       and n.name in ("is_trusted_level", "is_blacklisted")]
ns = {}
exec(compile(ast.Module(body=sets + fns, type_ignores=[]), "<m>", "exec"), ns)
is_blacklisted = ns["is_blacklisted"]

RVH_MAP = {'viridiplantae': 'Plant', 'fungi': 'Fungi', 'chordata': 'Animal',
           'invertebrate': 'Animal', 'metazoa': 'Animal', 'animal': 'Animal',
           'bacteria': 'Bacteria'}


def norm_c9(h):
    if h is None or str(h).strip() in ('', 'Unknown', 'None', 'nan'):
        return 'Unknown'
    h = str(h).strip()
    if h in ['Insecta', 'Arachnida', 'Aves', 'Human', 'Animal_other']:
        return 'Animal'
    if h == 'Oomycetes':
        return 'Protist'
    return h


def parse_rvh(row):
    h = str(row.get('pred|L1', 'Unknown')).strip().lower()
    ev = str(row.get('evidence', '')).strip().lower()
    if h == 'unknown' or ev == 'unclassified':
        return 'Unknown'
    return RVH_MAP.get(h, 'Unknown')


def parse_pb2(row):
    lin = (str(row.get('Host_NCBI_lineage', '')) + "|" +
           str(row.get('Host_GTDB_lineage', '')) + "|" +
           str(row.get('Host', ''))).lower()
    if 'bacteria' in lin: return 'Bacteria'
    if 'archaea' in lin: return 'Archaea'
    if 'streptophyta' in lin or 'viridiplantae' in lin or 'plant' in lin: return 'Plant'
    if 'fungi' in lin: return 'Fungi'
    if any(k in lin for k in ['metazoa', 'animal', 'chordata', 'arthropoda', 'insecta']):
        return 'Animal'
    return 'Unknown'


# C7 原生/藻名单
def load_c7(fn, namecol):
    d = defaultdict(set)
    with open(f"{C7}/{fn}", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[(r.get("Predicted_Host") or "").strip()].add((r.get(namecol) or "").strip())
    return d

fb = load_c7("family_host_probability.tsv", "Family")
gb = load_c7("genus_host_probability.tsv", "Genus")
pa_fam = fb.get("Protist", set()) | fb.get("Algae", set())
pa_gen = gb.get("Protist", set()) | gb.get("Algae", set())

for label, root in TREES.items():
    print("\n" + "=" * 78)
    print(f"[{label}] {root}")
    print("=" * 78)
    lvl = {}
    f9 = f"{root}/06_HostPrediction/C9_ICTV_result/classification_result.tsv"
    with open(f9, errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            lvl[r["contig_id"]] = (r.get("Determination_Level") or "").strip()

    p = f"{root}/06_HostPrediction/ensemble_host_summary.tsv"
    total = final_plant = fire = 0
    fire_final = Counter()
    fire_plant_method = Counter()
    fire_plant_evid = Counter()
    rescuable = Counter()
    plant_method = Counter()
    plant_lvl = Counter()
    plant_immune_pa = 0
    plant_notfire_pa_fam = Counter()
    fire_unknown = 0
    with open(p, errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            cid = row["contig_id"]
            lv = lvl.get(cid, "")
            ictv = norm_c9(row.get("Host_ICTV"))
            fam = (row.get("Family") or "").strip()
            gen = (row.get("Genus") or "").strip()
            fh = (row.get("Final_Host") or "").strip()
            meth = (row.get("Decision_Method") or "").strip()
            fired = (ictv == 'Plant' and is_blacklisted(fam, gen, lv))
            if fired:
                fire += 1
                fire_final[fh] += 1
                if fh == 'Plant':
                    fire_plant_method[meth] += 1
                    fire_plant_evid[(str(row.get('pred|L1', '')).strip(),
                                     str(row.get('evidence', '')).strip())] += 1
            if fh == 'Plant':
                final_plant += 1
                plant_method[meth] += 1
                plant_lvl[lv.split("(")[0]] += 1
                if fam in pa_fam or gen in pa_gen:
                    if not fired:
                        plant_immune_pa += 1
                        plant_notfire_pa_fam[(fam, gen, lv.split("(")[0])] += 1
    print(f"表 {total} 行 | Final_Host=Plant {final_plant} 行 | 黑名单开火 {fire} 行")
    print(f"\n开火行的最终归属: {dict(fire_final.most_common())}")
    print(f"\n开火但仍判 Plant 的 {sum(fire_plant_method.values())} 行 -> 决策路径:")
    for k, v in fire_plant_method.most_common():
        print(f"    {k:<22} {v}")
    print(f"    RVH 证据分布: {dict(fire_plant_evid.most_common(8))}")
    print(f"\nFinal_Host=Plant 的决策路径: {dict(plant_method.most_common())}")
    print(f"Final_Host=Plant 的判定层级: {dict(plant_lvl.most_common())}")
    print(f"\n* 未开火却带原生/藻科或属名的 Plant 行: {plant_immune_pa}")
    for (fm, g, l), n in plant_notfire_pa_fam.most_common(12):
        print(f"    {fm:<24} {g:<26} {l:<10} {n}")
