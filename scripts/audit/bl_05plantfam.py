"""这几个科是不是植物病毒科: 用自家 Plant.tsv 白名单 + 宿主黑名单 + 下游植物表三重核对."""
import csv
import os
from collections import Counter

FAMS7 = ["Hepaciviridae", "Pestiviridae", "Zimmerviridae", "Ambiguiviridae",
         "Autographiviridae", "Ourmiaviridae", "Parahypoviridae"]
# 版本漂移/置空行涉及的科
DRIFT = ["Mimiviridae", "Phycodnaviridae", "Marseilleviridae", "Pithoviridae", "Iridoviridae",
         "Schizomimiviridae", "Partitiviridae", "Retroviridae", "Astroviridae", "Orpheoviridae",
         "Epsomviridae", "Mycoalphaviridae", "Hydriviridae", "Botourmiaviridae", "Flaviviridae"]
ALL = FAMS7 + [f for f in DRIFT if f not in FAMS7]

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
print("Plant.tsv 存在:", os.path.exists(PLANT), os.path.getsize(PLANT) if os.path.exists(PLANT) else "")


def find_col(cols, *keys):
    for c in cols:
        lc = c.lower()
        if all(k in lc for k in keys):
            return c
    return None


def low(x):
    return (x or "").strip().lower()


# ── 1. Plant.tsv: 家里自己的植物病毒白名单 ────────────────────────────
plant_fam = Counter()
plant_gen = Counter()
plant_sp = set()
with open(PLANT, errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    cols = rd.fieldnames
    print("Plant.tsv 列:", cols)
    fc = find_col(cols, "family")
    gc = find_col(cols, "genus")
    for r in rd:
        if fc:
            v = (r.get(fc) or "").strip()
            if v:
                plant_fam[v] += 1
        if gc:
            v = (r.get(gc) or "").strip()
            if v:
                plant_gen[v] += 1

print("\n=== 这几个科在自家 Plant.tsv (201,488 行植物-病毒对) 里的出现情况 ===")
for k in ALL:
    n = plant_fam.get(k, 0)
    tag = "植物病毒库里' 有" if n else "植物病毒库里 无"
    print(f"  {k:<22} {n:>6} 行   {tag.replace(chr(39), '')}")

print("\nPlant.tsv 里出现最多的科 (前 15):")
for k, v in plant_fam.most_common(15):
    print(f"  {k:<26} {v}")

# ── 2. 宿主黑名单里有没有登记 ────────────────────────────────────────
RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
import ast
src = open(RUN, encoding="utf-8", errors="replace").read()
t = ast.parse(src)
keep = [n for n in t.body if isinstance(n, (ast.Import, ast.ImportFrom, ast.Assign, ast.FunctionDef, ast.AnnAssign))]
ns = {"__file__": RUN}
exec(compile(ast.Module(body=keep, type_ignores=[]), "<rhp>", "exec"), ns)
BF = set(ns.get("NON_PLANT_FAMILIES_FALLBACK", []) or [])
print(f"\n=== 宿主黑名单 NON_PLANT_FAMILIES_FALLBACK ({len(BF)} 条) 覆盖情况 ===")
for k in ALL:
    print(f"  {k:<22} {'在名单内' if k in BF else '不在名单内'}")

# ── 3. 下游植物表里实际有多少行落到这几个科 ──────────────────────────
DOWN = [
    ("枸杞 Lycium barbarum",
     "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/10_Reports/All_plant.viruses_info.tsv"),
    ("OneKP",
     "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"),
]
for label, p in DOWN:
    print(f"\n=== 下游植物表 {label}: {p}")
    if not os.path.exists(p):
        print("  不存在")
        continue
    with open(p, errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t")
        cols = rd.fieldnames
        rows = list(rd)
    print(f"  列: {cols}")
    print(f"  行数: {len(rows)}")
    hit = Counter()
    for r in rows:
        blob = " ".join(str(v) for v in r.values())
        for k in ALL:
            if k in blob:
                hit[k] += 1
    print("  命中的目标科:", dict(hit) if hit else "无")
