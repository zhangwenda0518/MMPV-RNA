"""排查(v2): 基于权威宿主概率表的"科-属不一致"。

判定逻辑 (客观、不自我循环):
  1. 科级: 从 family_host_probability.tsv 取每个科的 Predicted_Host
     - 植物科 := Predicted_Host == Plant
  2. 属级: 从 genus_host_probability.tsv 取每个属的 Predicted_Host
  3. 交叉: 对每个"植物科", 列出其下所有属, 标出 Predicted_Host != Plant 的属
     -> 这就是"科是植物、属不是植物"的漏网嫌疑
  4. 叠加黑名单标记: 已在黑名单的标 [已处理], 未在黑名单的标 [待核]
"""
import ast
import csv
from collections import defaultdict

CA = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"

# ── 黑名单 ──
RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
src = open(f"{RUN}/run_host_prediction.py").read()
tree = ast.parse(src)
set_assigns = [n for n in tree.body if isinstance(n, ast.Assign)
               and any(isinstance(t, ast.Name) and
                       t.id in ("NON_PLANT_GENERA", "NON_PLANT_FAMILIES_FALLBACK",
                                "TRUSTED_LEVELS") for t in n.targets)]
funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef)
         and n.name in ("is_trusted_level", "is_blacklisted")]
ns = {}
exec(compile(ast.Module(body=set_assigns + funcs, type_ignores=[]), "<m>", "exec"), ns)
is_blacklisted = ns["is_blacklisted"]
NON_PLANT_GENERA = ns.get("NON_PLANT_GENERA", set())
NON_PLANT_FAMILIES_FALLBACK = ns.get("NON_PLANT_FAMILIES_FALLBACK", set())


def load_prob(path, key):
    """{name: (Predicted_Host, P_Max, Confidence_Level, Total_Records, PPlant)}"""
    d = {}
    with open(path, errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            n = (row.get(key, "") or "").strip()
            if not n:
                continue
            try:
                pplant = float(row.get("P(Plant)", 0) or 0)
            except Exception:
                pplant = 0.0
            d[n] = {
                "host": (row.get("Predicted_Host", "") or "").strip(),
                "pmax": (row.get("P_Max", "") or "").strip(),
                "lvl": (row.get("Confidence_Level", "") or "").strip(),
                "n": (row.get("Total_Records", "") or "").strip(),
                "pplant": pplant,
            }
    return d


fam_p = load_prob(f"{CA}/family_host_probability.tsv", "Family")
gen_p = load_prob(f"{CA}/genus_host_probability.tsv", "Genus")
print(f"科概率表: {len(fam_p)} 条, 属概率表: {len(gen_p)} 条")

plant_fams = sorted([k for k, v in fam_p.items() if v["host"] == "Plant"])
print(f"\n权威判定为 Plant 的科 ({len(plant_fams)}):")
for k in plant_fams:
    print(f"  {k:<26} P={fam_p[k]['pmax']:<8} lvl={fam_p[k]['lvl']:<12} n={fam_p[k]['n']}")

# ── 科-属交叉 ──
# 建立 属->科 映射: 从 vmr.tsv / taxa.txt
gen_to_fam = {}
import os
for cand in ["/home/zhangwenda/database/virus-db/acvirus_db/vmr.tsv",
             "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"]:
    if not os.path.exists(cand):
        continue
    with open(cand, errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        try:
            gi = hdr.index("Genus")
            fi = hdr.index("Family")
        except ValueError:
            continue
        for row in rd:
            if len(row) <= max(gi, fi):
                continue
            g = row[gi].strip()
            fm = row[fi].strip()
            if g and fm and g not in gen_to_fam:
                gen_to_fam[g] = fm
    if gen_to_fam:
        break
print(f"\n属->科 映射: {len(gen_to_fam)} 条")

print("\n" + "=" * 78)
print("=== 植物科下, Predicted_Host != Plant 的属 ===")
print("=" * 78)
print(f"{'科':<22} {'属':<24} {'Pred_Host':<14} {'P(plant)':>8} {'黑名单':<8}")
print("-" * 78)
rows = []
for fm in plant_fams:
    for g, v in gen_p.items():
        if gen_to_fam.get(g) != fm:
            continue
        if v["host"] == "Plant":
            continue
        bl = "已处理" if is_blacklisted(fm, g, "Genus") else "待核"
        rows.append((fm, g, v["host"], v["pplant"], bl))

rows.sort(key=lambda x: (x[4] != "待核", -x[3]))
n_todo = sum(1 for r in rows if r[4] == "待核")
n_done = len(rows) - n_todo
for fm, g, h, pp, bl in rows:
    print(f"{fm:<22} {g:<24} {h:<14} {pp:>8.3f} {bl:<8}")
print(f"\n合计 {len(rows)} 条 (待核 {n_todo}, 已处理 {n_done})")
