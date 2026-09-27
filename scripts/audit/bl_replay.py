"""真跑验证: 用服务器现版 run_host_prediction.py 的级联函数重放两棵树的真实表格.
输入 = ensemble_host_summary.tsv (含 Host_ICTV/pred|L1/evidence/Host lineage/Class/Family/Genus)
     + C9 classification_result.tsv (取 Determination_Level)
输出 = 现版代码下的 Final_Host, 与树里历史 Final_Host 对比
"""
import ast
import csv
import os
from collections import Counter
from pathlib import Path

RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
K = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus"
G = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"

src = open(RUN, encoding="utf-8", errors="replace").read()
t = ast.parse(src)
keep = [n for n in t.body if isinstance(
    n, (ast.Import, ast.ImportFrom, ast.Assign, ast.FunctionDef, ast.AnnAssign))]
mod = ast.Module(body=keep, type_ignores=[])
ns = {"__file__": RUN}
exec(compile(mod, "<rhp>", "exec"), ns)
casc = ns["decision_tree_cascade"]
is_bl = ns["is_blacklisted"]
print("[0] 现版函数已加载:", casc.__name__, "/", len(ns.get("NON_PLANT_GENERA", [])), "属 /",
      len(ns.get("NON_PLANT_FAMILIES_FALLBACK", [])), "科 / 豁免", sorted(ns["TRUSTED_LEVELS"]))

# C7 原生/藻名单 (用于识别残留行)
import csv as _c7
PA_F, PA_G = set(), set()
for _fn, _col, _s in (("family_host_probability.tsv", "Family", PA_F),
                      ("genus_host_probability.tsv", "Genus", PA_G)):
    with open(f"/home/zhangwenda/MMPV-RNA/database/cross_analysis/{_fn}", errors="replace") as _f:
        for _r in _c7.DictReader(_f, delimiter="\t"):
            if (_r.get("Predicted_Host") or "").strip() in ("Protist", "Algae"):
                _s.add((_r.get(_col) or "").strip())
print(f"    C7 原生/藻: 科 {len(PA_F)} / 属 {len(PA_G)}")

for label, root in [("goji-Lycium", G), ("onekp", K)]:
    print("\n" + "=" * 78)
    print(f"[{label}]")
    print("=" * 78)
    lvl = {}
    with open(f"{root}/06_HostPrediction/C9_ICTV_result/classification_result.tsv", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            lvl[r["contig_id"]] = (r.get("Determination_Level") or "").strip()

    shift = Counter()
    fire_new = Counter()
    plant_new = Counter()
    resid_old = 0
    resid_new = Counter()
    n = 0
    n_plant_old = 0
    changed = 0
    with open(f"{root}/06_HostPrediction/ensemble_host_summary.tsv", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            n += 1
            cid = row["contig_id"]
            row["Determination_Level"] = lvl.get(cid, "")
            old = (row.get("Final_Host") or "").strip()
            new, meth = casc(row)
            if old != new:
                changed += 1
                shift[(old, new, meth)] += 1
            if old == "Plant":
                n_plant_old += 1
                plant_new[(new, meth)] += 1
            ictv = ns["normalize_c9"](row.get("Host_ICTV"))
            if ictv == "Plant" and is_bl(row.get("Family"), row.get("Genus"), row.get("Determination_Level")):
                fire_new[(old, new, meth)] += 1
            fam = (row.get("Family") or "").strip()
            gen = (row.get("Genus") or "").strip()
            if old == "Plant" and (fam in PA_F or gen in PA_G):
                resid_old += 1
                resid_new[(new, meth, row.get("Determination_Level", "").split("(")[0])] += 1
    print(f"表 {n} 行 | 历史 Final_Host=Plant {n_plant_old} 行 | 现版代码重放后变动 {changed} 行")
    print(f"\n变动明细 (旧 -> 新, 路径):")
    for k, v in shift.most_common(20):
        print(f"    {str(k):<58} {v}")
    print(f"\n历史 Plant 行重放后的归属:")
    for k, v in plant_new.most_common(12):
        print(f"    {str(k):<58} {v}")
    print(f"\n当前名单会开火的行 (旧->新):")
    for k, v in fire_new.most_common(12):
        print(f"    {str(k):<58} {v}")
    print(f"\n带原生/藻科或属名的历史 Plant 行: {resid_old} 行, 重放后归属:")
    for k, v in resid_new.most_common(12):
        print(f"    {str(k):<58} {v}")
