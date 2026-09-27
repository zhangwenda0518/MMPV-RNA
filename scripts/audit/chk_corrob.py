#!/usr/bin/env python3
"""量: 本行双名物种名佐证当前属的行数 —— 决定新函数是否把双名并入佐证集

R 侧判据 (harmonize_genus_species 同一套):
  Species 含空格 -> 首词; 有效需 virus$ 结尾 且非 viridae$/virinae$ 且非子阶元后缀
"""
import csv
import importlib.util
import re
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
            "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated")
PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")

TOOLS = ["ACVirus", "CAT", "diamond_lca", "genomad", "metabuli", "mmseqs", "VITAP"]
BIAS = {"ACVirus": 1.2, "VITAP": 1.1, "mmseqs": 1.0, "metabuli": 1.0,
        "CAT": 0.9, "genomad": 0.9, "diamond_lca": 0.8}
GW = {t: 64.0 * BIAS.get(t, 0.8) for t in TOOLS}
SUBRANK = ["viricotina", "viricetidae", "virineae", "virinae"]

spec = importlib.util.spec_from_file_location("ann", PIPE / "utils" / "annotate_nucleic_acid.py")
ann = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ann)
VMR_IDX = ann.load_vmr(VMR)

per_tool = defaultdict(dict)
for t in TOOLS:
    with open(BASE / ("standardized_%s.tsv" % t), newline="", encoding="utf-8",
              errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            per_tool[r.get("contig_id")][t] = r
with open(BASE / "final_integrated_classification.tsv", newline="", encoding="utf-8",
          errors="replace") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

B = ann.is_blank
V = lambda v: not B(v)


def sp_genus(sp):
    if not sp or " " not in sp:
        return None
    g = sp.split(" ")[0]
    if not re.search(r"virus$", g, re.IGNORECASE):
        return None
    if re.search(r"viridae$|virinae$", g, re.IGNORECASE):
        return None
    if any(re.search(re.escape(s) + "$", g, re.IGNORECASE) for s in SUBRANK):
        return None
    return g


fix_map, blank_ids = {}, set()
for r in rows:
    cid, F, G = r["contig_id"], r.get("Family", ""), r.get("Genus", "")
    if not V(F) or not V(G):
        continue
    tl = per_tool.get(cid, {})
    S = [t for t in tl if (tl[t].get("Family") or "") == F]
    if not S:
        continue
    cand = Counter()
    for t in S:
        g = tl[t].get("Genus")
        if V(g):
            cand[g] += GW[t]
    if G in cand:
        continue
    if cand:
        fix_map[cid] = cand.most_common(1)[0][0]
    else:
        blank_ids.add(cid)

print("校准集合: 替换 %d, 置空 %d" % (len(fix_map), len(blank_ids)))
print()

n_b = n_b_eq = n_f = n_f_eq = 0
ex_eq = []
for r in rows:
    cid, G = r["contig_id"], r.get("Genus", "")
    sg = sp_genus(r.get("Species", ""))
    hit = bool(sg) and sg.lower() == G.lower()
    if cid in blank_ids:
        n_b += 1
        if hit:
            n_b_eq += 1
            if len(ex_eq) < 6:
                ex_eq.append((cid, r.get("Family"), G, r.get("Species"), "置空"))
    elif cid in fix_map:
        n_f += 1
        if hit:
            n_f_eq += 1
            if len(ex_eq) < 6:
                ex_eq.append((cid, r.get("Family"), G, r.get("Species"), "替换"))

print("=== 双名佐证（本行 Species 首词 == 当前 Genus） ===")
print("  被置空行中双名佐证当前属 : %4d / %d" % (n_b_eq, n_b))
print("  被替换行中双名佐证当前属 : %4d / %d" % (n_f_eq, n_f))
print("  合计需豁免置空/替换      : %4d" % (n_b_eq + n_f_eq))
print()
if ex_eq:
    print("  样例:")
    for cid, F, G, sp, act in ex_eq:
        print("    [%s] Family=%s Genus=%s Species=%s" % (act, F, G, sp))

# VMR 里 Guapo partitivirus 到底有没有属
print()
print("=== VMR Species 键核查 ===")
lvl = VMR_IDX["lv"]
for nm in ["Guapo partitivirus", "Guapo_partitivirus"]:
    print("  Species 键 '%s' 在 VMR: %s" % (nm, nm in lvl.get("Species", {})))
gs = lvl.get("Genus", {})
print("  Genus 键 'Biavirus' 在 VMR: %s  -> %s" % ("Biavirus" in gs,
                                                  dict(gs["Biavirus"]) if "Biavirus" in gs else ""))
hits = [k for k in lvl.get("Species", {}) if k.lower().startswith("guapo")]
print("  VMR 中 guapo* 物种: %s" % hits[:8])
