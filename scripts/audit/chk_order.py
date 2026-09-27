#!/usr/bin/env python3
"""量清 科-属校准 与 harmonize_genus_species 的先后顺序风险

harmonize_genus_species 判据（照 R 源码复刻）：
  Species 含空格 -> 首词 genus_from_sp
  有效: 首词以 virus$ 结尾(忽略大小写) 且 不以 viridae$/virinae$ 结尾 且 不以子阶元后缀结尾
  且 要求 Genus 非 NA 且与首词不同 -> 才覆写（只覆写，不填空）
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
idx = ann.load_vmr(VMR)

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


def species_genus(sp):
    """复刻 harmonize_genus_species 的首词提取与有效性判定"""
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


# ---- 重算校准集合 ----
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

# ---- 关键交叉：置空行里有多少能靠双名救回 ----
n_blank_species = sum(1 for r in rows if r["contig_id"] in blank_ids
                      and species_genus(r.get("Species", "")))
n_fix_species = sum(1 for r in rows if r["contig_id"] in fix_map
                    and species_genus(r.get("Species", "")))
n_species_fire = sum(1 for r in rows
                     if species_genus(r.get("Species", ""))
                     and V(r.get("Genus", ""))
                     and species_genus(r.get("Species", "")).lower() != r["Genus"].lower())
n_guard_block = sum(1 for r in rows if B(r.get("Genus", ""))
                    and species_genus(r.get("Species", "")))
print()
print("=== 与 harmonize_genus_species 的交叉 ===")
print("  被置空且 Species 是合法双名            : %5d" % n_blank_species)
print("  被替换且 Species 是合法双名            : %5d" % n_fix_species)
print("  当前会触发 species->genus 覆写的行      : %5d" % n_species_fire)
print("  当前 Genus 为空但有合法双名(被守卫挡住)  : %5d" % n_guard_block)

# 顺序模拟：先 species 再 family  vs  先 family 再 species
print()
print("=== 两种顺序的终态差异（只看属这一列） ===")


def apply_species(g, sp):
    sg = species_genus(sp)
    if sg and V(g) and sg.lower() != g.lower():
        return sg
    return g


def apply_family(cid, g, blank=True):
    if cid in fix_map:
        return fix_map[cid]
    if blank and cid in blank_ids:
        return "NA"
    return g


o1 = 0   # 先 family 后 species
o2 = 0   # 先 species 后 family
diff = 0
for r in rows:
    cid, g0 = r["contig_id"], r.get("Genus", "")
    sp = r.get("Species", "")
    a = apply_species(apply_family(cid, g0), sp)          # family -> species
    b = apply_family(cid, apply_species(g0, sp))          # species -> family
    if a != b:
        diff += 1
    if a != g0:
        o1 += 1
    if b != g0:
        o2 += 1
print("  先 family 后 species : 属变化 %d 行" % o1)
print("  先 species 后 family : 属变化 %d 行" % o2)
print("  两种顺序终态不同的行   : %d" % diff)

# ---- VMR 里这条物种到底有没有属 ----
print()
print("=== VMR 核查: Guapo partitivirus ===")
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    cols = rd.fieldnames
    for r in rd:
        sp = r.get("Species") or r.get("Virus name(s)") or ""
        if "Guapo partitivirus" in sp:
            print("  命中: %s" % {k: r[k] for k in cols
                                 if k in ("Species", "Genus", "Family", "Genome", "Virus name(s)")})
print("  VMR 列名:", cols[:14])
