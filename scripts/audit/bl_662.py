"""还原被黑名单删掉的 662 条: 是什么、为什么被删、现行名单能拦到几分.
数据源:
  bak: onekp-virus/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv.bak_bl662_20260911
  当前: 同路径 All_plant.viruses_info.tsv
  属性: 06_HostPrediction/ensemble_host_summary.tsv + C9_ICTV_result/classification_result.tsv
"""
import ast
import csv
import os
from collections import Counter

RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
K = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus"
C7 = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"

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
NG = set(ns["NON_PLANT_GENERA"])
NF = set(ns["NON_PLANT_FAMILIES_FALLBACK"])


def load_c7(fn, col):
    d = {}
    with open(f"{C7}/{fn}", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[(r[col] or "").strip()] = (r.get("Predicted_Host") or "").strip()
    return d


fam_host = load_c7("family_host_probability.tsv", "Family")
gen_host = load_c7("genus_host_probability.tsv", "Genus")
PA = ("Protist", "Algae")


def ids(path):
    s = set()
    with open(path, errors="replace") as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            c = l.split("\t", 1)[0].strip()
            if c:
                s.add(c)
    return s


bak = f"{K}/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv.bak_bl662_20260911"
cur = f"{K}/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"
B, C = ids(bak), ids(cur)
gone = B - C
add = C - B
print("=" * 78)
print(f"[1] 备份 {len(B)} 条 / 当前 {len(C)} 条 / 备份有当前无 {len(gone)} 条 / 当前有备份无 {len(add)} 条")
print("=" * 78)

# 属性表
attr = {}
f9 = f"{K}/06_HostPrediction/C9_ICTV_result/classification_result.tsv"
with open(f9, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        attr[r["contig_id"]] = ((r.get("Family") or "").strip(),
                               (r.get("Genus") or "").strip(),
                               (r.get("Determination_Level") or "").strip(),
                               (r.get("Predicted_Host") or "").strip())
summ = {}
fS = f"{K}/06_HostPrediction/ensemble_host_summary.tsv"
with open(fS, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        summ[r["contig_id"]] = ((r.get("Host_ICTV") or "").strip(),
                               (r.get("Final_Host") or "").strip(),
                               (r.get("Decision_Method") or "").strip())

print("\n[2] 被删 662 条属性 (来自 C9 表):")
hit_c9 = Counter()
fam_c = Counter()
gen_c = Counter()
lvl_c = Counter()
host_c = Counter()
route = Counter()
for cid in sorted(gone):
    fm, g, lv, ph = attr.get(cid, ("<不在C9表>", "", "", ""))
    fam_c[fm] += 1
    gen_c[g] += 1
    lvl_c[lv.split("(")[0]] += 1
    host_c[ph] += 1
    hit_c9[ph] += 1
    if fm in NF and (not g or g in ("NA", "nan")):
        route["科兜底名单命中"] += 1
    elif g in NG:
        route["属黑名单命中"] += 1
    else:
        route["两条名单都不命中"] += 1

print(f"    C9 判定宿主: {dict(host_c.most_common())}")
print(f"    判定层级: {dict(lvl_c.most_common())}")
print(f"    删因路径: {dict(route.most_common())}")
print(f"\n    被删 top15 科:")
for k, v in fam_c.most_common(15):
    tag = ""
    if k in NF:
        tag += " [在科兜底名单]"
    if fam_host.get(k) in PA:
        tag += " [C7判原生/藻]"
    print(f"      {k:<26} {v:>4}{tag}")
print(f"\n    被删 top15 属:")
for k, v in gen_c.most_common(15):
    tag = ""
    if k in NG:
        tag += " [在属黑名单]"
    if gen_host.get(k) in PA:
        tag += " [C7判原生/藻]"
    print(f"      {k:<26} {v:>4}{tag}")

print("\n[3] 被删行的 summary 决策路径:")
m_c = Counter()
fh_c = Counter()
miss = 0
for cid in gone:
    if cid not in summ:
        miss += 1
        continue
    hi, fh, me = summ[cid]
    m_c[me] += 1
    fh_c[fh] += 1
print(f"    在 summary 中找不到的: {miss}")
print(f"    Final_Host: {dict(fh_c.most_common())}")
print(f"    Decision_Method: {dict(m_c.most_common())}")

print("\n[4] 当前仍在 Plant 名单里的残留 (对照):")
pf = f"{K}/06_HostPrediction/host_classified_fasta/Plant.classified.fasta"
cur_plant = set()
with open(pf, errors="replace") as f:
    for l in f:
        if l.startswith(">"):
            cur_plant.add(l[1:].split()[0])
print(f"    Plant.classified.fasta: {len(cur_plant)} 条")
print(f"    被删 662 条中仍在 fasta 的: {len(gone & cur_plant)}")
resid = []
for cid in cur_plant:
    fm, g, lv, ph = attr.get(cid, ("", "", "", ""))
    if fam_host.get(fm) in PA or gen_host.get(g) in PA:
        resid.append((fm, g, lv.split("(")[0], fam_host.get(fm, ""), gen_host.get(g, "")))
print(f"    fasta 中带原生/藻科或属名的: {len(resid)}")
rc = Counter(resid)
for (fm, g, lv, fh, gh), n in rc.most_common(20):
    print(f"      {fm:<22} {g:<26} {lv:<9} 科判{fh:<8} 属判{gh:<8} x{n}")
print(f"    其中判定层级分布: {dict(Counter(x[2] for x in resid).most_common())}")
