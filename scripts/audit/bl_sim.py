"""黑名单机制定量模拟 (只读, 不写结果).
问题: 靠"补原生/藻黑名单"能不能把混进植物表的行去掉?
做法:
  A. C7 概率表里判 Protist/Algae 的科属 vs 现有黑名单的覆盖差集
  B. 两棵树 Plant 行逐行跑服务器现版 is_blacklisted: 命中 / 漏网 + 漏网原因四分类
  C. 模拟补名单后的增量拦截与免疫行 (属/种级) 数量
  D. 植物表内"科属打架"证据: 科判植物、属判原生/藻 的行
"""
import ast
import csv
import os
from collections import Counter, defaultdict

RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
C7 = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"

TREES = {
    "goji-Lycium(未过黑名单)": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp(已过黑名单-662)": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}

# ---------- 1. 从服务器源码 AST 取黑名单三件套 ----------
src = open(RUN, encoding="utf-8", errors="replace").read()
tree = ast.parse(src)
sets = [n for n in tree.body if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id in (
            "NON_PLANT_GENERA", "NON_PLANT_FAMILIES_FALLBACK", "TRUSTED_LEVELS")
            for t in n.targets)]
fns = [n for n in tree.body if isinstance(n, ast.FunctionDef)
       and n.name in ("is_trusted_level", "is_blacklisted")]
ns = {}
exec(compile(ast.Module(body=sets + fns, type_ignores=[]), "<m>", "exec"), ns)
is_blacklisted = ns["is_blacklisted"]
is_trusted_level = ns["is_trusted_level"]
NG = set(ns["NON_PLANT_GENERA"])
NF = set(ns["NON_PLANT_FAMILIES_FALLBACK"])
TL = ns["TRUSTED_LEVELS"]
print("=" * 78)
print(f"[0] 服务器黑名单: 属 {len(NG)} 项 / 科兜底 {len(NF)} 项 / 豁免层级 {sorted(TL)}")
print("=" * 78)

# ---------- 2. C7 概率表: 按 Predicted_Host 分组的科/属名 ----------
def load_c7(fn, namecol):
    d = defaultdict(set)
    with open(f"{C7}/{fn}", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            h = (row.get("Predicted_Host") or "").strip()
            n = (row.get(namecol) or "").strip()
            if n:
                d[h].add(n)
    return d

fam_by_host = load_c7("family_host_probability.tsv", "Family")
gen_by_host = load_c7("genus_host_probability.tsv", "Genus")
PA = ("Protist", "Algae")          # 原生/藻
pa_fam = set()
pa_gen = set()
for h in fam_by_host:
    if h in PA:
        pa_fam |= fam_by_host[h]
for h in gen_by_host:
    if h in PA:
        pa_gen |= gen_by_host[h]
plant_fam = set(fam_by_host.get("Plant", set()))
plant_gen = set(gen_by_host.get("Plant", set()))

print(f"\n[1] C7 概率表宿主分布 (科 {sum(len(v) for v in fam_by_host.values())} 行 / 属 {sum(len(v) for v in gen_by_host.values())} 行)")
for h in sorted(fam_by_host, key=lambda x: -len(fam_by_host[x])):
    print(f"    科层 {h:<12} {len(fam_by_host[h]):>4}")
for h in sorted(gen_by_host, key=lambda x: -len(gen_by_host[x])):
    print(f"    属层 {h:<12} {len(gen_by_host[h]):>4}")

print(f"\n[2] C7 判原生/藻 与 现有黑名单 的覆盖差集")
miss_fam = sorted(pa_fam - NF)
hit_fam = sorted(pa_fam & NF)
miss_gen = sorted(pa_gen - NG)
hit_gen = sorted(pa_gen & NG)
print(f"    原生/藻科 {len(pa_fam)} 个: 已在科兜底 {len(hit_fam)}, 漏登 {len(miss_fam)}")
print(f"      已在: {hit_fam}")
print(f"      漏登: {miss_fam}")
print(f"    原生/藻属 {len(pa_gen)} 个: 已在属黑名单 {len(hit_gen)}, 漏登 {len(miss_gen)}")
print(f"      已在: {hit_gen}")
print(f"      漏登: {miss_gen}")
print(f"\n    同名冲突检查 (既判原生/藻 又判植物):")
print(f"      科: {sorted(pa_fam & plant_fam)}")
print(f"      属: {sorted(pa_gen & plant_gen)}")

# ---------- 3. 逐树模拟 ----------
for label, root in TREES.items():
    f9 = f"{root}/06_HostPrediction/C9_ICTV_result/classification_result.tsv"
    print("\n" + "=" * 78)
    print(f"[3] {label}")
    print(f"     {f9}")
    print("=" * 78)
    if not os.path.isfile(f9):
        print("    文件不存在")
        continue

    rows = []
    with open(f9, errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    plant = [r for r in rows if (r.get("Predicted_Host") or "").strip() == "Plant"]
    print(f"    全表 {len(rows)} 行, 判 Plant {len(plant)} 行")

    caught = 0
    reason = Counter()
    caught_det = Counter()
    immune_level = Counter()
    immune_with_pa = 0
    plant_fam_pa_genus = []
    r3_fam = Counter()
    r3_det = Counter()
    r2_gen = Counter()
    immune_pa = []
    for r in plant:
        fam = (r.get("Family") or "").strip()
        gen = (r.get("Genus") or "").strip()
        lv = (r.get("Determination_Level") or "").strip()
        trusted = is_trusted_level(lv)
        if trusted:
            immune_level[lv.split("(")[0]] += 1
            if fam in pa_fam or gen in pa_gen:
                immune_with_pa += 1
        if is_blacklisted(fam, gen, lv):
            caught += 1
            caught_det[(fam, gen)] += 1
            continue
        # 漏网原因
        if trusted:
            reason["R1 属/种级判定豁免 (任何名单都管不到)"] += 1
        elif gen and gen not in ("NA", "nan"):
            if fam in NF:
                reason["R5 有属名->科兜底被跳过 (科已在名单但属不在)"] += 1
            else:
                reason["R2 有属名且属不在名单 (含属名像植物)"] += 1
        elif fam and fam in NF:
            reason["R3b 无属+科在名单 (理论应已被拦, 异常)"] += 1
        elif fam:
            reason["R3 无属+科不在科兜底名单"] += 1
            r3_fam[fam] += 1
            r3_det[lv] += 1
        else:
            reason["R4 无属无科"] += 1
        if not trusted and gen and gen not in ("NA", "nan"):
            r2_gen[(fam, gen)] += 1
        if trusted and (fam in pa_fam or gen in pa_gen):
            immune_pa.append((fam, gen, lv, (r.get("Predicted_Host") or "").strip()))
        if fam in plant_fam and gen in pa_gen:
            plant_fam_pa_genus.append((fam, gen, lv))

    print(f"\n    现有黑名单命中: {caught} 行 ({caught/len(plant)*100:.2f}%)")
    print(f"    命中明细 top15:")
    for (fm, g), n in caught_det.most_common(15):
        print(f"      {fm:<24} {g:<22} {n}")
    print(f"\n    漏网 {len(plant)-caught} 行, 原因分类:")
    for k, v in reason.most_common():
        print(f"      {k:<48} {v}")
    print(f"\n    豁免层级分布: {dict(immune_level)}")
    print(f"    其中科或属命中原生/藻名字的豁免行: {immune_with_pa}")
    if immune_pa:
        print("      -> 明细 (科 / 属 / 判定层级):")
        for fm, g, lv, ph in immune_pa:
            print(f"         {fm:<22} {g:<24} {lv}")
    print(f"\n    R3 '无属且科不在名单' {sum(r3_fam.values())} 行的科分布:")
    for fm, n in r3_fam.most_common(25):
        mark = "  <== C7 判原生/藻科" if fm in pa_fam else ""
        print(f"      {fm:<28} {n:>4}{mark}")
    print(f"    R3 行判定层级: {dict(r3_det)}")
    print(f"\n    R2 '有属名但属不在名单' {sum(r2_gen.values())} 行的科-属分布:")
    for (fm, g), n in r2_gen.most_common(25):
        mark = ""
        if g in pa_gen:
            mark += "  <== C7 判原生/藻属"
        if g in plant_gen:
            mark += "  <== C7 也判植物属"
        print(f"      {fm:<24} {g:<24} {n:>4}{mark}")

    # D. 补名单增量模拟
    add_fam = pa_fam - NF
    add_gen = pa_gen - NG
    inc_fam = inc_gen = 0
    inc_fam_det = Counter()
    inc_gen_det = Counter()
    for r in plant:
        fam = (r.get("Family") or "").strip()
        gen = (r.get("Genus") or "").strip()
        lv = (r.get("Determination_Level") or "").strip()
        if is_blacklisted(fam, gen, lv):
            continue
        if is_trusted_level(lv):
            continue
        if gen and gen not in ("NA", "nan"):
            if gen in add_gen:
                inc_gen += 1
                inc_gen_det[(fam, gen)] += 1
        elif fam and fam in add_fam:
            inc_fam += 1
            inc_fam_det[fam] += 1
    print(f"\n    补名单增量 (仅对非豁免行生效):")
    print(f"      补 {len(add_gen)} 个原生/藻属 -> 多拦 {inc_gen} 行  {dict(inc_gen_det)}")
    print(f"      补 {len(add_fam)} 个原生/藻科到科兜底 -> 多拦 {inc_fam} 行  {dict(inc_fam_det)}")
    print(f"      合计可多拦 {inc_gen + inc_fam} 行 / 漏网 {len(plant)-caught} 行")

    # D2. 科属打架
    print(f"\n    植物表内'科判植物+属判原生/藻': {len(plant_fam_pa_genus)} 行")
    for (fm, g, lv), n in Counter(plant_fam_pa_genus).most_common(15):
        print(f"      {fm:<24} {g:<22} {lv}  x{n}")

print("\n" + "=" * 78)
print("[4] ensemble_host_summary.tsv 列 (判断能否回放 RVH/PB2 复活路径)")
print("=" * 78)
for label, root in TREES.items():
    p = f"{root}/06_HostPrediction/ensemble_host_summary.tsv"
    if os.path.isfile(p):
        with open(p, errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
        print(f"\n{label}: {len(hdr)} 列")
        print("   ", " | ".join(hdr))
