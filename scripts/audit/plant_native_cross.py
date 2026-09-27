"""原生/藻类与植物的科属交叉程度实测。

三个层面分别测:
  层面1 名称层面: 每个科/属在概率表里只有一个 Predicted_Host 标签 -> 是否真的单标签
  层面2 记录层面: 同一科/属的宿主记录里是否同时含 Plant 与 Protist/Algae
  层面3 层级层面: 同一科下不同属的宿主标签是否分属两个阵营 (需 taxa.txt 属->科映射)
  层面4 数据行层面: 全库 20,892 行里 Family 标签与 Genus 标签是否打架
只读。
"""
import csv
import os
from collections import Counter, defaultdict

CA = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
FULLHOST = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
            "RNA-Lycium_barbarum_out/06_HostPrediction/C9_ICTV_result/classification_result.tsv")
PLANT = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
         "RNA-Lycium_barbarum_out/06_HostPrediction/C9_ICTV_result/Plant.classified.tsv")
VMR = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"
NATIVE = {"Protist", "Algae"}


def rows_of(path):
    with open(path, errors="replace") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            yield {(k.strip() if k else k): (v.strip().strip('"') if isinstance(v, str) else v)
                   for k, v in d.items()}


def load_prob(level):
    out = {}
    for r in rows_of(f"{CA}/{level.lower()}_host_probability.tsv"):
        n = r.get(level, "")
        if n:
            out[n] = r
    return out


def fl(r, k):
    try:
        return float(r.get(k) or 0)
    except Exception:
        return 0.0


fam_p, gen_p = load_prob("Family"), load_prob("Genus")

print("=" * 92)
print("[层面1] 名称层面: 科/属在概率表里是否单标签")
print("=" * 92)
for lvl, tb in (("Family", fam_p), ("Genus", gen_p)):
    hosts = Counter(v["Predicted_Host"] for v in tb.values())
    print(f"  {lvl}: {len(tb)} 条, 出现的 Predicted_Host 种类 {len(hosts)} 种 (每条目唯一标签={len(hosts) > 0})")
nat_f = {k for k, v in fam_p.items() if v["Predicted_Host"] in NATIVE}
plt_f = {k for k, v in fam_p.items() if v["Predicted_Host"] == "Plant"}
nat_g = {k for k, v in gen_p.items() if v["Predicted_Host"] in NATIVE}
plt_g = {k for k, v in gen_p.items() if v["Predicted_Host"] == "Plant"}
print(f"  科名交集 (判Plant ∩ 判原生藻): {len(plt_f & nat_f)}  (空集即名称层不交叉)")
print(f"  属名交集 (判Plant ∩ 判原生藻): {len(plt_g & nat_g)}  (空集即名称层不交叉)")

print("\n" + "=" * 92)
print("[层面2] 记录层面: 同一科内是否同时含 Plant 与新原生/藻记录")
print("=" * 92)


def cross_list(tb, level):
    both, plant_first, nat_first = [], [], []
    for name, v in tb.items():
        pp = fl(v, "P(Plant)")
        pa = fl(v, "P(Protist)") + fl(v, "P(Algae)")
        pf = fl(v, "P(Fungi)")
        if pp > 0 and pa > 0:
            both.append((fl(v, "Total_Records"), name, v["Predicted_Host"], pp, pa, pf, fl(v, "Shannon_Entropy")))
        elif v["Predicted_Host"] == "Plant" and pa > 0:
            plant_first.append(name)
        elif v["Predicted_Host"] in NATIVE and pp > 0:
            nat_first.append(name)
    both.sort(reverse=True)
    return both, plant_first, nat_first


for lvl, tb in (("Family", fam_p), ("Genus", gen_p)):
    both, pf, nf = cross_list(tb, lvl)
    print(f"\n  --- {lvl} 级 ---")
    print(f"  同时含 Plant 与 原生/藻 记录: {len(both)} 个")
    print(f"  标签=Plant 但含原生/藻记录: {len(pf)} 个")
    print(f"  标签=原生/藻 但含 Plant 记录: {len(nf)} 个")
    if both:
        print(f"    {'记录数':>7} {'名称':<26} {'标签':<10} {'P(Plant)':>8} {'P(原生藻)':>9} {'P(Fungi)':>8} {'Shannon':>8}")
        for rec, name, h, pp, pa, pfn, sh in both[:30]:
            print(f"    {rec:>7} {name[:26]:<26} {h:<10} {pp:>8.3f} {pa:>9.3f} {pfn:>8.3f} {sh:>8.3f}")
    if nf:
        print(f"    标签=原生/藻却含植物记录: {', '.join(sorted(nf)[:20])}")

print("\n" + "=" * 92)
print("[层面3] 层级层面: 同一科下不同属的宿主标签是否分属两阵营")
print("=" * 92)
gen2fam = {}
if os.path.exists(VMR):
    with open(VMR, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter=",")
        tc = {c.strip().upper(): c for c in (rd.fieldnames or [])}
        cf, cg = tc.get("FAMILY"), tc.get("GENUS")
        for r in rd:
            g = (r.get(cg) or "").strip().strip('"').strip()
            f = (r.get(cf) or "").strip().strip('"').strip()
            if g and f and g not in gen2fam:
                gen2fam[g] = f
print(f"  属->科 映射: {len(gen2fam)} 属 (源 {VMR})")

fam2gen = defaultdict(list)
for g, f in gen2fam.items():
    fam2gen[f].append(g)

print("\n  A. 科级标签=Plant, 但科下有属级标签为原生/藻:")
n_a = 0
for fam in sorted(plt_f):
    bad = [g for g in fam2gen.get(fam, []) if g in nat_g]
    if bad:
        n_a += 1
        items = [g + "(" + gen_p[g]["Predicted_Host"] + ")" for g in bad[:6]]
        print(f"    {fam:<26} 下 {len(fam2gen[fam])} 属, 其中原生/藻 {len(bad)}: {', '.join(items)}")
print(f"    合计 {n_a} 个植物科")

print("\n  B. 科级标签=原生/藻, 但科下有属级标签为 Plant:")
n_b = 0
for fam in sorted(nat_f):
    bad = [g for g in fam2gen.get(fam, []) if g in plt_g]
    if bad:
        n_b += 1
        print(f"    {fam:<26} ({fam_p[fam]['Predicted_Host']}) 下 {len(fam2gen[fam])} 属, "
              f"其中植物属 {len(bad)}: {', '.join(bad[:6])}")
print(f"    合计 {n_b} 个原生藻科")

print("\n  C. 概率表里属级标签与科级标签不一致的属 (全表):")
n_c = 0
mismatch = []
for g, gf in gen2fam.items():
    if g not in gen_p or gf not in fam_p:
        continue
    gh, fh = gen_p[g]["Predicted_Host"], fam_p[gf]["Predicted_Host"]
    if gh != fh:
        n_c += 1
        mismatch.append((gf, g, fh, gh))
print(f"    可比对属(属级+科级都有条目): "
      f"{sum(1 for g, gf in gen2fam.items() if g in gen_p and gf in fam_p)}, 其中标签不一致 {n_c}")
cc = Counter((fh, gh) for _, _, fh, gh in mismatch)
print("    不一致组合 top12 (科标签 -> 属标签):")
for (fh, gh), v in cc.most_common(12):
    print(f"      {fh:<14} -> {gh:<14} {v:>5}")
print("    原生/藻相关的属级不一致:")
for gf, g, fh, gh in sorted(mismatch):
    if gh in NATIVE or fh in NATIVE:
        print(f"      {gf:<24} 科标签={fh:<12} | {g:<26} 属标签={gh}")

print("\n" + "=" * 92)
print("[层面4] 数据行层面: 全库 20,892 行 Family 标签与 Genus 标签是否打架")
print("=" * 92)
src = FULLHOST if os.path.exists(FULLHOST) else PLANT
combo = Counter()
detail = Counter()
nat_row = 0
for r in rows_of(src):
    fam, gen = r.get("Family", ""), r.get("Genus", "")
    fh = fam_p.get(fam, {}).get("Predicted_Host", "")
    gh = gen_p.get(gen, {}).get("Predicted_Host", "")
    if not fh or not gh:
        combo[("无完整标签", "")] += 1
        continue
    if fh == gh:
        combo[("一致", fh)] += 1
    else:
        combo[("不一致", f"{fh}->{gh}")] += 1
        if fh in NATIVE or gh in NATIVE:
            nat_row += 1
            detail[(fh, gh)] += 1
tot = sum(combo.values())
print(f"  源: {src}")
n_same = sum(v for (k, _), v in combo.items() if k == "一致")
n_diff = sum(v for (k, _), v in combo.items() if k == "不一致")
n_na = tot - n_same - n_diff
print(f"  总行 {tot}: 科属标签一致 {n_same} ({n_same / tot * 100:.1f}%) / "
      f"不一致 {n_diff} ({n_diff / tot * 100:.1f}%) / 缺一侧条目 {n_na}")
print(f"  不一致且涉及原生/藻的行: {nat_row} ({nat_row / tot * 100:.1f}%)")
print("  涉原生/藻的不一致组合 top15:")
for (fh, gh), v in detail.most_common(15):
    print(f"    {fh:<14} -> {gh:<14} {v:>6}")
