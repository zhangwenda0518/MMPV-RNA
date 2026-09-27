# -*- coding: utf-8 -*-
"""残余风险量化: 科属交叉一致性 + 未登记非植物科清单。

对每张表算三个量:
  R1 冲突行: Family 不在植物库、Genus 却在植物库里  -> 最危险 (非植物科挂植物属)
  R2 倒挂行: Family 在植物库、Genus 不在植物库里     -> 植物科挂陌生属, 需人工看
  R3 属缺行: Family 非空但 Genus 空                   -> 科兜底名单唯一生效的场景

用法: python3 bl_residual_risk.py
"""
import csv
import os
import sys
from collections import Counter, defaultdict

PLANT_DB = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
RUN_HOST = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"

TABLES = [
    ("goji_v62d(新)", "/tmp/v62d_goji/final_integrated_classification.tsv"),
    ("goji_原表", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
                  "/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("goji_植物表", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
                    "/10_Reports/All_plant.viruses_info.tsv"),
    ("onekp_v62d(新)", "/tmp/v62d_onekp/final_integrated_classification.tsv"),
    ("onekp_植物表", "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus"
                     "/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"),
]
RANK_SUFFIX = ("viricota", "viricotina", "viricetes", "viricetidae", "virales",
               "virineae", "virinae", "viridae")


def norm(x):
    v = (x or "").strip()
    return "" if v.lower() in ("", "na", "nan", "none", "null", "-") else v


def split_lineage(s):
    toks = [t.strip() for t in str(s or "").split(";") if t.strip()]
    fam = gen = ""
    for i, t in enumerate(toks):
        if t.endswith("viridae") and not t.endswith("virinae"):
            fam = t
            if i + 1 < len(toks) and not toks[i + 1].endswith(RANK_SUFFIX):
                gen = toks[i + 1]
            break
    return fam, gen


def col_of(cols, *keys):
    for c in cols:
        if all(k in c.lower() for k in keys):
            return c
    return None


def read_tsv(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t")
        return rd.fieldnames or [], list(rd)


def main():
    import ast
    tree = ast.parse(open(RUN_HOST, encoding="utf-8", errors="replace").read())
    fallback = set()
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "NON_PLANT_FAMILIES_FALLBACK" for t in node.targets):
            fallback = set(ast.literal_eval(node.value))

    _, prows = read_tsv(PLANT_DB)
    plant_fam, plant_gen = set(), set()
    for r in prows:
        f, g = split_lineage(r.get("Virus_lineage"))
        if f:
            plant_fam.add(f)
        if g:
            plant_gen.add(g)
    print(f"[植物库] 科 {len(plant_fam)} | 属 {len(plant_gen)} | 宿主名单科兜底 {len(fallback)}")

    gap_fam = defaultdict(Counter)   # 未登记科 -> 表名计数
    for label, path in TABLES:
        if not os.path.exists(path):
            print(f"\n[跳过] {label}: {path} 不存在")
            continue
        cols, rows = read_tsv(path)
        fc, gc = col_of(cols, "family"), col_of(cols, "genus")
        if not fc:
            print(f"\n[跳过] {label}: 无 Family 列 {cols[:6]}")
            continue
        r1 = r2 = r3 = 0
        r1s = []
        r2f = Counter()
        for r in rows:
            f, g = norm(r.get(fc)), norm(r.get(gc))
            if not f:
                continue
            if not g:
                r3 += 1
            elif f not in plant_fam and g in plant_gen:
                r1 += 1
                if len(r1s) < 6:
                    r1s.append((r.get("contig_id") or r.get("Virus_Name") or "?", f, g, norm(r.get("species") or "")))
            elif f in plant_fam and g not in plant_gen:
                r2 += 1
                r2f[f] += 1
            if f not in plant_fam and f not in fallback:
                gap_fam[f][label] += 1
        print("\n" + "=" * 78)
        print(f"{label}: {len(rows)} 行")
        print(f"  R1 危险行 (非植物科挂植物属): {r1}")
        for c, f, g, s in r1s:
            print(f"       {c[:44]:<46}{f:<20}{g:<18}{s[:24]}")
        print(f"  R2 倒挂行 (植物科挂陌生属): {r2}  科分布 {r2f.most_common(6)}")
        print(f"  R3 属缺行 (科兜底名单生效场景): {r3}")

    print("\n" + "=" * 78)
    print("未登记非植物科 (植物库无记录 + 宿主名单未登记), 按出现表计:")
    for f, c in sorted(gap_fam.items(), key=lambda kv: -sum(kv[1].values())):
        print(f"  {f:<24}{dict(c)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
