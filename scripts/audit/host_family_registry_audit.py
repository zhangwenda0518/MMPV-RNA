# -*- coding: utf-8 -*-
"""宿主分类名单完备性审计 (修正版: Plant.tsv 的科名要从 Virus_lineage 解析)。

背景: 05 分类表出现过 Family=Mimiviridae/Genus=Potyvirus 这类科属打架行,
      科是巨大病毒科、属是植物病毒属。要回答两件事:
        (1) 这些科是不是植物病毒科 -> 决定要不要登记进宿主非植物名单;
        (2) 除它们之外还有哪些科可能同样有问题 -> 一次扫全, 免得漏。

判定证据链 (三层, 逐层收窄):
  L1 自家植物库 Plant.tsv 的 Virus_lineage 里, 该科有没有植物病毒记录
  L2 宿主名单 (NON_PLANT_FAMILIES_FALLBACK / NON_PLANT_GENERA) 有没有登记
  L3 成品表 / 下游植物表里该科实际有多少行, 属是哪些
"""
import ast
import csv
import os
import sys
from collections import Counter, defaultdict

RUN_HOST = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
PLANT_DB = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

TREES = [
    ("goji", "/tmp/v62d_goji/final_integrated_classification.tsv",
     "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
     "/10_Reports/All_plant.viruses_info.tsv"),
    ("onekp", "/tmp/v62d_onekp/final_integrated_classification.tsv",
     "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus"
     "/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"),
]

TARGET_FAMS = ["Mimiviridae", "Marseilleviridae", "Pithoviridae"]

RANK_SUFFIX = ("viricota", "viricotina", "viricetes", "viricetidae", "virales",
               "virineae", "virinae", "viridae")


def norm(x):
    v = (x or "").strip()
    if v.lower() in ("", "na", "nan", "none", "null", "-"):
        return ""
    return v


def load_lists(path):
    """只取模块级字面量赋值, 不执行模块体 (模块体里有 Path/环境依赖)。"""
    tree = ast.parse(open(path, encoding="utf-8", errors="replace").read())
    out = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for tgt in node.targets:
            if isinstance(tgt, ast.Name):
                try:
                    out[tgt.id] = ast.literal_eval(node.value)
                except Exception:
                    pass
    return out


def col_of(cols, *keys):
    for c in cols:
        lc = c.lower()
        if all(k in lc for k in keys):
            return c
    return None


def read_tsv(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t")
        return rd.fieldnames or [], list(rd)


def split_lineage(s):
    """返回 (family, genus)。family=以 viridae 结尾的 token; genus=紧随其后的 token。"""
    toks = [t.strip() for t in str(s or "").split(";") if t.strip()]
    fam = gen = ""
    for i, t in enumerate(toks):
        if t.endswith("viridae") and not t.endswith("virinae"):
            fam = t
            if i + 1 < len(toks):
                nxt = toks[i + 1]
                if not nxt.endswith(RANK_SUFFIX):
                    gen = nxt
            break
    return fam, gen


def main():
    ns = load_lists(RUN_HOST)
    fallback = set(ns.get("NON_PLANT_FAMILIES_FALLBACK") or [])
    genera_bl = set(ns.get("NON_PLANT_GENERA") or [])
    print(f"[宿主名单] 科兜底 {len(fallback)} 条 | 属黑名单 {len(genera_bl)} 条")

    # ── L1: Plant.tsv 植物病毒权威库 (科名从 lineage 解析) ──
    pcols, prows = read_tsv(PLANT_DB)
    plant_fam = Counter()
    plant_fam_gen = defaultdict(Counter)
    for r in prows:
        f, g = split_lineage(r.get("Virus_lineage"))
        if f:
            plant_fam[f] += 1
            if g:
                plant_fam_gen[f][g] += 1
    plant_genera = set(g for c in plant_fam_gen.values() for g in c)
    print(f"\n[L1 植物库] {len(prows)} 行 | 解析出科 {len(plant_fam)} 个 | 属 {len(plant_genera)} 个")
    print("  Top15 科: %s" % plant_fam.most_common(15))

    # ── A. 用户点名的三科 ──
    print("\n" + "=" * 78)
    print("A. 专项核对: %s" % TARGET_FAMS)
    for f in TARGET_FAMS:
        print(f"  {f:<18} 植物库记录 {plant_fam.get(f,0):>5} 行 | 植物库属 {sorted(plant_fam_gen[f])[:5]}"
              f" | 科兜底 {'已登记' if f in fallback else '未登记'}")

    # ── A2. 现有 38 条科兜底名单自身审计 (两边看: 误杀 / 漏登) ──
    print("\n" + "=" * 78)
    print("A2. 现有科兜底名单的植物库交叉检查")
    in_db = [(f, plant_fam[f]) for f in sorted(fallback) if f in plant_fam]
    no_db = [f for f in sorted(fallback) if f not in plant_fam]
    print(f"  名单里在植物库有记录的科 {len(in_db)} 个 (整科否决会有误杀风险):")
    for f, n in in_db:
        print(f"    {f:<22} 植物库 {n:>5} 行 | 植物库属样本 {sorted(plant_fam_gen[f])[:4]}")
    print(f"  名单里植物库无记录的科 {len(no_db)} 个 (符合建表口径)")

    # ── B/C: 两棵树 05 表全量科扫描 ──
    for label, tax_path, plant_path in TREES:
        if not os.path.exists(tax_path):
            print(f"\n[跳过] {label}: {tax_path} 不存在")
            continue
        cols, rows = read_tsv(tax_path)
        fc, gc = col_of(cols, "family"), col_of(cols, "genus")
        fam_rows, fam_gen = Counter(), defaultdict(Counter)
        for r in rows:
            f, g = norm(r.get(fc)), norm(r.get(gc))
            if not f:
                continue
            fam_rows[f] += 1
            if g:
                fam_gen[f][g] += 1
        print("\n" + "=" * 78)
        print(f"B. {label} 05 表: {len(rows)} 行 | 非空科 {sum(fam_rows.values())} 行 | 科名 {len(fam_rows)} 个")
        unknown = sorted(((f, n) for f, n in fam_rows.items()
                          if f not in plant_fam and f not in fallback), key=lambda t: -t[1])
        print(f"  植物库有记录 {sum(1 for f in fam_rows if f in plant_fam)} 科 | "
              f"已登记非植物 {sum(1 for f in fam_rows if f in fallback)} 科 | "
              f"两者都不是 {len(unknown)} 科 / {sum(n for _, n in unknown)} 行")
        print(f"\n  C. 未在植物库、也未登记非植物的科 (按行数):")
        for f, n in unknown:
            top = ", ".join(f"{g}({c})" for g, c in fam_gen[f].most_common(3))
            hit = sum(1 for g in fam_gen[f] if g in genera_bl)
            print(f"    {f:<24}{n:>6}  {top[:56]:<58}{'属已部分命中' if hit else ''}")

        # ── D. 下游植物表 ──
        if not os.path.exists(plant_path):
            print(f"\n  [跳过 D] {label}: {plant_path} 不存在")
            continue
        pcols2, prows2 = read_tsv(plant_path)
        fc2, gc2 = col_of(pcols2, "family"), col_of(pcols2, "genus")
        pf, pg = Counter(), defaultdict(Counter)
        for r in prows2:
            f, g = norm(r.get(fc2)), norm(r.get(gc2))
            if f:
                pf[f] += 1
                if g:
                    pg[f][g] += 1
        print(f"\n  D. {label} 下游植物表: {len(prows2)} 行 | 科名 {len(pf)} 个")
        notin = sorted(((f, n) for f, n in pf.items() if f not in plant_fam), key=lambda t: -t[1])
        print(f"    植物库里没有记录的科 {len(notin)} 个 / {sum(n for _, n in notin)} 行:")
        for f, n in notin:
            tag = "已登记非植物" if f in fallback else "未登记"
            gens = ", ".join(g for g, _ in pg[f].most_common(3))
            print(f"      {f:<22}{n:>6}  {tag:<12} 属: {gens[:44]}")
        still = sorted(((f, n) for f, n in pf.items() if f in fallback), key=lambda t: -t[1])
        print(f"    仍在植物表里的\"已登记非植物科\" {len(still)} 科 / {sum(n for _, n in still)} 行"
              f" (科级名单没拦住, 因为属/种级判定绕过名单):")
        for f, n in still[:12]:
            gens = ", ".join(g for g, _ in pg[f].most_common(3))
            print(f"      {f:<22}{n:>6} 属: {gens[:44]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
