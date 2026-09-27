#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐级相容性检查：统计 final_integrated_classification.tsv 里"相邻等级跨树"的嵌合行。

判据：对每一对相邻等级 (H 高, L 低)，用参照库查 L 的定型父级，与行内 H 比较；不等即矛盾。
参照：/home/zhangwenda/database/taxonomy/rankedlineage.dmp
      列序 tax_id | tax_name | species | genus | family | order | class | phylum | kingdom | superkingdom
      （病毒部分：superkingdom 列 = Realm，kingdom 列 = Kingdom）
只读，不改任何文件。
"""
import csv
import os
import sys
from collections import Counter, defaultdict

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy"
INTEG = os.path.join(BASE, "Votus.integrated", "final_integrated_classification.tsv")

# 由低到高
RANKS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
# rankedlineage 列下标（split 后）
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
# 相邻对：(高, 低, 低等级在参照里的父级列)
PAIRS = [
    ("Realm", "Kingdom", "Realm"),
    ("Kingdom", "Phylum", "Kingdom"),
    ("Phylum", "Class", "Phylum"),
    ("Class", "Order", "Class"),
    ("Order", "Family", "Order"),
    ("Family", "Genus", "Family"),
    ("Genus", "Species", "Genus"),
]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() == "NA" or v.upper() == "N/A":
        return None
    return v


def low(v):
    return v.lower() if v else None


def parse_agree(v):
    """'3/5: CAT,diamond_lca,mmseqs' -> set of tools"""
    if not v or ":" not in v:
        return None
    s = {t.strip() for t in v.split(":", 1)[1].split(",") if t.strip()}
    return s or None


def main():
    # ---- 1. 读整合表 ----
    rows = []
    with open(INTEG, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = [c.strip().strip('"') for c in rd.fieldnames]
        for r in rd:
            d = {}
            for k, v in r.items():
                d[k.strip().strip('"')] = norm(v)
            rows.append(d)
    print("整合表行数: %d" % len(rows))
    print("列名: %s" % cols)

    # ---- 2. 收集需要查的 name ----
    need = set()
    for d in rows:
        for L in RANKS:
            v = low(d.get(L))
            if v:
                need.add(v)
    print("需要查参照的 name 数: %d" % len(need))

    # ---- 3. 扫 rankedlineage.dmp，只留有需要的 ----
    ref = {}
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            # NCBI dmp 每行以 "\t|" 结尾，末字段会带上尾巴竖线，必须剥掉
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in need and nm not in ref:
                ref[nm] = {L: (parts[REF_IDX[L]] or None) for L in RANKS}
    print("命中参照的 name 数: %d" % len(ref))
    print("--- 参照解析自检 ---")
    for probe in ["partitiviridae", "biavirus", "caudoviricetes", "kyanoviridae", "chlorovirus"]:
        r = ref.get(probe)
        if r:
            print("  %-18s family=%s order=%s class=%s phylum=%s kingdom=%s realm=%s" % (
                probe, r.get("Family"), r.get("Order"), r.get("Class"), r.get("Phylum"), r.get("Kingdom"), r.get("Realm")))
        else:
            print("  %-18s <未命中>" % probe)

    # ---- 4. 逐级相容性检查 ----
    no_conflict = 0
    conflicts = []            # (row, [矛盾边界描述])
    boundary_counter = Counter()
    pair_combo = Counter()    # (Family, Genus) 组合
    unmatched = defaultdict(int)
    n_self_inconsistent = 0

    data_rows = []
    for d in rows:
        cid = d.get("contig_id")
        if not cid:
            continue
        bad = []
        for hi, lo, refcol in PAIRS:
            vh = d.get(hi)
            vl = d.get(lo)
            if not vh or not vl:
                continue
            rl = ref.get(low(vl))
            if rl is None:
                unmatched[lo] += 1
                continue
            expect = rl.get(refcol)
            if not expect:
                continue
            if low(expect) != low(vh):
                bad.append((hi, lo, vh, vl, expect))
                boundary_counter["%s <-> %s" % (hi, lo)] += 1
        rec = {"cid": cid, "row": d, "bad": bad, "srcbad": []}
        for hi, lo, _refcol in PAIRS:
            sh = parse_agree(d.get(hi + "_agree"))
            sl = parse_agree(d.get(lo + "_agree"))
            if sh and sl and sh.isdisjoint(sl):
                rec["srcbad"].append((hi, lo))
        data_rows.append(rec)
        if bad:
            conflicts.append(rec)
            if any(b[1] == "Genus" for b in bad):
                pair_combo[(d.get("Family") or "NA", d.get("Genus") or "NA")] += 1
        else:
            no_conflict += 1

    n_rows = len(data_rows)
    n_bad = len(conflicts)
    print("")
    print("=" * 70)
    print("总数据行: %d" % n_rows)
    print("无矛盾行: %d (%.2f%%)" % (no_conflict, 100.0 * no_conflict / n_rows))
    print("有矛盾行: %d (%.2f%%)" % (n_bad, 100.0 * n_bad / n_rows))

    fam_gen = [r for r in conflicts if any(b[0] == "Family" and b[1] == "Genus" for b in r["bad"])]
    print("")
    print("口径A（仅科-属矛盾 / 最小切口）: %d 行" % len(fam_gen))
    print("口径A+（任意相邻等级矛盾）    : %d 行" % n_bad)
    multi = [r for r in conflicts if len(r["bad"]) > 1]
    print("其中含 2 处以上矛盾的行        : %d 行" % len(multi))

    print("")
    print("--- 各边界的矛盾计数 ---")
    for k, v in boundary_counter.most_common():
        print("  %-22s %6d" % (k, v))

    print("")
    print("--- 矛盾行的科-属组合 top 25 ---")
    for (f, g), c in pair_combo.most_common(25):
        print("  %-28s + %-30s %6d" % (f, g, c))

    print("")
    print("--- 无法在参照里判定的低等级计数（name 未命中）---")
    for k, v in sorted(unmatched.items(), key=lambda x: -x[1]):
        print("  %-12s %7d 次" % (k, v))

    print("")
    print("--- 矛盾行的 Genus_agree 分布 top 15 ---")
    ga = Counter((r["row"].get("Genus_agree") or "NA") for r in conflicts)
    for k, v in ga.most_common(15):
        print("  %-40s %6d" % (k, v))

    single = sum(1 for r in conflicts if (r["row"].get("Genus_agree") or "").startswith("1/1"))
    print("  矛盾行里 Genus_agree 为单票(1/1:*)的: %d 行" % single)

    print("")
    print("--- 边界组合签名 top 15 ---")
    sig = Counter()
    for r in conflicts:
        key = "+".join(sorted("%s<->%s" % (b[0], b[1]) for b in r["bad"]))
        sig[key] += 1
    for k, v in sig.most_common(15):
        print("  %-46s %6d" % (k, v))

    print("")
    print("--- 矛盾行的 primary_tool 分布 ---")
    for k, v in Counter((r["row"].get("primary_tool") or "NA") for r in conflicts).most_common(12):
        print("  %-16s %6d" % (k, v))

    print("")
    print("--- 若按口径A执行（空白 Genus+Species）的影响面 ---")
    n_gen_blank = 0
    n_sp_blank = 0
    for r in fam_gen:
        if r["row"].get("Genus"):
            n_gen_blank += 1
        if r["row"].get("Species"):
            n_sp_blank += 1
    print("  待置空 Genus 单元格: %d" % n_gen_blank)
    print("  待置空 Species 单元格: %d" % n_sp_blank)

    print("")
    print("--- 口径A 样例 20 行 ---")
    for r in fam_gen[:20]:
        d = r["row"]
        print("  %s" % d.get("contig_id"))
        print("      Realm=%s Kingdom=%s Phylum=%s Class=%s" % (d.get("Realm"), d.get("Kingdom"), d.get("Phylum"), d.get("Class")))
        print("      Order=%s Family=%s Genus=%s Species=%s" % (d.get("Order"), d.get("Family"), d.get("Genus"), d.get("Species")))
        print("      primary=%s | Fam_agree=%s | Gen_agree=%s" % (d.get("primary_tool"), d.get("Family_agree"), d.get("Genus_agree")))
        for b in r["bad"]:
            print("      >> 矛盾 %s(%s) vs %s(%s) 的参照父级=%s" % (b[0], b[2], b[1], b[3], b[4]))

    print("")
    print("--- 第二参照交叉验证（VMR MSL41 taxa.txt 的 Family+Genus）---")
    TAXA = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"
    gmap = defaultdict(set)
    gmap_seen = 0
    try:
        with open(TAXA, newline="", encoding="utf-8", errors="replace") as fh:
            trd = csv.DictReader(fh, delimiter=",")
            tcols = {c.strip().upper(): c for c in (trd.fieldnames or [])}
            print("  taxa.txt 表头: %s" % (trd.fieldnames,))
            cf, cg = tcols.get("FAMILY"), tcols.get("GENUS")
            for r in trd:
                g = norm(r.get(cg))
                f = norm(r.get(cf))
                if g and f:
                    gmap[g.lower()].add(f)
                    gmap_seen += 1
        print("  taxa.txt 可用 genus->family 映射: %d 条 / %d 个属" % (gmap_seen, len(gmap)))
    except Exception as e:
        print("  taxa.txt 读取失败: %s" % e)

    buck = Counter()
    detail_rows = []
    for r in fam_gen:
        d = r["row"]
        g = low(d.get("Genus"))
        f = d.get("Family")
        fams = gmap.get(g) if g else None
        if not fams:
            buck["①仅单参照可判（第二参照无此属）"] += 1
            detail_rows.append(("单参照", r))
        elif low(f) in {low(x) for x in fams}:
            buck["②第二参照与行内 Family 一致（版本漂移嫌疑，待复核）"] += 1
            detail_rows.append(("版本嫌疑", r))
        else:
            buck["③两参照一致反驳行内 Family（强判嵌合）"] += 1
            detail_rows.append(("强判", r))
    for k in sorted(buck):
        print("  %-52s %6d" % (k, buck[k]))

    print("")
    print("--- 两判据关系（科-属边界）---")
    R = {i for i, r in enumerate(data_rows)
         if any(b[0] == "Family" and b[1] == "Genus" for b in r["bad"])}
    S = {i for i, r in enumerate(data_rows)
         if any(b[0] == "Family" and b[1] == "Genus" for b in r["srcbad"])}
    print("  R 参照库跨树不相容（可判错）: %d" % len(R))
    print("  S 票源不相交（可疑）        : %d" % len(S))
    print("  R ∩ S 两者同时命中          : %d" % len(R & S))
    print("  R 但非 S（同源却内部矛盾）  : %d" % len(R - S))
    print("  S 但非 R（异源但值相容）    : %d" % len(S - R))
    print("  --- 非 R 的 S 行样例（异源但谱系相容，说明 S 不能单独当判错依据）---")
    shown = 0
    for i in sorted(S - R):
        d = data_rows[i]["row"]
        print("    %s | Family=%s (%s) Genus=%s (%s)" % (
            d.get("contig_id"), d.get("Family"), d.get("Family_agree"),
            d.get("Genus"), d.get("Genus_agree")))
        shown += 1
        if shown >= 6:
            break

    print("")
    print("--- Biavirus 行诊断 ---")
    bv_all = [r for r in data_rows if (r["row"].get("Genus") or "").lower() == "biavirus"]
    print("  Genus=Biavirus 行数: %d" % len(bv_all))
    print("  其 Family 取值分布:")
    for k, v in Counter((r["row"].get("Family") or "NA") for r in bv_all).most_common(12):
        print("     %-24s %5d" % (k, v))
    bsig = Counter()
    for r in bv_all:
        bsig["+".join(sorted("%s<->%s" % (b[0], b[1]) for b in r["bad"])) or "无矛盾"] += 1
    print("  其矛盾边界分布:")
    for k, v in bsig.most_common(12):
        print("     %-46s %5d" % (k, v))
    print("  前 6 行矛盾明细:")
    for r in bv_all[:6]:
        d = r["row"]
        print("    %s" % d.get("contig_id"))
        print("      F=%s G=%s S=%s | %s" % (d.get("Family"), d.get("Genus"), d.get("Species"), d.get("Genus_agree")))
        for b in r["bad"]:
            print("      >> %s(%s) vs %s(%s) 参照父级=%s" % (b[0], b[2], b[1], b[3], b[4]))

    print("")
    print("--- Biavirus 相关行确认 ---")
    bv = [r for r in data_rows if (r["row"].get("Genus") or "").lower() == "biavirus"]
    print("  整合表中 Genus=Biavirus 的行: %d" % len(bv))
    print("  其中被判为矛盾: %d" % sum(1 for r in bv if r["bad"]))
    for r in bv[:5]:
        d = r["row"]
        print("    %s -> Family=%s Genus=%s Species=%s | %s" % (
            d.get("contig_id"), d.get("Family"), d.get("Genus"), d.get("Species"), "矛盾" if r["bad"] else "相容"))


if __name__ == "__main__":
    main()
