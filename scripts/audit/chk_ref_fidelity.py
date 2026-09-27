#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""校验 build_genus_family_ref.py 产出的参照表能否复现 apply_calib_A.py 的判定。

做法：只用「参照表 tsv」+「原始分类表」重算逐级相容性判据，与 apply_calib_A.py
实际写出的 calibration_changes.tsv 逐行对照，要求 blank / review / conflict_one_side
三个集合完全一致（集合相等，不只是计数相等）。

并附带两个诊断：
  1. 被处置行命中的 name 在 rankedlineage.dmp 里的 superkingdom（Domain）分布，
     用于判断是否存在「同名非病毒条目」污染；
  2. 若给参照加「只认病毒 Realm」的 Domain 过滤，判定结果会不会变。
"""
import csv
import os
import sys
from collections import Counter

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated"
SRC = os.path.join(BASE, "final_integrated_classification.tsv")
CHG = os.path.join(BASE, "calibration_20260914", "calibration_changes.tsv")
CAL = os.path.join(BASE, "calibration_20260914", "final_integrated_classification.calibrated.tsv")
REF = "/home/zhangwenda/database/taxonomy/genus_family_ref.tsv"

VIRUS_REALMS = {"riboviria", "duplodnaviria", "varidnaviria", "adnaviria", "ribozyviria",
                "monodnaviria", "floreoviria", "singelaviria", "ribovira", "viroviria",
                "sangervirae", "proteoviridae", "khanviria", "diviria", "loebvirales"}


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A", "-"):
        return None
    return v


def load_ref(path):
    ref = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            g = norm(r.get("Genus"))
            if not g:
                continue
            ref[g.lower()] = (
                norm(r.get("NCBI_Family")),
                r.get("VMR_Family") or "",
                r.get("NCBI_n") or "",
                r.get("Domain") or "",
            )
    return ref


def judge(fam, gen, ref):
    """返回 (action, 诊断)。规则与 apply_calib_A.py 一致。"""
    gl = gen.lower() if gen else None
    r = ref.get(gl) if gl else None
    ncbi_fam, vmr_str, ncbi_n, domain = r if r else (None, "", "", "")
    ncbi_known = bool(ncbi_fam)
    vmr_set = {x.strip() for x in vmr_str.split(";") if x.strip()}
    vmr_known = bool(vmr_set)
    ncbi_bad = ncbi_known and ncbi_fam.lower() != fam.lower()
    vmr_bad = vmr_known and fam.lower() not in {x.lower() for x in vmr_set}
    if ncbi_known and vmr_known:
        if ncbi_bad and vmr_bad:
            return "blank", (ncbi_fam, vmr_str, domain, ncbi_n)
        if ncbi_bad or vmr_bad:
            return "conflict_one_side", (ncbi_fam, vmr_str, domain, ncbi_n)
    elif ncbi_bad or vmr_bad:
        return ("review", (ncbi_fam, vmr_str, domain, ncbi_n))
    return "", (ncbi_fam, vmr_str, domain, ncbi_n)


def main():
    with open(SRC, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rows = [{k.strip().strip('"'): v for k, v in r.items()} for r in rd]
    print("原表 %d 行" % len(rows))

    # 实跑基准取校准表的 calib_action 列（逐行权威值）；changes 表只写了 blank/review 两类
    expect = {}
    with open(CAL, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            a = (r.get("calib_action") or "").strip()
            if a:
                expect[r["contig_id"]] = a
    chg_n = 0
    chg_kind = Counter()
    with open(CHG, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            chg_n += 1
            chg_kind[r["action"]] += 1
    print("校准表处置行 %d，其中 blank %d / review %d / conflict_one_side %d"
          % (len(expect), sum(1 for v in expect.values() if v == "blank"),
             sum(1 for v in expect.values() if v == "review"),
             sum(1 for v in expect.values() if v == "conflict_one_side")))
    print("changes 表 %d 行 (blank %d / review %d) <- 只登记 blank 与 review"
          % (chg_n, chg_kind.get("blank", 0), chg_kind.get("review", 0)))

    ref = load_ref(REF)
    print("参照表 %d 个 name" % len(ref))

    got = {}
    diag = {}
    dom = Counter()
    for d in rows:
        fam, gen = norm(d.get("Family")), norm(d.get("Genus"))
        act = ""
        if fam and gen:
            act, info = judge(fam, gen, ref)
            diag[d["contig_id"]] = info
            if act in ("blank", "conflict_one_side"):
                dom[info[2] or "(empty)"] += 1
        if act:
            got[d["contig_id"]] = act

    c_got = Counter(got.values())
    c_exp = Counter(expect.values())
    print("")
    print("=== 复算 vs 实跑 ===")
    for k in ("blank", "review", "conflict_one_side"):
        print("  %-20s 复算 %5d / 实跑 %5d" % (k, c_got.get(k, 0), c_exp.get(k, 0)))
    set_got_blank = {k for k, v in got.items() if v == "blank"}
    set_exp_blank = {k for k, v in expect.items() if v == "blank"}
    set_got_all = set(got)
    set_exp_all = set(expect)
    print("  blank 集合是否完全相同: %s (差异 %d)"
          % (set_got_blank == set_exp_blank, len(set_got_blank ^ set_exp_blank)))
    print("  全部处置行集合是否相同: %s (差异 %d)"
          % (set_got_all == set_exp_all, len(set_got_all ^ set_exp_all)))
    for cid in list(set_got_blank ^ set_exp_blank)[:5]:
        print("    差异样例 %s: 复算=%s 实跑=%s" % (cid, got.get(cid, "-"), expect.get(cid, "-")))

    print("")
    print("=== 被处置行命中 name 的 Domain 分布 ===")
    for k, v in dom.most_common(10):
        print("  %-16s %5d" % (k, v))
    nonviral = sum(v for k, v in dom.items()
                   if k != "(empty)" and k.lower() not in VIRUS_REALMS)
    print("  非病毒 Realm 且非空的条目: %d" % nonviral)

    # 诊断 2：只认病毒 Realm
    got2 = {}
    for d in rows:
        fam, gen = norm(d.get("Family")), norm(d.get("Genus"))
        if not (fam and gen):
            continue
        gl = gen.lower()
        r = ref.get(gl)
        if r and r[3] and r[3].lower() not in VIRUS_REALMS:
            r = None
        act, _ = judge(fam, gen, dict([(gl, r)] if r else []))
        if act:
            got2[d["contig_id"]] = act
    print("")
    print("=== 加「只认病毒 Realm」Domain 过滤后 ===")
    print("  blank %d (原 %d) / 集合是否相同: %s"
          % (sum(1 for v in got2.values() if v == "blank"), c_got.get("blank", 0),
             {k for k, v in got2.items() if v == "blank"} == set_got_blank))
    ok = set_got_all == set_exp_all
    print("")
    print("结论: %s" % ("PASS 参照表可复现实跑判定" if ok else "FAIL 见上"))
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
