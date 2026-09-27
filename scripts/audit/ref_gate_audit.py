#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ref_gate_audit.py —— 参照表 genus->family 多义性实测 + 闸门「首行口径 vs 并集口径」差异量化

背景：load_genus_family_ref_dual() 里 unique(ref_dt[...], by = "g_l") 对每个属只留第一行
（该语义从 v6.3 线上原样继承）。若一个属在参照表里出现多个科，判据就只看其中一个，
可能出现「行内 Family 明明与参照表另一行相容，却被判为双侧反驳而置空」。

本脚本做两件事：
  A. 统计参照表里「一属多科」的规模（NCBI 侧 / VMR 侧）。
  B. 对闸门实际置空的行（species_containment_blanked.tsv）逐行回查：
     首行口径判为互斥，但并集口径（任一参照任一行相容即豁免）会放过多少行 —— 即潜在误伤上限。

用法：python3 ref_gate_audit.py <ref_tsv> <blanked_tsv> [<blanked_tsv2> ...]
"""
import csv
import sys
from collections import defaultdict

INVALID = {"-", "NA", "NA", "N/A"}


def norm(x):
    if x is None:
        return None
    x = x.replace('"', "").replace("*", "").strip()
    if x == "":
        return None
    if x.upper() in INVALID:
        return None
    return x


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t", quotechar='"')
        head = next(rd)
        rows = [r for r in rd if r]
    return head, rows


def main():
    ref_path = sys.argv[1]
    blanked = sys.argv[2:]

    head, rows = read_tsv(ref_path)
    idx = {n: i for i, n in enumerate(head)}
    print("参照表: %s" % ref_path)
    print("  行数 %d，列 %s" % (len(rows), ",".join(head)))

    ncbi = defaultdict(set)
    vmr = defaultdict(set)
    n_genus = 0
    n_dropped = 0
    first_ncbi = {}
    first_vmr = {}
    for r in rows:
        g = norm(r[idx["Genus"]])
        if g is None:
            n_dropped += 1
            continue
        g = g.lower()
        n_genus += 1
        nf = norm(r[idx["NCBI_Family"]])
        vf = norm(r[idx["VMR_Family"]])
        if nf is None and vf is None:
            n_dropped += 1
            continue
        if nf is not None:
            nf = nf.lower()
            ncbi[g].add(nf)
            first_ncbi.setdefault(g, nf)
        if vf is not None:
            vf = vf.lower().replace(" ", "")
            vmr[g].add(vf)
            first_vmr.setdefault(g, vf)

    print("  可用行 %d（丢弃含空属名或双参照皆空的 %d 行）" % (n_genus - n_dropped, n_dropped))
    print("  去重后属数: NCBI 侧 %d，VMR 侧 %d" % (len(ncbi), len(vmr)))

    multi_n = [g for g, s in ncbi.items() if len(s) > 1]
    multi_v = [g for g, s in vmr.items() if len(s) > 1]
    print("  [A] 一属多科：NCBI 侧 %d 个属（占 %.3f%%），VMR 侧 %d 个属（占 %.3f%%）"
          % (len(multi_n), 100.0 * len(multi_n) / max(1, len(ncbi)),
             len(multi_v), 100.0 * len(multi_v) / max(1, len(vmr))))
    for tag, lst, mp in (("NCBI", multi_n, ncbi), ("VMR", multi_v, vmr)):
        if lst:
            print("      %s 侧举例（最多 8 个）:" % tag)
            for g in sorted(lst, key=lambda x: -len(mp[x]))[:8]:
                print("        %-28s -> %s" % (g, ",".join(sorted(mp[g]))[:120]))

    for bp in blanked:
        head2, rows2 = read_tsv(bp)
        idx2 = {n: i for i, n in enumerate(head2)}
        n_row = 0
        n_genus_multi = 0
        n_union_exempt = 0
        n_ncbi_any = 0
        n_vmr_any = 0
        ex = []
        for r in rows2:
            fam = norm(r[idx2["Family"]])
            sp = norm(r[idx2["Species"]])
            if fam is None or sp is None or " " not in sp:
                continue
            n_row += 1
            fam_c = fam.lower().replace(" ", "")
            g = sp.split(" ")[0].lower()
            gs_n = ncbi.get(g, set())
            gs_v = vmr.get(g, set())
            if len(gs_n) > 1 or len(gs_v) > 1:
                n_genus_multi += 1
            ok_n = fam_c in gs_n
            ok_v = fam_c in gs_v
            if ok_n:
                n_ncbi_any += 1
            if ok_v:
                n_vmr_any += 1
            if ok_n or ok_v:
                n_union_exempt += 1
                if len(ex) < 12:
                    ex.append("%s | Family=%s | Species=%s | 该属 NCBI科=%s | VMR科=%s"
                              % (r[idx2["contig_id"]][:40], fam, sp,
                                 ",".join(sorted(gs_n))[:80], ",".join(sorted(gs_v))[:80]))
        print("")
        print("  [B] 置空清单 %s" % bp)
        print("      可回查行数 %d（双名法 Species + 有 Family）" % n_row)
        print("      其中该属在参照表中「一属多科」的行 %d" % n_genus_multi)
        print("      并集口径会豁免的行 %d（NCBI 侧任一科相容 %d / VMR 侧 %d）-> 首行口径潜在误伤上限 %.2f%%"
              % (n_union_exempt, n_ncbi_any, n_vmr_any,
                 100.0 * n_union_exempt / max(1, n_row)))
        for s in ex:
            print("        %s" % s)


if __name__ == "__main__":
    main()
