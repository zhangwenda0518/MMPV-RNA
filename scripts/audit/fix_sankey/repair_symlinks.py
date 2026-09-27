#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""修复数据集树内断链 symlink（2026-09-12 relink 遗留的旧绝对路径失效）。

用法:
  python3 repair_symlinks.py [--fix] [--ds SUBSTR]
    --ds 只处理名字包含 SUBSTR 的数据集（用于分片跑, 避免 onekp 巨型树超时）

规则 (按 basename):
  <s>_virus.all.candidate.fasta              -> 02a_Identification/<s>/<同名>
  <s>.virus.candidate_uniprot_filtered.fasta -> <link 所在目录>/uniprot/all_candidates_uniprot.fasta
  <s>.virus.candidate_filtered.fasta         -> <link 所在目录>/uniprot/all_candidates_uniprot.fasta
  <s>.virus.candidate_cdd_filtered.fasta     -> <link 所在目录>/cdd/all_candidates_cdd_filter.fasta
  其它: 同数据集树内同名文件唯一匹配 (generic, 用 find -name 定位)
输出: /tmp/symlink_repair_20260916/<ds>.tsv
"""
import os
import subprocess
import sys
import time
from pathlib import Path

BASE = Path("/home/zhangwenda/MMPV-paper")
ROOTS = [
    BASE / "goji-virome/02_novel_virus/RNA-Alternaria_alternata_out",
    BASE / "goji-virome/02_novel_virus/RNA-Aphis_gossypii_out",
    BASE / "goji-virome/02_novel_virus/RNA-Fusarium_nematophilum_out",
    BASE / "goji-virome/02_novel_virus/RNA-Lycium_amarum_out",
    BASE / "goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    BASE / "goji-virome/02_novel_virus/RNA-Lycium_chinense_out",
    BASE / "goji-virome/02_novel_virus/RNA-Lycium_ruthenicum_out",
    BASE / "onekp-virome/onekp-virus",
]
OUT = Path("/tmp/symlink_repair_20260916")


def find_broken(ds):
    """用 find -xtype l 列出断链（比 os.walk 快得多）"""
    r = subprocess.run(["find", str(ds), "-xtype", "l", "-print"],
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    return [Path(x) for x in r.stdout.decode("utf-8", "replace").splitlines() if x.strip()]


def rule_target(link, ds):
    parent = link.parent
    name = link.name
    s = parent.name
    if name.endswith("_virus.all.candidate.fasta"):
        cands = [(ds / "02a_Identification" / s / name, "02a_all_candidate")]
    elif name.endswith(".virus.candidate_cdd_filtered.fasta"):
        cands = [(parent / "cdd" / "all_candidates_cdd_filter.fasta", "cdd_filtered")]
    elif name.endswith(".virus.candidate_uniprot_filtered.fasta"):
        cands = [(parent / "uniprot" / "all_candidates_uniprot.fasta", "uniprot_filtered")]
    elif name.endswith(".virus.candidate_filtered.fasta"):
        cands = [(parent / "uniprot" / "all_candidates_uniprot.fasta", "strict_from_uniprot")]
    else:
        cands = []
    for c, rule in cands:
        if c.is_file():
            return rule, c
    return "fallback", None


def unique_by_name(ds, name):
    r = subprocess.run(["find", str(ds), "-name", name, "-print"],
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    hits = [Path(x) for x in r.stdout.decode("utf-8", "replace").splitlines() if x.strip()]
    return hits


def main():
    fix = "--fix" in sys.argv
    only = None
    if "--ds" in sys.argv:
        only = sys.argv[sys.argv.index("--ds") + 1]
    OUT.mkdir(parents=True, exist_ok=True)
    print("MODE=%s ONLY=%s" % ("FIX" if fix else "DRY-RUN", only))
    g = {"broken_before": 0, "fixed": 0, "unresolved": 0, "ambiguous": 0, "post_broken": 0}
    for ds in ROOTS:
        if only and only not in ds.name:
            continue
        if not ds.is_dir():
            print("SKIP missing ds: %s" % ds)
            continue
        t0 = time.time()
        broken = find_broken(ds)
        t_scan = time.time() - t0
        plan, unresolved = [], []
        for link in broken:
            rule, tgt = rule_target(link, ds)
            if tgt is None:
                unresolved.append((rule, link))
            else:
                plan.append((rule, link, tgt))
        for rule, link in unresolved:
            hits = unique_by_name(ds, link.name)
            if len(hits) == 1:
                plan.append(("generic", link, hits[0]))
            elif len(hits) > 1:
                plan.append(("ambiguous", link, None))
            else:
                plan.append(("unresolved", link, None))
        stat = {}
        for rule, link, tgt in plan:
            stat[rule] = stat.get(rule, 0) + 1
        print("\n##### %s  broken=%d  scan=%.1fs  rules=%s" % (ds.name, len(broken), t_scan, stat))
        g["broken_before"] += len(broken)
        g["unresolved"] += stat.get("unresolved", 0)
        g["ambiguous"] += stat.get("ambiguous", 0)
        man = OUT / ("%s.tsv" % ds.name)
        with open(man, "w") as fh:
            fh.write("link_relpath\told_target\trule\tnew_target\taction\n")
            for rule, link, tgt in plan:
                old = os.readlink(str(link))
                rel = os.path.relpath(str(link), str(ds))
                if tgt is None:
                    fh.write("%s\t%s\t%s\t\tUNRESOLVED\n" % (rel, old, rule))
                    print("   ! %s  rule=%s  old=%s" % (rel, rule, old))
                    continue
                newrel = os.path.relpath(str(tgt), str(link.parent))
                act = "DRY"
                if fix:
                    os.unlink(str(link))
                    os.symlink(newrel, str(link))
                    act = "FIXED"
                    g["fixed"] += 1
                fh.write("%s\t%s\t%s\t%s\t%s\n" % (rel, old, rule, newrel, act))
        print("   manifest -> %s" % man)
        if fix:
            after = find_broken(ds)
            g["post_broken"] += len(after)
            print("   post_check_broken=%d" % len(after))
    print("\n===== SUMMARY =====")
    for k, v in g.items():
        print("%s=%s" % (k, v))
    return 0


if __name__ == "__main__":
    sys.exit(main())
