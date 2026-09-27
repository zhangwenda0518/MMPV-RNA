#!/usr/bin/env python3
"""数据集级 v68(fill-fix1) 前后对照汇总。

用法:
  python3 verify_v68.py Alternaria <dataset_dir> [Aphis <dir> ...]

每数据集输出: combined/成品 的行数与 md5 (旧=备份, 新=线上)、gate 关键计数新旧对照、
成品 Genus 多词值计数新旧对照、以及 /tmp/refresh_v68/<label>/cmp_combined.txt 的方向行。
只读。
"""
import csv, hashlib, os, sys

STAMP = "20260916"
LOGROOT = "/tmp/refresh_v68"
KEY = ["rows_total", "family_filled", "genus_filled", "genus_ref_known",
       "dual_ref_conflict", "single_ref_conflict", "species_dual_ref_conflict",
       "species_single_ref_conflict"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def nrows(p):
    try:
        with open(p, "rb") as f:
            return sum(1 for _ in f) - 1
    except FileNotFoundError:
        return -1


def gate(p):
    d = {}
    try:
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                c = (r.get("check") or "").strip()
                if c:
                    d[c] = ((r.get("n") or "").strip(), (r.get("status") or "").strip())
    except FileNotFoundError:
        pass
    return d


def multiword(p):
    """返回 (多词行数, 不同值数)"""
    try:
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            n, s = 0, set()
            for r in csv.DictReader(f, delimiter="\t"):
                g = (r.get("Genus") or "").strip().strip('"').strip()
                if g and g.upper() != "NA" and " " in g:
                    n += 1
                    s.add(g)
            return n, len(s)
    except FileNotFoundError:
        return -1, -1


def show(label, d):
    INT = os.path.join(d, "05_Taxonomy", "Votus.integrated")
    CLS = os.path.join(d, "05_Taxonomy", "Votus.classed")
    P_NEW = os.path.join(INT, "final_integrated_classification.tsv")
    P_OLD = INT + ".bak_fillfix_" + STAMP + "/final_integrated_classification.tsv"
    C_NEW = os.path.join(CLS, "Votus_combined_taxonomy.tsv")
    C_OLD = os.path.join(CLS, "Votus_combined_taxonomy.tsv.bak_fillfix_" + STAMP)
    print("=" * 100)
    print("%s  (%s)" % (label, os.path.basename(d)))
    for tag, old, new in (("combined", C_OLD, C_NEW), ("成品", P_OLD, P_NEW)):
        ro, rn = nrows(old), nrows(new)
        mo = md5(old)[:12] if os.path.exists(old) else "NA"
        mn = md5(new)[:12] if os.path.exists(new) else "NA"
        print("  %-8s 行数 %6d -> %-6d (%+d)   md5 %s -> %s%s"
              % (tag, ro, rn, rn - ro, mo, mn, "  <== 相同" if mo == mn and mo != "NA" else ""))
    go = gate(INT + ".bak_fillfix_" + STAMP + "/taxonomy_gate_check.tsv")
    gn = gate(os.path.join(INT, "taxonomy_gate_check.tsv"))
    print("  gate     " + "  ".join(
        "%s %s->%s%s" % (k, go.get(k, ("-", ""))[0], gn.get(k, ("-", ""))[0],
                         "" if gn.get(k, ("-", ""))[1] == "PASS" else "(%s)" % gn.get(k, ("-", ""))[1])
        for k in ("dual_ref_conflict", "single_ref_conflict", "species_dual_ref_conflict", "species_single_ref_conflict")))
    print("  gate计数 " + "  ".join("%s %s->%s" % (k, go.get(k, ("-",))[0], gn.get(k, ("-",))[0])
                                  for k in ("rows_total", "family_filled", "genus_filled", "genus_ref_known")))
    po, pn = multiword(P_OLD), multiword(P_NEW)
    print("  Genus多词 行 %d -> %d ; 不同值 %d -> %d" % (po[0], pn[0], po[1], pn[1]))
    lg = os.path.join(LOGROOT, label, "cmp_combined.txt")
    for line in open(lg, encoding="utf-8", errors="replace"):
        if "改动行" in line or "单元格方向" in line or "各阶元" in line:
            print("  combined " + line.rstrip())


argv = sys.argv[1:]
for i in range(0, len(argv) - 1, 2):
    show(argv[i], argv[i + 1])
