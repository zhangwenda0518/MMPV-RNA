#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读探针：量化 R 产物中「元信息列与实际分类列口径不一致」的问题。
  A. *_agree 列在闸门之前算：统计「值已被置空但 agree 仍声称有工具支持」的格数
  B. completeness / confidence 也在闸门之前算：与实际非空阶元数对比
  C. cascade 残留的科种不匹配 30 条是什么
  D. Kingdom 层是否有占位值/后缀异常
用法: python3 probe_stale.py P1.tsv [P2.tsv ...]
"""
import csv, io, re, sys, os
from collections import defaultdict, Counter

VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
TAX = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
NA_TOKENS = {"", "na", "n/a", "nan", "null", "none", "undefined"}


def norm(v):
    v = (v or "").strip()
    return "" if v.lower() in NA_TOKENS else v


def load(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        raw = f.read()
    lines = [ln for ln in raw.splitlines() if not ln.startswith("#")]
    rd = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    return list(rd)


def load_vmr():
    out = defaultdict(set)
    with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t", quotechar='"')
        hdr = next(rd)
        idx = {n: hdr.index(n) for n in TAX}
        for row in rd:
            if len(row) <= idx["Species"]:
                continue
            chain = [(row[idx[r]] or "").strip() for r in TAX]
            if not chain[5]:
                continue
            for i, rk in enumerate(TAX):
                if chain[i]:
                    out[(rk, chain[i].lower())].add(chain[5])
    return out


def n_ag(agstr):
    m = re.match(r"\s*(\d+)\s*/", agstr or "")
    return int(m.group(1)) if m else 0


for path in sys.argv[1:]:
    if not os.path.isfile(path):
        print("缺失 %s" % path)
        continue
    rows = load(path)
    print("\n" + "=" * 78)
    print("FILE %s  (%d 行)" % (path, len(rows)))

    print("\n[A] 值已置空但 *_agree 仍声称有支持（说明 agree 是闸门前口径）")
    for r in TAX:
        a = r + "_agree"
        if a not in (rows[0].keys() if rows else []):
            continue
        bad = sum(1 for x in rows if not norm(x.get(r)) and n_ag(x.get(a)) > 0)
        if bad:
            print("    %-9s %5d 行" % (r, bad))

    print("\n[B] completeness / confidence 与实际非空阶元数不符")
    mism = Counter()
    for x in rows:
        c = x.get("completeness")
        try:
            c = int(float(c))
        except Exception:
            continue
        real = sum(1 for r in TAX if norm(x.get(r)))
        if c != real:
            mism[(c, real)] += 1
    tot = sum(mism.values())
    print("    不符行数 %d / %d" % (tot, len(rows)))
    for (c, real), n in mism.most_common(6):
        print("      completeness=%d 实际=%d : %d 行" % (c, real, n))

    print("\n[C] Kingdom 层取值分布（前 8）")
    kc = Counter(norm(x.get("Kingdom")) or "<空>" for x in rows)
    for k, n in kc.most_common(8):
        print("      %5d  %s" % (n, k))
    odd = [k for k in kc if k not in ("<空>",) and not k.endswith("virae")]
    print("    不以 -virae 结尾的 Kingdom 值: %s" % (odd[:8] if odd else "无"))

vmr = load_vmr()
print("\n" + "=" * 78)
print("[D] 科种不匹配残留明细（VMR 口径，最多 12 条）")
shown = 0
for path in sys.argv[1:]:
    if not os.path.isfile(path):
        continue
    for x in load(path):
        fam, sp = norm(x.get("Family")), norm(x.get("Species"))
        if not (fam and sp):
            continue
        f = vmr.get(("Species", sp.lower()))
        if f and fam not in f:
            print("      %-30s %-38s VMR科=%s" % (fam[:30], sp[:38], sorted(f)[:3]))
            shown += 1
            if shown >= 12:
                break
    if shown >= 12:
        break
