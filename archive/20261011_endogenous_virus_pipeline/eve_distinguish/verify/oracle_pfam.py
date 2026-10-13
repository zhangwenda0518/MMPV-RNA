#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oracle_pfam.py — 组件判定的独立 oracle (hmmscan/Pfam)
=====================================================
问题: `components` 列的 MP/CP/AP/RT/RH 来自**面板条目名**, `credit_components` 只是在
条目名的基础上加了"独立区间 + 长度封顶"两道闸。整条链路上没有任何一步真的看过这段
query 序列里有没有对应的蛋白结构域 —— 这是在用"参考序列的注释"推断"query 的注释"。

本脚本用完全独立的第二套方法复核: 把 query 的同源区按 HSP 读框翻译成肽段, 跑 hmmscan
对 Pfam-A, 再按结构域家族判断这段序列到底有没有 RT / RNase H / 天冬氨酸蛋白酶。

覆盖范围与已知局限（写清楚，别让 oracle 显得比它实际能证明的更强）:
  * RT  -> PF00078 RVT_1 / PF07727 RVT_2 / PF13456 RVT_3
  * RH  -> PF00075 RNase_H
  * AP  -> PF00077 RVP  (CANON 的 "AP" 是 Aspartic protease, 面板条目名 Short=PR)
  * MP / CP **没有可靠的单家族映射** (运动蛋白/外壳蛋白家族繁多), 本 oracle 不判定,
    在输出里记 "n/a" —— 缺判据不等于判据通过。

用法:
  python3 oracle_pfam.py emit --run-dir REAL -o peptides.fa
  hmmscan --domtblout pfam.domtbl -E 1e-3 --cpu 32 Pfam-A.hmm peptides.fa > /dev/null
  python3 oracle_pfam.py compare --s2 REAL/s2_domains.tsv --domtbl pfam.domtbl \
      --json oracle.json
"""
import argparse
import collections
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import null_model as nm            # noqa: E402  复用 read_hits / tsv_load
import s1_decay_scan as s1         # noqa: E402  read_fasta 在这里

# Pfam 家族 -> CANON 组件. 只收映射无歧义的三种; MP/CP 记 n/a.
PFAM2COMP = {
    "PF00078": "RT", "PF07727": "RT", "PF13456": "RT",
    "PF00075": "RH",
    "PF00077": "AP",
}
ORACLE_COVERS = {"RT", "RH", "AP"}

CODON = {}
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _b1 in enumerate(_BASES):
    for _j, _b2 in enumerate(_BASES):
        for _k, _b3 in enumerate(_BASES):
            CODON[_b1 + _b2 + _b3] = _AAS[_i * 16 + _j * 4 + _k]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def translate(nt):
    nt = nt.upper().replace("U", "T")
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def hsp_peptides(run_dir, out_fa):
    """每条候选 -> 同源区按读框翻译出的肽段 (与 s2 的 HSP 过滤门槛一致)."""
    seqs = s1.read_fasta(os.path.join(run_dir, "q.fa"))
    hits = nm.read_hits(os.path.join(run_dir, "panel_hits.tsv"))
    import s2_domain_scan as s2
    n_pep = 0
    with open(out_fa, "w", encoding="utf-8") as fo:
        for cid, seq in seqs.items():
            seen = set()
            for h in hits.get(cid, []):
                if h["length"] < s2.MIN_ALN or h["pident"] < s2.MIN_PID:
                    continue
                qs, qe = min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])
                if (qs, qe) in seen:
                    continue
                seen.add((qs, qe))
                if h["frame"].startswith("-"):
                    sub = seq[qs - 1:qe]
                    sub = sub.translate(COMP)[::-1]
                else:
                    sub = seq[qs - 1:qe]
                # 按 HSP 的读框偏移对齐
                off = abs(int(h["frame"])) - 1
                pep = translate(sub[off:])
                # 去终止符, 丢过短的
                pep = pep.split("*")[0] if "*" in pep else pep
                if len(pep) < 30:
                    continue
                fo.write(">%s|%d-%d|%s\n" % (cid, qs, qe, h["frame"]))
                for i in range(0, len(pep), 60):
                    fo.write(pep[i:i + 60] + "\n")
                n_pep += 1
    print("[oracle] 写出 %d 条同源区肽段 -> %s" % (n_pep, out_fa))
    return n_pep


def parse_domtbl(path):
    """hmmscan --domtblout -> {contig: {组件}}."""
    per = collections.defaultdict(set)
    accs = collections.Counter()
    if not os.path.exists(path):
        raise SystemExit("找不到 domtblout: %s" % path)
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 23:
                continue
            # domtblout 列: target_name accession tlen query_name qlen full_E ...
            # 注意 p[3] 是**完整 query 名**, 与 read_fasta 的键同源; 不要再 split 掉,
            # 否则键对不上, 对照脚本会给出全 0 的假阴性 (已踩过)。
            acc = p[1].split(".")[0]
            accs[acc] += 1
            comp = PFAM2COMP.get(acc)
            if comp:
                per[p[3]].add(comp)
    return per, accs


def control(panel_fa, domtbl):
    """阳性对照 —— **跑 compare 之前必须先看这个**.

    拿面板参考蛋白当阳性样本: 它们条目名里就写着 RT/RH/AP, hmmscan 理应在它们身上看到
    对应家族。检出率低说明 oracle 自身灵敏度不够, 这时 compare 报出的"未被证实"绝大
    多数是 oracle 假阴性, **不能**当成"s2 过度归属组件"的证据。

    实测 (Pfam-A v35, E<=1e-3, 2026-09-24): RT 50% (30/60), RH 0% (0/6), AP 0% (0/15)
    —— 即这套 oracle 目前**不足以裁决**组件判定, 其 compare 结果只能记为"未定论"。
    要让它可用: 换 Caulimoviridae 专属家族、改用 --cut_ga (gathering threshold)、
    或改用 CDD rpsblast。
    """
    seqs = s1.read_fasta(panel_fa)
    hits = collections.defaultdict(set)
    for line in open(domtbl, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) < 23:
            continue
        comp = PFAM2COMP.get(p[1].split(".")[0])
        if comp:
            hits[p[3]].add(comp)
    stat = collections.defaultdict(lambda: [0, 0])
    for name in seqs:
        labeled = set(c for c in name.split("|")[0].split("+")
                      if c in ("MP", "CP", "AP", "RT", "RH"))
        got = hits.get(name, set())
        for c in labeled & ORACLE_COVERS:
            stat[c][1] += 1
            if c in got:
                stat[c][0] += 1
    print("=" * 74)
    print("oracle 阳性对照: 有该标注的面板参考蛋白里, Pfam 能检出的比例")
    print("=" * 74)
    weak = []
    for c in sorted(ORACLE_COVERS):
        ok, tot = stat[c]
        rate = (100.0 * ok / tot) if tot else float("nan")
        print("  %-4s %3d / %3d = %s" % (c, ok, tot, ("%.0f%%" % rate) if tot else "n/a"))
        if tot and rate < 80:
            weak.append(c)
    if weak:
        print("")
        print("!! 对照检出率低于 80%% 的组件: %s" % ", ".join(sorted(weak)))
        print("   -> oracle 灵敏度不足, **不要**把 compare 的'未被证实'当证据。")
    else:
        print("")
        print("对照全部通过, compare 的结论可用。")
    return {c: list(stat[c]) for c in sorted(ORACLE_COVERS)}


def contig_of(qname):
    """候选肽段名 `<contig>|<qs>-<qe>|<frame>` -> contig id; 面板条目的名字不匹配则 None.

    parse_domtbl 一律以**完整 query 名**做键 (对面板条目是对的, 它们的名字里也有 `|`,
    split 会切错); 候选侧要用 contig 聚合, 所以在这里显式剥掉后缀。两种命名混用同一
    个字典时, 只有本函数能把它们分开。
    """
    m = re.match(r"^(.*)\|\d+-\d+\|[+-]?\d+$", qname)
    return m.group(1) if m else None


def compare(s2_path, domtbl, json_out):
    rows, _ = nm.tsv_load(s2_path)
    per_name, accs = parse_domtbl(domtbl)
    oracle = collections.defaultdict(set)
    for qname, comps in per_name.items():
        cid = contig_of(qname)
        if cid:
            oracle[cid] |= comps
    n = len(rows)
    # 逐候选: s2 声称的组件 vs oracle 实测的组件 (只在 oracle 覆盖的 3 种上比)
    agree = disagree_miss = disagree_extra = 0
    both = neither = 0
    details = []
    for cid, r in rows.items():
        claim = set(c for c in (r.get("components") or "-").split("+") if c)
        o = oracle.get(cid, set())
        for comp in sorted(ORACLE_COVERS):
            c_in, o_in = comp in claim, comp in o
            if c_in and o_in:
                agree += 1
            elif c_in and not o_in:
                disagree_miss += 1        # s2 说有, oracle 没看到
            elif not c_in and o_in:
                disagree_extra += 1       # oracle 看到, s2 没记
            else:
                neither += 1
        if (claim & ORACLE_COVERS) or o:
            details.append({"contig": cid, "claim": sorted(claim),
                            "oracle": sorted(o)})
    tot = agree + disagree_miss + disagree_extra + neither
    print("=" * 74)
    print("组件判定独立 oracle (hmmscan vs Pfam-A) —— 只覆盖 RT / RH / AP")
    print("=" * 74)
    print("候选 %d 条; 逐候选 × 逐个(已覆盖组件) 共 %d 个判定点" % (n, tot))
    print("  双方都给      : %5d  (%.1f%%)" % (agree, 100.0 * agree / max(tot, 1)))
    print("  s2 说有/oracle无: %5d  (%.1f%%)  <- 面板条目名过度归属的主要嫌疑"
          % (disagree_miss, 100.0 * disagree_miss / max(tot, 1)))
    print("  oracle有/s2无  : %5d  (%.1f%%)  <- 漏记"
          % (disagree_extra, 100.0 * disagree_extra / max(tot, 1)))
    print("  双方都无      : %5d  (%.1f%%)" % (neither, 100.0 * neither / max(tot, 1)))
    cons = agree / max(agree + disagree_miss, 1)
    print("\n阳性一致率 (s2 声称里有多少被 oracle 证实) = %.1f%%" % (100 * cons))
    print("Pfam 命中最多的家族: %s" % ", ".join(
        "%s×%d" % (a, c) for a, c in accs.most_common(8)))
    print("\n注意: MP / CP 无可靠单家族映射, 本 oracle 不判定, 不计入上表。")
    res = {"n_candidates": n, "agree": agree, "claim_only": disagree_miss,
           "oracle_only": disagree_extra, "neither": neither,
           "positive_concordance": round(cons, 4),
           "top_pfam": dict(accs.most_common(20))}
    if json_out:
        with open(json_out, "w", encoding="utf-8") as f:
            json.dump({"summary": res, "details": details}, f,
                      ensure_ascii=False, indent=1)
        print("[oracle] -> %s" % json_out)
    return res


def main():
    ap = argparse.ArgumentParser(description="组件判定独立 oracle")
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("emit")
    e.add_argument("--run-dir", required=True)
    e.add_argument("-o", "--out", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--s2", required=True)
    c.add_argument("--domtbl", required=True)
    c.add_argument("--json", default="", dest="json_out")
    k = sub.add_parser("control", help="阳性对照: 先确认 oracle 自己看得见这些结构域")
    k.add_argument("--panel", required=True, help="面板参考蛋白 fasta")
    k.add_argument("--domtbl", required=True, help="面板蛋白的 hmmscan domtblout")
    a = ap.parse_args()
    if a.cmd == "emit":
        hsp_peptides(a.run_dir, a.out)
    elif a.cmd == "control":
        control(a.panel, a.domtbl)
    else:
        compare(a.s2, a.domtbl, a.json_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
