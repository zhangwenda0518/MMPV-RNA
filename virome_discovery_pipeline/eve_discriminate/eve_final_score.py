#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_final_score.py — 综合多路证据，输出候选最终分层

整合的五路证据
--------------
1. **同源参照** (ref_classify): 与 41,933 个已知植物 EVE 元件的同源性
2. **转录组嵌合** (chimera): 两端宿主 + 中间病毒 = 整合旁证
3. **位点架构** (locus_architecture): 同一宿主 scaffold 上组件是否齐全
4. **多拷贝家族** (mmseqs): 是否有家族扩增
5. **BAM 覆盖** (bam_coverage): 是否有活跃病毒样覆盖

评分逻辑
--------
不给"是/否"的硬判定，而是给**证据计数 + 分层**:

  tier 1 (内源强证据): >=2 条独立的内源证据
  tier 2 (内源弱证据): 1 条内源证据
  tier 3 (无证据)    : 两侧都无
  tier 4 (偏外源)    : 有外源证据且无内源证据

**重要**: 这是"证据汇总"，不是"概率"。每个 tier 都附上具体是哪几条证据，
下游可自行按需过滤。

输入可缺省
----------
`--ref-classify / --chimera / --locus-arch` 都不是必填: 没给的通道按"无证据"计。
这样三种场景都能先出表 —— 有参考基因组 (三路齐全)、只有转录组 (三路全空,
结论请看 02_classify 的序列判别)、或只跑到一半。HTML 版见
`03_report/eve_final_report.py`。
"""
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

HDR = ("contig_id\tfinal_tier\tn_endo_evidence\tn_exo_evidence\t"
       "ref_class\teve_bits\tvirus_bits\tchimera_call\tchimera_mid_gap\t"
       "locus_arch\tlocus_ncomp\tmulticopy_size\tbam_mean_depth\t"
       "evidence_list\n")


def load_tsv(path, key="contig_id"):
    d = {}
    if not path or not Path(path).is_file():
        return d
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get(key):
                d[r[key]] = r
    return d


def load_multicopy(cluster_tsv):
    """mmseqs cluster -> {member: cluster_size}。"""
    size = Counter()
    with open(cluster_tsv, encoding="utf-8") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2:
                size[c[0]] += 1
    return size


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-classify", default="",
                    help="双库分类表; 不给则该通道记为空 (无参考基因组时常见)")
    ap.add_argument("--chimera", default="",
                    help="转录本嵌合表; 不给则该通道记为空")
    ap.add_argument("--locus-arch", default="",
                    help="位点架构表; 不给则该通道记为空")
    ap.add_argument("--multicopy-clusters", default="")
    ap.add_argument("--bam-coverage", default="")
    ap.add_argument("--query-fa", required=True,
                    help="候选 fasta (定义报告要列哪些候选, 不参与打分)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    if not any([a.ref_classify, a.chimera, a.locus_arch]):
        print("[warn] 三路内源通道都没给 (ref/chimera/locus): 所有候选都会落到 "
              "tier3_no_evidence。若要按序列判别出结论, 请另看 02_classify 的 "
              "ML / Δbit-score 报告。", file=sys.stderr)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ids = []
    with open(a.query_fa, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].strip()
                if h:
                    ids.append(h.split()[0])

    ref = load_tsv(a.ref_classify)
    chi = load_tsv(a.chimera)
    arch = load_tsv(a.locus_arch)
    multi = load_multicopy(a.multicopy_clusters) if a.multicopy_clusters else {}

    rows = []
    tiers = Counter()
    for cid in ids:
        r = ref.get(cid, {})
        c = chi.get(cid, {})
        g = arch.get(cid, {})
        endo, exo = [], []

        # 证据 1: 同源参照
        rcls = r.get("ref_class", "")
        if rcls == "endogenous_supported":
            endo.append(f"命中已知EVE元件(bits {r.get('eve_bits','')}, "
                        f"{r.get('eve_pident','')}%)")
        elif rcls == "exogenous_likely":
            exo.append(f"命中真实病毒(bits {r.get('virus_bits','')})")
        elif rcls == "both_families":
            endo.append("两侧均命中(家族内外兼有)")

        # 证据 2: 嵌合
        ccall = c.get("chimera_call", "")
        if ccall == "chimera_supported":
            endo.append(f"转录组嵌合(两端宿主, 间隔 {c.get('mid_gap','')}nt)")

        # 证据 3: 位点架构
        arch_v = g.get("locus_arch", "")
        if arch_v == "locus_full":
            endo.append(f"位点组件齐全({g.get('locus_ncomp','')}/5)")

        # 证据 4: 多拷贝
        msize = multi.get(cid, 0)
        if msize >= 3:
            endo.append(f"多拷贝家族({msize} 成员)")

        # 证据 5: BAM 覆盖 (仅当有数据时)
        # 注意: bam 是按样本聚合的, 这里只做提示

        n_endo, n_exo = len(endo), len(exo)
        if n_endo >= 2:
            tier = "tier1_endogenous_strong"
        elif n_endo == 1:
            tier = "tier2_endogenous_weak"
        elif n_exo >= 1:
            tier = "tier4_exogenous"
        else:
            tier = "tier3_no_evidence"
        tiers[tier] += 1

        rows.append((
            cid, tier, n_endo, n_exo,
            rcls, r.get("eve_bits", ""), r.get("virus_bits", ""),
            ccall, c.get("mid_gap", ""),
            arch_v, g.get("locus_ncomp", ""),
            msize, "",
            "; ".join(endo + exo) or "-",
        ))

    tsv = out / "final_tiers.tsv"
    with open(tsv, "w", encoding="utf-8", newline="\n") as fo:
        fo.write(HDR)
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")

    print(f"[done] {tsv}  ({len(rows)} 条)")
    print("\n分层分布:")
    for k, v in tiers.most_common():
        print(f"  {k:28s} {v:5d}  ({v/len(rows)*100:.1f}%)")

    # 各 tier 的证据构成
    print("\n各 tier 的证据组合:")
    for tier in ("tier1_endogenous_strong", "tier2_endogenous_weak",
                 "tier4_exogenous"):
        sub = [r for r in rows if r[1] == tier]
        if not sub:
            continue
        combos = Counter()
        for r in sub:
            ev = r[13]
            kinds = []
            for k in ("命中已知EVE元件", "转录组嵌合", "位点组件齐全",
                      "多拷贝家族", "命中真实病毒"):
                if k in ev:
                    kinds.append(k)
            combos["+".join(kinds)] += 1
        print(f"  {tier} (n={len(sub)}):")
        for k, v in combos.most_common(5):
            print(f"      {v:4d}  {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
