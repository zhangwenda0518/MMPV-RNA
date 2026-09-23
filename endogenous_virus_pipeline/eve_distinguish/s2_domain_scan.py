#!/usr/bin/env python3
"""S2 功能通道 v4: 基因座(locus)架构 + 前病毒规模判定.

v3 的做法: architecture_score = 非重叠病毒基因座数, 判定里"loci>=2 或 span>=1500"
  就认为是强候选. 这其实是拿"完整性"当 EVE 证据, 方向是反的 —— RNA-seq 里真 DNA 病毒
  本来就拼不全(环状 7.5kb 基因组靠短读段只能拿到局部), "不全" 不等于 EVE.
v4 的关键倒置: 真正有判别力的是 **单一位点上 MP+CP+AP+RT+RH 齐全** (前病毒规模),
  且必须与 S1 的密码子退化**联合**才有意义:
     齐全 + 密码子退化 -> 整合的前病毒 (EVE)
     齐全 + 密码子完好 -> 感染中的真病毒 (要留的)
     不全 + 完好       -> 病毒碎片 (信息不足, 不能据此判 EVE)
  所以本模块只负责把架构事实测准, 不独自下 EVE 结论 (判定在 s3).

- MP/CP/AP 为病毒特有基因, 转座子没有 -> 架构由 locus 驱动
- TE flag: baits FALSE_* (LTR 逆转录转座子域) bitscore>40

用法: python3 s2_domain_scan.py <query.fasta> <panel_hits.tsv> <baits_hits.tsv> <out.tsv>
"""
import sys
from collections import defaultdict

CANON = ["MP", "CP", "AP", "RT", "RH"]
MIN_ALN, MIN_PID, MIN_BITS = 50, 25, 40
PROVIRUS_SPAN = 1500     # 单个基因座达到此跨度才够"前病毒"规模
PROVIRUS_NCOMP = 4       # 单基因座上至少 4/5 个规范组件


def merge_loci(hsps):
    """按 query 坐标合并重叠 HSP -> 基因座列表. 负链 HSP qstart>qend, 先归一化."""
    hsps = sorted(hsps, key=lambda h: h["qs"])
    loci = []
    for h in hsps:
        if loci and h["qs"] <= loci[-1]["qe"] + 30:  # 30nt 内视为同一基因座
            L = loci[-1]
            L["qe"] = max(L["qe"], h["qe"])
            for c in h["comps"]:
                if c not in L["comps"]:
                    L["comps"] += "+" + c
            L["bits"] = max(L["bits"], h["bits"])
            L["detail"].append(h)
        else:
            loci.append({"qs": h["qs"], "qe": h["qe"], "comps": h["comps"],
                         "bits": h["bits"], "detail": [h]})
    return loci


def main():
    fa, panel_tsv, baits_tsv, out_tsv = sys.argv[1:5]
    ids = []
    name = None
    for line in open(fa):
        if line.startswith(">"):
            if name:
                ids.append(name)
            name = line[1:].split()[0]
    if name:
        ids.append(name)

    # ---- panel 命中 -> 每条 HSP (组件标签来自面板条目, 但计分按 locus) ----
    contig_hsps = defaultdict(list)
    with open(panel_tsv) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or line.startswith("qseqid"):
                continue
            try:
                aln, pident, bits = int(p[3]), float(p[2]), float(p[9])
                qs, qe, sstart, send = int(p[4]), int(p[5]), int(p[6]), int(p[7])
            except ValueError:
                continue
            if aln < MIN_ALN or pident < MIN_PID:
                continue
            if abs(send - sstart) + 1 < MIN_ALN:
                continue
            comps = "+".join(c for c in p[1].split("|")[0].split("+") if c in CANON)
            if not comps:
                continue
            lo, hi = min(qs, qe), max(qs, qe)
            contig_hsps[p[0]].append({"qs": lo, "qe": hi, "comps": comps,
                                      "bits": bits, "pident": pident, "aln": aln})

    # ---- baits 命中 -> TE flag / cauli 同源 ----
    te_bits, cauli_bits, te_label, cauli_label = (defaultdict(float), defaultdict(float),
                                                  defaultdict(str), defaultdict(str))
    with open(baits_tsv) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or line.startswith("qseqid"):
                continue
            try:
                bits, pident = float(p[9]), float(p[2])
            except ValueError:
                continue
            if pident < MIN_PID or bits < MIN_BITS:
                continue
            sid, label = p[1], p[11] if len(p) > 11 else ""
            if sid.startswith("FALSE_"):
                if bits > te_bits[p[0]]:
                    te_bits[p[0]], te_label[p[0]] = bits, label
            elif bits > cauli_bits[p[0]]:
                cauli_bits[p[0]], cauli_label[p[0]] = bits, label

    def locus_ncomp(L):
        return len([c for c in L["comps"].split("+") if c in CANON])

    with open(out_tsv, "w") as out:
        out.write("contig_id\tarchitecture_score\tmax_locus_span\tloci_detail\tcomponents\t"
                  "evidence_detail\tte_flag\tte_bitscore\tte_label\tcauli_bitscore\tcauli_label\t"
                  "best_locus_span\tbest_locus_ncomp\tbest_locus_comps\tlocus_completeness\t"
                  "provirus_scale\tmulti_locus_arch\n")
        n_scale = 0
        for cid in ids:
            loci = merge_loci(contig_hsps.get(cid, []))
            score = min(len(loci), 5)
            comps_all = "+".join(dict.fromkeys(
                c for L in loci for c in L["comps"].split("+") if c in CANON))
            loci_detail = ";".join(
                f"L{i+1}:{L['qs']}-{L['qe']}({L['comps']},bit{L['bits']:.0f})"
                for i, L in enumerate(loci))
            detail = ";".join(
                f"{h['comps']}:{h['pident']:.0f}%/{h['aln']}aa"
                for L in loci for h in sorted(L["detail"], key=lambda x: -x["bits"])[:2])
            # v4: 以"单一基因座"为单位 — richest locus 决定前病毒规模
            if loci:
                best = max(loci, key=lambda L: (locus_ncomp(L), L["qe"] - L["qs"]))
                best_span = best["qe"] - best["qs"]
                best_n = locus_ncomp(best)
                best_comps = best["comps"]
            else:
                best_span, best_n, best_comps = 0, 0, "-"
            completeness = round(best_n / len(CANON), 3)
            provirus = (best_span >= PROVIRUS_SPAN and best_n >= PROVIRUS_NCOMP)
            if provirus:
                n_scale += 1
            multi = "TRUE" if sum(1 for L in loci if locus_ncomp(L) >= 1) >= 2 else "FALSE"
            te = "TRUE" if te_bits.get(cid, 0) > MIN_BITS else "FALSE"
            out.write(f"{cid}\t{score}\t{max((L['qe'] - L['qs'] for L in loci), default=0)}\t"
                      f"{loci_detail or '-'}\t{comps_all or '-'}\t{detail or '-'}\t{te}\t"
                      f"{te_bits.get(cid, 0):.0f}\t{te_label.get(cid, '-')[:60]}\t"
                      f"{cauli_bits.get(cid, 0):.0f}\t{cauli_label.get(cid, '-')[:60]}\t"
                      f"{best_span}\t{best_n}\t{best_comps}\t{completeness}\t"
                      f"{'TRUE' if provirus else 'FALSE'}\t{multi}\n")
    print(f"[done] s2(v4 single-locus): {len(ids)} contigs, "
          f"前病毒规模(单基因座>={PROVIRUS_SPAN}nt 且>={PROVIRUS_NCOMP}/5 组件) {n_scale} 条 -> {out_tsv}")


if __name__ == "__main__":
    main()
