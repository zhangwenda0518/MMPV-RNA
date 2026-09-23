#!/usr/bin/env python3
"""S2b 位点归并: 把候选 contig 按寄主基因组 scaffold 聚合, 在位点层面算基因座架构.

为什么需要这一步: "单一位点上 MP+CP+AP+RT+RH 齐全" 里的"位点"是**寄主基因组位点**,
不是 contig. 一个感染中的病毒在 RNA-seq 里会被打成若干个 contig, 每个只带一两个基因;
而一个整合的前病毒在位点上是完整的, 只是被组装/注释切成碎片. 所以要把映射到同一个
寄主 scaffold 的候选 contig 合并起来再看架构.

同时这一步输出"候选 vs 自身寄主"的同源强度, 用来剔除把寄主基因当病毒的假阳性
(实测: 两条 27-29% 流行率的 DNA 候选就是寄主 HSP70, 靠 CDD 的 hsp70 域 PASS_VIRAL 混进来的).

输入 ( blastn outfmt 6 ):
  qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen slen
  (query = 候选 contig, subject = 该样本物种的 OneKP 转录组组装)

用法: python3 s2b_locus_scan.py <locus_blastn.tsv> <s2_domains.tsv> <s1_decay.tsv> <out.tsv> [host_map.tsv]
host_map.tsv: 三列 "样本号 <TAB> 4位代码 <TAB> 物种目录名" (run_all.sh 用的同一份).
  给了它, 寄主同源强度只统计**自身物种**的比对 —— 跨物种旁系同源会虚高, 不能算自身寄主.
输出: locus_architecture.tsv
"""
import sys
import csv
from collections import defaultdict

CANON = ["MP", "CP", "AP", "RT", "RH"]


def read_tsv(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(hdr, line.rstrip("\n").split("\t")))


def scaffold_code(s):
    """scaffold-IAJW-2026239-Amentotaxus_argotaenia -> IAJW (OneKP 四位物种代码)."""
    p = s.split("-")
    return p[1] if len(p) >= 4 and p[0] == "scaffold" else s


def union_len(ivs):
    """已排序区间的并集总长."""
    if not ivs:
        return 0
    tot, cs, ce = 0, ivs[0][0], ivs[0][1]
    for s, e in ivs[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            tot += ce - cs
            cs, ce = s, e
    return tot + (ce - cs)


def parse_blastn(path):
    """寄主比对 -> 每 contig 每物种的同源强度.

    返回 cid -> {code: (wpid, cov, aln_total, best_scaffold, best_key)}
      wpid  比对长度加权平均一致性
      cov   该物种全部 scaffold 对 query 的区间并集覆盖率 (并集而非单条: OneKP 转录组
            组装是碎的, 一个寄主基因常被打断在多个 scaffold 上, 单条覆盖率会严重低估)
    """
    agg = defaultdict(lambda: defaultdict(
        lambda: {"aln": 0, "pid_w": 0.0, "ivs": [], "qlen": 0, "best": (0.0, 0, "", 0.0)}))
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 14:
                continue
            q, s, pid, aln, bits = p[0], p[1], float(p[2]), int(p[3]), float(p[11])
            qs, qe, qlen = int(p[6]), int(p[7]), int(p[12])
            code = scaffold_code(s)
            a = agg[q][code]
            a["aln"] += aln
            a["pid_w"] += aln * pid
            a["ivs"].append((min(qs, qe), max(qs, qe)))
            a["qlen"] = qlen
            if (bits, aln) > (a["best"][3], a["best"][1]):
                a["best"] = (pid, aln, s, bits)
    out = {}
    for q, per_code in agg.items():
        out[q] = {}
        for code, a in per_code.items():
            cov = union_len(sorted(a["ivs"])) / max(a["qlen"], 1)
            out[q][code] = (round(a["pid_w"] / a["aln"], 2), round(cov, 3),
                            a["aln"], a["best"][2], a["best"][3])
    return out


def comps_of(s):
    return [c for c in (s or "").split("+") if c in CANON]


def main():
    argv = sys.argv[1:]
    if len(argv) < 4:
        print(__doc__)
        sys.exit(1)
    blast_tsv, s2_tsv, s1_tsv, out_tsv = argv[:4]
    host_map_tsv = argv[4] if len(argv) > 4 else ""

    sample2code = {}
    if host_map_tsv:
        with open(host_map_tsv) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[0] and not p[0].startswith("#"):
                    sample2code[p[0]] = p[1]

    host = parse_blastn(blast_tsv)
    s2 = {r["contig_id"]: r for r in read_tsv(s2_tsv)}
    s1 = {r["contig_id"]: r for r in read_tsv(s1_tsv)} if s1_tsv else {}

    # 逐 contig: 位点归属 + 自身寄主同源强度 (覆盖 s1 ∪ s2, 不留空洞)
    rows = {}
    for cid in set(s2) | set(s1):
        d2 = s2.get(cid, {})
        per_code = host.get(cid, {})
        sample = cid.split("_clean_")[0]
        own = sample2code.get(sample)
        # 自身物种: 映射表优先. 没有映射表时退化为"最佳单条命中的 bitscore 所在物种" ——
        # 跨物种旁系同源的比对更短更弱, 该启发式实测 20/23 正确; 但在旁系特别强的样本
        # (如 Hemerocallis HSP70 撞上 Borya 旁系) 会判错, 所以映射表仍是首选输入.
        if own and own in per_code:
            code = own
        elif per_code:
            code = max(per_code, key=lambda c: per_code[c][4])
        else:
            code = ""
        if code:
            wpid, cov, _, scaff, _ = per_code[code]
            xeno = max(((v[0], v[2], c) for c, v in per_code.items() if c != code),
                       default=(0.0, 0, ""))
        else:
            wpid, cov, scaff, xeno = "", 0.0, "-", (0.0, 0, "")
        rows[cid] = {
            "sample": sample,
            "host_scaffold": scaff,
            "host_code": code,
            "host_wpid": wpid,
            "host_cov": cov,
            "host_xeno": "%s/%dnt" % (xeno[2], xeno[1]) if xeno[2] else "-",
            "comps_contig": set(comps_of(d2.get("components", "-"))),
            "best_locus_comps": set(comps_of(d2.get("best_locus_comps", "-"))),
            "decay": s1.get(cid, {}).get("decay_class", ""),
            "te": d2.get("te_flag", "FALSE"),
        }

    # 位点级 (sample x scaffold) 汇总
    loci = defaultdict(list)
    for cid, r in rows.items():
        loci[(r["sample"], r["host_scaffold"])].append(cid)

    for key, members in loci.items():
        union = set()
        for m in members:
            union |= rows[m]["comps_contig"]
        n = len(union & set(CANON))
        decays = {rows[m]["decay"] for m in members}
        has_decay = any(d == "distributed_decay" for d in decays)
        has_break = any(d == "assembly_breakpoint" for d in decays)
        pol_only = bool(union) and not (union & {"MP", "CP", "AP"})
        if n >= 4:
            arch = "locus_full"
        elif n >= 2:
            arch = "locus_partial"
        elif n == 1:
            arch = "locus_single_comp"
        else:
            arch = "locus_no_cauli"
        if arch == "locus_full" and has_decay:
            hint = "provirus_integrated_signature"   # 齐全 + 密码子退化
        elif arch == "locus_full" and has_break:
            hint = "locus_full_but_assembled_breakpoint"
        elif arch == "locus_full":
            hint = "locus_full_intact"              # 齐全 + 完好 = 感染中的真病毒
        elif pol_only:
            hint = "pol_only_te_like"
        elif arch == "locus_no_cauli":
            hint = "host_mapped_no_cauli_domain"
        else:
            hint = "insufficient_info"
        for m in members:
            rows[m]["locus_arch"] = arch
            rows[m]["locus_ncomp"] = n
            rows[m]["locus_comps"] = "+".join(c for c in CANON if c in union) or "-"
            rows[m]["locus_members"] = len(members)
            rows[m]["locus_hint"] = hint

    with open(out_tsv, "w") as out:
        out.write("contig_id\tsample\thost_scaffold\thost_code\thost_wpid\thost_cov\thost_xeno\t"
                  "comps_contig\tbest_locus_comps\tdecay_class\tte_flag\tlocus_arch\tlocus_ncomp\t"
                  "locus_comps\tlocus_members\tlocus_hint\n")
        for cid in rows:
            r = rows[cid]
            out.write(f"{cid}\t{r['sample']}\t{r.get('host_scaffold','-')}\t{r.get('host_code','')}\t"
                      f"{r.get('host_wpid','')}\t{r.get('host_cov','')}\t{r.get('host_xeno','-')}\t"
                      f"{'+'.join(c for c in CANON if c in r.get('comps_contig',set())) or '-'}\t"
                      f"{'+'.join(c for c in CANON if c in r.get('best_locus_comps',set())) or '-'}\t"
                      f"{r.get('decay','')}\t{r.get('te','FALSE')}\t"
                      f"{r.get('locus_arch','unresolved')}\t{r.get('locus_ncomp','')}\t"
                      f"{r.get('locus_comps','-')}\t{r.get('locus_members','')}\t"
                      f"{r.get('locus_hint','no_host_assembly')}\n")
    print(f"[done] s2b -> {out_tsv}  (自身物种映射 {len(sample2code)} 条)")


if __name__ == "__main__":
    main()
