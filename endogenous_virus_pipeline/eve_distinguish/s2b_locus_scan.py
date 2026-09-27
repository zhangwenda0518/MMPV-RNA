#!/usr/bin/env python3
"""S2b 位点归并 [判定逻辑最后变更: 模块 v5.1]: 把候选 contig 按寄主基因组 scaffold 聚合, 在位点层面算基因座架构.

为什么需要这一步: "单一位点上 MP+CP+AP+RT+RH 齐全" 里的"位点"是**寄主基因组位点**,
不是 contig. 一个感染中的病毒在 RNA-seq 里会被打成若干个 contig, 每个只带一两个基因;
而一个整合的前病毒在位点上是完整的, 只是被组装/注释切成碎片. 所以要把映射到同一个
寄主 scaffold 的候选 contig 合并起来再看架构.

v5 修订 (2026-09-23): 两个对判定可信度影响很大的问题, 都在这一版修掉:

  1. **未映射候选不得合并成伪位点** (旧版核心缺陷). 旧代码用
     `loci[(sample, host_scaffold)]` 归组, 而没有任何寄主命中的候选 host_scaffold
     是 "-", 于是同一样本下**所有**未映射候选被塌缩成一个位点 —— 若干各带一两个
     组件的 contig 并起来凑够 >=4 组件, 就伪造出"单位点齐全", 进而伪造 EVE_STRONG
     / virus_candidate. 实测一次 1027 候选的运行: 787 条无寄主命中, 484 条落进
     多人伪位点, 39 个 locus_full 里 27 个出自伪位点. 现在未映射候选每条**自成一位点**,
     arch 记为 `no_host_locus`, 不参与位点齐全判定 (s3 的 locus_full 只看真位点).
  2. **"自身物种"口径必须可见**. 旧版在样本没有映射表条目时静默退化成"最佳 bitscore
     所在物种"的跨物种启发式, 但同源强度和 reason 文本都还写着"自身寄主物种".
     现在逐候选记录 host_scope:
       own_species         样本在映射表内且该物种有比对 -> 可当自身寄主证据
       cross_species_guess 样本无映射条目(表里没有) -> 命中来自哪个物种就算哪个, 只是旁证
       none                无任何寄主命中 -> 位点通道无输入
     s3 的寄主否决 reason 会带上这个口径, 下游可按 host_scope 过滤.

同时这一步输出"候选 vs 寄主"的同源强度, 用来剔除把寄主基因当病毒的假阳性
(实测: 两条 27-29% 流行率的 DNA 候选就是寄主 HSP70, 靠 CDD 的 hsp70 域 PASS_VIRAL 混进来的).

输入 ( blastn outfmt 6 ):
  qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen slen
  (query = 候选 contig, subject = 该样本物种的 OneKP 转录组组装)

用法: python3 s2b_locus_scan.py <locus_blastn.tsv> <s2_domains.tsv> <s1_decay.tsv> <out.tsv> [host_map.tsv]
      python3 s2b_locus_scan.py --emit-header      # 只打印表头 (run_all.sh 的 --no-host 分支复用)
host_map.tsv: 三列 "样本号 <TAB> 4位代码 <TAB> 物种目录名" (run_all.sh 用的同一份).
  给了它, 寄主同源强度只统计**自身物种**的比对 —— 跨物种旁系同源会虚高, 不能算自身寄主.
输出: locus_architecture.tsv
"""
import argparse
import sys
from collections import defaultdict

CANON = ["MP", "CP", "AP", "RT", "RH"]
BLASTN_MIN_COLS = 14          # qseqid..slen; 少于这个数说明 outfmt 配错了
NOHOST_LABEL = "__no_host__"  # 未寄主映射候选的位点组名 (每组只剩自己)


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
    """寄主比对 -> 每 contig 每物种的同源强度 + 解析统计.

    返回 (out, stats):
      out   cid -> {code: {wpid, cov, aln, best_scaffold, best_bits, n_scafs}}
              wpid      比对长度加权平均一致性
              cov       该物种全部 scaffold 对 query 的区间并集覆盖率 (并集而非单条:
                        OneKP 转录组组装是碎的, 一个寄主基因常被打断在多个 scaffold 上,
                        单条覆盖率会严重低估)
              n_scafs   该物种下命中了几条不同 scaffold。位点归并只取最佳那条 scaffold,
                        这个数用来暴露"同源区实际被打断在多条 scaffold 上"的漏并风险.
      stats dict(n_rows, n_bad): 行数与列数不足被跳过的行数。n_bad>0 而解析行数为 0
                        时基本就是 outfmt 配错了, 调用方必须报错而不是让全部候选静默
                        变成"无寄主映射"。
    """
    agg = defaultdict(lambda: defaultdict(
        lambda: {"aln": 0, "pid_w": 0.0, "ivs": [], "qlen": 0,
                 "best": (0.0, 0, "", 0.0), "scafs": set()}))
    stats = {"n_rows": 0, "n_bad": 0}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            stats["n_rows"] += 1
            p = line.rstrip("\n").split("\t")
            if len(p) < BLASTN_MIN_COLS:
                stats["n_bad"] += 1
                continue
            q, s, pid, aln, bits = p[0], p[1], float(p[2]), int(p[3]), float(p[11])
            qs, qe, qlen = int(p[6]), int(p[7]), int(p[12])
            code = scaffold_code(s)
            a = agg[q][code]
            a["aln"] += aln
            a["pid_w"] += aln * pid
            a["ivs"].append((min(qs, qe), max(qs, qe)))
            a["qlen"] = qlen
            a["scafs"].add(s)
            if (bits, aln) > (a["best"][3], a["best"][1]):
                a["best"] = (pid, aln, s, bits)
    out = {}
    for q, per_code in agg.items():
        out[q] = {}
        for code, a in per_code.items():
            cov = union_len(sorted(a["ivs"])) / max(a["qlen"], 1)
            out[q][code] = {"wpid": round(a["pid_w"] / a["aln"], 2),
                            "cov": round(cov, 3),
                            "aln": a["aln"],
                            "best_scaffold": a["best"][2],
                            "best_bits": a["best"][3],
                            "n_scafs": len(a["scafs"])}
    return out, stats


HDR = ("contig_id\tsample\tsample_flag\thost_scaffold\thost_code\thost_scope\thost_wpid\t"
       "host_cov\thost_xeno\tn_host_scaffolds\tcomps_contig\tbest_locus_comps\tdecay_class\t"
       "te_flag\tlocus_key\tlocus_arch\tlocus_ncomp\tlocus_comps\tlocus_members\tlocus_hint\n")


def comps_of(s):
    return [c for c in (s or "").split("+") if c in CANON]


def parse_sample(cid):
    """候选头 -> (样本号, 标记)。保守序列按 <样本号>_clean_<NODE...> 命名; 认不出这个
    标记时(如 NC_013134.1 / contig_22919)退回整个头当样本号 —— 每条自成一组, 不会和
    别的候选混在一起, 并把 no_clean_marker 标出来, 免得静默按错口径判寄主同源。"""
    if "_clean_" in cid:
        return cid.split("_clean_", 1)[0], "clean"
    return cid, "no_clean_marker"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('blast_tsv', nargs='?', help='candidate vs OneKP blastn outfmt 6')
    ap.add_argument('s2_tsv', nargs='?', help='s2_domain_scan.py 产出')
    ap.add_argument('s1_tsv', nargs='?', help='s1_decay_scan.py 产出 (可为空)')
    ap.add_argument('out_tsv', nargs='?', help='输出 locus_architecture.tsv')
    ap.add_argument('host_map_tsv', nargs='?', default='',
                    help='host_map.tsv: 样本号 <TAB> 4位代码 <TAB> 物种目录名')
    ap.add_argument('--emit-header', action='store_true',
                    help='只打印输出表头并退出 (run_all.sh 的 --no-host 分支复用)')
    a = ap.parse_args()

    # 表头只写进出文件: s3 按首行当列名读这份表, 少了表头就会把第一条候选当
    # 表头吃掉、后面整行列错位。--emit-header 是给 run_all.sh 的 --no-host 分支
    # 用的 (那次不跑 s2b), 只在 stdout 打一份, 两头共用同一个常量, 不会不一致。
    if a.emit_header:
        sys.stdout.write(HDR)
        return 0
    if not all([a.blast_tsv, a.s2_tsv, a.out_tsv]):
        ap.error('需要 blast_tsv / s2_tsv / out_tsv (或改用 --emit-header)')

    sample2code = {}
    if a.host_map_tsv:
        with open(a.host_map_tsv) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[0] and not p[0].startswith("#"):
                    sample2code[p[0]] = p[1]

    host, st = parse_blastn(a.blast_tsv)
    if st["n_bad"] and not host:
        sys.exit("[s2b] %s 有 %d 行但一行都解析不出来 (每行至少 %d 列: "
                 "qseqid sseqid pident length mismatch gaps qstart qend sstart send "
                 "evalue bitscore qlen slen)。八成是 blastn 的 -outfmt 与这里不一致; "
                 "修好再跑, 否则全部候选都会被当成'无寄主映射', 位点通道静默失效。"
                 % (a.blast_tsv, st["n_rows"], BLASTN_MIN_COLS))

    s2 = {r["contig_id"]: r for r in read_tsv(a.s2_tsv)}
    s1 = {r["contig_id"]: r for r in read_tsv(a.s1_tsv)} if a.s1_tsv else {}

    # 逐 contig: 位点归属 + 寄主同源强度 (覆盖 s1 ∪ s2, 不留空洞)
    # 输出按 contig_id 排序: 同一份输入必须得到逐字节相同的表。旧版按 set 迭代序遍历,
    # 两次跑出来的行序不同, 没法 diff/复现比对。
    rows = {}
    for cid in sorted(set(s2) | set(s1)):
        d2 = s2.get(cid, {})
        per_code = host.get(cid, {})
        sample, sflag = parse_sample(cid)
        own = sample2code.get(sample)
        # 自身物种: 映射表优先。没有映射表条目时退化为"最佳单条命中的 bitscore 所在物种" ——
        # 跨物种旁系同源的比对更短更弱, 该启发式实测 20/23 正确; 但在旁系特别强的样本
        # (如 Hemerocallis HSP70 撞上 Borya 旁系) 会判错, 所以映射表仍是首选输入,
        # 且这种候选一律带 host_scope=cross_species_guess, 不当自身寄主证据用。
        if own and own in per_code:
            code, scope = own, "own_species"
        elif per_code:
            code = max(per_code, key=lambda c: per_code[c]["best_bits"])
            scope = "cross_species_guess"
        else:
            code, scope = "", "none"
        if code:
            h = per_code[code]
            wpid, cov, scaff, n_scafs = (h["wpid"], h["cov"],
                                         h["best_scaffold"], h["n_scafs"])
            xeno = max(((v["wpid"], v["aln"], c) for c, v in per_code.items() if c != code),
                       default=(0.0, 0, ""))
        else:
            wpid, cov, scaff, n_scafs, xeno = "", 0.0, "-", 0, (0.0, 0, "")
        rows[cid] = {
            "sample": sample,
            "sample_flag": sflag,
            "host_scaffold": scaff,
            "host_code": code,
            "host_scope": scope,
            "host_wpid": wpid,
            "host_cov": cov,
            "host_xeno": "%s/%dnt" % (xeno[2], xeno[1]) if xeno[2] else "-",
            "n_host_scaffolds": n_scafs,
            "comps_contig": set(comps_of(d2.get("components", "-"))),
            "best_locus_comps": set(comps_of(d2.get("best_locus_comps", "-"))),
            "decay": s1.get(cid, {}).get("decay_class", ""),
            "te": d2.get("te_flag", "FALSE"),
        }

    # ---- 位点级汇总 ----
    # 真位点 = (样本, 寄主 scaffold)。未寄主映射的候选**不归组**: 旧版把同样本的它们
    # 全部塌缩进 ("-") 一个位点, 若干各带一两个组件的 contig 并起来就能凑出"单位点
    # 齐全", 这是伪造 EVE 判定的直接原因。现在每条自成一位点, arch=no_host_locus。
    loci = defaultdict(list)
    for cid, r in rows.items():
        if r["host_scaffold"] in ("-", ""):
            loci[(NOHOST_LABEL, cid)].append(cid)
        else:
            loci[(r["sample"], r["host_scaffold"])].append(cid)

    for key, members in loci.items():
        nohost = key[0] == NOHOST_LABEL
        union = set()
        for m in members:
            union |= rows[m]["comps_contig"]
        n = len(union & set(CANON))
        decays = {rows[m]["decay"] for m in members}
        has_decay = any(d == "distributed_decay" for d in decays)
        has_break = any(d == "assembly_breakpoint" for d in decays)
        pol_only = bool(union) and not (union & {"MP", "CP", "AP"})
        if nohost:
            # 没有寄主位点就没有"位点级"证据: 只描述这条 contig 自己带了什么组件,
            # 不与其他候选合成架构, 更不许当 locus_full
            arch, hint = "no_host_locus", "no_host_mapping"
        elif n >= 4:
            arch = "locus_full"
        elif n >= 2:
            arch = "locus_partial"
        elif n == 1:
            arch = "locus_single_comp"
        else:
            arch = "locus_no_cauli"
        if not nohost:
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
        lkey = ("%s|nohost:%s" % (rows[members[0]]["sample"], members[0])
                if nohost else "%s|%s" % (key[0], key[1]))
        for m in members:
            rows[m].update({
                "locus_arch": arch,
                "locus_ncomp": n,
                "locus_comps": "+".join(c for c in CANON if c in union) or "-",
                "locus_members": len(members),
                "locus_hint": hint,
                "locus_key": lkey,
            })

    with open(a.out_tsv, "w") as out:
        out.write(HDR)
        for cid in rows:
            r = rows[cid]
            out.write(f"{cid}\t{r['sample']}\t{r['sample_flag']}\t{r.get('host_scaffold','-')}\t"
                      f"{r.get('host_code','')}\t{r.get('host_scope','none')}\t"
                      f"{r.get('host_wpid','')}\t{r.get('host_cov','')}\t{r.get('host_xeno','-')}\t"
                      f"{r.get('n_host_scaffolds',0)}\t"
                      f"{'+'.join(c for c in CANON if c in r.get('comps_contig',set())) or '-'}\t"
                      f"{'+'.join(c for c in CANON if c in r.get('best_locus_comps',set())) or '-'}\t"
                      f"{r.get('decay','')}\t{r.get('te','FALSE')}\t{r.get('locus_key','')}\t"
                      f"{r.get('locus_arch','unresolved')}\t{r.get('locus_ncomp','')}\t"
                      f"{r.get('locus_comps','-')}\t{r.get('locus_members','')}\t"
                      f"{r.get('locus_hint','no_host_assembly')}\n")

    # ---- 运行摘要: 位点通道到底覆盖了多少候选, 一眼可见 ----
    n_real = len({(r["sample"], r["host_scaffold"]) for r in rows.values()
                  if r["host_scaffold"] not in ("-", "")})
    n_nohost = sum(1 for r in rows.values() if r["locus_arch"] == "no_host_locus")
    n_full = sum(1 for r in rows.values() if r["locus_arch"] == "locus_full")
    sc = defaultdict(int)
    for r in rows.values():
        sc[r["host_scope"]] += 1
    print(f"[done] s2b -> {a.out_tsv}  (自身物种映射 {len(sample2code)} 条)")
    print(f"  候选 {len(rows)} 条: 真实位点 {n_real} 个 / 无寄主映射 {n_nohost} 条; "
          f"locus_full {n_full} 条")
    print(f"  寄主比对范围: 自身物种 {sc['own_species']} 条, "
          f"跨物种启发式 {sc['cross_species_guess']} 条, 无命中 {sc['none']} 条")
    if st["n_bad"]:
        print(f"[warn] blastn 表有 {st['n_bad']}/{st['n_rows']} 行列数不足 "
              f"({BLASTN_MIN_COLS}) 被跳过", file=sys.stderr)
    n_flag = sum(1 for r in rows.values() if r["sample_flag"] == "no_clean_marker")
    if n_flag:
        print(f"[warn] {n_flag} 条候选头没有 '_clean_' 标记, 样本号取整个头, 各自成组; "
              f"自身物种口径对它们不可用", file=sys.stderr)
    if sc["own_species"] * 2 < len(rows):
        print(f"[warn] 只有 {sc['own_species']}/{len(rows)} 条候选能按自身物种比: host_map "
              f"覆盖的样本太少, 寄主同源大半是跨物种旁证 (host_scope 列可分辨)",
              file=sys.stderr)
    if n_nohost * 2 > len(rows):
        print(f"[warn] {n_nohost}/{len(rows)} 条候选无任何寄主命中: 位点级架构对多数候选"
              f"不可用, s3 对这些只能给 EVE_suspect / *_review, 不会判 EVE_STRONG",
              file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
