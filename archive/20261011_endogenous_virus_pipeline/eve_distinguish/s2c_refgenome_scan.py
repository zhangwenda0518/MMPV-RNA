#!/usr/bin/env python3
"""S2c 参考基因组通道 [判定逻辑最后变更: 模块 v6.3]: 候选 contig blastn 回**自身物种的
染色体级参考基因组**, 用整合信号区分内源 (EVE) 与真病毒.

与 S2b 的分工: S2b 比对的是 OneKP **转录组组装** (碎、无染色体上下文、只有映射表里的
物种), 回答"候选是不是寄主基因 / 位点上组件齐不齐"; S2c 比对的是**现成的参考基因组**,
回答"这条候选在寄主基因组里有没有整合证据". 两条证据正交, S2c 只在给过
`-G/--genome-manifest` 时运行; 没有参考基因组的物种自动跳过 (记 no_ref_genome),
**跳过不是"真病毒"的证据** —— 参考基因组缺失是常态, 证据是单向的:

  比对上且结构成立  -> 内源 (EVE) 方向的证据 (下面两种 call)
  比对不上/没基因组 -> 不定案, 只记列, 不改任何 verdict 方向

两种整合判据 (S3 采信强度不同, 见 s3_verdict.py 的 rg_adjust):

  判据一 两端宿主侧翼同一位点  contig 的两端各有一段高一致宿主比对, 两段落于**同一条
                              染色体的同一条链**, 间隔不超过 MAX_FLANK_GAP —— 这是
                              "基因组片段跨着整合位点被组装出来"的铁证 (前病毒).
                              两端都有宿主段但位点对不上 -> flank_locus_conflict
                              (更像组装嵌合, 只记列不定案).
  判据二 宿主外显子嵌合        contig 内部有宿主高一致段, 呈剪接样式 (同一染色体上
                              相邻段的 subject 间隔 >= SPLICE_GAP, 即"宿主内含子被
                              剪掉了"的转座录本) 或与任一端侧翼并存 -> 这是转座录的
                              EVE 与宿主外显子拼在一起的嵌合转录本. 注意这是 contig
                              级代理: read 级的 split-splice junction 才是金标准,
                              见 README "参考基因组通道" 一节.

另有一个移交信号: 自身参考基因组聚合覆盖 >=80% 且 >=90% 一致、又无 MP/CP/AP 病毒结构
基因 (host_dominant) —— 这就是寄主序列, 交给 S3 的寄主否决三分支 (与 S2b 同一闸门,
阈值一致). 全程遵守 v6 的教训: 带 MP/CP/AP 或上游有病毒注释的候选, 不由单一同源通道
单方面判死.

输入 refgenome_blastn.tsv: run_all.sh 生成, 每行 = 4位代码 + blastn outfmt 6 的 14 列
  (qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore
   qlen slen), query = 候选 contig, subject = 该物种参考基因组; 每个物种比对完追加
  一行 `#BLASTED	<代码>` 标记 (零命中也留标记, 区分"比过没打中"与"根本没比").
  行首 `#` 的其他行一律当注释跳过.

用法: python3 s2c_refgenome_scan.py <refgenome_blastn.tsv> <s2_domains.tsv> \
          <host_map.tsv> <out.tsv>
      python3 s2c_refgenome_scan.py --emit-header
host_map.tsv: 与 run_all.sh 同一份三列映射 (样本号 <TAB> 4位代码 <TAB> 物种目录名),
  用来把候选归属到"自身物种的参考基因组"; 样本不在映射表里 -> no_ref_genome
  (rg_scope=no_sample_map), 代码不在清单/比对失败 -> no_ref_genome (genome_absent).

输出 refgenome_evidence.tsv: 逐候选 rg_call 与结构特征, S3 以可选第 6 参读它.
"""
import argparse
import sys
from collections import defaultdict

from s2b_locus_scan import parse_sample, read_tsv   # 同一批候选, 同一份样本号口径

BLASTN_MIN_COLS = 14          # outfmt 6 的 14 列; 外加首列物种代码共 15
MIN_SEG_PID = 90.0            # 宿主段的最低一致性 (近期整合 99+; 古老整合可到 90 上下)
MIN_SEG_ALN = 100             # 单段低于这个长度不算结构证据 (防散在噪声 HSP)
END_FRAC = 0.05               # 端部窗口: query 前/后 5% 内有宿主段才算"端侧翼"
SPLICE_GAP = 500              # 相邻宿主段的 subject 间隔 >= 此值 = 剪接 (内含子) 样式
MAX_FLANK_GAP = 500_000       # 两端侧翼在同一条染色体上的最大间隔 (Cauli 基因组 ~8kb,
                              # 放宽到 500kb 容纳局部组装空缺; 再远就是两个位点)
HOSTDOM_PID, HOSTDOM_COV = 90.0, 0.80   # 与 s3_verdict.HOST_PIDENT/HOST_COV 同一闸门

RG_HDR = ("contig_id\tsample\tsample_flag\trg_code\trg_scope\trg_call\t"
          "flank_left\tflank_right\tflank_same_locus\tflank_gap\tsplice_pattern\t"
          "internal_host_segs\thost_seg_cov\tbest_host_pid\tviral_comps\t"
          "xeno_code_hit\tcall_reason\n")


def _hsp(code, p):
    """一行 blastn (带代码前缀) -> HSP dict; 坐标统一 min/max 归一化 (v4 的教训)."""
    qs, qe, ss, se = int(p[7]), int(p[8]), int(p[9]), int(p[10])
    return {"code": code, "sseqid": p[2], "pid": float(p[3]), "aln": int(p[4]),
            "qs": min(qs, qe), "qe": max(qs, qe),
            "sm": min(ss, se), "sM": max(ss, se),
            "strand": "+" if se >= ss else "-",
            "bits": float(p[12]), "qlen": int(p[13])}


def parse_refgenome_blastn(path):
    """比对表 -> (hits, blasted, stats).  hits[cid][code] = [HSP...];  blasted = 代码集."""
    hits = defaultdict(lambda: defaultdict(list))
    blasted = set()
    stats = {"n_rows": 0, "n_bad": 0}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                p = line.split("\t")
                if len(p) >= 2 and p[0] == "#BLASTED":
                    blasted.add(p[1].strip())
                continue
            stats["n_rows"] += 1
            p = line.split("\t")
            if len(p) < 1 + BLASTN_MIN_COLS:
                stats["n_bad"] += 1
                continue
            hits[p[1]][p[0]].append(_hsp(p[0], p))
    return hits, blasted, stats


def merge_segments(hsps):
    """qualifying HSP -> query 区间段. 同一染色体同链且 query 重叠/相接的 HSP 并成一段,
    段的 subject 位点取贡献碱基最多的那条 HSP (代表侧翼落点). 不同染色体/异链的
    HSP **不并** —— 那是转位/嵌合信号, 合并会把 conflict 洗成同一位点."""
    segs = []
    for h in sorted(hsps, key=lambda x: (x["qs"], x["qe"])):
        for s in segs:
            if (s["sseqid"] == h["sseqid"] and s["strand"] == h["strand"]
                    and h["qs"] <= s["qe"]):
                s["qe"] = max(s["qe"], h["qe"])
                s["ivs"].append((h["qs"], h["qe"]))
                if h["aln"] > s["rep_aln"]:
                    s.update(rep_aln=h["aln"], sm=h["sm"], sM=h["sM"], pid=h["pid"])
                break
        else:
            segs.append({"qs": h["qs"], "qe": h["qe"], "sseqid": h["sseqid"],
                         "strand": h["strand"], "sm": h["sm"], "sM": h["sM"],
                         "pid": h["pid"], "rep_aln": h["aln"],
                         "ivs": [(h["qs"], h["qe"])]})
    return segs


def classify(hsps, has_mp_cp_ap):
    """qualifying HSPs -> (call, info dict). 纯函数, 测试直接调.

    call 的判定顺序是有意的:
      host_dominant 先于一切 —— 整条 contig 基本就是寄主序列时, "两端侧翼"没有意义
      (侧翼就是它自己), 交 S3 寄主否决;
      然后判据一 (两端侧翼同一位点) 强于判据二 (外显子嵌合);
      conflict / weak 证据只记列, 不定案.
    """
    info = {"flank_left": "-", "flank_right": "-", "flank_same_locus": "FALSE",
            "flank_gap": "", "splice_pattern": "FALSE",
            "internal_host_segs": 0, "host_seg_cov": 0.0, "best_host_pid": 0.0}
    if not hsps:
        return "no_structure", info
    qlen = max(h["qlen"] for h in hsps)
    aln_tot = sum(h["aln"] for h in hsps)
    wp = sum(h["aln"] * h["pid"] for h in hsps) / max(aln_tot, 1)
    ivs = sorted((h["qs"], h["qe"]) for h in hsps)
    cov = _union(ivs) / max(qlen, 1)
    best_pid = max(h["pid"] for h in hsps)
    info.update(host_seg_cov=round(cov, 3), best_host_pid=round(best_pid, 2))

    segs = merge_segments(hsps)
    lwin = max(1, int(END_FRAC * qlen))
    rwin = qlen - lwin + 1
    lefts = [s for s in segs if s["qs"] <= lwin]
    rights = [s for s in segs if s["qe"] >= rwin]
    left = min(lefts, key=lambda s: s["qs"]) if lefts else None
    right = max(rights, key=lambda s: s["qe"]) if rights else None
    single_span = left is not None and right is not None and left is right

    if left is not None and not single_span:
        info["flank_left"] = _locus(left)
    if right is not None and not single_span:
        info["flank_right"] = _locus(right)

    both = left is not None and right is not None and not single_span
    if both:
        if (left["sseqid"] == right["sseqid"] and left["strand"] == right["strand"]):
            # gap = 两段之间的缺失碱基数 (1-bit 含端坐标差 -1), 即中间插段的量级
            gap = (right["sm"] - left["sM"] - 1) if left["strand"] == "+" \
                else (left["sm"] - right["sM"] - 1)
            info["flank_gap"] = gap
            if 0 <= gap <= MAX_FLANK_GAP:
                info["flank_same_locus"] = "TRUE"

    internal = [s for s in segs if s["qs"] > lwin and s["qe"] < rwin]
    info["internal_host_segs"] = len(internal)
    splice = False
    for a, b in zip(internal, internal[1:]):
        if a["sseqid"] == b["sseqid"] and a["strand"] == b["strand"]:
            sgap = (b["sm"] - a["sM"]) if a["strand"] == "+" else (a["sm"] - b["sM"])
            if sgap >= SPLICE_GAP:
                splice = True
    info["splice_pattern"] = "TRUE" if splice else "FALSE"

    if cov >= HOSTDOM_COV and wp >= HOSTDOM_PID and not has_mp_cp_ap:
        call = "host_dominant"
    elif both and info["flank_same_locus"] == "TRUE":
        call = "EVE_flank_confirmed"
    elif both:
        call = "flank_locus_conflict"
    elif splice or (internal and (left is not None or right is not None)):
        call = "EVE_splice_chimera"
    elif internal:
        call = "host_mix_weak"
    else:
        call = "no_structure"
    return call, info


def _union(ivs):
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


def _locus(s):
    return "%s:%d-%d(%s)" % (s["sseqid"], s["sm"], s["sM"], s["strand"])


REASON = {
    "EVE_flank_confirmed": "两端宿主侧翼落于同一条染色体同链 (判据一: 前病毒结构)",
    "flank_locus_conflict": "两端均有宿主段但位点对不上 (疑组装嵌合), 只记列不定案",
    "EVE_splice_chimera": "宿主外显子嵌合/剪接样式 (判据二, contig 级代理; read 级待验证)",
    "host_mix_weak": "内部有宿主高一致段但不满足剪接样式, 仅旁证",
    "host_dominant": "自身参考基因组聚合覆盖占优且无 MP/CP/AP -> 移交寄主否决",
    "no_structure": "有宿主命中但无侧翼/嵌合结构, 不定案",
    "no_hit": "自身参考基因组比对零命中, 不定案 (未比对上不是真病毒的证据)",
    "no_ref_genome": "无自身参考基因组 (清单未含或比对失败), 本通道跳过",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('blast_tsv', nargs='?', help='run_all.sh 生成的 refgenome_blastn.tsv')
    ap.add_argument('s2_tsv', nargs='?', help='s2_domain_scan.py 产出 (取病毒组件)')
    ap.add_argument('host_map_tsv', nargs='?',
                    help='host_map.tsv: 样本号 <TAB> 4位代码 <TAB> 物种目录名')
    ap.add_argument('out_tsv', nargs='?', help='输出 refgenome_evidence.tsv')
    ap.add_argument('--emit-header', action='store_true', help='只打印表头并退出')
    a = ap.parse_args()

    if a.emit_header:
        sys.stdout.write(RG_HDR)
        return 0
    if not all([a.blast_tsv, a.s2_tsv, a.host_map_tsv, a.out_tsv]):
        ap.error('需要 blast_tsv / s2_tsv / host_map_tsv / out_tsv (或 --emit-header)')

    sample2code = {}
    with open(a.host_map_tsv) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] and not p[0].startswith("#"):
                sample2code[p[0]] = p[1]

    hits, blasted, st = parse_refgenome_blastn(a.blast_tsv)
    if st["n_bad"] and not hits:
        sys.exit("[s2c] %s 有 %d 行但一行都解析不出来 (每行 = 物种代码 + %d 列 blastn "
                 "outfmt 6)。八成是 run_all.sh 的 -outfmt 与这里不一致, 修好再跑。"
                 % (a.blast_tsv, st["n_rows"], BLASTN_MIN_COLS))

    s2 = {r["contig_id"]: r for r in read_tsv(a.s2_tsv)}

    with open(a.out_tsv, "w") as out:
        out.write(RG_HDR)
        calls = defaultdict(int)
        n_map = 0
        for cid in sorted(set(s2) | set(hits)):
            # info 只在真正比对上时由 classify() 填; 其余分支保持空 dict,
            # 写行时统一落到默认值, 避免未初始化
            info = {}
            comps = (s2.get(cid, {}).get("components", "-") or "-")
            has_viral = any(c in comps for c in ("MP", "CP", "AP"))
            sample, sflag = parse_sample(cid)
            own = sample2code.get(sample)
            if own is None:
                call, scope, code = "no_ref_genome", "no_sample_map", "-"
            elif own not in blasted:
                call, scope, code = "no_ref_genome", "genome_absent", own
            else:
                code, scope = own, "own_genome"
                qual = [h for h in hits.get(cid, {}).get(own, [])
                        if h["pid"] >= MIN_SEG_PID and h["aln"] >= MIN_SEG_ALN]
                if not qual:
                    call = "no_hit"
                else:
                    call, info = classify(qual, has_viral)
            calls[call] += 1
            if scope == "own_genome":
                n_map += 1
            # 跨物种旁证: 其他基因组的最佳 qualifying 命中, 只展示不判定
            xeno = max(((h["code"], h["aln"], h["pid"], h["bits"])
                        for cs in hits.get(cid, {}).values() for h in cs
                        if h["code"] != code
                        and h["pid"] >= MIN_SEG_PID and h["aln"] >= MIN_SEG_ALN),
                       default=None, key=lambda t: (t[3], t[1]))
            reason = REASON[call]
            if call == "no_ref_genome" and scope == "genome_absent":
                reason += "; 代码 %s 不在可用清单" % code
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                cid, sample, sflag, code, scope, call,
                info.get("flank_left", "-"), info.get("flank_right", "-"),
                info.get("flank_same_locus", "FALSE"), info.get("flank_gap", ""),
                info.get("splice_pattern", "FALSE"), info.get("internal_host_segs", 0),
                info.get("host_seg_cov", ""), info.get("best_host_pid", ""),
                comps, ("%s/%dnt@%.1f%%" % xeno[:3]) if xeno else "-", reason))

    n_all = sum(calls.values())
    print(f"[done] s2c -> {a.out_tsv}  (映射表 {len(sample2code)} 条, "
          f"可用参考基因组 {len(blasted)} 个)")
    print(f"  候选 {n_all} 条: 自身基因组可判 {n_map} 条")
    for c in ("EVE_flank_confirmed", "flank_locus_conflict", "EVE_splice_chimera",
              "host_mix_weak", "host_dominant", "no_structure", "no_hit",
              "no_ref_genome"):
        if calls[c]:
            print(f"    {c:22s} {calls[c]}")
    n_norg = calls["no_ref_genome"]
    if n_all and n_norg == n_all:
        print("[warn] 所有候选都无可用参考基因组: 本通道整体空转 (清单为空或代码"
              "与 host_map 第二列不一致), 检查 refgenome_errors.log / genome_manifest",
              file=sys.stderr)
    elif n_norg * 2 > n_all:
        print(f"[warn] {n_norg}/{n_all} 条候选无参考基因组可判: 本通道只覆盖有现成"
              f"基因组的物种, 属预期; 这些候选不由本通道定案", file=sys.stderr)
    if st["n_bad"]:
        print(f"[warn] blastn 表有 {st['n_bad']}/{st['n_rows']} 行列数不足 "
              f"(1+{BLASTN_MIN_COLS}) 被跳过", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
