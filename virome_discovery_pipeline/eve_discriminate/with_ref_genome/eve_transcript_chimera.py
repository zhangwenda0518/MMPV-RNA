#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_transcript_chimera.py — 转录组嵌合检测 (无宿主基因组时的整合旁证)

原理
----
整合进宿主基因组的病毒片段, 其转录本两端会是**宿主外显子**, 中间是病毒序列。
这种"宿主-病毒-宿主"嵌合结构, 在没有宿主基因组时是**最接近整合位点的证据**。

判定 (基于候选 vs 宿主转录组 blastn 的 query 坐标):
  * 候选的**前 X%** 有宿主命中 (高一致性)  -> 左侧宿主段
  * 候选的**后 X%** 有宿主命中            -> 右侧宿主段
  * 中间是病毒序列 (无宿主命中)
  * 且左右两段不重叠 (中间确实有间隔)

与 eve_host_integration_scan.py 的区别
-----------------------------------
后者要求两端落在**参考基因组同一染色体同一位点** (真正的整合位点);
本脚本只能要求"两端都是宿主", 证据强度低一档 (转录本嵌合也可能是
反式剪接/组装错误), 但远强于没有任何证据。

输出: <out>/chimera.tsv
  contig_id  chimera_call  left_cov  right_cov  mid_gap  host_identity
             n_host_segs  reason
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

HDR = ("contig_id\tchimera_call\tleft_cov\tright_cov\tmid_gap\t"
       "host_identity\tn_host_segs\treason\n")

END_FRAC = 0.15        # 端部窗口: 前/后 15% 内有宿主命中才算"端侧翼"
MIN_HOST_PID = 85.0    # 宿主段的比对一致性下限
MIN_SEG_ALN = 100      # 单段最短长度 (nt)
# v7.5: 中间间隔下限 —— 实测发现 227 条"嵌合"里 92 条的间隔 <100nt,
# 那只是同一段宿主序列被 blastn 打断成多个 HSP 的空档, 不是"中间夹着病毒"。
# 真正的整合位点中间应该有病毒序列 (通常 >=500nt), 所以要求间隔至少这么多。
MIN_MID_GAP = 500      # 中间间隔下限 (nt)


def parse_blastn(path):
    """候选 vs 宿主转录组 blastn -> {qid: [hsp, ...]}。

    每 hsp: {qs, qe, pid, aln, qlen, sid}
    """
    hits = defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 14:
                continue
            try:
                pid = float(c[2])
                aln = int(c[3])
                qs, qe = int(c[6]), int(c[7])
                qlen = int(c[12])
            except ValueError:
                continue
            if pid < MIN_HOST_PID or aln < MIN_SEG_ALN:
                continue
            hits[c[0]].append({
                "qs": min(qs, qe), "qe": max(qs, qe),
                "pid": pid, "aln": aln, "qlen": qlen, "sid": c[1],
            })
    return hits


def classify_contig(hsps, end_frac=END_FRAC, min_mid_gap=MIN_MID_GAP):
    """单条候选的嵌合判定 -> (call, info dict)。

    call:
      chimera_supported  两端都有宿主段且中间有间隔
      host_flanked       只有一端有宿主段 (弱)
      no_host_structure  无端部宿主段
    """
    info = {"left_cov": 0.0, "right_cov": 0.0, "mid_gap": 0,
            "host_identity": 0.0, "n_host_segs": 0}
    if not hsps:
        return "no_host_structure", info
    qlen = max(h["qlen"] for h in hsps)
    if qlen < 200:
        return "no_host_structure", info
    info["n_host_segs"] = len(hsps)
    info["host_identity"] = round(max(h["pid"] for h in hsps), 2)
    lwin = max(1, int(end_frac * qlen))
    rwin = qlen - lwin + 1
    lefts = [h for h in hsps if h["qs"] <= lwin]
    rights = [h for h in hsps if h["qe"] >= rwin]
    # 左右段覆盖: 各自并集长度 / 窗口长度 (封顶 1.0)
    def cov(segs, lo, hi):
        if not segs:
            return 0.0
        ivs = sorted((max(s["qs"], lo), min(s["qe"], hi)) for s in segs)
        tot = 0
        cs, ce = ivs[0]
        for a, b in ivs[1:]:
            if a <= ce + 1:
                ce = max(ce, b)
            else:
                tot += ce - cs + 1
                cs, ce = a, b
        tot += ce - cs + 1
        return round(min(tot / max(hi - lo + 1, 1), 1.0), 3)
    info["left_cov"] = cov(lefts, 1, lwin)
    info["right_cov"] = cov(rights, rwin, qlen)
    # 中间间隔: 左段最右端 到 右段最左端
    if lefts and rights:
        lmax = max(h["qe"] for h in lefts)
        rmin = min(h["qs"] for h in rights)
        info["mid_gap"] = max(0, rmin - lmax)
    if (info["left_cov"] >= 0.5 and info["right_cov"] >= 0.5
            and info["mid_gap"] >= MIN_MID_GAP):
        return "chimera_supported", info
    if info["left_cov"] >= 0.5 or info["right_cov"] >= 0.5:
        return "host_flanked", info
    return "no_host_structure", info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blastn", required=True, help="候选 vs 宿主转录组 blastn")
    ap.add_argument("--query-fa", required=True, help="候选 fasta (取全部 ID)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--end-frac", type=float, default=END_FRAC)
    ap.add_argument("--min-mid-gap", type=int, default=MIN_MID_GAP,
                    help="中间间隔下限 (nt); 小于此值视为比对空档而非嵌合")
    a = ap.parse_args(argv)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ids = []
    with open(a.query_fa, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].strip()
                if h:
                    ids.append(h.split()[0])
    hits = parse_blastn(a.blastn)
    print(f"[load] 候选 {len(ids)} / 有宿主命中 {len(hits)}", file=sys.stderr)

    from collections import Counter
    cnt = Counter()
    rows = []
    for cid in ids:
        call, info = classify_contig(hits.get(cid, []), a.end_frac,
                                     a.min_mid_gap)
        cnt[call] += 1
        if call == "chimera_supported":
            why = (f"两端均有宿主段 (左 {info['left_cov']:.2f} / 右 "
                   f"{info['right_cov']:.2f}), 中间间隔 {info['mid_gap']}nt")
        elif call == "host_flanked":
            why = (f"仅一端有宿主段 (左 {info['left_cov']:.2f} / 右 "
                   f"{info['right_cov']:.2f})")
        else:
            why = "无端部宿主段"
        rows.append((cid, call, info["left_cov"], info["right_cov"],
                     info["mid_gap"], info["host_identity"],
                     info["n_host_segs"], why))

    tsv = out / "chimera.tsv"
    with open(tsv, "w", encoding="utf-8", newline="\n") as fo:
        fo.write(HDR)
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")

    print(f"\n[done] {tsv}  ({len(rows)} 条)")
    print("\n判定分布:")
    for k, v in cnt.most_common():
        print(f"  {k:22s} {v:5d}  ({v/len(rows)*100:.1f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
