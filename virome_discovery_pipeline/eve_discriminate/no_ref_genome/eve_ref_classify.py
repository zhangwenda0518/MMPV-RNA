#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_ref_classify.py — 候选 vs EVE 元件库 / 真实病毒库 的双侧判定

背景
----
机器学习实验证明: 序列内在特征**无法**区分内源与真实病毒
(硬负样本下 AUC 0.66)。可行的替代是**同源参照**:

  * 候选命中 **已知植物 EVE 元件库** (41,933 条, 有基因组坐标背书)
    -> 与"已经整合进宿主基因组的东西"同源 = 内源旁证
  * 候选命中 **真实病毒基因组库** -> 外源病毒证据
  * 两侧都命中 -> 该家族既有内源也有外源 (如 Caulimoviridae), 无法二分

判定规则 (双侧 bitscore 比较, 与 A 线 Stage2 同一哲学):
  EVE 侧最佳 bitscore >= 阈值 且 > 病毒侧  -> 内源旁证强
  病毒侧最佳 >= 阈值 且 > EVE 侧          -> 更像外源病毒
  两侧都 >= 阈值且接近                     -> 两者兼有 (需人工)
  都不过阈值                              -> 无法判定

输出: <out>/ref_classify.tsv
  contig_id  ref_class  eve_bits  eve_id  eve_pident  virus_bits  virus_id
  virus_pident  reason
"""
import argparse
import sys
from pathlib import Path

CLASS_EVE = "endogenous_supported"
CLASS_VIRUS = "exogenous_likely"
CLASS_BOTH = "both_families"
CLASS_NONE = "undetermined"

HDR = ("contig_id\tref_class\teve_bits\teve_id\teve_pident\t"
       "virus_bits\tvirus_id\tvirus_pident\treason\n")


def load_best(path):
    """blastn outfmt 6 -> {qid: (bits, sid, pident, aln, qlen)} 每 query 取最佳。"""
    best = {}
    if not path or not Path(path).is_file():
        return best
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            q, sid = c[0], c[1]
            try:
                bits = float(c[11])
                pid = float(c[2])
                aln = int(c[3])
                qlen = int(c[4])
            except ValueError:
                continue
            if q not in best or bits > best[q][0]:
                best[q] = (bits, sid, pid, aln, qlen)
    return best


def read_ids(fasta):
    ids = []
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].strip()
                if h:
                    ids.append(h.split()[0])
    return ids


def classify(eve, vir, min_bits, ratio):
    """单条候选判定 -> (class, reason)。"""
    e = eve.get("bits", 0.0)
    v = vir.get("bits", 0.0)
    if e >= min_bits and e > v:
        if v >= min_bits and v >= e * ratio:
            return CLASS_BOTH, (f"两侧均过阈值且接近 (EVE {e:.0f} vs 病毒 {v:.0f}): "
                                f"该家族既有内源也有外源, 无法二分")
        return CLASS_EVE, (f"命中已知植物 EVE 元件 (bits {e:.0f}), "
                           f"高于病毒侧 (bits {v:.0f})")
    if v >= min_bits and v > e:
        if e >= min_bits and e >= v * ratio:
            return CLASS_BOTH, (f"两侧均过阈值且接近 (EVE {e:.0f} vs 病毒 {v:.0f})")
        return CLASS_VIRUS, (f"命中真实病毒基因组 (bits {v:.0f}), "
                             f"高于 EVE 侧 (bits {e:.0f})")
    if e >= min_bits and v >= min_bits:
        return CLASS_BOTH, (f"两侧均过阈值 (EVE {e:.0f} vs 病毒 {v:.0f})")
    if e >= min_bits:
        return CLASS_EVE, f"仅命中 EVE 元件 (bits {e:.0f})"
    if v >= min_bits:
        return CLASS_VIRUS, f"仅命中真实病毒 (bits {v:.0f})"
    return CLASS_NONE, (f"两侧均未过阈值 (EVE {e:.0f} / 病毒 {v:.0f}, "
                        f"门槛 {min_bits:.0f})")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query-fa", required=True)
    ap.add_argument("--eve-hits", required=True, help="候选 vs EVE 元件 blastn")
    ap.add_argument("--virus-hits", required=True, help="候选 vs 真实病毒 blastn")
    ap.add_argument("--eve-hits-prot", default="",
                    help="可选: 候选 vs EVE 元件蛋白库 diamond blastx (v7.6 Δ 占优子标签)")
    ap.add_argument("--virus-hits-prot", default="",
                    help="可选: 候选 vs 真实病毒蛋白库 diamond blastx")
    ap.add_argument("--delta-dominant", type=float, default=10.0,
                    help="both 类的 Δ 占优阈值: Δ=EXO−EVE (蛋白级), |Δ|>=该值给占优子标签")
    ap.add_argument("--classify-on", choices=["nt", "prot"], default="nt",
                    help="分类依据: nt=核苷级 bits (原始口径); prot=蛋白级 bits (需给 "
                         "--*-hits-prot, 敏感度更高: 1027 候选 nt undetermined 788 -> "
                         "prot 170, both 5 -> 287)。给 prot 时 eve_bits/virus_bits 列"
                         "仍保持核苷级口径, 另加 classify_basis 列标注")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-bits", type=float, default=50.0,
                    help="两侧判定阈值 (与 A 线 Stage2 的 HOST_BS/VIRAL_BS 一致)")
    ap.add_argument("--ratio", type=float, default=0.8,
                    help="两侧接近的判据: 较弱侧 >= 较强侧 * 该值 -> both_families")
    a = ap.parse_args(argv)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ids = read_ids(a.query_fa)
    eve = load_best(a.eve_hits)
    vir = load_best(a.virus_hits)
    # v7.6: 蛋白级 bits (可选) -> Δ=EXO−EVE 的占优子标签。
    # 1027 候选实测 (2026-10-20): both_families 287 条按 Δ 干净三分 106/104/77,
    # 且 endogenous_supported 254 条 100% Δ<=-10、exogenous_likely 316 条 100% Δ>=+10,
    # 两套规则方向完全自洽 —— 子标签把 both 类从"人工定案"细化为"带倾向的人工定案"。
    eve_p = load_best(a.eve_hits_prot) if a.eve_hits_prot else {}
    vir_p = load_best(a.virus_hits_prot) if a.virus_hits_prot else {}
    print(f"[load] 候选 {len(ids)} / EVE 命中 {len(eve)} / 病毒命中 {len(vir)}"
          + (f" / 蛋白级 EVE {len(eve_p)} / 病毒 {len(vir_p)}" if eve_p or vir_p else ""),
          file=sys.stderr)

    rows = []
    from collections import Counter
    cnt = Counter()
    dom_cnt = Counter()
    use_prot = bool(eve_p or vir_p)
    classify_prot = use_prot and a.classify_on == "prot"
    hdr = HDR.rstrip("\n") + ("\tdelta_prot\tdominance\tclassify_basis\n"
                              if classify_prot else
                              ("\tdelta_prot\tdominance\n" if use_prot else "\n"))
    for cid in ids:
        e = eve.get(cid)
        v = vir.get(cid)
        er = {"bits": e[0], "sid": e[1], "pid": e[2]} if e else {}
        vr = {"bits": v[0], "sid": v[1], "pid": v[2]} if v else {}
        # --classify-on prot: 分类用蛋白级 bits (敏感度高), 表内 eve_bits/virus_bits
        # 列保持核苷级口径不变, basis 列写明分类依据
        if classify_prot:
            ep = eve_p.get(cid)
            vp = vir_p.get(cid)
            cls, why = classify({"bits": ep[0]} if ep else {},
                                {"bits": vp[0]} if vp else {},
                                a.min_bits, a.ratio)
        else:
            cls, why = classify(er, vr, a.min_bits, a.ratio)
        cnt[cls] += 1
        base = (cid, cls,
                round(er.get("bits", 0), 1), er.get("sid", "-"),
                round(er.get("pid", 0), 2),
                round(vr.get("bits", 0), 1), vr.get("sid", "-"),
                round(vr.get("pid", 0), 2), why)
        if use_prot:
            d = vir_p.get(cid, (0.0,))[0] - eve_p.get(cid, (0.0,))[0]
            if cls == CLASS_BOTH:
                dom = ("both_eve_dominant" if d <= -a.delta_dominant else
                       "both_exo_dominant" if d >= a.delta_dominant else "both_balanced")
            elif cls == CLASS_EVE:
                dom = "eve"
            elif cls == CLASS_VIRUS:
                dom = "exo"
            else:
                dom = "none"
            dom_cnt[dom] += 1
            rows.append(base + (round(d, 1), dom,
                                "prot" if classify_prot else "nt"))
        else:
            rows.append(base)

    tsv = out / "ref_classify.tsv"
    with open(tsv, "w", encoding="utf-8", newline="\n") as fo:
        fo.write(hdr)
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")

    print(f"\n[done] {tsv}  ({len(rows)} 条)")
    print("\n判定分布:")
    for k, v in cnt.most_common():
        print(f"  {k:24s} {v:5d}  ({v/len(rows)*100:.1f}%)")
    if use_prot:
        print("\nΔ 占优子标签分布 (蛋白级, Δ=EXO−EVE):")
        for k, v in dom_cnt.most_common():
            print(f"  {k:22s} {v:5d}")

    # 附加统计: 命中 EVE 的候选里, 有多少同时命中病毒 (家族重叠)
    both_hit = sum(1 for cid in ids if cid in eve and cid in vir)
    print(f"\n两侧都有命中: {both_hit} 条")
    print(f"仅 EVE 命中: {sum(1 for cid in ids if cid in eve and cid not in vir)} 条")
    print(f"仅病毒命中: {sum(1 for cid in ids if cid in vir and cid not in eve)} 条")
    print(f"两侧都无: {sum(1 for cid in ids if cid not in eve and cid not in vir)} 条")
    return 0


if __name__ == "__main__":
    sys.exit(main())
