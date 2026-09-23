#!/usr/bin/env python3
"""S3 判别 v4: 单一位点架构 (S2/S2b) x 密码子退化 (S1) -> verdict.

v3 的规则把"完整"当 EVE 证据: ORF 完整 + loci>=2 或 span>=1500 -> review_extend;
stops>=3 -> ancient_EVE. 两个方向都不对:
  * RNA-seq 里真 DNA 病毒本来就拼不全, "不全" 不能推出 EVE;
  * stops>=3 是绝对计数, 被组装 indel 污染 (v4 S1 已改用同批读框基线 + 分布剖面).
v4 判定表 (只对 Caulimoviridae / 面板可用的 DNA 候选):

  密码子退化 \\ 单一位点架构   齐全(>=4/5 组件, 或 S2b 位点齐全)   部分/单组件
  -------------------------  ----------------------------------  ---------------
  分布式退化 (化石签名)       EVE_STRONG_provirus (整合前病毒)     ancient_EVE / EVE_suspect
  完好 (编码受保留)           virus_candidate (感染中的真病毒)     virus_fragment_review
  组装断裂 (单点 indel)       assembly_breakpoint_review          assembly_breakpoint_review
  组成本底噪声                compositional_noise_review         compositional_noise_review

  注意"部分 + 退化"只给 EVE_suspect: 它既可能是化石碎片, 也可能是有 indel 的真病毒碎片,
  contig 级无法定案, 必须靠 S2b 位点或寄主基因组.
  非 Cauli (面板不适用): 只走结构通道, 退而不判 EVE.

输入: s1_decay.tsv s2_domains.tsv [s2b_locus_architecture.tsv] rescue_evidence_scored.tsv <out.tsv>
"""
import sys
import os
from collections import Counter

HOST_PIDENT, HOST_COV = 90.0, 0.80   # 与自身寄主物种 >=90% 一致且区间覆盖 >=80% -> 寄主序列
POL_COMPS = {"RT", "RH"}   # 纯 pol 区 (无 MP/CP/AP) 更像转座子/退化 pol


def read_tsv(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        d = {}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p:
                d[p[0]] = dict(zip(hdr, p))
        return d


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if len(args) < 5:
        print(__doc__)
        sys.exit(1)
    s1_tsv, s2_tsv, s2b_tsv, ev_tsv, out_tsv = args[0], args[1], args[2], args[3], args[4]
    s2b = read_tsv(s2b_tsv)
    s1 = read_tsv(s1_tsv)
    s2 = read_tsv(s2_tsv)

    ev = read_tsv(ev_tsv)

    out = open(out_tsv, "w")
    out.write("contig_id\ttax_family\tcategory\tcheckv_completeness\tlocus_ncomp\tlocus_completeness\t"
              "provirus_scale\tlocus_arch\tlocus_members\tdecay_class\tstop_enrichment\t"
              "premature_stops_region\tframe_switches\torf_max_fraction\tte_flag\thost_wpid\t"
              "host_cov\tverdict\tverdict_reason\n")
    verdicts = Counter()
    for cid, d1 in s1.items():
        d2 = s2.get(cid, {})
        db = s2b.get(cid, {})
        e = ev.get(cid, {})
        fam = e.get("tax_family", "")
        checkv = _f(e.get("checkv_completeness"))
        best_n = int(d2.get("best_locus_ncomp") or 0)
        comps_str = d2.get("components", "-") or "-"
        provirus = d2.get("provirus_scale", "FALSE") == "TRUE"
        locus_arch = db.get("locus_arch", "")
        locus_hint = db.get("locus_hint", "")
        locus_full = (locus_arch == "locus_full") or provirus
        decay = d1.get("decay_class", "") or "no_hsp"
        enrich = _f(d1.get("stop_enrichment"))
        stops = _f(d1.get("premature_stops_region"))
        stops = int(stops) if stops is not None else None
        switches = _f(d1.get("frame_switches"))
        switches = int(switches) if switches is not None else None
        te = d2.get("te_flag", "FALSE")
        te_bits = _f(d2.get("te_bitscore")) or 0
        cauli_bits = _f(d2.get("cauli_bitscore")) or 0
        has_mp_cp_ap = any(c in comps_str for c in ("MP", "CP", "AP"))
        pol_only = bool(comps_str not in ("-", "")) and not has_mp_cp_ap
        hpid, hcov = _f(db.get("host_wpid")) or 0.0, _f(db.get("host_cov")) or 0.0
        reason = []

        def emit(v):
            verdicts[v] += 1
            out.write(f"{cid}\t{fam}\t{e.get('category','')}\t{e.get('checkv_completeness','')}\t"
                      f"{best_n}\t{d2.get('locus_completeness','')}\t{provirus}\t{locus_arch}\t"
                      f"{db.get('locus_members','')}\t{decay}\t{enrich}\t{stops}\t{switches}\t"
                      f"{d1.get('orf_max_fraction','')}\t{te}\t{hpid}\t{hcov}\t"
                      f"{v}\t{'|'.join(reason)}\n")

        # ---- 寄主同源否决: 候选就是自身寄主的序列 (如 HSP70 假阳性) ----
        # OneKP 转录组组装是碎的, 所以用物种级区间并集覆盖而非单条比对覆盖
        if hpid >= HOST_PIDENT and hcov >= HOST_COV and not has_mp_cp_ap:
            reason.append(f"自身寄主物种聚合 {hpid:.1f}% 一致性 / {hcov:.2f} 区间覆盖, "
                          "且无 MP/CP/AP 病毒结构基因 -> 寄主序列 (含寄主基因与整合拷贝), 非感染病毒")
            emit("host_contamination_likely")

        # ---- 非 Cauli / 无面板命中: 只走结构通道, 退而不判 EVE ----
        elif fam != "Caulimoviridae" and comps_str in ("-", ""):
            if te == "TRUE" and te_bits >= 100:
                reason.append("TE bitscore 领先且无病毒结构基因")
                emit("EVE_LTR_TE")
            elif decay == "coding_intact":
                reason.append("结构完整/编码保留, 面板不适用")
                emit("structure_intact_review")
            elif decay == "distributed_decay":
                reason.append("分布式退化但无 Cauli 面板证据, 不下 EVE 结论")
                emit("structure_decay_review")
            else:
                reason.append("证据不足")
                emit("review")

        # ---- TE 否决 (仅适用于无 MP/CP/AP 的纯 pol 区) ----
        elif te == "TRUE" and te_bits > cauli_bits and pol_only:
            reason.append("纯 pol 区且 TE bitscore 领先 -> 转座子")
            emit("EVE_LTR_TE")
        elif fam == "Metaviridae":
            reason.append("Metaviridae 即 LTR 逆转录转座子科")
            emit("EVE_LTR_TE")

        else:
            # ---- 二维判定: 单一位点架构 x 密码子退化 ----
            if decay == "distributed_decay":
                if locus_full:
                    reason.append("单一位点架构齐全(>=4/5 组件)且密码子分布式退化 = 整合前病毒")
                    emit("EVE_STRONG_provirus")
                elif pol_only:
                    reason.append("纯 pol 区分布式退化")
                    emit("ancient_EVE")
                else:
                    reason.append("架构不全 + 分布式退化; contig 级无法定案, 需位点/寄主基因组")
                    emit("EVE_suspect")
            elif decay == "coding_intact":
                if locus_full:
                    cv = "" if checkv is None else f", CheckV 完整度 {checkv:.0f}%"
                    reason.append("单一位点架构齐全且密码子完好 = 感染中的真病毒 (非 EVE)" + cv)
                    if locus_hint == "locus_full_but_assembled_breakpoint":
                        reason.append("位点齐全但存在组装断裂, 复拼后可升级")
                    emit("virus_candidate")
                else:
                    reason.append("架构不全 + 编码完好: RNA-seq 里真 DNA 病毒本就拼不全, 不判 EVE")
                    emit("virus_fragment_review")
            elif decay == "assembly_breakpoint":
                reason.append("终止集中在单个断裂点 (同批读框基线已排除组成本底), "
                              "更像组装 indel, 建议 rescue 延伸后复判")
                emit("assembly_breakpoint_review")
            elif decay == "compositional_noise":
                reason.append("病毒读框不比其他读框更脏, 终止来自组成/组装噪声, 无退化证据")
                emit("compositional_noise_review")
            else:
                reason.append("同源区过短或模式不明")
                emit("review")
    out.close()
    print(f"[done] s3(v4) -> {out_tsv}")
    print("verdict 分布:")
    for v, n in verdicts.most_common():
        print(f"  {v:30s} {n}")


if __name__ == "__main__":
    main()
