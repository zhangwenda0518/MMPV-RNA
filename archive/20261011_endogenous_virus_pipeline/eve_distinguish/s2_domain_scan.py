#!/usr/bin/env python3
"""S2 功能通道 [判定逻辑最后变更: 模块 v6]: 基因座(locus)架构 + 前病毒规模判定.
模块 v6.2 只把 provirus / multi_locus 两个判据抽成函数以便测试, 输出逐字节不变.

v3 的做法: architecture_score = 非重叠病毒基因座数, 判定里"loci>=2 或 span>=1500"
  就认为是强候选. 这其实是拿"完整性"当 EVE 证据, 方向是反的 —— RNA-seq 里真 DNA 病毒
  本来就拼不全(环状 7.5kb 基因组靠短读段只能拿到局部), "不全" 不等于 EVE.
v4 的关键倒置: 真正有判别力的是 **单一位点上 MP+CP/AP/RT/RH 齐全** (前病毒规模),
  且必须与 S1 的密码子退化**联合**才有意义:
     齐全 + 密码子退化 -> 整合的前病毒 (EVE)
     齐全 + 密码子完好 -> 感染中的真病毒 (要留的)
     不全 + 完好       -> 病毒碎片 (信息不足, 不能据此判 EVE)
  所以本模块只负责把架构事实测准, 不独自下 EVE 结论 (判定在 s3).

- MP/CP/AP 为病毒特有基因, 转座子没有 -> 架构由 locus 驱动
- TE flag: baits FALSE_* (LTR 逆转录转座子域) bitscore>40
- multi_locus_arch (v5 去同义反复): 旧定义为"组件数>=1 的基因座 >=2", 而每条 locus
  必然由带命中的 HSP 构成, 组件数恒 >=1, 所以旧定义等价于"loci>=2", 是个永真式.
  现在改成有信息量的表述: **至少两个基因座且合计组件 >=3** —— 病毒基因确实分布在
  多个基因座上, 而不是把同一段序列打成两个相邻碎片.

v6 修订 (2026-09-23): **组件计数必须有独立证据**。旧版把一条 HSP 打中的面板条目名
里列出的组件全部记上, 于是"打中一个基因组多聚蛋白条目"就等于"5 个组件齐全"。那 5 个
组件是**参考序列**的事实 (多聚蛋白确实含 5 个蛋白), 不是 contig 的事实 —— 一段比对只
覆盖 contig 上的一小段, 不可能同时证明 5 个蛋白都存在。实测 1027 候选的运行里, 有条
506 nt 的 contig 被判 5 组件齐全并进了 EVE_STRONG_provirus, 而它根本没有 5×100 aa
的独立同源区。现在每个组件要同时满足两条才计入:

  (a) query 区间独占: 同一基因座里, 与该 HSP 的 query 区间重叠的其他 HSP 已占用过的
      区间不算数。两段重叠的比对落在同一段序列上, 不可能各自证明不同的蛋白。
      界阈值按 nt 算 (3 x COMP_AA_MIN): qs/qe 是 nt 而比对长度是 aa, 混单位比会
      让一截 200 nt 的重叠尾巴也过关, 等于没设这道闸。
  (b) 比对长度封顶: 一条 N aa 的比对最多计入 N//COMP_AA_MIN 个组件 (融合条目
      MP+CP / AP+RT / RT+RH 的名字带多个组件, 但那是条目的事实)。超限时按 CANON
      序取前几个, 并在 credit_hsps 列暴露"到底几条 HSP 撑起了几个组件"。

新输出列 credit_hsps = 真正贡献了至少一个组件的 HSP 条数。它可以小于 ncomp: 一条
600 aa 的多聚蛋白比对确实能同时证 4 个组件, 那不是拼出来的。所以判"拼"不能看
ncomp > credit_hsps, 而要看 credit_hsps 本身: ncomp=4 而 credit_hsps=1 只在比对
长度也够 (>=(4-1)*COMP_AA_MIN + 起跳那 100 aa) 时成立, 可用 credit_hsps 列把这类
行挑出来人工复核比对长度。

用法: python3 s2_domain_scan.py <query.fasta> <panel_hits.tsv> <baits_hits.tsv> <out.tsv>
"""
import sys
from collections import defaultdict

CANON = ["MP", "CP", "AP", "RT", "RH"]
MIN_ALN, MIN_PID, MIN_BITS = 50, 25, 40
PROVIRUS_SPAN = 1500     # 单个基因座达到此跨度才够"前病毒"规模
PROVIRUS_NCOMP = 4       # 单基因座上至少 4/5 个规范组件
MULTI_MIN_LOCI = 2       # multi_locus_arch: 至少这么多基因座
MULTI_MIN_COMPS = 3      #               且合计这么多不同组件
COMP_AA_MIN = 100        # 一个组件要被计入, 至少要这么长的独立比对支撑 (aa)
COMP_NT_MIN = COMP_AA_MIN * 3   # 同一门槛换算成 nt: blastx 的 qs/qe 是 nt, aln 是 aa

# 输出表头: 测试直接 import S2_HDR 引用, 避免两头各抄一份
S2_HDR = ("contig_id\tarchitecture_score\tmax_locus_span\tloci_detail\tcomponents\t"
          "evidence_detail\tte_flag\tte_bitscore\tte_label\tcauli_bitscore\tcauli_label\t"
          "best_locus_span\tbest_locus_ncomp\tbest_locus_comps\tlocus_completeness\t"
          "provirus_scale\tmulti_locus_arch\tcredit_hsps\n")


def locus_ncomp(L):
    """单个基因座覆盖了几个规范组件."""
    return len([c for c in L["comps"].split("+") if c in CANON])


def fmt_comps(cs):
    """组件集合 -> CANON 序字符串.

    直接把 set 抖进 f-string 会印出 `{'MP'}` 这种 repr, 而且 set 的迭代顺序不保证,
    同一份输入两次跑出的 loci_detail 可能不同 —— s2b 靠这列做位点分组, 行序可复现就没了。
    """
    return "+".join(c for c in CANON if c in cs) or "-"


def subtract_ivs(free, lo, hi):
    """从区间列表 free 里扣掉 [lo,hi)."""
    out = []
    for a, b in free:
        if hi <= a or lo >= b:
            out.append((a, b))
            continue
        if lo > a:
            out.append((a, min(lo, b)))
        if hi < b:
            out.append((max(hi, a), b))
    return out


def credit_components(hsps, comp_aa_min=COMP_AA_MIN):
    """一个基因座里哪些组件真有独立证据 -> (组件集合, 贡献了证据的 HSP 条数).

    按 bits 降序遍历 HSP, 已判过的 HSP 的 query 区间视为已占用; 之后落进同一段序列的
    比对不再证明新组件 (两段重叠的比对不可能各证一个不同蛋白)。每条 HSP 自己还受长度
    封顶: N aa 最多证 N//comp_aa_min 个组件, 超限时按 CANON 序取前几个 —— 融合条目
    名里的多个组件是**参考**的事实, 不是这段比对的事实。

    注意单位: qs/qe 是 nt (blastx 查询坐标), aln 是 aa, 所以"可用区间够放一个蛋白"
    要按 3x 折算成 nt 再比。混着比会让 200 nt 的重叠尾巴也过关, 等于没设这道闸。
    """
    nt_min = comp_aa_min * 3
    credited, n_used, occupied = set(), 0, []
    for h in sorted(hsps, key=lambda h: (-h["bits"], h["qs"])):
        free = [(h["qs"], h["qe"])]
        for lo, hi in occupied:
            free = subtract_ivs(free, lo, hi)
        avail = sum(b - a + 1 for a, b in free)
        labels = [c for c in CANON if c in h["comps"].split("+")]
        cap = max(1, h["aln"] // comp_aa_min)
        # 可用区间要够放一个蛋白; 融合条目超限时按 CANON 序截断
        got = labels[:cap] if avail >= nt_min else []
        if got:
            n_used += 1
            credited.update(got)
        occupied.append((h["qs"], h["qe"]))
    return credited, n_used


def locus_credits(loci):
    """逐基因座算 credit_components, 返回 [(组件集合, 贡献 HSP 数), ...] 与 loci 平行."""
    return [credit_components(L["detail"]) for L in loci]


def provirus_ok(span, ncomp):
    """前病毒规模判据: 单个基因座跨度够 + 同一位点上规范组件数够.

    抽成函数而不是留在 main 里的行内表达式 —— 测试要能直接验"这个判据不是永真式";
    照着源码重抄一遍表达式只能证明抄对了, 证明不了代码里用的就是它。
    """
    return span >= PROVIRUS_SPAN and ncomp >= PROVIRUS_NCOMP


def multi_locus_ok(n_loci, comps):
    """多基因座架构判据: >=MULTI_MIN_LOCI 个基因座 **且** 合计 >=MULTI_MIN_COMPS 组件.

    旧定义 ("组件数>=1 的基因座 >=2") 是永真式 —— 每条 locus 必由带命中的 HSP 构成,
    组件数恒 >=1, 等价于"loci>=2", 没有信息量。这里保持有信息量的表述, 并把
    "它不是永真式"这件事交给测试锁住。
    """
    return n_loci >= MULTI_MIN_LOCI and len(comps) >= MULTI_MIN_COMPS


def merge_loci(hsps):
    """按 query 坐标合并重叠 HSP -> 基因座列表. 负链 HSP qstart>qend, 先归一化."""
    hsps = sorted(hsps, key=lambda h: h["qs"])
    loci = []
    for h in hsps:
        if loci and h["qs"] <= loci[-1]["qe"] + 30:  # 30nt 内视为同一基因座
            L = loci[-1]
            L["qe"] = max(L["qe"], h["qe"])
            # 按整词并集, 不能按字符遍历: 组件名是两位/三位字母 (MP/AP/RH),
            # `for c in h["comps"]` 会把 "CP" 拆成 'C' 和 'P', 于是第二个组件
            # 永远并进来, best_locus_ncomp / provirus_scale 全部低估
            for c in h["comps"].split("+"):
                if c in CANON and c not in L["comps"].split("+"):
                    L["comps"] = (L["comps"] + "+" + c) if L["comps"] else c
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

    with open(out_tsv, "w") as out:
        out.write(S2_HDR)
        n_scale = 0
        for cid in ids:
            loci = merge_loci(contig_hsps.get(cid, []))
            score = min(len(loci), 5)
            # v6: 组件只认有独立 query 区间 + 足够比对长度的 (见 credit_components)
            creds = locus_credits(loci)
            allc = set()
            for cs, _ in creds:
                allc |= cs
            comps_all = "+".join(c for c in CANON if c in allc)
            n_credit_hsps = sum(n for _, n in creds)
            loci_detail = ";".join(
                f"L{i+1}:{L['qs']}-{L['qe']}({fmt_comps(comps)},bit{L['bits']:.0f})"
                for i, (L, (comps, _)) in enumerate(zip(loci, creds)))
            detail = ";".join(
                f"{h['comps']}:{h['pident']:.0f}%/{h['aln']}aa"
                for L in loci for h in sorted(L["detail"], key=lambda x: -x["bits"])[:3])
            # v4: 以"单一基因座"为单位 — richest locus 决定前病毒规模
            if loci:
                bi = max(range(len(loci)),
                         key=lambda i: (len(creds[i][0]), loci[i]["qe"] - loci[i]["qs"]))
                best = loci[bi]
                best_span = best["qe"] - best["qs"]
                best_comps_set = creds[bi][0]
                best_n = len(best_comps_set)
                best_comps = "+".join(c for c in CANON if c in best_comps_set) or "-"
            else:
                best_span, best_n, best_comps = 0, 0, "-"
            completeness = round(best_n / len(CANON), 3)
            provirus = provirus_ok(best_span, best_n)
            if provirus:
                n_scale += 1
            # v5: 多基因座要有信息量 —— >=2 个基因座且合计 >=3 个不同组件.
            # 旧定义 "组件>=1 的基因座 >=2" 是永真式(每条 locus 必有组件), 删掉.
            multi = ("TRUE" if multi_locus_ok(len(loci),
                                              set(comps_all.split("+"))) else "FALSE")
            te = "TRUE" if te_bits.get(cid, 0) > MIN_BITS else "FALSE"
            out.write(f"{cid}\t{score}\t{max((L['qe'] - L['qs'] for L in loci), default=0)}\t"
                      f"{loci_detail or '-'}\t{comps_all or '-'}\t{detail or '-'}\t{te}\t"
                      f"{te_bits.get(cid, 0):.0f}\t{te_label.get(cid, '-')[:60]}\t"
                      f"{cauli_bits.get(cid, 0):.0f}\t{cauli_label.get(cid, '-')[:60]}\t"
                      f"{best_span}\t{best_n}\t{best_comps}\t{completeness}\t"
                      f"{'TRUE' if provirus else 'FALSE'}\t{multi}\t{n_credit_hsps}\n")
    print(f"[done] s2(logic v6, credited components): {len(ids)} contigs, "
          f"前病毒规模(单基因座>={PROVIRUS_SPAN}nt 且>={PROVIRUS_NCOMP}/5 组件) {n_scale} 条 -> {out_tsv}")


if __name__ == "__main__":
    main()
