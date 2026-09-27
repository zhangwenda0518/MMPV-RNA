#!/usr/bin/env python3
"""S1 结构通道 [判定逻辑最后变更: 模块 v6.3]: 编码退化签名 + 同链读框基线.

v3 的缺陷 (2026-09-23 修订): premature_stops_region 用绝对计数 + 固定阈值 (stops>=3 判退化),
  会被组装 indel 污染 —— 一个单碱基 indel 就足以让其后全部密码子错读, 在远缘同源区
  (25-45% pid) 上造出十几个"提前终止", 与真正的化石退化无法区分.
v4 的做法: 不跟绝对阈值较劲, 而是把每个 contig **自身**当基线:

  * 同源区同一段序列、同一条链、另外两个读框的终止密度 = 本地 null.
    这个 null 是**读框基线(组成本底)**, 不是"同批组装错误率": 三个读框读的是同一段
    碱基, 终止密码子的期望密度由该段的碱基/密码子组成决定 (随机序列约 3/64 ≈ 4.7%),
    换一次组装、换一个样本, 它都还在同一个量级。它能扣掉的只是"这段序列本来就容易读出
    终止"这一项, 扣不掉组装 indel —— indel 只污染比对上的那一个读框, 另外两框根本不在
    读框上, 所以对"病毒读框被 indel 打断"这件事, 基线提供不了任何保护。v6 起把旧文档里
    "同批组装错误率基线"的说法改掉: 那个叫法让人以为 indel 已被扣除, 并没有。
  * enrichment = 病毒读框终止数 / (另外两框均数 + 1):
      - 远大于 1   病毒读框比噪声底**更脏** = 退化证据 (终止密码子远超基线)
      - 接近 1     病毒读框和随机读框一样脏 = 无退化证据, 终止来自组成/组装噪声
      - 明显小于 1 病毒读框比另外两框**更干净** = 编码仍受选择保留 (同源区确实在翻译)
    注意方向: 分子是"终止数", 越大越脏。旧版 docstring 把这句写反了 (说远大于 1
    "更干净"), 与代码行为相反, 照它读会把化石签名和完好编码搞混。
  * 停止分布剖面区分两种"接近 1":
      - max_stopfree_frac >= 0.25  -> 存在一段完整可读长区, 是**单个 indel 断裂**
      - head_stops >= 2           -> 终止从同源区最前端起就密集出现, 是**分布式退化**
    前者是组装噪声, 后者才是化石签名。head_stops 取 2 而不是 1: 单个终止落在最前
    15 个密码子里, 同样可能只是那里刚好有个 indel, 不足以当化石签名 (旧版因此把
    只有 2 个终止的 contig 判成"分布式退化")。

v5 另修一处归类错误: stops<=1 才认"编码完好"太窄。实测有 9/40 条 compositional_noise
的病毒读框比同链另两框干净一个数量级 (stops<=5 而另两框均数几十) —— 那是编码受保留
的强签名, 不是噪声。现在 enrich<=0.5 且 stops<=5 的判回 coding_intact。

v6.3 修一处系统性高估: stop_profile 的最长无终止段旧版按相邻终止的"起点差"计,
把上游终止密码子自身占的 3nt 也当成了无终止, 每个间隔高估 3nt —— msf 在 0.25 闸门
附近的条目可能被错放进 assembly_breakpoint。改为间隔扣 3nt、末段从最后一个终止的
结尾起算; 语义不变, 数值系统性变小 (最多差 3nt/span)。

输入: 候选 fasta + DIAMOND blastx 结果 (含 qframe, panel 与 baits 各跑一次合并传入).
输出: s1_decay.tsv (逐 contig).

用法: python3 s1_decay_scan.py <query.fasta> <diamond_hits.tsv> <out.tsv>
diamond_hits.tsv 列: qseqid sseqid pident length qstart qend sstart send evalue bitscore qframe subject_label
"""
import sys
from collections import defaultdict

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
ORF_MIN_NT = 300  # 100 aa

HEAD_CODONS = 15   # 同源区最前端 15 个密码子, 用于判断"终止是否从一开头就有"
MIN_SPAN = 90      # 同源区短于此不做基线判断 (噪声太大)

# v5 判定常数 (见 classify_decay):
CLEAN_ENRICH = 0.5     # 病毒读框比同链另两框干净到这个程度, 认为编码仍受选择保留
CLEAN_MAX_STOPS = 5    # 且绝对终止数也不多 (再多就说不清了, 归 compositional_noise)
MIN_HEAD_STOPS = 2     # 头部至少这么多个终止才认"分布式退化" (1 个更像恰好一个 indel)
MIN_DECAY_STOPS = 3    # 或终止数密集到这个程度, 才认"分布式退化"

# 输出表头: 测试直接 import S1_HDR 引用, 避免两头各抄一份
S1_HDR = ("contig_id\tlength\tdom_frame\tmain_orf_aa\torf_max_fraction\tcoding_density\t"
          "stops_dom_per_kb\tframe_switches\tn_hsp_frames\tpremature_stops_region\t"
          "hit_bitscore\thit_label\thsp_span\tspan_len\tstops_other_mean\tstop_enrichment\t"
          "head_stops\tfirst_stop_frac\tmax_stopfree_frac\tdecay_class\n")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def frame_intervals(seq, offset):
    """单读框 (offset=0,1,2) 的 stop-to-stop 区间 [(start_nt,end_nt_excl)], 序列全长坐标."""
    n = []
    i = offset
    while i + 3 <= len(seq):
        if seq[i:i + 3].upper() in STOPS:
            n.append(i)
        i += 3
    bounds = [-3] + n + [i]
    out = []
    for a, b in zip(bounds, bounds[1:]):
        s, e = a + 3, b
        if e - s >= ORF_MIN_NT:
            out.append((s, e))
    return out, len(n)


def orf_scan(seq):
    """6 框 ORF 扫描, 返回主读框统计."""
    L = len(seq)
    best = {"frame": None, "orf_nt": 0, "coding": 0}
    stats = {}
    for strand, sseq in (("+", seq), ("-", revcomp(seq))):
        for offset in range(3):
            ivs, n_stops = frame_intervals(sseq, offset)
            frame = f"{strand}{offset+1}"
            coding_nt = sum(e - s for s, e in ivs)
            stats[frame] = (coding_nt, n_stops, ivs)
            if coding_nt > best["coding"] or (
                    coding_nt == best["coding"] and ivs and
                    max(e - s for s, e in ivs) > best["orf_nt"]):
                best = {"frame": frame, "orf_nt": max(e - s for s, e in ivs) if ivs else 0,
                        "coding": coding_nt}
    dom = best["frame"]
    if dom is None:
        dom = "+1"
        best["orf_nt"] = 0
    coding_nt, n_stops, ivs = stats[dom]
    main_orf_nt = best["orf_nt"]
    return {
        "dom_frame": dom,
        "main_orf_aa": main_orf_nt // 3,
        "orf_max_fraction": round(main_orf_nt / L, 4) if L else 0,
        "coding_density": round(coding_nt / L, 4) if L else 0,
        "stops_dom_per_kb": round(n_stops * 1000 / L, 3) if L else 0,
    }


def scan_frame_stops(sseq, start, end, offset):
    """[start,end) 内指定读框的终止密码子位置列表 (坐标相对 start)."""
    pos = []
    i = start + offset
    while i + 3 <= end:
        if sseq[i:i + 3].upper() in STOPS:
            pos.append(i - start)
        i += 3
    return pos


def stop_profile(pos, span_nt):
    """终止密码子分布剖面 -> 判断 分布式退化 vs 单点断裂."""
    if not pos or span_nt <= 0:
        return {"head_stops": 0, "first_stop_frac": "", "max_stopfree_frac": ""}
    head_end = min(HEAD_CODONS * 3, span_nt)
    head_stops = sum(1 for p in pos if p < head_end)
    first_stop_frac = round(pos[0] / span_nt, 4)
    # 最长无终止区间 (含两端截断)。相邻终止之间要扣掉上游终止密码子自身占的 3nt
    # (pos 是终止的起点): v6.2 及之前的 b-a 把那 3nt 也算成"无终止", 每个间隔
    # 高估 3nt, msf 逼近 0.25 闸门的条目会被错放进 assembly_breakpoint。
    # 末段从最后一个终止的结尾 (pos[-1]+3) 算到同源区末尾; 扫描循环保证
    # pos[-1]+3 <= span_nt, 所以每段都 >= 0, 不需要额外守卫。
    runs = ([pos[0]]
            + [b - a - 3 for a, b in zip(pos, pos[1:])]
            + [span_nt - pos[-1] - 3])
    max_stopfree_frac = round(max(runs) / span_nt, 4)
    return {"head_stops": head_stops, "first_stop_frac": first_stop_frac,
            "max_stopfree_frac": max_stopfree_frac}


def classify_decay(stops, other_mean, prof, span_nt):
    """v5 分类: 编码保留 / 分布式退化(化石签名) / 组装断裂 / 组成本底噪声.

    enrichment = (stops+1)/(other_mean+1), 方向是"越大越脏":
      enrich <= CLEAN_ENRICH  -> 病毒读框明显比基线干净, 是编码受保留, 不是退化
      enrich < 1              -> 与基线无差, 终止来自序列组成/组装噪声
      enrich >= 1             -> 病毒读框确实更脏, 再看分布剖面定性
    """
    if stops is None or span_nt < MIN_SPAN:
        return "no_hsp"
    if stops <= 1:
        return "coding_intact"
    enrich = (stops + 1) / (other_mean + 1)
    if enrich < 1.0:
        # 病毒读框并不比另外两框更脏
        if stops <= CLEAN_MAX_STOPS and enrich <= CLEAN_ENRICH:
            # 比基线干净一个数量级 -> 编码仍在被选择保留 (同源区就是 ORF),
            # 旧版把这些一律归成 compositional_noise, 白丢了"编码完好"这一路证据
            return "coding_intact"
        return "compositional_noise"
    msf = prof["max_stopfree_frac"]
    head = prof["head_stops"]
    if msf != "" and msf >= 0.25:
        # 存在一段 >=1/4 同源区长度的完整无终止区 -> 单个 indel 断裂, 而非全段退化.
        # (断裂点两侧各是一段完整可读区, 所以最长无终止段必然占去大半)
        return "assembly_breakpoint"
    if head >= MIN_HEAD_STOPS or stops >= MIN_DECAY_STOPS:
        # 终止从最前端起就密集出现, 或数量密集且无任何长完整段 -> 全段分布式退化 = 化石签名.
        # 注意不能用 "head==0" 反推断裂: 实测 92 终止的铁证 contig 前端恰好无终止,
        # 但 msf 仅 0.059, 显然不是单点断裂.
        return "distributed_decay"
    return "assembly_breakpoint"


def hsp_metrics(hsps, seq):
    """读框漂移 + 同源区内 各读框终止分布 + 基线化退化分类."""
    empty = {"frame_switches": "", "n_hsp_frames": "", "premature_stops_region": "",
             "hit_bitscore": "", "hit_label": "", "hsp_span": "", "span_len": "",
             "stops_other_mean": "", "stop_enrichment": "", "head_stops": "",
             "first_stop_frac": "", "max_stopfree_frac": "", "decay_class": "no_hsp"}
    if not hsps:
        return empty
    # 按 subject bitscore 取最佳 target 的全部 HSP, 按 qstart 排序
    best_sid = max(hsps, key=lambda h: h["bits"])["sid"]
    sub = sorted([h for h in hsps if h["sid"] == best_sid], key=lambda h: h["qstart"])
    frames = [h["frame"] for h in sub]
    switches = sum(1 for a, b in zip(frames, frames[1:]) if a != b)
    qs = min(min(h["qstart"], h["qend"]) for h in sub)
    qe = max(max(h["qstart"], h["qend"]) for h in sub)
    from collections import Counter
    dom_frame = Counter(frames).most_common(1)[0][0]
    fnum = int(dom_frame)
    strand = "+" if fnum > 0 else "-"
    offset = abs(fnum) - 1
    L = len(seq)
    if strand == "+":
        sseq, scan_s, scan_e = seq, qs - 1, min(qe, len(seq))
    else:
        sseq = revcomp(seq)
        scan_s, scan_e = L - qe, L - (qs - 1)
        scan_e = min(scan_e, L)
    span_nt = max(0, scan_e - scan_s)
    # 三读框同链扫描: 病毒读框 + 另外两个读框 (后者 = 本地读框基线 / 组成本底 null)
    stops = scan_frame_stops(sseq, scan_s, scan_e, offset)
    other = []
    for off in range(3):
        if off == offset:
            continue
        other.append(len(scan_frame_stops(sseq, scan_s, scan_e, off)))
    other_mean = round(sum(other) / len(other), 2) if other else ""
    prof = stop_profile(stops, span_nt)
    enrich = round((len(stops) + 1) / (other_mean + 1), 3) if other_mean != "" else ""
    cls = classify_decay(len(stops), other_mean, prof, span_nt)
    return {"frame_switches": switches, "n_hsp_frames": len(set(frames)),
            "premature_stops_region": len(stops),
            "hit_bitscore": max(h["bits"] for h in sub),
            "hit_label": sub[0].get("label", "")[:80],
            "hsp_span": f"{qs}-{qe}", "span_len": span_nt,
            "stops_other_mean": other_mean, "stop_enrichment": enrich,
            "head_stops": prof["head_stops"],
            "first_stop_frac": prof["first_stop_frac"],
            "max_stopfree_frac": prof["max_stopfree_frac"],
            "decay_class": cls}


def main():
    fa, hits_tsv, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
    seqs = read_fasta(fa)
    hits = defaultdict(list)
    with open(hits_tsv) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or line.startswith("qseqid"):
                continue
            try:
                hits[p[0]].append({
                    "sid": p[1], "pident": float(p[2]), "qstart": int(p[4]),
                    "qend": int(p[5]), "evalue": float(p[8]), "bits": float(p[9]),
                    "frame": p[10], "label": p[11] if len(p) > 11 else ""})
            except (ValueError, IndexError):
                continue
    with open(out_tsv, "w") as out:
        out.write(S1_HDR)
        for cid, seq in seqs.items():
            r = orf_scan(seq)
            h = hsp_metrics(hits.get(cid, []), seq)
            out.write(f"{cid}\t{len(seq)}\t{r['dom_frame']}\t{r['main_orf_aa']}\t"
                      f"{r['orf_max_fraction']}\t{r['coding_density']}\t{r['stops_dom_per_kb']}\t"
                      f"{h['frame_switches']}\t{h['n_hsp_frames']}\t{h['premature_stops_region']}\t"
                      f"{h['hit_bitscore']}\t{h['hit_label']}\t{h['hsp_span']}\t{h['span_len']}\t"
                      f"{h['stops_other_mean']}\t{h['stop_enrichment']}\t{h['head_stops']}\t"
                      f"{h['first_stop_frac']}\t{h['max_stopfree_frac']}\t{h['decay_class']}\n")
    print(f"[done] s1(logic v6.3, frame baseline): {len(seqs)} contigs -> {out_tsv}")


if __name__ == "__main__":
    main()
