#!/usr/bin/env python3
"""S1 结构通道 v4: 编码退化签名 + 同批组装错误率基线.

v3 的缺陷 (2026-09-23 修订): premature_stops_region 用绝对计数 + 固定阈值 (stops>=3 判退化),
  会被组装 indel 污染 —— 一个单碱基 indel 就足以让其后全部密码子错读, 在远缘同源区
  (25-45% pid) 上造出十几个"提前终止", 与真正的化石退化无法区分.
v4 的做法: 不跟绝对阈值较劲, 而是把每个 contig **自身**当基线:

  * 同源区同一段序列、同一条链、另外两个读框的终止密度 = 本地 null.
    同一批组装、同一读段、同一碱基组成, 组装 indel 对三个读框的污染是同等的,
    所以另外两框的终止密度就代表了"这套组装的噪声底噪".
  * enrichment = 病毒读框终止数 / (另外两框均数 + 1): 远大于 1 说明病毒读框
    确实比噪声底更干净(编码受选择保留); 接近 1 说明该读框和随机无异 = 退化.
  * 停止分布剖面区分两种"接近 1":
      - max_stopfree_frac >= 0.25  -> 存在一段完整可读长区, 是**单个 indel 断裂**
      - head_stops >= 1            -> 终止从同源区最前端就有, 是**分布式退化**
    前者是组装噪声, 后者才是化石签名.

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
    # 最长无终止区间 (含两端截断)
    runs = [pos[0]] + [b - a for a, b in zip(pos, pos[1:])] + [span_nt - pos[-1]]
    max_stopfree_frac = round(max(runs) / span_nt, 4)
    return {"head_stops": head_stops, "first_stop_frac": first_stop_frac,
            "max_stopfree_frac": max_stopfree_frac}


def classify_decay(stops, other_mean, prof, span_nt):
    """v4 分类: 编码保留 / 分布式退化(化石签名) / 组装断裂 / 组成本底噪声."""
    if stops is None or span_nt < MIN_SPAN:
        return "no_hsp"
    if stops <= 1:
        return "coding_intact"
    enrich = (stops + 1) / (other_mean + 1)
    if enrich < 1.0:
        # 病毒读框并不比另外两框更脏 -> 终止来自序列组成/组装噪声, 不是退化证据
        return "compositional_noise"
    msf = prof["max_stopfree_frac"]
    head = prof["head_stops"]
    if msf != "" and msf >= 0.25:
        # 存在一段 >=1/4 同源区长度的完整无终止区 -> 单个 indel 断裂, 而非全段退化.
        # (断裂点两侧各是一段完整可读区, 所以最长无终止段必然占去大半)
        return "assembly_breakpoint"
    if head >= 1 or stops >= 3:
        # 终止从最前端就有, 或数量密集且无任何长完整段 -> 全段分布式退化 = 化石签名.
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
    # 三读框同链扫描: 病毒读框 + 另外两个读框 (后者 = 同批组装的本地噪声基线)
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
        out.write("contig_id\tlength\tdom_frame\tmain_orf_aa\torf_max_fraction\t"
                  "coding_density\tstops_dom_per_kb\tframe_switches\tn_hsp_frames\t"
                  "premature_stops_region\thit_bitscore\thit_label\t"
                  "hsp_span\tspan_len\tstops_other_mean\tstop_enrichment\thead_stops\t"
                  "first_stop_frac\tmax_stopfree_frac\tdecay_class\n")
        for cid, seq in seqs.items():
            r = orf_scan(seq)
            h = hsp_metrics(hits.get(cid, []), seq)
            out.write(f"{cid}\t{len(seq)}\t{r['dom_frame']}\t{r['main_orf_aa']}\t"
                      f"{r['orf_max_fraction']}\t{r['coding_density']}\t{r['stops_dom_per_kb']}\t"
                      f"{h['frame_switches']}\t{h['n_hsp_frames']}\t{h['premature_stops_region']}\t"
                      f"{h['hit_bitscore']}\t{h['hit_label']}\t{h['hsp_span']}\t{h['span_len']}\t"
                      f"{h['stops_other_mean']}\t{h['stop_enrichment']}\t{h['head_stops']}\t"
                      f"{h['first_stop_frac']}\t{h['max_stopfree_frac']}\t{h['decay_class']}\n")
    print(f"[done] s1(v4 baseline): {len(seqs)} contigs -> {out_tsv}")


if __name__ == "__main__":
    main()
