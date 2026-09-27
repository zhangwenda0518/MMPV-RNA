#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
null_model.py — 判别器的空模型特异性检验（免外部真值）
======================================================
要回答的问题: s1 的"化石退化"签名与 s2 的"组件齐全"判据都没有被校准过。我们只证明了
代码按定义执行，没证明这个定义能把**真信号**从**组成噪声**里分出来。空模型就是补这一步。

做法: 把候选序列洗牌（严格保持碱基组成、破坏密码子结构与读框连续性），拿洗牌序列
**重跑同一套脚本、同一套阈值**，然后比较判定分布:

    P(判据 | 洗牌) = 只看组成、没有任何编码结构时该判据的触发概率

这是"诱饵库"式的标准阴性对照。真实富集倍数 = P(判据|真实) / P(判据|洗牌):
接近 1 说明该判据只是组成噪声的代理；越大越有信息量。

两种洗牌粒度
------------
  contig  整条 contig 单核苷酸洗牌 —— 端到端空模型（要重跑 diamond，最接近"随机序列
          进这条流水线会得到什么"）。这是主口径。
  region  只洗同源区、保留 contig 两侧上下文与坐标 —— 细粒度诊断：同源区**内部**的
          信号是否超出组成噪声。不需要 diamond，秒级。

用法
----
  # 1) 端到端空模型
  python3 null_model.py emit --run-dir run/ --mode contig -o shuf.fa
  #   然后把 shuf.fa 当作候选跑一遍 run_all.sh（或手跑各 stage），得到 *_shuf.tsv
  python3 null_model.py compare --stage s2 --real run/s2_domains.tsv \
      --null shuf_out/s2_domains.tsv
  python3 null_model.py compare --stage s3 --real run/eve_distinguish_verdict.tsv \
      --null shuf_out/eve_distinguish_verdict.tsv

  # 2) 同源区细粒度诊断（无需外部工具）
  python3 null_model.py region --run-dir run/ --n-shuffle 20 --json region.json
"""
import argparse
import collections
import json
import os
import pathlib
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import s1_decay_scan as s1          # noqa: E402


def tsv_load(path, key="contig_id"):
    """读带表头的 TSV -> {key: {列名: 值}} + 列名列表."""
    rows, hdr = {}, []
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f):
            p = line.rstrip("\n").split("\t")
            if i == 0:
                hdr = p
                continue
            if not line.strip():
                continue
            p += [""] * (len(hdr) - len(p))
            rows[p[hdr.index(key)] if key in hdr else p[0]] = dict(zip(hdr, p))
    return rows, hdr


def read_hits(path):
    hits = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return hits
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or p[0] == "qseqid":
                continue
            try:
                hits[p[0]].append({
                    "sid": p[1], "pident": float(p[2]), "length": int(p[3]),
                    "qstart": int(p[4]), "qend": int(p[5]),
                    "evalue": float(p[8]), "bits": float(p[9]),
                    "frame": p[10], "label": p[11] if len(p) > 11 else ""})
            except (ValueError, IndexError):
                continue
    return hits


def read_many(paths):
    out = collections.defaultdict(list)
    for p in paths:
        for k, v in read_hits(p).items():
            out[k].extend(v)
    return out


# ───────────────────────── 洗牌 (三种零模型) ─────────────────────────
# 单一零模型不够: 单核苷酸洗牌只保住碱基总组成, 同时破坏了**密码子位置特异组成**。
# 而 s1 的判据正是"病毒读框 vs 同链另两框"—— 三框在同一段碱基上, 位置特异组成一变,
# 三框的终止期望密度就跟着变。所以要按"保住得越来越多"分三档, 看结论是否零模型无关。
NULL_MODES = ("mono", "frame", "di")


def shuffle_nt(seq, rng):
    """mono: 单核苷酸洗牌 —— 保住碱基总组成, 破坏密码子结构与位置特异组成."""
    a = list(seq)
    rng.shuffle(a)
    return "".join(a)


def shuffle_frame(seq, rng):
    """frame: **保读框**洗牌 —— 第 1/2/3 密码子位各自独立洗.

    这档对 s1 是最锋利的对照: 位置特异组成与总组成都被严格保住, 三框的终止期望密度
    完全一致。于是任何残余富集只能归因于**密码子顺序**本身, 而不是组成差异。
    """
    cols = [list(seq[i::3]) for i in range(3)]
    for c in cols:
        rng.shuffle(c)
    out = []
    for i in range(max(len(c) for c in cols)):
        for c in cols:
            if i < len(c):
                out.append(c[i])
    return "".join(out)


def shuffle_di(seq, rng):
    """di: 双核苷酸保持洗牌 (Altschul-Erikson 欧拉路径法).

    保住全部二核苷酸频率 -> 终止密码子(三核苷酸)的期望频率也被约束住, 比 mono 更严。
    做法: 以二核苷酸为边建多重图, 随机走一条欧拉路径 (Hierholzer)。

    注意 Hierholzer 的迭代写法: **栈就是路径**, 只有在无路可走时才把顶点收进 path,
    最后整体反转。曾经写错过 —— 边走边往 out 里塞顶点、回溯时不补, 于是侧环被接在
    了错误的位置, 走不完的边被落下, 长度凑巧还对得上, 检查不出来。这种错误只能靠
    "逐次试验断言二核苷酸多重集不变"来抓。
    """
    if len(seq) < 4:
        return shuffle_nt(seq, rng)
    adj = collections.defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        adj[a].append(b)
    for v in adj:
        rng.shuffle(adj[v])
    stack, path = [seq[0]], []
    while stack:
        v = stack[-1]
        if adj.get(v):
            stack.append(adj[v].pop())
        else:
            path.append(stack.pop())
    path.reverse()
    res = "".join(path)
    if len(res) != len(seq):
        # 原始序列本身就是一条欧拉路径, 所以理论上走不完不该发生。真发生就**明确告警**,
        # 绝不静默退回 mono —— 那会让一个不保二核苷酸的序列冒充 di 洗牌结果。
        sys.stderr.write("[null_model] WARN: di 洗牌未覆盖全长, 退回 mono (结论请改用 "
                         "mono 或 frame 零模型复核)\n")
        return shuffle_nt(seq, rng)
    return res


def shuffler(seq, rng, mode):
    return {"mono": shuffle_nt, "frame": shuffle_frame, "di": shuffle_di}[mode](seq, rng)


def emit(run_dir, mode, out_fa, seed, null_mode="mono"):
    """生成洗牌候选 fasta. contig=整条洗; region=只洗最佳同源区."""
    seqs = s1.read_fasta(os.path.join(run_dir, "q.fa"))
    rng = random.Random(seed)
    hits = read_many([os.path.join(run_dir, f)
                      for f in ("panel_hits.tsv", "baits_hits.tsv")]) \
        if mode == "region" else {}
    n_reg = 0
    with open(out_fa, "w", encoding="utf-8") as fo:
        for cid, seq in seqs.items():
            if mode == "contig":
                new = shuffler(seq, rng, null_mode)
            else:
                reg = align_region(hits.get(cid, []), seq)
                if reg is None:            # 没有可用同源区: 原样保留, 只洗有区间的那批
                    new = seq
                else:
                    sseq, s, e, _ = reg
                    sh = shuffler(sseq[s:e], rng, null_mode)
                    n_reg += 1
                    if sseq is seq:
                        new = seq[:s] + sh + seq[e:]
                    else:
                        L = len(seq)
                        new = seq[:L - e] + sh + seq[L - s:]
            assert len(new) == len(seq), "洗牌必须保长度"
            fo.write(f">{cid}\n")
            for i in range(0, len(new), 60):
                fo.write(new[i:i + 60] + "\n")
    extra = f", 其中 {n_reg} 条洗了同源区" if mode == "region" else ""
    print(f"[emit] mode={mode} null={null_mode}: {len(seqs)} 条候选 -> {out_fa}{extra}")


def align_region(hsps, seq):
    """最佳 subject 的 HSP 覆盖区间 -> (sseq, start, end, offset)."""
    if not hsps:
        return None
    best = max(hsps, key=lambda h: h["bits"])["sid"]
    sub = [h for h in hsps if h["sid"] == best]
    qs = min(min(h["qstart"], h["qend"]) for h in sub)
    qe = max(max(h["qstart"], h["qend"]) for h in sub)
    fnum = int(collections.Counter(h["frame"] for h in sub).most_common(1)[0][0])
    L = len(seq)
    if fnum > 0:
        sseq, s, e = seq, qs - 1, min(qe, L)
    else:
        sseq = s1.revcomp(seq)
        s, e = L - qe, min(L - (qs - 1), L)
    if e - s < s1.MIN_SPAN:
        return None
    return sseq, s, e, abs(fnum) - 1


# ─────────────────── 同源区细粒度空模型（纯计算） ───────────────────
def classify_region(sseq, s, e, offset):
    """在给定区间跑 s1 的判据 (与 hsp_metrics 内部逐行同源)."""
    stops = s1.scan_frame_stops(sseq, s, e, offset)
    other = [len(s1.scan_frame_stops(sseq, s, e, o))
             for o in range(3) if o != offset]
    om = sum(other) / len(other) if other else 0.0
    prof = s1.stop_profile(stops, e - s)
    return (s1.classify_decay(len(stops), om, prof, e - s),
            prof["head_stops"], prof["max_stopfree_frac"],
            round((len(stops) + 1) / (om + 1), 3), e - s)


def region_null(run_dir, k, seed, json_out, null_mode="mono", quiet=False):
    seqs = s1.read_fasta(os.path.join(run_dir, "q.fa"))
    hits = read_many([os.path.join(run_dir, f) for f in
                      ("panel_hits.tsv", "baits_hits.tsv")])
    rng = random.Random(seed)
    real_c, null_c = collections.Counter(), collections.Counter()
    rows, n_used, n_nohit, n_short = [], 0, 0, 0
    for cid, seq in seqs.items():
        reg = align_region(hits.get(cid, []), seq)
        if reg is None:
            if hits.get(cid):
                n_short += 1
            else:
                n_nohit += 1
            continue
        n_used += 1
        sseq, s, e, off = reg
        rc, rh, rm, re_, span = classify_region(sseq, s, e, off)
        real_c[rc] += 1
        region = sseq[s:e]
        nd = collections.Counter()
        for _ in range(k):
            sh = shuffler(region, rng, null_mode)
            nc = classify_region(sh, 0, len(sh), off)[0]
            nd[nc] += 1
            null_c[nc] += 1
        rows.append({"contig": cid, "span": span, "real": rc,
                     "real_enrich": re_, "real_msf": rm,
                     "p_decay_null": round(nd.get("distributed_decay", 0) / k, 3),
                     "p_intact_null": round(nd.get("coding_intact", 0) / k, 3)})
    agg = _agg(real_c, null_c, k, n_used, n_nohit, n_short)
    agg["null_mode"] = null_mode
    if not quiet:
        _report(agg)
    if json_out:
        with open(json_out, "w", encoding="utf-8") as f:
            json.dump({"agg": agg, "rows": rows}, f, ensure_ascii=False, indent=1)
        print(f"[null_model] 明细 -> {json_out}")
    return agg


KEYS = ["coding_intact", "compositional_noise", "assembly_breakpoint",
        "distributed_decay"]


def _agg(real_c, null_c, k, n_used, n_nohit, n_short):
    rt = sum(real_c.values()) or 1
    nt = sum(null_c.values()) or 1
    a = {"n_with_region": n_used, "n_no_hit": n_nohit, "n_region_short": n_short,
         "k": k, "real": {c: real_c.get(c, 0) for c in KEYS},
         "null": {c: null_c.get(c, 0) for c in KEYS},
         "real_frac": {c: round(real_c.get(c, 0) / rt, 4) for c in KEYS},
         "null_frac": {c: round(null_c.get(c, 0) / nt, 4) for c in KEYS}}
    for c in ("distributed_decay", "coding_intact"):
        rf, nf = a["real_frac"][c], a["null_frac"][c]
        a[c + "_excess"] = round(rf / nf, 2) if nf else None
    return a


def _report(a):
    print("=" * 74)
    print("s1 结构通道 · 同源区空模型 (零模型=%s)" % a.get("null_mode", "mono"))
    print("  mono=保碱基组成 | frame=保读框(位置特异组成也保住) | di=保双核苷酸")
    print("=" * 74)
    print(f"有可用同源区 {a['n_with_region']} 条"
          f"（无命中 {a['n_no_hit']}, 区间过短 {a['n_region_short']}），每条洗牌 {a['k']} 次\n")
    print(f"{'判定':<24}{'真实':>7}{'真实占比':>11}{'空模型':>9}{'空模型占比':>12}{'富集':>8}")
    for c in KEYS:
        rf, nf = a["real_frac"][c], a["null_frac"][c]
        r = ("%.2f" % (rf / nf)) if nf else "inf"
        print(f"{c:<24}{a['real'][c]:>7}{rf:>11.4f}{a['null'][c]:>9}{nf:>12.4f}{r:>8}")
    print()
    for c, lab in (("distributed_decay", "化石退化"), ("coding_intact", "编码保留")):
        e = a[c + "_excess"]
        print(f"{lab}: 真实/空模型 = {e if e is not None else 'inf'}"
              + ("   ← 与纯组成噪声无法区分" if e is not None and e < 1.5 else ""))
    print("\n读法: 富集≈1 表示该判定只是组成噪声的代理, 没有额外信息量。")


# ─────────────────── 端到端空模型: 比较两次运行的产物 ───────────────────
CMP = {
    "s1": ("decay_class", ["coding_intact", "compositional_noise",
                           "assembly_breakpoint", "distributed_decay", "no_hsp"]),
    "s2": ("provirus_scale", None),
    "s3": ("verdict", None),
}


def compare(stage, real_f, null_f, json_out):
    real, hdr = tsv_load(real_f)
    null, _ = tsv_load(null_f)
    common = sorted(set(real) & set(null))
    print("=" * 74)
    print(f"端到端空模型 · stage={stage}")
    print("=" * 74)
    print(f"真实 {len(real)} 条 / 洗牌 {len(null)} 条 / 可比对 {len(common)} 条\n")
    if stage == "s1":
        col, order = CMP["s1"]
        rc = collections.Counter(real[c].get(col, "") for c in common)
        nc = collections.Counter(null[c].get(col, "") for c in common)
        print(f"{'decay_class':<24}{'真实':>7}{'占比':>9}{'洗牌':>8}{'占比':>9}{'富集':>8}")
        for k in order:
            rf = rc.get(k, 0) / max(len(common), 1)
            nf = nc.get(k, 0) / max(len(common), 1)
            print(f"{k:<24}{rc.get(k, 0):>7}{rf:>9.3f}{nc.get(k, 0):>8}{nf:>9.3f}"
                  f"{(('%.2f' % (rf / nf)) if nf else 'inf'):>8}")
        res = {"stage": stage, "real": dict(rc), "null": dict(nc)}
    elif stage == "s2":
        def flag(cid):
            return (real[cid].get("best_locus_ncomp", "0"),
                    real[cid].get("best_locus_comps", "-"))
        rn = collections.Counter(int(real[c]["best_locus_ncomp"] or 0) for c in common)
        nn = collections.Counter(int(null[c]["best_locus_ncomp"] or 0) for c in common)
        print(f"{'best_locus_ncomp':<20}{'真实':>7}{'洗牌':>8}")
        for k in sorted(set(rn) | set(nn)):
            print(f"{k:<20}{rn.get(k, 0):>7}{nn.get(k, 0):>8}")
        rp = sum(1 for c in common if real[c]["provirus_scale"] == "TRUE")
        np_ = sum(1 for c in common if null[c]["provirus_scale"] == "TRUE")
        print(f"\nprovirus_scale TRUE: 真实 {rp} / 洗牌 {np_}"
              f"  (富集 {('%.2f' % (rp / np_)) if np_ else 'inf'})")
        print(f"任一组件 (ncomp>=1): 真实 {len(common) - rn.get(0, 0)}"
              f" / 洗牌 {len(common) - nn.get(0, 0)}")
        res = {"stage": stage, "ncomp_real": dict(rn), "ncomp_null": dict(nn),
               "provirus_real": rp, "provirus_null": np_}
    else:
        rv = collections.Counter(real[c]["verdict"] for c in common)
        nv = collections.Counter(null[c]["verdict"] for c in common)
        tot = max(len(common), 1)
        print(f"{'verdict':<30}{'真实':>7}{'占比':>9}{'洗牌':>8}{'占比':>9}{'富集':>8}")
        for k in sorted(set(rv) | set(nv)):
            rf, nf = rv.get(k, 0) / tot, nv.get(k, 0) / tot
            print(f"{k:<30}{rv.get(k, 0):>7}{rf:>9.3f}{nv.get(k, 0):>8}{nf:>9.3f}"
                  f"{(('%.2f' % (rf / nf)) if nf else 'inf'):>8}")
        res = {"stage": stage, "real": dict(rv), "null": dict(nv)}
    if json_out:
        with open(json_out, "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=False, indent=1)
        print(f"\n[null_model] -> {json_out}")
    return res




# ─────────────────── 已知答案对照 ───────────────────
STOPS = {"TAA", "TAG", "TGA"}
_BASES = "ACGT"


def _rand_no_stop(rng, n_codons):
    """造一段**无终止**的随机编码序列 (frame 1 干净)."""
    out = []
    for _ in range(n_codons):
        while True:
            c = "".join(rng.choice(_BASES) for _ in range(3))
            if c not in STOPS:
                out.append(c)
                break
    return "".join(out)


def _make_case(work, kind, n, rng, span=900):
    """造一个合成 run-dir: q.fa + panel_hits.tsv, 每条 contig 一个覆盖全长的 HSP.

    kind="orf"    同源区 frame 1 无终止 (已知编码) -> 判据必须报 coding_intact
    kind="random" 同源区纯随机, 三框无差别 (已知无信号) -> 判据不该报富集
    """
    d = pathlib.Path(work) / kind
    d.mkdir(parents=True, exist_ok=True)
    flank = "".join(rng.choice(_BASES) for _ in range(60))
    with open(d / "q.fa", "w", encoding="utf-8") as fq, \
            open(d / "panel_hits.tsv", "w", encoding="utf-8") as fh:
        for i in range(n):
            if kind == "orf":
                region = _rand_no_stop(rng, span // 3)
            else:
                region = "".join(rng.choice(_BASES) for _ in range(span))
            seq = flank + region + flank
            cid = "%s_%03d" % (kind, i)
            fq.write(">%s\n%s\n" % (cid, seq))
            # 一个 HSP: 覆盖 region 全长, frame +1, bits 高, 比对长度足
            fh.write("%s\tSID_x\t90\t%d\t%d\t%d\t1\t%d\t1e-50\t500\t1\tlabel\n"
                     % (cid, len(region) // 3, 61, 60 + len(region),
                        len(region) // 3))
    return str(d)


def control(work, n, k):
    """已知答案对照 —— 跑任何 compare/region 结论之前先过这一关.

    构造两类**答案由构造保证**的输入:
      orf     frame 1 无终止密码子 => "编码受保留"为真, 判据必须看得见 (富集显著 > 1)
      random  纯随机          => 无信号为真, 判据不该声称富集 (富集 ~ 1)

    仪器连这两种都分不开的话, 它对真实数据说的任何话都不算数。
    """
    rng = random.Random(20260924)
    print("=" * 74)
    print("null_model 已知答案对照 (or 编码 vs 纯随机)")
    print("=" * 74)
    res = {}
    for kind in ("orf", "random"):
        rd = _make_case(work, kind, n, rng)
        agg = region_null(rd, k, 20260924, "", null_mode="mono", quiet=True)
        rf, nf = agg["real_frac"]["coding_intact"], agg["null_frac"]["coding_intact"]
        exc = (rf / nf) if nf else float("inf")
        res[kind] = {"real_frac": rf, "null_frac": nf, "excess": exc}
        print("  %-7s coding_intact: 真实 %.3f / 空模型 %.3f -> 富集 %s"
              % (kind, rf, nf, ("%.2f" % exc) if nf else "inf"))
    ok_pos = res["orf"]["real_frac"] >= 0.5 and res["orf"]["excess"] >= 3.0
    ok_neg = res["random"]["excess"] <= 3.0
    print()
    print("  阳性对照 (已知编码必须被看见): %s" % ("通过" if ok_pos else "**失败**"))
    print("  阴性对照 (纯随机不该有富集):   %s" % ("通过" if ok_neg else "**失败**"))
    if not (ok_pos and ok_neg):
        print("")
        print("!! 仪器通不过自己的对照, 它对真实数据给出的结论一律不可用。")
        return 1
    print("")
    print("对照通过: 该仪器能分开'已知编码'与'已知随机', 可以读它的结论。")
    return 0

def main():
    ap = argparse.ArgumentParser(description="EVE 判别器空模型检验")
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("emit", help="生成洗牌候选 fasta")
    e.add_argument("--run-dir", required=True)
    e.add_argument("--mode", default="contig", choices=["contig", "region"])
    e.add_argument("-o", "--out", required=True)
    e.add_argument("--seed", type=int, default=20260924)
    e.add_argument("--null-mode", default="mono", choices=NULL_MODES)
    r = sub.add_parser("region", help="同源区细粒度空模型 (纯计算)")
    r.add_argument("--run-dir", required=True)
    r.add_argument("--n-shuffle", type=int, default=20)
    r.add_argument("--seed", type=int, default=20260924)
    r.add_argument("--json", default="", dest="json_out")
    r.add_argument("--null-mode", default="mono", choices=NULL_MODES)
    r.add_argument("--quiet", action="store_true", help="不打印明细 (给 control 用)")
    c = sub.add_parser("compare", help="比较真实 vs 洗牌两次运行的产物")
    c.add_argument("--stage", required=True, choices=["s1", "s2", "s3"])
    c.add_argument("--real", required=True)
    c.add_argument("--null", required=True)
    c.add_argument("--json", default="", dest="json_out")
    k = sub.add_parser("control", help="已知答案对照: 先确认仪器看得见已知信号")
    k.add_argument("--work", required=True, help="临时目录 (自己造合成输入)")
    k.add_argument("--n", type=int, default=60, help="每类造多少条")
    k.add_argument("--n-shuffle", type=int, default=20)
    a = ap.parse_args()
    if a.cmd == "emit":
        emit(a.run_dir, a.mode, a.out, a.seed, a.null_mode)
    elif a.cmd == "region":
        region_null(a.run_dir, a.n_shuffle, a.seed, a.json_out, a.null_mode, a.quiet)
    elif a.cmd == "control":
        return control(a.work, a.n, a.n_shuffle)
    else:
        compare(a.stage, a.real, a.null, a.json_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
