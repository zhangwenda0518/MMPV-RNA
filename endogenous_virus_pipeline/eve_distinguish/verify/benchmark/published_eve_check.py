#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
published_eve_check.py — 独立复算"已发表 EVE 集"的双侧判定，与上游结论对帐
================================================================================
背景: 上游 hi-fever 研究对 Liu et al. 2026 (Mol Plant) 公开的植物 EVE 集做过一次
独立验证, 结论是 **3,063 条全部不通过双侧判定、96.6% 归为 host_like, 成分是植物
自身的 RDR (RdRp 类) 基因 / 线粒体 ORF / RT 域蛋白**。那个结论是那份研究的对外材料
之一, 但它本身只是一次跑出来的数字 —— 按本项目一贯的纪律, **别人给的数字不能直接引用,
要么自己复算一遍, 要么把它标成未核实**。

本脚本做的就是复算: 用同一套库、同一套规则重新跑一遍, 然后与上游的逐条判定表对帐,
输出**一致率**。一致就引用该结论; 不一致就必须查出分歧在哪。

复算链路 (与 `eve_scan_core.verdict_loci` 同构):
  1. 输入 = 已发表植物 EVE 的**蛋白**序列 (上游已把核酸 6 框翻译好)
  2. diamond blastp vs Stage-2 双侧库 (`plant_virus_final.dmnd`)
  3. 按 `id2div_final.tsv` 把命中参考分成 viral / plant 两侧
  4. 取两侧最佳 bitscore, 套规则:
       viral_bs >= 50 且 viral_bs >= plant_bs  -> viral_supported
       plant_bs >= 50 且 plant_bs >  viral_bs  -> host_like
       否则                                     -> undetermined
  注意: **Stage-1 病毒参考库那一遍不参与判定** —— 上游文档也写明规则只落在双侧库上
  (病毒参考库那遍只是用来报"85% 有病毒命中")。所以这里也不跑它, 否则复算的不是同一个东西。

用法:
  python3 published_eve_check.py --fa mp_plant_eve_p.fa \\
      --pv-db plant_virus_final.dmnd --id2div id2div_final.tsv \\
      --theirs mp_eve_verdict.tsv --work WORKDIR [--threads 32]
"""
import argparse
import collections
import os
import subprocess
import sys

HOST_BS = 50
VIRAL_BS = 50


def load_id2div(path):
    """sseqid -> 'viral'|'plant'. 与 eve_scan_core.load_id2div 同口径 (制表符两列)."""
    m = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                m[p[0]] = p[1]
    return m


def run_blastp(fa, db, out, threads, diamond="diamond"):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print("  复用已有命中表 %s" % out)
        return
    cmd = [diamond, "blastp", "-q", fa, "-d", db, "-o", out,
           "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue",
           "bitscore", "stitle", "-e", "1e-5", "--max-target-seqs", "20",
           "--threads", str(threads)]
    print("  跑 diamond blastp: %s" % " ".join(cmd[:6] + ["..."]))
    subprocess.run(cmd, check=True)


def call_two_sided(hits_tsv, id2div):
    """-> {qseqid: (verdict, viral_bs, plant_bs, viral_hit, plant_hit)}"""
    best = collections.defaultdict(lambda: [0.0, "", 0.0, ""])   # v_bs, v_id, p_bs, p_id
    with open(hits_tsv, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                bs = float(p[5])
            except ValueError:
                continue
            q, sid = p[0], p[1]
            dv = id2div.get(sid, "?")
            b = best[q]
            if dv == "viral" and bs > b[0]:
                b[0], b[1] = bs, sid
            elif dv == "plant" and bs > b[2]:
                b[2], b[3] = bs, sid
    out = {}
    for q, (vb, vid, pb, pid) in best.items():
        if vb >= VIRAL_BS and vb >= pb:
            v = "viral_supported"
        elif pb >= HOST_BS and pb > vb:
            v = "host_like"
        else:
            v = "undetermined"
        out[q] = (v, vb, pb, vid, pid)
    return out


def load_theirs(path):
    """读上游判定表 -> {eve_id: (verdict, viral_bs, plant_bs)}"""
    out = {}
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        vi = hdr.index("verdict")
        ki = hdr.index("eve_id") if "eve_id" in hdr else 0
        bi = hdr.index("viral_bs") if "viral_bs" in hdr else None
        pi = hdr.index("plant_bs") if "plant_bs" in hdr else None
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= vi:
                continue
            vb = float(p[bi]) if bi is not None and p[bi] not in ("", "NA") else 0.0
            pb = float(p[pi]) if pi is not None and p[pi] not in ("", "NA") else 0.0
            out[p[ki]] = (p[vi], vb, pb)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fa", required=True, help="已发表植物 EVE 的蛋白 fasta")
    ap.add_argument("--pv-db", required=True)
    ap.add_argument("--id2div", required=True)
    ap.add_argument("--theirs", required=True, help="上游逐条判定表 (用于对帐)")
    ap.add_argument("--work", required=True)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)

    id2div = load_id2div(a.id2div)
    print("id2div 条目: %d" % len(id2div))
    hits = os.path.join(a.work, "pub_vs_pvfinal.tsv")
    run_blastp(a.fa, a.pv_db, hits, a.threads)
    mine = call_two_sided(hits, id2div)
    theirs = load_theirs(a.theirs)
    n_hit = sum(1 for _ in open(hits, encoding="utf-8"))
    print("blastp 命中行: %d; 复算覆盖序列: %d; 上游表: %d" % (n_hit, len(mine), len(theirs)))

    common = sorted(set(mine) & set(theirs))
    print("可对帐: %d 条" % len(common))
    mism = [q for q in common if mine[q][0] != theirs[q][0]]
    print("\n=== 判定一致率 ===")
    if common:
        print("  一致 %d / %d = %.2f%%" % (len(common) - len(mism), len(common),
                                          100.0 * (len(common) - len(mism)) / len(common)))
    cm = collections.Counter(mine[q][0] for q in common)
    ct = collections.Counter(theirs[q][0] for q in common)
    print("  复算分布: %s" % dict(cm))
    print("  上游分布: %s" % dict(ct))
    if mism:
        print("\n=== 分歧样例 (最多 5 条) ===")
        for q in mism[:5]:
            print("  %s\n    复算 %s (v=%.0f p=%.0f) / 上游 %s (v=%.0f p=%.0f)"
                  % (q[:60], mine[q][0], mine[q][1], mine[q][2],
                     theirs[q][0], theirs[q][1], theirs[q][2]))
    # 关键结论的复现检查
    n_vs_mine = sum(1 for q in common if mine[q][0] == "viral_supported")
    n_vs_theirs = sum(1 for q in common if theirs[q][0] == "viral_supported")
    print("\n=== 上游关键结论: '0 条通过双侧判定' ===")
    print("  复算 viral_supported: %d 条; 上游表: %d 条 -> %s"
          % (n_vs_mine, n_vs_theirs,
             "复现" if (n_vs_mine == 0 and n_vs_theirs == 0) else "**不一致, 需查**"))
    if a.json:
        import json
        with open(a.json, "w", encoding="utf-8", newline="\n") as fo:
            json.dump({"n_common": len(common), "n_mismatch": len(mism),
                       "agree_pct": round(100.0 * (len(common) - len(mism))
                                          / max(len(common), 1), 2),
                       "mine_dist": dict(cm), "theirs_dist": dict(ct),
                       "viral_supported_mine": n_vs_mine,
                       "viral_supported_theirs": n_vs_theirs}, fo,
                      ensure_ascii=False, indent=1)
        print("-> %s" % a.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
