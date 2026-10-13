#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tabulate.py — 基准集的逐类判定矩阵（判别器的灵敏度/特异性）
============================================================
读 s3 的 verdict 表 + 类别映射, 输出"类别 × verdict"矩阵, 并按键的严重程度标出
**方向性错误**:

  假阴性 (漏): 真的整合元件 (B1) 被判成病毒侧 (virus_candidate / virus_fragment_review)
               -> 把 EVE 当现行病毒留下, 下游会漏清
  假阳性 (误): 现存病毒 (A) 或寄主基因 (C) 被判成 EVE 侧
               (ancient_EVE / EVE_*.STRONG / EVE_LTR_TE -> MOVE_EVE; host_contamination -> REMOVE)
               -> 会删/移错数据, 是最贵的一类错

用法:
  python3 tabulate.py --verdict WORK/eve_distinguish_verdict.tsv \
                      --filter  WORK/dna_vs_eve_filter.tsv \
                      --map     BENCH/class_map.tsv [--json OUT]
"""
import argparse
import collections
import csv
import json
import os
import sys

EVE_SIDE = {"ancient_EVE", "EVE_STRONG_provirus", "EVE_LTR_TE", "EVE_suspect"}
VIRUS_SIDE = {"virus_candidate"}
HOST_SIDE = {"host_contamination_likely", "host_conflict_review",
             "host_homology_cross_species"}


def load_map(path):
    m = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            m[r["id"]] = r["class"]
    return m


def load_tsv(path, key):
    out = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[r[key]] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verdict", required=True)
    ap.add_argument("--filter", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    cmap = load_map(a.map)
    vd = load_tsv(a.verdict, "contig_id") if os.path.exists(a.verdict) else {}
    fl = load_tsv(a.filter, "contig_id") if os.path.exists(a.filter) else {}

    classes = ["A_extant_virus", "B1_eprv_clean", "B1b_eprv_decayed",
               "B2_eprv_activatable", "C_host_viral_domain"]
    per = collections.defaultdict(collections.Counter)
    acts = collections.defaultdict(collections.Counter)
    missing = collections.Counter()
    for cid, cls in cmap.items():
        if cid not in vd:
            missing[cls] += 1
            per[cls]["(无判定)"] += 1
            continue
        per[cls][vd[cid].get("verdict", "")] += 1
        acts[cls][fl.get(cid, {}).get("action", "")] += 1

    print("=" * 78)
    print("策展基准集 × 判别器 (寄主通道关闭, 与 sensitivity 同口径)")
    print("=" * 78)
    for cls in classes:
        if not per[cls]:
            continue
        tot = sum(per[cls].values())
        print("\n%s  (%d 条)" % (cls, tot))
        for v, n in per[cls].most_common():
            print("    %-30s %5d  %5.1f%%" % (v, n, 100.0 * n / max(tot, 1)))
        if acts[cls]:
            print("    action: " + "  ".join("%s=%d" % (k, v)
                                             for k, v in acts[cls].most_common()))

    # 方向性错误要分**两个严重度**报, 混在一起会把风险讲大:
    #   consequential  会真的动数据 (MOVE_EVE / REMOVE_host_contamination)
    #   flagged        只是"判成 EVE 侧但转人工复核" (EVE_suspect -> REVIEW), 不删不移
    print("\n" + "=" * 78)
    print("方向性错误 (分两个严重度)")
    print("=" * 78)
    err = {}
    for cls in classes:
        tot = sum(per[cls].values())
        if not tot:
            continue
        moved = acts[cls].get("MOVE_EVE", 0) + acts[cls].get("REMOVE_host_contamination", 0)
        flagged = sum(n for v, n in per[cls].items() if v in EVE_SIDE) - \
            acts[cls].get("MOVE_EVE", 0)
        if cls in ("B1_eprv_clean", "B1b_eprv_decayed"):
            print("  %-22s 灵敏度(真的 MOVE_EVE 掉): %3d / %3d = %5.1f%%"
                  % (cls, moved, tot, 100.0 * moved / tot))
            print("  %-22s 其余 %d 条留在病毒集里 (review/%s 等) —— 漏清, 但不删错"
                  % ("", tot - moved, "virus_fragment_review"))
            err[cls] = {"sensitivity_move": moved, "total": tot,
                        "rate": round(100.0 * moved / tot, 2)}
        elif cls in ("A_extant_virus", "C_host_viral_domain"):
            cons = acts[cls].get("MOVE_EVE", 0) + acts[cls].get("REMOVE_host_contamination", 0)
            print("  %-22s 会真动数据的假阳性: %3d / %3d = %5.1f%%   (MOVE_EVE+REMOVE)"
                  % (cls, cons, tot, 100.0 * cons / tot))
            print("  %-22s 仅判成 EVE 侧转人工(不删不移): %3d / %3d = %5.1f%%"
                  % ("", flagged, tot, 100.0 * flagged / tot))
            err[cls] = {"consequential_fp": cons, "flagged_only": flagged, "total": tot,
                        "consequential_rate": round(100.0 * cons / tot, 2)}
        else:
            print("  %-22s 标签本身有歧义 (可激活 EPRV), 只报分布不作对错" % cls)
    if missing:
        print("\n  注意: 有 %d 条基准序列没进 verdict 表 (%s)"
              % (sum(missing.values()), dict(missing)))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump({"per_class": {k: dict(v) for k, v in per.items()},
                       "actions": {k: dict(v) for k, v in acts.items()},
                       "errors": err}, f, ensure_ascii=False, indent=1)
        print("\n-> %s" % a.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
