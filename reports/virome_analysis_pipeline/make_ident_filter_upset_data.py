#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_ident_filter_upset_data.py — 为「7 轨道 x 过滤状态」堆叠 UpSet 汇总数据。

口径（2026-10-13 探针实测确认，见 _probe_ident_filter.py 输出）：
  集合（7）= 6 个有命中的病毒工具 + viroid_blastn
     viralverify / metabuli / blast / genomad / virbot / rdrpcatch  +  viroid_blastn
  逐样本产物：
     02a_Identification/{S}/{S}_virus.{tool}.result.id     工具命中 id（每行一个 contig）
     02a_Identification/{S}/{S}_virus.all.result.id        6 病毒工具并集（= 过滤链的输入）
     02a_Identification/{S}/{S}_viroids.result.id          viroid_blastn 命中 id
     02b_Filter/{S}/cdd/all_candidates_cdd_filter.fasta    CDD 过滤后保留（⊆ union）
     02b_Filter/{S}/uniprot/all_candidates_uniprot.fasta   Uniprot 过滤后保留（⊆ cdd_kept）
  实测链式嵌套（3 样本 + 全域均成立）：union ⊇ cdd_kept ⊇ retained，且 viroid ∩ union = 0。

状态（互斥）：
  retained  = 在 uniprot 保留集（最终保留）
  cdd_only  = 通过 CDD 但被 Uniprot 滤掉
  removed   = 被 CDD 滤掉（未进入 cdd_kept）
  viroid    = viroid_blastn 命中（不在 union 内，**不参与**病毒过滤链）
  nofilter  = 该样本无 02b uniprot 产物（无法判定过滤结果）

产出（聚合后极小，只落「签名 x 状态计数」）：
  ident_filter_upset_signature.tsv
  ident_filter_upset_meta.json
"""
import collections
import glob
import json
import os
import sys
from pathlib import Path

BASE = Path(os.path.expanduser(
    "~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"))
A = BASE / "02a_Identification"
B = BASE / "02b_Filter"
OUT = Path(os.path.expanduser("~/ident_filter_upset"))

TOOLS = ["viralverify", "metabuli", "blast", "genomad", "virbot", "rdrpcatch"]
VIROID = "viroid_blastn"
ALL7 = TOOLS + [VIROID]

STATUSES = ["retained", "cdd_only", "removed", "viroid", "nofilter"]


def read_id_file(path):
    out = set()
    with open(path) as f:
        for line in f:
            c = line.strip()
            if c:
                out.add(c)
    return out


def read_fasta_ids(path):
    out = set()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                out.add(line[1:].split()[0])
    return out


def main():
    samples = sorted(d.name for d in A.iterdir()
                     if d.is_dir() and glob.glob(str(A / d.name / "*virus.all.result.id")))
    print(f"样本数 = {len(samples)}")

    sig_counter = collections.defaultdict(collections.Counter)  # signature -> Counter(status)
    n_samples_no_filter = []
    sum_union = sum_viroid = sum_retained = sum_cdd = 0
    viroid_in_union = 0

    for i, s in enumerate(samples, 1):
        ad = A / s
        sets = {}
        ok = True
        for t in TOOLS:
            g = glob.glob(str(ad / f"*virus.{t}.result.id"))
            if not g:
                ok = False
                break
            sets[t] = read_id_file(g[0])
        if not ok:
            print(f"  [skip] {s} 缺工具产物", file=sys.stderr)
            continue
        gv = glob.glob(str(ad / "*viroids.result.id"))
        vset = read_id_file(gv[0]) if gv else set()
        union = read_id_file(glob.glob(str(ad / "*virus.all.result.id"))[0])

        # 自检：并集 = 各工具重算并集
        recomputed = set().union(*sets.values())
        if recomputed != union:
            print(f"  [warn] {s} 并集不一致: id文件={len(union)} 重算={len(recomputed)}",
                  file=sys.stderr)
            union = recomputed

        viroid_in_union += len(vset & union)

        up = B / s / "uniprot" / "all_candidates_uniprot.fasta"
        cdd = B / s / "cdd" / "all_candidates_cdd_filter.fasta"
        retained = read_fasta_ids(up) if up.is_file() else None
        cdd_kept = read_fasta_ids(cdd) if cdd.is_file() else set()
        if retained is None:
            n_samples_no_filter.append(s)
            retained = set()
            has_filter = False
        else:
            has_filter = True

        # 自检：嵌套
        if has_filter and not (retained <= cdd_kept <= union):
            print(f"  [warn] {s} 嵌套失败 retained⊆cdd⊆union: "
                  f"{retained <= cdd_kept} / {cdd_kept <= union}", file=sys.stderr)

        # 逐 contig 归类（union ∪ viroid）
        # ⚠ 状态优先级：**先看是否在 union 内**（= 进了病毒过滤链），只有 union 之外
        #    的 viroid 命中才算 "viroid"。实测全域 viroid∩union = 88 ≠ 0（早先只在
        #    3 个样本上抽验得 0，不能推广）→ 若把 viroid 判在前面，这 88 条会被
        #    错误地记成 "viroid"，其过滤结局丢失。
        for contig in (union | vset):
            bits = "".join("1" if contig in sets[t] else "0" for t in TOOLS)
            bits += "1" if contig in vset else "0"
            if contig in union:
                if not has_filter:
                    sig_counter[bits]["nofilter"] += 1
                elif contig in retained:
                    sig_counter[bits]["retained"] += 1
                elif contig in cdd_kept:
                    sig_counter[bits]["cdd_only"] += 1
                else:
                    sig_counter[bits]["removed"] += 1
            else:
                sig_counter[bits]["viroid"] += 1

        sum_union += len(union)
        sum_viroid += len(vset)
        sum_retained += len(retained)
        sum_cdd += len(cdd_kept)

        if i % 50 == 0:
            print(f"  ... {i}/{len(samples)}", flush=True)

    # ── 落盘 ──
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for sig, c in sig_counter.items():
        rows.append((sig, sum(c.values()), c["retained"], c["cdd_only"],
                     c["removed"], c["viroid"], c["nofilter"]))
    rows.sort(key=lambda r: (-r[1], r[0]))

    tsv = OUT / "ident_filter_upset_signature.tsv"
    with tsv.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("signature\tn_total\tn_retained\tn_cdd_only\tn_removed\tn_viroid\tn_nofilter\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    meta = dict(
        source_ident="02a_Identification/{S}/{S}_virus.{tool}.result.id",
        source_union="02a_Identification/{S}/{S}_virus.all.result.id",
        source_viroid="02a_Identification/{S}/{S}_viroids.result.id",
        source_cdd="02b_Filter/{S}/cdd/all_candidates_cdd_filter.fasta",
        source_uniprot="02b_Filter/{S}/uniprot/all_candidates_uniprot.fasta",
        tool_order=ALL7,
        n_samples=len(samples),
        n_samples_no_filter=len(n_samples_no_filter),
        samples_no_filter=n_samples_no_filter,
        n_signatures=len(rows),
        n_sequences=sum(r[1] for r in rows),
        sum_union=sum_union,
        sum_viroid=sum_viroid,
        sum_retained=sum_retained,
        sum_cdd_kept=sum_cdd,
        viroid_in_union=viroid_in_union,
        retained_pct_of_union=round(100.0 * sum_retained / sum_union, 2),
        status_legend=dict(
            retained="passed CDD + UniProt filter (final kept)",
            cdd_only="passed CDD, removed by UniProt",
            removed="removed by CDD",
            viroid="viroid_blastn hit (separate track, outside the virus filter chain)",
            nofilter="sample has no 02b UniProt output - filter outcome unknown",
        ),
    )
    with (OUT / "ident_filter_upset_meta.json").open("w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)

    print()
    print(f"✅ {tsv}")
    print(f"   签名数 {len(rows)}｜序列总数 {meta['n_sequences']:,}")
    print(f"   Σunion {sum_union:,}｜Σviroid {sum_viroid:,}｜Σretained {sum_retained:,}"
          f"（{meta['retained_pct_of_union']}%）｜Σcdd_kept {sum_cdd:,}")
    print(f"   viroid∩union = {viroid_in_union}（>0 → 类病毒轨道与病毒轨道确有重叠）")
    print(f"   无过滤产物样本 {len(n_samples_no_filter)} 个: {n_samples_no_filter}")
    tot = {k: sum(c[k] for c in sig_counter.values()) for k in STATUSES}
    print("   状态合计: " + "  ".join(f"{k}={v:,}" for k, v in tot.items()))
    print()
    print("   Top 12 签名:")
    for r in rows[:12]:
        print(f"     {r[0]}  total={r[1]:>7,}  ret={r[2]:>6,}  cdd_only={r[3]:>7,}"
              f"  removed={r[4]:>7,}  viroid={r[5]:>4,}  nofilter={r[6]:>5,}")


if __name__ == "__main__":
    main()
