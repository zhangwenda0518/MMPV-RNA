#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_bam_coverage.py — 用 BAM 覆盖特征区分"整合"与"活跃病毒"

背景
----
MMPV 项目的 BAM 是**reads 比对到 5,679 条参考病毒基因组**的结果 (不是比对到
候选 contig)。所以不能直接给候选 contig 算覆盖。

但仍有可用路径: 如果一条候选能通过 blastn 映射到某个参考病毒, 就可以用
**该参考病毒的 BAM 覆盖特征**来判断这条病毒在样本里是"活跃感染"还是
"低水平/整合样":

  * 覆盖均匀性 (CV): 活跃病毒沿基因组均匀; 整合/低水平常有断崖
  * 覆盖深度: 活跃病毒高; 整合样极低
  * 链偏好: 整合位点常有链偏好 (单向转录); 活跃病毒双链均衡
  * 跨样本一致性: 同物种多样本都低覆盖 -> 更像内源

输入:
  --bam-dir    BAM 目录 (含 <sample>.sorted.bam)
  --sample-map 样本号 -> 参考病毒 accession 的映射 (来自 blastn 结果)
  --ref-fa     参考病毒 fasta (取长度)
输出: <out>/bam_coverage.tsv
"""
import argparse
import subprocess
import sys
from pathlib import Path


def bam_depth(bam, region, samtools="samtools"):
    """samtools depth 一个区域 -> 每碱基深度列表。"""
    p = subprocess.run([samtools, "depth", "-r", region, bam],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return []
    out = []
    for line in p.stdout.splitlines():
        c = line.split("\t")
        if len(c) >= 3:
            try:
                out.append(int(c[2]))
            except ValueError:
                pass
    return out


def coverage_stats(depths):
    """深度列表 -> (均值, CV, 零覆盖比例, 最大/中位比)。"""
    if not depths:
        return {"mean": 0.0, "cv": 0.0, "zero_frac": 1.0, "max_med_ratio": 0.0}
    n = len(depths)
    mean = sum(depths) / n
    if mean == 0:
        return {"mean": 0.0, "cv": 0.0, "zero_frac": 1.0, "max_med_ratio": 0.0}
    var = sum((d - mean) ** 2 for d in depths) / n
    cv = (var ** 0.5) / mean
    zero = sum(1 for d in depths if d == 0) / n
    s = sorted(depths)
    med = s[n // 2] or 1
    return {"mean": round(mean, 3), "cv": round(cv, 3),
            "zero_frac": round(zero, 3), "max_med_ratio": round(s[-1] / med, 2)}


def strand_bias(bam, region, samtools="samtools"):
    """正/负链 reads 数 -> (forward, reverse, bias)。

    bias = |f-r|/(f+r); 接近 0 = 双链均衡 (活跃病毒), 接近 1 = 单链偏好 (整合样)。
    """
    p = subprocess.run([samtools, "view", "-c", "-f", "64", region, bam],
                       capture_output=True, text=True)
    f = int(p.stdout.strip() or 0) if p.returncode == 0 else 0
    p = subprocess.run([samtools, "view", "-c", "-f", "128", region, bam],
                       capture_output=True, text=True)
    r = int(p.stdout.strip() or 0) if p.returncode == 0 else 0
    tot = f + r
    bias = round(abs(f - r) / tot, 3) if tot else 0.0
    return f, r, bias


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam-dir", required=True)
    ap.add_argument("--sample-map", required=True,
                    help="两列 TSV: 样本号 <TAB> 参考病毒 accession")
    ap.add_argument("--out", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    bam_dir = Path(a.bam_dir)

    pairs = []
    with open(a.sample_map, encoding="utf-8") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2 and c[0] and c[1]:
                pairs.append((c[0], c[1]))
    if a.limit:
        pairs = pairs[:a.limit]
    print(f"[load] {len(pairs)} 个 样本-参考 对", file=sys.stderr)

    rows = []
    for i, (sample, acc) in enumerate(pairs, 1):
        bam = bam_dir / f"{sample}.sorted.bam"
        if not bam.is_file():
            continue
        depths = bam_depth(str(bam), acc, a.samtools)
        st = coverage_stats(depths)
        f, r, bias = strand_bias(str(bam), acc, a.samtools)
        rows.append((sample, acc, len(depths), st["mean"], st["cv"],
                     st["zero_frac"], st["max_med_ratio"], f, r, bias))
        if i % 20 == 0:
            print(f"  {i}/{len(pairs)}", file=sys.stderr)

    tsv = out / "bam_coverage.tsv"
    hdr = ("sample\tref_acc\tn_bases\tmean_depth\tcv\tzero_frac\t"
           "max_med_ratio\tforward\treverse\tstrand_bias\n")
    with open(tsv, "w", encoding="utf-8", newline="\n") as fo:
        fo.write(hdr)
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")
    print(f"\n[done] {tsv}  ({len(rows)} 行)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
