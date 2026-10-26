#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
regen_barbarum_detection_report.py — 宁夏枸杞 (Lycium barbarum) 已知病毒检测报告重生成

背景:
  known_virus_barbarum/01_detection/summary/analysis_report.txt 原为指向
  known_virus_all_v2 的 symlink, 内容是跨物种合并口径 (处理样本 1348 / 确诊 754 /
  分母隐含 1400), 与单物种 barbarum 完全不符。

本脚本复用 batch_virus_depth.py::generate_report_txt 的报告结构与统计口径
(平均 CPM / FPKM / 覆盖度 / Poisson Ratio / 群体检出率), 但:
  - 分母 = barbarum 权威样本集 (01_Assembly 目录对应的 run), 而非 v2 合并队列
  - 仅统计 best.summary 中实际检出的病毒
  - 额外补充 v2 版本没有的: 检出矩阵规模 / 每样本病毒数分布 / 覆盖度分级

用法:
  python regen_barbarum_detection_report.py \
      -v known_virus_barbarum/01_detection/summary/all_viruses.best.summary.tsv \
      -d known_virus_barbarum/01_detection/summary \
      --denominator-sources 02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly
"""
import argparse
import collections
import csv
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# 分子类型 -> 归类 (用于报告里的 DNA/RNA 分节)
DNA_HINTS = ("viroid", "mitovirus", "nucleovirus", "geminivirus", "tobamovirus",
             "caulimovirus", "circular ssDNA", "dsDNA")
RNA_HINTS = ("virus", "viroid", "rna")


def mol_class(name):
    n = name.lower()
    if "viroid" in n:
        return "Viroid"
    if any(h in n for h in ("DNA virus", "circular ssDNA", "dsDNA")):
        return "DNA virus"
    return "RNA virus"


def mean_of(rows, key):
    vals = []
    for r in rows:
        raw = (r.get(key) or "").strip()
        if not raw:
            continue
        try:
            vals.append(float(raw))
        except ValueError:
            continue
    return (sum(vals) / len(vals)) if vals else 0.0


def count_assembly_runs(asm_root):
    """从 novel_virus/01_Assembly 目录统计权威 run 数 (权威数据来源)。"""
    root = Path(asm_root)
    if not root.is_dir():
        return None
    runs = set()
    for p in root.iterdir():
        if not p.is_dir() or p.name.lower() in ("logs", "log"):
            continue
        # 目录名可能直接是 run, 也可能带前后缀
        tok = p.name
        for pref in ("sample_", "SRR", "ERR", "DRR", "CRR"):
            if tok.startswith(pref):
                break
        runs.add(tok)
    return len(runs)


def build_report(best_path, denominator, denom_desc, engine="SALMON"):
    with open(best_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise SystemExit(f"[FATAL] best.summary 为空: {best_path}")

    sample_col = "Sample" if "Sample" in rows[0] else list(rows[0].keys())[0]
    sp_col = "Adjusted_Species" if "Adjusted_Species" in rows[0] else list(rows[0].keys())[2]

    samples = sorted({(r.get(sample_col) or "").strip() for r in rows} - {""})
    by_virus = collections.defaultdict(list)
    for r in rows:
        by_virus[(r.get(sp_col) or "").strip()].append(r)

    # 每样本检出病毒数
    per_sample = collections.Counter()
    for r in rows:
        s = (r.get(sample_col) or "").strip()
        if s:
            per_sample[s] += 1

    out = []
    W = 80
    out.append("=" * W)
    out.append("宏病毒定量全景报告 (V42.4 Engine · Lycium barbarum 单物种口径)")
    out.append("=" * W)
    out.append("")
    out.append(f"  核心引擎: {engine} (极速伪比对)")
    out.append(f"  数据来源: 01_detection/summary/all_viruses.best.summary.tsv")
    out.append(f"  分母口径: {denominator} 个样本  <- {denom_desc}")
    out.append(f"  检出记录: {len(rows)} 条 | 涉及样本: {len(samples)} 个")
    out.append(f"  检出病毒: {len(by_virus)} 种")
    out.append(f"  打假策略: Poisson Ratio 建模 (覆盖比值阈值 >= 0.5, 已在best 阶段固化)")
    out.append("")

    # ---- 队列概览 ----
    out.append("【队列构成】")
    out.append("-" * 60)
    neg = sum(1 for s in samples if per_sample[s] == 0)
    multi = sum(1 for s in samples if per_sample[s] >= 2)
    out.append(f"  队列总样本(分母): {denominator}")
    out.append(f"  其中检出病毒  : {len(samples)} ({len(samples)/denominator*100:.1f}%)")
    out.append(f"  未检出任何病毒: {denominator - len(samples)} "
               f"({(denominator-len(samples))/denominator*100:.1f}%)")
    out.append(f"  单一感染      : {len(samples) - multi}")
    out.append(f"  合并感染(>=2) : {multi}")
    out.append("")

    # ---- 确诊群落概览 (逐病毒, 按检出数降序) ----
    out.append("【确诊群落概览】(按检出例数降序)")
    out.append("-" * 60)
    order = sorted(by_virus.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    for name, grp in order:
        n = len(grp)
        rate = n / denominator * 100
        mc = mol_class(name)
        out.append(f"[{mc}] {name}")
        out.append(f"    检出 {n} 例 (群体检出率 {rate:.1f}%)")
        out.append(f"    |- 平均 CPM  : {mean_of(grp,'Asm_CPM'):.2f}")
        out.append(f"    |- 平均 FPKM : {mean_of(grp,'Asm_FPKM'):.2f}")
        # ANI / Pi 仅传统比对 (coverM) 模式才有; 伪比对下整列为空, 不打印
        ani, pi = mean_of(grp, "Avg_Read_ANI"), mean_of(grp, "Avg_Pi")
        if any((r.get("Avg_Read_ANI") or "").strip() for r in grp):
            out.append(f"    |- 测定 ANI  : {ani:.2f}%")
            out.append(f"    |- 测定 Pi   : {pi:.5f}")
        out.append(f"    |- 覆盖度    : {mean_of(grp,'Rep_Coverage(%)'):.2f}%")
        out.append(f"    `- 泊松比值  : {mean_of(grp,'Poisson_Ratio'):.2f}")
        out.append("")

    # ---- 高流行率病毒汇总 (>=5%) ----
    hi = [(nm, len(g)) for nm, g in order if len(g) / denominator >= 0.05]
    out.append("【高流行率病毒 (检出率 >= 5%)】")
    out.append("-" * 60)
    if hi:
        for nm, n in hi:
            out.append(f"  {nm:56s} {n:5d} 例  {n/denominator*100:5.1f}%")
    else:
        out.append("  (无)")
    out.append("")

    # ---- 共感染结构 ----
    out.append("【共感染结构】")
    out.append("-" * 60)
    dist = collections.Counter(per_sample[s] for s in samples)
    for k in sorted(dist):
        out.append(f"  检出 {k} 种病毒的样本: {dist[k]:5d} ({dist[k]/len(samples)*100:.1f}%)")
    out.append("")

    # ---- 覆盖度质量分级 ----
    out.append("【覆盖度质量分级 (按代表序列平均覆盖度)】")
    out.append("-" * 60)
    buckets = collections.Counter()
    for nm, g in order:
        c = mean_of(g, "Rep_Coverage(%)")
        if c >= 90:
            buckets[">=90% (完整)"] += 1
        elif c >= 70:
            buckets["70-90% (高)"] += 1
        elif c >= 50:
            buckets["50-70% (中)"] += 1
        else:
            buckets["<50% (低/碎片)"] += 1
    for k in [">=90% (完整)", "70-90% (高)", "50-70% (中)", "<50% (低/碎片)"]:
        if buckets.get(k):
            out.append(f"  {k:20s}: {buckets[k]:3d} 种")
    out.append("")

    out.append("=" * W)
    out.append("注: 本报告分母为 Lycium barbarum 单物种队列, 与 known_virus_all_v2 的")
    out.append("    跨物种合并口径 (v2+v3, 1400样本) 不可直接比较。")
    out.append("=" * W)
    return "\n".join(out) + "\n", by_virus, samples, per_sample


def main():
    ap = argparse.ArgumentParser(description="重生成 L. barbarum 已知病毒检测报告")
    ap.add_argument("-v", "--virus-summary", required=True, help="all_viruses.best.summary.tsv")
    ap.add_argument("-d", "--out-dir", required=True, help="输出目录 (01_detection/summary)")
    ap.add_argument("--denominator", type=int, default=None, help="显式指定分母")
    ap.add_argument("--denominator-sources", default=None,
                    help="novel_virus/01_Assembly 目录, 用于自动统计权威 run 数")
    ap.add_argument("--engine", default="SALMON")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    vsum = Path(a.virus_summary)
    if not vsum.is_file():
        raise SystemExit(f"[FATAL] 不存在: {vsum}")

    n_asm = count_assembly_runs(a.denominator_sources) if a.denominator_sources else None
    if a.denominator:
        denom, desc = a.denominator, "显式指定"
    elif n_asm:
        denom, desc = n_asm, f"权威组装目录 {a.denominator_sources} (run 去重)"
    else:
        raise SystemExit("[FATAL] 需--denominator 或 --denominator-sources 之一")

    text, by_virus, samples, per_sample = build_report(vsum, denom, desc, a.engine)

    if a.dry_run:
        sys.stdout.write(text)
        print(f"\n[dry-run] 未写盘。分母={denom} 病毒={len(by_virus)} 检出样本={len(samples)}")
        return

    outdir = Path(a.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    dest = outdir / "analysis_report.txt"
    # 覆盖前若目标是 symlink, 先解除链接, 避免写穿到 v2
    if dest.is_symlink():
        dest.unlink()
        print(f"[guard] 已解除 symlink: {dest} -> 不写穿 v2")
    dest.write_text(text, encoding="utf-8")
    print(f"[OK] 写出 {dest} ({len(text.encode('utf-8'))} B)")
    print(f"     分母={denom} | 检出样本={len(samples)} | 病毒={len(by_virus)}")


if __name__ == "__main__":
    main()