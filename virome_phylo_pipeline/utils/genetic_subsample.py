#!/usr/bin/env python3
"""
utils/genetic_subsample.py — 遗传信息驱动的两段式降采样（FPS / clade）+ 敏感性对照
================================================================================
补齐 geo_analysis.smart_subsample 的遗传盲区：现有 unique_timeloc 按
(地点×时间窗) 取**最长**序列——**遗传上任意的**（可能留 5 条几乎相同的、
丢掉差异大的，而地理信号恰藏在遗传差异里）。本模块做两段式：

  段1（已有思想，保留）: 按 (地点 × 时间窗) 分层 —— 保护时空多样性
  段2（新增）          : cell 内 FPS 最远点采样 / 全局 clade 完全连锁聚类 —— 保护遗传多样性

方法学出处（逐条可追溯）
----------------------
* FPS（最远点采样，组内 max_reps 个遗传上最分散的代表；种子=含 N 最少）与
  clade（全局 SNP 距离 + 完全连锁聚类 cutoff，每 地点/时间/clade 留 1 代表）：
  git-repo/phymap-workflow/scripts/prep_metadata.py:335-506 prune_sequences 的
  两种方法（思想与判据）；实现为 numpy 向量化重写，不引入 pairsnp 依赖。
* SNP 距离口径：只统计**双端均为无歧义碱基**的位点差异（phymap 同款）；
  U 归一为 T（RNA 病毒比对常见）；N/gap 不计差异也不计相同。
* ⚠ 类病毒/小基因组尺度警告（borrow 报告 §批次3）：phymap 默认阈值
  （clade cutoff=5 SNP、min_snp_diff=2）按病原基因组级多样性设定；
  PSTVd 359bp 全基因组成对差异可能仅几十 SNP，直接照抄会"全并一簇"。
  阈值必须按数据重标定并写入 METHODS.md——本模块默认值沿用 phymap 原值
  但在 QC 输出中标注 'needs_calibration'。
* 敏感性对照（zika_Vietnam snakemake/py/subsample.py 的"全量 vs 子样本"框架）：
  对全量/子采样各算 S（分离位点数）、π（平均 p-distance）、成对 SNP 距离分位数，
  输出对照表，声明降采样对多样性估计的影响。

自包含性：numpy/pandas/Bio.SeqIO + 标准库；不动 geo_analysis 既有行为。
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from utils.virphy_bridge import LogCollector
except ImportError:  # 允许直接运行
    from .virphy_bridge import LogCollector  # type: ignore

_VALID_BASES = np.array([ord(c) for c in "ACGT"], dtype=np.uint8)


# ─────────────────────────────────────────────────────────────────────
# SNP 距离
# ─────────────────────────────────────────────────────────────────────

def read_fasta(path: str) -> List:
    from Bio import SeqIO
    recs = list(SeqIO.parse(path, "fasta"))
    if not recs:
        raise SystemExit(f"FASTA 无序列: {path}")
    return recs


def _encode(recs: List) -> Tuple[np.ndarray, List[str]]:
    """序列 → uint8 矩阵（U→T 归一, 大写），返回 (n × L) 矩阵与 id 列表。"""
    ids, arrs = [], []
    for r in recs:
        s = str(r.seq).upper().replace("U", "T")
        arrs.append(np.frombuffer(s.encode("ascii"), dtype=np.uint8))
        ids.append(r.id)
    lengths = {a.size for a in arrs}
    if len(lengths) > 1:
        raise SystemExit(f"序列未比对（长度不一: {sorted(lengths)}）——本模块要求已比对 FASTA")
    return np.vstack(arrs), ids


def snp_distance_matrix(recs: List) -> Tuple[List[str], np.ndarray]:
    """成对 SNP 距离矩阵（双端无歧义位点口径；phymap compute_snp_distance_matrix
    的 numpy 路径等价实现，向量化为一次性 n×n）。"""
    seqs, ids = _encode(recs)
    n = len(ids)
    valid = np.isin(seqs, _VALID_BASES)              # n × L
    d = np.zeros((n, n), dtype=np.int32)
    for i in range(n):
        both = valid[i] & valid[i + 1:]               # (n-i-1) × L
        diff = (seqs[i] != seqs[i + 1:]) & both
        d[i, i + 1:] = diff.sum(axis=1)
    d = d + d.T
    return ids, d


def n_counts(recs: List) -> Dict[str, int]:
    return {r.id: str(r.seq).upper().count("N") for r in recs}


# ─────────────────────────────────────────────────────────────────────
# FPS / clade
# ─────────────────────────────────────────────────────────────────────

def fps_select(group_ids: List[str], dist: np.ndarray, id_index: Dict[str, int],
               max_reps: int, min_snp_diff: int, ncount: Dict[str, int]) -> List[str]:
    """组内最远点采样（phymap prune_sequences/fps 判据）：
    种子 = 含 N 最少；每轮加入"与已选集合最小距离最大"的候选；
    候选最小距离 < min_snp_diff 时停止。"""
    if len(group_ids) <= max_reps:
        return list(group_ids)
    idx = [id_index[g] for g in group_ids]
    selected = [min(group_ids, key=lambda g: ncount.get(g, 0))]
    sel_idx = [id_index[selected[0]]]
    while len(selected) < max_reps:
        best_cand, best_min_d = None, -1
        for g in group_ids:
            if g in selected:
                continue
            gi = id_index[g]
            min_d = int(min(dist[gi][si] for si in sel_idx))
            if min_d > best_min_d:
                best_min_d, best_cand = min_d, g
        if best_cand is None or best_min_d < min_snp_diff:
            break
        selected.append(best_cand)
        sel_idx.append(id_index[best_cand])
    return selected


def clade_clusters(ids: List[str], dist: np.ndarray, cutoff: int) -> List[List[str]]:
    """完全连锁聚类（phymap prune_sequences/clade 判据）：
    簇间距离 = 跨簇任意两成员距离的**最大值**；最小簇间距离 ≤ cutoff 则合并。
    朴素 O(n³)，n>2000 时建议改用 fps（调用方守卫）。"""
    clusters = [[i] for i in range(len(ids))]
    while True:
        best = (None, None, float("inf"))
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                max_d = 0
                for i in clusters[a]:
                    row = dist[i]
                    for j in clusters[b]:
                        if row[j] > max_d:
                            max_d = row[j]
                if max_d < best[2]:
                    best = (a, b, max_d)
        a, b, min_max_d = best
        if a is None or min_max_d > cutoff:
            break
        clusters[a].extend(clusters[b])
        clusters.pop(b)
    return [[ids[i] for i in c] for c in clusters]


# ─────────────────────────────────────────────────────────────────────
# 分组
# ─────────────────────────────────────────────────────────────────────

def _date_group(date_str, resolution: str = "month") -> str:
    val = str(date_str or "").strip()
    if not val or val.lower() in ("nan", "none", "null"):
        return "unknown"
    val = val.replace("/", "-")
    parts = val.split("-")
    if not parts or len(parts[0]) != 4:
        return "unknown"
    year = parts[0]
    month = parts[1].zfill(2) if len(parts) >= 2 else "07"
    if resolution == "year":
        return year
    return f"{year}-{month}"


def _load_meta(metadata_csv: str) -> pd.DataFrame:
    df = pd.read_csv(metadata_csv, encoding="utf-8-sig")
    cols = {c.strip().lower(): c for c in df.columns}
    name_col = cols.get("name") or list(df.columns)[0]
    loc_col = cols.get("location")
    date_col = cols.get("date") or cols.get("collection_date")
    out = df[[name_col] + ([loc_col] if loc_col else []) + ([date_col] if date_col else [])].copy()
    out.columns = ["name"] + ([c for c in ("location", "date") if c]) [:len(out.columns) - 1]
    return out


# ─────────────────────────────────────────────────────────────────────
# 编排
# ─────────────────────────────────────────────────────────────────────

def genetic_subsample(
    fasta: str,
    metadata_csv: Optional[str] = None,
    output_fasta: str = "genetic_subsampled.fasta",
    mode: str = "two_stage",
    max_reps: int = 3,
    min_snp_diff: int = 2,
    clade_cutoff: int = 5,
    date_resolution: str = "month",
    log: Optional[LogCollector] = None,
) -> Dict:
    """两段式降采样主入口。

    mode:
      two_stage — (地点[×时间窗]) 分层 + cell 内 FPS（默认；段1保护时空，段2保护遗传）
      clade     — 全局完全连锁聚类 + 每(地点×时间×clade)留 1 代表（n>2000 自动降级 two_stage）
      none      — 仅按 ID 去重（等价不降采样）
    """
    log = log or LogCollector()
    recs = read_fasta(fasta)
    rec_map = {r.id: r for r in recs}
    meta = _load_meta(metadata_csv) if metadata_csv else pd.DataFrame(columns=["name"])
    meta_by_name = {str(r["name"]).strip(): r for _, r in meta.iterrows()}

    matched = [r.id for r in recs if r.id in meta_by_name] if metadata_csv else [r.id for r in recs]
    if metadata_csv and len(matched) < len(recs):
        log.warning(f"{len(recs) - len(matched)} 条 FASTA 序列无元数据行（不参与降采样）")
    if metadata_csv and not matched:
        raise SystemExit("FASTA 与元数据零交集")
    work_ids = matched

    ids, dist = snp_distance_matrix([rec_map[i] for i in work_ids])
    id_index = {g: i for i, g in enumerate(ids)}
    ncount = n_counts([rec_map[i] for i in work_ids])

    meta_of = lambda gid: meta_by_name.get(gid, {})  # noqa: E731
    retained: List[str] = []
    removed_rows: List[Dict] = []

    if mode == "clade" and len(work_ids) > 2000:
        log.warning(f"n={len(work_ids)} > 2000, clade 完全连锁 O(n³) 过慢 → 自动降级 two_stage")
        mode = "two_stage"

    if mode == "clade":
        clusters = clade_clusters(ids, dist, clade_cutoff)
        log.emit(f"clade: {len(ids)} 序列 → {len(clusters)} 个 clade (cutoff={clade_cutoff} SNP)")
        seq_to_clade = {m: ci for ci, c in enumerate(clusters) for m in c}
        groups: Dict[Tuple, List[str]] = {}
        for g in ids:
            m = meta_of(g)
            loc = str(m.get("location", "") or "").strip()
            dg = _date_group(m.get("date"), date_resolution)
            groups.setdefault((loc, dg, seq_to_clade[g]), []).append(g)
        for key, members in sorted(groups.items(), key=lambda kv: str(kv[0])):
            keep = min(members, key=lambda g: ncount.get(g, 0))
            retained.append(keep)
            for g in members:
                if g != keep:
                    removed_rows.append({"removed_id": g, "kept_id": keep, "group": str(key),
                                         "reason": f"clade member (cutoff={clade_cutoff})"})
    elif mode == "two_stage":
        groups: Dict[Tuple, List[str]] = {}
        for g in ids:
            m = meta_of(g)
            loc = str(m.get("location", "") or "").strip()
            dg = _date_group(m.get("date"), date_resolution) if "date" in m else "all"
            groups.setdefault((loc, dg), []).append(g)
        for key, members in sorted(groups.items(), key=lambda kv: str(kv[0])):
            if len(members) <= 1:
                retained.extend(members)
                continue
            sub_ids = [g for g in members]
            sub_dist = dist[np.ix_([id_index[g] for g in sub_ids], [id_index[g] for g in sub_ids])]
            sub_index = {g: i for i, g in enumerate(sub_ids)}
            keep = fps_select(sub_ids, sub_dist, sub_index, max_reps, min_snp_diff, ncount)
            retained.extend(keep)
            for g in members:
                if g not in keep:
                    removed_rows.append({"removed_id": g, "kept_id": keep[0], "group": str(key),
                                         "reason": f"FPS (max_reps={max_reps}, min_snp_diff={min_snp_diff})"})
    else:  # none
        retained = list(ids)

    # 落盘 FASTA（保持原文件顺序）
    from Bio import SeqIO
    keep_set = set(retained)
    kept_recs = [r for r in recs if r.id in keep_set]
    os.makedirs(os.path.dirname(os.path.abspath(output_fasta)), exist_ok=True)
    SeqIO.write(kept_recs, output_fasta, "fasta")

    qc = {
        "mode": mode,
        "n_input": len(recs),
        "n_matched": len(work_ids),
        "n_retained": len(kept_recs),
        "n_removed": len(removed_rows),
        "reduction_pct": round(100 * (1 - len(kept_recs) / max(1, len(recs))), 1),
        "params": {"max_reps": max_reps, "min_snp_diff": min_snp_diff,
                   "clade_cutoff": clade_cutoff, "date_resolution": date_resolution},
        "needs_calibration": True,
        "note": "阈值沿用 phymap 默认（病原基因组级多样性设定）；小基因组(如 PSTVd 359bp)须先按本数据成对 SNP 分布重标定",
    }
    out = {"success": True, "output_fasta": output_fasta, "qc": qc,
           "removed": removed_rows, "retained": retained}
    log.emit(f"genetic_subsample[{mode}]: {qc['n_input']} → {qc['n_retained']} "
             f"({qc['reduction_pct']}% 减少)")
    return out


# ─────────────────────────────────────────────────────────────────────
# 敏感性对照（zika_Vietnam 全量 vs 子样本框架）
# ─────────────────────────────────────────────────────────────────────

def compare_diversity(full_fasta: str, sub_fasta: str, out_csv: Optional[str] = None,
                      log: Optional[LogCollector] = None) -> pd.DataFrame:
    """全量 vs 降采样多样性对照：S(分离位点)/π(平均 p-distance)/成对 SNP 分位数。"""
    log = log or LogCollector()
    rows = []
    for label, path in (("full", full_fasta), ("subsampled", sub_fasta)):
        recs = read_fasta(path)
        ids, dist = snp_distance_matrix(recs)
        iu = np.triu_indices(len(ids), k=1)
        pairs = dist[iu].astype(float)
        seqs, _ = _encode(recs)
        variable_cols = int(sum(1 for c in range(seqs.shape[1])
                                if len(set(seqs[:, c].tolist())) > 1
                                and np.isin(seqs[:, c], _VALID_BASES).sum() > 0))
        l = seqs.shape[1]
        pi = float((pairs / l).mean()) if pairs.size else 0.0
        rows.append({
            "dataset": label, "n_seqs": len(ids), "segregating_sites": variable_cols,
            "pi_pdistance": round(pi, 6),
            "pair_snp_median": float(np.median(pairs)) if pairs.size else 0,
            "pair_snp_p95": float(np.percentile(pairs, 95)) if pairs.size else 0,
            "pair_snp_max": float(pairs.max()) if pairs.size else 0,
        })
    df = pd.DataFrame(rows)
    if out_csv:
        os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
        df.to_csv(out_csv, index=False)
    d_pi = (df.loc[1, "pi_pdistance"] - df.loc[0, "pi_pdistance"]) / max(df.loc[0, "pi_pdistance"], 1e-12)
    log.emit(f"敏感性: π 变化 {d_pi:+.1%}, S {df.loc[0,'segregating_sites']}→{df.loc[1,'segregating_sites']}")
    return df


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="遗传信息驱动降采样（两段式: 时空分层 + cell 内 FPS / 全局 clade）+ 敏感性对照")
    ap.add_argument("--fasta", required=True, help="已比对 FASTA")
    ap.add_argument("--metadata", default=None, help="name,location[,date] CSV")
    ap.add_argument("-o", "--output", default="genetic_subsampled.fasta")
    ap.add_argument("--mode", default="two_stage", choices=["two_stage", "clade", "none"])
    ap.add_argument("--max-reps", type=int, default=3)
    ap.add_argument("--min-snp-diff", type=int, default=2)
    ap.add_argument("--clade-cutoff", type=int, default=5)
    ap.add_argument("--date-resolution", default="month", choices=["year", "month"])
    ap.add_argument("--sensitivity", default=None,
                    help="对照报告输出 CSV 路径（给了则自动对全量 vs 子采样算多样性）")
    ap.add_argument("--qc-json", default=None, help="QC 输出路径（默认与 output 同目录）")
    args = ap.parse_args()
    log = LogCollector()
    r = genetic_subsample(args.fasta, args.metadata, args.output, mode=args.mode,
                          max_reps=args.max_reps, min_snp_diff=args.min_snp_diff,
                          clade_cutoff=args.clade_cutoff, date_resolution=args.date_resolution,
                          log=log)
    qc_path = args.qc_json or os.path.splitext(args.output)[0] + "_qc.json"
    with open(qc_path, "w", encoding="utf-8") as f:
        json.dump(r["qc"], f, ensure_ascii=False, indent=2)
    if args.sensitivity:
        compare_diversity(args.fasta, args.output, out_csv=args.sensitivity, log=log)
    print("GENETIC_SUBSAMPLE_DONE")


if __name__ == "__main__":
    main()
