#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""coexp_module_conservation.py — 跨宿主共表达模块保守性 (OrthoFinder 直系组)
=========================================================================
方法学依据: KingdomStress 预印本 (bioRxiv 10.64898/2026.04.17.719161) Fig 5 —
跨物种模块对经直系组 (orthogroup, OG) 映射后, 重叠系数
    OC = |共享 OG 数| / min(模块 A 的 OG 数, 模块 B 的 OG 数)
判保守: OC >= 0.5 且共享 OG >= 10 (论文阈值, 可调)。零模型: 在每个物种内
打乱 基因→OG 映射 (--permutations 次, 种子固定), 对候选对重算 OC 得经验 p 值
p = (1 + #{null_OC >= obs_OC}) / (P + 1)。(论文用 Wilcoxon 符号秩对全部
配对做总体检验; 本脚本输出逐对经验 p, 可在下游自行做 Wilcoxon。)

用途: 同一 MMPV 项目跑多个宿主 (各一份 ⑧ 输出) 时, 比较两宿主的
共表达模块 (coexp_modules.py 的 module_table_<组>.tsv) 哪些跨宿主保守。
单宿主项目用不到 — 独立脚本, 不在编排器链上。

输入:
  --modules-a/-b   module_table_<组>.tsv (Gene, Module, Degree)
  --orthogroups    OrthoFinder Orthogroups.tsv (首列 Orthogroup, 其余列=物种名;
                   单元格内基因以逗号分隔 — 兼容空格/逗号混排)
  --species-a/-b   物种在 Orthogroups.tsv 里的列名 (默认取 --name-a/-b)
输出 (--output_dir):
  conservation_pairs.tsv   Name_A, Module_A, N_OG_A, Name_B, Module_B, N_OG_B,
                           Shared_OG, OC, Conserved, P_Emp, P_Emp_BH,
                           Significant_BH, Null_OC_Median
                           (P_Emp_BH = 对候选对经验 p 的 BH-FDR, 2026-10-10 加)
  conservation_summary.tsv 参数 + 计数 (含未映射基因数)

用法:
  python coexp_module_conservation.py \
      --modules-a hostA/05_Coexp/modules/module_table_positive.tsv \
      --modules-b hostB/05_Coexp/modules/module_table_positive.tsv \
      --name-a hostA --name-b hostB \
      --orthogroups OrthoFinder/Results/Orthogroups.tsv \
      --species-a hostA --species-b hostB --output_dir conservation/
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR (NaN 透传)。与 enrich_analysis/viral_load_scan 同式。"""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order] * ok.sum() / (np.arange(ok.sum()) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(int(ok.sum()))
    res[order] = np.clip(ranked, 0, 1)
    out[ok] = res
    return out


def read_module_table(path: Path) -> dict:
    """module_table tsv → {module_id: [gene, ...]}。"""
    df = pd.read_csv(path, sep="\t")
    if not {"Gene", "Module"}.issubset(df.columns):
        sys.exit(f"[ERROR] 模块表需含 Gene/Module 列: {path}")
    out = {}
    for mid, grp in df.groupby("Module"):
        out[int(mid)] = sorted(set(grp["Gene"].astype(str)))
    return out


def read_orthogroups(path: Path, species_col: str) -> dict:
    """Orthogroups.tsv → {gene: og} (指定物种列; 缺列报错)。"""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "Orthogroup" not in df.columns:
        sys.exit(f"[ERROR] Orthogroups.tsv 缺 Orthogroup 列: {path}")
    if species_col not in df.columns:
        sys.exit(f"[ERROR] Orthogroups.tsv 无物种列 {species_col!r}; "
                 f"可用: {list(df.columns)[:8]}...")
    g2og = {}
    for og, cell in zip(df["Orthogroup"], df[species_col]):
        cell = str(cell).strip()
        if not cell:
            continue
        for g in cell.replace(",", " ").split():
            g2og[g] = str(og)
    return g2og


def module_og_sets(modules: dict, g2og: dict):
    """模块 → OG 集合; 返回 (每模块 OG 集, 基因→OG 索引数组素材, 未映射基因数)。

    索引素材: (genes_list, og_list) — 模块存为 genes_list 中的下标数组,
    供置换时用 og_perm[idx] 快速重算。"""
    genes_list, gene_idx = [], {}
    for m, genes in modules.items():
        for g in genes:
            if g in g2og and g not in gene_idx:
                gene_idx[g] = len(genes_list)
                genes_list.append(g)
    og_list = np.array([g2og[g] for g in genes_list], dtype=object)
    sets, idxs, unmapped = {}, {}, 0
    for m, genes in modules.items():
        idx = np.array([gene_idx[g] for g in genes if g in g2og], dtype=int)
        unmapped += len(genes) - len(idx)
        idxs[m] = idx
        ogs = {g2og[g] for g in genes if g in g2og}
        sets[m] = ogs
    return sets, idxs, og_list, unmapped


def pair_stats(sets_a: dict, sets_b: dict):
    """全部模块对的 (shared, oc); 任一方 OG 数为 0 的对记 nan。"""
    rows = []
    for ma, sa in sets_a.items():
        for mb, sb in sets_b.items():
            shared = len(sa & sb)
            denom = min(len(sa), len(sb))
            oc = shared / denom if denom else float("nan")
            rows.append((ma, mb, len(sa), len(sb), shared, oc))
    return rows


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description="跨宿主共表达模块保守性 (KingdomStress Fig5 OC + 置换零模型)")
    ap.add_argument("--modules-a", required=True, help="宿主 A 模块长表")
    ap.add_argument("--modules-b", required=True, help="宿主 B 模块长表")
    ap.add_argument("--name-a", required=True, help="宿主 A 名称 (输出列)")
    ap.add_argument("--name-b", required=True, help="宿主 B 名称")
    ap.add_argument("--orthogroups", required=True,
                    help="OrthoFinder Orthogroups.tsv")
    ap.add_argument("--species-a", help="A 的物种列名 (默认 = --name-a)")
    ap.add_argument("--species-b", help="B 的物种列名 (默认 = --name-b)")
    ap.add_argument("--oc-cutoff", type=float, default=0.5,
                    help="保守判定 OC 阈值 (论文 0.5)")
    ap.add_argument("--min-shared-og", type=int, default=10,
                    help="保守判定最少共享 OG (论文 10); 置换仅算 shared≥此值的候选对")
    ap.add_argument("--permutations", type=int, default=100,
                    help="物种内 基因→OG 置换次数")
    ap.add_argument("--seed", type=int, default=42, help="置换随机种子")
    ap.add_argument("--output_dir", required=True, help="输出目录")
    args = ap.parse_args(argv)

    sp_a = args.species_a or args.name_a
    sp_b = args.species_b or args.name_b
    mods_a = read_module_table(Path(args.modules_a))
    mods_b = read_module_table(Path(args.modules_b))
    print(f"宿主 {args.name_a}: {len(mods_a)} 模块; 宿主 {args.name_b}: "
          f"{len(mods_b)} 模块")

    g2og_a = read_orthogroups(Path(args.orthogroups), sp_a)
    g2og_b = read_orthogroups(Path(args.orthogroups), sp_b)
    sets_a, idxs_a, og_a, unm_a = module_og_sets(mods_a, g2og_a)
    sets_b, idxs_b, og_b, unm_b = module_og_sets(mods_b, g2og_b)
    print(f"OG 映射: A 未映射基因 {unm_a}, B 未映射基因 {unm_b} "
          f"(不在 Orthogroups.tsv 中, 已剔除)")

    stats = pair_stats(sets_a, sets_b)
    cand = [(ma, mb, shared, oc) for ma, mb, _, _, shared, oc in stats
            if shared == shared and shared >= args.min_shared_og]
    print(f"模块对: {len(stats)} (有效候选 shared≥{args.min_shared_og}: "
          f"{len(cand)})")

    # 置换零模型: 物种内打乱 基因→OG, 重算候选对 OC
    rng = np.random.default_rng(args.seed)
    null_ge = {k: 0 for k, _ in enumerate(cand)}
    null_med = {k: [] for k, _ in enumerate(cand)}
    if cand and args.permutations > 0:
        for _ in range(args.permutations):
            perm_a = rng.permutation(og_a)
            perm_b = rng.permutation(og_b)
            sa = {m: set(perm_a[idx].tolist()) for m, idx in idxs_a.items()
                  if len(idx)}
            sb = {m: set(perm_b[idx].tolist()) for m, idx in idxs_b.items()
                  if len(idx)}
            for k, (ma, mb, shared, oc) in enumerate(cand):
                qa, qb = sa.get(ma), sb.get(mb)
                denom = min(len(qa), len(qb)) if qa and qb else 0
                if not denom:
                    continue
                nulloc = len(qa & qb) / denom
                null_med[k].append(nulloc)
                if nulloc >= oc:
                    null_ge[k] += 1

    # 2026-10-10 修复: 逐对经验 p (P_Emp) 无跨对多重检验校正 → 模块对多时假阳
    # 膨胀。对全部候选对的经验 p 做 BH-FDR, 输出 P_Emp_BH 与 Significant_BH
    # (P_Emp_BH < 0.05); P_Emp 原值保留以兼容下游 Wilcoxon。
    p_emp_all = np.array(
        [(1 + null_ge[k]) / (args.permutations + 1) for k in range(len(cand))],
        dtype=float)
    p_emp_bh = bh_fdr(p_emp_all)

    rows = []
    for k, (ma, mb, shared, oc) in enumerate(cand):
        p_emp = float(p_emp_all[k])
        p_bh = float(p_emp_bh[k])
        med = float(np.median(null_med[k])) if null_med[k] else float("nan")
        rows.append({
            "Name_A": args.name_a, "Module_A": ma, "N_OG_A": len(sets_a[ma]),
            "Name_B": args.name_b, "Module_B": mb, "N_OG_B": len(sets_b[mb]),
            "Shared_OG": shared, "OC": round(oc, 4),
            "Conserved": bool(oc >= args.oc_cutoff and shared >= args.min_shared_og),
            "P_Emp": round(p_emp, 4), "P_Emp_BH": round(p_bh, 4),
            "Significant_BH": bool(p_bh < 0.05),
            "Null_OC_Median": round(med, 4),
        })
    pairs = pd.DataFrame(rows).sort_values(
        ["Conserved", "OC"], ascending=[False, False]) if rows else pd.DataFrame(
        columns=["Name_A", "Module_A", "N_OG_A", "Name_B", "Module_B", "N_OG_B",
                 "Shared_OG", "OC", "Conserved", "P_Emp", "P_Emp_BH",
                 "Significant_BH", "Null_OC_Median"])

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(out / "conservation_pairs.tsv", sep="\t", index=False)
    n_cons = int(pairs["Conserved"].sum()) if len(pairs) else 0
    pd.DataFrame([{
        "name_a": args.name_a, "name_b": args.name_b,
        "n_modules_a": len(mods_a), "n_modules_b": len(mods_b),
        "n_pairs": len(stats), "n_candidate_pairs": len(cand),
        "oc_cutoff": args.oc_cutoff, "min_shared_og": args.min_shared_og,
        "permutations": args.permutations, "seed": args.seed,
        "n_conserved": n_cons, "unmapped_genes_a": unm_a,
        "unmapped_genes_b": unm_b,
    }]).to_csv(out / "conservation_summary.tsv", sep="\t", index=False)
    print(f"保守模块对 (OC≥{args.oc_cutoff:g} 且共享 OG≥{args.min_shared_og}): "
          f"{n_cons}")
    print(f"输出: {out / 'conservation_pairs.tsv'} / conservation_summary.tsv")


if __name__ == "__main__":
    main()


def main_testable(argv):
    """测试入口: 显式传 argv 列表, 不动 sys.argv。"""
    main(argv)
