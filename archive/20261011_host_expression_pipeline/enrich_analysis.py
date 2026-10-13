#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""enrich_analysis.py — 功能富集分析 (⑧ 宿主表达模块, 07_Enrich 阶段)
=========================================================================
对 ⑧ 产出的基因集 (06_DE 显著 DEG / 05_Coexp Δz hub / 模块成员) 做
超几何过表达分析 (ORA) + BH-FDR。本脚本**不跑注释**——注释由用户事先
对宿主蛋白/CDS 跑一次 eggNOG-mapper (或 Mercator4 → 导出映射表), 本脚本
消费其标准输出, 因此物种无关、零重依赖 (scipy 超几何)。

注释输入 (二选一, 自动嗅探):
  1. eggNOG-mapper v2 `.emapper.annotations`: 提取 GOs 列 (GO:XXXX,…) 与
     KEGG_ko 列 (ko:KXXXX) 两类 term; namespace 分别记 GO/KEGG_KO
  2. 通用映射表: TSV 必含 gene,term 列, 可选 namespace,name 列
     (MapMan bin / Mercator4 导出 / 任意 自定义词汇均可)

基因集输入: --gene-set NAME=path (可重复); path 为每行一个基因的 txt,
或首列为基因名的 TSV。背景: --background (每行一个基因; 统计上必需,
编排器自动传 TPM 全基因)。

输出 (--output_dir, 07_Enrich/):
  enrich_<Set>.tsv    每基因集的全部检验 term (按 Padj 排序, Significant 列
                      标记 Padj<0.05; 集合为空/跳过时仅表头)
  enrich_summary.tsv  每基因集: 规模/可注释数/检验 term 数/显著数

统计: P(X>=k) = hypergeom.sf(k-1, N_BG, N_Term_BG, N_Set)。N_BG = 背景∩注释,
N_Set = 查询∩注释 (可注释空间; 2026-10-03 修复 —— 无注释基因不可能命中任何
term, 不该进入抽样空间)。该 N 口径与 R clusterProfiler::enricher 一致
(enricher 内部对传入 universe 做 `intersect(universe, TERM2GENE 基因)`),
故跨工具 p 值可比 (2026-10-10 明确)。term 背景规模过滤在
[--min-term-genes, --max-term-genes]; 检验族 = **背景内全部通过规模门槛的
term** (2026-10-10 修复: 旧实现只对"被查询命中的 term"做 BH, 检验族由查询
结果决定 → 假阳膨胀); BH 校正全部检验 term (k=0 的 term p=1 仍计入族)。
输出表同时给出 N_BG(可注释) 与 N_BG_Raw(背景原始规模) 便于跨工具核对口径。

用法:
  python enrich_analysis.py --annotation emapper.annotations \
      --gene-set DE=sets/DE.txt --gene-set module_0=sets/module_0.txt \
      --background background_genes.txt --output_dir 07_Enrich
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


# ── 注释解析 ──────────────────────────────────────────────────────────

def _looks_like_emapper(path: Path) -> bool:
    """头部若干行内存在 '#…query…' 表头行即为 emapper 格式
    (真文件前两行是 ##/### 版本注释, #query 在第 3 行)。"""
    with path.open(encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i >= 50:
                break
            s = line.strip()
            if s.startswith("#") and "query" in s.split("\t")[0].lower():
                return True
    return False


def parse_emapper(path: Path) -> pd.DataFrame:
    """eggNOG-mapper v2 → 长表 (gene, term, namespace, name='')。

    header 动态定位 (文件头有 2-3 行 # 注释), 取 GOs 与 KEGG_ko 两列。"""
    with path.open(encoding="utf-8", errors="replace") as fh:
        header_idx = None
        for i, line in enumerate(fh):
            if line.startswith("#") and "query" in line.split("\t")[0].lower():
                header_idx = i
                break
        if header_idx is None:
            sys.exit(f"[ERROR] emapper 注释找不到 #query 表头行: {path}")
    df = pd.read_csv(path, sep="\t", skiprows=header_idx, dtype=str)
    df.columns = [c.lstrip("#") for c in df.columns]
    qcol = "query" if "query" in df.columns else df.columns[0]
    rows = []
    for src, col, ns in (("GO", "GOs", "GO"), ("KEGG_KO", "KEGG_ko", "KEGG_KO")):
        if col not in df.columns:
            print(f"[WARN] emapper 无 {col} 列, 跳过 {src} 注释")
            continue
        sub = df[[qcol, col]].dropna()
        for gene, cell in zip(sub[qcol], sub[col]):
            for t in str(cell).replace(",", " ").split():
                t = t.strip()
                # 2026-10-10 修复: eggNOG 未注释用 "-" 占位, 旧实现把 "-" 当
                # 合法 term 建背景集合 (H06), 现跳过空串与 "-"。
                if not t or t == "-":
                    continue
                rows.append((str(gene), t, ns))
    return pd.DataFrame(rows, columns=["gene", "term", "namespace"])


def parse_generic_map(path: Path) -> pd.DataFrame:
    """通用 gene↔term 映射表 → 长表。"""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if not {"gene", "term"}.issubset(df.columns):
        sys.exit(f"[ERROR] 通用注释表需含 gene/term 列 (可选 namespace/name): {path}")
    df = df.rename(columns={c: c.strip().lower() for c in df.columns})
    out = df[["gene", "term"]].copy()
    out["namespace"] = df.get("namespace", "")
    out["name"] = df.get("name", "") if "name" in df.columns else ""
    out = out.dropna(subset=["gene", "term"])
    # 2026-10-10 (H06): 剔除 "-" / 空串占位符, 避免其被当合法 term 建背景集合。
    out["gene"] = out["gene"].astype(str).str.strip()
    out["term"] = out["term"].astype(str).str.strip()
    out = out[(out["gene"] != "") & (out["gene"] != "-")
              & (out["term"] != "") & (out["term"] != "-")]
    return out


def parse_annotations(path: Path) -> pd.DataFrame:
    """返回长表 gene,term,namespace[,name]; emapper 自动嗅探。"""
    path = Path(path)
    if not path.is_file():
        sys.exit(f"[ERROR] 注释文件不存在: {path}")
    df = parse_emapper(path) if _looks_like_emapper(path) \
        else parse_generic_map(path)
    df["gene"] = df["gene"].astype(str)
    df["term"] = df["term"].astype(str)
    if "name" not in df.columns:
        df["name"] = ""
    return df


# ── 基因集 / 背景 ─────────────────────────────────────────────────────

def read_gene_list(path: Path) -> list:
    """txt (每行一基因, 无表头) 或 TSV (首列基因, 有表头) → 去重保序列表。

    判别: 首行含 tab → TSV (首行是表头, 取首列数据行); 否则逐行即基因。"""
    p = Path(path)
    if not p.is_file():
        sys.exit(f"[ERROR] 基因集文件不存在: {p}")
    lines = [l.rstrip("\r\n") for l in p.open(encoding="utf-8") if l.strip()]
    if not lines:
        return []
    if "\t" in lines[0]:
        lines = lines[1:]                      # TSV: 跳过表头, 取首列
        lines = [l.split("\t", 1)[0] for l in lines]
    seen, out = set(), []
    for g in lines:
        g = g.strip()
        if g and g not in seen:
            seen.add(g)
            out.append(g)
    return out


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
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


# ── ORA 主流程 ────────────────────────────────────────────────────────

def run_ora(gene_sets: dict, annot: pd.DataFrame, background: list,
            min_term: int, max_term: int):
    """返回 (逐 term 长表, 汇总行列表)。基因集 <3 个可注释基因时记 skip。

    2026-10-03 修复 (统计空间): 超几何的 N (总体) 与 n (抽样) 必须与"成功项
    空间"一致。旧实现 N=len(bg_set) / n=len(查询∩背景) 把无注释基因也算进
    抽样空间 (它们不可能命中任何 term), 当查询集注释率高于背景时 P 值系统性
    偏小。现统一到可注释空间: universe = 背景∩注释, query = 查询∩注释,
    M=term 背景规模、k 不变; 输出列 N_BG 也随之 = len(universe)。
    """
    bg_set = set(background)
    annot_genes = set(annot["gene"])
    universe = bg_set & annot_genes
    term2genes = {}
    for g, t in zip(annot["gene"], annot["term"]):
        if g in bg_set:
            term2genes.setdefault(t, set()).add(g)
    term_meta = (annot.drop_duplicates("term")
                 .set_index("term")[["namespace", "name"]].to_dict("index"))
    bg_terms = {t: gs for t, gs in term2genes.items()
                if min_term <= len(gs) <= max_term}
    term_bg_size = {t: len(gs) for t, gs in term2genes.items()}

    rows, summary = [], []
    for name, genes in gene_sets.items():
        in_bg = [g for g in genes if g in bg_set]
        gene_set = set(in_bg) & annot_genes      # 可注释空间内的查询集
        n_annot = len(gene_set)
        if n_annot < 3:
            summary.append({"Set": name, "N_Genes": len(genes),
                            "N_In_BG": len(in_bg), "N_Annotated": n_annot,
                            "N_Terms_Tested": 0, "N_Sig_Padj05": 0,
                            "Note": "可注释基因 <3, 跳过"})
            continue
        n_set = len(gene_set)
        # 2026-10-10 修复 (H05): 检验族 = 背景内全部通过规模门槛的 term, 而不是
        # "被查询命中的 term" —— 否则检验族由查询结果决定, BH 假阳膨胀。
        # 未命中的 term k=0 → p=1, 仍计入族。与 clusterProfiler 一致 (全 term 检验)。
        tested = list(bg_terms)
        n_bg = len(universe)
        pvals, hits = [], []
        for t in tested:
            k = len(term2genes[t] & gene_set)
            p = float(hypergeom.sf(k - 1, n_bg, term_bg_size[t], n_set))
            pvals.append(p)
            hits.append((t, k))
        padj = bh_fdr(np.asarray(pvals)) if pvals else np.array([])
        n_sig = 0
        for (t, k), pv, pj in zip(hits, pvals, padj):
            if pv != pv:
                continue
            sig = pj < 0.05
            n_sig += int(sig)
            meta = term_meta.get(t, {})
            fold = (k / n_set) / (term_bg_size[t] / n_bg) if n_set else float("nan")
            rows.append({
                "Set": name, "Term": t,
                "Namespace": meta.get("namespace", ""),
                "Name": meta.get("name", ""),
                "N_Set": n_set, "N_Set_in_Term": k,
                "N_Term_BG": term_bg_size[t], "N_BG": n_bg,
                "N_BG_Raw": len(bg_set),
                "Fold": round(fold, 3) if fold == fold else "",
                "Pvalue": pv, "Padj": pj, "Significant": bool(sig),
            })
        summary.append({"Set": name, "N_Genes": len(genes),
                        "N_In_BG": len(in_bg), "N_Annotated": n_annot,
                        "N_Terms_Tested": len(tested), "N_Sig_Padj05": n_sig,
                        "N_BG": n_bg, "N_BG_Raw": len(bg_set),
                        "Note": ""})
    table = pd.DataFrame(rows, columns=[
        "Set", "Term", "Namespace", "Name", "N_Set", "N_Set_in_Term",
        "N_Term_BG", "N_BG", "N_BG_Raw", "Fold", "Pvalue", "Padj",
        "Significant"])
    return table, summary


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description="功能富集 ORA (超几何+BH; 消费 eggNOG-mapper/通用注释)")
    ap.add_argument("--annotation", required=True,
                    help="注释: eggNOG-mapper .emapper.annotations 或 "
                         "通用 gene/term 映射 TSV")
    ap.add_argument("--gene-set", action="append", required=True,
                    help="NAME=path 基因集 (可重复; txt 每行一基因或 TSV 首列)")
    ap.add_argument("--background", required=True,
                    help="背景基因 (每行一个; 通常 = TPM 矩阵全部基因)")
    ap.add_argument("--output_dir", required=True, help="输出目录 (07_Enrich)")
    ap.add_argument("--min-term-genes", type=int, default=5,
                    help="term 背景规模下限 (默认 5)")
    ap.add_argument("--max-term-genes", type=int, default=500,
                    help="term 背景规模上限 (默认 500)")
    args = ap.parse_args(argv)

    annot = parse_annotations(Path(args.annotation))
    background = [l.strip() for l in Path(args.background)
                  .open(encoding="utf-8") if l.strip()]
    gene_sets = {}
    for spec in args.gene_set:
        name, _, p = spec.partition("=")
        if not p:
            sys.exit(f"[ERROR] --gene-set 需要 NAME=path 形式, 得到: {spec}")
        gene_sets[name] = read_gene_list(Path(p))
    print(f"注释: {annot['gene'].nunique()} 基因 / {annot['term'].nunique()} term; "
          f"背景 {len(background)}; 基因集 {len(gene_sets)}: "
          f"{ {k: len(v) for k, v in gene_sets.items()} }")

    table, summary = run_ora(gene_sets, annot, background,
                             args.min_term_genes, args.max_term_genes)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for r in summary:                      # 每个集合都落文件 (含空/跳过 → 仅表头)
        sub = table[table["Set"] == r["Set"]].sort_values("Padj") \
            if len(table) else table
        sub.to_csv(out / f"enrich_{r['Set']}.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(out / "enrich_summary.tsv", sep="\t",
                                 index=False)
    n_tot = sum(r["N_Sig_Padj05"] for r in summary)
    for r in summary:
        print(f"  {r['Set']}: 基因 {r['N_Genes']} (可注释 {r['N_Annotated']}), "
              f"检验 {r['N_Terms_Tested']} term, 显著 {r['N_Sig_Padj05']} "
              f"{('[' + r['Note'] + ']') if r['Note'] else ''}")
    print(f"输出: {out / 'enrich_summary.tsv'} "
          f"(显著 term 总数 {n_tot}; 逐集见 enrich_<Set>.tsv)")


if __name__ == "__main__":
    main()


def main_testable(argv):
    """测试入口: 显式传 argv 列表, 不动 sys.argv。"""
    main(argv)
