#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""derive_virus_labels.py — 病毒标签推导 (⑧ 宿主表达模块, 04_Labels 阶段)
=========================================================================
把 ④ 已知病毒分析管线的样本×病毒检出表 (all_viruses.best.summary.tsv)
转成宿主表达分析所需的分组变量与剂量变量:

  - 逐病毒标签: 样本级阳性判定 (基因组覆盖度 >= 阈值, 默认 90%, 与
    goji-virome 口径一致) + 连续载量 (load 列, 供载量×表达相关扫描)
  - 样本级标签: Virus_any (+/-) 与逐病毒 0/1 阳性列 + 连续载量列
    (供 DESeq2 分组、逐病毒对比与载量扫描直接消费)
  - 元数据合并 (可选): Core14 表 Tissue / ScientificName / CollectionDate

输入列名自动解析 (可用 --*-col 覆盖):
  Sample  Adjusted_Species  Rep_Coverage(%)  Rep_MeanDepth  Asm_TPM  Rep_Reads

输出 (默认落在 04_Labels/):
  virus_labels.tsv    长表: Sample, Virus, Covered_pct, MeanDepth, Load, Reads, Positive
  sample_labels.tsv   宽表: Sample, Virus_any, N_virus, <病毒>_positive/_load, [Tissue, ...]
  virus_columns.tsv   列名映射与逐病毒计数: Virus, Column, N_positive, N_samples

用法:
  python derive_virus_labels.py \
      --detect <项目>/04_Analysis/01_Detection/summary/all_viruses.best.summary.tsv \
      --samples <08_HostExpression>/03_QC/qc_summary.tsv \
      --metadata <项目>/01_PublicData/02_Unified/Global_Unified_Metadata_Core14.csv \
      --output_dir <08_HostExpression>/04_Labels
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

# 逻辑列 → 候选实际列名 (按优先级; 覆盖 MMPV 检出表与常见变体)
COL_CANDIDATES: Dict[str, Sequence[str]] = {
    "sample": ("Sample", "sample", "Run", "run"),
    "virus":  ("Adjusted_Species", "Species_NCBI", "Species_ICTV",
               "Virus", "virus", "Species"),
    "cov":    ("Rep_Coverage(%)", "Rep_Coverage", "Coverage(%)",
               "Coverage", "cov", "Covered_pct"),
    "depth":  ("Rep_MeanDepth", "MeanDepth", "Depth", "depth"),
    "load":   ("Asm_TPM", "TPM", "Load", "load"),
    "reads":  ("Rep_Reads", "Asm_EM_Reads", "Reads", "Uniq_Reads"),
}


def sniff_read(path: Path) -> pd.DataFrame:
    """读取 TSV/CSV (按首行分隔符嗅探); 检出表统一按字符串读入。"""
    first = path.open("r", encoding="utf-8", errors="replace").readline()
    sep = "," if first.count(",") > first.count("\t") else "\t"
    return pd.read_csv(path, sep=sep, dtype=str, low_memory=False)


def read_detection(paths: Sequence[Path]) -> pd.DataFrame:
    """读取一个或多个检出表 (目录自动展开为其中的 *.tsv/*.csv) 并纵向合并。"""
    files: List[Path] = []
    for p in paths:
        if p.is_dir():
            files.extend(sorted(list(p.glob("*.tsv")) + list(p.glob("*.csv"))))
        else:
            files.append(p)
    if not files:
        raise FileNotFoundError(f"检出表为空: {[str(p) for p in paths]}")
    frames = []
    for f in files:
        df = sniff_read(f)
        df["__source_table"] = f.name
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def resolve_columns(df: pd.DataFrame,
                    overrides: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """把逻辑列解析为实际列名; 必需列缺失即报错 (信息含实际表头)。"""
    overrides = {k: v for k, v in (overrides or {}).items() if v}
    cols: Dict[str, str] = {}
    for logical, candidates in COL_CANDIDATES.items():
        if logical in overrides:
            if overrides[logical] not in df.columns:
                raise KeyError(f"指定列 {overrides[logical]!r} 不在表中; "
                               f"实际表头: {list(df.columns)}")
            cols[logical] = overrides[logical]
            continue
        hit = next((c for c in candidates if c in df.columns), None)
        if hit is not None:
            cols[logical] = hit
    for required in ("sample", "virus", "cov"):
        if required not in cols:
            raise KeyError(
                f"检出表缺少必需列 {required!r} (候选: "
                f"{COL_CANDIDATES[required]}); 实际表头: {list(df.columns)}; "
                f"可用 --{required}-col 显式指定")
    return cols


def _to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def build_virus_labels(df: pd.DataFrame, cols: Dict[str, str],
                       cov_threshold: float = 90.0,
                       depth_threshold: float = 0.0) -> pd.DataFrame:
    """逐 (Sample, Virus) 聚合 → 长表标签。

    同一 (Sample, Virus) 多行 (分节段/多参考) 时: 覆盖与深度取 max,
    载量取 max, reads 求和。阳性判定: Covered_pct >= cov_threshold
    (可选叠加 MeanDepth >= depth_threshold, 默认关闭)。
    """
    work = pd.DataFrame({
        "Sample": df[cols["sample"]].astype(str).str.strip(),
        "Virus": df[cols["virus"]].astype(str).str.strip(),
        "Covered_pct": _to_num(df[cols["cov"]]),
        "MeanDepth": _to_num(df[cols["depth"]]) if "depth" in cols else 0.0,
        "Load": _to_num(df[cols["load"]]) if "load" in cols else 0.0,
        "Reads": _to_num(df[cols["reads"]]) if "reads" in cols else 0.0,
    })
    work = work[(work["Sample"] != "") & (work["Virus"] != "")
                & (work["Virus"].str.lower() != "nan")]

    agg = (work.groupby(["Sample", "Virus"], as_index=False)
                .agg(Covered_pct=("Covered_pct", "max"),
                     MeanDepth=("MeanDepth", "max"),
                     Load=("Load", "max"),
                     Reads=("Reads", "sum")))
    positive = agg["Covered_pct"].fillna(0) >= cov_threshold
    if depth_threshold > 0:
        positive &= agg["MeanDepth"].fillna(0) >= depth_threshold
    agg["Positive"] = positive.astype(int)
    return agg


def sanitize_virus(name: str) -> str:
    """病毒名 → 列名安全形式 (空格/引号/斜杠等一律转下划线)。"""
    s = re.sub(r"[^0-9A-Za-z_.-]+", "_", str(name)).strip("_")
    return s or "virus"


def build_sample_labels(vlabels: pd.DataFrame,
                        samples: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """长表 → 样本级宽表 (Virus_any +/- , N_virus, 逐病毒阳性列与载量列)。

    samples: 全样本名单 (含未检出的阴性样本); 未提供的样本视为全阴性。
    """
    all_samples = sorted(set(samples) if samples else set(vlabels["Sample"]))
    viruses = sorted(vlabels["Virus"].unique())

    rows: List[Dict] = []
    by_sample = {s: g for s, g in vlabels.groupby("Sample")}
    for s in all_samples:
        g = by_sample.get(s)
        rec: Dict = {"Sample": s}
        n_pos = 0
        for v in viruses:
            col = sanitize_virus(v)
            if g is None:
                pos, load = 0, 0.0
            else:
                sub = g[g["Virus"] == v]
                if sub.empty:
                    pos, load = 0, 0.0
                else:
                    pos = int(sub["Positive"].iloc[0])
                    load = float(sub["Load"].iloc[0]) if pd.notna(sub["Load"].iloc[0]) else 0.0
            n_pos += pos
            rec[f"{col}_positive"] = pos
            rec[f"{col}_load"] = load
        rec["Virus_any"] = "+" if n_pos > 0 else "-"
        rec["N_virus"] = n_pos
        rows.append(rec)

    base = ["Sample", "Virus_any", "N_virus"]
    per_virus = [c for r in rows for c in r if c not in base]
    ordered = base + sorted(set(per_virus))
    out = pd.DataFrame(rows)[ordered]
    return out


def read_sample_list(path: Path) -> List[str]:
    """读样本名单: 支持每行一个, 或带表头的 TSV (取 sample/Sample 列或首列)。"""
    with path.open(encoding="utf-8") as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return []
    first = lines[0].split("\t")
    names: List[str] = []
    if any("sample" in c.lower() for c in first):
        idx = next(i for i, c in enumerate(first) if "sample" in c.lower())
        for l in lines[1:]:
            parts = l.split("\t")
            if idx < len(parts) and parts[idx].strip():
                names.append(parts[idx].strip())
    else:
        # 无表头: 整行或首列
        for l in lines:
            names.append(l.split("\t")[0].strip())
    return names


def merge_metadata(slabels: pd.DataFrame, meta_path: Optional[Path]) -> pd.DataFrame:
    """把 Core14 元数据的 Tissue/ScientificName/CollectionDate 并到样本表。"""
    if meta_path is None:
        return slabels
    meta = sniff_read(Path(meta_path))
    key = next((c for c in ("Run", "Sample", "sample", "run") if c in meta.columns),
               None)
    if key is None:
        print(f"[WARN] 元数据表无 Run/Sample 列, 跳过合并: {meta_path}",
              file=sys.stderr)
        return slabels
    keep = [key] + [c for c in ("Tissue", "ScientificName", "CollectionDate")
                    if c in meta.columns]
    meta = meta[keep].rename(columns={key: "Sample"}).drop_duplicates("Sample")
    merged = slabels.merge(meta, on="Sample", how="left")
    n_miss = merged["Tissue"].isna().sum() if "Tissue" in merged.columns else 0
    if n_miss:
        print(f"[INFO] {n_miss}/{len(merged)} 样本在元数据表中无 Tissue 标注",
              file=sys.stderr)
    return merged


def main() -> None:
    ap = argparse.ArgumentParser(description="由 ④ 检出表推导病毒标签")
    ap.add_argument("--detect", required=True, nargs="+",
                    help="检出表路径 (all_viruses.best.summary.tsv) 或所在目录; 可多个")
    ap.add_argument("--output_dir", required=True, help="输出目录 (04_Labels)")
    ap.add_argument("--samples", help="全样本名单 (如 03_QC/qc_summary.tsv); "
                                      "用于把未检出样本补为阴性")
    ap.add_argument("--metadata", help="Core14 元数据表 (csv/tsv), 合并 Tissue 等")
    ap.add_argument("--cov-threshold", type=float, default=90.0,
                    help="阳性判定: 基因组覆盖度阈值 (%%), 默认 90")
    ap.add_argument("--depth-threshold", type=float, default=0.0,
                    help="阳性判定附加深度阈值; 0=关闭 (默认)")
    for logical in ("sample", "virus", "cov", "depth", "load", "reads"):
        ap.add_argument(f"--{logical}-col", help=f"显式指定 {logical} 列名 (默认自动解析)")
    args = ap.parse_args()

    df = read_detection([Path(p) for p in args.detect])
    overrides = {k: getattr(args, f"{k}_col")
                 for k in ("sample", "virus", "cov", "depth", "load", "reads")}
    cols = resolve_columns(df, overrides)
    print(f"检出表: {len(df)} 行 × {len(df.columns)} 列  "
          f"(样本列={cols['sample']!r}, 病毒列={cols['virus']!r}, "
          f"覆盖列={cols['cov']!r})")

    vlabels = build_virus_labels(df, cols, args.cov_threshold, args.depth_threshold)
    n_virus = vlabels["Virus"].nunique()
    n_detected = vlabels["Sample"].nunique()
    print(f"逐病毒标签: {n_detected} 样本 × {n_virus} 病毒 "
          f"(阳性判定 cov>={args.cov_threshold}%)")

    samples = read_sample_list(Path(args.samples)) if args.samples else None
    if samples:
        extra = sorted(set(samples) - set(vlabels["Sample"]))
        if extra:
            print(f"补充 {len(extra)} 个未检出样本为全阴性 (来自 --samples)")
    slabels = build_sample_labels(vlabels, samples)
    slabels = merge_metadata(slabels, Path(args.metadata) if args.metadata else None)

    n_pos = int((slabels["Virus_any"] == "+").sum())
    print(f"样本分组: 阳性 {n_pos} / 阴性 {len(slabels) - n_pos}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    vlabels.to_csv(out_dir / "virus_labels.tsv", sep="\t", index=False)
    slabels.to_csv(out_dir / "sample_labels.tsv", sep="\t", index=False)

    stats = (vlabels.groupby("Virus", as_index=False)
                    .agg(N_positive=("Positive", "sum"),
                         N_samples=("Sample", "nunique")))
    stats["Column"] = stats["Virus"].map(sanitize_virus)
    stats.to_csv(out_dir / "virus_columns.tsv", sep="\t", index=False)

    print(f"输出: {out_dir / 'virus_labels.tsv'}   ({len(vlabels)} 行)")
    print(f"      {out_dir / 'sample_labels.tsv'}  ({len(slabels)} 样本)")
    print(f"      {out_dir / 'virus_columns.tsv'}")


if __name__ == "__main__":
    main()
