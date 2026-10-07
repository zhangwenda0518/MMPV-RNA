#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""quant_qc.py — 宿主定量 QC (⑧ 宿主表达模块, 03_QC 阶段)
=========================================================================
从 salmon quant 输出 (aux_info/meta_info.json) 提取每样本:
  num_processed / num_mapped / percent_mapped, 并按双阈值判定 pass/fail。

默认阈值:
  --min-mapped    50.0   比对率 >= 50%      (16S/垃圾样本通常 0-5%;
                          同时会剔掉宿主转录组占比低的样本 — 宏转录组
                          队列常见, 下游样本量不足时可 CLI 放宽)
  --min-processed 1e6    处理 reads >= 100 万

输出:
  qc_summary.tsv          每样本: num_processed, num_mapped, percent_mapped, pass, note
  qc_passed_samples.txt   通过样本名单 (每行一个 sample 名)

用法:
  python quant_qc.py --quant_dir <08_HostExpression>/02_Quant \
      --output_dir <08_HostExpression>/03_QC
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List


def find_sample_dirs(quant_dir: Path) -> List[Path]:
    """定位 salmon 输出目录下的逐样本子目录 ({quant_dir}/quant/<sample>/)。"""
    root = quant_dir / "quant"
    if not root.is_dir():
        return []
    return sorted((p for p in root.iterdir() if p.is_dir()),
                  key=lambda p: p.name)


def read_sample_qc(sample_dir: Path) -> Dict:
    """读取单样本 meta_info.json; 缺失/损坏时返回 note 标记。"""
    meta = sample_dir / "aux_info" / "meta_info.json"
    if not meta.is_file():
        return {"sample": sample_dir.name, "num_processed": None,
                "num_mapped": None, "percent_mapped": None,
                "note": "no_meta_info"}
    try:
        info = json.loads(meta.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        return {"sample": sample_dir.name, "num_processed": None,
                "num_mapped": None, "percent_mapped": None,
                "note": f"meta_info_unreadable: {exc}"}
    return {
        "sample": sample_dir.name,
        "num_processed": info.get("num_processed"),
        "num_mapped": info.get("num_mapped"),
        "percent_mapped": info.get("percent_mapped"),
        "note": "",
    }


def collect_qc(quant_dir: Path,
               min_mapped: float = 50.0,
               min_processed: float = 1e6) -> List[Dict]:
    """扫描全部样本, 返回逐样本 QC 记录 (含 pass 判定)。"""
    rows: List[Dict] = []
    for sdir in find_sample_dirs(quant_dir):
        row = read_sample_qc(sdir)
        n_proc = row["num_processed"]
        pct = row["percent_mapped"]
        if n_proc is None or pct is None:
            row["pass"] = False
        else:
            ok_proc = n_proc >= min_processed
            ok_map = pct >= min_mapped
            row["pass"] = bool(ok_proc and ok_map)
            if not ok_proc and not ok_map:
                row["note"] = "low_reads+low_mapped"
            elif not ok_proc:
                row["note"] = "low_reads"
            elif not ok_map:
                row["note"] = "low_mapped"
        rows.append(row)
    return rows


def write_qc(rows: List[Dict], output_dir: Path) -> Dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = output_dir / "qc_summary.tsv"
    passed = output_dir / "qc_passed_samples.txt"

    with summary.open("w", encoding="utf-8", newline="") as f:
        f.write("sample\tnum_processed\tnum_mapped\tpercent_mapped\tpass\tnote\n")
        for r in rows:
            f.write("\t".join([
                r["sample"],
                "" if r["num_processed"] is None else str(r["num_processed"]),
                "" if r["num_mapped"] is None else str(r["num_mapped"]),
                "" if r["percent_mapped"] is None else f"{r['percent_mapped']:.2f}",
                "1" if r["pass"] else "0",
                r["note"],
            ]) + "\n")

    with passed.open("w", encoding="utf-8") as f:
        for r in rows:
            if r["pass"]:
                f.write(r["sample"] + "\n")

    return {"summary": summary, "passed": passed}


def main() -> None:
    ap = argparse.ArgumentParser(description="宿主定量 QC (salmon meta_info)")
    ap.add_argument("--quant_dir", required=True,
                    help="salmon 定量输出目录 (02_Quant, 内含 quant/<sample>/)")
    ap.add_argument("--output_dir", required=True, help="输出目录 (03_QC)")
    ap.add_argument("--min-mapped", type=float, default=50.0,
                    help="percent_mapped 阈值 (%%), 默认 50")
    ap.add_argument("--min-processed", type=float, default=1e6,
                    help="处理 reads 数阈值, 默认 1e6")
    args = ap.parse_args()

    quant_dir = Path(args.quant_dir)
    rows = collect_qc(quant_dir, args.min_mapped, args.min_processed)
    if not rows:
        print(f"[ERROR] 未在 {quant_dir}/quant/ 找到任何样本输出; "
              f"请先跑 quant 阶段", file=sys.stderr)
        sys.exit(1)

    n_pass = sum(1 for r in rows if r["pass"])
    out = write_qc(rows, Path(args.output_dir))
    print(f"QC: {n_pass}/{len(rows)} 通过 (min_mapped={args.min_mapped}%, "
          f"min_processed={args.min_processed:g})")
    for r in rows:
        if not r["pass"]:
            pct = "-" if r["percent_mapped"] is None else f"{r['percent_mapped']:.1f}%"
            print(f"  [剔除] {r['sample']}  mapped={pct}  {r['note']}")
    print(f"输出: {out['summary']}")
    print(f"      {out['passed']}")


if __name__ == "__main__":
    main()
