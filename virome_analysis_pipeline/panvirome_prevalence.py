#!/usr/bin/env python3
"""
panvirome_prevalence.py — 病毒组序列级发生率汇总表生成器
========================================================
从 detect 阶段的 high_conf.summary.tsv 出发，合并 rescue 的属平均长度、
CheckV 完备性、以及 taxonomy 的完整分类等级，生成一张"序列级发生率总表"。

以 contig_id (vOTU 序列名) 为最小单元统计发生率，而不是按物种名聚合，
因为分类工具的物种标签可能不确定/不准。每一列都能回溯到原始文件。

输入（按宿主物种，路径自动拼接）:
  {base}/prevalence_{sp}_{detect_tag}/02_filtering/high_conf.summary.tsv
      检测阶段的高置信检出表 (auto_known_virus.py --stage detect --filter 输出)
  {base}/{sp_prefix}_{sp}_out/08_Rescue/Plant/rescue_report.tsv
      拯救报告 (提供 genus_avg_len / branch / method)
  {base}/{sp_prefix}_{sp}_out/07_Checkv/*/completeness.tsv
      CheckV 完备性 (aai_completeness / aai_confidence)
  {base}/{sp_prefix}_{sp}_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv
      完整分类等级 (Realm..Species)

输出:
  {out}  序列级发生率总表 (默认 prevalence_full_table.tsv)

用法示例:
  python panvirome_prevalence.py \
      --base ~/virus/data-2026/data-test \
      --species barbarum,ruthenicum,chinense,amarum \
      --detect-tag bowtie2 \
      --out prevalence_full_table.tsv
"""

import argparse
import os
import sys
import csv
import re
from collections import defaultdict
from pathlib import Path

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))


# ── 输出列 ──────────────────────────────────────────────
OUT_COLS = [
    "host_species", "contig_id", "length", "genus_avg_len",
    "n_samples", "prevalence_%", "checkv_completeness", "checkv_confidence",
    "avg_ANI", "Realm", "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="病毒组序列级发生率汇总表生成器",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--base", required=True,
                   help="项目根目录，含 {sp_prefix}_{sp}_out/ 与 prevalence_{sp}_{tag}/")
    p.add_argument("--species", default="barbarum,ruthenicum,chinense,amarum",
                   help="宿主物种列表，逗号分隔 (默认 4 种枸杞)")
    p.add_argument("--detect-tag", default="bowtie2",
                   help="检测输出目录后缀，即 prevalence_{sp}_{tag}/ (默认 bowtie2)")
    p.add_argument("--sp-prefix", default="RNA-Lycium",
                   help="discovery 输出目录前缀，即 {prefix}_{sp}_out/ (默认 RNA-Lycium)")
    p.add_argument("--out", default="prevalence_full_table.tsv",
                   help="输出 TSV 路径 (默认 prevalence_full_table.tsv)")
    p.add_argument("--min-prevalence", type=float, default=0.0,
                   help="最低发生率%%过滤，0 不过滤 (默认 0)")
    p.add_argument("--min-samples", type=int, default=0,
                   help="最低检出样本数过滤，0 不过滤 (默认 0)")
    p.add_argument("--sample-counts", default=None,
                   help="手动指定样本总数，格式 sp1:n1,sp2:n2 (默认从 detect 日志自动读)")
    p.add_argument("--checkv-hosts", default="Plant,Fungi,Animal,Unknown,Bacteria,Archaea,Algae,Protist,Mammalia",
                   help="CheckV 宿主目录列表，逗号分隔 (默认全部宿主)")
    return p.parse_args()


def clean(v):
    """清洗空值/占位符"""
    v = (v or "").strip()
    return "" if v.lower() in ("na", "nan", "none", "-", "") else v


def load_tsv_map(path, key_col):
    """读 TSV -> {key: row_dict}"""
    m = {}
    if not path.exists():
        return m
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            k = (row.get(key_col) or "").strip()
            if k:
                m[k] = row
    return m


def read_sample_total(detect_dir):
    """从 detect 日志读取样本总数 (总计发现样本: N 个)"""
    log = detect_dir / "detect.log"
    if not log.exists():
        return None
    txt = log.read_text(errors="ignore")
    m = re.search(r"总计发现样本[:\s]*(\d+)", txt)
    if m:
        return int(m.group(1))
    return None


def load_checkv(sp_out, hosts):
    """合并各宿主目录的 CheckV completeness.tsv -> {contig_id: (completeness, confidence)}"""
    m = {}
    cv_dir = sp_out / _D['d_checkv']
    if not cv_dir.exists():
        return m
    for host in hosts:
        f = cv_dir / host / "completeness.tsv"
        if not f.exists():
            continue
        with open(f, encoding="utf-8") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                cid = (row.get("contig_id") or "").strip()
                if cid:
                    m[cid] = (
                        (row.get("aai_completeness") or "").strip(),
                        (row.get("aai_confidence") or "").strip(),
                    )
    return m


def aggregate_high_conf(high_conf_path):
    """按 Rep_Accession 分组，返回 {acc: {length, n_samples, cov, depth, tpm, ani, label}}"""
    groups = defaultdict(lambda: {
        "n_samples": set(), "length": None, "cov": [], "depth": [],
        "tpm": [], "ani": [], "label": None,
    })
    with open(high_conf_path, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            acc = (row.get("Rep_Accession") or "").strip()
            if not acc:
                continue
            g = groups[acc]
            g["n_samples"].add(row.get("Sample", ""))
            if g["length"] is None:
                g["length"] = row.get("Rep_Length", "")
            if g["label"] is None:
                g["label"] = row.get("Adjusted_Species", "")
            for col, key in [("cov", "Rep_Coverage(%)"), ("depth", "Rep_MeanDepth"),
                             ("tpm", "Asm_TPM"), ("ani", "Avg_Read_ANI")]:
                v = row.get(key, "")
                if v not in ("", "NA", "nan"):
                    try:
                        g[col].append(float(v))
                    except ValueError:
                        pass
    return groups


def _mean(vals):
    return round(sum(vals) / len(vals), 2) if vals else ""


def main():
    args = parse_args()
    base = Path(args.base).expanduser()
    species = [s.strip() for s in args.species.split(",") if s.strip()]
    hosts = [h.strip() for h in args.checkv_hosts.split(",") if h.strip()]

    # 手动样本总数覆盖
    manual_counts = {}
    if args.sample_counts:
        for pair in args.sample_counts.split(","):
            if ":" in pair:
                k, v = pair.split(":", 1)
                manual_counts[k.strip()] = int(v.strip())

    all_rows = []
    for sp in species:
        detect_dir = base / f"prevalence_{sp}_{args.detect_tag}"
        sp_out = base / f"{args.sp_prefix}_{sp}_out"
        high_conf = detect_dir / _D['a_filter'] / "high_conf.summary.tsv"
        if not high_conf.exists():
            print(f"[SKIP] {sp}: high_conf 不存在 ({high_conf})")
            continue

        n_total = manual_counts.get(sp) or read_sample_total(detect_dir)
        if n_total is None:
            print(f"[WARN] {sp}: 无法从日志读取样本总数，请用 --sample-counts 指定")
            n_total = 0

        groups = aggregate_high_conf(high_conf)
        rr = load_tsv_map(sp_out / _D['d_rescue'] / "Plant" / "rescue_report.tsv", "contig_id")
        tx = load_tsv_map(sp_out / _D['d_taxonomy'] / "Votus.integrated" / "final_integrated_classification.tsv", "contig_id")
        cv = load_checkv(sp_out, hosts)

        n_rr = sum(1 for a in groups if a in rr)
        n_tx = sum(1 for a in groups if a in tx)
        n_cv = sum(1 for a in groups if a in cv)
        print(f"{sp}: {len(groups)} 序列 (总样本={n_total}) | rescue={n_rr} taxonomy={n_tx} checkv={n_cv}")

        for acc in sorted(groups, key=lambda a: -len(groups[a]["n_samples"])):
            g = groups[acc]
            n_samples = len(g["n_samples"])
            prev = round(n_samples / n_total * 100, 1) if n_total else 0.0
            if n_samples < args.min_samples or prev < args.min_prevalence:
                continue
            r = rr.get(acc, {})
            t = tx.get(acc, {})
            c = cv.get(acc, ("", ""))
            all_rows.append({
                "host_species": sp,
                "contig_id": acc,
                "length": g["length"] or "",
                "genus_avg_len": clean(r.get("genus_avg_len")),
                "n_samples": str(n_samples),
                "prevalence_%": f"{prev:.1f}",
                "checkv_completeness": c[0],
                "checkv_confidence": c[1],
                "avg_ANI": str(_mean(g["ani"])),
                "Realm": clean(t.get("Realm")),
                "Kingdom": clean(t.get("Kingdom")),
                "Phylum": clean(t.get("Phylum")),
                "Class": clean(t.get("Class")),
                "Order": clean(t.get("Order")),
                "Family": clean(t.get("Family")),
                "Genus": clean(t.get("Genus")),
                "Species": clean(t.get("Species")),
            })

    out = Path(args.out).expanduser()
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)
    print(f"\n完成: {out} ({len(all_rows)} 行, {len(species)} 物种)")


if __name__ == "__main__":
    main()
