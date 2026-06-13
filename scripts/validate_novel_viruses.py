#!/usr/bin/env python3
"""
validate_novel_viruses.py — 基于分类层级判断新病毒 v2.0
========================================================

读取 taxonomy + CD-HIT 已知关联, 按分类完整性分级, 无需 BLASTN。

分类规则:
  Species ≠ NA              → ★ known (已知病毒)
  Genus ≠ NA, Species = NA  → ★★ novel species (新种)
  Family ≠ NA, Genus = NA   → ★★ novel genus   (新属)
  Order ≠ NA, Family = NA   → ★★★ novel family  (新科)
  Class ≠ NA, Order = NA    → ★★★ novel order   (新目)
  全是 NA                    → ★★★ truly novel   (全新)

CD-HIT 已知标记的 centroids 直接归为 ★ known。

用法:
  python validate_novel_viruses.py \\
      -i centoids.fasta \\
      --taxonomy final_integrated_classification.tsv \\
      --cdhit-known known_association.tsv \\
      --host ensemble_host_summary.tsv \\
      -o 07_Validation/
"""

import sys, os, logging, argparse, json, time
from pathlib import Path
from collections import defaultdict

import polars as pl
from Bio import SeqIO

TAX_LEVELS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def setup_logger(level="INFO"):
    l = logging.getLogger("validate")
    l.setLevel(getattr(logging, level.upper(), logging.INFO))
    if not l.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%H:%M:%S'))
        l.addHandler(h)
    return l


def classify_by_taxonomy(row):
    """
    根据分类完整性分级:
      Species ≠ NA → known
      Genus ≠ NA, Species = NA → novel_species
      Family ≠ NA, Genus = NA → novel_genus
      Class ≠ NA, Family = NA → novel_family
      全是 NA → truly_novel
    """
    # CD-HIT 已标记
    if row.get("cdhit_known", False):
        return "★ known", "known", "CD-HIT"

    species = row.get("Species", None)
    genus = row.get("Genus", None)
    family = row.get("Family", None)
    order = row.get("Order", None)
    class_ = row.get("Class", None)

    def has_val(v):
        return v and v != "NA" and str(v).strip() != ""

    if has_val(species):
        return "★ known", "known", "Species"
    if has_val(genus):
        return "★★ novel_species", "novel_species", "Genus"
    if has_val(family):
        return "★★ novel_genus", "novel_genus", "Family"
    if has_val(order):
        return "★★★ novel_family", "novel_family", "Order"
    if has_val(class_):
        return "★★★ novel_order", "novel_order", "Class"
    return "★★★ truly_novel", "truly_novel", "No_hit"


def generate_html_report(stats, taxonomy_df, out_path):
    """生成含 Plotly.js 的 HTML 报告"""
    counts = stats["counts"]
    total = sum(counts.values())
    class_label = stats.get("class_label", "Class")

    # 宿主分布
    host_counts = stats.get("host_counts", {})

    # 分类层级分布
    level_counts = stats.get("level_counts", {})

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<title>病毒新颖性验证报告</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #f5f7fa; color: #333; padding: 20px; }}
h1 {{ text-align: center; color: #1A5276; margin: 20px 0 10px; font-size: 28px; }}
.subtitle {{ text-align: center; color: #7F8C8D; margin-bottom: 30px; font-size: 14px; }}
.cards {{ display: flex; gap: 20px; justify-content: center; flex-wrap: wrap; margin-bottom: 30px; }}
.card {{ background: #fff; border-radius: 12px; padding: 24px 32px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); text-align: center; min-width: 160px; }}
.card .num {{ font-size: 42px; font-weight: 700; }}
.card .label {{ font-size: 14px; color: #7F8C8D; margin-top: 6px; }}
.card.known .num {{ color: #27AE60; }}
.card.novel .num {{ color: #F39C12; }}
.card.truly .num {{ color: #E74C3C; }}
.chart-row {{ display: flex; gap: 20px; justify-content: center; flex-wrap: wrap; margin-bottom: 30px; }}
.chart-box {{ background: #fff; border-radius: 12px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); padding: 20px; flex: 1; min-width: 400px; max-width: 600px; }}
table {{ border-collapse: collapse; width: 100%; }}
table th, table td {{ padding: 10px 14px; border-bottom: 1px solid #eee; text-align: left; font-size: 13px; }}
table th {{ background: #f0f4f8; color: #1A5276; font-weight: 600; }}
.footer {{ text-align: center; color: #95A5A6; font-size: 12px; margin-top: 40px; padding: 20px; }}
</style>
</head>
<body>

<h1>病毒新颖性验证报告</h1>
<p class="subtitle">validate_novel_viruses.py v2.0 — 基于分类层级判断 | {time.strftime('%Y-%m-%d %H:%M:%S')}</p>

<div class="cards">
  <div class="card known">
    <div class="num">{counts.get("known", 0):,}</div>
    <div class="label">★ 已知病毒 (Species ≠ NA)</div>
  </div>
  <div class="card novel">
    <div class="num">{counts.get("novel_species", 0) + counts.get("novel_genus", 0):,}</div>
    <div class="label">★★ 新种/新属</div>
  </div>
  <div class="card truly">
    <div class="num">{counts.get("novel_family", 0) + counts.get("novel_order", 0) + counts.get("truly_novel", 0):,}</div>
    <div class="label">★★★ 新科以上</div>
  </div>
  <div class="card">
    <div class="num" style="color:#2980B9">{total:,}</div>
    <div class="label">总计 vOTU</div>
  </div>
  <div class="card">
    <div class="num" style="color:#1ABC9C">{stats.get("total_samples", 0):,}</div>
    <div class="label">样本数</div>
  </div>
</div>

<div class="chart-row">
  <div class="chart-box" id="pieChart"></div>
  <div class="chart-box" id="levelChart"></div>
</div>

<div class="chart-row">
  <div class="chart-box" id="hostChart"></div>
  <div class="chart-box" id="freqChart"></div>
</div>

<script>
Plotly.newPlot("pieChart", [{{
  type: "pie",
  labels: ["★ 已知 (Species)", "★★ 新种 (Genus)", "★★ 新属 (Family)", "★★★ 新科+(Order/Class)", "★★★ 全新 (无分类)"],
  values: [{counts.get("known", 0)}, {counts.get("novel_species", 0)}, {counts.get("novel_genus", 0)}, {counts.get("novel_family", 0) + counts.get("novel_order", 0)}, {counts.get("truly_novel", 0)}],
  marker: {{ colors: ["#27AE60", "#2ECC71", "#F39C12", "#E67E22", "#E74C3C"] }},
  hole: 0.4,
  textinfo: "label+percent",
  textfont: {{ size: 12 }}
}}, {{ title: "vOTU 分类层级分布", height: 450 }});

Plotly.newPlot("levelChart", [{{
  type: "bar",
  x: {json.dumps(list(level_counts.keys()))},
  y: {json.dumps(list(level_counts.values()))},
  marker: {{ color: "#3498DB" }}
}}, {{
  title: "各 {class_label} 的新颖性分布",
  xaxis: {{ title: "{class_label}" }},
  yaxis: {{ title: "vOTU 数量" }},
  height: 450,
}});

Plotly.newPlot("hostChart", [{{
  type: "bar",
  x: {json.dumps(list(host_counts.keys()))},
  y: {json.dumps(list(host_counts.values()))},
  marker: {{ color: "#8E44AD" }}
}}, {{
  title: "宿主分布",
  xaxis: {{ title: "宿主类别" }},
  yaxis: {{ title: "vOTU 数量" }},
  height: 400,
}});

Plotly.newPlot("freqChart", [{{
  type: "histogram",
  x: {json.dumps(stats.get("frequency_list", []))},
  xbins: {{ start: 1, size: 1 }},
  marker: {{ color: "#1ABC9C", line: {{ color: "#fff", width: 1 }} }}
}}, {{
  title: "病毒频率分布 (检出样本数)",
  xaxis: {{ title: "样本数", dtick: 1 }},
  yaxis: {{ title: "vOTU 数量" }},
  height: 400,
}});
</script>

<p class="footer">validate_novel_viruses.py v2.0 — 病毒宏基因组分析全流程框架</p>
</body>
</html>"""

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)


def main():
    p = argparse.ArgumentParser(description="基于分类层级判断病毒新颖性 v2.0")

    p.add_argument("-i", "--input", required=True, type=Path,
                   help="vOTU centroids FASTA (final_centroids.fasta)")
    p.add_argument("--taxonomy", required=True, type=Path,
                   help="分类注释 TSV (final_integrated_classification.tsv)")
    p.add_argument("--cdhit-known", type=Path,
                   help="CD-HIT known_association.tsv (已知 centroids 标记)")
    p.add_argument("--clusters-tsv", type=Path,
                   help="vclust_clusters.tsv (计算病毒在各样本中的频率)")
    p.add_argument("--host", type=Path,
                   help="宿主预测 TSV (ensemble_host_summary.tsv)")
    p.add_argument("--known-summary", type=Path,
                   help="已知病毒 summary TSV (可选, 补充丰度信息)")
    p.add_argument("--ref-info", type=Path,
                   help="参考元数据 TSV (可选, 补充物种名)")

    p.add_argument("-o", "--output-dir", default="07_Validation", type=Path,
                   help="输出目录 (默认 07_Validation)")
    p.add_argument("--log-level", default="INFO", help="日志级别")
    p.add_argument("--force", action="store_true", help="强制重新生成报告")

    args = p.parse_args()
    logger = setup_logger(args.log_level)

    os.makedirs(args.output_dir, exist_ok=True)
    work_dir = os.path.abspath(args.output_dir)

    # 断点续传: 输出文件已存在则跳过
    report_path = os.path.join(work_dir, "validation_report.html")
    if not args.force and os.path.exists(report_path) and os.path.getsize(report_path) > 100:
        logger.info("输出文件已存在, 跳过 (--force 强制重跑)")
        logger.info("  %s", report_path)
        return

    # 校验输入
    for path, name in [(args.input, "input"), (args.taxonomy, "taxonomy")]:
        if not os.path.exists(str(path)):
            logger.error(f"文件不存在: {name} = {path}")
            sys.exit(1)

    # ── Step 1: 加载数据 ──
    logger.info("=" * 60)
    logger.info("Step 1: 加载输入数据")

    logger.info(f"  加载 vOTU: {args.input}")
    votu_ids = set()
    votu_seqs = {}
    for rec in SeqIO.parse(str(args.input), "fasta"):
        votu_ids.add(rec.id)
        votu_seqs[rec.id] = str(rec.seq)
    logger.info(f"    {len(votu_ids):,} 条 vOTU")

    logger.info(f"  加载分类: {args.taxonomy}")
    taxonomy = pl.read_csv(str(args.taxonomy), separator="\t", null_values=["NA", "N/A", ""])
    logger.info(f"    {taxonomy.height:,} 行, {taxonomy.width} 列")

    # 加载 CD-HIT 已知
    cdhit_known_ids = set()
    cdhit_ref_map = {}
    if args.cdhit_known and os.path.exists(str(args.cdhit_known)):
        logger.info(f"  加载 CD-HIT known: {args.cdhit_known}")
        cdhit = pl.read_csv(str(args.cdhit_known), separator="\t", null_values=["NA", "N/A", ""])
        if "contig_id" in cdhit.columns and "ref_accession" in cdhit.columns:
            for row in cdhit.iter_rows(named=True):
                cdhit_known_ids.add(row["contig_id"])
                cdhit_ref_map[row["contig_id"]] = row["ref_accession"]
        logger.info(f"    CD-HIT known: {len(cdhit_known_ids):,} 条")
    else:
        logger.info("  无 CD-HIT known (全部通过 taxonomy 判断)")

    # 加载 clusters.tsv → 计算病毒频率 (每 cluster 去重样本数)
    cluster_frequency = {}  # {centroid_id: frequency}
    cluster_samples = {}    # {centroid_id: [sample_names]}
    if args.clusters_tsv and os.path.exists(str(args.clusters_tsv)):
        logger.info(f"  加载 clusters.tsv: {args.clusters_tsv}")
        clust_map = defaultdict(list)  # {cluster_id: [member_ids]}
        with open(str(args.clusters_tsv)) as f:
            f.readline()  # skip header
            for line in f:
                p = line.strip().split()
                if len(p) >= 2:
                    clust_map[p[1].strip()].append(p[0].strip())

        # 对每个 cluster, 提取所有成员 contig 的样本前缀
        for cid_name, members in clust_map.items():
            # 找代表 (第一个成员)
            rep = members[0]
            samples = set()
            for m in members:
                # 提取样本名: 第一个 _ 之前
                sample = m.split('_')[0]
                samples.add(sample)
            cluster_frequency[rep] = len(samples)
            cluster_samples[rep] = sorted(samples)
        logger.info(f"    {len(cluster_frequency):,} 个 cluster 有频率信息")

    # 加载宿主
    host_df = None
    if args.host and os.path.exists(str(args.host)):
        logger.info(f"  加载宿主: {args.host}")
        host_df = pl.read_csv(str(args.host), separator="\t", null_values=["NA", "N/A", ""])
        logger.info(f"    {host_df.height:,} 行")

    # 加载参考元数据 (可选, 补充物种名)
    ref_info = None
    if args.ref_info and os.path.exists(str(args.ref_info)):
        logger.info(f"  加载参考元数据: {args.ref_info}")
        ref_info = pl.read_csv(str(args.ref_info), separator="\t", null_values=["NA", "N/A", ""])

    # ── Step 2: 分类判断 ──
    logger.info("=" * 60)
    logger.info("Step 2: 分类层级判断")

    # 合并 classification + CD-HIT 标记
    classifications = []
    counts = {"known": 0, "novel_species": 0, "novel_genus": 0,
              "novel_family": 0, "novel_order": 0, "truly_novel": 0}

    for row in taxonomy.iter_rows(named=True):
        cid = row.get("contig_id", "")

        # CD-HIT 标记
        is_cdhit = cid in cdhit_known_ids

        # taxonomy 数据
        tax_data = {}
        for level in TAX_LEVELS:
            v = row.get(level, None)
            tax_data[level] = v

        label, cat, method = classify_by_taxonomy({
            "cdhit_known": is_cdhit,
            **tax_data,
        })

        # 如果 CD-HIT 标记了但 taxonomy 无 Species, 从 ref_info 补充
        ref_species = ""
        if is_cdhit and cid in cdhit_ref_map:
            ref_acc = cdhit_ref_map[cid]
            if ref_info is not None:
                species_col = None
                for c in ["Species_NCBI", "Species_ICTV", "Species", "description"]:
                    if c in ref_info.columns:
                        species_col = c
                        break
                if species_col:
                    match = ref_info.filter(pl.col("Accession").cast(pl.Utf8) == ref_acc)
                    if match.height > 0:
                        ref_species = str(match[species_col][0])

        freq = cluster_frequency.get(cid, 0)
        samp_list = ",".join(cluster_samples.get(cid, []))
        classifications.append({
            "contig_id": cid,
            "category": cat,
            "label": label,
            "method": method,
            **tax_data,
            "cdhit_known": is_cdhit,
            "frequency": freq,
            "samples": samp_list,
            "ref_species": ref_species,
            "sequence": votu_seqs.get(cid, ""),
        })
        counts[cat] = counts.get(cat, 0) + 1

    class_df = pl.DataFrame(classifications)

    # 统计
    known_total = counts["known"]
    novel_total = counts["novel_species"] + counts["novel_genus"]
    truly_total = counts["novel_family"] + counts["novel_order"] + counts["truly_novel"]
    total = sum(counts.values())

    logger.info(f"  ★ known:               {known_total:,} (CD-HIT + Species ≠ NA)")
    logger.info(f"  ★★ novel_species:      {counts['novel_species']:,} (Genus ≠ NA, Species = NA)")
    logger.info(f"  ★★ novel_genus:        {counts['novel_genus']:,} (Family ≠ NA, Genus = NA)")
    logger.info(f"  ★★★ novel_family+:     {counts['novel_family'] + counts['novel_order']:,} (Class/Order ≠ NA)")
    logger.info(f"  ★★★ truly_novel:       {counts['truly_novel']:,} (全部 NA)")

    # ── Step 3: 补充分类/宿主 ──
    logger.info("=" * 60)
    logger.info("Step 3: 补充宿主信息")

    if host_df is not None:
        host_cols = ["contig_id"]
        for c in ["Final_Host", "Decision_Method", "Host_ICTV"]:
            if c in host_df.columns:
                host_cols.append(c)
        class_df = class_df.join(host_df.select(host_cols), on="contig_id", how="left", suffix="_host")
        logger.info(f"  已补充宿主预测")

    # ── Step 4: 统计收集 ──
    # 宿主分布
    host_counts = {}
    if "Final_Host" in class_df.columns:
        for val in class_df["Final_Host"].to_list():
            v = val if val and val != "NA" else "Unknown"
            host_counts[v] = host_counts.get(v, 0) + 1
    else:
        host_counts = {"Unknown": total}

    # 分类层级分布 (Class 级别)
    level_counts = {}
    if "Class" in class_df.columns:
        for val in class_df["Class"].to_list():
            if val and val != "NA":
                level_counts[val] = level_counts.get(val, 0) + 1
    level_counts = dict(sorted(level_counts.items(), key=lambda x: -x[1])[:15])

    # 频率统计
    frequency_list = class_df["frequency"].to_list() if "frequency" in class_df.columns else []
    all_samples = set()
    for s in class_df["samples"].to_list() if "samples" in class_df.columns else []:
        if s:
            for x in s.split(","):
                if x.strip():
                    all_samples.add(x.strip())

    stats = {
        "counts": counts,
        "host_counts": host_counts,
        "level_counts": level_counts,
        "class_label": "Class",
        "frequency_list": frequency_list,
        "total_samples": len(all_samples),
    }

    # ── Step 5: 输出 ──
    logger.info("=" * 60)
    logger.info("Step 4: 输出文件")

    # novel_viruses.annotated.tsv (全部 vOTU)
    annot_out = os.path.join(work_dir, "novel_viruses.annotated.tsv")
    annot_cols = [c for c in class_df.columns if c != "sequence"]
    class_df.select(annot_cols).write_csv(annot_out, separator="\t")
    logger.info(f"  novel_viruses.annotated.tsv: {class_df.height} 行")

    # final_virus_catalog.fasta
    catalog_out = os.path.join(work_dir, "final_virus_catalog.fasta")
    with open(catalog_out, "w") as cf:
        for row in class_df.iter_rows(named=True):
            cid = row["contig_id"]
            cat = row["category"]
            seq = row.get("sequence", "")
            if seq:
                cf.write(f">{cid}|{cat}\n")
                for i in range(0, len(seq), 60):
                    cf.write(seq[i:i+60] + "\n")
    logger.info(f"  final_virus_catalog.fasta: {class_df.height} 条")

    # validation_report.html
    report_out = os.path.join(work_dir, "validation_report.html")
    generate_html_report(stats, class_df, report_out)
    logger.info(f"  validation_report.html: {report_out}")

    # ── 汇总 ──
    logger.info("=" * 60)
    logger.info("验证完成!")
    logger.info(f"  已知病毒:      {known_total:,} 个 (CD-HIT锚定 + Species确定)")
    logger.info(f"  新种/新属:     {novel_total:,} 个 (Genus已知但Species新)")
    logger.info(f"  新科及以上:    {truly_total:,} 个 (需要进一步验证)")
    if frequency_list:
        freq_1 = sum(1 for f in frequency_list if f == 1)
        freq_2_5 = sum(1 for f in frequency_list if 2 <= f <= 5)
        freq_6plus = sum(1 for f in frequency_list if f >= 6)
        logger.info(f"  病毒频率:      singletons={freq_1:,}, 2-5样本={freq_2_5:,}, ≥6样本={freq_6plus:,}")
        logger.info(f"  总样本数:      {len(all_samples):,}")
    logger.info(f"  输出目录:      {work_dir}")


if __name__ == "__main__":
    main()
