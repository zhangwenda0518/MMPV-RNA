#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_identfilter_data.py — 汇总「鉴定 / 过滤」双簇条形图的数据中间件。

数据源（只读，全部来自本地镜像 reports/_discovery_sci_20261012/_lt/10_Reports/）：
  ident_summary.tsv       11 列 = Sample + All_Candidate + 9 个病毒工具 + viroid_blastn；末行 TOTAL。
                          All_Candidate = 各**病毒**工具 id 的并集；viroid_blastn **不计入**
                          （2026-10-13 在 02a_Identification 逐样本 comm 验证：viroid id 全部落在并集外）。
  filter_per_sample.tsv   395 行 = Sample / Uniprot_Filtered / Cdd_Filtered。
                          ⚠ 该文件是 v18 本地新建（服务器无此文件，属遗留 R8）；
                            它替代了恒为 0 的 filter_summary.tsv。

产出：
  identfilter_data/ident_filter_summary.tsv   长表 group / track / category / value
  identfilter_data/ident_filter_meta.json     并集、保留率、零命中工具等口径元数据
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REP = ROOT / "reports" / "_discovery_sci_20261012" / "_lt" / "10_Reports"
OUT = Path(__file__).resolve().parent / "identfilter_data"

VIRUS_TOOLS = ["genomad", "blast", "metabuli", "virsorter2", "viralverify",
               "virhunter", "virbot", "viralm", "rdrpcatch"]
VIROID_TOOL = "viroid_blastn"


def read_total_row(path: Path):
    """返回 ident_summary.tsv 的表头与 TOTAL 行（dict）。"""
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        total = None
        n_sample = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "TOTAL":
                total = dict(zip(header, parts))
            elif parts[0] != "summary_plots":
                n_sample += 1
    if total is None:
        sys.exit(f"❌ {path} 里找不到 TOTAL 行")
    return header, total, n_sample


def main():
    ident_p = REP / "ident_summary.tsv"
    filt_p = REP / "filter_per_sample.tsv"
    for p in (ident_p, filt_p):
        if not p.is_file():
            sys.exit(f"❌ 缺少数据源: {p}")

    header, total, n_sample = read_total_row(ident_p)

    # ── 口径自检 1：表头必须含 viroid_blastn（旧本地副本只有 10 列，会静默漏掉类病毒轨道）──
    if VIROID_TOOL not in header:
        sys.exit(f"❌ ident_summary.tsv 缺列 {VIROID_TOOL} —— 这是 2026-10-13 前的旧版，"
                 f"请从服务器重新拉取（服务器版 11 列，md5 f9478309930ccbf5e9b1bd8dd281226e）")
    missing = [t for t in VIRUS_TOOLS if t not in header]
    if missing:
        sys.exit(f"❌ ident_summary.tsv 缺病毒工具列: {missing}")

    # ── 口径自检 2：并集必须落在「工具和」以内（并集 > 和 = 说明口径读错）──
    vals = {t: int(total[t]) for t in VIRUS_TOOLS}
    union = int(total["All_Candidate"])
    viroid = int(total[VIROID_TOOL])
    sum_tools = sum(vals.values())
    if union > sum_tools:
        sys.exit(f"❌ 并集 {union} > 各工具之和 {sum_tools} —— All_Candidate 口径异常")

    # ── 过滤侧：只取 Uniprot 保留总数（用户 2026-10-13 选定「只展示保留总数」）──
    with filt_p.open(encoding="utf-8") as fh:
        fh.readline()
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    uni = sum(int(r[1]) for r in rows)
    cdd = sum(int(r[2]) for r in rows)
    n_filt_sample = len(rows)

    # ── 长表 ──
    recs = []
    for t in VIRUS_TOOLS:
        if vals[t] > 0:
            recs.append(dict(group="Identification", track="virus", category=t, value=vals[t]))
    if viroid > 0:
        recs.append(dict(group="Identification", track="viroid", category=VIROID_TOOL, value=viroid))
    recs.append(dict(group="Filtering", track="virus", category="Retained (02b_Filter)", value=uni))
    recs.sort(key=lambda r: (r["group"] != "Identification", -r["value"]))

    OUT.mkdir(parents=True, exist_ok=True)
    tsv = OUT / "ident_filter_summary.tsv"
    with tsv.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("group\ttrack\tcategory\tvalue\n")
        for r in recs:
            fh.write(f"{r['group']}\t{r['track']}\t{r['category']}\t{r['value']}\n")

    zero_tools = [t for t in VIRUS_TOOLS if vals[t] == 0]
    meta = dict(
        source_ident="ident_summary.tsv (TOTAL 行)",
        source_filter="filter_per_sample.tsv (Uniprot_Filtered 求和)",
        n_samples_ident=n_sample,
        n_samples_filter=n_filt_sample,
        union_all_candidate=union,
        sum_virus_tools=sum_tools,
        overlap_virus_tools=sum_tools - union,
        viroid_blastn=viroid,
        viroid_in_union=False,
        retained_uniprot=uni,
        retained_pct_of_union=round(100.0 * uni / union, 2),
        cdd_filtered_not_shown=cdd,
        active_virus_tools=[t for t in VIRUS_TOOLS if vals[t] > 0],
        zero_virus_tools=zero_tools,
        per_tool_values=vals,
    )
    with (OUT / "ident_filter_meta.json").open("w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)

    print(f"✅ {tsv}")
    print(f"   鉴定并集 {union:,}｜各工具和 {sum_tools:,}（重叠 {sum_tools-union:,}）"
          f"｜类病毒 {viroid:,}（不在并集）")
    print(f"   过滤保留 {uni:,} = 并集的 {meta['retained_pct_of_union']}%")
    print(f"   零命中工具 {len(zero_tools)} 个: {', '.join(zero_tools) if zero_tools else '无'}")
    print(f"   样本数: 鉴定 {n_sample} / 过滤 {n_filt_sample}")


if __name__ == "__main__":
    main()
