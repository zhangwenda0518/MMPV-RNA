"""层1g: 剔除 ensemble_host_summary.tsv 中的幽灵行 (Final_Host=Plant 但无对应序列)。

根因:
  ensemble_host_summary.tsv 由分类表(--tax, 05_Taxonomy/final_integrated_classification.tsv)逐行生成;
  host_classified_fasta/*.fasta 由输入序列按 _resolve_host 分组。
  两者源集不同 -> summary 里可能有分类表有、输入序列没有的"幽灵行"。
  实测: Plant 幽灵行 = HE862273.1, MK032712.1
        (VITAP 分类, Order=Geplafuvirales, Family/Genus/Species=NA, Determination_Level=Order(via Order))
  这 2 条不在任何 host fasta (含 Unknown), 即输入序列中根本不存在。

目标: 使 06 报告 Plant 数 (20056) 与下游序列层 (20054) 一致。

策略: 动态计算 Plant 幽灵行 (Final_Host=='Plant' 且 contig_id 不在 Plant.classified.fasta), 备份后剔除。
     非 Plant 的幽灵行一并统计报告 (本次若有也会列出, 但只删 Plant? 不 : 全删以保一致)。
     注意: 仅对 Final_Host 明确的幽灵行操作; Unknown 桶的差值是源集差异, 不动。
"""
import csv
import os
import shutil
import time

R = "/home/zhangwenda/data-test/out/06_HostPrediction"
SP = f"{R}/ensemble_host_summary.tsv"
FA = f"{R}/host_classified_fasta/Plant.classified.fasta"
STAMP = time.strftime("%Y%m%d")

# 1) 读 Plant fasta 的 id 集合
fa_ids = set()
with open(FA, errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            fa_ids.add(line[1:].split()[0])
print(f"Plant.classified.fasta: {len(fa_ids)} 条")

# 2) 读 summary, 找出 Final_Host==Plant 但不在 fasta 的行
with open(SP, errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    hdr = next(rd)
    ci = hdr.index("contig_id")
    fi = hdr.index("Final_Host")
    rows = list(rd)

# 先扫全部 host 的幽灵行 (信息性)
import glob as _g
host_ids = {}
for fp in _g.glob(f"{R}/host_classified_fasta/*.classified.fasta"):
    hn = os.path.basename(fp).replace(".classified.fasta", "")
    ids = set()
    with open(fp, errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    host_ids[hn] = ids
print("各 host: summary 计数 vs fasta 计数")
from collections import Counter
sm = Counter(r[fi].strip() for r in rows if r and len(r) > fi)
for hn in sorted(host_ids):
    s, fa = sm.get(hn, 0), len(host_ids[hn])
    flag = "" if s == fa else f"  <== 差 {s-fa}"
    print(f"  {hn:<12} {s:>7} vs {fa:>7}{flag}")
print()

ghosts = [r for r in rows if r and len(r) > max(ci, fi)
          and r[fi].strip() == "Plant" and r[ci].strip() not in fa_ids]
print(f"summary 总行: {len(rows)}")
print(f"Plant 幽灵行 (summary有/fasta无): {len(ghosts)}")
for g in ghosts:
    print(f"   - {g[ci]}")

if not ghosts:
    print("无幽灵行, 无需处理")
    raise SystemExit

# 3) 备份 + 剔除
bak = f"{SP}.bak_ghost_{STAMP}"
shutil.copy2(SP, bak)
ghost_ids = set(g[ci].strip() for g in ghosts)
kept = [r for r in rows if not (r and len(r) > ci and r[ci].strip() in ghost_ids)]

with open(SP, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(hdr)
    w.writerows(kept)

print(f"\n[剔除] {len(rows)} -> {len(kept)} (删 {len(rows)-len(kept)})")
print(f"备份: {bak}")

# 4) 复核
with open(SP, errors="replace") as f:
    rd = csv.reader(f, delimiter="\t"); h = next(rd)
    ci2 = h.index("contig_id"); fi2 = h.index("Final_Host")
    n_plant = sum(1 for r in rd if len(r) > fi2 and r[fi2].strip() == "Plant")
print(f"复核: 现 Plant={n_plant} (fasta={len(fa_ids)}) {'✓' if n_plant==len(fa_ids) else '✗'}")
