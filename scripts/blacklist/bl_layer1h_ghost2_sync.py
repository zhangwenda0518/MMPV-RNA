"""层1h: 剔除 ensemble_host_summary.tsv 中剩余的非 Plant 幽灵行 (summary有/Fasta无)。

背景: 层1g 已剔除 Plant 幽灵行 2 条 (HE862273.1, MK032712.1)。
      全 host 扫描发现 Bacteria 也有 1 条幽灵行 MT374858.1 (Peduovirus)。
      07_CheckV 子目录 (由 host fasta 派生) 各 host 计数:
        Algae 46176 / Animal 155944 / Archaea 76 / Bacteria 67513 / Fungi 11641
        Mammalia 33612 / Plant 20054 / Protist 185751 / Unknown 13671  => 534438
      剔除全部幽灵行后 06 summary 应与 07 的 534438 对齐。

范围: 对除 Unknown 外的所有 host, 剔除 Final_Host==host 但 contig_id 不在该 host fasta 的行。
      (Unknown 是源集差异, 不动)
"""
import csv
import glob
import os
import shutil
import time
from collections import Counter

R = "/home/zhangwenda/data-test/out/06_HostPrediction"
SP = f"{R}/ensemble_host_summary.tsv"
STAMP = time.strftime("%Y%m%d")

# 各 host fasta id 集合
host_ids = {}
for fp in glob.glob(f"{R}/host_classified_fasta/*.classified.fasta"):
    hn = os.path.basename(fp).replace(".classified.fasta", "")
    ids = set()
    with open(fp, errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    host_ids[hn] = ids

with open(SP, errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    hdr = next(rd)
    ci = hdr.index("contig_id"); fi = hdr.index("Final_Host")
    rows = list(rd)

# 找幽灵行 (Unknown 除外)
ghost_ids = set()
print("幽灵行明细:")
for r in rows:
    if not r or len(r) <= max(ci, fi):
        continue
    h = r[fi].strip(); cid = r[ci].strip()
    if h == "Unknown":
        continue
    if h in host_ids and cid not in host_ids[h]:
        ghost_ids.add(cid)
        print(f"   - [{h}] {cid}")

if not ghost_ids:
    print("无幽灵行"); raise SystemExit

bak = f"{SP}.bak_ghost2_{STAMP}"
shutil.copy2(SP, bak)

kept = [r for r in rows if not (r and len(r) > ci and r[ci].strip() in ghost_ids)]
with open(SP, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(hdr); w.writerows(kept)

print(f"\n[剔除] {len(rows)} -> {len(kept)} (删 {len(rows)-len(kept)})")
print(f"备份: {bak}\n")

# 复核各 host
with open(SP, errors="replace") as f:
    rd = csv.reader(f, delimiter="\t"); h2 = next(rd); fi2 = h2.index("Final_Host")
    sm = Counter(r[fi2].strip() for r in rd if r and len(r) > fi2)
print("复核 (summary vs fasta):")
for hn in sorted(host_ids):
    s, fa = sm.get(hn, 0), len(host_ids[hn])
    print(f"  {hn:<12} {s:>7} vs {fa:>7} {'✓' if s==fa else '<== 差'}")
print(f"  合计 summary: {sum(sm.values())}  (07 total 目标 534438)")
