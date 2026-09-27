"""层1e: 清理 09b virus_validation 中的被删 contig。

发现: 报告 _collect 里 "Rescue: K/R/D" 和 "Detected: n/N" 读的是:
      _resolve_acv_dir(...)/virus_validation/rescue_evidence_scored.tsv
      实际路径: 09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv

目标:
  1. rescue_evidence_scored.tsv  (contig_id 首列, 9414 -> 8915)
  2. cdd_evidence_report.tsv     (seq 首列,      9414 -> 8915)

策略: 备份 -> 过滤 -> 校验
"""
import csv
import os
import shutil
import time

D = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/virus_validation"
KILL = set(l.strip() for l in open("/tmp/blacklist_kill_ids.txt") if l.strip())
STAMP = time.strftime("%Y%m%d")

print(f"662 清单: {len(KILL)} 条\n")

TARGETS = [
    ("rescue_evidence_scored.tsv", "contig_id"),
    ("cdd_evidence_report.tsv", "seq"),
]

for fn, idcol in TARGETS:
    p = f"{D}/{fn}"
    if not os.path.exists(p):
        print(f"[跳过] 不存在: {fn}")
        continue

    bak = f"{p}.bak_bl662_{STAMP}"
    shutil.copy2(p, bak)

    with open(p, errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        if idcol not in hdr:
            print(f"[跳过] {fn} 无 '{idcol}' 列 (表头: {hdr[:4]})")
            os.remove(bak)
            continue
        ci = hdr.index(idcol)
        kept, dropped = [], 0
        for r in rd:
            if not r:
                continue
            if r[ci].strip() in KILL:
                dropped += 1
            else:
                kept.append(r)

    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(hdr)
        w.writerows(kept)

    print(f"[过滤] {fn}")
    print(f"    {len(kept)+dropped} -> {len(kept)} (删 {dropped})")
    print(f"    备份: {bak}")

    with open(p, errors="replace") as f:
        rd = csv.reader(f, delimiter="\t"); h = next(rd); ci = h.index(idcol)
        left = sum(1 for r in rd if r and r[ci].strip() in KILL)
    print(f"    残留复核: {left} {'✓' if left == 0 else '✗ 异常!'}")
    print()
