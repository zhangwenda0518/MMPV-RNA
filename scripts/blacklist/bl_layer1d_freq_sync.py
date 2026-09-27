"""层1d: 清理 rescue_detection/frequency_table.tsv 中的被删 contig。

发现: 该表首列 contig_id, 9414 行 (= 09_Plant_virus 总数), 命中 662 清单 499 条。
      报告 read 它为 rescue_detection_summary 源之一。

附带: plots/ 目录若有基于该表的图, 仅日志提示 (不重绘, 保持数据一致即可)。
"""
import csv
import os
import shutil
import time

P = "/home/zhangwenda/data-test/out/09_Virome_Analysis/rescue_detection/frequency_table.tsv"
KILL = set(l.strip() for l in open("/tmp/blacklist_kill_ids.txt") if l.strip())
STAMP = time.strftime("%Y%m%d")

bak = f"{P}.bak_bl662_{STAMP}"
shutil.copy2(P, bak)

with open(P, errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    hdr = next(rd)
    ci = hdr.index("contig_id")
    kept, dropped = [], 0
    for r in rd:
        if not r:
            continue
        if r[ci].strip() in KILL:
            dropped += 1
        else:
            kept.append(r)

with open(P, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(hdr)
    w.writerows(kept)

print(f"[过滤] frequency_table.tsv")
print(f"    {len(kept)+dropped} -> {len(kept)} (删 {dropped})")
print(f"    备份: {bak}")

with open(P, errors="replace") as f:
    rd = csv.reader(f, delimiter="\t"); h = next(rd); ci = h.index("contig_id")
    left = sum(1 for r in rd if r and r[ci].strip() in KILL)
print(f"    残留复核: {left} {'✓' if left == 0 else '✗'}")
