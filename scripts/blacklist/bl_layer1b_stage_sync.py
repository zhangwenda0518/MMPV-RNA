"""层1b: 同步清理 06/07 stage 原始产物中的 662 条。

背景: 报告 _collect_hostprediction / _collect_checkv 直接读 stage 原始目录,
      与下游 host_classified_fasta 不同源, 需同步删才能让报告数字归位。

目标:
  1. 06_HostPrediction/ensemble_host_summary.tsv   (contig_id 列, 首列)
  2. 07_Checkv/Plant/completeness.tsv              (contig_id 列, 首列)

策略: 备份 → 过滤 → 校验行数
"""
import os
import shutil
import time

ROOT = "/home/zhangwenda/data-test/out"
KILL = set(l.strip() for l in open("/tmp/blacklist_kill_ids.txt") if l.strip())
STAMP = time.strftime("%Y%m%d")

TARGETS = [
    f"{ROOT}/06_HostPrediction/ensemble_host_summary.tsv",
    f"{ROOT}/07_Checkv/Plant/completeness.tsv",
]

print(f"662 清单: {len(KILL)} 条\n")

for p in TARGETS:
    if not os.path.exists(p):
        print(f"[跳过] 不存在: {p}")
        continue

    bak = f"{p}.bak_bl662_{STAMP}"
    shutil.copy2(p, bak)

    with open(p, errors="replace") as f:
        lines = f.readlines()

    hdr = lines[0]
    body = lines[1:]
    kept = [l for l in body if l.split("\t", 1)[0].strip() not in KILL]
    dropped = len(body) - len(kept)

    # ID 列判定: 首列是否就是 contig_id
    first_hdr = hdr.split("\t")[0].strip()
    if first_hdr not in ("contig_id", "Contig_ID", "qseqid", "ID"):
        print(f"[警告] {os.path.basename(p)} 首列是 '{first_hdr}', 非预期 ID 列, 跳过!")
        os.remove(bak)
        continue

    with open(p, "w", encoding="utf-8") as f:
        f.write(hdr)
        f.writelines(kept)

    print(f"[处理] {os.path.basename(p)}")
    print(f"    {len(body)} -> {len(kept)} (删除 {dropped})")
    print(f"    备份: {bak}")

    # 复核: 残留应为 0
    with open(p, errors="replace") as f:
        f.readline()
        left = sum(1 for l in f if l.split("\t", 1)[0].strip() in KILL)
    print(f"    残留复核: {left} {'✓' if left == 0 else '✗ 异常!'}")
    print()
