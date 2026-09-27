"""层1f: 清理 08_Rescue final_centroids.fasta 中的被删 contig。

发现: 报告 _generate_plant_virus_summary 从这两个 fasta 合并构建 plant_virus_summary:
  1. 08_Rescue/known/centroids/final_centroids.fasta   (免拯救, 199 -> 198)
  2. 08_Rescue/Plant/centroids/final_centroids.fasta   (rescued, 9215 -> 8717)

注: 08_Rescue/HQ_plant_viruses.fasta 已是 8915 (层1 已同步), 此处补齐其上游源。

策略: 备份 -> 过滤 (保留可选第二行) -> 校验
"""
import os
import shutil
import time

R = "/home/zhangwenda/data-test/out/08_Rescue"
KILL = set(l.strip() for l in open("/tmp/blacklist_kill_ids.txt") if l.strip())
STAMP = time.strftime("%Y%m%d")

print(f"662 清单: {len(KILL)} 条\n")

TARGETS = [
    (f"{R}/known/centroids/final_centroids.fasta", "免拯救"),
    (f"{R}/Plant/centroids/final_centroids.fasta", "rescued"),
]

for p, label in TARGETS:
    if not os.path.exists(p):
        print(f"[跳过] 不存在: {p}")
        continue

    bak = f"{p}.bak_bl662_{STAMP}"
    shutil.copy2(p, bak)

    kept_lines = []
    n_tot = n_drop = 0
    cur_keep = False
    with open(p, errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                n_tot += 1
                cid = line[1:].split()[0]
                cur_keep = cid not in KILL
                if not cur_keep:
                    n_drop += 1
            if cur_keep:
                kept_lines.append(line)

    with open(p, "w", encoding="utf-8") as f:
        f.writelines(kept_lines)

    n_kept = n_tot - n_drop
    print(f"[过滤] {label}  {os.path.basename(os.path.dirname(os.path.dirname(p)))}/final_centroids.fasta")
    print(f"    {n_tot} -> {n_kept} (删 {n_drop})")
    print(f"    备份: {bak}")

    # 校验: 序列行数应等于 header 数
    with open(p, errors="replace") as f:
        h = sum(1 for l in f if l.startswith(">"))
    print(f"    复核: 现 {h} 条 {'✓' if h == n_kept else '✗'}")
    print()
