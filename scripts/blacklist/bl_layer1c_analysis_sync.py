"""层1c: 同步清理 09_Virome_Analysis/all_plant_analysis 中的 662 条。

发现: 报告 _collect 里 "All plant: N seqs" 读的是:
      analysis/all_plant_analysis/All_plant.viruses_info.tsv

目标:
  1. All_plant.viruses_info.tsv          (contig_id 首列, 20716 -> 20054)
  2. plant_virus_cluster_summary.tsv     (cluster_id 首列, contig_id==cluster_id 逐簇)
  3. all_plant_viruses_genus_summary.tsv (从 1 重算 Genus 聚合)

策略: 备份 -> 过滤 -> 重算聚合 -> 校验
"""
import csv
import os
import shutil
import time
from collections import defaultdict

D = "/home/zhangwenda/data-test/out/09_Virome_Analysis/all_plant_analysis"
KILL = set(l.strip() for l in open("/tmp/blacklist_kill_ids.txt") if l.strip())
STAMP = time.strftime("%Y%m%d")

print(f"662 清单: {len(KILL)} 条\n")


def filt_tsv(path, id_col_name):
    bak = f"{path}.bak_bl662_{STAMP}"
    shutil.copy2(path, bak)
    with open(path, errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ci = hdr.index(id_col_name)
        kept, dropped = [], 0
        for r in rd:
            if not r:
                continue
            if r[ci].strip() in KILL:
                dropped += 1
            else:
                kept.append(r)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(hdr)
        w.writerows(kept)
    print(f"[过滤] {os.path.basename(path)}")
    print(f"    {len(kept)+dropped} -> {len(kept)} (删 {dropped})")
    print(f"    备份: {bak}")
    return hdr, kept


# ── 1. All_plant.viruses_info.tsv ──
hdr1, kept1 = filt_tsv(f"{D}/All_plant.viruses_info.tsv", "contig_id")
print()

# ── 2. plant_virus_cluster_summary.tsv ──
hdr2, kept2 = filt_tsv(f"{D}/plant_virus_cluster_summary.tsv", "cluster_id")
print()

# ── 3. 重算 genus_summary ──
gs_path = f"{D}/all_plant_viruses_genus_summary.tsv"
old_hdr = open(gs_path, errors="replace").readline().rstrip("\n").split("\t")
print(f"[重算] all_plant_viruses_genus_summary.tsv (旧表头: {old_hdr})")

# 找 Genus / n_contigs 列
gi = hdr1.index("Genus") if "Genus" in hdr1 else None
if gi is None:
    print("    [跳过] All_plant 无 Genus 列, 无法重算")
else:
    from collections import Counter
    cnt = Counter()
    for r in kept1:
        g = r[gi].strip() if len(r) > gi else ""
        cnt[g if g else "NA"] += 1
    bak = f"{gs_path}.bak_bl662_{STAMP}"
    shutil.copy2(gs_path, bak)
    with open(gs_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        # 保持原表头结构 (其余列若无法重算则置 0)
        w.writerow(old_hdr)
        for g, n in cnt.most_common():
            row = [g, str(n)] + ["0"] * (len(old_hdr) - 2)
            w.writerow(row)
    print(f"    属数 {len(cnt)}, 总 contig {sum(cnt.values())}")
    print(f"    备份: {bak}")

# ── 校验 ──
print("\n=== 校验 ===")
for p, col in [(f"{D}/All_plant.viruses_info.tsv", "contig_id"),
               (f"{D}/plant_virus_cluster_summary.tsv", "cluster_id")]:
    with open(p, errors="replace") as f:
        rd = csv.reader(f, delimiter="\t"); h = next(rd)
        ci = h.index(col)
        left = sum(1 for r in rd if r and r[ci].strip() in KILL)
    print(f"  {os.path.basename(p)}: 残留 {left} {'✓' if left == 0 else '✗'}")
