"""扫描 662 条待删 ID 在管线各层产物中的分布。"""
import os
import subprocess

KILL = "/tmp/blacklist_kill_ids.txt"
ids = set()
with open(KILL) as f:
    for l in f:
        l = l.strip()
        if l:
            ids.add(l)
print(f"待删 ID: {len(ids)} 条\n")

OUT = "/home/zhangwenda/data-test/out"

# 关键产物清单 (path, 是否 fasta)
targets = [
    (f"{OUT}/06_HostPrediction/host_classified_fasta/Plant.classified.fasta", "fasta"),
    (f"{OUT}/06_HostPrediction/C9_ICTV_result/Plant.classified.fasta", "fasta"),
    (f"{OUT}/06_HostPrediction/ensemble_host_summary.tsv", "tsv"),
]
# 07/08/09/10 下自动发现
for sub in ["07_Checkv", "08_Rescue", "09_Virome_Analysis", "09b_Analysis_Verify"]:
    d = os.path.join(OUT, sub)
    if not os.path.isdir(d):
        continue
    for root, dirs, files in os.walk(d):
        for fn in files:
            if fn.endswith((".fasta", ".fa", ".fna", ".tsv")):
                targets.append((os.path.join(root, fn),
                                "fasta" if fn.endswith((".fasta", ".fa", ".fna")) else "tsv"))

print(f"{'产物':<72} {'命中':>7}")
print("-" * 82)
hits = []
for path, kind in targets:
    if not os.path.exists(path):
        continue
    if kind == "fasta":
        n = 0
        with open(path, errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    if line[1:].split()[0] in ids:
                        n += 1
    else:
        # TSV: 按第一列精确匹配
        n = 0
        with open(path, errors="replace") as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                c = line.split("\t", 1)[0].strip().strip('"')
                if c in ids:
                    n += 1
    if n > 0:
        rel = path.replace(OUT + "/", "")
        hits.append((rel, n))
        print(f"{rel:<72} {n:>7}")

print("-" * 82)
print(f"命中产物的文件数: {len(hits)}")
