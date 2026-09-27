import csv
import os
import subprocess

BASE = "/home/zhangwenda/data-test/out"
V = f"{BASE}/09b_Analysis_Verify/virus_validation"

valid = set()
with open(f"{V}/rescue_evidence_scored.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        valid.add(r["contig_id"])
print(f"scored.tsv: {len(valid)} 条")

# 从备份重建 499 = 备份 txt 里有、scored 里没有
s = set()
for v in ["keep", "review", "drop"]:
    p = f"{V}/{v}_candidates.txt.bak_bl499_20260912"
    for l in open(p, errors="replace"):
        l = l.strip()
        if l and l.split("\t")[0] not in valid:
            s.add(l.split("\t")[0])
with open("/tmp/stale_499.txt", "w") as f:
    f.write("\n".join(sorted(s)) + "\n")
print(f"499 清单重建: {len(s)} 条")

# 限定范围检索: 排除 tmp/ 与 .bak*
dirs = ["06_HostPrediction", "07_Checkv", "08_Rescue",
        "09_Virome_Analysis", "09b_Analysis_Verify", "10_Reports"]
print("\n=== 残留检索 (排除 tmp/ 与 .bak*) ===")
for d in dirs:
    cmd = ["grep", "-rlFf", "/tmp/stale_499.txt", f"{BASE}/{d}",
           "--exclude=*.bak*", "--exclude-dir=tmp", "--exclude=*.log"]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=180).stdout.strip()
    except subprocess.TimeoutExpired:
        print(f"  {d}: [超时]")
        continue
    files = [x for x in out.split("\n") if x]
    print(f"  {d}: 命中 {len(files)} 个文件")
    for x in files[:8]:
        print(f"      {x.replace(BASE+'/','')}")
