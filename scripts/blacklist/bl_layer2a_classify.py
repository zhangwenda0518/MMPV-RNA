"""层2 步骤A: 重跑 ACVirus classify (基于删除后的 class_KEEP.fasta)。

清理旧 classify 产物 → 重新分类 → 得到正确的科清单。
"""
import os
import shutil
import subprocess
import sys
import time

HOME = os.path.expanduser("~")
ACV = f"{HOME}/MMPV-RNA/biosoft/ACVirus/cli.py"
DB = f"{HOME}/database/virus-db/acvirus_db"
ACVDIR = "/home/zhangwenda/data-test/out/09b_Analysis_Verify"
KEEP = f"{ACVDIR}/class_KEEP.fasta"
CLS = f"{ACVDIR}/acvirus_classify"

THREADS = 60

# 检查输入
n = sum(1 for l in open(KEEP) if l.startswith(">"))
print(f"class_KEEP.fasta: {n} 条 (期望 2322)")

# 清理旧 classify 产物 (保留目录)
print(f"\n清理旧产物: {CLS}")
for fn in os.listdir(CLS):
    p = os.path.join(CLS, fn)
    if os.path.isfile(p):
        os.remove(p)
print("  已清空")

# 重跑 classify
print(f"\n=== ACVirus classify (threads={THREADS}) ===")
t0 = time.time()
cmd = ["python3", ACV, "classify", "--contigs", KEEP, "--data_path", DB,
       "--out", CLS, "-t", str(THREADS)]
print("CMD:", " ".join(cmd))
r = subprocess.run(cmd, capture_output=True, text=True)
print(f"rc={r.returncode}  耗时 {time.time()-t0:.0f}s")
if r.returncode != 0:
    print("STDERR:", r.stderr[-3000:])
    sys.exit(1)
print("STDOUT 尾部:", r.stdout[-1000:])

res = f"{CLS}/final_result_with_confidence.tsv"
print(f"\n产物: {res}")
if os.path.exists(res):
    with open(res) as f:
        hdr = f.readline().rstrip().split("\t")
        print("表头:", hdr[:15])
        rows = sum(1 for _ in f)
    print(f"记录数: {rows}")
