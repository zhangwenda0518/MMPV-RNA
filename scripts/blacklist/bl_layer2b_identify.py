"""层2 步骤B: 重跑 ACVirus identify (新病毒鉴定)。

基于重跑后的 acvirus_classify 结果。
"""
import os
import subprocess
import sys
import time

HOME = os.path.expanduser("~")
ACV = f"{HOME}/MMPV-RNA/biosoft/ACVirus/cli.py"
DB = f"{HOME}/database/virus-db/acvirus_db"
ACVDIR = "/home/zhangwenda/data-test/out/09b_Analysis_Verify"
CLS = f"{ACVDIR}/acvirus_classify"
IDT = f"{ACVDIR}/acvirus_identify"

THREADS = 60
os.makedirs(IDT, exist_ok=True)
out_csv = f"{IDT}/identification_result_score.csv"

print(f"输入 classify: {CLS}/final_result_with_confidence.tsv")
print(f"输出: {out_csv}")
print(f"\n=== ACVirus identify (threads={THREADS}) ===")
t0 = time.time()
cmd = ["python3", ACV, "identify", "--classify_out", CLS, "--data_path", DB,
       "--out", out_csv, "--min_cluster", "3", "-t", str(THREADS)]
print("CMD:", " ".join(cmd))
r = subprocess.run(cmd, capture_output=True, text=True)
print(f"rc={r.returncode}  耗时 {time.time()-t0:.0f}s")
if r.returncode != 0:
    print("STDERR:", r.stderr[-3000:])
    sys.exit(1)
print("STDOUT 尾部:", r.stdout[-1200:])

if os.path.exists(out_csv):
    with open(out_csv, errors="replace") as f:
        rows = sum(1 for _ in f) - 1
    print(f"\n产物: {out_csv}")
    print(f"记录数: {rows}")
else:
    print("\n警告: 未生成产物")
