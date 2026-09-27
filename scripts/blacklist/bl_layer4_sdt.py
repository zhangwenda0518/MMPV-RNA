"""层4: 重算受影响属的 SDT 矩阵。

范围: Alphaorpheovirus(14) / Hatfieldvirus(2) / Tethysvirus(50)
参数对齐 virome_pipeline.py::run_sdt_matrix。
"""
import os
import subprocess
import sys
import time

SDT_SCRIPT = "/home/zhangwenda/MMPV-RNA/virome_analysis_pipeline/sdt_genus_matrix.py"
ACVDIR = "/home/zhangwenda/data-test/out/09b_Analysis_Verify"
SDT_OUT = f"{ACVDIR}/SDT_matrix"

GENERA = ["Tethysvirus", "Alphaorpheovirus", "Hatfieldvirus"]  # 按序列数降序
THREADS = 8

os.makedirs(SDT_OUT, exist_ok=True)

print("=== 层4 重算 SDT 属矩阵 ===")
results = []
for g in GENERA:
    od = os.path.join(SDT_OUT, g)
    os.makedirs(od, exist_ok=True)
    t0 = time.time()
    cmd = ["python3", SDT_SCRIPT,
           "--analysis_dir", ACVDIR, "--genus", g, "-o", od,
           "--threads", str(THREADS), "--short_labels",
           "--palette", "cividis", "--resume",
           "--other_count", "0", "--target_cap", "20"]
    print(f"\n[{g}] CMD: {' '.join(cmd)}", flush=True)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
        dt = time.time() - t0
        print(f"[{g}] rc={r.returncode}  {dt:.0f}s", flush=True)
        if r.returncode != 0:
            print(f"[{g}] STDERR: {r.stderr[-600:]}", flush=True)
        results.append((g, r.returncode, dt))
    except Exception as e:
        print(f"[{g}] 异常: {e}", flush=True)
        results.append((g, 99, 0))

print("\n=== SDT 汇总 ===")
for g, rc, dt in results:
    od = os.path.join(SDT_OUT, g)
    files = os.listdir(od) if os.path.isdir(od) else []
    mark = "OK" if rc == 0 else "FAIL"
    print(f"  {mark:<5} {g:<20} rc={rc}  {dt:.0f}s  产物: {files}")
