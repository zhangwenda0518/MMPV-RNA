"""层3: 重生成报告 (report_pipeline.py)。

在数据层删除 + 树重建 + SDT 重算后, 重新生成:
  - pipeline_report.html
  - taxonomy_sunburst.html
  - classification_sankey.html
  - 各类汇总 TSV
"""
import os
import subprocess
import sys
import time

HOME = os.path.expanduser("~")
REPORT = f"{HOME}/MMPV-RNA/virome_discovery_pipeline/report_pipeline.py"
ROOT = "/home/zhangwenda/data-test/out"

print(f"=== 层3 重生成报告 ===")
print(f"脚本: {REPORT}")
print(f"目标: {ROOT}")

# 记录旧产物时间戳 (对比用)
reports_dir = os.path.join(ROOT, "10_Reports")
before = {}
for fn in ["pipeline_report.html", "taxonomy_sunburst.html",
           "classification_sankey.html", "stage_summary.tsv"]:
    p = os.path.join(reports_dir, fn)
    if os.path.exists(p):
        before[fn] = (os.path.getmtime(p), os.path.getsize(p))

t0 = time.time()
cmd = [sys.executable, REPORT, "-o", ROOT]
print(f"\nCMD: {' '.join(cmd)}\n", flush=True)
r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
dt = time.time() - t0
print(f"rc={r.returncode}  耗时 {dt:.0f}s")
if r.stdout:
    print("STDOUT 尾部:")
    print("\n".join(r.stdout.strip().split("\n")[-30:]))
if r.returncode != 0 and r.stderr:
    print("STDERR 尾部:")
    print("\n".join(r.stderr.strip().split("\n")[-15:]))

print("\n=== 产物对比 ===")
for fn in ["pipeline_report.html", "taxonomy_sunburst.html",
           "classification_sankey.html", "stage_summary.tsv"]:
    p = os.path.join(reports_dir, fn)
    if os.path.exists(p):
        mt, sz = os.path.getmtime(p), os.path.getsize(p)
        old = before.get(fn)
        if old:
            chg = "已更新" if (mt > old[0] + 1) else "未变"
            print(f"  {fn:<30} {sz/1e6:.2f}MB  {chg}")
        else:
            print(f"  {fn:<30} {sz/1e6:.2f}MB  新建")
    else:
        print(f"  {fn:<30} 缺失!")
