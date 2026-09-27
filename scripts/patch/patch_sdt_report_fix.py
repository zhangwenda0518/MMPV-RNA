#!/usr/bin/env python3
"""
修 9b 后段的两个真 bug:

P5（SDT 全部失败：name '_sp' is not defined）
  _run_sdt_matrices() 是独立方法, 内部只 import 了 csv as _c,
  但用了 _sp.run(...)。_sp 只在 run_analysis_verify() 里 import,
  Python 局部 import 不跨方法 → 每个属都抛 NameError, 0 个属出矩阵。
  → 在 _run_sdt_matrices 内补 import subprocess as _sp。

P6（report 阶段 600s 超时）
  run_reports() 用 subprocess.run(..., timeout=600) 调 report_pipeline.py,
  该脚本遍历全流程产物 (含 8k+ 目录), 600s 不够 → 超时, 报告未生成。
  → 超时提升到 3600s, 并记录耗时。
"""
import shutil
import sys
from pathlib import Path

PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virome_pipeline.py")
BAK = PIPE.with_suffix(".py.bak_sdtfix_20260910")

if not PIPE.is_file():
    sys.exit(f"找不到: {PIPE}")

src = PIPE.read_text(encoding="utf-8")
orig = src
applied = []

# ── P5 ──
old5 = """        import csv as _c
        sdt_script = SCRIPT_DIR.parent / 'virome_analysis_pipeline' / 'sdt_genus_matrix.py'"""
new5 = """        import csv as _c
        import subprocess as _sp  # 独立方法作用域, 必须自行导入 (曾缺此报 NameError)
        sdt_script = SCRIPT_DIR.parent / 'virome_analysis_pipeline' / 'sdt_genus_matrix.py'"""
if old5 in src:
    src = src.replace(old5, new5, 1)
    applied.append("P5: _run_sdt_matrices 补 import subprocess as _sp")
elif "独立方法作用域" in src:
    applied.append("P5: 已修复")
else:
    sys.exit("P5 锚点未匹配")

# ── P6 ──
old6 = """        cmd = [sys.executable, str(report_script), "-o", str(self.d['root'])]
        self.log.info("  → %s", " ".join(cmd))
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)"""
new6 = """        cmd = [sys.executable, str(report_script), "-o", str(self.d['root'])]
        self.log.info("  → %s", " ".join(cmd))
        try:
            import time as _t
            _t0 = _t.time()
            # 报告需遍历全流程产物 (含数千目录), 600s 不够 → 3600s
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            self.log.info("  report_pipeline 耗时 %.1fs (rc=%d)", _t.time() - _t0, result.returncode)"""
if old6 in src:
    src = src.replace(old6, new6, 1)
    applied.append("P6: report_pipeline 超时 600→3600s + 耗时日志")
elif "timeout=3600" in src and "report_pipeline 耗时" in src:
    applied.append("P6: 已修复")
else:
    sys.exit("P6 锚点未匹配")

if src != orig:
    if not BAK.exists():
        shutil.copy2(PIPE, BAK)
        print(f"备份: {BAK}")
    PIPE.write_text(src, encoding="utf-8")
    print("[virome_pipeline.py]")
    for a in applied:
        print("  -", a)
else:
    print("无改动")
