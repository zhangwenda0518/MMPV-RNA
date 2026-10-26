# -*- coding: utf-8 -*-
"""
补修健壮性缺陷：阶段目录不存在时 collect_charts() 崩溃
=======================================================

端到端验证 standard 布局时暴露：
    FileNotFoundError: [Errno 2] No such file or directory: '/tmp/std_mini/04_PostAnalysis'
      at collect_charts() -> `for d in Path(post_dir).iterdir():`

这条路径一直存在，只是以前目录名总能被找到，从未触发。目录解析改为
"按磁盘探测"后，未跑该阶段的场景会真正走到这里 —— 于是暴露。

**这是独立的健壮性缺陷，不是本次整改引入的**：语义上「该阶段没跑」
应当是"跳过该段并继续出报告"，而不是整个报告崩掉。

修法（两处，最小改动）:
  1. collect_charts(): iterdir() 前加 is_dir() 守卫
  2. collect_dvg_charts(): 入口加同样的守卫（对称，防止同类问题）

原子：写盘最后 + 锚点断言。
"""
import os
import sys

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_REPO   = f"{_SERVER}/MMPV-RNA"


TARGET = f"{_REPO}/virome_analysis_pipeline/generate_pipeline_report.py"

with open(TARGET, "r", encoding="utf-8") as f:
    src = f.read()

orig = src
edits = []


def rep(old, new, n_expect, tag):
    global src
    n = src.count(old)
    if n != n_expect:
        print(f"FAIL[{tag}] 期望 {n_expect} 处，实际 {n} 处 -> 整体回滚")
        print("---- 锚点 ----")
        print(old[:300])
        sys.exit(2)
    src = src.replace(old, new, n_expect)
    edits.append(tag)


# ── 1. collect_charts: iterdir 前守卫 ─────────────────────────────────
OLD1 = '''    charts = {}
    vdir = Path(post_dir) / virus_acc
    if not vdir.exists():
        for d in Path(post_dir).iterdir():
'''
NEW1 = '''    charts = {}
    # 该阶段未跑时目录可能不存在 -> 降级为"无图"，而非抛 FileNotFoundError。
    if not Path(post_dir).is_dir():
        return charts
    vdir = Path(post_dir) / virus_acc
    if not vdir.exists():
        for d in Path(post_dir).iterdir():
'''
rep(OLD1, NEW1, 1, "1:collect_charts-guard")

# ── 2. collect_dvg_charts 已有守卫，无需改动 ────────────────────────
# 实测确认: `if not plots_dir.exists(): return charts`（L186）已覆盖。
# 当初想对称补一个入口守卫属于过度设计 —— 撤销，只保留必要的那一处。
CHECK_DVG = '''def collect_dvg_charts(dvg_dir, virus_acc):
    """Collect DVG charts from Stage 9 output."""
    charts = {}
    plots_dir = Path(dvg_dir) / "Summary_Analysis_Report" / "Virus_Specific_Plots"
    if not plots_dir.exists(): return charts'''
if src.count(CHECK_DVG) != 1:
    print("FAIL[2] collect_dvg_charts 的既有守卫形态与预期不符，需人工确认")
    sys.exit(2)
edits.append("2:collect_dvg_charts-already-guarded (no change)")

tmp = TARGET + ".tmp_fix9"
with open(tmp, "w", encoding="utf-8") as f:
    f.write(src)
os.replace(tmp, TARGET)

print("PATCHED OK ->", TARGET)
print("size:", len(orig), "->", len(src), f"(delta {len(src)-len(orig)})")
for e in edits:
    print("   ", e)