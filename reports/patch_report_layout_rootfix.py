# -*- coding: utf-8 -*-
"""
整改 MMPV 分析管线报告生成器的目录解析（根因修复，2026-10-09）
=============================================================

背景：`generate_pipeline_report.py` 虽然 import 了统一布局表 `mmpv_common/io_layout`，
但实际用法自相矛盾，且多处绕过它写死目录名：

  1. `_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))`
     —— 默认强行 legacy，而 `io_layout.DEFAULT_LAYOUT == "standard"`。
        同一套代码，管道产standard 大写目录名，报告脚本默认去找 legacy 小写名。
  2. `354`: `out / "08_similarity"`          写死
  3. `2103`: `out / _D['a_similarity']` 后再if 不到回退 `out / "08_similarity"`  写死
  4. `2406`: `_legacy_rep = out / "10_report"`  写死

后果（已实测）：standard 布局下报告脚本找不到目录；更糟的是 2406 那处会
**凭空创建 `10_report/`**（实测 152 MB / 1039 文件），且因"存在就沿用"的
回退逻辑，后续每次都往错目录里塞。

修法：不新增第二套目录表，而是复用 `io_layout` 既有的"按磁盘探测"范式
（同 `locate_centroids()` / `locate_taxonomy()`），由它统一探测三代命名：
  10阶段legacy / 9阶段 legacy / standard。
调用方只声明"我要哪个逻辑阶段"，不再关心磁盘上叫什么。

原子：写盘放最后 + 每处锚点计数断言，失败整体回滚。
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
        print(f"FAIL[{tag}] 期望 {n_expect} 处锚点，实际 {n} 处 -> 整体回滚，不写盘")
        print("---- 锚点 ----")
        print(old[:300])
        sys.exit(2)
    src = src.replace(old, new, n_expect)
    edits.append(tag)


# ── 1. import 统一布局的定位器；默认值改为遵循 DEFAULT_LAYOUT ──────────
OLD1 = '''from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))
'''
NEW1 = '''from mmpv_common.io_layout import layout_dirs as _layout_dirs
from mmpv_common.io_layout import locate_analysis_dir, locate_analysis_dirs
# 布局默认值统一交给 io_layout.DEFAULT_LAYOUT 决定（不再在此处另行默认 legacy，
# 否则与管道实际产出的目录名不一致）。需要强制 legacy 时设 MMPV_IO_LAYOUT=legacy。
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))
'''
rep(OLD1, NEW1, 1, "1:import+default")

# ── 2. 删掉自造的 _DIR_FALLBACK / _resolve_dir，改走统一探测 ──────────
#    io_layout.locate_analysis_dir() 按三代候选探测磁盘，覆盖了原表的全部回退，
#    且额外能吃下 standard 大写命名 —— 原表做不到。
OLD2 = '''_DIR_FALLBACK = {
    _D['a_post']: "06_post_analysis",
    _D['a_assembly']: "04_assembly",
    _D['a_extract']: "05_extraction",
}

def _resolve_dir(base_dir, path_segment):
    """Resolve directory path with new→old naming fallback."""
    d = Path(base_dir) / path_segment
    if d.is_dir():
        return d
    if path_segment in _DIR_FALLBACK:
        fb = Path(base_dir) / _DIR_FALLBACK[path_segment]
        if fb.is_dir():
            return fb
    return d
'''
NEW2 = '''def _resolve_dir(base_dir, stage_key, default_segment=None):
    """按逻辑阶段键解析真实目录：委托 io_layout 统一探测三代命名。

    stage_key 用逻辑名（post/assembly/extract/similarity/dvg/report...），
    不再传目录名字符串 —— 后者正是各处硬编码的根源。
    """
    return locate_analysis_dir(base_dir, stage_key, default_segment=default_segment)
'''
rep(OLD2, NEW2, 1, "2:resolve-dir")

# ── 3. Similarity 数据收集：写死的 08_similarity ──────────────────────
OLD3 = '''    vdir = _find_virus_dir(out / "08_similarity", virus_acc)
'''
NEW3 = '''    vdir = _find_virus_dir(locate_analysis_dir(out, "similarity"), virus_acc)
'''
rep(OLD3, NEW3, 1, "3:similarity-collect")

# ── 4. copy_results_to_report 里的 sim_dir 回退 ───────────────────────
OLD4 = '''    sim_dir = out / _D['a_similarity']
    if not sim_dir.is_dir():
        sim_dir = out / "08_similarity"      # 兼容 10 阶段旧产物
'''
NEW4 = '''    sim_dir = locate_analysis_dir(out, "similarity")
'''
rep(OLD4, NEW4, 1, "4:sim-dir")

# ── 5. report_dir：去掉写死的 10_report 优先（它会凭空创建目录）────────
OLD5 = '''    _legacy_rep = out / "10_report"
    report_dir = _legacy_rep if _legacy_rep.is_dir() else out / _D['a_report']
'''
NEW5 = '''    # 报告目录一律按逻辑键探测（10_report / 09_report / 09_Report 都覆盖）。
    # 注意：不能再"10_report 不存在就自己造一个" —— 那会凭空多出 152 MB。
    report_dir = locate_analysis_dir(out, "report", default_segment=_D['a_report'])
'''
rep(OLD5, NEW5, 1, "5:report-dir")

# ── 写盘（放最后） ──────────────────────────────────────────────────
tmp = TARGET + ".tmp_fix5"
with open(tmp, "w", encoding="utf-8") as f:
    f.write(src)
os.replace(tmp, TARGET)

print("PATCHED OK ->", TARGET)
print("edits:", ", ".join(edits))
print("size:", len(orig), "->", len(src), f"(delta {len(src)-len(orig)})")