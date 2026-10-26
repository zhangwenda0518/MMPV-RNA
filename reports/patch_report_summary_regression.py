# -*- coding: utf-8 -*-
"""
跟进修复：main() 的 summary 定位仍直接拼 `_D['a_filter']` 目录名字符串。

背景：上一批整改把 `_D` 的默认值从legacy 改成遵循 `io_layout.DEFAULT_LAYOUT`
(=standard)，于是 `_D['a_filter']` 从 `02_filtering` 变成 `02_Filtering`。
barbarum 产物是 legacy 小写命名，`(out / _D['a_filter'] / p).is_file()` 全部落空
→ `No summary file found` 直接退出。**这是上一批引入的回归。**

修法：main() 也改走 `locate_analysis_dir()`（按磁盘实际探测），不再拼目录名字符串。
这既修回归，又把main() 从"依赖布局名"变成"依赖磁盘事实"，与前三处一致。

原子：写盘放最后 + 锚点计数断言。
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


OLD = '''    summary_in = None
    for p in ["high_conf.summary.tsv", "all_viruses.best.summary.tsv", "all_viruses.summary.tsv"]:
        for loc in [out / _D['a_filter'] / p, out / _D['a_detect'] / "summary" / p]:
            if loc.is_file(): summary_in = loc; break
        if summary_in: break
'''
NEW = '''    summary_in = None
    # 按磁盘实际目录定位（不拼 _D['a_filter'] 这类目录名字符串）：
    # 产物可能是 standard 大写 / 9 阶段 legacy 小写 / 10 阶段旧命名，探测才可靠。
    _filter_dir = locate_analysis_dir(out, "filter")
    _detect_dir = locate_analysis_dir(out, "detect")
    for p in ["high_conf.summary.tsv", "all_viruses.best.summary.tsv", "all_viruses.summary.tsv"]:
        for loc in [_filter_dir / p, _detect_dir / "summary" / p]:
            if loc.is_file(): summary_in = loc; break
        if summary_in: break
'''
rep(OLD, NEW, 1, "main:summary-lookup")

# ensure_summary_plots / _collect_* 系列也用 _D 拼目录名，一并排查是否还有同类硬依赖
# 只做报告，不改动，先列出调用点供人工确认。

tmp = TARGET + ".tmp_fix6"
with open(tmp, "w", encoding="utf-8") as f:
    f.write(src)
os.replace(tmp, TARGET)

print("PATCHED OK ->", TARGET)
print("edits:", ", ".join(edits))
print("size:", len(orig), "->", len(src), f"(delta {len(src)-len(orig)})")