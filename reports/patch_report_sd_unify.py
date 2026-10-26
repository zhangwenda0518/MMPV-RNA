# -*- coding: utf-8 -*-
"""
整改第二批：消除 generate_pipeline_report.py 全部 23 处 `_D['a_*']` 目录名拼接。
=================================================================================

起因（我自己引入的回归）
------------------------
第一批整改把 `_D` 的默认值从 legacy 改成遵循 `io_layout.DEFAULT_LAYOUT`(standard)，
本意是消除"管道产standard 大写、报告脚本默认找 legacy 小写"的矛盾。但该文件有
**23 处** `out / _D['a_filter']` 这类**目录名字符串拼接**，分布在 6 个函数里：

    1215_collect_extract_stats        2006 ensure_summary_plots
    2010/2017/2018                   2020/2024
    2061..2076 copy_results_to_report 2079post_dir
    2107 dvg_plots                    2349 main() summary 定位
    2374..2376                        2380/2386                2400 report_dir

默认值一改，这 23 处全部从 `02_filtering` 变成 `02_Filtering`，在 legacy 产物上
全部落空 → `No summary file found` 直接退出（实测已复现）。

本批修法（零回归风险）
----------------------
**不改 `_D` 的默认值**（仍显式 legacy），只把「拼目录名字符串」换成「按磁盘探测」：
新增 `_sd(root, key)` 辅助函数，内部调`io_layout.locate_analysis_dir()`，
再把 23 处 `out / _D['a_x']` / `os.path.join(_D['a_x'], ...)` 统一替换。

这样同时达成两件事，且都不破坏既有行为：
  ① standard 大写布局能正确解析（原先完全失效）
  ② legacy 小写/ 10 阶段旧命名行为逐字节不变（探测会先命中 legacy 候选）

注意 `_D` 仍保留：它承载非 `a_*` 的键、以及「阶段未跑时该用哪个默认名」，
两者都需要布局表。

原子：写盘放最后 + 每处锚点计数断言，失败整体回滚。
"""
import os
import re
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

# ── 0. 回滚第一批对 _D 默认值的改动（恢复 legacy 默认，零回归）──────────
OLD0 = '''from mmpv_common.io_layout import locate_analysis_dir, locate_analysis_dirs
# 布局默认值统一交给 io_layout.DEFAULT_LAYOUT 决定（不再在此处另行默认 legacy，
# 否则与管道实际产出的目录名不一致）。需要强制 legacy 时设 MMPV_IO_LAYOUT=legacy。
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))
'''
NEW0 = '''from mmpv_common.io_layout import locate_analysis_dir, locate_analysis_dirs
# 注意: 这里仍显式默认 legacy（存量项目多为 legacy 小写命名）。
# 目录名的实际解析不依赖这个默认值 —— 全部走 _sd() 按磁盘探测（见下），
# 因此 standard 大写 / 10 阶段旧命名同样能解析。
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))


def _sd(root, stage_key):
    """按逻辑阶段键解析真实目录（委托 io_layout 探测三代命名）。

    逻辑键: detect/filter/variants/post/assembly/extract/similarity/dvg/report
    替代 `out / _D['a_filter']` —— 后者把**布局名**当目录名用，一旦布局默认值
    变化或产物是另一代命名就会静默落空（本文件曾因此报 No summary file found）。
    """
    return locate_analysis_dir(root, stage_key)
'''
if src.count(OLD0) != 1:
    print("FAIL[0:rollback] 未找到第一批的 _D 改动块 -> 整体回滚")
    sys.exit(2)
src = src.replace(OLD0, NEW0, 1)
edits.append("0:rollback-_D-default+add-_sd")

# ── 1~N. 把所有 _D['a_*'] 目录名拼接换成 _sd() ---------------------------
# 替换规则（按顺序应用，逐一记录命中数）：
#   _resolve_dir(<root>, _D['a_x'])   -> _sd(<root>, 'logical')
#   <root> / _D['a_x'] / <tail>       -> _sd(<root>, 'logical') / <tail>   (链式)
#   <root>/_D['a_x']                  -> _sd(<root>, 'logical')  (无空格)
#   os.path.join(_D['a_x'], <rest>)   -> str(_sd(<root>, 'logical') / <rest>)
#   _d(_D['a_x'])                     -> _sd(out, 'logical')      (内部小helper)
#
# <root> 有多种写法: 函数内是 out / out_dir, main() 里是 _P(args.dir)。
# **不能一律写成 out** —— 否则 main() 里会引用未定义变量。

TABLE = {                       # 布局键 -> 逻辑阶段键
    "a_detect": "detect",
    "a_filter": "filter",
    "a_variants": "variants",
    "a_post": "post",
    "a_assembly": "assembly",
    "a_extract": "extract",
    "a_similarity": "similarity",
    "a_dvg": "dvg",
    "a_report": "report",
}
ROOTS = ("_P(args.dir)", "out_dir", "out")


def _record(tag, rv, lk, n):
    if n:
        edits.append(f"{tag}:{rv}:{lk} x{n}")


# ══════════════════════════════════════════════════════════════════
# 逐行、作用域感知的转换
# ══════════════════════════════════════════════════════════════════
# 为什么不用全局正则：根变量 (out / out_dir / _P(args.dir)) 是**函数局部**的。
# 全局正则猜根变量会张冠李戴 —— 实测曾把 _collect_detection_summary() 里的
# `out` 换成 main() 的 _P(args.dir)，编译能过但运行 NameError。
# 这里改为逐行追踪 `def` 边界来判定当前函数，再据此选根变量。

REV = {v: k for k, v in TABLE.items()}      # logical -> layout key


def _root_for(fn_name):
    """按所在函数决定根表达式。"""
    if fn_name == "main":
        return "_P(args.dir)"
    return "out"


def _conv_line(line, fn_name):
    """把一行里的 _D['a_x'] 目录名拼接换成 _sd(root, logical)。返回 (新行, 是否改)。"""
    if "_D['a_" not in line:
        return line, False
    root = _root_for(fn_name)
    orig = line

    # (a) 已是 locate_analysis_dir(..., default_segment=_D['a_x']) 的正当用法
    line = re.sub(r"default_segment=_D\['a_\w+'\]", "@@KEEP@@", line)

    # (b) <root> / os.path.join(_D['a_x'], <rest>)  —— 前缀属于整个表达式
    for logical, lk in REV.items():
        line = re.sub(
            re.escape(root) + r" / os\.path\.join\(_D\['" + lk + r"'\],\s*(.+?)\)",
            lambda m: f"_sd({root}, '{logical}') / {m.group(1)}", line)

    # (c) os.path.join(_D['a_x'], <rest>) —— 无前缀（用作 str 路径）
    for logical, lk in REV.items():
        line = re.sub(
            r"os\.path\.join\(_D\['" + lk + r"'\],\s*(.+?)\)",
            lambda m: f"_sd({root}, '{logical}') / {m.group(1)}", line)

    # (d) _d(_D['a_x']) / _resolve_dir(<root>, _D['a_x'])
    for logical, lk in REV.items():
        line = line.replace(f"_d(_D['{lk}'])", f"_sd({root}, '{logical}')")
        line = re.sub(r"_resolve_dir\(\s*[^,]+,\s*_D\['" + lk + r"'\]\s*\)",
                      f"_sd({root}, '{logical}')", line)

    # (e) <root> / _D['a_x']  与 <root>/_D['a_x']（含链式尾巴）
    for logical, lk in REV.items():
        line = re.sub(
            re.escape(root) + r" ?/ ?_D\['" + lk + r"'\]((?: / [^\n,)]*)?)",
            lambda m: f"_sd({root}, '{logical}'){m.group(1) or ''}", line)

    # (f) 字典值 / 列表元素：NN:_D['a_x']
    for logical, lk in REV.items():
        line = re.sub(r"(\d+:)_D\['" + lk + r"'\]",
                      lambda m: f"{m.group(1)}_sd({root}, '{logical}').name", line)

    # (g) 裸 `_D['a_x']`（列表/元组元素，无 <root> / 前缀）
    #     用 .name 保持"目录名字符串"语义，与 (f) 一致。
    for logical, lk in REV.items():
        line = re.sub(r"_D\['" + lk + r"'\]",
                      lambda m: f"_sd({root}, '{logical}').name", line)

    line = line.replace("@@KEEP@@", "default_segment=_D['a_report']")
    return line, (line != orig)


_lines = src.split("\n")
_fn = "<module>"
for _i, _ln in enumerate(_lines):
    _m = re.match(r"\s*def\s+(\w+)\s*\(", _ln)
    if _m:
        _fn = _m.group(1)
    _new, _changed = _conv_line(_ln, _fn)
    if _changed:
        _lines[_i] = _new
        edits.append(f"conv:L{_i+1} in {_fn}()")
src = "\n".join(_lines)

# ── 特例 1 已被规则 (g) 覆盖 ─────────────────────────────────────────
# 早先需为 `stage_rel == _D['a_variants']` 单独改写；规则 (g)（裸 _D['a_x']）
# 已能把它变成 `_sd(<root>,'variants').name`，两侧都指向同一解析结果，
# 故不再需要特例。若将来该行回归为裸串，规则 (g) 会再次接住。
if "if stage_rel == _D['a_variants']:" in src:
    print("FAIL[special1] 规则 (g) 未覆盖 stage_rel 比较行")
    sys.exit(2)

# ── 作用域校验：args.dir 只应出现在 main() 内 ───────────────────────
# 兜底防线：逐行正则猜错根变量时, 编译能过但运行 NameError。
cur_fn = "<module>"
for ln, line in enumerate(src.splitlines(), 1):
    m = re.match(r"\s*def\s+(\w+)\s*\(", line)
    if m:
        cur_fn = m.group(1)
    if "args.dir" in line and cur_fn != "main":
        print(f"FAIL[scope] L{ln} 在函数 {cur_fn}() 内引用了 args.dir -> 根变量推断有误")
        print("   ", line.strip()[:120])
        sys.exit(2)

# ── 收尾校验：不应再有 a_* 目录名拼接 ─────────────────────────────────
# 唯一合法例外: locate_analysis_dir(..., default_segment=_D['a_x'])
# —— 那是"该阶段未跑时用哪个默认名"，是布局表的正当用途，不是路径拼接。
ALLOWED = {"_D['a_report']"}
bad = [m for m in re.findall(r"_D\['a_\w+'\]", src) if m not in ALLOWED]
if bad:
    print(f"FAIL[final] 仍有未转换的目录名拼接: {sorted(set(bad))}")
    for ln, line in enumerate(src.splitlines(), 1):
        found = re.findall(r"_D\['a_\w+'\]", line)
        if found and not all(m in ALLOWED for m in found):
            print(f"   L{ln}: {line.strip()[:110]}")
    sys.exit(2)

# ── 写盘（放最后） ──────────────────────────────────────────────────
tmp = TARGET + ".tmp_fix7"
with open(tmp, "w", encoding="utf-8") as f:
    f.write(src)
os.replace(tmp, TARGET)

print("PATCHED OK ->", TARGET)
print("size:", len(orig), "->", len(src), f"(delta {len(src)-len(orig)})")
print("edits (%d):" % len(edits))
for e in edits:
    print("   ", e)