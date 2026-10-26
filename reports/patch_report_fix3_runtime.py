# -*- coding: utf-8 -*-
"""
修复第三批遗留的两处运行期错误（端到端回归暴露）
=================================================

报错: `UnboundLocalError: local variable '_P' referenced before assignment`
  at main(): `for loc in [out / _sd(_P(args.dir), 'filter').name / p, ...]`

两个独立错误:
  1. **根变量选错**: main() 开头已有 `out = Path(args.dir)`，根应当是 `out`。
     我却写成了 `_P(args.dir)` —— 而 `_P` 只在下面的 `if args.ai_key:` 分支里
     才`from pathlib import Path as _P` 导入（L2381），在它之前引用即未绑定。
     → 修法: `_P(args.dir)` -> `out`，并把 `_P` 的导入提到 main() 开头（顺带修好
       本就存在的"局部导入"隐患：它在ai 分支内，却服务后续多处使用）。
  2. **多加了 .name**: 产出是 `out / _sd(...).name / p` —— 目录名后面又拼 `/ p`，
     语义变成"拿目录的名字字符串当目录用"。规则 (f)/(g) 的 `.name` 只适用于
     **字典值/字符串比较**场景（需要目录名字符串），不适用于路径拼接。
     → 修法: 路径拼接处去掉 `.name`。

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


# ── 1. main() 的 summary 定位：根用 out，去掉 .name ─────────────────────
OLD1 = """        for loc in [out / _sd(_P(args.dir), 'filter').name / p, out / _sd(_P(args.dir), 'detect').name / "summary" / p]:"""
NEW1 = """        for loc in [_sd(out, 'filter') / p, _sd(out, 'detect') / "summary" / p]:"""
rep(OLD1, NEW1, 1, "1:main-summary-root")

# ── 2. main() 里其余 _P(args.dir) 用法统一改成 out ─────────────────────
#    （_P 仅在 ai 分支里局部导入，服务于这些行必然踩未绑定）
import re
n_left = 0
while "_P(args.dir)" in src:
    src = src.replace("_P(args.dir)", "out", 1)
    n_left += 1
if n_left:
    edits.append(f"2:_P(args.dir)->out x{n_left}")

# ── 3. _P 的局部导入移到模块顶部（它服务多处，不是 ai 专属）────────────
OLD3 = """from mmpv_common.io_layout import layout_dirs as _layout_dirs"""
NEW3 = """from pathlib import Path
_P = Path          # 兼容旧代码里的 `_P(args.dir)` 写法(main 内曾依赖局部导入)

from mmpv_common.io_layout import layout_dirs as _layout_dirs"""
rep(OLD3, NEW3, 1, "3:module-level-_P")

# ── 4. 路径拼接处误加的 .name：`_sd(root, 'x').name / rest` -> `_sd(root, 'x') / rest`
before = src
src = re.sub(r"(_sd\([^()]*?\))\s*\.name(\s*/\s*)", r"\1\2", src)
if src != before:
    edits.append("4:drop-.name-before-path")

# ── 5. 列表内类型一致性：`_sd(x,'k').name` 与相邻 Path 元素混用 ────────
#    L2024: `for d in [_sd(out,'filter').name, _sd(out,'detect') / 'summary']:`
#    左边目录名字符串、右边 Path，混用会让后续 Path()/str() 行为不一致。
#    统一成 Path 语义。
OLD5 = """            for d in [_sd(out, 'filter').name,_sd(out, 'detect') / 'summary']:"""
NEW5 = """            for d in [_sd(out, 'filter'), _sd(out, 'detect') / 'summary']:"""
if src.count(OLD5) == 1:
    src = src.replace(OLD5, NEW5, 1)
    edits.append("5:list-type-consistency")
elif ".name" in src and "for d in [_sd" in src:
    print("FAIL[5] 列表混用行形态与预期不符，需人工确认")
    for ln, line in enumerate(src.splitlines(), 1):
        if "for d in [_sd" in line:
            print(f"   L{ln}: {line.strip()[:120]}")
    sys.exit(2)

# ── 6. main() 内ai 分支的局部 `from pathlib import Path as _P` 已冗余 ──
OLD6 = """        # Quick-load data needed for AI summary
        from pathlib import Path as _P
"""
if src.count(OLD6) == 1:
    src = src.replace(OLD6, """        # Quick-load data needed for AI summary
""", 1)
    edits.append("6:drop-redundant-local-_P-import")

# ── 收尾：不应再有 _sd(...).name / 形态 ──────────────────────────────
bad = re.findall(r"_sd\([^()]*?\)\s*\.name\s*/", src)
if bad:
    print(f"FAIL[final] 路径拼接处仍带 .name: {bad[:3]}")
    sys.exit(2)

tmp = TARGET + ".tmp_fix8"
with open(tmp, "w", encoding="utf-8") as f:
    f.write(src)
os.replace(tmp, TARGET)

print("PATCHED OK ->", TARGET)
print("size:", len(orig), "->", len(src), f"(delta {len(src)-len(orig)})")
for e in edits:
    print("   ", e)