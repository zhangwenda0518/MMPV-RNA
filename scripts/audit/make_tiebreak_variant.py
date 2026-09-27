#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 TIE_BREAK_ORDER 变体补丁，用于「平票兜底顺序」敏感性测试。

用法: python3 make_tiebreak_variant.py <源补丁> <输出补丁> "T1,T2,T3,..."
只替换唯一的 `TIE_BREAK_ORDER <- c(...)` 定义行，其余字节不动。
"""
import re
import sys

if len(sys.argv) < 4:
    print(__doc__)
    raise SystemExit(2)

src, out, order_arg = sys.argv[1], sys.argv[2], sys.argv[3]
order = [t.strip() for t in order_arg.split(",") if t.strip()]
if len(order) != 7:
    print("需要 7 个工具名，实际 %d" % len(order))
    raise SystemExit(2)

txt = open(src, encoding="utf-8", errors="replace").read()
new_line = "TIE_BREAK_ORDER <- c(%s)" % ", ".join('"%s"' % t for t in order)
txt2, n = re.subn(r"^TIE_BREAK_ORDER <- c\(.*\)$", new_line, txt, count=1, flags=re.M)
if n != 1:
    print("未找到 TIE_BREAK_ORDER 定义行 (n=%d)" % n)
    raise SystemExit(3)

with open(out, "w", encoding="utf-8", newline="\n") as f:
    f.write(txt2)
print("源: %s" % src)
print("出: %s" % out)
print("新: %s" % new_line)
