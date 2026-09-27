#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""文本级验证: 备份 vs 现装 的差异必须是"纯新增", 即除新增行外零改动。"""
import difflib
from pathlib import Path

D = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
BAK = D / "run_host_prediction.py.bak_plantwl_20260915"
CUR = D / "run_host_prediction.py"

a = BAK.read_text(encoding="utf-8").splitlines(keepends=True)
b = CUR.read_text(encoding="utf-8").splitlines(keepends=True)
sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
added, removed, changed = [], [], 0
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == "equal":
        continue
    if tag == "insert":
        added.append((i1 + 1, j2 - j1))
    elif tag == "delete":
        removed.append((i1 + 1, i2 - i1))
    else:
        changed += 1
        print("!! replace 块: 备份 %d-%d -> 现装 %d-%d" % (i1 + 1, i2, j1 + 1, j2))

print("备份 %d 行 / 现装 %d 行" % (len(a), len(b)))
print("新增块 %d 处, 共 %d 行: %s" % (len(added), sum(n for _, n in added), added))
print("删除块 %d 处, 共 %d 行: %s" % (len(removed), sum(n for _, n in removed), removed))
print("替换块 %d 处" % changed)

print("\n新增块位置上下文:")
for ln, n in added:
    print("  在备份第 %d 行之后插入 %d 行 ; 上一行=%r" % (ln - 1, n, a[ln - 2].rstrip()[:90]))
    print("    下一行=%r" % a[ln - 1].rstrip()[:90])

# 黑名单定义区必须逐字节相同
def block(lines, start_kw, end_kw):
    s = next(i for i, l in enumerate(lines) if start_kw in l)
    e = next(i for i, l in enumerate(lines) if end_kw in l and i > s)
    return lines[s:e]

for name, s_kw, e_kw in [("科黑名单区", "NON_PLANT_FAMILIES_FALLBACK = [", "FAMILY_FIRST_VETO"),
                         ("属黑名单区", "NON_PLANT_GENERA = {", "def ") ]:
    try:
        ba, bb = block(a, s_kw, e_kw), block(b, s_kw, e_kw)
        print("\n%s: 备份 %d 行 / 现装 %d 行 -> %s"
              % (name, len(ba), len(bb), "逐字节相同" if ba == bb else "**有差异**"))
    except StopIteration:
        print("\n%s: 未定位到 (跳过)" % name)

ok = not removed and changed == 0
print("\n%s" % ("PASS 差异全部为新增行, 无删改" if ok else "FAIL 存在删除/替换"))
