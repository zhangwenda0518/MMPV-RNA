"""统一的文本 diff 工具：输出 v6.3 与 v6.4 补丁之间的全部改动块。

用法: python difftool.py 旧文件 新文件 [上下文行数]
"""
import difflib
import sys

old_p, new_p = sys.argv[1], sys.argv[2]
ctx = int(sys.argv[3]) if len(sys.argv) > 3 else 4
out_p = sys.argv[4] if len(sys.argv) > 4 else None
out = open(out_p, "w", encoding="utf-8") if out_p else None


def emit(line=""):
    if out:
        out.write(line + "\n")
    else:
        print(line)

with open(old_p, encoding="utf-8") as fh:
    old = fh.read().splitlines()
with open(new_p, encoding="utf-8") as fh:
    new = fh.read().splitlines()

sm = difflib.SequenceMatcher(None, old, new, autojunk=False)
n_add = n_del = 0
blocks = 0
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == "equal":
        continue
    blocks += 1
    n_del += i2 - i1
    n_add += j2 - j1
    emit("=" * 78)
    emit("[%s] 旧 %d-%d (%d 行)  ->  新 %d-%d (%d 行)" % (tag, i1 + 1, i2, i2 - i1, j1 + 1, j2, j2 - j1))
    lo = max(0, i1 - ctx)
    hi = min(len(old), i2 + ctx)
    for k in range(lo, i1):
        emit("  %5d| %s" % (k + 1, old[k]))
    for k in range(i1, i2):
        emit("- %5d| %s" % (k + 1, old[k]))
    for k in range(j1, j2):
        emit("+ %5d| %s" % (k + 1, new[k]))
    for k in range(i2, hi):
        emit("  %5d| %s" % (k + 1, old[k]))
emit("=" * 78)
emit("合计改动块 %d 个；删除 %d 行，新增 %d 行；旧 %d 行 -> 新 %d 行" % (blocks, n_del, n_add, len(old), len(new)))
if out:
    out.close()
