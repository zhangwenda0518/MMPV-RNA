"""黑名单机制侦察: 只读, 不写任何结果.
目标:
  1. 服务器 run_host_prediction.py 的黑名单三件套 + 判定链 (RVH/PB2 映射)
  2. C7 概率表 (family/genus) 的列与原生/藻条目
  3. 两个结果树的 C9 分类表列名与行数
"""
import ast
import csv
import os
import re
from collections import Counter, defaultdict

RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
SRC = f"{RUN}/run_host_prediction.py"

src = open(SRC, encoding="utf-8", errors="replace").read()
print("=" * 78)
print(f"[1] run_host_prediction.py: {len(src.splitlines())} 行, {len(src)} 字节")
print("=" * 78)

tree = ast.parse(src)
names = []
for n in tree.body:
    if isinstance(n, ast.Assign):
        for t in n.targets:
            if isinstance(t, ast.Name):
                names.append(n.targets[0].id if isinstance(n.targets[0], ast.Name) else "?")
    elif isinstance(n, ast.FunctionDef):
        names.append(n.name + "()")
print("顶层名字:", ", ".join(names))

sets = {}
for n in tree.body:
    if isinstance(n, ast.Assign):
        for t in n.targets:
            if isinstance(t, ast.Name) and t.id in (
                    "NON_PLANT_GENERA", "NON_PLANT_FAMILIES_FALLBACK",
                    "TRUSTED_LEVELS", "RVH_MAP", "PB2_MAP", "HOST_ALIASES"):
                try:
                    sets[t.id] = ast.literal_eval(n.value)
                except Exception as e:
                    print(f"  [literal_eval 失败] {t.id}: {e}")

for k, v in sets.items():
    print(f"\n--- {k} ({len(v)} 项) ---")
    if isinstance(v, dict):
        for kk, vv in list(v.items()):
            print(f"    {kk!r} -> {vv!r}")
    else:
        print("   ", sorted(v))

# 追踪 h_ictv / rvh / pb2 相关的判定行
print("\n" + "=" * 78)
print("[2] 判定链关键行 (含 rvh / pb2 / ictv / blacklist 的行)")
print("=" * 78)
for i, line in enumerate(src.splitlines(), 1):
    low = line.lower()
    if any(k in low for k in ("viridiplantae", "parse_rvh", "pb2", "is_blacklisted",
                              "h_ictv", "final_host", "def decision", "def _resolve")):
        print(f"{i:5d}| {line.rstrip()}")

# C7 概率表
print("\n" + "=" * 78)
print("[3] C7 概率表")
print("=" * 78)
C7 = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
if os.path.isdir(C7):
    for fn in sorted(os.listdir(C7)):
        p = os.path.join(C7, fn)
        if not os.path.isfile(p):
            continue
        with open(p, errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
        print(f"\n--- {fn} ({os.path.getsize(p)}B) ---")
        print("  列:", " | ".join(f"{i+1}:{h}" for i, h in enumerate(hdr)))
else:
    print(f"  目录不存在: {C7}")

# C9 表列名 + Plant 行数
print("\n" + "=" * 78)
print("[4] 两棵树的 C9 分类表")
print("=" * 78)
TREES = {
    "goji-Lycium": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}
for label, root in TREES.items():
    c9 = f"{root}/06_HostPrediction/C9_ICTV_result"
    print(f"\n### {label}: {c9}")
    if not os.path.isdir(c9):
        print("  目录不存在")
        continue
    for fn in sorted(os.listdir(c9)):
        print("   -", fn, os.path.getsize(os.path.join(c9, fn)))
    f = f"{c9}/classification_result.tsv"
    if os.path.isfile(f):
        with open(f, errors="replace") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            hdr = rd.fieldnames
            cnt = Counter()
            lvl = Counter()
            for row in rd:
                cnt[(row.get("Predicted_Host") or "").strip()] += 1
                lvl[(row.get("Determination_Level") or "").strip().split("(")[0]] += 1
        print("  列:", " | ".join(hdr))
        print("  Predicted_Host:", dict(cnt.most_common()))
        print("  Determination_Level 主前缀:", dict(lvl.most_common()))
