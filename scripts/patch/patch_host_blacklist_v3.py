"""黑名单补齐 v3: 加入"植物科下的非植物属"漏网项 (单行 list 格式)。

依据: 权威宿主概率表交叉验证
  Plant 科(35个) → Pred_Host != Plant 的属 → 实际命中/零命中

新增属:
  Betapartitivirus  (Partitiviridae, Fungi, P(plant)=0.208, 实测 21 条)
  Botrexvirus       (Alphaflexiviridae, Fungi)
  Sclerodarnavirus  (Alphaflexiviridae, Fungi)
  Draflysatellite   (Alphasatellitidae, Insecta)
  Unirnavirus       (Amalgaviridae, Fungi)
  Zybavirus         (Amalgaviridae, Fungi)
  Betaendornavirus  (Endornaviridae, Fungi)
  Cryspovirus       (Partitiviridae, Protist)
  Gammapartitivirus (Partitiviridae, Fungi)
"""
import ast
import re
import shutil
import sys
import py_compile

TARGET = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
BAK = TARGET + ".bak_blacklist4_20260902"

NEW_GENERA = [
    "Betapartitivirus", "Botrexvirus", "Sclerodarnavirus", "Draflysatellite",
    "Unirnavirus", "Zybavirus", "Betaendornavirus", "Cryspovirus",
    "Gammapartitivirus",
]

src = open(TARGET, encoding="utf-8").read()
lines = src.split("\n")

target_idx = None
for i, l in enumerate(lines):
    if l.startswith("NON_PLANT_GENERA = ["):
        target_idx = i
        break
if target_idx is None:
    print("FAIL: 未找到 NON_PLANT_GENERA 单行定义")
    sys.exit(1)

line = lines[target_idx]
existing = ast.literal_eval(line[len("NON_PLANT_GENERA = "):])
print(f"现有属: {len(existing)}")

to_add = [g for g in NEW_GENERA if g not in existing]
print(f"待新增: {len(to_add)} -> {to_add}")
if not to_add:
    print("无需修改")
    sys.exit(0)

merged = existing + to_add
# 按原格式单行重建
lines[target_idx] = "NON_PLANT_GENERA = [" + ", ".join(f"'{g}'" for g in merged) + "]"

shutil.copy2(TARGET, BAK)
print(f"备份: {BAK}")
open(TARGET, "w", encoding="utf-8").write("\n".join(lines))

try:
    py_compile.compile(TARGET, doraise=True)
    print("py_compile: PASS")
except py_compile.PyCompileError as e:
    print(f"py_compile FAIL: {e}")
    shutil.copy2(BAK, TARGET)
    print("已回滚")
    sys.exit(1)

t = ast.parse(open(TARGET, encoding="utf-8").read())
for n in t.body:
    if isinstance(n, ast.Assign):
        for tg in n.targets:
            if isinstance(tg, ast.Name) and tg.id == "NON_PLANT_GENERA":
                v = ast.literal_eval(n.value)
                print(f"复核: NON_PLANT_GENERA 现有 {len(v)} 个")
                miss = [g for g in NEW_GENERA if g not in v]
                print(f"新增校验: {'全部到位' if not miss else '缺失 ' + str(miss)}")
