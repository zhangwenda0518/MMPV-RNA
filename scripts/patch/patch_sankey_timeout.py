#!/usr/bin/env python3
"""
修 report 的 Sankey 图超时 (classification_sankey.html 未生成):

report_pipeline.py 调用 stats/taxonomic_sankey.py 时 timeout=120s (两处),
但输入 final_integrated_classification.tsv 是 179MB, 解析需要更久 → 超时。

→ 两处 timeout=120 提到 1800s (仅限 sankey 调用, 不动 assembly/ident/cobra 的 60/120s)。
"""
import shutil
import sys
from pathlib import Path

RP = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/report_pipeline.py")
BAK = RP.with_suffix(".py.bak_sankey_20260910")

if not RP.is_file():
    sys.exit(f"找不到: {RP}")

lines = RP.read_text(encoding="utf-8").split("\n")
if not BAK.exists():
    shutil.copy2(RP, BAK)
    print(f"备份: {BAK}")

changed = 0
for i in (1013, 1035):  # 0-based: 1014行 / 1036行
    if i < len(lines) and "timeout=120" in lines[i]:
        lines[i] = lines[i].replace("timeout=120", "timeout=1800")
        print(f"  行{i+1}: timeout=120 → 1800")
        changed += 1
    else:
        print(f"  行{i+1}: 未匹配 (内容: {lines[i][:90] if i < len(lines) else 'EOF'})")

if changed:
    RP.write_text("\n".join(lines), encoding="utf-8")
    print(f"[report_pipeline.py] 改了 {changed} 处")
else:
    print("无改动")
