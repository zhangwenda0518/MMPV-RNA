#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 /tmp/plant_whitelist_block.py 插入 run_host_prediction.py (锚点=科级否决开关注释行之前)。

幂等: 若文件中已存在 'PLANT_FAMILIES_WHITELIST', 直接跳过。
用法: python3 patch_insert_whitelist.py [--check]
"""
import sys
from pathlib import Path

TARGET = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py")
BLOCK = Path("/tmp/plant_whitelist_block.py")
ANCHOR = "# 科级否决开关: True=科在非植物科名单即否决 Plant (科级优先);"

text = TARGET.read_text(encoding="utf-8")
if "PLANT_FAMILIES_WHITELIST" in text:
    print("[SKIP] 已存在白名单常量, 未改动")
    sys.exit(0)
if ANCHOR not in text:
    print("[FAIL] 锚点未找到: %s" % ANCHOR)
    sys.exit(1)

block = BLOCK.read_text(encoding="utf-8").rstrip("\n")
new = text.replace(ANCHOR, block + "\n\n" + ANCHOR, 1)
if "--check" in sys.argv:
    print("[CHECK] 将插入 %d 行, 文件 %d -> %d 行" % (block.count("\n") + 1,
                                                    text.count("\n") + 1, new.count("\n") + 1))
    sys.exit(0)
TARGET.write_text(new, encoding="utf-8")
print("[OK] 已插入 %d 行, 文件 %d -> %d 行" % (block.count("\n") + 1,
                                              text.count("\n") + 1, new.count("\n") + 1))
