#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""补丁2: 在 run_host_prediction.py 中加 WHITELIST_OVERRIDES_FAMILY_VETO 开关 (默认 False, 行为不变)。

锚点 A: FAMILY_FIRST_VETO = True 之后插入常量定义
锚点 B: is_blacklisted 里科级否决分支加白名单放行分支
幂等: 已存在 WHITELIST_OVERRIDES_FAMILY_VETO 则跳过。
"""
import sys
from pathlib import Path

TARGET = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py")
A = "FAMILY_FIRST_VETO = True\n"
A_ADD = """FAMILY_FIRST_VETO = True

# 白名单越权开关 (默认 False, 保持 2026-09-15 现状行为, 需显式决策后再开):
#   True  = 科命中非植物科名单时, 若判定层为属/种级 且 属在 PLANT_GENERA_WHITELIST,
#           则放行 Plant (即 "植物库认证的属级证据 胜过 科级否决")。
#   实测 (2026-09-15, 运行期替换 is_blacklisted 对拍):
#       goji  Plant 1029 -> 1042  (+13)
#       onekp Plant 19560 -> 19967 (+407)
#   被救回行都是 C9 依植物属给出 Plant 的嵌合行 (科=非植物科, 属=植物病毒属), 例:
#       Mimiviridae+Potyvirus、Orthoherpesviridae+Macluravirus、Barnaviridae+Sobemovirus、
#       Pithoviridae+Potexvirus、Nimaviridae+Allexivirus
#   代价: 若某行科名正确而属名是错的 (真为非植物病毒), 会引入 Plant 假阳性。
#   佐证: goji 44 行科级否决里 15 行、onekp 579 行里 463 行的属命中植物白名单。
#   明细: /tmp/w3_rescued_{goji,onekp}.tsv ; 审计 scripts/audit/HOST_REGISTRY_AUDIT_20260915.md
WHITELIST_OVERRIDES_FAMILY_VETO = False
"""

B = """    if FAMILY_FIRST_VETO and fam and fam in NON_PLANT_FAMILIES_FALLBACK:
        return True
"""
B_NEW = """    if FAMILY_FIRST_VETO and fam and fam in NON_PLANT_FAMILIES_FALLBACK:
        if WHITELIST_OVERRIDES_FAMILY_VETO:
            gen_w = str(genus).strip() if genus is not None else ''
            if is_trusted_level(det_level) and gen_w in PLANT_GENERA_WHITELIST:
                return False
        return True
"""

text = TARGET.read_text(encoding="utf-8")
if "WHITELIST_OVERRIDES_FAMILY_VETO" in text:
    print("[SKIP] 开关已存在")
    sys.exit(0)
for name, anchor in (("A", A), ("B", B)):
    if text.count(anchor) != 1:
        print("[FAIL] 锚点 %s 命中 %d 次 (需恰好 1)" % (name, text.count(anchor)))
        sys.exit(1)
text = text.replace(A, A_ADD, 1).replace(B, B_NEW, 1)
if "--check" in sys.argv:
    print("[CHECK] 锚点就位, 可打补丁")
    sys.exit(0)
TARGET.write_text(text, encoding="utf-8")
print("[OK] 已加入 WHITELIST_OVERRIDES_FAMILY_VETO (默认 False)")
