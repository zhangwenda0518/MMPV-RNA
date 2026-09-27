#!/usr/bin/env python3
"""从 patched 版派生 3 个参数组合变体 (顺序 × blank_orphan)，锚点必须命中，否则报错退出"""
from pathlib import Path

SRC = Path("/tmp/virus_classifier_analysis.patched.R")
CALLLINE = "  wide <- harmonize_family_genus(wide, stacked, tool_weights)"
CALLNB = "  wide <- harmonize_family_genus(wide, stacked, tool_weights, blank_orphan = FALSE)"
HEAD = "  # \u79d1-\u5c5e\u6821\u51c6\u653e\u5728\u6700\u540e"
BLOCK = (
    "  wide <- harmonize_genus_species(wide)\n"
    + HEAD + "\uff0c\u4fdd\u8bc1\u6700\u7ec8\u8868\u91cc Genus \u4e00\u5b9a\u80fd\u8ffd\u6eaf\u5230\u300c\u81ea\u62a5\u79d1 == \u5171\u8bc6\u79d1\u300d\u7684\u5de5\u5177\n"
    "  # blank_orphan=FALSE \u65f6\u53ea\u66ff\u6362\u4e0d\u7f6e\u7a7a\uff08\u4fdd\u7559\u5b64\u7acb\u5c5e\uff0c\u4e0e harmonize_genus_species \u7684\u4fdd\u5b88\u53d6\u5411\u4e00\u81f4\uff09\n"
    + CALLLINE + "\n"
)
BLOCK_FIRST = CALLLINE + "\n  wide <- harmonize_genus_species(wide)\n"
BLOCK_FIRST_NB = CALLNB + "\n  wide <- harmonize_genus_species(wide)\n"

src = SRC.read_text(encoding="utf-8")
assert src.count(CALLLINE) == 1, "调用行锚点未命中或重复"
assert src.count(BLOCK) == 1, "调用块锚点未命中"
assert "blank_orphan" in src, "源文件未见 blank_orphan"

variants = {
    "r_v_last_noblank.R": src.replace(CALLLINE, CALLNB),
    "r_v_first_blank.R": src.replace(BLOCK, BLOCK_FIRST),
    "r_v_first_noblank.R": src.replace(BLOCK, BLOCK_FIRST_NB),
}
for name, txt in variants.items():
    assert txt != src, name + " 替换未生效"
    p = Path("/tmp") / name
    p.write_text(txt, encoding="utf-8")
    lines = txt.splitlines()
    hits = [(i, l.strip()) for i, l in enumerate(lines)
            if "harmonize_family_genus(wide" in l or "harmonize_genus_species(wide)" in l]
    print("== %s (%d B) ==" % (name, len(txt.encode("utf-8"))))
    for i, l in hits:
        print("   %d: %s" % (i + 1, l))
print("\nOK: 3 个变体已生成")
