#!/bin/sh
set -e
P=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py
TS=20260830
cp -p "$P" "$P.bak_wlgenus_$TS"
python3 - <<'EOF'
import io, re
P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
src = io.open(P, encoding="utf-8").read()
old = "WHITELIST_OVERRIDES_FAMILY_VETO = False"
new = "WHITELIST_OVERRIDES_FAMILY_VETO = True"
n = src.count(old)
assert n == 1, "预期 1 处赋值, 实际 %d" % n
io.open(P, "w", encoding="utf-8", newline="").write(src.replace(old, new))
print("已替换 1 处: %s -> %s" % (old, new))
EOF
python3 -m py_compile "$P" && echo "py_compile OK"
md5sum "$P" "$P.bak_wlgenus_$TS"
grep -n "WHITELIST_OVERRIDES_FAMILY_VETO" "$P" | head -5
