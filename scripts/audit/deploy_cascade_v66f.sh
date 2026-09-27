#!/bin/bash
# v6.6f 共识引擎上线（一次性做完：校验 -> 备份 -> 安装 -> 复验）
# 只换 R 共识脚本，不动管线版本戳，不动任何已有产物。
# 用法: bash deploy_cascade_v66f.sh [补丁路径]     默认 /tmp/virus_classifier_analysis.R.cascade_v66f
set -eu
SRC=${1:-/tmp/virus_classifier_analysis.R.cascade_v66f}
PIPE_DIR=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline
TARGET=$PIPE_DIR/virus_classifier_analysis.R
TS=$(date +%Y%m%d_%H%M%S)
BAK=$TARGET.bak_cascade_v66f_$TS

echo "=== 0. 上线前校验 ==="
[ -f "$SRC" ] || { echo "补丁不存在: $SRC"; exit 1; }
Rscript -e "invisible(parse('$SRC')); cat('PARSE OK\n')" 2>&1 | tail -2
grep -n 'TAX_GATE_VERSION <-' "$SRC"
grep -n '^TIE_BREAK_ORDER' "$SRC"
echo "补丁 md5: $(md5sum "$SRC" | cut -d' ' -f1)"

echo
echo "=== 1. 当前线上脚本 ==="
md5sum "$TARGET"
grep -n 'TAX_GATE_VERSION <-' "$TARGET"

echo
echo "=== 2. 备份 ==="
cp -p "$TARGET" "$BAK"
echo "备份 -> $BAK"
md5sum "$BAK"

echo
echo "=== 3. 安装 ==="
cp "$SRC" "$TARGET"
tr -d '\r' < "$TARGET" > "$TARGET.lf" && mv "$TARGET.lf" "$TARGET"
chmod 664 "$TARGET"
md5sum "$TARGET"
grep -n 'TAX_GATE_VERSION <-' "$TARGET"

echo
echo "=== 4. 复验（parse + 版本戳） ==="
Rscript -e "invisible(parse('$TARGET')); cat('TARGET PARSE OK\n')" 2>&1 | tail -2

echo
echo "如需回滚: cp -p $BAK $TARGET"
echo "备份清单:"; ls -la "$TARGET"*bak_cascade_v66f_* 2>/dev/null || true
