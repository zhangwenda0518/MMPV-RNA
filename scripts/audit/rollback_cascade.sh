#!/bin/bash
# v6.6f 回滚：把线上 R 共识脚本还原到指定备份。
# 用法: bash rollback_cascade.sh [备份文件]    不给参数则列出可用备份并还原最近一个
set -eu
PIPE_DIR=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline
TARGET=$PIPE_DIR/virus_classifier_analysis.R
BAK=${1:-}

if [ -z "$BAK" ]; then
  ls -1t "$TARGET".bak_cascade_v66f_* 2>/dev/null | head -5 || true
  BAK=$(ls -1t "$TARGET".bak_cascade_v66f_* 2>/dev/null | head -1)
  [ -n "$BAK" ] || { echo "找不到任何 bak_cascade_v66f_* 备份"; exit 1; }
  echo "自动选定: $BAK"
fi
[ -f "$BAK" ] || { echo "备份不存在: $BAK"; exit 1; }

echo "回滚前: $(md5sum "$TARGET" | cut -d' ' -f1)"
cp -p "$BAK" "$TARGET"
tr -d '\r' < "$TARGET" > "$TARGET.lf" && mv "$TARGET.lf" "$TARGET"
echo "回滚后: $(md5sum "$TARGET" | cut -d' ' -f1)"
grep -n 'TAX_GATE_VERSION <-' "$TARGET"
Rscript -e "invisible(parse('$TARGET')); cat('PARSE OK\n')" 2>&1 | tail -2
