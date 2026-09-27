#!/bin/bash
# 只读侦察：全部物种项目清单 + 09 目录布局 + 已有的定点校准补丁脚本
echo "=== [1] ~/MMPV-paper 顶层（全部物种项目） ==="
ls -la ~/MMPV-paper/

echo
echo "=== [2] 各项目顶层目录布局 ==="
for d in ~/MMPV-paper/*/; do
  echo "--- $d"
  ls -1 "$d" 2>/dev/null | head -25
done

echo
echo "=== [3] 09 目录定位（各项目） ==="
for d in ~/MMPV-paper/*/; do
  find "$d" -maxdepth 3 -type d -name "09*" 2>/dev/null
done

echo
echo "=== [4] 已有定点校准/补丁脚本 ==="
find ~/MMPV-RNA/virome_discovery_pipeline -maxdepth 2 -name "*calib*" -o -maxdepth 2 -name "*patch*" 2>/dev/null | head -20
echo "--- /tmp 里的补丁脚本 ---"
ls -la /tmp/*calib* /tmp/*patch* /tmp/apply_* 2>/dev/null
echo "--- 项目内 scripts 目录 ---"
for d in ~/MMPV-paper/*/; do
  find "$d" -maxdepth 3 -type d -name "scripts" 2>/dev/null
done
