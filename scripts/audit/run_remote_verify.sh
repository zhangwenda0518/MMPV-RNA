#!/bin/bash
# 服务器侧验证批处理 (避免 PowerShell 引号/进程替换被破坏)
set -u
cd /home/zhangwenda/MMPV-RNA/virome_discovery_pipeline || exit 1
echo "########## 1. 文本级 diff: 备份 vs 现装 ##########"
python3 /tmp/diff_check_plantwl.py
echo
echo "########## 2. 白名单 vs Plant.tsv 独立复算 ##########"
python3 /tmp/verify_whitelist_vs_planttsv.py
