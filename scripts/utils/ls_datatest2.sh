#!/bin/bash
cd /home/zhangwenda/data-test
echo "== 未显示的目录/文件 =="
ls -la | grep -v '^total' | awk '{print $NF}' | while read f; do
  case "$f" in
    blacklist_bak_20260902|gbk_files|heiguo_deg|host_depletion_logs|onekp_analysis|onekp-virus|out10|out5|raw|smoke|spades-res|validation_20260905|vs2_test|.|..) ;;
    *) if [ -d "$f" ]; then echo "DIR  $f : $(du -sh --max-depth=0 "$f" 2>/dev/null | cut -f1) : $(ls "$f" | wc -l) entries : $(ls "$f" | head -3 | tr '\n' '|')"; else echo "FILE $f : $(du -sh "$f" 2>/dev/null | cut -f1)"; fi ;;
  esac
done
