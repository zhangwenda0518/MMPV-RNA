#!/bin/bash
# 只读：在 NCBI names.dmp 里查这 6 个名字到底是什么（taxid / 名字 / 分类级别）
NAMES=/home/zhangwenda/database/taxonomy/names.dmp
for n in Sylvanvirus Crucivirus Rimosavirus Klosneuvirus Hokovirus Indivirus; do
  echo "===== $n ====="
  grep -w "$n" "$NAMES" | head -6 | awk -F'|' '{gsub(/^ +| +$/,"",$1); gsub(/^ +| +$/,"",$2); gsub(/^ +| +$/,"",$4); print "  taxid="$1"  name="$2"  class="$4}'
  echo "  ---- 关联名字里含 genus/acronym 的 variant 记录 ----"
  grep -w "$n" "$NAMES" | awk -F'|' '{gsub(/^ +| +$/,"",$4)} $4=="genbank acronym"||$4=="equivalent name"||$4=="synonym"{c++} END{print "  等价名/缩写条数 =", c+0}'
done
echo
echo "===== 顺带看 Crucivirus 的 taxid 对应的分类链（若在 names.dmp 有 taxid） ====="
grep -w "Crucivirus" "$NAMES" | head -3 | awk -F'|' '{gsub(/^ +| +$/,"",$1); print $1}' | while read t; do
  echo "-- taxid $t 在 rankedlineage.dmp 的行:"
  grep -w "^$t" /home/zhangwenda/database/taxonomy/rankedlineage.dmp 2>/dev/null | head -2
done
