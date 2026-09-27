#!/bin/bash
# stage_breakdown.sh — 解析 VSI run.log, 输出各阶段累计耗时
# 用法: stage_breakdown.sh <run.log>
LOG="$1"
awk '
function ep(t){ sub(/\..*/, "", t); gsub(/[-:]/, " ", t); return mktime(t) }
/^[A-Za-z][A-Za-z0-9]*:$/ { st=$0; sub(/:$/, "", st) }
/Start time:/ { t0[st]=ep($2" "$3) }
/End time:/   { d=ep($2" "$3); if (st in t0) { dur[st]+=d-t0[st]; cnt[st]++; delete t0[st] } }
END {
  total=0
  for (k in dur) if (k!="extendOneScaffold" && k!="growScaffoldWithAssembly" && k!="checkCoverage" && k!="extendOneScaffold2") total+=dur[k]
  printf "%-30s %10s %8s %6s\n", "stage", "minutes", "share", "calls"
  for (k in dur) printf "%-30s %10.1f %7.1f%% %6d\n", k, dur[k]/60, 100*dur[k]/total, cnt[k]
  printf "%-30s %10.1f %7s\n", "TOTAL(leaf stages)", total/60, ""
}' "$LOG"
