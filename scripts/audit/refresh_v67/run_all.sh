#!/bin/bash
# Independent cross-verification runner (read-only on products; all temp under /tmp/indep_verify)
set -u
SC=/tmp/indep_chimera.py
LOGD=/tmp/indep_verify/logs
mkdir -p "$LOGD"
{
  for d in Alternaria amarum Aphis barbarum chinense Fusarium onekp ruthenicum; do
    p=$(cut -f2 /tmp/refresh_v67/$d/manifest.tsv)
    printf '%s\t%s\t%s\n' "$d.new" "$p" "$LOGD/$d.new.log"
    printf '%s\t%s\t%s\n' "$d.old" "/tmp/refresh_v67/$d/old.tsv" "$LOGD/$d.old.log"
  done
} > /tmp/indep_verify/joblist.tsv
cat /tmp/indep_verify/joblist.tsv | xargs -P 3 -n 3 bash -c 'python3 /tmp/indep_chimera.py "$1" > "$2" 2>&1; echo "[done] $0 rc=$?"'
echo "ALL_DONE"
