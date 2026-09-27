#!/bin/bash
# Runs detail.py once over all 16 tables (new + old), parsing dmp a single time.
set -u
ARGS=""
for d in Alternaria amarum Aphis barbarum chinense Fusarium onekp ruthenicum; do
  p=$(cut -f2 /tmp/refresh_v67/$d/manifest.tsv)
  ARGS="$ARGS $p /tmp/refresh_v67/$d/old.tsv"
done
python3 /tmp/indep_verify/detail.py /home/zhangwenda/database/taxonomy/rankedlineage.dmp $ARGS
echo DETAIL_ALL_DONE
