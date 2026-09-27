#!/bin/sh
# 服务器端跑 v6.4 规则层自测
cd /tmp || exit 1
tr -d '\r' < logic_selftest.R > logic_selftest_lf.R
echo "== patch md5 =="
md5sum /tmp/virus_classifier_analysis.R.cascade_v64
echo "== selftest md5 =="
md5sum /tmp/logic_selftest_lf.R
echo "== run =="
Rscript /tmp/logic_selftest_lf.R /tmp/virus_classifier_analysis.R.cascade_v64 2>&1 | tee /tmp/logic_selftest.out
echo "== Rscript exit = ${PIPESTATUS:-see-above} =="
