#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐行等价回归 (重活): 备份版 vs 现装版 decision_tree_cascade + is_blacklisted 全量比对。
结果同时写到 /tmp/regress_plantwl.log, 便于回捞。"""
import subprocess
import sys

with open("/tmp/regress_plantwl.log", "w", encoding="utf-8") as fh:
    p = subprocess.run([sys.executable, "-u", "/tmp/regress_plantwl_equivalence.py"],
                       stdout=fh, stderr=subprocess.STDOUT)
print("exit=%s" % p.returncode)
print(open("/tmp/regress_plantwl.log", encoding="utf-8").read()[-4000:])
