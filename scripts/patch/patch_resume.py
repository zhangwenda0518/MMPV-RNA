#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""修复 auto_known_virus.py Stage1 checkpoint 跳过时仍执行裸 --resume 的 bug"""
import shutil, time

FN = '/home/zhangwenda/MMPV-RNA/virome_analysis_pipeline/auto_known_virus.py'
bak = FN + '.bak_resume1_%s' % time.strftime('%Y%m%d')
shutil.copy2(FN, bak)
print('backup ->', bak)

src = open(FN).read()
old = '''        if not args.no_resume and not args.force:
            parts.append("--resume")
        if not run(" ".join(parts), log, "batch_virus_depth"):
            sys.exit(1)'''
new = '''        if parts:
            if not args.no_resume and not args.force:
                parts.append("--resume")
            if not run(" ".join(parts), log, "batch_virus_depth"):
                sys.exit(1)'''
assert src.count(old) == 1, 'pattern count=%d' % src.count(old)
open(FN, 'w').write(src.replace(old, new))
print('patched OK')

import py_compile
py_compile.compile(FN, doraise=True)
print('py_compile OK')
