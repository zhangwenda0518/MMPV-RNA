#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""修复 batch_virema_dvg.py 的 /1 /2 后缀位置 bug：
后缀必须插在 read 名 token 内（第一个空格之前），否则 ViReMa 取键时被丢弃，R1/R2 同名导致 KeyError。"""
import shutil, time

FN = '/home/zhangwenda/MMPV-RNA/virome_analysis_pipeline/batch_virema_dvg.py'
bak = FN + '.bak_virema_tag_%s' % time.strftime('%Y%m%d')
shutil.copy2(FN, bak)
print('backup ->', bak)

src = open(FN).read()
old = """                                if line.startswith('>'):
                                    # Ensure unique read names for ViReMa dict
                                    if not line.rstrip().endswith(tag):
                                        line = line.rstrip() + tag
                                out_f.write(line + '\\n')"""
new = """                                if line.startswith('>'):
                                    # Ensure unique read names for ViReMa dict.
                                    # Fix: tag must be inserted INSIDE the name token (before any
                                    # whitespace), otherwise ViReMa's MakeReadDict (Name.split()[0][1:])
                                    # drops it and R1/R2 with identical names cause KeyError.
                                    parts = line.split(None, 1)
                                    name_tok = parts[0]
                                    rest = (' ' + parts[1]) if len(parts) > 1 else ''
                                    if not name_tok.endswith(tag):
                                        name_tok = name_tok + tag
                                    line = name_tok + rest
                                out_f.write(line + '\\n')"""
assert src.count(old) == 1, 'pattern count=%d' % src.count(old)
open(FN, 'w').write(src.replace(old, new))
print('patched OK')

import py_compile
py_compile.compile(FN, doraise=True)
print('py_compile OK')
