#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""复现 batch_virema_dvg.py 的输入构建逻辑，检查 ERR2040134 走哪个分支"""
import re
from pathlib import Path

rd = Path('/home/zhangwenda/data-test/out/00b_HostDepletion')
sname = 'ERR2040134'
s_clean = sname.lower()
all_reads = [f for f in rd.rglob("*") if f.is_file()
             and any(ext in f.name.lower() for ext in ['.fq', '.fastq', '.fa', '.fasta', '.gz'])]
matched = [f for f in all_reads
           if re.search(r'\b' + re.escape(s_clean.replace('_', ' ')) + r'\b', f.name.lower().replace('_', ' ').replace('-', ' ').replace('.', ' '))
           or f.name.lower().startswith(s_clean + "_") or f.name.lower().startswith(s_clean + ".")]
matched = list(set(matched))
print('matched files:', [f.name for f in matched])

r1, r2 = None, None
for f in matched:
    nl = f.name.lower()
    if any(x in nl for x in ['_r2', '_2.', '.r2', '_2_']): r2 = f
    elif any(x in nl for x in ['_r1', '_1.', '.r1', '_1_']): r1 = f
    elif not r1: r1 = f
print('r1:', r1.name if r1 else None, ' r2:', r2.name if r2 else None)
is_single = not (r1 and r2)
print('is_single:', is_single)

is_fasta = bool(r1 and any(ext in r1.name.lower() for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']))
print('is_fasta:', is_fasta)
print('need_suffix:', (not is_single) and r1 and r2 and r1.exists() and r2.exists())

# 构建输入（与前 300 条），验证头名
import gzip
def _open_src(p):
    if str(p).endswith('.gz'):
        return gzip.open(p, 'rt')
    return open(p, 'r')

n = 0
with open('/tmp/repro_input.fa', 'w') as out_f:
    for src, tag in [(r1, '/1'), (r2, '/2')]:
        with _open_src(src) as in_f:
            for line in in_f:
                line = line.rstrip('\n\r')
                if line.startswith('>'):
                    if not line.rstrip().endswith(tag):
                        line = line.rstrip() + tag
                    out_f.write(line + '\n')
                    n += 1
                    if n >= 300: break
        if n >= 300: break
print('--- 复现输入头 4 行 ---')
for i, l in enumerate(open('/tmp/repro_input.fa')):
    print(l.rstrip())
    if i >= 3: break
