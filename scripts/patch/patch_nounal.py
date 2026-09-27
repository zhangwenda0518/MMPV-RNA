#!/usr/bin/env python3
"""补丁: bowtie2 加 --no-unal, 删两个 samtools view 中间层"""
p = '/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/Virseqimprover.py'
lines = open(p, encoding='utf-8').readlines()
Q = '"'
bowtie_hits = []
view_hits = []
for i, l in enumerate(lines):
    if Q + 'bowtie2 --local --threads ' + Q in l:
        bowtie_hits.append(i + 1)
        lines[i] = l.replace(Q + 'bowtie2 --local --threads ' + Q,
                             Q + 'bowtie2 --local --no-unal --threads ' + Q)
    if ' | samtools view -bS - | samtools view -h -F 0x04 -b - | ' in l:
        view_hits.append(i + 1)
        lines[i] = l.replace(' | samtools view -bS - | samtools view -h -F 0x04 -b - | ', ' | ')
print('bowtie hits:', bowtie_hits)
print('view hits:', view_hits)
assert len(bowtie_hits) == 2, f'expect 2 bowtie lines, got {len(bowtie_hits)}'
assert len(view_hits) == 2, f'expect 2 view lines, got {len(view_hits)}'
open(p, 'w', encoding='utf-8').writelines(lines)
print('PATCH_OK')
