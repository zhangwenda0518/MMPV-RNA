#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读侦察：第 3 条诉求（按分类补 DNA/RNA）的可落地性。
1) VMR 的 Genome 列取值分布
2) 05 表 Realm x Kingdom 交叉表（看能否只靠分类层级判定 DNA/RNA）
3) 属/种级 join 的可行性：同一属内 Genome 是否唯一
4) 管线现有代码里 genome_type 的用法
"""
import collections
import csv
import os

VMR = '/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv'
G = '/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out'
O = '/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus'
P = '/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline'


def rd(path):
    with open(path, 'r', encoding='utf-8', errors='replace', newline='') as f:
        r = csv.reader(f, delimiter='\t', quotechar='"')
        hdr = next(r)
        hdr = [h.strip() for h in hdr]
        for row in r:
            yield hdr, row


def norm(x):
    if x is None:
        return None
    x = x.strip().strip('"').replace('*', '').strip()
    if x == '' or x.upper() in ('NA', 'N/A', '-'):
        return None
    return x


print('=== [1] VMR Genome 列取值 ===')
genome_cnt = collections.Counter()
genus_genomes = collections.defaultdict(set)
sp_genomes = collections.defaultdict(set)
for hdr, row in rd(VMR):
    if len(row) < 27:
        continue
    g = norm(row[25])
    if g:
        genome_cnt[g] += 1
    gn = norm(row[15])
    sp = norm(row[17])
    if gn and g:
        genus_genomes[gn].add(g)
    if sp and g:
        sp_genomes[sp].add(g)
for k, v in genome_cnt.most_common():
    print('  %-14s %d' % (k, v))

print()
print('=== [2] 属级 join 可行性 ===')
amb = [g for g, s in genus_genomes.items() if len(s) > 1]
print('  VMR 属总数 %d，Genome 取值不唯一的属 %d' % (len(genus_genomes), len(amb)))
for g in sorted(amb)[:10]:
    print('    %s: %s' % (g, '|'.join(sorted(genus_genomes[g]))))
amba = [s for s, v in sp_genomes.items() if len(v) > 1]
print('  VMR 种总数 %d，Genome 取值不唯一的种 %d' % (len(sp_genomes), len(amba)))

print()
print('=== [3] 05 表 Realm x Kingdom 交叉表（前 20） ===')
for label, path in (('GOJI', os.path.join(G, '05_Taxonomy/Votus.integrated/final_integrated_classification.tsv')),
                    ('ONEKP', os.path.join(O, '05_Taxonomy/Votus.integrated/final_integrated_classification.tsv'))):
    if not os.path.isfile(path):
        print('  %s 不存在' % path)
        continue
    ct = collections.Counter()
    kg_none = collections.Counter()
    n = 0
    for hdr, row in rd(path):
        if len(row) < 12:
            continue
        n += 1
        realm = norm(row[4]) or '(空)'
        king = norm(row[5]) or '(空)'
        ct[(realm, king)] += 1
        if king == '(空)':
            kg_none[realm] += 1
    print('  --- %s 行数 %d ---' % (label, n))
    for (realm, king), c in ct.most_common(20):
        print('    %-28s %-16s %d' % (realm[:28], king[:16], c))
    print('    [Kingdom 缺失的 Realm 分布] %s' % dict(kg_none.most_common(8)))

print()
print('=== [4] 管线里 genome_type 相关用法 ===')
for fn in ('integrated_summary.py', 'virome_pipeline.py'):
    p = os.path.join(P, fn)
    if not os.path.isfile(p):
        print('  %s 不存在' % fn)
        continue
    print('  --- %s ---' % fn)
    with open(p, 'r', encoding='utf-8', errors='replace') as f:
        for i, line in enumerate(f, 1):
            low = line.lower()
            if any(k in low for k in ('genome_type', 'nucleic', 'moltype', 'molecule_type', 'dna_or_rna')):
                print('    %5d %s' % (i, line.rstrip()[:150]))
