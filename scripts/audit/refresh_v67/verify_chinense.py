#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核 chinense 门槛 v6.7 重跑：行数 / contig 集合 / 矛盾率 / 各阶元取值变化。

自写复算，不复用 /tmp/chimera_multi.py。判据与参考实现对齐：
对每对相邻阶元 (H 高, L 低)，用 rankedlineage.dmp 查 L 的定型父级，
低阶元有参照父级且行内高阶元有值时才比，不等即矛盾。
只读，不改动任何文件。
"""
import csv
from collections import Counter, OrderedDict

RANKED = '/home/zhangwenda/database/taxonomy/rankedlineage.dmp'
OLD = '/tmp/refresh_v67/chinense/old.tsv'
NEW = ('/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/'
       'RNA-Lycium_chinense_out/05_Taxonomy/Votus.integrated/'
       'final_integrated_classification.tsv')
RANKS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
REF_IDX = {'Species': 2, 'Genus': 3, 'Family': 4, 'Order': 5,
           'Class': 6, 'Phylum': 7, 'Kingdom': 8, 'Realm': 9}
PAIRS = [('Realm', 'Kingdom', 'Realm'), ('Kingdom', 'Phylum', 'Kingdom'),
         ('Phylum', 'Class', 'Phylum'), ('Class', 'Order', 'Class'),
         ('Order', 'Family', 'Order'), ('Family', 'Genus', 'Family'),
         ('Genus', 'Species', 'Genus')]


def norm(v):
    if v is None:
        return None
    v = v.strip()
    if len(v) >= 2 and v[0] == '"' and v[-1] == '"':
        v = v[1:-1].strip()
    while v.endswith('*'):
        v = v[:-1].strip()
    if not v or v.upper() in ('NA', 'N/A'):
        return None
    return v


def read(path):
    rows = OrderedDict()
    with open(path, newline='', encoding='utf-8', errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            cid = norm(r.get('contig_id'))
            rows[cid] = {k.strip(): norm(v) for k, v in r.items()}
    return rows


old = read(OLD)
new = read(NEW)
print('OLD rows=%d  NEW rows=%d' % (len(old), len(new)))
so, sn = set(old), set(new)
print('contig set equal: %s  only_old=%d  only_new=%d'
      % (so == sn, len(so - sn), len(sn - so)))

need = set()
for src in (old, new):
    for r in src.values():
        for L in RANKS:
            v = r.get(L)
            if v:
                need.add(v.lower())
ref = {}
with open(RANKED, encoding='utf-8', errors='replace') as fh:
    for line in fh:
        parts = [p.strip().strip('|').strip() for p in line.split('\t|\t')]
        if len(parts) < 10:
            continue
        nm = parts[1].lower()
        if nm in need and nm not in ref:
            ref[nm] = {L: parts[REF_IDX[L]] for L in RANKS}
print('ref name hit=%d / need=%d' % (len(ref), len(need)))


def rate(rows):
    bcount = Counter()
    n_conf = 0
    for r in rows.values():
        bad = 0
        for h, l, rc in PAIRS:
            hv, lv = r.get(h), r.get(l)
            if not hv or not lv:
                continue
            e = ref.get(lv.lower())
            if e is None:
                continue
            pv = (e.get(rc) or '').strip()
            if pv and pv != hv:
                bcount['%s-%s' % (h, l)] += 1
                bad = 1
        n_conf += bad
    return bcount, n_conf


for tag, rows in (('OLD_v66f', old), ('NEW_v67', new)):
    b, n = rate(rows)
    print('[%s] rows=%d conflict_rows=%d rate=%.2f%%'
          % (tag, len(rows), n, 100.0 * n / len(rows)))
    for h, l, _ in PAIRS:
        print('   pair %-14s %d' % ('%s-%s' % (h, l), b['%s-%s' % (h, l)]))

keys = [k for k in old if k in new]
chg, lost = Counter(), Counter()
for k in keys:
    o, nw = old[k], new[k]
    for L in RANKS:
        if (o.get(L) or '') != (nw.get(L) or ''):
            chg[L] += 1
    for L in ('Order', 'Family'):
        if o.get(L) and not nw.get(L):
            lost[L] += 1
        elif not o.get(L) and nw.get(L):
            lost[L + '_gain'] += 1
print('common contigs=%d' % len(keys))
print('changed rows per rank:')
for L in RANKS:
    print('   %-8s %d' % (L, chg[L]))
print('Order filled->empty=%d  Family filled->empty=%d'
      % (lost['Order'], lost['Family']))
print('Order empty->filled=%d  Family empty->filled=%d'
      % (lost['Order_gain'], lost['Family_gain']))
