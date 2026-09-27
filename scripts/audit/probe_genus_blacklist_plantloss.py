#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E1 核实: Miraophiovirus / Oleurovirus 是否被属黑名单误否决, 影响多少行。

只读脚本: 不写任何管线文件。用法: python3 /tmp/probe_genus_blacklist_plantloss.py
"""
import os
import re
import sys
from collections import Counter

sys.path.insert(0, '/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline')
import run_host_prediction as R  # noqa: E402

TARGETS = ['Miraophiovirus', 'Oleurovirus', 'Ophiovirus', 'Oleavirus', 'Geminivirus']

print('=' * 78)
print('[1] 名单状态直测')
print('=' * 78)
for fam, gen, det in [
    ('Aspiviridae', 'Miraophiovirus', 'Family(via Family)'),
    ('Aspiviridae', 'Miraophiovirus', 'Genus(via Genus)'),
    ('Aspiviridae', 'Miraophiovirus', 'Species(via Species)'),
    ('Geminiviridae', 'Oleurovirus', 'Family(via Family)'),
    ('Geminiviridae', 'Oleurovirus', 'Genus(via Genus)'),
]:
    print('  %-14s %-16s %-20s -> is_blacklisted=%s' % (
        fam, gen, det, R.is_blacklisted(fam, gen, det)))
print('  in NON_PLANT_GENERA :', {g: (g in R.NON_PLANT_GENERA) for g in TARGETS})
print('  in W3 属白名单      :', {g: (g in R.PLANT_GENERA_WHITELIST) for g in TARGETS})
print('  Aspiviridae in fam白名单:', 'Aspiviridae' in R.PLANT_FAMILIES_WHITELIST,
      '/ in 科黑名单:', 'Aspiviridae' in R.NON_PLANT_FAMILIES_FALLBACK)
print('  Geminiviridae in fam白名单:', 'Geminiviridae' in R.PLANT_FAMILIES_WHITELIST,
      '/ in 科黑名单:', 'Geminiviridae' in R.NON_PLANT_FAMILIES_FALLBACK)
print('  FAMILY_FIRST_VETO =', R.FAMILY_FIRST_VETO,
      '/ WHITELIST_OVERRIDES_FAMILY_VETO =', R.WHITELIST_OVERRIDES_FAMILY_VETO)

TREES = {
    'goji': '/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out',
    'onekp': '/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus',
}


def find_tables(root):
    """找 05 整合表与 C9 表。"""
    hits = {'05': [], 'c9': []}
    for dirpath, dirnames, filenames in os.walk(root):
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            if fn == 'final_integrated_classification.tsv':
                hits['05'].append(p)
            if fn == 'classification_result.tsv' and 'C9_ICTV_result' in dirpath:
                hits['c9'].append(p)
    return hits


for tree, root in TREES.items():
    print()
    print('=' * 78)
    print('[2] 树 %s' % tree)
    print('=' * 78)
    hits = find_tables(root)
    for k in ('05', 'c9'):
        print('  %s 表: %d 个' % (k, len(hits[k])))
        for p in hits[k][:4]:
            print('     ', p, os.path.getsize(p))

    for p in hits['05'][:2]:
        with open(p, encoding='utf-8', errors='replace') as fh:
            header = fh.readline().rstrip('\n').split('\t')
            idx = {c: i for i, c in enumerate(header)}
            fam_col = idx.get('Family')
            gen_col = idx.get('Genus')
            print('  [05] %s' % p)
            print('       Family列=%s Genus列=%s 共%d列' % (fam_col, gen_col, len(header)))
            cnt = Counter()
            rows = []
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < len(header):
                    parts += [''] * (len(header) - len(parts))
                fam = parts[fam_col] if fam_col is not None else ''
                gen = parts[gen_col] if gen_col is not None else ''
                if gen in TARGETS or fam in ('Aspiviridae', 'Geminiviridae'):
                    cnt[(fam, gen)] += 1
                    if gen in ('Miraophiovirus', 'Oleurovirus'):
                        rows.append(parts)
                elif not fam_col:
                    if any(t in line for t in ('Miraophiovirus', 'Oleurovirus')):
                        cnt[('REGEX_HIT', '?')] += 1
            print('       Family/Genus 组合分布:')
            for k2, v in cnt.most_common(30):
                print('         %-24s %-18s %d' % (k2[0], k2[1], v))
            if rows:
                print('       命中行样例 (前 3):')
                for r in rows[:3]:
                    print('         ', ' | '.join(
                        '%s=%s' % (header[i], r[i]) for i in range(min(len(header), 26))
                        if header[i] in ('contig_id', 'Virus_name', 'Family', 'Genus',
                                         'Determination_Level', 'Predicted_Host', 'Final_Host',
                                         'Host_ICTV', 'primary_tool', 'Family*_agree',
                                         'Genus*_agree', 'Confidence_Level', 'Decision_Method')))

    for p in hits['c9'][:2]:
        with open(p, encoding='utf-8', errors='replace') as fh:
            header = fh.readline().rstrip('\n').split('\t')
            idx = {c: i for i, c in enumerate(header)}
            ph = idx.get('Predicted_Host')
            dl = idx.get('Determination_Level')
            f_col = idx.get('Family')
            g_col = idx.get('Genus')
            print('  [C9] %s' % p)
            print('       Predicted_Host列=%s Determination_Level列=%s Family列=%s Genus列=%s'
                  % (ph, dl, f_col, g_col))
            cnt = Counter()
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < len(header):
                    parts += [''] * (len(header) - len(parts))
                gen = parts[g_col] if g_col is not None else ''
                fam = parts[f_col] if f_col is not None else ''
                if gen not in TARGETS and fam not in ('Aspiviridae', 'Geminiviridae'):
                    continue
                host = parts[ph] if ph is not None else ''
                det = (parts[dl] if dl is not None else '')
                det_lv = det.split('(')[0]
                cnt[(fam, gen, host, det_lv)] += 1
            print('       (Family, Genus, Predicted_Host, 判定层级) 分布:')
            for k2, v in cnt.most_common(30):
                print('         %-22s %-16s %-10s %-9s %d' % (k2[0], k2[1], k2[2], k2[3], v))
print()
print('DONE')
