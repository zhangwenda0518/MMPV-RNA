# -*- coding: utf-8 -*-
"""干跑验证：改后的收集器产出的条目里，非 barbarum 的样本级文件数应为 0。

不写盘，只调用 collect() 并对每个条目的源路径做统计。
"""
import os
import re
import sys

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


sys.path.insert(0, f"{_GOJI}/00_scripts")
import extract_barbarum_complete as E  # noqa: E402

runs, _ = E.load_runs(f"{_GOJI}/01_data")
E.RUNS.clear()
E.RUNS.update(runs)          # 模拟 main() 填充全局
print('barbarum RUNS =', len(E.RUNS))

RX = re.compile(r'^(?:sample_)?((?:SRR|CRR|ERR|DRR)\d+)(?=[._-]|$)')
V2 = f"{_GOJI}/03_known_virus/known_virus_all_v2"


def bad_in(src, kind):
    """统计单个条目（文件或整树）里非 barbarum 的样本级条目数。"""
    if kind == 'file':
        m = RX.match(os.path.basename(src))
        return 1 if (m and m.group(1) not in E.RUNS) else 0
    n = 0
    for dp, dns, fns in os.walk(src):
        for fn in fns:
            m = RX.match(fn)
            if m and m.group(1) not in E.RUNS:
                n += 1
        for dn in dns:
            m = RX.match(dn)
            if m and m.group(1) not in E.RUNS:
                n += 1
    return n


total_files = 0
all_bad = 0
for st in ['04_post_analysis', '07_capheine', '08_similarity',
           '09_dvg', '10_report', '03_variants']:
    src = os.path.join(V2, st)
    if not os.path.isdir(src):
        print('MISSING', st)
        continue
    ents = E.collect(st, src, E.RUNS)
    bad = 0
    tot = 0
    for e in ents:
        p = e['src']
        k = e['kind']
        if k in ('file', 'table_tsv', 'parquet'):
            tot += 1
            bad += bad_in(p, 'file')
        else:
            if os.path.isdir(p):
                for _dp, _dns, fns in os.walk(p):
                    tot += len(fns)
            else:
                tot += 1
            bad += bad_in(p, 'tree')
    total_files += tot
    all_bad += bad
    print('%s %-18s entries=%-6d files=%-7d mixed_in=%d'
          % ('OK  ' if bad == 0 else 'FAIL', st, len(ents), tot, bad))

print('total files =', total_files, '| total mixed_in =', all_bad)
print('RESULT:', 'PASS' if all_bad == 0 else 'FAIL')
