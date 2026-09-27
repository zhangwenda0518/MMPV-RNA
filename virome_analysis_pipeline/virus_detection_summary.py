#!/usr/bin/env python3
"""
virus_detection_summary.py — 病毒检测结果综合表（panvirome 风格）
从 best + high_conf + consensus(不含N) + assembly/extraction 拼一张 per-virus 综合表。

列:
  Virus / Reference / best_n / high_n / consensus_n / consensus_noN_n /
  assembly_n / ref_length / assembly_avg_len

用法:
  python virus_detection_summary.py \
    --base ~/virus/data-2026/known_virus_all_v2 \
    --out virus_detection_summary.tsv
"""
import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))

ACC_RE = re.compile(r'([A-Z]{1,4}_?\d{4,9}\.\d+)')

def extract_acc(text):
    """从目录名/文件名提取 RefSeq 风格 accession（如 OR489165.1）。"""
    m = ACC_RE.search(str(text))
    return m.group(1) if m else None

def seq_stats(fasta_path):
    """统计一个 fasta：返回 (条目数, 平均长度, 不含N条数)。"""
    lens = []
    noN = 0
    cur = []
    with open(fasta_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur:
                    s = ''.join(cur)
                    lens.append(len(s))
                    if 'N' not in s.upper():
                        noN += 1
                    cur = []
            else:
                cur.append(line)
        if cur:
            s = ''.join(cur)
            lens.append(len(s))
            if 'N' not in s.upper():
                noN += 1
    return len(lens), (sum(lens) / len(lens) if lens else 0.0), noN

def seq_n_pct(fasta_path):
    """返回最长序列的 N 含量百分比（consensus/full 通常单条）。"""
    cur = []
    best = ''
    try:
        with open(fasta_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if cur:
                        s = ''.join(cur)
                        if len(s) > len(best):
                            best = s
                        cur = []
                else:
                    cur.append(line)
            if cur:
                s = ''.join(cur)
                if len(s) > len(best):
                    best = s
    except Exception:
        return None
    if not best:
        return None
    return round(best.upper().count('N') / len(best) * 100, 1)


def main():
    ap = argparse.ArgumentParser(description='病毒检测结果综合表')
    ap.add_argument('--base', required=True, help='known_virus_all_v2 根目录')
    ap.add_argument('--out', default='virus_detection_summary.tsv')
    ap.add_argument('--ref_info', default='', help='参考库 final.cluster.ref_info.tsv (节段判定, 默认自动搜索)')
    args = ap.parse_args()

    # 节段映射: 物种 -> 段集合 (ref_info)
    seg_map = {}
    ref_info_candidates = [args.ref_info] if args.ref_info else []
    ref_info_candidates += [
        str(Path.cwd() / 'database' / 'final.cluster.ref_info.tsv'),
        str(Path.home() / 'plant_virus_db' / '3.final-ref-virus.db' / 'final.cluster.ref_info.tsv'),
    ]
    for rc in ref_info_candidates:
        rp = Path(rc)
        if rp.exists():
            try:
                import csv as _csv
                with open(rp, newline='', encoding='utf-8', errors='replace') as f:
                    for r in _csv.DictReader(f, delimiter='\t'):
                        sp = (r.get('Species_NCBI') or r.get('Species_ICTV') or '').strip()
                        seg = (r.get('Segment') or '').strip()
                        if sp and seg:
                            seg_map.setdefault(sp, set()).add(re.sub(r'\s+', '', seg.upper()))
                print(f'[ref_info] loaded from {rp} ({len(seg_map)} species with segment)')
                break
            except Exception as e:
                print(f'[warn] ref_info parse failed {rp}: {e}')
    base = Path(args.base).expanduser()

    best_p = base / _D['a_detect'] / 'summary' / 'all_viruses.best.summary.tsv'
    high_p = base / _D['a_filter'] / 'high_conf.summary.tsv'
    if not best_p.exists() or not high_p.exists():
        print(f'[ERROR] missing {best_p} or {high_p}', file=sys.stderr)
        sys.exit(1)

    best = pd.read_csv(best_p, sep='\t', dtype=str)
    high = pd.read_csv(high_p, sep='\t', dtype=str)
    best['_sample'] = best['Sample'].astype(str)
    high['_sample'] = high['Sample'].astype(str)

    # accession → 物种名 / 参考长度
    acc2sp = best.drop_duplicates('Rep_Accession').set_index('Rep_Accession')['Adjusted_Species'].to_dict()
    acc2len = {}
    for acc, rl in best[['Rep_Accession', 'Rep_Length']].drop_duplicates('Rep_Accession').values:
        try:
            acc2len[acc] = int(float(rl))
        except (TypeError, ValueError):
            pass

    rows = {}

    def ensure(acc):
        if acc not in rows:
            rows[acc] = {
                'Virus': str(acc2sp.get(acc, acc)),
                'Reference': acc,
                'is_segmented': len(seg_map.get(str(acc2sp.get(acc, acc)), set())) > 1,
                'best_n': 0, 'high_n': 0,
                'consensus_n': 0, 'consensus_noN_n': 0,
                'consensus_n0': 0, 'consensus_n1_5': 0, 'consensus_n_ge5': 0,
                'assembly_n': 0, 'extract_n': 0,
                'assembly_avg_len': None, 'ref_length': acc2len.get(acc),
                'note': '',
            }
        return rows[acc]

    # best / high 检出数（每病毒 unique sample 数）
    for acc, grp in best.groupby('Rep_Accession'):
        ensure(str(acc))['best_n'] = grp['_sample'].nunique()
    for acc, grp in high.groupby('Rep_Accession'):
        ensure(str(acc))['high_n'] = grp['_sample'].nunique()

    # consensus: 03_variants/virus-consensus/<病毒>/<样本>/<样本>.consensus.fasta
    cons_root = base / _D['a_variants'] / 'virus-consensus'
    if cons_root.is_dir():
        for vdir in cons_root.iterdir():
            if not vdir.is_dir():
                continue
            acc = extract_acc(vdir.name)
            if not acc:
                continue
            for f in vdir.rglob('*.consensus.fasta'):
                n, _avg, noN = seq_stats(f)
                if n == 0:
                    continue
                row = ensure(acc)
                row['consensus_n'] += n
                row['consensus_noN_n'] += noN
                p = seq_n_pct(f)
                if p is not None:
                    if p == 0:
                        row['consensus_n0'] += 1
                    elif p < 5:
                        row['consensus_n1_5'] += 1
                    else:
                        row['consensus_n_ge5'] += 1

    # assembly: 05_assembly/<病毒>/<样本>/ 子目录数 = 组装成功样本数
    asm_root = base / _D['a_assembly']
    if asm_root.is_dir():
        for vdir in asm_root.iterdir():
            if not vdir.is_dir():
                continue
            acc = extract_acc(vdir.name)
            if not acc:
                continue
            n_dirs = sum(1 for x in vdir.iterdir() if x.is_dir())
            if n_dirs:
                ensure(acc)['assembly_n'] = n_dirs

    # extraction: 06_extraction/<病毒>/*.full.fasta
    ext_root = base / _D['a_extract']
    if ext_root.is_dir():
        for vdir in ext_root.iterdir():
            if not vdir.is_dir():
                continue
            acc = extract_acc(vdir.name)
            if not acc:
                continue
            row = ensure(acc)
            lens = []
            bad = []
            for f in vdir.glob('*.full.fasta'):
                n, avg, _noN = seq_stats(f)
                if n == 0:
                    bad.append(f.name)
                    continue
                if avg < 50:  # 异常/空产物
                    bad.append(f'{f.name}({int(avg)}bp)')
                    continue
                row['extract_n'] += n
                lens.append(avg)
            if bad:
                row['note'] = '; '.join(bad)[:120]
            if lens:
                row['assembly_avg_len'] = round(sum(lens) / len(lens))

    # extraction 显式状态 → note 自动填充原因
    status_p = ext_root / 'extraction_status.tsv' if ext_root.is_dir() else None
    if status_p and status_p.exists():
        try:
            st = pd.read_csv(status_p, sep='\t', dtype=str)
            for virus, grp in st.groupby('virus'):
                acc = extract_acc(virus)
                if not acc or acc not in rows:
                    continue
                cnt = grp['status'].value_counts().to_dict()
                skipped = {k: v for k, v in cnt.items() if k != 'extracted'}
                if skipped:
                    reason = '; '.join(f'{k}({v})' for k, v in sorted(skipped.items()))
                    rows[acc]['note'] = reason
        except Exception as e:
            print(f'[warn] read extraction_status failed: {e}')

    # consensus N 分布 → note 补充
    for acc, row in rows.items():
        if row['consensus_n']:
            dist = f"consN0:{row['consensus_n0']} N1-5:{row['consensus_n1_5']} N>=5:{row['consensus_n_ge5']}"
            row['note'] = (row['note'] + ' | ' + dist).strip(' |') if row['note'] else dist

    cols = ['Virus', 'Reference', 'is_segmented', 'best_n', 'high_n', 'consensus_n', 'consensus_noN_n',
            'consensus_n0', 'consensus_n1_5', 'consensus_n_ge5',
            'assembly_n', 'extract_n', 'assembly_avg_len', 'ref_length', 'note']
    out_df = pd.DataFrame([rows[k] for k in sorted(rows, key=lambda k: -rows[k]['best_n'])], columns=cols)
    out_df.to_csv(args.out, sep='\t', index=False)
    print(f'[OK] {args.out} ({len(out_df)} viruses)')


if __name__ == '__main__':
    main()