#!/usr/bin/env python3
"""sync_sqn_from_csv.py v2 — GUI 编辑后的 CSV 一键回写 submission.sqn

用法 (服务器): python3 sync_sqn_from_csv.py <dataset> [csv_path]
  csv_path 默认 /tmp/submission_csv_<dataset>.csv (本地 sync_client.py 上传)
流程:
  1. 校验 CSV sequence_name 集合 == source_filled.src 保留集 (不一致即退出)
  2. 用 CSV 非占位值覆盖 source_filled.src 与 source_norm.src 对应行
  3. 按 CSV authors/gb-title 重建 template.sbt (首跑备份为 template_orig.sbt)
  4. 调 table2asn 重出 submission.sqn 并报告错误数
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

EMAIL = os.environ.get('NCBI_EMAIL', '').strip()   # 不硬编码邮箱; 未设置时留空 + 告警
CONTACT = ('Zhang', 'Wenda', '', 'Ningxia University',
           'College of Life Sciences', 'Yinchuan', 'Ningxia', 'China', '', '750021')
# 解包顺序: last, first, middle, affil, div, city, sub, country, street, postal

# 以下路径均可用环境变量覆盖, 默认值为 2026-10-26 ssh 实测的服务器现网位置:
#   BASE  → MMPV_SUBMIT_BASE (各数据集 *_out 目录的父目录)
#   MIX   → MMPV_SUBMIT_MIX  (mix 数据集输出目录)
#   T2ASN → MMPV_T2ASN       (table2asn 可执行文件绝对路径)
# 旧默认 '/home/zhangwenda/virus/data-2026/...' 已整体不存在 (目录被迁移), 故更正;
# 两份脚本原先的 MIX 前缀还不一致 (本文件多一段 data-2026/), 现统一。
BASE = Path(os.environ.get('MMPV_SUBMIT_BASE',
                           '/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus'))
MIX = Path(os.environ.get('MMPV_SUBMIT_MIX',
                          '/home/zhangwenda/virus/ningxiagouqi/11.merge_assembly/mix/out'))
DATASETS = {
    'barbarum':   BASE / 'RNA-Lycium_barbarum_out',
    'ruthenicum': BASE / 'RNA-Lycium_ruthenicum_out',
    'chinense':   BASE / 'RNA-Lycium_chinense_out',
    'amarum':     BASE / 'RNA-Lycium_amarum_out',
    'Fusarium':   BASE / 'RNA-Fusarium_nematophilum_out',
    'Alternaria': BASE / 'RNA-Alternaria_alternata_out',
    'Aphis':      BASE / 'RNA-Aphis_gossypii_out',
    'mix':        MIX,
}
HOST_MAP = {
    'barbarum': 'Lycium barbarum', 'ruthenicum': 'Lycium ruthenicum',
    'chinense': 'Lycium chinense', 'amarum': 'Lycium amarum',
    'Fusarium': 'Fusarium nematophilum', 'Alternaria': 'Alternaria alternata',
    'Aphis': 'Aphis gossypii', 'mix': 'Lycium spp. (mixed samples)',
}
# table2asn 可执行文件; 可用环境变量 MMPV_T2ASN 覆盖 (默认值与服务器现网一致)
T2ASN = os.environ.get('MMPV_T2ASN', '/home/zhangwenda/.pixi/bin/table2asn')

SRC_COLS = ['Sequence_ID', 'Organism', 'Isolate', 'Collection_date', 'geo_loc_name',
            'Lat_Lon', 'Bioproject', 'Biosample', 'SRA', 'Segment', 'Metagenomic',
            'Metagenome_source', 'Host']
# CSV 列名 -> src 列索引
CSV_TO_SRC = {'organism': 1, 'collection_date': 3, 'src-geo_loc_name': 4,
              'src-Lat_Lon': 5, 'bioproject': 6, 'biosample': 7, 'sra': 8,
              'src-Host': 12}
PLACEHOLDER_RE = re.compile(
    r'XXXX|YYYY|PRJNAXXXX|Country:Region|SAMNXXXXXXXX|XX\.\d+|^Last,\s*First$', re.I)


def norm_id(i):
    """norm 副本里的本地 ID 规则 (与 fix_sqn.py 保持一致)"""
    s = re.sub(r'[^A-Za-z0-9_]', '_', i)
    if len(s) > 45:
        s = s[:36] + '_' + hashlib.md5(s.encode()).hexdigest()[:8]
    return s


def esc(s):
    return str(s).replace('"', "'")


def is_placeholder(v):
    return not v or bool(PLACEHOLDER_RE.search(str(v)))


def parse_authors(text):
    out = []
    for part in (text or '').split(';'):
        part = part.strip()
        if not part or PLACEHOLDER_RE.search(part):
            continue
        if ',' in part:
            last, _, rest = part.partition(',')
            toks = rest.strip().split()
        else:
            toks = part.split()
            last, toks = (toks[0], toks[1:]) if toks else ('', [])
        out.append((last, toks[0] if toks else '', ' '.join(toks[1:])))
    return out or [('Zhang', 'Wenda', '')]


# ────────────────── template.sbt ──────────────────

def _nm_fields(a, d):
    return (f'{d}last "{esc(a[0])}",\n'
            f'{d}first "{esc(a[1])}",\n'
            f'{d}middle "{esc(a[2])}",\n'
            f'{d}initials "",\n'
            f'{d}suffix "",\n'
            f'{d}title ""')


def _entry(a, brace_ind):
    """names std 列表项: {\n  name name {...}\n}; 括号自带配平"""
    ni = brace_ind + '  '
    return (f'{brace_ind}{{\n'
            f'{ni}name name {{\n'
            + _nm_fields(a, ni + '  ') + '\n'
            + f'{ni}}}\n'
            + f'{brace_ind}}}')


def _affil(ind, email):
    lines = [f'{ind}affil std {{',
             f'{ind}  affil "{CONTACT[3]}",',
             f'{ind}  div "{CONTACT[4]}",',
             f'{ind}  city "{CONTACT[5]}",',
             f'{ind}  sub "{CONTACT[6]}",',
             f'{ind}  country "{CONTACT[7]}",',
             f'{ind}  street "{CONTACT[8]}",']
    if email:
        lines.append(f'{ind}  email "{EMAIL}",')
    lines.append(f'{ind}  postal-code "{CONTACT[9]}"')
    lines.append(f'{ind}}}')
    return '\n'.join(lines)


def _names_list(authors, list_open, entry_brace_ind):
    ent = ',\n'.join(_entry(a, entry_brace_ind) for a in authors)
    close_ind = ' ' * (len(entry_brace_ind) - 6)
    return (f'{list_open}names std {{\n'
            + ent + '\n'
            + f'{close_ind}}},' if False else
            f'{list_open}names std {{\n' + ent)


def build_sbt(name, authors, title):
    """结构与 template_orig.sbt 完全一致 (该版本经 table2asn 验证通过)"""
    cnm = ('      name name {\n' + _nm_fields(authors[0], '        ') + '\n      },')
    cit_entries = ',\n'.join(_entry(a, '        ') for a in authors)      # 缩进8 大括号
    pub_entries = ',\n'.join(_entry(a, '          ') for a in authors)    # 缩进10 大括号
    if not EMAIL:
        print('⚠️ 未设置 NCBI_EMAIL 环境变量: .sbt 联系邮箱留空, 提交前必须补全')
    return (
        'Submit-block ::= {\n'
        '  contact {\n'
        '    contact {\n'
        + cnm + '\n'
        + _affil('      ', email=bool(EMAIL)) + '\n'
        '    }\n'
        '  },\n'
        '  cit {\n'
        '    authors {\n'
        '      names std {\n'
        + cit_entries + '\n'
        '      },\n'
        + _affil('      ', email=False) + '\n'
        '    }\n'
        '  },\n'
        '  subtype new\n'
        '}\n'
        'Seqdesc ::= pub {\n'
        '  pub {\n'
        '    gen {\n'
        '      cit "unpublished",\n'
        '      authors {\n'
        '        names std {\n'
        + pub_entries + '\n'
        '        }\n'
        '      },\n'
        f'      title "{esc(title)}"\n'
        '    }\n'
        '  }\n'
        '}\n'
        + (('Seqdesc ::= user {\n'
            '  type str "Submission",\n'
            '  data {\n'
            '    {\n'
            '      label str "AdditionalComment",\n'
            f'      data str "ALT EMAIL:{EMAIL}"\n'
            '    }\n'
            '  }\n'
            '}\n') if EMAIL else '')
        + 'Seqdesc ::= user {\n'
        '  type str "Submission",\n'
        '  data {\n'
        '    {\n'
        '      label str "AdditionalComment",\n'
        f'      data str "Submission Title:{esc(title)}"\n'
        '    }\n'
        '  }\n'
        '}\n')


# ────────────────── 主流程 ──────────────────

def load_csv(csv_path):
    import csv
    with open(csv_path, encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    return {r['sequence_name'].strip(): r for r in rows}


def load_src(path):
    lines = path.read_text(encoding='utf-8').rstrip('\n').split('\n')
    return [ln.split('\t') for ln in lines]


def save_src(path, recs):
    path.write_text('\n'.join('\t'.join(r) for r in recs) + '\n', encoding='utf-8')


def pick_title(rows, ds):
    for r in rows.values():
        t = (r.get('gb-title') or '').strip()
        if t and not is_placeholder(t):
            return t
    return f'Plant virome of {HOST_MAP.get(ds, ds)} in Ningxia, China'


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    ds = sys.argv[1]
    csv_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(f'/tmp/submission_csv_{ds}.csv')
    root = DATASETS.get(ds)
    if not root:
        print(f'[ERROR] 未知数据集: {ds}（可选: {", ".join(DATASETS)}）')
        sys.exit(1)
    if not root.is_dir():
        hint = ''
        if ds == 'mix':
            hint = ('\n        mix 数据集目前无提交产物 (11.merge_assembly/ 仅 merge 中间文件);'
                    ' 先跑提交管线, 或用 MMPV_SUBMIT_MIX 指向实际输出目录。')
        print(f'[ERROR] 数据集目录不存在: {root}{hint}')
        sys.exit(1)
    sub = root / f'submission_{ds}_virome'
    sv, idn = sub / 'suvtk_submission', sub / 'id_norm'
    for need in (idn / 'source_filled.src', sv / 'submission.sqn'):
        if not need.exists():
            print(f'[ERROR] 缺少必需文件: {need} — 该数据集尚未完成提交准备, 先跑 fix_sqn.py')
            sys.exit(1)
    if not csv_path.is_file():
        print(f'[ERROR] CSV 不存在: {csv_path}')
        sys.exit(1)

    # 1. 读入并校验
    rows = load_csv(csv_path)
    filled = load_src(idn / 'source_filled.src')
    kept = {r[0] for r in filled[1:]}
    extra = set(rows) - kept
    missing = kept - set(rows)
    if extra or missing:
        print(f'[ERROR] CSV 与保留集不一致: 多余={sorted(extra)[:5]} 缺失={sorted(missing)[:5]}')
        sys.exit(1)

    # 2. 覆盖两份 src 的元数据列
    # CSV 按原始 sequence_name 键控; source_norm.src 的首列已被 norm_id 规范化
    # (点/横线→下划线), 只按原键反查会 100% 失配静默跳过 — 同时按规范化键兜底
    n_upd = 0
    n_miss = []
    rows_by_norm = {norm_id(k): v for k, v in rows.items()}
    for fn in ('source_filled.src', 'source_norm.src'):
        p = idn / fn
        recs = load_src(p)
        for r in recs[1:]:
            key = r[0][5:] if r[0].startswith('lcl|') else r[0]
            row = rows.get(key) or rows_by_norm.get(norm_id(key))
            if not row:
                n_miss.append(f'{fn}:{r[0]}')
                continue
            for col, idx in CSV_TO_SRC.items():
                v = (row.get(col) or '').strip()
                if v and not is_placeholder(v):
                    if r[idx] != v:
                        r[idx] = v
                        n_upd += 1
        save_src(p, recs)

    if n_miss:
        print(f'[WARN] {len(n_miss)} 行按原始 ID 与规范化 ID 均未命中 CSV, '
              f'元数据未更新 (前 5): {", ".join(n_miss[:5])}')

    # 3. 重建 template.sbt
    authors, seen = [], set()
    for r in rows.values():
        for a in parse_authors(r.get('authors') or ''):
            if a not in seen:
                seen.add(a)
                authors.append(a)
    title = pick_title(rows, ds)
    sbt = sv / 'template.sbt'
    orig_bak = sv / 'template_orig.sbt'
    if not orig_bak.exists() and sbt.exists():
        shutil.copy2(sbt, orig_bak)
    txt = build_sbt(ds, authors, title)
    assert txt.count('{') == txt.count('}'), 'build_sbt 括号不平衡!'
    sbt.write_text(txt, encoding='utf-8')

    # 4. table2asn 重出
    sqn = sv / 'submission.sqn'
    cmd = [T2ASN, '-i', str(idn / 'sequences_norm.fna'), '-o', str(sqn),
           '-t', str(sbt), '-f', str(idn / 'featuretable_norm.tbl'),
           '-src-file', str(idn / 'source_norm.src'),
           '-w', str(idn / 'comments_norm.cmt'), '-V', 'vb', '-a', 's']
    proc = subprocess.run(cmd, capture_output=True, text=True)
    err_n = len(re.findall(r'^\s+Error:', proc.stdout + proc.stderr, re.M))
    if proc.returncode != 0 or not sqn.exists() or sqn.stat().st_size < 10000:
        print(f'[FAIL] {ds}: table2asn rc={proc.returncode}, errors={err_n}')
        tail = (proc.stdout + proc.stderr).strip().splitlines()
        print('\n'.join(tail[-15:]))
        sys.exit(1)

    print(f'[OK] {ds}: .sqn {sqn.stat().st_size // 1024} KB | '
          f'行更新 {n_upd} | 作者 {len(authors)} 人 | 错误 {err_n}')
    print(f'     title: {title}')


if __name__ == '__main__':
    main()
