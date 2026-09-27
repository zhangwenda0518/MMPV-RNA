#!/usr/bin/env python3
# fix_sqn.py v2 — 提交前整理一条龙:
#   1. Core13 元数据填充 source.src (collection_date/geo/bioproject/biosample/sra/host)
#   2. Organism 用 suvtk taxonomy.tsv 的物种分类
#   3. 剔除参考序列 (NC_/PV_/HM_ 等 accession, 不可提交)
#   4. 序列 ID 加 lcl| 前缀 (table2asn 要求本地 ID 格式)
#   5. 重跑 suvtk table2asn → submission.sqn
import os
import re
import subprocess
import sys
from pathlib import Path
from datetime import datetime

from Bio import SeqIO

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): ③ 输出目录名随 MMPV_IO_LAYOUT 解析
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))

BASE = Path(os.environ.get('MMPV_SUBMIT_BASE', '/home/zhangwenda/virus/data-2026/data-test'))
MIX = Path(os.environ.get('MMPV_SUBMIT_MIX', '/home/zhangwenda/virus/ningxiagouqi/11.merge_assembly/mix/out'))
CORE13 = Path(os.environ.get('MMPV_SUBMIT_CORE13',
                            '/home/zhangwenda/virus/data-2026/sra_rna.data4/global_metadata/Global_Unified_Metadata_Core13.tsv'))

DATASETS = {
    'barbarum': BASE / 'RNA-Lycium_barbarum_out',
    'ruthenicum': BASE / 'RNA-Lycium_ruthenicum_out',
    'chinense': BASE / 'RNA-Lycium_chinense_out',
    'amarum': BASE / 'RNA-Lycium_amarum_out',
    'Fusarium': BASE / 'RNA-Fusarium_nematophilum_out',
    'Alternaria': BASE / 'RNA-Alternaria_alternata_out',
    'Aphis': BASE / 'RNA-Aphis_gossypii_out',
    'mix': MIX,
}

# 数据集 → 宿主 (用于无 SRA 本地序列的占位 + template.sbt 标题)
HOST_MAP = {
    'barbarum': 'Lycium barbarum',
    'ruthenicum': 'Lycium ruthenicum',
    'chinense': 'Lycium chinense',
    'amarum': 'Lycium amarum',
    'Fusarium': 'Fusarium nematophilum',
    'Alternaria': 'Alternaria alternata',
    'Aphis': 'Aphis gossypii',
    'mix': 'Lycium spp. (mixed samples)',
}

# 本地测序序列 (无 SRA 元数据) 的占位默认值, 提交前人工补
DEFAULT_DATE = '2024'
DEFAULT_GEO = 'China:Ningxia'

# NCBI 参考序列 accession 前缀 (提交时剔除)
# 兼容两种格式: RefSeq (NC_024458.1) 与 GenBank 传统 accession (HM590055.1, 无下划线)
ACC_RE = re.compile(r'^[A-Z]{2}_?\d+(\.\d+)?$')
SRA_RE = re.compile(r'([SCER]{2,3}R\d+)')


def load_core13():
    meta = {}
    with open(CORE13, encoding='utf-8-sig') as f:
        hdr = f.readline().strip().split('\t')
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < len(hdr):
                cols += [''] * (len(hdr) - len(cols))
            d = dict(zip(hdr, cols))
            run = d.get('Run', '').strip()
            if run:
                meta[run] = d
    return meta


def load_taxonomy(p):
    """{contig: 末段分类名}"""
    m = {}
    if not p.exists():
        return m
    with open(p) as f:
        hdr = f.readline().strip().split('\t')
        ci = next((i for i, h in enumerate(hdr) if h.lower() in ('contig', 'seq_id', 'sequence_id')), 0)
        ti = next((i for i, h in enumerate(hdr) if h.lower() in ('taxonomy', 'tax')), 1)
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) <= max(ci, ti):
                continue
            tax = cols[ti].strip()
            leaf = tax.split(';')[-1].strip() if tax else ''
            m[cols[ci]] = leaf
    return m


def extract_sra(contig):
    mm = SRA_RE.match(contig)
    return mm.group(1) if mm else ''


def fmt_date(s):
    if not s:
        return ''
    for fm in ('%Y-%m-%d', '%Y/%m/%d', '%d-%b-%Y', '%Y'):
        try:
            return datetime.strptime(s[:10], fm).strftime('%d-%b-%Y')
        except ValueError:
            pass
    return s


def fmt_geo(s):
    if not s:
        return ''
    parts = [p.strip().replace('_AI', '').strip() for p in s.split(',')]
    parts = [p for p in parts if p]
    if len(parts) >= 2:
        return parts[0] + ':' + parts[1]
    return parts[0] if parts else ''


def norm_id(i):
    s = re.sub(r'[^A-Za-z0-9_]', '_', i)
    if len(s) > 45:
        import hashlib
        s = s[:36] + '_' + hashlib.md5(s.encode()).hexdigest()[:8]
    return 'lcl|' + s


def keep_blocks(tbl, keep_ids):
    """按 >Feature 块过滤 tbl, 保留 keep_ids 的块"""
    blocks = []
    cur = None
    for line in open(tbl, encoding='utf-8'):
        if line.startswith('>Feature '):
            if cur is not None:
                blocks.append(cur)
            cur = [line]
        elif cur is not None:
            cur.append(line)
    if cur is not None:
        blocks.append(cur)
    out = []
    for b in blocks:
        hid = b[0][len('>Feature '):].strip()
        if hid in keep_ids:
            out.append(b)
    return out


def rewrite_sbt(sbt_path, name):
    """改 template.sbt: 真实邮箱 + 按数据集宿主的标题"""
    host = HOST_MAP.get(name, 'plant')
    title = f'Plant virome of {host} in Ningxia, China'
    txt = sbt_path.read_text(encoding='utf-8')
    txt = txt.replace('zhangwenda@example.com', 'zhangwenda05@163.com')
    txt = txt.replace('Plant virome of Lycium chinense in Ningxia, China', title)
    txt = txt.replace('Submission Title:None', f'Submission Title:{title}')
    sbt_path.write_text(txt, encoding='utf-8')


def process(name, out_root, meta_all):
    sub = out_root / f'submission_{name}_virome'
    sv = sub / 'suvtk_submission'
    feat = out_root / _D['d_rescue'] / 'suvtk.features_output'
    hypo = out_root / _D['d_rescue'] / 'analyze_hypothetical'
    fna = feat / 'reoriented_nucleotide_sequences.fna'
    tbl = hypo / 'featuretable_updated.tbl'
    src, cmt, sbt = sv / 'source.src', sv / 'comments.cmt', sv / 'template.sbt'
    tax_tsv = out_root / _D['d_rescue'] / 'suvtk.taxonomy_output' / 'taxonomy.tsv'

    if not all(p.exists() for p in [fna, tbl, src, cmt, sbt]):
        print(f'[{name}] 缺输入, 跳过')
        return

    tax_map = load_taxonomy(tax_tsv)
    recs = list(SeqIO.parse(str(fna), 'fasta'))

    # 分保留/参考
    kept, refs = [], []
    for r in recs:
        if ACC_RE.match(r.id):
            refs.append(r)
        else:
            kept.append(r)
    keep_ids = {r.id for r in kept}

    # 重建 source.src
    src_lines = ['Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\t'
                 'Lat_Lon\tBioproject\tBiosample\tSRA\tSegment\tMetagenomic\tMetagenome_source\tHost']
    n_sra_meta = 0
    for i, r in enumerate(kept, 1):
        sra = extract_sra(r.id)
        md = meta_all.get(sra, {}) if sra else {}
        if md:
            n_sra_meta += 1
        org = tax_map.get(r.id, '') or 'unclassified plant virus'
        date = fmt_date(md.get('CollectionDate', '')) or DEFAULT_DATE
        geo = fmt_geo(md.get('Location', '')) or DEFAULT_GEO
        bp = md.get('BioProject', '')
        bs = md.get('BioSample', '')
        host = md.get('ScientificName', '') or HOST_MAP.get(name, 'plant')
        iso = f'NXGQ{i:03d}'
        src_lines.append('\t'.join([
            r.id, org, iso, date, geo, '', bp, bs, sra, '',
            'TRUE', 'plant virome', host]))
    src_out = sub / 'id_norm' / 'source_filled.src'
    src_out.parent.mkdir(parents=True, exist_ok=True)
    src_out.write_text('\n'.join(src_lines) + '\n', encoding='utf-8')

    # 过滤并重建 cmt (保留序列行)
    cmt_lines = []
    with open(cmt, encoding='utf-8') as f:
        for line in f:
            c0 = line.split('\t', 1)[0]
            if c0 in keep_ids or line.startswith('Sequence_ID'):
                cmt_lines.append(line)
    cmt_out = sub / 'id_norm' / 'comments_filled.cmt'
    cmt_out.write_text(''.join(cmt_lines), encoding='utf-8')

    # 过滤 tbl
    tbl_blocks = keep_blocks(tbl, keep_ids)
    tbl_out = sub / 'id_norm' / 'featuretable_filled.tbl'
    with open(tbl_out, 'w', encoding='utf-8') as f:
        for b in tbl_blocks:
            f.writelines(b)

    # 过滤 fna
    fna_out = sub / 'id_norm' / 'sequences_filled.fna'
    with open(fna_out, 'w', encoding='utf-8') as f:
        SeqIO.write(kept, f, 'fasta')

    # ID 规范化 (lcl|) — 四个文件
    def nf_fasta(p, out):
        with open(p) as f, open(out, 'w') as g:
            for line in f:
                if line.startswith('>'):
                    parts = line[1:].split(None, 1)
                    head = norm_id(parts[0])
                    if len(parts) > 1:
                        head += ' ' + parts[1]
                    g.write('>' + head + '\n')
                else:
                    g.write(line)

    def nf_tbl(p, out):
        with open(p) as f, open(out, 'w') as g:
            for line in f:
                if line.startswith('>Feature '):
                    g.write('>Feature ' + norm_id(line[len('>Feature '):].strip()) + '\n')
                else:
                    g.write(line)

    def nf_tsv(p, out):
        with open(p) as f, open(out, 'w') as g:
            for i, line in enumerate(f):
                if i == 0 or not line.strip():
                    g.write(line)
                    continue
                cols = line.rstrip('\n').split('\t')
                cols[0] = norm_id(cols[0])
                g.write('\t'.join(cols) + '\n')

    nf_fasta(fna_out, sub / 'id_norm' / 'sequences_norm.fna')
    nf_tbl(tbl_out, sub / 'id_norm' / 'featuretable_norm.tbl')
    nf_tsv(src_out, sub / 'id_norm' / 'source_norm.src')
    nf_tsv(cmt_out, sub / 'id_norm' / 'comments_norm.cmt')

    # 改 template.sbt (真实邮箱 + 标题)
    rewrite_sbt(sbt, name)

    # table2asn
    r = subprocess.run(
        f'cd {sv} && suvtk table2asn -i {sub}/id_norm/sequences_norm.fna '
        f'-o {sv}/submission -s {sub}/id_norm/source_norm.src '
        f'-f {sub}/id_norm/featuretable_norm.tbl -t {sbt} -c {sub}/id_norm/comments_norm.cmt',
        shell=True, executable='/bin/bash', capture_output=True, text=True)
    sqn = sv / 'submission.sqn'
    if sqn.exists() and sqn.stat().st_size > 0:
        print(f'[{name}] ✓ .sqn 生成 ({sqn.stat().st_size // 1024} KB) | 保留 {len(kept)} 条, '
              f'剔参考 {len(refs)} 条, 元数据命中 {n_sra_meta}/{len(kept)}')
    else:
        print(f'[{name}] ✗ 失败 rc={r.returncode}')
        for ln in ((r.stdout or '') + (r.stderr or '')).splitlines()[-20:]:
            if 'BiopythonWarning' in ln or 'warnings.warn' in ln or 'site-packages' in ln:
                continue
            print('   ', ln)
    # 参考序列清单
    if refs:
        with open(sub / 'ref_sequences_removed.txt', 'w') as f:
            for r in refs:
                f.write(r.id + '\n')


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    meta_all = load_core13()
    print(f'Core13 元数据: {len(meta_all)} 条 Run')
    for name, out_root in DATASETS.items():
        if only and name != only:
            continue
        try:
            process(name, out_root, meta_all)
        except Exception as e:
            print(f'[{name}] 异常: {e}')


if __name__ == '__main__':
    main()
