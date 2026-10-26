#!/usr/bin/env python3
# fix_sqn.py v2 — 提交前整理一条龙:
#   1. Core13 元数据填充 source.src (collection_date/geo/bioproject/biosample/sra/host)
#   2. Organism 用 suvtk taxonomy.tsv 的物种分类
#   3. 剔除参考序列 (NC_/PV_/HM_ 等 accession, 不可提交)
#   4. 序列 ID 加 lcl| 前缀 (table2asn 要求本地 ID 格式)
#   5. 重跑 suvtk table2asn → submission.sqn
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from datetime import datetime

from Bio import SeqIO

# bash 解析: Linux= /bin/bash; Windows(Git Bash) 经 shutil.which 找到 bash.exe。
# Git 常装在带空格路径, 故用 [BASH, '-c', cmd] 列表式调用 (shell=True+executable 会断)。
BASH = shutil.which('bash') or '/bin/bash'

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): ③ 输出目录名随 MMPV_IO_LAYOUT 解析
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))

# 以下三个路径均可用环境变量覆盖, 默认值为 2026-10-26 ssh 实测的服务器现网位置:
#   BASE   → MMPV_SUBMIT_BASE   (各数据集 *_out 目录的父目录)
#   MIX    → MMPV_SUBMIT_MIX    (mix 数据集输出目录)
#   CORE13 → MMPV_SUBMIT_CORE13 (Global_Unified_Metadata_Core13.tsv)
# 实测依据: 7 个 RNA-*_out 及其 submission_<ds>_virome 均在 02_novel_virus/ 下;
#           Core13 在 01_data/sra_rna.data4/global_metadata/ 下 (159,722 B, 07-29)。
# 旧默认 '/home/zhangwenda/virus/data-2026/...' 已整体不存在 (目录被迁移), 故更正。
# 注意: mix 数据集在服务器上**尚无提交产物** (11.merge_assembly/ 只有 merge 中间文件,
#       mix/out 不存在、全盘无 submission_mix_virome) → 用 mix 前须先跑提交管线,
#       或用 MMPV_SUBMIT_MIX 指向实际输出目录。
BASE = Path(os.environ.get('MMPV_SUBMIT_BASE',
                           '/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus'))
MIX = Path(os.environ.get('MMPV_SUBMIT_MIX',
                          '/home/zhangwenda/virus/ningxiagouqi/11.merge_assembly/mix/out'))
CORE13 = Path(os.environ.get('MMPV_SUBMIT_CORE13',
                             '/home/zhangwenda/MMPV-paper/goji-virome/01_data/sra_rna.data4'
                             '/global_metadata/Global_Unified_Metadata_Core13.tsv'))

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

# 本地测序序列 (无 SRA 元数据) 不留默认值: 空串 + 告警, 提交前必须人工补
# (旧版默认 '2024'/'China:Ningxia' 属于把推测写成观测)
DEFAULT_DATE = ''
DEFAULT_GEO = ''

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
    """日期格式化为 DD-Mmm-YYYY

    只接受完整日期 (%Y-%m-%d / %Y/%m/%d / %d-%b-%Y); 年份 (YYYY) 等不完整
    输入原样返回, 不补全为 '01-Jan-YYYY' (凭空精确到日 = 把推测写成观测)。
    """
    if not s:
        return ''
    for fm in ('%Y-%m-%d', '%Y/%m/%d', '%d-%b-%Y'):
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
    """改 template.sbt: 联系邮箱 (NCBI_EMAIL 环境变量) + 按数据集宿主的标题

    不再硬编码任何邮箱; 环境变量取不到时把邮箱替换为空串并告警。
    """
    host = HOST_MAP.get(name, 'plant')
    title = f'Plant virome of {host} in Ningxia, China'
    email = os.environ.get('NCBI_EMAIL', '').strip()
    if not email:
        print(f'[{name}] ⚠️ 未设置 NCBI_EMAIL 环境变量: .sbt 联系邮箱留空, 提交前必须补全')
    txt = sbt_path.read_text(encoding='utf-8')
    txt = txt.replace('zhangwenda@example.com', email)
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
    n_missing_date = 0
    n_missing_geo = 0
    for i, r in enumerate(kept, 1):
        sra = extract_sra(r.id)
        md = meta_all.get(sra, {}) if sra else {}
        if md:
            n_sra_meta += 1
        org = tax_map.get(r.id, '') or 'unclassified plant virus'
        date = fmt_date(md.get('CollectionDate', '')) or DEFAULT_DATE
        if not date:
            n_missing_date += 1
        geo = fmt_geo(md.get('Location', '')) or DEFAULT_GEO
        if not geo:
            n_missing_geo += 1
        bp = md.get('BioProject', '')
        bs = md.get('BioSample', '')
        host = md.get('ScientificName', '') or HOST_MAP.get(name, 'plant')
        # isolate = (物种, 样本) 组哈希: 同物种同样本的多条 contig (分段病毒各片段)
        # 共享同一 isolate → NCBI/suvtk 按 Isolate 归组; 哈希与 FASTA 行序无关
        # (修复旧 NXGQ### 顺序漂移), 跨数据集合并也不重名
        iso = 'NXGQ' + hashlib.md5(f'{org}|{sra or r.id}'.encode('utf-8')).hexdigest()[:5].upper()
        src_lines.append('\t'.join([
            r.id, org, iso, date, geo, '', bp, bs, sra, '',
            'TRUE', 'plant virome', host]))
    if n_missing_date or n_missing_geo:
        print(f'[{name}] ⚠️ 元数据缺失: Collection_date 留空 {n_missing_date}/{len(kept)}, '
              f'geo_loc_name 留空 {n_missing_geo}/{len(kept)} — 提交前必须人工补全')
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
    try:
        r = subprocess.run(
            [BASH, '-c',
             f'cd {sv} && suvtk table2asn -i {sub}/id_norm/sequences_norm.fna '
             f'-o {sv}/submission -s {sub}/id_norm/source_norm.src '
             f'-f {sub}/id_norm/featuretable_norm.tbl -t {sbt} -c {sub}/id_norm/comments_norm.cmt'],
            capture_output=True, text=True)
    except OSError as e:
        print(f'[{name}] ✗ 无法执行 table2asn (shell/bash 未找到: {e})')
        return
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
        # 数据集可能尚未跑到提交阶段 (mix 当前就没有产物): 明确跳过, 不要报"异常"
        sub = out_root / f'submission_{name}_virome'
        if not out_root.is_dir():
            hint = ('  ← mix 需先跑提交管线, 或用 MMPV_SUBMIT_MIX 指向实际输出目录'
                    if name == 'mix' else '')
            print(f'[{name}] ⊘ 跳过: 数据集目录不存在 ({out_root}){hint}')
            continue
        if not sub.is_dir():
            print(f'[{name}] ⊘ 跳过: 无提交准备目录 ({sub})')
            continue
        try:
            process(name, out_root, meta_all)
        except Exception as e:
            print(f'[{name}] 异常: {e}')


if __name__ == '__main__':
    main()
