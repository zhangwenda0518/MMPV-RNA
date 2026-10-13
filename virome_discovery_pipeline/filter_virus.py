#!/usr/bin/env python3
"""
02b_Filter: UniProt virus contig filter (CDD/mmseqs 层已移除, 2026-10-11).
UniProt: Diamond BLASTX + virus TaxID + keyword rescue

Modes:
  -m filter  UniProt(filter): viral TaxID / 关键词抢救 / 比对不可信抢救 / no-hit 保留
  -m strict  UniProt(strict): 仅 viral TaxID
  -m raw     symlink only

Split mode (merge-then-split):
  --split-by-prefix   Split output by contig prefix (first _ segment) into sample dirs
  -o split_outdir/    Base output dir for per-sample results
"""

import argparse, os, sys, subprocess, time
from collections import defaultdict

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))

DEFAULTS = {
    'UNIPROT_DB':     os.environ.get('UNIPROT_DB', os.path.expanduser('~/database/uniport_db/uniref90/uniref90')),
    'VIRUS_TAXID':    os.environ.get('VIRUS_TAXID', os.path.expanduser('~/MMPV-RNA/database/uniprot/virus.taxid.txt')),
    'THREADS':        os.environ.get('THREADS', '120'),
}

def run(cmd, desc=''):
    # fail-fast: 工具非零退出 (OOM/DB 路径错等) 时下游"输出为空 → symlink 全放行"
    # 会把崩溃与"零命中"混为一谈, 过滤层静默失效 → 直接报错退出
    print(f'  [{desc}] {cmd}', flush=True)
    t0 = time.time()
    proc = subprocess.run(cmd, shell=True, executable='/bin/bash',
                          stderr=subprocess.PIPE, text=True)
    wall = time.time() - t0
    print(f'  [{desc}] done: wall={wall:.1f}s', flush=True)
    if proc.returncode != 0:
        tail = '\n'.join((proc.stderr or '').strip().split('\n')[-5:]) or '(stderr 为空)'
        print(f'  [{desc}] ERROR: exit {proc.returncode}; stderr tail:\n{tail}',
              file=sys.stderr, flush=True)
        sys.exit(1)

def symlink_force(src, dst):
    if os.path.islink(dst) or os.path.exists(dst):
        os.unlink(dst)
    os.symlink(os.path.abspath(src), dst)

# ===================== UniProt =====================

def run_uniprot(input_fa, outdir, uniprot_db, threads, virus_taxid_path=None, mode='filter',
                force=False):
    base = os.path.basename(input_fa).rsplit('.', 1)[0]
    out_fa = os.path.join(outdir, f'{base}_uniprot.fasta')

    # 断点续传: 已有过滤结果则跳过 (--force 时忽略既有产物强制重算)
    if not force and os.path.isfile(out_fa) and os.path.getsize(out_fa) > 0:
        print(f'  UniProt: [resume] 已有 {out_fa}, 跳过', flush=True)
        return out_fa

    dmnd_db = uniprot_db if uniprot_db.endswith('.dmnd') else uniprot_db + '.dmnd'
    if not uniprot_db or not os.path.exists(dmnd_db):
        print('  UniProt: DB not found, symlinking input', flush=True)
        symlink_force(input_fa, out_fa)
        return out_fa

    sub_dir = os.path.join(outdir, 'diamond_uniprot')
    os.makedirs(sub_dir, exist_ok=True)
    blast_out = os.path.join(sub_dir, 'blastx.tsv')

    run(f'diamond blastx -q {input_fa} --db {uniprot_db} --long-reads '
        f'-o {blast_out} -e 1e-3 --threads {min(int(threads), 64)} '
        f'--max-target-seqs 10 --block-size 30 --index-chunks 2 '
        f'--outfmt 6 qseqid sseqid pident length evalue bitscore stitle '
        f'> {os.path.join(sub_dir, "blastx.log")} 2>&1',
        'uniprot-diamond')

    if not os.path.isfile(blast_out) or os.path.getsize(blast_out) == 0:
        # 仅当 diamond 成功 (run() 非零已退出) 且输出确为空时才走"零命中"放行
        print('  UniProt: 0 hits (diamond 成功且输出为空), symlinking input', flush=True)
        symlink_force(input_fa, out_fa)
        return out_fa

    virus_taxids = set()
    if virus_taxid_path and os.path.isfile(virus_taxid_path):
        with open(virus_taxid_path) as f:
            for l in f:
                l = l.strip()
                if l and l.isdigit():
                    virus_taxids.add(l)

    # ── 后置过滤: Top-N 多数投票 (polars 加速版, 完全复刻 virus_identification16.py) ──
    import polars as pl
    TOP_N = 10
    blast_cols = ['qseqid','sseqid','pident','length','evalue','bitscore','stitle']
    df = pl.read_csv(blast_out, separator='\t', has_header=False, new_columns=blast_cols,
                     schema_overrides={'pident': pl.Float64, 'length': pl.Int64, 'bitscore': pl.Float64, 'evalue': pl.Float64})
    if 'evalue' in df.columns:
        df = df.filter(pl.col('evalue') <= 0.001)
    df = df.sort(['qseqid','bitscore','evalue'], descending=[False, True, False])
    df_topn = df.group_by('qseqid').head(TOP_N)
    df_topn = df_topn.with_columns(
        pl.col('stitle').cast(pl.Utf8).str.extract(r'(?:OX|TaxID)=(\d+)', 1).alias('taxid'))

    all_fasta_ids = set()
    with open(input_fa) as fin:
        for line in fin:
            if line.startswith('>'):
                all_fasta_ids.add(line[1:].split()[0])
    all_df_qids = set(df['qseqid'].unique().to_list())

    viral_keywords = ['phage', 'virus', 'virion', 'capsid', 'tail', 'head',
                      'portal', 'terminase', 'integrase', 'prophage']

    def has_viral_keyword(stitle):
        sl = str(stitle if stitle is not None else '').lower()
        return any(kw in sl for kw in viral_keywords)

    def _compute_passed(mode):
        passed = set()
        for (qid,), group in df_topn.group_by('qseqid'):
            group_valid = group.filter(~pl.col('taxid').is_null())
            any_viral = any(t in virus_taxids for t in group_valid['taxid'].to_list()) if len(group_valid) > 0 else False
            if mode == 'strict':
                if any_viral:
                    passed.add(qid)
            else:  # filter
                if any_viral:
                    passed.add(qid)                         # ① 已知病毒
                elif len(group_valid) == 0:
                    pass  # 无 taxid → 留给全局 no-hit
                elif any(has_viral_keyword(s) for s in group_valid['stitle'].to_list()):
                    passed.add(qid)                         # ② 关键词抢救
                else:
                    for row in group_valid.iter_rows(named=True):
                        l = row.get('length', 0) or 0
                        p = row.get('pident', 0) or 0
                        e = row.get('evalue', 1.0) or 1.0
                        if l < 50 or p < 30 or e > 1e-5:
                            passed.add(qid); break          # ③ 比对不可信抢救
        if mode != 'strict':
            passed |= (all_fasta_ids - all_df_qids)        # ④ 全局 no-hit 保留
        return passed

    passed = _compute_passed(mode)

    keep_count = 0
    total_count = 0
    with open(input_fa) as fin, open(out_fa, 'w') as fout:
        write = False
        for line in fin:
            if line.startswith('>'):
                total_count += 1
                contig = line[1:].split()[0]
                write = contig in passed
                if write: keep_count += 1
            if write:
                fout.write(line)

    # 输出 UniProt unclassified (blastx 无任何 hit) 的 contig ID 列表
    hit_contigs = set()
    with open(blast_out) as f:
        for line in f:
            p = line.split('\t')
            if p:
                hit_contigs.add(p[0])
    uncl_id = os.path.join(outdir, f'{base}_uniprot.unclassified.id')
    uncl_set = set()
    with open(input_fa) as fin:
        for line in fin:
            if line.startswith('>'):
                contig = line[1:].split()[0]
                if contig not in hit_contigs:
                    uncl_set.add(contig)
    with open(uncl_id, 'w') as uf:
        for c in sorted(uncl_set):
            uf.write(c + '\n')
    print(f'  UniProt unclassified ID → {uncl_id} ({len(uncl_set)} contigs)', flush=True)

    print(f'  UniProt: kept {keep_count}/{total_count} contigs', flush=True)
    return out_fa

# ===================== Split =====================

def _sample_of(contig, sample_names=None):
    """contig header → 样本名。sample_names 提供时按最长前缀匹配 (样本名含 _ 时
    不再被首段截断, 如样本 SRR1_clean 的 contig SRR1_clean_k141_1 正确归 SRR1_clean);
    未提供时维持旧行为 (首段 _ 前缀)。匹配要求样本名后跟 _ 或 句点, 防 SRR1 误配
    到 SRR10 的 contig; 无匹配时回退旧行为。"""
    if sample_names:
        best = None
        for s in sample_names:
            if (contig == s or contig.startswith(s + '_')
                    or contig.startswith(s + '.')):
                if best is None or len(s) > len(best):
                    best = s
        if best is not None:
            return best
    return contig.split('_', 1)[0]


def split_output(split_outdir, uniprot_fa, uniprot_dir, base, sample_names=None):
    """Split merged output by contig prefix → per-sample dirs with uniprot/.
    After split, remove merged uniprot/ dir."""
    from collections import defaultdict

    # 确定源文件: UniProt 结果
    if not uniprot_fa or not os.path.isfile(uniprot_fa):
        print('Split: 无有效过滤结果 (UniProt 未产出), 跳过拆分')
        return

    # Pass 1: collect sample → contig set
    sample_contigs = defaultdict(set)
    with open(uniprot_fa) as f:
        for line in f:
            if line.startswith('>'):
                contig = line[1:].split()[0]
                sample = _sample_of(contig, sample_names)
                sample_contigs[sample].add(contig)

    all_samples = sorted(sample_contigs.keys())
    print(f'\nSplit: {len(all_samples)} samples -> {split_outdir}')

    def _split_fasta(src_fa, subdir, outname):
        """把合并 fasta 按样本前缀拆到 split_outdir/<s>/subdir/outname.
        利用输入按样本分组的特性: 单文件句柄顺序写, 样本切换时才换文件."""
        cur_fh = None
        cur_sample = None
        try:
            with open(src_fa) as fin:
                for line in fin:
                    if line.startswith('>'):
                        contig = line[1:].split()[0]
                        sample = _sample_of(contig, sample_names)
                        # 样本切换 → 关旧开新 (输入按样本分组, 切换次数=样本数-1)
                        if sample != cur_sample:
                            if cur_fh:
                                cur_fh.close()
                            out_dir = os.path.join(split_outdir, sample, subdir)
                            os.makedirs(out_dir, exist_ok=True)
                            cur_fh = open(os.path.join(out_dir, outname), 'w')
                            cur_sample = sample
                    if cur_fh:
                        cur_fh.write(line)
        finally:
            if cur_fh:
                cur_fh.close()

    # UniProt 结果
    if uniprot_fa and os.path.isfile(uniprot_fa):
        _split_fasta(uniprot_fa, 'uniprot', f'{base}_uniprot.fasta')

    # 并行: fasta 拆分已完成 (串行一遍), blast/id 并行
    from concurrent.futures import ThreadPoolExecutor, as_completed

    # 拆分 blastx.tsv (UniProt diamond 结果)
    def _split_blast():
        blast_tsv = os.path.join(uniprot_dir, 'diamond_uniprot', 'blastx.tsv') if uniprot_dir else None
        if blast_tsv and os.path.isfile(blast_tsv):
            sample_blast = defaultdict(list)
            with open(blast_tsv) as f:
                for line in f:
                    contig = line.split('\t', 1)[0]
                    sample = _sample_of(contig, sample_names)
                    sample_blast[sample].append(line)
            for s in all_samples:
                if sample_blast[s]:
                    s_uni = os.path.join(split_outdir, s, 'uniprot', 'diamond_uniprot')
                    os.makedirs(s_uni, exist_ok=True)
                    with open(os.path.join(s_uni, 'blastx.tsv'), 'w') as f:
                        for l in sample_blast[s]:
                            f.write(l)

    # 拆分 unclassified ID 列表: 每样本一行到 <s>/uniprot/
    def _split_id(id_file, subdir, outname):
        if not id_file or not os.path.isfile(id_file):
            return
        sample_ids = defaultdict(list)
        with open(id_file) as f:
            for line in f:
                contig = line.strip()
                if not contig:
                    continue
                sample = _sample_of(contig, sample_names)
                sample_ids[sample].append(contig)
        for s, ids in sample_ids.items():
            s_dir = os.path.join(split_outdir, s, subdir)
            os.makedirs(s_dir, exist_ok=True)
            with open(os.path.join(s_dir, outname), 'w') as f:
                f.write('\n'.join(ids) + '\n')

    # UniProt unclassified.id 实际由 run_uniprot 生成, 用 glob 匹配避免文件名不一致
    uniprot_uncl_id = None
    if uniprot_dir and os.path.isdir(uniprot_dir):
        import glob as _glob
        _hits = _glob.glob(os.path.join(uniprot_dir, '*_uniprot.unclassified.id'))
        if _hits:
            uniprot_uncl_id = _hits[0]
        else:
            _hits = _glob.glob(os.path.join(uniprot_dir, '*.unclassified.id'))
            if _hits:
                uniprot_uncl_id = _hits[0]

    jobs = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        blast_tsv_check = os.path.join(uniprot_dir, 'diamond_uniprot', 'blastx.tsv') if uniprot_dir else None
        if blast_tsv_check and os.path.isfile(blast_tsv_check):
            jobs.append(ex.submit(_split_blast))
        if uniprot_uncl_id and os.path.isfile(uniprot_uncl_id):
            jobs.append(ex.submit(_split_id, uniprot_uncl_id, 'uniprot', 'uniprot.unclassified.id'))
        for j in as_completed(jobs):
            j.result()

    # COBRA 兼容 symlink: 分层输出 02a原始 → UniProt
    ident_dir = os.path.normpath(os.path.join(split_outdir, '..', _D['d_ident']))
    for s in all_samples:
        sdir = os.path.join(split_outdir, s)
        uni_out = os.path.join(sdir, 'uniprot', f'{base}_uniprot.fasta')
        # 02a 原始候选
        orig_fa = os.path.join(ident_dir, s, f'{s}_virus.all.candidate.fasta')
        if not os.path.isfile(orig_fa):
            orig_fa = os.path.join(ident_dir, s, f'{s}.virus.candidate.fasta')
        if os.path.isfile(orig_fa):
            symlink_force(orig_fa, os.path.join(sdir, f'{s}_virus.all.candidate.fasta'))
        # UniProt 结果 (仅 UniProt 真跑过, 无回退)
        if os.path.isfile(uni_out):
            symlink_force(uni_out, os.path.join(sdir, f'{s}.virus.candidate_uniprot_filtered.fasta'))
        # 综合: UniProt > 02a 原始
        final_fa = uni_out if os.path.isfile(uni_out) else orig_fa
        if os.path.isfile(final_fa):
            symlink_force(final_fa, os.path.join(sdir, f'{s}.virus.candidate_filtered.fasta'))

    # 拆分完成后删除合并的 uniprot/
    import shutil
    if uniprot_dir and os.path.isdir(uniprot_dir):
        shutil.rmtree(uniprot_dir, ignore_errors=True)
        print(f'  removed merged dir: {uniprot_dir}')

    print(f'Split complete: {len(all_samples)} samples')

# ===================== FASTA 工具 =====================

def _fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                ids.append(line[1:].split()[0] if len(line) > 1 else '')
    return ids


# ===================== COBRA 主链镜像 (恢复 filter→COBRA 主链契约) =====================

def build_cobra_mirror(outdir, base, filter_mode, split_by_prefix=False,
                       verdict_dir=None, mirror_root=None):
    """把 filter 最终候选集以 02a 兼容结构镜像到
    {d_filter}/candidates/{S}/{S}_virus.all.candidate.fasta, 供 COBRA 主链消费。

    镜像根统一为 filter 输出根的 candidates (2026-10-13 用户裁定):
      - merged 形态 (--split-by-prefix, outdir = filter 根): outdir/candidates
      - per-sample 形态 (单样本/--no-merge, outdir = filter 根/<S>): 上提一层到
        os.path.dirname(outdir)/candidates (旧行为错误落在 <outdir>/candidates
        即 d_filter/<S>/candidates, 与 COBRA 消费位 d_filter/candidates 错位)
    mirror_root 显式传入时优先于自动判定。

    口径 (2026-10-07 用户裁定): graph 只对 filter 后无证据 (unclassified 且无
    blastx hit) 的结果做过滤 — 该动作由 graph 阶段 (graph_propagation.py,
    本阶段之后运行) 的裁决完成, 本镜像原样搬运, 不做额外增删。

    背景: 原始设计 COBRA 吃 filter 结果 (编排器按 cdd>uniprot>raw 逐级选择;
    CDD 层已于 2026-10-11 移除, 现为 uniprot>raw);
    feb352f (2026-09-30) 因 02b 结构与 cobra raw 模式不兼容导致 8/8 崩溃, 被迫
    绕开 filter。本镜像以 cobra 认识的结构承载 filter 结果, 恢复主链契约且不
    复发崩溃。"""
    del verdict_dir  # 历史参数: 早期镜像曾在此做 graph_nonviral 剔除, 现口径收窄后不再需要
    entries = []
    persample = False   # 最终链文件直接位于 outdir → outdir 即样本目录 (per-sample 形态)
    if split_by_prefix:
        for d in sorted(os.listdir(outdir)):
            sdir = os.path.join(outdir, d)
            if not os.path.isdir(sdir):
                continue
            chain = os.path.join(sdir, f'{d}.virus.candidate_filtered.fasta')
            if os.path.isfile(chain):
                entries.append((d, chain))
    else:
        chain = os.path.join(outdir, f'{base}_filtered.fasta')
        alt = os.path.join(outdir, f'{base}_virus.all.candidate.fasta')
        pick = chain if os.path.exists(chain) else (alt if os.path.exists(alt) else None)
        if pick:
            entries.append((os.path.basename(os.path.normpath(outdir)), pick))
            persample = True
    if not entries:
        return
    if mirror_root is None:
        if persample:
            parent = os.path.dirname(os.path.normpath(outdir))
            if parent:
                mirror_root = os.path.join(parent, 'candidates')
            else:
                print(f'CobraMirror: WARN 无法从 {outdir} 推出 filter 根 (无父目录), '
                      f'镜像根回退 {outdir}/candidates', file=sys.stderr, flush=True)
                mirror_root = os.path.join(outdir, 'candidates')
        else:
            mirror_root = os.path.join(outdir, 'candidates')
    os.makedirs(mirror_root, exist_ok=True)
    tot = kept_tot = 0
    for s, chain in entries:
        mdir = os.path.join(mirror_root, s)
        os.makedirs(mdir, exist_ok=True)
        dst = os.path.join(mdir, f'{s}_virus.all.candidate.fasta')
        with open(chain, 'rb') as fi, open(dst, 'wb') as fo:
            fo.write(fi.read())
        kept_tot += len(_fasta_ids(dst))
        tot += 1
    print(f'CobraMirror: {tot} 样本镜像 -> {mirror_root} '
          f'(共 {kept_tot} 条最终候选; COBRA 主链恢复消费 filter 结果)', flush=True)


# ===================== Main =====================

def main():
    parser = argparse.ArgumentParser(description='UniProt virus contig filter')
    parser.add_argument('-i', '--input', help='Input fasta (use with -o for single sample)')
    parser.add_argument('--input-dir', help='Input directory of per-sample candidates (merges automatically)')
    parser.add_argument('-o', '--outdir', required=True)
    parser.add_argument('-m', '--filter-mode', choices=['filter','strict','raw'], default='filter')
    parser.add_argument('--skip-uniprot', action='store_true')
    parser.add_argument('--uniprot-db', default=DEFAULTS['UNIPROT_DB'])
    parser.add_argument('--virus-taxid', default=DEFAULTS['VIRUS_TAXID'])
    parser.add_argument('-t', '--threads', default=DEFAULTS['THREADS'])
    parser.add_argument('--split-by-prefix', action='store_true',
                        help='Split output by contig prefix (first _ segment) into sample dirs')
    parser.add_argument('--sample-names', default=None,
                        help='逗号分隔的已知样本名列表: 按最长前缀匹配确定 contig 归属样本 '
                             '(样本名含 _ 时不被首段截断); 缺省沿用首段 _ 截断旧行为')
    parser.add_argument('--force', action='store_true', default=False,
                        help='忽略既有产物强制重算 (跳过 resume 检查)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --input-dir: 合并目录下所有真病毒候选 (*_virus.all.candidate.fasta)
    if args.input_dir:
        import glob
        merged_fa = os.path.join(args.outdir, 'all_candidates.fasta')
        with open(merged_fa, 'w') as fout:
            for f in sorted(glob.glob(os.path.join(args.input_dir, '*', '*_virus.all.candidate.fasta'))):
                with open(f) as fin:
                    fout.write(fin.read())
        args.input = merged_fa
        if args.filter_mode != 'raw':
            args.split_by_prefix = True  # auto-split back

    if not args.input:
        print('ERROR: -i or --input-dir required')
        sys.exit(1)

    # Mode resolution
    if args.filter_mode == 'raw':
        args.skip_uniprot = True
    args.uniprot_mode = 'strict' if args.filter_mode == 'strict' else 'filter'

    print(f'Input:  {args.input}')
    print(f'Output: {args.outdir}')
    print(f'Mode:   {args.filter_mode}', flush=True)

    current_input = args.input
    base = os.path.basename(current_input).rsplit('.', 1)[0]

    # ---- UniProt ----
    uniprot_dir = os.path.join(args.outdir, 'uniprot')
    os.makedirs(uniprot_dir, exist_ok=True)

    if args.skip_uniprot:
        uniprot_out = None
        uniprot_out_link = os.path.join(uniprot_dir, f'{base}_uniprot.fasta')
        symlink_force(current_input, uniprot_out_link)
        print('UniProt: SKIP -> symlink')
    else:
        uniprot_out = run_uniprot(current_input, uniprot_dir, args.uniprot_db,
                                  args.threads, args.virus_taxid, args.uniprot_mode,
                                  force=args.force)

    # ---- Output ----
    # 确定实际运行的组件 (供 split 判断)
    uni_real = uniprot_out if not args.skip_uniprot else None      # UniProt 真跑过

    if args.split_by_prefix:
        sample_names = [s for s in (args.sample_names or '').split(',') if s]
        split_output(args.outdir, uni_real, uniprot_dir, base, sample_names)
    else:
        final_src = uni_real or current_input
        final_dst = os.path.join(args.outdir, f'{base}_filtered.fasta')
        symlink_force(final_src, final_dst)
        cobra_dst = os.path.join(args.outdir, f'{base}_virus.all.candidate.fasta')
        symlink_force(final_src, cobra_dst)
        print(f'Final: {final_dst} -> {final_src}')

    # ---- COBRA 主链镜像 (恢复 filter→COBRA 主链契约; 02a 兼容结构) ----
    build_cobra_mirror(args.outdir, base, args.filter_mode,
                       split_by_prefix=bool(args.split_by_prefix))

if __name__ == '__main__':
    main()
