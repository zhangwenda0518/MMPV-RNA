#!/usr/bin/env python3
"""
02b_Filter: UniProt + CDD tiered virus contig filter.
UniProt: Diamond BLASTX + virus TaxID + keyword rescue
CDD:     mmseqs translated search + 3-tier classification

Modes:
  -m filter  UniProt(filter) + CDD(filter) → viral+viral_plant+viral_bacteria
  -m strict  UniProt(strict) + CDD(strict)  → viral only
  -m raw     symlink only

Split mode (merge-then-split):
  --split-by-prefix   Split output by contig prefix (first _ segment) into sample dirs
  -o split_outdir/    Base output dir for per-sample results
"""

import argparse, os, sys, subprocess, csv, re, time
from collections import Counter, defaultdict

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))

DEFAULTS = {
    'CDD_DB':        os.path.expanduser('~/database/cdd/cdd-db/cdd_db'),
    'CDD_WHITELIST': os.path.expanduser('~/MMPV-RNA/database/cdd/cdd_virus_final_v4.txt'),
    'CDD_CLASSIFIED': os.path.expanduser('~/MMPV-RNA/database/cdd/cdd_classified_taxid.tsv'),
    'UNIPROT_DB':     os.environ.get('UNIPROT_DB', os.path.expanduser('~/database/uniport_db/uniref90/uniref90')),
    'VIRUS_TAXID':    os.environ.get('VIRUS_TAXID', os.path.expanduser('~/MMPV-RNA/database/uniprot/virus.taxid.txt')),
    'THREADS':        os.environ.get('THREADS', '120'),
}

def run(cmd, desc=''):
    print(f'  [{desc}] {cmd}', flush=True)
    t0 = time.time()
    ret = subprocess.run(cmd, shell=True, executable='/bin/bash').returncode
    wall = time.time() - t0
    print(f'  [{desc}] done: wall={wall:.1f}s', flush=True)
    if ret != 0:
        print(f'  [{desc}] WARNING: exit {ret}', flush=True)

def symlink_force(src, dst):
    if os.path.islink(dst) or os.path.exists(dst):
        os.unlink(dst)
    os.symlink(os.path.abspath(src), dst)

# ===================== UniProt =====================

def run_uniprot(input_fa, outdir, uniprot_db, threads, virus_taxid_path=None, mode='filter'):
    base = os.path.basename(input_fa).rsplit('.', 1)[0]
    out_fa = os.path.join(outdir, f'{base}_uniprot.fasta')

    # 断点续传: 已有过滤结果则跳过
    if os.path.isfile(out_fa) and os.path.getsize(out_fa) > 0:
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
        print('  UniProt: no hits, symlinking input', flush=True)
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

# ===================== CDD =====================

def run_cdd(input_fa, outdir, cdd_db, whitelist_path, classified_path, threads, mode, orig_base=None):
    base = orig_base or os.path.basename(input_fa).rsplit('.', 1)[0]
    out_fa = os.path.join(outdir, f'{base}_cdd_{mode}.fasta')

    whitelist = set(open(whitelist_path).read().strip().split('\n'))

    tier = {}
    cdd_cat = {}
    # 逆转录酶 name 正则：Root 类里被 category 漏掉的病毒结构域（RT_LTR/RNase_HI_RT_Ty3 等）
    RT_RE = re.compile(r"(?i)\brt_|_rt_|_rt\b|reverse transcriptase|rvt_|rtvl")
    with open(classified_path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            acc = r['accession']
            cat = r['category']
            title = r.get('title', '')
            cdd_cat[acc] = cat
            if acc not in whitelist: continue
            if cat == 'Virus': tier[acc] = 1
            elif cat == 'Virus_title': tier[acc] = 1
            elif cat == 'Root' and RT_RE.search(title): tier[acc] = 1  # Root 里的逆转录酶
            elif cat == 'Root': tier[acc] = 3
            else: tier[acc] = 2
    for acc in whitelist:
        if acc not in tier:
            tier[acc] = 2  # default empirical, Root handled above

    mmd = os.path.join(outdir, 'mmseqs')
    os.makedirs(mmd, exist_ok=True)
    qdb = os.path.join(mmd, 'qDB')
    rdb = os.path.join(mmd, 'rDB')
    tsv = os.path.join(mmd, 'hits.tsv')

    run(f'mmseqs createdb {input_fa} {qdb} -v 1 --threads {threads}', 'mmseqs-createdb')
    run(f'mmseqs search {qdb} {cdd_db} {rdb} {mmd}/tmp --search-type 2 -s 4 '
        f'-v 1 --threads {threads} --max-seqs 300', 'mmseqs-search')
    run(f'mmseqs convertalis {qdb} {cdd_db} {rdb} {tsv} '
        f'--format-output query,target,evalue --threads {threads} -v 1', 'mmseqs-convertalis')

    contig_tiers = defaultdict(lambda: Counter())
    contig_hosts = defaultdict(lambda: Counter())

    with open(tsv) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 3: continue
            q, t, ev = p[0], p[1], float(p[2])
            if ev > 1e-3: continue
            if t in whitelist:
                ti = tier.get(t, 3)
                contig_tiers[q][f'tier{ti}'] += 1
            else:
                contig_tiers[q]['non_viral'] += 1
                cat = cdd_cat.get(t, '')
                if cat == 'Plant': contig_hosts[q]['Plant'] += 1
                elif cat == 'Bacteria': contig_hosts[q]['Bacteria'] += 1
                elif cat == 'Fungi': contig_hosts[q]['Fungi'] += 1
                elif cat in ('Animal', 'Metazoa'): contig_hosts[q]['Animal'] += 1
                elif cat == 'Archaea': contig_hosts[q]['Archaea'] += 1

    def dominant_host(contig):
        ch = contig_hosts.get(contig, Counter())
        if not ch: return None
        total = sum(ch.values())
        top = ch.most_common(1)[0]
        return top[0] if top[1] > total * 0.5 else top[0]

    def callit(contig):
        ct = contig_tiers.get(contig, Counter())
        t1, t2, t3 = ct.get('tier1', 0), ct.get('tier2', 0), ct.get('tier3', 0)
        nv = ct.get('non_viral', 0)
        host = dominant_host(contig)
        # T3 (Root/universal) treated same as T2 → ambiguous
        if t1 >= 1:
            return f'viral_{host.lower()}' if host else 'viral'
        if t2 >= 1 or t3 >= 1:
            return f'ambiguous_{host.lower()}' if host else 'ambiguous'
        if host: return f'non_viral_{host.lower()}'
        if nv > 0: return 'non_viral'
        return 'unclassified'

    all_contigs = set()
    with open(input_fa) as f:
        for line in f:
            if line.startswith('>'):
                all_contigs.add(line[1:].split()[0])

    kept = set()
    calls = Counter()
    contig_call = {}
    for contig in sorted(all_contigs):
        call = callit(contig)
        contig_call[contig] = call
        calls[call] += 1
        # 保留规则: 命中白名单 (tier1/2/3 → viral*/ambiguous*) 或 CDD 无命中 (unclassified)
        if mode == 'strict':
            keep = (call == 'viral')  # strict: 仅纯病毒域 (无宿主污染)
        else:
            keep = call.startswith('viral') or call.startswith('ambiguous') or call == 'unclassified'
        if keep:
            kept.add(contig)

    keep_count = 0
    with open(input_fa) as fin, open(out_fa, 'w') as fout:
        write = False
        for line in fin:
            if line.startswith('>'):
                write = line[1:].split()[0] in kept
                if write: keep_count += 1
            if write: fout.write(line)

    summary = os.path.join(outdir, 'cdd_calls.tsv')
    with open(summary, 'w') as f:
        f.write('contig\tcall\n')
        for contig in sorted(all_contigs):
            f.write(f'{contig}\t{contig_call[contig]}\n')

    # 输出 CDD unclassified (无任何 CDD 命中) 的 contig ID 列表
    if mode == 'filter':
        uncl_id = os.path.join(outdir, f'{base}_cdd.unclassified.id')
        n_uncl = 0
        with open(uncl_id, 'w') as f:
            for contig in sorted(all_contigs):
                if contig_call[contig] == 'unclassified':
                    f.write(contig + '\n')
                    n_uncl += 1
        print(f'  CDD unclassified ID → {uncl_id} ({n_uncl} contigs)', flush=True)

    print(f'  CDD ({mode}):')
    for c, n in calls.most_common():
        if mode == 'strict':
            keep = (c == 'viral')
        else:
            keep = c.startswith('viral') or c.startswith('ambiguous') or c == 'unclassified'
        marker = ' KEPT' if keep else ''
        print(f'    {c}: {n}{marker}')
    print(f'  Kept: {keep_count}/{len(all_contigs)}')

    return out_fa

# ===================== Split =====================

def split_output(split_outdir, cdd_fa, uniprot_fa, cdd_dir, uniprot_dir, calls_tsv, base, mode='filter'):
    """Split merged output by contig prefix → per-sample dirs with cdd/ + uniprot/.
    After split, remove merged cdd/ and uniprot/ dirs."""
    from collections import defaultdict

    # 确定源文件: UniProt 结果优先, 否则 CDD 结果
    filtered_fa = uniprot_fa if uniprot_fa and os.path.isfile(uniprot_fa) else cdd_fa
    if not filtered_fa or not os.path.isfile(filtered_fa):
        print('Split: 无有效过滤结果 (CDD 和 UniProt 均未产出), 跳过拆分')
        return

    # Pass 1: collect sample → contig set
    sample_contigs = defaultdict(set)
    with open(filtered_fa) as f:
        for line in f:
            if line.startswith('>'):
                contig = line[1:].split()[0]
                sample = contig.split('_', 1)[0]
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
                        sample = contig.split('_', 1)[0]
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

    # CDD 结果 (仅当存在)
    if cdd_fa and os.path.isfile(cdd_fa):
        _split_fasta(cdd_fa, 'cdd', f'{base}_cdd_{mode}.fasta')
    # UniProt 结果
    if uniprot_fa and os.path.isfile(uniprot_fa):
        _split_fasta(uniprot_fa, 'uniprot', f'{base}_uniprot.fasta')

    # 拆分 cdd_calls.tsv
    def _split_calls():
        sample_calls = defaultdict(list)
        if calls_tsv and os.path.isfile(calls_tsv):
            with open(calls_tsv) as f:
                next(f)
                for line in f:
                    contig = line.split('\t', 1)[0]
                    sample = contig.split('_', 1)[0]
                    sample_calls[sample].append(line)
            for s in all_samples:
                if sample_calls[s]:
                    s_cdd = os.path.join(split_outdir, s, 'cdd')
                    os.makedirs(s_cdd, exist_ok=True)
                    with open(os.path.join(s_cdd, 'cdd_calls.tsv'), 'w') as f:
                        f.write('contig\tcall\n')
                        for l in sample_calls[s]:
                            f.write(l)

    # 拆分 blastx.tsv (UniProt diamond 结果)
    def _split_blast():
        blast_tsv = os.path.join(uniprot_dir, 'diamond_uniprot', 'blastx.tsv') if uniprot_dir else None
        if blast_tsv and os.path.isfile(blast_tsv):
            sample_blast = defaultdict(list)
            with open(blast_tsv) as f:
                for line in f:
                    contig = line.split('\t', 1)[0]
                    sample = contig.split('_', 1)[0]
                    sample_blast[sample].append(line)
            for s in all_samples:
                if sample_blast[s]:
                    s_uni = os.path.join(split_outdir, s, 'uniprot', 'diamond_uniprot')
                    os.makedirs(s_uni, exist_ok=True)
                    with open(os.path.join(s_uni, 'blastx.tsv'), 'w') as f:
                        for l in sample_blast[s]:
                            f.write(l)

    # 并行: fasta 拆分已完成 (CDD/UniProt 串行但每文件一遍), calls/blast/id 并行
    from concurrent.futures import ThreadPoolExecutor, as_completed

    # 拆分 unclassified ID 列表: 每样本一行到 <s>/cdd/ 或 <s>/uniprot/
    def _split_id(id_file, subdir, outname):
        if not id_file or not os.path.isfile(id_file):
            return
        sample_ids = defaultdict(list)
        with open(id_file) as f:
            for line in f:
                contig = line.strip()
                if not contig:
                    continue
                sample = contig.split('_', 1)[0]
                sample_ids[sample].append(contig)
        for s, ids in sample_ids.items():
            s_dir = os.path.join(split_outdir, s, subdir)
            os.makedirs(s_dir, exist_ok=True)
            with open(os.path.join(s_dir, outname), 'w') as f:
                f.write('\n'.join(ids) + '\n')

    cdd_uncl_id = os.path.join(cdd_dir, f'{base}_cdd.unclassified.id') if cdd_dir else None
    # UniProt unclassified.id 实际由 run_uniprot 生成, base 是 CDD 过滤文件名 (all_candidates_cdd_filter),
    # 用 glob 匹配避免文件名不一致
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
        if calls_tsv and os.path.isfile(calls_tsv):
            jobs.append(ex.submit(_split_calls))
        blast_tsv_check = os.path.join(uniprot_dir, 'diamond_uniprot', 'blastx.tsv') if uniprot_dir else None
        if blast_tsv_check and os.path.isfile(blast_tsv_check):
            jobs.append(ex.submit(_split_blast))
        if cdd_uncl_id and os.path.isfile(cdd_uncl_id):
            jobs.append(ex.submit(_split_id, cdd_uncl_id, 'cdd', 'cdd.unclassified.id'))
        if uniprot_uncl_id and os.path.isfile(uniprot_uncl_id):
            jobs.append(ex.submit(_split_id, uniprot_uncl_id, 'uniprot', 'uniprot.unclassified.id'))
        for j in as_completed(jobs):
            j.result()

    # COBRA 兼容 symlink: 分层输出 02a原始 → CDD → UniProt
    ident_dir = os.path.normpath(os.path.join(split_outdir, '..', _D['d_ident']))
    for s in all_samples:
        sdir = os.path.join(split_outdir, s)
        cdd_ok = cdd_fa and os.path.isfile(cdd_fa)
        uni_ok = uniprot_fa and os.path.isfile(uniprot_fa)
        cdd_out = os.path.join(sdir, 'cdd', f'{base}_cdd_{mode}.fasta') if cdd_ok else None
        uni_out = os.path.join(sdir, 'uniprot', f'{base}_uniprot.fasta') if uni_ok else None
        # 02a 原始候选
        orig_fa = os.path.join(ident_dir, s, f'{s}_virus.all.candidate.fasta')
        if not os.path.isfile(orig_fa):
            orig_fa = os.path.join(ident_dir, s, f'{s}.virus.candidate.fasta')
        if os.path.isfile(orig_fa):
            symlink_force(orig_fa, os.path.join(sdir, f'{s}_virus.all.candidate.fasta'))
        # CDD 结果 (仅 CDD 真跑过)
        if cdd_ok:
            symlink_force(cdd_out, os.path.join(sdir, f'{s}.virus.candidate_cdd_filtered.fasta'))
        # UniProt 结果 (仅 UniProt 真跑过, 无回退)
        if uni_ok:
            symlink_force(uni_out, os.path.join(sdir, f'{s}.virus.candidate_uniprot_filtered.fasta'))
        # 综合: UniProt > CDD > 02a 原始
        final_fa = uni_out or cdd_out or orig_fa
        if os.path.isfile(final_fa):
            symlink_force(final_fa, os.path.join(sdir, f'{s}.virus.candidate_filtered.fasta'))

    # 拆分完成后删除合并的 cdd/ 和 uniprot/
    import shutil
    for d in (cdd_dir, uniprot_dir):
        if d and os.path.isdir(d):
            shutil.rmtree(d, ignore_errors=True)
            print(f'  removed merged dir: {d}')

    print(f'Split complete: {len(all_samples)} samples')

# ===================== Main =====================

def main():
    parser = argparse.ArgumentParser(description='UniProt + CDD virus contig filter')
    parser.add_argument('-i', '--input', help='Input fasta (use with -o for single sample)')
    parser.add_argument('--input-dir', help='Input directory of per-sample candidates (merges automatically)')
    parser.add_argument('-o', '--outdir', required=True)
    parser.add_argument('-m', '--filter-mode', choices=['filter','strict','raw'], default='filter')
    parser.add_argument('--skip-uniprot', action='store_true')
    parser.add_argument('--skip-cdd', action='store_true')
    parser.add_argument('--uniprot-db', default=DEFAULTS['UNIPROT_DB'])
    parser.add_argument('--virus-taxid', default=DEFAULTS['VIRUS_TAXID'])
    parser.add_argument('--cdd-db', default=DEFAULTS['CDD_DB'])
    parser.add_argument('--cdd-whitelist', default=DEFAULTS['CDD_WHITELIST'])
    parser.add_argument('--cdd-classified', default=DEFAULTS['CDD_CLASSIFIED'])
    parser.add_argument('-t', '--threads', default=DEFAULTS['THREADS'])
    parser.add_argument('--split-by-prefix', action='store_true',
                        help='Split output by contig prefix (first _ segment) into sample dirs')
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
        if not args.skip_cdd:
            args.split_by_prefix = True  # auto-split back

    if not args.input:
        print('ERROR: -i or --input-dir required')
        sys.exit(1)

    # Mode resolution
    skip_uniprot_explicit = '--skip-uniprot' in sys.argv
    skip_cdd_explicit = '--skip-cdd' in sys.argv

    if skip_uniprot_explicit or skip_cdd_explicit:
        if not args.skip_uniprot:
            args.uniprot_mode = 'filter'
        if not args.skip_cdd:
            args.mode = 'filter'
    else:
        if args.filter_mode == 'raw':
            args.skip_uniprot = True
            args.skip_cdd = True
        elif args.filter_mode == 'filter':
            args.skip_uniprot = False
            args.skip_cdd = False
            args.uniprot_mode = 'filter'
            args.mode = 'filter'
        elif args.filter_mode == 'strict':
            args.skip_uniprot = False
            args.skip_cdd = False
            args.uniprot_mode = 'strict'
            args.mode = 'strict'

    print(f'Input:  {args.input}')
    print(f'Output: {args.outdir}')
    print(f'Mode:   {args.filter_mode}', flush=True)

    current_input = args.input

    # ---- CDD first (faster, reduces input for UniProt) ----
    cdd_dir = os.path.join(args.outdir, 'cdd')
    os.makedirs(cdd_dir, exist_ok=True)
    filter_fa_path = None
    calls_tsv_path = None
    base = os.path.basename(current_input).rsplit('.', 1)[0]

    if args.skip_cdd:
        filter_fa = os.path.join(cdd_dir, f'{base}_cdd_filter.fasta')
        strict_fa = os.path.join(cdd_dir, f'{base}_cdd_strict.fasta')
        symlink_force(current_input, filter_fa)
        symlink_force(current_input, strict_fa)
        print('CDD: SKIP -> symlink')
    else:
        filter_fa_path = os.path.join(cdd_dir, f'{base}_cdd_filter.fasta')
        calls_tsv_path = os.path.join(cdd_dir, 'cdd_calls.tsv')
        if args.mode in ('filter', 'both'):
            if os.path.isfile(filter_fa_path) and os.path.getsize(filter_fa_path) > 0:
                print(f'CDD: [resume] 已有 {filter_fa_path}, 跳过 CDD 搜索')
            else:
                run_cdd(current_input, cdd_dir, args.cdd_db,
                        args.cdd_whitelist, args.cdd_classified, args.threads, 'filter', orig_base=base)
        if args.mode in ('strict', 'both'):
            strict_fa = os.path.join(cdd_dir, f'{base}_cdd_strict.fasta')
            if os.path.isfile(strict_fa) and os.path.getsize(strict_fa) > 0:
                print(f'CDD: [resume] 已有 {strict_fa}, 跳过 strict CDD')
            else:
                run_cdd(current_input, cdd_dir, args.cdd_db,
                        args.cdd_whitelist, args.cdd_classified, args.threads, 'strict', orig_base=base)

    # ---- UniProt second (only on CDD-filtered contigs, much smaller) ----
    uniprot_dir = os.path.join(args.outdir, 'uniprot')
    os.makedirs(uniprot_dir, exist_ok=True)

    if args.skip_uniprot:
        uniprot_out = filter_fa_path or current_input
        uniprot_out_link = os.path.join(uniprot_dir, f'{base}_uniprot.fasta')
        symlink_force(uniprot_out, uniprot_out_link)
        print('UniProt: SKIP -> symlink')
    else:
        uni_input = filter_fa_path or current_input
        uniprot_out = run_uniprot(uni_input, uniprot_dir, args.uniprot_db,
                                  args.threads, args.virus_taxid, args.uniprot_mode)

    # ---- Output ----
    # 确定实际运行的组件 (供 split 判断)
    cdd_real = filter_fa_path if filter_fa_path else None          # CDD 真跑过
    uni_real = uniprot_out if not args.skip_uniprot else None      # UniProt 真跑过

    if args.split_by_prefix:
        split_output(args.outdir, cdd_real, uni_real,
                     cdd_dir, uniprot_dir, calls_tsv_path, base, mode=args.mode)
    else:
        final_src = uni_real or cdd_real or current_input
        final_dst = os.path.join(args.outdir, f'{base}_filtered.fasta')
        symlink_force(final_src, final_dst)
        cobra_dst = os.path.join(args.outdir, f'{base}_virus.all.candidate.fasta')
        symlink_force(final_src, cobra_dst)
        print(f'Final: {final_dst} -> {final_src}')

if __name__ == '__main__':
    main()
