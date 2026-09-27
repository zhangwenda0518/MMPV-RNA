#!/usr/bin/env python3
"""Batch viroid BLASTN best-hit analysis: all + per-virus + per-sample + per-species."""
import csv
from collections import defaultdict
from pathlib import Path


def load_taxonomy(tsv_path):
    tax = {}
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            tax[row['accession']] = dict(row)
    return tax


def load_fasta_lengths(fasta_path):
    lens = {}
    if not Path(fasta_path).is_file():
        return lens
    with open(fasta_path) as f:
        sid, seq = '', ''
        for line in f:
            if line.startswith('>'):
                if seq:
                    lens[sid.split()[0]] = len(seq)
                sid = line[1:].strip().split()[0]; seq = ''
            else:
                seq += line.strip()
        if seq:
            lens[sid.split()[0]] = len(seq)
    return lens


def process_sample(sample_dir, tax, ref_lens):
    sample = sample_dir.parent.name
    best_files = sorted(sample_dir.glob('*viroids.blastn.result.txt.best'))
    if not best_files:
        return []
    best = best_files[0]
    results = []
    with open(best) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue
            qid, sid = cols[0], cols[1]
            pid, alen = float(cols[2]), int(cols[3])
            ev, bs = cols[10], float(cols[11])
            t = tax.get(sid, {})
            results.append({
                'sample': sample, 'qid': qid, 'sid': sid,
                'pid': pid, 'alen': alen, 'evalue': ev, 'bitscore': bs,
                'rlen': ref_lens.get(sid, 0),
                'species': t.get('species', ''),
                'genus': t.get('genus', ''),
                'family': t.get('family', ''),
                'description': t.get('description', ''),
            })
    return results


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--ident-dir', required=True, help='02_Identification dir')
    p.add_argument('--taxonomy', required=True)
    p.add_argument('--ref-fasta', required=True)
    p.add_argument('-o', '--output', default='./viroid_analysis')
    args = p.parse_args()
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    tax = load_taxonomy(args.taxonomy)
    ref_lens = load_fasta_lengths(args.ref_fasta)

    all_results = []
    total_samples = 0
    for sd in sorted(Path(args.ident_dir).iterdir()):
        vd = sd / 'viroid_output'
        if not vd.is_dir():
            continue
        total_samples += 1
        r = process_sample(vd, tax, ref_lens)
        all_results.extend(r)
        if len(r) > 0:
            print('  {}: {} contigs'.format(sd.name, len(r)))
    print('Samples with viroid output: {}'.format(total_samples))
    print('Total viroid contigs: {}'.format(len(all_results)))
    if not all_results:
        print('No viroid hits found')
        return

    # ── 1. all_viroids.tsv (master table) ──
    with open(out / 'Viroid.all_info.tsv', 'w') as tf:
        tf.write('sample\tcontig_id\tref_acc\tref_len\tpident\talen\tevalue\tbitscore\tfamily\tgenus\tspecies\tdescription\n')
        for r in sorted(all_results, key=lambda x: -x['bitscore']):
            tf.write('{sample}\t{qid}\t{sid}\t{rlen}\t{pid:.1f}\t{alen}\t{evalue}\t{bitscore:.0f}\t{family}\t{genus}\t{species}\t{description}\n'.format(**r))

    # ── 2. viroid_species_summary.tsv ──
    sp_stats = defaultdict(lambda: {'n_contigs': 0, 'n_samples': set(), 'sum_pid': 0.0, 'rlen': 0})
    for r in all_results:
        sp = r['species'] or 'Unclassified'
        sp_stats[sp]['n_contigs'] += 1
        sp_stats[sp]['n_samples'].add(r['sample'])
        sp_stats[sp]['sum_pid'] += r['pid']
        if r['rlen'] > 0:
            sp_stats[sp]['rlen'] = r['rlen']

    with open(out / 'Viroid.species_info.tsv', 'w') as sf:
        sf.write('species\tgenus\tfamily\tref_len\tn_contigs\tn_samples\tavg_pident\n')
        for sp, info in sorted(sp_stats.items(), key=lambda x: -x[1]['n_contigs']):
            ge = ''; fa = ''
            for r in all_results:
                if (r['species'] or 'Unclassified') == sp:
                    ge = r['genus']; fa = r['family']; break
            avg_pid = info['sum_pid'] / info['n_contigs']
            sf.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.1f}\n'.format(
                sp, ge, fa, info['rlen'], info['n_contigs'], len(info['n_samples']), avg_pid))

    # ── 3. per-virus breakdown (by ref accession) ──
    virus_stats = defaultdict(lambda: {'n_contigs': 0, 'n_samples': set(), 'sum_pid': 0.0, 'sum_alen': 0, 'rlen': 0})
    for r in all_results:
        vid = r['sid']
        virus_stats[vid]['n_contigs'] += 1
        virus_stats[vid]['n_samples'].add(r['sample'])
        virus_stats[vid]['sum_pid'] += r['pid']
        virus_stats[vid]['sum_alen'] += r['alen']
        if r['rlen'] > 0:
            virus_stats[vid]['rlen'] = r['rlen']

    with open(out / 'Viroid.per_virus.tsv', 'w') as vf:
        vf.write('ref_acc\tspecies\tgenus\tfamily\tref_len\tn_contigs\tn_samples\tavg_pident\tavg_alen\tcoverage_pct\n')
        for vid, info in sorted(virus_stats.items(), key=lambda x: -x[1]['n_contigs']):
            sp = ''; ge = ''; fa = ''
            for r in all_results:
                if r['sid'] == vid:
                    sp = r['species']; ge = r['genus']; fa = r['family']; break
            avg_pid = info['sum_pid'] / info['n_contigs']
            avg_alen = info['sum_alen'] / info['n_contigs']
            cov_pct = 100 * sum(r['alen'] for r in all_results if r['sid'] == vid) / info['rlen'] if info['rlen'] else 0
            vf.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.0f}\t{:.1f}\n'.format(
                vid, sp, ge, fa, info['rlen'], info['n_contigs'], len(info['n_samples']),
                avg_pid, avg_alen, cov_pct))

    # ── 4. viroid contig FASTA ──
    viroid_ids = set(r['qid'] for r in all_results)
    from Bio import SeqIO
    n_written = 0
    with open(out / 'Viroid.contigs.fasta', 'w') as vf_out:
        for sd in sorted(Path(args.ident_dir).iterdir()):
            vd = sd / 'viroid_output'
            if not vd.is_dir():
                continue
            # look up taxonomy info for this sample's contigs
            contig_info = {}
            for r in all_results:
                if r['sample'] == sd.name:
                    contig_info[r['qid']] = '{}|{}|{}'.format(
                        r.get('species', '').replace(' ', '_'),
                        r.get('genus', '').replace(' ', '_'),
                        r.get('family', '').replace(' ', '_'))
            # find virus candidate FASTA
            for cand_name in ['{}.virus.candidate.fasta'.format(sd.name), 'virus.candidate.fasta', '{}_virus.candidate.fasta'.format(sd.name)]:
                cand_fa = sd / cand_name
                if cand_fa.is_file():
                    break
            if not cand_fa.is_file():
                continue
            for rec in SeqIO.parse(str(cand_fa), 'fasta'):
                if rec.id.split()[0] in viroid_ids:
                    info = contig_info.get(rec.id.split()[0], '')
                    if info:
                        rec.description = info
                    SeqIO.write(rec, vf_out, 'fasta')
                    n_written += 1
    print('Viroid contig sequences: {} wrote to viroid_contigs.fasta'.format(n_written))

    # ── 5. per-sample breakdown ──
    sample_stats = defaultdict(lambda: {'n_contigs': 0, 'species': set(), 'top_species': '', 'top_count': 0})
    for r in all_results:
        s = r['sample']
        sample_stats[s]['n_contigs'] += 1
        sample_stats[s]['species'].add(r['species'] or 'Unclassified')
    # determine top species per sample
    for s in sample_stats:
        sp_count = defaultdict(int)
        for r in all_results:
            if r['sample'] == s:
                sp_count[r['species'] or 'Unclassified'] += 1
        top_sp = max(sp_count, key=sp_count.get) if sp_count else ''
        sample_stats[s]['top_species'] = top_sp
        sample_stats[s]['top_count'] = sp_count[top_sp]

    with open(out / 'Viroid.per_sample.tsv', 'w') as smf:
        smf.write('sample\tn_viroid_contigs\tn_species\ttop_species\ttop_count\n')
        for s, info in sorted(sample_stats.items()):
            smf.write('{}\t{}\t{}\t{}\t{}\n'.format(
                s, info['n_contigs'], len(info['species']), info['top_species'], info['top_count']))

    # ── Terminal summary ──
    fams = defaultdict(int)
    gens = defaultdict(int)
    for r in all_results:
        fams[r['family'] or 'Unclassified'] += 1
        gens[r['genus'] or 'Unclassified'] += 1
    print('\nFamily distribution:')
    for fa, n in sorted(fams.items(), key=lambda x: -x[1]):
        print('  {}: {}'.format(fa, n))
    print('\nTop genera:')
    for ge, n in sorted(gens.items(), key=lambda x: -x[1])[:15]:
        print('  {}: {}'.format(ge, n))
    print('\nTop species:')
    for sp, info in sorted(sp_stats.items(), key=lambda x: -x[1]['n_contigs'])[:15]:
        print('  {}: ref={}bp, {} contigs, {} samples, avg_pid={:.1f}%'.format(
            sp, info['rlen'], info['n_contigs'], len(info['n_samples']),
            info['sum_pid'] / info['n_contigs']))
    print('\nOutput: {}'.format(out))


if __name__ == '__main__':
    main()
