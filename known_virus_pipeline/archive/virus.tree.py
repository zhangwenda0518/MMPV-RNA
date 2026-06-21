#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ACVirus Pangenomic & Evolutionary Analysis Pipeline
===================================================
Features upgraded (Latest):
1. Cladogram toggle: `--ignore_branch_length` for topological flattening.
2. Lineage strict deduplication: Auto CD-HIT redundant sequence removal.
3. Advanced Right-Aligned Tree Render (Nature/Cell Style Tips)
4. Centralized Label Tanglegram plotting
5. Diamond/BLASTP & Smart ModelTest-NG -> RAxML-NG -> IQ-TREE Engine
"""

import argparse
import os
import re
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ==========================================
# 模块 1：GenBank 解析器
# ==========================================
def parse_gb_to_baits(gb_file: Path, outdir: Path):
    baits = []
    print(f"[INFO] Parsing GenBank file: {gb_file.name}")
    for record in SeqIO.parse(gb_file, "genbank"):
        for feat in record.features:
            if feat.type == "CDS":
                gene_name = feat.qualifiers.get("gene", feat.qualifiers.get("product", ["unknown"]))[0]
                gene_name = "".join(c if c.isalnum() else "_" for c in gene_name)
                translation = feat.qualifiers.get("translation", [None])[0]
                if translation:
                    aa_record = SeqRecord(Seq(translation), id=gene_name, description=record.id)
                    bait_path = outdir / f"bait_{gene_name}.faa"
                    SeqIO.write([aa_record], bait_path, "fasta")
                    baits.append({'gene': gene_name, 'bait_path': bait_path})
    print(f"[INFO] Successfully extracted {len(baits)} genes from GenBank.")
    return baits


# ==========================================
# 模块 2：智能多模态抽样、智能长度聚类(CD-HIT)
# ==========================================
def get_sampled_accessions(df, target_genus, mode, target_species=None, add_outgroup=False, specific_outgroup=None, target_ratio=1.0, other_count=5, outgroup_count=3):
    target_info = df[df['Genus'].str.lower() == target_genus.lower()]
    if target_info.empty: raise ValueError(f"[ERROR] Cannot find Genus '{target_genus}' in taxa.txt!")
    t_family, t_order = target_info.iloc[0]['Family'], target_info.iloc[0]['Order']
    if pd.isna(t_order) or str(t_order).strip() == '': t_order, outgroup_level = target_info.iloc[0]['Class'], 'Class'
    else: outgroup_level = 'Order'

    sampled_accessions = set()
    
    if mode == 'phylogeny':
        print(f"[INFO] Mode: Phylogeny (Family-level). Target: {t_family} / {target_genus}")
        fam_df = df[df['Family'] == t_family]
        for genus, group in sorted(fam_df.groupby('Genus')):
            if str(genus).lower() == target_genus.lower():
                n_sample = max(1, int(len(group) * target_ratio))
                sampled = group.sample(n=n_sample)
                if target_ratio == 1.0: print(f"  - Keeping ALL {len(sampled)} sequences from target Genus: {genus}")
                else: print(f"  - Target Genus {genus}: Sampled {n_sample}/{len(group)} sequences")
            else:
                n_sample = min(len(group), other_count)
                sampled = group.sample(n=n_sample)
                print(f"  - Sampled {len(sampled)} sequences from Genus: {genus}")
            sampled_accessions.update(sampled['Virus GENBANK accession'].dropna().tolist())
            
        out_df = df[(df[outgroup_level] == t_order) & (df['Family'] != t_family)]
        out_families = out_df['Family'].unique()
        if len(out_families) > 0:
            for family, group in out_df.groupby('Family'):
                sampled_accessions.update(group.sample(n=min(len(group), outgroup_count))['Virus GENBANK accession'].dropna().tolist())
            print(f"  - Sampled outgroups from {len(out_families)} other families in {t_order} (up to {outgroup_count} per family)")
        else:
            print(f"  - [WARNING] No outgroup found in {outgroup_level} '{t_order}'.")
                
    elif mode == 'simple':
        print(f"[INFO] Mode: Simple (Genus-level). Target: {target_genus}")
        target_group = df[df['Genus'].str.lower() == target_genus.lower()]
        sampled_accessions.update(target_group['Virus GENBANK accession'].dropna().tolist())
        print(f"  - Keeping ALL {len(target_group)} sequences from target Genus: {target_genus}")
        
        out_df = df[(df['Family'] == t_family) & (df['Genus'].str.lower() != target_genus.lower())]
        if not out_df.empty:
            out_sample = out_df.sample(n=min(len(out_df), outgroup_count))
            sampled_accessions.update(out_sample['Virus GENBANK accession'].dropna().tolist())
            print(f"  - Sampled {len(out_sample)} sequences from other genera in family '{t_family}' as outgroups.")
            
    elif mode == 'lineage':
        if not target_species: raise ValueError("[ERROR] For 'lineage' mode, specify --species")
        print(f"[INFO] Mode: Lineage (Species-level). Target Species: {target_species}")
        target_group = df[df['Species'].str.lower() == target_species.lower()]
        if not target_group.empty: 
            sampled_accessions.update(target_group['Virus GENBANK accession'].dropna().tolist())
            print(f"  - Keeping ALL {len(target_group)} reference strains of species '{target_species}'.")

        if specific_outgroup:
            out_group = df[df['Species'].str.lower() == specific_outgroup.lower()]
            if not out_group.empty: 
                sampled_accessions.update(out_group.sample(n=1)['Virus GENBANK accession'].dropna().tolist())
                print(f"  - Using SPECIFIED outgroup: {specific_outgroup} (1 sequence).")
        elif add_outgroup:
            out_df = df[(df['Genus'].str.lower() == target_genus.lower()) & (df['Species'].str.lower() != target_species.lower())]
            if not out_df.empty: 
                sampled_accessions.update(out_df.sample(n=1)['Virus GENBANK accession'].dropna().tolist())
                print(f"  - AUTO sampled 1 sequence from a different species in '{target_genus}' as outgroup.")

    return [str(acc).split(';')[0].strip() for acc in sampled_accessions if pd.notna(acc)]

def sample_user_fasta(input_fasta, output_fasta, mode, max_n_ratio=0.05):
    valid_recs, dropped_recs = [], []
    for rec in SeqIO.parse(input_fasta, "fasta"):
        seq_str = str(rec.seq).upper()
        if not seq_str: continue
        n_ratio = seq_str.count('N') / len(seq_str)
        if mode == 'lineage' and n_ratio > max_n_ratio: dropped_recs.append((rec.id, n_ratio))
        else: valid_recs.append(rec)
        
    if dropped_recs:
        print(f"[WARNING] Dropped {len(dropped_recs)} user sequences due to excessive 'N's (> {max_n_ratio*100}%):")
        for drop_id, r in dropped_recs: print(f"    - {drop_id} (N ratio: {r:.2%})")
        
    print(f"[INFO] Kept {len(valid_recs)} valid user sequences.")
    SeqIO.write(valid_recs, output_fasta, "fasta")

def extract_and_merge(db_accs, db_fasta, user_sampled_fasta, final_fasta, threads):
    outdir = Path(final_fasta).parent
    db_ids_file, temp_db_seqs = outdir / "temp_db_ids.txt", outdir / "temp_db_seqs.fasta"
    
    with open(db_ids_file, "w") as f:
        for acc in db_accs: f.write(f"^{acc}(\\.[0-9]+)?\n")
        
    if db_accs: 
        subprocess.run(['seqkit', 'grep', '-j', str(threads), '-r', '-f', str(db_ids_file), str(db_fasta), '-o', str(temp_db_seqs)], check=True)
    else: 
        Path(temp_db_seqs).touch()
        
    with open(final_fasta, 'wb') as outfile:
        if user_sampled_fasta and Path(user_sampled_fasta).exists():
            with open(user_sampled_fasta, 'rb') as infile: outfile.write(infile.read())
        if Path(temp_db_seqs).exists():
            with open(temp_db_seqs, 'rb') as infile: outfile.write(infile.read())
            
    subprocess.run(['seqkit', 'rmdup', '-j', str(threads), '-s', '-i', str(final_fasta), '-o', str(final_fasta)+".tmp"], check=True, capture_output=True)
    os.replace(str(final_fasta)+".tmp", str(final_fasta))
    
    for tmp in [db_ids_file, temp_db_seqs, user_sampled_fasta]:
        if tmp and isinstance(tmp, Path) and tmp.exists(): tmp.unlink()

# --- [新增] CD-HIT 高级冗余去除引擎 ---
def calculate_adaptive_threshold(ref_length):
    if ref_length == 0: return 0.99
    if ref_length < 1000: return max(0.85, 1.0 - (5.0 / ref_length)) 
    elif ref_length < 10000: return max(0.90, 1.0 - (15.0 / ref_length))
    else: return max(0.95, 1.0 - (30.0 / ref_length))

def run_cdhit_and_parse(in_fasta: Path, out_fasta: Path, report_tsv: Path, threshold: float, threads: int):
    """运用 CD-HIT 生成聚类并且输出只包含代表序列的 Fasta 及 包含折叠信息的分析报告"""
    clstr_file = Path(str(out_fasta) + ".clstr")
    cmd = [
        "cd-hit-est", "-i", str(in_fasta), "-o", str(out_fasta),
        "-c", str(threshold), "-aL", "0.95", "-aS", "0.85",
        "-d", "0", "-M", "0", "-T", str(threads)
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        print(f"[ERROR] CD-HIT execution failed! Is cd-hit-est installed? Base Error: {e}")
        raise

    clusters = []
    with open(clstr_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                clusters.append({'id': line.strip().replace('>', ''), 'rep': None, 'rep_size': '', 'members': []})
            else:
                parts = line.strip().split(', >')
                str_size = parts[0].split('\t')[1]
                rest = parts[1].split('... ')
                samp_id = rest[0].split()[0]
                ident = rest[1].strip()
                if ident == '*':
                    clusters[-1]['rep'] = samp_id
                    clusters[-1]['rep_size'] = str_size
                    clusters[-1]['members'].append((samp_id, str_size, '100% (Ref)'))
                else:
                    m = re.search(r'([0-9\.]+\%)', ident)
                    ident_val = m.group(1) if m else ident
                    clusters[-1]['members'].append((samp_id, str_size, ident_val))
    
    with open(report_tsv, 'w', encoding='utf-8') as f:
        f.write(f"# Target CD-HIT Auto Threshold: {threshold*100:.2f}%\n")
        f.write("Cluster_ID\tRepresentative_Sample\tRep_Length\tCluster_Size\tCluster_Members_Details\n")
        for c in clusters:
            members_str = ", ".join([f"{m[0]}({m[2]})" for m in c['members'] if m[0] != c['rep']])
            if not members_str: members_str = "None"
            f.write(f"{c['id']}\t{c['rep']}\t{c['rep_size']}\t{len(c['members'])}\t{members_str}\n")
            
    if clstr_file.exists(): clstr_file.unlink()

# ==========================================
# 模块 3：同源靶标提取与共线性数据准备
# ==========================================
def predict_and_extract_markers(combined_fasta: Path, bait_protein: Path, out_faa: Path, out_fna: Path, threads: int, resume: bool, tool: str, evalue: str, max_target_seqs: int):
    if resume and out_faa.exists() and out_fna.exists() and out_faa.stat().st_size > 0: return True, []
    outdir, all_faa, all_fna = combined_fasta.parent, combined_fasta.parent/"temp_all_predicted.faa", combined_fasta.parent/"temp_all_predicted.fna"
    blast_tsv = outdir / "temp_marker_blast.tsv"
    
    if not (all_faa.exists() and all_fna.exists()):
        print("    -> Predicting ORFs (AA & NT) for all genomes...")
        subprocess.check_call(f'prodigal -i {combined_fasta} -a {all_faa} -d {all_fna} -p meta -q 1>/dev/null', shell=True)
        
    if tool == 'diamond':
        dmnd_db = outdir / "temp_marker_db.dmnd"
        print(f"    -> Mapping bait to genomes via Diamond (Ultra-Sensitive mode)...")
        if not dmnd_db.exists(): subprocess.check_call(f'diamond makedb --in {all_faa} -d {dmnd_db} -p {threads} --quiet', shell=True)
        subprocess.check_call(f'diamond blastp --ultra-sensitive -q {bait_protein} -d {dmnd_db} -p {threads} -o {blast_tsv} -k {max_target_seqs} --outfmt 6 qseqid sseqid pident length evalue --evalue {evalue} --quiet', shell=True)
    elif tool == 'blastp':
        blast_db = outdir / "temp_marker_blastdb"
        print(f"    -> Mapping bait to genomes via NCBI BLASTP (Deep Evolutionary Search)...")
        if not Path(str(blast_db)+".phr").exists(): subprocess.check_call(f'makeblastdb -in {all_faa} -dbtype prot -out {blast_db} -quiet', shell=True)
        subprocess.check_call(f'blastp -query {bait_protein} -db {blast_db} -num_threads {threads} -out {blast_tsv} -max_target_seqs {max_target_seqs} -outfmt "6 qseqid sseqid pident length evalue" -evalue {evalue}', shell=True)
    
    genome_best_hits = {}
    with open(blast_tsv, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            sseqid, ev = parts[1], float(parts[4])
            genome_id = sseqid.rsplit('_', 1)[0]
            if genome_id not in genome_best_hits or ev < genome_best_hits[genome_id]['evalue']:
                genome_best_hits[genome_id] = {'orf_id': sseqid, 'evalue': ev}
                
    orf_to_keep = {v['orf_id'] for v in genome_best_hits.values()}
    print(f"    -> Found homologous sequences in {len(orf_to_keep)} genomes.")

    if len(orf_to_keep) < 3: 
        print(f"    -> [WARNING] Too few sequences found ({len(orf_to_keep)}). Skipping gene.")
        return False, []

    if not (resume and out_faa.exists() and out_fna.exists() and out_faa.stat().st_size > 0):
        def extract_and_rename(in_file, out_file, desc):
            recs = [rec for rec in SeqIO.parse(in_file, "fasta") if rec.id in orf_to_keep]
            for r in recs: r.id, r.description = r.id.rsplit('_', 1)[0], desc
            SeqIO.write(recs, out_file, "fasta")
        extract_and_rename(all_faa, out_faa, "Marker_AA")
        extract_and_rename(all_fna, out_fna, "Marker_NT")
    
    if blast_tsv.exists(): blast_tsv.unlink()
    return True, list(orf_to_keep)

def process_prodigal_gff(input_file, output_file):
    seq_lengths, protein_counters, results = {}, {}, []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Sequence Data'):
                seq_id = [p for p in line.split(';') if 'seqhdr' in p][0].split('"')[1].split()[0]
                seq_lengths[seq_id] = int([p for p in line.split(';') if 'seqlen' in p][0].split('=')[1])
            elif not line.startswith('#') and line:
                fields = line.split('\t')
                if len(fields) >= 9 and fields[2] == 'CDS':
                    seq_id, start, end = fields[0], int(fields[3]), int(fields[4])
                    if start > end: start, end = end, start
                    protein_counters[seq_id] = protein_counters.get(seq_id, 0) + 1
                    results.append([seq_id, f"{seq_id}_{protein_counters[seq_id]}", start, end, seq_lengths.get(seq_id, 0)])
    pd.DataFrame(results, columns=['nucl_id', 'protein', 'start', 'end', 'nucl_length']).to_csv(output_file, sep='\t', index=False)

def prep_synteny_data(infile: Path, outdir: Path, threads: int, resume: bool):
    prefix = infile.stem
    gff, faa, dmnd, tab, tsv = [outdir / f"{prefix}{ext}" for ext in ['.gff', '_prot.faa', '.dmnd', '.tab', '.tsv']]
    if not (resume and tab.exists() and tab.stat().st_size > 0 and gff.exists()):
        subprocess.check_call(f'prodigal -i {infile} -a {faa} -f gff -p meta -o {gff} -q 1>/dev/null', shell=True)
        subprocess.check_call(f'diamond makedb --in {faa} -d {dmnd} -p {threads} --quiet', shell=True)
        subprocess.check_call(f'diamond blastp --more-sensitive -q {faa} -d {dmnd} -p {threads} -o {tab} --evalue 1e-3 --quiet', shell=True)
    if not tsv.exists() or not resume:
        process_prodigal_gff(gff, tsv)
        
    pos_df = pd.read_csv(tsv, sep='\t')
    blast_df = pd.read_csv(tab, sep='\t', header=None, usecols=[0,1,2], names=['protein1','protein2','similarity'])
    sim_df = blast_df[blast_df['protein1'] != blast_df['protein2']]
    return pos_df, sim_df


# ==========================================
# 模块 4：原生外群定根与智能建树引擎
# ==========================================
def find_outgroup_list(aln_file: Path, target_genus: str, taxa_meta: dict, mode: str, target_species: str = None, specific_outgroup: str = None):
    seq_ids = [rec.id for rec in SeqIO.parse(aln_file, "fasta")]
    target_family = None
    for acc, meta in taxa_meta.items():
        if str(meta.get('Genus')).lower() == target_genus.lower():
            target_family = meta.get('Family')
            break

    outgroup_nodes = []
    for seq_id in seq_ids:
        base_acc = seq_id.split('.')[0]
        meta = taxa_meta.get(base_acc, {})
        genus, family, species = meta.get('Genus', 'Unknown'), meta.get('Family', 'Unknown'), meta.get('Species', 'Unknown')
        if genus == 'Unknown' or family == 'Unknown': continue

        if mode == 'phylogeny' and family != target_family: outgroup_nodes.append(seq_id)
        elif mode == 'simple' and genus.lower() != target_genus.lower(): outgroup_nodes.append(seq_id)
        elif mode == 'lineage':
            if specific_outgroup and species.lower() == specific_outgroup.lower(): outgroup_nodes.append(seq_id)
            elif not specific_outgroup and species.lower() != target_species.lower() and genus.lower() == target_genus.lower(): outgroup_nodes.append(seq_id)

    return ",".join(outgroup_nodes) if outgroup_nodes else None

def parse_best_model(modeltest_out_file, criterion="BIC"):
    if not os.path.isfile(modeltest_out_file): return None
    best_model = None
    found_target = False
    with open(modeltest_out_file, 'r') as f:
        for line in f:
            clean_line = line.strip()
            if clean_line.startswith(f"Best model according to {criterion}"): found_target = True
            if found_target and clean_line.startswith("Model:"):
                best_model = clean_line.split()[1]
                break
    return best_model

def run_smart_tree(aln_file: Path, seq_type: str, threads: int, resume: bool, run_mt: bool, criterion: str, outgroups_str: str):
    mt_seq_type = 'nt' if 'Nucleotide' in seq_type or seq_type == 'NT' else 'aa'
    iqtree_out, raxml_out, raxml_best = Path(str(aln_file) + '.treefile'), Path(str(aln_file) + '.raxml.support'), Path(str(aln_file) + '.raxml.bestTree')

    if resume and (raxml_out.exists() or raxml_best.exists() or iqtree_out.exists()): 
        return raxml_out if raxml_out.exists() else (raxml_best if raxml_best.exists() else iqtree_out)

    use_iqtree_fallback = not run_mt

    if run_mt:
        print(f"    -> [ModelTest-NG] Selecting best model ({criterion}) for {seq_type}...")
        try:
            subprocess.run(f"modeltest-ng -i {aln_file} -d {mt_seq_type} -p {threads} -T raxml", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            best_model = parse_best_model(f"{aln_file}.out", criterion)
            
            if best_model:
                print(f"    -> [RAxML-NG] Building tree using model: {best_model}...")
                out_param = f"--outgroup {outgroups_str}" if outgroups_str else ""
                
                if outgroups_str: print(f"      -> [Rooting] Instructing RAxML-NG to root with outgroups: {outgroups_str}")
                subprocess.run(f"raxml-ng --all --msa {aln_file} --model {best_model} {out_param} --threads {threads} --bs-trees 1000 --redo", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                return raxml_out if raxml_out.exists() else raxml_best
            else:
                print("    -> [WARNING] ModelTest-NG parsing failed! Falling back to IQ-TREE.")
                use_iqtree_fallback = True
        except Exception:
            print("    -> [WARNING] ModelTest-NG/RAxML-NG execution crashed or not found! Falling back to IQ-TREE.")
            use_iqtree_fallback = True

    if use_iqtree_fallback:
        print(f"    -> [IQ-TREE] Building {seq_type} tree with MFP...")
        redo_flag = "" if resume else "-redo"
        out_param = f"-o {outgroups_str}" if outgroups_str else ""
        if outgroups_str: print(f"      -> [Rooting] Instructing IQ-TREE to root with outgroups: {outgroups_str}")
        subprocess.check_call(f'iqtree -s {aln_file} -m MFP {out_param} -B 1000 -T {threads} {redo_flag} -quiet', shell=True)
        return iqtree_out

def build_single_tree(fasta_file: Path, out_prefix: Path, seq_type: str, threads: int, resume: bool, run_mt: bool, criterion: str, target_genus: str, taxa_meta: dict, mode: str, target_species: str = None, specific_outgroup: str = None):
    mafft_out = Path(str(out_prefix) + '.mafft.fasta')
    if not (resume and mafft_out.exists() and mafft_out.stat().st_size > 0):
        print(f"    -> [MAFFT] Aligning {seq_type} sequences...")
        subprocess.check_call(f'mafft --auto --thread {threads} {fasta_file} > {mafft_out} 2>/dev/null', shell=True)
        
    outgroups_str = find_outgroup_list(mafft_out, target_genus, taxa_meta, mode, target_species, specific_outgroup)
    treefile = run_smart_tree(mafft_out, seq_type, threads, resume, run_mt, criterion, outgroups_str)
    return treefile, mafft_out


# ==========================================
# 模块 5：核心渲染引擎与排版修正
# ==========================================
# --- [新增] 枝长动态压平功能 ---
def flatten_branch_lengths(tree, default_length=1.0):
    """将整个拓扑树压平以便进行类似于 cladogram 的展示，极大地减少不均匀分支带来的文本混乱"""
    for clade in tree.find_clades():
        clade.branch_length = default_length

def draw_tree_layer(ax_tree, tree, terminals, leaf_y, taxa_meta, max_depth):
    def calc_coords(node, current_x):
        if node.is_terminal(): node.x, node.y = current_x + (node.branch_length or 0.0), leaf_y[node.name]; return node.y
        else:
            node.x = current_x + (node.branch_length or 0.0)
            node.y = sum([calc_coords(c, node.x) for c in node.clades]) / len(node.clades); return node.y
    calc_coords(tree.root, 0.0)

    align_x = max_depth * 1.15
    strip_x = max_depth * 3.5  

    blocks, curr_g, start_y = [], None, None
    for n in sorted(terminals, key=lambda n: n.y):
        g = taxa_meta.get(n.name.split('.')[0], {}).get('Genus', 'Unknown')
        if pd.isna(g): g = 'Unknown'
        if g != curr_g:
            if curr_g and curr_g != 'Unknown': blocks.append((curr_g, start_y, n.y - 1))
            curr_g, start_y = g, n.y
    if curr_g and curr_g != 'Unknown': blocks.append((curr_g, start_y, sorted(terminals, key=lambda n: n.y)[-1].y))

    unique_genera = list(set([b[0] for b in blocks]))
    cmap_bg = plt.get_cmap('Pastel1')
    gen_colors = {g: cmap_bg(i % 9) for i, g in enumerate(unique_genera)}

    for g, sy, ey in blocks:
        rect = patches.Rectangle((0, sy - 0.5), align_x, (ey - sy) + 1, facecolor=gen_colors[g], alpha=0.3, lw=0, zorder=0)
        ax_tree.add_patch(rect)

    def draw_clade(node):
        if not node.is_terminal():
            min_y, max_y = min(c.y for c in node.clades), max(c.y for c in node.clades)
            ax_tree.plot([node.x, node.x], [min_y, max_y], color='black', lw=1.2, zorder=2)
            
            support = node.confidence if hasattr(node, 'confidence') and node.confidence is not None else getattr(node, 'name', None)
            if support and node != tree.root:
                try:
                    val = float(str(support).split('/')[0])
                    if val >= 50: ax_tree.text(node.x - max_depth*0.01, node.y - 0.1, str(support), color='blue', fontsize=7, ha='right', va='bottom', zorder=4)
                except ValueError: pass

            for c in node.clades: 
                ax_tree.plot([node.x, c.x], [c.y, c.y], color='black', lw=1.2, zorder=2)
                draw_clade(c)
    draw_clade(tree.root)
    
    for node in terminals:
        base_acc = node.name.split('.')[0]
        sp = taxa_meta.get(base_acc, {}).get('Species', 'Unknown')
        g = taxa_meta.get(base_acc, {}).get('Genus', 'Unknown')
        label = f"{sp} ({node.name})" if sp != 'Unknown' and not pd.isna(sp) else node.name
        
        ax_tree.plot([node.x, align_x], [node.y, node.y], color='gray', linestyle=':', lw=1.0, zorder=2)
        dot_color = gen_colors.get(g, 'lightgray') if g != 'Unknown' else 'lightgray'
        ax_tree.scatter(align_x, node.y, color=dot_color, s=40, edgecolor='black', zorder=4)
        ax_tree.text(align_x + max_depth*0.05, node.y, label, va='center', ha='left', fontsize=10, zorder=5)

    for g, sy, ey in blocks:
        sy_line = sy - 0.3 if sy == ey else sy
        ey_line = ey + 0.3 if sy == ey else ey
        ax_tree.plot([strip_x, strip_x], [sy_line, ey_line], color='black', lw=2.0, zorder=5)
        ax_tree.text(strip_x + max_depth*0.05, (sy + ey) / 2.0, g, va='center', ha='left', fontsize=11, fontweight='bold', fontstyle='italic', zorder=5)

    ax_tree.set_xlim(0, max_depth * 4.8)


def plot_composite_figure(tree_file: Path, out_file: Path, taxa_meta: dict, aln_file: Path = None, seq_type: str = "AA", pos_df: pd.DataFrame = None, sim_df: pd.DataFrame = None, target_orfs: list = None, title: str = None, ignore_bl: bool = False):
    tree = Phylo.read(tree_file, 'newick')
    # [触发] 忽略枝长转换 Cladogram
    if ignore_bl: flatten_branch_lengths(tree)
    tree.ladderize(reverse=True)
    terminals = tree.get_terminals()
    leaf_y = {node.name: i for i, node in enumerate(terminals)}
    genome_order = [node.name for node in terminals]
    num_taxa = len(terminals)

    panels, width_ratios = 1, [3.8]
    do_synteny, do_msa = pos_df is not None and not pos_df.empty, aln_file is not None and aln_file.exists()

    if do_synteny: panels += 1; width_ratios.append(3.0)
    if do_msa: panels += 1; width_ratios.append(1.5)

    fig = plt.figure(figsize=(sum(width_ratios)*3.5, max(6, num_taxa * 0.4)))
    gs = fig.add_gridspec(1, panels, width_ratios=width_ratios, wspace=0.08)
    current_ax = 0
    
    ax_tree = fig.add_subplot(gs[current_ax]); current_ax += 1
    max_depth = max(tree.distance(tree.root, n) for n in terminals)
    draw_tree_layer(ax_tree, tree, terminals, leaf_y, taxa_meta, max_depth)
    ax_tree.set_ylim(num_taxa - 0.5, -0.5); ax_tree.axis('off')

    if do_synteny:
        ax_syn = fig.add_subplot(gs[current_ax]); current_ax += 1
        _pos_df, _sim_df = pos_df[pos_df['nucl_id'].isin(genome_order)].copy(), sim_df.copy()
        
        if target_orfs:
            _pos_df = _pos_df[_pos_df['protein'].isin(target_orfs)]
            _sim_df = _sim_df[(_sim_df['protein1'].isin(target_orfs)) & (_sim_df['protein2'].isin(target_orfs))]
            if title: ax_syn.set_title(f"Synteny Map ({title})", fontsize=14, pad=20, fontweight='bold')
        else:
            if title: ax_syn.set_title("Full-Genome Synteny Map", fontsize=14, pad=20, fontweight='bold')

        if not _pos_df.empty:
            nucl_max_len = _pos_df['nucl_length'].max()
            for i, nucl_id in enumerate(genome_order):
                length = pos_df[pos_df['nucl_id']==nucl_id]['nucl_length'].max() if not target_orfs else _pos_df[_pos_df['nucl_id']==nucl_id]['nucl_length'].max()
                if pd.notna(length): ax_syn.hlines(i, 0, length, lw=2, color='#2c3e50')
            
            for _, row in _pos_df.iterrows(): ax_syn.add_patch(patches.Rectangle((row['start'], leaf_y[row['nucl_id']] - 0.15), row['end']-row['start'], 0.3, facecolor='#2ecc71' if target_orfs else 'gray', alpha=0.8, zorder=2))
            protein_to_genome = _pos_df.set_index('protein')['nucl_id'].to_dict()
            sim_cmap = LinearSegmentedColormap.from_list('cg', ['#3498db', '#2ecc71', '#f1c40f', '#e67e22'])
            norm = plt.Normalize(0, 100)
            
            for _, row in _sim_df.iterrows():
                g1, g2 = protein_to_genome.get(row['protein1']), protein_to_genome.get(row['protein2'])
                if g1 and g2 and g1 != g2 and abs(leaf_y[g1] - leaf_y[g2]) == 1:
                    p1, p2 = _pos_df[_pos_df['protein']==row['protein1']].iloc[0], _pos_df[_pos_df['protein']==row['protein2']].iloc[0]
                    ax_syn.add_patch(plt.Polygon([[p1['start'], leaf_y[g1]], [p2['start'], leaf_y[g2]], [p2['end'], leaf_y[g2]], [p1['end'], leaf_y[g1]]], facecolor=sim_cmap(norm(row['similarity'])), alpha=0.3, edgecolor='none', zorder=1))
            
            ax_syn.set_xlim(-nucl_max_len*0.02, nucl_max_len*1.02)
            sm = plt.cm.ScalarMappable(cmap=sim_cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax_syn, orientation='horizontal', fraction=0.03, pad=0.02, aspect=40)
            cbar.set_label('Amino Acid Identity (%)', fontsize=10)
            
        ax_syn.set_ylim(num_taxa - 0.5, -0.5); ax_syn.axis('off')

    if do_msa:
        ax_msa = fig.add_subplot(gs[current_ax])
        if title: ax_msa.set_title(f"MSA Heatmap ({seq_type})", fontsize=14, pad=20, fontweight='bold')
        aln_dict = {rec.id.replace(':', '_').replace('|', '_').replace('-', '_'): str(rec.seq).upper() for rec in SeqIO.parse(aln_file, "fasta")}
        aln_len = len(next(iter(aln_dict.values())))
        if seq_type == "NT": chars, colors = "-ACGTN", ["#FFFFFF", "#BC8F8F", "#FF8247", "#FFEC8B", "#B0E2FF", "#E0E0E0"] 
        else: chars, colors = "-ACDEFGHIKLMNPQRSTVWYX", ["#FFFFFF"] + list(plt.get_cmap('tab20').colors) + ["#E0E0E0"]
        char_map = {c: i for i, c in enumerate(chars)}
        matrix = np.zeros((num_taxa, aln_len), dtype=int)
        for leaf_name, y in leaf_y.items():
            if leaf_name in aln_dict: matrix[y, :] = [char_map.get(c, len(chars)-1) for c in aln_dict[leaf_name]]
        ax_msa.imshow(matrix, aspect='auto', cmap=ListedColormap(colors), interpolation='nearest')
        ax_msa.set_ylim(num_taxa - 0.5, -0.5); ax_msa.axis('off')

    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()


def plot_tanglegram(tree_L_file: Path, tree_R_file: Path, out_file: Path, taxa_meta: dict, title_L="Amino Acid Tree", title_R="Nucleotide Tree", ignore_bl: bool = False):
    tree_L, tree_R = Phylo.read(tree_L_file, 'newick'), Phylo.read(tree_R_file, 'newick')
    # [触发] 忽略枝长转换 Cladogram
    if ignore_bl:
        flatten_branch_lengths(tree_L)
        flatten_branch_lengths(tree_R)
        
    tree_L.ladderize(reverse=True); tree_R.ladderize(reverse=True)
    terminals_L, terminals_R = tree_L.get_terminals(), tree_R.get_terminals()
    leaf_y_L = {node.name: i for i, node in enumerate(terminals_L)}
    leaf_y_R = {node.name: i for i, node in enumerate(terminals_R)}
    
    def assign_y(tree, leaf_y):
        def _assign(node):
            if node.is_terminal(): node.y = leaf_y[node.name]
            else: node.y = sum([_assign(c) for c in node.clades]) / len(node.clades)
            return node.y
        _assign(tree.root)
    assign_y(tree_L, leaf_y_L); assign_y(tree_R, leaf_y_R)

    def assign_x(node, curr_x):
        node.x = curr_x + (node.branch_length or 0.0)
        for c in node.clades: assign_x(c, node.x)
    assign_x(tree_L.root, 0.0); max_x_L = max(n.x for n in terminals_L)
    assign_x(tree_R.root, 0.0); max_x_R = max(n.x for n in terminals_R)

    gap = max_x_L * 2.5 
    offset_R = max_x_L + gap + max_x_R

    def flip_x_R(node):
        node.x = offset_R - node.x
        for c in node.clades: flip_x_R(c)
    flip_x_R(tree_R.root)

    fig, ax = plt.subplots(figsize=(18, max(6, len(terminals_L) * 0.4)))

    def get_blocks(terminals, leaf_y):
        sorted_t = sorted(terminals, key=lambda n: n.y)
        blocks, curr_g, start_y = [], None, None
        for n in sorted_t:
            g = taxa_meta.get(n.name.split('.')[0], {}).get('Genus', 'Unknown')
            if pd.isna(g): g = 'Unknown'
            if g != curr_g:
                if curr_g and curr_g != 'Unknown': blocks.append((curr_g, start_y, n.y - 1))
                curr_g, start_y = g, n.y
        if curr_g and curr_g != 'Unknown': blocks.append((curr_g, start_y, sorted_t[-1].y))
        return blocks
        
    blocks_L, blocks_R = get_blocks(terminals_L, leaf_y_L), get_blocks(terminals_R, leaf_y_R)
    unique_genera = list(set([b[0] for b in blocks_L + blocks_R]))
    cmap_bg = plt.get_cmap('Pastel1')
    gen_colors = {g: cmap_bg(i % 9) for i, g in enumerate(unique_genera)}
    
    for g, sy, ey in blocks_L: ax.add_patch(patches.Rectangle((0, sy-0.5), max_x_L*1.1, ey-sy+1, facecolor=gen_colors[g], alpha=0.3, lw=0, zorder=0))
    for g, sy, ey in blocks_R: ax.add_patch(patches.Rectangle((offset_R - max_x_R*1.1, sy-0.5), max_x_R*1.1, ey-sy+1, facecolor=gen_colors[g], alpha=0.3, lw=0, zorder=0))

    def draw_lines(node, is_left_tree=True):
        if not node.is_terminal():
            min_y, max_y = min(c.y for c in node.clades), max(c.y for c in node.clades)
            ax.plot([node.x, node.x], [min_y, max_y], color='black', lw=1.2, zorder=2)
            
            support = node.confidence if hasattr(node, 'confidence') and node.confidence is not None else getattr(node, 'name', None)
            if support and hasattr(tree_L, 'root') and node != tree_L.root and node != tree_R.root:
                try:
                    val = float(str(support).split('/')[0])
                    if val >= 50:
                        x_offs = -max_x_L*0.01 if is_left_tree else max_x_R*0.01
                        ax.text(node.x + x_offs, node.y - 0.1, str(support), color='blue', fontsize=7, ha=('right' if is_left_tree else 'left'), va='bottom', zorder=4)
                except ValueError: pass

            for c in node.clades:
                ax.plot([node.x, c.x], [c.y, c.y], color='black', lw=1.2, zorder=2)
                draw_lines(c, is_left_tree)
                
    draw_lines(tree_L.root, True); draw_lines(tree_R.root, False)

    align_L = max_x_L * 1.15
    align_R = offset_R - max_x_R * 1.15
    center_text_x = (align_L + align_R) / 2.0

    for n_L in terminals_L:
        n_R = next((n for n in terminals_R if n.name == n_L.name), None)
        base_acc = n_L.name.split('.')[0]
        sp = taxa_meta.get(base_acc, {}).get('Species', 'Unknown')
        g = taxa_meta.get(base_acc, {}).get('Genus', 'Unknown')
        label = f"{sp} ({n_L.name})" if sp != 'Unknown' and not pd.isna(sp) else n_L.name
        
        ax.plot([n_L.x, align_L], [n_L.y, n_L.y], ':', color='gray', lw=1.0)
        ax.plot([n_R.x, align_R], [n_R.y, n_R.y], ':', color='gray', lw=1.0) if n_R else None
        
        dot_color = gen_colors.get(g, 'lightgray') if g != 'Unknown' else 'lightgray'

        if n_R:
            ax.plot([align_L, align_R], [n_L.y, n_R.y], color='gray', alpha=0.3, lw=1.2, zorder=1)
            ax.scatter(align_R, n_R.y, color=dot_color, s=25, edgecolor='black', zorder=4)
            ax.text(center_text_x, (n_L.y + n_R.y)/2.0, label, va='center', ha='center', fontsize=9, zorder=5, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1.5))
        else:
            ax.text(align_L + max_x_L*0.05, n_L.y, label, va='center', ha='left', fontsize=9, zorder=5)

        ax.scatter(align_L, n_L.y, color=dot_color, s=25, edgecolor='black', zorder=4)

    ax.text(max_x_L/2, -1.5, title_L, ha='center', fontsize=14, fontweight='bold')
    ax.text(offset_R - max_x_R/2, -1.5, title_R, ha='center', fontsize=14, fontweight='bold')
    ax.set_ylim(len(terminals_L), -2); ax.axis('off')
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()


# ==========================================
# 主程序
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description='ACVirus Advanced Pangenomic & Evolutionary Pipeline\n(Nature Style Typography & CD-HIT Lineage Dedup)',
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='mode', required=True, help="Analysis Mode:")

    def add_common_args(p):
        group_base = p.add_argument_group('Basic Input Options')
        group_base.add_argument('--genus', required=True, help="Target viral Genus (Required)")
        group_base.add_argument('--contigs', required=False, default=None, help="Input FASTA file containing viral contigs (OPTIONAL)")
        group_base.add_argument('--db_taxa', required=True, help="Path to reference taxa.txt")
        group_base.add_argument('--db_fasta', required=True, help="Path to reference all_virus.fasta")
        group_base.add_argument('--outdir', required=True, help="Output directory path")
        
        group_tree = p.add_argument_group('Tree Building Strategies')
        group_tree.add_argument('--tree_strategy', choices=['trim', 'marker', 'full'], default='trim')
        group_tree.add_argument('--run_modeltest', action='store_true', help="Run ModelTest-NG -> RAxML-NG.")
        group_tree.add_argument('--criterion', choices=['BIC', 'AIC', 'AICc'], default='BIC')
        # [NEW] --ignore_branch_length 支持所有模式强制压缩至 Cladogram
        group_tree.add_argument('--ignore_branch_length', action='store_true', help="Draw tree as cladogram (ignoring real branch lengths)")
        
        group_bait = p.add_mutually_exclusive_group()
        group_bait.add_argument('--bait', type=str, help="FASTA file of marker PROTEIN.")
        group_bait.add_argument('--bait_gb', type=str, help="GenBank file. Parses and builds trees for ALL annotated CDS.")
        
        group_search = p.add_argument_group('Homology Search options')
        group_search.add_argument('--search_tool', choices=['diamond', 'blastp'], default='diamond')
        group_search.add_argument('--evalue', type=str, default='1e-3')
        group_search.add_argument('--max_target_seqs', type=int, default=10000)
        
        group_sample = p.add_argument_group('Sampling Control')
        group_sample.add_argument('--target_genus_ratio', type=float, default=1.0)
        group_sample.add_argument('--other_genus_count', type=int, default=5)
        group_sample.add_argument('--outgroup_count', type=int, default=3)
        
        group_perf = p.add_argument_group('Performance')
        group_perf.add_argument('--threads', type=int, default=8)
        group_perf.add_argument('--resume', action='store_true')

    phylo_parser = subparsers.add_parser('phylogeny', help='Macro-evolution: Reconstruct cross-family tree.')
    add_common_args(phylo_parser)

    simple_parser = subparsers.add_parser('simple', help='Genus-level tree: Target genus + Outgroups from other genera.')
    add_common_args(simple_parser)

    lineage_parser = subparsers.add_parser('lineage', help='Micro-evolution: Reconstruct intra-genus tree.')
    add_common_args(lineage_parser)
    lineage_parser.add_argument('--species', type=str, help="Target Species for Lineage Mode")
    lineage_parser.add_argument('--add_outgroup', action='store_true', help="Auto-sample outgroups from other species")
    lineage_parser.add_argument('--specific_outgroup', type=str, help="Specify EXACT Species name to use as outgroup")
    # 🔥 这里就是之前触发百分号转义崩溃的地方，现在改成了 %% Ns，百分百安全！
    lineage_parser.add_argument('--max_n_ratio', type=float, default=0.05, help="Drop user contigs with >5%% Ns")
    # [NEW] CD-HIT 参数，允许微小尺度建树前高度去冗
    lineage_parser.add_argument('--cdhit', action='store_true', help="Enable CD-HIT intelligent adaptive deduplication")
    lineage_parser.add_argument('--cdhit_threshold', type=float, default=None, help="CD-HIT identity threshold (0.85-1.0). Auto-calculated if neglected.")

    args = parser.parse_args()
    if args.tree_strategy == 'marker' and not (args.bait or args.bait_gb): parser.error("--bait or --bait_gb is REQUIRED when --tree_strategy is 'marker'")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== Starting {args.mode.upper()} Analysis for {args.genus} ===")

    df_taxa = pd.read_csv(args.db_taxa, dtype=str)
    df_taxa['Base_Acc'] = df_taxa['Virus GENBANK accession'].astype(str).apply(lambda x: x.split('.')[0].split(';')[0])
    taxa_meta = df_taxa.set_index('Base_Acc')[['Species', 'Genus', 'Family']].to_dict('index')
    
    final_fasta = outdir / f"{args.genus}_{args.mode}_combined.fasta"
    if not (args.resume and final_fasta.exists() and final_fasta.stat().st_size > 0):
        db_accessions = get_sampled_accessions(df_taxa, args.genus, args.mode, getattr(args, 'species', None), getattr(args, 'add_outgroup', False), getattr(args, 'specific_outgroup', None), args.target_genus_ratio, args.other_genus_count, args.outgroup_count)
        
        if args.contigs:
            user_temp_fasta = outdir / "temp_user_sampled.fasta"
            sample_user_fasta(args.contigs, user_temp_fasta, args.mode, getattr(args, 'max_n_ratio', 0.05))
            user_fasta_to_merge = user_temp_fasta
        else:
            user_fasta_to_merge = None
            
        extract_and_merge(db_accessions, args.db_fasta, user_fasta_to_merge, final_fasta, args.threads)

    # --- [NEW] Intercept and perform CD-HIT if lineage deduplication is requested ---
    if args.mode == 'lineage' and getattr(args, 'cdhit', False):
        dedup_fasta = outdir / f"{args.genus}_{args.mode}_combined_dedup.fasta"
        report_tsv = outdir / f"{args.genus}_Lineage_Clustering_Report.tsv"
        
        if not (args.resume and dedup_fasta.exists() and report_tsv.exists()):
            recs = list(SeqIO.parse(final_fasta, "fasta"))
            ref_len = np.mean([len(r.seq) for r in recs]) if recs else 0
            threshold = args.cdhit_threshold if args.cdhit_threshold else calculate_adaptive_threshold(ref_len)
            
            print(f"\n[INFO] Lineage Dedup: Running CD-HIT Clustering (Threshold: {threshold*100:.2f}%)")
            run_cdhit_and_parse(final_fasta, dedup_fasta, report_tsv, threshold, args.threads)
            print(f"       -> Deduplication Report generated: {report_tsv.name}")
            
        # Re-point the final fasta target to the deduplicated one
        final_fasta = dedup_fasta
    # -------------------------------------------------------------------------

    if args.tree_strategy == 'marker':
        baits = parse_gb_to_baits(Path(args.bait_gb), outdir) if args.bait_gb else [{'gene': 'Target_Protein', 'bait_path': Path(args.bait)}]
        print("\n[INFO] Preparing Background Synteny Database...")
        pos_df, sim_df = prep_synteny_data(final_fasta, outdir, args.threads, args.resume)

        for bait_info in baits:
            gene = bait_info['gene']
            print(f"\n>>> Processing Gene: {gene} <<<")
            out_faa, out_fna = outdir / f"{gene}_extracted_AA.fasta", outdir / f"{gene}_extracted_NT.fasta"
            
            success, target_orfs = predict_and_extract_markers(final_fasta, bait_info['bait_path'], out_faa, out_fna, args.threads, args.resume, args.search_tool, args.evalue, args.max_target_seqs)
            if not success: continue
                
            treefile_aa, aln_aa = build_single_tree(out_faa, outdir / f"Tree_{gene}_AA", "Amino Acid", args.threads, args.resume, args.run_modeltest, args.criterion, args.genus, taxa_meta, args.mode, getattr(args, 'species', None), getattr(args, 'specific_outgroup', None))
            treefile_nt, aln_nt = build_single_tree(out_fna, outdir / f"Tree_{gene}_NT", "Nucleotide", args.threads, args.resume, args.run_modeltest, args.criterion, args.genus, taxa_meta, args.mode, getattr(args, 'species', None), getattr(args, 'specific_outgroup', None))
            
            print(f"    -> [PLOT] Generating advanced Right-Aligned visualizations for {gene}...")
            # 引入参数 ignore_bl 解决绘制要求
            plot_tanglegram(treefile_aa, treefile_nt, outdir/f"{gene}_Tanglegram_AA_vs_NT.png", taxa_meta, title_L=f"{gene} (AA)", title_R=f"{gene} (NT)", ignore_bl=args.ignore_branch_length)
            plot_composite_figure(treefile_aa, outdir/f"{gene}_Tree_with_AA_MSA.png", taxa_meta, aln_file=aln_aa, seq_type="AA", title=gene, ignore_bl=args.ignore_branch_length)
            plot_composite_figure(treefile_nt, outdir/f"{gene}_Tree_with_NT_MSA.png", taxa_meta, aln_file=aln_nt, seq_type="NT", title=gene, ignore_bl=args.ignore_branch_length)
            plot_composite_figure(treefile_aa, outdir/f"{gene}_Tree_with_SingleGene_Synteny.png", taxa_meta, pos_df=pos_df, sim_df=sim_df, target_orfs=target_orfs, title=gene, ignore_bl=args.ignore_branch_length)
            plot_composite_figure(treefile_aa, outdir/f"{gene}_Tree_Synteny_MSA_Combo.png", taxa_meta, aln_file=aln_aa, seq_type="AA", pos_df=pos_df, sim_df=sim_df, target_orfs=target_orfs, title=gene, ignore_bl=args.ignore_branch_length)

    else:
        print("\n[INFO] Running Full-Genome/Trimmed Pipeline...")
        alignment = outdir / 'mafft.fasta'
        if not (args.resume and alignment.exists()):
            print(f"    -> [MAFFT] Aligning sequences...")
            subprocess.check_call(f'mafft --auto --thread {args.threads} {final_fasta} > {alignment} 2>/dev/null', shell=True)
            
        if args.tree_strategy == 'trim':
            trim_out = outdir / 'trim.fasta'
            if not (args.resume and trim_out.exists()):
                print(f"    -> [AliFilter] Trimming alignment...")
                subprocess.check_call(f'AliFilter -p {args.threads} -i {alignment} -o {trim_out}', shell=True)
            alignment = trim_out

        treefile, _ = build_single_tree(alignment, outdir / 'Full', "Nucleotide", args.threads, args.resume, args.run_modeltest, args.criterion, args.genus, taxa_meta, args.mode, getattr(args, 'species', None), getattr(args, 'specific_outgroup', None))
            
        print("\n[INFO] Preparing Background Synteny Database...")
        pos_df, sim_df = prep_synteny_data(final_fasta, outdir, args.threads, args.resume)
        
        print(f"    -> [PLOT] Generating advanced Right-Aligned visualizations...")
        plot_composite_figure(treefile, outdir/"Tree_with_NT_MSA.png", taxa_meta, aln_file=alignment, seq_type="NT", title="Full Genome", ignore_bl=args.ignore_branch_length)
        plot_composite_figure(treefile, outdir/"Tree_with_FullGenome_Synteny.png", taxa_meta, pos_df=pos_df, sim_df=sim_df, title="Full Genome", ignore_bl=args.ignore_branch_length)

    print(f"\n[SUCCESS] Analysis complete! Check {outdir.resolve()} for Nature/Cell style plots.")

if __name__ == '__main__':
    main()
