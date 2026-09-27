#!/usr/bin/env python3
"""known_virus_submit.py — cawlign codon-aware 参考引导 GenBank 提交

cawlign 比对 assembly → reference → refmap 格式 (查询插入相对参考被移除)
→ 参考CDS坐标直接映射到对齐列 → 提取assembly CDS → .tbl → .gb

优势: 零假ORF, 零blastp, curator注释直接继承
"""

import argparse, os, re, sys, csv, logging, subprocess, shutil, tempfile
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

def _setup_logger(out_dir):
    logger = logging.getLogger("KnownSubmit"); logger.setLevel(logging.DEBUG); logger.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)
    ch = logging.StreamHandler(); ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)
    return logger

def _run(cmd, log, step, check=False):
    log.info("[%s] %s", step, cmd[:160])
    try:
        subprocess.run(cmd, shell=True, check=check, executable='/bin/bash')
        log.info("[%s] OK", step); return True
    except subprocess.CalledProcessError as e:
        log.error("[%s] FAIL (%d)", step, e.returncode); return False

# ═══════ Collect ═══════

def collect_viruses(d, log):
    d = Path(d)
    fs = sorted(d.glob("*.full.fasta")) + sorted(d.glob("*.fasta"))
    if fs:
        m = re.search(r'([A-Z]{1,4}_?\d+\.\d+)', d.name)
        acc = m.group(1) if m else d.name
        return {acc: {"name": d.name, "dir": d, "fastas": fs, "n": len(fs)}}
    vs = {}
    for sd in sorted(d.iterdir()):
        if not sd.is_dir() or sd.name.startswith('.'): continue
        fs = sorted(sd.glob("*.full.fasta")) or sorted(sd.glob("*.fasta"))
        if not fs: continue
        m = re.search(r'([A-Z]{1,4}_?\d+\.\d+)', sd.name)
        acc = m.group(1) if m else sd.name
        vs[acc] = {"name": sd.name, "dir": sd, "fastas": fs, "n": len(fs)}
    log.info("collect %d viruses (%d seqs)", len(vs), sum(v["n"] for v in vs.values()))
    return vs

# ═══════ Ref DB ═══════

def load_ref_db(ref_db_dir, log):
    r = Path(ref_db_dir)
    ri = {}; 
    with open(r / "final.cluster.ref_info.tsv") as f:
        for row in csv.DictReader(f, delimiter='\t'):
            a = row.get("Accession","")
            if a: ri[a] = row
    rg = {}
    for rec in SeqIO.parse(r / "final.cluster.ref.fasta", "fasta-pearson"): rg[rec.id] = rec
    log.info("  ref: %d info, %d genomes", len(ri), len(rg))
    return ri, rg

# ═══════ NCBI curator CDS ═══════

def fetch_ncbi_cds(accessions, cache_dir, log):
    cd = Path(cache_dir); cd.mkdir(parents=True, exist_ok=True)
    products = {}
    for acc in accessions:
        gb = cd / f"{acc}.gb"
        if not gb.exists():
            try:
                from Bio import Entrez; Entrez.email = "zhangwenda@example.com"
                log.info("  fetch %s", acc)
                h = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
                with open(gb, "w") as f: f.write(h.read()); h.close()
            except Exception as e: log.warning("  %s fetch fail: %s", acc, e); continue
        try:
            cds = []
            for rec in SeqIO.parse(gb, "genbank"):
                for feat in rec.features:
                    if feat.type == "CDS":
                        cds.append((
                            int(feat.location.start) + 1,
                            int(feat.location.end),
                            "+" if feat.location.strand == 1 else "-",
                            feat.qualifiers.get("product", ["hypothetical protein"])[0],
                            feat.qualifiers.get("gene", [""])[0]
                        ))
            if cds:
                products[acc] = cds
                log.info("  %s: %d curated CDS", acc, len(cds))
        except Exception as e: log.warning("  %s parse fail: %s", acc, e)
    return products

# ═══════ Taxonomy ═══════

def gen_taxonomy(viruses, ref_info, out_dir, log):
    td = Path(out_dir) / "1_taxonomy"; td.mkdir(parents=True, exist_ok=True)
    tr, mr = [], []
    for acc, vi in sorted(viruses.items()):
        ri = ref_info.get(acc, {})
        sp = ri.get("Species_NCBI", acc); fm = ri.get("VMR_Family","")
        gn = ri.get("VMR_Genus",""); tx = ri.get("Taxid","")
        tp = ri.get("Topology","linear"); ml = ri.get("Molecule_type","ssRNA(+)")
        for fa in vi["fastas"]:
            c = fa.stem
            tr.append([c, f"Viruses;{fm};{gn};{sp}", "species", fm, gn, sp, tx, tp, ml])
            try:
                recs = list(SeqIO.parse(fa, "fasta-pearson"))
                if not recs: continue
                mr.append([c, tx, sp, fm, gn, tp, ml, str(fa.stat().st_size),
                           str(len(recs[0].seq))])
            except Exception:
                continue
    for fn, h, d in [
        ("taxonomy.tsv", ["contig","taxonomy","rank","family","genus","species","taxid","topology","molecule_type"], tr),
        ("miuvig_taxonomy.tsv", ["sample_id","taxid","species","family","genus","topology","molecule_type","file_size","seq_length"], mr)
    ]:
        with open(td / fn, "w", newline='') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(h); w.writerows(d)
    log.info("  taxonomy: %d records", len(tr))
    return td

# ═══════ Features: cawlign ═══════

CAWLIGN = shutil.which("cawlign") or os.path.expanduser("~/.pixi/bin/cawlign")

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def _translate(seq_str):
    s = str(seq_str).upper(); aa = []
    for i in range(0, len(s)-2, 3):
        codon = s[i:i+3]
        if len(codon) < 3: break
        aa.append(CODON_TABLE.get(codon, 'X'))
    return Seq(''.join(aa))

def gen_features(viruses, ref_genomes, ncbi_products, ref_info, out_dir, log):
    """cawlign refmap → extract CDS (并行 cawlign + 串行 CDS提取)"""
    fd = Path(out_dir) / "2_features"; fd.mkdir(parents=True, exist_ok=True)
    all_tbl, all_fna, all_faa = [], [], []
    ng, ns = 0, 0

    for acc, vi in tqdm(sorted(viruses.items()), desc="Features", unit="virus"):
        ref_rec = ref_genomes.get(acc)
        if not ref_rec: ns += vi["n"]; continue
        ref_len = len(ref_rec.seq)

        curator_cds = ncbi_products.get(acc, [])
        if not curator_cds:
            for fa in vi["fastas"]:
                try:
                    recs = list(SeqIO.parse(fa, "fasta-pearson"))
                    if recs:
                        rec = recs[0]
                        all_tbl.extend([f">Feature {rec.id}",
                            f"1\t{len(rec.seq)}\tsource",
                            f"\t\t\tmol_type\tgenomic RNA",
                            f"\t\t\torganism\t{ref_info.get(acc,{}).get('Species_NCBI',acc)}"])
                except Exception: pass
            log.info("  %s: viroid/no-CDS, %d source-only", acc, vi["n"])
            ns += vi["n"]
            continue

        ref_tmp = fd / "tmp" / f"{acc}_ref.fna"
        ref_tmp.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(ref_rec, str(ref_tmp), "fasta")

        # 准备 assembly 列表
        tasks = []
        for fa in vi["fastas"]:
            try:
                recs = list(SeqIO.parse(fa, "fasta-pearson"))
                if not recs: continue
                rec = recs[0]
            except Exception: continue
            asm_tmp = fd / "tmp" / f"{rec.id}.fna"
            if not asm_tmp.exists():
                SeqIO.write(rec, str(asm_tmp), "fasta")
            aln_out = fd / "tmp" / f"{rec.id}_aln.fasta"
            tasks.append((rec, fa, asm_tmp, aln_out))

        if not tasks: continue

        # 并行 cawlign
        def _do_cawlign(rec, fa_path, asm_tmp, aln_out):
            if aln_out.exists() and aln_out.stat().st_size > 500:
                return rec, fa_path, aln_out, True
            cmd = f"{CAWLIGN} -t nucleotide -r {ref_tmp} -f refmap -o {aln_out} {asm_tmp}"
            try:
                subprocess.run(cmd, shell=True, capture_output=True, timeout=120)
                ok = aln_out.exists() and aln_out.stat().st_size > 100
            except Exception:
                ok = False
            return rec, fa_path, aln_out, ok

        with ThreadPoolExecutor(max_workers=8) as ex:
            futures = {ex.submit(_do_cawlign, r, f, a, o): r.id for r, f, a, o in tasks}
            cawlign_results = {}
            for fut in as_completed(futures):
                rec, fa_path, aln_out, ok = fut.result()
                cawlign_results[rec.id] = (rec, fa_path, aln_out, ok)

        # 串行 CDS 提取 (很快, 不需要并行)
        for rec, fa_path, aln_out, ok in [cawlign_results[r.id] for r, f, a, o in tasks]:
            if not ok: ns += 1; continue

            try:
                aln_recs = list(SeqIO.parse(aln_out, "fasta-pearson"))
            except Exception: ns += 1; continue
            if not aln_recs: ns += 1; continue
            aln_seq = str(aln_recs[0].seq).upper()
            aln_len = len(aln_seq)
            if aln_len < ref_len * 0.5: ns += 1; continue

            tbl_block = [f">Feature {rec.id}"]
            gi = 0
            for (cds_start, cds_end, strand, product, gene_name) in curator_cds:
                col_s = cds_start - 1; col_e = cds_end - 1
                if col_s >= aln_len or col_e >= aln_len: continue
                bases = []
                for col in range(col_s, min(col_e + 1, aln_len)):
                    b = aln_seq[col]
                    if b != '-': bases.append(b)
                cds_str = ''.join(bases)
                if len(cds_str) < 30: continue
                rem = len(cds_str) % 3
                if rem: cds_str = cds_str[:len(cds_str) - rem]
                gi += 1; gn = f"gene_{gi:02d}"
                cds_seq = Seq(cds_str)
                pep = _translate(cds_str)
                lt = f"{acc.replace('.','_')}_{gi:02d}"
                asm_start, asm_end = _est_asm_coords(aln_seq, cds_start, cds_end)
                tbl_block += [
                    f"{asm_start}\t{asm_end}\tgene", f"\t\t\tgene\t{gn}",
                    f"{asm_start}\t{asm_end}\tCDS", f"\t\t\tproduct\t{product}",
                    f"\t\t\tlocus_tag\t{lt}"]
                all_faa.append(SeqRecord(pep, id=f"{rec.id}_{gn}", description=f"[{product}] [{acc}]"))
                all_fna.append(SeqRecord(cds_seq, id=f"{rec.id}_{gn}", description=f"[{product}] [{acc}]"))
                ng += 1
            all_tbl.extend(tbl_block)

    # 写 .tbl
    tbl_path = fd / "featuretable.tbl"
    with open(tbl_path, "w") as f:
        all_ids = []
        for v in viruses.values():
            for fa in v["fastas"]:
                try:
                    recs = list(SeqIO.parse(fa, "fasta-pearson"))
                    if recs: all_ids.append(recs[0].id)
                except Exception: pass
        f.write(">Feature " + " ".join(all_ids if all_ids else ["no_seqs"]))
        f.write("\n" + "\n".join(all_tbl) + "\n")
    faa_path = fd / "proteins.faa"; SeqIO.write(all_faa, faa_path, "fasta")
    fna_path = fd / "reoriented_nucleotide_sequences.fna"; SeqIO.write(all_fna, fna_path, "fasta")
    log.info("  .tbl: %d genes (%d skipped)", ng, ns)
    return tbl_path, faa_path, fna_path

def _est_asm_coords(aln_seq_str, ref_start, ref_end):
    """从 cawlign refmap 对齐估算 assembly 坐标 (不跑 blastn)
    
    refmap: 列 = 参考位置, '-' = assembly 上的 gap。
    通过计数 gap 之前的总碱基数估算 assembly 坐标。
    """
    # 计算 start: ref_start 列之前有多少非 gap 碱基
    pre_bases = sum(1 for c in aln_seq_str[:ref_start-1] if c != '-')
    in_bases = sum(1 for c in aln_seq_str[ref_start-1:ref_end] if c != '-')
    return pre_bases + 1, pre_bases + in_bases
    """blastn CDS → assembly 获取原始坐标，超时回退到线性估算"""
    tq = tempfile.mktemp(suffix=".fna"); to = tempfile.mktemp(suffix=".txt")
    seq_len = len(cds_seq_str)
    with open(tq, "w") as f: f.write(f">cds\n{cds_seq_str}\n")
    try:
        subprocess.run(
            f"blastn -task blastn -query {tq} -subject {asm_fasta} -out {to} "
            f"-outfmt '6 sstart send' -word_size 11 -dust no -evalue 1e-20",
            shell=True, capture_output=True, timeout=60)
    except (subprocess.TimeoutExpired, Exception):
        pass
    s, e = 1, seq_len
    if os.path.exists(to) and os.path.getsize(to) > 5:
        with open(to) as f:
            cols = f.readline().strip().split("\t")
            if len(cols) >= 2:
                s, e = int(cols[0]), int(cols[1])
                if s > e: s, e = e, s
    for f in [tq, to]:
        try: os.unlink(f)
        except: pass
    return s, e

# ═══════ .gb ═══════

def run_tbl2gb(tbl_path, tax_dir, feat_dir, log):
    gb_dir = Path(feat_dir) / "genbank"; gb_dir.mkdir(parents=True, exist_ok=True)
    script = None
    for d in [Path(__file__).parent.parent / "virome_analysis_pipeline" / "utils" / "tbl2gb.py",
              Path.home() / "MMPV-RNA" / "virome_analysis_pipeline" / "utils" / "tbl2gb.py"]:
        if d.exists(): script = d; break
    if not script: log.warning("tbl2gb.py not found"); return gb_dir
    full_fna = Path(feat_dir) / "full_assemblies.fna"
    _run(f"python {script} --tbl {tbl_path} --fasta {full_fna} "
         f"--taxonomy {Path(tax_dir)/'taxonomy.tsv'} -o {gb_dir}", log, "tbl2gb")
    gb = list(gb_dir.glob("*.gb"))
    if gb: log.info("  .gb: %d files", len(gb))
    return gb_dir

# ═══════ Metadata ═══════

def run_suvtk_sqn(tax_dir, feat_dir, meta_dir, output_dir, log):
    """suvtk comments + table2asn → .sqn"""
    out = Path(output_dir)
    
    # 生成正确的 template.sbt
    sbt = out / "template.sbt"
    with open(sbt, "w") as f:
        f.write('Submit-block ::= {\n'
          '  contact {\n'
          '    contact {\n'
          '      name name { last "Zhang", first "Wenda", middle "", initials "", suffix "", title "" },\n'
          '      affil std { affil "Ningxia University", div "College of Life Sciences",\n'
          '        city "Yinchuan", sub "Ningxia", country "China",\n'
          '        street "Ningxia University", email "zhangwenda@example.com", postal-code "750021" }\n'
          '    }\n'
          '  },\n'
          '  cit { authors { names std { { name name { last "Zhang", first "Wenda", middle "", initials "", suffix "", title "" } } },\n'
          '    affil std { affil "Ningxia University", div "College of Life Sciences",\n'
          '      city "Yinchuan", sub "Ningxia", country "China", street "Ningxia University", postal-code "750021" } } },\n'
          '  subtype new\n'
          '}\n'
          'Seqdesc ::= pub { pub { gen { cit "unpublished",\n'
          '    authors { names std { { name name { last "Zhang", first "Wenda", middle "", initials "", suffix "", title "" } } } },\n'
          '    title "Plant virome of Lycium chinense in Ningxia, China" } } }\n'
          'Seqdesc ::= user { type str "Submission",\n'
          '  data { { label str "AdditionalComment", data str "ALT EMAIL:zhangwenda@example.com" } } }\n')
    
    # miuvig + assembly
    meta_dir_p = Path(meta_dir)
    if not (meta_dir_p / "miuvig.tsv").exists():
        with open(meta_dir_p / "miuvig.tsv", "w") as f:
            f.write("MIUVIG_parameter\tvalue\n"
                    "viral_enrichment\trRNA_depletion\n"
                    "sequencing_platform\tIllumina_NovaSeq\n"
                    "sequencing_method\tmetatranscriptomic\n"
                    "assembly_software\tMEGAHIT;1.2.9;default parameters\n"
                    "assembly_method\treference_guided\n"
                    "quality_check_software\tCheckV\n"
                    "source_uvig\tmetatranscriptome (not viral targeted)\n")
    if not (meta_dir_p / "assembly.tsv").exists():
        with open(meta_dir_p / "assembly.tsv", "w") as f:
            f.write("Assembly_parameter\tvalue\n"
                    "Sequencing Technology\tIllumina NovaSeq\n"
                    "Assembly Method\treference guided\n")
    
    # miuvig_features.tsv
    tbl = Path(feat_dir) / "featuretable.tbl"
    n_cds, n_hypo = 0, 0
    if tbl.exists():
        with open(tbl) as f:
            for line in f:
                if "product" in line:
                    n_cds += 1
                    if "hypothetical" in line.lower(): n_hypo += 1
    mf = Path(feat_dir) / "miuvig_features.tsv"
    with open(mf, "w") as f:
        f.write(f"MIUVIG_parameter\tvalue\nn_CDS\t{n_cds}\nn_hypothetical\t{n_hypo}\n")
    
    # suvtk comments
    cmt_dir = out / "4_comments"
    miuvig_tax = Path(tax_dir) / "miuvig_taxonomy.tsv"
    cmd = (f"suvtk comments --taxonomy {miuvig_tax} --features {mf} "
           f"--miuvig {meta_dir_p/'miuvig.tsv'} --assembly {meta_dir_p/'assembly.tsv'} -o {cmt_dir}")
    _run(cmd, log, "suvtk comments")
    
    # suvtk table2asn (uses full_assemblies.fna, not CDS)  
    sqn_dir = out / "5_submission"; sqn_dir.mkdir(exist_ok=True)
    fna = Path(feat_dir) / "full_assemblies.fna"
    src = meta_dir_p / "source.clean.src"
    # 去注释的 source.src
    src_all = meta_dir_p / "source.src.template"
    if not src.exists():
        with open(src, "w") as sf:
            sf.write("Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\tLat_Lon\tBioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment\n")
            with open(src_all) as af:
                for line in af:
                    if not line.startswith("#") and "\t" in line:
                        sf.write(line)
    
    cmt_file = cmt_dir / "output.cmt" if (cmt_dir / "output.cmt").exists() else out / "4_comments.cmt"
    if not cmt_file.exists():
        cmt_file = out / "4_comments.cmt"
    cmd = (f"suvtk table2asn -i {fna} -f {tbl} -s {src} -t {sbt} -c {cmt_file} -o {sqn_dir}")
    _run(cmd, log, "suvtk table2asn", check=False)
    
    sqn_files = list(sqn_dir.glob("*.sqn"))
    if sqn_files: log.info("  .sqn: %d files", len(sqn_files))
    return sqn_files


def gen_metadata(viruses, tax_dir, out_dir, meta_file, log):
    md = Path(out_dir) / "3_metadata"; md.mkdir(parents=True, exist_ok=True)
    ml = {}
    if meta_file and os.path.exists(meta_file):
        try:
            import pandas as pd
            df = pd.read_csv(meta_file, sep=None, engine='python')
            for _, row in df.iterrows():
                rn = str(row.get('Run','')).strip()
                if not rn or rn.lower() in ('nan',''): continue
                loc = str(row.get('Location','')); ll = str(row.get('Lat_Lon',''))
                if (not ll or ll.lower() in ('nan','','not_provided')):
                    lat, lon = row.get('lat',''), row.get('lon','')
                    if lat and lon and str(lat).lower() not in ('nan','') and str(lon).lower() not in ('nan',''):
                        ll = f"{lat} N {lon} E" if float(lat)>=0 else f"{-float(lat)} S {abs(float(lon))} W"
                ml[rn] = {
                    'date': str(row.get('CollectionDate','')).split(' ')[0],
                    'geo': loc if loc and 'nan' not in str(loc).lower() else 'Country:Region',
                    'll': ll or 'XX.XX N XXX.XX E',
                    'bp': str(row.get('BioProject','PRJNAXXXXXX')),
                    'bs': str(row.get('BioSample','SAMNXXXXXXXX')),
                    'src': str(row.get('Tissue','plant virome'))}
            log.info("  meta: %d runs", len(ml))
        except Exception as e: log.warning("  meta fail: %s", e)

    tm = {}
    with open(Path(tax_dir)/"taxonomy.tsv") as f:
        h = f.readline().strip().split('\t')
        ci = h.index('contig') if 'contig' in h else 0
        ti = h.index('species') if 'species' in h else 3
        for l in f:
            c = l.strip().split('\t')
            if len(c) > max(ci,ti): tm[c[ci]] = c[ti]

    def _v(sra, k, dflt):
        if sra in ml:
            v = ml[sra].get(k, dflt)
            if v and 'nan' not in str(v).lower() and v != 'not_provided': return v
        return dflt

    al = []
    for acc, vi in sorted(viruses.items()):
        vl = []
        for fa in vi["fastas"]:
            c = fa.stem
            sra = (re.match(r'([SC]RR\d+)', c) or ["UNKNOWN"])[0]
            sra = sra.group(1) if hasattr(sra, 'group') else sra
            org = tm.get(c, acc)
            ln = "\t".join([c, org, f"{acc}_{sra}",
                _v(sra,'date','DD-Mmm-YYYY'), _v(sra,'geo','Country:Region'),
                _v(sra,'ll','XX.XX N XXX.XX E'), _v(sra,'bp','PRJNAXXXXXX'),
                _v(sra,'bs','SAMNXXXXXXXX'), sra, "TRUE",
                _v(sra,'src','plant virome'), ""])
            vl.append(ln); al.append(ln)
        with open(md / f"source_{acc}.src.template", "w") as f:
            f.write(f"# source.src - {acc}\nseqs: {len(vl)}\n")
            f.write("Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\tLat_Lon\t"
                    "Bioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment\n")
            f.write("\n".join(vl) + "\n")
    with open(md / "source.src.template", "w") as f:
        f.write(f"# source.src - all ({len(al)})\n")
        f.write("Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\tLat_Lon\t"
                "Bioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment\n")
        f.write("\n".join(al) + "\n")
    log.info("  source.src: %d records", len(al))
    
    # miuvig.tsv
    with open(md / "miuvig.tsv", "w") as f:
        f.write("MIUVIG_parameter\tvalue\n")
        f.write("viral_enrichment\trRNA_depletion\nsequencing_platform\tIllumina_NovaSeq\n"
                "sequencing_method\tmetatranscriptomic\nassembly_software\tMEGAHIT\n"
                "assembly_method\treference_guided\nquality_check_software\tCheckV\n"
                "source_uvig\tmetatranscriptome (not viral targeted)\n")
    # assembly.tsv
    with open(md / "assembly.tsv", "w") as f:
        f.write("Assembly_parameter\tvalue\n"
                "Sequencing Technology\tIllumina NovaSeq\n"
                "Assembly Method\treference guided\n"
                "Assembly Name\tMMPV-RNA v3.0\n"
                "Assembly Software\tMEGAHIT\n"
                "Coverage\tNOT_PROVIDED\n")
    return md

# ═══════ main ═══════

def main():
    p = argparse.ArgumentParser(description="known_virus_submit v3 — cawlign + NCBI curator")
    p.add_argument("--extraction-dir", required=True)
    p.add_argument("--ref-db", required=True)
    p.add_argument("--metadata")
    p.add_argument("--output", "-o", default="./gb_submission")
    p.add_argument("--threads", "-t", type=int, default=24)
    args = p.parse_args()
    out = Path(args.output).resolve(); log = _setup_logger(str(out))
    log.info("v3 cawlign + NCBI curator CDS")

    viruses = collect_viruses(args.extraction_dir, log)
    if not viruses: sys.exit(1)
    ref_info, ref_genomes = load_ref_db(args.ref_db, log)
    ncbi_cds = fetch_ncbi_cds(list(viruses.keys()), str(out / "ncbi_cache"), log)
    tax_dir = gen_taxonomy(viruses, ref_info, str(out), log)
    tbl, faa, fna = gen_features(viruses, ref_genomes, ncbi_cds, ref_info, str(out), log)

    feat_dir = out / "2_features"
    with open(str(feat_dir / "full_assemblies.fna"), "w") as ff:
        for v in viruses.values():
            for fa in v["fastas"]:
                try:
                    recs = list(SeqIO.parse(fa, "fasta-pearson"))
                    if recs: r = recs[0]
                    else: continue
                    ff.write(f">{r.id}\n{r.seq}\n")
                except Exception: continue
    run_tbl2gb(tbl, str(tax_dir), str(feat_dir), log)
    gen_metadata(viruses, str(tax_dir), str(out), args.metadata, log)
    run_suvtk_sqn(str(tax_dir), str(feat_dir), str(out / "3_metadata"), str(out), log)
    log.info("Done: %s", out)

if __name__ == "__main__":
    main()
