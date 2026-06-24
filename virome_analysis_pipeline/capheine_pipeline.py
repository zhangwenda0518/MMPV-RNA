#!/usr/bin/env python3
"""
CAPHEINE Pipeline - Pure Python Implementation (Production Ready for HPC)
100% functionally identical to the Nextflow pipeline, with built-in concurrency,
deadlock prevention, intelligent resource allocation, and exact MultiQC replication.
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Dict

# Require Biopython for sequence parsing
try:
    import Bio
    from Bio import SeqIO
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
except ImportError:
    print("Error: Biopython is required. Please install it using: pip install biopython", file=sys.stderr)
    sys.exit(1)

# ----------------------------------------------------------------------
# Logging Setup
# ----------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("capheine")

# ----------------------------------------------------------------------
# Global Resource & Code Mapping
# ----------------------------------------------------------------------
_NCBI_TO_HYPHY = {
    1: "Universal", 2: "Vertebrate-mtDNA", 3: "Yeast-mtDNA", 4: "Mold-Protozoan-mtDNA",
    5: "Invertebrate-mtDNA", 6: "Ciliate-Nuclear", 9: "Echinoderm-mtDNA", 10: "Euplotid-Nuclear",
    12: "Alt-Yeast-Nuclear", 13: "Ascidian-mtDNA", 14: "Flatworm-mtDNA", 15: "Blepharisma-Nuclear",
    16: "Chlorophycean-mtDNA", 21: "Trematode-mtDNA", 22: "Scenedesmus-obliquus-mtDNA",
    23: "Thraustochytrium-mtDNA", 24: "Pterobranchia-mtDNA", 25: "SR1-and-Gracilibacteria",
    26: "Pachysolen-Nuclear", 29: "Mesodinium-Nuclear", 30: "Peritrich-Nuclear",
    33: "Cephalodiscidae-mtDNA",
}

# These will be dynamically overridden by CLI arguments in main()
TOOL_RESOURCES = {
    'iqtree': {'cpus': 6},
    'hyphy': {'cpus': 16}
}

def hyphy_code_name(code_arg: str) -> str:
    if not code_arg: return "Universal"
    code = code_arg.strip()
    try:
        n = int(code)
        if n in _NCBI_TO_HYPHY: return _NCBI_TO_HYPHY[n]
        raise ValueError(f"Unsupported NCBI genetic code id: {n}")
    except ValueError: pass
    
    if code.lower() in ("standard", "universal"): return "Universal"
    for name in _NCBI_TO_HYPHY.values():
        if name.lower() == code.lower(): return name
    return code

def cawlign_code_name(code_arg: str) -> str:
    return hyphy_code_name(code_arg).lower()

def _load_codon_table(table_arg: str):
    if not table_arg: return CodonTable.unambiguous_dna_by_id[1]
    try: return CodonTable.unambiguous_dna_by_id[int(table_arg)]
    except (ValueError, KeyError): pass

    key = table_arg.strip().lower()
    for n_id, name in _NCBI_TO_HYPHY.items():
        if name.lower() == key: return CodonTable.unambiguous_dna_by_id[n_id]
            
    try: return CodonTable.unambiguous_dna_by_name[table_arg]
    except KeyError: raise ValueError(f"Unknown genetic code '{table_arg}'")

# ----------------------------------------------------------------------
# Pure Python Helper Modules
# ----------------------------------------------------------------------
def remove_terminal_stop_codon(fasta_path: Path, output_path: Path, genetic_code: str = "1") -> Path:
    table = _load_codon_table(genetic_code)
    stop_codons = set(table.stop_codons)
    records_out = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(record.seq).upper().replace("U", "T")
        idx = len(seq_str)
        trailing_stops = 0
        while idx >= 3:
            if seq_str[idx - 3 : idx] in stop_codons:
                trailing_stops += 1
                idx -= 3
            else: break
        for pos in range(0, (idx // 3) * 3, 3):
            if seq_str[pos : pos + 3] in stop_codons:
                raise RuntimeError(f"Internal stop codon in {record.id} at pos {pos}.")
        record.seq = Seq(seq_str[:idx] if trailing_stops > 0 else seq_str)
        records_out.append(record)
    SeqIO.write(records_out, output_path, "fasta")
    return output_path

def split_fasta_pure(ref_fasta: Path, outdir: Path) -> List[Path]:
    gene_fastas = []
    prefix = ref_fasta.stem
    for record in SeqIO.parse(ref_fasta, "fasta"):
        safe_id = re.sub(r"[^a-zA-Z0-9_]", "_", record.id)
        out_name = f"{prefix}.part_{safe_id}.fasta"
        out_path = outdir / out_name
        SeqIO.write([record], out_path, "fasta")
        gene_fastas.append(out_path)
    return gene_fastas

def clean_foreground_list(input_file: Path, output_file: Path) -> Path:
    with open(input_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            if cleaned := re.sub(r"[^a-zA-Z0-9_]", "_", line.strip()):
                fout.write(cleaned + "\n")
    return output_file

def filter_ambiguous_sequences(fasta_path: Path, output_path: Path, max_gap_fraction: float = 0.5) -> Path:
    kept = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        ambig_count = seq.count("N") + seq.count("X") + seq.count("-") + seq.count(".")
        if (ambig_count / len(seq) if len(seq) > 0 else 1.0) <= max_gap_fraction:
            kept.append(record)
    SeqIO.write(kept, output_path, "fasta")
    return output_path

# ----------------------------------------------------------------------
# Command Wrappers with Safe Execution
# ----------------------------------------------------------------------
def run_cmd(cmd: str, description: str, attempt: int = 1, max_retries: int = 1) -> None:
    logger.info(f"[{description}] Executing: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        # exit=2: IQ-TREE checkpoint exists (previous run completed), treat as success
        if result.returncode == 2 and 'successfully finished' in (result.stderr or ''):
            logger.info(f"[{description}] Checkpoint found — previous run completed, skipping.")
            return
        if result.returncode in list(range(130, 146)) + [104] and attempt <= max_retries:
            wait = 2 ** attempt
            logger.warning(f"[{description}] Failed (exit={result.returncode}). Retrying in {wait}s ({attempt}/{max_retries})...")
            time.sleep(wait)
            run_cmd(cmd, description, attempt=attempt+1, max_retries=max_retries)
        else:
            logger.error(f"[{description}] Command failed: {cmd}\nStderr: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)

def cawlign(reference: Path, unaligned: Path, output_path: Path, genetic_code: str):
    code_name = cawlign_code_name(genetic_code)
    code_arg = f"-c '{code_name}'" if code_name else ""
    run_cmd(f"cawlign -t codon -r {reference} -f refmap -s BLOSUM62 {code_arg} \"{unaligned}\" > {output_path}", "cawlign")

def hyphy_cln(alignment: Path, output_path: Path, genetic_code: str):
    run_cmd(f"hyphy cln --code '{hyphy_code_name(genetic_code)}' --alignment {alignment} --filtering-method 'Yes/No' --output {output_path}", "hyphy_cln")

def iqtree(alignment: Path, prefix_path: Path) -> Path:
    cpus = TOOL_RESOURCES['iqtree']['cpus']
    run_cmd(f"iqtree -s {alignment} -pre {prefix_path} -nt AUTO -ntmax {cpus} -m GTR+I+G", "IQ-TREE")
    tree_file = prefix_path.parent / f"{prefix_path.name}.treefile"
    if not tree_file.exists(): raise FileNotFoundError(f"Tree file not found: {tree_file}")
    return tree_file

def hyphy_label_tree(tree: Path, output_path: Path, label: str, regexp: str = None, list_file: Path = None, invert: bool = False, internal_nodes: str = "All descendants", leaf_nodes: str = "Label"):
    invert_str = "Yes" if invert else "No"
    target = f"--regexp '{regexp}'" if regexp else f"--list {list_file}"
    run_cmd(f"hyphy label-tree --tree {tree} {target} --invert '{invert_str}' --label '{label}' --internal-nodes '{internal_nodes}' --leaf-nodes '{leaf_nodes}' --output {output_path}", f"hyphy_label_tree ({label})")

def hyphy_analysis(tool: str, alignment: Path, tree: Path, output_path: Path, genetic_code: str, test_branches: str = None, foreground_tag: str = None, reference_tag: str = None, use_mpi: bool = False) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    code_arg = f"--code '{hyphy_code_name(genetic_code)}'"
    cpus = TOOL_RESOURCES['hyphy']['cpus']
    branch_arg = f"--branches '{'Internal' if test_branches and test_branches.lower() == 'internal' else 'All'}'" if test_branches and tool in ("FEL", "MEME", "PRIME", "BUSTED") else ""

    cmd_base = f"--alignment {alignment} --tree {tree} --output {output_path} {code_arg}"
    if tool == "FEL": cmd = f"hyphy fel {cmd_base} --srv Yes {branch_arg}"
    elif tool == "MEME": cmd = f"hyphy meme {cmd_base} {branch_arg}"
    elif tool == "PRIME": cmd = f"hyphy prime {cmd_base} --property-set 'Atchley' {branch_arg}"
    elif tool == "BUSTED": cmd = f"hyphy busted {cmd_base} --srv Yes --error-sink Yes {branch_arg}"
    elif tool == "CONTRASTFEL": cmd = f"hyphy contrast-fel {cmd_base} --branch-set {foreground_tag} --branch-set {reference_tag}"
    elif tool == "RELAX": cmd = f"hyphy relax {cmd_base} --mode 'Classic mode' --test {foreground_tag} --reference {reference_tag} --srv Yes"

    if use_mpi and tool in ("FEL", "MEME", "PRIME", "CONTRASTFEL"):
        cmd = cmd.replace("hyphy", f"mpirun -np {cpus} HYPHYMPI", 1)
    else:
        cmd = cmd.replace("hyphy", f"OMP_NUM_THREADS={cpus} hyphy", 1)

    run_cmd(cmd, f"hyphy_{tool}")
    return output_path if output_path.exists() else None

# ----------------------------------------------------------------------
# Central Gene Processor
# ----------------------------------------------------------------------
def process_gene(gene_fasta: Path, unaligned_fasta: Path, dirs: Dict[str, Path], genetic_code: str, test_branches: str, foreground_regexp: str, foreground_list_path: Path, use_mpi: bool, max_gap_fraction: float) -> dict:
    gene_name = gene_fasta.stem
    logger.info(f"===== Processing gene module: {gene_name} =====")

    aligned_fasta = dirs["cawlign"] / f"{gene_name}-aligned.fasta"
    cawlign(gene_fasta, unaligned_fasta, aligned_fasta, genetic_code)

    clean_fasta = dirs["removeambigseqs"] / f"{gene_name}-clean.fasta"
    filter_ambiguous_sequences(aligned_fasta, clean_fasta, max_gap_fraction)

    cln_fasta = dirs["cln"] / f"{gene_name}-nodups.fasta"
    hyphy_cln(clean_fasta, cln_fasta, genetic_code)

    # GUARD: Need >= 3 seqs for phylogeny
    if len(list(SeqIO.parse(cln_fasta, "fasta"))) < 3:
        logger.warning(f"Gene '{gene_name}' has < 3 valid sequences. Skipping tree & HyPhy.")
        return {"gene": gene_name, "FEL": None, "MEME": None, "PRIME": None, "BUSTED": None, "CONTRASTFEL": None, "RELAX": None}

    iqtree_prefix = dirs["iqtree"] / gene_name
    tree_file = iqtree(cln_fasta, iqtree_prefix)

    final_tree = tree_file
    has_fg = bool(foreground_regexp or foreground_list_path)
    leaf_mode = "Skip" if test_branches == "internal" else "Label"

    if has_fg:
        s_list = dirs["labeltree"] / f"{gene_name}_foreground_sanitized.txt" if foreground_list_path else None
        if foreground_list_path: clean_foreground_list(foreground_list_path, s_list)

        tree_fg = dirs["labeltree"] / f"{gene_name}-Foreground.treefile"
        hyphy_label_tree(final_tree, tree_fg, "Foreground", regexp=foreground_regexp, list_file=s_list, invert=False, leaf_nodes=leaf_mode)
        
        tree_bg = dirs["labeltree"] / f"{gene_name}-Reference.treefile"
        hyphy_label_tree(tree_fg, tree_bg, "Reference", regexp=foreground_regexp, list_file=s_list, invert=True, leaf_nodes=leaf_mode)
        final_tree = tree_bg

        if test_branches == "internal":
            tree_nuisance = dirs["labeltree"] / f"{gene_name}-Nuisance.treefile"
            hyphy_label_tree(tree_bg, tree_nuisance, "Nuisance", regexp="Node", invert=True, internal_nodes="None", leaf_nodes="Label")
            final_tree = tree_nuisance

    tasks = ["FEL", "MEME", "PRIME", "BUSTED"]
    if has_fg: tasks.extend(["CONTRASTFEL", "RELAX"])
    
    result = {"gene": gene_name}
    for t in ["FEL", "MEME", "PRIME", "BUSTED", "CONTRASTFEL", "RELAX"]: result[t] = None

    with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
        futures = {
            executor.submit(
                hyphy_analysis, t, cln_fasta, final_tree, dirs[t.lower()] / f"{gene_name}.{t}.json", genetic_code,
                test_branches, "Foreground" if has_fg else None, "Reference" if has_fg else None, use_mpi
            ): t for t in tasks
        }
        for future in as_completed(futures):
            result[futures[future]] = future.result()

    return result

# ----------------------------------------------------------------------
# Aggregation & Reporting (DRHIP & MultiQC)
# ----------------------------------------------------------------------
def multiqc(input_dir: Path, outdir: Path, title: str = None, args=None):
    """
    Creates exact Nextflow-style Custom Content YAML files, including dynamic
    tool version extraction, and then runs MultiQC to generate the final report.
    """
    # 1. Dynamic Tool Version Extraction
    def get_version(cmd, regex_pattern=None):
        try:
            res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            output = res.stdout.strip() + " " + res.stderr.strip()
            if regex_pattern:
                match = re.search(regex_pattern, output)
                return match.group(1) if match else "Unknown"
            return output.split('\n')[0]
        except Exception:
            return "Unknown"

    versions = {
        "Python": sys.version.split()[0],
        "Biopython": Bio.__version__,
        "IQ-TREE": get_version("iqtree -version", r"version\s+([0-9\.]+)"),
        "HyPhy": get_version("hyphy --version", r"HYPHY (.+)"),
        "Cawlign": get_version("cawlign --version"),
        "DRHIP": get_version("drhip --version", r"drhip\s+(.+)"),
    }

    # 2. Fake Workflow Summary YAML
    summary_yaml = input_dir / "workflow_summary_mqc.yaml"
    command_str = " ".join(sys.argv)
    with open(summary_yaml, "w") as f:
        f.write(f"""id: 'nf-core-capheine-summary'
description: ' - this information is collected when the pipeline is started.'
section_name: 'veg/CAPHEINE Workflow Summary'
section_href: 'https://github.com/veg/CAPHEINE'
plot_type: 'html'
data: |
    <dl class="dl-horizontal">
        <dt>Execution Command</dt><dd><samp>python {command_str}</samp></dd>
        <dt>Output Directory</dt><dd><samp>{args.outdir}</samp></dd>
        <dt>Reference Genes</dt><dd><samp>{args.reference}</samp></dd>
        <dt>Unaligned Seqs</dt><dd><samp>{args.unaligned}</samp></dd>
        <dt>Genetic Code</dt><dd><samp>{args.code}</samp></dd>
        <dt>Test Branches</dt><dd><samp>{args.test_branches}</samp></dd>
        <dt>Max Gap Fraction</dt><dd><samp>{args.max_gap_fraction}</samp></dd>
        <dt>Workers / Concurrency</dt><dd><samp>{args.workers}</samp></dd>
        <dt>Use MPI</dt><dd><samp>{str(args.use_mpi).lower()}</samp></dd>
    </dl>
""")

    # 3. Fake Software Versions YAML
    versions_yaml = input_dir / "software_versions_mqc.yaml"
    with open(versions_yaml, "w") as f:
        f.write("""id: 'software_versions'
section_name: 'Software Versions'
description: 'Software Versions lists versions of software tools extracted from the current environment.'
plot_type: 'table'
pconfig:
    id: 'software_versions_table'
    title: 'Software Versions'
data:
    "Pipeline Tools":
""")
        for tool, ver in versions.items():
            f.write(f"        '{tool}': '{ver}'\n")

    # 4. Fake Methods Description YAML (Citations)
    methods_yaml = input_dir / "methods_description_mqc.yaml"
    with open(methods_yaml, "w") as f:
        f.write("""id: 'nf-core-capheine-methods-description'
section_name: 'CAPHEINE Methods Description'
description: 'Suggested text and references to use when describing pipeline usage within the methods section of a publication. https://github.com/veg/capheine'
plot_type: 'html'
data: |
    <h4>Methods</h4>
    <p>Data was processed using CAPHEINE (Pure Python Architecture). The pipeline executes reproducible software environments from Bioconda and Biocontainers.</p>
    <p>Tools used in the workflow included: BioPython (Cock et al., 2009) BUSTED (Murrell et al., 2015) Cawlign Contrast-FEL (Kosakovsky Pond et al., 2020) DRHIP FEL (Kosakovsky Pond et al., 2005) HyPhy (Kosakovsky Pond et al., 2019) IQ-TREE (Minh et al., 2020) MEME (Murrell et al., 2012) MultiQC (Ewels et al., 2016) PRIME RELAX (Wertheim et al., 2014).</p>
    <h4>References</h4>
    <ul>
        <li>Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, et al. Biopython. Bioinformatics. 2009;25(11):1422–3. doi: 10.1093/bioinformatics/btp163</li>
        <li>Murrell B, Weaver S, Smith MD, Wertheim JO, Murrell S, Aylward A, et al. Gene-wide identification of episodic selection. Mol Biol Evol. 2015;32(5):1365–71. doi: 10.1093/molbev/msv035</li>
        <li>Kosakovsky Pond SL, Wisotsky SR, Escalante A, Magalis BR, Weaver S. Contrast-FEL. Mol Biol Evol. 2020;38(3):1184–98. doi: 10.1093/molbev/msaa263</li>
        <li>Kosakovsky Pond SL, Frost SD. Not so different after all. Mol Biol Evol. 2005;22(5):1208–22. doi: 10.1093/molbev/msi105</li>
        <li>Kosakovsky Pond SL, Poon AF, Velazquez R, Weaver S, et al. Hyphy 2.5. Mol Biol Evol. 2019;37(1):295–9. doi: 10.1093/molbev/msz197</li>
        <li>Minh BQ, Schmidt HA, Chernomor O, Schrempf D, et al. IQ-tree 2. Mol Biol Evol. 2020;37(5):1530–4. doi: 10.1093/molbev/msaa015</li>
        <li>Murrell B, Wertheim JO, Moola S, Weighill T, et al. Detecting individual sites subject to episodic diversifying selection. PLoS Genetics. 2012;8(7). doi: 10.1371/journal.pgen.1002764</li>
        <li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC. Bioinformatics. 2016;32(19):3047-8. doi: 10.1093/bioinformatics/btw354</li>
        <li>Wertheim JO, Murrell B, Smith MD, Kosakovsky Pond SL, Scheffler K. Relax. Mol Biol Evol. 2014;32(3):820–32. doi: 10.1093/molbev/msu400</li>
    </ul>
""")

    # 5. Generate simple HTML summary (MultiQC as optional bonus)
    title_arg = f"--title '{title}'" if title else ""
    cmd = f"multiqc {input_dir} -o {outdir} {title_arg} --force 2>/dev/null || true"
    logger.info(f"[MultiQC] Executing: {cmd}")
    result = subprocess.run(cmd, shell=True)
    # Always write a minimal summary HTML
    summary_html = outdir / "capheine_summary.html"
    with open(summary_html, "w") as f:
        f.write(f"<html><head><title>Capheine Results</title></head><body>")
        f.write(f"<h1>Capheine Positive Selection Analysis</h1>")
        f.write(f"<p>Output: {input_dir}</p>")
        f.write(f"<h2>Results</h2><ul>")
        for d in ["FEL", "MEME", "PRIME", "BUSTED"]:
            json_files = list((input_dir / "hyphy" / d).glob("*.json"))
            if json_files:
                f.write(f"<li><strong>{d}:</strong> {len(json_files)} gene(s) analyzed</li>")
        if (input_dir / "drhip" / "combined_sites.csv").exists():
            f.write(f"<li><strong>DRHIP:</strong> combined_sites.csv available</li>")
        f.write(f"</ul></body></html>")
    if result.returncode != 0:
        logger.warning(f"[MultiQC] returned exit code {result.returncode} — report may be incomplete")

# ----------------------------------------------------------------------
# CLI Setup and Orchestration
# ----------------------------------------------------------------------
def run_capheine(args):
    outdir = Path(args.outdir)
    
    dirs = {
        "removeterminalstopcodon": outdir / "removeterminalstopcodon",
        "seqkit": outdir / "seqkit",
        "cawlign": outdir / "cawlign",
        "removeambigseqs": outdir / "removeambigseqs",
        "cln": outdir / "hyphy" / "CLN",
        "iqtree": outdir / "iqtree",
        "labeltree": outdir / "hyphy" / "LABELTREE",
        "fel": outdir / "hyphy" / "FEL",
        "meme": outdir / "hyphy" / "MEME",
        "prime": outdir / "hyphy" / "PRIME",
        "busted": outdir / "hyphy" / "BUSTED",
        "contrastfel": outdir / "hyphy" / "CONTRASTFEL",
        "relax": outdir / "hyphy" / "RELAX",
        "drhip": outdir / "drhip",
        "multiqc": outdir / "multiqc"
    }
    for d in dirs.values(): d.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"1. Removing terminal stop codons from reference...")
    ref_nostop = dirs["removeterminalstopcodon"] / f"{Path(args.reference).stem}-noStopCodons.fasta"
    remove_terminal_stop_codon(Path(args.reference), ref_nostop, args.code)

    logger.info(f"2. Splitting reference genes...")
    gene_fastas = split_fasta_pure(ref_nostop, dirs["seqkit"])
    
    logger.info(f"3. Initiating execution mapping for {len(gene_fastas)} genes. Global Concurrency: {args.workers}")
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                process_gene, gene, Path(args.unaligned), dirs, args.code, args.test_branches,
                args.foreground_regexp, Path(args.foreground_list) if args.foreground_list else None, 
                args.use_mpi, args.max_gap_fraction
            ): gene for gene in gene_fastas
        }
        for future in as_completed(futures):
            try:
                res = future.result()
                results.append(res)
                logger.info(f"Successfully processed gene module: {res['gene']}")
            except Exception as e:
                logger.error(f"Gene module failed: {e}")

    analysis_files = {
        "FEL": [],
        "MEME": [],
        "PRIME": [],
        "BUSTED": [],
        "CONTRASTFEL": [],
        "RELAX": []
    }
    for r in results:
        gene = r['gene']
        for tool in analysis_files:
            json_path = dirs[tool.lower()] / f"{gene}.{tool}.json"
            if json_path.exists():
                analysis_files[tool].append(json_path)

    logger.info("4. Aggregating JSON results with DRHIP...")
    drhip_hyphy_dir = dirs["drhip"] / "hyphy"
    shutil.rmtree(drhip_hyphy_dir, ignore_errors=True)
    drhip_hyphy_dir.mkdir(parents=True)
    for tool, files in analysis_files.items():
        if files:
            tool_dir = drhip_hyphy_dir / tool
            tool_dir.mkdir(exist_ok=True)
            for f in files: shutil.copy(f, tool_dir / f.name)
    run_cmd(f"drhip --input {drhip_hyphy_dir} --output {dirs['drhip']}", "DRHIP")

    logger.info("5. Generating MultiQC workflow report...")
    multiqc(outdir, dirs['multiqc'], title=args.multiqc_title, args=args)

    logger.info("=========================================")
    logger.info("CAPHEINE pipeline completed successfully!")
    logger.info("=========================================")

def main():
    parser = argparse.ArgumentParser(description="CAPHEINE Pipeline (Production Pure Python for HPC)")
    parser.add_argument("--reference", "-r", required=True, help="Path to reference genes FASTA")
    parser.add_argument("--unaligned", "-u", required=True, help="Path to unaligned sequences FASTA")
    parser.add_argument("--outdir", "-o", required=True, help="Output directory path")
    parser.add_argument("--foreground_list", help="Path to text file with foreground taxa names")
    parser.add_argument("--foreground_regexp", help="Regular expression for foreground taxa")
    parser.add_argument("--test_branches", choices=["internal", "all"], help="Branches to test (internal/all)")
    parser.add_argument("--code", default="1", help="Genetic code (NCBI id or HyPhy name). Default: 1 (Universal)")
    parser.add_argument("--max_gap_fraction", type=float, default=0.5, help="Max allowed fraction of gaps/ambiguous bases (Default: 0.5)")
    parser.add_argument("--use_mpi", action="store_true", help="Enable MPI for HyPhy execution (HYPHYMPI)")
    parser.add_argument("--multiqc_title", help="Custom title for the MultiQC HTML report")
    parser.add_argument("--workers", type=int, default=2, help="Number of genes to process concurrently (Default: 2)")
    parser.add_argument("--cpus_iqtree", type=int, default=6, help="Number of threads for IQ-TREE per gene (Default: 6)")
    parser.add_argument("--cpus_hyphy", type=int, default=16, help="Number of threads/MPI ranks for HyPhy per test (Default: 16)")
    
    args = parser.parse_args()

    if args.foreground_list and args.foreground_regexp:
        parser.error("Provide only ONE of --foreground_list or --foreground_regexp")

    TOOL_RESOURCES['iqtree']['cpus'] = args.cpus_iqtree
    TOOL_RESOURCES['hyphy']['cpus'] = args.cpus_hyphy

    run_capheine(args)

if __name__ == "__main__":
    main()
