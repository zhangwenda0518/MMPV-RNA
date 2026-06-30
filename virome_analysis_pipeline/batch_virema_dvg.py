#!/usr/bin/env python3
"""
batch_virema_dvg.py — 高通量 DVG 与病毒重组全自动流水线
整合 virema/src 全部能力: ViReMa + Compiler + Visualize + DI-tector
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import colorlog
import pandas as pd
from tqdm import tqdm

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.lines import Line2D
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

SCRIPT_DIR = Path(__file__).resolve().parent

# ==========================================
# 1. 日志与系统工具
# ==========================================
def setup_logging(verbose: bool = False) -> logging.Logger:
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        "%(log_color)s%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
        log_colors={"DEBUG": "cyan", "INFO": "green", "WARNING": "yellow", "ERROR": "red", "CRITICAL": "bold_red"},
    ))
    log = colorlog.getLogger("virema_batch")
    log.handlers.clear()
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    return log

logger = setup_logging()

def safe_name(s: str, max_len: int = 100) -> str:
    s = str(s)
    s = re.sub(r'[^A-Za-z0-9\-.]', '_', s)
    s = re.sub(r'_+', '_', s)
    return s.strip('_.')[:max_len]


def build_runtime_env() -> dict:
    """Ensure Python is in PATH for subprocess (pattern from virema/src)."""
    env = os.environ.copy()
    py_bin = str(Path(sys.executable).resolve().parent)
    current_path = env.get("PATH", "")
    if py_bin not in current_path.split(os.pathsep):
        env["PATH"] = py_bin + os.pathsep + current_path
    return env


def run_cmd(cmd: list, log_path: str = None, check: bool = True):
    """Execute a command list with logging."""
    cmd_str = " ".join(str(c) for c in cmd)
    result = subprocess.run(cmd, capture_output=True, text=True, env=build_runtime_env())

    if log_path:
        lp = Path(log_path)
        lp.parent.mkdir(parents=True, exist_ok=True)
        with open(lp, "a", encoding="utf-8") as f:
            f.write(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] CMD: {cmd_str}\n")
            f.write(f"EXIT_CODE: {result.returncode}\n")
            if result.stdout: f.write(f"--- STDOUT ---\n{result.stdout.strip()}\n")
            if result.stderr: f.write(f"--- STDERR ---\n{result.stderr.strip()}\n")
            f.write("-" * 60 + "\n")

    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd_str, output=result.stdout, stderr=result.stderr)
    return result


# ==========================================
# 2. 弧线图可视化 (from virema/src/visualize.py)
# ==========================================
def _safe_filename(text):
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in text)


def plot_arc_diagram(bed_path: Path, chrom_name: str, output_dir: Path, genome_length=None):
    """Arc diagram: Deletion (red) / Duplication (blue) / Back-Splice (purple)."""
    if not HAS_MPL:
        logger.warning("matplotlib not available, skipping arc diagram")
        return
    if not bed_path.exists():
        return

    data = []
    with open(bed_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("track") or not line: continue
            cols = line.split("\t")
            try:
                start, stop = int(cols[1]), int(cols[2])
                data.append({
                    "Start": start, "Stop": stop, "Type": cols[3],
                    "Count": int(cols[4]) if len(cols) > 4 else 1,
                })
            except (ValueError, IndexError):
                continue
    if not data:
        return

    fig, ax = plt.subplots(figsize=(14, 6))
    max_count = max(d["Count"] for d in data)
    color_map = {"Deletion": "red", "Duplication": "blue", "Back-Splice": "purple", "Insertion": "green"}

    for row in data:
        start, stop = row["Start"], row["Stop"]
        center = (start + stop) / 2
        width = abs(stop - start)
        alpha = min(1.0, 0.1 + 0.9 * (np.log1p(row["Count"]) / np.log1p(max_count))) if max_count > 1 else 1.0
        color = next((c for k, c in color_map.items() if k in row["Type"]), "gray")

        arc = patches.Arc((center, 0), width, width, theta1=0, theta2=180,
                          edgecolor=color, alpha=alpha, linewidth=1.5)
        ax.add_patch(arc)

    limit = genome_length or max(d["Stop"] for d in data)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, max(1, max(abs(d["Stop"] - d["Start"]) for d in data) / 1.5))
    ax.set_xlabel("Genome Position (nt)"); ax.set_ylabel("Jump Distance")
    ax.set_title(f"Recombination Arc Diagram\n{chrom_name}")
    legend_elements = [Line2D([0], [0], color=c, lw=2, label=t) for t, c in color_map.items()]
    ax.legend(handles=legend_elements, loc="upper right", title="Event Type")

    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"arc_{_safe_filename(chrom_name)}.png", dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / f"arc_{_safe_filename(chrom_name)}.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_top_events(bed_path: Path, chrom_name: str, output_dir: Path, top_n=15):
    """Bar chart of top N recombination events."""
    if not HAS_MPL or not bed_path.exists():
        return

    data = []
    with open(bed_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("track") or not line: continue
            cols = line.split("\t")
            try:
                start, stop = int(cols[1]), int(cols[2])
                data.append({
                    "Label": f"{start}->{stop}\n({cols[3]})",
                    "Count": int(cols[4]) if len(cols) > 4 else 1,
                    "Type": cols[3],
                })
            except (ValueError, IndexError):
                continue
    if not data:
        return

    top = sorted(data, key=lambda x: x["Count"], reverse=True)[:top_n]
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ["red" if "Deletion" in t["Type"] else ("blue" if "Duplication" in t["Type"] else "gray") for t in top]
    bars = ax.bar([t["Label"] for t in top], [t["Count"] for t in top], color=colors)
    for bar in bars:
        ax.annotate(f'{bar.get_height()}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom')
    ax.set_ylabel("Read Count"); ax.set_title(f"Top {top_n} Recombination Events\n{chrom_name}")
    plt.xticks(rotation=45, ha='right'); plt.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"top_events_{_safe_filename(chrom_name)}.png", dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / f"top_events_{_safe_filename(chrom_name)}.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)


# ==========================================
# 3. 核心 Worker：ViReMa + Visualize
# ==========================================
def worker_virema(args_tuple):
    (sample, virus, r1_str, r2_str, is_single, ref_fa, out_dir_str,
     seed, mindel, defuzz, virema_path, threads, resume, log_file) = args_tuple

    out_dir = Path(out_dir_str)
    sam_basename = f"{sample}_{virus}_ViReMa.sam"
    output_tag = f"{sample}_{virus}"

    # Resume: skip if BED output already exists (task previously completed)
    bed_dir = out_dir / "bed_results"
    if resume and bed_dir.exists() and any(bed_dir.glob("*.bed")):
        return True

    # Only create directory when task actually starts running
    out_dir.mkdir(parents=True, exist_ok=True)

    r1, r2 = Path(r1_str) if r1_str else None, Path(r2_str) if r2_str else None
    is_fasta = bool(r1 and any(ext in r1.name.lower() for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']))
    fmt_ext = ".fa" if is_fasta else ".fq"

    tmp_dir = out_dir / "virema_sandbox"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_input = tmp_dir / f"input_virema{fmt_ext}"
    if tmp_input.exists(): tmp_input.unlink()

    try:
        # A) Prepare input: for paired R1+R2, add /1 /2 suffix to avoid ViReMa KeyError
        need_suffix = not is_single and r1 and r2 and r1.exists() and r2.exists()

        def _open_src(src_path):
            if str(src_path).endswith('.gz'):
                import gzip; return gzip.open(src_path, 'rt')
            return open(src_path, 'r')

        if is_single and r1 and r1.exists():
            with _open_src(r1) as f_in, open(tmp_input, 'w') as f_out:
                shutil.copyfileobj(f_in, f_out)
        elif need_suffix:
            with open(tmp_input, 'w') as out_f:
                for src, tag in [(r1, '/1'), (r2, '/2')]:
                    if src and src.exists():
                        with _open_src(src) as in_f:
                            for line in in_f:
                                line = line.rstrip('\n\r')
                                if line.startswith('>'):
                                    # Ensure unique read names for ViReMa dict
                                    if not line.rstrip().endswith(tag):
                                        line = line.rstrip() + tag
                                out_f.write(line + '\n')
        else:
            return True
        if not tmp_input.exists() or tmp_input.stat().st_size == 0:
            return True

        # B) ViReMa
        virema_cmd = [
            sys.executable, str(virema_path),
            str(ref_fa), str(tmp_input), sam_basename,
            "--Seed", str(seed), "--MicroInDel_Length", str(mindel),
            "--Defuzz", str(defuzz), "--p", str(threads),
            "-BED", "-Overwrite",
            "--Output_Tag", output_tag, "--Output_Dir", str(tmp_dir),
        ] + (["-Fasta"] if is_fasta else [])
        run_cmd(virema_cmd, log_path=log_file)

        # C) Collect SAM
        produced_sam = tmp_dir / sam_basename
        if produced_sam.exists():
            shutil.move(str(produced_sam), str(out_dir / sam_basename))

        # D) Collect BED outputs
        dest_bed = out_dir / "bed_results"
        dest_bed.mkdir(exist_ok=True)
        bed_dir = tmp_dir / "BED_Files"
        if bed_dir.exists():
            for pattern in ["*.bed", "*.bedgraph", "*.BEDPE"]:
                for f in bed_dir.glob(pattern):
                    shutil.move(str(f), str(dest_bed / f.name))
        for pattern in ["*_Results.txt", "*_Insertions.txt", "*_Substitutions.txt", "*_Micro*.txt"]:
            for f in tmp_dir.glob(pattern):
                shutil.move(str(f), str(dest_bed / f.name))

        return True
    except Exception:
        with open(log_file, "a") as lf:
            lf.write(f"\n[Python Exception] {sample}/{virus}: {traceback.format_exc()}\n")
        return False
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ==========================================
# 4. DI-tector 替代引擎
# ==========================================
def worker_ditector(args_tuple):
    """DI-tector_06.py worker: BWA-based DVG detection (alternative to ViReMa)."""
    (sample, virus, r1_str, r2_str, is_single, ref_fa, out_dir_str,
     ditector_path, threads, tag) = args_tuple

    out_dir = Path(out_dir_str)
    out_dir.mkdir(parents=True, exist_ok=True)

    r1, r2 = Path(r1_str) if r1_str else None, Path(r2_str) if r2_str else None
    is_fasta = bool(r1 and any(ext in r1.name.lower() for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']))

    tmp_dir = out_dir / "ditector_sandbox"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_input = tmp_dir / "input.fa"
    if tmp_input.exists(): tmp_input.unlink()

    try:
        # Concat R1+R2 (decompress .gz)
        for src in [r1, r2]:
            if src and src.exists():
                if str(src).endswith('.gz'):
                    import gzip
                    with gzip.open(src, 'rb') as f_in, open(tmp_input, 'ab') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                else:
                    with open(src, 'rb') as f_in, open(tmp_input, 'ab') as f_out:
                        shutil.copyfileobj(f_in, f_out)
        if tmp_input.stat().st_size == 0:
            return True

        # BWA index if needed
        bwa_idx = Path(ref_fa).with_suffix(".amb")
        if not bwa_idx.exists():
            run_cmd(["bwa", "index", str(ref_fa)], check=False)

        # DI-tector
        ditector_cmd = [
            sys.executable, str(ditector_path),
            str(ref_fa), str(tmp_input),
            "-o", str(tmp_dir), "-t", tag, "-x", str(threads),
            "-s", "15", "-m", "25",
        ] + (["-f"] if is_fasta else [])
        run_cmd(ditector_cmd, check=False)

        # Collect output
        for f in tmp_dir.glob(f"{tag}_counts.txt"):
            shutil.copy(str(f), str(out_dir / f.name))
        for f in tmp_dir.glob(f"{tag}_DVG*"):
            shutil.copy(str(f), str(out_dir / f.name))
        return True
    except Exception:
        return False
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ==========================================
# 5. 主控枢纽
# ==========================================
class BatchViReMaPipeline:
    def __init__(self, args):
        self.args = args
        self.out = Path(args.output_dir)
        self.d_refs = self.out / "individual_refs"
        self.d_results = self.out / "virema_results"
        self.d_all_beds = self.out / "all_gathered_beds"
        for d in [self.d_refs, self.d_results, self.d_all_beds]:
            d.mkdir(parents=True, exist_ok=True)
        self.df, self.ref_map, self.sample_files, self.mol_types = None, {}, {}, {}

    def run(self):
        engine = "DI-tector" if self.args.ditector else "ViReMa"
        logger.info(f"DVG Pipeline — {engine} + Visualize + Circos")

        self.check_dependencies()
        self.load_summary()
        self.load_molecule_types()
        self.extract_virus_fastas()
        self.find_reads()
        self.print_dataset_statistics()

        if self.args.ditector:
            self.run_ditector()
        else:
            self.run_virema()
            self.run_arc_diagrams()

        self.run_auto_r_report()
        logger.info(f"Done -> {self.out}")

    def check_dependencies(self):
        missing = []
        if self.args.ditector:
            for exe in ("bwa", "samtools"):
                if not shutil.which(exe): missing.append(exe)
        else:
            if not shutil.which("bowtie-build") and not shutil.which("bowtie2-build"):
                missing.append("bowtie-build or bowtie2-build")

        vp = Path(self.args.virema_script)
        if vp.is_dir():
            py_files = list(vp.glob("ViReMa*.py"))
            if py_files: self.args.virema_script = str(py_files[0])
            else: missing.append(f"ViReMa.py not found in {vp}")
        elif not vp.exists():
            missing.append(f"Script not found: {vp}")
        if missing:
            logger.error("Missing: " + ", ".join(missing)); sys.exit(1)

    def load_summary(self):
        sep = "," if Path(self.args.summary).suffix.lower() == ".csv" else "\t"
        df = pd.read_csv(self.args.summary, sep=sep, dtype=str)
        df['Sample'] = df['Sample'].apply(lambda x: re.sub(r'(?i)(_clean|_trimmed|_filtered|_val|_fastp)+$', '', str(x).strip()))
        for col in ['Rep_Accession', 'Accession', 'Virus']:
            if col in df.columns: df['Virus_Acc'] = df[col]; break
        for col in ['Taxonomy', 'Adjusted_Species', 'Species_NCBI']:
            if col in df.columns: df['Taxonomy'] = df[col]; break
        if 'Virus_Acc' not in df.columns: logger.error("No accession column"); sys.exit(1)
        if 'Taxonomy' not in df.columns: df['Taxonomy'] = "Unannotated"
        if self.args.min_cov > 0:
            cov_cols = ['Rep_Coverage(%)', 'Coverage(%)', 'Coverage', 'Cov']
            cov_col = next((c for c in cov_cols if c in df.columns), None)
            if cov_col:
                df[cov_col] = df[cov_col].astype(str).str.replace('%', '').apply(pd.to_numeric, errors='coerce')
                n = len(df); df = df[df[cov_col] >= self.args.min_cov]
                if n - len(df) > 0: logger.info(f"Coverage filter: removed {n - len(df)} records")
        self.df = df.drop_duplicates(subset=['Sample', 'Virus_Acc'])

    def load_molecule_types(self):
        ref_info = getattr(self.args, 'ref_info', None)
        if not ref_info or not Path(ref_info).exists(): return
        with open(ref_info) as f:
            header = f.readline().strip().split('\t')
            try: idx_acc = header.index('Accession')
            except ValueError: return
            idx_mol = next((header.index(c) for c in ['Molecule_Type2', 'Molecule_type'] if c in header), None)
            if idx_mol is None: return
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) > idx_mol: self.mol_types[cols[idx_acc]] = cols[idx_mol]

    @staticmethod
    def _needs_polyA(mol_type: str) -> bool:
        return bool(mol_type and 'ssRNA(-)' in mol_type.upper())

    def extract_virus_fastas(self):
        target_viruses = set(self.df['Virus_Acc'].unique())
        seq_buf, vid_cur = [], None
        def _save():
            if vid_cur in target_viruses and seq_buf:
                ref_fa = self.d_refs / f"{safe_name(vid_cur)}.fasta"
                if not ref_fa.exists():
                    seq = "".join(seq_buf)
                    mol = self.mol_types.get(vid_cur, self.mol_types.get(vid_cur.split('.')[0], ''))
                    if self._needs_polyA(mol):
                        seq += "A" * 160
                        logger.info(f"  polyA padded: {vid_cur} ({mol})")
                    with open(ref_fa, "w") as f: f.write(f">{vid_cur}\n{seq}\n")
                self.ref_map[vid_cur] = str(ref_fa)
        with open(self.args.reference, "r") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"): _save(); vid_cur = line[1:].split()[0]; seq_buf = []
                else: seq_buf.append(line)
        _save()

    def find_reads(self):
        rd = Path(self.args.reads_dir)
        all_reads = [f for f in rd.rglob("*") if f.is_file()
                     and any(ext in f.name.lower() for ext in ['.fq', '.fastq', '.fa', '.fasta', '.gz'])]
        for sname in self.df['Sample'].unique():
            s_clean = sname.lower()
            matched = [f for f in all_reads
                       if re.search(r'\b' + re.escape(s_clean.replace('_', ' ')) + r'\b', f.name.lower().replace('_', ' ').replace('-', ' ').replace('.', ' '))
                       or f.name.lower().startswith(s_clean + "_") or f.name.lower().startswith(s_clean + ".")]
            matched = list(set(matched))
            if not matched: continue
            r1, r2 = None, None
            for f in matched:
                nl = f.name.lower()
                if any(x in nl for x in ['_r2', '_2.', '.r2', '_2_']): r2 = f
                elif any(x in nl for x in ['_r1', '_1.', '.r1', '_1_']): r1 = f
                elif not r1: r1 = f
            self.sample_files[sname] = {'r1': str(r1), 'r2': str(r2), 'is_single': not (r1 and r2)}

    def print_dataset_statistics(self):
        pe = sum(1 for s in self.sample_files.values() if not s['is_single'])
        logger.info(f"\n  Samples: {len(self.sample_files)} | Viruses: {len(self.ref_map)} | PE:{pe} SE:{len(self.sample_files)-pe}\n")

    def run_virema(self):
        tasks = []
        for _, r in self.df.iterrows():
            s, v, tax = r["Sample"], r["Virus_Acc"], r["Taxonomy"]
            if s not in self.sample_files or v not in self.ref_map: continue
            sf = self.sample_files[s]
            od = self.d_results / safe_name(tax) / f"{safe_name(s)}_{safe_name(v)}"
            tasks.append((s, v, sf['r1'], sf['r2'], sf['is_single'], self.ref_map[v], str(od),
                          self.args.seed, self.args.mindel, self.args.defuzz,
                          self.args.virema_script, self.args.threads, self.args.resume,
                          str(od / "virema.log")))

        logger.info(f"ViReMa: {len(tasks)} tasks (jobs={self.args.jobs})")
        with tqdm(total=len(tasks), desc="ViReMa", position=0) as pbar:
            with ProcessPoolExecutor(max_workers=self.args.jobs) as ex:
                futures = {ex.submit(worker_virema, t): t for t in tasks}
                ok = fail = 0
                for fut in as_completed(futures):
                    s, v = futures[fut][:2]
                    try:
                        if fut.result(): ok += 1
                        else: fail += 1; logger.warning(f"FAIL: {s}/{v}")
                    except Exception as e:
                        fail += 1; logger.error(f"CRASH: {s}/{v} - {e}")
                    pbar.update(1); pbar.set_postfix_str(f"OK:{ok} FAIL:{fail}")
        logger.info(f"ViReMa: {ok} OK, {fail} failed")

    def run_ditector(self):
        ditector_path = Path(self.args.virema_script).parent / "DI-tector_06.py"
        if not ditector_path.exists():
            logger.error(f"DI-tector not found: {ditector_path}"); return

        tasks = []
        for _, r in self.df.iterrows():
            s, v, tax = r["Sample"], r["Virus_Acc"], r["Taxonomy"]
            if s not in self.sample_files or v not in self.ref_map: continue
            sf = self.sample_files[s]
            od = self.d_results / safe_name(tax) / f"{safe_name(s)}_{safe_name(v)}"
            od.mkdir(parents=True, exist_ok=True)
            tasks.append((s, v, sf['r1'], sf['r2'], sf['is_single'], self.ref_map[v], str(od),
                          str(ditector_path), self.args.threads, f"{s}_{v}"))

        logger.info(f"DI-tector: {len(tasks)} tasks (jobs={self.args.jobs})")
        with tqdm(total=len(tasks), desc="DI-tector", position=0) as pbar:
            with ProcessPoolExecutor(max_workers=self.args.jobs) as ex:
                futures = {ex.submit(worker_ditector, t): t for t in tasks}
                ok = fail = 0
                for fut in as_completed(futures):
                    s, v = futures[fut][:2]
                    try:
                        if fut.result(): ok += 1
                        else: fail += 1; logger.warning(f"FAIL: {s}/{v}")
                    except Exception as e:
                        fail += 1; logger.error(f"CRASH: {s}/{v} - {e}")
                    pbar.update(1); pbar.set_postfix_str(f"OK:{ok} FAIL:{fail}")
        logger.info(f"DI-tector: {ok} OK, {fail} failed")

    def run_arc_diagrams(self):
        """Generate arc diagrams + bar charts for each virus (from visualize.py)."""
        if not HAS_MPL:
            logger.warning("matplotlib not available, skipping arc diagrams")
            return
        logger.info("Generating arc diagrams...")
        for virus_dir in self.d_results.iterdir():
            if not virus_dir.is_dir(): continue
            for sample_dir in virus_dir.iterdir():
                if not sample_dir.is_dir(): continue
                bed_files = list((sample_dir / "bed_results").glob("*Recombination_Results.bed"))
                for bed_file in bed_files[:1]:  # one per sample
                    plots_dir = sample_dir / "plots"
                    plot_arc_diagram(bed_file, sample_dir.name, plots_dir)
                    plot_top_events(bed_file, sample_dir.name, plots_dir)

    def run_auto_r_report(self):
        r_script = SCRIPT_DIR / "virema_summary_report.R"
        if not r_script.exists():
            r_script = Path(self.args.virema_script).parent / "virema_summary_report.R"
        if not r_script.exists():
            logger.warning("virema_summary_report.R not found, skipping Circos")
            return

        logger.info("Running Circos report (R)...")
        for pattern in ["bed_results/*.bed", "bed_results/*.bedgraph", "bed_results/*.txt"]:
            for f in self.d_results.rglob(pattern):
                shutil.copy(f, self.d_all_beds / f.name)

        r_out = self.out / "Summary_Analysis_Report"
        run_cmd(["Rscript", str(r_script), "-i", str(self.d_all_beds), "-o", str(r_out), "--auto_annotate"],
                log_path=str(self.out / "r_summary_run.log"), check=False)
        logger.info(f"Circos report -> {r_out}")


if __name__ == "__main__":
    if sys.platform == "win32":
        import colorama; colorama.init()

    parser = argparse.ArgumentParser(description="Batch ViReMa/DI-tector DVG Pipeline")
    parser.add_argument("-s", "--summary", required=True)
    parser.add_argument("-r", "--reference", required=True)
    parser.add_argument("-d", "--reads_dir", required=True)
    parser.add_argument("-v", "--virema_script", required=True)
    parser.add_argument("--ref_info", default=None)
    parser.add_argument("-o", "--output_dir", default="./virema_out")

    eng = parser.add_argument_group("Engine")
    eng.add_argument("--ditector", action="store_true", help="Use DI-tector (BWA) instead of ViReMa")

    vr = parser.add_argument_group("ViReMa params")
    vr.add_argument("--min_cov", type=float, default=0.0)
    vr.add_argument("--seed", type=int, default=25)
    vr.add_argument("--mindel", type=int, default=15)
    vr.add_argument("--defuzz", type=int, choices=[0, 3, 5], default=0)

    ctl = parser.add_argument_group("Concurrency")
    ctl.add_argument("-j", "--jobs", type=int, default=4)
    ctl.add_argument("-t", "--threads", type=int, default=8)
    ctl.add_argument("--resume", action="store_true")

    args = parser.parse_args()
    BatchViReMaPipeline(args).run()
