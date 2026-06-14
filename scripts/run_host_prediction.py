#!/usr/bin/env python3
"""
02_run_host_prediction.py — 宿主预测及集成决策树 (ICTV > RNAVirHost > PhaBOX2)

Usage:
  python 02_run_host_prediction.py \
    -i vclust_centroids.fasta \
    --tax taxonomy_out/integrated/final_integrated_classification.tsv \
    -o host_out \
    -t 40 \
    --phabox-db ~/database/virus-db/phabox_db_v2_2
"""

import argparse
import os
import subprocess
import time
import sys
from datetime import datetime
from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent

def run_step(cmd, step_name):
    print(f"\n{'='*75}\n[{datetime.now().strftime('%H:%M:%S')}] {step_name}")
    head = ' '.join(cmd) if isinstance(cmd, list) else cmd
    print(f"[CMD] {head[:200]}\n{'='*75}")
    t0 = time.time()
    try:
        result = subprocess.run(cmd, shell=True if isinstance(cmd, str) else False,
                                check=True, capture_output=True, text=True)
        print(f"[OK] {step_name} — {time.time() - t0:.0f}s")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[FAIL] {step_name} — {time.time() - t0:.0f}s (rc={e.returncode})")
        if e.stderr:
            for l in e.stderr.strip().split('\n')[-6:]: print(f"  [stderr] {l}")
        return False

def check_file(path, min_size=10):
    return Path(path).is_file() and Path(path).stat().st_size > min_size

# ----------------- 数据准备与工具执行 -----------------

def run_tools(args):
    # 1. RNAVirHost
    rvh_csv = os.path.join(args.output_dir, "RVH_result", "result.csv")
    if not args.skip_rnavirhost:
        if not args.force and check_file(rvh_csv, 50):
            print(f"[SKIP] RNAVirHost — exists: {rvh_csv}")
        else:
            cmd1 = ["rnavirhost", "classify_order", "-i", args.input]
            run_step(cmd1, "RNAVirHost Classify Order")
            cmd2 = ["rnavirhost", "predict", "-i", args.input, "-o", os.path.join(args.output_dir, "RVH_result")]
            run_step(cmd2, "RNAVirHost Predict")

    # 2. PhaBOX2 CHERRY
    pb2_tsv = os.path.join(args.output_dir, "phabox2_output", "final_prediction", "cherry_prediction.tsv")
    if not args.skip_phabox:
        if not args.force and check_file(pb2_tsv, 50):
            print(f"[SKIP] PhaBOX2 — exists: {pb2_tsv}")
        else:
            cmd = ["phabox2", "--task", "cherry", "--dbdir", args.phabox_db,
                   "--outpth", os.path.join(args.output_dir, "phabox2_output"),
                   "--contigs", args.input, "--threads", str(args.threads), "--len", "500"]
            run_step(cmd, "PhaBOX2 CHERRY Host Prediction")

    # 3. ICTV C9 Lookup
    c9_tsv = os.path.join(args.output_dir, "C9_ICTV_result", "classification_result.tsv")
    if not args.skip_ictv:
        if not args.force and check_file(c9_tsv, 50):
            print(f"[SKIP] ICTV Lookup (C9) — exists: {c9_tsv}")
        elif check_file(args.tax):
            host_dir = os.path.join(args.output_dir, "C9_ICTV_result")
            os.makedirs(host_dir, exist_ok=True)
            cmd = f"python {SCRIPT_DIR}/C9_classify_contigs.py -i {args.tax} -f {args.input} --output_dir {host_dir} --prob_dir {args.prob_dir} --mode high"
            run_step(cmd, "ICTV Taxonomy Lookup (C9)")
        else:
            print("[WARN] Valid taxonomy TSV required for C9 ICTV lookup.")

    return rvh_csv, pb2_tsv, c9_tsv

# ----------------- 决策树与整合逻辑 -----------------

def normalize_c9(h):
    if pd.isna(h) or h in ['Unknown', 'None']: return 'Unknown'
    h = str(h)
    if h in ['Insecta', 'Arachnida', 'Aves', 'Human', 'Animal_other']: return 'Animal'
    if h in ['Oomycetes']: return 'Protist'
    return h

def parse_rvh(row):
    h = str(row.get('pred|L1', 'Unknown')).lower()
    ev = str(row.get('evidence', '')).lower()
    if h == 'unknown' or ev == 'unclassified' or pd.isna(h): return 'Unknown'
    if h == 'viridiplantae': return 'Plant'
    if h == 'fungi': return 'Fungi'
    if h in ['chordata', 'invertebrate', 'metazoa', 'animal']: return 'Animal'
    if h == 'bacteria': return 'Bacteria'
    return 'Unknown'

def parse_pb2(row):
    lineage = (str(row.get('Host_NCBI_lineage', '')) + "|" +
               str(row.get('Host_GTDB_lineage', '')) + "|" +
               str(row.get('Host', ''))).lower()
    if 'bacteria' in lineage: return 'Bacteria'
    if 'archaea' in lineage: return 'Archaea'
    if 'streptophyta' in lineage or 'viridiplantae' in lineage or 'plant' in lineage: return 'Plant'
    if 'fungi' in lineage: return 'Fungi'
    if any(k in lineage for k in ['metazoa', 'animal', 'chordata', 'arthropoda', 'insecta']): return 'Animal'
    return 'Unknown'

def decision_tree(row):
    tax_class = str(row.get('Class', '')).lower()
    if tax_class in ['leviviricetes', 'caudoviricetes', 'vidaverviricetes', 'faserviricetes']:
        return 'Bacteria', 'Rule_Class_Taxonomy'

    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))
    h_rvh  = parse_rvh(row)
    h_pb2  = parse_pb2(row)

    # 1. 全部未知
    if h_ictv == 'Unknown' and h_rvh == 'Unknown' and h_pb2 == 'Unknown':
        return 'Unknown', 'Unassigned'
    # 2. ICTV 与 RNAVirHost 一致
    if h_ictv == h_rvh and h_ictv != 'Unknown':
        return h_ictv, 'ICTV_RVH_Agree'
    # 3. 分歧时 PhaBOX2 决胜
    if h_pb2 != 'Unknown':
        if h_pb2 == h_ictv: return h_ictv, 'ICTV_PB2_Agree_Tiebreaker'
        if h_pb2 == h_rvh:  return h_rvh,  'RVH_PB2_Agree_Tiebreaker'
    # 4. 兜底: ICTV > RNAVirHost > PhaBOX2
    if h_ictv != 'Unknown': return h_ictv, 'ICTV_Preferred'
    if h_rvh  != 'Unknown': return h_rvh,  'RVH_Preferred'
    if h_pb2  != 'Unknown': return h_pb2,  'PB2_Preferred'
    return 'Unknown', 'Fallback'

def run_ensemble(args, rvh_csv, pb2_tsv, c9_tsv):
    print(f"\n{'='*75}\n[Ensemble] Running Decision Tree Integration\n{'='*75}")
    df = pd.read_csv(args.tax, sep='\t') if check_file(args.tax) else pd.DataFrame(columns=['contig_id', 'Class'])

    if check_file(c9_tsv):
        c9 = pd.read_csv(c9_tsv, sep='\t')[['contig_id', 'Predicted_Host']].rename(columns={'Predicted_Host': 'Host_ICTV'})
        df = df.merge(c9, on='contig_id', how='left')

    if check_file(rvh_csv):
        rvh = pd.read_csv(rvh_csv)
        id_col = rvh.columns[0]
        if 'Unnamed' in id_col or 'y|virus order' in rvh.columns:
            rvh.rename(columns={id_col: 'contig_id'}, inplace=True)
            use_cols = ['contig_id'] + [c for c in ['pred|L1', 'pred|L2', 'evidence'] if c in rvh.columns]
            df = df.merge(rvh[use_cols], on='contig_id', how='left')

    if check_file(pb2_tsv):
        pb2 = pd.read_csv(pb2_tsv, sep='\t').rename(columns={'Accession': 'contig_id'})
        use_cols = [c for c in ['contig_id', 'Host', 'Host_NCBI_lineage', 'Host_GTDB_lineage'] if c in pb2.columns]
        df = df.merge(pb2[use_cols], on='contig_id', how='left')

    results = df.apply(decision_tree, axis=1)
    df['Final_Host'], df['Decision_Method'] = zip(*results)

    out_tsv = os.path.join(args.output_dir, "ensemble_host_summary.tsv")
    df.to_csv(out_tsv, sep='\t', index=False)

    print("\n  [Summary]")
    for host, cnt in df['Final_Host'].value_counts().items():
        print(f"    {host:<15s}: {cnt:>6} contigs")

    print("\n  [Splitting FastA]")
    host_dict = dict(zip(df['contig_id'].astype(str), df['Final_Host']))
    out_fastas = {}
    with open(args.input, 'r') as f:
        curr_id, curr_seq = None, []
        for line in f:
            if line.startswith('>'):
                if curr_id:
                    out_fastas.setdefault(host_dict.get(curr_id, 'Unknown'), []).append(f">{curr_id}\n{''.join(curr_seq)}")
                curr_id = line[1:].split()[0]
                curr_seq = []
            else: curr_seq.append(line)
        if curr_id:
            out_fastas.setdefault(host_dict.get(curr_id, 'Unknown'), []).append(f">{curr_id}\n{''.join(curr_seq)}")

    fasta_dir = os.path.join(args.output_dir, "host_classified_fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    for host, seqs in out_fastas.items():
        with open(os.path.join(fasta_dir, f"{host}.classified.fasta"), 'w') as f: f.write("".join(seqs))
    print(f"    Saved into: {fasta_dir}/\n")

def main():
    p = argparse.ArgumentParser(description="Ensemble Host Prediction Workflow")
    p.add_argument("-i", "--input", required=True, help="Input FASTA file")
    p.add_argument("--tax", required=True, help="Integrated taxonomy TSV (from 01_run_taxonomy.py)")
    p.add_argument("-o", "--output-dir", default="host_out", help="Output directory")
    p.add_argument("-t", "--threads", type=int, default=40, help="Threads to use")
    p.add_argument("--phabox-db", default=os.path.expanduser("~/database/virus-db/phabox_db_v2_2"), help="PhaBOX2 DB")
    p.add_argument("--prob-dir", default="cross_analysis",
                   help="ICTV 宿主概率表目录 (C9_classify_contigs 使用, 默认: cross_analysis/)")

    p.add_argument("--mode", default="all",
                   choices=["all", "ICTV", "RNAVirHost", "PhaBOX2"],
                   help="运行模式: all(三工具集成) ICTV RNAVirHost PhaBOX2 (默认: all)")
    p.add_argument("--skip-rnavirhost", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--skip-phabox", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--skip-ictv", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("-f", "--force", action="store_true", help="Force re-run all stages")
    args = p.parse_args()

    # mode 自动设置 skip 标志 (--skip-* 显式传入时优先)
    if args.mode != "all":
        args.skip_rnavirhost = args.skip_rnavirhost or (args.mode != "RNAVirHost")
        args.skip_phabox     = args.skip_phabox     or (args.mode != "PhaBOX2")
        args.skip_ictv       = args.skip_ictv       or (args.mode != "ICTV")

    if not os.path.isfile(args.input): sys.exit(f"ERROR: FastA file not found: {args.input}")
    os.makedirs(args.output_dir, exist_ok=True)

    rvh_csv, pb2_tsv, c9_tsv = run_tools(args)
    run_ensemble(args, rvh_csv, pb2_tsv, c9_tsv)

if __name__ == "__main__":
    main()
