#!/usr/bin/env python3
"""
rdp5_bdp_plot.py — Wrapper for RDP5_RBDP_Rgrapher
==================================================
一键生成可发表断点分布图。

用法:
  python rdp5_bdp_plot.py -b breakpoint.csv -g orf.csv -p positions.csv -n "Virus Name"

需要 R 和 RDP5_RBDP_Rgrapher（同目录下）。
"""

import argparse, os, shutil, subprocess, sys, tempfile
from pathlib import Path
from datetime import datetime

SCRIPT_DIR = Path(__file__).parent
RGRAPHER_DIR = SCRIPT_DIR / "RDP5_RBDP_Rgrapher"
MAINFILE = RGRAPHER_DIR / "MainFile.R"
RSCRIPT = shutil.which("Rscript") or r"C:\Program Files\R\R-4.5.3\bin\Rscript.exe"


def parse_args():
    p = argparse.ArgumentParser(description="RDP5 Breakpoint Distribution Plot")
    p.add_argument("-b", "--breakpoint", required=True, help="RDP5 Breakpoint Distribution CSV")
    p.add_argument("-g", "--gene-map", required=True, help="ORF Coords CSV (*ORFCoords.csv)")
    p.add_argument("-p", "--positions", required=True, help="Breakpoint Positions CSV (*BreakpointPositions.csv)")
    p.add_argument("-n", "--name", default="Virus", help="Virus name for plot title")
    p.add_argument("--taxid", default="", help="NCBI TaxID (optional)")
    p.add_argument("-o", "--output", default=None, help="Output PDF path (default: <name>_breakpoint_plot.pdf)")
    p.add_argument("--rscript", default=RSCRIPT, help="Path to Rscript executable")
    return p.parse_args()


def check_r_packages():
    """Check if required R packages are installed."""
    required = ["ggplot2", "tidyverse", "scales", "ggfittext", "gridExtra"]
    for pkg in required:
        result = subprocess.run(
            [RSCRIPT, "-e", f"library({pkg})"],
            capture_output=True, text=True, encoding='utf-8', errors='replace'
        )
        if result.returncode != 0:
            print(f"Installing missing R package: {pkg}")
            subprocess.run(
                [RSCRIPT, "-e", f"install.packages('{pkg}', repos='https://cran.r-project.org')"],
                capture_output=True
            )


def run():
    args = parse_args()
    
    if not MAINFILE.exists():
        sys.exit(f"RDP5_RBDP_Rgrapher not found at {RGRAPHER_DIR}")
    
    # Check R
    try:
        subprocess.run([RSCRIPT, "--version"], capture_output=True, check=True)
    # 2026-09-15: 原为裸 except, 连 KeyboardInterrupt/SystemExit 都吞。
    except (OSError, subprocess.CalledProcessError) as _e:
        sys.exit(f"R not found ({_e}). Install R and ensure Rscript is in PATH.")
    
    check_r_packages()
    
    # Copy input files to temp (R can't handle unicode paths)
    work_dir = Path(tempfile.mkdtemp(prefix="rdp5_bdp_"))
    bdp_csv = work_dir / "breakpoint_data.csv"
    orf_csv = work_dir / "breakpoint_data.csvORFCoords.csv"
    pos_csv = work_dir / "breakpoint_data.csvBreakpointPositions.csv"
    
    shutil.copy(args.breakpoint, bdp_csv)
    shutil.copy(args.gene_map, orf_csv)
    shutil.copy(args.positions, pos_csv)
    
    # Output path (will copy from temp after R finishes)
    r_out_pdf = work_dir / "output.pdf"
    out_pdf = Path(args.output) if args.output else (SCRIPT_DIR / f"{args.name}_breakpoint_plot.pdf")
    
    # Read and modify MainFile.R
    with open(MAINFILE, 'r') as f:
        r_code = f.read()
    
    # Replace library(RCurl) - not needed
    r_code = r_code.replace("library(RCurl)\n", "")
    
    # Replace file paths
    r_code = r_code.replace(
        'breakpointData <- read.csv("example_files/Sarbecovirus.csv")',
        f'breakpointData <- read.csv("{bdp_csv.as_posix()}")'
    )
    r_code = r_code.replace(
        'geneMap <- read.csv("example_files/Sarbecovirus.csvORFCoords.csv")',
        f'geneMap <- read.csv("{orf_csv.as_posix()}")'
    )
    r_code = r_code.replace(
        'breakpointDotPos <- read.csv("example_files/Sarbecovirus.csvBreakpointPositions.csv")',
        f'breakpointDotPos <- read.csv("{pos_csv.as_posix()}")'
    )
    
    # Add missing "Bottom" column (RDP5 doesn't output this)
    r_code = r_code.replace(
        'bottom <- breakpointData$Bottom',
        'breakpointData$Bottom <- pmin(breakpointData$Recombination.breakpoint.number..200nt.win., breakpointData$Upper.99..CI, breakpointData$Lower.99..CI, breakpointData$Upper.95..CI, breakpointData$Lower.95..CI, na.rm=TRUE)\nbottom <- breakpointData$Bottom'
    )
    r_code = r_code.replace('virusName <- "Sarbecovirus"', f'virusName <- "{args.name}"')
    r_code = r_code.replace('taxId <- "694014"', f'taxId <- "{args.taxid}"')
    
    # Replace print(p) with ggsave
    r_code = r_code.replace(
        'print(p)',
        f'ggsave("{r_out_pdf.as_posix()}", plot=p, width=30, height=8, limitsize=FALSE)'
    )
    r_code = r_code.replace(
        'print("Graph generated successfully without overlapping gene names.")',
        'cat("OK\\n")'
    )
    r_script = work_dir / "run.R"
    with open(r_script, 'w') as f:
        f.write(r_code)
    
    # Run R
    print(f"Generating breakpoint plot for {args.name}...")
    result = subprocess.run([RSCRIPT, str(r_script)], capture_output=True, text=True,
                           encoding='latin-1', errors='replace')
    
    if r_out_pdf.exists():
        shutil.copy(r_out_pdf, out_pdf)
        print(f"✓ Saved: {out_pdf}  ({r_out_pdf.stat().st_size//1024} KB)")
    
    # Cleanup
    shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    run()
