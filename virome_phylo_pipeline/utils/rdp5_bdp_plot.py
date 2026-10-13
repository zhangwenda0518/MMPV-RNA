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
    # 2026-10-12: 画布尺寸原先在 ggsave 里**写死 width=30, height=8** —— 30 in 宽
    # 的图塞进 9.8 in 通栏只剩 0.33 倍，源图 12 pt 文字排到 4 pt，**无法排版**。
    # 现暴露为 CLI（默认仍 30/8 → 既有调用行为逐字节不变）。
    p.add_argument("--width", type=float, default=30.0,
                   help="Output canvas width in inches (default 30; 组合图排版建议 12–16)")
    p.add_argument("--height", type=float, default=8.0,
                   help="Output canvas height in inches (default 8)")
    # 2026-10-12: R 主题字号是**绝对 pt**（MainFile.R: fontSizeMultiplier <- 2 →
    # title 40 / axis.text 30 / axis.title 20 pt），画布一小就溢出被裁。
    # 暴露为 --font-scale（默认 None = 沿用 R 里写死的 2 → 行为逐字节不变）。
    # 组合图排版建议：画布 16 in 宽时用 ~1.0–1.1。
    p.add_argument("--font-scale", type=float, default=None,
                   help="覆盖 R 的 fontSizeMultiplier (默认沿用 R 内写死的 2)")
    p.add_argument("--gene-font-size", type=float, default=None,
                   help="覆盖 R 的 fontSizeForGeneMap (默认沿用 R 内写死的 6)")
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
    # 2026-10-12: taxid 缺省是空串 → 标题会渲染成 "Virus [TaxID: ]"（空括号很难看）。
    # 改成有 taxid 才拼方括号。
    r_code = r_code.replace(
        'ggtitle(paste0(virusName, " [TaxID: ", taxId, "]"))',
        'ggtitle(if (nzchar(taxId)) paste0(virusName, " [TaxID: ", taxId, "]") else virusName)')
    # 2026-10-12: 可选覆盖 R 主题字号倍率（默认 None → 不改，沿用 R 内写死的 2）
    if args.font_scale is not None:
        r_code = r_code.replace('fontSizeMultiplier <- 2',
                                f'fontSizeMultiplier <- {args.font_scale}')
    if args.gene_font_size is not None:
        r_code = r_code.replace('fontSizeForGeneMap <- 6',
                                f'fontSizeForGeneMap <- {args.gene_font_size}')

    # 2026-10-12: 配色统一到项目 OKABE 调色板。原 MainFile.R 用的是
    # #E85E5D(红) / #C2DBD5(浅青) / #5F6AB1(蓝紫) / #87B4A9(青灰) —— 都不在
    # 项目调色板内，与其它 21 个面板的色语义对不上。与字体/尺寸一样在运行时
    # 做字符串覆盖，**不改第三方源码**（RDP5_RBDP_Rgrapher 保持原样）。
    _RDP5_COLORS = {
        '"#E85E5D"': '"#D55E00"',   # 断点 RBN 主填充      → 朱红
        '"#C2DBD5"': '"#56B4E9"',   # 99% CI 上界填充      → 天蓝
        '"#5F6AB1"': '"#0072B2"',   # 99% CI 下界填充      → 蓝
        '"#87B4A9"': '"#999999"',   # 99% CI 边界线        → 灰
    }
    for _old, _new in _RDP5_COLORS.items():
        r_code = r_code.replace(_old, _new)
    # 基因轨道青色渐变 → OKABE 蓝系
    r_code = r_code.replace(
        'c("#043F3D", "#085F63", "#0B7F89", "#13A3A9", "#17C6C9")',
        'c("#56B4E9", "#0072B2")')
    r_code = r_code.replace('geneMapDataFrame$color <- "#085F63"',
                            'geneMapDataFrame$color <- "#0072B2"')

    # Replace print(p) with ggsave
    r_code = r_code.replace(
        'print(p)',
        f'ggsave("{r_out_pdf.as_posix()}", plot=p, width={args.width}, height={args.height}, limitsize=FALSE)'
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
