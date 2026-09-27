#!/usr/bin/env python3
"""
Memory-safe data loading utilities for large TSV/CSV files.

Provides fallback strategies when in-memory loading would exceed available RAM.
"""

import os
import sys
from pathlib import Path


def available_memory_mb():
    """Estimate available system memory in MB. Returns None on Windows without psutil."""
    try:
        import psutil
        return psutil.virtual_memory().available / (1024 * 1024)
    except ImportError:
        pass
    # Linux fallback via /proc/meminfo
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1024
    except Exception:
        pass
    return None


def estimate_file_memory_mb(path, sep="\t"):
    """Estimate memory (MB) needed to load a TSV/CSV fully in pandas/polars.
    
    Rough heuristic: file_size * 5 (parsed overhead for string columns).
    """
    if not Path(path).exists():
        return 0
    size_mb = Path(path).stat().st_size / (1024 * 1024)
    return size_mb * 5


def safe_read_tsv(path, sep="\t", max_memory_ratio=0.5, chunk_size=1000):
    """Read a TSV file, falling back to chunked reading if memory is tight.
    
    Args:
        path: Path to TSV/CSV file
        sep: Column separator (default: tab)
        max_memory_ratio: Max fraction of available memory to use (default: 0.5)
        chunk_size: Rows per chunk when falling back to chunked mode
    
    Returns:
        pandas DataFrame, or None if file not found
    
    Raises:
        MemoryError: if the file is too large to fit even with chunking
    """
    try:
        import pandas as pd
    except ImportError:
        print("[WARN] pandas not installed, cannot read TSV", file=sys.stderr)
        return None

    if not Path(path).exists():
        return None

    est_mb = estimate_file_memory_mb(path, sep)
    avail_mb = available_memory_mb()

    if avail_mb and est_mb > avail_mb * max_memory_ratio:
        print(f"[INFO] Large file detected ({est_mb:.0f}MB est, {avail_mb:.0f}MB avail). "
              f"Using chunked reading (chunk_size={chunk_size})...", file=sys.stderr)
        chunks = []
        for chunk in pd.read_csv(path, sep=sep, chunksize=chunk_size, low_memory=False):
            chunks.append(chunk)
        return pd.concat(chunks, ignore_index=True) if chunks else None

    try:
        return pd.read_csv(path, sep=sep, low_memory=False)
    except MemoryError:
        print(f"[WARN] MemoryError on {path}, retrying with chunked reading...", file=sys.stderr)
        chunks = []
        for chunk in pd.read_csv(path, sep=sep, chunksize=chunk_size, low_memory=False):
            chunks.append(chunk)
        return pd.concat(chunks, ignore_index=True) if chunks else None


def safe_read_polars(path, sep="\t", max_memory_ratio=0.5):
    """Read a TSV file with Polars, using lazy scanning if memory is tight.
    
    Args:
        path: Path to TSV/CSV file
        sep: Column separator (default: tab)
        max_memory_ratio: Max fraction of available memory to use (default: 0.5)
    
    Returns:
        polars DataFrame (lazy if memory is tight, eager otherwise), or None
    """
    try:
        import polars as pl
    except ImportError:
        print("[WARN] polars not installed, cannot read TSV", file=sys.stderr)
        return None

    if not Path(path).exists():
        return None

    est_mb = estimate_file_memory_mb(path, sep)
    avail_mb = available_memory_mb()

    if avail_mb and est_mb > avail_mb * max_memory_ratio:
        print(f"[INFO] Large file detected ({est_mb:.0f}MB est). Using Polars lazy scan.", file=sys.stderr)
        return pl.scan_csv(path, separator=sep)

    return pl.read_csv(path, separator=sep)


# ── Convenience: read summary TSV (auto-detects pandas vs polars) ──
def read_summary_safe(path, use_polars=False):
    """Read a pipeline summary TSV with memory safety.
    
    Args:
        path: Path to summary TSV
        use_polars: If True, prefer polars; otherwise use pandas
    
    Returns:
        DataFrame (pandas or polars), or None if file not found
    """
    if use_polars:
        return safe_read_polars(path)
    return safe_read_tsv(path)
