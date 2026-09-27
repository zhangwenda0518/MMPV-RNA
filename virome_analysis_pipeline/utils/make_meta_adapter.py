#!/usr/bin/env python3
"""make_meta_adapter.py - Build a Run+Taxonomy adapter table for virus_metadata_plot.py.

virus_metadata_plot.py expects the virus detection table to carry 'Run' and
'Taxonomy' columns, while filter_summary.py emits 'Sample' + species columns.
This adapter extracts (Sample, Taxonomy) from a high-confidence summary so the
metadata association stage can consume it directly.

Usage:
    python make_meta_adapter.py <high_conf.summary.tsv> <adapter_out.tsv>
"""
import sys

import pandas as pd


def main():
    if len(sys.argv) != 3:
        print("usage: make_meta_adapter.py <high_conf.summary.tsv> <adapter_out.tsv>")
        return 2
    src, dst = sys.argv[1], sys.argv[2]

    df = pd.read_csv(src, sep="\t", dtype=str)
    df.columns = [str(c).replace("\ufeff", "").strip() for c in df.columns]

    sample_col = next((c for c in df.columns if c.lower() == "sample"),
                      df.columns[0])
    tax_col = next(
        (c for c in df.columns
         if c in ("Adjusted_Species", "Species_NCBI", "Species_ICTV")),
        df.columns[1])

    out = df[[sample_col, tax_col]].rename(columns={sample_col: "Sample",
                                                    tax_col: "Taxonomy"})
    out.to_csv(dst, sep="\t", index=False)
    print(f"[meta_adapter] {len(out)} rows <- {src} ({sample_col}, {tax_col}) -> {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
