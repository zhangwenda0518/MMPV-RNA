"""Bridge module to call gsa_sra.search.py / gsa_sra.info.py from the GUI."""

import sys
import os
import importlib.util
from typing import Optional, List

_GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PROJECT_DIR = os.path.dirname(_GUI_DIR)
_PIPELINE_DIR = os.path.join(_PROJECT_DIR, "public_metadata_pipeline")
_SEARCH_PATH = os.path.join(_PIPELINE_DIR, "gsa_sra.search.py")


def _load_search_module():
    spec = importlib.util.spec_from_file_location("gsa_sra_search", _SEARCH_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def search_sra(query: str, source: str = "", outdir: str = "",
               detailed: bool = True, ncbi_api: str = None) -> dict:
    import pandas as pd
    # Normalize "All" to no filter
    src = None if (not source or source.strip().lower() == "all") else source
    try:
        mod = _load_search_module()
        engine = mod.SRAEngine(query, src, outdir or os.getcwd(),
                               detailed=detailed, ncbi_api=ncbi_api)
        df = engine.fetch_runinfo()
        if not df.empty:
            df["Database"] = "SRA"
            # Map NCBI runinfo size columns
            for src_col, dst_col in [("size_MB", "FileSize_MB"),
                                       ("bases", "Bases"),
                                       ("spots", "Spots")]:
                if src_col in df.columns:
                    df[dst_col] = df[src_col]
        return {"ok": True, "df": df if not df.empty else pd.DataFrame(), "error": ""}
    except Exception as e:
        return {"ok": False, "df": pd.DataFrame(), "error": str(e)}


def _extract_gsa_filesizes(outdir: str, df) -> None:
    """Extract file sizes for GSA runs from cached HTML and add to DataFrame."""
    import re, os
    if df is None or df.empty or "Run" not in df.columns:
        return
    cache_dir = os.path.join(outdir, "GSA_Results", "0_web_cache")
    if not os.path.isdir(cache_dir):
        return
    df["FileSize_MB"] = ""
    df["FileCount"] = ""
    for idx, row in df.iterrows():
        run = str(row["Run"]).strip()
        if not run:
            continue
        html_path = os.path.join(cache_dir, f"{run}.html")
        if not os.path.isfile(html_path):
            continue
        try:
            with open(html_path, "r", encoding="utf-8") as f:
                html = f.read()
            # Extract file sizes: pattern matches file name + size cell
            sizes = re.findall(
                r'<td[^>]*>\s*([\d.,]+\s*(?:bytes?|KB|MB|GB|TB|GiB|MiB|KiB)?)\s*</td>',
                html, re.I)
            if not sizes:
                # Try number patterns with file extensions nearby
                sizes = re.findall(
                    r'(?:fq|fastq|bam)\.(?:gz)?.*?</td>\s*<td[^>]*>\s*([\d.,]+\s*\S*)',
                    html, re.I)
            if sizes:
                df.at[idx, "FileSize_MB"] = " | ".join(sizes[:4])
                df.at[idx, "FileCount"] = str(len(sizes))
        except Exception:
            pass


def search_gsa(query: str, source: str = "", outdir: str = "",
               detailed: bool = True) -> dict:
    import pandas as pd
    try:
        mod = _load_search_module()
        engine = mod.GSAEngine(query, source or None, outdir or os.getcwd(),
                               detailed=detailed)
        df = engine.fetch_gsa()
        if not df.empty:
            df["Database"] = "GSA"
            _extract_gsa_filesizes(outdir, df)
        return {"ok": True, "df": df if not df.empty else pd.DataFrame(), "error": ""}
    except Exception as e:
        return {"ok": False, "df": pd.DataFrame(), "error": str(e)}


def search_both(query: str, source: str = "", outdir: str = "",
                detailed: bool = True, ncbi_api: str = None) -> dict:
    result = {"ok": True, "sra": None, "gsa": None}
    result["sra"] = search_sra(query, source, outdir, detailed, ncbi_api)
    result["gsa"] = search_gsa(query, source, outdir, detailed)
    if not result["sra"]["ok"] and not result["gsa"]["ok"]:
        result["ok"] = False
    return result


# ── gsa_sra.info.py deep extraction ──────────────────

def deep_extract(run_ids: List[str], outdir: str = "",
                 deepseek_api: str = None, model: str = "deepseek-chat",
                 ncbi_api: str = None) -> dict:
    """Run gsa_sra.info.py on a list of Run IDs via subprocess.
    Returns {'ok': bool, 'tsv_path': str, 'full_path': str, 'error': str}."""
    import subprocess, tempfile
    if not run_ids:
        return {"ok": False, "tsv_path": "", "full_path": "", "error": "No Run IDs"}

    # Write Run IDs to temp file
    fd, id_file = tempfile.mkstemp(suffix=".txt")
    with os.fdopen(fd, "w") as f:
        f.write("\n".join(run_ids))

    out = outdir or tempfile.mkdtemp(prefix="info_extract_")
    os.makedirs(out, exist_ok=True)

    script = os.path.join(_PIPELINE_DIR, "gsa_sra.info.py")
    cmd = [sys.executable, script, "-i", id_file, "-o", out, "-m", "api"]

    env = os.environ.copy()
    if deepseek_api:
        env["DEEPSEEK_API_KEY"] = deepseek_api
    if ncbi_api:
        env["NCBI_API_KEY"] = ncbi_api

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
            encoding="utf-8", errors="replace",
            cwd=_PIPELINE_DIR, env=env)
        if result.returncode != 0:
            return {"ok": False, "tsv_path": "", "full_path": "",
                    "error": result.stderr[-500:] or "Unknown error"}

        tsv = os.path.join(out, "Global_Unified_Metadata_Core13.tsv")
        full = os.path.join(out, "Global_Unified_Metadata_Full.tsv")
        if os.path.isfile(tsv):
            return {"ok": True, "tsv_path": tsv, "full_path": full, "error": ""}
        return {"ok": False, "tsv_path": "", "full_path": "",
                "error": f"Output not found: {tsv}"}
    except subprocess.TimeoutExpired:
        return {"ok": False, "tsv_path": "", "full_path": "",
                "error": "Deep extraction timed out (>10min)"}
    except Exception as e:
        return {"ok": False, "tsv_path": "", "full_path": "",
                "error": str(e)}
