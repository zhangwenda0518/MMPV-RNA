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
            # Map NCBI runinfo size_MB → GB
            if "size_MB" in df.columns:
                df["FileSize_GB"] = df["size_MB"].apply(
                    lambda x: f"{float(x)/1024:.1f}" if pd.notna(x) and str(x).replace('.','').isdigit() else "")
        return {"ok": True, "df": df if not df.empty else pd.DataFrame(), "error": ""}
    except Exception as e:
        return {"ok": False, "df": pd.DataFrame(), "error": str(e)}


def _extract_gsa_filesizes(outdir: str, df) -> None:
    """Extract GSA file sizes: preferred from Run Excel, fallback to HTML."""
    import re, os
    if df is None or df.empty or "Run" not in df.columns:
        return
    df["FileSize_GB"] = ""
    df["FileCount"] = ""

    # Best source: Run Excel files (1_xls_cache)
    xls_dir = os.path.join(outdir, "GSA_Results", "1_xls_cache")
    if os.path.isdir(xls_dir):
        import pandas as pd
        import re
        run_sizes = {}
        for fname in os.listdir(xls_dir):
            if not fname.endswith(".xlsx"):
                continue
            xls_path = os.path.join(xls_dir, fname)
            try:
                xls = pd.ExcelFile(xls_path)
                if "Run" not in xls.sheet_names:
                    continue
                df_run = pd.read_excel(xls, sheet_name="Run")
                for _, rrow in df_run.iterrows():
                    acc = str(rrow.get("Accession", "")).strip()
                    if not acc:
                        continue
                    # Extract file sizes from Read filename columns
                    # Format: "CRR1126132_f1.fq.gz (1529402668 bytes)"
                    total_bytes = 0
                    file_count = 0
                    for col in df_run.columns:
                        cl = col.strip().lower()
                        if "read filename" in cl and "md5" not in cl:
                            val = str(rrow[col]).strip() if pd.notna(rrow[col]) else ""
                            m = re.search(r'\((\d+)\s*bytes?\)', val, re.I)
                            if m:
                                total_bytes += int(m.group(1))
                                file_count += 1
                    if file_count > 0:
                        total_gb = total_bytes / 1073741824
                        label = f"{total_gb:.1f}"
                        run_sizes[acc] = (label, file_count)
            except Exception:
                pass

        for idx, row in df.iterrows():
            run = str(row["Run"]).strip()
            if run in run_sizes:
                label, n = run_sizes[run]
                df.at[idx, "FileSize_GB"] = label
                df.at[idx, "FileCount"] = str(n)
        if run_sizes:
            return  # Excel data is sufficient

    # Fallback: HTML cache
    cache_dir = os.path.join(outdir, "GSA_Results", "0_web_cache")
    if not os.path.isdir(cache_dir):
        return

    for idx, row in df.iterrows():
        run = str(row["Run"]).strip()
        if not run:
            continue
        # Try exact run ID, then try matching CRR prefix
        for fname in [f"{run}.html"] + [
            f for f in os.listdir(cache_dir)
            if f.startswith(run) and f.endswith(".html")
        ]:
            html_path = os.path.join(cache_dir, fname)
            if not os.path.isfile(html_path):
                continue
            try:
                with open(html_path, "r", encoding="utf-8") as f:
                    html = f.read()

                # Find file sizes in bytes: "filename.fastq.gz (1234567890 bytes)" or similar
                # Extract all (number bytes/GB/MB/KB) patterns
                byte_sizes = re.findall(
                    r'\((\d+)\s*bytes?\)', html, re.I)
                if not byte_sizes:
                    # Try GB/MB format
                    gb_sizes = re.findall(
                        r'\(([\d.]+)\s*GB?\)', html, re.I)
                    byte_sizes = [str(float(s) * 1073741824) for s in gb_sizes]
                if not byte_sizes:
                    mb_sizes = re.findall(
                        r'\(([\d.]+)\s*MB?\)', html, re.I)
                    byte_sizes = [str(float(s) * 1048576) for s in mb_sizes]

                if byte_sizes:
                    # Filter: only file sizes (not html content sizes)
                    # Usually file sizes are > 100000 bytes
                    size_bytes = [int(s) for s in byte_sizes
                                  if int(s) > 100000 and int(s) < 100000000000]
                    if len(size_bytes) >= 1:
                        total_bytes = sum(size_bytes)
                        total_mb = total_bytes / 1048576
                        # Check if paired (2 files roughly same size)
                        if len(size_bytes) == 2 and abs(size_bytes[0] - size_bytes[1]) < size_bytes[0] * 0.3:
                            df.at[idx, "FileSize_GB"] = f"{total_mb:.0f} MB (PE)"
                        else:
                            df.at[idx, "FileSize_GB"] = f"{total_mb:.0f} MB"
                        df.at[idx, "FileCount"] = str(len(size_bytes))
                        break
            except Exception:
                pass

        # Also try Excel cache as fallback
        if not df.at[idx, "FileSize_GB"]:
            xls_dir = os.path.join(outdir, "GSA_Results", "1_xls_cache")
            if os.path.isdir(xls_dir):
                try:
                    import pandas as pd
                    for xf in os.listdir(xls_dir):
                        if not xf.endswith(".xlsx"):
                            continue
                        xls = pd.ExcelFile(os.path.join(xls_dir, xf))
                        if "Run" not in xls.sheet_names:
                            continue
                        df_run = pd.read_excel(xls, sheet_name="Run")
                        for _, rrow in df_run.iterrows():
                            if str(rrow.get("Accession", "")).strip() == run:
                                dl_cols = [c for c in df_run.columns
                                           if "download read file" in c.strip().lower()
                                           and "md5" not in c.strip().lower()]
                                urls = [str(rrow[c]).strip() for c in dl_cols
                                        if pd.notna(rrow[c]) and str(rrow[c]).strip()]
                                if urls:
                                    df.at[idx, "FileSize_GB"] = f"{len(urls)} file(s)"
                                    df.at[idx, "FileCount"] = str(len(urls))
                    break
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
