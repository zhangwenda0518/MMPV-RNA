#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""本地核心文件哈希 + 与服务器清单对比 (服务器端 md5 清单: %TMP%/server_core_md5.txt)"""
import hashlib, os, sys
from pathlib import Path

R = Path(r"D:\桌面\延伸基因组\MMPV-RNA")
TMP = os.environ.get("TMP", r"C:\Users\17711\AppData\Local\Temp")
SRV = Path(TMP) / "server_core_md5.txt"
LOC = Path(TMP) / "local_core_md5.txt"

DIRS = ["virome_discovery_pipeline", "virome_analysis_pipeline", "public_metadata_pipeline",
        "virome_submission_pipeline", "virome_phylo_pipeline", "metadata_gui", "submission_gui"]
FILES = ["pipeline_config.yaml", "DATABASE_SETUP.md", "SOFTWARE_VERSIONS.txt", "README.md", "pixi.toml"]
EXTS = {".py", ".R", ".yaml", ".md", ".pl", ".toml", ".txt", ".sh"}
EXCL = ("__pycache__", "archive", ".bak", ".orig", "cache/", "phylo_results", "_scratch",
        "_archive", "codon_usage", ".log", "/build/", "/dist/", "node_modules",
        "server_checks", "_dev_smoke", "_test_", "results_PSTVd", "_backup_utils_flatten")

def kept(p: str) -> bool:
    return not any(x in p for x in EXCL)

rows = []
for d in DIRS:
    for p in (R / d).rglob("*"):
        if p.is_file() and p.suffix in EXTS and kept(str(p).replace("\\", "/")):
            rows.append(p)
for f in FILES:
    p = R / f
    if p.exists():
        rows.append(p)

loc = {}
for p in rows:
    h = hashlib.md5(p.read_bytes()).hexdigest()
    rel = str(p.relative_to(R)).replace("\\", "/")
    loc[rel] = h
    print(h, rel, file=open(LOC, "a", encoding="utf-8"))
LOC.write_text("\n".join(f"{h}  {r}" for r, h in sorted(loc.items())) + "\n", encoding="utf-8")

srv = {}
for line in SRV.read_text(encoding="utf-8", errors="ignore").splitlines():
    parts = line.strip().split(None, 1)
    if len(parts) == 2:
        srv[parts[1].lstrip("*").strip()] = parts[0]

common_diff = [r for r in sorted(set(loc) & set(srv)) if loc[r] != srv[r]]
only_loc = sorted(set(loc) - set(srv))
only_srv = sorted(set(srv) - set(loc))

print(f"local files: {len(loc)} | server files: {len(srv)} | common: {len(set(loc)&set(srv))}")
print(f"\n=== 内容不同 ({len(common_diff)}) ===")
for r in common_diff:
    print(" ", r)
print(f"\n=== 仅本地 ({len(only_loc)}) ===")
for r in only_loc:
    print(" ", r)
print(f"\n=== 仅服务器 ({len(only_srv)}) ===")
for r in only_srv:
    print(" ", r)
