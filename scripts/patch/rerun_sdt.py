#!/usr/bin/env python3
"""单独补跑 9b 的 SDT 属级矩阵 (P5 修复后)。
复用管线同样的筛选逻辑: classify 的 Genus 列 + keep_summary len_ratio。"""
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path("/home/zhangwenda/data-test/out")
ACV = ROOT / "09b_Analysis_Verify"
SDT_SCRIPT = Path("/home/zhangwenda/MMPV-RNA/virome_analysis_pipeline/sdt_genus_matrix.py")
CLS = ACV / "acvirus_classify" / "final_result_with_confidence.tsv"
KEEP_FA = ACV / "class_KEEP.fasta"
SDT_OUT = ACV / "SDT_matrix"
MIN_RATIO = 0.0
THREADS = 8
MAX_SEQS = 80  # 超过则跳过 (O(N^2) 组合控制, 与管线 --sdt-max-seqs 默认一致)

keep_ids = set()
with open(KEEP_FA) as f:
    for line in f:
        if line.startswith(">"):
            keep_ids.add(line[1:].strip())

ratio_map = {}
ksum = ACV / "keep_summary.tsv"
if ksum.is_file():
    for r in csv.DictReader(open(ksum), delimiter="\t"):
        try:
            ratio_map[r["contig_id"]] = float(r["len_ratio"] or 0)
        except Exception:
            pass

gen_of = {}
with open(CLS) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        cid = (row.get("Nucleotide", "") or "").strip()
        g = (row.get("Genus", "") or "").strip()
        if cid in keep_ids and g and g != "NA":
            gen_of.setdefault(g, []).append(cid)

todo = []
skip_big = []
for g, cids in sorted(gen_of.items(), key=lambda x: -len(x[1])):
    user_ok = [c for c in cids if ratio_map.get(c, 1.0) >= MIN_RATIO]
    if len(user_ok) < 2:
        continue
    if len(user_ok) > MAX_SEQS:
        skip_big.append((g, len(user_ok)))
        continue
    od = SDT_OUT / g
    if (od / f"SDT_{g}.pdf").exists() or (od / f"SDT_{g}.png").exists():
        continue
    todo.append((g, len(user_ok)))

print(f"待补跑属数: {len(todo)} / 总属数 {len(gen_of)}  (超上限 {MAX_SEQS} 跳过 {len(skip_big)} 个)", flush=True)
for g, n in skip_big:
    print(f"  [skip] {g}: {n} 条 → {n*(n-1)//2} 对", flush=True)
print(f"总 pair 数: {sum(n*(n-1)//2 for _, n in todo)}", flush=True)
ok = fail = 0
for i, (g, n) in enumerate(todo, 1):
    od = SDT_OUT / g
    od.mkdir(parents=True, exist_ok=True)
    try:
        r = subprocess.run(
            ["python3", str(SDT_SCRIPT), "--analysis_dir", str(ACV),
             "--genus", g, "-o", str(od), "--threads", str(THREADS),
             "--short_labels", "--palette", "cividis", "--resume",
             "--other_count", "0", "--target_cap", "20"],
            capture_output=True, text=True, timeout=7200)
        if r.returncode == 0:
            ok += 1
            print(f"[{i}/{len(todo)}] {g} ({n}) OK", flush=True)
        else:
            fail += 1
            print(f"[{i}/{len(todo)}] {g} ({n}) FAIL rc={r.returncode} {r.stderr[-150:]}", flush=True)
    except Exception as e:
        fail += 1
        print(f"[{i}/{len(todo)}] {g} ({n}) ERR {e}", flush=True)

print(f"\n完成: ok={ok} fail={fail}", flush=True)
