"""层2 步骤C: 清理旧树 + 按新 classify 结果重拆 KEEP + 重建受影响科的树。

逻辑对齐 virome_pipeline.py::run_analysis_verify 的建树段落:
  1. 从 acvirus_classify/final_result_with_confidence.tsv 读 Family / Nucleotide / Genus
  2. 用 genus_lens 做 len_ratio>=0.7 过滤
  3. 把 class_KEEP.fasta 拆成 acvirus_trees/<Fam>/keep_<Fam>.fasta
  4. 对受影响的 8 个科 (有树目录的) 重建树
"""
import csv
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO

HOME = os.path.expanduser("~")
ACV_SCRIPT = f"{HOME}/MMPV-RNA/virome_discovery_pipeline/utils/acvirus_tree_pro.py"
DB = f"{HOME}/database/virus-db/acvirus_db"
DB_TAXA = f"{DB}/taxa.txt"
DB_FASTA = f"{DB}/all_virus.fasta"
GENUS_LENS = f"{HOME}/database/virus-db/db/genus_lens"

ACVDIR = Path("/home/zhangwenda/data-test/out/09b_Analysis_Verify")
CLS = ACVDIR / "acvirus_classify"
TREE_ROOT = ACVDIR / "acvirus_trees"
KEEP_FA = ACVDIR / "class_KEEP.fasta"

# 受影响科 (有树目录的 8 个)
AFFECTED = ["Tombusviridae", "Partitiviridae", "Aspiviridae", "Caulimoviridae",
            "Metaviridae", "Geminiviridae", "Potyviridae", "Secoviridae"]
MIN_RATIO = 0.7
FAM_THREADS = 15
N_PARALLEL = 4

cls_res = CLS / "final_result_with_confidence.tsv"
print(f"classify 结果: {cls_res}")

# ── 1. 读 classify: contig -> Family / Genus ──
fam_of, genus_of = {}, {}
fams = set()
with open(cls_res, errors="replace") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        fam = (row.get("Family", "") or "").strip()
        cid = (row.get("Nucleotide", "") or row.get("Contig", "") or "").strip()
        gen = (row.get("Genus", "") or "").strip()
        if fam and fam != "NA":
            fams.add(fam)
            if cid:
                fam_of[cid] = fam
        if cid and gen:
            genus_of[cid] = gen
print(f"新 classify 科数: {len(fams)}")
print(f"受影响科在新结果中: {[f for f in AFFECTED if f in fams]}")

# ── 2. 属平均长度 ──
gl = {}
if os.path.exists(GENUS_LENS):
    for line in open(GENUS_LENS):
        p = line.rstrip().split("\t")
        if len(p) >= 2 and p[0].startswith("g__"):
            try:
                gl[p[0][3:]] = float(p[1])
            except Exception:
                pass
print(f"genus_lens 条目: {len(gl)}")

# ── 3. 拆 KEEP (仅受影响科, 其余科不动) ──
print(f"\n=== 重拆受影响科的 keep_<Fam>.fasta ===")
fam_recs = {}
n_drop = 0
n_total = 0
for rec in SeqIO.parse(str(KEEP_FA), "fasta"):
    n_total += 1
    f_ = fam_of.get(rec.id) or fam_of.get(rec.id.split()[0])
    if not f_ or f_ not in AFFECTED:
        continue
    g = genus_of.get(rec.id, "")
    glen = gl.get(g, 0.0) or 0.0
    if glen > 0:
        ratio = len(rec.seq) / glen
        if ratio < MIN_RATIO:
            n_drop += 1
            continue
    fam_recs.setdefault(f_, []).append(rec)

print(f"class_KEEP 总数: {n_total}")
print(f"受影响科保留: {sum(len(v) for v in fam_recs.values())} (len_ratio<{MIN_RATIO} 筛掉 {n_drop})")
for f_, recs in sorted(fam_recs.items(), key=lambda x: -len(x[1])):
    print(f"  {f_:<22} {len(recs)}")

# 写 keep_<Fam>.fasta (覆盖)
for f_ in AFFECTED:
    d = TREE_ROOT / f_
    d.mkdir(parents=True, exist_ok=True)
    recs = fam_recs.get(f_, [])
    with open(d / f"keep_{f_}.fasta", "w") as fo:
        for rc in recs:
            SeqIO.write(rc, fo, "fasta")

# ── 4. 清理旧树产物 (保留新写的 keep_<Fam>.fasta) ──
print(f"\n=== 清理旧树产物 ===")
for f_ in AFFECTED:
    d = TREE_ROOT / f_
    removed = []
    for child in d.iterdir():
        if child.name == f"keep_{f_}.fasta":
            continue
        if child.is_dir():
            shutil.rmtree(child)
            removed.append(child.name + "/")
        else:
            child.unlink()
            removed.append(child.name)
    print(f"  {f_:<22} 清理 {len(removed)} 项")

# ── 5. 并行建树 ──
print(f"\n=== 并行建树 ({N_PARALLEL} 路 x {FAM_THREADS} threads) ===")


def build_one(fam):
    fam_out = TREE_ROOT / fam
    fam_keep = fam_out / f"keep_{fam}.fasta"
    if not fam_keep.exists():
        return fam, 3, "no keep file"
    n = sum(1 for l in open(fam_keep) if l.startswith(">"))
    if n == 0:
        return fam, 3, "empty keep"
    if list(fam_out.glob("*.treefile")):
        return fam, 4, "tree exists"
    cmd = ["python3", ACV_SCRIPT,
           "--mode", "macro", "--target_name", fam, "--target_rank", "Family",
           "--contigs", str(fam_keep), "--db_taxa", DB_TAXA,
           "--db_fasta", DB_FASTA, "--outdir", str(fam_out),
           "--threads", str(FAM_THREADS)]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
        return fam, r.returncode, f"{time.time()-t0:.0f}s {(r.stderr or '')[-150:]}"
    except Exception as e:
        return fam, 99, str(e)[-200:]


todo = [f for f in AFFECTED if (TREE_ROOT / f / f"keep_{f}.fasta").exists()]
print(f"待建树: {len(todo)} 科 {todo}")

results = []
with ThreadPoolExecutor(max_workers=N_PARALLEL) as ex:
    futs = {ex.submit(build_one, f): f for f in todo}
    for fu in as_completed(futs):
        fam, rc, msg = fu.result()
        results.append((fam, rc, msg))
        print(f"  [{fam}] rc={rc}  {msg}", flush=True)

print("\n=== 建树汇总 ===")
ok = sum(1 for _, rc, _ in results if rc == 0)
for fam, rc, msg in sorted(results):
    mark = "OK" if rc == 0 else ("SKIP" if rc in (3, 4) else "FAIL")
    print(f"  {mark:<5} {fam:<22} rc={rc}  {msg[:80]}")
print(f"\n成功: {ok}/{len(results)}")
