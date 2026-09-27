#!/usr/bin/env python3
"""ACVirus 的 Biavirus 判定：直接复现它的蛋白比对，看命中谁"""
import csv
import subprocess
from collections import Counter
from pathlib import Path

DB = Path("/home/zhangwenda/database/virus-db/acvirus_db")
DT = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
STD = DT / "05_Taxonomy/Votus.integrated"
W = Path("/tmp/bia_check")
W.mkdir(exist_ok=True)


def sh(c):
    p = subprocess.run(c, shell=True, capture_output=True, text=True)
    return (p.stdout or "") + (p.stderr or "")


def clean(v):
    return (v or "").strip().strip('"')


def load(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fa(path, keep=None):
    out, name, buf = {}, None, []
    for ln in open(path, encoding="utf-8", errors="replace"):
        if ln.startswith(">"):
            if name:
                out[name] = "".join(buf)
            toks = ln[1:].split()
            name, buf = (toks[0] if toks else str(len(out))), []
        else:
            buf.append(ln.strip())
    if name:
        out[name] = "".join(buf)
    if keep:
        out = {k: v for k, v in out.items() if k in keep}
    return out


print("=" * 118)
print("1. ACVirus 本体在哪")
print("=" * 118)
print(sh("ls -d /home/zhangwenda/mambaforge/envs/*/bin/ACVirus /home/zhangwenda/*/ACVirus /opt/*/ACVirus 2>/dev/null | head; "
         "find /home/zhangwenda -maxdepth 4 -name 'ACVirus*' -not -path '*/database/*' 2>/dev/null | head -10"))

print("=" * 118)
print("2. 库里有没有 Biavirus / HG999358 / OQ968299 / Guapo")
print("=" * 118)
for f in ["vmr.tsv", "taxa.txt", "species.tsv", "processed_accessions_b.tsv",
          "processed_accessions_b.fa_names.tsv", "fixed_vmr_b.tsv", "bad_accessions_b.tsv"]:
    p = DB / f
    if not p.exists():
        continue
    for pat in ["Biavirus", "HG999358", "OQ968299", "Guapo", "Schizomimiviridae"]:
        n = sh("grep -c -e '%s' '%s' 2>/dev/null" % (pat, p)).strip().split("\n")[-1]
        if n and n != "0":
            print("  %-38s %-20s %s 行" % (f, pat, n))
print("  --- 命中的样例 ---")
for pat in ["Biavirus", "HG999358", "OQ968299", "Guapo"]:
    r = sh("grep -h -m2 -e '%s' %s/vmr.tsv %s/taxa.txt %s/species.tsv 2>/dev/null | head -3" % (pat, DB, DB, DB))
    if r.strip():
        print("  [%s]\n%s" % (pat, r.strip()[:600]))
print("  --- 蛋白库里搜（266MB, 只取计数与首条）---")
for pat in ["Biavirus", "HG999358", "OQ968299", "Guapo", "Prymnesium"]:
    c = sh("grep -c -e '%s' %s/all_virus.faa 2>/dev/null" % (pat, DB)).strip().split("\n")[-1]
    s = sh("grep -m1 -e '%s' %s/all_virus.faa 2>/dev/null" % (pat, DB)).strip()
    print("  %-14s 头部出现 %s 次   首条: %s" % (pat, c, s[:110]))
print("  --- 核酸库 ---")
for pat in ["HG999358", "OQ968299", "Guapo"]:
    c = sh("grep -c -e '%s' %s/all_virus.fna 2>/dev/null" % (pat, DB)).strip().split("\n")[-1]
    s = sh("grep -m1 -e '%s' %s/all_virus.fna 2>/dev/null" % (pat, DB)).strip()
    print("  %-14s 出现 %s 次   首条: %s" % (pat, c, s[:110]))

print("=" * 118)
print("3. 拿 13 个 contig 直接打 ACVirus 的 database.dmnd，看命中什么蛋白")
print("=" * 118)
plant = set()
for ln in open(DT / "09_Virome_Analysis/All_plant.viruses.fasta", encoding="utf-8", errors="replace"):
    if ln.startswith(">"):
        plant.add(ln[1:].split()[0])
final = {r["contig_id"]: r for r in load(STD / "final_integrated_classification.tsv")}
targets = [c for c, r in final.items() if c in plant and clean(r.get("Genus")) == "Biavirus"]
seqs = fa(DT / "09_Virome_Analysis/All_plant.viruses.fasta", set(targets))
qf = W / "q13.fa"
with open(qf, "w") as f:
    for k, v in seqs.items():
        f.write(">%s\n%s\n" % (k, v))
print("  查询: %d 条" % len(seqs))
DIAMOND = sh("which diamond 2>/dev/null").strip() or "/home/zhangwenda/biosoft/diamond/bin/bin/diamond"
print("  diamond: %s" % DIAMOND)
out = sh("%s blastx -q %s -d %s/database.dmnd -o %s/dia.tsv -f 6 qseqid sseqid pident length "
         "evalue bitscore stitle --max-target-seqs 10 --threads 16 2>&1 | tail -3" % (DIAMOND, qf, DB, W))
print("  运行输出: %s" % out.strip()[:300])
rows = [l.split("\t") for l in open(W / "dia.tsv") if l.strip()]
print("  命中 %d 条\n" % len(rows))
seen = set()
for q, s, pid, al, ev, bs, *st in sorted(rows, key=lambda r: (r[0], -float(r[5]))):
    title = st[0][:80] if st else ""
    if q not in seen:
        print("  %s" % q)
        seen.add(q)
    print("      bits=%-6s pid=%-6s alen=%-5s e=%-9s %-22s %s" % (bs, pid, al, ev, s[:22], title))
