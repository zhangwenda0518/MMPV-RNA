#!/usr/bin/env python3
"""A. ACVirus 的选科规则  B. 13 行全工具票型  C. 命中基因的注释"""
import csv
import re
import subprocess
from pathlib import Path

DB = Path("/home/zhangwenda/database/virus-db/acvirus_db")
DT = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
STD = DT / "05_Taxonomy/Votus.integrated"
SRC = [Path("/home/zhangwenda/biosoft/ACVirus"), Path("/home/zhangwenda/MMPV-RNA/biosoft/ACVirus")]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def sh(c):
    p = subprocess.run(c, shell=True, capture_output=True, text=True)
    return (p.stdout or "") + (p.stderr or "")


def clean(v):
    return (v or "").strip().strip('"')


def load(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


print("=" * 120)
print("A1. ACVirus 源码位置与文件")
print("=" * 120)
for s in SRC:
    print("\n[%s]" % s)
    print(sh("ls -la %s 2>&1 | head -25" % s))
print("\n所有 .py 文件:")
print(sh("find %s -name '*.py' -printf '%%s\\t%%p\\n' 2>/dev/null | sort -rn | head -20"
         % " ".join(str(s) for s in SRC)))

print("=" * 120)
print("A2. 选科规则：搜 min_coverage / lineage 相关代码")
print("=" * 120)
for s in SRC:
    r = sh("grep -rn -e 'taxon_min_coverage' -e 'min_coverage' -e 'coverage_threshold' "
           "-e 'def classify' -e 'lineage' --include='*.py' %s 2>/dev/null | head -40" % s)
    if r.strip():
        print("\n[%s]\n%s" % (s, r[:4000]))

print("=" * 120)
print("A3. taxon_min_coverage.csv 里 Biavirus / Schizomimiviridae / HG999358")
print("=" * 120)
print("表头:", sh("head -1 %s/taxon_min_coverage.csv" % DB).strip()[:400])
print("前 3 行:", sh("head -4 %s/taxon_min_coverage.csv | tail -3" % DB).strip()[:400])
for pat in ["Biavirus", "Schizomimiviridae", "HG999358", "Prymnesium", "Tupanvirus", "Lymphocystis"]:
    r = sh("grep -m2 -e '%s' %s/taxon_min_coverage.csv 2>/dev/null" % (pat, DB)).strip()
    print("  [%-20s] %s" % (pat, r[:200] if r else "(无)"))

print("=" * 120)
print("A4. 命中 subject id 怎么映射到 lineage（看 sseqid 解析）")
print("=" * 120)
for s in SRC:
    r = sh("grep -rn -e 'sseqid' -e 'qseqid' -e 'split(' --include='*.py' %s 2>/dev/null | head -30" % s)
    if r.strip():
        print("\n[%s]\n%s" % (s, r[:3000]))

print("=" * 120)
print("B. 13 行全工具票型")
print("=" * 120)
plant = set()
for ln in open(DT / "09_Virome_Analysis/All_plant.viruses.fasta", encoding="utf-8", errors="replace"):
    if ln.startswith(">"):
        plant.add(ln[1:].split()[0])
final = {r["contig_id"]: r for r in load(STD / "final_integrated_classification.tsv")}
tools = ["ACVirus", "CAT", "VITAP", "diamond_lca", "genomad", "metabuli", "mmseqs"]
tabs = {t: {r["contig_id"]: r for r in load(STD / ("standardized_%s.tsv" % t))} for t in tools}
rawmap = {}
rp = DT / "05_Taxonomy/Votus.classed/Votus_ACVirus_taxonomy.tsv"
if rp.exists():
    for i, ln in enumerate(open(rp, encoding="utf-8", errors="replace")):
        if i:
            ps = ln.rstrip("\n").split("\t")
            if len(ps) >= 3:
                rawmap[ps[0]] = ps[2]
targets = sorted(c for c in final if c in plant and clean(final[c].get("Genus")) == "Biavirus")
for c in targets:
    print("\n%s" % c)
    print("    ACVirus 原始 lineage: %s" % rawmap.get(c, "<无>")[:160])
    for t in tools:
        r = tabs[t].get(c)
        print("    %-12s %-28s %s" % (t, clean(r.get("Family")) if r else "-", clean(r.get("Genus")) if r else "-"))
    print("    共识: Family=%s Genus=%s agree=%s/%s" %
          (clean(final[c].get("Family")), clean(final[c].get("Genus")),
           clean(final[c].get("Family_agree")), clean(final[c].get("Genus_agree"))))

print("=" * 120)
print("C. 命中基因的注释（拿 KX643370 LCDV-3 的 feature table 试）")
print("=" * 120)
for acc in ["KX643370", "HG999358"]:
    r = sh('curl -s -g "%s/efetch.fcgi?db=nuccore&id=%s&rettype=ft&retmode=text" | head -c 200' % (EUTILS, acc))
    print("  [%s] ft 前 200 字节: %s" % (acc, r[:200]))
r = sh('curl -s -g "%s/efetch.fcgi?db=nuccore&id=KX643370&rettype=gb&retmode=text" > /tmp/bia_check/kx.gb; wc -c /tmp/bia_check/kx.gb' % EUTILS)
print("  gb 下载: %s" % r.strip())
gb = Path("/tmp/bia_check/kx.gb")
if gb.exists() and gb.stat().st_size > 10000:
    txt = gb.read_text(errors="replace")
    # 找覆盖 154904-156967 的 CDS
    for m in re.finditer(r'( {5}CDS {2,}complement\((\d+)\.\.(\d+)\)| {5}CDS {2,}(\d+)\.\.(\d+))(.{0,400}?)(?=\n {5}\S)', txt, re.S):
        s = int(m.group(2) or m.group(4) or 0)
        e = int(m.group(3) or m.group(5) or 0)
        if abs(s - 154904) < 3000 or abs(e - 156967) < 3000:
            prod = re.search(r'/product="([^"]+)"', m.group(6) or "")
            gene = re.search(r'/gene="([^"]+)"', m.group(6) or "")
            note = re.search(r'/note="([^"]+)"', m.group(6) or "")
            print("  KX643370 %d-%d product=%s gene=%s note=%s" %
                  (s, e, prod.group(1) if prod else "-", gene.group(1) if gene else "-",
                   (note.group(1)[:80] if note else "-")))
