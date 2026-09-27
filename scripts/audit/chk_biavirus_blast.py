#!/usr/bin/env python3
"""决定性比对：13 个 Biavirus contig 到底像谁
参考 A  HG999358  Biavirus raunefjordenense = Prymnesium kappa virus (1.4 Mb, dsDNA, 原生生物)
参考 B  Guapo partitivirus (NCBI, unclassified Partitiviridae, 非 ICTV)
三档灵敏度：megablast(默认) / blastn(11) / dc-megablast，全部关掉 dust 屏蔽
另测：参考 B 与参考 A 是否互相重叠（若是，两套参考指向同一序列）
"""
import csv
import json
import subprocess
from collections import Counter
from pathlib import Path

DT = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
WORK = Path("/tmp/bia_check")
WORK.mkdir(exist_ok=True)
BLASTN = "/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin/blastn"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def clean(v):
    return (v or "").strip().strip('"')


def sh(cmd):
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return p.stdout, p.stderr


def get(url):
    o, e = sh('curl -s -g "%s"' % url)
    if not o.strip():
        print("    [curl 空响应] err=%s" % e.strip()[:120])
    return o


def load(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fasta_seqs(text, tag):
    """key = tag|第一个 token（整行含空格会让 blast 解析出错）"""
    out = {}
    name = None
    buf = []
    for ln in text.split("\n"):
        if ln.startswith(">"):
            if name:
                out["%s|%s" % (tag, name)] = "".join(buf)
            toks = ln[1:].split()
            name, buf = (toks[0] if toks else "%s_%d" % (tag, len(out))), []
        elif ln.strip():
            buf.append(ln.strip())
    if name:
        out["%s|%s" % (tag, name)] = "".join(buf)
    return out


# ---------- 1. 目标 contig ----------
plant_ids = set()
with open(DT / "09_Virome_Analysis/All_plant.viruses.fasta", encoding="utf-8", errors="replace") as f:
    for ln in f:
        if ln.startswith(">"):
            plant_ids.add(ln[1:].split()[0].strip())
base = {r["contig_id"]: r for r in load(
    DT / "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")}
targets = [c for c, r in base.items() if c in plant_ids and clean(r.get("Genus")) == "Biavirus"]
allq = fasta_seqs(open(DT / "09_Virome_Analysis/All_plant.viruses.fasta", encoding="utf-8",
                       errors="replace").read(), "q")
qs = {k.split("|", 1)[1]: s for k, s in allq.items() if k.split("|", 1)[1] in targets}
print("目标 contig: %d 个\n" % len(qs))

# ---------- 2. 参考 A ----------
hq = get("%s/efetch.fcgi?db=nuccore&id=HG999358&rettype=fasta&retmode=text" % EUTILS)
refA = fasta_seqs(hq, "A_HG999358_Prymnesium_kappa_virus")
for k, v in refA.items():
    print("参考A: %s  长度 %d" % (k, len(v)))

# ---------- 3. 参考 B ----------
ids, titles = [], {}
for term in ("Guapo%20partitivirus", "txid2732505", "txid3078405"):
    js = get("%s/esearch.fcgi?db=nuccore&term=%s&retmax=30&retmode=json" % (EUTILS, term))
    try:
        d = json.loads(js)["esearchresult"]
        print("  esearch %-22s -> %s 条 %s" % (term, d.get("count"), d.get("idlist", [])[:8]))
        ids += d.get("idlist", [])
    except Exception as ex:
        print("  esearch %-22s -> 失败 %s | 响应头 %s" % (term, str(ex)[:50], js[:90].replace("\n", " ")))
ids = list(dict.fromkeys(ids))[:30]
blob = ""
if ids:
    sm = get("%s/esummary.fcgi?db=nuccore&id=%s&retmode=json" % (EUTILS, ",".join(ids)))
    try:
        for k, v in json.loads(sm)["result"].items():
            if k != "uids":
                titles[k] = "%s | %s | %s bp" % (v.get("accessionversion"), v.get("title", "")[:60], v.get("slen"))
                print("     %s" % titles[k])
    except Exception as ex:
        print("     esummary 失败 %s" % str(ex)[:60])
    blob = get("%s/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=text" % (EUTILS, ",".join(ids)))
refB = fasta_seqs(blob, "B_guapo")
print("参考B: %d 条" % len(refB))
for k, v in list(refB.items())[:6]:
    print("     %s 长度 %d" % (k, len(v)))
print()

# ---------- 4. 建库 + 三档比对 ----------
allref = {}
allref.update(refA)
allref.update(refB)
rf = WORK / "refs2.fa"
with open(rf, "w") as f:
    for k, v in allref.items():
        f.write(">%s\n%s\n" % (k, v))
sh("makeblastdb -in %s -dbtype nucl -out %s/db2 2>&1" % (rf, WORK))
qf = WORK / "q2.fa"
with open(qf, "w") as f:
    for k, v in qs.items():
        f.write(">%s\n%s\n" % (k, v))

TASKS = [("megablast", ""), ("blastn", "-task blastn -word_size 11"), ("dc-megablast", "-task dc-megablast")]
print("=" * 130)
print("contig -> 参考 (三档灵敏度, dust 关闭)")
print("=" * 130)
summary = {}
for tname, opt in TASKS:
    out, _ = sh("%s -query %s -db %s/db2 -outfmt '6 qseqid sseqid pident length qlen qcovs evalue bitscore' "
                "-max_target_seqs 3 -evalue 1e-3 -dust no %s" % (BLASTN, qf, WORK, opt))
    rows = [r.split("\t") for r in out.strip().split("\n") if r.strip()]
    print("\n--- %s (%d 条命中) ---" % (tname, len(rows)))
    for q, s, pid, al, ql, qcov, ev, bs in rows:
        tgt = "B" if s.startswith("B_") else "A"
        print("  %-50s -> %s %-40s pid=%-6s alen=%-6s qcov=%-4s e=%s" %
              (q[:50], tgt, s[2:44], pid, al, qcov, ev))
        summary.setdefault(q, []).append((tname, tgt, pid, qcov))
print()
for q in qs:
    got = summary.get(q, [])
    print("  %-50s %s" % (q[:50], got if got else "三档全无命中"))

# ---------- 5. 参考 B 对参考 A（是否同一序列） ----------
print("\n" + "=" * 130)
print("参考B(Guapo) 对 参考A(HG999358) 是否重叠")
print("=" * 130)
if refB:
    bf = WORK / "refB.fa"
    with open(bf, "w") as f:
        for k, v in refB.items():
            f.write(">%s\n%s\n" % (k, v))
    out, _ = sh("%s -query %s -db %s/db2 -outfmt '6 qseqid sseqid pident length qlen qcovs evalue' "
                "-max_target_seqs 3 -evalue 1e-3 -dust no -task blastn -word_size 11" % (BLASTN, bf, WORK))
    print(out.strip() or "  无重叠命中")
else:
    print("  参考B 未取到，跳过")

# ---------- 6. ACVirus 那 768 个 contig 的家族画像 ----------
print("\n" + "=" * 130)
print("ACVirus 判成 Biavirus/Schizomimiviridae 的 contig，长度分布与共识科")
print("=" * 130)
ac = load(DT / "05_Taxonomy/Votus.integrated/standardized_ACVirus.tsv")
b_ids = [r["contig_id"] for r in ac if "Biavirus" in (r.get("Genus") or "")
         or "Schizomimiviridae" in (r.get("Family") or "")]
con = Counter()
for r in load(DT / "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"):
    if r["contig_id"] in set(b_ids):
        con[clean(r.get("Family")) or "NA"] += 1
print("  contig 数: %d" % len(b_ids))
print("  共识科分布 top10: %s" % dict(con.most_common(10)))
lens = sorted(len(v) for v in qs.values())
print("  植物子集这 13 个的长度: %s" % lens)
