#!/usr/bin/env python3
"""
panvirome_novel_known.py — 新病毒 vs 已知病毒综合判定表生成器
================================================================
对每个宿主物种的 rescue centroids (vOTU)，综合三层证据 + 流程溯源信息，
判定每条序列是"已知病毒"还是"新病毒候选"，并附上判定依据。

三层证据:
  1. 核酸 identity  — blastn 对 nt 病毒库
  2. 蛋白 identity  — mmseqs translated-search 对 ICTV 蛋白库
  3. CDD 结构域     — mmseqs profile-search 对 CDD 库 (区分 RdRp/RT vs 宿主结构域)

流程溯源:
  - blast_*      02 鉴定阶段 diamond 对 RVDB 病毒库的 best-hit (pident/alnlen/evalue/target)
  - identify_tools 02 鉴定阶段哪些工具判成病毒 (blast/genomad/metabuli/rdrpcatch/viralverify/virbot)
  - classify_tool  05 分类阶段 primary_tool (ACVirus/CAT/metabuli/...)

输入(自动拼接):
  {base}/{sp_prefix}_{sp}_out/08_Rescue/Plant/centroids/final_centroids.fasta  (rescue vOTU)
  {base}/{sp_prefix}_{sp}_out/02a_Identification/                              (02 鉴定结果)
  {base}/{sp_prefix}_{sp}_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv (05 分类)

中间产物(复用, 存在则跳过重跑):
  {base}/all_rescue_centroids.fasta
  {base}/blastn_ntviruses_best_hit.tsv
  {base}/mmseqs_blastx_best_hit.tsv
  {base}/cdd_hits.tsv

用法:
  python panvirome_novel_known.py --base ~/virus/data-2026/data-test \
      --species barbarum,ruthenicum,chinense,amarum --detect-tag bowtie2
"""

import argparse
import csv
import os
import sys
import re
import subprocess
from collections import defaultdict
from pathlib import Path

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))


VIRAL_DOM = ["RdRp", "RT_", "RT_LTR", "RNase_HI_RT", "coat", "capsid", "viridae",
             "replicase", "helicase", "nsP", "movement", "Polyprotein"]
HOST_DOM = ["tubulin", "KOG", "PLN", "PTZ", "RING", "Ras", "kinase", "GEF",
            "Glycogen", "phosphorylase", "peptidase", "WD40", "A2M", "RecA", "ALDH"]

# 默认数据库/工具路径
DB = {
    "nt_viruses": "~/database/nt-db/nt_viruses/nt_viruses",
    "ictv_protein": "~/database/virus-db/ictv-db/ictv_nr_db/ictv_nr_db",
    "cdd_db": "~/database/cdd/cdd-db/cdd_db",
    "cdd_annot": "~/database/cdd/cddannot.dat",
}
TOOLS = {
    "blastn": "~/biosoft/ncbi-blast-2.13.0+/bin/blastn",
    "mmseqs": "~/biosoft/binary/mmseqs",
    "prodigal": "~/mambaforge/bin/prodigal",
}


def parse_args():
    p = argparse.ArgumentParser(
        description="新病毒 vs 已知病毒综合判定表生成器",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--base", required=True, help="项目根目录")
    p.add_argument("--species", default="barbarum,ruthenicum,chinense,amarum",
                   help="宿主物种列表, 逗号分隔")
    p.add_argument("--sp-prefix", default="RNA-Lycium",
                   help="discovery 输出目录前缀, 即 {prefix}_{sp}_out/")
    p.add_argument("--out", default=None, help="输出 TSV (默认 {base}/novel_vs_known_final.tsv)")
    p.add_argument("--threads", type=int, default=60, help="比对线程数")
    p.add_argument("--rerun-blastn", action="store_true", help="强制重跑 blastn")
    p.add_argument("--rerun-mmseqs", action="store_true", help="强制重跑 mmseqs 蛋白比对")
    p.add_argument("--rerun-cdd", action="store_true", help="强制重跑 CDD 注释")
    return p.parse_args()


def sh(cmd):
    print(f"  [RUN] {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def resolve(p):
    return str(Path(p).expanduser())


# ── 1. 合并 rescue centroids ──────────────────────────────
def merge_centroids(base, species, sp_prefix):
    all_fa = Path(base) / "all_rescue_centroids.fasta"
    if all_fa.exists():
        print(f"[1] 复用 centroids: {all_fa}")
        return all_fa
    parts = []
    for sp in species:
        fa = Path(base) / f"{sp_prefix}_{sp}_out" / _D['d_rescue'] / "Plant" / "centroids" / "final_centroids.fasta"
        if fa.exists():
            parts.append(str(fa))
    if not parts:
        raise SystemExit("[ERROR] 未找到任何 final_centroids.fasta")
    sh(f"cat {' '.join(parts)} > {all_fa}")
    print(f"[1] 合并 centroids: {all_fa} ({sum(1 for _ in open(all_fa) if _.startswith('>'))} 条)")
    return all_fa


# ── 2. blastn 对 nt 病毒库 (核酸 identity) ────────────────
def run_blastn(base, all_fa, threads, rerun):
    out = Path(base) / "blastn_ntviruses_best_hit.tsv"
    raw = Path(base) / "blastn_vs_ntviruses.tsv"
    if out.exists() and not rerun:
        print("[2] 复用 blastn 结果")
        return out
    nt = resolve(DB["nt_viruses"])
    sh(f"{resolve(TOOLS['blastn'])} -query {all_fa} -db {nt} "
        f"-outfmt '6 qseqid sseqid pident length qcovs stitle' "
        f"-max_target_seqs 3 -num_threads {threads} -evalue 1e-5 -out {raw}")
    best = {}
    with open(raw, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 6 and p[0] not in best:
                best[p[0]] = {"pident": p[2], "qcovs": p[4], "stitle": p[5]}
    all_ids = {l[1:].split()[0] for l in open(all_fa) if l.startswith(">")}
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["contig_id", "pident", "qcovs", "stitle"])
        for c in all_ids:
            b = best.get(c)
            w.writerow([c, b["pident"] if b else "", b["qcovs"] if b else "",
                        b["stitle"] if b else "No hit"])
    print(f"[2] blastn 完成: {out}")
    return out


# ── 3. mmseqs 蛋白比对 (蛋白 identity) ────────────────────
def run_mmseqs(base, all_fa, threads, rerun):
    out = Path(base) / "mmseqs_blastx_best_hit.tsv"
    raw = Path(base) / "mmseqs_blastx_ictv.tsv"
    if out.exists() and not rerun:
        print("[3] 复用 mmseqs 结果")
        return out
    prot = resolve(DB["ictv_protein"])
    mm = resolve(TOOLS["mmseqs"])
    tmp = Path(base) / "tmp_mmseqs"
    sh(f"{mm} easy-search {all_fa} {prot} {raw} {tmp} --search-type 2 -e 1e-3 "
        f"--max-seqs 3 --threads {threads} --format-output query,target,pident,qcov,tcov,evalue")
    best = {}
    with open(raw, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 6 and p[0] not in best:
                best[p[0]] = {"pident": p[2], "target": p[1]}
    # target accession -> 物种名 (用 ICTV accession2taxid + names.dmp)
    acc2tax = {}
    acc2tax_f = Path(resolve(DB["ictv_protein"])).parent.parent / "accession2taxid.tsv"
    with open(acc2tax_f, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] in {b["target"] for b in best.values()}:
                acc2tax[p[0]] = p[1]
    tax2name = {}
    names_f = Path(resolve(DB["ictv_protein"])).parent.parent / "ncbi_taxdump" / "names.dmp"
    with open(names_f, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("|")
            if len(p) >= 4 and p[3].strip() == "scientific name":
                tax2name[p[0].strip()] = p[1].strip()
    all_ids = {l[1:].split()[0] for l in open(all_fa) if l.startswith(">")}
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["contig_id", "pident", "target", "species"])
        for c in all_ids:
            b = best.get(c)
            sp_name = tax2name.get(acc2tax.get(b["target"], ""), b["target"]) if b else "No hit"
            w.writerow([c, b["pident"] if b else "", b["target"] if b else "", sp_name])
    print(f"[3] mmseqs 完成: {out}")
    return out


# ── 4. CDD 注释 ─────────────────────────────────────────
def run_cdd(base, all_fa, threads, rerun):
    out = Path(base) / "cdd_hits.tsv"
    if out.exists() and not rerun:
        print("[4] 复用 CDD 结果")
        return out
    orfs = Path(base) / "all_rescue_orfs.faa"
    prod = resolve(TOOLS["prodigal"])
    mm = resolve(TOOLS["mmseqs"])
    cdd_db = resolve(DB["cdd_db"])
    sh(f"{prod} -i {all_fa} -a {orfs} -p meta -q -o /dev/null")
    qdb = Path(base) / "queryDB_cdd"
    rdb = Path(base) / "resultDB_cdd"
    tmp = Path(base) / "tmp_cdd"
    sh(f"{mm} createdb {orfs} {qdb} && {mm} search {qdb} {cdd_db} {rdb} {tmp} "
        f"--threads {threads} -s 7.5 -e 1e-3 && "
        f"{mm} convertalis {qdb} {cdd_db} {rdb} {out} --format-output query,target,evalue,bits,pident")
    print(f"[4] CDD 完成: {out}")
    return out


# ── 读取 TSV 工具 ──────────────────────────────────────
def load_map(path, key="contig_id"):
    m = {}
    if not os.path.exists(path):
        return m
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            m[r[key]] = r
    return m


def cdd_to_contig(cdd_hits, cdd_annot):
    cdd2name = {}
    with open(cdd_annot, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[1] not in cdd2name:
                cdd2name[p[1]] = p[2]
    contig_cdd = defaultdict(list)
    with open(cdd_hits, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            q = p[0]
            m = re.match(r"^(.*)_(\d+)$", q)
            contig = m.group(1) if m else q
            contig_cdd[contig].append(cdd2name.get(p[1], p[1]))
    return contig_cdd


# ── 判定 ───────────────────────────────────────────────
def has_viral(cdd):
    return any(k.lower() in cdd.lower() for k in VIRAL_DOM)


def has_host(cdd):
    return any(k.lower() in cdd.lower() for k in HOST_DOM)


def judge(nt, aa, cdd):
    nt_p = float(nt.get("pident")) if nt.get("pident") else None
    nt_q = float(nt.get("qcovs")) if nt.get("qcovs") else None
    aa_p = float(aa.get("pident")) if aa.get("pident") else None
    aa_sp = aa.get("species", "")
    viral = has_viral(cdd)
    host = has_host(cdd)

    if nt_p is not None and nt_q is not None and nt_q >= 50:
        if nt_p >= 85:
            return "已知病毒"
        return "新种候选"
    if aa_p is not None:
        if any(k in aa_sp for k in ["phage", "Caudoviricetes", "Siphoviridae", "Myoviridae"]):
            return "非病毒(噬菌体)"
        if aa_p >= 70:
            return "已知病毒(蛋白)"
        if aa_p >= 50:
            return "新种候选(蛋白)"
        if viral:
            return "新属/远缘候选"
        return "远缘未定"
    if viral:
        return "远缘新病毒候选"
    if host:
        return "非病毒(宿主)"
    return "未定"


# ── 02 阶段 blast + 鉴定工具 ────────────────────────────
def extract_02(base, species, sp_prefix, contigs):
    sample_of = {}
    for c in contigs:
        m = re.match(r"^([A-Z0-9]+)_clean", c)
        sample_of[c] = m.group(1) if m else c.split("_")[0]
    sample2dir = {}
    for sp in species:
        d = Path(base) / f"{sp_prefix}_{sp}_out" / _D['d_ident']
        if not d.is_dir():
            continue
        for s in os.listdir(d):
            if (d / s).is_dir():
                sample2dir.setdefault(s, d)
    TOOLS02 = ["blast", "genomad", "metabuli", "rdrpcatch", "viralverify", "virbot"]
    blast_best = {}
    identify = defaultdict(set)
    cache = {}
    for c, sample in sample_of.items():
        d = sample2dir.get(sample)
        if not d:
            continue
        if sample not in cache:
            tc = {}
            for t in TOOLS02:
                f = d / sample / f"{sample}_virus.{t}.result.id"
                if f.exists():
                    tc[t] = f.read_text(encoding="utf-8", errors="ignore")
            vp = d / sample / "blast_output" / f"{sample}.vp.txt"
            vp_content = vp.read_text(encoding="utf-8", errors="ignore") if vp.exists() else ""
            cache[sample] = (tc, vp_content)
        tc, vp_content = cache[sample]
        for t, content in tc.items():
            if c in content:
                identify[c].add(t)
        if c not in blast_best:
            for line in vp_content.splitlines():
                p = line.split("\t")
                if len(p) >= 11 and p[0] == c:
                    blast_best[c] = (p[2], p[3], p[10], p[1])
                    break
    return blast_best, identify


def main():
    args = parse_args()
    base = Path(args.base).expanduser()
    species = [s.strip() for s in args.species.split(",") if s.strip()]
    out = Path(args.out) if args.out else base / "novel_vs_known_final.tsv"

    all_fa = merge_centroids(base, species, args.sp_prefix)
    nt_path = run_blastn(base, all_fa, args.threads, args.rerun_blastn)
    aa_path = run_mmseqs(base, all_fa, args.threads, args.rerun_mmseqs)
    cdd_path = run_cdd(base, all_fa, args.threads, args.rerun_cdd)

    nt_map = load_map(nt_path)
    aa_map = load_map(aa_path)
    contig_cdd = cdd_to_contig(cdd_path, resolve(DB["cdd_annot"]))
    contigs = list(nt_map.keys())

    blast_best, identify = extract_02(base, species, args.sp_prefix, contigs)

    # 05 分类工具
    classify_tool = {}
    for sp in species:
        f = base / f"{args.sp_prefix}_{sp}_out" / _D['d_taxonomy'] / "Votus.integrated" / "final_integrated_classification.tsv"
        if not f.exists():
            continue
        with open(f, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cid = (r.get("contig_id") or "").strip()
                if cid:
                    classify_tool[cid] = r.get("primary_tool", "")

    cols = ["contig_id", "nt_pident", "nt_qcovs", "aa_pident", "aa_species",
            "cdd_top", "final_judge", "blast_pident", "blast_alnlen", "blast_evalue",
            "blast_target", "identify_tools", "classify_tool"]
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for c in sorted(contigs):
            nt = nt_map.get(c, {})
            aa = aa_map.get(c, {})
            doms = contig_cdd.get(c, [])
            from collections import Counter
            cdd_top = ";".join(d for d, _ in Counter(doms).most_common(3))
            b = blast_best.get(c)
            w.writerow({
                "contig_id": c,
                "nt_pident": nt.get("pident", ""),
                "nt_qcovs": nt.get("qcovs", ""),
                "aa_pident": aa.get("pident", ""),
                "aa_species": aa.get("species", ""),
                "cdd_top": cdd_top,
                "final_judge": judge(nt, aa, cdd_top),
                "blast_pident": b[0] if b else "",
                "blast_alnlen": b[1] if b else "",
                "blast_evalue": b[2] if b else "",
                "blast_target": b[3] if b else "",
                "identify_tools": ",".join(sorted(identify.get(c, []))) or "-",
                "classify_tool": classify_tool.get(c, "-"),
            })
    print(f"\n完成: {out} ({len(contigs)} 条)")


if __name__ == "__main__":
    main()
