#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成逐级相容性约束用的「属 -> 定型科」参照表（供 virus_classifier_analysis.R 使用）。

两张参照：
  1. NCBI taxdump rankedlineage.dmp
     列序（\\t|\\t 分隔，10 列）：
       tax_id | tax_name | species | genus | family | order | class | phylum | kingdom | superkingdom
     病毒条目里 superkingdom 列放的是 Realm，kingdom 列放的是 Kingdom；每行只填严格高于自身等级的祖级。
     口径与 scripts/audit/apply_calib_A.py 的 ref1 完全一致：按文件顺序取该 name 的**首次出现**，
     family 列填了就用，没填就记为未知（后面即使再出现也不覆盖）。
  2. VMR MSL41 taxa.txt（逗号分隔，含 Virus GENBANK accession 这类带逗号的引号字段）
     Genus + Family 两列都非空才记一条 genus -> family；同一属对应多个科时用 ";" 连接。

输出 TSV（表头）：Genus  NCBI_Family  NCBI_n  VMR_Family  VMR_n  Domain
  NCBI_n  = 该 name 在 dmp 全部出现行里出现过的不同科数（>1 表示同名多科，供诊断，不参与判定）
  VMR_n   = 该属在 taxa.txt 里关联的科数
  Domain  = dmp 首次出现行的 superkingdom 值（病毒条目即 Realm），供诊断

默认剔除 Domain ∈ {Bacteria, Archaea, Eukaryota} 且不在 VMR 属名单里的条目：病毒分类管线里
出现的属不应该是细胞生物的属，留着只会引入同名污染。本数据集实测剔除前后判定结果完全一致。
加 --keep-all 可关掉该剔除。只读参照库，只写输出文件；不触碰任何管线产物。
"""
import csv
import hashlib
import os
import sys
from collections import defaultdict

DMP = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
VMR = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"
OUTDIR = "/home/zhangwenda/database/taxonomy"
OUT = os.path.join(OUTDIR, "genus_family_ref.tsv")
META = os.path.join(OUTDIR, "genus_family_ref.meta.txt")

# NCBI 分类体系里的三个细胞生物 superkingdom；这三个之外的 superkingdom（含 Realm 名与空值）都保留
CELLULAR = {"bacteria", "archaea", "eukaryota"}
KEEP_ALL = "--keep-all" in sys.argv


def md5_of(path, limit=None):
    h = hashlib.md5()
    n = 0
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            h.update(b)
            n += len(b)
            if limit is not None and n >= limit:
                break
    return h.hexdigest()


def norm(v):
    if v is None:
        return ""
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if v.upper() in ("NA", "N/A"):
        return ""
    return v


def load_ncbi():
    """返回 first_fam: name -> 首次出现的科(可能为 '')，以及 distinct: name -> set(科)"""
    first_fam = {}
    domain = {}
    distinct = defaultdict(set)
    n_lines = 0
    with open(DMP, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = [p.strip().strip("|").strip() for p in line.rstrip("\n").split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if not nm:
                continue
            n_lines += 1
            fam = parts[4]
            if nm not in first_fam:
                first_fam[nm] = fam
                domain[nm] = parts[9]
                if fam:
                    distinct[nm].add(fam)
            elif nm in distinct and fam:
                distinct[nm].add(fam)
    print("  dmp 有效行 %d / 唯一 name %d / 首次出现带科者 %d"
          % (n_lines, len(first_fam), sum(1 for v in first_fam.values() if v)))
    return first_fam, domain, distinct


def load_vmr():
    gen2fam = defaultdict(set)
    n = 0
    with open(VMR, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter=",")
        tc = {c.strip().upper(): c for c in (rd.fieldnames or [])}
        cf, cg = tc.get("FAMILY"), tc.get("GENUS")
        if not cf or not cg:
            raise SystemExit("ERROR: taxa.txt 表头异常 %s" % (rd.fieldnames,))
        for r in rd:
            g = norm(r.get(cg))
            f = norm(r.get(cf))
            if g and f:
                gen2fam[g.lower()].add(f)
                n += 1
    print("  VMR 映射 %d 条 / 属 %d 个" % (n, len(gen2fam)))
    return gen2fam


def main():
    for p in (DMP, VMR):
        if not os.path.exists(p):
            print("ERROR: 找不到参照库 %s" % p)
            return 1
    if not os.path.isdir(OUTDIR):
        print("ERROR: 输出目录不存在 %s" % OUTDIR)
        return 1

    print("读取 NCBI rankedlineage.dmp ...")
    first_fam, domain, distinct = load_ncbi()
    print("读取 VMR taxa.txt ...")
    gen2fam = load_vmr()

    # 名称展示用的大小写：优先 NCBI 原名，其次 VMR 原名
    names = set(first_fam) | set(gen2fam)
    n_ncbi_only = 0
    n_vmr_only = 0
    n_both = 0
    n_written = 0
    n_conflict = 0
    n_trimmed = 0
    domain_counter = defaultdict(int)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        fh.write("Genus\tNCBI_Family\tNCBI_n\tVMR_Family\tVMR_n\tDomain\n")
        for nm in sorted(names):
            nfam = first_fam.get(nm, "")
            vs = gen2fam.get(nm)
            if not nfam and not vs:
                continue
            dom = domain.get(nm, "")
            if KEEP_ALL is False and dom.lower() in CELLULAR and not vs:
                n_trimmed += 1
                continue
            vmr_str = ";".join(sorted(vs)) if vs else ""
            if nfam and vs:
                n_both += 1
                if nfam not in vs:
                    n_conflict += 1
            elif nfam:
                n_ncbi_only += 1
            else:
                n_vmr_only += 1
            fh.write("%s\t%s\t%d\t%s\t%d\t%s\n" % (
                nm, nfam, len(distinct.get(nm, ())), vmr_str, len(vs) if vs else 0,
                dom))
            n_written += 1
            if dom:
                domain_counter[dom] += 1

    print("")
    print("=== 写入 %s ===" % OUT)
    print("  数据行 %d (剔除细胞生物同名条目 %d)" % (n_written, n_trimmed))
    print("  仅 NCBI 有科 %d / 仅 VMR 有科 %d / 两边都有 %d" % (n_ncbi_only, n_vmr_only, n_both))
    print("  两边都有但科不一致 %d（多数应为分类学版本漂移）" % n_conflict)
    print("  大小 %.1f MB" % (os.path.getsize(OUT) / 1048576.0))
    print("  按 Domain(superkingdom) 统计前 12:")
    for k, v in sorted(domain_counter.items(), key=lambda x: -x[1])[:12]:
        print("    %-28s %7d" % (k, v))

    with open(META, "w", encoding="utf-8") as fh:
        fh.write("file: %s\n" % OUT)
        fh.write("rows: %d\n" % n_written)
        fh.write("source_ncbi: %s\n" % DMP)
        fh.write("source_ncbi_md5: %s\n" % md5_of(DMP))
        fh.write("source_vmr: %s\n" % VMR)
        fh.write("source_vmr_md5: %s\n" % md5_of(VMR))
        fh.write("generator: scripts/audit/build_genus_family_ref.py\n")
        fh.write("genus_ncbi_only: %d\n" % n_ncbi_only)
        fh.write("genus_vmr_only: %d\n" % n_vmr_only)
        fh.write("genus_both: %d\n" % n_both)
        fh.write("genus_both_family_mismatch: %d\n" % n_conflict)
        fh.write("trimmed_cellular: %d\n" % n_trimmed)
        fh.write("keep_all: %s\n" % KEEP_ALL)
    print("  provenance: %s" % META)
    return 0


if __name__ == "__main__":
    sys.exit(main())
