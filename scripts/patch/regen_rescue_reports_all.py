#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用修好的 _write_rescue_report 重新生成 8 个数据集的 rescue 报告（不重跑流程）。

修复点：
  1) 分支 C 的 Completo 直达 pass（BLASTN qcov>=98%）此前没有对应分支，全部被记成 fail；
  2) 分支 A 计数含 prepass 免拯救 contig，混进候选表的分母，导致「原因数 > 未拯救数」。
用法: python3 regen_rescue_reports.py [数据集名 ...]   (默认全部)
"""
import os
import shutil
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
import rescue_pipeline as rp  # noqa: E402

ROOTS = [
    Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus"),
    Path("/home/zhangwenda/MMPV-paper/onekp-virome"),
]


def find_genus_len():
    for p in (Path(os.path.expanduser("~/database/virus-db/db/genus_lens")),
              Path("/home/zhangwenda/MMPV-RNA/database/genus_lens")):
        if p.is_file():
            return str(p)
    return None


def count_fa(p):
    return sum(1 for _ in SeqIO.parse(str(p), "fasta")) if p and Path(p).is_file() else 0


def main():
    want = set(sys.argv[1:])
    genus_len = find_genus_len()
    print("genus_len:", genus_len)
    rows = []
    for root in ROOTS:
        for ds in sorted(root.iterdir()):
            if not ds.is_dir() or (want and ds.name not in want):
                continue
            out = ds / "08_Rescue" / "Plant"
            if not out.is_dir():
                continue
            cen = ds / "04_CLUSTER" / "rescue_centroids.fasta"
            tax = ds / "05_Taxonomy" / "Votus.integrated" / "final_integrated_classification.tsv"
            if not cen.is_file() or not tax.is_file():
                print("  跳过 %s (缺 centroids 或 taxonomy)" % ds.name)
                continue
            recs = list(SeqIO.parse(str(cen), "fasta"))
            fa = {k: out / ("branch_%s" % k) / ("branch%s_pass.fasta" % k.upper()) for k in "abcd"}
            cnt = {k: count_fa(fa[k]) for k in "abcd"}
            fa = {k: (str(v) if v.is_file() else None) for k, v in fa.items()}
            n_pre = 0
            pf = out / "branch_a_prepass.txt"
            if pf.is_file():
                n_pre = sum(1 for _ in open(pf))
            cnt["a"] += n_pre
            n_final = count_fa(out / "centroids" / "final_centroids.fasta")
            for name in ("rescue_report.tsv", "rescue_summary.md"):
                p = out / name
                if p.is_file() and not Path(str(p) + ".bak_v661").exists():
                    shutil.copy2(str(p), str(p) + ".bak_v661")
            print("\n##### %s  候选=%d A=%d(含prepass %d) B=%d C=%d D=%d final=%d  报告旧版已备份 .bak_v661"
                  % (ds.name, len(recs), cnt["a"], n_pre, cnt["b"], cnt["c"], cnt["d"], n_final))
            rp._write_rescue_report(
                out, recs, {},
                fa["a"], fa["b"], fa["c"], fa["d"],
                cnt["a"], cnt["b"], cnt["c"], cnt["d"], n_final,
                taxonomy_tsv=str(tax), genus_lens_path=genus_len, min_vsi_len=2000,
                prepass=n_pre)
            # 结果核对
            import csv
            from collections import Counter
            bc = Counter()
            with open(out / "rescue_report.tsv") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    bc[r["branch"]] += 1
            print("  新报告 branch 分布: %s" % dict(bc))
            rows.append((ds.name, len(recs), cnt, n_pre, n_final, dict(bc), n_pre))
    print("\n===== 汇总 =====")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
