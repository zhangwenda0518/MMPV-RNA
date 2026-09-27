#!/usr/bin/env python3
# 量化 fill_taxonomy_na 的实际改动：从 per-tool 产物重建 fill 之前的 combined，与磁盘上 fill 之后的逐行对齐比对
# 用法: python3 /tmp/measure_fill_impact.py [label ...]
import sys, os, glob, collections

sys.path.insert(0, "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
import virus_classifier as vc

RANK = list(vc.RANK_NAMES)
TF_SUFFIX = {
    "genomad": "_genomad_taxonomy.tsv",
    "metabuli": "_metabuli_taxonomy.tsv",
    "diamond_lca": "_diamond_lca_taxonomy.tsv",
    "CAT": "_CAT_taxonomy.tsv",
    "mmseqs": "_mmseqs_taxonomy.tsv",
    "VITAP": "_VITAP_taxonomy.tsv",
    "ACVirus": "_ACVirus_taxonomy.tsv",
    "vcontact3": "_vcontact3_taxonomy.tsv",
}
ROOTS = {
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
    "Alternaria": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Alternaria_alternata_out",
    "Aphis": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Aphis_gossypii_out",
    "Fusarium": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Fusarium_nematophilum_out",
    "amarum": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_amarum_out",
    "barbarum": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "chinense": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_chinense_out",
    "ruthenicum": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_ruthenicum_out",
}


def classed_dir(label):
    for sub in ("05_Taxonomy/Votus.classed", "05_Taxonomy"):
        d = os.path.join(ROOTS[label], sub)
        if glob.glob(os.path.join(d, "*_combined_taxonomy.tsv")):
            return d
    raise SystemExit("no combined found for " + label)


def gen_rebuilt(sample, outdir, tools):
    for tool in tools:
        tf = os.path.join(outdir, sample + TF_SUFFIX[tool])
        if not os.path.exists(tf):
            continue
        with open(tf) as f:
            lines = f.readlines()
        has_hdr = bool(lines) and lines[0].startswith("seq_name")
        for line in lines[int(has_hdr):]:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split('\t')
            if len(parts) >= 5:
                r = vc.lineage_to_ranks(parts[4])
            elif len(parts) >= 3:
                r = vc.lineage_to_ranks(parts[2])
            else:
                continue
            yield [parts[0], tool] + [r.get(rn, "NA") for rn in RANK]


def measure(label):
    d = classed_dir(label)
    comb = glob.glob(os.path.join(d, "*_combined_taxonomy.tsv"))[0]
    sample = os.path.basename(comb).split("_combined_taxonomy.tsv")[0]
    # 从磁盘文件恢复 tool 出现顺序
    tools, seen = [], set()
    with open(comb) as f:
        f.readline()
        for line in f:
            p = line.split('\t')
            if len(p) > 1 and p[1] not in seen:
                seen.add(p[1]); tools.append(p[1])
    print("=" * 70)
    print("[%s] sample=%s tools=%s" % (label, sample, ",".join(tools)))

    overwrite = collections.Counter()   # (rank, old, new) 旧值非 NA 被改
    filled = collections.Counter()      # (rank, new) 旧值 NA 被填
    ex = []
    n_rows = n_chg_rows = n_ov_rows = 0
    with open(comb) as f:
        next(f)
        for before, line in zip(gen_rebuilt(sample, d, tools), f):
            after = line.strip().split('\t')
            n_rows += 1
            row_chg = row_ov = False
            for i, rn in enumerate(RANK):
                o = before[2 + i]
                n = after[2 + i].strip().strip('"')
                if o == n:
                    continue
                row_chg = True
                if o in ("", "NA"):
                    filled[(rn, n)] += 1
                else:
                    overwrite[(rn, o, n)] += 1
                    row_ov = True
                    if len(ex) < 6:
                        ex.append((before[0], before[1], rn, o, n, " | ".join(before[2:])))
            if row_chg:
                n_chg_rows += 1
            if row_ov:
                n_ov_rows += 1
    print("  行数=%d  有改动行=%d  其中覆盖非NA行=%d" % (n_rows, n_chg_rows, n_ov_rows))
    print("  覆盖事件(非NA→不同值) 合计=%d  top15:" % sum(overwrite.values()))
    for (rn, o, n), c in overwrite.most_common(15):
        print("    %-8s %-34s -> %-34s x%d" % (rn, o[:34], n[:34], c))
    print("  回填事件(NA→值) 合计=%d  top15:" % sum(filled.values()))
    for (rn, n), c in filled.most_common(15):
        print("    %-8s -> %-40s x%d" % (rn, n[:40], c))
    print("  覆盖样例:")
    for e in ex:
        print("    %s [%s] %s: %s -> %s" % (e[0][:40], e[1], e[2], e[3][:30], e[4][:30]))
        print("       行内原值: %s" % e[5][:200])


if __name__ == "__main__":
    for lb in (sys.argv[1:] or list(ROOTS)):
        try:
            measure(lb)
        except Exception as exc:
            print("[%s] ERROR %s" % (lb, exc))
