#!/usr/bin/env python3
"""该 contig 按现行代码/脚本的完整分类结果。

修正点（前一版错的）：
  1. decide() 的语义是「8 层各出一票取多数派」，不是「按 LOOKUP_ORDER 命中即停」。
  2. resolve_rank 返回 None 有三种不同原因，之前一律标成「VMR 无此名」，现分开报：
       空白 / VMR 未收录该名 / 真跨大界（该层名下有 DNA 与 RNA 两派，按设计弃票）
"""
import csv
import importlib.util
from collections import Counter
from pathlib import Path

PREFIX = "CRR1126135_clean_NODE_68_length_2879_cov_8"
PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
FILES = [
    ("磁盘现行产物(8-12, 未重跑)", Path(
        "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
        "/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")),
    ("已部署脚本重跑(校准后)", Path("/tmp/rc_live/final_integrated_classification.tsv")),
]

spec = importlib.util.spec_from_file_location("ann", PIPE / "utils" / "annotate_nucleic_acid.py")
ann = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ann)
IDX = ann.load_vmr(VMR)


def why_none(v, rank):
    """resolve_rank 返回 None 的三种原因，分开报"""
    if ann.is_blank(v):
        return "空白"
    gc = IDX["lv"][rank].get(v.strip())
    if not gc:
        return "VMR 未收录该名"
    dna = sorted({k for k in gc if ann.to_broad(k) == "DNA"})
    rna = sorted({k for k in gc if ann.to_broad(k) == "RNA"})
    s = "真跨大界(按设计弃票): DNA侧=" + ("/".join(dna) if dna else "无")
    s += "; RNA侧=" + ("/".join(rna) if rna else "无")
    return s


for label, p in FILES:
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        hit = [r for r in csv.DictReader(f, delimiter="\t") if r["contig_id"].startswith(PREFIX)]
    if not hit:
        print("[%s] 未命中\n" % label)
        continue
    row = hit[0]
    print("=" * 78)
    print("%s" % label)
    print("=" * 78)
    print("  contig_id = %s" % row["contig_id"])
    print("  primary_tool = %s   completeness = %s   confidence = %s"
          % (row.get("primary_tool"), row.get("completeness"), row.get("confidence")))
    print("\n  ── 分类谱系（最终表里的值）──")
    for rank in ann.RANKS:
        print("      %-8s %s" % (rank, row.get(rank) or "NA"))
    print("\n  ── 核酸判定：8 层各一票，取多数派 ──")
    per = {}
    for rank in ann.RANKS:
        v = row.get(rank, "") or ""
        b = ann.resolve_rank(v, rank, IDX)
        if b:
            per[rank] = b
            gc = IDX["lv"][rank].get(v.strip(), {})
            det = ", ".join("%s x%d" % (k, n) for k, n in Counter(gc).most_common(3))
            print("      %-8s %-22s -> %-3s 票      [VMR 该层名下: %s]" % (rank, v or "NA", b, det))
        else:
            print("      %-8s %-22s -> 弃票   (%s)" % (rank, v or "NA", why_none(v, rank)))
    broad, via = ann.decide(row, IDX)
    tally = Counter(per.values())
    print("\n      计票: %s   -> 多数派 %s" % (dict(tally), broad if broad != "NA" else "无"))
    print("      Nucleic_acid = %s   (日志层级 = %s)" % (broad, via))
    print()
