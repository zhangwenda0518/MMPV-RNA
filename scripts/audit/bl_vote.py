"""抽查: 带原生/藻科名的 Plant 行, 其宿主判定究竟靠哪个字段.
若 Species 名是植物病毒种名 -> 宿主票来自 Species 表, 科字段是坏标签;
若 Species 名也是巨型病毒种名 -> 宿主票本身可疑.
"""
import csv
from collections import Counter

K = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus"
G = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
C7 = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"


def load(tbl, col):
    d = {}
    with open(f"{C7}/{tbl}", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[(r.get(col) or "").strip()] = ((r.get("Predicted_Host") or "").strip(),
                                             (r.get("Confidence_Level") or "").strip(),
                                             (r.get("Total_Records") or "").strip())
    return d


sp = load("species_host_probability.tsv", "Species")
gen = load("genus_host_probability.tsv", "Genus")
fam = load("family_host_probability.tsv", "Family")
PA = ("Protist", "Algae")

for label, root in (("onekp", K), ("goji", G)):
    print("\n" + "=" * 100)
    print(f"[{label}] 带原生/藻科或属名 且 Final_Host=Plant 的行: 逐条看宿主票从哪来")
    print("=" * 100)
    lvl = {}
    with open(f"{root}/06_HostPrediction/C9_ICTV_result/classification_result.tsv", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            lvl[r["contig_id"]] = r
    rows = []
    with open(f"{root}/06_HostPrediction/ensemble_host_summary.tsv", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if (r.get("Final_Host") or "").strip() != "Plant":
                continue
            fm = (r.get("Family") or "").strip()
            g = (r.get("Genus") or "").strip()
            if fam.get(fm, ("",))[0] in PA or gen.get(g, ("",))[0] in PA:
                rows.append(r)
    print(f"共 {len(rows)} 行")
    # 统计: species 名在 species 表命中时的宿主
    vote = Counter()
    for r in rows:
        c9 = lvl.get(r["contig_id"], {})
        spname = (c9.get("Species") or "").strip()
        lv = (c9.get("Determination_Level") or "").strip()
        real = lv.split("(")[0]
        if real == "Species":
            h = sp.get(spname, ("<species表未命中>",))[0]
        elif real == "Genus":
            h = gen.get((c9.get("Genus") or "").strip(), ("<genus表未命中>",))[0]
        elif real == "Family":
            h = fam.get((c9.get("Family") or "").strip(), ("<family表未命中>",))[0]
        else:
            h = "?"
        vote[(real, h)] += 1
    print("  宿主票来源层级 -> 该层级查表判定的宿主:")
    for k, v in vote.most_common():
        print(f"    {str(k):<40} {v}")

    # 明细抽样
    print("\n  抽样 12 行:")
    print(f"  {'Family':<20}{'Genus':<22}{'Species':<32}{'Species表判定':<14}{'层级'}")
    shown = 0
    for r in rows:
        c9 = lvl.get(r["contig_id"], {})
        spname = (c9.get("Species") or "").strip()
        lv = (c9.get("Determination_Level") or "").strip()
        if lv.split("(")[0] != "Species":
            continue
        print(f"  {(r.get('Family') or '')[:19]:<20}{(r.get('Genus') or '')[:21]:<22}"
              f"{spname[:31]:<32}{sp.get(spname, ('<未命中>',))[0]:<14}{lv}")
        shown += 1
        if shown >= 12:
            break
