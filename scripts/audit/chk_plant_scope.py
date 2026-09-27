#!/usr/bin/env python3
"""「这种问题」在植物病毒子集里有多少，以及每一行是什么性质

基线统一用 R 的两份输出（口径一致，不自己重算不变量）：
  before = /tmp/rc_base   现行代码（无科-属校准）
  after  = /tmp/rc_live   已部署代码（含科-属校准）
change 集 = before 属非空 且 after 属不同 = R 日志里的「替换 + 置空」

受影响行分四类（用数据化的植物病毒科集合判，不硬编码）：
  G1 外来污染      : 属的 VMR 科不在植物病毒科集合，且共识科在
  G2 共识科非植物  : 共识科本身不在植物病毒科集合
  G3 植物属被改    : 属的 VMR 科在植物病毒科集合，且与共识科不同
  G4 其他
"""
import csv
import importlib.util
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

DT = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
BEFORE = Path("/tmp/rc_base/final_integrated_classification.tsv")
AFTER = Path("/tmp/rc_live/final_integrated_classification.tsv")
PLANT_TAX = DT / "10_Reports/plant_final_taxonomy.tsv"
FASTAS = {
    "植物病毒 All_plant": DT / "09_Virome_Analysis/All_plant.viruses.fasta",
    "HQ 植物病毒": DT / "09_Virome_Analysis/HQ_plant_viruses.fasta",
}


def blank(v):
    return (not v) or str(v).strip() in ("", "NA", "na", "nan", "None", "-")


def clean(v):
    return (v or "").strip().strip('"')


def load(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fa_ids(p):
    ids = set()
    with open(p, encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith(">"):
                ids.add(ln[1:].split()[0].strip().strip('"'))
    return ids


spec = importlib.util.spec_from_file_location("ann", PIPE / "utils" / "annotate_nucleic_acid.py")
ann = importlib.util.module_from_spec(spec); spec.loader.exec_module(ann)
IDX = ann.load_vmr(VMR)

# 属 -> VMR 科
g2fam = {}
for r in load(VMR):
    g = clean(r.get("Genus"))
    if g:
        g2fam.setdefault(g.lower(), Counter())[clean(r.get("Family")) or "(空)"] += 1

# 本项目的植物病毒科集合（数据化：取植物病毒分类表里出现过的科）
plant_fams = {clean(r.get("Family")) for r in load(PLANT_TAX) if clean(r.get("Family"))}
plant_fams.discard("NA")
print("植物病毒科集合 (%d 个): %s\n" % (len(plant_fams), ", ".join(sorted(plant_fams))))

before = {r["contig_id"]: r for r in load(BEFORE)}
after = {r["contig_id"]: r for r in load(AFTER)}
print("before = %s (%d 行)" % (BEFORE, len(before)))
print("after  = %s (%d 行)\n" % (AFTER, len(after)))

subsets = [("全部行", None)]
for label, p in FASTAS.items():
    subsets.append((label, fa_ids(p)))
    print("%-18s ID 数 %d" % (label, len(subsets[-1][1])))
print()


def famof(gen):
    c = g2fam.get(gen.lower())
    return c.most_common(1)[0][0] if c else None


subset_rows = {}
det = defaultdict(lambda: {"G3": [], "G1": []})
for label, ids in subsets:
    rows = [before[c] for c in before if ids is None or c in ids]
    subset_rows[label] = rows
    n_conf_b = n_conf_a = n_rep = n_blank = 0
    groups = Counter()
    for r in rows:
        cid = r["contig_id"]
        bf = ann.resolve_rank(clean(r.get("Family")), "Family", IDX)
        bg = ann.resolve_rank(clean(r.get("Genus")), "Genus", IDX)
        if bf and bg and bf != bg:
            n_conf_b += 1
        a = after[cid]
        af = ann.resolve_rank(clean(a.get("Family")), "Family", IDX)
        ag = ann.resolve_rank(clean(a.get("Genus")), "Genus", IDX)
        if af and ag and af != ag:
            n_conf_a += 1
        g0, g1v = clean(r.get("Genus")), clean(a.get("Genus"))
        if blank(g0) or g0 == g1v:
            continue
        fam_c = clean(r.get("Family")) or "NA"
        vf = famof(g0)
        is_blanked = blank(g1v)
        if is_blanked:
            n_blank += 1
        else:
            n_rep += 1
        if vf is None:
            g = "G4 属未收录"
        elif fam_c not in plant_fams:
            g = "G2 共识科非植物科"
        elif vf in plant_fams:
            g = "G3 植物病毒属被改"
        else:
            g = "G1 外来污染属"
        groups[g] += 1
        rec = (cid, fam_c, g0, vf, "置空" if is_blanked else "替换为 " + g1v)
        if g == "G3 植物病毒属被改":
            det[label]["G3"].append(rec)
        elif g == "G1 外来污染属":
            det[label]["G1"].append(rec)
    print("=" * 92)
    print("子集: %s   (行数 %d)" % (label, len(rows)))
    print("  科属矛盾: before %d  ->  after %d" % (n_conf_b, n_conf_a))
    print("  受影响行: %d   (替换 %d / 置空 %d)" % (n_rep + n_blank, n_rep, n_blank))
    for k in sorted(groups):
        print("      %-20s %5d" % (k, groups[k]))
    # 汇总每类的 共识科 -> 属(VMR科) 组合
    combo = Counter()
    for r in rows:
        cid = r["contig_id"]
        g0 = clean(r.get("Genus")); g1v = clean(after[cid].get("Genus"))
        if blank(g0) or g0 == g1v:
            continue
        combo[("%s -> %s (VMR: %s)" % (clean(r.get("Family")) or "NA", g0,
                                       famof(g0) or "未收录"))] += 1
    print("  受影响行的 共识科 -> 属(VMR科) top12:")
    for k, v in combo.most_common(12):
        print("      %-62s %4d" % (k, v))
    print()

for label in [l for l, _ in subsets if l != "全部行"]:
    print("\n" + "=" * 92)
    print("【%s】G3 明细：共识科是植物病毒科，但属来自另一个植物病毒科，属被改掉" % label)
    print("=" * 92)
    for cid, fam_c, g0, vf, act in det[label]["G3"]:
        print("  %-46s %-20s -> %-20s (VMR: %-18s) [%s]"
              % (cid.replace("_clean_", "|"), fam_c, g0, vf, act))
    print("\n【%s】G1 明细前 20：属来自非植物病毒科，判为外来污染" % label)
    for cid, fam_c, g0, vf, act in det[label]["G1"][:20]:
        print("  %-46s %-20s -> %-20s (VMR: %-18s) [%s]"
              % (cid.replace("_clean_", "|"), fam_c, g0, vf, act))
