#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""追查：Species 与 Family 互相矛盾的行到底怎么回事（只读）

做法：用 VMR 的 Species -> Family 权威映射，去比对 8 个 All_plant 表里的 Family 列，
      统计冲突组合、涉及的判定工具，再抽具体 contig 追到上游 05 表。
"""
import csv, os, glob
from collections import Counter, defaultdict

ROOT = os.path.expanduser("~/MMPV-paper")
VMR = os.path.expanduser("~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")

sp2fam, sp2gen, sp2realm = {}, {}, {}
gen2fam = defaultdict(set)
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t"); next(rd)
    for r in rd:
        if len(r) < 27: continue
        sp, fam, ge = r[17].strip().strip('"'), r[13].strip().strip('"'), r[15].strip().strip('"')  # Species / Family / Genus（index 13 才是 Family，12 是 Suborder）
        if sp:
            sp2fam.setdefault(sp, fam); sp2gen.setdefault(sp, ge); sp2realm.setdefault(sp, r[3].strip().strip('"'))
        if ge and fam:
            gen2fam[ge].add(fam)

def rows(p):
    with open(p, newline="", encoding="utf-8", errors="surrogateescape") as f:
        return list(csv.reader(f, delimiter="\t"))

files = []
for pat in ("*/*/*_out/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
            "*/onekp-virus/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"):
    files += sorted(glob.glob(os.path.join(ROOT, pat)))

pairs = Counter(); tools = Counter(); srcs = Counter(); resc = Counter()
selfconsistent = 0   # 表内 Family 与表内 Genus 在 VMR 里对应的科一致（即 Family 跟自己同行的 Genus 自洽）
samples = []
fam_side = Counter(); exp_side = Counter()
n_match = n_mis = n_sp = 0
mixed_rank = 0

for q in files:
    r = rows(q); idx = {h: i for i, h in enumerate(r[0])}
    for x in r[1:]:
        if not x or not x[0] or len(x) < 13: continue
        sp = x[idx["Species"]].strip(); fam = x[idx["Family"]].strip()
        if sp not in sp2fam: continue
        n_sp += 1
        vf = sp2fam[sp]
        if fam == vf:
            n_match += 1
        elif fam and fam not in ("NA",):
            n_mis += 1
            ge = x[idx["Genus"]].strip()
            gf = gen2fam.get(ge, set())
            if len(gf) == 1 and fam in gf:
                selfconsistent += 1
            pairs[(fam, vf)] += 1
            fam_side[fam] += 1; exp_side[vf] += 1
            tools[x[idx["primary_tool"]]] += 1
            srcs[x[idx["source"]]] += 1
            resc[x[idx["rescued"]]] += 1
            if len(samples) < 12:
                samples.append((q.split("MMPV-paper/")[-1].split("/")[-3], x[idx["contig_id"]], sp, fam, vf,
                                x[idx["primary_tool"]], x[idx["rescued"]], x[idx["confidence"]],
                                x[idx["ref_cov"]], x[idx["n_samples"]],
                                x[idx["Realm"]], x[idx["Kingdom"]]))
        # 高等级列内部是否自相矛盾（Realm 与 Kingdom 不同域）
        rm, kd = x[idx["Realm"]].strip(), x[idx["Kingdom"]].strip()
        if (rm == "Riboviria" and kd in ("Bamfordvirae", "Loebvirae", "Shotokuvirae", "Trapavirae", "Sangervirae")) or \
           (rm in ("Varidnaviria", "Monodnaviria") and kd in ("Orthornavirae", "Pararnavirae")):
            mixed_rank += 1

print("有 VMR 物种可核对的样本行: %d" % n_sp)
print("  Family 与 VMR 一致: %d" % n_match)
print("  Family 与 VMR 冲突: %d  (占可核对行 %.1f%%)" % (n_mis, 100.0 * n_mis / max(n_sp, 1)))
print("  其中 表内 Family 与表内 Genus 在 VMR 里对应同一个科(即科随属、但与种冲突)的行: %d" % selfconsistent)
print("  高等级列自相矛盾(Realm 与 Kingdom 不同域)的行: %d" % mixed_rank)

print("\n=== 冲突组合 top 20 (表里 Family -> VMR 该物种的真实 Family) ===")
for (f1, f2), c in pairs.most_common(20):
    print("  %-22s -> %-24s %d" % (f1, f2, c))

print("\n=== 表内 Family 侧 top 10 ===")
for k, v in fam_side.most_common(10): print("  %-26s %d" % (k, v))
print("=== VMR 真实 Family 侧 top 10 ===")
for k, v in exp_side.most_common(10): print("  %-26s %d" % (k, v))
print("=== primary_tool 分布 ===", dict(tools))
print("=== source 分布 ===", dict(srcs))
print("=== rescued 分布 ===", dict(resc))

print("\n=== 样例 12 行 ===")
print("%-16s %-46s %-24s %-20s %-22s %-12s %-4s %-5s %-8s %-3s %-12s %-13s" % (
    "项目", "contig_id", "Species", "表内Family", "VMR真实Family", "tool", "resc", "conf", "ref_cov", "n", "Realm", "Kingdom"))
for s in samples:
    print("%-16s %-46s %-24s %-20s %-22s %-12s %-4s %-5s %-8s %-3s %-12s %-13s" % (
        s[0][:16], s[1][:46], s[2][:24], s[3][:20], s[4][:22], s[5][:12], s[6], s[7], s[8], s[9], s[10][:12], s[11][:13]))
