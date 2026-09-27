#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读预演：把"每阶元独立取加权票首"换成"共识优先，不一致/无值才按工具优先级兜底"，并去掉最完整工具行兜底。

复刻依据：virus_classifier_analysis.R v6.2（md5 3e6c0cd8）
  :30    TAX_LEVELS
  :39-40 RANK_DEPTH_WEIGHTS / TOOL_BIAS
  :107   is_valid_value_vec
  :112   species_quality_score_vec
  :138   clean_all_ranks（投票前过滤）
  :251   compute_tool_weights
  :511   build_consensus（第 3 步加权投票、第 4 步最完整工具行兜底）
  :981-1017 consensus_stats 自举（未加权计数，common_ids 全工具共有）
不写任何产物，只打印。
"""
import csv, os, re, sys
from collections import defaultdict

BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
        "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated")

TAX = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
RDW = {"Realm": 1, "Kingdom": 2, "Phylum": 4, "Class": 8,
       "Order": 16, "Family": 32, "Genus": 64, "Species": 128}
BIAS = {"ACVirus": 1.2, "VITAP": 1.1, "mmseqs": 1.0, "metabuli": 1.0,
        "CAT": 0.9, "genomad": 0.9, "diamond_lca": 0.8}
PRIORITY = ["ACVirus", "VITAP", "mmseqs", "metabuli", "CAT", "genomad", "diamond_lca"]
PRIO_RANK = {t: i for i, t in enumerate(PRIORITY)}
KNOWN_REALMS = {"riboviria", "monodnaviria", "duplodnaviria", "varidnaviria",
                "adnaviria", "ribozyviria"}
SUBRANK = ["viricotina", "viricetidae", "virineae", "virinae"]
INVALID = {"-", "NA", "na", "N/A", "no rank", "undefined", "unknown", "null",
           "default", "Unclassified"}


def valid(v):
    return v is not None and v != "" and v not in INVALID


def species_quality(s):
    v = s.lower()
    if re.search(r"(inae|idae|ales|icetes|viricota)$", v):
        return 0.0
    q = 1.0
    if re.search(r"\bsp\.?\b|\bcf\.?\b|\baff\.?\b", v):
        q *= 0.3
    if "unclassified" in v or "environmental" in v or "uncultured" in v:
        q *= 0.1
    if v.startswith("unplaced") or v.startswith("novel_"):
        q *= 0.2
    if re.match(r"^[A-Z][a-z]+virus\s+[a-z]", s):
        q *= 1.5
    if " " not in s:
        q *= 0.5
    return q


def clean_rank(rank, val):
    if not valid(val):
        return None
    v = val.lower()
    if rank == "Realm":
        return val if v in KNOWN_REALMS else None
    if rank == "Phylum":
        return val if v.endswith("viricota") else None
    if rank == "Class":
        return val if v.endswith("viricetes") else None
    if rank == "Order":
        return val if v.endswith("virales") else None
    if rank == "Family":
        return val if v.endswith("viridae") else None
    if rank == "Genus":
        return None if (v.endswith("viridae") or v.endswith("virinae")) else val
    if rank == "Species":
        return None if re.search(r"(inae|idae|ales|icetes|viricota)$", v) else val
    return val


def load_tool(path):
    out = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = r.get("contig_id")
            if not cid:
                continue
            row = {}
            for k in TAX:
                row[k] = clean_rank(k, (r.get(k) or "").strip())
            for k in list(row):
                if row[k] and any(row[k].lower().endswith(s) for s in SUBRANK):
                    row[k] = None
            cur = out.get(cid)
            if cur is None:
                out[cid] = row
            else:  # 复刻 R: unique(by=contig_id) 取首行
                pass
    return out


tools = []
for t in PRIORITY:
    p = os.path.join(BASE, "standardized_%s.tsv" % t)
    if os.path.exists(p):
        tools.append(t)
data = {t: load_tool(os.path.join(BASE, "standardized_%s.tsv" % t)) for t in tools}
print("=== 载入的工具 (%d) ===" % len(tools))
for t in tools:
    print("  %-12s %d contig" % (t, len(data[t])))

all_ids = set()
for t in tools:
    all_ids |= set(data[t])
print("  contig 并集 %d" % len(all_ids))

# ---- consensus_stats 自举（未加权计数） ----
cnt_tool = defaultdict(int)
for cid in all_ids:
    for t in tools:
        if cid in data[t]:
            cnt_tool[cid] += 1
n_tools = len(tools)
common = [c for c in all_ids if cnt_tool[c] == n_tools]
ids_w = common if common else [c for c in all_ids if cnt_tool[c] >= max(2, n_tools // 2)]
print("\n=== 权重自举样本 ===")
print("  common_ids(全工具共有) = %d ; 自举用 %d 条" % (len(common), len(ids_w)))

boot = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))  # cid->rank->taxon->count
for cid in ids_w:
    for t in tools:
        row = data[t].get(cid)
        if not row:
            continue
        for rank in TAX:
            if row[rank]:
                boot[cid][rank][row[rank]] += 1

first_seen = {}
for i, cid in enumerate(ids_w):
    for t in tools:
        row = data[t].get(cid)
        if not row:
            continue
        for rank in TAX:
            if row[rank]:
                first_seen.setdefault((cid, rank, row[rank]), len(first_seen))

boot_cons = {}
for cid in boot:
    for rank in boot[cid]:
        items = sorted(boot[cid][rank].items(), key=lambda kv: (-kv[1], first_seen.get((cid, rank, kv[0]), 0)))
        boot_cons[(cid, rank)] = items[0][0]

agree = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # tool->rank->[agreed,total]
for cid in boot:
    for rank in boot[cid]:
        cons = boot_cons[(cid, rank)]
        for t in tools:
            row = data[t].get(cid)
            if not row or not row[rank]:
                continue
            agree[t][rank][1] += 1
            if row[rank].lower() == cons.lower():
                agree[t][rank][0] += 1

print("\n=== 各工具在各阶元的自举一致率 (agreed/total) ===")
print("  %-12s %s" % ("tool", " ".join("%-9s" % r[:8] for r in TAX)))
for t in tools:
    line = []
    for r in TAX:
        a, b = agree[t][r]
        line.append("%-9s" % ("%.2f" % (a / b) if b else "n/a"))
    print("  %-12s %s" % (t, " ".join(line)))

W = {}
for t in tools:
    W[t] = {}
    for r in TAX:
        w = RDW[r] * BIAS.get(t, 0.8)
        a, b = agree[t][r]
        if b:
            rate = a / b
            if rate > 0:
                w *= rate
        W[t][r] = w

# ---- 每 contig 每阶元收集 (taxon, weight) ----
cells = defaultdict(lambda: defaultdict(list))
for cid in all_ids:
    for t in tools:
        row = data[t].get(cid)
        if not row:
            continue
        for r in TAX:
            v = row[r]
            if v:
                w = W[t][r]
                if r == "Species":
                    w *= species_quality(v)
                cells[cid][r].append((v, w, t))

MISSING = "<空>"
old_wide, new_strict, new_majority = {}, {}, {}
census = {r: defaultdict(int) for r in TAX}
census_tool = {r: defaultdict(int) for r in TAX}  # 兜底时用的是哪个工具
for cid in all_ids:
    old_wide[cid], new_strict[cid], new_majority[cid] = {}, {}, {}
    for r in TAX:
        recs = cells[cid][r]
        if not recs:
            census[r]["无值"] += 1
            old_wide[cid][r] = new_strict[cid][r] = new_majority[cid][r] = MISSING
            continue
        agg = defaultdict(float)
        first = {}
        for i, (v, w, t) in enumerate(recs):
            agg[v] += w
            first.setdefault(v, i)
        order = sorted(agg.items(), key=lambda kv: (-kv[1], first[kv[0]]))
        top = order[0][0]
        tot = sum(agg.values())
        ratio = agg[top] / tot

        # 兜底：工具优先级最高的、报了该阶元的工具
        fb = None
        for t in PRIORITY:
            for v, w, tt in recs:
                if tt == t:
                    fb = v
                    break
            if fb is not None:
                break

        old_wide[cid][r] = top
        if len(agg) == 1:
            census[r]["全票一致"] += 1
            new_strict[cid][r] = new_majority[cid][r] = top
        elif ratio > 0.5:
            census[r]["多数(>50%)"] += 1
            new_strict[cid][r] = fb
            new_majority[cid][r] = top
            census_tool[r][fb] += 1
        else:
            census[r]["分歧"] += 1
            new_strict[cid][r] = fb
            new_majority[cid][r] = fb
            census_tool[r][fb] += 1

print("\n=== 1) 各阶元的票型普查（每 contig 每阶元一个格子） ===")
print("  %-9s %9s %9s %11s %8s" % ("阶元", "全票一致", "多数>50%", "分歧", "无值"))
for r in TAX:
    c = census[r]
    print("  %-9s %9d %9d %11d %8d" % (r, c["全票一致"], c["多数(>50%)"], c["分歧"], c["无值"]))
tc = sum(census[TAX[0]].values())
print("  （阶元间行数不同属正常：低阶元报值的 contig 少）")

print("\n=== 2) 兜底时实际被选中的工具分布（分歧/多数格） ===")
for r in ("Family", "Genus", "Species"):
    tot = sum(census_tool[r].values())
    if tot:
        print("  %-9s %s" % (r, dict(sorted(census_tool[r].items(), key=lambda kv: -kv[1]))))

print("\n=== 3) 与现行「独立取票首」相比，新规则改了多少格 ===")
for name, new in (("严格全票(不一致即兜底)", new_strict), ("多数>50%", new_majority)):
    print("  [%s]" % name)
    print("    %-9s %9s %9s %9s" % ("阶元", "改变", "变空", "由空变有"))
    for r in TAX:
        ch = em = fl = 0
        for cid in all_ids:
            o, n = old_wide[cid][r], new[cid][r]
            if o != n:
                ch += 1
                if n == MISSING:
                    em += 1
                if o == MISSING:
                    fl += 1
        print("    %-9s %9d %9d %9d" % (r, ch, em, fl))

# ---- 第 4 步兜底死代码核验 ----
score_best = {}
for cid in all_ids:
    best, bs = None, -1
    for t in tools:
        row = data[t].get(cid)
        if not row:
            continue
        s = sum(1 for r in TAX if row[r])
        if s > bs:
            bs, best = s, row
    score_best[cid] = best

touched = 0
sample_touched = []
for cid in all_ids:
    for r in TAX:
        if not cells[cid][r]:                      # 票首为空
            b = score_best[cid]
            if b and b[r]:                          # 最完整工具行确实有值
                touched += 1
                if len(sample_touched) < 5:
                    sample_touched.append((cid, r, b[r]))
print("\n=== 4) 「最完整工具行」兜底(:537-546) 实际生效格数 ===")
print("  票首为空且最完整工具行有值的格子: %d" % touched)
for s in sample_touched:
    print("    %s %s -> %s" % s)

# ---- 样例：Mimiviridae 行 ----
shipped = os.path.join(BASE, "final_integrated_classification.tsv")
fam_mimi = []
if os.path.exists(shipped):
    with open(shipped, newline="", encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if (row.get("Family") or "").strip() == "Mimiviridae":
                fam_mimi.append(row["contig_id"])
print("\n=== 5) 现产物里 Family=Mimiviridae 的 contig，在新规则下属/种怎么变 ===")
print("  共 %d 条" % len(fam_mimi))
for cid in fam_mimi[:10]:
    o = old_wide.get(cid, {})
    s = new_strict.get(cid, {})
    m = new_majority.get(cid, {})
    print("  %s" % cid[:58])
    print("    现行票首   Genus=%-14s Species=%s" % (o.get("Genus"), o.get("Species")))
    print("    严格全票   Genus=%-14s Species=%s" % (s.get("Genus"), s.get("Species")))
    print("    多数50%%    Genus=%-14s Species=%s" % (m.get("Genus"), m.get("Species")))
