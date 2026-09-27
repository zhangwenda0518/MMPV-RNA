#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读预演：逐级共识 + 异议工具淘汰（大王 2026-08-30 口径）。

机制（严格照大王的说法）：
  从最粗阶元（Realm）开始逐级向下。
  每一级只看「当前还活着的工具」报的有效值：
    - 若活跃工具在这一级都没值（NA）-> 该级留空，active 不变（NA 不是反对票）
    - 若只有部分工具报值、且只有一个不同取值 -> 采用它，不淘汰任何人
      （大王：只有 ACVirus 有分类、其他都没分类，那就按 ACVirus）
    - 若出现多个不同取值 -> 取加权票首为共识；报了别的值的工具从 active 移除，
      后面更细的阶元不再看它们
  不引入 VMR 谱系约束（VMR 只用于事后评测"科属还打不打架"）。
  同时去掉「最完整工具行」兜底（已实测生效 0 格）。

两种淘汰门槛：
  P1 保守：只有票首占比 > 50%（有明确多数）时才淘汰异议者
  P2 激进：只要该级存在多个不同取值，就淘汰所有非票首工具
另测 P1b：被淘汰工具在"活跃工具全体无值"的阶元上可复活补值
"""
import csv, os, re
from collections import defaultdict

BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
        "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated")
VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"

TAX = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
RDW = {"Realm": 1, "Kingdom": 2, "Phylum": 4, "Class": 8,
       "Order": 16, "Family": 32, "Genus": 64, "Species": 128}
BIAS = {"ACVirus": 1.2, "VITAP": 1.1, "mmseqs": 1.0, "metabuli": 1.0,
        "CAT": 0.9, "genomad": 0.9, "diamond_lca": 0.8}
PRIORITY = ["ACVirus", "VITAP", "mmseqs", "metabuli", "CAT", "genomad", "diamond_lca"]
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
            if not cid or cid in out:
                continue
            row = {k: clean_rank(k, (r.get(k) or "").strip()) for k in TAX}
            for k in list(row):
                if row[k] and any(row[k].lower().endswith(s) for s in SUBRANK):
                    row[k] = None
            out[cid] = row
    return out


# ---------- VMR 谱系（仅用于评测） ----------
vmr_chain_by_species = defaultdict(set)
vmr_chain_by_genus = defaultdict(set)
vmr_chain_by_rank = defaultdict(set)
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t", quotechar='"')
    hdr = next(rd)
    idx = {n: hdr.index(n) for n in TAX}
    for row in rd:
        if len(row) <= idx["Species"]:
            continue
        chain = tuple((row[idx[r]] or "").strip() for r in TAX)
        if chain[7]:
            vmr_chain_by_species[chain[7].lower()].add(chain)
        if chain[6]:
            vmr_chain_by_genus[chain[6].lower()].add(chain)
        for i, rk in enumerate(TAX):
            if chain[i]:
                vmr_chain_by_rank[(rk, chain[i].lower())].add(chain)


def vmr_family_of(rank, value):
    chs = vmr_chain_by_rank.get((rank, value.lower()))
    if not chs:
        return None
    return {ch[5] for ch in chs if ch[5]}


# ---------- 载入 ----------
tools = [t for t in PRIORITY
         if os.path.exists(os.path.join(BASE, "standardized_%s.tsv" % t))]
data = {t: load_tool(os.path.join(BASE, "standardized_%s.tsv" % t)) for t in tools}
all_ids = set()
for t in tools:
    all_ids |= set(data[t])
print("工具 %d 个: %s" % (len(tools), ",".join(tools)))
print("contig 并集 %d" % len(all_ids))

# ---------- 权重（复刻自举） ----------
cnt = defaultdict(int)
for cid in all_ids:
    for t in tools:
        if cid in data[t]:
            cnt[cid] += 1
n_tools = len(tools)
common = [c for c in all_ids if cnt[c] == n_tools]
ids_w = common if common else [c for c in all_ids if cnt[c] >= max(2, n_tools // 2)]
boot = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
for cid in ids_w:
    for t in tools:
        row = data[t].get(cid)
        if not row:
            continue
        for r in TAX:
            if row[r]:
                boot[cid][r][row[r]] += 1
fs = {}
for cid in ids_w:
    for t in tools:
        row = data[t].get(cid)
        if not row:
            continue
        for r in TAX:
            if row[r]:
                fs.setdefault((cid, r, row[r]), len(fs))
bcons = {}
for cid in boot:
    for r in boot[cid]:
        bcons[(cid, r)] = sorted(boot[cid][r].items(),
                                 key=lambda kv: (-kv[1], fs.get((cid, r, kv[0]), 0)))[0][0]
agree = defaultdict(lambda: defaultdict(lambda: [0, 0]))
for cid in boot:
    for r in boot[cid]:
        c = bcons[(cid, r)]
        for t in tools:
            row = data[t].get(cid)
            if not row or not row[r]:
                continue
            agree[t][r][1] += 1
            if row[r].lower() == c.lower():
                agree[t][r][0] += 1
W = {t: {} for t in tools}
for t in tools:
    for r in TAX:
        w = RDW[r] * BIAS.get(t, 0.8)
        a, b = agree[t][r]
        if b and a / b > 0:
            w *= a / b
        W[t][r] = w
print("自举样本 %d 条 (全工具共有 %d)" % (len(ids_w), len(common)))

cells = defaultdict(lambda: defaultdict(list))
for cid in all_ids:
    for t in tools:
        row = data[t].get(cid)
        if not row:
            continue
        for r in TAX:
            if row[r]:
                w = W[t][r] * (species_quality(row[r]) if r == "Species" else 1.0)
                cells[cid][r].append((row[r], w, t))

# ---------- 现行基线 ----------
OLD = {}
for cid in all_ids:
    OLD[cid] = {}
    for r in TAX:
        recs = cells[cid][r]
        if not recs:
            continue
        agg, fir = defaultdict(float), {}
        for i, (v, w, t) in enumerate(recs):
            agg[v] += w
            fir.setdefault(v, i)
        OLD[cid][r] = sorted(agg.items(), key=lambda kv: (-kv[1], fir[kv[0]]))[0][0]


def topsum(pool):
    agg, fir = defaultdict(float), {}
    for i, (v, w, t) in enumerate(pool):
        agg[v] += w
        fir.setdefault(v, i)
    order = sorted(agg.items(), key=lambda kv: (-kv[1], fir[kv[0]]))
    return order[0][0], agg[order[0][0]] / sum(agg.values()), len(order)


# ---------- P1 / P2 / P1b ----------
RES = {"P1": {}, "P2": {}, "P1b": {}}
kill_log = {r: defaultdict(int) for r in TAX}     # 在哪个阶元被淘汰
solo = defaultdict(int)                            # 单工具独有值的格数
for cid in all_ids:
    for tag in RES:
        RES[tag][cid] = {}
    for tag, mode in (("P1", "majority"), ("P2", "any"), ("P1b", "revive")):
        active = set(tools)
        out = RES[tag][cid]
        for r in TAX:
            recs = [(v, w, t) for (v, w, t) in cells[cid][r] if t in active]
            if not recs:
                if mode == "revive":
                    pool = cells[cid][r]
                    if pool:
                        out[r] = topsum(pool)[0]
                continue
            seen = {v for v, w, t in recs}
            if len(seen) == 1:
                out[r] = recs[0][0]
                continue
            top, ratio, nd = topsum(recs)
            out[r] = top
            if mode == "majority" and ratio <= 0.5:
                continue
            if tag == "P1b":
                if ratio <= 0.5:
                    continue
            kill = sorted({t for v, w, t in recs if v != top},
                          key=lambda t: PRIORITY.index(t))
            if kill:
                kill_log[r][kill[0]] += 1
            active -= set(kill)
    # 单工具独有值统计
    for r in TAX:
        recs = cells[cid][r]
        if len(recs) == 1:
            solo[r] += 1

print("\n=== 各阶元「只有 1 个工具报值」的格数（大王：这种要按该工具取值） ===")
for r in TAX:
    tot = len([1 for cid in all_ids if cells[cid][r]])
    print("  %-9s %6d / %6d  (%s)" % (r, solo[r], tot,
                                      "%.1f%%" % (100.0 * solo[r] / tot) if tot else "n/a"))

print("\n=== 淘汰事件：在哪个阶元上首次被淘汰的工具 ===")
for r in TAX:
    if kill_log[r]:
        print("  %-9s %s" % (r, dict(sorted(kill_log[r].items(), key=lambda kv: -kv[1]))))


def metrics(tag, tab):
    f_sp = f_ge = n_sp = n_ge = 0
    for cid in all_ids:
        row = tab[cid]
        fam, sp, ge = row.get("Family"), row.get("Species"), row.get("Genus")
        if fam and sp:
            fs_ = vmr_family_of("Species", sp)
            if fs_:
                n_sp += 1
                if fam not in fs_:
                    f_sp += 1
        if fam and ge:
            fg_ = vmr_family_of("Genus", ge)
            if fg_:
                n_ge += 1
                if fam not in fg_:
                    f_ge += 1
    print("  %-6s Species侧科不一致 %4d /%6d (%.1f%%)   Genus侧科属不一致 %4d /%6d (%.1f%%)"
          % (tag, f_sp, n_sp, 100.0 * f_sp / n_sp if n_sp else 0,
             f_ge, n_ge, 100.0 * f_ge / n_ge if n_ge else 0))


print("\n=== 科属/科种不匹配（VMR 谱系校验） ===")
metrics("OLD", OLD)
for tag in ("P1", "P2", "P1b"):
    metrics(tag, RES[tag])

print("\n=== 与现行逐格对比 ===")
for tag in ("P1", "P2", "P1b"):
    tab = RES[tag]
    print("  [%s]" % tag)
    print("    %-9s %8s %8s %8s" % ("阶元", "改变", "变空", "由空变有"))
    for r in TAX:
        ch = em = fl = 0
        for cid in all_ids:
            o, n = OLD[cid].get(r), tab[cid].get(r)
            if (o or None) != (n or None):
                ch += 1
                if o and not n:
                    em += 1
                if n and not o:
                    fl += 1
        print("    %-9s %8d %8d %8d" % (r, ch, em, fl))

# 活跃工具枯竭
print("\n=== 淘汰后的活跃工具数分布（P1，取最后一个有值阶元处） ===")
dist = defaultdict(int)
for cid in all_ids:
    alive = set(tools)
    for r in TAX:
        recs = [(v, w, t) for (v, w, t) in cells[cid][r] if t in alive]
        if not recs:
            continue
        seen = {v for v, w, t in recs}
        if len(seen) == 1:
            continue
        top, ratio, nd = topsum(recs)
        if ratio > 0.5:
            alive -= {t for v, w, t in recs if v != top}
    dist[len(alive)] += 1
for k in sorted(dist):
    print("  剩 %d 个工具: %d contig" % (k, dist[k]))

print("\n=== 实例：门层 ACVirus 与多数不同 → 被淘汰，后续阶元变化（P1） ===")
shown = 0
for cid in sorted(all_ids):
    recs = cells[cid].get("Phylum") or []
    if len({v for v, w, t in recs}) < 2:
        continue
    top = topsum(recs)[0]
    if not any(v != top and t == "ACVirus" for v, w, t in recs):
        continue
    print("  %s" % cid[:56])
    for r in TAX:
        o = OLD[cid].get(r)
        n = RES["P1"][cid].get(r)
        if (o or None) != (n or None):
            print("    %-8s 旧 %-22s -> P1 %s" % (r, o or "<空>", n or "<空>"))
    shown += 1
    if shown >= 6:
        break
print("  共找到 %d 条" % sum(1 for cid in all_ids
                          if len({v for v, w, t in (cells[cid].get("Phylum") or [])}) > 1
                          and any(v != topsum(cells[cid]["Phylum"])[0] and t == "ACVirus"
                                  for v, w, t in cells[cid]["Phylum"])))
