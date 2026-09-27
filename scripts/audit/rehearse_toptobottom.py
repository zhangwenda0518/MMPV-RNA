#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读预演 B/C：把「从粗到细」纳入共识/兜底，量化"科属不匹配"能否消掉。

大王口径（2026-08-30）：每个级别全都应该是共识排序；都不一致或都没有时，
最后兜底走「从粗到细」的顺序（Realm -> Species），以免引起科属不匹配；
同时去掉「最完整工具行」兜底（已实测生效 0 格）。

本脚本对比四个方案（全部只读，不写任何产物）：
  OLD  : 现行 = 每阶元独立取加权票首（复刻 build_consensus 第 3 步）
  A    : 共识(全票一致|多数>50%) + 不一致/无值→按工具优先级取票
  B    : 现行票首 + 从粗到细谱系相容过滤（最小改动版）
  C    : 共识 + 不一致/无值→从粗到细（父级 VMR 投影）兜底

指标：主指标 = final 表口径「Species 在 VMR 里的 Family 与行内 Family 不同」的行数；
      辅指标 = 「Genus 在 VMR 里的 Family 与行内 Family 不同」的行数（科属不匹配）。
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
            row = {}
            for k in TAX:
                row[k] = clean_rank(k, (r.get(k) or "").strip())
            for k in list(row):
                if row[k] and any(row[k].lower().endswith(s) for s in SUBRANK):
                    row[k] = None
            out[cid] = row
    return out


# ---------------- VMR 谱系 ----------------
vmr_chain_by_species = defaultdict(set)   # species.lower -> set(chain tuple)
vmr_chain_by_genus = defaultdict(set)
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t", quotechar='"')
    hdr = next(rd)
    idx = {n: hdr.index(n) for n in ("Realm", "Kingdom", "Phylum", "Class",
                                     "Order", "Family", "Genus", "Species")}
    for row in rd:
        if len(row) <= idx["Species"]:
            continue
        chain = tuple((row[idx[r]] or "").strip() for r in TAX)
        if not any(chain):
            continue
        if chain[7]:
            vmr_chain_by_species[chain[7].lower()].add(chain)
        if chain[6]:
            vmr_chain_by_genus[chain[6].lower()].add(chain)
print("VMR: species 键 %d, genus 键 %d" % (len(vmr_chain_by_species), len(vmr_chain_by_genus)))

# ---------------- 载入工具 ----------------
tools = [t for t in PRIORITY
         if os.path.exists(os.path.join(BASE, "standardized_%s.tsv" % t))]
data = {t: load_tool(os.path.join(BASE, "standardized_%s.tsv" % t)) for t in tools}
all_ids = set()
for t in tools:
    all_ids |= set(data[t])
print("工具 %d 个, contig 并集 %d" % (len(tools), len(all_ids)))

# ---------------- 权重（复刻 consensus_stats 自举） ----------------
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
print("自举样本 %d 条 (全工具共有 %d)" % (len(ids_w), len(common)))

W = {t: {} for t in tools}
for t in tools:
    for r in TAX:
        w = RDW[r] * BIAS.get(t, 0.8)
        a, b = agree[t][r]
        if b and a / b > 0:
            w *= a / b
        W[t][r] = w

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

# ---------------- 谱系相容工具 ----------------
def compatible(rank, value, prev):
    """value 在 rank 层是否与已定粗阶元 prev 相容。无法校验（不在 VMR）返回 True。"""
    if not prev:
        return True
    if rank == "Species":
        chains = vmr_chain_by_species.get(value.lower())
    elif rank == "Genus":
        chains = vmr_chain_by_genus.get(value.lower())
    else:
        chains = None
    if not chains:
        # 非属/种，或 VMR 无此值：用 VMR 全表反查该阶元值所属链，取不到就宽容
        chains = vmr_chain_by_rank.get((rank, value.lower()))
        if not chains:
            return True
    for ch in chains:
        ok = True
        for r2, v2 in prev.items():
            if ch[TAX.index(r2)] and ch[TAX.index(r2)].lower() != v2.lower():
                ok = False
                break
        if ok:
            return True
    return False


vmr_chain_by_rank = defaultdict(set)
for s in vmr_chain_by_species.values():
    for ch in s:
        for i, r in enumerate(TAX):
            if ch[i]:
                vmr_chain_by_rank[(r, ch[i].lower())].add(ch)


def children_of(rank, value):
    """VMR 中该阶元值下面所有子阶元取值集合。"""
    chs = vmr_chain_by_rank.get((rank, value.lower()))
    if not chs:
        return None
    i = TAX.index(rank)
    return {ch[i + 1] for ch in chs if i + 1 < 8 and ch[i + 1]}


def vmr_family_of(rank, value):
    chs = vmr_chain_by_rank.get((rank, value.lower()))
    if not chs:
        return None
    fams = {ch[5] for ch in chs if ch[5]}
    return fams


# ---------------- 四种求值 ----------------
OLD, A, B, C = {}, {}, {}, {}
for cid in all_ids:
    OLD[cid] = {}
    # OLD
    for r in TAX:
        recs = cells[cid][r]
        if not recs:
            continue
        agg = defaultdict(float)
        fir = {}
        for i, (v, w, t) in enumerate(recs):
            agg[v] += w
            fir.setdefault(v, i)
        OLD[cid][r] = sorted(agg.items(), key=lambda kv: (-kv[1], fir[kv[0]]))[0][0]

    def summarize(r, pool):
        agg = defaultdict(float)
        fir = {}
        for i, (v, w, t) in enumerate(pool):
            agg[v] += w
            fir.setdefault(v, i)
        order = sorted(agg.items(), key=lambda kv: (-kv[1], fir[kv[0]]))
        top = order[0][0]
        tot = sum(agg.values())
        return top, (agg[top] / tot), len(agg)

    def fallback_tool(r):
        for t in PRIORITY:
            for v, w, tt in cells[cid][r]:
                if tt == t:
                    return v
        return None

    # A: 共识 + 工具优先级兜底
    A[cid] = {}
    for r in TAX:
        if not cells[cid][r]:
            continue
        top, ratio, nd = summarize(r, cells[cid][r])
        A[cid][r] = top if (nd == 1 or ratio > 0.5) else fallback_tool(r)

    # B: 现行票首 + 从粗到细相容过滤（最小改动）
    B[cid] = {}
    prev = {}
    for r in TAX:
        recs = cells[cid][r]
        if not recs:
            continue
        top, ratio, nd = summarize(r, recs)
        pool = [x for x in recs if compatible(r, x[0], prev)]
        if pool:
            v = summarize(r, pool)[0]
        else:
            v = None          # 全部与粗阶元冲突 -> 留空，不倒向越权值
        if v:
            B[cid][r] = v
            prev[r] = v
        else:
            prev[r] = top      # 保守：该阶元留空但用票首继续约束下游

    # C: 共识 + 从粗到细兜底
    C[cid] = {}
    prev = {}
    for r in TAX:
        recs = cells[cid][r]
        pool = [x for x in recs if compatible(r, x[0], prev)]
        v = None
        if pool:
            top, ratio, nd = summarize(r, pool)
            if nd == 1 or ratio > 0.5:
                v = top
        if v is None:
            # 从粗到细兜底：父级（最近的已定粗阶元）在 VMR 里的子阶元投影
            for pr in reversed(TAX[:TAX.index(r)]):
                if pr in prev:
                    kids = children_of(pr, prev[pr])
                    if kids:
                        if len(kids) == 1:
                            v = next(iter(kids))
                        elif pool:
                            inter = [x for x in pool if x[0] in kids]
                            if inter:
                                v = summarize(r, inter)[0]
                    break
        if v is None and pool:
            v = summarize(r, pool)[0]      # 最后仍无解时退回本阶元票首（与旧行为一致）
        if v:
            C[cid][r] = v
            prev[r] = v
        elif pool:
            prev[r] = summarize(r, pool)[0]


# ---------------- 指标 ----------------
def metrics(tag, tab):
    f_sp = f_ge = 0
    n_sp_known = n_ge_known = 0
    for cid in all_ids:
        row = tab.get(cid) or {}
        fam = row.get("Family")
        sp, ge = row.get("Species"), row.get("Genus")
        if fam and sp:
            fs_ = vmr_family_of("Species", sp)
            if fs_:
                n_sp_known += 1
                if fam not in fs_:
                    f_sp += 1
        if fam and ge:
            fg_ = vmr_family_of("Genus", ge)
            if fg_:
                n_ge_known += 1
                if fam not in fg_:
                    f_ge += 1
    print("  %-6s Species侧科不一致 %4d /%6d   Genus侧科属不一致 %4d /%6d"
          % (tag, f_sp, n_sp_known, f_ge, n_ge_known))
    return f_sp, f_ge


print("\n=== 科属/科种不匹配（VMR 谱系校验） ===")
print("  仅统计能在 VMR 里查到 Species / Genus 的行（新病毒无法校验，不计）")
mets = {}
mets["OLD"] = metrics("OLD", OLD)
mets["A"] = metrics("A", A)
mets["B"] = metrics("B", B)
mets["C"] = metrics("C", C)

print("\n=== 与现行产物逐格对比（仅密度<100%的阶元会变） ===")
print("  变空 = 有无 -> 无（信息损失）；消失的越权值 = 原值与已定粗阶元谱系冲突")
for tag, tab in (("A 共识+工具优先级", A), ("B 票首+从粗到细", B), ("C 共识+从粗到细", C)):
    print("  [%s]" % tag)
    print("    %-9s %8s %8s %8s" % ("阶元", "改变", "变空", "由空变有"))
    for r in TAX:
        ch = em = fl = 0
        for cid in all_ids:
            o = OLD[cid].get(r)
            n = tab[cid].get(r)
            if (o or None) != (n or None):
                ch += 1
                if o and not n:
                    em += 1
                if n and not o:
                    fl += 1
        print("    %-9s %8d %8d %8d" % (r, ch, em, fl))

print("\n=== 抽样：Family=Mimiviridae 且属/种越权（B/C 修正效果） ===")
shown = 0
for cid in sorted(all_ids):
    if (OLD[cid].get("Family") or "") != "Mimiviridae":
        continue
    og, os_ = OLD[cid].get("Genus"), OLD[cid].get("Species")
    if not og or og == os_:
        continue
    bad = vmr_family_of("Genus", og)
    if not bad or "Mimiviridae" in bad:
        continue
    print("  %s" % cid[:58])
    print("    旧  F=%s G=%s S=%s" % (OLD[cid].get("Family"), og, os_))
    print("    B   F=%s G=%s S=%s" % (B[cid].get("Family"), B[cid].get("Genus"), B[cid].get("Species")))
    print("    C   F=%s G=%s S=%s" % (C[cid].get("Family"), C[cid].get("Genus"), C[cid].get("Species")))
    shown += 1
    if shown >= 8:
        break

# 反向抽样：Genus 定得住、Family 是越权的
print("\n=== 抽样：Genus=Potyvirus 但 Family 非 Potyviridae ===")
shown = 0
for cid in sorted(all_ids):
    if (OLD[cid].get("Genus") or "") != "Potyvirus":
        continue
    f = OLD[cid].get("Family")
    if f == "Potyviridae":
        continue
    print("  %s" % cid[:58])
    print("    旧  F=%s G=%s S=%s" % (f, "Potyvirus", OLD[cid].get("Species")))
    print("    B   F=%s G=%s S=%s" % (B[cid].get("Family"), B[cid].get("Genus"), B[cid].get("Species")))
    print("    C   F=%s G=%s S=%s" % (C[cid].get("Family"), C[cid].get("Genus"), C[cid].get("Species")))
    shown += 1
    if shown >= 6:
        break
