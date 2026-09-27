#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐阶元投票实录：把 cascade 在单个 contig 上的每一步还原出来，并全量验证。

语义（对照补丁 v6.7 的 build_cascade_winners）：
    Realm -> Species 逐级递降；每级只在「仍活跃」的工具之间按权重计票；
    该级活跃工具都没报值 -> 留空，活跃集不变（无值不是反对票）；
    票首达到淘汰阈值且票首唯一 -> 采用票首，并把该级报了其它取值的工具移出活跃集；
    否则 -> 采用票首，不淘汰任何人；被淘汰者不复活；
    平票 -> 先看「下一（更细）级的一致性」（同值工具数 / 支持工具数），仍分不出则按 TIE_ORDER 兜底。

RULE 开关：gt = v6.6（share > 门槛）；ge = v6.7（票首唯一 且 share >= 门槛）。

先做全量对账：仿真结果必须与 <OUT>/final_integrated_classification.tsv 逐格一致，
一致才可以用它解释「科属为什么不匹配」。
"""
import csv
import glob
import os
from collections import Counter, OrderedDict

OUT = os.environ.get("DIAG_OUT", "/tmp/pregate_cascade")
LEVELS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
SHARE = 0.5
# 淘汰判据开关：gt = v6.6（share > 门槛）；ge = v6.7（票首唯一 且 share >= 门槛）
RULE = os.environ.get("MMPV_RULE", "gt")


def triggers(share, n_top):
    if RULE == "ge":
        return share >= SHARE and n_top == 1
    return share > SHARE


EMPTY = {"", "na", "n/a", "nan", "-", "none", "null"}
# 与补丁 TIE_BREAK_ORDER 同口径：平票兜底顺序
TIE_ORDER = ["ACVirus", "CAT", "VITAP", "diamond_lca", "genomad", "metabuli", "mmseqs"]


def norm(x):
    x = (x or "").strip().strip('"').strip()
    return "" if x.lower() in EMPTY else x


# ---------- 权重 ----------
w = {}
tools_order = []
with open(os.path.join(OUT, "tool_weights.tsv"), newline="", encoding="utf-8") as f:
    rd = csv.DictReader(f, delimiter="\t")
    ranks_in_file = [c for c in rd.fieldnames if c != "Tool"]
    for r in rd:
        tool = norm(r["Tool"])
        tools_order.append(tool)
        for lv in ranks_in_file:
            w[(tool, lv)] = float(r.get(lv) or 0)
print("权重表列: %s" % ", ".join(ranks_in_file))
print("工具顺序（补丁里由 TIE_BREAK_ORDER 显式排定）: %s" % ", ".join(tools_order))
# 平票兜底序 = TIE_BREAK_ORDER（未列入的工具按权重表原序追加）
RANK = {t: i for i, t in enumerate([t for t in TIE_ORDER if t in tools_order]
                                   + [t for t in tools_order if t not in TIE_ORDER])}
print("平票兜底序: %s%s" % (", ".join(sorted(RANK, key=RANK.get)),
                            "" if list(sorted(RANK, key=RANK.get)) == tools_order else "  （与权重表行序不同！）"))


# ---------- 平票裁决（与补丁同口径）----------
def tie_scores(tops, votes, per_tool, lv):
    """平票候选打分：下一级同值工具数 / 支持工具数；次键为工具兜底序。"""
    nxt = LEVELS[LEVELS.index(lv) + 1] if lv != LEVELS[-1] else None
    rows = []
    for v in tops:
        ts = votes[v]
        nv = [x for x in (per_tool.get(t, {}).get(nxt, "") for t in ts) if x] if nxt else []
        modal = max(Counter(nv).values()) if nv else 0
        rows.append(((modal / len(ts) if ts else 0.0), min(RANK[t] for t in ts), v))
    rows.sort(key=lambda r: (-r[0], r[1]))
    return rows


def resolve(votes, per_tool, lv):
    """返回 (票首值, 平票归属)：None 非平票 / "next" 下一级一致性定 / "order" 工具顺序兜底。"""
    wsum = {v: sum(w[(t, lv)] for t in ts) for v, ts in votes.items()}
    top = max(wsum.values())
    tops = [v for v, s in wsum.items() if abs(s - top) < 1e-9]
    if len(tops) == 1:
        return tops[0], None
    rows = tie_scores(tops, votes, per_tool, lv)
    kind = "next" if (rows[0][0] > 0 and (len(rows) == 1 or rows[0][0] > rows[1][0])) else "order"
    return rows[0][2], kind

# ---------- 各工具逐阶元取值 ----------
vals = {}          # contig -> tool -> {rank: value}
for p in sorted(glob.glob(os.path.join(OUT, "standardized_*.tsv"))):
    tool = os.path.basename(p)[len("standardized_"):-len(".tsv")]
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = norm(r["contig_id"])
            d = vals.setdefault(cid, {})
            d[tool] = {lv: norm(r.get(lv)) for lv in LEVELS}
print("contig 数 %d" % len(vals))

# ---------- 仿真 ----------
sim = {}
traces = {}
n_cell_prune = 0      # 有异议被淘汰的格次（与补丁日志「淘汰格次」同口径）
n_event_prune = 0     # (contig, 工具) 淘汰事件（与补丁日志「淘汰工具投票权」同口径）
for cid, per_tool in vals.items():
    active = list(tools_order)
    out = {}
    tr = []
    for lv in LEVELS:
        votes = OrderedDict()
        for t in active:
            v = per_tool.get(t, {}).get(lv, "")
            if v:
                votes.setdefault(v, []).append(t)
        if not votes:
            out[lv] = ""
            tr.append((lv, {}, None, 0.0, "该级活跃工具全无值 -> 留空，活跃集不变", []))
            continue
        tot = sum(sum(w[(t, lv)] for t in ts) for ts in votes.values())
        # v6.6 口径：先算票首；平票则先看下一级一致性，再按 TIE_ORDER 兜底
        best_v, tie_kind = resolve(votes, per_tool, lv)
        best = (best_v, votes[best_v])
        bw = sum(w[(t, lv)] for t in best[1])
        share = bw / tot if tot > 0 else 0.0
        out[lv] = best[0]
        killed = []
        n_top = len([v for v in votes if abs(sum(w[(t, lv)] for t in votes[v]) - max(sum(w[(t, lv)] for t in votes[x]) for x in votes)) < 1e-9])
        if triggers(share, n_top):
            killed = [t for v, ts in votes.items() if v != best[0] for t in ts]
            active = [t for t in active if t not in killed]
            n_cell_prune += 1
            n_event_prune += len(killed)
            why = "票首占比 %.1f%% 达阈值且票首唯一 -> 采用并淘汰异议工具 %s" % (share * 100, ",".join(killed) or "∅")
        else:
            why = "票首占比 %.1f%% 未达阈值或票首并列 -> 采用但不淘汰（避免误杀）" % (share * 100)
        if tie_kind == "next":
            why = "平票(下一级一致性定) " + why
        elif tie_kind == "order":
            why = "平票(工具顺序兜底) " + why
        tr.append((lv, {v: sorted(ts) for v, ts in votes.items()}, best[0], share, why, killed))
    sim[cid] = out
    traces[cid] = tr

# ---------- 全量对账 ----------
actual = {}
with open(os.path.join(OUT, "final_integrated_classification.tsv"), newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        actual[norm(r["contig_id"])] = {lv: norm(r.get(lv)) for lv in LEVELS}

diff = Counter()
for cid in actual:
    s = sim.get(cid, {})
    for lv in LEVELS:
        if s.get(lv, "") != actual[cid][lv]:
            diff[lv] += 1
print()
print("仿真 vs 实际（OUT=%s）逐格差异：" % OUT)
if not diff:
    print("  全部 8 个阶元 0 差异，仿真与 R 实现一致")
else:
    for lv in LEVELS:
        print("  %-8s %6d" % (lv, diff[lv]))

print("\n淘汰判据 RULE=%s（gt=v6.6 share>门槛；ge=v6.7 唯一票首且 share>=门槛）" % RULE)
print("对照补丁日志：本仿真 淘汰格次 %d / 淘汰工具投票权 %d 次" % (n_cell_prune, n_event_prune))
stat = Counter()
per_rank = Counter(); per_rank_elim = Counter(); per_rank_tie = Counter()
tie_win_order = Counter(); tie_win_next = Counter()
for cid, tr in traces.items():
    for lv, votes, win, share, why, killed in tr:
        if not votes:
            stat["该级无任何活跃工具报值（留空）"] += 1
            continue
        per_rank[lv] += 1
        top = max(sum(w[(t, lv)] for t in ts) for ts in votes.values())
        n_top = sum(1 for ts in votes.values() if abs(sum(w[(t, lv)] for t in ts) - top) < 1e-9)
        if n_top > 1:
            stat["票首平票（多值同分）"] += 1
            per_rank_tie[lv] += 1
            _, tie_kind = resolve(votes, vals[cid], lv)
            stat["  平票由下一级一致性定"] += 1 if tie_kind == "next" else 0
            stat["  平票由工具顺序兜底"] += 1 if tie_kind == "order" else 0
            win_tool = min(votes[win], key=lambda t: RANK[t])
            (tie_win_next if tie_kind == "next" else tie_win_order)[win_tool] += 1
        if killed:
            stat["票首达阈值且唯一 -> 采用并淘汰异议工具"] += 1
            per_rank_elim[lv] += 1
        else:
            stat["票首未达阈值或并列 -> 采用但不淘汰任何人"] += 1
print()
print("投票决定全局统计（%d 个 contig × 8 阶元）" % len(traces))
for k, v in stat.most_common():
    print("  %-36s %7d" % (k, v))
print("  按阶元 计票格次 / 不淘汰格次 / 平票格次")
for lv in LEVELS:
    print("    %-9s %6d %6d %6d" % (lv, per_rank[lv], per_rank[lv] - per_rank_elim[lv], per_rank_tie[lv]))
print("  平票格次归属：下一级一致性定 %d / 工具顺序兜底 %d"
      % (stat["  平票由下一级一致性定"], stat["  平票由工具顺序兜底"]))
print()
print("平票赢家的值由哪个工具提供（按票首值归属到兜底序最早的支持工具）")
print("  %-12s %8s %8s %8s" % ("工具", "顺序定", "下一级定", "合计"))
for t in sorted(tie_win_order.keys() | tie_win_next.keys(), key=lambda t: RANK[t]):
    a, b = tie_win_order[t], tie_win_next[t]
    print("  %-12s %8d %8d %8d" % (t, a, b, a + b))
print("  %-12s %8d %8d %8d" % ("合计", sum(tie_win_order.values()), sum(tie_win_next.values()),
                                sum(tie_win_order.values()) + sum(tie_win_next.values())))

# ---------- 单 contig 实录 ----------
def dump(cid, title):
    print()
    print("=" * 78)
    print("[%s] %s" % (title, cid))
    print("=" * 78)
    print("  %-9s %-46s %-11s %s" % ("阶元", "活跃工具的取值(权重)", "票首占比", "决定"))
    for lv, votes, win, share, why, killed in traces.get(cid, []):
        s = " ".join("%s[%s]" % (v, "+".join("%.1f" % w[(t, lv)] for t in ts))
                     for v, ts in sorted(votes.items(), key=lambda kv: -sum(w[(t, lv)] for t in kv[1])))
        print("  %-9s %-46s %-11s %s" % (lv, s[:46], ("%.1f%%" % (share * 100)) if votes else "-", why))
    for lv in ("Family", "Genus"):
        print("  实际落定 %-7s %s" % (lv, actual.get(cid, {}).get(lv, "")))

# ---------- 199 行分解 ----------
mism = []
fam_gen_pairs = 0
for cid in vals:
    per_tool = vals[cid]
    fam = actual.get(cid, {}).get("Family", "")
    gen = actual.get(cid, {}).get("Genus", "")
    if not fam or not gen:
        continue
    fam_gen_pairs += 1
    fam_tools = [t for t in tools_order if per_tool.get(t, {}).get("Family", "") == fam]
    gen_tools = [t for t in tools_order if per_tool.get(t, {}).get("Genus", "") == gen]
    if set(fam_tools) & set(gen_tools):
        continue
    mism.append(cid)

print()
print("科属都非空的行数（分母） %d；其中工具层不自洽 %d（%.2f%%）"
      % (fam_gen_pairs, len(mism), 100.0 * len(mism) / fam_gen_pairs if fam_gen_pairs else 0.0))

# 这 199 行的属，是 cascade 投票定的，还是 cascade 之后被「由种提属」改的？
after = Counter()
for cid in mism:
    sg, ag = sim[cid].get("Genus", ""), actual[cid]["Genus"]
    sp = actual[cid].get("Species", "")
    first = sp.split(" ")[0] if " " in sp else ""
    if sg == ag:
        after["A 属就是 cascade 投票结果：不匹配是投票本身造成的"] += 1
    else:
        after["B 属在 cascade 之后被改写"] += 1
        if first and first.lower() == ag.lower():
            after["  其中 属 == Species 首词（由种提属补出来的）"] += 1
        else:
            after["  其中 属另有来源"] += 1
for k, v in after.most_common():
    print("  %-52s %5d" % (k, v))
buckets = Counter()
gen_elim = Counter()
for cid in mism:
    tr = {lv: item for lv, item in zip(LEVELS, traces[cid])}
    lv, votes, win, share, why, killed = tr["Family"]
    fam = actual[cid]["Family"]
    g_item = tr["Genus"]
    gen_win = g_item[2]
    gen_tools = [t for v, ts in g_item[1].items() if v == gen_win for t in ts]
    killed_before = set()
    for l in LEVELS[:LEVELS.index("Genus")]:
        killed_before |= set(tr[l][5])
    if not votes:
        buckets["A 科级活跃工具全无值(科是更早阶元或填充来的)"] += 1
    elif triggers(share, len([v for v in votes if abs(sum(w[(t, lv)] for t in votes[v]) - max(sum(w[(t, lv)] for t in votes[x]) for x in votes)) < 1e-9])):
        buckets["B 科级票首达阈值且唯一并淘汰了异议工具"] += 1
    else:
        buckets["C 科级票首未达阈值或并列 -> 不淘汰任何工具，属级重新投票"] += 1
    if any(t in killed_before for t in gen_tools):
        gen_elim["属级赢家中有工具在更粗阶元已被淘汰(不该发生)"] += 1
    else:
        gen_elim["属级赢家全部在属级仍是活跃工具(符合语义)"] += 1
    # 属级赢家在科级报了什么
    if gen_tools:
        fv = {vals[cid].get(t, {}).get("Family", "") for t in gen_tools}
        if fv == {fam}:
            gen_elim["属级赢家在科级报的就是共识科"] += 1
        elif "" in fv and len(fv) == 1:
            gen_elim["属级赢家在科级沉默"] += 1
        else:
            gen_elim["属级赢家在科级报了别的科/部分沉默"] += 1
    else:
        gen_elim["属级赢家没有工具支持(属是填充或由种提属来的)"] += 1

# ---------- 单 contig 实录 ----------
for pref, title in (("CRR1440126_clean_NODE_914", "各阶元赢家不同工具"),
                    ("CRR1126135_clean_NODE_29_length", "沉默逃逸"),
                    ("CRR1126135_clean_NODE_499", "属无任何工具支持")):
    hit = [c for c in vals if c.startswith(pref)]
    if hit:
        dump(hit[0], title)

print()
for k, v in buckets.most_common():
    print("  %-46s %5d" % (k, v))
for k, v in gen_elim.most_common():
    print("  %-46s %5d" % (k, v))
