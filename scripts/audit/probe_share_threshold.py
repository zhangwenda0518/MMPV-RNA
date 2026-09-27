#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""量「唯一票首但不过半」的格次，并测算两种放宽门槛的代价。

背景：cascade 的淘汰判据是「票首占比 > CASCADE_MIN_SHARE(0.5) 才淘汰异议工具」。
7 家工具等权时票首占比只能落在 1/7 的格点上，3 票领跑（42.9%）或并列（<=50%）都踢不了人。
本脚本对每个 (contig, 阶元) 计票格次分类：

    E 淘汰类   share > 0.5                     -> 现状会踢人
    T 并列类   票首有 n_top >= 2 个值同分       -> 永不 > 0.5，走两级平票裁决
    U 卡线类   唯一票首且 share <= 0.5          -> 票势其实明确，却被门槛拦下

并跑三种规则各一遍全量仿真，比较最终值变化（同一套引擎，只换淘汰判据）：
    A 现状      share > 0.50
    B 降门槛    share > 0.40
    C 唯一票首  票首唯一即定音并淘汰异议（并列仍走两级裁决）

用法：DIAG_OUT=/tmp/rc66_count python3 probe_share_threshold.py
"""
import csv
import glob
import os
from collections import Counter, OrderedDict

OUT = os.environ.get("DIAG_OUT", "/tmp/rc66_count")
LEVELS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
EMPTY = {"", "na", "n/a", "nan", "-", "none", "null"}
TIE_ORDER = ["ACVirus", "CAT", "VITAP", "diamond_lca", "genomad", "metabuli", "mmseqs"]


def norm(x):
    x = (x or "").strip().strip('"').strip()
    return "" if x.lower() in EMPTY else x


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
RANK = {t: i for i, t in enumerate([t for t in TIE_ORDER if t in tools_order]
                                   + [t for t in tools_order if t not in TIE_ORDER])}


def tie_scores(tops, votes, per_tool, lv):
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
    wsum = {v: sum(w[(t, lv)] for t in ts) for v, ts in votes.items()}
    top = max(wsum.values())
    tops = [v for v, s in wsum.items() if abs(s - top) < 1e-9]
    if len(tops) == 1:
        return tops[0], None
    rows = tie_scores(tops, votes, per_tool, lv)
    kind = "next" if (rows[0][0] > 0 and (len(rows) == 1 or rows[0][0] > rows[1][0])) else "order"
    return rows[0][2], kind


vals = {}
for p in sorted(glob.glob(os.path.join(OUT, "standardized_*.tsv"))):
    tool = os.path.basename(p)[len("standardized_"):-len(".tsv")]
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = norm(r["contig_id"])
            vals.setdefault(cid, {})[tool] = {lv: norm(r.get(lv)) for lv in LEVELS}
print("输入 %s：contig %d，工具 %s" % (OUT, len(vals), ",".join(tools_order)))


def simulate(rule, min_share=0.5):
    """返回 (逐格结果, 淘汰格次, 淘汰工具事件数)。"""
    sim = {}
    n_cell = 0
    n_ev = 0
    for cid, per_tool in vals.items():
        active = list(tools_order)
        out = {}
        for lv in LEVELS:
            votes = OrderedDict()
            for t in active:
                v = per_tool.get(t, {}).get(lv, "")
                if v:
                    votes.setdefault(v, []).append(t)
            if not votes:
                out[lv] = ""
                continue
            tot = sum(sum(w[(t, lv)] for t in ts) for ts in votes.values())
            best_v, _ = resolve(votes, per_tool, lv)
            bw = sum(w[(t, lv)] for t in votes[best_v])
            share = bw / tot if tot > 0 else 0.0
            tops = [v for v in votes if abs(sum(w[(t, lv)] for t in votes[v]) - bw) < 1e-9]
            out[lv] = best_v
            if rule == "gshare":
                hit = share > min_share
            elif rule == "unique":
                hit = len(tops) == 1
            elif rule == "unique_ge":
                hit = len(tops) == 1 and share >= min_share
            else:
                hit = share > min_share and len(tops) == 1
            if hit:
                killed = [t for v, ts in votes.items() if v != best_v for t in ts]
                if killed:
                    n_cell += 1
                    n_ev += len(killed)
                active = [t for t in active if t not in killed]
        sim[cid] = out
    return sim, n_cell, n_ev


# ---------- 1. 计票格次分类（基线口径）----------
cls = Counter()
per_rank = Counter()
per_rank_cls = {lv: Counter() for lv in LEVELS}
patterns_u = Counter()
patterns_u_rank = Counter()
u_cases = []
n_prune = 0
for cid, per_tool in vals.items():
    active = list(tools_order)
    for lv in LEVELS:
        votes = OrderedDict()
        for t in active:
            v = per_tool.get(t, {}).get(lv, "")
            if v:
                votes.setdefault(v, []).append(t)
        if not votes:
            continue
        tot = sum(sum(w[(t, lv)] for t in ts) for ts in votes.values())
        best_v, _ = resolve(votes, per_tool, lv)
        bw = sum(w[(t, lv)] for t in votes[best_v])
        share = bw / tot if tot > 0 else 0.0
        tops = [v for v in votes if abs(sum(w[(t, lv)] for t in votes[v]) - bw) < 1e-9]
        per_rank[lv] += 1
        dissent = [t for v, ts in votes.items() if v != best_v for t in ts]
        bool_trigger = share > 0.5 and bool(dissent)
        if len(tops) > 1:
            per_rank_cls[lv]["T 并列票首"] += 1
            cls["T 并列票首（n_top>=2，数学上永不 > 0.5，走两级平票裁决）"] += 1
        elif share > 0.5 and not dissent:
            per_rank_cls[lv]["E0 全票一致"] += 1
            cls["E0 全票一致（share=100%，无人可踢）"] += 1
        elif share > 0.5:
            per_rank_cls[lv]["E1 触发淘汰"] += 1
            cls["E1 过半且有异议（现状真的踢人）"] += 1
        else:
            per_rank_cls[lv]["U 唯一票首但不过半"] += 1
            cls["U 唯一票首但不过半（现状不淘汰）"] += 1
            if abs(share - 0.5) < 1e-9:
                cls["   其中 share 恰好 50%（唯一最高但卡线）"] += 1
                patterns_u_rank[(lv, "精确卡线 50%")] += 1
            else:
                cls["   其中 share 严格小于 50%"] += 1
            n = sum(len(ts) for ts in votes.values())
            pat = ":".join(str(x) for x in sorted((len(ts) for ts in votes.values()), reverse=True))
            patterns_u["%s (有值工具 %d)" % (pat, n)] += 1
            patterns_u_rank[(lv, pat)] += 1
            u_cases.append((cid, lv, pat, share, best_v))
        if bool_trigger:
            n_prune += len(dissent)
            active = [t for t in active if t not in set(dissent)]

tot_cells = sum(per_rank.values())
print()
print("计票格次总计 %d（对照补丁日志「阶元计票 144771 格次」）" % tot_cells)
for k, v in cls.most_common():
    print("  %-46s %7d  %5.2f%%" % (k, v, 100.0 * v / tot_cells))
n_need = cls["E1 过半且有异议（现状真的踢人）"] + cls["T 并列票首（n_top>=2，数学上永不 > 0.5，走两级平票裁决）"] + cls["U 唯一票首但不过半（现状不淘汰）"]
print("  合计计票 %d；其中「真需要裁决」的格次（有异议）%d：过半踢人 %d + 并列 %d + 唯一票首卡线 %d"
      % (tot_cells, n_need, cls["E1 过半且有异议（现状真的踢人）"],
         cls["T 并列票首（n_top>=2，数学上永不 > 0.5，走两级平票裁决）"],
         cls["U 唯一票首但不过半（现状不淘汰）"]))
print("  工具淘汰事件合计 %d（对照补丁日志「淘汰工具投票权 8523 格次」；单位是 contig×工具 对，不是格次）" % n_prune)
print("  按阶元 计票 / 全票一致 / 过半踢人 / 并列 / 唯一票首卡线")
for lv in LEVELS:
    c = per_rank_cls[lv]
    print("    %-9s %6d %6d %6d %6d %6d" % (lv, per_rank[lv], c["E0 全票一致"], c["E1 触发淘汰"],
                                            c["T 并列票首"], c["U 唯一票首但不过半"]))

nU = sum(1 for _ in u_cases)
print()
print("U 类票型分布（唯一票首但不过半，共 %d 格次）" % nU)
for k, v in patterns_u.most_common(20):
    print("  %-28s %6d  %5.1f%%" % (k, v, 100.0 * v / max(nU, 1)))
print()
print("U 类按阶元 × 票型（前 20）")
for (lv, pat), v in patterns_u_rank.most_common(20):
    print("  %-9s %-24s %6d" % (lv, pat, v))

# ---------- 2. 三种规则的全量仿真对照 ----------
base, cell_a, ev_a = simulate("gshare", 0.5)
alt_b, _, _ = simulate("gshare", 0.4)
alt_c, _, _ = simulate("unique", 0.5)
alt_d, _, _ = simulate("gshare_uniq", 0.4)
alt_e, cell_e, ev_e = simulate("unique_ge", 0.5)
print()
print("淘汰判据的格次级代价（与补丁日志同口径：淘汰格次 / 淘汰工具投票权）")
print("  A 现状 share > 0.5         淘汰格次 %5d, 工具投票权 %5d  (对照 v6.6 日志 7519 / 8523)" % (cell_a, ev_a))
print("  E 唯一票首且 share >= 0.5  淘汰格次 %5d, 工具投票权 %5d  (对照 v6.7 日志 7698 / 8995)" % (cell_e, ev_e))


def cmp_sim(x, y, tag):
    per = Counter()
    contig_hit = set()
    for cid in y:
        for lv in LEVELS:
            if x.get(cid, {}).get(lv, "") != y[cid][lv]:
                per[lv] += 1
                contig_hit.add(cid)
    print("  %-26s 变动格次 %s | 涉及 contig %d" %
          (tag, " ".join("%s=%d" % (lv, per[lv]) for lv in LEVELS), len(contig_hit)))
    return per, contig_hit


print()
print("规则放宽后的最终值变动（相对 A 现状；Genus/Species 含复刻未跟的后置闸门，仅参考）")
cmp_sim(base, alt_b, "B 门槛降到 0.40")
cmp_sim(base, alt_c, "C 唯一票首即定音")
cmp_sim(base, alt_d, "D 门槛 0.40 且票首唯一（并列不踢）")
cmp_sim(base, alt_e, "E 唯一票首且 share>=0.5（只救 694 卡线格次）")

# ---------- 3. 抽样 ----------
print()
print("U 类抽样（前 10 条，看票势到底明不明确）")
for cid, lv, pat, share, win in u_cases[:10]:
    print("  %-42s %-9s %-12s 票首=%-24s %.1f%%" % (cid[:42], lv, pat, win[:24], share * 100))
