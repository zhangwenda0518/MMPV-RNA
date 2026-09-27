#!/usr/bin/env python3
"""对成品对比中的 值->NA 丢失格给出机械判定, 不靠个案叙述。

三级判据(依次):
  A. 祖先链不相容: 丢失值在参照库中的任一较粗祖先(rank-1..Realm)与成品该 contig 的新值不同
     -> 由产物层「阶元相容」既有规则清空。
  B. 占位值取胜: 该阶元在新 combined 里取最多数/并列最多的是 environmental/unclassified 类占位串
     (is_valid_value_vec L181-184 不过滤多词占位串, 产物层 L918 才清空且不回落)
     -> 由既有投票 + 占位清理规则清空。
  C. 其余 -> 未解释, 强制打印。

用法:
  python3 explain_lost_cells.py --prod-old O.tsv --prod-new N.tsv [--comb-new C.tsv] [--label X] [--show]
"""
import csv, re, argparse

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PH = re.compile("environmental|uncultured|unclassified|unidentified|unassigned|not assigned|no hit|unknown|undefined", re.I)
NAV = {"", "NA", "N/A", "-", "no rank", "undefined", "unknown", "null", "default", "Unclassified"}
NAV_LC = {x.lower() for x in NAV}


def is_val(s):
    """同 R 脚本 is_valid_value_vec (L181-184): NA 类取值不算投票, 也不算异议。"""
    return bool(s) and s.lower() not in NAV_LC

ap = argparse.ArgumentParser()
ap.add_argument("--prod-old", required=True)
ap.add_argument("--prod-new", required=True)
ap.add_argument("--comb-new", default="")
ap.add_argument("--label", default="")
ap.add_argument("--show", action="store_true")
a = ap.parse_args()


def norm(v):
    return (v or "").strip().strip('"').strip()


def load(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t")
        return {norm(r.get("contig_id")): r for r in rd}


old, new = load(a.prod_old), load(a.prod_new)
lost = []
for k, ro in old.items():
    rn = new.get(k)
    if not rn:
        continue
    for rk in RANKS:
        vo, vn = norm(ro.get(rk)), norm(rn.get(rk))
        if vo not in NAV and vn in NAV:
            lost.append((k, rk, vo))

need = {v.lower() for _, _, v in lost}
ref = {}
if need:
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in need and nm not in ref:
                ref[nm] = parts

comb = {}
if a.comb_new:
    needk = {k for k, _, _ in lost}
    with open(a.comb_new, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        ccols = next(rd)
        for r in rd:
            if r and r[0].strip('"') in needk:
                comb.setdefault(r[0].strip('"'), []).append(r)
    cidx = {c: i for i, c in enumerate(ccols)}

nA = nA2 = nB = nB2 = nC = 0
print("== %s  值->NA %d 格" % (a.label, len(lost)))
for k, rk, v in lost:
    ri = RANKS.index(rk)
    parts = ref.get(v.lower())
    # A. 整条祖先链比对
    mismatch = None
    if parts:
        for anc in RANKS[:ri]:
            ref_anc = (parts[REF_IDX[anc]] or "").strip()
            new_anc = norm(new[k].get(anc))
            if ref_anc and new_anc and ref_anc != new_anc:
                mismatch = (anc, new_anc, anc, ref_anc)
                break
    if mismatch:
        verdict = "A 祖先链不相容 (新%s=%s, 该值参照%s=%s)" % mismatch
        nA += 1
    elif parts and not any((parts[REF_IDX[anc]] or "").strip() for anc in RANKS[:ri]) \
            and any(is_val(norm(new[k].get(anc))) for anc in RANKS[:ri]):
        # A2: 参照库该名一条较粗阶元都没记 (如 NCBI Pandoravirus 无科级), 无法机械核验
        verdict = "A2 参照库该名无祖先记录, 无法机械核验 (新%s=%s)" % (
            RANKS[ri - 1], norm(new[k].get(RANKS[ri - 1])) or "NA")
        nA2 += 1
    else:
        # B. 占位值取胜
        vals = {}
        if cidx.get(rk) is not None:
            for r in comb.get(k, []):
                tv = norm(r[cidx[rk]]) if len(r) > cidx[rk] else ""
                vals[tv] = vals.get(tv, 0) + 1
        top = sorted(vals.items(), key=lambda kv: -kv[1])
        topn = top[0][1] if top else 0
        ph = [x for x in top if PH.search(x[0]) and x[1] >= topn]
        dist = ", ".join("%s x%d" % (x[0][:34], x[1]) for x in top[:5]) or "无 combined"
        # C. 逐级淘汰: 该值唯一来源工具在更粗阶元报异议, 失去后续投票权
        prune = None
        if cidx.get("tool") is not None:
            for r in comb.get(k, []):
                if len(r) <= cidx[rk] or norm(r[cidx[rk]]) != v:
                    continue
                for anc in RANKS[:ri]:
                    ai = cidx.get(anc)
                    if ai is None or len(r) <= ai:
                        continue
                    tv, pv = norm(r[ai]), norm(new[k].get(anc))
                    if is_val(tv) and is_val(pv) and tv != pv:
                        prune = (norm(r[cidx["tool"]]), anc, tv, pv)
                        break
                if prune:
                    break
        if ph:
            verdict = "B 占位值多数/并列取胜后清空 (该阶元取值: %s)" % dist
            nB += 1
        elif prune:
            verdict = "C 逐级淘汰 (唯一来源工具 %s 在 %s 报异议 %s, 成品 %s)" % prune
            nB2 += 1
        else:
            verdict = "D !! 未解释 (该阶元取值: %s)" % dist
            nC += 1
    if a.show or verdict.startswith("D"):
        print("   %s %s %s: %s -> NA | %s" % ("D " if verdict.startswith("D") else "v ", k[:44], rk, v, verdict))
print("   汇总: A 祖先链不相容 %d + A2 参照库无祖先 %d + B 占位值取胜 %d + C 逐级淘汰 %d = 机械解释 %d 格 / D 未解释 %d 格"
      % (nA, nA2, nB, nB2, nA + nA2 + nB + nB2, nC))
