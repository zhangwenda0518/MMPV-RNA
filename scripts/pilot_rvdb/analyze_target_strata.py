#!/usr/bin/env python3
"""路径②内部按 blast 强度分层, 看 HMM 命中是否向强 blast 一侧收敛。"""
import csv
from collections import defaultdict

V3 = "/home/zhangwenda/rvdb_hmm_pilot/rvdb_hmm_per_contig_v3.tsv"
S = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv"

blast = {}
with open(S) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        blast[r["contig_id"]] = (r.get("aa_pident", ""), r.get("aa_qcov", ""),
                                 r.get("nt_pident", ""), r.get("nt_qcovs", ""), r.get("aa_species", ""))

rows = []
with open(V3) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(r)


def fn(x):
    try:
        return float(x)
    except Exception:
        return 0.0


def band(ap, aq, np_, nq):
    # blast 强度分层: 蛋白优先, 退化到核酸
    if ap >= 85 and aq >= 85:
        return "aa>=85 (已知种)"
    if ap >= 40 and aq >= 50:
        return "aa 40-85 (新种候选)"
    if np_ >= 85 and nq >= 85:
        return "nt>=85 (核酸强)"
    if np_ >= 40 and nq >= 50:
        return "nt 40-85 (核酸弱)"
    return "其他"


tgt = [r for r in rows if r["bucket"] == "no_viral_domain_but_blast"]
print(f"路径② 共 {len(tgt)}\n")

g = defaultdict(lambda: defaultdict(int))
for r in tgt:
    ap, aq, np_, nq, sp = blast.get(r["contig_id"], ("", "", "", "", ""))
    b = band(fn(ap), fn(aq), fn(np_), fn(nq))
    gg = g[b]
    gg["n"] += 1
    ev = fn(r["best_evalue"]) if r["best_evalue"] else None
    cov = fn(r["best_cov"]) if r["best_cov"] else 0.0
    if ev is not None and ev <= 1e-5:
        gg["e5"] += 1
        if cov >= 0.5:
            gg["e5c50"] += 1
    if r["n_hits"] and int(r["n_hits"]) > 0:
        gg["any"] += 1

print(f"{'blast 分层':<24}{'n':<7}{'有HMM命中%':<12}{'E5%':<8}{'E5&cov50%':<11}")
for b in ["aa>=85 (已知种)", "aa 40-85 (新种候选)", "nt>=85 (核酸强)", "nt 40-85 (核酸弱)", "其他"]:
    gg = g.get(b)
    if not gg:
        continue
    n = gg["n"]
    print(f"{b:<24}{n:<7}{100*gg['any']/n:<12.1f}{100*gg['e5']/n:<8.1f}{100*gg['e5c50']/n:<11.1f}")

# 对照: DROP 里也看 aa_pident 分布(应为 0)
print("\n=== 对照 DROP 的 blast 强度(应几乎无命中) ===")
dr = [r for r in rows if r["bucket"] == "DROP"]
g2 = defaultdict(int)
for r in dr[:2000]:
    ap, aq, np_, nq, sp = blast.get(r["contig_id"], ("", "", "", "", ""))
    g2[band(fn(ap), fn(aq), fn(np_), fn(nq))] += 1
for k, v in g2.items():
    print(f"  {k}: {v}")
