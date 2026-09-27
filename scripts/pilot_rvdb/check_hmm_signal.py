"""检查 HMM 证据对 09b KEEP/REVIEW/DROP 的区分力。

数据源:
  A. 09b/virus_validation/rescue_evidence_scored.tsv  -> verdict (真值)
  B. 07_Checkv/Plant/completeness.tsv                 -> CheckV 的 HMM/AAI 列
     (hmm_num_hits, hmm_completeness_lower/upper, aai_num_hits, aai_completeness, kmer_freq)

问题: 现有裁决完全基于 CDD 域 + blast + 长度门槛, 未用 HMM。
      若把 HMM 加进来, 能不能把 REVIEW/DROP 里的"漏网病毒"捞回来?
"""
import csv
import statistics as st

B = "/home/zhangwenda/data-test/out"
S = f"{B}/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv"
C = f"{B}/07_Checkv/Plant/completeness.tsv"

verdict = {}
cdd_evid = {}
with open(S) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        verdict[r["contig_id"]] = r["verdict"]
        cdd_evid[r["contig_id"]] = r.get("cdd_evidence", "")
print(f"scored.tsv: {len(verdict)} 条")

cv = {}
with open(C) as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        cv[r["contig_id"]] = r
print(f"07 Plant completeness.tsv: {len(cv)} 条")

inter = set(verdict) & set(cv)
print(f"交集(scored ∩ 07Plant): {len(inter)} / scored {len(verdict)}\n")


def fnum(x):
    try:
        return float(x)
    except Exception:
        return None


print("=== 各 verdict 的 HMM / AAI 信号 ===")
hdr = f"{'verdict':<8}{'n':<7}{'有HMM命中%':<12}{'HMM命中中位':<12}{'有AAI命中%':<12}{'最大HMM完整度':<14}{'有kmer%':<10}"
print(hdr)
for v in ["KEEP", "REVIEW", "DROP"]:
    ids = [i for i in inter if verdict[i] == v]
    if not ids:
        print(f"{v:<8}0")
        continue
    hmm_pos = sum(1 for i in ids if (fnum(cv[i].get("hmm_num_hits")) or 0) > 0)
    hmm_n = [fnum(cv[i].get("hmm_num_hits")) or 0 for i in ids]
    aai_pos = sum(1 for i in ids if (fnum(cv[i].get("aai_num_hits")) or 0) > 0)
    hmm_c = [fnum(cv[i].get("hmm_completeness_upper")) or 0 for i in ids]
    kmer_pos = sum(1 for i in ids if (fnum(cv[i].get("kmer_freq")) or 0) > 0)
    n = len(ids)
    print(f"{v:<8}{n:<7}{100*hmm_pos/n:<12.1f}{st.median(hmm_n):<12.1f}"
          f"{100*aai_pos/n:<12.1f}{max(hmm_c):<14.1f}{100*kmer_pos/n:<10.1f}")

print("\n=== 关键: REVIEW 里被 CDD 判'非病毒域'的, HMM 有命中的有多少 ===")
for v in ["REVIEW", "DROP"]:
    ids = [i for i in inter if verdict[i] == v]
    novir = [i for i in ids if cdd_evid[i] not in ("PASS_CORE", "PASS_CORE_FAM", "PASS_VIRAL")]
    novir_hmm = [i for i in novir if (fnum(cv[i].get("hmm_num_hits")) or 0) > 0]
    print(f"  {v}: 非病毒域 {len(novir)} 条, 其中 HMM 有命中 {len(novir_hmm)} 条 "
          f"({100*len(novir_hmm)/len(novir):.1f}%)" if novir else f"  {v}: 无")

print("\n=== KEEP 里非病毒域的(DROP侧对照) ===")
ids = [i for i in inter if verdict[i] == "KEEP"]
novir = [i for i in ids if cdd_evid[i] not in ("PASS_CORE", "PASS_CORE_FAM", "PASS_VIRAL")]
print(f"  KEEP 总数 {len(ids)}, 非病毒域 {len(novir)}")

print("\n=== cdd_evidence × HMM 交叉 (交集内) ===")
from collections import Counter
cc = Counter()
for i in inter:
    e = cdd_evid[i] or "?"
    h = "HMM+" if (fnum(cv[i].get("hmm_num_hits")) or 0) > 0 else "HMM-"
    cc[(verdict[i], e, h)] += 1
for k in sorted(cc, key=lambda x: (x[0] != "KEEP", x[0], x[1], x[2])):
    print(f"  {k[0]:<7}{k[1]:<14}{k[2]:<6}{cc[k]}")
