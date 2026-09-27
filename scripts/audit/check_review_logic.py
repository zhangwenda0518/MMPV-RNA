"""复盘 09b 现有 REVIEW 判定的确切构成。

integrate_rescue_evidence.py 的 REVIEW 有三条不同触发路径, 都写进了 fp_reason。
直接从 rescue_evidence_scored.tsv 里统计, 不用推断。
"""
import csv
import statistics as st

S = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv"

rows = []
with open(S) as f:
    rd = csv.DictReader(f, delimiter="\t")
    cols = rd.fieldnames
    for r in rd:
        rows.append(r)
print("scored.tsv 列:", cols)
print("总条数:", len(rows))

rev = [r for r in rows if r["verdict"] == "REVIEW"]
print(f"\nREVIEW 总数: {len(rev)}\n")

print("=== REVIEW 按 fp_reason 归类 ===")
from collections import Counter
c = Counter()
for r in rev:
    fr = (r.get("fp_reason") or "").strip()
    key = fr.split(":")[0].strip() if fr else "(空)"
    c[key] += 1
for k, v in c.most_common():
    print(f"  {v:<6} {k}")

print("\n=== REVIEW 按 cdd_evidence × fp_reason ===")
c2 = Counter()
for r in rev:
    fr = (r.get("fp_reason") or "").strip()
    key = fr.split(":")[0].strip() if fr else "(空)"
    c2[(r.get("cdd_evidence", ""), key)] += 1
for k in sorted(c2, key=lambda x: (x[0], x[1])):
    print(f"  {k[0]:<14}{k[1]:<28}{c2[k]}")

print("\n=== 三条触发路径的关键量 ===")


def fnum(x, d=0.0):
    try:
        return float(x)
    except Exception:
        return d


# 路径1: PASS_VIRAL 但过短
p1 = [r for r in rev if r.get("cdd_evidence") == "PASS_VIRAL"]
# 路径2/3: 无病毒域
p23 = [r for r in rev if r.get("cdd_evidence") not in ("PASS_CORE", "PASS_CORE_FAM", "PASS_VIRAL")]
p2 = [r for r in p23 if (r.get("fp_reason") or "").startswith("no_viral_domain_but_blast")]
p3 = [r for r in p23 if (r.get("fp_reason") or "").startswith("short_fragment")]
print(f"  路径1 PASS_VIRAL 但过短      : {len(p1)}")
if p1:
    ratios = [fnum(r["length"]) / fnum(r["genus_avg_len"]) for r in p1 if fnum(r["genus_avg_len"]) > 0]
    print(f"        长度/属均长 中位 {st.median(ratios):.3f}" if ratios else "")
print(f"  路径2 无病毒域 + 够长 + blast : {len(p2)}")
if p2:
    aa = [fnum(r["aa_pident"]) for r in p2 if fnum(r["aa_pident"]) > 0]
    nt = [fnum(r["nt_pident"]) for r in p2 if fnum(r["nt_pident"]) > 0]
    print(f"        aa_pident 中位 {st.median(aa):.1f}" if aa else "        无 aa 命中")
    print(f"        nt_pident 中位 {st.median(nt):.1f}" if nt else "        无 nt 命中")
print(f"  路径3 无病毒域 + 过短         : {len(p3)}")
print(f"  三路径合计: {len(p1)+len(p2)+len(p3)}")
