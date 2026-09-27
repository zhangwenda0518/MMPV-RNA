#!/usr/bin/env python3
"""量化 R 脚本 KNOWN_REALMS 白名单 (L33 + L215) 造成的 Realm 静默清空。

R 脚本 L215: dt[!tolower(Realm) %in% tolower(KNOWN_REALMS), Realm := NA]
KNOWN_REALMS 只有 6 个旧名, 缺 ICTV 2025 新立的 4 个 realm (Monodnaviria 拆分而来):
  Efunaviria / Floreoviria / Pleomoviria / Volvereviria

用法:
  python3 probe_realm_whitelist.py --comb C.tsv --prod P.tsv [--label X]
输出:
  1. combined 中出现的全部有效 Realm 值 + 计数 + 是否在白名单
  2. 成品 Realm=NA 但 combined 里有非白名单 realm 的 contig 数 (可恢复格)
"""
import csv, argparse

OLD = {"riboviria", "monodnaviria", "duplodnaviria", "varidnaviria", "adnaviria", "ribozyviria"}
NEW2025 = {"efunaviria", "floreoviria", "pleomoviria", "volvereviria"}
NAV = {"", "NA", "N/A", "-", "no rank", "undefined", "unknown", "null", "unclassified"}
NAV_LC = {x.lower() for x in NAV}

ap = argparse.ArgumentParser()
ap.add_argument("--comb", required=True)
ap.add_argument("--prod", required=True)
ap.add_argument("--label", default="")
a = ap.parse_args()


def norm(v):
    return (v or "").strip().strip('"').strip()


inv, combrealm = {}, {}
with open(a.comb, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    cols = next(rd)
    ci = {c: i for i, c in enumerate(cols)}
    ri, ki = ci["Realm"], ci[cols[0]]
    for r in rd:
        if not r or len(r) <= ri:
            continue
        v = norm(r[ri])
        if not v or v.lower() in NAV_LC:
            continue
        inv[v] = inv.get(v, 0) + 1
        combrealm.setdefault(norm(r[ki]), set()).add(v)

prod_na, prod_tot = 0, 0
recover_detail = {}
with open(a.prod, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    cols = next(rd)
    ki, ri = cols.index("contig_id"), cols.index("Realm")
    for r in rd:
        if not r or len(r) <= ri:
            continue
        prod_tot += 1
        pv = norm(r[ri])
        if pv and pv.lower() not in NAV_LC:
            continue
        prod_na += 1
        for cv in combrealm.get(norm(r[ki]), ()):
            if cv.lower() not in OLD:
                recover_detail[cv] = recover_detail.get(cv, 0) + 1

print("== %s 成品 %d 行, Realm=NA %d 行" % (a.label, prod_tot, prod_na))
print("   combined 有效 Realm 清单:")
for v, n in sorted(inv.items(), key=lambda kv: -kv[1]):
    tag = "旧白名单" if v.lower() in OLD else ("ICTV2025新增(被清空)" if v.lower() in NEW2025 else "未知")
    print("     %-16s %7d  %s" % (v[:16], n, tag))
tot = sum(recover_detail.values())
print("   成品 Realm=NA 且 combined 里有非旧白名单 realm 的 (格/值): %d %s"
      % (tot, dict(sorted(recover_detail.items(), key=lambda kv: -kv[1])[:6])))
