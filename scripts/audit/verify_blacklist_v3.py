"""验证 v3: 用真实 Determination_Level 判定 (与管线同口径)。

同时对比: 662 旧清单 vs 新黑名单在当前数据的实际剔除范围。
"""
import ast
import csv
from collections import Counter

TARGET = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
PLANT_FA = ("/home/zhangwenda/data-test/out/06_HostPrediction/"
            "host_classified_fasta/Plant.classified.fasta")
CLS = ("/home/zhangwenda/data-test/out/06_HostPrediction/"
       "C9_ICTV_result/classification_result.tsv")

src = open(TARGET, encoding="utf-8").read()
t = ast.parse(src)
sets = [n for n in t.body if isinstance(n, ast.Assign)
        and any(isinstance(tg, ast.Name) and
                tg.id in ("NON_PLANT_GENERA", "NON_PLANT_FAMILIES_FALLBACK",
                          "TRUSTED_LEVELS") for tg in n.targets)]
fns = [n for n in t.body if isinstance(n, ast.FunctionDef)
       and n.name in ("is_trusted_level", "is_blacklisted")]
ns = {}
exec(compile(ast.Module(body=sets + fns, type_ignores=[]), "<m>", "exec"), ns)
is_blacklisted = ns["is_blacklisted"]

info = {}
with open(CLS, errors="replace") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        info[row["contig_id"]] = (
            (row.get("Genus", "") or "").strip(),
            (row.get("Family", "") or "").strip(),
            (row.get("Determination_Level", "") or "").strip())

kill = set(l.strip() for l in open("/tmp/blacklist_kill_ids.txt") if l.strip())

total = 0
hits = Counter()
in662 = 0
for l in open(PLANT_FA):
    if not l.startswith(">"):
        continue
    total += 1
    cid = l[1:].split()[0]
    if cid in kill:
        in662 += 1
    gen, fm, lv = info.get(cid, ("", "", ""))
    if is_blacklisted(fm, gen, lv):
        hits[(fm, gen, lv)] += 1

n = sum(hits.values())
print(f"Plant.classified.fasta: {total} 条")
print(f"其中命中旧 662 清单: {in662}")
print(f"v3 黑名单命中(真实Det_Level口径): {n}")
print("\n按科-属-判定层级:")
for (fm, gen, lv), c in hits.most_common(15):
    print(f"  {fm:<22} {gen:<24} {lv:<14} {c}")
