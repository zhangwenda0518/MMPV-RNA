"""实测: v2 属级黑名单作用于 ensemble Plant.classified.fasta(20716) 的效果。

关注两点:
  1. 剔除多少条
  2. 剔除的是不是全是真垃圾 (属在黑名单, 且该属确非植物)
  3. 有没有误杀真植物病毒 (属是真植物属却因科兜底被拦)
"""
import ast
import csv
import collections

RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
src = open(f"{RUN}/run_host_prediction.py").read()
tree = ast.parse(src)

set_assigns = [n for n in tree.body if isinstance(n, ast.Assign)
               and any(isinstance(t, ast.Name) and
                       t.id in ("NON_PLANT_GENERA", "NON_PLANT_FAMILIES_FALLBACK",
                                "TRUSTED_LEVELS") for t in n.targets)]
funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef)
         and n.name in ("is_trusted_level", "is_blacklisted")]
ns = {}
exec(compile(ast.Module(body=set_assigns + funcs, type_ignores=[]), "<m>", "exec"), ns)
is_blacklisted = ns["is_blacklisted"]
NG = ns["NON_PLANT_GENERA"]

# 权威库属清单
auth = set()
with open("/tmp/auth_plant_genera.txt", errors="replace") as f:
    for l in f:
        l = l.strip()
        if l:
            auth.add(l)

# C9 全量: contig -> (genus, family, level, conf)
info = {}
with open("/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result/"
          "classification_result.tsv", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        info[row["contig_id"]] = (
            row.get("Genus", "").strip(), row.get("Family", "").strip(),
            row.get("Determination_Level", "").strip(),
            row.get("Confidence_Level", "").strip())

# ensemble Plant.fasta
fa = []
with open("/home/zhangwenda/data-test/out/06_HostPrediction/host_classified_fasta/"
          "Plant.classified.fasta") as f:
    for l in f:
        if l.startswith(">"):
            fa.append(l[1:].split()[0])

killed = collections.Counter()
kept = collections.Counter()
mimi_kept = 0
kill_in_auth = collections.Counter()   # 被删的属恰在权威库中 = 疑似误杀
for cid in fa:
    g, fm, lv, cl = info.get(cid, ("?", "?", "?", "?"))
    if is_blacklisted(fm, g, lv):
        killed[(fm, g, cl)] += 1
        if g in auth and g != "NA":
            kill_in_auth[(fm, g)] += 1
    else:
        kept[(fm, g)] += 1
        if fm == "Mimiviridae":
            mimi_kept += 1

print(f"ensemble Plant.fasta 总数: {len(fa)}")
print(f"被黑名单剔除: {sum(killed.values())} ({sum(killed.values())/len(fa)*100:.2f}%)")
print(f"保留: {len(fa)-sum(killed.values())}")
print(f"\nMimiviridae 名下保留: {mimi_kept} 条")
print(f"\n=== 剔除明细 (科/属/置信度) ===")
for k, v in killed.most_common(30):
    print(f"  {k[0]:<18} {k[1]:<20} {k[2]:<22} {v}")
print(f"\n=== 疑似误杀 (被删属在权威库中) ===")
if kill_in_auth:
    for k, v in kill_in_auth.most_common(20):
        print(f"  {k[0]:<18} {k[1]:<20} {v}")
else:
    print("  无")
