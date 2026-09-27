"""排查: 植物科下的非植物属 (黑名单之外漏掉的)。

判定逻辑:
  1. 从 C9/classify 结果收集 (Family, Genus) 组合
  2. 用黑名单 is_blacklisted 标记已知非植物属
  3. 对每个"植物科", 列出其下所有属, 标出:
     - 已在黑名单 -> 已处理
     - 未在黑名单但形态可疑 -> 待核查
  4. 用权威植物病毒属库 (/tmp/auth_plant_genera.txt) 交叉验证

输出重点: 植物病毒科中, 属不在权威植物库、也不在黑名单的组合。
"""
import ast
import csv
import os
from collections import defaultdict

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
NON_PLANT_GENERA = ns.get("NON_PLANT_GENERA", set())
NON_PLANT_FAMILIES_FALLBACK = ns.get("NON_PLANT_FAMILIES_FALLBACK", set())
print(f"黑名单: 属 {len(NON_PLANT_GENERA)} 个, 科兜底 {len(NON_PLANT_FAMILIES_FALLBACK)} 个")

# 权威植物病毒属库
auth = set()
if os.path.exists("/tmp/auth_plant_genera.txt"):
    for l in open("/tmp/auth_plant_genera.txt", errors="replace"):
        l = l.strip()
        if l:
            auth.add(l)
print(f"权威植物属库: {len(auth)} 条")

# 数据源: C9 classify 全量
CLS = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result/classification_result.tsv"
fam_genus = defaultdict(lambda: defaultdict(int))
fam_total = defaultdict(int)
with open(CLS, errors="replace") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        fam = (row.get("Family", "") or "").strip()
        gen = (row.get("Genus", "") or "").strip()
        if not fam or fam == "NA":
            continue
        fam_total[fam] += 1
        fam_genus[fam][gen or "NA"] += 1

print(f"\n数据源: {CLS}")
print(f"科数: {len(fam_total)}\n")

# ── 判定什么是"植物科" ──
# 用 Plant 判定结果: C9 里 Final_Host=Plant 的科
# 这里改用 ensemble 结果里的 Plant.classified.fasta 反查科
import collections
plant_fams = collections.Counter()
info = {}
with open(CLS, errors="replace") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        info[row["contig_id"]] = (row.get("Family", "").strip(), row.get("Genus", "").strip())
with open("/home/zhangwenda/data-test/out/06_HostPrediction/host_classified_fasta/"
          "Plant.classified.fasta") as f:
    for l in f:
        if l.startswith(">"):
            cid = l[1:].split()[0]
            fm = info.get(cid, ("", ""))[0]
            if fm:
                plant_fams[fm] += 1

print(f"=== Plant.classified.fasta 中的科 (共 {len(plant_fams)} 个) ===")
for fm, n in plant_fams.most_common():
    print(f"  {fm:<28} {n}")

# ── 核心: 植物科下的属, 找非植物属 ──
print("\n" + "=" * 70)
print("=== 排查: 植物科下, 非植物属 (不在黑名单 + 不在权威植物库) ===")
print("=" * 70)

suspect = []
for fm in plant_fams:
    for gen, n in fam_genus.get(fm, {}).items():
        in_bl = is_blacklisted(fm, gen, "Family(via Family)")
        in_auth = gen in auth
        if not in_bl and not in_auth and gen != "NA":
            suspect.append((fm, gen, n))
        elif gen == "NA" and not is_blacklisted(fm, "NA", "Family(via Family)"):
            suspect.append((fm, gen, n))

suspect.sort(key=lambda x: -x[2])
print(f"\n可疑组合数: {len(suspect)}\n")
print(f"{'科':<24} {'属':<24} {'条数':>6}")
print("-" * 60)
for fm, gen, n in suspect[:60]:
    print(f"{fm:<24} {gen:<24} {n:>6}")
