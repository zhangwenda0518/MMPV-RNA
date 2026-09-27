import csv
import collections
import glob
import os
import re

src = open("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py").read()

# 用 AST 解析集合字面量, 避免文本解析歧义
import ast
tree = ast.parse(src)
found = {}
for node in ast.walk(tree):
    if isinstance(node, ast.Assign):
        for t in node.targets:
            if isinstance(t, ast.Name) and t.id in (
                    "NON_PLANT_GENERA", "NON_PLANT_FAMILIES_FALLBACK"):
                found[t.id] = ast.literal_eval(node.value)

NON_PLANT_GENERA = found["NON_PLANT_GENERA"]
FALLBACK = found["NON_PLANT_FAMILIES_FALLBACK"]
TRUSTED = {"genus", "species", "species*"}

print(f"AST 加载: 属黑名单 {len(NON_PLANT_GENERA)} / 科兜底 {len(FALLBACK)}")
assert "Biavirus" in NON_PLANT_GENERA, "Biavirus 缺失!"
assert "Tobaliviridae" in FALLBACK, "Tobaliviridae 缺失!"
assert "Fabavirus" not in NON_PLANT_GENERA, "Fabavirus 误入黑名单!"

def is_blacklisted(fam, gen, det):
    if det is None:
        return False
    if str(det).split("(")[0].strip().lower() in TRUSTED:
        return False
    gen = (gen or "").strip()
    fam = (fam or "").strip()
    if gen and gen not in ("NA", "nan"):
        return gen in NON_PLANT_GENERA
    return bool(fam and fam in FALLBACK)

base_dir = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"
killed = collections.Counter()
kept_mimi = 0
plant_total = 0
for tsv in glob.glob(os.path.join(base_dir, "*.classified.tsv")):
    with open(tsv, errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if (row.get("Predicted_Host", "") or "").strip() != "Plant":
                continue
            plant_total += 1
            fam = (row.get("Family", "") or "").strip()
            gen = (row.get("Genus", "") or "").strip()
            if is_blacklisted(fam, gen, row.get("Determination_Level")):
                killed[(fam, gen)] += 1
            elif fam == "Mimiviridae":
                kept_mimi += 1

tot_k = sum(killed.values())
print(f"\nC9 判 Plant 总数: {plant_total}")
print(f"被否决: {tot_k} ({tot_k/plant_total*100:.2f}%) | 保留: {plant_total-tot_k}")
print(f"\nMimiviridae 判 Plant 保留数: {kept_mimi} (期望 244 = 全部保留)")
print(f"\n=== 被否决明细 top20 ===")
for (fm, g), n in killed.most_common(20):
    print(f"  {fm:<22} {g:<22} {n}")
