import csv
import collections
import glob
import os

# 检查 C9 里判 Plant 的 Family 字段实际值 (是否有空格/大小写问题)
base_dir = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"
fam_vals = collections.Counter()
mimi_examples = []
for tsv in glob.glob(os.path.join(base_dir, "*.classified.tsv")):
    with open(tsv, errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if (row.get("Predicted_Host", "") or "").strip() != "Plant":
                continue
            fam = row.get("Family", "")
            fam_vals[fam] += 1
            if "Mimi" in str(fam) or "Tobali" in str(fam):
                if len(mimi_examples) < 5:
                    mimi_examples.append((repr(fam), row.get("Genus"), row.get("Determination_Level")))

print("=== 判 Plant 的 Family 字段原始值 (含引号) top20 ===")
for k, v in fam_vals.most_common(20):
    print(f"  {k!r}: {v}")

print("\n=== Mamiviridae/Tobaliviridae 样例 ===")
for e in mimi_examples:
    print(f"  Family={e[0]} Genus={e[1]!r} Level={e[2]!r}")

# 检查源码里黑名单解析结果
src = open("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py").read()
def extract_set(name):
    key = name + " = "
    i = src.index(key) + len(key)
    j = src.index("}", i)
    body = src[i:j]
    return {tok.strip().strip("'\"") for tok in body.split(",") if tok.strip()}

S = extract_set("NON_PLANT_FAMILIES")
print(f"\n解析出 {len(S)} 个科, 其中含 Mimiviridae? {'Mimiviridae' in S}")
print(f"含 Tobaliviridae? {'Tobaliviridae' in S}")
# 看有没有脏 token
dirty = [t for t in S if not (t.endswith("viridae") or t.endswith("virinae"))]
print(f"非 viridae 结尾的 token ({len(dirty)}): {sorted(dirty)[:20]}")
