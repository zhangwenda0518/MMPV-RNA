"""生成黑名单剔除清单: 从 ensemble Plant.classified.fasta 中识别需删除的 contig。

输出:
  /tmp/blacklist_kill_ids.txt    待删 contig_id 清单 (每行一个)
  /tmp/blacklist_kill_detail.tsv 明细 (contig_id / genus / family / conf / reason)
"""
import ast
import csv

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

info = {}
with open("/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result/"
          "classification_result.tsv", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        info[row["contig_id"]] = (
            row.get("Genus", "").strip(), row.get("Family", "").strip(),
            row.get("Determination_Level", "").strip(),
            row.get("Confidence_Level", "").strip())

fa = []
with open("/home/zhangwenda/data-test/out/06_HostPrediction/host_classified_fasta/"
          "Plant.classified.fasta") as f:
    for l in f:
        if l.startswith(">"):
            fa.append(l[1:].split()[0])

kills = []
with open("/tmp/blacklist_kill_ids.txt", "w") as fo, \
     open("/tmp/blacklist_kill_detail.tsv", "w") as fd:
    fd.write("contig_id\tfamily\tgenus\tconfidence\tdetermination\n")
    for cid in fa:
        g, fm, lv, cl = info.get(cid, ("?", "?", "?", "?"))
        if is_blacklisted(fm, g, lv):
            kills.append(cid)
            fo.write(cid + "\n")
            fd.write(f"{cid}\t{fm}\t{g}\t{cl}\t{lv}\n")

print(f"待删清单: {len(kills)} 条")
print(f"  → /tmp/blacklist_kill_ids.txt")
print(f"  → /tmp/blacklist_kill_detail.tsv")
