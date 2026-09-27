"""核对黑名单清单条数 + is_blacklisted 门控原文 + plant_virus_db 白名单路径是否存在."""
import ast
import os

RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
src = open(RUN, encoding="utf-8", errors="replace").read()
t = ast.parse(src)
keep = [n for n in t.body if isinstance(
    n, (ast.Import, ast.ImportFrom, ast.Assign, ast.FunctionDef, ast.AnnAssign))]
ns = {"__file__": RUN}
exec(compile(ast.Module(body=keep, type_ignores=[]), "<rhp>", "exec"), ns)

G = ns["NON_PLANT_GENERA"]
F = ns["NON_PLANT_FAMILIES_FALLBACK"]
print("NON_PLANT_GENERA =", len(G), "条")
print("NON_PLANT_FAMILIES_FALLBACK =", len(F), "条")
print("TRUSTED_LEVELS =", sorted(ns["TRUSTED_LEVELS"]))
print("PHAGE_CLASSES =", sorted(ns["PHAGE_CLASSES"]))
print("RVH_MAP 键数 =", len(ns.get("RVH_MAP", {})))
for fn in ("is_trusted_level", "is_blacklisted", "normalize_c9"):
    node = next(n for n in t.body if isinstance(n, ast.FunctionDef) and n.name == fn)
    print("\n--- %s ---" % fn)
    print(ast.get_source_segment(src, node).rstrip())

v3 = ["Betapartitivirus", "Gammapartitivirus", "Cryspovirus"]
print("\n[v3 补入的三个代表是否在属名单]", [(x, x in G) for x in v3])

MISS_G = ["Alphacedratvirus", "Bacilladnavirus", "Bacillarnavirus", "Diatodnavirus", "Dinornavirus",
          "Eimeriavirus", "Fadolivirus", "Heliosvirus", "Keisodnavirus", "Leishmaniavirus",
          "Losannavirus", "Marseillevirus", "Mazarbulvirus", "Medusavirus", "Mimivirus",
          "Mimoreovirus", "Moumouvirus", "Nimphelosvirus", "Oceanusvirus", "Pandoravirus",
          "Phaeovirus", "Phialvirus", "Pithovirus", "Protobacilladnavirus", "Raphidovirus",
          "Shilevirus", "Sogarnavirus", "Sputnikvirus", "Stormycovirus", "Styxivirus",
          "Tetrivirus", "Yaravirus"]
MISS_F = ["Alvernaviridae", "Bacilladnaviridae", "Circoviridae", "Leishbuviridae", "Marnaviridae",
          "Ouroboviridae", "Phypoliviridae", "Pseudototiviridae", "Sputniviroviridae", "Yaraviridae"]
print("[仍不在属名单的漏登属] %d/%d:" % (sum(x not in G for x in MISS_G), len(MISS_G)),
      [x for x in MISS_G if x not in G])
print("[仍不在科名单的漏登科] %d/%d:" % (sum(x not in F for x in MISS_F), len(MISS_F)),
      [x for x in MISS_F if x not in F])

print("\n=== plant_virus_db 白名单路径核对 ===")
cands = [
    "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv",
    "/home/zhangwenda/MMPV-RNA/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv",
    "/home/zhangwenda/database/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv",
]
for p in cands:
    if os.path.exists(p):
        with open(p, errors="replace") as f:
            n = sum(1 for _ in f)
        print("  存在:", p, "| 行数(含表头):", n, "| 字节:", os.path.getsize(p))
        with open(p, errors="replace") as f:
            print("  表头:", f.readline().rstrip()[:160])
    else:
        print("  不存在:", p)
d = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify"
if os.path.isdir(d):
    print("  上级目录", d, "内容:", sorted(os.listdir(d))[:25])
    dc = os.path.join(d, "classified_clean")
    if os.path.isdir(dc):
        print("  classified_clean 内容:")
        for fn in sorted(os.listdir(dc))[:25]:
            fp = os.path.join(dc, fn)
            print("    %-34s %s" % (fn, os.path.getsize(fp) if os.path.isfile(fp) else "<dir>"))

print("\n=== C7 概率表行数 ===")
C7 = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
for fn in ["order_host_probability.tsv", "family_host_probability.tsv",
           "genus_host_probability.tsv", "species_host_probability.tsv"]:
    fp = os.path.join(C7, fn)
    if os.path.exists(fp):
        with open(fp, errors="replace") as f:
            n = sum(1 for _ in f) - 1
        print("  %-34s %d 行" % (fn, n))
    else:
        print("  %-34s 不存在" % fn)
