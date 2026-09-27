"""从 R 端落盘的 tool_weights.tsv 反推每个工具各阶元的 rate。

base(Tool,Rank) = RANK_DEPTH_WEIGHTS[Rank] * TOOL_BIAS[Tool]
weight = base * max(rate, RATE_FLOOR)   (若有 rate)，否则 weight = base
implied = weight / base
  implied == RATE_FLOOR(0.25) -> rate <= 0.25（被下限截断，v6.4 修复对象：rate 极低或为 0）
  implied == 1                -> rate == 1 或该阶元无 rate 记录（两者无法区分，需另证）
  其它                        -> rate == implied
"""
import sys

RANK_DEPTH_WEIGHTS = {"Realm": 1, "Kingdom": 2, "Phylum": 4, "Class": 8,
                      "Order": 16, "Family": 32, "Genus": 64, "Species": 128}
TOOL_BIAS = {"ACVirus": 1.2, "VITAP": 1.1, "mmseqs": 1.0, "metabuli": 1.0, "CAT": 0.9,
             "genomad": 0.9, "diamond_lca": 0.8, "vcontact3": 0.7, "contigtax": 0.6,
             "BASTA": 0.6, "PhaGCN3": 0.8}
FLOOR = 0.25

for path in sys.argv[1:]:
    with open(path, encoding="utf-8") as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    head, body = rows[0], rows[1:]
    print("=" * 78)
    print("FILE %s" % path)
    print("  %-12s %s" % ("Tool", "  ".join("%9s" % h for h in head[1:])))
    floor_hits = []
    for r in body:
        tool = r[0]
        bias = TOOL_BIAS.get(tool, 0.8)
        cells = []
        for h, v in zip(head[1:], r[1:]):
            base = RANK_DEPTH_WEIGHTS[h] * bias
            imp = float(v) / base
            if abs(imp - FLOOR) < 1e-6:
                cells.append("%9s" % "<=%.2f" % FLOOR)
                floor_hits.append((tool, h))
            else:
                cells.append("%9.3f" % imp)
        print("  %-12s %s" % (tool, "  ".join(cells)))
    print("  [下限生效格] %d 个: %s" % (len(floor_hits), ", ".join("%s/%s" % t for t in floor_hits)))
