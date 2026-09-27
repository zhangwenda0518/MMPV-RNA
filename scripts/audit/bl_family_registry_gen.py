# -*- coding: utf-8 -*-
"""生成"非植物科登记清单": 植物库零记录 + 在结果表中出现的科, 逐科带计数。

产物:
  /tmp/non_plant_family_registry.tsv  逐科证据 (科/各表行数/代表属/植物库记录数/是否已登记)
  /tmp/non_plant_family_registry.txt  同上, 便于 sticky 阅读
"""
import ast, csv, os, sys
from collections import Counter, defaultdict

PLANT_DB = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
RUN_HOST = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
TABLES = [
    ("goji05", "/tmp/v62d_goji/final_integrated_classification.tsv"),
    ("onekp05", "/tmp/v62d_onekp/final_integrated_classification.tsv"),
    ("gojiPL", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
               "/10_Reports/All_plant.viruses_info.tsv"),
    ("onekpPL", "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus"
                "/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"),
]
RANK_SUFFIX = ("viricota", "viricotina", "viricetes", "viricetidae", "virales",
               "virineae", "virinae", "viridae")


def norm(x):
    v = (x or "").strip()
    return "" if v.lower() in ("", "na", "nan", "none", "null", "-") else v


def split_lineage(s):
    toks = [t.strip() for t in str(s or "").split(";") if t.strip()]
    for i, t in enumerate(toks):
        if t.endswith("viridae") and not t.endswith("virinae"):
            gen = ""
            if i + 1 < len(toks) and not toks[i + 1].endswith(RANK_SUFFIX):
                gen = toks[i + 1]
            return t, gen
    return "", ""


def col_of(cols, *keys):
    for c in cols:
        if all(k in c.lower() for k in keys):
            return c
    return None


def read_tsv(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t")
        return rd.fieldnames or [], list(rd)


def main():
    tree = ast.parse(open(RUN_HOST, encoding="utf-8", errors="replace").read())
    registered = set()
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "NON_PLANT_FAMILIES_FALLBACK" for t in node.targets):
            registered = set(ast.literal_eval(node.value))

    _, prows = read_tsv(PLANT_DB)
    plant_fam = Counter()
    for r in prows:
        f, _ = split_lineage(r.get("Virus_lineage"))
        if f:
            plant_fam[f] += 1

    fam_cnt = defaultdict(Counter)
    fam_gen = defaultdict(Counter)
    for label, path in TABLES:
        if not os.path.exists(path):
            print(f"[skip] {label}", file=sys.stderr)
            continue
        cols, rows = read_tsv(path)
        fc, gc = col_of(cols, "family"), col_of(cols, "genus")
        for r in rows:
            f, g = norm(r.get(fc)), norm(r.get(gc))
            if not f:
                continue
            fam_cnt[f][label] += 1
            if g:
                fam_gen[f][g] += 1

    cand = sorted(f for f in fam_cnt if f not in plant_fam and f not in registered)
    with open("/tmp/non_plant_family_registry.tsv", "w", encoding="utf-8", newline="") as fh:
        fh.write("Family\tgoji05\tonekp05\tgojiPL\tonekpPL\t植物库记录\t代表属(Top3)\n")
        for f in cand:
            c = fam_cnt[f]
            fh.write("\t".join([f, str(c.get("goji05", 0)), str(c.get("onekp05", 0)),
                                str(c.get("gojiPL", 0)), str(c.get("onekpPL", 0)),
                                str(plant_fam.get(f, 0)),
                                ", ".join(f"{g}({n})" for g, n in fam_gen[f].most_common(3))]) + "\n")
    print(f"候选科 {len(cand)} 个; 已登记 {len(registered)} 个; 植物库科 {len(plant_fam)} 个")
    inPL = [f for f in cand if fam_cnt[f].get("gojiPL") or fam_cnt[f].get("onekpPL")]
    print(f"其中出现在下游植物表的 {len(inPL)} 个: {inPL}")
    print(f"仅出现在 05 表的 {len(cand)-len(inPL)} 个")
    # 生成可直接粘贴的 Python 列表片段
    with open("/tmp/non_plant_family_registry.txt", "w", encoding="utf-8") as fh:
        for f in cand:
            c = fam_cnt[f]
            fh.write(f"    '{f}',  # 05:g{c.get('goji05',0)}/o{c.get('onekp05',0)}"
                     f" plant:g{c.get('gojiPL',0)}/o{c.get('onekpPL',0)}"
                     f" | {', '.join(g for g, _ in fam_gen[f].most_common(2))}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
