"""05 分类选择层: 科属不相容现状实测 (亲代 vs 校准后).

判据同 virus_classifier_analysis.R::enforce_rank_containment:
  属 -> 定型科 (NCBI rankedlineage.dmp / VMR MSL41 taxa.txt)
  两参照都收录该属, 且两参照给出的科都与行内 Family 不同 -> both_bad (会被置空)
  单侧不同 -> one_side (保守不动)
"""
import csv
import os
from collections import Counter

REF = "/home/zhangwenda/database/taxonomy/genus_family_ref.tsv"
GOJI = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated"


def norm(x):
    x = (x or "").replace('"', "").replace("*", "").strip()
    if x == "" or x.upper() in ("NA", "N/A", "-"):
        return None
    return x


def famset(x):
    x = norm(x)
    if x is None:
        return set()
    out = set()
    for part in x.replace(",", ";").split(";"):
        part = part.replace(" ", "").lower()
        if part:
            out.add(part)
    return out


ref = {}
with open(REF, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        g = norm(r.get("Genus"))
        if not g:
            continue
        gl = g.lower()
        if gl not in ref:
            ref[gl] = (famset(r.get("NCBI_Family")), famset(r.get("VMR_Family")))
print(f"参照表属数: {len(ref)}")


def scan(path, label):
    print("\n" + "=" * 88)
    print(f"[{label}] {os.path.basename(path)}")
    print("=" * 88)
    rows = list(csv.DictReader(open(path, errors="replace"), delimiter="\t"))
    print("  行数:", len(rows))
    st = Counter()
    acts = Counter()
    for r in rows:
        fam = norm(r.get("Family"))
        gen = norm(r.get("Genus"))
        if "calib_action" in r:
            acts[(r.get("calib_action") or "").strip()] += 1
        if gen is None:
            st["Genus 为空"] += 1
            continue
        if fam is None:
            st["Family 为空, Genus 有值"] += 1
            continue
        ncbi_f, vmr_f = ref.get(gen.lower(), (set(), set()))
        famc = fam.replace(" ", "").lower()
        ncbi_hit = bool(ncbi_f)
        vmr_hit = bool(vmr_f)
        if not ncbi_hit and not vmr_hit:
            st["属无跨参照记录(无法判)"] += 1
            continue
        ncbi_bad = ncbi_hit and famc not in ncbi_f
        vmr_bad = vmr_hit and famc not in vmr_f
        if ncbi_bad and vmr_bad:
            st["both_bad 双侧不相容"] += 1
        elif ncbi_bad or vmr_bad:
            st["one_side 单侧不相容"] += 1
        else:
            st["相容"] += 1
    for k, v in st.most_common():
        print(f"    {k:<26} {v}")
    if acts:
        print("  calib_action 分布:", dict(acts))
    return st


scan(f"{GOJI}/final_integrated_classification.tsv", "亲代表 (08-12, 原表未改)")
scan(f"{GOJI}/calibration_20260914/final_integrated_classification.calibrated.tsv", "校准表 (09-14)")
