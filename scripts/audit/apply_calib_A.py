#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""口径A校准补丁：按逐级相容性判据修 final_integrated_classification.tsv 的科-属嵌合行。

判据（双重参照取交集，只有两把尺子都判不相容才动手）：
  参照1 = NCBI taxdump rankedlineage.dmp（病毒部分 superkingdom 列=Realm, kingdom 列=Kingdom）
  参照2 = VMR MSL41 acvirus_db/taxa.txt（Comma 分隔，Family/Genus 定型列）
  对每一行：查 Genus 在参照里的定型科，与行内 Family 比较。
    - 两参照都已知且都不相容      -> action=blank（置空 Genus+Species，不兜底）
    - 仅一个参照已知且判不相容    -> action=review（只标记，不置空，交人工过）
    - 任一参照判相容 / 两参照都未知 -> 不动
另加不依赖参照库的票源判据（相邻等级投票工具集合不相交）作为标记列，只标记不处置。

只读原表，输出到独立目录，原表零改动。
"""
import csv
import hashlib
import os
import sys
from collections import Counter, defaultdict

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated"
SRC = os.path.join(BASE, "final_integrated_classification.tsv")
OUTDIR = os.path.join(BASE, "calibration_20260914")
OUT = os.path.join(OUTDIR, "final_integrated_classification.calibrated.tsv")
CHG = os.path.join(OUTDIR, "calibration_changes.tsv")
REP = os.path.join(OUTDIR, "calibration_report.md")

NCBI = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
VMR = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"

CALIB_COLS = ["calib_action", "calib_prev_genus", "calib_prev_species",
              "calib_evidence", "calib_flag_src", "calib_flag_single_ref"]

RANKS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PAIRS = [
    ("Realm", "Kingdom", "Realm"),
    ("Kingdom", "Phylum", "Kingdom"),
    ("Phylum", "Class", "Phylum"),
    ("Class", "Order", "Class"),
    ("Order", "Family", "Order"),
    ("Family", "Genus", "Family"),
    ("Genus", "Species", "Genus"),
]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def low(v):
    return v.lower() if v else None


def gv(d, key):
    """取归一化后的值（只用于比较，不用于输出）"""
    return norm(d.get(key))


def parse_agree(v):
    """'3/5: CAT,diamond_lca,mmseqs' -> set of tools"""
    if not v or ":" not in v:
        return None
    s = {t.strip() for t in v.split(":", 1)[1].split(",") if t.strip()}
    return s or None


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    if not os.path.isfile(SRC):
        print("ERROR: 找不到原表 %s" % SRC)
        return 1

    with open(SRC, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        orig_cols = [c.strip().strip('"') for c in rd.fieldnames]
        if "calib_action" in orig_cols:
            print("ERROR: 输入表已含 calib_action 列，说明喂进来的是校准后产物。")
            print("       本补丁只接受原始表，重跑请从原始表重新生成（保证幂等）。")
            return 1
        rows = []
        for r in rd:
            rows.append({k.strip().strip('"'): (v if v is not None else "") for k, v in r.items()})
    print("原表: %s" % SRC)
    print("  行数 %d / 列数 %d" % (len(rows), len(orig_cols)))
    src_md5 = md5_of(SRC)
    print("  md5 %s" % src_md5)

    # ---- 参照1：NCBI rankedlineage.dmp ----
    need = set()
    for d in rows:
        for L in RANKS:
            v = low(gv(d, L))
            if v:
                need.add(v)
    ref1 = {}
    with open(NCBI, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in need and nm not in ref1:
                ref1[nm] = {L: (parts[REF_IDX[L]] or None) for L in RANKS}
    print("参照1 NCBI: 需要 %d 个 name，命中 %d" % (len(need), len(ref1)))

    # ---- 参照2：VMR MSL41 taxa.txt ----
    ref2 = defaultdict(set)
    n_vmr = 0
    with open(VMR, newline="", encoding="utf-8", errors="replace") as fh:
        rd2 = csv.DictReader(fh, delimiter=",")
        tc = {c.strip().upper(): c for c in (rd2.fieldnames or [])}
        cf, cg = tc.get("FAMILY"), tc.get("GENUS")
        if not cf or not cg:
            print("ERROR: taxa.txt 表头异常 %s" % (rd2.fieldnames,))
            return 1
        for r in rd2:
            g = norm(r.get(cg))
            f = norm(r.get(cf))
            if g and f:
                ref2[g.lower()].add(f)
                n_vmr += 1
    print("参照2 VMR : %d 条映射 / %d 个属" % (n_vmr, len(ref2)))

    # ---- 逐行判定 ----
    act_counter = Counter()
    review_kind = Counter()
    boundary_counter = Counter()
    src_counter = Counter()
    singles = []
    out_rows = []
    chg_rows = []
    for d in rows:
        cid = d.get("contig_id")
        if not cid:
            continue
        fam = gv(d, "Family")
        gen = gv(d, "Genus")
        sp = gv(d, "Species")
        raw_gen = d.get("Genus") or ""
        raw_sp = d.get("Species") or ""
        action = ""
        evidence = ""
        single_ref = ""
        if fam and gen:
            gl = low(gen)
            r1 = ref1.get(gl)
            ncbi_fam = r1.get("Family") if r1 else None
            vmr_set = ref2.get(gl)
            ncbi_known = bool(ncbi_fam)
            vmr_known = bool(vmr_set)
            ncbi_bad = ncbi_known and low(ncbi_fam) != low(fam)
            vmr_bad = vmr_known and low(fam) not in {low(x) for x in vmr_set}
            ev = "Family=%s(%s)|Genus=%s(%s)|NCBI定型科=%s|VMR定型科=%s" % (
                fam, d.get("Family_agree") or "NA", gen, d.get("Genus_agree") or "NA",
                ncbi_fam or "-", ";".join(sorted(vmr_set)) if vmr_set else "-")
            if ncbi_known and vmr_known:
                if ncbi_bad and vmr_bad:
                    action = "blank"
                    evidence = ev
                elif ncbi_bad or vmr_bad:
                    action = "conflict_one_side"  # 一参照相容一参照不相容，保守不动
                    evidence = ev
            elif ncbi_bad or vmr_bad:
                action = "review"
                single_ref = "ncbi_only" if ncbi_bad else "vmr_only"
                evidence = ev
                singles.append((cid, fam, gen, ncbi_fam,
                                ";".join(sorted(vmr_set)) if vmr_set else "-", single_ref))

        # 票源判据（标记用，不处置）
        srcbad = []
        for hi, lo, _rc in PAIRS:
            sh = parse_agree(d.get(hi + "_agree"))
            sl = parse_agree(d.get(lo + "_agree"))
            if sh and sl and sh.isdisjoint(sl):
                srcbad.append("%s<->%s" % (hi, lo))
                src_counter["%s<->%s" % (hi, lo)] += 1
        flag_src = ";".join(srcbad)

        # 相邻边界矛盾（参照库全等级，诊断用）
        for hi, lo, refcol in PAIRS:
            vh, vl = gv(d, hi), gv(d, lo)
            if not vh or not vl:
                continue
            rl = ref1.get(low(vl))
            if not rl:
                continue
            exp = rl.get(refcol)
            if exp and low(exp) != low(vh):
                boundary_counter["%s<->%s" % (hi, lo)] += 1

        act_counter[action or "none"] += 1
        if action == "review":
            review_kind[single_ref] += 1
        new_d = dict(d)
        if action == "blank":
            new_d["Genus"] = ""
            new_d["Species"] = ""
        new_d["calib_action"] = action
        new_d["calib_prev_genus"] = raw_gen
        new_d["calib_prev_species"] = raw_sp
        new_d["calib_evidence"] = evidence
        new_d["calib_flag_src"] = flag_src
        new_d["calib_flag_single_ref"] = single_ref
        out_rows.append(new_d)
        if action in ("blank", "review"):
            chg_rows.append([
                cid, action, "Family<->Genus",
                raw_gen, raw_sp,
                "" if action == "blank" else raw_gen,
                "" if action == "blank" else raw_sp,
                fam or "", d.get("Order") or "", d.get("Class") or "",
                d.get("primary_tool") or "",
                d.get("Family_agree") or "NA", d.get("Genus_agree") or "NA",
                d.get("Species_agree") or "NA", evidence,
            ])

    os.makedirs(OUTDIR, exist_ok=True)
    all_cols = list(orig_cols) + CALIB_COLS
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", quotechar='"', quoting=csv.QUOTE_ALL)
        w.writerow(all_cols)
        for d in out_rows:
            w.writerow(["" if d.get(c) is None else d.get(c) for c in all_cols])
    with open(CHG, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", quotechar='"', quoting=csv.QUOTE_ALL)
        w.writerow(["contig_id", "action", "boundary", "prev_genus", "prev_species",
                    "new_genus", "new_species", "Family", "Order", "Class", "primary_tool",
                    "Family_agree", "Genus_agree", "Species_agree", "evidence"])
        for r in chg_rows:
            w.writerow(r)

    print("")
    print("=== 判定结果（%d 行）===" % len(out_rows))
    for k, v in act_counter.most_common():
        print("  %-20s %6d" % (k, v))
    print("  review 组成:")
    for k, v in review_kind.most_common():
        print("    %-14s %6d" % (k, v))
    print("")
    print("=== 参照库全等级边界矛盾（诊断，未处置）===")
    for k, v in boundary_counter.most_common():
        print("  %-22s %6d" % (k, v))
    print("")
    print("=== 票源不相交（诊断，未处置）===")
    for k, v in src_counter.most_common():
        print("  %-22s %6d" % (k, v))

    # ---- 复验：重读输出，科-属矛盾应只剩 review 行 ----
    remain = 0
    remain_ids = []
    with open(OUT, newline="", encoding="utf-8", errors="replace") as fh:
        rd3 = csv.DictReader(fh, delimiter="\t")
        for r in rd3:
            d = {k.strip().strip('"'): norm(v) for k, v in r.items()}
            fam, gen = d.get("Family"), d.get("Genus")
            if not fam or not gen:
                continue
            rl = ref1.get(low(gen))
            if not rl:
                continue
            exp = rl.get("Family")
            if exp and low(exp) != low(fam):
                remain += 1
                remain_ids.append((d.get("contig_id"), d.get("calib_action")))
    print("")
    print("=== 复验：输出表残留科-属矛盾 %d 行（应等于 review 行数）===" % remain)
    for cid, act in remain_ids[:30]:
        print("   %s  action=%s" % (cid, act))
    extra = [x for x in remain_ids if x[1] != "review"]
    if extra:
        print("   !! 异常：非 review 行残留 %d 条" % len(extra))
        for cid, act in extra[:10]:
            print("     %s action=%s" % (cid, act))
        return 2

    # ---- 报告 ----
    import datetime
    with open(REP, "w", newline="", encoding="utf-8") as fh:
        fh.write("# 口径A 校准报告\n\n")
        fh.write("生成时间: %s\n\n" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        fh.write("## 输入（零改动）\n\n")
        fh.write("- 文件: `%s`\n" % SRC)
        fh.write("- md5: `%s`\n" % src_md5)
        fh.write("- 行数: %d，列数: %d\n\n" % (len(rows), len(orig_cols)))
        fh.write("## 参照\n\n")
        fh.write("- NCBI taxdump `%s`（mtime %s, %d bytes, md5 `%s`）\n" % (
            NCBI, datetime.datetime.fromtimestamp(os.path.getmtime(NCBI)).strftime("%Y-%m-%d %H:%M"),
            os.path.getsize(NCBI), md5_of(NCBI)))
        fh.write("- VMR MSL41 `%s`（mtime %s, %d bytes, md5 `%s`）\n\n" % (
            VMR, datetime.datetime.fromtimestamp(os.path.getmtime(VMR)).strftime("%Y-%m-%d %H:%M"),
            os.path.getsize(VMR), md5_of(VMR)))
        fh.write("## 判定\n\n| action | 行数 | 处置 |\n|---|---|---|\n")
        fh.write("| blank | %d | Genus+Species 置空 |\n" % act_counter.get("blank", 0))
        fh.write("| review | %d | 只标记，人工过 |\n" % act_counter.get("review", 0))
        fh.write("| conflict_one_side | %d | 一参照相容一参照不相容，保守不动 |\n" % act_counter.get("conflict_one_side", 0))
        fh.write("| none | %d | 不动 |\n\n" % act_counter.get("none", 0))
        fh.write("### review 两个亚类（性质不同，均不置空）\n\n")
        fh.write("- `ncbi_only` %d 行：NCBI 判不相容，VMR MSL41 里无此属名。其中相当一部分的 Genus 列填的本身就不是属名"
                 "（如 `Human papillomavirus`、`uncultured partitivirus`、`Organic Lake phycodnavirus`），属于名称层级混入，需逐条看。\n"
                 % review_kind.get("ncbi_only", 0))
        fh.write("- `vmr_only` %d 行：VMR MSL41 判不相容，NCBI 无此属（新属）。这批是**分类学版本漂移**，"
                 "不是嵌合：MSL41 已把原 Mimiviridae 拆出新科 Hydriviridae，工具的库还停在老科名。"
                 "行内 Family 是过时的科名，属级主张未必错，因此不置空。\n\n" % review_kind.get("vmr_only", 0))
        fh.write("### conflict_one_side %d 行\n\n" % act_counter.get("conflict_one_side", 0))
        fh.write("两参照都有该属，但一个相容一个不相容。判据要求双参照一致才动手，故保守不动。\n\n")
        fh.write("## 输出\n\n")
        fh.write("- 校准后表: `%s`\n" % OUT)
        fh.write("- 变更清单: `%s`\n" % CHG)
        fh.write("- 原表未修改（md5 见上），回滚 = 不用这个新文件即可\n\n")
        fh.write("## review 行明细（%d 行）\n\n" % len(singles))
        fh.write("| contig_id | kind | Family | Genus | NCBI定型科 | VMR定型科 |\n|---|---|---|---|---|---|\n")
        for cid, fam, gen, n1, n2, kind in singles:
            fh.write("| %s | %s | %s | %s | %s | %s |\n" % (cid, kind, fam, gen, n1 or "-", n2 or "-"))
        fh.write("\n## 残留诊断（未处置，供后续决策）\n\n")
        fh.write("参照库全等级边界矛盾:\n\n| 边界 | 行数 |\n|---|---|\n")
        for k, v in boundary_counter.most_common():
            fh.write("| %s | %d |\n" % (k, v))
        fh.write("\n票源不相交:\n\n| 边界 | 行数 |\n|---|---|\n")
        for k, v in src_counter.most_common():
            fh.write("| %s | %d |\n" % (k, v))
    print("")
    print("输出: %s" % OUT)
    print("清单: %s" % CHG)
    print("报告: %s" % REP)
    return 0


if __name__ == "__main__":
    sys.exit(main())
