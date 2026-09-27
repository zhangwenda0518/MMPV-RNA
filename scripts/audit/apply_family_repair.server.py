#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""科-属校准（以属为准修正科）离线补丁。

这条补丁是 05 阶段管线内 R 层 `repair_family_from_genus` 的离线孪生实现，规则完全一致，
用途是把"闸门上线之前"已经产出的 final_integrated_classification.tsv 一次性对齐到新口径。

背景
    R 共识里的 `harmonize_family_genus` 假定"共识科"一定正确，于是用科去约束属。
    当科层被少数派工具带偏时（例：少数工具报 Mimiviridae，多数工具报 Potyviridae 且属层一致
    报 Potyvirus），这个方向性假定会把正确的属替换掉或置空 —— 这就是"科错属对"行被改坏的机制。

规则（与 virome_discovery_pipeline/virus_classifier_analysis.R 的 repair_family_from_genus 一致）
    参照 = genus_family_ref.tsv，同一条属名同时挂 NCBI 与 ICTV/VMR 两套体系的定型科。
    A `rewrite_family`：行内 Family 非空，该 Genus 被两参照都收录且给出同一个科 F2，且 F2 != 行内 Family
                        -> Family := F2（用参照里的正常大小写写法）
    B `fill_family`   ：行内 Family 为空，参照能给出定型科 -> 补上
    只有 A 要求双参照一致；单侧参照（只 NCBI 或只 VMR 收录该属）不参与改写，避免用一把尺子推翻主张。
    只改 Family，绝不改 Genus / Species / 其他列。

安全性
    * 原表零改动：只读输入表，产物写独立目录。
    * 幂等：重跑同一输入得到同一输出；对已修过的表再跑，A/B 命中数归零。
    * 自检 fail-loud：输出表里"属的双参照定型科 != 行内 Family"的行数必须为 0，否则返回非零。
    * 台账 family_repair_log.tsv 逐行记录 contig_id / 旧科 / 新科 / 属 / 动作 / 双参照证据 / 属名歧义计数。
"""
import argparse
import csv
import hashlib
import os
import sys
from collections import Counter


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm(v):
    """归一化：去空白、去引号、去尾部星号；空的/NA 归 None。"""
    if v is None:
        return None
    v = str(v).strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A", "-"):
        return None
    return v


def low(v):
    return v.lower() if v else None


def load_ref(path):
    """genus_family_ref.tsv -> {genus_lower: dict(ncbi_raw, vmr_raw, ncbi_n, vmr_n)}

    key 取属名首现（与 R 层 unique(by=g_l) 同口径）。"""
    ref = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = {c.strip().strip('"'): c for c in (rd.fieldnames or [])}
        cg = cols.get("Genus")
        cf1 = cols.get("NCBI_Family")
        cf2 = cols.get("VMR_Family")
        if not (cg and cf1 and cf2):
            raise SystemExit("ERROR: 参照表缺列，需要 Genus / NCBI_Family / VMR_Family，实到 %s" % rd.fieldnames)
        n1c, n2c = cols.get("NCBI_n"), cols.get("VMR_n")
        for r in rd:
            g = norm(r.get(cg))
            if not g:
                continue
            gl = g.lower()
            if gl in ref:
                continue
            ncbi = norm(r.get(cf1))
            vmr = norm(r.get(cf2))
            if not ncbi and not vmr:
                continue
            ref[gl] = {
                "genus_raw": g,
                "ncbi_raw": ncbi,
                "vmr_raw": vmr,
                "ncbi_n": norm(r.get(n1c)) if n1c else None,
                "vmr_n": norm(r.get(n2c)) if n2c else None,
            }
    return ref


def target_of(entry):
    """返回 (定型科原样写法, 是否双参照一致)。两参照不一致且都有值时返回 (None, False)。"""
    if entry is None:
        return None, False
    n1, n2 = entry["ncbi_raw"], entry["vmr_raw"]
    if n1 and n2:
        if n1.lower() == n2.lower():
            return n1, True
        return None, False
    if n1:
        return n1, False
    if n2:
        return n2, False
    return None, False


def main():
    ap = argparse.ArgumentParser(description="科-属校准（以属为准修正科）离线补丁")
    ap.add_argument("--input", required=True, help="输入 final_integrated_classification.tsv")
    ap.add_argument("--ref", default=os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv"),
                    help="属-科参照表 genus_family_ref.tsv")
    ap.add_argument("--outdir", required=True, help="输出目录（独立目录，原表不动）")
    args = ap.parse_args()

    if not os.path.isfile(args.input):
        print("ERROR: 找不到输入表 %s" % args.input)
        return 1
    if not os.path.isfile(args.ref):
        print("ERROR: 找不到参照表 %s" % args.ref)
        return 1

    ref = load_ref(args.ref)
    print("参照表: %s" % args.ref)
    print("  属数 %d（双参照一致 %d）" % (len(ref), sum(1 for v in ref.values() if v["ncbi_raw"] and v["vmr_raw"]
                                                        and v["ncbi_raw"].lower() == v["vmr_raw"].lower())))

    with open(args.input, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = [c.strip().strip('"') for c in (rd.fieldnames or [])]
        if "Family" not in cols or "Genus" not in cols:
            print("ERROR: 输入表缺 Family / Genus 列")
            return 1
        rows = [{k.strip().strip('"'): ("" if v is None else v) for k, v in r.items()} for r in rd]

    src_md5 = md5_of(args.input)
    print("输入表: %s" % args.input)
    print("  行数 %d / 列数 %d / md5 %s" % (len(rows), len(cols), src_md5))

    log_rows = []
    act = Counter()
    homonym = Counter()
    fam_before = Counter()
    for d in rows:
        cid = d.get("contig_id") or ""
        fam = norm(d.get("Family"))
        gen = norm(d.get("Genus"))
        if not gen:
            act["no_genus"] += 1
            continue
        entry = ref.get(gen.lower())
        tgt, dual = target_of(entry)
        if not tgt:
            act["ref_unknown_or_ambiguous"] += 1
            continue
        if fam and low(fam) != low(tgt):
            if not dual:
                act["single_side_skip"] += 1
                continue
            action = "rewrite_family"
        elif not fam:
            action = "fill_family"
        else:
            act["already_consistent"] += 1
            continue

        d["Family"] = tgt
        act[action] += 1
        fam_before[fam or "(空)"] += 1
        nn, vn = entry["ncbi_n"], entry["vmr_n"]
        amb = ""
        try:
            if nn and int(nn) > 1:
                amb = "NCBI_n=%s" % nn
        except ValueError:
            pass
        if amb:
            homonym[amb] += 1
        log_rows.append([
            cid, fam or "", tgt, gen, action,
            "NCBI=%s;VMR=%s" % (entry["ncbi_raw"] or "-", entry["vmr_raw"] or "-"),
            nn or "", vn or "", amb,
        ])

    os.makedirs(args.outdir, exist_ok=True)
    out_tsv = os.path.join(args.outdir, "final_integrated_classification.repaired.tsv")
    out_log = os.path.join(args.outdir, "family_repair_log.tsv")
    out_rep = os.path.join(args.outdir, "family_repair_report.md")

    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)
        w.writerow(cols)
        for d in rows:
            w.writerow(["" if d.get(c) is None else d.get(c) for c in cols])

    with open(out_log, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)
        w.writerow(["contig_id", "old_family", "new_family", "genus", "action",
                    "evidence", "genus_ncbi_n", "genus_vmr_n", "genus_homonym_flag"])
        for r in log_rows:
            w.writerow(r)

    # 自检：输出表里双参照定型科与行内 Family 冲突的行数必须为 0
    remain = []
    with open(out_tsv, newline="", encoding="utf-8", errors="replace") as fh:
        rd2 = csv.DictReader(fh, delimiter="\t")
        for r in rd2:
            fam = norm(r.get("Family"))
            gen = norm(r.get("Genus"))
            if not fam or not gen:
                continue
            entry = ref.get(gen.lower())
            tgt, dual = target_of(entry)
            if dual and tgt and low(fam) != low(tgt):
                remain.append((r.get("contig_id"), fam, gen, tgt))

    print("")
    print("=== 处置统计（%d 行）===" % len(rows))
    for k, v in act.most_common():
        print("  %-26s %6d" % (k, v))
    if fam_before:
        print("  被改写前的科名分布（前 15）:")
        for k, v in fam_before.most_common(15):
            print("    %-28s %6d" % (k, v))
    print("")
    print("=== 自检 ===")
    print("  输出表残留「双参照定型科 != 行内 Family」行数 = %d" % len(remain))
    for cid, fam, gen, tgt in remain[:10]:
        print("     %s  Family=%s Genus=%s 应为 %s" % (cid, fam, gen, tgt))
    ok = (len(remain) == 0)

    import datetime
    with open(out_rep, "w", newline="", encoding="utf-8") as fh:
        fh.write("# 科-属校准（以属为准修正科）报告\n\n")
        fh.write("生成时间: %s\n\n" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        fh.write("## 输入（零改动）\n\n")
        fh.write("- 表: `%s`\n- 行数: %d / 列数: %d\n- md5: `%s`\n" % (args.input, len(rows), len(cols), src_md5))
        fh.write("- 参照: `%s`（属 %d 条，双参照一致 %d 条，md5 `%s`）\n\n"
                 % (args.ref, len(ref), sum(1 for v in ref.values() if v["ncbi_raw"] and v["vmr_raw"]
                                            and v["ncbi_raw"].lower() == v["vmr_raw"].lower()), md5_of(args.ref)))
        fh.write("## 处置\n\n| action | 行数 |\n|---|---|\n")
        for k, v in act.most_common():
            fh.write("| %s | %d |\n" % (k, v))
        fh.write("\n改写涉及的原科名（前 30）：\n\n| 原科名 | 行数 |\n|---|---|\n")
        for k, v in fam_before.most_common(30):
            fh.write("| %s | %d |\n" % (k, v))
        if homonym:
            fh.write("\n属名歧义提醒（NCBI 同名多科，%d 行）：\n\n" % sum(homonym.values()))
            for k, v in homonym.most_common():
                fh.write("- %s × %d\n" % (k, v))
        fh.write("\n## 自检\n\n输出表残留双参照科-属冲突 **%d** 行（期望 0）→ %s\n\n"
                 % (len(remain), "PASS" if ok else "FAIL"))
        fh.write("## 输出\n\n- 修复后表: `%s`\n- 台账: `%s`\n- 原表未改动，回滚 = 不使用本目录产物\n" % (out_tsv, out_log))

    print("")
    print("输出: %s" % out_tsv)
    print("台账: %s" % out_log)
    print("报告: %s" % out_rep)
    print("结论: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
