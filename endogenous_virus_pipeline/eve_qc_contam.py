#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eve_qc_contam.py — 测序污染筛查（对一次已跑完的筛查结果事后追加）
==================================================================
为什么需要它：候选里混进 phiX174（Illumina 建库的 spike-in）或其它大肠杆菌噬菌体序列时，
它们会打中病毒参考库里的噬菌体条目，然后被当成"病毒信号"一路走到汇总表。上游同类研究
实测 Phage 层占元件总数 **28%**（11,855/41,887），其质控记录直接写明"Phage 类含 phiX
污染"——也就是说这一层不做筛查，王国级的组成结论会被污染抬高。我们的管线此前没有这一步。

为什么做成**独立脚本**而不是又一个 stage：判定只需要已有的命中表，不需要重跑比对，所以
任何一次历史结果都能事后补筛；也避免改动主编排的阶段注册与断点语义。上游研究用的也是
这种"事后 QC pass"的形式（`eve_qc_clean.py`），这里保持了同样的形状。

判定口径（与上游一致，按参考 accession 前缀 + 标题关键词）：
  命中 phiX174 / 大肠杆菌噬菌体的参考序列 = 几乎必然是 spike-in 或细菌污染。
  **其它噬菌体不在本脚本范围内** —— 它们更可能是组装里带的细菌序列，属于 ICTV
  genome-type 分类那一层的职责（"Phage (bacterial)" 类），不在这里混为一谈。

用法
----
  python3 eve_qc_contam.py --outdir <筛查输出根目录> [--frac 0.5] [--json F] [--md F]

产物（默认写在 --outdir 下）
  qc_contamination.tsv   逐 contig 的污染命中占比与最强污染参考
  QC_CONTAMINATION.md    可直接引用的质控记录（含建议剔除清单）
"""
import argparse
import collections
import json
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from eve_scan_core import (CONTAM_HDR, STAGE_DIRS, contamination_scan,
                          is_contaminant_hit, logmsg)   # noqa: E402

# 两种输入布局 —— 移植时踩到的差异: 上游那套 (eve_kingdom) 的逐基因组目录里**没有**
# Stage-1 原始命中, 只有 `<NAME>.s1_best.tsv`(每位点的最佳命中)。于是"污染占比"的分母
# 从"该 contig 的全部命中"变成"该基因组的全部位点", 问题也相应变成更能直接回答的那个:
# **报出来的位点里有多少是 phiX**。代价是拿不到标题通道(只有 accessions), 所以那一档
# 只能靠 accession 前缀判 —— 这一点必须写在产物里, 不能悄悄降级。
LAYOUT_OURS = "ours"
LAYOUT_KINGDOM = "kingdom"


def detect_layout(outdir):
    root = Path(outdir)
    if (root / STAGE_DIRS["loci"]).is_dir():
        return LAYOUT_OURS
    if any(root.glob("*/*.s1_best.tsv")):
        return LAYOUT_KINGDOM
    raise SystemExit("既没有 %s/ 也没有 */<NAME>.s1_best.tsv, 认不出这是哪套输出: %s"
                     % (STAGE_DIRS["loci"], outdir))


def _viral_supported(name_dir, name):
    """该基因组被判 viral_supported 的位点集合 (从 <NAME>.s2_verdict.tsv 读).

    为什么要做这一步: `s1_best.tsv` 里是**全部**位点, 而真正进结论的只有
    viral_supported 那一批。拿"全部位点"当分母会把污染占比算小 (大量未通过的位点被
    算进分母), 拿"污染位点数"直接对全库报又会算大 (那些位点可能根本没进结论)。
    实测 Abies_alba 有 540,547 个位点、其中 149,372 个最佳命中是 phiX, 但它的
    viral_supported 只有一小部分 —— 不区分就会得出"全库七成是 phiX"这种错结论。
    """
    vp = Path(name_dir) / ("%s.s2_verdict.tsv" % name)
    if not vp.is_file():
        return None
    out = set()
    with open(vp, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        if "locus" not in hdr or "verdict" not in hdr:
            return None
        li, vi = hdr.index("locus"), hdr.index("verdict")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(li, vi) and p[vi] == "viral_supported":
                out.add(p[li])
    return out


def scan_kingdom(outdir, frac_thr=0.5, quiet=False):
    """按 `<NAME>/<NAME>.s1_best.tsv` 扫: 每位点看最佳命中的参考是不是污染源。

    返回 (逐位点行, 汇总, 超阈值行). 有 `s2_verdict.tsv` 时**同时**给两个口径:
      `n_hits` / `n_contam_hits`  全部位点口径
      `n_vs` / `n_contam_vs`      **viral_supported 口径**（决定性的那个）
    `contam_frac` 用 viral_supported 口径算 (没有该表时退回全位点口径)。
    """
    root = Path(outdir)
    rows = []
    n_genomes = 0
    per_genome = collections.Counter()
    top_refs = collections.Counter()
    top_refs_vs = collections.Counter()
    for best in sorted(root.glob("*/*.s1_best.tsv")):
        name = best.parent.name
        vs = _viral_supported(best.parent, name)
        n_loc = n_vs = 0
        hit = []
        hit_vs = []
        with open(best, encoding="utf-8") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                p += [""] * (len(hdr) - len(p))
                r = dict(zip(hdr, p))
                n_loc += 1
                locus = r.get("locus", "")
                in_vs = (vs is None) or (locus in vs)
                if vs is not None and in_vs:
                    n_vs += 1
                if is_contaminant_hit(r.get("ref_sseqid", ""), ""):
                    hit.append(r)
                    if in_vs:
                        hit_vs.append(r)
        n_genomes += 1
        if not hit:
            continue
        for r in hit:
            top_refs[r.get("ref_sseqid", "")] += 1
        for r in hit_vs:
            top_refs_vs[r.get("ref_sseqid", "")] += 1
        # 决定性口径: 污染位点里有多少真的进了结论
        denom = n_vs if vs is not None else n_loc
        frac = len(hit_vs) / denom if denom else 0.0
        if frac >= frac_thr:
            per_genome[name] = len(hit_vs)
        rows.append({"genome": name, "query_id": "(locus-level)",
                     "n_hits": str(n_loc), "n_contam_hits": str(len(hit)),
                     "contam_frac": "%.3f" % frac,
                     "n_vs": str(denom), "n_contam_vs": str(len(hit_vs)),
                     "top_contam_ref": ";".join(sorted(
                         {h.get("ref_sseqid", "") for h in hit})[:5]),
                     "top_contam_title": "(kingdom 布局无标题通道)"})
    rows.sort(key=lambda r: (r["genome"], r["query_id"]))
    flag = [r for r in rows if float(r["contam_frac"]) >= frac_thr]
    summary = {
        "layout": LAYOUT_KINGDOM,
        "n_genomes_scanned": n_genomes,
        "n_contigs_with_contam_hits": len(rows),
        "n_contigs_flagged": len(flag),
        "frac_threshold": frac_thr,
        "n_genomes_with_flagged": len(per_genome),
        "total_loci": sum(int(r["n_hits"]) for r in rows),
        "total_contam_loci": sum(int(r["n_contam_hits"]) for r in rows),
        "total_viral_supported": sum(int(r["n_vs"]) for r in rows),
        "total_contam_viral_supported": sum(int(r["n_contam_vs"]) for r in rows),
        "top_contam_refs": dict(top_refs.most_common(10)),
        "top_contam_refs_viral_supported": dict(top_refs_vs.most_common(10)),
    }
    if not quiet:
        logmsg("污染筛查(kingdom 布局): 扫了 %d 个基因组; 有污染命中的 %d 个。"
               "决定性口径 —— viral_supported 位点里污染占 %d / %d = %.2f%%"
               % (n_genomes, len(rows), summary["total_contam_viral_supported"],
                  summary["total_viral_supported"],
                  100.0 * summary["total_contam_viral_supported"]
                  / max(summary["total_viral_supported"], 1)))
    return rows, summary, flag


def scan_outdir(outdir, frac_thr=0.5, quiet=False):
    """扫 <outdir>/01_Loci/*/ 下的 Stage-1 原始命中 -> (逐 contig 行, 汇总)."""
    root = Path(outdir)
    loci_root = root / STAGE_DIRS["loci"]
    if not loci_root.is_dir():
        raise SystemExit("无 %s/ 目录, 先跑筛查阶段: %s"
                         % (STAGE_DIRS["loci"], outdir))
    rows = []
    per_genome = collections.Counter()
    top_refs = collections.Counter()
    n_genomes = 0
    for d in sorted(loci_root.iterdir()):
        if not d.is_dir():
            continue
        raw = d / ("%s.s1_raw.tsv" % d.name)
        if not raw.is_file() or raw.stat().st_size == 0:
            continue
        n_genomes += 1
        tmp = d / "_qc_contam.tsv"
        n_q, n_hit = contamination_scan(raw, tmp)
        if not n_hit:
            tmp.unlink(missing_ok=True)
            continue
        with open(tmp, encoding="utf-8") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                p = line.rstrip("\n").split("\t")
                p += [""] * (len(hdr) - len(p))
                r = dict(zip(hdr, p))
                r["genome"] = d.name
                rows.append(r)
                if float(r["contam_frac"]) >= frac_thr:
                    per_genome[d.name] += 1
                top_refs[r["top_contam_ref"]] += 1
        tmp.unlink(missing_ok=True)
    rows.sort(key=lambda r: (r["genome"], r["query_id"]))
    flag = [r for r in rows if float(r["contam_frac"]) >= frac_thr]
    summary = {
        "n_genomes_scanned": n_genomes,
        "n_contigs_with_contam_hits": len(rows),
        "n_contigs_flagged": len(flag),
        "frac_threshold": frac_thr,
        "n_genomes_with_flagged": len(per_genome),
        "top_contam_refs": dict(top_refs.most_common(10)),
    }
    if not quiet:
        logmsg("污染筛查: 扫了 %d 个基因组; 有污染命中的 contig %d 条, 其中 "
               "占比 >=%.2f 的 %d 条 (分布在 %d 个基因组)"
               % (n_genomes, len(rows), frac_thr, len(flag), len(per_genome)))
    return rows, summary, flag


def write_outputs(outdir, rows, summary, flag, json_path="", md_path=""):
    root = Path(outdir)
    tsv = root / "qc_contamination.tsv"
    with open(tsv, "w", encoding="utf-8", newline="\n") as fo:
        fo.write("genome\t" + CONTAM_HDR.rstrip("\n")
                 + "\tn_viral_supported\tn_contam_viral_supported\n")
        for r in rows:
            fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                r["genome"], r["query_id"], r["n_hits"], r["n_contam_hits"],
                r["contam_frac"], r["top_contam_ref"], r["top_contam_title"],
                r.get("n_vs", ""), r.get("n_contam_vs", "")))
    md = root / "QC_CONTAMINATION.md"
    with open(md, "w", encoding="utf-8", newline="\n") as fo:
        fo.write("# 测序污染筛查 (phiX174 / 大肠杆菌噬菌体)\n\n")
        lay = summary.get("layout", LAYOUT_OURS)
        fo.write("输入布局: **%s**%s\n\n" % (lay, (
            "" if lay == LAYOUT_OURS else
            "（该布局只有每位点的最佳命中，故占比的分母是**该基因组位点总数**，"
            "且**没有标题通道** —— 只用 accession 前缀判，检出率会低于 ours 布局）")))
        fo.write("判定: 命中参考的 accession 以 %s 开头, 或标题含 %s。\n\n"
                 % ("/".join(("NP_0406", "NP_0407", "YP_51237")),
                    ", ".join(("phix", "phi x", "escherichia phage"))))
        fo.write("| 项目 | 数值 |\n|---|---|\n")
        fo.write("| 扫描的基因组 | %d |\n" % summary["n_genomes_scanned"])
        fo.write("| 有污染命中的 contig | %d |\n" % summary["n_contigs_with_contam_hits"])
        fo.write("| 其中污染占比 >= %.2f 的 contig | %d |\n"
                 % (summary["frac_threshold"], summary["n_contigs_flagged"]))
        fo.write("| 涉及基因组 | %d |\n" % summary["n_genomes_with_flagged"])
        if "total_viral_supported" in summary:
            vs = summary["total_viral_supported"]
            cv = summary["total_contam_viral_supported"]
            fo.write("| **viral_supported 位点合计（仅这些受污染基因组）** | **%d** |\n" % vs)
            fo.write("| **其中最佳命中是污染源的** | **%d（%.2f%%）** |\n"
                     % (cv, 100.0 * cv / max(vs, 1)))
            fo.write("\n> 注意分母：上表只统计**出现过污染命中的那些基因组**的 viral_supported "
                     "位点。要看它对全库的影响，需再除以全库 viral_supported 总数"
                     "（本次为 207,415，见 `EVE_DB/DB_SUMMARY.md`）。\n")
        fo.write("\n## 出现最多的污染参考 (Top 10)\n\n| 参考 | contig 数 |\n|---|---|\n")
        for ref, n in sorted(summary["top_contam_refs"].items()):
            fo.write("| %s | %d |\n" % (ref, n))
        fo.write("\n## 建议剔除清单 (contam_frac >= %.2f)\n\n"
                 % summary["frac_threshold"])
        fo.write("逐 contig 明细见 `qc_contamination.tsv`。下游按 `genome` + `query_id` "
                 "过滤即可；阈值可按需调整（本表未替你下结论）。\n\n")
        if flag:
            fo.write("| genome | contig | contam_frac | top_ref |\n|---|---|---|---|\n")
            for r in flag:
                fo.write("| %s | %s | %s | %s |\n" % (
                    r["genome"], r["query_id"], r["contam_frac"], r["top_contam_ref"]))
        else:
            fo.write("（无）\n")
    if json_path:
        with open(json_path, "w", encoding="utf-8", newline="\n") as fo:
            json.dump(summary, fo, ensure_ascii=False, indent=1)
    return tsv, md


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="对一次已跑完的 EVE 筛查结果事后做测序污染筛查")
    ap.add_argument("--outdir", required=True, help="筛查输出根目录（含 01_Loci/）")
    ap.add_argument("--layout", default="auto",
                    choices=["auto", LAYOUT_OURS, LAYOUT_KINGDOM],
                    help="ours=本管线布局(01_Loci/<NAME>/<NAME>.s1_raw.tsv, 占比按"
                         "该 contig 全部命中算); kingdom=上游布局(<NAME>/<NAME>.s1_best.tsv, "
                         "占比按该基因组位点数算); auto=自动识别")
    ap.add_argument("--frac", type=float, default=0.5,
                    help="污染命中占比达到多少算'建议剔除'（默认 0.5；只影响标记，"
                         "明细表恒含全部有污染命中的条目）")
    ap.add_argument("--json", default="", help="汇总 JSON 输出路径")
    ap.add_argument("--md", default="", help="质控记录输出路径（默认 <outdir>/QC_CONTAMINATION.md）")
    a = ap.parse_args(argv)
    layout = a.layout if a.layout != "auto" else detect_layout(a.outdir)
    logmsg("输入布局: %s" % layout)
    if layout == LAYOUT_OURS:
        rows, summary, flag = scan_outdir(a.outdir, a.frac)
    else:
        rows, summary, flag = scan_kingdom(a.outdir, a.frac)
    tsv, md = write_outputs(a.outdir, rows, summary, flag, a.json, a.md)
    logmsg("产物: %s / %s" % (tsv, md))
    return 0


if __name__ == "__main__":
    sys.exit(main())
