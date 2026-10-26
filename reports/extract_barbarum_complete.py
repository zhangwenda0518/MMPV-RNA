#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_barbarum_complete.py

把 `03_known_virus/known_virus_all_v2` 与 `known_virus_all_v3` 中
**宁夏枸杞 (Lycium barbarum)** 的样本数据，完整镜像到
`03_known_virus/known_virus_barbarum/`（纯新增，不改任何上游文件）。

口径
----
物种标签与 `00_scripts/split_known_virus_by_species.py` 完全一致：
`01_data/{gsa.tab, ncbi.DNA.tab, ncbi.RNA.tab}` 的**列头**即分组名。
barbarum 分组 = {Lycium barbarum-RNA, RNA-Lycium barbarum,
                 TRANSCRIPTOMIC SINGLE CELL-L. b}
`ncbi.tab` 是 DNA+RNA 的并集，读入但不额外贡献分组（union 去重，避免重复计数）。

切分规则
--------
1. per-sample 文件/目录：名字以 `<Run>` 开头（Run 后紧跟 `.` 或 `_`），
   或 `sample_<Run>...` 形式 -> 命中 barbarum 才复制。
2. 表 (.tsv/.csv)：按 `Sample`/`sample` 列过滤后写出；
   若表头没有该列则整体复制。
3. parquet：按 `Sample` 列过滤重写（需要 pyarrow，缺失时整体复制）。
4. 病毒级 / 全局产物（病毒注释、比对、索引、图、报告…）：整体复制，
   无法按样本切分（README 会标注）。

合并策略
--------
**v2 优先**：先复制 v2 的 barbarum 子集，再用 v3 补 v2 缺失
（目标路径已存在即跳过）。两版的 stage 目录名一致
（05_assembly / 06_extraction 用下划线，01/03 用点号）。

用法
----
    python3 extract_barbarum_complete.py --dry-run
    python3 extract_barbarum_complete.py
    python3 extract_barbarum_complete.py --link                 # 硬链接（省空间）
    python3 extract_barbarum_complete.py --stages 05 06         # 只做指定 stage
"""
import argparse
import csv
import json
import os
import re
import shutil
import sys
from collections import OrderedDict

HOME = os.path.expanduser("~")
GOJI = os.path.join(HOME, "MMPV-paper", "goji-virome")

SPECIES_MAP = {
    "Lycium barbarum-RNA":             "Lycium_barbarum",
    "RNA-Lycium barbarum":             "Lycium_barbarum",
    "TRANSCRIPTOMIC SINGLE CELL-L. b": "Lycium_barbarum",
}
TARGET = "Lycium_barbarum"
TAB_FILES = ["gsa.tab", "ncbi.DNA.tab", "ncbi.RNA.tab", "ncbi.tab"]

RUN_RE = re.compile(r"^([A-Za-z]{3}\d+)(?=[._]|$)")
SAMPLE_RE = re.compile(r"^sample_([A-Za-z]{3}\d+)")
STAGE_RE = re.compile(r"^\d{2}_[a-z_]+$")

RUNS = set()          # barbarum run 集合（main 中填充）
PLANNED = set()       # 已计划/已写入的目标相对路径（用于 v2 优先）


# ------------------------------------------------------------------ tab 读取
def read_tab(path):
    """读转置表：第 1 行是分组名，每列是该分组的 run 列表。"""
    groups = OrderedDict()
    with open(path, encoding="utf-8-sig", newline="") as f:
        rows = [r.rstrip("\n").rstrip("\r").split("\t") for r in f]
    if not rows:
        return groups
    for ci, g in enumerate(rows[0]):
        g = g.strip()
        if not g:
            continue
        groups[g] = [r[ci].strip() for r in rows[1:]
                     if ci < len(r) and r[ci].strip()]
    return groups


def load_runs(data_dir):
    runs, per_file = set(), OrderedDict()
    for tf in TAB_FILES:
        p = os.path.join(data_dir, tf)
        if not os.path.isfile(p):
            per_file[tf] = {"present": False}
            continue
        groups = read_tab(p)
        hit_groups = [g for g in groups if SPECIES_MAP.get(g) == TARGET]
        got = []
        for g in hit_groups:
            got.extend(groups[g])
        per_file[tf] = {"present": True, "groups": hit_groups,
                        "n": len(set(got))}
        runs.update(got)
    return runs, per_file


def find_run(name):
    m = RUN_RE.match(name)
    if m and m.group(1) in RUNS:
        return m.group(1)
    m = SAMPLE_RE.match(name)
    if m and m.group(1) in RUNS:
        return m.group(1)
    return None


# ------------------------------------------------------------ 各 stage 收集器
def _file(src, rel):
    return dict(src=src, rel=rel, kind="file")


def _tree(src, rel):
    return dict(src=src, rel=rel, kind="tree")


def _table(src, rel, col="Sample"):
    return dict(src=src, rel=rel, kind="table_tsv", col=col)


def c_01(src, runs):
    """01_detection"""
    E = []
    for sub in ("bam", "stat"):                       # 全是 per-sample
        d = os.path.join(src, sub)
        if os.path.isdir(d):
            for fn in sorted(os.listdir(d)):
                fp = os.path.join(d, fn)
                if os.path.isfile(fp) and find_run(fn):
                    E.append(_file(fp, "%s/%s" % (sub, fn)))
    d = os.path.join(src, "depth_plots")              # 目录级 per-sample
    if os.path.isdir(d):
        for dn in sorted(os.listdir(d)):
            dp = os.path.join(d, dn)
            if os.path.isdir(dp) and find_run(dn):
                E.append(_tree(dp, "depth_plots/%s" % dn))
    d = os.path.join(src, "logs")                     # sample_<Run>.log + 全局
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            fp = os.path.join(d, fn)
            if not os.path.isfile(fp):
                continue
            if find_run(fn) or fn.startswith(("index_build", "pipeline_",
                                              "run")):
                E.append(_file(fp, "logs/%s" % fn))
    d = os.path.join(src, "summary")                  # 表按 Sample 过滤
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            fp = os.path.join(d, fn)
            if os.path.isdir(fp):
                E.append(_tree(fp, "summary/%s" % fn))
            elif fn.endswith(".tsv"):
                E.append(_table(fp, "summary/%s" % fn, "Sample"))
            else:
                E.append(_file(fp, "summary/%s" % fn))
    d = os.path.join(src, "batches")                  # parquet 过滤
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            fp = os.path.join(d, fn)
            if fn.endswith(".parquet"):
                E.append(dict(src=fp, rel="batches/%s" % fn,
                              kind="parquet", col="Sample"))
            else:
                E.append(_file(fp, "batches/%s" % fn))
    for sub in ("index", "plots", "sample_distribution",
                "coabundance", "pseudo_quant"):       # 整体
        d = os.path.join(src, sub)
        if os.path.isdir(d):
            E.append(_tree(d, sub))
    for fn in sorted(os.listdir(src)):                # 顶层 run*.log / run.ok
        fp = os.path.join(src, fn)
        if os.path.isfile(fp) and (fn.startswith("run") or fn.endswith(".log")):
            E.append(_file(fp, fn))
    return E


def c_02(src, runs):
    """02_filtering：表按 Sample 过滤，其余整体"""
    E = []
    for x in sorted(os.listdir(src)):
        p = os.path.join(src, x)
        if os.path.isdir(p):
            E.append(_tree(p, x))
        elif x.endswith(".tsv"):
            E.append(_table(p, x, "Sample"))
        else:
            E.append(_file(p, x))
    return E


def c_03(src, runs):
    """03_variants"""
    E = []
    for sub in ("virus-variants", "virus-consensus",
                "virus-SnpEff", "virus-SNPGenie"):
        d = os.path.join(src, sub)
        if not os.path.isdir(d):
            continue
        for v in sorted(os.listdir(d)):
            vd = os.path.join(d, v)
            if not os.path.isdir(vd):
                continue
            for s in sorted(os.listdir(vd)):
                sd = os.path.join(vd, s)
                if not find_run(s):
                    continue
                rel = "%s/%s/%s" % (sub, v, s)
                E.append(_tree(sd, rel) if os.path.isdir(sd)
                         else _file(sd, rel))
    for sub in ("virus-bam", "virus_reads"):          # 文件级 per-sample
        d = os.path.join(src, sub)
        if os.path.isdir(d):
            for fn in sorted(os.listdir(d)):
                if find_run(fn):
                    E.append(_file(os.path.join(d, fn), "%s/%s" % (sub, fn)))
    d = os.path.join(src, "logs")                     # per-sample + 全局
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            fp = os.path.join(d, fn)
            if os.path.isfile(fp) and (find_run(fn)
                                       or fn.startswith(("run", "bowtie2", "index"))):
                E.append(_file(fp, "logs/%s" % fn))
    for sub in ("virus-annotations", "virus-fasta", "virus-index"):  # 病毒级
        d = os.path.join(src, sub)
        if os.path.isdir(d):
            E.append(_tree(d, sub))
    d = os.path.join(src, "summary")
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            fp = os.path.join(d, fn)
            if os.path.isdir(fp):
                E.append(_tree(fp, "summary/%s" % fn))
            elif fn.endswith(".tsv"):
                E.append(_table(fp, "summary/%s" % fn, "Sample"))
            else:
                E.append(_file(fp, "summary/%s" % fn))
    for fn in sorted(os.listdir(src)):                # 顶层 run*.log
        fp = os.path.join(src, fn)
        if os.path.isfile(fp) and fn.startswith("run"):
            E.append(_file(fp, fn))
    return E


def c_05(src, runs):
    """05_assembly：<virus>/<Run>_<acc>/ 整树 + <Run>_<acc>_assembly.log"""
    E = []
    for v in sorted(os.listdir(src)):
        vd = os.path.join(src, v)
        if not os.path.isdir(vd):
            E.append(_file(vd, v))
            continue
        for s in sorted(os.listdir(vd)):
            sd = os.path.join(vd, s)
            if not find_run(s):
                continue
            rel = "%s/%s" % (v, s)
            E.append(_tree(sd, rel) if os.path.isdir(sd) else _file(sd, rel))
    return E


def c_06(src, runs):
    """06_extraction：<virus>/<Run>.<acc>.full.fasta + 表过滤"""
    E = []
    for x in sorted(os.listdir(src)):
        p = os.path.join(src, x)
        if os.path.isdir(p):
            for fn in sorted(os.listdir(p)):
                if find_run(fn):
                    E.append(_file(os.path.join(p, fn), "%s/%s" % (x, fn)))
        elif x.endswith(".tsv"):
            E.append(_table(p, x, "sample"))
        else:
            E.append(_file(p, x))
    return E


def c_default(src, runs):
    """07/08/09/10：病毒级/全局，整体复制"""
    E = []
    for x in sorted(os.listdir(src)):
        p = os.path.join(src, x)
        E.append(_tree(p, x) if os.path.isdir(p) else _file(p, x))
    return E


COLLECTORS = {"01": c_01, "02": c_02, "03": c_03, "05": c_05, "06": c_06}


def collect(stage, src, runs):
    fn = COLLECTORS.get(stage[:2], c_default)
    return fn(src, runs)


# ---------------------------------------------------------------- 写盘 / 估算
def cell_matches(v):
    """表单元格是否属于 barbarum：可能是纯 Run，也可能是 <Run>_<acc> / <Run>.<acc>"""
    if v in RUNS:
        return True
    head = v.split("_")[0].split(".")[0]
    return head in RUNS


def filter_table(src, dst, col, est):
    try:
        with open(src, encoding="utf-8-sig", newline="") as f:
            rows = list(csv.reader(f, delimiter="\t"))
    except Exception:
        if est:
            return 0, 0
        shutil.copy2(src, dst)
        return 0, 0
    if not rows:
        return 0, 0
    hdr = rows[0]
    ci = hdr.index(col) if col in hdr else None
    if ci is None:
        if not est:
            shutil.copy2(src, dst)
        return len(rows) - 1, 0
    kept = [r for r in rows[1:] if ci < len(r) and cell_matches(r[ci])]
    if not est:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with open(dst, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(hdr)
            w.writerows(kept)
    return len(rows) - 1, len(kept)


def filter_parquet(src, dst, col):
    try:
        import pyarrow.parquet as pq
        import pyarrow.compute as pc
    except Exception:
        return None
    try:
        t = pq.read_table(src)
        if col not in t.column_names:
            return None
        mask = pc.is_in(t[col], value_set=__import__("pyarrow").array(sorted(RUNS)))
        out = t.filter(mask)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        pq.write_table(out, dst)
        return t.num_rows, out.num_rows
    except Exception:
        return None


def copy_tree(src, dst, mode):
    for root, _dirs, files in os.walk(src):
        rel = os.path.relpath(root, src)
        d = dst if rel == "." else os.path.join(dst, rel)
        os.makedirs(d, exist_ok=True)
        for f in files:
            s = os.path.join(root, f)
            t = os.path.join(d, f)
            if os.path.lexists(t):
                continue
            link_file(s, t, mode)


def link_file(s, t, mode):
    """按 mode 落地单个文件：copy / symlink（绝对路径）/ hardlink（失败回落 copy）"""
    if mode == "symlink":
        os.symlink(os.path.abspath(s), t)
    elif mode == "hardlink":
        try:
            os.link(s, t)
        except OSError:
            shutil.copy2(s, t)
    else:
        shutil.copy2(s, t)


def tree_stat(path):
    total = 0
    n = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
                n += 1
            except OSError:
                pass
    return total, n


def handle(ent, out_root, opts, stat, src_name, stage):
    dst = os.path.join(out_root, ent["rel"])
    if os.path.lexists(dst) or ent["rel"] in PLANNED:
        stat["skipped_exists"] += 1
        return
    kind = ent["kind"]
    src = ent["src"]
    size = 0
    nfile = 0
    if kind == "file":
        size = os.path.getsize(src)
        nfile = 1
        if not opts["dry"]:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            link_file(src, dst, opts["mode"])
    elif kind == "tree":
        if opts["dry"]:
            size, nfile = tree_stat(src)
        else:
            copy_tree(src, dst, opts["mode"])
            size, nfile = tree_stat(dst)
    elif kind == "table_tsv":
        size = os.path.getsize(src)
        tot, kept = filter_table(src, dst, ent["col"], opts["dry"])
        nfile = 1
        stat["rows_in"] += tot
        stat["rows_kept"] += kept
    elif kind == "parquet":
        size = os.path.getsize(src)
        nfile = 1
        if not opts["dry"]:
            r = filter_parquet(src, dst, ent["col"])
            if r is None:
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy2(src, dst)
        else:
            pass
    PLANNED.add(ent["rel"])
    key = "%s/%s" % (src_name, stage)
    st = stat["by_stage"].setdefault(key, {"entries": 0, "files": 0, "bytes": 0})
    st["entries"] += 1
    st["files"] += nfile
    st["bytes"] += size
    stat["entries"] += 1
    stat["files"] += nfile
    stat["bytes"] += size


# ---------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=GOJI)
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--sources", nargs="*",
                    default=["known_virus_all_v2", "known_virus_all_v3"])
    ap.add_argument("--stages", nargs="*", default=None,
                    help="只处理这些 stage 前缀，如: --stages 05 06")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--mode", choices=["copy", "symlink", "hardlink"],
                    default="copy",
                    help="copy=真实复制（默认）；symlink=符号链接；hardlink=硬链接")
    a = ap.parse_args()

    global RUNS
    root = os.path.abspath(a.root)
    data_dir = a.data_dir or os.path.join(root, "01_data")
    out = a.out or os.path.join(root, "03_known_virus", "known_virus_barbarum")
    opts = {"dry": a.dry_run, "mode": a.mode}

    RUNS, per_file = load_runs(data_dir)
    print("[info] barbarum runs = %d" % len(RUNS))
    for tf, info in per_file.items():
        print("       %-16s %s" % (tf, info))

    stat = {"mode": "dry-run" if a.dry_run else a.mode,
            "out": out, "runs_total": len(RUNS), "tab_files": per_file,
            "sources": a.sources, "entries": 0, "files": 0, "bytes": 0,
            "skipped_exists": 0, "rows_in": 0, "rows_kept": 0,
            "by_stage": OrderedDict()}

    if not a.dry_run:
        os.makedirs(out, exist_ok=True)

    for s in a.sources:
        sroot = os.path.join(root, "03_known_virus", s)
        if not os.path.isdir(sroot):
            print("[warn] 源不存在: %s" % sroot)
            continue
        stages = [d for d in sorted(os.listdir(sroot)) if STAGE_RE.match(d)]
        if a.stages:
            stages = [d for d in stages
                      if any(d.startswith(x) for x in a.stages)]
        print("\n[source] %s -> %d 个 stage: %s"
              % (s, len(stages), ", ".join(stages)))
        for st in stages:
            src = os.path.join(sroot, st)
            entries = collect(st, src, RUNS)
            for e in entries:                       # 补上 stage 前缀，还原 01..10 层级
                e["rel"] = "%s/%s" % (st, e["rel"])
            before = stat["entries"]
            for e in entries:
                handle(e, out, opts, stat, s, st)
            print("  %-22s 条目 %5d" % (st, stat["entries"] - before))
        inp = os.path.join(sroot, "input")
        if os.path.isdir(inp):
            entries = [_file(os.path.join(inp, fn), "input/%s" % fn)
                       for fn in sorted(os.listdir(inp))
                       if os.path.isfile(os.path.join(inp, fn))
                       and find_run(fn)]
            before = stat["entries"]
            for e in entries:
                handle(e, out, opts, stat, s, "input")
            print("  %-22s 条目 %5d" % ("input", stat["entries"] - before))

    mb = stat["bytes"] / 1024.0 / 1024.0
    print("\n[done] 条目 %d / 文件 %d / 合计 %.1f MB (%.2f GB)"
          % (stat["entries"], stat["files"], mb, mb / 1024.0))
    print("       跳过已存在 %d / 表行 %d -> %d"
          % (stat["skipped_exists"], stat["rows_in"], stat["rows_kept"]))

    if not a.dry_run:
        with open(os.path.join(out, "_extract_report.json"), "w",
                  encoding="utf-8") as f:
            json.dump(stat, f, ensure_ascii=False, indent=2)
        print("[write] %s" % os.path.join(out, "_extract_report.json"))


if __name__ == "__main__":
    main()
