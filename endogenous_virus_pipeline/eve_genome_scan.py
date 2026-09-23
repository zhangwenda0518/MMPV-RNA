#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eve_genome_scan.py — 单基因组 EVE 筛查 worker (断点续传)
=========================================================
对单个宿主参考基因组跑 Stage1→2→3 + 汇总. 由 eve_screen.py 编排器并行调用,
也可独立运行 (单基因组重跑/续跑, 不必启动整个批次):

  python eve_genome_scan.py -g genome.fa -n NAME -o OUTDIR -t 40 \
      --ref-db virus_ref.pep.dmnd --pv-db plant_virus.dmnd \
      --id2div id2div.tsv --rvdb-db RVDB.dmnd [--fast] [--blastn-viroid] \
      [--stages 1,2,3] [--force] [--cleanup]

断点续传语义: 阶段产物文件存在即视为完成 (空文件 = 已完成且零结果),
重发同一命令即从断点继续; --force 删除旧产物全量重跑.
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from eve_scan_core import (EveConfig, STAGE_DIRS, check_tools, cleanup_stage_files,
                           logmsg, safe_dir, stage1_discover, stage2_verdict,
                           stage3_annotate, stage_viroid, summarize_genome)

SUMMARY_SUFFIX = "_eve_summary.tsv"


def setup_logger(log_path=None, level="INFO"):
    logger = logging.getLogger("EVEScan")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, level.upper(), logging.INFO))
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_path, encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def parse_stages(spec):
    """'1,2,3' / 'discover,verdict' -> {'discover','verdict','annotate'} 子集."""
    alias = {"1": "discover", "2": "verdict", "3": "annotate",
             "discover": "discover", "verdict": "verdict", "annotate": "annotate"}
    out = set()
    for tok in str(spec).split(","):
        tok = tok.strip()
        if tok in ("all", ""):
            return {"discover", "verdict", "annotate"}
        if tok not in alias:
            raise SystemExit(f"未知阶段: {tok} (可选: 1/discover, 2/verdict, "
                             f"3/annotate, all)")
        out.add(alias[tok])
    return out


def genome_done(outdir, name):
    return (Path(outdir) / STAGE_DIRS["summary"] / f"{name}{SUMMARY_SUFFIX}").exists()


def resolve_genome_fa(genome, outdir, name, stages):
    """跳过 Stage1 时推断可读的基因组 fasta.

    discover 不在 stages 时 Stage1 的解压产物不存在, 而 samtools/diamond
    读不了 .tar.gz/.gz — 压缩输入必须改用 Stage1 落盘的
    01_Loci/<NAME>/<NAME>.fna (须先跑过 --stage discover 解包).
    """
    if "discover" in stages:
        return None
    g = Path(genome)
    unpacked = (Path(outdir) / STAGE_DIRS["loci"] / name / f"{name}.fna")
    if str(genome).endswith((".tar.gz", ".gz")):
        if unpacked.is_file():
            return unpacked
        raise FileNotFoundError(
            f"输入是压缩/打包基因组 ({genome}), 但未找到 Stage1 解压产物 "
            f"{unpacked}; 请先跑 --stage discover 解包, 或改用裸 fasta 输入")
    if g.is_file():
        return g
    if unpacked.is_file():
        return unpacked
    return g      # 原路径不存在: 交给后续 is_file 检查报统一错误


def run_genome(job):
    """进程池/串行 worker: 单基因组全流程 (断点续传 + 故障隔离).

    job 键: name, genome, outdir, threads, stages(set), fast, overlap,
            viroid, skip_s3, cleanup, force, cfg(EveConfig), log_path
    返回 dict: name/status/elapsed/counts[/error]
    """
    name = job["name"]
    t0 = time.time()
    outdir = safe_dir(job["outdir"])
    cfg = job["cfg"]
    log = setup_logger(job.get("log_path"))
    log.info("=== %s 开始 (stages=%s) ===", name, ",".join(sorted(job["stages"])))
    try:
        if not job["force"] and genome_done(outdir, name):
            log.info("已完成 (summary 存在), 跳过; --force 可重跑")
            return {"name": name, "status": "skip", "elapsed": 0}

        stages = job["stages"]
        genome_fa = None
        if "discover" in stages:
            genome_fa, _ = stage1_discover(cfg, job["genome"], name, outdir,
                                           job["threads"], job["fast"],
                                           job["overlap"], job["viroid"],
                                           log, force=job["force"])
        else:
            # 跳过 Stage1: 原路径优先, 压缩/打包输入改用 Stage1 解压产物
            genome_fa = resolve_genome_fa(job["genome"], outdir, name, stages)
        if "verdict" in stages:
            stage2_verdict(cfg, name, outdir, job["threads"], log,
                           force=job["force"])
        if "annotate" in stages and not job.get("skip_s3"):
            if genome_fa is None or not Path(genome_fa).is_file():
                raise FileNotFoundError(f"Stage3 需要基因组 fasta: {genome_fa}")
            stage3_annotate(cfg, name, outdir, job["threads"], genome_fa,
                            log, force=job["force"])
        if job["viroid"] and genome_fa and Path(genome_fa).is_file():
            stage_viroid(cfg, name, outdir, genome_fa, job["threads"], log,
                         force=job["force"])

        summ = summarize_genome(name, outdir)
        freed = cleanup_stage_files(outdir, name) if job.get("cleanup") else 0

        counts = {}
        vp = (Path(outdir) / STAGE_DIRS["verdict"] / name
              / f"{name}.s2_verdict.tsv")
        if vp.exists():
            with open(vp, encoding="utf-8") as fh:
                fh.readline()
                for line in fh:
                    parts = line.rstrip("\n").split("\t", 2)
                    if len(parts) > 1:
                        counts[parts[1]] = counts.get(parts[1], 0) + 1
        log.info("=== %s 完成: %s (freed %.1f GB) ===", name, summ,
                 freed / 1073741824)
        return {"name": name, "status": "done",
                "elapsed": round(time.time() - t0), "counts": counts,
                "freed_gb": round(freed / 1073741824, 1)}
    except Exception as exc:
        log.error("=== %s 失败: %s: %s ===", name, type(exc).__name__, exc)
        return {"name": name, "status": "error",
                "elapsed": round(time.time() - t0),
                "error": f"{type(exc).__name__}: {exc}"[:300]}


def build_parser():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="单基因组 EVE 三阶段筛查 worker")
    p.add_argument("-g", "--genome", required=True,
                   help="宿主基因组 FASTA (支持 .fa/.fasta/.fna/.gz/.tar.gz)")
    p.add_argument("-n", "--name", default=None,
                   help="基因组名 (默认: 文件名去后缀)")
    p.add_argument("-o", "--outdir", required=True, help="输出根目录")
    p.add_argument("-t", "--threads", type=int, default=40)
    p.add_argument("--stage", "--stages", default="1,2,3",
                   help="阶段: 1/discover,2/verdict,3/annotate 逗号组合 或 all")
    p.add_argument("--fast", action="store_true",
                   help="Stage1 用 --fast (约快6倍, 位点数约-54%%)")
    p.add_argument("--overlap", action="store_true",
                   help="切块步长减半 (25kb 重叠, 修跨界截断; 重点基因组复扫用)")
    p.add_argument("--blastn-viroid", action="store_true",
                   help="附加类病毒 blastn 层 (默认关)")
    p.add_argument("--skip-s3", action="store_true", help="跳过 RVDB 层")
    p.add_argument("--cleanup", action="store_true",
                   help="完成后删除已被下游产物取代的中间文件")
    p.add_argument("--force", action="store_true", help="忽略断点全量重跑")
    p.add_argument("--log", default=None, help="日志文件路径 (默认仅 stderr)")
    # 数据库
    g = p.add_argument_group("数据库 (单跑时用 CLI 指定; 批量跑由编排器按 "
                             "pipeline_config.yaml 注入)")
    g.add_argument("--ref-db", default=None, help="Stage1 病毒参考蛋白 diamond 库")
    g.add_argument("--pv-db", default=None, help="Stage2 植物/病毒拆分库")
    g.add_argument("--id2div", default=None, help="Stage2 sseqid->viral/plant 表")
    g.add_argument("--rvdb-db", default=None, help="Stage3 RVDB diamond 库")
    g.add_argument("--viroids-db", default=None, help="类病毒 fasta (blastn 层)")
    # 工具
    g = p.add_argument_group("工具")
    g.add_argument("--diamond", default=None, help="diamond 可执行文件 (默认 PATH)")
    g.add_argument("--samtools", default=None,
                   help="samtools 可执行文件 (默认 PATH)")
    # 参数
    g = p.add_argument_group("运行参数")
    g.add_argument("--window", type=int, default=50000)
    g.add_argument("--merge-distance", type=int, default=300)
    g.add_argument("--evalue", default="1e-5")
    g.add_argument("--host-bs", type=float, default=50)
    g.add_argument("--viral-bs", type=float, default=50)
    g.add_argument("--cmd-timeout", type=int, default=0,
                   help="单条外部命令超时 (秒); 0=不限 (防工具挂死占核)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    log = setup_logger(args.log)
    name = args.name or Path(args.genome).name
    for suf in (".fna", ".fasta", ".fa", ".fasta.gz", ".fa.gz", ".fna.gz",
                ".tar.gz", ".gz"):
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    from eve_scan_core import clean_name
    cfg = EveConfig(ref_dmnd=args.ref_db or "", pv_dmnd=args.pv_db or "",
                    id2div=args.id2div or "", rvdb_dmnd=args.rvdb_db or "",
                    viroids_fa=args.viroids_db or "", diamond=args.diamond or "",
                    samtools=args.samtools or "", window=args.window,
                    merge_d=args.merge_distance, evalue=args.evalue,
                    host_bs=args.host_bs, viral_bs=args.viral_bs,
                    cmd_timeout=args.cmd_timeout)
    missing = check_tools(cfg, need_viroid=args.blastn_viroid)
    if missing:
        log.error("缺少工具/数据库: %s", "; ".join(missing))
        return 1
    job = {"name": clean_name(name), "genome": str(Path(args.genome)),
           "outdir": str(safe_dir(args.outdir)), "threads": args.threads,
           "stages": parse_stages(args.stage), "fast": args.fast,
           "overlap": args.overlap, "viroid": args.blastn_viroid,
           "skip_s3": args.skip_s3, "cleanup": args.cleanup,
           "force": args.force, "cfg": cfg, "log_path": args.log}
    r = run_genome(job)
    if r["status"] == "error":
        log.error("FAILED %s: %s", r["name"], r.get("error", ""))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
