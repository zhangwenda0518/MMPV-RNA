#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eve_screen.py — 内源性病毒元件 (EVE) 筛查编排器
================================================
对宿主参考基因组做三阶段 EVE 筛查 (移植自 hi-fever eve_kingdom.py):

  discover  (1) 正向发现: 滑窗 → diamond blastx vs 病毒参考蛋白 → 位点合并 → 抽序列
  verdict   (2) 双侧判定: vs 植物/病毒拆分库 → viral_supported / host_like / undetermined
  annotate  (3) RVDB 深度归属: viral_supported + undetermined vs RVDB
  merge          跨基因组汇总 (kingdom_summary / family_by_genome)

用法:
  # 单基因组
  python eve_screen.py -g genome.fa -n Solanum_lycopersicum -o out/ -t 40

  # 批量并行 (batch.tsv 每行 NAME<TAB>/path/to/genome.fa; NAME 必须唯一,
  # 重名会共用 01_Loci/<NAME>/ 产物, 开跑前直接报错中止)
  python eve_screen.py -B batch.tsv -o out/ -t 40 -J 5 --fast
  # -J 并行基因组数, -t 每基因组线程数; 总核数 = J x t (256 核推荐 -J5 -t40)
  # 批量路径可直接填 .tar.gz/.gz 基因组包 (Stage1 自动解包)

  # 断点续传: 以 04_Summary/<NAME>_eve_summary.tsv 为完成标记, 已完成基因组自动跳过
  python eve_screen.py -B batch.tsv -o out/ -t 40 -J 5

  # 汇总与干跑
  python eve_screen.py --merge -o out/
  python eve_screen.py -B batch.tsv --list

常用选项:
  --fast         Stage1 用 --fast (约快6倍, 位点数约-54%, 千种全扫推荐)
  --overlap      切块步长减半 (25kb 重叠, 修跨界截断; 重点基因组复扫用)
  --skip-s3      跳过 RVDB 层
  --blastn-viroid  附加类病毒 blastn 层 (默认关)
  --cleanup      每基因组完成后删除已被下游产物取代的中间文件
  --force        忽略断点全量重跑
  --no-merge     批量完成后不自动汇总 (默认自动 merge)
  --cmd-timeout N 单条外部命令超时秒数; 0=不限 (防 diamond/seqkit 挂死占核)

抽序列: samtools faidx (需 samtools >= 1.11, 提供 faidx 的 --region-file
  /--fai-idx/--output; 版本不够会在开跑前体检时报错); 可用 --samtools 指定
  可执行文件, 或写在 pipeline_config.yaml 的 tools: 块 (MMPV_SAMTOOLS).

数据库路径: CLI > pipeline_config.yaml (--profile) > 环境变量
  (MMPV_VIRUS_DB / MMPV_DB_ROOT / MMPV_EVE_PV_DB / MMPV_EVE_RVDB_DB)
"""

import argparse
import logging
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from eve_scan_core import (EveConfig, check_tools, clean_name, merge_all,
                           safe_dir)
from eve_genome_scan import parse_stages, run_genome, setup_logger

STAGE_HELP = {
    "discover":  "01: 正向发现 (滑窗 → blastx vs 病毒参考蛋白 → 位点合并 → 抽序列)",
    "verdict":   "02: 双侧判定 (viral_supported / host_like / undetermined)",
    "annotate":  "03: RVDB 深度归属 (viral_supported + undetermined)",
    "merge":     "04: 跨基因组汇总 (kingdom_summary / family_by_genome)",
}

# pipeline_config.yaml 读取键: (yaml 键, argparse dest)
DB_KEYS = [
    ("eve_ref_db", "ref_db"),
    ("eve_pv_db", "pv_db"),
    ("eve_id2div", "id2div"),
    ("eve_rvdb_db", "rvdb_db"),
    ("viroids_db", "viroids_db"),
]
EVE_PARAM_KEYS = ["window", "merge_distance", "evalue", "host_bs", "viral_bs",
                  "cmd_timeout"]
# CLI 长选项 -> args dest (判断 CLI 是否显式覆盖 yaml 的 eve: 参数)
_EVE_FLAGS = {"--window": "window", "--merge-distance": "merge_distance",
              "--evalue": "evalue", "--host-bs": "host_bs",
              "--viral-bs": "viral_bs", "--cmd-timeout": "cmd_timeout"}
# 运行参数 CLI -> dest (含短选项, 兼容 --flag value 与 --flag=value)
_RT_FLAGS = {"-t": "threads", "--threads": "threads",
             "-J": "jobs", "--jobs": "jobs"}


def _expand_env_config(obj):
    """递归展开 ${VAR} / ${VAR:-默认} (与 virome_pipeline.py 语义一致)."""
    pat = re.compile(r"\$\{([A-Za-z_]\w*)(?::-((?:[^{}]|\{[^{}]*\})*))?\}")

    def _s(v):
        def _sub(m):
            name, default = m.group(1), m.group(2)
            val = os.environ.get(name)
            if val:
                return val
            return default if default is not None else ""
        for _ in range(10):
            new = pat.sub(_sub, v)
            if new == v:
                break
            v = new
        return v

    if isinstance(obj, str):
        return _s(obj)
    if isinstance(obj, dict):
        return {k: _expand_env_config(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env_config(v) for v in obj]
    return obj


def _cli_flags(argv, flag_map):
    """扫描 argv, 返回 CLI 显式指定的 dest 集合 (--flag value / --flag=value)."""
    return {dest for a in argv for flag, dest in flag_map.items()
            if a == flag or a.startswith(flag + "=")}


def _load_config(args, argv=None):
    """加载 pipeline_config.yaml: CLI 未显式指定 > profile > 环境变量默认."""
    argv = list(sys.argv[1:]) if argv is None else list(argv)
    config_path = args.config
    if not config_path:
        for p in [SCRIPT_DIR.parent / "pipeline_config.yaml",
                  SCRIPT_DIR / "pipeline_config.yaml",
                  Path.cwd() / "pipeline_config.yaml"]:
            if p.is_file():
                config_path = str(p)
                break
    if not config_path:
        return
    try:
        import yaml
    except ImportError:
        print("[WARN] PyYAML 未安装, 跳过配置文件加载 (pip install pyyaml)")
        return
    with open(config_path, encoding="utf-8") as cf:
        config = _expand_env_config(yaml.safe_load(cf))
    profiles = config.get("profiles", {})
    profile = profiles.get(args.profile, profiles.get("default", {}))
    if not profile:
        print(f"[WARN] profile '{args.profile}' 未找到, 仅用 CLI 参数")
        return

    db = profile.get("databases", {})
    for ykey, field in DB_KEYS:
        cli_val = getattr(args, field, None)
        if not cli_val and db.get(ykey):
            setattr(args, field, db[ykey])
    tools = profile.get("tools", {})
    if not args.diamond and tools.get("diamond"):
        args.diamond = os.path.expanduser(tools["diamond"])
    if not args.samtools and tools.get("samtools"):
        args.samtools = os.path.expanduser(tools["samtools"])
    eve = profile.get("eve", {})
    # CLI 显式指定的 EVE 参数优先, 其余才用 yaml 值 (与 runtime 判定同思路)
    _eve_cli = _cli_flags(argv, _EVE_FLAGS)
    for key in EVE_PARAM_KEYS:
        if key in eve and key not in _eve_cli:
            setattr(args, key, eve[key])
    rt = profile.get("runtime", {})
    _cli = _cli_flags(argv, _RT_FLAGS)
    for key in ("threads", "jobs"):
        if key in rt and key not in _cli:
            setattr(args, key, int(rt[key]))


def build_parser():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="内源性病毒元件 (EVE) 三阶段筛查编排器")
    g = p.add_argument_group("输入")
    g.add_argument("-g", "--genome", help="单基因组 FASTA (.fa/.fasta/.fna/.gz/.tar.gz)")
    g.add_argument("-n", "--name", help="基因组名 (默认: 文件名去后缀)")
    g.add_argument("-B", "--batch", help="批量 TSV: 每行 NAME<TAB>基因组路径")
    g.add_argument("--list", action="store_true",
                   help="仅列出任务 (name\\tpath) 不执行")

    g = p.add_argument_group("运行")
    g.add_argument("-o", "--outdir", default="eve_results", help="输出根目录")
    g.add_argument("-t", "--threads", type=int, default=40)
    g.add_argument("-J", "--jobs", type=int, default=1,
                   help="并行基因组数 (总核数 = jobs x threads)")
    g.add_argument("--stage", default="all",
                   help="阶段: discover/verdict/annotate 逗号组合 或 all")
    g.add_argument("--merge", action="store_true",
                   help="只跑跨基因组汇总 (不筛查)")
    g.add_argument("--no-merge", action="store_true",
                   help="批量完成后不自动汇总")
    g.add_argument("--force", action="store_true",
                   help="忽略断点全量重跑 (已完成基因组也重跑)")
    g.add_argument("--fast", action="store_true",
                   help="Stage1 用 --fast (约快6倍, 位点数约-54%%)")
    g.add_argument("--overlap", action="store_true",
                   help="切块 25kb 重叠 (修跨界截断; 重点基因组复扫用)")
    g.add_argument("--skip-s3", action="store_true", help="跳过 RVDB 层")
    g.add_argument("--blastn-viroid", action="store_true",
                   help="附加类病毒 blastn 层 (默认关)")
    g.add_argument("--cleanup", action="store_true",
                   help="完成后删除已无用中间文件 (chunks/raw/hits)")
    g.add_argument("--config", default=None,
                   help="YAML 配置文件 (默认: 自动查找 pipeline_config.yaml)")
    g.add_argument("--profile", default="default",
                   help="配置预设 (default/plant, 可选见 pipeline_config.yaml)")

    g = p.add_argument_group("数据库 (默认: 从 pipeline_config.yaml 读取)")
    g.add_argument("--ref-db", default=None, help="Stage1 病毒参考蛋白 diamond 库")
    g.add_argument("--pv-db", default=None, help="Stage2 植物/病毒拆分库")
    g.add_argument("--id2div", default=None, help="Stage2 sseqid->viral/plant 表")
    g.add_argument("--rvdb-db", default=None, help="Stage3 RVDB diamond 库")
    g.add_argument("--viroids-db", default=None, help="类病毒 fasta (blastn 层)")

    g = p.add_argument_group("工具")
    g.add_argument("--diamond", default=None, help="diamond 可执行文件 (默认 PATH)")
    g.add_argument("--samtools", default=None,
                   help="samtools 可执行文件 (默认 PATH)")

    g = p.add_argument_group("EVE 参数")
    g.add_argument("--window", type=int, default=50000, help="滑窗宽度 (nt)")
    g.add_argument("--merge-distance", type=int, default=300,
                   help="位点合并距离 (nt)")
    g.add_argument("--evalue", default="1e-5")
    g.add_argument("--host-bs", type=float, default=50,
                   help="Stage2 植物侧判定阈值 (bitscore)")
    g.add_argument("--viral-bs", type=float, default=50,
                   help="Stage2 病毒侧判定阈值 (bitscore)")
    g.add_argument("--cmd-timeout", type=int, default=0,
                   help="单条外部命令超时 (秒); 0=不限 (防工具挂死占核)")
    return p


def build_jobs(args, stages, log):
    """从 -g / -B 构建任务列表."""
    spec_common = dict(threads=args.threads, stages=stages, fast=args.fast,
                       overlap=args.overlap, viroid=args.blastn_viroid,
                       skip_s3=args.skip_s3, cleanup=args.cleanup,
                       force=args.force)
    todo = []
    if args.batch:
        with open(args.batch, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh, 1):
                c = line.rstrip("\n").split("\t")
                if len(c) != 2 or c[0].startswith("#"):
                    continue
                todo.append(dict(spec_common, name=clean_name(c[0].strip()),
                                 genome=os.path.expanduser(c[1].strip()),
                                 _src=f"{args.batch}:{lineno}"))
    elif args.genome:
        name = args.name or Path(args.genome).name
        for suf in (".fna", ".fasta", ".fa", ".fasta.gz", ".fa.gz",
                    ".fna.gz", ".tar.gz", ".gz"):
            if name.endswith(suf):
                name = name[: -len(suf)]
                break
        todo.append(dict(spec_common, name=clean_name(name),
                         genome=str(Path(args.genome)), _src="-g"))
    else:
        return []
    # 重名检查: clean_name 后重名会共用 01_Loci/<NAME>/ 产物 (互相覆盖且
    # summary 对不上), 必须在开跑前挡住而不是静默覆盖
    seen = {}
    for t in todo:
        seen.setdefault(t["name"], []).append(t["_src"])
    dups = {n: srcs for n, srcs in seen.items() if len(srcs) > 1}
    if dups:
        for n, srcs in sorted(dups.items()):
            log.error("批次重名 '%s' (来自 %s) — 请改 TSV 里的名字", n,
                      ", ".join(srcs))
        raise SystemExit(f"批次中有 {len(dups)} 个重名, 已中止")
    outdir = safe_dir(args.outdir)
    cfg = EveConfig(ref_dmnd=args.ref_db or "", pv_dmnd=args.pv_db or "",
                    id2div=args.id2div or "", rvdb_dmnd=args.rvdb_db or "",
                    viroids_fa=args.viroids_db or "",
                    diamond=args.diamond or "", samtools=args.samtools or "",
                    window=args.window, merge_d=args.merge_distance,
                    evalue=args.evalue, host_bs=args.host_bs,
                    viral_bs=args.viral_bs, cmd_timeout=args.cmd_timeout)
    jobs = []
    for t in todo:
        jobs.append(dict(t, outdir=str(outdir), cfg=cfg,
                         log_path=str(Path(outdir) / "logs" / f"{t['name']}.log")))
    return jobs


def main(argv=None):
    args = build_parser().parse_args(argv)
    log = setup_logger(Path(safe_dir(args.outdir)) / "logs" / "eve_screen.log")
    _load_config(args, argv)

    if args.merge:
        n_g, n_l = merge_all(safe_dir(args.outdir))
        log.info("MERGE COMPLETE: %d genomes, %d loci -> %s", n_g, n_l,
                 safe_dir(args.outdir))
        return 0

    stages = parse_stages(args.stage)
    jobs = build_jobs(args, stages, log)
    if args.list:
        for j in jobs:
            print(f"{j['name']}\t{j['genome']}")
        return 0
    if not jobs:
        log.error("need -g/-B/--merge")
        return 1

    cfg = jobs[0]["cfg"]
    missing = check_tools(cfg, need_viroid=args.blastn_viroid)
    if missing:
        log.error("缺少工具/数据库: %s", "; ".join(missing))
        return 1

    log.info("批次: %d 个基因组, jobs=%d, threads=%d, stages=%s",
             len(jobs), args.jobs, args.threads, ",".join(sorted(stages)))
    t0 = time.time()
    results = []
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_genome, j): j["name"] for j in jobs}
            for i, fu in enumerate(as_completed(futs), 1):
                r = fu.result()
                results.append(r)
                _log_result(log, i, len(jobs), r)
    else:
        for i, j in enumerate(jobs, 1):
            r = run_genome(j)
            results.append(r)
            _log_result(log, i, len(jobs), r)

    n_err = sum(1 for r in results if r["status"] == "error")
    n_skip = sum(1 for r in results if r["status"] == "skip")
    log.info("BATCH COMPLETE: %d total, %d done, %d skipped, %d error "
             "(%.0fs)", len(results), len(results) - n_err - n_skip, n_skip,
             n_err, time.time() - t0)
    for r in results:
        if r["status"] == "error":
            log.error("  FAILED %s: %s", r["name"], r.get("error", ""))

    if not args.no_merge and stages == {"discover", "verdict", "annotate"}:
        n_g, n_l = merge_all(safe_dir(args.outdir))
        log.info("MERGE COMPLETE: %d genomes, %d loci", n_g, n_l)
    return 1 if n_err else 0


def _log_result(log, i, total, r):
    extra = ""
    if r["status"] == "done":
        extra = f" counts={r.get('counts', {})}"
    elif r["status"] == "error":
        extra = f" ERROR {r.get('error', '')}"
    log.info("[%d/%d] %s: %s (%ss)%s", i, total, r["name"], r["status"],
             r.get("elapsed"), extra)


if __name__ == "__main__":
    sys.exit(main())
