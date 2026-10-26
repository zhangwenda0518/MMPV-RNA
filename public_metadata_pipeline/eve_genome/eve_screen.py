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
import hashlib
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from eve_scan_core import (BATCH_MANIFEST_NAME, EVALUE, FLANK_MERGE_D, FLANK_NT,
                           HOST_BS, MERGE_D, VIRAL_BS, WINDOW, EveConfig,
                           check_tools, clean_name, merge_all, safe_dir)
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
    g.add_argument("-o", "--outdir", default="01_Scan",
                   help="输出根目录 (默认 01_Scan; 编排器传 out/eve/01_Scan)")
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
    g.add_argument("--allow-duplicate-assembly", action="store_true",
                   help="放行 assembly 指纹重复的批次登记 (默认中止; 放行也记录在案)")
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
    g.add_argument("--window", type=int, default=WINDOW, help="滑窗宽度 (nt)")
    g.add_argument("--merge-distance", type=int, default=MERGE_D,
                   help="位点合并距离 (nt)")
    g.add_argument("--evalue", default=EVALUE)
    g.add_argument("--host-bs", type=float, default=HOST_BS,
                   help="Stage2 植物侧判定阈值 (bitscore)")
    g.add_argument("--viral-bs", type=float, default=VIRAL_BS,
                   help="Stage2 病毒侧判定阈值 (bitscore)")
    g.add_argument("--cmd-timeout", type=int, default=0,
                   help="单条外部命令超时 (秒); 0=不限 (防工具挂死占核)")
    g.add_argument("--flank-nt", type=int, default=FLANK_NT,
                   help="侧翼上下文长度 (nt, 两侧各扩; 0=不抽)")
    g.add_argument("--flank-merge-d", type=int, default=FLANK_MERGE_D,
                   help="抽侧翼前先合并位点的距离 (nt)")
    return p


def assembly_fingerprint(path, tail=65536):
    """基因组文件的 assembly 指纹 (P1-5b): sha256(头 64KB + 尾 64KB), 大文件不读中间.

    目的只有一个: 抓"同一文件在批次里换名登记两次" —— 这种重复会把元件聚类里的
    跨基因组支持数 (n_genomes) 抬高 (审查 R02)。小文件 (<=64KB) 整文件哈希
    (mode=full); 大文件只读头尾 (mode=head_tail_64k), 上 TB 级批量不能在开跑前
    做全量哈希。**保证**: 同一字节序列必得同一指纹; **不承诺**: 同一基因组重新
    打包/重压缩后字节变了就抓不到 (那要解开包才算得出)。文件读不到返回
    ("", "missing") —— 输入存在性由 Stage1 负责, 这里不做第二套校验。
    """
    p = Path(path)
    try:
        size = p.stat().st_size
    except OSError:
        return "", "missing"
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        head = fh.read(tail)
        h.update(head)
        if size > len(head):
            fh.seek(max(0, size - tail))
            h.update(fh.read(tail))
    return h.hexdigest(), ("full" if size <= tail else "head_tail_64k")


def write_batch_manifest(outdir, todo):
    """批次清单 (溯源 + 元件层 n_assemblies 去重的依据) -> <outdir>/batch_manifest.tsv."""
    p = Path(outdir) / BATCH_MANIFEST_NAME
    p.parent.mkdir(parents=True, exist_ok=True)   # safe_dir 只归一化不建目录
    with open(p, "w", encoding="utf-8", newline="\n") as fo:
        fo.write("name\tgenome_path\tsize_bytes\tfingerprint\tfingerprint_mode\n")
        for t in sorted(todo, key=lambda x: x["name"]):
            try:
                size = os.path.getsize(t["genome"])
            except OSError:
                size = ""
            fo.write(f"{t['name']}\t{t['genome']}\t{size}\t"
                     f"{t.get('fingerprint', '')}\t{t.get('fp_mode', '')}\n")
    return p


def build_jobs(args, stages, log, fingerprint=True):
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
    # P1-5b: assembly 指纹去重 —— 同一文件换名登记两次会把元件层的跨基因组
    # 支持数抬高, 且没有任何其他机制能发现它。指纹重复默认中止, 显式
    # --allow-duplicate-assembly 才放行 (放行也记录在案, n_assemblies 列可见)。
    n_unique_asm = len(todo)
    if fingerprint:
        for t in todo:
            t["fingerprint"], t["fp_mode"] = assembly_fingerprint(t["genome"])
        by_fp = {}
        for t in todo:
            if t["fingerprint"]:
                by_fp.setdefault(t["fingerprint"], []).append(t)
        dup_asm = {fp: ts for fp, ts in by_fp.items()
                   if len({x["name"] for x in ts}) > 1}
        if dup_asm:
            for fp, ts in sorted(dup_asm.items()):
                names = ", ".join(f"{x['name']} <- {x['genome']}" for x in ts)
                log.error("批次 assembly 指纹重复 (%s…): %s",
                          fp[:12], names)
            if not getattr(args, "allow_duplicate_assembly", False):
                raise SystemExit(
                    f"批次中有 {len(dup_asm)} 组重复 assembly 登记 (同一文件多个名字), "
                    f"已中止; 确属有意请加 --allow-duplicate-assembly")
        n_unique_asm = len({t["fingerprint"] or f"name:{t['name']}" for t in todo})
        log.info("assembly 指纹: %d 个登记名 / %d 个唯一 assembly", len(todo),
                 n_unique_asm)
    outdir = safe_dir(args.outdir)
    if fingerprint:
        write_batch_manifest(outdir, todo)
    cfg = EveConfig(ref_dmnd=args.ref_db or "", pv_dmnd=args.pv_db or "",
                    id2div=args.id2div or "", rvdb_dmnd=args.rvdb_db or "",
                    viroids_fa=args.viroids_db or "",
                    diamond=args.diamond or "", samtools=args.samtools or "",
                    window=args.window, merge_d=args.merge_distance,
                    evalue=args.evalue, host_bs=args.host_bs,
                    viral_bs=args.viral_bs, cmd_timeout=args.cmd_timeout,
                    flank_nt=args.flank_nt, flank_merge_d=args.flank_merge_d)
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
    # --list 是干跑: 不为它做指纹 (可能上 TB 的哈希只为打印一行名单)
    jobs = build_jobs(args, stages, log, fingerprint=not args.list)
    if args.list:
        for j in jobs:
            print(f"{j['name']}\t{j['genome']}")
        return 0
    if not jobs:
        log.error("need -g/-B/--merge")
        return 1

    cfg = jobs[0]["cfg"]
    missing = check_tools(cfg, need_viroid=args.blastn_viroid,
                          skip_s3=args.skip_s3)
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
                # worker 硬死 (OOM killer / 段错误) 会让 fu.result() 抛
                # BrokenProcessPool; 那是整个进程池的状态, 从此所有 future 都取不到
                # 结果. 逐 future 兜住, 至少把它记成 error 行而不是让整批消失.
                try:
                    r = fu.result()
                except Exception as exc:
                    name = futs[fu]
                    r = {"name": name, "status": "error", "elapsed": 0,
                         "error": f"worker 死亡: {type(exc).__name__}"[:300]}
                    log.error("worker 死亡, 剩余 future 可能都取不到结果: %s: %s",
                              name, type(exc).__name__)
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
