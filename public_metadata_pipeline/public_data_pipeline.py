#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Public Data Pipeline - 公共数据获取管道 v3.2
=============================================
九阶段: search -> info -> prescreen -> down -> convert -> plot -> hostref -> hostdb -> report
全部参数命令行指定。宿主参考基因组默认复用 --species/--taxid (检索物种即宿主物种)。
stage 直接调用工作脚本 (download_host_genome.py / build_hostbase.py), 不嵌套子管道。

prescreen 阶段 (v3.2 新增)
--------------------------
调 stat_prescreen.py, 用 NCBI STAT / CNCB GSA 的**预计算分类结果**做两件事,
**不下载任何测序数据**:

  1) 植物物种鉴定 —— 用实测主要物种对拍元信息 ScientificName, 揪出记错的样本
     (避免元信息错误一路带到下游)
  2) 病毒预鉴定   —— 按病毒占比给出 decision, 供 down 阶段过滤
     (避免把无病毒样本也下载下来, 浪费带宽与磁盘)

默认**不改变** down 的行为 (--prescreen-filter off); 显式指定才启用过滤。

prescreen 之后的**病毒层面画像** (v3.2, 默认开)
-----------------------------------------------
样本级问题回答完, 病毒本身也是结果。prescreen 跑完自动调 stat_virus_profile.py,
把 run×taxon 明细上卷成病毒层面的画像, 落到 <prescreen>/virus_profile/:

  * 有哪些病毒 —— 病毒名录 + 阶元 + 大类 (viroid/virus/phage/...)
  * 各占多少   —— 样本内丰度 (mean/median/max)
  * 多少样本里出现 —— 跨样本检出率 (prevalence), 分 leaf(末端节点) / all 两套口径
  * 图表       —— 检出率条形图 / 检出率-丰度气泡图 / 样本×病毒热图 / 大类构成 / 组合图

纯附加产物, 不改 sample_decision.tsv, 也不影响 down。
用 --no-virus-profile 关闭, --no-virus-profile-plot 只出表不绘图。

用法:
  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage search info plot

  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage all

  # 启用预筛过滤 (safe: 保留 DOWNLOAD+REVIEW+LOCAL_CHECK, 不漏未知样本)
  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --stage prescreen down --prescreen-filter safe
"""

import os
import sys
import time
import csv
import json
import re
import argparse
import shlex
from typing import Dict, List, Optional, Tuple
from pathlib import Path

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import (build_meta_dirs, normalize_layout_env,
                                   dir_name as layout_dir_name)
from mmpv_common.stage_args import parse_multi_stage

from utils.pipeline_utils import UI, Checkpoint, run_cmd

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


# ==========================================
# Public Data Pipeline
# ==========================================
class PublicDataPipeline:
    STAGES = ['search', 'info', 'prescreen', 'down', 'convert',
              'plot', 'hostref', 'hostdb', 'eve', 'report']

    STAGE_DESC = {
        'search':    'GSA + SRA 双引擎物种检索',
        'info':      '元数据深度解析 & 文献溯源',
        'prescreen': '病毒预筛 & 物种鉴定 (预计算分类, 不下载数据)',
        'down':      '高通量数据下载 (NGDC + NCBI)',
        'convert':   'SRA → FASTQ.GZ 转换',
        'plot':      'SCI 级六图可视化',
        'hostref':   '宿主基因组下载 (datasets/ngd/FTP/gget 四通道; --host-fasta 可跳过)',
        'hostdb':    '宿主索引构建 (Kraken2/Bowtie2/HISAT2/Minimap2)',
        'eve':       '基因组 EVE 分析 (宿主库构建后; 默认关, --eve-scan 开启)',
        'report':    'HTML 报告 + sample_handoff.csv 交接清单',
    }

    # prescreen 过滤模式: 决定 down 阶段保留哪些 decision 的样本
    PRESCREEN_KEEP = {
        'safe':   {'DOWNLOAD', 'REVIEW', 'LOCAL_CHECK'},  # 不漏「无分类数据」的未知样本
        'strict': {'DOWNLOAD', 'REVIEW'},                 # 只下确有病毒信号的
    }

    # hostdb 建库默认参数 (与 build_host_pipeline.py 历史行为对齐)
    HOST_K2_LIBS = 'archaea,bacteria,plasmid,fungi,protozoa,UniVec'
    HOST_SEQ_TYPES = 'dna-short,rna-short,nanopore,pacbio'

    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.bin_dir = args.bin_dir or os.path.dirname(os.path.abspath(__file__))
        self.work_dir = os.path.abspath(args.work_dir)
        self.ckpt = Checkpoint(self.work_dir)
        self.log_dir = os.path.join(self.work_dir, 'logs')
        # 布局解析 (CLI > MMPV_IO_LAYOUT > standard) 并写回环境供子进程继承
        self.layout = normalize_layout_env(getattr(args, 'io_layout', None))
        self.dirs = build_meta_dirs(Path(self.work_dir), self.layout)
        for _name, _path in self.dirs.items():
            os.makedirs(_path, exist_ok=True)
        # AI 溯源: 若 work_dir 之上能找到项目 90_Handoff, 让 search/info 子进程把
        # AI 调用记进 <项目>/90_Handoff/state/ai_log.jsonl (与平台 AI 层同账本)
        self.ai_state_dir = self._find_ai_state_dir()
        if self.ai_state_dir:
            os.environ["MMPV_AI_STATE_DIR"] = self.ai_state_dir

    def _find_ai_state_dir(self) -> Optional[str]:
        p = Path(self.work_dir)
        for cand in (p, *list(p.parents)[:4]):
            if (cand / "90_Handoff").is_dir():
                return str(cand / "90_Handoff" / "state")
        return None

    def _bin(self, script: str) -> str:
        p = os.path.join(self.bin_dir, script)
        if not os.path.isfile(p) and not script.endswith('.py'):
            p = os.path.join(self.bin_dir, script + '.py')
        return p

    def _py(self) -> str:
        # 裸 `python` 可能落到无依赖的解释器, 统一用当前解释器 (与 plot/report 阶段一致)
        return shlex.quote(sys.executable or 'python')

    def _secrets(self) -> List[str]:
        """Return list of API key strings to mask in logs."""
        s = []
        if self.args.deepseek_api:
            s.append(self.args.deepseek_api)
        if self.args.ncbi_api:
            s.append(self.args.ncbi_api)
        return s

    def _api_cli(self) -> str:
        parts = []
        if self.args.deepseek_api:
            parts.append(f'--deepseek-api {shlex.quote(self.args.deepseek_api)}')
            parts.append(f'--deepseek-model {shlex.quote(self.args.deepseek_model)}')
        if self.args.ncbi_api:
            parts.append(f'--ncbi-api {shlex.quote(self.args.ncbi_api)}')
        return ' '.join(parts)

    def _extract_runs(self) -> Tuple[Optional[str], int]:
        merged = os.path.join(self.dirs['search'], 'SRA_GSA_Merged_Final.csv')
        sra_list = os.path.join(self.dirs['info'], 'sra.list')
        if os.path.isfile(merged):
            runs = []
            with open(merged, 'r', encoding='utf-8-sig') as f:
                for row in csv.DictReader(f):
                    r = row.get('Run', '').strip()
                    if r and r != 'Run':
                        runs.append(r)
            with open(sra_list, 'w') as f:
                f.write('\n'.join(runs))
            return sra_list, len(runs)
        return None, 0

    def run_search(self) -> bool:
        s = 'search'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        detailed = '--detailed' if self.args.detailed else ''
        search_ai = ''
        if self.args.deepseek_search:
            search_ai = '--deepseek-search' + (' --deepseek-fetch' if self.args.deepseek_fetch else '')
        cmd = (
            f'{self._py()} {shlex.quote(self._bin("gsa_sra.search.py"))} '
            f'--query {shlex.quote(self.args.species)} '
            f'--source {shlex.quote(self.args.source_type)} '
            f'{detailed} {self._api_cli()} {search_ai} '
            f'--outdir {shlex.quote(self.dirs["search"])}'
        )
        rc = run_cmd(cmd, s, self.log_dir, secrets=self._secrets())
        merged = os.path.join(self.dirs['search'], 'SRA_GSA_Merged_Final.csv')
        if rc == 0 and os.path.isfile(merged):
            UI.ok(f"output: {merged}")
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

    def run_info(self) -> bool:
        s = 'info'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        sra_list, n = self._extract_runs()
        if not sra_list:
            UI.err("missing search output, run search first")
            self.ckpt.mark_fail(s)
            return False
        UI.ok(f"extracted {n} runs -> {sra_list}")
        mode = 'both' if self.args.deepseek_api else 'local'
        if mode == 'local':
            UI.warn("未提供 --deepseek-api, 将跳过 AI 深度清洗 (mode=local)")
        # search 阶段的 --deepseek-search 检索回填只落在 Merged CSV, info 重抓元数据,
        # 不显式贯通就会丢掉 —— 这里把表交给 info (仅补空值, 已有值不动)
        merged = os.path.join(self.dirs['search'], 'SRA_GSA_Merged_Final.csv')
        backfill_cli = ''
        if os.path.isfile(merged):
            backfill_cli = f'--search-backfill {shlex.quote(merged)}'
            UI.ok("检索回填表并入 Core14 (--search-backfill; 仅补空值)")
        cmd = (
            f'{self._py()} {shlex.quote(self._bin("gsa_sra.info.py"))} '
            f'--input {shlex.quote(sra_list)} --mode {shlex.quote(mode)} --fill-date '
            f'{self._api_cli()} {backfill_cli} --threads {self.args.threads} '
            f'--outdir {shlex.quote(self.dirs["info"])}'
        )
        rc = run_cmd(cmd, s, self.log_dir, secrets=self._secrets())
        if rc == 0:
            UI.ok(f"output: {self.dirs['info']}/Global_Unified_Metadata_Core14.csv")
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

    def _prescreen_filtered_list(self, mode: str) -> Optional[str]:
        """按 sample_decision.tsv 的 decision 列生成过滤后的 run 列表。

        mode: safe   -> 保留 DOWNLOAD + REVIEW + LOCAL_CHECK (不漏未知样本)
              strict -> 只保留 DOWNLOAD + REVIEW (最大化省资源)

        返回 None 表示无法过滤 (预筛产物不存在或 mode 非法),
        调用方应回退到完整列表 —— 绝不因为过滤失败而少下样本。
        """
        dec = os.path.join(self.dirs['prescreen'], 'sample_decision.tsv')
        if not os.path.isfile(dec):
            UI.warn(f"未找到预筛主表 {dec} —— 过滤跳过, 将下载全部样本")
            UI.warn("  (需要先运行 prescreen 阶段; 或去掉 --prescreen-filter)")
            return None
        keep = self.PRESCREEN_KEEP.get(mode)
        if not keep:
            UI.err(f"未知 --prescreen-filter 模式: {mode}")
            return None

        runs, tally = [], {}
        with open(dec, 'r', encoding='utf-8') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                d = (row.get('decision') or '').strip()
                tally[d] = tally.get(d, 0) + 1
                if d in keep:
                    runs.append((row.get('acc') or '').strip())
        runs = [r for r in runs if r]

        out = os.path.join(self.dirs['prescreen'], f'prescreen_filtered_{mode}.list')
        with open(out, 'w', encoding='utf-8') as f:
            f.write('\n'.join(sorted(runs)) + ('\n' if runs else ''))
        UI.ok(f"预筛过滤 [{mode}]: 保留 {len(runs)} runs -> {out}")
        UI.info("  决策分布: " + ", ".join(f"{k}={v}" for k, v in sorted(tally.items())))
        dropped = sum(v for k, v in tally.items() if k not in keep)
        if dropped:
            UI.warn(f"  {dropped} 个样本按预筛结论过滤掉, 不会下载")
        return out

    def run_prescreen(self) -> bool:
        """病毒预筛 & 物种鉴定: 调 stat_prescreen.py, **不下载任何测序数据**。

        输入: info 产出的 sra.list + Global_Unified_Metadata_Core14.csv
        产出: sample_decision.tsv (主表) / prescreen_run.list (下载清单) / ...
        产物落在 self.dirs['prescreen'], 供 down 阶段 (--prescreen-filter) 消费。
        """
        s = 'prescreen'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)

        sra_list = os.path.join(self.dirs['info'], 'sra.list')
        if not os.path.isfile(sra_list):
            # 与 run_down 一致: 缺 sra.list 时从 search 产物重建
            _, n = self._extract_runs()
            sra_list = os.path.join(self.dirs['info'], 'sra.list')
            if n:
                UI.info(f"从 search 产物重建 sra.list ({n} runs)")
        if not os.path.isfile(sra_list):
            UI.err(f"missing sra list: {sra_list} (先运行 search/info 阶段)")
            self.ckpt.mark_fail(s, "missing sra.list")
            return False

        outdir = self.dirs['prescreen']
        os.makedirs(outdir, exist_ok=True)

        # 物种对拍基准: 优先 Core14 的逐 run ScientificName, 退回 --species
        core_csv = os.path.join(self.dirs['info'], 'Global_Unified_Metadata_Core14.csv')
        core_arg = f'--core-csv {shlex.quote(core_csv)}' if os.path.isfile(core_csv) else ''
        if not core_arg:
            UI.warn("未找到 Core14 元数据 —— 物种对拍退化为全局 --species 基准")

        thr_arg = ''
        if self.args.prescreen_virus_pct_ident is not None:
            thr_arg = f'--virus-pct-ident {self.args.prescreen_virus_pct_ident}'
        elif self.args.prescreen_virus_pct is not None:
            thr_arg = f'--virus-pct {self.args.prescreen_virus_pct}'

        taxdump_arg = f'--taxdump {shlex.quote(self.args.taxdump)}' if self.args.taxdump else ''
        gsa_arg = '--no-gsa' if self.args.prescreen_no_gsa else ''

        cmd = (
            f'{self._py()} {shlex.quote(self._bin("stat_prescreen.py"))} '
            f'--runs {shlex.quote(sra_list)} '
            f'{core_arg} '
            f'--expect-species {shlex.quote(self.args.species)} '
            f'{thr_arg} {taxdump_arg} {gsa_arg} '
            f'-j {self.args.prescreen_jobs} '
            f'--gsa-jobs {self.args.prescreen_gsa_jobs} '
            f'-o {shlex.quote(outdir)}'
        )
        UI.info("预筛只查上游预计算分类结果, 不下载测序数据")
        rc = run_cmd(cmd, s, self.log_dir, timeout=7200)
        dec = os.path.join(outdir, 'sample_decision.tsv')
        if rc == 0 and os.path.isfile(dec):
            UI.ok(f"output: {dec}")
            self._run_virus_profile(outdir)
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

    def _run_virus_profile(self, prescreen_dir: str) -> bool:
        """prescreen 之后跑**病毒层面画像** (名录 / 丰度 / 跨样本检出率 + 图)。

        病毒本身也是结果: 项目里有哪些病毒、各占多少、在多少样本里出现,
        与"这条 run 要不要下载"是两件独立的事。本步骤纯附加 ——
        只往 <prescreen>/virus_profile/ 写新文件, 不改 sample_decision.tsv,
        也不影响 down 阶段; 绘图依赖缺失时自动降级为只出表。
        """
        if not getattr(self.args, 'virus_profile', True):
            return True
        taxa = os.path.join(prescreen_dir, 'stat_viral_taxa.tsv')
        if not os.path.isfile(taxa):
            UI.warn("  未找到 stat_viral_taxa.tsv, 跳过病毒层面画像")
            return True

        outdir = os.path.join(prescreen_dir, 'virus_profile')
        plot_arg = '--plot' if getattr(self.args, 'virus_profile_plot', True) else ''
        if plot_arg:
            import importlib.util as _ilu
            missing = [m for m in ('numpy', 'pandas', 'matplotlib', 'seaborn')
                       if _ilu.find_spec(m) is None]
            if missing:
                UI.warn(f"  本解释器缺少 {'/'.join(missing)} —— 画像只出表不出图; "
                        f"要图请换用装了 matplotlib/seaborn 的解释器运行本管线")
        # 用**当前解释器**跑, 不用裸 `python`: run_cmd 经 shell 启动, Windows 下
        # `python` 会按 PATH 解析到另一个解释器 (实测父进程 C:\Python312, 子进程
        # 却落到无 matplotlib 的 3.13.12), 导致画像静默出不了图。用 sys.executable
        # 可保证子进程与父进程能力一致, 上面那次依赖探测也才代表子进程的真实情况。
        py = shlex.quote(sys.executable or 'python')
        cmd = (
            f'{py} {shlex.quote(self._bin("stat_virus_profile.py"))} '
            f'-i {shlex.quote(prescreen_dir)} '
            f'-o {shlex.quote(outdir)} '
            f'--min-samples {self.args.virus_profile_min_samples} '
            f'--top-n {self.args.virus_profile_top_n} '
            f'--prevalence-by {self.args.virus_profile_by} '
            f'{plot_arg}'
        )
        UI.info("病毒层面画像: 有哪些病毒 / 各占多少 / 在多少样本里出现")
        rc = run_cmd(cmd, 'prescreen', self.log_dir, timeout=1800)
        if rc == 0:
            UI.ok(f"output: {os.path.join(outdir, 'virus_prevalence.tsv')}")
        else:
            UI.warn(f"  病毒层面画像失败 (exit={rc}), 不影响 prescreen 主产物")
        return True

    def run_down(self) -> bool:
        s = 'down'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        sra_list = os.path.join(self.dirs['info'], 'sra.list')
        if not os.path.isfile(sra_list):
            _, _ = self._extract_runs()
            sra_list = os.path.join(self.dirs['info'], 'sra.list')
        if not os.path.isfile(sra_list):
            UI.err(f"missing sra list: {sra_list}")
            self.ckpt.mark_fail(s)
            return False

        # 预筛过滤: 默认 off, 不改变既有行为; 显式指定才按 sample_decision.tsv 筛
        if self.args.prescreen_filter != 'off':
            filtered = self._prescreen_filtered_list(self.args.prescreen_filter)
            if filtered:
                sra_list = filtered

        skip_arg = f'--skip-list {shlex.quote(self.args.skip_list)}' if self.args.skip_list else ''
        links_file = os.path.join(self.dirs['search'], 'download_links.csv')
        links_arg = f'--links {shlex.quote(links_file)}' if os.path.isfile(links_file) else ''
        cmd = (
            f'{self._py()} {shlex.quote(self._bin("gsa_sra.down.py"))} '
            f'--list {shlex.quote(sra_list)} '
            f'--ngdc-method {shlex.quote(self.args.ngdc_method)} '
            f'--ngdc-concurrency {self.args.ngdc_concurrency} '
            f'--prefetch-concurrency {self.args.prefetch_concurrency} '
            f'{skip_arg} '
            f'{links_arg} '
            f'--output {shlex.quote(self.dirs["down"])}'
        )
        UI.warn("data download may take hours to days")
        rc = run_cmd(cmd, s, self.log_dir, timeout=86400 * 7)
        if rc == 0:
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

    def run_convert(self) -> bool:
        s = 'convert'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)

        down_dir = self.dirs['down']
        sra_files = list(Path(down_dir).rglob('*.sra'))
        if not sra_files:
            UI.ok("no .sra files found, nothing to convert")
            self.ckpt.mark_done(s)
            return True

        UI.info(f"found {len(sra_files)} .sra file(s)")
        from sra2fastx import convert_all
        success, failures = convert_all(
            down_dir,
            jobs=self.args.convert_jobs,
            threads=self.args.convert_threads,
            progress=True,
            backend=self.args.convert_backend,
            verify=self.args.convert_verify,
        )
        if failures:
            UI.err(f"{len(failures)} conversion(s) failed")
            self.ckpt.mark_fail(s, f"{len(failures)} failures")
            return False
        UI.ok(f"{success} file(s) converted successfully")
        self.ckpt.mark_done(s)
        return True

    def run_plot(self) -> bool:
        s = 'plot'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)

        # 收集可用数据源: (标签, csv路径, 输出子目录)
        sources = []
        search_csv = os.path.join(self.dirs['search'], 'SRA_GSA_Merged_Final.csv')
        info_csv = os.path.join(self.dirs['info'], 'Global_Unified_Metadata_Core14.csv')
        if os.path.isfile(search_csv):
            sources.append(('search', search_csv, 'from_search'))
        if os.path.isfile(info_csv):
            sources.append(('info', info_csv, 'from_info'))
        if not sources:
            UI.err(f"missing both search and info output")
            self.ckpt.mark_fail(s)
            return False

        all_ok = True
        # 与画像阶段同理: 用 sys.executable, 裸 `python` 在 Windows 会解析到别的解释器
        py = shlex.quote(sys.executable or 'python')
        for label, csv_path, subdir in sources:
            outdir = os.path.join(self.dirs['plot'], subdir)
            os.makedirs(outdir, exist_ok=True)
            UI.info(f"绘制 [{label}] → {outdir}")
            cmd = (
                f'{py} {shlex.quote(self._bin("gsa_sra.plot.py"))} '
                f'--input {shlex.quote(csv_path)} '
                f'--outdir {shlex.quote(outdir)}'
            )
            rc = run_cmd(cmd, s, self.log_dir)
            if rc != 0:
                UI.err(f"plot [{label}] failed")
                all_ok = False

        if all_ok:
            UI.ok(f"output: {self.dirs['plot']}/from_search/ + from_info/")
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, "one or more plot sources failed")
        return False

    def run_hostref(self) -> bool:
        """宿主基因组下载: 直接调 download_host_genome.py (四通道自动回退); --host-fasta 时跳过."""
        s = 'hostref'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        if self.args.host_fasta:
            if not os.path.isfile(self.args.host_fasta):
                UI.err(f"host fasta not found: {self.args.host_fasta}")
                self.ckpt.mark_fail(s, "host fasta missing")
                return False
            UI.ok(f"user genome provided, skip download: {self.args.host_fasta}")
            self.ckpt.mark_done(s)
            return True
        host_species = self.args.host_species or self.args.species
        host_taxid = self.args.host_taxid or self.args.taxid
        genome_dir = os.path.join(self.work_dir, layout_dir_name('h_genome', self.layout))
        ncbi = f'--ncbi-api {shlex.quote(self.args.ncbi_api)}' if self.args.ncbi_api else ''
        cmd = (
            f'{self._py()} {shlex.quote(self._bin("download_host_genome.py"))} '
            f'--species {shlex.quote(str(host_species))} '
            f'--taxid {host_taxid} '
            f'--outdir {shlex.quote(genome_dir)} '
            f'--source {shlex.quote(self.args.host_download_source)} '
            f'{ncbi}'
        )
        rc = run_cmd(cmd, s, self.log_dir, secrets=self._secrets(), timeout=7200)
        fasta = os.path.join(genome_dir, 'all.genome.uniq.fasta')
        if rc == 0 and os.path.isfile(fasta):
            UI.ok(f"output: {fasta}")
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

    def run_hostdb(self) -> bool:
        """宿主索引构建: 直接调 build_hostbase.py (FASTA 来自 --host-fasta 或 hostref 产物)."""
        s = 'hostdb'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        host_taxid = self.args.host_taxid or self.args.taxid
        fasta = self.args.host_fasta or os.path.join(
            self.work_dir, layout_dir_name('h_genome', self.layout), 'all.genome.uniq.fasta')
        if not os.path.isfile(fasta):
            UI.err(f"missing genome fasta: {fasta} (先运行 hostref 阶段下载, 或 --host-fasta 指定)")
            self.ckpt.mark_fail(s, "missing genome FASTA")
            return False
        hostdb_dir = os.path.join(self.work_dir, layout_dir_name('h_hostdb', self.layout))
        cmd = (
            f'{self._py()} {shlex.quote(self._bin("build_hostbase.py"))} '
            f'--tool {shlex.quote(self.args.hostdb_tools)} '
            f'--input {shlex.quote(fasta)} '
            f'--output {shlex.quote(hostdb_dir)} '
            f'--threads {self.args.threads} '
            f'--seq-type {shlex.quote(self.HOST_SEQ_TYPES)} '
            f'--k2-libs {shlex.quote(self.HOST_K2_LIBS)} '
            f'--taxid {host_taxid}'
        )
        UI.warn("Kraken2 建库可能耗时数小时 (checkpoint 支持断点与跳过)")
        rc = run_cmd(cmd, s, self.log_dir, secrets=self._secrets(), timeout=86400 * 3)
        if rc == 0 and os.path.isdir(hostdb_dir):
            UI.ok(f"output: {hostdb_dir}/")
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

    def run_eve(self) -> bool:
        """[eve] 基因组 EVE 分析 (接在 hostdb 之后).

        掩码仍在 eve-mask 之外单独控制: 这里只做"发现 + 汇总", 不改已建好的宿主库。
        产物落 <host_reference>/eve/ (04_Summary/eve_report.html + 可选元件库/蛋白库),
        report 阶段会把 EVE 明细内嵌到 Pipeline_Summary_Report.html 的 S7；
        独立 eve_report.html 仅保留为可选导出。
        """
        s = 'eve'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        if not getattr(self.args, 'eve_scan', False):
            UI.warn("EVE 分析未启用 (--eve-scan): 跳过")
            self.ckpt.mark_done(s)
            return True

        genome_dir = os.path.join(self.work_dir,
                                  layout_dir_name('h_genome', self.layout))
        fasta = os.path.join(genome_dir, 'all.genome.uniq.fasta')
        if not os.path.isfile(fasta):
            UI.err(f"missing genome fasta: {fasta} (先跑 hostref, 或 --host-fasta)")
            self.ckpt.mark_fail(s, 'missing genome FASTA')
            return False
        eve_dir = os.path.join(genome_dir, 'eve')
        os.makedirs(eve_dir, exist_ok=True)
        species = (getattr(self.args, 'host_species', None) or self.args.species)
        name = re.sub(r'[^0-9A-Za-z_.-]+', '_', species).strip('_') or 'host'
        threads = str(getattr(self.args, 'threads', 30))

        screen = self._bin_multi('eve_genome', 'eve_screen.py')
        if not os.path.isfile(screen):
            UI.err(f"eve_screen.py 不存在: {screen}")
            self.ckpt.mark_fail(s, 'missing eve_screen.py')
            return False
        cmd = (f'{self._py()} {shlex.quote(screen)} -g {shlex.quote(fasta)} '
               f'-n {shlex.quote(name)} --outdir {shlex.quote(eve_dir)} '
               f'--stage 1,2,3 --fast -t {threads}')
        UI.warn("EVE 扫描耗时较长 (滑窗 blastx); 结果仅作分析/参考库")
        if run_cmd(cmd, s, self.log_dir, secrets=self._secrets(),
                   timeout=86400 * 2) != 0:
            self.ckpt.mark_fail(s, 'eve_screen failed')
            return False
        run_cmd(f'{self._py()} {shlex.quote(screen)} --merge -o {shlex.quote(eve_dir)}',
                s, self.log_dir, secrets=self._secrets(), timeout=86400)

        manifest = {'eve_root': eve_dir, 'scan_html': None,
                    'consensus_library': None, 'protein_library': None}
        if getattr(self.args, 'eve_elements', False):
            elem = self._bin_multi('eve_genome', 'eve_element_db.py')
            if os.path.isfile(elem):
                ecmd = (f'{self._py()} {shlex.quote(elem)} -o {shlex.quote(eve_dir)} '
                        f'--stages {shlex.quote(self.args.eve_element_stages)} '
                        f'-t {threads}')
                if getattr(self.args, 'eve_diamond', ''):
                    ecmd += f' --diamond {shlex.quote(self.args.eve_diamond)}'
                if getattr(self.args, 'eve_ictv_prot', ''):
                    ecmd += f' --ictv-prot {shlex.quote(self.args.eve_ictv_prot)}'
                run_cmd(ecmd, s, self.log_dir, secrets=self._secrets(),
                        timeout=86400 * 3)
                lib = os.path.join(eve_dir, '05_Elements',
                                   'eve_consensus_library.fa')
                if os.path.isfile(lib):
                    manifest['consensus_library'] = lib
                    prot = self._bin_multi('eve_genome', 'eve_prot_lib.py')
                    faa = lib.replace('.fa', '.sixframe.faa')
                    if os.path.isfile(prot) and run_cmd(
                            f'{self._py()} {shlex.quote(prot)} --in {shlex.quote(lib)} '
                            f'--out {shlex.quote(faa)}', s, self.log_dir,
                            secrets=self._secrets(), timeout=86400) == 0:
                        manifest['protein_library'] = faa

        rep = self._bin_multi('eve_genome', 'eve_scan_report.py')
        if os.path.isfile(rep):
            html = os.path.join(eve_dir, '04_Summary', 'eve_report.html')
            run_cmd(f'{self._py()} {shlex.quote(rep)} -o {shlex.quote(eve_dir)} '
                    f'--out {shlex.quote(html)}', s, self.log_dir,
                    secrets=self._secrets(), timeout=3600)
            if os.path.isfile(html):
                manifest['scan_html'] = html
        with open(os.path.join(eve_dir, 'eve_manifest.json'), 'w',
                  encoding='utf-8') as fh:
            json.dump(manifest, fh, ensure_ascii=False, indent=2)
        UI.ok(f"EVE 报告: {manifest['scan_html'] or eve_dir}")
        if manifest['protein_library']:
            UI.ok(f"EVE 蛋白库 (供发现管线 --eve-faa): {manifest['protein_library']}")
        self.ckpt.mark_done(s)
        return True

    def _bin_multi(self, subdir: str, script: str) -> str:
        return os.path.join(self.bin_dir, subdir, script)

    def run_report(self) -> bool:
        s = 'report'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        # sys.executable: 与 plot/画像阶段同理, 不用裸 `python` (PATH 可能指向别的解释器)
        py = shlex.quote(sys.executable or 'python')
        # 配置了 DeepSeek 时把 key/model 透传给报告生成器, 报告里的 AI 摘要才有内容;
        # 未配置则不传 (generate_report.py 的 --ai-key 默认空, 不会误触发)
        ai_arg = ''
        if self.args.deepseek_api:
            ai_arg = (f' --ai-key {shlex.quote(self.args.deepseek_api)}'
                      f' --ai-provider deepseek'
                      f' --ai-model {shlex.quote(self.args.deepseek_model)}')
        cmd = f'{py} {shlex.quote(self._bin("generate_report.py"))} -d {shlex.quote(self.work_dir)}{ai_arg}'
        rc = run_cmd(cmd, s, self.log_dir, secrets=self._secrets())
        html = os.path.join(self.work_dir, 'Pipeline_Summary_Report.html')
        handoff = os.path.join(self.work_dir, 'sample_handoff.csv')
        if rc == 0 and os.path.isfile(html):
            UI.ok(f"report:  {html}")
            if os.path.isfile(handoff):
                UI.ok(f"handoff: {handoff}")
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Public Data Pipeline v3.2 - 公共数据获取管道",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage search info plot

  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage all

  # 预筛 + 过滤下载 (safe: 不漏「无分类数据」的未知样本)
  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --stage prescreen down --prescreen-filter safe \\
      --prescreen-virus-pct-ident 0.01
        """
    )
    parser.add_argument('--species', required=True, help='物种拉丁学名')
    parser.add_argument('--taxid', type=int, required=True, help='NCBI Taxonomy ID')
    parser.add_argument('--stage', nargs='+', required=True,
                        help=f'执行阶段 (空格或逗号分隔, 如: --stage search info 或 '
                             f'--stage search,info). 可选: {PublicDataPipeline.STAGES + ["all"]}')
    parser.add_argument('--deepseek-api', help='DeepSeek API Key')
    parser.add_argument('--deepseek-model', default='deepseek-v4-flash', help='DeepSeek 模型')
    parser.add_argument('--deepseek-search', action='store_true',
                        help='search 阶段启用 DeepSeek 网页检索补全缺失 metadata '
                             '(采样时间/采样地/组织/生长年限; ¥5/千次检索)')
    parser.add_argument('--deepseek-fetch', action='store_true',
                        help='检索证据不足 (rerank<1.5) 时追加 Fetch 抓取网页正文')
    parser.add_argument('--ncbi-api', help='NCBI API Key')
    parser.add_argument('--work-dir', default='./public_data_pipeline_output', help='输出根目录')
    parser.add_argument('--io-layout', choices=['legacy', 'standard'], default=None,
                        help='I/O 目录布局: standard=统一编号 (默认); legacy=v3.0 现行名 '
                             '(存量项目/checkpoint 兼容; doc/IO_LAYOUT_DESIGN.md); '
                             '亦可 MMPV_IO_LAYOUT 环境变量')
    parser.add_argument('--bin-dir', help='脚本目录 (默认自动检测)')
    parser.add_argument('--source-type', default='TRANSCRIPTOMIC', help='测序类型')
    parser.add_argument('--detailed', action='store_true', default=True, help='详细模式')
    parser.add_argument('--no-detailed', action='store_false', dest='detailed', help='关闭详细模式')
    # ---- prescreen 阶段 (病毒预筛 & 物种鉴定) ----
    parser.add_argument('--prescreen-filter', default='off', choices=['off', 'safe', 'strict'],
                        help='down 是否按预筛结论过滤下载清单: off=不过滤(默认, 行为不变); '
                             'safe=保留 DOWNLOAD+REVIEW+LOCAL_CHECK(不漏无分类数据的样本); '
                             'strict=只保留 DOWNLOAD+REVIEW(最省资源)')
    parser.add_argument('--prescreen-jobs', type=int, default=6, help='prescreen NCBI 并发数')
    parser.add_argument('--prescreen-gsa-jobs', type=int, default=3,
                        help='prescreen GSA 并发数 (NGDC 有反爬, 勿调高)')
    parser.add_argument('--prescreen-virus-pct-ident', type=float, default=None,
                        help='prescreen 归一化病毒占比阈值(%%) —— 病毒命中/已鉴定命中, '
                             '跨 NCBI/GSA 可比, 推荐用这个; 同时驱动 decision 列')
    parser.add_argument('--prescreen-virus-pct', type=float, default=None,
                        help='prescreen 旧口径病毒占比阈值(%%)')
    parser.add_argument('--prescreen-no-gsa', action='store_true',
                        help='prescreen 只跑 NCBI STAT, 跳过 GSA 通路')
    # ---- prescreen 之后的病毒层面画像 (名录/丰度/检出率 + 图) ----
    parser.add_argument('--virus-profile', dest='virus_profile', action='store_true',
                        default=True,
                        help='prescreen 后生成病毒层面画像 (默认开; 纯附加产物, '
                             '不改 sample_decision.tsv, 也不影响 down)')
    parser.add_argument('--no-virus-profile', dest='virus_profile', action='store_false',
                        help='关闭病毒层面画像')
    parser.add_argument('--virus-profile-plot', dest='virus_profile_plot',
                        action='store_true', default=True,
                        help='画像生成图表 (默认开; 缺 matplotlib/seaborn 时自动降级为只出表)')
    parser.add_argument('--no-virus-profile-plot', dest='virus_profile_plot',
                        action='store_false', help='画像只出表, 不绘图')
    parser.add_argument('--virus-profile-min-samples', type=int, default=2,
                        help='画像收录的样本数下限 (默认 2, 滤掉单例噪声; 表仍保留全部)')
    parser.add_argument('--virus-profile-top-n', type=int, default=25,
                        help='画像条形图展示的病毒数 (默认 25)')
    parser.add_argument('--virus-profile-by', default='leaf', choices=['leaf', 'all'],
                        help='画像绘图口径: leaf=末端节点(默认, 父子不并列); all=任意层级')
    parser.add_argument('--taxdump', help='NCBI taxdump 目录(含 nodes.dmp), 用于精确病毒子树判定')
    parser.add_argument('--ngdc-method', default='aria2c', choices=['aria2c', 'requests', 'wget'])
    parser.add_argument('--ngdc-concurrency', type=int, default=5)
    parser.add_argument('--prefetch-concurrency', type=int, default=3)
    parser.add_argument('--skip-list', help='已下载样本跳过列表 (传给 gsa_sra.down.py)')
    parser.add_argument('--threads', type=int, default=4, help='CPU 线程数')
    parser.add_argument('--convert-jobs', type=int, default=2, help='SRA 转换并行任务数')
    parser.add_argument('--convert-threads', type=int, default=6, help='每个转换任务的线程数')
    parser.add_argument('--convert-backend', default='auto', choices=['auto', 'sracha', 'sratools'],
                        help='SRA→FASTQ 后端: auto=sracha 优先、失败回退 sra-tools (默认)')
    parser.add_argument('--convert-verify', action='store_true',
                        help='sracha: 解码后额外核对 quality 数值分布 (更严格, 明显变慢)')
    parser.add_argument('--host-species', help='宿主物种拉丁学名 (默认复用 --species)')
    parser.add_argument('--host-taxid', type=int, help='宿主 NCBI TaxID (默认复用 --taxid)')
    parser.add_argument('--host-fasta', help='已有宿主基因组 FASTA (hostref 阶段跳过下载)')
    parser.add_argument('--host-download-source', default='auto',
                        choices=['auto', 'datasets', 'ngd', 'ftp', 'gget'],
                        help='宿主基因组下载通道 (默认: auto 自动回退; gget 走 Ensembl)')
    parser.add_argument('--eve-scan', action='store_true',
                        help='hostdb 后做基因组 EVE 分析 (默认关; 产物并入 HTML 报告)')
    parser.add_argument('--eve-elements', action='store_true',
                        help='EVE 扫描后继续跑元件层 (聚类/consensus/CDD) 并产蛋白库')
    parser.add_argument('--eve-element-stages', default='collect,cluster,consensus,cdd',
                        help='EVE 元件层阶段 (default: collect,cluster,consensus,cdd)')
    parser.add_argument('--eve-diamond', default='', help='EVE 元件层 DIAMOND 可执行文件')
    parser.add_argument('--eve-ictv-prot', default='', help='ICTV 六框蛋白库 (EVE 定科用)')
    parser.add_argument('--hostdb-tools', default='kraken2,bowtie2,hisat2,minimap2',
                        help='hostref 建库工具 (默认: kraken2,bowtie2,hisat2,minimap2)')
    parser.add_argument('--force', action='store_true', help='强制重新执行 (忽略已完成状态)')
    parser.add_argument('--dry-run', action='store_true', help='仅预览')

    args = parser.parse_args()
    # [multi-stage] 空格/逗号统一拆分 + 校验 (合法值与旧 choices 一致)
    args.stage = parse_multi_stage(args.stage, tuple(PublicDataPipeline.STAGES) + ('all',))
    if not args.deepseek_api:
        args.deepseek_api = os.environ.get('DEEPSEEK_API_KEY') or os.environ.get('MMPV_AI_API_KEY', '')
    if not args.ncbi_api:
        args.ncbi_api = os.environ.get('NCBI_API_KEY', '')
    if not args.bin_dir:
        args.bin_dir = os.path.dirname(os.path.abspath(__file__))

    UI.banner("Public Data Pipeline v3.2")

    if 'all' in args.stage:
        stages_to_run = list(PublicDataPipeline.STAGES)
    else:
        stages_to_run = [s for s in args.stage if s in PublicDataPipeline.STAGES]

    ckpt = Checkpoint(os.path.abspath(args.work_dir))
    if args.force:
        ckpt.reset()

    print(f"  species:      {args.species} (taxid: {args.taxid})")
    print(f"  source type:  {args.source_type}")
    print(f"  work dir:     {os.path.abspath(args.work_dir)}")
    print(f"  threads:      {args.threads}")
    print(f"  convert:      jobs={args.convert_jobs}, threads={args.convert_threads}")
    print(f"  deepseek:     {'configured' if args.deepseek_api else 'NOT set'}")
    print(f"  ncbi api:     {'configured' if args.ncbi_api else 'NOT set'}")
    print(f"  stage plan:   {' -> '.join(stages_to_run)}")
    print(f"\n{ckpt.summary(PublicDataPipeline.STAGES)}")

    if args.dry_run:
        print(f"\n{UI.C['cyan']}dry-run mode.{UI.C['reset']}")
        for s in stages_to_run:
            print(f"  [{'SKIP' if ckpt.is_done(s) else 'RUN'}] {s}: {PublicDataPipeline.STAGE_DESC[s]}")
        return

    p = PublicDataPipeline(args)
    funcs = {'search': p.run_search, 'info': p.run_info, 'prescreen': p.run_prescreen,
             'down': p.run_down, 'convert': p.run_convert, 'plot': p.run_plot,
             'hostref': p.run_hostref, 'hostdb': p.run_hostdb,
             'eve': p.run_eve, 'report': p.run_report}
    failed = []
    t0 = time.time()
    for s in stages_to_run:
        if ckpt.is_done(s) and not args.force:
            UI.stage(PublicDataPipeline.STAGE_DESC[s], 'skip')
            continue
        try:
            if not funcs[s]():
                failed.append(s)
                UI.err(f"stage [{s}] failed")
                break
        except KeyboardInterrupt:
            UI.warn("interrupted")
            failed.append(s)
            break
        except Exception as e:
            UI.err(f"stage [{s}] exception: {e}")
            failed.append(s)
            break

    elapsed = time.time() - t0
    print(f"\n{UI.C['purple']}{UI.C['bold']}{'=' * 60}{UI.C['reset']}")
    if failed:
        print(f"  {UI.C['red']}[FAILED] stages: {', '.join(failed)}{UI.C['reset']}")
        print(f"  re-run: --stage {' '.join(failed)}")
        sys.exit(1)
    else:
        print(f"  {UI.C['green']}[SUCCESS] elapsed: {elapsed / 60:.1f} min{UI.C['reset']}")
        for name, path in p.dirs.items():
            if os.path.isdir(path):
                print(f"    {name}: {path}")


if __name__ == '__main__':
    main()
