#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Public Data Pipeline - 公共数据获取管道 v3.1
=============================================
八阶段: search -> info -> down -> convert -> plot -> hostref -> hostdb -> report
全部参数命令行指定。宿主参考基因组默认复用 --species/--taxid (检索物种即宿主物种)。
stage 直接调用工作脚本 (download_host_genome.py / build_hostbase.py), 不嵌套子管道。

用法:
  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage search info plot

  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage all
"""

import os
import sys
import time
import csv
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

from utils.pipeline_utils import UI, Checkpoint, run_cmd

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


# ==========================================
# Public Data Pipeline
# ==========================================
class PublicDataPipeline:
    STAGES = ['search', 'info', 'down', 'convert', 'plot', 'hostref', 'hostdb', 'report']

    STAGE_DESC = {
        'search':  'GSA + SRA 双引擎物种检索',
        'info':    '元数据深度解析 & 文献溯源',
        'down':    '高通量数据下载 (NGDC + NCBI)',
        'convert': 'SRA → FASTQ.GZ 转换',
        'plot':    'SCI 级六图可视化',
        'hostref': '宿主基因组下载 (datasets/ngd/FTP/gget 四通道; --host-fasta 可跳过)',
        'hostdb':  '宿主索引构建 (Kraken2/Bowtie2/HISAT2/Minimap2)',
        'report':  'HTML 报告 + sample_handoff.csv 交接清单',
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
        # 布局解析 (CLI > MMPV_IO_LAYOUT > legacy) 并写回环境供子进程继承
        self.layout = normalize_layout_env(getattr(args, 'io_layout', None))
        self.dirs = build_meta_dirs(Path(self.work_dir), self.layout)
        for _name, _path in self.dirs.items():
            os.makedirs(_path, exist_ok=True)

    def _bin(self, script: str) -> str:
        p = os.path.join(self.bin_dir, script)
        if not os.path.isfile(p) and not script.endswith('.py'):
            p = os.path.join(self.bin_dir, script + '.py')
        return p

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
        cmd = (
            f'python {shlex.quote(self._bin("gsa_sra.search.py"))} '
            f'--query {shlex.quote(self.args.species)} '
            f'--source {shlex.quote(self.args.source_type)} '
            f'{detailed} {self._api_cli()} '
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
        cmd = (
            f'python {shlex.quote(self._bin("gsa_sra.info.py"))} '
            f'--input {shlex.quote(sra_list)} --mode {shlex.quote(mode)} --fill-date '
            f'{self._api_cli()} --threads {self.args.threads} '
            f'--outdir {shlex.quote(self.dirs["info"])}'
        )
        rc = run_cmd(cmd, s, self.log_dir, secrets=self._secrets())
        if rc == 0:
            UI.ok(f"output: {self.dirs['info']}/Global_Unified_Metadata_Core14.csv")
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

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
        skip_arg = f'--skip-list {shlex.quote(self.args.skip_list)}' if self.args.skip_list else ''
        links_file = os.path.join(self.dirs['search'], 'download_links.csv')
        links_arg = f'--links {shlex.quote(links_file)}' if os.path.isfile(links_file) else ''
        cmd = (
            f'python {shlex.quote(self._bin("gsa_sra.down.py"))} '
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
        for label, csv_path, subdir in sources:
            outdir = os.path.join(self.dirs['plot'], subdir)
            os.makedirs(outdir, exist_ok=True)
            UI.info(f"绘制 [{label}] → {outdir}")
            cmd = (
                f'python {shlex.quote(self._bin("gsa_sra.plot.py"))} '
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
            f'python {shlex.quote(self._bin("download_host_genome.py"))} '
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
            f'python {shlex.quote(self._bin("build_hostbase.py"))} '
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

    def run_report(self) -> bool:
        s = 'report'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        cmd = f'python {shlex.quote(self._bin("generate_report.py"))} -d {shlex.quote(self.work_dir)}'
        rc = run_cmd(cmd, s, self.log_dir)
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
        description="Public Data Pipeline v3.1 - 公共数据获取管道",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage search info plot

  python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage all
        """
    )
    parser.add_argument('--species', required=True, help='物种拉丁学名')
    parser.add_argument('--taxid', type=int, required=True, help='NCBI Taxonomy ID')
    parser.add_argument('--stage', nargs='+', required=True,
                        choices=PublicDataPipeline.STAGES + ['all'],
                        help=f'执行阶段. 可选: {PublicDataPipeline.STAGES}')
    parser.add_argument('--deepseek-api', help='DeepSeek API Key')
    parser.add_argument('--deepseek-model', default='deepseek-v4-flash', help='DeepSeek 模型')
    parser.add_argument('--ncbi-api', help='NCBI API Key')
    parser.add_argument('--work-dir', default='./public_data_pipeline_output', help='输出根目录')
    parser.add_argument('--io-layout', choices=['legacy', 'standard'], default=None,
                        help='I/O 目录布局: legacy=v3.0 现行名 (默认); standard=统一编号'
                             '(doc/IO_LAYOUT_DESIGN.md); 亦可 MMPV_IO_LAYOUT 环境变量')
    parser.add_argument('--bin-dir', help='脚本目录 (默认自动检测)')
    parser.add_argument('--source-type', default='TRANSCRIPTOMIC', help='测序类型')
    parser.add_argument('--detailed', action='store_true', default=True, help='详细模式')
    parser.add_argument('--no-detailed', action='store_false', dest='detailed', help='关闭详细模式')
    parser.add_argument('--ngdc-method', default='aria2c', choices=['aria2c', 'requests', 'wget'])
    parser.add_argument('--ngdc-concurrency', type=int, default=5)
    parser.add_argument('--prefetch-concurrency', type=int, default=3)
    parser.add_argument('--skip-list', help='已下载样本跳过列表 (传给 gsa_sra.down.py)')
    parser.add_argument('--threads', type=int, default=4, help='CPU 线程数')
    parser.add_argument('--convert-jobs', type=int, default=2, help='SRA 转换并行任务数')
    parser.add_argument('--convert-threads', type=int, default=6, help='每个 fasterq-dump 的线程数')
    parser.add_argument('--host-species', help='宿主物种拉丁学名 (默认复用 --species)')
    parser.add_argument('--host-taxid', type=int, help='宿主 NCBI TaxID (默认复用 --taxid)')
    parser.add_argument('--host-fasta', help='已有宿主基因组 FASTA (hostref 阶段跳过下载)')
    parser.add_argument('--host-download-source', default='auto',
                        choices=['auto', 'datasets', 'ngd', 'ftp', 'gget'],
                        help='宿主基因组下载通道 (默认: auto 自动回退; gget 走 Ensembl)')
    parser.add_argument('--hostdb-tools', default='kraken2,bowtie2,hisat2,minimap2',
                        help='hostref 建库工具 (默认: kraken2,bowtie2,hisat2,minimap2)')
    parser.add_argument('--force', action='store_true', help='强制重新执行 (忽略已完成状态)')
    parser.add_argument('--dry-run', action='store_true', help='仅预览')

    args = parser.parse_args()
    if not args.deepseek_api:
        args.deepseek_api = os.environ.get('DEEPSEEK_API_KEY', '')
    if not args.ncbi_api:
        args.ncbi_api = os.environ.get('NCBI_API_KEY', '')
    if not args.bin_dir:
        args.bin_dir = os.path.dirname(os.path.abspath(__file__))

    UI.banner("Public Data Pipeline v3.1")

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
    funcs = {'search': p.run_search, 'info': p.run_info, 'down': p.run_down,
             'convert': p.run_convert, 'plot': p.run_plot, 'hostref': p.run_hostref,
             'hostdb': p.run_hostdb, 'report': p.run_report}
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
