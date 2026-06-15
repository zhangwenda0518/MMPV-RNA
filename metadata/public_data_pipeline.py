#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Public Data Pipeline - 公共数据获取管道 v3.1
=============================================
四阶段: search -> info -> down -> plot
全部参数命令行指定。

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

from pipeline_utils import UI, Checkpoint, run_cmd

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


# ==========================================
# Public Data Pipeline
# ==========================================
class PublicDataPipeline:
    STAGES = ['search', 'info', 'down', 'plot']

    STAGE_DESC = {
        'search': 'GSA + SRA 双引擎物种检索',
        'info':   '元数据深度解析 & 文献溯源',
        'down':   '高通量数据下载 (NGDC + NCBI)',
        'plot':   'SCI 级六图可视化',
    }

    def __init__(self, args):
        self.args = args
        self.bin_dir = args.bin_dir or os.path.dirname(os.path.abspath(__file__))
        self.work_dir = os.path.abspath(args.work_dir)
        self.ckpt = Checkpoint(self.work_dir)
        self.log_dir = os.path.join(self.work_dir, 'logs')
        self.dirs = {}
        for name in ['search', 'plot', 'info', 'down']:
            self.dirs[name] = os.path.join(self.work_dir, name)
            os.makedirs(self.dirs[name], exist_ok=True)

    def _bin(self, script):
        p = os.path.join(self.bin_dir, script)
        if not os.path.isfile(p) and not script.endswith('.py'):
            p = os.path.join(self.bin_dir, script + '.py')
        return p

    def _secrets(self):
        """Return list of API key strings to mask in logs."""
        s = []
        if self.args.deepseek_api:
            s.append(self.args.deepseek_api)
        if self.args.ncbi_api:
            s.append(self.args.ncbi_api)
        return s

    def _api_cli(self):
        parts = []
        if self.args.deepseek_api:
            parts.append(f'--deepseek-api {shlex.quote(self.args.deepseek_api)}')
            parts.append(f'--deepseek-model {shlex.quote(self.args.deepseek_model)}')
        if self.args.ncbi_api:
            parts.append(f'--ncbi-api {shlex.quote(self.args.ncbi_api)}')
        return ' '.join(parts)

    def _extract_runs(self):
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

    def run_search(self):
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

    def run_info(self):
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
        cmd = (
            f'python {shlex.quote(self._bin("gsa_sra.info.py"))} '
            f'--input {shlex.quote(sra_list)} --mode {shlex.quote(mode)} --fill-date '
            f'{self._api_cli()} --threads {self.args.threads} '
            f'--outdir {shlex.quote(self.dirs["info"])}'
        )
        rc = run_cmd(cmd, s, self.log_dir, secrets=self._secrets())
        if rc == 0:
            UI.ok(f"output: {self.dirs['info']}/Global_Unified_Metadata_Core13.csv")
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

    def run_down(self):
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
        cmd = (
            f'python {shlex.quote(self._bin("gsa_sra.down.py"))} '
            f'--list {shlex.quote(sra_list)} '
            f'--ngdc-method {shlex.quote(self.args.ngdc_method)} '
            f'--ngdc-concurrency {self.args.ngdc_concurrency} '
            f'--prefetch-concurrency {self.args.prefetch_concurrency} '
            f'{skip_arg} '
            f'--output {shlex.quote(self.dirs["down"])}'
        )
        UI.warn("data download may take hours to days")
        rc = run_cmd(cmd, s, self.log_dir, timeout=86400 * 7)
        if rc == 0:
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False

    def run_plot(self):
        s = 'plot'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)
        merged = os.path.join(self.dirs['search'], 'SRA_GSA_Merged_Final.csv')
        if not os.path.isfile(merged):
            UI.err(f"missing search output: {merged}")
            self.ckpt.mark_fail(s)
            return False
        cmd = (
            f'python {shlex.quote(self._bin("gsa_sra.plot.py"))} '
            f'--input {shlex.quote(merged)} '
            f'--outdir {shlex.quote(self.dirs["plot"])}'
        )
        rc = run_cmd(cmd, s, self.log_dir)
        if rc == 0:
            UI.ok(f"output: {self.dirs['plot']}/Combined_Landscape_Full.pdf")
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
    parser.add_argument('--bin-dir', help='脚本目录 (默认自动检测)')
    parser.add_argument('--source-type', default='TRANSCRIPTOMIC', help='测序类型')
    parser.add_argument('--detailed', action='store_true', default=True, help='详细模式')
    parser.add_argument('--no-detailed', action='store_false', dest='detailed', help='关闭详细模式')
    parser.add_argument('--ngdc-method', default='aria2c', choices=['aria2c', 'requests', 'wget'])
    parser.add_argument('--ngdc-concurrency', type=int, default=5)
    parser.add_argument('--prefetch-concurrency', type=int, default=3)
    parser.add_argument('--skip-list', help='已下载样本跳过列表 (传给 gsa_sra.down.py)')
    parser.add_argument('--threads', type=int, default=4, help='CPU 线程数')
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
    funcs = {'search': p.run_search, 'info': p.run_info, 'down': p.run_down, 'plot': p.run_plot}
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
