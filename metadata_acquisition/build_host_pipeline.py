#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build Host Pipeline - 宿主数据库构建管道 v3.1
===============================================
两阶段: genome-down -> hostdb
全部参数命令行指定。

用法:
  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --ncbi-api "xxx" --stage genome-down hostdb

  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --ncbi-api "xxx" --taxonomy-dir ~/database/taxonomy --stage all

  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --genome-fasta /path/to/existing.fasta --stage hostdb
"""

import os
import sys
import time
import json
import gzip
import shutil
import argparse
import subprocess
import shlex

from pipeline_utils import UI, Checkpoint, run_cmd

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


# ==========================================
# Build Host Pipeline
# ==========================================
class BuildHostPipeline:
    STAGES = ['genome-down', 'hostdb']

    STAGE_DESC = {
        'genome-down': '宿主参考基因组下载 (NCBI datasets)',
        'hostdb':      '宿主竞争性数据库构建 (K2/BT2/HT2/MM2)',
    }

    def __init__(self, args):
        self.args = args
        self.bin_dir = args.bin_dir or os.path.dirname(os.path.abspath(__file__))
        self.work_dir = os.path.abspath(args.work_dir)
        self.ckpt = Checkpoint(self.work_dir)
        self.log_dir = os.path.join(self.work_dir, 'logs')
        self.dirs = {
            'genome': os.path.join(self.work_dir, 'genome'),
            'hostdb': os.path.join(self.work_dir, 'hostdb'),
        }
        for d in self.dirs.values():
            os.makedirs(d, exist_ok=True)

    def _bin(self, script):
        p = os.path.join(self.bin_dir, script)
        if not os.path.isfile(p) and not script.endswith('.py'):
            p = os.path.join(self.bin_dir, script + '.py')
        return p

    def _secrets(self):
        """Return list of API key strings to mask in logs."""
        s = []
        if self.args.ncbi_api:
            s.append(self.args.ncbi_api)
        return s

    # ── genome-down ────────────────────────────────
    def run_genome_down(self):
        s = 'genome-down'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)

        genome_dir = self.dirs['genome']
        genome_zip = os.path.join(genome_dir, 'genome_down.zip')

        # 如果用户直接提供了 FASTA，跳过下载
        if self.args.genome_fasta and os.path.isfile(self.args.genome_fasta):
            UI.ok(f"using provided genome: {self.args.genome_fasta}")
            self._genome_fasta = self.args.genome_fasta
            self.ckpt.mark_done(s)
            return True

        if not os.path.isfile(genome_zip):
            if shutil.which('datasets'):
                cmd = (
                    f'datasets download genome taxon {shlex.quote(self.args.species)} '
                    f'--filename {shlex.quote(genome_zip)} --include genome,gff3,seq-report'
                )
                if self.args.ncbi_api:
                    cmd += f' --api-key {shlex.quote(self.args.ncbi_api)}'
                rc = run_cmd(cmd, s, self.log_dir, timeout=3600, secrets=self._secrets())
                if rc != 0:
                    self.ckpt.mark_fail(s, "datasets CLI failed")
                    return False
                UI.ok("datasets download complete")
            else:
                UI.err("datasets CLI not found. install: conda install -c conda-forge ncbi-datasets-cli")
                self.ckpt.mark_fail(s, "datasets CLI missing")
                return False
        else:
            UI.ok("genome zip exists, skip download")

        # unzip — use list form (no shell=True) to avoid injection
        extracted = os.path.join(genome_dir, 'extracted')
        if not os.path.isdir(extracted) or not os.listdir(extracted):
            os.makedirs(extracted, exist_ok=True)
            UI.info("extracting genome...")
            try:
                subprocess.run(['unzip', '-o', genome_zip, '-d', extracted],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                import zipfile
                try:
                    with zipfile.ZipFile(genome_zip, 'r') as zf:
                        zf.extractall(extracted)
                except Exception as e:
                    self.ckpt.mark_fail(s, f"unzip failed: {e}")
                    return False
            UI.ok("extraction complete")

        # collect and merge FASTA
        fasta_exts = {'.fna', '.fasta', '.fa', '.fna.gz', '.fasta.gz', '.fa.gz'}
        fasta_files = []
        for root, _, files in os.walk(extracted):
            for f in files:
                ext = os.path.splitext(f)[1]
                if ext in fasta_exts or (f.endswith('.gz') and os.path.splitext(f[:-3])[1] in fasta_exts):
                    fasta_files.append(os.path.join(root, f))
        if not fasta_files:
            self.ckpt.mark_fail(s, "no FASTA files found")
            return False

        merged_fasta = os.path.join(genome_dir, 'all.genome.uniq.fasta')
        UI.info(f"merging {len(fasta_files)} fasta files -> {merged_fasta}")

        seen_ids, dup_count = set(), 0
        with open(merged_fasta, 'w', encoding='utf-8') as out:
            out.write(f"# Merged host genome: {self.args.species}\n")
            out.write(f"# Created: {datetime.now().isoformat()}\n#\n")
            for fa in fasta_files:
                try:
                    if fa.endswith('.gz'):
                        with gzip.open(fa, 'rt', encoding='utf-8') as inf:
                            content = inf.read()
                    else:
                        with open(fa, 'r', encoding='utf-8') as inf:
                            content = inf.read()
                except Exception:
                    continue
                for line in content.split('\n'):
                    if line.startswith('>'):
                        seq_id = line[1:].strip().split()[0]
                        if seq_id in seen_ids:
                            suffix = 1
                            while f"{seq_id}_dup{suffix}" in seen_ids:
                                suffix += 1
                            new_id = f"{seq_id}_dup{suffix}"
                            seen_ids.add(new_id)
                            rest = line[1:].strip()[len(seq_id):]
                            out.write(f">{new_id} {rest}\n".rstrip() + '\n')
                            dup_count += 1
                        else:
                            seen_ids.add(seq_id)
                            out.write(line + '\n')
                    else:
                        out.write(line + '\n')

        file_mb = os.path.getsize(merged_fasta) / (1024 * 1024)
        UI.ok(f"merge complete: {len(seen_ids)} sequences, {file_mb:.1f} MB")
        if dup_count:
            UI.warn(f"renamed {dup_count} duplicates")

        self._genome_fasta = merged_fasta
        self.ckpt.mark_done(s)
        return True

    # ── hostdb ─────────────────────────────────────
    def run_hostdb(self):
        s = 'hostdb'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)

        # 查找 FASTA
        genome_fasta = getattr(self, '_genome_fasta', None)
        if not genome_fasta or not os.path.isfile(genome_fasta):
            genome_fasta = self.args.genome_fasta or os.path.join(self.dirs['genome'], 'all.genome.uniq.fasta')
        if not os.path.isfile(genome_fasta):
            UI.err(f"missing genome fasta: {genome_fasta}")
            self.ckpt.mark_fail(s, "missing genome FASTA")
            return False

        cmd_parts = [
            f'python {shlex.quote(self._bin("build_hostbase.py"))}',
            f'--tool {shlex.quote(self.args.hostdb_tools)}',
            f'--input {shlex.quote(genome_fasta)}',
            f'--output {shlex.quote(self.dirs["hostdb"])}',
            f'--threads {self.args.threads}',
            f'--seq-type {shlex.quote(self.args.seq_types)}',
            f'--k2-libs {shlex.quote(self.args.k2_libs)}',
            f'--taxid {self.args.taxid}',
        ]
        if self.args.taxonomy_dir:
            if os.path.isfile(os.path.join(self.args.taxonomy_dir, 'names.dmp')):
                cmd_parts.append(f'--taxonomy {shlex.quote(self.args.taxonomy_dir)}')
                UI.ok(f"taxonomy validated: {self.args.taxonomy_dir}")
            else:
                UI.err(f"taxonomy dir missing names.dmp: {self.args.taxonomy_dir}")
                self.ckpt.mark_fail(s, "invalid taxonomy dir")
                return False

        UI.warn("Kraken2 build may take several hours")
        rc = run_cmd(' '.join(cmd_parts), s, self.log_dir, timeout=86400 * 3)
        if rc == 0:
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Build Host Pipeline v3.1 - 宿主数据库构建管道",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --ncbi-api "xxx" --stage genome-down hostdb

  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --ncbi-api "xxx" --taxonomy-dir ~/database/taxonomy --stage all

  # 使用已有基因组 FASTA 直接建库 (跳过 genome-down)
  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --genome-fasta /path/to/genome.fasta --stage hostdb
        """
    )
    parser.add_argument('--species', required=True, help='物种拉丁学名')
    parser.add_argument('--taxid', type=int, required=True, help='NCBI Taxonomy ID')
    parser.add_argument('--stage', nargs='+', required=True,
                        choices=BuildHostPipeline.STAGES + ['all'],
                        help=f'执行阶段. 可选: {BuildHostPipeline.STAGES}')
    parser.add_argument('--ncbi-api', help='NCBI API Key')
    parser.add_argument('--work-dir', default='./build_host_pipeline_output', help='输出根目录')
    parser.add_argument('--bin-dir', help='脚本目录 (默认自动检测)')
    parser.add_argument('--genome-fasta', help='已有基因组 FASTA 路径 (跳过 genome-down)')
    parser.add_argument('--taxonomy-dir', help='Kraken2 本地 Taxonomy 目录')
    parser.add_argument('--hostdb-tools', default='kraken2,bowtie2,hisat2,minimap2',
                        help='建库工具 (default: kraken2,bowtie2,hisat2,minimap2)')
    parser.add_argument('--seq-types', default='dna-short,rna-short,nanopore,pacbio',
                        help='Minimap2 序列类型 (default: dna-short,rna-short,nanopore,pacbio)')
    parser.add_argument('--k2-libs', default='archaea,bacteria,plasmid,fungi,protozoa,UniVec',
                        help='Kraken2 标准库. none=跳过 (default: archaea,bacteria,...)')
    parser.add_argument('--threads', type=int, default=30, help='CPU 线程数 (default: 30)')
    parser.add_argument('--force', action='store_true', help='强制重新执行 (忽略已完成状态)')
    parser.add_argument('--dry-run', action='store_true', help='仅预览')

    args = parser.parse_args()
    if not args.ncbi_api:
        args.ncbi_api = os.environ.get('NCBI_API_KEY', '')
    if not args.bin_dir:
        args.bin_dir = os.path.dirname(os.path.abspath(__file__))

    UI.banner("Build Host Pipeline v3.1")

    if 'all' in args.stage:
        stages_to_run = list(BuildHostPipeline.STAGES)
    else:
        stages_to_run = [s for s in BuildHostPipeline.STAGES if s in args.stage]

    ckpt = Checkpoint(os.path.abspath(args.work_dir))
    if args.force:
        ckpt.reset()

    print(f"  species:      {args.species} (taxid: {args.taxid})")
    print(f"  work dir:     {os.path.abspath(args.work_dir)}")
    print(f"  ncbi api:     {'configured' if args.ncbi_api else 'NOT set'}")
    print(f"  genome fasta: {args.genome_fasta or '(from datasets CLI)'}")
    print(f"  hostdb tools: {args.hostdb_tools}")
    print(f"  seq types:    {args.seq_types}")
    print(f"  k2 libs:      {args.k2_libs}")
    print(f"  threads:      {args.threads}")
    print(f"  taxonomy dir: {args.taxonomy_dir or '(auto-download)'}")
    print(f"  stage plan:   {' -> '.join(stages_to_run)}")
    print(f"\n{ckpt.summary(BuildHostPipeline.STAGES)}")

    if args.dry_run:
        print(f"\n{UI.C['cyan']}dry-run mode.{UI.C['reset']}")
        for s in stages_to_run:
            print(f"  [{'SKIP' if ckpt.is_done(s) else 'RUN'}] {s}: {BuildHostPipeline.STAGE_DESC[s]}")
        return

    p = BuildHostPipeline(args)
    funcs = {'genome-down': p.run_genome_down, 'hostdb': p.run_hostdb}
    failed = []
    t0 = time.time()
    for s in stages_to_run:
        if ckpt.is_done(s) and not args.force:
            UI.stage(BuildHostPipeline.STAGE_DESC[s], 'skip')
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
