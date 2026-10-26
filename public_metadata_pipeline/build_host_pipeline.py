#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build Host Pipeline - 宿主数据库构建管道 v3.1
===============================================
三阶段: genome-down -> [eve-mask] -> hostdb
`eve-mask`（内源病毒样区域屏蔽）**默认关闭**：不给 `--eve-mask` 或任何掩码源时，
阶段计划里不会出现它，建库行为与 v3.1 完全一致。
全部参数命令行指定。

eve-mask 的上下游衔接:
  · 上游: genome-down 产出 `<work>/genome/all.genome.uniq.fasta`；
          掩码源来自 `eve_screen.py`（归档管线）的 `*.s2_verdict.tsv` / `kingdom_loci_all.tsv.gz`，
          或任意 BED；也可 `--eve-scan` 在现场跑一遍归档扫描生成判定表。
  · 输出: `<work>/genome/all.genome.masked.fasta` + `.report.tsv` + `.meta.json`
  · 下游: 掩码开启时 hostdb 默认落到 `<work>/hostdb_evemask/`（不覆盖原始 `hostdb/`），
          并把 `db_provenance.json` 写进库目录；`data_preprocessing*.py --host_db <该目录>` 直接可用
          （find_host_db 认 kraken2/ 与 bowtie2|hisat2|minimap2/ 子目录）。

用法:
  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --ncbi-api "xxx" --stage genome-down hostdb

  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --ncbi-api "xxx" --taxonomy-dir ~/database/taxonomy --stage all

  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --genome-fasta /path/to/existing.fasta --stage hostdb

  # EVE 掩码（可选臂）：用已有判定表掩码后建库，落到 hostdb_evemask/
  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --genome-fasta genome/all.genome.uniq.fasta \\
      --eve-verdict 06_EVE/kingdom_loci_all.tsv.gz --eve-pad 50 \\
      --stage eve-mask hostdb
"""

import os
import re
import sys
import time
import json
import gzip
import shutil
import argparse
import subprocess
import shlex
from datetime import datetime

from utils.pipeline_utils import UI, Checkpoint, run_cmd

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


# ==========================================
# Build Host Pipeline
# ==========================================
class BuildHostPipeline:
    STAGES = ['genome-down', 'eve-mask', 'hostdb', 'eve-report']

    STAGE_DESC = {
        'genome-down': '宿主参考基因组下载 (NCBI datasets)',
        'eve-mask':    '内源病毒样区域屏蔽 (可选, 默认关闭)',
        'hostdb':      '宿主竞争性数据库构建 (K2/BT2/HT2/MM2)',
        'eve-report':  '基因组 EVE 汇总报告 (扫描 HTML + 可选元件库/蛋白库)',
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
                    f'--filename {shlex.quote(genome_zip)} --reference --include genome,gff3,seq-report'
                )
                if self.args.ncbi_api:
                    cmd += f' --api-key {shlex.quote(self.args.ncbi_api)}'
                rc = run_cmd(cmd, s, self.log_dir, timeout=3600, secrets=self._secrets())
                if rc != 0 and self.args.download_source in ('auto', 'ngd', 'ftp') and self.args.taxid:
                    # 多通道回退: 交由 download_host_genome.py 走 ncbi-genome-download / FTP 直连
                    UI.warn("datasets failed, fallback to download_host_genome.py multi-channel...")
                    sub = (f'python {shlex.quote(self._bin("download_host_genome.py"))} '
                           f'--species {shlex.quote(self.args.species)} '
                           f'--outdir {shlex.quote(genome_dir)} '
                           f'--source {self.args.download_source} '
                           f'--taxid {self.args.taxid}')
                    if self.args.ncbi_api:
                        sub += f' --ncbi-api {shlex.quote(self.args.ncbi_api)}'
                    rc2 = run_cmd(sub, s, self.log_dir, timeout=7200, secrets=self._secrets())
                    merged = os.path.join(genome_dir, 'all.genome.uniq.fasta')
                    if rc2 == 0 and os.path.isfile(merged):
                        self._genome_fasta = merged
                        self.ckpt.mark_done(s)
                        return True
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
        # FASTA 规范不允许 # 注释行 (bowtie2-build 会拒绝), 元数据写旁车 .meta.txt
        meta_path = merged_fasta + '.meta.txt'
        with open(meta_path, 'w', encoding='utf-8') as meta:
            meta.write(f"Merged host genome: {self.args.species}\n")
            meta.write(f"Created: {datetime.now().isoformat()}\n")
        with open(merged_fasta, 'w', encoding='utf-8') as out:
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

    # ── eve-mask (可选, 默认关闭) ────────────────────
    def _eve_mask_enabled(self):
        """默认关闭；显式 --eve-mask 或给出任一掩码源时才启用。"""
        return eve_mask_enabled(self.args)

    def _raw_genome_fasta(self):
        p = getattr(self, '_genome_fasta', None)
        if p and os.path.isfile(p):
            return p
        p = self.args.genome_fasta
        if p and os.path.isfile(p):
            return p
        return os.path.join(self.dirs['genome'], 'all.genome.uniq.fasta')

    def _masked_fasta_path(self):
        return os.path.abspath(
            self.args.eve_masked_fasta
            or os.path.join(self.dirs['genome'], 'all.genome.masked.fasta'))

    def _masked_hostdb_dir(self):
        return os.path.abspath(
            self.args.masked_hostdb_dir
            or os.path.join(self.work_dir, 'hostdb_evemask'))

    @staticmethod
    def _fingerprint(path):
        """size + md5（分块）——写进 provenance，保证「这个库对应哪个 FASTA」可追溯。"""
        import hashlib
        h = hashlib.md5()
        size = 0
        with open(path, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
                size += len(chunk)
        return {'path': os.path.abspath(path), 'size_bytes': size, 'md5': h.hexdigest()}

    @staticmethod
    def _mask_advice(frac):
        """按 doc/16-host-eve-masking.md §5 的决策表给一句话建议。"""
        if frac < 0.001:
            return 'masked_frac < 0.1%：建议不掩码，写进 Limitations 即可'
        if frac <= 0.01:
            return 'masked_frac 0.1%–1%：先做 contig 级证据层（P1），掩码库当敏感性分析'
        return 'masked_frac > 1%：建议启用掩码库（P2）并在 Methods 声明'

    def resolve_inputs(self):
        """hostdb 阶段用哪个 FASTA / 写到哪个目录。

        EVE 掩码关闭（默认）→ 与旧行为完全一致：原始基因组 + `<work>/hostdb`。
        掩码开启 → 掩码基因组 + `<work>/hostdb_evemask`（可用 --masked-hostdb-dir 覆盖）。
        """
        masked = getattr(self, '_masked_genome', None)
        if not masked and self.args.masked_genome:
            masked = os.path.abspath(self.args.masked_genome)
        if not masked and self._eve_mask_enabled():
            cand = self._masked_fasta_path()
            if os.path.isfile(cand):
                masked = cand
        if self._eve_mask_enabled():
            if not masked or not os.path.isfile(masked):
                raise FileNotFoundError(
                    f"EVE 掩码已启用但找不到掩码后的基因组: {self._masked_fasta_path()}"
                    f"（先跑 `--stage eve-mask`，或用 --masked-genome 指定现成文件）")
            return masked, self._masked_hostdb_dir(), True
        return self._raw_genome_fasta(), self.dirs['hostdb'], False

    def _find_verdict_tables(self, outdir):
        """在 eve_screen 输出目录里找判定表：优先跨基因组合并表，退化为单基因组表。"""
        hits = []
        for pat in ('kingdom_loci_all.tsv.gz', 'kingdom_loci_all.tsv'):
            for root, _, files in os.walk(outdir):
                if pat in files:
                    hits.append(os.path.join(root, pat))
        if hits:
            return sorted(set(hits))
        for root, _, files in os.walk(outdir):
            for f in files:
                if f.endswith('.s2_verdict.tsv'):
                    hits.append(os.path.join(root, f))
        return sorted(hits)

    def _run_eve_scan(self, genome_fasta):
        """可选：调归档管线 eve_screen.py 现场扫一遍，返回判定表路径列表。

        归档口径是「方法不成熟」（阈值/污染-退化区分未验证），这里只当**掩码源生成器**用，
        不当 EVE 判定结论。
        """
        repo_root = os.path.dirname(os.path.abspath(self.bin_dir))
        if self.args.eve_scan_script:
            candidates = [self.args.eve_scan_script]
        else:
            # 现行位置: 本管线自己的 eve_genome/；旧归档仅作兜底
            candidates = [
                os.path.join(self.bin_dir, 'eve_genome', 'eve_screen.py'),
                os.path.join(repo_root, 'archive', '20261011_endogenous_virus_pipeline',
                             'eve_screen.py'),
            ]
        script = next((p for p in candidates if os.path.isfile(p)), candidates[0])
        if not os.path.isfile(script):
            UI.err("eve_screen.py 不存在，已尝试: "
                   + " , ".join(candidates) + "（用 --eve-scan-script 指定）")
            return []
        outdir = os.path.join(self.work_dir, 'eve')
        name = re.sub(r'[^0-9A-Za-z_.-]+', '_', self.args.species).strip('_') or 'host'
        UI.warn("EVE 现场扫描 (--eve-scan)：耗时较长；结果仅作掩码源，不作结论")
        rc = run_cmd(
            f'python {shlex.quote(script)} -g {shlex.quote(genome_fasta)} '
            f'-n {shlex.quote(name)} --outdir {shlex.quote(outdir)} '
            f'--stage 1,2,3 --fast -t {self.args.threads}',
            'eve-mask', self.log_dir, timeout=86400 * 2)
        if rc != 0:
            UI.err("eve_screen 扫描失败（可先手工跑 `--stage` 子步骤排查）")
            return []
        run_cmd(f'python {shlex.quote(script)} --merge -o {shlex.quote(outdir)}',
                'eve-mask', self.log_dir, timeout=86400)
        tables = self._find_verdict_tables(outdir)
        if tables:
            UI.ok(f"eve_screen 判定表: {len(tables)} 个（{os.path.dirname(tables[0])}）")
        else:
            UI.warn(f"eve_screen 输出里没找到判定表: {outdir}")
        return tables

    def run_eve_mask(self):
        s = 'eve-mask'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)

        if not self._eve_mask_enabled():
            UI.warn("EVE 掩码未启用（默认关闭）：跳过，hostdb 使用原始宿主基因组")
            self.ckpt.mark_done(s)
            return True

        # 1) 现成的掩码基因组：只登记，不重算
        if self.args.masked_genome:
            fasta = os.path.abspath(self.args.masked_genome)
            if not os.path.isfile(fasta):
                UI.err(f"missing masked genome: {fasta}")
                self.ckpt.mark_fail(s, "missing masked genome")
                return False
            self._masked_genome = fasta
            meta = {
                'schema': 'mmpv-eve-mask/1',
                'created_at': datetime.now().isoformat(timespec='seconds'),
                'mode': 'provided',
                'masked_genome': self._fingerprint(fasta),
                'mask_sources': {'provided': fasta},
                'params': {},
                'totals': {},
                'advice': '外部提供的掩码基因组：掩码比例未知，请自行确认来源',
            }
            meta_path = fasta + '.meta.json'
            with open(meta_path, 'w', encoding='utf-8') as fh:
                json.dump(meta, fh, ensure_ascii=False, indent=2)
            UI.ok(f"使用外部掩码基因组: {fasta}")
            self.ckpt.mark_done(s)
            return True

        src = self._raw_genome_fasta()
        if not os.path.isfile(src):
            UI.err(f"missing genome fasta: {src}")
            self.ckpt.mark_fail(s, "missing genome FASTA")
            return False

        bed = [os.path.abspath(p) for p in (self.args.eve_bed or [])]
        verdict_tables = [os.path.abspath(p) for p in (self.args.eve_verdict or [])]
        if self.args.eve_scan and not verdict_tables:
            verdict_tables = self._run_eve_scan(src)
        if not bed and not verdict_tables:
            UI.err("eve-mask 需要掩码源: --eve-bed / --eve-verdict / --masked-genome（或 --eve-scan）")
            self.ckpt.mark_fail(s, "no mask source")
            return False
        missing = [p for p in list(bed) + list(verdict_tables) if not os.path.isfile(p)]
        if missing:
            UI.err(f"mask source 不存在: {missing}")
            self.ckpt.mark_fail(s, "missing mask source files")
            return False

        # 2) 解析区间（BED 0 基半开 + verdict 1 基闭区间 → 统一成 0 基半开）
        try:
            from mask_host_genome import parse_bed, parse_verdict_tsv, mask_genome
        except ImportError:  # 以包/别处 cwd 导入时，按文件路径加载同目录模块
            import importlib.util
            _mp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mask_host_genome.py')
            _spec = importlib.util.spec_from_file_location('mask_host_genome', _mp)
            _mod = importlib.util.module_from_spec(_spec)
            _spec.loader.exec_module(_mod)
            parse_bed, parse_verdict_tsv, mask_genome = \
                _mod.parse_bed, _mod.parse_verdict_tsv, _mod.mask_genome

        vals = [v.strip() for v in (self.args.eve_verdict_values or 'viral_supported').split(',')
                if v.strip()]
        ivs = {}
        if bed:
            b_ivs, n_bed, n_bed_skip = parse_bed(bed)
            UI.info(f"BED 区间: {n_bed} 段 (跳过 {n_bed_skip} 行)")
            for contig, segs in b_ivs.items():
                ivs.setdefault(contig, []).extend(segs)
        if verdict_tables:
            v_ivs, n_v, n_v_skip = parse_verdict_tsv(verdict_tables, vals)
            UI.info(f"EVE 判定表: {n_v} 位点 (verdict∈{{{','.join(vals)}}}, 跳过 {n_v_skip})")
            for contig, segs in v_ivs.items():
                ivs.setdefault(contig, []).extend(segs)
        if not ivs:
            UI.err("掩码源里没有可用区间：检查 verdict/BED 的 contig 名是否与 FASTA 头一致、"
                   "verdict 值是否在 --eve-verdict-values 里")
            self.ckpt.mark_fail(s, "empty mask intervals")
            return False

        # 3) 掩码（长度不变）+ 报告 + 溯源
        out = self._masked_fasta_path()
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        report_path = out + '.report.tsv'
        rows = []
        summary = mask_genome(
            src, ivs, out, rows, os.path.splitext(os.path.basename(src))[0],
            pad=self.args.eve_pad, merge_gap=self.args.eve_merge_gap,
            mode=self.args.eve_mask_mode, line_width=60)
        with open(report_path, 'w', encoding='utf-8') as fh:
            fh.write("genome\tcontig\tlength\tn_loci\tmasked_bp\tmasked_frac\n")
            for row in rows:
                fh.write("\t".join(row) + "\n")

        frac = float(summary.get('masked_frac', 0.0))
        if int(summary.get('n_loci', 0)) == 0 or int(summary.get('n_masked_bp', 0)) == 0:
            # 区间解析出来了、但一条也没落到 FASTA 的 contig 上（contig 名不匹配是最常见原因）
            # → 宁可报错，也不要留下一个「其实没掩码」的 masked FASTA 冒充掩码库。
            try:
                os.remove(out)
            except OSError:
                pass
            UI.err("掩码区间没有落在任何 contig 上（verdict/BED 的 contig 名与 FASTA 头是否一致？）"
                   " → 拒绝生成掩码库")
            self.ckpt.mark_fail(s, "no interval matched any contig")
            return False
        meta = {
            'schema': 'mmpv-eve-mask/1',
            'created_at': datetime.now().isoformat(timespec='seconds'),
            'mode': 'scan' if self.args.eve_scan else 'tables',
            'source_genome': dict(self._fingerprint(src),
                                  contigs=int(summary.get('n_contigs', 0))),
            'masked_genome': self._fingerprint(out),
            'mask_sources': {
                'bed': bed,
                'verdict': verdict_tables,
                'verdict_values': vals,
                'eve_scan': bool(self.args.eve_scan),
            },
            'params': {
                'pad': self.args.eve_pad,
                'merge_gap': self.args.eve_merge_gap,
                'mask_mode': self.args.eve_mask_mode,
                'line_width': 60,
                'intervals_after_parse': sum(len(v) for v in ivs.values()),
                'contigs_with_intervals': len(ivs),
            },
            'totals': {
                'bp': int(summary.get('n_bp', 0)),
                'masked_bp': int(summary.get('n_masked_bp', 0)),
                'loci': int(summary.get('n_loci', 0)),
                'masked_frac': round(frac, 6),
            },
            'report': os.path.abspath(report_path),
            'advice': self._mask_advice(frac),
        }
        with open(out + '.meta.json', 'w', encoding='utf-8') as fh:
            json.dump(meta, fh, ensure_ascii=False, indent=2)

        self._masked_genome = out
        UI.ok(f"掩码完成: {int(summary['n_masked_bp']):,} bp ({frac * 100:.3f}%) / "
              f"{int(summary['n_loci'])} 段 → {out}")
        UI.warn(meta['advice'])
        self.ckpt.mark_done(s)
        return True

    def write_db_provenance(self, db_dir, genome_fasta, masked):
        """把「这个宿主库是哪个 FASTA、掩没掩、用哪些工具/参数建的」写进库目录。"""
        meta = {
            'schema': 'mmpv-hostdb-provenance/1',
            'created_at': datetime.now().isoformat(timespec='seconds'),
            'pipeline': 'build_host_pipeline.py v3.1',
            'species': self.args.species,
            'taxid': self.args.taxid,
            'genome': self._fingerprint(genome_fasta),
            'masked': bool(masked),
            'eve_mask_meta': None,
            'hostdb_tools': self.args.hostdb_tools,
            'seq_types': self.args.seq_types,
            'k2_libs': self.args.k2_libs,
            'taxonomy_dir': self.args.taxonomy_dir,
            'threads': self.args.threads,
            'db_dir': os.path.abspath(db_dir),
        }
        mask_meta_path = genome_fasta + '.meta.json'
        if masked and os.path.isfile(mask_meta_path):
            try:
                with open(mask_meta_path, 'r', encoding='utf-8') as fh:
                    meta['eve_mask_meta'] = json.load(fh)
            except Exception as e:
                meta['eve_mask_meta'] = {'error': str(e)}
        out = os.path.join(db_dir, 'db_provenance.json')
        os.makedirs(os.path.dirname(os.path.abspath(out)) or '.', exist_ok=True)
        with open(out, 'w', encoding='utf-8') as fh:
            json.dump(meta, fh, ensure_ascii=False, indent=2)
        return out

    # ── hostdb ─────────────────────────────────────
    def run_eve_report(self):
        """[eve-report] 基因组 EVE 汇总: 扫描 HTML (+ 可选元件库/蛋白库).

        位置在 hostdb 之后: 掩码发生在 hostdb 之前 (保证宿主库不含 EVE 序列),
        这里只做"分析结果汇总", 不改动已建好的宿主库。
        """
        s = 'eve-report'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)

        eve_dir = os.path.join(self.work_dir, 'eve')
        if not os.path.isdir(eve_dir):
            UI.warn("未发现 EVE 扫描产物 (--eve-scan 未启用或 eve-mask 跳过): "
                    "跳过 EVE 报告")
            self.ckpt.mark_done(s)
            return True

        manifest = {'eve_root': eve_dir, 'scan_html': None, 'summary': None,
                    'consensus_library': None, 'protein_library': None}

        # 可选: 元件层 (聚类/consensus/CDD) + 六框蛋白库
        if getattr(self.args, 'eve_elements', False):
            elem_script = os.path.join(self.bin_dir, 'eve_genome', 'eve_element_db.py')
            if not os.path.isfile(elem_script):
                UI.err(f"eve_element_db.py 不存在: {elem_script}")
                self.ckpt.mark_fail(s, 'missing eve_element_db.py')
                return False
            cmd = (f'python {shlex.quote(elem_script)} -o {shlex.quote(eve_dir)} '
                   f'--stages {shlex.quote(self.args.eve_element_stages)} '
                   f'-t {self.args.threads}')
            if self.args.eve_diamond:
                cmd += f' --diamond {shlex.quote(self.args.eve_diamond)}'
            if self.args.eve_ictv_prot:
                cmd += f' --ictv-prot {shlex.quote(self.args.eve_ictv_prot)}'
            if run_cmd(cmd, s, self.log_dir, timeout=86400 * 3) != 0:
                UI.warn("元件层失败: 扫描报告仍会生成, 但不产出 EVE 参考库")
            lib = os.path.join(eve_dir, '05_Elements', 'eve_consensus_library.fa')
            if os.path.isfile(lib):
                manifest['consensus_library'] = lib
                prot_script = os.path.join(self.bin_dir, 'eve_genome', 'eve_prot_lib.py')
                faa = os.path.join(eve_dir, '05_Elements',
                                   'eve_consensus_library.sixframe.faa')
                if os.path.isfile(prot_script) and run_cmd(
                        f'python {shlex.quote(prot_script)} --in {shlex.quote(lib)} '
                        f'--out {shlex.quote(faa)}', s, self.log_dir,
                        timeout=86400) == 0 and os.path.isfile(faa):
                    manifest['protein_library'] = faa

        # 扫描 HTML
        rep_script = os.path.join(self.bin_dir, 'eve_genome', 'eve_scan_report.py')
        if not os.path.isfile(rep_script):
            UI.err(f"eve_scan_report.py 不存在: {rep_script}")
            self.ckpt.mark_fail(s, 'missing eve_scan_report.py')
            return False
        html = os.path.join(eve_dir, '04_Summary', 'eve_report.html')
        if run_cmd(f'python {shlex.quote(rep_script)} -o {shlex.quote(eve_dir)} '
                   f'--out {shlex.quote(html)}', s, self.log_dir, timeout=3600) != 0:
            UI.err("EVE HTML 报告生成失败")
            self.ckpt.mark_fail(s, 'eve report failed')
            return False
        if os.path.isfile(html):
            manifest['scan_html'] = html
        summ = os.path.join(eve_dir, '04_Summary', 'kingdom_summary.tsv')
        if os.path.isfile(summ):
            manifest['summary'] = summ

        man_path = os.path.join(eve_dir, 'eve_manifest.json')
        with open(man_path, 'w', encoding='utf-8') as fh:
            json.dump(manifest, fh, ensure_ascii=False, indent=2)
        UI.ok(f"EVE 报告: {html}")
        if manifest['protein_library']:
            UI.ok(f"EVE 蛋白参考库 (供发现管线 Δbit-score): "
                  f"{manifest['protein_library']}")
        self.dirs['eve'] = eve_dir
        self.ckpt.mark_done(s)
        return True

    def run_hostdb(self):
        s = 'hostdb'
        UI.stage(self.STAGE_DESC[s], 'start')
        self.ckpt.mark_start(s)

        # 查找 FASTA 与库目录（EVE 掩码关闭时 = 原始基因组 + <work>/hostdb）
        try:
            genome_fasta, hostdb_dir, masked = self.resolve_inputs()
        except FileNotFoundError as e:
            UI.err(str(e))
            self.ckpt.mark_fail(s, "missing genome FASTA")
            return False
        if not os.path.isfile(genome_fasta):
            UI.err(f"missing genome fasta: {genome_fasta}")
            self.ckpt.mark_fail(s, "missing genome FASTA")
            return False
        os.makedirs(hostdb_dir, exist_ok=True)
        if masked:
            UI.warn(f"使用 EVE 掩码后的宿主基因组建库 → {hostdb_dir}")
        else:
            UI.info(f"使用原始宿主基因组建库 → {hostdb_dir}")

        cmd_parts = [
            f'python {shlex.quote(self._bin("build_hostbase.py"))}',
            f'--tool {shlex.quote(self.args.hostdb_tools)}',
            f'--input {shlex.quote(genome_fasta)}',
            f'--output {shlex.quote(hostdb_dir)}',
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
            try:
                prov = self.write_db_provenance(hostdb_dir, genome_fasta, masked)
                UI.ok(f"db provenance: {prov}")
            except Exception as e:
                UI.warn(f"provenance 写入失败（不影响建库）: {e}")
            if masked:
                UI.ok(f"masked hostdb: {hostdb_dir}")
                self._hostdb_hint = hostdb_dir
            self.ckpt.mark_done(s)
            return True
        self.ckpt.mark_fail(s, f"exit={rc}")
        return False


def build_parser():
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

  # EVE 掩码臂（可选，默认关闭）：掩码后建库，落到 hostdb_evemask/
  python build_host_pipeline.py --species "Lycium barbarum" --taxid 112863 \\
      --genome-fasta genome/all.genome.uniq.fasta \\
      --eve-verdict 06_EVE/kingdom_loci_all.tsv.gz --eve-pad 50 \\
      --stage eve-mask hostdb
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
    parser.add_argument('--download-source', default='auto',
                        choices=['auto', 'datasets', 'ngd', 'ftp'],
                        help='基因组下载通道: auto=datasets 失败自动回退 ncbi-genome-download / FTP 直连')
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

    g = parser.add_argument_group('EVE/EPRV 掩码 (可选, 默认关闭)')
    g.add_argument('--eve-mask', action='store_true',
                   help='启用 eve-mask 阶段；也可由任一掩码源 (--eve-bed/--eve-verdict/'
                        '--masked-genome/--eve-scan) 自动启用')
    g.add_argument('--eve-bed', nargs='*', default=[],
                   help='掩码区间 BED (0 基半开: contig start end)，可多个取并集')
    g.add_argument('--eve-verdict', nargs='*', default=[],
                   help='EVE 判定表 (*.s2_verdict.tsv / kingdom_loci_all.tsv.gz)，可多个')
    g.add_argument('--eve-verdict-values', default='viral_supported',
                   help='判定表中要掩码的 verdict，逗号分隔 (default: viral_supported；'
                        'viral_supported,undetermined = 激进臂)')
    g.add_argument('--eve-scan', action='store_true',
                   help='没给判定表时，现场调归档 eve_screen.py 扫描生成掩码源 (慢，默认关)')
    g.add_argument('--eve-scan-script', help='eve_screen.py 路径 (默认自动定位: archive/20261011_* 或顶层 endogenous_virus_pipeline/)')
    g.add_argument('--eve-pad', type=int, default=0,
                   help='区间两侧外扩 bp (default: 0；建议 25-50)')
    g.add_argument('--eve-merge-gap', type=int, default=0,
                   help='相邻区间合并的最大间隙 bp (default: 0)')
    g.add_argument('--eve-mask-mode', choices=['N', 'lower'], default='N',
                   help='N=替换成 N (建库首选, default); lower=软掩码')
    g.add_argument('--eve-masked-fasta',
                   help='掩码后 FASTA 输出 (默认 <work>/genome/all.genome.masked.fasta)')
    g.add_argument('--masked-genome',
                   help='已掩码好的宿主 FASTA：跳过掩码计算，只登记 provenance')
    g.add_argument('--masked-hostdb-dir',
                   help='掩码库输出目录 (默认 <work>/hostdb_evemask，不覆盖原始 hostdb/)')
    g = parser.add_argument_group('EVE 分析汇总 (eve-report 阶段)')
    g.add_argument('--eve-elements', action='store_true',
                   help='扫描后继续跑元件层 (聚类/consensus/CDD)，产出 EVE 参考库；'
                        '并自动六框翻译成蛋白库供发现管线 Δbit-score 使用')
    g.add_argument('--eve-element-stages', default='collect,cluster,consensus,cdd',
                   help='元件层阶段 (default: collect,cluster,consensus,cdd)')
    g.add_argument('--eve-diamond', default='',
                   help='元件层 DIYMOND 可执行文件 (Other_viral 重分类用)')
    g.add_argument('--eve-ictv-prot', default='',
                   help='ICTV 六框蛋白库 (Other_viral 重分类用)')
    return parser


def eve_mask_enabled(args):
    """eve-mask 默认关闭：显式 --eve-mask 或给出任一掩码源才启用。"""
    return bool(getattr(args, 'eve_mask', False)
                or getattr(args, 'masked_genome', None)
                or getattr(args, 'eve_bed', None)
                or getattr(args, 'eve_verdict', None)
                or getattr(args, 'eve_scan', False))


def plan_stages(args):
    """决定本轮跑哪些 stage；eve-mask 未启用时从 `--stage all` 里移除（不写 checkpoint）。"""
    if 'all' in args.stage:
        stages = list(BuildHostPipeline.STAGES)
    else:
        stages = [s for s in args.stage if s in BuildHostPipeline.STAGES]
    explicit = 'all' not in args.stage
    if 'eve-mask' in stages and not eve_mask_enabled(args):
        if explicit:
            raise SystemExit(
                "[ERROR] --stage eve-mask 需要掩码源: --eve-mask / --eve-bed / "
                "--eve-verdict / --masked-genome / --eve-scan")
        stages.remove('eve-mask')
    return stages


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.ncbi_api:
        args.ncbi_api = os.environ.get('NCBI_API_KEY', '')
    if not args.bin_dir:
        args.bin_dir = os.path.dirname(os.path.abspath(__file__))

    UI.banner("Build Host Pipeline v3.1")

    stages_to_run = plan_stages(args)

    ckpt = Checkpoint(os.path.abspath(args.work_dir))
    if args.force:
        ckpt.reset()

    print(f"  species:      {args.species} (taxid: {args.taxid})")
    print(f"  work dir:     {os.path.abspath(args.work_dir)}")
    print(f"  ncbi api:     {'configured' if args.ncbi_api else 'NOT set'}")
    print(f"  genome fasta: {args.genome_fasta or '(from datasets CLI)'}")
    eve_on = eve_mask_enabled(args)
    if eve_on:
        srcs = []
        if args.eve_bed:
            srcs.append(f'BED×{len(args.eve_bed)}')
        if args.eve_verdict:
            srcs.append(f'verdict×{len(args.eve_verdict)}')
        if args.eve_scan:
            srcs.append('eve-scan')
        if args.masked_genome:
            srcs.append('masked-genome')
        print(f"  eve mask:     ON ({', '.join(srcs) or '--eve-mask'})")
        print(f"  mask values:  {args.eve_verdict_values}  pad={args.eve_pad} "
              f"gap={args.eve_merge_gap} mode={args.eve_mask_mode}")
    else:
        print(f"  eve mask:     OFF (默认关闭；加 --eve-mask 或给掩码源打开)")
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
    funcs = {'genome-down': p.run_genome_down,
             'eve-mask': p.run_eve_mask,
             'hostdb': p.run_hostdb,
             'eve-report': p.run_eve_report}
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
        hint = getattr(p, '_hostdb_hint', None)
        if hint:
            print(f"    masked hostdb: {hint}")
            print(f"    downstream:    --host_db {hint}")
            print(f"    下游用法:      python data_preprocessing_pipeline/data_preprocessing.py "
                  f"--host_db {hint} ...  (preprocess_unified.py 同理)")


if __name__ == '__main__':
    main()
