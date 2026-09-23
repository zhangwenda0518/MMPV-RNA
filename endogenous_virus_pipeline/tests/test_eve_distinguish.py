#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_eve_distinguish.py — eve_distinguish 子模块装配测试
=========================================================
不需要真实 diamond / blastn (blastn 根本不会被调用, 因为寄主通道关闭)。

运行: python -m unittest discover -s endogenous_virus_pipeline/tests -v
  或: python endogenous_virus_pipeline/tests/test_eve_distinguish.py

锁两层:
  1. s4_filter.py 的 verdict -> action 映射。s3 的 12 种 verdict 必须全在 ACTIONS 表里,
     否则下游按列过滤会漏; 表外的新值必须按 REVIEW 兜底并返回非零, 不能静默丢行。
  2. run_all.sh 作为独立后运行脚本的装配: 参数解析 / 默认 OUT / 缺 -A 时报错 /
     模块面板复用 / 五个 python 阶段的文件交接 / -D 按发现管道目录约定定位输入。
     (这一层原先写在 virome_pipeline.py 的 run_eve_distinguish 里, EVE 判别退出主流程后
      搬到了脚本自身, 所以在这里重锁。)
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]                      # endogenous_virus_pipeline/
MOD = ROOT / 'eve_distinguish'
sys.path.insert(0, str(MOD))

from s4_filter import ACTIONS, CARRY, OUT_HDR  # noqa: E402

RUN_ALL = MOD / 'run_all.sh'
PY = sys.executable


def find_bash():
    """真的 bash 可执行文件。

    Windows 上直接调 'bash' 会被 CreateProcess 解析到 System32\\bash.exe (WSL 存根,
    这里跑不起来, 只会吐一段 UTF-16 的 wsl 提示然后 rc=1), 所以要显式找 Git 自带的那个。
    """
    if os.name != 'nt':
        return 'bash'
    cands = []
    which = shutil.which('bash')
    if which:
        cands.append(which)
    for base in (os.environ.get('ProgramFiles', r'C:\Program Files'),
                 os.environ.get('ProgramFiles(x86)', r'C:\Program Files (x86)'),
                 r'D:\Program Files'):
        cands += [str(Path(base) / 'Git' / 'usr' / 'bin' / 'bash.exe'),
                  str(Path(base) / 'Git' / 'bin' / 'bash.exe')]
    for c in cands:
        # WSL 存根也在 PATH 里, 见到 system32 就跳过
        if c and Path(c).is_file() and 'system32' not in c.lower():
            return c
    return 'bash'


BASH = find_bash()
# 传进 bash 的 PYTHON: Windows 上 sys.executable 带反斜杠, 给命令名更稳
PY_FOR_BASH = PY if os.name != 'nt' else 'python'

# 假 diamond: makedb 只按 -d 名字落一个 .dmnd 空文件, blastx 只把 -o 指向的文件建成空文件。
# 这样 run_all.sh 的建库/比对两步"成功"但没有任何命中, 正好用来验后续 python 阶段的交接。
DIAMOND_SHIM = """#!/bin/sh
out=""; db=""; prev=""
for a in "$@"; do
    [ "$prev" = "-o" ] && out="$a"
    [ "$prev" = "-d" ] && db="$a"
    prev="$a"
done
[ -n "$out" ] && : > "$out"
[ -n "$db" ] && : > "$db.dmnd"
exit 0
"""

# s3_verdict.py 写出的表头 (逐列对齐, s4 靠 DictReader 取列), 兼作列名 -> 下标索引
S3_HDR = ("contig_id\ttax_family\tcategory\tcheckv_completeness\tlocus_ncomp\tlocus_completeness\t"
          "provirus_scale\tlocus_arch\tlocus_members\tdecay_class\tstop_enrichment\t"
          "premature_stops_region\tframe_switches\torf_max_fraction\tte_flag\thost_wpid\t"
          "host_cov\tverdict\tverdict_reason\n")
S3_COL = {name: i for i, name in enumerate(S3_HDR.rstrip('\n').split('\t'))}
VERDICT_COL = S3_COL['verdict']

# s2b_locus_scan.py / run_all.sh 的 --no-host 分支写出的 s2b 表头
S2B_HDR = ('contig_id\tsample\thost_scaffold\thost_code\thost_wpid\thost_cov\thost_xeno\t'
           'comps_contig\tbest_locus_comps\tdecay_class\tte_flag\tlocus_arch\tlocus_ncomp\t'
           'locus_comps\tlocus_members\tlocus_hint\n').rstrip('\n').split('\t')

QUERY_FA = (
    ">cand_A_length_1200_cov_10.0\n" + ("ATGAAACCCGGGTTTAAACCCGGGTTTAAA" * 40) + "\n"
    ">cand_B_length_900_cov_5.0\n" + ("ATGTTTAAACCCGGGAAATTTCCCGGGAAA" * 30) + "\n"
)
EVIDENCE_TSV = ("contig_id\ttax_family\tcategory\tcheckv_completeness\n"
                "cand_A_length_1200_cov_10.0\tCaulimoviridae\tDNA_virus\t62.5\n"
                "cand_B_length_900_cov_5.0\t\tRNA_virus\t\n")


def write(path, text, mode=0o644):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding='utf-8')
    if mode:
        os.chmod(p, mode)
    return p


def _posix(p):
    """Windows 路径与 bash 拼出来的混合分隔符路径, 统一成正斜杠再比。"""
    return str(p).replace('\\', '/')


def verdict_rows(verdicts):
    """按 s3 的表头造若干行; verdict 列在最后一列 (verdict_reason 之前)。"""
    rows = [S3_HDR]
    for i, v in enumerate(verdicts):
        rows.append(f"cand_{i}\tCaulimoviridae\tDNA_virus\t50\t3\t0.6\tFALSE\t"
                    f"locus_full\tm1;m2\tdistributed_decay\t2.0\t12\t0\t0.4\tFALSE\t"
                    f"0\t0\t{v}\treason\n")
    return ''.join(rows)


def run(args, env=None, cwd=None):
    e = dict(os.environ)
    e['PYTHON'] = PY_FOR_BASH
    e['PYTHONUNBUFFERED'] = '1'
    e.update(env or {})
    return subprocess.run([BASH, str(RUN_ALL)] + args, capture_output=True, text=True,
                          encoding='utf-8', errors='replace', env=e,
                          cwd=str(cwd) if cwd else None)


class TestS4Filter(unittest.TestCase):
    """verdict -> action 映射: 12 种全覆盖 + 未知值兜底。"""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='eve_s4_'))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _filter(self, verdicts):
        src = write(self.tmp / f'verdict_{len(verdicts)}.tsv', verdict_rows(verdicts))
        out = self.tmp / 'f.tsv'
        r = subprocess.run([PY, str(MOD / 's4_filter.py'), str(src), '-o', str(out)],
                           capture_output=True, text=True, encoding='utf-8', errors='replace')
        return r, out

    def test_all_known_verdicts_map(self):
        # 只喂已知 verdict: 必须干净退出, 且一种不漏
        r, out = self._filter(sorted(ACTIONS))
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = out.read_text(encoding='utf-8').splitlines()
        self.assertEqual(lines[0].split('\t'), OUT_HDR)
        self.assertEqual(len(lines) - 1, len(ACTIONS))
        acts = [ln.split('\t')[-2] for ln in lines[1:]]
        self.assertEqual(set(acts), {'MOVE_EVE', 'KEEP_virus', 'REVIEW', 'REMOVE_host_contamination'})
        # 6 个 *_review 变体 + review 本身 = 7 行 REVIEW
        self.assertEqual(acts.count('REVIEW'), 7)

    def test_carry_columns_preserved(self):
        r, out = self._filter(sorted(ACTIONS))
        self.assertEqual(r.returncode, 0, r.stderr)
        v = sorted(ACTIONS)[1]
        src = [ln.split('\t') for ln in (self.tmp / f'verdict_{len(ACTIONS)}.tsv').read_text().splitlines()]
        row = [x for x in src[1:] if x[VERDICT_COL] == v][0]
        got = [ln.split('\t') for ln in out.read_text().splitlines() if ln.split('\t')[CARRY.index('verdict')] == v][0]
        for i, col in enumerate(CARRY):
            self.assertEqual(got[i], row[S3_COL[col]], f'列 {col} 应原样搬运')

    def test_unknown_verdict_returns_nonzero(self):
        # 表外新值: 按 REVIEW 兜底写全, 但必须返回非零, 让调用方知道 s3 加了新 verdict
        r, out = self._filter(sorted(ACTIONS) + ['brand_new_verdict_from_future_s3'])
        self.assertEqual(r.returncode, 2, '表外 verdict 必须报非零, 不能静默丢行')
        self.assertIn('brand_new_verdict_from_future_s3', r.stderr)
        lines = out.read_text(encoding='utf-8').splitlines()
        self.assertEqual(len(lines) - 1, len(ACTIONS) + 1)
        self.assertEqual(lines[-1].split('\t')[-2], 'REVIEW')


class TestRunAllCLI(unittest.TestCase):
    """参数校验: 缺输入/缺寄主根目录必须立刻失败, 不能带着跑。"""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='eve_cli_'))
        self.query = write(self.tmp / 'cand.fasta', QUERY_FA)
        self.evidence = write(self.tmp / 'rescue_evidence_scored.tsv', EVIDENCE_TSV)
        self.asm = self.tmp / 'assemblies'
        self.asm.mkdir()
        self.shim = write(self.tmp / 'diamond', DIAMOND_SHIM, mode=0o755)
        # 走到工具检查之前的用例无所谓; 真跑起来的用例要有个能用的 diamond
        self.env = {'DIAMOND': str(self.shim)}

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_help(self):
        r = run(['--help'])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn('--asm-root', r.stdout)
        self.assertIn('--no-host', r.stdout)
        self.assertIn('dna_vs_eve_filter.tsv', r.stdout)

    def test_missing_query(self):
        r = run([], env=self.env)
        self.assertEqual(r.returncode, 1)
        self.assertIn('缺少候选 fasta', r.stderr)

    def test_missing_evidence(self):
        r = run(['-q', str(self.query)], env=self.env)
        self.assertEqual(r.returncode, 1)
        self.assertIn('证据表', r.stderr)

    def test_missing_asm_root_is_fatal(self):
        # 寄主通道是唯一能识别"候选其实就是寄主基因"的手段, 缺 -A 必须报错而不是静默降级
        r = run(['-q', str(self.query), '-e', str(self.evidence)], env=self.env)
        self.assertEqual(r.returncode, 1)
        self.assertIn('asm-root', r.stderr)
        self.assertIn('host_contamination_likely', r.stderr)

    def test_empty_asm_root_dir_is_fatal(self):
        r = run(['-q', str(self.query), '-e', str(self.evidence),
                 '-A', str(self.tmp / 'nope')], env=self.env)
        self.assertEqual(r.returncode, 1)
        self.assertIn('asm-root', r.stderr)

    def test_unknown_arg(self):
        r = run(['--bogus'], env=self.env)
        self.assertEqual(r.returncode, 1)
        self.assertIn('未知参数', r.stderr)

    def test_no_host_allows_missing_asm_root(self):
        # --no-host 显式关掉寄主通道时才允许没有 -A
        r = run(['-q', str(self.query), '-e', str(self.evidence), '--no-host',
                 '-o', str(self.tmp / 'o1')], env=self.env)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)

    def test_diamond_crash_is_reported_not_silent(self):
        # set -e + 2>/dev/null 会让被 OOM killer 干掉的 diamond 静默退出; 必须留下 FAILED 行
        bad = write(self.tmp / 'diamond_crash', '#!/bin/sh\nexit 137\n', mode=0o755)
        r = run(['-q', str(self.query), '-e', str(self.evidence), '--no-host',
                 '-o', str(self.tmp / 'o2')], env={**self.env, 'DIAMOND': str(bad)})
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('FAILED', r.stderr)
        self.assertIn('DIAMOND', r.stderr)
        # 崩溃时不能留下半成品产物冒充成功
        self.assertFalse((self.tmp / 'o2' / 'dna_vs_eve_filter.tsv').exists())

    def test_python3_missing_falls_back_to_python(self):
        # Windows 上常只装 python; 等跑到第三个阶段才报 "python3: command not found" 太晚
        if shutil.which('python3'):
            self.skipTest('本机 PATH 里有 python3, 降级分支触发不了')
        r = run(['-q', str(self.query), '-e', str(self.evidence), '--no-host',
                 '-o', str(self.tmp / 'o3')],
                env={**self.env, 'PYTHON': 'python3'})
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        self.assertIn('无 python3, 改用 python', r.stderr + r.stdout)


class TestRunAllPlumbing(unittest.TestCase):
    """端到端装配 (假 diamond, 无命中): 五个 python 阶段的文件交接必须全通。"""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='eve_e2e_'))
        self.query = write(self.tmp / 'cand.fasta', QUERY_FA)
        self.evidence = write(self.tmp / 'rescue_evidence_scored.tsv', EVIDENCE_TSV)
        self.shim = write(self.tmp / 'diamond', DIAMOND_SHIM, mode=0o755)
        self.env = {'DIAMOND': str(self.shim)}

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, extra, out=None):
        o = out or self.tmp / 'out'
        r = run(['-q', str(self.query), '-e', str(self.evidence), '--no-host', '-o', str(o)]
                + extra, env=self.env)
        self.assertEqual(r.returncode, 0, r.stderr + '\n--- stdout ---\n' + r.stdout)
        return o

    def test_stage_files_handoff(self):
        out = self._run([])
        for name in ('panel.fasta', 'panel_hits.tsv', 'baits_hits.tsv', 'all_hits.tsv',
                     's1_decay.tsv', 's2_domains.tsv', 'locus_architecture.tsv',
                     'eve_distinguish_verdict.tsv', 'dna_vs_eve_filter.tsv'):
            self.assertTrue((out / name).is_file(), f'缺少 {name}')

    def test_locus_table_header_complete_when_no_host(self):
        # --no-host 时 s2b 表由 printf 拼出来; 重定向只写在最后一行会让前几列漏进 stdout,
        # 表头只剩 3 列 (s3 靠列名取值, 少一列就静默变 0)
        out = self._run([])
        hdr = (out / 'locus_architecture.tsv').read_text().splitlines()[0].split('\t')
        self.assertEqual(hdr, S2B_HDR)

    def test_module_panel_reused_not_refetched(self):
        # 模块自带面板必须被复用 (软链, 或 Windows 无符号链接权限时复制),
        # 绝不能重新联网抓取 —— build_panel.py 要 NCBI 可达
        out = self._run([])
        p = out / 'panel.fasta'
        self.assertTrue(p.is_file())
        self.assertTrue(p.is_symlink() or p.read_bytes() == (MOD / 'panel.fasta').read_bytes(),
                        '面板应复用模块产物')

    def test_filter_table_shape(self):
        out = self._run([])
        lines = (out / 'dna_vs_eve_filter.tsv').read_text(encoding='utf-8').splitlines()
        self.assertEqual(lines[0].split('\t'), OUT_HDR)
        # 无命中时 s1 给 no_hsp -> s3 走非 Cauli 分支 -> review -> REVIEW
        self.assertEqual(len(lines) - 1, QUERY_FA.count('>'))
        for ln in lines[1:]:
            self.assertEqual(ln.split('\t')[-2], 'REVIEW')

    def test_discovery_dir_resolves_inputs_and_default_out(self):
        # -D 指向一次发现管道输出: query/evidence 按阶段目录约定定位, OUT 默认 08b_EVE_Distinguish
        disc = self.tmp / 'onekp-virus'
        write(disc / '10_Reports' / 'eve_candidates.fasta', QUERY_FA)
        write(disc / '10_Reports' / 'rescue_evidence_scored.tsv', EVIDENCE_TSV)
        r = run(['-D', str(disc), '--no-host'], env=self.env)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        # bash 把 -D 的 Windows 路径和后面的 /10_Reports/... 拼成一串混合分隔符, 归一化后再比
        printed = [l for l in r.stdout.splitlines() if l.startswith('[run_all] QUERY=')]
        self.assertTrue(printed, r.stdout)
        self.assertIn(_posix(disc / '10_Reports' / 'eve_candidates.fasta'), _posix(printed[0]))
        out = disc / '08b_EVE_Distinguish'
        self.assertTrue((out / 'dna_vs_eve_filter.tsv').is_file(),
                        '-D 时 OUT 应默认 <DISCOVERY>/08b_EVE_Distinguish')

    def test_explicit_out_beats_discovery_default(self):
        disc = self.tmp / 'onekp-virus'
        write(disc / '10_Reports' / 'eve_candidates.fasta', QUERY_FA)
        write(disc / '10_Reports' / 'rescue_evidence_scored.tsv', EVIDENCE_TSV)
        mine = self.tmp / 'mine'
        r = run(['-D', str(disc), '--no-host', '-o', str(mine)], env=self.env)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        self.assertTrue((mine / 'dna_vs_eve_filter.tsv').is_file())
        self.assertFalse((disc / '08b_EVE_Distinguish').exists())

    def test_missing_baits_library_still_completes(self):
        # 假 diamond 不认 makedb 的 -d, baits.dmnd 建不出来 -> 不能因为 cat 缺文件把整个运行
        # 带崩 (真实场景: 模块只带了 panel 没带 baits, 或 makedb 被跳过)
        out = self._run(['--baits', str(self.tmp / 'no_such_baits.fa')])
        self.assertTrue((out / 'all_hits.tsv').is_file())
        self.assertTrue((out / 'dna_vs_eve_filter.tsv').is_file())


if __name__ == '__main__':
    unittest.main(verbosity=2)
