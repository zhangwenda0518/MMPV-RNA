#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_eve_element_db.py — 元件级后处理层 (eve_element_db.py) 测试
================================================================
离线: 不需要真 mmseqs / mafft / RepeatMasker —— 三个假工具是写进临时 bin/ 的小程序
(Windows 落 .cmd 启动器, POSIX 落 sh 包装), 借助 PATH 被解析, 与
test_eve_distinguish.py 里假 diamond 的做法同一套思路。

三层:
  1. 纯逻辑 (进程内): locus 解析 / 元件区归并 / 候选筛选 / consensus 输入装配 /
     比对列投票 / 聚类聚合 / 批量表解析 / RepeatMasker .out 解析。
  2. 端到端 (subprocess 跑真 CLI): 假工具跑完 collect→cluster→regions→consensus→rm,
     逐条核对表头常量、产物文件、[done] 行、以及工具收到的参数 (防"悄悄改了参数")。
     输入不是手抄的 summary, 而是交给 eve_scan_core.summarize_genome() 生成的 ——
     抄一份表头就只能证明抄对了, 证明不了上游真这么写。
  3. 失败可见性: 工具崩掉必须留 [fail] 行且不留半成品表; 需要额外依赖的阶段缺参数
     必须在开跑前报出来 (源脚本 eve_kingdom_post.py 是三处静默跳过, 这里锁住"不许静默")。

运行: cd endogenous_virus_pipeline && python -m pytest tests/ -q
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
sys.path.insert(0, str(ROOT))

import eve_element_db as ed                 # noqa: E402
from eve_scan_core import summarize_genome  # noqa: E402

MOD = ROOT / "eve_element_db.py"
PY = sys.executable


# ─────────────────────── 假外部工具 ───────────────────────
MMSEQS_SHIM = r'''"""假 mmseqs: 只实现元件层用到的四条子命令, 调用记录写进 $EVE_FAKE_LOG."""
import os, pathlib, sys

a = sys.argv[1:]
if os.environ.get("EVE_FAKE_LOG"):
    with open(os.environ["EVE_FAKE_LOG"], "a", encoding="utf-8") as fh:
        fh.write("mmseqs " + " ".join(a) + "\n")
cmd = a[0] if a else ""
if cmd == "easy-linclust":
    prefix, tmp = a[2], a[3]
    pathlib.Path(tmp).mkdir(parents=True, exist_ok=True)
    src = os.environ.get("EVE_FAKE_CLUSTERS", "")
    body = pathlib.Path(src).read_text(encoding="utf-8") if src else ""
    with open(prefix + "_cluster.tsv", "w", encoding="utf-8", newline="\n") as fo:
        fo.write(body)
elif cmd == "createdb":
    pathlib.Path(a[2]).mkdir(parents=True, exist_ok=True)
elif cmd == "search":
    pathlib.Path(a[3]).mkdir(parents=True, exist_ok=True)
elif cmd == "convertalis":
    with open(a[4], "w", encoding="utf-8", newline="\n") as fo:
        fo.write(os.environ.get("EVE_FAKE_DOMAINS", ""))
else:
    sys.stderr.write("假 mmseqs 不认识的子命令: %r\n" % cmd)
    sys.exit(2)
'''

MAFFT_SHIM = r'''"""假 mafft: 把每条序列右侧补 '-' 到等长后写出 (只求列数一致, 好验投票逻辑)."""
import os, sys

a = sys.argv[1:]
inp, out, i = None, None, 0
while i < len(a):
    if a[i] == "-o":
        out = a[i + 1]; i += 2; continue
    if a[i] == "--thread":
        i += 2; continue
    if a[i].startswith("-"):
        i += 1; continue
    inp = a[i]; i += 1
if os.environ.get("EVE_FAKE_LOG"):
    with open(os.environ["EVE_FAKE_LOG"], "a", encoding="utf-8") as fh:
        fh.write("mafft " + " ".join(a) + "\n")
recs, name, buf = [], None, []
with open(inp, encoding="utf-8") as fh:
    for line in fh:
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
if name:
    recs.append((name, "".join(buf)))
w = max((len(s) for _, s in recs), default=0)
with open(out, "w", encoding="utf-8", newline="\n") as fo:
    for n, s in recs:
        fo.write(">%s\n%s\n" % (n, s + "-" * (w - len(s))))
'''

# -lib 里的每个成员给每个基因组造一条命中; 两行表头与一行说明必须被解析器跳过
RM_SHIM = r'''"""假 RepeatMasker: 按 -lib 成员为每条 contig 造一条 .out 命中."""
import os, pathlib, sys

a = sys.argv[1:]
lib, rdir, genome, i = "", "", "", 0
while i < len(a):
    if a[i] in ("-lib", "-dir", "-pa"):
        val = a[i + 1]
        if a[i] == "-lib":
            lib = val
        elif a[i] == "-dir":
            rdir = val
        i += 2; continue
    if a[i].startswith("-"):
        i += 1; continue
    genome = a[i]; i += 1
if os.environ.get("EVE_FAKE_LOG"):
    with open(os.environ["EVE_FAKE_LOG"], "a", encoding="utf-8") as fh:
        fh.write("RepeatMasker " + " ".join(a) + "\n")
names = []
with open(lib, encoding="utf-8") as fh:
    for line in fh:
        if line.startswith(">"):
            names.append(line[1:].strip().split()[0])
contig = "chr1"
with open(genome, encoding="utf-8") as fh:
    for line in fh:
        if line.startswith(">"):
            contig = line[1:].strip().split()[0]; break
pathlib.Path(rdir).mkdir(parents=True, exist_ok=True)
out = pathlib.Path(rdir) / (pathlib.Path(genome).name + ".out")
with open(out, "w", encoding="utf-8", newline="\n") as fo:
    fo.write("   SW   perc perc perc  query     position in query              "
             "matching                    repeat        position in repeat\n")
    fo.write("score   div. del. ins.  sequence  begin    end          (left)   "
             "repeat                      class/family begin  end    (left)   ID\n\n")
    fo.write("There were no repetitive sequences detected in the query sequence\n")
    for k, nm in enumerate(names):
        b = 100 + k * 1000
        # (left) 列按真 .out 的写法 '(0)' (无内空格), 一行 15 列
        fo.write("  1234   5.0  0.0  0.0  %-8s  %8d %8d (0) + %s  Test/Test "
                 "   1   200  (0)   %d\n" % (contig, b, b + 199, nm, k + 1))
'''

CRASH_SHIM = 'import sys\nsys.stderr.write("假工具: 被 OOM killer 干掉了\\n")\nsys.exit(137)\n'


def write_shim(bindir, name, code):
    """在 bindir 里落一个可跨平台的假工具, 返回可执行文件路径.

    Windows 上 CreateProcess 跑不了裸 .py (不会走文件关联), PATH 查找也只认 PATHEXT
    里的后缀 —— 所以两边都落: <name>.py 干活, <name>.cmd / <name> 只负责转交。
    """
    bindir = Path(bindir)
    bindir.mkdir(parents=True, exist_ok=True)
    py = bindir / f"_{name}_shim.py"
    with open(py, "w", encoding="utf-8", newline="\n") as fo:
        fo.write(code)
    if os.name == "nt":
        exe = bindir / f"{name}.cmd"
        with open(exe, "w", encoding="utf-8", newline="\n") as fo:
            fo.write(f'@echo off\r\n"{sys.executable}" "{py}" %*\r\n')
    else:
        exe = bindir / name
        with open(exe, "w", encoding="utf-8", newline="\n") as fo:
            fo.write(f'#!/bin/sh\nexec "{sys.executable}" "{py}" "$@"\n')
        os.chmod(exe, 0o755)
    return exe


# ─────────────────────── 合成输入 ───────────────────────
# 表头由 eve_scan_core.summarize_genome() 生成 (见 build_screen_output), 不在这里抄
G1_ROWS = [
    # locus, verdict, ref_family, ref_sseqid, ref_bitscore, rvdb_stitle, rvdb_bitscore
    ("chr1:100-300", "viral_supported", "Caulimoviridae", "REF_A", "120",
     "Cauliflower mosaic virus", "200"),
    ("chr1:5000-5300", "viral_supported", "Caulimoviridae", "REF_A", "90", "", ""),
    ("chr1:20000-20200", "host_like", "Caulimoviridae", "REF_A", "300", "", ""),
    # 证据不足 (max(10, 空) = 10 < min_bs 50): collect 不要, regions 要
    ("chr1:30000-30200", "viral_supported", "Caulimoviridae", "REF_A", "10", "", ""),
    ("chr2:100-200", "viral_supported", "Geminiviridae", "REF_B", "150", "", ""),
]
G2_ROWS = [
    ("scaf_1:1000-1200", "viral_supported", "Caulimoviridae", "REF_A", "130", "", ""),
    ("scaf_1:60000-60200", "viral_supported", "Geminiviridae", "REF_B", "20",
     "Tomato yellow leaf curl virus", "80"),
]
CORE = "ATGAAACCCGGGTTTAAACC"
SEQS = {
    "chr1:100-300": CORE + "AAAA" + "G" * 30,
    "chr1:5000-5300": CORE + "AAAA" + "G" * 20,
    "chr1:20000-20200": "CCCC" * 12,
    "chr1:30000-30200": CORE + "TTTT",
    "chr2:100-200": "TTTTGGGGCCCCAAAATTTT" * 3,
    "scaf_1:1000-1200": CORE + "AAAA" + "G" * 10,
    "scaf_1:60000-60200": "TTTTGGGGCCCCAAAATTTT" * 4,
}
# 假 mmseqs 的聚类结果: C000000 = Caulimoviridae 3 簇成员, C000001 = Geminiviridae 2 个
FAKE_CLUSTERS = (
    "G1__chr1:100-300\tG1__chr1:100-300\n"
    "G1__chr1:100-300\tG1__chr1:5000-5300\n"
    "G1__chr1:100-300\tG2__scaf_1:1000-1200\n"
    "G1__chr2:100-200\tG1__chr2:100-200\n"
    "G1__chr2:100-200\tG2__scaf_1:60000-60200\n"
)
FAKE_DOMAINS = "C000000_Caulimoviridae\tPF00078\t1e-30\t99.5\t2\t120\t5\t123\t150\n"


def build_screen_output(out, *, underscore_headers=False):
    """造一份"跑完的 eve_screen 输出", 尽量走真上游代码.

    每个基因组: 写 01_Loci/<N>/<N>.loci.bed + .s1_best.tsv, 02_Verdict 的 s2_verdict,
    03_RVDB 的 s3_rvdb, 再调 summarize_genome() 生成 04_Summary/<N>_eve_summary.tsv ——
    本层的输入契约因此与上游实现绑定, 上游改列名这里会立刻失败。
    """
    out = Path(out)
    for name, rows in (("G1", G1_ROWS), ("G2", G2_ROWS)):
        loci_d = out / "01_Loci" / name
        verd_d = out / "02_Verdict" / name
        rvdb_d = out / "03_RVDB" / name
        for d in (loci_d, verd_d, rvdb_d):
            d.mkdir(parents=True, exist_ok=True)
        with open(loci_d / f"{name}.loci.bed", "w", encoding="utf-8",
                  newline="\n") as fo:
            for locus, *_rest in rows:
                contig, span = locus.rsplit(":", 1)
                s, e = span.split("-")
                fo.write(f"{contig}\t{int(s) - 1}\t{e}\t{locus}\n")
        with open(loci_d / f"{name}.s1_best.tsv", "w", encoding="utf-8",
                  newline="\n") as fo:
            fo.write("locus\tref_sseqid\tfamily\tref_evalue\tref_bitscore\tref_qcov\n")
            for locus, _v, fam, sid, bs, _st, _rb in rows:
                key = ed.locus_key(locus) if underscore_headers else locus
                fo.write(f"{key}\t{sid}\t{fam}\t1e-20\t{bs}\t80\n")
        with open(verd_d / f"{name}.s2_verdict.tsv", "w", encoding="utf-8",
                  newline="\n") as fo:
            fo.write("locus\tverdict\tbest_viral_id\tbest_viral_bs\t"
                     "best_plant_id\tbest_plant_bs\n")
            for locus, v, _f, sid, bs, _st, _rb in rows:
                key = ed.locus_key(locus) if underscore_headers else locus
                fo.write(f"{key}\t{v}\t{sid}\t{bs}\tHOST1\t5\n")
        with open(rvdb_d / f"{name}.s3_rvdb.tsv", "w", encoding="utf-8",
                  newline="\n") as fo:
            for locus, _v, _f, sid, _bs, st, rb in rows:
                key = ed.locus_key(locus) if underscore_headers else locus
                if rb:
                    # diamond OUTFMT2: query sseqid pident length qstart qend
                    # sstart send evalue bitscore stitle
                    fo.write(f"{key}\tRVDB_{sid}\t0.9\t300\t1\t300\t10\t310\t"
                             f"1e-40\t{rb}\t{st}\n")
        summarize_genome(name, out)
    # loci.fa: Stage1 的 samtools faidx 产物 (头即区域串)
    for name, rows in (("G1", G1_ROWS), ("G2", G2_ROWS)):
        p = out / "01_Loci" / name / f"{name}.loci.fa"
        with open(p, "w", encoding="utf-8", newline="\n") as fo:
            for locus, *_rest in rows:
                head = ed.locus_key(locus) if underscore_headers else locus
                fo.write(f">{head}\n{SEQS[locus]}\n")
    return out


def build_batch(tmp, out):
    """rm 阶段的批量表 + 两个小基因组文件 (假 RepeatMasker 只读头一行取 contig 名)."""
    lines = []
    for name in ("G1", "G2"):
        fa = Path(tmp) / f"{name}.fa"
        contig = "chr1" if name == "G1" else "scaf_1"
        with open(fa, "w", encoding="utf-8", newline="\n") as fo:
            fo.write(f">{contig}\nACGTACGTACGT\n")
        lines.append(f"{name}\t{fa}\n")
    batch = Path(tmp) / "batch.tsv"
    with open(batch, "w", encoding="utf-8", newline="\n") as fo:
        fo.writelines(["# name\tgenome\n"] + lines)
    return batch


class FakeToolFixture:
    """一次性的临时环境: 假工具 bin/ + 合成输入输出根目录 + 可控的调用记录."""

    def __init__(self, prefix="eve_elem_"):
        self.tmp = Path(tempfile.mkdtemp(prefix=prefix))
        self.bindir = self.tmp / "bin"
        self.mmseqs = write_shim(self.bindir, "mmseqs", MMSEQS_SHIM)
        self.mafft = write_shim(self.bindir, "mafft", MAFFT_SHIM)
        self.rm = write_shim(self.bindir, "RepeatMasker", RM_SHIM)
        self.crash = write_shim(self.bindir, "crash_tool", CRASH_SHIM)
        self.clusters = self.tmp / "fake_clusters.tsv"
        with open(self.clusters, "w", encoding="utf-8", newline="\n") as fo:
            fo.write(FAKE_CLUSTERS)
        self.log = self.tmp / "calls.log"
        self.out = build_screen_output(self.tmp / "screen")
        self.batch = build_batch(self.tmp, self.out)

    def env(self):
        e = dict(os.environ)
        e["PATH"] = str(self.bindir) + os.pathsep + e.get("PATH", "")
        e["PYTHONUTF8"] = "1"
        e["PYTHONIOENCODING"] = "utf-8"
        e["EVE_FAKE_LOG"] = str(self.log)
        e["EVE_FAKE_CLUSTERS"] = str(self.clusters)
        e["EVE_FAKE_DOMAINS"] = FAKE_DOMAINS
        return e

    def run(self, *args, env_extra=None):
        e = self.env()
        if env_extra:
            e.update(env_extra)
        return subprocess.run([PY, str(MOD), "-o", str(self.out), *args],
                              capture_output=True, text=True, encoding="utf-8",
                              errors="replace", env=e)

    def calls(self):
        if not self.log.is_file():
            return []
        return [l for l in self.log.read_text(encoding="utf-8").splitlines() if l]

    def cleanup(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


# ─────────────────────── 1) 纯逻辑 ───────────────────────
class TestLocusParsing(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(ed.locus_coords("chr1:99199-100700"),
                         ("chr1", 99199, 100700))

    def test_contig_with_colon_keeps_left_part_whole(self):
        # 只 rsplit(":", 1): 名字里带 ':' 的 contig 不能被截断
        self.assertEqual(ed.locus_coords("scaf:win1:10-20"),
                         ("scaf:win1", 10, 20))

    def test_malformed(self):
        for bad in ("", "chr1", "chr1:", "chr1:10", "chr1:10-", "chr1:a-b",
                    ":10-20", "chr1:10-20-30"):
            self.assertIsNone(ed.locus_coords(bad), bad)


class TestMergeRegions(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(ed.merge_regions([], 5000), [])

    def test_merge_within_and_beyond(self):
        pts = [(100, 300, "A"), (5000, 5300, "A"), (30000, 30200, "A")]
        got = ed.merge_regions(pts, 5000)
        self.assertEqual(got, [(100, 5300, ["A"]), (30000, 30200, ["A"])])

    def test_boundary_merges_exactly_at_dist(self):
        # s <= ce + er_dist 是"含端点": 间隔正好 5000 要并, 5001 就不并
        self.assertEqual(ed.merge_regions([(0, 100, "A"), (5100, 5200, "A")], 5000),
                         [(0, 5200, ["A"])])
        self.assertEqual(len(ed.merge_regions([(0, 100, "A"), (5101, 5200, "A")],
                                              5000)), 2)

    def test_families_union_sorted_and_dedup(self):
        got = ed.merge_regions([(0, 10, "B"), (5, 20, "A"), (30, 40, "B")], 5000)
        self.assertEqual(got, [(0, 40, ["A", "B"])])

    def test_creates_no_region_for_empty_family_but_keeps_it(self):
        # 源脚本口径: ref_family 为空串时它照样算一个"家族" (保留, 不替上游过滤)
        self.assertEqual(ed.merge_regions([(0, 10, "")], 5000), [(0, 10, [""])])

    def test_unsorted_input_still_merges_in_coordinate_order(self):
        pts = [(30000, 30200, "A"), (100, 300, "A")]
        self.assertEqual([r[:2] for r in ed.merge_regions(pts, 5000)],
                         [(100, 300), (30000, 30200)])


class TestCandidateFilter(unittest.TestCase):
    """viral_supported 行筛选: 按列名取值, 不按列下标 (列错位是静默错数据的经典来源)."""

    HDR = ["locus", "verdict", "ref_family", "ref_sseqid", "ref_bitscore",
           "ref_qcov", "best_viral_sp", "best_viral_bs", "best_plant_sp",
           "best_plant_bs", "rvdb_stitle", "rvdb_bitscore"]

    def _row(self, **kw):
        base = dict(zip(self.HDR, ["chr1:1-9", "viral_supported", "Caulimoviridae",
                                   "REF_A", "120", "80", "REF_A", "120", "HOST",
                                   "5", "", ""]))
        base.update(kw)
        return [base[c] for c in self.HDR]

    def test_verdict_filter(self):
        rows = [self._row(), self._row(verdict="host_like"),
                self._row(verdict="undetermined")]
        got = ed.candidates_from_rows("G1", self.HDR, rows, {"viral_supported"}, 50)
        self.assertEqual(len(got), 1)
        self.assertEqual(got[0][3], "viral_supported")

    def test_evidence_bs_takes_max_of_two_channels(self):
        rows = [self._row(ref_bitscore="30", rvdb_bitscore="210")]
        got = ed.candidates_from_rows("G1", self.HDR, rows, {"viral_supported"}, 50)
        self.assertEqual(got[0][5], 210.0)
        self.assertEqual(len(ed.candidates_from_rows("G1", self.HDR, rows,
                                                     {"viral_supported"}, 300)), 0)

    def test_min_bs_none_keeps_low_evidence(self):
        # regions 阶段的口径: 只看 verdict, 不设证据阈值
        rows = [self._row(ref_bitscore="0.0")]
        self.assertEqual(len(ed.candidates_from_rows("G1", self.HDR, rows,
                                                     {"viral_supported"}, None)), 1)
        self.assertEqual(len(ed.candidates_from_rows("G1", self.HDR, rows,
                                                     {"viral_supported"}, 50)), 0)

    def test_bad_numbers_do_not_crash(self):
        rows = [self._row(ref_bitscore="NA", rvdb_bitscore="")]
        got = ed.candidates_from_rows("G1", self.HDR, rows, {"viral_supported"}, None)
        self.assertEqual(got[0][5], 0.0)

    def test_missing_columns_and_short_rows_tolerated(self):
        # 少列/短行不得错位: 取不到的列按空值处理, 该行照旧按 verdict 判
        hdr = ["locus", "verdict"]
        rows = [["chr1:1-9", "viral_supported"], ["chr1:9-9"]]
        got = ed.candidates_from_rows("G1", hdr, rows, {"viral_supported"}, None)
        self.assertEqual(got, [("G1", "chr1:1-9", "", "viral_supported", "", 0.0)])

    def test_legacy_13_column_summary(self):
        # 旧布局 (首列多一个 genome) 按列名读照样对
        hdr = ["genome"] + self.HDR
        row = ["G1"] + self._row()
        got = ed.candidates_from_rows("G1", hdr, [row], {"viral_supported"}, 50)
        self.assertEqual(got[0][1], "chr1:1-9")
        self.assertEqual(got[0][2], "Caulimoviridae")

    def test_empty_locus_skipped(self):
        rows = [self._row(locus="")]
        self.assertEqual(ed.candidates_from_rows("G1", self.HDR, rows,
                                                 {"viral_supported"}, None), [])


class TestConsensusInput(unittest.TestCase):
    def setUp(self):
        self.seqs = {"m1": "A" * 10, "m2": "A" * 20, "m3": "A" * 30}

    def test_rep_picks_longest_and_tie_takes_smallest_name(self):
        self.assertEqual(ed.pick_representative(["m2", "m3", "m1"], self.seqs), "m3")
        tie = {"b": "AA", "a": "CC"}
        self.assertEqual(ed.pick_representative(["b", "a"], tie), "a")

    def test_rep_when_no_sequence_available(self):
        self.assertIsNone(ed.pick_representative(["zz"], self.seqs))
        mode, payload = ed.consensus_input(["zz"], self.seqs, False)
        self.assertEqual((mode, payload), ("rep", []))

    def test_mafft_mode_by_member_count(self):
        mode, payload = ed.consensus_input(["m1", "m2"], self.seqs, True)
        self.assertEqual(mode, "mafft")
        self.assertEqual(payload, {"m1": "A" * 10, "m2": "A" * 20})
        # 单成员 / 超过 50 成员: 退回最长成员 (与源脚本 2..50 的判据一致)
        self.assertEqual(ed.consensus_input(["m1"], self.seqs, True)[0], "rep")
        many = {f"m{i}": "A" * i for i in range(1, 60)}
        self.assertEqual(ed.consensus_input(list(many), many, True)[0], "rep")

    def test_mafft_falls_back_when_sequences_vanished(self):
        # 成员数够但序列不在 all_candidates.fa 里: 源脚本会拿空集合去跑 mafft
        # (mafft 报错), 这里退回最长成员
        mode, payload = ed.consensus_input(["m1", "gone1", "gone2"], self.seqs, True)
        self.assertEqual(mode, "rep")
        self.assertEqual(payload, [("m1", "A" * 10)])

    def test_mafft_needs_two_available_sequences(self):
        seqs = {"keep": "AAAA", "gone": ""}
        mode, _p = ed.consensus_input(["keep", "gone"], seqs, True)
        self.assertEqual(mode, "rep")


class TestMsaConsensus(unittest.TestCase):
    def test_frequency_rule_and_uppercase(self):
        # 第 1 列 3 个 A vs 1 个 C: 3/4 = 0.75 >= 0.4 取 A; 小写也统一成大写
        recs = [("a", "A"), ("b", "a"), ("c", "A"), ("d", "C")]
        self.assertEqual(ed.msa_consensus(recs), "A")

    def test_below_frequency_threshold_becomes_n(self):
        # 第 2 列 5 条各不相同: 众数 1/5 = 0.2 < 0.4 -> n
        recs = [("1", "AA"), ("2", "AC"), ("3", "AG"), ("4", "AT"), ("5", "AN")]
        self.assertEqual(ed.msa_consensus(recs), "An")
        # 0.4 是含端点的阈值: 2/5 = 0.4 恰好取众数
        boundary = [("1", "AT"), ("2", "AT"), ("3", "AC"), ("4", "AG"), ("5", "AA")]
        self.assertEqual(ed.msa_consensus(boundary), "AT")

    def test_gap_majority_becomes_n(self):
        # 空位占多数 (3/5) 时即使空位自己过了频率阈值也写 n (不确定位)
        recs = [("a", "AT"), ("b", "AG"), ("c", "A-"), ("d", "C-"), ("e", "G-")]
        self.assertEqual(ed.msa_consensus(recs), "An")

    def test_column_with_fewer_than_two_non_gap_dropped(self):
        # 第 2 列只有 1 条非空位: 整列丢弃 (不产出也不占位)
        recs = [("a", "A-"), ("b", "A-"), ("c", "A-"), ("d", "AT")]
        self.assertEqual(ed.msa_consensus(recs), "A")

    def test_empty(self):
        self.assertEqual(ed.msa_consensus([]), "")


class TestAggregateClusters(unittest.TestCase):
    def _meta(self):
        return {
            "g1__L1": ["g1", "L1", "Caulimoviridae", "viral_supported", "REF_A", "120"],
            "g2__L2": ["g2", "L2", "Caulimoviridae", "viral_supported", "REF_A", "90"],
            "g3__L3": ["g3", "L3", "Geminiviridae", "viral_supported", "REF_B", "300"],
        }

    def test_rows_and_ids(self):
        pairs = [("g1__L1", "g1__L1"), ("g1__L1", "g2__L2"), ("g3__L3", "g3__L3")]
        rows, members = ed.aggregate_clusters(self._meta(), pairs)
        self.assertEqual([r[0] for r in rows], ["C000000", "C000001"])
        self.assertEqual(rows[0][2], 2)                     # n_members
        self.assertEqual(rows[0][3], 2)                     # n_genomes
        self.assertEqual(rows[0][4], "Caulimoviridae")      # 家族投票
        self.assertEqual(rows[0][5], "")                    # 无次家族 -> 非嵌合体
        self.assertEqual(rows[0][6], "g1__L1")              # 最佳成员按 evidence_bs
        self.assertEqual(rows[0][7], 120.0)
        self.assertEqual(members["C000000"], ["g1__L1", "g2__L2"])

    def test_members_sorted_and_tie_order_stable(self):
        # 两个等大簇: 并列时按 rep 名排序 -> 元件 ID 不随 mmseqs 行序漂
        meta = self._meta()
        meta["g9__L9"] = ["g9", "L9", "Caulimoviridae", "viral_supported", "R", "10"]
        pairs_a = [("g3__L3", "g3__L3"), ("g3__L3", "g9__L9")]
        pairs_b = [("g3__L3", "g9__L9"), ("g3__L3", "g3__L3")]
        self.assertEqual(ed.aggregate_clusters(meta, pairs_a)[0],
                         ed.aggregate_clusters(meta, pairs_b)[0])
        _rows, mem = ed.aggregate_clusters(meta, pairs_b)
        self.assertEqual(mem["C000000"], ["g3__L3", "g9__L9"])

    def test_chimera_flag_needs_minor_family_at_threshold(self):
        def build(minor):
            meta = {f"g{i}__L": [f"g{i}", "L", "Caulimoviridae", "viral_supported",
                                 "R", "10"] for i in range(5 - minor)}
            for i in range(minor):
                meta[f"gm{i}__L"] = [f"gm{i}", "L", "Geminiviridae",
                                     "viral_supported", "R", "10"]
            pairs = [("g0__L", m) for m in sorted(meta)]
            return ed.aggregate_clusters(meta, pairs)[0]
        # 5 个成员里 1 个次家族: 1 < max(2, 5//5) = 2 -> 不标嵌合体
        self.assertEqual(build(1)[0][5], "")
        # 5 个成员里 2 个次家族: 2 >= 2 -> 标 (且不改主家族投票)
        rows = build(2)
        self.assertEqual(rows[0][5], "possible_chimera")
        self.assertEqual(rows[0][4], "Caulimoviridae")

    def test_unknown_member_does_not_crash(self):
        rows, _mem = ed.aggregate_clusters(self._meta(), [("ghost", "ghost")])
        self.assertEqual(rows[0][4], "Other_viral")
        self.assertEqual(rows[0][3], 0)

    def test_empty_input(self):
        self.assertEqual(ed.aggregate_clusters({}, []), ([], {}))


class TestBatchAndRmParse(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="eve_elem_unit_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write(self, name, text):
        p = self.tmp / name
        with open(p, "w", encoding="utf-8", newline="\n") as fo:
            fo.write(text)
        return p

    def test_parse_batch_filters_and_cleans(self):
        p = self._write("b.tsv", "# 注释\n"
                                  "G1\t~/g1.fa\n"
                                  "with/slash\t/tmp/g2.fa\n"
                                  "G1\n"                       # 只有一列: 跳过
                                  "\n")
        got = ed.parse_batch(p)
        self.assertEqual(len(got), 2)
        self.assertEqual(got[0][0], "G1")
        self.assertEqual(got[1][0], "with_slash")        # clean_name 挡住路径逃逸
        self.assertNotIn("~", got[0][1])                 # expanduser 生效

    def test_pending_rm_jobs_skips_existing(self):
        rm = self.tmp / "rm"
        (rm / "G1").mkdir(parents=True)
        (rm / "G1" / "G1.out").write_text("done", encoding="utf-8", newline="\n")
        got = ed.pending_rm_jobs([("G1", "/g1.fa"), ("G2", "/g2.fa")], rm)
        self.assertEqual(got, [("G2", "/g2.fa")])

    def test_build_lib_fam_map_keeps_underscore_in_family(self):
        m = ed.build_lib_fam_map(["C000000_Caulimoviridae",
                                  "C000001_Caulimoviridae_endo_",
                                  "C000002"])
        self.assertEqual(m["C000000"], "Caulimoviridae")
        self.assertEqual(m["C000001"], "Caulimoviridae_endo_")
        self.assertEqual(m["C000002"], "C000002")

    def test_parse_rm_out_lines_on_real_repeatmasker_format(self):
        """真 .out 的列序回归: 0 SW score, 1-3 perc, 4 query 名, 5/6 query begin/end,
        7 left, 8 链, 9 repeat, 10 class/family, 11-13 repeat 坐标, 14 ID.

        源脚本 (eve_kingdom_post.py:387-391) 用 c[4]/c[5] 当 begin/end —— c[4] 是
        query 序列名, 对真实数据行永远不是数字, 于是**每一行命中都被当表头丢掉**,
        family_quant.tsv 只剩表头而退出码是 0。这里用一段真 .out 头+数据行锁住。
        """
        text = (
            "   SW   perc perc perc  query     position in query              "
            "matching                    repeat        position in repeat\n"
            "score   div. del. ins.  sequence  begin    end          (left)   "
            "repeat                      class/family begin  end    (left)   ID\n"
            "\n"
            " 2334    3.8  0.5  0.0  chr4         50350    50657 (18534399) + "
            "C000000_Caulimoviridae   LTR/Gypsy     1745   2167 (1383)    1\n"
            " 10602   3.9  2.6  0.0  chr4         51025    52324 (18532732) C "
            "C000000_Caulimoviridae   LTR/Gypsy    (324)   3226   1893    2\n"
            "There were no repetitive sequences detected in the query sequence\n")
        got = list(ed.parse_rm_out_lines(text.splitlines()))
        self.assertEqual(got, [(50350, 50657, "C000000_Caulimoviridae"),
                               (51025, 52324, "C000000_Caulimoviridae")])

    def test_parse_rm_out_lines_strips_hash_and_ignores_short_rows(self):
        rows = ["  1 2 3 4 q 10 20 (0) + name#class/1  c  1 2 (0) 1",
                "  1 2 3 q 10 20 30",           # 7 列 (< 11): 丢掉
                "  1 2 3 4 q 10 20 (0) + nm c"]  # 11 列: 够解析 (尾部 ID 列可缺)
        got = list(ed.parse_rm_out_lines(rows))
        self.assertEqual(got, [(10, 20, "name"), (10, 20, "nm")])

    def test_padded_left_column_is_skipped_not_misattributed(self):
        # 万一某版本把 (left) 写成 '(     0)': split() 会多切一列, 硬编 c[9] 就取到
        # strand 本身, 会把命计算到一个叫 '+' 的家族头上 —— 宁可整行跳过
        row = "  1 2 3 4 q 10 20 (     0) + nm c 1 2 (0) 1"
        self.assertEqual(list(ed.parse_rm_out_lines([row])), [])

    def test_accumulate_quant_counts_copies_and_bp(self):
        lib = {"C000000": "Caulimoviridae", "C000001": "Geminiviridae"}
        hits = [(50350, 50657, "C000000_Caulimoviridae"),
                (51025, 52324, "C000000_Caulimoviridae"),
                (100, 150, "C000001_Geminiviridae")]
        q = ed.accumulate_quant(hits, lib, "G1")
        self.assertEqual(q[("Caulimoviridae", "G1")], [2, 308 + 1300])
        self.assertEqual(q[("Geminiviridae", "G1")], [1, 51])

    def test_accumulate_quant_unknown_id_falls_back_to_repeat_name(self):
        # 文库头不认识时至少留下名字本身, 不要静默归到某个家族名下
        q = ed.accumulate_quant([(1, 10, "ODD_name")], {}, "G1")
        self.assertEqual(list(q), [("ODD_name", "G1")])


class TestBookkeeping(unittest.TestCase):
    def test_stage_outputs_cover_every_stage(self):
        self.assertEqual(set(ed.STAGE_OUTPUTS), set(ed.ALL_STAGES))

    def test_parse_stages(self):
        self.assertEqual(ed.parse_stages("collect,cluster"),
                         {"collect", "cluster"})
        self.assertEqual(ed.parse_stages("1,2"), {"collect", "cluster"})
        self.assertEqual(ed.parse_stages("all"), set(ed.ALL_STAGES))
        self.assertEqual(ed.parse_stages(""), set(ed.ALL_STAGES))
        self.assertEqual(ed.parse_stages(",".join(ed.DEFAULT_STAGES)),
                         {"collect", "cluster", "regions", "consensus"})
        with self.assertRaises(SystemExit):
            ed.parse_stages("nonsense")

    def test_default_stages_need_no_extra_db(self):
        # 默认阶段必须是不依赖额外参考库的那四个 (文档与 CLI 帮助都这么写)
        self.assertEqual(set(ed.DEFAULT_STAGES),
                         {"collect", "cluster", "regions", "consensus"})

    def test_safe_tag_vs_clean_name(self):
        self.assertEqual(ed.safe_tag("Caulimoviridae(endo)"),
                         "Caulimoviridae_endo_")
        self.assertEqual(ed.safe_tag("A:B/C"), "A_B_C")
        # 差别只在空名/点开头: clean_name 要拿它拼目录名, 以 asm_ 兜底;
        # safe_tag 只做字符替换 (. 是安全字符, 不动)
        self.assertEqual(ed.safe_tag(""), "")
        self.assertEqual(ed.clean_name(""), "asm_")
        self.assertEqual(ed.safe_tag(".hidden"), ".hidden")
        self.assertEqual(ed.clean_name(".hidden"), "asm_.hidden")

    def test_writers_emit_lf_only(self):
        tmp = Path(tempfile.mkdtemp(prefix="eve_elem_io_"))
        try:
            t = tmp / "t.tsv"
            ed.write_table(t, ed.CAND_META_HDR, [("a", "b")])
            f = tmp / "t.fa"
            ed.write_fa(f, [("s1", "ACGTACGT" * 20)])
            for p in (t, f):
                self.assertNotIn(b"\r", p.read_bytes(), p)
            self.assertEqual(t.read_bytes().splitlines()[0].decode(), ed.CAND_META_HDR)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


# ─────────────────────── 2) 端到端 (假工具) ───────────────────────
class TestElementLayerEndToEnd(unittest.TestCase):
    """假 mmseqs/mafft/RepeatMasker 跑完 collect→cluster→regions→consensus→rm."""

    @classmethod
    def setUpClass(cls):
        cls.fx = FakeToolFixture()
        cls.res = cls.fx.run("--stages", "collect,cluster,regions,consensus",
                             "--rm", "-B", str(cls.fx.batch), "-t", "2", "-J", "2")
        cls.elem = cls.fx.out / ed.ELEM_DIR
        cls.calls = cls.fx.calls()

    @classmethod
    def tearDownClass(cls):
        cls.fx.cleanup()

    def read(self, name):
        return (self.elem / name).read_text(encoding="utf-8").splitlines()

    def test_run_succeeds_and_prints_done_per_stage(self):
        self.assertEqual(self.res.returncode, 0, self.res.stderr)
        for stage in ("collect", "cluster", "regions", "consensus", "rm"):
            self.assertIn(f"[done] {stage}:", self.res.stderr, stage)
        self.assertIn("[done] all:", self.res.stderr)

    def test_outputs_exist(self):
        for stage, files in ed.STAGE_OUTPUTS.items():
            if stage in ("cdd", "phylo"):
                continue
            for f in files:
                self.assertTrue((self.elem / f).is_file(), f"{stage}/{f}")

    def test_header_constants_match_written_tables(self):
        """表头常量必须与写盘的第一行逐列一致 (常量在别处再抄一份就会悄悄漂)."""
        for fname, hdr in (("candidates_meta.tsv", ed.CAND_META_HDR),
                           ("element_table.tsv", ed.ELEMENT_TABLE_HDR),
                           ("element_regions.tsv", ed.ELEMENT_REGIONS_HDR),
                           ("family_quant.tsv", ed.FAMILY_QUANT_HDR)):
            self.assertEqual(self.read(fname)[0].split("\t"), hdr.split("\t"),
                             fname)

    def test_collect_counts_and_meta(self):
        # G1: 3 条 viral_supported 过阈值 (host_like 与 bs=10 的两条被挡掉), G2: 2 条
        meta = self.read("candidates_meta.tsv")[1:]
        self.assertEqual(len(meta), 5, meta)
        self.assertEqual([r.split("\t")[0] for r in meta], ["G1"] * 3 + ["G2"] * 2)
        first = meta[0].split("\t")
        self.assertEqual(first[1], "chr1:100-300")
        self.assertEqual(first[2], "Caulimoviridae")
        self.assertEqual(first[5], "200.0")      # max(ref 120, rvdb 200)
        # 序列收集齐全 (ID 带基因组前缀, 防跨基因组撞名)
        fa = (self.elem / "all_candidates.fa").read_text(encoding="utf-8")
        self.assertIn(">G1__chr1:100-300\n", fa)
        self.assertIn(">G2__scaf_1:60000-60200\n", fa)

    def test_regions_apply_no_bitscore_threshold(self):
        rows = [r.split("\t") for r in self.read("element_regions.tsv")[1:]]
        # G1 的 chr1 三个位点里, 100-300 与 5000-5300 差 4700 <= 5000 并成一个区
        g1 = [r for r in rows if r[0] == "G1"]
        self.assertEqual([(r[2], r[3], r[4]) for r in g1],
                         [("chr1", "100", "5300"), ("chr1", "30000", "30200"),
                          ("chr2", "100", "200")])
        # bs=10 的那条不在 collect 里 (证据阈值) 却在元件区表里 (只看 verdict):
        # 两个阶段的取行口径确实不同, 不是同一个过滤器抄了两遍
        collect_loci = [r.split("\t")[1] for r in
                        self.read("candidates_meta.tsv")[1:]]
        self.assertNotIn("chr1:30000-30200", collect_loci)
        self.assertNotIn(">G1__chr1:30000-30200",
                         self.read("all_candidates.fa"))

    def test_region_ids_are_unique_per_genome(self):
        ids = [r.split("\t")[1] for r in self.read("element_regions.tsv")[1:]]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(ids[0], "G1__ER00001")

    def test_element_table(self):
        rows = [r.split("\t") for r in self.read("element_table.tsv")[1:]]
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0][:5], ["C000000", "G1__chr1:100-300", "3", "2",
                                       "Caulimoviridae"])
        self.assertEqual(rows[0][5], "")                      # 无次家族
        self.assertEqual(rows[0][6], "G1__chr1:100-300")      # 最佳 evidence_bs
        self.assertEqual(rows[0][7], "200.0")
        self.assertEqual(rows[1][0], "C000001")
        self.assertEqual(rows[1][4], "Geminiviridae")

    def test_member_files_written_sorted(self):
        self.assertEqual(self.read("_members_C000000.txt"),
                         ["G1__chr1:100-300", "G1__chr1:5000-5300",
                          "G2__scaf_1:1000-1200"])

    def test_consensus_library_uses_longest_member(self):
        lines = self.read("eve_consensus_library.fa")
        self.assertEqual(lines[0], ">C000000_Caulimoviridae")
        self.assertEqual(lines[1], SEQS["chr1:100-300"])
        self.assertEqual(lines[2], ">C000001_Geminiviridae")
        self.assertEqual(lines[3], SEQS["scaf_1:60000-60200"])
        self.assertIn("(longest-member)", self.res.stderr)

    def test_family_quant_counts_hits_per_family_and_genome(self):
        """RepeatMasker .out 的列序回归 (源脚本在这里取错一列 -> 表恒空)."""
        rows = [r.split("\t") for r in self.read("family_quant.tsv")[1:]]
        self.assertEqual(sorted((r[0], r[1], r[2], r[3]) for r in rows),
                         [("Caulimoviridae", "G1", "1", "200"),
                          ("Caulimoviridae", "G2", "1", "200"),
                          ("Geminiviridae", "G1", "1", "200"),
                          ("Geminiviridae", "G2", "1", "200")])

    def test_tool_arguments_are_not_retuned(self):
        """算法参数是继承来的, 不许悄悄改: 从调用记录里核对."""
        linclust = [c for c in self.calls if c.startswith("mmseqs easy-linclust")]
        self.assertTrue(linclust, self.calls)
        for want in ("--min-seq-id 0.8", "-c 0.8", "--cov-mode 1"):
            self.assertIn(want, linclust[0])
        rm = [c for c in self.calls if c.startswith("RepeatMasker")]
        self.assertTrue(rm)
        for want in ("-no_is", "-norna", "-lib", "-pa 2"):
            self.assertIn(want, rm[0])

    def test_second_run_skips_everything_without_calling_tools(self):
        n_before = len(self.fx.calls())
        r = self.fx.run("--stages", "collect,cluster,regions,consensus",
                        "--rm", "-B", str(self.fx.batch))
        self.assertEqual(r.returncode, 0, r.stderr)
        for stage in ("collect", "cluster", "regions", "consensus", "rm"):
            self.assertIn(f"[done] {stage}: 产物已存在, 跳过", r.stderr, stage)
        self.assertEqual(len(self.fx.calls()), n_before, "断点跳过时不得再调外部工具")

    def test_force_reruns(self):
        n_before = len(self.fx.calls())
        r = self.fx.run("--stages", "cluster", "--force")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("[done] cluster:", r.stderr)
        self.assertNotIn("跳过", r.stderr)
        self.assertGreater(len(self.fx.calls()), n_before)

    def test_mafft_consensus_mode(self):
        r = self.fx.run("--stages", "consensus", "--force", "--mafft-consensus",
                        "--mmseqs", str(self.fx.mmseqs))
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("(mafft)", r.stderr)
        lines = self.read("eve_consensus_library.fa")
        self.assertEqual(lines[0], ">C000000_Caulimoviridae")
        # 三条成员序列都以公共核心开头 -> 共识串必须以它开头
        self.assertTrue(lines[1].startswith(CORE), lines[1])
        self.assertTrue(any(c.startswith("mafft ") for c in self.fx.calls()))

    def test_cdd_stage_uses_translated_search_and_field_spec(self):
        db = self.fx.tmp / "cdd_db"
        db.mkdir(exist_ok=True)
        r = self.fx.run("--stages", "cdd", "--cdd-db", str(db))
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("[done] cdd: 1 条结构域命中", r.stderr)
        self.assertEqual((self.elem / "consensus_domains.tsv").read_text(
            encoding="utf-8"), FAKE_DOMAINS)
        search = [c for c in self.fx.calls() if c.startswith("mmseqs search")]
        self.assertIn("--search-type 2", search[-1])
        self.assertIn("-s 4.0", search[-1])
        self.assertIn("--max-seqs 300", search[-1])
        conv = [c for c in self.fx.calls() if c.startswith("mmseqs convertalis")]
        self.assertIn(ed.CDD_CONVERT_FIELDS, conv[-1])

    def test_locus_key_join_tolerates_rewritten_fasta_headers(self):
        """抽序列工具把 ':' 改写成 '_' 时, 位点序列不许被静默丢掉."""
        fx = FakeToolFixture(prefix="eve_elem_uf_")
        try:
            build_screen_output(fx.out, underscore_headers=True)
            r = fx.run("--stages", "collect", "--force")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("[done] collect: 5 个位点", r.stderr)
        finally:
            fx.cleanup()


class TestDeterminism(unittest.TestCase):
    """同一份输入跑两遍, 所有表逐字节一致 (本仓库的硬要求)."""

    TABLES = ("candidates_meta.tsv", "cluster_members.tsv", "element_table.tsv",
              "element_regions.tsv", "eve_consensus_library.fa", "family_quant.tsv",
              "all_candidates.fa", "rm/G1/G1.out")

    def test_byte_identical_on_rerun(self):
        fx = FakeToolFixture(prefix="eve_elem_det_")
        try:
            args = ("--stages", "collect,cluster,regions,consensus", "--rm",
                    "-B", str(fx.batch), "-t", "2", "-J", "2")
            r1 = fx.run(*args)
            self.assertEqual(r1.returncode, 0, r1.stderr)
            elem = fx.out / ed.ELEM_DIR
            snap = {t: (elem / t).read_bytes() for t in self.TABLES}
            self.assertGreater(len(snap["candidates_meta.tsv"]), 0)
            r2 = fx.run(*args, "--force")
            self.assertEqual(r2.returncode, 0, r2.stderr)
            for t in self.TABLES:
                self.assertEqual((elem / t).read_bytes(), snap[t], t)
        finally:
            fx.cleanup()


# ─────────────────────── 3) 失败必须可见 ───────────────────────
class TestFailuresAreVisible(unittest.TestCase):
    def setUp(self):
        self.fx = FakeToolFixture(prefix="eve_elem_fail_")

    def tearDown(self):
        self.fx.cleanup()

    def test_crashing_tool_fails_loudly_and_writes_no_table(self):
        r = self.fx.run("--stages", "collect,cluster",
                        "--mmseqs", str(self.fx.crash))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("[done] collect:", r.stderr)           # 前一个阶段是好的
        self.assertIn("[fail] cluster:", r.stderr)
        self.assertIn("cluster", r.stderr)
        elem = self.fx.out / ed.ELEM_DIR
        self.assertTrue((elem / "candidates_meta.tsv").is_file())
        self.assertFalse((elem / "element_table.tsv").exists(),
                         "工具崩了不许留半成品表冒充成功")

    def test_mmseqs_missing_cluster_output_is_fatal(self):
        # 退出码 0 但没产出 _clu_cluster.tsv: 必须指名报错, 不是 FileNotFoundError
        fake = write_shim(self.fx.tmp / "bin2", "mmseqs",
                          "import sys\nsys.exit(0)\n")
        r = self.fx.run("--stages", "collect,cluster", "--mmseqs", str(fake))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("[fail] cluster:", r.stderr)
        self.assertIn("_clu_cluster.tsv", r.stderr)

    def test_rm_without_batch_is_fatal_not_silent(self):
        # 源脚本: rm 在 stages 里但没 -g 就静默什么都不做, 退出码 0
        r = self.fx.run("--stages", "rm")
        self.assertEqual(r.returncode, 1)
        self.assertIn("--batch", r.stderr)
        self.assertNotIn("[done] rm:", r.stderr)

    def test_cdd_without_db_is_fatal_not_silent(self):
        r = self.fx.run("--stages", "cdd")
        self.assertEqual(r.returncode, 1)
        self.assertIn("--cdd-db", r.stderr)
        self.assertNotIn("[done] cdd:", r.stderr)

    def test_run_phylo_without_rt_ref_is_fatal_not_silent(self):
        r = self.fx.run("--stages", "phylo", "--run-phylo")
        self.assertEqual(r.returncode, 1)
        self.assertIn("--rt-ref", r.stderr)

    def test_stages_all_names_missing_pieces(self):
        # all 会连带 rm/cdd/phylo: 缺的东西必须开跑前一次性报出来
        # (假 mmseqs / RepeatMasker 已经在 PATH 上, 所以这两样不算缺)
        r = self.fx.run("--stages", "all")
        self.assertEqual(r.returncode, 1)
        for want in ("--cdd-db", "--batch"):
            self.assertIn(want, r.stderr)
        self.assertNotIn("mmseqs=", r.stderr)

    def test_outdir_without_summary_is_fatal(self):
        empty = Path(tempfile.mkdtemp(prefix="eve_elem_empty_"))
        try:
            r = subprocess.run([PY, str(MOD), "-o", str(empty)],
                               capture_output=True, text=True, encoding="utf-8",
                               errors="replace", env=self.fx.env())
            self.assertEqual(r.returncode, 1)
            self.assertIn("[fail] collect:", r.stderr)
            self.assertIn("04_Summary", r.stderr)
        finally:
            shutil.rmtree(empty, ignore_errors=True)

    def test_no_candidates_is_fatal(self):
        # 把 verdict 换成人人都不是的值 -> 收集为空必须报错, 不写空表当成功
        r = self.fx.run("--stages", "collect", "--verdicts", "nonexistent_verdict")
        self.assertEqual(r.returncode, 1)
        self.assertIn("[fail] collect:", r.stderr)
        self.assertIn("没有任何 viral_supported 位点", r.stderr)


class TestCli(unittest.TestCase):
    def test_help(self):
        r = subprocess.run([PY, str(MOD), "--help"], capture_output=True,
                           text=True, encoding="utf-8", errors="replace")
        self.assertEqual(r.returncode, 0, r.stderr)
        for want in ("--stages", "--outdir", "--cdd-db", "--rt-ref", "--rm-bin",
                     "--force", "collect", "phylo"):
            self.assertIn(want, r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
