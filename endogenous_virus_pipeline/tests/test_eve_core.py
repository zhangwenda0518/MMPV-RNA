#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_eve_core.py — eve_scan_core 纯逻辑单元测试 (无需 diamond/seqkit/samtools)

运行: python -m unittest discover -s endogenous_virus_pipeline/tests -v
  或: python endogenous_virus_pipeline/tests/test_eve_core.py

用合成的 diamond blastx 输出锁定行为: 坐标还原 / 位点合并 / 最佳命中归属 /
双侧判定 / 单基因组汇总 / 跨基因组 merge / 名字清洗 / seqkit 命令构造 /
工具体检 / BED→samtools 区域文件换算 / 命令超时 / 批次重名 / 压缩输入推断.
"""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eve_scan_core import (EveConfig, assign_best_hits, bed_to_regions,
                           check_extracted, check_tools, clean_name,
                           extract_cmd, family_of, load_id2div, locus_key,
                           merge_all, merge_loci, restore_coordinates, run,
                           samtools_too_old, sliding_cmd, summarize_genome,
                           verdict_loci)
from eve_genome_scan import parse_stages, resolve_genome_fa


def write(path, text):
    Path(path).write_text(text, encoding="utf-8")


class _Log:
    """run()/build_jobs 只需要 log.info / log.error 的哑日志."""
    def info(self, *a, **k):
        pass

    def error(self, *a, **k):
        pass


class TestCleanName(unittest.TestCase):
    def test_sanitize(self):
        self.assertEqual(clean_name("Solanum lycopersicum (v6)"), "Solanum_lycopersicum__v6_")
        # 以 '.' 开头会被加 asm_ 前缀 (防隐藏/相对路径目录名)
        self.assertEqual(clean_name("../../etc/passwd"), "asm_.._.._etc_passwd")
        self.assertEqual(clean_name(""), "asm_")
        self.assertEqual(clean_name(".hidden"), "asm_.hidden")
        # 关键: 不允许残留 '..' 路径逃逸风险之外, 仅保留安全字符
        self.assertNotIn("/", clean_name("a/b"))


class TestRestoreCoordinates(unittest.TestCase):
    """滑窗 query 坐标 -> 基因组坐标 (seqkit sliding 命名 <contig>_sliding:<start>-<end>)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.raw = Path(self.tmp.name) / "s1_raw.tsv"
        self.bed = Path(self.tmp.name) / "s1_hits.bed"

    def tearDown(self):
        self.tmp.cleanup()

    def test_two_chunks(self):
        # OUTFMT1: qseqid qlen sseqid slen stitle evalue bitscore qcovhsp pident length qstart qend
        write(self.raw,
              "chr1_sliding:0-50000\t50000\tNC_001\t900\tCauliflower mosaic virus MP\t"
              "1e-30\t120\t80\t60.0\t400\t100\t500\n"
              "chr1_sliding:50000-100000\t50000\tNC_002\t300\tTobacco mosaic virus CP\t"
              "1e-20\t90\t70\t55.0\t300\t1\t300\n"
              "chr2_sliding:0-50000\t50000\tNC_003\t300\tBadnavirus RT\t"
              "1e-15\t80\t60\t50.0\t300\t50\t250\n")
        n = restore_coordinates(self.raw, self.bed)
        self.assertEqual(n, 3)
        lines = self.bed.read_text().splitlines()
        # 块起点 0 + qstart-1 = 98 (0-based bed start = gs-1)
        self.assertTrue(lines[0].startswith("chr1\t98\t499\tchr1_sliding"))
        # 块起点 50000: 50000+1-1=50000 → bed start 49999
        self.assertTrue(lines[1].startswith("chr1\t49999\t50299\tchr1_sliding"))
        # 块起点 0, qstart=50: gs=49 → bed start 48
        self.assertTrue(lines[2].startswith("chr2\t48\t249\tchr2_sliding"))

    def test_empty(self):
        write(self.raw, "")
        self.assertEqual(restore_coordinates(self.raw, self.bed), 0)
        self.assertEqual(self.bed.read_text(), "")


class TestMergeLoci(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.hits = Path(self.tmp.name) / "hits.bed"
        self.loci = Path(self.tmp.name) / "loci.bed"

    def tearDown(self):
        self.tmp.cleanup()

    def test_merge_within_distance(self):
        # 98-499 与 600-700: 600 <= 499+300 → 合并为 98-700
        # 49999-50299 与 98-700 间隔 >300 → 独立位点
        write(self.hits,
              "chr1\t98\t499\tx\nchr1\t600\t700\tx\nchr1\t49999\t50299\tx\n"
              "chr1\t50500\t50700\tx\n")   # 50500 <= 50299+300 → 并入上位点
        n = merge_loci(self.hits, self.loci, merge_d=300)
        self.assertEqual(n, 2)
        lines = self.loci.read_text().splitlines()
        self.assertEqual(lines[0], "chr1\t98\t700\tchr1:99-700")
        self.assertEqual(lines[1], "chr1\t49999\t50700\tchr1:50000-50700")

    def test_nearby_not_merged(self):
        write(self.hits, "chr1\t0\t100\tx\nchr1\t401\t500\tx\n")  # 401 > 100+300
        self.assertEqual(merge_loci(self.hits, self.loci, merge_d=300), 2)

    def test_empty(self):
        write(self.hits, "")
        self.assertEqual(merge_loci(self.hits, self.loci), 0)
        self.assertEqual(self.loci.read_text(), "")


class TestAssignBestHits(unittest.TestCase):
    """每位点取 bitscore 最大的 ref 命中 (含 bed+raw 混合列解析)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.loci = Path(self.tmp.name) / "loci.bed"
        self.hits = Path(self.tmp.name) / "hits.bed"
        self.best = Path(self.tmp.name) / "s1_best.tsv"
        write(self.loci, "chr1\t98\t700\tchr1:99-700\n")
        # hits.bed 第4列起为原始 raw 行 (15 列):
        # qseqid qlen sseqid slen stitle evalue bitscore qcovhsp pident length qstart qend
        # c[5]=sseqid c[7]=stitle c[8]=evalue c[9]=bitscore c[10]=qcovhsp
        write(self.hits,
              "chr1\t98\t499\t"
              "chunk\t50000\tNC_001\t900\tCauliflower mosaic virus movement protein\t"
              "1e-30\t80\t70\t60.0\t400\t100\t500\n"
              "chr1\t300\t700\t"
              "chunk\t50000\tNC_002\t900\tBadnavirus RT\t"
              "1e-40\t120\t90\t65.0\t400\t350\t750\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_best_bitscore(self):
        n = assign_best_hits(self.hits, self.loci, self.best)
        self.assertEqual(n, 1)
        lines = self.best.read_text().splitlines()
        self.assertEqual(lines[0],
                         "locus\tref_sseqid\tfamily\tref_evalue\tref_bitscore\tref_qcov")
        f = lines[1].split("\t")
        self.assertEqual(f[0], "chr1:99-700")
        self.assertEqual(f[1], "NC_002")          # bitscore 120 胜出
        self.assertEqual(f[2], "Caulimoviridae")  # badna -> Caulimoviridae
        self.assertEqual(f[3], "1e-40")
        self.assertEqual(f[4], "120")
        self.assertEqual(f[5], "90")


class TestVerdict(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.raw = Path(self.tmp.name) / "s2_raw.tsv"
        self.out = Path(self.tmp.name) / "s2_verdict.tsv"
        self.id2div = Path(self.tmp.name) / "id2div.tsv"
        write(self.id2div, "SP_VIRAL_A\tviral\nSP_PLANT_A\tplant\nSP_VIRAL_B\tviral\nSP_PLANT_B\tplant\n")
        # OUTFMT2: qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle
        write(self.raw,
              # A: viral 120 >= plant 40 → viral_supported
              "L1\tSP_VIRAL_A\t40\t100\t1\t100\t1\t100\t1e-20\t120\tv\n"
              "L1\tSP_PLANT_A\t30\t80\t1\t80\t1\t80\t1e-5\t40\tp\n"
              # B: plant 100 > viral 30 且过 host_bs → host_like
              "L2\tSP_VIRAL_B\t20\t60\t1\t60\t1\t60\t1e-3\t30\tv\n"
              "L2\tSP_PLANT_B\t60\t120\t1\t120\t1\t120\t1e-30\t100\tp\n"
              # C: 两侧都不过阈值 → undetermined
              "L3\tSP_VIRAL_A\t20\t40\t1\t40\t1\t40\t1e-2\t20\tv\n"
              "L3\tSP_PLANT_A\t20\t40\t1\t40\t1\t40\t1e-2\t20\tp\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_two_sided(self):
        div = load_id2div(self.id2div)
        nv, nh, nu = verdict_loci(self.raw, div, self.out, 50, 50)
        self.assertEqual((nv, nh, nu), (1, 1, 1))
        rows = {l.split("\t")[0]: l.split("\t") for l in
                self.out.read_text().splitlines()[1:]}
        self.assertEqual(rows["L1"][1], "viral_supported")
        self.assertEqual(rows["L1"][3], "120")
        self.assertEqual(rows["L2"][1], "host_like")
        self.assertEqual(rows["L3"][1], "undetermined")


class TestFamilyOf(unittest.TestCase):
    def test_keywords(self):
        # 属/科级词根命中
        self.assertEqual(family_of("Cauliflower mosaic virus isolate X"), "Caulimoviridae")
        self.assertEqual(family_of("Nanovirus sp."), "Nanoviridae")
        self.assertEqual(family_of("Tobamovirus sp."), "Virgaviridae")
        self.assertEqual(family_of("Potyvirus sp."), "Potyviridae")
        self.assertEqual(family_of("unknown mycovirus"), "Mycovirus")
        self.assertEqual(family_of("zzz no keyword"), "Other_viral")

    def test_known_species_level_gaps(self):
        # 已知局限 (见 FAMILY_KEYWORDS 注释): 种名级同义词不展开, 落 Other_viral
        # — 这是 "轻量映射" 的预期行为, 精确分类走 Stage3 RVDB
        self.assertEqual(family_of("Tobacco mosaic virus"), "Other_viral")
        self.assertEqual(family_of("Bean yellow dwarf virus"), "Other_viral")


class TestSummarizeAndMerge(unittest.TestCase):
    """模拟两个基因组的阶段产物 → 单基因组 summary + 跨基因组 merge."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "out"
        self._make_genome("G1", loci=[("chr1:99-700", "viral_supported"),
                                      ("chr1:50000-50700", "host_like")])
        self._make_genome("G2", loci=[("chr2:1-300", "undetermined")])

    def _make_genome(self, name, loci):
        base = self.out
        ld = base / "01_Loci" / name
        vd = base / "02_Verdict" / name
        rd = base / "03_RVDB" / name
        for d in (ld, vd, rd):
            d.mkdir(parents=True, exist_ok=True)
        bed = "".join(f"{loc.split(':')[0]}\t0\t100\t{loc}\n" for loc, _ in loci)
        write(ld / f"{name}.loci.bed", bed)
        write(ld / f"{name}.s1_best.tsv",
              "locus\tref_sseqid\tfamily\tref_evalue\tref_bitscore\tref_qcov\n"
              + "".join(f"{loc}\tNC_x\tGeminiviridae\t1e-9\t95\t80\n"
                        for loc, _ in loci))
        write(vd / f"{name}.s2_verdict.tsv",
              "locus\tverdict\tbest_viral_id\tbest_viral_bs\tbest_plant_id\tbest_plant_bs\n"
              + "".join(f"{loc}\t{v}\tSP_V\t120\tSP_P\t10\n" for loc, v in loci))
        write(rd / f"{name}.s3_rvdb.tsv",
              # OUTFMT2 11 列, c[9]=bitscore c[10]=stitle
              "".join(f"{loc}\tRVDB_x\t40\t100\t1\t100\t1\t100\t1e-9\t70\t"
                      f"Geminiviridae sp. [RVDB]\n" for loc, v in loci))

    def tearDown(self):
        self.tmp.cleanup()

    def test_summary_rows(self):
        summ = summarize_genome("G1", self.out)
        lines = summ.read_text().splitlines()
        self.assertEqual(len(lines), 3)  # header + 2 loci
        f = lines[1].split("\t")
        self.assertEqual(f[:3], ["G1", "chr1:99-700", "viral_supported"])
        self.assertEqual(f[3], "Geminiviridae")
        self.assertEqual(f[7], "SP_V")
        self.assertEqual(f[12], "70")  # rvdb_bitscore

    def test_merge_all(self):
        summarize_genome("G1", self.out)
        summarize_genome("G2", self.out)
        n_g, n_l = merge_all(self.out)
        self.assertEqual((n_g, n_l), (2, 3))
        ks = (self.out / "kingdom_summary.tsv").read_text().splitlines()
        self.assertEqual(ks[0], "genome\tviral_supported\thost_like\tundetermined")
        self.assertEqual(ks[1], "G1\t1\t1\t0")
        self.assertEqual(ks[2], "G2\t0\t0\t1")
        fam = (self.out / "family_by_genome.tsv").read_text().splitlines()
        self.assertEqual(fam[0], "family\tG1\tG2")
        self.assertEqual(fam[1], "Geminiviridae\t1\t0")  # 仅 viral_supported 计入


class TestLocusKeyJoin(unittest.TestCase):
    """loci.bed 写 chr1:50-699, 而抽序列工具/diamond 可能把 ':' 改写成 '_'.
    两侧都过 locus_key 后 join, 否则 Stage2 判定被静默丢弃 (全变 undetermined,
    Stage3 候选为 0) — 这是离线 mock E2E 抓到的真实回归."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "out"
        ld = self.out / "01_Loci" / "G"
        vd = self.out / "02_Verdict" / "G"
        rd = self.out / "03_RVDB" / "G"
        for d in (ld, vd, rd):
            d.mkdir(parents=True, exist_ok=True)
        write(ld / "G.loci.bed", "chr1\t49\t699\tchr1:50-699\n")
        write(ld / "G.s1_best.tsv",
              "locus\tref_sseqid\tfamily\tref_evalue\tref_bitscore\tref_qcov\n"
              "chr1:50-699\tNC_x\tGeminiviridae\t1e-9\t95\t80\n")
        # 下游键全部被 sanitize 成 '_' (工具行为不确定, 必须两边都兼容)
        write(vd / "G.s2_verdict.tsv",
              "locus\tverdict\tbest_viral_id\tbest_viral_bs\tbest_plant_id\tbest_plant_bs\n"
              "chr1_50-699\tviral_supported\tSP_V\t120\tSP_P\t10\n")
        write(rd / "G.s3_rvdb.tsv",
              "chr1_50-699\tRVDB_x\t40\t100\t1\t100\t1\t100\t1e-9\t70\t"
              "Geminiviridae sp. [RVDB]\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_key_canonical(self):
        self.assertEqual(locus_key("chr1:50-699"), locus_key("chr1_50-699"))

    def test_summary_joins_sanitized_keys(self):
        summ = summarize_genome("G", self.out)
        f = summ.read_text().splitlines()[1].split("\t")
        self.assertEqual(f[1], "chr1:50-699")   # 展示保持 loci.bed 原名
        self.assertEqual(f[2], "viral_supported")
        self.assertEqual(f[7], "SP_V")
        self.assertEqual(f[11], "Geminiviridae sp. [RVDB]")
        self.assertEqual(f[12], "70")


class TestParseStages(unittest.TestCase):
    def test_aliases(self):
        self.assertEqual(parse_stages("1,2,3"), {"discover", "verdict", "annotate"})
        self.assertEqual(parse_stages("discover,annotate"), {"discover", "annotate"})
        self.assertEqual(parse_stages("all"), {"discover", "verdict", "annotate"})
        with self.assertRaises(SystemExit):
            parse_stages("9")


class TestSlidingCmd(unittest.TestCase):
    """seqkit sliding: -g 是布尔开关 --greedy, 基因组是位置参数.

    曾把 '基因组路径' 塞给 -g 当值, 靠 pflag 漏值碰巧跑通 — 意图全丢且
    路径以 '-' 开头即失败. 这里锁定显式写法.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cfg = EveConfig(window=50000)
        self.g = Path(self.tmp.name) / "genome.fa"
        self.chunks = Path(self.tmp.name) / "chunks.fna"

    def tearDown(self):
        self.tmp.cleanup()

    def test_explicit_flags_and_positional_input(self):
        cmd = sliding_cmd(self.cfg, self.g, self.chunks, 50000)
        self.assertEqual(cmd, ["seqkit", "sliding", "-s", "50000",
                               "-W", "50000", "--greedy",
                               "-o", str(self.chunks), str(self.g)])

    def test_overlap_halves_step_only(self):
        cmd = sliding_cmd(self.cfg, self.g, self.chunks, 25000)
        self.assertEqual(cmd[cmd.index("-s") + 1], "25000")
        self.assertEqual(cmd[cmd.index("-W") + 1], "50000")

    def test_greedy_takes_no_value(self):
        # --greedy 后不能再跟"值": 基因组必须直接作为位置参数在末尾
        cmd = sliding_cmd(self.cfg, self.g, self.chunks, 50000)
        i = cmd.index("--greedy")
        self.assertEqual(cmd[i + 1], "-o")
        self.assertEqual(cmd[-1], str(self.g))
        self.assertNotIn("-g", cmd)   # 短选项 -g 不再直接出现于命令中


class TestCheckTools(unittest.TestCase):
    """工具体检: samtools 必查 (抽序列只用它); blastn 仅类病毒层需要."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.db = Path(self.tmp.name) / "db"
        self.db.mkdir()
        self.cfg = EveConfig(ref_dmnd=str(self.db / "r.dmnd"),
                             pv_dmnd=str(self.db / "p.dmnd"),
                             id2div=str(self.db / "d.tsv"),
                             rvdb_dmnd=str(self.db / "v.dmnd"))
        for p in ("r.dmnd", "p.dmnd", "d.tsv", "v.dmnd"):
            write(self.db / p, "x")

    def tearDown(self):
        self.tmp.cleanup()

    def test_all_ready(self):
        with mock.patch("eve_scan_core.shutil.which", return_value="/x"):
            self.assertEqual(check_tools(self.cfg), [])

    def test_missing_samtools_reported(self):
        with mock.patch("eve_scan_core.shutil.which",
                        side_effect=lambda n: None if n == "samtools" else "/x"):
            self.assertEqual(check_tools(self.cfg), ["samtools=samtools"])

    def test_configured_samtools_path_must_exist(self):
        self.cfg.samtools = str(self.db / "nope-samtools")

        def _which(n):
            # 已配置但文件不存在的路径才算缺; 裸工具名一律当作找得到
            if n.startswith(str(self.db)) and not Path(n).exists():
                return None
            return "/x"

        with mock.patch("eve_scan_core.shutil.which", side_effect=_which):
            self.assertIn(f"samtools={self.cfg.samtools}",
                          check_tools(self.cfg))
            write(self.cfg.samtools, "x")   # 路径存在即可, 不要求 PATH 有
            self.assertEqual(check_tools(self.cfg), [])

    def test_samtools_too_old_blocked_before_run(self):
        # 版本不够要在开跑前报出来, 不能等每个基因组 Stage1 才炸
        with mock.patch("eve_scan_core.shutil.which", return_value="/x"), \
                mock.patch("eve_scan_core.samtools_too_old", return_value="1.9"):
            bad = check_tools(self.cfg)
        self.assertTrue(any("samtools" in m and "1.9" in m for m in bad), bad)

    def test_blastn_required_for_viroid_layer(self):
        with mock.patch("eve_scan_core.shutil.which", return_value="/x"):
            self.assertEqual(check_tools(self.cfg), [])
            self.cfg.viroids_fa = str(self.db / "viroids.fa")
            write(self.cfg.viroids_fa, "x")
            self.assertEqual(check_tools(self.cfg, need_viroid=True), [])
        with mock.patch("eve_scan_core.shutil.which",
                        side_effect=lambda n: None if n == "blastn" else "/x"):
            self.assertIn("blastn=blastn",
                          check_tools(self.cfg, need_viroid=True))
            self.assertNotIn("blastn=blastn", check_tools(self.cfg))

    def test_missing_diamond_binary_reported(self):
        self.cfg.diamond = str(self.db / "no-diamond")
        with mock.patch("eve_scan_core.shutil.which", return_value=None):
            self.assertIn("diamond=" + self.cfg.diamond,
                          check_tools(self.cfg))

    def test_missing_db_reported(self):
        self.cfg.rvdb_dmnd = str(self.db / "gone.dmnd")
        with mock.patch("eve_scan_core.shutil.which", return_value="/x"):
            self.assertIn(f"rvdb_dmnd={self.cfg.rvdb_dmnd}",
                          check_tools(self.cfg))


class TestSamtoolsVersionProbe(unittest.TestCase):
    """samtools 版本探测: 读第一行的版本号; 探测失败不拦任务 (交给 samtools 兜底)."""

    def _patch_version(self, stdout):
        return mock.patch("eve_scan_core.subprocess.run",
                          return_value=subprocess.CompletedProcess(
                              [], 0, stdout, ""))

    def test_current_version_accepted(self):
        with self._patch_version("samtools 1.21\nUsing htslib 1.21\n"):
            self.assertIsNone(samtools_too_old("/x/samtools"))

    def test_old_version_reported(self):
        with self._patch_version("samtools 1.9\nUsing htslib 1.9\n"):
            self.assertEqual(samtools_too_old("/x/samtools"), "1.9")

    def test_htslib_line_not_mistaken_for_samtools(self):
        # 版本号取自第一行; 不能抓到 Copyright 里的年份或 htslib 版本
        with self._patch_version("samtools 1.10\nUsing htslib 1.10\n"
                                 "Copyright (C) 2024 Genome Research Ltd.\n"):
            self.assertEqual(samtools_too_old("/x/samtools"), "1.10")

    def test_probe_failure_does_not_block(self):
        with mock.patch("eve_scan_core.subprocess.run",
                        side_effect=OSError("no exec")):
            self.assertIsNone(samtools_too_old("/x/samtools"))
        with self._patch_version("?? 不是版本行 ??\n"):
            self.assertIsNone(samtools_too_old("/x/samtools"))


class TestBedToRegions(unittest.TestCase):
    """BED(0基半开) → samtools 区域文件(1基闭区间): 坐标差一位全批位点静默错位."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.bed = Path(self.tmp.name) / "loci.bed"
        self.regions = Path(self.tmp.name) / "loci.fa.regions.txt"

    def tearDown(self):
        self.tmp.cleanup()

    def test_one_based_inclusive_conversion(self):
        write(self.bed, "chr1\t49\t699\tchr1:50-699\t0\t+\n"
                        "chr2\t0\t300\tchr2:1-300\t0\t-\n")
        self.assertEqual(bed_to_regions(self.bed, self.regions), 2)
        self.assertEqual(self.regions.read_text(encoding="utf-8"),
                         "chr1:50-699\nchr2:1-300\n")

    def test_name_column_not_used(self):
        # 区域名取自坐标本身; 第4列叫什么都不影响 (samtools 以区域串作输出头)
        write(self.bed, "chr1\t10\t20\tWHATEVER\n")
        self.assertEqual(bed_to_regions(self.bed, self.regions), 1)
        self.assertEqual(self.regions.read_text(encoding="utf-8"),
                         "chr1:11-20\n")

    def test_skips_comments_short_and_inverted(self):
        write(self.bed, "# comment\ntrack name=x\n"
                        "chr1\t100\n"                  # 列数不足
                        "chr1\t500\t400\tbad\n"        # 起止倒挂
                        "chr1\t400\t400\tempty\n"      # 零长度
                        "chr1\tabc\t200\tnonnum\n"     # 非整数
                        "chr3\t999\t1200\tok\n")
        self.assertEqual(bed_to_regions(self.bed, self.regions), 1)
        self.assertEqual(self.regions.read_text(encoding="utf-8"),
                         "chr3:1000-1200\n")

    def test_empty_bed_yields_no_regions(self):
        write(self.bed, "")
        self.assertEqual(bed_to_regions(self.bed, self.regions), 0)
        self.assertEqual(self.regions.read_text(encoding="utf-8"), "")


class TestExtractCmd(unittest.TestCase):
    """samtools faidx 命令: 区域走文件, 索引写输出目录."""

    def setUp(self):
        self.cfg = EveConfig(ref_dmnd="x", pv_dmnd="x", id2div="x",
                             rvdb_dmnd="x")

    def test_cmd_shape(self):
        with mock.patch("eve_scan_core.shutil.which",
                        return_value="/usr/bin/samtools"):
            cmd = extract_cmd(self.cfg, "/db/g.fa", "/out/o.genome.fai",
                              "/out/loci.fa.regions.txt", "/out/loci.fa")
        self.assertEqual(cmd[0], "/usr/bin/samtools")   # 未配置 = PATH 解析结果
        self.assertEqual(cmd[1], "faidx")
        self.assertEqual(cmd[2], "/db/g.fa")
        self.assertIn("--fai-idx", cmd)
        self.assertEqual(cmd[cmd.index("--fai-idx") + 1], "/out/o.genome.fai")
        # 必须用正式选项名 --fai-idx: --fai 只靠无歧义前缀匹配才碰巧能用
        self.assertNotIn("--fai", cmd)
        self.assertEqual(cmd[cmd.index("--region-file") + 1],
                         "/out/loci.fa.regions.txt")
        self.assertEqual(cmd[cmd.index("--output") + 1], "/out/loci.fa")
        # 区域不当命令行参数: 万级位点会撞 ARG_MAX
        self.assertNotIn("chr1:50-699", cmd)

    def test_configured_samtools_used(self):
        with tempfile.TemporaryDirectory() as td:
            exe = Path(td) / "samtools"
            write(exe, "x")
            self.cfg.samtools = str(exe)
            cmd = extract_cmd(self.cfg, "/db/g.fa", "/out/f.fai", "/out/r.txt",
                              "/out/loci.fa")
            self.assertEqual(cmd[0], str(exe))


class TestCheckExtracted(unittest.TestCase):
    """抽序列产物核对: 越界/截断区域必须先炸, 不能静默流到下游.

    真 samtools 1.21 对这些情况 exit 仍是 0 (只在 stderr 写一句
    [faidx] Zero length/Truncated sequence), run() 的退出码检查拦不住.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.regions = Path(self.tmp.name) / "r.regions.txt"
        self.out = Path(self.tmp.name) / "loci.fa"

    def tearDown(self):
        self.tmp.cleanup()

    def _fa(self, recs):
        """按真 samtools 的形态写输出: 头=区域串, 每条记录一行不折行."""
        write(self.out, "".join(f">{h}\n{s}\n" for h, s in recs))

    def test_matching_passes(self):
        write(self.regions, "chr1:1-10\nchr1:50-699\nchr2:5-6\n")
        self._fa([("chr1:1-10", "A" * 10),
                  ("chr1:50-699", "C" * 650),
                  ("chr2:5-6", "G" * 2)])
        self.assertEqual(check_extracted(self.regions, self.out), 3)

    def test_wrapped_sequences_counted_by_length(self):
        write(self.regions, "chr1:1-10\n")
        write(self.out, ">chr1:1-10\n" + "AAAAA\nAAAAA\n")   # 折不折行只看累计长度
        self.assertEqual(check_extracted(self.regions, self.out), 1)

    def test_underscore_headers_accepted(self):
        # 下游 join 用 locus_key() 归一化 ':'/'_', 核对也得用同一把尺子
        write(self.regions, "chr1:50-699\n")
        self._fa([("chr1_50-699", "C" * 650)])
        self.assertEqual(check_extracted(self.regions, self.out), 1)

    def test_zero_length_region_raises(self):
        write(self.regions, "chr1:12346-19346\nchr1:1-10\n")
        self._fa([("chr1:12346-19346", ""),        # 真 samtools 的越界产物
                  ("chr1:1-10", "A" * 10)])
        with self.assertRaises(RuntimeError) as ctx:
            check_extracted(self.regions, self.out)
        msg = str(ctx.exception)
        self.assertIn("chr1:12346-19346", msg)
        self.assertIn("0bp != 7001bp", msg)
        self.assertIn("1/2", msg)

    def test_truncated_region_raises(self):
        write(self.regions, "chr1:2302-3000\n")
        self._fa([("chr1:2302-3000", "A")])        # 染色体只有 2302bp
        with self.assertRaises(RuntimeError) as ctx:
            check_extracted(self.regions, self.out)
        self.assertIn("1bp != 699bp", str(ctx.exception))

    def test_missing_record_raises(self):
        write(self.regions, "chr1:1-10\nchr1:50-699\n")
        self._fa([("chr1:1-10", "A" * 10)])        # 第二条整个丢了
        with self.assertRaises(RuntimeError) as ctx:
            check_extracted(self.regions, self.out)
        self.assertIn("chr1:50-699 无记录", str(ctx.exception))

    def test_error_message_truncated_to_five(self):
        write(self.regions, "".join(f"chr1:{i + 1}-{i + 2}\n" for i in range(8)))
        self._fa([])                               # 全部无记录
        with self.assertRaises(RuntimeError) as ctx:
            check_extracted(self.regions, self.out)
        self.assertIn("...", str(ctx.exception))


class TestRunTimeout(unittest.TestCase):
    """run() 超时: 挂死的外部命令必须杀掉并报错, 不能白占并行核."""

    def test_timeout_kills_and_raises(self):
        cmd = [sys.executable, "-c", "import time; time.sleep(30)"]
        with self.assertRaises(RuntimeError) as ctx:
            run(cmd, _Log(), timeout=1)
        self.assertIn("TIMEOUT", str(ctx.exception))

    def test_success_returns_normally(self):
        cmd = [sys.executable, "-c", "print('ok')"]
        run(cmd, _Log(), timeout=30)

    def test_nonzero_exit_raises(self):
        cmd = [sys.executable, "-c", "import sys; sys.exit(3)"]
        with self.assertRaises(RuntimeError) as ctx:
            run(cmd, _Log())
        self.assertIn("exit 3", str(ctx.exception))


class TestResolveGenomeFa(unittest.TestCase):
    """跳过 Stage1 时推断 fasta: 压缩输入必须落在 Stage1 解压产物上."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "out"
        self.g = self.out / "01_Loci" / "G" / "G.fna"
        self.g.parent.mkdir(parents=True)
        all_stages = {"discover", "verdict", "annotate"}
        self.no_discover = {"verdict", "annotate"}

    def tearDown(self):
        self.tmp.cleanup()

    def test_discover_stage_returns_none(self):
        self.assertIsNone(resolve_genome_fa("/x/y.fa", self.out, "G",
                                            {"discover", "verdict"}))

    def test_plain_fasta_preferred(self):
        raw = self.out / "raw.fa"
        write(raw, ">a\nACGT\n")
        self.assertEqual(resolve_genome_fa(str(raw), self.out, "G",
                                           self.no_discover), raw)

    def test_archive_uses_unpacked_fna(self):
        write(self.g, ">a\nACGT\n")
        self.assertEqual(resolve_genome_fa("/x/G.tar.gz", self.out, "G",
                                           self.no_discover), self.g)
        self.assertEqual(resolve_genome_fa("/x/G.fa.gz", self.out, "G",
                                           self.no_discover), self.g)

    def test_archive_without_unpacked_raises(self):
        with self.assertRaises(FileNotFoundError) as ctx:
            resolve_genome_fa("/x/G.tar.gz", self.out, "G", self.no_discover)
        self.assertIn("stage discover", str(ctx.exception))

    def test_missing_plain_falls_back_to_unpacked(self):
        write(self.g, ">a\nACGT\n")
        self.assertEqual(resolve_genome_fa("/x/gone.fa", self.out, "G",
                                           self.no_discover), self.g)

    def test_missing_plain_without_unpacked_returns_original(self):
        p = self.out / "gone.fa"
        self.assertEqual(resolve_genome_fa(str(p), self.out, "G",
                                           self.no_discover), p)


class TestBatchJobs(unittest.TestCase):
    """eve_screen.build_jobs: batch.tsv 解析 + 重名必须开跑前挡住."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        from eve_screen import build_jobs, build_parser
        self.build_jobs, self.build_parser = build_jobs, build_parser

    def tearDown(self):
        self.tmp.cleanup()

    def _args(self, lines):
        tsv = Path(self.tmp.name) / "batch.tsv"
        write(tsv, "".join(lines))
        return self.build_parser().parse_args(
            ["-B", str(tsv), "-o", str(Path(self.tmp.name) / "out")])

    _STAGES = {"discover", "verdict", "annotate"}

    def test_crlf_comments_and_columns(self):
        args = self._args(["A\t/p/a.fa\r\n", "# 注释\r\n", "坏行\r\n",
                           "B\t/p/b.fa\n"])
        jobs = self.build_jobs(args, self._STAGES, _Log())
        self.assertEqual([j["name"] for j in jobs], ["A", "B"])
        self.assertEqual(jobs[0]["genome"], "/p/a.fa")   # \r 被 strip 掉

    def test_duplicate_names_abort(self):
        args = self._args(["A(1)\t/p/a.fa\n", "A_1_\t/p/b.fa\n", "A\t/p/c.fa\n"])
        with self.assertRaises(SystemExit) as ctx:
            self.build_jobs(args, self._STAGES, _Log())
        self.assertIn("重名", str(ctx.exception))

    def test_distinct_names_ok(self):
        args = self._args(["A\t/p/a.fa\n", "B\t/p/b.fa\n"])
        jobs = self.build_jobs(args, self._STAGES, _Log())
        self.assertEqual(len(jobs), 2)
        self.assertEqual({j["name"] for j in jobs}, {"A", "B"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
