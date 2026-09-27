#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_eve_refgenome_scan.py — S2c 参考基因组通道 + S3 rg 单向调整的判定测试
==========================================================================
v6.3 新通道的三层锁定 (夹具 helpers 复用 test_eve_distinguish_logic, 不另抄一份):

  s2c classify 的结构判定:
    两端侧翼同一位点 -> EVE_flank_confirmed / 位点对不上或间隔超限 -> flank_locus_conflict
    内部宿主段剪接样式 (含 与端侧翼并存) -> EVE_splice_chimera / 单独内部段 -> host_mix_weak
    host_dominant 先于一切 (整条 contig 就是寄主序列时, "侧翼"是它自己)
    低一致/短段 HSP 不构成结构证据 (no_hit)
  s2c 装配:
    no_sample_map / genome_absent (无 #BLASTED 标记) 与 no_hit 三种"未判"分得开;
    列数不足整表报错而不是静默变 no_hit; 同一输入两次运行逐字节相同 (行序可 diff)
  s3 rg_adjust 单向性:
    只把 review 档往 EVE 升 / 矛盾转 host_conflict_review; no_hit 与未给 s2c
    绝不改 verdict ("没比对上"不是真病毒的证据); VERDICT_HDR 带 rg_call/rg_scope;
    s4 的 CARRY 把 rg 口径带给下游.

运行: python -m unittest discover -s endogenous_virus_pipeline/tests -v
  或: python endogenous_virus_pipeline/tests/test_eve_refgenome_scan.py
"""
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]                      # endogenous_virus_pipeline/
MOD = ROOT / 'eve_distinguish'
sys.path.insert(0, str(MOD))
sys.path.insert(0, str(HERE))               # 复用兄弟测试文件的夹具 helpers

import s2c_refgenome_scan as s2c            # noqa: E402
from s1_decay_scan import S1_HDR            # noqa: E402
from s2b_locus_scan import HDR as S2B_HDR   # noqa: E402
from s2_domain_scan import S2_HDR           # noqa: E402
from s3_verdict import VERDICT_HDR, VERDICT_COLS, rg_adjust  # noqa: E402
from s4_filter import CARRY, OUT_HDR        # noqa: E402
from test_eve_distinguish_logic import (    # noqa: E402
    EV_HDR, run_py, read_tsv_rows, s1_row, s2_row, write)


# ---------------------------------------------------------------- 夹具
def rg_row(code, qseqid, sseqid, pid="99.0", aln="400", qs="1", qe="400",
           ss="100001", se="100400", qlen="5000", bits="800", slen="100000000"):
    """refgenome_blastn.tsv 的一行: 代码 + outfmt 6 的 14 列 (mismatch/gaps 置 0)."""
    return "\t".join([code, qseqid, sseqid, pid, aln, "0", "0", qs, qe,
                      ss, se, "1e-100", bits, qlen, slen]) + "\n"


def hsp(qs, qe, sseqid="chr1", sm=100001, sM=None, strand="+", pid=99.0,
        aln=None, qlen=5000, bits=800):
    """直接构造 classify() 的输入 (run_all 落表后的等价形态)."""
    if aln is None:
        aln = qe - qs + 1
    if sM is None:
        sM = sm + aln - 1
    return {"code": "GBCQ", "sseqid": sseqid, "pid": pid, "aln": aln,
            "qs": min(qs, qe), "qe": max(qs, qe),
            "sm": min(sm, sM), "sM": max(sm, sM),
            "strand": strand, "bits": bits, "qlen": qlen}


def s2b_row(cid, arch="no_host_locus", scope="none", wpid="", cov="",
            ncomp=""):
    cols = S2B_HDR.rstrip("\n").split("\t")
    vals = {"contig_id": cid, "sample": cid.split("_clean_")[0],
            "sample_flag": "clean", "host_scaffold": "-", "host_code": "",
            "host_scope": scope, "host_wpid": wpid, "host_cov": cov,
            "host_xeno": "-", "n_host_scaffolds": "0", "comps_contig": "-",
            "best_locus_comps": "-", "decay_class": "", "te_flag": "FALSE",
            "locus_key": "k", "locus_arch": arch, "locus_ncomp": ncomp,
            "locus_comps": "-", "locus_members": "1",
            "locus_hint": "no_host_mapping" if arch == "no_host_locus" else "x"}
    return "\t".join(vals[c] for c in cols) + "\n"


def ev_row(cid, fam="Caulimoviridae", aa_pid="", aa_spec="", nt_pid="", nt_spec=""):
    return (f"{cid}\t{fam}\tvirus\t60\t{aa_pid}\t{aa_spec}\t{nt_pid}\t{nt_spec}\n")


RG_HDR_LINE = s2c.RG_HDR


# ══════════════════════════════════════════════════ S2c classify (纯函数)
class TestClassifyFlanks(unittest.TestCase):
    """判据一: 两端宿主侧翼。同染色体同链 + 间隔有限 = 前病毒; 其余一律不定案."""

    def test_flank_pair_same_locus_is_confirmed(self):
        hsps = [hsp(1, 400, sm=100001),            # 左端 400nt, chr1:+
                hsp(4601, 5000, sm=104601)]        # 右端 400nt, 间隔 4200nt
        call, info = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertEqual(call, "EVE_flank_confirmed")
        self.assertEqual(info["flank_same_locus"], "TRUE")
        self.assertEqual(info["flank_gap"], 4200)
        self.assertTrue(info["flank_left"].startswith("chr1:"))
        self.assertIn("(+)", info["flank_right"])

    def test_minus_strand_contig_also_confirmed(self):
        # 候选拼成了参考基因组的负链: 两条侧翼 HSP 都应为 '-', 间隔按负链序数
        hsps = [hsp(1, 400, sm=104601, sM=105000, strand="-"),
                hsp(4601, 5000, sm=100001, sM=100400, strand="-")]
        call, info = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertEqual(call, "EVE_flank_confirmed")
        self.assertEqual(info["flank_gap"], 4200)

    def test_flanks_on_different_chromosomes_are_conflict(self):
        hsps = [hsp(1, 400, sseqid="chr1", sm=100001),
                hsp(4601, 5000, sseqid="chr2", sm=104601)]
        call, info = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertEqual(call, "flank_locus_conflict")
        self.assertEqual(info["flank_same_locus"], "FALSE")

    def test_flanks_too_far_apart_are_conflict(self):
        # 间隔 > MAX_FLANK_GAP: 同一条染色体, 但不是"围着同一个整合位点"的距离
        hsps = [hsp(1, 400, sm=100001),
                hsp(4601, 5000, sm=100001 + s2c.MAX_FLANK_GAP + 4200)]
        call, _ = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertEqual(call, "flank_locus_conflict")

    def test_overlapping_strands_conflict(self):
        # 同一染色体但异链: 不是一条链上跨插段的两个侧翼
        hsps = [hsp(1, 400, sm=100001, strand="+"),
                hsp(4601, 5000, sm=104601, strand="-")]
        call, _ = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertEqual(call, "flank_locus_conflict")


class TestClassifyChimeraAndHostdom(unittest.TestCase):
    """判据二 (外显子嵌合) 与 host_dominant 的优先级."""

    def test_internal_splice_pattern_is_chimera(self):
        # 两段内部宿主段, 同染色体同链, subject 间隔 >= SPLICE_GAP = 剪接样式
        hsps = [hsp(1000, 1300, sm=500001),
                hsp(2000, 2300, sm=600001)]
        call, info = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertEqual(call, "EVE_splice_chimera")
        self.assertEqual(info["splice_pattern"], "TRUE")
        self.assertEqual(info["internal_host_segs"], 2)

    def test_internal_plus_end_flank_is_chimera(self):
        hsps = [hsp(1, 300, sm=100001),            # 端侧翼
                hsp(2000, 2300, sm=500001)]        # 内部段 (无第二条内部段可比剪接)
        call, _ = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertEqual(call, "EVE_splice_chimera")

    def test_single_internal_segment_is_only_weak(self):
        # 一条内部段不满足任何样式: 组装嵌合和真实嵌合转录分不开, 只记列
        call, _ = s2c.classify([hsp(2000, 2300, sm=500001)], has_mp_cp_ap=True)
        self.assertEqual(call, "host_mix_weak")

    def test_host_dominant_beats_everything(self):
        # 整条 contig 96% 被一段 99% 一致的宿主段盖住且无病毒组件: 这是寄主序列,
        # 即便它同时"两端触边" (single_span, 侧翼就是它自己)
        hsps = [hsp(1, 4800, sm=100001, pid=99.0)]
        call, info = s2c.classify(hsps, has_mp_cp_ap=False)
        self.assertEqual(call, "host_dominant")
        self.assertEqual(info["host_seg_cov"], 0.96)

    def test_host_dominant_requires_no_viral_comps(self):
        # 同样的宿主覆盖, 但带 MP/CP/AP: 不能判死, 落回嵌合/结构通道
        hsps = [hsp(1, 4800, sm=100001, pid=99.0)]
        call, _ = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertNotEqual(call, "host_dominant")

    def test_low_pidentity_hits_are_not_evidence(self):
        # MIN_SEG_PID 以下的命中不该由 classify 处理 —— 装配层过滤 (no_hit);
        # 这里锁的是: 即便混进来, 也决不能凑成"侧翼同位点"
        hsps = [hsp(1, 400, pid=85.0), hsp(4601, 5000, pid=85.0)]
        call, info = s2c.classify(hsps, has_mp_cp_ap=True)
        self.assertIn(call, ("flank_locus_conflict", "no_structure"))
        self.assertEqual(info["flank_same_locus"], "FALSE")

    def test_merge_keeps_different_chromosomes_separate(self):
        # query 上重叠的两条 HSP 分属不同染色体: 不许并成一段 (会把转位洗成同位点)
        hsps = [hsp(1, 400, sseqid="chr1", sm=100001),
                hsp(1, 400, sseqid="chr2", sm=200001)]
        segs = s2c.merge_segments(hsps)
        self.assertEqual(len(segs), 2)
        self.assertEqual({s["sseqid"] for s in segs}, {"chr1", "chr2"})


# ══════════════════════════════════════════════════ S2c 装配 (文件级)
class TestS2cAssembly(unittest.TestCase):
    def _fixture(self, d):
        # 三种归属口径的样本: ERR1 有比对过的参考基因组 / ERR2 的代码没比过 /
        # ERR9 压根不在映射表里
        host_map = write(d / "host_map.tsv",
                         "ERR1\tGBCQ\tGBCQ-Aerva_persica\n"
                         "ERR2\tZZZZ\tZZZZ-nomap\n")
        s2 = write(d / "s2_domains.tsv",
                   S2_HDR +
                   s2_row("ERR1_clean_N1", "MP+CP", "2", "MP+CP") +
                   s2_row("ERR2_clean_N1", "MP", "1", "MP") +
                   s2_row("ERR9_clean_N1", "MP", "1", "MP"))
        return host_map, s2

    def test_calls_and_scope_split(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            host_map, s2 = self._fixture(d)
            blast = write(d / "refgenome_blastn.tsv",
                          # ERR1 的 GBCQ 参考基因组: 两条侧翼同一位点
                          rg_row("GBCQ", "ERR1_clean_N1", "chr1", qs="1", qe="400",
                                 ss="100001", se="100400") +
                          rg_row("GBCQ", "ERR1_clean_N1", "chr1", qs="4601", qe="5000",
                                 ss="104601", se="105000") +
                          "#BLASTED\tGBCQ\n" +
                          # ERR2 的代码 ZZZZ 没有任何标记行 -> genome_absent
                          "#BLASTED\tXXXX\n")
            rc = run_py("s2c_refgenome_scan.py", blast, s2, host_map,
                        d / "refgenome_evidence.tsv")
            self.assertEqual(rc.returncode, 0, rc.stderr)
            rows = {r["contig_id"]: r
                    for r in read_tsv_rows(d / "refgenome_evidence.tsv")}
            self.assertEqual(rows["ERR1_clean_N1"]["rg_call"], "EVE_flank_confirmed")
            self.assertEqual(rows["ERR1_clean_N1"]["rg_scope"], "own_genome")
            self.assertEqual(rows["ERR1_clean_N1"]["flank_gap"], "4200")
            # 在映射表里但代码没比过: genome_absent; 不在映射表: no_sample_map
            self.assertEqual(rows["ERR2_clean_N1"]["rg_call"], "no_ref_genome")
            self.assertEqual(rows["ERR2_clean_N1"]["rg_scope"], "genome_absent")
            self.assertEqual(rows["ERR9_clean_N1"]["rg_scope"], "no_sample_map")
            # 比过但零命中 (GBCQ 有标记、候选无行): 不与 genome_absent 混
            blast2 = write(d / "refgenome_blastn_zero.tsv",
                           "#BLASTED\tGBCQ\n")
            rc2 = run_py("s2c_refgenome_scan.py", blast2, s2, host_map,
                         d / "rg_zero.tsv")
            self.assertEqual(rc2.returncode, 0, rc2.stderr)
            rows2 = {r["contig_id"]: r for r in read_tsv_rows(d / "rg_zero.tsv")}
            self.assertEqual(rows2["ERR1_clean_N1"]["rg_call"], "no_hit")
            self.assertEqual(rows2["ERR1_clean_N1"]["rg_scope"], "own_genome")
            self.assertEqual(rows2["ERR2_clean_N1"]["rg_scope"], "genome_absent")
            self.assertEqual(rows2["ERR9_clean_N1"]["rg_scope"], "no_sample_map")

    def test_deterministic_row_order(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            host_map, s2 = self._fixture(d)
            lines = [rg_row("GBCQ", "ERR1_clean_N1", "chr1", qs="1", qe="400"),
                     rg_row("GBCQ", "ERR1_clean_N1", "chr1", qs="4601", qe="5000",
                            ss="104601", se="105000"),
                     "#BLASTED\tGBCQ\n"]
            for i, cid in enumerate([f"ERR1_clean_N{i}" for i in range(2, 12)]):
                lines.append(rg_row("GBCQ", cid, "chr3", qs="10", qe="409",
                                    ss="1", se="400"))
            blast = write(d / "b.tsv", "".join(lines))
            out1, out2 = d / "o1.tsv", d / "o2.tsv"
            self.assertEqual(run_py("s2c_refgenome_scan.py", blast, s2, host_map,
                                    out1).returncode, 0)
            self.assertEqual(run_py("s2c_refgenome_scan.py", blast, s2, host_map,
                                    out2).returncode, 0)
            self.assertEqual(out1.read_bytes(), out2.read_bytes())

    def test_bad_columns_abort_instead_of_silent_no_hit(self):
        # 列数不对的表若被当成"零命中", 全部候选静默变 no_hit —— 必须报错拦下
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            host_map, s2 = self._fixture(d)
            blast = write(d / "b.tsv", "GBCQ\tERR1_clean_N1\tchr1\t99.0\n")
            rc = run_py("s2c_refgenome_scan.py", blast, s2, host_map,
                        d / "o.tsv")
            self.assertNotEqual(rc.returncode, 0)
            self.assertIn("outfmt", rc.stderr)

    def test_emit_header_matches_module_constant(self):
        rc = run_py("s2c_refgenome_scan.py", "--emit-header")
        self.assertEqual(rc.returncode, 0)
        self.assertEqual(rc.stdout, RG_HDR_LINE)


# ══════════════════════════════════════════════════ S3 rg_adjust (单向性)
class TestRgAdjust(unittest.TestCase):
    def test_flank_upgrades_review_only(self):
        self.assertEqual(rg_adjust("review", "EVE_flank_confirmed", False, False),
                         "EVE_STRONG_provirus")
        self.assertEqual(rg_adjust("EVE_suspect", "EVE_flank_confirmed", False, False),
                         "EVE_STRONG_provirus")
        self.assertEqual(rg_adjust("compositional_noise_review",
                                   "EVE_flank_confirmed", False, False),
                         "EVE_STRONG_provirus")

    def test_splice_upgrades_to_suspect_not_strong(self):
        # 判据二是 contig 级代理: 最多升到 EVE_suspect, 不冒充前病毒铁证
        self.assertEqual(rg_adjust("review", "EVE_splice_chimera", False, False),
                         "EVE_suspect")

    def test_flank_vs_virus_candidate_is_conflict(self):
        # 侧翼证据说"整合的", 二维表说"感染中的真病毒": 矛盾, 人工定案
        self.assertEqual(rg_adjust("virus_candidate", "EVE_flank_confirmed",
                                   False, False),
                         "host_conflict_review")

    def test_flank_vs_remove_is_conflict(self):
        self.assertEqual(rg_adjust("host_contamination_likely", "EVE_flank_confirmed",
                                   False, False),
                         "host_conflict_review")

    def test_host_dominant_follows_three_gate(self):
        self.assertEqual(rg_adjust("review", "host_dominant", False, False),
                         "host_contamination_likely")
        self.assertEqual(rg_adjust("host_homology_cross_species", "host_dominant",
                                   False, False),
                         "host_contamination_likely")
        # 上游有病毒信号: 一票转矛盾复核, 不删数据 (v6 的教训)
        self.assertEqual(rg_adjust("review", "host_dominant", True, False),
                         "host_conflict_review")
        self.assertEqual(rg_adjust("review", "host_dominant", False, True),
                         "host_conflict_review")

    def test_no_hit_and_unknown_calls_never_change_verdict(self):
        for v in ("review", "EVE_STRONG_provirus", "virus_candidate",
                  "host_contamination_likely"):
            for c in ("no_hit", "no_structure", "host_mix_weak",
                      "flank_locus_conflict", "no_ref_genome", ""):
                self.assertEqual(rg_adjust(v, c, False, False), v, (v, c))

    def test_strong_verdict_not_downgraded(self):
        self.assertEqual(rg_adjust("EVE_STRONG_provirus", "EVE_splice_chimera",
                                   False, False),
                         "EVE_STRONG_provirus")


class TestS3Integration(unittest.TestCase):
    """s3 带/不带第 6 参的端到端: rg 列存在、调整生效、不给 s2c 时行为不变."""

    def _run_s3(self, d, rg_text=None):
        # s3.read_tsv 把每个文件的首行当表头: 四张输入表都必须带各自的 HDR,
        # 否则整表错位 (测过: 漏表头时 flank 升级会误判到错误的候选上)
        s1 = write(d / "s1.tsv",
                   S1_HDR +
                   s1_row("R1", "compositional_noise") +
                   s1_row("R2", "compositional_noise") +
                   s1_row("R3", "coding_intact") +
                   s1_row("R4", "compositional_noise"))
        s2 = write(d / "s2.tsv",
                   S2_HDR +
                   s2_row("R1", "MP", "1", "MP") +
                   s2_row("R2", "MP", "1", "MP") +
                   s2_row("R3", "MP+CP+AP+RT", "4", "MP+CP+AP+RT") +
                   s2_row("R4", "MP", "1", "MP"))
        s2b = write(d / "s2b.tsv",
                    S2B_HDR +
                    s2b_row("R1") + s2b_row("R2") +
                    s2b_row("R3", arch="locus_full", scope="own_species",
                            wpid="95.0", cov="0.9", ncomp="4") +
                    s2b_row("R4"))
        ev = write(d / "ev.tsv",
                   EV_HDR +
                   ev_row("R1") + ev_row("R2") + ev_row("R3") + ev_row("R4",
                                                                         fam="-"))
        rg = None
        if rg_text is not None:
            rg = write(d / "rg.tsv", rg_text)
        args = ["s3_verdict.py", s1, s2, s2b, ev, d / "verdict.tsv"]
        if rg is not None:
            args.append(rg)
        rc = run_py(*args)
        self.assertEqual(rc.returncode, 0, rc.stderr)
        return {r["contig_id"]: r for r in read_tsv_rows(d / "verdict.tsv")}

    def test_without_rg_table_columns_still_present(self):
        with tempfile.TemporaryDirectory() as td:
            rows = self._run_s3(Path(td))
            self.assertEqual(VERDICT_COLS[-2:], ["verdict", "verdict_reason"])
            self.assertIn("rg_call", VERDICT_COLS)
            self.assertIn("rg_scope", VERDICT_COLS)
            self.assertEqual(rows["R1"]["rg_call"], "")
            self.assertEqual(rows["R1"]["verdict"], "compositional_noise_review")

    def test_flank_confirmed_upgrades_review(self):
        rg_text = (RG_HDR_LINE +
                   "R1\tERR1\tclean\tGBCQ\town_genome\tEVE_flank_confirmed\t"
                   "chr1:100001-100400(+)\tchr1:104601-105000(+)\tTRUE\t4200\t"
                   "FALSE\t0\t0.16\t99.0\tMP\t-\t两端宿主侧翼落于同一条染色体同链\n")
        with tempfile.TemporaryDirectory() as td:
            rows = self._run_s3(Path(td), rg_text)
            self.assertEqual(rows["R1"]["verdict"], "EVE_STRONG_provirus")
            self.assertIn("EVE_flank_confirmed", rows["R1"]["verdict_reason"])
            self.assertIn("s2c 调整", rows["R1"]["verdict_reason"])
            # 对照: 没有 rg 行的 R4 保持原判
            self.assertEqual(rows["R4"]["verdict"], "compositional_noise_review")

    def test_host_dominant_with_upstream_family_is_conflict(self):
        rg_text = (RG_HDR_LINE +
                   "R2\tERR1\tclean\tGBCQ\town_genome\thost_dominant\t-\t-\t"
                   "FALSE\t\tFALSE\t0\t0.92\t95.0\t-\t-\t聚合覆盖占优\n")
        with tempfile.TemporaryDirectory() as td:
            rows = self._run_s3(Path(td), rg_text)
            # R2 上游 tax_family=Caulimoviridae: 寄主否决不许一票定生死
            self.assertEqual(rows["R2"]["verdict"], "host_conflict_review")

    def test_virus_candidate_plus_flank_is_conflict(self):
        rg_text = (RG_HDR_LINE +
                   "R3\tERR1\tclean\tGBCQ\town_genome\tEVE_flank_confirmed\t"
                   "chr1:100001-100400(+)\tchr1:104601-105000(+)\tTRUE\t4200\t"
                   "FALSE\t0\t0.16\t99.0\tMP\t-\t两端宿主侧翼\n")
        with tempfile.TemporaryDirectory() as td:
            rows = self._run_s3(Path(td), rg_text)
            self.assertEqual(rows["R3"]["verdict"], "host_conflict_review")

    def test_s4_carries_rg_columns_before_verdict(self):
        self.assertIn("rg_call", CARRY)
        self.assertIn("rg_scope", CARRY)
        self.assertLess(OUT_HDR.index("rg_call"), OUT_HDR.index("verdict"))


if __name__ == '__main__':
    unittest.main()
