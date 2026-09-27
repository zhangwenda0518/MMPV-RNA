#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_eve_qc_contam.py — 测序污染筛查的测试
==========================================
筛查的价值全在"能不能把 phiX 认出来、又不会把真病毒一起误杀"，所以测试盯三件事:
  1. 判定本身: phiX/大肠杆菌噬菌体参考必须命中; 真病毒参考必须不命中;
  2. 统计口径: `contam_frac` 是**该 contig 全部命中里**污染命中的占比, 不是绝对条数
     (单看条数会把"打中一次 phiX 的长 contig"和"整条都是 phiX"混为一谈);
  3. 可复现: 同一份输入跑两次, 产物逐字节相同。
"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import eve_qc_contam as qc          # noqa: E402
from eve_scan_core import (CONTAM_HDR, STAGE_DIRS, contamination_scan,
                          is_contaminant_hit)   # noqa: E402

PHIX = "NP_040704.1"
PHIX_TITLE = "DNA replication initiation [Escherichia phage phiX174]"
CAMV = "sp|P03542.1|CAPSD_CAMVS"
CAMV_TITLE = "Capsid protein [Cauliflower mosaic virus]"


def raw_row(q, sid, title, bits):
    """造一行 OUTFMT1 (12 列). 列序见 eve_scan_core.OUTFMT1。"""
    return "\t".join([q, sid, "90", "100", "1", "300", "1", "300",
                      "1e-50", str(bits), "1", title])


class TestContaminantPredicate(unittest.TestCase):
    def test_phix_accessions_are_flagged(self):
        for acc in ("NP_040704.1", "NP_040703.1", "NP_040708.1", "YP_512373.1"):
            self.assertTrue(is_contaminant_hit(acc, ""), acc)

    def test_real_virus_reference_is_not_flagged(self):
        self.assertFalse(is_contaminant_hit(CAMV, CAMV_TITLE))
        self.assertFalse(is_contaminant_hit("YP_009551903.1",
                                           "RdRp [plant virus]"))

    def test_title_keyword_channel(self):
        """accession 不是特征段时, 标题里的 phiX 也能命中 (第二道闸)."""
        self.assertTrue(is_contaminant_hit("XP_999999.1",
                                          "hypothetical protein [phiX174]"))
        self.assertFalse(is_contaminant_hit("XP_999999.1", "hypothetical protein"))

    def test_other_phages_are_out_of_scope(self):
        """其它噬菌体**不在**本脚本范围 —— 它们更可能是组装带的细菌序列, 归
        ICTV genome-type 的 'Phage (bacterial)' 层管, 不在这里混为一谈。"""
        self.assertFalse(is_contaminant_hit("YP_000001.1",
                                           "tail fiber [Bacillus phage SPP1]"))


class TestContaminationScan(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _write(self, rows):
        f = self.d / "raw.tsv"
        f.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
        return f

    def test_frac_is_share_of_all_hits_not_count(self):
        """C1 三条命中里 1 条 phiX -> frac=0.333; C2 两条里 2 条 -> 1.0。
        绝对条数相同(都是打中 phiX 的 contig), 只有占比能把两者分开。"""
        raw = self._write([
            raw_row("C1", PHIX, PHIX_TITLE, 500),
            raw_row("C1", CAMV, CAMV_TITLE, 400),
            raw_row("C1", CAMV, CAMV_TITLE, 300),
            raw_row("C2", PHIX, PHIX_TITLE, 500),
            raw_row("C2", PHIX, PHIX_TITLE, 450),
            raw_row("C3", CAMV, CAMV_TITLE, 500),
        ])
        out = self.d / "c.tsv"
        n_q, n_hit = contamination_scan(raw, out)
        self.assertEqual((n_q, n_hit), (3, 2))   # 3 个 query, 2 个有污染命中
        rows = dict((p[0], p) for p in
                    (l.rstrip("\n").split("\t")
                     for l in out.read_text(encoding="utf-8").splitlines()[1:]))
        self.assertEqual(rows["C1"][3], "0.333")
        self.assertEqual(rows["C2"][3], "1.000")
        self.assertNotIn("C3", rows, "干净的 contig 不该出现在表里")

    def test_top_ref_is_by_bitscore(self):
        """top_contam_ref 记录污染命中里 bitscore 最强的那条参考."""
        raw = self._write([
            raw_row("C1", "NP_040703.1", "head morphogenesis [phiX174]", 100),
            raw_row("C1", "NP_040708.1", "major spike protein [phiX174]", 900),
        ])
        out = self.d / "c.tsv"
        contamination_scan(raw, out)
        row = out.read_text(encoding="utf-8").splitlines()[1].split("\t")
        self.assertEqual(row[4], "NP_040708.1")

    def test_header_matches_writer(self):
        self.assertEqual(len(CONTAM_HDR.rstrip("\n").split("\t")), 6)
        raw = self._write([raw_row("C1", PHIX, PHIX_TITLE, 500)])
        out = self.d / "c.tsv"
        contamination_scan(raw, out)
        got = out.read_text(encoding="utf-8").splitlines()[0]
        self.assertEqual(got, CONTAM_HDR.rstrip("\n"))

    def test_deterministic(self):
        raw = self._write([raw_row("C%d" % i, PHIX, PHIX_TITLE, 100 + i)
                           for i in range(20)])
        a, b = self.d / "a.tsv", self.d / "b.tsv"
        contamination_scan(raw, a)
        contamination_scan(raw, b)
        self.assertEqual(a.read_bytes(), b.read_bytes())


class TestScanOutdirCli(unittest.TestCase):
    """端到端: 造一个假的筛查输出目录, 跑 CLI, 看产物与阈值行为."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "out"
        for g, rows in (("G_dirty", [raw_row("c1", PHIX, PHIX_TITLE, 900)]),
                        ("G_mixed", [raw_row("c2", PHIX, PHIX_TITLE, 500),
                                     raw_row("c2", CAMV, CAMV_TITLE, 400)]),
                        ("G_clean", [raw_row("c3", CAMV, CAMV_TITLE, 500)])):
            d = self.out / STAGE_DIRS["loci"] / g
            d.mkdir(parents=True)
            (d / ("%s.s1_raw.tsv" % g)).write_text(
                "\n".join(rows) + "\n", encoding="utf-8", newline="\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_cli_writes_both_products(self):
        rc = qc.main(["--outdir", str(self.out), "--frac", "0.5",
                      "--json", str(self.out / "s.json")])
        self.assertEqual(rc, 0)
        tsv = self.out / "qc_contamination.tsv"
        md = self.out / "QC_CONTAMINATION.md"
        self.assertTrue(tsv.is_file() and md.is_file())
        lines = tsv.read_text(encoding="utf-8").splitlines()
        # 表头 = genome + CONTAM_HDR + 两个 viral_supported 口径列
        self.assertEqual(lines[0], "genome\t" + CONTAM_HDR.rstrip("\n")
                         + "\tn_viral_supported\tn_contam_viral_supported")
        self.assertEqual(len(lines) - 1, 2, "应有 2 条有污染命中的 contig")
        # ours 布局拿不到位点级 verdict, 这两列留空 (不编造数字)
        for ln in lines[1:]:
            self.assertEqual(ln.split("\t")[7:9], ["", ""])
        body = md.read_text(encoding="utf-8")
        self.assertIn("G_dirty", body)      # frac=1.0 -> 建议剔除
        self.assertIn("G_mixed", body)      # frac=0.5 -> 达到阈值
        self.assertNotIn("G_clean", body.split("建议剔除清单")[1])

    def test_threshold_only_affects_the_flagged_list(self):
        """阈值只影响"建议剔除"，明细表恒含全部有污染命中的 contig。"""
        qc.main(["--outdir", str(self.out), "--frac", "0.99"])
        md = (self.out / "QC_CONTAMINATION.md").read_text(encoding="utf-8")
        flagged = md.split("建议剔除清单")[1]
        self.assertIn("G_dirty", flagged)
        self.assertNotIn("G_mixed", flagged)     # 0.5 < 0.99 不标
        tsv = (self.out / "qc_contamination.tsv").read_text(encoding="utf-8")
        self.assertIn("G_mixed", tsv)            # 但明细里还在

    def test_missing_loci_dir_is_a_loud_error(self):
        empty = Path(self.tmp.name) / "nothing"
        empty.mkdir()
        with self.assertRaises(SystemExit):
            qc.scan_outdir(empty)


if __name__ == "__main__":
    unittest.main(verbosity=2)
