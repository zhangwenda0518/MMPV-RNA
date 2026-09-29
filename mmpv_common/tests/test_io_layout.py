#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_io_layout.py — mmpv_common.io_layout 单元测试
====================================================
直接运行:  python mmpv_common/tests/test_io_layout.py
或 pytest: python -m pytest mmpv_common/tests/test_io_layout.py -q

核心不变量:
  1. legacy 布局与 v3.0 硬编码目录名逐键一致 (改名即测试失败, 防回归)
  2. standard 布局键集与 legacy 完全一致 (只许改名, 不许丢键)
  3. 布局解析优先级 CLI > 环境变量 > legacy; 非法值报错
  4. 三个构造器 (discovery/analysis/meta) 两布局键集一致
  5. 边界定位 locate_centroids / locate_taxonomy 命中两布局与三代旧名
  6. handoff 台账幂等重写
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from mmpv_common import io_layout as iol


class TestLegacyFrozen(unittest.TestCase):
    """legacy 布局冻结测试: 改任何一行都应红 (保护既有 checkpoint 兼容)。"""

    LEGACY_EXPECT = {
        "d_clean": "00a_CleanData",
        "d_hostdep": "00b_HostDepletion",
        "d_bbnorm": "00c_BBnorm",
        "d_asm": "01_Assembly",
        "d_ident": "02a_Identification",
        "d_filter": "02b_Filter",
        "d_cobra": "03a_COBRA",
        "d_merge": "03b_MergeSamples",
        "d_cluster": "04_CLUSTER",
        "d_centroids": "04_CLUSTER/4_centroids",
        "d_taxonomy": "05_Taxonomy",
        "d_host_pred": "06_HostPrediction",
        "d_checkv": "07_Checkv",
        "d_rescue": "08_Rescue",
        "d_analysis": "09a_Virome_Analysis",
        "d_verify": "09b_Analysis_Verify",
        "d_reports": "10_Reports",
        "a_detect": "01_detection",
        "a_filter": "02_filtering",
        "a_variants": "03_variants",
        "a_post": "04_post_analysis",
        "a_assembly": "05_assembly",
        "a_extract": "06_extraction",
        "a_similarity": "07_similarity",
        "a_dvg": "08_dvg",
        "a_report": "09_report",
        "m_search": "search",
        "m_info": "info",
        "m_plot": "plot",
        "m_down": "down",
        "h_genome": "host_reference/genome",
        "h_hostdb": "host_reference/hostdb",
        "p_root": "phylo_results",
        "ph_data": "data",
        "ph_phylogeny": "phylogeny",
        "ph_recomb": "recomb",
        "ph_popgen": "popgen",
        "ph_select": "select",
        "ph_time": "time",
        "ph_geography": "geography",
        "ph_report": "report",
        "e_loci": "01_Loci",
        "e_verdict": "02_Verdict",
        "e_rvdb": "03_RVDB",
        "e_summary": "04_Summary",
        "x_handoff": "90_Handoff",
    }

    def test_legacy_frozen(self):
        actual = iol.layout_dirs("legacy")
        for key, expected in self.LEGACY_EXPECT.items():
            self.assertEqual(actual.get(key), expected,
                             f"legacy 键 {key} 被改动: {actual.get(key)!r} != {expected!r}")

    def test_key_sets_identical(self):
        self.assertEqual(set(iol.layout_dirs("legacy")),
                         set(iol.layout_dirs("standard")),
                         "legacy 与 standard 键集必须一致 (只许改名, 不许丢键)")


class TestResolution(unittest.TestCase):
    def setUp(self):
        os.environ.pop(iol.LAYOUT_ENV, None)

    def tearDown(self):
        os.environ.pop(iol.LAYOUT_ENV, None)

    def test_default_legacy(self):
        self.assertEqual(iol.resolve_layout(), "legacy")

    def test_env_override(self):
        os.environ[iol.LAYOUT_ENV] = "standard"
        self.assertEqual(iol.resolve_layout(), "standard")

    def test_explicit_beats_env(self):
        os.environ[iol.LAYOUT_ENV] = "standard"
        self.assertEqual(iol.resolve_layout("legacy"), "legacy")

    def test_invalid_raises(self):
        with self.assertRaises(ValueError):
            iol.resolve_layout("nope")

    def test_normalize_writes_env(self):
        os.environ[iol.LAYOUT_ENV] = "legacy"
        self.assertEqual(iol.normalize_layout_env("standard"), "standard")
        self.assertEqual(os.environ[iol.LAYOUT_ENV], "standard")


class TestBuilders(unittest.TestCase):
    def test_discovery_legacy_matches_v30_literal(self):
        from pathlib import Path as P
        out, raw = P("/proj"), P("/reads")
        d = iol.build_discovery_dirs(out, raw, "legacy")
        self.assertEqual(d['clean'], out / '00a_CleanData')
        self.assertEqual(d['hostdep'], out / '00b_HostDepletion')
        self.assertEqual(d['centroids'], out / '04_CLUSTER' / '4_centroids')
        self.assertEqual(d['_centroids_v1'], out / '04_CLUSTER' / '04_centroids')
        self.assertEqual(d['_analysis_dirs'], [out / '09a_Virome_Analysis',
                                               out / '09_Virome_Analysis'])
        self.assertEqual(d['raw'], raw)

    def test_discovery_standard_new_names(self):
        out, raw = Path("/proj"), Path("/reads")
        d = iol.build_discovery_dirs(out, raw, "standard")
        # standard: out = ③ 自己的根 <项目>/03_Discovery, 根内独立编号
        self.assertEqual(d['clean'], out / '01_CleanData')
        self.assertEqual(d['hostdep'], out / '02_HostDepleted')
        self.assertEqual(d['asm'], out / '01_Assembly')
        self.assertEqual(d['centroids'], out / '06_CLUSTER' / 'centroids')
        self.assertEqual(d['taxonomy'], out / '07_Taxonomy')
        self.assertEqual(d['rescue_dir'], out / '10_Rescue')
        self.assertEqual(d['reports'], out / '13_Report')

    def test_analysis_dirs(self):
        out = Path("/an")
        legacy = iol.build_analysis_dirs(out, "legacy")
        std = iol.build_analysis_dirs(out, "standard")
        self.assertEqual(set(legacy), set(std))
        self.assertEqual(legacy['extract'], out / '06_extraction')
        self.assertEqual(std['extract'], out / '06_Extraction')
        self.assertEqual(legacy['detect'], out / '01_detection')
        self.assertEqual(std['detect'], out / '01_Detection')
        self.assertEqual(std['report'], out / '09_Report')

    def test_meta_dirs(self):
        w = Path("/meta")
        legacy = iol.build_meta_dirs(w, "legacy")
        std = iol.build_meta_dirs(w, "standard")
        self.assertEqual(legacy['down'], str(w / 'down'))
        self.assertEqual(std['down'], str(w / '03_RawData'))
        self.assertEqual(std['search'], str(w / '01_Search'))
        self.assertEqual(iol.dir_name('h_genome', 'standard'), '05_HostRef/genome')

    def test_phylo_modules_numbered_by_output_order(self):
        """standard: ⑤ per-virus 模块按输出顺序编号; legacy 保持原名。"""
        w = Path("/vir")
        legacy = iol.build_phylo_dirs(w, "legacy")
        std = iol.build_phylo_dirs(w, "standard")
        self.assertEqual(legacy['data'], w / 'data')
        self.assertEqual(legacy['phylogeny'], w / 'phylogeny')
        self.assertEqual([p.name for p in
                          (std[m] for m in iol.PHYLO_MODULES)],
                         ['01_data', '02_phylogeny', '03_popgen', '04_recomb',
                          '05_select', '06_time', '07_geography', '08_report'])
        self.assertEqual(std['popgen'], w / '03_popgen')
        self.assertEqual(std['recomb'], w / '04_recomb')
        self.assertEqual(std['phylogeny'], w / '02_phylogeny')
        # ph_dir 助手: 返回目录名字符串, 可直接嵌入 os.path.join
        self.assertEqual(iol.ph_dir('time', 'legacy'), 'time')
        self.assertEqual(iol.ph_dir('time', 'standard'), '06_time')
        with self.assertRaises(KeyError):
            iol.ph_dir('nope')

    def test_eve_dirs(self):
        out = Path("/eve")
        legacy = iol.build_eve_dirs(out, "legacy")
        std = iol.build_eve_dirs(out, "standard")
        # EVE 内部编号两布局一致 (本就符合编号惯例); standard 输出根约定 06_EVE
        self.assertEqual(legacy['loci'], out / '01_Loci')
        self.assertEqual(std['summary'], out / '04_Summary')
        self.assertEqual(set(legacy), set(std))


class TestLocateAndManifest(unittest.TestCase):
    def _mk(self, root: Path, centroids: Path, taxonomy: Path):
        centroids.parent.mkdir(parents=True, exist_ok=True)
        taxonomy.parent.mkdir(parents=True, exist_ok=True)
        centroids.write_text(">a\nACGT\n", encoding="utf-8")
        taxonomy.write_text("Accession\n", encoding="utf-8")

    def test_locate_legacy(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._mk(root,
                     root / "04_CLUSTER/4_centroids/final_centroids.fasta",
                     root / "05_Taxonomy/integrated/final_integrated_classification.tsv")
            self.assertEqual(iol.locate_centroids(root),
                             root / "04_CLUSTER/4_centroids/final_centroids.fasta")
            tax = iol.locate_taxonomy(root)
            self.assertIsInstance(tax, list) and self.assertTrue(tax)
            self.assertIn("integrated", str(tax[0]))

    def test_locate_legacy_old_centroid_names(self):
        for old in ("04_centroids", "4.centroids"):
            with tempfile.TemporaryDirectory() as td:
                root = Path(td)
                self._mk(root,
                         root / "04_CLUSTER" / old / "final_centroids.fasta",
                         root / "05_Taxonomy/integrated/final_integrated_classification.tsv")
                self.assertIsNotNone(iol.locate_centroids(root), old)

    def test_locate_standard(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._mk(root,
                     root / "06_CLUSTER/centroids/final_centroids.fasta",
                     root / "07_Taxonomy/integrated/final_integrated_classification.tsv")
            self.assertEqual(iol.locate_centroids(root),
                             root / "06_CLUSTER/centroids/final_centroids.fasta")
            tax = iol.locate_taxonomy(root)
            self.assertIsInstance(tax, list) and self.assertTrue(tax)
            self.assertEqual(tax[0], root / "07_Taxonomy/integrated/final_integrated_classification.tsv")

    def test_locate_taxonomy_per_sample_files(self):
        """taxonomy 真实形态 = 逐样本 {sample}.integrated/; locate_taxonomy 返回全部"""
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            tax = root / "07_Taxonomy"
            for s in ("S1", "S2"):
                f = tax / f"{s}.integrated" / "final_integrated_classification.tsv"
                f.parent.mkdir(parents=True, exist_ok=True)
                f.write_text("contig_id\n", encoding="utf-8")
            got = iol.locate_taxonomy(root)
            assert isinstance(got, list) and len(got) == 2, got

    def test_manifest_idempotent(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            m1 = iol.write_boundary(root, "discovery->analysis",
                                    {"centroids": "/a.fasta"}, "t1")
            m2 = iol.write_boundary(root, "discovery->analysis",
                                    {"centroids": "/b.fasta", "taxonomy": "/t.tsv"}, "t2")
            self.assertEqual(m1, m2)
            lines = m2.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(lines[0].split("\t"),
                             ["boundary", "key", "path", "produced_by", "produced_at"])
            data = {tuple(l.split("\t")[:2]): l.split("\t") for l in lines[1:]}
            self.assertEqual(len(data), 2)                      # 幂等: 无重复行
            self.assertEqual(data[("discovery->analysis", "centroids")][2], "/b.fasta")

    def test_manifest_lands_at_project_root(self):
        """传入管线独立根时, 台账应上提到项目根的 90_Handoff/。"""
        with tempfile.TemporaryDirectory() as td:
            proj = Path(td)
            droot = proj / "03_Discovery"
            droot.mkdir()
            m = iol.write_boundary(droot, "discovery->analysis",
                                   {"centroids": "/a.fasta"}, "t")
            self.assertEqual(m, proj / "90_Handoff" / "handoff_manifest.tsv")
            self.assertTrue(m.is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
