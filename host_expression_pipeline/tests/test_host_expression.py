#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_host_expression.py — ⑧ 宿主表达模块单元测试
=====================================================
直接运行:  python host_expression_pipeline/tests/test_host_expression.py
或 pytest: python -m pytest host_expression_pipeline/tests/test_host_expression.py -q

覆盖纯逻辑核心 (不需要 salmon/gffread/DESeq2 等外部工具):
  1. derive_virus_labels: 列名解析 / 聚合 / 阳性判定 / 宽表 / 样本名单
  2. quant_qc: meta_info 解析与双阈值判定 (含缺失文件)
  3. viral_load_scan: BH-FDR 已知值 + Spearman 与 scipy 逐值一致
  4. host_expression 编排器: clean reads / 基因组 / GFF / Core14 / 检出表 定位
  5. host_expression_report: 合成产物 HTML 生成 + 缺段降级
"""

import json
import sys
import tempfile
import unittest
from pathlib import Path

PIPELINE_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = PIPELINE_DIR.parent
sys.path.insert(0, str(PIPELINE_DIR))
sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pandas as pd

import derive_virus_labels as dvl
import host_expression_report as her
import host_expression as hexp
import quant_qc
import viral_load_scan as vls


class TestDeriveLabels(unittest.TestCase):
    def _df(self):
        return pd.DataFrame({
            "Sample": ["S1", "S1", "S1", "S1", "S2", "S2"],
            "Adjusted_Species": ["GCVA", "GCVA", "PSTVd", "PSTVd", "GCVA", "PSTVd"],
            "Rep_Coverage(%)": [95.0, 100.0, 50.0, 40.0, 10.0, 99.0],
            "Rep_MeanDepth": [5.0, 8.0, 1.0, 1.5, 0.2, 3.0],
            "Asm_TPM": [12.0, 15.0, 2.0, 3.0, 0.5, 7.0],
            "Rep_Reads": [100, 150, 20, 30, 5, 60],
        })

    def test_resolve_columns(self):
        cols = dvl.resolve_columns(self._df())
        self.assertEqual(cols["sample"], "Sample")
        self.assertEqual(cols["virus"], "Adjusted_Species")
        self.assertEqual(cols["cov"], "Rep_Coverage(%)")
        self.assertEqual(cols["load"], "Asm_TPM")

    def test_resolve_columns_override_and_missing(self):
        df = self._df().rename(columns={"Rep_Coverage(%)": "mycov"})
        cols = dvl.resolve_columns(df, {"cov": "mycov"})
        self.assertEqual(cols["cov"], "mycov")
        with self.assertRaises(KeyError):
            dvl.resolve_columns(pd.DataFrame({"foo": [1]}))   # 缺 Sample 等

    def test_build_virus_labels_aggregation_and_threshold(self):
        vlab = dvl.build_virus_labels(self._df(), dvl.resolve_columns(self._df()),
                                      cov_threshold=90.0)
        # (S1, GCVA): cov max=100, depth max=8, load max=15, reads sum=250, 阳性
        row = vlab[(vlab["Sample"] == "S1") & (vlab["Virus"] == "GCVA")].iloc[0]
        self.assertEqual(row["Covered_pct"], 100.0)
        self.assertEqual(row["MeanDepth"], 8.0)
        self.assertEqual(row["Load"], 15.0)
        self.assertEqual(row["Reads"], 250)
        self.assertEqual(int(row["Positive"]), 1)
        # (S1, PSTVd): cov max=50 < 90 → 阴性; (S2, GCVA): 10 → 阴性
        self.assertEqual(int(vlab[(vlab["Sample"] == "S1")
                                  & (vlab["Virus"] == "PSTVd")]["Positive"].iloc[0]), 0)
        self.assertEqual(int(vlab[(vlab["Sample"] == "S2")
                                  & (vlab["Virus"] == "GCVA")]["Positive"].iloc[0]), 0)
        # (S2, PSTVd): 99 → 阳性
        self.assertEqual(int(vlab[(vlab["Sample"] == "S2")
                                  & (vlab["Virus"] == "PSTVd")]["Positive"].iloc[0]), 1)

    def test_build_virus_labels_depth_threshold(self):
        vlab = dvl.build_virus_labels(self._df(), dvl.resolve_columns(self._df()),
                                      cov_threshold=90.0, depth_threshold=4.0)
        # (S2, PSTVd): cov 99 但 depth 3 < 4 → 阴性
        self.assertEqual(int(vlab[(vlab["Sample"] == "S2")
                                  & (vlab["Virus"] == "PSTVd")]["Positive"].iloc[0]), 0)

    def test_sample_labels_wide_and_undetected_samples(self):
        vlab = dvl.build_virus_labels(self._df(), dvl.resolve_columns(self._df()))
        slab = dvl.build_sample_labels(vlab, samples=["S1", "S2", "S3"])
        self.assertEqual(len(slab), 3)
        s1 = slab[slab["Sample"] == "S1"].iloc[0]
        self.assertEqual(s1["Virus_any"], "+")
        self.assertEqual(int(s1["N_virus"]), 1)
        self.assertEqual(int(s1["GCVA_positive"]), 1)
        self.assertEqual(int(s1["PSTVd_positive"]), 0)
        self.assertAlmostEqual(float(s1["PSTVd_load"]), 3.0)
        # S3 未出现在检出表中 → 全阴性
        s3 = slab[slab["Sample"] == "S3"].iloc[0]
        self.assertEqual(s3["Virus_any"], "-")
        self.assertEqual(int(s3["N_virus"]), 0)

    def test_sanitize_virus(self):
        self.assertEqual(dvl.sanitize_virus("Cytorhabdovirus sp. 'lycii'"),
                         "Cytorhabdovirus_sp._lycii")
        self.assertTrue(dvl.sanitize_virus("   "))          # 空名兜底

    def test_read_sample_list(self):
        with tempfile.TemporaryDirectory() as td:
            f = Path(td) / "qc_summary.tsv"
            f.write_text("sample\tpercent_mapped\nS1\t88\nS2\t70\n",
                         encoding="utf-8")
            self.assertEqual(dvl.read_sample_list(f), ["S1", "S2"])
            f2 = Path(td) / "plain.txt"
            f2.write_text("S1\nS2\n\n", encoding="utf-8")
            self.assertEqual(dvl.read_sample_list(f2), ["S1", "S2"])


class TestQuantQC(unittest.TestCase):
    def _meta(self, sample_dir: Path, processed: int, mapped: int):
        aux = sample_dir / "aux_info"
        aux.mkdir(parents=True, exist_ok=True)
        (aux / "meta_info.json").write_text(json.dumps({
            "num_processed": processed, "num_mapped": mapped,
            "percent_mapped": 100.0 * mapped / processed,
        }), encoding="utf-8")

    def test_collect_qc_thresholds_and_missing(self):
        with tempfile.TemporaryDirectory() as td:
            quant = Path(td)
            (quant / "quant").mkdir(parents=True)
            self._meta(quant / "quant" / "good", 2_000_000, 1_400_000)
            self._meta(quant / "quant" / "lowmap", 2_000_000, 4_000)      # 0.2%
            self._meta(quant / "quant" / "fewreads", 100_000, 90_000)     # 9e4 < 1e6
            (quant / "quant" / "broken").mkdir()                          # 无 meta_info
            rows = quant_qc.collect_qc(quant, min_mapped=10.0,
                                       min_processed=1e6)
            by = {r["sample"]: r for r in rows}
            self.assertTrue(by["good"]["pass"])
            self.assertFalse(by["lowmap"]["pass"])
            self.assertEqual(by["lowmap"]["note"], "low_mapped")
            self.assertFalse(by["fewreads"]["pass"])
            self.assertEqual(by["fewreads"]["note"], "low_reads")
            self.assertFalse(by["broken"]["pass"])
            self.assertEqual(by["broken"]["note"], "no_meta_info")

    def test_write_qc_outputs(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "qc"
            rows = [{"sample": "A", "num_processed": 10, "num_mapped": 9,
                     "percent_mapped": 90.0, "pass": True, "note": ""},
                    {"sample": "B", "num_processed": 10, "num_mapped": 1,
                     "percent_mapped": 10.0, "pass": False, "note": "low_mapped"}]
            paths = quant_qc.write_qc(rows, out)
            self.assertEqual(paths["passed"].read_text(encoding="utf-8").strip(), "A")
            self.assertTrue(paths["summary"].is_file())


class TestLoadScanMath(unittest.TestCase):
    def test_bh_fdr_known_values(self):
        p = np.array([0.01, 0.04, 0.03, 0.002])
        q = vls.bh_fdr(p)
        np.testing.assert_allclose(q, [0.02, 0.04, 0.04, 0.008], atol=1e-12)

    def test_bh_fdr_monotone_and_nan_safe(self):
        p = np.array([0.001, 0.01, 0.02, 0.2, np.nan])
        q = vls.bh_fdr(p)
        self.assertTrue(np.isnan(q[-1]))
        valid = q[:-1]
        self.assertTrue(np.all(valid >= p[:-1]))
        self.assertTrue(np.all(np.diff(valid) >= 0))

    def test_spearman_matches_scipy(self):
        from scipy.stats import spearmanr
        rng = np.random.default_rng(42)
        load = rng.normal(size=30)
        X = rng.normal(size=(5, 30))
        X[3] = 1.0                      # 零方差行 → NaN
        X[2] = X[1] + 0.01 * rng.normal(size=30)   # 近似共线行
        rho, p = vls.spearman_scan(load, X)
        for i in range(5):
            if np.isnan(rho[i]):
                continue
            ref = spearmanr(load, X[i])
            self.assertAlmostEqual(rho[i], ref.statistic, places=10)
            self.assertAlmostEqual(p[i], ref.pvalue, places=10)
        self.assertTrue(np.isnan(rho[3]))


class TestOrchestratorLocate(unittest.TestCase):
    def test_locate_standard_tree(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            (p / "02_Preprocessing/01_CleanData/3.clumpify").mkdir(parents=True)
            (p / "02_Preprocessing/01_CleanData/3.clumpify/S1_1.fa.gz").write_text("")
            (p / "04_Analysis/01_Detection/summary").mkdir(parents=True)
            (p / "04_Analysis/01_Detection/summary/all_viruses.best.summary.tsv"
             ).write_text("Sample\n")
            (p / "01_PublicData/05_HostRef/genome").mkdir(parents=True)
            (p / "01_PublicData/05_HostRef/genome/all.genome.uniq.fasta"
             ).write_text(">c\nACGT\n")
            (p / "01_PublicData/05_HostRef/genome/host.gff3").write_text("##gff\n")
            (p / "01_PublicData/02_Unified").mkdir(parents=True)
            (p / "01_PublicData/02_Unified/Global_Unified_Metadata_Core14.csv"
             ).write_text("Run\n")

            self.assertEqual(hexp.locate_clean_reads(p).name, "3.clumpify")
            self.assertEqual(hexp.locate_host_genome(p).name, "all.genome.uniq.fasta")
            self.assertEqual(hexp.locate_host_gff(p).name, "host.gff3")
            self.assertEqual(hexp.locate_metadata_core14(p).name,
                             "Global_Unified_Metadata_Core14.csv")
            self.assertEqual(hexp.locate_detect_table(p).name,
                             "all_viruses.best.summary.tsv")

    def test_locate_legacy_tree(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            (p / "00a_CleanData/2.fasta").mkdir(parents=True)
            (p / "00a_CleanData/2.fasta/S1_1.fa.gz").write_text("")
            self.assertEqual(hexp.locate_clean_reads(p).name, "2.fasta")


class TestSalmonEngineSuffixes(unittest.TestCase):
    """抽取引擎必须识别 ② clean data 的 FASTA 命名 ({sample}_1.fa.gz)。"""

    def test_detect_fasta_pe_samples(self):
        import salmon_auto_pipeline as sap
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            for name in ("S1_1.fa.gz", "S1_2.fa.gz", "S2_1.fa.gz", "S2_2.fa.gz"):
                (d / name).write_text(">r\nACGT\n", encoding="utf-8")
            samples = sap.detect_samples(d)
            self.assertEqual({s.name for s in samples}, {"S1", "S2"})
            for s in samples:
                self.assertEqual(s.layout, "PE")
                self.assertEqual(len(s.read1), 1)
                self.assertEqual(len(s.read2), 1)

    def test_detect_single_end_fasta(self):
        import salmon_auto_pipeline as sap
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / "S3.fa.gz").write_text(">r\nACGT\n", encoding="utf-8")
            samples = sap.detect_samples(d)
            self.assertEqual([s.name for s in samples], ["S3"])
            self.assertEqual(samples[0].layout, "SE")


class TestBuildCoexpMatrices(unittest.TestCase):
    """切矩阵 + 闸门 (ruthenicum 失衡教训的工程化)。"""

    def _setup(self, root: Path, n_pos=35, n_neg=35, n_genes=10):
        import build_coexp_matrices as bcm
        samples = [f"P{i:02d}" for i in range(n_pos)] + \
                  [f"N{i:02d}" for i in range(n_neg)]
        tpm = pd.DataFrame(
            np.arange(n_genes * len(samples), dtype=float).reshape(n_genes, -1),
            index=[f"g{i}" for i in range(n_genes)], columns=samples)
        tpm_path = root / "TPM_matrix.tsv"
        tpm.to_csv(tpm_path, sep="\t")
        lab = pd.DataFrame({"Sample": samples,
                            "Virus_any": ["+"] * n_pos + ["-"] * n_neg})
        lab_path = root / "sample_labels.tsv"
        lab.to_csv(lab_path, sep="\t", index=False)
        return bcm, tpm_path, lab_path

    def test_split_ok(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            bcm, tpm_path, lab_path = self._setup(root)
            out = root / "matrices"
            bcm.main_testable(["--tpm", str(tpm_path), "--labels", str(lab_path),
                               "--output_dir", str(out)])
            pos = pd.read_csv(out / "positive.tsv", sep="\t", index_col=0)
            neg = pd.read_csv(out / "negative.tsv", sep="\t", index_col=0)
            self.assertEqual(pos.shape, (10, 35))
            self.assertEqual(neg.shape, (10, 35))
            self.assertEqual(list(pos.columns)[:2], ["P00", "P01"])

    def test_gate_small_group(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            bcm, tpm_path, lab_path = self._setup(root, n_pos=5, n_neg=35)
            with self.assertRaises(SystemExit) as ctx:
                bcm.main_testable(["--tpm", str(tpm_path), "--labels", str(lab_path),
                                   "--output_dir", str(root / "matrices")])
            self.assertEqual(ctx.exception.code, 2)

    def test_gate_imbalance_skip_mode(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            bcm, tpm_path, lab_path = self._setup(root, n_pos=40, n_neg=35)
            out = root / "matrices"
            with self.assertRaises(SystemExit) as ctx:   # skip 模式: exit(0)
                bcm.main_testable(["--tpm", str(tpm_path), "--labels", str(lab_path),
                                   "--output_dir", str(out),
                                   "--max-imbalance", "1.0",
                                   "--on-gate-fail", "skip"])
            self.assertEqual(ctx.exception.code, 0)
            summary = pd.read_csv(out / "group_summary.tsv", sep="\t")
            self.assertEqual(summary.iloc[0]["status"], "SKIPPED")

    def _write_gcn(self, root: Path, name: str, genes: list, zmap):
        """zmap: {gene: {target: z}}; 未给的填 0.5; 自边填 9.0 (模拟 TEA-GCN)。"""
        d = root / name
        d.mkdir(parents=True)
        for g in genes:
            zvals = [9.0 if t == g else zmap.get(g, {}).get(t, 0.5) for t in genes]
            pd.DataFrame({"Target": genes, "Co-exp_Str": [0.8] * len(genes),
                          "Co-exp_Str_MR": [1.0] * len(genes),
                          "zScore(Co-exp_Str_MR)": zvals}).to_csv(
                d / g, sep="\t", index=False)
        return d

    def test_delta_z_and_self_edge(self):
        import diff_coexp as dc
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            rng = np.random.default_rng(7)
            n_bg = 60
            genes = [f"b{i:02d}" for i in range(n_bg)] + ["g1", "g2", "g3", "g4"]
            zp = {g: {t: round(float(rng.normal(0, 0.2)), 3)
                      for t in genes if t != g} for g in genes}
            zn = {g: {t: round(float(rng.normal(0, 0.2)), 3)
                      for t in genes if t != g} for g in genes}
            # 植入: g1-g2 增强, g3-g4 减弱 (背景 ~N(0,0.2), 极端值稀有)
            for a, b in (("g1", "g2"), ("g2", "g1")):
                zp[a][b] = 8.0; zn[a][b] = 0.0
            for a, b in (("g3", "g4"), ("g4", "g3")):
                zp[a][b] = -8.0; zn[a][b] = 0.0
            dpos = self._write_gcn(root, "pos", genes, zp)
            dneg = self._write_gcn(root, "neg", genes, zn)
            out = root / "out"
            dc.main_testable(["--pos-gcn", str(dpos), "--neg-gcn", str(dneg),
                              "--output_dir", str(out), "--top-n", "5"])
            sig = pd.read_csv(out / "diff_coexp_significant.tsv", sep="\t")
            self.assertTrue((sig["Gene"] != sig["Target"]).all())   # 自边已剔除
            m = ((sig["Gene"] == "g1") & (sig["Target"] == "g2") |
                 (sig["Gene"] == "g2") & (sig["Target"] == "g1") |
                 (sig["Gene"] == "g3") & (sig["Target"] == "g4") |
                 (sig["Gene"] == "g4") & (sig["Target"] == "g3"))
            planted = sig[m]
            self.assertEqual(len(planted), 4)                        # 4 条植入有向行
            self.assertEqual(set(planted["Delta_z"].round(6)), {8.0, -8.0})
            summary = pd.read_csv(out / "diff_coexp_summary.tsv", sep="\t")
            self.assertGreaterEqual(int(summary.iloc[0]["n_unique_pairs"]), 2)

    def test_row_order_mismatch_aligned_by_name(self):
        """行集不同/行序不同 → 按 Target 名对齐, 植入信号仍正确检出。"""
        import diff_coexp as dc
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            rng = np.random.default_rng(7)
            n_bg = 60
            genes = [f"b{i:02d}" for i in range(n_bg)] + ["g1", "g2"]
            zp = {g: {t: round(float(rng.normal(0, 0.2)), 3)
                      for t in genes if t != g} for g in genes}
            zn = {g: {t: round(float(rng.normal(0, 0.2)), 3)
                      for t in genes if t != g} for g in genes}
            for a, b in (("g1", "g2"), ("g2", "g1")):
                zp[a][b] = 8.0
                zn[a][b] = 0.0
            dpos = self._write_gcn(root, "pos", genes, zp)
            dneg = self._write_gcn(root, "neg", list(reversed(genes)), zn)
            out = root / "out"
            dc.main_testable(["--pos-gcn", str(dpos), "--neg-gcn", str(dneg),
                              "--output_dir", str(out), "--top-n", "5"])
            sig = pd.read_csv(out / "diff_coexp_significant.tsv", sep="\t")
            m = ((sig["Gene"] == "g1") & (sig["Target"] == "g2") |
                 (sig["Gene"] == "g2") & (sig["Target"] == "g1"))
            self.assertEqual(int(sig[m]["Delta_z"].round(6).eq(8.0).sum()), 2)

    def test_row_subset_aligned_by_name(self):
        """行数不同 (单侧 NaN 过滤) → 交集对齐, 单侧缺失边不计入。"""
        import diff_coexp as dc
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            genes = ["g1", "g2", "g3", "g4", "g5"]
            zp = {g: {t: 0.1 for t in genes if t != g} for g in genes}
            zn = {g: {t: 0.1 for t in genes if t != g} for g in genes}
            zp["g1"]["g2"] = 6.0
            zn["g1"]["g2"] = 0.0
            zn["g1"].pop("g3")                       # neg 缺 g3 行 (NaN 过滤)
            dpos = self._write_gcn(root, "pos", genes, zp)
            dneg = self._write_gcn(root, "neg", genes, zn)
            # neg 的 g1 文件删掉 g3 行, 其余文件保持全行
            p = dneg / "g1"
            rows = p.read_text(encoding="utf-8").splitlines()
            keep = [r for r in rows if not r.startswith("g3\t") or
                    rows.index(r) == 0]
            p.write_text("\n".join(keep) + "\n", encoding="utf-8")
            out = root / "out"
            dc.main_testable(["--pos-gcn", str(dpos), "--neg-gcn", str(dneg),
                              "--output_dir", str(out), "--top-n", "5"])
            sig = pd.read_csv(out / "diff_coexp_significant.tsv", sep="\t")
            m = (sig["Gene"] == "g1") & (sig["Target"] == "g2")
            self.assertEqual(float(sig[m]["Delta_z"].iloc[0]), 6.0)


class TestGenericGrouping(unittest.TestCase):
    """通用性: 任意二分组列/水平 (非病毒场景, 如 Condition: control/treated)。"""

    def test_bcm_custom_column_and_levels(self):
        import build_coexp_matrices as bcm
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            n_c, n_t = 12, 12
            samples = [f"C{i:02d}" for i in range(n_c)] + \
                      [f"T{i:02d}" for i in range(n_t)]
            tpm = pd.DataFrame(np.random.default_rng(3).normal(
                50, 5, size=(8, len(samples))),
                index=[f"g{i}" for i in range(8)], columns=samples)
            tpm_path = root / "TPM_matrix.tsv"
            tpm.to_csv(tpm_path, sep="\t")
            lab = pd.DataFrame({"Sample": samples,
                                "Condition": ["control"] * n_c + ["treated"] * n_t})
            lab_path = root / "labels.tsv"
            lab.to_csv(lab_path, sep="\t", index=False)
            out = root / "matrices"
            bcm.main_testable([
                "--tpm", str(tpm_path), "--labels", str(lab_path),
                "--output_dir", str(out),
                "--group-col", "Condition",
                "--pos-level", "treated", "--neg-level", "control",
                "--min-group", "5"])
            pos = pd.read_csv(out / "positive.tsv", sep="\t", index_col=0)
            self.assertEqual(list(pos.columns), [f"T{i:02d}" for i in range(n_t)])
            summary = pd.read_csv(out / "group_summary.tsv", sep="\t")
            self.assertEqual(summary.iloc[0]["group_col"], "Condition")

    def test_loadscan_strata_col(self):
        """组内扫描分层列可换 (非 Tissue, 如 Batch)。"""
        import viral_load_scan as vls
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            rng = np.random.default_rng(5)
            samples = [f"S{i:02d}" for i in range(12)]
            genes = [f"g{i}" for i in range(30)]
            tpm = pd.DataFrame(np.abs(rng.normal(50, 10, size=(30, 12))),
                               index=genes, columns=samples)
            tpm_path = root / "TPM_matrix.tsv"
            tpm.to_csv(tpm_path, sep="\t")
            # 6 阳性 (带载量), 分属 batch a/b
            rows = []
            for i, s in enumerate(samples):
                pos = i % 2 == 0            # 阳性交错分布: 每 batch 各 3 阳
                rows.append({"Sample": s, "Virus": "VX",
                             "Load": 100.0 + i if pos else 0.0,
                             "Positive": 1 if pos else 0})
            vlab_path = root / "virus_labels.tsv"
            pd.DataFrame(rows).to_csv(vlab_path, sep="\t", index=False)
            meta_path = root / "sample_labels.tsv"
            pd.DataFrame({"Sample": samples,
                          "Batch": ["a"] * 6 + ["b"] * 6}).to_csv(
                meta_path, sep="\t", index=False)
            out = root / "scan"
            vls.main_testable(["--tpm", str(tpm_path), "--labels", str(vlab_path),
                               "--meta", str(meta_path), "--output_dir", str(out),
                               "--strata-col", "Batch",
                               "--min-samples", "4", "--min-positives", "2"])
            summary = pd.read_csv(out / "load_scan_summary.tsv", sep="\t")
            self.assertEqual(set(summary["stratum"]), {"all", "a", "b"})
            self.assertTrue((out / "VX__a_genome_scan.tsv").is_file())
            self.assertTrue((out / "VX__b_genome_scan.tsv").is_file())


class TestModulesSubstepWiring(unittest.TestCase):
    """编排器 _run_modules_substep 真实调用: 网络定位/断点/产物。"""

    def _project(self, root: Path):
        import shutil
        cdir = root / "05_Coexp"
        dpos, dneg = _planted_networks(root / "_net")
        for grp, src in (("teagcn_positive", dpos), ("teagcn_negative", dneg)):
            dst = cdir / grp / "Completed_GCN" / "netA"
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
        (root / "06_DE").mkdir(exist_ok=True)
        pd.DataFrame({"log2FoldChange": [2.0] * 6, "padj": [0.01] * 6},
                     index=["g1", "g2", "g3", "g4", "g5", "g6"]).to_csv(
            root / "06_DE" / "deseq2_degs.tsv", sep="\t", index_label="gene")
        return cdir

    def test_run_and_checkpoint(self):
        import argparse
        import logging
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            cdir = self._project(root)
            args = argparse.Namespace(force=False, module_deg=None,
                                      module_z=1.7, module_resolution=1.0)
            dirs = {"coexp": cdir, "de": root / "06_DE",
                    "labels": root / "04_Labels", "quant": root / "02_Quant"}
            log = logging.getLogger("test_wiring")
            log.addHandler(logging.NullHandler())
            hexp._run_modules_substep(args, dirs, log)
            s = cdir / "modules" / "module_summary.tsv"
            self.assertTrue(s.is_file())
            summ = pd.read_csv(s, sep="\t")
            self.assertEqual(set(summ["Group"]), {"positive"})  # 阴性无模块
            ck = cdir / ".checkpoints" / "coexp_modules.ok"
            self.assertTrue(ck.is_file())
            s.unlink()
            hexp._run_modules_substep(args, dirs, log)   # 断点在但产物缺 → 重跑
            self.assertTrue(s.is_file())


class TestCoexpVirusSpec(unittest.TestCase):
    """--virus 指定病毒: coexp 分组切 <病毒>_positive 1 vs 0 (默认不变)。"""

    @staticmethod
    def _labels(root: Path) -> Path:
        lab = root / "sample_labels.tsv"
        lab.write_text(
            "Sample\tVirus_any\tGCVA_positive\tPSTVd_positive\n"
            "S1\t+\t1\t0\nS2\t-\t0\t1\n", encoding="utf-8")
        return lab

    @staticmethod
    def _args(virus):
        import argparse
        return argparse.Namespace(virus=virus, group_col="Virus_any",
                                  pos_level="+", neg_level="-")

    def test_default_unchanged_and_virus_switch(self):
        with tempfile.TemporaryDirectory() as td:
            lab = self._labels(Path(td))
            self.assertEqual(hexp._coexp_group_spec(self._args(None), lab),
                             ("Virus_any", "+", "-"))
            self.assertEqual(hexp._coexp_group_spec(self._args("GCVA"), lab),
                             ("GCVA_positive", "1", "0"))

    def test_unknown_virus_lists_available(self):
        with tempfile.TemporaryDirectory() as td:
            lab = self._labels(Path(td))
            with self.assertRaises(FileNotFoundError) as ctx:
                hexp._coexp_group_spec(self._args("NOPE"), lab)
            self.assertIn("GCVA_positive", str(ctx.exception))
            self.assertIn("PSTVd_positive", str(ctx.exception))


class TestQCFilterWiring(unittest.TestCase):
    """QC 名单 (--samples) 在 LoadScan 与模块×载量关联中真实生效。"""

    def test_loadscan_samples_filter(self):
        import viral_load_scan as vls
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            n = 10
            load = np.arange(1, n + 1, dtype=float)
            samples = [f"S{i:02d}" for i in range(1, n + 1)]
            rng = np.random.default_rng(5)
            tpm = pd.DataFrame({s: [10.0, 20.0] for s in samples},
                               index=["g1", "g2"])
            tpm.loc["g1"] = load * 10 + rng.normal(0, 0.5, n)   # 植入正相关
            tpm.to_csv(root / "TPM.tsv", sep="\t")
            pd.DataFrame({"Sample": samples, "Virus": "GCVA",
                          "Positive": 1, "Load": load}).to_csv(
                root / "virus_labels.tsv", sep="\t", index=False)
            (root / "qc_passed.txt").write_text(
                "\n".join(samples[:6]) + "\n", encoding="utf-8")
            out = root / "scan"
            vls.main_testable(["--tpm", str(root / "TPM.tsv"),
                               "--labels", str(root / "virus_labels.tsv"),
                               "--samples", str(root / "qc_passed.txt"),
                               "--min-samples", "4", "--min-positives", "2",
                               "--output_dir", str(out)])
            summ = pd.read_csv(out / "load_scan_summary.tsv", sep="\t")
            self.assertEqual(int(summ.iloc[0]["n_samples"]), 6)  # QC 过滤生效

    def test_modules_assoc_samples_filter(self):
        import coexp_modules as cm
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            dpos, _ = _planted_networks(root)
            pd.DataFrame({"log2FoldChange": [2.0] * 6, "padj": [0.01] * 6},
                         index=["g1", "g2", "g3", "g4", "g5", "g6"]).to_csv(
                root / "degs.tsv", sep="\t", index_label="gene")
            n_bg, n = 40, 12
            load = np.arange(1, n + 1, dtype=float)
            samples = [f"S{i:02d}" for i in range(1, n + 1)]
            tpm = pd.DataFrame({s: [10.0] * (n_bg + 6) for s in samples},
                               index=[f"b{i:02d}" for i in range(n_bg)] +
                                     ["g1", "g2", "g3", "g4", "g5", "g6"])
            for k, g in enumerate(["g1", "g2", "g3"]):
                rng = np.random.default_rng(11 + k)
                tpm.loc[g] = load * 10 + rng.normal(0, 0.5, n)
            tpm.to_csv(root / "TPM.tsv", sep="\t")
            pd.DataFrame({"Sample": samples, "Virus": "GCVA",
                          "Positive": 1, "Load": load}).to_csv(
                root / "virus_labels.tsv", sep="\t", index=False)
            (root / "qc_passed.txt").write_text(
                "\n".join(samples[:8]) + "\n", encoding="utf-8")
            out = root / "modules"
            cm.main_testable(["--gcn-dir", f"positive={dpos}",
                              "--deg", str(root / "degs.tsv"),
                              "--tpm", str(root / "TPM.tsv"),
                              "--virus-labels", str(root / "virus_labels.tsv"),
                              "--samples", str(root / "qc_passed.txt"),
                              "--output_dir", str(out)])
            assoc = pd.read_csv(out / "module_load_assoc.tsv", sep="\t")
            self.assertEqual(set(assoc["N_Samples"]), {8})      # QC 过滤生效
            self.assertGreater(float(assoc["Rho"].max()), 0.9)  # 信号保留

    def test_coexp_dryrun_tolerates_mp_on_windows(self):
        """dry-run 下 mp 入口 fork 不可用只警告, 不中断预览。"""
        import argparse
        import logging
        import os as _os
        with tempfile.TemporaryDirectory() as td:
            args = argparse.Namespace(
                dry_run=True, stage=["all"], virus=None,
                group_col="Virus_any", pos_level="+", neg_level="-",
                coexp_min_group=30, coexp_max_imbalance=3.0,
                coexp_ks=5, coexp_ke=30, teagcn_entry="mp",
                teagcn_python=None, coexp_workers=4, coexp_qc_gate="warn",
                coexp_top_n=20, module_z=1.7, module_resolution=1.0,
                rank_workers=8)
            dirs = {"coexp": Path(td) / "05_Coexp"}
            log = logging.getLogger("test_dryrun_mp")
            log.addHandler(logging.NullHandler())
            hexp.stage_coexp(args, dirs, log)                   # 不应 raise
            if not hasattr(_os, "fork"):                        # 非 dry-run 才硬失败
                args2 = argparse.Namespace(**{**vars(args), "dry_run": False})
                with self.assertRaises(EnvironmentError):
                    hexp.stage_coexp(args2, dirs, log)


class TestEnrichAnalysis(unittest.TestCase):
    """07_Enrich: 注释解析 / ORA 统计 / 基因集收集。"""

    @staticmethod
    def _synthetic(root: Path):
        bg = [f"g{i:03d}" for i in range(300)]
        t1, t2 = set(bg[:30]), set(bg[30:60])
        # 2026-10-03: ORA 统一到可注释空间 (N=背景∩注释=60, n=查询∩注释=10) 后,
        # 旧植入量 (t1 8/t2 2) 的 padj≈0.06 不再显著 — 提高 t1 植入量到 9,
        # 保持"植入显著 / decoy 不显著"的原意 (无注释基因仍留 2 个以覆盖本修复)。
        de = list(t1)[:9] + list(t2)[:1] + ["g290", "g291"]
        (root / "bg.txt").write_text("\n".join(bg) + "\n", encoding="utf-8")
        (root / "DE.txt").write_text("\n".join(de) + "\n", encoding="utf-8")
        (root / "tiny.txt").write_text("g001\n", encoding="utf-8")
        rows = ["## eggNOG 5.0.0", "### rules",
                "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\t"
                "max_annot_lvl\tCOG_category\tDescription\tPreferred_name\t"
                "GOs\tEC\tKEGG_ko"]
        for g in bg:
            gos = "GO:0001111" if g in t1 else \
                  ("GO:0002222" if g in t2 else "")
            ko = "ko:K00001" if g in t1 else ""
            rows.append(f"{g}\t-\t0\t0\t-\t-\t-\t-\t-\t{gos}\t-\t{ko}")
        (root / "anno.emapper.annotations").write_text(
            "\n".join(rows) + "\n", encoding="utf-8")
        return t1, t2

    def test_emapper_ora_planted_and_decoy(self):
        import enrich_analysis as ea
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._synthetic(root)
            out = root / "out"
            ea.main_testable(["--annotation", str(root / "anno.emapper.annotations"),
                              "--gene-set", f"DE={root / 'DE.txt'}",
                              "--gene-set", f"tiny={root / 'tiny.txt'}",
                              "--background", str(root / "bg.txt"),
                              "--output_dir", str(out)])
            s = pd.read_csv(out / "enrich_summary.tsv", sep="\t")
            de_row = s[s["Set"] == "DE"].iloc[0]
            self.assertEqual(int(de_row["N_Genes"]), 12)       # 首基因不得被吞
            self.assertEqual(int(de_row["N_Sig_Padj05"]), 2)   # GO + KEGG_KO
            t = pd.read_csv(out / "enrich_DE.tsv", sep="\t")
            self.assertTrue(bool(t[t["Term"] == "GO:0001111"]["Significant"].iloc[0]))
            self.assertFalse(bool(t[t["Term"] == "GO:0002222"]["Significant"].iloc[0]))
            # 统计空间回归 (2026-10-03 修复): N=可注释背景=60, n=可注释查询=10
            from scipy.stats import hypergeom
            self.assertAlmostEqual(
                float(t[t["Term"] == "GO:0001111"]["Pvalue"].iloc[0]),
                float(hypergeom.sf(9 - 1, 60, 30, 10)), places=10)
            self.assertTrue((out / "enrich_tiny.tsv").is_file())  # 跳过也落文件

    def test_generic_map_and_first_gene_preserved(self):
        import enrich_analysis as ea
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            t1, t2 = self._synthetic(root)
            gm = ["gene\tterm\tnamespace\tname"] + \
                 [f"{g}\tMAPMAN_1\tmapman\tphotosynthesis" for g in sorted(t1)] + \
                 [f"{g}\tMAPMAN_2\tmapman\tstress" for g in sorted(t2)]
            (root / "generic.tsv").write_text("\n".join(gm) + "\n",
                                              encoding="utf-8")
            out = root / "out"
            ea.main_testable(["--annotation", str(root / "generic.tsv"),
                              "--gene-set", f"DE={root / 'DE.txt'}",
                              "--background", str(root / "bg.txt"),
                              "--output_dir", str(out)])
            t = pd.read_csv(out / "enrich_DE.tsv", sep="\t")
            row = t[t["Term"] == "MAPMAN_1"].iloc[0]
            self.assertTrue(bool(row["Significant"]))
            self.assertEqual(row["Name"], "photosynthesis")

    def test_read_gene_list_formats(self):
        import enrich_analysis as ea
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "plain.txt").write_text("gA\ngB\ngA\n", encoding="utf-8")
            self.assertEqual(ea.read_gene_list(root / "plain.txt"),
                             ["gA", "gB"])                     # 去重保序, 无表头
            (root / "tab.tsv").write_text("gene\tx\ngA\t1\ngB\t2\n",
                                          encoding="utf-8")
            self.assertEqual(ea.read_gene_list(root / "tab.tsv"), ["gA", "gB"])

    def test_collect_enrich_sets(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            de = root / "06_DE"
            co = root / "05_Coexp"
            de.mkdir(parents=True)
            (co / "modules").mkdir(parents=True)
            (de / "deseq2_degs.tsv").write_text(
                "gene\tlog2FoldChange\tpadj\ndA\t2\t0.01\ndB\t1.5\t0.02\n"
                "dC\t-1.2\t0.03\n", encoding="utf-8")
            (co / "diff_coexp_gene_aggregate.tsv").write_text(
                "Gene\tn_sig_pairs\nhA\t5\nhB\t3\n", encoding="utf-8")
            (co / "modules" / "module_table_positive.tsv").write_text(
                "Gene\tModule\tDegree\nm1\t0\t2\nm2\t0\t1\nm3\t1\t1\n",
                encoding="utf-8")
            edir = root / "07_Enrich"
            found = hexp._collect_enrich_sets(
                {"de": de, "coexp": co}, edir, hub_n=50)
            self.assertEqual(set(found), {"DE", "DE_up", "DE_down", "dz_hub",
                                          "module_positive_0",
                                          "module_positive_1"})
            self.assertEqual((edir / "sets" / "DE.txt").read_text(
                encoding="utf-8").split(), ["dA", "dB", "dC"])
            self.assertEqual((edir / "sets" / "DE_up.txt").read_text(
                encoding="utf-8").split(), ["dA", "dB"])
            self.assertEqual((edir / "sets" / "DE_down.txt").read_text(
                encoding="utf-8").split(), ["dC"])
            self.assertEqual((edir / "sets" / "dz_hub.txt").read_text(
                encoding="utf-8").split(), ["hA", "hB"])
            self.assertEqual((edir / "sets" / "module_positive_1.txt")
                             .read_text(encoding="utf-8").split(), ["m3"])


class TestRankTransformMp(unittest.TestCase):
    """rank_transform_mp: 输出与 vendored 串行 Rank_transform 逐字节一致。"""

    def test_byte_identical_with_serial(self):
        import os
        import random
        import hashlib
        import shutil
        import subprocess
        rng = random.Random(7)
        n = 40
        genes = [f"g{i:03d}" for i in range(n)]
        shuffled = genes[:]
        rng.shuffle(shuffled)
        with tempfile.TemporaryDirectory() as td:
            root0 = Path(td) / "work0"
            raw = root0 / "RAW_GCN" / "PCC" / "k2"
            raw.mkdir(parents=True)
            for src in shuffled:
                lines = [f"# {src}"]
                for tgt in shuffled:
                    if tgt == src:
                        lines.append(f"{tgt}\t{1.0:.5f}\t1")
                        continue
                    if rng.random() < 0.05:          # 单向缺失 → NaN
                        continue
                    lines.append(f"{tgt}\t{rng.uniform(-1, 1):.8f}"
                                 f"\t{rng.randint(2, n)}")
                (raw / src).write_text("\n".join(lines) + "\n",
                                       encoding="utf-8")
            outs = []
            for tag in ("serial", "mp"):
                root = Path(td) / f"work_{tag}"
                shutil.copytree(root0, root)
                if tag == "serial":
                    env = dict(os.environ,
                               PYTHONPATH=str(PIPELINE_DIR / "teagcn"))
                    r = subprocess.run(
                        [sys.executable, "main/Rank_transform.py",
                         "-o", str(root)],
                        cwd=str(PIPELINE_DIR / "teagcn"), env=env,
                        capture_output=True, text=True, timeout=300)
                else:
                    r = subprocess.run(
                        [sys.executable,
                         str(PIPELINE_DIR / "rank_transform_mp.py"),
                         "-o", str(root), "-w", "3"],
                        capture_output=True, text=True, timeout=300)
                self.assertEqual(r.returncode, 0, r.stderr[-300:])
                out = root / "Completed_GCN" / "k2"
                outs.append({p.name: hashlib.md5(p.read_bytes()).hexdigest()
                             for p in out.iterdir() if p.is_file()})
            self.assertEqual(len(outs[0]), n)
            self.assertEqual(outs[0], outs[1])           # 逐文件 md5 一致


class TestExpandStages(unittest.TestCase):
    """--stage 展开: 'all' 默认含 loadscan, --no-with-loadscan 才排除。"""

    def test_all_includes_loadscan_by_default(self):
        st = hexp._expand_stages(["all"], with_loadscan=True)
        self.assertIn("loadscan", st)
        self.assertEqual(st[0], "hostref")               # 保序

    def test_loadscan_disabled_explicitly(self):
        self.assertNotIn("loadscan", hexp._expand_stages(["all"], False))
        # 显式点名永远生效
        self.assertIn("loadscan", hexp._expand_stages(["all", "loadscan"], False))
        self.assertEqual(hexp._expand_stages(["loadscan"], False), ["loadscan"])

    def test_dedupe_preserves_order(self):
        # 去重保首现位置: 先点名的 coexp 保持最前, all 只补未出现的阶段
        st = hexp._expand_stages(["coexp", "all", "de"], with_loadscan=False)
        self.assertEqual(st, ["coexp", "hostref", "quant", "qc", "labels",
                              "de", "enrich", "report"])


class TestHostExpressionPlots(unittest.TestCase):
    """host_expression_plots: 合成产物出图 (Agg 后端, 只验文件产出)。"""

    @staticmethod
    def _synthetic(root: Path):
        cdir = root / "05_Coexp"
        (cdir / "modules").mkdir(parents=True)
        (root / "07_Enrich").mkdir(parents=True)
        pd.DataFrame({"Gene": [f"g{i}" for i in range(15)],
                      "n_sig_pairs": range(15, 0, -1),
                      "mean_delta_z": [0.5] * 15,
                      "max_abs_delta_z": [2.0] * 15}).to_csv(
            cdir / "diff_coexp_gene_aggregate.tsv", sep="	", index=False)
        rows = []
        for i, (g, m) in enumerate([("gA", 0), ("gB", 0), ("gC", 1),
                                    ("gD", 1), ("gE", 2)]):
            rows.append({"Group": "positive", "Module": m, "N_Genes": 2,
                         "Virus": "GCVA", "N_Samples": 20,
                         "Rho": -0.5 + 0.2 * i,
                         "Pvalue": 1e-6, "Padj": 10.0 ** (-4 - i)})
        pd.DataFrame(rows).to_csv(
            cdir / "modules" / "module_load_assoc.tsv", sep="	", index=False)
        terms = [("GO:0001", "resp A"), ("GO:0002", "resp B"),
                 ("GO:0003", "resp C"), ("GO:0004", "resp D")]
        for set_name in ("DE_up", "DE_down"):
            pd.DataFrame([{"Set": set_name, "Term": t, "Namespace": "GO",
                           "Name": n, "N_Set": 10, "N_Set_in_Term": 5,
                           "N_Term_BG": 20, "N_BG": 100, "Fold": 2.5 - 0.2 * i,
                           "Pvalue": 1e-5, "Padj": 0.005 - 1e-4 * i,
                           "Significant": True}
                          for i, (t, n) in enumerate(terms)]).to_csv(
                root / "07_Enrich" / f"enrich_{set_name}.tsv", sep="	",
                index=False)
        # 模块网络: 2 模块, 环形边
        genes = [f"m{i}g{j}" for i in range(2) for j in range(10)]
        pd.DataFrame({"Gene": genes,
                      "Module": [i // 10 for i in range(20)],
                      "Degree": [5] * 20}).to_csv(
            cdir / "modules" / "module_table_positive.tsv", sep="	",
            index=False)
        e = []
        for i in range(2):
            mem = [f"m{i}g{j}" for j in range(10)]
            for a in range(10):
                e.append((mem[a], mem[(a + 1) % 10], 2.5, i, i))
        pd.DataFrame(e, columns=["Gene1", "Gene2", "Z",
                                 "Module1", "Module2"]).to_csv(
            cdir / "modules" / "module_edges_positive.tsv", sep="	",
            index=False)
        # fig01 (Δz 分布) 数据源
        rng = np.random.default_rng(3)
        dz = np.concatenate([rng.normal(-3.5, 0.3, 200),
                             rng.normal(3.5, 0.3, 200)])
        pd.DataFrame({"Delta_z": dz}).to_csv(
            cdir / "diff_coexp_candidates_sampled.tsv", sep="	", index=False)
        pd.DataFrame([{"n_genes": 100, "top_n": 20, "q_low": 0.5,
                       "q_high": 99.5, "threshold_low": -2.9,
                       "threshold_high": 2.9, "n_significant_rows": 400,
                       "n_unique_pairs": 380, "n_enhanced": 200,
                       "n_weakened": 200, "n_genes_involved": 90}]).to_csv(
            cdir / "diff_coexp_summary.tsv", sep="	", index=False)

    def test_all_figs_render(self):
        import host_expression_plots as hep
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._synthetic(root)
            out = root / "figs"
            hep.main_testable(["--root", str(root), "--out-dir", str(out)])
            made = {f.name for f in out.iterdir() if f.suffix == ".png"}
            self.assertIn("fig01_delta_z_distribution.png", made)
            self.assertIn("fig02_module_load.png", made)
            self.assertIn("fig03_enrichment.png", made)
            self.assertIn("fig04_module_network_positive.png", made)


class TestReportSmoke(unittest.TestCase):
    def test_build_report_with_and_without_artifacts(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            html = her.build_report(root)          # 全缺 → 降级
            self.assertIn("MMPV 宿主表达报告", html)
            self.assertIn("（无此段产物）", html)

            (root / "03_QC").mkdir()
            (root / "03_QC/qc_summary.tsv").write_text(
                "sample\tnum_processed\tnum_mapped\tpercent_mapped\tpass\tnote\n"
                "S1\t2000000\t1400000\t70.00\t1\t\n"
                "S2\t500\t0\t0.00\t0\tlow_mapped\n", encoding="utf-8")
            (root / "04_Labels").mkdir()
            (root / "04_Labels/sample_labels.tsv").write_text(
                "Sample\tVirus_any\tN_virus\tGCVA_positive\n"
                "S1\t+\t1\t1\nS2\t-\t0\t0\n", encoding="utf-8")
            (root / "06_DE").mkdir()
            (root / "06_DE/deseq2_summary.tsv").write_text(
                "contrast\tn_deg\tcontrast\tn_up\tn_down\n"
                "Virus_any: + vs -\t3\t0\t3\n", encoding="utf-8")
            html = her.build_report(root)
            self.assertIn("class='fail'", html)     # QC 不通过行高亮
            self.assertIn("class='card'", html)
            self.assertIn("1 / 1", html)            # 阳性/阴性卡片


# ═══════════ KingdomStress 增量 (P1/P2/P3): 网络 z-QC / 模块 / 保守性 ═══════════

def _write_gcn_dir(root: Path, name: str, genes: list, zmap):
    """zmap: {gene: {target: z}}; 未给的填 0.5; 自边填 9.0 (模拟 TEA-GCN 逐基因文件)。"""
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    for g in genes:
        zvals = [9.0 if t == g else zmap.get(g, {}).get(t, 0.5) for t in genes]
        pd.DataFrame({"Target": genes, "Co-exp_Str": [0.8] * len(genes),
                      "Co-exp_Str_MR": [1.0] * len(genes),
                      "zScore(Co-exp_Str_MR)": zvals}).to_csv(
            d / g, sep="\t", index=False)
    return d


def _planted_networks(root: Path, n_bg: int = 40):
    """40 背景基因 + 两个 DEG 三角团 {g1,g2,g3}/{g4,g5,g6}; 返回 (dpos, dneg)。"""
    rng = np.random.default_rng(7)
    genes = [f"b{i:02d}" for i in range(n_bg)] + \
            ["g1", "g2", "g3", "g4", "g5", "g6"]
    zp = {g: {t: round(float(rng.normal(0, 0.2)), 3)
              for t in genes if t != g} for g in genes}
    zn = {g: {t: round(float(rng.normal(0, 0.2)), 3)
              for t in genes if t != g} for g in genes}
    for a, b in (("g1", "g2"), ("g2", "g3"), ("g1", "g3"),
                 ("g4", "g5"), ("g5", "g6"), ("g4", "g6")):
        zp[a][b] = 3.0
        zp[b][a] = 3.0
    dpos = _write_gcn_dir(root, "pos", genes, zp)
    dneg = _write_gcn_dir(root, "neg", genes, zn)
    return dpos, dneg


class TestNetworkQC(unittest.TestCase):
    """P1: 网络 z 分布 QC (KingdomStress Fig S8 物种剔除规则的工程化)。"""

    def test_status_pass_fail_and_columns(self):
        import coexp_network_qc as cnq
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            dpos, dneg = _planted_networks(root)
            out = root / "coexp"
            cnq.main_testable(["--gcn-dir", f"positive={dpos}",
                               "--gcn-dir", f"negative={dneg}",
                               "--output_dir", str(out)])      # warn 模式 exit 0
            q = pd.read_csv(out / "network_qc.tsv", sep="\t")
            self.assertEqual(list(q["Group"]), ["positive", "negative"])
            pos = q[q["Group"] == "positive"].iloc[0]
            neg = q[q["Group"] == "negative"].iloc[0]
            self.assertEqual(pos["Status"], "PASS")            # 植入 z=3.0
            self.assertEqual(float(pos["Max_z"]), 3.0)
            self.assertEqual(int(pos["N_z_ge_1.7"]), 12)       # 6 边 × 双向
            self.assertEqual(neg["Status"], "FAIL")            # 纯噪声 max<1.5
            self.assertLess(float(neg["Max_z"]), 1.5)

    def test_on_fail_error_exit(self):
        import coexp_network_qc as cnq
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            _, dneg = _planted_networks(root)
            with self.assertRaises(SystemExit) as ctx:
                cnq.main_testable(["--gcn-dir", f"negative={dneg}",
                                   "--output_dir", str(root / "o"),
                                   "--on-fail", "error"])
            self.assertEqual(ctx.exception.code, 3)


class TestCoexpModules(unittest.TestCase):
    """P2: DEG 限定模块检测 (z>=1.7 → Louvain) + 模块×载量关联。"""

    @staticmethod
    def _write_degs(root: Path):
        degs = ["g1", "g2", "g3", "g4", "g5", "g6"]
        pd.DataFrame({"log2FoldChange": [2.0] * 6, "padj": [0.01] * 6},
                     index=degs).to_csv(root / "degs.tsv", sep="\t",
                                        index_label="gene")
        return degs

    def test_module_recovery_and_deg_restriction(self):
        import coexp_modules as cm
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            dpos, dneg = _planted_networks(root)
            degs = self._write_degs(root)
            out = root / "modules"
            cm.main_testable(["--gcn-dir", f"positive={dpos}",
                              "--gcn-dir", f"negative={dneg}",
                              "--deg", str(root / "degs.tsv"),
                              "--output_dir", str(out)])
            s = pd.read_csv(out / "module_summary.tsv", sep="\t")
            self.assertEqual(len(s), 2)                        # 两个植入团
            self.assertTrue((s[s["Group"] == "positive"]["N_Genes"] == 3).all())
            t = pd.read_csv(out / "module_table_positive.tsv", sep="\t")
            mod = dict(zip(t["Gene"], t["Module"]))
            self.assertEqual(set(mod), set(degs))              # 非 DEG 不得入表
            self.assertEqual(mod["g1"], mod["g2"])
            self.assertEqual(mod["g2"], mod["g3"])
            self.assertEqual(mod["g4"], mod["g5"])
            self.assertNotEqual(mod["g1"], mod["g4"])          # 两团分开
            # negative 全噪声 → 无 z>=1.7 边 → 无模块行
            self.assertEqual(len(s[s["Group"] == "negative"]), 0)
            # 边表导出: 植入对在列, 同模块, 权重=z
            e = pd.read_csv(out / "module_edges_positive.tsv", sep="	")
            planted_e = e[((e["Gene1"] == "g1") & (e["Gene2"] == "g2")) |
                          ((e["Gene1"] == "g2") & (e["Gene1"].eq("g1")))]
            self.assertGreaterEqual(len(planted_e), 1)
            self.assertEqual(float(planted_e.iloc[0]["Z"]), 3.0)
            self.assertEqual(planted_e.iloc[0]["Module1"],
                             planted_e.iloc[0]["Module2"])

    def test_qc_gate_skips_low_z_group(self):
        import coexp_modules as cm
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            dpos, dneg = _planted_networks(root)
            self._write_degs(root)
            qc = root / "coexp"
            qc.mkdir()
            pd.DataFrame({"Group": ["positive", "negative"],
                          "Max_z": [3.0, 0.6],
                          "Status": ["PASS", "FAIL"]}).to_csv(
                qc / "network_qc.tsv", sep="\t", index=False)
            out = root / "modules"
            cm.main_testable(["--gcn-dir", f"positive={dpos}",
                              "--gcn-dir", f"negative={dneg}",
                              "--deg", str(root / "degs.tsv"),
                              "--network-qc", str(qc / "network_qc.tsv"),
                              "--output_dir", str(out)])
            s = pd.read_csv(out / "module_summary.tsv", sep="\t")
            self.assertEqual(set(s["Group"]), {"positive"})    # FAIL 组被跳过

    def test_load_association(self):
        import coexp_modules as cm
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            dpos, _ = _planted_networks(root)
            self._write_degs(root)
            n_bg, n = 40, 12
            load = np.arange(1, n + 1, dtype=float)
            samples = [f"S{i:02d}" for i in range(1, n + 1)]
            tpm = pd.DataFrame({s: [10.0] * (n_bg + 6) for s in samples},
                               index=[f"b{i:02d}" for i in range(n_bg)] +
                                     ["g1", "g2", "g3", "g4", "g5", "g6"])
            for k, g in enumerate(["g1", "g2", "g3"]):         # 植入: 随载量上升
                rng = np.random.default_rng(11 + k)
                tpm.loc[g] = load * 10 + rng.normal(0, 0.5, n)
            for g in ["g4", "g5", "g6"]:                       # 对照: 与载量无关
                rng = np.random.default_rng(21 + ord(g[-1]))
                tpm.loc[g] = 10.0 + rng.normal(0, 1.0, n)
            tpm.to_csv(root / "TPM.tsv", sep="\t")
            pd.DataFrame({"Sample": samples, "Virus": "GCVA",
                          "Positive": 1, "Load": load}).to_csv(
                root / "virus_labels.tsv", sep="\t", index=False)
            out = root / "modules"
            cm.main_testable(["--gcn-dir", f"positive={dpos}",
                              "--deg", str(root / "degs.tsv"),
                              "--tpm", str(root / "TPM.tsv"),
                              "--virus-labels", str(root / "virus_labels.tsv"),
                              "--output_dir", str(out)])
            assoc = pd.read_csv(out / "module_load_assoc.tsv", sep="\t")
            self.assertEqual(len(assoc), 2)                    # 2 模块 × 1 病毒
            self.assertIn("Padj", assoc.columns)
            self.assertGreater(float(assoc["Rho"].max()), 0.9)  # 植入信号


class TestModuleConservation(unittest.TestCase):
    """P3: 跨宿主模块保守性 (OC + 物种内置换零模型)。"""

    @staticmethod
    def _write_modules(path: Path, rows):
        pd.DataFrame(rows, columns=["Gene", "Module", "Degree"]).to_csv(
            path, sep="\t", index=False)

    def test_conserved_pair_and_null(self):
        import coexp_module_conservation as cmc
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._write_modules(root / "ma.tsv",
                                [("g1", 0, 2), ("g2", 0, 2), ("g3", 0, 2),
                                 ("g4", 1, 1), ("g5", 1, 1), ("g6", 1, 1)])
            self._write_modules(root / "mb.tsv",
                                [("h1", 0, 2), ("h2", 0, 2), ("h3", 0, 2),
                                 ("h4", 1, 1), ("h5", 1, 1), ("h6", 1, 1)])
            og = [(f"OG{i:05d}", f"g{i}, h{i}") for i in range(1, 4)]  # 团0 完全重叠
            og += [(f"OGA{i:05d}", f"g{i}") for i in range(4, 7)]      # 团1 各自独立
            og += [(f"OGB{i:05d}", f"h{i}") for i in range(4, 7)]
            pd.DataFrame(og, columns=["Orthogroup", "hostA_hostB"]).to_csv(
                root / "Orthogroups.tsv", sep="\t", index=False)
            out = root / "cons"
            cmc.main_testable(["--modules-a", str(root / "ma.tsv"),
                               "--modules-b", str(root / "mb.tsv"),
                               "--name-a", "hostA", "--name-b", "hostB",
                               "--orthogroups", str(root / "Orthogroups.tsv"),
                               "--species-a", "hostA_hostB",
                               "--species-b", "hostA_hostB",
                               "--min-shared-og", "3", "--permutations", "50",
                               "--output_dir", str(out)])
            pairs = pd.read_csv(out / "conservation_pairs.tsv", sep="\t")
            self.assertEqual(len(pairs), 1)                    # 弱对被门槛滤掉
            row = pairs.iloc[0]
            self.assertEqual(int(row["Shared_OG"]), 3)
            self.assertEqual(float(row["OC"]), 1.0)
            self.assertTrue(bool(row["Conserved"]))
            self.assertLess(float(row["P_Emp"]), 0.1)          # 完美重叠 vs 置换
            summ = pd.read_csv(out / "conservation_summary.tsv", sep="\t")
            self.assertEqual(int(summ.iloc[0]["n_conserved"]), 1)

    def test_unmapped_genes_reported(self):
        import coexp_module_conservation as cmc
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._write_modules(root / "ma.tsv", [("g1", 0, 1), ("gX", 0, 1)])
            self._write_modules(root / "mb.tsv", [("h1", 0, 1)])
            pd.DataFrame([("OG00001", "g1, h1")],
                         columns=["Orthogroup", "sp"]).to_csv(
                root / "Orthogroups.tsv", sep="\t", index=False)
            out = root / "cons"
            cmc.main_testable(["--modules-a", str(root / "ma.tsv"),
                               "--modules-b", str(root / "mb.tsv"),
                               "--name-a", "A", "--name-b", "B",
                               "--orthogroups", str(root / "Orthogroups.tsv"),
                               "--species-a", "sp", "--species-b", "sp",
                               "--min-shared-og", "1", "--permutations", "10",
                               "--output_dir", str(out)])
            summ = pd.read_csv(out / "conservation_summary.tsv", sep="\t")
            self.assertEqual(int(summ.iloc[0]["unmapped_genes_a"]), 1)
            pairs = pd.read_csv(out / "conservation_pairs.tsv", sep="\t")
            self.assertEqual(float(pairs.iloc[0]["OC"]), 1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
