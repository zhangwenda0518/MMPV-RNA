#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_graph_propagation.py — graph_propagation.py 纯逻辑单元测试
====================================================================
直接运行:  python virome_discovery_pipeline/tests/test_graph_propagation.py
或 pytest: python -m pytest virome_discovery_pipeline/tests/test_graph_propagation.py -q

只测纯函数与 tmp 目录 I/O, 不依赖 bwa-mem2/samtools/coverm 等外部工具。
"""

import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "virome_discovery_pipeline"))

import graph_propagation as gp  # noqa: E402


def _mk_fasta(path, records):
    """records: [(id, seq)] → 写 FASTA, 返回 {id: (len, gc)} 期望值。"""
    with open(path, 'w') as fh:
        for cid, seq in records:
            fh.write(f">{cid} desc\n{seq}\n")
    gc = {cid: (len(s), sum(1 for ch in s.upper() if ch in 'GC') / len(s))
          for cid, s in records}
    return gc


class TestFastaRead(unittest.TestCase):
    def test_lengths_gc(self):
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "a.fasta"
            _mk_fasta(fa, [("c1", "ACGTACGTAC"), ("c2", "GGCC")])
            got = gp.read_fasta_lengths_gc(fa)
            self.assertEqual(got["c1"][0], 10)
            self.assertAlmostEqual(got["c1"][1], 0.5)
            self.assertEqual(got["c2"][0], 4)
            self.assertAlmostEqual(got["c2"][1], 1.0)

    def test_gz_input(self):
        import gzip
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "a.fasta.gz"
            with gzip.open(fa, 'wt') as fh:
                fh.write(">x\nAC\n")
            self.assertEqual(gp.read_fasta_lengths_gc(fa)["x"][0], 2)


class TestSeedVotes(unittest.TestCase):
    def test_votes_and_clean_variant(self):
        with tempfile.TemporaryDirectory() as td:
            ident = Path(td)
            sdir = ident / "ERR1"
            sdir.mkdir()
            (sdir / "ERR1_virus.genomad.result.id").write_text("A\nB\n")
            (sdir / "ERR1_virus.virbot.result.id").write_text("A\n")
            (sdir / "ERR1_viroids.result.id").write_text("D\n")
            votes, viroid, n, hit = gp.load_seed_votes(
                ident, ["ERR1_clean", "ERR1"])
            self.assertEqual(n, "ERR1")
            self.assertEqual(votes["A"], 2)
            self.assertEqual(votes["B"], 1)
            self.assertNotIn("C", votes)
            self.assertEqual(viroid, {"D"})

    def test_missing_returns_empty(self):
        with tempfile.TemporaryDirectory() as td:
            votes, viroid, n, hit = gp.load_seed_votes(Path(td), ["NOPE"])
            self.assertIsNone(n)
            self.assertEqual(votes, {})
            self.assertEqual(viroid, set())

    def test_all_union_file_not_counted_as_tool(self):
        """all.result.id 是 ensemble 并集派生文件, 不得按独立工具计票 —
        计入会虚增 1 票, 使 1 票弱阳性被错误保护成 viral_seed (246 mock 实测)。"""
        with tempfile.TemporaryDirectory() as td:
            sdir = Path(td) / "S1"
            sdir.mkdir()
            (sdir / "S1_virus.genomad.result.id").write_text("A\n")
            (sdir / "S1_virus.virbot.result.id").write_text("A\nB\n")
            (sdir / "S1_virus.all.result.id").write_text("A\nB\nC\n")
            votes, _, n, _ = gp.load_seed_votes(Path(td), ["S1"])
            self.assertEqual(votes["A"], 2)   # 两工具 + 并集不计
            self.assertEqual(votes["B"], 1)   # 1 票弱阳性不被并集虚增
            self.assertNotIn("C", votes)

    def test_uniprot_filtered_ids_not_counted(self):
        """{N}_virus.{tool}.uniprot_filtered.id 不匹配 *.result.id, 不应计票。"""
        with tempfile.TemporaryDirectory() as td:
            sdir = Path(td) / "S1"
            sdir.mkdir()
            (sdir / "S1_virus.genomad.result.id").write_text("A\n")
            (sdir / "S1_virus.genomad.uniprot_filtered.id").write_text("B\n")
            votes, _, n, _ = gp.load_seed_votes(Path(td), ["S1"])
            self.assertIn("A", votes)
            self.assertNotIn("B", votes)


class TestAssignSeeds(unittest.TestCase):
    def test_three_classes(self):
        lengths = {"A": 1500, "B": 1500, "C": 1500, "D": 800, "E": 500}
        votes = {"A": 2, "B": 1}
        seeds = gp.assign_seeds(lengths, lengths, votes, {"D"},
                                seed_min_votes=1, nonviral_min_len=1000)
        self.assertEqual(seeds["A"], gp.V_VIRAL_SEED)
        self.assertEqual(seeds["B"], gp.V_VIRAL_SEED)
        self.assertEqual(seeds["C"], gp.V_NONVIRAL_SEED)
        self.assertEqual(seeds["D"], gp.V_UNKNOWN)   # viroid 阳性: 两边都不站
        self.assertEqual(seeds["E"], gp.V_UNKNOWN)   # 短 contig 未被工具检查

    def test_seed_min_votes_two_weak_positive_unknown(self):
        """1 票弱阳性必须归 unknown (交图裁决), 不得钉成 nonviral_seed —
        否则 filter 降级与 COBRA 拦截两个消费点永远无法触发。"""
        lengths = {"A": 1500, "B": 1500, "C": 1500}
        votes = {"A": 2, "B": 1}
        seeds = gp.assign_seeds(lengths, lengths, votes, set(),
                                seed_min_votes=2, nonviral_min_len=1000)
        self.assertEqual(seeds["A"], gp.V_VIRAL_SEED)
        self.assertEqual(seeds["B"], gp.V_UNKNOWN)        # 弱阳性 → 图裁决
        self.assertEqual(seeds["C"], gp.V_NONVIRAL_SEED)  # 0 票且被检查过 → 非病毒种子


class TestCountCrossLinks(unittest.TestCase):
    @staticmethod
    def _sam(flag, rname, rnext):
        return f"read1\t{flag}\t{rname}\t100\t60\t50M\t{rnext}\t200\t0\tACGT\t*\n"

    def test_cross_contig_proper(self):
        lines = [self._sam(99, "A", "B"), self._sam(99, "B", "A"),
                 self._sam(99, "A", "=")]      # 同 contig ('='), 不建边
        links = gp.count_cross_links(lines, require_proper=True)
        self.assertEqual(links[("A", "B")], 2)
        self.assertEqual(len(links), 1)

    def test_filters(self):
        lines = [
            self._sam(73, "A", "*"),            # mate unmapped (0x8)
            self._sam(99 | 0x100, "A", "B"),    # secondary
            self._sam(99 | 0x800, "A", "B"),    # supplementary
            self._sam(65, "A", "B"),            # 非 proper
            self._sam(99, "A", "B"),            # 合格
        ]
        self.assertEqual(gp.count_cross_links(lines, require_proper=True)
                         .get(("A", "B")), 1)
        got = gp.count_cross_links(lines, require_proper=False)
        self.assertEqual(got.get(("A", "B")), 2)  # 放宽后 65 与 99 都收

    def test_non_proper_cross_contig_is_the_real_default(self):
        """bwa 的 0x2 仅同 contig 置位 — 跨 contig pair 永远 non-proper (flag 65),
        require_proper=False (生产默认) 必须能收进来 (246 实测: proper 跨 contig=0)。"""
        lines = [self._sam(65, "A", "B"), self._sam(129, "B", "A")]  # both-mapped 非proper
        self.assertEqual(gp.count_cross_links(lines, require_proper=False)
                         .get(("A", "B")), 1)
        self.assertEqual(gp.count_cross_links(lines, require_proper=True), {})

    def test_header_and_malformed_ignored(self):
        links = gp.count_cross_links(["@SQ\tSN:A\n", "", "short\tline\n"])
        self.assertEqual(links, {})


class TestPruneEdges(unittest.TestCase):
    def setUp(self):
        self.gc = {"v": 0.5, "host": 0.2, "mid": 0.4}
        self.depth = {"v": 10.0, "host": 1.0, "mid": 1.0, "zero": 0.0}

    def test_weight_threshold(self):
        links = {("v", "mid"): 2, ("v", "host"): 3}
        edges = gp.prune_edges(links, self.gc, self.depth, min_link_pairs=3,
                               gc_diff=1.0, depth_ratio=1e9, has_depth=False)
        self.assertEqual(edges, {("v", "host"): 3})

    def test_gc_diff(self):
        links = {("v", "host"): 5, ("v", "mid"): 5}
        edges = gp.prune_edges(links, self.gc, self.depth, min_link_pairs=1,
                               gc_diff=0.20, depth_ratio=1e9, has_depth=False)
        self.assertEqual(set(edges), {("v", "mid")})   # |0.5-0.2|=0.3>0.2 断

    def test_depth_ratio_and_zero(self):
        links = {("v", "mid"): 5, ("v", "host"): 5, ("v", "zero"): 5}
        # v-mid 10:1 (≤20 保留); v-host 100:1 (>20 断); v-zero 一端 0 深度 → 断
        depth = {"v": 10.0, "mid": 1.0, "host": 0.1, "zero": 0.0}
        edges = gp.prune_edges(links, self.gc, depth, min_link_pairs=1,
                               gc_diff=1.0, depth_ratio=20.0, has_depth=True)
        self.assertEqual(set(edges), {("v", "mid")})

    def test_no_depth_table_skips_rule(self):
        links = {("v", "host"): 5}
        edges = gp.prune_edges(links, self.gc, {}, min_link_pairs=1,
                               gc_diff=1.0, depth_ratio=20.0, has_depth=False)
        self.assertEqual(set(edges), {("v", "host")})


class TestPropagate(unittest.TestCase):
    @staticmethod
    def _run(edges, contigs, seeds, **kw):
        return gp.propagate(edges, contigs, seeds, **kw)

    def test_viral_dominant(self):
        edges = {("v1", "u1"): 5, ("u1", "u2"): 5}
        contigs = ["v1", "u1", "u2", "iso"]
        seeds = {"v1": gp.V_VIRAL_SEED, "iso": gp.V_UNKNOWN}
        verdicts, comm, stats = self._run(edges, contigs, seeds)
        self.assertEqual(verdicts["v1"], gp.V_VIRAL_SEED)   # 种子不可变
        self.assertEqual(verdicts["u1"], gp.V_GRAPH_VIRAL)
        self.assertEqual(verdicts["u2"], gp.V_GRAPH_VIRAL)
        self.assertEqual(verdicts["iso"], gp.V_UNKNOWN)     # 孤立节点

    def test_nonviral_dominant_boundary_10x(self):
        # N=10, V=1: 10 >= 10*1 恰好达标 → nonviral-dominant; viral 种子不被改标
        contigs = ["n0"] + [f"n{i}" for i in range(10)] + ["v1"]
        edges = {("n0", f"n{i}"): 5 for i in range(10)}
        edges[("n0", "v1")] = 5
        seeds = {c: gp.V_NONVIRAL_SEED for c in contigs if c != "v1"}
        seeds["v1"] = gp.V_VIRAL_SEED
        verdicts, _, stats = self._run(edges, contigs, seeds)
        self.assertEqual(verdicts["v1"], gp.V_VIRAL_SEED)
        self.assertEqual(stats[0]['rule'], 'nonviral-dominant')

    def test_conflicting(self):
        edges = {("v1", "n1"): 5, ("n1", "u1"): 5}
        seeds = {"v1": gp.V_VIRAL_SEED, "n1": gp.V_NONVIRAL_SEED}
        verdicts, _, stats = self._run(edges, ["v1", "n1", "u1"], seeds)
        self.assertEqual(verdicts["u1"], gp.V_CONFLICTING)
        self.assertEqual(stats[0]['rule'], 'conflicting')

    def test_no_seeds_community(self):
        edges = {("a", "b"): 5}
        verdicts, _, stats = self._run(edges, ["a", "b"], {})
        self.assertEqual(verdicts["a"], gp.V_UNKNOWN)
        self.assertEqual(stats[0]['rule'], 'no-seeds')

    def test_all_unknown_without_seeds(self):
        verdicts, _, _ = self._run({}, ["x"], {})
        self.assertEqual(verdicts, {"x": gp.V_UNKNOWN})


class TestProcessToolGraphRows(unittest.TestCase):
    CFG = {'min_link_pairs': 1, 'gc_diff': 0.2, 'depth_ratio': 20.0,
           'require_proper': True, 'resolution': 1.0, 'seed_min_votes': 1,
           'nonviral_min_len': 1000, 'majority_ratio': 10.0, 'max_reads': 0,
           'keep_bam': False}

    def test_rows_contract_and_no_reads_degradation(self):
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "S1_megahit.contig.fasta"
            _mk_fasta(fa, [("S1_v1", "ACGT" * 300), ("S1_h1", "AC" * 600)])
            seeds = {"S1_v1": gp.V_VIRAL_SEED, "S1_h1": gp.V_NONVIRAL_SEED}
            rows, demoted, summary = gp.process_tool_graph(
                fa, "S1", "megahit", gp.read_fasta_lengths_gc(fa), seeds,
                {"S1_v1": 1}, Path(td) / "tmp", 1, None, None, self.CFG)
            self.assertEqual(set(rows[0]), set(gp.VERDICT_COLUMNS))   # 列契约
            self.assertEqual(demoted, set())                          # 无 reads 不降级
            by = {r['contig']: r for r in rows}
            self.assertEqual(by['S1_v1']['verdict'], gp.V_VIRAL_SEED)
            self.assertEqual(by['S1_h1']['verdict'], gp.V_NONVIRAL_SEED)
            self.assertEqual(len(summary), 2)                         # 两个孤立社区


class TestMapMinLen(unittest.TestCase):
    def test_write_filtered_reference(self):
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "a.fasta"
            _mk_fasta(fa, [("long1", "ACGT" * 400), ("short1", "ACGT" * 100),
                           ("long2", "ACGT" * 500)])
            lengths = gp.read_fasta_lengths_gc(fa)
            dst, kept = gp.write_filtered_reference(fa, lengths, 1000, Path(td))
            text = Path(dst).read_text()
            self.assertEqual(kept, 2)
            self.assertIn(">long1", text)
            self.assertIn(">long2", text)
            self.assertNotIn("short1", text)

    def test_cfg_carries_map_min_len(self):
        cfg = {'map_min_len': 1000}
        self.assertEqual(cfg.get('map_min_len', 0), 1000)


class TestDerivedCandidates(unittest.TestCase):
    def test_subtract_and_copy(self):
        with tempfile.TemporaryDirectory() as td:
            ident = Path(td) / "ident"
            sdir = ident / "S1"
            sdir.mkdir(parents=True)
            src = sdir / "S1_virus.all.candidate.fasta"
            with open(src, 'w') as fh:
                fh.write(">keep1 d\nACGT\n>drop d\nTTTT\n>keep2 d\nGGGG\n")
            root = Path(td) / "cand"
            _, dst, removed = gp.write_derived_candidates(
                ident, "S1", {"drop"}, root)
            self.assertEqual(removed, 1)
            text = Path(dst).read_text()
            self.assertIn(">keep1", text)
            self.assertNotIn(">drop", text)
            self.assertIn(">keep2", text)
            # 无降级 → 逐字节一致
            _, dst2, removed2 = gp.write_derived_candidates(ident, "S1", set(), root)
            self.assertEqual(removed2, 0)
            self.assertEqual(Path(dst2).read_bytes(), src.read_bytes())

    def test_fallback_source_name(self):
        with tempfile.TemporaryDirectory() as td:
            ident = Path(td)
            sdir = ident / "S1"
            sdir.mkdir()
            (sdir / "S1.virus.candidate.fasta").write_text(">a\nAC\n")
            _, dst, _ = gp.write_derived_candidates(ident, "S1", set(),
                                                    Path(td) / "cand")
            self.assertEqual(Path(dst).name, "S1_virus.all.candidate.fasta")

    def test_missing_returns_none(self):
        with tempfile.TemporaryDirectory() as td:
            src, dst, removed = gp.write_derived_candidates(
                Path(td), "NOPE", set(), Path(td) / "cand")
            self.assertIsNone(src)
            self.assertIsNone(dst)


class TestComputeAdjustments(unittest.TestCase):
    def test_demote_only_zero_evidence(self):
        kept = {"a", "b", "c"}
        cands = {"a", "b", "c", "d"}
        cdd_uncls = {"b", "c"}
        blastx = {"c"}                       # c 有 hit → 正面证据保护
        verdicts = {"b": gp.V_GRAPH_NONVIRAL, "c": gp.V_GRAPH_NONVIRAL,
                    "d": gp.V_GRAPH_VIRAL}
        final, demoted, rescued = gp.compute_adjustments(
            kept, cands, cdd_uncls, blastx, verdicts)
        self.assertEqual(demoted, {"b"})
        self.assertEqual(rescued, {"d"})
        self.assertEqual(final, {"a", "c", "d"})

    def test_rescue_only_dropped(self):
        final, demoted, rescued = gp.compute_adjustments(
            {"a"}, {"a"}, set(), set(), {"a": gp.V_GRAPH_VIRAL})
        self.assertEqual(rescued, set())
        self.assertEqual(final, {"a"})

    def test_switches(self):
        kept, cands = {"a", "b"}, {"a", "b", "c"}
        verdicts = {"b": gp.V_GRAPH_NONVIRAL, "c": gp.V_GRAPH_VIRAL}
        final, demoted, rescued = gp.compute_adjustments(
            kept, cands, {"b"}, set(), verdicts,
            do_demote=False, do_rescue=False)
        self.assertEqual((demoted, rescued, final), (set(), set(), kept))
        final2, demoted2, _ = gp.compute_adjustments(
            kept, cands, {"b"}, set(), verdicts, do_demote=False)
        self.assertEqual(demoted2, set())
        self.assertEqual(final2, {"a", "b", "c"})


class TestTruncateReads(unittest.TestCase):
    def test_fasta_prefix(self):
        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "r.fa"
            src.write_text(">a\nAC\n>b\nGT\n>c\nAA\n")
            dst = gp.truncate_reads(src, Path(td) / "t.fa", 2)
            text = dst.read_text()
            self.assertIn(">a", text)
            self.assertIn(">b", text)
            self.assertNotIn(">c", text)

    def test_fastq_four_line_groups(self):
        # 序列行以 @ 开头也不能破坏 4 行组逻辑
        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "r.fq"
            rec = "@r{i}\n{seq}\n+\n{qual}\n"
            src.write_text("".join(
                rec.format(i=i, seq="@" + "A" * 9, qual="I" * 10)
                for i in range(3)))
            dst = gp.truncate_reads(src, Path(td) / "t.fq", 1)
            self.assertEqual(dst.read_text().count("@r"), 1)

    def test_gz(self):
        import gzip
        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "r.fa.gz"
            with gzip.open(src, 'wt') as fh:
                fh.write(">a\nAC\n>b\nGT\n")
            dst = gp.truncate_reads(src, Path(td) / "t.fa", 1)
            self.assertIn(">a", dst.read_text())
            self.assertNotIn(">b", dst.read_text())


class TestReadPairMatching(unittest.TestCase):
    def test_pe_and_variants(self):
        with tempfile.TemporaryDirectory() as td:
            rd = Path(td)
            (rd / "S1_clean_1.fa.gz").write_text("")
            (rd / "S1_clean_2.fa.gz").write_text("")
            r1, r2 = gp.find_read_pair(rd, "S1_clean", ["S1_clean", "S1"])
            self.assertEqual(r1.name, "S1_clean_1.fa.gz")
            self.assertEqual(r2.name, "S1_clean_2.fa.gz")
            (rd / "S2_1.fa.gz").write_text("")
            r1, r2 = gp.find_read_pair(rd, "S2")
            self.assertIsNotNone(r1)
            self.assertIsNone(r2)
            self.assertEqual(gp.find_read_pair(rd, "S3"), (None, None))

    def test_contig_file_discovery(self):
        with tempfile.TemporaryDirectory() as td:
            sdir = Path(td) / "S1"
            sdir.mkdir()
            (sdir / "S1_megahit.contig.fasta").write_text(">a\nAC\n")
            (sdir / "S1_rnaviralspades.contig.fasta").write_text(">b\nAC\n")
            (sdir / "S1.all_tools.refineC.merged.fasta").write_text(">c\nAC\n")
            sub = sdir / "S1_megahit_refineC"
            sub.mkdir()
            (sub / "S1_megahit_refineC.split.fasta").write_text(">d\nAC\n")
            files = gp.find_contig_files(sdir)
            self.assertEqual([f.name for f in files],
                             ["S1_megahit.contig.fasta",
                              "S1_rnaviralspades.contig.fasta"])
            self.assertEqual(gp.tool_from_contig_file(files[0], "S1"), "megahit")


class TestFilterAdjudicationEnd2End(unittest.TestCase):
    """adjudicate_sample 三场景回归 (新主链 filter→graph, 合成布局)。"""

    def _mk_layout(self, td):
        out_root = td / "d_graph"
        sdir = td / "d_filter" / "ERR1"
        (sdir / "uniprot" / "diamond_uniprot").mkdir(parents=True, exist_ok=True)
        (sdir / "cdd").mkdir(exist_ok=True)
        (out_root / "ERR1").mkdir(parents=True, exist_ok=True)
        # graph verdicts
        (out_root / "ERR1" / "ERR1.graph_verdict.tsv").write_text(
            "contig\tsource_tool\tverdict\nb\tmegahit\tgraph_nonviral\n"
            "d\tmegahit\tgraph_viral\n")
        # ident 候选 (全量: a,b,c,d)
        ident = td / "ident"
        (ident / "ERR1").mkdir(parents=True, exist_ok=True)
        (ident / "ERR1" / "ERR1_virus.all.candidate.fasta").write_text(
            ">a\nAC\n>b\nAC\n>c\nAC\n>d\nAC\n")
        # filter 最终链 (uniprot 真实文件: a,b,c; b 为零证据, d 已被丢)
        chain = sdir / "ERR1.virus.candidate_filtered.fasta"
        chain.write_text(">a\nAC\n>b\nAC\n>c\nAC\n")
        (sdir / "cdd" / "cdd_calls.tsv").write_text(
            "contig\tcall\na\tviral\nb\tunclassified\nc\tunclassified\nd\tnon_viral_plant\n")
        (sdir / "uniprot" / "diamond_uniprot" / "blastx.tsv").write_text("a\tvirus\t100\n")
        return out_root, sdir, ident

    def test_demote_and_rescue(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            out_root, sdir, ident = self._mk_layout(Path(td))
            r = gp.adjudicate_sample("ERR1", str(out_root), str(sdir.parent),
                                     str(ident))
            self.assertEqual(r, (3, 1, 1))           # final=3, 降级 1 (b), 救回 1 (d)
            fin = out_root / "ERR1" / "ERR1.virus.candidate.final.fasta"
            self.assertEqual([l[1:] for l in fin.read_text().splitlines()
                              if l.startswith(">")], ["a", "c", "d"])
            mirror = sdir.parent / "candidates" / "ERR1" / \
                "ERR1_virus.all.candidate.fasta"
            self.assertEqual([l[1:] for l in mirror.read_text().splitlines()
                              if l.startswith(">")], ["a", "c", "d"])
            audit = out_root / "ERR1" / "ERR1.graph_adjust_calls.tsv"
            self.assertIn("b\tunclassified\tno\tgraph_nonviral\tdemoted",
                          audit.read_text())

    def test_switches_off_passthrough(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            out_root, sdir, ident = self._mk_layout(Path(td))
            r = gp.adjudicate_sample("ERR1", str(out_root), str(sdir.parent),
                                     str(ident), demote=False, rescue=False)
            self.assertEqual(r, (3, 0, 0))           # 无增删, 镜像=链
            mirror = sdir.parent / "candidates" / "ERR1" / \
                "ERR1_virus.all.candidate.fasta"
            self.assertEqual([l[1:] for l in mirror.read_text().splitlines()
                              if l.startswith(">")], ["a", "b", "c"])

    def test_rescue_needs_candidate_source(self):
        """ident 候选源缺失时救回放弃 (链中没有 d), 仅保留降级"""
        import tempfile
        import shutil
        with tempfile.TemporaryDirectory() as td:
            out_root, sdir, ident = self._mk_layout(Path(td))
            shutil.rmtree(ident)
            r = gp.adjudicate_sample("ERR1", str(out_root), str(sdir.parent),
                                     str(ident))
            self.assertEqual(r[2], 0)                # 救回放弃
            fin = out_root / "ERR1" / "ERR1.virus.candidate.final.fasta"
            self.assertEqual([l[1:] for l in fin.read_text().splitlines()
                              if l.startswith(">")], ["a", "c"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
