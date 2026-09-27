#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_eve_invariants.py — 判别器产物不变量
=========================================
这个文件把此前散落在仓库根目录 `eve_verify_*.py` 里的一次性核查固化成测试。原先是
"我跑过、数字在 README 里"的**叙述性证据**; 现在改成任何人都能重跑的**可复现证据**。

两层:
  1. 判定逻辑层 (总是跑): 用合成小表锁住各条不变量的**谓词**本身 —— 表长得对才谈得上
     产物对。离线、秒级。
  2. 真实产物层 (需要环境变量): 设 `EVE_VERIFY_DIR=<跑完的产物目录>` 时, 对真实表逐条
     检查同一批不变量。没设就 skip —— 单测不该依赖服务器上的数据。

用法:
  python -m pytest tests/test_eve_invariants.py -q
  EVE_VERIFY_DIR=/path/to/run python -m pytest tests/test_eve_invariants.py -q

为什么需要第二层: 第一层只能证明"如果表长成这样, 检查会通过"。真实数据里出现过
`locus_full` 挂在 `no_host_locus` 上(伪位点)、`best_locus_comps` 出现 `RT+A+P` 这种
按字符撕裂的垃圾值(334 条)、`locus_ncomp` 与 `locus_comps` 对不上 —— 这些只有拿真实
产物跑才抓得到。
"""
import csv
import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "eve_distinguish"))
import s2_domain_scan as s2        # noqa: E402
import s4_filter as s4             # noqa: E402

CANON = set(s2.CANON)
VERDICT_DIR = os.environ.get("EVE_VERIFY_DIR", "")


def load_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def comps_of(s):
    return set(c for c in (s or "").replace(";", "+").split("+") if c and c != "-")


# ─────────────────────── 1) 判定逻辑层（离线） ───────────────────────
class TestComponentTokens(unittest.TestCase):
    """组件列只允许出现 CANON 里的整词 —— 按字符撕裂会产出 'A'/'P'/'C' 这种垃圾."""

    def test_canon_tokens_only(self):
        self.assertEqual(comps_of("MP+CP+AP+RT+RH"), CANON)
        self.assertTrue(comps_of("MP+CP+AP+RT+RH") <= CANON)
        self.assertEqual(comps_of("-"), set())
        # 撕裂过的值: 里面有词落在 CANON 外 (A/P/C/H 单字母)
        for bad in ["RT+A+P", "MP+++C", "RT+H+A+P"]:
            self.assertFalse(comps_of(bad) <= CANON,
                             "%r 应该被判为含非法词" % bad)
            self.assertGreater(len(comps_of(bad) - CANON), 0)

    def test_fmt_comps_roundtrip(self):
        for cs in (set(), {"MP"}, {"MP", "CP"}, CANON):
            self.assertEqual(comps_of(s2.fmt_comps(cs)), cs)


class TestLocusInvariants(unittest.TestCase):
    """位点表的结构性不变量（谓词层）."""

    def _row(self, **kw):
        base = {"locus_arch": "no_host_locus", "host_scope": "none",
                "n_host_scaffolds": "0", "locus_key": "S|nohost:1",
                "locus_ncomp": "0", "locus_comps": "-", "locus_members": "1"}
        base.update(kw)
        return base

    def _check(self, r):
        errs = []
        if r["locus_arch"] == "locus_full" and int(r["n_host_scaffolds"] or 0) < 1:
            errs.append("locus_full 却没有寄主 scaffold")
        if r["locus_arch"] != "no_host_locus" and "nohost" in r["locus_key"]:
            errs.append("非 no_host_locus 的行带着 nohost 伪位点键")
        if int(r["locus_ncomp"] or 0) > len(CANON):
            errs.append("locus_ncomp 超过 5")
        if not comps_of(r["locus_comps"]) <= CANON:
            errs.append("locus_comps 含非 CANON 词")
        if int(r["locus_ncomp"] or 0) != len(comps_of(r["locus_comps"])):
            errs.append("locus_ncomp 与 locus_comps 不一致")
        return errs

    def test_good_rows_pass(self):
        self.assertEqual(self._check(self._row()), [])
        self.assertEqual(self._check(self._row(
            locus_arch="locus_full", host_scope="own_species",
            n_host_scaffolds="2", locus_key="S|scaf1",
            locus_ncomp="4", locus_comps="MP+CP+RT+RH")), [])

    def test_bad_rows_are_caught(self):
        self.assertTrue(self._check(self._row(locus_arch="locus_full")))
        self.assertTrue(self._check(self._row(locus_arch="partial",
                                             locus_key="S|nohost:3")))
        self.assertTrue(self._check(self._row(locus_ncomp="4",
                                             locus_comps="RT+A+P")))
        self.assertTrue(self._check(self._row(locus_ncomp="6",
                                             locus_comps="MP+CP+AP+RT+RH")))


class TestProvirusAndMultiInvariants(unittest.TestCase):
    """provirus_scale / multi_locus_arch 必须由它们自己的定义推出, 不许永真.

    直接调 s2 里那两个判据函数, 而不是在测试里把表达式重抄一遍 —— 重抄只能证明抄对了,
    证明不了代码里用的就是那个判据。
    """

    def test_provirus_definition(self):
        self.assertTrue(s2.provirus_ok(s2.PROVIRUS_SPAN, s2.PROVIRUS_NCOMP))
        self.assertTrue(s2.provirus_ok(s2.PROVIRUS_SPAN + 500,
                                       s2.PROVIRUS_NCOMP + 1))
        self.assertFalse(s2.provirus_ok(s2.PROVIRUS_SPAN - 1, s2.PROVIRUS_NCOMP))
        self.assertFalse(s2.provirus_ok(s2.PROVIRUS_SPAN, s2.PROVIRUS_NCOMP - 1))
        self.assertFalse(s2.provirus_ok(0, 0))

    def test_multi_locus_not_tautological(self):
        """旧定义 '组件>=1 的位点>=2' 是永真式; 新定义要求 >=2 位点且合计 >=3 组件."""
        self.assertTrue(s2.multi_locus_ok(2, {"MP", "CP", "RT"}))
        self.assertTrue(s2.multi_locus_ok(3, CANON))
        # 位点够但组件不够 -> 不成立 (这正是旧定义抓不到的情形)
        self.assertFalse(s2.multi_locus_ok(2, {"MP", "CP"}))
        self.assertFalse(s2.multi_locus_ok(5, {"MP"}))
        # 组件够但只有一个位点 -> 也不是"多基因座"
        self.assertFalse(s2.multi_locus_ok(1, CANON))
        # 空输入不得为真 (永真式的特征: 空集也为真)
        self.assertFalse(s2.multi_locus_ok(0, set()))

    def test_credit_hsps_is_bounded(self):
        """credit_hsps 不能超过实际 HSP 条数, 且 ncomp>=1 时必须 >=1."""
        hsps = [{"qs": 1, "qe": 400, "comps": "MP+CP", "bits": 100, "aln": 133},
                {"qs": 500, "qe": 900, "comps": "RT+RH", "bits": 90, "aln": 133}]
        credited, n_used = s2.credit_components(hsps)
        self.assertLessEqual(n_used, len(hsps))
        self.assertGreaterEqual(n_used, 1)
        self.assertTrue(credited)


class TestSourceLineEndings(unittest.TestCase):
    """EVE 模块的源文件必须是 LF —— 这是跨平台能不能跑的问题, 不是风格问题.

    回归: 工作区里的 CRLF 会让 `set -euo pipefail` 变成 `set -euo pipefail\\r`,
    Linux 上的 bash 判它 "invalid option name", run_all.sh 直接跑不起来。
    Windows 的 MSYS bash 容忍 CRLF, 所以本地跑得好好的, 一上服务器就炸 ——
    这个测试是本地唯一能拦住它的关卡。

    成因通常是脚本里 `read_text()` -> `write_text()` 回环: Windows 上 read_text 走
    universal newlines 把 CRLF 读成 \\n, write_text 又把 \\n 写回 CRLF, 一次回环
    就把整个文件翻过去。改文件请用能保住行尾的方式, 或改完跑本测试。
    """

    SUFFIXES = (".py", ".sh", ".tsv")
    MODULE = Path(__file__).resolve().parents[1] / "eve_distinguish"
    ROOT = Path(__file__).resolve().parents[1]

    def _files(self):
        for base in (self.MODULE, self.ROOT):
            for p in base.rglob("*"):
                if (p.is_file() and p.suffix in self.SUFFIXES
                        and ".pytest_cache" not in str(p)
                        and "__pycache__" not in str(p)):
                    yield p

    def test_no_crlf_in_module_sources(self):
        bad = [str(p) for p in self._files() if b"\r\n" in p.read_bytes()]
        self.assertEqual(bad, [], "这些文件带 CRLF, shell 脚本在 Linux 上会挂: %s" % bad)

    def test_run_all_is_posix_lf(self):
        f = self.MODULE / "run_all.sh"
        b = f.read_bytes()
        self.assertNotIn(b"\r\n", b)
        self.assertTrue(b.startswith(b"#!/bin/bash\n"), "shebang 行尾必须是 LF")


class TestS4MappingInvariants(unittest.TestCase):
    """action 只允许出现在声明集合里; 只有 host_contamination_likely 允许 REMOVE."""

    DECLARED = {"REMOVE_host_contamination", "MOVE_EVE", "KEEP_virus", "REVIEW"}

    def test_actions_declared(self):
        self.assertTrue(s4.ACTIONS, "ACTIONS 不能为空")
        for verdict, (act, reason) in s4.ACTIONS.items():
            self.assertIn(act, self.DECLARED, "%s -> %s" % (verdict, act))
            self.assertTrue(reason.strip(), "%s 缺 action_reason" % verdict)

    def test_remove_only_for_host_contamination(self):
        for verdict, (act, _) in s4.ACTIONS.items():
            if act == "REMOVE_host_contamination":
                self.assertEqual(verdict, "host_contamination_likely")

    def test_every_verdict_has_an_action(self):
        """s3 产出的每个 verdict 名都必须在 ACTIONS 里有映射 (否则按 REVIEW 兜底)."""
        produced = {"review", "EVE_LTR_TE", "EVE_STRONG_provirus", "ancient_EVE",
                    "virus_candidate", "virus_fragment_review", "EVE_suspect",
                    "assembly_breakpoint_review", "compositional_noise_review",
                    "structure_intact_review", "structure_decay_review",
                    "host_contamination_likely", "host_homology_cross_species",
                    "host_conflict_review"}
        self.assertEqual(produced - set(s4.ACTIONS), set(),
                         "有 verdict 没有 action 映射")


# ─────────────────────── 2) 真实产物层（需 EVE_VERIFY_DIR） ───────────────────────
@unittest.skipUnless(VERDICT_DIR, "未设 EVE_VERIFY_DIR, 跳过真实产物不变量检查")
class TestRealOutputsInvariants(unittest.TestCase):
    """对一份真实跑完的产物目录逐条检查不变量.

    跑法: EVE_VERIFY_DIR=/path/to/run python -m pytest tests/test_eve_invariants.py
    注意: 该目录里的派生表必须是用**当前代码**生成的。实测踩过坑: 复用上一次运行留下的
    s2_domains.tsv / locus_architecture.tsv 会得到与文档对不上的数字 (旧 s2b 的 ncomp
    偏高, EVE_STRONG_provirus 4 -> 11)。要么全量重跑, 要么从原始命中重算派生表。
    """

    @classmethod
    def setUpClass(cls):
        cls.d = Path(VERDICT_DIR)
        need = ["s1_decay.tsv", "s2_domains.tsv", "locus_architecture.tsv",
                "eve_distinguish_verdict.tsv", "dna_vs_eve_filter.tsv"]
        missing = [f for f in need if not (cls.d / f).exists()]
        if missing:
            raise unittest.SkipTest("EVE_VERIFY_DIR 缺少: %s" % ", ".join(missing))
        cls.s1 = load_tsv(cls.d / "s1_decay.tsv")
        cls.s2 = load_tsv(cls.d / "s2_domains.tsv")
        cls.la = load_tsv(cls.d / "locus_architecture.tsv")
        cls.vd = load_tsv(cls.d / "eve_distinguish_verdict.tsv")
        cls.fl = load_tsv(cls.d / "dna_vs_eve_filter.tsv")

    def test_s2_component_tokens_are_canon(self):
        bad = [r["contig_id"] for r in self.s2
               if not comps_of(r["components"]) <= CANON]
        self.assertEqual(bad, [], "components 列出现非 CANON 词: %s" % bad[:5])

    def test_s2_no_set_repr_leaked(self):
        """loci_detail 曾经把 set 直接抖进 f-string, 印出 {'MP', 'CP'}."""
        bad = [r["contig_id"] for r in self.s2
               if "{" in r["loci_detail"] or "'" in r["loci_detail"]]
        self.assertEqual(bad, [], "loci_detail 含 set repr: %s" % bad[:5])

    def test_credit_hsps_consistent(self):
        bad = [r["contig_id"] for r in self.s2
               if int(r.get("credit_hsps") or 0) > 0
               and int(r.get("best_locus_ncomp") or 0) == 0]
        self.assertEqual(bad, [], "有 credit_hsps 却 zero 组件: %s" % bad[:5])

    def test_locus_full_is_attached_to_real_locus(self):
        bad = [(r["contig_id"], r["locus_arch"], r["n_host_scaffolds"])
               for r in self.la
               if r["locus_arch"] == "locus_full"
               and int(r["n_host_scaffolds"] or 0) < 1]
        self.assertEqual(bad, [], "locus_full 没有寄主 scaffold: %s" % bad[:5])

    def test_no_pseudo_locus_pooling(self):
        bad = [(r["contig_id"], r["locus_key"]) for r in self.la
               if r["locus_arch"] != "no_host_locus" and "nohost" in r["locus_key"]]
        self.assertEqual(bad, [], "非 no_host_locus 用了伪位点键: %s" % bad[:5])

    def test_locus_ncomp_matches_comps(self):
        bad = [(r["contig_id"], r["locus_ncomp"], r["locus_comps"]) for r in self.la
               if int(r["locus_ncomp"] or 0) != len(comps_of(r["locus_comps"]))
               and r["locus_arch"] != "no_host_locus"]
        self.assertEqual(bad, [], "位点组件数与组件列不符: %s" % bad[:5])

    def test_evs_strong_only_on_locus_full(self):
        by = {r["contig_id"]: r for r in self.la}
        bad = [r["contig_id"] for r in self.vd
               if r["verdict"] == "EVE_STRONG_provirus"
               and by.get(r["contig_id"], {}).get("locus_arch") != "locus_full"]
        self.assertEqual(bad, [], "EVE_STRONG_provirus 未挂在 locus_full 上: %s" % bad[:5])

    def test_provirus_scale_column_definition(self):
        bad = [r["contig_id"] for r in self.s2
               if r["provirus_scale"] == "TRUE"
               and (int(r["best_locus_span"] or 0) < s2.PROVIRUS_SPAN
                    or int(r["best_locus_ncomp"] or 0) < s2.PROVIRUS_NCOMP)]
        self.assertEqual(bad, [], "provirus_scale 违反自身定义: %s" % bad[:5])

    def test_filter_table_matches_verdict_table(self):
        self.assertEqual(len(self.fl), len(self.vd),
                         "action 表与 verdict 表行数不一致")
        self.assertEqual({r["contig_id"] for r in self.fl},
                         {r["contig_id"] for r in self.vd})

    def test_only_host_contamination_is_removed(self):
        bad = [r["contig_id"] for r in self.fl
               if r["action"] == "REMOVE_host_contamination"
               and r["verdict"] != "host_contamination_likely"]
        self.assertEqual(bad, [], "非寄主污染被 REMOVE: %s" % bad[:5])

    def test_s1_decay_class_declared(self):
        allowed = {"coding_intact", "distributed_decay", "assembly_breakpoint",
                   "compositional_noise", "no_hsp"}
        bad = sorted({r["decay_class"] for r in self.s1} - allowed)
        self.assertEqual(bad, [], "decay_class 出现未声明取值: %s" % bad)


if __name__ == "__main__":
    unittest.main(verbosity=2)
