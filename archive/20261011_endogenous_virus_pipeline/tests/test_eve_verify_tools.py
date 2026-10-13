#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_eve_verify_tools.py — **校准层自身的**测试
===============================================
`verify/` 下的工具有 900+ 行, 它们负责给判别侧"定罪"或"洗白"。这个文件管的是:
**谁来看住这些工具？**

四类测试, 对应校准层实测踩过的四类失效:

1. 已知答案对照 (`TestNullModelControl`)
   仪器瞎 → 报出完全虚假的数字。oracle 初版就是这么报出"88% 未获证实"的。
   做法: 喂答案由构造保证的输入, 断言仪器分得开。
     * 已知编码 (frame 1 无终止)  -> 判据必须看得见
     * 纯随机                   -> 判据不该声称富集

2. 干预生效性 (`TestSensitivityIntervention`)
   把"我没改成"报成"这个参数无影响"。`COMP_AA_MIN` 初版就是这么报出 0% 的
   (默认参数在 def 时求值, 改全局到不了代码路径)。
   做法: 检查**必须随之移动的下游量**, 而不是被改的变量本身。

3. 零模型稳健性 (`TestNullModelShufflers`)
   结论只在一种零模型下成立就写进论文。做法: 三种洗牌器逐条验性质
   (保长度 / 保组成 / 保读框 / 保双核苷酸), 且都不得复用同一份碱基池以外的东西。

4. 独立实现差分 (`TestBruteForceDifferential`)
   把"我读过代码、看着对"换成"两个独立实现互相印证"。
   暴力版用与线上完全不同的机制重算同一件事, 在真实输入上断言逐条一致。

跑法:
  python -m pytest tests/test_eve_verify_tools.py -q            # 离线部分
  EVE_VERIFY_RUN_DIR=/path/to/run python -m pytest tests/test_eve_verify_tools.py -q
       # 追加真实数据上的差分与多零模型 (需要 q.fa + panel_hits.tsv)
"""
import collections
import os
import random
import sys
import unittest
from pathlib import Path

TESTS = Path(__file__).resolve().parent
PKG = TESTS.parent
sys.path.insert(0, str(PKG))
sys.path.insert(0, str(PKG / "eve_distinguish"))
# verify/ 的位置随打包方式而变: 仓库里是 <pkg>/eve_distinguish/verify/, 服务器上验证时
# 常把测试与工具一起拍平在 <dir>/ 与 <dir>/verify/。
# 顺序要紧: `insert(0, ...)` 是后插的优先, 所以 ***/verify 放在最后插 —— 工具的
# `MOD = HERE/..` 是相对**工具自己所在目录**算的, 万一旁边还留着一份拍平的旧副本,
# 让它抢先就会指到上一层去找阶段脚本, 直接 FileNotFoundError (这个坑踩过)。
for _cand in (PKG, PKG / "eve_distinguish", TESTS,
              PKG / "eve_distinguish" / "verify", PKG / "verify", TESTS / "verify"):
    if _cand.is_dir():
        sys.path.insert(0, str(_cand))

import null_model as nm            # noqa: E402
import sensitivity as sens         # noqa: E402
import s2_domain_scan as s2        # noqa: E402
import s1_decay_scan as s1         # noqa: E402

RUN_DIR = os.environ.get("EVE_VERIFY_RUN_DIR", "")


# ─────────── 1) 已知答案对照: 仪器必须看得见已知信号 ───────────
class TestNullModelControl(unittest.TestCase):
    """`null_model control` 的已知答案对照必须通过 —— 它是判据结论的前置闸门."""

    @classmethod
    def setUpClass(cls):
        import tempfile
        cls.tmp = tempfile.mkdtemp(prefix="nmctl_")
        cls.rc = nm.control(cls.tmp, n=40, k=10)

    def test_control_passes(self):
        self.assertEqual(self.rc, 0, "null_model 通不过自己的已知答案对照")

    def test_sees_known_coding_signal(self):
        """已知无终止的 ORF 必须被报成 coding_intact (否则仪器瞎)."""
        import tempfile
        d = tempfile.mkdtemp(prefix="nmctl_orf_")
        rd = nm._make_case(d, "orf", 40, random.Random(7))
        agg = nm.region_null(rd, 10, 7, "", null_mode="mono", quiet=True)
        self.assertGreaterEqual(agg["real_frac"]["coding_intact"], 0.5,
                                "已知编码序列没被认出来")
        nf = agg["null_frac"]["coding_intact"]
        self.assertTrue(nf == 0 or agg["real_frac"]["coding_intact"] / nf >= 3.0,
                        "已知编码的相对富集不显著")

    def test_random_has_no_enrichment(self):
        """纯随机不该被声称有编码富集."""
        import tempfile
        d = tempfile.mkdtemp(prefix="nmctl_rnd_")
        rd = nm._make_case(d, "random", 40, random.Random(11))
        agg = nm.region_null(rd, 10, 11, "", null_mode="mono", quiet=True)
        self.assertLessEqual(agg["real_frac"]["coding_intact"], 0.2,
                             "纯随机序列被判出了编码保留")


# ─────────── 2) 干预生效性: 看下游量, 不看变量本身 ───────────
class TestSensitivityIntervention(unittest.TestCase):
    """注入必须改变**下游诊断量** —— 否则扫描会静默给出"该参数无影响"的假结论.

    回归: 初版自检只断言"模块全局被改成了新值", 它通过了, 而结论依然错 ——
    `credit_components(hsps, comp_aa_min=COMP_AA_MIN)` 的默认参数在 def 时求值。
    """

    def test_override_reaches_def_time_default(self):
        """`=CONST` 形式的默认参数也必须被覆盖 (负向后顾不能吃掉 `default=CONST`)."""
        m_small = sens.load_module("s2_domain_scan.py", {"COMP_AA_MIN": 10})
        m_big = sens.load_module("s2_domain_scan.py", {"COMP_AA_MIN": 1000})
        # 同一条 HSP: 短比对在小门槛下能记功, 在大门槛下不能
        h = [{"qs": 1, "qe": 400, "comps": "MP+CP", "bits": 100, "aln": 120}]
        self.assertTrue(m_small.credit_components(h)[0],
                        "COMP_AA_MIN=10 时 120aa 的比对应该能记功")
        self.assertFalse(m_big.credit_components(h)[0],
                         "COMP_AA_MIN=1000 时同一条比对不该记功 —— 默认参数没被覆盖")

    def test_diagnostic_moves_when_param_is_binding(self):
        """诊断量定义必须对参数敏感: provirus TRUE 数随 PROVIRUS_NCOMP 单调不增."""
        import tempfile
        tmp = tempfile.mkdtemp(prefix="diag_")
        rows = [{"provirus_scale": "TRUE", "best_locus_span": "2000",
                 "best_locus_ncomp": "5"},
                {"provirus_scale": "TRUE", "best_locus_span": "1600",
                 "best_locus_ncomp": "4"},
                {"provirus_scale": "FALSE", "best_locus_span": "100",
                 "best_locus_ncomp": "1"}]
        f = Path(tmp) / "s2.tsv"
        hdr = "contig_id\tbest_locus_span\tbest_locus_ncomp\tprovirus_scale\n"
        f.write_text(hdr + "".join(
            "%s\t%s\t%s\t%s\n" % (i, r["best_locus_span"], r["best_locus_ncomp"],
                                  r["provirus_scale"]) for i, r in enumerate(rows)),
            encoding="utf-8")
        self.assertEqual(sens._diag_s2(str(f), "provirus_n"), 2)

    def test_control_flags_unreachable_value(self):
        """注入 NCOMP>组件数 后 provirus_ok 必须恒假 —— control 的判据基础.

        不能用默认常量的模块测这件事: 默认 NCOMP=4 时 ncomp=6 当然 >=4。必须用
        **注入过的**模块 —— 这正是 control 里"不可能值"能成立的前提。
        """
        m = sens.load_module("s2_domain_scan.py", {"PROVIRUS_NCOMP": 6})
        self.assertFalse(m.provirus_ok(10 ** 9, 5),
                         "NCOMP=6 时 5 个组件不该算达标")
        self.assertFalse(m.provirus_ok(10 ** 9, 4))
        self.assertTrue(m.provirus_ok(10 ** 9, 6))
        self.assertFalse(m.provirus_ok(m.PROVIRUS_SPAN - 1, 6))

    def test_blind_params_are_declared(self):
        """只在寄主分支起作用的参数必须被显式标成不可测, 不能报"无影响"."""
        self.assertIn("AA_VIRAL_PID", sens.BLIND_NO_HOST)
        stage, diag = sens.DIAGNOSTICS["AA_VIRAL_PID"]
        self.assertEqual((stage, diag), ("s3", "host_branch"))


# ─────────── 3) 三种零模型洗牌器的性质 ───────────
class TestNullModelShufflers(unittest.TestCase):
    """三种零模型的**性质**必须各自成立, 否则"零模型无关"这句话就不成立."""

    def setUp(self):
        self.rng = random.Random(2026)
        self.seq = "".join(self.rng.choice("ACGT") for _ in range(1200))

    def test_all_preserve_length_and_composition(self):
        for m in nm.NULL_MODES:
            out = nm.shuffler(self.seq, self.rng, m)
            self.assertEqual(len(out), len(self.seq), m)
            self.assertEqual(sorted(out), sorted(self.seq), "%s 改变了碱基组成" % m)

    def test_frame_preserves_per_position_composition(self):
        """frame 洗牌必须保住**位置特异**组成 —— 这是它比 mono 更严的全部理由."""
        out = nm.shuffle_frame(self.seq, self.rng)
        for off in range(3):
            self.assertEqual(sorted(out[off::3]), sorted(self.seq[off::3]),
                             "第 %d 密码子位的组成变了" % (off + 1))

    def test_frame_actually_reorders(self):
        out = nm.shuffle_frame(self.seq, self.rng)
        self.assertNotEqual(out, self.seq)

    def test_di_preserves_dinucleotide_composition(self):
        """di 洗牌必须保住**全部二核苷酸频率** (欧拉路径法)."""
        out = nm.shuffle_di(self.seq, self.rng)
        self.assertEqual(len(out), len(self.seq))
        self.assertEqual(collections.Counter(zip(out, out[1:])),
                         collections.Counter(zip(self.seq, self.seq[1:])))

    def test_di_falls_back_on_short_input(self):
        """太短的序列退回 mono, 但绝不能返回长度不对的东西."""
        for s in ("A", "AC", "ACG"):
            out = nm.shuffle_di(s, self.rng)
            self.assertEqual(len(out), len(s))

    def test_modes_are_distinct(self):
        """三种洗牌得到的结果不该全同 —— 否则"多零模型"是假的."""
        outs = {m: nm.shuffler(self.seq, random.Random(5), m) for m in nm.NULL_MODES}
        self.assertGreater(len(set(outs.values())), 1)


# ─────────── 4) 独立实现差分 ───────────
def merge_loci_bruteforce(hsps, gap=30):
    """暴力版位点合并: 反复把任意一对可合区间并掉, 直到不动点.

    与线上版机制不同 —— 线上版是"排序后单趟扩展", 这里是"反复扫描到收敛",
    不共用任何区间算术, 所以两者一致才说明结果不依赖实现方式。
    """
    ivs = sorted((h["qs"], h["qe"]) for h in hsps)
    changed = True
    while changed:
        changed = False
        out = []
        for s, e in ivs:
            if out and s <= out[-1][1] + gap:
                if e > out[-1][1]:
                    out[-1] = (out[-1][0], e)
                changed = True
            else:
                out.append((s, e))
        ivs = out
    return ivs


def credit_components_bruteforce(hsps, comp_aa_min=100):
    """暴力版组件记功: 用**被占用的核苷酸位置集合**判断可用区间, 不做区间减法.

    线上版用 `subtract_ivs` 做区间算术 + 求和; 这里把占用展开成位置集合再取差集。
    机制完全不同, 因而是一个真正的独立参照。
    """
    nt_min = comp_aa_min * 3
    credited, occupied, n_used = set(), set(), 0
    for h in sorted(hsps, key=lambda x: (-x["bits"], x["qs"])):
        avail = len(set(range(h["qs"], h["qe"] + 1)) - occupied)
        labels = [c for c in s2.CANON if c in h["comps"].split("+")]
        cap = max(1, h["aln"] // comp_aa_min)
        got = labels[:cap] if avail >= nt_min else []
        if got:
            n_used += 1
            credited.update(got)
        occupied |= set(range(h["qs"], h["qe"] + 1))
    return credited, n_used


class TestBruteForceDifferential(unittest.TestCase):
    """线上实现 vs 暴力实现, 在**随机构造的输入**上逐例比对 (离线)."""

    def _rand_hsps(self, rng, n):
        out = []
        for _ in range(n):
            qs = rng.randint(1, 4000)
            aln = rng.choice([40, 60, 99, 100, 101, 200, 300, 450])
            qe = qs + aln * 3 - 1
            comps = "+".join(rng.sample(s2.CANON, rng.randint(1, 5)))
            out.append({"qs": qs, "qe": qe, "comps": comps,
                        "bits": rng.randint(40, 600), "aln": aln})
        return out

    def test_merge_loci_matches_bruteforce(self):
        rng = random.Random(4242)
        for _ in range(300):
            hsps = self._rand_hsps(rng, rng.randint(1, 9))
            got = [(L["qs"], L["qe"]) for L in s2.merge_loci(hsps)]
            want = merge_loci_bruteforce(hsps)
            self.assertEqual(got, want, "合并结果与暴力版不一致: %r" % hsps)

    def test_credit_components_matches_bruteforce(self):
        rng = random.Random(99)
        for _ in range(300):
            hsps = self._rand_hsps(rng, rng.randint(1, 8))
            g_comp, g_n = s2.credit_components(hsps)
            w_comp, w_n = credit_components_bruteforce(hsps)
            self.assertEqual(g_comp, w_comp, "组件集合不一致: %r" % hsps)
            self.assertEqual(g_n, w_n, "贡献 HSP 条数不一致: %r" % hsps)

    def test_credit_components_respects_threshold(self):
        """跨门槛行为: 恰好 nt_min 可用则记功, 差一个碱基则不记."""
        nt = s2.COMP_NT_MIN
        hsps = [{"qs": 1, "qe": nt, "comps": "MP", "bits": 100, "aln": 400}]
        self.assertTrue(s2.credit_components(hsps)[0], "恰好够长应记功")
        hsps2 = [{"qs": 1, "qe": nt - 1, "comps": "MP", "bits": 100, "aln": 400}]
        self.assertFalse(s2.credit_components(hsps2)[0], "差一个碱基不该记功")


class TestClassifyDecaySpecDifferential(unittest.TestCase):
    """`classify_decay` 的 docstring 就是规格, 枚举网格逐点比对实现与规格.

    规格 (见 s1_decay_scan.classify_decay docstring):
      span < MIN_SPAN                     -> no_hsp
      stops <= 1                          -> coding_intact
      enrich < 1:
          stops <= CLEAN_MAX_STOPS 且 enrich <= CLEAN_ENRICH -> coding_intact
          否则                                              -> compositional_noise
      enrich >= 1:
          msf >= 0.25                     -> assembly_breakpoint
          head >= MIN_HEAD_STOPS 或 stops >= MIN_DECAY_STOPS -> distributed_decay
          否则                             -> assembly_breakpoint
    """

    @staticmethod
    def spec(stops, other_mean, head, msf, span):
        if stops is None or span < s1.MIN_SPAN:
            return "no_hsp"
        if stops <= 1:
            return "coding_intact"
        enrich = (stops + 1) / (other_mean + 1)
        if enrich < 1.0:
            if stops <= s1.CLEAN_MAX_STOPS and enrich <= s1.CLEAN_ENRICH:
                return "coding_intact"
            return "compositional_noise"
        if msf != "" and msf >= 0.25:
            return "assembly_breakpoint"
        if head >= s1.MIN_HEAD_STOPS or stops >= s1.MIN_DECAY_STOPS:
            return "distributed_decay"
        return "assembly_breakpoint"

    def test_grid_matches_spec(self):
        grid_stops = [0, 1, 2, 3, 5, 8, 20]
        grid_other = [0.0, 1.0, 2.0, 5.0, 20.0]
        grid_head = [0, 1, 2, 3]
        grid_msf = ["", 0.05, 0.24, 0.25, 0.9]
        grid_span = [0, 89, 90, 300]
        n = 0
        for stops in grid_stops:
            for other in grid_other:
                for head in grid_head:
                    for msf in grid_msf:
                        for span in grid_span:
                            prof = {"head_stops": head, "max_stopfree_frac": msf,
                                    "first_stop_frac": 0}
                            got = s1.classify_decay(stops, other, prof, span)
                            want = self.spec(stops, other, head, msf, span)
                            self.assertEqual(got, want,
                                             "stops=%s other=%s head=%s msf=%s span=%s"
                                             % (stops, other, head, msf, span))
                            n += 1
        self.assertGreater(n, 1000, "网格太小, 覆盖不足")


# ─────────── 真实数据上的差分 (需 EVE_VERIFY_RUN_DIR) ───────────
@unittest.skipUnless(RUN_DIR, "未设 EVE_VERIFY_RUN_DIR, 跳过真实数据差分")
class TestRealDataDifferential(unittest.TestCase):
    """在真实候选上跑同样的差分与多零模型 —— 随机输入覆盖不到真实 HSP 分布."""

    @classmethod
    def setUpClass(cls):
        seqs = s1.read_fasta(os.path.join(RUN_DIR, "q.fa"))
        hits = nm.read_hits(os.path.join(RUN_DIR, "panel_hits.tsv"))
        cls.cases = []
        for cid, seq in seqs.items():
            hs = hits.get(cid, [])
            if not hs:
                continue
            raw = []
            for h in hs:
                if h["length"] < s2.MIN_ALN or h["pident"] < s2.MIN_PID:
                    continue
                comps = "+".join(c for c in h["label"].split("|")[0].split("+")
                                 if c in s2.CANON)
                if not comps:
                    continue
                lo, hi = min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])
                raw.append({"qs": lo, "qe": hi, "comps": comps,
                            "bits": h["bits"], "aln": h["length"]})
            if raw:
                cls.cases.append((cid, raw))

    def test_bruteforce_agrees_on_real_candidates(self):
        self.assertGreater(len(self.cases), 100, "可比对的候选太少")
        for cid, raw in self.cases:
            self.assertEqual([(L["qs"], L["qe"]) for L in s2.merge_loci(raw)],
                             merge_loci_bruteforce(raw), cid)
            self.assertEqual(s2.credit_components(raw),
                             credit_components_bruteforce(raw), cid)

    def test_multi_null_is_run_instrumented(self):
        """多零模型在真实数据上要能跑出三组数 (数值见 verify/README.md)."""
        for m in nm.NULL_MODES:
            agg = nm.region_null(RUN_DIR, 5, 20260924, "", null_mode=m, quiet=True)
            self.assertIn("coding_intact", agg["real_frac"])
            self.assertEqual(agg["null_mode"], m)
            self.assertGreater(agg["n_with_region"], 100,
                               "%s 零模型下可比对候选太少" % m)


if __name__ == "__main__":
    unittest.main(verbosity=2)
