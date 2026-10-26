#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_eve_distinguish_logic.py — eve_distinguish 四个判定阶段的逻辑测试
=====================================================================
与 test_eve_distinguish.py 分工: 那个锁 s4 映射表与 run_all.sh 的装配; 这个锁
S1/S2/S2b/S3 里**真正参与判定的计算**, 每条都对应一次线上踩过的坑:

  s1  classify_decay 的方向 (enrichment 越大越脏) / 阈值 /
       未寄主映射候选伪造 locus_full (伪位点) /
       run_all.sh --no-host 的表头由 python 单点生成 (bash 不许再 printf 一份)
  s2  merge_loci 的 30nt 间隔 / multi_locus_arch 不再永真
  s2b 未映射候选不合并 / 同 scaffold 才合并 / host_scope 口径 / 行序可复现 /
       blastn 表列数不对要报错而不是让全部候选静默变 no_host_locus
  s3  provirus_scale 列写 TRUE/FALSE (不是 Python 的 True/False) /
       no_host_locus 不得判 EVE_STRONG_provirus / 寄主否决 reason 写清比较口径

运行: python -m unittest discover -s endogenous_virus_pipeline/tests -v
  或: python endogenous_virus_pipeline/tests/test_eve_distinguish_logic.py
"""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]                      # endogenous_virus_pipeline/
MOD = ROOT / 'eve_distinguish'
sys.path.insert(0, str(MOD))

import s1_decay_scan as s1                   # noqa: E402
import s2_domain_scan as s2                  # noqa: E402
import s2b_locus_scan as s2b                 # noqa: E402
from s1_decay_scan import S1_HDR             # noqa: E402
from s2_domain_scan import S2_HDR            # noqa: E402
from s3_verdict import VERDICT_COLS, VERDICT_HDR   # noqa: E402
from s4_filter import CARRY, OUT_HDR         # noqa: E402

PY = sys.executable

# ---------------------------------------------------------------- 表头/夹具
S2B_HDR = s2b.HDR                            # 由被测模块自己定义, 不许两头各写一份
# 上游证据表: s3 读 tax_family + aa_pident/aa_species (v6 的"上游病毒信号"两条通道)
EV_HDR = ("contig_id\ttax_family\tcategory\tcheckv_completeness\t"
          "aa_pident\taa_species\tnt_pident\tnt_species\n")

# blastn outfmt 6 的 14 列 (qlen/slen 必须有, 少一列 parse_blastn 就整行丢掉)
BN_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gaps",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]


def write(path, text, mode=0o644):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding='utf-8')
    if mode:
        os.chmod(p, mode)
    return p


def s1_row(cid, cls, stops="", other="", enrich="", head="", msf="", span="600"):
    return (f"{cid}\t600\t+1\t200\t0.5\t0.5\t10\t0\t2\t{stops}\t500\thit\t1-{span}\t"
            f"{span}\t{other}\t{enrich}\t{head}\t0.1\t{msf}\t{cls}\n")


def s2_row(cid, comps, best_n, best_comps, provirus="FALSE", multi="FALSE",
           te="FALSE", te_bits="0", cauli_bits="200", credit=""):
    # 最后一列 credit_hsps = 真贡献了组件的 HSP 条数, 缺省取组件数 (单组件单 HSP 的
    # 正常情况)。s3 不读它, 但表头列数必须对齐, 少一列整行就串位
    return (f"{cid}\t2\t2000\tL1:1-2000({comps})\t{comps}\tdetail\t{te}\t{te_bits}\t-\t"
            f"{cauli_bits}\t-\t1999\t{best_n}\t{best_comps}\t0.8\t{provirus}\t{multi}\t"
            f"{credit or best_n}\n")


def bn_row(qseqid, sseqid, pid="99.0", aln="500", qs="10", qe="509",
           bits="900", qlen="1000", slen="5000", sstart="100", send="1000"):
    d = dict(qseqid=qseqid, sseqid=sseqid, pident=pid, length=aln, mismatch="5", gaps="0",
             qstart=qs, qend=qe, sstart=sstart, send=send, evalue="1e-100",
             bitscore=bits, qlen=qlen, slen=slen)
    return "\t".join(d[c] for c in BN_COLS) + "\n"


def run_py(script, *args):
    return subprocess.run([PY, str(MOD / script), *[str(a) for a in args]],
                          capture_output=True, text=True,
                          encoding='utf-8', errors='replace')


def read_tsv_rows(path):
    with open(path, newline='') as f:
        lines = f.read().splitlines()
    hdr = lines[0].split('\t')
    return [dict(zip(hdr, ln.split('\t'))) for ln in lines[1:]]


# ══════════════════════════════════════════════════════════════════ S1
class TestS1ClassifyDecay(unittest.TestCase):
    """classify_decay: enrichment 的方向是"越大越脏", 阈值按 v5 定。"""

    def prof(self, head=0, msf=0.1):
        return {"head_stops": head, "first_stop_frac": 0.5, "max_stopfree_frac": msf}

    def test_span_too_short_is_no_hsp(self):
        self.assertEqual(s1.classify_decay(5, 1.0, self.prof(), 80), "no_hsp")
        self.assertEqual(s1.classify_decay(None, None, self.prof(), 500), "no_hsp")

    def test_zero_or_one_stop_is_coding_intact(self):
        for stops in (0, 1):
            self.assertEqual(s1.classify_decay(stops, 1.0, self.prof(head=1), 600),
                             "coding_intact")

    def test_viral_frame_much_cleaner_than_baseline_is_coding_intact(self):
        # 病毒读框 2 个终止而同链另两框均 36 个: 编码受保留, 不是"组成噪声".
        # 旧版这里一律给 compositional_noise, 把"编码完好"这条证据白丢了
        self.assertEqual(s1.classify_decay(2, 36.0, self.prof(head=0, msf=0.1), 600),
                         "coding_intact")

    def test_noisy_frames_give_compositional_noise(self):
        # 两框一样脏 (20 vs 40 均数) -> 终止来自组成/组装噪声, 无退化证据
        self.assertEqual(s1.classify_decay(20, 40.0, self.prof(head=2), 600),
                         "compositional_noise")

    def test_single_indel_gives_assembly_breakpoint(self):
        # 长无终止段占 90% -> 单点断裂, 不是化石退化
        self.assertEqual(s1.classify_decay(5, 0.5, self.prof(head=0, msf=0.9), 600),
                         "assembly_breakpoint")

    def test_many_stops_no_long_clean_run_is_distributed_decay(self):
        self.assertEqual(s1.classify_decay(12, 1.0, self.prof(head=3, msf=0.05), 600),
                         "distributed_decay")
        self.assertEqual(s1.classify_decay(3, 1.0, self.prof(head=0, msf=0.1), 600),
                         "distributed_decay")

    def test_single_head_stop_does_not_prove_distributed_decay(self):
        # 只 2 个终止, 其中 1 个落在最前 15 个密码子: 旧版判"分布式退化"(化石签名),
        # v5 要求 head_stops>=2 —— 单个头部终止同样可能只是那里恰好有个 indel
        self.assertEqual(s1.classify_decay(2, 1.0, self.prof(head=1, msf=0.1), 600),
                         "assembly_breakpoint")
        self.assertEqual(s1.classify_decay(2, 1.0, self.prof(head=0, msf=0.1), 600),
                         "assembly_breakpoint")

    def test_two_head_stops_is_distributed_decay(self):
        self.assertEqual(s1.classify_decay(4, 1.0, self.prof(head=2, msf=0.1), 600),
                         "distributed_decay")


class TestS1Helpers(unittest.TestCase):
    def test_stop_profile_counts_head_and_longest_clean_run(self):
        pos = [3, 30, 60]        # span 90: head_end = min(45, 90) = 45
        p = s1.stop_profile(pos, 90)
        self.assertEqual(p["head_stops"], 2)          # 3 和 30 都在前 45nt 内
        self.assertEqual(p["first_stop_frac"], round(3 / 90, 4))
        # 各段无终止长度: [3, 30-3-3, 60-30-3, 90-60-3] = [3, 24, 27, 27]
        # (v6.3 起相邻终止之间扣掉终止密码子自身 3nt, 末段从最后一个终止的结尾起算)
        self.assertEqual(p["max_stopfree_frac"], round(27 / 90, 4))

    def test_stop_profile_zeroes_on_stop_free_edges(self):
        # 终止紧贴两端: 首段/末段被压到 0, 不再凭空多出 3nt
        self.assertEqual(s1.stop_profile([0], 90)["max_stopfree_frac"],
                         round(87 / 90, 4))           # 只剩末段 90-0-3
        self.assertEqual(s1.stop_profile([87], 90)["max_stopfree_frac"],
                         round(87 / 90, 4))           # 只剩首段 87
        # 相邻终止 (间隔恰为 3nt): 中间段长度为 0, 不为负
        p = s1.stop_profile([0, 3, 6], 90)
        self.assertEqual(p["max_stopfree_frac"], round(81 / 90, 4))

    def test_stop_profile_empty(self):
        self.assertEqual(s1.stop_profile([], 90)["head_stops"], 0)
        self.assertEqual(s1.stop_profile([], 0)["max_stopfree_frac"], "")

    def test_frame_intervals_are_stop_to_stop(self):
        long_ = "ATG" * 120 + "TAA" + "ATG" * 120   # 唯一终止在第 360 位之后
        ivs, n = s1.frame_intervals(long_, 0)
        self.assertEqual(n, 1)
        self.assertEqual(ivs, [(0, 360), (363, len(long_))])

    def test_frame_intervals_drops_short_intervals(self):
        seq = "ATG" * 10 + "TAA" + "ATG" * 10      # 两段都只有 30nt < ORF_MIN_NT
        ivs, n = s1.frame_intervals(seq, 0)
        self.assertEqual(n, 1)
        self.assertEqual(ivs, [])

    def test_long_clean_orf_wins_over_random_backbone(self):
        # 随机序列每个读框都密布终止, 最长无终止段很短; 插进去的干净 ORF 必须被
        # 认成主读框, 且 orf_max_fraction 明显抬升
        import random
        rng = random.Random(20260923)
        base = list("".join(rng.choice("ACGT") for _ in range(1500)))
        base[996:999] = list("TAA")          # 插入点前放一个终止: ORF 正好从 ATG 起
        base = "".join(base)
        orf = "ATG" + "AAA" * 300 + "TAA"     # 903nt 干净 ORF = 301 个密码子
        plain = s1.orf_scan(base)
        with_orf = s1.orf_scan(base[:999] + orf + base[999:])
        self.assertGreater(with_orf["main_orf_aa"], 290)
        self.assertLess(plain["main_orf_aa"], 200)     # 随机序列没有长 ORF
        self.assertGreater(with_orf["orf_max_fraction"], plain["orf_max_fraction"])
        # 精确核对 +1 读框: 从 999 起必须有一段 903nt 的干净 ORF
        ivs, n_stops = s1.frame_intervals(base[:999] + orf + base[999:], 0)
        self.assertIn((999, 1902), ivs)

    # ---- P0-1 真值用例: 期望值从 blastx 坐标契约推出, 不照实现重抄公式 ----
    # 契约: qstart/qend 是首个/末个已对齐密码子的核苷酸坐标 (1 基), 所以
    #   正链锚点 = qstart-1, 负链锚点 = L-qe, 都已在主读框栅格上;
    # 双工具实测 (NCBI BLAST / DIAMOND, 2026-10-11 审查): qstart-1 ≡ |qframe|-1 (mod 3)。
    _ORF = "ATG" + "GCT" * 219          # 660nt, 不含任何终止密码子

    def _plus_probe(self, f):
        """同一条无终止编码段放进 +f 读框: 前置 f-1 个填充碱基。"""
        seq = "N" * (f - 1) + self._ORF
        hsps = [{"sid": "s1", "pident": 35.0, "qstart": f, "qend": f - 1 + len(self._ORF),
                 "evalue": 1e-50, "bits": 400.0, "frame": f"+{f}", "label": "lab"}]
        return seq, hsps

    def _minus_probe(self, k):
        """无终止编码段放进 -k 读框: revcomp 后前置 k-1 个填充碱基。

        revcomp(seq)[k-1:] 从编码段第一个密码子开始, 对应 plus 坐标
        qstart = L-(k-1), qend = qstart-659 (负链 qstart > qend)。
        """
        seq = s1.revcomp(self._ORF) + "N" * (k - 1)
        L = len(seq)
        qstart = L - (k - 1)
        hsps = [{"sid": "s1", "pident": 35.0, "qstart": qstart, "qend": qstart - 659,
                 "evalue": 1e-50, "bits": 400.0, "frame": f"-{k}", "label": "lab"}]
        return seq, hsps

    def test_hsp_metrics_six_frame_clean_orf_has_zero_stops(self):
        """P0-1 回归: 6 个读框里的同一条无终止编码段, 提前终止必须都是 0。

        旧版把 |qframe|-1 在锚点上又加了一遍: ±2/±3 框错位读框, 干净编码段被
        读出 35 个"终止"判成 distributed_decay (真值 0 / coding_intact)。
        """
        for f in (1, 2, 3):
            seq, hsps = self._plus_probe(f)
            m = s1.hsp_metrics(hsps, seq)
            self.assertEqual(m["premature_stops_region"], 0, f"+{f}")
            self.assertEqual(m["decay_class"], "coding_intact", f"+{f}")
            self.assertEqual(m["span_len"], len(self._ORF), f"+{f}")
        for k in (1, 2, 3):
            seq, hsps = self._minus_probe(k)
            m = s1.hsp_metrics(hsps, seq)
            self.assertEqual(m["premature_stops_region"], 0, f"-{k}")
            self.assertEqual(m["decay_class"], "coding_intact", f"-{k}")
            self.assertEqual(m["span_len"], len(self._ORF), f"-{k}")

    def test_hsp_metrics_stops_are_counted_on_aligned_frame_grid(self):
        """+2 读框的终止只在"1 基位置 2 起步"的栅格上数: 手工构造, 逐个可核算。

        seq = N + TAATAGTGA + A*150: +2 读框 (1 基 2,5,8,...) 的前三个密码子是
        TAA/TAG/TGA -> 恰 3 个终止; 相对锚点错后 1/2 nt 的两个基线相位无终止。
        旧版 (offset 双重应用) 会去数 +1/+3 框, 终止数读成 0。
        """
        seq = "N" + "TAATAGTGA" + "A" * 150
        hsps = [{"sid": "s1", "pident": 35.0, "qstart": 2, "qend": 160,
                 "evalue": 1e-50, "bits": 400.0, "frame": "+2", "label": "lab"}]
        m = s1.hsp_metrics(hsps, seq)
        self.assertEqual(m["premature_stops_region"], 3)
        self.assertEqual(m["stops_other_mean"], 0.0)
        self.assertEqual(m["head_stops"], 3)
        self.assertEqual(m["span_len"], 159)

    def test_hsp_metrics_minus_strand_span_and_stops(self):
        # 锁 minus 链坐标换算: span_len 必须等于 qe-qs+1; 病毒读框的终止数必须
        # 等于在 revcomp(seq)[L-qe : L-(qs-1)] 区段上、以区段起点为锚点
        # (相位 0) 独立数出来的结果; 基线 = 锚点相位 +1/+2 的两个读框。
        seq = ("ACGTTCGAAGCCTTAAGGGATCCATCGATTACGGATCCTTAAGGCCTTAAGGGCC"
               "TTAAGGCCTTAAGGGATCCATCGATTACGGATCCTTAAGGCCAATTGGCCAATT"
               "TTAAGGCCTTAAGGGATCCATCGATTACGGATCCTTAAGGCCTTAAGGGCCAAT" * 3)
        L = len(seq)
        qs, qe = 40, 320
        hsps = [{"sid": "s1", "pident": 30.0, "qstart": qs, "qend": qe,
                 "evalue": 1e-50, "bits": 400.0, "frame": "-2", "label": "lab"}]
        m = s1.hsp_metrics(hsps, seq)
        self.assertEqual(m["span_len"], qe - qs + 1)
        region = s1.revcomp(seq)[L - qe: L - (qs - 1)]
        self.assertEqual(len(region), m["span_len"])
        # 独立复算 (按定义, 不引用实现): 锚点 = 区段起点, 相位 0/1/2 各数一遍
        def stops_at(ph):
            return sum(1 for i in range(ph, len(region) - 2, 3)
                       if region[i:i + 3].upper() in s1.STOPS)
        self.assertEqual(m["premature_stops_region"], stops_at(0))
        self.assertAlmostEqual(m["stops_other_mean"],
                               round((stops_at(1) + stops_at(2)) / 2, 2), places=2)


# ══════════════════════════════════════════════════════════════════ S2
class TestS2MergeLoci(unittest.TestCase):
    def hsp(self, qs, qe, comps="RT", bits=100):
        return {"qs": qs, "qe": qe, "comps": comps, "bits": bits,
                "pident": 40.0, "aln": 200, "detail": [self.hsp.__dict__] if False else []}

    def test_adjacent_hsps_merge_within_30nt(self):
        got = s2.merge_loci([self.hsp(100, 300), self.hsp(320, 500)])  # 320 <= 300+30
        self.assertEqual(len(got), 1)
        self.assertEqual(got[0]["qs"], 100)
        self.assertEqual(got[0]["qe"], 500)

    def test_distant_hsps_stay_separate(self):
        got = s2.merge_loci([self.hsp(100, 300), self.hsp(400, 600)])
        self.assertEqual(len(got), 2)

    def test_components_union_in_merged_locus(self):
        got = s2.merge_loci([self.hsp(100, 300, "MP"), self.hsp(320, 500, "CP")])
        self.assertEqual(got[0]["comps"], "MP+CP")
        self.assertEqual(s2.locus_ncomp(got[0]), 2)

    def test_fmt_comps_is_canon_ordered_and_set_free(self):
        # 直接把 set 抖进 f-string 会印出 "{'MP'}" 这种 repr, 且 set 迭代顺序不保证,
        # 同一份输入两次跑出的 loci_detail 可能不同 —— 必须走 fmt_comps
        self.assertEqual(s2.fmt_comps({"MP", "CP"}), "MP+CP")
        self.assertEqual(s2.fmt_comps({"CP", "MP", "AP"}), "MP+CP+AP")
        self.assertEqual(s2.fmt_comps(set()), "-")
        self.assertNotIn("{", s2.fmt_comps({"MP"}))


class TestS2EndToEnd(unittest.TestCase):
    """跑真 main: 锁 provirus_scale 与 multi_locus_arch 的新定义。"""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='eve_s2_'))

    def _run(self, panel_rows):
        # 输出按 fasta 头里的 ID 走, 所以头必须是 panel 行里的 qseqid
        fa = write(self.tmp / 'q.fasta', ">c1\n" + "ATGCCCGGGTTTAAA" * 200 + "\n")
        panel = write(self.tmp / 'panel.tsv', "".join(panel_rows))
        baits = write(self.tmp / 'baits.tsv', "")
        out = self.tmp / 's2.tsv'
        r = run_py('s2_domain_scan.py', fa, panel, baits, out)
        self.assertEqual(r.returncode, 0, r.stderr)
        return read_tsv_rows(out)[0]

    def _hsp(self, qs, qe, comps, bits=500, aln=None):
        # aln 默认按 nt 跨度折算成 aa 长度: s2 的 COMP_AA_MIN 以 aa 计, 一条
        # 200 aa 的比对最多证 2 个组件 (v6 起 component credit 要独立比对长度)
        if aln is None:
            aln = max(1, (qe - qs + 1) // 3)
        # sseqid 形如 "MP+CP|Cauli|NC_x": 组件标签取 | 前 + 分隔的部分
        return (f"c1\t{comps}|panel|NC_1\t30.0\t{aln}\t{qs}\t{qe}\t100\t900\t"
                f"1e-20\t{bits}\t+1\tlabel\n")

    def test_provirus_scale_needs_span_and_ncomp(self):
        r = self._run([self._hsp(1, 1999, "MP+CP+AP+RT")])
        self.assertEqual(r['provirus_scale'], 'TRUE')
        self.assertEqual(r['best_locus_ncomp'], '4')
        r2 = self._run([self._hsp(1, 599, "MP+CP+AP+RT")])   # 跨度不够
        self.assertEqual(r2['provirus_scale'], 'FALSE')

    def test_one_short_hsp_cannot_credit_four_components(self):
        # v6 核心回归: 一条 120 aa 的比对打中融合条目 MP+CP+AP+RT, 只证 1 个组件,
        # 哪怕 nt 跨度看着够。旧版这里 ncomp=4 + span 够 -> provirus_scale=TRUE
        r = self._run([self._hsp(1, 1999, "MP+CP+AP+RT", aln=120)])
        self.assertEqual(r['best_locus_ncomp'], '1')
        self.assertEqual(r['provirus_scale'], 'FALSE')
        self.assertEqual(r['credit_hsps'], '1')
    def test_overlapping_hsps_do_not_stack_components(self):
        # 两段比对落在同一段 query 上: 不可能各自证明不同蛋白。
        # 第二条的可用区间 = 900-700 = 200nt < 300nt(COMP_AA_MIN*3) -> 不计
        r = self._run([self._hsp(100, 700, "MP", aln=200, bits=500),
                       self._hsp(300, 900, "CP", aln=200, bits=400)])
        self.assertEqual(r['best_locus_ncomp'], '1')
        self.assertEqual(r['credit_hsps'], '1')
        # 拉开到互不重叠 -> 两个组件都站得住 (两个基因座, 合计 MP+CP)
        r2 = self._run([self._hsp(100, 700, "MP", aln=200, bits=500),
                        self._hsp(1200, 1800, "CP", aln=200, bits=400)])
        self.assertEqual(r2['components'], 'MP+CP')
        self.assertEqual(r2['credit_hsps'], '2')
        self.assertEqual(r2['best_locus_ncomp'], '1')     # 单基因座仍只 1 个
        self.assertIn('L1:100-700', r2['loci_detail'])
        self.assertIn('L2:1200-1800', r2['loci_detail'])

    def test_long_single_hsp_credits_only_its_own_length(self):
        # 一条 200 aa 的比对打中融合条目 MP+CP+AP+RT: 封顶 aln//100 = 2 个组件,
        # 名字里列 4 个不顶用 (那是参考序列的事实)。v6 起 ncomp 不再等于条目里的标签数
        r = self._run([self._hsp(1, 1999, "MP+CP+AP+RT", aln=200)])
        self.assertEqual(r['best_locus_ncomp'], '2')
        self.assertEqual(r['provirus_scale'], 'FALSE')
        # 但比对本身够长 (600 aa) 时, 4 个组件是这条比对能撑起来的
        r2 = self._run([self._hsp(1, 1999, "MP+CP+AP+RT", aln=600)])
        self.assertEqual(r2['best_locus_ncomp'], '4')
        self.assertEqual(r2['provirus_scale'], 'TRUE')

    def test_available_interval_caps_components_not_all_or_nothing(self):
        """P1-4 回归: 可用区间闸门必须按'剩余空闲能放几个完整蛋白'折算.

        审查实测场景: 300 aa 融合比对 (cap=3), 但它先被一条重叠比对占掉大半,
        只剩 600 nt 空闲 —— 3 个组件需要 900 nt, 旧版 avail>=nt_min 一票放行
        照记 3 个; v7 起只能记 2 个。
        """
        # 第一条占 1-2799; 第二条 2700-3399 (300aa, cap=3),
        # 扣掉占用后空闲 = 2800-3399 = 600nt -> n_fit = 2
        hsps = [
            {"qs": 1, "qe": 2799, "comps": "MP", "bits": 500, "aln": 200},
            {"qs": 2700, "qe": 3399, "comps": "MP+CP+AP", "bits": 400, "aln": 300},
        ]
        credited, n_used = s2.credit_components(hsps)
        self.assertEqual(credited, {"MP", "CP"})       # 600nt 只放得下 2 个蛋白
        self.assertEqual(n_used, 2)
        # 对照: 空闲够 3 个蛋白 (900nt) 时, 300aa 融合比对记满 3 个
        hsps2 = [
            {"qs": 1, "qe": 2399, "comps": "MP", "bits": 500, "aln": 200},
            {"qs": 2300, "qe": 3599, "comps": "MP+CP+AP", "bits": 400, "aln": 300},
        ]
        credited2, _ = s2.credit_components(hsps2)
        self.assertEqual(credited2, {"MP", "CP", "AP"})

    def test_multi_locus_arch_not_tautological(self):
        # 单基因座 4 组件: 旧定义也给 TRUE(永真式), 新定义要求 >=2 个基因座
        r = self._run([self._hsp(1, 1999, "MP+CP+AP+RT")])
        self.assertEqual(r['multi_locus_arch'], 'FALSE')
        # 两个远隔基因座, 合计 4 个不同组件 -> TRUE (每个基因座 >=600nt, 按剩余
        # 空闲折算够放 2 个蛋白; v7 P1-4 起 491nt 的短基因座只撑得起 1 个组件)
        r2 = self._run([self._hsp(10, 1699, "MP+CP", aln=200),
                        self._hsp(3000, 4699, "AP+RT", aln=200)])
        self.assertEqual(r2['multi_locus_arch'], 'TRUE')
        self.assertEqual(r2['components'], 'MP+CP+AP+RT')
        self.assertEqual(len(r2['loci_detail'].split(';')), 2)
        # 两个远隔基因座但只有 2 个组件 -> 不够 (相邻碎片不算"多基因座架构")
        r3 = self._run([self._hsp(10, 500, "MP", aln=200),
                        self._hsp(3000, 3500, "CP", aln=200)])
        self.assertEqual(r3['multi_locus_arch'], 'FALSE')

    def test_loci_detail_is_plain_canon_ordered_text(self):
        # 回归: loci_detail 里的组件曾是 credit_components 返回的 set 原样插值,
        # 印成 `({'AP', 'MP'},bit440)`。s2b 按这列做位点分组, 不许带 repr/括号
        r = self._run([self._hsp(100, 700, "MP+CP", aln=200, bits=500),
                       self._hsp(1200, 1800, "AP", aln=200, bits=400)])
        d = r['loci_detail']
        self.assertIn("L1:100-700(MP+CP,bit500)", d)
        self.assertIn("L2:1200-1800(AP,bit400)", d)
        for bad in ("{", "}", "'", "set(", "("):
            self.assertNotIn(bad, d.replace("(MP+CP,bit", "").replace("(AP,bit", ""))
        # 与 components 列口径一致: 都是 CANON 序
        self.assertEqual(r['components'], "MP+CP+AP")


# ══════════════════════════════════════════════════════════════════ S2b
class TestS2bHelpers(unittest.TestCase):
    def test_parse_sample_clean_marker(self):
        self.assertEqual(s2b.parse_sample("ERR2040310_clean_NODE_2_length_100"),
                         ("ERR2040310", "clean"))

    def test_parse_sample_without_clean_marker(self):
        # NC_/contig_ 头: 退回整头当样本号, 打标记, 绝不猜一个前缀
        self.assertEqual(s2b.parse_sample("NC_013134.1"),
                         ("NC_013134.1", "no_clean_marker"))
        self.assertEqual(s2b.parse_sample("contig_22919"),
                         ("contig_22919", "no_clean_marker"))

    def test_scaffold_code(self):
        self.assertEqual(s2b.scaffold_code("scaffold-IAJW-2026239-Amentotaxus_argotaenia"),
                         "IAJW")
        self.assertEqual(s2b.scaffold_code("plain_id"), "plain_id")

    def test_union_len(self):
        self.assertEqual(s2b.union_len([(1, 10), (5, 20), (40, 50)]), 29)
        self.assertEqual(s2b.union_len([]), 0)


class TestS2bLocusGrouping(unittest.TestCase):
    """位点归并: 只有真位点能合并; 未寄主映射候选一条一个位点。"""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='eve_s2b_'))
        self.s2 = write(self.tmp / 's2.tsv', S2_HDR +
                        s2_row("S1_clean_1", "MP+CP", 2, "MP+CP") +
                        s2_row("S1_clean_2", "AP+RT", 2, "AP+RT") +
                        s2_row("S1_clean_3", "MP", 1, "MP") +
                        s2_row("S1_clean_4", "MP", 1, "MP") +
                        s2_row("S2_clean_1", "CP+AP", 2, "CP+AP") +
                        s2_row("NC_013134.1", "MP+CP+AP+RT+RH", 5, "MP+CP+AP+RT+RH"))
        self.s1 = write(self.tmp / 's1.tsv', S1_HDR +
                        "".join(s1_row(c, "distributed_decay") for c in
                                ("S1_clean_1", "S1_clean_2", "S1_clean_3",
                                 "S1_clean_4", "S2_clean_1", "NC_013134.1")))
        self.host_map = write(self.tmp / 'host_map.tsv',
                              "S1\tIAJW\tDirA\nS2\tGBCQ\tDirB\n")

    def _run(self, blastn_text, with_map=True):
        bl = write(self.tmp / 'locus_blastn.tsv', blastn_text)
        out = self.tmp / 'locus_arch.tsv'
        args = [bl, self.s2, self.s1, out]
        if with_map:
            args += [self.host_map]
        r = run_py('s2b_locus_scan.py', *args)
        return r, {row['contig_id']: row for row in read_tsv_rows(out)}

    def test_unmapped_candidates_are_not_pooled(self):
        # 一条 blastn 都不给: 6 条候选全部无寄主映射. 旧版会把同一样本的它们合并,
        # S1_clean_1(MP+CP) 与 S1_clean_2(AP+RT) 并起来就是 4 组件 -> locus_full
        r, rows = self._run("")
        self.assertEqual(r.returncode, 0, r.stderr)
        for cid, row in rows.items():
            self.assertEqual(row['locus_arch'], 'no_host_locus', cid)
            self.assertEqual(row['locus_members'], '1', cid)
            self.assertEqual(row['host_scope'], 'none', cid)
            self.assertIn('nohost:', row['locus_key'], cid)
        # 组件只描述自己, 不与他人合成
        self.assertEqual(rows['S1_clean_1']['locus_comps'], 'MP+CP')
        self.assertEqual(rows['NC_013134.1']['locus_comps'], 'MP+CP+AP+RT+RH')

    def test_same_scaffold_is_pooled(self):
        # 两条候选同样本、同 scaffold -> 合并成一位点, 组件并集 4 -> locus_full
        bl = (bn_row("S1_clean_1", "scaffold-IAJW-1-x") +
              bn_row("S1_clean_2", "scaffold-IAJW-1-x"))
        r, rows = self._run(bl)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(rows['S1_clean_1']['locus_arch'], 'locus_full')
        self.assertEqual(rows['S1_clean_2']['locus_arch'], 'locus_full')
        self.assertEqual(rows['S1_clean_1']['locus_members'], '2')
        self.assertEqual(rows['S1_clean_1']['locus_comps'], 'MP+CP+AP+RT')   # 按 CANON 序输出
        self.assertEqual(rows['S1_clean_1']['locus_ncomp'], '4')
        self.assertEqual(rows['S1_clean_1']['locus_hint'], 'provirus_integrated_signature')
        self.assertEqual(rows['S1_clean_1']['host_scope'], 'own_species')
        self.assertEqual(rows['S1_clean_1']['host_code'], 'IAJW')

    def test_different_scaffolds_are_not_pooled(self):
        bl = (bn_row("S1_clean_1", "scaffold-IAJW-1-x") +
              bn_row("S1_clean_2", "scaffold-IAJW-9-z"))
        r, rows = self._run(bl)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotEqual(rows['S1_clean_1']['locus_key'],
                            rows['S1_clean_2']['locus_key'])
        for cid in ('S1_clean_1', 'S1_clean_2'):
            self.assertEqual(rows[cid]['locus_arch'], 'locus_partial')
            self.assertEqual(rows[cid]['locus_members'], '1')

    def test_different_samples_are_not_pooled(self):
        bl = (bn_row("S1_clean_1", "scaffold-IAJW-1-x") +
              bn_row("S2_clean_1", "scaffold-IAJW-1-x"))     # 同一 scaffold, 不同样本
        r, rows = self._run(bl)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotEqual(rows['S1_clean_1']['locus_key'],
                            rows['S2_clean_1']['locus_key'])

    def test_no_clean_marker_candidate_gets_own_scope_and_group(self):
        bl = rc = bn_row("NC_013134.1", "scaffold-IAJW-1-x")
        r, rows = self._run(bl)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(rows['NC_013134.1']['sample_flag'], 'no_clean_marker')
        self.assertEqual(rows['NC_013134.1']['locus_members'], '1')
        # 样本号不在映射表里 -> 跨物种启发式, 不是自身物种口径
        self.assertEqual(rows['NC_013134.1']['host_scope'], 'cross_species_guess')

    def test_host_scope_own_species_requires_map_entry(self):
        bl = (bn_row("S1_clean_1", "scaffold-IAJW-1-x") +
              bn_row("S2_clean_1", "scaffold-IAJW-1-x"))
        r, rows = self._run(bl, with_map=False)   # 没给映射表: 全部跨物种启发式
        self.assertEqual(r.returncode, 0, r.stderr)
        # 只断言真有寄主命中的那些: 无命中的候选 scope 只能是 none
        for row in rows.values():
            if row['host_scaffold'] != '-':
                self.assertEqual(row['host_scope'], 'cross_species_guess', row['contig_id'])
        self.assertTrue(any(row['host_scaffold'] != '-' for row in rows.values()))

    def test_n_host_scaffolds_counts_split_hits(self):
        # 同一条候选的同源区被打断在两条 scaffold 上: 位点只取最佳那条, 但要能看见
        bl = (bn_row("S1_clean_1", "scaffold-IAJW-1-x", bits="900") +
              bn_row("S1_clean_1", "scaffold-IAJW-2-y", bits="100"))
        r, rows = self._run(bl)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(rows['S1_clean_1']['n_host_scaffolds'], '2')
        self.assertEqual(rows['S1_clean_1']['host_scaffold'], 'scaffold-IAJW-1-x')

    def test_output_is_byte_identical_across_runs(self):
        bl = "".join(bn_row(c, "scaffold-IAJW-1-x") for c in
                     ("S2_clean_1", "S1_clean_1", "S1_clean_2", "NC_013134.1"))
        r1, _ = self._run(bl)
        self.assertEqual(r1.returncode, 0, r1.stderr)
        first = (self.tmp / 'locus_arch.tsv').read_bytes()
        r2, _ = self._run(bl)
        self.assertEqual(r2.returncode, 0, r2.stderr)
        # 旧版按 set 迭代序遍历, 两次跑行序不同, 没法 diff; 现在按 contig_id 排序
        self.assertEqual((self.tmp / 'locus_arch.tsv').read_bytes(), first)

    def test_unparsable_blastn_table_is_fatal(self):
        # 行数很多但都是 12 列: 八成是 -outfmt 配错了. 必须报错, 不许让全部候选
        # 静默变成 no_host_locus 还照样往下判
        bad = write(self.tmp / 'bad.tsv',
                    "\n".join("\t".join(x for x in ("q", "s", "99", "500", "1", "500",
                                                    "1", "500", "1e-5", "900", "100", "1000"))
                              for _ in range(5)) + "\n")
        r = run_py('s2b_locus_scan.py', bad, self.s2, self.s1, self.tmp / 'o.tsv')
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('outfmt', r.stderr)

    def test_header_only_columns(self):
        r = run_py('s2b_locus_scan.py', '--emit-header')
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout, S2B_HDR)
        self.assertEqual('host_scope', S2B_HDR.split('\n')[0].split('\t')[5])
        self.assertIn('locus_key', S2B_HDR)
        # v7.4: 距离约束归并的审计列
        for col in ('subj_start', 'subj_end', 'locus_span'):
            self.assertIn(col, S2B_HDR, col)

    def test_same_scaffold_far_apart_not_pooled(self):
        """v7.4 回归: 同一 scaffold 上相隔很远的两个候选不得并成一个位点.

        旧版只按 scaffold 名分组, 于是同一 scaffold 上相隔几十 kb 的两个病毒
        片段会把组件并集凑成"齐全" (伪造 locus_full)。现在按 subject 坐标聚类,
        间隔 > LOCUS_MAX_GAP (10kb) 就是两个独立位点。
        """
        bl = (bn_row("S1_clean_1", "scaffold-IAJW-1-x", sstart="1000", send="2000") +
              bn_row("S1_clean_2", "scaffold-IAJW-1-x", sstart="90000", send="91000"))
        r, rows = self._run(bl)
        self.assertEqual(r.returncode, 0, r.stderr)
        # 两条各自成位点: 组件不合并, 都是 partial (各 2 组件)
        self.assertNotEqual(rows['S1_clean_1']['locus_key'],
                            rows['S1_clean_2']['locus_key'])
        self.assertEqual(rows['S1_clean_1']['locus_members'], '1')
        self.assertEqual(rows['S1_clean_2']['locus_members'], '1')
        self.assertEqual(rows['S1_clean_1']['locus_arch'], 'locus_partial')

    def test_same_scaffold_nearby_pooled(self):
        """同一 scaffold 上相邻 (间隔 < 10kb) 的候选仍要并成一个位点."""
        bl = (bn_row("S1_clean_1", "scaffold-IAJW-1-x", sstart="1000", send="2000") +
              bn_row("S1_clean_2", "scaffold-IAJW-1-x", sstart="3000", send="4000"))
        r, rows = self._run(bl)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(rows['S1_clean_1']['locus_key'],
                         rows['S1_clean_2']['locus_key'])
        self.assertEqual(rows['S1_clean_1']['locus_members'], '2')
        self.assertEqual(rows['S1_clean_1']['locus_arch'], 'locus_full')
        self.assertEqual(rows['S1_clean_1']['locus_span'], '3001')

    def test_locus_key_distinguishes_clusters_on_same_scaffold(self):
        """两簇在同一 scaffold 上时 locus_key 必须不同 (带簇序号)."""
        bl = (bn_row("S1_clean_1", "scaffold-IAJW-1-x", sstart="1000", send="2000") +
              bn_row("S1_clean_2", "scaffold-IAJW-1-x", sstart="2000", send="3000") +
              bn_row("S1_clean_3", "scaffold-IAJW-1-x", sstart="90000", send="91000") +
              bn_row("S1_clean_4", "scaffold-IAJW-1-x", sstart="91000", send="92000"))
        r, rows = self._run(bl)
        self.assertEqual(r.returncode, 0, r.stderr)
        k1 = rows['S1_clean_1']['locus_key']
        k3 = rows['S1_clean_3']['locus_key']
        self.assertNotEqual(k1, k3)
        self.assertTrue(k1.endswith('#0'))
        self.assertTrue(k3.endswith('#1'))


# ══════════════════════════════════════════════════════════════════ S3
def s2b_row(cid, arch, ncomp, comps, members, hint, scope="none", code="",
            scaff="-", wpid="0", cov="0", decay="", flag="clean", sample="S1",
            subj_start="", subj_end="", locus_span="0"):
    """按 s2b.HDR 的列拼一行 (列名必须与 s2b_locus_scan.HDR 完全一致).

    v7.4: 新增 subj_start / subj_end / locus_span 三列 (距离约束归并的审计列)。
    这里按列名填, 不按位置 —— HDR 加列时只有这一处要改。
    """
    hdr = S2B_HDR.split("\n")[0].split("\t")
    vals = {c: "" for c in hdr}
    vals.update({
        "contig_id": cid, "sample": sample, "sample_flag": flag,
        "host_scaffold": scaff, "host_code": code, "host_scope": scope,
        "host_wpid": wpid, "host_cov": cov, "host_xeno": "-",
        "n_host_scaffolds": "1",
        "subj_start": subj_start, "subj_end": subj_end,
        "comps_contig": comps, "best_locus_comps": comps,
        "decay_class": decay, "te_flag": "FALSE",
        "locus_key": (f"{sample}|nohost:{cid}" if scaff in ("-", "")
                      else f"{sample}|{scaff}#0"),
        "locus_arch": arch, "locus_ncomp": str(ncomp), "locus_comps": comps,
        "locus_members": str(members), "locus_span": locus_span,
        "locus_hint": hint,
    })
    return "\t".join(vals[c] for c in hdr) + "\n"


class TestS3Verdict(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix='eve_s3_'))

        def ev_row(cid, fam, aa_pid="", aa_spec="", cat="DNA_virus", cv="",
                   nt_pid="", nt_spec=""):
            return (f"{cid}\t{fam}\t{cat}\t{cv}\t{aa_pid}\t{aa_spec}\t"
                    f"{nt_pid}\t{nt_spec}\n")
        self.s1 = write(self.tmp / 's1.tsv', S1_HDR +
                        s1_row("A_distributed", "distributed_decay", stops=12, other=1,
                               enrich=6.5, head=3, msf=0.05) +
                        s1_row("A_distributed2", "distributed_decay", stops=12, other=1,
                               enrich=6.5, head=3, msf=0.05) +
                        s1_row("A_provirus_nolocus", "distributed_decay", stops=12,
                               other=1, enrich=6.5, head=3, msf=0.05) +
                        s1_row("B_hostown", "coding_intact", stops=0) +
                        s1_row("B_conflict", "coding_intact", stops=0) +
                        s1_row("B_conflict_aa", "coding_intact", stops=0) +
                        # 上游只给了核苷酸通道 (blastn 71.1% 打到病毒参考株), 科和蛋白
                        # 通道都空 —— v6.1 之前这样的候选会被判成寄主序列删掉
                        s1_row("B_ntonly", "coding_intact", stops=0) +
                        # 反向对照: nt_pident 有值但没有物种名, 不算记录
                        s1_row("B_ntnospecies", "coding_intact", stops=0) +
                        s1_row("B_coding", "coding_intact", stops=0) +
                        s1_row("C_clean_frame", "coding_intact", stops=2, other=36,
                               enrich=0.081, head=0, msf=0.1) +
                        s1_row("C_hostcross", "coding_intact", stops=2, other=36,
                               enrich=0.081, head=0, msf=0.1) +
                        s1_row("D_review", "compositional_noise", stops=20, other=40,
                               enrich=0.512, head=2, msf=0.1))
        self.s2 = write(self.tmp / 's2.tsv', S2_HDR +
                        s2_row("A_distributed", "MP+CP+AP+RT", 4, "MP+CP+AP+RT") +
                        s2_row("A_distributed2", "MP+CP+AP", 3, "MP+CP+AP") +
                        # s2 说前病毒规模, 但 s2b 说没有寄主位点: 后者赢 (v6)
                        s2_row("A_provirus_nolocus", "MP+CP+AP+RT", 4, "MP+CP+AP+RT",
                               provirus="TRUE") +
                        # 寄主假阳性没有病毒结构基因, 只有 pol 区
                        s2_row("B_hostown", "RT", 1, "RT") +
                        s2_row("B_conflict", "RT", 1, "RT") +
                        s2_row("B_conflict_aa", "RT", 1, "RT") +
                        s2_row("B_ntonly", "RT", 1, "RT") +
                        s2_row("B_ntnospecies", "RT", 1, "RT") +
                        s2_row("B_coding", "MP+CP+AP+RT", 4, "MP+CP+AP+RT") +
                        s2_row("C_clean_frame", "MP+CP", 2, "MP+CP") +
                        s2_row("C_hostcross", "RT", 1, "RT") +
                        s2_row("D_review", "RT", 1, "RT"))
        # A/ A2/ A_provirus_nolocus: 同一套组件与退化, 差别只在有没有寄主位点
        self.s2b = write(self.tmp / 's2b.tsv', S2B_HDR +
                         s2b_row("A_distributed", "locus_full", 4, "MP+CP+AP+RT", 1,
                                 "provirus_integrated_signature", scope="own_species",
                                 code="IAJW", scaff="scaffold-IAJW-1",
                                 decay="distributed_decay") +
                         s2b_row("A_distributed2", "no_host_locus", 3, "MP+CP+AP", 1,
                                 "no_host_mapping", decay="distributed_decay") +
                         s2b_row("A_provirus_nolocus", "no_host_locus", 4, "MP+CP+AP+RT", 1,
                                 "no_host_mapping", decay="distributed_decay") +
                         s2b_row("B_hostown", "no_host_locus", 1, "RT", 1,
                                 "no_host_mapping", scope="own_species", code="IAJW",
                                 scaff="scaffold-IAJW-2", wpid="99.0", cov="0.9",
                                 decay="coding_intact") +
                         # 与 B_hostown 同样的同源强度, 但上游给了病毒科: 两通道矛盾
                         s2b_row("B_conflict", "no_host_locus", 1, "RT", 1,
                                 "no_host_mapping", scope="own_species", code="IAJW",
                                 scaff="scaffold-IAJW-3", wpid="98.5", cov="0.88",
                                 decay="coding_intact") +
                         # 上游只给了蛋白级命中 (aa_pident 97.5%), 没给科: 同样矛盾
                         s2b_row("B_conflict_aa", "no_host_locus", 1, "RT", 1,
                                 "no_host_mapping", scope="own_species", code="IAJW",
                                 scaff="scaffold-IAJW-4", wpid="99.2", cov="0.91",
                                 decay="coding_intact") +
                         # 与 B_hostown / B_conflict_aa 同强度的寄主命中, 差别只在上游
                         # 记录的形态 (科 / 蛋白 / 核苷酸 / 什么都没有)
                         s2b_row("B_ntonly", "no_host_locus", 1, "RT", 1,
                                 "no_host_mapping", scope="own_species", code="IAJW",
                                 scaff="scaffold-IAJW-5", wpid="99.0", cov="0.9",
                                 decay="coding_intact") +
                         s2b_row("B_ntnospecies", "no_host_locus", 1, "RT", 1,
                                 "no_host_mapping", scope="own_species", code="IAJW",
                                 scaff="scaffold-IAJW-6", wpid="99.0", cov="0.9",
                                 decay="coding_intact") +
                         s2b_row("B_coding", "locus_full", 4, "MP+CP+AP+RT", 2,
                                 "locus_full_intact", scope="own_species", code="IAJW",
                                 scaff="scaffold-IAJW-2", decay="coding_intact") +
                         s2b_row("C_clean_frame", "no_host_locus", 2, "MP+CP", 1,
                                 "no_host_mapping", decay="coding_intact") +
                         s2b_row("C_hostcross", "no_host_locus", 1, "RT", 1,
                                 "no_host_mapping", scope="cross_species_guess",
                                 code="GBCQ", scaff="scaffold-GBCQ-1",
                                 wpid="97.0", cov="0.85", decay="coding_intact") +
                         s2b_row("D_review", "no_host_locus", 1, "RT", 1,
                                 "no_host_mapping", decay="compositional_noise"))
        # B_conflict 上游带病毒科注释; B_conflict_aa 上游没给科, 但最佳蛋白命中
        # 97.5% 一致到某个已知蛋白 —— 同样是"病毒侧有证据", 也走矛盾复核
        self.ev = write(self.tmp / 'evidence.tsv', EV_HDR +
                        ev_row("A_distributed", "Caulimoviridae", cv="60") +
                        ev_row("A_distributed2", "Caulimoviridae", cv="60") +
                        ev_row("A_provirus_nolocus", "Caulimoviridae", cv="60") +
                        ev_row("B_hostown", "") +
                        ev_row("B_conflict", "Betaflexiviridae") +
                        ev_row("B_coding", "Caulimoviridae", cv="70") +
                        ev_row("C_clean_frame", "Caulimoviridae", cv="55") +
                        ev_row("C_hostcross", "") +
                        ev_row("B_conflict_aa", "", aa_pid="97.5",
                               aa_spec="Tobacco_virus_sp.|sp|P03542.1") +
                        ev_row("D_review", "") +
                        # B_ntonly: 只有 blastn 命中, 无科无蛋白 —— 同样是病毒侧证据
                        ev_row("B_ntonly", "", nt_pid="71.1",
                               nt_spec="Begomovirus_sp.|gi|PP728250.1") +
                        # B_ntnospecies: pident 有值但物种空, 不算记录 -> 仍走 REMOVE
                        ev_row("B_ntnospecies", "", nt_pid="71.1", nt_spec=""))

    def _run(self):
        out = self.tmp / 'verdict.tsv'
        r = run_py('s3_verdict.py', self.s1, self.s2, self.s2b, self.ev, out)
        self.assertEqual(r.returncode, 0, r.stderr)
        return {row['contig_id']: row for row in read_tsv_rows(out)}

    def test_no_host_locus_cannot_be_provirus(self):
        # 核心回归: 未寄主映射候选即使组件接近齐全, 也不许判 EVE_STRONG_provirus
        rows = self._run()
        self.assertEqual(rows['A_distributed2']['verdict'], 'EVE_suspect')
        self.assertNotEqual(rows['A_distributed2']['verdict'], 'EVE_STRONG_provirus')
        self.assertEqual(rows['A_distributed']['verdict'], 'EVE_STRONG_provirus')

    def test_provirus_column_literal(self):
        rows = self._run()
        # 旧版写 Python 的 True/False, 下游按 "TRUE" 过滤会静默漏行
        for row in rows.values():
            self.assertIn(row['provirus_scale'], ('TRUE', 'FALSE'), row['contig_id'])

    def test_coding_intact_with_full_locus_is_virus_candidate(self):
        rows = self._run()
        self.assertEqual(rows['B_coding']['verdict'], 'virus_candidate')

    def test_clean_viral_frame_without_locus_is_only_fragment(self):
        # 病毒框比同链另两框干净一个数量级 -> 编码保留; 但没有寄主位点时只能给
        # virus_fragment_review, 不能升级成 virus_candidate
        rows = self._run()
        self.assertEqual(rows['C_clean_frame']['verdict'], 'virus_fragment_review')

    def test_host_veto_branches_by_evidence_scope(self):
        # v6 三条闸: 旧的"寄主同源 -> host_contamination_likely(REMOVE)" 在 1027 条真实
        # 候选上被打错 48 次, 现在按口径分流, 只有最硬的那种才允许删数据
        rows = self._run()
        # 自身物种 + 上游无病毒信号 -> 唯一能删数据的分支
        self.assertEqual(rows['B_hostown']['verdict'], 'host_contamination_likely')
        self.assertIn('自身寄主物种', rows['B_hostown']['verdict_reason'])
        self.assertIn('IAJW', rows['B_hostown']['verdict_reason'])
        # 跨物种旁证口径: 同源强度再高也不算自身寄主证据
        self.assertEqual(rows['C_hostcross']['verdict'], 'host_homology_cross_species')
        self.assertIn('跨物种旁证', rows['C_hostcross']['verdict_reason'])
        self.assertIn('GBCQ', rows['C_hostcross']['verdict_reason'])
        # 自身物种但上游给了病毒科 -> 两通道矛盾, contig 级不仲裁
        self.assertEqual(rows['B_conflict']['verdict'], 'host_conflict_review')
        self.assertIn('Betaflexiviridae', rows['B_conflict']['verdict_reason'])
        # 上游只给了蛋白级命中 (>=95% 一致到某个已知蛋白) 也一样: 没给科不等于没病毒信号
        self.assertEqual(rows['B_conflict_aa']['verdict'], 'host_conflict_review')
        self.assertIn('97.5%', rows['B_conflict_aa']['verdict_reason'])
        # v6.1: 只给了核苷酸命中 (blastn 71.1% 到某个病毒参考株) 也一样。这条专治
        # 48 条旧 REMOVE 里那 2 条只带 nt 记录的候选 (71.1/72.1% vs PP728250.1)
        self.assertEqual(rows['B_ntonly']['verdict'], 'host_conflict_review')
        self.assertIn('核苷酸命中 71.1%', rows['B_ntonly']['verdict_reason'])
        self.assertIn('PP728250.1', rows['B_ntonly']['verdict_reason'])
        # 但 pident 有值而物种名为空不算记录: 同强度寄主命中仍走删数据那条分支
        self.assertEqual(rows['B_ntnospecies']['verdict'], 'host_contamination_likely')
        # 唯一能删数据的分支: 自身物种 + 上游三条同源通道都没有记录
        self.assertEqual(rows['B_hostown']['verdict'], 'host_contamination_likely')
        # 判词只承诺它真正查过的东西, 不说成"上游无任何病毒信号"
        self.assertIn('三条同源通道', rows['B_hostown']['verdict_reason'])
        self.assertIn('均无记录', rows['B_hostown']['verdict_reason'])
        self.assertNotIn('上游无病毒科注释', rows['B_hostown']['verdict_reason'])
        # host_scope 列本身可分辨
        self.assertEqual(rows['B_hostown']['host_scope'], 'own_species')
        self.assertEqual(rows['C_hostcross']['host_scope'], 'cross_species_guess')
        self.assertEqual(rows['A_distributed']['host_scope'], 'own_species')
        self.assertEqual(rows['A_distributed2']['host_scope'], 'none')

    def test_host_veto_not_fed_by_provirus_scale(self):
        # v6: 纵向集成判定不再吃 s2 的 provirus_scale —— 那 6 条 EVE_STRONG 的 s2b 位点
        # 其实是 no_host_locus, 靠旧版 `locus_arch or provirus_scale` 混过去的
        rows = self._run()
        self.assertEqual(rows['A_provirus_nolocus']['provirus_scale'], 'TRUE')
        self.assertEqual(rows['A_provirus_nolocus']['locus_arch'], 'no_host_locus')
        self.assertEqual(rows['A_provirus_nolocus']['verdict'], 'EVE_suspect')
        self.assertNotEqual(rows['A_provirus_nolocus']['verdict'], 'EVE_STRONG_provirus')

    def test_missing_evidence_row_cannot_reach_remove(self):
        """P0-2 回归: evidence 表缺行时, 候选不得滑进唯一删数据的 REMOVE 分支.

        B_hostown 本该是 host_contamination_likely (自身物种强同源 + 上游三通道全空),
        但把它的 evidence 行抽掉后, "三通道全空"就不再是被评估出来的事实 —— 旧版
        ev.get(cid, {}) 取空照样判 REMOVE; v7 起必须 missing_evidence (REVIEW)。
        """
        kept = [ln for ln in self.ev.read_text().splitlines()[1:]
                if not ln.startswith("B_hostown")]
        # splitlines() 会剥掉行尾换行, 拼回去时必须补 (否则整个表只剩一行)
        ev = write(self.tmp / 'ev_partial.tsv', EV_HDR + "\n".join(kept) + "\n")
        out = self.tmp / 'verdict_partial.tsv'
        r = run_py('s3_verdict.py', self.s1, self.s2, self.s2b, ev, out)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {row['contig_id']: row for row in read_tsv_rows(out)}
        self.assertEqual(rows['B_hostown']['verdict'], 'missing_evidence')
        self.assertIn('从未评估', rows['B_hostown']['verdict_reason'])
        self.assertIn('P0-2', rows['B_hostown']['verdict_reason'])
        # evidence 在场的候选照常走各自分支, 不受闸门影响
        self.assertEqual(rows['B_conflict']['verdict'], 'host_conflict_review')
        self.assertEqual(rows['B_coding']['verdict'], 'virus_candidate')
        # 溯源指纹: 双方 sha256 + 覆盖账目落在 verdict 表旁
        prov = (self.tmp / 'verdict_partial.tsv.provenance.tsv').read_text()
        self.assertIn('evidence_sha256\t', prov)
        self.assertIn('query_sha256\t', prov)
        self.assertIn('missing_evidence_ids\tB_hostown', prov)
        self.assertIn('extra_evidence_ids\t', prov)
        # stdout 报账 (失败不静默)
        self.assertIn('[gate] evidence 覆盖 11/12', r.stdout)
        self.assertIn('missing_evidence', r.stdout)

    def test_panel_not_applicable_blocks_remove(self):
        """v7.1 回归 (审稿意见): 面板不适用的候选不得 REMOVE.

        同一条"自身物种强同源 + 上游三通道全空"的候选, 如果它在 Cauli 面板上
        零命中 (no_panel), "无 MP/CP/AP" 就是"没测到"而非"没有"; 旧版只加注脚
        仍判 host_contamination_likely, 用单科面板给非 Cauli 候选下寄主结论。
        """
        # D_review 在夹具里没有 MP/CP/AP (只有 RT), 且无上游病毒信号;
        # 把它改造成"自身物种 99%/0.9 + 面板零命中"的候选
        s2 = write(self.tmp / 's2_nopanel.tsv', S2_HDR +
                   s2_row("D_review", "-", 0, "-", cauli_bits="0", credit="0"))
        s2b = write(self.tmp / 's2b_nopanel.tsv', S2B_HDR +
                    s2b_row("D_review", "no_host_locus", 0, "-", 1,
                            "no_host_mapping", scope="own_species", code="IAJW",
                            scaff="scaffold-IAJW-9", wpid="99.0", cov="0.90",
                            decay="coding_intact"))
        out = self.tmp / 'verdict_nopanel.tsv'
        r = run_py('s3_verdict.py', self.s1, s2, s2b, self.ev, out)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {row['contig_id']: row for row in read_tsv_rows(out)}
        self.assertEqual(rows['D_review']['verdict'], 'host_homology_panel_na')
        self.assertIn('面板不适用', rows['D_review']['verdict_reason'])
        self.assertNotEqual(rows['D_review']['verdict'], 'host_contamination_likely')
        # 对照: 面板适用 (有组件命中) 时同强度寄主同源仍走 REMOVE 分支
        s2b_ok = write(self.tmp / 's2b_panel.tsv', S2B_HDR +
                       s2b_row("B_hostown", "no_host_locus", 1, "RT", 1,
                               "no_host_mapping", scope="own_species", code="IAJW",
                               scaff="scaffold-IAJW-2", wpid="99.0", cov="0.9",
                               decay="coding_intact"))
        out2 = self.tmp / 'verdict_panel.tsv'
        r2 = run_py('s3_verdict.py', self.s1, self.s2, s2b_ok, self.ev, out2)
        self.assertEqual(r2.returncode, 0, r2.stderr)
        rows2 = {row['contig_id']: row for row in read_tsv_rows(out2)}
        self.assertEqual(rows2['B_hostown']['verdict'], 'host_contamination_likely')

    def test_evidence_tier_column(self):
        """v7.1: evidence_tier 按可用证据的最强形态分层, 与 verdict 解耦."""
        rows = self._run()
        # A 档: s2c 侧翼证据 (本夹具没给 -G, 所以全部无 A 档)
        self.assertTrue(all(r['evidence_tier'] != 'A_genomic_flank'
                            for r in rows.values()))
        # B 档: 位点齐全 (locus_full)
        self.assertEqual(rows['A_distributed']['evidence_tier'], 'B_locus_architecture')
        self.assertEqual(rows['B_coding']['evidence_tier'], 'B_locus_architecture')
        # C 档: 只有退化通道, 没有位点
        self.assertEqual(rows['A_distributed2']['evidence_tier'], 'C_decay_only')
        # D 档: 两者都没有
        self.assertEqual(rows['D_review']['evidence_tier'], 'D_none')
        self.assertIn('evidence_tier', VERDICT_HDR)
        # 列位置: evidence_tier 必须在 verdict 之前 (s4 按列名搬, 但人读要顺)
        self.assertLess(VERDICT_COLS.index('evidence_tier'),
                        VERDICT_COLS.index('verdict'))

    def test_recon_contradiction_downgrades_to_review(self):
        """v7.3: 重建层报 0 事件而 s1 判 distributed_decay -> 转 REVIEW.

        s1 在原始序列上数终止, 会被组装 indel 的错读框污染; 重建层每个 X 都
        核对过密码子, 是独立通道。两者矛盾时不许各自坚持。
        """
        recon = write(self.tmp / 'recon.tsv',
                      "locus\tbackend\tn_frameshift\tn_stop\tn_intron\taln_aa\t"
                      "identity\tcoverage\tearliest_event_frac\t"
                      "local_flank_identity\tstatus\n"
                      "A_distributed\tminiprot+genewise\t0\t0\t0\t400\t0.9\t0.9\t\t\tok\n")
        out = self.tmp / 'v_recon.tsv'
        r = run_py('s3_verdict.py', self.s1, self.s2, self.s2b, self.ev, out,
                   "", str(recon))
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {row['contig_id']: row for row in read_tsv_rows(out)}
        self.assertEqual(rows['A_distributed']['verdict'], 'review')
        self.assertIn('重建层 0 终止', rows['A_distributed']['verdict_reason'])
        self.assertEqual(rows['A_distributed']['recon_stop'], '0')
        self.assertEqual(rows['A_distributed']['recon_status'], 'ok')

    def test_recon_confirms_decay_and_keeps_verdict(self):
        """重建层确认有事件 -> 支持退化判定, verdict 不变."""
        recon = write(self.tmp / 'recon2.tsv',
                      "locus\tbackend\tn_frameshift\tn_stop\tn_intron\taln_aa\t"
                      "identity\tcoverage\tearliest_event_frac\t"
                      "local_flank_identity\tstatus\n"
                      "A_distributed\tminiprot+genewise\t1\t8\t0\t400\t0.9\t0.9\t0.2\t0.8\tok\n")
        out = self.tmp / 'v_recon2.tsv'
        r = run_py('s3_verdict.py', self.s1, self.s2, self.s2b, self.ev, out,
                   "", str(recon))
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {row['contig_id']: row for row in read_tsv_rows(out)}
        self.assertEqual(rows['A_distributed']['verdict'], 'EVE_STRONG_provirus')
        self.assertIn('重建层确认 9 个事件', rows['A_distributed']['verdict_reason'])
        self.assertEqual(rows['A_distributed']['recon_stop'], '8')
        self.assertEqual(rows['A_distributed']['recon_frameshift'], '1')

    def test_recon_filtered_blocks_eve_move(self):
        """重建层过滤器判突变不可靠 -> EVE 类降为 REVIEW."""
        recon = write(self.tmp / 'recon3.tsv',
                      "locus\tbackend\tn_frameshift\tn_stop\tn_intron\taln_aa\t"
                      "identity\tcoverage\tearliest_event_frac\t"
                      "local_flank_identity\tstatus\n"
                      "A_distributed\tminiprot+genewise\t0\t3\t0\t400\t0.9\t0.9\t0.95\t0.2\tfiltered\n")
        out = self.tmp / 'v_recon3.tsv'
        r = run_py('s3_verdict.py', self.s1, self.s2, self.s2b, self.ev, out,
                   "", str(recon))
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {row['contig_id']: row for row in read_tsv_rows(out)}
        self.assertEqual(rows['A_distributed']['verdict'], 'review')
        self.assertIn('重建层过滤', rows['A_distributed']['verdict_reason'])
        self.assertEqual(rows['A_distributed']['recon_status'], 'filtered')


# ══════════════════════════════════════════════════════════════════ S4
class TestS4CarriesHostScope(unittest.TestCase):
    def test_host_scope_is_carried(self):
        self.assertIn('host_scope', CARRY)
        self.assertIn('host_scope', OUT_HDR)
        self.assertTrue(OUT_HDR.index('host_scope') < OUT_HDR.index('verdict'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
