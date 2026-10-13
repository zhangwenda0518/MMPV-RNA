# -*- coding: utf-8 -*-
"""⑤ 群体遗传统计的数值对拍测试 (2026-10-13, 修 10-10 P0-3/P0-5/P0-6)。

原则: 统计公式一律调用 scikit-allel / 仓库既有 bh_fdr, 本文件只验证
「实现 == 库的参考值」与「旧 bug 场景不再复现」, 不重新实现任何公式。
"""
import json
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_PHYLO = _HERE.parent
_REPO = _PHYLO.parent
for p in (str(_PHYLO), str(_REPO)):
    if p not in sys.path:
        sys.path.insert(0, p)

# ⑤ 的 utils/ 与其他管线 (③④) 的同名包冲突: 全量 pytest 时若其他测试先导入了
# 别的 utils, 这里必须清缓存并把 ⑤ 目录插到最前, 保证 utils.dataset_config
# 解析到 ⑤ 自己的包。
for _m in [m for m in list(sys.modules) if m == "utils" or m.startswith("utils.")]:
    del sys.modules[_m]
sys.path.insert(0, str(_PHYLO))

from popgen_analysis import PopGenAnalyzer, _ac_matrix  # noqa: E402
from selection_analysis import parse_fel_meme  # noqa: E402


def _write_inputs(tmp_path, seqs, pops):
    fa = tmp_path / "aln.fasta"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()), encoding="ascii")
    meta = tmp_path / "meta.csv"
    meta.write_text("name,location\n"
                    + "".join(f"{k},{v}\n" for k, v in pops.items()), encoding="ascii")
    return str(fa), str(meta)


# 组 A 的 A_gap 在两个多态位点上都是 gap → 触发旧 fst_one「剔 gap 后标签错位」的场景
SEQS_GAP = {
    "A1": "AAAAAACCGG",
    "A2": "AAAAAACGCG",
    "A_gap": "AAAAAA--GG",
    "B1": "TTTTTTCCGG",
    "B2": "TTTTTTCGCG",
    "B2b": "TTTTTTCCGG",
}
POPS = {"A1": "A", "A2": "A", "A_gap": "A", "B1": "B", "B2": "B", "B2b": "B"}


def test_tajima_d_matches_scikit_allel(tmp_path):
    """P0-3: Tajima's D 必须等于 scikit-allel 参考值 (绝对量纲标准式)。

    旧手写实现分子混用每位点 pi 与绝对 theta_w → 符号可翻转;
    现实现直接调 allel.tajima_d, 以库输出为对照。
    """
    import allel

    fa, meta = _write_inputs(tmp_path, SEQS_GAP, POPS)
    pg = PopGenAnalyzer(fa, meta, name_col="name")
    d = pg.diversity()

    ac = _ac_matrix(list(SEQS_GAP.values()))
    ref = float(allel.stats.diversity.tajima_d(ac, min_sites=1))
    assert d["tajima_D"] == pytest.approx(round(ref, 4), abs=1e-9)
    # 反回归自检: 在本 fixture 上, 旧式分子 (每位点 pi − θW) 与新式分子
    # (k = pi·L − θW) 符号相反 —— 即旧实现会翻转 Tajima's D 的解读方向。
    old_numerator = d["pi"] - d["theta_w"]
    new_numerator = d["pi"] * d["L"] - d["theta_w"]
    assert old_numerator * new_numerator < 0


def test_hd_bounded_and_sequences_used_unified(tmp_path):
    """P0-4 回归锁: Hd 恒在 [0,1], n 与单倍型频数同源。"""
    fa, meta = _write_inputs(tmp_path, SEQS_GAP, POPS)
    pg = PopGenAnalyzer(fa, meta, name_col="name")
    d = pg.diversity()
    assert 0.0 <= d["Hd"] <= 1.0
    assert d["n"] == len(pg.sequences_used)
    assert d["haplotypes"] >= 1


def test_fst_gap_label_sync(tmp_path):
    """P0-5: Fst 管线输出 == 测试内独立构组计数矩阵直调 allel 的参考值。

    旧 fst_one 列内剔 gap 后群体标签未同步 → 含 gap 位点标签错位。
    新实现按样本计数矩阵聚合; 这里在测试侧**独立**重建各组计数
    (不经 _ac_matrix/PopGen 内部路径), 直调 allel.hudson_fst 对拍。
    """
    from collections import Counter

    import allel

    fa, meta = _write_inputs(tmp_path, SEQS_GAP, POPS)
    pg = PopGenAnalyzer(fa, meta, name_col="name")
    f = pg.fst("location", min_n=1, n_perm=99)
    assert f["fst_site"] is not None and -0.5 <= f["fst_site"] <= 1.0

    L = len(SEQS_GAP["A1"])
    alleles = "ACGT"
    group_ac = {}
    for g in ("A", "B"):
        ac = [[0] * 4 for _ in range(L)]
        for k, s in SEQS_GAP.items():
            if POPS[k] != g:
                continue
            for j, ch in enumerate(s):
                if ch in alleles:
                    ac[j][alleles.index(ch)] += 1
        group_ac[g] = ac
    num, den = allel.stats.fst.hudson_fst(group_ac["A"], group_ac["B"])
    ok = [i for i, d in enumerate(den) if d == d and d > 0]
    ref = sum(num[i] for i in ok) / sum(den[i] for i in ok)
    assert f["fst_site"] == pytest.approx(ref, abs=5e-4)

    fp = pg.fst_pairwise("location", min_n=1)
    m = fp["matrix"]
    assert m[0][0] == 0.0 and m[1][1] == 0.0
    assert m[0][1] == m[1][0]
    # 全局 (两组) 与两两一致 —— 自洽性
    assert m[0][1] == pytest.approx(f["fst_site"], abs=5e-4)


def test_fel_meme_q_gates_positivity(tmp_path):
    """P0-6: 阳性判定必须由 BH q < fdr 门控, --fdr 必须真实控制结果。

    BH 的性质是 q ≤ p, 所以收紧 --fdr 才能区分: fdr=0.01 时,
    4 个 p=0.001 (q=0.005) 阳性; 16 个 p=0.04 (p<0.05 但 q=0.04 ≥ 0.01)
    旧实现 (只看 p<0.05) 全部误判阳性, 新实现必须排除。
    """
    n_sites = 20
    pvals = [0.001] * 4 + [0.04] * (n_sites - 4)
    # HyPhy JSON 的 content 行是数组 (与 headers 按位 zip), 不是 dict
    rows = [[p, 2.0, 1.0, 10.0] for p in pvals]
    payload = {"MLE": {"headers": [["p-value"], ["beta"], ["alpha"], ["LRT"]],
                       "content": {"part": rows}}}
    jf = tmp_path / "g.FEL.json"
    jf.write_text(json.dumps(payload), encoding="utf-8")

    r = parse_fel_meme(str(jf), "g", fdr=0.01, method="FEL")
    assert r["n_positive"] == 4
    assert all(h["q_value"] < 0.01 for h in r["hits"])
    assert all(h["beta"] > h["alpha"] for h in r["hits"])

    # MEME 分支: 同样 q 门控 (无 beta>alpha 附加条件)
    r2 = parse_fel_meme(str(jf), "g", fdr=0.01, method="MEME")
    assert r2["n_positive"] == 4
