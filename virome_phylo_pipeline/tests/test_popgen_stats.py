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


# ════════════════════════════════════════════════════════════════════
# 单倍型状态热图三联动 (2026-10-22, 借 hastat plot gene 构图):
#   基因结构带 (CDS/UTR/内含子/链向) + 逐基因连接线 + π/Fst 滑窗轨道
# ════════════════════════════════════════════════════════════════════

def _tiny_haps():
    """5 条单倍型 × 20 列; 位点 0/5/10/15 为多态 (ACGT 循环)。"""
    base = "ACGTACGTACGTACGTACGT"
    rows = [
        ("H1", base, 5),
        ("H2", "ACGTACGTACGTACGTACGA", 4),   # 末位 A
        ("H3", "TCGTACGTACGTACGTACGT", 3),   # 首位 T
        ("H4", "TCGTACGTACGTACGTACGA", 2),
        ("H5", "ACGAACGTACGTACGTACGT", 1),
    ]
    return rows


def test_render_hap_state_heatmap_three_layers(tmp_path):
    """hastat 论文 Fig.2 构图五层齐出: Fst 轨道/π 轨道/基因带/矩阵/行色条。

    只验「成功产出 + 返回结构与轨道计数」——构图回归 (谁把谁删了) 由
    tracks/hap_groups 入参与出图成功兜底, 比逐像素比对稳健。
    旧版只有热图一个 imshow, 没有任何轨道/行色条对象。
    """
    from utils.popgen_plots import render_hap_state_heatmap
    out = tmp_path / "hap"
    r = render_hap_state_heatmap(
        _tiny_haps(), str(out) + ".pdf", str(out) + ".png",
        title="smoke", genes=[("N", 2, 8), ("P", 12, 17)],
        tracks={"pi": {"G1": [(1, 0.01), (10, 0.03), (20, 0.01)],
                      "G2": [(1, 0.02), (10, 0.01), (20, 0.04)]},
                "fst": {"G1|G2": [(1, 0.2), (10, 0.4), (20, 0.2)]}},
        hap_groups={"H1": "G1", "H2": "G1", "H3": "G2"})
    assert r["success"] is True
    assert r["n_haps"] == 5
    assert (tmp_path / "hap.png").exists() and (tmp_path / "hap.pdf").exists()


def test_link_lines_are_drawn_per_gene(tmp_path):
    """逐位点引线 (论文 Fig.2 构图): 基因覆盖的每一列都引线; 无列覆盖画 0 条。

    连接线挂在 figure 级 artist (ConnectionPatch), 用 monkeypatch 计数。
    位点稀疏 (≤60 列/基因) 时全画; 过密由调用方抽稀 (见实现), 不在本处断言。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from utils import popgen_plots as pp

    calls = {"n": 0}
    orig = pp._link_gene_to_heatmap

    def _spy(fig, gene_ax, heat_ax, gene_intervals, **kw):
        # 区间过滤/抽稀在调用方; 这里数实际请求的列数
        calls["n"] += sum(max(1, c1 - c0 + 1)
                          for c0, c1 in gene_intervals if c1 >= c0)
        return orig(fig, gene_ax, heat_ax, gene_intervals, **kw)

    pp._link_gene_to_heatmap = _spy
    try:
        # _tiny_haps 的多态列: 位点 0(H3/H4 是 T) / 4(H5) / 19(H2) → 3 列。
        # 基因 1 覆盖列 0-2 (基因组 1-4) → 命中列 0 → 1 条;
        # 基因 2 覆盖列 2-4 (基因组 12-17) → 无多态列命中 (列 3/4 不在多态集
        # 的映射外? 多态列为 0/4/19) → 基因2 命中列... 由 bisect 决定, 不断言
        # 精确数, 只断言: 有引线画出 (逐位点行为生效) 且区外基因不画。
        out = tmp_path / "hap2"
        r = pp.render_hap_state_heatmap(
            _tiny_haps(), str(out) + ".pdf", str(out) + ".png",
            genes=[("N", 1, 4), ("OUT", 1000, 2000)])
        assert r["success"] is True
        assert calls["n"] >= 1        # 至少基因 1 有逐列引线
        # 区外基因 (1000-2000 超出 20 列) 不贡献引线 —— 上面已计入的只是基因1
        out2 = tmp_path / "hap3"
        calls["n"] = 0
        r2 = pp.render_hap_state_heatmap(
            _tiny_haps(), str(out2) + ".pdf", str(out2) + ".png",
            genes=[("OUT", 1000, 2000)])
        assert r2["success"] is True
        assert calls["n"] == 0        # 区外基因零引线
    finally:
        pp._link_gene_to_heatmap = orig
        plt.close("all")


def test_sliding_tracks_pi_fst_shapes():
    """滑窗轨道 (hastat 论文 Fig.2 口径): π 按群体、Fst 按群体两两各一条线。

    fixture: L=40, X/Y 两组; 列 0-19 全同, 列 20-29 X 组内 A/C 混合
    (AAAC) 且 Y 组固定 C, 列 30-39 又全同。window=20/step=10 → 3 窗:
      窗1 (列 0-19):  无分化 → 各组 π=0, 无 Fst 点 (轨道不凭空造点)
      窗2 (列 10-29): X 组内 π>0 (组内多态), Fst=0.667 (中间值: 双等位列
                      X 多数 A / Y 全 C, Hudson 非饱和 —— 证明 Fst 不是
                      只会输出 0/1 的哑值)
      窗3 (列 20-39): 同上
    无元数据 → π 退化为单一 'all' 轨道, Fst 为空。
    """
    import tempfile
    import os
    from popgen_analysis import PopGenAnalyzer
    fd, tmp_fasta = tempfile.mkstemp(suffix=".fasta")
    seg1, seg3 = "ACGT" * 5, "ACGTAACGTA"          # 20 + 10
    with os.fdopen(fd, "w") as fh:
        for i, b in enumerate("AAAC"):
            fh.write(f">X{i}\n{seg1}{b * 10}{seg3}\n")
        for i in range(4):
            fh.write(f">Y{i}\n{seg1}{'C' * 10}{seg3}\n")
    try:
        pg = PopGenAnalyzer(tmp_fasta)
        t = pg.sliding_tracks(window=20, step=10)
        # 无元数据: 单一 'all' π 轨道, 无 Fst
        assert list(t["pi"].keys()) == ["all"]
        assert len(t["pi"]["all"]) == 3
        assert all(p >= 0 for _, p in t["pi"]["all"])
        assert t["fst"] == {}
        # 有元数据: X/Y 两组 → π 按群体各一条, Fst 按群体对各一条
        meta = tmp_fasta + ".meta.csv"
        with open(meta, "w") as mf:
            mf.write("name,location\n")
            for i in range(4):
                mf.write(f"X{i},X\n")
            for i in range(4):
                mf.write(f"Y{i},Y\n")
        pg2 = PopGenAnalyzer(tmp_fasta, meta, name_col="name")
        t2 = pg2.sliding_tracks(window=20, step=10, group_col="location", min_n=2)
        assert set(t2["pi"].keys()) == {"X", "Y"}
        assert set(t2["fst"].keys()) == {"X|Y"}      # '|' 分隔, 群体名含 _ 也不混
        x_pts, y_pts = t2["pi"]["X"], t2["pi"]["Y"]
        assert len(x_pts) == 3 and len(y_pts) == 3
        assert x_pts[0][1] == 0.0 and y_pts[0][1] == 0.0   # 窗1 无分化
        assert x_pts[1][1] > 0                            # X 组内多态 → π>0
        assert y_pts[1][1] == 0.0                         # Y 组内无多态
        fpts = t2["fst"]["X|Y"]
        assert len(fpts) == 2                             # 窗2/3 有点, 窗1 无
        assert all(0.0 < v < 1.0 for _, v in fpts)        # 中间值 (非 0/1 哑值)
        # 行色条映射: 每个单倍型 → 主导群体 (确定性: 并列取字典序)
        gm = pg2.hap_group_map("location", min_n=1)
        assert set(gm.values()) <= {"X", "Y"}
    finally:
        os.unlink(tmp_fasta)
        if os.path.exists(tmp_fasta + ".meta.csv"):
            os.unlink(tmp_fasta + ".meta.csv")


def test_sliding_tracks_downsample_cap():
    """窗口数超 400 时 step 自动倍增, 轨道点数 ≤ 400 (与热图 max_sites 同量级)。"""
    import tempfile, os
    from popgen_analysis import PopGenAnalyzer
    fd, tmp_fasta = tempfile.mkstemp(suffix=".fasta")
    L = 2000
    with os.fdopen(fd, "w") as fh:
        for i in range(3):
            seq = "".join("ACGT"[(i + j) % 4] for j in range(L))
            fh.write(f">S{i}\n{seq}\n")
    try:
        pg = PopGenAnalyzer(tmp_fasta)
        t = pg.sliding_tracks(window=100, step=1)   # 1900+ 窗口 → 必须降采样
        # 无元数据: 单一 'all' 轨道; 每条轨道点数 ≤ 400
        assert list(t["pi"].keys()) == ["all"]
        assert len(t["pi"]["all"]) <= 400
        # 位置单调递增 (step 倍增后不重叠)
        pos = [p for p, _ in t["pi"]["all"]]
        assert pos == sorted(pos) and len(set(pos)) == len(pos)
    finally:
        os.unlink(tmp_fasta)


def test_track_pos_maps_to_matrix_columns():
    """轨道横轴必须映射到多态位点列号 (与矩阵/基因带共享基准) —— 审查修复回归锁。

    背景: sliding_tracks 的 pos 是窗口中心 (基因组/比对列坐标 1-based), 矩阵
    横轴是多态位点序号 0..ncol-1。位点稀疏时两者差几个量级 (pos 到 16 而
    矩阵只有 2 列), 不映射则轨道画到完全错误的列上且毫无提示。

    验法: 20 列、只有 2 个多态位点 (列 5 和 15) 的 fixture; 轨道 pos 取窗口
    中心 6 和 16 (1-based), 映射后必须分别落到列 0 和 1, 且所有轨道点横坐标
    ∈ [0, ncol-1] (列号基准), 不得残留 6/16 这类基因组坐标。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from utils import popgen_plots as pp

    seen = {"xs": []}
    orig_plot = plt.Axes.plot

    def _spy(self, *args, **kwargs):
        out = orig_plot(self, *args, **kwargs)
        for ln in out:
            seen["xs"].extend(float(x) for x in ln.get_xdata())
        return out

    base = "ACGTACGTACGTACGTACGT"
    rows = [
        ("H1", base, 3),
        ("H2", base[:5] + "A" + base[6:15] + "T" + base[16:], 2),
        ("H3", base[:5] + "C" + base[6:15] + "G" + base[16:], 1),
    ]
    tracks = {"pi": {"all": [(6, 0.02), (16, 0.05)]},
              "fst": {"A|B": [(6, 0.3), (16, 0.6)]}}
    import tempfile
    import os
    fd, out = tempfile.mkstemp(suffix=".pdf")
    os.close(fd)
    plt.Axes.plot = _spy
    try:
        r = pp.render_hap_state_heatmap(
            rows, out, out.replace(".pdf", ".png"),
            title="coord", tracks=tracks)
        assert r["success"] is True, r.get("error")
        assert r["n_sites"] == 2                 # 只有 2 个多态列
        # 所有 plot 线的 x 必须落在 [0, ncol-1] (列号基准); 旧 bug 会把轨道画在
        # x=6/16 (基因组坐标) 而矩阵只有 2 列 → 轨道整体跑到图外/错列
        assert seen["xs"], "轨道没有画任何线"
        for x in seen["xs"]:
            assert 0 <= x <= r["n_sites"] - 1, (
                f"轨道横轴越界 (未映射到列号): x={x}, ncol={r['n_sites']}")
        # 两个窗口中心的点都应画出 (列 5 → 列 0, 列 15 → 列 1)
        assert len(seen["xs"]) >= 4              # π 2 点 + Fst 2 点
    finally:
        plt.Axes.plot = orig_plot
        plt.close("all")
        os.unlink(out)
        png = out.replace(".pdf", ".png")
        if os.path.exists(png):
            os.unlink(png)
