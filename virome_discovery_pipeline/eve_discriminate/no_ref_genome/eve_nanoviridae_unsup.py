#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_nanoviridae_unsup.py — Nanoviridae 无监督降解打分（v3.1）

v2 → v3 改了什么
----------------
1. **修正窗口不独立**。v2 把同一序列的重叠窗口当独立事件做二项检验，
   p 值偏乐观（重叠窗口高度相关）。v3 用 **长度匹配的 Monte Carlo 零分布**：
   生成与候选等长的伪序列 M 次，得到"序列级统计量"的零分布，再算 p_final。
   窗口相关性被零分布吸收。批量打分可加 `--fdr` 做 BH 校正。

2. **零分布改为单条参考 ORF 平铺**。曾试过单类 k-NN 新颖度和"随机拼接多条 ORF"，
   两者都失败：拼接混合不同基因组的 GC，零方差被撑大，已知 EVE 检出从 4/5 掉到 1/5。
   平铺保持单一基因组的组成。

3. **末端终止密码子必须去掉**。pyrodigal 的 `gene.sequence()` 带末端 TAA/TAG/TGA；
   平铺时每 ~400 bp 插一个终止，零分布自己变成"降解序列"（orf 0.09–0.17），
   所有真 EVE 都会被判成病毒。

4. **加入合成降解基准**（`nanoviridae_bench.py`）：
   按 1%/2%/5%/10% 引入终止密码子和移码，量化检出下限。

注意：`(z_stop - z_orf + |z_gc|)/sqrt(3)` 的三项等权仍是先验设定，不是学出来的。
p 值对任何固定组合都有效，但 `stat` 本身不能解释为效应量。

参考零分布
----------
pyrodigal-gv 在 193 条 Nanoviridae 参考基因组上预测的 197 条 ORF（>=250 bp）。
（v1 直接用整条基因组是错的：只有 26/193 条 orf_max_frac>=0.7。）

特征
----
  stop6 : 六框最小终止密度（每千密码子）
  orf   : 最长无终止区间占比（六框）
  gc    : GC 含量

用法
----
  python3 eve_nanoviridae_unsup.py --candidates cand.fa --out score.tsv
  python3 eve_nanoviridae_unsup.py --candidates cand.fa --out s.tsv --fdr
"""
import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

# 旧实现来自 eve_score_with_deploy.py（server 上的 /home/zhangwenda/eve_ml）。
# 收敛到本脚本内，避免对服务器绝对路径的依赖；[1]=gc, [3]=orf_max_frac 的
# 索引口径与 v3.1 实验完全一致（本脚本只消费这两维）。
def read_fasta(path):
    out, name, buf = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def nine_dim(seq):
    """旧 9 维特征（保持 v3.1 口径）。

    布局: [L, gc, gc_skew, orf_max_frac, stop_density/1000, 0, atg, at, ent]
    本脚本只用到 [1] 与 [3]，其余位保留以维持索引不变。
    """
    s = seq.upper()
    L = len(s)
    a, t = s.count("A"), s.count("T")
    g, c = s.count("G"), s.count("C")
    gc = (g + c) / L if L else 0.0
    best = 0
    if L >= 60:
        for st in (s, s.translate(_C6)[::-1]):
            for off in range(3):
                nc = (L - off) // 3
                run = br = 0
                for i in range(nc):
                    if st[off + 3 * i: off + 3 * i + 3] in STOPS:
                        br = max(br, run)
                        run = 0
                    else:
                        run += 1
                best = max(best, br, run)
    of = min(1.0, best * 3 / L) if L else 0.0
    sp = 1000.0 * sum(1 for i in range(0, L - 2, 3)
                      if s[i:i + 3] in STOPS) / max(1, L // 3)
    atg = sum(1 for i in range(0, L - 3, 3)
              if s[i:i + 3] == "ATG") / max(1, L // 3)
    ent = 0.0
    if L >= 2:
        cnt = {}
        for j in range(L - 1):
            k2 = s[j:j + 2]
            cnt[k2] = cnt.get(k2, 0) + 1
        tot = sum(cnt.values())
        ent = -sum((v / tot) * math.log2(v / tot)
                   for v in cnt.values() if v)
    return [L, gc, (g - c) / (g + c) if (g + c) else 0.0, of, sp / 100.0, 0.0,
            atg, (a + t) / L if L else 0.0, ent]


# 参考资产默认放 database/eve_ml/；也可用 MMPV_EVE_ML_DIR 指到别处
_ML_DIR = Path(os.environ.get(
    "MMPV_EVE_ML_DIR",
    str(Path(__file__).resolve().parents[3] / "database" / "eve_ml")))
ICTV = _ML_DIR / "ictv_by_family.fa"
ORF_REF = _ML_DIR / "nan_orf_ref.fa"
CACHE = _ML_DIR / "nanoviridae_ref_v3.json"
STOPS = ("TAA", "TAG", "TGA")
_C6 = str.maketrans("ACGTN", "TGCAN")
WINDOWS = (400, 250)
STEP_FRAC = 0.5          # 序列内窗口步长 w/2
K_NN = 10                # 近邻数
M_MC = 400               # Monte Carlo 重复次数
SEED = 20261022


# ---------------------------------------------------------------- 特征
def six_stop(seq):
    s = seq.upper()
    L = len(s)
    best = float("inf")
    for strand in (s, s.translate(_C6)[::-1]):
        for off in range(3):
            nc = (L - off) // 3
            if nc < 5:
                continue
            n = sum(1 for i in range(nc)
                    if strand[off + 3 * i:off + 3 * i + 3] in STOPS)
            best = min(best, 1000.0 * n / nc)
    return best if best != float("inf") else 0.0


def feat3(seq):
    f = nine_dim(seq)
    return np.array([six_stop(seq), f[3], f[1]], dtype=np.float64)


def windows_of(seq, w, step):
    return [seq[i:i + w] for i in range(0, len(seq) - w + 1, step)]


# ---------------------------------------------------------------- 参考
def ensure_orf_ref(path):
    p = Path(path)
    if p.exists():
        return
    if not Path(ICTV).is_file():
        raise SystemExit(
            "缺少 Nanoviridae 参考 ORF 库: %s; 二选一: (1) --ref 指定现成 ORF fasta; "
            "(2) 把 ictv_by_family.fa 放到 %s (或设 MMPV_EVE_ML_DIR) 后由本脚本自动构建。"
            % (p, Path(ICTV).parent))
    import pyrodigal_gv
    ref = {k: v for k, v in read_fasta(ICTV).items()
           if k.startswith("Nanoviridae|")}
    gf = pyrodigal_gv.ViralGeneFinder(meta=True)
    with open(p, "w", encoding="utf-8") as fo:
        for k, s in ref.items():
            for g in gf.find_genes(s):
                nt = g.sequence()
                if len(nt) >= 250:
                    fo.write(f">{k}|orf{g.begin}\n{nt}\n")


def build_ref(ref_path, cache_path):
    if cache_path.exists():
        z = json.loads(cache_path.read_text(encoding="utf-8"))
        return {int(k): {kk: np.array(vv) for kk, vv in v.items()}
                for k, v in z.items()}
    ensure_orf_ref(ref_path)
    seqs = read_fasta(ref_path)
    out = {}
    for w in WINDOWS:
        ws = []
        for s in seqs.values():
            ws.extend(windows_of(s, w, max(25, int(w * STEP_FRAC))))
        X = np.vstack([feat3(s) for s in ws])
        med = np.median(X, 0)
        q1, q3 = np.percentile(X, [25, 75], axis=0)
        sd = np.maximum((q3 - q1) / 1.349, np.array([0.5, 0.005, 0.005]))
        out[w] = dict(X=X, med=med, sd=sd)
        print(f"[ref] w={w}: {len(ws)} 窗口, stop p50={np.median(X[:, 0]):.1f}, "
              f"orf p50={np.median(X[:, 1]):.3f}, gc p50={np.median(X[:, 2]):.3f}",
              flush=True)
    cache_path.write_text(json.dumps(
        {str(w): {k: v.tolist() for k, v in d.items()}
         for w, d in out.items()}), encoding="utf-8")
    return out


def window_comp(seq, w, ref):
    z = (feat3(seq) - ref["med"]) / ref["sd"]
    return float((z[0] - z[1] + abs(z[2])) / np.sqrt(3.0)), z


def seq_stat(seq, w, ref):
    """序列级统计量：各窗口 z 的均值再合成。

    用均值而不是"最差窗口"：EVE 的组成/降解沿着整条序列一致，
    最差窗口会被单个测序错误或短重复带偏。
    """
    ws = windows_of(seq, w, max(25, w // 2)) or [seq[:w]]
    Z = np.vstack([(feat3(s) - ref["med"]) / ref["sd"] for s in ws])
    mz = Z.mean(0)
    comp = float((mz[0] - mz[1] + abs(mz[2])) / np.sqrt(3.0))
    per = (Z[:, 0] - Z[:, 1] + np.abs(Z[:, 2])) / np.sqrt(3.0)
    return dict(stat=comp, max_stat=float(per.max()), n_windows=len(ws),
                z=mz, per_comp=per)


# ---------------------------------------------------------------- MC 零分布
def random_tile(pool, L, rng):
    """随机取一条参考 ORF，平铺到长度 L。

    为什么不是随机拼接多条 ORF：拼接会把不同基因组的 GC 混在一起，
    零分布的 GC 方差被撑大，正好吃掉 EVE 最主要的组成信号（见 v3.0 实测：
    已知 EVE 检出从 4/5 掉到 1/5）。平铺保持单一基因组的组成。
    """
    s = pool[rng.integers(len(pool))]
    s = s.upper()
    # pyrodigal 的 gene.sequence() 带末端终止密码子；不去掉的话平铺会
    # 每 ~400 bp 周期性插入一个终止，零分布自己变成"降解序列"。
    while len(s) >= 3 and s[-3:] in STOPS:
        s = s[:-3]
    s = s[:len(s) // 3 * 3]        # 按密码子对齐
    if len(s) < 3:
        return s * (L // 3 + 1)
    if len(s) >= L:
        off = int(rng.integers(0, len(s) - L + 1))
        return s[off:off + L]
    reps = L // len(s) + 1
    return (s * reps)[:L]


def mc_null(pool, L, w, ref, m=M_MC, seed=SEED):
    rng = np.random.default_rng(seed)
    stats = []
    for _ in range(m):
        s = random_tile(pool, L, rng)
        stats.append(seq_stat(s, w, ref)["stat"])
    return np.array(stats)


# ---------------------------------------------------------------- BH
def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        prev = min(prev, p[i] * n / (rank + 1))
        q[i] = prev
    return q


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--ids", default=None)
    ap.add_argument("--ref", default=str(ORF_REF))
    ap.add_argument("--null-cache", default=str(CACHE))
    ap.add_argument("--mc", type=int, default=M_MC)
    ap.add_argument("--fdr", action="store_true",
                    help="批内 BH-FDR 校正（候选多时用）")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    refs = build_ref(Path(a.ref), Path(a.null_cache))
    pool = list(read_fasta(a.ref).values())
    seqs = read_fasta(a.candidates)
    if a.ids:
        keep = {l.strip() for l in open(a.ids, encoding="utf-8") if l.strip()}
        seqs = {k: v for k, v in seqs.items() if k in keep}
    print(f"[load] 候选 {len(seqs)} 条, 参考 ORF {len(pool)} 条", flush=True)

    rows, null_cache = [], {}
    for cid, s in seqs.items():
        L = len(s)
        w = next((x for x in WINDOWS if L >= x), None)
        if w is None:
            print(f"[skip] {cid}: {L} bp < 250", flush=True)
            continue
        if L not in null_cache:
            null_cache[L] = mc_null(pool, L, w, refs[w], m=a.mc)
        obs = seq_stat(s, w, refs[w])
        null = null_cache[L]
        p_final = (1.0 + float((null >= obs["stat"]).sum())) / (1.0 + len(null))
        p_max = (1.0 + float((null >= obs["max_stat"]).sum())) / (1.0 + len(null))
        # 稳健性：去掉最异常的 2 个窗口后结论是否改变
        per = np.sort(obs["per_comp"])
        if len(per) >= 4:
            drop2 = float(per[:-2].mean())
            loo_p = (1.0 + float((null >= drop2).sum())) / (1.0 + len(null))
        else:
            loo_p = p_final
        rows.append(dict(
            contig_id=cid, len_nt=L, window=w, n_windows=obs["n_windows"],
            stop_density=round(six_stop(s), 2),
            orf_max_frac=round(nine_dim(s)[3], 4),
            gc=round(nine_dim(s)[1], 4),
            z_stop=round(float(obs["z"][0]), 3),
            z_orf=round(float(obs["z"][1]), 3),
            z_gc=round(float(obs["z"][2]), 3),
            stat=round(obs["stat"], 4),
            p_final=float(f"{p_final:.3g}"),
            p_mosaic=float(f"{p_max:.3g}"),
            p_drop2=float(f"{loo_p:.3g}")))

    if a.fdr and rows:
        q = bh_fdr([r["p_final"] for r in rows])
        for r, qi in zip(rows, q):
            r["q_value"] = float(f"{qi:.3g}")
    else:
        for r in rows:
            r["q_value"] = ""

    for r in rows:
        key = r["q_value"] if r["q_value"] != "" else r["p_final"]
        r["call"] = ("EVE_like" if key < 0.05 else
                     "ambiguous" if key < 0.20 else "virus_like")
        if r["call"] != "EVE_like" and r["p_mosaic"] < 0.05:
            r["call"] = "mosaic_EVE"

    cols = ["contig_id", "len_nt", "window", "n_windows",
            "stop_density", "orf_max_frac", "gc",
            "z_stop", "z_orf", "z_gc", "stat",
            "p_final", "q_value", "p_mosaic", "p_drop2", "call"]
    rows.sort(key=lambda r: r["p_final"])
    with open(a.out, "w", encoding="utf-8", newline="\n") as fo:
        w_ = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w_.writeheader()
        for r in rows:
            w_.writerow({c: r.get(c, "") for c in cols})

    print(f"\n{'id':<42}{'L':>6}{'stop':>7}{'orf':>7}{'gc':>7}"
          f"{'stat':>8}{'p_final':>10}{'p_mosaic':>10}{'p_drop2':>9}  call")
    for r in rows:
        print(f"{r['contig_id'][:40]:<42}{r['len_nt']:>6}"
              f"{r['stop_density']:>7.1f}{r['orf_max_frac']:>7.3f}"
              f"{r['gc']:>7.3f}{r['stat']:>8.2f}"
              f"{r['p_final']:>10.3g}{r['p_mosaic']:>10.3g}"
              f"{r['p_drop2']:>9.3g}  {r['call']}")
    print(f"\n[done] {a.out}  (MC={a.mc}, k={K_NN})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
