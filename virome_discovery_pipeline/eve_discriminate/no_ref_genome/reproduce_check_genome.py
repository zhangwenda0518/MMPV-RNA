#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""reproduce_check.py — genome bundle 自检

从包内 eve_score_genome.py 导入特征实现（独立于训练代码），
对 sample/ 的每条序列重算分数，与 expected.tsv 逐条比对（容差 1e-6）。
"""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from eve_score_genome import base_feats, kmer_vec, read_fasta  # noqa

TOL = 1e-6


def verify_manifest():
    mf = HERE / "MANIFEST.sha256"
    if not mf.exists():
        print("[skip] 无 MANIFEST")
        return True
    bad = []
    for line in mf.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        h, rel = line.split("  ", 1)
        p = HERE / rel
        if not p.exists():
            bad.append(f"缺失 {rel}")
        elif hashlib.sha256(p.read_bytes()).hexdigest() != h:
            bad.append(f"哈希不符 {rel}")
    if bad:
        print("[FAIL] MANIFEST:")
        for b in bad:
            print("   ", b)
        return False
    print("[ok] MANIFEST 全部一致")
    return True


def auc_np(y, s):
    o = np.argsort(s, kind="mergesort")
    so, rank = s[o], np.empty(len(s))
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and so[j + 1] == so[j]:
            j += 1
        rank[o[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    return float((rank[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def main(argv=None):
    global HERE
    ap = argparse.ArgumentParser(
        description="genome bundle 自检: 用包内特征实现重算 sample/ 分数, "
                    "与 expected.tsv 逐条比对并复核 16 科 AUC")
    ap.add_argument("--dir", default=str(HERE),
                    help="bundle 目录 (含 sample/ models/ index.json)")
    a = ap.parse_args(argv)
    HERE = Path(a.dir).resolve()
    ok = verify_manifest()
    seqs = read_fasta(HERE / "sample" / "sample.fa")
    exp = list(csv.DictReader(open(HERE / "sample" / "expected.tsv",
                                   encoding="utf-8"), delimiter="\t"))
    index = {r["family"]: r for r in
             json.loads((HERE / "index.json").read_text(encoding="utf-8"))}
    models = {}
    for p in (HERE / "models").glob("*.npz"):
        z = np.load(p)
        models[p.stem] = dict(w=z["w"], mu=z["mu"], sd=z["sd"],
                              mean_A=z["mean_A"], ym=float(z["ym"][0]),
                              dim=int(z["dim"][0]))
    print(f"[load] 样本 {len(seqs)} 条, 期望 {len(exp)} 条, 模型 {len(models)} 个")

    diffs, per_fam = [], {}
    for r in exp:
        s = seqs.get(r["contig_id"])
        if s is None:
            print(f"[FAIL] 样本缺失 {r['contig_id']}")
            ok = False
            continue
        m = models[r["family"]]
        g7, _ = base_feats(s)
        X = np.concatenate([np.array(g7[1:], dtype=np.float32),
                            kmer_vec(s, 5)]).astype(np.float32)
        if len(X) != m["dim"]:
            print(f"[FAIL] {r['contig_id']}: 特征 {len(X)} ≠ 模型 {m['dim']}")
            ok = False
            continue
        A = (X - m["mu"]) / m["sd"]
        A[np.isnan(A)] = 0.0
        sc = float((A - m["mean_A"]) @ m["w"] + m["ym"])
        diffs.append(abs(sc - float(r["score"])))
        per_fam.setdefault(r["family"], []).append(
            (1 if r["label"] == "EVE" else 0, sc))

    md = max(diffs) if diffs else float("inf")
    print(f"[score] max|diff| = {md:.3e}  (容差 {TOL:.0e})")
    if md > TOL:
        ok = False

    print(f"\n{'family':<22}{'n':>6}{'AUC':>9}{'阈值':>9}")
    for fam, rs in per_fam.items():
        y = np.array([a for a, _ in rs])
        s = np.array([b for _, b in rs])
        a = auc_np(y, s) if len(set(y.tolist())) == 2 else float("nan")
        print(f"{fam:<22}{len(rs):>6}{a:>9.4f}{index[fam]['threshold']:>9.4f}")
        if not (a > 0.90):
            print(f"   [FAIL] {fam} AUC {a:.4f} <= 0.90")
            ok = False

    print("\n" + ("PASS：特征实现与模型可复现" if ok
                  else "FAIL：复现不一致，检查 FEATURE_SPEC 是否被改动"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
