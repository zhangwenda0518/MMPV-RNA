#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_dbitscore.py — Brait et al. 的 Δbit-score 双库竞争法，接入本管道

方法（EVE_threshold.R, Brait et al. 2023 bioRxiv 556789）：
  每条 contig 对两个蛋白库各做一次 blastx，取最优 hit：
    EXO_BITSCORE = 对"外源活病毒蛋白库"的最优 bitscore
    EVE_BITSCORE = 对"已知 EVE 参考库"的最优 bitscore
  ΔBITSCORE = EXO_BITSCORE − EVE_BITSCORE
  Δ > 0 且 max 分高 → 真病毒；Δ < 0 → EVE 转录本

本管道的库：
  EXO 库 = ICTV 六框翻译（reclass/ictv_sixframe.faa，6.86M 蛋白）
  EVE 库 = 41,933 簇共识库的六框翻译（prot_check/our_eve_sixframe.faa）
  基准评测时 EVE 库需剔除测试 EVE 自己的簇（防自命中泄漏）；
  应用模式用全库（与 Brait 原文一致——他们的 EVE 参考也来自同一批蚊子）。

用法：
  diamond blastx 需先建好两个 .dmnd（本脚本调 diamond makedb/blastx）
"""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np


def run(cmd, log=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[cmd fail] {cmd}\n{r.stderr[:500]}", flush=True)
        raise SystemExit(1)
    if log:
        Path(log).write_text(r.stdout, encoding="utf-8")
    return r.stdout


def read_fasta(p):
    out, name, buf = {}, None, []
    for line in open(p, encoding="utf-8", errors="replace"):
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


def diamond_db(faa, dmnd, threads, diamond):
    if Path(dmnd).exists():
        return
    run(f'"{diamond}" makedb --in {faa} -d {dmnd} --threads {threads}')


def best_bits(query_fa, dmnd, threads, diamond, out_tsv, sensitivity="--very-sensitive"):
    run(f'"{diamond}" blastx -d {dmnd} -q {query_fa} -o {out_tsv} '
        f'--outfmt 6 qseqid bitscore pident length --max-target-seqs 1 '
        f'{sensitivity} --threads {threads} --masking 0')
    best = {}
    for line in open(out_tsv, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        q, b = p[0], float(p[1])
        if q not in best or b > best[q]:
            best[q] = b
    return best


def auc_np(y, s):
    y, s = np.asarray(y), np.asarray(s, float)
    if len(set(y.tolist())) < 2:
        return float("nan")
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries", required=True, help="候选核酸 FASTA")
    ap.add_argument("--labels", default=None,
                    help="TSV: seq_id/label/family（有标签时算 AUC）")
    ap.add_argument("--label-col", default="label")
    ap.add_argument("--id-col", default="seq_id")
    ap.add_argument("--fam-col", default="family")
    ap.add_argument("--exo-dmnd", required=True)
    ap.add_argument("--eve-faa", required=True)
    ap.add_argument("--eve-dmnd", required=True)
    ap.add_argument("--exclude-clusters", default=None,
                    help="从 EVE 库剔除的簇 ID 文件（每行一个，如 C000001）")
    ap.add_argument("--work", default="dbitscore_work")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--diamond", default="diamond")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    eve_faa, eve_dmnd = Path(a.eve_faa), Path(a.eve_dmnd)

    # 1) 按需剔除自命中簇
    if a.exclude_clusters:
        excl = {l.strip() for l in open(a.exclude_clusters, encoding="utf-8")
                if l.strip()}
        filtered = work / "eve_db_filtered.faa"
        n_in = n_out = 0
        with open(eve_faa, encoding="utf-8", errors="replace") as fi, \
                open(filtered, "w", encoding="utf-8") as fo:
            keep = False
            for line in fi:
                if line.startswith(">"):
                    cid = line[1:].strip().split("|")[0]
                    keep = cid not in excl
                    n_in += 1
                    n_out += not keep
                if keep:
                    fo.write(line)
        print(f"[eve-db] 剔除 {n_out} 条自命中簇蛋白，"
              f"保留 {n_in - n_out}/{n_in}", flush=True)
        eve_faa = filtered
        eve_dmnd = work / "eve_ho.dmnd"

    diamond = a.diamond
    diamond_db(a.eve_faa if not a.exclude_clusters else eve_faa, eve_dmnd,
               a.threads, diamond)
    if not Path(a.exo_dmnd).exists() and Path(a.exo_dmnd).with_suffix("").exists():
        pass
    # EXO dmnd 若不存在则从 fasta 建（--exo-faa 由调用者保证）
    if not Path(a.exo_dmnd).exists():
        faa = Path(a.exo_dmnd).with_suffix(".faa")
        if faa.exists():
            diamond_db(faa, a.exo_dmnd, a.threads, diamond)
        else:
            raise SystemExit(f"EXO 库不存在: {a.exo_dmnd}")

    # 2) 双 blastx
    print("[blastx] vs EXO 库 ...", flush=True)
    exo = best_bits(a.queries, a.exo_dmnd, a.threads, diamond,
                    work / "exo_hits.tsv")
    print(f"  {len(exo)} 条有 EXO 命中", flush=True)
    print("[blastx] vs EVE 库 ...", flush=True)
    eve = best_bits(a.queries, eve_dmnd, a.threads, diamond,
                    work / "eve_hits.tsv")
    print(f"  {len(eve)} 条有 EVE 命中", flush=True)

    # 3) 合并
    rows = []
    qids = list(read_fasta(a.queries))
    labels = {}
    if a.labels:
        for r in csv.DictReader(open(a.labels, encoding="utf-8"),
                                delimiter="\t"):
            labels[r[a.id_col]] = r
    for q in qids:
        e = exo.get(q, 0.0)
        v = eve.get(q, 0.0)
        row = dict(contig_id=q, exo_bitscore=round(e, 1),
                   eve_bitscore=round(v, 1),
                   delta_bitscore=round(e - v, 1),
                   max_bitscore=round(max(e, v), 1))
        if q in labels:
            lr = labels[q]
            row["label"] = lr.get(a.label_col, "")
            row["family"] = lr.get(a.fam_col, "")
        rows.append(row)

    cols = ["contig_id", "label", "family", "exo_bitscore", "eve_bitscore",
            "delta_bitscore", "max_bitscore"]
    with open(a.out, "w", encoding="utf-8", newline="\n") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})

    # 4) 汇总
    out = {}
    if a.labels and rows and rows[0].get("label"):
        by = {}
        for r in rows:
            by.setdefault((r.get("family"), r.get("label")), []).append(r)
        for (fam, lab), rs in sorted(by.items()):
            d = [r["delta_bitscore"] for r in rs]
            print(f"  {fam:<20}{lab:<14} n={len(rs):<5} "
                  f"Δ中位={np.median(d):>8.1f}  Δ>0: "
                  f"{sum(1 for x in d if x > 0) / len(d) * 100:.1f}%", flush=True)
        for fam in sorted({r.get("family") for r in rows if r.get("family")}):
            rs = [r for r in rows if r.get("family") == fam
                  and r.get("label") in ("EVE", "VIRUS")]
            if len(rs) < 20:
                continue
            y = [1 if r["label"] == "EVE" else 0 for r in rs]
            for key, colname in (("delta", "delta_bitscore"),
                                 ("max", "max_bitscore")):
                s = [r[colname] for r in rs]
                a_ = auc_np(y, s)
                out.setdefault(fam, {})[f"auc_{key}"] = round(a_, 4)
                print(f"  [AUC] {fam} {key}: {a_:.4f}", flush=True)
        y_all = [1 if r.get("label") == "EVE" else 0 for r in rows
                 if r.get("label") in ("EVE", "VIRUS")]
        s_all = [r["delta_bitscore"] for r in rows
                 if r.get("label") in ("EVE", "VIRUS")]
        if len(set(y_all)) == 2:
            a_ = auc_np(y_all, s_all)
            out["overall"] = {"auc_delta": round(a_, 4)}
            print(f"  [AUC] overall Δbit-score: {a_:.4f}", flush=True)
    json.dump(out, open(a.out.replace(".tsv", ".json"), "w"), indent=2)
    print(f"[done] {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
