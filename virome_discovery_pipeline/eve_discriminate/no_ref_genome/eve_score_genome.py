#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_score_genome.py — deploy_genome 模型打分（1,031 维，genome-primary 规格）

与 eve_score_with_deploy.py（旧 1,033 维规格）的区别：
  1. 全局特征 7 维（去掉 coding_density；orf_max_frac 只作诊断列）
  2. stop_density 取六框中"编码最完整"那一框的终止密度
  3. 二核苷酸熵按 4 字母码 16 种二核苷酸（与训练一致）
  4. k-mer = 5-mer(1024)，合计 1031 维
  5. 模型文件来自 deploy_genome/models/*.npz（16 科，阳性=完整 EVE，
     阴性=同科基因组等长窗口）

输出列：contig_id, family, species, len_nt, score, threshold, call,
        orf_max_frac(诊断), stop_density, gc
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

STOPS3 = (48, 50, 56)          # TAA / TAG / TGA（A=0 C=1 G=2 T=3）
LUT = np.full(256, 4, dtype=np.uint8)
for _i, _c in enumerate("ACGT"):
    LUT[ord(_c)] = _i
    LUT[ord(_c.lower())] = _i
CCODE = np.array([3, 2, 1, 0, 4], dtype=np.uint8)


def _codes(seq):
    return LUT[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]


def sixframe(seq):
    """(orf_max_frac, stop_density_bestframe)：6 框取编码最完整框的终止密度。"""
    c = _codes(seq)
    L = len(c)
    if L < 60:
        return 0.0, 1.0
    best_orf, best_stop, best_coding = 0, 1.0, -1.0
    for strand in (c, CCODE[c][::-1]):
        for off in range(3):
            n = (L - off) // 3
            if n < 5:
                continue
            m = strand[off:off + 3 * n].reshape(n, 3)
            code = m[:, 0].astype(np.int32) * 16 + m[:, 1] * 4 + m[:, 2]
            stops = np.isin(code, STOPS3)
            ns = int(stops.sum())
            coding = (n - ns) / n
            if coding > best_coding:
                best_coding, best_stop = coding, ns / n
            idx = np.flatnonzero(stops)
            if idx.size:
                longest = int(np.max(np.diff(np.concatenate([[-1], idx, [n]])) - 1))
            else:
                longest = n
            best_orf = max(best_orf, longest)
    return best_orf * 3 / L, best_stop


def base_feats(seq):
    """7 维全局特征。"""
    s = seq.upper()
    L = len(s)
    if L == 0:
        return [0.0] * 7, 0.0
    c = _codes(s)
    g = int((c == 2).sum())
    cc = int((c == 1).sum())
    a = int((c == 0).sum())
    t = int((c == 3).sum())
    gc = (g + cc) / L
    skew = (g - cc) / (g + cc) if (g + cc) else 0.0
    _, stop = sixframe(s)
    if L >= 6:
        m = c[0:3 * (L // 3)].reshape(-1, 3)
        atg = float(np.isin(m[:, 0] * 16 + m[:, 1] * 4 + m[:, 2], (14,)).mean()) \
            if m.shape[0] else 0.0
    else:
        atg = 0.0
    if L >= 2:
        code2 = c[:-1].astype(np.int32) * 4 + c[1:]
        cnt = np.bincount(code2[code2 < 16], minlength=16).astype(np.float64)
        pr = cnt / cnt.sum()
        pr = pr[pr > 0]
        ent = float(-(pr * np.log2(pr)).sum())
    else:
        ent = 0.0
    return [float(L), gc, skew, stop, atg, (a + t) / L, ent], 0.0


def kmer_vec(seq, k=5):
    c = _codes(seq)
    n = len(c) - k + 1
    D = 4 ** k
    if n <= 0:
        return np.zeros(D, dtype=np.float32)
    code = np.zeros(n, dtype=np.int64)
    for j in range(k):
        code = code * 4 + c[j:j + n]
    bad = (c[:n] == 4).astype(np.int64)
    for j in range(1, k):
        bad = bad | (c[j:j + n] == 4)
    code = code[~bad.astype(bool)]
    if not len(code):
        return np.zeros(D, dtype=np.float32)
    v = np.bincount(code, minlength=D).astype(np.float32)
    return v / v.sum()


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--models", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    seqs = read_fasta(a.candidates)
    meta = {}
    with open(a.meta, encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        flds = rd.fieldnames or []
        low = {c: (c or "").strip().lower() for c in flds}
        idcol = next((c for c in flds
                      if low[c] in ("contig_id", "nucleotide", "id", "seqid")),
                     None)
        famcol = next((c for c in flds if low[c] in ("family", "family_like")),
                      None)
        spccol = next((c for c in flds
                       if low[c] in ("species", "host", "host_species")), None)
        if idcol is None or famcol is None:
            raise SystemExit(f"meta 缺 ID/family 列: {flds[:8]}")
        for r in rd:
            meta[(r.get(idcol) or "").strip()] = dict(
                family=(r.get(famcol) or "").strip(),
                species=(r.get(spccol) or "").strip() if spccol else "")

    index = {r["family"]: r
             for r in json.loads(Path(a.index).read_text(encoding="utf-8"))}
    models = {}
    for p in sorted(Path(a.models).glob("*.npz")):
        z = np.load(p)
        models[p.stem] = dict(w=z["w"], mu=z["mu"], sd=z["sd"],
                              mean_A=z["mean_A"], ym=float(z["ym"][0]),
                              dim=int(z["dim"][0]))
    print(f"[load] 候选 {len(seqs)}, meta {len(meta)}, 模型 {len(models)}",
          flush=True)

    rows, skip = [], {"no_family": 0, "no_model": 0}
    for cid, s in seqs.items():
        fam = meta.get(cid, {}).get("family", "")
        if not fam:
            skip["no_family"] += 1
            continue
        m = models.get(fam)
        if m is None:
            skip["no_model"] += 1
            continue
        g7, orf_diag = base_feats(s)
        # deploy_genome 特征顺序（见 deploy_genome/train.log）：
        # gc, gc_skew, stop_density, atg_frac, at_frac, dinuc_entropy + 5-mer
        # len_nt 不进模型（长度已按设计匹配）
        X = np.concatenate([np.array(g7[1:], dtype=np.float32),
                            kmer_vec(s, 5)]).astype(np.float32)
        if len(X) != m["dim"]:
            print(f"[skip] {cid}: 特征 {len(X)} 维 ≠ 模型 {m['dim']} 维",
                  flush=True)
            continue
        A = (X - m["mu"]) / m["sd"]
        A[np.isnan(A)] = 0.0
        sc = float((A - m["mean_A"]) @ m["w"] + m["ym"])
        thr = float(index.get(fam, {}).get("threshold", 0.5))
        rows.append(dict(contig_id=cid, family=fam,
                         species=meta.get(cid, {}).get("species", ""),
                         len_nt=len(s), score=round(sc, 4),
                         threshold=round(thr, 4),
                         call="EVE_like" if sc >= thr else "virus_like",
                         orf_max_frac=round(orf_diag, 4),
                         stop_density=round(g7[3], 4), gc=round(g7[1], 4)))

    cols = ["contig_id", "family", "species", "len_nt", "score", "threshold",
            "call", "orf_max_frac", "stop_density", "gc"]
    rows.sort(key=lambda r: -r["score"])
    with open(a.out, "w", encoding="utf-8", newline="\n") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"[out] 打分 {len(rows)} 条, 跳过 {skip}", flush=True)
    by = {}
    for r in rows:
        by.setdefault(r["family"], []).append(r)
    print(f"\n{'family':<22}{'n':>6}{'n_EVE':>7}{'ratio':>8}{'median':>9}")
    for fam, rs in sorted(by.items(), key=lambda kv: -len(kv[1])):
        k = sum(1 for r in rs if r["call"] == "EVE_like")
        med = float(np.median([r["score"] for r in rs]))
        print(f"{fam:<22}{len(rs):>6}{k:>7}{k / len(rs) * 100:>7.1f}%{med:>9.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
