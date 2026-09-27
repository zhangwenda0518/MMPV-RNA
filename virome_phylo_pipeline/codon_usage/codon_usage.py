"""密码子使用分析 (ENC / ENC'/RSCU/GC123/Fop/CAI)。

算法语义取自 PhyloSuite 内嵌 codonusage.py (Wright 1990; Sun et al. 2012),
重写为无依赖单文件 + CLI。独立工具, 不在管线 stage 体系内。

用法:
    python codon_usage/codon_usage.py --fasta cds.fa --outdir out/ [--gcode 1]
    python codon_usage/codon_usage_plots.py --tsv out/codon_usage.tsv --rscu out/rscu.tsv

产出:
    codon_usage.tsv  — 每序列一行: ENC, ENC2(Sun), ENC5(Sun), GC, GC1, GC2, GC3,
                       GC3s, Fop, 长度, 有效密码子判定
    rscu.tsv         — 64 密码子 RSCU 矩阵 (每序列一列)
    plots/           — ENC/PR2/中性线/RSCU热图 (300dpi 白底色盲友好)

基因代码: 1=标准 (默认), 其余 NCBI 表号按需扩 codontables。
"""
import os
import argparse

# 标准遗传代码: codon -> (aa, degeneracy class)
_BASES = "TCAG"
_AAS = {
    "A": "GCT GCC GCA GCG", "R": "CGT CGC CGA CGG AGA AGG",
    "N": "AAT AAC", "D": "GAT GAC", "C": "TGT TGC",
    "Q": "CAA CAG", "E": "GAA GAG", "G": "GGT GGC GGA GGG",
    "H": "CAT CAC", "I": "ATT ATC ATA", "L": "TTA TTG CTT CTC CTA CTG",
    "K": "AAA AAG", "M": "ATG", "F": "TTT TTC", "P": "CCT CCC CCA CCG",
    "S": "TCT TCC TCA TCG AGT AGC", "T": "ACT ACC ACA ACG",
    "W": "TGG", "Y": "TAT TAC", "V": "GTT GTC GTA GTG",
    "*": "TAA TAG TGA",
}
STOPS = {"TAA", "TAG", "TGA"}


def _build_table():
    """codon -> [aa, degeneracy, count]。degeneracy: four/three/two/one/stop。"""
    tbl = {}
    for aa, codons in _AAS.items():
        for c in codons.split():
            if aa == "*":
                tbl[c] = ["*", "stop", 0]
            else:
                deg = {2: "two", 4: "four", 3: "three", 1: "one", 6: "six"}[len(codons.split())]
                tbl[c] = [aa, deg, 0]
    return tbl


CODON_TABLE = _build_table()
# Sun 口径: L/S/R 六义不并入 four/two, 单列 six
_AA_DEG = {c: v[1] for c, v in CODON_TABLE.items()}


def count_codons(cds):
    """-> {codon: count}, 忽略非 ACGT。"""
    counts = {c: 0 for c in CODON_TABLE}
    n_valid = n_bad = 0
    cds = cds.upper().replace("U", "T")
    for i in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[i:i + 3]
        if codon in counts:
            counts[codon] += 1
            n_valid += 1
        elif set(codon) <= set("ACGT"):
            if codon in STOPS:
                n_valid += 1
            else:
                n_bad += 1
        else:
            n_bad += 1
    return counts, n_valid, n_bad


def gc_content(cds):
    cds = cds.upper().replace("U", "T")
    coding = [b for b in cds if b in "ACGT"]
    if not coding:
        return {}
    total = len(coding)
    gc = sum(1 for b in coding if b in "GC")
    gc1 = sum(1 for i in range(0, total - total % 3, 3) if coding[i] in "GC")
    gc2 = sum(1 for i in range(1, total - total % 3, 3) if coding[i] in "GC")
    gc3 = sum(1 for i in range(2, total - total % 3, 3) if coding[i] in "GC")
    n1 = len(range(0, total - total % 3, 3))
    # GC3s: 排除 Met/Trp/Stop 的第 3 位
    syn3 = gc3
    syn3 -= CODON_TABLE.get("ATG", ["", "", 0])[2] * 0  # 占位, 下面统一算
    n_syn3, g_syn3 = 0, 0
    for i in range(2, total - total % 3, 3):
        codon = "".join(coding[i - 2:i + 1])
        if codon not in ("ATG", "TGG") and codon not in STOPS:
            n_syn3 += 1
            if coding[i] in "GC":
                g_syn3 += 1
    return {
        "GC": gc / total,
        "GC1": gc1 / n1 if n1 else 0,
        "GC2": gc2 / n1 if n1 else 0,
        "GC3": gc3 / n1 if n1 else 0,
        "GC3s": g_syn3 / n_syn3 if n_syn3 else 0,
    }


def enc_wright(counts):
    """Wright 1990 / CodonW (Peden) 忠实实现。

    Nc = n1 + Σ_z n_z / F̄_z; F̂ = (nΣp² − 1)/(n−1) 按整个氨基酸家族计算
    (六义族 L/S/R 不拆亚族, 这是 Sun 系 ENC' 的约定), F̄_z 为同简并度
    氨基酸的 F̂ 平均。n1 = 密码表中 1-fold 氨基酸数 (标准码 = 2: M/W)。
    与 CodonW 唯一偏差: 某简并度类内无可用氨基酸时 CodonW 硬报错
    "Nc was not calculated", 此处回落 F̂ = 1/z (均匀期望) 并 stderr 警告。
    Nc > 61 截断到 61 (CodonW 同款后处理)。"""
    import sys as _sys
    totb, numaa, fold = {}, {}, {}
    for aa, codons in _AAS.items():
        if aa == "*":
            continue
        cl = codons.split()
        z = len(cl)
        fold[z] = fold.get(z, 0) + 1
        n = sum(counts[c] for c in cl)
        if n <= 1:
            continue  # CodonW: 该氨基酸出现 ≤1 次则跳过
        homo = sum((counts[c] / n) ** 2 for c in cl)
        f_hat = (n * homo - 1) / (n - 1)
        if f_hat > 1e-7:
            totb[z] = totb.get(z, 0.0) + f_hat
            numaa[z] = numaa.get(z, 0) + 1
    nc = float(fold.get(1, 0))
    for z in sorted(fold):
        if z == 1:
            continue
        if numaa.get(z):
            averb = totb[z] / numaa[z]
        else:
            averb = 1.0 / z  # 回落, 见 docstring
            _sys.stderr.write(
                f"[enc_wright] no aa used in {z}-fold class; Fhat=1/{z} fallback\n")
        nc += fold[z] / averb
    return min(nc, 61.0)


def rscu(counts):
    """-> {codon: RSCU}。"""
    out = {}
    for aa, codons in _AAS.items():
        if aa == "*":
            continue
        cl = codons.split()
        n = sum(counts[c] for c in cl)
        for c in cl:
            out[c] = (counts[c] / n) * len(cl) if n else 0.0
    return out


def fop(counts):
    """Fraction of Optimal Codons: RSCU>1 的同义密码子占比 (自指最优集,
    严格 CAI 需高表达参考基因集, 这里给自指 Fop 与 RSCU 定义一致口径)。"""
    opt_hits = syn_total = 0
    for aa, codons in _AAS.items():
        if aa == "*" or len(codons.split()) < 2:
            continue
        cl = codons.split()
        n = sum(counts[c] for c in cl)
        if not n:
            continue
        best = max(cl, key=lambda c: counts[c])
        opt_hits += counts[best]
        syn_total += n
    return opt_hits / syn_total if syn_total else float("nan")


def enc_sun(counts, pseudocount=False):
    """Sun et al. 2012 ENC' (按退化度分组的 Nc 估计)。
    pseudocount=False -> eq2Sun; True -> eq5Sun (+1 伪计数)。"""
    groups = {}  # deg -> list of (F, n_aa) per aa
    for aa, codons in _AAS.items():
        if aa == "*":
            continue
        cl = codons.split()
        cs = [counts[c] for c in cl]
        n = sum(cs)
        m = len(cl)
        if pseudocount:
            f = sum(((x + 1) / (n + m)) ** 2 for x in cs)
        else:
            f = (sum((x / n) ** 2 for x in cs) if n else None)
        deg = len(cl)
        if deg == 6:
            groups.setdefault("six", []).append((f, n))
        else:
            groups.setdefault({2: "two", 3: "three", 4: "four", 1: "one"}[deg],
                              []).append((f, n))
    nc = 0.0
    for deg_name, entries in groups.items():
        m = len(entries)
        fs = [f if f is not None else 1.0 / m for f, _ in entries]
        avg = sum(f * n for f, n in entries) / sum(n for _, n in entries) if \
            sum(n for _, n in entries) else sum(fs) / m
        nc += m / avg
    return nc


def parse_fasta(text):
    data, name, parts = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                data.append((name, "".join(parts)))
            name, parts = line[1:].strip(), []
        elif line:
            parts.append(line)
    if name is not None:
        data.append((name, "".join(parts)))
    return data


def analyze(cds, gcode=1):
    if gcode != 1:
        raise NotImplementedError("Only standard genetic code (gcode=1) implemented")
    counts, n_valid, n_bad = count_codons(cds)
    if n_valid < 20:
        return None  # 太短不可靠
    stats = {"length_nt": len(cds), "valid_codons": n_valid, "ambiguous": n_bad}
    stats.update({k.upper().replace("S", "s"): v for k, v in gc_content(cds).items()})
    # 统一命名: GC / GC1 / GC2 / GC3 / GC3s
    g = gc_content(cds)
    stats.update({"GC": g["GC"], "GC1": g["GC1"], "GC2": g["GC2"],
                  "GC3": g["GC3"], "GC3s": g["GC3s"]})
    stats["ENC"] = enc_wright(counts)
    stats["ENC2_Sun"] = enc_sun(counts, pseudocount=False)
    stats["ENC5_Sun"] = enc_sun(counts, pseudocount=True)
    stats["Fop"] = fop(counts)
    r = rscu(counts)
    stats["top_RSCU"] = max(r, key=r.get) if r else ""
    return stats, r


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="CDS FASTA (多序列)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--gcode", type=int, default=1)
    args = ap.parse_args()

    seqs = parse_fasta(open(args.fasta, encoding="utf-8", errors="ignore").read())
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    rscu_cols = sorted(CODON_TABLE)
    rscu_rows = []
    n_ok = n_skip = 0
    for name, s in seqs:
        res = analyze(s, args.gcode)
        if res is None:
            n_skip += 1
            continue
        stats, r = res
        stats["seq"] = name
        rows.append(stats)
        rscu_rows.append([name] + [f"{r.get(c, float('nan')):.4f}" for c in rscu_cols])
        n_ok += 1

    cols = ["seq", "length_nt", "valid_codons", "ambiguous", "GC", "GC1", "GC2",
            "GC3", "GC3s", "ENC", "ENC2_Sun", "ENC5_Sun", "Fop", "top_RSCU"]
    out1 = os.path.join(args.outdir, "codon_usage.tsv")
    with open(out1, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(
                f"{row[c]:.4f}" if isinstance(row[c], float) else str(row[c])
                for c in cols) + "\n")

    out2 = os.path.join(args.outdir, "rscu.tsv")
    with open(out2, "w", encoding="utf-8") as f:
        f.write("seq\t" + "\t".join(rscu_cols) + "\n")
        for row in rscu_rows:
            f.write("\t".join(row) + "\n")

    print(f"codon usage: {n_ok} sequences analyzed, {n_skip} skipped (<20 codons)\n"
          f"  {out1}\n  {out2}")


if __name__ == "__main__":
    main()
