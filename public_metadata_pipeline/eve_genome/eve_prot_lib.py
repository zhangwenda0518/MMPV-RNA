#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_prot_lib.py — 把 01 阶段的核苷酸元件库翻成蛋白参考库（六框）。

为什么需要: `02`（发现管线的 EVE 判别）里 Δbit-score 的 `--eve-faa` 要的是
蛋白库（diamond makedb + blastx），而 `eve_element_db.py` 的
`eve_consensus_library.fa` 是核苷酸。这里按六框翻译、按终止符切 ORF，
只保留 >= min_aa 的片段，与 ICTV 六框库同一口径。

用法:
  python eve_prot_lib.py --in eve_consensus_library.fa \
      --out eve_consensus_library.sixframe.faa [--min-aa 30]
"""

import argparse
import sys

_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf)
    return seqs


def sixframe_faa(fa_in, out_faa, min_aa=30):
    seqs = read_fasta(fa_in)
    n = 0
    with open(out_faa, "w", encoding="utf-8", newline="\n") as out:
        for sid, s in seqs.items():
            s = s.upper()
            for strand, seq in (("+", s), ("-", s.translate(_COMP)[::-1])):
                for frame in range(3):
                    aa, start = [], frame
                    for i in range(frame, len(seq) - 2, 3):
                        a = _CODON.get(seq[i:i + 3], "X")
                        if a == "*":
                            if len(aa) >= min_aa:
                                n += 1
                                out.write(f">{sid}|{strand}{frame}|{start}-{i}\n"
                                          + "".join(aa) + "\n")
                            aa, start = [], i + 3
                        else:
                            aa.append(a)
                    if len(aa) >= min_aa:
                        n += 1
                        out.write(f">{sid}|{strand}{frame}|{start}-{len(seq)}\n"
                                  + "".join(aa) + "\n")
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="fa_in", required=True,
                    help="核苷酸元件库 (eve_consensus_library.fa)")
    ap.add_argument("--out", required=True, help="输出蛋白 FASTA")
    ap.add_argument("--min-aa", type=int, default=30,
                    help="保留的最短 ORF 氨基酸数")
    a = ap.parse_args(argv)
    n = sixframe_faa(a.fa_in, a.out, a.min_aa)
    print(f"[done] {n} ORFs -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
