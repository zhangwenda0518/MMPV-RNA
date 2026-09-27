"""TIGER v2 位点速率评估 (Cummins et al. 2018, BMC Evolutionary Biology)。

按 PhyloSuite 内嵌 TIGER 2.0 源码 (src/tiger_{index,rate,output}.py) 语义复刻,
去掉 pickle 中间文件 / PyQt queue 依赖, 纯内存 + CLI。

用法:
    python -m utils.tiger_sites --aln input.aln --outdir out/ \
        [--bins 10] [--exclude 9,10] [--mask] [--unknowns ?,-]

产出:
    {prefix}.tiger.rates.csv  — Site,Rate,Bin 三列 (1-based 位点)
    {prefix}.tiger.fas        — 低速率 bin 剔除/掩蔽后的比对 (X 掩蔽或直接删除)
    {prefix}.tiger.hist.txt   — bin 直方图

原理: 位点速率 = 1 - 与其它非恒定位点 pattern 的平均相似度。
速率高 = 与其它位点相似 = 进化慢; 恒定位点硬编码 rate=1.0。
bin 编号与原版一致: bin 1 = 最高速率, bin N = 最低速率。
默认排除 bin 9,10 (最低速率 20% = 高变异/饱和位点), 与
PhyloSuite 主调用默认 [9,10] 一致。注意: bin 边界为双闭区间,
边界值归编号更小的 bin (原版 get_bin 同款行为)。
"""
import os
import argparse

# 模块级默认 (历史坑: 只靠 run() 内 global 赋值, 直接调用 site_pattern 会 NameError)
unknown_chars = set("?,-")


def parse_fasta(text):
    """-> [(name, seq)] 保序。"""
    data = []
    name = None
    seq_parts = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                data.append((name, "".join(seq_parts)))
            name = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    if name is not None:
        data.append((name, "".join(seq_parts)))
    return data


def check_aln(seqs):
    base = len(seqs[0][1])
    for name, s in seqs[1:]:
        if len(s) != base:
            raise ValueError(
                f"Sequences not aligned: {name} len={len(s)} != {base}")
    return base


def site_pattern(site):
    """列字符 -> pattern 字符串 (同字符位置聚组, '|' 分隔)。
    未知字符的位置从 pattern 中剔除 (与真实 TIGER CLI 一致;
    PhyloSuite 移植版未实现 -u, 此处按原版 TIGER 语义补全)。"""
    pat = {}
    order = []
    for i, base in enumerate(site):
        if base in unknown_chars:
            continue  # 剔除, 不参与打分
        try:
            pat[base].append(i)
        except KeyError:
            pat[base] = [i]
            order.append(base)
    return "|".join(",".join(str(x) for x in pat[o]) for o in order)


def set_pattern(p):
    return [frozenset(int(y) for y in x.split(",")) for x in p.split("|")]


def score(a, b):
    """b 的每组在 a 中找到超集的比例。"""
    found = 0
    for set_b in b:
        for set_a in a:
            if set_b.issubset(set_a):
                found += 1
                break
    return found / len(b)


def site_rate(a_str, uniq_pats, self_key):
    """TIGER site_rate: 1 - mean(1 - score(self, other)) over 非恒定位点。"""
    patA = set_pattern(a_str)
    pat_rates = []
    dividand = 0
    for p, info in uniq_pats.items():
        if "|" not in p:
            continue  # 恒定位点不参与
        pr = 1.0 - score(patA, set_pattern(p))
        reps = info["count"]
        if p == self_key:
            reps -= 1
        for _ in range(reps):
            pat_rates.append(pr)
            dividand += 1
    if dividand == 0:
        return 1.0
    return 1.0 - (sum(pat_rates) / dividand)


def bin_rates(rates_by_pattern, bin_no):
    """速率分 bin。注意编号方向与原版一致: bin 1 = 最高速率,
    bin_no = 最低速率 (原版 get_bin 的反向编号, 已数值验证)。"""
    all_r = [info["rate"] for info in rates_by_pattern.values()]
    upper, lower = max(all_r), min(all_r)
    step = (upper - lower) / bin_no
    if step == 0:
        raise ValueError("Data too homogeneous to bin")
    divs = [lower + step * i for i in range(bin_no + 1)]
    divs[-1] = upper
    bins = {}
    for k, info in rates_by_pattern.items():
        r = info["rate"]
        # 双闭区间, 先匹配优先 (与原版 get_bin 一致)
        # 编号: i=0 是最低速率区间 -> bin N; i=N-1 最高速率 -> bin 1
        # 即 bin 1 = 最高速率 (保守位点), bin N = 最低速率 (变异位点)
        assigned = bin_no
        for i in range(len(divs) - 1):
            if divs[i] <= r <= divs[i + 1]:
                assigned = (len(divs) - 1) - i
                break
        bins[k] = assigned
    return bins


def map_bins_to_positions(rates_by_pattern, bins):
    bin_map = {}
    for k, info in rates_by_pattern.items():
        for pos in info["sites"]:
            bin_map[pos] = bins[k]
    return bin_map


def run(aln_path, outdir, bins=10, exclude="9,10", mask=True, unknowns="?,-"):
    global unknown_chars
    unknown_chars = set(unknowns.split(",")) if unknowns else set()

    seqs = parse_fasta(open(aln_path, encoding="utf-8", errors="ignore").read())
    aln_len = check_aln(seqs)
    if len(seqs) < 3:
        raise ValueError("Need >=3 sequences for site-rate estimation")

    # 1. patterns
    pats = []
    for x in range(aln_len):
        col = [s[1][x] for s in seqs]
        pats.append(site_pattern(col))
    uniq = {}
    for xi, p in enumerate(pats):
        if p in uniq:
            uniq[p]["count"] += 1
            uniq[p]["sites"].append(xi)
        else:
            uniq[p] = {"count": 1, "sites": [xi]}

    # 2. rates
    for k in uniq:
        if "|" not in k:
            uniq[k]["rate"] = 1.0
        else:
            uniq[k]["rate"] = site_rate(k, uniq, k)

    # 3. bins
    bins_map = bin_rates(uniq, int(bins))
    bin_map = map_bins_to_positions(uniq, bins_map)

    # 4. outputs
    os.makedirs(outdir, exist_ok=True)
    prefix = os.path.splitext(os.path.basename(aln_path))[0]
    excl = set(int(b) for b in exclude.split(",")) if exclude else set()

    rates_csv = os.path.join(outdir, f"{prefix}.tiger.rates.csv")
    with open(rates_csv, "w") as f:
        f.write("Site,Rate,Bin\n")
        for pos in range(aln_len):
            f.write(f"{pos + 1},{uniq[pats[pos]]['rate']:.6f},{bin_map[pos]}\n")

    masked_fas = os.path.join(outdir, f"{prefix}.tiger.fas")
    with open(masked_fas, "w") as f:
        for name, s in seqs:
            chars = []
            for pos, base in enumerate(s):
                if bin_map[pos] in excl:
                    chars.append("X" if mask else "")
                else:
                    chars.append(base)
            f.write(f">{name}\n{''.join(chars)}\n")

    hist = os.path.join(outdir, f"{prefix}.tiger.hist.txt")
    with open(hist, "w") as f:
        counts = {}
        for pos in range(aln_len):
            b = bin_map[pos]
            counts[b] = counts.get(b, 0) + 1
        total = aln_len
        for b in sorted(counts):
            bar = "#" * max(1, int(counts[b] / total * 50))
            f.write(f"bin {b:2d}: {counts[b]:6d} ({counts[b]/total:5.1%}) {bar}\n")
        f.write(f"\nexcluded bins: {sorted(excl)} "
                f"({sum(c for b, c in counts.items() if b in excl)} sites)\n")

    n_excl = sum(1 for pos in range(aln_len) if bin_map[pos] in excl)
    return rates_csv, masked_fas, hist, aln_len, n_excl


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True, help="比对好的 FASTA")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bins", type=int, default=10)
    ap.add_argument("--exclude", default="9,10",
                    help="要剔除/掩蔽的 bin (默认 9,10 = 最低速率 20%, 与原版一致)")
    ap.add_argument("--mask", action="store_true", default=True,
                    help="低速率位点以 X 掩蔽 (默认); 否则直接删除")
    ap.add_argument("--no-mask", dest="mask", action="store_false")
    ap.add_argument("--unknowns", default="?,-", help="未知字符集")
    args = ap.parse_args()

    rates_csv, masked_fas, hist, n, n_excl = run(
        args.aln, args.outdir, args.bins, args.exclude, args.mask, args.unknowns)
    print(f"TIGER site rates: {n} sites, {n_excl} excluded "
          f"(bins {args.exclude})\n  {rates_csv}\n  {masked_fas}\n  {hist}")


if __name__ == "__main__":
    main()
