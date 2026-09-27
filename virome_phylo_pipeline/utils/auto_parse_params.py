#!/usr/bin/env python3
"""自动解析合并 log 关键参数 (供监控 v3 完成时调用)

用法: python -m utils.auto_parse_params <work> <merged>
"""
import sys, os, csv, json

# 脚本方式运行时 (python utils/auto_parse_params.py) sys.path 指向 utils/,
# utils.ess 不可见; 补根路径保证两种方式都能 import
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# 历史坑: 输入 merged.log 已被 logcombiner -burnin 烧过一次, 再硬烧 500 点是双重 burnin。
# 改为只跳过表头行 (i==0), 不再二次烧样本。
BURNIN = 0


def read_col(log_file, col):
    vals = []
    with open(log_file, encoding='utf-8', errors='replace') as f:
        i = -1
        for line in f:
            if line.startswith('#'):
                continue
            i += 1
            if i == 0 or i <= BURNIN:
                continue
            parts = line.strip().split('\t')
            try:
                vals.append(float(parts[col]))
            except (ValueError, IndexError):
                pass
    return vals


def stats(vals):
    n = len(vals)
    if n < 100:
        return None
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / (n - 1)
    if var == 0:
        ess = n
    else:
        # Tracer 同算法 (Geyer monotone sequence); 旧版 AR(1) 近似严重低估
        from utils.ess import calculate_ess
        ess = calculate_ess(vals)
    vs = sorted(vals)
    # 95% HPD (Tracer 口径, 最短区间) — 旧版用 2.5%/97.5% 分位数是等尾区间, 不是 HPD
    from utils.ess import calculate_95hpd
    hpd_lo, hpd_hi = calculate_95hpd(vals)
    return {'mean': mean, 'hpd_lower': hpd_lo, 'hpd_upper': hpd_hi,
            'ess': ess, 'n': n}


def main():
    if len(sys.argv) < 3:
        print("用法: python -m utils.auto_parse_params <work> <merged>")
        return 1
    work, merged = sys.argv[1], sys.argv[2]
    log_file = f'{merged}/merged.log'
    if not os.path.exists(log_file):
        print("merged.log 不存在")
        return 0

    hdr = None
    with open(log_file, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                continue
            hdr = line.strip().split('\t')
            break
    if not hdr:
        return 0

    result = {}
    for col, label in [('treeModel.rootHeight', 'tmrca'), ('ucld.mean', 'clock_rate'),
                       ('ucld.stdev', 'ucld_stdev')]:
        if col in hdr:
            s = stats(read_col(log_file, hdr.index(col)))
            if s:
                result[label] = s
    # skyline popSize 均值
    pop_vals = []
    for j in range(1, 8):
        col = f'skyline.popSize{j}'
        if col in hdr:
            s = stats(read_col(log_file, hdr.index(col)))
            if s:
                pop_vals.append(s['mean'])
    if pop_vals:
        result['popSize_mean'] = sum(pop_vals) / len(pop_vals)
        result['popSize_series'] = pop_vals

    out = f'{merged}/parameter_summary.json'
    with open(out, 'w', encoding='utf-8') as _jf:
        json.dump(result, _jf, indent=2)
    print(f"参数摘要: {out}")
    return 0


if __name__ == '__main__':
    sys.exit(main())