"""geo_subsampler.py — 序列地理/区域亚采样（自研，等效 VirPhyKit GeoSubsampler）

三种模式:
  equal_sampling=True  : 按序列名解析出的地点等量抽样（每地点抽 min 数量）
  region=<关键词>      : 抽走指定区域内 N 条（输出被抽走的 + remaining）
  else                : 从全部序列随机抽 N 条

默认不覆盖原文件（VirPhyKit 原版 region 模式会覆盖; 这里 remaining 输出到 output_dir）。
"""
from __future__ import annotations

import os
import random
import re
from collections import OrderedDict

# 2026-09-15: 默认种子 (固定值)。旧版用全局未播种的 random.sample, 每次跑
# 抽到的序列都不同 -> 下游 Mantel / 树地理结果随运行漂移, 且无法复现。
SUBSAMPLE_DEFAULT_SEED = 20260915


def _read_fasta(fasta_file: str) -> "OrderedDict[str, str]":
    """读 FASTA → OrderedDict {header: seq}。"""
    dict_fas = OrderedDict()
    with open(fasta_file, 'r', encoding='utf-8', errors='replace') as f:
        line = f.readline()
        while line != '':
            while not line.startswith('>') and line != '':
                line = f.readline()
            if line == '':
                break
            fas_name = line.strip()
            fas_seq = ''
            line = f.readline()
            while not line.startswith('>') and line != '':
                fas_seq += re.sub(r'\s', '', line)
                line = f.readline()
            dict_fas[fas_name] = fas_seq
    return dict_fas


def _site_of(name: str) -> str:
    """从序列名提取地点（匹配下划线/结尾前面的字母数字段）。"""
    m = re.search(r'(?:^|>|_)([A-Za-z0-9]+)(?:_|$)', name)
    if not m:
        raise ValueError(f'Cannot parse site from sequence name: {name}')
    return m.group(1)


def _write_fasta(path: str, seqs):
    with open(path, 'w', encoding='utf-8') as f:
        # 与原版一致：join 不补尾部换行
        f.write('\n'.join(f'{name}\n{seq}' for name, seq in seqs))


def subsample_fasta(fasta_file: str, num_seqs: int, region: str = None,
                    output_dir: str = None, equal_sampling: bool = False,
                    output_file_name: str = 'extract.fas',
                    overwrite_original: bool = False,
                    seed: int = None) -> dict:
    """亚采样。

    返回 {'extract': 输出路径, 'remaining': remaining 路径(region 模式) 或 None,
          'n': 抽取数, 'seed': 实际使用的随机种子}

    可复现性 (2026-09-15 修复): 旧版三处都调全局未播种的 `random.sample`,
    同一份输入两次跑抽出的序列不同, 下游 Mantel / 树地理结果随之漂移。
    现改用独立的 `random.Random(seed)`; seed 缺省取固定常量
    SUBSAMPLE_DEFAULT_SEED, 实际使用的值随返回值一并返回便于追溯。
    """
    used_seed = SUBSAMPLE_DEFAULT_SEED if seed is None else int(seed)
    rng = random.Random(used_seed)
    dict_fas = _read_fasta(fasta_file)
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(fasta_file)) or '.'
    os.makedirs(output_dir, exist_ok=True)
    extract_file = os.path.join(output_dir, output_file_name)

    if equal_sampling:
        site_dict = OrderedDict()
        for name, seq in dict_fas.items():
            site = _site_of(name)
            site_dict.setdefault(site, []).append((name, seq))
        min_seq_num = min(len(seqs) for seqs in site_dict.values())
        if min_seq_num == 0:
            raise ValueError('Some sites have zero sequences.')
        sampled = []
        for seqs in site_dict.values():
            sampled.extend(rng.sample(seqs, min_seq_num))
        _write_fasta(extract_file, sampled)
        stats = {site: min_seq_num for site in site_dict}
        return {'extract': extract_file, 'remaining': None, 'n': len(sampled),
                'per_site': stats, 'seed': used_seed}

    elif region:
        region_seqs = [(n, s) for n, s in dict_fas.items() if region in n]
        total_region = len(region_seqs)
        if num_seqs > total_region:
            raise ValueError(f'Requested to remove {num_seqs} sequences, but only '
                             f'{total_region} sequences with region {region} are available.')
        if num_seqs == total_region:
            raise ValueError('Cannot remove all sequences, please backup the original file.')
        to_remove = rng.sample(region_seqs, num_seqs)
        remove_ids = {n for n, _ in to_remove}
        remaining = [(n, s) for n, s in dict_fas.items() if n not in remove_ids]
        _write_fasta(extract_file, to_remove)
        remaining_path = None
        if overwrite_original:
            _write_fasta(fasta_file, remaining)          # 原版行为
        else:
            remaining_path = os.path.join(output_dir, 'remaining.fas')
            _write_fasta(remaining_path, remaining)
        return {'extract': extract_file, 'remaining': remaining_path, 'n': num_seqs,
                'seed': used_seed}

    else:
        if num_seqs > len(dict_fas):
            raise ValueError(f'Requested {num_seqs} sequences, but only {len(dict_fas)} available.')
        sampled = rng.sample(list(dict_fas.items()), num_seqs)
        _write_fasta(extract_file, sampled)
        return {'extract': extract_file, 'remaining': None, 'n': num_seqs,
                'seed': used_seed}


def main():
    import argparse
    ap = argparse.ArgumentParser(description='序列地理/区域亚采样')
    ap.add_argument('fasta', help='输入 FASTA')
    ap.add_argument('-n', '--num', type=int, default=10, help='抽取数量')
    ap.add_argument('--region', default=None, help='区域关键词（抽走该区域 N 条）')
    ap.add_argument('--equal', action='store_true', help='按地点等量抽样')
    ap.add_argument('-o', '--out', default=None, help='输出目录')
    ap.add_argument('--seed', type=int, default=None,
                    help=f'随机种子 (缺省 {SUBSAMPLE_DEFAULT_SEED}, 固定值以保证可复现)')
    args = ap.parse_args()
    r = subsample_fasta(args.fasta, args.num, region=args.region,
                        output_dir=args.out, equal_sampling=args.equal,
                        seed=args.seed)
    print(r)


if __name__ == '__main__':
    main()