"""seqid_renamer.py — 序列 ID 重命名（自研，等效 VirPhyKit SeqIDRenamer）

输入 fasta + 映射（TSV 或 dict），输出重命名后的 fasta（Renamed_ 前缀）。
"""
from __future__ import annotations

import os
from typing import Dict, Optional


def rename_sequences(seq_file: str, rename_file: Optional[str] = None,
                     output_dir: str = '.',
                     rename_map: Optional[Dict[str, str]] = None,
                     output_prefix: str = 'Renamed_') -> str:
    """重命名 FASTA 序列头。

    参数
    ----
    seq_file : 输入 FASTA
    rename_file : 映射 TSV（旧名\\t新名，可选，与 rename_map 二选一）
    output_dir : 输出目录
    rename_map : 直接给映射 dict（可选）
    output_prefix : 输出文件名前缀

    返回输出文件路径；未命中映射的头保持原样。
    """
    if rename_map is None:
        rename_map = {}
        if rename_file:
            with open(rename_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) == 2:
                        rename_map[parts[0].strip()] = parts[1].strip()
    if not rename_map:
        raise ValueError('No rename mapping provided (rename_file or rename_map)')

    os.makedirs(output_dir, exist_ok=True)
    input_name = os.path.basename(seq_file)
    output_file = os.path.join(output_dir, f'{output_prefix}{input_name}')

    renamed = 0
    with open(seq_file, 'r', encoding='utf-8', errors='replace') as fin, \
            open(output_file, 'w', encoding='utf-8') as fout:
        for line in fin:
            if line.startswith('>'):
                seq_name = line[1:].strip()
                new_name = rename_map.get(seq_name)
                if new_name:
                    fout.write(f'>{new_name}\n')
                    renamed += 1
                else:
                    fout.write(line if line.endswith('\n') else line + '\n')
            else:
                fout.write(line if line.endswith('\n') else line + '\n')
    print(f'[seqid_renamer] {renamed}/{_count_seqs(seq_file)} heads renamed → {output_file}')
    return output_file


def _count_seqs(fasta_path: str) -> int:
    n = 0
    with open(fasta_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


def main():
    import argparse
    ap = argparse.ArgumentParser(description='序列 ID 重命名')
    ap.add_argument('fasta', help='输入 FASTA')
    ap.add_argument('--map', help='映射 TSV (旧名\\t新名)')
    ap.add_argument('--out', default='.', help='输出目录')
    args = ap.parse_args()
    if not args.map:
        raise SystemExit('需 --map 映射文件')
    rename_sequences(args.fasta, args.map, args.out)


if __name__ == '__main__':
    main()