#!/usr/bin/env python3
import sys

def fasta_to_phy_no_trunc(input_fasta, output_phy):
    records = {}
    current_id = None
    current_seq = []

    with open(input_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    records[current_id] = ''.join(current_seq)
                # 保留完整的 ID（仍然只取第一个空格前的部分作为名称，防止描述混入）
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.replace(' ', ''))
        if current_id is not None:
            records[current_id] = ''.join(current_seq)

    if not records:
        print("错误：文件中没有找到任何序列。")
        return

    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError("序列长度不一致，请确认输入文件已经是多重比对结果。")
    n_tax = len(records)
    n_char = lengths.pop()

    with open(output_phy, 'w') as out:
        out.write(f"{n_tax} {n_char}\n")
        for seq_id, seq in records.items():
            # 只用一个空格分隔完整 ID 和序列，不做截断、不定宽对齐
            out.write(f"{seq_id} {seq}\n")

    print(f"转换完成：{n_tax} 条序列，长度 {n_char}，已写入 {output_phy}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python fa2phy.py <输入.fasta> <输出.phy>")
        sys.exit(1)
    fasta_to_phy_no_trunc(sys.argv[1], sys.argv[2])
