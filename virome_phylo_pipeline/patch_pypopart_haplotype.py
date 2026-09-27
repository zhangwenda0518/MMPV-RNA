#!/usr/bin/env python3
"""patch_pypopart_haplotype.py — 修复 pypopart 单倍型距离崩溃（根因）

问题:
  identify_haplotypes_from_alignment 用 remove_gaps() 后的序列作单倍型代表序列,
  去 gap 后序列不等长 → hamming_distance 长度检查必炸 (病毒比对有 indel)。

修复:
  单倍型分组的 key 仍用去 gap 序列 (逻辑不变), 但代表序列保留完整比对 (等长含 gap)。
  pypopart 的 hamming_distance(ignore_gaps=True) 与 π 计算 (a!='-' and b!='-')
  均为 pairwise gap 处理, 设计意图就是序列含 gap。

用法:  python patch_pypopart_haplotype.py [pypopart_src_dir]
"""
import os
import sys

# 默认按仓库定位（<repo>/biosoft/pypopart/src/pypopart），不写死本机绝对路径
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DEFAULT_SRC = os.path.join(_REPO, "biosoft", "pypopart", "src", "pypopart")


def main():
    # 历史坑: 默认 /tmp 路径会被服务器清理; 改为项目内 pypopart 真实路径
    src_dir = sys.argv[1] if len(sys.argv) > 1 else _DEFAULT_SRC
    target = os.path.join(src_dir, "core", "haplotype.py")
    if not os.path.exists(target):
        print(f"ERROR: {target} 不存在")
        sys.exit(1)

    src = open(target).read()
    old = """        if key not in haplotype_map:
            haplotype_map[key] = []
            sequence_map[key] = ungapped
"""
    new = """        if key not in haplotype_map:
            haplotype_map[key] = []
            # FIX(2026-08-28): 代表序列保留完整比对(等长含gap)而非 remove_gaps 后的序列。
            # 原实现去 gap 后不等长, hamming_distance 长度检查必炸; 而 hamming_distance
            # 与 π 计算均已有 pairwise gap 处理 (ignore_gaps / a!='-'), 设计意图即序列含 gap。
            sequence_map[key] = seq
"""
    if old not in src:
        print("OK: 已 patch 过或结构不同, 检查现状:")
        idx = src.find("sequence_map[key]")
        print(src[max(0, idx-200):idx+80])
        return
    src = src.replace(old, new)
    open(target, "w").write(src)
    print(f"PATCHED: {target}")

if __name__ == "__main__":
    main()
