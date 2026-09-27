#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare.py — 基准集去重 + 类别映射
==================================
`build_benchmark.sh` 用多条 NCBI 查询拼接出四类, 同一条记录可能被两条查询同时命中
(例如某条 EPRV 既在 "endogenous[Title] AND Caulimoviridae" 里, 也在 "endogenous
pararetrovirus" 里)。直接跑会让这些序列在分母里被数两次, 逐类比率就偏了。

本脚本按登录号去重, 输出:
  bench.fa          合并后的 fasta (id 前缀加类别, 便于下游按前缀取数)
  class_map.tsv     id <TAB> 类别 <TAB> 原始标题
"""
import argparse
import os
import re
import sys

CLASSES = ("A_extant_virus", "B1_eprv_clean", "B1b_eprv_decayed",
           "B2_eprv_activatable", "C_host_viral_domain")

# 类 -> 输入文件通配. 用通配而不是"合并后的单个文件": 合并是 build_benchmark.sh 的
# 最后一节, 前面任何一条 efetch 出问题都会让带 set -e 的脚本提前退出, 于是合并文件
# 不存在、prepare 静默跳过整类 —— 踩过: B1b 整类被跳过, 结果看起来"和上次一模一样",
# 差点当成"加了退化型 EPRV 后数字没变"。通配直接吃原始分片, 不依赖前置步骤。
FILES = {
    "A_extant_virus": "A_extant_virus.fna",
    "B1_eprv_clean": "B1_eprv_clean.fna",
    "B1b_eprv_decayed": "B1b_eprv_decayed_beet*.fna",
    "B2_eprv_activatable": "B2_eprv_activatable_all.fna",
    "C_host_viral_domain": "C_host_viral_domain.fna",
}


def read_fa(path):
    """-> [(id, desc, seq)]"""
    out, cur = [], None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    out.append(cur)
                head = line[1:].rstrip("\n")
                parts = head.split(None, 1)
                cur = [parts[0], parts[1] if len(parts) > 1 else "", []]
            elif cur:
                cur[2].append(line.strip())
    if cur:
        out.append(cur)
    return [(i, d, "".join(s)) for i, d, s in out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True, help="build_benchmark.sh 的输出目录")
    ap.add_argument("--out", required=True, help="输出目录")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    seen, rows, n_dup = {}, [], 0
    import glob as _glob
    for cls in CLASSES:
        pat = os.path.join(a.bench, FILES[cls])
        parts = sorted(_glob.glob(pat))
        if not parts:
            # 缺整类必须**响亮**报错, 不能静默跳过 —— 静默跳过会让结果看起来
            # "和上次一样", 极易被误读成"加了这一类后结论没变"。
            raise SystemExit("缺基准分片: %s (类 %s 无法评估)" % (pat, cls))
        for f in parts:
            for rid, desc, seq in read_fa(f):
                base = re.sub(r"\.\d+$", "", rid)   # NC_003378.1 -> NC_003378
                if base in seen:
                    n_dup += 1
                    continue
                seen[base] = cls
                new_id = "%s|%s" % (cls, base)
                rows.append((new_id, cls, rid, desc, seq))

    with open(os.path.join(a.out, "bench.fa"), "w", encoding="utf-8") as fo, \
            open(os.path.join(a.out, "class_map.tsv"), "w", encoding="utf-8") as fm:
        fm.write("id\tclass\tacc\tdesc\n")
        for new_id, cls, rid, desc, seq in rows:
            fo.write(">%s %s\n" % (new_id, desc[:120]))
            for i in range(0, len(seq), 60):
                fo.write(seq[i:i + 60] + "\n")
            fm.write("%s\t%s\t%s\t%s\n" % (new_id, cls, rid, desc[:120]))

    by = {}
    for _, cls, _, _, _ in rows:
        by[cls] = by.get(cls, 0) + 1
    print("去重后基准集: %d 条 (合并时丢掉重复 %d 条)" % (len(rows), n_dup))
    for c in CLASSES:
        if c in by:
            print("  %-24s %5d" % (c, by[c]))
    print("-> %s/bench.fa 与 class_map.tsv" % a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
