#!/usr/bin/env python3
"""
run_bbnorm.py — 独立 BBNorm 覆盖度归一化脚本
=============================================
自动检测 PE/SE, 逐对(或逐个)运行 bbnorm.sh 进行 k-mer 覆盖度归一化。

用法:
  python run_bbnorm.py -i 00b_HostDepletion/ -o 00c_BBnorm/ -t 16

  # 与 virome_pipeline 配合:
  python run_bbnorm.py -i out/00b_HostDepletion/ -o out/00c_BBnorm/ -t 16

参数:
  target=70    — 目标覆盖度 (高于此值会被降采样)
  mindepth=2   — 最低 k-mer 深度阈值 (低于此值被过滤)
  prefilter=t  — 内存优化模式
"""

import argparse, os, sys, subprocess, glob
from pathlib import Path

def find_reads(input_dir):
    """自动检测 PE/SE reads, 返回 (pe_pairs, se_files)"""
    d = Path(input_dir)
    r1_pats = ["*_R1*.fastq.gz", "*_R1*.fq.gz", "*_1.fastq.gz", "*_1.fq.gz", "*_1.fa.gz"]
    r2_pats = ["*_R2*.fastq.gz", "*_R2*.fq.gz", "*_2.fastq.gz", "*_2.fq.gz", "*_2.fa.gz"]
    all_r1 = set(); all_r2 = set()
    for pat in r1_pats:
        for f in d.glob(pat): all_r1.add(f)
    for pat in r2_pats:
        for f in d.glob(pat): all_r2.add(f)

    # PE 配对: 按 stem 匹配
    def _stem(p, tag):
        name = p.name
        for t in tag:
            idx = name.find(t)
            if idx >= 0: return name[:idx]
        return name.rsplit(".", 1)[0]

    r1_by_stem = {}
    for f in all_r1:
        s = _stem(f, ["_R1", "_1."])
        # 去 _R1/ _1 后缀再匹配 R2
        for t in ["_R1", "_1"]:
            if t in Path(f).name:
                s = Path(f).name.split(t)[0]
                break
        r1_by_stem.setdefault(s, []).append(f)

    r2_by_stem = {}
    for f in all_r2:
        for t in ["_R2", "_2"]:
            if t in Path(f).name:
                s = Path(f).name.split(t)[0]
                r2_by_stem.setdefault(s, []).append(f)
                break

    pe_pairs = []
    for stem in sorted(r1_by_stem):
        if stem in r2_by_stem:
            pe_pairs.append((r1_by_stem[stem][0], r2_by_stem[stem][0]))
        else:
            pe_pairs.append((r1_by_stem[stem][0], None))  # SE

    # 剩下的 R2 无匹配 R1 → SE fallback
    for stem in sorted(r2_by_stem):
        if stem not in r1_by_stem:
            pe_pairs.append((r2_by_stem[stem][0], None))

    se_files = []
    # 也检查非 _R1/_R2 文件
    for pat in ["*.fastq.gz", "*.fq.gz", "*.fa.gz"]:
        for f in d.glob(pat):
            fn = f.name
            if "_R1" not in fn and "_R2" not in fn and "_1." not in fn and "_2." not in fn:
                if f not in all_r1 and f not in all_r2:
                    se_files.append(f)
    se_files = sorted(set(se_files))

    return pe_pairs, se_files

def main():
    p = argparse.ArgumentParser(description="BBNorm 覆盖度归一化 (自动 PE/SE 识别)")
    p.add_argument("-i", "--input", required=True, help="输入 reads 目录")
    p.add_argument("-o", "--output", required=True, help="输出目录")
    p.add_argument("-t", "--threads", type=int, default=16, help="线程 (默认 16)")
    p.add_argument("--target", type=int, default=70, help="目标覆盖度 (默认 70)")
    p.add_argument("--mindepth", type=int, default=2, help="最低 k-mer 深度 (默认 2)")
    p.add_argument("--keep-se-as-pe", action="store_true", help="SE 文件也按 PE 强制配对 (默认: 单独处理)")
    args = p.parse_args()

    inp = Path(args.input)
    if not inp.is_dir():
        sys.exit(f"输入目录不存在: {inp}")

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    pe_pairs, se_files = find_reads(inp)
    n_pe = len([p for p in pe_pairs if p[1] is not None])
    n_se = len([p for p in pe_pairs if p[1] is None]) + len(se_files)

    print(f"检测: {n_pe} PE 对, {n_se} SE 文件")
    print(f"参数: target={args.target}, mindepth={args.mindepth}, prefilter=t, threads={args.threads}")
    print(f"输入: {inp}")
    print(f"输出: {out}")

    # PE pairs
    done, total = 0, 0
    for r1, r2 in pe_pairs:
        stem = Path(r1).name
        for t in ["_R1", "_1"]:
            if t in stem:
                stem = stem.split(t)[0]
                break
        nr1 = out / f"{stem}_norm_R1.fq.gz"
        nr2 = out / f"{stem}_norm_R2.fq.gz"
        if nr1.is_file() and nr1.stat().st_size > 0:
            done += 1; total += 1
            continue
        total += 1
        status = "PE" if r2 else "SE"
        print(f"  [{total}/{n_pe+n_se}] {status} {stem} ...", end=" ", flush=True)
        cmd = [f"bbnorm.sh", f"in1={r1}", f"out1={nr1}",
               f"target={args.target}", f"mindepth={args.mindepth}",
               f"prefilter=t", f"threads={args.threads}"]
        if r2:
            cmd.insert(2, f"in2={r2}")
            cmd.insert(4, f"out2={nr2}")
        try:
            subprocess.run(" ".join(cmd), shell=True, capture_output=True, check=False)
            if nr1.is_file() and nr1.stat().st_size > 0:
                print("OK")
            else:
                print("EMPTY")
        except: print("FAILED")

    # SE files
    for se in se_files:
        stem = se.name.rsplit(".", 1)[0]
        nr1 = out / f"{stem}_norm_SE.fq.gz"
        if nr1.is_file() and nr1.stat().st_size > 0:
            done += 1; total += 1
            continue
        total += 1
        print(f"  [{total}/{n_pe+n_se}] SE {stem} ...", end=" ", flush=True)
        cmd = [f"bbnorm.sh", f"in={se}", f"out={nr1}",
               f"target={args.target}", f"mindepth={args.mindepth}",
               f"prefilter=t", f"threads={args.threads}"]
        try:
            subprocess.run(" ".join(cmd), shell=True, capture_output=True, check=False)
            print("OK" if nr1.is_file() else "EMPTY")
        except: print("FAILED")

    output_files = list(out.glob("*_norm_*.fq.gz"))
    print(f"\n完成: {len(output_files)} 个归一化文件 → {out}")

if __name__ == "__main__":
    main()
