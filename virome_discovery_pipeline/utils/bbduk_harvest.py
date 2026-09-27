#!/usr/bin/env python3
"""
bbduk_harvest.py — 用 BBDuk k-mer 匹配一趟捞取与 tip 序列重叠的 reads,
替代 Virseqimprover.py 中 salmon quant + filterbyname.sh 的两趟全量扫描。

分工约定 (与上游默认行为对齐):
  默认走 BBDuk; 需要覆盖度过滤语义时 (-readFrac > 0, 即 salmon --quasiCoverage)
  回退 salmon 路径 — BBDuk 没有 reads 级覆盖度权重概念, 语义不等价。

实测依据 (2026-09-04, Alternaria contig_116, 38.7M 对 reads, 10 线程):
  召回 k=25 100.000% (0 漏捞, 以 salmon+filterbyname 输出为真值), 配对完整,
  62-72s vs 234s (约 3.3x); reads 可直接读 .gz (BBTools 内部解压;
  不能经 bash 进程替换喂流 — FASTA 读取器要求真实文件)。
"""

import os
import shutil
import subprocess

# JVM 参数: 关断言 (BBTools wrapper 默认 -ea, 热循环纯开销) + 固定堆 (ref 极小,
# 默认自动抓 85% 物理内存毫无必要)
_BBDUK_JAVA_ARGS = ("-da", "-Xmx8g")


def _safe_path(p, what):
    """规范化并校验路径, 拒绝含父目录组件的显式穿越写法, 返回绝对路径。"""
    if not p:
        return ""
    norm = os.path.normpath(p)
    if os.pardir in norm.split(os.sep) or os.pardir in p.replace("\\", "/").split("/"):
        raise ValueError("Illegal path component in %s: %s" % (what, p))
    return os.path.abspath(norm)


def bbduk_harvest(read1, read2, work_dir, threads, k=25, bbduk_bin=None):
    """
    在 work_dir 下产出 tmp/mapped_reads_{1,2}.fastq 供 SPAdes 使用。

    read1/read2 : reads 文件 (.gz 或明文, FASTA/FASTQ 均可); read2 为空表示单端
    work_dir    : scaffold-truncated 目录 (需已含 scaffold-start-end.fasta)
    k           : k-mer 长度 (默认 25; 实测 25/31 召回均 100%, 20 噪音过多)

    返回 out1 (单端) 或 (out1, out2) (双端)。
    """
    r1 = _safe_path(read1, "read1")
    r2 = _safe_path(read2, "read2")
    work = _safe_path(work_dir, "work_dir")
    binp = _safe_path(bbduk_bin or shutil.which("bbduk.sh") or "", "bbduk_bin")

    if not os.path.isfile(r1):
        raise FileNotFoundError("read1 not found: %s" % r1)
    if r2 and not os.path.isfile(r2):
        raise FileNotFoundError("read2 not found: %s" % r2)
    if not binp or not os.path.isfile(binp):
        raise FileNotFoundError("bbduk.sh not found: %s" % (binp or "(PATH)"))
    ref = _safe_path(os.path.join(work, "scaffold-start-end.fasta"), "ref")
    fa = _safe_path(os.path.join(work, "scaffold.fasta"), "fasta")
    bed = _safe_path(os.path.join(work, "scaffold-start-end.bed"), "bed")
    for label, pth in (("fasta", fa), ("bed", bed)):
        if not os.path.isfile(pth):
            raise FileNotFoundError("%s not found: %s" % (label, pth))

    # tip 序列每次重建 (scaffold 每轮都在变; 等价于原 runAlignment 内的 bedtools getfasta)
    bedtools = shutil.which("bedtools")
    if not bedtools:
        raise FileNotFoundError("bedtools not found in PATH")
    sub = subprocess.run(
        [bedtools, "getfasta", "-fi", fa, "-bed", bed, "-fo", ref],
        shell=False, check=False, capture_output=True, text=True,
    )
    if sub.returncode != 0 or not os.path.isfile(ref):
        tail = "\n".join((sub.stderr or sub.stdout or "").splitlines()[-5:])
        raise RuntimeError("bedtools getfasta failed (rc=%d):\n%s" % (sub.returncode, tail))

    tmp = _safe_path(os.path.join(work, "tmp"), "tmp")
    if os.path.isdir(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp, exist_ok=True)

    out1 = _safe_path(os.path.join(tmp, "mapped_reads_1.fastq"), "out1")
    out2 = _safe_path(os.path.join(tmp, "mapped_reads_2.fastq"), "out2")

    # hdist: ref kmer 的容错位数 (1 = 允许 1 错配锚定, 捞回带测序错误的边缘 reads)。
    # 默认 0; 通过环境变量 VSI_BBDUK_HDIST 覆盖 (A/B 实验与调参入口)。
    try:
        hdist = int(os.environ.get("VSI_BBDUK_HDIST", "0"))
    except ValueError:
        hdist = 0

    # outm/outm2 = 匹配 ref 的 reads (等价 salmon+filterbyname 的 mapped 集合);
    # 注意不能用 out=/out2= (那是"未匹配"输出, 会把全量 reads 喂给 SPAdes)
    proc = subprocess.run(
        [
            binp,
            *_BBDUK_JAVA_ARGS,
            "in=%s" % r1,
            "outm=%s" % out1,
            "ref=%s" % ref,
            "k=%d" % int(k),
            "hdist=%d" % hdist,
            "rcomp=t",
            "threads=%d" % int(threads),
            "overwrite=t",
            *(("in2=%s" % r2, "outm2=%s" % out2) if r2 else ()),
        ],
        shell=False,
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr or proc.stdout or "").splitlines()[-8:])
        raise RuntimeError("bbduk failed (rc=%d):\n%s" % (proc.returncode, tail))
    for ln in (proc.stderr or "").splitlines():
        if any(key in ln for key in ("Input:", "Total Removed:", "Result:", "Time:")):
            print("[bbduk] %s" % ln.strip())

    return (out1, out2) if r2 else out1


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 5:
        print("Usage: python bbduk_harvest.py <read1> <read2|-> <work_dir> <threads> [k] [bbduk.sh]")
        sys.exit(1)
    r2arg = None if sys.argv[2] == "-" else sys.argv[2]
    karg = int(sys.argv[5]) if len(sys.argv) > 5 else 25
    binarg = sys.argv[6] if len(sys.argv) > 6 else None
    print(bbduk_harvest(sys.argv[1], r2arg, sys.argv[3], int(sys.argv[4]), karg, binarg))
