#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Patch Virseqimprover.py: minibwa 优先比对器，bowtie2 保留兜底。

改动点:
1. 新增全局变量 mapperMode (默认 minibwa) + minibwaBin
2. parseArguments 新增 -mapper 参数 (minibwa|bowtie2)
3. 比对命令生成改为 dual-path: minibwa 主路径 + 产物缺失自动回退 bowtie2
4. 产物名保持 bowtie2-mapped.bam / samtools-coverage.txt 不变 (下游 L912/L961 依赖)
"""
import re
import sys
import shutil

TARGET = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/Virseqimprover.py"

with open(TARGET, 'r', encoding='utf-8') as f:
    src = f.read()

orig = src

# ---------- 1) 全局变量 ----------
anchor_g = 'genus_avg_len = 0\n'
new_g = ('genus_avg_len = 0\n'
         '# 比对器: minibwa (默认, 实测比 bowtie2 --local 快 9.4x 且召回 +1.3%)\n'
         '# bowtie2 保留为兜底: minibwa 不可用或产物缺失时自动回退\n'
         'mapperMode = "minibwa"\n'
         'minibwaBin = os.path.expanduser("~/bin/minibwa")\n'
         'bowtie2Bin = "bowtie2"\n')
assert src.count(anchor_g) == 1, "全局变量锚点不唯一"
src = src.replace(anchor_g, new_g)

# ---------- 2) global 声明 ----------
anchor_gl = '    global alignerMode, bbdukBin, bbdukK\n'
new_gl = ('    global alignerMode, bbdukBin, bbdukK\n'
          '    global mapperMode, minibwaBin, bowtie2Bin\n')
assert src.count(anchor_gl) == 1, "global 锚点不唯一"
src = src.replace(anchor_gl, new_gl)

# ---------- 3) 参数解析 ----------
anchor_p = '''                    elif args[i] == "-checkv_db":
                        checkv_db = os.path.abspath(args[i + 1])'''
new_p = '''                    elif args[i] == "-mapper":
                        mapperMode = args[i + 1].lower()
                        if mapperMode not in ("minibwa", "bowtie2"):
                            print("Invalid -mapper: " + args[i + 1] + " (choose minibwa or bowtie2)")
                            return
                    elif args[i] == "-minibwa":
                        minibwaBin = os.path.abspath(args[i + 1])
                    elif args[i] == "-bowtie2":
                        bowtie2Bin = os.path.abspath(args[i + 1])
                    elif args[i] == "-checkv_db":
                        checkv_db = os.path.abspath(args[i + 1])'''
assert src.count(anchor_p) == 1, "参数解析锚点不唯一"
src = src.replace(anchor_p, new_p)

# ---------- 4) 帮助文本 ----------
anchor_h = '  -aligner <str>         Tip-reads harvester: bbduk (default, faster) or salmon;'
new_h = ('  -mapper <str>          Read mapper for extension: minibwa (default, ~9x faster)\n'
         '                         or bowtie2 (legacy). bowtie2 auto-fallback on failure.\n'
         '  -minibwa <path>        Path to minibwa binary (default: ~/bin/minibwa)\n'
         '  -bowtie2 <path>        Path to bowtie2 binary (default: bowtie2 in PATH)\n'
         '  -aligner <str>         Tip-reads harvester: bbduk (default, faster) or salmon;')
assert src.count(anchor_h) == 1, "帮助文本锚点不唯一"
src = src.replace(anchor_h, new_h)

# ---------- 5) 比对命令块 ----------
old_block = '''    if len(read2) == 0:
        cmd = str("cd " + outputDir + "\\n" \\
                  + "bowtie2-build --threads " + str(threads) + " scaffold.fasta bowtie2-index\\n" \\
                  + "bowtie2 --local --no-unal --threads " + str(threads) + bowtie2FastaFlag + " -x bowtie2-index " \\
                  + "-U " + read1 \\
                  + " | " \\
                  + "samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\\n" \\
                  + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\\n")
    else:
        cmd = str("cd " + outputDir + "\\n" \\
                  + "bowtie2-build --threads " + str(threads) + " scaffold.fasta bowtie2-index\\n" \\
                  + "bowtie2 --local --no-unal --threads " + str(threads) + bowtie2FastaFlag + " -x bowtie2-index " \\
                  + "-1 " + read1 \\
                  + " -2 " + read2 \\
                  + " | " \\
                  + "samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\\n" \\
                  + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = "bash " + outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)
'''

new_block = '''    # ---------- 比对命令: minibwa 优先, bowtie2 兜底 ----------
    # 产物名固定为 bowtie2-mapped.bam / samtools-coverage.txt (下游依赖, 不可改)
    paired = (len(read2) > 0)

    if paired:
        mbwaInput = read1 + " " + read2
        bt2Input = "-1 " + read1 + " -2 " + read2
    else:
        mbwaInput = read1
        bt2Input = "-U " + read1

    # minibwa 路径: index -> map -u -> sort -> depth
    mbwaCmd = (
        "cd " + outputDir + "\\n"
        + "set -o pipefail\\n"
        + minibwaBin + " index -t " + str(threads) + " scaffold.fasta > /dev/null 2>&1\\n"
        + minibwaBin + " map -u -t " + str(threads) + " scaffold.fasta " + mbwaInput
        + " | samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\\n"
        + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\\n"
    )

    # bowtie2 路径 (兜底, 也是 -mapper bowtie2 时的主路径)
    bt2Cmd = (
        "cd " + outputDir + "\\n"
        + "set -o pipefail\\n"
        + bowtie2Bin + "-build --threads " + str(threads) + " scaffold.fasta bowtie2-index\\n"
        + bowtie2Bin + " --local --no-unal --threads " + str(threads) + bowtie2FastaFlag
        + " -x bowtie2-index " + bt2Input
        + " | samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\\n"
        + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\\n"
    )

    def _run_mapper(scriptBody, tag):
        """写 run.sh 并执行, 返回 (ok, tail_of_output)"""
        sf = open(outputDir + "/run.sh", 'w')
        sf.write('#!/bin/bash\\n')
        sf.write(scriptBody)
        sf.close()
        try:
            out = subprocess.check_output("bash " + outputDir + "/run.sh",
                                          shell=True, stderr=subprocess.STDOUT)
            return True, out.decode('utf-8', 'replace')[-1500:]
        except subprocess.CalledProcessError as e:
            tail = (e.output or b'').decode('utf-8', 'replace')[-1500:]
            return False, tail

    alignerUsed = mapperMode
    runOk = False

    if mapperMode == "minibwa":
        if not os.path.exists(minibwaBin):
            print("[mapper] minibwa not found at " + minibwaBin + ", fallback to bowtie2")
            alignerUsed = "bowtie2"
        else:
            print('[mapper] using minibwa: ' + minibwaBin)
            runOk, tail = _run_mapper(mbwaCmd, "minibwa")
            if not runOk:
                print("[mapper] minibwa failed, fallback to bowtie2. tail:\\n" + tail)
                alignerUsed = "bowtie2"

    if alignerUsed == "bowtie2":
        # 清理 minibwa 可能留下的半成品, 避免下游读到坏产物
        for junk in ("bowtie2-mapped.bam", "samtools-coverage.txt"):
            p = os.path.join(outputDir, junk)
            if os.path.exists(p) and mapperMode == "minibwa":
                try: os.remove(p)
                except OSError: pass
        print('[mapper] using bowtie2: ' + bowtie2Bin)
        runOk, tail = _run_mapper(bt2Cmd, "bowtie2")
        if not runOk:
            print("[mapper] bowtie2 also failed. tail:\\n" + tail)
            raise RuntimeError("Both minibwa and bowtie2 failed for " + outputDir)

    # 产物存在性校验 (兜底机制的最后一道闸)
    if not os.path.exists(os.path.join(outputDir, "samtools-coverage.txt")):
        raise RuntimeError("Alignment product missing after " + alignerUsed + ": " + outputDir)

    print('[mapper] done with ' + alignerUsed)
'''

assert src.count(old_block) == 1, "比对命令块锚点不唯一 (可能已有改动)"
src = src.replace(old_block, new_block)

if src == orig:
    print("NO CHANGE APPLIED")
    sys.exit(1)

# 语法检查
import py_compile
with open(TARGET, 'w', encoding='utf-8') as f:
    f.write(src)
try:
    py_compile.compile(TARGET, doraise=True)
    print("PATCH OK, syntax OK")
except py_compile.PyCompileError as e:
    print("SYNTAX ERROR!", e)
    sys.exit(1)
