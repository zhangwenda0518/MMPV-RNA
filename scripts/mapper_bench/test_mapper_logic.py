#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
直接验证 VSI 的 mapper 决策逻辑 (不跑 SPAdes 全流程)
从 Virseqimprover.py 提取比对段的真实行为
"""
import os, sys, subprocess, shutil, tempfile

os.chdir("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")

# 模拟 VSI 的全局变量与决策逻辑 (与源码一致)
threads = 10
mapperMode = "minibwa"
minibwaBin = os.path.expanduser("~/bin/minibwa")
bowtie2Bin = "bowtie2"
outputDir = "/tmp/vsi_mapper_logic_test"
read1 = "/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_1.fa.gz"
read2 = "/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_2.fa.gz"
bowtie2FastaFlag = " -f"

shutil.rmtree(outputDir, ignore_errors=True)
os.makedirs(outputDir)
shutil.copy("/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/scaffolds/ERR2040118_clean_NODE_150_length_3726_cov_29.718708.fasta",
            os.path.join(outputDir, "scaffold.fasta"))

# ---- 以下代码块从补丁后的 Virseqimprover.py 原样复制 ----
paired = (len(read2) > 0)
if paired:
    mbwaInput = read1 + " " + read2
    bt2Input = "-1 " + read1 + " -2 " + read2
else:
    mbwaInput = read1
    bt2Input = "-U " + read1

mbwaCmd = (
    "cd " + outputDir + "\n"
    + "set -o pipefail\n"
    + minibwaBin + " index -t " + str(threads) + " scaffold.fasta > /dev/null 2>&1\n"
    + minibwaBin + " map -u -t " + str(threads) + " scaffold.fasta " + mbwaInput
    + " | samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\n"
    + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\n"
)

bt2Cmd = (
    "cd " + outputDir + "\n"
    + "set -o pipefail\n"
    + bowtie2Bin + "-build --threads " + str(threads) + " scaffold.fasta bowtie2-index\n"
    + bowtie2Bin + " --local --no-unal --threads " + str(threads) + bowtie2FastaFlag
    + " -x bowtie2-index " + bt2Input
    + " | samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\n"
    + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\n"
)

def _run_mapper(scriptBody, tag):
    sf = open(outputDir + "/run.sh", 'w')
    sf.write('#!/bin/bash\n')
    sf.write(scriptBody)
    sf.close()
    try:
        out = subprocess.check_output("bash " + outputDir + "/run.sh",
                                      shell=True, stderr=subprocess.STDOUT)
        return True, out.decode('utf-8', 'replace')[-800:]
    except subprocess.CalledProcessError as e:
        tail = (e.output or b'').decode('utf-8', 'replace')[-800:]
        return False, tail

alignerUsed = mapperMode
if mapperMode == "minibwa":
    if not os.path.exists(minibwaBin):
        print("[mapper] minibwa not found -> bowtie2")
        alignerUsed = "bowtie2"
    else:
        print("[mapper] using minibwa")
        ok, tail = _run_mapper(mbwaCmd, "minibwa")
        if not ok:
            print("[mapper] minibwa FAILED:\n" + tail)
            alignerUsed = "bowtie2"
        else:
            print("[mapper] minibwa OK")

# ---- 结果校验 ----
print("")
print("=" * 50)
cov = os.path.join(outputDir, "samtools-coverage.txt")
bam = os.path.join(outputDir, "bowtie2-mapped.bam")
print("alignerUsed: " + alignerUsed)
print("samtools-coverage.txt exists: " + str(os.path.exists(cov)))
print("bowtie2-mapped.bam exists:    " + str(os.path.exists(bam)))
if os.path.exists(bam):
    r = subprocess.check_output("samtools view -c " + bam, shell=True).decode().strip()
    print("BAM records: " + r)
if os.path.exists(cov):
    with open(cov) as f:
        lines = f.readlines()
    dsum = sum(int(l.split('\t')[2]) for l in lines if len(l.split('\t')) >= 3)
    print("coverage lines: " + str(len(lines)) + ", depth sum: " + str(dsum))
print("=" * 50)
