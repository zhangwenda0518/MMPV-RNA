#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""验证 bowtie2 兜底路径真能产出 (补丁中 alignerUsed=='bowtie2' 分支)"""
import os, subprocess, shutil

outputDir = "/tmp/vsi_fallback_full"
shutil.rmtree(outputDir, ignore_errors=True)
os.makedirs(outputDir)
shutil.copy("/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/scaffolds/"
            "ERR2040118_clean_NODE_150_length_3726_cov_29.718708.fasta",
            outputDir + "/scaffold.fasta")

threads = 10
r1 = "/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_1.fa.gz"
r2 = "/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_2.fa.gz"
bowtie2FastaFlag = " -f"

bt2Cmd = (
    "cd " + outputDir + "\n"
    + "set -o pipefail\n"
    + "bowtie2-build --threads " + str(threads) + " scaffold.fasta bowtie2-index\n"
    + "bowtie2 --local --no-unal --threads " + str(threads) + bowtie2FastaFlag
    + " -x bowtie2-index -1 " + r1 + " -2 " + r2
    + " | samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\n"
    + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\n"
)

open(outputDir + "/run.sh", "w").write("#!/bin/bash\n" + bt2Cmd)
print("--- 执行 bowtie2 兜底路径 ---")
out = subprocess.run("bash " + outputDir + "/run.sh", shell=True,
                     capture_output=True, text=True)
print("returncode:", out.returncode)
if out.returncode != 0:
    print("stderr tail:", out.stderr[-800:])

bam = outputDir + "/bowtie2-mapped.bam"
cov = outputDir + "/samtools-coverage.txt"
print("bam exists:", os.path.exists(bam), "| cov exists:", os.path.exists(cov))
if os.path.exists(bam):
    print("records:", subprocess.check_output("samtools view -c " + bam,
                                              shell=True).decode().strip())
if os.path.exists(cov):
    dsum = sum(int(l.split('\t')[2]) for l in open(cov) if len(l.split('\t')) >= 3)
    print("depth sum:", dsum)
