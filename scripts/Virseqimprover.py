import subprocess
import os
import sys
import math
import datetime

# Required inputs
read1 = ""
read2 = "" # Automatically transfer to single-end reads if not specified.
outputDir = ""
scaffold = ""

spadesKmerlen = "default"
minOverlapCircular = 5000
minIdentityCircular = 95
salmonReadFraction = 0
minSuspiciousLen = 1000
threads = 16
salmonBin = "salmon"
checkv_db = ""
run_counter = 1
checkv_triggered = False


def printHelp():
    helpText = """
Virseqimprover - Iteratively extend viral genome scaffolds using read mapping and assembly.

Usage:
  python Virseqimprover.py -1 <read1.fastq> -scaffold <scaffold.fasta> -o <output_dir> [options]

Required arguments:
  -1 <file>              Input read file 1 (FASTQ or FASTA, .gz supported)
  -scaffold <file>       Input scaffold file (FASTA format)
  -o <dir>               Output directory

Optional arguments:
  -2 <file>              Input read file 2 for paired-end reads
  -t, --threads <int>    Number of threads (default: 16)
  -salmon <path>         Path to salmon binary (default: salmon)
  -checkv_db <dir>       Path to CheckV database. If provided, enables auto-stop when completeness > 90%
  -spadeskmer <str>      SPAdes k-mer length (default: auto)
  -minOverlapCircular <int>   Minimum overlap for circular detection (default: 5000)
  -minIdentityCircular <int>  Minimum identity (%) for circular detection (default: 95)
  -readFrac <float>      Salmon read fraction for quasiCoverage (default: 0)
  -minSuspiciousLen <int>     Minimum length for suspicious region detection (default: 1000)
  -h, --help             Show this help message and exit

Example:
  python Virseqimprover.py -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \\
      -scaffold cluster_20.ref.fasta -o output_dir -t 32 -checkv_db ~/database/virus-db/checkv-db-v1.6
"""
    print(helpText)
    sys.exit(0)


def parseArguments(args):
    global outputDir, read1, read2, scaffold
    global spadesKmerlen, minOverlapCircular, minIdentityCircular, salmonReadFraction, minSuspiciousLen
    global threads, salmonBin, checkv_db

    if len(args) == 0:
        printHelp()
    else:
        skip_next = False
        for i in range(len(args)):
            if skip_next:
                skip_next = False
                continue
            if args[i] == "-h" or args[i] == "--help":
                printHelp()
            elif args[i][0] == "-":
                if (i+1) >= len(args):
                    print("Missing argument after " + args[i] + ".")
                    return
                else:
                    if args[i] == "-o":
                        outputDir = os.path.abspath(args[i + 1])
                    elif args[i] == "-1":
                        read1 = os.path.abspath(args[i + 1])
                    elif args[i] == "-2":
                        read2 = os.path.abspath(args[i + 1])
                    elif args[i] == "-scaffold":
                        scaffold = os.path.abspath(args[i + 1])
                    elif args[i] == "-spadeskmer":
                        spadesKmerlen = args[i + 1]
                    elif args[i] == "-minOverlapCircular":
                        minOverlapCircular = int(args[i + 1])
                    elif args[i] == "-minIdentityCircular":
                        minIdentityCircular = int(args[i + 1])
                    elif args[i] == "-readFrac":
                        salmonReadFraction = float(args[i + 1])
                    elif args[i] == "-minSuspiciousLen":
                        minSuspiciousLen = int(args[i + 1])
                    elif args[i] == "-t" or args[i] == "--threads":
                        threads = int(args[i + 1])
                    elif args[i] == "-salmon":
                        salmonBin = args[i + 1]
                    elif args[i] == "-checkv_db":
                        checkv_db = os.path.abspath(args[i + 1])
                    else:
                        print("Invalid argument: " + args[i])
                        return
                skip_next = True


def reverse_complement(seq):
    """
    生成 DNA 序列的反向互补链
    """
    mapping = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return seq.translate(mapping)[::-1]


def read_fasta(file_path):
    """
    极简 FASTA 读取器，将序列合并为一条大写字符串
    """
    seq = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    seq.append(line.strip())
    except Exception:
        pass
    return "".join(seq).upper()


def verify_tandem_repeat(seq, delta):
    """
    核心校验引擎：
    检查序列两端新增的 delta 长度，是否与相邻的 delta 长度序列内容一致（包含反向互补校验）
    """
    if delta <= 0 or len(seq) < 2 * delta:
        return False

    def is_similar(s1, s2):
        if len(s1) != len(s2) or len(s1) == 0:
            return False
        mismatches = sum(1 for a, b in zip(s1, s2) if a != b)
        max_mismatches = max(1, int(len(s1) * 0.05))
        return mismatches <= max_mismatches

    # 提取 3' 端的两个相邻 block
    end3_new = seq[-delta:]
    end3_adj = seq[-2*delta : -delta]

    # 提取 5' 端的两个相邻 block
    end5_new = seq[:delta]
    end5_adj = seq[delta : 2*delta]

    # 1. 检查 3' 端是否套娃（正向匹配 or 反向互补匹配）
    if is_similar(end3_new, end3_adj) or is_similar(end3_new, reverse_complement(end3_adj)):
        return True

    # 2. 检查 5' 端是否套娃
    if is_similar(end5_new, end5_adj) or is_similar(end5_new, reverse_complement(end5_adj)):
        return True

    return False


def getReadLen():
    print("getReadLen:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    cmd = ""
    if read1.endswith(".gz"):
        cmd = "gzip -dc " + read1 + \
        " | awk 'NR%4 == 2 {lenSum+=length($0); readCount++;} END {print lenSum/readCount}'"
    else:
        cmd = "awk 'NR%4 == 2 {lenSum+=length($0); readCount++;} END {print lenSum/readCount}' " + read1

    os.makedirs(outputDir, exist_ok=True)
    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = "bash " + outputDir + "/run.sh"
    reader = subprocess.check_output(cmd, shell=True)
    reader = reader.decode()
    str_val = reader.replace("\n", "")
    global avgReadLen
    avgReadLen = float(str_val)

    if avgReadLen == 0:
        print("Could not extract average read length from read file.")
        return
    print("Estimated average read length: " + str(avgReadLen))

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def createBed():
    print("createBed:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    scaffoldFile = outputDir + "/scaffold-truncated/tmp/spades-res/scaffold.fasta"
    if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        cmd = "cd " + outputDir + "/scaffold-truncated\n" + "rm -f scaffold.fasta*\n" + "cd tmp/spades-res\n" + "mv scaffold.fasta " + outputDir + "/scaffold-truncated\n" + "mv scaffold.fasta.fai " + outputDir + "/scaffold-truncated\n"

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    br = open(outputDir + "/scaffold-truncated/scaffold.fasta.fai")
    bw = open(outputDir + "/scaffold-truncated/scaffold-start-end.bed",'w')
    bwOutLog = open(outputDir + "/output-log.txt",'a')
    str_val = br.readline()
    results = []
    results = str_val.split("\t")
    scaffoldLength = int(results[1].strip())

    if scaffoldLength >= avgReadLen * 2:
        scaffoldId = results[0].strip()
        bw.write(str(scaffoldId) + "\t" + str(0) + "\t" + str(math.ceil(avgReadLen * 1.5)) + "\n")
        bw.write(str(scaffoldId) + "\t" + str(math.ceil(scaffoldLength - avgReadLen * 1.5)) + "\t" + str(scaffoldLength) + "\n")
        bwOutLog.write("Trying to grow scaffold " + scaffoldId + " with length " + str(scaffoldLength) + "\n")
        if scaffoldLength > 300000:
            bwOutLog.write("Length of " + scaffoldId + \
                           " is already greater than 300kbp, so stop extending this one.\n")

    br.close()
    bw.close()
    bwOutLog.close()

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def runAlignment():
    print("runAlignment:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    cmd = ""
    if len(read2) == 0:
        cmd = str("cd " + outputDir + "/scaffold-truncated\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta\n" \
                  + "rm -rf salmon-index\n" \
                  + "rm -rf salmon-res\n" \
                  + "rm -f salmon-mapped.sam\n" \
                  + salmonBin + " index -t scaffold-start-end.fasta -i salmon-index\n" \
                  + salmonBin + " quant -i salmon-index -l A " \
                  + "-r " + read1 + " -o salmon-res --writeMappings -p " + str(threads) + " --quasiCoverage " \
                  + str(salmonReadFraction) \
                  + " | samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam\n")
    else:
        cmd = str("cd " + outputDir + "/scaffold-truncated\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta\n" \
                  + "rm -rf salmon-index\n" \
                  + "rm -rf salmon-res\n" \
                  + "rm -f salmon-mapped.sam\n" \
                  + salmonBin + " index -t scaffold-start-end.fasta -i salmon-index\n" \
                  + salmonBin + " quant -i salmon-index -l A " \
                  + "-1 " + read1 + " -2 " + read2 + " -o salmon-res --writeMappings -p " + str(threads) + " --quasiCoverage " \
                  + str(salmonReadFraction) \
                  + "| samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\nset -e\nset -o pipefail\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = "bash " + outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def getMappedReads():
    print("getMappedReads:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    cmd = ""
    if len(read2) == 0:
        cmd = str("cd " + outputDir + "/scaffold-truncated\n" \
                  + "rm -rf tmp\n" \
                  + "mkdir tmp\n" \
                  + "bash filterbyname.sh in=" + read1 \
                  + " out=tmp/mapped_reads_1.fastq names=" \
                  + "salmon-mapped.sam include=t\n")
    else:
        cmd = str("cd " + outputDir + "/scaffold-truncated\n" \
                  + "rm -rf tmp\n" \
                  + "mkdir tmp\n" \
                  + "bash filterbyname.sh in=" + read1 + " in2=" + read2 \
                  + " out=tmp/mapped_reads_1.fastq out2=tmp/mapped_reads_2.fastq names=" \
                  + "salmon-mapped.sam include=t\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = "bash " + outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def runSpades():
    print("runSpades:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    cmd = ""
    if len(read2) == 0:
        if spadesKmerlen == "default":
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp\n" \
                      + "spades.py --rnaviral --threads " + str(threads) + " -o " \
                      + "spades-res -s mapped_reads_1.fastq " \
                      + "--trusted-contigs ../scaffold.fasta" \
                      + " --only-assembler\n" \
                      + "cd spades-res\n" \
                      + "samtools faidx scaffolds.fasta\n")
        else:
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp\n" \
                      + "spades.py --rnaviral --threads " + str(threads) + " -o " \
                      + "spades-res -s mapped_reads_1.fastq " \
                      + "--trusted-contigs ../scaffold.fasta -k " + str(spadesKmerlen) \
                      + " --only-assembler\n" \
                      + "cd spades-res\n" \
                      + "samtools faidx scaffolds.fasta\n")
    else:
        if spadesKmerlen == "default":
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp\n" \
                      + "spades.py --rnaviral --threads " + str(threads) + " -o " \
                      + "spades-res -1 mapped_reads_1.fastq -2 mapped_reads_2.fastq " \
                      + "--trusted-contigs ../scaffold.fasta" \
                      + " --only-assembler\n" \
                      + "cd spades-res\n" \
                      + "samtools faidx scaffolds.fasta\n")
        else:
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp\n" \
                      + "spades.py --rnaviral --threads " + str(threads) + " -o " \
                      + "spades-res -1 mapped_reads_1.fastq -2 mapped_reads_2.fastq " \
                      + "--trusted-contigs ../scaffold.fasta -k " + str(spadesKmerlen) \
                      + " --only-assembler\n" \
                      + "cd spades-res\n" \
                      + "samtools faidx scaffolds.fasta\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = "bash " + outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def getScaffoldFromScaffolds():
    print("getScaffoldFromScaffolds:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    maxScaffoldLength = 0
    maxScaffoldId = ""

    scaffoldFile = outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta"
    if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        br = open(outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta.fai")
        str_val = br.readline()
        str_val.strip()
        results = []
        results = str_val.split("\t")
        maxScaffoldId = results[0].strip()
        maxScaffoldLength = int(results[1].strip())
        br.close()
        if (maxScaffoldId != "") and (maxScaffoldLength != 0):
            cmd = str("cd " + outputDir + "/scaffold-truncated/tmp/spades-res\n" \
                      + "bash filterbyname.sh " \
                      + "in=scaffolds.fasta " \
                      + "out=scaffold.fasta names=" + maxScaffoldId \
                      + " include=t\n" \
                      + "samtools faidx scaffold.fasta\n")

            shellFileWriter = open(outputDir + "/run.sh",'w')
            shellFileWriter.write('#'+"!/bin/bash\n")
            shellFileWriter.write(cmd)
            shellFileWriter.close()

            cmd = "bash " + outputDir + "/run.sh"
            subprocess.check_output(cmd, shell=True)
    else:
        br = open(outputDir + "/scaffold-truncated/scaffold.fasta.fai")
        str_val = br.readline()
        str_val.strip()
        results = []
        results = str_val.split("\t")
        maxScaffoldLength = int(results[1].strip())
        br.close()

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))
    return maxScaffoldLength


def getScaffoldLenFromTruncatedExtend():
    print("getScaffoldLenFromTruncatedExtend:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    scaffoldLength = 0
    scaffoldFile = str(outputDir + "/scaffold-truncated/scaffold.fasta")
    if (os.path.exists(scaffoldFile)) == True and (os.path.isdir(scaffoldFile) == False):
        br = open(outputDir + "/scaffold-truncated/scaffold.fasta.fai")
        str_val = br.readline()
        str_val.strip()
        results = []
        results = str_val.split("\t")
        scaffoldLength = int(results[1].strip())
        br.close()
    else:
        br = open(outputDir + "/scaffold.fasta.fai")
        str_val = br.readline()
        str_val.strip()
        results = []
        results = str_val.split("\t")
        scaffoldLength = int(results[1].strip())
        br.close()

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))
    return scaffoldLength


def updateCurrentScaffold():
    print("updateCurrentScaffold:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    scaffoldFile = str(outputDir + "/scaffold-truncated/scaffold.fasta")
    if (os.path.exists(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        cmd = str("cd " + outputDir + "\n" \
                  + "rm -f scaffold.fasta*\n" \
                  + "cd scaffold-truncated\n" \
                  + "cp scaffold.fasta " + outputDir + "\n" \
                  + "cp scaffold.fasta.fai " + outputDir + "\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def createBedForTruncatedScaffold(truncatedLen):
    print("createBedForTruncatedScaffold:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    br = open(outputDir + "/scaffold.fasta.fai")
    bw = open(outputDir + "/scaffold-truncated.bed",'w')
    str_val = br.readline()
    results = str_val.split("\t")
    scaffoldLength = int(results[1].strip())
    if scaffoldLength >= truncatedLen:
        scaffoldId = results[0].strip()
        bw.write(scaffoldId + "\t" + str(truncatedLen) + "\t" + str(scaffoldLength - truncatedLen) + "\n")
    br.close()
    bw.close()

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def growScaffoldWithAssembly():
    global checkv_triggered
    print("growScaffoldWithAssembly:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    iteration = 1
    extendContig = True
    prevLength = 0
    delta_history = []
    while(extendContig):
        currentLength = getScaffoldFromScaffolds()
        if currentLength > 300000:
            extendContig = False
        elif currentLength > prevLength:
            if prevLength > 0:
                delta = currentLength - prevLength
                delta_history.append(delta)

                # 连续3次延伸相同长度 → 触发序列内容校验
                if len(delta_history) >= 3 and delta_history[-1] == delta_history[-2] == delta_history[-3]:
                    fasta_path = outputDir + "/scaffold-truncated/tmp/spades-res/scaffold.fasta"
                    if not os.path.exists(fasta_path):
                        fasta_path = outputDir + "/scaffold-truncated/scaffold.fasta"

                    curr_seq = read_fasta(fasta_path)

                    if verify_tandem_repeat(curr_seq, delta):
                        bwLog = open(outputDir + "/output-log.txt", 'a')
                        bwLog.write("\n[BLOCKED] Rigorous tandem repeat loop detected (+" + str(delta) + "bp sequence confirmed!). Stopping assembly.\n")
                        bwLog.close()
                        extendContig = False
                        break

            prevLength = currentLength
            createBed()
            runAlignment()
            getMappedReads()
            runSpades()

            scaffoldFile = outputDir + "/scaffold-truncated/tmp/spades-res/scaffolds.fasta"
            if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
                iteration += 1
                extendContig = True
            else:
                extendContig = False

            # 每5轮在 growing scaffold 上运行 CheckV，防止过度延伸
            if checkv_db != "" and iteration % 5 == 0:
                grow_fasta = outputDir + "/scaffold-truncated/tmp/spades-res/scaffold.fasta"
                if not os.path.exists(grow_fasta):
                    grow_fasta = outputDir + "/scaffold-truncated/scaffold.fasta"
                if os.path.exists(grow_fasta):
                    is_complete, exp_len = runCheckV(grow_fasta)
                    if is_complete:
                        checkv_triggered = True
                        extendContig = False
                        break
        else:
            extendContig = False

        if extendContig and iteration > 2000:
            extendContig = False

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def getTruncatedScaffoldAndExtend(currentLength):
    print("getTruncatedScaffoldAndExtend:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    truncation_lengths = [300, 500, 700, 1000, 1300, 1500, 1700, 2000]

    for trunc_len in truncation_lengths:
        createBedForTruncatedScaffold(trunc_len)
        cmd = str("cd " + outputDir + "\n" \
                  + "rm -rf scaffold-truncated\n" \
                  + "mkdir scaffold-truncated\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta\n" \
                  + "cd scaffold-truncated\n" \
                  + "samtools faidx scaffold.fasta\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

        growScaffoldWithAssembly()
        lengthFromGrowingTruncatedScaffold = getScaffoldLenFromTruncatedExtend()
        if lengthFromGrowingTruncatedScaffold > currentLength:
            return

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def checkCircularity():
    print("checkCircularity:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    file_list = os.listdir(outputDir)
    has_doc = False
    for file in file_list:
        path = os.path.join(outputDir, file)
        if "blastn-" in path:
            has_doc = True
            break
    if has_doc:
        cmd = "cd " + outputDir + "\n" + "rm -rf blastn-*\n"

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    # 加长探针：至少 300bp，验证前后序列是否一致
    br = open(outputDir + "/scaffold.fasta.fai")
    bw1 = open(outputDir + "/blastn-subject-1stround.bed",'w')
    bw2 = open(outputDir + "/blastn-query-1stround.bed",'w')
    str_val = br.readline()
    results = str_val.split("\t")
    scaffoldLength = int(results[1].strip())

    check_len = max(int(avgReadLen), 300)
    if check_len > scaffoldLength // 3:
        check_len = scaffoldLength // 3

    if scaffoldLength > (check_len * 2):
        scaffoldId = results[0].strip()
        bw1.write(str(scaffoldId) + "\t" + str(0) + "\t" + str(scaffoldLength - check_len) + "\n")
        bw2.write(str(scaffoldId) + "\t" + str(scaffoldLength - check_len) + "\t" + str(scaffoldLength) + "\n")
    br.close()
    bw1.close()
    bw2.close()

    cmd = str("cd " + outputDir + "\n" \
              + "bedtools getfasta -fi scaffold.fasta -bed blastn-subject-1stround.bed -fo blastn-subject-1stround.fasta\n" \
              + "bedtools getfasta -fi scaffold.fasta -bed blastn-query-1stround.bed -fo blastn-query-1stround.fasta\n" \
              + "samtools faidx blastn-subject-1stround.fasta\n" \
              + "samtools faidx blastn-query-1stround.fasta\n" \
              + "makeblastdb -in blastn-subject-1stround.fasta -dbtype nucl\n" \
              + "blastn -query blastn-query-1stround.fasta -db blastn-subject-1stround.fasta -num_threads " + str(threads) + " -outfmt '7' -out blastn-res-1stround.txt\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = "bash " + outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    # parse blastn result
    isCircular = False
    subjectStart = 0
    queryStart = 0
    br = open(outputDir + "/blastn-res-1stround.txt")
    bwCircularOutputLog = open(outputDir + "/circularity-output-log.txt",'a')
    bwOutputLog = open(outputDir + "/output-log.txt",'a')
    str_val = br.readline()
    results = []
    while(str_val != "" and isCircular == False):
        if str_val[0] != "#":
            str_val = str_val.strip()
            results = str_val.split("\t")
            percentIden = float(results[2].strip())
            alignmentLen = int(results[3].strip())
            subjectStart = int(results[8].strip())
            subjectEnd = int(results[9].strip())
            queryStart = int(results[6].strip())

            # 要求比对长度达到探针长度的 90%
            minAlignmentLength = int(check_len * 0.90)

            if (percentIden >= minIdentityCircular) and (alignmentLen >= minAlignmentLength):
                bwCircularOutputLog.write("Scaffold seems TRULY circular (Flanking sequences match!). Scaffold position " + str(results[6]) \
                                          + " to " + str(results[7]) + " mapped to position " + str(results[8]) + " to " + str(results[9]) \
                                          + " with " + str(percentIden) + "% identity and " + str(alignmentLen) + " alignment length.\n")
                bwOutputLog.write("Scaffold seems TRULY circular (Flanking sequences match!). Scaffold position " + str(results[6]) \
                                  + " to " + str(results[7]) + " mapped to position " + str(results[8]) + " to " + str(results[9]) \
                                  + " with " + str(percentIden) + "% identity and " + str(alignmentLen) + " alignment length.\n")
                isCircular = True
            else:
                if alignmentLen >= int(avgReadLen * 0.95):
                    bwCircularOutputLog.write("Ignored internal repeat (Flanking sequences mismatch). Alignment length " + str(alignmentLen) + " is less than required " + str(minAlignmentLength) + ".\n")
                    bwOutputLog.write("Ignored internal repeat (Flanking sequences mismatch). Alignment length " + str(alignmentLen) + " is less than required " + str(minAlignmentLength) + ".\n")

        str_val = br.readline()
    br.close()
    bwOutputLog.close()
    bwCircularOutputLog.close()

    # blastn subject query 2nd round
    if isCircular:
        br = open(outputDir + "/scaffold.fasta.fai")
        bw1 = open(outputDir + "/blastn-subject-2ndround.bed",'w')
        bw2 = open(outputDir + "/blastn-query-2ndround.bed",'w')
        str_val = br.readline()
        results = str_val.split("\t")
        scaffoldLength = int(results[1].strip())

        if(scaffoldLength > (check_len * 2)):
            scaffoldId = results[0].strip()
            bw1.write(str(scaffoldId) + "\t" + str(0) + "\t" + str(subjectStart) + "\n")
            bw2.write(str(scaffoldId) + "\t" + str(scaffoldLength - (check_len * 2) - subjectStart) + "\t" + str(scaffoldLength - (check_len * 2)) + "\n")
        br.close()
        bw1.close()
        bw2.close()

        cmd = str("cd " + outputDir + "\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed blastn-subject-2ndround.bed -fo blastn-subject-2ndround.fasta\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed blastn-query-2ndround.bed -fo blastn-query-2ndround.fasta\n" \
                  + "samtools faidx blastn-subject-2ndround.fasta\n" \
                  + "samtools faidx blastn-query-2ndround.fasta\n" \
                  + "makeblastdb -in blastn-subject-2ndround.fasta -dbtype nucl\n" \
                  + "blastn -query blastn-query-2ndround.fasta -db blastn-subject-2ndround.fasta -num_threads " + str(threads) + " -outfmt '7' -out blastn-res-2ndround.txt\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))
    return isCircular


def runAlignmentGetCoverage():
    print("runAlignmentGetCoverage:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    cmd = ""

    file_list = os.listdir(outputDir)
    has_bowtie2 = False
    has_samtools = False
    for file in file_list:
        path = os.path.join(outputDir, file)
        if "bowtie2" in path:
            has_bowtie2 = True
            break
    if has_bowtie2:
        cmd = "cd " + outputDir + "\n" + "rm -rf bowtie2*\n"

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    for file in file_list:
        path = os.path.join(outputDir, file)
        if "samtools" in path:
            has_samtools = True
            break
    if has_samtools:
        cmd = "cd " + outputDir + "\n" + "rm -rf samtools*\n"

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    isFasta = read1.replace('.gz','').endswith(('.fa','.fasta','.fna','.ffn','.frn'))
    bowtie2FastaFlag = " -f" if isFasta else ""

    if len(read2) == 0:
        cmd = str("cd " + outputDir + "\n" \
                  + "bowtie2-build --threads " + str(threads) + " scaffold.fasta bowtie2-index\n" \
                  + "bowtie2 --threads " + str(threads) + bowtie2FastaFlag + " -x bowtie2-index " \
                  + "-U " + read1 \
                  + " | samtools view -bS - | samtools view -h -F 0x04 -b - | " \
                  + "samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\n" \
                  + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\n")
    else:
        cmd = str("cd " + outputDir + "\n" \
                  + "bowtie2-build --threads " + str(threads) + " scaffold.fasta bowtie2-index\n" \
                  + "bowtie2 --threads " + str(threads) + bowtie2FastaFlag + " -x bowtie2-index " \
                  + "-1 " + read1 \
                  + " -2 " + read2 \
                  + " | samtools view -bS - | samtools view -h -F 0x04 -b - | " \
                  + "samtools sort -@ " + str(threads) + " - -o bowtie2-mapped.bam\n" \
                  + "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = "bash " + outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def percentile(values, percentiles):
    percentileResults = []
    temp = values
    values = values.sort()
    for percentile in percentiles:
        index = int(math.ceil((percentile / 100.00) * len(temp)))
        percentileResults.append(temp[index - 1])
    return percentileResults


def readCoverageGetPercentile():
    coverages = []
    br = open(outputDir + "/samtools-coverage.txt")
    str_val = br.readline()
    results = []
    while(str_val != ""):
        str_val = str_val.strip()
        results = str_val.split("\t")
        coverage = int(results[2])
        if coverage != 0:
            coverages.append(coverage)
        str_val = br.readline()
    br.close()

    percentiles = []
    percentiles.append(15.00)
    percentiles.append(85.00)
    percentileResults = percentile(coverages, percentiles)

    return percentileResults


def writeCoverageQuantile(quantile15Percent, quantile85Percent):
    print("writeCoverageQuantile:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    hasSuspiciousRegion = False
    scaffoldLength = 0
    scaffoldId = ""
    suspiciousStarts = []
    suspiciousEnds = []

    br = open(outputDir + "/scaffold.fasta.fai")
    str_val = br.readline()
    str_val = str_val.strip()
    results = []
    results = str_val.split("\t")
    scaffoldLength = int(results[1].strip())
    scaffoldId = results[0].strip()
    br.close()

    str_val = ""
    results = []
    coverage = 0
    startBase = 0
    currentBase = 0
    lowHigh = False

    br = open(outputDir + "/samtools-coverage.txt")
    str_val = br.readline()
    while(str_val != ""):
        str_val = str_val.strip()
        results = str_val.split("\t")
        coverage = int(results[2])
        currentBase = int(results[1])
        if (coverage >= quantile15Percent) and (coverage <= quantile85Percent):
            if lowHigh:
                if (currentBase - startBase) > minSuspiciousLen:
                    suspiciousStarts.append(startBase)
                    suspiciousEnds.append(currentBase-1)
            lowHigh = False
        else:
            if (startBase == 0) or (lowHigh == False):
                startBase = currentBase
            lowHigh = True
        str_val = br.readline()
    br.close()

    if len(suspiciousStarts) == 0:
        hasSuspiciousRegion = False
        bwCoverageOutputLog = open(outputDir + "/suspicious-regions.log",'w')
        bwCoverageOutputLog.write(scaffoldId + " has no suspicious region\n")
        bwCoverageOutputLog.close()
    else:
        hasSuspiciousRegion = True
        maxStartBase = 1
        maxLength = 0
        start = 0
        bwCoverageOutputLog = open(outputDir + "/suspicious-regions.log",'w')
        bwCoverageOutputLog.write(scaffoldId + " has " + str(len(suspiciousStarts)) + " suspicious regions\n")

        for i in range(len(suspiciousStarts)):
            start = suspiciousStarts[i]
            if maxLength == 0:
                maxLength = start - 1
            else:
                if (start - suspiciousEnds[i-1]) > maxLength:
                    maxLength = start - suspiciousEnds[i-1]
                    maxStartBase = suspiciousEnds[i-1] + 1

        if (scaffoldLength - suspiciousEnds[len(suspiciousEnds)-1]) > maxLength:
            maxLength = scaffoldLength - suspiciousEnds[len(suspiciousEnds)-1]
            maxStartBase = suspiciousEnds[len(suspiciousEnds)-1] + 1

        bwCoverageOutputLog.close()

        bw = open(outputDir + "/longest-non-suspicious.bed",'w')
        bw.write(scaffoldId + "\t" + str(maxStartBase) + "\t" + str(maxStartBase + maxLength - 1) + "\n")
        bw.close()

        cmd = str("cd " + outputDir + "\n" \
                  + "rm -rf longest-non-suspicious.fasta*\n" \
                  + "bedtools getfasta -fi scaffold.fasta -bed longest-non-suspicious.bed -fo longest-non-suspicious.fasta\n" \
                  + "samtools faidx longest-non-suspicious.fasta\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))
    return hasSuspiciousRegion


def updateScaffoldWithLongest():
    print("updateScaffoldWithLongest:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    scaffoldFile = outputDir + "/scaffold.fasta"
    if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        cmd = str("cd " + outputDir + "\n" \
                  + "rm -f scaffold.fasta*\n" \
                  + "cp longest-non-suspicious.fasta scaffold.fasta\n" \
                  + "samtools faidx scaffold.fasta\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def copyScaffoldForAssembly():
    print("copyScaffoldForAssembly:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    scaffoldFile = outputDir + "/scaffold.fasta"
    if (os.path.isfile(scaffoldFile) == True) and (os.path.isdir(scaffoldFile) == False):
        cmd = str("cd " + outputDir + "\n" \
                    + "rm -rf scaffold-truncated\n" \
                    + "mkdir scaffold-truncated\n" \
                    + "cp scaffold.fasta scaffold-truncated/scaffold.fasta\n" \
                    + "cd scaffold-truncated\n" \
                    + "samtools faidx scaffold.fasta\n")

        shellFileWriter = open(outputDir + "/run.sh",'w')
        shellFileWriter.write('#'+"!/bin/bash\n")
        shellFileWriter.write(cmd)
        shellFileWriter.close()

        cmd = "bash " + outputDir + "/run.sh"
        subprocess.check_output(cmd, shell=True)

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


def runCheckV(fasta_file):
    print("runCheckV:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    checkv_outdir = os.path.join(outputDir, "checkv_tmp")
    cmd = f"rm -rf {checkv_outdir} && checkv completeness {fasta_file} {checkv_outdir} -t {threads} -d {checkv_db}"

    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(f"[WARNING] CheckV execution failed. Skipping completeness check. Error: {e.output.decode()}")
        return False, 0

    tsv_file = os.path.join(checkv_outdir, "completeness.tsv")
    if not os.path.exists(tsv_file):
        print("[WARNING] CheckV completeness.tsv not found.")
        return False, 0

    is_complete = False
    expected_length = 0
    with open(tsv_file, 'r') as f:
        lines = f.readlines()
        if len(lines) > 1:
            header = lines[0].strip().split('\t')
            data = lines[1].strip().split('\t')

            try:
                idx_id = header.index("contig_id")
                idx_len = header.index("contig_length")
                idx_exp_len = header.index("aai_expected_length")
                idx_comp = header.index("aai_completeness")
            except ValueError:
                idx_id, idx_len, idx_exp_len, idx_comp = 0, 1, 3, 4

            # 提取预期长度
            if len(data) > idx_exp_len:
                exp_str = data[idx_exp_len]
                if exp_str != "NA" and exp_str != "Not-determined":
                    try:
                        expected_length = float(exp_str)
                    except ValueError:
                        expected_length = 0

            # 提取当前长度
            contig_length = 0
            if len(data) > idx_len:
                try:
                    contig_length = int(data[idx_len])
                except ValueError:
                    contig_length = 0

            # 提取完整度
            comp_val = 0
            if len(data) > idx_comp:
                comp_str = data[idx_comp]
                if comp_str != "NA" and comp_str != "Not-determined":
                    try:
                        comp_val = float(comp_str)
                    except ValueError:
                        comp_val = 0

            # 停摆条件1: 完整度 >= 90%
            if comp_val >= 90.0:
                is_complete = True
                msg = (
                    "\n" + "="*70 + "\n"
                    + "[SUCCESS] Assembly Complete! CheckV Completeness >= 90%\n"
                    + f"  - Contig ID:       {data[idx_id]}\n"
                    + f"  - Contig Length:   {contig_length} bp\n"
                    + f"  - Expected Length: {expected_length:.0f} bp\n"
                    + f"  - Completeness:    {comp_val}%\n"
                    + "="*70 + "\n"
                )
                print(msg)
                bwOutLog = open(os.path.join(outputDir, "output-log.txt"), 'a')
                bwOutLog.write(f"\nAssembly stopped: CheckV completeness reached {comp_val}%.\n")
                bwOutLog.write(f"Contig Length: {contig_length} bp, Expected Length: {expected_length:.0f} bp.\n")
                bwOutLog.close()

            # 停摆条件2: 长度 >= 2x 预期长度 (硬性上限，防止极端过度延伸)
            elif expected_length > 0 and contig_length >= expected_length * 2:
                is_complete = True
                msg = (
                    "\n" + "="*70 + "\n"
                    + "[WARNING] Assembly Stopped! Contig length exceeded 2x expected length\n"
                    + f"  - Contig ID:       {data[idx_id]}\n"
                    + f"  - Contig Length:   {contig_length} bp\n"
                    + f"  - Expected Length: {expected_length:.0f} bp\n"
                    + f"  - Completeness:    {comp_val}%\n"
                    + "="*70 + "\n"
                )
                print(msg)
                bwOutLog = open(os.path.join(outputDir, "output-log.txt"), 'a')
                bwOutLog.write(f"\nAssembly stopped: Contig length ({contig_length} bp) >= 2x expected length ({expected_length:.0f} bp).\n")
                bwOutLog.close()

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))
    return is_complete, expected_length


def checkCoverage():
    global run_counter
    print("checkCoverage:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    runAlignmentGetCoverage()
    percentiles = readCoverageGetPercentile()
    hasSuspiciousRegion = writeCoverageQuantile(percentiles[0], percentiles[1])

    bwOutLog = open(outputDir + "/output-log.txt",'a')
    is_complete = False

    if hasSuspiciousRegion:
        updateScaffoldWithLongest()
        bwOutLog.write("Scaffold got updated for suspicious region.\n")
        bwOutLog.close()
        copyScaffoldForAssembly()
        growScaffoldWithAssembly()
    else:
        bwOutLog.write("Scaffold didn't get updated as there is no suspicious region.\n")

        scaffoldFile = outputDir + "/scaffold.fasta"
        if os.path.isfile(scaffoldFile):
            backup_name = "run" + str(run_counter) + ".fasta"
            backup_path = os.path.join(outputDir, backup_name)
            cmd = str("cd " + outputDir + "\n" \
                      + "cp scaffold.fasta " + backup_name + "\n")

            shellFileWriter = open(outputDir + "/run.sh",'w')
            shellFileWriter.write('#'+"!/bin/bash\n")
            shellFileWriter.write(cmd)
            shellFileWriter.close()

            cmd = "bash " + outputDir + "/run.sh"
            subprocess.check_output(cmd, shell=True)

            bwOutLog.write("Saved current optimal scaffold as " + backup_name + " before truncation.\n")
            run_counter += 1

            # 触发 CheckV 完整度检测
            if checkv_db != "":
                bwOutLog.close()
                is_complete, _ = runCheckV(backup_path)
                bwOutLog = open(outputDir + "/output-log.txt",'a')

        bwOutLog.close()

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))
    return hasSuspiciousRegion, is_complete


def extendOneScaffold():
    print("extendOneScaffold:")
    print('Start time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))

    iteration = 1
    cmd = str("cd " + outputDir + "\n" \
              + "cp " + scaffold + " scaffold.fasta\n" \
              + "samtools faidx scaffold.fasta\n")

    shellFileWriter = open(outputDir + "/run.sh",'w')
    shellFileWriter.write('#'+"!/bin/bash\n")
    shellFileWriter.write(cmd)
    shellFileWriter.close()

    cmd = "bash " + outputDir + "/run.sh"
    subprocess.check_output(cmd, shell=True)

    extendContig = True
    prevLength = 0
    needUpdate = False

    while(extendContig):
        currentLength = getScaffoldLenFromTruncatedExtend()
        if currentLength > 300000:
            extendContig = False
        elif currentLength > prevLength:
            prevLength = currentLength
            updateCurrentScaffold()
            if checkCircularity():
                extendContig = False
            else:
                needUpdate, isComplete = checkCoverage()

                if isComplete:
                    extendContig = False
                    print("\n[INFO] Stopping extension process because CheckV confirmed genome completeness.\n")
                    break

                if needUpdate:
                    newLength = getScaffoldLenFromTruncatedExtend()
                    updateCurrentScaffold()
                    prevLength = newLength
                    getTruncatedScaffoldAndExtend(newLength)
                else:
                    getTruncatedScaffoldAndExtend(currentLength)

                # 检查 growScaffoldWithAssembly 内部 CheckV 是否触发停止
                if checkv_triggered:
                    checkv_triggered = False
                    extendContig = False
                    print("\n[INFO] Stopping extension process because CheckV triggered in growing scaffold.\n")
                    break

                iteration += 1
        else:
            extendContig = False

        if extendContig and iteration > 2000:
            extendContig = False

    print('End time: {0}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))


args = sys.argv[1:]
parseArguments(args)
print("Finished parsing input arguments")
getReadLen()
print("Started growing scaffold: {0}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))
extendOneScaffold()
print("Finished growing scaffold: {0}".format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')))
