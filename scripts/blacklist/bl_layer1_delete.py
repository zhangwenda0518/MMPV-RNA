"""层1: 从各产物 FASTA 中剔除黑名单 contig。

策略: 流式读写, 只跳过 header first-token 在 kill 集合中的序列。
每文件处理后立即报告前后条数, 便于核对。
"""
import os
import shutil

KILL = "/tmp/blacklist_kill_ids.txt"
OUT = "/home/zhangwenda/data-test/out"

kill = set()
with open(KILL) as f:
    for l in f:
        l = l.strip()
        if l:
            kill.add(l)
print(f"剔除集: {len(kill)} 条\n")

TARGETS = [
    "06_HostPrediction/host_classified_fasta/Plant.classified.fasta",
    "07_Checkv/Plant.fasta",
    "08_Rescue/Plant/input_centroids.fasta",
    "08_Rescue/HQ_plant_viruses.fasta",
    "09_Virome_Analysis/HQ_plant_viruses.fasta",
    "09b_Analysis_Verify/class_KEEP.fasta",
    "09b_Analysis_Verify/acvirus_classify/contigs.fna",
    "09b_Analysis_Verify/virus_validation/keep_candidates.fasta",
    "09b_Analysis_Verify/virus_validation/review_candidates.fasta",
]

print(f"{'文件':<60} {'前':>7} {'后':>7} {'删':>6}")
print("-" * 85)
for rel in TARGETS:
    path = os.path.join(OUT, rel)
    if not os.path.isfile(path):
        print(f"{rel:<60} {'缺失':>7}")
        continue
    tmp = path + ".tmp_kill"
    n_before = n_after = 0
    with open(path, errors="replace") as fin, open(tmp, "w") as fout:
        keep = True
        for line in fin:
            if line.startswith(">"):
                cid = line[1:].split()[0]
                n_before += 1
                keep = cid not in kill
                if keep:
                    n_after += 1
                    fout.write(line)
            elif keep:
                fout.write(line)
    os.replace(tmp, path)
    print(f"{rel:<60} {n_before:>7} {n_after:>7} {n_before-n_after:>6}")

print("\n层1 完成。")
