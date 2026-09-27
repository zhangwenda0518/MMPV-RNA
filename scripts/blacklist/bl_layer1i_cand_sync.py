"""清理 09b virus_validation 候选文件中残留的 499 条 (与 scored.tsv 对齐)。

背景: scored.tsv 已从 9414 清到 8915 (KEEP2592/REVIEW2614/DROP3709)。
      但派生的候选文件未同步:
        keep_candidates.txt   2636 (多44)
        review_candidates.txt 3045 (多431)
        drop_candidates.txt   3733 (多24)
        drop_candidates.fasta 3733 (多24)  <- keep/review fasta 已同步
      清理 = 剔除不在 scored.tsv 中的 id。顺带核验 class_KEEP.fasta 是否 <= scored KEEP。
"""
import csv
import os
import shutil
import time

V = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/virus_validation"
ACV = "/home/zhangwenda/data-test/out/09b_Analysis_Verify"
STAMP = time.strftime("%Y%m%d")

# scored 真值 (按文件顺序)
scored_ids = []
verdict = {}
with open(f"{V}/rescue_evidence_scored.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        cid = r["contig_id"]
        scored_ids.append(cid)
        verdict[cid] = r["verdict"]
valid = set(scored_ids)
print(f"scored.tsv: {len(valid)} 条 (KEEP={sum(1 for v in verdict.values() if v=='KEEP')} "
      f"REVIEW={sum(1 for v in verdict.values() if v=='REVIEW')} "
      f"DROP={sum(1 for v in verdict.values() if v=='DROP')})")

# ── 1. 核验 class_KEEP.fasta ──
keep_ids = {c for c, v in verdict.items() if v == "KEEP"}
ck = f"{ACV}/class_KEEP.fasta"
ck_ids = set()
with open(ck, errors="replace") as f:
    for l in f:
        if l.startswith(">"):
            ck_ids.add(l[1:].split()[0])
print(f"\nclass_KEEP.fasta: {len(ck_ids)} 条; 全部∈scored KEEP? {ck_ids <= keep_ids} "
      f"(越界 {len(ck_ids - keep_ids)}); 被长度筛掉 {len(keep_ids) - len(ck_ids)}")

# ── 2. 清理 3 个 txt ──
print("\n=== 清理 txt ===")
for v in ["keep", "review", "drop"]:
    p = f"{V}/{v}_candidates.txt"
    with open(p, errors="replace") as f:
        lines = [l for l in f.read().split("\n") if l.strip()]
    keep_lines = [l for l in lines if l.split("\t")[0].strip() in valid]
    if len(keep_lines) != len(lines):
        shutil.copy2(p, f"{p}.bak_bl499_{STAMP}")
        with open(p, "w") as f:
            f.write("\n".join(keep_lines) + "\n")
    print(f"  {v}_candidates.txt: {len(lines)} -> {len(keep_lines)}  ({'已清' if len(keep_lines)!=len(lines) else '无需'})")

# ── 3. 清理 drop_candidates.fasta ──
print("\n=== 清理 fasta ===")
for v in ["keep", "review", "drop"]:
    p = f"{V}/{v}_candidates.fasta"
    if not os.path.exists(p):
        print(f"  {v}_candidates.fasta: 缺"); continue
    out_blocks, cur_hdr, cur_seq, n_in, n_out = [], None, [], 0, 0
    with open(p, errors="replace") as f:
        for l in f:
            if l.startswith(">"):
                if cur_hdr is not None:
                    out_blocks.append((cur_hdr, cur_seq))
                cur_hdr, cur_seq = l.rstrip("\n"), []
            else:
                cur_seq.append(l.rstrip("\n"))
        if cur_hdr is not None:
            out_blocks.append((cur_hdr, cur_seq))
    n_in = len(out_blocks)
    kept = [(h, s) for h, s in out_blocks if h[1:].split()[0] in valid]
    n_out = len(kept)
    if n_out != n_in:
        shutil.copy2(p, f"{p}.bak_bl499_{STAMP}")
        with open(p, "w") as f:
            for h, s in kept:
                f.write(h + "\n" + "\n".join(s) + "\n")
    print(f"  {v}_candidates.fasta: {n_in} -> {n_out}  ({'已清' if n_out!=n_in else '无需'})")

# ── 4. 复核 ──
print("\n=== 复核 ===")
for v in ["keep", "review", "drop"]:
    t = sum(1 for _ in open(f"{V}/{v}_candidates.txt"))
    fa = sum(1 for l in open(f"{V}/{v}_candidates.fasta", errors="replace") if l.startswith(">"))
    truth = sum(1 for x in verdict.values() if x == v.upper())
    print(f"  {v}: txt={t} fasta={fa} scored={truth} "
          f"{'OK' if t==truth else 'txt不符'} / {'OK' if fa<=truth else 'fasta超'}")
