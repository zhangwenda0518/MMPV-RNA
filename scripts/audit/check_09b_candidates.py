import csv
V = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/virus_validation"

sc = {}
with open(f"{V}/rescue_evidence_scored.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        sc[r["contig_id"]] = r["verdict"]
print(f"scored.tsv 真值: {len(sc)} 条")


def rd_txt(p):
    s = set()
    with open(p, errors="replace") as f:
        for l in f:
            l = l.strip()
            if l and l != "contig_id":
                s.add(l.split("\t")[0])
    return s


def rd_fa(p):
    s = set()
    with open(p, errors="replace") as f:
        for l in f:
            if l.startswith(">"):
                s.add(l[1:].split()[0])
    return s


for v in ["KEEP", "REVIEW", "DROP"]:
    low = v.lower()
    truth = {k for k, x in sc.items() if x == v}
    t = rd_txt(f"{V}/{low}_candidates.txt")
    fa = rd_fa(f"{V}/{low}_candidates.fasta")
    print(f"{v:<7} scored={len(truth):<5} txt={len(t):<5}(多{len(t-truth)}/少{len(truth-t)}) "
          f"fasta={len(fa):<5}(多{len(fa-truth)}/少{len(truth-fa)})")

print()
allbad = set()
for v in ["keep", "review", "drop"]:
    t = rd_txt(f"{V}/{v}_candidates.txt")
    allbad |= {k for k in t if sc.get(k) != v.upper()}
print(f"三 txt 合计与 scored 不符的 id: {len(allbad)}")
# 这些 id 是否还在 scored 里(改判) 还是已被彻底删除
in_sc = sum(1 for k in allbad if k in sc)
print(f"  其中仍在 scored 里(verdict 变了): {in_sc}")
print(f"  已被删除(不在 scored): {len(allbad)-in_sc}")
