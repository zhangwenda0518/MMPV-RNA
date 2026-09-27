import collections, os, sys

path = sys.argv[1]
DS = sys.argv[2] if len(sys.argv) > 2 else "RNA-Lycium_barbarum_out"
rows = []
with open(path) as f:
    f.readline()
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 5:
            rows.append(p)

print("total rows:", len(rows))

def top(rel):
    return rel.split("/")[0]

print("\n== by top dir ==")
for k, v in collections.Counter(top(r[0]) for r in rows).most_common():
    print(f"{v:6d}  {k}")

print("\n== by rule ==")
for k, v in collections.Counter(r[2] for r in rows).most_common():
    print(f"{v:6d}  {k}")

print("\n== rules within 02b_Filter (by relative depth-3 sample dir) ==")
sub = [r for r in rows if top(r[0]) == "02b_Filter"]
print("  02b_Filter broken:", len(sub))
for k, v in collections.Counter(r[2] for r in sub).most_common():
    print(f"{v:6d}  {k}")

print("\n== 02b_Filter links NOT matching the 4 rules ==")
n = 0
for r in sub:
    if r[2] in ("02a_all_candidate", "cdd_filtered", "uniprot_filtered", "strict_from_uniprot"):
        continue
    print(f"  {r[0]}\n     rule={r[2]}  old={r[1]}  new={r[3]}")
    n += 1
    if n >= 25:
        break
print("  ... shown", n)

print("\n== basename patterns in 02b_Filter ==")
pat = collections.Counter()
for r in sub:
    b = r[0].split("/")[-1]
    for suf in ("_virus.all.candidate.fasta", ".virus.candidate_cdd_filtered.fasta",
                ".virus.candidate_uniprot_filtered.fasta", ".virus.candidate_filtered.fasta"):
        if b.endswith(suf):
            pat[suf] += 1
            break
    else:
        pat["OTHER:" + b] += 1
for k, v in pat.most_common(20):
    print(f"{v:6d}  {k}")

print("\n== non-02b_Filter top dirs x basename-ish ==")
other = [r for r in rows if top(r[0]) != "02b_Filter"]
for k, v in collections.Counter(top(r[0]) for r in other).most_common(15):
    print(f"{v:6d}  {k}")
print("\n-- sample 3 old targets of non-02b --")
for r in other[:3]:
    print(f"  {r[0]}\n     old={r[1]}\n     new={r[3]}")
