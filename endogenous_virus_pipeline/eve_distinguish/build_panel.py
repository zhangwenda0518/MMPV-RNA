#!/usr/bin/env python3
"""S2 参考面板构建 v2: 从 NCBI protein 库拉取 Caulimoviridae (txid151433) RefSeq 蛋白,
defline 标题自带功能标签 (coat protein / movement protein / reverse transcriptase ...),
按关键词标注组件 (MP/CP/AP/RT/RH/GAG/VAP).
输出: panel.fasta, panel_report.tsv. 一次性运行; 需要联网."""
import os, sys, time, csv, re
from collections import Counter
from Bio import Entrez

EMAIL = "zhangwenda@example.com"
TAXID = 186534          # Caulimoviridae (NCBI taxonomy 实测)
RETMAX = 5000
PER_COMP_CAP = 60       # 每组件最多入库蛋白数, 控制面板大小

COMPONENT_RULES = [
    ("MP",  ["movement protein", "cell-to-cell transport"]),
    ("CP",  ["coat protein", "capsid protein"]),
    ("AP",  ["aspartic protease", "aspartyl protease", "aspartic proteinase"]),
    ("RT",  ["reverse transcriptase"]),
    ("RH",  ["ribonuclease h", "rnase h"]),
    ("GAG", ["gag-pol", "gag pol", "gag protein", "gag "]),
    ("VAP", ["virion-associated", "virion associated"]),
]
SKIP = ["hypothetical", "unknown", "partial", "fragment", "uncharacterized"]

COMP_ORDER = ["MP", "CP", "AP", "RT", "RH", "GAG", "VAP"]


def label(title):
    p = title.lower()
    if any(s in p for s in SKIP):
        return []
    tags = [c for c, kws in COMPONENT_RULES if any(k in p for k in kws)]
    # 顺序规整
    return [c for c in COMP_ORDER if c in tags]


def main(outdir):
    Entrez.email = EMAIL
    # 1) Swiss-Prot 优先 (标注最全), 不足再补 GenBank 全量
    ids = []
    for term in (f"txid{TAXID}[Organism] AND srcdb_swiss_prot[PROP]",
                 f"txid{TAXID}[Organism]"):
        try:
            h = Entrez.esearch(db="protein", retmax=RETMAX, term=term)
            batch = Entrez.read(h)["IdList"]
            h.close()
        except Exception as e:
            sys.stderr.write(f"[warn] esearch({term}): {e}\n")
            batch = []
        new = [i for i in batch if i not in ids]
        ids.extend(new)
        print(f"esearch[{term.split('AND')[-1].strip() or 'all'}] +{len(new)}")
        if ids and len(ids) >= RETMAX // 2:
            break
        time.sleep(1)
    print(f"蛋白 id 总数 {len(ids)}")
    # 2) 分批拉 fasta
    recs = []  # (defline, seq)
    B = 800
    for i in range(0, len(ids), B):
        batch = ids[i:i + B]
        for att in range(3):
            try:
                h = Entrez.efetch(db="protein", id=",".join(batch),
                                  rettype="fasta", retmode="text")
                txt = h.read(); h.close()
                break
            except Exception as e:
                sys.stderr.write(f"[warn] efetch {i}: {e}\n"); time.sleep(5)
        else:
            continue
        cur, buf = None, []
        for line in txt.splitlines():
            if line.startswith(">"):
                if cur:
                    recs.append((cur, "".join(buf)))
                cur, buf = line[1:], []
            else:
                buf.append(line.strip())
        if cur:
            recs.append((cur, "".join(buf)))
        time.sleep(0.4)
    print(f"efetch 取回 {len(recs)} 条")
    # 3) 标注 + 按组件限量
    fa_path = os.path.join(outdir, "panel.fasta")
    tsv_path = os.path.join(outdir, "panel_report.tsv")
    comp_count = Counter()
    rows = []
    for defline, seq in recs:
        # defline 形如 "YP_009123.1 reverse transcriptase [X virus]"
        m = re.match(r"(\S+)\s+(.*)", defline)
        if not m:
            continue
        pid, title = m.group(1), m.group(2)
        # 去掉尾部 [species] 前取标签
        lab_part = re.sub(r"\s*\[[^\]]*\]\s*$", "", title)
        tags = label(lab_part + " " + title.lower())
        if not tags:
            continue
        species = re.search(r"\[([^\]]+)\]\s*$", title)
        species = species.group(1) if species else "?"
        if any(comp_count[c] >= PER_COMP_CAP for c in tags) and len(tags) == 1:
            continue
        for c in tags:
            comp_count[c] += 1
        rows.append((pid, title, "+".join(tags), len(seq), species))
    with open(fa_path, "w") as fa, open(tsv_path, "w", newline="") as tsv:
        w = csv.writer(tsv, delimiter="\t")
        w.writerow(["protein_id", "title", "components", "prot_len", "species"])
        for pid, title, comp, ln, sp in rows:
            fa.write(f">{comp}|{pid}|{title}\n{seq_of(recs, pid)}\n")
            w.writerow([pid, title, comp, ln, sp])
    print(f"[done] 面板入库 {len(rows)} 条 -> {fa_path}")
    print("组件覆盖:", {k: comp_count[k] for k in COMP_ORDER if comp_count[k]})
    missing = [c for c in ("MP", "CP", "AP", "RT", "RH") if comp_count[c] < 5]
    if missing:
        print(f"[warn] 组件参考蛋白过少: {missing}")


_seq_cache = {}
def seq_of(recs, pid):
    if not _seq_cache:
        for dl, sq in recs:
            _seq_cache[dl.split()[0]] = sq
    return _seq_cache.get(pid, "")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
