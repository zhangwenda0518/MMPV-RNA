#!/usr/bin/env python3
"""Caulimoviridae 参考面板构建 [面板自身 v3, 独立于模块版本]

v3: bare "polyprotein" 兜底标注 + --supplement 属级补充模式。
  v2 面板 (222 条) 的组件标签全靠标题关键词, 而 Badnavirus 等属的 GenBank 蛋白大量
  只标 "polyprotein" —— 全部被 v2 的"无关键词即丢弃"规则扔掉。实测 v2 面板里
  Badnavirus 只有 1 条 VAP (VAP 还不在 s2 的 CANON 五组件里), 而该属是科内最大属、
  已报道植物 EVE 的主力类群 (Banana streak 系)。v3 两处改动:

    1. 标题不含任何组件关键词但含 "polyprotein" 的干净条目 (SKIP 词已在 label()
       里挡掉), 按科的 pol 区规范标 AP+RT+RH —— 蛋白酶+逆转录酶+RNase H 是该科
       pol ORF 的通用内容; MP/CP 有独立 ORF/条目, 不会因此虚标。
    2. --supplement 模式: **不重建面板** —— NCBI 内容漂移会让全量重建出来的面板
       和已提交的版本逐条不同, 全部下游数字跟着搬家。补充模式读入既有
       panel.fasta + panel_report.tsv, 抓科内蛋白, 跳过已收录 accession 与已覆盖
       物种, 只为面板里没有的物种按属级配额补条目 (Badnavirus 25 种、其他属 8 种;
       每物种每个组件组至多 1 条、取最长), 追加写入。既有条目一条不动, 已跑过的
       命中只增不减。

用法
----
  python3 build_panel.py <outdir>                 # 全量重建 (联网; v3 规则)
  python3 build_panel.py <outdir> --supplement    # 在 outdir 既有面板上补充 (联网)

输出: <outdir>/panel.fasta, panel_report.tsv;
      --supplement 另出 panel_genus_coverage.tsv (逐物种属级归属, 新旧条目都在)。
"""
import os
import re
import sys
import time
import csv
from collections import Counter, defaultdict

from Bio import Entrez

EMAIL = "zhangwenda@example.com"
TAXID = 186534          # Caulimoviridae (NCBI taxonomy 实测; 旧 docstring 写的
                        # txid151433 是过时数字, 186534 才是该科)
RETMAX = 30000
PER_COMP_CAP = 60       # 每组件最多入库蛋白数 (仅全量重建模式用)

# --supplement 的属级物种配额。Badnavirus 是主目标给大配额; 其他属只补存在感,
# 防止某个 GenBank 投放多的属把面板撑爆。
SUPPLEMENT_SPECIES_CAP = {"Badnavirus": 25}
SUPPLEMENT_SPECIES_CAP_DEFAULT = 8
PER_SPECIES_TAGSETS = 4     # 每物种最多几个不同组件组 (MP / CP / AP+RT+RH 通常 3 个)

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
# bare polyprotein 的兜底标注: 该科 pol ORF 的通用组件内容 (见模块 docstring)
POLYPROTEIN_COMPS = ["AP", "RT", "RH"]


def label(title):
    """None = 命中 SKIP 词 (硬跳过); [] = 无组件关键词 (v3 里可走 polyprotein 兜底)."""
    p = title.lower()
    if any(s in p for s in SKIP):
        return None
    tags = [c for c, kws in COMPONENT_RULES if any(k in p for k in kws)]
    return [c for c in COMP_ORDER if c in tags]


def tag_entry(title):
    """v3 标注: 关键词优先; 干净标题的 bare polyprotein 兜底标 AP+RT+RH."""
    tags = label(title)
    if tags is None:
        return None
    if tags:
        return tags
    if "polyprotein" in title.lower():
        return list(POLYPROTEIN_COMPS)
    return None


def _retry(fn, n=3, wait=5):
    for att in range(n):
        try:
            return fn()
        except Exception as e:
            sys.stderr.write(f"[warn] NCBI 第 {att + 1} 次失败: {e}\n")
            time.sleep(wait)
    return None


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def accession_of(pid):
    """sp|P03542.1|CAPSD_CAMVS -> P03542.1; 其余原样 (GenBank defline 首词即版本号)."""
    return pid.split("|")[1] if pid.startswith("sp|") else pid


def fetch_family_proteins():
    """全科蛋白元数据 -> [(pid, title)]; pid 为 defline 首词 (accession 或 sp|acc|NAME).

    只解析头行不搬序列: 科内蛋白上万条, 序列只对最终选中的几百条二次抓取。
    """
    h = Entrez.esearch(db="protein", retmax=RETMAX, term=f"txid{TAXID}[Organism]")
    ids = Entrez.read(h)["IdList"]
    h.close()
    print(f"esearch: 科内蛋白 {len(ids)} 条")
    recs, seen = [], set()
    for batch in _chunks(ids, 800):
        txt = _retry(lambda b=batch: Entrez.efetch(
            db="protein", id=",".join(b), rettype="fasta", retmode="text").read())
        if not txt:
            continue
        for line in txt.splitlines():
            if not line.startswith(">"):
                continue
            m = re.match(r"(\S+)\s+(.*)", line[1:])
            if m and m.group(1) not in seen:
                seen.add(m.group(1))
                recs.append((m.group(1), m.group(2).strip()))
        time.sleep(0.4)
    return recs


def fetch_seqs(pids):
    """按 accession 批量取序列 -> {accession: seq}.

    返回头首词对 Swiss-Prot 是 sp|ACC|NAME, 归一化回 ACC 再做键, 与请求侧
    accession_of() 对得上。
    """
    seqs = {}

    def key_of(tok):
        return tok.split("|")[1] if tok.startswith("sp|") else tok

    for batch in _chunks(pids, 800):
        txt = _retry(lambda b=batch: Entrez.efetch(
            db="protein", id=",".join(b), rettype="fasta", retmode="text").read())
        if not txt:
            continue
        acc, buf = None, []
        for line in txt.splitlines():
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc, buf = key_of(line[1:].split()[0]), []
            else:
                buf.append(line.strip())
        if acc:
            seqs[acc] = "".join(buf)
        time.sleep(0.4)
    return seqs


def taxid_of(accessions):
    """esummary protein: accession -> taxid."""
    out = {}
    for batch in _chunks(sorted(accessions), 200):
        summ = _retry(lambda b=batch: Entrez.read(
            Entrez.esummary(db="protein", id=",".join(b))))
        if not summ:
            continue
        for d in summ:
            acc = d.get("AccessionVersion") or d.get("Caption") or ""
            if d.get("TaxId"):
                # biopython 1.87 的 IntegerElement.__str__ 返回 repr 风格字符串
                # ("IntegerElement(3051985, ...)"), 必须先 int() 再 str()
                out[acc] = str(int(d["TaxId"]))
        time.sleep(0.34)
    return out


def genus_of_taxids(taxids):
    """efetch taxonomy: taxid -> 属名 (LineageEx 里 rank=genus 的节点)."""
    out = {}
    for batch in _chunks(sorted(taxids), 100):
        recs = _retry(lambda b=batch: Entrez.read(
            Entrez.efetch(db="taxonomy", id=",".join(b))))
        if not recs:
            continue
        for r in recs:
            genus = "?"
            for node in r.get("LineageEx", []):
                if node.get("Rank") == "genus":
                    genus = node.get("ScientificName", "?")
            out[str(int(r["TaxId"]))] = genus
    return out


def genus_of_species_names(names):
    """物种名 -> 属名 (esearch taxonomy 逐个查; 查不到记 '?')."""
    out = {}
    for n in sorted(names):
        h = _retry(lambda: Entrez.esearch(db="taxonomy", retmax=1,
                                          term=f"{n}[SCIN]"))
        ids = Entrez.read(h)["IdList"] if h else []
        if ids:
            g = genus_of_taxids([ids[0]])
            out[n] = g.get(ids[0], "?")
        else:
            out[n] = "?"
        time.sleep(0.34)
    return out


def read_report(tsv_path):
    with open(tsv_path, newline="", encoding="utf-8") as f:
        return list(csv.reader(f, delimiter="\t"))


def supplement(outdir):
    """在既有面板上补属级覆盖: 只为面板没有的物种加条目, 既有条目一条不动."""
    fa_path = os.path.join(outdir, "panel.fasta")
    tsv_path = os.path.join(outdir, "panel_report.tsv")
    rows = read_report(tsv_path)
    hdr, body = rows[0], rows[1:]
    assert hdr[:5] == ["protein_id", "title", "components", "prot_len", "species"], hdr
    n_fasta = sum(1 for line in open(fa_path, encoding="utf-8") if line.startswith(">"))
    if n_fasta != len(body):
        sys.exit(f"[supplement] panel.fasta ({n_fasta} 条) 与 panel_report.tsv "
                 f"({len(body)} 行) 对不上, 先核对再补")
    existing_pids = {r[0] for r in body}
    existing_species = {r[4] for r in body}
    print(f"[supplement] 既有面板 {len(body)} 条, {len(existing_species)} 个物种")

    recs = fetch_family_proteins()
    cand = {}
    for pid, title in recs:
        if pid in existing_pids or pid in cand:
            continue
        sp_m = re.search(r"\[([^\]]+)\]\s*$", title)
        if not sp_m:
            continue                        # 无 [物种] 的不补, 口径不明
        species = sp_m.group(1)
        if species in existing_species:
            continue                        # 只补面板没有的物种
        tags = tag_entry(title)
        if not tags:
            continue
        cand[pid] = {"pid": pid, "title": title, "tags": tags,
                     "comps": "+".join(tags), "species": species}
    print(f"[supplement] 新物种候选条目 {len(cand)} 条 "
          f"({len({e['species'] for e in cand.values()})} 个物种)")
    if not cand:
        print("[supplement] 没有需要新增的物种, 面板保持不变")
        return

    # 序列先取齐 (候选只有几百到一千多条), 属级选择按"最长序列"取代表
    seqs = fetch_seqs([accession_of(pid) for pid in cand])
    cand = {p: e for p, e in cand.items() if seqs.get(accession_of(p))}
    for e in cand.values():
        e["seq"] = seqs[accession_of(e["pid"])]

    t2s = taxid_of({accession_of(e["pid"]) for e in cand.values()})
    for e in cand.values():
        e["taxid"] = t2s.get(accession_of(e["pid"]), "")
    g2g = genus_of_taxids({e["taxid"] for e in cand.values() if e["taxid"]})
    for e in cand.values():
        e["genus"] = g2g.get(e["taxid"], "?") if e["taxid"] else "?"

    # 组内选择: 每物种每组件组取最长; 每物种至多 PER_SPECIES_TAGSETS 个组件组;
    # 属级物种配额 (Badnavirus 大配额, 其他属小配额); 物种按名字排序保证确定性.
    groups = defaultdict(lambda: defaultdict(list))
    for e in cand.values():
        groups[(e["genus"], e["species"])][e["comps"]].append(e)
    quota = defaultdict(int)
    chosen = []
    for (genus, species), by_comps in sorted(groups.items()):
        picked = sorted((max(v, key=lambda e: len(e["seq"])) for v in by_comps.values()),
                        key=lambda e: e["comps"])[:PER_SPECIES_TAGSETS]
        cap = SUPPLEMENT_SPECIES_CAP.get(genus, SUPPLEMENT_SPECIES_CAP_DEFAULT)
        if quota[genus] >= cap:
            continue
        quota[genus] += 1
        chosen.append((genus, species, picked))
    print("[supplement] 属级新增物种数: "
          + ", ".join(f"{g} {quota[g]}" for g in sorted(quota, key=lambda x: -quota[x])))
    if not chosen:
        print("[supplement] 配额内没有可新增的物种, 面板保持不变")
        return

    # 追加写 fasta + 报告
    new_rows = []
    for genus, species, picked in chosen:
        for e in sorted(picked, key=lambda e: (e["comps"], e["pid"])):
            with open(fa_path, "a", encoding="utf-8", newline="\n") as fa:
                fa.write(f">{e['comps']}|{e['pid']}|{e['title']}\n{e['seq']}\n")
            new_rows.append([e["pid"], e["title"], e["comps"], str(len(e["seq"])),
                             species])
    with open(tsv_path, "a", encoding="utf-8", newline="") as tsv:
        w = csv.writer(tsv, delimiter="\t", lineterminator="\n")
        w.writerows(new_rows)
    print(f"[supplement] 面板 {len(body)} -> {len(body) + len(new_rows)} 条")

    # 属级覆盖表: 既有物种 + 新增物种都给出属 (旧物种逐名查 taxonomy, 失败记 '?')
    old_g = genus_of_species_names(existing_species - {"?"})
    with open(os.path.join(outdir, "panel_genus_coverage.tsv"), "w",
              encoding="utf-8", newline="") as cov:
        w = csv.writer(cov, delimiter="\t", lineterminator="\n")
        w.writerow(["genus", "species", "source", "components"])
        for sp in sorted(old_g):
            comps = sorted({r[2] for r in body if r[4] == sp})
            w.writerow([old_g[sp], sp, "existing", "|".join(comps)])
        for genus, species, picked in chosen:
            w.writerow([genus, species, "added",
                        "|".join(sorted({e["comps"] for e in picked}))])
    n_added = Counter()
    for genus, species, picked in chosen:
        n_added[genus] += len(picked)
    n_old = Counter()
    for r in body:
        n_old[old_g.get(r[4], "?")] += 1
    print("[supplement] 属级条目数 existing+added:")
    for g in sorted(set(n_old) | set(n_added), key=lambda g: -(n_old[g] + n_added[g])):
        print(f"  {g:20s} {n_old[g]:4d} + {n_added[g]:3d}")


def rebuild(outdir):
    """全量重建 (v3 规则)。注意: 结果与已提交的 v2 面板不可逐条对齐, 主要用于新部署."""
    fa_path = os.path.join(outdir, "panel.fasta")
    tsv_path = os.path.join(outdir, "panel_report.tsv")
    recs = fetch_family_proteins()
    print(f"efetch 元数据 {len(recs)} 条, 按标题标注组件")
    comp_count = Counter()
    rows = []
    for pid, title in recs:
        lab_part = re.sub(r"\s*\[[^\]]*\]\s*$", "", title)
        title_all = lab_part + " " + title.lower()
        tags = tag_entry(title_all)
        if not tags:
            continue
        species_m = re.search(r"\[([^\]]+)\]\s*$", title)
        species = species_m.group(1) if species_m else "?"
        # v2 的组件上限只拦单组件条目 (保住 AP+RT 这类融合条目)。v3 的 polyprotein
        # 兜底条目都是 AP+RT+RH 三组件, 若同样放行会把 AP/RT/RH 各灌到几百条
        # (实测 448/468/438), 全是同质 pol 区 —— 兜底条目改按"任一组件满额即拦"。
        is_fallback = not label(title_all)
        if is_fallback and any(comp_count[c] >= PER_COMP_CAP for c in tags):
            continue
        if any(comp_count[c] >= PER_COMP_CAP for c in tags) and len(tags) == 1:
            continue
        for c in tags:
            comp_count[c] += 1
        rows.append((pid, title, "+".join(tags), species))
    print(f"标注通过 {len(rows)} 条, 取序列")
    seqs = fetch_seqs([accession_of(pid) for pid, *_ in rows])
    written = 0
    with open(fa_path, "w", encoding="utf-8", newline="\n") as fa, \
            open(tsv_path, "w", encoding="utf-8", newline="") as tsv:
        w = csv.writer(tsv, delimiter="\t", lineterminator="\n")
        w.writerow(["protein_id", "title", "components", "prot_len", "species"])
        for pid, title, comp, sp in rows:
            seq = seqs.get(accession_of(pid))
            if not seq:
                continue
            fa.write(f">{comp}|{pid}|{title}\n{seq}\n")
            w.writerow([pid, title, comp, len(seq), sp])
            written += 1
    print(f"[done] 面板入库 {written} 条 -> {fa_path}")
    print("组件覆盖:", {k: comp_count[k] for k in COMP_ORDER if comp_count[k]})
    missing = [c for c in ("MP", "CP", "AP", "RT", "RH") if comp_count[c] < 5]
    if missing:
        print(f"[warn] 组件参考蛋白过少: {missing}")


def main():
    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print(__doc__)
        return
    Entrez.email = EMAIL
    flags = [a for a in sys.argv[1:] if a.startswith("-")]
    positional = [a for a in sys.argv[1:] if not a.startswith("-")]
    outdir = positional[0] if positional else "."
    os.makedirs(outdir, exist_ok=True)
    if "--supplement" in flags:
        supplement(outdir)
    else:
        rebuild(outdir)


if __name__ == "__main__":
    main()
