#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""09 层补 DNA/RNA 列（追加式，原字节零改动）
用法:  python3 dnarna_apply.py            # 预演，只写 /tmp/dnarna_preview/
       python3 dnarna_apply.py --apply    # 落地：备份 + 原地追加 3 列

追加列（原列一行不改，只在行尾加）:
  Genome_Type      DNA / RNA
  Genome_Structure VMR Genome 原文（如 ssRNA(+)、dsDNA-RT、ssDNA）；分类兜底行为空
  Genome_Source    species / genus / tax / none（取值出处，便于追溯）

取值优先级: VMR Species 精确匹配 > VMR Genus 唯一取值 > 分类规则(Kingdom 先于 Realm) > none
"""
import csv, os, glob, sys, shutil, hashlib, datetime
from collections import Counter, defaultdict

APPLY = "--apply" in sys.argv
STAMP = "20260830"
VMR = os.path.expanduser("~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
ROOT = os.path.expanduser("~/MMPV-paper")
OUTDIR = "/tmp/dnarna_preview"
LOG = "/tmp/dnarna_apply_%s.log" % STAMP

TARGETS = [
    ("09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv", "full"),
    ("09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv", "full"),
    ("09_Virome_Analysis/all_plant_analysis/all_plant_viruses_genus_summary.tsv", "genus_only"),
]
NEW_COLS = ["Genome_Type", "Genome_Structure", "Genome_Source"]

# ---------- VMR 字典 ----------
sp2g, g2gs = {}, defaultdict(set)
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t"); next(rd)
    for r in rd:
        if len(r) < 27:
            continue
        g, genus, sp = r[25].strip().strip('"'), r[15].strip().strip('"'), r[17].strip().strip('"')
        if sp and g: sp2g.setdefault(sp, g)
        if genus and g: g2gs[genus].add(g)

def split_genome(g):
    """返回 (DNA/RNA, VMR Genome 原文)。原文比自造缩写更可追溯，不做二次解释。"""
    if not g:
        return "", ""
    gu = g.upper()
    t = "DNA" if "DNA" in gu else ("RNA" if "RNA" in gu else "")
    return t, g

# 分类规则：Kingdom 优先（Riboviria 含 Orthornavirae=RNA 与 Pararnavirae=RT-DNA，Realm 单独判会错）
KINGDOM_RULE = {"Orthornavirae": "RNA", "Pararnavirae": "DNA", "Loebvirae": "DNA", "Sangervirae": "DNA",
                "Shotokuvirae": "DNA", "Trapavirae": "DNA", "Bamfordvirae": "DNA", "Helvetiavirae": "DNA"}
REALM_RULE = {"Ribozyviria": "RNA", "Monodnaviria": "DNA", "Varidnaviria": "DNA",
              "Duplodnaviria": "DNA", "Adnaviria": "DNA"}

def tax_call(f, idx):
    def gv(name):
        return f[idx[name]].strip().strip('"') if name in idx and idx[name] < len(f) else ""
    k, rm = gv("Kingdom"), gv("Realm")
    if k in KINGDOM_RULE: return KINGDOM_RULE[k]
    if rm in REALM_RULE: return REALM_RULE[rm]
    return ""

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:12]

files = []
for rel, kind in TARGETS:
    for pat in (os.path.join(ROOT, "*", "02_novel_virus", "*_out", rel),
                os.path.join(ROOT, "*", "onekp-virus", rel)):
        for p in sorted(glob.glob(pat)):
            files.append((p, kind))

log_lines = ["# 09 层 DNA/RNA 补列记录  %s  mode=%s" % (
    datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "APPLY" if APPLY else "PREVIEW")]
print("%d 个目标文件  mode=%s\n" % (len(files), "APPLY" if APPLY else "PREVIEW"))
grand = Counter()

for p, kind in files:
    with open(p, encoding="utf-8", errors="surrogateescape", newline="") as f:
        raw = f.read()
    lines = raw.split("\n")
    hdr = next(csv.reader([lines[0]], delimiter="\t"))
    idx = {h.strip().strip('"'): i for i, h in enumerate(hdr)}
    if "Genome_Type" in idx:
        print("跳过(已有列): %s" % p); continue
    ncol = len(hdr)
    has_sp = "Species" in idx and kind == "full"
    has_ge = "Genus" in idx
    stat = Counter()
    newlines = [lines[0] + "\t" + "\t".join(NEW_COLS)]
    for line in lines[1:]:
        if not line.strip("\r"):
            newlines.append(line); continue
        try:
            f_ = next(csv.reader([line], delimiter="\t"))
        except Exception:
            f_ = line.split("\t")
        if len(f_) != ncol:
            stat["cols_mismatch"] += 1
        sp = f_[idx["Species"]].strip().strip('"') if has_sp and idx["Species"] < len(f_) else ""
        ge = f_[idx["Genus"]].strip().strip('"') if has_ge and idx["Genus"] < len(f_) else ""
        gt = gs = ""; src = "none"
        if sp and sp in sp2g:
            gt, gs = split_genome(sp2g[sp]); src = "species"
        elif ge and len(g2gs.get(ge, ())) == 1:
            gt, gs = split_genome(next(iter(g2gs[ge]))); src = "genus"
        elif ge and len(g2gs.get(ge, ())) > 1:
            src = "genus_ambiguous"
        if not gt and kind == "full":
            t = tax_call(f_, idx)
            if t:
                gt = t; src = "tax" if src == "none" else src + "+tax"
        stat[src] += 1
        grand[src] += 1
        newlines.append(line + "\t" + gt + "\t" + gs + "\t" + src)
    if stat["cols_mismatch"]:
        print("  !! 字段数不一致行=%d (见 %s)" % (stat["cols_mismatch"], p))
    out = "\n".join(newlines)

    if APPLY:
        bak = p + ".bak_dnarna_" + STAMP
        if not os.path.exists(bak):
            shutil.copy2(p, bak)
        before = md5(bak)
        with open(p, "w", encoding="utf-8", errors="surrogateescape", newline="") as f:
            f.write(out)
        after = md5(p)
        # 复核：行数一致 + 原列逐行逐字段等于备份
        with open(bak, encoding="utf-8", errors="surrogateescape", newline="") as f:
            old = f.read().split("\n")
        with open(p, encoding="utf-8", errors="surrogateescape", newline="") as f:
            new = f.read().split("\n")
        assert len(old) == len(new), "行数不一致 %s" % p
        bad = 0
        for a, b in zip(old[1:], new[1:]):
            if not a.strip("\r"):      # 末尾空行 / 空行，不参与前缀校验
                continue
            if not b.startswith(a + "\t"):
                bad += 1
        assert bad == 0, "原列被改动 %d 行 %s" % (bad, p)
        cnt = Counter(b.rsplit("\t", 3)[-3] for b in new[1:] if b.strip("\r"))
        rec = ("%s\n    backup=%s  md5=%s -> %s  行=%d  列 %d->%d\n    species=%d genus=%d tax=%d ambiguous=%d none=%d | DNA=%d RNA=%d | 原列校验=PASS" % (
            p, bak, before, after, len(new), ncol, ncol + 3,
            stat["species"], stat["genus"], stat["tax"] + stat["species+tax"] + stat["genus+tax"],
            stat["genus_ambiguous"], stat["none"], cnt["DNA"], cnt["RNA"]))
    else:
        os.makedirs(OUTDIR, exist_ok=True)
        dst = os.path.join(OUTDIR, p.split("MMPV-paper/")[-1].replace("/", "__"))
        with open(dst, "w", encoding="utf-8", errors="surrogateescape", newline="") as f:
            f.write(out)
        cnt = Counter(b.rsplit("\t", 3)[-3] for b in newlines[1:] if b.strip("\r"))
        rec = ("%s\n    -> %s  行=%d  列 %d->%d\n    species=%d genus=%d tax=%d ambiguous=%d none=%d | DNA=%d RNA=%d" % (
            p, dst, len(newlines), ncol, ncol + 3,
            stat["species"], stat["genus"], stat["tax"] + stat["species+tax"] + stat["genus+tax"],
            stat["genus_ambiguous"], stat["none"], cnt["DNA"], cnt["RNA"]))
    print(rec)
    log_lines.append(rec)

tot = sum(grand.values())
log_lines.append("\n合计 行=%d | species=%d(%.1f%%) genus=%d(%.1f%%) tax兜底=%d(%.1f%%) 属级歧义=%d 无值=%d(%.1f%%)" % (
    tot, grand["species"], 100 * grand["species"] / tot, grand["genus"], 100 * grand["genus"] / tot,
    grand["tax"] + grand["species+tax"] + grand["genus+tax"],
    100 * (grand["tax"] + grand["species+tax"] + grand["genus+tax"]) / tot,
    grand["genus_ambiguous"], grand["none"], 100 * grand["none"] / tot))
print("\n" + log_lines[-1])
with open(LOG, "a", encoding="utf-8") as f:
    f.write("\n".join(log_lines) + "\n\n")
print("\n记录: %s" % LOG)
