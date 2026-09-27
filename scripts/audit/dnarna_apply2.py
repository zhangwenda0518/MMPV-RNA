#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""09 层补 DNA/RNA 列（字节级追加，原列零改动）v2

用法:  python3 dnarna_apply2.py            # 预演，只写 /tmp/dnarna_preview/
       python3 dnarna_apply2.py --apply    # 落地：备份 + 原地追加 3 列

追加列: Genome_Type (DNA/RNA) | Genome_Structure (VMR Genome 原文) | Genome_Source (species/genus/tax/none)
取值优先级: VMR Species 精确匹配 > VMR Genus 唯一取值 > 分类规则(Kingdom 优先于 Realm) > none

v2 相对 v1 的修正:
  1. 按字节处理，兼容 LF / CRLF / CR 三种行尾：追加内容一律插在行尾 \\r 之前，绝不落在 \\r 之后
  2. 校验改为双层：字节层(重算期望内容逐字节比对) + 解析层(独立 csv 解析后逐行逐字段比对原列)
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

def genome_fields(g):
    """DNA/RNA + VMR Genome 原文"""
    if not g:
        return "", ""
    gu = g.upper()
    return ("DNA" if "DNA" in gu else ("RNA" if "RNA" in gu else "")), g

KINGDOM_RULE = {"Orthornavirae": "RNA", "Pararnavirae": "DNA", "Loebvirae": "DNA", "Sangervirae": "DNA",
                "Shotokuvirae": "DNA", "Trapavirae": "DNA", "Bamfordvirae": "DNA", "Helvetiavirae": "DNA"}
REALM_RULE = {"Ribozyviria": "RNA", "Monodnaviria": "DNA", "Varidnaviria": "DNA",
              "Duplodnaviria": "DNA", "Adnaviria": "DNA"}

def tax_call(f, idx):
    def gv(n):
        return f[idx[n]].strip().strip('"') if n in idx and idx[n] < len(f) else ""
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

def annotate(raw, kind, stat):
    """按字节给每行追加 3 列，返回新字节串；行尾 \\r 保持在最后"""
    parts = raw.split(b"\n")
    out = []
    hdr = None
    for i, part in enumerate(parts):
        body, eol = (part[:-1], b"\r") if part.endswith(b"\r") else (part, b"")
        if not body.strip():
            out.append(part); continue
        text = body.decode("utf-8", "surrogateescape")
        try:
            f_ = next(csv.reader([text], delimiter="\t"))
        except Exception:
            f_ = text.split("\t")
        if hdr is None:
            hdr = [h.strip().strip('"') for h in f_]
            out.append(body + b"\t" + "\t".join(NEW_COLS).encode("utf-8") + eol)
            continue
        idx = {h: j for j, h in enumerate(hdr)}
        ncol = len(hdr)
        if len(f_) != ncol:
            stat["cols_mismatch"] += 1
        sp = f_[idx["Species"]].strip().strip('"') if kind == "full" and "Species" in idx else ""
        ge = f_[idx["Genus"]].strip().strip('"') if "Genus" in idx else ""
        gt = gs = ""; src = "none"
        if sp and sp in sp2g:
            gt, gs = genome_fields(sp2g[sp]); src = "species"
        elif ge and len(g2gs.get(ge, ())) == 1:
            gt, gs = genome_fields(next(iter(g2gs[ge]))); src = "genus"
        elif ge and len(g2gs.get(ge, ())) > 1:
            src = "genus_ambiguous"
        if not gt and kind == "full":
            t = tax_call(f_, idx)
            if t:
                gt = t; src = "tax" if src == "none" else src + "+tax"
        stat[src] += 1
        out.append(body + b"\t" + ("%s\t%s\t%s" % (gt, gs, src)).encode("utf-8") + eol)
    return b"\n".join(out), hdr

files = []
for rel, kind in TARGETS:
    for pat in (os.path.join(ROOT, "*", "02_novel_virus", "*_out", rel), os.path.join(ROOT, "*", "onekp-virus", rel)):
        files += [(p, kind) for p in sorted(glob.glob(pat))]

log_lines = ["# 09 层 DNA/RNA 补列记录 %s mode=%s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                                                "APPLY" if APPLY else "PREVIEW")]
print("%d 个目标文件  mode=%s\n" % (len(files), "APPLY" if APPLY else "PREVIEW"))
grand = Counter()
ok = 0

for p, kind in files:
    bak = p + ".bak_dnarna_" + STAMP
    src = bak if (APPLY and os.path.exists(bak)) else p
    raw = open(src, "rb").read()
    eol_kind = "CRLF" if b"\r\n" in raw else ("CR" if b"\r" in raw else "LF")
    stat = Counter()
    new_raw, hdr = annotate(raw, kind, stat)
    ncol = len(hdr)
    cnt = Counter(r.split(b"\t")[-3].decode("utf-8", "replace") for r in new_raw.split(b"\n")[1:] if r.strip(b"\r"))

    if APPLY:
        if not os.path.exists(bak):
            shutil.copy2(p, bak)
        before = md5(bak)
        with open(p, "wb") as f:
            f.write(new_raw)
        after = md5(p)
        # 校验 1：字节层，从备份重算期望内容逐字节比对
        exp, _ = annotate(open(bak, "rb").read(), kind, Counter())
        byte_ok = (exp == open(p, "rb").read())
        # 校验 2：解析层，独立 csv 解析，原列逐行逐字段一致
        def rows(q):
            with open(q, newline="", encoding="utf-8", errors="surrogateescape") as f:
                return list(csv.reader(f, delimiter="\t"))
        ro, rn = rows(bak), rows(p)
        parse_ok = (len(ro) == len(rn)) and all(a == b[:len(a)] for a, b in zip(ro, rn)) \
            and len(rn[0]) == ncol + 3 and rn[0][-3:] == NEW_COLS
        ok += 1 if (byte_ok and parse_ok) else 0
        rec = ("%s\n    eol=%s backup=%s md5 %s->%s 行=%d 列 %d->%d\n"
               "    species=%d genus=%d tax=%d ambiguous=%d none=%d | DNA=%d RNA=%d | 字段数异常行=%d\n"
               "    校验: 字节层=%s 解析层=%s" % (
                   p, eol_kind, bak, before, after, len(ro), ncol, ncol + 3,
                   stat["species"], stat["genus"], stat["tax"] + stat["species+tax"] + stat["genus+tax"],
                   stat["genus_ambiguous"], stat["none"], cnt["DNA"], cnt["RNA"], stat["cols_mismatch"],
                   "PASS" if byte_ok else "FAIL", "PASS" if parse_ok else "FAIL"))
    else:
        os.makedirs(OUTDIR, exist_ok=True)
        dst = os.path.join(OUTDIR, p.split("MMPV-paper/")[-1].replace("/", "__"))
        open(dst, "wb").write(new_raw)
        rec = ("%s\n    eol=%s -> %s 列 %d->%d\n    species=%d genus=%d tax=%d ambiguous=%d none=%d | DNA=%d RNA=%d | 字段数异常行=%d" % (
            p, eol_kind, dst, ncol, ncol + 3,
            stat["species"], stat["genus"], stat["tax"] + stat["species+tax"] + stat["genus+tax"],
            stat["genus_ambiguous"], stat["none"], cnt["DNA"], cnt["RNA"], stat["cols_mismatch"]))
    grand.update(stat)
    print(rec)
    log_lines.append(rec)

tot = sum(grand.values())
s = "\n合计 行=%d | species=%d(%.1f%%) genus=%d(%.1f%%) tax兜底=%d(%.1f%%) 属级歧义=%d 无值=%d(%.1f%%) | 双层校验全通过文件=%d/%d" % (
    tot, grand["species"], 100 * grand["species"] / tot, grand["genus"], 100 * grand["genus"] / tot,
    grand["tax"] + grand["species+tax"] + grand["genus+tax"],
    100 * (grand["tax"] + grand["species+tax"] + grand["genus+tax"]) / tot,
    grand["genus_ambiguous"], grand["none"], 100 * grand["none"] / tot, ok, len(files))
print(s); log_lines.append(s)
with open(LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines) + "\n")
print("\n记录: %s" % LOG)
