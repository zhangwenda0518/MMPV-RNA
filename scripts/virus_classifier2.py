#!/usr/bin/env python3
# virus_classifier2.py — 病毒分类整合脚本 v4.2
# 支持: genomad, metabuli, CAT, diamond_lca, VITAP, mmseqs, ACVirus, vcontact3, PhaGCN3
# 输出: 8 级 combined_taxonomy.tsv (Realm Kingdom Phylum Class Order Family Genus Species)

import os, sys, argparse, subprocess, glob, time, copy, re
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
try: from tqdm import tqdm; HAS_TQDM = True
except ImportError: HAS_TQDM = False; tqdm = None
try: import psutil; HAS_PSUTIL = True
except ImportError: HAS_PSUTIL = False; psutil = None

RANK_NAMES = ["realm","kingdom","phylum","class","order","family","genus","species"]
HEADER = ["seq_name","tool","Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]

def safe_print(msg):
    if HAS_TQDM: tqdm.write(msg)
    else: print(msg)

def is_file_valid(path, min_size=1):
    return os.path.exists(path) and os.path.getsize(path) > min_size

def _sample_memory_peak(pid, stop_event, result_holder):
    if not HAS_PSUTIL: return
    peak = 0
    try:
        proc = psutil.Process(pid)
        while not stop_event.is_set():
            try:
                rss = proc.memory_info().rss
                for c in proc.children(recursive=True):
                    try: rss += c.memory_info().rss
                    except: pass
                if rss > peak: peak = rss
            except: break
            stop_event.wait(0.5)
    except: pass
    result_holder["peak_rss"] = peak

def run_command(cmd, log_file=None):
    try:
        if log_file:
            with open(log_file,'w') as f:
                subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)
        else:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
        return True, ""
    except subprocess.CalledProcessError as e:
        return False, f"code={e.returncode}"
    except Exception as e:
        return False, str(e)

# ==========================================================
# lineage → 8 级 rank 映射
# ==========================================================

def lineage_to_ranks(lineage_str):
    KNOWN_REALMS = {"Riboviria","Monodnaviria","Duplodnaviria","Varidnaviria","Adnaviria","Ribozyviria"}
    KNOWN_KINGDOMS = {"Orthornavirae","Shotokuvirae","Heunggongvirae","Lenarviricota"}
    SUBRANKS = ("viricotina","viricetidae","virineae","virinae","viricotina")
    skip = ("default","unplaced","unclassified","novel_subfamily","novel_order","cellular","root",
            "viruses","DNA viruses","dsDNA viruses","ssDNA viruses","RNA viruses","ssRNA viruses","dsRNA viruses",
            "unclassified phages")
    raw = []
    for p in lineage_str.split(";"):
        p = p.strip()
        if not p or p == "Viruses": continue
        p = p.replace("unclassified ","").replace("Unclassified ","")
        if raw and p == raw[-1]: continue  # 去相邻重复 ("Riboviria;unclassified Riboviria")
        raw.append(p)
    # 先过滤再丢弃 NA; 保留过滤后的有效段用于锚定 (避免 NA 占位推偏索引)
    parts_filtered = [p for p in raw if p != "-" and not any(w in p.lower() for w in skip)
                      and not any(p.endswith(s) for s in SUBRANKS)]
    ranks = {r:"NA" for r in RANK_NAMES}
    if not parts_filtered: return ranks
    # vcontact3: raw 有 8 段 (含 "-" 占位), 用 raw 按位置直接映射
    if len(raw) == 8:
        for i, rn in enumerate(RANK_NAMES):
            if raw[i] != "-": ranks[rn] = raw[i]
        return ranks
    parts = parts_filtered
    # 寻找 anchor: 已知realm/kingdom > -viricota(phylum) > -viricetes(class) > -idae(family) > -ales(order) > -virus(genus)
    anchor = None; rp = None
    for part_idx, (suffix_or_set, rank_name) in enumerate([
        (KNOWN_REALMS, "realm"), (KNOWN_KINGDOMS, "kingdom"),
        ("viricota","phylum"), ("viricetes","class"), ("idae","family"),
        ("ales","order"), ("viridae","family"), ("virus","genus")
    ]):
        if isinstance(suffix_or_set, set):
            for i, p in enumerate(parts):
                if p in suffix_or_set:
                    anchor = i; rp = RANK_NAMES.index(rank_name); break
        else:
            for i, p in enumerate(parts):
                if p.endswith(suffix_or_set):
                    anchor = i; rp = RANK_NAMES.index(rank_name); break
        if anchor is not None: break
    valid_parts = [p for p in parts if p!="NA"]
    # 特判: 只有2个有效段(Realm + species名)
    if len(valid_parts)==2 and valid_parts[0] in KNOWN_REALMS:
        ranks = {r:"NA" for r in RANK_NAMES}
        ranks["realm"] = valid_parts[0]
        ranks["species"] = valid_parts[1]
        return ranks
    if anchor is not None:
        for o, rn in enumerate(RANK_NAMES):
            ix = anchor + (o - rp)
            if 0 <= ix < len(parts) and parts[ix]!="NA":
                ranks[rn] = parts[ix]
    else:
        sub = parts[-6:] if len(parts)>=6 else parts
        tr = RANK_NAMES[-len(sub):]
        for i, p in enumerate(sub):
            if i < len(tr) and p!="NA": ranks[tr[i]] = p
    return ranks

# ==========================================================
# 分类工具
# ==========================================================

_DO = None
def _ensure_diamond_blastx(inp, s, out, uniprot_db, th):
    global _DO
    if _DO and _DO[1]: return _DO[0]
    d = os.path.join(out, "diamond_classify_output"); os.makedirs(d, exist_ok=True)
    bo = os.path.join(d, f"{s}_classify_blast.tsv")
    if is_file_valid(bo, 1000): _DO = (bo, True); return bo
    safe_print("  [diamond] 共享 blastx...")
    cmd = (f"diamond blastx -q '{inp}' --db '{uniprot_db}' --threads {th} "
           f"--more-sensitive --top 10 -e 0.001 "
           f"--outfmt 6 qseqid sseqid pident length mismatch gapopen "
           f"qstart qend sstart send evalue bitscore staxids -o '{bo}'")
    ok, _ = run_command(cmd, os.path.join(d, "diamond.log"))
    _DO = (bo, ok); return bo if ok else None

def classify_genomad(inp, s, out, db, th):
    d = os.path.join(out, "genomad_annotate_output"); os.makedirs(d, exist_ok=True)
    r = os.path.join(out, f"{s}_genomad_taxonomy.tsv")
    stem = os.path.splitext(os.path.basename(inp))[0]
    exp = os.path.join(d, f"{stem}_annotate", f"{stem}_taxonomy.tsv")
    if is_file_valid(r,10): return r
    cmd = f"genomad annotate --cleanup --full-ictv-lineage --lenient-taxonomy --threads {th} '{inp}' '{d}' '{db}'"
    ok, _ = run_command(cmd, os.path.join(out, "genomad_annotate.log"))
    if ok and os.path.exists(exp): os.system(f"cp '{exp}' '{r}' 2>/dev/null")
    return r

def classify_metabuli(inp, s, out, db, th):
    d = os.path.join(out, "metabuli_output"); os.makedirs(d, exist_ok=True)
    r = os.path.join(out, f"{s}_metabuli_taxonomy.tsv")
    if is_file_valid(r,10): return r
    ct = os.path.join(d, f"{s}_classifications.tsv")
    if not is_file_valid(ct,10):
        ok,_ = run_command(f"metabuli classify --seq-mode 3 --threads {th} '{inp}' '{db}' '{d}' '{s}'", os.path.join(out,"metabuli.log"))
        if not ok: return r
    if not is_file_valid(ct,10): return r
    tf = os.path.join(d, f"{s}_taxids.txt")
    os.system(f"awk '$1==1{{print $3}}' '{ct}' | sort -u > '{tf}' 2>/dev/null")
    if is_file_valid(tf,1):
        os.system(f"taxonkit lineage '{tf}' | awk -F'\\t' 'FNR==NR{{lin[$1]=$2; next}} $1==1 && $3 in lin && lin[$3] ~ /^Viruses;/{{print $2\"\\t1\\t1.0000\\t\"$3\"\\t\"lin[$3]}}' - '{ct}' > '{r}' 2>/dev/null")
    return r

def classify_cat(inp, s, out, cat_db, cat_tax, th):
    cat_dir = os.path.join(out, "cat_output"); os.makedirs(cat_dir, exist_ok=True)
    r = os.path.join(out, f"{s}_CAT_taxonomy.tsv")
    if is_file_valid(r, 10): return r
    # Step 1: CAT contigs
    ok, _ = run_command(
        f"CAT_pack contigs -c '{inp}' -o '{cat_dir}/CAT_output' -d '{cat_db}' -t '{cat_tax}' --nproc {th}",
        os.path.join(cat_dir, "CAT.log"))
    cf = os.path.join(cat_dir, "CAT_output.contig2classification.txt")
    if not ok or not is_file_valid(cf, 10): return r
    # Step 2: CAT add_names (去掉 --only_official)
    nf = os.path.join(cat_dir, "CAT_output.contig2classification.named.txt")
    ok2, _ = run_command(
        f"CAT_pack add_names -i '{cf}' -o '{nf}' -t '{cat_tax}' --exclude_scores",
        os.path.join(cat_dir, "CAT_add_names.log"))
    if not ok2 or not is_file_valid(nf, 10): return r
    # Step 3: 解析 "full lineage names" 列 (col 6), rank 格式: Name (rank)
    RANK_MAP = {"realm":"realm","kingdom":"kingdom","phylum":"phylum","class":"class",
                "order":"order","family":"family","genus":"genus","species":"species",
                "superkingdom":"realm","subfamily":"family","subgenus":"genus"}
    with open(nf) as f, open(r, 'w') as fo:
        fo.write("seq_name\ttaxid\tlineage\n")
        hdr = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7: continue
            # col 3: lineage (1;10239;...;185954*), col 4: lineage scores, col 5+: full lineage names (每 rank 一列)
            if "10239" not in parts[3]: continue
            seq_id = parts[0]
            taxid = parts[3].split(";")[-1].rstrip('*')
            # parts[5]=root(no rank), parts[6]=Viruses(acellular root), parts[7]=Riboviria(realm), ...
            rank_vals = {}
            for seg in parts[7:]:  # 跳过 root + Viruses
                m = re.match(r'(.+) \((.+)\)', seg)
                if m:
                    name, rk = m.group(1).strip(), m.group(2).strip().lower()
                    rk = RANK_MAP.get(rk, rk)
                    if name and "no rank" not in rk:
                        rank_vals[rk] = name
            # 按标准顺序输出: realm,kingdom,phylum,class,order,family,genus,species
            ranked = [rank_vals.get(rn, "-") for rn in RANK_NAMES]
            if any(v != "-" for v in ranked):
                fo.write(seq_id + "\t" + taxid + "\tViruses;" + ";".join(ranked) + "\n")
    return r


def classify_diamond_lca(inp, s, out, uniprot_db, th):
    d = os.path.join(out, "diamond_output"); os.makedirs(d, exist_ok=True)
    r = os.path.join(out, f"{s}_diamond_lca_taxonomy.tsv")
    if is_file_valid(r,10): return r
    lr = os.path.join(d, f"{s}_diamond_lca_raw.tsv")
    cmd = (f"diamond blastx --range-culling --top 10 -F 15 "
           f"-q '{inp}' --db '{uniprot_db}' --threads {th} "
           f"--outfmt 102 --include-lineage -o '{lr}'")
    ok, _ = run_command(cmd, os.path.join(d,"diamond_lca.log"))
    if ok and is_file_valid(lr,10):
        os.system(f"awk -F'\\t' '$4 ~ /^Viruses;/{{print $1\"\\t\"$2\"\\t\"$4}}' '{lr}' > '{r}' 2>/dev/null")
    return r

# ==========================================================
# 后处理: VITAP/mmseqs/ACVirus/vContact3/PhaGCN3 → standard
# ==========================================================

def postproc_mmseqs(inp, s, out):
    raw = os.path.join(out, "mmseqs_results", f"{s}_lca.tsv")
    if not is_file_valid(raw,10):
        for a in [f"{s}_taxonomy.tsv", f"{s}_lca.tsv"]:
            p = os.path.join(out, "mmseqs_results", a)
            if is_file_valid(p,10): raw = p; break
    r = os.path.join(out, f"{s}_mmseqs_taxonomy.tsv")
    if not is_file_valid(raw,10): return r
    with open(raw) as f, open(r,'w') as fo:
        fo.write("seq_name\ttaxid\tlineage\n")
        for ln in f:
            ps = ln.strip().split('\t')
            if len(ps)<9: continue
            sid, tid, lin = ps[0], ps[1], ps[8]
            if "cellular" in lin.lower() and "viruses" not in lin.lower() and "viria" not in lin.lower(): continue
            lp = []
            for p in lin.split(";"):
                p = p.strip()
                if not p: continue
                if p.startswith("-_"): p = p[2:]  # -_Riboviria → Riboviria (realm)
                elif len(p)>2 and p[1]=='_': p = p[2:]  # k_/p_/c_/o_/f_/g_/s_ 前缀
                if p and p != "Viruses": lp.append(p)
            fo.write(sid + "\t" + tid + "\t" + "Viruses;" + ";".join(lp) + "\n")
    return r

def postproc_vitap(inp, s, out):
    raw = os.path.join(out, "VITAP_results", f"{s}.vitap", "all_lineages.tsv")
    r = os.path.join(out, f"{s}_VITAP_taxonomy.tsv")
    if not is_file_valid(raw,10): return r
    seen = set()
    with open(raw) as f, open(r,'w') as fo:
        fo.write("seq_name\ttaxid\tlineage\n")
        next(f)
        for ln in f:
            ps = ln.strip().split('\t')
            if len(ps)<2: continue
            sid, lin = ps[0], ps[1]
            if sid in seen: continue; seen.add(sid)
            lps = [p for p in reversed(lin.split(";")) if p and p!="-"]
            fo.write(sid + "\t\t" + ("Viruses;"+";".join(lps) if lps else "Viruses") + "\n")
    return r

def postproc_acvirus(inp, s, out):
    raw = os.path.join(out, "ACVirus_results", f"{s}.acvirus", "final_result.tsv")
    r = os.path.join(out, f"{s}_ACVirus_taxonomy.tsv")
    if not is_file_valid(raw,10): return r
    with open(raw) as f, open(r,'w') as fo:
        fo.write("seq_name\ttaxid\tlineage\n")
        hdr = f.readline().strip().split('\t')
        for ln in f:
            ps = ln.strip().split('\t')
            if len(ps)<len(hdr): continue
            row = dict(zip(hdr, ps))
            sid = row.get('Nucleotide', ps[0])
            ranks = []
            for rk in ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]:
                v = row.get(rk,"").strip()
                if v and v not in ("-","NA","","no support"): ranks.append(v)
            fo.write(sid + "\t\t" + ("Viruses;"+";".join(ranks) if ranks else "Viruses") + "\n")
    return r

def postproc_vcontact3(inp, s, out):
    od = os.path.join(out, "vcontact3_results")
    # vContact3 标准输出: genome_by_genome_overview.csv (在输出根目录)
    raw = os.path.join(od, "genome_by_genome_overview.csv")
    if not is_file_valid(raw,10):
        cs = glob.glob(os.path.join(od,"**","*overview*.csv"), recursive=True) + \
             glob.glob(os.path.join(od,"**","final_assignments.csv"), recursive=True)
        raw = cs[0] if cs else raw
    r = os.path.join(out, f"{s}_vcontact3_taxonomy.tsv")
    if not is_file_valid(raw,10): return r
    with open(raw) as f, open(r,'w') as fo:
        fo.write("seq_name\ttaxid\tlineage\n")
        h = f.readline().strip().strip('#')
        hdr = [x.lower().strip('"') for x in h.split(',')]
        ic = next((i for i,x in enumerate(hdr) if x in ('genome','bin_name','contig','contig_id')),0)
        rc = next((i for i,x in enumerate(hdr) if x=='reference'),-1)
        rm = {}
        for rk in ["realm","kingdom","phylum","class","order","family","genus","species"]:
            for i,x in enumerate(hdr):
                if x==f"{rk}_prediction": rm[rk]=i; break
        skip_low = ("-","na","","unclassified","unknown","default")
        for ln in f:
            ps = ln.strip().split(',')
            if len(ps)<len(hdr): continue
            if rc>=0 and len(ps)>rc and ps[rc].strip().lower()!='false': continue
            sid = ps[ic].strip() if ic<len(ps) else ps[0]
            rv = []
            for rk in ["realm","kingdom","phylum","class","order","family","genus","species"]:
                idx = rm.get(rk)
                if idx is not None and idx<len(ps):
                    v = ps[idx].strip()
                    if v and v.lower() not in skip_low \
                       and "unplaced" not in v.lower() \
                       and not any(v.lower().startswith(w) for w in ("novel_genus","novel_subfamily","novel_family","novel_order")):
                        rv.append(v)
                    else: rv.append("-")
                else: rv.append("-")
            if any(v!="-" for v in rv):
                fo.write(sid + "\t\t" + "Viruses;" + ";".join(rv) + "\n")
    return r

def postproc_phagcn3(inp, s, out):
    od = os.path.join(out, "PhaGCN3_results"); os.makedirs(od, exist_ok=True)
    raw = os.path.join(od, f"{s}.phagcn3.csv")
    r = os.path.join(out, f"{s}_PhaGCN3_taxonomy.tsv")
    if not is_file_valid(raw,10): return r
    with open(raw) as f, open(r,'w') as fo:
        fo.write("seq_name\ttaxid\tlineage\n")
        hdr = f.readline().strip().split(',')
        for ln in f:
            ps = ln.strip().split(',')
            if len(ps)<2: continue
            sid = ps[0]
            found = {}
            for col in ps[1:]:
                for m in re.finditer(r'([rkpcofgs]);([^;]*)', col):
                    found[m.group(1)] = m.group(2)
            ranks = [found[k] for k in "rkpcofgs" if k in found and found[k] and found[k]!="unclassified"]
            fo.write(sid + "\t\t" + ("Viruses;"+";".join(ranks) if ranks else "Viruses") + "\n")
    return r

# ==========================================================
# 合并输出
# ==========================================================

def merge_taxonomy_results(sample, output_dir, tools_ran):
    combined = os.path.join(output_dir, f"{sample}_combined_taxonomy.tsv")
    rows = []
    tf_map = {
        "genomad": os.path.join(output_dir, f"{sample}_genomad_taxonomy.tsv"),
        "metabuli": os.path.join(output_dir, f"{sample}_metabuli_taxonomy.tsv"),
        "diamond_lca": os.path.join(output_dir, f"{sample}_diamond_lca_taxonomy.tsv"),
        "CAT": os.path.join(output_dir, f"{sample}_CAT_taxonomy.tsv"),
        "mmseqs": os.path.join(output_dir, f"{sample}_mmseqs_taxonomy.tsv"),
        "VITAP": os.path.join(output_dir, f"{sample}_VITAP_taxonomy.tsv"),
        "ACVirus": os.path.join(output_dir, f"{sample}_ACVirus_taxonomy.tsv"),
        "vcontact3": os.path.join(output_dir, f"{sample}_vcontact3_taxonomy.tsv"),
        "PhaGCN3": os.path.join(output_dir, f"{sample}_PhaGCN3_taxonomy.tsv"),
    }
    for tool in tools_ran:
        tf = tf_map.get(tool)
        if not tf or not os.path.exists(tf): continue
        with open(tf) as f: lines = f.readlines()
        has_hdr = bool(lines) and lines[0].startswith("seq_name")
        for line in lines[int(has_hdr):]:
            line = line.strip()
            if not line or line.startswith("#"): continue
            parts = line.split('\t')
            if len(parts) >= 5:
                r = lineage_to_ranks(parts[4])
                rows.append([parts[0], tool] + [r.get(rn,"NA") for rn in RANK_NAMES])
            elif len(parts) >= 3:
                r = lineage_to_ranks(parts[2])
                rows.append([parts[0], tool] + [r.get(rn,"NA") for rn in RANK_NAMES])
    with open(combined, 'w') as f:
        f.write('\t'.join(HEADER) + '\n')
        for row in rows:
            line = '\t'.join(row)
            if line.strip(): f.write(line + '\n')
    safe_print(f"  [合并] {len(rows)} 条 -> {os.path.basename(combined)}")
    return combined

# ==========================================================
# taxonkit 回填 NA
# ==========================================================

def save_resource_summary(out_dir, sample, metrics):
    usage_file = os.path.join(out_dir, f"{sample}_resource_usage.tsv")
    # 读已有记录, 只更新当前运行的工具 (避免覆盖其他工具)
    existing = {}
    if os.path.exists(usage_file):
        with open(usage_file) as f:
            for line in f.readlines()[1:]:
                ps = line.strip().split('\t')
                if len(ps) >= 5: existing[ps[1]] = ps
    for tool, m in metrics.items():
        if tool.startswith("_"): continue
        if m.get('wall_time_sec',0) > 0:  # 只保留实际运行的工具
            existing[tool] = [sample, tool,
                f"{m.get('wall_time_sec',0):.1f}",
                f"{m.get('cpu_time_sec',0):.1f}",
                f"{m.get('peak_rss_mb',0):.1f}", "OK"]
    with open(usage_file, 'w') as f:
        f.write("sample\ttool\twall_sec\tcpu_sec\tmem_mb\tstatus\n")
        for tool in sorted(existing):
            f.write('\t'.join(existing[tool]) + '\n')
    safe_print(f"  资源消耗: {os.path.basename(usage_file)}")


def fill_taxonomy_na(tsv_path, output_path):
    if not is_file_valid(tsv_path, 100): return
    with open(tsv_path) as f:
        lines = [l.rstrip('\r\n') for l in f.readlines()]
    hdr = lines[0].strip().split('\t')
    rank_cols = list(RANK_NAMES)
    ci = {}
    for r in rank_cols:
        for i, h in enumerate(hdr):
            if h.lower()==r.lower(): ci[r]=i; break
    skip_vals = {"NA","N/A","-","no rank","unknown","Unclassified","","default"}
    subranks = ("viricotina","viricetidae","virineae","virinae")
    skip_words = ("unplaced","novel_subfamily","novel_order","novel_genus","cellular","root",
                  "viruses","DNA viruses","dsDNA viruses","RNA viruses","unclassified phages")

    def _ok(v):
        if not v or v in skip_vals: return False
        return not any(w in v.lower() for w in skip_words)

    unames = set()
    to_fill = []
    for i, line in enumerate(lines[1:], 1):
        line = line.rstrip('\n')
        if not line.strip(): continue
        ps = line.split('\t')
        if len(ps)<len(hdr): continue
        if not any(not _ok(ps[ci[r]].strip()) for r in rank_cols if r in ci and ci[r]<len(ps)): continue
        dn = None
        for r in reversed(rank_cols):
            idx = ci.get(r)
            if idx is not None and idx<len(ps):
                if _ok(ps[idx].strip()): dn=ps[idx].strip(); break
        if dn: unames.add(dn); to_fill.append((i, dn))
        else: to_fill.append((i, None))
    if not unames:
        safe_print("  [fill] 无需回填 (0 NA)"); return
    safe_print(f"  [fill] {len(to_fill)} 行需回填, {len(unames)} 个唯一名, 批量 taxonkit...")
    r1 = subprocess.run("taxonkit name2taxid 2>/dev/null", input='\n'.join(unames), shell=True, capture_output=True, text=True)
    n2t = {}
    for ln in r1.stdout.strip().split('\n'):
        ps = ln.split('\t')
        if len(ps)>=2 and ps[1].isdigit(): n2t[ps[0]] = ps[1]
    tids = list(set(n2t.values()))
    if not tids:
        safe_print("  [fill] taxonkit name2taxid 无结果"); return
    r2 = subprocess.run("taxonkit lineage 2>/dev/null", input='\n'.join(tids), shell=True, capture_output=True, text=True)
    t2l = {}
    for ln in r2.stdout.strip().split('\n'):
        ps = ln.split('\t')
        if len(ps)>=2: t2l[ps[0]] = ps[1]
    n2r = {}
    for name, tid in n2t.items():
        l = t2l.get(tid,"")
        if not l: continue
        lp = []
        for p in l.split(";"):
            p = p.strip()
            if not p or p == "Viruses": continue
            p = p.replace("unclassified ","").replace("Unclassified ","")
            if any(w in p.lower() for w in skip_words): continue
            if any(p.endswith(s) for s in subranks): continue
            if lp and p == lp[-1]: continue
            lp.append(p)
        fill = {rn:"" for rn in rank_cols}
        if len(lp)==8:
            for i, rn in enumerate(rank_cols):
                if lp[i]!="NA" and lp[i] not in skip_vals: fill[rn]=lp[i]
        else:
            pi = None; rp = None
            KR = {"Riboviria","Monodnaviria","Duplodnaviria","Varidnaviria","Adnaviria","Ribozyviria"}
            KK = {"Orthornavirae","Shotokuvirae","Heunggongvirae","Lenarviricota"}
            for item, rn in [(KR,"realm"),(KK,"kingdom"),("viricota","phylum"),("viricetes","class"),
                              ("idae","family"),("ales","order"),("viridae","family"),("virus","genus")]:
                if isinstance(item, set):
                    for i, p in enumerate(lp):
                        if p in item: pi=i; rp=rank_cols.index(rn); break
                else:
                    for i, p in enumerate(lp):
                        if p.endswith(item): pi=i; rp=rank_cols.index(rn); break
                if pi is not None: break
            if pi is not None:
                for o, rn in enumerate(rank_cols):
                    ix = pi + (o - rp)
                    if 0<=ix<len(lp): fill[rn]=lp[ix]
            else:
                tail = lp[-6:] if len(lp)>=6 else lp
                tr = ["phylum","class","order","family","genus","species"][-len(tail):]
                for i, p in enumerate(tail):
                    if i<len(tr): fill[tr[i]]=p
        n2r[name]=fill
    fc = 0
    for row_idx, dn in to_fill:
        if not dn or dn not in n2r: continue
        fill = n2r[dn]
        ps = lines[row_idx].rstrip('\n').split('\t')
        for r in rank_cols:
            idx = ci.get(r)
            if idx is not None and idx<len(ps):
                if not _ok(ps[idx].strip()):
                    nv = fill.get(r,"")
                    if nv and nv != dn: ps[idx]=nv; fc+=1  # 不同名才填 (防止 genus→species 复制)
        lines[row_idx] = '\t'.join(ps)
    non_empty = [l for l in lines if l.strip()]
    with open(output_path, 'w') as f:
        f.write('\n'.join(non_empty) + '\n')
    safe_print(f"  [fill] {fc} 个 rank 回填完成")

# ==========================================================
# 主类
# ==========================================================

class VirusClassifier:
    def __init__(self, args, quiet_console=False, db_paths=None):
        self.args = args
        self.genomes = args.genomes
        self.sample = args.sample
        self.tools = args.tools
        self.threads = args.threads if hasattr(args,'threads') and args.threads else 20
        self.quiet_console = quiet_console
        self.db_paths = db_paths or {}
        self.output_dir = Path(args.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.resource_metrics = {}

    def print_progress(self, msg):
        if self.quiet_console: safe_print(msg)
        else: print(msg)

    def _run_cmd_with_resources(self, cmd, tool_name):
        m = {"wall_time_sec":0,"peak_rss_mb":None,"cpu_time_sec":None}
        try:
            ws = time.perf_counter()
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            se = threading.Event(); mr = {"peak_rss":0}
            mt = threading.Thread(target=_sample_memory_peak, args=(proc.pid, se, mr), daemon=True) if HAS_PSUTIL else None
            if mt: mt.start()
            stdout, stderr = proc.communicate()
            se.set()
            if mt: mt.join(timeout=1)
            m["wall_time_sec"] = round(time.perf_counter()-ws, 2)
            if HAS_PSUTIL and mr["peak_rss"]>0: m["peak_rss_mb"] = round(mr["peak_rss"]/(1024*1024), 2)
            self.resource_metrics[tool_name] = m
            if proc.returncode != 0:
                err = stderr[:500] if stderr else ""
                return False, f"rc={proc.returncode} {err}"
            return True, ""
        except Exception as e:
            self.resource_metrics[tool_name] = m
            return False, str(e)

    def run_classify(self, tool_name, func, *args):
        self.print_progress(f" [{tool_name}]...")
        ws = time.perf_counter()
        cs = psutil.Process().cpu_times() if HAS_PSUTIL else None
        se = threading.Event(); mr = {"peak_rss":0}
        mt = threading.Thread(target=_sample_memory_peak, args=(os.getpid(), se, mr), daemon=True) if HAS_PSUTIL else None
        if mt: mt.start()
        result = func(*args)
        se.set()
        if mt: mt.join(timeout=1)
        wall = round(time.perf_counter()-ws,2)
        m = {"wall_time_sec":wall}
        if HAS_PSUTIL and mr["peak_rss"]>0: m["peak_rss_mb"]=round(mr["peak_rss"]/(1024*1024),2)
        if HAS_PSUTIL and cs:
            try:
                ce = psutil.Process().cpu_times()
                m["cpu_time_sec"]=round((ce.user-cs.user)+(ce.system-cs.system),2)
            except: pass
        self.resource_metrics[tool_name]=m
        info = f"{wall:.0f}s"
        if m.get("peak_rss_mb"): info+=f", {m['peak_rss_mb']:.0f}MB"
        if m.get("cpu_time_sec"): info+=f", CPU:{m['cpu_time_sec']:.0f}s"
        self.print_progress(f"   {tool_name} ({info})")
        return result

    def run_vc_analysis(self):
        tools_ran = []
        t0 = time.time()
        out = str(self.output_dir)
        inp = self.genomes; s = self.sample
        uniprot = self.db_paths.get("uniprot","")

        def _skip(tool, db_path, tax_out, run_fn):
            if is_file_valid(tax_out,10) and not self.args.force:
                safe_print(f"  [{tool}] 已有结果, 跳过"); tools_ran.append(tool); return
            if db_path and not os.path.exists(db_path):
                safe_print(f"  [{tool}] DB 跳过: {db_path}"); return
            self.run_classify(tool, run_fn); tools_ran.append(tool)

        # virootaxonomy
        _skip("genomad", self.db_paths.get("genomad", os.path.expanduser("~/database/virus-db/genomad_db")),
              os.path.join(out, f"{s}_genomad_taxonomy.tsv"),
              lambda: classify_genomad(inp, s, out, self.db_paths.get("genomad", os.path.expanduser("~/database/virus-db/genomad_db")), self.threads))
        _skip("metabuli", self.db_paths.get("metabuli", os.path.expanduser("~/database/virus-db/RVDB-v31/RVDB_viroids.metabuli_db")),
              os.path.join(out, f"{s}_metabuli_taxonomy.tsv"),
              lambda: classify_metabuli(inp, s, out, self.db_paths.get("metabuli", os.path.expanduser("~/database/virus-db/RVDB-v31/RVDB_viroids.metabuli_db")), self.threads))
        _skip("diamond_lca", uniprot,
              os.path.join(out, f"{s}_diamond_lca_taxonomy.tsv"),
              lambda: classify_diamond_lca(inp, s, out, uniprot, self.threads))
        cat_db = self.db_paths.get("cat", os.path.expanduser("~/database/virus-db/RVDB-30/CAT-db/db"))
        cat_tax = self.db_paths.get("cat_tax", os.path.expanduser("~/database/virus-db/RVDB-30/CAT-db/tax"))
        _skip("CAT", cat_db,
              os.path.join(out, f"{s}_CAT_taxonomy.tsv"),
              lambda: classify_cat(inp, s, out, cat_db, cat_tax, self.threads))

        # ncbi-lca + ictv-network
        def _run_ski(tool, db, tax_out, pre_fn, post_fn):
            if is_file_valid(tax_out,10) and not self.args.force:
                safe_print(f"  [{tool}] 已有结果, 跳过"); tools_ran.append(tool); return
            if db and not os.path.exists(db):
                safe_print(f"  [{tool}] DB 跳过: {db}"); return
            def _r():
                pre_fn()
                post_fn(inp, s, out)
            self.run_classify(tool, _r); tools_ran.append(tool)

        vdb = self.db_paths.get("VITAP", os.path.expanduser("~/database/virus-db/vitap-db/VMR-MSL40_DB"))
        _run_ski("VITAP", vdb, os.path.join(out, f"{s}_VITAP_taxonomy.tsv"),
                 lambda: (Path(out,"VITAP_results").mkdir(exist_ok=True),
                          os.system(f"VITAP assignment -i {inp} -d {vdb} -p {self.threads} -o {Path(out,'VITAP_results')}/{s}.vitap > /dev/null 2>&1")),
                 postproc_vitap)

        mdb = self.db_paths.get("mmseqs", os.path.expanduser("~/database/virus-db/RVDB-30/RVDB.mmseqs"))
        if not os.path.exists(mdb):
            for alt in ["RVDB-30/RVDB.mmseqs", "RVDB-v31/RVDB.mmseqs_db"]:
                p = os.path.join(os.path.expanduser("~/database/virus-db"), alt)
                if os.path.exists(p): mdb = p; break
        _run_ski("mmseqs", mdb, os.path.join(out, f"{s}_mmseqs_taxonomy.tsv"),
                 lambda: (Path(out,"mmseqs_results").mkdir(exist_ok=True),
                          (Path(out,"mmseqs_results")/"tmp").mkdir(exist_ok=True),
                          os.system(f"mmseqs easy-taxonomy {inp} {mdb} {Path(out,'mmseqs_results')}/{s} {Path(out,'mmseqs_results')}/tmp --blacklist '' --tax-lineage 1 --threads {self.threads} --split-memory-limit 80G > /dev/null 2>&1")),
                 postproc_mmseqs)

        adb = self.db_paths.get("ACVirus", os.path.expanduser("~/database/virus-db/acvirus_db"))
        _run_ski("ACVirus", adb, os.path.join(out, f"{s}_ACVirus_taxonomy.tsv"),
                 lambda: (Path(out,"ACVirus_results").mkdir(exist_ok=True),
                          os.system(f"ACVirus classify --contig {inp} --data_path {adb} --out {Path(out,'ACVirus_results')}/{s}.acvirus > /dev/null 2>&1")),
                 postproc_acvirus)

        cdb = self.db_paths.get("vcontact3", os.path.expanduser("~/database/virus-db/vConTACT3_db"))
        _run_ski("vcontact3", cdb, os.path.join(out, f"{s}_vcontact3_taxonomy.tsv"),
                 lambda: os.system(f"vcontact3 run --nucleotide {inp} --output {Path(out,'vcontact3_results')} --db-version 232 --db-path {cdb} --threads {self.threads} --pyrodigal-gv --db-domain eukaryotes --export-all --keep-fna --keep-temp --exports cytoscape graphml profiles completeness centroids > /dev/null 2>&1"),
                 postproc_vcontact3)

        _run_ski("PhaGCN3", None, os.path.join(out, f"{s}_PhaGCN3_taxonomy.tsv"),
                 lambda: None,
                 lambda i,s2,o: (Path(o,"PhaGCN3_results").mkdir(exist_ok=True),
                                 postproc_phagcn3(i,s2,o) if is_file_valid(os.path.join(o,"PhaGCN3_results",f"{s2}.phagcn3.csv"),10) else None))

        if len(tools_ran)>=1:
            merged = merge_taxonomy_results(s, out, tools_ran)
            safe_print("  回填空缺 rank...")
            fill_taxonomy_na(merged, merged)
        wall_total = time.time() - t0
        save_resource_summary(out, s, self.resource_metrics)
        safe_print(f"\n[{s}] {len(tools_ran)}/{len(self.tools)} 工具, {wall_total:.0f}s")
        return True

# ==========================================================
# 入口
# ==========================================================

def process_single_wrapper(args_bundle):
    a, dp = args_bundle
    try: return VirusClassifier(a, quiet_console=True, db_paths=dp).run_vc_analysis()
    except Exception as e: print(f"\n致命: {a.sample}: {e}"); return False

def main():
    p = argparse.ArgumentParser(description="病毒分类整合脚本 v4.2 — 8级 taxonomy", formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="工具: genomad,metabuli,diamond_lca,VITAP,mmseqs,ACVirus,vcontact3,PhaGCN3,all")
    p.add_argument('-g','--genomes', help='FASTA')
    p.add_argument('-s','--sample', help='样品名')
    p.add_argument('-i','--input-dir', help='输入目录')
    p.add_argument('-e','--ext', default='.fasta', help='扩展名')
    p.add_argument('--remove-suffix', help='去后缀')
    p.add_argument('-j','--jobs', type=int, default=1, help='并行数')
    p.add_argument('-t','--tools', default='all', help='工具')
    p.add_argument('-o','--output-dir', default='./classify_output')
    p.add_argument('-p','--threads', type=int, default=20, help='线程')
    p.add_argument('-f','--force', action='store_true')
    p.add_argument('--db-dir', default=os.path.expanduser('~/database/virus-db'))
    p.add_argument('--genomad-db'); p.add_argument('--metabuli-db')
    p.add_argument('--cat-db'); p.add_argument('--cat-tax')
    p.add_argument('--uniprot-db')
    p.add_argument('--vitap-db'); p.add_argument('--mmseqs-db')
    p.add_argument('--acvirus-db'); p.add_argument('--vcontact3-db')
    args = p.parse_args()

    if not args.input_dir and (not args.genomes or not args.sample):
        p.error("需要 -i 或 -g + -s")

    all_tools = ["genomad","metabuli","CAT","diamond_lca","VITAP","mmseqs","ACVirus","vcontact3","PhaGCN3"]
    args.tools = all_tools if args.tools.lower()=='all' else [t.strip() for t in args.tools.split(',') if t.strip() in all_tools]

    db_paths = {
        "genomad": args.genomad_db or os.path.join(args.db_dir,"genomad_db"),
        "metabuli": args.metabuli_db or os.path.join(args.db_dir,"RVDB-v31","RVDB_viroids.metabuli_db"),
        "uniprot": args.uniprot_db or "",
        "cat": args.cat_db or os.path.join(args.db_dir,"RVDB-30","CAT-db","db"),
        "cat_tax": args.cat_tax or os.path.join(args.db_dir,"RVDB-30","CAT-db","tax"),
        "VITAP": args.vitap_db or os.path.join(args.db_dir,"vitap-db","VMR-MSL40_DB"),
        "mmseqs": args.mmseqs_db or os.path.join(args.db_dir,"RVDB-30","RVDB.mmseqs"),
        "ACVirus": args.acvirus_db or os.path.join(args.db_dir,"acvirus_db"),
        "vcontact3": args.vcontact3_db or os.path.join(args.db_dir,"vConTACT3_db"),
    }

    if args.input_dir:
        ip = Path(args.input_dir)
        if not ip.exists(): sys.exit(f"目录不存在: {ip}")
        files = list(ip.glob(f"*{args.ext}"))
        if not files: sys.exit(f"未找到 *{args.ext}")
        bo = Path(args.output_dir)
        tasks = []; skipped = []
        for f in files:
            sn = f.name
            if args.remove_suffix: sn = sn.replace(args.remove_suffix,'')
            elif args.ext: sn = sn.replace(args.ext,'')
            sf = bo / f"{sn}.virus_classed" / f"{sn}_combined_taxonomy.tsv"
            if sf.exists() and not args.force: skipped.append(sn)
            else: tasks.append((f, sn))
        print(f"批量: {len(files)} 文件, 跳过 {len(skipped)}, 需处理 {len(tasks)}")
        success = len(skipped)
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futures = {}
            for f, sn in tasks:
                la = copy.copy(args)
                la.genomes = str(f.absolute())
                la.sample = sn
                la.output_dir = str(bo / f"{sn}.virus_classed")
                futures[ex.submit(process_single_wrapper, (la, db_paths))] = sn
            it = as_completed(futures)
            if HAS_TQDM: it = tqdm(it, total=len(tasks), desc="进度", unit="样本")
            for fu in it:
                sn = futures[fu]
                if fu.result(): success += 1
        print(f"\n完成: {success}/{len(files)}")
    else:
        la = copy.copy(args)
        la.output_dir = os.path.join(args.output_dir, f"{args.sample}.virus_classed")
        VirusClassifier(la, quiet_console=False, db_paths=db_paths).run_vc_analysis()

if __name__ == "__main__":
    import threading
    main()
