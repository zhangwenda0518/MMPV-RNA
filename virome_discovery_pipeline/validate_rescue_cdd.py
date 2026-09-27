#!/usr/bin/env python3
"""
validate_rescue_cdd.py — rescue 后序列的 CDD 结构域验证（假病毒辅助去除）

用途：
    对 rescue 延伸后的序列（或 HQ_plant_viruses.fasta 全集）做 CDD 结构域验证，
    结合分类信息判断：
      1) 是否含有病毒结构域（T1 白名单命中）
      2) 是否含有其所属属的核心结构域（对照属级核心域表）
    输出证据等级，供假病毒去除 / 人工复核。

证据等级（EVIDENCE）：
    PASS_CORE       有 T1 病毒域 + 命中属级核心域      → 强证据，保留
    PASS_CORE_FAM   有 T1 病毒域 + 仅命中科级核心域    → 强证据（科级），保留
    PASS_VIRAL      有 T1 病毒域（属/科核心均未命中）  → 弱-中证据，保留但标注（可能为节段/片段）
    AMBIGUOUS       仅 T2 实证宿主蛋白或 T3 Root 白名单命中 → 待人工复核
    REVIEW          有分类但无 T1/T2 命中（延伸过头/短片段/CDD 盲区）→ 需人工复核
    NON_VIRAL       无任何白名单命中且无分类           → 疑似假阳性（宿主），建议去除

科级兜底：分类只有科（或属不在属级表中）时，用科级核心域表判定；
报告同时输出属级命中 (core_hits) 与科级命中 (fam_core_hits) 供对比。

v2 超家族归并集成（2026-08）：
    CDD 多标签效应会使同一功能家族的多模型命中被重复计数（实测膨胀 ~6×），
    且跨界同源超家族（RT/解旋酶/pol/HSP70）存在 viral/宿主双版本模型，
    制造虚假外源域信号（超家族双重投影）。
    本版在 convertalis 增加 qstart/qend（核酸坐标），统一聚类成位点后按成员构成分类：
      viral_only    仅病毒域
      dual          ext+病毒域同位点 = 双重投影伪影嫌疑
      proximal_dual ext 紧邻病毒域(±300nt同frame) = 区间错开伪影
      confirmed_capture HSP70 家族位点（文献支持真捕获）
      pure_clean    干净的真实外源捕获候选
    报告新增列：n_sites/site构成/suprafam_flag/pure_ext_domains。
    原 PASS_* 判定语义不变（向后兼容）；REVIEW 行若含 pure_clean 位点会在 note 列提示。

用法：
    python validate_rescue_cdd.py --input rescue.fasta --outdir out
                                  [--taxonomy consensus_tax.tsv]
                                  [--genus-table genera_core_aux_species_v3.tsv]
                                  [--threads 30] [--mode filter|report]
"""
import argparse, csv, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

# 默认路径（可被 CLI 覆盖）
DEF_CDD_DB      = Path(os.path.expanduser('~/database/cdd/cdd-db/cdd_db'))
DEF_CLASSIFIED  = Path(os.path.expanduser('~/MMPV-RNA/database/cdd/cdd_classified_taxid.tsv'))
DEF_WHITELIST   = Path(os.path.expanduser('~/MMPV-RNA/database/cdd/cdd_virus_final_v4.txt'))
DEF_GENUS_TABLE = Path(os.path.expanduser('~/plant_virus_db/3.final-ref-virus.db/genera_core_aux_species_v3.tsv'))
DEF_FAMILY_TABLE = Path(os.path.expanduser('~/plant_virus_db/3.final-ref-virus.db/family_core_aux_species_v3.tsv'))

def run(cmd, desc="", check=False):
    print(f"  [CMD] {cmd}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and r.returncode != 0:
        raise RuntimeError(f"{desc} failed:\n{r.stderr[-2000:]}")
    return r

# ---------- 1. 加载 CDD 分类 ----------
def load_classified(path):
    """cdd_classified_taxid.tsv → {accession: category}"""
    cat, title = {}, {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            cat[r['accession']] = r['category']
            title[r['accession']] = r.get('title', '')
    return cat, title

def load_whitelist(path):
    with open(path) as f:
        return {l.strip() for l in f if l.strip()}

# ---------- 2. 加载属级 / 科级核心域表 ----------
def load_genus_table(path):
    """genera_core_aux_species_v3.tsv → {genus: {core100:[], core90:[], aux:[]}}"""
    table = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            g = r['genus']
            table[g] = {
                'core100': set(x for x in r.get('core100_domains','').split(';') if x),
                'core90':  set(x for x in r.get('core90_domains','').split(';') if x),
                'aux':     set(x for x in r.get('aux_domains','').split(';') if x),
                'core_funcs': set(x for x in r.get('core_functions','').split(',') if x),
            }
    return table

def load_family_table(path):
    """family_core_aux_species_v3.tsv → {family: {core100:[], core90:[], aux:[]}} (列名兼容 family/n_species)"""
    table = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            fam = r.get('family', '') or r.get('Family', '')
            if not fam:
                continue
            table[fam] = {
                'core100': set(x for x in r.get('core100_domains','').split(';') if x),
                'core90':  set(x for x in r.get('core90_domains','').split(';') if x),
                'aux':     set(x for x in r.get('aux_domains','').split(';') if x),
                'core_funcs': set(x for x in r.get('core_functions','').split(',') if x),
            }
    return table

# ---------- 3. 加载分类信息（taxonomy 阶段产物） ----------
def load_taxonomy(path):
    """consensus_tax.tsv → {contig_id: {lineage, genus, family, ...}} (自动探测列名)"""
    tax = {}
    if not path or not Path(path).is_file():
        return tax
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            return tax
        # 识别 ID 列
        id_col = None
        for cand in ['contig_id', 'id', 'query', 'sequence_id', 'Sequence_ID', 'vOTU']:
            if cand in reader.fieldnames:
                id_col = cand; break
        if not id_col:
            id_col = reader.fieldnames[0]
        genus_col = next((c for c in ['genus', 'Genus', 'VMR_Genus', 'best_genus'] if c in reader.fieldnames), None)
        fam_col   = next((c for c in ['family', 'Family', 'VMR_Family', 'best_family'] if c in reader.fieldnames), None)
        for r in reader:
            row = dict(r)
            tax[row[id_col]] = {
                'genus':  row.get(genus_col, '') if genus_col else '',
                'family': row.get(fam_col, '') if fam_col else '',
                'lineage': row.get('lineage', row.get('Lineage', '')),
            }
    return tax

# ---------- 4. CDD 搜索 ----------
def run_cdd_search(input_fa, outdir, cdd_db, threads):
    """mmseqs 搜索并输出 hits.tsv (query,target,evalue)"""
    outdir.mkdir(parents=True, exist_ok=True)
    qdb = outdir / "queryDB"
    rdb = outdir / "resDB"
    tmp = outdir / "tmp"
    mmseqs = shutil.which('mmseqs') or 'mmseqs'
    run(f"{mmseqs} createdb {input_fa} {qdb} -v 1", "createdb", check=True)
    run(f"{mmseqs} search {qdb} {cdd_db} {rdb} {tmp} --search-type 2 -s 4 -v 1 --threads {threads} --max-seqs 300",
        "search", check=True)
    hits = outdir / "cdd_hits.tsv"
    run(f"{mmseqs} convertalis {qdb} {cdd_db} {rdb} {hits} --format-output query,target,evalue,qstart,qend -v 1",
        "convertalis", check=True)
    return hits

def parse_hits(hits_path, evalue=1e-3):
    """兼容 3 列(旧缓存)与 5 列(qstart,qend 核酸坐标)格式。
    返回: seq_hits {seq:set(acc)}, seq_hit_pos {seq:[(acc,ev,qs,qe)]}
    """
    seq_hits, seq_hit_pos = defaultdict(set), defaultdict(list)
    with open(hits_path) as f:
        header = f.readline()
        has_pos = len(header.rstrip('\n').split('\t')) >= 5
        f.seek(0)
        for l in f:
            p = l.strip().split('\t')
            if len(p) < 3: continue
            q, t, ev = p[0], p[1], float(p[2])
            if ev > evalue: continue
            seq_hits[q].add(t)
            if has_pos and len(p) >= 5:
                qs, qe = int(p[3]), int(p[4])
                seq_hit_pos[q].append((t, ev, min(qs, qe), max(qs, qe)))
    return seq_hits, seq_hit_pos

# ---------- 4b. 超家族归并（v2） ----------
def fam_of(acc, title):
    """ext 模型功能家族关键词映射"""
    t = (title.get(acc, '') or '').lower()
    if 'helicase' in t: return 'HELICASE'
    if 'reverse transcriptase' in t or t.startswith('rt_'): return 'RT'
    if 'hsp70' in t or 'heat shock' in t: return 'HSP70'
    if 'triphosphatase' in t or 'itp' in t or 'ditp' in t: return 'ITPASE'
    if 'methyltransferase' in t: return 'MTASE'
    if 'protease' in t: return 'PROTEASE'
    if 'polysaccharide' in t or 'caps_synth' in t: return 'CAPSULE'
    if 'polymerase' in t: return 'POLYMERASE'
    if 'photosystem' in t or 'psba' in t: return 'PHOTOSYS'
    return None

def is_ext(acc, cat, whitelist):
    return acc in whitelist and cat.get(acc, '') not in ('Virus', 'Root', '')

def is_viral(acc, cat, whitelist):
    return acc in whitelist and cat.get(acc) == 'Virus'

def cluster_sites(hit_pos):
    """同 reading frame 内, 与最强 E-value 锚区间重叠>=50%(较短者) 聚类。
    hit_pos: [(acc,ev,s,e)] 核酸坐标; 返回 [{s,e,members}]
    """
    by_frame = defaultdict(list)
    for acc, ev, s, e in hit_pos:
        by_frame[(s - 1) % 3].append((acc, ev, s, e))
    sites = []
    for fr, hs in by_frame.items():
        clusters = []
        for h in sorted(hs, key=lambda x: x[1]):
            acc, ev, s, e = h
            joined = False
            for cl in clusters:
                ov = min(e, cl['e']) - max(s, cl['s'])
                if ov > 0 and ov / min(e - s, cl['e'] - cl['s']) >= 0.5:
                    cl['members'].append(h); joined = True; break
            if not joined:
                clusters.append({'s': s, 'e': e, 'members': [h]})
        sites.extend(clusters)
    return sites

def classify_sites(hit_pos, cat, whitelist, title, prox=300):
    """位点三级分类; 返回 dict(n_viral_only,n_dual,n_proximal,n_confirmed,n_pure,
                              pure_fams set, has_suprafam bool) 或 None(无坐标)"""
    if not hit_pos:
        return None
    sites = []
    for cl in cluster_sites(hit_pos):
        ms = cl['members']
        ext_m = [m for m in ms if is_ext(m[0], cat, whitelist)]
        vir_m = [m for m in ms if is_viral(m[0], cat, whitelist)]
        sites.append({'s': cl['s'], 'e': cl['e'], 'ext': ext_m, 'vir': vir_m})
    # viral-only 区间索引(供 proximal 检测)
    vir_iv = defaultdict(list)
    for st in sites:
        if st['vir'] and not st['ext']:
            vir_iv[(st['s'] - 1) % 3].append((st['s'], st['e']))
    out = {'n_viral_only': 0, 'n_dual': 0, 'n_proximal': 0, 'n_confirmed': 0, 'n_pure': 0,
           'pure_fams': set(), 'dual_fams': set()}
    for st in sites:
        if not st['ext']:
            if st['vir']: out['n_viral_only'] += 1
            continue
        efams = {(fam_of(m[0], title) or f"DOM_{m[0]}") for m in st['ext']}
        if st['vir']:
            if 'HSP70' in efams:
                out['n_confirmed'] += 1   # 文献支持的真捕获
            else:
                out['n_dual'] += 1
                out['dual_fams'].update(efams)
        else:
            fr = (st['s'] - 1) % 3
            near = any(not (st['e'] < vs - prox or st['s'] > ve + prox) for vs, ve in vir_iv.get(fr, []))
            if near:
                out['n_proximal'] += 1
                out['dual_fams'].update(efams)
            else:
                out['n_pure'] += 1
                out['pure_fams'].update(efams)
    out['has_suprafam'] = (out['n_dual'] + out['n_proximal']) > 0
    return out

# ---------- 5. 判定 ----------
def classify_domain(acc, cat, whitelist):
    if acc in whitelist:
        return cat.get(acc, 'unclassified')
    return None  # 白名单外

def evaluate(seq_hits, cat, whitelist, genus_table, family_table, seq_tax,
             seq_hit_pos=None, title=None):
    """对每条序列输出判定行（v2: 附超家族归并位点统计，不改变原判定语义）

    判定优先级：
      1) 属级分类存在且属在表中 → 用属级核心域
      2) 科级分类存在且科在表中 → 用科级核心域（兜底）
      3) 属和科都命中 → PASS_CORE_FAM 标记（仅科级核心命中）
    """
    rows = []
    for q, doms in sorted(seq_hits.items()):
        # 白名单命中 → 分类
        wl = {a for a in doms if a in whitelist}
        t1 = {a for a in wl if cat.get(a) == 'Virus'}
        t2 = {a for a in wl if cat.get(a) not in ('Virus', 'Root')}
        t3 = {a for a in wl if cat.get(a) == 'Root'}
        # 分类信息
        info = seq_tax.get(q, {})
        genus = info.get('genus', '') or ''
        fam = info.get('family', '')
        has_tax = bool(genus or fam)
        # 属级核心域命中
        core_hit, aux_hit = '', ''
        hit_core = set()
        if genus and genus in genus_table:
            gt = genus_table[genus]
            hit_core = t1 & (gt['core100'] | gt['core90'])
            hit_aux = t1 & gt['aux']
            core_hit = ';'.join(sorted(hit_core))[:300]
            aux_hit = ';'.join(sorted(hit_aux))[:300]
        # 科级核心域命中（兜底，属级未命中时用）
        fam_core_hit, fam_aux_hit = '', ''
        fam_hit_core = set()
        if fam and fam in family_table:
            ft = family_table[fam]
            fam_hit_core = t1 & (ft['core100'] | ft['core90'])
            fam_hit_aux = t1 & ft['aux']
            fam_core_hit = ';'.join(sorted(fam_hit_core))[:300]
            fam_aux_hit = ';'.join(sorted(fam_hit_aux))[:300]
        # 证据等级：属核心优先，科核心兜底
        if t1 and hit_core:
            evid = 'PASS_CORE'          # 属级核心命中
        elif t1 and fam_hit_core:
            evid = 'PASS_CORE_FAM'      # 仅科级核心命中
        elif t1:
            evid = 'PASS_VIRAL'
        elif t2:
            evid = 'AMBIGUOUS'
        elif has_tax:
            evid = 'REVIEW'             # 有分类但无 T1/T2
        else:
            evid = 'NON_VIRAL'
        rows.append({
            'seq': q, 'genus': genus, 'family': fam,
            'n_t1': len(t1), 'n_t2': len(t2), 'n_t3': len(t3),
            't1_domains': ';'.join(sorted(t1))[:400],
            'core_hits': core_hit,
            'aux_hits': aux_hit,
            'fam_core_hits': fam_core_hit,
            'fam_aux_hits': fam_aux_hit,
            'evidence': evid,
            'suprafam_flag': '', 'site_summary': '', 'pure_ext_domains': '',
        })
        # v2: 位点级归并统计
        sc = classify_sites((seq_hit_pos or {}).get(q, []), cat, whitelist, title or {})
        if sc is not None:
            r = rows[-1]
            r['suprafam_flag'] = 'YES' if sc['has_suprafam'] else 'NO'
            r['site_summary'] = (f"vir={sc['n_viral_only']},dual={sc['n_dual']},"
                                 f"prox={sc['n_proximal']},hsp70cap={sc['n_confirmed']},pure={sc['n_pure']}")
            r['pure_ext_domains'] = ';'.join(sorted(sc['pure_fams']))[:200]
    return rows

# ---------- 主流程 ----------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', required=True, help='rescue 后序列 fasta')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--taxonomy', default=None, help='分类信息 TSV (contig_id → genus/family)')
    ap.add_argument('--genus-table', default=str(DEF_GENUS_TABLE))
    ap.add_argument('--family-table', default=str(DEF_FAMILY_TABLE), help='科级核心域表 (兜底, 属级未覆盖时用)')
    ap.add_argument('--cdd-db', default=str(DEF_CDD_DB))
    ap.add_argument('--classified', default=str(DEF_CLASSIFIED))
    ap.add_argument('--whitelist', default=str(DEF_WHITELIST))
    ap.add_argument('--threads', type=int, default=30)
    ap.add_argument('--evalue', type=float, default=1e-3)
    ap.add_argument('--mode', choices=['filter', 'report'], default='report',
                    help='filter: 同时输出 保留/去除 列表; report: 仅报告')
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("[1/5] 加载 CDD 分类 + 白名单 ...")
    cat, title = load_classified(args.classified)
    whitelist = load_whitelist(args.whitelist)
    print(f"      分类 {len(cat)} 条, 白名单 {len(whitelist)} 条")

    print("[2/5] 加载属级/科级核心域表 ...")
    genus_table = load_genus_table(args.genus_table)
    family_table = load_family_table(args.family_table) if args.family_table and Path(args.family_table).is_file() else {}
    print(f"      {len(genus_table)} 个属, {len(family_table)} 个科")

    print("[3/5] 加载分类信息 ...")
    seq_tax = load_taxonomy(args.taxonomy)
    print(f"      {len(seq_tax)} 条序列有分类")

    print("[4/5] CDD 搜索 ...")
    hits_file = outdir / 'cdd_hits.tsv'
    if hits_file.is_file():
        print("      复用已有 hits: " + str(hits_file))
    else:
        hits_file = run_cdd_search(args.input, outdir, args.cdd_db, args.threads)
    seq_hits, seq_hit_pos = parse_hits(hits_file, args.evalue)
    if not seq_hit_pos:
        print("      [警告] hits 无坐标列(旧缓存), 超家族归并降级为不可用")
    print(f"      {len(seq_hits)} 条序列有 CDD 命中")

    print("[5/5] 判定 (含超家族归并) ...")
    rows = evaluate(seq_hits, cat, whitelist, genus_table, family_table, seq_tax,
                    seq_hit_pos=seq_hit_pos, title=title)

    # 补充：无 CDD 命中的序列也进入报告 (证据=NON_VIRAL 或 REVIEW)
    from Bio import SeqIO
    seen = {r['seq'] for r in rows}
    for rec in SeqIO.parse(args.input, 'fasta'):
        q = rec.id
        if q in seen:
            continue
        info = seq_tax.get(q, {})
        has_tax = bool(info.get('genus', '') or info.get('family', ''))
        rows.append({
            'seq': q,
            'genus': info.get('genus', ''),
            'family': info.get('family', ''),
            'n_t1': 0, 'n_t2': 0, 'n_t3': 0,
            't1_domains': '', 'core_hits': '', 'aux_hits': '',
            'fam_core_hits': '', 'fam_aux_hits': '',
            'evidence': 'REVIEW' if has_tax else 'NON_VIRAL',
            'suprafam_flag': '', 'site_summary': '', 'pure_ext_domains': '',
        })

    # 报告 TSV
    report = outdir / 'cdd_evidence_report.tsv'
    with open(report, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['seq','genus','family','n_t1','n_t2','n_t3',
                                          't1_domains','core_hits','aux_hits',
                                          'fam_core_hits','fam_aux_hits','evidence',
                                          'suprafam_flag','site_summary','pure_ext_domains'],
                           delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\n报告: {report} ({len(rows)} 条)")

    # 汇总
    from collections import Counter
    cnt = Counter(r['evidence'] for r in rows)
    print("\n=== 证据等级汇总 ===")
    for k in ['PASS_CORE', 'PASS_CORE_FAM', 'PASS_VIRAL', 'AMBIGUOUS', 'REVIEW', 'NON_VIRAL']:
        print(f"  {k:<14} {cnt.get(k,0):>5}")
    total = len(rows)
    if total:
        core = cnt.get('PASS_CORE', 0) + cnt.get('PASS_CORE_FAM', 0) + cnt.get('PASS_VIRAL', 0)
        print(f"\n  保留建议: {core} ({100*core/total:.1f}%)  复核: {cnt.get('AMBIGUOUS',0)+cnt.get('REVIEW',0)}  去除建议: {cnt.get('NON_VIRAL',0)}")
    # v2 超家族归并汇总
    sf = [r for r in rows if r.get('site_summary')]
    if sf:
        import re as _re
        agg = Counter()
        for r in sf:
            for kv in r['site_summary'].split(','):
                k, v = kv.split('=')
                agg[k] += int(v)
        n_sf = sum(1 for r in rows if r.get('suprafam_flag') == 'YES')
        print(f"\n=== 超家族归并位点统计 ({len(sf)} 条序列, {sum(agg.values())} 个位点) ===")
        for k in ['vir', 'dual', 'prox', 'hsp70cap', 'pure']:
            print(f"  {k:<9} {agg.get(k,0):>6}")
        print(f"  含双重投影伪影嫌疑的序列: {n_sf}")
        pure_rows = [r for r in rows if r.get('pure_ext_domains')]
        if pure_rows:
            print(f"  真实捕获候选(pure_clean): {len(pure_rows)} 条 → 见报告 pure_ext_domains 列")

    if args.mode == 'filter':
        keep = outdir / 'keep_candidates.txt'
        review = outdir / 'review_candidates.txt'
        drop = outdir / 'drop_candidates.txt'
        with open(keep, 'w') as f:
            for r in rows:
                if r['evidence'] in ('PASS_CORE', 'PASS_CORE_FAM', 'PASS_VIRAL'):
                    f.write(r['seq'] + '\n')
        with open(review, 'w') as f:
            for r in rows:
                if r['evidence'] in ('AMBIGUOUS', 'REVIEW'):
                    f.write(r['seq'] + '\n')
        with open(drop, 'w') as f:
            for r in rows:
                if r['evidence'] == 'NON_VIRAL':
                    f.write(r['seq'] + '\n')
        print(f"  保留列表: {keep}  复核列表: {review}  去除列表: {drop}")

if __name__ == '__main__':
    main()
