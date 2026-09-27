#!/usr/bin/env python3
"""
gen_final_judgement.py — 按物种生成 final_judgement_table.tsv (七工具分类 + blast 身份)
供 integrate_rescue_evidence.py 五层证据整合使用.

final_judgement_table.tsv 11 列:
  contig_id, host_species, length, prevalence_%,
  nt_pident, nt_qcovs, nt_species,          <- blastn vs 核苷酸病毒库
  aa_pident, aa_species,                    <- blastx/mmseqs vs 蛋白病毒库
  cdd_top, category

blast 库默认:
  nt:  ncbi-virus/ncbi_virus.nucl.fasta   (171K 核苷酸病毒, 含新株)
  aa:  ncbi-virus/ncbi_virus.prot.fasta   (367K 蛋白病毒, 配套)
可用其它库覆盖. category 按 CDD 结构域分类表 + blast + 长度推导.

重要: category 是【展示性标签】(供 report 可视化 / 论文统计/人为审查),
【绝不参与】integrate_rescue_evidence 的 verdict(KEEP/REVIEW/DROP) 裁决.
假病毒去除的裁决只看: blast 命中 + CDD 病毒域证据 + 跨样本检测频次.
cdd_classified_taxid.tsv 判定的是结构域保守于哪类生物(Plant/Animal/...),
只是“线索”(指示可能的宿主/Contamination来源), 不是宿主身份的硬证据;
若某序列无病毒 CDD 域 + blast 弱, 即使不用该表也会归入非保留类,
其 verdict 不受 category 影响.

用法:
  python gen_final_judgement.py --fasta HQ_plant_viruses.fasta \
      --freq frequency_table.tsv --cdd cdd_hits.tsv \
      --cdd-evidence cdd_evidence_report.tsv \
      --sample barbarum --out final_judgement_table.tsv \
      [--nt-db ...] [--aa-db ...] [--threads 32]
"""
import argparse, csv, os, subprocess, sys, tempfile
from pathlib import Path
from collections import defaultdict

# Known virus判定阈值 (与 integrate_rescue_evidence 对齐)
NT_KNOWN = 85.0
NT_RELATED = 40.0
AA_RELATED = 40.0
# query coverage (nt_qcovs) 阈值: 物种级(已知种)对齐主流85%, 属/新种放宽到50%
NT_QCOV_KNOWN = 85.0
NT_QCOV_RELATED = 50.0
SHORT_LEN = 500      # short fragment / 低覆盖 阈值
DEF_CDD_CLASSIFIED = os.environ.get('MMPV_CDD_TAXID_TSV', str(
    Path(__file__).resolve().parent.parent / 'database' / 'cdd' / 'cdd_classified_taxid.tsv'))

VIRUS_KEYWORDS = ['phage', 'phiv', 'bacteriophage']          # 噬菌体
INSECT_KEYWORDS = ['insect', 'aphis', 'arbovirus', 'mosquito', 'bee', 'wasp', 'aphid', 'louse', 'moth', 'butterfly', 'planthopper', 'leafhopper']
FUNGAL_KEYWORDS = ['fungus', 'fungi', 'mycovirus', 'ascomycete', 'basidiomycete', 'hypovirus', 'totivirus', 'partitivirus']


def run(cmd, desc="", check=True):
    print(f"  ▶ {desc}: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    r = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True)
    if check and r.returncode != 0:
        sys.exit(f"[{desc}] 失败: {r.stderr[-2000:] if r.stderr else r.stdout[-2000:]}")
    return r


def safe_float(x, default=0.0):
    try:
        return float(x) if x not in ('', None) else default
    except (ValueError, TypeError):
        return default


def parse_length_from_id(cid):
    """contig_2246_length_... 或 NODE_..._length_... → int, 失败返回 ''"""
    import re
    m = re.search(r'length[=_]([0-9]+)', cid)
    return int(m.group(1)) if m else ''


def classify_host(species_name):
    """按 blast 命中物种名判断来源类型 → plant/phage/insect/fungus/other"""
    s = (species_name or '').lower()
    if not s:
        return 'unknown'
    if any(k in s for k in VIRUS_KEYWORDS):
        return 'phage'
    if any(k in s for k in INSECT_KEYWORDS):
        return 'insect'
    if any(k in s for k in FUNGAL_KEYWORDS):
        return 'fungus'
    return 'plant'   # 其余病毒视为植物/未分类病毒


def load_fasta_ids(fasta):
    ids, lens = [], {}
    with open(fasta) as f:
        cid = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                cid = line[1:].split()[0]
                ids.append(cid)
            elif cid and line:
                lens[cid] = lens.get(cid, 0) + len(line)
    for cid in ids:
        lens.setdefault(cid, '')
    return ids, lens


def load_freq(path):
    """frequency_table.tsv → {contig_id: prevalence_%, ...}"""
    d = {}
    if not path or not Path(path).is_file():
        return d
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            d[r['contig_id']] = r.get('prevalence_%', '')
    return d


def load_cdd_top(path):
    """cdd_hits.tsv → {contig_id: top domain(s) 分号拼接}. 取每序列前3个 distinct 域"""
    d = defaultdict(list)
    if not path or not Path(path).is_file():
        return {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            q, dom = r.get('qseqid', r.get('contig_id', '')), r.get('sacc', r.get('cd_accession', ''))
            # cdd_hits 是 <queryid, cdaccession, evalue> 三列
            if not dom and len([x for x in r]) >= 2:
                vals = list(r.values())
                q, dom = vals[0], vals[1]
            if q and dom and dom not in d[q]:
                d[q].append(dom)
    return {q: ';'.join(doms[:5]) for q, doms in d.items()}


def load_cdd_evidence(path):
    """cdd_evidence_report.tsv → {contig_id: evidence}"""
    d = {}
    if not path or not Path(path).is_file():
        return d
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if 'seq' in r and 'evidence' in r:
                d[r['seq']] = r['evidence']
    return d


def blastn_query(fasta, nt_db, threads, out_tsv):
    """blastn vs 核苷酸库 → 写 out_tsv (outfmt6). 返回 {q: (pident, qcov, sseqid, stitle)}"""
    r = run(['blastn', '-query', str(fasta), '-db', str(nt_db),
             '-out', out_tsv, '-outfmt', '6 qseqid sseqid pident qcovs evalue bitscore stitle',
             '-num_threads', str(threads), '-max_target_seqs', '1',
             '-evalue', '1e-3'],
            'blastn')
    best = {}
    with open(out_tsv) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 7:
                continue
            q = p[0]
            if q not in best:
                best[q] = (float(p[2]), float(p[3]), p[1], p[6])
    return best


def mmseqs_query(fasta, aa_db, threads, workdir):
    """diamond blastx vs 蛋白病毒库 → 返回 {q: (aa_pident, sseqid, subject_cov)}.
    aa_db 为蛋白 fasta; 内部用 diamond 的 .dmnd 索引.
    subject_cov = 比对对象(subject 蛋白)被查询覆盖的比例(合并同蛋白所有 HSP 区间/蛋白全长).
    判断"这条 contig 是否编码了某个病毒蛋白"用 subject 覆盖度."""
    work = Path(workdir); work.mkdir(parents=True, exist_ok=True)
    db_dmnd = str(aa_db) + '.dmnd'
    if not Path(db_dmnd).is_file():
        run(['diamond', 'makedb', '--in', str(aa_db), '--db', db_dmnd], 'diamond makedb (蛋白库索引)')
    out_tsv = work / 'blastx.m8'
    # 输出含 subject 坐标(sstart/send)和全长(slen), 多个 target 便于评估
    run(['diamond', 'blastx', '--query', str(fasta), '--db', db_dmnd,
         '--out', str(out_tsv), '--evalue', '1e-3',
         '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'sstart', 'send', 'slen', 'evalue', 'bitscore',
         '--max-target-seqs', '5', '--max-hsps', '5', '--threads', str(threads), '--sensitive'],
        'diamond blastx')
    # 按 subject 蛋白分组合并 HSP 区间, 算该蛋白被覆盖比例; 取 pident 最高的 HSP 的 identity
    from collections import defaultdict
    per_q = defaultdict(lambda: defaultdict(list))  # q → target → [(sstart,send),...]
    pid_map = defaultdict(dict)  # q → target → max_pident
    with open(out_tsv) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 7:
                continue
            q, t = p[0], p[1]
            try:
                pid = float(p[2])
                ss, se = int(float(p[3])), int(float(p[4]))
                slen = int(float(p[5]))
                if ss > se: ss, se = se, ss
            except Exception:
                continue
            per_q[q][t].append((ss, se, slen))
            if t not in pid_map[q] or pid > pid_map[q][t]:
                pid_map[q][t] = pid
    # 对每个 query, 在它的 target 里取"subject覆盖占比最大"的代表
    best = {}
    for q, tars in per_q.items():
        best_t, best_cov, best_pid = None, -1, 0.0
        for t, hsp in tars.items():
            slen = hsp[0][2]
            ivs = sorted((a, b) for a, b, _ in hsp)
            merged = []
            for a, b in ivs:
                if merged and a <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            cov = sum(b - a for a, b in merged)
            cov_frac = cov / slen if slen else 0.0
            if cov_frac > best_cov:  # 选覆盖该蛋白比例最大的
                best_t, best_cov, best_pid = t, cov_frac, pid_map[q][t]
        if best_t:
            best[q] = (best_pid, best_t, round(100 * best_cov, 1))
    return best



def load_cdd_classified(path):
    """cdd_classified_taxid.tsv → {accession: category}. category∈{Virus,Plant,Animal,Fungi,Bacteria,Archaea,Eukaryota,Protist,Root,...}"""
    cat = {}
    if path and Path(path).is_file():
        with open(path) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                a = r.get('accession', '').strip()
                if a:
                    cat[a] = r.get('category', '').strip()
    return cat


def decide_category(cdd_cat, evidence, nt_pi, nt_qc, aa_pi, nt_sp, length, cdd_domains):
    """category: 基于 cdd 域分类(病毒/植物/动物/真菌/细菌) + blast + 长度.
    与 archive 口径对齐: Contamination类由 cdd 域类型主导, blast 决定病毒已知/新种/Distantly related."""
    # 该序列各 cdd 域的分类 (cdd_classified_taxid 的 category 列)
    dom_cats = set()
    for d in cdd_domains:
        c = cdd_cat.get(d, '')
        if c:
            dom_cats.add(c)
    virus_dom = 'Virus' in dom_cats
    plant_dom = 'Plant' in dom_cats
    anim_dom  = 'Animal' in dom_cats or 'Arthropoda' in dom_cats
    fung_dom  = 'Fungi' in dom_cats
    bact_dom  = 'Bacteria' in dom_cats or 'Archaea' in dom_cats
    euk_dom   = 'Eukaryota' in dom_cats or 'Protist' in dom_cats
    host_dom  = plant_dom or anim_dom or fung_dom or bact_dom or euk_dom
    host = classify_host(nt_sp)

    # 1. Short fragment / low coverage (no viral domain)
    if not virus_dom and length and int(length) < SHORT_LEN:
        if host_dom:
            return 'Host gene contamination'
        return 'Short fragment / low coverage'

    # 2. 病毒已知型 (blast nt 命中 + query coverage 门槛)
    #    已知种需 qcov>=85% (物种级), 新种候选需 qcov>=50% (跨种/属级)
    if nt_pi >= NT_KNOWN and nt_qc >= NT_QCOV_KNOWN:
        return 'Known plant virus'
    if nt_pi >= NT_RELATED and nt_qc >= NT_QCOV_RELATED:
        return '新种候选'

    # 3. 噬菌体 (blast 命中噬菌体, 或 cdd 域与噬菌体相关)
    if host == 'phage':
        return '噬菌体Contamination'

    # 4. 真菌/昆虫 (cdd 域类型 + blast 命中辅助)
    if fung_dom or host == 'fungus':
        return '真菌病毒/Contamination'
    if anim_dom or host == 'insect':
        return 'Insect virus/Contamination'

    # 5. Distantly related (蛋白级命中)
    if aa_pi >= AA_RELATED:
        if virus_dom:
            return 'Distantly-related virus candidate'
        return 'Distantly-related (undetermined)'

    # 6. Host gene contamination (有宿主域, no viral domain, blast 弱)
    if host_dom and not virus_dom:
        return 'Host gene contamination'

    # 7. 兜底
    if virus_dom:
        return 'Distantly-related virus candidate'
    return 'Undetermined'


def main():
    ap = argparse.ArgumentParser(description='生成 final_judgement_table.tsv')
    ap.add_argument('--fasta', required=True, help='HQ_plant_viruses.fasta (全序列, 含known参考)')
    ap.add_argument('--freq', default='', help='rescue_detection/frequency_table.tsv')
    ap.add_argument('--cdd', default='', help='rescue_validation/cdd_hits.tsv')
    ap.add_argument('--cdd-evidence', default='', help='rescue_validation/cdd_evidence_report.tsv')
    ap.add_argument('--sample', default='', help='物种名 (host_species)')
    ap.add_argument('--out', default='final_judgement_table.tsv')
    ap.add_argument('--nt-db', default=os.environ.get(
        'MMPV_NT_DB', str(Path.home() / 'database' / 'virus-db' / 'ncbi-virus' / 'ncbi_virus.nucl.fasta')),
        help='核酸参考库 (默认 $MMPV_NT_DB 或 ~/database/virus-db/ncbi-virus/ncbi_virus.nucl.fasta)')
    ap.add_argument('--aa-db', default=os.environ.get(
        'MMPV_AA_DB', str(Path.home() / 'database' / 'virus-db' / 'ncbi-virus' / 'ncbi_virus.prot.fasta')),
        help='蛋白参考库 (默认 $MMPV_AA_DB 或 ~/database/virus-db/ncbi-virus/ncbi_virus.prot.fasta)')
    ap.add_argument('--install-mode', action='store_true', help='只建 blast/mmseqs 索引, 不跑 query')
    ap.add_argument('--threads', type=int, default=32)
    args = ap.parse_args()

    import tempfile as _tmp
    workdir = Path(_tmp.mkdtemp(prefix='judgement_'))
    print(f"[gen_final_judgement] sample={args.sample or '?'} fasta={args.fasta}")
    print(f"[gen_final_judgement] 工作目录: {workdir}")

    # ── blast 索引 (若缺则建) ──
    nt_db, aa_db = Path(args.nt_db), Path(args.aa_db)
    if not nt_db.is_file() or not aa_db.is_file():
        sys.exit(f'blast库不存在: {nt_db} / {aa_db}')

    # blastn 索引 (.nin)
    nt_idx = Path(str(nt_db) + '.nin')
    if nt_idx.is_file():
        print('  ✓ blastn 索引已存在, 跳过 makeblastdb')
    else:
        print(f'  建 blastn 索引 ({nt_db.name}) ...')
        run(['makeblastdb', '-in', str(nt_db), '-dbtype', 'nucl',
             '-parse_seqids', '-title', nt_db.stem], 'makeblastdb nucl')

    # diamond blastx 蛋白索引 (target 库 <name>.dmnd)
    diamond_dmnd = str(aa_db) + '.dmnd'
    if Path(diamond_dmnd).is_file():
        print('  ✓ diamond blastx 索引已存在, 跳过')
    else:
        print(f'  建 diamond blastx 索引 ({aa_db.name}) ...')
        run(['diamond', 'makedb', '--in', str(aa_db), '--db', diamond_dmnd], 'diamond makedb')

    if args.install_mode:
        print('[install-mode] 索引已建好, 退出.')
        return

    # ── 读 fasta IDs ──
    ids, lens = load_fasta_ids(args.fasta)
    if not ids:
        sys.exit('fasta 无序列')

    # ── blastn query ──
    print(f'  blastn query ({len(ids)} 条 × {nt_db.stem}) ...')
    blast_res = blastn_query(args.fasta, nt_db, args.threads,
                              str(workdir / 'blastn.m8'))

    # ── diamond blastx query ──
    print(f'  diamond blastx query ({len(ids)} 条 × {aa_db.stem}) ...')
    aa_res = mmseqs_query(args.fasta, aa_db, args.threads, str(workdir / 'aa_query'))

    # ── 读辅助文件 ──
    freq = load_freq(args.freq)
    cdd_top = load_cdd_top(args.cdd)
    cdd_evid = load_cdd_evidence(args.cdd_evidence)
    cdd_cat = load_cdd_classified(DEF_CDD_CLASSIFIED)

    # ── 写 final_judgement_table.tsv ──
    out = Path(args.out)
    with open(out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['contig_id', 'host_species', 'length', 'prevalence_%',
                    'nt_pident', 'nt_qcovs', 'nt_species',
                    'aa_pident', 'aa_qcov', 'aa_species', 'cdd_top', 'category'])
        for cid in ids:
            if cid in blast_res:
                nt_pi, nt_qc, nt_sp = blast_res[cid][0], blast_res[cid][1], blast_res[cid][2]
            else:
                nt_pi, nt_qc, nt_sp = 0.0, 0.0, ''
            aa_pi, aa_sp, aa_qc = aa_res.get(cid, (0.0, '', 0.0)) if cid in aa_res else (0.0, '', 0.0)
            length = lens.get(cid, '')
            evidence = cdd_evid.get(cid, '')
            cat = decide_category(cdd_cat, evidence, nt_pi, nt_qc, aa_pi, nt_sp, length,
                                  [d for d in cdd_top.get(cid, '').split(';') if d])
            w.writerow([
                cid, args.sample, length,
                freq.get(cid, ''),
                f'{nt_pi:.1f}' if nt_pi else '',
                f'{nt_qc:.1f}' if nt_qc else '',
                nt_sp,
                f'{aa_pi:.1f}' if aa_pi else '',
                f'{aa_qc:.1f}' if aa_qc else '',
                aa_sp,
                cdd_top.get(cid, ''),
                cat,
            ])
    print(f'✓ final_judgement_table.tsv: {len(ids)} 条 → {out}')
    # 汇总
    if out.is_file():
        from collections import Counter
        cats = Counter()
        with open(out) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                cats[r['category']] += 1
        print('  category 分布:')
        for k, v in cats.most_common():
            print(f'    {k:<14} {v}')


if __name__ == '__main__':
    main()
