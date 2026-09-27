#!/usr/bin/env python3
"""hmm_ct3_evidence.py -- CT3 病毒 HMM 库第三路证据 (pyhmmer) + 模型名分级

背景与定位:
    9b (analysis_verify) 的证据链是「一道门 + 两次补票」:
      门    = CDD 结构域验证 (validate_rescue_cdd.py, curated 模型, 特异性高覆盖窄)
      补票1 = blast 身份注释 (known / new)
      补票2 = HMM 全长 profile 库 (覆盖面宽, 能碰到 CDD 没有的类群)
    本脚本承担补票2, 取代早期的 RVDB-prot 路径 (hmm_rvdb_evidence.py)。

    核酸 contig 先做六框翻译成 ORF, 再用 pyhmmer hmmscan 打 CT3 的五个库:
      RDRP_HMMs / DNA_rep_HMMs / Virion_HMMs / Useful_Annotation_HMMs / phrogs_for_ct
    每个库单独扫一遍, 命中按 contig 取并集 (不物理合并 HMM: CT3 是二进制 .h3m,
    RVDB 是已 press 的文本, 合并需转文本重压且混两套命名体系, 得不偿失)。

为什么不能「命中即病毒」:
    CT3 不是纯病毒库。phrogs_for_ct 是噬菌体专属; Useful_Annotation 里塞了大量
    cdd_cluster (COG/PHA 噬菌体结构蛋白) 与 baits_missed (噬菌体尾部/门户);
    DNA_rep 里有 circular_genetic_element / plasmid 这类非病毒环状元件。
    因此必须按「模型名」把命中分成证据等级, 而不是简单计数。

模型命名体系 (实测):
    <cluster_id>/<SOURCE>:<ACC>-<description>
      cluster_id 三种形态:
        ① 分类单元名   Alphaflexiviridae.1 / Geminiviridae_cluster / Varsani_CRESSviridae_cluster
        ② 支系簇名     negarna_cluster.1 / reo_cluster.3 / some_caudo_cluster / phrog_505
        ③ 家族簇名     f.0186.1 / cdd_cluster.379 / baits_missed.110 / VOG4586
      description 是功能注释: RNA-dependent-RNA-polymerase / Terminase / Coat-protein ...

证据等级 (六桶, 高 -> 低):
    VIRAL_FAMILY   cluster_id 带 ICTV 分类单元后缀 (viridae/virinae/virales/virus/viroid/satellit)
    VIRAL_CLUSTER  cluster_id 是已知病毒支系簇 (negarna/reo/tombus/gemini/unclassified_RNA/varsani_ ...)
    VIRAL_FUNCTION description 是病毒专属功能 (RdRp/replicase/coat/capsid/movement/polyprotein ...)
    AMBIGUOUS      VOG 这类跨噬菌体/真核病毒无法判定的模型 (既不当作病毒证据也不当作噬菌体)
    PHAGE          噬菌体 (phrogs 库 / phrog_ seeker caudo phage podo sipho myo micro ino 等簇)
    OTHER          其余 (plasmid / 宿主样功能簇 / 未知)

    三类 VIRAL_* 合起来 = 非噬菌体病毒证据。下游 integrate_rescue_evidence.py 用
    --hmm-min-specific 在这个计数上设阈值; 用 --hmm-require-family 收紧到只认 VIRAL_FAMILY;
    需要放宽时用 --hmm-count-ambiguous 把 AMBIGUOUS 也算进来。

    RDRP_HMMs 整库就是病毒 RdRp 库, 库内非噬菌体模型一律归 VIRAL_FUNCTION,
    不受 description 写法影响。

与 CDD 的分工:
    CDD 是 curated 结构域 (窄而准); CT3 是数据库自建全长 profile (宽而杂)。互补不替代。

用法:
    python hmm_ct3_evidence.py --fasta HQ_plant_viruses.fasta \
        --ct3-dir ~/database/virus-db/ct3_DBs/hmmscan_DBs/v3.1.1 \
        --outdir 09b_Analysis_Verify/virus_validation --threads 60

输出 (outdir):
    hmm_ct3_hits.tsv      每 contig 一行汇总 (等级计数 / 最佳病毒命中 / 库来源)
    hmm_ct3_hits_raw.tsv  每个显著性命中一行 (ORF 级, 含核酸坐标与证据等级)
"""
import argparse
import csv
import re
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

# site-packages 里有个 pyproject.toml, 会让 Biopython 误报 source-tree 警告, 噪音屏蔽
warnings.filterwarnings('ignore', message='.*importing Biopython from inside the source tree.*')

from Bio.Seq import Seq

# pyhmmer 只在真正扫描时必需; 模块级分级规则 (classify_ct3_model) 可离线单测
try:
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextSequence
    from pyhmmer.plan7 import HMMFile
    from pyhmmer.hmmer import hmmscan
    _PYHMMER_ERR = None
except ImportError as _e:
    pyhmmer = None
    _PYHMMER_ERR = _e


# ── CT3 默认库与分级规则 ────────────────────────────────────────────────
DEFAULT_CT3_DIR = '~/database/virus-db/ct3_DBs/hmmscan_DBs/v3.1.1'
DEFAULT_DBS = ['RDRP_HMMs.h3m', 'DNA_rep_HMMs.h3m', 'Virion_HMMs.h3m',
               'Useful_Annotation_HMMs.h3m', 'phrogs_for_ct.h3m']

# 整库即为噬菌体, 库内模型一律归 PHAGE
PHAGE_DBS = {'phrogs_for_ct.h3m'}

# 整库即为病毒功能库: 库内非噬菌体模型一律归 VIRAL_FUNCTION
# (RDRP 库收录的就是病毒 RdRp, description 写法五花八门, 不靠关键词碰运气)
VIRAL_FUNCTION_DBS = {'RDRP_HMMs.h3m'}

# 证据等级 (高 -> 低); 用于挑「本 contig 最强的一类」
TIER_ORDER = ['VIRAL_FAMILY', 'VIRAL_CLUSTER', 'VIRAL_FUNCTION',
              'AMBIGUOUS', 'PHAGE', 'OTHER']
VIRAL_TIERS = ('VIRAL_FAMILY', 'VIRAL_CLUSTER', 'VIRAL_FUNCTION')

# ① 噬菌体支系簇 token (cluster_id 子串, 小写匹配)
PHAGE_TOKENS = (
    'phrog_', 'seeker_phage', 'seeker',
    'caudo',                      # some_caudo_cluster / unclassified_caudo_cluster
    'phage',                      # unclassified_dsDNA_phage_cluster
    'crass',                      # crassvirales / crass_cluster
    'ackermann', 'herell',
    'podo_cluster', 'sipho_cluster', 'myo_cluster', 'tecti_cluster',
    'levi_cluster', 'corti_cluster', 'halo_cluster', 'fusello_cluster',
    'microviridae_cluster', 'micro_cluster', 'ino_cluster',
    'some_myo', 'some_sipho', 'some_podo', 'some_tecti', 'some_levi',
)

# ② ICTV 分类单元后缀 (cluster_id)
VIRAL_TAXON_RE = re.compile(
    r'(viridae|virinae|virales|viricetes|viricota|viroid|satellit)', re.I)

# ③ 已知病毒支系簇 token (cluster_id 子串, 小写匹配)
VIRAL_CLUSTER_TOKENS = (
    # RNA 病毒
    'negarna', 'narna', 'tombus', 'sobemo', 'luteo', 'toti', 'partiti',
    'clostero', 'poty', 'sequi', 'flavi', 'toga', 'corona', 'arteri', 'ifla',
    'bunya', 'phlebo', 'orthomyxo', 'rhabdo', 'filo', 'borna', 'beny',
    'reo_cluster', 'astro_cluster', 'calici', 'picorna', 'virga', 'bromo',
    'secovi', 'tymo', 'cucomo', 'polero', 'carmo', 'furo', 'noda', 'lavida',
    'alphaflexi', 'betaflexi', 'gammalexi', 'deltaflexi',
    # 未分类 / 环境病毒簇
    'unclassified_rna', 'unclassified_dsrna', 'unclassified_ssrna',
    'unclassified_cress', 'unclassified_ssdna', 'varsani_', 'cress',
    'all_no_hit_virus', 'shth',
    # DNA 病毒 (含大型 dsDNA)
    'gemin', 'nano', 'circo', 'baculo', 'nudi', 'irido', 'asco', 'pox',
    'mimi', 'phycodna', 'marseille', 'pitho', 'nima', 'herpes', 'adeno',
    'adinto', 'adoma', 'papilloma', 'polyoma', 'parvo', 'anello', 'ligamen',
    'genomoviridae', 'smacoviridae', 'satellites', 'hydrosa', 'endo',
)

# ④ 病毒专属功能注释 (description 子串, 正则)
# ④ 病毒专属功能注释 (description 子串, 正则)
VIRAL_FUNC_RE = re.compile(
    r'RNA[- ]dependent[- ]RNA[- ]polymerase|RNA[- ]directed[- ]RNA[- ]polymerase|'
    r'RNA[- ]dependent[- ]DNA[- ]polymerase|'
    r'replicase|replication[- ]protein|replication[- ]associated[- ]protein|'
    r'replicative[- ]helicase|large[- ]T[- ]antigen|\bT[- ]antigen|\bE1\b|\bRep[- ]protein|'
    r'\bRdRp\b|'
    r'\bcoat\b|\bcapsid\b|nucleocapsid|\bnucleoprotein\b|hexon|penton|adenain|'
    r'movement[- ]protein|polyprotein|genome[- ]linked|VPg|'
    r'silencing[- ]suppressor|triple[- ]gene[- ]block|helper[- ]component|'
    r'virion[- ]protein|readthrough[- ]protein',
    re.I)

# ⑤ 噬菌体专属功能注释 (description 子串, 正则) — 只收噬菌体排他性词条
PHAGE_FUNC_RE = re.compile(
    r'terminase|portal[- ]protein|baseplate|'
    r'tail[- ](?:tube|fiber|fibre|sheath|assembly|adaptor|length|tape|spike|protein|tip)|'
    r'tape[- ]measure|dna[- ]packaging|packaging[- ](?:atpase|protein|terminase)|'
    r'tegument|head[- ](?:to[- ]tail|protein|assembly|morphogenesis)|'
    r'holin|endolysin|lysozyme|muramidase|virion[- ]structural|'
    r'prohead|scaffold[- ]protein|neck[- ]protein|'
    r'phage[- ](?:tail|protein|structural)',
    re.I)


def split_model(model_name):
    """<cluster_id>/<SOURCE>:<ACC>-<description> -> (cluster_id, desc)

    注意: 不能拿第一个 '-' 去切说明文字。有些条目写成
    `PFAM:PHA00497;pol;RNA-dependent-RNA-polymerase`, 第一个 '-' 落在
    「RNA-dependent」里面, 会把描述切成 `dependent-RNA-polymerase` 而丢掉主语,
    导致 RdRp 被误归 OTHER。统一保留 '/' 之后的整段作为可搜索文本。
    """
    head, _, rest = model_name.partition('/')
    return head, rest


def classify_ct3_model(db_name, model_name):
    """把一条 CT3 命中按模型名分到证据等级。返回 (tier, cluster_id, desc)。

    判定顺序即优先级: 噬菌体 -> 分类单元 -> 病毒支系簇 -> 病毒功能注释 -> 其他。
    噬菌体必须最先排除, 否则 microviridae_cluster 会被 viridae 规则误判为病毒家族。
    """
    head, desc = split_model(model_name)
    hl, dl = head.lower(), desc.lower()

    # ── 1. 噬菌体 ──
    if db_name in PHAGE_DBS:
        return 'PHAGE', head, desc
    if hl.startswith('vog'):
        # VOG (Viral Orthologous Groups) 跨噬菌体/真核病毒, 无法判定, 单独归 AMBIGUOUS
        return 'AMBIGUOUS', head, desc
    if any(t in hl for t in PHAGE_TOKENS):
        return 'PHAGE', head, desc
    if PHAGE_FUNC_RE.search(desc):
        return 'PHAGE', head, desc

    # ── 2. 命名的病毒分类单元 ──
    if VIRAL_TAXON_RE.search(head):
        return 'VIRAL_FAMILY', head, desc

    # ── 3. 已知病毒支系簇 ──
    if any(t in hl for t in VIRAL_CLUSTER_TOKENS):
        return 'VIRAL_CLUSTER', head, desc

    # ── 4. 病毒专属功能注释 ──
    if VIRAL_FUNC_RE.search(desc):
        return 'VIRAL_FUNCTION', head, desc

    # ── 5. 整库即病毒功能库的兜底 (RDRP) ──
    if db_name in VIRAL_FUNCTION_DBS:
        return 'VIRAL_FUNCTION', head, desc

    return 'OTHER', head, desc


# ---------- 六框 ORF ----------
def six_frame_orfs(contig_id, seq, min_aa=30):
    """返回 [(orf_id, prot, nt_start, nt_end, strand, frame)] (1-based frame)"""
    seq = seq.upper()
    rc = str(Seq(seq).reverse_complement())
    orfs = []
    for strand, s in (('+', seq), ('-', rc)):
        for f in range(3):
            sub = s[f:]
            ncod = len(sub) // 3
            if ncod == 0:
                continue
            aa = str(Seq(sub[:ncod * 3]).translate())
            start = 0
            for i in range(len(aa) + 1):
                if i == len(aa) or aa[i] == '*':
                    frag = aa[start:i]
                    if len(frag) >= min_aa:
                        orfs.append((
                            '%s_ORF%d_%s%d' % (contig_id, len(orfs) + 1, strand, f + 1),
                            frag,
                            f + start * 3,
                            f + i * 3,
                            strand,
                            f + 1,
                        ))
                    start = i + 1
    return orfs


def _s(x):
    """hmmscan 返回的 name 可能是 str 也可能是 bytes, 统一成 str"""
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def read_fasta(path):
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)


def resolve_dbs(ct3_dir, dbs):
    """把裸文件名或路径解析成存在的 .h3m 绝对路径列表"""
    out = []
    for d in dbs:
        p = Path(d)
        if not p.is_file():
            p = Path(ct3_dir).expanduser() / d
        if not p.is_file():
            sys.exit('CT3 库不存在: %s (在 %s 下也没找到)' % (d, ct3_dir))
        out.append(p)
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True, help='核酸 contig FASTA')
    ap.add_argument('--ct3-dir', default=DEFAULT_CT3_DIR,
                    help='CT3 hmmscan 库目录 (默认 %s)' % DEFAULT_CT3_DIR)
    ap.add_argument('--dbs', default=None,
                    help='逗号分隔的库文件名/路径 (默认五库全部: %s)' % ','.join(DEFAULT_DBS))
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--threads', type=int, default=60)
    ap.add_argument('--min-aa', type=int, default=30, help='ORF 最小氨基酸长度')
    ap.add_argument('--evalue', type=float, default=1e-3, help='序列级 E-value 阈值 (扫描)')
    ap.add_argument('--dom-evalue', type=float, default=1e-3, help='域级 E-value 阈值 (扫描)')
    ap.add_argument('--strict-evalue', type=float, default=1e-5,
                    help='汇总时的强命中阈值 (等级计数与门控用)')
    ap.add_argument('--top-models', type=int, default=5, help='每 contig 保留前 N 个模型名')
    args = ap.parse_args()

    if pyhmmer is None:
        sys.exit('pyhmmer 不可用: %s (pip install pyhmmer)' % _PYHMMER_ERR)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dbs = resolve_dbs(args.ct3_dir, (args.dbs.split(',') if args.dbs else DEFAULT_DBS))
    print('[ct3-hmm] 库 (%d): %s' % (len(dbs), ', '.join(p.name for p in dbs)))

    # ── 六框翻译 (一次, 五个库共用) ──
    contigs = list(read_fasta(args.fasta))
    abc = Alphabet.amino()
    queries = []
    meta = {}
    contig_of = {}
    for cid, seq in contigs:
        for oid, prot, ns, ne, strand, frame in six_frame_orfs(cid, seq, args.min_aa):
            queries.append(TextSequence(name=oid.encode(), sequence=prot).digitize(abc))
            meta[oid] = (cid, ns, ne, strand, frame, len(prot))
            contig_of[oid] = cid
    print('[ct3-hmm] contigs=%d, ORFs(>=%daa)=%d' % (len(contigs), args.min_aa, len(queries)))
    if not queries:
        print('[ct3-hmm] 无 ORF, 退出')
        return

    # ── 逐库 pyhmmer 搜索, 命中并集 ──
    # raw 列: 0 contig_id 1 orf_id 2 strand 3 frame 4 nt_start 5 nt_end 6 aa_len
    #         7 db 8 model_name 9 cluster_id 10 tier 11 desc 12 evalue 13 score
    raw_rows = []
    hit_by_contig = defaultdict(list)
    for db in dbs:
        n_db = 0
        with HMMFile(str(db)) as hf:
            for hits in hmmscan(queries, hf, cpus=args.threads,
                                E=args.evalue, domE=args.dom_evalue):
                qname = _s(hits.query.name)
                cid = contig_of.get(qname)
                if cid is None:
                    continue
                for hit in hits:
                    model = _s(hit.name)
                    tier, cluster, desc = classify_ct3_model(db.name, model)
                    row = (
                        cid, qname, meta[qname][3], meta[qname][4],
                        meta[qname][1], meta[qname][2], meta[qname][5],
                        db.name, model, cluster, tier, desc,
                        '%.3g' % hit.evalue, '%.1f' % hit.score,
                    )
                    raw_rows.append(row)
                    hit_by_contig[cid].append(row)
                    n_db += 1
        print('[ct3-hmm]   %-28s 显著命中 %d' % (db.name, n_db))

    print('[ct3-hmm] 合计显著命中 %d 条, 覆盖 contig %d 个' % (len(raw_rows), len(hit_by_contig)))

    # ── 输出 raw ──
    raw_head = ['contig_id', 'orf_id', 'strand', 'frame', 'nt_start', 'nt_end',
                'aa_len', 'db', 'model_name', 'cluster_id', 'tier', 'desc',
                'evalue', 'score']
    with open(outdir / 'hmm_ct3_hits_raw.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(raw_head)
        for r in sorted(raw_rows, key=lambda x: (x[0], float(x[12]))):
            w.writerow(r)

    # ── 汇总 per contig ──
    sum_head = ['contig_id', 'n_ct3_hit', 'n_ct3_model', 'n_dbs_hit', 'dbs_hit',
                'best_model', 'best_db', 'best_tier', 'best_evalue', 'best_score',
                'ct3_hmm_evidence', 'ct3_viral_evidence',
                'ct3_n_viral', 'ct3_n_family', 'ct3_n_cluster', 'ct3_n_function',
                'ct3_n_ambiguous', 'ct3_n_phage', 'ct3_n_other', 'ct3_tier_counts',
                'ct3_best_viral_model', 'ct3_best_viral_cluster', 'ct3_best_viral_tier',
                'ct3_best_viral_db', 'ct3_best_viral_desc',
                'ct3_best_viral_evalue', 'ct3_best_viral_score',
                'ct3_best_family_model', 'ct3_best_family_cluster',
                'ct3_best_family_evalue', 'ct3_best_family_score',
                'ct3_top_models']
    with open(outdir / 'hmm_ct3_hits.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(sum_head)
        for cid, _ in contigs:
            rows = hit_by_contig.get(cid)
            if not rows:
                w.writerow([cid, 0, 0, 0, '', '', '', '', '', '',
                            'CT3_NO_HIT', 'CT3_NO_VIRAL',
                            0, 0, 0, 0, 0, 0, 0, '',
                            '', '', '', '', '', '', '',
                            '', '', '', '', ''])
                continue

            rows_sorted = sorted(rows, key=lambda x: float(x[12]))
            best = rows_sorted[0]

            # 强命中 (E <= strict) 的等级计数
            strict = [r for r in rows_sorted if float(r[12]) <= args.strict_evalue]
            cc = Counter(r[10] for r in strict)
            viral = [r for r in strict if r[10] in VIRAL_TIERS]
            family = [r for r in strict if r[10] == 'VIRAL_FAMILY']

            def pick(bucket, fields):
                """取该桶内 E-value 最小的命中, 按 fields 抽出对应列"""
                if not bucket:
                    return [''] * len(fields)
                b = min(bucket, key=lambda x: float(x[12]))
                return [b[i] for i in fields]

            # viral 桶: model(8) cluster(9) tier(10) db(7) desc(11) evalue(12) score(13)
            v = pick(viral, [8, 9, 10, 7, 11, 12, 13])
            # family 桶: model(8) cluster(9) evalue(12) score(13)
            fm = pick(family, [8, 9, 12, 13])
            counts = ';'.join('%s:%d' % (k, cc[k]) for k in TIER_ORDER if cc.get(k))
            top_models = []
            for r in rows_sorted:
                if r[8] not in top_models:
                    top_models.append(r[8])
                if len(top_models) >= args.top_models:
                    break
            w.writerow([cid, len(rows), len(set(r[8] for r in rows)),
                        len(set(r[7] for r in rows)), ';'.join(sorted(set(r[7] for r in rows))),
                        best[8], best[7], best[10], best[12], best[13],
                        'CT3_HIT',
                        'CT3_VIRAL' if viral else 'CT3_NO_VIRAL',
                        len(viral), len(family), cc.get('VIRAL_CLUSTER', 0),
                        cc.get('VIRAL_FUNCTION', 0), cc.get('AMBIGUOUS', 0),
                        cc.get('PHAGE', 0), cc.get('OTHER', 0),
                        counts] + v + fm + [';'.join(top_models)])
    n_hit_c = sum(1 for cid, _ in contigs if hit_by_contig.get(cid))
    n_viral_c = sum(1 for cid, _ in contigs
                    if any(r[10] in VIRAL_TIERS and float(r[12]) <= args.strict_evalue
                           for r in hit_by_contig.get(cid, [])))
    print('[ct3-hmm] 完成: %d/%d contig 有命中, %d/%d 有非噬菌体病毒证据 '
          '(E<=%.0e) → %s' % (n_hit_c, len(contigs), n_viral_c, len(contigs),
                              args.strict_evalue, outdir))

    # 分级统计, 供审计 (判定规则是否合理看这里)
    tc = Counter()
    for r in raw_rows:
        if float(r[12]) <= args.strict_evalue:
            tc[(r[10], r[7])] += 1
    print('[ct3-hmm] 强命中等级 x 库 分布:')
    for tier in TIER_ORDER:
        tot = sum(v for (t, _), v in tc.items() if t == tier)
        if not tot:
            continue
        det = ', '.join('%s=%d' % (db, v) for (t, db), v in sorted(tc.items())
                        if t == tier)
        print('    %-15s %5d  (%s)' % (tier, tot, det))


if __name__ == '__main__':
    main()
