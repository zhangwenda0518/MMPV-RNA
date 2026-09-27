#!/usr/bin/env python3
"""
integrate_rescue_evidence.py — rescue 病毒五层证据链整合 + 30/30/40 加权评分
（virome_discovery_pipeline 组件，接在 validate_rescue_cdd.py 之后）

证据链:
  ① 病毒属性 = category (七工具交叉, final_judgement_table)
  ② 质量评估 = checkv_completeness (prevalence_full_table)
  ③ 域验证   = CDD T1/T2/T3 (validate_rescue_cdd 输出)          ← 唯一的「门」
  ④ 身份注释 = nt_species (blastn) + aa_species (blastx)         ← 补票 1
  ⑤ 新病毒判定 = identity + 30/30/40 加权评分
  ⑥ 全长 profile 补充证据 = CT3 病毒 HMM 库 (pyhmmer, hmm_ct3_evidence.py)
                ← 补票 2。默认参与裁决: 无 CDD 病毒域但 CT3 有非噬菌体病毒命中时
                  把 DROP 抬到 REVIEW。用 --no-hmm-rescue 可退回「只注释不裁决」。

  一句话: 一道门 (CDD 结构域) 加两次补票 (blast 身份 / HMM 全长 profile)。

加权评分 (对齐 ViralQuest 30/30/40):
  score = 0.3*blastn + 0.3*blastx + 0.4*domain
  domain = max(cdd, hmm)   ← 域证据层: CDD 与 CT3-HMM 取 max, 同为 profile 证据
                              不叠加 (防双重计分), CDD 盲区由 HMM 补上
  最终判定: CDD 病毒域(T1)为第一判据, blast 为辅助身份
  HMM 进分默认开 (--no-hmm-in-score 回退: 域证据层只看 CDD)

用法:
  python integrate_rescue_evidence.py \
      --judgement final_judgement_table.tsv \
      --cdd-report validate_out/cdd_evidence_report.tsv \
      --prevalence prevalence_full_table.tsv \
      --hmm-report validate_out/hmm_ct3_hits.tsv \
      --outdir .
"""
import argparse, csv, os
from collections import Counter

# CDD 证据等级 → 分数 (30/30/40 中的 40% 权重层)
CDD_SCORE = {
    'PASS_CORE': 100.0,
    'PASS_CORE_FAM': 90.0,
    'PASS_VIRAL': 80.0,
    'AMBIGUOUS': 40.0,
    'REVIEW': 20.0,
    'NON_VIRAL': 0.0,
}

# CT3-HMM 最佳病毒命中等级 → base 分 (与 CDD_SCORE 同尺度; 进域证据层取 max)
# 定位: 全长 profile 命中特异性低于 curated CDD 域, 各档比 CDD 同义档低一档:
#   VIRAL_FAMILY ≈ CDD PASS_VIRAL (75<80), VIRAL_FUNCTION ≈ AMBIGUOUS/REVIEW 之间
HMM_TIER_SCORE = {
    'VIRAL_FAMILY': 75.0,
    'VIRAL_CLUSTER': 65.0,
    'VIRAL_FUNCTION': 50.0,
    'AMBIGUOUS': 20.0,
}
HMM_COUNT_BONUS = 5.0   # 每多 1 条口径内病毒命中 +5
HMM_BONUS_CAP = 25.0    # 数量加成封顶 +5 条 (防长 contig 多 ORF 轰分)

# 新病毒判定阈值
NT_KNOWN = 85.0    # >=85% 已知种
NT_RELATED = 40.0  # 40-85% Closely related/新种候选
AA_RELATED = 40.0  # 蛋白级Distantly related阈值


def safe_float(x, default=0.0):
    try:
        return float(x) if x not in ('', None) else default
    except (ValueError, TypeError):
        return default


def safe_int(x, default=0):
    try:
        return int(float(x)) if x not in ('', None) else default
    except (ValueError, TypeError):
        return default


def load_judgement(path):
    rows = {}
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            rows[r['contig_id']] = r
    return rows


def load_cdd_report(path):
    evid = {}
    if not os.path.exists(path):
        return evid
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            evid[r['seq']] = r
    return evid


def load_prevalence(path):
    cv = {}
    if not path or not os.path.exists(path):
        return cv
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            cid = r.get('contig_id', '')
            if cid:
                cv[cid] = {
                    'completeness': r.get('checkv_completeness', ''),
                    'confidence': r.get('checkv_confidence', ''),
                    'prevalence': r.get('prevalence_%', ''),
                    'n_samples': r.get('n_samples', ''),
                    'avg_ANI': r.get('avg_ANI', ''),
                    'genus_avg_len': r.get('genus_avg_len', ''),
                    'genus': r.get('Genus', ''),
                    'family': r.get('Family', ''),
                }
    return cv


def load_fasta(path):
    """读 fasta → {seq_id: sequence}"""
    seqs = {}
    if not path or not os.path.exists(path):
        return seqs
    from Bio import SeqIO
    for rec in SeqIO.parse(path, 'fasta'):
        seqs[rec.id] = str(rec.seq)
    return seqs


def load_taxonomy(path):
    """读 consensus_tax (final_integrated_classification.tsv) → contig_id → {genus, family, species, primary_tool}"""
    tax = {}
    if not path or not os.path.exists(path):
        return tax
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            cid = r.get('contig_id', '').strip().strip('"')
            if not cid:
                continue
            tax[cid] = {
                'genus': (r.get('Genus', '') or '').strip(),
                'family': (r.get('Family', '') or '').strip(),
                'species': (r.get('Species', '') or '').strip(),
                'primary_tool': (r.get('primary_tool', '') or '').strip(),
            }
    return tax


def load_hmm_report(path):
    """读 hmm_ct3_hits.tsv (CT3-HMM 全长 profile 探针汇总) → contig_id → 整行

    兼容旧 RVDB 产物 (hmm_rvdb_hits.tsv, rvdb_* 列名): 取值时两套列名都试。
    """
    hmm = {}
    if not path or not os.path.exists(path):
        return hmm
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            cid = r.get('contig_id', '')
            if cid:
                hmm[cid] = r
    return hmm


def hm_get(row, ct3_key, old_key, default=''):
    """先取 ct3_* 列, 没有就退回旧 rvdb_* 列"""
    v = row.get(ct3_key)
    if v in (None, ''):
        v = row.get(old_key)
    return default if v in (None, '') else v


def hmm_evidence_score(best_tier, n_specific):
    """CT3-HMM → 0-100 分: 最佳病毒命中等级 base + 饱和数量加成.

    只有 VIRAL_* 等级给分 (AMBIGUOUS/PHAGE/OTHER 不给); n_specific 用裁决口径
    计数 (随 --hmm-require-family / --hmm-count-ambiguous 联动), 无有效命中得 0.
    """
    base = HMM_TIER_SCORE.get((best_tier or '').strip().upper(), 0.0)
    if base <= 0.0 or n_specific <= 0:
        return 0.0
    bonus = min(HMM_COUNT_BONUS * max(n_specific - 1, 0), HMM_BONUS_CAP)
    return min(base + bonus, 100.0)


def genus_of_species(species_name):
    """从物种名提取属名 (首个词, 去掉 sp./virus 等)"""
    if not species_name:
        return ''
    # 取第一个词 (属名)
    first = species_name.strip().split()[0]
    return first


def parse_length_from_id(contig_id):
    """从 contig_id 解析序列长度 (格式 ..._length_XXX_...)"""
    import re
    m = re.search(r'length_(\d+)', contig_id)
    return m.group(1) if m else ''


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--judgement', required=True, help='final_judgement_table.tsv (七工具分类 + blast)')
    ap.add_argument('--cdd-report', required=True, help='validate_rescue_cdd 输出的 cdd_evidence_report.tsv')
    ap.add_argument('--prevalence', default=None, help='prevalence_full_table.tsv (checkv 质量, 可选)')
    ap.add_argument('--taxonomy', default=None, help='consensus_tax (final_integrated_classification.tsv, 可选)')
    ap.add_argument('--fasta', default=None, help='输入序列 fasta (用于输出 keep/review/drop 的 fasta, 可选)')
    ap.add_argument('--genus-lens', default=None, help='属平均长度表 genus_lens (属→平均基因组长度, 用于core domain validation的长度门槛)')
    ap.add_argument('--outdir', default='.', help='输出目录')
    ap.add_argument('--hmm-report', default=None,
                    help='hmm_ct3_hits.tsv (CT3-HMM 全长 profile 探针汇总; 可选, 缺则第三路列留空)')
    ap.add_argument('--hmm-rescue', dest='hmm_rescue', action='store_true', default=True,
                    help='第三路探针 rescue (默认开): 无 CDD 病毒域但 CT3-HMM 有病毒命中时 DROP→REVIEW')
    ap.add_argument('--no-hmm-rescue', dest='hmm_rescue', action='store_false',
                    help='关闭第三路探针 rescue (退回只注释不裁决)')
    ap.add_argument('--hmm-promote-keep', dest='hmm_promote_keep', action='store_true', default=True,
                    help='第三路探针提升 (默认开): CT3-HMM 有非噬菌体病毒命中时 DROP/REVIEW 直接提升为 KEEP')
    ap.add_argument('--no-hmm-promote-keep', dest='hmm_promote_keep', action='store_false',
                    help='关闭 CT3-HMM 升 KEEP (退回 DROP→REVIEW 的保守口径)')
    ap.add_argument('--hmm-promote-min', type=int, default=1,
                    help='触发升 KEEP 所需的最少非噬菌体病毒命中数 (默认 1)')
    ap.add_argument('--hmm-min-specific', type=int, default=1,
                    help='触发 rescue 所需的最少病毒命中数 (默认 1), 口径见 --hmm-require-family')
    ap.add_argument('--hmm-require-family', action='store_true',
                    help='收紧口径: 只认命中 ICTV 分类单元的模型 (ct3_n_family), 而非全部非噬菌体病毒命中')
    ap.add_argument('--hmm-count-ambiguous', action='store_true',
                    help='放宽口径: 把 VOG 等无法判定的 AMBIGUOUS 命中也算进病毒证据')
    ap.add_argument('--hmm-in-score', dest='hmm_in_score', action='store_true', default=True,
                    help='HMM 进分 (默认开): 域证据层取 max(CDD, CT3-HMM), CDD 盲区由 HMM 补; '
                         '30/30/40 权重不变, 同为 profile 证据取 max 防双重计分')
    ap.add_argument('--no-hmm-in-score', dest='hmm_in_score', action='store_false',
                    help='关闭 HMM 进分 (域证据层只看 CDD, 回退旧口径; hmm_score 列仍输出供审计)')
    args = ap.parse_args()

    # 默认 genus_lens 路径 (优先 v2 表: total 已按节段属取最短完整段组)
    if not args.genus_lens:
        import os as _os
        for cand in [_os.path.expanduser("~/database/virus-db/db/genus_lens_v2.tsv"),
                     _os.path.expanduser("~/database/virus-db/db/genus_lens"),
                     _os.path.expanduser("~/MMPV-RNA/database/genus_lens_v2.tsv"),
                     _os.path.expanduser("~/MMPV-RNA/database/genus_lens")]:
            if _os.path.isfile(cand):
                args.genus_lens = cand
                break
    genus_len_map = {}
    if args.genus_lens and os.path.isfile(args.genus_lens):
        with open(args.genus_lens) as _gl:
            for _l in _gl:
                _l = _l.strip()
                if not _l or _l.startswith('genus') or '\t' not in _l:
                    continue
                _parts = _l.split('\t')      # v2 表多列, 只取前两列
                if len(_parts) < 2:
                    continue
                _g = _parts[0].replace('g__', '').strip()
                try:
                    genus_len_map[_g] = float(_parts[1])
                except Exception:
                    pass
        print(f'  属平均长度表: {args.genus_lens} ({len(genus_len_map)} 个属)')

    cdd_evidence = load_cdd_report(args.cdd_report)
    jt_rows = load_judgement(args.judgement)
    checkv = load_prevalence(args.prevalence)
    tax = load_taxonomy(args.taxonomy)
    seqs = load_fasta(args.fasta)
    hmm = load_hmm_report(args.hmm_report)
    print(f'CDD 判定 {len(cdd_evidence)} | 分类 {len(jt_rows)} | checkv {len(checkv)} | taxonomy {len(tax)} | fasta {len(seqs)} | ct3-hmm {len(hmm)}')

    # 第三路探针的裁决口径 (计数器在扫描阶段已按等级拆好, 这里只决定怎么算)
    gate_field = 'ct3_n_family' if args.hmm_require_family else 'ct3_n_viral'
    gate_field_old = 'rvdb_n_specific'
    gate_desc = ('CT3-HMM 命中 ICTV 分类单元' if args.hmm_require_family
                 else 'CT3-HMM 非噬菌体病毒命中')
    if args.hmm_count_ambiguous:
        gate_desc += ' + AMBIGUOUS'

    rows = []
    for cid, jr in jt_rows.items():
        cdd = cdd_evidence.get(cid, {})
        cdd_evid = cdd.get('evidence', 'NON_VIRAL')
        cdd_score = CDD_SCORE.get(cdd_evid, 0.0)

        nt_pident = safe_float(jr.get('nt_pident', ''))
        nt_qcovs = safe_float(jr.get('nt_qcovs', ''))
        nt_species = jr.get('nt_species', '')
        aa_pident = safe_float(jr.get('aa_pident', ''))
        aa_qcov = safe_float(jr.get('aa_qcov', ''))
        aa_species = jr.get('aa_species', '')
        cdd_top = jr.get('cdd_top', '')
        category = jr.get('category', '')
        length = jr.get('length', '')
        # length 列若空则从 contig_id 解析
        if not length:
            length = parse_length_from_id(cid)

        cv = checkv.get(cid, {})
        completeness = safe_float(cv.get('completeness', ''))
        confidence = cv.get('confidence', '')
        genus_avg_len = cv.get('genus_avg_len', '')
        n_samples = cv.get('n_samples', '')
        prevalence = cv.get('prevalence', '')

        # taxonomy 分类 (consensus_tax)
        tx = tax.get(cid, {})
        tax_genus = tx.get('genus', '')
        tax_family = tx.get('family', '')
        tax_species = tx.get('species', '')
        # blast 命中物种的属名
        blast_genus = genus_of_species(nt_species)

        # ── 第三路探针: CT3 病毒 HMM 库全长 profile 证据 (--hmm-report 缺则全空) ──
        hm = hmm.get(cid, {})
        hm_evidence = hm_get(hm, 'ct3_hmm_evidence', 'rvdb_hmm_evidence')
        hm_viral = hm_get(hm, 'ct3_viral_evidence', 'rvdb_viral_evidence')
        hm_n_viral = safe_int(hm_get(hm, 'ct3_n_viral', 'rvdb_n_viral', 0))
        hm_n_family = safe_int(hm_get(hm, 'ct3_n_family', 'rvdb_n_specific', 0))
        hm_n_ambiguous = safe_int(hm_get(hm, 'ct3_n_ambiguous', '', 0))
        hm_best_model = hm_get(hm, 'ct3_best_viral_model', 'rvdb_best_viral_fam')
        hm_best_cluster = hm_get(hm, 'ct3_best_viral_cluster', 'rvdb_best_viral_lca')
        hm_best_tier = hm_get(hm, 'ct3_best_viral_tier', 'rvdb_best_viral_class')
        hm_best_db = hm_get(hm, 'ct3_best_viral_db', '')
        hm_best_desc = hm_get(hm, 'ct3_best_viral_desc', '')
        hm_best_ev = hm_get(hm, 'ct3_best_viral_evalue', 'rvdb_best_viral_evalue')
        hm_best_score = hm_get(hm, 'ct3_best_viral_score', 'rvdb_best_viral_score')
        hm_classes = hm_get(hm, 'ct3_tier_counts', 'rvdb_class_counts')
        # 裁决用计数
        hm_n_specific = hm_n_family if args.hmm_require_family else hm_n_viral
        if args.hmm_count_ambiguous and not args.hmm_require_family:
            hm_n_specific += hm_n_ambiguous

        # 加权评分 (30/30/40); 域证据层 = max(CDD, HMM):
        #   CDD 盲区 (NON_VIRAL/REVIEW 但 CT3-HMM 病毒命中强) 由 HMM 补上,
        #   CDD 命中时取 max 不叠加 — 两者同源 (profile 证据), 相加会双重计分。
        #   --no-hmm-in-score 时 hmm_score 置 0, domain 退回纯 CDD (旧口径)。
        blastn_score = (nt_pident * 0.7 + nt_qcovs * 0.3) if nt_pident > 0 else 0.0
        blastx_score = aa_pident if aa_pident > 0 else 0.0
        hmm_score = hmm_evidence_score(hm_best_tier, hm_n_specific) if args.hmm_in_score else 0.0
        domain_score = max(cdd_score, hmm_score)
        total_score = 0.3 * blastn_score + 0.3 * blastx_score + 0.4 * domain_score

        # ═ 判定框架 (三层职责 + 长度门槛) ═
        # 1. CDD 结构域 → 判是不是病毒
        # 2. blastx/blastn 蛋白/核酸命中 → 判已知还是新病毒 (仅对"是病毒"的)
        # 3. CDD 核心域(科/属) → 核验 taxonomy 分类是否正确
        # 长度门槛: core domain validation只在 长度>=属平均长度×20% 时起作用 (short fragmentmay lack core domain)
        has_tax_class = bool(tax_family or tax_genus or tax_species)
        # 长度门槛: 查该属平均长度; 够长才做core domain validation
        genus_avg = genus_len_map.get(tax_genus, 0.0)
        length_f = float(length) if length else 0.0
        long_enough = True  # 默认够长
        if genus_avg and length_f:
            long_enough = length_f >= genus_avg * 0.20
        # ── 1. CDD 判是不是病毒 ──
        is_viral_dom = cdd_evid in ('PASS_CORE', 'PASS_CORE_FAM', 'PASS_VIRAL')  # 有病毒结构域证据
        has_virus_core = cdd_evid in ('PASS_CORE', 'PASS_CORE_FAM')              # 命中 taxonomy 属/科核心域(分类被支持)
        # ── 2. blast 判已知/新病毒 (blastx 蛋白命中为主, blastn 辅; 覆盖率用 subject 覆盖) ──
        # 覆盖率: blastx=aa_qcov(subject覆盖), blastn=nt_qcovs
        bx_known = (aa_pident >= 85 and aa_qcov >= 85) or (nt_pident >= 85 and nt_qcovs >= 85)
        bx_new   = (aa_pident >= 40 and aa_qcov >= 50) or (nt_pident >= 40 and nt_qcovs >= 50)
        # ── 判定 ──
        if is_viral_dom:
            # CDD 说是病毒结构域; 用 blast 判已知/新病毒
            if has_virus_core:
                verdict = 'KEEP'
            else:
                # PASS_VIRAL: 有病毒域但未命中 tax 属核心域
                if long_enough:
                    verdict = 'KEEP'
                else:
                    verdict = 'REVIEW'
            if bx_known:
                novelty = 'Known virus(blast蛋白)'
            elif bx_new:
                novelty = 'Novel candidate (blast protein)'
            else:
                novelty = 'Novel candidate (domain-only evidence)'
            if not has_virus_core:
                if long_enough:
                    fp_reason = 'tax_mismatch: 有病毒域但CDD未命中taxenom属核心域(分类未被CDD支持)'
                else:
                    fp_reason = 'short_fragment: short fragment (<20% genus mean) excluded from core domain validation, hits retained for review'
            else:
                fp_reason = ''
        else:
            # CDD 无病毒结构域 → 可能是宿主Contamination/short fragment
            if long_enough:
                # 够长 + no viral domain: 看 blast
                if bx_new:
                    verdict = 'REVIEW'  # no viral domain但有blast, 可疑
                    novelty = 'Novel candidate (blast only, no viral domain)'
                    fp_reason = 'no_viral_domain_but_blast: CDDno viral domain but blast hit, pending'
                else:
                    verdict = 'DROP'  # 假阳性(宿主Contamination/无证据)
                    _tax = f"{tax_family or '?'}/{tax_genus or '?'}/{tax_species or '?'}" if has_tax_class else '无分类/NA'
                    novelty = 'Host / non-viral'
                    fp_reason = f'FALSE_POSITIVE: taxonomy 分类为 {_tax} 但 CDD 无病毒核心域且 blast 无有效证据(aa={aa_pident:.0f}%/{aa_qcov:.0f}%cov)'
            else:
                # short fragment(长度<属均值20%) → may lack core domain, 不武断, keep REVIEW
                verdict = 'REVIEW'
                novelty = 'Short fragment - pending'
                fp_reason = 'short_fragment: length below 20% of genus mean, may lack core domain, pending'

        # ── 第三路探针提升 (默认开, --no-hmm-promote-keep 关闭) ──
        # 口径: "比对上就 KEEP"。CDD 无病毒域但 CT3-HMM 有非噬菌体病毒命中时,
        # 把 DROP/REVIEW 直接提升为 KEEP。
        hmm_rescue = ''
        hmm_promoted = False
        if args.hmm_promote_keep and hm_n_specific >= args.hmm_promote_min:
            if verdict in ('DROP', 'REVIEW'):
                _prev_v = verdict
                verdict = 'KEEP'
                hmm_rescue = f'{_prev_v}->KEEP'
                novelty = 'Novel candidate (CT3-HMM viral profile)'
                fp_reason = ('hmm_promote: CDD 无病毒域, 但 %s 有 %d 条命中 '
                             '(best=%s [%s/%s] E=%s)' % (gate_desc, hm_n_specific,
                                                         hm_best_model, hm_best_tier,
                                                         hm_best_cluster, hm_best_ev))
                hmm_promoted = True

        # ── 第三路探针 rescue (DROP→REVIEW; 仅在未升 KEEP 时生效) ──
        if (not hmm_promoted) and args.hmm_rescue and (not is_viral_dom) and hm_n_specific >= args.hmm_min_specific:
            if verdict == 'DROP':
                verdict = 'REVIEW'
                hmm_rescue = 'DROP->REVIEW'
                novelty = 'Novel candidate (CT3-HMM viral profile, no CDD domain)'
                fp_reason = ('hmm_rescue: CDD 无病毒域, 但 %s 有 %d 条命中 '
                             '(best=%s [%s/%s] E=%s)' % (gate_desc, hm_n_specific,
                                                         hm_best_model, hm_best_tier,
                                                         hm_best_cluster, hm_best_ev))
            elif verdict == 'REVIEW':
                hmm_rescue = 'REVIEW(kept)'

        rows.append({
            'contig_id': cid, 'length': length, 'category': category,
            'cdd_evidence': cdd_evid, 'n_t1': cdd.get('n_t1', ''), 'cdd_top': cdd_top,
            'nt_pident': f'{nt_pident:.1f}' if nt_pident else '',
            'nt_qcovs': f'{nt_qcovs:.1f}' if nt_qcovs else '',
            'nt_species': nt_species,
            'aa_pident': f'{aa_pident:.1f}' if aa_pident else '',
            'aa_qcov': f'{aa_qcov:.1f}' if aa_qcov else '',
            'aa_species': aa_species,
            'tax_genus': tax_genus, 'tax_family': tax_family, 'tax_species': tax_species,
            'genus_avg_len': f'{genus_avg:.1f}' if genus_avg else (genus_avg_len or ''),
            'n_samples': n_samples, 'prevalence_%': prevalence,
            'checkv_completeness': f'{completeness:.1f}' if completeness else '',
            'checkv_confidence': confidence,
            'blastn_score': f'{blastn_score:.1f}', 'blastx_score': f'{blastx_score:.1f}',
            'cdd_score': f'{cdd_score:.1f}', 'hmm_score': f'{hmm_score:.1f}',
            'domain_score': f'{domain_score:.1f}', 'total_score': f'{total_score:.1f}',
            'novelty': novelty, 'verdict': verdict, 'fp_reason': fp_reason,
            'ct3_hmm_evidence': hm_evidence, 'ct3_viral_evidence': hm_viral,
            'ct3_n_viral': hm_n_viral, 'ct3_n_family': hm_n_family,
            'ct3_n_ambiguous': hm_n_ambiguous,
            'ct3_best_viral_model': hm_best_model, 'ct3_best_viral_cluster': hm_best_cluster,
            'ct3_best_viral_tier': hm_best_tier, 'ct3_best_viral_db': hm_best_db,
            'ct3_best_viral_desc': hm_best_desc,
            'ct3_best_viral_evalue': hm_best_ev, 'ct3_best_viral_score': hm_best_score,
            'ct3_tier_counts': hm_classes, 'hmm_rescue': hmm_rescue,
        })

    rows.sort(key=lambda x: -float(x['total_score']))

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, 'rescue_evidence_scored.tsv')
    fieldnames = ['contig_id','length','category','cdd_evidence','n_t1','cdd_top',
                  'nt_pident','nt_qcovs','nt_species','aa_pident','aa_qcov','aa_species',
                  'tax_genus','tax_family','tax_species','genus_avg_len',
                  'n_samples','prevalence_%',
                  'checkv_completeness','checkv_confidence',
                  'blastn_score','blastx_score','cdd_score','hmm_score','domain_score','total_score',
                  'novelty','verdict','fp_reason',
                  'ct3_hmm_evidence','ct3_viral_evidence','ct3_n_viral','ct3_n_family',
                  'ct3_n_ambiguous',
                  'ct3_best_viral_model','ct3_best_viral_cluster','ct3_best_viral_tier',
                  'ct3_best_viral_db','ct3_best_viral_desc',
                  'ct3_best_viral_evalue','ct3_best_viral_score','ct3_tier_counts',
                  'hmm_rescue']
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f'输出: {out} ({len(rows)} 条)')

    verdict_cnt = Counter(r['verdict'] for r in rows)
    novelty_cnt = Counter(r['novelty'] for r in rows)
    print('\n=== 最终判定 ===')
    for k in ['KEEP', 'REVIEW', 'DROP']:
        print(f'  {k}: {verdict_cnt.get(k, 0)}')
    n_rescued = sum(1 for r in rows if r['hmm_rescue'] == 'DROP->REVIEW')
    print(f'  其中 CT3-HMM rescue 抬升 DROP->REVIEW: {n_rescued}')
    n_promoted = sum(1 for r in rows if r['hmm_rescue'].endswith('->KEEP'))
    print(f'  其中 CT3-HMM 比对命中提升 KEEP: {n_promoted}')
    print('\n=== 新病毒标注 ===')
    for k, v in novelty_cnt.most_common():
        print(f'  {k}: {v}')

    # ---------- 输出 keep/review/drop 的 txt + fasta ----------
    for verdict in ['KEEP', 'REVIEW', 'DROP']:
        sub = [r for r in rows if r['verdict'] == verdict]
        ids = [r['contig_id'] for r in sub]
        # txt 列表
        txt_path = os.path.join(args.outdir, f'{verdict.lower()}_candidates.txt')
        with open(txt_path, 'w') as f:
            for cid in ids:
                f.write(cid + '\n')
        # fasta (若有序列)
        if seqs:
            fa_path = os.path.join(args.outdir, f'{verdict.lower()}_candidates.fasta')
            n_written = 0
            with open(fa_path, 'w') as f:
                for cid in ids:
                    if cid in seqs:
                        f.write(f'>{cid}\n{seqs[cid]}\n')
                        n_written += 1
            print(f'  {verdict}: {len(ids)} 条 → {verdict.lower()}_candidates.txt (+ {n_written} 条 fasta)')

    # ---------- 第三路探针 (CT3-HMM) rescue 候选清单 ----------
    # 无论是否开启 --hmm-rescue 都产出, 供人工过目: 非 KEEP 且带病毒命中的 contig
    if args.hmm_report:
        cand = [r for r in rows
                if r['verdict'] != 'KEEP'
                and int(r.get('ct3_n_viral', 0) or 0) >= 1]
        cfields = ['contig_id', 'verdict', 'hmm_would_change', 'cdd_evidence', 'category',
                   'length', 'aa_pident', 'aa_qcov', 'aa_species', 'nt_pident', 'nt_qcovs', 'nt_species',
                   'ct3_n_viral', 'ct3_n_family', 'ct3_n_ambiguous',
                   'ct3_best_viral_model', 'ct3_best_viral_cluster', 'ct3_best_viral_tier',
                   'ct3_best_viral_db', 'ct3_best_viral_evalue', 'ct3_best_viral_score',
                   'ct3_tier_counts']
        cpath = os.path.join(args.outdir, 'hmm_ct3_rescue_candidates.tsv')
        with open(cpath, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cfields, delimiter='\t')
            w.writeheader()
            for r in sorted(cand, key=lambda x: -int(x.get('ct3_n_viral', 0) or 0)):
                row = {k: r.get(k, '') for k in cfields}
                if r.get('hmm_rescue') == 'DROP->REVIEW':
                    row['hmm_would_change'] = 'applied(DROP->REVIEW)'
                elif r['verdict'] == 'DROP':
                    row['hmm_would_change'] = 'DROP->REVIEW'
                else:
                    row['hmm_would_change'] = ''
                w.writerow(row)
        print(f'\nCT3-HMM rescue 候选 (非 KEEP 且有病毒命中): {len(cand)} 条 → {cpath}')


if __name__ == '__main__':
    main()
