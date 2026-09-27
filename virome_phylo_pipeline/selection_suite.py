#!/usr/bin/env python3
"""selection_suite.py — 选择分析套件（通用, 所有基因/所有病毒）

子命令:
  fubar        贝叶斯位点级正选择 (HyPhy FUBAR) + 与 FEL/MEME 交叉表
  relax        选择强度放松/收紧 (HyPhy RELAX, 需分组)
  contrastfel  宿主分组位点级对比 (HyPhy CONTRASTFEL, 需分组)
  timetree     正选择位点 × BEAST 时间树            (待实现)
  network      正选择等位基因 × 单倍型网络          (待实现)
  structure    正选择位点 × 蛋白结构/结构域映射     (待实现)

注: codeml M7/M8 子命令已于 2026-08-30 归档止损 (正选择证据链唯一归属 CAPHEINE)。
    唯一实现 codeml_bridge.py 在 _archive_codeml_20260830/, 已不在 utils/ 下;
    如需复活, 从归档目录取回并重新接入 main() 分发, 不要再引用 utils.codeml_bridge。

用法:
  python selection_suite fubar \
      --capheine-dir <virus>/capheine/ \
      [--genes G,L] [--outdir selection_out] [--threads 8]

输入 (CAPHEINE 标准目录结构):
  cawlign/ref_cds-noStopCodons.part_{GENE}-aligned.fasta  密码子比对
  iqtree/ref_cds-noStopCodons.part_{GENE}.treefile        ML 树
"""
import argparse, json, os, subprocess, sys, shutil
from pathlib import Path


def find_hyphy():
    for cmd in ['hyphy']:
        p = shutil.which(cmd)
        if p:
            return p
    return None


def discover_partitions(capheine_dir):
    """自动发现基因分区: 返回 {gene: {'aln': path, 'tree': path}}

    优先 CLN/*-nodups.fasta (CAPHEINE 去重后输入, 与 iqtree 树匹配)
    回退 cawlign/*-aligned.fasta
    """
    base = Path(capheine_dir)
    tree_dir = base / 'iqtree'
    genes = {}
    if not tree_dir.exists():
        return genes
    # 候选比对目录 (优先级): hyphy/CLN (CAPHEINE 去重输入) > CLN > cawlign
    aln_dirs = []
    for cand in [base / 'hyphy' / 'CLN', base / 'CLN', base / 'cawlign']:
        if cand.exists():
            aln_dirs.append(cand)
    for aln_dir in aln_dirs:
        for aln in sorted(aln_dir.glob('*part_*.fasta')):
            name = aln.name
            gene = None
            if '-nodups.fasta' in name:
                gene = name.split('.part_')[-1].replace('-nodups.fasta', '')
            elif '-aligned.fasta' in name:
                gene = name.split('.part_')[-1].replace('-aligned.fasta', '')
            elif '.part_' in name and name.endswith('.fasta'):
                gene = name.split('.part_')[-1].replace('.fasta', '')
            if not gene or gene in genes:
                continue  # 已有更高优先级目录的条目
            tree = tree_dir / f'ref_cds-noStopCodons.part_{gene}.treefile'
            if tree.exists():
                genes[gene] = {'aln': str(aln), 'tree': str(tree)}
    return genes


def run_hyphy(method, aln, tree, out_json, hyphy, extra=None):
    """运行 HyPhy 方法

    2026-09-15 修复:
      · 未捕获 TimeoutExpired → 超时直接抛异常中断整条流程;
      · text=True 未指定 encoding → Linux LANG=C 下 locale 为 ASCII,
        HyPhy 输出含非 ASCII 即 UnicodeDecodeError。
    超时返回 returncode=124 的 CompletedProcess, 保持调用方契约不变。
    """
    cmd = [hyphy, method, '--alignment', aln, '--tree', tree,
           '--output', out_json]
    if extra:
        cmd += extra
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=3600)
    except subprocess.TimeoutExpired:
        r = subprocess.CompletedProcess(cmd, 124, '',
                                        f'{method} 超时 (>3600s)')
    return r


# ---------------- FUBAR ----------------
def parse_fubar(json_path, gene):
    d = json.load(open(json_path))
    mle = d.get('MLE', {})
    headers = [h[0] for h in mle.get('headers', [])]
    content = mle.get('content', {})
    sites = []
    for part in content.values():
        for row in part:
            sites.append(dict(zip(headers, row)))
    # FUBAR: 正选择 = Prob[alpha<beta] >= 0.9 (HyPhy 推荐阈值)
    # 注意: headers 6 列但数据行 8 列 (尾部 0,0), zip 取前 6; 且只用正选择后验列
    pos_idx = None
    for i, h0 in enumerate(headers):
        if 'positive selection' in h0.lower() or 'prob[alpha<beta]' in h0.lower():
            pos_idx = i
            break
    hits = []
    for si, s in enumerate(sites, 1):
        if pos_idx is None:
            break
        pv = s.get(headers[pos_idx]) if isinstance(s, dict) else None
        if isinstance(pv, (int, float)) and pv >= 0.9:
            hits.append({'site': si, 'alpha': s.get('alpha'),
                         'beta': s.get('beta'), 'prob': pv})
    return {'gene': gene, 'n_sites': len(sites), 'n_positive': len(hits),
            'hits': hits}


def cmd_fubar(args):
    hyphy = find_hyphy()
    if not hyphy:
        sys.exit('hyphy 未找到')
    genes = discover_partitions(args.capheine_dir)
    if args.genes:
        genes = {g: v for g, v in genes.items() if g in args.genes.split(',')}
    if not genes:
        sys.exit(f'未发现基因分区: {args.capheine_dir}')
    outdir = Path(args.outdir) / 'fubar'
    outdir.mkdir(parents=True, exist_ok=True)

    print(f'=== FUBAR (贝叶斯位点级) === 基因: {list(genes.keys())}')
    results = {}
    for gene, paths in genes.items():
        out_json = outdir / f'{gene}.FUBAR.json'
        if not out_json.exists():
            r = run_hyphy('fubar', paths['aln'], paths['tree'],
                          str(out_json), hyphy)
            if r.returncode != 0:
                print(f'  {gene}: 失败 ({r.stderr[-300:]})')
                continue
        try:
            res = parse_fubar(out_json, gene)
            results[gene] = res
            print(f"  {gene}: {res['n_positive']}/{res['n_sites']} 正选择位点 "
                  f"(Prob≥0.9)" + (f" → site {[h['site'] for h in res['hits']]}" if res['hits'] else ''))
        except Exception as e:
            print(f'  {gene}: 解析失败 {e}')

    # 与 FEL/MEME 交叉表 (若已有 selection_results.json)
    cross = {'fubar': results}
    sel_json = Path(args.outdir) / 'selection_results.json'
    if sel_json.exists():
        sel = json.load(open(sel_json))
        for method in ['fel', 'meme']:
            cross[method] = {}
            for m in sel.get(method, []):
                cross[method][m['gene']] = {h['site'] for h in m['hits']}
        print('\n=== FUBAR × FEL × MEME 交叉 (site 编号) ===')
        for gene in results:
            fu = {h['site'] for h in results[gene]['hits']}
            fe = cross['fel'].get(gene, set())
            me = cross['meme'].get(gene, set())
            triple = fu & fe & me
            print(f'  {gene}: FUBAR={sorted(fu)} FEL={sorted(fe)} MEME={sorted(me)} '
                  f'→ 三方法共检出={sorted(triple)}')

    with open(outdir / 'fubar_summary.json', 'w') as f:
        json.dump(cross, f, indent=2, default=str)
    print(f'\n结果: {outdir}/fubar_summary.json')
    return 0


# ---------------- 宿主分组 (RELAX / CONTRASTFEL) ----------------
def load_host_groups(metadata_csv, test_key='ruthenicum', ref_key='barbarum'):
    """元数据 → {run: group}, group ∈ {test, ref}"""
    import csv as _csv
    groups = {}
    with open(metadata_csv, encoding='utf-8-sig') as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = _csv.Sniffer().sniff(sample, delimiters='\t,')
        except Exception:
            dialect = _csv.excel_tab
        for row in _csv.DictReader(f, dialect=dialect):
            run = (row.get('Run') or '').strip()
            sp = (row.get('ScientificName') or '').strip().lower()
            if test_key in sp:
                groups[run] = 'test'
            elif ref_key in sp:
                groups[run] = 'ref'
    return groups


def host_tips(aln_fasta, groups):
    """比对 → {tip_id: group|None}"""
    from Bio import SeqIO
    import re as _re
    tips = {}
    for rec in SeqIO.parse(aln_fasta, 'fasta'):
        # 2026-09-15 修复: 旧正则 (CRR|SRR) 漏掉 SRA 的 ERR/DRR 编号 →
        # 欧/美来源样本全归 None, 分组数不足后被静默跳过 (整块分析无声消失)。
        # 与 mask_recombination.py 的 (CRR|SRR|ERR|DRR) 口径对齐。
        m = _re.match(r'(CRR|SRR|ERR|DRR)\d+', rec.id)
        g = groups.get(m.group(0), None) if m else None
        tips[rec.id] = g
    return tips


def annotate_tree(tree_path, tips, out_path):
    """树 tip 加 {group} 标注 (HyPhy 语法): name{test}/name{reference}"""
    import re as _re
    txt = Path(tree_path).read_text().strip()
    def _ann(m):
        name = m.group(1)
        g = tips.get(name)
        if g:
            return f'{name}{{{g}}}:'
        return f'{name}:'
    # tip 名 = 冒号前无括号的标识符
    txt2 = _re.sub(r'([^(),:]+?)\s*:', _ann, txt)
    Path(out_path).write_text(txt2)
    return out_path


def run_group_method(method, aln, tree, test_tips, ref_tips, out_json, hyphy):
    """RELAX / CONTRASTFEL: 宿主分组对比 (树文件 {group} 标注)"""
    import tempfile
    tips = {t: 'test' for t in test_tips}
    tips.update({t: 'reference' for t in ref_tips})
    tree_ann = annotate_tree(tree, tips, f'{out_json}.ann.tree')
    if method == 'contrast-fel':
        # CONTRASTFEL 用 --branch-set (test 组 vs 背景), 不用 --test/--reference
        cmd = [hyphy, method, '--alignment', aln, '--tree', tree_ann,
               '--branch-set', 'test', '--output', out_json]
    else:
        cmd = [hyphy, method, '--alignment', aln, '--tree', tree_ann,
               '--test', 'test', '--reference', 'reference',
               '--output', out_json]
    # 2026-09-15: 补 encoding + 捕获 TimeoutExpired (同 run_hyphy)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=3600)
    except subprocess.TimeoutExpired:
        r = subprocess.CompletedProcess(cmd, 124, '',
                                        f'{method} 超时 (>3600s)')
    return r


def cmd_relax(args):
    return _cmd_grouped('relax', args)


def cmd_contrastfel(args):
    return _cmd_grouped('contrast-fel', args)


def _cmd_grouped(method, args):
    hyphy = find_hyphy()
    if not hyphy:
        sys.exit('hyphy 未找到')
    if not args.metadata:
        sys.exit(f'{method} 需要 --metadata (Global_Unified_Metadata CSV)')
    genes = discover_partitions(args.capheine_dir)
    if args.genes:
        genes = {g: v for g, v in genes.items() if g in args.genes.split(',')}
    if not genes:
        sys.exit(f'未发现基因分区: {args.capheine_dir}')
    groups = load_host_groups(args.metadata, args.test_host, args.ref_host)
    outdir = Path(args.outdir) / ('relax' if method == 'relax' else 'contrastfel')
    outdir.mkdir(parents=True, exist_ok=True)
    print(f'=== {method.upper()} (宿主对比: {args.test_host} vs {args.ref_host}) ===')
    print(f'宿主分组: test={sum(1 for v in groups.values() if v=="test")} '
          f'ref={sum(1 for v in groups.values() if v=="ref")}')
    results = {}
    for gene, paths in genes.items():
        tips = host_tips(paths['aln'], groups)
        test_tips = [t for t, g in tips.items() if g == 'test']
        ref_tips = [t for t, g in tips.items() if g == 'ref']
        if not test_tips or not ref_tips:
            print(f'  {gene}: 分组不足 (test={len(test_tips)} ref={len(ref_tips)})')
            continue
        out_json = outdir / f'{gene}.{method.upper()}.json'
        if not out_json.exists():
            r = run_group_method(method, paths['aln'], paths['tree'],
                                 test_tips, ref_tips, str(out_json), hyphy)
            if r.returncode != 0:
                print(f'  {gene}: 运行失败 ({r.stderr[-300:]})')
                continue
        try:
            d = json.load(open(out_json))
            tr = d.get('test results', {})
            if method == 'relax':
                k = tr.get('relaxation or intensification parameter')
                pv = tr.get('p-value')
                verdict = ('intensified' if k and k > 1 else 'relaxed') if k else '?'
                print(f'  {gene}: k={k} p={pv} {verdict}')
                results[gene] = {'k': k, 'p_value': pv, 'verdict': verdict}
            else:
                print(f'  {gene}: {json.dumps(tr)[:200]}')
                results[gene] = {'test_results': tr}
        except Exception as e:
            print(f'  {gene}: 解析失败 {e}')
    with open(outdir / 'summary.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f'\n结果: {outdir}/summary.json')
    return 0


# ---------------- 子命令分发 ----------------
def main():
    ap = argparse.ArgumentParser(description='选择分析套件')
    sub = ap.add_subparsers(dest='command', required=True)

    for name, help_ in [('fubar', 'FUBAR 位点级'), ('relax', 'RELAX 选择强度'),
                        ('contrastfel', 'CONTRASTFEL 分组对比'),
                        ('codeml', 'M7/M8 LRT (已归档, 未接入)'),
                        ('timetree', '位点×时间树 (待实现)'), ('network', '位点×网络 (待实现)'),
                        ('structure', '位点×结构 (待实现)')]:
        p = sub.add_parser(name, help=help_)
        p.add_argument('--capheine-dir', required=True, help='CAPHEINE 病毒输出目录 (如 <virus>/capheine/)')
        p.add_argument('--genes', default=None, help='基因列表 (逗号分隔, 默认全部)')
        p.add_argument('--outdir', default='selection_out')
        p.add_argument('--threads', default='8')
        p.add_argument('--metadata', default=None, help='宿主元数据 CSV (relax/contrastfel 需要)')
        p.add_argument('--test-host', default='ruthenicum', help='测试组宿主 (ScientificName 子串, 默认黑果枸杞)')
        p.add_argument('--ref-host', default='barbarum', help='参考组宿主 (默认枸杞)')

    args = ap.parse_args()
    if args.command == 'fubar':
        return cmd_fubar(args)
    if args.command == 'relax':
        return cmd_relax(args)
    if args.command == 'contrastfel':
        return cmd_contrastfel(args)
    # 其余子命令 (codeml/timetree/network/structure) 未接入分发。
    # codeml 追加归档说明: 2026-08-30 决策 — 正选择证据链以 CAPHEINE 为主;
    # 实现见 _archive_codeml_20260830/codeml_bridge.py (不在 utils/ 下, 勿再 import utils.codeml_bridge)。
    print(f'[待实现] {args.command} — 逐步开发中')
    return 0


# cmd_codeml 已于 2026-08-30 随 codeml 链路一并归档移除。
# 旧版此函数是死代码 (main() 不分发 codeml, 永远走不到), 且从已删除的 utils.codeml_bridge
# 导入, 一旦有人接回分发即 ImportError。实现保留在 _archive_codeml_20260830/codeml_bridge.py。


if __name__ == '__main__':
    sys.exit(main())
