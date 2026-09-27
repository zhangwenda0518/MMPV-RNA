#!/usr/bin/env python3
"""selection_analysis.py — 选择压力分析套件（通用, 任何病毒）

子命令:
  parse        解析已有 CAPHEINE 结果 (BUSTED/FEL/MEME)
  fubar        HyPhy FUBAR 位点级交叉验证 (贝叶斯)
  relax        HyPhy RELAX 选择强度 (需宿主/分组)
  contrastfel  HyPhy CONTRASTFEL 分组位点对比 (需宿主/分组)
  m7m8         codeml M7/M8 LRT 金标准 (PAML)
  time-tree    正选择位点 × BEAST 时间树 (图5式)
  network      正选择位点 × MJN 网络 (图6式)
  structure    正选择位点 × 蛋白域/结构 (图7式)

通用输入: CDS 分区比对 (part_G.fasta 等) + 基因坐标 CSV + 树 + 分组元数据
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ============================================================
# 工具函数
# ============================================================

def find_hyphy():
    # 2026-08-30 审计修复: 旧版用 Unix `which` 子进程，Windows 无 which → 永远 None。
    # 改用 shutil.which 跨平台。
    import shutil
    for cmd in ['hyphy']:
        p = shutil.which(cmd)
        if p:
            return p
    return None


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR 校正"""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    qs = [0.0] * n
    for rank, i in enumerate(order, 1):
        qs[i] = pvals[i] * n / rank
    for k in range(n - 2, -1, -1):
        qs[order[k]] = min(qs[order[k]], qs[order[k + 1]])
    return qs


def load_genes(genes_csv):
    """基因注释 CSV → {gene: (start, end)}"""
    import csv
    out = {}
    if genes_csv and os.path.exists(genes_csv):
        with open(genes_csv) as f:
            for row in csv.DictReader(f):
                if row['gene'] not in out:
                    out[row['gene']] = (int(row['start']), int(row['end']))
    return out


def gene_to_genome(hits, gene, gene_map):
    """位点(密码子编号) → 基因组坐标"""
    gstart = gene_map.get(gene, (None, None))[0]
    if gstart is None:
        return hits
    for h in hits:
        h['genome_pos'] = gstart + (h['site'] - 1) * 3
    return hits


def load_group_map(group_csv, seq_column='name', group_column='host'):
    """分组元数据 → {序列名: 组名}"""
    import csv
    out = {}
    if not group_csv or not os.path.exists(group_csv):
        return out
    with open(group_csv) as f:
        for row in csv.DictReader(f):
            nm = row.get(seq_column, '')
            g = row.get(group_column, '')
            if nm and g:
                out[nm] = g
    return out


def _strip_internal_node_labels(nwk):
    """删掉 newick 内部节点的标签 (如 IQ-TREE 的 `82/97` 支持值), 保留枝叶名与枝长。

    为什么需要 (2026-09-16): `build_group_tree` 的产物直接喂 HyPhy
    (`run_hyphy_method`), 而 HyPhy 的分组标注用 `{group}` 花括号。花括号紧跟一个
    既有节点标签 (`...){group}82/97:1`) 时 HyPhy 的解析行为不可预期 —— 主树在
    2026-09-16 补了分支支持值后, 这种"既有标签"必然出现 (补之前没有)。
    分组标注只需要**拓扑 + 枝长**, 所以这里先把内部标签清掉, 恢复成 HyPhy
    一直被验证过的那个输入形状。
    """
    spans = []
    i, n = 0, len(nwk)
    while i < n:
        c = nwk[i]
        if c in '(),:;':
            i += 1
            continue
        j = i
        while j < n and nwk[j] not in '(),:;':
            j += 1
        # 紧跟在 `)` 之后的 token = 内部节点标签 (跳过空白判断)
        k = i - 1
        while k >= 0 and nwk[k] in ' \t\r\n':
            k -= 1
        if k >= 0 and nwk[k] == ')':
            spans.append((i, j))
        i = j
    if not spans:
        return nwk
    out, prev = [], 0
    for a, b in spans:
        out.append(nwk[prev:a])
        prev = b
    out.append(nwk[prev:])
    return ''.join(out)


def _annotate_clade_text(nwk, tip_names, clade_name):
    """newick 文本级 clade 标注: 目标叶 LCA 的右括号后插 {clade_name}

    栈法: 记录每个叶的祖先 '(' 位置 → 公共前缀 = LCA

    ⚠️ 2026-09-16 修 (P0: 死循环): 原实现把 `,` `:` `;` 这些**分隔符**也丢进 else
    分支, 而那个分支先用 `while nwk[j] not in '(),:;'` 找 token 结束位 —— 首字符
    就是分隔符时该循环立刻为假, `j == i`, 于是 `i = j` **原地踏步 → 对任何 newick
    都永久挂死**(第 2 个字符必然是 `(`/逗号/冒号)。调用链
    `build_group_tree → _annotate_clade_text` 属于 selection_analysis 的
    RELAX/CONTRASTFEL 分组标注 → 整条选择压力分析**永久卡住不返回**。
    """
    stack = []
    leaf_paths = {}
    i, n = 0, len(nwk)
    while i < n:
        c = nwk[i]
        if c == '(':
            stack.append(i)
            i += 1
        elif c == ')':
            if stack:
                stack.pop()
            i += 1
        elif c in ',:;':
            i += 1                          # 分隔符: 必须自己前进 (否则死循环)
        else:
            j = i
            while j < n and nwk[j] not in '(),:;':
                j += 1
            name = nwk[i:j].strip()
            # 紧跟在 `)` 后面的名字是**内部节点标签**(如 IQ-TREE 的 `82/97` 支持值),
            # 不是叶子 —— 主树补支持值后这类 token 必然出现, 必须排除,
            # 否则会被当成"叶子"塞进 leaf_paths (污染 LCA 推断)。
            if name and nwk[i - 1:i] != ')':
                leaf_paths[name] = list(stack)
            i = j
    paths = [leaf_paths[t] for t in tip_names if t in leaf_paths]
    if not paths:
        return nwk, False
    common = []
    for k in range(min(len(p) for p in paths)):
        if all(p[k] == paths[0][k] for p in paths):
            common.append(paths[0][k])
        else:
            break
    if not common:
        return nwk, False
    start = common[-1]
    depth = 0
    end = None
    for i in range(start, n):
        if nwk[i] == '(':
            depth += 1
        elif nwk[i] == ')':
            depth -= 1
            if depth == 0:
                end = i
                break
    if end is None:
        return nwk, False
    return nwk[:end + 1] + '{' + clade_name + '}' + nwk[end + 1:], True


def build_group_tree(tree_path, group_map, out_path, group_names):
    """构建带分组标注的树: 目标组叶的 LCA 子树标注 {group_<组名>}

    HyPhy RELAX/CONTRASTFEL 的 --test 用花括号分支标注

    2026-09-16: 标注前先清掉内部节点标签 —— 主树补了 IQ-TREE 支持值后,
    内部节点会带 `82/97` 这类标签; 与 HyPhy 的 `{group}` 花括号混在一个
    节点上解析行为不可预期 (见 `_strip_internal_node_labels` 注释)。
    """
    from Bio import Phylo
    nwk = _strip_internal_node_labels(Path(tree_path).read_text().strip())
    tree = Phylo.read(tree_path, 'newick')
    test_names = []
    n_renamed = 0
    for gname in group_names:
        leaves = []
        for leaf in tree.get_terminals():
            nm = leaf.name
            base = nm.split('_OR489165')[0] if '_OR' in nm else nm.split('/')[0]
            g = group_map.get(nm) or group_map.get(base) or group_map.get(nm.split('/')[-1])
            if g == gname:
                leaves.append(leaf.name)
        if not leaves:
            continue
        clade_name = f'group_{gname.replace(" ", "_")}'
        nwk, ok = _annotate_clade_text(nwk, leaves, clade_name)
        if ok:
            test_names.append(clade_name)
            n_renamed += len(leaves)
    Path(out_path).write_text(nwk + '\n')
    test_re = '|'.join(test_names) if test_names else None
    return out_path, test_re, n_renamed


# ============================================================
# 子命令: FUBAR / RELAX / CONTRASTFEL (HyPhy 执行)
# ============================================================

def run_hyphy_method(method, aligned_fasta, tree_path, outdir, test_re=None):
    """执行 HyPhy 单方法, 返回输出 JSON 路径"""
    hyphy = find_hyphy()
    if not hyphy:
        print('[错误] 未找到 hyphy')
        return None
    os.makedirs(outdir, exist_ok=True)
    out_json = os.path.join(outdir, f'{method}.json')
    cmd = [hyphy, method, '--alignment', aligned_fasta,
           '--tree', tree_path, '--output', out_json]
    if test_re and method in ('relax', 'contrastfel'):
        cmd += ['--test', test_re]
    print(f'[run] {" ".join(cmd)}')
    # 2026-09-15 修复: 旧版未捕获 TimeoutExpired (超时抛异常中断整条流程),
    # 且 text=True 未指定 encoding (Linux LANG=C 下 HyPhy 非 ASCII 输出会
    # UnicodeDecodeError)。两处一并补上。
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=3600)
    except subprocess.TimeoutExpired:
        print(f'[警告] {method} 超时 (>3600s), 跳过')
        return None
    if r.returncode != 0:
        print(f'[警告] {method} 退出码 {r.returncode}: {(r.stderr or "")[-600:]}')
        return None
    return out_json if os.path.exists(out_json) else None


def parse_fubar(json_path, gene):
    """FUBAR JSON → 位点表 (posterior prob)"""
    with open(json_path, encoding='utf-8', errors='replace') as _f:
        d = json.load(_f)
    mle = d['MLE']
    headers = [h[0] for h in mle['headers']]
    content = mle['content']
    sites = []
    for part in content.values():
        for row in part:
            sites.append(dict(zip(headers, row)))
    hits = []
    # 2026-09-15 修复: 正选择 = dN>dS 即 beta>alpha, 对应后验列 Prob[alpha<beta]。
    # 旧实现找的是 'Prob[alpha>beta]' (纯化选择列), 且元组两个元素完全相同
    # → 把纯化选择位点当正选择上报, 与 selection_suite 的 Prob[alpha<beta] 结论相反。
    # 现按名称优先取正选择列; 名称缺失时仅在"恰好只有一个 prob 列"时兜底,
    # 避免在多列情况下误取纯化选择列。
    for si, s in enumerate(sites, 1):
        pp = None
        for k, v in s.items():
            if not isinstance(v, (int, float)):
                continue
            kl = k.lower()
            if 'positive selection' in kl or 'prob[alpha<beta]' in kl:
                pp = v
                break
        if pp is None:
            _probs = [(k, v) for k, v in s.items()
                      if 'prob' in k.lower() and isinstance(v, (int, float))]
            if len(_probs) == 1:
                pp = _probs[0][1]
        if pp is not None and pp >= 0.9:
            hits.append({'site': si, 'posterior_prob': pp, 'method': 'FUBAR'})
    return {'gene': gene, 'n_sites': len(sites), 'n_positive': len(hits), 'hits': hits}


def parse_relax_json(json_path, gene):
    """RELAX JSON → k + p-value"""
    d = json.load(open(json_path))
    tr = d.get('test results', {})
    return {'gene': gene, 'k': tr.get('relaxation or intensification parameter'),
            'p_value': tr.get('p-value'),
            'verdict': ('intensified' if (tr.get('relaxation or intensification parameter') or 0) > 1
                        else 'relaxed')}


def parse_contrastfel_json(json_path, gene):
    """CONTRASTFEL JSON → 位点级差异 (两组 β 差)"""
    d = json.load(open(json_path))
    mle = d['MLE']
    headers = [h[0] for h in mle['headers']]
    content = mle['content']
    sites = []
    for part in content.values():
        for row in part:
            sites.append(dict(zip(headers, row)))
    hits = []
    for si, s in enumerate(sites, 1):
        p = s.get('p-value')
        if isinstance(p, (int, float)) and p < 0.05:
            hits.append({'site': si, 'p_value': p, 'method': 'CONTRASTFEL'})
    return {'gene': gene, 'n_sites': len(sites), 'n_positive': len(hits), 'hits': hits}


def cmd_hyphy(args):
    """fubar / relax / contrastfel 子命令"""
    aligned_dir = Path(args.aligned_dir)
    tree_path = args.tree
    if tree_path != 'auto' and not os.path.exists(tree_path):
        print(f'[错误] 树不存在: {tree_path}')
        return 1
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 分组映射 (relax/contrastfel 需要, 循环内逐基因标注树)
    group_map = None
    test_groups = None
    if args.method in ('relax', 'contrastfel') and args.group_csv:
        group_map = load_group_map(args.group_csv, args.seq_col, args.group_col)
        test_groups = args.test_groups.split(',') if args.test_groups else None
        if not group_map:
            print('[警告] 分组映射为空')

    # 逐基因执行 (支持 cds_prep 输出 part_G.fas + 树 part_G.treefile)
    results = {}
    gene_map = load_genes(args.genes)
    fas_files = sorted(aligned_dir.glob('part_*.fas')) + \
                sorted(aligned_dir.glob('*-aligned.fasta')) + \
                sorted(aligned_dir.glob('*.part_*.fasta'))
    seen_files = set()
    for fas in fas_files:
        if str(fas) in seen_files:
            continue
        seen_files.add(str(fas))
        gene = None
        for g in ['N', 'P', 'P4', 'M', 'G', 'L']:
            if f'.part_{g}.' in fas.name or f'part_{g}.' in fas.name:
                gene = g
                break
        if not gene:
            gene = fas.stem.split('.')[-1]
        # 树: 同前缀 .treefile (若存在)
        tree_use = tree_path
        pref = str(fas)[:-4] if fas.suffix == '.fas' else str(fas).rsplit('.', 1)[0]
        cand_tree = f'{pref}.treefile'
        if args.tree == 'auto' and os.path.exists(cand_tree):
            tree_use = cand_tree
        elif args.tree == 'auto':
            print(f'[跳过] {gene}: 无树文件 ({cand_tree})')
            continue
        # 分组标注树 (逐基因)
        test_re = None
        if group_map and test_groups:
            tagged_tree = str(outdir / f'tagged_{gene}.nwk')
            tree_use, test_re, n_renamed = build_group_tree(
                tree_use, group_map, tagged_tree, test_groups)
            print(f'[分组] {gene}: 标注 {n_renamed} 叶, test={test_re}')
            if n_renamed == 0:
                print('  [警告] 叶标注失败 — 检查序列名与分组 CSV')
                from Bio import Phylo
                t = Phylo.read(tree_use, 'newick')
                print('  示例叶:', [x.name for x in t.get_terminals()[:3]])
        out_json = run_hyphy_method(args.method, str(fas), tree_use,
                                    str(outdir), test_re)
        if not out_json:
            continue
        if args.method == 'fubar':
            res = parse_fubar(out_json, gene)
            gene_to_genome(res['hits'], gene, gene_map)
        elif args.method == 'relax':
            res = parse_relax_json(out_json, gene)
        else:
            res = parse_contrastfel_json(out_json, gene)
            gene_to_genome(res['hits'], gene, gene_map)
        results[gene] = res
        # 打印
        if args.method == 'fubar':
            print(f"{gene}: {res['n_positive']}/{res['n_sites']} 正选择位点 "
                  f"(pp≥0.9)" + (f" → {[h['site'] for h in res['hits'][:15]]}" if res['hits'] else ''))
        elif args.method == 'relax':
            print(f"{gene}: k={res.get('k')} p={res.get('p_value')} ({res.get('verdict')})")
        else:
            print(f"{gene}: {res['n_positive']} 位点差异"
                  + (f" → {[h['site'] for h in res['hits'][:15]]}" if res['hits'] else ''))

    with open(outdir / f'{args.method}_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f'\n结果: {outdir}/{args.method}_results.json')
    return 0


# ============================================================
# 子命令: parse (解析已有 CAPHEINE 结果)
# ============================================================

def parse_busted(path, gene):
    d = json.load(open(path))
    tr = d.get('test results', {})
    out = {'gene': gene, 'lrt': tr.get('LRT'), 'p_value': tr.get('p-value')}
    try:
        fits = d['fits']
        unconstrained = fits.get('MG94xREV with separate rates for branch sets', {})
        rd = unconstrained.get('Rate Distributions', {})
        bg = rd.get('non-synonymous/synonymous rate ratio for *background*')
        test = rd.get('non-synonymous/synonymous rate ratio for *test*')
        out['omega_background'] = bg[0][0] if bg else None
        out['omega_test'] = test[0][0] if test else None
    except Exception:
        pass
    out['verdict'] = ('positive' if (out['p_value'] is not None and out['p_value'] < 0.05)
                      else 'no-signal')
    return out


def parse_fel_meme(path, gene, fdr=0.1, method='FEL'):
    d = json.load(open(path))
    mle = d['MLE']
    headers = [h[0] for h in mle['headers']]
    content = mle['content']
    sites = []
    for part in content.values():
        for row in part:
            sites.append(dict(zip(headers, row)))
    pvals = [s['p-value'] for s in sites
             if isinstance(s.get('p-value'), (int, float))]
    qs = bh_fdr(pvals) if pvals else []
    qi = 0
    hits = []
    for si, s in enumerate(sites, 1):
        p = s.get('p-value')
        if not isinstance(p, (int, float)):
            continue
        q = qs[qi]; qi += 1
        if method == 'FEL':
            positive = (s.get('beta', 0) > s.get('alpha', 0)) and p < 0.05
        else:
            positive = p < 0.05
        if positive:
            hits.append({'site': si, 'alpha': s.get('alpha'),
                         'beta': s.get('beta'), 'p_value': p, 'q_value': q,
                         'lrt': s.get('LRT')})
    return {'gene': gene, 'n_sites': len(sites), 'n_positive': len(hits),
            'hits': hits, 'method': method}


def cmd_parse(args):
    from collections import OrderedDict
    base = Path(args.capheine_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gene_map = load_genes(args.genes)
    results = OrderedDict()

    bd = (base / 'hyphy' / 'BUSTED') if (base / 'hyphy' / 'BUSTED').exists() else base / 'BUSTED'
    busted = []
    if bd.exists():
        for f in sorted(bd.glob('*.BUSTED.json')):
            gene = f.name.replace('.BUSTED.json', '').split('.part_')[-1]
            busted.append(parse_busted(f, gene))
    results['busted'] = busted

    fd = (base / 'hyphy' / 'FEL') if (base / 'hyphy' / 'FEL').exists() else base / 'FEL'
    fel = []
    if fd.exists():
        for f in sorted(fd.glob('*.FEL.json')):
            gene = f.name.replace('.FEL.json', '').split('.part_')[-1]
            fel.append(parse_fel_meme(f, gene, args.fdr, 'FEL'))
    results['fel'] = fel

    md = (base / 'hyphy' / 'MEME') if (base / 'hyphy' / 'MEME').exists() else base / 'MEME'
    meme = []
    if md.exists():
        for f in sorted(md.glob('*.MEME.json')):
            gene = f.name.replace('.MEME.json', '').split('.part_')[-1]
            meme.append(parse_fel_meme(f, gene, args.fdr, 'MEME'))
    results['meme'] = meme

    rd = (base / 'hyphy' / 'RELAX') if (base / 'hyphy' / 'RELAX').exists() else base / 'RELAX'
    relax = []
    if rd.exists():
        for f in sorted(rd.glob('*.RELAX.json')):
            gene = f.name.replace('.RELAX.json', '').split('.part_')[-1]
            relax.append(parse_relax_json(f, gene))
    results['relax'] = relax

    cd = (base / 'hyphy' / 'CONTRASTFEL') if (base / 'hyphy' / 'CONTRASTFEL').exists() else base / 'CONTRASTFEL'
    contrast = []
    if cd.exists():
        for f in sorted(cd.glob('*.CONTRASTFEL.json')):
            gene = f.name.replace('.CONTRASTFEL.json', '').split('.part_')[-1]
            contrast.append(parse_contrastfel_json(f, gene))
    results['contrastfel'] = contrast

    print('=== BUSTED (基因级正选择) ===')
    print(f"{'基因':<6}{'LRT':>10}{'p-value':>12}{'ω背景':>10}{'ω测试':>10}  判定")
    for b in busted:
        print(f"{b['gene']:<6}{str(b.get('lrt')):>10}{str(b.get('p_value')):>12}"
              f"{str(b.get('omega_background')):>10}{str(b.get('omega_test')):>10}  {b['verdict']}")

    print('\n=== FEL (位点级, p<0.05) ===')
    for f_ in fel:
        gstart = gene_map.get(f_['gene'], (None, None))[0]
        pos_str = ''
        if gstart:
            pos_str = ' → 基因组: ' + ','.join(str(gstart + (h['site'] - 1) * 3) for h in f_['hits'][:15])
        print(f"{f_['gene']}: {f_['n_positive']}/{f_['n_sites']} 正选择位点 "
              + (f" → site {[h['site'] for h in f_['hits'][:15]]}{pos_str}" if f_['hits'] else ''))

    print('\n=== MEME (episodic 位点) ===')
    for m_ in meme:
        gstart = gene_map.get(m_['gene'], (None, None))[0]
        pos_str = ''
        if gstart:
            pos_str = ' → 基因组: ' + ','.join(str(gstart + (h['site'] - 1) * 3) for h in m_['hits'][:10])
        if m_['hits']:
            print(f"{m_['gene']}: {m_['n_positive']} 位点 (site {[h['site'] for h in m_['hits'][:10]]}{pos_str})")
        else:
            print(f"{m_['gene']}: 0")

    print('\n=== RELAX ===')
    for r_ in relax:
        print(f"{r_['gene']}: k={r_.get('k')} p={r_.get('p_value')} {r_.get('verdict', '')}")

    print('\n=== CONTRASTFEL ===')
    if contrast:
        for c_ in contrast:
            print(f"{c_['gene']}: {json.dumps(c_['test_results'])[:200]}")
    else:
        print('无结果')

    with open(outdir / 'selection_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f'\n结果: {outdir}/selection_results.json')
    return 0


# ============================================================
# 主入口
# ============================================================

def main():
    ap = argparse.ArgumentParser(description='选择压力分析套件')
    sub = ap.add_subparsers(dest='cmd', required=True)

    # parse (已有 CAPHEINE 结果)
    p = sub.add_parser('parse', help='解析已有 CAPHEINE 结果')
    p.add_argument('--capheine-dir', required=True)
    p.add_argument('--genes', default=None)
    p.add_argument('--fdr', type=float, default=0.1)
    p.add_argument('--outdir', default='selection_out')

    # hyphy 方法
    p = sub.add_parser('run', help='运行 HyPhy 方法 (fubar/relax/contrastfel)')
    p.add_argument('--method', required=True, choices=['fubar', 'relax', 'contrastfel'])
    p.add_argument('--aligned-dir', required=True, help='CDS 分区比对目录')
    p.add_argument('--tree', default='auto', help='树 (newick/nhx); auto=用 cds_prep 每基因树')
    p.add_argument('--genes', default=None, help='基因坐标 CSV')
    p.add_argument('--group-csv', default=None, help='分组元数据 CSV')
    p.add_argument('--seq-col', default='name')
    p.add_argument('--group-col', default='host')
    p.add_argument('--test-groups', default=None, help='测试组 (逗号分隔, 如 Lycium ruthenicum)')
    p.add_argument('--outdir', default='selection_out')

    # 占位 (后续实现)
    for name in ['m7m8', 'time-tree', 'network', 'structure']:
        p = sub.add_parser(name, help=f'{name} (规划中)')
        p.add_argument('--input', required=True)
        p.add_argument('--outdir', default='selection_out')

    args = ap.parse_args()

    if args.cmd == 'parse':
        return cmd_parse(args)
    elif args.cmd == 'run':
        return cmd_hyphy(args)
    else:
        print(f'[{args.cmd}] 模块规划中 — 后续步骤实现')
        return 0


if __name__ == '__main__':
    sys.exit(main())
