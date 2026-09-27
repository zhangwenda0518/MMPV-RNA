#!/usr/bin/env python3
"""dnasp_bridge.py — DnaSP 6 (dnasp.py) 分析桥接

调用 dnasp.py (DnaSP 6 完整重实现) 补充:
  - Fu & Li's D*/F*   (中性检验缺口)
  - Rm 最小重组数      (重组检测)
  - Ka/Ks (Nei-Gojobori) (正选择整体 ω)
  - HKA / MK / Fay&Wu / SFS / LD / tstv 等

用法:
  from utils.dnasp_bridge import DnaSPAnalyzer
  d = DnaSPAnalyzer(fasta)
  d.full_report(out_dir)     # 全套 DnaSP 分析
  d.neutrality()             # Tajima/Fu&Li/Fu's Fs/R2
  d.recombination()          # Rm
  d.kaks()                   # Ka/Ks

CLI:
  python -m utils.dnasp_bridge --fasta X --out dir
"""
import argparse
import os
import sys
import json
from collections import Counter
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dnasp  # DnaSP 6 reimplementation


class DnaSPAnalyzer:
    """DnaSP 6 分析封装 (dnasp.py 引擎)"""

    def __init__(self, fasta_file, pop_file=None):
        self.aln = dnasp.load_alignment(Path(fasta_file))
        self.pop_assignments = None
        if pop_file:
            self.pop_assignments = dnasp.load_pop_file(Path(pop_file))

    def run(self, analyses, out_dir=None):
        """运行指定分析集合, 返回结果字典"""
        res = dnasp.run_analysis(
            self.aln,
            analyses=set(analyses),
            pop_assignments=self.pop_assignments,
        )
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
            import json
            # 序列化可 JSON 的结果
            def _conv(o):
                if hasattr(o, '__dict__'):
                    return {k: _conv(v) for k, v in vars(o).items()
                            if not k.startswith('_')}
                if isinstance(o, (list, tuple)):
                    return [_conv(x) for x in o]
                if isinstance(o, dict):
                    return {k: _conv(v) for k, v in o.items()}
                return o
            try:
                json.dump(_conv(res), open(f'{out_dir}/dnasp_results.json', 'w'),
                          indent=2, default=str)
            except Exception as e:
                print(f'[dnasp_bridge] 警告: dnasp_results.json 写入失败: {e}', file=sys.stderr)
        return res

    def kaks_by_gene(self, genes_csv, out_json=None, reference=None):
        """按基因 Ka/Ks (通用: 任意病毒 + 基因注释 CSV)

        genes_csv: 基因坐标 (gene,start,end[,type]) 参考坐标 1-based
        reference: 参考序列名 (坐标映射基准); None=第一条
        注意: 建议传入含参考序列的原始比对 (非无gap过滤版)
        """
        # 读基因坐标
        import csv as _csv
        genes = []
        with open(genes_csv) as f:
            for row in _csv.DictReader(f):
                genes.append((row['gene'], int(row['start']), int(row['end'])))
        seen = set()
        genes_u = []
        for g in genes:
            if g[0] not in seen:
                seen.add(g[0])
                genes_u.append(g)

        seqs = self._raw_seqs()
        # 选参考序列 (指定名或第一条)
        # 2026-09-15 修复: 旧版 `reference in nm` 是子串匹配 → `CRR1` 命中 `CRR10`;
        # 且显式指定的参考找不到时**静默回落**到第一条, 坐标基准悄悄改变。
        # 现走统一三级匹配, 显式指定找不到即报错。
        try:
            from utils.seq_match import resolve_reference as _resolve_ref
        except ImportError:                                  # 直接以顶层模块导入时
            from seq_match import resolve_reference as _resolve_ref
        ref_name, ref_seq = _resolve_ref(
            seqs, reference, who="dnasp_bridge.kaks_by_gene")

        # 比对列 → 参考坐标映射
        ref_pos = {}
        seen_n = 0
        for i, c in enumerate(ref_seq):
            if c != '-':
                seen_n += 1
                ref_pos[seen_n] = i

        results = []
        for g, start, end in genes_u:
            if start not in ref_pos or end not in ref_pos:
                continue
            cs, ce = ref_pos[start], ref_pos[end] + 1
            gene_seqs = [s[cs:ce] for _, s in seqs]
            n_codons_full = (ce - cs) // 3 * 3
            gene_seqs_c = [s[:n_codons_full] for s in gene_seqs]
            try:
                r = dnasp.compute_ka_ks(gene_seqs_c)
                if r.omega is not None:
                    verdict = ('positive' if r.omega > 1 else
                               ('purifying' if r.omega < 0.7 else 'near-neutral'))
                    results.append({
                        'gene': g, 'start': start, 'end': end,
                        # cd_length 报实际参与计算的密码子长度 (非 3 倍数尾部被截)
                        'cd_length': n_codons_full,
                        'n_codons': r.n_codons, 'S_sites': round(r.S_sites, 2),
                        'N_sites': round(r.N_sites, 2),
                        # 2026-09-15 (审查 P2-18) 修复: 旧版 `round(r.Ka,5) if r.Ka else None`
                        # —— Ka/Ks/omega 恰为 0.0 时被当成"缺失"写成 None, 下游
                        # 看到的是"没算出来"而不是"算出来等于 0"。改用 is not None。
                        'ka': round(r.Ka, 5) if r.Ka is not None else None,
                        'ks': round(r.Ks, 5) if r.Ks is not None else None,
                        'omega': round(r.omega, 4) if r.omega is not None else None,
                        'verdict': verdict,
                    })
            except Exception as e:
                # 单基因失败不再静默: 记录告警 + 失败基因, 便于排查
                print(f'[dnasp_bridge] 警告: {g} Ka/Ks 计算失败: {e}', file=sys.stderr)
                results.append({'gene': g, 'start': start, 'end': end,
                                'error': str(e)})
        if out_json:
            json.dump(results, open(out_json, 'w'), indent=2)
        return results

    def _raw_seqs(self):
        """从 dnasp Alignment 提取原始序列 [(id, seq)]"""
        return list(zip(self.aln.names, self.aln.seqs))

    def recombination_profile(self, genes_csv=None, out_json=None):
        """重组信号分析: Rm + 断点分布 + 基因映射 (通用)

        返回: {rm, n_pairs, breakpoints_by_gene, density_per_kb,
               hotspot_windows (100bp top)}
        """
        res = dnasp.run_analysis(self.aln, analyses={'recombination'})
        recomb = res['recombination']
        pairs = recomb.incompatible_pairs
        n_pairs = recomb.n_incompatible_pairs
        rm = recomb.Rm

        # 断点 = 不兼容位点对中点
        breakpoints = [(a + b) / 2 for a, b in pairs]

        # 100bp 窗口热点
        bins = Counter(int(bp // 100) * 100 for bp in breakpoints)
        hotspots = [{'window_start': w, 'count': c}
                    for w, c in sorted(bins.items(), key=lambda x: -x[1])[:10]]

        result = {'rm': rm, 'n_incompatible_pairs': n_pairs,
                  'hotspot_windows': hotspots}

        # 基因映射 (若提供基因注释)
        if genes_csv:
            genes = []
            with open(genes_csv) as f:
                for row in __import__('csv').DictReader(f):
                    # 历史坑: `or True` 恒真, type 白名单是死代码, 所有行都会被纳入
                    # (含 ncRNA/intergenic 等), 断点可能被错误归到非编码区。
                    if row.get('type', 'CDS') in ('CDS', 'gene', ''):
                        genes.append((row['gene'], int(row['start']), int(row['end'])))
            seen = set()
            genes_u = []
            for g in genes:
                if g[0] not in seen:
                    seen.add(g[0])
                    genes_u.append(g)
            gene_hits = Counter()
            for a, b in pairs:
                mid = (a + b) / 2
                for g, start, end in genes_u:
                    if start <= mid <= end:
                        gene_hits[g] += 1
                        break
                else:
                    gene_hits['intergenic'] += 1
            total = max(sum(gene_hits.values()), 1)
            result['breakpoints_by_gene'] = {
                g: {'count': c, 'pct': round(c / total * 100, 1)}
                for g, c in gene_hits.most_common()
            }
            result['density_per_kb'] = {
                g: round(c / (end - start + 1) * 1000, 1)
                for g, start, end in genes_u
                for g2, c in gene_hits.items() if g2 == g
            }

        if out_json:
            json.dump(result, open(out_json, 'w'), indent=2, default=str)
        return result

    def neutrality(self, out_dir=None):
        """中性检验: Tajima's D + Fu & Li's D*/F* + Fu's Fs + R2"""
        return self.run(['polymorphism', 'fufs'], out_dir)

    def recombination(self, out_dir=None):
        """重组: Rm 最小重组数"""
        return self.run(['recombination'], out_dir)

    def kaks(self, out_dir=None):
        """Ka/Ks (Nei-Gojobori)"""
        return self.run(['kaks'], out_dir)

    def ld(self, out_dir=None):
        """连锁不平衡: D/D'/R²/ZnS"""
        return self.run(['ld'], out_dir)

    def divergence(self, fasta2, out_dir=None):
        """群体间分化: Dxy/Da (需要第二群体)"""
        aln2 = dnasp.load_alignment(Path(fasta2))
        return dnasp.run_analysis(self.aln, aln2=aln2, analyses={'divergence'})

    def full_report(self, out_dir):
        """全套 DnaSP 分析 + TSV + Markdown 报告 + 图"""
        analyses = {'polymorphism', 'ld', 'recombination', 'popsize',
                    'indel', 'kaks', 'fufs', 'sfs', 'tstv', 'codon', 'faywu'}
        res = self.run(analyses, out_dir)
        dnasp.make_figures(Path(out_dir), res)
        dnasp.write_reproducibility(Path(out_dir), res, analyses)
        return res


def main():
    ap = argparse.ArgumentParser(description='DnaSP 6 (dnasp.py) 分析')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--pop-file', default=None, help='群体分组文件 (seq→pop)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--analyses', default='polymorphism,fufs,recombination,kaks,ld',
                    help='分析集合 (逗号分隔)')
    ap.add_argument('--cli-report', action='store_true',
                    help='用 dnasp.py 原 CLI 输出完整报告 (TSV/Markdown/图)')
    ap.add_argument('--genes', default=None, help='基因坐标 CSV (gene,start,end) 用于按基因 Ka/Ks + 断点映射')
    ap.add_argument('--kaks-by-gene', action='store_true', help='按基因算 Ka/Ks')
    ap.add_argument('--reference', default=None, help='参考序列名 (坐标映射基准)')
    ap.add_argument('--recomb-profile', action='store_true', help='重组断点分析')
    args = ap.parse_args()

    if args.cli_report:
        # 用 dnasp.py 原生 CLI 出报告
        # 历史坑: 硬编码 'python3' 在 Windows 无 python3.exe → FileNotFoundError 被上层
        # 吞成 warning，DnaSP 补充检验静默无产出。统一 sys.executable。
        cmd = [sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dnasp.py'),
               '--input', args.fasta, '--analysis', args.analyses, '--output', args.out]
        if args.pop_file:
            cmd += ['--pop-file', args.pop_file]
        import subprocess
        # 2026-09-15: 补退码检查 —— 原来既不查 returncode 也不读 stderr,
        # dnasp 失败时静默 return, 调用方以为分析成功。
        _r = subprocess.run(cmd, capture_output=True, text=True,
                            encoding='utf-8', errors='replace')
        if _r.returncode != 0:
            raise SystemExit(f"[dnasp] 子进程退出码 {_r.returncode}: "
                             f"{(_r.stderr or '')[-500:]}")
        return

    d = DnaSPAnalyzer(args.fasta, args.pop_file)

    # 按基因 Ka/Ks (通用)
    if args.kaks_by_gene:
        if not args.genes:
            print('需 --genes 基因坐标 CSV')
            return 1
        try:
            kaks = d.kaks_by_gene(args.genes, out_json=f'{args.out}/kaks_by_gene.json',
                                  reference=args.reference)
        except ValueError as e:
            # 显式指定的参考找不到时拒绝静默换基准 (2026-09-15)
            print(f'[错误] {e}', file=sys.stderr)
            return 1
        print('=== 按基因 Ka/Ks ===')
        for r in kaks:
            # 2026-09-15 修复: 单基因失败时 kaks_by_gene 会追加 {'error': ...} 条目,
            # 旧版无条件读 r['omega'] → KeyError 直接崩掉整个 CLI。
            if r.get('error'):
                print(f"{r.get('gene', '?')}: 计算失败 ({r['error']})")
                continue
            print(f"{r['gene']}: ω={r['omega']} ({r['verdict']}) Ka={r['ka']} Ks={r['ks']}")
        return

    # 重组断点分析 (通用)
    if args.recomb_profile:
        rp = d.recombination_profile(args.genes, out_json=f'{args.out}/recombination_profile.json')
        print(f"Rm = {rp['rm']}, 不兼容位点对 = {rp['n_incompatible_pairs']}")
        if 'breakpoints_by_gene' in rp:
            print('断点按基因:', rp['breakpoints_by_gene'])
            print('断点密度/kb:', rp['density_per_kb'])
        print('断点热点窗口:', rp['hotspot_windows'][:5])
        return

    analyses = set(a.strip() for a in args.analyses.split(','))
    res = d.run(analyses, args.out)
    # 简要输出关键统计
    for key in ['polymorphism', 'recombination', 'kaks', 'fufs']:
        if key in res:
            print(f"[{key}] {str(res[key])[:200]}")


if __name__ == '__main__':
    main()
