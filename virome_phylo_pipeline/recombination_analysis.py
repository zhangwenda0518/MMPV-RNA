#!/usr/bin/env python3
"""recombination_analysis.py — 重组分析统一入口（RDP5 全流程）

串联: RDP5 检测 → IQ-TREE/SimPlot 验证 → 断点分布图 → dnasp Rm 交叉验证

用法:
  python -m recombination_analysis \
      --fasta alignment.fas --prefix my_virus --outdir out \
      --genes genes.csv [--rdp5-script ~/MMPV-RNA/biosoft/rdp5/run_rdp5.sh]

阶段 (可跳过):
  run       RDP5 检测 (run_rdp5.sh, 输出 <prefix>.csv / RecIDTests / 断点分布)
  validate  rdp5_validate.py (IQ-TREE 拓扑切换 + SimPlot + 树图)
  plot      rdp5_bdp_plot.py (断点分布图, 需 R + RDP5_RBDP_Rgrapher)
  cross     dnasp Rm 断点 vs RDP5 断点交叉验证 (仅统计/报告)

输出:
  <outdir>/recombination_report.md  完整报告
  <outdir>/rdp5_events.json         事件表 (含基因映射/方法p值/可靠性)
"""
import argparse, json, os, subprocess, sys
from pathlib import Path

UTILS = Path(__file__).parent


def run_rdp5(rdp5_script, fasta, prefix, outdir):
    """RDP5 检测: bash run_rdp5.sh <fasta> <prefix> <outdir>"""
    os.makedirs(outdir, exist_ok=True)
    print(f"[run] RDP5: {rdp5_script} {fasta} {prefix} {outdir}")
    try:
        r = subprocess.run(['bash', str(rdp5_script), str(fasta), prefix, str(outdir)],
                           capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=3600)
    except subprocess.TimeoutExpired:
        print('[警告] RDP5 超时 (>3600s)')
        return None
    print(r.stdout[-2000:] if r.stdout else '')
    if r.returncode != 0:
        print(f'[警告] RDP5 退出码 {r.returncode}: {r.stderr[-800:]}')
        return None
    # RDP5 输出文件 (主 CSV + RecIDTests + 断点分布)
    csv = Path(outdir) / f'{prefix}.csv'
    return csv if csv.exists() else None


def validate(csv_path, fasta, outdir, skip_trees=False, threads='8'):
    """验证: rdp5_validate.py"""
    print(f"[validate] 事件解析 + IQ-TREE 拓扑验证 + SimPlot")
    cmd = [sys.executable, str(UTILS / 'rdp5_validate.py'),
           '-r', str(csv_path), '-a', str(fasta), '-o', str(outdir),
           '-t', threads]
    if skip_trees:
        cmd.append('--skip-trees')
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=7200)
    except subprocess.TimeoutExpired:
        print('[警告] 验证超时 (>7200s)')
        return None
    print(r.stdout[-2500:] if r.stdout else '')
    if r.returncode != 0:
        print(f'[警告] 验证退出码 {r.returncode}: {r.stderr[-800:]}')
    return outdir


def plot(bdp_csv, genes_csv, outdir, name):
    """断点分布图: rdp5_bdp_plot.py
    需要 -b/-p/-g 三个 required 参数 + -n(标题) + -o(PDF)。
    历史坑: 旧调用用 -t (非法参数) 且缺 -p, 必然 argparse 失败, 管线却只警告继续。
    """
    bdp = Path(bdp_csv)
    if not bdp.exists():
        print(f'[plot] 跳过: 无断点分布 CSV {bdp}')
        return None
    print(f"[plot] 断点分布图: {bdp.name}")
    # -p positions: 断点位置 CSV (RDP5 *BreakpointPositions.csv)
    pos_csv = None
    for pat in ('*BreakpointPositions*.csv', '*breakpoint*position*.csv', '*positions*.csv'):
        cands = list(bdp.parent.glob(pat))
        if cands:
            pos_csv = cands[0]
            break
    if pos_csv is None:
        print(f'[plot] 跳过: 未找到断点位置 CSV (BreakpointPositions.csv), rdp5_bdp_plot 需要 -p')
        return None
    # -g gene-map: ORF 坐标 CSV (*ORFCoords.csv, RDP5 产物)
    orf_csv = None
    for pat in ('*ORFCoords*.csv', '*orf*coords*.csv'):
        cands = list(bdp.parent.glob(pat))
        if cands:
            orf_csv = cands[0]
            break
    if orf_csv is None:
        print(f'[plot] 跳过: 未找到 ORF 坐标 CSV (ORFCoords.csv), rdp5_bdp_plot 需要 -g')
        return None
    out_pdf = Path(outdir) / f'{name}_breakpoint_plot.pdf'
    cmd = [sys.executable, str(UTILS / 'rdp5_bdp_plot.py'),
           '-b', str(bdp), '-p', str(pos_csv), '-g', str(orf_csv),
           '-n', f'{name} RDP5 Breakpoints', '-o', str(out_pdf)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=600)
    except subprocess.TimeoutExpired:
        print('[警告] 断点图超时 (>600s)')
        return None
    print(r.stdout[-800:] if r.stdout else '')
    if r.returncode != 0:
        print(f'[警告] 断点图退出码 {r.returncode}: {r.stderr[-500:]}')
        return None
    return str(out_pdf) if out_pdf.exists() else None


def parse_events(csv_path):
    """复用 rdp5_validate 的解析器提取事件摘要"""
    sys.path.insert(0, str(UTILS))
    try:
        from utils.rdp5_validate import parse_rdp5_csv
        return parse_rdp5_csv(csv_path)
    except Exception as e:
        print(f'[cross] 事件解析失败: {e}')
        return []


def dnasp_rm_cross(fasta, events, genes_csv):
    """dnasp Rm 与 RDP5 断点交叉 (统计性)

    2026-09-15 (审查 P2-24) 修复: genes_csv 参数旧版完全未使用 (调用方传了
    基因坐标表却没有任何效果)。现在用它把重合窗口注释到基因名 —— 报告里能直接
    看出"重组热点落在哪个基因", 而不是只有一串裸窗口坐标。
    """
    sys.path.insert(0, str(UTILS))
    try:
        import dnasp
        from pathlib import Path as P
        from collections import Counter
        aln = dnasp.load_alignment(P(fasta))
        res = dnasp.run_analysis(aln, analyses={'recombination'})
        rm = res['recombination'].Rm
        n_pairs = res['recombination'].n_incompatible_pairs
        pairs = res['recombination'].incompatible_pairs
        rm_bins = Counter(int((a + b) / 2 // 100) * 100 for a, b in pairs)
        rm_hot = [w for w, c in rm_bins.most_common(5)]
        # RDP5 断点
        rdp_hot = []
        for ev in events:
            rdp_hot.append((ev['bp_start'] // 100) * 100)
            rdp_hot.append((ev['bp_end'] // 100) * 100)
        overlap = sorted(set(rm_hot) & set(rdp_hot))

        # 重合窗口 → 基因注释 (genes_csv: gene,start,end; 1-based 含端)
        overlap_genes = {}
        genes = []
        if genes_csv and os.path.exists(str(genes_csv)):
            try:
                import csv as _csv
                with open(genes_csv, encoding='utf-8-sig') as _f:
                    for row in _csv.DictReader(_f):
                        try:
                            genes.append((str(row.get('gene') or row.get('name') or '?'),
                                          int(row['start']), int(row['end'])))
                        except (KeyError, TypeError, ValueError):
                            continue
            except OSError as _e:
                print(f'[cross] genes_csv 读取失败: {_e}')
        for w in overlap:
            hit = [g for g, s, e in genes if not (e < w or s > w + 99)]
            overlap_genes[w] = hit or ['(intergenic)']

        return {'rm': rm, 'n_incompatible_pairs': n_pairs,
                'rm_hotspot_windows': rm_hot,
                'rdp5_breakpoint_windows': sorted(set(rdp_hot)),
                'overlap_windows': overlap,
                'overlap_genes': overlap_genes}
    except Exception as e:
        return {'error': str(e)}


def main():
    ap = argparse.ArgumentParser(description='RDP5 重组分析统一入口')
    ap.add_argument('--fasta', required=True, help='比对 FASTA')
    ap.add_argument('--prefix', required=True, help='输出前缀 (样本/病毒名)')
    ap.add_argument('--outdir', default='recombination_out', help='输出目录')
    ap.add_argument('--genes', default=None, help='基因注释 CSV (gene,start,end)')
    ap.add_argument('--rdp5-script', default=None,
                    help='run_rdp5.sh 路径 (默认无, 需显式提供或已有结果)')
    ap.add_argument('--rdp5-csv', default=None,
                    help='已有 RDP5 结果 CSV (跳过 run, 直接解析)')
    ap.add_argument('--stage', default='all',
                    choices=['all', 'run', 'validate', 'plot', 'cross'],
                    help='执行阶段')
    ap.add_argument('--skip-trees', action='store_true', help='验证跳过 IQ-TREE')
    ap.add_argument('--threads', default='8', help='IQ-TREE 线程')
    args = ap.parse_args()

    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix
    csv_path = Path(args.rdp5_csv).expanduser() if args.rdp5_csv else outdir / f'{prefix}.csv'
    bdp_csv = outdir / f'{prefix} breakpoint distribution.csv'
    stages = [args.stage] if args.stage != 'all' else ['run', 'validate', 'plot', 'cross']

    if 'run' in stages:
        rdp5 = os.path.expanduser(args.rdp5_script) if args.rdp5_script else None
        if not rdp5 or not os.path.exists(rdp5):
            print(f'[run] RDP5 脚本不存在/未提供: {rdp5} — 若已有结果可用 --stage validate/plot/cross')
        else:
            csv_path = run_rdp5(rdp5, args.fasta, prefix, outdir) or csv_path

    # 事件解析 (csv 存在即解析, 供 validate/cross 共用)
    events = []
    if csv_path.exists():
        events = parse_events(csv_path)
        print(f'[parse] 解析到 {len(events)} 个重组事件')
    elif 'validate' in stages or 'cross' in stages:
        print(f'[parse] 无 {csv_path} — 先运行 RDP5 (或提供已有结果)')

    if 'validate' in stages:
        if csv_path.exists():
            validate(csv_path, args.fasta, outdir,
                     skip_trees=args.skip_trees, threads=args.threads)
        else:
            print(f'[validate] 跳过: 无 {csv_path} (先跑 RDP5)')

    if 'plot' in stages:
        plot(bdp_csv, args.genes, outdir, prefix)

    report = [f'# 重组分析报告: {prefix}\n',
              f'## RDP5 事件 ({len(events)})']
    if events:
        for ev in events:
            # 历史坑: methods 为空 dict 时 min() 崩 (所有方法 p 值 NS/空)
            best_p = min(ev['methods'].values()) if ev.get('methods') else float('nan')
            report.append(
                f"- **Event {ev['num']}** {ev['bp_start']}-{ev['bp_end']}: "
                f"重组子 {ev['recombinant']} ({ev['minor_parent']}/{ev['major_parent']}), "
                f"{ev['n_methods']} 方法 (best p={best_p:.2e})")
    else:
        report.append('- 无事件或未解析')

    if 'cross' in stages:
        cross = dnasp_rm_cross(args.fasta, events, args.genes)
        report.append(f"\n## dnasp Rm 交叉验证")
        if 'error' in cross:
            report.append(f"- {cross['error']}")
        else:
            report.append(f"- dnasp Rm = {cross['rm']} ({cross['n_incompatible_pairs']} 不兼容对)")
            report.append(f"- dnasp 断点热点窗口: {cross['rm_hotspot_windows']}")
            report.append(f"- RDP5 断点窗口: {cross['rdp5_breakpoint_windows']}")
            report.append(f"- **重合窗口: {cross['overlap_windows']}**")
            _og = cross.get('overlap_genes') or {}
            if _og:
                report.append("- 重合窗口 → 基因:")
                for _w, _gs in sorted(_og.items()):
                    report.append(f"  - {_w}-{_w + 99}: {', '.join(_gs)}")
        with open(outdir / 'rdp5_events.json', 'w') as f:
            json.dump([{'num': e['num'], 'bp_start': e['bp_start'], 'bp_end': e['bp_end'],
                        'recombinant': e['recombinant'], 'minor_parent': e['minor_parent'],
                        'major_parent': e['major_parent'], 'methods': e['methods']}
                       for e in events], f, indent=2, default=str)

    (outdir / 'recombination_report.md').write_text('\n'.join(report))
    print(f'\n报告: {outdir}/recombination_report.md')
    print('\n'.join(report))
    return 0


if __name__ == '__main__':
    sys.exit(main())
