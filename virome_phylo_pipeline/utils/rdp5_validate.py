#!/usr/bin/env python3
"""
rdp5_validate.py — Universal RDP5 recombination validator
=========================================================
一键完成：事件解析 → IQ-TREE ML树 → 亲本切换验证 → SimPlot → 树图
"""

import argparse, os, re, subprocess, shutil, sys
from pathlib import Path
from Bio import SeqIO, Phylo

def parse_args():
    p = argparse.ArgumentParser(description="RDP5 recombination validation — one click pipeline")
    p.add_argument("-r", "--rdp5-csv", required=True, help="RDP5 output CSV")
    p.add_argument("-a", "--alignment", required=True, help="Aligned FASTA")
    p.add_argument("-o", "--outdir", default="./rdp5_validation", help="Output dir")
    p.add_argument("--iqtree", default=None, help="Path to iqtree2 (auto-detect)")
    p.add_argument("-t", "--threads", default="1", help="IQ-TREE threads")
    p.add_argument("--quiet", action="store_true", help="Less output")
    p.add_argument("--skip-trees", action="store_true", help="Skip IQ-TREE")
    p.add_argument("--skip-simplots", action="store_true", help="Skip SimPlot graphs")
    p.add_argument("--skip-treefigs", action="store_true", help="Skip tree figures")
    p.add_argument("--clean", action="store_true", help="Remove intermediate IQ-TREE files")
    return p.parse_args()

# ── IQ-TREE auto-detect ──
def find_iqtree():
    for cmd in ["iqtree2", "iqtree"]:
        if shutil.which(cmd): return shutil.which(cmd)
    for p in Path(__file__).parent.glob("iqtree/**/iqtree2.exe"):
        return str(p)
    sys.exit("IQ-TREE not found. Use --iqtree or: conda install -c bioconda iqtree")

# ── RDP5 CSV parser ──
def parse_rdp5_csv(csv_path):
    events = []
    with open(csv_path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('Table key','~','*','$','^','Minor','Major','Unknown','NS')):
                continue
            m = re.match(r'\s*(\d+)\s*,\s*(\d+)([~\^\*\$]*)\s*,\s*(\d+)\s*,\s*(\d+)', line)
            if not m: continue
            fields = [x.strip() for x in line.split(',')]
            if len(fields) <= 10: continue
            
            rec_id = fields[8].split('\n')[0].strip().rstrip(',').lstrip('^')
            def get_parent(field):
                if field.startswith('Unknown'):
                    m2 = re.search(r'\(([^)]+)\)', field); return m2.group(1) if m2 else ''
                return field.split('\n')[0].strip().rstrip(',')
            
            minor_id = get_parent(fields[9])
            major_id = get_parent(fields[10])
            
            method_names = ['RDP','GENECONV','Bootscan','Maxchi','Chimaera','SiScan','PhylPro','LARD','3Seq']
            pvals = {}
            for i, mn in enumerate(method_names):
                idx = 11 + i
                if idx < len(fields) and fields[idx] and fields[idx] != 'NS':
                    try: pvals[mn] = float(fields[idx])
                    except ValueError: pass
            
            events.append({
                'num': int(m.group(1)), 'markers': m.group(3),
                'bp_start': int(m.group(4)), 'bp_end': int(m.group(5)),
                'recombinant': rec_id, 'minor_parent': minor_id, 'major_parent': major_id,
                'methods': pvals, 'n_methods': len(pvals)
            })
    return events

# ── Fuzzy sequence finder ──
def find_seq(seq_id, records):
    # 历史坑: 子串包含判断会把 CRR1 匹配到 CRR10; 先精确匹配, 再 token, 最后才子串
    for rec in records:
        if rec.id == seq_id:
            return rec
    toks = seq_id.split('_')
    for rec in records:
        if toks[0] in rec.id and toks[-1] in rec.id:
            return rec
    for rec in records:
        if seq_id in rec.id:
            return rec
    return None

# ── Region utils ──
def split_regions(triplet, bp_start, bp_end):
    aln_len = len(triplet[0].seq)
    bp_start = max(1, min(bp_start, aln_len))
    bp_end = max(1, min(bp_end, aln_len))
    if bp_start >= bp_end:
        a = [str(s.seq)[bp_end:bp_start-1] for s in triplet]
        b = [str(s.seq)[bp_start-1:] + str(s.seq)[:bp_end] for s in triplet]
        return a, b, f"nonrecomb_{bp_end+1}_{bp_start-1}", f"recomb_{bp_start}_{aln_len}+1_{bp_end}"
    else:
        a = [str(s.seq)[:bp_start-1] + str(s.seq)[bp_end:] for s in triplet]
        b = [str(s.seq)[bp_start-1:bp_end] for s in triplet]
        return a, b, f"nonrecomb_1_{bp_start-1}+{bp_end+1}_{aln_len}", f"recomb_{bp_start}_{bp_end}"

def valid_region(seqs):
    return any(sum(1 for c in s if c not in '-Nn?') > 0 for s in seqs)

# ── SimPlot generator ──
def make_simplot(rec, minor, major, bp, parents, outpath):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    
    rs = np.array(list(str(rec.seq).upper()))
    ms = np.array(list(str(minor.seq).upper())) if minor else None
    js = np.array(list(str(major.seq).upper())) if major else None
    
    win, step = 200, 20
    nw = (len(rs) - win) // step
    sim_minor, sim_major, pos = [], [], []
    for i in range(nw):
        s, e = i*step, i*step+win
        v = (rs[s:e]!='-')
        if ms is not None: v &= (ms[s:e]!='-')
        if js is not None: v &= (js[s:e]!='-')
        if v.sum() > 0:
            if minor: sim_minor.append((rs[s:e][v]==ms[s:e][v]).sum()/v.sum()*100)
            if major: sim_major.append((rs[s:e][v]==js[s:e][v]).sum()/v.sum()*100)
        else:
            if minor: sim_minor.append(np.nan)
            if major: sim_major.append(np.nan)
        pos.append(s+win//2)
    
    fig, ax = plt.subplots(figsize=(12, 4))
    minor_label = f"vs Minor ({parents[0][1].id[:15]}...)" if parents else "Minor"
    major_label = f"vs Major ({parents[1][1].id[:15]}...)" if len(parents)>1 else "Major"
    if minor: ax.plot(pos, sim_minor, 'r-', label=minor_label, lw=1.5)
    if major: ax.plot(pos, sim_major, 'b-', label=major_label, lw=1.5)
    ax.axvspan(bp[0], bp[1], alpha=0.15, color='orange', label=f'BP {bp[0]}-{bp[1]}')
    ax.set_xlabel('Genome Position (bp)'); ax.set_ylabel('Similarity (%)')
    ax.set_title(f'SimPlot: {rec.id[:30]}...')
    ax.legend(loc='lower right'); ax.set_ylim(50, 105); ax.grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(str(outpath), dpi=150); plt.close()

# ── Tree figure generator ──
def make_tree_figure(tree_a_path, tree_b_path, rec_id, minor_id, major_id, bp, outpath):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    
    COLORS = {'recombinant': '#E74C3C', 'major': '#27AE60', 'minor': '#3498DB', 'outgroup': '#95A5A6'}
    
    def leaf_color(label):
        if rec_id in (label or ''): return COLORS['recombinant']
        if major_id in (label or ''): return COLORS['major']
        if minor_id and minor_id in (label or ''): return COLORS['minor']
        return COLORS['outgroup']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for idx, (tp, title, ax) in enumerate([
        (tree_a_path, f'Non-recombinant\n(1-{bp[0]-1} + {bp[1]+1}-end)', axes[0]),
        (tree_b_path, f'Recombinant ({bp[0]}-{bp[1]})', axes[1])
    ]):
        if not os.path.exists(tp): continue
        tree = Phylo.read(tp, 'newick'); tree.ladderize()
        for leaf in tree.get_terminals(): leaf.color = leaf_color(leaf.name)
        Phylo.draw(tree, axes=ax, do_show=False,
                   label_func=lambda n: (n.name or '')[:25],
                   label_colors=lambda n: getattr(n, 'color', 'black'))
        ax.set_title(title, fontsize=11, fontweight='bold'); ax.axis('off')
    
    legend_els = [
        Rectangle((0,0),1,1, facecolor=COLORS['recombinant'], label='Recombinant'),
        Rectangle((0,0),1,1, facecolor=COLORS['major'], label='Major Parent'),
        Rectangle((0,0),1,1, facecolor=COLORS['minor'], label='Minor Parent'),
        Rectangle((0,0),1,1, facecolor=COLORS['outgroup'], label='Reference'),
    ]
    fig.legend(handles=legend_els, loc='lower center', ncol=4, fontsize=9, bbox_to_anchor=(0.5,-0.02))
    fig.suptitle(f'{rec_id[:35]}... — ML Tree Validation', fontsize=13, fontweight='bold')
    plt.tight_layout(rect=[0,0.05,1,0.93])
    plt.savefig(str(outpath), dpi=200, bbox_inches='tight'); plt.close()

# ── MAIN ──
def main():
    args = parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    iqtree = args.iqtree or find_iqtree()
    
    events = parse_rdp5_csv(args.rdp5_csv)
    if not events:
        print("No events found."); return
    
    records = list(SeqIO.parse(args.alignment, "fasta"))
    print(f"Alignment: {len(records)} seqs  |  Events: {len(events)}  |  IQ-TREE: {iqtree}")
    
    report = [f"# RDP5 Recombination Validation Report\n"
              f"**Input**: {args.rdp5_csv}, {len(records)} sequences\n"]
    
    passed = failed = 0
    for evt in events:
        n = evt['num']; d = outdir / f"event_{n:02d}"; d.mkdir(exist_ok=True)
        rec_id = evt['recombinant']; minor_id = evt['minor_parent']; major_id = evt['major_parent']
        bp = (evt['bp_start'], evt['bp_end'])
        
        print(f"\nEvent {n}: {rec_id}  BP={bp}  {'~' if evt['markers'] else ''}")
        
        rec = find_seq(rec_id, records)
        if not rec: print("  ⚠ Rec not found"); continue
        parents = [(t, s) for t, sid in [('major',major_id),('minor',minor_id)] if sid and (s:=find_seq(sid,records))]
        if not parents: print("  ⚠ No parents"); continue
        
        triplet = [rec] + [p[1] for p in parents]
        
        # ── IQ-TREE ──
        topology = "UNDETERMINED"
        if not args.skip_trees:
            ra, rb, ral, rbl = split_regions(triplet, bp[0], bp[1])
            for rname, rseqs, rlabel in [('region_a',ra,ral),('region_b',rb,rbl)]:
                if not valid_region(rseqs): continue
                fa = d / f"{rname}.fasta"
                with open(fa, 'w') as f:
                    f.write(f">{rec_id}\n{rseqs[0]}\n")
                    for j, (pt, ps) in enumerate(parents):
                        f.write(f">{ps.id}[{pt}]\n{rseqs[j+1]}\n")
                
                cmd = [iqtree, "-s", str(fa), "-m", "GTR+G", "-nt", args.threads,
                       "--prefix", str(d/rname), "-redo", "-quiet"]
                try:
                    r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=300)
                except subprocess.TimeoutExpired:
                    print(f'[validate] 警告: IQ-TREE 超时 ({rname}), 跳过该区域', file=sys.stderr)
                    continue
                tf = str(d/rname) + ".treefile"
                if r.returncode == 0 and os.path.exists(tf):
                    shutil.copy(tf, d / f"{rname}_{rlabel}.nwk")
                    if not args.quiet: print(f"  ✓ {rname} tree")
            
            ta = d / "region_a.treefile"; tb = d / "region_b.treefile"
            if os.path.exists(ta) and os.path.exists(tb):
                try:
                    t_a = Phylo.read(ta, "newick"); t_b = Phylo.read(tb, "newick")
                    rl = ml = nil = None
                    for lf in t_a.get_terminals():
                        n = lf.name or ''
                        if rec_id in n: rl = n
                        elif major_id in n: ml = n
                        elif minor_id and minor_id in n: nil = n
                    if rl and ml and nil:
                        pa = t_a.distance(rl, ml) < t_a.distance(rl, nil)
                        pb = t_b.distance(rl, nil) < t_b.distance(rl, ml)
                        topology = "PASSED" if (pa and pb) else "FAILED"
                        print(f"  Topology: {topology}  (non-recomb ✓={pa} recomb ✓={pb})")
                except Exception as e: print(f"  ⚠ Topo: {e}")
        
        # ── SimPlot ──
        if not args.skip_simplots and len(parents) >= 1:
            try:
                minor = next((p[1] for p in parents if p[0]=='minor'), None)
                major = next((p[1] for p in parents if p[0]=='major'), None)
                make_simplot(rec, minor, major, bp, parents, d / "simplot.png")
                print(f"  ✓ SimPlot")
            except Exception as e: print(f"  ⚠ SimPlot: {e}")
        
        # ── Tree Figure ──
        ta = d / "region_a.treefile"; tb = d / "region_b.treefile"
        if not args.skip_treefigs and os.path.exists(ta) and os.path.exists(tb):
            try:
                make_tree_figure(ta, tb, rec_id, minor_id, major_id, bp, d / "tree_validation.png")
                print(f"  ✓ Tree figure")
            except Exception as e: print(f"  ⚠ TreeFig: {e}")
        
        # ── Cleanup ──
        if args.clean:
            for pattern in ["*.fasta", "*.bionj", "*.ckp.gz", "*.iqtree", "*.log", "*.mldist", "*.model.gz", "*.uniqueseq.phy"]:
                for f in d.glob(pattern):
                    # 2026-09-15: 原为裸 except (会吞 KeyboardInterrupt), 且静默。
                    try:
                        f.unlink()
                    except OSError as _e:
                        print(f"  [cleanup] 跳过 {f.name}: {_e}")
        
        if topology == "PASSED": passed += 1
        elif topology == "FAILED": failed += 1
        
        report.append(f"### Event {n}: {rec_id}\n"
                      f"- Markers: `{evt['markers'] or 'none'}`, BP: {bp[0]}-{bp[1]}\n"
                      f"- Minor: {minor_id or '?'}  Major: {major_id or '?'}\n"
                      f"- Methods: {evt['n_methods']} ({', '.join(evt['methods'].keys())})\n"
                      f"- Topology: **{topology}**\n")
    
    print(f"\nPASSED: {passed}  FAILED: {failed}")
    (outdir / "validation_report.md").write_text('\n'.join(report))

if __name__ == "__main__":
    main()
