#!/usr/bin/env python3
"""
splitstree_bridge.py -- SplitsTree6 分裂网络桥接工具
将比对 FASTA 转为完整 .stree6 workflow (数据+算法链嵌入), 通过 workflow-run CLI
执行 Neighbor Net 分析, 输出 Nexus 网络。

依赖:
  - splitstree6-tools 6.9.5+ (~/biosoft/splitstree6-tools/tools/workflow-run)
  - Xvfb (JavaFX 22+ 移除 Monocle, 需虚拟显示)
  - python3 + Bio (仅用于读 fasta; 纯文本解析亦可)

用法:
  python3 splitstree_bridge.py -i alignment.fasta -o out_dir [--network neighbornet|median_joining|parsimony] [--max-taxa 200]

管线集成: 在 align (MAFFT) 之后运行; 产物 out_dir/neighbornet.nexus
"""

import argparse
import os
import shutil
import subprocess
import sys
import time

ST6_TOOLS = os.path.expanduser('~/biosoft/splitstree6-tools/tools')
XVFB = shutil.which('Xvfb') or os.path.expanduser('~/mambaforge/bin/Xvfb')

# SplitsTree6 算法真名 (jar 类名), 来自 st6_alglist.sh 实测
DIST_METHODS = {
    'p': 'PDistance',
    'hamming': 'HammingDistance',
    'jc': 'JukesCantorDistance',
    'k2p': 'K2PDistance',
}
NETWORKS = {
    'neighbornet': ('Neighbor Net', 'DISTANCES'),
    'median_joining': ('Median Joining', 'CHARACTERS'),
    'parsimony': ('Parsimony Splits', 'CHARACTERS'),
}


def read_fasta(path):
    """极简 FASTA 解析, 避免依赖 Biopython."""
    labels, seqs = [], []
    cur = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                cur = line[1:].split()[0].replace(' ', '_')[:32]
                labels.append(cur)
                seqs.append([])
            elif line and cur is not None:
                seqs[-1].append(line.upper())
    return labels, [''.join(s) for s in seqs]


def pdistance_matrix(seqs):
    n = len(seqs)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            diff = same = 0
            for a, b in zip(seqs[i], seqs[j]):
                if a in '-?N' or b in '-?N':
                    continue
                if a != b:
                    diff += 1
                same += 1
            D[i][j] = D[j][i] = diff / same if same else 0.0
    return D


def sanitize(label):
    return ''.join(c if c.isalnum() or c == '_' else '_' for c in label)[:32]


def build_stree6(labels, seqs, distances, net_alg, dist_name):
    """生成完整 .stree6: TAXA -> Characters -> PDistance -> Distances -> Network."""
    n, L = len(labels), max(len(s) for s in seqs)
    labels = [sanitize(x) for x in labels]
    # 去重 (NEXUS label 冲突)
    seen = {}
    for k, lb in enumerate(labels):
        if lb in seen:
            seen[lb] += 1
            labels[k] = '%s_%d' % (lb, seen[lb])
        else:
            seen[lb] = 0
    D = pdistance_matrix(seqs) if distances is None else distances
    out = []
    w = out.append
    w('#NEXUS [SplitsTree6]\n')
    w('\nBEGIN SPLITSTREE6;\n\tDIMENSIONS nDataNodes=6 nAlgorithms=4;\n')
    w("\tPROGRAM version='SplitsTree App (version 6.9.5)';\n\tWORKFLOW creationDate='%d';\nEND; [SPLITSTREE6]\n" % int(time.time() * 1000))
    # Input Taxa
    w('\nBEGIN TAXA;\n\tTITLE \'Input Taxa\';\n\tLINK ALGORITHM = \'Input Data Loader\';\n\tDIMENSIONS ntax=%d;\n\tTAXLABELS\n' % n)
    for i, lb in enumerate(labels, 1):
        w("\t[%d] '%s'\n" % (i, lb))
    w(';\nEND; [TAXA]\n')
    w('\nBEGIN ALGORITHM;\n\tTITLE \'Taxa Filter\';\n\tLINK ALGORITHM = \'Input Data Loader\';\n\tNAME \'Taxa Filter\';\n\tOPTIONS\n\t;\nEND; [ALGORITHM]\n')
    # Working Taxa
    w('\nBEGIN TAXA;\n\tTITLE \'Working Taxa\';\n\tLINK ALGORITHM = \'Taxa Filter\';\n\tDIMENSIONS ntax=%d;\n\tTAXLABELS\n' % n)
    for i, lb in enumerate(labels, 1):
        w("\t[%d] '%s'\n" % (i, lb))
    w(';\nEND; [TAXA]\n')
    # Input Characters
    w('\nBEGIN CHARACTERS;\n\tTITLE \'Input Characters\';\n\tLINK ALGORITHM = \'Input Data Loader\';\n\tDIMENSIONS ntax=%d nchar=%d;\n\tFORMAT\n\t\tdatatype=DNA missing=? gap=- labels=left;\n\tMATRIX\n' % (n, L))
    for lb, s in zip(labels, seqs):
        w('%s %s\n' % (lb, s))
    w(';\nEND; [CHARACTERS]\n')
    w('\nBEGIN ALGORITHM;\n\tTITLE \'Input Data Filter\';\n\tLINK CHARACTERS = \'Input Characters\';\n\tNAME \'Characters Taxa Filter\';\n\tOPTIONS\n\t;\nEND; [ALGORITHM]\n')
    # Working Characters
    w('\nBEGIN CHARACTERS;\n\tTITLE \'Working Characters\';\n\tLINK ALGORITHM = \'Input Data Filter\';\n\tDIMENSIONS ntax=%d nchar=%d;\n\tFORMAT\n\t\tdatatype=DNA missing=? gap=- labels=left;\n\tMATRIX\n' % (n, L))
    for lb, s in zip(labels, seqs):
        w('%s %s\n' % (lb, s))
    w(';\nEND; [CHARACTERS]\n')

    if net_alg == 'Neighbor Net':
        w('\nBEGIN ALGORITHM;\n\tTITLE \'%s\';\n\tLINK CHARACTERS = \'Working Characters\';\n\tNAME \'%s\';\n\tOPTIONS\n\t;\nEND; [ALGORITHM]\n' % (dist_name, dist_name))
        w('\nBEGIN DISTANCES;\n\tTITLE \'Distances\';\n\tLINK ALGORITHM = \'%s\';\n\tDIMENSIONS ntax=%d;\n\tFORMAT labels=left diagonal triangle=Both;\n\tMATRIX\n' % (dist_name, n))
        for i, lb in enumerate(labels, 1):
            w('[%d] \'%s\' %s\n' % (i, lb, ' '.join('%.6f' % x for x in D[i - 1])))
        w(';\nEND; [DISTANCES]\n')
        w('\nBEGIN ALGORITHM;\n\tTITLE \'Neighbor Net\';\n\tLINK DISTANCES = \'Distances\';\n\tNAME \'Neighbor Net\';\n\tOPTIONS\n\t;\nEND; [ALGORITHM]\n')
        w('\nBEGIN SPLITS;\n\tTITLE \'Splits\';\n\tLINK ALGORITHM = \'Neighbor Net\';\n\tDIMENSIONS ntax=%d nsplits=0;\n\tCYCLE\n\t%s\n;\nEND; [SPLITS]\n' % (n, ' '.join(str(i + 1) for i in range(n))))
    else:
        # characters2network: Median Joining / Parsimony Splits
        w('\nBEGIN ALGORITHM;\n\tTITLE \'%s\';\n\tLINK CHARACTERS = \'Working Characters\';\n\tNAME \'%s\';\n\tOPTIONS\n\t;\nEND; [ALGORITHM]\n' % (net_alg, net_alg))
        w('\nBEGIN NETWORK;\n\tTITLE \'Network\';\n\tLINK ALGORITHM = \'%s\';\n\tDIMENSIONS nVertices=0 nEdges=0;\n\tTYPE HaplotypeNetwork;\nEND; [NETWORK]\n' % net_alg)
    return ''.join(out)


def run_workflow(stree6_path, out_nexus, timeout=1800):
    """启动 Xvfb, 运行 workflow-run, 输出 nexus."""
    if not os.path.exists(os.path.join(ST6_TOOLS, 'workflow-run')):
        raise FileNotFoundError('workflow-run not found in %s' % ST6_TOOLS)
    xpid = None
    env = dict(os.environ, DISPLAY=':99')
    if subprocess.run(['bash', '-c', 'pgrep -f "Xvfb :99" >/dev/null'], check=False).returncode != 0:
        if not XVFB:
            raise FileNotFoundError('Xvfb not found')
        xpid = subprocess.Popen([XVFB, ':99', '-screen', '0', '1024x768x24'],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        time.sleep(2)
    try:
        cmd = [os.path.join(ST6_TOOLS, 'workflow-run'),
               '-w', stree6_path, '-i', stree6_path, '-o', out_nexus,
               '-t', '%ds' % timeout]
        r = subprocess.run(cmd, env=env, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=timeout + 60)
        if not os.path.exists(out_nexus) or os.path.getsize(out_nexus) < 100:
            sys.stderr.write(r.stdout[-2000:] + '\n' + r.stderr[-2000:])
            raise RuntimeError('workflow-run failed (rc=%d)' % r.returncode)
        return r.returncode
    finally:
        if xpid:
            xpid.terminate()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--input', required=True, help='aligned FASTA')
    ap.add_argument('-o', '--outdir', required=True, help='output directory')
    ap.add_argument('-n', '--network', default='neighbornet', choices=list(NETWORKS))
    ap.add_argument('-d', '--distance', default='p', choices=list(DIST_METHODS))
    ap.add_argument('--max-taxa', type=int, default=200, help='taxa cap; subsample beyond this')
    ap.add_argument('--min-taxa', type=int, default=4)
    args = ap.parse_args()

    labels, seqs = read_fasta(args.input)
    if len(labels) > args.max_taxa:
        step = len(labels) / args.max_taxa
        idx = [int(i * step) for i in range(args.max_taxa)]
        labels = [labels[i] for i in idx]
        seqs = [seqs[i] for i in idx]
        sys.stderr.write('subsampled to %d taxa\n' % len(labels))
    if len(labels) < args.min_taxa:
        sys.stderr.write('skip: only %d taxa (<%d)\n' % (len(labels), args.min_taxa))
        return 1

    os.makedirs(args.outdir, exist_ok=True)
    net_alg, _ = NETWORKS[args.network]
    dist_name = DIST_METHODS[args.distance]
    D = pdistance_matrix(seqs) if args.network == 'neighbornet' else None
    stree6 = os.path.join(args.outdir, '%s.stree6' % args.network)
    with open(stree6, 'w', encoding='utf-8') as f:
        f.write(build_stree6(labels, seqs, D, net_alg, dist_name))

    out_nexus = os.path.join(args.outdir, '%s.nexus' % args.network)
    t0 = time.time()
    run_workflow(stree6, out_nexus)
    sys.stderr.write('%s done: %s (%.1fs)\n' % (args.network, out_nexus, time.time() - t0))
    return 0


if __name__ == '__main__':
    sys.exit(main())
