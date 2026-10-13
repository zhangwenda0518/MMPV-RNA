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
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request

ST6_TOOLS = os.path.expanduser('~/biosoft/splitstree6-tools/tools')
XVFB = shutil.which('Xvfb') or os.path.expanduser('~/mambaforge/bin/Xvfb')
WEBSERVER = os.path.join(ST6_TOOLS, 'webserver')

# ⚠ webserver 的 HTTP context 名用**下划线**: /draw_distances, /draw_sequences。
#   `/drawDistances` 会被 com.sun.net.httpserver 的**前缀匹配**吃掉 → 落到
#   `/draw` → 报 `Unknown parameter: 'matrix'` (2026-10-12 第 8 轮实测)。
_EP_DRAW_DISTANCES = '/splitstree/draw_distances'
_DIST_BLOCK_RE = re.compile(r'BEGIN\s+DISTANCES\s*;.*?END\s*;\s*\[DISTANCES\]',
                            re.S | re.I)

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


# ──────────────────────────────────────────────────────────────────────────────
# 原生布局坐标 (SplitsTree6 自己的 SplitNetworkLayout)
# ──────────────────────────────────────────────────────────────────────────────
def _ensure_xvfb(display=':99'):
    """Xvfb 没起就起一个 (JavaFX 22+ 无 Monocle, 必须有显示)。返回 Popen 或 None。"""
    if subprocess.run(['bash', '-c', 'pgrep -f "Xvfb %s" >/dev/null' % display],
                      check=False).returncode == 0:
        return None
    if not XVFB:
        raise FileNotFoundError('Xvfb not found')
    p = subprocess.Popen([XVFB, display, '-screen', '0', '1600x1200x24'],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    time.sleep(3)
    return p


def _ensure_webserver(port=8001, display=':99'):
    """SplitsTreeServer 没起就起一个。返回 Popen 或 None。"""
    if subprocess.run(['bash', '-c', 'pgrep -f "SplitsTreeServer" >/dev/null'],
                      check=False).returncode == 0:
        return None
    if not os.path.exists(WEBSERVER):
        raise FileNotFoundError('webserver not found: %s' % WEBSERVER)
    env = dict(os.environ, DISPLAY=display)
    p = subprocess.Popen([WEBSERVER, '-p', str(port)], env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    time.sleep(6)
    return p


def extract_distances(src, out_nex):
    """从 .stree6 里抽出 DISTANCES 块 → 独立 nexus (webserver 只吃 DISTANCES)。"""
    text = open(src, encoding='utf-8', errors='replace').read()
    m = _DIST_BLOCK_RE.search(text)
    if not m:
        return None
    with open(out_nex, 'w', encoding='utf-8') as f:
        f.write('#NEXUS\n' + m.group(0) + '\n')
    return out_nex


def _http_post(port, path, params, timeout):
    """POST 表单。⚠ 必须用 POST: DISTANCES 块 ~190 KB, 塞进 GET query 会被
    sun.net.httpserver 的请求头上限截断 → `ConnectionResetError: [Errno 104]`。"""
    data = urllib.parse.urlencode(params).encode('utf-8')
    req = urllib.request.Request(
        'http://127.0.0.1:%d%s' % (port, path), data=data,
        headers={'Content-Type': 'application/x-www-form-urlencoded'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode('utf-8', 'replace')


def _radial_score(txt):
    """tip 半径的 max/中位 —— 越小越"辐射/星形"。用来在两种布局实现里挑纯 EqualAngle 那个。

    为什么需要: webserver 同一请求会返回**两种实现** (见 native_coords 文档),
    它们在几何上差别极大 —— 实测本数据 449 节点版 max/中位=1.8 (辐射), 514 节点版=4.3 (彗星)。
    选"节点数最多"会挑到后者。**没有标签的节点不参与统计** (只算 tip)。
    """
    lines = txt.rstrip('\n').split('\n')
    try:
        n = int(lines[0].split()[0])
    except (IndexError, ValueError):
        return 0.0
    pts, rs = [], []
    for i in range(1, min(1 + n, len(lines))):
        f = lines[i].split('\t')
        try:
            x, y = float(f[1]), float(f[2])
        except (IndexError, ValueError):
            continue
        pts.append((x, y))
        if len(f) > 3 and f[3].strip():
            rs.append([x, y])
    if len(pts) < 3 or not rs:
        return 0.0
    cx = 0.5 * (min(p[0] for p in pts) + max(p[0] for p in pts))
    cy = 0.5 * (min(p[1] for p in pts) + max(p[1] for p in pts))
    rr = sorted(math.hypot(p[0] - cx, p[1] - cy) for p in rs)
    med = rr[len(rr) // 2]
    return (rr[-1] / med) if med > 0 else 0.0


def native_coords(src, out_coords, width=2886, height=1496, algorithm='nnet',
                  layout='radial', port=8001, attempts=6, timeout=1800,
                  select='radial', verbose=False):
    """取 SplitsTree6 自己的 `SplitNetworkLayout` 坐标 (纯文本)。

    ⚠ 实测 (2026-10-12 第 8/9 轮) —— **同一输入会随机返回两种实现之一**:
        实现 A: 449 节点 / 625 边, bbox 2170x2164 (近正方, 辐射/星芒感)
                = `EqualAngle.apply()` 返回 true (全部 split 判为 circular)
        实现 B: 514 节点 / 829 边, bbox 2406x1486 (宽扁网状, 与 SplitsTree GUI 一致)
                = EqualAngle 判 false → 追加 `ConvexHull.apply()` 重建
      两者**都用满全部 271 个 split** (实测), 差别只在布局算法本身。
      每种实现的坐标集本身**完全固定** (跨进程 md5 一致), 随机性来自 NeighborNet
      的 cycle 判定。故这里连发 `attempts` 次, 再按 `select` 挑:

        select='radial'    (默认, 第 9 轮起) —— 取 **tip 半径 max/中位最小** 者, 即纯
                           EqualAngle 辐射布局。用户 2026-10-12 定夺要这个 (要"辐射感")。
        select='max_nodes' (第 8 轮旧行为) —— 取节点数最多者 = ConvexHull 宽扁版。

    src      : .stree6 或含 DISTANCES 的 nexus
    out_coords: 输出坐标文本
    返回 out_coords; 失败抛异常。
    """
    _ensure_xvfb()
    _ensure_webserver(port=port)
    work = os.path.join(os.path.dirname(os.path.abspath(out_coords)) or '.', '_distances.nex')
    nex = src if re.search(r'BEGIN\s+DISTANCES', open(src, errors='ignore').read(), re.I) \
        else extract_distances(src, work)
    if not nex:
        raise RuntimeError('no DISTANCES block in %s' % src)
    matrix = open(nex, encoding='utf-8', errors='replace').read()
    cands, seen = [], []
    for k in range(max(1, attempts)):
        txt = _http_post(port, _EP_DRAW_DISTANCES, {
            'matrix': matrix, 'algorithm': algorithm, 'layout': layout,
            'output': 'coordinates', 'width': int(width), 'height': int(height),
        }, timeout)
        head = txt.split('\n', 1)[0].split()
        if len(head) != 2:
            raise RuntimeError('unexpected webserver reply: %r' % txt[:200])
        n_n = int(head[0])
        sc = _radial_score(txt)
        seen.append((n_n, round(sc, 2)))
        cands.append((n_n, sc, txt))
        if verbose:
            sys.stderr.write('  [native_coords] try %d: %d nodes, tip r max/med=%.2f\n'
                             % (k + 1, n_n, sc))
    if select == 'max_nodes':
        best = max(cands, key=lambda c: (c[0], -c[1]))
    else:  # 'radial'
        best = min(cands, key=lambda c: (c[1] or 1e9, c[0]))
    with open(out_coords, 'w', encoding='utf-8') as f:
        f.write(best[2])
    sys.stderr.write('  [native_coords] %d attempts (nodes, tip r max/med)=%s -> kept %d nodes '
                     '(%s)\n' % (len(seen), seen, best[0], select))
    return out_coords


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--input', required=True, help='aligned FASTA')
    ap.add_argument('-o', '--outdir', required=True, help='output directory')
    ap.add_argument('-n', '--network', default='neighbornet', choices=list(NETWORKS))
    ap.add_argument('-d', '--distance', default='p', choices=list(DIST_METHODS))
    ap.add_argument('--max-taxa', type=int, default=200, help='taxa cap; subsample beyond this')
    ap.add_argument('--min-taxa', type=int, default=4)
    ap.add_argument('--native-coords', default=None,
                    help='额外输出 SplitsTree6 原生布局坐标 (webserver /draw_distances)')
    ap.add_argument('--native-size', default='2886x1496',
                    help='原生布局画布 WxH (默认 2886x1496, 与 SplitsTree GUI 观感一致)')
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

    if args.native_coords:
        try:
            w, h = (int(v) for v in str(args.native_size).lower().split('x'))
        except Exception:  # noqa: BLE001
            w, h = 2886, 1496
        os.makedirs(os.path.dirname(os.path.abspath(args.native_coords)), exist_ok=True)
        try:
            native_coords(stree6, args.native_coords, width=w, height=h,
                          algorithm='nnet' if args.network == 'neighbornet' else args.network,
                          verbose=True)
        except Exception as e:  # noqa: BLE001 — 坐标失败不影响 nexus 产物
            sys.stderr.write('[warn] native_coords failed: %s\n' % e)
    return 0


if __name__ == '__main__':
    sys.exit(main())
