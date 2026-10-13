#!/usr/bin/env python3
"""从 SplitsTree6 加权坐标 -> 不加权 (拓扑) 布局.

原理 (源码 EqualAngle):
  assignAnglesToEdges : graph.setAngle(e, angles[getSplit(e)])   <- 角度是 **per-split**
  assignCoordinatesToNodesRec:
      p_w = p_v + weight_e * (cos(angle_s), sin(angle_s))
      weight_e = useWeights ? getWeight(e) : 1
  => 拓扑布局 = 同一套角度, 把每条边权重换成 1 重放一次 DFS。

角度获取: **不能**按 split 把各边单位向量求和 —— 同一 split 的平行边在文件里
source/target 朝向不一致, 正负相消 (实测位移塌成 ~0.001)。改为在 DFS 中直接用
**原始坐标的位移方向** atan2(p_w - p_v) 反推该 split 的角度 (同一次遍历顺序 =>
方向必然一致), 并用「原权重重放 == 原坐标」自检。
"""
import math
import sys


def load(p):
    L = open(p).read().rstrip('\n').split('\n')
    n, e = map(int, L[0].split())
    nodes, labels, edges = {}, {}, []
    for i in range(1, 1 + n):
        f = L[i].split('\t')
        nodes[int(f[0])] = (float(f[1]), float(f[2]))
        if len(f) > 3 and f[3].strip():
            labels[int(f[0])] = f[3].strip()
    for i in range(1 + n, 1 + n + e):
        f = L[i].split('\t')
        edges.append((int(f[0]), int(f[1]), int(f[2]), float(f[3]), int(f[4])))
    return nodes, labels, edges


def walk(nodes, edges, root, weight_fn):
    """重放 assignCoordinatesToNodesRec。返回 (pos, theta)。

    ⚠ splitsInPath **每路径回溯**; nodesVisited 全局只挡重复递归;
      pos[w] 在检查 visited **之前**写入 (后到路径会覆盖) —— 全部照抄源码。
    """
    adj = {}
    for eid, u, v, w, s in edges:
        adj.setdefault(u, []).append((v, s, w))
        adj.setdefault(v, []).append((u, s, w))
    pos = {root: (0.0, 0.0)}
    theta, ln_of, conflict = {}, {}, []
    visited, path = set(), set()
    sys.setrecursionlimit(200000)

    def rec(v):
        if v in visited:
            return
        visited.add(v)
        for w_, s, ow in adj.get(v, []):
            if s < 0 or s in path:
                continue
            d = (nodes[w_][0] - nodes[v][0], nodes[w_][1] - nodes[v][1])
            th = math.atan2(d[1], d[0])
            ln = math.hypot(d[0], d[1])
            if s in theta:                       # 同 split 多次遍历 -> 必须同向
                dd = abs((th - theta[s] + math.pi) % (2 * math.pi) - math.pi)
                if dd > 1e-3:
                    conflict.append((s, theta[s], th))
            else:
                theta[s], ln_of[s] = th, ln
            wt = weight_fn(s, ow)
            pv = pos[v]
            pos[w_] = (pv[0] + wt * math.cos(theta[s]), pv[1] + wt * math.sin(theta[s]))
            path.add(s)
            rec(w_)
            path.discard(s)

    rec(root)
    return pos, theta, ln_of, conflict


def main():
    src, dst = sys.argv[1], sys.argv[2]
    nodes, labels, edges = load(src)
    root = min(nodes, key=lambda k: nodes[k][0] ** 2 + nodes[k][1] ** 2)

    pos0, theta, ln_of, conflict = walk(nodes, edges, root, lambda s, ow: ow)

    # 坐标 = scale * raw  (LayoutUtils.normalize(maintainAspectRatio=true) 只做均匀缩放)
    # 用最小二乘拟合 scale:  min Σ |scale*pos0 - orig|²
    num = sum(pos0[k][0] * nodes[k][0] + pos0[k][1] * nodes[k][1]
              for k in nodes if k in pos0)
    den = sum(pos0[k][0] ** 2 + pos0[k][1] ** 2 for k in nodes if k in pos0)
    scale = num / den if den else 1.0
    res = [math.hypot(scale * pos0[k][0] - nodes[k][0], scale * pos0[k][1] - nodes[k][1])
           for k in nodes if k in pos0]
    ratios = sorted(ln_of[s] / ow for _, _, _, ow, s in edges if s in ln_of and ow > 0)
    print('root=%d %s' % (root, labels.get(root, '-')))
    print('自检: 覆盖 %d/%d  拟合 scale=%.1f  最大残差 %.3f' % (len(pos0), len(nodes), scale, max(res)))
    print('  |d|/weight 比值: min=%.1f med=%.1f max=%.1f  (应恒定)'
          % (ratios[0], ratios[len(ratios) // 2], ratios[-1]))
    if max(res) > 0.5:
        worst = max((k for k in nodes if k in pos0),
                    key=lambda k: math.hypot(scale * pos0[k][0] - nodes[k][0],
                                             scale * pos0[k][1] - nodes[k][1]))
        print('!! 残差过大 @ node %d orig=%s fit=%s'
              % (worst, nodes[worst], (scale * pos0[worst][0], scale * pos0[worst][1])))
        print('!! 自检未通过, 不产出')
        return 1

    pos1, _, _, _ = walk(nodes, edges, root, lambda s, ow: 1.0)
    xs = [v[0] for v in pos1.values()]
    ys = [v[1] for v in pos1.values()]
    tr = sorted(math.hypot(pos1[k][0], pos1[k][1]) for k in labels)
    import statistics
    print('不加权: %d 节点 (标签 %d)  bbox %.0f x %.0f  tipR med=%.1f max=%.1f ratio=%.2f'
          % (len(pos1), sum(1 for k in labels if k in pos1),
             max(xs) - min(xs), max(ys) - min(ys),
             statistics.median(tr), tr[-1], tr[-1] / statistics.median(tr)))

    with open(dst, 'w') as f:
        f.write('%d\t%d\n' % (len(pos1), len(edges)))
        for k in sorted(pos1):
            lab = '\t' + labels[k] if k in labels else ''
            f.write('%d\t%.1f\t%.1f%s\n' % (k, pos1[k][0], pos1[k][1], lab))
        for eid, u, v, w, s in edges:
            f.write('%d\t%d\t%d\t%.8f\t%d\n' % (eid, u, v, w, s))
    print('wrote', dst)
    return 0


if __name__ == '__main__':
    sys.exit(main())
