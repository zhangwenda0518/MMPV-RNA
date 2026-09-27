#!/usr/bin/env python3
"""
tree_ops.py — 纯 Newick 树操作（无 ETE/Bio.Phylo 依赖）

递归下降解析 newick → 树节点，提供:
  parse_newick      newick 字符串 → Node (递归下降, 支持枝长/支持值/注释)
  to_newick         Node → newick 字符串
  find_mrca         给定一组 taxa → 它们的共同祖先 (MRCA) 节点
  annotate_mrca     给 MRCA 节点加注释 → 返回带注释 newick
  extract_subtree   提取 MRCA 子树 → newick
  genetic_distance  树上两个 taxa 的枝长距离
  get_leaf_names    所有叶名

参考: YR-MPE streamlined_ete 功能 (但用更可靠的递归下降解析, 无 token+深度遍历 bug)

用法:
  from utils.tree_ops import parse_newick, find_mrca, annotate_mrca, genetic_distance
  tree = parse_newick('(A:0.1,(B:0.2,C:0.3)0.9:0.1);')
  mrca = find_mrca(tree, {'B', 'C'})
  d = genetic_distance(tree, 'A', 'B')
"""

from dataclasses import dataclass, field
from typing import List, Optional, Set


@dataclass
class Node:
    name: str = ""
    length: Optional[float] = None
    children: List["Node"] = field(default_factory=list)
    support: str = ""          # bootstrap 值 / 注释 (如 "100", "&label=xxx")
    parent: Optional["Node"] = None


# ─────────────────────────── 解析 ───────────────────────────

def parse_newick(newick: str) -> Node:
    """递归下降解析 newick → Node 树"""
    s = newick.strip()
    if not s or s == ';':
        raise ValueError('Newick 内容为空')
    if not s.endswith(';'):
        s += ';'
    i = 0
    n = len(s)

    def skip_ws():
        nonlocal i
        while i < n and s[i].isspace():
            i += 1

    def read_label() -> str:
        """读一个标签：单引号包裹（'' 转义）或裸读到分隔符；引号本身不并入名字"""
        nonlocal i
        if i < n and s[i] == "'":
            i += 1
            buf = []
            while i < n:
                if s[i] == "'":
                    if i + 1 < n and s[i + 1] == "'":
                        buf.append("'")
                        i += 2
                        continue
                    i += 1
                    break
                buf.append(s[i])
                i += 1
            # 历史坑: 引号名 'AS_KC430335_2010.6685' 若连引号一起当名字，与
            # Meta 表/另一棵树的裸名字对不上 → 拓扑比较整片错位且看不出来。
            return ''.join(buf)
        j = i
        while j < n and s[j] not in ',):;[':
            j += 1
        label = s[i:j].strip()
        i = j
        return label

    def parse_node() -> Node:
        nonlocal i
        skip_ws()
        if i >= n:
            raise ValueError(f'newick 语法错误 @{i}: 输入提前结束（括号不配对？）')
        node = Node()
        if s[i] == '(':
            i += 1  # '('
            # 解析子节点
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                skip_ws()
                if i < n and s[i] == ',':
                    i += 1
                    continue
                if i < n and s[i] == ')':
                    i += 1
                    break
                raise ValueError(f'newick 语法错误 @{i}: 期望 , 或 )')
            # `)` 与枝长之间的部分：支持值 / 注释 / 内部节点名，可能叠加
            # （如 `Ningxia[&posterior=0.9]`），循环吃到分隔符为止。
            skip_ws()
            while i < n:
                if s[i] == '[':
                    j = s.find(']', i)
                    if j == -1:
                        node.support = (node.support or '') + s[i:]
                        i = n
                    else:
                        node.support = (node.support or '') + s[i:j + 1]
                        i = j + 1
                elif s[i] in ':,);':
                    break
                else:
                    j = i
                    # 历史坑: 只读整数 → MrBayes/TreeAnnotator 小数后验概率 (如 0.9) 停在 '.',
                    # 不命中终止符 → ValueError。数字解析含小数点/负号/指数。
                    while j < n and (s[j].isdigit() or s[j] in '.-+eE'):
                        j += 1
                    if j > i and (j >= n or s[j] in ':,);'):
                        try:
                            float(s[i:j])
                            node.support = s[i:j]
                            i = j
                            skip_ws()
                            continue
                        except ValueError:
                            pass
                    # 历史坑: 只认数字/[&…] → 带名内部节点 `(A,B)Ningxia` 的名字
                    # 留在原地，父层会报「期望 , 或 )」；名字读到分隔符为止。
                    before = i
                    node.name = read_label()
                    if i == before:
                        break
                skip_ws()
        else:
            # 叶名（引号名剥引号；裸名读到 , ) : ; 或 BEAST 注释起点 [）
            node.name = read_label()
            # 历史坑: tip 注释 A[&location=China] 被并入叶名 → 跨树比对名字错位; 归入 support
            if i < n and s[i] == '[':
                j2 = s.find(']', i)
                if j2 == -1:
                    node.support = s[i:]
                    i = n
                else:
                    node.support = s[i:j2 + 1]
                    i = j2 + 1
        # 枝长
        skip_ws()
        if i < n and s[i] == ':':
            i += 1
            skip_ws()
            j = i
            while j < n and (s[j].isdigit() or s[j] in '.-+eE'):
                j += 1
            if j > i:
                try:
                    node.length = float(s[i:j])
                except ValueError:
                    pass
                i = j
        skip_ws()
        return node

    # BEAST/Nexus 树前缀注释（如 `[&R] (((…`）：不剥掉会被当叶名，整棵树解析错位
    skip_ws()
    while i < n and s[i] == '[':
        j = s.find(']', i)
        if j == -1:
            raise ValueError(f'newick 前缀注释未闭合 @{i}: {s[i:i + 30]!r}')
        i = j + 1
        skip_ws()

    root = parse_node()
    # 解析不完就报错，不许「截断成子树」返回：括号不配对 / 顶层多个逗号分隔子树
    # 时，旧实现把剩余分支静默丢掉 —— 下游（拓扑比较 / MRCA / 遗传距离）会拿着
    # 一棵少了样品的树照常出结果，且看不出来。同口径见平台 phylogeo.parse_newick。
    skip_ws()
    if i < n and s[i] == ';':        # 收尾分号（本函数会自动补，故此处必然在）
        i += 1
        skip_ws()
    if i < n:
        raise ValueError(
            f"Newick 解析未走完：位置 {i} 处还剩 {n - i} 个字符没解析 "
            f"({s[i:i + 30]!r}) —— 多半是括号不配对或顶层有多个逗号分隔的子树；"
            '继续下去会把剩余分支静默丢掉，故直接报错')
    return root


def _format_name(name: str) -> str:
    """名字含分隔符/空格/引号时加单引号（内部 ' 翻倍转义），保证输出能再解析"""
    if name and any(ch in name for ch in "(),:;[]' \t"):
        return "'" + name.replace("'", "''") + "'"
    return name


def to_newick(node: Node) -> str:
    """Node → newick 字符串"""
    if node.children:
        inner = ','.join(to_newick(c) for c in node.children)
        s = f'({inner})'
        if node.name:
            s += _format_name(node.name)
        if node.support:
            s += node.support
    else:
        s = _format_name(node.name)
    if node.length is not None:
        # repr 而非 %g：%g 只保留 6 位有效数字，往返会把枝长静默改掉
        # （如 0.000478047 → 0.000478），下游建树/定年会用到这些数值。
        s += f':{node.length!r}'
    return s


# ─────────────────────────── 树操作 ───────────────────────────

def get_leaf_names(node: Node) -> Set[str]:
    if not node.children:
        return {node.name} if node.name else set()
    out = set()
    for c in node.children:
        out |= get_leaf_names(c)
    return out


def find_mrca(node: Node, taxa: Set[str]) -> Optional[Node]:
    """找包含全部 taxa 的最小共同祖先 (MRCA)"""
    taxa = set(taxa)

    def dfs(n: Node):
        if not n.children:
            return {n.name} if n.name in taxa else set()
        found = set()
        for c in n.children:
            found |= dfs(c)
        if taxa.issubset(found) and len(found) > 0:
            # 这是第一个(最浅)包含全部 taxa 的节点
            return found
        return found

    # 递归找 MRCA: 后序遍历, 第一个其子树包含全部 taxa 的节点
    result = None

    def postorder(n: Node):
        nonlocal result
        if result is not None:
            return set()
        if not n.children:
            return {n.name} if n.name in taxa else set()
        found = set()
        for c in n.children:
            found |= postorder(c)
        if result is None and taxa and taxa.issubset(found):
            result = n
        return found

    postorder(node)
    return result


def annotate_mrca(newick: str, taxa: Set[str], label: str = "[&MRCA]") -> str:
    """给 MRCA 节点加注释, 返回带注释 newick"""
    tree = parse_newick(newick)
    mrca = find_mrca(tree, taxa)
    if mrca is None:
        return newick
    if not mrca.support:
        mrca.support = label
    elif label not in mrca.support:
        mrca.support += label
    return to_newick(tree) + ';'


def extract_subtree(newick: str, taxa: Set[str]) -> str:
    """提取 MRCA 子树 → newick"""
    tree = parse_newick(newick)
    mrca = find_mrca(tree, taxa)
    if mrca is None:
        return ''
    return to_newick(mrca) + ';'


def _path_to_root(node: Node):
    path = []
    cur = node
    while cur is not None:
        path.append(cur)
        cur = cur.parent
    return path


def _find_leaf(node: Node, name: str) -> Optional[Node]:
    if not node.children:
        return node if node.name == name else None
    for c in node.children:
        r = _find_leaf(c, name)
        if r:
            return r
    return None


def genetic_distance(newick: str, taxa1: str, taxa2: str) -> Optional[float]:
    """树上两个 taxa 的枝长距离 (patristic distance)"""
    tree = parse_newick(newick)
    n1 = _find_leaf(tree, taxa1)
    n2 = _find_leaf(tree, taxa2)
    if n1 is None or n2 is None:
        return None
    p1 = _path_to_root(n1)
    p2 = _path_to_root(n2)
    # LCA: 两条路径第一个共同节点
    # 2026-09-15 (审查 P2-19) 修复: 旧版写 `if a in p2` —— list.__contains__ 走
    # 元素 __eq__; Node 是 dataclass, 默认按字段比较, 两个**结构相同**但位置不同的
    # 节点会判等 → LCA 可能取到错误的节点 (且逐个比较代价 O(n²))。
    # 改用 id() 集合按对象身份判定。
    p2_ids = {id(x) for x in p2}
    lca = None
    for a in p1:
        if id(a) in p2_ids:
            lca = a
            break
    if lca is None:
        return None
    dist = 0.0

    def sum_len(node, stop):
        d = 0.0
        cur = node
        while cur is not None and cur is not stop:
            d += cur.length or 0.0
            cur = cur.parent
        return d

    return sum_len(n1, lca) + sum_len(n2, lca)


def sort_topology(newick: str, ascending: bool = True) -> str:
    """按子树叶数排序拓扑 (ascending=True 叶少的 clade 在前)"""
    tree = parse_newick(newick)

    def rec_sort(node: Node):
        for c in node.children:
            rec_sort(c)
        node.children.sort(key=lambda c: len(get_leaf_names(c)),
                           reverse=not ascending)

    rec_sort(tree)
    return to_newick(tree) + ';'


def _norm_cdf_not_used():
    pass


if __name__ == '__main__':
    import sys
    nwk = sys.argv[1] if len(sys.argv) > 1 else '(A:0.1,(B:0.2,C:0.3)0.9:0.1);'
    t = parse_newick(nwk)
    print('叶名:', sorted(get_leaf_names(t)))
    mrca = find_mrca(t, {'B', 'C'})
    print('MRCA(B,C):', to_newick(mrca) if mrca else None)
    print('距离(A,B):', genetic_distance(nwk, 'A', 'B'))
    print('排序:', sort_topology(nwk))
