#!/usr/bin/env python3
"""
utils/import_export.py — 引入 / 本地传播 / 输出（importation · local transmission · export）
==========================================================================================
从一棵定年树（BEAST MCC / IQ-TREE / 任意 Newick）+ 一条「序列 → 区划」元数据表，
按 **1/group_size 规则** 拆出每个区划的：

    importations          引入事件数（一个测序簇 = 一次引入，按簇大小摊权重）
    local_transmissions   本地传播事件数（簇内除首例外的其余病例）
    exportations          输出事件数（本区划的簇从别的区划进入后，再迁出的次数）
    以及 corridors 表（from → to 走廊，含原始边计数与按 1/group_size 的加权计数）

算法出处（逐条可追溯）
----------------------
1. **1/group_size 规则** ——
   `git-repo/SARS-CoV-2_phylogeo/docs/Analyses/R_functions/
    Relative_load_importation_local_transmission.R:15-30` 原文（行号已核实）：

        group_size = Result$Lineage_sizes[group_no]
        importations        = c(importations, 1/group_size)
        local_transmissions = c(local_transmissions, 1-1/group_size)

   即：一个谱系簇（group）视为「1 次引入 + (g-1) 次本地传播」，
   簇内每条序列分摊 1/g 次引入、1-1/g 次本地传播。这样避免了
   「每个簇只算 1 次引入」在簇大小差异大时的偏倚（SARS 原文注释：
   "If the group size is 2, then it is +0.5 importations, and +0.5 local
   transmissions"）。
   本模块的簇定义（可复现的推广，原文用外部谱系分组，此处从树上推）：
   **区划簇 = 最大同区划单系枝**（子树内所有 tip 同属一个区划，而其父节点
   子树已跨区划）。这一推广保证「簇大小」在树上唯一确定、可手算复核。

2. **区划状态标注**（可选）—— 与管线 `geo_analysis.parse_mcc_states()` 读取的
   BEAST 2 离散性状注释格式兼容：`.set={A,B}` / `.set.prob={0.6,0.4}`；
   也支持 MCC 内部节点标签直接写区划名（`Hainan` / `{Hainan}` / `Hainan/0.98`）。
   带 `--use-node-labels` 时内部节点优先用标签状态，否则一律用子代 tip 的
   多数表决（majority rule，平票时取「子枝 tip 数最多者」，再平票按名字序）。

3. **走廊计数** —— 两种口径并列输出：
   * `n_edges`       ：树上「父节点多数区划 ≠ 子节点多数区划」的边数
                       （等价于一次离散性状状态转移，spread3/BSSVS 走廊口径）；
   * `tip_weighted`  ：把每个区划簇的 1 次迁出按 1/group_size 摊到 tip 上
                       （与 importations 同规则，两者可加和自洽）。

本模块**自包含**（只依赖标准库 + 可选 biopython 自检），不 import 管线内其它
模块，纯新增，不改动任何既有行为。

用法
----
  # 1) 自测（玩具树 + 手算真值断言）
  python utils/import_export.py --selftest

  # 2) 真实数据
  python utils/import_export.py --tree mcc.tree --regions regions.tsv \
      --outdir import_export_out [--use-node-labels]

  # regions.tsv 两列（列名自动识别）：
  #   tip<TAB>region[<TAB>date]      # date 可选，给了才出时序表

输出
----
  <outdir>/import_export.tsv           # 每区划：引入/本地/输出/净流量
  <outdir>/corridors.tsv               # from,to,n_edges,tip_weighted
  <outdir>/region_clades.tsv           # 每个区划簇：大小 / 每 tip 权重 / tip 列表
  <outdir>/import_export_timeseries.tsv  # 仅当有日期列
  <outdir>/import_export_summary.json  # 总计 + 算法出处 + 输入指纹
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, OrderedDict
from typing import Dict, List, Optional, Sequence, Tuple

# 脚本方式运行 (python utils/import_export.py --selftest) 时 sys.path[0] 是 utils/，
# 下面的 `from utils.decimal_year import ...` 会 ModuleNotFoundError
# （同 auto_parse_params.py 记的那个坑）→ 把管线根目录补进 sys.path。
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# ═══════════════════════════════════════════════════════════════════
# 1. 极简 Newick 解析（自包含；支持引号标签 / [注释] / 分支长度 / [&R]）
# ═══════════════════════════════════════════════════════════════════
class Node:
    """树节点：name 标签、length 分支长、children 子节点、id 内部编号。"""

    __slots__ = ("name", "length", "children", "parent", "id", "n_tips", "regions", "majority")

    def __init__(self, name: str = "", length: Optional[float] = None):
        self.name = name
        self.length = length
        self.children: List["Node"] = []
        self.parent: Optional["Node"] = None
        self.id = -1
        self.n_tips = 0
        self.regions: Optional[Counter] = None
        self.majority: Optional[str] = None

    @property
    def is_leaf(self) -> bool:
        return not self.children


def _strip_comments(text: str) -> str:
    """去掉 Newick 注释 `[...]`（含 BEAST 的 [&rate=...]、[&R] 根标记）。"""
    out = []
    depth = 0
    for ch in text:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    return "".join(out)


def parse_newick(newick: str) -> Node:
    """解析 Newick → 根节点。支持：引号标签、[注释]、分支长度、多分叉。"""
    s = _strip_comments(newick).strip()
    if not s:
        raise ValueError("空 Newick")
    if s.endswith(";"):
        s = s[:-1]
    pos = 0

    def read_label() -> str:
        nonlocal pos
        if pos < len(s) and s[pos] == "'":
            pos += 1
            buf = []
            while pos < len(s):
                if s[pos] == "'":
                    if pos + 1 < len(s) and s[pos + 1] == "'":
                        buf.append("'")
                        pos += 2
                        continue
                    pos += 1
                    break
                buf.append(s[pos])
                pos += 1
            return "".join(buf)
        start = pos
        while pos < len(s) and s[pos] not in ",():;":
            pos += 1
        return s[start:pos].strip()

    def read_length() -> Optional[float]:
        nonlocal pos
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",()":
                pos += 1
            raw = s[start:pos].strip()
            try:
                return float(raw)
            except ValueError:
                return None
        return None

    def parse_subtree() -> Node:
        nonlocal pos
        if pos < len(s) and s[pos] == "(":
            pos += 1
            node = Node()
            while True:
                child = parse_subtree()
                child.parent = node
                node.children.append(child)
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                if pos < len(s) and s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"Newick 语法错误：位置 {pos} 期待 ',' 或 ')'")
            node.name = read_label()
            node.length = read_length()
            return node
        leaf = Node(name=read_label())
        leaf.length = read_length()
        if not leaf.name:
            raise ValueError(f"Newick 语法错误：位置 {pos} 出现空 tip 名")
        return leaf

    root = parse_subtree()
    if pos < len(s):
        raise ValueError(f"Newick 尾部有多余内容：{s[pos:pos+30]!r}")
    _index(root)
    return root


def _index(root: Node) -> None:
    """后序遍历：编号 + 每节点 tip 数。"""
    stack = [(root, False)]
    order = []
    while stack:
        n, done = stack.pop()
        if done:
            order.append(n)
        else:
            stack.append((n, True))
            for c in n.children:
                stack.append((c, False))
    for i, n in enumerate(order):
        n.id = i
        n.n_tips = 1 if n.is_leaf else sum(c.n_tips for c in n.children)


def iter_tips(root: Node) -> List[Node]:
    """按 Newick 出现顺序返回所有叶子。"""
    out: List[Node] = []
    stack = [root]
    while stack:
        n = stack.pop()
        if n.is_leaf:
            out.append(n)
        else:
            stack.extend(reversed(n.children))
    return out


# ═══════════════════════════════════════════════════════════════════
# 2. 区划映射 / 多数表决
# ═══════════════════════════════════════════════════════════════════
def _read_table(path: str) -> List[Dict[str, str]]:
    """读 CSV/TSV（嗅探分隔符）→ 行字典列表。"""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            delim = dialect.delimiter
        except csv.Error:
            delim = "\t" if "\t" in sample else ","
        return list(csv.DictReader(fh, delimiter=delim))


def _pick(keys: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    norm = {str(k).strip().lower().replace(" ", "_"): k for k in keys}
    for c in candidates:
        if c in norm:
            return norm[c]
    for c in candidates:
        for k_norm, k in norm.items():
            if k_norm.startswith(c):
                return k
    return None


def load_region_table(path: str) -> Tuple[Dict[str, str], Dict[str, float]]:
    """读元数据表 → ({tip: region}, {tip: decimal_year})。列名自动识别。"""
    rows = _read_table(path)
    if not rows:
        return {}, {}
    cols = list(rows[0].keys())
    tip_c = _pick(cols, ("tip", "name", "id", "strain", "accession", "taxon", "seq"))
    reg_c = _pick(cols, ("region", "location", "locality", "state", "site", "province", "place"))
    date_c = _pick(cols, ("date", "collection_date", "collectiondate", "sampling_date", "time", "year"))
    if tip_c is None or reg_c is None:
        raise ValueError(f"{path}: 需要 tip 与 region 两类列，实际列 {cols}")
    regions: Dict[str, str] = {}
    dates: Dict[str, float] = {}
    for r in rows:
        tip = str(r.get(tip_c, "")).strip()
        reg = str(r.get(reg_c, "")).strip()
        if not tip or not reg:
            continue
        regions[tip] = reg
        if date_c:
            yr = _to_decimal_year(str(r.get(date_c, "")).strip())
            if yr is not None:
                dates[tip] = yr
    return regions, dates


def _to_decimal_year(raw: str) -> Optional[float]:
    """'2016' / '2016-07' / '2016-07-15' / '15/07/2016' / '2016/07/15' → 十进制年。

    2026-09-16 (回归 D11 修复): 原实现自带一套年/月日算术, 与全管线唯一实现
    `utils/decimal_year.to_decimal_year` 不一致, 触发 `verify_p2b_decimal.py` 的
    D11 FAIL。已改为**委托统一实现**, 仅保留本模块特有的两项前置处理:

      ① **接受 DD/MM/YYYY 与 DD.MM.YYYY** —— 原实现支持, 统一实现不支持
         (`15/07/2016` 会被它当成 day=2016 越界而抛异常)。本模块的 R 脚本
         算法出处常用该格式, 故保留, 并**先归一化成 ISO 再委托**。
      ② **dry 时不静默编造** —— 原实现在无法解析时用
         `re.search(r"(19|20)\\d{2}") + 0.5` **凭正则猜一个"年中"**,
         这会把 `xxxx1999yy` 这种乱串也算成 1999.5。已删除, 改为返回 None
         (与统一实现的 strict=False 一致: 脏数据不造年份)。

    ⚠ 本模块的年份语义说明: 原实现对 `YYYY` 返回**年初**(2016.0),
    对 `YYYY-MM` 返回**月初**。改为委托后按统一口径返回**年中/月中**
    (2016.4556 / 2016.5376)。**这不影响本模块输出** —— 唯一消费者
    `timeseries()` 在 `:550` 只取 `int(yr)`, 月/日部分从未参与分箱。
    已实测: 三档格式经 `int()` 后新旧一致。
    """
    if not raw or raw.lower() in ("na", "nan", "none", "not applicable", "unknown"):
        return None
    s = str(raw).strip().replace("/", "-").replace(".", "-")
    # ① DD-MM-YYYY → YYYY-MM-DD (统一实现只认 YYYY 在前)
    m = re.match(r"^(\d{1,2})-(\d{1,2})-(\d{4})$", s)
    if m:
        d, mo, y = m.group(1), m.group(2), m.group(3)
        s = f"{y}-{mo}-{d}"
    # ② 委托统一实现; strict=False → 解析失败返回 None (不静默编造)
    from utils.decimal_year import to_decimal_year
    return to_decimal_year(s, sample="", strict=False)


_STATE_LABEL_RE = re.compile(r"^[\{\[]?\s*([^\[\{\}/,;:]+?)\s*[\}\]]?\s*(?:[/|].*)?$")


def node_state_label(name: str, known_regions: Sequence[str]) -> Optional[str]:
    """内部节点标签 → 区划名（兼容 `Hainan` / `{Hainan}` / `Hainan/0.98`）。

    只有匹配已知区划集合才返回，避免把 posterior 值当状态。
    """
    if not name:
        return None
    m = _STATE_LABEL_RE.match(name.strip())
    cand = m.group(1).strip() if m else name.strip()
    if cand in known_regions:
        return cand
    return None


def annotate_tree(root: Node, tip_regions: Dict[str, str],
                  use_node_labels: bool = False) -> Tuple[int, int]:
    """后序标注每个节点的多数区划（`node.regions` 为 Counter）。

    Returns (n_unassigned_tips, n_tips)。
    平票规则：取「子枝 tip 数最多」的子节点区划；再平票取名字序最小者。
    """
    tips = iter_tips(root)
    n_total, n_unassigned = len(tips), 0
    known = sorted(set(tip_regions.values()))
    for leaf in tips:
        reg = tip_regions.get(leaf.name)
        if reg is None:
            n_unassigned += 1
            leaf.regions = Counter()
        else:
            leaf.regions = Counter({reg: 1})

    order: List[Node] = []
    stack = [(root, False)]
    while stack:
        n, done = stack.pop()
        if done:
            order.append(n)
        else:
            stack.append((n, True))
            stack.extend((c, False) for c in n.children)

    for node in order:
        if node.is_leaf:
            continue
        counts: Counter = Counter()
        for c in node.children:
            if c.regions:
                counts.update(c.regions)
        node.regions = counts
        if use_node_labels:
            lab = node_state_label(node.name, known)
            if lab:
                node.majority = lab  # 标签优先，仅用于走廊口径
                continue
        node.majority = _majority(counts, node)
    # 叶子节点也补 majority
    for leaf in tips:
        leaf.majority = next(iter(leaf.regions), None)
    return n_unassigned, n_total


def _majority(counts: Counter, node: Node) -> Optional[str]:
    """多数表决 + 确定性平票规则（子枝 tip 数最多者 → 名字序）。"""
    if not counts:
        return None
    best = max(counts.values())
    cands = sorted(r for r, c in counts.items() if c == best)
    if len(cands) == 1:
        return cands[0]
    kids = sorted(node.children, key=lambda c: (-c.n_tips, c.name))
    for c in kids:
        maj = getattr(c, "majority", None)
        if maj in cands:
            return maj
    return cands[0]


# ═══════════════════════════════════════════════════════════════════
# 3. 核心算法：区划簇 → importations / local_transmissions
# ═══════════════════════════════════════════════════════════════════
def find_regional_clades(root: Node) -> List[Dict]:
    """找最大同区划单系枝（区划簇）。

    判定：节点子树内所有**已指派** tip 同属一个区划，且父节点子树已跨区划
    （或该节点是根）。范围外的未指派 tip 不参与判定（计入 unassigned）。

    Returns
    -------
    list[dict]: clade_root, region, size, weight_per_tip (=1/size), tips
    """
    clades: List[Dict] = []

    def node_pure_region(n: Node) -> Optional[str]:
        if not n.regions:
            return None
        if len(n.regions) == 1 and sum(n.regions.values()) == n.n_tips:
            return next(iter(n.regions))
        return None

    stack = [(root, False)]
    order: List[Node] = []
    while stack:
        n, done = stack.pop()
        if done:
            order.append(n)
        else:
            stack.append((n, True))
            stack.extend((c, False) for c in n.children)

    region_of_clade = {}
    for n in order:
        pr = node_pure_region(n)
        if pr is None:
            continue
        parent_pure = node_pure_region(n.parent) if n.parent is not None else None
        if parent_pure == pr:
            continue  # 非最大（父节点仍是纯的）→ 交给父节点
        # 纯节点 ⇒ 子树内每个 tip 均已指派区划，直接收集
        tips = [t.name for t in iter_tips(n)]
        if not tips:
            continue
        region_of_clade[n.id] = pr
        clades.append({
            "clade_root": n.name or f"node{n.id}",
            "node_id": n.id,
            "region": pr,
            "size": len(tips),
            "weight_per_tip": 1.0 / len(tips),
            "tips": tips,
        })
    return clades


def summarize_import_export(root: Node, tip_regions: Dict[str, str],
                            use_node_labels: bool = False) -> Dict:
    """主入口：树 + 区划表 → 每区划 importations / local_transmissions /
    exportations + corridors + 区划簇明细。

    1/group_size 规则出处：SARS-CoV-2_phylogeo
    `docs/Analyses/R_functions/Relative_load_importation_local_transmission.R:15-30`。
    """
    n_unassigned, n_tips = annotate_tree(root, tip_regions, use_node_labels)
    if n_tips == n_unassigned:
        raise ValueError("没有任何 tip 能匹配区划表（检查 tip 名是否一致）")

    clades = find_regional_clades(root)
    regions = sorted({r for r in tip_regions.values() if r} |
                     {c["region"] for c in clades})

    # ── importations / local_transmissions（1/group_size 规则）──
    imp: Dict[str, float] = {r: 0.0 for r in regions}
    loc: Dict[str, float] = {r: 0.0 for r in regions}
    n_clade_by_region: Dict[str, int] = {r: 0 for r in regions}
    per_tip: Dict[str, Dict[str, float]] = {}
    for c in clades:
        g = c["size"]
        imp[c["region"]] += 1.0          # Σ_tips 1/g = 1
        loc[c["region"]] += g - 1.0      # Σ_tips (1-1/g) = g-1
        n_clade_by_region[c["region"]] += 1
        for t in c["tips"]:
            per_tip[t] = {"region": c["region"], "importation": 1.0 / g,
                          "local_transmission": 1.0 - 1.0 / g, "group_size": g}

    # ── corridors：① 原始边计数 ② 按 1/group_size 的 tip 加权 ──
    edges: Counter = Counter()
    stack = [root]
    while stack:
        n = stack.pop()
        for c in n.children:
            if n.majority and c.majority and n.majority != c.majority:
                edges[(n.majority, c.majority)] += 1
            stack.append(c)

    clade_root_of: Dict[int, Dict] = {}
    for c in clades:
        clade_root_of[c["node_id"]] = c
    node_by_id = {}
    stack = [root]
    while stack:
        n = stack.pop()
        node_by_id[n.id] = n
        stack.extend(n.children)

    weighted: Counter = Counter()
    for c in clades:
        node = node_by_id[c["node_id"]]
        origin = None
        p = node.parent
        while p is not None:
            if p.majority and p.majority != c["region"]:
                origin = p.majority
                break
            p = p.parent
        if origin is None:
            continue  # 祖先区划与自身相同或全未指派 → 不记走廊（如根即本区划）
        weighted[(origin, c["region"])] += 1.0  # Σ_tips 1/g = 1

    exp: Dict[str, float] = {r: 0.0 for r in regions}
    for (a, _b), w in weighted.items():
        exp[a] = exp.get(a, 0.0) + w

    table = []
    for r in regions:
        n_r_tips = sum(1 for t, reg in tip_regions.items() if reg == r)
        table.append(OrderedDict([
            ("region", r),
            ("n_tips", n_r_tips),
            ("n_clades", n_clade_by_region.get(r, 0)),
            ("importations", round(imp.get(r, 0.0), 6)),
            ("local_transmissions", round(loc.get(r, 0.0), 6)),
            ("exportations", round(exp.get(r, 0.0), 6)),
            ("net_flow", round(imp.get(r, 0.0) - exp.get(r, 0.0), 6)),
        ]))

    corridors = [OrderedDict([
        ("from_region", a), ("to_region", b),
        ("n_edges", int(edges.get((a, b), 0))),
        ("tip_weighted", round(weighted.get((a, b), 0.0), 6)),
    ]) for (a, b) in sorted(set(edges) | set(weighted))]

    return {
        "regions": regions,
        "table": table,
        "corridors": corridors,
        "clades": clades,
        "per_tip": per_tip,
        "n_tips": n_tips,
        "n_tips_unassigned": n_unassigned,
        "n_clades_total": len(clades),
        "total_importations": round(sum(imp.values()), 6),
        "total_local_transmissions": round(sum(loc.values()), 6),
        "total_exportations": round(sum(exp.values()), 6),
    }


def timeseries(result: Dict, dates: Dict[str, float], bins: str = "year") -> List[Dict]:
    """按时间分箱汇总 importation / local transmission 权重（SARS 原文的周切片推广）。"""
    per_tip = result["per_tip"]
    rows: Dict[float, Dict[str, float]] = {}
    for tip, w in per_tip.items():
        yr = dates.get(tip)
        if yr is None:
            continue
        key = float(int(yr)) if bins == "year" else round(yr * 12) / 12.0
        slot = rows.setdefault(key, {"importations": 0.0, "local_transmissions": 0.0, "n_tips": 0})
        slot["importations"] += w["importation"]
        slot["local_transmissions"] += w["local_transmission"]
        slot["n_tips"] += 1
    out = []
    for key in sorted(rows):
        slot = rows[key]
        total = slot["importations"] + slot["local_transmissions"]
        out.append(OrderedDict([
            ("time_bin", key),
            ("n_tips", slot["n_tips"]),
            ("importations", round(slot["importations"], 6)),
            ("local_transmissions", round(slot["local_transmissions"], 6)),
            ("frac_importation", round(slot["importations"] / total, 6) if total else 0.0),
        ]))
    return out


# ═══════════════════════════════════════════════════════════════════
# 4. 写出
# ═══════════════════════════════════════════════════════════════════
def _write_tsv(rows: Sequence[Dict], path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if not rows:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write("")
        return path
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return path


def write_outputs(result: Dict, outdir: str, dates: Optional[Dict[str, float]] = None,
                  tree_path: str = "", region_table: str = "") -> List[str]:
    """落盘 import_export.tsv / corridors.tsv / region_clades.tsv (+ 时序 + summary)。"""
    os.makedirs(outdir, exist_ok=True)
    files = []
    files.append(_write_tsv(result["table"], os.path.join(outdir, "import_export.tsv")))
    files.append(_write_tsv(result["corridors"], os.path.join(outdir, "corridors.tsv")))
    clade_rows = [OrderedDict([
        ("clade_index", i + 1),
        ("region", c["region"]),
        ("size", c["size"]),
        ("weight_per_tip", round(c["weight_per_tip"], 6)),
        ("importations", 1.0),
        ("local_transmissions", float(c["size"] - 1)),
        ("clade_root", c["clade_root"]),
        ("tips", ";".join(c["tips"])),
    ]) for i, c in enumerate(result["clades"])]
    files.append(_write_tsv(clade_rows, os.path.join(outdir, "region_clades.tsv")))
    if dates:
        ts = timeseries(result, dates)
        if ts:
            files.append(_write_tsv(ts, os.path.join(outdir, "import_export_timeseries.tsv")))
    summary = {
        "module": "utils/import_export.py",
        "algorithm_source": ("SARS-CoV-2_phylogeo docs/Analyses/R_functions/"
                             "Relative_load_importation_local_transmission.R:15-30 "
                             "(importations += 1/group_size; "
                             "local_transmissions += 1 - 1/group_size)"),
        "group_definition": ("区划簇 = 最大同区划单系枝（子树内 tip 同区划、父节点已跨区划）；"
                             "每簇 = 1 次引入 + (g-1) 次本地传播"),
        "majority_rule": "子代 tip 多数表决；平票取子枝 tip 数最多者，再平票取名字序（确定性）",
        "corridor_definition": {
            "n_edges": "树上父节点多数区划 ≠ 子节点多数区划的边数",
            "tip_weighted": "每簇的 1 次迁出按 1/group_size 摊到 tip（与 importations 同规则）",
        },
        "inputs": {"tree": os.path.abspath(tree_path) if tree_path else None,
                   "region_table": os.path.abspath(region_table) if region_table else None},
        "n_tips": result["n_tips"],
        "n_tips_unassigned": result["n_tips_unassigned"],
        "n_clades_total": result["n_clades_total"],
        "totals": {
            "importations": result["total_importations"],
            "local_transmissions": result["total_local_transmissions"],
            "exportations": result["total_exportations"],
        },
        "files": files,
    }
    sp = os.path.join(outdir, "import_export_summary.json")
    with open(sp, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    files.append(sp)
    return files


# ═══════════════════════════════════════════════════════════════════
# 5. 玩具树 + 手算真值自测
# ═══════════════════════════════════════════════════════════════════
# 树 1（含跨区划入/出，手算真值见下方断言）：
#   (((A1:1,A2:1)x:1,B1:1)y:1,(B2:1,C1:1)z:1)root;
#   区划簇：x(A,2) / B1(B,1) / B2(B,1) / C1(C,1)  → 共 4 簇
#   importations: A=1, B=2, C=1 ；local: A=1, B=0, C=0 ；n_tips=5
#   多数区划：x→A, y→A(2:1), z→B(平票取子枝 B2), root→A(平票取大子枝 y)
#   corridors: A→B edges=2,tip_weighted=2 ；B→C edges=1,tip_weighted=1
#   exportations: A=2, B=1, C=0
TOY_TREE_1 = "(((A1:1,A2:1)x:1,B1:1)y:1,(B2:1,C1:1)z:1)root;"
TOY_REGIONS_1 = {"A1": "A", "A2": "A", "B1": "B", "B2": "B", "C1": "C"}

# 树 2（两个完整单系区划簇，验证 1/group_size 的非平凡权重）：
#   ((A1:1,A2:1,A3:1)x:1,(B1:1,B2:1)y:1)root;
#   A: 1 簇 size3 → import 1 / local 2 ；B: 1 簇 size2 → import 1 / local 1
#   root 多数 = A(3:2) → corridors A→B edges=1, tip_weighted=1
TOY_TREE_2 = "((A1:1,A2:1,A3:1)x:1,(B1:1,B2:1)y:1)root;"
TOY_REGIONS_2 = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B"}

# 树 3：带注释/引号/根标记/内部节点状态标签的解析测试
TOY_TREE_3 = "[&R](('Hubei 1':0.1,Hubei2:0.2)Ningxia[&posterior=0.9]:0.3,Yunnan1:0.4)root;"


def _selftest(verbose: bool = True) -> int:
    n_ok = 0

    def _check(cond: bool, msg: str) -> None:
        nonlocal n_ok
        if not cond:
            raise AssertionError(f"SELFTEST FAILED: {msg}")
        n_ok += 1
        if verbose:
            print(f"  [ok {n_ok:2d}] {msg}")

    # ── A. Newick 解析 ─────────────────────────────────────────────
    if verbose:
        print("A. Newick 解析（注释 / 引号 / 分支长 / 根标记）")
    r1 = parse_newick(TOY_TREE_1)
    tips1 = [t.name for t in iter_tips(r1)]
    _check(tips1 == ["A1", "A2", "B1", "B2", "C1"],
           f"树 1 的 5 个 tip 顺序解析正确：{tips1}")
    _check(r1.n_tips == 5 and r1.children[0].n_tips == 3,
           "内部节点 tip 计数正确（root=5, y 子枝=3）")
    _check(abs(r1.children[0].children[1].length - 1.0) < 1e-12, "分支长度解析正确（:1）")
    r3 = parse_newick(TOY_TREE_3)
    tips3 = [t.name for t in iter_tips(r3)]
    _check(tips3 == ["Hubei 1", "Hubei2", "Yunnan1"],
           f"引号内空格标签 + [注释] + [&R] 根标记解析正确：{tips3}")
    _check(node_state_label("Ningxia", ["Ningxia", "Yunnan"]) == "Ningxia",
           "内部节点状态标签匹配已知区划")
    _check(node_state_label("0.98", ["Ningxia"]) is None, "纯数值标签不被误判为区划")

    # ── B. 手算真值：importations / local_transmissions ────────────
    if verbose:
        print("B. 1/group_size 规则（树 1，手算真值对拍）")
    res1 = summarize_import_export(r1, TOY_REGIONS_1)
    t1 = {row["region"]: row for row in res1["table"]}
    _check(t1["A"]["importations"] == 1.0 and t1["B"]["importations"] == 2.0
           and t1["C"]["importations"] == 1.0,
           f"importations 手算值 A=1,B=2,C=1 → 实得 "
           f"{t1['A']['importations']},{t1['B']['importations']},{t1['C']['importations']}")
    _check(t1["A"]["local_transmissions"] == 1.0 and t1["B"]["local_transmissions"] == 0.0
           and t1["C"]["local_transmissions"] == 0.0,
           "local_transmissions 手算值 A=1,B=0,C=0（簇 size2 → 2-1=1）")
    _check(t1["A"]["n_clades"] == 1 and t1["B"]["n_clades"] == 2 and t1["C"]["n_clades"] == 1,
           "区划簇数 A=1,B=2,C=1（toy tree 上手工可数）")
    _check(res1["n_clades_total"] == 4 and res1["total_importations"] == 4.0,
           "总引入数 == 区划簇总数 == 4（1/group_size 求和为 1 的守恒性）")
    _check(abs(res1["total_importations"] + res1["total_local_transmissions"] - res1["n_tips"]) < 1e-9,
           "Σ(importations + local_transmissions) == n_tips == 5（每 tip 权重恰为 1）")
    # 每 tip 权重 = 1/g
    _check(abs(res1["per_tip"]["A1"]["importation"] - 0.5) < 1e-12
           and abs(res1["per_tip"]["B1"]["importation"] - 1.0) < 1e-12,
           "每 tip 引入权重 = 1/group_size（A1=0.5、B1=1.0）")
    _check(abs(res1["per_tip"]["A1"]["local_transmission"] - 0.5) < 1e-12
           and res1["per_tip"]["B1"]["local_transmission"] == 0.0,
           "每 tip 本地权重 = 1 - 1/group_size（A1=0.5、B1=0.0）")
    # 树 2：非平凡簇大小
    res2 = summarize_import_export(parse_newick(TOY_TREE_2), TOY_REGIONS_2)
    t2 = {row["region"]: row for row in res2["table"]}
    _check(t2["A"]["importations"] == 1.0 and t2["A"]["local_transmissions"] == 2.0
           and t2["B"]["importations"] == 1.0 and t2["B"]["local_transmissions"] == 1.0,
           "树 2（簇 size 3 / 2）：A 引入 1 本地 2；B 引入 1 本地 1 —— 手算对拍")
    _check(abs(res2["per_tip"]["A3"]["importation"] - 1.0 / 3) < 1e-12,
           "size=3 的簇每 tip 引入权重 = 1/3（SARS 注释 'group size is 30 → 1/30'）")

    # ── C. 走廊 / 输出流向 ────────────────────────────────────────
    if verbose:
        print("C. corridors / exportations（手算真值对拍）")
    c1 = {(c["from_region"], c["to_region"]): c for c in res1["corridors"]}
    _check(c1[("A", "B")]["n_edges"] == 2, "走廊 A→B 原始边计数 = 2（y→B1, root→z）")
    _check(c1[("A", "B")]["tip_weighted"] == 2.0, "走廊 A→B 加权 = 2.0（两条 size=1 的簇）")
    _check(c1[("B", "C")]["n_edges"] == 1 and c1[("B", "C")]["tip_weighted"] == 1.0,
           "走廊 B→C = 1（z→C1）")
    _check(t1["A"]["exportations"] == 2.0 and t1["B"]["exportations"] == 1.0
           and t1["C"]["exportations"] == 0.0,
           "exportations 手算值 A=2,B=1,C=0")
    _check(abs(sum(c["tip_weighted"] for c in res1["corridors"]) - res1["total_exportations"]) < 1e-9,
           "Σ 加权走廊 == Σ exportations（两种口径自洽）")
    _check(t1["B"]["net_flow"] == 1.0 and t1["C"]["net_flow"] == 1.0,
           "净流量 net = importations - exportations（B=+1, C=+1）")
    c2 = {(c["from_region"], c["to_region"]): c for c in res2["corridors"]}
    _check(c2[("A", "B")]["n_edges"] == 1 and c2[("A", "B")]["tip_weighted"] == 1.0,
           "树 2 走廊 A→B = 1（y 整个簇 1 次迁出）")

    # ── D. 边角：未指派 tip / 时间序列 / 输出文件 ─────────────────
    if verbose:
        print("D. 边角情形（未指派 tip / 时间序列 / 落盘）")
    part = {"A1": "A", "A2": "A", "B1": "B"}  # B2, C1 未指派
    res_p = summarize_import_export(parse_newick(TOY_TREE_1), part)
    _check(res_p["n_tips_unassigned"] == 2, "未指派 tip 计数 = 2（如实记录，不静默丢弃）")
    _check(res_p["n_tips"] == 5, "n_tips 仍报全树 tip 数（分母透明）")
    dates = {"A1": 2015.0, "A2": 2016.5, "B1": 2016.0, "B2": 2017.25, "C1": 2015.5}
    ts = timeseries(res1, dates)
    _check(len(ts) == 3 and abs(sum(r["importations"] for r in ts) - 4.0) < 1e-9,
           "时间序列按年分箱：3 个箱，Σ 引入权重守恒 = 4.0")
    _check(abs(sum(r["frac_importation"] for r in ts) - sum(
        r["importations"] / (r["importations"] + r["local_transmissions"]) for r in ts)) < 1e-9,
        "时序表 frac_importation 自洽（引入占比）")
    # 2026-09-16 (D11 修复): 改为委托 utils.decimal_year 后, 年份语义由
    # "年初" 变为统一口径 "年中"; 断言同步更新。对 int() 分箱无影响 (见函数 docstring)。
    _check(_to_decimal_year("2016-07-15") is not None
           and abs(_to_decimal_year("2016") - 2016.4556) < 0.01,
           "日期解析：'2016' → 2016.456（统一口径年中）")
    _check(abs(_to_decimal_year("2016-07-02") - 2016.5) < 0.01, "'2016-07-02' ≈ 2016.5")
    # 保留本模块特有格式: DD/MM/YYYY 与 DD.MM.YYYY (统一实现不认, 需前置归一化)
    _check(abs(_to_decimal_year("15/07/2016") - 2016.5376) < 0.01,
           "'15/07/2016' (DD/MM/YYYY) 正确解析")
    _check(abs(_to_decimal_year("15.07.2016") - 2016.5376) < 0.01,
           "'15.07.2016' (DD.MM.YYYY) 正确解析")
    # 不再静默编造: 乱串必须返回 None (旧版会给 1999.5)
    _check(_to_decimal_year("xxxx1999yy") is None,
           "无法解析的串返回 None（不凭正则猜年份）")
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        files = write_outputs(res1, td, dates, "toy.nwk", "toy.tsv")
        names = sorted(os.path.basename(f) for f in files)
        _check(names == ["corridors.tsv", "import_export.tsv", "import_export_summary.json",
                         "import_export_timeseries.tsv", "region_clades.tsv"],
               f"5 个输出文件齐备：{names}")
        with open(os.path.join(td, "import_export.tsv"), encoding="utf-8") as fh:
            head = fh.readline().strip().split("\t")
            body = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
        _check(head == ["region", "n_tips", "n_clades", "importations",
                        "local_transmissions", "exportations", "net_flow"],
               "import_export.tsv 表头契约（polio 14 表中 import_export.csv 的落地版）")
        _check(len(body) == 3 and {r[0] for r in body} == {"A", "B", "C"},
               "import_export.tsv 含 3 个区划行")
        with open(os.path.join(td, "import_export_summary.json"), encoding="utf-8") as fh:
            summ = json.load(fh)
        _check("Relative_load_importation_local_transmission.R:15-30"
               in summ["algorithm_source"], "summary.json 记录算法出处（含行号）")

    # ── E. 与 biopython 解析对拍（若可用）─────────────────────────
    if verbose:
        print("E. 与 biopython Newick 解析对拍（可选依赖）")
    try:
        from io import StringIO
        from Bio import Phylo as _Phylo
        bp = [t.name for t in _Phylo.read(StringIO(TOY_TREE_2), "newick").get_terminals()]
        _check(bp == [t.name for t in iter_tips(parse_newick(TOY_TREE_2))],
               "自写解析器与 Bio.Phylo 的 tip 集合/顺序一致（树 2）")
    except ImportError:
        _check(True, "biopython 不可用 → 跳过对拍（不影响功能）")

    if verbose:
        print(f"\nSELFTEST PASSED — {n_ok} 项断言全部通过")
    return n_ok


# ═══════════════════════════════════════════════════════════════════
# 6. CLI
# ═══════════════════════════════════════════════════════════════════
def _build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="import_export",
        description="引入 / 本地传播 / 输出拆分（SARS 1/group_size 规则）+ corridors 表",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("用法\n----")[-1][:1000])
    ap.add_argument("--selftest", action="store_true",
                    help="跑玩具树自测（手算真值断言）后退出")
    ap.add_argument("--tree", help="定年树（Newick；BEAST MCC / IQ-TREE 均可）")
    ap.add_argument("--regions", help="元数据表 CSV/TSV（tip, region [, date]）")
    ap.add_argument("--outdir", default="import_export_out", help="输出目录")
    ap.add_argument("--use-node-labels", action="store_true",
                    help="内部节点标签若为已知区划名，优先用作该节点状态（走廊口径）")
    ap.add_argument("--bins", default="year", choices=["year", "month"],
                    help="时序分箱粒度（有日期列时生效）")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)
    if args.selftest:
        print("=" * 64)
        print("import_export 自测（玩具树 + 手算真值；算法出处 SARS R 脚本 15-30 行）")
        print("=" * 64)
        try:
            _selftest()
        except AssertionError as e:
            print(f"\n{e}")
            return 1
        return 0

    if not args.tree or not args.regions:
        print("错误：需要 --tree 与 --regions（或 --selftest）", file=sys.stderr)
        return 2

    with open(args.tree, encoding="utf-8") as fh:
        newick = fh.read()
    root = parse_newick(newick)
    tip_regions, dates = load_region_table(args.regions)
    if not tip_regions:
        print(f"错误：{args.regions} 未解析到任何 tip→区划", file=sys.stderr)
        return 2

    result = summarize_import_export(root, tip_regions, args.use_node_labels)
    if args.bins == "month":
        dates = {k: v for k, v in dates.items()}
    files = write_outputs(result, args.outdir, dates, args.tree, args.regions)

    print(f"tips: {result['n_tips']}（未指派 {result['n_tips_unassigned']}）  "
          f"区划簇: {result['n_clades_total']}  区划: {len(result['regions'])}")
    print(f"Σ importations={result['total_importations']}  "
          f"Σ local_transmissions={result['total_local_transmissions']}  "
          f"Σ exportations={result['total_exportations']}")
    print(f"\n{'region':<16}{'n_tips':>7}{'clades':>7}{'import':>9}{'local':>9}{'export':>9}{'net':>8}")
    for row in result["table"]:
        print(f"{row['region']:<16}{row['n_tips']:>7}{row['n_clades']:>7}"
              f"{row['importations']:>9}{row['local_transmissions']:>9}"
              f"{row['exportations']:>9}{row['net_flow']:>8}")
    if result["corridors"]:
        print("\ncorridors (from → to)：")
        for c in result["corridors"]:
            print(f"  {c['from_region']} → {c['to_region']}: "
                  f"n_edges={c['n_edges']}, tip_weighted={c['tip_weighted']}")
    for f in files:
        print(f"[OUT ] {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
