#!/usr/bin/env python3
"""
phylo_cross.py — 交叉验证工具 (参考 YR-MPE 能力, 规范化对照)

1. run_raxml_ng : RAxML-NG ML 建树 (与 IQ-TREE 交叉验证)
2. run_mrbayes  : MrBayes 贝叶斯建树 (工具缺失时优雅跳过)
3. run_phylobayes: PhyloBayes-MPI (工具缺失时优雅跳过)
4. run_lsd2     : LSD2 最小二乘定年 (工具缺失时优雅跳过, BEAST 定年交叉验证)
5. export_icytree: 树 → NEXUS (IcyTree 可打开)

用法:
  python -m utils.phylo_cross raxml --aln mafft.aln.fasta --outdir xv_out
  python -m utils.phylo_cross lsd2 --tree iqtree.treefile --meta dates.csv
  python -m utils.phylo_cross icytree --tree iqtree.treefile --out icytree.nex
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


def _which(tool):
    return shutil.which(tool)


_TREE_STMT_RE = re.compile(r'\btree\s+\S+\s*=\s*(.+?);', re.S)


def _extract_tree_stmt(content):
    """从 NEXUS / Newick 文本里取 (树语句, Translate 块)。

    历史坑: 旧实现只认「以 `(` 开头、以 `;` 结尾」的单行；BEAST2 的 .tre 里树语句是
    `tree TREE1 = [&R] (…);` → 一行都匹配不上，最后落到 `content` 兜底，把**整个
    NEXUS 文件**当成 newick 写进输出（还返回 success）——IcyTree 打不开，且看不出来。
    另外 BEAST2 默认带 Translate 块把叶名缩写成 1..N，丢了它叶名就成数字。

    返回值里的 newick **不带**结尾分号（由调用方补一个），两条路径口径必须一致：
    历史坑: 正则路径剥了分号、兜底路径留着分号，写出去变成 `TREE tree = (…)\\nEND;`
    ——NEXUS 阅读器按 ';' 切命令，`END` 被并进树命令，整段 TREES 块解析出 0 棵树
    （Bio.Nexus 读回 trees=0 translate=0，且不报错）。
    """
    m = _TREE_STMT_RE.search(content)
    if m:
        tr = re.search(r'\bTranslate\b(.*?);', content, re.S | re.I)
        return m.group(1).strip().rstrip(';').strip(), \
            ('Translate' + tr.group(1) + ';') if tr else ''
    for line in content.splitlines():          # 纯 Newick（可能带 [&R] 前缀）
        line = line.strip()
        if line.startswith('(') and line.endswith(';'):
            return line[:-1].strip(), ''
    txt = content.strip()
    if txt.startswith('(') and txt.endswith(';'):
        return txt[:-1].strip(), ''
    return '', ''


def run_raxml_ng(aln, outdir="xv_out", threads=8, model="GTR+G", log=None):
    """RAxML-NG ML 建树 (与 IQ-TREE 交叉验证)"""
    def emit(m):
        if log:
            log(m)
    bin_ = _which("raxml-ng")
    if not bin_:
        return {"success": False, "error": "raxml-ng 未安装"}
    os.makedirs(outdir, exist_ok=True)
    prefix = os.path.join(outdir, "raxml")
    emit(f"[raxml-ng] 建树 (model={model}, threads={threads})...")
    cmd = [bin_, "--all", "--msa", aln, "--model", model, "--prefix", prefix,
           "--threads", str(threads), "--bs-trees", "100", "--seed", "42"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=7200)
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "超时"}
    tree = prefix + ".raxml.bestTree"
    if r.returncode == 0 and os.path.exists(tree):
        emit(f"[raxml-ng] → {tree}")
        return {"success": True, "tree": tree,
                "bootstrap": prefix + ".raxml.support",
                "log": prefix + ".raxml.log"}
    return {"success": False, "error": (r.stderr or r.stdout)[-300:]}


def run_mrbayes(aln, outdir="mrbayes_out", chain=1000000, log=None):
    """MrBayes 贝叶斯建树 (工具缺失时优雅跳过)"""
    def emit(m):
        if log:
            log(m)
    bin_ = _which("mb") or _which("mrbayes")
    if not bin_:
        return {"success": False, "error": "MrBayes 未安装 (conda install -c bioconda mrbayes)"}
    os.makedirs(outdir, exist_ok=True)
    nexus = os.path.join(outdir, "input.nex")
    # 转 NEXUS
    try:
        from Bio import SeqIO
        with open(nexus, "w") as f:
            SeqIO.write(SeqIO.parse(aln, "fasta"), f, "nexus")
    except Exception as e:
        return {"success": False, "error": f"NEXUS 转换失败: {e}"}
    with open(nexus, "a") as f:
        f.write(f"\nbegin mrbayes;\n  lset nst=6 rates=invgamma;\n  mcmc ngen={chain} nruns=2 nchains=4 samplefreq=100;\n  sump burnin=25%;\n  sumt burnin=25%;\nend;\n")
    emit(f"[mrbayes] 运行 (chain={chain})...")
    try:
        r = subprocess.run([bin_, nexus], capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=7200, cwd=outdir)
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "超时"}
    if r.returncode == 0:
        return {"success": True, "nexus": nexus, "output": outdir}
    return {"success": False, "error": (r.stderr or r.stdout)[-300:]}


def run_phylobayes(aln, outdir="phylobayes_out", chain=1000, log=None):
    """PhyloBayes-MPI (工具缺失时优雅跳过)"""
    def emit(m):
        if log:
            log(m)
    bin_ = _which("pb_mpi") or _which("phylobayes")
    if not bin_:
        return {"success": False, "error": "PhyloBayes 未安装"}
    os.makedirs(outdir, exist_ok=True)
    emit(f"[phylobayes] 运行 (chain={chain})...")
    try:
        r = subprocess.run([bin_, "-d", aln, "-cat", "-gtr"], capture_output=True,
                           text=True, encoding='utf-8', errors='replace', timeout=7200, cwd=outdir)
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "超时"}
    if r.returncode == 0:
        return {"success": True, "output": outdir}
    return {"success": False, "error": (r.stderr or r.stdout)[-300:]}


def run_lsd2(tree, meta_csv, outdir="lsd2_out", log=None):
    """LSD2 最小二乘定年 (BEAST 定年的快速交叉验证; 工具缺失时优雅跳过)"""
    def emit(m):
        if log:
            log(m)
    bin_ = _which("lsd2")
    if not bin_:
        return {"success": False, "error": "LSD2 未安装 (https://github.com/tothuhien/lsd2)"}
    os.makedirs(outdir, exist_ok=True)
    # 提取 dates (name,date)
    dates_file = os.path.join(outdir, "dates.txt")
    import csv as _csv
    rows = []
    with open(meta_csv, encoding="utf-8-sig") as f:
        for row in _csv.DictReader(f):
            name = row.get("name") or row.get("Name")
            date = row.get("date") or row.get("Date")
            if name and date:
                y = date[:4]
                if y.isdigit():
                    rows.append(f"{name} {y}")
    with open(dates_file, "w") as f:
        f.write("\n".join(rows) + "\n")
    emit(f"[lsd2] 定年 ({len(rows)} 样本日期)...")
    try:
        r = subprocess.run([bin_, "-i", tree, "-d", dates_file, "-o", outdir, "-c", "-r", "a"],
                           capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=3600)
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "超时"}
    if r.returncode == 0:
        return {"success": True, "output": outdir}
    return {"success": False, "error": (r.stderr or r.stdout)[-300:]}


def compare_topology(newick1: str, newick2: str):
    """两树拓扑一致性 (clade 集合比较, 用 tree_ops)
    返回 {n_clades1, n_clades2, shared, pct, rf_like}。"""
    from utils.tree_ops import parse_newick, get_leaf_names

    def clades(newick):
        root = parse_newick(newick)
        total = set(get_leaf_names(root))
        result = set()

        def walk(n):
            if not n.children:
                return
            leaves = set(get_leaf_names(n))
            if 1 < len(leaves) < len(total):
                result.add(frozenset(leaves))
            for c in n.children:
                walk(c)
        walk(root)
        return result

    c1 = clades(newick1)
    c2 = clades(newick2)
    shared = len(c1 & c2)
    union = len(c1 | c2)
    pct = 100.0 * shared / union if union else 0.0
    # RF 距离 = 两树独有的 clade 数 (对称差 / 2 ≈ 拆分距离)
    rf = (len(c1) + len(c2) - 2 * shared) / 2.0
    return {"n_clades1": len(c1), "n_clades2": len(c2),
            "shared": shared, "pct": round(pct, 1), "rf_like": rf}


def export_icytree(tree, out_path=None, log=None):
    """树 → NEXUS (IcyTree 浏览器可打开: https://icytree.org)"""
    def emit(m):
        if log:
            log(m)
    content = Path(tree).read_text(encoding="utf-8", errors="replace").strip()
    newick, tr_block = _extract_tree_stmt(content)
    if not newick:
        return {"success": False,
                "error": "无法从树文件提取 newick（既没有 `tree NAME = …;` 语句，"
                         "也不是 `(…);` 形式的纯 Newick）：%s" % tree}
    out_path = out_path or (str(Path(tree)) + ".icytree.nex")
    # newline="\n": Windows 上默认会把 \n 写成 \r\n；NEXUS/Newick 是要跨平台
    # 交换给 IcyTree / Bio.Nexus 的，别掺 CR。
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("#NEXUS\nBEGIN TREES;\n")
        if tr_block:                       # BEAST2 的数字叶名靠 Translate 还原
            f.write("  " + tr_block + "\n")
        f.write(f"  TREE tree = {newick};\nEND;\n")
    emit(f"[icytree] → {out_path} (用 https://icytree.org 打开)")
    return {"success": True, "nexus": out_path}


def main():
    ap = argparse.ArgumentParser(description="交叉验证工具 (RAxML-NG/MrBayes/PhyloBayes/LSD2/IcyTree)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("raxml", help="RAxML-NG ML 建树")
    p1.add_argument("--aln", required=True)
    p1.add_argument("--outdir", default="xv_out")
    p1.add_argument("--threads", type=int, default=8)
    p1.add_argument("--model", default="GTR+G")
    p1.add_argument("--ref-tree", default=None, help="参考树 (如 IQ-TREE 树) → 拓扑比较")

    p2 = sub.add_parser("mrbayes", help="MrBayes 贝叶斯")
    p2.add_argument("--aln", required=True)
    p2.add_argument("--chain", type=int, default=1000000)

    p3 = sub.add_parser("phylobayes", help="PhyloBayes")
    p3.add_argument("--aln", required=True)

    p4 = sub.add_parser("lsd2", help="LSD2 定年")
    p4.add_argument("--tree", required=True)
    p4.add_argument("--meta", required=True)

    p5 = sub.add_parser("icytree", help="树 → IcyTree NEXUS")
    p5.add_argument("--tree", required=True)
    p5.add_argument("--out", default=None)

    args = ap.parse_args()
    if args.cmd == "raxml":
        r = run_raxml_ng(args.aln, args.outdir, args.threads, args.model, log=print)
        # 拓扑比较: RAxML-NG 树 vs 参考树 (如 IQ-TREE)
        if r.get("success") and args.ref_tree and os.path.exists(args.ref_tree):
            try:
                ref_newick = open(args.ref_tree).read().strip()
                alt_newick = open(r["tree"]).read().strip()
                cmp = compare_topology(ref_newick, alt_newick)
                print(f"[拓扑比较] 参考树 clades={cmp['n_clades1']}, RAxML-NG clades={cmp['n_clades2']}")
                print(f"[拓扑比较] 共享 clade: {cmp['shared']} ({cmp['pct']}%), RF≈{cmp['rf_like']}")
                r["topology"] = cmp
            except Exception as e:
                print(f"[拓扑比较] 跳过: {e}")
    elif args.cmd == "mrbayes":
        r = run_mrbayes(args.aln, "mrbayes_out", args.chain, log=print)
    elif args.cmd == "phylobayes":
        r = run_phylobayes(args.aln, log=print)
    elif args.cmd == "lsd2":
        r = run_lsd2(args.tree, args.meta, log=print)
    elif args.cmd == "icytree":
        r = export_icytree(args.tree, args.out, log=print)
    if r.get("success"):
        print(f"✓ {args.cmd} 完成")
        return 0
    print(f"✗ {args.cmd}: {r.get('error')}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
