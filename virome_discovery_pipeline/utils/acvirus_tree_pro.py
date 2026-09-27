#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ACVirus Tree-Pro (Streamlined Edition) — 按科建树 + 属级分色 + 共线性双面板

v2 新增 (移植自 vfam_trees, 离线版):
  - LCA 内部节点分类注释 (基于 VMR taxa 表的 ranked lineage, 无需联网)
  - 分类学引导定根 (MAD / 中点法兜底)
  - 属级 HLS 配色 (亚科色相带, colorsys 标准库)
  - 输出: phylogeny.annotated.nwk / phylogeny.phyloxml (rooted, rerootable=false)
          / tree_auspice.json (Nextstrain v2, 可拖入 auspice.us) / tree_metadata.tsv
  无新增外部依赖 (全部 python 标准库 + 既有 biopython)。
"""
from __future__ import annotations
import argparse
import colorsys
import copy
import io
import json
import os
import re
import shlex
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from matplotlib.colors import LinearSegmentedColormap, Normalize

COMMAND_LOG_NAME = "tree_pro_commands.jsonl"
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

def find_executable(names): 
    for name in names:
        path = shutil.which(name)
        if path: return path
    return None

def run_cmd(cmd, label, log_dir, stdout_file=None):
    cmd = [str(x) for x in cmd]
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / COMMAND_LOG_NAME
    start = datetime.now().isoformat(timespec="seconds")
    stdout_handle = stdout_file.open("w") if stdout_file else subprocess.PIPE
    print(f"[RUN] {label} -> {' '.join(cmd[:6])}{' ...' if len(cmd) > 6 else ''}")
    res = subprocess.run(cmd, stdout=stdout_handle, stderr=subprocess.PIPE, text=True)
    if stdout_file: stdout_handle.close()
    with log_file.open("a", encoding="utf-8") as h:
        h.write(json.dumps({"label": label, "cmd": shlex.join(cmd), "start": start,
                            "returncode": res.returncode, "stderr": res.stderr}) + "\n")
    if res.returncode != 0:
        print(f"[ERROR] {label} failed:\n{res.stderr}", file=sys.stderr)
        raise subprocess.CalledProcessError(res.returncode, cmd, stderr=res.stderr)
    return res

def get_db_lengths(db_fasta, outdir):
    """seqkit fx2tab 扫参考库 → {基础acc(无版本号): 序列长度}; 结果缓存到 outdir/db_lengths.tsv"""
    cache = outdir / "db_lengths.tsv"
    if not cache.exists():
        seqkit = find_executable(["seqkit"])
        if not seqkit:
            raise FileNotFoundError("seqkit required!")
        run_cmd([seqkit, "fx2tab", "-i", "-l", str(db_fasta)], "seqkit fx2tab (lengths)",
                outdir, stdout_file=cache)
    lens = {}
    with cache.open() as h:
        for line in h:
            parts = line.rstrip("\n").rsplit("\t", 1)
            if len(parts) != 2:
                continue
            try:
                L = float(parts[1])
            except ValueError:
                continue
            base = parts[0].split()[0].split(".")[0].strip()
            if base and base not in lens:
                lens[base] = L
    return lens


def get_hierarchical_accessions(df_taxa, target_rank, target_name, mode,
                                target_ratio=1.0, other_count=3,
                                acc_lengths=None, length_ratio=0.0,
                                length_ratio_fallback=0.0):
    """同科内智能抽样 (macro/genus/lineage), 防止全量抓取;
    可选抽样前长度预筛 (两段式): 先按科内属中位数 ±length_ratio 严格窗筛,
    属内全越界的属放宽 ±length_ratio_fallback 兜底窗重试, 仍空则保底最接近中位数 1 条"""
    df = df_taxa.copy()
    for col in ['Realm','Order','Class','Family','Genus','Species']:
        if col in df.columns:
            df[col] = df[col].fillna("").astype(str).str.strip()
    sampled_accs = set()
    target_info = df[df[target_rank].str.lower() == target_name.lower()]
    if target_info.empty:
        raise ValueError(f"Cannot find {target_rank} '{target_name}' in database taxonomy!")
    t_family = target_info.iloc[0].get('Family','')
    t_genus  = target_info.iloc[0].get('Genus','')
    print(f"\n[SAMPLING] Mode: {mode.upper()} | Target: {target_rank}={target_name} (Family: {t_family})")
    df['Clean_Acc'] = df['Virus GENBANK accession'].apply(lambda x: str(x).split('.')[0].split(';')[0].strip() if pd.notna(x) else '')
    # 长度预筛 (抽样前, 保证各属抽样名额从合格池里出), 两段式:
    #   严格窗 ±length_ratio → 属内全越界的属放宽兜底窗 ±length_ratio_fallback 重试
    #   → 仍空的属保底最接近中位数的 1 条 (中位数取科内该属全部候选)
    if acc_lengths and length_ratio > 0:
        fam_mask = (df['Family'].str.lower() == str(t_family).lower()) if 'Family' in df.columns else pd.Series(True, index=df.index)
        fam_df = df[fam_mask]
        genus_lens = defaultdict(list)
        for g, acc in zip(fam_df['Genus'].astype(str), fam_df['Clean_Acc']):
            L = acc_lengths.get(str(acc))
            if L:
                genus_lens[g].append(L)
        med_of = {g: statistics.median(v) for g, v in genus_lens.items()}
        fb_ratio = (length_ratio_fallback
                    if length_ratio_fallback and length_ratio_fallback > length_ratio else None)

        def _window_pass(ratio, only_genera=None):
            keep_idx, drop_idx = [], []
            for idx, g, acc in zip(fam_df.index, fam_df['Genus'].astype(str), fam_df['Clean_Acc']):
                if only_genera is not None and g not in only_genera:
                    continue
                med = med_of.get(g)
                L = acc_lengths.get(str(acc))
                if not med or not L or (1 - ratio) * med <= L <= (1 + ratio) * med:
                    keep_idx.append(idx)
                else:
                    drop_idx.append(idx)
            return keep_idx, drop_idx

        keep_idx, _ = _window_pass(length_ratio)
        kept_g = set(fam_df.loc[keep_idx, 'Genus'].astype(str))
        empty_g = [g for g in med_of if g not in kept_g]
        n_fallback = n_floor = 0
        if empty_g and fb_ratio:  # 兜底窗: 全越界的属放宽重试
            fidx, _ = _window_pass(fb_ratio, only_genera=set(empty_g))
            keep_idx.extend(fidx)
            n_fallback = len(fidx)
            kept_g = set(fam_df.loc[keep_idx, 'Genus'].astype(str))
            empty_g = [g for g in med_of if g not in kept_g]
        for g in empty_g:  # 兜底仍空 → 每属保底最接近中位数 1 条
            sub = fam_df[fam_df['Genus'].astype(str) == g]
            best_idx, best_d = None, float('inf')
            for idx, acc in zip(sub.index, sub['Clean_Acc']):
                L = acc_lengths.get(str(acc))
                if L:
                    d = abs(L - med_of[g])
                    if d < best_d:
                        best_d, best_idx = d, idx
            if best_idx is not None:
                keep_idx.append(best_idx)
                n_floor += 1
        n_drop = len(fam_df) - len(keep_idx)
        df = df.loc[keep_idx]
        if n_drop:
            msg = (f"[LENGTH] 抽样前长度预筛 ±{length_ratio * 100:g}% (属中位数): "
                   f"剔除 {n_drop} 条, 可选候选 {len(df)} 条")
            if n_fallback:
                msg += f"; {n_fallback} 条经 ±{fb_ratio * 100:g}% 兜底窗保留"
            if n_floor:
                msg += f"; {n_floor} 属全越界保底 1 条"
            print(msg)
    if mode == "macro":
        fam_df = df[df['Family'].str.lower()==t_family.lower()]
        for genus, group in fam_df.groupby('Genus'):
            if genus == "" or pd.isna(genus): continue
            if str(genus).lower()==str(t_genus).lower():
                n = max(1, int(len(group)*target_ratio))
                sampled = group.sample(n=n, random_state=42)
            else:
                sampled = group.sample(n=min(len(group), other_count), random_state=42)
            sampled_accs.update(sampled['Virus GENBANK accession'].tolist())
    elif mode == "genus":
        target_group = df[df['Genus'].str.lower()==str(t_genus).lower()]
        sampled_accs.update(target_group['Virus GENBANK accession'].tolist())
        fam_other = df[(df['Family'].str.lower()==t_family.lower()) & (df['Genus'].str.lower()!=str(t_genus).lower())]
        if not fam_other.empty:
            for _, group in fam_other.groupby('Genus'):
                sampled_accs.update(group.sample(n=min(len(group), other_count), random_state=42)['Virus GENBANK accession'].tolist())
    elif mode == "lineage":
        sp_name = target_name if target_rank=='Species' else target_info.iloc[0].get('Species','')
        target_group = df[df['Species'].str.lower()==str(sp_name).lower()]
        sampled_accs.update(target_group['Virus GENBANK accession'].tolist())
        same_genus_other = df[(df['Genus'].str.lower()==str(t_genus).lower()) & (df['Species'].str.lower()!=str(sp_name).lower())]
        if not same_genus_other.empty:
            for _, group in same_genus_other.groupby('Species'):
                sampled_accs.update(group.sample(n=1, random_state=42)['Virus GENBANK accession'].tolist())
    final_accs = [str(acc).split(';')[0].strip() for acc in sampled_accs if pd.notna(acc) and str(acc).strip()]
    _tax_cols = [c for c in ('Realm','Kingdom','Phylum','Class','Order','Family',
                             'Subfamily','Genus','Species') if c in df.columns]
    taxa_meta = df.set_index('Clean_Acc')[_tax_cols].to_dict('index')
    # 补充: custom 模式外, 每个入库 accession 的 ranked lineage (root→tip)
    taxa_meta = {k: {c: v for c, v in row.items() if pd.notna(v)} for k, row in taxa_meta.items()}
    print(f"  -> Sampled {len(final_accs)} representative genomes from database.")
    return final_accs, taxa_meta

def parse_best_model(model_out_file):
    if not Path(model_out_file).exists(): return None
    with open(model_out_file,'r') as f:
        for line in f:
            if line.strip().startswith("Best model according to"):
                for next_line in f:
                    if next_line.strip().startswith("Model:"):
                        return next_line.strip().split()[1]
    return None

def run_phylogeny_engine(aln_file, outdir, threads, engine_preference="iqtree", bootstrap=1000,
                          rooting="unrest", aln_fasta=None):
    """建树引擎: (默认) IQ-TREE + UNREST 非可逆模型 → 天然有根树 (无外群定根)
       也可用 MFP (ModelFinder 可逆模型, 无根) + 后续定根."""
    tree_file = outdir / "phylogeny.treefile"
    modeltest = find_executable(["modeltest-ng"]); raxml = find_executable(["raxml-ng"])
    iqtree = find_executable(["iqtree2","iqtree","iqtree3"])
    fasttree = find_executable(["FastTree","FastTreeDbl","fasttree"])
    # FastTree 优先: 秒级出树, SH-like 支持值 (0-100) 与 bootstrap 同尺度
    if engine_preference == "fasttree" or (engine_preference == "auto" and fasttree):
        if not fasttree: raise FileNotFoundError("--engine fasttree requested but not found!")
        # 喂 FASTA: trimAl 的交错 phylip 后续块顶格无名字区, FastTree 会误判
        ft_in = aln_fasta if (aln_fasta and Path(aln_fasta).exists()) else aln_file
        print(f"[TREE] FastTree GTR+Gamma ({Path(ft_in).name}, SH-like supports, unrooted)")
        run_cmd([fasttree,"-nt","-gtr","-gamma","-quiet",str(ft_in)],"FastTree",outdir,
                stdout_file=tree_file)
        return tree_file
    # 大比对 (>80 序列) 跳过 ModelTest+RAxML-NG (bootstrap 极慢易超时), 直接 IQ-TREE
    try:
        with open(aln_file) as _f: n_seq_aln = int(_f.readline().split()[0])
    except Exception:
        n_seq_aln = 0
    use_raxml = (engine_preference == "raxml") or (
        engine_preference == "auto" and bool(modeltest and raxml) and n_seq_aln <= 80)
    if use_raxml:
        print("[TREE] ModelTest-NG -> RAxML-NG (注意: RAxML 产无根树)")
        try:
            mt = outdir/"modeltest"
            run_cmd([modeltest,"-i",str(aln_file),"-d","nt","-p",str(threads),"-o",str(mt)],"ModelTest-NG",outdir)
            best_model = parse_best_model(str(mt)+".out") or "GTR+G"
            print(f"  -> Model: {best_model}")
            rp = outdir/"raxml"
            run_cmd([raxml,"--all","--msa",str(aln_file),"--model",best_model,"--threads",str(threads),
                     "--bs-trees",str(bootstrap),"--prefix",str(rp),"--redo"],"RAxML-NG",outdir)
            bt = Path(str(rp)+".raxml.support")
            if not bt.exists(): bt = Path(str(rp)+".raxml.bestTree")
            if bt.exists():
                shutil.copy2(bt, tree_file); return tree_file
        except Exception as e:
            print(f"[WARNING] RAxML failed ({e}), IQ-TREE fallback.")
    if not iqtree: raise FileNotFoundError("No RAxML/IQ-TREE!")
    # 默认: UNREST 非可逆模型 → 天然有根树 (无外群定根)
    if rooting == "unrest":
        print("[TREE] IQ-TREE UNREST (non-reversible) -> rooted tree (no outgroup needed)")
        ip = outdir/"iqtree_run_unrest"
        run_cmd([iqtree,"-s",str(aln_file),"-m","UNREST","--model-joint","UNREST","-B",str(bootstrap),
                 "-T",str(threads),"--prefix",str(ip),"-redo"],"IQ-TREE UNREST",outdir)
        bt = Path(str(ip)+".treefile")
    else:
        print("[TREE] IQ-TREE ModelFinder (unrooted, MFP)")
        ip = outdir/"iqtree_run"
        run_cmd([iqtree,"-s",str(aln_file),"-m","MFP","-B",str(bootstrap),"-T",str(threads),
                 "--prefix",str(ip),"-redo"],"IQ-TREE",outdir)
        bt = Path(str(ip)+".treefile")
    shutil.copy2(bt, tree_file); return tree_file

def parse_prodigal_gff(gff_file, fasta_file):
    fasta_lens = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_file,"fasta")}
    rows=[]; pc=defaultdict(int)
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p=line.strip().split("\t")
            if len(p)>=9 and p[2]=="CDS":
                seqid,start,end = p[0].strip(),int(p[3]),int(p[4])
                if start>end: start,end=end,start
                pc[seqid]+=1
                rows.append({"nucl_id":seqid,"protein":f"{seqid}_{pc[seqid]}","start":start,"end":end,
                             "nucl_length":fasta_lens.get(seqid,0),"strand":p[6] if p[6] in['+','-'] else '+'})
    return pd.DataFrame(rows)

def prepare_synteny_and_diamond(fasta_file, outdir, threads):
    prodigal=find_executable(["prodigal"]); diamond=find_executable(["diamond"])
    if not (prodigal and diamond): raise FileNotFoundError("prodigal/diamond required")
    faa=outdir/"all_proteins.faa"; gff=outdir/"all_proteins.gff"
    dmnd=outdir/"proteins.dmnd"; dt=outdir/"diamond_all_vs_all.tsv"
    run_cmd([prodigal,"-i",str(fasta_file),"-a",str(faa),"-f","gff","-p","meta","-o",str(gff),"-q"],"Prodigal",outdir)
    pos_df=parse_prodigal_gff(gff,fasta_file)
    run_cmd([diamond,"makedb","--in",str(faa),"-d",str(dmnd),"--threads",str(threads),"--quiet"],"DIAMOND makedb",outdir)
    run_cmd([diamond,"blastp","-q",str(faa),"-d",str(dmnd),"-p",str(threads),"-o",str(dt),
             "--evalue","1e-5","--outfmt","6","qseqid","sseqid","pident","length","bitscore","--quiet"],"DIAMOND blastp",outdir)
    sim=pd.read_csv(dt,sep="\t",header=None,names=["protein1","protein2","identity","length","bitscore"])
    sim=sim[sim["protein1"]!=sim["protein2"]].copy()
    return pos_df, sim

def _root_path_len(tree, node):
    d=0.0; cur=node
    while cur is not tree.root:
        d += (cur.branch_length or 0.0)
        cur = _parent_of(tree, cur)
    return d

def _parent_of(tree, node):
    for clade in tree.get_nonterminals():
        if node in clade.clades:
            return clade
    return tree.root

def _dist_between(tree, a, b):
    # 双节点到根路径长之和近似
    return _root_path_len(tree,a) + _root_path_len(tree,b)

def _root_at_midpoint(tree, far1=None, far2=None):
    """mid 定根: 直径(最远两片叶), 选第二个为 outgroup 重根定根"""
    try:
        terms = tree.get_terminals()
        if len(terms) < 2:
            return tree
        if far1 is None:
            far1 = max(terms, key=lambda t: _root_path_len(tree, t))
        if far2 is None:
            far2 = max(terms, key=lambda t: _dist_between(tree, far1, t))
        if far1.name != far2.name:
            tree.root_with_outgroup(far2)
        return tree
    except Exception as e:
        return tree

# ============================================================================
# 分类注释 / 定根 / 配色 (移植自 vfam_trees v1.3.x, 离线 VMR taxa 表驱动)
# ============================================================================

# NCBI/ICTV rank 深度 (越大越特异); "no rank"/"clade" 视为 0
_RANK_DEPTH = {
    "realm": 1, "subrealm": 2, "kingdom": 3, "phylum": 4, "subphylum": 5,
    "class": 6, "order": 7, "family": 8, "subfamily": 9,
    "genus": 10, "subgenus": 11, "species": 12, "subspecies": 13, "strain": 14,
}

# VMR taxa 表列名 → lineage rank (root→tip 顺序)
_VMR_RANK_COLUMNS = [
    ("realm", "Realm"), ("kingdom", "Kingdom"), ("phylum", "Phylum"),
    ("class", "Class"), ("order", "Order"), ("family", "Family"),
    ("subfamily", "Subfamily"), ("genus", "Genus"), ("species", "Species"),
]

# ICTV 学名后缀 → rank (注释时兜底推断; 单词且以 virus 结尾视为属)
_ICTV_RANK_SUFFIXES = [
    ("viricotina", "subphylum"), ("viricetes", "class"), ("viricota", "phylum"),
    ("virinae", "subfamily"), ("viridae", "family"), ("virales", "order"),
    ("virae", "kingdom"), ("viria", "realm"), ("vira", "subrealm"),
]

_PLACEHOLDER_VALUES = {"", "na", "n/a", "none", "nan", "-"}
# 大树保护: taxonomy 扫描/MAD 为 O(n^2)+, 超限自动降级
_TAXO_LEAF_CAP = 1200
_MAD_LEAF_CAP = 600


def _infer_rank(name):
    lower = name.lower()
    for suffix, rank in _ICTV_RANK_SUFFIXES:
        if lower.endswith(suffix):
            return rank
    if lower.endswith("virus") and " " not in name.strip():
        return "genus"
    return ""


def build_lineage_index(taxa_meta):
    """taxa_meta (无版本号 acc → VMR 分类列) → lineage 索引; 样本 contig 无条目"""
    idx = {}
    for acc, row in (taxa_meta or {}).items():
        ranked = []
        for rank, col in _VMR_RANK_COLUMNS:
            val = str(row.get(col, "") or "").strip()
            if val.lower() in _PLACEHOLDER_VALUES:
                continue
            ranked.append({"name": val, "rank": rank})
        if ranked:
            idx[str(acc)] = {
                "lineage": [e["name"] for e in ranked],
                "lineage_ranked": ranked,
                "species": str(row.get("Species", "") or "").strip(),
                "genus": str(row.get("Genus", "") or "").strip(),
                "subfamily": str(row.get("Subfamily", "") or "").strip(),
                "family": str(row.get("Family", "") or "").strip(),
                "order": str(row.get("Order", "") or "").strip(),
                "class": str(row.get("Class", "") or "").strip(),
                "realm": str(row.get("Realm", "") or "").strip(),
            }
    return idx


def _leaf_lineage(name, lineage_index):
    """树叶名 (可能带版本号, 或是 contig) → lineage_ranked 列表; 无记录 → []"""
    base = str(name).split(".")[0].strip()
    return (lineage_index.get(base) or {}).get("lineage_ranked", [])


def normalize_support(tree):
    """支持值统一到 0-100 整数 (FastTree SH-like 为 [0,1] 需 ×100). 返回支持类型标签"""
    confs = [c.confidence for c in tree.find_clades() if c.confidence is not None]
    if not confs:
        return ""
    rescale = max(confs) <= 1.0
    for c in tree.find_clades():
        if c.confidence is not None:
            c.confidence = round(c.confidence * 100) if rescale else round(c.confidence)
    return "SH-like" if rescale else "bootstrap/UFBoot"


def root_tree(tree, lineage_index, method="taxonomy"):
    """按优先级定根: taxonomy 引导 → MAD → midpoint. 返回实际使用的方法名."""
    if method in ("none", "engine"):
        return method
    n_leaves = len(tree.get_terminals())
    if method == "taxonomy" and lineage_index:
        sid_lin = {t.name: _leaf_lineage(t.name, lineage_index)
                   for t in tree.get_terminals()}
        sid_lin = {k: v for k, v in sid_lin.items() if v}
        if sid_lin:
            if n_leaves <= _TAXO_LEAF_CAP:
                if _taxonomy_guided_root(tree, sid_lin):
                    return "taxonomy"
            else:
                print(f"[ROOT] {n_leaves} leaves > {_TAXO_LEAF_CAP}, "
                      f"taxonomy 扫描跳过 → MAD/midpoint")
    if method in ("taxonomy", "mad"):
        if n_leaves <= _MAD_LEAF_CAP:
            if _mad_root(tree):
                return "mad"
        else:
            print(f"[ROOT] MAD 跳过 ({n_leaves} leaves > {_MAD_LEAF_CAP}) → midpoint")
    try:
        tree.root_at_midpoint()
        return "midpoint"
    except Exception as e:
        print(f"[ROOT] midpoint 失败 ({e}), 保留原根")
        return "engine"


def _taxonomy_guided_root(tree, sid_to_lineage):
    """逐枝试根, 取所有内部节点平均 LCA 特异性最高者 (vfam_trees 同款)"""
    branches = list(tree.find_clades())
    if len(branches) < 3:
        return False
    best_score, best_clade, best_blen = -1.0, None, float("inf")
    for clade in branches:
        if clade is tree.root:
            continue
        score = _score_rooting(tree, clade, sid_to_lineage)
        blen = clade.branch_length or 0.0
        if score > best_score or (score == best_score and blen < best_blen):
            best_score, best_clade, best_blen = score, clade, blen
    if best_clade is None:
        return False
    try:
        # 带拆分长度重根: 在边中点干净切开, 避免 Bio.Phylo 吸收父节点产生多岔根
        tree.root_with_outgroup(best_clade,
                                outgroup_branch_length=max(best_blen / 2.0, 1e-9))
        print(f"[ROOT] taxonomy 引导定根完成 (mean specificity={best_score:.3f})")
        return True
    except Exception as e:
        print(f"[ROOT] taxonomy 定根失败 ({e})")
        return False


def _score_rooting(tree, clade, sid_to_lineage):
    """临时在 clade 边内重根, 计算平均 LCA 特异性后还原 (deepcopy 隔离)"""
    try:
        target_leaves = frozenset(c.name for c in clade.get_terminals())
        tmp_tree = copy.deepcopy(tree)
        tmp_clade = _find_clade_by_leaves(tmp_tree, target_leaves)
        if tmp_clade is None:
            return -1.0
        tmp_tree.root_with_outgroup(
            tmp_clade, outgroup_branch_length=max((tmp_clade.branch_length or 0.0) / 2.0, 1e-9))
        return _mean_lca_specificity(tmp_tree.root, sid_to_lineage)
    except Exception:
        return -1.0


def _find_clade_by_leaves(tree, target_leaves):
    for clade in tree.find_clades():
        if frozenset(c.name for c in clade.get_terminals()) == target_leaves:
            return clade
    return None


def _mean_lca_specificity(root, sid_to_lineage):
    """按终端数加权的平均 LCA 谱系深度 (谱系覆盖 <50% 的节点不计)"""
    weighted_sum, total_weight = 0.0, 0
    for clade in root.find_clades():
        if clade.is_terminal():
            continue
        terminals = clade.get_terminals()
        leaf_lins = [sid_to_lineage[lf.name] for lf in terminals
                     if lf.name in sid_to_lineage]
        if not leaf_lins or len(leaf_lins) * 2 < len(terminals):
            continue
        weighted_sum += len(_find_lca_lineage(leaf_lins)) * len(terminals)
        total_weight += len(terminals)
    return weighted_sum / total_weight if total_weight else 0.0


def _find_lca_lineage(lineages):
    """多条谱系的公共前缀 (即 LCA 谱系)"""
    if not lineages:
        return []
    common = lineages[0][:]
    for lin in lineages[1:]:
        trimmed = []
        for a, b in zip(common, lin):
            if a == b:
                trimmed.append(a)
            else:
                break
        common = trimmed
    return common


def _mad_root(tree):
    """最小祖先偏差 (MAD, Tria et al. 2017) 定根 — 纯 python 移植自 vfam_trees"""
    import collections
    terminals = tree.get_terminals()
    n = len(terminals)
    if n < 3:
        return False
    adj = {}
    for clade in tree.find_clades():
        for child in clade.clades:
            edge = child.branch_length or 0.0
            adj.setdefault(id(clade), []).append((child, edge))
            adj.setdefault(id(child), []).append((clade, edge))
    pdist = {}
    for start in terminals:
        queue = collections.deque([(start, 0.0)])
        visited = {id(start)}
        while queue:
            node, d = queue.popleft()
            if node.is_terminal() and node is not start:
                pdist[(start.name, node.name)] = d
                pdist[(node.name, start.name)] = d
            for nbr, edge in adj.get(id(node), []):
                if id(nbr) not in visited:
                    visited.add(id(nbr))
                    queue.append((nbr, d + edge))
    subdist = {}

    def _build_subdist(clade):
        if clade.is_terminal():
            subdist[id(clade)] = {clade.name: 0.0}
            return subdist[id(clade)]
        result = {}
        for child in clade.clades:
            edge = child.branch_length or 0.0
            for name, d in _build_subdist(child).items():
                result[name] = d + edge
        subdist[id(clade)] = result
        return result

    _build_subdist(tree.root)
    all_names = {t.name for t in terminals}
    n_pairs = n * (n - 1) / 2
    best_score, best_clade, best_x = float("inf"), None, 0.0
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        L = clade.branch_length
        if not L or L <= 0:
            continue
        B = subdist[id(clade)]
        A_names = all_names - B.keys()
        if not A_names:
            continue
        j0 = next(iter(B))
        a = {ni: pdist[(ni, j0)] - L - B[j0] for ni in A_names}
        A_list, B_list = list(A_names), list(B)
        sum_num = sum_den = 0.0
        for ni in A_list:
            for nj in B_list:
                d_ij = a[ni] + L + B[nj]
                if d_ij <= 0:
                    continue
                inv_d2 = 1.0 / (d_ij * d_ij)
                sum_num += (a[ni] - B[nj] - L) * inv_d2
                sum_den += inv_d2
        if sum_den == 0:
            continue
        x = max(0.0, min(L, -sum_num / (2.0 * sum_den)))
        mad = 0.0
        for ni in A_list:
            for nj in B_list:
                d_ij = a[ni] + L + B[nj]
                if d_ij <= 0:
                    continue
                rho = (a[ni] + x - B[nj] - (L - x)) / d_ij
                mad += rho * rho
        for p in range(len(A_list)):
            for q in range(p + 1, len(A_list)):
                d_ij = pdist.get((A_list[p], A_list[q]), 0.0)
                if d_ij <= 0:
                    continue
                mad += ((a[A_list[p]] - a[A_list[q]]) / d_ij) ** 2
        for p in range(len(B_list)):
            for q in range(p + 1, len(B_list)):
                d_ij = pdist.get((B_list[p], B_list[q]), 0.0)
                if d_ij <= 0:
                    continue
                mad += ((B[B_list[p]] - B[B_list[q]]) / d_ij) ** 2
        score = mad / n_pairs
        if score < best_score:
            best_score, best_clade, best_x = score, clade, x
    if best_clade is None:
        return False
    B_leaves = frozenset(c.name for c in best_clade.get_terminals())
    orig_len = best_clade.branch_length
    try:
        # 按最优分割点 x 在边内干净重根 (outgroup 侧 = orig_len - x)
        split = max(min(best_x, orig_len - 1e-9), 1e-9) if (orig_len and orig_len > 0) else 1e-9
        tree.root_with_outgroup(best_clade, outgroup_branch_length=orig_len - split)
        for child in tree.root.clades:
            child_leaves = frozenset(c.name for c in child.get_terminals())
            child.branch_length = (orig_len - split) if child_leaves == B_leaves else split
        print(f"[ROOT] MAD 定根完成 (score={best_score:.6f})")
        return True
    except Exception as e:
        print(f"[ROOT] MAD 定根失败 ({e})")
        return False


def annotate_tree(tree, lineage_index, lca_min_rank="none"):
    """内部节点 LCA 分类注释 + 冠群去重标签. 纯离线, 用 build_lineage_index 的谱系."""
    sid_ranked, sid_lin = {}, {}
    for t in tree.get_terminals():
        ranked = _leaf_lineage(t.name, lineage_index)
        if ranked:
            sid_ranked[t.name] = ranked
            sid_lin[t.name] = [e["name"] for e in ranked]
    if not sid_lin:
        print("[NOTE] 无 lineage 数据 (custom 模式或全部为 contig), 跳过内部节点注释")
        return {}
    if lca_min_rank and lca_min_rank != "none":
        min_depth = _RANK_DEPTH.get(lca_min_rank.lower(), 0)
        if min_depth:
            before = len(sid_lin)
            sid_lin = {k: v for k, v in sid_lin.items()
                       if any(_RANK_DEPTH.get(e["rank"], 0) >= min_depth
                              for e in sid_ranked[k])}
            if before - len(sid_lin):
                print(f"[NOTE] lca_min_rank={lca_min_rank}: "
                      f"排除 {before - len(sid_lin)} 条谱系过浅叶子")
    name_to_rank = {}
    for ranked in sid_ranked.values():
        for e in ranked:
            if e["rank"] != "no rank" and e["name"] not in name_to_rank:
                name_to_rank[e["name"]] = e["rank"]
    _annotate_internal_nodes(tree.root, sid_lin, name_to_rank)
    _keep_deepest_labels(tree)
    return name_to_rank


def _annotate_internal_nodes(clade, sid_to_lineage, name_to_rank):
    """递归: 内部节点名 = 其叶子的 LCA 分类单元 (谱系覆盖 <50% 不注释)"""
    if clade.is_terminal():
        return
    for child in clade.clades:
        _annotate_internal_nodes(child, sid_to_lineage, name_to_rank)
    terminals = clade.get_terminals()
    leaf_lins = [sid_to_lineage[lf.name] for lf in terminals
                 if lf.name in sid_to_lineage]
    if not leaf_lins or len(leaf_lins) * 2 < len(terminals):
        return
    lca_lineage = _find_lca_lineage(leaf_lins)
    if lca_lineage:
        taxon = lca_lineage[-1]
        clade.name = taxon
        clade._taxonomy_rank = name_to_rank.get(taxon) or _infer_rank(taxon)


def _keep_deepest_labels(tree):
    """同一分类单元只保留在最大 (冠群) 节点, 小节点上的重复标签清除"""
    seen = set()
    internal = []
    for idx, c in enumerate(tree.find_clades(order="preorder")):
        if not c.is_terminal() and c.name:
            internal.append((-len(c.get_terminals()), idx, c))
    internal.sort(key=lambda t: (t[0], t[1]))
    for _, _, clade in internal:
        taxon = clade.name or ""
        if not taxon:
            continue
        if taxon in seen:
            clade.name = ""
        else:
            seen.add(taxon)


def assign_genus_colors(tree, lineage_index):
    """属级 HLS 配色 (移植 vfam_trees colors.py):
    亚科各占一个色相带, 带内属按明度区分; 无亚科 → 属均分全色相轮.
    返回 (genus→hex, leaf_name→hex, subfamily→[genus...])"""
    leaf_taxa = {}
    for t in tree.get_terminals():
        base = str(t.name).split(".")[0].strip()
        info = lineage_index.get(base) or {}
        leaf_taxa[t.name] = {"genus": info.get("genus", ""),
                             "subfamily": info.get("subfamily", "")}
    subfamily_genera, unclassified = {}, set()
    for taxa in leaf_taxa.values():
        g = taxa["genus"]
        if not g:
            continue
        sf = taxa["subfamily"]
        if sf:
            subfamily_genera.setdefault(sf, set()).add(g)
        else:
            unclassified.add(g)
    n_groups = len(subfamily_genera) + (1 if unclassified else 0)
    genus_to_color = {}
    subfamily_to_genera = {}
    if n_groups:
        groups = sorted(subfamily_genera.keys())
        if unclassified:
            groups.append("")
        hue_step = 1.0 / n_groups
        no_subfamilies = len(subfamily_genera) <= 1
        for gi, subfamily in enumerate(groups):
            base_hue = gi * hue_step
            genera = (sorted(subfamily_genera[subfamily]) if subfamily
                      else sorted(unclassified))
            subfamily_to_genera[subfamily or "(unclassified)"] = genera
            for j, genus in enumerate(genera):
                if no_subfamilies:
                    hue = j / max(len(genera), 1)
                    r, g_, b = colorsys.hls_to_rgb(hue, 0.45, 0.80)
                else:
                    lightness = (0.35 + 0.25 * (j / max(len(genera) - 1, 1))
                                 if len(genera) > 1 else 0.45)
                    r, g_, b = colorsys.hls_to_rgb(base_hue, lightness, 0.80)
                genus_to_color[genus] = "#{:02x}{:02x}{:02x}".format(
                    int(r * 255), int(g_ * 255), int(b * 255))
    sample_gray, default_gray = "#95a5a6", "#888888"
    leaf_colors = {}
    for name, taxa in leaf_taxa.items():
        leaf_colors[name] = genus_to_color.get(taxa["genus"], default_gray)
        if not taxa["genus"]:
            leaf_colors[name] = sample_gray
    return genus_to_color, leaf_colors, subfamily_to_genera


def _safe_out_path(outdir, name):
    """规范化输出文件名并强制限定在 outdir 内 (拒绝 ../ 越界与绝对路径注入)"""
    base = Path(outdir).resolve()
    p = (base / name).resolve()
    if p.parent != base:
        raise ValueError(f"输出路径越界 (仅允许 outdir 下固定文件名): {name}")
    return p


def write_phyloxml(tree, outdir, title="", leaf_colors=None):
    """PhyloXML: 保留内部分类标签 + 支持值 + 叶色 + rooted 属性.
    (纯写出, 不解析任何外部 XML)"""
    tree.rooted = True
    if leaf_colors:
        for t in tree.get_terminals():
            try:
                t.color = leaf_colors.get(t.name, "#888888")
            except Exception:
                pass
    buf = io.StringIO()
    Phylo.write([tree], buf, "phyloxml")
    xml = buf.getvalue()
    if "rerootable" not in xml:  # Bio.Phylo 不支持该属性, 写出后注入
        xml = re.sub(r'<phylogeny(?=[\s>])', '<phylogeny rerootable="false"', xml, count=1)
    path = _safe_out_path(outdir, "phylogeny.phyloxml")
    path.write_text(xml, encoding="utf-8")
    print(f"[OUT] PhyloXML: {path}")


def write_auspice_json(tree, outdir, lineage_index, genus_colors, title=""):
    """Auspice v2 (Nextstrain) JSON — 可直接拖入 https://auspice.us 交互查看"""
    def build_node(clade, parent_div):
        div = parent_div + (clade.branch_length or 0.0)
        is_leaf = clade.is_terminal()
        name = clade.name if (is_leaf and clade.name) else (clade.name or "")
        if not is_leaf and not name and clade is tree.root:
            name = "ROOT"
        node = {"name": name,
                "node_attrs": {"div": round(div, 6)},
                "branch_attrs": {}}
        if is_leaf:
            info = lineage_index.get(str(clade.name or "").split(".")[0], {}) or {}
            for key in ("realm", "class", "order", "family", "genus", "species"):
                if info.get(key):
                    node["node_attrs"][key] = {"value": info[key]}
            node["node_attrs"]["is_sample"] = {"value": not bool(info)}
        else:
            labels = {}
            if clade.name:
                labels["taxon"] = clade.name
            if clade.confidence is not None:
                labels["support"] = str(int(round(clade.confidence)))
            if labels:
                node["branch_attrs"]["labels"] = labels
        children = [build_node(c, div) for c in clade.clades]
        if children:
            node["children"] = children
        return node

    scale = [{"value": g, "color": c} for g, c in sorted(genus_colors.items())]
    doc = {
        "version": "v2",
        "meta": {
            "title": title or "phylogeny",
            "updated": datetime.now().strftime("%Y-%m-%d"),
            "colorings": [
                {"key": "genus", "title": "Genus", "type": "categorical", "scale": scale},
                {"key": "family", "title": "Family", "type": "categorical"},
                {"key": "species", "title": "Species", "type": "categorical"},
                {"key": "is_sample", "title": "Sample contig", "type": "boolean"},
            ],
            "filters": {"dimensions": ["genus", "family", "species"]},
            "display_defaults": {"layout": "rect", "distance_measure": "div",
                                 "color_by": "genus", "branch_label": "support"},
            "panels": ["tree"],
        },
        "tree": build_node(tree.root, 0.0),
    }
    path = _safe_out_path(outdir, "tree_auspice.json")
    with path.open("w", encoding="utf-8") as h:
        json.dump(doc, h, ensure_ascii=False)
    print(f"[OUT] Auspice JSON: {path} (拖入 auspice.us 即可交互查看)")


def write_annotated_nwk(tree, outdir):
    """注释树 Newick。Newick 每节点仅一个标签: Bio.Phylo 会把节点名与支持值
    拼成一个 token (如 'Betavirus82.00'), 故写出前清掉内部节点支持值 —
    支持值完整保留在 phylogeny.phyloxml / tree_auspice.json / 原始 treefile 中"""
    out = copy.deepcopy(tree)
    for c in out.get_nonterminals():
        c.confidence = None
    path = outdir / "phylogeny.annotated.nwk"
    Phylo.write(out, str(path), "newick")
    print(f"[OUT] Annotated Newick: {path} (内部标签=LCA 分类单元; 支持值见 phyloxml/auspice)")


def write_metadata_tsv(tree, outdir, lineage_index, leaf_colors):
    """每叶子一行的元数据表: 分类谱系 + 颜色 + 是否样本 contig"""
    path = _safe_out_path(outdir, "tree_metadata.tsv")
    with path.open("w", encoding="utf-8") as h:
        h.write("leaf_id\tis_sample\tspecies\tgenus\tsubfamily\tfamily\torder\t"
                "class\trealm\tlineage\tcolor\n")
        for t in tree.get_terminals():
            info = lineage_index.get(str(t.name or "").split(".")[0], {}) or {}
            lineage = ";".join(info.get("lineage", []))
            h.write("\t".join([
                t.name or "", "yes" if not info else "no",
                info.get("species", ""), info.get("genus", ""),
                info.get("subfamily", ""), info.get("family", ""),
                info.get("order", ""), info.get("class", ""),
                info.get("realm", ""), lineage,
                leaf_colors.get(t.name, "#888888"),
            ]) + "\n")
    print(f"[OUT] Metadata TSV: {path}")


def draw_tree_synteny_figure(tree, pos_df, sim_df, taxa_meta, out_prefix,
                             genus_colors=None, draw_internal_labels=True,
                             ignore_branch_length=False, min_identity=30.0, show_support=70,
                             keep_tree_order=False):
    # tree 为已定根/注释/归一化支持值的 Bio.Phylo 树对象 (由 main 传入)
    if ignore_branch_length:
        for c in tree.find_clades(): c.branch_length=1.0
    # mid 定根: UNREST 非可逆模型已天然定根, 不再手动重根; 仅非 UNREST 引擎调用时
    # _root_at_midpoint(tree) 已移除 - UNREST 建树直接得到有根树
    try:
        if not keep_tree_order:
            tree.ladderize(reverse=True)
    except Exception:
        pass
    terminals=tree.get_terminals(); num_taxa=len(terminals)
    genome_order=[n.name for n in terminals]
    y_coords={name:i for i,name in enumerate(genome_order)}
    fig=plt.figure(figsize=(19,max(6,num_taxa*0.45)),dpi=300)
    gs=fig.add_gridspec(1,3,width_ratios=[0.09,0.40,0.51],wspace=0.05)
    ax_legend=fig.add_subplot(gs[0])
    ax_tree=fig.add_subplot(gs[1])
    # 图例列（独立轴, 仅绘图例色块, 不画树）
    def assign_pos(clade,curr_x):
        clade.x=curr_x+(clade.branch_length or 0.0)
        if clade.is_terminal():
            clade.y=y_coords[clade.name]
        else:
            for child in clade.clades: assign_pos(child,clade.x)
            clade.y=np.mean([c.y for c in clade.clades])
    assign_pos(tree.root,0.0)
    max_depth=max(n.x for n in terminals)
    # 属级背景色块
    blocks=[]; curr_g=None; start_y=None
    for n in terminals:
        base_acc=n.name.split('.')[0]
        g=str(taxa_meta.get(base_acc,{}).get('Genus','Sample_Contig'))
        if g!=curr_g:
            if curr_g is not None: blocks.append((curr_g,start_y,y_coords[n.name]-1))
            curr_g,start_y=g,y_coords[n.name]
    if curr_g is not None: blocks.append((curr_g,start_y,num_taxa-1))
    cmap_bg=plt.get_cmap('Pastel1')
    gen_colors={g:cmap_bg(i%9) for i,(g,_,_) in enumerate(blocks)}
    if genus_colors:  # 属级 HLS 配色优先 (与 PhyloXML/Auspice/图例同源)
        gen_colors={g:genus_colors.get(g,'#eeeeee') for g,_,_ in blocks}
    # 分色: 我们的样本标红色斜体, 参考按属色块
    sample_colors=gen_colors
    for g,sy,ey in blocks:
        rect=patches.Rectangle((0,sy-0.4),max_depth*1.05,(ey-sy)+0.8,
                 facecolor=gen_colors.get(g,'#eeeeee'),alpha=0.35,lw=0,zorder=0)
        ax_tree.add_patch(rect)
    def draw_clade_lines(clade):
        if not clade.is_terminal():
            ys=[c.y for c in clade.clades]
            ax_tree.plot([clade.x,clade.x],[min(ys),max(ys)],color="#2c3e50",lw=1.2)
            for c in clade.clades:
                ax_tree.plot([clade.x,c.x],[c.y,c.y],color="#2c3e50",lw=1.2)
                # 标注支持值 (bootstrap/SH-like) - 仅显示 >= show_support
                supp = getattr(c, 'confidence', None)
                if supp is not None:
                    try: supp = float(supp)
                    except Exception: supp = None
                if supp is not None:
                    if supp < 1.0: supp *= 100.0  # FastTree SH-like 0-1 → 0-100 尺度
                    if supp >= float(show_support):
                        ax_tree.text(clade.x + (c.x-clade.x)/2, c.y, str(int(round(supp))),
                                     ha='center', va='center', fontsize=6.5, color='#111111',
                                     fontweight='bold', zorder=8, bbox=dict(boxstyle='round,pad=0.12',
                                     fc='white', ec='none', alpha=0.75))
                draw_clade_lines(c)
    draw_clade_lines(tree.root)
    # 内部节点 LCA 分类标签 (斜体灰字, ≥3 叶的节点才标, 避免杂乱)
    if draw_internal_labels:
        for c in tree.get_nonterminals():
            nm=(getattr(c,'name','') or '').strip()
            if not nm:
                continue
            if len(c.get_terminals())<3:
                continue
            ax_tree.text(c.x,c.y-0.18,nm,fontsize=6.5,style='italic',
                         color='#566573',ha='center',va='bottom',zorder=9,
                         bbox=dict(boxstyle='round,pad=0.1',fc='white',ec='none',alpha=0.65))
    align_x=max_depth*1.08
    for n in terminals:
        base_acc=n.name.split('.')[0]
        sp=taxa_meta.get(base_acc,{}).get('Species','')
        label=f"{sp} ({n.name})" if sp else n.name
        ax_tree.plot([n.x,align_x],[n.y,n.y],color="#bdc3c7",linestyle=":",lw=0.8)
        is_sample = base_acc not in taxa_meta
        ax_tree.text(align_x+max_depth*0.02,n.y,label,va="center",ha="left",fontsize=8,
                     fontweight="bold" if is_sample else "normal", color="#c0392b" if is_sample else "#2c3e50")
    ax_tree.set_xlim(0,max_depth*1.6); ax_tree.set_ylim(num_taxa-0.5,-0.5)
    ax_tree.axis("off"); ax_tree.set_title("Maximum-Likelihood Phylogeny",loc="left",fontsize=11,fontweight="bold",pad=12)
    # 属图例 - 独立列 (ax_legend), 顶部色块+标签竖直排列, 完全不接触树
    y_leg = 0.92
    for g,_,_ in blocks[:12]:
        ax_legend.add_patch(patches.Rectangle((0.15,y_leg), 0.35, 0.045,
                          facecolor=gen_colors.get(g,'#ccc'), edgecolor='#555',linewidth=0.6, transform=ax_legend.transAxes))
        ax_legend.text(0.55, y_leg+0.018, g, fontsize=8, va='center', ha='left', transform=ax_legend.transAxes)
        y_leg -= 0.07
    if blocks:
        ax_legend.text(0.15, 0.97, 'Genus', fontsize=9, fontweight='bold', transform=ax_legend.transAxes, va='top')
    ax_legend.axis('off')
    ax_legend.set_ylim(-0.1, 1.2)
    # 右面板 synteny
    ax_syn=fig.add_subplot(gs[2])
    max_len=pos_df["nucl_length"].max() if not pos_df.empty else 10000
    pos_df_sub=pos_df[pos_df["nucl_id"].isin(genome_order)].copy()
    ymax=0
    for genome,y in y_coords.items():
        gl=pos_df[pos_df["nucl_id"]==genome]["nucl_length"].max()
        if pd.notna(gl):
            ymax=max(ymax,gl)
            ax_syn.hlines(y,0,gl,color="#34495e",lw=2,zorder=2)
    for row in pos_df_sub.itertuples():
        y=y_coords[row.nucl_id]
        w=max(row.end-row.start,1.0)
        rect=patches.FancyBboxPatch((row.start,y-0.15),w,0.3,boxstyle="round,pad=0.02",facecolor="#3498db",edgecolor="#2980b9",lw=0.5,zorder=4)
        ax_syn.add_patch(rect)
    prot_to_info=pos_df_sub.set_index("protein").to_dict("index")
    sim_cmap=LinearSegmentedColormap.from_list("sim",["#d5e8f7","#3498db","#2ecc71","#e67e22"])
    norm=Normalize(vmin=min_identity,vmax=100.0)
    for row in sim_df.itertuples():
        if row.identity<min_identity: continue
        p1,p2=str(row.protein1),str(row.protein2)
        if p1 in prot_to_info and p2 in prot_to_info:
            i1,i2=prot_to_info[p1],prot_to_info[p2]; g1,g2=i1["nucl_id"],i2["nucl_id"]
            if abs(y_coords[g1]-y_coords[g2])==1:
                y1,y2=y_coords[g1],y_coords[g2]
                s1=a if False else i1["start"]; e1=i1["end"]
                y1u = y1+0.15 if y1<y2 else y1-0.15
                y2u = y2-0.15 if y1<y2 else y2+0.15
                poly=patches.Polygon([(i1["start"],y1u),(i1["end"],y1u),(i2["end"],y2u),(i2["start"],y2u)],
                      closed=True,facecolor=sim_cmap(norm(row.identity)),alpha=0.35,edgecolor="none",zorder=3)
                ax_syn.add_patch(poly)
    ax_syn.set_xlim(-max_len*0.02,max_len*1.05); ax_syn.set_ylim(num_taxa-0.5,-0.5)
    ax_syn.set_xlabel("Genome Coordinates (bp)",fontsize=9,fontweight="bold")
    ax_syn.spines['top'].set_visible(False); ax_syn.spines['right'].set_visible(False); ax_syn.spines['left'].set_visible(False)
    ax_syn.set_yticks([]); ax_syn.set_title("Comparative Synteny & Homology",loc="left",fontsize=11,fontweight="bold",pad=12)
    cbar_ax=ax_syn.inset_axes([0.72,1.02,0.25,0.03])
    cbar=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=sim_cmap),cax=cbar_ax,orientation="horizontal")
    cbar.set_label("Protein Identity (%)",fontsize=7); cbar.ax.tick_params(labelsize=6)
    fig.savefig(f"{out_prefix}.png",dpi=300,bbox_inches="tight")
    fig.savefig(f"{out_prefix}.pdf",bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT] Figures: {out_prefix}.png / .pdf")

def main():
    parser=argparse.ArgumentParser(prog="acvirus_tree_pro",description="ACVirus Tree-Pro: Phylogeny & Synteny",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mode",choices=["macro","genus","lineage","custom"],default="macro")
    parser.add_argument("--target_name",type=str,help="Target taxon name (e.g. Rhabdoviridae)")
    parser.add_argument("--target_rank",choices=["Family","Genus","Species"],default="Family")
    parser.add_argument("--contigs",type=Path,help="Input query viral contigs FASTA")
    parser.add_argument("--db_taxa",type=Path,required=True)
    parser.add_argument("--db_fasta",type=Path,required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    parser.add_argument("--target_ratio",type=float,default=1.0)
    parser.add_argument("--other_count",type=int,default=3)
    parser.add_argument("--engine",choices=["auto","fasttree","raxml","iqtree"],default="auto")
    parser.add_argument("--threads",type=int,default=16)
    parser.add_argument("--bootstrap",type=int,default=1000)
    parser.add_argument("--min_identity",type=float,default=30.0)
    parser.add_argument("--ignore_branch_length",action="store_true")
    parser.add_argument("--show_support",type=float,default=70.0,help="Minimum branch support displayed (default 70)")
    parser.add_argument("--keep_tree_order",action="store_true",help="Preserve tree leaf order (no ladderizing)")
    parser.add_argument("--rooting",choices=["unrest","mfp"],default="mfp",help="Tree rooting: unrest (UNREST non-reversible, rooted) or mfp (ModelFinder, default unrooted)")
    parser.add_argument("--root_method",choices=["taxonomy","mad","midpoint","engine","none"],default=None,
                        help="建树后二次定根: taxonomy=分类学引导(默认, MAD/中点兜底); mad=最小祖先偏差; "
                             "midpoint=中点; engine=保留建树引擎的原根(unrest 默认); none=不定根")
    parser.add_argument("--no_annotate",action="store_true",help="跳过 LCA 内部节点分类注释")
    parser.add_argument("--lca_min_rank",choices=["none","family","genus","species"],default="none",
                        help="谱系深度不足该阶元的叶子不参与 LCA 投票")
    parser.add_argument("--length_ratio",type=float,default=0.3,
                        help="抽样前参考长度预筛(严格窗): 剔除偏离科内属中位数 ±ratio 的候选 (0=关闭, 默认 0.3)")
    parser.add_argument("--length_ratio_fallback",type=float,default=0.5,
                        help="长度预筛兜底窗: 属内严格窗全空时放宽 ±ratio 重试, 仍空则保底最接近中位数 1 条 (默认 0.5, 0=跳过兜底直接保底)")
    args=parser.parse_args()
    outdir=args.outdir.resolve(); outdir.mkdir(parents=True,exist_ok=True)
    df_taxa=pd.read_csv(args.db_taxa)
    if args.mode!="custom":
        if not args.target_name: parser.error("--target_name required!")
        acc_lengths = get_db_lengths(args.db_fasta, outdir) if args.length_ratio > 0 else None
        sampled_accs,taxa_meta=get_hierarchical_accessions(df_taxa,args.target_rank,args.target_name,args.mode,args.target_ratio,args.other_count,
                                                           acc_lengths=acc_lengths,length_ratio=args.length_ratio,
                                                           length_ratio_fallback=args.length_ratio_fallback)
    else:
        sampled_accs,taxa_meta=[],{}
    seqkit=find_executable(["seqkit"])
    if not seqkit: raise FileNotFoundError("seqkit required!")
    final_fasta=outdir/"analysis_sequences.fasta"
    acc_file=outdir/"sampled_accessions.txt"
    with acc_file.open("w") as f:
        for acc in sampled_accs: f.write(f"^{acc}(\\.[0-9]+)?\n")
    db_extracted=outdir/"db_extracted.fasta"
    if sampled_accs:
        run_cmd([seqkit,"grep","-r","-f",str(acc_file),str(args.db_fasta),"-o",str(db_extracted)],"Extract DB",outdir)
    # 合并 contigs + DB 序列: 保证全局唯一 ID
    # (自研 contig ID 可能与 DB 参考同名, 如 BK061149.1, FastTree 会因重名崩溃)
    _seen_ids = set()

    def _read_fasta_records(_path):
        _rid, _buf = None, []
        with _path.open() as _h:
            for _ln in _h:
                if _ln.startswith(">"):
                    if _rid is not None:
                        yield _rid, _buf
                    _rid = _ln[1:].split()[0].strip()
                    _buf = []
                else:
                    _buf.append(_ln)
            if _rid is not None:
                yield _rid, _buf

    def _uniq(_base, _mark_prefix=False):
        # 默认原名; 冲突时 (或显式要求标记) 加 sample| 前缀, 再冲突加 __2/__3
        _cand = f"sample|{_base}" if (_mark_prefix or _base in _seen_ids) else _base
        if _cand in _seen_ids:
            _i = 2
            while f"{_cand}__{_i}" in _seen_ids:
                _i += 1
            _cand = f"{_cand}__{_i}"
        _seen_ids.add(_cand)
        return _cand

    with final_fasta.open("w") as out_f:
        # DB 先写, 优先占住原始名 (参考序列保持可辨识)
        if db_extracted.exists():
            for _rid, _buf in _read_fasta_records(db_extracted):
                out_f.write(f">{_uniq(_rid)}\n")
                out_f.writelines(_buf)
        if args.contigs and args.contigs.exists():
            for _rid, _buf in _read_fasta_records(args.contigs):
                out_f.write(f">{_uniq(_rid)}\n")
                out_f.writelines(_buf)
    mafft=find_executable(["mafft"]); trimal=find_executable(["trimal"])
    aln_out=outdir/"alignment.mafft"; trimmed=outdir/"alignment.trimmed.phy"
    run_cmd([mafft,"--auto","--thread",str(args.threads),str(final_fasta)],"MAFFT",outdir,stdout_file=aln_out)
    # 统一大写 (soft-mask 小写符号会让 trimAl 报 symbol not defined)
    upper_aln=outdir/"alignment.upper.fasta"
    with aln_out.open() as _fi, upper_aln.open("w") as _fo:
        for _line in _fi: _fo.write(_line if _line.startswith(">") else _line.upper())
    aln_out=upper_aln
    run_cmd([trimal,"-in",str(aln_out),"-out",str(trimmed),"-gt","0.8","-st","0.005","-phylip"],"trimAl",outdir)
    # 保护: trim 后序列数不足则跳过 (trimal 删光全 gap 序列后不会写出文件)
    def _count_phylip(p):
        try:
            with p.open() as _f:
                return int(_f.readline().split()[0])
        except Exception:
            return 0
    n_after=_count_phylip(trimmed)
    if n_after<3:
        # 退出码 3 = 跳过 (trimal 删空后无产物; 需与真失败区分)
        print(f"[SKIP] {outdir.name}: only {n_after} sequences after trimming, skip tree")
        sys.exit(3)
    tree_file=run_phylogeny_engine(trimmed,outdir,args.threads,args.engine,args.bootstrap,rooting=getattr(args,'rooting','mfp'),aln_fasta=aln_out)

    # ── 定根 + LCA 注释 + 属级配色 + 多格式输出 (移植自 vfam_trees) ──
    tree=Phylo.read(tree_file,"newick")
    support_type=normalize_support(tree)
    if support_type:
        print(f"[TREE] 支持值归一化: {support_type} (0-100)")
    lineage_index=build_lineage_index(taxa_meta)
    root_method=args.root_method
    if root_method is None:  # UNREST 树自带根, 默认保留; 其余 taxonomy 引导
        root_method="engine" if getattr(args,'rooting','mfp')=="unrest" else "taxonomy"
    used_root=root_tree(tree,lineage_index,root_method)
    print(f"[ROOT] 定根方法: {used_root}")
    genus_colors,leaf_colors,_=assign_genus_colors(tree,lineage_index)
    if not args.no_annotate:
        annotate_tree(tree,lineage_index,args.lca_min_rank)
    title=f"{args.target_name} phylogeny" if args.target_name else outdir.name
    write_annotated_nwk(tree,outdir)
    write_phyloxml(tree,outdir,title=title,leaf_colors=leaf_colors)
    write_auspice_json(tree,outdir,lineage_index,genus_colors,title=title)
    write_metadata_tsv(tree,outdir,lineage_index,leaf_colors)

    pos_df,sim_df=prepare_synteny_and_diamond(final_fasta,outdir,args.threads)
    draw_tree_synteny_figure(tree,pos_df,sim_df,taxa_meta,out_prefix=outdir/"Tree_Synteny_Composite",
                             genus_colors=genus_colors,draw_internal_labels=not args.no_annotate,
                             ignore_branch_length=args.ignore_branch_length,min_identity=args.min_identity,
                             show_support=getattr(args,'show_support',70.0),keep_tree_order=getattr(args,'keep_tree_order',False))
    print("\n[SUCCESS] ACVirus Tree-Pro completed!")

if __name__=="__main__":
    main()
