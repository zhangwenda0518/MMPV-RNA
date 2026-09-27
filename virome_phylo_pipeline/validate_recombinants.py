#!/usr/bin/env python3
"""
validate_recombinants.py — 重组分析管线：验证层调度器 + 汇总表生成
==================================================================

功能：
  1. 读取 filter_recombinants.py 的高置信度事件 CSV
  2. 自动调度三个验证工具:
     - VirusRecom  (谱系模式或单序列模式)
     - IQ-TREE      (重组区域子比对建树)
     - Simplot++    (相似性曲线配置生成)
  3. 生成发表级汇总表

输入:
  - filter_recombinants.py 输出的 *_high_confidence.csv
  - 原始多序列比对 FASTA

输出:
  validation/
  ├── run_all.sh / run_all.ps1    一键执行脚本
  ├── master_summary.tsv          发表级汇总表
  ├── Seq_A/                      每个重组体一个子目录
  │   ├── region_100_500.fasta    重组区域 ± margin 子比对
  │   ├── iqtree/                 IQ-TREE 结果
  │   ├── virusrecom/             VirusRecom 结果
  │   └── simplot_config.txt      Simplot++ 配置
  ├── Seq_B/
  └── ...

用法:
  # 完整运行（自动执行 IQ-TREE + VirusRecom）
  python validate_recombinants.py -e high_confidence.csv -a alignment.fasta --run

  # 仅生成脚本和配置（不实际运行工具）
  python validate_recombinants.py -e high_confidence.csv -a alignment.fasta

  # 指定边距、线程数
  python validate_recombinants.py -e high_confidence.csv -a alignment.fasta --margin 300 -t 20 --run
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# ============================================================================
#  FASTA I/O
# ============================================================================

def read_fasta(fpath: str) -> dict:
    """返回 {seq_id: sequence}。"""
    seqs = {}
    cur_id, cur_seq = None, []
    with open(fpath) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            elif line:
                cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)
    return seqs


def write_fasta(seqs: dict, fpath: str):
    os.makedirs(os.path.dirname(fpath) or ".", exist_ok=True)
    with open(fpath, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")


# ============================================================================
#  序列名模糊匹配
# ============================================================================

def find_matching_ids(target: str, all_ids: Set[str]) -> List[str]:
    """在 FASTA 序列名中模糊匹配重组体/亲本名。"""
    target_clean = target.strip().strip('"\'')
    # 精确匹配
    if target_clean in all_ids:
        return [target_clean]
    # 包含关系
    matches = [sid for sid in all_ids if target_clean in sid or sid in target_clean]
    if matches:
        return matches
    # 前缀匹配（常见于从 accession 中取前几个字符的情况）
    matches = [sid for sid in all_ids if sid.startswith(target_clean.split("_")[0])]
    return matches


# ============================================================================
#  子比对提取
# ============================================================================

def extract_region_alignment(
    alignment: dict,
    start: int,
    end: int,
    margin: int = 200,
    recombinant_id: Optional[str] = None,
    parent1_id: Optional[str] = None,
    parent2_id: Optional[str] = None,
) -> Tuple[dict, int, int]:
    """从全比对中截取重组区域 ± margin，只保留相关序列。

    坐标基于比对位置（1-indexed），自动跳过大片gap。
    返回 (子比对, 实际起始位点, 实际终止位点)。
    """
    all_ids = set(alignment.keys())
    target_ids = set()
    for tid in [recombinant_id, parent1_id, parent2_id]:
        if tid:
            matches = find_matching_ids(tid, all_ids)
            target_ids.update(matches)

    # 计算比对长度（取第一条序列的长度）
    ref_len = len(next(iter(alignment.values())))
    region_start = max(1, start - margin)
    region_end = min(ref_len, end + margin)

    # 对每条序列截取区域
    sub_aln = {}
    seq_ids_to_include = target_ids if target_ids else set(alignment.keys())
    for sid in seq_ids_to_include:
        if sid in alignment:
            full_seq = alignment[sid]
            sub_aln[sid] = full_seq[region_start - 1:region_end]

    # 去除全 gap 列
    if sub_aln:
        seq_list = list(sub_aln.values())
        n_seqs = len(seq_list)
        n_cols = len(seq_list[0])
        keep_cols = []
        for j in range(n_cols):
            all_gap = all(seq_list[i][j] in "-" for i in range(n_seqs))
            if not all_gap:
                keep_cols.append(j)
        if len(keep_cols) < n_cols:
            for sid in sub_aln:
                sub_aln[sid] = "".join(sub_aln[sid][j] for j in keep_cols)

    return sub_aln, region_start, region_end


# ============================================================================
#  IQ-TREE 建树
# ============================================================================

def generate_iqtree_job(event: dict, sub_aln: dict, outdir: str,
                        margin: int, threads: int = 4) -> Optional[str]:
    """为重组区域生成并可选执行 IQ-TREE 建树。

    返回 IQ-TREE 命令字符串，如果成功执行则同时写 Newick 文件。
    """
    if len(sub_aln) < 4:
        return None  # 不足 4 条序列不建树

    iqtree_dir = os.path.join(outdir, "iqtree")
    os.makedirs(iqtree_dir, exist_ok=True)
    fasta_path = os.path.join(iqtree_dir, "region.fasta")
    write_fasta(sub_aln, fasta_path)

    cmd = (
        f"iqtree -s {fasta_path} -m MFP -bb 1000 "
        f"-nt {threads} --quiet -pre {iqtree_dir}/tree"
    )
    return cmd


# ============================================================================
#  VirusRecom 配置生成
# ============================================================================

def generate_virusrecom_job(event: dict, alignment: dict, outdir: str,
                            margin: int) -> Optional[str]:
    """为重组事件生成 VirusRecom 验证配置。

    自动聚类谱系，生成 mapping 文件和命令。
    """
    recombinant = event["recombinant"]
    all_ids = set(alignment.keys())
    rec_matches = find_matching_ids(recombinant, all_ids)
    if not rec_matches:
        # 用完整比对，所有序列作为参考
        rec_matches = [recombinant]

    vr_dir = os.path.join(outdir, "virusrecom")
    os.makedirs(vr_dir, exist_ok=True)

    # 构建子比对（重组体 + 所有其他序列，完整长度或区域）
    # VirusRecom 适合用全长比对
    sub_aln = {}
    for sid in alignment:
        sub_aln[sid] = alignment[sid]

    # 自动谱系：用成对距离聚类
    seq_ids = list(sub_aln.keys())
    n = len(seq_ids)

    # 简单聚类：p-distance < 0.05 归为一组
    distances = {}
    for i in range(n):
        for j in range(i + 1, n):
            s1, s2 = sub_aln[seq_ids[i]], sub_aln[seq_ids[j]]
            diffs = sum(1 for a, b in zip(s1, s2) if a != "-" and b != "-" and a != b)
            total = sum(1 for a, b in zip(s1, s2) if a != "-" and b != "-")
            distances[(i, j)] = diffs / max(total, 1)

    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        px, py = find(x), find(y)
        if px != py: parent[py] = px

    for (i, j), d in sorted(distances.items(), key=lambda x: x[1]):
        if d <= 0.05:
            union(i, j)

    clusters = defaultdict(list)
    for i, sid in enumerate(seq_ids):
        clusters[find(i)].append(sid)

    mapping = {}
    for cid, members in clusters.items():
        label = members[0][:20].replace(" ", "_") if len(members) == 1 else f"lineage_{cid+1}"
        for sid in members:
            mapping[sid] = label

    # 写 mapping 文件
    map_path = os.path.join(vr_dir, "lineage_map.txt")
    with open(map_path, "w") as fh:
        fh.write("#sequence name\tlineage name\n")
        for sid, lin in sorted(mapping.items()):
            fh.write(f"{sid}\t{lin}\n")

    # 写全长比对
    aln_path = os.path.join(vr_dir, "alignment.fasta")
    write_fasta(sub_aln, aln_path)

    # 确定 query 谱系名（重组体所属的谱系）
    query_lineage = mapping.get(rec_matches[0], rec_matches[0])

    n_lineages = len(set(mapping.values()))
    if n_lineages == len(seq_ids):
        # 单序列模式
        cmd = (
            f"virusrecom -a {aln_path} -q {query_lineage} "
            f"-map {map_path} -g n -m a -w 800 -s 100 "
            f"-cp 0.7 -mr 6000 -le r -o {vr_dir}/output"
        )
    else:
        cmd = (
            f"virusrecom -a {aln_path} -q {query_lineage} "
            f"-map {map_path} -g n -m p -w 100 -s 20 "
            f"-o {vr_dir}/output"
        )
    return cmd


# ============================================================================
#  Simplot++ 配置生成
# ============================================================================

def generate_simplot_config(event: dict, alignment: dict, outdir: str,
                            margin: int) -> Optional[str]:
    """生成 Simplot++ 的配置文件（query 序列 + 参考序列列表）。

    Simplot++ 命令行用法:
      simplot -i alignment.fasta -q query_id -r ref1,ref2,ref3 -o output
    """
    recombinant = event["recombinant"]
    parent1 = event.get("parent1", "")
    parent2 = event.get("parent2", "")

    all_ids = set(alignment.keys())
    rec_matches = find_matching_ids(recombinant, all_ids)

    if not rec_matches:
        return None

    query_id = rec_matches[0]
    ref_ids = []
    for p in [parent1, parent2]:
        if p:
            matches = find_matching_ids(p, all_ids)
            if matches:
                ref_ids.append(matches[0])

    # 如果没有亲本信息，用所有其他序列作为参考
    if not ref_ids:
        ref_ids = [sid for sid in all_ids if sid != query_id][:10]

    simplot_dir = os.path.join(outdir, "simplot")
    os.makedirs(simplot_dir, exist_ok=True)

    # 写子比对
    sub_ids = [query_id] + ref_ids
    sub_aln = {sid: alignment[sid] for sid in sub_ids if sid in alignment}
    aln_path = os.path.join(simplot_dir, "query_and_refs.fasta")
    write_fasta(sub_aln, aln_path)

    # 写配置（simplot++ 用命令行参数，这里生成命令）
    ref_str = ",".join(ref_ids[:20])  # 限制参考序列数
    out_prefix = os.path.join(simplot_dir, f"{query_id}_simplot")
    cmd = f"simplot -i {aln_path} -q {query_id} -r {ref_str} -o {out_prefix}"
    # 历史坑: docstring 与 master_summary 都声称有 simplot_config.txt, 但从未落盘 → 汇总表指向不存在的文件
    with open(os.path.join(simplot_dir, "simplot_config.txt"), "w", encoding="utf-8") as _sf:
        _sf.write(cmd + "\n")
    return cmd


# ============================================================================
#  汇总表生成
# ============================================================================

def generate_master_summary(
    events: List[dict],
    validation_dir: str,
    meta: dict,
) -> str:
    """生成发表级汇总表 TSV。

    列:
      recombinant, start, end, length,
      num_methods, methods, best_pvalue,
      parent1, parent2, parents_consistent,
      iqtree_result, virusrecom_result, simplot_result,
      notes
    """
    tsv_path = os.path.join(validation_dir, "master_summary.tsv")
    columns = [
        "recombinant", "breakpoint_start", "breakpoint_end", "region_length",
        "num_methods", "detection_methods", "best_pvalue",
        "major_parent", "minor_parent", "parents_consistent",
        "iqtree_tree", "virusrecom_report", "simplot_command",
        "validation_status", "notes",
    ]

    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for e in events:
            rec_name = e["recombinant"]
            val_dir = os.path.join(validation_dir, _safe_dirname(rec_name))

            methods_raw = e.get("methods_str", "") or ",".join(e.get("methods", []))
            p1 = e.get("parent1", "") or e.get("major_parent", "")
            p2 = e.get("parent2", "") or e.get("minor_parent", "")

            row = {
                "recombinant":          rec_name,
                "breakpoint_start":     int(float(e.get("start", 0))),
                "breakpoint_end":       int(float(e.get("end", 0))),
                "region_length":        int(float(e.get("length", 0))),
                "num_methods":          int(float(e.get("num_methods", 0))),
                "detection_methods":    methods_raw,
                "best_pvalue":          _fmt_pvalue(e.get("best_pvalue", "")),
                "major_parent":         p1,
                "minor_parent":         p2,
                "parents_consistent":   str(e.get("parents_consistent", "")),
                "iqtree_tree":          os.path.join(val_dir, "iqtree", "tree.treefile"),
                "virusrecom_report":    os.path.join(val_dir, "virusrecom", "output",
                                                     "Possible_recombination_event_conciseness.txt"),
                "simplot_command":      os.path.join(val_dir, "simplot_config.txt"),
                "validation_status":    "pending",
                "notes":                "",
            }
            writer.writerow(row)

    return tsv_path


def _fmt_pvalue(val) -> str:
    """格式化 p 值，处理字符串和浮点数。"""
    try:
        return f"{float(val):.2e}"
    except (ValueError, TypeError):
        return str(val)


def _safe_dirname(name: str) -> str:
    """将序列名转为安全的目录名。"""
    safe = name.strip().replace(" ", "_").replace("/", "_").replace("\\", "_")
    safe = "".join(c for c in safe if c.isalnum() or c in "._-")
    return safe[:60]


# ============================================================================
#  主流程
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="重组分析验证层调度器 + 汇总表生成",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
示例:
  # 生成脚本和配置（不执行）
  python validate_recombinants.py -e high_confidence.csv -a alignment.fasta

  # 完整运行
  python validate_recombinants.py -e high_confidence.csv -a alignment.fasta --run -t 20
        """),
    )
    parser.add_argument("-e", "--events", required=True,
                        help="filter_recombinants.py 输出的高置信度事件 CSV")
    parser.add_argument("-a", "--alignment", required=True,
                        help="原始多序列比对 FASTA")
    parser.add_argument("-o", "--outdir", default="validation",
                        help="输出目录 (默认: validation)")
    parser.add_argument("--margin", type=int, default=200,
                        help="重组区域两侧扩展 bp 数 (默认: 200)")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="线程数 (默认: 4)")
    parser.add_argument("--min-seqs-for-tree", type=int, default=4,
                        help="子比对最少序列数才建树 (默认: 4)")
    parser.add_argument("--run", action="store_true",
                        help="实际执行 IQ-TREE 和 VirusRecom（需要工具已安装）")
    parser.add_argument("--run-script-name", default="run_all",
                        help="生成的一键执行脚本名前缀")
    parser.add_argument("--fusion-model", default="",
                        help="融合分类器模型文件路径 (.pkl)")
    args = parser.parse_args()

    # --- 读取输入 ---
    if not os.path.exists(args.events):
        print(f"错误: 找不到事件文件 {args.events}")
        sys.exit(1)
    if not os.path.exists(args.alignment):
        print(f"错误: 找不到比对文件 {args.alignment}")
        sys.exit(1)

    print(f"读取事件文件: {args.events}")
    with open(args.events, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        events = list(reader)

    if not events:
        print("警告: 事件文件中没有记录")
        sys.exit(0)

    # 历史坑: docstring 声称吃 filter_recombinants 的 *_high_confidence.csv (仓库中不存在),
    # 实际调用方喂的是原始 RDP5 CSV — DictReader 裸读后 event["start"] KeyError 必崩。
    # 统一规范化为下游期望的字段 (兼容两种命名), 复用 rdp5_validate 的解析口径。
    try:
        from utils.rdp5_validate import parse_rdp5_csv
        parsed = parse_rdp5_csv(args.events)
        if parsed:
            events = parsed
            print("  [i] 已用 parse_rdp5_csv 解析原始 RDP5 CSV (无表头格式)")
    except Exception as _pe:
        print(f"  [!] parse_rdp5_csv 不可用 ({_pe}), 按 DictReader 字段继续")
    for e in events:
        e.setdefault("start", e.get("bp_start", 0))
        e.setdefault("end", e.get("bp_end", 0))
        e.setdefault("parent1", e.get("major_parent", ""))
        e.setdefault("parent2", e.get("minor_parent", ""))
        e.setdefault("num_methods", e.get("n_methods", 0))
        if not e.get("methods_str"):
            _m = e.get("methods")
            e["methods_str"] = ",".join(_m.keys()) if isinstance(_m, dict) else (",".join(_m) if isinstance(_m, list) else "")
        if not e.get("best_pvalue"):
            _m = e.get("methods")
            if isinstance(_m, dict) and _m:
                e["best_pvalue"] = min(_m.values())

    print(f"  共 {len(events)} 个高置信度重组事件")

    print(f"读取比对文件: {args.alignment}")
    alignment = read_fasta(args.alignment)
    print(f"  共 {len(alignment)} 条序列, 比对长度 {len(next(iter(alignment.values())))} bp")

    # --- 创建输出目录 ---
    os.makedirs(args.outdir, exist_ok=True)

    # --- 可选：亲本精炼 ---
    refined_path = None
    try:
        from parent_refiner import ParentRefiner
        print(f"\n{'='*60}")
        print(f"  亲本精炼 (Hamming 距离)")
        print(f"{'='*60}")
        refiner = ParentRefiner(alignment, bp_margin=args.margin)
        refined_events = refiner.refine(events)
        refined_path = os.path.join(args.outdir, "refined_events.csv")
        # 写精炼后的 CSV
        if refined_events:
            import csv as _csv
            with open(refined_path, "w", newline="", encoding="utf-8") as fh:
                w = _csv.DictWriter(fh, fieldnames=list(refined_events[0].keys()), extrasaction="ignore")
                w.writeheader()
                w.writerows(refined_events)
        n_modified = sum(1 for r in refined_events if r.get("parents_modified"))
        print(f"  {n_modified}/{len(refined_events)} 个事件的亲本被修正")

        # 导出 .rdp5ML
        rdp5ml_path = os.path.join(args.outdir, "events.rdp5ML")
        try:
            from rdp5ml_writer import events_to_rdp5ml
            events_to_rdp5ml(refined_events, rdp5ml_path)
            print(f"  .rdp5ML 已导出: {rdp5ml_path}")
            print(f"  可用 RDP5CL.exe -rdp5ml {rdp5ml_path} 重新分析")
        except ImportError:
            print(f"  (rdp5ml_writer 未找到，跳过 .rdp5ML 导出)")

        # 用精炼后的事件继续验证
        events_for_validation = refined_events
    except ImportError:
        print(f"  (parent_refiner 未找到，跳过亲本精炼)")
        events_for_validation = events

    # --- 可选：融合分类器 ---
    fusion_model = getattr(args, 'fusion_model', None) or os.environ.get("FUSION_MODEL", "")
    if fusion_model and os.path.exists(fusion_model):
        print(f"\n{'='*60}")
        print(f"  融合分类器评分")
        print(f"{'='*60}")
        try:
            from fusion_classifier import FusionClassifier
            fc = FusionClassifier.load(fusion_model)
            scored = fc.predict(events_for_validation)
            n_approved = sum(1 for s in scored if s.get("fusion_prediction") == "recombinant")
            print(f"  {n_approved}/{len(scored)} 个事件通过融合分类器 (模型={fc.model_type})")
            events_for_validation = scored
        except ImportError:
            print(f"  (fusion_classifier 未找到，跳过)")
        except Exception as e:
            print(f"  融合分类器错误: {e}")

    # --- 收集所有待执行命令 ---
    all_commands: List[Tuple[str, str]] = []  # (label, cmd)
    jobs_summary = []

    for i, event in enumerate(events_for_validation):
        rec_name = event["recombinant"]
        parent1 = event.get("parent1", "")
        parent2 = event.get("parent2", "")
        start = int(event["start"])
        end = int(event["end"])
        methods = event.get("methods_str", "")

        safe_name = _safe_dirname(rec_name)
        event_dir = os.path.join(args.outdir, safe_name)
        os.makedirs(event_dir, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"  [{i+1}/{len(events_for_validation)}] {rec_name}  ({start}–{end} bp, {event.get('num_methods','?')} 方法)")
        print(f"{'='*60}")

        # --- 1. 提取子比对 ---
        sub_aln, aln_start, aln_end = extract_region_alignment(
            alignment, start, end, args.margin, rec_name, parent1, parent2
        )
        region_path = os.path.join(event_dir, f"region_{aln_start}_{aln_end}.fasta")
        write_fasta(sub_aln, region_path)
        print(f"  子比对: {len(sub_aln)} 条序列 → {region_path}")

        # --- 2. IQ-TREE ---
        iqtree_cmd = None
        if len(sub_aln) >= args.min_seqs_for_tree:
            iqtree_cmd = generate_iqtree_job(event, sub_aln, event_dir, args.margin, args.threads)
            if iqtree_cmd:
                all_commands.append((f"{safe_name}/IQ-TREE", iqtree_cmd))
                print(f"  IQ-TREE: 已生成命令")

        # --- 3. VirusRecom ---
        vr_cmd = generate_virusrecom_job(event, alignment, event_dir, args.margin)
        if vr_cmd:
            all_commands.append((f"{safe_name}/VirusRecom", vr_cmd))
            print(f"  VirusRecom: 已生成命令")

        # --- 4. Simplot++ ---
        sp_cmd = generate_simplot_config(event, alignment, event_dir, args.margin)
        if sp_cmd:
            all_commands.append((f"{safe_name}/Simplot++", sp_cmd))
            print(f"  Simplot++: 已生成命令")

        jobs_summary.append({
            "recombinant": rec_name,
            "event_dir": event_dir,
            "sub_aln": region_path,
            "iqtree_cmd": iqtree_cmd,
            "vr_cmd": vr_cmd,
            "sp_cmd": sp_cmd,
        })

    # --- 生成一键执行脚本 ---
    _write_run_scripts(all_commands, args.outdir, args.run_script_name)

    # --- 生成汇总表 ---
    meta = {
        "alignment": os.path.abspath(args.alignment),
        "margin": args.margin,
        "threads": args.threads,
    }
    summary_path = generate_master_summary(events_for_validation, args.outdir, meta)
    print(f"\n汇总表: {summary_path}")

    # --- 可选：执行 ---
    if args.run:
        print(f"\n{'='*60}")
        print(f"  执行验证工具（共 {len(all_commands)} 个任务）")
        print(f"{'='*60}")
        for label, cmd in all_commands:
            print(f"\n--- {label} ---")
            print(f"  {cmd}")
            try:
                result = subprocess.run(cmd, shell=True, capture_output=True,
                                       text=True, encoding='utf-8', errors='replace', timeout=1800)
                if result.returncode == 0:
                    print(f"  完成")
                else:
                    print(f"  失败 (code={result.returncode}): {result.stderr[:200]}")
            except subprocess.TimeoutExpired:
                print(f"  超时（跳过）")
            except FileNotFoundError:
                print(f"  工具未安装（跳过）")

    # --- 最终输出 ---
    print(f"\n{'='*60}")
    print(f"  验证层准备完毕")
    print(f"{'='*60}")
    print(f"  输出目录:  {os.path.abspath(args.outdir)}")
    print(f"  事件数:    {len(events_for_validation)}")
    print(f"  任务数:    {len(all_commands)}")
    if not args.run:
        print(f"  → 运行 {args.run_script_name}.sh 或 {args.run_script_name}.ps1 执行所有验证任务")


def _write_run_scripts(commands: List[Tuple[str, str]], outdir: str, name: str):
    """生成 Linux 和 Windows 双版本一键执行脚本。"""

    # Linux/Mac
    sh_path = os.path.join(outdir, f"{name}.sh")
    with open(sh_path, "w") as fh:
        fh.write("#!/bin/bash\n")
        fh.write(f"# 重组分析验证层一键执行脚本\n")
        fh.write(f"# 生成时间: {__import__('datetime').datetime.now()}\n")
        fh.write(f"# 任务数: {len(commands)}\n\n")
        fh.write("set -e\n\n")
        for i, (label, cmd) in enumerate(commands):
            fh.write(f"echo '[{i+1}/{len(commands)}] {label}'\n")
            fh.write(f"{cmd}\n\n")
        fh.write(f"echo '\\n全部验证任务完成'\n")
    os.chmod(sh_path, 0o755) if hasattr(os, "chmod") else None
    print(f"  一键脚本: {sh_path}")

    # Windows PowerShell
    ps_path = os.path.join(outdir, f"{name}.ps1")
    with open(ps_path, "w", encoding="utf-8") as fh:
        fh.write(f"# 重组分析验证层一键执行脚本 (PowerShell)\n")
        fh.write(f"# 生成时间: {__import__('datetime').datetime.now()}\n\n")
        fh.write("$ErrorActionPreference = 'Continue'\n\n")
        for i, (label, cmd) in enumerate(commands):
            fh.write(f"Write-Host '[{i+1}/{len(commands)}] {label}' -ForegroundColor Cyan\n")
            fh.write(f"cmd /c \"{cmd}\"\n")
            fh.write(f"if ($LASTEXITCODE -ne 0) {{ Write-Host '  警告: 返回码=' $LASTEXITCODE -ForegroundColor Yellow }}\n\n")
        fh.write(f"Write-Host '\\n全部验证任务完成' -ForegroundColor Green\n")
    print(f"  一键脚本: {ps_path}")


if __name__ == "__main__":
    main()
