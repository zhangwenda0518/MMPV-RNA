#!/usr/bin/env python3
"""
geo_analysis.py — 地理分析扩展：智能子抽样 + RRT + TempMig
==========================================================
论文级地理分析模块:

  1. smart_subsample()      — 按时间+地点分层去冗余 (替代盲目 CD-HIT)
  2. run_rrt()              — RRT 方法一: 遗传距离矩阵法 (组间距离随机化; 无需 BEAST)
  3. run_rrt_tree()         — RRT 方法二: VirPhyKit 树注释法 (.set/.set.prob);
                              ⚠ 口径与 run_rrt **不同**, 二者结论不可互换引用
  4. run_tempmig()          — Temporal Migration Tracker (时序迁移追踪)
  + collect_randomized_mccs() / _dump_rrt_tree_summary()  — 方法二的接线辅助
  + inspect_tree_annotations() / select_root_states()     — 输入树诊断 / 根状态定位
  + run_pathways()          — 事件级传播路径表 + episode 驻留折叠 + 三级置信 + LTL
                              (utils.transmission_paths; 2026-09-16 接入, 方法出处见其 docstring)
  + run_dispersion()        — 扩散统计族 (速度/wavefront/扩散系数) + 坐标重标注置换零模型
                              (utils.diffusion_stats; MCC 描述层, 不替代 BSSVS 推断)

📌 方法二的输入树 (2026-09-16 更正): `<键>.set` / `<键>.set.prob` 由 **TreeAnnotator**
   汇总后验样本时自动生成, 不是 BEAST XML 选项、也无需自己算后验 → 直接喂
   phylogeo 的 `merged/mcc.tree` 即可 (它本就是 treeannotator 产物)。细节见
   `parse_mcc_states` docstring 与 `inspect_tree_annotations` docstring。

CLI (2026-09-16):
  --rrt-tree <真实MCC树> --rrt-randomized <随机化树目录/glob>   # 启用方法二
  --no-rrt                                                    # 只跑方法二时免去做方法一

参考文献:
  - RRT: 检验地理聚类是否显著 (VirPhyKit RRT 模块)
  - TempMig: 追踪不同时间段迁移事件 (VirPhyKit TempMig 模块)
  - 16 papers: 多数用 Mantel + BSSVS，RRT和TempMig是补充
"""

import csv, os, sys, re, random
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from Bio import SeqIO, AlignIO, Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ═══════════════════════════════════════════════════════════════════
# 1. Smart Subsampling (时空分层去冗余)
# ═══════════════════════════════════════════════════════════════════

def smart_subsample(
    fasta_file: str,
    metadata_csv: str,
    output_fasta: str,
    time_window_years: float = 1.0,
    method: str = "unique_timeloc",
    location_column: str = "location",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    智能子抽样：保留唯一(时间窗口, 地点)组合的代表序列。

    策略:
      - unique_timeloc: 同一地点+时间窗口内只保留一条（最长的）
      - stratified: 每个地点等量抽样（VirPhyKit GeoSubsampler 模式）
      - none: 仅去完全重复序列（同ID）

    与 CD-HIT 的区别:
      CD-HIT: 只看序列相似性 → 可能删掉不同地点的相同序列（丢失迁移信号）
      smart: 看(时间, 地点) → 保留时空多样性

    Returns
    -------
    dict: {n_original, n_kept, n_removed, per_location, output_fasta}
    """
    log = log or LogCollector()

    # Load sequences
    seqs = list(SeqIO.parse(fasta_file, "fasta"))

    # Load metadata
    meta = {}
    with open(metadata_csv, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('name', '').strip()
            date = row.get('date', '').strip()
            loc = row.get(location_column, '').strip()
            if name:
                meta[name] = {"date": date, "location": loc}

    # Parse decimal years
    def _to_year(d):
        try:
            return float(str(d).split('-')[0])
        except (ValueError, IndexError):
            return 0

    # Group by (time_window, location)
    groups = defaultdict(list)
    for s in seqs:
        m = meta.get(s.id, {})
        loc = m.get("location", "Unknown")
        yr = _to_year(m.get("date", "0"))
        time_bin = int(yr / time_window_years) if yr > 0 else -1
        key = (loc, time_bin)
        seq_len = len(str(s.seq))
        groups[key].append((s, seq_len, yr))

    # Keep longest per group, record what was removed
    kept = []
    removed_list = []
    removed_count = 0
    per_location = defaultdict(lambda: {"kept": 0, "removed": 0, "groups": 0})

    for (loc, time_bin), group_seqs in groups.items():
        group_seqs.sort(key=lambda x: -x[1])  # sort by length desc
        best = group_seqs[0]
        kept.append(best[0])
        per_location[loc]["kept"] += 1
        per_location[loc]["groups"] += 1
        time_label = f"{time_bin * time_window_years:.0f}-{(time_bin + 1) * time_window_years:.0f}" if time_bin >= 0 else "unknown"

        if len(group_seqs) > 1:
            removed_count += len(group_seqs) - 1
            per_location[loc]["removed"] += len(group_seqs) - 1
            for s, slen, yr in group_seqs[1:]:
                removed_list.append({
                    "removed_id": s.id,
                    "kept_id": best[0].id,
                    "location": loc,
                    "time_bin": time_label,
                    "removed_len": slen,
                    "kept_len": best[1],
                    "reason": f"duplicate (same location={loc}, time≈{time_label})",
                })

    # Write output FASTA
    SeqIO.write(kept, output_fasta, "fasta")

    # Write redundancy summary CSV
    summary_csv = os.path.splitext(output_fasta)[0] + "_redundancy_report.csv"
    with open(summary_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["removed_id", "kept_id", "location",
                                                "time_bin", "removed_len", "kept_len", "reason"])
        writer.writeheader()
        writer.writerows(removed_list)

    # Per-location summary
    loc_csv = os.path.splitext(output_fasta)[0] + "_per_location.csv"
    with open(loc_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["location", "groups", "kept", "removed", "redundancy_rate"])
        for loc in sorted(per_location.keys()):
            p = per_location[loc]
            rate = p["removed"] / (p["kept"] + p["removed"]) * 100 if (p["kept"] + p["removed"]) > 0 else 0
            writer.writerow([loc, p["groups"], p["kept"], p["removed"], f"{rate:.1f}%"])

    result = {
        "n_original": len(seqs),
        "n_kept": len(kept),
        "n_removed": removed_count,
        "per_location": dict(per_location),
        "output_fasta": output_fasta,
        "redundancy_csv": summary_csv,
        "per_location_csv": loc_csv,
    }

    log.emit(f"Smart subsample: {len(seqs)} → {len(kept)} ({removed_count} removed)")
    log.emit(f"  Groups (location × time): {len(groups)}")
    log.emit(f"  Redundancy report: {summary_csv}")
    log.emit(f"  Per-location: {loc_csv}")

    # Print per-location summary
    log.emit(f"\n  {'Location':<20s} {'Groups':>7s} {'Kept':>6s} {'Removed':>8s} {'Redundancy':>10s}")
    log.emit(f"  {'─'*20} {'─'*7} {'─'*6} {'─'*8} {'─'*10}")
    for loc in sorted(per_location.keys()):
        p = per_location[loc]
        rate = p["removed"] / (p["kept"] + p["removed"]) * 100 if (p["kept"] + p["removed"]) > 0 else 0
        log.emit(f"  {loc:<20s} {p['groups']:>7d} {p['kept']:>6d} {p['removed']:>8d} {rate:>9.1f}%")

    return result


# ═══════════════════════════════════════════════════════════════════
# 2. RRT — Region Randomization Test
# ═══════════════════════════════════════════════════════════════════

# 2026-09-15: 置换检验默认种子 (固定值, 保证同一输入两次运行得到同一 p 值)。
# 需要做 seed 敏感性分析时可显式传 run_rrt(..., seed=<其它值>)。
RRT_DEFAULT_SEED = 20260915

def run_rrt(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    n_randomizations: int = 1000,
    location_column: str = "location",
    log: Optional[LogCollector] = None,
    seed: Optional[int] = None,
) -> Dict:
    """
    Region Randomization Test: 检验遗传距离的地理聚类是否显著。

    方法:
      1. 计算真实的组间遗传距离中位数 (between-location)
      2. 将地点标签随机打乱 N 次，每次计算随机组间距离
      3. p = (随机距离 ≤ 真实距离的次数) / N

    可复现性 (2026-09-15 修复): 旧版直接用全局未播种的 `random.shuffle`,
    每次运行 p 值都不同 —— 同一份数据两次跑出两个 p, 无法写进论文。
    现改用独立的 `random.Random(seed)`, seed 默认取固定常量
    `RRT_DEFAULT_SEED`, 并把实际使用的 seed 写进返回值便于追溯。
    注意: 置换检验自身的 Monte-Carlo 误差仍在 (N 次打乱的离散性),
    所以判"显著/不显著"仍应结合效应量看, 不要只看 p 的小数位。

    与 Mantel 检验的区别:
      Mantel: 遗传距离矩阵 vs 地理距离矩阵的线性相关
      RRT: 直接比较组间遗传距离的分布 → 更直接、更robust

    Returns
    -------
    dict: {success, p_value, observed_between, randomized_mean, significant, plot}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "p_value": None, "observed_between": None,
              "randomized_between": [], "significant": None, "plot": None, "error": None}

    try:
        from Bio.Phylo.TreeConstruction import DistanceCalculator

        # Load data
        aln = AlignIO.read(fasta_file, "fasta")
        meta = {}
        with open(metadata_csv, encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                name = row.get('name', '').strip()
                loc = row.get(location_column, '').strip()
                if name and loc:
                    meta[name] = loc

        # Map sequence indices to locations
        seq_names = [r.id for r in aln]
        locations = [meta.get(n, "Unknown") for n in seq_names]
        unique_locs = sorted(set(locations))
        n_seqs = len(seq_names)

        # Compute observed between-group distances
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)

        observed = _compute_between_group_median(dm, locations, unique_locs)
        result["observed_between"] = observed

        # Randomizations
        # 2026-09-15: 独立 Random 实例 + 固定种子 -> 可复现
        used_seed = RRT_DEFAULT_SEED if seed is None else int(seed)
        rng = random.Random(used_seed)
        result["seed"] = used_seed
        result["n_randomizations"] = n_randomizations
        randomized = []
        for i in range(n_randomizations):
            shuffled = locations[:]
            rng.shuffle(shuffled)
            rnd = _compute_between_group_median(dm, shuffled, unique_locs)
            randomized.append(rnd)

        randomized = np.array(randomized)
        result["randomized_between"] = randomized.tolist()
        result["randomized_mean"] = float(randomized.mean())

        # p-value: proportion of randomized ≤ observed
        # 2026-09-15: 改用标准置换公式 (r+1)/(N+1), 与 utils/virphy_bridge.py 的
        # Mantel 检验口径一致。旧版 np.mean(randomized <= observed) 可能返回
        # 精确的 0 —— 置换检验里 p=0 不成立 (至少有一个"置换"是原排列),
        # 写进论文会被审稿人质疑, 也低估了 Monte-Carlo 误差。
        n_le = int(np.sum(randomized <= observed))
        p_value = (n_le + 1) / (n_randomizations + 1)
        result["n_perm_le_observed"] = n_le
        result["p_value"] = float(p_value)
        result["significant"] = p_value < 0.05

        # Plot
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(randomized, bins=min(50, n_randomizations//10),
                color='#b0b0b0', edgecolor='#666', alpha=0.7, label='Randomized')
        ax.axvline(observed, color='#E64B35', linewidth=3, linestyle='--',
                  label=f'Observed = {observed:.4f}')
        ax.axvline(np.percentile(randomized, 95), color='#4DBBD5', linewidth=2,
                  linestyle=':', label=f'95th percentile')
        ax.set_xlabel('Between-group median distance', fontsize=13)
        ax.set_ylabel('Frequency', fontsize=13)
        ax.set_title(f'Region Randomization Test (p = {p_value:.4f})', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10, framealpha=0.9)

        rrt_plot = os.path.join(output_dir, "rrt_results.pdf")
        fig.savefig(rrt_plot, dpi=300, format='pdf', bbox_inches='tight')
        plt.close(fig)
        result["plot"] = rrt_plot

        log.emit(f"RRT: observed={observed:.4f}, p={p_value:.4f} "
                f"({'significant' if p_value < 0.05 else 'not significant'})")
        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"RRT failed: {e}")
        return result


def _compute_between_group_median(dm, labels, unique_locs):
    """Compute median pairwise distance between different location groups."""
    n = len(labels)
    between_dists = []
    for i in range(n):
        for j in range(i + 1, n):
            if labels[i] != labels[j] and labels[i] in unique_locs and labels[j] in unique_locs:
                between_dists.append(dm[i, j])
    return float(np.median(between_dists)) if between_dists else 0.0


# ═══════════════════════════════════════════════════════════════════
# 3. TempMig — Temporal Migration Tracker
# ═══════════════════════════════════════════════════════════════════

def run_tempmig(
    mcc_tree_file: str,
    metadata_csv: str,
    output_dir: str,
    time_bins: int = 5,
    location_column: str = "location",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    Temporal Migration Tracker: 统计不同时间段的迁移事件。

    从 MCC 树的祖先节点提取地点变迁事件，按时间窗口统计。

    方法:
      1. 用 BEAST MCC 树的节点注释提取每个内节点的 (时间, 祖先地点, 后代地点)
      2. 将时间分桶，统计每个窗口的迁移事件数
      3. 绘制时序迁移柱状图

    与 BSSVS 的互补:
      BSSVS: 回答 "哪条路线有迁移" (空间)
      TempMig: 回答 "什么时候迁移多" (时间)

    Returns
    -------
    dict: {success, migration_events, time_bins, plot, summary_csv}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "migration_events": [], "error": None}

    try:
        # 格式自动检测 (NEXUS / Newick) — 历史坑: 写死 newick 会拒绝 treeannotator 的 #NEXUS 输出
        with open(mcc_tree_file, encoding='utf-8', errors='replace') as f:
            head = f.read(1024)
        fmt = 'nexus' if head.lstrip().upper().startswith('#NEXUS') else 'newick'
        tree = Phylo.read(mcc_tree_file, fmt)

        # 统一注释解析 (与 utils.tempmig_full.parse_beast_annotation 一致,
        # 支持 location={X=prob} / Location.set+Location.set.prob / Location="X" 三种格式)
        from utils.tempmig_full import parse_beast_annotation

        # Load location metadata
        tip_locations = {}
        with open(metadata_csv, encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                name = row.get('name', '').strip()
                loc = row.get(location_column, '').strip()
                if name and loc:
                    tip_locations[name] = loc

        # Also extract from tree annotations
        node_locations = {}  # id(clade) → location
        node_heights = {}    # id(clade) → height

        for clade in tree.find_clades():
            cid = id(clade)
            if clade.is_terminal():
                node_locations[cid] = tip_locations.get(clade.name, "Unknown")
                node_heights[cid] = 0.0
            else:
                info = parse_beast_annotation(clade)
                loc = info.get('location')
                if loc:
                    node_locations[cid] = max(loc, key=loc.get)
                h = info.get('height')
                if h is not None:
                    node_heights[cid] = h
                # Default from branch length (仅当无 height 注释)
                if cid not in node_heights:
                    node_heights[cid] = clade.branch_length or 0

        # Extract migration events: parent → child location change
        root_height = max(node_heights.values()) if node_heights else 100
        events = []

        for clade in tree.find_clades():
            if clade.is_terminal():
                continue
            parent_loc = node_locations.get(id(clade), "Unknown")
            for child in clade.clades:
                child_loc = node_locations.get(id(child), "Unknown")
                if parent_loc != child_loc and parent_loc != "Unknown" and child_loc != "Unknown":
                    # Time = root_height - node_height (root=oldest, tips=present)
                    time_ago = root_height - node_heights.get(id(clade), 0)
                    events.append({
                        "from": parent_loc,
                        "to": child_loc,
                        "time": round(time_ago, 2),
                    })

        result["migration_events"] = events

        if not events:
            log.warning("No migration events extracted (tree may lack location annotations)")
            result["error"] = "No migration events"
            return result

        # Time binning
        times = [e["time"] for e in events]
        min_t, max_t = min(times), max(times)
        bins = np.linspace(max_t, min_t, time_bins + 1)  # reverse: oldest→recent
        bin_labels = [f"{bins[i+1]:.0f}–{bins[i]:.0f}" for i in range(time_bins)]

        counts = np.zeros(time_bins, dtype=int)
        for e in events:
            for i in range(time_bins):
                if bins[i+1] <= e["time"] < bins[i]:
                    counts[i] += 1
                    break

        # Plot
        fig, ax = plt.subplots(figsize=(10, 5))
        colors = plt.cm.YlOrRd(np.linspace(0.3, 0.9, time_bins))
        ax.bar(range(time_bins), counts, color=colors, edgecolor='white', linewidth=1.5)
        ax.set_xticks(range(time_bins))
        ax.set_xticklabels(bin_labels, rotation=45, ha='right')
        ax.set_xlabel('Time (years before present)', fontsize=13)
        ax.set_ylabel('Number of migration events', fontsize=13)
        ax.set_title('Temporal Migration Tracker', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plot_path = os.path.join(output_dir, "tempmig.pdf")
        fig.savefig(plot_path, dpi=300, format='pdf', bbox_inches='tight')
        plt.close(fig)
        result["plot"] = plot_path

        # Summary CSV
        summary = [{"time_bin": bl, "n_events": int(c)} for bl, c in zip(bin_labels, counts)]
        csv_path = os.path.join(output_dir, "tempmig_summary.csv")
        pd.DataFrame(summary).to_csv(csv_path, index=False)
        result["summary_csv"] = csv_path

        log.emit(f"TempMig: {len(events)} migration events across {time_bins} time bins")
        for s in summary:
            log.emit(f"  {s['time_bin']}: {s['n_events']} events")

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"TempMig failed: {e}")
        return result


def _parse_location_set(v: str) -> Dict[str, float]:
    """Parse BEAST location annotation like {Ningxia=0.8, Beijing=0.2}"""
    v = v.strip('{}')
    result = {}
    for item in v.split(','):
        item = item.strip()
        if '=' in item:
            loc, prob = item.rsplit('=', 1)
            try:
                result[loc.strip()] = float(prob.strip())
            except ValueError:
                pass
    return result


def parse_mcc_states(tree_content: str, prefix: Optional[str] = None
                     ) -> Tuple[List[str], List[float]]:
    """解析 BEAST MCC 树的离散性状注释 (VirPhyKit RRT 格式)

    **关键认知 (2026-09-16 更正)** —— 这两个注释键既不是 BEAST XML 的选项、也不是
    下游自己算的后验，而是 **TreeAnnotator 自动生成**的：

      TreeAnnotator 汇总后验树样本时，若某节点属性在样本里是**字符串(离散)**，
      就对该节点做频次统计并写出三个键 (BEAST 1.x `dr.app.tools.TreeAnnotator$CladeSystem`
      常量池含 `.set`/`.set.prob`/`.prob`；BEAST 2 同源 `TreeAnnotator.java`:
      `isDiscrete = v[i] instanceof String` → `annotateModeAttribute` + `annotateFrequencyAttribute`
      → `node.setMetaData(label + ".set", ...)` / `(label + ".set.prob", ...)`)：

        <属性>=<众数>                   ← 点状态, tie 时会拼成 "A+B"
        <属性>.prob=<众数占比>           ← 众数概率 (tie 时含计数修正)
        <属性>.set={A,B,...}            ← 样本中出现过的状态
        <属性>.set.prob={0.6,0.3,...}   ← 各自的后验频率 (归一化到 1)

      实证: 用本地 BEAST1 jar 跑
        java -cp beast1.jar dr.app.tools.TreeAnnotator -burnin 0 -heights median in.trees out.tree
      输入节点属性名 `location.states`(字符串) 时，输出 MCC 里确实出现
        location.states.set={"Beijing","Ningxia"}, location.states.set.prob={0.5,0.5}

    因此 `parse_mcc_states` 的正则**与键名无关**：`Region.set` / `location.states.set` /
    `max.set` 都能吃 —— 喂它 treeannotator 产出的 MCC 即可，无需改 BEAST XML。

    取第几个？**根节点**。Newick 里根的 `[&...]` 注释块排在整棵树最后，所以按
    注释块从后往前找「同时含 `<前缀>.set` 与 `<前缀>.set.prob` 且前缀一致」的第一块。
    (2026-09-16 之前是全文正则取最后一个匹配 —— 对单性状树恰好等价，但若树里同时有
     两个离散性状 / 值里含 `.set` 字样就会串味，故改为按块解析。)

    Parameters
    ----------
    tree_content : str
        树文件全文 (MCC 树; 逐样本 `.trees` 只有点状态, 解析不出)
    prefix : str, optional
        显式指定性状前缀 (如 `location.states`)。给定时优先用它。

    Returns
    -------
    (states, probs): 状态名列表 + 概率列表; 解析不到返回 ([], [])
    """
    states, probs, _ = select_root_states(tree_content, prefer_prefix=prefix)
    return states, probs


# ── RRT(树注释法) 接线辅助 (2026-09-16 新增) ─────────────────────────────────
_ANNOT_BLOCK_RE = re.compile(r'\[&([^\]]*)\]')
_ANNOT_KEY_RE = re.compile(r'([A-Za-z_][\w.]*)\s*=')
# 前缀可选: 兼容 `set={...}`(裸) 与 `Region.set={...}`
_ANNOT_SET_RE = re.compile(r'(?:([A-Za-z_][\w.]*)\.)?set=\{\s*([^}]*?)\s*\}')
_ANNOT_SETPROB_RE = re.compile(r'(?:([A-Za-z_][\w.]*)\.)?set\.prob=\{\s*([^}]*?)\s*\}')
_ANNOT_QUOTED_RE = re.compile(r'(?:^|,)\s*([A-Za-z_][\w.]*)\s*=\s*"')
_TREE_LINE_RE = re.compile(r'(?m)^\s*(?:tree|utree)\s')


def _block_state_pairs(block: str) -> Dict[str, Tuple[List[str], List[float]]]:
    """一个 `[&...]` 注释块里所有「前缀一致且长度自洽」的 set/set.prob 对。"""
    sets = {m.group(1) or '': m.group(2) for m in _ANNOT_SET_RE.finditer(block)}
    probs = {m.group(1) or '': m.group(2) for m in _ANNOT_SETPROB_RE.finditer(block)}
    out: Dict[str, Tuple[List[str], List[float]]] = {}
    for k, raw_states in sets.items():
        if k not in probs:
            continue
        states = [s.strip().strip('"').strip("'")
                  for s in raw_states.split(',') if s.strip()]
        try:
            pv = [float(p.strip()) for p in probs[k].split(',') if p.strip()]
        except ValueError:
            continue
        if not states or len(states) != len(pv):
            continue
        out[k] = (states, pv)
    return out


def select_root_states(tree_content: str, prefer_prefix: Optional[str] = None
                       ) -> Tuple[List[str], List[float], Dict]:
    """定位**根节点**的状态后验向量, 并给出可追溯的诊断信息。

    Returns
    -------
    (states, probs, info)
        info: {n_blocks, n_trees, block_index, prefix, candidates, n_states, tail_block}
              `prefix=None` 表示没解析到。
    """
    blocks = _ANNOT_BLOCK_RE.findall(tree_content)
    info: Dict = {"n_blocks": len(blocks),
                  "n_trees": len(_TREE_LINE_RE.findall(tree_content)),
                  "block_index": None, "prefix": None,
                  "candidates": [], "n_states": 0, "tail_block": False}

    # 显式指定前缀 → 从后往前找第一个含该前缀的注释块 (根的注释块在最末)
    if prefer_prefix is not None:
        for idx in range(len(blocks) - 1, -1, -1):
            pairs = _block_state_pairs(blocks[idx])
            if prefer_prefix in pairs:
                states, probs = pairs[prefer_prefix]
                info.update({"block_index": idx, "prefix": prefer_prefix,
                             "candidates": sorted(pairs), "n_states": len(states),
                             "tail_block": idx == len(blocks) - 1})
                return states, probs, info

    # 默认: 最后一个含 set/set.prob 对的注释块 = 根节点
    for idx in range(len(blocks) - 1, -1, -1):
        pairs = _block_state_pairs(blocks[idx])
        if not pairs:
            continue
        # 同一块里若有多个性状前缀, 取状态数最多的 (确定性 tie-break: 前缀字典序)
        key = sorted(pairs, key=lambda k: (-len(pairs[k][0]), k))[0]
        states, probs = pairs[key]
        info.update({"block_index": idx, "prefix": key,
                     "candidates": sorted(pairs), "n_states": len(states),
                     "tail_block": idx == len(blocks) - 1})
        return states, probs, info
    return [], [], info


def inspect_tree_annotations(tree_path: str) -> Dict:
    """诊断一棵树里**实际存在**的节点注释标签名, 并在缺 `*.set` 时给出**正确**处方。

    为什么需要它 (2026-09-16 更正)：早期版本这里写着"要把 BEAST XML 的 tree logger 改成
    输出 `列名.set`/`列名.set.prob`，再重跑 BEAST" —— **这条指引是错的**。经查证：
    `.set`/`.set.prob` 由 **TreeAnnotator** 在汇总后验样本时自动生成（见
    `parse_mcc_states` docstring 的机制说明与实证），BEAST 的 logTree 只负责把
    **逐样本点状态**写进 `.trees`。所以"树里只有点状态"的唯一含义是：
    **你手上的是逐样本 `.trees`（或 treeannotator 之前的目标树），不是 MCC。**

    三类常见输入与处方：
      · MCC 树 (treeannotator 产出)  → 有 `*.set`/`*.set.prob`, ok=True
      · 逐样本 `.trees`              → 只有 `X="状态"` 点状态, ok=False → 跑 treeannotator
      · 无注释树 (IQ-TREE 树等)       → 无 `[&...]` 块, ok=False → 选错文件

    Returns
    -------
    dict: {n_annot_blocks, n_trees, tags, set_tags, point_state_tags,
           has_set, has_set_prob, ok, hint}
    """
    try:
        with open(tree_path, encoding='utf-8', errors='replace') as f:
            txt = f.read()
    except Exception as e:                                    # noqa: BLE001
        return {"n_annot_blocks": 0, "n_trees": 0, "tags": [], "set_tags": [],
                "point_state_tags": [], "has_set": False, "has_set_prob": False,
                "ok": False, "hint": f"读取失败: {e}"}
    blocks = _ANNOT_BLOCK_RE.findall(txt)
    tags = set()
    point_state_tags = set()
    for b in blocks:
        tags.update(_ANNOT_KEY_RE.findall(b))
        point_state_tags.update(_ANNOT_QUOTED_RE.findall(b))
    set_tags = sorted(t for t in tags if t == 'set' or t.endswith('.set')
                      or t == 'set.prob' or t.endswith('.set.prob'))
    has_set = any(t == 'set' or t.endswith('.set') for t in tags)
    has_prob = any(t == 'set.prob' or t.endswith('.set.prob') for t in tags)
    ok = has_set and has_prob
    n_trees = len(_TREE_LINE_RE.findall(txt))
    hint = ""
    if not ok:
        if not blocks:
            hint = ("树里没有任何 [&...] 节点注释 → 多半选错文件了（要 MCC 树，"
                    "不是未注释的 IQ-TREE 树）")
        elif has_set or has_prob:
            hint = ("只有 `*.set` 或只有 `*.set.prob`，两者不成对 → 疑似 treeannotator "
                    "输出被截断/被文本工具改写过，重新生成一次 MCC")
        elif point_state_tags:
            hint = (
                "树里是**逐样本点状态**（%s=\"X\"）而非 MCC → 这是 BEAST 的 `.trees`"
                "（或 treeannotator 之前的输入树）。处方：对它跑 TreeAnnotator 得到 MCC，"
                "`.set`/`.set.prob` 会自动出现（无需改 BEAST XML）：\n"
                "  treeannotator -burnin <N> -heights median <你的>.trees mcc.tree"
                % "/".join(sorted(point_state_tags)[:4]))
            if n_trees > 1:
                hint += "\n  (该文件含 %d 棵取样树 → 确认是后验 `.trees`)" % n_trees
        elif tags:
            hint = ("树有注释但不是离散性状（实际: %s）→ 该分析没记录地点性状，"
                    "需确认 phylogeo BEAST 的 logTree 里挂了 "
                    "<trait name=\"location.states\" tag=\"Location\">"
                    % ", ".join(sorted(tags)[:8]))
        else:
            hint = "有 [&...] 块但未解析出键名 → 树格式异常"
    return {"n_annot_blocks": len(blocks), "n_trees": n_trees, "tags": sorted(tags),
            "set_tags": set_tags, "point_state_tags": sorted(point_state_tags),
            "has_set": has_set, "has_set_prob": has_prob, "ok": ok, "hint": hint}



def _natural_key(path: str):
    """自然排序键: `..._random2.tre` 排在 `..._random10.tre` 之前。

    上游示例命名是 `PVS_random1..20.tre`, 纯字典序会得到 1,10,11,…,2 的顺序。
    由于 `run_rrt_tree` 用列表序号写 `Random{i}` 行标签, 顺序错了会让
    "第 i 次随机化"与文件名错位 —— 不影响判定, 但影响表格可读性与追溯。
    """
    parts = re.split(r'(\d+)', os.path.basename(path))
    return [(1, int(t)) if t.isdigit() else (0, t.lower()) for t in parts if t != '']


def collect_randomized_mccs(spec: str, exclude: Optional[str] = None) -> List[str]:
    """收集区域随机化 MCC 树 (RRT 树注释法的输入)。

    spec 支持三种写法:
      · 目录 → 优先 `*random*` 命名的树 (上游 `PVS_random1..20.tre` 风格);
               该目录若没有这种命名, 退化为目录下全部 `*.tre/*.tree/*.trees/*.nwk`
      · 含 `*` `?` `[` 的 glob → 直接展开
      · 单个文件路径 → 就那一棵

    exclude: 真实 MCC 树路径; 若它恰好落在同一目录里会被剔除 (否则会被当成"随机化树",
             把真实后验混进随机化分布, 直接拉低通过率)。

    返回自然序 (见 `_natural_key`) 的绝对/相对原样路径列表。
    """
    import glob as _glob
    spec = os.path.expanduser(str(spec or ''))
    files: List[str] = []
    if os.path.isdir(spec):
        for pat in ("*random*.tre*", "*random*.tree*", "*random*.nwk"):
            hit = _glob.glob(os.path.join(spec, pat))
            if hit:
                files = hit
                break
        if not files:
            for pat in ("*.tre", "*.tree", "*.trees", "*.nwk"):
                files += _glob.glob(os.path.join(spec, pat))
            files = list(set(files))
        files = sorted(files, key=_natural_key)
    elif any(ch in spec for ch in "*?["):
        files = sorted(_glob.glob(spec), key=_natural_key)
    elif os.path.exists(spec):
        files = [spec]
    if exclude:
        ex = os.path.abspath(exclude)
        files = [f for f in files if os.path.abspath(f) != ex]
    return files


def _dump_rrt_tree_summary(res: Dict, outdir: str, log: Optional[LogCollector] = None) -> str:
    """把 `run_rrt_tree` 的结论落成 `rrt_tree_summary.json` (供 pipeline 汇总/论文取数)。

    只挑结论性字段, 不写 `randomized_probs` 的完整数组以外的东西 (那也已经很短)。
    """
    import json
    log = log or LogCollector()
    keys = ("success", "passed", "skipped", "error",
            "real_max_prob", "target_state", "max_random_prob",
            "margin", "second_best", "tie_states", "stable",
            "n_randomized", "n_skipped", "n_missing_target",
            "annotation_prefix", "randomized_probs", "table", "plot")
    payload = {k: res.get(k) for k in keys}
    p = os.path.join(outdir, "rrt_tree_summary.json")
    os.makedirs(outdir, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    if res.get("success"):
        log.emit(f"树注释法 RRT 结论: {'PASSED' if res.get('passed') else 'FAILED'}"
                 f" (target={res.get('target_state')}, margin={res.get('margin')}) → {p}")
    return p


def run_rrt_tree(
    original_mcc: str,
    randomized_mccs: List[str],
    output_dir: str,
    log: Optional[LogCollector] = None,
) -> Dict:
    """VirPhyKit 风格 RRT (区域随机化检验, 树注释法)

    复刻 VirPhyKit RRT 模块的核心逻辑:
      1. 解析真实 MCC 树的 .set/.set.prob 注释 (根状态后验概率)
      2. 找最大后验概率的 target 状态及其概率
      3. 解析 N 棵"区域随机化"的 MCC 树 (地理状态随机打乱后重跑 BEAST)
      4. 判定: `真实 max_prob > max(随机树中**同一个 target 状态**的后验)`
         → PASSED (地理信号显著)

    ⚠️ 第 4 条的口径务必看清（2026-09-16 更正，原来这里写成"所有随机化 max_prob"是**错的**）：
       · **本函数（口径 1，正确）**：真实树取 argmax 状态 → 在每棵随机树里**只读该同名状态**
         的概率，与真实树该状态的概率比。等价于 VirPhyKit
         `function_rrt.py:91` 的 `[data.get(target_country, 0.0) for data in all_data]`。
       · **口径 2（错误，别用）**：各棵随机树各自取自己的 `max` 再与真实 max 比 ——
         随机树里任一状态偶然集中就会被当成"更显著"。姊妹平台实测该口径会误判
         （原始 0.3958 落在随机树各自 max 的区间内）。
       代码 `:551-556` 实现的是**口径 1**（`rp[rs.index(target_state)]`），已实测确认。

    ✅ 接线状态（2026-09-16 **已接入**，原为"正确的死代码"）：
       · `geo_analysis.main()` 新增 `--rrt-tree` / `--rrt-randomized`（+ `--no-rrt` 跳过距离矩阵法），
         随机化树由 `collect_randomized_mccs()` 按自然序收集（上游示例 `PVS_random1..20.tre` 风格）。
       · 结论落 `rrt_tree_summary.json`，PDF/CSV 落 `<outdir>/rrt_tree/`。
       · `phylo_pipeline.py::run_stage_geo` 按 `--rrt_randomized_dir`（或 config `phylogeo.rrt_randomized_dir`）
         显式转发；**默认关闭**——因为它需要 N 棵随机化 MCC 树（上游口径 N=20，即 ~21 次
         完整 BEAST 重跑），且会改变地理显著性结论，不能默认替用户开启。
       · 随机化树**必须自己准备**（VirPhyKit 原版同样是外部输入，`main_rrt.py` 只做 GUI 选择）。
         生成方法：把地点标签在序列间打乱 → 重跑 phylogeo BEAST → treeannotator 取 MCC，
         重复 N 次。本函数只负责"拿到这些树之后的判定"。
       · 平局/多状态/缺 target 的稳健性诊断已补（见下方 3b），但**判定语义仍严格照上游**
         （argmax 取第一个）—— 诊断只 warn，不改变结论，这样才能声称"与 VirPhyKit 同口径"。
       · **注释键不是阻塞**（2026-09-16 更正，早期文档曾误记为"要改 BEAST XML 的 tree
         logger"）：`.set`/`.set.prob` 是 **TreeAnnotator 汇总后验样本时自动生成**的
         （机制+实证见 `parse_mcc_states` docstring），所以本管线 phylogeo 产出的
         `merged/mcc.tree` 本来就带 `location.states.set` / `location.states.set.prob`。
         正则与键名无关 —— `Region.set` / `max.set` / `location.states.set` 都能吃。
         若树里只有点状态，说明拿到的是逐样本 `.trees`，跑一次 treeannotator 即可
         （`inspect_tree_annotations` 会给出确切命令）。

    与 run_rrt (距离矩阵法) 的区别:
      run_rrt       基于遗传距离矩阵的组间距离随机化 (自写, 不需 BEAST)
      run_rrt_tree  基于 BEAST MCC 树离散性状后验概率 (VirPhyKit 原版)

    Returns
    -------
    dict: {success, passed, real_max_prob, target_state, randomized_probs, table, plot}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "passed": None, "real_max_prob": None,
              "target_state": None, "randomized_probs": [],
              "table": None, "plot": None, "error": None}

    # 0. 空随机化集守卫 (2026-09-16 接线时补)
    #    若 randomized_mccs 为空, 下游 max(random_target_probs) 落到 `0.0` 兜底,
    #    于是 real_max_prob > 0 恒成立 → 静默 PASS。这是**假阳性**, 必须拦住。
    if not randomized_mccs:
        result["error"] = ("未提供随机化 MCC 树 → 无法判定 (空列表会让 max_random=0 "
                           "而假性 PASS)。上游 RRT 需 1 棵真实 + N 棵区域随机化树 "
                           "(VirPhyKit Example 用 N=20)")
        log.error(result["error"])
        return result

    try:
        # 1. 真实 MCC 树
        with open(original_mcc) as f:
            states, probs, ann = select_root_states(f.read())
        if not states:
            # 只报"无注释"会让人无从下手 —— 附上实际标签诊断 + **正确**处方
            # (2026-09-16 更正: 早期这里的建议是"改 BEAST XML", 错的; 见
            #  inspect_tree_annotations docstring)。
            diag = inspect_tree_annotations(original_mcc)
            result["annotation_diag"] = diag
            result["error"] = ("真实 MCC 树无 `<键>.set` / `<键>.set.prob` 注释 → "
                               "解析不到根状态后验。诊断: "
                               + (diag.get("hint") or "未知"))
            log.error(result["error"])
            if diag.get("tags"):
                log.error(f"  实际标签: {', '.join(diag['tags'][:12])}"
                          f" (共 {diag['n_annot_blocks']} 个注释块"
                          f"{', %d 棵树' % diag['n_trees'] if diag.get('n_trees') else ''})")
            return result

        real_prefix = ann.get("prefix")
        result["annotation_prefix"] = real_prefix
        real_max_prob = max(probs)
        target_idx = probs.index(real_max_prob)
        target_state = states[target_idx]

        # 2. 随机化 MCC 树
        random_target_probs = []
        random_rows = []
        n_skipped = 0            # 树里没解析出 .set/.set.prob 注释 → 跳过
        n_missing_target = 0     # 树里没有 target 状态 → 按上游口径记 0.0
        for i, mcc in enumerate(randomized_mccs, 1):
            with open(mcc) as f:
                # 与真实树读**同一性状前缀** (树里若含多个离散性状, 防止串味)
                rs, rp = parse_mcc_states(f.read(), prefix=real_prefix)

            if not rs:
                n_skipped += 1
                log.warning(f"  随机化树无 .set/.set.prob 注释, 跳过: {os.path.basename(mcc)}")
                continue
            row = {"Replicates": f"Random{i}"}
            for s, p in zip(rs, rp):
                row[s] = p
            random_rows.append(row)
            # target 状态在随机树中的概率（上游 :91 `data.get(target_country, 0.0)` 口径）
            if target_state in rs:
                random_target_probs.append(rp[rs.index(target_state)])
            else:
                random_target_probs.append(0.0)
                n_missing_target += 1

        # 3. 判定 (VirPhyKit 原逻辑)
        max_random = max(random_target_probs) if random_target_probs else 0.0
        passed = real_max_prob > max_random

        # 3b. 平局/稳健性诊断 (2026-09-16 接线时补)
        # 口径仍与上游一致: argmax 取**第一个**（上游 `probs.index(max_prob)`）。
        # 但平局会让 target_state 依赖状态在注释里的顺序 —— 必须报出来, 不能沉默。
        tie_states = [s for s, p in zip(states, probs)
                      if abs(p - real_max_prob) < 1e-12]
        _sorted_p = sorted(probs, reverse=True)
        second_best = _sorted_p[1] if len(_sorted_p) > 1 else 0.0
        margin = real_max_prob - max_random
        if len(tie_states) > 1:
            log.warning(f"  ⚠ target 状态平局: {tie_states} 同为 {real_max_prob:.4g}; "
                        f"按上游口径取第一个 {target_state} —— 结论可能随状态顺序变化")
        if second_best and (real_max_prob - second_best) < 0.02:
            log.warning(f"  ⚠ target 与次高状态仅差 {real_max_prob - second_best:.4g} "
                        f"({target_state} vs 次高 {second_best:.4g}) → argmax 不稳固")
        if n_missing_target:
            log.warning(f"  ⚠ {n_missing_target}/{len(random_rows)} 棵随机化树里没有 target "
                        f"状态 {target_state}, 按上游口径记 0.0")
        if margin <= 0.02:
            log.warning(f"  ⚠ 与随机化最大值的差距仅 {margin:.4g} → 判定处在临界")
        # 根节点的注释块理应排在整棵树最末 (Newick 结构决定)。若不是, 说明该文件里
        # 根注释不在尾部 (多树文件/异常写法) —— 我们取到的是"最后一个带 set 的块",
        # 不一定是根, 必须报出来而不是沉默。
        if not ann.get("tail_block"):
            log.warning(f"  ⚠ 根状态取自第 {ann.get('block_index')} 个注释块 (共 "
                        f"{ann.get('n_blocks')} 块), 不在文件末尾 → 可能不是根节点, "
                        f"请确认输入是单棵 MCC 树")
        if ann.get("candidates") and len(ann["candidates"]) > 1:
            log.warning(f"  ⚠ 该注释块含多个离散性状前缀 {ann['candidates']}, "
                        f"已取状态数最多者 `{real_prefix}` —— 若不对请用 prefix= 指定")

        # 4. 表 (VirPhyKit generate_table 风格: Random → 空行 → Real → Min → Max)
        df = pd.DataFrame(random_rows)
        real_row = {"Replicates": "Real"}
        for s, p in zip(states, probs):
            real_row[s] = p
        # 空行分隔: 上游 :80-81 用 `[""] + [""]*(n-1)` —— 是**空字符串** `""` 而不是 NaN,
        # 落盘成 `""`、读回为 ''。2026-09-16 对拍时发现原先用 `pd.DataFrame([{}])`
        # 写的是真空(读回 NaN), 语义等价但字节不同; 照抄上游后 CSV 可与上游逐字节比对。
        _empty_row = pd.DataFrame([[""] + [""] * (len(df.columns) - 1)], columns=df.columns)
        df = pd.concat([df, _empty_row,
                        pd.DataFrame([real_row])], ignore_index=True)
        if len(random_rows) > 0:
            min_row = {"Replicates": "Min"}
            max_row = {"Replicates": "Max"}
            num_cols = df.columns[1:]
            for c in num_cols:
                min_row[c] = df[c].iloc[:-2].min()               # 只取随机行
                max_row[c] = df[c].iloc[:-2].max()
            df = pd.concat([df, pd.DataFrame([min_row]), pd.DataFrame([max_row])],
                           ignore_index=True)
        table_path = os.path.join(output_dir, "rrt_tree_table.csv")
        df.to_csv(table_path, index=False, encoding='utf-8-sig')

        # 5. 图 (VirPhyKit plot_graph_from_csv 风格)
        if len(random_rows) > 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            cols = df.columns[1:]
            real_vals = df[df["Replicates"] == "Real"].iloc[0, 1:]
            min_vals = df[df["Replicates"] == "Min"].iloc[0, 1:]
            max_vals = df[df["Replicates"] == "Max"].iloc[0, 1:]
            ax.plot(cols, real_vals, label="Real", marker='o', color='#5BD4D8')
            ax.plot(cols, min_vals, label="Min", marker='o',
                    color='#FFDBB3', linestyle='--')
            ax.plot(cols, max_vals, label="Max", marker='o',
                    color='#FE6993', linestyle='--')
            ax.set_title("Region Randomization Test (VirPhyKit)")
            ax.set_xlabel("Region")
            ax.set_ylabel("Probability")
            ax.legend()
            ax.grid(True)
            plot_path = os.path.join(output_dir, "rrt_tree.pdf")
            fig.savefig(plot_path, dpi=300, format='pdf', bbox_inches='tight')
            plt.close(fig)
            result["plot"] = plot_path

        result.update({
            "success": True, "passed": passed,
            "real_max_prob": real_max_prob, "target_state": target_state,
            "randomized_probs": random_target_probs,
            "max_random_prob": max_random,
            "table": table_path,
            # ── 稳健性诊断 (2026-09-16 接线时补; 不改变判定语义, 只暴露风险) ──
            "n_randomized": len(random_rows),
            "n_skipped": n_skipped,
            "n_missing_target": n_missing_target,
            "tie_states": tie_states,
            "second_best": second_best,
            "margin": margin,
            "stable": bool(len(tie_states) == 1 and margin > 0.02),
        })
        log.emit(f"VirPhyKit RRT: 真实 max_prob={real_max_prob:.4f} "
                 f"({target_state}) vs 随机化 max={max_random:.4f} → "
                 f"{'PASSED' if passed else 'FAILED'} "
                 f"[有效随机化 {len(random_rows)}/{len(randomized_mccs)} 棵"
                 f"{', 跳过 %d' % n_skipped if n_skipped else ''}]")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"run_rrt_tree failed: {e}")
        return result


def main():
    """CLI 入口: 智能子抽样 + RRT 距离矩阵法 + (可选) RRT 树注释法 + (可选) TempMig"""
    import argparse
    ap = argparse.ArgumentParser(description="地理分析扩展 (子抽样 + RRT + TempMig)")
    ap.add_argument("--fasta", required=True, help="比对 fasta")
    ap.add_argument("--metadata", required=True, help="name,location CSV")
    ap.add_argument("--mcc-tree", default=None, help="MCC 树 (TempMig 可选)")
    ap.add_argument("-o", "--outdir", default="geo_analysis", help="输出目录")
    ap.add_argument("--no-subsample", action="store_true")
    # ── VirPhyKit 口径 RRT (树注释法); 2026-09-16 接入 ──
    ap.add_argument("--rrt-tree", default=None,
                    help="RRT 树注释法: 真实 MCC 树 (需 .set/.set.prob 注释)")
    ap.add_argument("--rrt-randomized", default=None,
                    help="RRT 树注释法: 随机化 MCC 树 目录/glob (上游 Example 用 20 棵)")
    ap.add_argument("--no-rrt", action="store_true",
                    help="跳过距离矩阵法 RRT (只想跑树注释法 RRT 时用)")
    # ── 事件级传播路径表 + 扩散统计 (2026-09-16 接入) ──
    ap.add_argument("--pathways", action="store_true",
                    help="事件级传播路径表/episode/LTL (复用 --mcc-tree; 输出 pathways/)")
    ap.add_argument("--dispersion", action="store_true",
                    help="扩散统计族+置换零模型 (依赖 --pathways 产物; 输出 dispersion/)")
    ap.add_argument("--direct-snp", type=float, default=2.0,
                    help="三级分类 direct 期望替换数阈值 (需按数据标定, 见 METHODS)")
    ap.add_argument("--indirect-snp", type=float, default=5.0)
    ap.add_argument("--censor-years", type=float, default=0.5, help="LTL 持久性删失窗(年)")
    ap.add_argument("--n-perm", type=int, default=200, help="置换零模型次数")
    ap.add_argument("--calibrate", default="manual", choices=["auto", "manual"],
                    help="auto=按本数据 E[S] 分布 q25/q75 自动标定置信阈值")
    # ── 遗传信息驱动降采样 (2026-09-16 接入; phymap FPS/clade 判据, numpy 重写) ──
    ap.add_argument("--genetic-subsample", action="store_true",
                    help="两段式降采样: (地点×时间窗)分层 + cell 内 FPS (输出 extended/genetic_subsampled.fasta)")
    ap.add_argument("--gs-mode", default="two_stage", choices=["two_stage", "clade", "none"])
    ap.add_argument("--gs-max-reps", type=int, default=3, help="cell 内最多保留代表数 (FPS)")
    ap.add_argument("--gs-min-snp", type=int, default=2, help="FPS 提前停止的最小 SNP 差")
    ap.add_argument("--gs-cutoff", type=int, default=5, help="clade 完全连锁 cutoff (SNP)")
    ap.add_argument("--gs-sensitivity", action="store_true",
                    help="附带给 全量 vs 子采样 的多样性对照 CSV (S/π/成对 SNP 分位数)")
    # ── 发表级传播图/动画 (2026-09-16 接入; 需 --pathways 先行或同批) ──
    ap.add_argument("--map", action="store_true",
                    help="静态传播图 PNG+PDF+SVG (默认纯矢量合规底图; --coords/geo_resolver 提供坐标)")
    ap.add_argument("--gif", action="store_true", help="附带事件级时间揭示动画 GIF")
    ap.add_argument("--gif-trail", type=int, default=-1,
                    help="-1=累积模式(默认,路线永久保留); 正整数=残影模式(只高亮最近 N 条)")
    ap.add_argument("--tree-map", action="store_true",
                    help="附带树-图联动双面板 (phymapr Tree|Map 图版替代; 需 --mcc-tree+--metadata)")
    ap.add_argument("--coords", default=None, help="location,lat,lon CSV (缺省走离线 geo_resolver)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    log = LogCollector()
    fa = args.fasta
    if not args.no_subsample:
        ss = smart_subsample(args.fasta, args.metadata,
                             os.path.join(args.outdir, "subsampled.fasta"), log=log)
        if ss.get("success") and ss.get("output_fasta"):
            fa = ss["output_fasta"]
    if not args.no_rrt:
        run_rrt(fa, args.metadata, args.outdir, log=log)
    else:
        log.emit("跳过距离矩阵法 RRT (--no-rrt)")
    if args.mcc_tree and os.path.exists(args.mcc_tree):
        run_tempmig(args.mcc_tree, args.metadata, args.outdir, log=log)

    # ── 事件级传播路径表 + 扩散统计 (+ 传播图/动画) (旁挂子分析, 不影响上游产物) ──
    if args.pathways or args.dispersion or args.map:
        if not (args.mcc_tree and os.path.exists(args.mcc_tree)):
            log.error("--pathways/--dispersion/--map 需要有效的 --mcc-tree → 跳过")
        else:
            from utils.transmission_paths import run_pathways
            pw_dir = os.path.join(args.outdir, "pathways")
            pr = run_pathways(args.mcc_tree, args.metadata, pw_dir, fasta=fa,
                              direct_snp=args.direct_snp, indirect_snp=args.indirect_snp,
                              censor_years=args.censor_years, calibrate=args.calibrate,
                              log=log)
            if pr.get("success") and args.dispersion:
                from utils.diffusion_stats import (find_posterior_trees,
                                                   markov_jumps_posterior, run_dispersion)
                run_dispersion(pr["branch_table"],
                               output_dir=os.path.join(args.outdir, "dispersion"),
                               n_perm=args.n_perm, log=log)
                _pt = find_posterior_trees(args.mcc_tree)
                if _pt:
                    mj = markov_jumps_posterior(_pt, log=log)
                    if mj.get("success"):
                        with open(os.path.join(args.outdir, "dispersion",
                                               "markov_jumps.json"), "w",
                                  encoding="utf-8") as _f:
                            json.dump(mj, _f, ensure_ascii=False, indent=2)
            # 发表级传播图 (+可选动画/树图联动) — 吃本批 pathways 产物
            if args.map or args.tree_map:
                if not pr.get("success"):
                    log.error("--map/--tree-map 依赖 --pathways 产物 → 跳过")
                else:
                    from utils.transmission_map import (draw_animation, draw_static_map,
                                                        draw_tree_map, resolve_coordinates)
                    import pandas as _pd
                    bt = _pd.read_csv(pr["branch_table"])
                    locs = sorted(set(bt["start_location"].dropna())
                                  | set(bt["end_location"].dropna()))
                    coords = resolve_coordinates(locs, args.coords, log=log)
                    missing = [l for l in locs if coords.get(l) is None]
                    if missing:
                        log.error(f"--map/--tree-map 跳过: 地点缺坐标 {missing}")
                    else:
                        if args.map:
                            draw_static_map(bt, coords, os.path.join(args.outdir, "transmission_map.png"),
                                            log=log)
                        if args.gif:
                            draw_animation(bt, coords,
                                           os.path.join(args.outdir, "transmission_map.gif"),
                                           trail=args.gif_trail,
                                           log=log)
                        if args.tree_map:
                            pw_csv = (pr["pathways"] if isinstance(pr.get("pathways"), str)
                                      else os.path.join(pw_dir, "pathways.csv"))
                            draw_tree_map(args.mcc_tree, args.metadata, coords,
                                          os.path.join(args.outdir, "transmission_tree_map.png"),
                                          pathways_csv=pw_csv, log=log)

    # ── 遗传信息驱动降采样 (独立于 --pathways; 输出 extended/) ──
    if args.genetic_subsample:
        from utils.genetic_subsample import compare_diversity, genetic_subsample
        gs_dir = os.path.join(args.outdir, "extended")
        os.makedirs(gs_dir, exist_ok=True)
        gs_fa = os.path.join(gs_dir, "genetic_subsampled.fasta")
        gr = genetic_subsample(fa, args.metadata, gs_fa, mode=args.gs_mode,
                               max_reps=args.gs_max_reps, min_snp_diff=args.gs_min_snp,
                               clade_cutoff=args.gs_cutoff, log=log)
        if args.gs_sensitivity:
            compare_diversity(fa, gs_fa,
                              out_csv=os.path.join(gs_dir, "genetic_subsample_sensitivity.csv"),
                              log=log)

    # ── RRT 树注释法 (VirPhyKit 同口径): 显式给真实树 + 随机化树才跑 ──
    if args.rrt_tree or args.rrt_randomized:
        if not args.rrt_tree:
            log.error("给了 --rrt-randomized 但没有 --rrt-tree (真实 MCC 树) → 跳过树注释法 RRT")
        elif not os.path.exists(args.rrt_tree):
            log.error(f"--rrt-tree 不存在: {args.rrt_tree} → 跳过树注释法 RRT")
        else:
            trees = collect_randomized_mccs(args.rrt_randomized,
                                            exclude=args.rrt_tree) if args.rrt_randomized else []
            if not trees:
                log.error(f"--rrt-randomized 未收集到任何树: {args.rrt_randomized} "
                          f"→ 跳过树注释法 RRT (空集会导致假性 PASS)")
            else:
                log.emit(f"树注释法 RRT: 真实 1 棵 + 随机化 {len(trees)} 棵")
                rrt_dir = os.path.join(args.outdir, "rrt_tree")
                rr = run_rrt_tree(args.rrt_tree, trees, rrt_dir, log=log)
                _dump_rrt_tree_summary(rr, args.outdir, log)
    print("GEO_ANALYSIS_DONE")


if __name__ == "__main__":
    main()
