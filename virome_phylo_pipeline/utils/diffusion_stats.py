#!/usr/bin/env python3
"""
utils/diffusion_stats.py — 扩散统计族（wavefront/速度/扩散系数）+ 置换零模型
==========================================================================
对 transmission_paths 产出的 branch_table（MCC 逐分支时空表）计算 seraphim
`spreadStatistics`（wnv_north_america R_script_analyses.r §6）的** MCC 描述版**：

  1. 分支扩散速度   — 大圆距离 / 时限（km/yr），均值 + 时长加权均值
  2. 扩散系数       — displacement² / (4·duration)（km²/yr，RRW 扩散参数 D 的经验对应）
  3. wavefront 距离时间序列 — 每个时间切片上"活跃谱系位置集合"的最大两两大圆距离
     （谱系位置沿分支线性插值，phylomovie pts_and_hpds_at_time.R:23-98 同法）

零模型（⚠ 与 seraphim 的 RRW 参数化零模型**不同**）：
  seraphim 零模型需 posterior log 的 sigma 矩阵重模拟（simulatorRRW1）；
  本模块用**路径端点配对置换**（打乱 start→end 配对、保留计数与时长），
  与管线已有 RRT/TEMPMIG 随机化检验同一精神：回答"观测的起终点配对是否强于随机"。
  两者都只服务于 MCC **描述层**，不作推断层（推断层是 BSSVS BF）。

自包含性：numpy/pandas + 标准库；坐标同 transmission_map（用户 CSV 或离线词典）。
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from utils.virphy_bridge import LogCollector
    from utils.transmission_map import haversine_km, resolve_coordinates
except ImportError:  # 允许直接运行
    from .virphy_bridge import LogCollector  # type: ignore
    from .transmission_map import haversine_km, resolve_coordinates  # type: ignore


# ─────────────────────────────────────────────────────────────────────
# 核心统计
# ─────────────────────────────────────────────────────────────────────

def pathway_velocity_table(pw: pd.DataFrame, coords: Dict) -> pd.DataFrame:
    """每条传播事件的位移/速度/扩散系数。"""
    rows = []
    for _, r in pw.iterrows():
        la0, lo0 = coords[r["start_location"]]
        la1, lo1 = coords[r["end_location"]]
        d = haversine_km(la0, lo0, la1, lo1)
        dur = float(r["duration"]) if pd.notna(r["duration"]) and r["duration"] > 0 else np.nan
        rows.append({
            "branch_id": r.get("branch_id", ""),
            "from": r["start_location"], "to": r["end_location"],
            "start_age": round(float(r["start_age"]), 4),
            "duration_yr": (round(dur, 4) if pd.notna(dur) else None),
            "displacement_km": round(d, 2),
            "velocity_km_yr": (round(d / dur, 3) if pd.notna(dur) else None),
            "diffusion_coeff_km2_yr": (round(d * d / (4 * dur), 2) if pd.notna(dur) else None),
        })
    return pd.DataFrame(rows)


def wavefront_series(
    branch_table: pd.DataFrame, coords: Dict, n_slices: int = 200,
) -> pd.DataFrame:
    """wavefront 距离时间序列（MCC 描述版）。

    对时间切片 t（根年龄→0），活跃谱系 = 跨越 t 的分支；谱系位置沿其分支
    在起止坐标间线性插值（phylomovie 同法）；wavefront = 活跃谱系位置两两
    大圆距离的最大值。root_distance = 活跃谱系到根位置的最大距离。
    """
    bt = branch_table.dropna(subset=["start_age", "end_age"])
    root_loc = None
    # 根 = age 最大的节点：branch_table 里 start_age 最大的 start_location
    if len(bt):
        root_loc = bt.loc[bt["start_age"].idxmax(), "start_location"]
    root_xy = coords.get(root_loc) if root_loc else None

    a_max = float(bt["start_age"].max())
    a_min = float(bt["end_age"].min())
    rows = []
    for t in np.linspace(a_max, a_min, n_slices):
        pts: List[Tuple[float, float]] = []
        for _, br in bt.iterrows():
            if br["start_age"] >= t >= br["end_age"]:
                span = br["start_age"] - br["end_age"]
                frac = 0.0 if span <= 0 else (br["start_age"] - t) / span
                p0, p1 = coords.get(br["start_location"]), coords.get(br["end_location"])
                if p0 and p1:
                    pts.append((p0[0] + frac * (p1[0] - p0[0]),
                                p0[1] + frac * (p1[1] - p0[1])))
        if len(pts) < 2:
            continue
        # 两两大圆距离最大值（点数 = 活跃谱系数，通常 < 50，O(n²) 可接受）
        wf = 0.0
        rd = 0.0
        for i in range(len(pts)):
            if root_xy:
                rd = max(rd, haversine_km(pts[i][0], pts[i][1], root_xy[0], root_xy[1]))
            for j in range(i + 1, len(pts)):
                wf = max(wf, haversine_km(pts[i][0], pts[i][1], pts[j][0], pts[j][1]))
        rows.append({"age": round(float(t), 4), "n_lineages": len(pts),
                     "wavefront_km": round(wf, 2),
                     "root_distance_km": round(rd, 2)})
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────
# 置换零模型
# ─────────────────────────────────────────────────────────────────────

def permutation_test(
    pw: pd.DataFrame, coords: Dict, n_perm: int = 200, seed: int = 42,
) -> Dict:
    """坐标重标注置换零模型：打乱"地点名 → 坐标"的分配（网络结构与时长不变）。

    ⚠ 设计说明：初版用"终点洗牌"，在单一源头网络（全部事件共享起点）下位移
    均值在置换下恒等（同一距离多重集），统计量无检验力——真实 gcva32 数据
    实测暴露（perm_mean == obs）。改为对**坐标指派**置换（RRT 随机化同精神）：
    零假设 = "观测传播网络的空间排布不强于随机地理指派"。
    统计量：平均位移 km 与平均速度 km/yr；单侧 p = P(perm >= obs)。
    RNG 种子显式记录（管线纪律：结果须可复现）。"""
    rng = np.random.default_rng(seed)
    locs = sorted({l for l in list(pw["start_location"]) + list(pw["end_location"])
                   if coords.get(l)})
    if len(locs) < 3:
        return {"n_perm": n_perm, "seed": seed, "note": "地点数<3，置换无意义，跳过",
                "p_displacement": None, "p_velocity": None}
    coord_pool = [coords[l] for l in locs]
    dur = pw["duration"].astype(float).fillna(0).values
    starts = pw["start_location"].tolist()
    ends = pw["end_location"].tolist()

    def stat(mapping: Dict[str, Tuple[float, float]]) -> Tuple[float, float]:
        ds, vs = [], []
        for s, e, d in zip(starts, ends, dur):
            if mapping.get(s) is None or mapping.get(e) is None:
                continue
            dist = haversine_km(mapping[s][0], mapping[s][1], mapping[e][0], mapping[e][1])
            ds.append(dist)
            if d > 0:
                vs.append(dist / d)
        return (float(np.mean(ds)) if ds else 0.0,
                float(np.mean(vs)) if vs else 0.0)

    obs_d, obs_v = stat(coords)
    perm_d, perm_v = [], []
    for _ in range(n_perm):
        pool = list(coord_pool)
        rng.shuffle(pool)
        mapping = dict(zip(locs, pool))
        d, v = stat(mapping)
        perm_d.append(d)
        perm_v.append(v)
    perm_d, perm_v = np.array(perm_d), np.array(perm_v)
    return {
        "n_perm": n_perm, "seed": seed,
        "null_scheme": "coordinate_relabeling",
        "observed": {"mean_displacement_km": round(obs_d, 2),
                     "mean_velocity_km_yr": round(obs_v, 4)},
        "perm_mean_displacement_km": round(float(perm_d.mean()), 2),
        "perm_mean_velocity_km_yr": round(float(perm_v.mean()), 4),
        "p_displacement": round(float((perm_d >= obs_d).mean()), 4),
        "p_velocity": round(float((perm_v >= obs_v).mean()), 4),
        "note": "坐标重标注置换（描述层零模型）；与 seraphim RRW 参数化零模型不同，不替代 BSSVS 推断",
    }


# ─────────────────────────────────────────────────────────────────────
# 编排
# ─────────────────────────────────────────────────────────────────────

def run_dispersion(
    branch_table_csv: str,
    coords_csv: Optional[str] = None,
    output_dir: str = "dispersion",
    n_slices: int = 200,
    n_perm: int = 200,
    seed: int = 42,
    log: Optional[LogCollector] = None,
) -> Dict:
    """扩散统计全套产物：velocity_table.csv / wavefront_series.csv / dispersion_summary.json"""
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    bt = pd.read_csv(branch_table_csv)
    locs = sorted(set(bt["start_location"].dropna()) | set(bt["end_location"].dropna()))
    coords = resolve_coordinates(locs, coords_csv, log=log)
    missing = [l for l in locs if coords.get(l) is None]
    if missing:
        # 2026-09-17 修: 缺坐标地点相关事件跳行+警告 (真实 PSTVd MCC 含 Unknown 节点)
        log.warning(f"地点缺坐标, 相关事件跳过: {missing}")
        bt = bt[(~bt["start_location"].isin(missing)) & (~bt["end_location"].isin(missing))]
        if bt.empty:
            raise SystemExit(f"全部事件都涉及缺坐标地点: {missing}")

    pw = bt[bt.get("location_change", pd.Series(dtype=bool)) == True]  # noqa: E712
    if pw.empty:
        raise SystemExit("branch_table 无地点变化分支")

    vt = pathway_velocity_table(pw, coords)
    vt_path = os.path.join(output_dir, "velocity_table.csv")
    vt.to_csv(vt_path, index=False)

    ws = wavefront_series(bt, coords, n_slices=n_slices)
    ws_path = os.path.join(output_dir, "wavefront_series.csv")
    ws.to_csv(ws_path, index=False)

    vt_valid = vt.dropna(subset=["velocity_km_yr"])
    summary = {
        "n_pathways": len(pw),
        "mean_velocity_km_yr": round(float(vt_valid["velocity_km_yr"].mean()), 4) if len(vt_valid) else None,
        "weighted_velocity_km_yr": (
            round(float(np.average(vt_valid["velocity_km_yr"],
                                   weights=vt_valid["duration_yr"].clip(lower=1e-9))), 4)
            if len(vt_valid) else None),
        "median_diffusion_coeff_km2_yr": (round(float(vt_valid["diffusion_coeff_km2_yr"].median()), 2)
                                          if len(vt_valid) else None),
        "max_wavefront_km": (round(float(ws["wavefront_km"].max()), 2) if len(ws) else None),
        "max_root_distance_km": (round(float(ws["root_distance_km"].max()), 2) if len(ws) else None),
    }
    null = permutation_test(pw, coords, n_perm=n_perm, seed=seed)
    summary["permutation_null"] = null
    summary_path = os.path.join(output_dir, "dispersion_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    log.emit(f"dispersion: {len(pw)} 事件, mean velocity = {summary['mean_velocity_km_yr']} km/yr, "
             f"max wavefront = {summary['max_wavefront_km']} km, "
             f"perm p(disp) = {null['p_displacement']}")
    return {"success": True, "velocity_table": vt_path, "wavefront_series": ws_path,
            "summary": summary_path, "summary_data": summary}


def main():
    import argparse
    ap = argparse.ArgumentParser(description="扩散统计族 + 置换零模型（MCC 描述层）")
    ap.add_argument("--branch-table", required=True, help="transmission_paths 的 branch_table.csv")
    ap.add_argument("--coords", default=None, help="location,lat,lon CSV（缺省走离线词典）")
    ap.add_argument("-o", "--outdir", default="dispersion")
    ap.add_argument("--n-slices", type=int, default=200)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    r = run_dispersion(args.branch_table, args.coords, args.outdir,
                       n_slices=args.n_slices, n_perm=args.n_perm, seed=args.seed)
    if not r["success"]:
        raise SystemExit("dispersion failed")
    print("DISPERSION_DONE")


if __name__ == "__main__":
    main()


# ─────────────────────────────────────────────────────────────────────
# 后验 Markov jumps (SpreaD3 同名统计的纯 Python 实现; 2026-09-17 新增)
# ─────────────────────────────────────────────────────────────────────

def markov_jumps_posterior(
    trees_file: str,
    trait: str = "Location",
    burnin_frac: float = 0.1,
    thin: int = 1,
    max_trees: int = 500,
    log: Optional[LogCollector] = None,
) -> Dict:
    """后验树上逐棵统计地点变迁边数 (Markov jumps) → 均值/中位/95% 区间。

    ⚠ 解析要点 (真实 GCVA 后验树实测): BEAST1 裸树行的节点注解在文本中按
    **后序 (postorder, 自左向右)** 出现, 且 Bio.Phylo NewickIO 在"节点注释+
    分支注释并存"时只保留一个 → 内部节点 Location 丢失 (实测 209 注解只解析
    出 1, 逐树 jump 恒 0)。故本函数用 Bio.Phylo 取**拓扑**, 用原始文本注解按
    postorder 对齐取**状态**。

    动机: MCC 单树上的 location-change 分支数是后验的**单个保守实现**
    (GCVA 全量 MCC=12 vs 后验均值≈20); 论文里"总迁移量"应引用后验均值+95%CI。
    """
    import re
    from io import StringIO
    from Bio import Phylo
    log = log or LogCollector()
    pat = re.compile(trait + r'="([^"]+)"')
    lines = [l for l in open(trees_file, encoding="utf-8", errors="replace")
             if l.startswith("tree ")]
    lines = lines[int(len(lines) * burnin_frac)::thin][:max_trees]
    counts, skipped = [], 0
    for l in lines:
        body = l.split("= ", 1)[1].strip() if "= " in l else l.strip()
        try:
            tr = Phylo.read(StringIO(body), "newick")
            raw = pat.findall(body)
            clades = list(tr.find_clades(order="postorder"))
            if len(raw) != len(clades):
                skipped += 1
                continue
            loc = {id(c): raw[i] for i, c in enumerate(clades)}
            n = 0
            for c in clades:
                for ch in c.clades:
                    a, b = loc.get(id(c)), loc.get(id(ch))
                    if a and b and a != b:
                        n += 1
            counts.append(n)
        except Exception:  # noqa: BLE001
            skipped += 1
    if not counts:
        return {"success": False, "error": "无可用树 (注解对齐全部失败)"}
    a = np.array(counts)
    out = {
        "success": True,
        "n_trees_used": len(counts), "n_trees_skipped": skipped,
        "burnin_frac": burnin_frac, "thin": thin, "trait": trait,
        "mean": round(float(a.mean()), 2),
        "median": int(np.median(a)),
        "hpd95": [int(np.percentile(a, 2.5)), int(np.percentile(a, 97.5))],
        "note": "后验 Markov jumps (每次 MCMC 采样的地点变迁边数); "
                "MCC 单树计数是其单个保守实现",
    }
    log.emit(f"Markov jumps: mean={out['mean']} median={out['median']} "
             f"95%CI={out['hpd95']} ({out['n_trees_used']} trees)")
    return out


def find_posterior_trees(mcc_tree: str) -> Optional[str]:
    """从 MCC 树路径探测同批后验树文件 (merged_sampled / merged / combined)。"""
    import glob as _glob
    d = os.path.dirname(os.path.abspath(mcc_tree))
    for pat in ("merged_sampled.trees", "merged.trees",
                "../combined_thinned.trees", "*.trees"):
        hits = sorted(_glob.glob(os.path.join(d, pat)))
        for h in hits:
            if os.path.abspath(h) != os.path.abspath(mcc_tree):
                return h
    return None
