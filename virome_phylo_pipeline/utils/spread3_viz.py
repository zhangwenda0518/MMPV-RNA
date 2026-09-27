#!/usr/bin/env python3
"""
spread3_viz.py — SpreaD3 交互式系统地理学可视化
================================================
生成 SpreaD3 兼容 JSON + 独立 HTML 地图。

SpreaD3 格式参考: Bielejec et al. (2016) SpreaD3: Interactive Visualization
of Spatiotemporal History and Trait Evolutionary Processes.
Mol Biol Evol 33:2167-2169.

功能:
  1. 从 BEAST log 提取迁移率矩阵 → SpreaD3 JSON
  2. 从 MCC 树提取节点位置+时间 → 时空散点
  3. 生成独立 HTML 交互式地图 (D3.js/Leaflet)

用法:
  python spread3_viz.py --log phylogeo_beast1.log --mcc mcc.tree \
      --locations Beijing,Neimenggu,Ningxia --output spread3.html
"""

import csv, json, os, sys, re, math
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector
# 2026-09-15 (审查 P2-9): BF 阈值/发散截断/公式 全管线唯一定义
from utils.virphy_bridge import (DEFAULT_BF_THRESHOLD, BF_INFINITE_CAP,
                                 bayes_factor, find_indicator_columns)


# ═══════════════════════════════════════════════════════════════════
# SpreaD3 JSON Generator
# ═══════════════════════════════════════════════════════════════════

# 常用地点坐标 (lat, lon)。**只收真坐标**, 缺失一律 fail-loud (见 generate_spread3_json
# 的 E5 说明): 历史版本对未收录地点用字符和编坐标, 地图上凭空多一个点。
_DEFAULT_LAT_LON = {
    "Beijing": (39.9, 116.4),
    "Ningxia": (37.3, 106.0),
    "Neimenggu": (44.0, 113.0),
    "China": (35.0, 105.0),
    "Japan": (36.0, 138.0),
    "Korea": (37.0, 127.5),
    "India": (20.0, 78.0),
    "Europe": (50.0, 10.0),
    "USA": (40.0, -100.0),
    "South America": (-15.0, -60.0),
    "Africa": (0.0, 25.0),
    "Other": (30.0, 80.0),
}

def generate_spread3_json(
    log_file: str,
    mcc_tree_file: str,
    locations: List[str],
    output_dir: str,
    burnin_pct: float = 10.0,
    bf_threshold: float = DEFAULT_BF_THRESHOLD,
    lat_lon_map: Optional[Dict[str, Tuple[float, float]]] = None,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    Generate SpreaD3-compatible JSON from BEAST output.

    Parameters
    ----------
    log_file : str — BEAST .log file with location indicators
    mcc_tree_file : str — MCC tree (Newick with location annotations)
    locations : List[str] — ordered list of discrete location states
    output_dir : str
    burnin_pct : float
    bf_threshold : float — BF threshold for significant routes
    lat_lon_map : dict — {location: (lat, lon)} for map placement
    log : LogCollector

    Returns
    -------
    dict: {success, json_path, n_routes, n_significant}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "json_path": None, "n_routes": 0,
              "n_significant": 0, "coord_missing": [], "error": None}

    try:
        n_loc = len(locations)

        # ── 前置守卫 (2026-09-15, 审查 P2-11) ───────────────────────────
        # 历史坑: log 缺失时 _extract_migration_matrix 只 warning + 返回 [],
        # 上层照样把 n_routes=0 的 JSON 落盘并置 success=True —— "0 条迁移路线"
        # 被当成成功产物, 调用方无法与"真的没有显著路线"区分。
        if not log_file or not os.path.exists(log_file):
            result["error"] = f"BEAST log not found: {log_file}"
            log.error(result["error"])
            return result
        if n_loc < 2:
            result["error"] = f"Need >=2 locations for migration extraction, got {n_loc}"
            log.error(result["error"])
            return result

        # ── 地点坐标 (E5, 2026-09-16: fail-loud, 绝不编坐标) ────────────
        # 历史坑: 未收录地点用 `(30 + sum(ord(ch))%20, 100 + sum(ord(ch))%35)`
        # **编**一组坐标 —— 落点在中国中部附近, 与真实位置毫无关系, 且地图上
        # 看不出来是编的 (静默造假)。现在按 ①调用方表 → ②常用内置表 →
        # ③geo_resolver (中国省/市离线表 + 缓存) 三级解析; 仍失败 → 该地点
        # latitude/longitude = None, 并进 coord_missing 由调用方/地图显式列出。
        explicit = lat_lon_map is not None
        lookup = dict(_DEFAULT_LAT_LON)
        if explicit:
            lookup.update(lat_lon_map)
        try:
            from utils.geo_resolver import resolve_location
        except ImportError:
            resolve_location = None
        coords, missing = {}, []
        for loc in locations:
            c = lookup.get(loc)
            src = "caller" if explicit and loc in (lat_lon_map or {}) else "builtin"
            if c is None and resolve_location is not None:
                c = resolve_location(loc)
                src = "geo_resolver" if c else src
            if c is None:
                missing.append(loc)
            else:
                coords[loc] = (float(c[0]), float(c[1]))
                log.emit(f"  坐标 {loc}: ({c[0]:.2f}, {c[1]:.2f}) ← {src}")
        if missing:
            log.warning(f"  {len(missing)}/{n_loc} 个地点无坐标 → 地图上不画点 "
                        f"(不编造): {', '.join(missing[:8])}"
                        + ("…" if len(missing) > 8 else ""))
        result["coord_missing"] = missing
        if not coords:
            result["error"] = (
                "全部地点都解析不到坐标 → 地图无点可画。修复: 用 lat_lon_map 传坐标表, "
                "或在 locations 里用已知地名 (如 " + ", ".join(sorted(_DEFAULT_LAT_LON)[:6])
                + " …)")
            log.error(result["error"])
            return result

        # ── Build nodes (无坐标的显式写 null, 不静默换默认值) ──
        nodes = []
        for i, loc in enumerate(locations):
            c = coords.get(loc)
            nodes.append({
                "id": loc,
                "label": loc,
                "latitude": round(c[0], 4) if c else None,
                "longitude": round(c[1], 4) if c else None,
            })

        # ── Extract migration rates from BEAST log ──
        migrations = _extract_migration_matrix(log_file, locations, burnin_pct, bf_threshold, log)
        result["n_routes"] = len(migrations)
        result["n_significant"] = sum(1 for m in migrations if m.get("significant"))

        if not migrations:
            # 2026-09-15 (审查 P2-11): 0 路线不再静默算成功。
            result["error"] = ("log 中未解析到任何 BSSVS 迁移路线 (indicator/rate 列缺失 "
                               "或列名与 locations 不匹配) — SpreaD3 JSON 未生成")
            log.error(result["error"])
            return result

        # ── Extract MCC tree events ──
        mcc_events = []
        if mcc_tree_file and os.path.exists(mcc_tree_file):
            mcc_events = _parse_mcc_tree_events(mcc_tree_file, locations, log)
            log.emit(f"MCC tree: {len(mcc_events)} annotated nodes")
        elif not mcc_tree_file:
            # 2026-08-30 审计修复 B6: 空值旧行为是静默跳过 → 时空散点维度整体缺失无人知。
            log.warning("mcc_tree_file 为空: SpreaD3 时空散点 (tree events) 未生成, 仅迁移路线可用")
        else:
            # 2026-09-16: 给了路径但文件不存在时旧行为**两个分支都不进** → 静默无散点。
            log.warning(f"mcc_tree_file 不存在: {mcc_tree_file} → 时空散点未生成")

        # ── Build SpreaD3 JSON ──
        spread3 = {
            "nodes": nodes,
            "links": [
                {
                    # fail-loud: from_idx/to_idx 由 _extract_migration_matrix 保证存在
                    "source": m["from_idx"],
                    "target": m["to_idx"],
                    "bf": m["bf"],
                    "posterior": m["posterior"],
                    "rate": m.get("rate", 1.0),
                }
                for m in migrations
            ],
            "mcc_events": mcc_events,
            "metadata": {
                "title": "Bayesian Phylogeography — Migration Routes",
                "bf_threshold": bf_threshold,
                "n_locations": n_loc,
                "n_routes": result["n_routes"],
                "n_significant": result["n_significant"],
                "coord_missing": missing,
                "generated": datetime.now().isoformat(),
                "source_log": os.path.basename(log_file),
                "units": {"height": "时间 (与 BEAST log 的 treeModel.rootHeight 同单位)",
                          "rate": "迁移率 (BSSVS 后验均值, 见 log 中的 location.rates)"},
            }
        }

        json_path = os.path.join(output_dir, "spread3.json")
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(spread3, f, indent=2, ensure_ascii=False)
        result["json_path"] = json_path

        log.emit(f"SpreaD3 JSON: {json_path}")
        log.emit(f"  Routes: {result['n_routes']} total, {result['n_significant']} significant")
        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"SpreaD3 JSON failed: {e}")
        return result


def _extract_migration_matrix(
    log_file: str,
    locations: List[str],
    burnin_pct: float,
    bf_threshold: float,
    log: LogCollector,
) -> List[Dict]:
    """Extract migration rates and BFs from BEAST log."""

    if not os.path.exists(log_file):
        log.warning(f"Log file not found: {log_file}")
        return []

    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
        burnin_idx = int(len(df) * burnin_pct / 100.0)
        df = df.iloc[burnin_idx:]

        n_loc = len(locations)
        n_pairs = n_loc * (n_loc - 1)
        migrations = []

        # Find indicator and rate columns (列检测口径全管线统一, P2-27 补)
        ind_cols = find_indicator_columns(df.columns)
        rate_cols = [c for c in df.columns
                     if 'rates' in c.lower() and 'indicator' not in c.lower()
                     and 'location' in c.lower()]

        # 2026-08-30 审计修复 (BUG-2): 旧版纯位置索引假设列序==参数维度序，但 BEAST 1.10
        # 实际输出命名列 location.indicators.<from>.<to>，列序与参数顺序可能不一致
        # (字母序) → 路线张冠李戴。与 phylogeo_bridge.extract_migration_bf 统一：
        # 命名列解析优先，位置索引仅作 fallback。
        def _pos_pair(idx):
            from_i = idx // (n_loc - 1)
            remainder = idx % (n_loc - 1)
            to_i = remainder if remainder < from_i else remainder + 1
            return from_i, to_i

        def _numeric_pair(col):
            """location.indicators<k> (1-based) → 声明序第 k 对。

            2026-09-15 (审查 P2-27 补): 本管线 beast1_bridge 写的是**单个向量参数**
            `location.indicators` (dimension=n_pairs), log 列名因此是
            `location.indicators<k>` (数字), 命名正则命不中 → 旧版只能落到
            "出现顺序" 回退, 字母序 logger 下整体错位。与 beast_parser /
            phylogeo_bridge 的 idx_map 同口径: k-1 → row-major 跳对角。
            """
            m = re.search(r'(\d+)\s*$', str(col))
            if not m:
                return None
            k = int(m.group(1)) - 1
            if not (0 <= k < n_pairs):
                return None
            return _pos_pair(k)

        def _col_route(col, idx):
            """→ (src_name, dst_name, from_i, to_i); 无法映射时返回 None。

            2026-09-15 (审查 P2-11 补, 重大): 必须**同时**给出地点索引 ——
            下游 links 构建读的是 m["from_idx"]/m["to_idx"], 而旧版两个分支都只写
            from/to 名字 → generate_spread3_json 每次都 KeyError 'from_idx' 静默失败
            (实测: 该函数从未产出过 spread3.json; 备份版同样如此)。
            """
            m_col = re.match(r'location\.indicators\.(\w+)\.(\w+)', col,
                             re.IGNORECASE)
            if m_col:
                src, dst = m_col.group(1), m_col.group(2)
                if src in locations and dst in locations:
                    return src, dst, locations.index(src), locations.index(dst)
                log.warning(f"  跳过列名地点不在 locations 里的 indicator 列: {col}")
                return None
            fi, ti = _numeric_pair(col) or _pos_pair(idx)
            if not (0 <= fi < n_loc and 0 <= ti < n_loc) or fi == ti:
                log.warning(f"  跳过无法映射到地点对的 indicator 列: {col}")
                return None
            return locations[fi], locations[ti], fi, ti

        def _rate_col_for(src, dst):
            """按 (from, to) 找 rate 列, 大小写不敏感 (log 可能写 Location.rates)。"""
            want = f'location.rates.{src}.{dst}'.lower()
            for _c in rate_cols:
                if _c.lower() == want:
                    return _c
            return None

        # 命名列直接按名字索引，无名字的才走数字后缀 / 位置映射
        named_cols = [(c, i) for i, c in enumerate(ind_cols)
                      if re.match(r'location\.indicators\.\w+\.\w+', c,
                                  re.IGNORECASE)]
        col_routes = []
        if len(named_cols) == len(ind_cols) and len(ind_cols) == n_pairs:
            for c, _ in named_cols:
                r = _col_route(c, 0)
                if r is None:
                    continue
                src, dst, fi, ti = r
                rate_col = _rate_col_for(src, dst)
                mean_rate = float(df[rate_col].mean()) if rate_col else 0.0
                posterior = float(df[c].mean())
                col_routes.append((src, dst, fi, ti, posterior, mean_rate))
        else:
            for idx, c in enumerate(ind_cols[:n_pairs]):
                r = _col_route(c, idx)
                if r is None:
                    continue
                src_name, dst_name, fi, ti = r
                posterior = float(df[c].mean())
                mean_rate = float(df[rate_cols[idx]].mean()) if idx < len(rate_cols) else 0.0
                col_routes.append((src_name, dst_name, fi, ti, posterior, mean_rate))

        for src_name, dst_name, fi, ti, posterior, mean_rate in col_routes:
            # BF: posterior odds / prior odds (全管线唯一定义, 见 virphy_bridge)
            bf = bayes_factor(posterior)

            sig = bf >= bf_threshold

            migrations.append({
                "from": src_name,
                "to": dst_name,
                # 2026-09-15 (P2-11 补): 下游 links 直接读这两个键, 缺一即整个 JSON
                # 生成失败。必须写出 (不要改成 m.get 静默成 null)。
                "from_idx": fi,
                "to_idx": ti,
                "bf": round(bf, 2),
                "posterior": round(posterior, 4),
                "rate": round(mean_rate, 6),
                "significant": sig,
            })

        migrations.sort(key=lambda x: -x["bf"])
        return migrations

    except Exception as e:
        log.warning(f"Migration matrix extraction failed: {e}")
        return []


def _split_annot_parts(comment: str) -> List[str]:
    """按**顶层**逗号切 BEAST 注释块 —— 不切开 `{...}` 里的逗号。

    `[&name.set={"A","B"},name.set.prob={0.6,0.4},height=1.0]` 这种块里，
    朴素 `comment.split(',')` 会把集合/概率向量切碎 → 必须做括号深度感知切分。
    """
    parts, buf, depth, inq = [], [], 0, False
    for ch in comment:
        if ch == '"':
            inq = not inq
            buf.append(ch)
        elif inq:
            buf.append(ch)
        elif ch in '{[':
            depth += 1
            buf.append(ch)
        elif ch in '}]':
            depth -= 1
            buf.append(ch)
        elif ch == ',' and depth == 0:
            parts.append(''.join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        parts.append(''.join(buf))
    return [p.strip() for p in parts if p.strip()]


def _annot_dicts(comment: str):
    """把一个注释块拆成 (flat, sets, setprobs)。

    flat      : {'height': '1.0', 'max': 'JAP', 'posterior': '0.98', ...}
    sets      : {'max': ['JAP','SWM'], 'location.states': [...]}
    setprobs  : {'max': [0.9986,0.0014], ...}
    """
    flat, sets, setprobs = {}, {}, {}
    for part in _split_annot_parts(comment):
        if '=' not in part:
            continue
        key, val = part.split('=', 1)
        key, val = key.strip().strip('"'), val.strip()
        if key.endswith('.set.prob'):
            pf = key[:-len('.set.prob')]
            body = val.strip('{}')
            try:
                setprobs[pf] = [float(x) for x in body.split(',') if x.strip()]
            except ValueError:
                pass
        elif key.endswith('.set'):
            pf = key[:-len('.set')]
            body = val.strip('{}')
            sets[pf] = [x.strip().strip('"').strip("'")
                        for x in body.split(',') if x.strip()]
        else:
            flat[key] = val.strip('"')
    return flat, sets, setprobs


def _parse_mcc_tree_events(
    tree_file: str,
    locations: List[str],
    log: LogCollector,
) -> List[Dict]:
    """Parse MCC tree with location annotations into SpreaD3 events.

    2026-09-16 修复（原实现在**真实 MCC 树上恒返回 0 事件且不报错**，三重缺陷）：

    ① **读法**：旧版 `Phylo.read(f, "newick")` 读 NEXUS 格式的 `.tre/.tree`
       （BEAST/treeannotator 输出的标准格式，含 `#NEXUS` + `Begin trees`）
       会抛 `ValueError: There are multiple trees in this file`，被 `except`
       吞掉 → 直接 `return []`。现在按文件头探测格式（NEXUS 优先，取第一棵）。

    ② **节点过滤**：旧版 `if not clade.name: continue` —— BEAST Newick 里
       **内部节点无名**，而祖先状态恰恰标在内部节点上 → 全部被跳过。
       现在不因无名而跳过（无名节点合成 `node<N>` 标签）。

    ③ **注释键**：旧版只认 `location.states`/`location`/`Location` 三种**点状态**键。
       真实 MCC 树的离散性状是 `<前缀>.set={...}` + `<前缀>.set.prob={...}`
       （由 treeannotator 自动生成；前缀名因分析而异，实测有 `Region`、`max`、
       `location.states`）→ 现在按 `*.set`/`*.set.prob` **取集合内后验最大的已知地点**，
       并保留对点状态键的兼容（值落在 `locations` 里才认）。

    ④ 注释存在但一个地点都没解出来时给 warning（旧版 0 事件是静默的）。
    """
    try:
        blocks = None
        try:
            with open(tree_file, encoding='utf-8', errors='replace') as f:
                head = f.read(4096)
            fmt = "nexus" if '#NEXUS' in head.upper() else "newick"
        except Exception:
            fmt = "newick"
        from Bio import Phylo
        it = Phylo.parse(tree_file, fmt)
        tree = next(iter(it), None)
        del it
        if tree is None:
            log.warning(f"MCC 树里没有可解析的树: {os.path.basename(tree_file)}")
            return []
    except Exception as e:                                     # noqa: BLE001
        log.warning(f"MCC 树解析失败 ({os.path.basename(tree_file)}): {e}")
        return []

    loc_set = set(locations)
    events = []
    n_annot = 0
    n_idx = 0
    for clade in tree.find_clades():
        n_idx += 1
        height = clade.branch_length or 0
        posterior = None
        location, state_prob = None, None

        comment = getattr(clade, 'comment', None)
        if comment:
            n_annot += 1
            flat, sets, setprobs = _annot_dicts(comment.strip('[]'))
            # height / posterior
            try:
                if 'height' in flat:
                    height = float(flat['height'])
            except ValueError:
                pass
            try:
                if 'posterior' in flat:
                    posterior = float(flat['posterior'])
            except ValueError:
                pass

            # ① 集合法：<前缀>.set + <前缀>.set.prob → 取已知地点里后验最大的
            for pf, states in sets.items():
                probs = setprobs.get(pf) or []
                cand = [(s, probs[i] if i < len(probs) else 0.0)
                        for i, s in enumerate(states) if s in loc_set]
                if cand:
                    location, state_prob = max(cand, key=lambda x: x[1])
                    break
            # ② 点状态法（兼容旧格式）：值正好是已知地点
            if location is None:
                for key, val in flat.items():
                    if key in ('height', 'posterior') or val.endswith('.prob'):
                        continue
                    v = val
                    if v in loc_set:
                        location = v
                    elif '+' in v:                # 平局被 treeannotator 拼成 "A+B"
                        hit = [x for x in v.split('+') if x in loc_set]
                        if len(hit) == 1:
                            location = hit[0]
                    if location:
                        try:
                            state_prob = float(flat.get(key + '.prob', ''))
                        except (TypeError, ValueError):
                            state_prob = None
                        break

        if location:
            events.append({
                "node": (clade.name or ("node%d" % n_idx))[:30],
                "location": location,
                "location_idx": locations.index(location),
                "height": round(height, 6),
                "posterior": round(posterior if posterior is not None
                                   else (state_prob if state_prob is not None else 1.0), 4),
                "state_prob": (round(state_prob, 4) if state_prob is not None else None),
            })

    if n_annot and not events:
        log.warning(f"MCC 树有 {n_annot} 个带注释节点，但没解出任何已知地点 "
                    f"(locations={locations[:6]}{'...' if len(locations) > 6 else ''}) "
                    f"→ 检查地点编码是否与注释里的状态名一致: "
                    f"{os.path.basename(tree_file)}")
    return events



# ═══════════════════════════════════════════════════════════════════
# Standalone HTML Interactive Map
# ═══════════════════════════════════════════════════════════════════

_MAP_CSS = """
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
         background: #f5f6fa; color: #2c3e50; }
  .wrap { max-width: 1180px; margin: 0 auto; padding: 18px; }
  h1 { font-size: 19px; margin-bottom: 2px; }
  .meta { color: #888; font-size: 12px; margin-bottom: 12px; }
  .card { background: #fff; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,.07);
          padding: 14px; margin-bottom: 14px; }
  .banner { background: #fdf6e3; border-left: 3px solid #f39c12; color: #7f5a04;
            font-size: 12px; line-height: 1.7; padding: 8px 12px; border-radius: 4px;
            margin-bottom: 10px; }
  .banner.err { background: #fdecea; border-color: #e74c3c; color: #922b21; }
  table { width: 100%; border-collapse: collapse; font-size: 12px; }
  th { background: #2c3e50; color: #fff; text-align: left; padding: 6px 10px; }
  td { padding: 5px 10px; border-bottom: 1px solid #eee; font-family: monospace; }
  .sig { color: #e74c3c; font-weight: 700; }
  .legend { display: flex; gap: 16px; flex-wrap: wrap; font-size: 12px; color: #555;
            margin-top: 8px; }
  .swatch { display: inline-block; width: 22px; height: 0; border-top-width: 3px;
            border-top-style: solid; vertical-align: middle; margin-right: 5px; }
  .foot { color: #999; font-size: 11px; text-align: center; margin: 18px 0 6px; }
  svg { display: block; width: 100%; height: auto; background: #fff; }
"""

# 地图配色 (与 phylogeo_figures.NATURE_PALETTE 同族)
_MAP_COLORS = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
               '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85']


def _esc(s) -> str:
    import html as _h
    return _h.escape(str(s), quote=True)


def _offline_svg_html(nodes, links, all_links, metadata, title) -> str:
    """离线自绘等距圆柱 (equirectangular) 地图 —— 纯 SVG, 零 JS / 零 CDN / 零底图。

    E4 (2026-09-16): 旧版默认走 Leaflet + CartoDB 瓦片, 离线/内网直接白屏 (仅有一个
    客户端 JS 兜底分支)。现在离线自绘是**默认**, 且在服务端生成 —— 连 JS 都不需要,
    打开即见。经纬网 + 弧线箭头 + BF 标度 + 缺坐标清单全部画进去。
    """
    pts = {n['id']: (float(n['latitude']), float(n['longitude']))
           for n in nodes if n.get('latitude') is not None}
    missing = metadata.get('coord_missing') or \
        [n['label'] for n in nodes if n.get('latitude') is None]
    thr = metadata.get('bf_threshold', 5)

    W, H = 980, 520
    pad = 56
    if pts:
        lats = [p[0] for p in pts.values()]
        lons = [p[1] for p in pts.values()]
        lat0, lat1 = min(lats), max(lats)
        lon0, lon1 = min(lons), max(lons)
        # 单点/极小范围: 给一个最小视野, 避免除零
        if lat1 - lat0 < 2:
            lat0, lat1 = lat0 - 4, lat1 + 4
        if lon1 - lon0 < 2:
            lon0, lon1 = lon0 - 4, lon1 + 4
        lat0, lat1 = lat0 - 2, lat1 + 2
        lon0, lon1 = lon0 - 2, lon1 + 2
    else:  # pragma: no cover - 上面已保证 coords 非空
        lat0, lat1, lon0, lon1 = -60, 80, -170, 170

    def X(lon):
        return pad + (lon - lon0) / (lon1 - lon0) * (W - 2 * pad)

    def Y(lat):
        return H - pad - (lat - lat0) / (lat1 - lat0) * (H - 2 * pad)

    svg = [f'<svg viewBox="0 0 {W} {H}" preserveAspectRatio="xMidYMid meet" '
           f'role="img" aria-label="{_esc(title)}">']
    svg.append('<defs>'
               '<marker id="ah-sig" viewBox="0 0 10 10" refX="9" refY="5" '
               'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
               f'<path d="M 0 0 L 10 5 L 0 10 z" fill="#e74c3c"/></marker>'
               '<marker id="ah-ns" viewBox="0 0 10 10" refX="9" refY="5" '
               'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
               f'<path d="M 0 0 L 10 5 L 0 10 z" fill="#95a5a6"/></marker>'
               '</defs>')
    svg.append(f'<rect x="0" y="0" width="{W}" height="{H}" fill="#fbfcfd" '
               f'stroke="#e8e8ef"/>')

    # 经纬网 (5° 起, 自适应步长)
    step = 5 if (lat1 - lat0) <= 40 else 10
    g = (int(lat0 / step) + 1) * step
    while g < lat1:
        y = Y(g)
        svg.append(f'<line x1="{pad}" y1="{y:.1f}" x2="{W-pad}" y2="{y:.1f}" '
                   f'stroke="#eef1f4"/><text x="{pad-6}" y="{y+3:.1f}" font-size="9" '
                   f'fill="#aab" text-anchor="end">{g}°N</text>')
        g += step
    g = (int(lon0 / step) + 1) * step
    while g < lon1:
        x = X(g)
        svg.append(f'<line x1="{x:.1f}" y1="{pad}" x2="{x:.1f}" y2="{H-pad}" '
                   f'stroke="#eef1f4"/><text x="{x:.1f}" y="{H-pad+14}" font-size="9" '
                   f'fill="#aab" text-anchor="middle">{g}°E</text>')
        g += step

    # 弧线 (先画非显著, 显著在上)
    bf_max = max([float(l.get('bf', 0) or 0) for l in all_links] or [1.0]) or 1.0
    for l in sorted(all_links, key=lambda x: float(x.get('bf', 0) or 0)):
        s, t = l.get('source'), l.get('target')
        if s not in pts or t not in pts or s == t:
            continue
        (la1, lo1), (la2, lo2) = pts[s], pts[t]
        x1, y1, x2, y2 = X(lo1), Y(la1), X(lo2), Y(la2)
        bf = float(l.get('bf', 0) or 0)
        sig = bf >= thr
        w = 1.0 + 4.0 * (bf / bf_max) ** 0.5
        # 轻微弓形, 避免双向重叠
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        dx, dy = x2 - x1, y2 - y1
        bow = 0.12
        cx, cy = mx - dy * bow, my + dx * bow
        names = {n['id']: n['label'] for n in nodes}
        tip = (f"{_esc(names.get(s, s))} → {_esc(names.get(t, t))}  "
               f"BF={bf:.1f}  posterior={float(l.get('posterior', 0) or 0):.3f}  "
               f"rate={float(l.get('rate', 0) or 0):.4g}")
        # 2026-09-16 兼容修复: f-string 表达式含转义引号是 3.12+ 语法,
        # 服务器 Python 3.10 报 SyntaxError → 提到表达式外 (语义不变)
        dash_attr = "" if sig else "stroke-dasharray='5,4'"
        svg.append(f'<path d="M {x1:.1f} {y1:.1f} Q {cx:.1f} {cy:.1f} {x2:.1f} {y2:.1f}" '
                   f'fill="none" stroke="{"#e74c3c" if sig else "#95a5a6"}" '
                   f'stroke-width="{w:.2f}" opacity="{0.85 if sig else 0.45}" '
                   f'{dash_attr} '
                   f'marker-end="url(#{"ah-sig" if sig else "ah-ns"})">'
                   f'<title>{tip}</title></path>')

    # 节点 + 标签 (带避让: 标签在点上方, 边缘处内收)
    for i, n in enumerate(nodes):
        if n['id'] not in pts:
            continue
        la, lo = pts[n['id']]
        x, y = X(lo), Y(la)
        col = _MAP_COLORS[i % len(_MAP_COLORS)]
        anchor = 'middle'
        lx = x
        if x < pad + 40:
            anchor, lx = 'start', x + 10
        elif x > W - pad - 40:
            anchor, lx = 'end', x - 10
        svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{col}" '
                   f'stroke="#fff" stroke-width="2"><title>{_esc(n["label"])} '
                   f'({la:.2f}, {lo:.2f})</title></circle>')
        svg.append(f'<text x="{lx:.1f}" y="{y-11:.1f}" font-size="11" fill="#2c3e50" '
                   f'text-anchor="{anchor}" font-weight="600">{_esc(n["label"])}</text>')
    svg.append(f'<text x="{pad}" y="{H-8}" font-size="10" fill="#9aa">'
               f'equirectangular projection · 自绘 SVG, 无底图 / 无 CDN / 无 JS</text>')
    svg.append('</svg>')

    # 路线表
    rows = []
    for l in sorted(all_links, key=lambda x: -float(x.get('bf', 0) or 0))[:15]:
        s, t = l.get('source'), l.get('target')
        names = {n['id']: n['label'] for n in nodes}
        bf = float(l.get('bf', 0) or 0)
        rows.append(f'<tr><td>{_esc(names.get(s, s))}</td><td>{_esc(names.get(t, t))}</td>'
                    f'<td class="{"sig" if bf >= thr else ""}">{bf:.1f}</td>'
                    f'<td>{float(l.get("posterior", 0) or 0):.3f}</td>'
                    f'<td>{float(l.get("rate", 0) or 0):.4g}</td></tr>')
    more = (f'<p class="meta" style="margin-top:6px">共 {len(all_links)} 条路线, '
            f'表中按 BF 取前 15。</p>') if len(all_links) > 15 else ''

    banner = ''
    if missing:
        banner = (f'<div class="banner">⚠ {len(missing)} 个地点没有坐标, '
                  f'<b>未在地图上画点</b>（不编造坐标）：{_esc("、".join(missing))}<br>'
                  f'修复: 调用方传 lat_lon_map, 或用 geo_resolver 认识的地名。</div>')
    if not pts:
        banner += '<div class="banner err">没有任何地点可画 (全部缺坐标)。</div>'

    return f'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{_esc(title)} — Phylogeographic Map (offline)</title>
<style>{_MAP_CSS}</style>
</head>
<body>
<div class="wrap">
  <h1>{_esc(title)}</h1>
  <p class="meta">{metadata.get('n_locations', '?')} locations ·
     {metadata.get('n_significant', '?')} significant routes (BF≥{thr}) ·
     {metadata.get('n_routes', '?')} total · 源: {_esc(metadata.get('source_log', '?'))}</p>
  {banner}
  <div class="card">{''.join(svg)}
    <div class="legend">
      <span><span class="swatch" style="border-color:#e74c3c"></span>显著路线 (BF≥{thr})</span>
      <span><span class="swatch" style="border-color:#95a5a6;border-top-style:dashed"></span>非显著 (虚线)</span>
      <span>线宽 ∝ √BF (最大 BF={bf_max:.1f})</span>
      <span>悬停看 BF / 后验 / 速率</span>
    </div>
  </div>
  <div class="card">
    <table><thead><tr><th>From</th><th>To</th><th>BF</th><th>Posterior</th><th>Rate</th></tr></thead>
    <tbody>{''.join(rows) or '<tr><td colspan="5">无路线</td></tr>'}</tbody></table>
    {more}
  </div>
  <p class="foot">MMPV-RNA virome_phylo_pipeline · 单文件离线地图 ·
     {metadata.get('generated', '')}</p>
</div>
</body>
</html>'''


def _leaflet_html(nodes, links, all_links, metadata, title,
                  tile_url: Optional[str] = None,
                  tile_attribution: str = '') -> str:
    """在线版 (Leaflet)。⚠ 底图合规 (2026-09-16 整改):

    默认 ``tile_url=None`` → **不加载任何第三方瓦片** (仅点位/连线矢量)。
    旧版硬编码 CartoDB/OSM 瓦片属无审图号底图, 不符合中国发表要求 (borrow 报告 §4.2)。
    需要底图时显式传合规源 XYZ 模板, 例如天地图:
        tile_url='https://t{s}.tianditu.gov.cn/DataServer?T=vec_w&x={x}&y={y}&l={z}&tk=<你的密钥>'
        (天地图 {s} 为 0-7 子域; 高德/腾讯同理自备密钥)
    离线场景请用 offline=True (自绘 SVG, 零外链)。
    """
    nodes_js = json.dumps(nodes, ensure_ascii=False)
    links_js = json.dumps(links, ensure_ascii=False)
    all_links_js = json.dumps(all_links, ensure_ascii=False)
    if tile_url:
        tile_js = (f"L.tileLayer('{tile_url}', "
                   f"{{attribution: '{tile_attribution}', maxZoom: 12}}).addTo(map);")
    else:
        tile_js = ("// 底图合规默认: 不加载第三方瓦片 (无审图号底图不可用于中国发表)。\n"
                   "// 需要底图请经 generate_interactive_map_html(..., tile_url=<合规XYZ模板>) 传入。")
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title} — Phylogeographic Map</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; }}
  #map {{ height: 100vh; width: 100%; }}
  .info-panel {{
    position: absolute; top: 10px; right: 10px; z-index: 1000;
    background: rgba(255,255,255,0.95); padding: 15px 20px;
    border-radius: 8px; box-shadow: 0 2px 12px rgba(0,0,0,0.15);
    max-width: 360px; font-size: 13px;
  }}
  .info-panel h2 {{ font-size: 16px; margin-bottom: 8px; color: #2c3e50; }}
  .info-panel table {{ width: 100%; border-collapse: collapse; }}
  .info-panel th, .info-panel td {{ padding: 4px 8px; text-align: left; border-bottom: 1px solid #eee; }}
  .info-panel th {{ font-size: 11px; color: #888; font-weight: 600; }}
  .info-panel .sig {{ color: #e74c3c; font-weight: bold; }}
  .legend {{
    position: absolute; bottom: 30px; left: 10px; z-index: 1000;
    background: rgba(255,255,255,0.9); padding: 10px 15px;
    border-radius: 6px; font-size: 12px;
  }}
</style>
</head>
<body>
<div id="map"></div>
<div class="info-panel" id="info">
  <h2>{title}</h2>
  <p style="color:#888;font-size:12px;">
    {metadata.get('n_locations', '?')} locations &middot;
    {metadata.get('n_significant', '?')} significant routes &middot;
    BF threshold: {metadata.get('bf_threshold', '?')}
  </p>
  <table id="routeTable">
    <tr><th>From</th><th>To</th><th>BF</th></tr>
  </table>
</div>
<div class="legend">
  <div>● Location &nbsp; → Migration route (width ∝ BF)</div>
</div>

<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script>
const nodes = {nodes_js};
const links = {links_js};
const allLinks = {all_links_js};
const usable = nodes.filter(n => n.latitude !== null && n.longitude !== null);
const pt = {{}};
usable.forEach(n => {{ pt[n.id] = n; }});

if (typeof L === 'undefined') {{
  document.addEventListener('DOMContentLoaded', () => {{
    document.getElementById('map').innerHTML =
      '<div style="padding:16px;color:#c0392b;font-weight:600">' +
      '⚠ Leaflet CDN 不可用 (离线) —— 请改用 offline=True 生成自绘 SVG 地图 ' +
      '(generate_interactive_map_html(..., offline=True))</div>';
  }});
}} else {{
const centerLat = usable.reduce((a,b) => a + b.latitude, 0) / (usable.length || 1);
const centerLon = usable.reduce((a,b) => a + b.longitude, 0) / (usable.length || 1);
const map = L.map('map').setView([centerLat, centerLon], 5);
{tile_js}

const colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf'];
const locColors = {{}};
usable.forEach((n, i) => {{ locColors[n.id] = colors[i % colors.length]; }});

usable.forEach(n => {{
  L.circleMarker([n.latitude, n.longitude], {{
    radius: 10, fillColor: locColors[n.id], color: '#333', weight: 2,
    opacity: 1, fillOpacity: 0.85
  }}).addTo(map).bindPopup(`<b>${{n.label}}</b><br>Lat: ${{n.latitude.toFixed(2)}}<br>Lon: ${{n.longitude.toFixed(2)}}`);
}});

const maxBF = Math.max(...links.map(l => l.bf || 0), 1);
links.forEach(l => {{
  const src = pt[nodes[l.source] && nodes[l.source].id];
  const tgt = pt[nodes[l.target] && nodes[l.target].id];
  if (!src || !tgt) return;
  const bf = l.bf || 1;
  const color = bf >= {metadata.get('bf_threshold', 5)} ? '#e74c3c' : '#95a5a6';
  const line = L.polyline([[src.latitude, src.longitude], [tgt.latitude, tgt.longitude]], {{
    color, weight: Math.max(1, Math.log(bf) * 2),
    opacity: Math.min(0.9, 0.3 + bf / maxBF * 0.6),
    dashArray: bf >= {metadata.get('bf_threshold', 5)} ? null : '5,5'
  }}).addTo(map);
  line.bindPopup(`${{src.label}} → ${{tgt.label}}<br>BF: ${{bf.toFixed(1)}}<br>Posterior: ${{(l.posterior||0).toFixed(3)}}`);
}});

const table = document.getElementById('routeTable');
allLinks.slice(0, 15).forEach(l => {{
  const src = nodes[l.source], tgt = nodes[l.target];
  if (!src || !tgt) return;
  const bf = l.bf || 0;
  const cls = bf >= {metadata.get('bf_threshold', 5)} ? 'sig' : '';
  const row = table.insertRow();
  row.innerHTML = `<td>${{src.label}}</td><td>${{tgt.label}}</td><td class="${{cls}}">${{bf.toFixed(1)}}</td>`;
}});

if (usable.length) map.fitBounds(L.latLngBounds(usable.map(n => [n.latitude, n.longitude])).pad(0.2));
}} // end Leaflet-available branch
</script>
</body>
</html>'''


def generate_interactive_map_html(
    spread3_json_path: str,
    output_path: str,
    title: str = "Virus Phylogeography",
    offline: bool = True,
    tile_url: Optional[str] = None,
    tile_attribution: str = '',
    log: Optional[LogCollector] = None,
) -> str:
    """
    SpreaD3 JSON → 独立 HTML 地图。

    offline=True (默认): 服务端自绘等距圆柱投影 SVG —— 无 CDN / 无瓦片 / 无 JS,
                         离线必可见, 并在写盘后做「零外链」自检 (strict)。
    offline=False:        Leaflet 在线版。⚠ 底图合规 (2026-09-16): 默认 tile_url=None
                         **不加载任何第三方瓦片** (仅矢量点位/连线); 旧版硬编码
                         CartoDB/OSM 已移除——无审图号底图不符合中国发表要求。
                         需要底图时传合规 XYZ 模板 (天地图/高德/腾讯, 自备密钥), 如:
                         tile_url='https://t{s}.tianditu.gov.cn/DataServer?T=vec_w&x={x}&y={y}&l={z}&tk=<密钥>'

    Returns path to HTML file.
    """
    log = log or LogCollector()

    with open(spread3_json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    nodes = data["nodes"]
    links = data["links"]
    metadata = data.get("metadata", {})

    # Only show significant links
    sig_links = [l for l in links if l.get("bf", 0) >= metadata.get("bf_threshold", 5)]
    all_links = sorted(links, key=lambda x: -x.get("bf", 0))

    if offline:
        html = _offline_svg_html(nodes, all_links, all_links, metadata, title)
        try:
            from utils.self_contained_html import write_html
            write_html(html, output_path, log=log, strict=True,
                       what="离线系统地理地图")
        except ImportError:
            os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html)
            log.warning("self_contained_html 不可用, 未做外链自检")
        if metadata.get("coord_missing"):
            log.warning(f"  地图缺坐标地点 (未画点): {', '.join(metadata['coord_missing'][:8])}")
        log.emit(f"Interactive map (offline SVG): {output_path}")
        return output_path

    html = _leaflet_html(nodes, sig_links if sig_links else all_links[:10],
                         all_links, metadata, title,
                         tile_url=tile_url, tile_attribution=tile_attribution)
    try:
        from utils.self_contained_html import write_html
        write_html(html, output_path, log=log, what="Leaflet 地图")
    except ImportError:
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
    log.emit(f"Interactive map (Leaflet, 需联网): {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════
# Result Dashboard
# ═══════════════════════════════════════════════════════════════════

def build_phylogeo_report(
    spread3_json: str,
    summary_csv: str,
    output_dir: str,
    log: Optional[LogCollector] = None,
) -> str:
    """
    Generate a comprehensive HTML report combining all analysis results.

    Returns path to report HTML.
    """
    log = log or LogCollector()

    # Read data
    with open(spread3_json, 'r', encoding='utf-8') as f:
        spread3 = json.load(f)

    summary_rows = []
    if summary_csv and os.path.exists(summary_csv):
        with open(summary_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            summary_rows = list(reader)

    metadata = spread3.get("metadata", {})
    nodes = spread3.get("nodes", [])
    links = spread3.get("links", [])
    sig_links = [l for l in links if l.get("bf", 0) >= metadata.get("bf_threshold", 5)]

    # Build tables
    # 2026-09-16 (E5 联动): 坐标不可解析时 latitude/longitude = None, 旧的
    # f"{None:.1f}" 会直接 TypeError 崩掉整份报告。坐标缺失必须显式可见 (项目铁律)。
    _coord_missing = set(metadata.get("coord_missing") or [])

    def _fmt_coord(v):
        return f"{v:.1f}" if isinstance(v, (int, float)) else "—"

    loc_rows = ""
    for n in nodes:
        absent = (n.get('latitude') is None or n.get('longitude') is None)
        note = ' <span style="color:#e67e22">(缺坐标, 未画点)</span>' if absent else ""
        loc_rows += (f"<tr><td>{_esc(n['label'])}{note}</td>"
                     f"<td>{_fmt_coord(n.get('latitude'))}</td>"
                     f"<td>{_fmt_coord(n.get('longitude'))}</td></tr>")

    coord_note = ""
    if _coord_missing:
        coord_note = (f'<p class="meta" style="color:#e67e22">⚠ {len(_coord_missing)} 个地点解析不到经纬度'
                      f' (未画点): {_esc(", ".join(sorted(_coord_missing)))} — '
                      f'补法: 传 lat_lon_map 或往 utils/geo_resolver 内置表加条目</p>')

    route_rows = ""
    for l in sorted(links, key=lambda x: -x.get("bf", 0))[:20]:
        src = nodes[l["source"]]["label"] if l["source"] < len(nodes) else "?"
        tgt = nodes[l["target"]]["label"] if l["target"] < len(nodes) else "?"
        bf = l.get("bf", 0)
        sig = "✅" if bf >= metadata.get("bf_threshold", 5) else ""
        route_rows += (f"<tr><td>{_esc(src)}</td><td>{_esc(tgt)}</td>"
                       f"<td>{bf:.1f}</td><td>{sig}</td></tr>")

    summary_html = ""
    if summary_rows:
        cols = list(summary_rows[0].keys())
        summary_html += "<tr>" + "".join(f"<th>{c}</th>" for c in cols) + "</tr>"
        for row in summary_rows:
            summary_html += "<tr>" + "".join(f"<td>{row.get(c,'')}</td>" for c in cols) + "</tr>"

    report = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Phylogeography Report — {metadata.get("title", "BEAST Analysis")}</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; background: #f5f6fa; color: #2c3e50; }}
  .container {{ max-width: 1100px; margin: 0 auto; padding: 30px; }}
  h1 {{ font-size: 24px; margin-bottom: 5px; }}
  h2 {{ font-size: 18px; margin: 25px 0 12px; color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 5px; }}
  .meta {{ color: #888; font-size: 14px; margin-bottom: 20px; }}
  table {{ width: 100%; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 4px rgba(0,0,0,0.08); margin-bottom: 20px; }}
  th {{ background: #34495e; color: white; padding: 10px 12px; text-align: left; font-weight: 600; font-size: 13px; }}
  td {{ padding: 8px 12px; border-bottom: 1px solid #ecf0f1; font-size: 13px; }}
  tr:hover {{ background: #f8f9fa; }}
  .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-bottom: 25px; }}
  .stat {{ background: white; padding: 18px; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.08); text-align: center; }}
  .stat .value {{ font-size: 28px; font-weight: 700; color: #3498db; }}
  .stat .label {{ font-size: 12px; color: #888; margin-top: 5px; }}
  .footer {{ text-align: center; color: #bbb; font-size: 12px; margin-top: 40px; }}
</style>
</head>
<body>
<div class="container">
<h1>📊 Bayesian Phylogeography Report</h1>
<p class="meta">{metadata.get("title", "")} — Generated {metadata.get("generated", "")}</p>

<div class="stats">
  <div class="stat"><div class="value">{metadata.get("n_locations", "?")}</div><div class="label">Locations</div></div>
  <div class="stat"><div class="value">{metadata.get("n_routes", "?")}</div><div class="label">Total Routes</div></div>
  <div class="stat"><div class="value" style="color:#e74c3c">{metadata.get("n_significant", "?")}</div><div class="label">Significant (BF≥{metadata.get("bf_threshold", "?")})</div></div>
  <div class="stat"><div class="value">{len(nodes)}</div><div class="label">Discrete States</div></div>
</div>

<h2>🗺 Locations</h2>
{coord_note}
<table><tr><th>Location</th><th>Latitude</th><th>Longitude</th></tr>{loc_rows}</table>

<h2>✈ Migration Routes (top 20 by BF)</h2>
<table><tr><th>From</th><th>To</th><th>Bayes Factor</th><th>Significant</th></tr>{route_rows}</table>
'''

    if summary_rows:
        report += f'<h2>📋 Analysis Summary</h2><table>{summary_html}</table>'

    report += f'''
<div class="footer">
  MMPV-RNA virome_phylo_pipeline &middot; BEAST Phylogeography &middot; {datetime.now().strftime("%Y-%m-%d %H:%M")}
</div>
</div>
</body>
</html>'''

    report_path = os.path.join(output_dir, "phylogeo_report.html")
    # 2026-09-16: 旧版直接 open(), output_dir 不存在时 FileNotFoundError 炸掉整条
    # 可视化收尾 (调用方通常已在, 但独立调用/新目录必然踩到)。与 generate_spread3_json 一致。
    os.makedirs(output_dir, exist_ok=True)
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)

    log.emit(f"Report: {report_path}")
    return report_path
