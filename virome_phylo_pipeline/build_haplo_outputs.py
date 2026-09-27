#!/usr/bin/env python3
"""build_haplo_outputs.py — pypopart 全套绘图/分析输出整合

输出 (每套数据):
  1. 静态发表级网络图 (SVG/PDF/PNG)  — create_publication_figure
  2. 交互式 HTML 网络                  — InteractiveNetworkPlotter
  3. 地理地图可视化                    — geo-visualize (lat/lon)
  4. 网络统计                          — pypopart analyze

用法:
  python build_haplo_outputs.py --virus PSTVd|GCVA
"""
import argparse
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.seq_ids import base_sample_id
from utils.dataset_config import (dataset as get_dataset, colors as get_colors,
                                  provinces as get_provinces, env as get_env,
                                  repo as get_repo, known_dataset_names)

# 优先取 datasets.yaml 的 env.pypopart_path（默认已是 {repo}/biosoft/pypopart/src）；
# 兜底同样按仓库定位，不写死本机绝对路径
PYPOPART_PATH = get_env('pypopart_path') or os.path.join(get_repo(), 'biosoft', 'pypopart', 'src')
if os.path.isdir(PYPOPART_PATH) and PYPOPART_PATH not in sys.path:
    sys.path.insert(0, PYPOPART_PATH)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pypopart.io import load_alignment
from pypopart.algorithms.mjn import MedianJoiningNetwork
from pypopart.visualization.static_plot import create_publication_figure
from pypopart.visualization.interactive_plot import InteractiveNetworkPlotter

# 采样点/配色从 datasets.yaml 读取（新病毒自行扩展）
def _prov_coords():
    return get_provinces() or {}


def _loc_colors():
    return get_colors('location') or {}


def _host_colors():
    return get_colors('host') or {}


def _host_meta_path(meta_path):
    """宿主值的来源表。

    2026-09-16 第三轮: `host` 已从地理表独立成 `data/host.csv` (宿主通道),
    故 `--color-by host` **不能**再读 `cfg['meta']` (那里已无 host 列)。
    优先同目录 `host.csv`; 没有则退回 `meta_path` 本身 (兼容旧版地理表)。
    """
    cand = os.path.join(os.path.dirname(meta_path), "host.csv")
    return cand if os.path.exists(cand) else meta_path


def _host_col_of(path):
    """该表的 host 列名 (无则 None)。"""
    try:
        with open(path, encoding='utf-8-sig', errors='replace') as f:
            for c in (csv.DictReader(f).fieldnames or []):
                if (c or '').strip().lower() in ('host', 'host_species',
                                                 'host_scientific_name'):
                    return c
    except OSError:
        return None
    return None


def host_of(host):
    for h in _host_colors():
        if h in host:
            return h
    return 'Unknown'


def province_of(loc):
    for p in _prov_coords():
        if p in loc:
            return p
    return 'Unknown'


def _haplo_cfg(virus):
    """取数据集单倍型配置 (fasta/meta/out/recombinants)。"""
    d = get_dataset(virus)
    if d is None or 'haplo' not in d:
        raise SystemExit(f"数据集 {virus} 未配置 haplo 段 (检查 datasets.yaml)")
    return d['haplo']


def make_geo_meta(meta_path, out_path):
    """生成带 lat/lon 的 metadata (geo-visualize 输入)"""
    provs = _prov_coords()
    rows = []
    with open(meta_path, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            prov = province_of(row.get('location', ''))
            lat, lon = provs.get(prov, provs['Unknown'])
            # 2026-09-15: 与 data_collector / popgen 统一口径 (见 utils/seq_ids.py)
            sid = base_sample_id(row.get('name', ''))
            rows.append({'id': sid, 'latitude': lat, 'longitude': lon,
                         'population': prov, 'location': row.get('location', '')})
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['id', 'latitude', 'longitude', 'population', 'location'])
        w.writeheader()
        w.writerows(rows)
    return out_path


def build(virus, color_by='location'):
    cfg = _haplo_cfg(virus)
    os.makedirs(cfg['out'], exist_ok=True)
    print(f"=== {virus} (着色: {color_by}) ===")

    # 1. 构建网络
    aln = load_alignment(cfg['fasta'])
    network = MedianJoiningNetwork(distance_method='hamming', epsilon=0).build_network(aln)
    print(f"网络: {network.num_nodes} 节点, {network.num_edges} 边")

    # 1.5 设置 population (地点或宿主)
    colors_map = _loc_colors() if color_by == 'location' else _host_colors()
    # host 已独立成宿主通道 (data/host.csv) → 按着色维度选对来源表
    meta_loc = {}
    src_meta = cfg['meta'] if color_by == 'location' else _host_meta_path(cfg['meta'])
    if color_by == 'host' and _host_col_of(src_meta) is None:
        raise SystemExit(
            f"[build_haplo_outputs] --color-by host 需要宿主通道, 但 {src_meta} "
            f"没有 host 列。\n  host 已从地理表独立出来 → 请确认存在 "
            f"{os.path.join(os.path.dirname(cfg['meta']), 'host.csv')} "
            f"(由 phylo_pipeline 的通道拆分/治理产出)。")
    print(f"  metadata 来源: {src_meta}")
    with open(src_meta, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            # 2026-09-15: 与 data_collector / popgen 统一口径 (见 utils/seq_ids.py)
            sid = base_sample_id(row.get('name', ''))
            if color_by == 'location':
                meta_loc[sid] = province_of(row.get('location', ''))
            else:
                meta_loc[sid] = host_of(row.get('host', ''))
    for hap_id in network.nodes:
        hap = network.get_haplotype(hap_id)
        if hap is None:
            continue
        for sid in hap.sample_ids:
            parts = sid.split('/')
            base = base_sample_id(sid)
            pop = meta_loc.get(base, 'Unknown')
            hap.add_sample(sid, pop)

    # 2. 静态发表级图 (按 color_by 着色, 重组子三角标记)
    print("静态发表级图...")
    base = f"{cfg['out']}/{virus.lower()}_network_{color_by}_pub"
    import networkx as nx
    layout = nx.spring_layout(network._graph, seed=42)
    fig, ax = create_publication_figure(
        network, population_colors=colors_map,
        filename=base + '.png', layout=layout)
    # 叠加重组子三角标记 (含重组子的单倍型节点用红色三角描边)
    recombinants = set(cfg.get('recombinants', []) or [])
    if recombinants:
        recomb_nodes = []
        for hap_id in network.nodes:
            hap = network.get_haplotype(hap_id)
            if hap is None:
                continue
            for sid in hap.sample_ids:
                parts = sid.split('/')
                base_id = base_sample_id(sid)
                if base_id in recombinants:
                    recomb_nodes.append(hap_id)
                    break
        if recomb_nodes:
            nx.draw_networkx_nodes(
                network._graph, layout, nodelist=recomb_nodes,
                node_shape='^', node_color='none',
                edgecolors='red', linewidths=2.5, ax=ax)
            print(f"  重组子三角标记: {len(recomb_nodes)} 节点")
    fig.savefig(base + '.svg', bbox_inches='tight')
    fig.savefig(base + '.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"  → {base}.{{png,svg,pdf}}")

    # 3. 交互式 HTML
    print("交互式 HTML...")
    plotter = InteractiveNetworkPlotter(network)
    plotter.plot()
    plotter.add_population_legend(colors_map)
    plotter.save_html(f"{cfg['out']}/{virus.lower()}_network_{color_by}_interactive.html")
    print(f"  → {virus.lower()}_network_{color_by}_interactive.html")

    # 4. 地理地图 (matplotlib 省坐标 + 网络叠加)
    print("地理地图...")
    geo_meta = make_geo_meta(cfg['meta'], f"{cfg['out']}/geo_meta.csv")
    try:
        provs = _prov_coords()
        loc_colors = _loc_colors()
        fig, ax = plt.subplots(figsize=(9, 7))
        # 省坐标
        for prov, (lat, lon) in provs.items():
            ax.scatter(lon, lat, s=400, c=loc_colors.get(prov, '#999999'),
                       alpha=0.3, zorder=1)
            ax.text(lon, lat + 0.8, prov, ha='center', fontsize=8)
        # 网络节点 → 省份坐标 (单倍型所属省, 抖动避免重叠)
        import numpy as np
        rng = np.random.RandomState(7)
        node_pop = {}
        for n in network.nodes:
            hap = network.get_haplotype(n)
            if hap is None:
                continue
            fi = hap.get_frequency_info()
            pop = 'Unknown'
            for p in provs:
                if p in str(fi.by_population):
                    pop = p
                    break
            node_pop[n] = pop
        for n, pop in node_pop.items():
            lat, lon = provs.get(pop, provs['Unknown'])
            j_lat = lat + rng.uniform(-1.2, 1.2)
            j_lon = lon + rng.uniform(-1.2, 1.2)
            ax.scatter(j_lon, j_lat, s=40, c=loc_colors.get(pop, '#999999'),
                       alpha=0.7, edgecolors='white', linewidths=0.3, zorder=2)
        # 边 (单倍型连接)
        for e in network.edges:
            a, b = e[0], e[1]
            if a in node_pop and b in node_pop:
                la1, lo1 = provs.get(node_pop[a], provs['Unknown'])
                la2, lo2 = provs.get(node_pop[b], provs['Unknown'])
                ax.plot([lo1, lo2], [la1, la2], color='#bbbbbb', lw=0.4, alpha=0.5, zorder=0)
        ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
        ax.set_title(f'{virus} haplotype network over sampling provinces', fontsize=11)
        ax.set_ylim(20, 45); ax.set_xlim(95, 130)
        fig.tight_layout()
        fig.savefig(f"{cfg['out']}/{virus.lower()}_network_map.png", dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  → {virus.lower()}_network_map.png")
    except Exception as e:
        print(f"  地理地图失败: {e}")

    # 5. 网络统计
    print("网络统计...")
    stats = network.calculate_stats() if hasattr(network, 'calculate_stats') else None
    summary = {
        'nodes': network.num_nodes, 'edges': network.num_edges,
        'median_vectors': len(network.median_vector_ids) if hasattr(network, 'median_vector_ids') else None,
        'connected_components': len(network.get_connected_components()) if hasattr(network, 'get_connected_components') else None,
    }
    # 2026-09-15: 原写法把句柄丢给 json.dump, 永不 close (Windows 上延迟落盘/占用
    # 句柄), 且未指定编码。
    _stats_path = f"{cfg['out']}/network_stats.json"
    with open(_stats_path, 'w', encoding='utf-8') as _sf:
        json.dump(summary, _sf, indent=2)
    print(f"  → network_stats.json: {summary}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--virus', choices=known_dataset_names(), required=True)
    ap.add_argument('--color-by', choices=['location', 'host'], default='location',
                    help='着色方式: location(地点) 或 host(宿主)')
    args = ap.parse_args()
    build(args.virus.upper(), color_by=args.color_by)


if __name__ == '__main__':
    main()
