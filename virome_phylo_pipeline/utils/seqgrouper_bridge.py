#!/usr/bin/env python3
"""
seqgrouper_bridge.py — 序列分组（自研，等效 VirPhyKit SeqGrouper）
=================================================================
按元数据（location/host）对序列分组。不依赖 VirPhyKit。
"""

import os, csv, logging
from typing import Dict, List, Optional



class LogCollector:
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.messages = []
    def emit(self, msg):
        self.messages.append(msg); self.logger.info(msg)
    def warning(self, msg):
        self.messages.append(f"WARN: {msg}"); self.logger.warning(msg)  # 历史坑: 告警不入 messages → 下游读 log.messages 丢告警
    def error(self, msg):
        self.messages.append(f"ERROR: {msg}"); self.logger.error(msg)


def run_seqgrouper(
    metadata_csv: str,
    output_dir: str,
    group_by: str = "location",
    mapping_file: Optional[str] = None,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    桥接 VirPhyKit SeqGrouper 分组逻辑。

    复用:
    - function_group.py parse_single_record 的字段名 ('Geo Location', 'Host', 'Collection Date')
    - function_group.py preview_groups 的分组/映射表逻辑

    Parameters
    ----------
    metadata_csv : str
    output_dir : str
    group_by : str — "location" 或 "host"
    mapping_file : str, optional — 分组映射表
    log : LogCollector, optional

    Returns
    -------
    dict: {success, groups, pie_chart, summary_csv}
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)

    try:
        # 读取元数据
        with open(metadata_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        log.emit(f"SeqGrouper (VirPhyKit): {len(rows)} records")

        # 确定分组字段 (VirPhyKit 原版字段名)
        field_map = {
            'location': 'Geo Location',
            'host': 'Host',
            'date': 'Collection Date',
        }
        field = field_map.get(group_by, group_by)
        if field not in (rows[0].keys() if rows else []):
            # 尝试小写匹配
            alt = {k.lower(): k for k in (rows[0].keys() if rows else [])}
            if field.lower() in alt:
                field = alt[field.lower()]
            else:
                return {"success": False, "error": f"Column '{field}' not found",
                        "groups": {}, "pie_chart": None}

        # 加载映射表 (VirPhyKit preview_groups 逻辑)
        group_map = {}
        if mapping_file and os.path.exists(mapping_file):
            with open(mapping_file) as f:
                for line in f:
                    if line.strip() and '\t' in line:
                        group, val = line.strip().split('\t', 1)
                        group_map[val.lower().strip()] = group
            log.emit(f"  Mapping table: {len(group_map)} entries")

        # 分组 (VirPhyKit 原版: geo_loc_name 保留国家:省份层级)
        groups = {}
        for row in rows:
            val = row.get(field, '').strip()
            if not val or val in ('Not_Provided', 'N/A', '', 'Unknown'):
                groups.setdefault('Unknown', 0)
                groups['Unknown'] += 1
                continue

            # geo_loc_name 简化: 保留 "国家:省份" 或完整的 "国家:省份,城市"
            if group_by == 'location':
                val = val.replace('"', '').strip()
                # "China: Ningxia, Yinchuan" → "China: Ningxia"
                parts = val.split(':')
                if len(parts) >= 2:
                    province = parts[1].split(',')[0].strip()
                    val = f"{parts[0].strip()}:{province}"
                else:
                    val = parts[0].split(',')[0].strip()

            if group_by == 'host':
                # 简化宿主名: "Solanum tuberosum L." → "S. tuberosum"
                val = val.strip()

            # 应用映射表
            if group_map:
                label = group_map.get(val.lower(), val)
            else:
                label = val

            groups[label] = groups.get(label, 0) + 1

        groups = dict(sorted(groups.items(), key=lambda x: -x[1]))
        log.emit(f"  Groups: {len(groups)}")

        # 饼图 (VirPhyKit 风格)
        labels = list(groups.keys())
        sizes = list(groups.values())
        if len(labels) > 20:
            threshold = sum(sizes) * 0.02
            keep = [(l, s) for l, s in zip(labels, sizes) if s >= threshold]
            other_s = sum(s for _, s in zip(labels, sizes) if s < threshold)
            labels, sizes = zip(*keep) if keep else ([], [])
            labels, sizes = list(labels), list(sizes)
            if other_s > 0:
                labels.append('Other'); sizes.append(other_s)

        fig, ax = plt.subplots(figsize=(10, 8))
        colors = plt.get_cmap('tab20', len(labels))(range(len(labels)))
        wedges, texts, autotexts = ax.pie(
            sizes, labels=None, autopct='%1.1f%%', startangle=90,
            colors=colors, pctdistance=0.85)
        ax.legend(wedges, [f'{l} ({s})' for l, s in zip(labels, sizes)],
                 title=group_by.title(), loc='center left',
                 bbox_to_anchor=(1, 0.5), fontsize=9)
        ax.set_title(f'Sequence Distribution by {group_by.title()}',
                    fontsize=14, fontweight='bold')

        pie_path = os.path.join(output_dir, f"seqgroup_{group_by}.pdf")
        fig.savefig(pie_path, bbox_inches='tight', dpi=150, format='pdf')
        plt.close(fig)

        summary_path = os.path.join(output_dir, f"seqgroup_{group_by}.csv")
        with open(summary_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow([group_by.title(), 'Count', 'Percentage'])
            total = sum(sizes)
            for l, s in zip(labels, sizes):
                w.writerow([l, s, f'{s*100/total:.1f}%'])

        return {
            "success": True,
            "groups": groups,
            "pie_chart": pie_path,
            "summary_csv": summary_path,
        }

    except Exception as e:
        log.error(f"SeqGrouper failed: {e}")
        return {"success": False, "error": str(e), "groups": {}, "pie_chart": None}
