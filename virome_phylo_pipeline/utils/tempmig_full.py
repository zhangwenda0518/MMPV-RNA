#!/usr/bin/env python3
"""tempmig_full.py — VirPhyKit MOT (Migration Over Time) 完整复刻

复刻 VirPhyKit 的 MOT + Get_categories 核心算法 (Ola Brynildsrud global L4 scripts):
  1. 解析 BEAST MCC 树的节点注释 (location={X=prob}, height=时间)
  2. 沿树遍历, 检测每个分支的"父地点 → 子地点"变化 (迁移事件)
  3. 用节点 height + branch length 把迁移事件映射到日历年份区间
  4. 输出迁移矩阵 (地点对 × 年份) + 时序迁移图

与 geo_analysis.run_tempmig 的区别:
  run_tempmig      简化版: 时间分桶 + 柱状图 (相对时间)
  run_tempmig_full 完整版: 迁移矩阵 × 日历年份 (VirPhyKit MOT 原版)

用法:
  python -m utils.tempmig_full --mcc mcc_annotated.tre --dates phylo_dates.csv \
      --outdir out [--trait location]
"""
import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict

from Bio import Phylo


def parse_beast_annotation(clade):
    """从 BEAST 1.x clade 注释提取 location + height。支持两种注释格式:
    1. location={Ningxia=0.394, Beijing=0.289}   (BEAUti 后处理/BEAST2)
    2. Location.set={"Ningxia",...} + Location.set.prob={0.9,...} (treeannotator 原生)
    3. Location="Ningxia"                         (确定性状态)
    """
    info = {}
    if not hasattr(clade, 'comment') or not clade.comment:
        return info
    # 只剥方括号与首个 '&'：BEAST 注释形如
    # [&height=…,type.set={A,B},type.set.prob={0.6,0.4}]，strip('{}') 会把
    # **最后一个键的收尾 '}' 也剥掉** → 多状态节点全部解析失败
    # （实测 MTT_MCC.tre 341 节点只成功 185，成功的都是 prob=1.0 的单状态节点）。
    comment = clade.comment.strip()
    if comment.startswith('['):
        comment = comment[1:]
    if comment.endswith(']'):
        comment = comment[:-1]
    if comment.startswith('&'):
        comment = comment[1:]
    # 格式1: location={Ningxia=0.394, Beijing=0.289}
    m = re.search(r'location=\{([^}]*)\}', comment)
    if m:
        locs = {}
        for item in m.group(1).split(','):
            if '=' in item:
                k, v = item.split('=', 1)
                try:
                    locs[k.strip()] = float(v.strip())
                except ValueError:
                    pass
        if locs:
            info['location'] = locs
    # 格式2: Location.set={"X","Y"} + Location.set.prob={0.9,0.1}
    if 'location' not in info:
        m_set = re.search(r'Location\.set=\{([^}]*)\}', comment)
        m_prob = re.search(r'Location\.set\.prob=\{([^}]*)\}', comment)
        if m_set:
            names = [x.strip().strip('"') for x in m_set.group(1).split(',') if x.strip()]
            if m_prob:
                probs = []
                for x in m_prob.group(1).split(','):
                    try:
                        probs.append(float(x.strip()))
                    except ValueError:
                        probs.append(0.0)
                if len(probs) == len(names):
                    info['location'] = dict(zip(names, probs))
            elif len(names) == 1:
                info['location'] = {names[0]: 1.0}
    # 格式2b: 任意 `<X>.set` / `<X>.set.prob`（VirPhyKit function_rrt 的通用做法）
    #   · 示例里的性状名各不相同（max / Region / type / location），写死必然漏
    #   · 也处理 `max.set={...}` 这种「最可能集合」
    if 'location' not in info:
        for _m in re.finditer(r'(?:^|[,&])\s*([A-Za-z][\w%.\-]*)\.set=\{([^}]*)\}',
                              comment):
            _name = _m.group(1)
            _names = [x.strip().strip('"') for x in _m.group(2).split(',') if x.strip()]
            _pm = re.search(r'(?:^|[,&])\s*' + re.escape(_name)
                            + r'\.set\.prob=\{([^}]*)\}', comment)
            if _pm:
                _probs = []
                for _x in _pm.group(1).split(','):
                    try:
                        _probs.append(float(_x.strip()))
                    except ValueError:
                        _probs.append(0.0)
                if len(_probs) == len(_names):
                    info['location'] = dict(zip(_names, _probs))
                    break
            if len(_names) == 1:
                info['location'] = {_names[0]: 1.0}
                break
    # 格式3: Location="Ningxia" (确定性状态)
    if 'location' not in info:
        m_loc = re.search(r'Location="([^"]+)"', comment)
        if m_loc:
            info['location'] = {m_loc.group(1): 1.0}
    # height=0.207
    m = re.search(r'height=([\d.eE-]+)', comment)
    if m:
        info['height'] = float(m.group(1))
    return info


def tip_year(name, dates_map):
    """tip 采样年份 (支持多种命名格式)"""
    if name in dates_map:
        return dates_map[name]
    # ⚠️ 必须在 `(\d{4})$` 之前：BEAST 常把十进制年写在名字末段
    # （CCD_AB622861_2008.58197），`(\d{4})$` 会抓到小数尾巴 `8197`，
    # 于是「最新采样年」变成 9863、整条年份轴报废。
    for _tok in reversed((name or '').split('_')):
        if re.fullmatch(r'\d{4}(\.\d+)?', _tok):
            return int(float(_tok))
    # 提取 accession: SRR12805583_OR489165.1 → SRR12805583
    m = re.search(r'^(CRR\d+|SRR\d+|ERR\d+|DRR\d+)', name)
    if m and m.group(1) in dates_map:
        return dates_map[m.group(1)]
    # VirNA 格式 /2021-06-01
    m = re.search(r'/(\d{4})[-/]', name)
    if m:
        return int(m.group(1))
    m = re.search(r'(\d{4})$', name)
    if m:
        return int(m.group(1))
    return None


def run_tempmig_full(mcc_tree_file, dates_csv, output_dir, trait='location'):
    """VirPhyKit MOT 完整复刻: 迁移矩阵 × 日历年份"""
    os.makedirs(output_dir, exist_ok=True)
    result = {'success': False, 'error': None}

    try:
        # 1. 采样日期映射
        dates_map = {}
        if dates_csv and os.path.exists(dates_csv):
            with open(dates_csv, encoding='utf-8-sig') as f:
                for row in csv.DictReader(f):
                    name = row.get('name', '')
                    date = row.get('date', '') or row.get('CollectionDate', '')
                    if name and date:
                        m = re.search(r'(\d{4})', date)
                        if m:
                            dates_map[name] = int(m.group(1))

        # 2. 解析树 (自动检测 NEXUS / Newick)
        # 教训: 旧版写死 'nexus', 喂 .newick 文件会解析出空树且静默返回
        # (exit 0 但零产物)。现在按文件头选择解析器, 并在错误信息中报告格式。
        try:
            with open(mcc_tree_file, encoding='utf-8', errors='replace') as _f:
                _head = _f.read(4096)
        except OSError as e:
            result['error'] = f'无法读取 MCC 文件: {e}'
            return result
        tree_fmt = 'nexus' if _head.lstrip().upper().startswith('#NEXUS') else 'newick'
        trees = list(Phylo.parse(mcc_tree_file, tree_fmt))
        if not trees:
            result['error'] = f'无树 (文件格式: {tree_fmt})'
            return result
        tree = trees[0]

        # 3. 提取节点信息 (location 从注释 + height 累积枝长)
        node_info = {}
        root = tree.root

        def process(clade, parent_clade=None, cum_height=0.0):
            info = parse_beast_annotation(clade)
            loc = None
            if info.get('location'):
                loc = max(info['location'], key=info['location'].get)
            # height = BEAST 节点年龄 (root=树高, tip=0)。VirPhyKit Get_categories 用
            # treeio read.beast 的节点 height (年龄)。旧版误用从根累积枝长
            # (root=0, tip=树高)，与原版互为补，导致整张年份轴镜像倒置。
            # 优先用 BEAST 注释里的 height=，无注释时用 root_height - 累积枝长换算。
            height = info.get('height')
            if height is None:
                height = None  # 树高需先归一化，在 process 后补算 (见下方 fixup)
            node_info[id(clade)] = {'clade': clade, 'height': height,
                                    'cum_height': cum_height,
                                    'location': loc, 'parent': parent_clade}
            for c in clade.clades:
                process(c, clade, cum_height + (c.branch_length or 0.0))

        process(root)

        # 树高 = 最大累积枝长 (tip 最深)
        root_height = max(v['cum_height'] for v in node_info.values())
        # 无 BEAST height 注释的节点: age = root_height - 累积枝长 (root=树高, tip=0)
        for v in node_info.values():
            if v['height'] is None:
                v['height'] = root_height - v['cum_height']

        # 完整性检查: BEAST 注释在 Newick 里可能被解析器丢失 (location 覆盖率骤降),
        # 此时迁移事件会被严重低估。覆盖率 < 60% 时警告, 建议改用 NEXUS 格式 MCC。
        n_loc = sum(1 for v in node_info.values() if v['location'])
        if n_loc < len(node_info) * 0.6:
            print(f'[tempmig_full] 警告: 仅 {n_loc}/{len(node_info)} 节点解析出 location '
                  f'(格式: {tree_fmt})。Newick 输入可能丢失 BEAST 注释, 请改用 '
                  f'NEXUS 格式 MCC 树 (treeannotator 输出) 以保证迁移事件完整。',
                  file=sys.stderr)

        # 4. 采样年份
        tip_years = []
        for c in tree.get_terminals():
            y = tip_year(c.name, dates_map)
            if y:
                tip_years.append(y)
        if not tip_years:
            result['error'] = '无法从 tip 名或 dates_csv 提取采样年份'
            return result
        most_current_year = max(tip_years)
        root_year = most_current_year - int(math.ceil(root_height))
        print(f'树高: {root_height:.2f}, 最新采样年: {most_current_year}, 根年: {root_year}')

        # 5. 迁移矩阵 (VirPhyKit Get_categories 复刻)
        unique_locs = sorted({v['location'] for v in node_info.values()
                              if v['location']})
        row_keys = [f'{a}_to_{b}' for a in unique_locs for b in unique_locs]
        col_keys = list(range(root_year, most_current_year + 1))
        calendar = {r: {c: 0 for c in col_keys} for r in row_keys}

        def find_placement(height, length):
            end = int(math.ceil(root_height) - math.ceil(height))
            start = int(end - math.floor(length))
            return start, end

        # 6. 逐节点检测迁移 (VirPhyKit 核心)
        for cid, info in node_info.items():
            clade = info['clade']
            if clade is root:
                continue
            node_height = info['height']
            node_loc = info['location']
            if node_height is None or node_loc is None:
                continue
            node_length = clade.branch_length or 0.0
            parent_clade = info['parent']
            if parent_clade is None:
                continue
            parent_loc = node_info[id(parent_clade)]['location']
            if parent_loc is None or parent_loc == node_loc:
                continue  # 无迁移

            # 迁移: 父地点 → 子地点, 映射年份区间
            start, end = find_placement(node_height, node_length)
            start += root_year
            end += root_year
            row = f'{parent_loc}_to_{node_loc}'
            for z in range(start, end + 1):
                try:
                    calendar[row][z] += 1
                except KeyError:
                    pass

        # 7. 输出迁移矩阵 CSV
        mat_path = os.path.join(output_dir, 'migration_over_time.csv')
        with open(mat_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['route'] + col_keys)
            for r in row_keys:
                w.writerow([r] + [calendar[r][c] for c in col_keys])

        # 8. 汇总: 每年总迁移事件
        yearly_total = {c: sum(calendar[r][c] for r in row_keys) for c in col_keys}
        summary = [{'year': c, 'n_migrations': yearly_total[c]} for c in col_keys]
        sum_path = os.path.join(output_dir, 'migration_yearly_summary.csv')
        with open(sum_path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['year', 'n_migrations'])
            w.writeheader()
            w.writerows(summary)

        # 9. 图 (时序迁移)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        years = [s['year'] for s in summary]
        counts = [s['n_migrations'] for s in summary]
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.bar(years, counts, color='#5BD4D8', edgecolor='white')
        ax.set_xlabel('Year', fontsize=13)
        ax.set_ylabel('Migration events', fontsize=13)
        ax.set_title('Migration Over Time (VirPhyKit MOT)', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plot_path = os.path.join(output_dir, 'migration_over_time.pdf')
        fig.savefig(plot_path, dpi=300, format='pdf', bbox_inches='tight')
        plt.close(fig)

        # 10. 热门迁移路线
        top_routes = []
        for r in row_keys:
            total = sum(calendar[r][c] for c in col_keys)
            if total > 0:
                top_routes.append((r, total))
        top_routes.sort(key=lambda x: -x[1])

        result.update({
            'success': True,
            'root_height': root_height, 'root_year': root_year,
            'most_current_year': most_current_year,
            'n_locations': len(unique_locs),
            'matrix': mat_path, 'summary': sum_path, 'plot': plot_path,
            'top_routes': top_routes[:10],
            'total_migrations': sum(yearly_total.values()),
        })
        print(f'\n总迁移事件: {sum(yearly_total.values())}')
        print(f'地点数: {len(unique_locs)}: {unique_locs}')
        print(f'热门路线 (top 5):')
        for r, t in top_routes[:5]:
            print(f'  {r}: {t} 次')

        return result

    except Exception as e:
        import traceback
        result['error'] = str(e)
        traceback.print_exc()
        return result


def main():
    ap = argparse.ArgumentParser(description='VirPhyKit MOT 迁移时序复刻')
    ap.add_argument('--mcc', required=True, help='BEAST MCC 注释树')
    ap.add_argument('--dates', default=None, help='采样日期 CSV (name,date)')
    ap.add_argument('--outdir', default='tempmig_full_out')
    args = ap.parse_args()
    result = run_tempmig_full(args.mcc, args.dates, args.outdir)
    if not result.get('success'):
        err = result.get('error') or '未知错误'
        print(f'[tempmig_full] 失败: {err}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
