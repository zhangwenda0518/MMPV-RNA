#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
修复 rescue_report 统计口径 bug。

问题: "过短 (<2000bp) 12265 条, 占比 110.2%"
根因: failed = total - rescued 用了两个不同口径的数
  - total = len(centroids_records)          = 20517  (原始 centroids 条数)
  - rescued = cnt_a+cnt_b+cnt_c+cnt_d       = 9392   (各分支产物文件条数, 含组组装新contig)
  - failed  = 20517 - 9392 = 11125
  但 fail_reasons 的累加发生在遍历原始 centroids 的循环里, 实际 total = 20467
  故 12265/11125 = 110.2%

修复: failed 改用 fail_reasons 累加和 (与分子同口径);
      占比分母统一为 total = len(centroids_records), 语义清晰;
      分支表新增"未拯救"行用同一分母。
"""
import sys, py_compile

TARGET = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/rescue_pipeline.py"

with open(TARGET, 'r', encoding='utf-8') as f:
    src = f.read()
orig = src

# ---------- 1) failed 的定义 ----------
old1 = '''    rescued = cnt_a + cnt_b + cnt_c + cnt_d
    total = len(centroids_records)
    failed = total - rescued
'''
new1 = '''    rescued = cnt_a + cnt_b + cnt_c + cnt_d
    total = len(centroids_records)
    # 注意: rescued 数的是各分支产物文件里的序列 (分支C的组组装会引入新 contig),
    # 与 total (原始 centroids 条数) 口径不同, 不可直接相减。
    # 未拯救数应取 fail_reasons 的累加和 —— 它在遍历原始 centroids 时逐条累加,
    # 与失败原因占比的分子同口径。
    failed = sum(fail_reasons.values())
    n_centroids_fail = failed  # 语义别名: 原始终端中未拯救的条数
'''
assert src.count(old1) == 1, "failed 定义锚点不唯一"
src = src.replace(old1, new1)

# ---------- 2) 分支统计表: 补分母说明 + 修正未拯救占比 ----------
old2 = '''        mf.write(f"| **合计** | | **{rescued}** | **{rescued/total*100:.1f}%** |\\n")
        mf.write(f"| 未拯救 | | {failed} | {failed/total*100:.1f}% |\\n\\n")'''
new2 = '''        mf.write(f"| **分支合计** | | **{rescued}** | — |\\n")
        mf.write(f"| 原始 centroids 中未拯救 | | {failed} | {failed/total*100:.1f}% |\\n\\n")
        mf.write(f"> 说明: 分支拯救数为各分支产物序列数(分支 C 的组组装会产出新 contig), "
                 f"与原始 centroids 总数 {total} 非同口径, 故分支合计不给出占比。"
                 f"未拯救占比以原始 centroids 为分母。\\n\\n")'''
assert src.count(old2) == 1, "合计行锚点不唯一"
src = src.replace(old2, new2)

# ---------- 3) 失败原因表: 分母明确 ----------
old3 = '''        if fail_reasons:
            mf.write(f"## 未拯救原因分布\\n\\n")
            mf.write(f"| 原因 | 数量 | 占比 |\\n")
            mf.write(f"|------|------|------|\\n")
            for reason, count in sorted(fail_reasons.items(), key=lambda x: -x[1]):
                label = reason_labels.get(reason, reason)
                mf.write(f"| {label} | {count} | {count/failed*100:.1f}% |\\n")
            mf.write("\\n")'''
new3 = '''        if fail_reasons:
            mf.write(f"## 未拯救原因分布\\n\\n")
            mf.write(f"| 原因 | 数量 | 占未拯救 |\\n")
            mf.write(f"|------|------|---------|\\n")
            _reasons_sum = sum(fail_reasons.values())
            for reason, count in sorted(fail_reasons.items(), key=lambda x: -x[1]):
                label = reason_labels.get(reason, reason)
                _pct = (count / _reasons_sum * 100) if _reasons_sum else 0.0
                mf.write(f"| {label} | {count} | {_pct:.1f}% |\\n")
            mf.write(f"| **合计** | **{_reasons_sum}** | **100.0%** |\\n")
            mf.write("\\n")'''
assert src.count(old3) == 1, "失败原因表锚点不唯一"
src = src.replace(old3, new3)

if src == orig:
    print("NO CHANGE APPLIED")
    sys.exit(1)

with open(TARGET, 'w', encoding='utf-8') as f:
    f.write(src)
try:
    py_compile.compile(TARGET, doraise=True)
    print("PATCH OK, syntax OK")
except py_compile.PyCompileError as e:
    print("SYNTAX ERROR!", e)
    sys.exit(1)
