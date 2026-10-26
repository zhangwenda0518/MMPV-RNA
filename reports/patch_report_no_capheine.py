# -*- coding: utf-8 -*-
"""generate_pipeline_report.py：对齐当前 9 阶段已知病毒管道（无 capheine）。

背景
----
auto_known_virus.py 第 19 行注明：
    "positive selection (capheine) analysis migrated to virome_phylo_pipeline/"
且 --stage 可选值已无 capheine（detect/filter/variants/full/extract/post/
similarity/dvg/report）。但报告生成器仍按 10 阶段硬编码：
  * report_dir = out / "10_report"          -> 新产物目录是 09_report
  * sim_dir    = out / "08_similarity"      -> 新产物目录是 07_similarity
  * S7 Capheine 导航段 / 流程节点无条件输出  -> 无 capheine 时是空段

改法（向后兼容：存在旧目录时仍走旧路径，v2 之类 10 阶段产物不受影响）
----
1. generate_html() 内计算 HAS_CAPHEINE = (out/"07_capheine").is_dir()
2. 无 capheine 时跳过 S7（用 sid 判断 + continue，S8/S9 的 sn 不受影响）
3. 流程节点的 S7 节点改为条件输出
4. copy_results_to_report(): sim_dir 布局优先，回退 08_similarity
5. main(): report_dir 已有 10_report 则沿用，否则用当前布局名 _D['a_report']
"""
import io
import os

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_REPO   = f"{_SERVER}/MMPV-RNA"


P = f"{_REPO}/virome_analysis_pipeline/generate_pipeline_report.py"
src = io.open(P, encoding='utf-8').read()

# ---------------------------------------------------------------- 1) HAS_CAPHEINE
OLD1 = '''def generate_html(out_dir, out_html, viruses, ai_html="", stage_summaries=None):
    """Generate interactive HTML report with sidebar navigation."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out = Path(out_dir)

    def _d(seg):
'''
NEW1 = '''def generate_html(out_dir, out_html, viruses, ai_html="", stage_summaries=None):
    """Generate interactive HTML report with sidebar navigation."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out = Path(out_dir)
    # 当前已知病毒管道已不含 capheine(迁移至 virome_phylo_pipeline)：
    # 无 07_capheine 产物时不输出 S7 段(导航 + 流程节点)，其余阶段编号保持不变。
    HAS_CAPHEINE = (out / "07_capheine").is_dir()

    def _d(seg):
'''
c1 = src.count(OLD1)
assert c1 == 1, 'OLD1 count = %d (expect 1)' % c1
src = src.replace(OLD1, NEW1, 1)
print('OK[1]: HAS_CAPHEINE 已定义')

# ---------------------------------------------------------------- 2) 跳过 S7
OLD2 = '''    for i, (sid, sname, sdesc, has_sub) in enumerate(stage_cfg):
'''
NEW2 = '''    for i, (sid, sname, sdesc, has_sub) in enumerate(stage_cfg):
        if sid == "S7" and not HAS_CAPHEINE:
            continue          # 无 capheine 产物 -> 跳过 S7；S8/S9 的 sn 不变
'''
c2 = src.count(OLD2)
assert c2 == 1, 'OLD2 count = %d (expect 1)' % c2
src = src.replace(OLD2, NEW2, 1)
print('OK[2]: 无 capheine 时跳过 S7')

# ---------------------------------------------------------------- 3) 流程节点
OLD3 = '''        <div class="flow-node"><b>S7</b><small>Selection</small></div><span class="flow-arrow">&#8594;</span>
'''
NEW3 = '''        {('<div class="flow-node"><b>S7</b><small>Selection</small></div><span class="flow-arrow">&#8594;</span>' if HAS_CAPHEINE else '')}
'''
c3 = src.count(OLD3)
assert c3 == 1, 'OLD3 count = %d (expect 1)' % c3
src = src.replace(OLD3, NEW3, 1)
print('OK[3]: 流程节点 S7 改为条件输出')

# ---------------------------------------------------------------- 4) sim_dir
OLD4 = '''    sim_dir = out / "08_similarity"
'''
NEW4 = '''    sim_dir = out / _D['a_similarity']
    if not sim_dir.is_dir():
        sim_dir = out / "08_similarity"      # 兼容 10 阶段旧产物
'''
c4 = src.count(OLD4)
assert c4 == 1, 'OLD4 count = %d (expect 1)' % c4
src = src.replace(OLD4, NEW4, 1)
print('OK[4]: sim_dir 改为布局感知')

# ---------------------------------------------------------------- 5) report_dir
OLD5 = '''    report_dir = out / "10_report"
'''
NEW5 = '''    _legacy_rep = out / "10_report"
    report_dir = _legacy_rep if _legacy_rep.is_dir() else out / _D['a_report']
'''
c5 = src.count(OLD5)
assert c5 == 1, 'OLD5 count = %d (expect 1)' % c5
src = src.replace(OLD5, NEW5, 1)
print('OK[5]: report_dir 改为布局感知（兼容旧 10_report）')

io.open(P, 'w', encoding='utf-8').write(src)
print('DONE')
