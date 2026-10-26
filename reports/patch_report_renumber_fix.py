# -*- coding: utf-8 -*-
"""generate_pipeline_report.py 补丁 3：修正重编号的两处遗漏。

补丁 2 的 bug：remap 把 Similarity 的 sid 改成 "S7" 后，循环里
`if sid == "S7" and not HAS_CAPHEINE: continue` 会把 Similarity 也一起跳过。
改为按下标判断（capheine 是第 7 项，下标 6），与 sid 无关。

另外流程条(flow-strip)是静态 HTML，S8/S9 标签需同样条件化。
"""
import io
import os

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_REPO   = f"{_SERVER}/MMPV-RNA"


P = f"{_REPO}/virome_analysis_pipeline/generate_pipeline_report.py"
src = io.open(P, encoding='utf-8').read()

# ---------------------------------------------- 8) 跳过条件改按下标
OLD1 = '''        if sid == "S7" and not HAS_CAPHEINE:
            continue          # 无 capheine 产物 -> 跳过 S7；S8/S9 的 sn 不变
'''
NEW1 = '''        if i == 6 and not HAS_CAPHEINE:
            # 无 capheine 产物 -> 跳过 Capheine 段(第 7 项, 下标 6)。
            # 按下标而非 sid 判断: 重编号后 Similarity 的 sid 也变成 S7。
            continue
'''
c1 = src.count(OLD1)
assert c1 == 1, 'OLD1 count = %d (expect 1)' % c1
src = src.replace(OLD1, NEW1, 1)
print('OK[8]: 跳过条件改为按下标 i == 6')

# ---------------------------------------------- 9) 流程条标签条件化
OLD2 = '''        <div class="flow-node"><b>S8</b><small>Similarity</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S9</b><small>DVG</small></div>
'''
NEW2 = '''        <div class="flow-node"><b>{'S8' if HAS_CAPHEINE else 'S7'}</b><small>Similarity</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>{'S9' if HAS_CAPHEINE else 'S8'}</b><small>DVG</small></div>
'''
c2 = src.count(OLD2)
assert c2 == 1, 'OLD2 count = %d (expect 1)' % c2
src = src.replace(OLD2, NEW2, 1)
print('OK[9]: 流程条 S8/S9 标签条件化')

io.open(P, 'w', encoding='utf-8').write(src)
print('DONE')
