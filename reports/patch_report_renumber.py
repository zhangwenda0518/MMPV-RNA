# -*- coding: utf-8 -*-
"""generate_pipeline_report.py 补丁 2：无 capheine 时重编号显示（S8->S7, S9->S8）。

补丁 1 已让 S7 段不输出，但导航/流程/目录名仍是 S8/S9，编号从 S6 跳到 S8。
本补丁在有/无 capheine 两种情况下都能给出连续编号：
  * 有 07_capheine（v2 等 10 阶段旧产物）：S7 Capheine / S8 Similarity / S9 DVG（不变）
  * 无 07_capheine（当前 9 阶段管道）：S7 Similarity / S8 DVG

注意：内容分支用的是 enumerate 得到的 sn（未受影响），只改显示用的 sid
和 copy_results_to_report 里的目录名。
"""
import io
import os

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_REPO   = f"{_SERVER}/MMPV-RNA"


P = f"{_REPO}/virome_analysis_pipeline/generate_pipeline_report.py"
src = io.open(P, encoding='utf-8').read()

# ---------------------------------------------- 6) stage_cfg 显示编号重排
OLD6 = '''        ("S9","DVG","ViReMa DVG detection: Circos 4-track recombination plots + arc diagrams.",True),
    ]
'''
NEW6 = '''        ("S9","DVG","ViReMa DVG detection: Circos 4-track recombination plots + arc diagrams.",True),
    ]
    if not HAS_CAPHEINE:
        # 无 capheine：显示编号前移一位，避免出现 S6 -> S8 的跳号
        # (内容分支用 enumerate 的 sn，不受影响)
        stage_cfg = [
            ("S7", n, d, h) if s == "S8" else ("S8", n, d, h) if s == "S9" else (s, n, d, h)
            for (s, n, d, h) in stage_cfg
        ]
'''
c6 = src.count(OLD6)
assert c6 == 1, 'OLD6 count = %d (expect 1)' % c6
src = src.replace(OLD6, NEW6, 1)
print('OK[6]: stage_cfg 显示编号条件重排')

# ---------------------------------------------- 7a) copy 段定义目录名
OLD7 = '''    out = Path(out_dir); rep = Path(report_dir)
    rep.mkdir(parents=True, exist_ok=True)
    copied = []; skipped = []
'''
NEW7 = '''    out = Path(out_dir); rep = Path(report_dir)
    rep.mkdir(parents=True, exist_ok=True)
    # 有 07_capheine(10 阶段旧产物) 时用 S8/S9, 否则前移为 S7/S8
    _has_cap = (out / "07_capheine").is_dir()
    _sim_seg = "S8_Similarity" if _has_cap else "S7_Similarity"
    _dvg_seg = "S9_DVG" if _has_cap else "S8_DVG"
    copied = []; skipped = []
'''
c7 = src.count(OLD7)
assert c7 == 1, 'OLD7 count = %d (expect 1)' % c7
src = src.replace(OLD7, NEW7, 1)
print('OK[7a]: 复制目录名改为条件')

# ---------------------------------------------- 7b) similarity 目标目录
OLD8 = '''                for f in vdir.rglob("*heatmap*"):
                    _cp(str(f.relative_to(out)), f"S8_Similarity/{vdir.name}/{f.name}")
                for f in vdir.rglob("*composite*"):
                    _cp(str(f.relative_to(out)), f"S8_Similarity/{vdir.name}/{f.name}")
                for f in vdir.rglob("*similarity*"):
                    _cp(str(f.relative_to(out)), f"S8_Similarity/{vdir.name}/{f.name}")
'''
NEW8 = '''                for f in vdir.rglob("*heatmap*"):
                    _cp(str(f.relative_to(out)), f"{_sim_seg}/{vdir.name}/{f.name}")
                for f in vdir.rglob("*composite*"):
                    _cp(str(f.relative_to(out)), f"{_sim_seg}/{vdir.name}/{f.name}")
                for f in vdir.rglob("*similarity*"):
                    _cp(str(f.relative_to(out)), f"{_sim_seg}/{vdir.name}/{f.name}")
'''
c8 = src.count(OLD8)
assert c8 == 1, 'OLD8 count = %d (expect 1)' % c8
src = src.replace(OLD8, NEW8, 1)
print('OK[7b]: similarity 目标目录名条件化')

# ---------------------------------------------- 7c) DVG 目标目录
OLD9 = '''                _cp(str(vdir.relative_to(out)), f"S9_DVG/{vdir.name}")
'''
NEW9 = '''                _cp(str(vdir.relative_to(out)), f"{_dvg_seg}/{vdir.name}")
'''
c9 = src.count(OLD9)
assert c9 == 1, 'OLD9 count = %d (expect 1)' % c9
src = src.replace(OLD9, NEW9, 1)
print('OK[7c]: DVG 目标目录名条件化')

io.open(P, 'w', encoding='utf-8').write(src)
print('DONE')
