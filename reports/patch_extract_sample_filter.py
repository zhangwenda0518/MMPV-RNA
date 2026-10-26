# -*- coding: utf-8 -*-
"""extract_barbarum_complete.py：c_default 改为样本级过滤。

背景
----
c_default() 对 04/07/08/09/10 一律"整树复制"，导致非 barbarum 样本被整棵
带入镜像。实测（v2 源，按 SRA run 号判定）：
    09_dvg/virema_results        790 文件 + 84 目录
    09_dvg/all_gathered_beds     711 文件
    03_variants/summary/individual_virus_reports   96 文件
    10_report/S6_Post_hoc          2 文件
合计约 1,943（与审计发现的混入文件数吻合）。

方案
----
新增通用递归收集器：
  * 名字形如 SRA run（SRR/CRR/ERR/DRR + 数字）的条目 -> 必须在 RUNS 内，
    否则连同子树一并丢弃；
  * 不形如 run 的条目（病毒级 / 全局 / 汇总）-> 一律保留；
  * 目录下若存在形如 run 的后代则继续细分，否则整树复制（保持效率）。

注意：run 号正则必须是 SRR|CRR|ERR|DRR 前缀。早期审计用 [A-Za-z]{3}\\d+
会把 "Fig01_VAF_Spectrum.pdf" 误判成样本名。
"""
import io
import os

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


P = f"{_GOJI}/00_scripts/extract_barbarum_complete.py"
src = io.open(P, encoding='utf-8').read()

# ---------------------------------------------------------------- 1) 插入通用工具
OLD1 = '''def _table(src, rel, col="Sample"):
    return dict(src=src, rel=rel, kind="table_tsv", col=col)


def c_01(src, runs):
'''

NEW1 = '''def _table(src, rel, col="Sample"):
    return dict(src=src, rel=rel, kind="table_tsv", col=col)


# ------------------------------------------- 样本级过滤（通用，供 04/07-10 用）
# c_default 原先对 04/07/09/10 整树复制，把非 barbarum 样本整棵带进镜像
# (09_dvg 1501 + individual_virus_reports 96 + 10_report/S6_Post_hoc 2)。
# 规则: 名字形如 SRA run 的条目必须属于 RUNS，否则连同子树丢弃；
#       不形如 run 的（病毒级 / 全局 / 汇总）一律保留。
# 注意: 必须是 SRR|CRR|ERR|DRR 前缀，用 [A-Za-z]{3}\\d+ 会把
#       "Fig01_VAF_Spectrum.pdf" 之类误判成样本名。
_RUN_SHAPE_RE = re.compile(r"^(?:sample_)?((?:SRR|CRR|ERR|DRR)\\d+)(?=[._-]|$)")


def _run_shape(name):
    """名字是否形如 SRA run；是则返回 run 号，否则 None。"""
    m = _RUN_SHAPE_RE.match(name)
    return m.group(1) if m else None


def _keep_sample(name):
    """非样本名 -> 保留；样本名 -> 仅当属于 barbarum 才保留。"""
    r = _run_shape(name)
    return True if r is None else (r in RUNS)


def _has_run_descendant(d, depth=0):
    """子树中是否存在形如 SRA run 的条目（决定该目录是否需要继续细分）。"""
    if depth > 8:
        return False
    try:
        names = sorted(os.listdir(d))
    except OSError:
        return False
    if any(_run_shape(x) for x in names):
        return True
    for x in names:
        p = os.path.join(d, x)
        if os.path.isdir(p) and _has_run_descendant(p, depth + 1):
            return True
    return False


def _walk_samples(src, rel_root, E):
    """递归收集：样本级条目按 barbarum 过滤，其余整树 / 单文件保留。"""
    for x in sorted(os.listdir(src)):
        p = os.path.join(src, x)
        r = "%s/%s" % (rel_root, x) if rel_root else x
        if not _keep_sample(x):
            continue
        if os.path.isdir(p):
            if _has_run_descendant(p):
                _walk_samples(p, r, E)
            else:
                E.append(_tree(p, r))
        else:
            E.append(_file(p, r))


def c_01(src, runs):
'''

c1 = src.count(OLD1)
assert c1 == 1, 'OLD1 (helpers anchor) count = %d (expect 1)' % c1
src = src.replace(OLD1, NEW1, 1)
print('OK[1]: 通用样本级过滤工具已插入')

# ---------------------------------------------------------------- 2) 重写 c_default
OLD2 = '''def c_default(src, runs):
    """07/08/09/10：病毒级/全局，整体复制"""
    E = []
    for x in sorted(os.listdir(src)):
        p = os.path.join(src, x)
        E.append(_tree(p, x) if os.path.isdir(p) else _file(p, x))
    return E
'''

NEW2 = '''def c_default(src, runs):
    """04/07/08/09/10：病毒级/全局整体保留，样本级条目按 barbarum 过滤"""
    E = []
    _walk_samples(src, "", E)
    return E
'''

c2 = src.count(OLD2)
assert c2 == 1, 'OLD2 (c_default) count = %d (expect 1)' % c2
src = src.replace(OLD2, NEW2, 1)
print('OK[2]: c_default 已改为样本级过滤')

# ---------------------------------------------------------------- 3) summary 子目录
# 只改 c_03 那处：c_01 的 summary 行带注释 "# 表按 Sample 过滤"，锚点不同，
# 且实测 01_detection/summary 里 0 个 run 形条目（全是汇总/图），无需过滤。
# 03_variants/summary/individual_virus_reports 则 442/442 全是样本文件。
OLD3 = '''    d = os.path.join(src, "summary")
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            fp = os.path.join(d, fn)
            if os.path.isdir(fp):
                E.append(_tree(fp, "summary/%s" % fn))
            elif fn.endswith(".tsv"):
                E.append(_table(fp, "summary/%s" % fn, "Sample"))
            else:
                E.append(_file(fp, "summary/%s" % fn))
'''

NEW3 = '''    d = os.path.join(src, "summary")
    if os.path.isdir(d):
        for fn in sorted(os.listdir(d)):
            fp = os.path.join(d, fn)
            if os.path.isdir(fp):
                if _has_run_descendant(fp):
                    _walk_samples(fp, "summary/%s" % fn, E)
                else:
                    E.append(_tree(fp, "summary/%s" % fn))
            elif fn.endswith(".tsv"):
                E.append(_table(fp, "summary/%s" % fn, "Sample"))
            else:
                E.append(_file(fp, "summary/%s" % fn))
'''

c3 = src.count(OLD3)
assert c3 == 1, 'OLD3 (c_03 summary block) count = %d (expect 1)' % c3
src = src.replace(OLD3, NEW3, 1)
print('OK[3]: c_03 summary 子目录已改为样本级过滤（命中 %d 处）' % c3)

io.open(P, 'w', encoding='utf-8').write(src)
print('DONE: extract_barbarum_complete.py 已更新')
