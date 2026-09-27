import numpy as np

# 历史坑修复 (2026-08-30 审计实锤): 旧版 from utils.ess import calculate_ess as _ess_geyer
# 之后又用本地 def _ess_geyer(x) 原地遮蔽了 import → bulk/tail ESS 实际仍在用旧口径
# (有偏 acov / 配对从 lag1 起)，与主 ESS 同链两个数。现真正统一到 utils.ess。
from utils.ess import calculate_ess as _ess_geyer  # noqa: F811


def ess_bulk(x):
    """bulk-ESS: 复用 utils.ess.calculate_ess (Tracer 口径), 与主 ESS 同源可比。"""
    return _ess_geyer(np.asarray(x, float))


def ess_tail(x, probs=(0.05, 0.95)):
    """tail-ESS (Vehtari et al. 2021): 对尾部概率分位数的指示序列分别算 ESS, 取最小。
    检查分布尾部 (HPD 端点) 的采样充分性, 比 bulk 保守。"""
    x = np.asarray(x, float)
    worst = float("inf")
    for q in probs:
        thr = np.quantile(x, q)
        ess_q = _ess_geyer((x <= thr).astype(float))
        worst = min(worst, ess_q)
    return worst


def ess_evolution(x, n_points=10, probs=(0.05, 0.95)):
    """ESS 随链长增长曲线: 分段累积采样计算 bulk/tail ESS。
    返回 (fracs, bulk_list, tail_list)。曲线趋平且过 200 = 链长足够。"""
    x = np.asarray(x, float)
    n = len(x)
    fracs = np.linspace(0.1, 1.0, n_points)
    bulk_l, tail_l = [], []
    for f in fracs:
        seg = x[: max(int(n * f), 100)]
        bulk_l.append(ess_bulk(seg))
        tail_l.append(ess_tail(seg, probs))
    return fracs, bulk_l, tail_l
