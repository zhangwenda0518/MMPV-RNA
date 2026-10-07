#!/usr/bin/env python3
"""sample_floor.py — MCMC 样本下限的统一登记处 (2026-10-03)

历史问题: 「样本太少」的下限阈值散落在多处、口径不一, 而且大多是**静默降级** ——
退化输入 (链未跑完 / logcombiner burnin 超界 → 近空 merged.log) 产出的不是报错,
而是"看似正常的空结果":

    · utils/ess.py::calculate_ess        < 4    → 直接返回样本数 (假 ESS: 1 行日志 → 全员 ESS=1)
    · utils/beast_parser.py              < 10   → 该参数静默不出值
    · utils/auto_parse_params.py         < 100  → 参数汇总为空 JSON (1 字节级), 无告警
    · utils/beast_bridge.py               (闸门)  < 10 → 拒出 ESS 表 (2026-10-03 新增)
    · utils/virphy_bridge.py BSP R 脚本   (闸门)  < 10 → 拒出 skyline 图, 退出码 2 (2026-10-03 新增)

统一口径 (本模块即唯一来源; 调用方不要各写各的):
    MIN_MCMC_SAMPLES          = 10   统计/出图下限 (burn-in 之后可用样本数)
    MIN_PARAM_SUMMARY_SAMPLES = 100  参数汇总 (mean/HPD/ESS) 下限

为什么不做成"硬失败": 链条 pending/短跑在工程上是常态, 这里只保证
**要么有结果, 要么有明确告警**, 不允许静默产出空结果。
"""

MIN_MCMC_SAMPLES = 10
MIN_PARAM_SUMMARY_SAMPLES = 100


def require_min_samples(n, what, floor=MIN_MCMC_SAMPLES, log=None):
    """检查样本数是否达到下限。

    Parameters
    ----------
    n : int            实际可用样本数 (burn-in 之后)
    what : str         步骤名 (进告警措辞, 例: "[ess] BEAST 合并日志")
    floor : int        下限 (默认 MIN_MCMC_SAMPLES)
    log : LogCollector, optional   给了就写它的 error(); 没给就 print

    Returns
    -------
    bool : True=达标; False=不足 (已打统一措辞的告警)
    """
    if n >= floor:
        return True
    msg = (f"{what}: 仅 {n} 个可用样本 (<{floor}) — 链未跑完 / 合并 burnin 超界? "
           f"该步结果不可用, 已拒出")
    if log is not None:
        for attr in ("error", "warning", "emit"):
            fn = getattr(log, attr, None)
            if callable(fn):
                fn(msg)
                break
    else:
        print(f"  ⚠️ {msg}")
    return False
