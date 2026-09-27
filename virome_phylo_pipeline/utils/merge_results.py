#!/usr/bin/env python3
"""一键合并工具: 检查完成 → logcombiner → treeannotator → 参数解析

用法:
  python merge_results.py --virus GCVA        # 按已知病毒名
  python merge_results.py --work-dir <path>   # 或指定工作目录
  python merge_results.py --status            # 只查看进度

正常流程: run_phylogeo.py 提交 → BEAST 跑完自动退出 → 本命令合并
"""
import argparse, json, os, re, subprocess, sys, time
from pathlib import Path
import sys as _sys

_here = os.path.dirname(os.path.abspath(__file__))
_sys.path.insert(0, _here)
_sys.path.insert(0, os.path.join(_here, 'utils'))
# 历史坑: 缺根目录路径 → 独立运行 python utils/merge_results.py 时
# from utils.dataset_config 直接 ModuleNotFoundError
_sys.path.insert(0, os.path.dirname(_here))
from utils.dataset_config import dataset as get_dataset, known_dataset_names


def known_work_dir(name):
    """已知数据集的工作目录 = dir/work (从 datasets.yaml 读)。"""
    d = get_dataset(name)
    if d is None:
        return None
    return f'{d["dir"]}/{d["work"]}'


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_task(work_dir):
    sp = f'{work_dir}/run_status.json'
    if not os.path.exists(sp):
        return None
    return json.load(open(sp))


def read_states(work_dir, prefixes):
    """读取各链 state (从 BEAST 主 log 最后数字行)"""
    states = []
    for p in prefixes:
        lf = f'{work_dir}/{p}.log'
        if not os.path.exists(lf):
            states.append(0)
            continue
        last = 0
        with open(lf, errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                try:
                    last = int(float(parts[0]))
                except (ValueError, IndexError):
                    pass
        states.append(last)
    return states


def _pid_alive(pid) -> bool:
    """进程是否存活。POSIX 用 os.kill(pid,0); 无法判定时返回 True (不误判为已死)。"""
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return True
    if pid <= 0:
        return True
    if os.name == 'nt':
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return True
    return True


def _count_running(st) -> int:
    """本任务仍在运行的 BEAST 进程数。

    2026-09-15 修复: 旧版固定用 `ps aux | grep beast | grep java | wc -l` 统计
    **全机器** java beast 进程 → 并行跑第二个数据集时该值永不为 0, all_done 永假,
    前台等待空转到 stale_minutes 才放弃 (退回 pending)。
    现优先用 run_status.json 里记录的 PID 逐个探测; 无 PID (旧任务) 才回退全局统计。
    """
    pids = [p for p in (st.get('pids') or []) if p]
    if pids:
        return sum(1 for pid in pids if _pid_alive(pid))
    try:
        r = subprocess.run('ps aux | grep "[b]east" | grep java | wc -l',
                           shell=True, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=30)
        return int(r.stdout.strip() or 0)
    except Exception as e:
        log(f"⚠️ 进程计数失败 ({e}); 以 state 判定为准")
        return 0


def status(work_dir, verbose=True):
    st = load_task(work_dir)
    if not st:
        log(f"⚠️ {work_dir} 无 run_status.json (可能是手动提交的任务)")
        return None
    prefixes = st['prefixes']
    states = read_states(work_dir, prefixes)
    cl = st['chain_length']
    n_proc = _count_running(st)
    if verbose:
        for p, s in zip(prefixes, states):
            pct = s / cl * 100
            print(f"  {p}: {s/1e6:.1f}M/{cl/1e6:.0f}M ({pct:.0f}%)")
        log(f"进程={n_proc}, 完成={sum(1 for s in states if s >= cl)}/{len(prefixes)}")
    return states, n_proc


def _pids_verifiable(st) -> bool:
    """能否用 PID 可靠判断链是否存活 (POSIX 且 run_status.json 里记了 pid)。"""
    return os.name != 'nt' and bool([p for p in (st.get('pids') or []) if p])


def _logs_quiet(work_dir, prefixes, quiet_sec: float = None) -> bool:
    """所有链日志在 quiet_sec 秒内都没有被追加写入。

    2026-09-15 修复 (P1-C): PID 心跳在"链被提交到另一台机器/容器"时不可见
    (os.kill 抛 ProcessLookupError -> 被判定已死), 会让 all_done 在链仍在
    运行时提前为真, 从而合并**不完整**的 .trees。产物本身不会骗人: 只要还有
    日志在被追加, 就绝不算完成。故在所有 PID 都探测不到时追加这条静默期判据。

    quiet_sec 默认 600s (10 min): 必须大于 BEAST 相邻两条 log 行的最大间隔,
    否则慢链 (logEvery 间隔 > quiet_sec) 会被误判为已完成。可用环境变量
    VP_MERGE_QUIET_SEC 覆盖 (例如把日志间隔很小的任务调成 120)。
    """
    import time as _time
    if quiet_sec is None:
        try:
            quiet_sec = float(os.environ.get('VP_MERGE_QUIET_SEC', '600'))
        except (TypeError, ValueError):
            quiet_sec = 600.0
    now = _time.time()
    for p in prefixes:
        fp = os.path.join(work_dir, f'{p}.log')
        try:
            if now - os.path.getmtime(fp) < quiet_sec:
                return False
        except OSError:
            return False        # 日志还没出现 -> 肯定没跑完
    return True


def all_done(work_dir, st):
    states, n_proc = status(work_dir, verbose=False)
    # 历史坑: BEAST 实际最后一行 state = 最大 logEvery 倍数 ≤ chain_length。
    # chain_length 非 logEvery 整数倍时 s < chain_length 恒成立 → 合并永久卡住。
    # 改判: 进程退出且 state 达到理论最后一行 (chain_length 向下取 logEvery 倍数)。
    # 2026-09-15 修复: log_every 可能为 None (XML 解析失败时 run_status.json 写 null),
    # 原 st.get('log_every', 5000) 在"键存在但值为 None"时返回 None → None % le 崩溃。
    le = st.get('log_every') or 5000
    try:
        le = int(le)
    except (TypeError, ValueError):
        le = 5000
    if le <= 0:
        le = 5000
    cl = int(st['chain_length'])
    target = cl - (cl % le) or cl
    prefixes = st.get('prefixes') or []
    # ① 进程层面: 能可靠探活且确实还在跑 -> 明确未完成
    if _pids_verifiable(st):
        if n_proc > 0:
            return False
    else:
        log("  (PID 不可验证: 非 POSIX 或旧任务缺 pids; 以日志静默期判定)")
    # ② 产物层面: 日志仍在被追加 -> 无论 PID 探测结果如何都判定未完成
    if prefixes and not _logs_quiet(work_dir, prefixes):
        log("  (链日志仍在更新 -> 判定未完成, 继续等待)")
        return False
    return all(s >= target for s in states)


def merge(work_dir, st, background=False):
    merged = f'{work_dir}/merged'
    os.makedirs(merged, exist_ok=True)
    prefixes = st['prefixes']
    # burnin 口径统一为链长 10% (与 run_phylogeo.MERGE_BURNIN_FRAC 一致; 旧 1% 已废弃)。
    # 默认 500000 是旧 1% 时代遗产, 改为 10% 兜底, 兼容缺 merge_burnin 字段的旧 run_status.json。
    burnin = st.get('merge_burnin') or max(int(st.get('chain_length', 5_000_000) * 0.10), 1)
    log(f"合并 (burnin={burnin})...")
    logs = [f'{work_dir}/{p}.log' for p in prefixes]
    trees = [f'{work_dir}/{p}.trees' for p in prefixes]

    # 1. log 合并 (快)
    # 历史坑: 此调用无 try/except, logcombiner 异常向上抛被管线吞 → 合并静默失败。
    # 且失败不阻断 trees 合并 (Python 合并不依赖 merged.log)。
    try:
        r = subprocess.run(['logcombiner', '-burnin', str(burnin)] + logs + [f'{merged}/merged.log'],
                           capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=3600)
    except (subprocess.TimeoutExpired, OSError) as _le:
        r = None
        log(f"  ⚠️ logcombiner 执行异常: {_le}")
    if r is not None and r.returncode == 0 and os.path.exists(f'{merged}/merged.log'):
        log(f"  ✅ merged.log ({os.path.getsize(f'{merged}/merged.log')/1e6:.1f}MB)")
    elif r is not None:
        log(f"  ❌ log 合并失败: {r.stderr[-200:]}")

    # 2. trees 合并 (Python: 按 STATE 过滤 burnin + #NEXUS 头)
    #    不用 logcombiner: BEAST1 LogCombiner 对注释树 (-burnin 按树数解释超界时)
    #    会异常慢/卡死 (实测 19MB 输入 8 分钟只写 5 字节)。treeannotator 也要求
    #    #NEXUS 头, 而裸 tree 行文件会报 "No trees"。
    log(f"  trees 合并 (Python)...")
    try:
        n_total = 0
        merged_trees = f'{merged}/merged.trees'
        with open(merged_trees, 'w') as fo:
            fo.write('#NEXUS\nBEGIN TREES;\n')
            for tp in trees:
                with open(tp, errors='replace') as f:
                    for line in f:
                        if line.startswith('tree STATE_'):
                            m = re.match(r'tree STATE_(\d+)', line)
                            state = int(m.group(1)) if m else 0
                            if state >= burnin:
                                fo.write(line)
                                n_total += 1
            fo.write('END;\n')
        log(f"  ✅ merged.trees ({n_total} trees, "
            f"{os.path.getsize(merged_trees)/1e6:.1f}MB)")

        # 3. TreeAnnotator (树多时抽样控制计算量: BEAST1 单线程, 大样本集很慢)
        annot_src = merged_trees
        if n_total > 3000:
            sampled = f'{merged}/merged_sampled.trees'
            step = max(1, n_total // 2000)
            n_sampled = 0
            with open(merged_trees) as fi, open(sampled, 'w') as fso:
                fso.write('#NEXUS\nBEGIN TREES;\n')
                for i, line in enumerate(fi):
                    if line.startswith('tree STATE_') and i % step == 0:
                        fso.write(line)
                        n_sampled += 1
                fso.write('END;\n')
            annot_src = sampled
            log(f"  → 抽样 {n_sampled} 棵树用于 treeannotator")
        r3 = subprocess.run(['treeannotator', '-burnin', '0', '-heights', 'median',
                             annot_src, f'{merged}/mcc.tree'],
                            capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=1800)
        if r3.returncode == 0 and os.path.exists(f'{merged}/mcc.tree'):
            log(f"  ✅ mcc.tree ({os.path.getsize(f'{merged}/mcc.tree')/1e6:.2f}MB)")
        else:
            log(f"  ❌ treeannotator 失败: {(r3.stderr or r3.stdout)[-200:]}")
    except Exception as e:
        log(f"  ❌ trees 合并异常: {e}")

    # 4. 参数解析 (若有)
    ap = f'{os.path.dirname(os.path.abspath(__file__))}/auto_parse_params.py'
    if os.path.exists(ap):
        # 用 -m utils.auto_parse_params (模块方式), 保证 utils.ess 可导入;
        # 脚本方式 python utils/xxx.py 会因 sys.path 指向 utils/ 而 import 失败
        subprocess.run([sys.executable, '-m', 'utils.auto_parse_params', work_dir, merged],
                       cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        if os.path.exists(f'{merged}/parameter_summary.json'):
            log(f"  ✅ parameter_summary.json")
    # 4. MCMC 收敛诊断图 (ESS 柱状 + trace/density/ACF + 双链散点)
    try:
        from utils.mcmc_diag_plot import plot_diagnostics, plot_chain_scatter
        diag_pdf = f'{merged}/mcmc_diagnostics.pdf'
        ess_values, trace_pdf, evo_pdf = plot_diagnostics(
            f'{merged}/merged.log', diag_pdf, burnin_frac=0.1)
        log("  ✅ mcmc_diagnostics.pdf (bulk/tail ESS: " +
            ", ".join(f"{p}={eb:.0f}/{et:.0f}" for p, _, _, eb, et in ess_values[:4]) + "…)")
        # 双链散点 (取第一条单链 log 对比合并 log)
        single_logs = [l for l in logs if os.path.basename(l) != 'merged.log']
        if single_logs:
            rs, _ = plot_chain_scatter(
                single_logs[0], f'{merged}/merged.log',
                f'{merged}/chain_scatter.pdf')
            log("  ✅ chain_scatter.pdf (r: " +
                ", ".join(f"{p}={r:.3f}" for p, r in rs[:3]) + "…)")
    except Exception as e:
        log(f"  ⚠️ mcmc diagnostics plot skipped: {e}")
    st['completed_at'] = time.strftime('%Y-%m-%d %H:%M:%S')
    json.dump(st, open(f'{work_dir}/run_status.json', 'w'), indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--virus', default=None, help='GCVA / PSTVd')
    ap.add_argument('--work-dir', default=None, help='或直接指定工作目录')
    ap.add_argument('--status', action='store_true', help='只查看进度')
    ap.add_argument('--background', action='store_true', help='trees 合并后台执行')
    args = ap.parse_args()

    work_dir = args.work_dir
    if args.virus:
        key = args.virus.upper()
        wd = known_work_dir(key)
        if wd is None:
            avail = ', '.join(known_dataset_names()) or '(无配置)'
            raise SystemExit(f"未知病毒: {args.virus} (可选: {avail})")
        work_dir = wd
    if not work_dir:
        raise SystemExit("需 --virus 或 --work-dir")

    st = load_task(work_dir)
    if args.status:
        if st:
            status(work_dir)
        else:
            log(f"无状态文件, 手动检查目录: {work_dir}")
        return

    if not st:
        raise SystemExit(f"无 run_status.json, 无法合并: {work_dir}")

    if not all_done(work_dir, st):
        states, n_proc = status(work_dir)
        log(f"❌ 链未全部完成 ({n_proc} 进程在跑), 等 BEAST 自然退出后再合并")
        sys.exit(1)

    merge(work_dir, st, background=args.background)


if __name__ == '__main__':
    main()
