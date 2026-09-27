"""进程守护 (process_guard)

解决"停了主命令、后台子进程(spades/kraken2/seqtk/worker)还在跑"的问题。

原理: 收到 SIGINT / SIGTERM 时, 一次性递归收集本进程的所有子孙进程 PID,
先 SIGTERM 优雅退出、宽限几秒、再 SIGKILL 兜底, 然后主进程退出。
- 一次性收集全部 PID (而不是分两次)可避免"子进程被过继给 init 后逃逸"的问题。
- 只处理 SIGINT/SIGTERM, 不碰 SIGHUP → 仍兼容 nohup(断开终端后台继续运行)。

用法: 入口脚本 main() 开头调用
    import process_guard
    process_guard.install(logger)   # logger 可选
停止整条流水线: 直接 kill <主进程PID> 即可, 子孙会被连带清理。
"""
import os
import sys
import signal
import time

_installed = False


def _iter_descendants(root_pid):
    """从 /proc 读取 ppid 关系, 返回 root_pid 的所有子孙 PID (BFS)。"""
    parent_of = {}
    try:
        entries = os.listdir('/proc')
    except OSError:
        return []
    for name in entries:
        if not name.isdigit():
            continue
        try:
            with open('/proc/%s/status' % name) as fh:
                for line in fh:
                    if line.startswith('PPid:'):
                        parent_of[int(name)] = int(line.split()[1])
                        break
        except (IOError, OSError, ValueError):
            continue
    child_of = {}
    for pid, ppid in parent_of.items():
        child_of.setdefault(ppid, []).append(pid)
    result, stack = [], [root_pid]
    while stack:
        p = stack.pop()
        for c in child_of.get(p, []):
            result.append(c)
            stack.append(c)
    return result


def install(logger=None, grace=2):
    """安装信号处理器: 停止时连带杀掉所有子孙进程。grace = SIGTERM 到 SIGKILL 的宽限秒数。"""
    global _installed
    if _installed:
        return
    _installed = True

    def _log(msg):
        if logger is not None:
            try:
                logger.warning(msg)
                return
            except Exception:
                pass
        sys.stderr.write(msg + "\n")
        sys.stderr.flush()

    def _handler(signum, frame):
        # 一次性快照全部子孙 PID (避免子进程被过继给 init 后逃逸)
        snapshot = _iter_descendants(os.getpid())
        _log("[process_guard] 收到信号 %d, 正在终止 %d 个子孙进程..." % (signum, len(snapshot)))
        for pid in snapshot:
            try:
                os.kill(pid, signal.SIGTERM)
            except OSError:
                pass
        time.sleep(grace)
        # SIGKILL 兜底: 对 (快照 ∪ 当前仍存在的子孙) 反复清理, 直到全部确认死亡或超过 10 轮
        for _ in range(10):
            alive = []
            for pid in set(snapshot) | set(_iter_descendants(os.getpid())):
                try:
                    os.kill(pid, 0)  # 探测存活 (不实际发信号)
                    alive.append(pid)
                except OSError:
                    pass
            if not alive:
                break
            for pid in alive:
                try:
                    os.kill(pid, signal.SIGKILL)
                except OSError:
                    pass
            time.sleep(0.3)
        _log("[process_guard] 子进程已清理, 退出。")
        os._exit(1)

    for s in (signal.SIGINT, signal.SIGTERM):
        try:
            signal.signal(s, _handler)
        except (ValueError, OSError):
            pass
