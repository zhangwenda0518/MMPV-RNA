#!/usr/bin/env python3
"""
Process Guard — ensures clean termination of the entire process tree on Ctrl+C.

When the orchestrator (auto_known_virus.py) spawns subprocesses that themselves
spawn grandchildren (megahit, bwa, samtools, etc.), a plain Ctrl+C only kills
the parent. This module installs signal handlers that recursively terminate the
entire process tree, preventing orphaned CPU-hungry background processes.

Usage:
    from process_guard import ProcessGuard

    guard = ProcessGuard()
    guard.install()           # install signal handlers + start tracking

    # ... run pipeline stages, subprocesses auto-tracked ...

    guard.cleanup()           # optional: remove handlers when done

Or as a context manager:
    with ProcessGuard() as guard:
        guard.run(['python', 'batch_virus_depth.py', ...])
"""

import os
import sys
import signal
import subprocess
import threading
import time
from pathlib import Path
from typing import Optional, List, Set

# ── Platform detection ──
IS_WINDOWS = sys.platform == "win32"
IS_LINUX = sys.platform.startswith("linux")

# ── Try psutil for reliable cross-platform process tree walking ──
try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


def _find_child_pids(pid: int) -> Set[int]:
    """Recursively find all child PIDs for a given parent PID.
    
    On Windows without psutil: uses wmic (slower but reliable).
    On Linux without psutil: walks /proc.
    """
    children = set()

    if HAS_PSUTIL:
        try:
            parent = psutil.Process(pid)
            for child in parent.children(recursive=True):
                try:
                    children.add(child.pid)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        return children

    # ── Fallback: OS-specific process tree walking ──
    if IS_WINDOWS:
        try:
            result = subprocess.run(
                ["wmic", "process", "where", f"(ParentProcessId={pid})", "get", "ProcessId"],
                capture_output=True, text=True, timeout=5
            )
            for line in result.stdout.strip().split("\n")[1:]:
                line = line.strip()
                if line.isdigit():
                    child_pid = int(line)
                    children.add(child_pid)
                    children.update(_find_child_pids(child_pid))
        except Exception:
            pass

    elif IS_LINUX:
        try:
            for proc_dir in Path("/proc").iterdir():
                if not proc_dir.name.isdigit():
                    continue
                stat_file = proc_dir / "stat"
                if not stat_file.exists():
                    continue
                try:
                    stat = stat_file.read_text()
                    # Format: pid (comm) state ppid ...
                    ppid_str = stat.split(") ", 1)[1].split()[1] if ") " in stat else ""
                    if ppid_str.isdigit() and int(ppid_str) == pid:
                        child_pid = int(proc_dir.name)
                        children.add(child_pid)
                        children.update(_find_child_pids(child_pid))
                except Exception:
                    pass
        except Exception:
            pass

    return children


def _kill_process_tree(pid: int, sig: int = signal.SIGTERM, timeout: float = 5.0):
    """Kill a process and all its descendants.
    
    Sends SIGTERM first, then SIGKILL after timeout if processes remain.
    """
    if IS_WINDOWS:
        # taskkill /T kills the entire tree natively
        try:
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(pid)],
                capture_output=True, timeout=timeout
            )
        except Exception:
            pass
        return

    # ── POSIX: find children and kill recursively ──
    children = _find_child_pids(pid)
    children.add(pid)

    # Phase 1: SIGTERM (graceful)
    for child_pid in sorted(children, reverse=True):
        try:
            os.kill(child_pid, sig)
        except (ProcessLookupError, PermissionError):
            pass

    # Phase 2: wait briefly, then SIGKILL stubborn survivors
    deadline = time.time() + timeout
    while time.time() < deadline:
        alive = set()
        for child_pid in children:
            try:
                os.kill(child_pid, 0)  # signal 0 = check if process exists
                alive.add(child_pid)
            except (ProcessLookupError, PermissionError):
                pass
        if not alive:
            break
        time.sleep(0.3)

    for child_pid in alive:
        try:
            os.kill(child_pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass


class ProcessGuard:
    """Track and clean up spawned subprocesses on abnormal termination.
    
    Features:
    - Installs SIGINT/SIGTERM handlers that kill the entire process tree
    - Auto-tracks subprocesses spawned via guard.run()
    - Works as a context manager for scoped cleanup
    - Cross-platform (Windows taskkill, POSIX os.kill + /proc walk)
    """

    def __init__(self, logger=None):
        self._logger = logger
        self._spawned_pids: Set[int] = set()
        self._lock = threading.Lock()
        self._original_sigint = None
        self._original_sigterm = None
        self._installed = False
        self._shutting_down = False

    def install(self):
        """Install signal handlers. Call once at pipeline startup."""
        if self._installed:
            return

        self._original_sigint = signal.signal(signal.SIGINT, self._on_signal)
        self._original_sigterm = signal.signal(signal.SIGTERM, self._on_signal)
        self._installed = True

        if self._logger:
            self._logger.info("[Guard] Process guard installed (SIGINT/SIGTERM → kill tree)")

    def uninstall(self):
        """Restore original signal handlers."""
        if not self._installed:
            return
        if self._original_sigint:
            signal.signal(signal.SIGINT, self._original_sigint)
        if self._original_sigterm:
            signal.signal(signal.SIGTERM, self._original_sigterm)
        self._installed = False

    def _on_signal(self, signum, frame):
        """Handle termination signal: kill all tracked children + self."""
        if self._shutting_down:
            return  # prevent re-entrant signals during cleanup
        self._shutting_down = True

        sig_name = signal.Signals(signum).name
        if self._logger:
            self._logger.warning("[Guard] Received %s — terminating entire process tree...", sig_name)
        else:
            print(f"\n[Guard] Received {sig_name} — terminating entire process tree...", file=sys.stderr)

        # Kill all tracked children
        with self._lock:
            pids = list(self._spawned_pids)

        for child_pid in pids:
            try:
                _kill_process_tree(child_pid)
            except Exception:
                pass

        # Also kill our own process group as a fallback
        try:
            if IS_WINDOWS:
                _kill_process_tree(os.getpid())
            else:
                os.killpg(os.getpgid(0), signal.SIGKILL)
        except Exception:
            pass

        # Restore original handler and re-raise
        self.uninstall()
        if self._logger:
            self._logger.warning("[Guard] Cleanup complete. Exiting.")
        sys.exit(128 + signum)

    def track_pid(self, pid: int):
        """Explicitly register a PID for cleanup tracking."""
        with self._lock:
            self._spawned_pids.add(pid)

    def run(self, cmd, **kwargs) -> subprocess.CompletedProcess:
        """Run a subprocess with automatic process tree tracking.
        
        Args:
            cmd: Command to run (str for shell=True, list otherwise)
            **kwargs: Passed to subprocess.run()
        
        Returns:
            subprocess.CompletedProcess instance
        
        On Windows: uses taskkill /T for cleanup.
        On Linux: tracks PID and recursively kills children.
        
        If the guard receives a signal during execution, the subprocess
        and all its descendants will be terminated.
        """
        is_shell = kwargs.get("shell", isinstance(cmd, str))
        proc = None

        try:
            # Start the process so we can capture its PID
            if is_shell:
                proc = subprocess.Popen(
                    cmd,
                    shell=True,
                    stdout=kwargs.get("stdout", subprocess.PIPE),
                    stderr=kwargs.get("stderr", subprocess.PIPE),
                    cwd=kwargs.get("cwd"),
                    text=kwargs.get("text", True),
                )
            else:
                proc = subprocess.Popen(
                    cmd,
                    stdout=kwargs.get("stdout", subprocess.PIPE),
                    stderr=kwargs.get("stderr", subprocess.PIPE),
                    cwd=kwargs.get("cwd"),
                    text=kwargs.get("text", True),
                )

            if proc.pid:
                self.track_pid(proc.pid)

            # Wait for completion
            stdout, stderr = proc.communicate(
                timeout=kwargs.get("timeout")
            )

            # Build CompletedProcess-like result
            return subprocess.CompletedProcess(
                args=cmd,
                returncode=proc.returncode,
                stdout=stdout or "",
                stderr=stderr or "",
            )

        except subprocess.TimeoutExpired:
            if proc and proc.pid:
                _kill_process_tree(proc.pid)
            proc.wait()
            raise
        except KeyboardInterrupt:
            if proc and proc.pid:
                _kill_process_tree(proc.pid)
            raise
        finally:
            if proc and proc.pid:
                with self._lock:
                    self._spawned_pids.discard(proc.pid)

    def cleanup(self):
        """Kill all tracked processes and uninstall handlers.
        
        Call at normal pipeline completion or in a finally block.
        """
        with self._lock:
            pids = list(self._spawned_pids)

        for child_pid in pids:
            try:
                _kill_process_tree(child_pid)
            except Exception:
                pass

        with self._lock:
            self._spawned_pids.clear()

        self.uninstall()

    def __enter__(self):
        self.install()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False  # don't suppress exceptions


# ── Convenience: monkey-patch subprocess.run for global tracking ──
_ORIGINAL_SUBPROCESS_RUN = subprocess.run
_GLOBAL_GUARD: Optional[ProcessGuard] = None


def _patched_subprocess_run(cmd, **kwargs):
    """Replacement for subprocess.run() that goes through the guard."""
    if _GLOBAL_GUARD and _GLOBAL_GUARD._installed:
        # Convert check=False semantics since our guard.run always returns
        check = kwargs.pop("check", False)
        result = _GLOBAL_GUARD.run(cmd, **kwargs)
        if check and result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode, cmd,
                output=result.stdout, stderr=result.stderr
            )
        return result
    return _ORIGINAL_SUBPROCESS_RUN(cmd, **kwargs)


def install_global_guard(logger=None) -> ProcessGuard:
    """Install a global process guard that intercepts all subprocess.run() calls.
    
    Returns:
        ProcessGuard instance (already installed).
        The guard is accessible as the global _GLOBAL_GUARD.
    """
    global _GLOBAL_GUARD
    if _GLOBAL_GUARD is not None:
        return _GLOBAL_GUARD

    _GLOBAL_GUARD = ProcessGuard(logger=logger)
    _GLOBAL_GUARD.install()
    subprocess.run = _patched_subprocess_run
    return _GLOBAL_GUARD


def uninstall_global_guard():
    """Restore original subprocess.run and remove the global guard."""
    global _GLOBAL_GUARD
    if _GLOBAL_GUARD:
        _GLOBAL_GUARD.cleanup()
        _GLOBAL_GUARD = None
    subprocess.run = _ORIGINAL_SUBPROCESS_RUN
