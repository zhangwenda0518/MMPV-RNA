#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline Utilities — shared UI, Checkpoint, and shell execution.
Used by public_data_pipeline.py and build_host_pipeline.py.
"""

import os
import sys
import json
import shutil
import subprocess
from datetime import datetime
from typing import Any, Dict, List, Optional, Iterable

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


# ==========================================
# UI
# ==========================================
class UI:
    C = {
        'cyan': '\033[96m', 'green': '\033[92m', 'yellow': '\033[93m',
        'red': '\033[91m', 'purple': '\033[95m', 'gray': '\033[90m',
        'bold': '\033[1m', 'reset': '\033[0m',
    }

    @staticmethod
    def banner(title: str) -> None:
        print(f"""{UI.C['purple']}{UI.C['bold']}
  =============================================================
    {title}
  ============================================================={UI.C['reset']}""")

    @classmethod
    def stage(cls, name: str, status: str = "start") -> None:
        sym = {'start': '[>>>]', 'skip': '[---]', 'done': '[OK]', 'fail': '[FAIL]'}
        color = {'start': cls.C['cyan'], 'skip': cls.C['gray'],
                 'done': cls.C['green'], 'fail': cls.C['red']}
        s, c = sym.get(status, '[>>>]'), color.get(status, cls.C['cyan'])
        print(f"\n{c}{cls.C['bold']}{'=' * 60}{cls.C['reset']}")
        print(f"{c}{cls.C['bold']} {s} {name}{cls.C['reset']}")
        print(f"{c}{cls.C['bold']}{'=' * 60}{cls.C['reset']}\n")

    @classmethod
    def ok(cls, msg: str) -> None:
        print(f"  {cls.C['green']}[OK]{cls.C['reset']} {msg}")

    @classmethod
    def warn(cls, msg: str) -> None:
        print(f"  {cls.C['yellow']}[WARN]{cls.C['reset']} {msg}")

    @classmethod
    def err(cls, msg: str) -> None:
        print(f"  {cls.C['red']}[ERR]{cls.C['reset']} {msg}")

    @classmethod
    def info(cls, msg: str) -> None:
        print(f"  {cls.C['gray']}> {msg}{cls.C['reset']}")


# ==========================================
# Checkpoint
# ==========================================
class Checkpoint:
    """JSON-based checkpoint manager for pipelined stage resumption."""

    def __init__(self, work_dir: str) -> None:
        self.dir: str = os.path.join(work_dir, '.checkpoints')
        os.makedirs(self.dir, exist_ok=True)
        self.file: str = os.path.join(self.dir, 'state.json')
        self.state: Dict[str, Any] = self._load()

    def _load(self) -> Dict[str, Any]:
        if not os.path.exists(self.file):
            return {'stages': {}}
        try:
            with open(self.file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            if not isinstance(data, dict):
                raise ValueError(f"根节点类型异常: {type(data).__name__}")
            if not isinstance(data.get('stages'), dict):
                data['stages'] = {}
            return data
        except (json.JSONDecodeError, OSError, ValueError) as e:
            # 上次写入中途被杀 / 磁盘错误 → 留档坏文件后重新初始化, 绝不阻塞管线
            backup = f"{self.file}.corrupt-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
            try:
                os.replace(self.file, backup)
                kept = f"已留档: {backup}"
            except OSError:
                kept = f"留档失败, 已忽略: {self.file}"
            print(f"  [WARNING] checkpoint 状态文件损坏或不可读 ({e}); "
                  f"{kept}; 已重新初始化为空状态")
            return {'stages': {}}

    def _save(self) -> None:
        # 原子写: 同目录临时文件 (带 pid 防并发碰撞) → os.replace, 不留半截 JSON
        tmp = f"{self.file}.tmp.{os.getpid()}"
        try:
            with open(tmp, 'w', encoding='utf-8') as f:
                json.dump(self.state, f, indent=2, ensure_ascii=False)
                f.flush()
                os.fsync(f.fileno())
            os.replace(tmp, self.file)
        except BaseException:
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except OSError:
                pass
            raise

    def is_done(self, stage: str) -> bool:
        return self.state.get('stages', {}).get(stage, {}).get('status') == 'done'

    def mark_start(self, stage: str) -> None:
        self.state['stages'][stage] = {'status': 'running', 'started': datetime.now().isoformat()}
        self._save()

    def mark_done(self, stage: str) -> None:
        self.state['stages'][stage] = {'status': 'done', 'completed': datetime.now().isoformat()}
        self._save()

    def mark_fail(self, stage: str, err: str = '') -> None:
        self.state['stages'][stage] = {'status': 'failed', 'error': str(err)[:200]}
        self._save()

    def reset(self) -> None:
        self.state['stages'] = {}
        self._save()

    def summary(self, stage_order: List[str]) -> str:
        lines = []
        for s in stage_order:
            info = self.state.get('stages', {}).get(s, {})
            status = info.get('status', 'pending')
            mark = {'done': '[OK]', 'running': '[..]', 'failed': '[FAIL]'}.get(status, '[  ]')
            lines.append(f"  {mark} {s}: {status}")
        return '\n'.join(lines)


# ==========================================
# Shell execution
# ==========================================
def run_cmd(
    cmd: str,
    stage_name: str,
    log_dir: str,
    timeout: int = 86400,
    secrets: Optional[Iterable[str]] = None,
) -> int:
    """Execute a shell command with output logging.

    Args:
        cmd: Full shell command string.
        stage_name: Label for log file naming.
        log_dir: Directory for log output.
        timeout: Max wall-clock seconds.
        secrets: Iterable of strings to mask in the log file.
    """
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{stage_name}.log")

    log_cmd = cmd
    if secrets:
        for s in secrets:
            if s:
                log_cmd = log_cmd.replace(s, '***')

    bash = shutil.which('bash') or '/bin/bash' if sys.platform != 'win32' else None
    # Windows: shlex.quote produces POSIX single quotes; cmd.exe needs double quotes
    if sys.platform == 'win32' and bash is None:
        cmd = cmd.replace("'", '"')
    with open(log_file, 'w', encoding='utf-8') as lf:
        lf.write(f"# Stage: {stage_name}\n"
                 f"# Cmd: {log_cmd}\n"
                 f"# {datetime.now().isoformat()}\n"
                 f"{'=' * 60}\n\n")
        proc = subprocess.Popen(cmd, shell=True, executable=bash,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                encoding='utf-8', errors='replace', bufsize=1,
                                env={**os.environ, 'PYTHONIOENCODING': 'utf-8'})
        for line in iter(proc.stdout.readline, ''):
            line_out = line
            if secrets:
                for s in secrets:
                    if s:
                        line_out = line_out.replace(s, '***')
            lf.write(line_out)
            s = line.rstrip()
            if s:
                if len(s) > 150:
                    s = s[:150] + "..."
                print(f"    {UI.C['gray']}{s}{UI.C['reset']}")
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            return -1
    return proc.returncode
