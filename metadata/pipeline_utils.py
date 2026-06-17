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
    def banner(title):
        print(f"""{UI.C['purple']}{UI.C['bold']}
  =============================================================
    {title}
  ============================================================={UI.C['reset']}""")

    @classmethod
    def stage(cls, name, status="start"):
        sym = {'start': '[>>>]', 'skip': '[---]', 'done': '[OK]', 'fail': '[FAIL]'}
        color = {'start': cls.C['cyan'], 'skip': cls.C['gray'],
                 'done': cls.C['green'], 'fail': cls.C['red']}
        s, c = sym.get(status, '[>>>]'), color.get(status, cls.C['cyan'])
        print(f"\n{c}{cls.C['bold']}{'=' * 60}{cls.C['reset']}")
        print(f"{c}{cls.C['bold']} {s} {name}{cls.C['reset']}")
        print(f"{c}{cls.C['bold']}{'=' * 60}{cls.C['reset']}\n")

    @classmethod
    def ok(cls, msg):
        print(f"  {cls.C['green']}[OK]{cls.C['reset']} {msg}")

    @classmethod
    def warn(cls, msg):
        print(f"  {cls.C['yellow']}[WARN]{cls.C['reset']} {msg}")

    @classmethod
    def err(cls, msg):
        print(f"  {cls.C['red']}[ERR]{cls.C['reset']} {msg}")

    @classmethod
    def info(cls, msg):
        print(f"  {cls.C['gray']}> {msg}{cls.C['reset']}")


# ==========================================
# Checkpoint
# ==========================================
class Checkpoint:
    def __init__(self, work_dir):
        self.dir = os.path.join(work_dir, '.checkpoints')
        os.makedirs(self.dir, exist_ok=True)
        self.file = os.path.join(self.dir, 'state.json')
        self.state = self._load()

    def _load(self):
        if os.path.exists(self.file):
            with open(self.file, 'r') as f:
                return json.load(f)
        return {'stages': {}}

    def _save(self):
        with open(self.file, 'w') as f:
            json.dump(self.state, f, indent=2, ensure_ascii=False)

    def is_done(self, stage):
        return self.state.get('stages', {}).get(stage, {}).get('status') == 'done'

    def mark_start(self, stage):
        self.state['stages'][stage] = {'status': 'running', 'started': datetime.now().isoformat()}
        self._save()

    def mark_done(self, stage):
        self.state['stages'][stage] = {'status': 'done', 'completed': datetime.now().isoformat()}
        self._save()

    def mark_fail(self, stage, err=''):
        self.state['stages'][stage] = {'status': 'failed', 'error': str(err)[:200]}
        self._save()

    def reset(self):
        self.state['stages'] = {}
        self._save()

    def summary(self, stage_order):
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
def run_cmd(cmd, stage_name, log_dir, timeout=86400, secrets=None):
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
