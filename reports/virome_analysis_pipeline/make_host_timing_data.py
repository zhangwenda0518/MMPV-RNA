#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_host_timing_data.py —— 生成图③面板 b1「宿主预测器耗时」数据表。

数据源（服务器结果树，只读）：
  <base>/06_HostPrediction/host.log
      由 virome_discovery_pipeline/run_host_prediction.py::run_step() 写入，
      每步一行 `[OK] <step_name> — <n>s`。

口径（2026-10-14 实测，逐条核对时间戳）：
  宿主预测只跑过 **一次**（2026-08-12 01:08:32 → 01:19:50，`--mode all`，全量输入）。
  实测时间戳首尾相接，无缺口：
      01:08:32 RNAVirHost Classify Order      → 139 s
      01:10:51 RNAVirHost Predict             → 214 s   (01:10:51 = 01:08:32 + 139 ✓)
      01:14:25 PhaBOX2 CHERRY Host Prediction → 322 s   (01:14:25 = 01:10:51 + 214 ✓)
      01:19:47 ICTV Taxonomy Lookup (C9)      →   2 s   (01:19:47 = 01:14:25 + 322 ✓)
  ⚠ 2026-09-15 只重跑了 **ensemble 决策树**（`ensemble_host_summary.tsv` 被重写、
    留 `.bak_kill_20260915`），**三路预测器产物 mtime 仍是 2026-08-12**：
      RVH_result/  2026-08-12 01:14 ｜ phabox2_output/ 2026-08-12 01:19 ｜ C9_ICTV_result/ 01:19
    ⇒ **本次耗时对当前数据仍然有效**，不是陈旧值。

  ⚠⚠ 三个口径必须踩对：
    ① **RVH 是两步**（classify_order + predict）→ 图上 RVH 必须画**两步之和 353 s**，
       只画 predict 的 214 s 会低估 39%。
    ② **Rule 没有独立耗时**。它 = `Rule_Class_Taxonomy`，是 ensemble 决策树里的
       **分类级回退**，在 run_host_prediction.py **进程内**完成、不走 `run_step()`
       → host.log 里**不存在**它的 `[OK]` 行。
       从时间戳可反推：ICTV 于 01:19:47 起、耗时 2 s → 01:19:49 结束；
       日志末行 `[01:19:50] 宿主预测完成` → **ensemble 步 ≈ 1 s**。
       但这是**时间戳反推**而非日志实测值 → 本表 `wall_sec` 记 **NA**，
       并在图上以文字注明，**不画成 0 高柱**（对数轴上 0 无定义、且会被读成"实测"）。
    ③ `host.log` 的 `[Summary]` 段落（Plant 1092 / Unknown 453）是 **2026-08-12 那次**
       ensemble 的输出；当前报告口径来自 2026-09-15 重跑的 `ensemble_host_summary.tsv`
       （Plant 1063 / Unknown 482）→ **两者不同批**，耗时表**只取 `[OK]` 行、不取 Summary**。

用法：
  python make_host_timing_data.py                       # 用本地 upset_data/host.log
  python make_host_timing_data.py --fetch               # 先从 246 拉最新 host.log 再解析
输出：
  reports/virome_analysis_pipeline/upset_data/host_timing.tsv
      tool  wall_sec  steps  note
"""
import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "upset_data")
LOG = os.path.join(DATA, "host.log")
OUT = os.path.join(DATA, "host_timing.tsv")

HOST = "zhangwenda@202.119.189.246"
REMOTE_LOG = ("~/MMPV-paper/goji-virome/02_novel_virus/"
              "RNA-Lycium_barbarum_out/06_HostPrediction/host.log")

# step_name（host.log 里的原文）→ 图上的预测器名
STEP2TOOL = {
    "RNAVirHost Classify Order": "RVH",
    "RNAVirHost Predict": "RVH",
    "PhaBOX2 CHERRY Host Prediction": "PB2",
    "ICTV Taxonomy Lookup (C9)": "ICTV",
}
# 图上展示顺序 = 耗时降序（Rule 恒垫底，它没有实测值）
DISPLAY = ["RVH", "PB2", "ICTV", "Rule"]
LABEL = {
    "ICTV": "ICTV (C9 lookup)",
    "RVH": "RNAVirHost",
    "PB2": "PhaBOX2 (Cherry)",
    "Rule": "Rule (class fallback)",
}
NOTE = {
    "ICTV": "C9_ICTV_result/",
    "RVH": "RVH_result/",
    "PB2": "phabox2_output/",
    "Rule": "in-process, no separate wall time",
}

OK_RE = re.compile(r"^\[OK\]\s+(.+?)\s+—\s+(\d+)s\s*$")


def fetch():
    """从 246 拉最新 host.log（只读；本地 CRLF 转 LF）。"""
    print(f"[fetch] {HOST}:{REMOTE_LOG}")
    p = subprocess.run(["ssh", "-o", "ConnectTimeout=15", HOST, f"cat {REMOTE_LOG}"],
                       capture_output=True, text=True, encoding="utf-8")
    if p.returncode != 0:
        sys.exit(f"[ERR] ssh 失败: {p.stderr}")
    txt = p.stdout.replace("\r\n", "\n").replace("\r", "\n")
    os.makedirs(DATA, exist_ok=True)
    with open(LOG, "w", encoding="utf-8", newline="\n") as f:
        f.write(txt)
    print(f"[fetch] 写入 {LOG}（{len(txt)} 字节）")


def parse():
    if not os.path.isfile(LOG):
        sys.exit(f"[ERR] 缺 {LOG}；先跑 --fetch")
    per_tool = {}
    seen = []
    with open(LOG, encoding="utf-8") as f:
        for line in f:
            m = OK_RE.match(line.strip())
            if not m:
                continue
            step, sec = m.group(1), int(m.group(2))
            seen.append((step, sec))
            tool = STEP2TOOL.get(step)
            if tool is None:
                print(f"[WARN] 未登记的 step（已忽略）: {step!r} — {sec}s")
                continue
            per_tool.setdefault(tool, []).append((step, sec))

    if not per_tool:
        sys.exit("[ERR] host.log 里没有可识别的 [OK] 行")

    # 自检：时间戳首尾相接（每步起始 = 上一步起始 + 上一步耗时）
    ts_re = re.compile(r"^\[(\d{2}):(\d{2}):(\d{2})\]\s+(.+)$")
    stamps = []
    with open(LOG, encoding="utf-8") as f:
        for line in f:
            m = ts_re.match(line.strip())
            if m and m.group(4).strip() in STEP2TOOL:
                h, mi, s = (int(m.group(i)) for i in (1, 2, 3))
                stamps.append((h * 3600 + mi * 60 + s, m.group(4).strip()))
    for (t0, s0), (t1, s1) in zip(stamps, stamps[1:]):
        sec0 = dict(seen).get(s0)
        if sec0 is not None and t0 + sec0 != t1:
            print(f"[WARN] 时间戳不连续: {s0} 起 {t0}s + {sec0}s ≠ {s1} 起 {t1}s "
                  f"（差 {t1 - t0 - sec0}s）")

    rows = []
    for tool in DISPLAY:
        if tool == "Rule":
            rows.append((LABEL[tool], "NA", NOTE[tool]))
            continue
        steps = per_tool.get(tool)
        if not steps:
            print(f"[WARN] 无耗时记录: {tool}")
            continue
        total = sum(s for _, s in steps)
        detail = " + ".join(f"{n} {s}s" for n, s in steps)
        rows.append((LABEL[tool], str(total), f"{detail}｜{NOTE[tool]}"))

    os.makedirs(DATA, exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("tool\twall_sec\tsteps\n")
        for lab, sec, note in rows:
            f.write(f"{lab}\t{sec}\t{note}\n")
    print(f"[out] {OUT}")
    for lab, sec, note in rows:
        print(f"       {lab:26s} {sec:>5s} s   {note}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true", help="先从 246 拉最新 host.log")
    a = ap.parse_args()
    if a.fetch:
        fetch()
    parse()


if __name__ == "__main__":
    main()
