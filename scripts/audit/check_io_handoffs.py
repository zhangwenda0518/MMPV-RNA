#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_io_handoffs.py — MMPV 管线间 I/O 衔接审计
================================================
扫描六大管线 + GUI + mmpv_common 的 Python 源码, 找出**功能层**硬编码的
legacy 目录名 (00a_CleanData / 01_detection / 08_Rescue ...), 防止
standard 布局 (doc/IO_LAYOUT_DESIGN.md) 下上下游衔接断裂。

判定分级:
  FAIL  精确字面量   字符串常量整串等于 legacy 目录名 (路径构件),
                    或相对子路径以 legacy 目录名开头 ('01_detection/xxx')
  WARN 长文本包含   帮助文本/HTML 模板/日志文案中提到目录名 (可接受,
                    但 standard 语境下建议措辞两布局中性)
白名单 (按路径前缀):
  mmpv_common/io_layout.py      注册表本体 (legacy 名的唯一定义点)
  mmpv_common/tests/            冻结测试 (legacy 名的守卫)
  (phylo_results/ 一次性诊断脚本 175+2 个已于 2026-09-28 归档至
   archive/20260928_phylo_results_diag/; 目录现仅存审计 .md 与结果数据)
  virome_analysis_pipeline/_debug_stage1.py   调试脚本
  (utils/auto_known_virus.py 已于 2026-09-28 执行 20260830 审计决定,
   归档为 archive/auto_known_virus.orphan_utils_20260915.py; archive/ 本身不扫描)
  virome_discovery_pipeline/integrated_summary.py     审计专用 (仅 scripts/audit 引用)
  endogenous_virus_pipeline/tests/      legacy fixture 测试

用法:
  python scripts/audit/check_io_handoffs.py           # 打印报告, FAIL 即退出码 1
  python scripts/audit/check_io_handoffs.py --quiet   # 只打印 FAIL/WARN 摘要
"""

import argparse
import ast
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

LEGACY_DIRS = [
    "00a_CleanData", "00b_HostDepletion", "00c_BBnorm", "01_Assembly",
    "02a_Identification", "02b_Filter", "03a_COBRA", "03b_MergeSamples",
    "04_CLUSTER", "05_Taxonomy", "06_HostPrediction", "07_Checkv",
    "08_Rescue", "09a_Virome_Analysis", "09_Virome_Analysis",
    "09b_Analysis_Verify", "10_Reports",
    "01_detection", "02_filtering", "03_variants", "04_post_analysis",
    "05_assembly", "06_extraction", "07_similarity", "08_dvg", "09_report",
    "09_dvg",
]
SCOPES = [
    "data_preprocessing_pipeline", "virome_discovery_pipeline",
    "virome_analysis_pipeline", "public_metadata_pipeline",
    "virome_phylo_pipeline", "endogenous_virus_pipeline",
    "metadata_gui", "submission_gui", "mmpv_common",
]
WHITELIST = (
    "mmpv_common/io_layout.py", "mmpv_common/tests/",
    "phylo_results/", "_debug_stage1",
    "integrated_summary.py",
    "endogenous_virus_pipeline/tests/",
)


def docstring_lines(tree):
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef,
                             ast.AsyncFunctionDef, ast.ClassDef)):
            if node.body and isinstance(node.body[0], ast.Expr) \
                    and isinstance(node.body[0].value, ast.Constant):
                c = node.body[0]
                out.update(range(c.lineno, (c.end_lineno or c.lineno) + 1))
    return out


def scan(quiet: bool):
    fails, warns = [], []
    for scope in SCOPES:
        for p in sorted(Path(REPO, scope).rglob("*.py")):
            rel = str(p.relative_to(REPO)).replace("\\", "/")
            if any(x in rel for x in ("pycache", "/archive/", "/build/",
                                      "/dist/", "server_checks")):
                continue
            if any(w in rel for w in WHITELIST):
                continue
            try:
                tree = ast.parse(p.read_text(encoding="utf-8", errors="replace"))
            except SyntaxError:
                continue
            skip = docstring_lines(tree)
            for node in ast.walk(tree):
                if not (isinstance(node, ast.Constant)
                        and isinstance(node.value, str)):
                    continue
                if node.lineno in skip:
                    continue
                val = node.value
                exact = [n for n in LEGACY_DIRS if val == n]
                prefix = [n for n in LEGACY_DIRS if val.startswith(n + "/")]
                contains = [n for n in LEGACY_DIRS if n in val]
                if exact or prefix:
                    fails.append((rel, node.lineno, exact or prefix,
                                  val.strip()[:70]))
                elif contains:
                    warns.append((rel, node.lineno, contains,
                                  val.strip()[:70]))
    return fails, warns


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--quiet", action="store_true", help="只打印 FAIL 行")
    args = ap.parse_args()

    fails, warns = scan(args.quiet)

    if not args.quiet:
        print(f"== MMPV I/O 衔接审计 (legacy 目录名功能残留扫描) ==")
        print(f"FAIL (功能硬编码, standard 布局会断裂): {len(fails)}")
        for rel, ln, names, val in fails:
            print(f"  FAIL {rel}:{ln} {names} {val!r}")
        print(f"WARN (长文本提及, 可接受): {len(warns)}")
        for rel, ln, names, val in warns:
            print(f"  warn {rel}:{ln} {names[:1]}")
    else:
        for rel, ln, names, val in fails:
            print(f"FAIL {rel}:{ln} {names}")

    print(f"-- 结论: {'PASS' if not fails else 'FAIL'} "
          f"({len(fails)} fail / {len(warns)} warn)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
