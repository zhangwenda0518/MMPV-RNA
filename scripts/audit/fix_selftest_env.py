#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""给 logic_selftest 的环境注入「由环境变量派生的常量」。

背景：自测框架只从补丁里提取 function 定义 + WANT_CONST 名单里的常量赋值；
v6.5 引入的 VOTE_WEIGHT_MODE 是 Sys.getenv(...)，不在提取范围，导致第三节
compute_tool_weights 直接 "object not found"。这里在提取段之后补默认值，
让自测能跑到底（不动补丁本身）。
"""
import io

P = "/tmp/logic_selftest_lf.R"
s = io.open(P, encoding="utf-8").read()

if "VOTE_WEIGHT_MODE" in s:
    print("已注入过，跳过")
    raise SystemExit(0)

old = '"CONSENSUS_MODE", "CASCADE_MIN_SHARE", "RATE_FLOOR")'
new = '"CONSENSUS_MODE", "CASCADE_MIN_SHARE", "RATE_FLOOR", "TIE_BREAK_ORDER")'
assert s.count(old) == 1, "WANT_CONST 锚点不唯一"
s = s.replace(old, new)

anchor = 'cat(sprintf("[提取] 补丁=%s\\n", SRC))'
assert s.count(anchor) == 1, "锚点不唯一"
inject = (
    "# v6.5/v6.6 补：环境变量派生的常量不在提取范围（Sys.getenv 不是 function 定义）\n"
    "if (!exists(\"VOTE_WEIGHT_MODE\")) {\n"
    "  .vw <- Sys.getenv(\"MMPV_VOTE_WEIGHT\", unset = \"count\")\n"
    "  VOTE_WEIGHT_MODE <- if (.vw %in% c(\"count\", \"weighted\")) .vw else \"count\"\n"
    "}\n"
    "if (!exists(\"TIE_BREAK_ORDER\")) {\n"
    "  TIE_BREAK_ORDER <- c(\"ACVirus\", \"CAT\", \"VITAP\", \"diamond_lca\", \"genomad\", \"metabuli\", \"mmseqs\")\n"
    "}\n"
    "cat(sprintf(\"[环境] VOTE_WEIGHT_MODE=%s | TIE_BREAK_ORDER=%s | TAX_GATE_VERSION=%s\\n\",\n"
    "            VOTE_WEIGHT_MODE, paste(TIE_BREAK_ORDER, collapse = \">\"), TAX_GATE_VERSION))\n\n"
)
s = s.replace(anchor, inject + anchor, 1)

io.open(P, "w", encoding="utf-8", newline="\n").write(s)
print("done, len=%d" % len(s))
