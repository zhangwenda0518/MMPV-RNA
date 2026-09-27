#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 v6.9「平票不猜」变体补丁：靠工具顺序兜底的平票格一律落空。

以 v6.6（当前主线行为）为基。设计要点：
  * 只影响「出值」这一步，绝不碰 win（win 还要供淘汰判定使用，动它会连带改变逐级淘汰本身）
  * 落空按阶元隔离（tie_ord_lvl 每层重置），不跨层累积
  * 逐级淘汰、淘汰计数、平票计数与 v6.6 完全一致；唯一差异是顺序兜底档的格子变空

七处改动：
  1/2. 新增 tie_ord_cells（总计）与 tie_ord_lvl（本层）
  3.   记录本层「下一级一致性没分辨出来、只能靠顺序」的格子
  4.   出值后把这些格子的 Taxon 置 NA
  5.   TAX_GATE_VERSION 6.6 -> 6.9
  6.   版本注补 6.9 行
  7.   日志把「工具顺序兜底」改为「顺序兜底已置空」
用法: python3 make_tieblank_variant.py <源补丁 v6.6> <输出>
"""
import sys

if len(sys.argv) < 3:
    print(__doc__)
    raise SystemExit(2)

SRC, OUT = sys.argv[1], sys.argv[2]
txt = open(SRC, encoding="utf-8", errors="replace").read()

EDITS = [
    (
        "  n_tie <- 0L; n_tie_next <- 0L; n_tie_order <- 0L\n",
        "  n_tie <- 0L; n_tie_next <- 0L; n_tie_order <- 0L\n"
        "  tie_ord_cells <- character(0)   # v6.9: 顺序兜底档的格（全阶元合计，仅用于日志）\n",
    ),
    (
        "    lvl <- TAX_LEVELS[i]\n",
        "    lvl <- TAX_LEVELS[i]\n"
        "    tie_ord_lvl <- character(0)   # v6.9: 本层顺序兜底档（每层重置，不跨层）\n",
    ),
    (
        "      n_tie_order <- n_tie_order + sum(!(dec$pick > 0 & dec$pick > dec$next_pick))\n",
        "      n_tie_order <- n_tie_order + sum(!(dec$pick > 0 & dec$pick > dec$next_pick))\n"
        "      # v6.9: 这一档没有证据，只有排序约定；记下来，出值时落空\n"
        "      tie_ord_lvl <- dec[!(pick > 0 & pick > next_pick), contig_id]\n"
        "      tie_ord_cells <- c(tie_ord_cells, tie_ord_lvl)\n",
    ),
    (
        "    pieces[[lvl]] <- win[, .(contig_id, Rank = lvl, Taxon = win_tax)]\n",
        "    pieces[[lvl]] <- win[, .(contig_id, Rank = lvl, Taxon = win_tax)]\n"
        "    # v6.9: 平票且下一级一致性分不出来 -> 该阶元不猜，置空（宁空不猜）\n"
        "    #       只动 pieces，不动 win（win 还要供下面的淘汰判定用）\n"
        "    if (length(tie_ord_lvl) > 0L) pieces[[lvl]][contig_id %in% tie_ord_lvl, Taxon := NA_character_]\n",
    ),
    (
        'TAX_GATE_VERSION <- "6.6"',
        'TAX_GATE_VERSION <- "6.9"',
    ),
    (
        "# 6.6 = 平票裁决两级化：先用下一（更细）级的一致性分辨，否则按 TIE_BREAK_ORDER 显式兜底；\n",
        "# 6.6 = 平票裁决两级化：先用下一（更细）级的一致性分辨，否则按 TIE_BREAK_ORDER 显式兜底；\n"
        "# 6.9 = 上面那条兜底取消：下一级一致性也分不出来的平票格，该阶元置空（不再用排序约定硬选）\n",
    ),
    (
        '"逐级淘汰: 阶元计票 %d 格次, 淘汰工具投票权 %d 格次 (票首门槛 %.2f); 平票 %d 格次: 下一级一致性定 %d, 工具顺序兜底 %d",',
        '"逐级淘汰: 阶元计票 %d 格次, 淘汰工具投票权 %d 格次 (票首门槛 %.2f); 平票 %d 格次: 下一级一致性定 %d, 顺序兜底已置空 %d",',
    ),
]

for old, new in EDITS:
    n = txt.count(old)
    if n != 1:
        print("!! 命中 %d 次，期望 1 次：%s" % (n, old.strip()[:60]))
        raise SystemExit(3)
    txt = txt.replace(old, new)

with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write(txt)
print("源: %s" % SRC)
print("出: %s" % OUT)
print("改动 %d 处，全部命中一次" % len(EDITS))
