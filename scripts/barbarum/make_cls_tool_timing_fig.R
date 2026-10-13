#!/usr/bin/env Rscript
# make_cls_tool_timing_fig.R —— 分类工具「全库耗时」条形图 **预渲染**（图② 左上角）
#
# 数据：upset_data/cls_tool_timing.tsv
#   来源（服务器）= ~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/
#                   05_Taxonomy/Votus.classed/Votus_resource_usage.tsv
#   （列 = sample tool wall_sec cpu_sec mem_mb status；7 行 = 7 个分类工具）
#
# ⚠⚠ **口径必须写清**（与 IdentFilter 的耗时图**不是一回事**）：
#   · 本图 = **一次性把全部 21,630 条 vOTU 交给该工具分类的总 wall time**（全库一次跑完）；
#   · IdentFilter 的 `tool_timing.tsv` = **单样本平均耗时**（395 样本）；
#   两者不可比，图注里必须声明，避免被误读成"每条 vOTU 花多少秒"。
#
# ⭐⭐ 行序（用户 2026-10-13）：**与图①/图② 右下点阵的 7 行严格同序** =
#   图①/图② 的集合行序由 ComplexUpset 按 **set size 降序**排（`sort_sets` 默认 TRUE）。
#   → 本脚本**不能**硬编码顺序，必须**按同一份矩阵算出各工具命中数**再降序排，
#     否则耗时图与点阵行错位（眼睛上下对读会读错）。
#
# ⚠⚠ 踩过的两个"改行序"静默坑（IdentFilter 同款，务必别重犯）：
#   ① `as.numeric(factor)` 给的是 **levels 序号**，不是轴标签序 → 必须显式 `ypos` 连续坐标；
#   ② `scale_y_continuous` 的 y **向上增大** → `ypos` 要**反相**（`N+1-match`）才能让第 1 名在最上。
#
# 用法：
#   Rscript make_cls_tool_timing_fig.R                    # 默认 → upset_data/cls_tool_timing_fig.png
#   Rscript make_cls_tool_timing_fig.R --w 2.9 --h 1.5    # 指定目标框英寸

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

script_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  getwd()
}
SD <- script_dir()
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d) { i <- which(args == k); if (length(i)) args[i + 1] else d }

out <- get_arg("--out", "")
# 目标框英寸（默认 = 图② 左上实测空区；见下）
# ⚠⚠ 默认尺寸必须与 fig_upset_classification.R 的 TIM_L/TIM_R/TIM_T/TIM_B 框**实测一致**
#   （2026-10-13 用户「时间的图位置往下些，往右边些」→ 框改为 2.3667 × 1.2100 in）。
#   ⚠ 改了宿主框却没重跑本脚本 → **预渲 PNG 不随参数更新**，会静默嵌旧尺寸图。
W_IN <- as.numeric(get_arg("--w", "2.3667"))
H_IN <- as.numeric(get_arg("--h", "1.2100"))
DPI  <- as.numeric(get_arg("--dpi", "600"))
RANK <- get_arg("--rank", "Species")
TOP  <- as.integer(get_arg("--top", "18"))

# ⚠⚠⚠ 默认输出名必须**随 RANK 变化**：两个 rank 的行序不同（Species / Family 末两位互换），
#   若都写 `cls_tool_timing_fig.png` 会**互相覆盖**（实测：`--rank Family` 仍写无后缀名，零报错）。
#   Species → `cls_tool_timing_fig.png`（无后缀，保持向后兼容）；
#   Family  → `cls_tool_timing_fig_Family.png`。
if (!nzchar(out)) {
  out <- file.path(SD, "upset_data",
                   if (identical(RANK, "Species")) "cls_tool_timing_fig.png"
                   else sprintf("cls_tool_timing_fig_%s.png", RANK))
  cat(sprintf("未指定 --out → 按 rank 推导：%s\n", basename(out)))
}

# ── 1. 算图② 点阵的行序（set size 降序）──
mtx <- fread(file.path(SD, "upset_data", sprintf("agree_matrix_%s.tsv", RANK)), sep = "\t")
TOOLS7 <- c("diamond_lca", "mmseqs", "CAT", "metabuli", "genomad", "VITAP", "ACVirus")
tools <- TOOLS7[TOOLS7 %in% names(mtx)]
# set size = 命中该工具的 vOTU 数（与 ComplexUpset 左侧 set-size 面板同口径）
set_size <- vapply(tools, function(t) sum(mtx[[t]]), numeric(1))
ORDER <- names(sort(set_size, decreasing = TRUE))   # 第 1 名在最上（矩阵首行）

# ── 2. 读耗时 ──
tdf <- fread(file.path(SD, "upset_data", "cls_tool_timing.tsv"), sep = "\t")
if (!setequal(tdf$tool, tools)) {
  stop(sprintf("耗时表的工具集合 %s 与矩阵 %s 不一致 → 请复查口径",
               paste(sort(tdf$tool), collapse="/"), paste(sort(tools), collapse="/")))
}
# 行序自检：ORDER 必须与矩阵列名集合一致
if (!setequal(ORDER, tools)) stop("点阵行序与工具集合不一致")

N <- length(ORDER)
tdf <- tdf[match(ORDER, tool)]
# ⚠ ypos 反相：ORDER[1] → N（最上）
tdf[, ypos := N + 1L - match(tool, ORDER)]
stopifnot(all(sort(tdf$ypos) == seq_len(N)))
cat(sprintf("行序（上→下）= %s\n", paste(ORDER, collapse=" > ")))
cat("各工具 set size / 耗时：\n")
print(data.frame(tool = ORDER, set_size = set_size[ORDER], wall_sec = tdf$wall_sec))

# ── 3. 画图 ──
FLOOR_S <- 10                       # x 轴下限（对数轴不能到 0）
lab_at <- setNames(N + 1L - seq_along(ORDER), ORDER)

p <- ggplot(tdf) +
  # ⚠ 用 geom_rect 而非 geom_col：横条从 FLOOR_S 起算，避免对数轴下界截断。
  geom_rect(aes(xmin = FLOOR_S, xmax = wall_sec,
                ymin = ypos - 0.30, ymax = ypos + 0.30),
            fill = "#0072B2") +
  geom_text(aes(x = wall_sec, y = ypos, label = sprintf("%.0f", wall_sec)),
            hjust = -0.18, size = 1.85, colour = "#222222") +
  scale_x_log10(limits = c(FLOOR_S, 1600), expand = c(0, 0),
                breaks = c(10, 30, 100, 300, 1000),
                labels = c("10", "30", "100", "300", "1000"),
                name = "Total runtime, whole library (s, log)") +
  scale_y_continuous(breaks = unname(lab_at), labels = names(lab_at),
                     limits = c(0.4, N + 0.6), expand = c(0, 0)) +
  labs(y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
    axis.text.y  = element_text(size = 6.4),
    axis.text.x  = element_text(size = 6.0),
    axis.title.x = element_text(size = 6.4, margin = margin(t = 1)),
    plot.margin  = margin(1, 8, 1, 1)
  )

ggsave(out, p, width = W_IN, height = H_IN, dpi = DPI, bg = "white")
cat(sprintf("耗时图预渲 -> %s（%d 工具，%.1f–%.1f s，%.4f × %.4f in @%.0f dpi）\n",
            out, N, min(tdf$wall_sec), max(tdf$wall_sec), W_IN, H_IN, DPI))
