#!/usr/bin/env Rscript
# ── 图件四格式导出 (SVG / PNG / PDF / PNG@1200dpi) — 默认开启 ──────────
# 2026-10-20: 所有出图脚本统一口径; MMPV_FIG_EXPORT=0 可整体关闭。
# 钩子只拦截 ggsave(); pdf()/png() 直写的老代码请改用 save_figure()。
local({
  .find_repo <- function(start) {
    d <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(d, "mmpv_common"))) return(d)
      p <- dirname(d)
      if (identical(p, d)) return(NA_character_)
      d <- p
    }
  }
  .here <- tryCatch({
    a <- commandArgs(FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
  }, error = function(e) getwd())
  .root <- .find_repo(.here)
  if (is.na(.root)) .root <- .find_repo(getwd())
  if (!is.na(.root)) {
    source(file.path(.root, "mmpv_common", "figure_export.R"), local = FALSE)
    install_ggsave_hook()
  } else {
    message("[figure_export] 未找到 mmpv_common/, 跳过四格式导出")
  }
})
# ──────────────────────────────────────────────────────────────────────
# make_tool_timing_fig.R —— 图⑧ 左上格「7 轨道单样本平均耗时」条形图 **预渲染**
#
# ⭐⭐ 为什么必须独立成脚本、独立进程渲染：
#   在 fig_ident_filter_stacked_all.R 内直接 `ggsave(timing_plot)` 会给整张主图引入
#   **18,173 px 的子像素重栅格噪声**（set-size 面板 / 右下点阵 / 交集柱全部变），
#   虽然 Δink 仅 −273 px（= 纯字形相位差，不是版式位移），但**破坏了「零回归」判据**。
#   逐项排除结果（详见 fig_ident_filter_stacked_all.R 内的注释）：
#     ✗ 换 device（png()/ggsave）  ✗ 加 dev.off()  ✗ 换渲染尺寸
#     ✓ 把 timing_plot 换成 theme_void() 玩具图 → 差异归 0
#   → 根因 = `theme_minimal` 图（含 axis.text / panel.grid）的 `ggplot_build()` 会改写
#     **进程级字形度量缓存**，而 ComplexUpset 的标签/格子尺寸依赖同一缓存。
#   → 故：本脚本**单独**出 PNG，主脚本只 `readPNG` 贴图，主进程**绝不构建**该 ggplot。
#
# 数据：identfilter_data/tool_timing.tsv（由 make_tool_timing_data.py 生成）
# 尺寸：与主脚本目标框一致（845 × 590 px @300dpi = 2.8167 × 1.9667 in），dpi = 600（2× 过采样）
# 用法：
#   Rscript make_tool_timing_fig.R                 # 用默认路径
#   Rscript make_tool_timing_fig.R --out <png> --tsv <tsv>
#
# ⭐⭐ 行序（用户 2026-10-13 第 8 次修订）：**必须与主图右下点阵的 7 行严格同序** =
#   `TRACKS` 顺序（viralverify / metabuli / blast / genomad / virbot / rdrpcatch / viroid_blastn，
#   即按 union 内命中数列降序）。理由：两图上下对读时，同一轨道在同一水平线上，眼睛不用来回找。
#   ⚠ 旧的「按平均耗时降序」（genomad 最上）会让两图错位 → 已改。
#   ⚠ 但**条长 = 耗时**，所以本图不再是单调递减 —— 这正是「命中多 ≠ 耗时高」的可见证据
#     （如 blast 命中 95,573 却只花 9.2 s；genomad 命中 11,526 却花 235.2 s = 全图最贵）。
#
# ⚠⚠ 本轮（第 8 次修订）实测踩到的**两个"改行序"专属静默坑**，务必别重犯：
#   ① **`as.numeric(factor)` 定位 ≠ 轴标签顺序**。首版把行序改成 `factor(levels=rev(...))`
#      后仍沿用 `y = as.numeric(tool)` 画条，而 `scale_y_continuous(breaks=1:7,
#      labels=TRACKS_ORDER)` 是按**原序**贴标签 → 两套方向相反 → **条长与标签错位整行**
#      （实测图面显示 viralverify = 7.5 s，真值 30.0 s），**零报错、零警告**，
#      必须逐行对着 tsv 核才看得出。→ 正解：显式 `ypos` 连续坐标，条与标签共用同一 ypos。
#   ② **`identical(sort(int), as.numeric(1:7))` 恒 FALSE** —— 整数 vs 双精度类型不同。
#      守卫要用 `all(a == b)`。（首版此守卫直接 `stop()`，反倒帮忙拦住了①的错位图。）
#   ③ 结论纪律：**改行序（y 轴类）改动后必须渲染出来逐行核对「标签 ↔ 条长」**，
#      不能只看「7 行都在、顺序对」就收工 —— 错位正是发生在"顺序对、配对错"的时候。

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

script_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  getwd()
}
SD <- script_dir()
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d) { i <- which(args == k); if (length(i)) args[i + 1] else d }

out    <- get_arg("--out",    file.path(SD, "identfilter_data", "tool_timing_fig.png"))
# ⭐ 目标框 = 主脚本 `fig_ident_filter_stacked_all.R` 内以画布 px 实测的 TIM 框。
#   ⚠⚠ **第 13 次修订（用户 2026-10-13 21:55 红框「再次调整大小，符合红色区域的大小」）**：
#      耗时图由 `x50..1162 / y424..1381`（3.7067 × 3.1900 in，比 1.162）
#      → **`x51..1156 / y545..1494`（3.6833 × 3.1633 in，比 1.164）**。
#      ⚠⚠⚠ 第 12 次修订的红框→画布映射 **Y 偏移错 ≈110 px**（矩阵 7 行棕点标定实测）
#      → 图整体偏高 121 px。本版改用矩阵点阵锚点（残差 0 px）后**整块下移 121 px**。
#   ⚠⚠ **框比变了源图必须重渲**（预渲 PNG 不随参数更新；且比不符会被拉伸变形）：
#      `Rscript make_tool_timing_fig.R --w 3.6833 --h 3.1633`（或直接改下面默认值）。
#   ⚠ 默认值必须与 TIM 框一致：wr = 1105/300 = 3.6833、hr = 949/300 = 3.1633。
W_IN   <- as.numeric(get_arg("--w", "3.6833"))
H_IN   <- as.numeric(get_arg("--h", "3.1633"))
DPI    <- as.numeric(get_arg("--dpi", "600"))

BASE <- get_arg("--base", file.path(SD, "identfilter_data"))
tdf <- data.table::fread(file.path(BASE, "tool_timing.tsv"))

# ⭐⭐ 行序 = 与主图点阵同序（见文件头）。顺序写死在下方 TRACKS_ORDER，
#   并**自检**它与 tsv 里的 tool 集合完全一致（少一个/多一个都报错，避免静默漏轨道）。
TRACKS_ORDER <- c("viralverify", "metabuli", "blast", "genomad",
                  "virbot", "rdrpcatch", "viroid_blastn")
if (!setequal(tdf$tool, TRACKS_ORDER)) {
  stop(sprintf("tsv 的 tool 集合 %s 与主图点阵行序 %s 不一致 → 请复查口径",
               paste(sort(tdf$tool), collapse = "/"), paste(TRACKS_ORDER, collapse = "/")))
}
#> 注意：`factor(levels = )` 在 ggplot 的**离散 y 轴**上是从**下往上**排的，
#  即以 levels[1] 为最底。而主图点阵的自上而下顺序是 viralverify→…→viroid_blastn，
#  故这里 levels 取 **rev**() → viralverify 落在最上，与点阵首行对齐。
tdf[, tool := factor(tool, levels = rev(TRACKS_ORDER))]
if (anyNA(tdf$tool)) stop("tool 因子化后出现 NA → TRACKS_ORDER 与 tsv 不匹配")
# ⚠⚠⚠ **不要用 `as.numeric(tool)` 定位** —— 实测踩坑：`as.numeric(factor)` 给的是
#   **levels 顺序的序号**（此处置为 rev 序），而 `scale_y_continuous(breaks=1:7,
#   labels=TRACKS_ORDER)` 又是按**原序**配 break → 两套方向相反，
#   结果「标签贴 A、条长却是 B」，实测整表错位一行（viralverify 显示成 7.5 s 的真值
#   实为 viroid_blastn）。**零报错**、肉眼要逐行对才看得出。
#   → 正解：把**画布纵向位置 ypos 显式写死**（1 = 最下），条与标签共用同一个 ypos，
#     并把 y 轴改成**连续轴 + breaks = ypos 的整数**，彻底断开对 factor 序号的依赖。
# ⚠⚠⚠ **`ypos` 的取向陷阱**（实测踩坑第 2 次）：`scale_y_continuous` 的 y **向上增大**，
#   故 ypos = 7 在**最上**、ypos = 1 在**最下**。若直接 `ypos := match(tool, ORDER)`
#   （viralverify→1），viralverify 会落在**最下**、viroid_blastn 落在**最上** —— 正好上下颠倒
#   （实测就是如此：打印说"上→下 = viralverify > …"，图上却是 viroid_blastn 在最上）。
#   → 正解：**反相** `ypos := length(ORDER) + 1 - match(...)`，让 ORDER[1] 落在最上。
N_TRACKS <- length(TRACKS_ORDER)
tdf[, ypos := N_TRACKS + 1L - match(tool, TRACKS_ORDER)]  # ORDER[1] → 7（最上）
# ⚠ `identical()` 对 1L:7L（整数）vs 1:7（双精度）返回 FALSE（类型不同）→
#   这里用 `all(... == ...)` 比数值，避免无谓的类型陷阱。
stopifnot(all(sort(tdf$ypos) == seq_along(TRACKS_ORDER)))  # 1..7 各一次
# 自检：最上那条（ypos = N）必须等于 TRACKS_ORDER[1]（与点阵首行同轨道）。
top_tool <- as.character(tdf$tool[which.max(tdf$ypos)])
if (top_tool != TRACKS_ORDER[1]) {
  stop(sprintf("行序自检失败：最上那条是 %s，应为 %s（TRACKS_ORDER[1]）",
               top_tool, TRACKS_ORDER[1]))
}
cat(sprintf("行序自检 OK：最上 = %s（= TRACKS_ORDER[1]）\n", top_tool))
# 轴标签：ypos 的整数位 → 工具名。ypos 已反相，故 break N 对 ORDER[1]。
lab_at <- setNames(N_TRACKS + 1L - seq_along(TRACKS_ORDER), TRACKS_ORDER)
FLOOR_S <- 3
tdf[, xmin := FLOOR_S]

p <- ggplot(tdf) +
  geom_rect(aes(xmin = xmin, xmax = mean_sec,
                ymin = ypos - 0.31, ymax = ypos + 0.31),
            fill = "#0072B2") +
  geom_text(aes(x = mean_sec, y = ypos,
                label = sprintf("%.1f", mean_sec)),
            hjust = -0.16, size = 1.85, colour = "#222222") +
  scale_x_log10(limits = c(FLOOR_S, 900), expand = c(0, 0),
                breaks = c(3, 10, 30, 100, 300),
                labels = c("3", "10", "30", "100", "300"),
                name = "Mean runtime per sample (s, log scale)") +
  # ⚠ y 轴用**连续轴**：breaks = ypos 的整数位、labels 由 `lab_at`（name→ypos 反相映射）。
  #   ypos 已反相 → break N 配 ORDER[1]，与"最上"语义一致。
  scale_y_continuous(breaks = unname(lab_at), labels = names(lab_at),
                     limits = c(0.5, N_TRACKS + 0.5),
                     expand = c(0, 0)) +
  labs(y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
    axis.text.y  = element_text(size = 6.6),
    axis.text.x  = element_text(size = 6.2),
    axis.title.x = element_text(size = 6.6, margin = margin(t = 1)),
    plot.margin  = margin(1, 7, 1, 1)
  )

ggsave(out, p, width = W_IN, height = H_IN, dpi = DPI, bg = "white")
# ⚠ 打印时按**画布自上而下**的顺序（= rev(TRACKS_ORDER)），便于与主图点阵目视逐行核对。
cat(sprintf("耗时图预渲 → %s（%d 轨道，%.1f–%.1f s，%.4f × %.4f in @%.0f dpi）\n",
            out, nrow(tdf), min(tdf$mean_sec), max(tdf$mean_sec), W_IN, H_IN, DPI))
cat(sprintf("  行序（上→下）= %s\n", paste(TRACKS_ORDER, collapse = " > ")))
