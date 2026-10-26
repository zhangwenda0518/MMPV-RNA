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
# fig_host_prediction_b_panel.R — 宿主预测 panel b 独立图（2026-10-14 第 2 版）
#
# 版式 = fig_upset_classification.R 已迭代 7 轮的成熟形态：
#   UpSet 全幅主体 + 两个预渲 inset（wall-time → 主体左上空胞；Venn → 交集柱上方天空区）。
#   ⚠ inset 必须走「ggsave PNG → rasterGrob → patchworkGrob 上叠」：
#     patchwork::inset_element 作用在 patchwork 上会静默误算 npc（分类图实测框缩一半），勿用。
#
# 修复记录（对照首版 .bak）：
#   ① `encode_sets = FALSE`（README 硬规则 18）：缺它时 set-size 面板集合名被编码成
#      "1 4 2 3" 数字乱码，零报错。
#   ② `set_sizes` 显式挂 `scale_y_continuous`：缺它时 ComplexUpset 默认给
#      `scale_y_reverse()` → 柱从右往左倒着长 + y 轴标题是 "group"。
#   ③ 交集柱 `counts = FALSE` + 自绘 annotate：内置柱顶文字会把最大柱的数字压进柱内
#     （分类图踩坑记录）；且 min_size=1 时 8/4/4/3 有并列 → 自绘标注必须去掉并列
#      （min_size = 5，保留 11 个互不相同的大小，脚本内自检兜底）。
#   ④ themes key 必须 === base_annotations 的元素名 "vOTUs"；写 "intersections" 静默失效。
#   ⑤ 消灭左上空胞与柱图天空空白：分别嵌入 wall-time 与 Venn inset（见 FRAMES）。
#
# 数据：host_timing.tsv / host_predictor_matrix.tsv
# 用法：
#   Rscript fig_host_prediction_b_panel.R                 # 成品（PNG+PDF）
#   Rscript fig_host_prediction_b_panel.R --base-only     # 只出无 inset 底图（量坐标用）
#   Rscript fig_host_prediction_b_panel.R --frames "..."  # 临时覆盖两处 inset 框（px @300dpi）

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexUpset)
  library(grid)
  library(ggVennDiagram)
  library(png)
})

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(f) else "."
}
source(file.path(script_dir(), "upset_common.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d) {
  i <- which(args == k)
  if (length(i)) args[i + 1] else d
}
datadir <- get_arg("--data", file.path(script_dir(), "upset_data"))
outdir  <- get_arg("--outdir", "reports/_discovery_sci_20261012/_lt/10_Reports/S10_SCI_Discovery")
# 布尔开关：共享 dest 时先注册者的 default 生效 → 用「只查存在性」的写法
base_only <- ("--base-only" %in% args)
# ⭐ v22.2（2026-10-14）配色开关（与 Q 图同款，用户确认「方案A 身份色」）：
want_pred_color <- !("--plain-setsize" %in% args)
quant_scheme <- get_arg("--quant-scheme", "viridis")

vfmt <- function(v) format(v, big.mark = ",", scientific = FALSE)

# ════════════════════════════════════════════════════════════════════════
# inset ①：宿主预测器耗时（log 水平条形）→ 叠在 UpSet 左上空胞
# ════════════════════════════════════════════════════════════════════════
tt <- fread(file.path(datadir, "host_timing.tsv"), sep = "\t")
tt[, short := fifelse(grepl("ICTV", tool), "ICTV",
  fifelse(grepl("RVH|RNAVirHost", tool), "RVH",
    fifelse(grepl("PB2|PhaBOX", tool), "PB2",
      fifelse(grepl("Rule", tool), "Rule", "Other"))))]
# ⭐ user 2026-10-14：Rule 行整行删掉（无柱空行），"in-process" 信息由副标题带过
tt[, tool := factor(tool, levels = rev(c(
  "RNAVirHost", "PhaBOX2 (Cherry)", "ICTV (C9 lookup)")))]
tt[, short := factor(short, levels = c("RVH", "PB2", "ICTV", "Rule"))]

timing_colors <- HOSTPRED_COLORS[levels(tt$short)]

pb1 <- ggplot(tt, aes(y = tool, x = wall_sec, fill = short)) +
  geom_col(data = tt[!is.na(wall_sec)], width = 0.6, color = "black", linewidth = 0.25) +
  geom_text(data = tt[!is.na(wall_sec)],
    aes(label = sprintf("%s s", format(wall_sec, big.mark = ","))),
    hjust = -0.15, size = 2.4) +
  scale_fill_manual(values = timing_colors, guide = "none") +
  scale_x_log10(expand = expansion(mult = c(0.01, 0.18)),
    breaks = c(1, 10, 100, 500),
    labels = c("1", "10", "100", "500")) +
  labs(
    title = "Predictor wall time",
    subtitle = expression(paste("Single full-set run (", log[10], " scale)")),
    x = "Wall time (s)", y = NULL
  ) +
  theme_upset(8) +
  theme(
    axis.text.y = element_text(size = 7.5, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(2, 4, 2, 4)
  )

# ════════════════════════════════════════════════════════════════════════
# inset ②：四预测器 Venn → 叠在交集柱上方天空区
# ════════════════════════════════════════════════════════════════════════
hp <- fread(file.path(datadir, "host_predictor_matrix.tsv"), sep = "\t")
preds <- c("ICTV", "RVH", "PB2", "Rule")
hp[, (preds) := lapply(.SD, as.logical), .SDcols = preds]

vsets <- lapply(preds, function(p) which(hp[[p]]))
names(vsets) <- preds
venn_union <- length(unique(unlist(vsets)))
cat(sprintf("Venn: union=%s | sets=%s\n",
    format(venn_union, big.mark = ","),
    paste(sprintf("%s %s", preds, format(vapply(vsets, length, numeric(1)), big.mark = ",")),
          collapse = " / ")))

# 用 ggVennDiagram 内部 API 取区域标签坐标
rdat <- as.data.frame(
  ggVennDiagram:::process_data(ggVennDiagram::Venn(vsets), shape_id = "401f")$regionLabel)
VENN_LABEL_MIN <- 50
rlab <- rdat[rdat$count >= VENN_LABEL_MIN, c("X", "Y", "count")]

pb2 <- ggVennDiagram::ggVennDiagram(
    vsets, shape_id = "401f",
    label = "none",
    set_color    = "#333333",
    set_size     = 2.4,
    edge_size    = 0.4) +
  scale_fill_gradient(low = "#F2F2F2", high = "#F2F2F2", guide = "none") +
  geom_text(data = rlab, aes(x = X, y = Y, label = vfmt(count)),
            size = 1.85, colour = "white", fontface = "bold") +
  geom_text(data = rlab, aes(x = X, y = Y, label = vfmt(count)),
            size = 1.85, colour = "#333333") +
  labs(title = "Predictor overlap (Venn)") +
  theme_void(base_size = 8) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0, margin = margin(2,0,2,4)),
    plot.margin = margin(2, 4, 2, 4)
  )

# ════════════════════════════════════════════════════════════════════════
# UpSet 主体（分类图同款修复：encode_sets / set_sizes / 自绘柱顶数字）
# ════════════════════════════════════════════════════════════════════════
hp[, n_pred := rowSums(.SD), .SDcols = preds]
cat(sprintf("UpSet universe: %s / %d\n", format(nrow(hp), big.mark = ","), 20892))

# 柱顶自绘标注要求「显示的交集大小互不相同」（并列时 ComplexUpset 先后与表序可能
# 不一致 → 数字整体错位，零报错；分类图同款自检）。
# 实测分布：… 49 / 8 / 4 / 4 / 3 → min_size=5 恰好去掉并列 4 与尾巴噪声，留 11 根。
sig <- hp[, do.call(paste0, lapply(.SD, as.integer)), .SDcols = preds]
sizes <- sort(table(sig), decreasing = TRUE)
min_size <- 5L
kept <- sizes[sizes >= min_size]
n_show <- length(kept)
if (!n_show || anyDuplicated(as.integer(kept))) {
  stop("显示的交集大小存在并列或为空 → 柱顶自绘标注会错位，请调 min_size")
}
cat(sprintf("显示交集 %d 个（size >= %d）：%s\n", n_show, min_size,
    paste(as.integer(kept), collapse = " / ")))
print(table(hp$n_pred))

pb_body <- upset(
  as.data.frame(hp), preds,
  name = "Host predictor",
  # ⭐ 与分类图第 7 次修订同款组合：矩阵行标签 element_blank（去重复 + 消 205px 槽）
  #   + width_ratio 0.312 补偿 → 点阵面板位置与分类图一致，set-size 面板最宽。
  width_ratio = 0.312,
  min_size = min_size,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  # ⭐ README 硬规则 18：缺 FALSE 时集合名被编码成整数 id，set-size 面板变 "1 4 2 3"。
  encode_sets = FALSE,
  # ⭐ 显式挂 y 标度：压掉默认 scale_y_reverse()（否则柱从右往左倒着长），
  #   并把默认 y 轴标题 "group" 换成语义名；集合名只由最左 set-size 面板显示。
  # ⭐⭐ v22.2（2026-10-14 用户「左侧条形图配色调整」）：set-size 条按**预测器身份**上色。
  #   ⚠⚠⚠ 只能用「颜色向量」：`aes(fill = ...)` 三种写法在 ComplexUpset 1.3.3 全报错
  #     （`object 'x'/'set' not found`、`Cannot add <ggproto> objects together`，见
  #      `upset_common.R::set_size_colors()` 注释里的实测清单）。
  #   ✅ `geom_bar(fill = <向量>)` 生效；向量按**传入 `preds` 顺序**逐条对应，
  #      **与渲染顺序无关**（Q 图默认 descending 下已实测确认）。
  #   ⚠ 宽度保持 0.6（默认）—— Q 图实测写成 0.62 会带来 3,855 px 差异。
  #   ⚠ 色源用 `HOSTPRED_COLORS`：与**左上耗时 inset 同色**（同一预测器同色），
  #      这是选「身份色」方案的核心收益。
  set_sizes = upset_set_size(
      geom = geom_bar(
        fill = if (want_pred_color) set_size_colors(preds, HOSTPRED_COLORS) else "grey35",
        width = 0.6)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    xlab("Set size") +
    theme_upset(8) +
    theme(axis.title.y = element_blank()),
  base_annotations = list(
    "vOTUs" = intersection_size(
      # ⭐ counts=FALSE：内置柱顶文字会把最大柱数字压进柱内（分类图踩坑）；改自绘。
      counts = FALSE,
      mapping = aes(fill = as.factor(n_pred))
    ) +
    annotate("text",
             x = seq_len(n_show),
             y = as.integer(kept),
             label = vfmt(as.integer(kept)),
             vjust = -0.5, size = 2.4, colour = "#222222") +
    scale_fill_manual(
      # ⭐ v22.2：改用**三张图统一的数量色阶**（`upset_common.R::QUANT_COLORS`）。
      #   ⚠ 旧口径（1=#F0E442 / 2=#56B4E9 / 3=#0072B2 / 4=#D55E00）与左侧预测器
      #     身份色撞了 2 个完全同色（#0072B2=ICTV、#D55E00=RVH）→ 故换掉。
      #   ⚠ 档位取 1..4，与 Q 图的 1..7 共用同一色表 → **同值同色**（跨图一致）。
      #   `--quant-scheme degree` 回到旧四色（零回归对拍用）。
      values = if (quant_scheme == "degree")
                 c("1" = "#F0E442", "2" = "#56B4E9", "3" = "#0072B2", "4" = "#D55E00")
               else if (quant_scheme == "blues")
                 QUANT_BLUES[1:4]
               else stats::setNames(quant_colors(1:4), as.character(1:4)),
      name = "Predictors", drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    ylab("vOTUs") +
    theme_upset(8) +
    theme(
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
  ),
  stripes = upset_stripes(mapping = aes(color = "black", fill = "grey92")),
  matrix = intersection_matrix(geom = geom_point(shape = 21, size = 2.8, stroke = 0.8)) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey88"), guide = "none"),
  themes = upset_modify_themes(list(
    # ⚠⚠ key 必须 === base_annotations 元素名 "vOTUs"；写 "intersections" 静默失效。
    "vOTUs" = theme_upset(8) +
      theme(
        axis.text.y = element_text(size = 7.5, color = "black"),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(4, 4, 2, 4)
      ),
    "overall_sizes" = theme_upset(8) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "black"),
        panel.grid.major.x = element_blank()
      ),
    # 拼回内置两段再叠加（直接整体替换会丢内置 blank → 底部冒出挤团的轴文字，实测坑）
    "intersections_matrix" = Reduce(`+`, ComplexUpset::upset_themes$intersections_matrix) +
      theme(axis.text.y = element_blank())
  ))
)

p <- pb_body +
  plot_annotation(
    title = "b  Host-predictor evidence",
    subtitle = paste0(
      "Wall time (left inset), Venn overlap (top-right inset), and UpSet intersections - ",
      "n = 20,443 vOTUs with >= 1 predictor hit; Rule = taxonomy-class fallback"
    ),
    theme = theme(
      plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.8, color = "#444444", hjust = 0, lineheight = 1.2),
      plot.margin = margin(6, 8, 6, 8)
    )
  )

# ════════════════════════════════════════════════════════════════════════
# 输出：底图（--base-only 量坐标用）/ 成品（两 inset 按 FRAMES 叠加）
# ════════════════════════════════════════════════════════════════════════
W_IN <- 9.32; H_IN <- 6.0          # 画布（与分类图一致）
W_PX <- round(W_IN * 300); H_PX <- round(H_IN * 300)

# 两处 inset 框（px @300dpi，画布 2796×1800）：
#   ⚠ 必须用「当前版本 --base-only 底图」重测，禁沿用旧版坐标（README 硬规则）。
#   实测工具：reports/virome_analysis_pipeline/_probe_bpanel_frames.py（2026-10-14，本版底图）。
#   ⭐ user 2026-10-14 13:xx（红/蓝框标注）：**wall time → 左空胞（红框）**、
#      **Venn → 交集柱右上天空区（蓝框）**。此前两框内容相反，此处已对调。
#   框① wall-time（左空胞）：x [30, 870]（右缘避开柱面板 y 轴刻度文字 x~895 起）、
#      y [300, 1100]（上界在副标题下、下界在 set-size 面板顶 y1142 之上）。
#   框② Venn（右上天空区）：x [1600, 2400]（左界在次高柱 4190 顶数字右侧净空外）、
#      y [210, 750]（上界在题头下、下界在次高柱顶 765 之上）。
#   落点自检：两框内底图墨迹均 = 0 px（_check 见 README 交付纪律）。
FR_WALLTIME <- c(30, 300, 870, 1100)     # L, T, R, B（红框＝左空胞 840x800px，2.80×2.67 in）
FR_VENN     <- c(1600, 210, 2400, 750)   # （蓝框＝交集柱上空区 800x540px，2.67×1.80 in）

frames <- get_arg("--frames", "")
if (nzchar(frames)) {
  v <- as.integer(strsplit(frames, "[;, ]+")[[1]])
  FR_WALLTIME <- v[1:4]; FR_VENN <- v[5:8]
}

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
base_png <- file.path(outdir, "Fig_GojiVirome_HostPrediction_BPanel_base.png")

if (base_only) {
  ggsave(base_png, p, width = W_IN, height = H_IN, dpi = 300, bg = "white")
  cat("Base only:", base_png, "\n")
  quit(save = "no", status = 0)
}

# inset 框 → npc
mk_ov <- function(fr, grob) {
  list(grob = grob,
       cx = (fr[1] + fr[3]) / 2 / W_PX, cy = 1 - (fr[2] + fr[4]) / 2 / H_PX,
       w = (fr[3] - fr[1]) / W_PX,      h = (fr[4] - fr[2]) / H_PX)
}
mk_grob <- function(plot, fr) {
  w_in <- (fr[3] - fr[1]) / 300; h_in <- (fr[4] - fr[2]) / 300
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, plot, width = w_in, height = h_in, dpi = 300, bg = "white")
  grid::rasterGrob(png::readPNG(tmp), interpolate = TRUE,
    width = unit(w_in, "in"), height = unit(h_in, "in"))
}

overlays <- list(
  mk_ov(FR_WALLTIME, mk_grob(pb1, FR_WALLTIME)),
  mk_ov(FR_VENN,     mk_grob(pb2, FR_VENN))
)

# ⭐ 不用 inset_element：作用在 patchwork 上会静默误算 npc（分类图实测，框缩一半）。
#   唯一零 reflow 画法 = patchworkGrob 整图绘制后按 npc 叠 rasterGrob。
one <- function(path, dev_open) {
  dev_open(path)
  grid::grid.newpage()
  grid::grid.draw(patchwork:::patchworkGrob(p))
  for (ov in overlays) {
    grid::pushViewport(grid::viewport(x = ov$cx, y = ov$cy,
                                      width = ov$w, height = ov$h, just = "center"))
    grid::grid.draw(ov$grob)
    grid::popViewport()
  }
  grDevices::dev.off()
  invisible(path)
}
name <- "Fig_GojiVirome_HostPrediction_BPanel"
one(file.path(outdir, paste0(name, ".png")),
    function(f) grDevices::png(f, width = W_IN, height = H_IN, units = "in", res = 300, bg = "white"))
one(file.path(outdir, paste0(name, ".pdf")),
    function(f) grDevices::pdf(f, width = W_IN, height = H_IN, bg = "white"))
cat("Written:", file.path(outdir, paste0(name, ".png")), "\n")
