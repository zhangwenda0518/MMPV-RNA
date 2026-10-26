#!/usr/bin/env Rscript
# fig_host_prediction.R — 图③ Host Prediction 图（补充 runtime + Venn 版，v2）
#
# 数据：06_HostPrediction/。
#   面板 a  Host distribution（host_distribution.tsv）
#   面板 b\u2081 Predictor wall time（host_timing.tsv）
#   面板 b\u2082 Predictor overlap Venn（host_predictor_matrix.tsv）
#   面板 b\u2083 Predictor overlap UpSet（host_predictor_matrix.tsv）
#   面板 c  Decision_Method 构成（host_decision_method.tsv）
#   面板 d  ICTV 置信度（host_confidence.tsv）
#
# 口径：面板 b\u2083 的 universe = ensemble_host_summary 全表 20,892 条；
#       四预测器任一命中者入图（20,443 条），四路全空者不入 UpSet/Venn。
#       \u26a0 RuleClassTaxonomy 是"分类级回退"，非独立序列比对。
#
# 修订（v2, 2026-10-14）：
#   \u2460 面板 b 拆为 b\u2081/b\u2082/b\u2083：参考 Fig_GojiVirome_IdentFilterStackedAll_norm 布局，
#     b\u2081=\u8017\u65f6（log 水平条形）、b\u2082=\u56db\u96c6 Venn、b\u2083=\u539f UpSet。
#   \u2461 布局改为左列 a/c/d + 右列 b\u2081/b\u2082/b\u2083（三行 \u00d7 两列），H 从 9.4 \u2192 10.5 in。
#   \u2462 Venn 用 ggVennDiagram::shape_id="601" + 平填充 #F2F2F2 + 白描边计数标签；
#     \u26a0 区域计数 <50 的 6 个区不标数字，免重叠（阈值在脚本中显式设）。
#   \u2463 耗时来源：host.log（2026-08-12 单次全量跑），RVH=139+214=353 s；
#     Rule 为 in-process \u2192 wall_sec=NA，图上以文字注记、不画零高柱。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexUpset)
  library(grid)
  library(ggVennDiagram)
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
outdir <- get_arg("--outdir", "reports/_discovery_sci_20261012/_lt/10_Reports/S10_SCI_Discovery")

# ---- 辅助：Venn 计数格式化（千分位、无科学计数） ----
vfmt <- function(v) format(v, big.mark = ",", scientific = FALSE)

# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# 面板 a：宿主分布（完全采用管线 host_distribution.tsv 的值）
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
hd <- fread(file.path(datadir, "host_distribution.tsv"), sep = "\t")
hd[, Host := factor(Host, levels = Host[order(-Count)])]
HOSTCOL <- setNames(OKABE_ITO[c(1, 2, 3, 4, 5, 6, 7, 8, 9)], levels(hd$Host))
pa <- ggplot(hd, aes(x = Host, y = Count, fill = Host)) +
  geom_col(width = 0.74, color = "black", linewidth = 0.25) +
  geom_text(aes(label = format(Count, big.mark = ",")), vjust = -0.45, size = 2.5) +
  scale_fill_manual(values = HOSTCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = "a  Host prediction distribution",
    subtitle = sprintf("Ensemble host call per vOTU (n = %s)", format(sum(hd$Count), big.mark = ",")),
    x = NULL, y = "vOTUs"
  ) +
  theme_upset(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, color = "black"))

# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# 面板 c：Decision_Method 构成
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
dm <- fread(file.path(datadir, "host_decision_method.tsv"), sep = "\t")
setnames(dm, c("Method", "Count"))
dm[, src := fifelse(grepl("ICTV", Method), "ICTV",
  fifelse(grepl("RVH", Method), "RVH",
    fifelse(grepl("PB2", Method), "PB2",
      fifelse(grepl("Rule", Method), "Rule", "None"))))]
dm[, Method := factor(Method, levels = Method[order(-Count)])]
dm[, MethodLab := gsub("_", "\n", as.character(Method))]
dm[, MethodLab := factor(MethodLab, levels = MethodLab[order(-Count)])]
SRCCOL <- c("ICTV" = "#0072B2", "RVH" = "#D55E00", "PB2" = "#009E73",
            "Rule" = "#E69F00", "None" = "#999999")
pc <- ggplot(dm, aes(x = MethodLab, y = Count, fill = src)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(aes(label = format(Count, big.mark = ",")), vjust = -0.45, size = 2.4) +
  scale_fill_manual(values = SRCCOL, name = "Winning source") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title = "c  Decision method",
    subtitle = sprintf("Which evidence resolved the call (total %s)", format(sum(dm$Count), big.mark = ",")),
    x = NULL, y = "vOTUs"
  ) +
  theme_upset(8) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5.4, color = "black", lineheight = 1.05),
    legend.position = "bottom", legend.key.size = unit(0.28, "cm")
  )

# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# 面板 d：ICTV 置信度
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
cf <- fread(file.path(datadir, "host_confidence.tsv"), sep = "\t")
setnames(cf, c("Level", "Count", "Percent"))
cf[, Level := factor(Level, levels = Level[order(-Count)])]
CFCOL <- c("#0072B2", "#E69F00", "#D55E00", "#CC79A7", "#999999")
pd <- ggplot(cf, aes(x = Level, y = Count, fill = Level)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%s\n(%s)", format(Count, big.mark = ","), Percent)),
    vjust = -0.35, size = 2.3) +
  scale_fill_manual(values = CFCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(
    title = "d  ICTV host-call confidence",
    subtitle = "Confidence level of the ICTV taxonomy lookup",
    x = NULL, y = "vOTUs"
  ) +
  theme_upset(8) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 6, color = "black"))

# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# 面板 b\u2081：宿主预测器耗时（log 水平条形）
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
tt <- fread(file.path(datadir, "host_timing.tsv"), sep = "\t")
# 提取短名（ICTV / RVH / PB2 / Rule）用于配色
tt[, short := fifelse(grepl("ICTV", tool), "ICTV",
  fifelse(grepl("RVH|RNAVirHost", tool), "RVH",
    fifelse(grepl("PB2|PhaBOX", tool), "PB2",
      fifelse(grepl("Rule", tool), "Rule", "Other"))))]
# 顺序 = 耗时降序；⭐ user 2026-10-14：Rule 行整行删掉（无柱空行），信息走副标题
tt[, tool := factor(tool, levels = rev(c(
  "RNAVirHost", "PhaBOX2 (Cherry)", "ICTV (C9 lookup)")))]
tt[, short := factor(short, levels = c("RVH", "PB2", "ICTV", "Rule"))]

# 配色映射
timing_colors <- HOSTPRED_COLORS[levels(tt$short)]

pb1 <- ggplot(tt, aes(y = tool, x = wall_sec, fill = short)) +
  geom_col(data = tt[!is.na(wall_sec)], width = 0.6, color = "black", linewidth = 0.25) +
  geom_text(data = tt[!is.na(wall_sec)],
    aes(label = sprintf("%s s", format(wall_sec, big.mark = ","))),
    hjust = -0.15, size = 2.4) +
  scale_fill_manual(values = timing_colors, guide = "none") +
  # ⚠ 右侧 expand 0.24：给 "353 s"/"322 s" 柱端数字留位（0.18 时会被裁掉）
  scale_x_log10(expand = expansion(mult = c(0.01, 0.24)),
    breaks = c(1, 10, 100, 500),
    labels = c("1", "10", "100", "500")) +
  labs(
    title = "Predictor wall time",
    # ⚠ 副标题不能长：b 列只有半幅宽，长副标题会溢出压到右侧 Venn 的标题（实测）。
    # ⚠ 副标题不能长：b 列只有半幅宽，长副标题会溢出压到右侧 Venn 的标题（实测）。
    #   "Rule in-process" 信息已由总副标题覆盖，这里不重复。
    subtitle = expression(paste("Single full-set run (", log[10], " scale)")),
    x = "Wall time (s)", y = NULL
  ) +
  theme_upset(8) +
  theme(
    axis.text.y = element_text(size = 7.5, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# 面板 b\u2082：四预测器 Venn（ggVennDiagram, shape_id="601", 平填充）
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
hp <- fread(file.path(datadir, "host_predictor_matrix.tsv"), sep = "\t")
preds <- c("ICTV", "RVH", "PB2", "Rule")
hp[, (preds) := lapply(.SD, as.logical), .SDcols = preds]

# 按预测器名取 contig index 集合
vsets <- lapply(preds, function(p) which(hp[[p]]))
names(vsets) <- preds
venn_n <- vapply(vsets, length, numeric(1))
venn_union <- length(unique(unlist(vsets)))
cat(sprintf("Venn \u53e3\u5f84\uff1a\u56db\u9884\u6d4b\u5668\u4efb\u4e00\u547d\u4e2d %s\uff5c\u5404\u96c6\u5408 %s\n",
    format(venn_union, big.mark = ","),
    paste(sprintf("%s %s", preds, format(venn_n, big.mark = ",")), collapse = " / ")))

# 用 ggVennDiagram 内部 API 取区域标签坐标
rdat <- as.data.frame(
  ggVennDiagram:::process_data(ggVennDiagram::Venn(vsets), shape_id = "401f")$regionLabel)
# 阈值：只标 count ≥100 的区域。⚠ 半幅宽度下 Venn 很小，中心 82/260/220 三标签会
# 相互叠压（实测）→ 主图抬到 100 舍弃 82；独立 b 面板全幅用 50。
VENN_LABEL_MIN <- 100
rlab <- rdat[rdat$count >= VENN_LABEL_MIN, c("X", "Y", "count")]
cat(sprintf("Venn \u533a\u57df\uff1a%d \u4e2a\uff08\u6807 %d \u4e2a\uff0c\u9608 \u2265%d\uff1b\u5176\u4f59 %d \u4e2a\u4e0d\u6807\uff09\n",
    nrow(rdat), nrow(rlab), VENN_LABEL_MIN, nrow(rdat) - nrow(rlab)))

# 设置色：HOSTPRED_COLORS 顺序 \u2192 与 Venn \u692d\u5706\u6807\u7b7e\u8272\u5bf9\u5e94
tool_venn_colors <- HOSTPRED_COLORS[preds]

pb2 <- ggVennDiagram::ggVennDiagram(
    vsets, shape_id = "401f",
    label = "none",
    set_color    = "#333333",
    set_size     = 2.4,
    edge_size    = 0.4) +
  # 平填充 #F2F2F2（禁渐变，免撞色；见 MEMORY 硬规则）
  scale_fill_gradient(low = "#F2F2F2", high = "#F2F2F2", guide = "none") +
  # 白描边 halo + 正文标签（半幅下用 1.7 号）
  geom_text(data = rlab, aes(x = X, y = Y, label = vfmt(count)),
            size = 1.7, colour = "white", fontface = "bold") +
  geom_text(data = rlab, aes(x = X, y = Y, label = vfmt(count)),
            size = 1.7, colour = "#333333") +
  labs(
    title = "Predictor overlap (Venn)",
    subtitle = sprintf("ICTV / RVH / PB2 / Rule - n = %s with >= 1 hit",
      format(venn_union, big.mark = ","))
  ) +
  theme_void(base_size = 8) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0, margin = margin(4,0,2,4)),
    plot.subtitle = element_text(size = 6.8, color = "#444444", hjust = 0, margin = margin(0,0,2,4)),
    plot.margin = margin(4, 4, 4, 4)
  )

# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# 面板 b\u2083：预测器 UpSet（原面板 b 的主体）
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
hp[, n_pred := rowSums(.SD), .SDcols = preds]
cat(sprintf("\u5bbf\u4e3b\u9884\u9884\u6d4b\u5668\u547d\u4e2d %s / \u5168\u8868 %d\n", format(nrow(hp), big.mark = ","), 20892))
print(table(hp$n_pred))

# ── 柱顶自绘标注（同 fig_host_prediction_b_panel.R 的修复 ③）：counts=TRUE 的内置
#    文字会把最大柱数字压进柱内（分类图踩坑）；且需「显示大小互不相同」，
#    min_size=5 去掉并列 4 与尾巴噪声（8/4/4/3），留 12 根。
sig_b <- hp[, do.call(paste0, lapply(.SD, as.integer)), .SDcols = preds]
sizes_b <- sort(table(sig_b), decreasing = TRUE)
min_size_b <- 5L
kept_b <- sizes_b[sizes_b >= min_size_b]
n_show_b <- length(kept_b)
if (!n_show_b || anyDuplicated(as.integer(kept_b))) {
  stop("b 面板显示的交集大小存在并列或为空 → 柱顶自绘标注会错位，请调 min_size_b")
}
cat(sprintf("b 面板显示交集 %d 个（size >= %d）\n", n_show_b, min_size_b))

pb_body <- upset(
  as.data.frame(hp), preds,
  name = "Host predictor",
  # ⭐ 半幅宽度下左列收窄到 0.24（独立 b 面板全幅用 0.312）：左上空胞过大会留白突兀。
  width_ratio = 0.24,
  min_size = min_size_b,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  # ⭐ README 硬规则 18：缺 FALSE 时集合名被编码成整数 id（set-size 面板变 "1 4 2 3"）。
  encode_sets = FALSE,
  # ⭐ 显式挂 y 标度压掉默认 scale_y_reverse()（否则柱从右往左倒着长 + "group" 标题）。
  set_sizes = upset_set_size() +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    xlab("Set size") +
    theme_upset(8) +
    theme(axis.title.y = element_blank()),
  base_annotations = list(
    "vOTUs" = intersection_size(
      counts = FALSE,
      mapping = aes(fill = as.factor(n_pred))
    ) +
    # ⭐ 半幅宽度柱距仅 ~56px：数字用 1.9 号才不相互打架（全幅独立面板用 2.4）。
    annotate("text",
             x = seq_len(n_show_b),
             y = as.integer(kept_b),
             label = format(as.integer(kept_b), big.mark = ",", scientific = FALSE),
             vjust = -0.5, size = 1.9, colour = "#222222") +
    scale_fill_manual(
      values = c("1" = "#F0E442", "2" = "#56B4E9", "3" = "#0072B2", "4" = "#D55E00"),
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
    # 拼回内置两段再叠加（整体替换会丢内置 blank → 底部冒出挤团的轴文字，实测坑）；
    # 行标签置 blank：集合名只由最左 set-size 面板显示一次。
    "intersections_matrix" = Reduce(`+`, ComplexUpset::upset_themes$intersections_matrix) +
      theme(axis.text.y = element_blank())
  ))
)

# 题头条（ComplexUpset 返回 patchwork，labs 会错位）
title_grob <- function(title, subtitle) {
  gTree(children = gList(
    textGrob(title, x = 0, hjust = 0, y = 0.66, gp = gpar(fontsize = 9.5)),
    textGrob(subtitle, x = 0, hjust = 0, y = 0.24,
      gp = gpar(fontsize = 6.6, col = "#444444"))
  ))
}


# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# 组装（左列 a/c/d，右列 b\u2081/b\u2082/b\u2083）
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# b 面板内部：左上耗时 + 右上 Venn + 下方 UpSet（Venn 在半幅下需要更宽才不挤）
b_top  <- pb1 + pb2 + plot_layout(widths = c(0.42, 0.58))
b_body <- b_top / pb_body + plot_layout(heights = c(0.55, 1.0))

b_panel <- wrap_elements(title_grob(
  "b  Host-predictor evidence",
  sprintf("Wall time, Venn overlap, and UpSet intersections - n = %s with >= 1 hit",
    format(nrow(hp), big.mark = ","))
)) / b_body +
  plot_layout(heights = c(0.08, 1))

# 主布局：2 x 2
layout <- (pa | pc) / (b_panel | pd) + plot_layout(heights = c(0.85, 1.15))

# 加总标题
layout <- layout +
  plot_annotation(
    title = "Host prediction of Lycium barbarum virome vOTUs",
    subtitle = paste0(
      "Ensemble of four predictors; 20,892 vOTUs in the host-prediction table ",
      "(1,063 Plant-host, report denominator), \n",
      "20,443 with >= 1 predictor hit. Rule = taxonomy-class fallback (not an independent sequence hit)."
    ),
    theme = theme(
      plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.8, color = "#444444", hjust = 0, lineheight = 1.2),
      plot.margin = margin(6, 8, 6, 8)
    )
  )

# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# 输出
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
save_upset(layout, outdir, "Fig_GojiVirome_HostPrediction", width = 9.32, height = 9.8)
cat("Written:", file.path(outdir, "Fig_GojiVirome_HostPrediction.png"), "\n")
