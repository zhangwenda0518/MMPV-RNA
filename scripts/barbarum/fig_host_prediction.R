#!/usr/bin/env Rscript
# fig_host_prediction.R — 图③ Host Prediction 图（补充）
#
# 数据：06_HostPrediction/。
#   面板 a  Host distribution（host_distribution.tsv）—— 与管线报告同源同值。
#   面板 b  预测器重叠 UpSet（ICTV / RVH / PB2 / Rule 四个证据源谁投了票）
#   面板 c  Decision_Method 构成（合成时由哪个证据主导）
#   面板 d  ICTV 置信度（High / Medium / Low / Unknown）
#
# 口径：面板 b 的 universe = ensemble_host_summary 全表 20,892 条；
#       四预测器任一命中者入图（20,443 条），四路全空者不入 UpSet。
#       ⚠ RuleClassTaxonomy 是"分类级回退"，非独立序列比对 → 面板 b/c/d 均标注。
#
# 来源：数据 = 项目 06_HostPrediction 管线产物；绘图代码为本次独立编写（PlotCase 简洁风格）。
# 说明：PlotCase MCP（127.0.0.1:17890）连接检查未通过，故未经案例库改编。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexUpset)
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
datadir <- get_arg("--data", "scripts/barbarum/upset_data")
outdir <- get_arg("--outdir", "reports/_discovery_sci_20261012/_lt/10_Reports/S10_SCI_Discovery")

# ---- 面板 a：宿主分布（完全采用管线 host_distribution.tsv 的值）----
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

# ---- 面板 b：预测器 UpSet ----
hp <- fread(file.path(datadir, "host_predictor_matrix.tsv"), sep = "\t")
preds <- c("ICTV", "RVH", "PB2", "Rule")
hp[, (preds) := lapply(.SD, as.logical), .SDcols = preds]
hp[, n_pred := rowSums(.SD), .SDcols = preds]
cat(sprintf("宿主预测器命中 %s / 全表 %d\n", format(nrow(hp), big.mark = ","), 20892))
print(table(hp$n_pred))

# ComplexUpset 返回 patchwork，labs() 会错位 → 用 grid 画一个独立题头条，
# 再与 UpSet 主体垂直拼接，保证标题稳定落在图上方。
library(grid)
title_grob <- function(title, subtitle) {
  gTree(children = gList(
    textGrob(title, x = 0, hjust = 0, y = 0.66, gp = gpar(fontsize = 9.5, fontface = "bold")),
    textGrob(subtitle, x = 0, hjust = 0, y = 0.24,
      gp = gpar(fontsize = 6.6, col = "#444444"))
  ))
}
pb_body <- upset(
  as.data.frame(hp), preds,
  name = "Host predictor",
  width_ratio = 0.30,
  min_size = 1,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "vOTUs" = intersection_size(
      counts = TRUE,
      mapping = aes(fill = as.factor(n_pred)),
      text = list(vjust = -0.5, size = 2.4)
    ) +
    scale_fill_manual(
      values = c("1" = "#F0E442", "2" = "#56B4E9", "3" = "#0072B2", "4" = "#D55E00"),
      name = "Predictors", drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
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
    "intersections" = theme_upset(8) +
      theme(
        axis.text.y = element_text(size = 7.5, color = "black"),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(4, 4, 2, 4)
      ),
    "overall_sizes" = theme_upset(8) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "black"),
        panel.grid.major.x = element_blank()
      )
  ))
)
pb <- wrap_elements(title_grob(
  "b  Host-predictor overlap",
  sprintf("ICTV lookup / RNAVirHost / PhaBOX2-Cherry / class-rule - n = %s with >= 1 hit",
    format(nrow(hp), big.mark = ","))
)) / pb_body +
  plot_layout(heights = c(0.11, 1))

# ---- 面板 c：Decision_Method 构成 ----
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

# ---- 面板 d：ICTV 置信度 ----
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

# ---- 组装 ----
layout <- (pa | pc) / (pb | pd) +
  plot_layout(heights = c(0.85, 1.15)) +
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

save_upset(layout, outdir, "Fig_GojiVirome_HostPrediction", width = 9.32, height = 9.4)
cat("Written:", file.path(outdir, "Fig_GojiVirome_HostPrediction.png"), "\n")
