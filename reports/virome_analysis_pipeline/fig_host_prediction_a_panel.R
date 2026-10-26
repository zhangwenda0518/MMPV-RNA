#!/usr/bin/env Rscript
# fig_host_prediction_a_panel.R — 宿主预测 panel a 独立图（2026-10-14 第 1 版）
#
# 来源：fig_host_prediction.R 的面板 a（Host prediction distribution，柱状图）。
# 用途：主图 v22 新增面板 O（与 BPanel = 面板 b 同批整合进 Fig_GojiVirome_Discovery）。
#
# 口径与 fig_host_prediction.R 逐字同源（host_distribution.tsv，OKABE_ITO 前 9 色）：
#   · 类目排序 = 计数降序（factor levels 按 -Count）；
#   · 柱顶数字千分位；x 标签 35° 防叠。
# ⚠ 独立图的标题块由 plot_annotation 出（同 fig_host_prediction_b_panel.R 的形态），
#   ggplot 内不再重复 title/subtitle —— 主图嵌入时会整块裁掉标题块（裁切线在
#   make_discovery_sci_figs.py::panel_O_hostclass 按像素自检）。
#
# 数据：upset_data/host_distribution.tsv
# 用法：
#   Rscript fig_host_prediction_a_panel.R [--data DIR] [--outdir DIR]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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

vfmt <- function(v) format(v, big.mark = ",", scientific = FALSE)

# ════════════════════════════════════════════════════════════════════════
# 面板 a：宿主分布（完全采用管线 host_distribution.tsv 的值）
# ════════════════════════════════════════════════════════════════════════
hd <- fread(file.path(datadir, "host_distribution.tsv"), sep = "\t")
stopifnot(nrow(hd) > 0, all(hd$Count > 0))
hd[, Host := factor(Host, levels = Host[order(-Count)])]
HOSTCOL <- setNames(OKABE_ITO[c(1, 2, 3, 4, 5, 6, 7, 8, 9)], levels(hd$Host))
cat(sprintf("Host distribution: n = %s | classes = %s\n",
    format(sum(hd$Count), big.mark = ","),
    paste(sprintf("%s %s", levels(hd$Host), vfmt(hd$Count[order(hd$Host)])),
          collapse = " / ")))

pa <- ggplot(hd, aes(x = Host, y = Count, fill = Host)) +
  geom_col(width = 0.74, color = "black", linewidth = 0.25) +
  geom_text(aes(label = vfmt(Count)), vjust = -0.45, size = 3.4) +
  scale_fill_manual(values = HOSTCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "vOTUs") +
  theme_upset(8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    panel.grid.major.x = element_blank()
  )

p <- pa +
  plot_annotation(
    title = "a  Host prediction distribution",
    subtitle = paste0(
      "Ensemble host call per vOTU (n = ", vfmt(sum(hd$Count)),
      "); colour = predicted host class"
    ),
    theme = theme(
      plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.8, color = "#444444", hjust = 0, lineheight = 1.2),
      plot.margin = margin(6, 8, 6, 8)
    )
  )

# 输出（画布与 BPanel 一致：9.32 × 6.0 in → 2796 × 1800 px @300dpi）
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
W_IN <- 9.32; H_IN <- 6.0
save_upset(p, outdir, "Fig_GojiVirome_HostPrediction_APanel",
           width = W_IN, height = H_IN)
cat("Written:", file.path(outdir, "Fig_GojiVirome_HostPrediction_APanel.png"), "\n")
