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
# panvirome_oncoprint.R — 病毒准种突变瀑布图 (ComplexHeatmap OncoPrint)
# 用法: Rscript panvirome_oncoprint.R --input Matrix_06.csv --output Figure_5.pdf

options(warn = -1)
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(grid)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "输入 CSV: Matrix_06_Ultimate_OncoPrint_State.csv"),
  make_option(c("-o", "--output"), type = "character", default = "Figure_5_Ultimate_Quasispecies_OncoPrint.pdf",
              help = "输出 PDF 路径"),
  make_option(c("-t", "--title"), type = "character", default = "Figure 5: Mutational Waterfall Landscape",
              help = "图标题"),
  make_option(c("-w", "--width"), type = "numeric", default = 16, help = "PDF 宽度 (英寸)"),
  make_option(c("-H", "--height"), type = "numeric", default = 10, help = "PDF 高度 (英寸)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# 突变类型配色
col_map <- c(
  "Fixed_Missense"     = "#E31A1C", "Major_Missense"  = "#FC4E2A", "Minor_Missense"     = "#FD8D3C",
  "Fixed_Synonymous"   = "#1F78B4", "Major_Synonymous" = "#41B6C4", "Minor_Synonymous"   = "#A1DAB4",
  "Fixed_Regulatory"   = "#33A02C", "Major_Regulatory" = "#74C476", "Minor_Regulatory"   = "#C7E9C0"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(1, "pt"), h - unit(1, "pt"), gp = gpar(fill = "#F5F5F5", col = NA))
  },
  Fixed_Missense   = function(x, y, w, h) { grid.rect(x, y, w, h,        gp = gpar(fill = col_map["Fixed_Missense"],   col = NA)) },
  Major_Missense   = function(x, y, w, h) { grid.rect(x, y, w, h * 0.75, gp = gpar(fill = col_map["Major_Missense"],   col = NA)) },
  Minor_Missense   = function(x, y, w, h) { grid.rect(x, y, w, h * 0.4,  gp = gpar(fill = col_map["Minor_Missense"],   col = NA)) },
  Fixed_Synonymous = function(x, y, w, h) { grid.rect(x, y, w, h,        gp = gpar(fill = col_map["Fixed_Synonymous"], col = NA)) },
  Major_Synonymous = function(x, y, w, h) { grid.rect(x, y, w, h * 0.75, gp = gpar(fill = col_map["Major_Synonymous"], col = NA)) },
  Minor_Synonymous = function(x, y, w, h) { grid.rect(x, y, w, h * 0.4,  gp = gpar(fill = col_map["Minor_Synonymous"], col = NA)) },
  Fixed_Regulatory = function(x, y, w, h) { grid.rect(x, y, w, h,        gp = gpar(fill = col_map["Fixed_Regulatory"], col = NA)) },
  Major_Regulatory = function(x, y, w, h) { grid.rect(x, y, w, h * 0.75, gp = gpar(fill = col_map["Major_Regulatory"], col = NA)) },
  Minor_Regulatory = function(x, y, w, h) { grid.rect(x, y, w, h * 0.4,  gp = gpar(fill = col_map["Minor_Regulatory"], col = NA)) }
)

cat(sprintf("  输入: %s\n", opt$input))
onco_mat <- as.matrix(read.csv(opt$input, row.names = 1, check.names = FALSE))
cat(sprintf("  矩阵: %d 变异 x %d 样本\n", nrow(onco_mat), ncol(onco_mat)))

ht <- oncoPrint(onco_mat,
  alter_fun = alter_fun, col = col_map,
  remove_empty_columns = TRUE, remove_empty_rows = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 9),
  pct_side = "right", row_names_side = "left",
  column_title = opt$title,
  heatmap_legend_param = list(title = "Fixation State", at = names(col_map))
)
# 2026-10-20: 四格式统一 (SVG/PNG/PDF/PNG@1200dpi); 每格式各重画一次
draw_figure_formats(function() draw(ht), opt$output,
                    width = opt$width, height = opt$height, bg = "white")
cat(sprintf("  → %s\n", opt$output))
