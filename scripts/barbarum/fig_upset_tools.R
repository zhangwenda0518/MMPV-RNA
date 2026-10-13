#!/usr/bin/env Rscript
# fig_upset_tools.R — 图① 分类工具重叠 UpSet（枸杞 discovery 侧 05_Taxonomy）
#
# 数据：05_Taxonomy/Votus.integrated/intersection_matrix.tsv
#       （21,630 条 vOTU × 7 个分类工具，由 virus_classifier_analysis.R 导出）
# 图形：ComplexUpset —— 交集大小柱（按 degree 着色）+ 成员点阵（工具 × 交集）
#
# 版式决策：7 工具 → 最多 127 个非空交集，全画不可读。
#   → 只画 **按交集大小降序的前 TOP_N 个**，其余归并，
#     并在副标题如实写明"已显示前 N 个交集（共 M 个）"。
#
# 来源：数据 = 项目 05_Taxonomy 管线产物；绘图代码为本次独立编写（PlotCase 简洁风格）。
# 说明：PlotCase MCP（127.0.0.1:17890）连接检查未通过，故未经案例库改编。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
matrix_tsv <- get_arg("--matrix", "reports/_discovery_sci_20261012/data/05_Taxonomy/intersection_matrix.tsv")
outdir <- get_arg("--outdir", "reports/_discovery_sci_20261012/_lt/10_Reports/S10_SCI_Discovery")
top_n <- as.integer(get_arg("--top", "30"))

m <- fread(matrix_tsv, sep = "\t", header = TRUE)
tools <- TOOLS7[TOOLS7 %in% names(m)]
m[, (tools) := lapply(.SD, as.logical), .SDcols = tools]
m[, count := rowSums(.SD), .SDcols = tools]

cat(sprintf("总计 %s 条 vOTU x %d 个分类工具\n", format(nrow(m), big.mark = ","), length(tools)))
print(table(m$count))

# 交集签名 → 大小；取前 TOP_N 作显示阈值
sig <- m[, do.call(paste0, lapply(.SD, as.integer)), .SDcols = tools]
sizes <- sort(table(sig), decreasing = TRUE)
n_all <- length(sizes)
n_show <- min(top_n, n_all)
min_size <- as.integer(sizes[n_show])
cat(sprintf("非空交集共 %d 个；显示前 %d 个（size >= %d），其余 %d 个不显示\n",
  n_all, n_show, min_size, n_all - n_show))

p <- upset(
  as.data.frame(m), tools,
  name = "Classification tools",
  width_ratio = 0.26,
  min_size = min_size,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Contigs" = intersection_size(
      counts = TRUE,
      mapping = aes(fill = as.factor(count)),
      text = list(vjust = -0.5, size = 2.4)
    ) +
    scale_fill_manual(values = DEGREE_COLORS, name = "Tools overlapped", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    ylab("Contigs") +
    theme_upset(8) +
    theme(
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
  ),
  stripes = upset_stripes(mapping = aes(color = "black", fill = "grey92")),
  matrix = intersection_matrix(geom = geom_point(shape = 21, size = 2.6, stroke = 0.8)) +
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

p <- p +
  patchwork::plot_annotation(
    title = "Taxonomic-classification tool overlap",
    subtitle = sprintf(
      "Lycium barbarum virome vOTUs (n = %s) classified by %d tools - showing top %d of %d intersections; colour = number of tools sharing the vOTU",
      format(nrow(m), big.mark = ","), length(tools), n_show, n_all
    ),
    theme = theme(
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.5, color = "#444444", hjust = 0),
      plot.margin = margin(6, 8, 4, 8)
    )
  )

save_upset(p, outdir, "Fig_GojiVirome_UpsetTools", width = 9.32, height = 6.0)
cat("Written:", file.path(outdir, "Fig_GojiVirome_UpsetTools.png"), "\n")
