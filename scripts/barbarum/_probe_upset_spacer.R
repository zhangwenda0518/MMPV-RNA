#!/usr/bin/env Rscript
# _probe_upset_spacer.R — 实测 ComplexUpset 生成的 patchwork 里各格子的 npc 位置。
# 目的：确定 (行1,列1) spacer 格在图①/图② 里的确切位置与尺寸，
#       供「补工具耗时条」时以英寸钉框（不写轴分数）。
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ComplexUpset)
})
SD <- "scripts/barbarum"
source(file.path(SD, "upset_common.R"))
m <- fread(file.path(SD, "upset_data/agree_matrix_Species.tsv"), sep = "\t")
tools <- TOOLS7[TOOLS7 %in% names(m)]
m[, (tools) := lapply(.SD, as.logical), .SDcols = tools]
m[, count := rowSums(.SD), .SDcols = tools]
sig <- m[, do.call(paste0, lapply(.SD, as.integer)), .SDcols = tools]
sizes <- sort(table(sig), decreasing = TRUE); min_size <- as.integer(sizes[18])
p <- upset(as.data.frame(m), tools, name="Classification tools", width_ratio=0.26,
  min_size=min_size, encode_sets=FALSE,
  sort_intersections="descending", sort_intersections_by="cardinality",
  base_annotations=list("Agreed vOTUs"=intersection_size(counts=TRUE)+ylab("n")),
  stripes=upset_stripes(mapping=aes(color="black", fill="grey92")),
  matrix=intersection_matrix(geom=geom_point(shape=21,size=2.6,stroke=0.8)))
cat("class:", class(p), "\n")
cat("patches$plots 数:", length(p$patches$plots), "\n")
for (i in seq_along(p$patches$plots)) {
  x <- p$patches$plots[[i]]
  cat(sprintf("  [%d] class=%s  inherits_spacer=%s\n", i,
      paste(class(x), collapse="/"), inherits(x, "spacer")))
}
# gtable 各格位置
g <- patchwork:::patchworkGrob(p)
cat("\n--- patchworkGrob 的 layout ---\n")
print(g$layout[, c("t","l","b","r","name")])
cat("\n--- widths ---\n"); print(g$widths)
cat("\n--- heights ---\n"); print(g$heights)
