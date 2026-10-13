# upset_common.R — 枸杞 discovery 侧 UpSet 系列：共用配色与 ggplot 主题
# 配色与 virome_discovery_pipeline/utils/sci_palette.py 的 OKABE_ITO 同源（顺序固定，勿改）。

OKABE_ITO <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#000000", # black
  "#999999"  # grey
)

# 7 个分类工具的固定配色（供三张图统一，按 OKABE_ITO 顺序取前 7）
TOOLS7 <- c("diamond_lca", "mmseqs", "CAT", "metabuli", "genomad", "VITAP", "ACVirus")
TOOL_COLORS <- setNames(OKABE_ITO[1:7], TOOLS7)

# degree（重叠工具数）1..7 的连续配色：浅→深取自同色系，末档用灰示意"全票共识"
DEGREE_COLORS <- c(
  "1" = "#F0E442", "2" = "#B8DF8E", "3" = "#7CCD7C", "4" = "#56B4E9",
  "5" = "#449B99", "6" = "#0072B2", "7" = "#D55E00"
)

# 宿主预测器配色（图③）
HOSTPRED_COLORS <- c(
  "ICTV" = "#0072B2",
  "RVH" = "#D55E00",
  "PB2" = "#009E73",
  "Rule" = "#E69F00"
)

theme_upset <- function(base_size = 8) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.5, color = "black"),
      plot.title = ggplot2::element_text(size = base_size + 1.5, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, color = "#444444", hjust = 0),
      legend.title = ggplot2::element_text(size = base_size - 0.5),
      legend.text = ggplot2::element_text(size = base_size - 1),
      plot.margin = ggplot2::margin(6, 8, 6, 8)
    )
}

# 保存 PNG + PDF（源画布 = 排版宽度，避免字号被压小）
save_upset <- function(p, outdir, name, width, height) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  ggplot2::ggsave(file.path(outdir, paste0(name, ".png")), p,
    width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(outdir, paste0(name, ".pdf")), p,
    width = width, height = height, bg = "white")
  invisible(file.path(outdir, paste0(name, ".png")))
}
