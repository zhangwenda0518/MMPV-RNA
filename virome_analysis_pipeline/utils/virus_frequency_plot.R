#!/usr/bin/env Rscript

suppressWarnings({
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(optparse)
    if (!require("viridis", quietly = TRUE)) install.packages("viridis", dependencies = TRUE, repos="http://cran.rstudio.com/")
    library(viridis)
  })
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "all_viruses.best.summary.tsv", help = "输入的最优筛选结果文件"),
  make_option(c("-o", "--output"), type = "character", default = "virus_analysis_plots", help = "输出目录及文件前缀"),
  make_option(c("-w", "--width"), type = "numeric", default = 12, help = "图片宽度(英寸)"),
  make_option(c("-e", "--height"), type = "numeric", default = 8, help = "图片高度(英寸)"),
  make_option(c("-m", "--modes"), type = "character", default = "all", help = "绘图模式: MeanDepth, FPKM, RPM, TPM, 或 all"),
  make_option(c("-p", "--point-size"), type = "numeric", default = 3, help = "散点大小"),
  make_option(c("-d", "--dpi"), type = "numeric", default = 300, help = "图片分辨率"),
  make_option(c("-t", "--theme"), type = "character", default = "bw", help = "ggplot主题"),
  make_option(c("--log10-transform"), type = "logical", default = FALSE, action = "store_true", help = "使用log10对数坐标"),
  make_option(c("--multi-plot"), type = "logical", default = FALSE, action = "store_true", help = "生成多指标组合面板图"),
  make_option(c("--format"), type = "character", default = "pdf", help = "输出图片格式 (png, pdf)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$input)) stop(sprintf("错误: 输入文件不存在: %s", opt$input))
cat(sprintf("✅ 正在读取核心汇总数据: %s\n", opt$input))

data <- if(grepl("\\.csv$", opt$input)) read.csv(opt$input, check.names=F) else read.delim(opt$input, check.names=F)

# 智能生成展示名称 (融合 Taxonomy 和 Accession)
if ("taxonomy" %in% colnames(data)) {
  data$Display_Name <- ifelse(
    data$taxonomy != "Unannotated" & data$taxonomy != "-",
    paste0(data$taxonomy, "\n(", data$Virus, ")"),
    data$Virus
  )
} else {
  data$Display_Name <- data$Virus
}

# 文本自动换行处理 (防止 Taxonomy 过长)
data$Display_Name <- sapply(data$Display_Name, function(x) paste(strwrap(x, width = 40), collapse = "\n"))

required_columns <- list("MeanDepth" = c("Display_Name", "MeanDepth"), "FPKM" = c("Display_Name", "FPKM"), "RPM" = c("Display_Name", "RPM"), "TPM" = c("Display_Name", "TPM"))
modes <- if(opt$modes == "all") names(required_columns) else strsplit(trimws(opt$modes), ",")[[1]]
available_modes <- modes[sapply(modes, function(m) all(required_columns[[m]] %in% colnames(data)))]

if (length(available_modes) == 0) stop("❌ 错误: 未找到指定的分析指标列。")

output_dir <- dirname(opt$output)
if (output_dir != "." && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 动态判断是否需要翻转坐标轴 (种类太多时翻转能保证标签可读性)
unique_viruses <- length(unique(data$Display_Name))
should_flip <- unique_viruses > 5

prepare_data <- function(data, value_column) {
  data <- data[!is.na(data[[value_column]]), ]
  # 按照中位数排序，让图表更具阶梯美感
  value_medians <- aggregate(data[[value_column]] ~ Display_Name, data, median)
  data$Display_Name <- factor(data$Display_Name, levels = value_medians[order(value_medians[[2]]), "Display_Name"])
  return(data)
}

theme_func <- switch(opt$theme, "classic"=theme_classic, "minimal"=theme_minimal, "bw"=theme_bw, theme_bw)

create_boxplot <- function(data, value_column, title_suffix = "") {
  plot_data <- prepare_data(data, value_column)
  y_label <- value_column
  y_trans <- NULL
  
  if (opt$`log10-transform`) {
    min_pos <- min(plot_data[[value_column]][plot_data[[value_column]] > 0], na.rm = TRUE)
    plot_data[[value_column]][plot_data[[value_column]] <= 0] <- min_pos / 10
    y_trans <- scale_y_log10()
    y_label <- paste0(y_label, " (log10)")
  }
  
  p <- ggplot(plot_data, aes(x = Display_Name, y = .data[[value_column]])) +
    geom_boxplot(aes(fill = Display_Name), alpha = 0.6, outlier.shape = NA, width = 0.6) +
    geom_point(aes(color = Display_Name), position = position_jitter(width = 0.2, height = 0), size = opt$`point-size`, alpha = 0.8) +
    scale_fill_viridis_d(option = "turbo") + 
    scale_color_viridis_d(option = "turbo") +
    labs(title = paste("Virus Abundance", title_suffix), x = "Identified Virus (Taxonomy & Accession)", y = y_label) +
    theme_func(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  if (should_flip) {
    p <- p + coord_flip() + theme(axis.text.y = element_text(size=10, face="italic"))
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10, face="italic"))
  }
  
  if (!is.null(y_trans)) p <- p + y_trans + annotation_logticks(sides = ifelse(should_flip, "b", "l"))
  return(p)
}

# 1. 逐一生成单指标分布图
for (mode in available_modes) {
  p <- create_boxplot(data, mode, paste("-", mode, "Distribution"))
  fn <- sprintf("%s_%s_boxplot%s.%s", opt$output, tolower(mode), ifelse(opt$`log10-transform`, "_log10", ""), opt$format)
  # 如果翻转了坐标轴，高度可以随着种类数量动态伸缩
  dynamic_height <- ifelse(should_flip, max(opt$height, unique_viruses * 0.8), opt$height)
  ggsave(fn, plot=p, width=opt$width, height=dynamic_height, dpi=opt$dpi, bg="white")
  fn_pdf <- sub("\\.[^.]+$", ".pdf", fn); if (fn != fn_pdf) ggsave(fn_pdf, plot=p, width=opt$width, height=dynamic_height, dpi=opt$dpi, bg="white")
  fn_png <- sub("\\.[^.]+$", ".png", fn); if (fn != fn_png) ggsave(fn_png, plot=p, width=opt$width, height=dynamic_height, dpi=opt$dpi, bg="white")
  cat(sprintf("📊 成功生成图表: %s\n", fn))
}

# 2. 生成多指标综合面板图
if (opt$`multi-plot` && length(available_modes) > 1) {
  plot_data_long <- data %>% select(Display_Name, Sample, all_of(available_modes)) %>%
    pivot_longer(cols = all_of(available_modes), names_to = "Metric", values_to = "Value") %>% filter(!is.na(Value))
  plot_data_long$Metric <- factor(plot_data_long$Metric, levels = available_modes)
  
  if (opt$`log10-transform`) {
    for (metric in unique(plot_data_long$Metric)) {
      md <- plot_data_long[plot_data_long$Metric == metric, ]
      min_p <- min(md$Value[md$Value > 0], na.rm=TRUE)
      plot_data_long$Value[plot_data_long$Metric == metric & plot_data_long$Value <= 0] <- min_p / 10
    }
  }
  
  # 按照平均丰度(以第一个Metric基准)排序
  base_metric <- available_modes[1]
  medians <- aggregate(Value ~ Display_Name, data=plot_data_long[plot_data_long$Metric==base_metric,], median)
  plot_data_long$Display_Name <- factor(plot_data_long$Display_Name, levels=medians[order(medians$Value), "Display_Name"])
  
  p_facet <- ggplot(plot_data_long, aes(x=Display_Name, y=Value)) +
    geom_boxplot(aes(fill=Display_Name), alpha=0.6, outlier.shape=NA) +
    geom_point(aes(color=Display_Name), position=position_jitter(width=0.2, height=0), alpha=0.6) +
    facet_wrap(~ Metric, scales="free_x", ncol=length(available_modes)) +
    scale_fill_viridis_d(option="turbo") + scale_color_viridis_d(option="turbo") +
    labs(title="Comprehensive Multi-metric Virus Abundance", x="Taxonomy (Accession)", y="Value") +
    theme_func(base_size=13) +
    theme(
      plot.title=element_text(hjust=0.5, face="bold"), 
      axis.text.y=element_text(size=10, face="italic"),
      strip.text = element_text(face="bold", size=12),
      legend.position="none"
    ) + coord_flip() # 综合面板图强制横向展示，更利于多指标对比
  
  if (opt$`log10-transform`) p_facet <- p_facet + scale_y_log10()
  
  fn_f <- sprintf("%s_multi_metrics%s.%s", opt$output, ifelse(opt$`log10-transform`, "_log10", ""), opt$format)
  dynamic_height_facet <- max(opt$height, unique_viruses * 0.8)
  ggsave(fn_f, plot=p_facet, width=opt$width * 1.5, height=dynamic_height_facet, dpi=opt$dpi, bg="white")
  fn_f_pdf <- sub("\\.[^.]+$", ".pdf", fn_f); if (fn_f != fn_f_pdf) ggsave(fn_f_pdf, plot=p_facet, width=opt$width * 1.5, height=dynamic_height_facet, dpi=opt$dpi, bg="white")
  fn_f_png <- sub("\\.[^.]+$", ".png", fn_f); if (fn_f != fn_f_png) ggsave(fn_f_png, plot=p_facet, width=opt$width * 1.5, height=dynamic_height_facet, dpi=opt$dpi, bg="white")
  cat(sprintf("📈 成功生成综合面板图: %s\n", fn_f))
}
cat("🎉 所有可视化任务完成！\n")
