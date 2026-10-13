#!/usr/bin/env Rscript
# fig_ident_filter.R — 枸杞 discovery：「鉴定 / 过滤」双簇分组条形图
#
# 来源声明：PlotCase 案例库连接检查**未通过**（2026-10-13 复查 127.0.0.1:17890 无监听、
#           curl exit 7）→ 本图为**独立编写**，未经案例库改编。
#
# 口径（用户 2026-10-13 选定「鉴定与过滤双簇」+「过滤只展示保留总数」）：
#   左簇 鉴定 = 6 个有命中的病毒工具（viralverify / metabuli / blast / genomad / virbot /
#              rdrpcatch）+ 1 个类病毒轨道（viroid_blastn，* 标记，**不计入并集**）
#   右簇 过滤 = Uniprot 过滤后保留总数（02b_Filter）
#   参考线   = 鉴定并集 All_Candidate 556,773（= 过滤的输入）
#   y 轴     = log10（否则 604 / 2,229 / 3,219 在线性轴上是 0.1–0.6% 高度，肉眼不可见）
#
# ⚠ 两个已踩的坑（见报告 §9 规则 42–43）：
#   1) geom_col() 在 scale_y_log10() 下 ymin=0 → -Inf → 柱子被**整批静默移除**（只有 warning）。
#      修法 = 改用 geom_rect(ymin = 轴下限)。
#   2) 副标题里的 U+00B7 / U+2192 触发 mbcsToSbcs 字形转换失败 → 渲染成点。
#      修法 = 图内**只留 ASCII**。

script_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  # ⚠ 含中文路径下 normalizePath() 会把路径编坏 → 只用 dirname()
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  getwd()
}
SD <- script_dir()
source(file.path(SD, "upset_common.R"))

suppressPackageStartupMessages({ library(ggplot2); library(scales) })

DATA   <- file.path(SD, "identfilter_data")
REPDIR <- file.path(SD, "..", "..", "reports", "_discovery_sci_20261012", "_lt", "10_Reports")
OUTDIR <- file.path(SD, "..", "..", "reports", "_discovery_sci_20261012", "_lt",
                    "10_Reports", "S10_SCI_Discovery")

d <- read.delim(file.path(DATA, "ident_filter_summary.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE)

# 并集（= 过滤输入）直接从权威表取 TOTAL 行，避免在 R 里重算口径
tot <- read.delim(file.path(REPDIR, "ident_summary.tsv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
UNION <- as.numeric(tot$All_Candidate[tot$Sample == "TOTAL"])
stopifnot(length(UNION) == 1, !is.na(UNION), UNION == 556773)
RETAINED <- as.numeric(d$value[d$group == "Filtering"])

# ── x 轴：左簇 7 根 + 间隔位 + 右簇 1 根（用数值坐标，便于 geom_rect 与定位）──
id_cats  <- d$category[d$group == "Identification"]
id_cats  <- id_cats[order(-d$value[d$group == "Identification"])]
filt_cat <- d$category[d$group == "Filtering"]
xlevs    <- c(id_cats, "__gap__", filt_cat)
d$xn     <- match(d$category, xlevs)

# ⚠ 只有类病毒轨道加 "*"（首版误给全部鉴定工具都加了 *）
xlab <- id_cats
xlab[xlab == "viroid_blastn"] <- "viroid_blastn*"
xlab <- c(xlab, "", gsub(" \\(", "\n(", filt_cat))

# 颜色：鉴定（病毒工具）蓝 / 类病毒橙 / 过滤绿 —— 全部取自 OKABE_ITO
d$fillc <- ifelse(d$group == "Filtering", "#009E73",
           ifelse(d$track == "viroid",    "#E69F00", "#0072B2"))
# ⚠ 数值标签必须放进柱内：柱顶 452,904 与并集线 556,773 在对数轴上只差 0.09 decade，
#    标签放柱外会被并集虚线**横穿**（首版实测碰撞，见报告 §9 规则 44）。
#    白字配蓝/绿柱、黑字配橙柱（对比度：蓝 5.0 / 绿 3.2 / 橙黑 9.4）。
d$labcol <- ifelse(d$track == "viroid", "black", "white")

VIROID <- subset(d, track == "viroid")
GAPX   <- which(xlevs == "__gap__")
W      <- 0.34
YFLOOR <- 100
YTOP   <- 1.7e6
IDMID  <- mean(seq_along(id_cats))
FILTX  <- length(xlevs)

p <- ggplot(d) +
  geom_rect(aes(xmin = xn - W, xmax = xn + W, ymin = YFLOOR, ymax = value),
            fill = d$fillc, colour = NA) +
  # 类病毒轨道：橙色实柱 + 虚线描边，强调「独立轨道、不在并集内」
  geom_rect(data = VIROID, aes(xmin = xn - W, xmax = xn + W, ymin = YFLOOR, ymax = value),
            fill = NA, colour = "#7A4F00", linetype = "dashed", linewidth = 0.45) +
  geom_text(aes(x = xn, y = value, label = comma(value), colour = labcol),
            vjust = 1.55, size = 2.45, show.legend = FALSE) +
  scale_colour_identity() +
  geom_hline(yintercept = UNION, linetype = "dashed", colour = "#D55E00", linewidth = 0.55) +
  # ⚠ vjust=0 让文字**底边**落在该 y 上，否则线会横穿文字（首版即如此）
  annotate("text", x = FILTX + 0.6, y = UNION * 1.15, hjust = 1, vjust = 0, size = 2.7,
           colour = "#D55E00",
           label = paste0("union = ", comma(UNION), "  (filter input)")) +
  # 保留率放进绿柱内部（白字，折两行以适配柱宽），放柱外会与柱顶数值标签相撞
  annotate("text", x = FILTX, y = RETAINED * 0.24, size = 2.5, colour = "white", lineheight = 1.05,
           label = paste0(sprintf("%.1f", 100 * RETAINED / UNION), "%\nof union")) +
  geom_vline(xintercept = GAPX, linetype = "dotted", colour = "#999999", linewidth = 0.45) +
  annotate("text", x = IDMID, y = YTOP, label = "Identification", fontface = "bold", size = 3.2) +
  annotate("text", x = FILTX, y = YTOP, label = "Filtering", fontface = "bold", size = 3.2) +
  scale_y_log10(
    limits = c(YFLOOR, YTOP * 1.25),
    breaks = c(1e2, 1e3, 1e4, 1e5, 1e6),
    labels = label_number(big.mark = ",", accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_x_continuous(breaks = seq_along(xlevs), labels = xlab,
                     limits = c(0.4, FILTX + 0.7),
                     expand = expansion(mult = 0)) +
  labs(
    x = NULL,
    y = expression("sequences  (log"[10]*" scale)"),
    title = "Identification and filtering of Lycium barbarum virome candidates",
    subtitle = paste0(
      "438 samples; 10 identification tools (7 with hits) -> union ", comma(UNION),
      " viral sequences -> Uniprot filter retained ", comma(RETAINED),
      " (", sprintf("%.1f", 100 * RETAINED / UNION), "% of union)\n",
      "* viroid_blastn is a separate viroid track and is NOT part of the union;  ",
      "3 of 10 tools returned 0 hits (virsorter2 / virhunter / viralm)"
    )
  ) +
  theme_upset(base_size = 9) +
  theme(
    axis.text.x = element_text(size = 7.2, colour = "black"),
    panel.grid.major.x = element_blank()
  )

save_upset(p, OUTDIR, "Fig_GojiVirome_IdentFilter", width = 9.32, height = 5.4)

cat("--- 作图数值 ---\n")
print(d[order(d$group, -d$value), c("group", "track", "category", "value")], row.names = FALSE)
cat("union =", UNION, " retained =", RETAINED,
    " pct =", round(100 * RETAINED / UNION, 2), "%\n")
