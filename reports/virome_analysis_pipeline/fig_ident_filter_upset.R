#!/usr/bin/env Rscript
# fig_ident_filter_upset.R — 图④「鉴定 x 过滤」堆叠 UpSet（枸杞 discovery 侧）
#
# 集合（7）= 6 个有命中的病毒鉴定工具 + 1 条独立类病毒轨道：
#   viralverify / metabuli / blast / genomad / virbot / rdrpcatch / viroid_blastn
# 交集柱 = **堆叠**，段 = 过滤结局（柱内分段，用户 2026-10-13 选定）：
#   Retained            通过 CDD + UniProt 两道过滤，最终保留
#   Removed by UniProt  通过 CDD、被 UniProt 滤掉
#   Removed by CDD      被 CDD 滤掉（未进入 cdd_kept）
#   Viroid track        viroid_blastn 命中且**不在**病毒过滤链内（独立轨道）
#
# 数据：identfilter_data/ident_filter_upset_signature.tsv（由 make_ident_filter_upset_data.py
#       在服务器上从 02a_Identification / 02b_Filter 逐样本汇总，395 个样本）。
# 口径（实测，见报告 §9）：
#   Σunion 556,950 ⊇ Σcdd_kept 473,523 ⊇ Σretained 94,826（= union 的 17.03%）
#   viroid∩union = 88 ≠ 0（类病毒轨道与病毒轨道**确有**重叠，非完全互斥）
#   5 个样本无 02b 产物（27 条序列）→ 过滤结局未知，**剔除不入图**，在副标题声明
#
# 版式：64 个非空交集，全画不可读 → 只画按大小降序的**前 TOP_N 个**，
#       并在副标题如实写明 "top N of M"。
#
# 来源：数据 = 项目 discovery 管线产物；绘图代码为本次独立编写（PlotCase 简洁风格）。
# 说明：PlotCase MCP（127.0.0.1:17890）需用户自行启动，故未经案例库改编。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexUpset)
  library(grid)
})
script_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  # ⚠ 含中文路径下 normalizePath() 会把路径编坏 → 只用 dirname()
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  getwd()
}
SD <- script_dir()
source(file.path(SD, "upset_common.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d) {
  i <- which(args == k)
  if (length(i)) args[i + 1] else d
}
datadir <- get_arg("--data", file.path(SD, "identfilter_data"))
outdir <- get_arg("--outdir",
  file.path(SD, "..", "..", "reports", "_discovery_sci_20261012", "_lt",
            "10_Reports", "S10_SCI_Discovery"))
top_n <- as.integer(get_arg("--top", "20"))

TRACKS <- c("viralverify", "metabuli", "blast", "genomad", "virbot", "rdrpcatch",
            "viroid_blastn")
STATUS_LEVELS <- c("Retained", "Removed by UniProt", "Removed by CDD", "Viroid track")
STATUS_COLORS <- c(
  "Retained"           = "#009E73",  # 绿：通过两道过滤
  "Removed by UniProt" = "#BFBFBF",  # 浅灰：CDD 过、UniProt 滤掉
  "Removed by CDD"     = "#7A7A7A",  # 深灰：CDD 即滤掉
  "Viroid track"       = "#E69F00"   # 橙：独立类病毒轨道（不在病毒过滤链内）
)

# ⚠ signature 是 7 位 0/1 串（如 "0100000"）→ fread 会当整数读、**吃掉前导零**，
#   必须强制按字符读，否则 nchar 校验失败 / 位序整体错位。
sig <- fread(file.path(datadir, "ident_filter_upset_signature.tsv"), sep = "\t",
             colClasses = list(character = "signature"))
stopifnot(all(nchar(sig$signature) == length(TRACKS)))
stopifnot(all(c("signature", "n_total", "n_retained", "n_cdd_only", "n_removed",
                "n_viroid", "n_nofilter") %in% names(sig)))

n_all <- nrow(sig)
n_nofilter <- sum(sig$n_nofilter)
n_plotted <- sum(sig$n_total) - n_nofilter
n_ret <- sum(sig$n_retained)
n_union <- sum(sig$n_total) - sum(sig$n_viroid) - n_nofilter  # union 内（不含 viroid-only）

cat(sprintf("签名数 %d｜序列（入图）%s｜其中 retained %s（%.2f%% of union）\n",
            n_all, format(n_plotted, big.mark = ","), format(n_ret, big.mark = ","),
            100 * n_ret / n_union))
cat(sprintf("剔除 nofilter = %d 条（5 个样本无 02b 产物）\n", n_nofilter))

# ── 展开：签名 -> 逐序列行（ComplexUpset 需要逐元素表）──
exp <- sig[, {
  st <- c(
    rep("Retained", n_retained),
    rep("Removed by UniProt", n_cdd_only),
    rep("Removed by CDD", n_removed),
    rep("Viroid track", n_viroid)
  )
  .(status = st)
}, by = signature]
for (i in seq_along(TRACKS)) {
  exp[[TRACKS[i]]] <- substr(exp$signature, i, i) == "1"
}
exp[, status := factor(status, levels = STATUS_LEVELS)]
stopifnot(nrow(exp) == n_plotted)
cat(sprintf("展开行数 %s\n", format(nrow(exp), big.mark = ",")))

# ── 显示阈值：前 TOP_N 个交集 ──
sizes <- sort(sig$n_total, decreasing = TRUE)
n_show <- min(top_n, n_all)
min_size <- as.integer(sizes[n_show])
cum_shown <- sum(sizes[seq_len(n_show)])
cat(sprintf("显示前 %d/%d 个交集（size >= %d），覆盖 %s 条 = %.2f%%\n",
            n_show, n_all, min_size, format(cum_shown, big.mark = ","),
            100 * cum_shown / sum(sig$n_total)))

# ── 柱顶标注用：前 TOP_N 个交集（顺序须与 ComplexUpset 的 descending-cardinality 一致）──
# ⚠ 自检：前 TOP_N 个交集大小必须互不相同 —— 有并列时 ComplexUpset 的先后与
#   本表顺序可能不一致，柱顶 % 会整体错位。并列时应改用「按交集键查表」而非按序号。
sig_top <- head(sig, n_show)
if (anyDuplicated(sig_top$n_total)) {
  stop(sprintf("前 %d 个交集存在并列大小（%s）→ 柱顶标注会错位，请调 --top 或改查表实现",
               n_show, paste(unique(sig_top$n_total[duplicated(sig_top$n_total)]), collapse = ", ")))
}
sig_top[, bar_h := n_total - n_nofilter]
sig_top[, pct := 100 * n_retained / pmax(bar_h, 1)]
sig_top[, pct_lab := ifelse(n_viroid > 0, "viroid",
                            sprintf("%.0f%%", pct))]
cat("柱顶标注:\n")
print(sig_top[, .(signature, bar_h, pct_lab)])

df <- as.data.frame(exp[, c(TRACKS, "status"), with = FALSE])

# 集合大小写进行标签（替代 ComplexUpset 默认那个只显示裸数字的左侧面板）
# ⚠ 标签必须**单行**：两行标签的高度超过点阵行距 → 会与上一行互相压字。
set_sizes_n <- vapply(seq_along(TRACKS), function(i) {
  sum(sig$n_total[substr(sig$signature, i, i) == "1"])
}, numeric(1))
names(set_sizes_n) <- TRACKS
LABELS <- setNames(
  sprintf("%s (%s)", TRACKS, format(set_sizes_n, big.mark = ",")), TRACKS)
setnames(df, TRACKS, unname(LABELS[TRACKS]))
cat("集合大小:\n"); print(set_sizes_n)

# ⚠ ComplexUpset 会把 intersect 给的顺序**上下翻转**再画 → 传 rev() 才能让
#   viralverify（最大集）落在最上面。
SET_ORDER <- rev(TRACKS)

p_body <- upset(
  df, unname(LABELS[SET_ORDER]),
  name = "Identification track",
  width_ratio = 0.30,
  min_size = min_size,
  sort_sets = FALSE,
  set_sizes = FALSE,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Sequences" = intersection_size(
      counts = FALSE,
      mapping = aes(fill = status)
    ) +
      # 每根柱顶标「最终保留率」—— 全图的定量结论：工具共识越高，保留率越高
      # （6 工具交集 98.5% vs 单工具 viralverify 13.0%）
      # ⚠ x 位置 = 交集在 descending-cardinality 序下的 1..N；实测前 20 个交集
      #   大小互不相同，ComplexUpset 的排序与 sig 表顺序一致，故可直接用序号定位。
      annotate("text",
               x = seq_len(nrow(sig_top)), y = sig_top$bar_h,
               label = sig_top$pct_lab,
               vjust = -0.45, size = 2.05, colour = "#1A5C3A") +
      scale_fill_manual(values = STATUS_COLORS, name = NULL, drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
      ylab("Sequences") +
      theme_upset(8) +
      theme(
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.26, "cm"),
        legend.text = element_text(size = 7)
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
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  ))
)

title_grob <- function(title, subtitle) {
  gTree(children = gList(
    textGrob(title, x = 0, hjust = 0, y = 0.72,
             gp = gpar(fontsize = 10.5, fontface = "bold")),
    textGrob(subtitle, x = 0, hjust = 0, y = 0.22,
             gp = gpar(fontsize = 6.6, col = "#444444", lineheight = 1.15))
  ))
}
n_cdd <- n_ret + sum(sig$n_cdd_only)
sub_txt <- paste0(
  "Lycium barbarum virome: which identification tool(s) found each contig, and what survived the filter chain\n",
  "395 samples; union ", format(n_union, big.mark = ","),
  " viral candidates -> CDD ", format(n_cdd, big.mark = ","),
  " kept -> UniProt ", format(n_ret, big.mark = ","),
  " kept (", sprintf("%.1f", 100 * n_ret / n_union), "% of union). ",
  "Showing top ", n_show, " of ", n_all, " intersections (",
  sprintf("%.2f", 100 * cum_shown / sum(sig$n_total)), "% of sequences).\n",
  "viroid_blastn is a separate viroid track; it overlaps the viral tools in 88 contigs. ",
  "27 contigs from 5 samples lacking 02b output are excluded (filter outcome unknown)."
)

layout <- wrap_elements(title_grob(
  "Identification x filtering of L. barbarum virome candidates", sub_txt
)) / p_body +
  plot_layout(heights = c(0.16, 1))

save_upset(layout, outdir, "Fig_GojiVirome_UpsetIdentFilter", width = 9.32, height = 6.8)
cat("Written:", file.path(outdir, "Fig_GojiVirome_UpsetIdentFilter.png"), "\n")
