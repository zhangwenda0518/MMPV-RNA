#!/usr/bin/env Rscript
# ess_coda_synth.R — 由 ess_convergence_fig.py 调用的 coda 对拍后端 (2026-10-14)
#
# 输入: <indir>/s*.csv   每文件一列 `x`（单条序列）
# 输出: <out_tsv>        列: file / ess_coda
#
# 用法: Rscript ess_coda_synth.R <indir> <out_tsv>
#
# 用途
# ====
# 在**已知答案的合成序列**上，用 coda（R 权威库）独立算 ESS，与本实现
# （utils/ess.py，Geyer IPS）对拍。
#
# ⚠⚠ 口径：coda 用的是 **AR 谱密度**（`spectrum0.ar`），我们用的是 **Geyer
#   initial positive sequence** —— 这是**两种不同算法**。所以本对拍只能作
#   「量级一致」的旁证，**不是位级判据**。真正的硬证据是：
#     ① 答案已知的金标准（i.i.d. 的 ESS 应 = n）；
#     ② 同算法的独立重写（utils/ess_convergence_fig.ess_ips_independent）。
#
# ⚠ 白噪声的特殊性：coda 对 i.i.d. 序列**精确返回 n**，因此它同时是
#   「序列实际长度」的第三方验证手段（曾用它证实 AR(1) 序列长度为 26503
#   而非生成长度 27003）。

suppressPackageStartupMessages(library(coda))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript ess_coda_synth.R <indir> <out_tsv>\n")
  quit(status = 2)
}
indir <- args[1]
out_tsv <- args[2]

fs <- sort(list.files(indir, pattern = "[.]csv$", full.names = TRUE))
if (length(fs) < 1) {
  cat("⚠ 无输入序列\n")
  quit(status = 3)
}

files <- character(0)
vals <- numeric(0)
for (f in fs) {
  d <- read.csv(f, check.names = FALSE)
  v <- as.numeric(d[[1]])
  files <- c(files, basename(f))
  vals <- c(vals, effectiveSize(as.mcmc(v)))
}

out <- data.frame(file = files, ess_coda = vals, stringsAsFactors = FALSE)
write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("written:", nrow(out), "rows ->", out_tsv, "\n")
