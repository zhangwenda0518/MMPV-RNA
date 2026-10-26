#!/usr/bin/env Rscript
# ess_rhat.R — 由 ess_convergence_fig.py 调用的 R 后端 (2026-10-14)
#
# 输入: <indir>/chain*.csv   每链一个 CSV（首行参数名，无行号）
# 输出: <out_tsv>            列: param / rhat / upper / coda_ess
#
# 用法: Rscript ess_rhat.R <indir> <out_tsv>
#
# 口径
# ====
# gelman.diag(autoburnin = FALSE, multivariate = FALSE) —— 经典 Gelman-Rubin，
#   psrf[, 1] = R-hat 点估计，psrf[, 2] = 97.5% 上界。
#   ⚠ autoburnin = FALSE 是刻意的：burn-in 已由 Python 侧按 10% 统一裁掉，
#     若让 coda 再自动丢弃前 50%，两条口径会打架且样本量对不上。
#   ⚠ 该判据**不经过任何 ESS 算法** —— 这正是它的价值：ESS 低既可能是链未收敛，
#     也可能是 logEvery 过大导致样本太少；R-hat 只看链间/链内方差比，能把
#     「链未收敛」这一成因单独拎出来（见 utils/convergence_gate.py 的语义边界）。
#
# ⚠⚠ 坑（2026-10-14 实测，勿改回）
#   `sapply(ml, effectiveSize)` 返回的是 **参数 × 链 矩阵**（16×3），
#   逐参数跨链求和必须用 `rowSums`；用 `colSums` 得到的是 3 个数（逐链跨参数求和）
#   → `data.frame` 直接报 "arguments imply differing number of rows: 16, 3"。
#
# ⚠ 单参数退化：若某参数在各链内方差为 0（常量列），gelman.diag 会报错。
#   本脚本对每列先做退化检查，退化列直接剔除而不是让整脚本崩掉。

suppressPackageStartupMessages(library(coda))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript ess_rhat.R <indir> <out_tsv>\n")
  quit(status = 2)
}
indir <- args[1]
out_tsv <- args[2]

fs <- sort(list.files(indir, pattern = "[.]csv$", full.names = TRUE))
if (length(fs) < 2) {
  cat("⚠ 少于 2 条链，无法计算 R-hat\n")
  quit(status = 3)
}

chains <- lapply(fs, function(f) as.mcmc(read.csv(f, check.names = FALSE)))
ml <- mcmc.list(chains)

# 退化列剔除：链内方差为 0 的列会让 gelman.diag 报错
mat <- do.call(cbind, lapply(chains, as.matrix))
vars <- apply(mat, 2, var)
keep <- names(vars)[is.finite(vars) & vars > 0]
if (length(keep) < length(vars)) {
  cat("⚠ 剔除", length(vars) - length(keep), "个零方差列\n")
}
if (length(keep) < 1) {
  cat("⚠ 无有效列\n")
  quit(status = 4)
}
ml <- mcmc.list(lapply(chains, function(ch) as.mcmc(as.matrix(ch)[, keep, drop = FALSE])))

gd <- gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)
ps <- as.matrix(gd$psrf)
pn <- rownames(ps)
if (is.null(pn)) pn <- keep

es <- effectiveSize(ml)
es_vals <- as.numeric(es[pn])
es_vals[is.na(es_vals)] <- NA_real_

out <- data.frame(param = pn,
                  rhat = as.numeric(ps[, 1]),
                  upper = as.numeric(ps[, 2]),
                  coda_ess = es_vals,
                  stringsAsFactors = FALSE)

write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("written:", nrow(out), "rows ->", out_tsv, "\n")
cat("min R-hat =", round(min(out$rhat, na.rm = TRUE), 4),
    "(", out$param[which.min(out$rhat)], ")\n")
cat("max R-hat =", round(max(out$rhat, na.rm = TRUE), 4),
    "(", out$param[which.max(out$rhat)], ")\n")
cat("n > 1.1  :", sum(out$rhat > 1.1, na.rm = TRUE), "/", nrow(out), "\n")
cat("n > 1.01 :", sum(out$rhat > 1.01, na.rm = TRUE), "/", nrow(out), "\n")
