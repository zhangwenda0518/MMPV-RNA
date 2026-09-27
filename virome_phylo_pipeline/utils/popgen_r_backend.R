#!/usr/bin/env Rscript
# popgen_r_backend.R — 调用 pegas 计算标准群体遗传统计 (DnaSP 口径)
# 用法: Rscript popgen_r_backend.R <fasta> <json输出>
suppressMessages({ library(pegas); library(ape) })
args <- commandArgs(trailingOnly = TRUE)
fasta <- args[1]
out <- args[2]

aln <- read.dna(fasta, format = "fasta", as.matrix = TRUE)
n <- nrow(aln); L <- ncol(aln)

# 单倍型计数由 python 端完成 (pegas haplotype 对含 gap 数据的 index 格式不稳定)
nh <- n
Hd <- 1.0

pi_val <- nuc.div(aln)
seg <- seg.sites(aln)
S <- length(seg)
a1 <- sum(1 / (1:(n - 1)))
theta_w <- S / (a1 * L)
k <- mean(dist.dna(aln, "N", pairwise.deletion = TRUE))

tt <- tryCatch(tajima.test(aln), error = function(e) NULL)
D <- if (!is.null(tt)) tt$D else NaN
p_beta <- if (!is.null(tt)) tt$Pval.beta else NaN

res <- list(n = n, L = L, haplotypes = nh, Hd = round(Hd, 6),
            pi = round(pi_val, 6), S = S, theta_w = round(theta_w, 6),
            k = round(k, 4), tajima_D = round(D, 4), tajima_p_beta = round(p_beta, 4))
write(jsonlite::toJSON(res, auto_unbox = TRUE), out)
cat("OK:", out, "\n")
