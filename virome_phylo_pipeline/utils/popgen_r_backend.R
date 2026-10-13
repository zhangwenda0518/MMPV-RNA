#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════
#  PLACEHOLDER — 未接线，nh/Hd 为占位值，勿用于发表
# ──────────────────────────────────────────────────────────────────────────
#  2026-10-12 审查确认: 本脚本在全仓库**无调用点** (仅 pixi.toml 注释与孤儿
#  审计清单提到它), 且下面的 `nh <- n; Hd <- 1.0` 是**占位假值** ——
#  真实单倍型数 / 单倍型多样性由 Python 端 pypopart 路径计算 (见
#  popgen_analysis.py; 本文件注释亦写明"单倍型计数由 python 端完成")。
#  谁若直接跑本脚本并把 JSON 里的 haplotypes/Hd 拿去做统计或进图,
#  得到的是编造结果。保留文件仅为记录 pegas 口径实现; 启用前必须先把
#  nh/Hd 换成真算 (pegas::haplotype / hap.div 需先解决含 gap 数据的 index 问题)。
# ══════════════════════════════════════════════════════════════════════════
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
