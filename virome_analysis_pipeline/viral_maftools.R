#!/usr/bin/env Rscript

# =========================================================
# 脚本名称：master_viral_maftools.R (The Masterpiece Edition v20)
# 功能：修复 NA% 错误，新增单图无缝拼接拼接，以及唯一的全基因组坐标棒棒糖总图
# =========================================================

suppressPackageStartupMessages({
  if (!require("optparse", quietly = TRUE)) { install.packages("optparse", repos="https://mirrors.westlake.edu.cn/CRAN/", quiet=TRUE); library(optparse) }
  if (!require("data.table", quietly = TRUE)) { install.packages("data.table", repos="https://mirrors.westlake.edu.cn/CRAN/", quiet=TRUE); library(data.table) }
  if (!require("maftools", quietly = TRUE)) { BiocManager::install("maftools", ask=FALSE, quiet=TRUE); library(maftools) }
  if (!require("ggplot2", quietly = TRUE)) { install.packages("ggplot2", repos="https://mirrors.westlake.edu.cn/CRAN/", quiet=TRUE); library(ggplot2) }
  if (!require("ggrepel", quietly = TRUE)) { install.packages("ggrepel", repos="https://mirrors.westlake.edu.cn/CRAN/", quiet=TRUE); library(ggrepel) }
  if (!require("patchwork", quietly = TRUE)) { install.packages("patchwork", repos="https://mirrors.westlake.edu.cn/CRAN/", quiet=TRUE); library(patchwork) }
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入路径"),
  make_option(c("-o", "--output"), type = "character", default = "maftools_out", help = "输出目录前缀"),
  make_option(c("-g", "--gene"), type = "character", default = NULL, help = "棒棒糖图指定目标"),
  make_option(c("--top"), type = "integer", default = 15, help = "Oncoplot 最大频数"),
  make_option(c("--min_mut"), type = "integer", default = 0, help = "低负荷剔除阈值")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) quit(status = 1)

out_dir <- opt$output
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
output_prefix <- file.path(out_dir, "maftools")

tcga_native_non_syn <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Targeted_Region", "5'Flank", "3'Flank", "IGR", "3'UTR", "5'UTR", "Intron")
custom_colors <- c("Missense_Mutation"="#33A02C", "Frame_Shift_Del"="#E31A1C", "Frame_Shift_Ins"="#FF7F00", "Nonsense_Mutation"="#1F78B4", "Nonstop_Mutation"="#A6CEE3", "Translation_Start_Site"="#000000", "In_Frame_Ins"="#FDBF6F", "In_Frame_Del"="#CAB2D6", "Silent"="#CCCCCC", "Splice_Site"="#6A3D9A", "5'Flank"="#B15928", "3'Flank"="#E7298A", "IGR"="#E7298A", "3'UTR"="#B15928", "5'UTR"="#E7298A", "Intron"="#A6CEE3", "Targeted_Region"="#808080")

# =================模块：推断引擎=================
get_ncbi_lengths <- function(acc) {
  if(is.na(acc) || acc == "auto" || acc == "Viral_Consensus") return(NULL)
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=", acc, "&rettype=ft&retmode=text")
  lines <- tryCatch(readLines(url, warn=FALSE, n=2000), error=function(e) NULL)
  if(is.null(lines) || length(lines) == 0) return(NULL)
  genes <- list(); curr_len <- NA
  for(line in lines) {
    if(grepl("^[<>]?[0-9]+\t[>0-9]+\tCDS", line)) {
      pts <- strsplit(line, "\t")[[1]]; s <- as.numeric(gsub("[^0-9]", "", pts[1])); e <- as.numeric(gsub("[^0-9]", "", pts[2]))
      if(!is.na(s) && !is.na(e)) curr_len <- floor((abs(e - s) + 1)/3)
    } else if(grepl("^\t\t\tgene\t", line) && !is.na(curr_len)) {
      g <- sub("^\t\t\tgene\t", "", line); genes[[g]] <- curr_len; curr_len <- NA
    } else if(grepl("^\t\t\tproduct\t", line) && !is.na(curr_len)) {
      g <- sub("^\t\t\tproduct\t", "", line)
      if(is.null(genes[[g]])) genes[[g]] <- curr_len
    }
  }
  if(length(genes)==0) return(NULL)
  data.frame(HGNC = names(genes), protein.length = unlist(genes), stringsAsFactors = FALSE)
}

get_empirical_lengths <- function(maf_dt) {
  dt <- copy(maf_dt)
  dt[, aa_pos := suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", Protein_Change)))]
  gl <- dt[!is.na(aa_pos), .(max_pos = max(aa_pos, na.rm=TRUE)), by = Hugo_Symbol]
  if(nrow(gl) == 0) return(NULL)
  gl[, est_len := floor(max_pos * 1.1) + 10]
  data.frame(HGNC = gl$Hugo_Symbol, protein.length = gl$est_len, stringsAsFactors = FALSE)
}

# =================模块：高定单基因棒棒糖引擎=================
draw_viral_lollipop <- function(maf_dt, n_samples, gene, prot_len_table, custom_colors) {
    muts <- copy(maf_dt[Hugo_Symbol == gene])
    target_class <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Silent", "Nonstop_Mutation", "Translation_Start_Site")
    muts <- muts[Variant_Classification %in% target_class]
    if(nrow(muts) == 0) return(NULL)
    
    muts[, aa_pos := suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", Protein_Change)))]
    muts <- muts[!is.na(aa_pos)]
    if(nrow(muts) == 0) return(NULL)
    
    # 修复 NA%: 参数 n_samples 安全传入，不再受限于 summary 对象错位！
    mutated_samples <- length(unique(muts$Tumor_Sample_Barcode))
    mut_rate_str <- sprintf("%.2f%%", (mutated_samples / n_samples) * 100)
    
    counts <- muts[, .N, by = .(aa_pos, Variant_Classification, Protein_Change)]
    max_n <- max(counts$N, na.rm = TRUE)
    counts[, pt_size := ifelse(Variant_Classification == "Silent", 1.0, 2.5 + (N/max_n) * 2.5)]
    counts[, pt_alpha := ifelse(Variant_Classification == "Silent", 0.2, 0.4 + (N/max_n) * 0.6)]
    counts[, line_alpha := ifelse(Variant_Classification == "Silent", 0.1, 0.3 + (N/max_n) * 0.4)]
    
    setorder(counts, -N); top_hits <- head(counts[Variant_Classification != "Silent"], 10)
    len <- prot_len_table[prot_len_table$HGNC == gene, "protein.length"][1]
    if(is.na(len) || len < max(counts$aa_pos, na.rm=TRUE)) len <- max(counts$aa_pos, na.rm=TRUE) + 10
    
    y_max <- max_n; y_limit <- y_max * 1.3
    rect_y_max <- 0; rect_y_min <- -(max(1, y_max * 0.08))
    plot_title <- bquote(italic(.(gene)) ~ " :[Cohort M-Rate: " ~ .(mut_rate_str) ~ "]")
    
    p <- ggplot(counts, aes(x = aa_pos, y = N, color = Variant_Classification)) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        annotate("rect", xmin = 0, xmax = len, ymin = rect_y_min, ymax = rect_y_max, fill = "#8F9CA3", color = "black", linewidth=0.5) +
        annotate("text", x = len/2, y = rect_y_min/2, label = paste(gene, "Domain"), size = 3.5, fontface = "italic", color="black") +
        geom_segment(aes(x = aa_pos, xend = aa_pos, y = 0, yend = N, alpha = line_alpha), color = "gray50", linewidth=0.5) +
        geom_point(aes(size = pt_size, alpha = pt_alpha)) +
        geom_text_repel(data = top_hits, aes(label = Protein_Change), size = 3, fontface="bold", direction = "y", nudge_y = y_max * 0.1, box.padding = 0.5, segment.size = 0.4, segment.color="gray40", min.segment.length = 0, show.legend = FALSE) +
        scale_x_continuous(limits = c(0, len), expand = c(0.01, 0.01)) + scale_y_continuous(limits = c(rect_y_min * 1.5, y_limit), breaks = scales::pretty_breaks(n = 4)) +
        scale_color_manual(values = custom_colors) + scale_size_identity() + scale_alpha_identity() + theme_classic() +
        theme(legend.position = "bottom", legend.title = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_line(color="black"), plot.title = element_text(hjust = 0, size=14), panel.grid.major.y = element_line(color="gray95", linetype="dashed")) +
        labs(x = "Amino Acid Position", y = "Mutation Count", title = plot_title)
    return(p)
}

# =================模块：全基因组统合棒棒糖引擎=================
draw_wholegenome_lollipop <- function(maf_dt, custom_colors) {
    muts <- copy(maf_dt[!is.na(Start_Position)])
    target_class <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Silent", "Nonstop_Mutation", "Translation_Start_Site")
    muts <- muts[Variant_Classification %in% target_class]
    if(nrow(muts) == 0) return(NULL)
    
    # 提取所有基因区块的基因组范围（底座地基）
    gene_bounds <- muts[Hugo_Symbol != "Unknown" & Hugo_Symbol != "NA", .(start = min(Start_Position, na.rm=T), end = max(Start_Position, na.rm=T)), by = Hugo_Symbol]
    
    counts <- muts[, .N, by = .(Start_Position, Variant_Classification, Hugo_Symbol, Protein_Change)]
    max_n <- max(counts$N, na.rm = TRUE)
    
    counts[, pt_size := ifelse(Variant_Classification == "Silent", 1.0, 2.0 + (N/max_n) * 2.0)]
    counts[, pt_alpha := ifelse(Variant_Classification == "Silent", 0.2, 0.4 + (N/max_n) * 0.6)]
    counts[, line_alpha := ifelse(Variant_Classification == "Silent", 0.1, 0.3 + (N/max_n) * 0.4)]
    
    setorder(counts, -N)
    top_hits <- head(counts[Variant_Classification != "Silent"], 12)
    # 为全基因组坐标标记，带上基因名前缀
    top_hits[, label_txt := paste0(Hugo_Symbol, ": ", Protein_Change)]
    top_hits[is.na(Protein_Change) | Protein_Change=="NA", label_txt := paste0(Hugo_Symbol, ": ", Start_Position, "bp")]
    
    y_max <- max_n; y_limit <- y_max * 1.3
    rect_y_max <- 0; rect_y_min <- -(max(1, y_max * 0.08))
    
    p <- ggplot() +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        # 批量绘制所有基因地基区块
        geom_rect(data = gene_bounds, aes(xmin = start, xmax = end, ymin = rect_y_min, ymax = rect_y_max), fill = "#8F9CA3", color = "black", linewidth=0.5, alpha = 0.8) +
        geom_text(data = gene_bounds, aes(x = (start+end)/2, y = rect_y_min/2, label = Hugo_Symbol), size = 3.5, fontface = "bold.italic", angle=15, color="black") +
        # 整体突变标枪
        geom_segment(data = counts, aes(x = Start_Position, xend = Start_Position, y = 0, yend = N, color = Variant_Classification, alpha = line_alpha), linewidth=0.4) +
        geom_point(data = counts, aes(x = Start_Position, y = N, color = Variant_Classification, size = pt_size, alpha = pt_alpha)) +
        # 牵引出全基因组里最耀眼的 12 个突变王
        geom_text_repel(data = top_hits, aes(x = Start_Position, y = N, label = label_txt), size = 3.5, fontface="bold", direction = "y", nudge_y = y_max * 0.1, box.padding = 0.5, segment.size = 0.4, segment.color="gray40", min.segment.length = 0) +
        scale_color_manual(values = custom_colors) + scale_size_identity() + scale_alpha_identity() + theme_classic() +
        theme(legend.position = "bottom", legend.title = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_line(color="black"), plot.title = element_text(hjust = 0.5, size=16, face="bold"), panel.grid.major.y = element_line(color="gray95", linetype="dashed")) +
        labs(x = "Whole Genomic Position (bp)", y = "Mutation Count in Cohort", title = "Holo-Genome Polyprotein Mutation Hubs")
    return(p)
}

# =================模块：数据聚合与安全清洗=================
mega_maf_dt <- NULL
if (dir.exists(opt$input)) {
  maf_files <- list.files(path = opt$input, pattern = "\\.maf$", full.names = TRUE, recursive = TRUE)
  if (length(maf_files) == 0) stop("\n[异常] 未找到任何 .maf 文件。")
  cat(sprintf("\n[探针] 发现 %d 个标本数据，开启装填...\n", length(maf_files)))
  dt_list <- list()
  for (f in maf_files) {
    tmp_dt <- tryCatch({ data.table::fread(file = f, colClasses = "character", stringsAsFactors = FALSE, showProgress = FALSE) }, error = function(e) NULL)
    if (!is.null(tmp_dt) && nrow(tmp_dt) > 0) dt_list[[basename(f)]] <- tmp_dt
  }
  mega_maf_dt <- data.table::rbindlist(dt_list, fill = TRUE)
} else if (file.exists(opt$input)) {
  mega_maf_dt <- data.table::fread(file = opt$input, colClasses = "character", stringsAsFactors = FALSE, showProgress = FALSE)
}

suppressWarnings({ mega_maf_dt[, Start_Position := as.numeric(Start_Position)]; mega_maf_dt[, End_Position := as.numeric(End_Position)] })
true_accession <- mega_maf_dt$NCBI_Build[1]
mega_maf_dt[, NCBI_Build := "Viral_Consensus"]; mega_maf_dt[, Center := "GCVA_Consensus"]

mega_maf_dt[Reference_Allele == "TRUE", Reference_Allele := "T"]; mega_maf_dt[Tumor_Seq_Allele1 == "TRUE", Tumor_Seq_Allele1 := "T"]; mega_maf_dt[Tumor_Seq_Allele2 == "TRUE", Tumor_Seq_Allele2 := "T"]
mega_maf_dt[Reference_Allele == "FALSE", Reference_Allele := "F"]; mega_maf_dt[Tumor_Seq_Allele1 == "FALSE", Tumor_Seq_Allele1 := "F"]; mega_maf_dt[Tumor_Seq_Allele2 == "FALSE", Tumor_Seq_Allele2 := "F"]

mega_maf_dt[, Variant_Classification := sub("&.*", "", Variant_Classification)]
mega_maf_dt[, Variant_Classification := data.table::fcase( Variant_Classification == "missense_variant", "Missense_Mutation", Variant_Classification == "frameshift_variant" & Variant_Type == "DEL", "Frame_Shift_Del", Variant_Classification == "frameshift_variant" & Variant_Type == "INS", "Frame_Shift_Ins", Variant_Classification == "frameshift_variant", "Frame_Shift_Del", Variant_Classification == "stop_gained", "Nonsense_Mutation", Variant_Classification == "stop_lost", "Nonstop_Mutation", Variant_Classification == "start_lost", "Translation_Start_Site", Variant_Classification == "inframe_insertion", "In_Frame_Ins", Variant_Classification == "inframe_deletion", "In_Frame_Del", Variant_Classification %in% c("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant"), "Splice_Site", Variant_Classification %in% c("protein_altering_variant", "coding_sequence_variant"), "Missense_Mutation", Variant_Classification %in% c("disruptive_inframe_deletion", "conservative_inframe_deletion"), "In_Frame_Del", Variant_Classification %in% c("disruptive_inframe_insertion", "conservative_inframe_insertion"), "In_Frame_Ins", Variant_Classification == "synonymous_variant", "Silent", Variant_Classification == "upstream_gene_variant", "5'Flank", Variant_Classification == "downstream_gene_variant", "3'Flank", Variant_Classification == "intergenic_region", "IGR", Variant_Classification == "3_prime_UTR_variant", "3'UTR", Variant_Classification == "5_prime_UTR_variant", "5'UTR", Variant_Classification == "intragenic_variant", "Intron", default = "Targeted_Region" )]

if (opt$min_mut > 0) mega_maf_dt <- mega_maf_dt[Tumor_Sample_Barcode %in% mega_maf_dt[, .N, by = Tumor_Sample_Barcode][N >= opt$min_mut, Tumor_Sample_Barcode]]

custom_prot_dat <- get_ncbi_lengths(true_accession)
if(is.null(custom_prot_dat)) custom_prot_dat <- get_empirical_lengths(mega_maf_dt)

cat("[预处理完成] 交送统计中 ...\n")
viral_maf <- tryCatch({ read.maf(maf = mega_maf_dt, vc_nonSyn = tcga_native_non_syn, verbose = FALSE) }, error = function(e) stop("\n[挂载失败] ", e$message))

summary_df <- as.data.frame(viral_maf@summary)
n_samples <- as.numeric(summary_df[summary_df$ID == "Samples", "summary"])
if(is.na(n_samples) || length(n_samples)==0) n_samples <- 1
cat(sprintf("   => 项目矩阵装甲达成 | N=%d | 有效突变总量=%d\n", n_samples, nrow(viral_maf@data)))

# ================== 渲染流水线 ==================
safe_plot <- function(filepath, width, height, call_expr) {
  pdf(filepath, width = width, height = height)
  res <- tryCatch({ suppressWarnings(eval(call_expr)); TRUE }, error = function(e) { cat("   -[拦截]: ", e$message, "\n"); FALSE })
  while(dev.cur() > 1) dev.off(); if(!res && file.exists(filepath)) file.remove(filepath)
  # PNG version
  png_path <- sub("\\.pdf$", ".png", filepath)
  png(png_path, width = width, height = height, units = "in", res = 300)
  res2 <- tryCatch({ suppressWarnings(eval(call_expr)); TRUE }, error = function(e) { cat("   -[PNG拦截]: ", e$message, "\n"); FALSE })
  while(dev.cur() > 1) dev.off(); if(!res2 && file.exists(png_path)) file.remove(png_path)
}

cat("\n================ 基因空间绘图起跑 ================\n")
cat(" -> [保留经典模型] 绘制瀑布图与突变面板...\n")
safe_plot(paste0(output_prefix, "_01_mafSummary_TCGA.pdf"), 11, 8, quote({ plotmafSummary(maf = viral_maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, color = custom_colors, textSize = 0.8) }))
safe_plot(paste0(output_prefix, "_06_TiTv_Summary.pdf"), 10, 7, quote({ plotTiTv(res = titv(maf = viral_maf, plot = FALSE, useSyn = TRUE)) }))
safe_plot(paste0(output_prefix, "_02_Oncoplot.pdf"), max(10, min(14, n_samples * 1.5)), if(n_samples>1) max(8, min(15, 0.4*opt$top)) else 8, quote({ oncoplot(maf=viral_maf, top=opt$top, fontSize=0.8, colors=custom_colors, showTumorSampleBarcodes=(n_samples>1 && n_samples<50), titleText = sprintf("Mutational Cohort N=%d", n_samples)) }))

gene_sum <- getGeneSummary(viral_maf)$Hugo_Symbol
target_genes <- if(is.null(opt$gene)) head(gene_sum, 5) else if(tolower(opt$gene)=="all") gene_sum else intersect(unlist(strsplit(opt$gene, ",")), gene_sum)

lollipops_list <- list()
if (length(target_genes) > 0) {
    for (g in target_genes) {
        cat(" -> 独立解析: 生成 [", g, "] 棒棒糖图 (已修复 NA%)\n")
        p <- draw_viral_lollipop(mega_maf_dt, n_samples, g, custom_prot_dat, custom_colors)
        if(!is.null(p)) {
            lollipops_list[[g]] <- p
            safe_plot(paste0(output_prefix, sprintf("_03_Lollipop_%s.pdf", g)), 10, 5, quote({ print(p) }))
        }
    }
}

# 【新增功能 1】：无缝拼接所有的棒棒糖图
if (length(lollipops_list) > 1) {
    cat(" -> 矩阵融合: 正在使用 Patchwork 将所有基因无缝拼接为长图鉴...\n")
    safe_plot(paste0(output_prefix, "_03b_Lollipop_All_Stitched.pdf"), 12, 4 * length(lollipops_list), quote({
        # 将列表全部拼图，一列展示，并且统一收集所有图底下的图例放在最下面
        print(wrap_plots(lollipops_list, ncol = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom"))
    }))
}

# 【新增功能 2】：史诗级的全基因组唯一坐标系总图
cat(" -> 时空展开: 刻画真正的全基因组坐标系棒棒糖图矩阵 ...\n")
holo_p <- draw_wholegenome_lollipop(mega_maf_dt, custom_colors)
if(!is.null(holo_p)) {
    safe_plot(paste0(output_prefix, "_03c_WholeGenome_Lollipop_Map.pdf"), 16, 6, quote({ print(holo_p) }))
}

cat("\n[全线杀青] 无懈可击的出版级图库已构建完毕！\n\n")
