#!/usr/bin/env Rscript
# ==============================================================================
# panvirome_suite.R — 泛病毒组全维图谱自动构建管线
# ==============================================================================
# Modules: Mutation Density | TiTv | Oncoplot | Lollipop | dN/dS | Sliding π
# Usage: Rscript panvirome_suite.R --merged-dir Merged_Results_All_Viruses/

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(zoo)
  library(ggpubr)
  library(fs)
  library(maftools)
})

option_list <- list(
  make_option(c("-d", "--merged-dir"), type="character", default="Merged_Results_All_Viruses",
              help="Aggregator 输出目录 [default: %default]"),
  make_option(c("-v", "--virus"), type="character", default=NULL,
              help="仅处理指定病毒 (默认: 全部)")
)
opt <- parse_args(OptionParser(option_list=option_list))

BASE_DIR <- opt[["merged-dir"]]
if (!dir.exists(BASE_DIR)) stop(sprintf("目录不存在: %s", BASE_DIR))

# TCGA 标准非同义定义
tcga_native_non_syn <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                         "Translation_Start_Site", "Frame_Shift_Ins", "Frame_Shift_Del",
                         "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Targeted_Region",
                         "5'Flank", "3'Flank", "IGR", "3'UTR", "5'UTR", "Intron")

custom_colors <- c("Missense_Mutation"="#33A02C", "Frame_Shift_Del"="#E31A1C",
                   "Frame_Shift_Ins"="#FF7F00", "Nonsense_Mutation"="#1F78B4",
                   "Nonstop_Mutation"="#A6CEE3", "Translation_Start_Site"="#000000",
                   "In_Frame_Ins"="#FDBF6F", "In_Frame_Del"="#CAB2D6", "Silent"="#CCCCCC",
                   "Splice_Site"="#6A3D9A", "5'Flank"="#B15928", "3'Flank"="#E7298A",
                   "IGR"="#E7298A", "3'UTR"="#B15928", "5'UTR"="#E7298A", "Intron"="#A6CEE3")

get_empirical_lengths <- function(maf_dt) {
  dt <- copy(maf_dt)
  dt[, aa_pos := suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", Protein_Change)))]
  gl <- dt[!is.na(aa_pos), .(max_pos = max(aa_pos, na.rm=TRUE)), by=Hugo_Symbol]
  if (nrow(gl) == 0) return(NULL)
  gl[, est_len := floor(max_pos * 1.1) + 10]
  data.frame(HGNC=gl$Hugo_Symbol, protein.length=gl$est_len, stringsAsFactors=FALSE)
}

draw_viral_lollipop <- function(maf_dt, n_samples, gene, prot_len_table) {
  muts <- copy(maf_dt[Hugo_Symbol == gene])
  target_class <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
                    "In_Frame_Del","In_Frame_Ins","Silent","Nonstop_Mutation","Translation_Start_Site")
  muts <- muts[Variant_Classification %in% target_class]
  if (nrow(muts) == 0) return(NULL)
  muts[, aa_pos := suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", Protein_Change)))]
  muts <- muts[!is.na(aa_pos)]
  if (nrow(muts) == 0) return(NULL)

  mutated_samples <- length(unique(muts$Tumor_Sample_Barcode))
  mut_rate_str <- sprintf("%.2f%%", (mutated_samples / n_samples) * 100)

  counts <- muts[, .N, by=.(aa_pos, Variant_Classification, Protein_Change)]
  max_n <- max(counts$N, na.rm=TRUE)
  counts[, pt_size := ifelse(Variant_Classification=="Silent", 1.0, 2.5 + (N/max_n)*2.5)]
  counts[, pt_alpha := ifelse(Variant_Classification=="Silent", 0.2, 0.4 + (N/max_n)*0.6)]
  counts[, line_alpha := ifelse(Variant_Classification=="Silent", 0.1, 0.3 + (N/max_n)*0.4)]
  setorder(counts, -N)
  top_hits <- head(counts[Variant_Classification != "Silent"], 10)

  len <- if (!is.null(prot_len_table) && gene %in% prot_len_table$HGNC)
    prot_len_table[prot_len_table$HGNC==gene, "protein.length"][1] else max(counts$aa_pos, na.rm=TRUE)+10

  y_limit <- max_n * 1.3
  rect_y_max <- 0; rect_y_min <- -(max(1, max_n * 0.08))

  ggplot(counts, aes(x=aa_pos, y=N, color=Variant_Classification)) +
    geom_hline(yintercept=0, color="black", linewidth=0.5) +
    annotate("rect", xmin=0, xmax=len, ymin=rect_y_min, ymax=rect_y_max,
             fill="#8F9CA3", color="black", linewidth=0.5) +
    annotate("text", x=len/2, y=rect_y_min/2, label=paste(gene,"Domain"),
             size=3.5, fontface="italic", color="black") +
    geom_segment(aes(x=aa_pos, xend=aa_pos, y=0, yend=N, alpha=line_alpha),
                 color="gray50", linewidth=0.5) +
    geom_point(aes(size=pt_size, alpha=pt_alpha)) +
    geom_text_repel(data=top_hits, aes(label=Protein_Change), size=3, fontface="bold",
                    direction="y", nudge_y=max_n*0.1, box.padding=0.5,
                    segment.color="gray40", min.segment.length=0, show.legend=FALSE) +
    scale_x_continuous(limits=c(0,len), expand=c(0.01,0.01)) +
    scale_y_continuous(limits=c(rect_y_min*1.5, y_limit), breaks=scales::pretty_breaks(n=4)) +
    scale_color_manual(values=custom_colors) + scale_size_identity() + scale_alpha_identity() +
    theme_classic() +
    theme(legend.position="bottom", legend.title=element_blank(),
          axis.line.x=element_blank(), plot.title=element_text(hjust=0, size=14),
          panel.grid.major.y=element_line(color="gray95", linetype="dashed")) +
    labs(x="Amino Acid Position", y="Mutation Count",
         title=bquote(italic(.(gene))~":[Cohort M-Rate:"~.(mut_rate_str)~"]"))
}

draw_wholegenome_lollipop <- function(maf_dt, abs_genes_df=NULL) {
  muts <- copy(maf_dt[!is.na(Start_Position)])
  target_class <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
                    "In_Frame_Del","In_Frame_Ins","Silent","Nonstop_Mutation","Translation_Start_Site")
  muts <- muts[Variant_Classification %in% target_class]
  if (nrow(muts) == 0) return(NULL)

  if (!is.null(abs_genes_df) && nrow(abs_genes_df) > 0) {
    gene_bounds <- data.frame(Hugo_Symbol=abs_genes_df$gene, start=abs_genes_df$start, end=abs_genes_df$end)
  } else {
    gene_bounds <- muts[Hugo_Symbol!="Unknown" & Hugo_Symbol!="NA",
                        .(start=min(Start_Position, na.rm=T), end=max(Start_Position, na.rm=T)), by=Hugo_Symbol]
  }

  counts <- muts[, .N, by=.(Start_Position, Variant_Classification, Hugo_Symbol, Protein_Change)]
  max_n <- max(counts$N, na.rm=TRUE)
  counts[, pt_size := ifelse(Variant_Classification=="Silent", 1.0, 2.0 + (N/max_n)*2.0)]
  counts[, pt_alpha := ifelse(Variant_Classification=="Silent", 0.2, 0.4 + (N/max_n)*0.6)]
  counts[, line_alpha := ifelse(Variant_Classification=="Silent", 0.1, 0.3 + (N/max_n)*0.4)]
  setorder(counts, -N)
  top_hits <- head(counts[Variant_Classification!="Silent"], 15)
  top_hits[, label_txt := paste0(Hugo_Symbol, ":", Protein_Change)]
  top_hits[is.na(Protein_Change)|Protein_Change=="NA", label_txt := paste0(Hugo_Symbol,":",Start_Position)]

  rect_y_max <- 0; rect_y_min <- -(max(1, max_n * 0.08))

  ggplot() +
    geom_hline(yintercept=0, color="black", linewidth=0.5) +
    geom_rect(data=gene_bounds, aes(xmin=start, xmax=end, ymin=rect_y_min, ymax=rect_y_max, fill=Hugo_Symbol),
              color="black", linewidth=0.5, alpha=0.6) +
    geom_text(data=gene_bounds, aes(x=(start+end)/2, y=rect_y_min/2, label=Hugo_Symbol),
              size=3, fontface="bold.italic", angle=0, color="black") +
    geom_segment(data=counts, aes(x=Start_Position, xend=Start_Position, y=0, yend=N,
                 color=Variant_Classification, alpha=line_alpha), linewidth=0.5) +
    geom_point(data=counts, aes(x=Start_Position, y=N, color=Variant_Classification,
               size=pt_size, alpha=pt_alpha)) +
    geom_text_repel(data=top_hits, aes(x=Start_Position, y=N, label=label_txt),
                    size=3, fontface="bold", direction="y", nudge_y=max_n*0.1,
                    box.padding=0.5, segment.color="gray40", min.segment.length=0) +
    scale_color_manual(values=custom_colors) + scale_size_identity() + scale_alpha_identity() +
    theme_classic() +
    theme(legend.position="bottom", legend.title=element_blank(),
          axis.line.x=element_blank(), plot.title=element_text(hjust=0.5, size=16, face="bold")) +
    guides(fill="none") +
    labs(x="Whole Genomic Position (bp)", y="Mutation Count in Cohort",
         title="Holo-Genome Architectural Polyprotein Mutation Hubs")
}

safe_plot <- function(filepath, width, height, call_expr) {
  pdf(filepath, width=width, height=height)
  res <- tryCatch({ suppressWarnings(eval(call_expr)); TRUE },
                  error=function(e) { cat("   -[skip]:", e$message, "\n"); FALSE })
  while (dev.cur() > 1) dev.off()
  if (!res && file.exists(filepath)) file.remove(filepath)
}

# ── 主循环 ──
virus_dirs <- fs::dir_ls(BASE_DIR, type="directory")
if (!is.null(opt[["virus"]])) {
  virus_dirs <- virus_dirs[basename(virus_dirs) == opt[["virus"]]]
}
cat(sprintf("共 %d 个病毒\n", length(virus_dirs)))

for (v_dir in virus_dirs) {
  virus_name <- basename(v_dir)
  cat(sprintf("\n▶ %s\n", virus_name))

  safe_read <- function(fp, ...) {
    if (file.exists(fp)) tryCatch(read.csv(fp, ...), error=function(e) NULL) else NULL
  }

  genes_df <- safe_read(file.path(v_dir, "gene_coordinates.csv"))
  is_viroid <- is.null(genes_df) || nrow(genes_df) == 0

  # ── I: 突变密度 ──
  df_ivar <- safe_read(file.path(v_dir, "merged_allele_freq.csv"))
  if (!is.null(df_ivar) && nrow(df_ivar) > 0) {
    p1 <- ggplot(df_ivar, aes(x=POS)) +
      geom_density(fill="steelblue", alpha=0.5, color="darkblue") +
      theme_classic() + labs(x="Genomic Position", y="Mutation Density", title="Global Mutational Landscape")
    if (!is_viroid) {
      y_est <- max(density(df_ivar$POS)$y, na.rm=TRUE)
      p1 <- p1 +
        geom_rect(data=genes_df, inherit.aes=FALSE,
                  aes(xmin=start, xmax=end, ymin=-y_est*0.05, ymax=-y_est*0.01, fill=gene), alpha=0.8) +
        geom_text(data=genes_df, inherit.aes=FALSE,
                  aes(x=(start+end)/2, y=-y_est*0.03, label=gene), size=3, color="white", fontface="bold") +
        theme(legend.position="none")
    }
    safe_plot(file.path(v_dir, "Plot1_Mutation_Density.pdf"), 10, 4, quote(print(p1)))
    cat("  Plot1 Density\n")
  }

  # ── II: Maftools ──
  if (!is_viroid) {
    maf_file <- file.path(v_dir, "merged_snpeff.maf")
    if (file.exists(maf_file) && file.info(maf_file)$size > 100) {
      cat("  MAF found, running maftools...\n")
      mega_maf_dt <- data.table::fread(file=maf_file, colClasses="character", showProgress=FALSE)
      if (nrow(mega_maf_dt) > 0) {
        suppressWarnings({
          mega_maf_dt[, Start_Position := as.numeric(Start_Position)]
          mega_maf_dt[, End_Position := as.numeric(End_Position)]
        })
        mega_maf_dt[Reference_Allele=="TRUE", Reference_Allele:="T"]
        mega_maf_dt[Tumor_Seq_Allele1=="TRUE", Tumor_Seq_Allele1:="T"]
        mega_maf_dt[Reference_Allele=="FALSE", Reference_Allele:="F"]
        mega_maf_dt[Tumor_Seq_Allele1=="FALSE", Tumor_Seq_Allele1:="F"]

        mega_maf_dt[, Variant_Classification := sub("&.*", "", Variant_Classification)]
        mega_maf_dt[, Variant_Classification := data.table::fcase(
          Variant_Classification=="missense_variant", "Missense_Mutation",
          Variant_Classification=="frameshift_variant" & Variant_Type=="DEL", "Frame_Shift_Del",
          Variant_Classification=="frameshift_variant" & Variant_Type=="INS", "Frame_Shift_Ins",
          Variant_Classification=="frameshift_variant", "Frame_Shift_Del",
          Variant_Classification=="stop_gained", "Nonsense_Mutation",
          Variant_Classification=="stop_lost", "Nonstop_Mutation",
          Variant_Classification=="start_lost", "Translation_Start_Site",
          Variant_Classification=="inframe_insertion", "In_Frame_Ins",
          Variant_Classification=="inframe_deletion", "In_Frame_Del",
          Variant_Classification %in% c("synonymous_variant","silent_variant"), "Silent",
          default="Missense_Mutation"
        )]

        viral_maf_obj <- tryCatch({
          read.maf(maf=mega_maf_dt, vc_nonSyn=tcga_native_non_syn, verbose=FALSE)
        }, error=function(e) NULL)

        if (!is.null(viral_maf_obj)) {
          n_samples <- as.numeric(as.data.frame(viral_maf_obj@summary)[1,"summary"])
          if (length(n_samples)==0 || is.na(n_samples)) n_samples <- 1

          safe_plot(file.path(v_dir,"Plot2_Maftools_Dashboard.pdf"), 11, 8,
                    quote(plotmafSummary(maf=viral_maf_obj, dashboard=TRUE, color=custom_colors)))
          safe_plot(file.path(v_dir,"Plot2_Maftools_TiTv.pdf"), 10, 7,
                    quote(plotTiTv(res=titv(maf=viral_maf_obj, plot=FALSE, useSyn=TRUE))))
          safe_plot(file.path(v_dir,"Plot2_Maftools_Oncoplot.pdf"),
                    max(10, min(14, n_samples*1.5)), 9,
                    quote(oncoplot(maf=viral_maf_obj, top=15, colors=custom_colors)))

          gene_levels <- getGeneSummary(viral_maf_obj)$Hugo_Symbol
          target_genes <- head(gene_levels, 6)
          emp_lens <- get_empirical_lengths(mega_maf_dt)
          lolli_list <- list()
          for (g in target_genes) {
            p_l <- draw_viral_lollipop(mega_maf_dt, n_samples, g, emp_lens)
            if (!is.null(p_l)) lolli_list[[g]] <- p_l
          }
          if (length(lolli_list) > 1) {
            safe_plot(file.path(v_dir,"Plot3_Lollipops_Stitched.pdf"), 12, max(4, 4*length(lolli_list)), quote({
              print(wrap_plots(lolli_list, ncol=1) + plot_layout(guides="collect") &
                      theme(legend.position="bottom"))
            }))
          }

          holo_p <- draw_wholegenome_lollipop(mega_maf_dt, abs_genes_df=genes_df)
          if (!is.null(holo_p))
            safe_plot(file.path(v_dir,"Plot3_WholeGenome_Lollipop.pdf"), 16, 6, quote(print(holo_p)))

          cat("  Plot2 Dashboard/TiTv/Oncoplot + Plot3 Lollipops\n")
        }
      }
    }
  }

  # ── III: 选择压力 ──
  if (!is_viroid) {
    df_prod <- safe_read(file.path(v_dir, "merged_product_results.csv"))
    if (!is.null(df_prod) && nrow(df_prod) > 0) {
      df_prod$piN <- as.numeric(df_prod$piN); df_prod$piS <- as.numeric(df_prod$piS)
      df_prod$Ratio <- df_prod$piN / (df_prod$piS + 1e-5)
      df_prod <- df_prod[!is.na(df_prod$Ratio) & df_prod$product %in% genes_df$gene,]
      if (nrow(df_prod) > 3) {
        p3 <- ggplot(df_prod, aes(x=product, y=Ratio, fill=product)) +
          geom_boxplot(outlier.shape=NA, alpha=0.7) + geom_jitter(width=0.2, alpha=0.5, color="darkgray") +
          geom_hline(yintercept=1, linetype="dashed", color="red") + theme_bw() +
          stat_compare_means(method="kruskal.test", label.y=max(df_prod$Ratio, na.rm=TRUE)*1.1) +
          labs(x="Viral Gene", y=expression(pi[N]/pi[S]~Ratio), title="Positive/Purifying Selection Profiling") +
          theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
        safe_plot(file.path(v_dir,"Plot4_Selection_Pressure.pdf"),
                  max(6, length(unique(df_prod$product))*1.2), 5, quote(print(p3)))
        cat("  Plot4 Selection\n")
      }
    }
  }

  # ── IV: 滑动 π ──
  df_site <- safe_read(file.path(v_dir, "merged_site_pi.csv"))
  if (!is.null(df_site) && nrow(df_site) > 50) {
    df_site <- df_site[order(df_site$site),]
    df_site$rolling_pi <- rollmean(df_site$pi, k=100, fill=NA, align="center")
    p4 <- ggplot() +
      geom_area(data=df_site, aes(x=site, y=rolling_pi), fill="indianred", alpha=0.6) +
      geom_line(data=df_site, aes(x=site, y=rolling_pi), color="darkred") +
      theme_classic() + labs(x="Genomic Position (bp)", y=expression(Mean~Diversity~(pi)),
                              title="Sliding-window Nucleotide Diversity (100bp)")
    if (!is_viroid && exists("genes_df") && nrow(genes_df) > 0) {
      my_pi <- max(df_site$rolling_pi, na.rm=TRUE)
      if (is.na(my_pi) || my_pi <= 0) my_pi <- 0.01
      p4 <- p4 +
        geom_rect(data=genes_df, inherit.aes=FALSE,
                  aes(xmin=start, xmax=end, ymin=-my_pi*0.05, ymax=-my_pi*0.01, fill=gene), alpha=0.8) +
        geom_text(data=genes_df, inherit.aes=FALSE,
                  aes(x=(start+end)/2, y=-my_pi*0.03, label=gene), size=3, color="white") +
        theme(legend.position="none")
    }
    safe_plot(file.path(v_dir, "Plot5_Diversity_Sliding_Window.pdf"), 10, 4, quote(print(p4)))
    cat("  Plot5 Sliding π\n")
  }
}

cat("\nDone.\n")
