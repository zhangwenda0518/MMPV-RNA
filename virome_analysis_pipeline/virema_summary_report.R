#!/usr/bin/env Rscript
# ==============================================================================
# ViReMa 顶刊多维可视化渲染引擎 (优化版)
# [优化]: +弧线图 +柱状图 +add_transparency +本地GB缓存 +空数据保护
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(dplyr)
  library(tidyr); library(stringr); library(circlize); library(httr)
})

option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default=".", help="Directory containing .bed/.bedgraph files"),
  make_option(c("-o", "--out_dir"), type="character", default="./ViReMa_Report", help="Output directory"),
  make_option(c("-g", "--auto_annotate"), action="store_true", default=FALSE, help="Auto-download GenBank files via NCBI"),
  make_option(c("-c", "--gb_cache"), type="character", default=NULL, help="Local GenBank cache directory (skip NCBI if found)")
)
opt <- parse_args(OptionParser(usage = "%prog [options]", option_list=option_list))

inp_dir <- opt$input_dir; out_dir <- opt$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ind_dir <- file.path(out_dir, "Virus_Specific_Plots"); dir.create(ind_dir, showWarnings = FALSE, recursive = TRUE)

message("==========================================================")
message("ViReMa Visualization Suite (Optimized)")
message("==========================================================")

# =======================================================
# 0. 工具函数
# =======================================================
safe_plot <- function(expr, fb_msg="") {
  tryCatch({ suppressWarnings(eval(expr, envir=parent.frame())); TRUE },
           error = function(e) { message("  [!] skip: ", fb_msg, " (", conditionMessage(e), ")"); FALSE })
}

add_transparency <- function(hex_col, alpha_val) {
  if(!all(grepl("^#[0-9A-Fa-f]{6}$", hex_col))) {
    hex_col <- sub("^#([0-9A-Fa-f]{6}).*", "#\\1", hex_col)
  }
  rgb_vals <- col2rgb(hex_col) / 255
  rgb(rgb_vals[1,], rgb_vals[2,], rgb_vals[3,], alpha_val)
}

getGradColor <- function(vals, max_v, c_low, c_high) {
    if(max_v <= 1) return(rep(c_high, length(vals)))
    pal <- colorRamp(c(c_low, c_high))
    norm_v <- (vals - 1) / (max_v - 1); norm_v[norm_v < 0] <- 0; norm_v[norm_v > 1] <- 1
    rgb_cols <- pal(norm_v)
    rgb(rgb_cols[,1], rgb_cols[,2], rgb_cols[,3], maxColorValue=255)
}

calc_coverage_fast <- function(df, max_len, step = 100){
  bins <- seq(0, max_len + step, by = step); cov_array <- numeric(length(bins))
  if(nrow(df)==0) return(data.frame(V1=bins, n=cov_array))
  for(i in 1:nrow(df)){
    s_bin <- max(1, floor(as.numeric(df$V2[i])/step)+1); e_bin <- max(1, floor(as.numeric(df$V3[i])/step)+1)
    if(s_bin > e_bin) { tmp <- s_bin; s_bin <- e_bin; e_bin <- tmp }
    cov_array[s_bin:e_bin] <- cov_array[s_bin:e_bin] + as.numeric(df$V5[i])
  }
  return(data.frame(V1 = bins, n = cov_array))
}

# =======================================================
# 1. 读取核心数据 (空数据优雅退出)
# =======================================================
bed_files <- list.files(inp_dir, pattern = "Recombination_Results\\.(bed|txt)$", full.names = T, ignore.case = T)
if(length(bed_files)==0) {
  message("No BED files found. Exiting gracefully.")
  write.csv(data.frame(Note="No recombination events detected"), file.path(out_dir,"Empty_Report.csv"), row.names=FALSE)
  quit(save="no", status=0)
}

filelist <- list()
for(i in seq_along(bed_files)){
  tmp <- tryCatch(read.table(bed_files[i], skip=1, sep="\t", stringsAsFactors=F, quote="", fill=T), error=function(e) NULL)
  if(!is.null(tmp) && nrow(tmp)>0) { tmp$File <- basename(bed_files[i]); filelist[[i]] <- tmp }
}
boundDF <- bind_rows(filelist) %>% filter(!is.na(V1) & !grepl("track",V1,ignore.case=T) & !is.na(V2) & !is.na(V3))
if(nrow(boundDF)==0) {
  message("No valid recombination events found. Exiting gracefully.")
  write.csv(data.frame(Note="All BED files parsed but no valid events"), file.path(out_dir,"Empty_Report.csv"), row.names=FALSE)
  quit(save="no", status=0)
}

boundDF$V2 <- as.numeric(boundDF$V2); boundDF$V3 <- as.numeric(boundDF$V3); boundDF$V5 <- as.numeric(boundDF$V5)
boundDF <- boundDF %>% group_by(V1, V2, V3, V4) %>% mutate(V11=n(), V5=sum(V5,na.rm=T)) %>% ungroup()
write.csv(boundDF, file.path(out_dir, "Aggregated_Sequence_Information.csv"), row.names = FALSE)

# =======================================================
# 2. 覆盖度 + 微小突变 (空数据保护)
# =======================================================
bg_files <- list.files(inp_dir, pattern = "\\.bedgraph$", full.names = T, ignore.case = T)
bgDF <- data.frame(); has_bg <- FALSE
if(length(bg_files)>0){
    bg_list <- list()
    for(i in seq_along(bg_files)){
        tmp <- tryCatch(read.table(bg_files[i], skip=1, sep="\t", stringsAsFactors=F, fill=T), error=function(e) NULL)
        if(!is.null(tmp) && nrow(tmp)>=4) {
            tmp <- tmp[,1:4]; colnames(tmp) <- c("Ref","Start","End","Value")
            tmp$Type <- ifelse(grepl("Conservation|Deletion",basename(bg_files[i]),ignore.case=T),"Deletion","Duplication")
            bg_list[[i]] <- tmp
        }
    }
    bgDF <- bind_rows(bg_list)
    if(nrow(bgDF)>0){
        bgDF$Start <- as.numeric(bgDF$Start); bgDF$End <- as.numeric(bgDF$End); bgDF$Value <- as.numeric(bgDF$Value)
        bgDF <- bgDF %>% group_by(Ref,Start,End,Type) %>% summarise(Value=sum(Value,na.rm=T),.groups="drop")
        has_bg <- TRUE
    }
}

micro_files <- list.files(inp_dir, pattern = "Micro.*\\.(bed|txt)$", full.names = T, ignore.case = T)
microDF <- data.frame(); has_micro <- FALSE
if(length(micro_files)>0) {
    mlist <- list()
    for(i in seq_along(micro_files)){
        tmp <- tryCatch(read.table(micro_files[i], skip=1, sep="\t", stringsAsFactors=F, fill=T), error=function(e) NULL)
        if(!is.null(tmp) && ncol(tmp)>=3){
            tmp_df <- tmp[,1:3]; colnames(tmp_df) <- c("Ref","Start","End")
            tmp_df$Sample <- basename(micro_files[i]); mlist[[i]] <- tmp_df
        }
    }
    microDF <- bind_rows(mlist) %>% filter(!grepl("track",Ref,ignore.case=T))
    if(nrow(microDF)>0) {
        microDF$Start <- as.numeric(microDF$Start); microDF$End <- as.numeric(microDF$End)
        microDF <- microDF %>% group_by(Ref,Start,End) %>% summarise(Prevalence=n_distinct(Sample),.groups="drop")
        has_micro <- TRUE
    }
}

# =======================================================
# 3. NCBI 注释 (本地 .gb 缓存)
# =======================================================
parse_genbank_features <- function(gb_text, ref_id) {
  gb_lines <- strsplit(gb_text, "\n")[[1]]
  res <- data.frame(RefSeq=character(), Start=numeric(), End=numeric(), Name=character(), stringsAsFactors=F)
  feat_regex <- "^ {5}(CDS|gene|5'UTR|3'UTR|mature_peptide)\\s+(.+)"
  attr_regex <- "^ {21}/(gene|product)=\"([^\"]+)\""
  curr_type <- NULL; curr_start <- NA; curr_end <- NA; curr_name <- NULL
  for (line in gb_lines) {
    if (grepl("^ORIGIN", line)) break
    if (grepl("^ {5}\\w", line)) {
      if (!is.null(curr_type) && !is.na(curr_start) && !is.na(curr_end))
        res <- rbind(res, data.frame(RefSeq=ref_id, Start=curr_start, End=curr_end,
                                      Name=ifelse(is.null(curr_name),curr_type,curr_name)))
      curr_type <- NULL; curr_start <- NA; curr_end <- NA; curr_name <- NULL
      if (grepl(feat_regex, line)) {
        matches <- regmatches(line, regexec(feat_regex, line))[[1]]; curr_type <- matches[2]
        nums <- as.numeric(unlist(regmatches(matches[3], gregexpr("\\d+", matches[3]))))
        if (length(nums)>=2) { curr_start <- min(nums); curr_end <- max(nums) }
      }
    } else if (!is.null(curr_type) && grepl(attr_regex, line)) {
      if (is.null(curr_name)) curr_name <- regmatches(line, regexec(attr_regex, line))[[1]][3]
    }
  }
  if (!is.null(curr_type) && !is.na(curr_start) && !is.na(curr_end))
    res <- rbind(res, data.frame(RefSeq=ref_id, Start=curr_start, End=curr_end,
                                  Name=ifelse(is.null(curr_name),curr_type,curr_name)))
  if(nrow(res)>0) {
    res <- res[!duplicated(res[,c("Start","End")]), ]
    npg_palette <- c("4DBBD5","00A087","3C5488","F39B7F","8491B4","91D1C2","DC0000","7E6148")
    res$Color <- npg_palette[(0:(nrow(res)-1)) %% length(npg_palette) + 1]
  }
  return(res)
}

anno_global <- NULL; virus_names_dict <- list()

if (opt$auto_annotate) {
  message("NCBI annotations (local .gb cache)...")
  gb_dir <- if (!is.null(opt$gb_cache)) opt$gb_cache else file.path(out_dir, ".gb_cache")
  dir.create(gb_dir, showWarnings=FALSE, recursive=TRUE)
  all_anno <- data.frame()

  for (ref in unique(boundDF$V1)) {
    acc_clean <- str_squish(strsplit(as.character(ref), " ")[[1]][1])
    gb_file <- file.path(gb_dir, paste0(acc_clean, ".gb"))
    gb_text <- NULL

    if (file.exists(gb_file) && file.info(gb_file)$size > 500) {
      gb_text <- paste(readLines(gb_file, warn=FALSE), collapse="\n")
    } else {
      req <- tryCatch(GET(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=gb&retmode=text", acc_clean), timeout(15)), error=function(e) NULL)
      if (!is.null(req) && status_code(req)==200) {
        gb_text <- content(req, "text", encoding="UTF-8")
        if (!grepl("Error|Item not found", substr(gb_text,1,100))) writeLines(gb_text, gb_file)
      }
    }

    if (!is.null(gb_text)) {
      gb_lines <- strsplit(gb_text, "\n")[[1]]
      for(line in gb_lines) if(grepl("^SOURCE",line)) { virus_names_dict[[ref]] <- str_squish(gsub("SOURCE","",line)); break }
      feat_df <- parse_genbank_features(gb_text, ref)
      if(nrow(feat_df)>0) all_anno <- rbind(all_anno, feat_df)
    }
  }
  if(nrow(all_anno)>0) anno_global <- all_anno
}

# =======================================================
# 4. 样本统计矩阵
# =======================================================
message("Building sample statistics...")
idDF <- boundDF %>% select(V1,V2,V3,V4,V5,V11,File) %>%
  mutate(id=paste0(V1,V2,V3,V4,collapse=",")) %>% group_by(id,File) %>% slice(1) %>% ungroup()
uniqueDF <- idDF %>% select(File, V11) %>% filter(V11==1) %>% group_by(File) %>% summarise(`Unique SV Events`=sum(V11))
filesDF <- boundDF %>% group_by(File) %>% slice(1) %>% ungroup() %>% bind_rows(filelist) %>%
  group_by(File) %>% summarise(`Total SV-Chimeric Reads`=sum(as.numeric(V5),na.rm=T), `Total DVG Jump Counts`=n())
individualFiles <- merge(filesDF, uniqueDF, by="File", all.x=TRUE); individualFiles[is.na(individualFiles)] <- 0
write.csv(individualFiles, file.path(out_dir, "Matrix_Sample_Wise_Recombination_Statistics.csv"), row.names=FALSE)

# =======================================================
# 5. 弧线图 + 柱状图 (from visualize.py)
# =======================================================
plot_arc_diagram <- function(df, prefix, genome_len=NULL) {
  if(nrow(df)==0) return()
  color_map <- c("Deletion"="#E74C3C","Duplication"="#2980B9","Back-Splice"="#8E44AD","Insertion"="#27AE60","Splice"="#E67E22")
  max_cnt <- max(df$V5, na.rm=TRUE)

  p <- ggplot() + geom_hline(yintercept=0, color="grey60", linewidth=0.5) +
    theme_minimal(base_size=12) + labs(x="Genome Position (nt)", y="Jump Distance",
      title=paste0("Recombination Arc Diagram\n", basename(prefix)))

  for(i in seq_len(nrow(df))) {
    s <- df$V2[i]; e <- df$V3[i]
    clr <- "grey50"
    for(k in names(color_map)) if(grepl(k, df$V4[i], ignore.case=TRUE)) { clr <- color_map[k]; break }
    alpha <- if(max_cnt<=1) 0.8 else max(0.15, 0.1 + 0.9*log1p(df$V5[i])/log1p(max_cnt))
    curv <- ifelse(s < e, 0.3, -0.3)
    p <- p + annotate("curve", x=s, xend=e, y=0, yend=0, curvature=curv, color=clr, alpha=alpha, linewidth=1.2)
  }
  x_lim <- if(!is.null(genome_len)) max(c(df$V2,df$V3,genome_len)) else max(c(df$V2,df$V3))
  p <- p + xlim(0, x_lim) + scale_color_manual(values=color_map) +
    theme(legend.position=c(0.85,0.85), legend.background=element_rect(fill="white",color="grey80"))
  ggsave(paste0(prefix,"_Arc_Diagram.png"), plot=p, width=14, height=6, dpi=300)
}

plot_top_events <- function(df, prefix, top_n=15) {
  if(nrow(df)==0) return()
  top <- df %>% arrange(desc(V5)) %>% head(top_n) %>%
    mutate(Label=paste0(V2,"->",V3,"\n(",V4,")"),
           Color=ifelse(grepl("Deletion",V4,ignore.case=TRUE),"#E74C3C",
                 ifelse(grepl("Duplication",V4,ignore.case=TRUE),"#2980B9","grey50")))

  p <- ggplot(top, aes(x=reorder(Label,V5), y=V5, fill=Color)) +
    geom_bar(stat="identity") + coord_flip() + scale_fill_identity() +
    geom_text(aes(label=V5), hjust=-0.2, size=3) +
    labs(x=NULL, y="Read Count", title=paste0("Top ",top_n," Recombination Events\n",basename(prefix))) +
    theme_minimal(base_size=12) + ylim(0, max(top$V5)*1.15)
  ggsave(paste0(prefix,"_Top_Events.png"), plot=p, width=12, height=6, dpi=300)
}

# =======================================================
# 5b. 断点核苷酸频率图 (from ViReMaShiny)
# =======================================================
plot_breakpoint_nt_freq <- function(df, prefix) {
  # Need V9 (donor_seq) and V10 (acceptor_seq) from ViReMa BED
  if (ncol(df) < 10) {
    message("  BED has <10 cols, skipping nucleotide frequency plot")
    return()
  }
  df_nt <- df %>% select(V2, V3, V4, V9, V10) %>% filter(!is.na(V9) & nchar(V9) > 0 & !is.na(V10) & nchar(V10) > 0)
  if (nrow(df_nt) == 0) return()

  # Donor site frequencies (V9 - sequence at donor)
  donor_pos <- df_nt %>% mutate(seq_clean = str_replace_all(V9, "\\|", "")) %>%
    mutate(chars = str_split(seq_clean, "")) %>% unnest(chars) %>%
    group_by(V2, V3, V4, V9) %>% mutate(pos = row_number() - 25) %>% ungroup() %>%
    count(chars, pos) %>% filter(chars %in% c("A", "C", "G", "T"))

  if (nrow(donor_pos) > 0) {
    p_donor <- ggplot(donor_pos, aes(x=pos, y=n, color=chars)) +
      geom_line(linewidth=0.8) + geom_vline(xintercept=0, linetype="dashed", color="grey50") +
      scale_color_manual(values=c(A="#E41A1C", C="#377EB8", G="#4DAF4A", T="#984EA3")) +
      ylab("Count") + xlab("Position relative to Donor (bp)") +
      theme_classic(base_size=14) + ggtitle(paste0("Nucleotide Frequency at Donor Sites\n", basename(prefix)))
    ggsave(paste0(prefix, "_Donor_NT_Freq.png"), plot=p_donor, width=10, height=5, dpi=300)
  }

  # Acceptor site frequencies (V10)
  accept_pos <- df_nt %>% mutate(seq_clean = str_replace_all(V10, "\\|", "")) %>%
    mutate(chars = str_split(seq_clean, "")) %>% unnest(chars) %>%
    group_by(V2, V3, V4, V10) %>% mutate(pos = row_number() - 25) %>% ungroup() %>%
    count(chars, pos) %>% filter(chars %in% c("A", "C", "G", "T"))

  if (nrow(accept_pos) > 0) {
    p_accept <- ggplot(accept_pos, aes(x=pos, y=n, color=chars)) +
      geom_line(linewidth=0.8) + geom_vline(xintercept=0, linetype="dashed", color="grey50") +
      scale_color_manual(values=c(A="#E41A1C", C="#377EB8", G="#4DAF4A", T="#984EA3")) +
      ylab("Count") + xlab("Position relative to Acceptor (bp)") +
      theme_classic(base_size=14) + ggtitle(paste0("Nucleotide Frequency at Acceptor Sites\n", basename(prefix)))
    ggsave(paste0(prefix, "_Acceptor_NT_Freq.png"), plot=p_accept, width=10, height=5, dpi=300)
  }
}

# =======================================================
# 6. 主绘图引擎 (单体病毒)
# =======================================================
plot_engine <- function(df, prefix_path, anno_df=NULL, virus_id="Global", virus_name="Multiple Viruses") {

  xLims <- df %>% group_by(V1) %>% summarise(max_pos=max(c(V2,V3),na.rm=T)) %>% as.data.frame()
  row.names(xLims) <- xLims$V1
  lim_matrix <- cbind(rep(0, nrow(xLims)), xLims$max_pos); row.names(lim_matrix) <- xLims$V1
  b1 <- df %>% select(V1, V2); b2 <- df %>% select(V1, V3)

  # CSV matrix
  report_df <- df %>%
    select(Virus_Acc=V1, Event_Type=V4, Donor_Start=V2, Acceptor_End=V3, Recombination_Reads=V5, Prevalence_in_Cohort=V11) %>%
    distinct(Donor_Start, Acceptor_End, .keep_all=TRUE) %>% arrange(desc(Prevalence_in_Cohort), desc(Recombination_Reads))
  write.csv(report_df, paste0(prefix_path,"_SV_Profile_Matrix.csv"), row.names=FALSE)

  # Arc + Top Events + NT Frequency
  safe_plot(quote(plot_arc_diagram(df, prefix_path)), "arc diagram")
  safe_plot(quote(plot_top_events(df, prefix_path)), "top events")
  safe_plot(quote(plot_breakpoint_nt_freq(df, prefix_path)), "nt frequency")

  # Scatterplot
  fill_limit <- c(1, max(df$V11, na.rm=TRUE))
  safe_plot(quote({
    p_scat <- ggplot(df, aes(x=V3, y=V2, color=V11)) + geom_point(aes(size=V5), alpha=0.7) +
      scale_size_continuous("Read Count", range=c(1,4)) + ylab("Donor (Start)") + xlab("Acceptor (End)") +
      theme_classic(base_size=14) + geom_rug(col=rgb(0.5,0.6,0.7,alpha=0.2)) + scale_color_viridis_c("Isolates", limits=fill_limit)
    ggsave(paste0(prefix_path,"_Scatterplot.pdf"), plot=p_scat, width=8, height=6)
    ggsave(paste0(prefix_path,"_Scatterplot.png"), plot=p_scat, width=8, height=6, dpi=300)
  }), "scatterplot")

  # Coverage
  mDel <- df %>% filter(V2 < V3) %>% distinct(V2, V3, V4, .keep_all=TRUE)
  if(nrow(mDel)>0) safe_plot(quote({
    p_del <- ggplot(calc_coverage_fast(mDel, max(mDel$V3,na.rm=T)), aes(x=V1, y=n)) +
      geom_area(fill="#3498DB", alpha=0.7) + geom_step() + ylab("Reads") + xlab("Position") +
      theme_classic(base_size=14) + ggtitle("Deletions Coverage")
    ggsave(paste0(prefix_path,"_Deletions.pdf"), plot=p_del, width=8, height=3)
    ggsave(paste0(prefix_path,"_Deletions.png"), plot=p_del, width=8, height=3, dpi=300)
  }))

  mDup <- df %>% filter(V2 > V3) %>% distinct(V2, V3, V4, .keep_all=TRUE)
  if(nrow(mDup)>0) safe_plot(quote({
    p_dup <- ggplot(calc_coverage_fast(mDup, max(mDup$V2,na.rm=T)), aes(x=V1, y=n)) +
      geom_area(fill="#E74C3C", alpha=0.7) + geom_step() + ylab("Reads") + xlab("Position") +
      theme_classic(base_size=14) + ggtitle("Duplications Coverage")
    ggsave(paste0(prefix_path,"_Duplications.pdf"), plot=p_dup, width=8, height=3)
    ggsave(paste0(prefix_path,"_Duplications.png"), plot=p_dup, width=8, height=3, dpi=300)
  }))

  # Circos
  n_samples_detected <- length(unique(df$File))
  genome_len <- max(lim_matrix[,2])
  if(genome_len <= 0 || nrow(b1) == 0) { message(sprintf("  Skip Circos: %s (no coords)", basename(prefix_path))); return() }

  line_widths <- ifelse(df$V5 > 1, log10(df$V5)*1.5, 0.5); line_widths[is.na(line_widths)] <- 0.5
  max_samples <- max(df$V11, na.rm=TRUE)

  df$EventColor <- "#000000"
  del_idx <- grepl("Deletion",df$V4,ignore.case=T); dup_idx <- grepl("Duplication",df$V4,ignore.case=T)
  ins_idx <- grepl("Ins",df$V4,ignore.case=T); oth_idx <- !(del_idx|dup_idx|ins_idx)
  if(any(del_idx)) df$EventColor[del_idx] <- getGradColor(df$V11[del_idx],max_samples,"#AED6F1","#154360")
  if(any(dup_idx)) df$EventColor[dup_idx] <- getGradColor(df$V11[dup_idx],max_samples,"#F5B7B1","#7B241C")
  if(any(ins_idx)) df$EventColor[ins_idx] <- getGradColor(df$V11[ins_idx],max_samples,"#ABEBC6","#145A32")
  if(any(oth_idx)) df$EventColor[oth_idx] <- getGradColor(df$V11[oth_idx],max_samples,"#F9E79F","#7D6608")

  dynamic_alpha <- 0.7 - (df$V11/max_samples)*0.4
  link_cols <- add_transparency(df$EventColor, dynamic_alpha)
  max_micro_prev <- if(has_micro && nrow(microDF)>0) max(microDF$Prevalence,na.rm=T) else 1

  sub_title <- sprintf("Virus: %s    |    Acc: %s    |    Identified in %d Samples",
                        virus_name, strsplit(virus_id," ")[[1]][1], n_samples_detected)

  draw_circos_core <- function() {
      par(mar=c(4,2,4,2), oma=c(2,0,1,0))
      circos.clear()
      circos.par("start.degree"=90, "gap.degree"=1, "cell.padding"=c(0,0,0,0), "points.overflow.warning"=F)
      circos.initialize(sectors=df$V1, xlim=lim_matrix)

      circos.track(ylim=c(0,1), track.height=0.05, bg.border=NA, panel.fun=function(x,y) {
          circos.axis(labels.cex=0.5, direction="outside",
                      major.at=seq(0,max(lim_matrix[CELL_META$sector.index,2]),by=max(1,max(lim_matrix[CELL_META$sector.index,2])%/%10)),
                      labels.facing="clockwise")
          if(has_micro && nrow(microDF)>0){
              s_micro <- microDF[microDF$Ref==CELL_META$sector.index,]
              if(nrow(s_micro)>0){
                  snp_colors <- getGradColor(s_micro$Prevalence,max_micro_prev,"#D5D8DC","#8E44AD")
                  circos.segments(x0=s_micro$Start,y0=0,x1=s_micro$End,y1=1,col=snp_colors,lwd=1)
              }
          }
      })

      circos.track(ylim=c(0,1), track.height=0.07, bg.border=NA, panel.fun=function(x,y) {
          if(!is.null(anno_df) && nrow(anno_df)>0){
              tgt_anno <- anno_df[anno_df$RefSeq==CELL_META$sector.index,]
              if(nrow(tgt_anno)>0){
                  for(i in 1:nrow(tgt_anno)){
                     circos.rect(xleft=tgt_anno$Start[i], ybottom=0.1, xright=tgt_anno$End[i], ytop=0.9,
                                 col=paste0("#",tgt_anno$Color[i],"B3"), border="white", lwd=0.5)
                     circos.text(x=(tgt_anno$Start[i]+tgt_anno$End[i])/2, y=0.5, labels=tgt_anno$Name[i],
                                 cex=0.55, facing="clockwise", niceFacing=T, col="white", font=2)
                  }
              }
          }
      })

      circos.track(ylim=c(-1,1), track.height=0.15, bg.border="grey95", panel.fun=function(x,y) {
          circos.lines(x=c(0,max(lim_matrix[CELL_META$sector.index,2])), y=c(0,0), lty=3, col="grey60")
          if(has_bg && nrow(bgDF)>0){
              s_bg <- bgDF[bgDF$Ref==CELL_META$sector.index,]
              if(nrow(s_bg)>0){
                  dup_data <- s_bg[s_bg$Type=="Duplication",]
                  if(nrow(dup_data)>0){
                      dup_max <- max(dup_data$Value,na.rm=T)
                      for(i in 1:nrow(dup_data)){
                          norm_val <- dup_data$Value[i]/dup_max
                          circos.rect(xleft=dup_data$Start[i], ybottom=0, xright=dup_data$End[i],
                                      ytop=norm_val, col=add_transparency("#C0392B",1-norm_val*0.9), border=NA)
                      }
                  }
                  del_data <- s_bg[s_bg$Type=="Deletion",]
                  if(nrow(del_data)>0){
                      del_max <- max(del_data$Value,na.rm=T)
                      for(i in 1:nrow(del_data)){
                          norm_val <- -(del_data$Value[i]/del_max)
                          circos.rect(xleft=del_data$Start[i], ybottom=norm_val, xright=del_data$End[i],
                                      ytop=0, col=add_transparency("#2980B9",1-abs(norm_val)*0.9), border=NA)
                      }
                  }
              }
          }
      })

      circos.genomicLink(b1, b2, col=link_cols, lwd=line_widths, arr.type="triangle", directional=1, border=NA)
      title(main="Recombination Mutational Landscape", cex.main=1.4, font.main=2, line=1.5)
      mtext(sub_title, side=3, line=0.2, cex=0.9, font=2, col="grey30")
      legend("bottomleft", inset=c(0,-0.05), title="Track Encoding",
             legend=c("Outer: MicroInDel","Mid: Duplication Peak","Mid: Deletion Peak","Inner: Deletion Jump","Inner: Duplication Jump"),
             fill=c("#8E44AD","#C0392B","#2980B9","#3498DB","#E74C3C"), border=NA, cex=0.75, bg="white", box.col=NA, xpd=TRUE)
      legend("bottomright", inset=c(0,-0.05), title="Event Prevalence",
             legend=c("Hotspot: Shared across samples","Rare: Singleton Event"),
             fill=c("#154360","#AED6F1"), border=NA, cex=0.75, bg="white", box.col=NA, xpd=TRUE)
  }

  safe_plot(quote({ pdf(paste0(prefix_path,"_Circos_4Track.pdf"),width=10,height=10); draw_circos_core(); dev.off() }), "Circos PDF")
  safe_plot(quote({ png(paste0(prefix_path,"_Circos_4Track.png"),width=3000,height=3000,res=300); draw_circos_core(); dev.off() }), "Circos PNG")
}

# =======================================================
# 7. 引擎启动
# =======================================================
message("Generating per-virus reports...")
unique_viruses <- unique(boundDF$V1)
for(v in unique_viruses) {
    sub_df <- boundDF %>% filter(V1==v)
    clean_v_name <- str_replace_all(v, "[^A-Za-z0-9_.-]", "_")
    v_full_name <- "Unannotated Virus Strain"
    if(exists("virus_names_dict") && !is.null(virus_names_dict[[v]])) v_full_name <- virus_names_dict[[v]]
    message(sprintf("   [%s] %s", clean_v_name, v_full_name))
    # Per-virus subdirectory: {Virus_Specific_Plots}/{Acc}/
    virus_dir <- file.path(ind_dir, clean_v_name)
    dir.create(virus_dir, showWarnings=FALSE, recursive=TRUE)
    plot_engine(sub_df, file.path(virus_dir, clean_v_name), anno_global, v, v_full_name)
}

message("==========================================================")
message("All reports -> ", out_dir)
message("==========================================================")
