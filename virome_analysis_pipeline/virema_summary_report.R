#!/usr/bin/env Rscript
# ==============================================================================
# ViReMa 顶刊多维可视化渲染引擎 (抛弃鸡肋全局图·单体深耕矩阵版)
# [优化]: 砍掉无用的 Global All-Virus，专注于单病毒的极限解析
# [新增]: 为每一个分离出的病毒，配套生成属于它的重组断点清单与流行统计 CSV 表
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(dplyr)
  library(tidyr); library(stringr); library(circlize); library(httr)
})

option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default=".", help="Directory containing .bed/.bedgraph files"),
  make_option(c("-o", "--out_dir"), type="character", default="./ViReMa_Report", help="Output directory"),
  make_option(c("-g", "--auto_annotate"), action="store_true", default=FALSE, help="Auto-download GenBank files via NCBI")
)
opt <- parse_args(OptionParser(usage = "%prog [options]", option_list=option_list))

inp_dir <- opt$input_dir; out_dir <- opt$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ind_dir <- file.path(out_dir, "Virus_Specific_Plots"); dir.create(ind_dir, showWarnings = FALSE, recursive = TRUE)

message("==========================================================")
message("🚀 激活 ViReMa 终极画匠中枢 (单体深耕·纯净输出版)")
message("==========================================================")

# =======================================================
# 1. 读取核心数据 (Central Web)
# =======================================================
bed_files <- list.files(inp_dir, pattern = "Recombination_Results\\.(bed|txt)$", full.names = T, ignore.case = T)
if(length(bed_files)==0) stop("❌ 未发现结果文件，终止！")

filelist <- list()
for(i in seq_along(bed_files)){
  tmp <- tryCatch(read.table(bed_files[i], skip=1, sep="\t", stringsAsFactors=F, quote="", fill=T), error=function(e) NULL)
  if(!is.null(tmp) && nrow(tmp)>0) { tmp$File <- basename(bed_files[i]); filelist[[i]] <- tmp }
}
boundDF <- bind_rows(filelist) %>% filter(!is.na(V1) & !grepl("track", V1, ignore.case=T) & !is.na(V2) & !is.na(V3))
if(nrow(boundDF)==0) { message("⚠️ 未检测到有效事件！"); quit(save="no", status=0) }

boundDF$V2 <- as.numeric(boundDF$V2); boundDF$V3 <- as.numeric(boundDF$V3); boundDF$V5 <- as.numeric(boundDF$V5)
boundDF <- boundDF %>% group_by(V1, V2, V3, V4) %>% mutate(V11 = n(), V5 = sum(V5, na.rm=T)) %>% ungroup()
write.csv(boundDF, file.path(out_dir, "Aggregated_Sequence_Information.csv"), row.names = FALSE)

# =======================================================
# 2. 地形覆盖度与微小突变
# =======================================================
bg_files <- list.files(inp_dir, pattern = "\\.bedgraph$", full.names = T, ignore.case = T)
bgDF <- data.frame()
if(length(bg_files)>0){
    bg_list <- list()
    for(i in seq_along(bg_files)){
        tmp <- tryCatch(read.table(bg_files[i], skip=1, sep="\t", stringsAsFactors=F, fill=T), error=function(e) NULL)
        if(!is.null(tmp) && nrow(tmp)>=4) {
            tmp <- tmp[, 1:4]; colnames(tmp) <- c("Ref", "Start", "End", "Value")
            tmp$Type <- ifelse(grepl("Conservation|Deletion", basename(bg_files[i]), ignore.case=T), "Deletion", "Duplication")
            bg_list[[i]] <- tmp
        }
    }
    bgDF <- bind_rows(bg_list)
    if(nrow(bgDF)>0){
        bgDF$Start <- as.numeric(bgDF$Start); bgDF$End <- as.numeric(bgDF$End); bgDF$Value <- as.numeric(bgDF$Value)
        bgDF <- bgDF %>% group_by(Ref, Start, End, Type) %>% summarise(Value=sum(Value, na.rm=T), .groups="drop")
    }
}

micro_files <- list.files(inp_dir, pattern = "Micro.*\\.(bed|txt)$", full.names = T, ignore.case = T)
microDF <- data.frame()
if(length(micro_files)>0) {
    mlist <- list()
    for(i in seq_along(micro_files)){
        tmp <- tryCatch(read.table(micro_files[i], skip=1, sep="\t", stringsAsFactors=F, fill=T), error=function(e) NULL)
        if(!is.null(tmp) && ncol(tmp)>=3){ 
            tmp_df <- tmp[, 1:3]; colnames(tmp_df) <- c("Ref", "Start", "End")
            tmp_df$Sample <- basename(micro_files[i]); mlist[[i]] <- tmp_df
        }
    }
    microDF <- bind_rows(mlist) %>% filter(!grepl("track", Ref, ignore.case=T))
    if(nrow(microDF)>0) {
        microDF$Start <- as.numeric(microDF$Start); microDF$End <- as.numeric(microDF$End)
        microDF <- microDF %>% group_by(Ref, Start, End) %>% summarise(Prevalence = n_distinct(Sample), .groups="drop")
    }
}

# =======================================================
# 3. NCBI 外骨骼标注与“病毒识名器”
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
      if (!is.null(curr_type) && !is.na(curr_start) && !is.na(curr_end)) res <- rbind(res, data.frame(RefSeq=ref_id, Start=curr_start, End=curr_end, Name=ifelse(is.null(curr_name), curr_type, curr_name)))
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
  if (!is.null(curr_type) && !is.na(curr_start) && !is.na(curr_end)) res <- rbind(res, data.frame(RefSeq=ref_id, Start=curr_start, End=curr_end, Name=ifelse(is.null(curr_name), curr_type, curr_name)))
  if(nrow(res)>0) {
    res <- res[!duplicated(res[, c("Start", "End")]), ]
    npg_palette <- c("4DBBD5", "00A087", "3C5488", "F39B7F", "8491B4", "91D1C2", "DC0000", "7E6148")
    res$Color <- npg_palette[(0:(nrow(res)-1)) %% length(npg_palette) + 1]
  }
  return(res)
}

anno_global <- NULL
virus_names_dict <- list()

if (opt$auto_annotate) {
  message("\n🌐 后台连接 NCBI 获取基因骨架与物种学名...")
  all_anno <- data.frame()
  for (ref in unique(boundDF$V1)) {
    acc_clean <- str_squish(strsplit(as.character(ref), " ")[[1]][1])
    req <- tryCatch(GET(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=gb&retmode=text", acc_clean), timeout(8)), error = function(e) NULL)
    if (!is.null(req) && status_code(req) == 200) {
      gb_text <- content(req, "text", encoding="UTF-8")
      
      gb_lines <- strsplit(gb_text, "\n")[[1]]
      for(line in gb_lines) {
          if (grepl("^SOURCE", line)) {
              virus_names_dict[[ref]] <- str_squish(gsub("SOURCE", "", line))
              break
          }
      }
      
      feat_df <- parse_genbank_features(gb_text, ref) 
      if(nrow(feat_df)>0) all_anno <- rbind(all_anno, feat_df)
    } 
  }
  if (nrow(all_anno)>0) anno_global <- all_anno
}

# =======================================================
# 5. 生成极其详尽的核心指标大表与精读信息汇总
# =======================================================
message("📊 脱水完成，构筑全景数据矩阵字典...")
# 表格1: 每个独立样本里的总读段指标统计，适合放到论文附件的数据总览表
idDF <- boundDF %>% 
  select(V1, V2, V3, V4, V5, V11, File) %>% 
  mutate(id = paste0(V1, V2, V3, V4, collapse = ",")) %>% 
  group_by(id, File) %>% slice(1) %>% ungroup()
uniqueDF <- idDF %>% select(File, V11) %>% filter(V11 == 1) %>% group_by(File) %>% summarise(`Unique SV Events` = sum(V11))
filesDF <- boundDF %>% group_by(File) %>% slice(1) %>% ungroup() %>% bind_rows(filelist) %>% group_by(File) %>% summarise(`Total SV-Chimeric Reads` = sum(as.numeric(V5), na.rm=T), `Total DVG Jump Counts` = n())
individualFiles <- merge(filesDF, uniqueDF, by = "File", all.x = TRUE); individualFiles[is.na(individualFiles)] <- 0
write.csv(individualFiles, file.path(out_dir, "Matrix_Sample_Wise_Recombination_Statistics.csv"), row.names = FALSE)


# =======================================================
# 6. 单体专属透射测绘与出图引擎
# =======================================================

getGradColor <- function(vals, max_v, c_low, c_high) {
    if(max_v <= 1) return(rep(c_high, length(vals)))
    pal <- colorRamp(c(c_low, c_high))
    norm_v <- (vals - 1) / (max_v - 1); norm_v[norm_v < 0] <- 0; norm_v[norm_v > 1] <- 1
    rgb_cols <- pal(norm_v)
    rgb(rgb_cols[,1], rgb_cols[,2], rgb_cols[,3], maxColorValue=255)
}

calc_coverage_fast <- function(df, max_len, step = 100){
  bins <- seq(0, max_len + step, by = step); cov_array <- numeric(length(bins))
  for(i in 1:nrow(df)){
    s_bin <- max(1, floor(as.numeric(df$V2[i])/step)+1); e_bin <- max(1, floor(as.numeric(df$V3[i])/step)+1)
    if(s_bin > e_bin) { tmp <- s_bin; s_bin <- e_bin; e_bin <- tmp }
    cov_array[s_bin:e_bin] <- cov_array[s_bin:e_bin] + as.numeric(df$V5[i])
  }
  return(data.frame(V1 = bins, n = cov_array))
}

plot_engine <- function(df, prefix_path, anno_df=NULL, virus_id="Global", virus_name="Multiple Viruses Detected") {
  
  xLims <- df %>% group_by(V1) %>% summarise(max_pos = max(c(V2, V3), na.rm=T)) %>% as.data.frame()
  row.names(xLims) <- xLims$V1
  lim_matrix <- cbind(rep(0, nrow(xLims)), xLims$max_pos); row.names(lim_matrix) <- xLims$V1
  b1 <- df %>% select(V1, V2); b2 <- df %>% select(V1, V3)
  
  # ------ 输出单病毒专属大表 (Virus Specific Report Matrix) ------
  # 去重：整理出这个病毒基因组上，所有确切发生过重组的位置
  virus_specific_csv_path <- paste0(prefix_path, "_SV_Profile_Matrix.csv")
  report_df <- df %>% 
    select(Virus_Acc=V1, Event_Type=V4, Donor_Start=V2, Acceptor_End=V3, Recombination_Reads=V5, Prevalence_in_Cohort=V11) %>% 
    distinct(Donor_Start, Acceptor_End, .keep_all = TRUE) %>% 
    arrange(desc(Prevalence_in_Cohort), desc(Recombination_Reads))
  write.csv(report_df, virus_specific_csv_path, row.names = FALSE)

  # ------ 出图：两仪散点图与覆盖度山峰辅助视图 ------ 
  fill_limit <- c(1, max(df$V11, na.rm = TRUE))
  p_scat <- ggplot(df, aes(x = V3, y = V2, color = V11)) + geom_point(aes(size = V5), alpha = 0.7) +
    scale_size_continuous("Read Count", range = c(1, 4)) + ylab("Donor (Start)") + xlab("Acceptor (End)") + 
    theme_classic(base_size = 14) + geom_rug(col = rgb(0.5, 0.6, 0.7, alpha = 0.2)) + scale_color_viridis_c("Isolates", limits = fill_limit)
  ggsave(paste0(prefix_path, "_Scatterplot.pdf"), plot = p_scat, width = 8, height = 6)

  mDel <- df %>% filter(V2 < V3) %>% distinct(V2, V3, V4, .keep_all = TRUE)
  if(nrow(mDel) > 0){
    p_del <- ggplot(calc_coverage_fast(mDel, max(mDel$V3, na.rm=T)), aes(x = V1, y = n)) + geom_area(fill = "#3498DB", alpha=0.7) + geom_step() + ylab("Reads Coverage") + xlab("Position") + theme_classic(base_size=14) + ggtitle("Deletions Coverage")
    ggsave(paste0(prefix_path, "_Deletions.pdf"), plot=p_del, width=8, height=3)
  }
  
  mDup <- df %>% filter(V2 > V3) %>% distinct(V2, V3, V4, .keep_all = TRUE)
  if(nrow(mDup) > 0){
    p_dup <- ggplot(calc_coverage_fast(mDup, max(mDup$V2, na.rm=T)), aes(x = V1, y = n)) + geom_area(fill = "#E74C3C", alpha=0.7) + geom_step() + ylab("Reads Coverage") + xlab("Position") + theme_classic(base_size=14) + ggtitle("Duplications Coverage")
    ggsave(paste0(prefix_path, "_Duplications.pdf"), plot=p_dup, width=8, height=3)
  }

  # ------ 核心绘图：Circos 超维展示 ------
  line_widths <- ifelse(df$V5 > 1, log10(df$V5) * 1.5, 0.5); line_widths[is.na(line_widths)] <- 0.5
  max_samples <- max(df$V11, na.rm=T)
  n_samples_detected <- length(unique(df$File))
  
  df$EventColor <- "#000000" 
  del_idx <- grepl("Deletion", df$V4, ignore.case=T); dup_idx <- grepl("Duplication", df$V4, ignore.case=T)
  ins_idx <- grepl("Ins", df$V4, ignore.case=T); oth_idx <- !(del_idx | dup_idx | ins_idx)
  
  if(any(del_idx)) df$EventColor[del_idx] <- getGradColor(df$V11[del_idx], max_samples, "#AED6F1", "#154360") 
  if(any(dup_idx)) df$EventColor[dup_idx] <- getGradColor(df$V11[dup_idx], max_samples, "#F5B7B1", "#7B241C") 
  if(any(ins_idx)) df$EventColor[ins_idx] <- getGradColor(df$V11[ins_idx], max_samples, "#ABEBC6", "#145A32") 
  if(any(oth_idx)) df$EventColor[oth_idx] <- getGradColor(df$V11[oth_idx], max_samples, "#F9E79F", "#7D6608") 
  
  dynamic_alpha <- 0.7 - (df$V11 / max_samples) * 0.4
  link_cols <- add_transparency(df$EventColor, dynamic_alpha)
  max_micro_prev <- if(nrow(microDF)>0) max(microDF$Prevalence, na.rm=T) else 1

  sub_title <- sprintf("Virus: %s    |    Acc: %s    |    Identified in %d Samples", virus_name, strsplit(virus_id, " ")[[1]][1], n_samples_detected)

  draw_circos_core <- function() {
      par(mar = c(4, 2, 4, 2), oma = c(2, 0, 1, 0))
      circos.clear()
      circos.par("start.degree"=90, "gap.degree"=1, "cell.padding"=c(0,0,0,0), "points.overflow.warning"=F)
      circos.initialize(sectors = df$V1, xlim = lim_matrix)
      
      circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(x, y) {
          circos.axis(labels.cex = 0.5, direction = "outside", major.at=seq(0, max(lim_matrix[CELL_META$sector.index, 2]), by=max(lim_matrix[CELL_META$sector.index, 2])%/%10), labels.facing="clockwise")
          sector.index <- CELL_META$sector.index
          if(nrow(microDF) > 0){
              s_micro <- microDF[microDF$Ref == sector.index, ]
              if(nrow(s_micro) > 0){
                  snp_colors <- getGradColor(s_micro$Prevalence, max_micro_prev, "#D5D8DC", "#8E44AD")
                  circos.segments(x0 = s_micro$Start, y0 = 0, x1 = s_micro$End, y1 = 1, col = snp_colors, lwd = 1)
              }
          }
      })
      
      circos.track(ylim = c(0, 1), track.height = 0.07, bg.border = NA, panel.fun = function(x, y) {
          sector.index <- CELL_META$sector.index
          if(!is.null(anno_df)){
              tgt_anno <- anno_df[anno_df$RefSeq == sector.index, ]
              if(nrow(tgt_anno) > 0){
                  for(i in 1:nrow(tgt_anno)){
                     circos.rect(xleft=tgt_anno$Start[i], ybottom=0.1, xright=tgt_anno$End[i], ytop=0.9, 
                                 col=paste0("#",tgt_anno$Color[i], "B3"), border="white", lwd=0.5)
                     circos.text(x=(tgt_anno$Start[i]+tgt_anno$End[i])/2, y=0.5, labels=tgt_anno$Name[i], 
                                 cex=0.55, facing="clockwise", niceFacing=T, col="white", font=2)
                  }
              }
          }
      })

      circos.track(ylim = c(-1, 1), track.height = 0.15, bg.border = "grey95", panel.fun = function(x, y) {
          circos.lines(x = c(0, max(lim_matrix[CELL_META$sector.index, 2])), y=c(0,0), lty=3, col="grey60")
          sector.index <- CELL_META$sector.index
          if(nrow(bgDF) > 0){
              s_bg <- bgDF[bgDF$Ref == sector.index, ]
              if(nrow(s_bg) > 0){
                  dup_data <- s_bg[s_bg$Type == "Duplication", ]
                  if(nrow(dup_data) > 0){
                      dup_max <- max(dup_data$Value, na.rm=T)
                      for(i in 1:nrow(dup_data)){
                          norm_val <- dup_data$Value[i] / dup_max
                          peak_col <- add_transparency("#C0392B", 1 - norm_val * 0.9)
                          circos.rect(xleft=dup_data$Start[i], ybottom=0, xright=dup_data$End[i], ytop=norm_val, col=peak_col, border=NA)
                      }
                  }
                  del_data <- s_bg[s_bg$Type == "Deletion", ]
                  if(nrow(del_data) > 0){
                      del_max <- max(del_data$Value, na.rm=T)
                      for(i in 1:nrow(del_data)){
                          norm_val <- -(del_data$Value[i] / del_max)
                          peak_col <- add_transparency("#2980B9", 1 - abs(norm_val) * 0.9)
                          circos.rect(xleft=del_data$Start[i], ybottom=norm_val, xright=del_data$End[i], ytop=0, col=peak_col, border=NA)
                      }
                  }
              }
          }
      })

      circos.genomicLink(b1, b2, col = link_cols, lwd = line_widths, arr.type = "triangle", directional = 1, border = NA)
      
      title(main = "Recombination Mutational Landscape", cex.main=1.4, font.main=2, line=1.5)
      mtext(sub_title, side=3, line=0.2, cex=0.9, font=2, col="grey30")
      
      legend("bottomleft", inset=c(0, -0.05), title="Track Encoding",
             legend=c("Outer: SNP / MicroInDel", "Mid: Duplication Peak", "Mid: Deletion Peak", "Inner: Deletion Jump", "Inner: Duplication Jump"),
             fill=c("#8E44AD", "#C0392B", "#2980B9", "#3498DB", "#E74C3C"), border=NA, cex=0.75, bg="white", box.col=NA, xpd=TRUE)
      
      legend("bottomright", inset=c(0, -0.05), title="Event Prevalence (Color Intensity)",
             legend=c("Hotspot: Highly Shared across samples", "Rare: Singleton Event"),
             fill=c("#154360", "#AED6F1"), border=NA, cex=0.75, bg="white", box.col=NA, xpd=TRUE)
  }

  pdf(paste0(prefix_path, "_Circos_4Track.pdf"), width=10, height=10); draw_circos_core(); dev.off()
  png(paste0(prefix_path, "_Circos_4Track.png"), width=3000, height=3000, res=300); draw_circos_core(); dev.off()
}

# ------------------------------------------------------------------------------
# 引擎启动 
# ------------------------------------------------------------------------------
message("🔬 剥夺全局伪图，专注生成每一株病毒的极致雷达投影与重组数据库...")
unique_viruses <- unique(boundDF$V1)
for(v in unique_viruses) {
    sub_df <- boundDF %>% filter(V1 == v)
    clean_v_name <- str_replace_all(v, "[^A-Za-z0-9_.-]", "_")
    
    v_full_name <- "Unannotated Virus Strain"
    if(exists("virus_names_dict") && !is.null(virus_names_dict[[v]])) {
        v_full_name <- virus_names_dict[[v]]
    }
    message(sprintf("   -> [单体刻录中] > %s", clean_v_name))
    plot_engine(sub_df, file.path(ind_dir, clean_v_name), anno_global, v, v_full_name)
}

message("==========================================================")
message("✨ (Matrix Active) 极净产出完成！现已为每个靶向病毒附赠独立重组矩阵！")
message("==========================================================")
