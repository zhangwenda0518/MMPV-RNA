#!/usr/bin/env Rscript

# ==============================================================================
# 病毒分类软件结果整合分析脚本 v4.0
# 兼容: virus_classifier2.py 的 combined_taxonomy.tsv + 各工具独立输出
# 新增: genomad, metabuli, CAT, diamond_lca, contigtax, BASTA 解析器
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(VennDiagram)
  library(ggVennDiagram)
  library(grid)
  library(gridExtra)
  library(stringr)
  library(RColorBrewer)
  library(scales)
  library(cowplot)
  library(patchwork)
  library(parallelly)
})

TAX_LEVELS <- c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# virus_classifier2.py 的输出列名 (小写)
TAX_V2_COLS <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

setup_threads <- function(cores = NULL) {
  avail <- availableCores()
  use <- if (is.null(cores)) max(2, min(avail - 1, 8)) else min(cores, avail)
  setDTthreads(use)
  cat(sprintf("[System] %d 核心, 使用 %d\n", avail, use))
}

log_msg <- function(level, msg, ...) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time(), "%H:%M:%S"), level, sprintf(msg, ...)))
}

# ==============================================================================
# 标准化
# ==============================================================================
standardize_dt <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  if (!"contig_id" %in% names(dt)) {
    id_guess <- grep("^(contig_id|seq_name|query|genome)$", names(dt), ignore.case=TRUE, value=TRUE)
    if (length(id_guess) > 0) setnames(dt, id_guess[1], "contig_id") else return(NULL)
  }
  dt[, contig_id := as.character(contig_id)]
  for (col in TAX_LEVELS) {
    if (!col %in% names(dt)) dt[, (col) := NA_character_]
    dt[, (col) := as.character(get(col))]
    dt[get(col) %in% c("", "-", "NA", "na", "N/A", "no rank", "undefined", "unknown", "null", "default", "Unclassified")
       | grepl("^Unclassified\\.", get(col))
       | grepl("^unplaced", get(col)), (col) := NA]
  }
  return(dt[, .SD, .SDcols = c("contig_id", TAX_LEVELS)])
}

# ==============================================================================
# 新增解析器: virus_classifier2.py combined_taxonomy.tsv
# ==============================================================================
parse_combined_v2 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 combined_taxonomy: %s", basename(file))
  dt <- fread(file, sep="\t")
  if (nrow(dt) == 0) return(NULL)

  # 按 tool 列拆分
  setnames(dt, names(dt), tolower(names(dt)))
  if (!"tool" %in% names(dt) || !"seq_name" %in% names(dt)) {
    log_msg("WARN", "combined_taxonomy 缺少 tool/seq_name 列")
    return(NULL)
  }
  setnames(dt, "seq_name", "contig_id")

  # v2 列名映射到标准 TAX_LEVELS (已用 tolower, 用小写匹配)
  col_map <- c(
    "superkingdom"="Realm", "realm"="Realm", "kingdom"="Kingdom",
    "phylum"="Phylum", "class"="Class", "order"="Order",
    "family"="Family", "genus"="Genus", "species"="Species"
  )
  for (v2col in names(col_map)) {
    if (v2col %in% names(dt)) setnames(dt, v2col, col_map[v2col])
  }

  # 按 tool 拆分为列表
  tool_list <- split(dt, by="tool")
  result <- lapply(names(tool_list), function(tn) {
    standardize_dt(tool_list[[tn]])
  })
  names(result) <- names(tool_list)
  return(result)
}

# ==============================================================================
# 各工具独立解析器 (兼容旧版)
# ==============================================================================
parse_mmseqs <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 MMseqs2: %s", basename(file))
  dt <- fread(file, header=FALSE, sep="\t", select=c(1,3,4,9),
              col.names=c("contig_id","rank","name","full_taxonomy"), fill=TRUE)
  dt <- dt[!rank %in% c("acellular root","root")]
  dt <- dt[grepl("Viruses|viria", full_taxonomy, ignore.case=TRUE)]
  if (nrow(dt)==0) return(NULL)

  dt[, (TAX_LEVELS) := NA_character_]
  parts <- tstrsplit(dt$full_taxonomy, ";", fixed=TRUE)
  rank_patterns <- list("^[kd]_"="Kingdom", "^p_"="Phylum", "^c_"="Class",
                        "^o_"="Order", "^f_"="Family", "^g_"="Genus", "^s_"="Species")
  valid_realms <- c("Adnaviria","Duplodnaviria","Monodnaviria","Riboviria","Ribozyviria","Varidnaviria")

  for (col_vec in parts) {
    col_vec <- trimws(col_vec)
    for (pat in names(rank_patterns)) {
      idx <- grep(pat, col_vec); if(length(idx)>0) dt[idx, (rank_patterns[[pat]]) := sub(pat,"",col_vec[idx])]
    }
    idx_r <- grep("^r_", col_vec); if(length(idx_r)>0) dt[idx_r, Realm := sub("^r_","",col_vec[idx_r])]
    idx_dash <- grep("^-_", col_vec)
    if (length(idx_dash)>0) {
      potentials <- sub("^-_","",col_vec[idx_dash])
      is_valid <- potentials %in% valid_realms
      rows <- idx_dash[is_valid]
      if (length(rows)>0) { rows2 <- rows[is.na(dt$Realm[rows])]; if(length(rows2)>0) dt[rows2, Realm := sub("^-_","",col_vec[rows2])] }
    }
  }
  rank_map <- c("species"="Species","genus"="Genus","family"="Family","order"="Order","class"="Class","phylum"="Phylum","kingdom"="Kingdom")
  for (r in names(rank_map)) dt[rank==r & !is.na(name) & name!="", (rank_map[[r]]) := name]
  return(standardize_dt(dt))
}

parse_vcontact3 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 vContact3: %s", basename(file))
  dt <- fread(file); if (nrow(dt)==0) return(NULL)
  setnames(dt, tolower(names(dt)))
  if ("reference" %in% names(dt)) dt <- dt[grepl("false", tolower(reference))]
  id_col <- grep("genome|bin", names(dt), value=TRUE)[1]; if(is.na(id_col)) return(NULL)
  setnames(dt, id_col, "contig_id"); dt <- dt[contig_id!="default"]
  for (lvl in TAX_LEVELS) {
    target <- grep(paste0("^",tolower(lvl),"_prediction"), names(dt), value=TRUE)[1]
    if (!is.na(target)) dt[, (lvl) := get(target)]
  }
  return(standardize_dt(dt))
}

parse_vitap <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 VITAP: %s", basename(file))
  hd <- readLines(file, n=20); skip_n <- grep("^Genome_ID", hd)[1]
  skip_n <- if(is.na(skip_n)) 0 else skip_n-1
  dt <- fread(file, skip=skip_n, fill=TRUE); if (nrow(dt)==0) return(NULL)
  setnames(dt, 1:2, c("contig_id", "lineage"))
  parts <- tstrsplit(dt$lineage, ";", fixed=TRUE)
  dt[, (TAX_LEVELS) := NA_character_]
  if (length(parts)>=1) dt[, Species := parts[[1]]]
  if (length(parts)>=2) dt[, Genus := parts[[2]]]
  if (length(parts)>=3) dt[, Family := parts[[3]]]
  if (length(parts)>=4) dt[, Order := parts[[4]]]
  if (length(parts)>=5) dt[, Class := parts[[5]]]
  if (length(parts)>=6) dt[, Phylum := parts[[6]]]
  if (length(parts)>=7) dt[, Kingdom := parts[[7]]]
  if (length(parts)>=8) dt[, Realm := parts[[8]]]
  return(standardize_dt(dt))
}

parse_phagcn3 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 PhaGCN3: %s", basename(file))
  dt <- fread(file, sep=",", fill=TRUE); if (nrow(dt)==0) return(NULL)
  if (!"contig_name" %in% names(dt)) setnames(dt, 1, "contig_name")
  setnames(dt, "contig_name", "contig_id")
  dt[, (TAX_LEVELS) := NA_character_]
  if ("prediction" %in% names(dt)) dt[, Family := prediction]
  char_cols <- names(dt)[sapply(dt, is.character)]
  char_cols <- setdiff(char_cols, c("contig_id", "prediction"))
  patterns <- list("r;"="Realm","k;"="Kingdom","p;"="Phylum","c;"="Class","o;"="Order","f;"="Family","g;"="Genus","s;"="Species")
  for (col in char_cols) {
    vv <- dt[[col]]; if (!any(grepl(";", vv[1:min(100,length(vv))]))) next
    for (pat in names(patterns)) {
      idx <- grep(pat, vv, fixed=TRUE)
      if (length(idx)>0) dt[idx, (patterns[[pat]]) := sub(paste0(".*",pat),"",vv[idx])]
    }
  }
  return(standardize_dt(dt))
}

parse_acvirus <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 ACVirus: %s", basename(file))
  dt <- fread(file, sep="\t", fill=TRUE); if (nrow(dt)==0) return(NULL)
  names(dt) <- paste0(toupper(substr(names(dt),1,1)), substr(names(dt),2,nchar(names(dt))))
  setnames(dt, 1, "contig_id")
  return(standardize_dt(dt))
}

# ==============================================================================
# 可视化 (保留原有全部函数)
# ==============================================================================
plot_classification_counts <- function(data_list, output_dir) {
  log_msg("VIS", "绘制分类数目条形图...")
  counts <- data.frame(Tool=names(data_list), Count=sapply(data_list, nrow))
  counts$Percentage <- round(counts$Count/sum(counts$Count)*100, 1)
  p <- ggplot(counts, aes(x=reorder(Tool,-Count), y=Count, fill=Tool)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%d\n(%.1f%%)", Count, Percentage)), vjust=-0.5, size=4) +
    labs(title="各工具分类数目", x="工具", y="序列数") +
    theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
    scale_fill_brewer(palette="Set2") + ylim(0, max(counts$Count)*1.15)
  ggsave(file.path(output_dir, "classification_counts.pdf"), p, width=max(8, nrow(counts)*1.2), height=6)
  return(counts)
}

plot_venn_diagrams <- function(data_list, output_dir, max_sets=5) {
  if (length(data_list) < 2) return()
  # >5 tools: 按计数取 top N
  if (length(data_list) > max_sets) {
    cnts <- sapply(data_list, nrow)
    keep <- names(sort(cnts, decreasing=TRUE)[1:max_sets])
    data_list <- data_list[keep]
    log_msg("INFO", "Venn图: 超过%d个工具, 只绘制top %d", max_sets, max_sets)
  }
  log_msg("VIS", "绘制 Venn 图...")
  id_lists <- lapply(data_list, function(x) unique(x$contig_id))
  colors <- brewer.pal(min(length(data_list), 8), "Set2")[1:length(data_list)]
  tryCatch({
    venn.plot <- venn.diagram(x=id_lists, filename=NULL, category.names=names(data_list),
      fill=colors, cex=1.2, cat.cex=1.2, cat.fontface="bold",
      main="病毒分类结果 Venn 图")
    pdf(file.path(output_dir, "venn_diagram.pdf"), width=10, height=10)
    grid.draw(venn.plot); dev.off()
  }, error=function(e) log_msg("WARN", "Venn: %s", e$message))
}

plot_venn_by_tax_level <- function(data_list, output_dir) {
  log_msg("VIS", "绘制各分类等级 Venn 图...")
  venn_plots <- list()
  for (level in TAX_LEVELS) {
    id_lists <- list()
    for (tn in names(data_list)) {
      ids <- data_list[[tn]][!is.na(get(level)) & get(level)!="", unique(contig_id)]
      if (length(ids)>0) id_lists[[tn]] <- ids
    }
    if (length(id_lists)>=2) {
      p <- ggVennDiagram(id_lists, label_alpha=0) +
        scale_fill_gradient(low="#F4FAFE", high="#4981BF") +
        labs(title=level) + theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.position="none")
      venn_plots[[level]] <- p
    }
  }
  if (length(venn_plots)>0) {
    combined <- wrap_plots(venn_plots, ncol=4, nrow=2) + plot_annotation(title="各分类等级 Venn 图")
    ggsave(file.path(output_dir, "combined_tax_level_venn.pdf"), combined, width=20, height=12)
  }
}

analyze_consistency_optimized <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return(NULL)
  log_msg("VIS", "分析一致性...")
  fll <- lapply(names(data_list), function(nm) {
    dt <- data_list[[nm]][contig_id %in% common_ids, .SD, .SDcols=c("contig_id",TAX_LEVELS)]; dt[, Tool:=nm]; dt
  })
  big_dt <- rbindlist(fll)
  results <- list()
  for (level in TAX_LEVELS) {
    vd <- big_dt[, .(contig_id, Tool, Taxon=get(level))][!is.na(Taxon) & Taxon!=""]
    if (nrow(vd)==0) next
    counts <- vd[, .N, by=.(contig_id, Taxon)]; setorder(counts, contig_id, -N, Taxon)
    cons <- unique(counts, by="contig_id")[, .(contig_id, Consensus=Taxon)]
    md <- merge(vd, cons, by="contig_id")
    md[, is_agreed := (!is.na(Taxon) & Taxon==Consensus)]
    stats <- md[!is.na(Taxon), .(Agreement_Count=sum(is_agreed), Total_Classified=.N), by=Tool]
    if (nrow(stats)>0) {
      stats[, `:=`(Taxonomic_Level=level, Agreement_Rate=Agreement_Count/Total_Classified)]
      results[[level]] <- stats
    }
  }
  ad <- rbindlist(results)
  if (nrow(ad)>0) {
    ad$Taxonomic_Level <- factor(ad$Taxonomic_Level, levels=TAX_LEVELS)
    p <- ggplot(ad, aes(x=Taxonomic_Level, y=Agreement_Rate, fill=Tool)) +
      geom_bar(stat="identity", position=position_dodge(0.9)) +
      geom_text(aes(label=sprintf("%.2f", Agreement_Rate)), position=position_dodge(0.9), vjust=-0.5, size=3) +
      labs(title="各工具与共识的一致率", x="分类等级", y="一致率") +
      theme_minimal() + scale_fill_brewer(palette="Set2") + scale_y_continuous(labels=percent, limits=c(0,1.1))
    ggsave(file.path(output_dir, "agreement_rates.pdf"), p, width=12, height=8)
    fwrite(ad, file.path(output_dir, "agreement_stats.tsv"), sep="\t")
  }
  return(ad)
}

analyze_detailed_consistency <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return(NULL)
  log_msg("INFO", "详细一致性分析...")
  tool_names <- names(data_list)
  consistency_summary <- data.frame()
  for (level in TAX_LEVELS) {
    long_list <- lapply(tool_names, function(tn) {
      dt <- data_list[[tn]][contig_id %in% common_ids, .(contig_id, taxon=get(level))]; dt[, Tool:=tn]
      dt[taxon %in% c("","NA","na","N/A","no rank","undefined","unknown","null"), taxon := NA_character_]; dt
    })
    big_dt <- rbindlist(long_list); valid_dt <- big_dt[!is.na(taxon)]
    if (nrow(valid_dt)>0) {
      counts <- valid_dt[, .N, by=.(contig_id, taxon)]; setorder(counts, contig_id, -N, taxon)
      consensus_dt <- counts[, .SD[1], by=contig_id, .SDcols="taxon"]; setnames(consensus_dt, "taxon", "consensus_taxon")
      stats_dt <- valid_dt[, .(uniq_cnt=uniqueN(taxon), all_vals_str=paste(sort(unique(taxon)), collapse="|")), by=contig_id]
      res_dt <- merge(stats_dt, consensus_dt, by="contig_id", all=TRUE)
      res_dt[, `:=`(is_consistent=(uniq_cnt==1), consistency_status=ifelse(uniq_cnt==1,"consistent","inconsistent"))]
    } else {
      res_dt <- data.table(contig_id=common_ids, uniq_cnt=0, all_vals_str=NA_character_,
                           consensus_taxon=NA_character_, is_consistent=NA, consistency_status="all_NA")
    }
    wide_dt <- dcast(big_dt, contig_id ~ Tool, value.var="taxon")
    final_comp <- merge(wide_dt, res_dt, by="contig_id", all.x=TRUE)
    consistent_out <- final_comp[is_consistent==TRUE, .(contig_id, consensus_taxon)]
    if (nrow(consistent_out)>0) fwrite(consistent_out, file.path(output_dir, paste0("consistent_",tolower(level),".tsv")), sep="\t")
    cols_order <- c("contig_id", tool_names, "all_vals_str", "consensus_taxon", "is_consistent")
    fwrite(final_comp[, ..cols_order], file.path(output_dir, paste0("comparison_",tolower(level),".tsv")), sep="\t")
    for (tn in tool_names) {
      if (tn %in% names(final_comp)) {
        cnt <- sum(!is.na(final_comp[[tn]]) & final_comp[[tn]]!="")
        consistency_summary <- rbind(consistency_summary, data.frame(Tool=tn, Level=level, Classified=cnt))
      }
    }
  }
  if (nrow(consistency_summary)>0) {
    p1 <- ggplot(consistency_summary, aes(x=Level, y=Classified, fill=Tool)) +
      geom_bar(stat="identity", position=position_dodge(0.8)) +
      labs(title="各工具各等级分类数", x="等级", y="序列数") + theme_minimal() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_fill_brewer(palette="Set2")
    p2 <- ggplot(consistency_summary, aes(x=Tool, y=Classified, fill=Level)) +
      geom_bar(stat="identity", position=position_stack()) +
      labs(title="分类数堆叠", x="工具", y="序列数") + theme_minimal() + scale_fill_brewer(palette="Set3")
    ggsave(file.path(output_dir, "classification_by_level.pdf"), p1, width=12, height=8)
    ggsave(file.path(output_dir, "classification_stacked.pdf"), p2, width=10, height=8)
    fwrite(consistency_summary, file.path(output_dir, "consistency_summary.tsv"), sep="\t")
  }
}

generate_grouped_consistent_files <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return()
  log_msg("INFO", "按分类单元分组...")
  tool_names <- names(data_list)
  overall_summary <- data.frame()
  for (level in TAX_LEVELS) {
    long_list <- lapply(tool_names, function(tn) {
      dt <- data_list[[tn]][contig_id %in% common_ids, .(contig_id, taxon=get(level))]
      dt[taxon %in% c("","NA","na","N/A"), taxon:=NA_character_]; dt
    })
    valid_dt <- rbindlist(long_list)[!is.na(taxon)]
    if (nrow(valid_dt)==0) next
    check_dt <- valid_dt[, .(uniq_cnt=uniqueN(taxon)), by=contig_id]
    consistent_ids <- check_dt[uniq_cnt==1, contig_id]
    if (length(consistent_ids)==0) next
    consistent_data <- valid_dt[contig_id %in% consistent_ids, .SD[1], by=contig_id]
    setnames(consistent_data, "taxon", "consensus_taxon")
    level_dir <- file.path(output_dir, paste0("consistent_",tolower(level)))
    dir.create(level_dir, recursive=TRUE, showWarnings=FALSE)
    summary_dt <- consistent_data[, .(Count=.N), by=consensus_taxon]; setorder(summary_dt, -Count)
    fwrite(summary_dt, file.path(level_dir, "summary.tsv"), sep="\t")
    for (taxon_name in unique(consistent_data$consensus_taxon)) {
      safe <- gsub("[^[:alnum:]_]", "_", taxon_name); safe <- gsub("_{2,}","_", gsub("^_|_$","", safe))
      sub_dt <- consistent_data[consensus_taxon==taxon_name, .(contig_id)]
      fwrite(sub_dt, file.path(level_dir, paste0(safe,".tsv")), sep="\t")
    }
    overall_summary <- rbind(overall_summary, data.frame(
      Taxonomic_Level=level, Unique_Taxa=length(unique(consistent_data$consensus_taxon)),
      Total_Sequences=nrow(consistent_data)))
  }
  if (nrow(overall_summary)>0) {
    fwrite(overall_summary, file.path(output_dir, "consistent_summary.tsv"), sep="\t")
    p1 <- ggplot(overall_summary, aes(x=Taxonomic_Level, y=Unique_Taxa, fill=Taxonomic_Level)) +
      geom_bar(stat="identity") + geom_text(aes(label=Unique_Taxa), vjust=-0.5) +
      labs(title="一致分类单元数", x="等级") + theme_minimal() + theme(legend.position="none") + scale_fill_brewer(palette="Set2")
    p2 <- ggplot(overall_summary, aes(x=Taxonomic_Level, y=Total_Sequences, fill=Taxonomic_Level)) +
      geom_bar(stat="identity") + geom_text(aes(label=Total_Sequences), vjust=-0.5) +
      labs(title="一致序列数", x="等级") + theme_minimal() + theme(legend.position="none") + scale_fill_brewer(palette="Set3")
    ggsave(file.path(output_dir, "consistent_summary.pdf"), plot_grid(p1, p2, ncol=2, labels="AUTO"), width=14, height=7)
  }
}

# ==============================================================================
# 主流程
# ==============================================================================
option_list <- list(
  make_option("--combined", type="character", default=NULL, help="virus_classifier2.py 的 combined_taxonomy.tsv"),
  make_option("--mmseqs", type="character", default=NULL),
  make_option("--vcontact3", type="character", default=NULL),
  make_option("--vitap", type="character", default=NULL),
  make_option("--phagcn3", type="character", default=NULL),
  make_option("--acvirus", type="character", default=NULL),
  make_option("--genomad", type="character", default=NULL, help="genomad taxonomy.tsv"),
  make_option("--metabuli", type="character", default=NULL, help="metabuli taxonomy.tsv"),
  make_option("--cat", type="character", default=NULL, help="CAT taxonomy.tsv"),
  make_option("--diamond_lca", type="character", default=NULL),
  make_option("--contigtax", type="character", default=NULL),
  make_option("--basta", type="character", default=NULL),
  make_option(c("-o","--output"), type="character", default="."),
  make_option("--cores", type="integer", default=NULL)
)

opt <- parse_args(OptionParser(option_list=option_list))
if (!dir.exists(opt$output)) dir.create(opt$output, recursive=TRUE)
setup_threads(opt$cores)

# 优先从 combined_taxonomy.tsv 解析所有工具
all_data <- list()
if (!is.null(opt$combined) && file.exists(opt$combined)) {
  combined_result <- parse_combined_v2(opt$combined)
  if (!is.null(combined_result)) all_data <- combined_result
}

# 补充独立工具文件 (覆盖/追加)
parsers <- list(
  mmseqs=parse_mmseqs, vcontact3=parse_vcontact3, vitap=parse_vitap,
  phagcn3=parse_phagcn3, acvirus=parse_acvirus
)
for (tool in names(parsers)) {
  f <- opt[[tool]]
  if (!is.null(f) && file.exists(f)) {
    tryCatch({
      res <- parsers[[tool]](f)
      if (!is.null(res) && nrow(res)>0) all_data[[tool]] <- res
    }, error=function(e) log_msg("ERROR", "%s: %s", tool, e$message))
  }
}

if (length(all_data) < 1) stop("无有效数据, 检查 --combined 或 --mmseqs/--vitap 等参数。")
log_msg("INFO", "共加载 %d 个工具数据: %s", length(all_data), paste(names(all_data), collapse=", "))

# 输出标准化表
for (tn in names(all_data)) {
  fwrite(unique(all_data[[tn]], by="contig_id"), file.path(opt$output, paste0("standardized_", tn, ".tsv")), sep="\t", na="NA")
}

# 交叉分析
all_ids <- unique(unlist(lapply(all_data, function(d) d$contig_id)))
n_total <- length(all_ids)
cat(sprintf("\n总序列数: %d\n", n_total))

count_stats <- plot_classification_counts(all_data, opt$output)
mat <- matrix(0, nrow=n_total, ncol=length(all_data), dimnames=list(all_ids, names(all_data)))
for (tn in names(all_data)) mat[all_data[[tn]]$contig_id, tn] <- 1
mat_dt <- as.data.table(mat, keep.rownames="contig_id")
mat_dt[, count := rowSums(.SD), .SDcols=names(all_data)]
fwrite(mat_dt, file.path(opt$output, "intersection_matrix.tsv"), sep="\t")

common_ids <- rownames(mat)[rowSums(mat)==length(all_data)]
if (length(common_ids)>0) writeLines(common_ids, file.path(opt$output, "common_ids.txt"))

if (length(all_data)>=2) {
  plot_venn_diagrams(all_data, opt$output)
  plot_venn_by_tax_level(all_data, opt$output)
}
if (length(common_ids)>0) {
  analyze_consistency_optimized(all_data, common_ids, opt$output)
  analyze_detailed_consistency(all_data, common_ids, opt$output)
  generate_grouped_consistent_files(all_data, common_ids, opt$output)
}

# 最终整合 (priority-based, Species 优先)
log_msg("INFO", "最终整合...")
long_dt <- rbindlist(lapply(names(all_data), function(n) {
  d <- copy(all_data[[n]]); d[, Tool:=n]; d
}), use.names=TRUE, fill=TRUE)
long_dt[, Completeness := rowSums(!is.na(.SD) & .SD!=""), .SDcols=TAX_LEVELS]
# 有效种名判定: 非空 + 非 -inae/-idae 结尾 + 非 sp./cf./aff.
valid_sp <- function(x) {
  if (is.na(x) || x=="") return(FALSE)
  if (grepl("inae$|idae$|ales$|icetes$", x, ignore.case=TRUE)) return(FALSE)
  if (grepl("\\bsp\\.?\\b|\\bcf\\.?\\b|\\baff\\.?\\b", x, ignore.case=TRUE)) return(FALSE)
  return(TRUE)
}
long_dt[, HasSpecies := valid_sp(Species)]
prio_map <- c("vcontact3"=1,"vitap"=2,"acvirus"=3,"phagcn3"=4,"mmseqs"=5,
              "genomad"=6,"CAT"=7,"metabuli"=8,"diamond_lca"=9,"contigtax"=10,"BASTA"=11)
long_dt[, Priority := prio_map[Tool]]; long_dt[is.na(Priority), Priority := 99]
# 排序: 有合法种名的优先, 其次填满度高, 再其次工具优先级
setorder(long_dt, contig_id, -HasSpecies, -Completeness, Priority)
final_result <- unique(long_dt, by="contig_id")
tools_agg <- long_dt[, .(All_Tools=paste(sort(unique(Tool)), collapse=",")), by=contig_id]
final_result <- merge(final_result, tools_agg, by="contig_id", sort=FALSE)
final_result[, Determination_Method := paste0(Tool, " (Score:", Completeness, ", Tools:", All_Tools, ")")]
out_cols <- c("contig_id", TAX_LEVELS, "Determination_Method")
fwrite(final_result[, ..out_cols], file.path(opt$output, "final_integrated_classification.tsv"), sep="\t", na="NA")

# 报告
sink(file.path(opt$output, "analysis_summary.txt"))
cat("=== 病毒分类分析报告 ===\n\n")
cat(sprintf("时间: %s\n", Sys.time()))
cat(sprintf("总序列: %d\n", n_total))
cat(sprintf("最终分类: %d\n", nrow(final_result)))
cat(sprintf("工具数: %d\n", length(all_data)))
cat("\n各工具分类数:\n"); print(count_stats)
cat("\n来源分布:\n"); print(table(final_result$Tool))
sink()
log_msg("INFO", "完成! 结果: %s", opt$output)
