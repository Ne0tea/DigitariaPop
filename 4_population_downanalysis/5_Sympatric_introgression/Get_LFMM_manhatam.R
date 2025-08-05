suppressPackageStartupMessages({
  library(optparse)
  library(qqman)
  library(qvalue)
})

process_files <- function(file_paths, ped_file){
  data_list <- lapply(file_paths, function(f){
    df <- read.table(f, header = TRUE,sep = ',')
    df[, 2]  # 提取第二列
  })
  p_matrix <- do.call(cbind, data_list)
  row_means <- rowMeans(p_matrix, na.rm = TRUE)
  
  base_df <- read.table(file_paths[1], header = TRUE)
  ped_name_df <- read.table(ped_file,header = FALSE)
  n_rows <- nrow(base_df)
  data.frame(
    CHR = ped_name_df[,1],
    SNP = paste(ped_name_df[,1], ped_name_df[,4], sep = "_"),
    BP = ped_name_df[,4],
    P = row_means,
    SNP_loc = paste(ped_name_df[,1], ped_name_df[,4], sep = "_")
  )
}

option_list <- list(
  make_option(c("-f", "--files"), type = "character",
              help = "输入文件路径", metavar = "FILE..."),
  make_option(c("-o", "--output"), type = "character",
              default = "result.csv",
              help = "输出文件路径 [默认 %default]"),
  make_option(c("-m", "--map"), type = "character",
              default = "result.csv",
              help = "指定lfmm输入的map文件")
  )

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

if(is.null(args$files)){
  cat("错误：必须指定输入文件\n")
  print_help(parser)
  quit(status = 1)
}
# file_list <- unlist(strsplit(args$files, " "))
file_list <- Sys.glob(args$files)
valid_files <- file_list[file.exists(file_list)]
if(length(valid_files) == 0){
  cat("未找到有效输入文件\n")
  quit(status = 2)
}

ped_file <- args$map
pvalue_dataframe <- process_files(valid_files, ped_file)
q <- qvalue(pvalue_dataframe$P)
q_value <- qvalue(pvalue_dataframe$P,fdr.level=0.05)
pvalue_dataframe$qvalues <- q_value$qvalues

pdf(paste0(args$output,'.pdf'), width=6, height=4)
manhattan(pvalue_dataframe,col = c('#00afb9'),annotatePval=0.01,cex.axis = 1.2,cex.lab=1.4)
dev.off()

write.table(pvalue_dataframe, args$output, row.names = FALSE,quote = FALSE, sep = '\t')