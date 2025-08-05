##If you want to build a costume GO database
if (FALSE) {
  library(clusterProfiler)
  library(openxlsx)
  library(dplyr)
  library(stringr)
  library(AnnotationForge)
  options(stringsAsFactors = F)
  #########读取蛋白组注释文件#########
  egg_f  <-  "E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Herbicide_analysis/Dsan_herbicide_DE/Dsan_eggNOG2.emapper.annotations.xlsx"
  egg  <-  read.xlsx(egg_f, sep = "\t")
  egg[egg==""]<-NA#将空行变成NA以便去除
  
  #########提取蛋白名称(Query)与egg注释信息(seed_ortholog)#########
  gene_info <- egg %>% dplyr::select(GID = query, EGGannot = seed_ortholog) %>% na.omit()
  #########提取蛋白名称(Query)与GO注释信息(GO_terms)#########
  gterms <- egg %>% dplyr::select(GID = query,GO_terms = GOs) %>% na.omit()
  gene2go <- data.frame(GID = character(),#由于此时存在一蛋白对应多个GO，因此将其拆成一对一的多列储存进新的dataframe中
                        GO = character(),
                        EVIDENCE = character())
  for (row in 1:nrow(gterms)) {
    gene_terms <- str_split(gterms[row,"GO_terms"], ",", simplify = FALSE)[[1]]  
    gene_id <- gterms[row, "GID"][[1]]
    tmp <- data.frame(GID = rep(gene_id, length(gene_terms)),
                      GO = gene_terms,
                      EVIDENCE = rep("IEA", length(gene_terms)))
    gene2go <- rbind(gene2go, tmp)}  
  #########建库保存#########
  gene2go <- gene2go[gene2go$GO!='-', ]
  genus = "Digitaria sanguinalis"
  species = "Ds"
  makeOrgPackage(go=gene2go,gene_info=gene_info, 
                 version="0.1",#指定该GO富集注释数据库的版本号
                 maintainer = 'Y. J. H. <yujiehuang@zju.edu.com>',#指定该数据库维护者信息
                 author = 'Y. J. H. <yujiehuang@zju.edu.com>',#指定该数据库作者信息
                 outputDir = "E:/Bio_analysis/Database", #数据库的保存路径
                 tax_id = "121769", #查询物种的Taxonomy https://www.ncbi.nlm.nih.gov/taxonomy
                 genus = genus,
                 species = species,
                 goTable = "go")
  # zymo_ZM4.orgdb <- str_c("org.", str_to_upper(str_sub(genus, 1, 1)) , species, ".eg.db", sep = "")#将数据库命名为genus首字母+species全称
}
parse_ratio <- function(ratio_str) {
parts <- as.numeric(unlist(strsplit(ratio_str, "/")))
  if (length(parts) == 2 && !any(is.na(parts))) {
    return(parts[1] / parts[2])
  } else {
    warning("Invalid ratio format")
    return(NA)
  }
}
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(dplyr)
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA

# install.packages("E:/Bio_analysis/Database/org.DDs.eg.db", repos=NULL, type="sources")#安装自建库
library(org.DDs.eg.db)#读取自建库
GO_database <- 'org.DDs.eg.db'

all_files <- list.files("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/Significant", recursive = TRUE, full.names = TRUE)
all_result <- data.frame()
for (signi_file in all_files) {
  cur_file <- read.table(signi_file,header = FALSE)
  genes <- paste0(cur_file$V1, '.mRNA1')

  ###################################GO富集分析###################################
  GO <- enrichGO(genes,
                OrgDb = GO_database,
                keyType = "GID", #设定读取的gene ID类型
                ont = "BP", #(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                # pvalueCutoff = 0.05, #设定p值阈值
                minGSSize = 1,
                #qvalueCutoff = 0.05,
                # pAdjustMethod = 'BH',
                readable = FALSE)

  df <- as.data.frame(GO)[,c('Description','GeneRatio','BgRatio','FoldEnrichment', 'qvalue','p.adjust','Count')]
        # %>%
        # arrange(p.adjust) %>%
        # head(5)
  if (nrow(df) == 0) {
    next
  }
  df$group <- basename(signi_file)
  all_result <- rbind(all_result, df)
  # df$GeneRatio <- sapply(df$GeneRatio, parse_ratio)
  # ggplot(df, aes(x = GeneRatio, y = reorder(str_wrap(Description, width = 20), FoldEnrichment, decreasing=TRUE))) +
  #   geom_point(aes(size = FoldEnrichment, color = -log10(p.adjust)),alpha = 0.8) +
  #   scale_size_continuous(name = "Fold Enrichment",range = c(3, 8), limits = c(5,60),breaks = scales::pretty_breaks(n = 4)) +
  #   scale_color_gradient(name = "-log10(p.adjust)",low = "blue",high = "red",guide = guide_colorbar(reverse = TRUE)) +
  #   # labs(x = "Fold Enrichment",y = "GO Term",title = "GO Enrichment Analysis") +
  #   theme_bw()
}
