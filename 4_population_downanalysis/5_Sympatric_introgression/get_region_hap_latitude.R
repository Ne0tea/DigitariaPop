library(ComplexHeatmap)
library(circlize)

args <- commandArgs(trailingOnly = TRUE)  # 获取命令行参数

###Import file 
prefix_012 <- 'intergenic.nonsyn.haplotype.012'
# prefix_012 <- args[1]
region_haplotype_df <- read.table(prefix_012)
region_haplotype_df <- region_haplotype_df[,-1]

meta_file <- '../Sympatric_material2.txt'
# meta_file <- args[2]
meta_df <- read.table(meta_file,header=T)

indv_file <- paste0(prefix_012,'.indv')
indv_df <- read.table(indv_file)

pos_file <- paste0(prefix_012,'.pos')
pos_df <- read.table(pos_file)

rownames(region_haplotype_df) <- indv_df$V1
colnames(region_haplotype_df) <- pos_df$V2

meta_df <- meta_df[meta_df$ID %in% indv_df$V1,]
meta_df <- meta_df[order(meta_df$Lat, decreasing = T),,drop=FALSE]
index_La <- meta_df$ID
region_haplotype_df <- region_haplotype_df[index_La, ]
rownames(meta_df) <- meta_df$ID
indv_Lat <- meta_df[index_La,'Lat',drop=FALSE]

mat <- as.matrix(region_haplotype_df)

la_range=c(32,36,40)
colors <- list('-1'='#006d77','0'='#edf6f9','1'='#ffddd2','2'='#e29578')
col_fun = colorRamp2(la_range, c("#006d77", "white", "#e29578"))
herbicide_anno_value = HeatmapAnnotation(HR = anno_points(indv_Lat,
                                                          ylim = c(32, 40),
                                                          axis_param = list(
                                                            side = "top",
                                                            at = la_range
                                                            # labels = c("zero", "half", "one")
                                                          )),
                                         HR_heatmap = anno_simple(indv_Lat$Lat,
                                                                  col = col_fun),
                                         which = 'row'
)


ht <- Heatmap(
  # mat_factor,
  mat,
  name = "Genotype",
  col = colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE, 
  # show_column_dend = TRUE,
  # column_split = column_split,
  show_row_names = TRUE,
  row_names_side = "left",
  show_column_names = TRUE,
  show_column_dend = FALSE,
  # top_annotation = ha_col,
  right_annotation = herbicide_anno_value,
  column_title = "Haplotype Heatmap",
  row_names_gp = gpar(fontsize = 4),
  heatmap_legend_param = list(
    title = "Genotype",
    at = c("-1", "0", "1", "2"),
    labels = c("./.", "0/0", "0/1", "1/1")
  )
)
# ht
pdf(paste0(prefix_012,'_haplotype_vis.pdf'),width = 7,height = 4)
draw(ht, ht_gap = unit(1, "mm"))
dev.off()