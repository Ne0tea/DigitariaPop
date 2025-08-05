library(qqman)
##USAGE: Rscrpt manhattan_draw_get_significant GR90.assoc.ps GR50.assoc.ps
args <- commandArgs(trailingOnly = TRUE)
emmax_data_file <- args[1]
# emmax_data_file <- 'F:/GWAS/emmax/Dsan_all.merge.gatk3.vcf.filtered.QC.Inbreeding.LD.emmax.GR90.assoc.ps'
emmax_data_df <- read.table(emmax_data_file, header = TRUE)
colnames(emmax_data_df)[1] <- 'CHR'
colnames(emmax_data_df)[2] <- 'SNP'
colnames(emmax_data_df)[3] <- 'BP'
colnames(emmax_data_df)[12] <- 'P'


pdf(paste0(emmax_data_file,'.pdf'), width=12, height=4)
manhattan(emmax_data_df,col = c('#006060', '#B8DADA','#FFBF3D'))
dev.off()

tiff(paste0(emmax_data_file,'.tiff'), width=2000, height=1000)
manhattan(emmax_data_df,col = c('#006060', '#B8DADA','#FFBF3D'))
dev.off()

pdf(paste0(emmax_data_file,'.qq.pdf'), width=5, height=5)
qq(emmax_data_df$P)
dev.off()

tiff(paste0(emmax_data_file,'.qq.tiff'), width=1500, height=1500)
qq(emmax_data_df$P)
dev.off()

