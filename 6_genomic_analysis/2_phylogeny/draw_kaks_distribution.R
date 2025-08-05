# 安装并加载必要的包
rm(list=ls())
library(ggridges)
library(ggplot2)
library(dplyr)

kaks_dir<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/Phylogeny/KaKs_result'

file_list<-dir(kaks_dir)
data_df <- read.table(paste0(kaks_dir,'/',file_list[1]),header = TRUE,stringsAsFactors = FALSE)
data_df$Type <- gsub('_all-results.txt','',file_list[1])

for (i in 2:length(file_list)) {
  cur_df <- read.table(paste0(kaks_dir,'/',file_list[i]),header = TRUE,stringsAsFactors = FALSE)
  cur_df$Type <- gsub('_all-results.txt','',file_list[i])
  data_df <- rbind(data_df, cur_df)
  }

data_df$Type <- factor(data_df$Type,levels = c('Dsan_subC_Drad','Dsan_subD_DmilD','Dsan_subE_DmilE',
                                               'DmilD_DmilE','Dsan_subD_subE','Dsan_subE_DmilD',
                                               'Dsan_subC_subD','Dsan_subC_DmilD','Dsan_subD_Drad',
                                               'Dsan_subC_subE','Dsan_subC_DmilE','Dsan_subE_Drad',
                                               'Dsan_subC_DexilA','Dsan_subD_DexilA','Dsan_subE_DexilA',
                                               'DexiA_Svir','DexiB_Svir', 
                                               'Dsan_subC_Ecur','Dsan_subD_Ecur','Dsan_subE_Ecur',
                                               'Dsan_subC_Osat','Dsan_subD_Osat','Dsan_subE_Osat')) 
# 绘制山脊图
density_data <- na.omit(data_df) %>%
       group_by(Type) %>%
       do({
             density_vals <- density(.$Ks)
             max_y <- max(density_vals$y)
             maxY_x <- density_vals$x[which.max(density_vals$y)]
             data.frame(DvT = maxY_x/2/6.5*1000000000, max_y = max_y, maxY_x = maxY_x)
         })
data_df <- na.omit(data_df)
ggplot(data_df, aes(x = Ks, y = Type, fill = Type)) +
  geom_density_ridges(scale = 1, size = 0.5, alpha = 0.7) +
  theme_ridges() +
  theme(axis.title.y = element_blank())+
  coord_cartesian(xlim = c(0, 1))+
  labs(x = "Ks (synonymous substitution rate)")+
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw()

