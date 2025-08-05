# 加载必要的包
library(ggplot2)
library(dplyr)

# ltr_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/Phylogeny/Dsan_chr.V2.LTR.repeatout.out'
ltr_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/Phylogeny/DZ2_chr.LTR.repeatout.out'
# ltr_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/Phylogeny/YZGJ2_chr.LTR.repeatout.out'

chr_sub_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/0_Dsan_chr_trans'
data_df<-read.table(ltr_file,header = FALSE,stringsAsFactors = FALSE)
chr_sub_df <- read.table(chr_sub_file,header = FALSE,stringsAsFactors = FALSE)
colnames(data_df) <- c('SW','div','del','ins','chr','begin','end','class')
colnames(chr_sub_df) <- c('chr','chrsub','sub')
subC_list<-chr_sub_df[chr_sub_df$sub=='SubC','chr']
subD_list<-chr_sub_df[chr_sub_df$sub=='SubD','chr']
subE_list<-chr_sub_df[chr_sub_df$sub=='SubE','chr']

# 根据第五列的值创建一个新列
data_df$chr <- gsub("_RagTag", "", data_df$chr)
data_df <- data_df %>%
  mutate(
    subg = case_when(
      chr %in% subC_list ~ "SubC",
      chr %in% subD_list ~ "SubD",
      chr %in% subE_list ~ "SubE",
      TRUE ~ "Other"  # 处理可能的其他值
    )
  )
data_df<- data_df[data_df$subg!="Other",]
df_gypsy <- data_df %>%
  filter(grepl("Gypsy", class, ignore.case = TRUE))

df_copia <- data_df %>%
  filter(grepl("Copia", class, ignore.case = TRUE))
sub_color_set=c('SubC'='#003049','SubD'='#d62828','SubE'='#f77f00')
ggplot(df_gypsy, aes(x = div, fill = subg)) +
  geom_density(alpha = 0.5,size=0.35) +
  scale_fill_manual(values=sub_color_set)+
  labs(x = "Divergence(%)", y = "Density") +
  theme_bw()
ggplot(df_copia, aes(x = div, fill = subg)) +
  geom_density(alpha = 0.5,size=0.35) +
  scale_fill_manual(values=sub_color_set)+
  labs(x = "Divergence(%)", y = "Density") +
  theme_bw()