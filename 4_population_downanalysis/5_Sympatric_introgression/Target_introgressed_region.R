# 加载必要的库
library(ggplot2)
library(dplyr)
Chr_subg_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/0_Dsan_chr_trans'
# window_fd_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/Introgressed_50k_window_pi_mean.csv'
window_fd_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/Introgressed_10k_window_fd_mean.csv'
Chr_subg_df <- read.table(Chr_subg_file, sep = "\t",header = FALSE)
colnames(Chr_subg_df) <- c("scaffold", "chromosome", "Subgenome")

window_fd_df <- read.csv(window_fd_file)
merged_data <- merge(window_fd_df,Chr_subg_df,by = 'scaffold')
# subg_color<-c('SubC'='#003049','SubD'='#d62828','SubE'='#f77f00')

# 指定要绘制的 scaffold 和范围
selected_scaffold <- "Chr04"  # 替换为你想要的 scaffold
start_range <- 51*1e6
end_range <- 57*1e6
group_order <- c('Group_1','Group_4','Group_6','Group_12','Group_2','Group_7','Group_3','Group_8')

# 筛选数据
filtered_data <- merged_data %>%
  filter(scaffold == selected_scaffold & start >= start_range & end <= end_range)
filtered_data[is.na(filtered_data)] <- 0

filtered_data$group <- factor(filtered_data$group, levels = group_order)
# filtered_data$dxyminor <- filtered_data$dxy_allopatric_C4 - filtered_data$dxy_sympatric_C4
# 创建散点图，并根据 group 分页
ggplot(filtered_data, ) +
  geom_point(aes(x = start, y = fd), color='grey') +
  # geom_point(aes(x = start, y = pi_sympatric), color='grey') +
  # geom_point(aes(x = start, y = dxy_allopatric_C4), color='blue') +
  annotate("rect",xmin = 54200000, xmax = 54600000, ymin = -Inf, ymax = Inf,
            fill = "#780000", alpha = 0.3) +
  labs(x = paste(selected_scaffold,'genomic position'),
       y = "fd") +
  scale_x_continuous(
    breaks = seq(ceiling(start_range / 1e6) * 1e6, floor(end_range / 1e6) * 1e6 + 1, by = 2000000),  # 设置刻度位置
    labels = function(x) paste0(x / 1e6, " Mb")
  ) +
  facet_grid(group ~ .) +  # 根据 group 分面
  theme_minimal()

ggsave('E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/pi_cold_Chr04.pdf',height = 6,width = 10)
