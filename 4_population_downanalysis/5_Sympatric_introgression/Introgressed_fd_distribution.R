# 加载必要的库
library(ggplot2)
Chr_subg_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/0_Dsan_chr_trans'
window_fd_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/Introgressed_100k_window_fd_mean.csv'
Chr_subg_df <- read.table(Chr_subg_file, sep = "\t",header = FALSE)
colnames(Chr_subg_df) <- c("scaffold", "chromosome", "Subgenome")

window_fd_df <- read.csv(window_fd_file)
merged_data <- merge(window_fd_df,Chr_subg_df,by = 'scaffold')
subg_color<-c('SubC'='#003049','SubD'='#d62828','SubE'='#f77f00')

for (group in unique(merged_data$group)) {
  cur_df <- merged_data[merged_data$group==group,]
  cur_plot<-ggplot(cur_df, aes(x = start, y = fd)) +
    # geom_line(aes(color = Subgenome), alpha = 0.5) + 
    geom_point(aes(color = Subgenome), size = 0.25) + 
    scale_color_manual(values = subg_color) +
    labs(x = "Subgenome",
         y = "fd") +
    facet_grid(chromosome ~ Subgenome) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(dirname(window_fd_file),'/',group,'_fd_distribution.pdf'),plot=cur_plot,height = 6, width = 12)
}

