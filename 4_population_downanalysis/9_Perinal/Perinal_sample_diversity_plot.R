library(dplyr)
library(reshape2)
library(entropy)
library(ggplot2)

div_scor <- function(datamelt){
  tmp<-aggregate(data=datamelt,fraction~refpop,FUN=sum)
  tmp$cumfrac<-tmp$fraction/sum(tmp$fraction)
  
  N<-length(unique(tmp$refpop))
  divscore<-(entropy.empirical(unlist(tmp$cumfrac)))/(entropy.empirical(rep(1/N,N)))
  return(divscore)
}

# args <- commandArgs(TRUE)
# proportion_file <- args[1]
proportion_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_nonsyn_499/Run2/GATK3_non_499_2.9.meanQ'
pop_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_499.txt'
perinal_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Perinal_analysis/Perinal_group_material_loc.txt'

perinal_df <- read.table(perinal_file,header=TRUE,stringsAsFactors=FALSE)
pop_data <- read.table(pop_file,stringsAsFactors = FALSE)
colnames(pop_data) <- 'ID'
data<-read.delim(proportion_file,sep="", header=FALSE,stringsAsFactors=FALSE)
data <- cbind(pop_data, data)
data$reference <- 'pop'
datamelt<-melt(data,id.vars=c("reference","ID"),variable.name="refpop",value.name="fraction")
datamelt <- datamelt[datamelt$ID %in% perinal_df$ID,]
datamelt <- merge(datamelt, perinal_df, all.x = TRUE)
datamelt$Batch[datamelt$Year %in% c(2013, 2015)] <- "Decade"
datamelt$Batch[datamelt$Year == 2023] <- "Current"
datamelt$Class <- paste(datamelt$Group, datamelt$Batch, sep = "_") 

sample_div <- by(
  data = datamelt,
  INDICES = datamelt$Class,
  FUN = function(df) div_scor(df)
)
div_result <- stack(sample_div)
div_plot_result <- div_result %>%
  mutate(
    group = sub("_(Current|Decade)$", "", ind),
    type = ifelse(grepl("Current$", ind), "Current", "Decade")
  )
div_plot_result$type<-factor(div_plot_result$type, levels = c("Decade","Current"))
ggplot(div_plot_result, aes(x = group, y = values, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ group, ncol = 3, scales ='free') +
  scale_fill_manual(values = c("Current" = "#E73847", "Decade" = "#A8DADB")) +
  theme(
    axis.line.x = element_blank(),  # 隐藏X轴线
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave('E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Perinal_analysis/Perinal_diversity_change.pdf',width = 6,height = 6)
# write.table(results,stdout(),quote=FALSE,sep="\t",row.names=FALSE)