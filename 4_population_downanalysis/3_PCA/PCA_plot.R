library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggrepel)

pca_data <- read.table("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/pca/GATK3_C5.eigenvec",header = F)
pca_data<-pca_data[,-1]
colnames(pca_data)[1]<-'name'
pca_data$name <- sub("_1$", "", pca_data$name)
eigval <- read.table("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/pca/GATK3_C5.eigenval", header = F)
pcs <- paste0("PC", 1:nrow(eigval))
eigval[nrow(eigval),1] <- 0
percentage <- eigval$V1/sum(eigval$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")
eigval_df$variance <- as.numeric(eigval_df$variance)
eigval_df$proportion <- as.numeric(eigval_df$proportion)
pc1_proportion <- paste0(round(eigval_df[1,3],2),"%")
pc2_proportion <- paste0(round(eigval_df[2,3],2),"%")

sample <- read.table("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/pca/sample_info.txt", header = F,sep='\t')
colnames(sample)<-c('name','Type','Year','Cluster')

data <- left_join(pca_data,sample,by="name")
data$Type[is.na(data$V2)] <- 'None'
colnames(data) <- c("Sample","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                    "PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20",
                    "PC21","PC22","PC23","PC24","PC25","PC26","PC27","PC28","PC29","PC30",
                    "PC31","PC32","PC33","PC34","PC35","PC36","PC37","PC38","PC39","PC40",
                    "PC41","PC42","PC43","PC44","PC45","PC46","PC47","PC48","PC49","PC50","Type","Year","Cluster")
# colnames(data) <- c("Sample","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Type","Year","Cluster")
# Important_sp<-c('CC254A','CC259A','CC262A','CC061A','CC142A','CC082A','CC027A','CC223A',
#                 'X55A','KM03A','LP23A','LP08A','LP25A','LP17A','LY08A','X25A', 'X58A', 
#                 'X21A', 'LW11A', 'KM19A', 'LY01A', 'KM20A', 'LZ06A', 'X59A', 'X22A', 
#                 'KM01A', 'KM06A', 'KM07A', 'X26A', 'LY07A', 'X20A', 'X09A', 'X24A', 
#                 'KM18A', 'LZ05A', 'X64A', 'KM11A', 'KM14A', 'LW03A', 'LW02A', 'X19A', 
#                 'X34A', 'LY10A', 'LZ04A', 'KM16A', 'X45A',
#                 #south america accession
#                 'CC019A', 'CC020A', 'CC021A', 'CC022A', 'CC023A', 'CC085A', 'CC088A', 
#                 'CC092A', 'CC098A', 'CC099A', 'CC103A', 'CC104A', 'CC116A', 'CC125A', 
#                 'CC127A', 'CC130A', 'CC131A', 'CC132A', 'CC163A', 'CC217A', 'CC218A', 
#                 'CC222A', 'CC223A', 'CC230A', 'CC269A', 'CC323A', 'CC332A', 'CC334A')
Continent_color=c('C1'='#264653','C2'='#287271','C3-1'='#2a9d8f','C3-2'='#e9c46a',
                  'C4'='#f4a261','C5'='#e76f51','Outlier'='#000000')
Type_color=c('C5-NE'='#005f73','C5-S'='#ae2012','C5-E1'='#94d2bd','C5-E2'='#e9d8a6',
             'Admix-C4'='#3d405b','Admix-E12'='#81b29a','Admix-E2S'='#f4f1de','Admix-E1S'='#e07a5f',
             'Outlier'='#000000')
year_color=c('2013'='#264653','2015'='#2a9d8f','2016'='#2a9d8f',
             '2018'='#2a9d8f','2019'='#2a9d8f','2021'='#e76f51','2023'='#9b2226')
data <- na.omit(data)
p <- ggplot(data,aes(PC2,PC3))+
  geom_point(aes(color=Type), size=2)+
  # geom_point(aes(color=Cluster),size=1)+
  # geom_text_repel(aes(label=Sample),direction="y",size=3,nudge_x = 0.02,max.overlaps =20 )+
  # stat_ellipse(aes(color=Type),level = 0.95, show.legend = FALSE, linewidth=1)+
  scale_color_manual(values =Type_color)+
  theme_classic()+
  labs(x=paste0("PC1(",pc1_proportion,")"),y=paste0("PC2(",pc2_proportion,")"))
p

p <- ggplot(data,aes(PC1,PC2))+
  # geom_point(aes(color=Type), size=3)+
  geom_point(aes(color=Cluster),size=1)+
  # geom_text_repel(aes(label=Sample),direction="y",size=3,nudge_x = 0.02,max.overlaps =20 )+
  # stat_ellipse(aes(color=Type),level = 0.95, show.legend = FALSE, linewidth=1)+
  # scale_color_manual(values =Type_color)+
  scale_color_manual(values =Continent_color)+
  # scale_x_continuous(limits = c(-0.1, 0.05))+
  # scale_y_continuous(limits = c(-0.05, 0.1))+
  theme_classic()+
  labs(x=paste0("PC1(",pc1_proportion,")"),y=paste0("PC2(",pc2_proportion,")"))
p

p <- ggplot(data,aes(PC1,PC2))+
  # geom_point(aes(color=Type), size=3)+
  geom_point(aes(color=Year),size=1)+
  # geom_text_repel(aes(label=Sample),direction="y",size=3,nudge_x = 0.02,max.overlaps =20 )+
  # stat_ellipse(aes(color=Type),level = 0.95, show.legend = FALSE, linewidth=1)+
  scale_color_manual(values =year_color)+
  # scale_color_manual(values =Continent_color)+
  # scale_x_continuous(limits = c(-0.1, 0.05))+
  # scale_y_continuous(limits = c(-0.05, 0.1))+
  theme(panel.grid = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.title = element_blank(),legend.key = element_blank(),
        axis.text = element_text(colour = "black", size=8),
        axis.title = element_text(color="black",size = 10),
        legend.text = element_text(colour = "black", size=10),
        legend.position = c(0.3,0.5),
        legend.direction = "horizontal")+
  labs(x=paste0("PC1(",pc1_proportion,")"),y=paste0("PC2(",pc2_proportion,")"))
p