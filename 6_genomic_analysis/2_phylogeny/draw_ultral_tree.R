# 读取mcmctree软件分析得到的树文件
library(ggplot2)
library(ggtree)
library(treeio)
library(ggh4x)
library(ggtreeExtra)

phy_col=c("BAM"="#b71515","ORY"="#e97a0c","POO"="#ffde0a","PAN"="#023e7d","CHL"="#a3cef1","OUT"="#a9a29c")
file <- "E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/Phylogeny/Ultral_tree_Digitaria.tre"
phy_file <- "E:/Bio_analysis/Weed_genome_project/Digitaria/Geomic_analysis/Phylogeny/Ultral_sp_phy.txt"

phy_df<-read.table(phy_file)
colnames(phy_df)<-c('label','Type')
phy_df$yanse <- unlist(phy_col[as.character(phy_df$Type)])
# phy_df <- phy_df[!(phy_df$label %in% to_drop),]
tr <- read.mcmctree(file)
phy_tree<-full_join(tr, phy_df, by = 'label')


p <- ggtree(phy_tree)+ 
  geom_range(range="`0.95`", center="reltime", alpha=0.4, color="#003049", size=2)+
  geom_tiplab(aes(label=label),fontface="italic")+
  geom_tippoint()+
  theme_tree2()+
  # scale_x_continuous(breaks = c(seq(0,0.8,by=0.2)),
  #                    limits = c(-1.5,3),
  #                    labels = c(seq(0.8,0,by=-0.2)))+
  guides(x=guide_axis_truncated(trunc_lower = 0,
                                trunc_upper = 0.8))+
  labs(x="Age (million years)")+
  theme(axis.title.x = element_text(hjust=0.2))
p$root.edge<-0.2

colnames(phy_df)<-c('label1','Type1','yanse2')
final_plot<-p +
  geom_fruit(data=phy_df,geom=geom_tile,
             mapping=aes(x=1, y=label1, fill=Type1),
             offset=0.12,pwidth = 0.01)+
  scale_fill_manual(
    name="Clade",
    values=phy_col,
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=3
    )
  )
final_plot