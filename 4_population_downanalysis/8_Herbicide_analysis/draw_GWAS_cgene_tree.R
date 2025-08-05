library(ggtree)
library(ggplot2)

Eco_type_color=c('C5-NE'='#005f73','C5-S'='#ae2012','C5-E1'='#94d2bd','C5-E2'='#e9d8a6',
                 'Admix-C4'='#3d405b','Admix-E12'='#81b29a','Admix-E2S'='#f4f1de','Admix-E1S'='#e07a5f',
                 'C4'='#bb3e03','C3-1'='#7e1671','C3-2'='#f1939c'
                 )

# args <- commandArgs(trailingOnly = TRUE)
# tree_file <- args[1]
# class_file <- args[2]
# output_dir <- args[3]

tree_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Herbicide_analysis/Chr05.3145.JTT.contree'
class_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Herbicide_analysis/Material_Ecotype.list'

tree <- read.tree(tree_file)
class_df <- read.table(class_file)
Ecotype_dic <- setNames(class_df$V2, class_df$V1)

gene_col<-c()
for (i in tree$tip.label) {
  gene_col<-append(gene_col,Eco_type_color[Ecotype_dic[i]])
}

tregraph_data<-ggtree(tree, layout="circular",  size=0.8 )
tregraph <- tregraph_data+
  geom_tiplab(size=1, geom = 'text', align = FALSE,
              label.size=1,fontface=1) +
  geom_tippoint(size=1)+
  # geom_nodelab(geom="text")+
  theme(legend.position='none')
  # geom_rootedge(rootedge = 0.01)+
  # xlim(NA, max(tregraph_data$data$x)*1.8)
tregraph
# ggsave("annotated_tree.pdf", width = 10, height = 8)
# ggsave("annotated_tree.png", width = 10, height = 8)