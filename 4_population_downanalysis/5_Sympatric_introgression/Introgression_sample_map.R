library(ggplot2)
library(scatterpie)
library(tidyr)
library(dplyr)
library(maps)
library(ggsci)
library(scales)
library(ggspatial)
library(sf)
material_data <- read.table("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Sympatric_material2.txt",
                            fill = TRUE,header =TRUE,stringsAsFactors = F,sep = '\t')
Sym_C4_data <- material_data[material_data$Ecotype == "C4", ]
Sym_C5_data <- material_data[material_data$Ecotype != "C4", ]

Sym_C4_data$Sp<-'Dc'
Sym_C5_data$Sp<-'Ds'

sample<-rbind(Sym_C4_data,Sym_C5_data)
sample$Lon<-as.numeric(sample$Lon)
sample$Lat<-as.numeric(sample$Lat)
sample$Alt<-as.numeric(sample$Alt)
sample$Ecotype<-as.character(sample$Ecotype)

Eco_color=c('C5-NE'='#005f73','C5-S'='#ae2012','C5-E1'='#94d2bd','C5-E2'='#e9d8a6',
             'Admix-C4'='#3d405b','Admix-E12'='#81b29a','Admix-E2S'='#f4f1de','Admix-E1S'='#e07a5f',
             'OUT'='#000000',
             'C4'='#000000','C1'='#000000','C2'='#000000','C3-1'='#000000','C3-2'='#000000')

Sp_shape=c('Ds'=21,'Dc'=24)

china_pro <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json")

sample$Ecotype<-factor(sample$Ecotype,levels = unique(sample$Ecotype))
region_extent <- extent(105, 124, 28, 42)
clim_data <- crop(worldclim[[6]], region_extent)

full_sample_plot<-ggplot() +
  geom_spatraster(data = clim_data, maxcell=5e+10) +
  scale_fill_distiller(palette = "RdBu", direction = -1,limits=c(-15,-3),na.value = "white" ) +
  new_scale("fill") +
  geom_sf(data = china_pro, fill = NA, color = "grey", linewidth = 0.3) +  # 省界叠加
  coord_sf(xlim=c(108,122),ylim=c(32,40), expand = TRUE)+
  geom_point(data = sample[sample$Sp=='Ds',], aes(x = Lon, y = Lat, fill=Ecotype),shape=21,size=2,alpha = 0.7)+
  geom_point(data = sample[sample$Sp=='Dc',], aes(x = Lon, y = Lat, fill=Ecotype),shape=24,size=2,alpha = 0.7)+
  scale_fill_manual(values = Eco_color,
                    # labels = paste0('G', 1:3), 
                    name = 'Type')+
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tl",#调整指北针位置，
                         style = north_arrow_fancy_orienteering(
                           fill = c("grey40", "white"),
                           line_col = "grey20")) +
  guides(radius = guide_legend(order = 1,nrow = 1), 
         fill = guide_legend(order = 2,nrow = 2),
         value = guide_legend(order = 3,nrow = 1, direction='horizontal')) +
  labs(x = "", y = "")+
  theme_bw()+
  theme(
    legend.position = c(0.5, -0.1),
    # legend.box.just = "left", 
    legend.box = 'vertical')
# full_sample_plot
ggsave('Sympatric_MT_sample_loc2.pdf',plot=full_sample_plot, width=6,height=4,
       path='E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric')
