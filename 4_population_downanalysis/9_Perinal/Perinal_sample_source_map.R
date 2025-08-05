library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
library(maps)
library(scales)
library(sf)
library(ggspatial)
library(PieGlyph)

material_data<-read.table("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Perinal_analysis/Perinal_group_material_loc.txt",
                   fill = TRUE,header =TRUE,stringsAsFactors = F,sep = '\t')
material_center_data<-read.table("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Perinal_analysis/Perinal_group_center_loc.txt",
                    fill = TRUE,header =TRUE,stringsAsFactors = F,sep = '\t')
material_data$Lon<-as.numeric(material_data$Lon)
material_data$Lat<-as.numeric(material_data$Lat)
material_data$Year<-as.character(material_data$Year)

material_center_data$Lat<-as.numeric(material_center_data$Lat)
material_center_data$Lon<-as.numeric(material_center_data$Lon)
material_center_data$Count<-as.numeric(material_center_data$Count)
material_center_data$Current<-as.numeric(material_center_data$Current)
material_center_data$Decade<-as.numeric(material_center_data$Decade)

Type_color<-c('Current'='#E73847','Decade'='#A8DADB')

year_color=c('2013'='#1D3557','2015'='#A8DADB','2016'='#2a9d8f',
             '2018'='#2a9d8f','2019'='#2a9d8f','2021'='#e76f51','2023'='#E73847')

china_pro <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json")

Chinses_map<-ggplot(data = china_pro) + 
  geom_sf(color='black',fill='white',linewidth=0.02) +
  coord_sf(xlim=c(108,122),ylim=c(32,40), expand = FALSE)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5))

full_sample_plot<-Chinses_map+
  geom_point(data = material_data, aes(x = Lon, y = Lat,group=NA, color=Year),
             size=2,alpha = 0.7, shape = 16)+
  scale_colour_manual(values=year_color)+
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tl",#调整指北针位置，
                         style = north_arrow_fancy_orienteering(
                           fill = c("grey40", "white"),
                           line_col = "grey20"))+
  theme(
    legend.position = c(0.8, 0.2),
    legend.box = 'horizonral')
full_sample_plot
# ggsave('Perinal_material_loc.pdf',path='E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Perinal_analysis')

Group_plot<-Chinses_map + 
  geom_pie_glyph(data=material_center_data,aes(x = Lon, y = Lat,radius = Count/6),
                 color='white',alpha=0.8,
                 slices = c('Current','Decade'),
                 inherit.aes = FALSE) +
  geom_text_repel(data=material_center_data,aes(x = Lon + 0.1, y = Lat - 0.3, label = Group)) +
  scale_radius_continuous(range = c(0.25, 0.6),
                          labels = function(x) {6*x}, 
                          name = 'Size',
                          unit = 'cm') +
  scale_fill_manual(values = Type_color,
                    # labels = paste0('G', 1:3), 
                    name = 'Record') +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tl",#调整指北针位置，
                         style = north_arrow_fancy_orienteering(
                           fill = c("grey40", "white"),
                           line_col = "grey20"))+
  guides(radius = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  theme(
    legend.position = c(0.8, 0.2),
    legend.box.just = "right", 
    legend.box = 'horizontal')
Group_plot
# ggsave('Perinal_group_center_loc.pdf',path='E:/Bio_analysis/Weed_genome_project/Digitaria/Population/Perinal_analysis')