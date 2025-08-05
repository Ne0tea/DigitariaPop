library(ggplot2)
library(PieGlyph)
library(tidyr)
library(tidyterra)
library(dplyr)
library(maps)
library(ggsci)
library(scales)
library(geodata)
library(ggrepel)
library(raster)
library(ggspatial)
library(ggnewscale)
library(sf)

hap_colors <- c('-1'='#006d77','0'='#edf6f9','1'='#ffddd2','2'='#c1121f')
hap_type_colors <- c('H1'='#edf6f9', 'H4'='#e63946', 'H3'='#efcac4')

###Import file
setwd("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Temp")
# target_pos <- '54406537'
target_pos_set <- c('54402930','54406026','54407007','54407074','54407753','54407771')
# prefix_012 <- 'Chr04_54Mb.nonsyn.haplotype.012'
prefix_012 <- 'Chr04_54Mb_all_C45.nosym.haplotype.012'
# prefix_012 <- args[1]
region_haplotype_df <- read.table(prefix_012)
region_haplotype_df <- region_haplotype_df[,-1]

# meta_file <- '../Sympatric_material2.txt'
meta_file <- '../../Dsan_all_material.txt'
# meta_file <- args[2]
meta_df <- read.table(meta_file,header=T)

Sympatric_group_df <- read.table("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Sympatric_group_center2.txt",
                                 fill = TRUE,header =TRUE,stringsAsFactors = F,sep = '\t',check.names = FALSE)

indv_file <- paste0(prefix_012,'.indv')
indv_df <- read.table(indv_file)

pos_file <- paste0(prefix_012,'.pos')
pos_df <- read.table(pos_file)

rownames(region_haplotype_df) <- indv_df$V1
colnames(region_haplotype_df) <- pos_df$V2
region_haplotype_df <- region_haplotype_df[target_pos_set]

region_haplotype_df[region_haplotype_df == -1] <- 0
region_haplotype_df[region_haplotype_df == 1] <- 2
haplotypes <- apply(region_haplotype_df, 1, paste, collapse = "|")
meterial_haptype <- data.frame(
  Material = rownames(region_haplotype_df),
  Haplotype = haplotypes,
  HaplotypeID = paste0("H", match(haplotypes, unique(haplotypes))),
  stringsAsFactors = FALSE
)
haplotype_counts <- meterial_haptype %>% 
  count(HaplotypeID, Haplotype, name = "Count")
###HaplotypeID get HapID that you want
meterial_haptype <- meterial_haptype[meterial_haptype$HaplotypeID %in% names(hap_type_colors),]

meta_df <- meta_df[meta_df$ID %in% indv_df$V1,]
# meta_df$hap <- target_var_df[match(meta_df$ID, rownames(target_var_df)), target_pos]
# meta_df$hap <- as.character(meta_df$hap)
meta_df$haptype <- meterial_haptype[match(meta_df$ID, rownames(meterial_haptype)), 'HaplotypeID']
meta_df <- na.omit(meta_df)
write.table(meta_df,file = 'Dsan_material_hap.txt', quote=FALSE, row.names = FALSE)
###if want haplotype distribution in phylogeny
if (FALSE) {
  meta_hap_count_df <- meta_df %>%
    group_by(Ecotype, haptype) %>%
    summarise(
      Count = n()
    ) %>%
    filter(Ecotype %in% c('C4', 'C5-NE', 'C5-E1', 'C5-E2', 'C5-S'))
  meta_hap_count_df$Ecotype <- factor(meta_hap_count_df$Ecotype, levels = c('C4', 'C5-NE', 'C5-E1', 'C5-E2', 'C5-S'))
  haptype_plot <- ggplot(meta_hap_count_df, aes(x = Ecotype, y = Count, fill = haptype)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = hap_type_colors) + 
    scale_y_continuous(labels = scales::percent) +  # y轴显示为百分比
    labs(y = "Percentage", fill = "Haptype") +
    theme_minimal()
  ggsave('Sympatric_Chr04_OsRZ3_Haptype_All_sample_hapdistribution.pdf',plot=haptype_plot, width=5.15,height=4,
         path='E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric')
}


# china_pro <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json")
# worldclim <- worldclim_global(name="worldclim", var="bio", res=2.5, path='Bio19')

region_extent <- extent(105, 124, 28, 42)
clim_data <- crop(worldclim[[6]], region_extent)
# bbox <- st_bbox(china_pro)
# clim_data <- crop(worldclim[[6]], bbox)

meta_df <- meta_df[meta_df$Cluster!='C4',]
full_sample_plot<-ggplot() +
  geom_spatraster(data = clim_data, maxcell=5e+10) +
  # scale_fill_distiller(palette = "RdBu", direction = -1,limits=c(-27,15),na.value = "white" ) +
  scale_fill_distiller(palette = "RdBu", direction = -1,limits=c(-15,-3),na.value = "white" ) +
  new_scale("fill") +
  geom_sf(data = china_pro, fill = NA, color = "grey", linewidth = 0.3) +  # 省界叠加
  # coord_sf(xlim=c(102,130),ylim=c(18,47), expand = TRUE)+
  coord_sf(xlim=c(108,122),ylim=c(32,40), expand = TRUE)+
  geom_point(data = meta_df, aes(x = Lon, y = Lat, fill=haptype),
             shape=21,size=2,alpha = 0.8, stroke = 0.5)+
  # geom_point(data = sample[sample$Sp=='Ds',], aes(x = Lon, y = Lat, fill=Ecotype),shape=21,size=2,alpha = 0.7)+
  scale_fill_manual(values = hap_type_colors,
                    # labels = paste0('G', 1:3), 
                    name = 'Haplotype')+
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tl",#调整指北针位置，
                         style = north_arrow_fancy_orienteering(
                           fill = c("grey40", "white"),
                           line_col = "grey20")) +
  guides(fill = guide_legend(order = 2,nrow = 1)) +
  labs(x = "", y = "")+
  theme_bw()+
  theme(
    legend.position = c(0.1, -0.1),
    # legend.box.just = "left", 
    legend.box = 'vertical')
full_sample_plot
# ggsave('Sympatric_Chr04_OsRZ3_Haptype_sample_loc.pdf',plot=full_sample_plot, width=6,height=4,
#        path='E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric')
ggsave('Sympatric_Chr04_OsRZ3_Haptype_All_sample_zo_loc.pdf',plot=full_sample_plot, width=6,height=4,
       path='E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric')

meta_group_df <- meta_df %>%
  pivot_wider(
    id_cols = "Cluster",
    names_from = "haptype",
    values_from = "haptype",
    values_fn = length,
    values_fill = 0
  )

meta_group_df <- merge(meta_group_df, Sympatric_group_df)
# meta_group_df$validnumber <- meta_group_df[names(hap_type_colors)]
meta_group_df$validnumber <- apply(meta_group_df[names(hap_type_colors)], 1,sum)
# if (all(c("0", "2") %in% names(meta_group_df))) {
  # if (all(c("0","1", "2") %in% names(meta_group_df))) {
  #   var_set <- c("0","1", "2")
  # } else{
  #   var_set <- c("0","2")
  # }
  full_sample_plot <- ggplot() +
    geom_spatraster(data = clim_data, maxcell=5e+10) +
    scale_fill_distiller(palette = "RdBu", direction = -1,limits=c(-15,-3),na.value = "white" ) +
    new_scale("fill") +
    geom_sf(data = china_pro, fill = NA, color = "grey", linewidth = 0.3) +  # 省界叠加
    coord_sf(xlim=c(108,122),ylim=c(32,40), expand = TRUE)+
    geom_pie_glyph(data=meta_group_df,aes(x=Center_Lon, y = Center_Lat,radius = validnumber),
                   color='black',linewidth=0.25,
                   slices = names(hap_type_colors),
                   inherit.aes = FALSE)+
    geom_text_repel(data=meta_group_df,aes(x=Center_Lon, y = Center_Lat+1, label = Cluster), 
                    size = 3)+
    scale_radius_continuous(range = c(0.25, 0.75),
                            labels = function(x) {x}, 
                            name = 'Count',
                            unit = 'cm') +
    scale_fill_manual(values = hap_type_colors,
                      # labels = paste0('G', 1:3), 
                      name = 'Hap')+
    annotation_scale(location = "br") +
    annotation_north_arrow(location = "tl",#调整指北针位置，
                           style = north_arrow_fancy_orienteering(
                             fill = c("grey40", "white"),
                             line_col = "grey20"))+
    guides(radius = guide_legend(order = 1,nrow = 1), 
           fill = guide_legend(order = 2,nrow = 1)) +
    labs(x = "", y = "")+
    theme_bw()+
    theme(
      legend.position = c(0.5, -0.1),
      # legend.box.just = "left", 
      legend.box = 'vertical')
  # full_sample_plot

  ggsave('Sympatric_Chr04_OsRZ3_Haptype_cluster_loc.pdf',plot=full_sample_plot, width=6,height=4,
         path='E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric')
# }
