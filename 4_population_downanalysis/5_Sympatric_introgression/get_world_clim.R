library(geodata)

group_site_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Sympatric_group_center2_subset.txt'
site <- read.table(group_site_file,header = TRUE,sep = '\t')
# site <- site[c(3,2,4,5,1),]
site <- site[,c(3,2)]
colnames(site)<-c('longtitude','latitude')

for (i in c("tmin", "tmax", "tavg", "prec", "wind", "vapr", "bio")) {
  MAPPP<-worldclim_global(i, res=10 ,path=i)
  bio <- extract(MAPPP, site)
  bio <- cbind(site, bio)
  write.table(bio, paste0('E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/',i,'.txt'), 
              sep = '\t', row.names = TRUE, quote = FALSE)
}

# C5_character_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Sympatric_material2.txt'
# site <- read.table(C5_character_file,header = TRUE,sep = '\t')
# # colnames(site)<-c('id','longtitude','latitude')

# for (i in c("bio", "tmin", "tmax", "tavg", "prec", "wind", "vapr")) {
#   MAPPP <- worldclim_global(i, res=10 ,path=i)
#   bio <- extract(MAPPP, site[,c(4,5)])
#   bio <- cbind(site, bio)
#   write.table(bio, paste0('Sympatric_material2_character_check_',i,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
# } 

# C5_character_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Dsan_all_material.txt'
# site <- read.table(C5_character_file,header = FALSE, sep = '\t')
# # colnames(site)<-c('id','longtitude','latitude')

# for (i in c("bio")) {
#   MAPPP <- worldclim_global(i, res=10 ,path=i)
#   bio <- extract(MAPPP, site[,c(4,5)])
#   bio <- cbind(site, bio)
#   write.table(bio, paste0('Dsan_all_material_bio_character_',i,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
# } 