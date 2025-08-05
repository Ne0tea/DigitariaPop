library(ecodist)
library(dplyr)
library(factoextra)

PSIG<- read.csv("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/PSIG_matrix.csv",header = TRUE) %>% as.matrix
diag(PSIG) <- 0
PSIG.edist <- ecodist::distance(PSIG, "euclidean")
geodistancematrix<- read.csv("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/geodis_matrix.csv",header = TRUE) %>% as.matrix
latitude<- read.csv("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/latdis_matrix.csv",header = TRUE) %>% as.matrix
longitude<- read.csv("E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/Introgressed_window/londis_matrix.csv",header = TRUE) %>% as.matrix

pre_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/prec.txt'
srad_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/vapr.txt'
wind_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/wind.txt'
bio_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/bio.txt'

tmax_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/tmax.txt'
tmin_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/tmin.txt'
tavg_file<-'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric/tavg.txt'
tmax_df<-read.table(tmax_file,header = TRUE)[,c(-1,-2,-3)]
tmin_df<-read.table(tmin_file,header = TRUE)[,c(-1,-2,-3)]
tavg_df<-read.table(tavg_file,header = TRUE)[,c(-1,-2,-3)]
tempraure_combined <- do.call(cbind, list(tmax_df, tmin_df, tavg_df))
rownames(tempraure_combined) <- NULL

all_env<-read.table(bio_file,header = TRUE)
all_env_matrix<-all_env[,c(-1,-2,-3)]
all_env_matrix.edist<-ecodist::distance(all_env_matrix, "manhattan")

for (i in c(pre_file,srad_file,wind_file)) {
  print(basename(i))
  cur_matrix<-read.table(i,header = TRUE)
  cur_matrix<-cur_matrix[,c(-1,-2,-3)]
  cur_matrix.edist<-ecodist::distance(cur_matrix, "euclidean")
  abund_temp <- mantel(cur_matrix.edist ~ PSIG.edist, mrank = FALSE, nperm = 9999)
  print(abund_temp)
  abund_temp <- mantel(cur_matrix.edist ~ PSIG.edist + all_env_matrix.edist, mrank = FALSE, nperm = 9999)
  print(abund_temp)
}

matrices <- list(latitude, longitude, geodistancematrix, tempraure_combined)
for (mat in matrices) {
  cur_matrix.edist<-ecodist::distance(mat, "euclidean")
  abund_temp <- mantel(cur_matrix.edist ~ PSIG.edist, mrank = FALSE,nperm = 9999)
  print(abund_temp)
  abund_temp <- mantel(cur_matrix.edist ~ PSIG.edist + all_env_matrix.edist, mrank = FALSE, nperm = 9999)
  print(abund_temp)
}

cur_matrix.edist<-ecodist::distance(geodistancematrix, "euclidean")
mantel(all_env_matrix.edist ~ PSIG.edist, mrank = FALSE, nperm = 9999)

# Mantel test
mantel(all91envir, PSIG)
mantel(geodistancematrix, PSIG)
mantel(longitude, PSIG)
mantel(latitude, PSIG)
mantel(pre, PSIG)
mantel(tem, PSIG)
mantel(srad, PSIG)
mantel(wind, PSIG)
mantel(all91envir, geodistancematrix)
# partial Mantel test
mantel.partial( geodistancematrix,PSIG, Fst)
mantel.partial( geodistancematrix,PSIG, all91envir)
mantel.partial(all91envir, PSIG, geodistancematrix)
mantel.partial(all91envir, PSIG, Fst)
mantel.partial(longitude, PSIG, geodistancematrix)
mantel.partial(longitude, PSIG, Fst)
mantel.partial(longitude, PSIG, all91envir)
mantel.partial(latitude, PSIG, geodistancematrix)
mantel.partial(latitude, PSIG, Fst)
mantel.partial(pre, PSIG, geodistancematrix)
mantel.partial(pre, PSIG, Fst)
mantel.partial(tem, PSIG, geodistancematrix)
mantel.partial(tem, PSIG, Fst)
mantel.partial(srad, PSIG, geodistancematrix)
mantel.partial(srad, PSIG, Fst)
mantel.partial(wind, PSIG, geodistancematrix)
mantel.partial(wind, PSIG, Fst)

### The correlation between genetic distance (FST) and environmental/geographic variables
mantel(geodistancematrix, Fst)
mantel(all91envir, Fst)
mantel(longitude, Fst)
mantel(latitude, Fst)
mantel(pre, Fst)
mantel(tem, Fst)
mantel(srad, Fst)
mantel(wind, Fst)
# partial Mantel test
mantel.partial(geodistancematrix,Fst, all91envir)
mantel.partial(all91envir, Fst, geodistancematrix)
mantel.partial(longitude, Fst, geodistancematrix)
mantel.partial(latitude, Fst, geodistancematrix)
mantel.partial(pre, Fst, geodistancematrix)
mantel.partial(tem, Fst, geodistancematrix)
mantel.partial(srad, Fst, geodistancematrix)
mantel.partial(wind, Fst, geodistancematrix)
