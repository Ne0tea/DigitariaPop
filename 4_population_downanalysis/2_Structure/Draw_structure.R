library(pophelper)
library(gridExtra)
library(stringr)
cluster_color <- c("#001219","#005f73","#0a9396","#94d2bd","#848e76","#e9d8a6","#ee9b00","#ca6702","#f1939c","#7e1671","#bb3e03","#ae2012","#9b2226")

# fstructure_path <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_syn_555//'
# ind_label_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/Dsan_555_sample.list'
# pop_group_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/Dsan_555_pop_test.txt'

fstructure_path <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_nonsyn_499/RUN2/'
ind_label_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_499.txt'
pop_group_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_non499_pop.txt'

# fstructure_path <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_SubE/'
# ind_label_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_Sub_sample.list'
# pop_group_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_Sub_allpop.txt'


ffiles <- list.files(path=fstructure_path, pattern = "Q$", full.names=T)
ffiles <- str_sort(ffiles, numeric = TRUE)
flist <- readQ(files=ffiles)

# otherwise add your own individual labels
inds <- read.delim(ind_label_file,header=FALSE,stringsAsFactors=F)
# if all runs are equal length, add indlab to all runs
if(length(unique(sapply(flist,nrow)))==1) flist <- lapply(flist,"rownames<-",inds$V1)

# grp_ord <- c('C3-2', 'C4', 'Admix-C4', 'C5-NE','C5-E1','Admix-E1S','Admix-E12', 'C5-E2', 'Admix-E2S','C5-S','Outlier')
grp_ord <- c('Admix-C4', 'C5-NE','C5-E1','Admix-E1S','Admix-E12', 'C5-E2', 'Admix-E2S','C5-S','Outlier')
pops <- read.delim(pop_group_file, header=F,stringsAsFactors=F,fill = TRUE)
colnames(pops) <- c('Loc')
pops$Grpindx <- letters[match(pops$Loc, grp_ord)]

plotQ(flist[c(2,3,4,5,6,7,8)],imgoutput="join",panelspacer = 0.3,
          clustercol=cluster_color[c(1,2,3,4,5,6,7,9,10)],
          showlegend=T,legendkeysize=3,legendtextsize=4,
          returnplot=T,exportplot=T,height=1, width=12,
          # subsetgrp=c('C3-2', 'C4', 'Admix-C4', 'C5-NE','C5-E1','Admix-E1S','Admix-E12', 'C5-E2', 'Admix-E2S','C5-S'),
          # subsetgrp=c('C5-NE','C5-E1','Admix-E1S','Admix-E12', 'C5-E2', 'Admix-E2S','C5-S'),
          # showyaxis=T, showticks=F,
          selgrp = 'Grpindx',
          grplab=pops,ordergrp=T,grplabsize=1.5,linesize=0.5,pointsize=1,
          basesize=5,showindlab=FALSE,useindlab=T,sortind = "all",sharedindlab = FALSE,indlabsize = 2,
          splab=paste0("K=",sapply(flist[c(2,3,4,5,6,7,8)],ncol)),
          outputfilename="Dsan_all_pop",imgtype="pdf",exportpath=fstructure_path)

