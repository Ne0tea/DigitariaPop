library(pophelper)
library(gridExtra)
cluster_color <- c("#001219","#005f73","#0a9396","#94d2bd","#848e76","#e9d8a6","#ee9b00","#ca6702","#bb3e03","#ae2012","#9b2226")

fstructure_path <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/GATK3_555'
ind_label_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/Dsan_555_sample.list'
pop_group_file <- 'E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/P_structure/Dsan_555_pop.list'

ffiles <- list.files(path=fstructure_path, pattern = "Q$", full.names=T)
flist <- readQ(files=ffiles)
# otherwise add your own individual labels
inds <- read.delim(ind_label_file,header=FALSE,stringsAsFactors=F)
# if all runs are equal length, add indlab to all runs
if(length(unique(sapply(flist,nrow)))==1) flist <- lapply(flist,"rownames<-",inds$V1)

sample<-read.table('E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression/Sympatric_material.txt',
                   fill = TRUE,header =TRUE,stringsAsFactors = F,sep = '\t',check.names = FALSE)
sample <- sample %>%
  mutate(Group = ifelse(Ecotype == "C4", "C4", Cluster))

filtered_list <- lapply(flist, function(df) {
  matched_rows <- rownames(df) %in% paste0(sample$ID, "_1")
  df[matched_rows, , drop = FALSE]
})
target_sample<-sub("_1$", "",rownames(filtered_list[1]$GATK3_555_1.2))
result <- sample[sample$ID %in% target_sample, c("ID", "Group")]
pops <- data.frame(result[match(target_sample, result$ID), "Group"],stringsAsFactors = FALSE)
colnames(pops)<-"Group"
pops$Group<-as.character(pops$Group)
plotQ(filtered_list[c(1,2,3,4,5)],imgoutput="join",panelspacer = 0.3,
      clustercol=cluster_color[c(1,2,4,5,6,7,9,11)],
      showlegend=T,legendkeysize=3,legendtextsize=4,
      returnplot=T,exportplot=T,height=1, width=12,
      # showyaxis=T, showticks=F,
      grplab=pops,ordergrp=T,grplabsize=1.5,linesize=0.5,pointsize=1,
      basesize=5,showindlab=FALSE,useindlab=T,sortind = "all",sharedindlab = FALSE,indlabsize = 2,
      splab=paste0("K=",sapply(filtered_list[c(1,2,3,4,5)],ncol)),
      outputfilename="GATK3_555",imgtype="pdf",exportpath='E:/Bio_analysis/Weed_genome_project/Digitaria/Population/03pop_structure/Introgression')