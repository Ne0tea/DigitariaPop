
#Eco_color=c('C5-NE'='#005f73','C5-S'='#ae2012','C4'='#848e76','Admix'='#f1939c',
#            'C5-E1'='#94d2bd','C5-SE'='#ee9b00','C5-E2'='#e9d8a6')
Eco_color=c('C5-NE'='#005f73','C5-S'='#ae2012','C4'='#848e76',
            'C5-E1'='#94d2bd','C5-E2'='#e9d8a6')

args <- commandArgs(trailingOnly = TRUE)
cairo_pdf(args[length(args)],height = 6.25,width = 10)

first_df<-read.table(args[1],header=T)
plot(first_df$year,first_df$Ne_median, log=c("xy"), type="n", xaxt = "n", yaxt = "n", 
     xlab="Years before present (μ=6.5×10^-9)", ylab="Effective Population Size (NE)",
     xlim=c(10,1000000),ylim=c(10,10^8), axes = FALSE)
axis(1, at = 10^(1:6), cex.axis = 0.6, labels = expression(10^1, 10^2, 10^3, 10^4, 10^5, 10^6))
axis(2, at = 10^(1:8), cex.axis = 0.6, labels = expression(10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8))

lty<-c()
fill_color<-c()
paint_color<-c()
for (pop in args[seq(length(args)-1)]) {
  cur_df<-read.table(pop,header=T)
  print(pop)
  cur_color<-Eco_color[which(args == pop)]
  polygon(c(cur_df$year, rev(cur_df$year)), 
          c(cur_df$Ne_2.5., rev(cur_df$Ne_97.5.)),
          col = paste0(cur_color,'4D'), border = NA)
  lines(cur_df$year, cur_df$Ne_median, type = "s", col = cur_color, lwd = 2)
  fill_color<-c(fill_color,paste0(cur_color,'4D'))
  paint_color<-c(paint_color,cur_color)
  lty<-c(lty,1)
}

legend("topright",legend = args[seq(length(args)-1)],
      lty=rep(1,length(args)-1),
       fill=fill_color,col=paint_color,cex=0.8)
box()
grid()
dev.off()
