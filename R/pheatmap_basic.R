## usage: Rscript pheatmap_basic.R [expression data with header]

library(pheatmap)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

exp<-read.delim(args[1], sep=" ", header=T)

z<-exp[,-1:-2]
rownames(z)<-exp[,1]
zz<-as.matrix(z)


#pheatmap(zz, scale = "row", cluster_cols = F, cluster_rows = F, annotation_row = anno)

outfile=paste0(args[1],".png")
# set defined cell width and height
pheatmap(zz, scale = "row", cluster_cols = F, cluster_rows = TRUE, cellheight = 10, cellwidth = 10, file = outfile,  annotation_legend=F)
# gaps_row=seq(2,538,by=2),
