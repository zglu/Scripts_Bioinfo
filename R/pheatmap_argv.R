## usage: Rscript pheatmap_basic.R [ids] [expression data with header]

library(pheatmap)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

ids<-read.table(args[1], header=F)
# with id column
exp<-read.delim(args[2], sep=" ", header=T)

ids$match<-match(ids$V1, exp$id, nomatch=0)
exp<-exp[ids$match,]
# if add group names in the second column
#z<-exp[,-1:-2]
z<-exp[,-1]
rownames(z)<-exp[,1]
zz<-as.matrix(z)


#pheatmap(zz, scale = "row", cluster_cols = F, cluster_rows = F, annotation_row = anno)

outfile=paste0(args[1],".heatmap.png")
# set defined cell width and height
pheatmap(zz, scale = "row", cluster_cols = F, cluster_rows = TRUE, cellheight = 10, cellwidth = 10, main=args[1], file = outfile,  annotation_legend=F)
# gaps_row=seq(2,538,by=2),
