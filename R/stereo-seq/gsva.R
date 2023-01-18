library(GSVA)
library(msigdbi)
library(Seurat)

normdata<-readRDS("data.rds")
avgexpr<-AverageExpression(normdata, assays = "SCT", group.by = "ident")
markers <- msigdbi::read.gmt('~/spatial/analysis-test/MOSTA/markers_whole_embryo.gmt')

markers.idx<-list()
for (i in 1:length(markers[[1]])){
  tmp<-unlist(markers[[1]][i])
  markers.idx[[markers[[3]][i]]]<-tmp[which(tmp!="")]
}


markers.score<-gsva(expr=as.matrix(avgexpr$SCT), gset.idx.list=markers.idx, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=T, parallel.sz=0)

write.csv(markers.score)
