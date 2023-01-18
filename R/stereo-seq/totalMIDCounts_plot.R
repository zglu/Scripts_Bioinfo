# Rscript totalMIDCounts.R [rds] [gene ids]
args<-commandArgs(T)

library(Seurat)
SeuObj<-readRDS(args[1])
C = SeuObj@assays$RNA@counts
sum_exp<-as.matrix(rowSums(C))

idlist<-read.table(args[2], header=F)
marker_exp<-sum_exp[which(rownames(sum_exp) %in% idlist$V1),]
par(mar=c(11,4,4,4))
barplot(marker_exp, las=2, ylab="Total MIDCounts")
