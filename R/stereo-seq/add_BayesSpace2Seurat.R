#Rscript add_BayesSpace2Seurat.R [Seurat rds] [BayesSpace cluster file]

args<-commandArgs(T)
if (length(args)<2) {
  stop("Usage: Rscript add_BayesSpace2Seurat.R [Seurat rds] [BayesSpace cluster file]")
}

library(Seurat)

SeuObj<-readRDS(args[1])
bayes<-read.table(args[2], header=T, sep=" ")



SeuObj@meta.data$bayesspace<-bayes$spatial.cluster
SeuObj@meta.data$bayesspace <- as.factor(SeuObj@meta.data$bayesspace)
saveRDS(SeuObj, paste0(args[1], "_BayesSpace.rds"))
