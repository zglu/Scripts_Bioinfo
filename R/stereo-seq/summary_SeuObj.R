args<-commandArgs(T)
suppressMessages(library(Seurat))
SeuObj<-readRDS(args[1])
message("Summary of Gene counts:")
summary(SeuObj@meta.data$nFeature_RNA)

message("Summary of MID counts:")
summary(SeuObj@meta.data$nCount_RNA)
