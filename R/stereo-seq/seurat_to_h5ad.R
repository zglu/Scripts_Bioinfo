# Seurat to h5ad on the RNA assay

args<-commandArgs(T)

library(Seurat)
library(SeuratDisk)

SeuObj<-readRDS(args[1])

DefaultAssay(SeuObj)<-'RNA'

# set data to counts; for TACCO
SeuObj@assays$RNA@data<-SeuObj@assays$RNA@counts

seu<-DietSeurat(SeuObj, counts = TRUE, data=TRUE, scale.data = FALSE, features=rownames(SeuObj), dimreducs=c("umap","spatial"), assays = "RNA")

i <- sapply(seu@meta.data, is.factor)
seu@meta.data[i] <- lapply(seu@meta.data[i], as.character)

SaveH5Seurat(seu, filename = "srt.h5seurat", overwrite = TRUE)
Convert("srt.h5seurat", paste0(args[1], "_RNA.h5ad"), assay="RNA", overwrite = TRUE)