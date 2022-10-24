# Rscript 1_integrate_sct.R [data1_SCT.rds] [data2_SCT.rds]

library(Seurat)
library(dplyr)
args<-commandArgs(T)

# each data has been normalised using sctransform: Rscript 1_single-sct-clustering.R
Data1<-readRDS(args[1])
Data2<-readRDS(args[2])

int_list<-list(Data1, Data2)
features <- SelectIntegrationFeatures(object.list = int_list, nfeatures = 3000)

int_list <- PrepSCTIntegration(object.list = int_list, anchor.features = features)

## integration using RPCA

int_list <- lapply(X = int_list, FUN = RunPCA, features = features)

int.anchors <- FindIntegrationAnchors(object.list = int_list, normalization.method = "SCT",
    anchor.features = features, reduction = "rpca")# , dims = 1:30, k.anchor = 20)
int.combined.sct <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT")#, dims = 1:30)

int.combined.sct <- RunPCA(int.combined.sct, verbose = FALSE)
int.combined.sct <- RunUMAP(int.combined.sct, reduction = "pca", dims = 1:30)
int.combined.sct <- FindNeighbors(int.combined.sct, reduction = "pca", dims = 1:30)
int.combined.sct <- FindClusters(int.combined.sct, resolution = 0.5)

saveRDS(int.combined.sct, file="integrated_SCT-RPCA.rds")
