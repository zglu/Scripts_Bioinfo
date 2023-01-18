# Rscript i1_integrate_SCT-CCA.R [data1_SCT.rds] [data2_SCT.rds]

args<-commandArgs(T)

if (length(args)<2) {
  stop("Usage: Rscript i1_integrate_SCT-CCA.R [data1_SCT.rds] [data2_SCT.rds]")
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(future)

# each data has been normalised using sctransform: Rscript 1_single-sct-clustering.R
Data1<-readRDS(args[1])
Data2<-readRDS(args[2])

int_list<-list(Data1, Data2)
features <- SelectIntegrationFeatures(object.list = int_list, nfeatures = 3000)

int_list <- PrepSCTIntegration(object.list = int_list, anchor.features = features)

## integration using CCA	



#plan("multiprocess", workers = 4)
#options(future.globals.maxSize = 500000 * 1024^2)
int.anchors <- FindIntegrationAnchors(object.list = int_list, normalization.method = "SCT",
    anchor.features = features)# dims = 1:30, k.anchor = 20)
int.combined.sct <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT")#, dims = 1:30)

int.combined.sct <- RunPCA(int.combined.sct, verbose = FALSE)
int.combined.sct <- RunUMAP(int.combined.sct, reduction = "pca", dims = 1:30)
int.combined.sct <- FindNeighbors(int.combined.sct, reduction = "pca", dims = 1:30)
int.combined.sct <- FindClusters(int.combined.sct, resolution = 0.5)

saveRDS(int.combined.sct, file="integrated_SCT-CCA.rds")

# plot dimplot
pdf("integrated_cca_dimplot.pdf", width=10, height=5)
DimPlot(int.combined.sct, reduction="umap", label = TRUE)+coord_fixed()
DimPlot(int.combined.sct, reduction="umap", split.by="orig.ident", label=TRUE)
dev.off()
