library(Seurat)
library(dplyr)

args<-commandArgs(T)

# each data was done ScaleData before: Rscript 2_single-lognorm-clustering.R
Fil1<-readRDS(args[1])
Fil2<-readRDS(args[2])

# integration using RPCA
int_list<-list(Fil1, Fil2)
features <- SelectIntegrationFeatures(object.list = list(Fil1, Fil2))








# You can increase the strength of alignment by increasing the k.anchor parameter, which is set to 5 by default. Increasing this parameter to 20 will assist in aligning these populations.
Fil.anchors <- FindIntegrationAnchors(object.list = list(Fil1, Fil2), anchor.features = features) ##------ k.anchor = 20
# this command creates an 'integrated' data assay
Fil.combined <- IntegrateData(anchorset = Fil.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Fil.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
Fil.combined <- ScaleData(Fil.combined, verbose = FALSE)
Fil.combined <- RunPCA(Fil.combined, npcs = 30, verbose = FALSE)
Fil.combined <- RunUMAP(Fil.combined, reduction = "pca", dims = 1:30)
Fil.combined <- FindNeighbors(Fil.combined, reduction = "pca", dims = 1:30)
Fil.combined <- FindClusters(Fil.combined, resolution = 1)

saveRDS(Fil.combined, file="integrated_Log-CCA.rds")
