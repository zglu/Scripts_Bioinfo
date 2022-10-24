library(Seurat)
library(dplyr)

Data1<-Read10X("~/zl3/Inqueries/Cheng_Sj_scRNA/_newCellRangerMapping2021-12/CellRanger_mito-4genes/F16/filtered_feature_bc_matrix/", gene.column=1) # default is col2: gene names. col1 is gene ids
Data2<-Read10X("~/zl3/Inqueries/Cheng_Sj_scRNA/_newCellRangerMapping2021-12/CellRanger_mito-4genes/F26/filtered_feature_bc_matrix/", gene.column=1) # default is col2: gene names. col1 is gene ids

Data1<-CreateSeuratObject(Data1, project = "F16", min.cells = 3, min.features = 200)
Data2<-CreateSeuratObject(Data2, project = "F26", min.cells = 3, min.features = 200)

Data1[["percent.mt"]] <- PercentageFeatureSet(object = Data1, pattern = "^Sj-")
Data2[["percent.mt"]] <- PercentageFeatureSet(object = Data2, pattern = "^Sj-")

Fil1 <- subset(x = Data1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 2.5)
Fil2 <- subset(x = Data2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 2.5)

Fil1 <- NormalizeData(Fil1, verbose = FALSE)
Fil2 <- NormalizeData(Fil2, verbose = FALSE)

Fil1 <- FindVariableFeatures(Fil1, selection.method = "vst", nfeatures = 2000, verbose = F)
Fil2 <- FindVariableFeatures(Fil2, selection.method = "vst", nfeatures = 2000, verbose = F)

# integration using RPCA

features <- SelectIntegrationFeatures(object.list = list(Fil1, Fil2))

##------- run PCA on common features
Fil1 <- ScaleData(Fil1, features = features) #single sample scale using all genes
Fil1 <- RunPCA(Fil1, features = features)

Fil2 <- ScaleData(Fil2, features = features)
Fil2 <- RunPCA(Fil2, features = features)
##-------

# You can increase the strength of alignment by increasing the k.anchor parameter, which is set to 5 by default. Increasing this parameter to 20 will assist in aligning these populations.
Fil.anchors <- FindIntegrationAnchors(object.list = list(Fil1, Fil2), anchor.features = features, reduction = "rpca") ##------ k.anchor = 20
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

saveRDS(Fil.combined, file="integrated_Log-RPCA.rds")
