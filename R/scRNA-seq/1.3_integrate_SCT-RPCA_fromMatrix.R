# Rscript 1.3_integrate_SCT-RPCA_fromMatrix.R

library(Seurat)
library(dplyr)

# each data
Data1<-Read10X("~/zl3/Inqueries/Cheng_Sj_scRNA/_newCellRangerMapping2021-12/CellRanger_mito-4genes/F16/filtered_feature_bc_matrix/", gene.column=1) # default is col2: gene names. col1 is gene ids
Data2<-Read10X("~/zl3/Inqueries/Cheng_Sj_scRNA/_newCellRangerMapping2021-12/CellRanger_mito-4genes/F26/filtered_feature_bc_matrix/", gene.column=1) # default is col2: gene names. col1 is gene ids

Data1<-CreateSeuratObject(Data1, project = "F16", min.cells = 3, min.features = 200)
Data2<-CreateSeuratObject(Data2, project = "F26", min.cells = 3, min.features = 200)

Data1[["percent.mt"]] <- PercentageFeatureSet(object = Data1, pattern = "^Sj-")
Data2[["percent.mt"]] <- PercentageFeatureSet(object = Data2, pattern = "^Sj-")

Fil1 <- subset(x = Data1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 2.5)
Fil2 <- subset(x = Data2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 2.5)

int_list<-list(Fil1, Fil2)
int_list <- lapply(X = int_list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = "percent.mt")

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
