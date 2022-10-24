# Rscript 1_sct-clustering.R [project name eg. F16]

library(Seurat)
library(dplyr)
args<-commandArgs(T)

rawData<-Read10X("~/zl3/Inqueries/Cheng_Sj_scRNA/_newCellRangerMapping2021-12/CellRanger_mito-4genes/F26/filtered_feature_bc_matrix/", gene.column=1) # default is col2: gene names. col1 is gene ids
rawData<-CreateSeuratObject(rawData, project = args[1], min.cells = 3, min.features = 200)
rawData[["percent.mt"]] <- PercentageFeatureSet(object = rawData, pattern = "^Sj-")
Fil <- subset(x = rawData, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 2.5)

Fil <- NormalizeData(Fil, verbose=FALSE) %>%
    FindVariableFeatures(Fil, selection.method = "vst", nfeatures = 2000, verbose = F) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters()

saveRDS(Fil, file=paste0(args[1],"_Fil_PC30.rds"))
