# Rscript 1_sct-clustering.R [project name eg. F16]

library(Seurat)
library(dplyr)
args<-commandArgs(T)

rawData<-Read10X("~/zl3/Inqueries/Cheng_Sj_scRNA/_newCellRangerMapping2021-12/CellRanger_mito-4genes/F26/filtered_feature_bc_matrix/", gene.column=1) # default is col2: gene names. col1 is gene ids
rawData<-CreateSeuratObject(rawData, project = args[1], min.cells = 3, min.features = 200)
rawData[["percent.mt"]] <- PercentageFeatureSet(object = rawData, pattern = "^Sj-")
Fil_sct <- subset(x = rawData, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 2.5)

Fil_sct <- SCTransform(Fil_sct, method = "glmGamPoi", vars.to.regress = "percent.mt") %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters()

saveRDS(Fil_sct, file=paste0(args[1],"_Fil_SCT_PC30.rds"))
