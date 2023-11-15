args <- commandArgs(trailingOnly = TRUE)


library(Seurat)
library(reticulate)
library(Matrix)
library(SeuratDisk)

rds <- readRDS(args[1])
#rds <- UpdateSeuratObject(rds)

anndata <- import("anndata")
assay <- DefaultAssay(object = rds)
exprs <- GetAssayData(object = rds, slot = "data", assay = assay)

col_data <- rds[[]]
for(name in names(col_data)) {
  if(class(col_data[[name]])=='data.frame') {
    col_data[[name]] <- NULL
  }
}
col_data$ident <- Idents(object = rds)
row_data <- rds[[assay]][[]]
if (length(row_data) == 0) {
  row_data['tmp'] = 'tmp'
}
adata <- anndata$AnnData(X = t(exprs), obs = col_data, var = row_data)

# dimreds
for (dr in Seurat:::FilterObjects(object = rds, classes.keep = "DimReduc")) {
  adata$obsm$setdefault(dr, Embeddings(object = rds[[dr]]))
}

h5ad_path=paste0(args[1], ".h5ad")
adata$write(h5ad_path)

h5_seurat_path <- paste0(tools::file_path_sans_ext(h5ad_path), '.h5Seurat')
SaveH5Seurat(rds, filename = h5_seurat_path, overwrite = TRUE)
Convert(h5_seurat_path, dest = "h5ad", overwrite = TRUE)
unlink(h5_seurat_path)