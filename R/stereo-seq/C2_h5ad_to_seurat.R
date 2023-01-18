library(Seurat)
library(SeuratDisk)

args<-commandArgs(T)

Convert(args[1], dest = "h5seurat", overwrite = TRUE)

# Rscript annh5ad2rds.R --infile [h5ad file] --outfile [rds file]
