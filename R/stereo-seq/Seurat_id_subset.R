# subsetting using selection ids

args<-commandArgs(T)

if (length(args)<2) {
  cat("Usage: Rscript Seurat_id_subset.R [rds] [id.txt]\n")
  q()
}

seuFile<-args[1]
selIds<-args[2]

SeuObj<-readRDS(seuFile)

SeuObj$cellNames<-rownames(SeuObj@meta.data)
toKeep.ids<-read.table(selIds, header=F) # 169 cells
#toKeep.ids<-newsc$cellNames[which(!(newsc$cellNames %in% toExcl.ids$V1))]
SeuObj.sel<-subset(SeuObj, subset=cellNames %in% toKeep.ids$V1)

saveRDS(newsc.rest, file=paste0(seuFile, "_sel.rds"))