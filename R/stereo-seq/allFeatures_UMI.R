# subset raw matrix with provided features; bin the matrix and plot nUMI
# output: spatial plot showing nUMI of cells expressing provided features
# eg. Rscript ../../allFeatures_UMI.R A202212220002.lasso.D1_VGP_T7.gem 50 transitional.ids

args<-commandArgs(T)

if (length(args)<3) {
  stop("Usage: Rscript allFeatures_UMI.R [bin1 matrix] [bin size] [feature list]")
}

library(Matrix)
library(data.table)
library(Seurat)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(gridExtra)

alldata<-fread(args[1], header=T, check.names=FALSE)
colnames(alldata)[1:4]<-c("geneID", "x", "y", "MIDCounts")

message("Original data looks like this:")
head(alldata, n=5)

#          geneID     x     y MIDCounts
#1: 0610005C13Rik 14467 17608         1
#2: 0610005C13Rik 13642 15853         1
#3: 0610005C13Rik 12873 18219         2

ids<-read.table(args[3],header=F)
ngene<-nrow(ids)
data<-alldata[which(alldata$geneID %in% ids$V1),]


binsize <- args[2]
binsize <- as.numeric(binsize)

message("Calculating bins:")
data$x = floor(data$x/binsize)
data$y = floor(data$y/binsize)
data$cellID<-paste(data$x,data$y,sep="_")

gene=unique(data$geneID) 
cell=unique(data$cellID) 
gene_idx=c(1:length(gene)) # 1:23389
cell_idx=c(1:length(cell)) # 1:5775
names(gene_idx)=gene
names(cell_idx)=cell

mat=sparseMatrix(i=gene_idx[data$geneID],j=cell_idx[data$cellID],x=data$MIDCounts)
rownames(mat)=gene
colnames(mat)=cell

message("The gene x cell count matrix:")
mat[1:5, 1:5]

dim(mat)

message("Creating Seurat Object:")
SeuObj<-CreateSeuratObject(counts = mat) # default assay: RNA
SeuObj[["percent.mt"]] <- PercentageFeatureSet(SeuObj, pattern = "^mt-|^MT-|^Mt-")
SeuObj@meta.data$coord_x=sub(rownames(SeuObj@meta.data),pattern = "_.*",replacement = "")
SeuObj@meta.data$coord_y=sub(rownames(SeuObj@meta.data),pattern = ".*_",replacement = "")
SeuObj@meta.data$coord_x=sub(SeuObj@meta.data$coord_x,pattern = "X",replacement = "")
SeuObj@meta.data$coord_x=as.integer(SeuObj@meta.data$coord_x)
SeuObj@meta.data$coord_y=as.integer(SeuObj@meta.data$coord_y)
SeuObj@meta.data$orig.ident <- paste0("bin",binsize)
 
SeuObj

head(SeuObj@meta.data)

# add coord info as embedding
sp.embedding<-SeuObj@meta.data[,c("coord_x", "coord_y")]
colnames(sp.embedding)<-c("Spatial_1", "Spatial_2") # prefix with number
sp.embedding<-as.matrix(sp.embedding) # needs to be a matrix
SeuObj$spatial<-CreateDimReducObject(embeddings=sp.embedding, key='Spatial_', assay='RNA')

#saveRDS(SeuObj, paste0(args[1], "_BIN", binsize, ".rds"))
#message("RDS file with raw counts saved.")


message("Summary of nFeature")
summary(SeuObj@meta.data$nFeature_RNA)

message("Summary of nCount")
summary(SeuObj@meta.data$nCount_RNA)


pdf(paste0(args[1], "_", args[3], "_bin", binsize, ".pdf"), width=5, height=5)
p10 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nCount_RNA))+
    geom_tile()+
    theme_void()+
    coord_fixed()+
    scale_fill_gradient(low = "#ededed", high = "red")+
    labs(fill="nUMI",x=NULL,y=NULL,title = NULL)+
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+    
    theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))

grid.arrange(p10,ncol = 1 ,top=paste0(args[3], " (n=", ngene, ")"))

dev.off()

message("UMI plot exported.")

