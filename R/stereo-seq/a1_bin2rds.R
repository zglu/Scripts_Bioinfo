# bin the matrix and save to RDS file
# Rscript a1_bin2rds.R [bin1 matrix gem/tsv] [bin size] [orig.ident]
# output: pdf with QC plots + binned RDS file

args<-commandArgs(T)

if (length(args)<3) {
  stop("Usage: Rscript a1_bin2rds.R [bin1 matrix] [bin size] [orig.ident]")
}

library(Matrix)
library(data.table)
library(Seurat)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(gridExtra)

data<-fread(args[1], header=T, check.names=FALSE)
colnames(data)[1:4]<-c("geneID", "x", "y", "MIDCounts")

message("Original data looks like this:")
head(data, n=5)

#          geneID     x     y MIDCounts
#1: 0610005C13Rik 14467 17608         1
#2: 0610005C13Rik 13642 15853         1
#3: 0610005C13Rik 12873 18219         2

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
SeuObj@meta.data$orig.ident <- args[3]
 
SeuObj

head(SeuObj@meta.data)

# add coord info as embedding
sp.embedding<-SeuObj@meta.data[,c("coord_x", "coord_y")]
colnames(sp.embedding)<-c("Spatial_1", "Spatial_2") # prefix with number
sp.embedding<-as.matrix(sp.embedding) # needs to be a matrix
SeuObj$spatial<-CreateDimReducObject(embeddings=sp.embedding, key='Spatial_', assay='RNA')

saveRDS(SeuObj, paste0(args[1], "_BIN", binsize, ".rds"))
message("RDS file with raw counts saved.")


message("Summary of nFeature")
summary(SeuObj@meta.data$nFeature_RNA)

message("Summary of nCount")
summary(SeuObj@meta.data$nCount_RNA)


pdf(paste0(args[1], "_BIN", binsize, ".rds_QC.pdf"), width=10, height=5)
#p1<-FeatureScatter(SeuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1.2<-ggplot(SeuObj@meta.data,aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
#  	stat_smooth(method=lm) +
#  	scale_x_log10() + 
#  	scale_y_log10() + 
#  	geom_vline(xintercept = 500) +
#  	geom_hline(yintercept = 250) +
  	theme_classic() 

# complexity
SeuObj$log10GenesPerUMI <- log10(SeuObj$nFeature_RNA) / log10(SeuObj$nCount_RNA)
p2<-ggplot(SeuObj@meta.data, aes(x=log10GenesPerUMI)) +
  	geom_density(alpha = 0.4, fill="khaki1") + labs(title = "log10 Genes per UMI")+
  	theme_classic() #+
  	#geom_vline(xintercept = 0.8)

#p3<-FeatureScatter(SeuObj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
p3.2<-ggplot(SeuObj@meta.data,aes(x=nFeature_RNA, y=percent.mt, color=nCount_RNA)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
#  	scale_x_log10() + 
#  	scale_y_log10() + 
#  	geom_vline(xintercept = 250) +
#  	geom_hline(yintercept = 2.5) +
  	theme_classic()

p4<-VlnPlot(SeuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)

median <- median(SeuObj@meta.data$nCount_RNA)
p5 <- ggplot(SeuObj@meta.data,aes(orig.ident,nCount_RNA))+
    geom_violin() + 
    labs(y="nUMI",x=NULL) +
    theme(panel.grid.major=element_line(colour=NA))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 15))+   
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),plot.title = element_text(hjust = 1,vjust=-40,size=14) )+
    theme(axis.text.x = element_text(size = 15,color="black"),axis.text.y = element_text(size = 15,color="black")) +
    geom_boxplot(width=0.2)+ annotate("text",label=as.character(median),y=as.numeric(median),x=1.1,hjust = 0,color='red',size=5)

median <- median(SeuObj@meta.data$nFeature_RNA)
p6 <- ggplot(SeuObj@meta.data,aes(orig.ident,nFeature_RNA))+
    geom_violin() + 
    labs(y="nGene",x=NULL) +
    theme(panel.grid.major=element_line(colour=NA))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 15))+   
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),plot.title = element_text(hjust = 1,vjust=-40,size=14) )+
    theme(axis.text.x = element_text(size = 15,color="black"),axis.text.y = element_text(size = 15,color="black")) +
    geom_boxplot(width=0.2)+ annotate("text",label=as.character(median),y=as.numeric(median),x=1.1,hjust = 0,color='red',size=5)


p7 <- ggplot(SeuObj@meta.data, aes(x=nFeature_RNA)) + 
      geom_density(alpha = 0.4, fill="lightblue")+labs(title = "distribution of nGene")+ 
      theme_classic()#+
      #theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+    
      #theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))

p8 <- ggplot(SeuObj@meta.data, aes(x=nCount_RNA)) + 
      geom_density(alpha = 0.4, fill="bisque")+labs(title = "distribution of nUMI")+
      theme_classic()#+
      #theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+    
      #theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))


p9 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nFeature_RNA))+
  geom_tile()+
  theme_void()+
  coord_fixed()+
  scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+
  labs(fill="nGene",x=NULL,y=NULL,title = NULL)+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+ 
  theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))

p10 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nCount_RNA))+
    geom_tile()+
    theme_void()+
    coord_fixed()+
    scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+
    labs(fill="nUMI",x=NULL,y=NULL,title = NULL)+
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+    
    theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))

grid.arrange(p9,p10,ncol = 2 ,top=paste0("Spatial distribution of nGene and nUMI (bin",binsize,")"))
grid.arrange(p6,p5,ncol = 3, top="nGene and nUMI per cell (bin unit)")
grid.arrange(p7,p8,p2,ncol = 3)
#print(p4)
grid.arrange(p1.2,p3.2,ncol = 2, top="Scatter plot of each cell (bin unit)")

C = SeuObj@assays$RNA@counts
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[10:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.05, las = 1, main="Top expressed genes", xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[10:1], horizontal = TRUE)

dev.off()

message("QC plots exported.")

