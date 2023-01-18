# Rscript b1_cellbin2rds.R [cellbin matrix] [orig.ident]

args<-commandArgs(T)

if (length(args)<2) {
  stop("Usage: Rscript b1_cellbin2rds.R [cellbin matrix] [orig.ident]")
}


library(Matrix)
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(rearrr)

cellbin<-fread(args[1], header=T)
cellbin<-cellbin[which(cellbin$label>0),]
colnames(cellbin)[1:5]<-c("geneID", "x", "y", "MIDCounts", "cell")

#geneID	x	y	MIDCounts	cell
#Peg10	19875	9300	1	170934
#Rps27a	19875	9300	1	170934
#Vbp1	19875	9301	2	170934

cellbin$cell<-as.character(cellbin$cell)

message("Cell bin data looks like this:")
head(cellbin)

message("Merging counts and getting cell coordinates:")
#cb.counts<-cellbin %>% group_by(cell, geneID) %>% summarise(counts=sum(MIDCounts))
#cb.coords<-cellbin %>% group_by(cell) %>% summarise(coord_x=floor(mean(x)), coord_y=floor(mean(y)))
cb.coords<-cellbin %>% group_by(cell) %>% summarise(coord_x=round(midrange(x)), coord_y=round(midrange(y)))

#cb.counts<-as.data.frame(cb.counts)
cb.coords<-as.data.frame(cb.coords)
rownames(cb.coords)<-cb.coords$cell
cb.coords<-cb.coords[,-1]

#      coord_x coord_y
#10000    7894   11167

#mat<-with(cb.counts, tapply(counts, list(geneID, cell), FUN=sum))
##mat<-with(cellbin, tapply(MIDCounts, list(geneID, cell), FUN=sum))
#mat[is.na(mat)] = 0
##write.csv(mat, file="genexcell_matrix.csv")

gene=unique(cellbin$geneID)
cell<-unique(cellbin$cell)
gene_idx=c(1:length(gene))
cell_idx=c(1:length(cell))
names(gene_idx)=gene
names(cell_idx)=cell

mat<-sparseMatrix(i=gene_idx[cellbin$geneID],j=cell_idx[cellbin$cell],x=cellbin$MIDCounts)
rownames(mat)=gene
colnames(mat)=cell

dim(mat)

message("Creating Seurat object:")
SeuObj<-CreateSeuratObject(counts = mat)
SeuObj[["percent.mt"]] <- PercentageFeatureSet(SeuObj, pattern = "^mt-|^MT-|^Mt-")
# add coordinates to meta data
SeuObj@meta.data<-merge(SeuObj@meta.data, cb.coords, by="row.names")
rownames(SeuObj@meta.data)<-SeuObj@meta.data$Row.names
SeuObj@meta.data<-SeuObj@meta.data[,-1]

SeuObj@meta.data$orig.ident <- args[2]

SeuObj

head(SeuObj@meta.data)

sp.embedding<-SeuObj@meta.data[,c("coord_x", "coord_y")]
colnames(sp.embedding)<-c("Spatial_1", "Spatial_2") # prefix with number
sp.embedding<-as.matrix(sp.embedding) # needs to be a matrix
SeuObj$spatial<-CreateDimReducObject(embeddings=sp.embedding, key='Spatial_', assay='RNA')

saveRDS(SeuObj, paste0(args[1], ".rds"))
message("RDS file with raw counts saved.")


message("Summary of nFeature")
summary(SeuObj@meta.data$nFeature_RNA)

message("Summary of nCount")
summary(SeuObj@meta.data$nCount_RNA)


pdf(paste0(args[1], "_QC.pdf"), width=8, height=4)
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

p9<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = nFeature_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nGene") #+scale_color_viridis_c(option="turbo")

#p9 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nFeature_RNA))+
#  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+
#  labs(fill="nGene",x=NULL,y=NULL,title = NULL)

p10<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = nCount_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nUMI") #+scale_color_viridis_c(option="turbo")

#p10 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nCount_RNA))+
#    geom_tile()+  theme_void()+ coord_fixed()+ scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+
#    labs(fill="read counts",x=NULL,y=NULL,title = "spatial distribution of read counts")

grid.arrange(p9,p10,ncol = 2, top="Spatial distribution of nGene and nUMI")
grid.arrange(p6,p5,ncol = 2, top="nGene and nUMI per cell")
print(p4)
grid.arrange(p1.2,p3.2,ncol = 2)
grid.arrange(p7,p8,p2,ncol = 3)
dev.off()

message("QC plots exported.")


