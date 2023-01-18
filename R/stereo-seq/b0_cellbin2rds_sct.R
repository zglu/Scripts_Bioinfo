# Rscript b0_cellbin2rds_sct.R [cellbin matrix] [orig.ident]

args<-commandArgs(T)

if (length(args)<2) {
  stop("Usage: Rscript b1_cellbin2rds.R [cellbin matrix] [orig.ident]")
}


library(Matrix)
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(rearrr)

cellbin<-fread(args[1], header=T)

#          geneID     x     y MIDCount  label tag
#1: 0610005C13Rik 12592 12373        1      0   0
#2: 0610005C13Rik  5475  6814        1  18812   0

message("Cell bin data looks like this:")
head(cellbin)

cellbin<-cellbin[which(cellbin$label>0),]
colnames(cellbin)[1:5]<-c("geneID", "x", "y", "MIDCounts", "cell")

cellbin$cell<-as.character(cellbin$cell)

# get DNB count for each cell
cellbin$DNB<-paste0(cellbin$x, "_", cellbin$y)
DNB_count<-aggregate(DNB~cell,cellbin,function(x) length(unique(x)))
rownames(DNB_count)<-DNB_count$cell
DNB_count$cell<-NULL

message("Merging counts and getting cell coordinates:")
# set center as cell coordinates: (min+max)/2
#cb.coords<-cellbin %>% group_by(cell) %>% summarise(coord_x=floor(mean(x)), coord_y=floor(mean(y)))
cb.coords<-cellbin %>% group_by(cell) %>% summarise(coord_x=round(midrange(x)), coord_y=round(midrange(y)))
cb.coords<-as.data.frame(cb.coords)
rownames(cb.coords)<-cb.coords$cell
cb.coords<-cb.coords[,-1]

cb.coords<-cbind(cb.coords, DNB_count)

#   coord_x coord_y  DNB
#1     7928   10106 1337
#10    7844   10366 2088

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

# add spatial embedding
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


#p1<-FeatureScatter(SeuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1.2<-ggplot(SeuObj@meta.data,aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
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
  	theme_classic() 

#p4<-VlnPlot(SeuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)

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

p8 <- ggplot(SeuObj@meta.data, aes(x=nCount_RNA)) + 
      geom_density(alpha = 0.4, fill="bisque")+labs(title = "distribution of nUMI")+
      theme_classic()#+

p9<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = nFeature_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nGene") #+scale_color_viridis_c(option="turbo")

#p9 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nFeature_RNA))+
#  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+
#  labs(fill="nGene",x=NULL,y=NULL,title = NULL)

p10<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = nCount_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nUMI") #+scale_color_viridis_c(option="turbo")

#p10 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nCount_RNA))+
#    geom_tile()+  theme_void()+ coord_fixed()+ scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+
#    labs(fill="read counts",x=NULL,y=NULL,title = "spatial distribution of read counts")

myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"),
     brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))


#SeuObj<-subset(SeuObj, subset = nFeature_RNA > 50 & nCount_RNA > 100 & percent.mt < 20)

message("SCTransform-RunPCA-RunUMAP-FindNeighbors-FindClusters:")
#SeuObj<-SCTransform(SeuObj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
SeuObj<-SCTransform(SeuObj, method = "glmGamPoi", verbose = FALSE)
SeuObj <- RunPCA(SeuObj, assay = "SCT", verbose = FALSE)
SeuObj <- RunUMAP(SeuObj, reduction = "pca", dims = 1:30)
SeuObj <- FindNeighbors(SeuObj, reduction = "pca", dims = 1:30)
SeuObj <- FindClusters(SeuObj, res=c(0.2, 0.4, 0.6, 1), verbose = FALSE)
SeuObj <- FindClusters(SeuObj, verbose = FALSE)

# add coord info as embedding
sp.embedding<-SeuObj@meta.data[,c("coord_x", "coord_y")]
colnames(sp.embedding)<-c("Spatial_1", "Spatial_2") # prefix with number
sp.embedding<-as.matrix(sp.embedding) # needs to be a matrix
SeuObj$spatial<-CreateDimReducObject(embeddings=sp.embedding, key='Spatial_', assay='SCT')

saveRDS(SeuObj, paste0(args[1], "_SCT.rds"))
message("RDS with sct counts saved.")

DefaultAssay(SeuObj) <- "RNA"
#deg <- FindAllMarkers(SeuObj,only.pos = TRUE)#, min.pct = 0.25, logfc.threshold = 0.5) # method: wilcox
#deg<-subset(deg, deg$p_val<1e-20)
#write.csv(deg,paste0(args[1], "_SCT.rds_allMarkers.csv"))

deg_roc <- FindAllMarkers(SeuObj,only.pos = TRUE, test.use="roc") # method: roc
write.csv(deg_roc,paste0(args[1], "_SCT.rds_allMarkers_roc.csv"))

#deg.top <- deg %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
deg.top <- deg_roc %>% group_by(cluster) %>% top_n(n = 3, wt = myAUC)
deg.topmarkers<-as.vector(deg.top$gene)

pdf(paste0(args[1], "_SCT.rds_CLUSTERS.pdf"),width=20,height=10)

pp3<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = seurat_clusters), shape=".")+theme_void()+coord_fixed()+scale_color_manual(values = myCol)

p0.2 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.2))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.2")

p0.4 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.4))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.4")

p0.6 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.6))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.6")

p1.0 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.1))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res1")

pp4 <- DimPlot(SeuObj, reduction = "umap", cols=myCol, label = TRUE) + NoLegend() + coord_fixed()
pp7<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols=myCol)
pp8<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nFeature_RNA", feature2 = "percent.mt", cols=myCol)
#pp9<-VlnPlot(SeuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)

grid.arrange(p9,p10,ncol = 2, top="Spatial distribution of nGene and nUMI")
grid.arrange(p6,p5,ncol = 3, top="nGene and nUMI per cell")
grid.arrange(p7,p8,p2,ncol = 3)
grid.arrange(p1.2,p3.2,ncol = 2)

#C = SeuObj@assays$RNA@counts
#C@x = C@x/rep.int(colSums(C), diff(C@p))
#most_expressed <- order(Matrix::rowSums(C), decreasing = T)[10:1]
#boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.05, las = 1, main="Top expressed genes", xlab = "% total count per cell",
#    col = (scales::hue_pal())(20)[10:1], horizontal = TRUE)

grid.arrange(pp4,pp3,ncol = 2, top="Seurat clustering")
#print(pp5)
#print(pp9)
grid.arrange(pp7,pp8,ncol = 2, top="Scatter plot grouped by clusters")
grid.arrange(p0.2,p0.4,ncol = 2, top="Seurat other resolutions")
grid.arrange(p0.6,p1.0,ncol = 2)
dev.off()

message("QC and clustering plots exported.")


p <- lapply(deg.topmarkers, function(x) {
tmexp<-as.matrix(SeuObj@assays$SCT@data[x,])  #@counts - unnormalised counts; @data normlised data (log or sct)
ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = tmexp[,1]), shape=".")+
    theme_void()+coord_fixed()+
    scale_color_gradient(low = "#fcf5eb", high = "#800080")+
    labs(title=x, color="")
})

ggsave(
   filename = paste0(args[1],"_SCT.rds_FEATURES.pdf"), 
   plot = marrangeGrob(p, nrow=3, ncol=3), 
   width = 15, height = 9
)

message("Feature spatial plots exported.")
