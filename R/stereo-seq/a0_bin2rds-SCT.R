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
library(dplyr)
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


myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
      pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
			pal_jama("default")(7),pal_jco("default")(10),
			pal_locuszoom("default")(7),pal_startrek("uniform")(7),
			pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
			pal_simpsons("springfield")(16),
			pal_gsea("default")(12)))

#myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"),
#     brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))

#saveRDS(SeuObj, paste0(args[1], "_BIN", binsize, ".rds"))

message("Summary of nFeature")
summary(SeuObj@meta.data$nFeature_RNA)

message("Summary of nCount")
summary(SeuObj@meta.data$nCount_RNA)

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

# filtering empty bins
SeuObj<-subset(SeuObj, subset = nFeature_RNA > 0 & nCount_RNA > 0)


message("Log-normalize on RNA counts: ")
SeuObj %>%  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = rownames(SeuObj)) -> SeuObj

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

#DimPlot(SeuObj, reduction="umap")
#DimPlot(SeuObj, reduction="spatial")

saveRDS(SeuObj, paste0(args[1], "_BIN", binsize, "_SCT.rds"))
message("RDS with sct data saved.")


message("Number of cells per cluster:")
table(SeuObj@active.ident)

DefaultAssay(SeuObj) <- "RNA"
#deg <- FindAllMarkers(SeuObj,only.pos = TRUE)#, min.pct = 0.25, logfc.threshold = 0.5) # method: wilcox
#deg<-subset(deg, deg$p_val<1e-20)
#write.csv(deg,paste0(args[1], "_BIN", binsize, "_SCT.rds_allMarkers.csv"))

deg_roc <- FindAllMarkers(SeuObj,only.pos = TRUE, test.use="roc") # method: roc
write.csv(deg_roc,paste0(args[1],"_BIN", binsize, "_SCT.rds_allMarkers_roc",".csv"))

message("Marker genes exported.")

#deg.top <- deg %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
deg.top <- deg_roc %>% group_by(cluster) %>% top_n(n = 3, wt = myAUC)
deg.topmarkers<-as.vector(deg.top$gene)

## plotting
pdf(paste0(args[1], "_BIN", binsize, "_SCT.rds_CLUSTERS.pdf"),width=10,height=5)

pp3 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= seurat_clusters))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat clusters")

pp4 <- DimPlot(SeuObj, reduction = "umap", cols=myCol, label = TRUE) + NoLegend() + coord_fixed()

p0.2 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.2))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.2")

p0.4 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.4))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.4")

p0.6 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.6))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.6")

p1.0 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.1))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res1")

#pp5 <- ElbowPlot(SeuObj)
pp7<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols=myCol)
pp8<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nFeature_RNA", feature2 = "percent.mt", cols=myCol)
pp9<-VlnPlot(SeuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)

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

grid.arrange(pp4,pp3,ncol = 2, top="Seurat clustering")
#print(pp5)
#print(pp9)
grid.arrange(pp7,pp8,ncol = 2, top="Scatter plot grouped by clusters")
grid.arrange(p0.2,p0.4,ncol = 2, top="Seurat other resolutions")
grid.arrange(p0.6,p1.0,ncol = 2)
dev.off()

message("QC and Cluster plots exported.")

##### top markers plots
p <- lapply(deg.topmarkers, function(x) {
tmexp<-as.matrix(SeuObj@assays$SCT@data[x,]) #@counts - unnormalised counts; @data normlised data (log or sct)
ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= tmexp[,1]))+
  geom_tile()+
  theme_void()+
  coord_fixed()+
  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
  labs(title=x, fill="") # sct-normalised data
})

ggsave(
   filename = paste0(args[1],"_BIN", binsize,"_SCT.rds_FEATURES.pdf"),
   plot = marrangeGrob(p, nrow=3, ncol=3),
   width = 10, height = 5
)

message("Feature spatial plots exported.")

## dotplot
pdf(paste0(args[1], "_BIN", binsize, "_SCT.rds_DOTPLOT.pdf"),width=10,height=5)
#pheat <- DoHeatmap(SeuObj,features = deg.top$gene)
pp6<-DotPlot(SeuObj, features = deg.topmarkers, cols = c("lightgrey", "red")) + NoLegend() + RotatedAxis()
print(pp6)
#print(pheat)
dev.off()

message("Dotplot plots exported.")

