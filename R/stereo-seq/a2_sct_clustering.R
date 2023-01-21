# Rscript a2_sct_clustering.R [rds file]
# output: SCtransformed RDS + cluster plot + feature spatial plot

library(Seurat)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(RColorBrewer)
library(dplyr)

args<-commandArgs(T)

myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
      pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
			pal_jama("default")(7),pal_jco("default")(10),
			pal_locuszoom("default")(7),pal_startrek("uniform")(7),
			pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
			pal_simpsons("springfield")(16),
			pal_gsea("default")(12)))

message("Reading RDS:")
SeuObj<-readRDS(args[1])

# write files for SpatialDE, HotSpot, SCENIC, etc


# Filtering
# cell-level
#SeuObj <- subset(SeuObj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &  nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 3)

# gene-level: keep genes with count >0 in at least 10 cells
#counts <- GetAssayData(object = SeuObj, slot = "counts")
#nonzero <- counts>0
#keep_genes<-Matrix::rowSums(nonzero) >=10
#filtered_counts<-counts[keep_genes,]
#SeuObj_fil<-CreateSeuratObject(filtered_counts, meta.data = SeuObj@meta.data)
#saveRDS(SeuObj_fil, file=paste0(args[1], "_filtered.rds"))
#SeuObj<-readRDS(file=paste0(args[1], "_filtered.rds"))


message("Log-normalize on RNA counts: ")
SeuObj %>%  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = rownames(SeuObj)) -> SeuObj

message("SCTransform-RunPCA-RunUMAP-FindNeighbors-FindClusters:")

#SeuObj %>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE) %>%
#  RunPCA(verbose = FALSE,assay="SCT") %>%
#  RunUMAP( dims = 1:30, verbose = FALSE)%>%
#  FindNeighbors(dims = 1:30, verbose = FALSE)%>%
#  FindClusters(verbose = FALSE) -> SeuObj # res=2

# use glmGamPoi package to substantially improves the speed of the learning procedure
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("glmGamPoi")
SeuObj<-SCTransform(SeuObj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
SeuObj <- RunPCA(SeuObj, assay = "SCT", verbose = FALSE)
SeuObj <- RunUMAP(SeuObj, reduction = "pca", dims = 1:30)
SeuObj <- FindNeighbors(SeuObj, reduction = "pca", dims = 1:30)
SeuObj <- FindClusters(SeuObj, verbose = FALSE)

# save data
saveRDS(SeuObj,paste0(args[1],"_SCT.rds"))
message("RDS with sct counts saved.")

SeuObj

message("Number of cells per cluster:")
table(SeuObj@active.ident)

DefaultAssay(SeuObj) <- "RNA"
#deg <- FindAllMarkers(SeuObj,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) # method: wilcox
#deg<-subset(deg, deg$p_val<1e-20)
#write.csv(deg,paste0(args[1],"_SCT.rds_allMarkers.csv"))

deg_roc <- FindAllMarkers(SeuObj,only.pos = TRUE, test.use="roc") # method: roc
write.csv(deg_roc,paste0(args[1],"_SCT.rds_allMarkers_roc",".csv"))

message("Marker genes exported.")

deg.top <- deg %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
deg.topmarkers<-as.vector(deg.top$gene)

## plotting 
pdf(paste0(args[1], "_SCT.rds_CLUSTERS.pdf"),width=10,height=5)

# spatial feature plot
p1<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nFeature_RNA))+
  geom_tile()+
  theme_void()+
  coord_fixed()+
  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
  labs(fill="nGene",x=NULL,y=NULL,title = NULL)+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+ 
  theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))

p2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nCount_RNA))+
  geom_tile()+
  theme_void()+
  coord_fixed()+
  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
  labs(fill="nUMI",x=NULL,y=NULL,title = NULL)+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+ 
  theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))


p3 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= seurat_clusters))+ 
  geom_tile()+
  theme_void()+
  coord_fixed()+ 
  scale_fill_manual(values = myCol)+
  labs(fill="Seurat clusters")

p4 <- DimPlot(SeuObj, reduction = "umap", cols=myCol, label = TRUE) + NoLegend() + coord_fixed()

#p3 <- DoHeatmap(SeuObj,features = deg.top$gene)

p5 <- ElbowPlot(SeuObj)
p7<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols=myCol)
p8<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nFeature_RNA", feature2 = "percent.mt", cols=myCol)
p9<-VlnPlot(SeuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)

grid.arrange(p4,p3,ncol = 2, top="Seurat clustering")
#print(p5)
print(p9)
grid.arrange(p7,p8,ncol = 2, top="Scatter plot grouped by clusters")
grid.arrange(p1,p2,ncol = 2, top="Spatial distribution of nGene and nUMI")
dev.off()

message("Cluster plots exported.")


pdf(paste0(args[1], "_SCT.rds_DOTPLOT.pdf"),width=10,height=5)
p6<-DotPlot(SeuObj, features = deg.topmarkers, cols = c("lightgrey", "red")) + NoLegend() + RotatedAxis()
print(p6)
dev.off()

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
   filename = paste0(args[1],"_SCT.rds_FEATURES.pdf"), 
   plot = marrangeGrob(p, nrow=3, ncol=3), 
   width = 10, height = 5
)

message("Feature spatial plots exported.")
