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
#SeuObj <- subset(SeuObj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &  nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 10)


message("Log-normalize on RNA counts: ")
SeuObj %>%  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() -> SeuObj #features = rownames(SeuObj)

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

message("Dim:")
dim(SeuObj@assays$RNA@counts)

message("Number of cells per cluster:")
table(SeuObj@active.ident)

DefaultAssay(SeuObj) <- "RNA"
SeuObj %>%  NormalizeData()%>%
  FindVariableFeatures()%>%
  ScaleData(features = rownames(SeuObj))

deg <- FindAllMarkers(SeuObj,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) # method: wilcox
deg<-subset(deg, deg$p_val<1e-20)
write.csv(deg,paste0(args[1],"_SCT.rds_allMarkers.csv"))

deg_roc <- FindAllMarkers(SeuObj,only.pos = TRUE, test.use="roc") # method: roc
write.csv(deg_roc,paste0(args[1],"_SCT.rds_allMarkers_roc",".csv"))

message("Marker genes exported.")

deg.top <- deg %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
deg.topmarkers<-as.vector(deg.top$gene)

## plotting 
pdf(paste0(args[1], "_SCT.rds_CLUSTERS.pdf"),width=10,height=5)

# spatial feature plot
#p1<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nFeature_RNA))+
#  geom_tile()+
#  theme_void()+
#  coord_fixed()+
#  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
#  labs(fill="gene counts",x=NULL,y=NULL,title = "nFeature_RNA")+
#  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+ 
#  theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))
#
#p2<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nCount_RNA))+
#  geom_tile()+
#  theme_void()+
#  coord_fixed()+
#  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
#  labs(fill="read counts",x=NULL,y=NULL,title = "nCount_RNA")+
#  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+ 
#  theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))


# spatial cluster plot
p3<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = seurat_clusters), shape=".")+theme_void()+coord_fixed()+scale_color_manual(values = myCol) # 

#p3 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= seurat_clusters))+ 
#  geom_tile()+
#  theme_void()+
#  coord_fixed()+ 
#  scale_fill_manual(values = myCol)+
#  labs(fill="clusters")

p4 <- DimPlot(SeuObj, reduction = "umap", cols=myCol, label = TRUE) + NoLegend() + coord_fixed()


p5 <- ElbowPlot(SeuObj)
p6<-DotPlot(SeuObj, features = deg.topmarkers, cols = c("lightgrey", "red")) + NoLegend()+RotatedAxis()
p7<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols=myCol)
p8<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nFeature_RNA", feature2 = "percent.mt", cols=myCol)
p9<-VlnPlot(SeuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)


grid.arrange(p4,p3,ncol = 2)
print(p6)
print(p5)
grid.arrange(p7,p8,ncol = 2)
print(p9)
dev.off()

message("Cluster plots exported.")

##### top markers plots
#p <- lapply(deg.topmarkers, function(x) {
#tmexp<-as.matrix(SeuObj@assays$SCT@counts[x,])  
#ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= tmexp[,1]))+ 
#  geom_tile()+
#  theme_void()+
#  coord_fixed()+ 
#  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
#  labs(title=x, fill="Counts")
#})

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
