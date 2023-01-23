# calculate module scores 

library(RColorBrewer)
library(Seurat)
library(gridExtra)
library(ggplot2)

args<-commandArgs(T)

SeuObj<-readRDS(args[1]) #("Integrated_SCT-RPCA.rds")

markers<-read.csv("Zeisel2015_markers.csv", header=T)

#   Astrocyte CA1_Pyramidal Endothelial Ependymal Interneuron  Microglia
#1       Aass       Abhd17c       Abcb1     Abca4  Ac009690.3      Abca9
#2      Abcb9          Abi1       Abcc6    Acad10  Ac073111.5      Abcb1

allgenes<-rownames(SeuObj@assays$integrated@data) #SCT / integrated / RNA

Astrocyte<-markers$Astrocyte
CA1_Pyramidal<-markers$CA1_Pyramidal
Endothelial<-markers$Endothelial
Ependymal<-markers$Ependymal
Interneuron<-markers$Interneuron
Microglia<-markers$Microglia
Mural<-markers$Mural
Oligodendrocyte<-markers$Oligodendrocyte
S1_Pyramidal<-markers$S1_Pyramidal

Astrocyte<-intersect(Astrocyte, allgenes)
CA1_Pyramidal<-intersect(CA1_Pyramidal, allgenes)
Endothelial<-intersect(Endothelial, allgenes)
Ependymal<-intersect(Ependymal, allgenes)
Interneuron<-intersect(Interneuron, allgenes)
Microglia<-intersect(Microglia, allgenes)
Mural<-intersect(Mural, allgenes)
Oligodendrocyte<-intersect(Oligodendrocyte, allgenes)
S1_Pyramidal<-intersect(S1_Pyramidal, allgenes)

SeuObj<-AddModuleScore(SeuObj, features=list(Astrocyte), name="Astrocyte")
SeuObj<-AddModuleScore(SeuObj, features=list(CA1_Pyramidal), name="CA1_Pyramidal")
SeuObj<-AddModuleScore(SeuObj, features=list(Endothelial), name="Endothelial")
SeuObj<-AddModuleScore(SeuObj, features=list(Ependymal), name="Ependymal")
SeuObj<-AddModuleScore(SeuObj, features=list(Interneuron), name="Interneuron")
SeuObj<-AddModuleScore(SeuObj, features=list(Microglia), name="Microglia")
SeuObj<-AddModuleScore(SeuObj, features=list(Mural), name="Mural")
SeuObj<-AddModuleScore(SeuObj, features=list(Oligodendrocyte), name="Oligodendrocyte")
SeuObj<-AddModuleScore(SeuObj, features=list(S1_Pyramidal), name="S1_Pyramidal")


write.csv(SeuObj@meta.data, file=paste0(args[1], "_markerScores.csv"))

p1<-FeaturePlot(SeuObj, features = 'Astrocyte1', label = TRUE, repel = TRUE) + coord_fixed() + scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlGnBu"))
p2<-FeaturePlot(SeuObj, features = 'CA1_Pyramidal1', label = TRUE, repel = TRUE) + coord_fixed()+ scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlGnBu"))
p3<-FeaturePlot(SeuObj, features = 'Endothelial1', label = TRUE, repel = TRUE) + coord_fixed()+ scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlGnBu"))
p4<-FeaturePlot(SeuObj, features = 'Ependymal1', label = TRUE, repel = TRUE) + coord_fixed()+ scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlGnBu"))
p5<-FeaturePlot(SeuObj, features = 'Interneuron1', label = TRUE, repel = TRUE) + coord_fixed()+ scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlGnBu"))
p6<-FeaturePlot(SeuObj, features = 'Microglia1', label = TRUE, repel = TRUE) + coord_fixed()+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds"))
p7<-FeaturePlot(SeuObj, features = 'Mural1', label = TRUE, repel = TRUE) + coord_fixed()+ scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlGnBu"))
p8<-FeaturePlot(SeuObj, features = 'Oligodendrocyte1', label = TRUE, repel = TRUE) + coord_fixed() +scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlGnBu"))
p9<-FeaturePlot(SeuObj, features = 'S1_Pyramidal1', label = TRUE, repel = TRUE) +coord_fixed() + scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlGnBu"))

pdf(paste0(args[1], "_markerScores.pdf", width=18, height=12))
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3)
dev.off()

# get No. of markers in assay
#print(args[1])
message("Astrocyte markers present in assay: ", length(Astrocyte))
message("CA1_Pyramidal markers present in assay: ", length(CA1_Pyramidal))
message("Endothelial markers present in assay: ", length(Endothelial))
message("Ependymal markers present in assay: ", length(Ependymal))
message("Interneuron markers present in assay: ", length(Interneuron))
message("Microglia markers present in assay: ", length(Microglia))
message("Mural markers present in assay: ", length(Mural))
message("Oligodendrocyte markers present in assay: ", length(Oligodendrocyte))
message("S1_Pyramidal markers present in assay: ", length(S1_Pyramidal))


# find markers
DefaultAssay(SeuObj)<-"RNA"

SeuObj%>% NormalizeData() --> SeuObj

c13.markers<-FindMarkers(SeuObj, ident.1=13, only.pos=TRUE, test.use="roc")
c7.markers<-FindMarkers(SeuObj, ident.1=7, only.pos=TRUE, test.use="roc", logfc.threshold = 0.5, min.pct = 0.25)

#subcluster 13 in the whole
subcluster13<-FindSubCluster(SeuObj, "13", graph.name="integrated_snn", subcluster.name = "sub13", resolution = 0.3, algorithm = 1)
head(subcluster13@meta.data)
DimPlot(subcluster13, group.by = "sub13", label = TRUE)+coord_fixed()


# take cluster 13 as a subset
c13<-subset(SeuObj, idents=13)
c13<-FindNeighbors(c13, dims=1:20) #k.param = 5
c13<-FindClusters(c13, res=0.15) #0.15
table(c13@active.ident)

c13.spatial <- ggplot(c13@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= seurat_clusters))+geom_tile()+ theme_void()+ coord_fixed() +labs(fill="C13 clusters")

pdf("c13_clusters.pdf")
DimPlot(c13, label=TRUE)+coord_fixed()
dev.off()

c13.allmarkers<-FindAllMarkers(c13, only.pos=TRUE, test.use="roc")


# average expression
avgexpr<-AverageExpression(c13, assays="SCT", slot="data")

AverageExpression(c13, features=c("Ncf2", "Fgd3"), slot="data", assays="SCT")

# feature plot
pdf("c13_featureplot.pdf", width=15, height=10)
# default slot RNA
FeaturePlot(c13, features=c("Ncf2", "Fgd3", "Fgd2","Rdx","Rcbtb1","Ctnnbip1","Il10ra","Pld4","Fbln7", "Lrrc7","Tjp1","Ctsd", "Dock2", "Ccl3", "Itgb2"), ncol=3, slot="data") & scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"), values=c(0, 60)) 
#FeaturePlot(c13, features=c("Ncf2", "Fgd3", "Fgd2","Rdx","Rcbtb1","Ctnnbip1","Il10ra","Pld4","Fbln7", "Lrrc7","Tjp1","Ctsd", "Dock2", "Ccl3", "Itgb2"), ncol=3, slot="data") & scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdPu"))) # & coord_fixed()


## feature spatial plot

myFeatures<-c("Ncf2", "Fgd2", "Il10ra", "Pld4", "Ctsd", "Dock2", "Ccl3", "Itgb2")
p <- lapply(myFeatures, function(x) {
tmexp<-as.matrix(c13@assays$RNA@data[x,])  #@counts - unnormalised counts; @data normlised data (log or sct)
ggplot(c13@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = tmexp[,1]), shape=".")+
    theme_void()+coord_fixed()+
    scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+
    labs(title=x, color="")
})

ggsave(
   filename = paste0("C13","_featureSpatial.pdf"), 
   plot = marrangeGrob(p, nrow=3, ncol=3), 
   width = 15, height = 9
)


## SCINA to identify microglia
SeuObj<-readRDS("Integrated_SCT-RPCA.rds")
DefaultAssay(SeuObj)<-"RNA"

SeuObj%>% NormalizeData() --> SeuObj

geneSets<-preprocess.signatures("Zeisel2015_markers.csv")
exprMatrix <- as.matrix(Seurat::GetAssayData(SeuObj))

SeuObj.scina = SCINA(exp = exprMatrix, signatures = geneSets, rm_overlap = FALSE, allow_unknown = TRUE)
SeuObj@meta.data$SCINA <- SeuObj.scina$cell_labels

table(SeuObj@meta.data$seurat_clusters[which(SeuObj@meta.data$SCINA=="Microglia")])

microglia<-subset(SeuObj, subset=SCINA=="Microglia")
microglia<-FindNeighbors(microglia, dims=1:20) #k.param = 5
microglia<-FindClusters(microglia, res=0.15)
saveRDS(microglia, "scina_Microglia_fromIntegrated.RDS")

microglia.allmarkers<-FindAllMarkers(microglia, only.pos=TRUE, test.use="roc")


myCol<-unique(c(brewer.pal(8, "Accent"), brewer.pal(12, "Set3"),brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), 
     brewer.pal(8, "Set1")))

pdf(paste0("Integrated", "_SCINA_cellTypes.pdf"), width=15, height=9)
DimPlot(SeuObj, reduction = "umap", cols=myCol, group.by='SCINA', label = FALSE) + coord_fixed()
DimPlot(microglia, label=FALSE)+coord_fixed()
dev.off()


