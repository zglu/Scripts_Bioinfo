# Rscript scina.R [SCT rds] [marker tsv/gmt]


args<-commandArgs(T)

if (length(args)<2) {
  stop("Usage: Rscript scina.R [SCT rds] [marker tsv/gmt]")
}

suppressMessages(library(Seurat))
suppressMessages(require(dplyr))
library(ggplot2)
library(ggsci)
suppressMessages(library(gridExtra))

myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
      pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
			pal_jama("default")(7),pal_jco("default")(10),
			pal_locuszoom("default")(7),pal_startrek("uniform")(7),
			pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
			pal_simpsons("springfield")(16),
			pal_gsea("default")(12)))


message("Reading RDS:")
SeuObj<-readRDS(args[1])


# change parameters
#SeuObj <- RunUMAP(SeuObj, dims = 1:20, verbose = F)
#SeuObj <- FindNeighbors(SeuObj, dims = 1:20, verbose = FALSE)
SeuObj <- FindClusters(SeuObj, verbose = FALSE, res=0.25)#, graph.name = "SCT_snn"


# use markers from tsv or gmt
library(msigdbi)
markers <- msigdbi::read.gmt(args[2])

# use markers from csv 

#Endothelial.cells	Fibroblasts	Hepatocytes
#Mmrn2	Nrxn1	Fabp1
#Ptprb	Carmn	Apoa2
#Cldn5	Gli3	Mup20

library(SCINA)
#markers<-preprocess.signatures(args[2])

#load("creLPS_exprMatrix.Rdata") #-->SeuObj.exprMatrix
SeuObj.exprMatrix <- as.matrix(Seurat::GetAssayData(SeuObj))
# rm_overlap=TRUE
SeuObj.scina = SCINA(exp = SeuObj.exprMatrix, signatures = markers, rm_overlap = FALSE, allow_unknown = TRUE)#signatures=markers$genesets
SeuObj@meta.data$SCINA <- SeuObj.scina$cell_labels

write.csv(SeuObj@meta.data, file=paste0(args[1],"_scina.csv"))

##### plot
pdf(paste0(args[1], "_SCINA.pdf"), width=8, height=4)

dim1<-DimPlot(SeuObj, reduction = "umap", cols=myCol, label = TRUE) + NoLegend() + coord_fixed()
seuratSpatial<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= seurat_clusters))+ 
  geom_tile()+
  theme_void()+
  coord_fixed()+ 
  scale_fill_manual(values = myCol)+
  labs(fill="clusters")

dim2<-DimPlot(SeuObj, reduction = "umap", cols=myCol, group.by='SCINA', label = FALSE) + coord_fixed()

scinaSpatial<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCINA))+ 
  geom_tile()+
  theme_void()+
  coord_fixed()+ 
  scale_fill_manual(values = myCol)+
  labs(fill="SCINA clusters")

scina_seurat<-ggplot(SeuObj@meta.data, aes(x=factor(seurat_clusters), fill=factor(SCINA)))+geom_bar()+theme_classic()+
	scale_fill_manual(values = myCol)+
	ggtitle("Proportion of predicted cell types in WT-PBS")+
	labs(x="Seurat cluster", y="Cell number", fill="SCINA")


# SCINA annotation plot
scina_heatmap<-plotheat.SCINA(SeuObj.exprMatrix, SeuObj.scina, markers)

grid.arrange(dim2,scinaSpatial,ncol = 2)
grid.arrange(dim1,seuratSpatial,ncol = 2)
grid.arrange(scina_seurat,ncol = 2)
print(scina_heatmap)

dev.off()
message("SCINA plots exported.")

