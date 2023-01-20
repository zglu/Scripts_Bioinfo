# calculate module scores using Seurat addModuleScore, AUCell, and UCell 

args<-commandArgs(T)

if (length(args)<3) {
  stop("Usage: Rscript Seurat_markerScores.R [data.RDS] [markers.csv] [assay(SCT/RNA/integrated)]")
}

suppressPackageStartupMessages({
  library(RColorBrewer)
  library(Seurat)
  library(gridExtra)
  library(ggplot2)
  library(SCINA)
  library(AUCell)
  library(UCell)
  library(dplyr)
  library(SingleCellExperiment)
})

myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"),
     brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))

SeuObj<-readRDS(args[1])

# AUCell and UCell are ranking-based, independent of normalisation
#DefaultAssay(SeuObj)<- args[3]
DefaultAssay(SeuObj)<- "RNA" #or SCT

#### AUCell ####
message("Cell assignment using AUCell:")
geneSets<-preprocess.signatures(args[2])
exprMatrix <- as.matrix(Seurat::GetAssayData(SeuObj))
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
#save(cells_rankings, file="cells_rankings.RData")
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
#save(cells_AUC, file="cells_AUC.RData")

# add scores to metadata
auc.scores<-DataFrame(t(assay(cells_AUC, "AUC")))
colnames(auc.scores)<-paste0(colnames(auc.scores), "_AUC")
SeuObj@meta.data<-cbind(SeuObj@meta.data, auc.scores)

#names(geneSets)
p2<-lapply(rownames(cells_AUC), function(i) {
    FeaturePlot(SeuObj, reduction="spatial", features = paste0(i, "_AUC"), label = FALSE, repel = TRUE) + coord_fixed()+theme_void()+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+labs(color="AUCellScore")
})

# error with ggplot2
#names<-paste0(rownames(cells_AUC), "_AUC")
#p2<-lapply(names, function(i) {
#    ggplot(SeuObj@meta.data, aes(coord_x, coord_y, fill=i))+geom_tile() + coord_fixed()+theme_void()+scale_fill_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+labs(fill="AUCellScore")
#})

ggsave(
   filename = paste0(args[1], "_AUCellScores.pdf"),
   plot = marrangeGrob(p2, nrow=2, ncol=4),
   width = 25, height = 15
)

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE)
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)

# summarise assignment
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
assignments <- assignmentTable %>% dplyr::group_by(cell) %>% 
  dplyr::summarize(AUCellType = paste(geneSet, collapse = "/")) %>%
  dplyr::group_by(AUCellType) %>%
  dplyr::mutate(n = length(unique(cell))) # %>% dplyr::filter(n > 5)
# add to meta data, NA if no match
SeuObj@meta.data$AUCellType<-assignments$AUCellType[match(rownames(SeuObj@meta.data),assignments$cell)]

pdf(paste0(args[1], "_AUCell_thresholds_assignments.pdf"), width=15, height=9)
set.seed(333)
par(mfrow=c(3,4))
thresholds <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE)
#p.aucell<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= AUCellType))+geom_tile()+ theme_void()+ coord_fixed()#+ scale_fill_manual(values = myCol)
#print(p.aucell)
dev.off()

#### UCell ####
message("Calculating UCell scores:")
SeuObj<-AddModuleScore_UCell(SeuObj, features=geneSets) # name=NULL
#p.ucell<-FeaturePlot(SeuObj, reduction="spatial", features=names(geneSets), ncol=4, coord.fixed=TRUE, pt.size=0.01)#keep.scale=all/feature/NULL
p3<-lapply(names(geneSets), function(i) {
    FeaturePlot(SeuObj, reduction="spatial", features = paste0(i, "_UCell"), label = FALSE, repel = TRUE) + coord_fixed()+theme_void()+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "RdPu"))+labs(color="UCellScore")
})

ggsave(
   filename = paste0(args[1], "_UCellScores.pdf"),
   plot = marrangeGrob(p3, nrow=2, ncol=4),
   width = 25, height = 15
)

#write.csv(SeuObj@meta.data, file=paste0(args[1], "_sigScores.csv"))
saveRDS(SeuObj, file=paste0(args[1], "_sigScores.rds"))

#### SCINA ####
message("Assignment using SCINA:")
SeuObj.scina = SCINA(exp = exprMatrix, signatures = geneSets, rm_overlap = FALSE, allow_unknown = TRUE)#signatures=markers$genesets
SeuObj@meta.data$SCINA <- SeuObj.scina$cell_labels
SeuObj@meta.data$scina_prob <- SeuObj.scina$probabilities

pdf(paste0(args[1], "_SCINA_cellTypes.pdf"), width=15, height=9)
ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCINA))+ 
   geom_tile()+ theme_void()+ coord_fixed()+ 
   scale_fill_manual(values = myCol)+
   labs(fill="SCINA clusters")
dev.off()

#### Seurat addModuleScore ####
message("Calculating scores using addModuleScore:")
# Seurat calculate the average expression of each cluster, thus a normalisation will needed
# SCT returns 3000 HVGs
DefaultAssay(SeuObj)<- args[3] #SCT / integrated

#allgenes<-rownames(SeuObj@assays$integrated@data) #SCT / integrated / RNA
allgenes<-rownames(SeuObj)

markers<-read.csv(args[2], header=T)
#   Astrocyte CA1_Pyramidal Endothelial Ependymal Interneuron  Microglia
#1       Aass       Abhd17c       Abcb1     Abca4  Ac009690.3      Abca9
#2      Abcb9          Abi1       Abcc6    Acad10  Ac073111.5      Abcb1

p1 <- lapply(colnames(markers), function(i) {
	m <- markers[, i]
	m <- intersect(m, allgenes)
	#print(i)
	SeuObj<-AddModuleScore(SeuObj, features=list(m), name=i)
	FeaturePlot(SeuObj, reduction="spatial", features = paste0(i, "1"), label = FALSE, repel = TRUE) + coord_fixed()+theme_void()+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+labs(color="SeuratScore")
})

ggsave(
   filename = paste0(args[1], "_SeuratScores.pdf"),
   plot = marrangeGrob(p1, nrow=2, ncol=4),
   width = 25, height = 15
)



### combine output pdfs
library(qpdf)
qpdf::pdf_combine(input = c(paste0(args[1], "_SeuratScores.pdf"), paste0(args[1], "_UCellScores.pdf"), paste0(args[1], "_AUCellScores.pdf"),paste0(args[1], "_AUCell_thresholds_assignments.pdf")), output = paste0(args[1], "_sigScores.pdf"))


