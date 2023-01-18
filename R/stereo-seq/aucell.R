# AUCell scoring method is ranking-based, AUCell is independent of the gene expression units and the normalization procedure
# Rscript aucell.R [seurat rds] [markers.csv]

args<-commandArgs(T)

if (length(args)<3) {
  stop("Usage: Rscript aucell.R [data.RDS] [markers.csv] [assay(SCT/RNA/integrated)]")
}

library(SCINA)
library(AUCell)
library(Seurat)
library(UCell)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggsci)

myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
      pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
			pal_jama("default")(7),pal_jco("default")(10),
			pal_locuszoom("default")(7),pal_startrek("uniform")(7),
			pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
			pal_simpsons("springfield")(16),
			pal_gsea("default")(12)))

geneSets<-preprocess.signatures(args[2])


SeuObj<-readRDS(args[1])
exprMatrix <- as.matrix(Seurat::GetAssayData(SeuObj))

# Build gene-expression rankings for each cell
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
#save(cells_rankings, file="cells_rankings.RData")

# Calculate enrichment for the gene signatures (AUC); by default using top 5%
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
#save(cells_AUC, file="cells_AUC.RData")

#               cells
#gene sets          168_380     241_234     240_229    173_365     208_295
#  Melanoma      0.03565703 0.008788077 0.022593670 0.01900125 0.031114542
#  Keratinocytes 0.54224874 0.370610889 0.491580458 0.46749226 0.534168995
#  CAFs          0.22316239 0.073846154 0.097435897 0.12341880 0.056495726

set.seed(333)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

#cells_assignment$Melanoma$aucThr$thresholds

# get cells assigned using automatic thresholds
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)

assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

#     cell  geneSet
#1 340_121 Melanoma
#2 312_167 Melanoma
#3 169_350 Melanoma


# one cell can be assigned to multiple cell types
# to summarize the assignments
assignments <- assignmentTable %>% dplyr::group_by(cell) %>% 
  dplyr::summarize(AUCellType = paste(geneSet, collapse = "/")) %>%
  dplyr::group_by(AUCellType) %>%
  dplyr::mutate(n = length(unique(cell))) # %>% dplyr::filter(n > 5)

head(assignments)
#  cell    AUCellType                 n
#  <chr>   <chr>                  <int>
#1 154_336 CAFs/Blood.vessels        85
#2 154_338 Immune.cells            2247
#3 154_339 Immune.cells/Pericytes   162

# add to meta data, NA if no match
SeuObj@meta.data$AUCellType<-assignments$AUCellType[match(rownames(SeuObj@meta.data),assignments$cell)]
p.aucell<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= AUCellType))+
	geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)


##### UCell scoring
#UCell.scores <- ScoreSignatures_UCell(exprMatrix, features=geneSets)
#head(UCell.scores)
##        Melanoma_UCell Keratinocytes_UCell CAFs_UCell Immune.cells_UCell
##211_266     0.00000000           0.2902941 0.17227451         0.00000000
##288_181     0.03400000           0.3423922 0.00000000         0.00000000
#
#melted <- reshape2::melt(UCell.scores)
#colnames(melted) <- c("Cell", "Signature", "UCell_score")
#p <- ggplot(melted, aes(x = Signature, y = UCell_score)) + geom_violin(aes(fill = Signature),
#    scale = "width") + geom_boxplot(width = 0.1, outlier.size = 0) + theme_bw() +
#    theme(axis.text.x = element_blank())
#p

# directly work on Seurat Object
SeuObj<-AddModuleScore_UCell(SeuObj, features=geneSets, name=NULL)
head(SeuObj@meta.data)
#        SCT_snn_res.0.8   Melanoma Keratinocytes       CAFs Immune.cells
#211_266               0 0.00000000     0.2902941 0.17227451   0.00000000
#288_181               6 0.03400000     0.3423922 0.00000000   0.00000000

pdf("UCell_sig.pdf", width=25, height=15)
p.ucell<-FeaturePlot(SeuObj, reduction="spatial", features=names(geneset), ncol=4, coord.fixed=TRUE, pt.size=0.01)#keep.scale=all/feature/NULL
dev.off()
