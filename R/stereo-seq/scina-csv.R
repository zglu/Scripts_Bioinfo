library(SCINA)

guill.allDEGs<-read.csv("~/spatial/analysis-test/VIB_mouse_liver/Guilliams_allDEGs.csv")

#Endothelial cells	Mmrn2	Ptprb	Cldn5
#Stromal cells	Nrxn1	Carmn	Gli3
#Hepatocytes	Fabp1	Apoa2	Mup20

guill.subDEGs<-guill.allDEGs[,c(1,3:52)]

newMarkers<-t(as.matrix(guill.subDEGs))

write.csv(newMarkers, file="Guilliams_top50Markers.csv", col.names=FALSE)
####### remove first row!!!


#Endothelial.cells	Fibroblasts	Hepatocytes
#Mmrn2	Nrxn1	Fabp1
#Ptprb	Carmn	Apoa2
#Cldn5	Gli3	Mup20

markers<-preprocess.signatures("Guilliams_top50Markers.csv")

##################

library(Seurat)
SeuObj<-readRDS(args[1])
#SeuObj <- FindClusters(SeuObj, verbose = FALSE, res=0.25)

SeuObj.exprMatrix <- as.matrix(Seurat::GetAssayData(SeuObj))

# save expr matrix as Rdata
save(SeuObj.exprMatrix, file="creLPS_exprMatrix.Rdata")

####
# load matrix objective
load("creLPS_exprMatrix.Rdata")
ls() #SeuObj.exprMatrix
