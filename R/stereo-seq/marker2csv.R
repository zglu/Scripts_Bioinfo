#Seurat markers to signature csv
# Find markers

args<-commandArgs(T)

if (length(args)<4) {
  stop("Usage: Rscript marker2csv.R [data.RDS] [group2test] [roc/wilcox/t] [No. top markers]")
}

library(Seurat)
library(dplyr)

scData<-args[1]
testGrp<-args[2]
testMethod<-args[3]
topGenes<-args[4]


sc1<-readRDS(scData)
sc1<-NormalizeData(sc1, assay = "RNA")

if (testGrp!=0){
  sc1<-SetIdent(sc1, value=sc1[[testGrp]])  
}

deg<-FindAllMarkers(sc1, assay="RNA", only.pos = TRUE, test.use = testMethod)


if (testGrp=="roc"){
  deg.top <- deg %>% group_by(cluster) %>% top_n(n = topGenes, wt = myAUC)
} else {
  deg.top <- deg %>% group_by(cluster) %>% top_n(n = topGenes, wt = avg_log2FC)
}

test<-deg.top[, c("cluster", "gene")] %>% group_by(cluster)

# add number for each marker within each group
test %>% group_by(cluster) %>% mutate(id=row_number()) -> test

# split gene column into multiple columns based on cluster
library(reshape2)
casted<-dcast(test, cluster ~ id, value.var="gene")
# remove colnames
colnames(casted)<-NULL
# transform  and save
write.csv(t(casted), file=paste0(scData, "_", testGrp, "_", testMethod, "_top", topGenes, ".csv"), quote=F, row.names = F, col.names = F, na="")

#############
# export as gene cluter list
test<-deg.top[, c("cluster", "gene")] %>% group_by(cluster)
geneCl<-aggregate(cluster~gene,data=test,FUN = function(x) paste0(x,collapse = ','))
#Actb   immune 2
write.table(geneCl, file=paste0(scData, "_", testGrp, "_", testMethod, "_top", topGenes, ".txt"), sep="\t", quote=F, row.names = F, col.names = F, na="")
