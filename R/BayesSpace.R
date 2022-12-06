# script to run BayesSpace for clustering

args<-commandArgs(T)
if (length(args)<2) {
  stop("Usage: Rscript BayesSpace.R [rds] [n_clusters]")
}

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

SeuObj<-readRDS(args[1])
n_clusters<-as.numeric(args[2])
n_pcs<-10

colData<-SeuObj@meta.data[,c("coord_x", "coord_y")]
colnames(colData)<-c("row", "col")

sce<-SingleCellExperiment(assays=list(counts=as(SeuObj@assays$RNA@counts, "dgCMatrix")), rowData=rownames(SeuObj@assays$RNA@counts), colData=colData)

set.seed(102)

# keeping the top n.PCs principal components

sce <- spatialPreprocess(sce, n.PCs=n_pcs, n.HVGs=2000, log.normalize=TRUE) #platform="ST"

sce <- qTune(sce, qs=seq(3, 15), d=n_pcs)

set.seed(149)
sce <- spatialCluster(sce, q=n_clusters, d=n_pcs, init.method="mclust", model="t", gamma=2, nrep=2000, burn.in=100, save.chain=TRUE)

write.table(colData(sce), file=paste0("bayesSpace_clusters_", n_clusters, ".txt"), quote=F)

#clusterPlot(sce)

# plot using ggplot
library(ggplot2)
library(RColorBrewer)
data<-read.delim(paste0("bayesSpace_clusters_", n_clusters, ".txt"), header=T, sep=" ")

myCol<-unique(c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"),
     brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))

p<-ggplot(data, aes(row, col, fill=as.factor(spatial.cluster)))+geom_tile()+theme_void()+scale_fill_manual(values = myCol)+labs(fill="BayesSpace clusters")+coord_fixed()

pdf(paste0("BayesSpace_clusters", n_clusters, ".pdf"))
print(p)
qPlot(sce)
dev.off()


# enhanced resolution
sce.enhanced <- spatialEnhance(sce, q=n_clusters, d=n_pcs, model="t", gamma=2, jitter_prior=0.3, jitter_scale=3.5, nrep=2000, burn.in=100, save.chain=TRUE)

write.table(colData(sce.enhanced), file="bayesSpace_clusters_enhanced.txt")

pdf(file=paste0("bayesSpace_clusters_enhanced_", n_clusters, ".pdf"))
clusterPlot(sce.enhanced)
dev.off()


### add bayesspace to Seurat objective
#SeuObj<-readRDS("683TL_F2_RGP_L20220920006.bin1.txt_BIN50_SCT.rds")
#bayes<-read.table("683TL_F2_RGP_bayesSpace_clusters_8.txt", header=T, sep=" ")
#SeuObj@meta.data$bayesspace<-bayes$spatial.cluster
#SeuObj@meta.data$bayesspace <- as.factor(SeuObj@meta.data$bayesspace)
#saveRDS(SeuObj, "_bayesspace.rds")
