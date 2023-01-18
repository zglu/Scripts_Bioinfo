#Rscript [cellType.R] [SCT rds] [feature] [name]

args<-commandArgs(T)


library(Seurat)
library(ggplot2)

SeuObj<-readRDS(args[1])
goi<-as.character(args[2])

DefaultAssay(SeuObj) <- "RNA"
#This is a bit of a misnomer in the plot: you're plotting the embedding values for PC1, not expression of a feature. These values are allowed to be negative (look at a DimPlot with reduciton = 'pca') so it's no concern that the violin plot show negative values.

pdf(paste0(args[2], "_", args[3],"_vln.pdf"), width=8, height=4)

# violin plot using ggplot2

vln_df<-data.frame(GOI =SeuObj[["RNA"]]@data[ goi, ], cluster = SeuObj$seurat_clusters)

#                     GOI cluster
#AAACCCACAATACCTG-1_1   0       1
#AAACCCACACCTCAGG-1_1   0       9

noise <- rnorm(n = length(x = vln_df[, "GOI"])) / 100000
vln_df$GOI<-vln_df$GOI+noise
ggplot(vln_df, aes(x=cluster, y=GOI, fill=cluster)) + geom_violin(adjust=1,trim=TRUE,scale="width") + geom_jitter()+theme_classic()+labs(title=goi,x="Cluster", y="Expression Level")+NoLegend()

dev.off()

