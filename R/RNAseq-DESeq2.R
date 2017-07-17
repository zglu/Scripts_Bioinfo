rm(list=ls()) #remove all objects
ls()
library("DESeq2")
#install.packages( pkgs= "gplots" )
library(gplots)
dir()
sessionInfo()
#Read data: counts and gene length for all genes
rd <-read.delim("Counts_ordered_ed.txt", sep = " ", row.names="GeneID")
head(rd)

############## define Condition, Counts and Coldata; DGE
dr.cond<- as.factor(substr(colnames(rd[,1:23]), 1,2))

dr.counts<-rd[,1:23]
dr.counts<-as.matrix(dr.counts)
head(dr.counts)

dr.coldata<-data.frame(row.names=colnames(rd[,1:23]), dr.cond)
dr.coldata

dds<-DESeqDataSetFromMatrix(countData=dr.counts, colData=dr.coldata, design=~dr.cond)
dds
dds<-DESeq(dds)

###get separate results###
results(dds, contrast=c("dr.cond","bO", "sO"))
## turn outlier cutoff by set off cooksCutoff
ress<-results(dds, cooksCutoff=FALSE, contrast=c("dr.cond","bO","sO")); 
ress<-ress[order(ress$padj),]; 
write.table(ress, file="ovDGE2.txt")

# check data
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plot(res$baseMean+1, -log10(res$pvalue), cex=.4)
##############
res<-results(dds)
res
mcols(res, use.names=TRUE)
summary(res)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "GeneID"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")
################

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

##PCA plot get percentage info
plotPCA(rld, intgroup=c("dr.cond"), ntop=9224)  # + returnData=TRUE return values for self-plotting; , col=colopca doesn't work

pca_data<-plotPCA(rld, intgroup=c("dr.cond"), ntop=9224, returnData = TRUE)
## ggplot for pca plotting
#p<-ggplot(pca_data, aes(x=PC1, y=PC2, fill=group))
#p+geom_point(size=3)+scale_fill_manual(values=c("red2","blue4","darkgoldenrod1","deepskyblue4","lightcoral","blue","wheat2","deepskyblue1"))

plot(pca_data$PC1, pca_data$PC2, col=coloa, xlab="PC1", ylab="PC2", pch=19, cex.lab=1.2)
legend('right',col=colo3,legend=c("bM","sM","bT","sT","bF","sF","bO","sO"), pch=19)

## on mac col works for plotPCA
colopca<-c("red2","blue4","darkgoldenrod1","deepskyblue4","lightcoral","blue","wheat2","deepskyblue1")
colo3<-c("blue4","blue","deepskyblue4","deepskyblue1","red2","lightcoral","darkgoldenrod1","wheat2")
coloa<-c(rep("blue4",3), rep("blue",3), rep("deepskyblue4",3),rep("deepskyblue1",3),
         rep("red2",3),rep("lightcoral",3),rep("darkgoldenrod1",3),rep("wheat2",2))


library(RColorBrewer)
hmcol2 <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(dr.cond))])
#levels: bF bM bO bT sF sM sO sT
macolo<-c("red2","blue4","darkgoldenrod1","deepskyblue4","lightcoral","blue","wheat2","deepskyblue1")

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=macolo[dr.cond], RowSideColors=macolo[dr.cond],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# PCA plot
#plotPCA(rld)

## heatmap of count matrix
heatmap.2(counts(dds,normalized=TRUE)[1:100,], col = hmcol2,
          Rowv = FALSE, Colv = FALSE, scale="row", labRow=NA,
          dendrogram="none", trace="none", margin=c(10,6))

## normalized read counts
rc.norm<-counts(dds, normalized=TRUE)
head(rc.norm)
rc.bMm<- rowMeans(rc.norm[, 1:3]); write.table(rc.bMm, file="bM.txt")
## calc row sd
library(matrixStats)
rc.bMsd<-rowSds(rc.norm[,1:3])


### top 35 genes
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           +            trace="none", dendrogram="column", 
           +            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           +            ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
             +                colData(rld)$treatment ] )

