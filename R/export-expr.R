rm(list=ls())
library("edgeR")
library("DESeq2")
# read counts table
reads <- read.delim("~/Documents/v7-meta/SmMeta1.txt", sep=" ", row.names="id") # last column is the gene length
ncol <- ncol(reads)
groups <- as.factor(substr(colnames(reads[,1:(ncol-1)]), 1,3))
counts <- reads[, 1:(ncol-1)]

#################### edgeR #########################
# build dgelist object
dgelist <- DGEList(counts = counts, group = groups)
glength <- reads[, ncol(reads)]

# library size
barplot(dgelist$samples$lib.size*1e-6, names=colnames(dgelist$counts),
        ylab="Library size (millions)")

# MDS plot
plotMDS(dgelist, top=2000, main = "500 / lLFC") #top=2000
plotMDS(dgelist, method="bcv", main="500 / BCV") 

# normalisation
dgelist<-calcNormFactors(dgelist)

# export RPKM after normalisation
rpkm_norm<-rpkm(dgelist, glength, normalized.lib.sizes = TRUE, log=FALSE)
barplot(rpkm_norm["Smp_070360",])
write.table(rpkm_norm, file="rpkm_norm_edgeR.txt", quote=F, col.names=T, row.names=T)
cat('RPKM values written to "rpkm_norm.txt"')

#################### DESeq 2 #########################
counts_matrix <- as.matrix(counts)
coldata<-data.frame(row.names=colnames(counts), groups)
dds<-DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~groups)
dds
dds<-DESeq(dds)

#PCA
rld <- rlogTransformation(dds)
plotPCA(rld, intgroup=c("groups"), ntop=10742)

## diff_exp export with comparisons
dge1<-results(dds, cooksCutoff=FALSE, contrast=c("groups","egg","bOv"))
dge1<-dge1[order(dge1$padj),]; 
#write.table(dge1, file="ovDGE2.txt", quote=F)

# check data
res<-results(dds)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plot(res$baseMean+1, -log10(res$pvalue), cex=.4)

## export counts
## normalized read counts
rc_norm<-counts(dds, normalized=TRUE)
barplot(rc_norm["Smp_070360",])
write.table(nc_norm, file="counts_norm_deseq2.txt", quote=F, col.names=T, row.names=T)

## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "GeneID"
head(resdata)
## Write results
#write.csv(resdata, file="deseq2_counts_norm.csv")

