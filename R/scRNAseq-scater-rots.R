## Use Scater for single cell analysis (normalisation)
library(scater, quietly = TRUE)
rd_raw<-read.delim("~/Documents/scRNAseq/1_RT_plate/reads_all_trim_1.tab", sep="\t", row.names="id")
rd<-rd_raw[,7:24]
inf<-read.delim("~/Documents/scRNAseq/1_RT_plate/cells_info_1.txt")
pd<-new("AnnotatedDataFrame", data = inf)
rownames(pd) <- pd$Cell
rd_sceset <- newSCESet(countData = rd, phenoData = pd)
keep_feature <- rowSums(counts(rd_sceset) > 0) > 0  #exprs()
rd_sceset <- rd_sceset[keep_feature,]
rd_sceset <- calculateQCMetrics(rd_sceset, feature_controls = 1:50)
scater_gui(rd_sceset) # will open browser
counts(rd_sceset)[1:10, 1:6]
exprs(rd_sceset)[1:3, 1:6]

#get_exprs(rd_sceset, exprs_values = "counts")

###### plotting ######
## QC plot
plotQC(rd_sceset, type = "exprs-freq-vs-mean")

# Expression distribution
plot(rd_sceset, block1 = "Number", block2 = "Type",
     colour_by = "Type", nfeatures = 300, exprs_values = "exprs")

# expression plot  #show_median = TRUE, show_violin = FALSE,
plotExpression(rd_sceset, rownames(rd)[1:6],
               x = "Type", exprs_values = "counts", colour = "Type")

## PCA
plotPCA(rd_sceset, colour_by = "Type", theme_size = 15) # ntop = 1000,

# detect outliers
rd_sceset <- plotPCA(rd_sceset, pca_data_input = "pdata", 
                     detect_outliers = TRUE, return_SCESet = TRUE)

## plot QC
plotQC(rd_sceset, type = "find-pcs", variable = "Type",
       plot_type = "pairs-pcs")
###### end plotting ######

## Normalisation
rd_sceset <- normaliseExprs(rd_sceset, exprs_values = "counts")
#rd_sceset <- normalise(rd_sceset)

### now accessible
# norm_counts norm_exprs norm_cpm norm_tpm norm_cpm
norm_counts(rd_sceset["Smp_095350",])

rd_normexprs<-data.frame(norm_exprs(rd_sceset))
plot(rd_normexprs[,"S1_1A"], rd_normexprs[,"S1_1B"])

# expression plot
plotExpression(rd_sceset, "Smp_095350", show_median = TRUE, show_violin = FALSE, #1:6,
               x = "Type", exprs_values = "norm_exprs", colour = "Type")

plotHighestExprs(rd_sceset, col_by_variable = "total_features", n = 50, exprs_values = "counts")

## Use ROTS for differential gene expression analysis
### should take normalised data as input
library("ROTS")
input <- rd_normexprs[, c(2,3,9,10)] # Oi vs Om
# groups <- as.factor(substr(colnames(rd), 1,2)) 
# groups = c(rep("Oi",2), rep("Om",2))
results <- ROTS(data=input, groups = groups, log = TRUE, B=100, K=500, seed=NULL)
#names(results)
summary(results, fdr = 0.05)
results$logfc["Smp_008900"]
write.table(results$logfc, file = "rots-results.txt")
write.table(results$FDR, file = "rots-results-fdr.txt")
plot(results, fdr = 0.05, type = "volcano")
plot(results, fdr = 0.05, type = "heatmap")

## Volcano plot
ovol<-read.delim("~/Documents/scRNAseq/1_RT_plate/rots_Om-Oi_logfc_fdr.txt", sep = "\t", row.names = 'geneid')
ovol$threshold=as.factor(abs(ovol$log2FC)>0.585 & ovol$FDR<0.05)
p1<-ggplot(data=ovol, aes(x=log2fc, y=-log10(fdr), colour=ovol$threshold))+geom_point(alpha=0.4, size=1.5)+
  +     xlab(expression(paste("log"[2],"FC")))+ylab(expression(paste("-log"[10],"FDR")))+ggtitle("Om/Oi")+
  +     scale_colour_manual(values=c("TRUE"="red","FALSE"="black"))+theme(legend.position = "none")
plot(p1)