# cell type annotation using reference with RCTD, CelliD, and SingleR
# by default is cell-level annotation
# Rscript Seurat_refAnno.R [scRNAseq_rds] [spatial_rds] [doublet/full/multi][name prefix]
# https://github.com/dmcable/spacexr/tree/master/vignettes

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(RColorBrewer)
  library(gridExtra)
  library(ggplot2)
  library(CelliD)
  library(SingleR)
})
  
args<-commandArgs(T)

if (length(args)<3) {
  stop("Usage: Rscript RCTD.R [scRNAseq_rds] [spatial_rds] [doublet/full/multi] [name perfix]")
}


myCol<-unique(c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"),brewer.pal(12, "Paired"), 
                brewer.pal(8, "Accent"), brewer.pal(8, "Dark2")))

# single-cell data
sc.data<-readRDS(args[1])
sp.data<-readRDS(args[2])


DefaultAssay(sc.data)<-"RNA"

message("Use seurat_clusters as cell type names")

####### SingleR on clusters or single cells ############ http://bioconductor.org/books/devel/SingleRBook/
message("SingleR based on similarities (Spearman correlation) in expression profiles")

# SingleR() expects reference datasets to be normalized and log-transformed (for marker detection)
sc.data<-NormalizeData(sc.data)
sp.data<-NormalizeData(sp.data, assay = "RNA")
# identifies marker genes from the reference and uses them to compute assignment scores (based on the Spearman correlation across markers) for each cell in the test dataset against each label in the reference. 
# The label with the highest score is the assigned to the test cell, possibly with further fine-tuning to resolve closely related labels.
sc.labels<-sc.data@meta.data$seurat_clusters
sc.expr<-as.matrix(sc.data@assays$RNA@data)
sp.expr<-as.matrix(sp.data@assays$RNA@data)

# de.methods: classic (sorts genes by logFC and take top de.n; ) / "wilcox" / 
# de.n: top n DE genes; classic: de.n=500 * (2/3) ^ log2(N), N=unique labels; otherwise default=10
singleR.classic<-SingleR(test = sp.expr, ref = sc.expr, labels = sc.labels) #clusters=seurat_clusters for cluster-level annotation
singleR.wilcox<-SingleR(test=sp.expr, ref=sc.expr, labels = sc.labels, de.method = "wilcox", de.n = 20) #doesn't work； aggr.ref=TRUE

#add labels to metadata
sp.data@meta.data$SingleR_classic<-singleR.classic$pruned.labels[match(rownames(sp.data@meta.data), rownames(singleR.classic))]
table(sp.data$SingleR_classic)
sp.data@meta.data$SingleR_wilcox<-singleR.wilcox$pruned.labels[match(rownames(sp.data@meta.data), rownames(singleR.wilcox))]
table(sp.data$SingleR_wilcox)

### SingleR can also take custom markers
# or use custom markers; note that by importing the markers.csv any space and special character were changed to .
roc_markers<-SCINA::preprocess.signatures("scRNA-seq/top20_roc_markers.csv")
# change marker names back
names(roc_markers)<-gsub("\\.", " ", names(roc_markers)) # . back to space
names(roc_markers)<-gsub("   ", " & ", names(roc_markers))

singleR.roc<-SingleR(test=sp.expr, ref = sc.expr, labels = sc.labels, genes = roc_markers)
sp.data$SingleR_roc<-singleR.roc$pruned.labels[match(rownames(sp.data@meta.data), rownames(singleR.roc))]

####### CelliD ###########
message("CelliD based on hypergeomtric test")

#sc data should also be log-normed; run Multiple Correspondence Analysis (MCA)
sc.data<-ScaleData(sc.data, features=rownames(sc.data))
sc.data<-RunMCA(sc.data)
# use group signatures; eg. cluster signatures
sc.groupsig<-GetGroupGeneSet(sc.data, dims=1:50, n.features = 200, group.by = "seurat_clusters")

# or use cell signatures
#sc.cellsig<-GetGroupGeneSet(sc.data, dims=1:50, n.features = 200)

# run matching from cluster signatures to query cells
sp.data<-NormalizeData(sp.data)
sp.data<-ScaleData(sp.data, features=rownames(sp.data))
sp.data<-RunMCA(sp.data)

hgt.pergrp<-RunCellHGT(sp.data, pathways = sc.groupsig, dims = 1:50, n.features = 200)
pergrp.pred <- rownames(hgt.pergrp)[apply(hgt.pergrp, 2, which.max)]
pergrp.pred.signif <- ifelse(apply(hgt.pergrp, 2, max)>2, yes = pergrp.pred, "unassigned")

# Save cell type predictions as metadata within the Seurat object
sp.data$CelliD_sc <- pergrp.pred.signif
table(sp.data$CelliD_sc)



####### RCTD/spacexr ########## 
# uses raw counts; Counts should be untransformed count-level data.
message("RCTD/spacexr")

sc.counts<-sc.data@assays$RNA@counts
meta_data<-sc.data@meta.data
# use labels from seurat_clusters
cell_types<-meta_data$seurat_clusters; names(cell_types)<-rownames(meta_data)
nUMI<-meta_data$nCount_RNA; names(nUMI)<-rownames(meta_data)

reference<-Reference(sc.counts, cell_types, nUMI)

# spatial data
sp.counts<-sp.data@assays$RNA@counts
coords<-sp.data@meta.data[, c("coord_x", "coord_y")]
sp.nUMI<-colSums(sp.counts)
puck<-SpatialRNA(coords, sp.counts, sp.nUMI)

hist(log(puck@nUMI,2))

barcodes <- colnames(puck@counts)
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 

# create RCTD 
message(paste0("RCTD mode: ", args[3]))

myRCTD<-create.RCTD(puck, reference, max_cores = 1) # for parallel processing, the number of cores used. 1 means no parallel processing.
myRCTD <- run.RCTD(myRCTD, doublet_mode = args[3]) 
# ‘doublet’ (at most 1-2 cell types per pixel), 
# ‘full’ (no restrictions on number of cell types), 
# ‘multi’ (finitely many cell types per pixel, e.g. 3 or 4)

#saveRDS(myRCTD, file = paste0(args[4],"_RCTD_", args[3], ".rds"))

results<-myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA

#resultsdir <- paste0(args[4], '_RCTD_Plots')
#if(!dir.exists(resultsdir))
#  dir.create(resultsdir)

#RCTD results are stored in the @results field. Of particular interest is @results$weights, 
#a data frame of cell type weights for each pixel (for full mode [i.e. doublet_mode = ‘full’]). 
#Alternatively This section will generate various plots which can be found in resultsdir. 
#The results of ‘doublet_mode=“doublet”’ are stored in @results$results_df and @results$weights_doublet, 
#the weights of each cell type. More specifically, the results_df object contains one column per pixel (barcodes as rownames).
#
#`spot_class`, a factor variable representing RCTD’s classification in doublet mode: 
#  “singlet” (1 cell type on pixel), “doublet_certain” (2 cell types on pixel), 
#  “doublet_uncertain” (2 cell types on pixel, but only confident of 1), 
#  “reject” (no prediction given for pixel).
#Next, the first_type column gives the first cell type predicted on the bead (for all spot_class conditions except “reject”).
#The second_type column gives the second cell type predicted on the bead for doublet spot_class conditions (not a confident prediction for “doublet_uncertain”).

if (args[3]=='full'){
#  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
#  plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 
  results@all_weights[1:3, ]
} else if (args[3]=='doublet'){
#  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
  doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 
#  plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names)
  table(results$results_df$spot_class)
  results$weights[1:3, ]
} else {
  results$sub_weights[1:3,]
}

## combine max weight with meta data
norm_weights = normalize_weights(results$weights); norm_weights<-as.data.frame(norm_weights)
noReject<-rownames(subset(results$results_df, subset=spot_class!="reject"))
norm_weights<-norm_weights[which(rownames(norm_weights) %in% noReject),]

maxcol<-colnames(norm_weights)[max.col(norm_weights, ties.method = "first")]
maxweight<-cbind(rownames(norm_weights), maxcol); 
maxweight<-as.data.frame(maxweight)

sp.data@meta.data$RCTD_maxweight<-maxweight$maxcol[match(rownames(sp.data@meta.data), maxweight$V1)]
sp.data@meta.data$RCTD_spotclass<-results$results_df$spot_class[match(rownames(sp.data@meta.data), rownames(results$results_df))]

saveRDS(sp.data, file=paste0(args[2], "_refAnno.rds"))

pdf(paste0(args[2], "_singleR_CelliD_RCTD-maxweight.pdf"), width=15, height=15)
# SingleR plots
#pHeat<-plotScoreHeatmap(singleR.classic)
#pDelta<-plotDeltaDistribution(singleR.classic, ncol=3)
#grid.arrange(pHeat, pDelta, ncol=2)
ggplot(sp.data@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SingleR_classic))+geom_tile()+
  theme_void()+ coord_fixed()+scale_fill_manual(values = myCol)
ggplot(sp.data@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), color= SingleR_classic))+geom_point(size=0.05)+
  theme_void()+ coord_fixed()+scale_color_manual(values = myCol)

ggplot(sp.data@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SingleR_wilcox))+geom_tile()+
  theme_void()+ coord_fixed()+scale_fill_manual(values = myCol)
ggplot(sp.data@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), color= SingleR_wilcox))+geom_point(size=0.05)+
  theme_void()+ coord_fixed()+scale_color_manual(values = myCol)

## SingleR diagnostics
plotScoreHeatmap(SingleR_classic)
plotDeltaDistribution(SingleR_classic)


# CelliD perdictions
DimPlot(sp.data, reduction="spatial", group.by = "CelliD_sc", cols = myCol, pt.size = 0.001)+coord_fixed()+theme_void()



# RCTD max weight plots
ggplot(sp.data@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= RCTD_maxweight))+geom_tile()+
  theme_void()+ coord_fixed()+scale_fill_manual(values = myCol)
ggplot(sp.data@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), color= RCTD_maxweight))+geom_point(size=0.05)+
  theme_void()+ coord_fixed()+scale_color_manual(values = myCol)
dev.off()
