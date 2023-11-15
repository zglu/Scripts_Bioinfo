# do DE on specific group
# check which marker set do the clusters found as DEG?

args<-commandArgs(T)

if (length(args)<4) {
  cat("Usage: Rscript groupDE_overlap.R [normalised_rds] [group name] [wilcox/roc/t] [sig table]\n")
  q()
}


suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
})

myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"),
                brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))

spData<-args[1]
grp<-as.character(args[2]) #eg. SCT_snn_res.0.4 or other metadata colnames
testMethod<-as.character(args[3])
sigtable<-args[4]

#Aatk         Mesenchymal-like
#Ackr1            Blood_vessels

SeuObj<-readRDS(spData)

SeuObj<-NormalizeData(SeuObj)

DefaultAssay(SeuObj) <- "RNA"
SeuObj<-SetIdent(SeuObj, value=SeuObj@meta.data[[grp]])

deg <- FindAllMarkers(SeuObj,only.pos = TRUE, test.use = testMethod)#, min.pct = 0.25, logfc.threshold = 0.5) # method: wilcox
if (testMethod=="wilcox" || testMethod=="t") {
  deg<-subset(deg, deg$p_val<1e-10)  
}

sigMarkers<-read.delim(sigtable, header=F, sep="\t")
deg_comb<-merge(deg, sigMarkers, all.x=TRUE, by.x="gene", by.y="V1")
deg_comb<-deg_comb[order(deg_comb$cluster, -deg_comb$avg_log2FC),]
if (testMethod=="roc") {
  colnames(deg_comb)[9]<-"cell_type_signature"
  deg_comb<-deg_comb[!is.na(deg_comb$cell_type_signature),]
  poverlap<-ggplot(deg_comb)+geom_jitter(aes(x=cluster, y=myAUC, color=cell_type_signature), position=position_jitter(0.2))+theme_classic()+labs(title=grp)
} else {
  colnames(deg_comb)[8]<-"cell_type_signature"
  deg_comb<-deg_comb[!is.na(deg_comb$cell_type_signature),]
  poverlap<-ggplot(deg_comb)+geom_jitter(aes(x=cluster, y=-log10(p_val_adj), color=cell_type_signature), position=position_jitter(0.2))+theme_classic()+labs(title=grp)
}

library(openxlsx)
write.xlsx(deg_comb,paste0(spData, "_allMarkers_", grp, "-", testMethod, ".xlsx"))

pdf(paste0(spData, "_allMarkers_", grp, "-", testMethod, ".pdf"), width=6, height=4)
print(poverlap)
dev.off()


