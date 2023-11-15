# compare two columns in the meta data using sankey plot

args<-commandArgs(T)

if (length(args)<3) {
  cat("Usage: Rscript meta_sankey.R [Seurat rds] [metadata col1] [metadata col2]\n")
  q()
}

seu_obj<-args[1]
col1<-args[2]
col2<-args[3]

library(Seurat)
library(highcharter)
library(ggplot2)
library(ggrepel)
library(ggalluvial)

SeuObj<-readRDS(seu_obj)
#SeuObj$RCTD_maxweight[which(is.na(SeuObj$RCTD_maxweight))]<-"unassigned"
meta.sel<-SeuObj@meta.data[, c(col1, col2)]
sankeyform<-data_to_sankey(meta.sel)

p<-ggplot(as.data.frame(sankeyform), aes(y=weight, axis1=from, axis2=to))+geom_flow(aes(fill=from, na.rm=FALSE), width=1/3)+
  geom_stratum(width=1/3)+geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=3, color="blue")+theme_void()+guides(fill="none")+
  scale_x_discrete(limits = c(col1, col2))

pdf(paste0(seu_obj, "_", col1, "-", col2, ".pdf"), width=6, height=8)
print(p)
dev.off()
