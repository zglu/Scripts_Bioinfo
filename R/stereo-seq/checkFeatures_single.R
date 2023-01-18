#Rscript checkFeatures.R [sct rds] [id]
# for i in `cat markers.ids`; do Rscript ../checkFeatures_single.R SS200000745BL_A2_merged_lasso_bin1.gem_BIN20.rds_nGene15_nCount30.RDS_SCT.rds $i; done


args<-commandArgs(T)

if (length(args)<2) {
  stop("Usage: Rscript checkFeatures.R [rds] [id]")
}

suppressMessages(library(Seurat))
library(ggplot2)
#suppressMessages(library(gridExtra))

SeuObj<-readRDS(args[1])

##### top markers plots

goi<-as.character(args[2])
## raw MID counts
#tmexp<-as.matrix(SeuObj@assays$RNA@counts[goi,]) #@counts - unnormalised counts; @data normlised data (log or sct)  

## SCT normalised data
tmexp<-as.matrix(SeuObj@assays$SCT@data[goi,]) #@counts - unnormalised counts; @data normlised data (log or sct)  

## SCT unnormalised data
#tmexp<-as.matrix(SeuObj@assays$SCT@counts[goi,]) #@counts - unnormalised counts; @data normlised data (log or sct)  

pdf(paste0(args[2],"_sctData.pdf"))

p<-ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= tmexp[,1]))+ 
  geom_tile()+
  theme_void()+
  coord_fixed()+ 
  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
  labs(title=goi, fill="") #Counts sct-normalised data

print(p)
dev.off()
