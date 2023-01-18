#Rscript checkFeatures.R [sct rds] [id]
# for i in `cat markers.ids`; do Rscript ../checkFeatures_single.R SS200000745BL_A2_merged_lasso_bin1.gem_BIN20.rds_nGene15_nCount30.RDS_SCT.rds $i; done


args<-commandArgs(T)

if (length(args)<2) {
  stop("Usage: Rscript checkFeatures.R [rds] [id]")
}

suppressMessages(library(Seurat))
library(ggplot2)
suppressMessages(library(dplyr))
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

suppressMessages(library(plotly))

f <- list(
  family = "Arial",
  size = 14,
  color = 'blue'
)

p<-ggplot(SeuObj@meta.data,aes(coord_x, coord_y, fill= tmexp[,1], text=paste0("data: ",tmexp[,1])))+ 
  geom_tile()+
  theme_minimal()+
  coord_fixed()+ 
  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
  labs(title=goi, fill="sct data") #Counts sct-normalised data

figplotly1<-ggplotly(p, tooltip=c("coord_x", "coord_y", "text")) %>% layout(font=f)

htmlwidgets::saveWidget(as_widget(figplotly1), paste0(args[2],"_spatialExp.html"), title=paste0(args[2], " Spatial Expression in ", args[1]))
