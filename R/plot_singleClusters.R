#Rscript plot_singleClusters.R [SCT rds]

args<-commandArgs(T)

library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggsci)

myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
                pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
                pal_jama("default")(7),pal_jco("default")(10),
                pal_locuszoom("default")(7),pal_startrek("uniform")(7),
                pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
                pal_simpsons("springfield")(16),
                pal_gsea("default")(12)))


SeuObj<-readRDS(args[1])

nClusters<-as.character(0:(max(as.numeric(SeuObj@meta.data$seurat_clusters))-1))
ps<- lapply(nClusters, function(x) {
  clusterdata<-subset(SeuObj@meta.data, subset=seurat_clusters==x)
  ggplot(clusterdata,aes(as.numeric(coord_x), as.numeric(coord_y), fill= seurat_clusters))+ 
    geom_tile()+
    theme_void()+
    coord_fixed()+ 
    labs(title=paste0("Cluster ", x), fill="")
})

ggsave(
  filename = paste0(args[1],"_singleClusters.pdf"), 
  plot = marrangeGrob(ps, nrow=3, ncol=3), 
  width = 15, height = 9
)
