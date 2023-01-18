library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

args<-commandArgs(T)


myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"),
     brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))


SeuObj<-readRDS(args[1])

pdf(paste0(args[1], "_SCT.rds_CLUSTERS.pdf"),width=20,height=10)

p1.2<-ggplot(SeuObj@meta.data,aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	theme_classic() 

# complexity
SeuObj$log10GenesPerUMI <- log10(SeuObj$nFeature_RNA) / log10(SeuObj$nCount_RNA)
p2<-ggplot(SeuObj@meta.data, aes(x=log10GenesPerUMI)) +
  	geom_density(alpha = 0.4, fill="khaki1") + labs(title = "log10 Genes per UMI")+
  	theme_classic() #+

p3.2<-ggplot(SeuObj@meta.data,aes(x=nFeature_RNA, y=percent.mt, color=nCount_RNA)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	theme_classic() 


median <- median(SeuObj@meta.data$nCount_RNA)
p5 <- ggplot(SeuObj@meta.data,aes(orig.ident,nCount_RNA))+
    geom_violin() + 
    labs(y="nUMI",x=NULL) +
    theme(panel.grid.major=element_line(colour=NA))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 15))+   
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),plot.title = element_text(hjust = 1,vjust=-40,size=14) )+
    theme(axis.text.x = element_text(size = 15,color="black"),axis.text.y = element_text(size = 15,color="black")) +
    geom_boxplot(width=0.2)+ annotate("text",label=as.character(median),y=as.numeric(median),x=1.1,hjust = 0,color='red',size=5)

median <- median(SeuObj@meta.data$nFeature_RNA)
p6 <- ggplot(SeuObj@meta.data,aes(orig.ident,nFeature_RNA))+
    geom_violin() + 
    labs(y="nGene",x=NULL) +
    theme(panel.grid.major=element_line(colour=NA))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 15))+   
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),plot.title = element_text(hjust = 1,vjust=-40,size=14) )+
    theme(axis.text.x = element_text(size = 15,color="black"),axis.text.y = element_text(size = 15,color="black")) +
    geom_boxplot(width=0.2)+ annotate("text",label=as.character(median),y=as.numeric(median),x=1.1,hjust = 0,color='red',size=5)


p7 <- ggplot(SeuObj@meta.data, aes(x=nFeature_RNA)) + 
      geom_density(alpha = 0.4, fill="lightblue")+labs(title = "distribution of nGene")+ 
      theme_classic()#+

p8 <- ggplot(SeuObj@meta.data, aes(x=nCount_RNA)) + 
      geom_density(alpha = 0.4, fill="bisque")+labs(title = "distribution of nUMI")+
      theme_classic()#+

p9<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = nFeature_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nGene") #+scale_color_viridis_c(option="turbo")
p10<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = nCount_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nUMI") #+scale_color_viridis_c(option="turbo")

pp3<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = seurat_clusters), shape=".")+theme_void()+coord_fixed()+scale_color_manual(values = myCol)

p0.2 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.2))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.2")

p0.4 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.4))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.4")

p0.6 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.0.6))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res0.6")

p1.0 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= SCT_snn_res.1))+
  geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat res1")

pp4 <- DimPlot(SeuObj, reduction = "umap", cols=myCol, label = TRUE) + NoLegend() + coord_fixed()
pp7<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols=myCol)
pp8<-FeatureScatter(SeuObj, group.by = "ident",feature1 = "nFeature_RNA", feature2 = "percent.mt", cols=myCol)

grid.arrange(p9,p10,ncol = 2, top="Spatial distribution of nGene and nUMI")
grid.arrange(p6,p5,ncol = 3, top="nGene and nUMI per cell")
grid.arrange(p7,p8,p2,ncol = 3)
grid.arrange(p1.2,p3.2,ncol = 2)
grid.arrange(pp4,pp3,ncol = 2, top="Seurat clustering")
grid.arrange(pp7,pp8,ncol = 2, top="Scatter plot grouped by clusters")
grid.arrange(p0.2,p0.4,ncol = 2, top="Seurat other resolutions")
grid.arrange(p0.6,p1.0,ncol = 2)
dev.off()

