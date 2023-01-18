# QC plots after binning or cellbin merging
# Rscript a1.2_qc_after_binning.R [rds with raw matrix] [nFeature cutoff] [nCount cutoff]
# output: pdf with QC plots 


args<-commandArgs(T)

if (length(args)<3) {
  stop("Usage: Rscript a1.3_qc_filter_after_binning.R [bined matrix] [nFeature cutoff] [nCount cutoff]")
}

suppressMessages(library(Seurat))
library(ggplot2)
library(ggsci)
library(gridExtra)
library(Matrix)

SeuObj0<-readRDS(args[1])

message("Original matrix:")
SeuObj0

message("Summary of nGene and nMID:")
summary(SeuObj0@meta.data$nFeature_RNA)
summary(SeuObj0@meta.data$nCount_RNA)

nFeatureCut<-as.numeric(args[2])
nCountCut<-as.numeric(args[3])

SeuObj<-subset(SeuObj0, subset = nFeature_RNA > nFeatureCut & nCount_RNA > nCountCut)

message(">>>")
message(paste0("Filtering conditions: nFeature > ", nFeatureCut, " and nCount > ", nCountCut))

message(">>>")
message("Filtered matrix:")
SeuObj

message("Summary of nGene and nMID:")
summary(SeuObj@meta.data$nFeature_RNA)
summary(SeuObj@meta.data$nCount_RNA)


saveRDS(SeuObj, file=paste0(args[1],"_nGene", args[2], "_nCount", args[3], ".RDS"))

pdf(paste0(args[1], "_nGene", args[2], "_nCount", args[3], "_QC.pdf"), width=8, height=4)
#p1<-FeatureScatter(SeuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1.2<-ggplot(SeuObj@meta.data,aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
#  	stat_smooth(method=lm) +
#  	scale_x_log10() + 
#  	scale_y_log10() + 
#  	geom_vline(xintercept = 500) +
#  	geom_hline(yintercept = 250) +
  	theme_classic() 

# complexity
SeuObj$log10GenesPerUMI <- log10(SeuObj$nFeature_RNA) / log10(SeuObj$nCount_RNA)
p2<-ggplot(SeuObj@meta.data, aes(x=log10GenesPerUMI)) +
  	geom_density(alpha = 0.4, fill="khaki1") + labs(title = "log10 Genes per UMI")+
  	theme_classic() #+
  	#geom_vline(xintercept = 0.8)

#p3<-FeatureScatter(SeuObj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
p3.2<-ggplot(SeuObj@meta.data,aes(x=nFeature_RNA, y=percent.mt, color=nCount_RNA)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
#  	scale_x_log10() + 
#  	scale_y_log10() + 
#  	geom_vline(xintercept = 250) +
#  	geom_hline(yintercept = 2.5) +
  	theme_classic()

p4<-VlnPlot(SeuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1)

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
      theme_classic()

p8 <- ggplot(SeuObj@meta.data, aes(x=nCount_RNA)) + 
      geom_density(alpha = 0.4, fill="bisque")+labs(title = "distribution of nUMI")+
      theme_classic()

p9 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nFeature_RNA))+
  geom_tile()+
  theme_void()+
  coord_fixed()+
  scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+
  labs(fill="nGene",x=NULL,y=NULL,title = NULL)+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+ 
  theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))

p10 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= nCount_RNA))+
    geom_tile()+
    theme_void()+
    coord_fixed()+
    scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+
    labs(fill="nUMI",x=NULL,y=NULL,title = NULL)+
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15) )+    
    theme(plot.title = element_text(hjust =0.5,vjust = 0.5,size = 15))

grid.arrange(p9,p10,ncol = 2, top="Spatial distribution of nGene and nUMI")
grid.arrange(p6,p5,ncol = 2, top="nGene and nUMI per cell (bin unit)")
grid.arrange(p7,p8,p2,ncol = 3)
#print(p4)
grid.arrange(p1.2,p3.2,ncol = 2)

C = SeuObj@assays$RNA@counts
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[10:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.05, las = 1, main="Top expressed genes", xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[10:1], horizontal = TRUE)


dev.off()

message("QC plots exported.")

