#Rscript checkFeatures.R [sct rds] [id list]


args<-commandArgs(T)

if (length(args)<2) {
  stop("Usage: Rscript checkFeatures.R [rds] [id list]")
}

suppressMessages(library(Seurat))
library(ggplot2)
suppressMessages(library(gridExtra))

SeuObj<-readRDS(args[1])

##### top markers plots

idlist<-read.table(args[2], header=F)
myFeatures<-as.vector(idlist$V1)

data<-SeuObj@assays$SCT@data
subdata<-data[which(rownames(data) %in% idlist$V1),]

p <- lapply(myFeatures, function(x) {
tmexp<-as.matrix(subdata[x,]) #@counts - unnormalised counts; @data normlised data (log or sct)  
ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= tmexp[,1]))+ 
  geom_tile()+
  theme_void()+
  coord_fixed()+ 
  scale_fill_gradient(low = "#fcf5eb", high = "#800080")+
  labs(title=x, fill="") # sct-normalised data
})

ggsave(
   filename = paste0(args[1],"_checkFEATURES.pdf"), 
   plot = marrangeGrob(p, nrow=3, ncol=3), 
   width = 15, height = 9
)

message("Feature spatial plots exported.")
