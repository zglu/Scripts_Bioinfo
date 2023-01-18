#Rscript distribution_plotly.R [cellbin RDS file]

args<-commandArgs(T)

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

SeuObj<-readRDS(args[1])

suppressMessages(library(plotly))


p1<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y, text=paste('nCount_RNA:', nCount_RNA)))+geom_point(aes(color = nFeature_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nGene") #+scale_color_viridis_c(option="turbo")
p2<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y, text=paste('nFeature_RNA:', nFeature_RNA)))+geom_point(aes(color = nCount_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nUMI") #+scale_color_viridis_c(option="turbo")


figplotly1<-ggplotly(p1)
figplotly2<-ggplotly(p2)

# save as HtmlWigdet
htmlwidgets::saveWidget(as_widget(figplotly1), paste0(args[1],"_nGene_spatial.html"), title=paste0(args[1], "_nGene"))
htmlwidgets::saveWidget(as_widget(figplotly2), paste0(args[1],"_nCount_spatial.html"), title=paste0(args[1], "_nCount"))
