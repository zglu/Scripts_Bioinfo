

args<-commandArgs(T)

if (length(args)<2) {
  cat("Usage: Rscript distribution_plotly.R [Seurat RDS] [bin|cellbin]\n")
  q()
}

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

SeuObj<-readRDS(args[1])

suppressMessages(library(plotly))


#figplotly1<-plot_ly(data=SeuObj@meta.data, x=~coord_x, y=~coord_y, text = ~paste('nGene', nFeature_RNA, "\n", 'nUMI', nCount_RNA), color = ~nFeature_RNA)
#figplotly2<-plot_ly(data=SeuObj@meta.data, x=~coord_x, y=~coord_y, text = ~paste('nUMI', nCount_RNA, '\n', 'nGene', nFeature_RNA), color = ~nCount_RNA)

if (args[2]=="cellbin") {
  p1<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y, text=paste('nCount_RNA:', nCount_RNA)))+geom_point(aes(color = nFeature_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nGene") #+scale_color_viridis_c(option="turbo")
  p2<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y, text=paste('nFeature_RNA:', nFeature_RNA)))+geom_point(aes(color = nCount_RNA), shape=".")+theme_void()+coord_fixed()+scale_color_gradient(low = "#fcf5eb", high = "#800080")+labs(fill="nUMI") #+scale_color_viridis_c(option="turbo")
} else {
  p1<-ggplot(SeuObj@meta.data,aes(coord_x, coord_y, fill= nFeature_RNA, text=paste('nCount_RNA:', nCount_RNA)))+ geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+labs(fill="nGene",x="coord_x",y="coord_y",title = NULL) #as.numeric(coord_x)
  p2<-ggplot(SeuObj@meta.data,aes(coord_x, coord_y, fill= nCount_RNA, text=paste('nFeature_RNA:', nFeature_RNA)))+ geom_tile()+ theme_void()+ coord_fixed()+ scale_fill_gradient(low = "palegoldenrod", high = "magenta4")+labs(fill="nUMI",x="coord_x",y="coord_y",title = NULL)
}


figplotly1<-ggplotly(p1)
figplotly2<-ggplotly(p2)



# save as HtmlWigdet
htmlwidgets::saveWidget(as_widget(figplotly1), paste0(args[1],"_nGene_spatial.html"), title=paste0(args[1], "_nGene"))
htmlwidgets::saveWidget(as_widget(figplotly2), paste0(args[1],"_nCount_spatial.html"), title=paste0(args[1], "_nMID"))
