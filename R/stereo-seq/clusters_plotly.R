#Rscript clusters_plotly_bin.R [sct RDS file]

args<-commandArgs(T)

if (length(args)<3) {
  cat("Usage: Rscript clusters_plotly.R [rds] [bin|cellbin]\n")
  q()
}

library(ggsci)
myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
      pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
			pal_jama("default")(7),pal_jco("default")(10),
			pal_locuszoom("default")(7),pal_startrek("uniform")(7),
			pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
			pal_simpsons("springfield")(16),
			pal_gsea("default")(12)))

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

SeuObj<-readRDS(args[1])

suppressMessages(library(plotly))

f <- list(
  family = "Arial")

if (args[2]=="bin") {
  p3 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= seurat_clusters))+ 
        geom_tile()+theme_void()+coord_fixed()+ scale_fill_manual(values = myCol)+labs(fill="Seurat clusters")
} else {
  p3<-ggplot(SeuObj@meta.data, aes(coord_x, coord_y))+geom_point(aes(color = seurat_clusters), shape=".")+theme_void()+coord_fixed()+scale_color_manual(values = myCol)
}

figplotly1<-ggplotly(p3) %>% layout(font=f)

# save as HtmlWigdet
htmlwidgets::saveWidget(as_widget(figplotly1), paste0(args[1],"_SeuratClusters.html"), title=paste0(args[1], "_SeuratClusters"))
