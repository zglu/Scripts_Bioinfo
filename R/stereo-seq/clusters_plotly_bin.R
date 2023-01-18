#Rscript clusters_plotly_bin.R [sct RDS file]

args<-commandArgs(T)

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

p3 <- ggplot(SeuObj@meta.data,aes(as.numeric(coord_x), as.numeric(coord_y), fill= seurat_clusters))+ 
  geom_tile()+
  theme_void()+
  coord_fixed()+ 
  scale_fill_manual(values = myCol)+
  labs(fill="Seurat clusters")


figplotly1<-ggplotly(p3) %>% layout(font=f)

# save as HtmlWigdet
htmlwidgets::saveWidget(as_widget(figplotly1), paste0(args[1],"_SeuratClusters.html"), title=paste0(args[1], "_SeuratClusters"))
