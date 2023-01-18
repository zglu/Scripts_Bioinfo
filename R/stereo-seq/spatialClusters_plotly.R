# Rscript spatialClusters_plotly.R [SCT rds file]

args<-commandArgs(T)

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggsci))
suppressMessages(library(plotly))

myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
      pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
			pal_jama("default")(7),pal_jco("default")(10),
			pal_locuszoom("default")(7),pal_startrek("uniform")(7),
			pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
			pal_simpsons("springfield")(16),
			pal_gsea("default")(12)))


SeuObj<-readRDS(args[1])

# change res
SeuObj <- FindClusters(SeuObj, res=0.5, verbose = FALSE)

p1 <- ggplot(SeuObj@meta.data,aes(coord_x, coord_y, fill= seurat_clusters))+ 
  geom_tile()+
  theme_void()+
  coord_fixed()+ 
  scale_fill_manual(values = myCol)+
  labs(fill="Seurat clusters")



figplotly1<-ggplotly(p1)

# save as HtmlWigdet
htmlwidgets::saveWidget(as_widget(figplotly1), paste0(args[1],"_Seurat_spatial.html"), title=paste0(args[1], "_Seurat"))
