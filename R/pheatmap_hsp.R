library(pheatmap)
library(RColorBrewer)

# Rscript pheatmap_hsp.R [Male/Female]

# in this script the 2nd column in the expression file is the group names
args = commandArgs(trailingOnly=TRUE)
exp<-read.delim(paste0("hsp_expression_", args[1], ".txt"), sep=" ",header=T)

z<-exp[,-1:-2]
rownames(z)<-exp[,1]
zz<-as.matrix(z)

anno<-data.frame(row.names=exp$id, Group=exp$group)

#newCols <- colorRampPalette(grDevices::rainbow(length(unique(anno$group))))
#annoCol <- newCols(length(unique(anno$group)))
#names(annoCol) <- unique(anno$group)
#annoCol <- list(category = annoCol)

# male
ma.annoCol<-list(Group=c(Chaperonin_CCT="orange", DNAj_Hsp40="red", Hsp70="blue", Hsp90="darkgrey", `Co-chaperones`="cyan", Ubiquination="green", Stress_granules="purple"))

# female
fe.annoCol<-list(Group=c(Chaperonin_CCT="orange", DNAj_Hsp40="red", Hsp60_Chaperonin="green", Hsp70="blue", Hsp90="darkgrey", `Co-chaperones`="cyan"))

if (args[1] == "Male") {
	annoCol <- ma.annoCol
} else {
	annoCol <- fe.annoCol
}

#pheatmap(zz, scale = "row", cluster_cols = F, cluster_rows = F, annotation_row = anno)

# set defined cell width and height
pheatmap(zz, scale = "row", cluster_cols = F, cluster_rows = F, annotation_row = anno, annotation_colors = annoCol, cellheight = 10, cellwidth = 10, main=paste0(args[1], " HSPs"), file = paste0(args[1], "_hsp_expr_ZT.png"))
