library(ggplot2)
library(ggrepel)
clstall<-read.delim("combined_clusters.txt", sep=" ", header = F, col.names=c("chr", "block", "func", "name", "genes", "fdr", "start", "end", "chrlen", "source", "label"))
#eg SM_V7_1 0 PF00169 PH 3 0.0450167 2212203 4703829 88881357
clst<-clstall[which(clstall$genes>=3),]

#pdf(file=paste0(args[1], "_geneCounts.pdf"), width=11.3, height=7.7)
# facet_grid to 
# ggrepel: https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
# geom_bar plot gene bars, add text and geom_rect to plot chromosome frame

pp<-ggplot(data=clst, aes(x=start/1000000, y=genes, fill=source)) + geom_bar(stat="identity", position="dodge", width=0.3) + geom_text_repel(data=clst[which(clst$label=="yes"),], mapping=aes(color=source, label=paste0(name, " (", genes, ")")), size=3.3, direction="both", max.iter=6000)  + geom_rect(mapping=aes(xmin=0, xmax=chrlen/1000000, ymin=0, ymax=25), fill="white", color="black", alpha=0, size=0.3) + facet_grid(rows=vars(chr)) 

pp<-pp+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + labs(x="Coordinate (Mb)", y="Genes") + ylim(0, 25) + scale_y_continuous(breaks=c(0, 10, 20)) + scale_fill_manual(labels=c("Gene family", "Pfam"),values=c("blue", "red")) + scale_color_manual(labels=c("Gene family", "Pfam"),values=c("blue", "red")) + theme(legend.position="bottom") + theme(legend.text = element_text(size=9.5, face="bold")) # + theme_classic() / theme_bw()

pp

#dev.off()
