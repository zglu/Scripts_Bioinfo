siggo<-read.delim("~/Documents/Inqueries/Carmen_sc_functions/sigGOterms.txt", sep="\t", header=T)
# go_name	-log10fdr	cell_type	color
# motor activity	2.27572413	Actin2+ muscle	#1B5E20
# actin cytoskeleton	1.698970004	Actin2+ muscle	#1B5E20
# supramolecular fiber	3.346787486	Positional muscle	#00c853
# motor activity	2.142667504	Positional muscle	#00c853
siggo001<-siggo[which(siggo$X.log10fdr>2),]
# reverse the row order
siggo001<- siggo001[seq(dim(siggo001)[1],1),]
siggo001$color<-as.character(siggo001$color)
# add space to the left
par(mar=c(5,18,4,1)+.5)
# plot with  axis label for y
# useful: names.arg, horiz, las, xpd, xex.names
barplot(siggo001$X.log10fdr, names.arg=siggo001$go_name, col=siggo001$color, horiz = T, las=2, xlim=c(0,10), xpd=F, border="black", cex.names=0.78, xlab=expression(paste("-log"[10], "FDR")), axes=F)
# add x axis and label, paralell to axis
axis(1, at=seq(0,10,by=5), las=0, labels=c("0", "5", "10"))
# add legend, define the shap, size, and reverse the order
legend("topright", legend=unique(rev(siggo001$cell_type)), col=unique(rev(siggo001$color)), pch=15, bty="n", cex=0.85, pt.cex=1.3)

#library(ggplot2)
#siggo001$go_name<-factor(siggo001$go_name, levels=unique(siggo001$go_name))
#ggplot(siggo001, aes(factor(go_name), X.log10fdr, fill=color))+geom_bar(stat="identity", position="dodge")+coord_flip()+theme_bw()
