## USAGE: Rscript plotExpression.R [gene ids] 

## expression data
args = commandArgs(trailingOnly=TRUE)

genExp <- read.delim("/Users/zlu/zl3/Kate-Circadian/_Kate_cycling_genes_Feb_18_Zhigang/Adam_cyclingGenes_201912/exprData_202005/plot-both-gender/meanExpr-cycling-mixed.txt", sep=" ")

# a list of gene ids as input
ids<-read.table(args[1], header=F)

ids$match<-match(ids$V1, genExp$id, nomatch=0)
y<-genExp[ids$match,]
z<-y[,-1:-2]
rownames(z)<-y[,1]
#colnames(z)<-as.character(seq(0,44,by=4))

## calculate relative expression values (value/group mean)
for (i in 1:nrow(z)) {
	idmean<-rowMeans(z[i,]);
	for (j in 1:ncol(z)) {
		z[i,j] = z[i,j]/idmean # avoid using rowMeans(z[i,]) here as the value changes after assigning each z[i,j]

	}
} 

z$id<-y[,1]

z$gender<-y[,2]

colnames(z)<-c(seq(0,44,by=4), "id", "gender")

## convert multiple columns into rows
## https://uc-r.github.io/tidyr

library(tidyr)
long.z <- z %>% gather(Time, Exp, 1:12) # ZT0:ZT44

library(ggplot2)
#ggplot(data = long.z, aes(x=Time, y=Exp, group=id)) + geom_point() + geom_line(aes(linetype=id, color=id))
pdf(paste0(args[1],".pdf"),width=3,height=3)

# normal line
p<-ggplot(data = long.z, aes(x=as.numeric(Time), y=Exp, group=id)) + geom_line(aes(color=gender)) + scale_color_manual(values=c("Female"="red", "Male"="blue"))

# smooth line
#p<-ggplot(data = long.z, aes(x=as.numeric(Time), y=Exp, group=id)) + geom_point(aes(color=id))
#p<-p + geom_smooth(method="loess", aes(color=id))

# label theme etc
p<-p + xlab("Circadian Time") + ylab("Relative Expression") + ggtitle(args[1]) + theme_bw() + theme(legend.position="none") + scale_x_discrete(limits=c(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44), labels=c("8","12","16","20","0","4","8","12","16","20","0","4")) + theme(axis.title.y = element_text(size = rel(0.8))) + theme(axis.title.x = element_text(size = rel(0.8)))
print(p)
dev.off()
