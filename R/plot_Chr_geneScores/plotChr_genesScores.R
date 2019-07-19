# Plot DEGs on chromosomes
library(ggplot2)
library(ggrepel)
args = commandArgs(trailingOnly=TRUE)
mydata<-read.delim(args[1], sep=" ", header = F, col.names=c("chr", "start", "end", "name", "score","strand", "chrlen"))
#eg SM_V7_1 11326094 11402473 Smp_124210  5.192      + 88881357

maxgene<-max(abs(mydata$score))

pdf(file=paste0(args[1], "_plotScore.pdf"), width=11.3, height=7.7)
# for each DEG plot a bar according to the abs value of logFC (score in the input), 
# and color it according to up- (>0 red) or down-regulation (<0 blue).
# take your chose which genes to label. I used score>=6 here
pp<-ggplot(data=mydata, aes(x=start/1000000, y=abs(mydata$score), fill=score>0)) + geom_bar(stat="identity", position="identity", width=0.6)+ scale_fill_manual(name="Score",values=c("blue", "red"), labels=c("< 0", "> 0")) + geom_text_repel(data=mydata[which(mydata$score>=6),], mapping=aes(x=start/1000000, y=score,label=paste0(name, " (", score, ")")), size=2, color="black")  + geom_rect(mapping=aes(xmin=0, xmax=chrlen/1000000, ymin=0, ymax=maxgene), fill="white", color="red", alpha=0.01, size=0.5) + facet_grid(rows=vars(chr))  
pp<-pp+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + labs(x="Coordinate (Mb)", y="Score") + theme_classic() # + theme(legend.position="none"))
pp

dev.off()
