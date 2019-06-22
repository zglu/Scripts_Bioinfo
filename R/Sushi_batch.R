#source("https://bioconductor.org/biocLite.R")
#biocLite("Sushi")

# usage: Rscript sushi_single_gene.R [file with gene id]
# will need: the gene coord bed file, and features (cd/utr/exon) bed file

library(Sushi)
### or sebsetting to a separate file
genebed<-read.delim("Sm_v7.3.bed", sep="\t", header=F, col.names=c("chr", "start", "end", "id", "score", "strand"))
combibed<-read.delim("Combined.bed", sep="\t", header=F, col.names=c("chrom", "start", "end", "gene", "score", "strand", "type", "col"))

args = commandArgs(trailingOnly=TRUE)
ids<-read.table(args[1], header=F)
myGenes<-as.character(ids$V1)


myGenes$match<-match(myGenes, genebed$id, nomatch=0)
y<-genebed[myGenes$match,]

pdf(paste0(args[1], ".pdf"), width=8.3, height=11.7)
par(mfrow=c(3,1))
i<-1;
for (i in 1:nrow(y)) {
chr=as.character(y[i,"chr"])
regstart<-ifelse(y[i, "strand"]=="+", y[i, "start"]-500, y[i, "start"]-200)
regend<-ifelse(y[i, "strand"]=="-", y[i, "end"]+500, y[i, "end"]+200)

bed.sub<-combibed[which(combibed[,"chrom"] == chr & combibed[,"start"] >= regstart & combibed[,"end"] <= regend),]
write.table(bed.sub, file="subbed", row.names = F, quote = F)
subbed<-read.delim("subbed", sep=" ", header = T)

plotGenes(subbed, chr, regstart, regend, types=subbed$type, maxrows=50,bheight=0.12,plotgenetype="box",bentline=FALSE,col=as.character(subbed$col), labeloffset=.3,fontsize=0.8,arrowlength = 0.025,labeltext=TRUE)

labelgenome(chr, regstart,regend,side=1,scipen=10,n=3,scale="Mb",line=.18,chromline=.8,scaleline=0.8)
}
dev.off()
