#source("https://bioconductor.org/biocLite.R")
#biocLite("Sushi")

# usage: Rscript sushi_single_gene.R [gene id]
# will need: the gene coord bed file, and features (cd/utr/exon) bed file

library(Sushi)
genebed<-read.delim("Sm_v7.3.bed", sep="\t", header=F, col.names=c("chr", "start", "end", "id", "score", "strand"))
combibed<-read.delim("Combined.bed", sep="\t", header=F, col.names=c("chrom", "start", "end", "gene", "score", "strand", "type", "col"))


args = commandArgs(trailingOnly=TRUE)
gene<-args[1]
shift<-as.numeric(args[2])

gene$match<-match(gene, genebed$id, nomatch=0)
y<-genebed[gene$match,]

chr=as.character(y$chr)
regstart<-ifelse(y$strand=="+", y$start-shift, y$start)
regend<-ifelse(y$strand=="-", y$end+shift, y$end)

bed.sub<-combibed[which(combibed[,"chrom"] == chr & combibed[,"start"] >= regstart & combibed[,"end"] <= regend),]
write.table(bed.sub, file="subbed", row.names = F, quote = F)
subbed<-read.delim("subbed", sep=" ", header = T)

pdf(paste0(args[1],"_sushi.pdf"), width=6, height=9)
par(mfrow=c(2,1))
plotGenes(subbed, chr, regstart, regend, types=subbed$type, maxrows=20,bheight=0.15,plotgenetype="box",bentline=FALSE,col=as.character(subbed$col), labeloffset=.3,fontsize=0.6,arrowlength = 0.025,labeltext=TRUE)
#plotBed(isobed, chr, regstart, regend, row="auto", wiggle=0.01)

labelgenome(chr, regstart,regend,side=1,scipen=10,n=3,scale="Mb",line=.18,chromline=.8,scaleline=0.8)
dev.off()
