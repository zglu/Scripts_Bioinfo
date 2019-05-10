#source("https://bioconductor.org/biocLite.R")
#biocLite("Sushi")

library(Sushi)
## standard plotting
bdata<-read.delim("~/Desktop/sushi.bed", sep=" ", header = T)
chrom            = "chr15"
chromstart       = 72998000
chromend         = 73020000

plotGenes(bdata, chrom, chromstart,chromend, maxrows=1,height=0.5,plotgenetype="arrow",bentline=FALSE,col="blue", labeloffset=.4,fontsize=1.2,arrowlength = 0.025)
labelgenome(chrom, chromstart,chromend,side=1,scipen=20,n=3,scale="Mb",line=.18,chromline=.5,scaleline=0.5)


## test data: plot one gene
bdata<-read.delim("~/Desktop/batch-sushi/sub.bed", sep=" ", header = T)
chrom = "SM_V7_2" # although defined, sushi will plot genes in the range in all chromosomes
chromstart = 11022560
chromend = 11202431

plotGenes(bdata, chrom=chrom,chromstart = chromstart,chromend = chromend, types=bdata$type, maxrows=50,bheight=0.08,plotgenetype="box",bentline=FALSE,col="brown", labeloffset=.2,fontsize=0.9,arrowlength = 0.025,labeltext=TRUE) # "arrow"
labelgenome(chrom, chromstart,chromend,side=1,scipen=10,n=5,scale="Mb",line=.18,chromline=.3,scaleline=0.8)


### or sebsetting to a separate file
# cat rattPac-all.gff | grep CDS | sed 's/;/ /'| awk '{print $1 " " $4 " " $5 " " $10 " " "." " " $7 " " $3}'|sed 's/Parent=//' | sed 's/+/1/g'>ratt.bed 
# can also plot utr
# pay attention to same genes on different contigs
bed<-read.delim("sub.bed", sep=" ", header=T) # can use whole genome annotation from different sources: combine them into one file
# cat rattPac-all.gff | grep gene | awk '{print $9 " " $1 " " $4-600 " " $5+600}'|sed 's/ID=//'|sort > ratt-ids
grange<-read.delim("sub-Smps.txt", sep=" ", header=F) # with start end and chr
pdf("batch-sushi.pdf",width = 8.3, height = 11.7);
par(mfrow=c(4,2));
i<-1;
for (i in 1:nrow(grange)) {
chr = as.character(grange[i, "V2"])
chrstart = grange[i, "V3"]
chrend = grange[i, "V4"]
bed.sub <- bed[which(bed[,"chrom"] == chr & bed[,"start"] >= chrstart & bed[,"end"] <= chrend),]
write.table(bed.sub, file = "subbed", quote = F, row.names = F)
subbed<-read.delim("subbed", sep=" ", header=T)
plotGenes(subbed, chr, chrstart, chrend, types=subbed$type, maxrows=50,bheight=0.08,plotgenetype="box",bentline=FALSE,col="brown", labeloffset=.2,fontsize=0.9,arrowlength = 0.025,labeltext=TRUE)
# plot using bed.sub causes an error: only defined on a data frame with all numeric variables
labelgenome(chr, chrstart,chrend,side=1,scipen=10,n=3,scale="Mb",line=.18,chromline=.8,scaleline=0.8)
}
dev.off()

## plotBed
bed1<-read.delim("~/Documents/SCHISTO/V7/beds/grunau-test.bed", sep=" ", header=T)
bed2<-read.table("~/Documents/SCHISTO/V7/beds/isoseq-test.bed", sep=" ", header=F)
bed3<-read.delim("~/Documents/SCHISTO/V7/beds/gff-test.bed", sep=" ", header=T)
chrom = "SM_V7_1" # although defined, sushi will plot genes in the range in all chromosomes
chromstart = 1185000
chromend = 1255000

par(mfrow=c(3,1));
plotBed(bed1, chrom, chromstart, chromend, row  = "auto",wiggle=0.01)
plotBed(bed2, chrom, chromstart, chromend, row  = "auto",wiggle=0.01)# row = 1, wiggle = 0.001
plotGenes(bed3, chrom=chrom,chromstart = chromstart,chromend = chromend, types=bdata$type, maxrows=50,bheight=0.08,plotgenetype="box",bentline=FALSE,col="brown", labeloffset=.2,fontsize=0.9,arrowlength = 0.025,labeltext=TRUE)
labelgenome(chrom,chromstart,chromend,n=2,scale="Mb")

## bam bed files
bambed<-read.delim("/Users/zl3/Documents/SCHISTO/V7/bams/mergedIsoseq.bed", sep=" ", header=T)
chr = "SM_V7_3"
chrstart = 35968584
chrend = 35977794

plotBed(bambed, chr, chrstart, chrend, row  = "auto",wiggle=0.01)
labelgenome(chr,chrstart,chrend,n=2,scale="Kb")
