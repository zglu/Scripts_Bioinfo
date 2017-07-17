rm(list=ls()) #remove all objects
ls()
library("edgeR")
dir()
sessionInfo()
#Read data
rd <-read.delim("Counts_ordered_ed.txt", sep=" ", row.names="GeneID")
head(rd)
cor.test(rd[,"sM1"], rd[,"sF1"]) #pearson
cor(rd[,1:6],method="spearman")
#Pearson correlation is a measure of linearity while 
#Spearman correlation is for measuring monotonicity
plot(rd$sM1, rd$sF1); abline(lm(rd$sM1 ~ rd$sF1), col="blue") # pch=20

rd.grp<- as.factor(substr(colnames(rd[,1:24]), 1,2))
rdlist<-DGEList(counts=rd[,1:24], group=rd.grp); 

##### BASIC EXPLORATION ####
# check MDS of all samples
cola2<-c(rep("blue4",3), rep("blue",3), rep("deepskyblue4",3),rep("deepskyblue1",3),
         rep("red2",3),rep("lightcoral",3),rep("darkgoldenrod1",3),rep("wheat2",3))
## MDS plot
plotMDS(rdlist, col=cola2, main = "500 / lLFC") #top=2000
plotMDS(rdlist,col=cola2, method="bcv", main="500 / BCV") 

###lib size accross all samples
barplot(rdlist$samples$lib.size*1e-6, names=colnames(rdlist$counts), col=cola2, 
        ylab="Library size (millions)")

## overview of sample relationship: sub-group
#pairs(log2(rdlist$counts[, 1:6]), pch='.', lower.panel=NULL)
pairs(log2(rdlist$counts[, c("sM1","sM2","sM3","sF1","sF2","sF3")]), 
      pch = '.', col = "blue", lower.panel=NULL) ##change
#or
plot(rd[1:10, 1:4], pch='.', lower.panel=NULL)

#################################### 
#### get RPK (reads per kilo bases)
rpk<-(rd[,1:23]/rd$length)*1000  #### rpk>3.333 as one transcript
write.table(rpk, file="~/Copy/seq/RPK-raw.dat")
#### get raw rpkm table for filtering
glength<-rd[,25]  #glength<-rd$length
rdlist.rpkm<-rpkm(rdlist, glength, normalized.lib.sizes=FALSE, log=FALSE)
# rpkm= [(reads/gene_length)*1000]/library_size *1000000
# normalization factors not calculated thus not difference
write.table(rdlist.rpkm, file="rpkm_Raw.txt")
## IN SHELL FILTER RPKM2 in at least one group and make a new read table

######### set sub-group (without sO3) ########
rd2 <-read.delim("readcount+length_rpkm2-filtered.table", row.names="gene_ID")
dr.grp<- as.factor(substr(colnames(rd2[,1:23]), 1,2)) #substr first letters of each column
dr<-DGEList(counts=rd2[,1:23], group = dr.grp) # genes = rawdata[,1])
dr.length=rd2[,25] ###
dim(dr)

## MDS using dots
par(mfrow=c(1,2))
plotMDS(dr, top=500,labels= '.', cex=5, col=cola2, main = "lLFC")
legend('right',col=colo3,legend=c("bM","sM","bT","sT","bF","sF","bO","sO"), pch=15)

plotMDS(dr, top=9224,labels= '.', cex=5, col=cola2, method="bcv", main = "BCV")
legend('right',col=colo3,legend=c("bM","sM","bT","sT","bF","sF","bO","sO"), pch=15)

#summary(dr.rpkm>2)
hist(cpm(dr)[,"bO1"], plot=FALSE) #breaks=c(0,2,3,5,10,20000)

### filtration is not needed
dr$samples

dr<-calcNormFactors(dr)
dr$samples

## CPM and rpkm makes the same
cpm<-cpm(dr, normalized.lib.sizes = TRUE)

dr.rpkm<-rpkm(dr, dr.length, normalized.lib.sizes = TRUE, log=FALSE) ####
##### dot plot #####
library(ggplot2)
sample=dr.grp;sample<-factor(sample, levels(sample)[c(2,6,4,8,1,5,3,7)])
dd<-data.frame(RPKM=rpkmx["Smp_000020.1",], Sample=sample)
g<-ggplot(dd, aes(Sample, RPKM))+geom_point(colour=colo2,size=3)+geom_errorbar(stat="hline", yintercept="mean", width=0.6, aes(ymax=..y.., ymin=..y..))
g+guides(colour=FALSE)+labs(title="Smp_000020", x="", y="Reads per kilobase per million")  +theme_bw()+theme(axis.text=element_text(size=12))

##### barplots ######
mean.bM <- rowMeans(dr.rpkm[, c(1, 2, 3)]); write.table(mean.bM, file="bM.txt")
library(matrixStats)
sd.bM<-rowSds(dr.rpkm[,1:3]); write.table(sd.bM, file="sd-bM.txt", quote=F, col.names=F, row.names=F)
sem.bM<-sd.bM/sqrt(3) ###20170505 SEM(x) = SD(x)/sqrt(length(x)) here bM has 3 records
## ... mean rpkm was calculated
mRPKM<-read.delim("~/Copy/seq/rpkmMean_2filter_allNorm1", sep = " ", header = T, row.names="gene_ID")
mr<-as.matrix(mRPKM)
mRPKM2<-read.delim("~/Copy/seq/rpkmMean2_prod", sep="\t",header = T, row.names="GeneID")
mrp<-as.matrix(mRPKM2)

colo3<-c("blue4","blue","deepskyblue4","deepskyblue1","red2","lightcoral","darkgoldenrod1","wheat2")
barplot(mr["Smp_179600.1", ], col = colo3, las = 2) 

idx<-"Smp_123300.1"
barplot(mr[idx, ], col = colo3, space = 0, 
        main = idx, 
        ylab="Average expression (RPKM)")
# relative
barplot(mr[idx, ]/max(mr[idx,]), col = colo3, space = 0, 
        main = paste(idx,"[SmFst]"), 
        ylab="Relative expression")
#title(idx, line=3) # adjust title position
## return scaled matrix
scaled<-t(apply(mr, 1, function(x)x/max(x)))

## a file per plot
GI<-scan(file="~/Copy/iloop//014-3-3.txt", what=character())
i<-1;
for (i in 1:length(GI)) {
  pdf(paste(GI[i],".pdf", sep=""));
  barplot(mr[GI[i],], col = colo3, space = 0, main=GI[i], ylab="Mean of normalized RPKM");
  mtext(mrp[GI[i], 9],side=1,line=3,at=c(0,0),adj=0, col="purple", font=3)
  dev.off()
}
## several plots in one page ******
GI<-scan(file="~/Copy/iloop//appendix.genes", what=character());
i<-1; pdf('appendix.pdf', width = 8.3, height = 11.7); #A4 size
par(mfrow=c(4,2)); ## can be removed so one plot per page (3,3)
for (i in 1:length(GI)){
  mb<-barplot(mr[GI[i], ], col = colo3, space=0, main=GI[i], xpd=F, 
          ylab="Average expression (RPKM)");
  mtext(mrp[GI[i], 9],side=1,line=3,at=c(0,0),adj=0, col="purple", font=3)};
dev.off()
###### add error bars    ######
sdRPKM<-read.delim("~/Dropbox/port/rpkmSD_norm2.txt", header=T, sep=" ", row.names = "GeneID")
msd<-as.matrix(sdRPKM)
GI<-scan(file="~/Copy/iloop//appendix.genes", what=character());
i<-1; pdf('appendix.pdf', width = 8.3, height = 11.7); #A4 size
par(mfrow=c(4,2)); ## can be removed so one plot per page (3,3)
for (i in 1:length(GI)){
  mb<-barplot(mr[GI[i], ], col = colo3, space=0, main=GI[i], xpd=F, 
              ylab="Average expression (RPKM)");
  segments(mb,mr[GI[i], ]-msd[GI[i], ], mb, mr[GI[i], ]+msd[GI[i], ], lwd = 1.5)
  arrows(mb,mr[GI[i], ]-msd[GI[i], ], mb, mr[GI[i], ]+msd[GI[i], ], lwd = 1.5, angle = 90, code = 3, length = 0.05)  
  mtext(mrp[GI[i], 9],side=1,line=3,at=c(0,0),adj=0, col="purple", font=3)};
dev.off()

#colo<-as.numeric(dr$samples$group) #### USEFUL ###
colo2<-c(rep("blue4",3), rep("blue",3), rep("deepskyblue4",3),rep("deepskyblue1",3),
         rep("red2",3),rep("lightcoral",3),rep("yellow",3),rep("wheat2",2))

## MDS plot
par(mfrow=c(1,2))
plotMDS(dr, col=colo2, main = "500 / lLFC") #top=2000
plotMDS(dr,col=colo2, method="bcv", main="500 / BCV") 

###lib size accross all samples before and after filtering
par(mfrow=c(2,1))
barplot(rdlist$samples$lib.size*1e-6, names=colnames(rdlist$counts), col=cola2,
        ylab="Library size (millions) original")
barplot(dr$samples$lib.size*1e-6, names=colnames(dr$counts), col=colo2,
        ylab="Library size (millions) after filtration")

#### 1.1 ANOVA-like** ### for housekeeper
design<-model.matrix(~dr.grp) ##first as intercept = (~1+grp)
rownames(design)<-colnames(dr); design
dr<-estimateGLMCommonDisp(dr, design, verbose=TRUE)
dr<-estimateGLMTrendedDisp(dr, design)
dr<-estimateGLMTagwiseDisp(dr, design)
plotBCV(dr)

fit<-glmFit(dr, design)
colnames(design)
lrta<-glmLRT(fit, coef=2:8) # all cmp; null hypothesis is all group means being equal
topTags(lrta, n=20)
write.table(topTags(lrta, n=Inf), file="glm_anova.txt")
lrta["Smp_079600.1"]
summary(dt<-decideTestsDGE(lrta))  ##ERROR because of multi-comparisons

# 1.2 specific group constrast(s): construct manually
design<-model.matrix(~0+dr.grp); colnames(design)
rownames(design)<-colnames(dr); design
dr<-estimateGLMCommonDisp(dr, design, verbose=TRUE)
dr<-estimateGLMTrendedDisp(dr, design)
dr<-estimateGLMTagwiseDisp(dr, design)

fit<-glmFit(dr, design)
con<-makeContrasts(Ov=dr.grpbO-dr.grpsO, Te=dr.grpbT-dr.grpsT, 
                   Fe=dr.grpbF-dr.grpsF, Ma=dr.grpbM-dr.grpsM, levels=design); con
lrt.Ov <- glmLRT(fit, con[,"Ov"])
topTags(lrt.Ov, n=10, sort.by=c("logFC"))

## bT vs sT; bO vs sO
te.con <- makeContrasts(dr.grpbT-dr.grpsT, levels=design)
ov.con<-makeContrasts(dr.grpbO-dr.grpsO, levels=design)
ma.con<-makeContrasts(dr.grpbM-dr.grpsM, levels=design)
fe.con<-makeContrasts(dr.grpbF-dr.grpsF, levels=design)

lrto<-glmLRT(fit, contrast=ov.con)
lrtt<-glmLRT(fit, contrast=te.con)
lrtm<-glmLRT(fit, contrast=ma.con)
lrtf<-glmLRT(fit, contrast=fe.con)

topTags(lrto, n=10, sort.by=c("logFC"))
topTags(lrtt, n=10, sort.by=c("logFC"))
lrto["Smp_079600.1"]
summary(dt<-decideTestsDGE(lrtt, p=0.05, adjust = "BH"))
detags<-rownames(dr)[as.logical(dt)]
plotSmear(lrtt, de.tags=detags, ylab="logFC: bT/sT")  ## adjust p above
abline(h=c(-0.585,0.585),col="blue", lty=3)

## Worms vs Organs
OvW.con<-makeContrasts((dr.grpbT+dr.grpsT+dr.grpbO+dr.grpsO)/4
                       -(dr.grpbM+dr.grpsM+dr.grpbF+dr.grpsF)/4, levels = design)
lrt.ow<-glmLRT(fit, contrast = OvW.con)
topTags(lrt.ow)
#options(scipen=6) #show only 6 dicimals
write.table(topTags(lrt.ow, n=Inf), file = "glm_ORGANvsWORM.txt")

## Males-sF cluster
MasF<-makeContrasts((dr.grpbM+dr.grpsM+dr.grpsF)/3 - (dr.grpbT+dr.grpsT+dr.grpbF+dr.grpbO+dr.grpsO)/5, 
                     levels=design)
lrt.MasF<-glmLRT(fit, contrast=MasF)
topTags(lrt.MasF)

### sF-males VS bF-males
masf2<-makeContrasts((dr.grpbM+dr.grpsM)/2 - dr.grpsF, levels=design)
lrt.masf2<-glmLRT(fit, contrast=masf2)

mabf2<-makeContrasts((dr.grpbM+dr.grpsM)/2 - dr.grpbF, levels=design)
lrt.mabf2<-glmLRT(fit, contrast=mabf2)

## bO-bF, bT-bM contracts
bOF<-makeContrasts(dr.grpbO-dr.grpbF, levels=design)
bTM<-makeContrasts(dr.grpbT-dr.grpbM, levels=design)
lrt.bOF<-glmLRT(fit, contrast=bOF)
lrt.bTM<-glmLRT(fit, contrast=bTM)

## gender effect
gender.con<-makeContrasts((dr.grpbT+dr.grpsT+dr.grpbM+dr.grpsM)/4
                       -(dr.grpbO+dr.grpsO+dr.grpbF+dr.grpsF)/4, levels = design)
lrt.gender<-glmLRT(fit, contrast=gender.con)
write.table(topTags(lrt.gender, n=Inf), file="MaTe-FeOv_contrast.txt")

bMF.con<-makeContrasts(dr.grpbM-dr.grpbF, levels=design)
sMF.con<-makeContrasts(dr.grpsM-dr.grpsF, levels=design)
bTO.con<-makeContrasts(dr.grpbT-dr.grpbO, levels=design)
sTO.con<-makeContrasts(dr.grpsT-dr.grpsO, levels=design)

lrt.bMF<-glmLRT(fit, contrast=bMF.con)
lrt.sMF<-glmLRT(fit, contrast=sMF.con)
lrt.bTO<-glmLRT(fit, contrast=bTO.con)
lrt.sTO<-glmLRT(fit, contrast=sTO.con)

#### 2. GLM more groups #####
Pairing<-factor(substring(colnames(dr), 1,1)) #b-sex/s-sex
Gender<-factor(c(rep("M", 12), rep("F", 11))) ## without sO3 Male/Female
#t<-c("W","W","W","W","W","W","O","O","O","O","O","O")
#Tissue<-factor(c(rep(t,2)))
Tissue<-factor(c("W","W","W","W","W","W","O","O","O","O","O","O",
                 "W","W","W","W","W","W","O","O","O","O","O")) #Worm/Organ
#design<-model.matrix(~Pairing+Gender+Tissue+Pairing:Gender+Pairing:Tissue+Gender:Tissue+Pairing:Gender:Tissue)
#design<-model.matrix(~Pairing * Gender * Tissue)  #equals last
design<-model.matrix(~Pairing + Gender + Tissue)  ## 
#data.frame(Pairing, Gender, Tissue)
rownames(design)<-colnames(dr); design

### Continue
dr<-estimateGLMCommonDisp(dr, design, verbose=TRUE)
dr<-estimateGLMTrendedDisp(dr, design)
dr<-estimateGLMTagwiseDisp(dr, design)
plotBCV(dr)
#sqrt(d$common.dispersion)  

plotMeanVar(dr, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)

### see separate pairs
par(mfrow=c(1,2))
plotSmear(dr, pair=c("sT","bT"), ylim=c(-10,10))
abline(h=c(-0.585, 0.585), col = "blue")
plotSmear(dr, pair=c("sO","bO"), ylim=c(-10,10))
abline(h=c(-0.585, 0.585), col = "blue")
###

fit<-glmFit(dr, design)
colnames(design)

lrt<-glmLRT(fit, coef=2:4) #select according to colnames(mm)
#lrti<-glmLRT(fit, coef=8) # genes affected by all 3 factors(interaction)
lrtp<-glmLRT(fit, coef=2) # Pairing
lrtg<-glmLRT(fit, coef=3) # Gender
lrtt<-glmLRT(fit, coef=4) # Tissue
#FDR <- p.adjust(lrt$table$PValue, method="BH")
topTags(lrt, n=20, sort.by=c("logFC")) #FDR
lrt["Smp_179600.1"]  ## check single gene
write.table(topTags(lrtp, n=Inf), file= "glm_Pairings.1.txt", sep="\t")

## optional: see cpm of top 10 
colnames(design)
o<-order(lrt$table$FDR)
rpkm(dr, dr.length)[o[1:10],]

#######---------#########

colo2<-c(rep("blue4",3), rep("blue",3), rep("deepskyblue4",3),rep("deepskyblue1",3),
         rep("red2",3),rep("lightcoral",3),rep("yellow",3),rep("wheat2",2))
barplot(dr.rpkm["Smp_179600.1", colnames(dr$counts)], col = colo2, las = 2)

## loop for plotting
G<-scan(file="1.1.txt", what=character())
i<-1;
while (i<=length(G)) {
       pdf(paste(G[i],".pdf", sep=""));
       barplot(dr.rpkm[G[i], colnames(dr$counts)], col = colo2, las = 2, main=G[i]); 
       dev.off();
       i<-i+1
}
dev.new()
## all pdf in a file
i<-1;
pdf('1.2.pdf')
for (i in 1:length(G)){
  barplot(dr.rpkm[G[i], colnames(dr$counts)], col = colo2, main=G[i])}
  dev.off()

dr.rpkm["Smp_179600.1",colnames(d$counts)]

## rpkm distribution
hist(dr.rpkm[,"bM1"], breaks = c(0, 2, 5, 10, 300000), plot = FALSE) # sample cpm histogram 

#### export glm and plot
#summary(dt<-decideTestsDGE(lrt)) #Error because of multi-cmps
TT<-topTags(lrt, n=Inf)
write.table(TT, file="allDEG_glm_Pairing-Gender-Tissue_rmdup.txt")
#write.table(topTags(lrt, n=Inf), file="file.txt)
## plot
detags<-rownames(d)[as.logical(dt)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1,1),col="blue")


###########################################################
#### SPECIFIC GROUP ###
rd2 <-read.delim("1.2.2_readcount+length_rpkm2-filtered.table", row.names="gene_ID")
grpt<- as.factor(substr(colnames(rd2[,7:12]), 1,2))
dt<-DGEList(counts=rd2[,7:12], group=grpt)# genes=rawdata[,1]) ####SAMPLE ORDER#####
#ds<-DGEList(counts=rawdata[, c("sF1", "sF2", "sF3", "sO1", "sO2")], group=grp2)
dt$samples; dim(dt)

pairs(log2(dt$counts), pch = '.', col = "blue", lower.panel=NULL) ##overview of sample relationship

################ Pairwise DGE #################
#Filtering
keep <- rowSums(rpkm(dt, dr.length)>2) >= 3 
dt <- dt[keep,]
table(keep)
dim(dt)
dt$samples$lib.size <- colSums(dt$counts) #re-calc lib size

#Normalization
dt <- calcNormFactors(dt)
dt$samples
#dt.rpkm<-rpkm(dt, dr.length, normalized.lib.sizes = TRUE, log=FALSE)
colx <- as.numeric(dt$samples$group)
plotMDS(dt, col = colx) # method="bcv", default method is leading Log Fold Change(lLFC)
#### other methods ## BCV "biological coefficient of variation"
par(mfrow=c(2,2))
plotMDS(dt,col=colx, main="500 / lLFC")
plotMDS(dt,col=colx, method="bcv", main="500 / BCV") 
plotMDS(dt,col=colx, top=2000, main="2000 / lLFC")
plotMDS(dt,col=colx, top=2000, method="bcv", main="2000 / BCV")

summary(dt$counts)
summary(cpm(dt))

### GET NORMALIZED+Filtered CPM ##
# compute CPM or RPKM
CPMn<-cpm(ds, normalized.lib.sizes=TRUE, log=FALSE) #cpm default norm
#cpm(d, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
write.table(CPMn, file="cpm_so12-bo13_5filtered_normalized.txt")
mean.s1 <- rowMeans(CPMn[, c(1, 2)])  #means of sTe "sT1"
mean.s2 <- rowMeans(CPMn[, c(3, 4, 5)])  #means of bTe
write.table(mean.s1, file = "cpm_sOv_mean.txt")
write.table(mean.s2, file = "cpm_bOv_mean.txt")

#Estimating the dispersions
ds <- estimateCommonDisp(ds) #, verbose=TRUE)
ds <-estimateTagwiseDisp(ds)
plotBCV(ds)

plotMeanVar(ds, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)

#Differential expression
et<-exactTest(ds, pair=c("sO", "bO"))
write.table(topTags(et, n=Inf), file = "pw_bO-bF_1.txt")
topTags(et, n=20, sort.by=c("logFC"))

# hist of logFC
FDR <- p.adjust(lrt$table$PValue, method="BH")
hist(et$table[,"logFC"], plot=FALSE)

## check with cpm
cps<-cpm(ds, normalized.lib.sizes=TRUE)  ##normalized
summary(cps)
barplot(cps["Smp_179600.1", colnames(ds$counts)], las = 2)
cps["Smp_179600.1", colnames(ds$counts)]
hist(cps[,"bO1"], breaks = c(0, 2, 5, 10, 300000), plot = FALSE)

TT<-topTags(et, n=Inf)
write.table(TT, file="OvariesPW_cpm5.txt") #!CHANGE

summary(dt<-decideTestsDGE(et, p=0.001, adjust ="BH"))
#number of DE genes at 5% FDR; -1 down-reg, 1 up-reg
detags<-rownames(ds)[as.logical(dt)]
plotSmear(et, de.tags=detags)  ## adjust p above
abline(h=c(-0.585,0.585),col="blue")
abline(h=c(-1,1),col="blue")

##########################################
### import all logFC and FDR for plot
x<- read.table ("OvaryPW_s12_cpm5_comp.txt", header = T, sep = " ")
head(x)
logFC<-x[,2]
FDR<-x[,5]

plot(logFC, -log10(FDR), xlab = 'log2 fold change', ylab = '-log10 FDR', pch = '.', col = 'blue')
abline(v=c(-0.585,0.585),col='red')
abline(h=3, col = 'red')
## histogram ##
hist(logFC, plot = FALSE)
hist(logFC, breaks = c(-14, -12, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, -0.585,
                        0, 0.585, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14), 
                        freq = TRUE, ylim = c(1, 2000), border = 'blue')
hist(logFC, breaks = 1)
## set x axis intervals, use frequency/density
hist(FDR, plot = FALSE)
hist(FDR, breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), freq = TRUE)
plot(logFC, pch = '.', col = 'blue')

##### extract logFC 0.585 && FDR 0.001 ###
y<-read.table("ovDEG_s12-b13_3299.txt", sep='\t'); head(y)
logFC_2 <- y[,2]; head(logFC_2)
FDR_2 <- y[,3]; head(FDR_2)
plot(logFC_2, -log10(FDR_2), pch = 20, col = 'blue')
hist(FDR_2, plot = FALSE)
hist(logFC_2, plot = FALSE)
hist(logFC_2, breaks = 20, border = 'blue')
hist(logFC_2, breaks = c(-14,-13,-12, -11,-10,-9,-8,-7,
    -6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), border = 'blue')

