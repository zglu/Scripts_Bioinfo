library(pathview)

# ARGS: gene list and pathway id
args = commandArgs(trailingOnly=TRUE)

degk<-read.delim(args[1], sep=" ", header = F)# avoid duplicate K numbers
rownames(degk)<-degk$V1
degk[,1] <- NULL
head(degk)

# id        log2FC
#K00006 -0.7641832
#K00026 -0.7620600
#K00030 -0.7097674
#K00108  1.2265815

# have to use as.matrix same.layer=F gives out K numbers; =T(default) gives abbr names; kegg.native=T gives map; =F gives network with K
pv.out1 <- pathview(gene.data = as.matrix(degk)[, 1], pathway.id = args[2], species = "ko", kegg.native = T, gene.idtype = "KEGG", same.layer=T, out.suffix = args[1], limit = c(0,24), bins = 12, mid = list(gene="yellow"))

#batch ko mapping
#GI<-scan(file="../ref7.1/ko.list.txt", what=character());
#i<-1;
#for (i in 1:length(GI)) {
#	pv.out.ko <- pathview(gene.data = as.matrix(degk)[, 1], pathway.id = GI[i], species = "ko", out.suffix = "feCycling", kegg.native = T, gene.idtype = "KEGG")
#	}
#  dev.off()

## multi-samples
# multideg<-read.delim("~/Dropbox/port/KEGG/Pairing-K-all_ed1130", sep=" ", row.names="ID")
# pv.outs<-pathview(gene.data = as.matrix(multideg)[, c(3,4)], pathway.id = "04215", species = "ko", out.suffix = "FeOvDEG", kegg.native = T, gene.idtype = "KEGG") #, same.layer=F ==K
# 
# GI<-scan(file="~/Dropbox/Scripts/ko.list", what=character());
# i<-1;
# for (i in 1:length(GI)) {
#   pv.out.kos <- pathview(gene.data = as.matrix(multideg)[, c(1,3)], pathway.id = GI[i], species = "ko", out.suffix = "MaFeDEG", kegg.native = T, gene.idtype = "KEGG")
# }
# dev.off()
