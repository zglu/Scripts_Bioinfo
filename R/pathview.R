### R code from vignette source 'pathview.Rnw'

###################################################
### code chunk number 1: synopsis1 (eval = FALSE)
###################################################
## library(pathview)
## data(gse16873.d)
## pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
##                    species = "hsa", out.suffix = "gse16873")


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("pathview")


###################################################
### code chunk number 3: <install (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("Rgraphviz", "png", "KEGGgraph", "org.Hs.eg.db"))


###################################################
### code chunk number 4: install (eval = FALSE)
###################################################
## install.packages("pathview",repos="http://R-Forge.R-project.org")


###################################################
### code chunk number 5: install (eval = FALSE)
###################################################
## install.packages("/your/local/directory/pathview_1.0.0.tar.gz", 
##     repos = NULL, type = "source")


###################################################
### code chunk number 6: <install (eval = FALSE)
###################################################
## install.packages("/your/local/directory/XML_3.95-0.2.zip", repos = NULL)

############ self experience: first XML then pathview #########
# sudo apt-get update
# sudo apt-get install r-cran-xml

library(pathview)
setwd("~/Documents/KEGG-re/")
degk<-read.delim("~/Documents/KEGG-re/fedeg-k_1.txt", sep="\t", row.names="ID") # avoid duplicate K numbers
head(degk)

#           log2FC
#K00006 -0.7641832
#K00026 -0.7620600
#K00030 -0.7097674
#K00108  1.2265815

degsmp<-read.delim("~/Documents/KEGG-re/fedeg-smp.txt", sep="\t", row.names = "Smp")

# have to use as.matrix same.layer=F gives out K numbers; =T(default) gives abbr names; kegg.native=T gives map; =F gives network with K
pv.out1 <- pathview(gene.data = as.matrix(degk)[, 1], pathway.id = "00020", species = "ko", out.suffix = "maDEG", kegg.native = T, gene.idtype = "KEGG", same.layer=F) # can also try kegg.native=F
pv.out2 <- pathview(gene.data = as.matrix(degsmp)[, 1], pathway.id = "04114", species = "smm", out.suffix = "DEG-smm", kegg.native = T, gene.idtype = "KEGG") # can also try kegg.native=F

#batch ko mapping
GI<-scan(file="~/Dropbox/Scripts/ko.list", what=character());
i<-1;
for (i in 1:length(GI)) {
	pv.out.ko <- pathview(gene.data = as.matrix(degk)[, 1], pathway.id = GI[i], species = "ko", out.suffix = "feDEG", kegg.native = T, gene.idtype = "KEGG")
	}
  dev.off()

# batch smm pathways
GI2<-scan(file="~/Dropbox/Scripts/smm.list", what=character());
i<-1;
for (i in 1:length(GI)) {
  pv.out.smm <- pathview(gene.data = as.matrix(degsmp)[, 1], pathway.id = GI2[i], species = "smm", out.suffix = "feDEG", kegg.native = T, gene.idtype = "KEGG")
  }
  dev.off()  
  
## multi-samples
multideg<-read.delim("~/Dropbox/port/KEGG/Pairing-K-all_ed1130", sep=" ", row.names="ID")
pv.outs<-pathview(gene.data = as.matrix(multideg)[, c(3,4)], pathway.id = "04215", species = "ko", out.suffix = "FeOvDEG", kegg.native = T, gene.idtype = "KEGG") #, same.layer=F ==K

GI<-scan(file="~/Dropbox/Scripts/ko.list", what=character());
i<-1;
for (i in 1:length(GI)) {
  pv.out.kos <- pathview(gene.data = as.matrix(multideg)[, c(1,3)], pathway.id = GI[i], species = "ko", out.suffix = "MaFeDEG", kegg.native = T, gene.idtype = "KEGG")
}
dev.off()
