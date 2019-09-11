## usage: Rscript topGO.R [go annotation with scores]

library("topGO")

## useful discussions:
# https://www.biostars.org/p/92374/

# set parameters
ndSize <- 10 
goCat <- 'BP'
nTerms <-50 # number of terms to write out (or only show significant ones)

# reference file containing go terms and scores
args = commandArgs(trailingOnly=TRUE)
refGO <- read.table(file=args[1], sep=" ", stringsAsFactor=F) 
# 
# Smp_000030 0.05 GO:0000502,GO:0005488,GO:0030234,GO:0042176
# Smp_000040 0.08 GO:0003777,GO:0005515,GO:0005871
#
names(refGO) = c('id', 'score', 'go')
ref.vec = strsplit(refGO$go, split=',', fixed=T)
names(ref.vec) <- refGO$id

geneList<-refGO$score
names(geneList) <- refGO$id

# select genes of interest: can only do cutoff based on scores 
# see: https://support.bioconductor.org/p/54268/
mySel<-function (score) {
	return(score>=0.85)
}

myGOdata <- new("topGOdata",
		description = "GO enrichment",
		ontology = goCat,
		allGenes = geneList, # specifies all annotated genes
		geneSel = mySel, # your genes of selection
		annot = annFUN.gene2GO,
		gene2GO = ref.vec,
		nodeSize = ndSize, # minimal number of genes for a term
)

myGOdata

resultClassicFisher <- runTest(myGOdata, algorithm="classic", statistic="Fisher") # classic algorithm doesn't consider GO hierarchy
resultWeight01Fisher <- runTest(myGOdata,algorithm="weight01",statistic="Fisher") # 
resultWeight01KS <- runTest(myGOdata,algorithm="weight01",statistic="ks") # 
#resultElim <- runTest(myGOdata,algorithm="elim",statistic="Fisher") # even less true positives

#resultClassicFisher
resultWeight01Fisher
resultWeight01KS
#resultElim

allRes <- GenTable(
  myGOdata,
  classicFisher = resultClassicFisher,
  weight01Fisher = resultWeight01Fisher, 
  weight01KS = resultWeight01KS, 
  #elim = resultElim, 
  orderBy = "weight01KS",
  ranksOf = "weight01KS",
  topNodes = nTerms # check:--> resultClassic --> "711 GO terms scored"
)

# add column with genes; add to table
# see the topGodata class: https://www.rdocumentation.org/packages/topGO/versions/2.24.0/topics/topGOdata-class
allRes$genes <- sapply(allRes$GO.ID, function(x)
    {
      genes<-genesInTerm(myGOdata, x) 
      genes[[1]][genes[[1]] %in% sigGenes(myGOdata)] # significant genes
    })

# convert the gene list to a character vector
allRes$genes <-vapply(allRes$genes, paste, collapse = ",", character(1L))

# only write significant terms
allRes<-subset(allRes, as.numeric(allRes[,"topGO"])<0.05 | grepl("e-", allRes[,"topGO"]))
#allRes<-allRes[ which(as.numeric(allRes[,"topGO"])<0.05 | grepl("e-", allRes[,"topGO"])),]

outfile <- paste("topgoKS_", args[1], "_", goCat, "_", ndSize, ".txt", sep="")
write.table(allRes, outfile,sep="\t", quote=F, row.names = F) #

## visualise
# showSigOfNodes(myGOdata, score(resultTopgo), firstSigNodes = 5, useInfo ='all')
# printGraph(myGOdata, resultTopgo, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

## retrive genes in a GO term
#allGO = genesInTerm(myGOdata)
#allGO["GO:0006732"]

## all annotated genes in significantly enriched terms
#sel.go <- names(score(resultTopgo))[1:5]
#selcTerm <- allRes$GO.ID[which(allRes$topGO<0.05)]
#selcGenes <- genesInTerm(myGOdata, whichGO=selcTerm)
#selcGenes
