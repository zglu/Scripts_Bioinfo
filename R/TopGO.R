## usage: Rscript topGO.R [gene ids]

library("topGO")

# this is the file containing gene IDs
args = commandArgs(trailingOnly=TRUE)
ids <- read.table(args[1], header=F)
#ids <- read.table("pzq.genes", header=F)
myGenes <- as.character(ids$V1)

ndSize <- 10 
goCat <- 'BP'
nTerms <-50 # number of terms to write out (or only show significant ones)

# GO annotation file
refGO <- read.table(file="femaleGORef-7.2.txt", sep=" ", stringsAsFactor=F) 
# 
# Smp_000030 GO:0000502,GO:0005488,GO:0030234,GO:0042176
# Smp_000040 GO:0003777,GO:0005515,GO:0005871
#
names(refGO) = c('id', 'go')
ref.vec = strsplit(refGO$go, split=',', fixed=T)
names(ref.vec) <- refGO$id

allAnnotated <- refGO$id
geneList <- factor(as.integer(allAnnotated %in% myGenes))
names(geneList) <- allAnnotated

myGOdata <- new("topGOdata",
		description = "topGO",
		ontology = goCat,
		allGenes = geneList, # specifies all annotated genes and of which those are yours
		annot = annFUN.gene2GO,
		gene2GO = ref.vec,
		nodeSize = ndSize, # can change this!!
)

myGOdata

resultClassic <- runTest(myGOdata, algorithm="classic", statistic="Fisher") # classic algorithm doesn't consider GO hierarchy
resultTopgo <- runTest(myGOdata,algorithm="weight01",statistic="Fisher") # 
#resultElim <- runTest(myGOdata,algorithm="elim",statistic="Fisher") # even less true positives

#resultClassic
resultTopgo
#resultElim

allRes <- GenTable(
  myGOdata,
  classic = resultClassic,
  topGO = resultTopgo, 
  #elim = resultElim, 
  orderBy = "topGO",
  ranksOf = "topGO",
  topNodes = nTerms # check:--> resultClassic --> "711 GO terms scored"
)

# add column with genes; add to table
allRes$genes <- sapply(allRes$GO.ID, function(x)
    {
      genes<-genesInTerm(myGOdata, x) 
      genes[[1]][genes[[1]] %in% myGenes]
    })
#allRes$genes[which(allRes$topGO<0.05)]

# convert the gene list to a character vector
allRes$genes <-vapply(allRes$genes, paste, collapse = ",", character(1L))
# only write significant terms
allRes<-subset(allRes, as.numeric(allRes[,"topGO"])<0.05 | grepl("<", allRes[,"topGO"]))
#allRes<-allRes[ which(as.numeric(allRes[,"topGO"])<0.05 | grepl("e-", allRes[,"topGO"])),]

outfile <- paste("topgo_", args[1], "_", goCat, "_", ndSize, ".txt", sep="")
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
