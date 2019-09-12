args = commandArgs(trailingOnly=TRUE)
fishertable <- read.delim(paste0(args[1],"-fisherTable.txt"), sep=" ", header=F)
## should only test those with mapped gene > 0 if with all genes
# fishertable <- fishertable[which(fishertable$V2 > 0),]
no.tests<-nrow(fishertable)
expotable<-data.frame()
for (i in 1:no.tests) {
  fisherp<-fisher.test(matrix(c(fishertable[i,'V3'],fishertable[i,'V4'],fishertable[i,'V5'],fishertable[i,'V6']),2,2), alternative ="greater")$p.value  # testing for over-representation (rather than under-representation) the alternative parameter is set to "greater"
  newline<-cbind(fishertable[i,], round(fisherp, digits = 4)) # 
  expotable<-rbind(expotable, newline)
}
colnames(expotable)<-c("accession", "term", "mapped-test", "mapped-all", "unmapped-test", "unmapped-all", "genes", "p_value")
## p-value adjustment for multiple comparisons ?p.adjust
### control of the family-wise error rate: bonferroni, holm, hochberg, hommel
expotable$Bonferroni <- round(p.adjust(expotable$p_value, method = "bonferroni"), digits=4) # p-values are multiplied by the number of comparisons (no.test)
expotable$Holm <- round(p.adjust(expotable$p_value, method = "holm"), digits=4)
### control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses. Less stringent 
expotable$FDR <- round(p.adjust(expotable$p_value, method = "BH"), digits=4) # Benjamini & Hochberg
#expotable$FDR<-NULL # drop a column
expotable<-expotable[order(expotable$FDR),]
expotable<-expotable[, c(1:6,8:11, 7)]
expotable<-subset(expotable, as.numeric(expotable[,"FDR"])<0.05 | grepl("e-", expotable[,"FDR"]))
write.table(expotable, file=paste0("fisher_", args[1], "_Results.txt"), row.names = F, sep="\t", quote = F)
