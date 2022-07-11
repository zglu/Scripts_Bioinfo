# R script to calculate distance between gene A start to the end coord of previous gene

################## method 1 ######################
smv9.genes<-read.delim("all_genes_ordered.txt", sep="\t", header=F)
colnames(smv9.genes)<-c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# gene list sorted by chromosome and coordinate (-1,1 -nk4)
#SM_V9_1	Liftoff	gene	68427	68783	.	-	.	ID=Smp_329140
#SM_V9_1	Liftoff	gene	83503	131262	.	-	.	ID=Smp_315690
#SM_V9_1	Liftoff	gene	216287	251096	.	+	.	ID=Smp_317470

contiglen<-read.delim("chr-length.txt", header=T, sep="\t", row.names = "chr")
#chr	length
#SM_V9_1	87984036
#SM_V9_2	45716228

# output the calculation to a file
sink("Calculated_dist.txt")
sink(stdout(), type = "message")

# ***need to calculate and add first and last gene manually***
for (row in 2:nrow(smv9.genes)) {
  # + strand compared to the previous gene
  if (smv9.genes[row, "strand"] == "+") {
    # if the former gene is on the same chr
    if (smv9.genes[row, "chr"] == smv9.genes[row-1, "chr"]) {
      message(smv9.genes[row, "chr"], "\t", smv9.genes[row-1, "end"]+1, "\t", smv9.genes[row, "start"]-1, "\t", "+", "\t", smv9.genes[row, "attributes"])
  } else {
    # if the previous gene is not on the same chr, it means that the current gene is the first (take to chr start 1)
      message(smv9.genes[row, "chr"], "\t", 1, "\t", smv9.genes[row, "start"]-1, "\t", "+", "\t", smv9.genes[row, "attributes"])
    }
  }  else { # - strand: compared to the next gene
    if (smv9.genes[row, "chr"] == smv9.genes[row+1, "chr"]) {
      message(smv9.genes[row, "chr"], "\t", smv9.genes[row, "end"]+1, "\t", smv9.genes[row+1, "start"]-1, "\t", "-", "\t", smv9.genes[row, "attributes"])
    } else {
      # if the previous gene is not on the same chr, it means that the current gene is the last (take to the chr end)
      message(smv9.genes[row, "chr"], "\t", smv9.genes[row, "end"]+1, "\t", contiglen[as.character(smv9.genes[row, "chr"]), "length"], "\t", "-", "\t", smv9.genes[row, "attributes"])
    }
  }
}
sink()

# manually add the first gene
# SM_V9_1   68784   83502   -   ID=Smp_329140

##################################################
################## method 2 ######################

# use dplyr lag function
# use data.table shift function

