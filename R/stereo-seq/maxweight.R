# read csv results with composition info; take the col with max value


args<-commandArgs(T)

results<-read.csv(args[1], row.names = "X")

#              CAFs    NCSC.like     immune.1    immune.2 keratinocytes.1 keratinocytes.2 keratinocytes.3 melanoma.pigmented
#38444 0.0006538332 0.0009147615 0.0017283390 0.002543489     0.155293495      0.07990093     0.004968082       1.391697e-02
#30499 0.0002460700 0.0000249803 0.0034521020 0.101345545 

maxcol<-colnames(results)[max.col(results, ties.method = "first")]
maxweight<-cbind(rownames(results), maxcol); 
maxweight<-as.data.frame(maxweight)
table(maxweight$maxcol)

