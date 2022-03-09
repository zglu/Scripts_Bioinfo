# combine gene information and export to a .md file each gene for Hugo

PROD <-read.delim("all_prod.txt", sep="\t", row.names="id")
pro <-as.matrix(PROD)

#id	product
#Smp_000020	Protein O-GlcNAcase

PSIDA <- read.delim("all_psid.txt", sep="\t", row.names="id")
psid <- as.matrix(PSIDA)

#id	psid
#Smp_000020
#Smp_000030
#Smp_000040	Smp_122080

CDOM <-read.delim("all_cdd.txt", sep="\t", header = T, row.names="id")
cdd <-as.matrix(CDOM)

#id	cdd
#Smp_000020	"NAGidase"
#Smp_000030	"PC_rep", "RPN1"

KEGG <-read.delim("all_kegg.txt", sep="\t", header = T, row.names="id")
kk <-as.matrix(KEGG)

#id	kegg	ko
#Smp_000020	K15719	Insulin resistance
#Smp_000030	K03028	Epstein-Barr virus infection, Proteasome
#Smp_000040	K10407	Salmonella infection
#Smp_000050

RPKM <-read.delim("rpkmAllMean.txt", sep=" ", header = T, row.names="id")
rpkmm <-as.matrix(RPKM)

#id bsTe bsOv Egg Mir Spo Cerc mCer fCer Som3 Som24 mSo1 fSo1 mSo2 fSo2 mSo3 fSo3 d18.ssMa d18.ssFe d21.ssMa d21.ssFe d28.ssMa d28.ssFe d35.ssMa d35.ssFe d38.ssMa d38.ssFe d42plus.ssMa d42plus.ssFe d21.bsMa d21.bsFe d28.bsMa d28.bsFe d35.bsMa d35.bsFe d38.bsMa d38.bsFe d42plus.bsMa d42plus.bsFe
#Smp_329140 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.119418603501414

COOR <-read.delim("all_coords.txt", sep="\t", header=T, row.names="id")
coord<-as.matrix(COOR)

#id	chr	start	end
#Smp_000020	SM_V7_1	46644182	46672976

GI<-readLines("all.genes") #file("stdin"),1

#Smp_331390
#Smp_331400
#...

i<-1;
for (i in 1:length(GI)) {
  fileConn<-file(paste(GI[i],".md", sep=""));
  writeLines(c("+++", paste("title = ", "\"", GI[i], "\"", sep=""), 
               paste("description = ", "\"", pro[GI[i],1], "\"", sep=""),
               paste("psid = ", "\"", psid[GI[i],1], "\"", sep=""),
               paste("pType = ", "\"", "Gene", "\"", sep=""), 
               paste("chr = ", "\"", coord[GI[i],1], "\"", sep=""), 
               paste("gstart = ", "\"", coord[GI[i],2], "\"", sep=""), 
               paste("gend = ", "\"", coord[GI[i],3], "\"", sep=""), 
               paste("cdd = ", "\"", cdd[GI[i],1], "\"", sep=""), 
               paste("kegg = ", "\"", kk[GI[i],1], "\"", sep=""), 
               paste("pathway = ", "\"", kk[GI[i],2], "\"", sep=""), 
	       paste0("bsTe = ", rpkmm[GI[i],1]),
	       paste0("bsOv = ", rpkmm[GI[i],2]),
	       paste0("Egg = ", rpkmm[GI[i],3]),
	       paste0("Mir = ", rpkmm[GI[i],4]),
	       paste0("Spo = ", rpkmm[GI[i],5]),
	       paste0("Cerc = ", rpkmm[GI[i],6]),
	       paste0("mCer = ", rpkmm[GI[i],7]),
	       paste0("fCer = ", rpkmm[GI[i],8]),
	       paste0("Som3 = ", rpkmm[GI[i],9]),
	       paste0("Som24 = ", rpkmm[GI[i],10]),
	       paste0("mSo1 = ", rpkmm[GI[i],11]),
	       paste0("fSo1 = ", rpkmm[GI[i],12]),
	       paste0("mSo2 = ", rpkmm[GI[i],13]),
	       paste0("fSo2 = ", rpkmm[GI[i],14]),
	       paste0("mSo3 = ", rpkmm[GI[i],15]),
	       paste0("fSo3 = ", rpkmm[GI[i],16]),
	       paste0("d18_ssMa = ", rpkmm[GI[i],17]),
	       paste0("d18_ssFe = ", rpkmm[GI[i],18]),
	       paste0("d21_ssMa = ", rpkmm[GI[i],19]),
	       paste0("d21_ssFe = ", rpkmm[GI[i],20]),
	       paste0("d28_ssMa = ", rpkmm[GI[i],21]),
	       paste0("d28_ssFe = ", rpkmm[GI[i],22]),
	       paste0("d35_ssMa = ", rpkmm[GI[i],23]),
	       paste0("d35_ssFe = ", rpkmm[GI[i],24]),
	       paste0("d38_ssMa = ", rpkmm[GI[i],25]),
	       paste0("d38_ssFe = ", rpkmm[GI[i],26]),
	       paste0("d42plus_ssMa = ", rpkmm[GI[i],27]),
	       paste0("d42plus_ssFe = ", rpkmm[GI[i],28]),
	       paste0("d21_bsMa = ", rpkmm[GI[i],29]),
	       paste0("d21_bsFe = ", rpkmm[GI[i],30]),
	       paste0("d28_bsMa = ", rpkmm[GI[i],31]),
	       paste0("d28_bsFe = ", rpkmm[GI[i],32]),
	       paste0("d35_bsMa = ", rpkmm[GI[i],33]),
	       paste0("d35_bsFe = ", rpkmm[GI[i],34]),
	       paste0("d38_bsMa = ", rpkmm[GI[i],35]),
	       paste0("d38_bsFe = ", rpkmm[GI[i],36]),
	       paste0("d42plus_bsMa = ", rpkmm[GI[i],37]),
	       paste0("d42plus_bsFe = ", rpkmm[GI[i],38]),
               "+++"), fileConn)
}
  close(fileConn)

#formatC(pro[GI[i],4], width=6, format="d", flag="0")
# paste0("abbr = ", "\"", sprintf("%06d", as.numeric(pro[GI[i],4])), "\""),
# paste numbers and keep 000 at the beginning
