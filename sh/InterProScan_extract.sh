#!/bin/bash

## Extract GO terms, Pfam, and InterPro accessions from the output gff3 of InterProScan
## Usage ./InterProScan_extract.sh [annotation.gff3] [gene/mRNA]
## zl3 28.05.2019

if [ $2 == 'mRNA' ] 
  then
    ### extract GO terms
    cat $1 | grep GO | awk '{print $9 "\t" $11}'| sed 's/;/ /g'| awk '{print $2 "\t" $4}'| sed 's/=/ /g'| awk '{print $2 "\t" $4}'| tr -d \"|sort | uniq | sed 's/,/ /g'|awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | sort | uniq | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}'| tr '\t' ','| sed 's/,/ /' > ips_GO_mRNA.txt
    
    ### extract InterPro accessions
    cat $1 | grep 'InterPro:'| sed 's/;/ /g'| awk '{print $10 "\t" $NF}'| sed 's/,/ /' | tr -d \" | sed 's/=/ /; s/:/ /'| awk '{print $2 " " $4}'| sort | uniq| awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | tr '\t' ','| sed 's/,/ /' > ips_InterPro_mRNA.txt
    
    ### Extract pfam
    cat $1 | grep Name=PF |awk '{ for(i=9; i<NF; i++) printf "%s",$i OFS; if(NF) printf "%s",$NF; printf ORS}'| sed 's/ /-/g'|sed 's/;/ /g'| awk '{print $2 "\t" $6}'| sed 's/-/ /; s/=/ /g'| awk '{print $2 "\t" $5}'| sort | uniq |grep PF > pfam1.txt
    cat $1 | grep Name=PF |awk '{ for(i=9; i<NF; i++) printf "%s",$i OFS; if(NF) printf "%s",$NF; printf ORS}'| sed 's/ /-/g'|sed 's/;/ /g'| awk '{print $2 "\t" $5}'| sed 's/-/ /; s/=/ /g'| awk '{print $2 "\t" $5}'| sort | uniq| grep PF > pfam2.txt
    cat pfam1.txt pfam2.txt | sort | uniq | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | tr '\t' ','| sed 's/,/ /' > ips_Pfam_mRNA.txt
    rm -f pfam1.txt pfam2.txt
    echo "Exported at mRNA level!"
elif [ $2 == 'gene' ]
  then
    ### extract GO terms
    cat $1 | grep GO | awk '{print $9 "\t" $11}'| sed 's/;/ /g'| awk '{print $2 "\t" $4}' | sed 's/=/ /g'| awk '{print $2 "\t" $4}'| tr -d \"|sed 's/,/ /g' | awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | sed 's/\./ /'| awk '{print $1 "\t" $3}'|sort -u | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | tr '\t' ','| sed 's/,/ /' > ips_GO_gene.txt

    ### extract InterPro accessions
    cat $1 | grep 'InterPro:'| sed 's/;/ /g'| awk '{print $10 "\t" $NF}'| sed 's/,/ /' | tr -d \" | sed 's/=/ /; s/:/ /'| awk '{print $2 " " $4}'|sed 's/\./ /'|awk '{print $1 "\t" $3}'| sort -u | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | tr '\t' ','| sed 's/,/ /' > ips_InterPro_gene.txt

    ### Extract pfam
    cat $1 | grep Name=PF |awk '{ for(i=9; i<NF; i++) printf "%s",$i OFS; if(NF) printf "%s",$NF; printf ORS}'| sed 's/ /-/g'|sed 's/;/ /g'| awk '{print $2 "\t" $6}'| sed 's/-/ /; s/=/ /g'| awk '{print $2 "\t" $5}'| sort | uniq |grep PF > pfam1.txt
    cat $1 | grep Name=PF |awk '{ for(i=9; i<NF; i++) printf "%s",$i OFS; if(NF) printf "%s",$NF; printf ORS}'| sed 's/ /-/g'|sed 's/;/ /g'| awk '{print $2 "\t" $5}'| sed 's/-/ /; s/=/ /g'| awk '{print $2 "\t" $5}'| sort | uniq| grep PF > pfam2.txt
    cat pfam1.txt pfam2.txt | sed 's/\./ /'|awk '{print $1 "\t" $3}'| sort -u | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | tr '\t' ','| sed 's/,/ /' > ips_Pfam_gene.txt
    rm -f pfam1.txt pfam2.txt
    echo "Exported at gene level!"
else
    echo "Need to provide the gff3 file and level (gene/mRNA)"
fi
