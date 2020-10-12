## genome files 
- @ /lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/split
### split genome
perl ~zl3/scripts/split.pl ../../Smansoni_v7.masked.masked.fa

### summary 
- bsub.py 4 summaryATGC summarizeACGTcontent.pl /lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Smansoni_v7.masked.masked.fa
- result in summaryATGC.o

### make chr.lst with content for all splits
- /lustre/scratch118/infgen/team133/alt/SCHISTO/augustus_V7_masked_star_prioritised_iso/split/SM_V7_1.fa  /lustre/scratch118/infgen/team133/alt/SCHISTO/V7_STAR_FINAL/hints_pri.gff       1       88881357
- cat summaryATGC.o |grep 'SM_V7_' |awk '{print "/lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/split/" $3 ".fa" "\t" "/lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/hints/hints_pri.gff" "\t" "1" "\t" $1}'> chr.lst

### chrom size
- cat summaryATGC.o | awk '{print $3 " " $1}'|grep 'SM_V7'>CHROMSIZE

## bam files and hints file 
- @ /lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/hints
### merge rnaseq and isoseq bam files
- samtools merge merged_rnaseq.bam 0HR.bam 12HR.bam 24HR.bam 3HR.bam 48HR.bam 6HR.bam 72HR.bam BF.bam BM.bam CERC.bam D13.bam D17.bam D21.bam D28.bam D35.bam D6.bam MIRACIDIA.bam SPOROCYST.bam
- COPY ISOSEQ-BAMS >cp /nfs/repository/working_area/SHISTO/V7/Isoseq_masked/Sman_iso_*.bam isoseq_bam/
- bsub.py 6 merge_isoseq samtools merge ../merged_isoseq.bam Sman_iso_fem_V7.sorted.bam Sman_iso_mal_V7.sorted.bam Sman_iso_som_V7.sorted.bam Sman_iso_sporo_nopoly_V7.sorted.bam

### sort merged bam and index them (or sort each bam file beforehand)
- samtools-1.3 sort --threads 8 -o merged_isoseq_sorted.bam merged_isoseq.bam
- samtools index merged_isoseq_sorted.bam (/lustre/scratch118/infgen/team133/alt/SCHISTO/V7_STAR_FINAL/merged_ALL.sorted.bam)

### bam2hints
- bsub.py 6 bam2hints_rnaseq /lustre/scratch118/infgen/archive/ss34/SCHISTO/augustus/augustus-3.2.2/bin/bam2hints --intronsonly --in=/lustre/scratch118/infgen/team133/alt/SCHISTO/V7_STAR_FINAL/merged_ALL.sorted.bam --out=intron_rnaseq.gff
- bsub.py 6 bam2hints_isoseq /lustre/scratch118/infgen/archive/ss34/SCHISTO/augustus/augustus-3.2.2/bin/bam2hints --intronsonly --in=/lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/hints/merged_isoseq_sorted.bam --out=intron_isoseq.gff
- change 'pri=4;' to 'pri=40;'in the last column (default pri=4)
- cat intron_isoseq.gff | sed 's/pri=4;/pri=40;/g'>intro_isoseq_ed.gff

### bam2wig
- bsub.py 12 -q normal bam2wig_rnaseq bam2wig.py -i /lustre/scratch118/infgen/team133/alt/SCHISTO/V7_STAR_FINAL/merged_ALL.sorted.bam -s /lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/split/CHROMSIZE -o wig_out_rnaseq
- bsub.py 12 bam2wig_isoseq bam2wig.py -i /lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/hints/merged_isoseq_sorted.bam -s /lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/split/CHROMSIZE -o wig_out_isoseq

### wig2hints with priority settings
#### ./wig2hints_rnaseq 
- cat wig_out_rnaseq.wig | /software/grit/augustus-3.2.2/scripts/wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --radius=4.5 --pri=4 --strand="." > exonparts_rnaseq.gff
- chmod u+x wig2hints_rnaseq
- bsub.py 8 wig2hints_rna ./wig2hints_rnaseq
#### ./wig2hints_isoseq
- cat wig_out_isoseq.wig | /software/grit/augustus-3.2.2/scripts/wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --radius=4.5 --pri=40 --strand="." > exonparts_isoseq.gff
- chmod u+x wig2hints_isoseq
- bsub.py 6 wig2hints_iso ./wig2hints_isoseq

### make complete hints.gff
- cat intron_rnaseq.gff intron_isoseq_ed.gff exonparts_rnaseq.gff exonparts_isoseq.gff > hints_pri.gff

## run Augustus
- @ /lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/augustus

### Create Augustus job list in folder /jobs:
- createAugustusJoblist.pl --sequences=../split/chr.lst --wrap="#" --overlap=50000 --chunksize=3000000 --outputdir=/lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/Augustus_Run_plusIsoseq2/augustus/ --joblist=jobs.lst --jobprefix=SM_V7_aug_ --command "/lustre/scratch118/infgen/archive/ss34/SCHISTO/augustus/augustus-3.2.2/bin/augustus --species=schistosoma2 --UTR=1 --nc=0 --gff3=on --alternatives-from-evidence=1 --extrinsicCfgFile=/lustre/scratch118/infgen/archive/ss34/SCHISTO/augustus/augustus-3.2.2/config/extrinsic/extrinsic.M.RM.E.W.cfg --AUGUSTUS_CONFIG_PATH=/lustre/scratch118/infgen/archive/ss34/SCHISTO/augustus/augustus-3.2.2/config/"

### run jobs /jobs
- for i in `cat jobs.lst`; do bsub.py 2 -q long $i ./$i; done

### generate gff files
- cat *gff | join_aug_pred.pl > all.gff
- grep SM_V7_ZW all.gff | sort -k1,1 -k4,4 -V > aug_ZW.sort.gff
- for i in {1..7}; do grep SM_V7_$i all.gff | sort -k1,1 -k4,4 -V > aug_C$i.sort.gff; done
- bgzip aug_C1.sort.gff
- tabix -p gff aug_C1.sort.gff.gz

### gffcompare with ratt
