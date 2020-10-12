## From Sushi.R manual
Other popular file formats such as BAM and GFF are not explicitly supported by Sushi. However, data stored in these formats can be easily converted to BED format using common command line tools such as the bedtools software suite available at https://github.com/arq5x/bedtools2. Some examples taken from the bedtools are shown below.
Convert BAM alignments to BED format.
bamToBed -i reads.bam > reads.bed
Convert BAM alignments to BED format using edit distance (NM) as the BED score.
bamToBed -i reads.bam -ed > reads.bed
Convert BAM alignments to BEDPE format.
bamToBed -i reads.bam -bedpe > reads.bedpe
These BED files can easily be read into R for use with Sushi using the following R command:
> read.table(file="reads.bed",sep="\t")

## bam to bed
bamToBed -i /lustre/scratch118/infgen/team133/zl3/smansoni/alt_v7/IsoSeq_bams/Sman_iso_som_V7.sorted.bam -ed > iso_som.bed
