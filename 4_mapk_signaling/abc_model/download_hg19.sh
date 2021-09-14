# Download hg19 from UCSC and prepare bwa indexes

# Go to working dir
mkdir hg19
cd hg19

wget 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz’ 
wget 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes’ 
gunzip hg19.fa.gz 

# Make a hg19.main.fa file
srun --pty bash
module load samtools
samtools faidx hg19.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > hg19.main.fa

# And a bwa index
module load bwa/0.7.15 
bwa index hg19.main.fa hg19.main

exit
