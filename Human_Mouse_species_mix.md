## fastq to bam
### fastqファイルをbamファイルに変換する
java -jar picard.jar FastqToSam \
#後日DDBJのURLを記入する
F1=3p133Lib2_S1_L001_R1_001.fastq \
F2=3p133Lib2_S1_L001_R2_001.fastq \
O=fqtobam.bam \
SM=HM

## bam sort
### 並びを整える
java -jar picard.jar SortSam \
I=fqtobam.bam \
O=sorted.bam \
SORT_ORDER=queryname

## Cell barcode using Drop-seq_tools-2.3.0
### bamファイルにcellバーコードのタグをくっつける
TagBamWithReadSequenceExtended \
TMP_DIR=./TMP \
INPUT=sorted.bam \
OUTPUT=unaligned_tagged_Cell.bam \
SUMMARY=unaligned_tagged_Cellular.bam_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1

## Molecular barcode
### bamファイルにmolecularバーコードのタグをくっつけ
TagBamWithReadSequenceExtended \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_Cell.bam \
OUTPUT=unaligned_tagged_CellMolecular.bam \
SUMMARY=unaligned_tagged_Molecular.bam_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1

## FilterBam:
FilterBam \
TMP_DIR=./TMP \
-m 10g \
TAG_REJECT=XQ \
INPUT=unaligned_tagged_CellMolecular.bam \
OUTPUT=unaligned_tagged_filtered.bam

## TrimStartingSequence
TrimStartingSequence \
TMP_DIR=./TMP \
-m 10g \
INPUT=unaligned_tagged_filtered.bam \
OUTPUT=unaligned_tagged_trimmed_smart.bam \
OUTPUT_SUMMARY=adapter_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5

## PolyATrimmer
PolyATrimmer \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_trimmed_smart.bam \
OUTPUT=unaligned_mc_tagged_polyA_filtered.bam \
OUTPUT_SUMMARY=polyA_trimming_report.txt \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=true

## SamToFastq
java -jar ../picard-2.23.1/picard.jar SamToFastq \
INPUT=unaligned_mc_tagged_polyA_filtered.bam \
FASTQ=unaligned_mc_tagged_polyA_filtered.fastq

## Removing sequences less than 20 bp
### 20bp以下の配列を削除
cat unaligned_mc_tagged_polyA_filtered.fastq | seqkit seq -m 20 > unaligned_mc_tagged_polyA_20bpfiltered.fastq

## Buiding database on Hisat2
### human and mouse transgenes reference metadata was obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63269
hisat2-build hg19_mm10_transgenes.fasta hg19_mm10

## Hisat alignment
hisat2 \
-p 12 \
-q \
-x hg19_mm10 \
-U unaligned_mc_tagged_polyA_20bpfiltered.fastq \
-S aligned.sam \
-t

## SAM validation
java -Xmx4g -jar ../picard-2.23.1/picard.jar ValidateSamFile I=aligned2.sam MODE=SUMMARY

## Removing rows that include information of Hisat2
sed '256d' aligned.sam > aligned2.sam

## samtobam/bam sort
java -jar picard.jar SortSam \
I=aligned2.sam \
O=aligned.sorted.bam \
SO=queryname

## MergeBamAlignment
java -jar picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=hg19_mm10_transgenes.fasta \
UNMAPPED_BAM=unaligned_mc_tagged_polyA_filtered.bam \
ALIGNED_BAM=aligned.sorted.bam \
OUTPUT=merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false

## TagReadWithGeneFunction
TagReadWithGeneFunction \
TMP_DIR=./TMP \
I=merged.bam \
O=hisat2_gene_exon_tagged.bam \
ANNOTATIONS_FILE=hg19_mm10_transgenes.refFlat

## DetectBeadSubstitutionErrors - Detecting and repairing substitution errors in cell barcodes
DetectBeadSubstitutionErrors \
TMP_DIR=./TMP \
I=hisat2_gene_exon_tagged.bam \
O=my_clean_subtitution.bam \
OUTPUT_REPORT=my_clean.substitution_report.txt

## DetectBeadSynthesisErrors
DetectBeadSynthesisErrors \
TMP_DIR=./TMP \
I=my_clean_subtitution.bam \
O=my_clean.bam \
REPORT=my_clean.indel_report.txt \
OUTPUT_STATS=my.synthesis_stats.txt \
SUMMARY=my.synthesis_stats.summary.txt \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

## Cell Selection
BAMTagHistogram \
TMP_DIR=./TMP \
I=my_clean.bam \
O=out_cell_readcounts.txt.gz \
TAG=XC

## on the R
a=read.table("out_cell_readcounts.txt.gz", header=F, stringsAsFactors=F) \
x=cumsum(a$V1) \
x=x/max(x) \
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", \
ylab="cumulative fraction of reads", xlim=c(1,1500))

## DigitalExpression of total cells
DigitalExpression \
TMP_DIR=./TMP \
I=my_clean.bam \
O=500.dge.txt.gz \
SUMMARY=500.dge.summary.txt \
NUM_CORE_BARCODES=500
#後日GEOをいれる

## FilterBam into Human/Mouse
### For human
FilterBam \
TMP_DIR=./TMP \
I=my_clean.bam \
O=my_clean_HUMAN.bam \
REF_SOFT_MATCHED_RETAINED=HUMAN

### For mouse
FilterBam \
TMP_DIR=./TMP \
I=my_clean.bam \
O=my_clean_MOUSE.bam \
REF_SOFT_MATCHED_RETAINED=MOUSE

## DigitalExpression of each specie
### HUMAN
DigitalExpression \
TMP_DIR=./TMP \
I=my_clean_HUMAN.bam \
O=HUMAN.dge.txt.gz \
SUMMARY=HUMAN.dge.summary.txt \
CELL_BC_FILE=500list.txt
#後日GEOをいれる

### MOUSE
DigitalExpression \
TMP_DIR=./TMP \
I=my_clean_MOUSE.bam \
O=MOUSE.dge.txt.gz \
SUMMARY=MOUSE.dge.summary.txt \
CELL_BC_FILE=500list.txt
#後日GEOをいれる

## SP-mix plotting on the R
human=read.table("HUMAN.dge.summary.txt", header=T, stringsAsFactors=F) \
mouse=read.table("MOUSE.dge.summary.txt", header=T, stringsAsFactors=F) \
human.bc=human[order(human$CELL_BARCODE),] \
mouse.bc=mouse[order(mouse$CELL_BARCODE),] \
spmix=cbind(human.bc,mouse.bc) \
colnames(spmix)=c("A","B","C","D","E","F","G","H") \
spmix=transform(spmix, I=spmix$C/(spmix$C+spmix$G)) \
spmix=transform(spmix, J=spmix$G/(spmix$C+spmix$G)) \
human.plot=spmix[spmix$I>=0.9,] \
mouse.plot=spmix[spmix$J>=0.9,] \
mix.plot=spmix[spmix$I<0.9&spmix$J<0.9,]

## Plotting a figure
library(ggplot2) \
library(extrafont) \
extrafont::loadfonts(device = "win") \
#https://rstudio-pubs-static.s3.amazonaws.com/459388_c0f0abce1612426c8aec91b678da8670.html

fig=ggplot(NULL)+ \
    geom_point(aes(mix.plot$C,mix.plot$G), col="#E41A1C", alpha=0.5, size=4, shape=17)+ \
    geom_point(aes(human.plot$C,human.plot$G), col="#377EB8", alpha=0.5, size=4, shape=19)+ \
    geom_point(aes(mouse.plot$C,mouse.plot$G), col="#4DAF4A", alpha=0.5, size=4, shape=15)+ \
    ylim(0,10000)+xlim(0,25000)+ \
    xlab("human transcripts (UMIs)") + ylab("mouse transcripts (UMIs)")+ \
    theme_classic()+ \
    theme(text=element_text(size=18, family="Arial"))+ \
    theme(axis.ticks=element_line(colour = "black"),　axis.text=element_text(colour = "black")) \
fig
