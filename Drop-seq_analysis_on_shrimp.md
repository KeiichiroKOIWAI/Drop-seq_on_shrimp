## CreateSequenceDictionary
### The first file needed is the sequence dictionary. This is a list of the contigs in the fastq file and their lengths.
java -jar ../picard-2.23.1/picard.jar CreateSequenceDictionary \
REFERENCE=Mj_hem_tr.fasta \
OUTPUT=Mj_hem_tr.fasta.dict \
SPECIES=Mj

## ConvertToRefFlat
### The next file is the refFlat file, which is generated using the sequence dictionary generated above.
../Drop-seq_tools-2.3.0/ConvertToRefFlat \
ANNOTATIONS_FILE=Mj_evigene.gtf \
SEQUENCE_DICTIONARY=Mj_evigene.dict \
OUTPUT=Mj_evigene.refFlat

## fastq to bam
java -jar ../../picard-2.23.1/picard.jar FastqToSam \
F1=scRNA_hemX_1.fastq.gz \
F2=scRNA_hemX_2.fastq.gz \
O=fqtobamX.bam \
SM=shrimpX

##bam sort
java -jar ../../picard-2.23.1/picard.jar SortSam \
I=fqtobamX.bam \
O=sortedX.bam \
SORT_ORDER=queryname

# Drop-seq tool
## Cell barcode
../../Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
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
../../Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
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
../../Drop-seq_tools-2.3.0/FilterBam \
TMP_DIR=./TMP \
TAG_REJECT=XQ \
INPUT=unaligned_tagged_CellMolecular.bam \
OUTPUT=unaligned_tagged_filtered.bam \
SUMMARY=Filter.summary.txt

## TrimStartingSequence
../../Drop-seq_tools-2.3.0/TrimStartingSequence \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_filtered.bam \
OUTPUT=unaligned_tagged_trimmed_smart.bam \
OUTPUT_SUMMARY=adapter_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5

## PolyATrimmer
../../Drop-seq_tools-2.3.0/PolyATrimmer \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_trimmed_smart.bam \
OUTPUT=unaligned_mc_tagged_polyA_filtered.bam \
OUTPUT_SUMMARY=polyA_trimming_report.txt \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=true

## SamToFastq
java -jar ../../picard-2.23.1/picard.jar SamToFastq \
INPUT=unaligned_mc_tagged_polyA_filtered.bam \
FASTQ=unaligned_mc_tagged_polyA_filtered.fastq

## seqkit stats
seqkit stat -j 6 unaligned_mc_tagged_polyA_filtered.fastq

## Remove sewuence less than 20bp
cat unaligned_mc_tagged_polyA_filtered.fastq | seqkit seq -j 12 -m 20 > unaligned_mc_tagged_polyA_20bpfiltered.fastq

## seqkit stats
seqkit stat -j 6 unaligned_mc_tagged_polyA_20bpfiltered.fastq

## Bowtie2 index
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) [ref](https://www.nature.com/articles/nmeth.1923) var. 2.4.1  
bowtie2-build Mj_hem_tr.fasta Mj_hem

## Bowtie2 mapping
bowtie2 \
-p 8 \
-q \
-N 0 \
--local \
-x ../bt2index/Mj_hem \
-U ../Mj_hem/shrimp2/unaligned_mc_tagged_polyA_20bpfiltered.fastq \
-S bt_aligned_evigene2.sam \
-t

## SAM validation
java -jar -Xms500m ../picard-2.23.1/picard.jar ValidateSamFile \
I=bt_aligned_evigene1.sam \
MODE=SUMMARY

## Remove extra raws
sed '30988d' bt_aligned_evigene1.sam > bt_aligned_evigene1_rm.sam

## samtobam/bam sort
java -jar -Xms500m ../picard-2.23.1/picard.jar SortSam \
I=bt_aligned_evigene1_rm.sam \
O=bt_aligned.sorted1.bam \
SO=queryname

## MergeBamAlignment
java -jar -Xms8000m ../picard-2.23.1/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=Mj_evigene.fasta \
UNMAPPED_BAM=../Mj_hem/shrimp1/unaligned_mc_tagged_polyA_filtered.bam \
ALIGNED_BAM=bt_aligned.sorted1.bam \
OUTPUT=bt_merged1.bam \
INCLUDE_SECONDARY_ALIGNMENTS=TRUE \
PAIRED_RUN=false

## TagReadWithGeneFunction
../Drop-seq_tools-2.3.0/TagReadWithGeneFunction \
TMP_DIR=./TMP \
I=bt_merged1.bam \
O=gene_exon_tagged1.bam \
ANNOTATIONS_FILE=Mj_evigene.refFlat

## DetectBeadSubstitutionErrors - Detecting and repairing substitution errors in cell barcodes
../Drop-seq_tools-2.3.0/DetectBeadSubstitutionErrors \
TMP_DIR=./TMP \
I=gene_exon_tagged1.bam \
O=bt_my_clean_subtitution1.bam \
OUTPUT_REPORT=bt_my_clean.substitution_report1.txt

## DetectBeadSynthesisErrors
../Drop-seq_tools-2.3.0/DetectBeadSynthesisErrors \
TMP_DIR=./TMP \
I=bt_my_clean_subtitution1.bam \
O=bt_my_clean1.bam \
REPORT=bt_my_clean.indel_report1.txt \
OUTPUT_STATS=bt_my_synthesis_stats1.txt \
SUMMARY=bt_my_synthesis_stats.summary1.txt \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

## Cell Selection
../Drop-seq_tools-2.3.0/BAMTagHistogram \
TMP_DIR=./TMP \
I=bt_my_clean1.bam \
O=out_cell_readcounts1.txt.gz \
TAG=XC

## On the R software
a=read.table("out_cell_readcounts3.txt.gz", header=F, stringsAsFactors=F) \
x=cumsum(a$V1) \
x=x/max(x) \
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", \
ylab="cumulative fraction of reads", xlim=c(1,3000))

## Extract cells
../Drop-seq_tools-2.3.0/DigitalExpression \
TMP_DIR=./TMP \
I=bt_my_clean3.bam \
O=1500.dge.txt.gz \
SUMMARY=1500.dge2.summary.txt \
NUM_CORE_BARCODES=1500
