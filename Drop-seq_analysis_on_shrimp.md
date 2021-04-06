## Creation of Sequence Dictionary

The first file needed is the sequence dictionary. This is a list of the contigs in the fastq file and their lengths.

```
java -jar ../Drop-seq/picard-2.23.1/picard.jar CreateSequenceDictionary \
REFERENCE=Mj_hem_tr.fasta \
OUTPUT=Mj_hem_tr.fasta.dict \
SPECIES=Mj
```

## ConvertToRefFlat

The next file is the refFlat file, which is generated using the sequence dictionary generated above.

```
../Drop-seq/Drop-seq_tools-2.3.0/ConvertToRefFlat \
ANNOTATIONS_FILE=Mj_hem_tr.gtf \
SEQUENCE_DICTIONARY=Mj_hem_tr.fasta.dict \
OUTPUT=Mj_hem_tr.refFlat
```

## fastq to bam

Convert  FASTQ files to BAM files

```
## shrimp1
java -jar ../Drop-seq/picard-2.23.1/picard.jar FastqToSam \
F1=../../../mnt/e/RawData_NGS/Mjaponicus/200706_Drop-seq/scRNA_hem1_1.fastq.gz \
F2=../../../mnt/e/RawData_NGS/Mjaponicus/200706_Drop-seq/scRNA_hem1_2.fastq.gz \
O=fqtobam1.bam \
SORT_ORDER=queryname \
SM=shrimp1
```

```
## shrimp2
java -jar ../Drop-seq/picard-2.23.1/picard.jar FastqToSam \
F1=../../../mnt/e/RawData_NGS/Mjaponicus/200706_Drop-seq/scRAN_hem2_1.fastq.gz \
F2=../../../mnt/e/RawData_NGS/Mjaponicus/200706_Drop-seq/scRAN_hem2_2.fastq.gz \
O=fqtobam2.bam \
SORT_ORDER=queryname \
SM=shrimp2
```

```
## shrimp3
java -jar ../Drop-seq/picard-2.23.1/picard.jar FastqToSam \
F1=../../../mnt/e/RawData_NGS/Mjaponicus/200706_Drop-seq/scRNA_hem3_1.fastq.gz \
F2=../../../mnt/e/RawData_NGS/Mjaponicus/200706_Drop-seq/scRNA_hem3_2.fastq.gz \
O=fqtobam3.bam \
SORT_ORDER=queryname \
SM=shrimp3
```

## Cell barcode

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
TMP_DIR=./TMP \
INPUT=fqtobam1.bam \
OUTPUT=unaligned_tagged_Cell_shrimp1.bam \
SUMMARY=unaligned_tagged_Cellular_shrimp1_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
TMP_DIR=./TMP \
INPUT=fqtobam2.bam \
OUTPUT=unaligned_tagged_Cell_shrimp2.bam \
SUMMARY=unaligned_tagged_Cellular_shrimp2_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
TMP_DIR=./TMP \
INPUT=fqtobam3.bam \
OUTPUT=unaligned_tagged_Cell_shrimp3.bam \
SUMMARY=unaligned_tagged_Cellular_shrimp3_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1
```

## Molecular Barcode

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_Cell_shrimp1.bam \
OUTPUT=unaligned_tagged_CellMolecular_shrimp1.bam \
SUMMARY=unaligned_tagged_CellMolecular_shrimp1_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_Cell_shrimp2.bam \
OUTPUT=unaligned_tagged_CellMolecular_shrimp2.bam \
SUMMARY=unaligned_tagged_CellMolecular_shrimp2_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/TagBamWithReadSequenceExtended \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_Cell_shrimp3.bam \
OUTPUT=unaligned_tagged_CellMolecular_shrimp3.bam \
SUMMARY=unaligned_tagged_CellMolecular_shrimp3_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1
```

## FilterBam

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/FilterBam \
TMP_DIR=./TMP \
TAG_REJECT=XQ \
INPUT=unaligned_tagged_CellMolecular_shrimp1.bam \
OUTPUT=unaligned_tagged_filtered_shrimp1.bam \
SUMMARY=Filter_shrimp1_summary.txt
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/FilterBam \
TMP_DIR=./TMP \
TAG_REJECT=XQ \
INPUT=unaligned_tagged_CellMolecular_shrimp2.bam \
OUTPUT=unaligned_tagged_filtered_shrimp2.bam \
SUMMARY=Filter_shrimp2_summary.txt
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/FilterBam \
TMP_DIR=./TMP \
TAG_REJECT=XQ \
INPUT=unaligned_tagged_CellMolecular_shrimp3.bam \
OUTPUT=unaligned_tagged_filtered_shrimp3.bam \
SUMMARY=Filter_shrimp3_summary.txt
```

## TrimStartingSequence

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/TrimStartingSequence \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_filtered_shrimp1.bam \
OUTPUT=unaligned_tagged_trimmed_smart_shrimp1.bam \
OUTPUT_SUMMARY=adapter_trimming_report_shrimp1.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/TrimStartingSequence \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_filtered_shrimp2.bam \
OUTPUT=unaligned_tagged_trimmed_smart_shrimp2.bam \
OUTPUT_SUMMARY=adapter_trimming_report_shrimp2.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/TrimStartingSequence \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_filtered_shrimp3.bam \
OUTPUT=unaligned_tagged_trimmed_smart_shrimp3.bam \
OUTPUT_SUMMARY=adapter_trimming_report_shrimp3.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5
```

## PolyATrimmer

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/PolyATrimmer \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_trimmed_smart_shrimp1.bam \
OUTPUT=unaligned_mc_tagged_polyA_filtered_shrimp1.bam \
OUTPUT_SUMMARY=polyA_trimming_report_shrimp1.txt \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=true
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/PolyATrimmer \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_trimmed_smart_shrimp2.bam \
OUTPUT=unaligned_mc_tagged_polyA_filtered_shrimp2.bam \
OUTPUT_SUMMARY=polyA_trimming_report_shrimp2.txt \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=true
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/PolyATrimmer \
TMP_DIR=./TMP \
INPUT=unaligned_tagged_trimmed_smart_shrimp3.bam \
OUTPUT=unaligned_mc_tagged_polyA_filtered_shrimp3.bam \
OUTPUT_SUMMARY=polyA_trimming_report_shrimp3.txt \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=true
```

## SamToFastq

```
## shrimp1
java -Xmx10g -jar ../Drop-seq/picard-2.23.1/picard.jar SamToFastq \
INPUT=unaligned_mc_tagged_polyA_filtered_shrimp1.bam \
FASTQ=unaligned_mc_tagged_polyA_filtered_shrimp1.fastq
```

```
## shrimp2
java -Xmx10g -jar ../Drop-seq/picard-2.23.1/picard.jar SamToFastq \
INPUT=unaligned_mc_tagged_polyA_filtered_shrimp2.bam \
FASTQ=unaligned_mc_tagged_polyA_filtered_shrimp2.fastq
```

```
## shrimp3
java -Xmx10g -jar ../Drop-seq/picard-2.23.1/picard.jar SamToFastq \
INPUT=unaligned_mc_tagged_polyA_filtered_shrimp3.bam \
FASTQ=unaligned_mc_tagged_polyA_filtered_shrimp3.fastq
```

## Trim Galore
[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) var. 0.6.6

Quality trimming

```
## shrimp1
trim_galore unaligned_mc_tagged_polyA_filtered_shrimp1.fastq
```

```
## shrimp2
trim_galore unaligned_mc_tagged_polyA_filtered_shrimp2.fastq
```

```
## shrimp3
trim_galore unaligned_mc_tagged_polyA_filtered_shrimp3.fastq
```

## STAR
[STAR](https://github.com/alexdobin/STAR) [ref](https://academic.oup.com/bioinformatics/article/29/1/15/272537) var. 2.7.8a

## Construction of a database

```
## construct a database
STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir Mj_hem_tr.fasta \
--genomeFastaFiles ../metadata/Mj_hem_tr.fasta \
--sjdbGTFfile ../metadata/Mj_hem_tr.gtf \
--sjdbOverhang 100 \
--genomeSAindexNbases 11 \
--genomeChrBinNbits 18
```

## mapping

```
## shrimp1
STAR \
--runThreadN 10 \
--genomeDir STAR/Mj_hem_tr.fasta \
--readFilesIn unaligned_mc_tagged_polyA_filtered_shrimp1_trimmed.fq \
--genomeLoad NoSharedMemory \
--outFilterScoreMinOverLread 0.5 \
--outFilterMatchNminOverLread 0.5 \
--outFilterMultimapNmax 20 \
--outFileNamePrefix STAR/shrimp1
```

```
## shrimp2
STAR \
--runThreadN 10 \
--genomeDir STAR/Mj_hem_tr.fasta \
--readFilesIn unaligned_mc_tagged_polyA_filtered_shrimp2_trimmed.fq \
--genomeLoad NoSharedMemory \
--outFilterScoreMinOverLread 0.5 \
--outFilterMatchNminOverLread 0.5 \
--outFilterMultimapNmax 20 \
--outFileNamePrefix STAR/shrimp2
```

```
## shrimp3
STAR \
--runThreadN 10 \
--genomeDir STAR/Mj_hem_tr.fasta \
--readFilesIn unaligned_mc_tagged_polyA_filtered_shrimp3_trimmed.fq \
--genomeLoad NoSharedMemory \
--outFilterScoreMinOverLread 0.5 \
--outFilterMatchNminOverLread 0.5 \
--outFilterMultimapNmax 20 \
--outFileNamePrefix STAR/shrimp3
```

## samtobam/bam sort

```
## shrimp1
java -jar ../Drop-seq/picard-2.23.1/picard.jar \
SortSam \
I=STAR/shrimp1Aligned.out.sam \
O=aligned.sorted_shrimp1.bam \
SO=queryname
```

```
## shrimp2
java -jar ../Drop-seq/picard-2.23.1/picard.jar \
SortSam \
I=STAR/shrimp2Aligned.out.sam \
O=aligned.sorted_shrimp2.bam \
SO=queryname
```

```
## shrimp3
java -jar ../Drop-seq/picard-2.23.1/picard.jar \
SortSam \
I=STAR/shrimp3Aligned.out.sam \
O=aligned.sorted_shrimp3.bam \
SO=queryname
```

## MergeBamAlignment

```
## shrimp1
java -jar -Xms8000m ../Drop-seq/picard-2.23.1/picard.jar \
MergeBamAlignment \
REFERENCE_SEQUENCE=metadata/Mj_hem_tr.fasta \
UNMAPPED_BAM=unaligned_mc_tagged_polyA_filtered_shrimp1.bam \
ALIGNED_BAM=aligned.sorted_shrimp1.bam \
OUTPUT=merged_shrimp1.bam \
INCLUDE_SECONDARY_ALIGNMENTS=TRUE \
PAIRED_RUN=false
```

```
## shrimp2
java -jar -Xms8000m ../Drop-seq/picard-2.23.1/picard.jar \
MergeBamAlignment \
REFERENCE_SEQUENCE=metadata/Mj_hem_tr.fasta \
UNMAPPED_BAM=unaligned_mc_tagged_polyA_filtered_shrimp2.bam \
ALIGNED_BAM=aligned.sorted_shrimp2.bam \
OUTPUT=merged_shrimp2.bam \
INCLUDE_SECONDARY_ALIGNMENTS=TRUE \
PAIRED_RUN=false
```

```
## shrimp3
java -jar -Xms8000m ../Drop-seq/picard-2.23.1/picard.jar \
MergeBamAlignment \
REFERENCE_SEQUENCE=metadata/Mj_hem_tr.fasta \
UNMAPPED_BAM=unaligned_mc_tagged_polyA_filtered_shrimp3.bam \
ALIGNED_BAM=aligned.sorted_shrimp3.bam \
OUTPUT=merged_shrimp3.bam \
INCLUDE_SECONDARY_ALIGNMENTS=TRUE \
PAIRED_RUN=false
```

## TagReadWithGeneFunction

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/TagReadWithGeneFunction \
TMP_DIR=./TMP \
I=merged_shrimp1.bam \
O=gene_exon_tagged_shrimp1.bam \
ANNOTATIONS_FILE=metadata/Mj_hem_tr.refFlat
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/TagReadWithGeneFunction \
TMP_DIR=./TMP \
I=merged_shrimp2.bam \
O=gene_exon_tagged_shrimp2.bam \
ANNOTATIONS_FILE=metadata/Mj_hem_tr.refFlat
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/TagReadWithGeneFunction \
TMP_DIR=./TMP \
I=merged_shrimp3.bam \
O=gene_exon_tagged_shrimp3.bam \
ANNOTATIONS_FILE=metadata/Mj_hem_tr.refFlat
```

## DetectBeadSubstitutionErrors - Detecting and repairing substitution errors in cell barcodes

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/DetectBeadSubstitutionErrors \
TMP_DIR=./TMP \
I=gene_exon_tagged_shrimp1.bam \
O=clean_substitution_shrimp1.bam \
OUTPUT_REPORT=clean.substitution_report_shrimp1.txt
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/DetectBeadSubstitutionErrors \
TMP_DIR=./TMP \
I=gene_exon_tagged_shrimp2.bam \
O=clean_substitution_shrimp2.bam \
OUTPUT_REPORT=clean.substitution_report_shrimp2.txt
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/DetectBeadSubstitutionErrors \
TMP_DIR=./TMP \
I=gene_exon_tagged_shrimp3.bam \
O=clean_substitution_shrimp3.bam \
OUTPUT_REPORT=clean.substitution_report_shrimp3.txt
```

## DetectBeadSynthesisErrors

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/DetectBeadSynthesisErrors \
TMP_DIR=./TMP \
I=clean_substitution_shrimp1.bam \
O=clean_shrimp1.bam \
REPORT=clean_report_shrimp1.txt \
OUTPUT_STATS=synthesis_stats_shrimp1.txt \
SUMMARY=synthesis_stats_summary_shrimp1.txt \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/DetectBeadSynthesisErrors \
TMP_DIR=./TMP \
I=clean_substitution_shrimp2.bam \
O=clean_shrimp2.bam \
REPORT=clean_report_shrimp2.txt \
OUTPUT_STATS=synthesis_stats_shrimp2.txt \
SUMMARY=synthesis_stats_summary_shrimp2.txt \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/DetectBeadSynthesisErrors \
TMP_DIR=./TMP \
I=clean_substitution_shrimp3.bam \
O=clean_shrimp3.bam \
REPORT=clean_report_shrimp3.txt \
OUTPUT_STATS=synthesis_stats_shrimp3.txt \
SUMMARY=synthesis_stats_summary_shrimp3.txt \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
```

## Cell Selection

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/BAMTagHistogram \
TMP_DIR=./TMP \
I=clean_shrimp1.bam \
O=out_cell_readcounts_shrimp1.txt.gz \
READ_MQ=0 \
TAG=XC
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/BAMTagHistogram \
TMP_DIR=./TMP \
I=clean_shrimp2.bam \
O=out_cell_readcounts_shrimp2.txt.gz \
READ_MQ=0 \
TAG=XC
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/BAMTagHistogram \
TMP_DIR=./TMP \
I=clean_shrimp3.bam \
O=out_cell_readcounts_shrimp3.txt.gz \
READ_MQ=0 \
TAG=XC
```

## Extract cells

```
## shrimp1
../Drop-seq/Drop-seq_tools-2.3.0/DigitalExpression \
TMP_DIR=./TMP \
I=clean_shrimp1.bam \
O=shrimp1.txt.gz \
SUMMARY=shrimp1_summary.txt \
NUM_CORE_BARCODES=3000 \
READ_MQ=0
```

```
## shrimp2
../Drop-seq/Drop-seq_tools-2.3.0/DigitalExpression \
TMP_DIR=./TMP \
I=clean_shrimp2.bam \
O=shrimp2.txt.gz \
SUMMARY=shrimp2_summary.txt \
NUM_CORE_BARCODES=3000 \
READ_MQ=0
```

```
## shrimp3
../Drop-seq/Drop-seq_tools-2.3.0/DigitalExpression \
TMP_DIR=./TMP \
I=clean_shrimp3.bam \
O=shrimp3.txt.gz \
SUMMARY=shrimp3_summary.txt \
NUM_CORE_BARCODES=3000 \
READ_MQ=0
```
