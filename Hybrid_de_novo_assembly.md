## Get short read sequence of shrimp hemocyte
### URL記入PRJDB4895
unpigz -c PRJDB4895_L.fastq.gz > PRJDB4895_L.fastq  
unpigz -c PRJDB4895_R.fastq.gz > PRJDB4895_R.fastq  

## Counting kmer for TALC
[jellyfish software](http://www.genome.umd.edu/jellyfish.html#Release)  
jellyfish count -m 21 -t 12 -s 10000000 -o jellyfish21 PRJDB4895_L.fastq PRJDB4895_R.fastq  

## jellyfish dump
jellyfish dump -c jellyfish21 > jellyfish21.dump

## Long read sequencing correction by TALC
[TALC software](https://gitlab.igh.cnrs.fr/lbroseus/TALC) [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.01.10.901728v3)  
talc "ONT_basecalled_sequence.fastq" --SRCounts  jellyfish21.dump -k 21 -o talc -t 10  

## Hybrid assenbly ONT/miseq by RNAspades
rnaspades.py \
--threads 12 \
--memory 28 \
--pe1-1 PRJDB4895_L.fastq.gz \
--pe1-2 PRJDB4895_R.fastq.gz \
--ss rf \
--fl-rna talc/talc.fa \
-o rnaSPAdes

## Seqkit stats
seqkit stats -a rnaSPAdes/transcripts.fasta
file                          format  type  num_seqs     sum_len  min_len  avg_len  max_len   Q1   Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)
rnaSPAdes/transcripts.fasta  FASTA   DNA     39,757  44,545,567       92  1,120.4   29,031  296  487  1,342        0  2,320       0       0

## launch docker Trinity
docker run -v $mnt/h/tri:/usr/local/data -it trinityrnaseq/trinityrnaseq

## Hybrid assembly ONT/miseq by Trinity
Trinity \
--seqType fq \
--max_memory 20G \
--left PRJDB4895_L.fastq.gz \
--right PRJDB4895_R.fastq.gz \
--CPU 10 \
--min_contig_length 200 \
--SS_lib_type RF \
--long_reads talc/talc.fa \
--output trinity

## Seqkit stats
seqkit stats -a trinity/Trinity.fasta
file                    format  type  num_seqs     sum_len  min_len  avg_len  max_len   Q1   Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)
trinity/Trinity.fasta  FASTA   DNA     69,126  77,715,789      185  1,124.3   20,174  300  533  1,391        0  2,204       0       0

## Trinity stats
bin/trinityrnaseq/util/TrinityStats.pl \
data/trinity/Trinity.fasta

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  44904
Total trinity transcripts:      69126
Percent GC: 42.76

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 6034
        Contig N20: 4458
        Contig N30: 3500
        Contig N40: 2806
        Contig N50: 2204

        Median contig length: 533
        Average contig: 1124.26
        Total assembled bases: 77715789


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 5491
        Contig N20: 3891
        Contig N30: 2942
        Contig N40: 2220
        Contig N50: 1643

        Median contig length: 387
        Average contig: 831.48
        Total assembled bases: 37336692

## cat all fasta
cat rnaSPAdes/transcripts.fasta trinity/Trinity.fasta > Mj_hem_raw.fasta

## Removing XXXXXXXXXXX by XXXXXXXXXX
perl ../Evigene/scripts/prot/tr2aacds4.pl \
-cdnaseq \
Mj_hem_raw.fasta \
-NCPU=10  -MAXMEM=10000 \
-MINCDS=90 \
-logfile

seqkit stats -a okayset/Mj_hem_raw.okay.mrna
#file                          format  type  num_seqs     sum_len  min_len  avg_len  max_len   Q1     Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)
#okayset/Mj_hem_raw.okay.mrna  FASTA   DNA     30,986  59,567,966      232  1,922.4   29,031  659  1,371  2,530        0  2,902       0       0

## Renaming of headers
seqtk seq -C /Mj_hem_raw.okay.mrna > seqs.fa.tmp
seqtk rename seqs.fa.tmp Mj- > Mj_hem_tr.fasta