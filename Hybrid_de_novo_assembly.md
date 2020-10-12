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

## Hybrid assenbly ONT/miseq by rnaSPAdes
[rnaSPAdes](https://cab.spbu.ru/software/rnaspades/)  
rnaspades.py \
--threads 12 \
--memory 28 \
--pe1-1 PRJDB4895_L.fastq.gz \
--pe1-2 PRJDB4895_R.fastq.gz \
--ss rf \
--fl-rna talc/talc.fa \
-o rnaSPAdes

## Seqkit stats
[seqkit](https://bioinf.shenwei.me/seqkit/)  
seqkit stats -a rnaSPAdes/transcripts.fasta

## launch docker Trinity
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)  
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

## Trinity stats
bin/trinityrnaseq/util/TrinityStats.pl \
data/trinity/Trinity.fasta

## cat all fasta
cat rnaSPAdes/transcripts.fasta trinity/Trinity.fasta > Mj_hem_raw.fasta

## Remove duplicates and find coding contigs by EvidentialGene
[EvidentialGene] (http://arthropods.eugenes.org/EvidentialGene/)  
perl ../Evigene/scripts/prot/tr2aacds4.pl \
-cdnaseq \
Mj_hem_raw.fasta \
-NCPU=10  -MAXMEM=10000 \
-MINCDS=90 \
-logfile

seqkit stats -a okayset/Mj_hem_raw.okay.mrna

## Renaming of headers
[seqtk](https://github.com/lh3/seqtk)  
seqtk seq -C /Mj_hem_raw.okay.mrna > seqs.fa.tmp
seqtk rename seqs.fa.tmp Mj- > Mj_hem_tr.fasta
