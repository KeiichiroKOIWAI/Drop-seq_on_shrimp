## Get short read sequence of shrimp hemocyte
In this study we used this raw sequences from kuruma shrimop hemocytes  
DDBJ Sequence Read Archive (DRA) accession number [DRR004781](https://ddbj.nig.ac.jp/DRASearch/submission?acc=DRA004781)  

```
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004781/DRX056819/DRR062578_1.fastq.bz2
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004781/DRX056819/DRR062578_2.fastq.bz2  
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004781/DRX056819/DRR062579_1.fastq.bz2  
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004781/DRX056819/DRR062579_2.fastq.bz2  
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004781/DRX056819/DRR062580_1.fastq.bz2  
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004781/DRX056819/DRR062580_2.fastq.bz2`
```

tar all files  
```
cat DRR062578_1.fastq DRR062579_1.fastq DRR062580_1.fastq > PRJDB4895_L.fastq
cat DRR062578_2.fastq DRR062579_2.fastq DRR062580_2.fastq > PRJDB4895_R.fastq
```

## Counting kmer
[jellyfish software](http://www.genome.umd.edu/jellyfish.html#Release) var. 2.3.0  
```
jellyfish count \
-m 21 \
-t 12 \
-s 10000000 \
-o jellyfish21 \
PRJDB4895_L.fastq PRJDB4895_R.fastq  
```

## jellyfish dump
```
jellyfish dump -c jellyfish21 > jellyfish21.dump
```

## Long read sequencing correction by TALC
[TALC](https://gitlab.igh.cnrs.fr/lbroseus/TALC) [ref](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btaa634/5872522?redirectedFrom=fulltext) var. 1.01  
```
talc "ONT_basecalled_sequence.fastq" \
--SRCounts  jellyfish21.dump \
-k 21 \
-o talc \
-t 10  
```

## Hybrid assenbly ONT/miseq by rnaSPAdes
[rnaSPAdes](https://cab.spbu.ru/software/rnaspades/) [ref](https://cab.spbu.ru/software/rnaspades/) var. 3.14.1  
```
rnaspades.py \
--threads 12 \
--memory 28 \
--pe1-1 PRJDB4895_L.fastq.gz \
--pe1-2 PRJDB4895_R.fastq.gz \
--ss rf \
--fl-rna talc/talc.fa \
-o rnaSPAdes
```

## Seqkit stats
[seqkit](https://bioinf.shenwei.me/seqkit/)  
```
seqkit stats -a rnaSPAdes/transcripts.fasta
```

## launch docker Trinity
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) [ref](https://www.nature.com/articles/nbt.1883) var. 2.10.0  
```
docker run -v $mnt/h/tri:/usr/local/data -it trinityrnaseq/trinityrnaseq
```

## Hybrid assembly ONT/miseq by Trinity
```
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
```

## Seqkit stats
```
seqkit stats -a trinity/Trinity.fasta
```

## Trinity stats
```
bin/trinityrnaseq/util/TrinityStats.pl \
data/trinity/Trinity.fasta
```
## cat all fasta
```
cat rnaSPAdes/transcripts.fasta trinity/Trinity.fasta > Mj_hem_raw.fasta
```

## Remove duplicates and find coding contigs by EvidentialGene
[EvidentialGene] (http://arthropods.eugenes.org/EvidentialGene/) var. 2022.01.20  
```
perl ../Evigene/scripts/prot/tr2aacds4.pl \
-cdnaseq \
Mj_hem_raw.fasta \
-NCPU=10  -MAXMEM=10000 \
-MINCDS=90 \
-logfile
```

## Seqkit stats
```
seqkit stats -a okayset/Mj_hem_raw.okay.mrna
```

## Renaming of headers
[seqtk](https://github.com/lh3/seqtk)  
```
seqtk seq -C /Mj_hem_raw.okay.mrna > seqs.fa.tmp
seqtk rename seqs.fa.tmp Mj- > Mj_hem_tr.fasta
```

## BUSCO
[BUSCO](https://busco.ezlab.org/) [ref](https://link.springer.com/protocol/10.1007%2F978-1-4939-9173-0_14) var.5.0.0
Checking the quality of assembled transcripts.
```
## BUSCO
busco \
-m transcriptome \
-i metadata/Mj_hem_tr.fasta \
-o BUSCO_result \
-l arthropoda_odb10 \
-c 4 \
-f
```
