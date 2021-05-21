# Drop-seq_on_shrimp
These are commads for Drop-seq analyzing of hemocytes of kuruma shrimp (*Marsupenaues japonicus*). The detailed scientific explanation is published under BioRxiv.  
Here, I uploaded commads used on this paper.

The data have been deposited with links to BioProject accession number PRJDB10658 in the DDBJ BioProject database.

## Basecalling of ONT RNA direct sequencing
[command link](https://github.com/KeiichiroKOIWAI/Drop-seq_on_shrimp/blob/main/linux_commandlines/ONT_basecalling_using_guppy.md)  

The raw sequence data of newly sequenced M. japonicus transcriptomic reads were archived in DDBJ Sequence Read Archive (DRA) of the DNA Data Bank of Japan : MinION mRNA direct sequencing: DRA010948. 

## Hybrid assembly of ONT long seqence and illumina short read sequence
[command link](https://github.com/KeiichiroKOIWAI/Drop-seq_on_shrimp/blob/main/linux_commandlines/Hybrid_de_novo_assembly.md)  

Two raw sequence data were used in this assembly, long reads:DRA010948 and short reads: DRA004781. 

The assembled transcripts were deposited [here](https://github.com/KeiichiroKOIWAI/Drop-seq_on_shrimp/blob/main/metadata.tar.gz).

## Command lines for Drop-seq analysis of shrimp hemocytes
[command link](https://github.com/KeiichiroKOIWAI/Drop-seq_on_shrimp/blob/main/linux_commandlines/Drop-seq_analysis_on_shrimp.md)   

The raw sequence data used in this study were archived in DDBJ Sequence Read Archive (DRA) of the DNA Data Bank of Japan: Drop-seq shrimp rep1: DRA010950, shrimp rep2: DRA010951, shrimp rep3: DRA010952. 

The obtained data were archived in Genomic Expression Archive of the DNA Data Bank of Japan: E-GEAD-403.

## Command lines for Seurat and monocle3 of shrimp hemocytes
Until UMAP reductuion and clustering by Seurat. [command link](https://github.com/KeiichiroKOIWAI/Drop-seq_on_shrimp/blob/main/R_commandlines/Seurat.md)  
Until pseudotime analysis by Seurat and Monocle3. [command link](https://github.com/KeiichiroKOIWAI/Drop-seq_on_shrimp/blob/main/R_commandlines/monocle.md)

## Command lines for Drop-seq analysis of human/mouse experiment
[command link](https://github.com/KeiichiroKOIWAI/Drop-seq_on_shrimp/blob/main/linux_commandlines/Human_Mouse_species_mix.md)  

The raw sequence data of newly sequenced  human/mouse cell transcriptomic reads were archived in DDBJ Sequence Read Archive (DRA) of the DNA Data Bank of Japan: mixture sample of HEK293 and 3T3 cells: DRA010949. 
