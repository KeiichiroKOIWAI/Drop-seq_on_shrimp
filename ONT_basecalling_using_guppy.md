## Get guppy software
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.6.1_linux64.tar.gz

## Untar files
tar xzvf ont-guppy_3.6.1_linux64.tar.gz

## Basecalling
/content/ont-guppy/bin/guppy_basecaller \
--input_path /"path_to_fast5_files" \
--flowcell FLO-MIN106 --kit SQK-RNA002 \
--save_path /"path_to_save_folder" \
--gpu_runners_per_device 20 \
-x auto \
--u_substitution on

## cat all fastaq files

## remove duplicates
[pyBioTools](https://github.com/a-slide/pyBioTools)  
python3 \
from pyBioTools import Fastq \
from pyBioTools.common import jhelp \
Fastq.Filter ("./fastq", "./ONT_fastq.fastq", remove_duplicates=True) \
exit()

## U>T
cat ONT_fastq.fastq | seqkit replace -p 'U' -r 'T' -s > ONT_final.fastq

## NanoPlot
[NanoPlot](https://github.com/wdecoster/NanoPlot)  
NanoPlot --fastq ONT_final.fastq --loglength -t 8 -o OUTDIR
