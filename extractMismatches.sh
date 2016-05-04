#!/bin/bash

PROJECTPATH=/scratch/cdw2854/tgirtDNA
DATAPATH=$PROJECTPATH/rawData
THREADS=24

for FQ1 in $DATAPATH/*R1_001.fastq.gz
do
    FQ2=${FQ1%R1_001.fastq.gz}R2_001.fastq.gz
    echo python mismatchExtraction.py \
        --fq1=$FQ1  \
        --fq2=$FQ2 \
        --outputPath=$PROJECTPATH \
        --threads=$THREADS
done | grep -v 'Error\|try'
