#!/bin/bash

datapath=/scratch/02727/cdw2854/TGIRT_plasma_dna/syntheticOligo/rawData
resultpath=$datapath
suffix=fastq.gz

for fq1 in $datapath/*R1_001.$suffix
do
    samplename=$(echo $(basename $fq1) | cut -d'_' -f1)
    fq2=${fq1%R1_001.$suffix}R2_001.$suffix
    echo python readClusterPairs.py \
            --cutoff=6 \
            --threads=12 \
            --idxBase=13 \
            --vote=0.9 \
            --barcodeCutOff=33 \
            --loglikThreshold=4.5 \
            --seqError=0.02 \
            --outputprefix=$resultpath/$samplename-ErrorFree \
            --fastq1=$fq1 --fastq2=$fq2 --retainN
done | grep -v 'try\|ErrorFree.Error'

