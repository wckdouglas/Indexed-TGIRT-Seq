#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/TGIRT_plasma_dna/syntheticOligo
FASTQPATH=$PROJECTPATH/trimmedFastq
BAMPATH=$PROJECTPATH/bamFile
BAMDATAPATH=$BAMPATH/bamData
EXECUTABLE=/work/02727/cdw2854/src/samData/bin/collectDataFromSam
INDEX=/corral-repl/utexas/2013lambowitz/Ref/syntheticOligos/celSpikeinLong/spikein.fa
THREADS=12
SUFFIX=1P.fastq.gz

for fq1 in $FASTQPATH/*$SUFFIX
do
    sample=$(basename ${fq1%_$SUFFIX})
    fq2=$FASTQPATH/$sample'_2P.fastq.gz'
    echo bowtie2 -p $THREADS --no-mixed \
            --very-sensitive-local --norc \
            -x $INDEX -1 $fq1 -2 $fq2 \
        \| samtools view -b@ $THREADS -F 4 -F 256 -F 2048 -  \
        \| samtools sort -@ $THREADS -O bam -T $BAMPATH/$sample - \
        \> $BAMPATH/$sample.bam 
done
