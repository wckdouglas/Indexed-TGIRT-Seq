#!/bin/bash

DATAPATH=/scratch/02727/cdw2854/TGIRT_plasma_dna/syntheticOligo/bamFile
RESULTPATH=$DATAPATH/mismatchData
REFFASTA=/corral-repl/utexas/2013lambowitz/Ref/syntheticOligos/celSpikeinLong/spikein.fa
EXEC=/work/02727/cdw2854/src/filterSamFile/bin/pileup2bed
#EXEC=/work/02727/cdw2854/src/filterSamFile/bin/parsePileup
SAMTOOLS=/opt/apps/samtools/0.1.19/samtools
QUALITY=0
MINCOV=0
#MINCOV=100
#QUALITY=30
DEPTH=10000000
THREADS=12

for BAM in $DATAPATH/*bam
do
    SAMPLE=$(basename ${BAM%.bam})
#    echo $SAMTOOLS mpileup \
#        -d $DEPTH \
#        -f $REFFASTA \
#        -q 15 -Q $QUALITY $BAM \
#        \| $EXEC - $QUALITY $MINCOV \
#        \| awk \'\$7!=\"-\"\' \
#        \> $RESULTPATH/$SAMPLE.cpp.tsv
    echo time python pileupBam.py \
        --refFasta=$REFFASTA \
        --qualThresh=30 --threads=$THREADS \
        --depth=$DEPTH --bamfile=$BAM \
        \> $RESULTPATH/$SAMPLE.tsv
done

