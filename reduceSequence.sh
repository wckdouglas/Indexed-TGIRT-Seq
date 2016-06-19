#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/jurkatCells
DATA_PATH=$PROJECT_PATH/rawData/completeData
RESULT_PATH=$PROJECT_PATH/splitted
LOG_PATH=$RESULT_PATH/logs
SUFFIX=_R1_001.fastq.gz
mkdir -p  $RESULT_PATH $LOG_PATH

for FQ1 in `ls $DATA_PATH/*${SUFFIX}`
do
	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
	FQ2=${FQ1/R1/R2}
	echo python readClusterPairs.py --outputprefix=$RESULT_PATH/$SAMPLE_NAME-errorFree \
	    --fastq1=$FQ1 --fastq2=$FQ2 --idxBase=13 --barcodeCutOff=30 \
	    --cutoff=0 --threads=12 --constant_region=CATCG 
done



