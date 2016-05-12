#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/TGIRT_DNA
DATA_PATH=$PROJECT_PATH/rawData
RESULT_PATH=$PROJECT_PATH/splitted
SUFFIX=_R1_001.fastq.gz
mkdir -p logs

for FQ1 in `ls $DATA_PATH/*${SUFFIX}`
do
	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
	FQ2=${FQ1/R1/R2}
	echo python readClusterPairs.py -o $RESULT_PATH/$SAMPLE_NAME-errorFree \
									-1 $FQ1 -2 $FQ2 -x 13 -q 30 \
									-m 6 -t 12 \
		\&\> logs/${SAMPLE_NAME}.log
done



