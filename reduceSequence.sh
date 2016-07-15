#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/jurkatCells
DATA_PATH=$PROJECT_PATH/rawData
RESULT_PATH=$DATA_PATH/splitted
SUFFIX=_R1_001.fastq.gz
PROGRAM=read_cluster_dict.py
mkdir -p  $RESULT_PATH

for FQ1 in `ls $DATA_PATH/*${SUFFIX}`
do
	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
	FQ2=${FQ1/R1/R2}
	echo $(which python) $PROGRAM \
		--outputprefix ${RESULT_PATH}/${SAMPLE_NAME}-errorFree-double-BC \
	    --fastq1 ${FQ1} \
		--fastq2 ${FQ2} \
		--idxBase 13 \
		--barcodeCutOff 30 \
	    --cutoff=0 \
		--constant_region=CATCG \
		--threads 24
done

#for FQ1 in `ls $DATA_PATH/*${SUFFIX}`
#do
#	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
#	FQ2=${FQ1/R1/R2}
#	echo $(which python) $PROGRAM --outputprefix=$RESULT_PATH/$SAMPLE_NAME-errorFree \
#	    --fastq1=$FQ1 --fastq2=$FQ2 --idxBase=13 --barcodeCutOff=30 \
#	    --cutoff=0 -l CATCG -r GAGTGTAGTGCATATGAGCACTGTCGAT
#done
