#!/bin/bash

PROJECT_PATH=$WORK/Data/NGS/JA16493
PROJECT_PATH=$WORK/Data/NGS/JA16468/combined
DATA_PATH=$PROJECT_PATH
RESULT_PATH=$PROJECT_PATH/splitted
SUFFIX=_R1_001.fastq.gz
PROGRAM=r1_index_cluster.py
THREADS=12
mkdir -p  $RESULT_PATH

for FQ1 in `ls $DATA_PATH/*${SUFFIX}`
do
	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
	FQ2=${FQ1/R1/R2}
	echo $(which python)  $PROGRAM \
		--outputprefix ${RESULT_PATH}/${SAMPLE_NAME}-errorFree-cython \
	    --fastq1 ${FQ1} \
		--fastq2 ${FQ2} \
		--idxBase 13 \
		--barcodeCutOff 30 \
	    --cutoff 0 \
		--constant_region CATCG \
		--threads $THREADS
done

#PROGRAM=double_index_cluster.py
#for FQ1 in `ls $DATA_PATH/DB*${SUFFIX}`
#do
#	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
#	FQ2=${FQ1/R1/R2}
#	echo $(which python) $PROGRAM \
#		--outputprefix=$RESULT_PATH/$SAMPLE_NAME-errorFree-double-BC \
#	    --fastq1=$FQ1 --fastq2=$FQ2 \
#		--idxBase=13 --barcodeCutOff=30 \
#	    --cutoff 0 -l CATCG -r GAGTGTAGTGCATATGAGCACTGTCGAT \
#		--threads $THREADS
#done
