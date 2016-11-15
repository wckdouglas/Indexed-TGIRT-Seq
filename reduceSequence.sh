#!/bin/bash

#PROJECT_PATH=$WORK/Data/NGS/JA16493
#PROJECT_PATH=$WORK/Data/NGS/JA16468/combined
PROJECT_PATH=$WORK/Data/NGS/SA16172/JA16594
DATA_PATH=$PROJECT_PATH
RESULT_PATH=$DATA_PATH/splitted
SUFFIX=_R1_001.fastq.gz
PROGRAM=read_cluster_pairs.py
THREADS=12
mkdir -p  $RESULT_PATH

for FQ1 in `ls $DATA_PATH/*${SUFFIX}`
do
	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
	FQ2=${FQ1/R1/R2}
	echo $(which python)  $PROGRAM \
		--outputprefix ${RESULT_PATH}/${SAMPLE_NAME}-errorFree-95 \
	    --fastq1 ${FQ1} \
		--fastq2 ${FQ2} \
		--cutoff 0 \
		--idxBase 13 \
		--barcodeCutOff 4 \
		--constant_region CATCG \
		--threads $THREADS \
		--mismatch 1 \
		--read read1 \
		\&\> ${RESULT_PATH}/${SAMPLE_NAME}.log
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
