#!/bin/bash

PROJECT_PATH=/corral-repl/utexas/2013lambowitz/Data/JA16381
PROJECT_PATH=/stor/work/Lambowitz/Data/NGS/JA16381
DATA_PATH=$PROJECT_PATH/combined
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
		--outputprefix ${RESULT_PATH}/${SAMPLE_NAME}-errorFree \
	    --fastq1 ${FQ1} \
		--fastq2 ${FQ2} \
		--idxBase 13 \
		--barcodeCutOff 30 \
	    --cutoff 0 \
		--constant_region CATCG \
		--threads $THREADS
done

PROGRAM=double_index_cluster.py
for FQ1 in `ls $DATA_PATH/DB*${SUFFIX}`
do
	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
	FQ2=${FQ1/R1/R2}
<<<<<<< HEAD
	echo $(which python) $PROGRAM \
		--outputprefix=$RESULT_PATH/$SAMPLE_NAME-errorFree-double-BC \
	    --fastq1=$FQ1 --fastq2=$FQ2 \
		--idxBase=13 --barcodeCutOff=30 \
	    --cutoff 0 -l CATCG -r GAGTGTAGTGCATATGAGCACTGTCGAT \
=======
	echo $(which python) $PROGRAM --outputprefix=$RESULT_PATH/$SAMPLE_NAME-errorFree-double-BC \
	    --fastq1=$FQ1 --fastq2=$FQ2 --idxBase=13 --barcodeCutOff=30 \
	    --cutoff=0 -l CATCG -r GAGTGTAGTGCATATGAGCACTGTCGAT \
>>>>>>> ee0f5b525b8e059e3f9bdb3d4ee607aaa0f5d9a6
		--threads $THREADS
done
