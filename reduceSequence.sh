#!/bin/bash

#PROJECT_PATH=$WORK/Data/NGS/JA16493
#PROJECT_PATH=$WORK/Data/NGS/JA16468/combined
#PROJECT_PATH=$WORK/Data/NGS/SA16172/JA16594
PROJECT_PATH=$DATA/JA16715
DATA_PATH=$PROJECT_PATH
RESULT_PATH=$DATA_PATH/splitted
SUFFIX=_R1_001.fastq.gz
PROGRAM=read_cluster_pairs.py
THREADS=12
mkdir -p  $RESULT_PATH

for FQ1 in `ls $DATA_PATH/*${SUFFIX}`
do
	SAMPLE_NAME=$(basename ${FQ1%$SUFFIX})
	if [[ $SAMPLE_NAME == cdk12* ]]
	then
		PRIMING_SITE=CGCCTTCGATATTGCTTCTTCGGTTTCATGGTGTTG
	elif [[ $SAMPLE_NAME == mdm4* ]]
	then
		PRIMING_SITE=CCCGTCTCGTGGTCTTTTCTCACATAAGCT
	fi

	FQ2=${FQ1/R1/R2}
	echo $(which python)  $PROGRAM \
		--outputprefix ${RESULT_PATH}/${SAMPLE_NAME}-errorFree \
	    --fastq1 ${FQ1} \
		--fastq2 ${FQ2} \
		--cutoff 0 \
		--idxBase 15 \
		--barcodeCutOff 20 \
		--constant_region $PRIMING_SITE \
		--threads $THREADS \
		--mismatch 3 \
		--read read2 \
		--fraction 0.66 \
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
