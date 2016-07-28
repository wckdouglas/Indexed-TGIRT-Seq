#!/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sys import stderr
import numpy as np
import sys
import argparse
import glob
import gzip
import time
import os
from itertools import izip, imap, product
from multiprocessing import Pool
from cluster_reads import *
from collections import defaultdict
import h5py
programname = os.path.basename(sys.argv[0]).split('.')[0]

def getOptions():
    '''
    reading input
    '''
    descriptions = 'Clustering fastq reads to fasta reads with the first $idx_base bases as cDNA-synthesis barcode. ' +\
                'Concensus bases are called only when the fraction of reads that contain the concensus base exceed some threshold. '+ \
                'Quality scores are generated by the average score for the bases that matched concensus base. '
    parser = argparse.ArgumentParser(description=descriptions)
    parser.add_argument('-o', '--outputprefix', required=True,
        help='Paired end Fastq files with R1_001.fastq.gz as suffix for read1, and R2_001.fastq.gz as suffix for read2')
    parser.add_argument('-1', '--fastq1', required=True,
        help='Paired end Fastq file 1 with four line/record')
    parser.add_argument('-2', '--fastq2',required=True,
        help='Paired end Fastq file 2 with four line/record')
    parser.add_argument('-m', '--cutoff', type=int,default=4,
        help="minimum read count for each read cluster (default: 4)")
    parser.add_argument("-x", "--idxBase", type=int, default=13,
        help="how many base in 5' end as index? (default: 13)")
    parser.add_argument('-q', '--barcodeCutOff', type=int, default=30,
        help="Average base calling quality for barcode sequence (default=30)")
    parser.add_argument("-c", "--constant_region", default='CATCG',
            help="Constant sequence after tags (default: CATCG ,e.g. Douglas's index-R1R)")
    args = parser.parse_args()
    return args

def readClustering(read1, read2, barcode_dict, idx_base, barcode_cut_off, constant, constant_length, hamming_threshold, usable_seq):
    """
    generate read cluster with a dictionary object and seqRecord class.
    index of the dictionary is the barcode extracted from first /idx_bases/ of read 1
    """
    idLeft, seqLeft, qualLeft = read1
    idRight, seqRight, qualRight = read2
    assert idLeft.split(' ')[0] == idRight.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(idRight, idLeft)
    barcode = seqLeft[:idx_base]
    constant_region = seqLeft[idx_base:usable_seq]
    barcodeQualmean = int(np.mean(map(ord,qualLeft[:idx_base])) - 33)
    if ('N' not in barcode \
            and barcodeQualmean > barcode_cut_off \
            and hammingDistance(constant, constant_region) <= hamming_threshold): #\
            #and not any(pattern in barcode for pattern in ['AAAAA','CCCCC','TTTTT','GGGGG']):
        seqLeft = seqLeft[usable_seq:]
        qualLeft = qualLeft[usable_seq:]
        barcode_dict[barcode].append([seqLeft,seqRight,qualLeft, qualRight])
        return 0
    return 1

def recordsToDict(outputprefix, inFastq1, inFastq2, idx_base, barcode_cut_off, constant, barcode_dict):
    read_num, discarded_sequence_count = 0, 0
    constant_length = len(constant)
    hamming_threshold = float(1)/constant_length
    usable_seq = idx_base + constant_length
    with gzip.open(inFastq1,'rb') as fq1, gzip.open(inFastq2,'rb') as fq2:
        for read1,read2 in izip(FastqGeneralIterator(fq1),FastqGeneralIterator(fq2)):
            discarded_sequence_count += readClustering(read1,read2,barcode_dict, idx_base, barcode_cut_off,
                    constant, constant_length, hamming_threshold, usable_seq)
            read_num += 1
            if read_num % 1000000 == 0:
                stderr.write('[%s] Parsed: %i sequence\n' %(programname,read_num))
    barcode_count = len(barcode_dict.keys())
    stderr.write('[%s] Extracted: %i barcode group\n' %(programname,barcode_count) +\
                 '[%s] discarded: %i sequences\n' %(programname, discarded_sequence_count) +\
                 '[%s] Parsed:    %i seqeucnes\n' %(programname, read_num))
    return barcode_dict, read_num, barcode_count

def clustering(outputprefix, inFastq1, inFastq2, idx_base, min_family_member_count, barcode_cut_off, constant):
    barcode_dict = defaultdict(list)
    barcode_dict, read_num, barcode_count = recordsToDict(outputprefix, inFastq1, inFastq2, idx_base, barcode_cut_off, constant, barcode_dict)
    barcode_member_counts = map(lambda index: len(barcode_dict[index]), barcode_dict.keys())
    p = plotBCdistribution(barcode_member_counts, outputprefix)
    dictToJson(barcode_dict, outputprefix)
    barcode_dict.clear()
    #output_cluster_count, read1File, read2File = writingAndClusteringReads(outputprefix, min_family_member_count, barcode_dict, barcode_count)
    output_cluster_count, read1File, read2File = writingAndClusteringReads(outputprefix, min_family_member_count, barcode_count)
    # all done!
    stderr.write('[%s] Finished writing error free reads\n' %programname)
    stderr.write('[%s] [Summary]                        \n' %programname)
    stderr.write('[%s] read1:                     %s\n' %(programname, read1File))
    stderr.write('[%s] read2:                     %s\n' %(programname, read2File))
    stderr.write('[%s] output clusters:           %i\n' %(programname, output_cluster_count))
    stderr.write('[%s] Percentage retained:       %.3f\n' %(programname, float(output_cluster_count)/read_num * 100))
    return 0

def main(args):
    """
    main function:
        controlling work flow
        1. generate read clusters by reading from fq1 and fq2
        2. obtain concensus sequence from read clusters
        3. writing concensus sequence to files
    """
    start = time.time()
    outputprefix = args.outputprefix
    inFastq1 = args.fastq1
    inFastq2 = args.fastq2
    idx_base = args.idxBase
    min_family_member_count = args.cutoff
    barcode_cut_off = args.barcodeCutOff
    constant = args.constant_region

    #print out parameters
    stderr.write('[%s] [Parameters] \n' %(programname))
    stderr.write('[%s] indexed bases:                     %i\n' %(programname,idx_base))
    stderr.write('[%s] minimum coverage:                  %i\n' %(programname,min_family_member_count))
    stderr.write('[%s] outputPrefix:                      %s\n' %(programname,outputprefix))
    stderr.write('[%s] using constant regions:   %s\n' %(programname,constant))

    # divide reads into subclusters
    clustering(outputprefix, inFastq1, inFastq2, idx_base, min_family_member_count, barcode_cut_off, constant)
    stderr.write('[%s] time lapsed:      %2.3f min\n' %(programname, np.true_divide(time.time()-start,60)))
    return 0

if __name__ == '__main__':
    args = getOptions()
    main(args)
