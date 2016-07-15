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
    descriptions = 'Clustering fastq reads to fasta reads with the first $IDXBASE bases as cDNA-synthesis barcode. ' +\
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
    parser.add_argument("-t", "--threads", default=1, type=int,
            help="Threads to use (default: 1)")
    args = parser.parse_args()
    return args

def concensusPairs(table):
    """ given a pair of reads as defined as the class: seqRecord
    return concensus sequence and mean quality of the pairs,
        as well as the number of reads that supports the concnesus pairs
    see function: concensusSeq, calculateConcensusBase
    """
    # get concensus left reads first
    sequenceLeft, qualityLeft = concensusSeq(table[0], table[2], range(len(table[0][0])))
    assert len(sequenceLeft) == len(qualityLeft), 'Wrong concensus sequence and quality!'
    # get concensus right reads first
    sequenceRight, qualityRight = concensusSeq(table[1], table[3], range(len(table[1][0])))
    assert len(sequenceRight) == len(qualityRight), 'Wrong concensus sequence and quality!'
    return sequenceLeft, qualityLeft, sequenceRight, qualityRight

def errorFreeReads(args):
    """
    main function for getting concensus sequences from read clusters.
    return  a pair of concensus reads with a 4-line fastq format
    see functions: 1. filterRead,
                  2. concensusPairs,
                  3. calculateConcensusBase
    """
    # skip if not enough sequences to perform voting
    index, h5_file, minReadCount = args
    with h5py.File(h5_file,'r') as h5:
        table = h5['barcodes'][index.upper()]
        member_count = table.shape[1]
        if member_count >= minReadCount:
            sequenceLeft, qualityLeft, sequenceRight, qualityRight = concensusPairs(table)
            leftRecord = '%s_%i_readCluster\n%s\n+\n%s\n' %(index, member_count, sequenceLeft, qualityLeft)
            rightRecord = '%s_%i_readCluster\n%s\n+\n%s\n' %(index, member_count, sequenceRight, qualityRight)
            return leftRecord, rightRecord

def readClustering(read1, read2, barcodeDict, idxBase, barcodeCutOff, constant, constant_length, hamming_threshold, usable_seq):
    """
    generate read cluster with a dictionary object and seqRecord class.
    index of the dictionary is the barcode extracted from first /idxBases/ of read 1
    """
    idLeft, seqLeft, qualLeft = read1
    idRight, seqRight, qualRight = read2
    assert idLeft.split(' ')[0] == idRight.split(' ')[0], 'Wrongly splitted files!! %s\n%s' %(idRight, idLeft)
    barcode = seqLeft[:idxBase]
    constant_region = seqLeft[idxBase:usable_seq]
    barcodeQualmean = int(np.mean(map(ord,qualLeft[:idxBase])) - 33)
    if ('N' not in barcode \
            and barcodeQualmean > barcodeCutOff \
            and not any(pattern in barcode for pattern in ['AAAAA','CCCCC','TTTTT','GGGGG']) \
            and hammingDistance(constant, constant_region) <= hamming_threshold):
        seqLeft = seqLeft[usable_seq:]
        qualLeft = qualLeft[usable_seq:]
        barcodeDict[barcode].addRecord(seqRight, qualRight, seqLeft, qualLeft)
        return 0
    return 1

def dictToh5File(barcodeDict, h5_file, barcode_file):
    """
    converting sequence dict to a h5 file for minimizing memory use
    """
    with h5py.File(h5_file,'w') as h5, open(barcode_file,'w') as bar_file:
        group = h5.create_group('barcodes')
        for index, seq_record in barcodeDict.iteritems():
            df = np.array([seq_record.seqListLeft,
                           seq_record.seqListRight,
                           seq_record.qualListLeft,
                           seq_record.qualListRight],
                           dtype='S256')
            table = group.create_dataset(index, data = df)
            bar_file.write(index + '\n')
    print 'Finished writting %s'  %h5_file
    return 0

def recordsToDict(outputprefix, inFastq1, inFastq2, idxBase, barcodeCutOff, constant):
    barcodeDict = defaultdict(seqRecord)
    read_num = 0
    constant_length = len(constant)
    hamming_threshold = float(1)/constant_length
    usable_seq = idxBase + constant_length
    discarded_sequence_count = 0
    with gzip.open(inFastq1,'rb') as fq1, gzip.open(inFastq2,'rb') as fq2:
        for read1,read2 in izip(FastqGeneralIterator(fq1),FastqGeneralIterator(fq2)):
            discarded_sequence_count += readClustering(read1,read2,barcodeDict, idxBase, barcodeCutOff,
                    constant, constant_length, hamming_threshold, usable_seq)
            read_num += 1
            if read_num % 1000000 == 0:
                stderr.write('[%s] Parsed: %i sequence\n' %(programname,read_num))
    stderr.write('[%s] Extracted: %i barcode sequences, discarded %i sequences\n' %(programname,len(barcodeDict.keys()), discarded_sequence_count))
    return barcodeDict, read_num

def writingAndClusteringReads(outputprefix, minReadCount, h5_file, threads, barcode_file):
    # From index library, generate error free reads
    # using multicore to process read clusters
    counter = 0
    output_cluster_count = 0
    read1File = outputprefix + '_R1_001.fastq.gz'
    read2File = outputprefix + '_R2_001.fastq.gz'
    with gzip.open(read1File,'wb') as read1, gzip.open(read2File,'wb') as read2, open(barcode_file,'r') as bar_file:
        args = ((str(index).strip(), h5_file, minReadCount) for index in bar_file)
        pool = Pool(threads)
        processes = pool.imap_unordered(errorFreeReads, args)
        #processes = imap(errorFreeReads, args)
        for p in processes:
            counter += 1
            if counter % 100000 == 0:
                stderr.write('[%s] Processed %i read clusters.\n' %(programname, counter))
            if p != None:
                leftRecord, rightRecord = p
                read1.write('@cluster%i_%s' %(output_cluster_count, leftRecord))
                read2.write('@cluster%i_%s' %(output_cluster_count, rightRecord))
                output_cluster_count += 1
    pool.close()
    pool.join()
    return output_cluster_count, read1File, read2File

def clustering(outputprefix, inFastq1, inFastq2, idxBase, minReadCount, barcodeCutOff, constant, threads):
    h5_file = outputprefix + '.h5'
    barcode_file = outputprefix + '.txt'
    barcodeDict, read_num = recordsToDict(outputprefix, inFastq1, inFastq2, idxBase, barcodeCutOff, constant)
    barcodeCount = map(lambda x: barcodeDict[x].member_count, barcodeDict.keys())
    p = plotBCdistribution(barcodeCount, outputprefix)
    dictToh5File(barcodeDict, h5_file, barcode_file)
    barcodeDict.clear()
    output_cluster_count, read1File, read2File = writingAndClusteringReads(outputprefix, minReadCount, h5_file, threads, barcode_file)
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
    idxBase = args.idxBase
    minReadCount = args.cutoff
    barcodeCutOff = args.barcodeCutOff
    constant = args.constant_region
    threads = args.threads

    #print out parameters
    stderr.write('[%s] [Parameters] \n' %(programname))
    stderr.write('[%s] indexed bases:                     %i\n' %(programname,idxBase))
    stderr.write('[%s] Threads:                           %i\n' %(programname,threads))
    stderr.write('[%s] minimum coverage:                  %i\n' %(programname,minReadCount))
    stderr.write('[%s] outputPrefix:                      %s\n' %(programname,outputprefix))
    stderr.write('[%s] using constant regions:   %s\n' %(programname,constant))

    # divide reads into subclusters
    clustering(outputprefix, inFastq1, inFastq2, idxBase, minReadCount, barcodeCutOff, constant, threads)
    stderr.write('[%s] time lapsed:      %2.3f min\n' %(programname, np.true_divide(time.time()-start,60)))
    return 0

if __name__ == '__main__':
    args = getOptions()
    main(args)